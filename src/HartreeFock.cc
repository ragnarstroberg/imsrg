
#include "HartreeFock.hh"
#include "ModelSpace.hh"

using namespace std;

template<class OPERATOR>
HartreeFock<OPERATOR>::HartreeFock(OPERATOR& hbare) : Hbare(hbare)
{
   tolerance = 1e-10;
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();

   C         = arma::mat(norbits,norbits,arma::fill::eye);
   Vab       = arma::mat(norbits,norbits,arma::fill::zeros);
   V3ab       = arma::mat(norbits,norbits,arma::fill::zeros);
   H         = arma::mat(norbits,norbits);
   Vmon      = arma::mat(nKets,nKets);
   Vmon_exch = arma::mat(nKets,nKets);
   prev_energies = arma::vec(norbits,arma::fill::zeros);

   t = Hbare.OneBody;
   energies = t.diag();
   BuildMonopoleV();
   cout << "About to build Vmon3." << endl;
   BuildMonopoleV3();
   cout << "Done building Vmon3." << endl;
   UpdateDensityMatrix();
   UpdateH();
}
template HartreeFock<Operator>::HartreeFock(Operator& );
template HartreeFock<Operator3>::HartreeFock(Operator3& );

template<class OPERATOR>
void HartreeFock<OPERATOR>::Solve()
{
   iterations = 0; // counter so we don't go on forever
   int maxiter = 1000;

   while( (! CheckConvergence()) && iterations < maxiter)
   {
      Diagonalize();  // Diagonalize the Fock matrix
      ReorderCoefficients(); // Reorder C so that new ordering corresponds to original ordering
      UpdateDensityMatrix();  // 1 body density matrix, used in UpdateH()
      UpdateH();  // Update the Fock matrix
      ++iterations;
   }

   ReorderCoefficients(); // Reorder C so that new ordering corresponds to original ordering
   if (iterations==maxiter)
   {
      cerr << "Warning: Hartree-Fock calculation didn't converge after " << maxiter << " iterations." << endl;
   }
   CalcEHF();
}
template void HartreeFock<Operator>::Solve();
template void HartreeFock<Operator3>::Solve();


//*******************************************************
// Calculate the HF energy.
//*******************************************************
template<class OPERATOR>
void HartreeFock<OPERATOR>::CalcEHF()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   EHF = 0;
   int norbits = Hbare.GetModelSpace()->GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      int jfactor = ms->GetOrbit(i).j2 +1;
      for (int j=0;j<norbits;j++)
      {
//         EHF += 0.5 * rho(i,j) * jfactor * (H(i,j)+t(i,j));
         EHF +=  rho(i,j) * jfactor * (t(i,j)+0.5*Vab(i,j)+1./6*V3ab(i,j));
      }
   }
}
template void HartreeFock<Operator>::CalcEHF();
template void HartreeFock<Operator3>::CalcEHF();


///**********************************************************************************
// Diagonalize() -- [See Suhonen eq. 4.85]
// Diagonalize <a|H|b> and put the
// eigenvectors in C(i,alpha) = <alpha|i>
// and eigenvalues in the vector energies.
// Save the last vector of energies to check
// for convergence.
// Different channels are diagonalized independently.
// This guarantees that J,Tz, and parity remain good. 
///***********************************************************************************
template<class OPERATOR>
void HartreeFock<OPERATOR>::Diagonalize()
{
   prev_energies = energies;
   vector<int> orbit_list;
   arma::mat H_ch;
   arma::mat C_ch;
   arma::vec E_ch;
   int norbits = Hbare.GetModelSpace()->GetNumberOrbits();
   int Jmax1 = Hbare.GetModelSpace()->OneBodyJmax;
   for (int p = 0; p<=1;p++)
   {
      for (int Tz = -1; Tz<=1; Tz+=2)
      {
          for (int J=0;J<=Jmax1;J++)
          {

            // Find all the SP orbits that have quantum numbers J,p,Tz
            // and store them in a list
             orbit_list.resize(0);
             for (int a=0;a<norbits;a++)
             {
                Orbit &oa = Hbare.GetModelSpace()->GetOrbit(a);
                if (oa.j2==J and oa.tz2==Tz and (oa.l%2)==p)
                {
                   orbit_list.push_back(a);
                }
             }
             // Now create submatrices corresponding to just these orbits
             int norb = orbit_list.size();
             if (norb < 1) continue;
             H_ch = arma::mat(norb,norb,arma::fill::zeros);
             C_ch = arma::mat(norb,norb,arma::fill::zeros);
             E_ch = arma::vec(norb,arma::fill::zeros);
             for (int a=0;a<norb;a++)
             {
                for (int b=0;b<norb;b++)
                {
                   H_ch(a,b) = H(orbit_list[a],orbit_list[b]);
                }
             }
             // Diagonalize the submatrix
             bool success = false;
             int diag_tries = 0;
             while ( not success and diag_tries<5 )
             {
                success = arma::eig_sym(E_ch, C_ch, H_ch);
                ++diag_tries;
             }
             if (not success)
             {
                cout << "Hartree-Fock: Failed to diagonalize the submatrix with J=" << J << " Tz=" << Tz << " parity = " << p
                     << " on iteration # " << iterations << ". The submatrix looks like:" << endl;
                H_ch.print();
             }
             // Update the full overlap matrix C and energy vector
             for (int a=0;a<norb;a++)
             {
                energies( orbit_list[a] ) = E_ch(a);
                for (int b=0;b<norb;b++)
                {
                   C(orbit_list[a],orbit_list[b]) = C_ch(a,b);
                }
             }

          } // for J...
      } // For Tz ...
   } // For p...
}
template void HartreeFock<Operator>::Diagonalize();
template void HartreeFock<Operator3>::Diagonalize();

//**************************************************
// BuildMonopoleV()
// Construct the unnormalized monople Hamiltonian
// <ab|V_mon|cd> = Sum_J (2J+1) <ab|V|cd>_J.
//         
//**************************************************
template<class OPERATOR>
void HartreeFock<OPERATOR>::BuildMonopoleV()
{
   Vmon.zeros();
   Vmon_exch.zeros();
   int nKets = Hbare.GetModelSpace()->GetNumberKets();
   for (int ibra=0;ibra<nKets;++ibra)
   {
      Ket & bra = Hbare.GetModelSpace()->GetKet(ibra);
      int a = bra.p;
      int b = bra.q;
      Orbit & oa = Hbare.GetModelSpace()->GetOrbit(a);
      Orbit & ob = Hbare.GetModelSpace()->GetOrbit(b);
      double norm = (oa.j2+1)*(ob.j2+1);
      for (int iket=ibra;iket<nKets;++iket)
      {
         Ket & ket = Hbare.GetModelSpace()->GetKet(iket);
         int c = ket.p;
         int d = ket.q;
         Vmon(ibra,iket)       = Hbare.GetTBMEmonopole(a,b,c,d)*norm;
         Vmon_exch(ibra,iket)  = Hbare.GetTBMEmonopole(a,b,d,c)*norm;
      }
   }
}
template void HartreeFock<Operator>::BuildMonopoleV();
template void HartreeFock<Operator3>::BuildMonopoleV();



// This is specific to 3 body interaction
//
// <iab|Vmon3|jcd> = sum_Jia,J (2J+1)sum_tia,T (2T+1)/2 <iab|V3(Jia,Jia,J,tia,tia,T)|jcd>
//
template<>
void HartreeFock<Operator>::BuildMonopoleV3()
{}

template<>
void HartreeFock<Operator3>::BuildMonopoleV3()
{
   ModelSpace * modelspace = Hbare.GetModelSpace();
   int norbits = modelspace->GetNumberOrbits();
   for (int i=0; i<norbits; ++i)
   {
     Orbit& oi = modelspace->GetOrbit(i);
     for (int a=0; a<norbits; ++a)
     {
       Orbit& oa = modelspace->GetOrbit(a);
       for (int b=0; b<norbits; ++b)
       {
         Orbit& ob = modelspace->GetOrbit(b);
         if ( 2*(oi.n+oa.n+ob.n)+oi.l+oa.l+ob.l > Hbare.E3max ) continue;
         for (int j=0; j<norbits; ++j)
         {
           Orbit& oj = modelspace->GetOrbit(j);
           if (oi.j2 != oj.j2 or oi.tz2 != oj.tz2 or oi.l != oj.l) continue;
           for (int c=0; c<norbits; ++c)
           {
             Orbit& oc = modelspace->GetOrbit(c);
             if (oa.j2 != oc.j2 or oa.tz2 != oc.tz2 or oa.l != oc.l) continue;
             for (int d=0; d<norbits; ++d)
             {
               Orbit& od = modelspace->GetOrbit(d);
               if ( 2*(oj.n+oc.n+od.n)+oj.l+oc.l+od.l > Hbare.E3max ) continue;
               if (ob.j2 != od.j2 or ob.tz2 != od.tz2 or ob.l != od.l) continue;
               if ( (oi.l+oa.l+ob.l+oj.l+oc.l+od.l)%2 >0) continue;
               unsigned long long orbit_index_pn =  i*1000000000000000LL
                                                  + a*1000000000000LL
                                                  + b*1000000000
                                                  + j*1000000
                                                  + c*1000
                                                  + d;
              double v = 0;
             int Jmin = max(max(1,oi.j2-oa.j2-ob.j2), oj.j2-oc.j2-od.j2);
             int Jmax = min(oi.j2+oa.j2+ob.j2, oj.j2+oc.j2+od.j2);
              for (int J=Jmin; J<=Jmax; J+=2)
              {
                int j2min = max( max(abs(oi.j2-oa.j2), abs(J-ob.j2)),max(abs(oj.j2-oc.j2), abs(J-od.j2) )) /2;
                int j2max = min( min(oi.j2+oa.j2, J+ob.j2), min(oj.j2+oc.j2, J+od.j2) )/2;
                for ( int j2=j2min; j2<j2max; ++j2)
                {
                  for (int T=1; T<=3; T+=2)
                  for ( int t2 = (T-1)/2; t2<=1; ++t2)
                  {
                     {
                       v += Hbare.GetThreeBodyME(j2,j2,J,t2,t2,T,i,a,b,j,c,d);
                     }
                  }
                }
              }
              Vmon3[orbit_index_pn] = v;
             }
           }
         }
       }
     }
   }

}




//**************************************************************************
// 1-body density matrix 
// <i|rho|j> = Sum_beta <i|beta> <j|beta>
// where beta runs over HF orbits in the core (i.e. below the fermi surface)
//**************************************************************************
template<class OPERATOR>
void HartreeFock<OPERATOR>::UpdateDensityMatrix()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();

   // Pcore is a projector onto core orbits.
   arma::mat Pcore = arma::mat(norbits,norbits,arma::fill::zeros);
   for (int &beta : ms->holes)
   {
       Pcore(beta,beta) = 1;
   }

   rho = C * Pcore * C.t();
}
template void HartreeFock<Operator>::UpdateDensityMatrix();
template void HartreeFock<Operator3>::UpdateDensityMatrix();


//*********************************************************************
// UpdateH() -- [See Suhonen eq 4.85]
// <a|H|b> = <a|t|b>  +  Sum_ij <i|rho|j> <ai|V_mon|bj>
// * H is the fock matrix, to be diagonalized
// * t is the kinetic energy
// * rho is the density matrix defined in UpdateDensityMatrix()
// * V_mon is the monopole component of the 2-body interaction.
//*********************************************************************
template <>
void HartreeFock<Operator>::UpdateH()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();
   int bra, ket;
   Vab.zeros();

   for (int a=0;a<norbits;a++)
   {
      Orbit& oa = ms->GetOrbit(a);
      for (int b=a;b<norbits;b++)
      {
         Orbit& ob = ms->GetOrbit(b);
         if (oa.j2 != ob.j2 or oa.tz2 != ob.tz2 or oa.l != ob.l)   continue;
         for (int i=0;i<norbits;i++)
         {
            Orbit& oi = ms->GetOrbit(i);
            bra = ms->GetKetIndex(min(a,i),max(a,i));
            for (int j=0;j<norbits;j++)
            {
               Orbit& oj = ms->GetOrbit(j);
               if (oi.j2 != oj.j2 or oi.tz2 != oj.tz2 or oi.l != oj.l)   continue;
               ket = ms->GetKetIndex(min(b,j),max(b,j));

               if (a>i xor b>j)
                  Vab(a,b) += rho(i,j)*Vmon_exch(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>
               else
                  Vab(a,b) += rho(i,j)*Vmon(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>
           }
         }
         Vab(a,b) /= (oa.j2+1);
         Vab(b,a) = Vab(a,b);  // Hermitian & real => symmetric
      }
   }
   H = t + Vab;
}



template<>
void HartreeFock<Operator3>::UpdateH()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();
   int bra, ket;
   Vab.zeros();
   V3ab.zeros();

   for (int a=0;a<norbits;a++)
   {
      Orbit& oa = ms->GetOrbit(a);
      for (int b=a;b<norbits;b++)
      {
         Orbit& ob = ms->GetOrbit(b);
         if (oa.j2 != ob.j2 or oa.tz2 != ob.tz2 or oa.l != ob.l)   continue;
         for (int i=0;i<norbits;i++)
         {
            Orbit& oi = ms->GetOrbit(i);
            bra = ms->GetKetIndex(min(a,i),max(a,i));
            for (int j=0;j<norbits;j++)
            {
               Orbit& oj = ms->GetOrbit(j);
               if (oi.j2 != oj.j2 or oi.tz2 != oj.tz2 or oi.l != oj.l)   continue;
               ket = ms->GetKetIndex(min(b,j),max(b,j));
               // 2body term <ai|V|bj>
               if (a>i xor b>j)
                  Vab(a,b) += rho(i,j)*Vmon_exch(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>
               else
                  Vab(a,b) += rho(i,j)*Vmon(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>

               // 3body term  <aik|V|bjl> <i|rho|j> <k|rho|l>
               for (int k=0;k<norbits;++k)
               {
                 Orbit& ok = ms->GetOrbit(k);
                 if ( 2*(oa.n+oi.n+ok.n)+oa.l+oi.l+ok.l > Hbare.E3max ) continue;
                 for (int l=0;l<norbits;++l)
                 {
                   Orbit& ol = ms->GetOrbit(l);
                   if ( 2*(ob.n+oj.n+ol.n)+ob.l+oj.l+ol.l > Hbare.E3max ) continue;
                   if (ok.j2 != ol.j2 or ok.tz2 != ol.tz2 or ok.l != ol.l) continue;
                   if ( (oa.l+oi.l+ok.l+ob.l+oj.l+ol.l)%2 >0) continue;
                   unsigned long long orbit_index_pn =  a*1000000000000000LL
                                                      + i*1000000000000LL
                                                      + k*1000000000
                                                      + b*1000000
                                                      + j*1000
                                                      + l;
                   V3ab(a,b) += rho(i,j) * rho(k,l) * Vmon3[orbit_index_pn];
//                   Vab(a,b) += rho(i,j) * rho(k,l) * Vmon3[orbit_index_pn];
                 }
               }
           }
         }
         Vab(a,b) /= (oa.j2+1);
         V3ab(a,b) /= (oa.j2+1);
         Vab(b,a) = Vab(a,b);  // Hermitian & real => symmetric
         V3ab(b,a) = V3ab(a,b);  // Hermitian & real => symmetric
      }
   }
//   H = t + Vab;
   H = t + Vab + 0.5*V3ab;
}



//********************************************************
// Check for convergence using difference in s.p. energies
// between iterations.
//********************************************************
template<class OPERATOR>
bool HartreeFock<OPERATOR>::CheckConvergence()
{
   arma::vec de = energies - prev_energies;
   double ediff = sqrt(arma::dot(de,de)) / energies.size();
   return (ediff < tolerance);
}
template bool HartreeFock<Operator>::CheckConvergence();
template bool HartreeFock<Operator3>::CheckConvergence();



template<class OPERATOR>
void HartreeFock<OPERATOR>::PrintOrbits()
{
  ModelSpace * ms = Hbare.GetModelSpace();
  int norbits = ms->GetNumberOrbits();
  for (int a=0;a<norbits;a++)
  {
     Orbit & oa = ms->GetOrbit(a);
     cout << a << ": E= " << oa.spe << " L=" << oa.l
          << " J=" << oa.j2 << " N=" << oa.n << " Tz=" << oa.tz2
          << " ph=" << oa.ph << " io=" << oa.io << endl;
  }
}
template void HartreeFock<Operator>::PrintOrbits();
template void HartreeFock<Operator3>::PrintOrbits();

//**********************************************************************
// Eigenvectors/values come out of the diagonalization energy-ordered.
// We want them ordered corresponding to the input ordering, i.e. we want
// the matrix C to be maximal and positive along the diagonal.
//**********************************************************************
template<class OPERATOR>
void HartreeFock<OPERATOR>::ReorderCoefficients()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   arma::mat C_tmp = C;
   arma::vec e_tmp = energies;
   for (int i=0;i<norbits;i++)
   {
      float fmax = 0.0;
      int kmax;
      int sign = 1;
      for (int k=0;k<norbits;k++)
      {
         if ( abs(C_tmp(i,k)) > fmax )
         {
            fmax = abs(C_tmp(i,k));
            kmax = k;
            if (C_tmp(i,k) < 0) sign = -1;
         }
      }
      // make sure we have a positive coefficient for the biggest term
      C.col(kmax) = C_tmp.col(i) * sign;
      energies(kmax) = e_tmp(i);
   }
}
template void HartreeFock<Operator>::ReorderCoefficients();
template void HartreeFock<Operator3>::ReorderCoefficients();




//**************************************************************************
// Takes in an operator expressed in the basis of the original Hamiltonian,
// and returns that operator in the Hartree-Fock basis.
//**************************************************************************
//Operator HartreeFock::TransformToHFBasis( Operator& OpIn)
template<class OPERATOR>
OPERATOR HartreeFock<OPERATOR>::TransformToHFBasis( OPERATOR& OpIn)
{
   OPERATOR OpHF = OpIn;

   // Easy part:
   //Update the one-body part by multiplying by the matrix C(i,a) = <i|a>.
   OpHF.OneBody = C.t() * OpIn.OneBody * C;

   // Medium part:
   //Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz. Most of the effort here is in constructing D.
   int nchan = OpIn.GetModelSpace()->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = OpIn.GetModelSpace()->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      if (npq<1) continue;

      arma::mat D     = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
      arma::mat Dexch = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ba> = <ji|ab>

      // loop over all possible original basis configurations <pq| in this J,p,Tz channel.
      // and all possible HF configurations |p'q'> in this J,p,Tz channel                                    
      // bra is in the original basis, ket is in the HF basis                                              
      // i and j are the indices of the matrix D for this channel                    
      for (int i=0; i<npq; ++i)    
      {
         Ket & bra = tbc.GetKet(i);   
         for (int j=0; j<npq; ++j)    
         {
            Ket & ket = tbc.GetKet(j); // 
            double normfactor = sqrt((1.0+ket.delta_pq())/(1.0+bra.delta_pq()));
            D(i,j) = C(bra.p,ket.p) * C(bra.q,ket.q) * normfactor;
            Dexch(i,j) = C(bra.p,ket.q) * C(bra.q,ket.p) * ket.Phase(tbc.J) * (1-ket.delta_pq()) * normfactor;
         }
      }

     // Do all the matrix multiplication in one expression so Armadillo can do optimizations.
     OpHF.TwoBody[ch]   = D.t()      * OpIn.TwoBody[ch] * D
                        + Dexch.t()  * OpIn.TwoBody[ch] * D
                        + D.t()      * OpIn.TwoBody[ch] * Dexch
                        + Dexch.t()  * OpIn.TwoBody[ch] * Dexch;

   }
   return OpHF;
}
template Operator HartreeFock<Operator>::TransformToHFBasis(Operator&);
template Operator3 HartreeFock<Operator3>::TransformToHFBasis(Operator3&);






