
#include "HartreeFock.hh"
#include "ModelSpace.hh"

#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif


using namespace std;

HartreeFock::HartreeFock(Operator& hbare) : Hbare(hbare)
{
   tolerance = 1e-10;
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();

   C         = arma::mat(norbits,norbits,arma::fill::eye);
   Vij       = arma::mat(norbits,norbits,arma::fill::zeros);
   V3ij       = arma::mat(norbits,norbits,arma::fill::zeros);
   F         = arma::mat(norbits,norbits);
   Vmon      = arma::mat(nKets,nKets);
   Vmon_exch = arma::mat(nKets,nKets);
   prev_energies = arma::vec(norbits,arma::fill::zeros);

   t = Hbare.OneBody;
   energies = t.diag();
   BuildMonopoleV();
   if (hbare.GetParticleRank()>2)
   {
      BuildMonopoleV3();
   }
   UpdateDensityMatrix();
   UpdateF();
}

void HartreeFock::Solve()
{
   iterations = 0; // counter so we don't go on forever
   int maxiter = 1000;

   for (iterations=0; iterations<maxiter; ++iterations)
   {
      Diagonalize();  // Diagonalize the Fock matrix
//      ReorderCoefficients(); // Reorder C so that new ordering corresponds to original ordering
      UpdateDensityMatrix();  // 1 body density matrix, used in UpdateF()
      UpdateF();  // Update the Fock matrix
      if ( CheckConvergence() ) break;
   }

   if (iterations==maxiter)
   {
      cerr << "Warning: Hartree-Fock calculation didn't converge after " << maxiter << " iterations." << endl;
   }

   ReorderCoefficients(); // Reorder C so that new ordering corresponds to original ordering
   CalcEHF();
// Just to check. Remove this...
//   Diagonalize();
//   UpdateDensityMatrix();
//   UpdateH();
//   CalcEHF();
}


/// Calculate the HF energy.
/// \f{eqnarray*} E_{HF} &=& \sum_{\alpha} t_{\alpha\alpha} 
///                    + \frac{1}{2}\sum_{\alpha\beta} V_{\alpha\beta\alpha\beta}
///                    + \frac{1}{6}\sum_{\alpha\beta\gamma} V_{\alpha\beta\gamma\alpha\beta\gamma} \\
///    &=& \sum_{ij} (2j_i+1) \rho_{ij} ( t_{ij} +\frac{1}{2}\tilde{V}^{(2)}_{ij} + \frac{1}{6}\tilde{V}^{(3)}_{ij} )
/// \f}
/// Where the matrices \f{eqnarray*}
///  \tilde{V}^{(2)}_{ij} &=& \sum_{ab} \rho_{ab}\bar{V}^{(2)}_{iajb} \\
///  \tilde{V}^{(2)}_{ij} &=& \sum_{abcd} \rho_{ab}\rho_{cd} \bar{V}^{(3)}_{iacjbd} \\
/// \f}
/// have already been calculated by UpdateF().
///
void HartreeFock::CalcEHF()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   EHF = 0;
   double e3hf = 0;
   int norbits = Hbare.GetModelSpace()->GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      int jfactor = ms->GetOrbit(i).j2 +1;
      for (int j=0;j<norbits;j++)
      {
//         EHF += 0.5 * rho(i,j) * jfactor * (H(i,j)+t(i,j));
         EHF +=  rho(i,j) * jfactor * (t(i,j)+0.5*Vij(i,j)+1./6*V3ij(i,j));
         e3hf +=  rho(i,j) * jfactor * (1./6*V3ij(i,j));
      }
   }
   cout << "e3hf = " << e3hf << endl;
}


/// [See Suhonen eq. 4.85]
/// Diagonalize the fock matrix \f$ <a|F|b> \f$ and put the
/// eigenvectors in \f$C(i,\alpha) = <i|\alpha> \f$
/// and eigenvalues in the vector energies.
/// Save the last vector of energies to check
/// for convergence.
/// Different channels are diagonalized independently.
/// This guarantees that J,Tz, and \f$ \pi \f$ remain good. 
void HartreeFock::Diagonalize()
{
   prev_energies = energies;
   vector<int> orbit_list;
   arma::mat F_ch;
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
             F_ch = arma::mat(norb,norb,arma::fill::zeros);
             C_ch = arma::mat(norb,norb,arma::fill::zeros);
             E_ch = arma::vec(norb,arma::fill::zeros);
             for (int a=0;a<norb;a++)
             {
                for (int b=0;b<norb;b++)
                {
                   F_ch(a,b) = F(orbit_list[a],orbit_list[b]);
//                   H_ch(a,b) = H(orbit_list[a],orbit_list[b]);
                }
             }
             // Diagonalize the submatrix
             bool success = false;
             int diag_tries = 0;
             while ( not success and diag_tries<5 )
             {
                success = arma::eig_sym(E_ch, C_ch, F_ch);
                ++diag_tries;
             }
             if (not success)
             {
                cout << "Hartree-Fock: Failed to diagonalize the submatrix with J=" << J << " Tz=" << Tz << " parity = " << p
                     << " on iteration # " << iterations << ". The submatrix looks like:" << endl;
                F_ch.print();
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


//*********************************************************************
/// Construct an unnormalized two-body monopole interaction
/// \f[ \langle ab | \bar{V}^{(2)} | cd \rangle = \sqrt{(1+\delta_{ab})(1+\delta_{cd})} \sum_{J} (2J+1) \langle ab | V^{(2)} | cd \rangle_{J} \f]
/// This method utilizes the operator method  Operator::GetTBMEmonopole() 
///
//*********************************************************************
void HartreeFock::BuildMonopoleV()
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



//*********************************************************************
/// Construct an unnormalized three-body monopole interaction
/// \f[ \langle iab | \bar{V}^{(3)} | jcd \rangle =
///     \sum\limits_{J,J_{12}}\sum_{Tt_{12}}(2J+1)(2T+1) 
///       \langle (ia)J_{12}t_{12};b JT| V^{(3)} | (jc)J_{12}t_{12}; d JT\rangle \f]
///
//*********************************************************************
void HartreeFock::BuildMonopoleV3()
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




//*********************************************************************
/// one-body density matrix 
/// \f$ <i|\rho|j> = \sum_{\beta} <i|\beta> <\beta|j> \f$
/// where beta runs over HF orbits in the core (i.e. below the fermi surface)
//*********************************************************************
void HartreeFock::UpdateDensityMatrix()
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


//*********************************************************************
///  [See Suhonen eq 4.85]
/// \f[ F_{ij} = t_{ij}  +  \frac{1}{2j_i+1}\sum_{ab} \rho_{ab} \bar{V}^{(2)}_{iajb}
///               + \frac{1}{2(2j_i+1)}\sum_{abcd}\rho_{ab} \rho_{cd} \bar{V}^{(3)}_{iacjbd}  \f]
/// * \f$ F \f$ is the Fock matrix, to be diagonalized
/// * \f$ t \f$ is the kinetic energy
/// * \f$\rho\f$ is the density matrix defined in UpdateDensityMatrix()
/// * \f$ \bar{V}^{(2)} \f$ is the monopole component of the 2-body interaction defined in BuildMonopoleV().
/// * \f$ \bar{V}^{(3)} \f$ is the monopole component of the 3-body interaction devined in BuildMonopoleV3().
//*********************************************************************
//void HartreeFock::UpdateH()
void HartreeFock::UpdateF()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();
   int bra, ket;
   Vij.zeros();
   V3ij.zeros();

   for (int i=0;i<norbits;i++)
   {
      Orbit& oi = ms->GetOrbit(i);
      for (int j=i;j<norbits;j++)
      {
         Orbit& oj = ms->GetOrbit(j);
         if (oi.j2 != oj.j2 or oi.tz2 != oj.tz2 or oi.l != oj.l)   continue;
         for (int a=0;a<norbits;++a)
         {
            Orbit& oa = ms->GetOrbit(a);
            bra = ms->GetKetIndex(min(i,a),max(i,a));
            for (int b=0;b<norbits;b++)
            {
               Orbit& ob = ms->GetOrbit(b);
               if (oa.j2 != ob.j2 or oa.tz2 != ob.tz2 or oa.l != ob.l)   continue;
               ket = ms->GetKetIndex(min(j,b),max(j,b));
               // 2body term <ai|V|bj>
               if (a>i xor b>j)
                  Vij(i,j) += rho(a,b)*Vmon_exch(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>
               else
                  Vij(i,j) += rho(a,b)*Vmon(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>

               if (Hbare.GetParticleRank()<3) continue;
               // 3body term  <aik|V|bjl> <i|rho|j> <k|rho|l>
               for (int c=0;c<norbits;++c)
               {
                 Orbit& oc = ms->GetOrbit(c);
                 if ( 2*(oi.n+oa.n+oc.n)+oi.l+oa.l+oc.l > Hbare.E3max ) continue;
                 for (int d=0;d<norbits;++d)
                 {
                   Orbit& od = ms->GetOrbit(d);
                   if ( 2*(oj.n+ob.n+od.n)+oj.l+ob.l+od.l > Hbare.E3max ) continue;
                   if (oc.j2 != od.j2 or oc.tz2 != od.tz2 or oc.l != od.l) continue;
                   if ( (oi.l+oa.l+oc.l+oj.l+ob.l+od.l)%2 >0) continue;
                   unsigned long long orbit_index_pn =  i*1000000000000000LL
                                                      + a*1000000000000LL
                                                      + c*1000000000
                                                      + j*1000000
                                                      + b*1000
                                                      + d;
                   V3ij(i,j) += rho(a,b) * rho(c,d) * Vmon3[orbit_index_pn];
                 }
               }
           }
         }
         Vij(i,j) /= (oi.j2+1);
         V3ij(i,j) /= (oi.j2+1)*2;
         Vij(j,i) = Vij(i,j);  // Hermitian & real => symmetric
         V3ij(j,i) = V3ij(i,j);  // Hermitian & real => symmetric
      }
   }
   F = t + Vij + 0.5*V3ij;
}



//********************************************************
/// Check for convergence using difference in s.p. energies
/// between iterations.
/// Converged when
/// \f[ \delta_{e} \equiv \sqrt{ \sum_{i}(e_{i}^{(n)}-e_{i}^{(n-1)})^2} < \textrm{tolerance} \f]
/// where \f$ e_{i}^{(n)} \f$ is the \f$ i \f$th eigenvalue of the Fock matrix after \f$ n \f$ iterations.
///
//********************************************************
bool HartreeFock::CheckConvergence()
{
   arma::vec de = energies - prev_energies;
   double ediff = sqrt(arma::dot(de,de)) / energies.size();
   return (ediff < tolerance);
}



void HartreeFock::PrintOrbits()
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

//**********************************************************************
/// Eigenvectors/values come out of the diagonalization energy-ordered.
/// We want them ordered corresponding to the input ordering, i.e. we want
/// the matrix C to be maximal and positive along the diagonal.
//**********************************************************************
void HartreeFock::ReorderCoefficients()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();

   int nswaps = 10;

   // first, reorder them so we still know the quantum numbers
   while (nswaps>0) // loop until we don't have to make any more swaps
   {
     nswaps = 0;
     for (int i=0;i<C.n_rows;++i) // loop through rows -> original basis states
     {
        arma::rowvec row = C.row(i);
        arma::uword imax; 
        double maxval = abs(row).max(imax);
        if (imax == i) continue;
        if (maxval > abs(C(imax,imax) ))
        {
           C.swap_cols(i,imax);
           energies.swap_rows(i,imax);
           ++nswaps;
        }
     }
   }
   // Make sure the diagonal terms are positive. (For easier comparison).
   for (int i=0;i<C.n_rows;++i) // loop through original basis states
   {
      if ( C(i,i) < 0 )
      {
         C.col(i) *= -1;
      }
   }

}



// where \f[ D(J)_{abcd} \equiv \sqrt{ \frac{1+\delta_{cd}} {1+\delta_{ab}} } \left( C_{ac} C_{bd} -(-1)^{j_c+j_d-J}C_{ad}C_{bc} \right) \f]
//**************************************************************************
/// Takes in an operator expressed in the basis of the original Hamiltonian,
/// and returns that operator in the Hartree-Fock basis.
/// \f[ t_{HF} = C^{-1} t_{HO} C \f]
/// \f[ V_{HF}^{J} = D^{-1} V^{J}_{HO} D \f]
/// \f[ \langle \alpha \beta \gamma | V^{(3)}_{HF} | \delta \epsilon \phi \rangle_{J_{\alpha\beta}J_{\delta\epsilon}J}
///    = \sum_{abcdef} \langle \alpha\beta\gamma | abc \rangle \langle \delta \epsilon \phi | def \rangle
///         \langle abc | V^{(3)}_{HO} | def \rangle_{J_{\alpha\beta}J_{\delta\epsilon}J} \f]
/// The matrix \f$ D \f$ is defined as
/// \f[ D_{ab\alpha\beta} \equiv \sqrt{ \frac{1+\delta_{ab}} {1+\delta_{\alpha\beta}} }  C_{a\alpha} C_{b\beta} \f]
/// The factor in the square root is due to the fact that we're using normalized TBME's.
/// Since only kets with \f$ a\leq b\f$ are stored, we can use the antisymmetry of the TBME's and define
/// \f[ D(J)_{ab\alpha\beta} \equiv \sqrt{ \frac{1+\delta_{ab}} {1+\delta_{\alpha\beta}} }
///      \left( C_{a\alpha} C_{b\beta} -(1-\delta_{ab})(-1)^{j_a+j_b-J} C_{b\alpha}C_{a\beta}\right) \f]
///
//**************************************************************************
Operator HartreeFock::TransformToHFBasis( Operator& OpIn)
{
   Operator OpHF = OpIn;

   if (OpIn.rank_J + OpIn.rank_T + OpIn.parity > 0)
   {
      cout << "Warning, HF transformation hasn't been checked for tensor operators" << endl;
   }

   cout << "Transform one body" << endl;
   // Easy part:
   //Update the one-body part by multiplying by the matrix C(i,a) = <i|a>
   // where |i> is the original basis and |a> is the HF basis
   OpHF.OneBody = C.t() * OpIn.OneBody * C;

   cout << "Transform two body" << endl;
   // Medium-difficulty part:
   //Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz. Most of the effort here is in constructing D.
   int nchan = OpIn.GetModelSpace()->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = OpIn.GetModelSpace()->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      if (npq<1) continue;

      arma::mat D     = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>

      // loop over all possible original basis configurations <pq| in this J,p,Tz channel.
      // and all possible HF configurations |p'q'> in this J,p,Tz channel                                    
      // bra is in the original basis, ket is in the HF basis                                              
      // i and j are the indices of the matrix D for this channel                    
      for (int i=0; i<npq; ++i)    
      {
         Ket & bra = tbc.GetKet(i);   
         for (int j=0; j<npq; ++j)    
         {
            Ket & ket = tbc.GetKet(j); 
            D(i,j) = C(bra.p,ket.p) * C(bra.q,ket.q);
            if (bra.p!=bra.q)
            {
               D(i,j) += C(bra.q,ket.p) * C(bra.p,ket.q) * bra.Phase(tbc.J);
            }
            if (bra.p==bra.q)    D(i,j) *= SQRT2;
            if (ket.p==ket.q)    D(i,j) /= SQRT2;
         }
      }

     auto& IN  =  OpIn.TwoBody[ch].at(ch);
     auto& OUT =  OpHF.TwoBody[ch].at(ch);
     OUT  =    D.t() * IN * D;
   }

   if (OpIn.GetParticleRank() < 3) return OpHF;

   ModelSpace *modelspace = OpIn.GetModelSpace();
   int norbits = modelspace->GetNumberOrbits();
   cout << "Transform three body" << endl;
   for (auto& it_ORB : OpHF.ThreeBody) // loop through orbit labels in the HF basis
   {
      // Double indices are in the transformed basis,
      // equivalent to greek letters in the notes.
      // Single indices are in the original basis
      int aa,bb,cc,dd,ee,ff;
      int a,b,c,d,e,f;
      OpHF.GetOrbitsFromThreeBodyIndex(it_ORB.first, aa,bb,cc,dd,ee,ff);
         Orbit& oa = modelspace->GetOrbit(aa);
         Orbit& ob = modelspace->GetOrbit(bb);
         Orbit& oc = modelspace->GetOrbit(cc);
         Orbit& od = modelspace->GetOrbit(dd);
         Orbit& oe = modelspace->GetOrbit(ee);
         Orbit& of = modelspace->GetOrbit(ff);
//        int Jmin = max(1, max(oa.j2-ob.j2-oc.j2, od.j2-oe.j2-of.j2) )/2;
//        int Jmax = min(oa.j2+ob.j2+oc.j2, od.j2+oe.j2+of.j2)/2;
        int Jmin = max(1, max(oa.j2-ob.j2-oc.j2, od.j2-oe.j2-of.j2) );
        int Jmax = min(oa.j2+ob.j2+oc.j2, od.j2+oe.j2+of.j2);


      for( int a=0;a<norbits;++a)
      {      
         if (C(a,aa) == 0) continue;
      for( int b=0;b<norbits;++b)
      {      
         if (C(b,bb) == 0) continue;
      for( int c=0;c<norbits;++c)
      {      
         if (C(c,cc) == 0) continue;
      for( int d=0;d<norbits;++d)
      {      
         if (C(d,dd) == 0) continue;
      for( int e=0;e<norbits;++e)
      {      
         if (C(e,ee) == 0) continue;
      for( int f=0;f<norbits;++f)
      {      
         if (C(f,ff) == 0) continue;
         double coeff = C(a,aa)*C(b,bb)*C(c,cc)*C(d,dd)*C(e,ee)*C(f,ff);

          for (int J=Jmin; J<=Jmax; J+=2)
          {
             int Jindex = (J-Jmin)/2;
      int Jabmin = max(abs(oa.j2-ob.j2),abs(J-oc.j2))/2;
      int Jabmax = min(oa.j2+ob.j2,J+oc.j2)/2;
      int Jdemin = max(abs(od.j2-oe.j2),abs(J-of.j2))/2;
      int Jdemax = min(od.j2+oe.j2,J+of.j2)/2;
        for (int Jab=Jabmin; Jab<=Jabmax; ++Jab)
        {
        for (int Jde=Jdemin; Jde<=Jdemax; ++Jde)
        {
          int J2index = (Jabmax-Jabmin+1)*(Jde-Jdemin) + (Jab-Jabmin);
          for (int tab=0;tab<=1;++tab)
          {
          for (int tde=0;tde<=1;++tde)
          {
           int Tmax = 2*min(tab,tde)+1;
           for (int T=1;T<=Tmax; T+=2)
           {
            int JTindex = J2index*5 + 2*tab + tde + (T-1)/2;
             cout << "first index: " << Jindex << " (" << it_ORB.second.size() << ")"<< endl;
             cout << "second index: " << JTindex << " (" << it_ORB.second.at(Jindex).size() << ")" << endl;
             cout << "J2index = " << J2index  << "  T= " << T << endl;
             (it_ORB).second.at(Jindex).at(JTindex) += OpIn.GetThreeBodyME(Jab,Jde,J,tab,tde,T,a,b,c,d,e,f) * coeff;
          }
          }
          }
          }
        }
        
//      for (auto& it_orb : OpIn.ThreeBody) // loop through orbit labels in the original basis
//      {
//         OpHF.GetOrbitsFromThreeBodyIndex(it_orb.first, a,b,c,d,e,f);
//         Orbit& oa = modelspace->GetOrbit(a);
//         Orbit& ob = modelspace->GetOrbit(b);
//         Orbit& od = modelspace->GetOrbit(d);
//         Orbit& oe = modelspace->GetOrbit(e);
//         double coeff = C(a,aa)*C(b,bb)*C(c,cc)*C(d,dd)*C(e,ee)*C(f,ff);
//         if (abs(coeff) < 1e-6) continue;
//         double Cabc = C(a,aa)*C(b,bb)*C(c,cc);
//         double Cdef = C(d,dd)*C(e,ee)*C(f,ff);
//         double Cbac = a==b ? 0 : C(b,aa)*C(a,bb)*C(c,cc);
//         double Cedf = d==e ? 0 : C(e,dd)*C(d,ee)*C(f,ff);

//         if ( (abs(Cabc)+abs(Cdef)+abs(Cbac)+abs(Cedf)) < 1e-6) continue;

         // loop over all the J,T couplings. The variable names here are terrible. Sorry.
//         for( auto it_J=it_ORB.second.begin(), it_j=it_orb.second.begin(); it_J!=it_ORB.second.end() and it_j!=it_orb.second.end(); ++it_J, ++it_j)
//         {
//            for ( auto it_JT=it_J->begin(), it_jt=it_j->begin(); it_JT!=it_J->end() and it_jt!=it_j->end(); ++it_JT, ++it_jt)
//            {
//         for( auto it_J=it_ORB.second.begin(), it_j=it_orb.second.begin(); it_J!=it_ORB.second.end() and it_j!=it_orb.second.end(); ++it_J, ++it_j)
//         {
//            for ( auto it_JT=it_J->begin(), it_jt=it_j->begin(); it_JT!=it_J->end() and it_jt!=it_j->end(); ++it_JT, ++it_jt)
//            {
//               (*it_JT) += (*it_jt) * coeff;
//               double coeff = (Cabc + modelspace->phase(oa.j2+ob.j2+1+Jab+tab)*Cbac);
//               double coeff *= (Cdef + modelspace->phase(od.j2+oe.j2+1+Jde+tde)*Cefd);
//               (*it_JT) += (*it_jt) * coeff;
//            }
//         }

      }
      }
      }
      }
      }
      }
      }

   }




   return OpHF;
}






