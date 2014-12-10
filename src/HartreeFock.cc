
#include "HartreeFock.hh"
#include "ModelSpace.hh"

using namespace std;

//HartreeFock::HartreeFock(Operator *hbare, Operator *hhf)
//HartreeFock::HartreeFock(Operator *hbare=NULL)
HartreeFock::HartreeFock(Operator& hbare) : Hbare(hbare)
{
//   if (hbare ==NULL) cout << "Ah! Null operator" << endl;
//   Hbare = hbare;
   tolerance = 1e-10;
   ediff = 1.0;
   //ModelSpace * ms = Hbare->GetModelSpace();
   ModelSpace * ms = Hbare.GetModelSpace();
   int noits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();

   C    = arma::mat(noits,noits,arma::fill::eye);
   Vab  = arma::mat(noits,noits);
   H    = arma::mat(noits,noits);
   Vmon = arma::mat(2*nKets,2*nKets);
   //Vmon = arma::mat(nKets,nKets);
   prev_energies = arma::vec(noits,arma::fill::zeros);

   //t = Hbare->OneBody;
   t = Hbare.OneBody;
   energies = t.diag();
   UpdateDensityMatrix();
   BuildMonopoleV();
}



void HartreeFock::Solve()
{
//   int ncore = Hbare->GetModelSpace()->nCore;
   //int noits = Hbare->GetModelSpace()->GetNumberOrbits();
   int noits = Hbare.GetModelSpace()->GetNumberOrbits();
   int iter = 0; // counter so we don't go on forever
   //ModelSpace *ms = Hbare->GetModelSpace();
   ModelSpace *ms = Hbare.GetModelSpace();

/*
   C(0,0) = sqrt(2./3);
   C(0,10) = sqrt(1./3);
   C(10,0) = -sqrt(1./3);
   C(10,10) = sqrt(2./3);
*/

//   C.swap_cols(0,10);
//   C.swap_cols(1,11);
   UpdateH();
//   return;
//   cout << "Input t" << endl;
//   t.print();
//   cout << "Input Vab" << endl;
//   Vab.print();
//   cout << "Input H" << endl;
//   H.print();
/*
   cout << "Input Vmon" << endl;
   Vmon.print();
   cout << "Input CHC" << endl;
   (C.t() * H * C).print();
   cout << "Input CtC" << endl;
   (C.t() * t * C).print();
   cout << "Input CVabC" << endl;
   (C.t() * Vab * C).print();
*/
   while( (! CheckConvergence()) && iter < 10000)
   {
   
//      UpdateHFOrbits();
      UpdateDensityMatrix();
      UpdateH();
/*
   cout << "-----------------------------------------------------------------" << endl;
   cout << "Input H " << iter << endl;
   H.print();
   cout << "Input t" << endl;
   t.print();
   cout << "Input Vab" << endl;
   Vab.print();
*/

      Diagonalize2();
//      cout << "Energies:" << endl;
//      energies.print();
//      cout << "Eigenvectors:" << endl;
//      C.print();

//      cout << iter << ":  " << ediff << endl;
      iter++;
   }
//   cout << "output H" << endl;
//   H.print();
//   cout << "output t" << endl;
//   t.print();
//   cout << "Output Vab" << endl;
//   Vab.print();
//   cout << "Cia:" << endl;
//   C.print();
/*
   cout << "Output CHC" << endl;
   (C.t() * H * C).print();
   cout << "Output CtC" << endl;
   (C.t() * t * C).print();
   cout << "Output CVabC" << endl;
   (C.t() * Vab * C).print();
*/
   if (iter==10000)
   {
      cerr << "HartreeFock calculation didn't converge after 10000 iterations" << endl;
   }
   // Armadillo returns the oits in order of increasing eigenvalue.
   // Reorder it so that it coincides with the original modelspace ordering.
   ReorderCoefficients();
   // Calculate the HF energy.
   CalcEHF();
/*
    for (int a=0;a<noits;a++)
  {
      Orbit & oa = ms->GetOrbit(a);
      for (int b=0;b<noits;b++)
      {
         Orbit & ob = ms->GetOrbit(b);
         if (oa.j2 != ob.j2 or oa.l != ob.l or oa.tz2 != ob.tz2) continue;
         
         cout << "H("<<a<<","<<b<<") = t(a,b) + V(a,b) = " << t(a,b) << " + " << Vab(a,b) << " = " << H(a,b) << endl;
      }
  }
*/
//   UpdateHFOrbits();

//   cout << "Transformed: " << endl;
//   cout << "T: " << endl;
//   ( C.t() * t * C ).print();
//   cout << "Vab: " << endl;
//   ( C.t() * Vab * C ).print();
//   cout << "H: " << endl;
//   ( C.t() * H * C ).print();

//  cout << "Energies: " << endl;
//  energies.print();
}

void HartreeFock::CalcEHF()
{
   //ModelSpace * ms = Hbare->GetModelSpace();
   ModelSpace * ms = Hbare.GetModelSpace();
   EHF = 0;
   //int noits = Hbare->GetModelSpace()->GetNumberOrbits();
   int noits = Hbare.GetModelSpace()->GetNumberOrbits();
   for (int i=0;i<noits;i++)
   {
      for (int j=0;j<noits;j++)
      {
         EHF += 0.5 * rho(i,j) * (ms->GetOrbit(i).j2+1) * (H(i,j)+t(i,j));
      }
   }

}


///************************************************
// Diagonalize() -- [See Suhonen eq. 4.85]
// Diagonalize <a|H|b> and put the
// eigenvectors in C(i,alpha) = <alpha|i>
// and eigenvalues in the vector energies.
// Save the last vector of energies to check
// for convergence.
///************************************************
void HartreeFock::Diagonalize()
{
   prev_energies = energies;
//   bool success = arma::eigs_sym(energies, C, H);
/*   int noits = Hbare->GetModelSpace()->GetNumberOrbits();
   for (int a=0;a<noits;a++)
   {
      float normalization = arma::norm(C.col(a));
      C.col(a) /= normalization;
   }
*/
}

///***********************************************************************************
// Try alternate strategy, where different channels are diagonalized independently.
// This guarantees that J,Tz, and parity remain good. In the first method,
// it seems that numerical instabilities can lead to mixing of those quantum numbers.
// Also, this way is faster, if less elegant.
///***********************************************************************************
void HartreeFock::Diagonalize2()
{
   prev_energies = energies;
   vector<int> oit_list;
   arma::mat H_ch;
   arma::mat C_ch;
   arma::vec E_ch;
   bool success;
   //int noits = Hbare->GetModelSpace()->GetNumberOrbits();
   int noits = Hbare.GetModelSpace()->GetNumberOrbits();
   for (int p = 0; p<=1;p++)
   {
      for (int Tz = -1; Tz<=1; Tz+=2)
      {
          for (int J=0;J<JMAX;J++)
          {

            // Find all the SP oits that have quantum numbers J,p,Tz
            // and store them in a list
             oit_list.resize(0);
             for (int a=0;a<noits;a++)
             {
                //Orbit &oa = Hbare->GetModelSpace()->GetOrbit(a);
                Orbit &oa = Hbare.GetModelSpace()->GetOrbit(a);
                if (oa.j2==J and oa.tz2==Tz and (oa.l%2)==p)
                {
                   oit_list.push_back(a);
                }
             }
             // Now create submatrices corresponding to just these oits
             int norb = oit_list.size();
             if (norb < 1) continue;
             H_ch = arma::mat(norb,norb,arma::fill::zeros);
             C_ch = arma::mat(norb,norb,arma::fill::zeros);
             E_ch = arma::vec(norb,arma::fill::zeros);
             for (int a=0;a<norb;a++)
             {
                for (int b=0;b<norb;b++)
                {
                   H_ch(a,b) = H(oit_list[a],oit_list[b]);
                }
             }
             // Diagonalize the submatrix
             success = arma::eig_sym(E_ch, C_ch, H_ch);
             // Update the full overlap matrix C and energy vector
             for (int a=0;a<norb;a++)
             {
                energies( oit_list[a] ) = E_ch(a);
                for (int b=0;b<norb;b++)
                {
                   C(oit_list[a],oit_list[b]) = C_ch(a,b);
                }
             }

          } // for J...
      } // For Tz ...
   } // For p...
}


//**************************************************
// BuildMonopoleV()
// Construct the unnormalized monople Hamiltonian
// <ab|V_mon|cd> = Sum_J (2J+1) <ab|V|cd>_J.
//   To facilitate looping, the matrix has 4 blocks.
//   If a<=b and c<=d, the matrix looks like:
//
//         [ <ab|V|cd> ... <ab|V|dc> ]
//    V  = |     :             :     |
//         [ <ba|V|cd> ... <ba|V|dc> ]
//         
//**************************************************

/*
void HartreeFock::BuildMonopoleV()
{
   Vmon.zeros();
   int nKets = Hbare->GetModelSpace()->GetNumberKets();
   int nchan = Hbare->GetModelSpace()->GetNumberTwoBodyChannels();
   
   for (int ch=0;ch<nchan;++ch) // loop over J,p,Tz channels
   {
      //TwoBodyChannel *tbc = Hbare->GetTwoBodyChannel(i);
      TwoBodyChannel& tbc = Hbare->GetModelSpace()->GetTwoBodyChannel(ch);
//      cout << "J,p,t = " << tbc.J << " " << tbc.parity << " " << tbc.Tz << endl;
//      Hbare->TwoBody[ch].print();
      int J = tbc.J;
      int npq = tbc.GetNumberKets();
      for (int a=0;a<npq;++a)
      {
         Ket & bra = tbc.GetKet(a);
         int ibra = tbc.GetKetIndex(a);
         for (int b=0;b<npq;++b)
         {
            Ket & ket  = tbc.GetKet(b);
            int iket = tbc.GetKetIndex(b);
            //double tbme = (2*J+1)*tbc.GetTBME(bra,ket);
//            if (a==0 and b==5)
//              cout << "J = " << J << ". Adding " << Hbare->TwoBody[ch](a,b) << endl;
            //double tbme = (2*J+1)* Hbare->TwoBody[ch](a,b);
            double tbme = (2*J+1)* sqrt((1+bra.delta_pq())*(1+ket.delta_pq())) *Hbare->TwoBody[ch](a,b);
            Vmon(ibra,iket) += tbme;
            Vmon(ibra+nKets,iket) += bra.Phase(J)*tbme;
            Vmon(ibra,iket+nKets) += ket.Phase(J)*tbme;
            Vmon(ibra+nKets,iket+nKets) += bra.Phase(J)*ket.Phase(J)*tbme;
         }
      }
   }
}
*/

//**************************************************
// BuildMonopoleV()
// Construct the unnormalized monople Hamiltonian
// <ab|V_mon|cd> = Sum_J (2J+1) <ab|V|cd>_J.
//  I calculate and store each term, rather than
//  repeatedly calculating things.
//   To facilitate looping, the matrix has 4 blocks.
//   If a<=b and c<=d, the matrix looks like:
//
//         [ <ab|V|cd> ... <ab|V|dc> ]
//    V  = |     :             :     |
//         [ <ba|V|cd> ... <ba|V|dc> ]
//         
//**************************************************
void HartreeFock::BuildMonopoleV()
{
   Vmon.zeros();
   int nKets = Hbare.GetModelSpace()->GetNumberKets();
   for (int ibra=0;ibra<nKets;++ibra)
   {
      Ket & bra = Hbare.GetModelSpace()->GetKet(ibra);
      int a = bra.p;
      int b = bra.q;
      Orbit & oa = Hbare.GetModelSpace()->GetOrbit(a);
      Orbit & ob = Hbare.GetModelSpace()->GetOrbit(b);
      double norm = (oa.j2+1)*(ob.j2+1);
      for (int iket=0;iket<nKets;++iket)
      {
         Ket & ket = Hbare.GetModelSpace()->GetKet(iket);
         int c = ket.p;
         int d = ket.q;
         Vmon(ibra,iket)             = Hbare.GetTBMEmonopole(a,b,c,d) * norm;
         Vmon(ibra+nKets,iket)       = Hbare.GetTBMEmonopole(b,a,c,d) * norm;
         Vmon(ibra,iket+nKets)       = Hbare.GetTBMEmonopole(a,b,d,c) * norm;
         Vmon(ibra+nKets,iket+nKets) = Hbare.GetTBMEmonopole(b,a,d,c) * norm;
      }
   }
   
}



//**************************************************************************
// 1-body density matrix 
// <i|rho|j> = Sum_beta <i|beta> <j|beta>
// where beta runs over HF oits in the core (i.e. below the fermi surface)
//**************************************************************************
void HartreeFock::UpdateDensityMatrix()
{
   //ModelSpace * ms = Hbare->GetModelSpace();
   ModelSpace * ms = Hbare.GetModelSpace();
   int noits = ms->GetNumberOrbits();

   // Pcore is a projector onto core oits.
   arma::mat Pcore = arma::mat(noits,noits,arma::fill::zeros);
   for (int &beta : ms->holes)
   {
       Pcore(beta,beta) = 1;
   }

   rho = C * Pcore * C.t();
}


//*********************************************************************
// UpdateH() -- [See Suhonen eq 4.85]
// <a|H|b> = <a|t|b>  +  Sum_ij <i|rho|j> <ai|V_mon|bj>
// * H is the fock matrix, to be diagonalized
// * t is the kinetic energy
// * rho is the density matrix defined in UpdateDensityMatrix()
// * V_mon is the monopole component of the 2-body interaction.
//*********************************************************************
void HartreeFock::UpdateH()
{
   //ModelSpace * ms = Hbare->GetModelSpace();
   ModelSpace * ms = Hbare.GetModelSpace();
   int noits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();
   int bra, ket;
   Vab.zeros();
   for (int a=0;a<noits;a++)
   {
      for (int b=a;b<noits;b++)
      {
         if (ms->GetOrbit(a).j2  != ms->GetOrbit(b).j2)  continue;
         if (ms->GetOrbit(a).tz2 != ms->GetOrbit(b).tz2) continue;
         if (ms->GetOrbit(a).l   != ms->GetOrbit(b).l)   continue;
         for (int i=0;i<noits;i++)
         {
            // The monopoles are listed for |ab> fist with a<=b, then for a>=b
            // so if a>b, add nKets.
            bra = ms->GetKetIndex(min(a,i),max(a,i));
//            bra = ms->GetKetIndex(a,i);
            if (a>i) bra += nKets;
            for (int j=0;j<noits;j++)
            {
               ket = ms->GetKetIndex(min(b,j),max(b,j));
//               ket = ms->GetKetIndex(b,j);
               if (b>j) ket += nKets;

               Vab(a,b) += rho(i,j)*Vmon(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>
               //Vab(a,b) += rho(i,j)*Vmon(bra,ket); // <i|rho|j> * <ai|Vmon|bj>

           }
         }
         Vab(a,b) /= (ms->GetOrbit(a).j2+1); // divide by 2ja+1
         Vab(b,a) = Vab(a,b);  // Hermitian & real => symmetric
      }
   }
   H = t + Vab;

}


bool HartreeFock::CheckConvergence()
{
   arma::vec de = energies - prev_energies;
   ediff = sqrt(arma::dot(de,de)) / energies.size();
//   cout << "ediff = " << ediff << endl;
   return (ediff < tolerance);
}




// Some more thought should go into how to do this properly.
void HartreeFock::UpdateHFOrbits()
{
   //ModelSpace * ms = Hbare->GetModelSpace();
   ModelSpace * ms = Hbare.GetModelSpace();
   int noits = ms->GetNumberOrbits();
   for (int a=0;a<noits;a++) // loop over HF basis
   {
      float L = 0;
      float J = 0;
      float N = 0;
      float Tz = 0;
      for (int i=0;i<noits;i++) // loop over original basis
      {
         Orbit & oi = ms->GetOrbit(i);
         float C2 = C(i,a)*C(i,a);  // Overlap |<a|i>|^2
         L  += C2*oi.l;
         J  += C2*oi.j2;
         N  += C2*oi.n;
         Tz += C2*oi.tz2;
      }
      int indx = ms->Index1(round(N), round(L), round(J), round(Tz));
      int ph = ms->GetOrbit(indx).ph;
      int io = ms->GetOrbit(indx).io;
      Orbit &oa = ms->GetOrbit(a);
      //oa.Set(round(N),round(L),round(J),round(Tz),Hvq,energies[a]);
      oa.Set(round(N),round(L),round(J),round(Tz),ph,io,energies[a]);
   }
}



void HartreeFock::PrintOrbits()
{
  //ModelSpace * ms = Hbare->GetModelSpace();
  ModelSpace * ms = Hbare.GetModelSpace();
  int noits = ms->GetNumberOrbits();
  for (int a=0;a<noits;a++)
  {
     Orbit & oa = ms->GetOrbit(a);
     cout << a << ": E= " << oa.spe << " L=" << oa.l
          << " J=" << oa.j2 << " N=" << oa.n << " Tz=" << oa.tz2
          << " ph=" << oa.ph << " io=" << oa.io << endl;
  }

}

void HartreeFock::ReorderCoefficients()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int noits = ms->GetNumberOrbits();
   arma::mat C_tmp = C;
   arma::vec e_tmp = energies;
   for (int i=0;i<noits;i++)
   {
      float fmax = 0.0;
      int kmax;
      int sign = 1;
      for (int k=0;k<noits;k++)
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


Operator HartreeFock::TransformToHFBasis( Operator& OpIn)
{
   Operator OpHF = OpIn;

   //Update the one-body part by multiplying by the matrix C(i,a) = <i|a>.
   OpHF.OneBody = C.t() * OpIn.OneBody * C;

   //Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz
   int nchan = OpIn.GetModelSpace()->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = OpIn.GetModelSpace()->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      if (npq<1) continue;
      arma::mat D     = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
      arma::mat Dexch = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ba> = <ji|ab>
      for (int i=0; i<npq; i++)       // loop over all possible original basis configurations <pq| in this J,p,Tz channel
      {                               // i and j are the indices of the small matrix for this channel
         Ket & bra = tbc.GetKet(i); // bra is in the original basis
         for (int j=0; j<npq; j++)       // loop over all possible HF configurations |p'q'> in this J,p,Tz channel
         {
            Ket & ket = tbc.GetKet(j); // ket is in the HF basis
            double normfactor = sqrt((1.0+ket.delta_pq())/(1.0+bra.delta_pq()));
            D(i,j) = C(bra.p,ket.p) * C(bra.q,ket.q) * normfactor;
            Dexch(i,j) = C(bra.p,ket.q) * C(bra.q,ket.p) * ket.Phase(tbc.J) * (1-ket.delta_pq()) * normfactor;
         }
      }

     OpHF.TwoBody[ch]   = D.t()      * OpIn.TwoBody[ch] * D;
     OpHF.TwoBody[ch]  += Dexch.t()  * OpIn.TwoBody[ch] * D;
     OpHF.TwoBody[ch]  += D.t()      * OpIn.TwoBody[ch] * Dexch;
     OpHF.TwoBody[ch]  += Dexch.t()  * OpIn.TwoBody[ch] * Dexch;

   }
   return OpHF;
}

