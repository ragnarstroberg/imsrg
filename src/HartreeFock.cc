
#include "HartreeFock.hh"
#include "ModelSpace.hh"

using namespace std;

//HartreeFock::HartreeFock(Operator *hbare, Operator *hhf)
HartreeFock::HartreeFock(Operator *hbare)
{
   Hbare = hbare;
   tolerance = 1e-10;
   ediff = 1.0;
   ModelSpace * ms = Hbare->GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();

   C    = arma::mat(norbits,norbits,arma::fill::eye);
   Vab  = arma::mat(norbits,norbits);
   H    = arma::mat(norbits,norbits);
   Vmon = arma::mat(2*nKets,2*nKets);
   prev_energies = arma::vec(norbits,arma::fill::zeros);

   t = Hbare->OneBody;
   energies = t.diag();
   UpdateDensityMatrix();
   BuildMonopoleV();
}



void HartreeFock::Solve()
{
//   int ncore = Hbare->GetModelSpace()->nCore;
   int norbits = Hbare->GetModelSpace()->GetNumberOrbits();
   int iter = 0; // counter so we don't go on forever
   ModelSpace *ms = Hbare->GetModelSpace();

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
/*
   cout << "Input H" << endl;
   H.print();
   cout << "Input t" << endl;
   t.print();
   cout << "Input Vab" << endl;
   Vab.print();
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
//      for (int a=0;a<norbits;a++)
//      {
//         if (Hbare->GetModelSpace()->GetOrbit(a)->hvq>0) continue;
//         cout << Hbare->GetModelSpace()->GetOrbit(a)->spe << " " ;
//      }
//      cout << endl;

      Diagonalize2();
//      cout << "Energies:" << endl;
//      energies.print();
//      cout << "Eigenvectors:" << endl;
//      C.print();

//      cout << iter << ":  " << ediff << endl;
      iter++;
   }
/*
   cout << "output H" << endl;
   H.print();
   cout << "output t" << endl;
   t.print();
   cout << "Output Vab" << endl;
   Vab.print();
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
   // Armadillo returns the orbits in order of increasing eigenvalue.
   // Reorder it so that it coincides with the original modelspace ordering.
   ReorderCoefficients();
   // Calculate the HF energy.
   CalcEHF();
/*
    for (int a=0;a<norbits;a++)
  {
      Orbit* oa = ms->GetOrbit(a);
      for (int b=0;b<norbits;b++)
      {
         Orbit* ob = ms->GetOrbit(b);
         if (oa->j2 != ob->j2 or oa->l != ob->l or oa->tz2 != ob->tz2) continue;
         
         cout << "H("<<a<<","<<b<<") = t(a,b) + V(a,b) = " << t(a,b) << " + " << Vab(a,b) << " = " << H(a,b) << endl;
      }
  }
*/
//   UpdateHFOrbits();
   
}

void HartreeFock::CalcEHF()
{
   ModelSpace * ms = Hbare->GetModelSpace();
   EHF = 0;
   int norbits = Hbare->GetModelSpace()->GetNumberOrbits();
   int nKets = Hbare->GetModelSpace()->GetNumberKets();
   for (int i=0;i<norbits;i++)
   {
      for (int j=0;j<norbits;j++)
      {
         EHF += 0.5 * rho(i,j) * (ms->GetOrbit(i)->j2+1) * (H(i,j)+t(i,j));
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
/*   int norbits = Hbare->GetModelSpace()->GetNumberOrbits();
   for (int a=0;a<norbits;a++)
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
   vector<int> orbit_list;
   arma::mat H_ch;
   arma::mat C_ch;
   arma::vec E_ch;
   bool success;
   int norbits = Hbare->GetModelSpace()->GetNumberOrbits();
   for (int p = 0; p<=1;p++)
   {
      for (int Tz = -1; Tz<=1; Tz+=2)
      {
          for (int J=0;J<JMAX;J++)
          {

            // Find all the SP orbits that have quantum numbers J,p,Tz
            // and store them in a list
             orbit_list.resize(0);
             for (int a=0;a<norbits;a++)
             {
                Orbit *orba = Hbare->GetModelSpace()->GetOrbit(a);
                if (orba->j2==J and orba->tz2==Tz and (orba->l%2)==p)
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
             success = arma::eig_sym(E_ch, C_ch, H_ch);
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

void HartreeFock::BuildMonopoleV()
{
   int nKets = Hbare->GetModelSpace()->GetNumberKets();
   Vmon.zeros();
   int nchan = Hbare->GetNumberTwoBodyChannels();
   
   for (int i=0;i<nchan;++i) // loop over J,p,Tz channels
   {
      TwoBodyChannel *tbc = Hbare->GetTwoBodyChannel(i);
      int J = tbc->J;
      int npq = tbc->GetNumberKets();
      for (int a=0;a<npq;++a)
      {
         Ket * bra = tbc->GetKet(a);
         int ibra = tbc->GetKetIndex(a);
         for (int b=0;b<npq;++b)
         {
            Ket * ket  = tbc->GetKet(b);
            int iket = tbc->GetKetIndex(b);
            Vmon(ibra,iket) += (2*J+1)*tbc->GetTBME(bra,ket);
            Vmon(ibra+nKets,iket) += (2*J+1)*bra->Phase(J)*tbc->GetTBME(bra,ket);
            Vmon(ibra,iket+nKets) += (2*J+1)*ket->Phase(J)*tbc->GetTBME(bra,ket);
            Vmon(ibra+nKets,iket+nKets) += (2*J+1)*bra->Phase(J)*ket->Phase(J)*tbc->GetTBME(bra,ket);
         }
      }
   }
}



//**************************************************************************
// 1-body density matrix 
// <i|rho|j> = Sum_beta <i|beta> <j|beta>
// where beta runs over HF orbits in the core (i.e. below the fermi surface)
//**************************************************************************
void HartreeFock::UpdateDensityMatrix()
{
   ModelSpace * ms = Hbare->GetModelSpace();
   int norbits = ms->GetNumberOrbits();

   // Pcore is a projector onto core orbits.
   arma::mat Pcore = arma::mat(norbits,norbits,arma::fill::eye);
   for (int beta=0;beta<norbits;beta++)
   {
      if (ms->GetOrbit(beta)->hvq > 0) // hvq==0 means hole state
      {
         Pcore(beta,beta) = 0;
      }
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
   ModelSpace * ms = Hbare->GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();
   int bra, ket;
//   Vab  = arma::mat(norbits,norbits,arma::fill::zeros);
   Vab.zeros();
   for (int a=0;a<norbits;a++)
   {
      for (int b=a;b<norbits;b++)
      {
         if (ms->GetOrbit(a)->j2  != ms->GetOrbit(b)->j2)  continue;
         if (ms->GetOrbit(a)->tz2 != ms->GetOrbit(b)->tz2) continue;
         if (ms->GetOrbit(a)->l   != ms->GetOrbit(b)->l)   continue;
         for (int i=0;i<norbits;i++)
         {
            // The monopoles are listed for |ab> fist with a<=b, then for a>=b
            // so if a>b, add nKets.
            bra = ms->GetKetIndex(min(a,i),max(a,i));
            if (a>i) bra += nKets;
            for (int j=0;j<norbits;j++)
            {
               ket = ms->GetKetIndex(min(b,j),max(b,j));
               if (b>j) ket += nKets;

               Vab(a,b) += rho(i,j)*Vmon(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>

           }
         }
         Vab(a,b) /= (ms->GetOrbit(a)->j2+1); // divide by 2ja+1
         Vab(b,a) = Vab(a,b);  // Hermitian & real => symmetric
      }
   }
   H = t + Vab;

}


bool HartreeFock::CheckConvergence()
{
   arma::vec de = energies - prev_energies;
   ediff = sqrt(arma::dot(de,de)) / energies.size();
   cout << "ediff = " << ediff << endl;
   return (ediff < tolerance);
}




// Some more thought should go into how to do this properly.
void HartreeFock::UpdateHFOrbits()
{
   ModelSpace * ms = Hbare->GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   for (int a=0;a<norbits;a++) // loop over HF basis
   {
      float L = 0;
      float J = 0;
      float N = 0;
      float Tz = 0;
      for (int i=0;i<norbits;i++) // loop over original basis
      {
         Orbit * orbi = ms->GetOrbit(i);
         float C2 = C(i,a)*C(i,a);  // Overlap |<a|i>|^2
         L  += C2*orbi->l;
         J  += C2*orbi->j2;
         N  += C2*orbi->n;
         Tz += C2*orbi->tz2;
      }
      int indx = ms->Index1(round(N), round(L), round(J), round(Tz));
      int Hvq = ms->GetOrbit(indx)->hvq;
      Orbit *orba = Hbare->GetModelSpace()->GetOrbit(a);
      orba->Set(round(N),round(L),round(J),round(Tz),Hvq,energies[a]);
   }
}



void HartreeFock::PrintOrbits()
{
  ModelSpace * ms = Hbare->GetModelSpace();
  int norbits = ms->GetNumberOrbits();
  for (int a=0;a<norbits;a++)
  {
     Orbit * orba = ms->GetOrbit(a);
     cout << a << ": E= " << orba->spe << " L=" << orba->l
          << " J=" << orba->j2 << " N=" << orba->n << " Tz=" << orba->tz2
          << " Hvq=" << orba->hvq << endl;
  }

}


void HartreeFock::ReorderCoefficients()
{
   ModelSpace * ms = Hbare->GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   arma::mat C_tmp = C;
   arma::vec e_tmp = energies;
   for (int i=0;i<norbits;i++)
   {
      float fmax = 0.0;
      int kmax;
      for (int k=0;k<norbits;k++)
      {
         if ( abs(C_tmp(i,k)) > fmax )
         {
            fmax = abs(C_tmp(i,k));
            kmax = k;
         }
      }
      C.col(kmax) = C_tmp.col(i);
      energies(kmax) = e_tmp(i);
   }
}


void HartreeFock::TransformToHFBasis(Operator *OpIn, Operator *OpHF)
{
   // Maybe put in a check that they should have the same modelspace?

   cout << "Fock matrix:" << endl;
//   (C.t() * H * C).print();

   ModelSpace * ms = OpIn->GetModelSpace();
   int norbits = ms->GetNumberOrbits();

   //Update the one-body part by multiplying by the matrix C(i,a) = <i|a>.

   OpHF->OneBody = C.t() * OpIn->OneBody * C;

   //Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz

   int nchan = OpIn->GetNumberTwoBodyChannels();
   for (int i=0;i<nchan;++i)
   {
      TwoBodyChannel *chin = OpIn->GetTwoBodyChannel(i);
      int npq = chin->GetNumberKets();
      if (npq<1) continue;
      TwoBodyChannel *chhf = OpHF->GetTwoBodyChannel(i);
      arma::mat D     = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
      arma::mat Dexch = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ba> = <ji|ab>
      for (int i=0; i<npq; i++)       // loop over all possible original basis configurations <pq| in this J,p,Tz channel
      {                               // i and j are the indices of the small matrix for this channel
         Ket * bra = chin->GetKet(i); // bra is in the original basis
         for (int j=0; j<npq; j++)       // loop over all possible HF configurations |pq> in this J,p,Tz channel
         {
            Ket * ket = chin->GetKet(j); // ket is in the HF basis
            D(i,j) = C(bra->p,ket->p) * C(bra->q,ket->q) * 1./(1.0+bra->delta_pq());
            Dexch(i,j) = C(bra->p,ket->q) * C(bra->q,ket->p) * ket->Phase(chin->J) * 1./(1.0+bra->delta_pq());
         }
      }
     chhf->TBME  = D.t()     * chin->TBME * D;
     chhf->TBME += Dexch.t() * chin->TBME * D;
     chhf->TBME += D.t()     * chin->TBME * Dexch;
     chhf->TBME += Dexch.t() * chin->TBME * Dexch;
   }
}

