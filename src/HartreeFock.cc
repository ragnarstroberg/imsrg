
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
   Vmon = arma::mat(2*nKets,2*nKets);
   prev_energies = arma::vec(norbits,arma::fill::zeros);

   t = Hbare->OneBody;
   energies = t.diag();
   UpdateDensityMatrix();
   BuildMonopoleV();
}



void HartreeFock::Solve()
{
   int ncore = Hbare->GetModelSpace()->nCore;
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
// <ab|V_mon|cd> = Sum_J (2J+1) <ab|V|cd>_J
//**************************************************

void HartreeFock::BuildMonopoleV()
{
   ModelSpace * ms = Hbare->GetModelSpace();
   int nKets = ms->GetNumberKets();
   //Vmon = arma::mat(2*nKets,2*nKets,arma::fill::zeros);
//   Vmon = arma::sp_mat(2*nKets,2*nKets);
   cout << "nKets = " << nKets << endl;
   Vmon.zeros();

   int ketcheck = ms->GetKetIndex(0,11);
   int bracheck = ms->GetKetIndex(1,10);
   
   
   for (int a=0;a<nKets;a++)
   {
      Ket* bra = Hbare->GetModelSpace()->GetKet(a%nKets);
      int par = bra->parity;
      int Tz = bra->Tz;
      for (int b=a;b<nKets;b++)
      {
         Ket* ket = Hbare->GetModelSpace()->GetKet(b%nKets);
         if (ket->parity != bra->parity or ket->Tz != bra->Tz) continue; // Interaction conserves T_z, parity

         int jmin = max(bra->Jmin, ket->Jmin);
         int jmax = min(bra->Jmax, ket->Jmax);
//         cout << "jmin = " << jmin << "  jmax = " << jmax << endl;
         int jstep = max(bra->Jstep, ket->Jstep);
         for (int J=jmin;J<=jmax;J+=jstep)
         {
            TwoBodyChannel * channel = Hbare->GetTwoBodyChannel(J,par,Tz);

            Vmon(a,b)             += (2*J+1)*channel->GetTBME(bra,ket);
            Vmon(a+nKets,b)       += (2*J+1)*bra->Phase(J)*channel->GetTBME(bra,ket);
            Vmon(a,b+nKets)       += (2*J+1)*ket->Phase(J)*channel->GetTBME(bra,ket);
            Vmon(a+nKets,b+nKets) += (2*J+1)*bra->Phase(J)*ket->Phase(J)*channel->GetTBME(bra,ket);
            //if (a==ketcheck and b==bracheck)
            if ((a==ketcheck and b==bracheck) or (a==bracheck and b==ketcheck))
            {
               cout << "KetCheck: J = " << J << " par = " << par << " Tz = " << Tz << " TBME = " << channel->GetTBME(bra,ket) << " Vmon(a,b) = " << Vmon(a,b) << " (" << Vmon(b,a) << ")" << endl;
            }
//            cout << "Building monopole: " << a << " " << b << " " << J << endl;
         }
         Vmon(b,a) = Vmon(a,b);
         Vmon(b,a+nKets) = Vmon(a+nKets,b);
         Vmon(b+nKets,a) = Vmon(a,b+nKets);
         Vmon(b+nKets,a+nKets) = Vmon(a+nKets,b+nKets);
      }
//      cout << a << " " << bra->a << " " << bra->b << endl;
   }

//   cout << "Vmon = " << endl;
//   Vmon.print();
//
//   Vmon.submat(0,0,6,6).print();
//   cout << "<0  0| Vmon |0   0> = " << Vmon(ms->GetKetIndex(0,0), ms->GetKetIndex(0,0)) << endl;
//   cout << "<0  0| Vmon |10  0> = " << Vmon(ms->GetKetIndex(0,0), ms->GetKetIndex(0,10)+nKets) << endl;
//   cout << "<0 10| Vmon |10 10> = " << Vmon(ms->GetKetIndex(0,10), ms->GetKetIndex(10,10)) << endl;
//   cout << "<0  0| Vmon |10 10> = " << Vmon(ms->GetKetIndex(0,0), ms->GetKetIndex(10,10)) << endl;
//   cout << "<0 10| Vmon |10  0> = " << Vmon(ms->GetKetIndex(0,0), ms->GetKetIndex(0,10)+nKets) << endl;
//   cout << "<0  1| Vmon |10  1> = " << Vmon(ms->GetKetIndex(0,1), ms->GetKetIndex(1,10)+nKets) << endl;
//   cout << "<0  1| Vmon |10 11> = " << Vmon(ms->GetKetIndex(0,1), ms->GetKetIndex(10,11)) << endl;
//   cout << "<0 11| Vmon |10 11> = " << Vmon(ms->GetKetIndex(0,11), ms->GetKetIndex(10,11)) << endl;
//   cout << "<0 11| Vmon |10  1> = " << Vmon(ms->GetKetIndex(0,11), ms->GetKetIndex(1,10)+nKets) << endl;


//   cout << "<5  1| Vmon |5  1> = " << Vmon(ms->GetKetIndex(1,5)+nKets, ms->GetKetIndex(1,5)+nKets) << endl;
//   cout << "<5 11| Vmon |5 11> = " << Vmon(ms->GetKetIndex(1,10), ms->GetKetIndex(1,10)) << endl;
//   cout << "<5  1| Vmon |5 11> = " << Vmon(ms->GetKetIndex(1,5)+nKets, ms->GetKetIndex(5,11)) << endl;
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

//   cout << "rho: " << endl;
//   rho.print();

}


//*********************************************************************
// UpdateH() -- [See Suhonen eq 4.85]
// <a|H|b> = <a|t|b>  +  Sum_ij <i|rho|j> <ai|V_mon|bj>
// Where H is the 1-body Hamiltonian to be diagonalized
// t is the kinetic energy
// rho is the density matrix defined in UpdateDensityMatrix()
// V_mon is the monopole 2-body interaction.
//*********************************************************************
void HartreeFock::UpdateH()
{
   ModelSpace * ms = Hbare->GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();
   int bra, ket;
   //Vab.zeros();
   Vab  = arma::mat(norbits,norbits,arma::fill::zeros);

   for (int a=0;a<norbits;a++)
   {
      int ja = ms->GetOrbit(a)->j2;
      for (int b=a;b<norbits;b++)
      {
         if (ms->GetOrbit(a)->j2  != ms->GetOrbit(b)->j2) continue;
         if (ms->GetOrbit(a)->tz2 != ms->GetOrbit(b)->tz2) continue;
         if (ms->GetOrbit(a)->l   != ms->GetOrbit(b)->l) continue;
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

               //Vab(a,b) += rho(i,j)*Vmon(bra,ket);
               Vab(a,b) += rho(i,j)*Vmon(min(bra,ket),max(bra,ket));

//               if (a==0 and b==10 and abs(rho(i,j))>1e-6)
               if (0)//abs(rho(i,j))>1e-5)
               {
                  cout << "rho(" << i << "," << j << ") = " << rho(i,j)
                       << "   < " << a << " " << i << " | Vmon | "
                       <<  b << " " << j << " > = "  << Vmon(min(bra,ket),max(bra,ket))
//                       <<  b << " " << j << " > = "  << Vmon(bra,ket)
                       << "  Vab = " << Vab(a,b)/(ja+1) << endl;
               }
           }
         }
         Vab(a,b) /= (ja+1);
         Vab(b,a) = Vab(a,b);
      }
   }
   H = t + Vab;
/*   cout << "tab: " << endl;
   t.print();
   cout << "Vab: " << endl;
   Vab.print();
   cout << "Hab: " << endl;
   H.print();
*/
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

   ModelSpace * ms = OpIn->GetModelSpace();
   int norbits = ms->GetNumberOrbits();

   //Update the one-body part by multiplying by the matrix C(i,a) = <i|a>.

   OpHF->OneBody = C.t() * OpIn->OneBody * C;

   //Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz

   for (int J=0; J<JMAX; J++)
   {
      for (int parity = 0; parity <=1; parity++)
      {
         for (int Tz = -1; Tz<=1; Tz++)
         {
            TwoBodyChannel *chin = OpIn->GetTwoBodyChannel(J,parity,Tz);
            int npq = chin->Number_pq_configs;
            if (npq<1) continue;
            TwoBodyChannel *chhf = OpHF->GetTwoBodyChannel(J,parity,Tz);
            arma::mat D     = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
            arma::mat Dexch = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ba> = <ji|ab>
            for (int i=0; i<npq; i++) // loop over all possible original basis configurations <pq| in this J,p,Tz channel
            {   // i and j are the indices of the small matrix for this channel
               Ket * bra = chin->GetKet(i); // bra is in the original basis
               for (int j=0; j<npq; j++) // loop over all possible HF configurations |pq> in this J,p,Tz channel
               {
                  Ket * ket = chin->GetKet(j); // ket is in the HF basis
                  D(i,j) = C(bra->p,ket->p) * C(bra->q,ket->q) * 1./(1.0+bra->delta_pq());

                  Dexch(i,j) = C(bra->p,ket->q) * C(bra->q,ket->p) * ket->Phase(J) * 1./(1.0+bra->delta_pq());

               }
            }

           chhf->TBME  = D.t()     * chin->TBME * D;
           chhf->TBME += Dexch.t() * chin->TBME * D;
           chhf->TBME += D.t()     * chin->TBME * Dexch;
           chhf->TBME += Dexch.t() * chin->TBME * Dexch;
         }
      }
   }
}

