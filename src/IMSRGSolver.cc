
#include "IMSRGSolver.hh"



// Constructor
IMSRGSolver::IMSRGSolver(Operator H_in)
{
   method = "BCH";
   generator = "white";
   s = 0;
   //ds = 10;
   ds = 0.1;
   ds_max = 1.0;
   smax  = 2.0;
   i_full_BCH = 5;
   norm_domega = 0.1;
   H_0 = H_in;
   H_s = H_in;
   Eta = H_in;
   Eta.EraseZeroBody();
   Eta.EraseOneBody();
   Eta.EraseTwoBody();
   Eta.SetAntiHermitian();
   modelspace = H_0.GetModelSpace();
   Omega = H_s;
   Omega.EraseZeroBody();
   Omega.EraseOneBody();
   Omega.EraseTwoBody();
   Omega.SetAntiHermitian();
   dOmega = Omega;
}



void IMSRGSolver::Solve()
{
   // If we have a flow output file, open it up and write to it here.
   ofstream flowf;
   if (flowfile != "")
      flowf.open(flowfile,ofstream::out);
   cout << " i     s       E0       ||H_1||      ||H_2||        ||Omega||     || Eta||    ||dOmega||     " << endl;

   UpdateEta();
   dOmega = Eta * ds; // Here's the Euler step.
//   double norm1 = dOmega.Norm();
   Operator H_last;
   Operator Omega_last;
   double norm_eta = Eta.Norm();
   double norm_eta_last = norm_eta;

   for (istep=0;s<smax;++istep)
   {
      // Write details of the flow
      WriteFlowStatus(flowf);
      WriteFlowStatus(cout);

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
      Omega = dOmega.BCH_Product( Omega ); 

      // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
      if (istep%i_full_BCH == i_full_BCH-1)
      {
         H_s = H_0.BCH_Transform( Omega );   
      }
      else
      {
         H_s = H_s.BCH_Transform( dOmega );  // less accurate, but converges with fewer commutators, since ||dOmega|| < ||Omega||
      }

      UpdateEta();   // generator for this step
      norm_eta = Eta.Norm();

      // choose ds such that ||dOmega|| = norm_domega (a fixed value)
      ds = min(norm_domega / norm_eta, ds_max); 
      if (ds == ds_max) norm_domega /=2;

      dOmega = Eta * ds; // Here's the Euler step.
      s += ds;
      norm_eta_last = norm_eta;
      H_last = H_s;
      Omega_last = Omega;

   }
   // if the last calculation of H_s was the quick way,
   // do it again the more accurate way.
   if (istep%i_full_BCH != i_full_BCH-1)
   {
      H_s = H_0.BCH_Transform( Omega ); 
   }

   if (flowfile != "")
      flowf.close();

}



Operator IMSRGSolver::Transform(Operator& OpIn)
{
   return OpIn.BCH_Transform(Omega);
}


void IMSRGSolver::UpdateOmega()
{
   Omega = dOmega.BCH_Product(Omega);
}



// Decide whether it's better to do the full
// transformation at each step, or just update
void IMSRGSolver::UpdateH()
{
   if (int(s/ds)%10 ==9)
      H_s = H_0.BCH_Transform( Omega );
   else
      H_s = H_s.BCH_Transform( dOmega );
}




void IMSRGSolver::UpdateEta()
{
   
   if (generator == "wegner")
   {
     ConstructGenerator_Wegner();
   } 
   else if (generator == "white")
   {
     ConstructGenerator_White();
   } 
   else if (generator == "shell-model")
   {
     ConstructGenerator_ShellModel();
   }

}



// Epstein-Nesbet energy denominators for White-type generators
double IMSRGSolver::GetEpsteinNesbet1bDenominator(int i, int j)
{
   double denominator = H_s.OneBody(i,i) - H_s.OneBody(j,j) - H_s.GetTBMEmonopole(j,i,j,i);
   if (abs(denominator ) < 0.01)
      cout << "1b denominator " << i << "," << j << " = " << denominator << endl;;
   return denominator;
}

double IMSRGSolver::GetEpsteinNesbet2bDenominator(int ch, int ibra, int iket)
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket * bra = tbc.GetKet(ibra);
   Ket * ket = tbc.GetKet(iket);
   int i = bra->p;
   int j = bra->q;
   int a = ket->p;
   int b = ket->q;
   double denominator = H_s.GetTBMEmonopole(i,j,i,j); // pp'pp'
   denominator       += H_s.GetTBMEmonopole(a,b,a,b); // hh'hh'
   denominator       -= H_s.GetTBMEmonopole(i,a,i,a); // phph
   denominator       -= H_s.GetTBMEmonopole(i,b,i,b); // ph'ph'
   denominator       -= H_s.GetTBMEmonopole(j,a,j,a); // p'hp'h
   denominator       -= H_s.GetTBMEmonopole(j,b,j,b); // p'h'p'h'

   denominator += H_s.OneBody(i,i)+ H_s.OneBody(j,j) - H_s.OneBody(a,a) - H_s.OneBody(b,b);

   if (abs(denominator ) < 0.01)
      cout << "2b denominator "  << ch << " " << ibra << "," << iket << " = " << denominator << endl;
   return denominator;
}



void IMSRGSolver::ConstructGenerator_Wegner()
{
   H_diag = H_s;
   H_diag.ZeroBody = 0;
   for (int &a : modelspace->holes)
   {
      for (int &b : modelspace->valence)
      {
         H_diag.OneBody(a,b) =0;
         H_diag.OneBody(b,a) =0;
      }
   }

   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {  // Note, should also decouple the v and q spaces
      // This is wrong. The projection operator should be different.
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      H_diag.TwoBody[ch] = (tbc.Proj_hh*H_diag.TwoBody[ch] + tbc.Proj_pp*H_diag.TwoBody[ch]);
   }

   Eta = H_diag.Commutator(H_s);
}



void IMSRGSolver::ConstructGenerator_White()
{
   // One body piece -- eliminate ph bits
   for ( int &i : modelspace->particles)
   {
      Orbit *oi = modelspace->GetOrbit(i);
      for (int &a : modelspace->holes)
      {
         Orbit *oa = modelspace->GetOrbit(a);
         double denominator = GetEpsteinNesbet1bDenominator(i,a);
         Eta.OneBody(i,a) = H_s.OneBody(i,a)/denominator;
         Eta.OneBody(a,i) = - Eta.OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   // This could likely be sped up by constructing and storing the monopole matrix
   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      for (int& ibra : tbc.KetIndex_pp)
      {
         for (int& iket : tbc.KetIndex_hh)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);

            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}







void IMSRGSolver::ConstructGenerator_ShellModel()
{
   ConstructGenerator_White(); // Start with the White generator

   // One body piece -- make sure the valence one-body part is diagonal
   for ( int &i : modelspace->valence)
   {
      Orbit *oi = modelspace->GetOrbit(i);
      for (int &j : modelspace->particles)
      {
         if (i==j) continue;
         Orbit *oj = modelspace->GetOrbit(j);
         double denominator = GetEpsteinNesbet1bDenominator(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
         if (abs(Eta.OneBody(i,j)) > 1e5)
         {
            cout << "Eta OneBody(" << i << "," << j << ") =  "<< Eta.OneBody(i,j) << endl;
         }
      }
   
   }
   // Two body piece -- eliminate ppvh and pqvv  ( vv'hh' was already accounted for with White )

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);


      for (int& ibra : tbc.KetIndex_pp)
      {
         for (int& iket : tbc.KetIndex_vh)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
            if (abs(Eta.TwoBody[ch](ibra,iket)) > 1e5)
               cout << "Eta TwoBody_ppvh(" << ibra << "," << iket << ") =  "<< Eta.TwoBody[ch](ibra,iket) << endl;
         }
      }

      for (int& ibra : tbc.KetIndex_vv)
      {
         for (int& iket : tbc.KetIndex_qq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
            if (abs(Eta.TwoBody[ch](ibra,iket)) > 1e5)
               cout << "Eta TwoBody_vvqq(" << ibra << "," << iket << ") =  "<< Eta.TwoBody[ch](ibra,iket) << endl;
         }
      }

      for (int& ibra : tbc.KetIndex_vv)
      {
         for (int& iket : tbc.KetIndex_vq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
            if (abs(Eta.TwoBody[ch](ibra,iket)) > 1e5)
            {
               Ket * bra = tbc.GetKet(ibra);
               Ket * ket = tbc.GetKet(iket);
               cout << "Eta TwoBody_vvvq(" << bra->p << "," << bra->q << "," << ket->p << "," << ket->q << ") (J="<<tbc.J << ") =  "<< Eta.TwoBody[ch](ibra,iket)
                    << ",  denominator = " << denominator << endl;
               
            }
         }
      }


    }
}




void IMSRGSolver::WriteFlowStatus(ostream& f)
{
//      cout << istep << "      " << s << "      " << H_s.ZeroBody << "     " << H_s.OneBodyNorm() << "    " << H_s.TwoBodyNorm() << "     " << Omega.Norm() << "     "  << Eta.Norm() << "   "  << dOmega.Norm() << endl;
      if ( f.good() )
         f << istep << "      " << s << "      " << H_s.ZeroBody << "     " << H_s.OneBodyNorm() << "    " << H_s.TwoBodyNorm() << "     " << Omega.Norm() << "     " << Eta.Norm() << "   " << dOmega.Norm() << endl;

}

