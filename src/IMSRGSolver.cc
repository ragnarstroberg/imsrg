
#include "IMSRGSolver.hh"



IMSRGSolver::IMSRGSolver(Operator H_in)
{
   method = "BCH";
   generator = "white";
   s = 0;
   ds = 0.1;
   smax  = 5.0;
   H_0 = H_in;
   H_s = H_in;
}



void IMSRGSolver::Solve()
{

   int imax = 50;
   for (int istep=0;i<imax;++istep)
   {
      UpdateEta();
      UpdateOmega();
      UpdateH();


   }

/*
  Do stuff...
*/

}



Operator IMSRGSolver::Transform(Operator& OpIn)
{
   return OpIn.BCH_Transform(Omega);
}


void IMSRGSolver::UpdateOmega()
{
   Omega = dOmega.BCH_Product(Omega);
}


void IMSRGSolver::UpdateH()
{
   H_s = H_s.BCH_Transform( dOmega );
}

void IMSRGSolver::UpdateEta()
{
   if (generator == "white")
   {
      H_diag = H_s;
      H_diag.ZeroBody = 0;
      for (int &a : modelspace->hole)
      {
         for (int &b : modelspace->valence)
         {
            H_diag(a,b) =0;
            H_diag(b,a) =0;
         }
      }

      for (int ch=0;ch<modelspace->GetNumberTowBodyChannels();++ch)
      {  // Note, should also decouple the v and q spaces
         H_diag.TwoBody[ch] = (P_hh*H_diag.TwoBody[ch] + P_pp*H_diag.TwoBody[ch]);
      }

      Eta = H_diag.Commutator(H_s);
   }

   dOmega = Eta * ds;
}






