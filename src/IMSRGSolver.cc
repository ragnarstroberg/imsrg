
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
   Eta = H_in;
   Eta.EraseZeroBody();
   Eta.EraseOneBody();
   Eta.EraseTwoBody();
   Eta.SetAntiHermitian();
   modelspace = H_0.GetModelSpace();
}



void IMSRGSolver::Solve()
{


   cout << "Norm of Hs = " << H_s.Norm() << endl;
   H_s.PrintTwoBody();

   UpdateEta();
   cout << "Norm of Eta = " << Eta.Norm() << endl;
   cout << "Eta ZeroBody:" << Eta.ZeroBody << endl;
   cout << "Eta OneBody: " << endl;
   Eta.OneBody.print();
   cout << "Eta TwoBody:" << endl;
   Eta.PrintTwoBody();

   cout << "=============== COMMUTATORS ==============" << endl;
   cout << "comm_110= " << Eta.comm110(H_s) << endl;
   cout << "comm_220= " << Eta.comm220(H_s) << endl;
   cout << "comm_111= " << endl;
   Eta.comm111(H_s).print();
//   cout << "comm_211= " << endl;
//   Eta.comm211(H_s).print();
//   cout << "comm_121= " << endl;
//   (Eta.comm121(H_s)).print();
   cout << "comm_121_both= " << endl;
   (Eta.comm121(H_s)).print();
   //(Eta.comm211(H_s)+ Eta.comm121(H_s)).print();
   cout << "norm = " << arma::norm(Eta.comm121(H_s),"frob" ) << endl;
   //cout << "norm = " << arma::norm(Eta.comm211(H_s)+Eta.comm121(H_s),"frob" ) << endl;
//   cout << "comm_221= " << endl;
//   Eta.comm221(H_s).print();
   cout << endl;
   arma::mat X2Y2_1 = Eta.comm221(H_s);
   cout << "comm_221 = " << endl;
   X2Y2_1.print();


  // Testing 2-body commutators
   Operator commOp = Eta;
   commOp.EraseOneBody();
   commOp.EraseTwoBody();
//   cout << " Eta2 H1: " << endl;
   //arma::mat X2Y1 = Eta.comm212(H_s,0);
//   Eta.comm212(H_s,commOp);
//   cout << "comm_212 ch0 = " << endl;
   //X2Y1.print();
//   commOp.TwoBody[0].print();
//   cout << endl;
//   cout << " H2 Eta1: " << endl;
   //arma::mat X1Y2 = Eta.comm122(H_s,0);
//   commOp.EraseTwoBody();

   Eta.comm122(H_s,commOp);
   cout << "comm_122 ch0 = " << endl;
   commOp.TwoBody[0].print();
   //X1Y2.print();
//   cout << "comm_212_both ch0 = " << endl;
//   (X2Y1 + X1Y2).print();

   //arma::mat X2Y2_2 = Eta.comm222_ph(H_s,0);
   commOp.EraseTwoBody();
   Eta.comm222_pp_hh(H_s,commOp);
   cout << "comm_222_pp_hh, ch=0 = " << endl;
   commOp.TwoBody[0].print();
   //X2Y2_2.print();
   //Eta.UpdateCrossCoupled();
   //cout << "Eta Cross-Coupled: " << endl;
   //Eta.TwoBody_CC[0].print();

  // Testing Cross-coupled nastiness...

   cout << endl << " Updating Hs Cross Coupled..." << endl;
   H_s.UpdateCrossCoupled();
   cout << endl << " Updating Eta Cross Coupled..." << endl;
   Eta.UpdateCrossCoupled();
   cout << "Hs Cross-Coupled left: " << endl;
   //H_s.TwoBody_CC[0].print();
      TwoBodyChannel_CC& tbc_print = modelspace->GetTwoBodyChannel_CC(120);
      cout << "J,p,t = " << tbc_print.J << " " << tbc_print.parity << " " << tbc_print.Tz << endl;
      for (int ik=0;ik<tbc_print.GetNumberKets();++ik)
      {
         Ket * ket = tbc_print.GetKet(ik);
         cout << "| " << ket->p << " " << ket->q << " >  ";
      }
      cout << endl;
   H_s.TwoBody_CC_left[151].print();
   cout << "Hs Cross-Coupled right: " << endl;
   H_s.TwoBody_CC_right[151].print();
   cout << " odd parity:" << endl;
   TwoBodyChannel_CC& tbc_print2 = modelspace->GetTwoBodyChannel_CC(150);
      cout << "J,p,t = " << tbc_print2.J << " " << tbc_print2.parity << " " << tbc_print2.Tz << endl;
      for (int ik=0;ik<tbc_print2.GetNumberKets();++ik)
      {
         Ket * ket = tbc_print2.GetKet(ik);
         cout << "| " << ket->p << " " << ket->q << " >  ";
      }
      cout << endl;
   cout << "Left" << endl;
   H_s.TwoBody_CC_left[151].print();
   cout << "Right" << endl;
   H_s.TwoBody_CC_right[151].print();
//   H_s.PrintTwoBody();
   cout << "Eta Cross-Coupled: " << endl;
   Eta.TwoBody_CC_left[151].print();
   cout << endl;
   Eta.TwoBody_CC_right[151].print();
//   Eta.PrintTwoBody();
   commOp.EraseTwoBody();
   Eta.comm222_ph(H_s,commOp);
   cout << "comm_222_ph, ch=0 = " << endl;
   commOp.TwoBody[0].print();
   commOp.PrintTwoBody();
   commOp.EraseTwoBody();
//   Eta.comm222_ph_slow(H_s,commOp);
//   cout << "THe slow way ... " << endl;
//   commOp.TwoBody[0].print();

//   TwoBodyChannel& tbc = Eta.GetModelSpace()->GetTwoBodyChannel(1);
//   cout << "J = " << tbc.J << "  Tz = " << tbc.Tz << "  P = " << tbc.parity << endl;
//   arma::mat X2Y2_2 = Eta.comm222_ph(H_s,1);
//   for (int ch=0;ch<Eta.GetModelSpace()->GetNumberTwoBodyChannels();++ch)
//   {
//    TwoBodyChannel& tbc = Eta.GetModelSpace()->GetTwoBodyChannel(ch);
//    cout <<"ch = " << ch << "  J = " << tbc.J << "  Tz = " << tbc.Tz << "  P = " << tbc.parity << endl;
//    arma::mat X2Y2_2 = Eta.comm222_ph(H_s,ch);
//    X2Y2_2.print();
//   }
//   cout << "comm_222_ph, ch=1 = " << endl;
//   X2Y2_2.print();


   return;
   int imax = 50;
   for (int istep=0;istep<imax;++istep)
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
   if (generator == "wegner")
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
   } // if wegner



   if (generator == "white")
   {
      cout << "Starting the one body part of ETA " << endl;
      // One body piece -- eliminate ph bits
      for ( int &i : modelspace->particles)
      {
         Orbit *oi = modelspace->GetOrbit(i);
         for (int &a : modelspace->holes)
         {
            Orbit *oa = modelspace->GetOrbit(a);
            if (oi->j2 != oa->j2 or oi->l != oa->l or oi->tz2 != oa->tz2) continue;
            double denominator = H_s.OneBody(i,i);
            denominator -= H_s.OneBody(a,a);
            denominator -= H_s.GetTBMEmonopole(a,i,a,i);
            Eta.OneBody(i,a) = H_s.OneBody(i,a)/denominator;
            Eta.OneBody(a,i) = - Eta.OneBody(i,a);
         }
      
      }
      cout << "Done with the one body part of ETA " << endl;
      // Two body piece -- eliminate pp'hh' bits
      // This could likely be sped up by constructing and storing the monopole matrix
      int nchan = modelspace->GetNumberTwoBodyChannels();
      for (int ch=0;ch<nchan;++ch)
      {
         TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
         if (tbc.GetNumberKets() < 1) continue;
         for (int& ibra : tbc.KetIndex_pp)
         {
            Ket * bra = tbc.GetKet(ibra);
            int i = bra->p;
            int j = bra->q;
            for (int& iket : tbc.KetIndex_hh)
            {
               Ket * ket = tbc.GetKet(iket);
               int a = ket->p;
               int b = ket->q;
               double denominator = H_s.GetTBMEmonopole(i,j,i,j); // pp'pp'
               denominator       += H_s.GetTBMEmonopole(a,b,a,b); // hh'hh'
               denominator       -= H_s.GetTBMEmonopole(i,a,i,a); // phph
               denominator       -= H_s.GetTBMEmonopole(i,b,i,b); // ph'ph'
               denominator       -= H_s.GetTBMEmonopole(j,a,j,a); // p'hp'h
               denominator       -= H_s.GetTBMEmonopole(j,b,j,b); // p'h'p'h'
               if (ch==0)
               {
                  cout << "ijab = " << i << "," << j << "," << a << "," << b << "  denom_A = " << denominator << endl;
                  cout << "TBMEmonopole(i,j,i,j) = " <<  H_s.GetTBMEmonopole(i,j,i,j) << endl;
                  cout << "TBMEmonopole(a,b,a,b) = " <<  H_s.GetTBMEmonopole(a,b,a,b) << endl;
                  cout << "TBMEmonopole(i,a,i,a) = " <<  H_s.GetTBMEmonopole(i,a,i,a) << endl;
                  cout << "TBMEmonopole(i,b,i,b) = " <<  H_s.GetTBMEmonopole(i,b,i,b) << endl;
                  cout << "TBMEmonopole(j,a,j,a) = " <<  H_s.GetTBMEmonopole(j,a,j,a) << endl;
                  cout << "TBMEmonopole(j,b,j,b) = " <<  H_s.GetTBMEmonopole(j,b,j,b) << endl;
                  cout << "      denom = " << denominator + H_s.OneBody(i,i)+ H_s.OneBody(j,j) - H_s.OneBody(a,a) - H_s.OneBody(b,b)
                       << "   numerator = " << H_s.TwoBody[ch](ibra,iket) << endl;
               }
               denominator += H_s.OneBody(i,i)+ H_s.OneBody(j,j) - H_s.OneBody(a,a) - H_s.OneBody(b,b);
               Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket)/denominator;
               Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ;
            }
         }
         
      }

   } // if white
//   ds = 1.0;
   dOmega = Eta * ds;
}






