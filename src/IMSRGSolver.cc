
#include "IMSRGSolver.hh"
#include <iomanip>

#ifndef NO_ODE
#include <boost/numeric/odeint.hpp>
#endif


IMSRGSolver::IMSRGSolver()
#ifndef NO_ODE
    :ode_monitor(*this)
#endif
{
   method = "BCH";
   generator = "white";
   s = 0;
   ds = 0.1;
   ds_max = 0.5;
   smax  = 2.0;
   norm_domega = 0.1;
   flowfile = "";
}

// Constructor
IMSRGSolver::IMSRGSolver(const Operator &H_in)
   : H_0(H_in), H_s(H_in), Eta(H_in), Omega(H_in) ,dOmega(H_in)
#ifndef NO_ODE
    ,ode_monitor(*this)
#endif
{
   method = "BCH";
   generator = "white";
   s = 0;
   ds = 0.1;
   ds_max = 0.5;
   smax  = 2.0;
   norm_domega = 0.1;
   flowfile = "";
   modelspace = H_in.GetModelSpace();
   Eta.Erase();
   Eta.SetAntiHermitian();
   Omega.Erase();
   Omega.SetAntiHermitian();
   dOmega.Erase();
   dOmega.SetAntiHermitian();
}

void IMSRGSolver::SetHin(const Operator & H_in)
{
   modelspace = H_in.GetModelSpace();
   H_0 = Operator(H_in);
   H_s = H_0;
   Eta = H_0;
   Eta.Erase();
   Eta.SetAntiHermitian();
   Omega = Eta;
   dOmega = Eta;  
}

void IMSRGSolver::Reset()
{
   s=0;
   Eta.Erase();
   Omega.Erase();
   dOmega.Erase();
}

void IMSRGSolver::Solve()
{
   // If we have a flow output file, open it up and write to it here.
   ofstream flowf;
   if (flowfile != "")
      flowf.open(flowfile,ofstream::out);
   WriteFlowStatusHeader(cout);

   istep = 0;
   UpdateEta();

    // Write details of the flow
   WriteFlowStatus(flowf);
   WriteFlowStatus(cout);

   for (istep=1;s<smax;++istep)
   {

      double norm_eta = Eta.Norm();
      double norm_omega = Omega.Norm();
      // ds should never be more than 1, as this is over-rotating
      //ds = min(norm_domega / norm_eta / (norm_omega+1.0e-9), ds_max); 
      ds = min(min(norm_domega/norm_eta, norm_domega / norm_eta / (norm_omega+1.0e-9)), ds_max); 
      if (ds == ds_max) norm_domega /=2;
      if (s+ds > smax) ds = smax-s;
      s += ds;
      dOmega = Eta * ds; // Here's the Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
//      cout << "Updating Omega" << endl;
      Omega = dOmega.BCH_Product( Omega ); 

      // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
//      cout << "Udating H" << endl;
      H_s = H_0.BCH_Transform( Omega );
        
      UpdateEta();

      // Write details of the flow
      WriteFlowStatus(flowf);
      WriteFlowStatus(cout);

   }

   if (flowfile != "")
      flowf.close();

}


#ifndef NO_ODE

void IMSRGSolver::Solve_ode()
{
   if (flowfile != "")
   {
     ofstream flowf;
     flowf.open(flowfile,ofstream::out);
     WriteFlowStatus(flowf);
     flowf.close();
   }
   WriteFlowStatusHeader(cout);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
   runge_kutta4<Operator, double, Operator, double, vector_space_algebra> stepper;
   auto system = std::bind( &IMSRGSolver::ODE_systemH, *this, pl::_1, pl::_2, pl::_3);
   auto monitor = ode_monitor;
   size_t steps = integrate_const(stepper, system, H_s, s, smax, ds, monitor);
   monitor.report();
}


void IMSRGSolver::ODE_systemH(const Operator& x, Operator& dxdt, const double t)
{
   H_s = x;
   s = t;
   UpdateEta();
   dxdt = Eta.Commutator(x);
   WriteFlowStatus(cout);
   if (flowfile != "")
   {
     ofstream flowf;
     flowf.open(flowfile,ofstream::app);
     WriteFlowStatus(flowf);
     flowf.close();
   }
}


void IMSRGSolver::Solve_ode_magnus()
{
   if (flowfile != "")
   {
     ofstream flowf;
     flowf.open(flowfile,ofstream::out);
     WriteFlowStatus(flowf);
     flowf.close();
   }
   WriteFlowStatusHeader(cout);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
   runge_kutta4<Operator, double, Operator, double, vector_space_algebra> stepper;
   auto system = std::bind( &IMSRGSolver::ODE_systemOmega, *this, pl::_1, pl::_2, pl::_3);
   auto monitor = ode_monitor;
   size_t steps = integrate_const(stepper, system, Omega, s, smax, ds, monitor);
   monitor.report();
}

void IMSRGSolver::ODE_systemOmega(const Operator& x, Operator& dxdt, const double t)
{
   s = t;
   Omega = x;
   H_s = H_0.BCH_Transform(Omega);
   UpdateEta();
   dxdt = Eta - 0.5*Omega.Commutator(Eta);
   WriteFlowStatus(cout);
   if (flowfile != "")
   {
     ofstream flowf;
     flowf.open(flowfile,ofstream::app);
     WriteFlowStatus(flowf);
     flowf.close();
   }
}

#endif


// Returns exp(Omega) OpIn exp(-Omega)
Operator IMSRGSolver::Transform(Operator& OpIn)
{
   return OpIn.BCH_Transform( Omega );
}

// Returns exp(-Omega) OpIn exp(Omega)
Operator IMSRGSolver::InverseTransform(Operator& OpIn)
{
   return OpIn.BCH_Transform( -Omega );
}

void IMSRGSolver::UpdateEta()
{
   Eta.EraseOneBody();
   Eta.EraseTwoBody();
   if (generator == "wegner")
   {
     ConstructGenerator_Wegner();
   } 
   else if (generator == "white")
   {
     ConstructGenerator_White();
   } 
   else if (generator == "atan")
   {
     ConstructGenerator_Atan();
   } 
   else if (generator == "shell-model")
   {
     ConstructGenerator_ShellModel();
   }
   else if (generator == "shell-model-atan")
   {
     ConstructGenerator_ShellModel_Atan();
   }
   else if (generator == "hartree-fock")
   {
     ConstructGenerator_HartreeFock();
   }
   else
   {
      cout << "Error. Unkown generator: " << generator << endl;
   }

}



// Epstein-Nesbet energy denominators for White-type generators
// i=particle, j=hole
//double IMSRGSolver::GetEpsteinNesbet1bDenominator(int i, int j) 
double IMSRGSolver::Get1bDenominator_ph(int i, int j) 
{
//   double denominator = H_s.OneBody(i,i) - H_s.OneBody(j,j) - H_s.GetTBMEmonopole(j,i,j,i);
   double denominator = H_s.OneBody(i,i) - H_s.OneBody(j,j) - H_s.GetTBMEmonopole(i,j,i,j);
   return denominator;
}



// This could likely be sped up by constructing and storing the monopole matrix
// bra=pp'  ket = hh'
//double IMSRGSolver::GetEpsteinNesbet2bDenominator(int ch, int ibra, int iket) 
double IMSRGSolver::Get2bDenominator_pphh(int ch, int ibra, int iket) 
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int i = bra.p;
   int j = bra.q;
   int a = ket.p;
   int b = ket.q;
   double denominator = H_s.OneBody(i,i)+ H_s.OneBody(j,j) - H_s.OneBody(a,a) - H_s.OneBody(b,b);

   denominator       += H_s.GetTBMEmonopole(i,j,i,j); // pp'pp'
   denominator       += H_s.GetTBMEmonopole(a,b,a,b); // hh'hh'
   denominator       -= H_s.GetTBMEmonopole(i,a,i,a); // phph
   denominator       -= H_s.GetTBMEmonopole(i,b,i,b); // ph'ph'
   denominator       -= H_s.GetTBMEmonopole(j,a,j,a); // p'hp'h
   denominator       -= H_s.GetTBMEmonopole(j,b,j,b); // p'h'p'h'

   return denominator;
}
// i=p j=p'
double IMSRGSolver::Get1bDenominator_pp(int i, int j)
{
   double denominator = H_s.OneBody(i,i) - H_s.OneBody(j,j);
   if (abs(denominator) < 1e-6) return 1e0;
   return denominator;
}

// bra = pp'  ket=hp
double IMSRGSolver::Get2bDenominator_pphp(int ch, int ibra, int iket)
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int p1 = bra.p;
   int p2 = bra.q;
   int h = ket.p;
   int p3 = ket.q;
   double denominator = H_s.OneBody(p1,p1)+ H_s.OneBody(p2,p2) - H_s.OneBody(p3,p3) - H_s.OneBody(h,h);

   denominator       += H_s.GetTBMEmonopole(p1,p2,p1,p2); // phph
   denominator       -= H_s.GetTBMEmonopole(p1,h,p1,h); // pp'pp'
   denominator       -= H_s.GetTBMEmonopole(p2,h,p2,h); // hh'hh'

   return denominator;
}
// bra = pq  ket=vv'
double IMSRGSolver::Get2bDenominator_pppp(int ch, int ibra, int iket)
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int p1 = bra.p;
   int p2 = bra.q;
   int p3 = ket.p;
   int p4 = ket.q;
   double denominator = H_s.OneBody(p1,p1)+ H_s.OneBody(p2,p2) - H_s.OneBody(p3,p3) - H_s.OneBody(p4,p4);

   denominator       += H_s.GetTBMEmonopole(p1,p2,p1,p2); // pp'pp'
   denominator       -= H_s.GetTBMEmonopole(p3,p4,p3,p4); // hh'hh'

   return denominator;
}



// I haven't used this, so I don't know if it's right.
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
   //   H_diag.TwoBody[ch] = (tbc.Proj_hh*H_diag.TwoBody[ch] + tbc.Proj_pp*H_diag.TwoBody[ch]);
      H_diag.TwoBody[ch].at(ch) = (tbc.Proj_hh + tbc.Proj_pp) * H_diag.TwoBody[ch].at(ch);
   }

   Eta = H_diag.Commutator(H_s);
}



void IMSRGSolver::ConstructGenerator_White()
{
   // One body piece -- eliminate ph bits
   for ( int &i : modelspace->particles)
   {
      for (int &a : modelspace->holes)
      {
         //double denominator = GetEpsteinNesbet1bDenominator(i,a);
         double denominator = Get1bDenominator_ph(i,a);
         Eta.OneBody(i,a) = H_s.OneBody(i,a)/denominator;
         Eta.OneBody(a,i) = - Eta.OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits. Note that the hh'hp pieces are accounted
   // for in the normal ordered one-body part.
   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto& ETA2 = Eta.TwoBody[ch].at(ch);
      auto& H2 = H_s.TwoBody[ch].at(ch);
      for (unsigned int& ibra : tbc.KetIndex_pp)
      {
         for (unsigned int& iket : tbc.KetIndex_hh)
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) =  H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}


void IMSRGSolver::ConstructGenerator_Atan()
{
   // One body piece -- eliminate ph bits
   for ( int &i : modelspace->particles)
   {
      for (int &a : modelspace->holes)
      {
         //double denominator = GetEpsteinNesbet1bDenominator(i,a);
         double denominator = Get1bDenominator_ph(i,a);
         Eta.OneBody(i,a) = 0.5*atan(2*H_s.OneBody(i,a)/denominator);
         Eta.OneBody(a,i) = - Eta.OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta.nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta.TwoBody[ch].at(ch);
      arma::mat& H2 = H_s.TwoBody[ch].at(ch);
      for (unsigned int& ibra : tbc.KetIndex_pp)
      {
         for (unsigned int& iket : tbc.KetIndex_hh)
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);

            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}






void IMSRGSolver::ConstructGenerator_ShellModel()
{
   // One body piece -- make sure the valence one-body part is diagonal
   // no excitations out of the core

   Eta.Erase();

   for (int &i : modelspace->particles )
//   for (int &i : modelspace->valence )
   {
      for (int &j : modelspace->holes )
      {
         double denominator = Get1bDenominator_ph(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
//         if (abs(Eta.OneBody(i,j)) > 10)
//         {
//           cout <<"Eta(" << i << "," << j << ") = " << Eta.OneBody(i,j) << "  -- " << H_s.OneBody(i,j) << " / " << denominator << endl;
//         }
      }
      for (int &j : modelspace->particles )
      {
         if (i==j) continue;
         double denominator = Get1bDenominator_pp(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
         if (abs(Eta.OneBody(i,j)) > 10)
         {
           cout << "Eta(" << i << "," << j << ") = " << Eta.OneBody(i,j) << "  -- " << H_s.OneBody(i,j) << " / " << denominator
                << "    ii " << H_s.OneBody(i,i) << "  jj " << H_s.OneBody(j,j) << endl;
         }
      }
   }


   // Two body piece 
   // we want to drive to zero any terms that take 
   //   |hh> to |pp>
   //   |vh> to |pp>
   //   |vv> to |qq> or |vq>

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& ETA2 = Eta.TwoBody[ch].at(ch);
      auto& H2 = H_s.TwoBody[ch].at(ch);

      // Decouple vv states from pq states
      for (int& iket : tbc.KetIndex_vv)
      {
         // < qq' | vv' >
         for (int& ibra : tbc.KetIndex_particleq_particleq) 
         {
            double denominator = Get2bDenominator_pppp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vq | vv' > 
         for (int& ibra : tbc.KetIndex_v_particleq) 
         {
            double denominator = Get2bDenominator_pppp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }


      // Decouple hh states

      for (int& iket : tbc.KetIndex_holeq_holeq)
      {
         // < qq' | hh' >
         for (int& ibra : tbc.KetIndex_particleq_particleq)
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         // < vq | hh' >
         for (int& ibra : tbc.KetIndex_v_particleq)
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         // < vv | hh' >
         for (int& ibra : tbc.KetIndex_vv)
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple vh states

      for (int& iket : tbc.KetIndex_v_holeq)
      {
         // < qq | vh >
         for (int& ibra : tbc.KetIndex_particleq_particleq)
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         // < vq | vh >
         for (int& ibra : tbc.KetIndex_v_particleq)
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         // < vv | vh >
         for (int& ibra : tbc.KetIndex_vv)
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

    }
}



void IMSRGSolver::ConstructGenerator_ShellModel_Atan()
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( int &i : modelspace->valence)
   {
//      for (int j=0; j<modelspace->GetNumberOrbits(); ++j)
      for (int &j : modelspace->holes)
      {
         double denominator = Get1bDenominator_ph(i,j);
         Eta.OneBody(i,j) = 0.5*atan(2*H_s.OneBody(i,j)/denominator);
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
      for (int &j : modelspace->particles)
      {
         if (i==j) continue;
//         double denominator = GetEpsteinNesbet1bDenominator(i,j);
         double denominator = Get1bDenominator_pp(i,j);
         Eta.OneBody(i,j) = 0.5*atan(2*H_s.OneBody(i,j)/denominator);
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
   
   }
   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta.TwoBody[ch].at(ch);
      arma::mat& H2 =  H_s.TwoBody[ch].at(ch);
      // Decouple vv from qq and qv

      for (int& ibra : tbc.KetIndex_vv)
      {
         for (int& iket : tbc.KetIndex_particleq_particleq) 
         {
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator_pppp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
//            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_holeq_holeq) 
         {
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
//            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_particleq) 
         {
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator_pppp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
//            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_holeq) 
         {
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
//            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }


      }


      // Decouple hh states

      for (int& iket : tbc.KetIndex_holeq_holeq)
      {
         for (int& ibra : tbc.KetIndex_particleq_particleq)
         {
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
//            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& ibra : tbc.KetIndex_v_particleq)
         {
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
//            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

      // Decouple vh states

      for (int& iket : tbc.KetIndex_v_holeq)
      {
         for (int& ibra : tbc.KetIndex_particleq_particleq)
         {
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
//            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& ibra : tbc.KetIndex_v_particleq)
         {
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
//            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}

void IMSRGSolver::ConstructGenerator_HartreeFock()
{
   Eta.SetParticleRank(1);
   // One body piece -- eliminate ph bits
   for ( int &a : modelspace->holes)
   {
      for (int &b : modelspace->holes)
      {
         // Note that for the hole-hole case, the denominator
         // has a minus sign relative to the pp case.
         if (a==b) continue;
         double denominator = -Get1bDenominator_pp(a,b);
         Eta.OneBody(a,b) = H_s.OneBody(a,b)/denominator;
         Eta.OneBody(b,a) = - Eta.OneBody(a,b);
      }
   }
   for ( int &i : modelspace->particles)
   {
      for (int &a : modelspace->holes)
      {
         double denominator = Get1bDenominator_ph(i,a);
         Eta.OneBody(i,a) = H_s.OneBody(i,a)/denominator;
         Eta.OneBody(a,i) = - Eta.OneBody(i,a);
      }
   }  
   for ( int &i : modelspace->particles)
   {
      for (int &j : modelspace->particles)
      {
         if (i==j) continue;
         double denominator = Get1bDenominator_pp(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
   }  
}



// count number of equations to be solved
int IMSRGSolver::GetSystemDimension()
{
   int dim = 1; // zero-body part

   // one-body part
   int N = H_0.OneBody.n_cols;
   dim += N*(N+1)/2;
   for( auto channels : H_0.TwoBody )
   {
     for( auto mtx : channels )
     {
        N = mtx.second.n_cols;
        dim += N*(N+1)/2;
     }
   }
   return dim;
}



void IMSRGSolver::WriteFlowStatus(ostream& f)
{
   if ( f.good() )
   {
      int fwidth = 16;
      int fprecision = 9;
      f.setf(ios::fixed);
      //f.width(20);
//      f.precision(10);
      f << setw(5) << istep
        << setw(fwidth) << setprecision(3) << s
        << setw(fwidth) << setprecision(fprecision) << H_s.ZeroBody 
        << setw(fwidth) << setprecision(fprecision) << H_s.OneBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << H_s.TwoBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << Omega.Norm()
        << setw(fwidth) << setprecision(fprecision) << Eta.OneBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << Eta.TwoBodyNorm()
//        << setw(fwidth) << setprecision(fprecision) << dOmega.Norm()
        << endl;
//      f << istep << "      " << s << "      " << H_s.ZeroBody << "     " << H_s.OneBodyNorm() << "    " << H_s.TwoBodyNorm() << "     " << Omega.Norm() << "     " << Eta.OneBodyNorm() << "    " << Eta.TwoBodyNorm() << "   " << dOmega.Norm() << endl;
   }

}

void IMSRGSolver::WriteFlowStatusHeader(ostream& f)
{
   if ( f.good() )
   {
      int fwidth = 16;
      int fprecision = 9;
      f.setf(ios::fixed);
      f << setw(5) << "i"
        << setw(fwidth) << setprecision(3) << "s"
        << setw(fwidth) << setprecision(fprecision) << "E0"
        << setw(fwidth) << setprecision(fprecision) << "||H_1||" 
        << setw(fwidth) << setprecision(fprecision) << "||H_2||" 
        << setw(fwidth) << setprecision(fprecision) << "||Omega||" 
        << setw(fwidth) << setprecision(fprecision) << "||Eta_1||" 
        << setw(fwidth) << setprecision(fprecision) << "||Eta_2||" 
//        << setw(fwidth) << setprecision(fprecision) << "||dOmega||" 
        << endl;
      f << "-----------------------------------------------------------------------------------------------------------------------" << endl;
   }

}
