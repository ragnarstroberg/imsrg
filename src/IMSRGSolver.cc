
#include "IMSRGSolver.hh"
#include <boost/numeric/odeint.hpp>



// Constructor
IMSRGSolver::IMSRGSolver(const Operator &H_in)
   : ode_monitor(*this)
{
   method = "BCH";
   generator = "white";
   s = 0;
   ds = 0.1;
   ds_max = 0.5;
   smax  = 2.0;
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
   flowfile = "";
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
      ds = min(norm_domega / norm_eta / (norm_omega+1.0e-9), ds_max); 
      if (ds == ds_max) norm_domega /=2;
      if (s+ds > smax) ds = smax-s;
      s += ds;
      dOmega = Eta * ds; // Here's the Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
      Omega = dOmega.BCH_Product( Omega ); 

      // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
      H_s = H_0.BCH_Transform( Omega );
        
      UpdateEta();

      // Write details of the flow
      WriteFlowStatus(flowf);
      WriteFlowStatus(cout);

   }

   if (flowfile != "")
      flowf.close();

}




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



// Returns exp(Omega) OpIn exp(-Omega)
Operator IMSRGSolver::Transform(Operator& OpIn)
{
   return OpIn.BCH_Transform( Omega );
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
   else
   {
      cout << "Error. Unkown generator: " << generator << endl;
   }

}



// Epstein-Nesbet energy denominators for White-type generators
double IMSRGSolver::GetEpsteinNesbet1bDenominator(int i, int j)
{
   double denominator = H_s.OneBody(i,i) - H_s.OneBody(j,j) - H_s.GetTBMEmonopole(j,i,j,i);
   return denominator;
}



// This could likely be sped up by constructing and storing the monopole matrix
double IMSRGSolver::GetEpsteinNesbet2bDenominator(int ch, int ibra, int iket)
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int i = bra.p;
   int j = bra.q;
   int a = ket.p;
   int b = ket.q;
   double denominator = H_s.GetTBMEmonopole(i,j,i,j); // pp'pp'
   denominator       += H_s.GetTBMEmonopole(a,b,a,b); // hh'hh'
   denominator       -= H_s.GetTBMEmonopole(i,a,i,a); // phph
   denominator       -= H_s.GetTBMEmonopole(i,b,i,b); // ph'ph'
   denominator       -= H_s.GetTBMEmonopole(j,a,j,a); // p'hp'h
   denominator       -= H_s.GetTBMEmonopole(j,b,j,b); // p'h'p'h'

   denominator += H_s.OneBody(i,i)+ H_s.OneBody(j,j) - H_s.OneBody(a,a) - H_s.OneBody(b,b);

//   if (abs(denominator ) < 0.01)
//      cout << "2b denominator "  << ch << " " << ibra << "," << iket << " = " << denominator << endl;
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
      Orbit &oi = modelspace->GetOrbit(i);
      for (int &a : modelspace->holes)
      {
         Orbit &oa = modelspace->GetOrbit(a);
         double denominator = GetEpsteinNesbet1bDenominator(i,a);
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
      for (int& ibra : tbc.KetIndex_hh)
      {
         for (int& iket : tbc.KetIndex_pp)
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
   // One body piece -- make sure the valence one-body part is diagonal
   for ( int &i : modelspace->valence)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      for (int j=0; j<modelspace->GetNumberOrbits(); ++j)
      {
         if (i==j) continue;
         Orbit &oj = modelspace->GetOrbit(j);
         double denominator = GetEpsteinNesbet1bDenominator(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
   
   }
   // Two body piece -- eliminate ppvh and pqvv  ( vv'hh' was already accounted for with White )

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      // Decouple vv from qq and qv

      for (int& ibra : tbc.KetIndex_vv)
      {
         for (int& iket : tbc.KetIndex_particleq_particleq) 
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_holeq_holeq) 
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_particleq) 
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_holeq) 
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }


      }


      // Decouple hh states

      for (int& ibra : tbc.KetIndex_holeq_holeq)
      {
         for (int& iket : tbc.KetIndex_particleq_particleq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_particleq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

      // Decouple vh states

      for (int& ibra : tbc.KetIndex_v_holeq)
      {
         for (int& iket : tbc.KetIndex_particleq_particleq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_particleq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}




void IMSRGSolver::ConstructGenerator_Atan()
{
   // One body piece -- eliminate ph bits
   double maxnum = 0;
   double maxdenom = 0;
   int maxi = -1;
   int maxa = -1;
   for ( int &i : modelspace->particles)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      for (int &a : modelspace->holes)
      {
         Orbit &oa = modelspace->GetOrbit(a);
         double denominator = GetEpsteinNesbet1bDenominator(i,a);
         Eta.OneBody(i,a) = 0.5*atan(2*H_s.OneBody(i,a)/denominator);
         Eta.OneBody(a,i) = - Eta.OneBody(i,a);
         if ( abs(H_s.OneBody(i,a)) > abs(maxnum) ) 
         {
           maxnum = H_s.OneBody(i,a);
           maxdenom = denominator;
           maxi = i;
           maxa = a;
         }
      }
   }
//   cout << "Maximum one-body term: f(" << maxi << "," << maxa << ") = " << maxnum << " / " << maxdenom << " = " << 0.5*atan(2*maxnum/maxdenom) << endl;

   maxnum = 0;
   maxdenom = 0;
   maxi = -1;
   maxa = -1;
   int maxch = -1;
   // Two body piece -- eliminate pp'hh' bits
   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      for (int& ibra : tbc.KetIndex_pp)
      {
         for (int& iket : tbc.KetIndex_hh)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);

            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
            if ( abs(H_s.TwoBody[ch](ibra,iket)) > abs(maxnum) ) 
            {
              maxnum = H_s.TwoBody[ch](ibra,iket);
              maxdenom = denominator;
              maxi = ibra;
              maxa = iket;
              maxch = ch;
            }
         }
      }
    }
   TwoBodyChannel& tbcmax = modelspace->GetTwoBodyChannel(maxch);
   Ket & bra  = tbcmax.GetKet(maxi);
   Ket & ket  = tbcmax.GetKet(maxa);
//   cout << "Maximum two-body term: < " << bra.p << " " << bra.q << " | V | " << ket.p << " " << ket.q << " >"  << "(J=" << tbcmax.J << ") = "
//        << maxnum << " / " << maxdenom << " = " << 0.5*atan(2*maxnum/maxdenom) << endl;
}



void IMSRGSolver::ConstructGenerator_ShellModel_Atan()
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( int &i : modelspace->valence)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      for (int j=0; j<modelspace->GetNumberOrbits(); ++j)
      {
         if (i==j) continue;
         Orbit &oj = modelspace->GetOrbit(j);
         double denominator = GetEpsteinNesbet1bDenominator(i,j);
         Eta.OneBody(i,j) = 0.5*atan(2*H_s.OneBody(i,j)/denominator);
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
   
   }
   // Two body piece -- eliminate ppvh and pqvv  ( vv'hh' was already accounted for with White )

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      // Decouple vv from qq and qv

      for (int& ibra : tbc.KetIndex_vv)
      {
         for (int& iket : tbc.KetIndex_particleq_particleq) 
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_holeq_holeq) 
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_particleq) 
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_holeq) 
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = 0.5*atan(2*H_s.TwoBody[ch](ibra,iket) / denominator);
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }


      }


      // Decouple hh states

      for (int& ibra : tbc.KetIndex_holeq_holeq)
      {
         for (int& iket : tbc.KetIndex_particleq_particleq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_particleq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

      // Decouple vh states

      for (int& ibra : tbc.KetIndex_v_holeq)
      {
         for (int& iket : tbc.KetIndex_particleq_particleq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }

         for (int& iket : tbc.KetIndex_v_particleq)
         {
            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            Eta.TwoBody[ch](ibra,iket) = H_s.TwoBody[ch](ibra,iket) / denominator;
            Eta.TwoBody[ch](iket,ibra) = - Eta.TwoBody[ch](ibra,iket) ; // Eta needs to be antisymmetric
         }
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
   for( arma::mat& mtx : H_0.TwoBody )
   {
      N = mtx.n_cols;
      dim += N*(N+1)/2;
   }
   return dim;
}



void IMSRGSolver::WriteFlowStatus(ostream& f)
{
   if ( f.good() )
   {
      int fwidth = 16;
      int fprecision = 9;
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
