
#include "IMSRGSolver.hh"
#include <iomanip>

#ifndef NO_ODE
#include <boost/numeric/odeint.hpp>
#endif

IMSRGSolver::~IMSRGSolver()
{
//   cout << "In IMSRGSolver destructor." << endl;
}

IMSRGSolver::IMSRGSolver()
#ifndef NO_ODE
    :ode_monitor(*this)
#endif
{
   method = "BCH";
   s = 0;
   ds = 0.1;
   ds_max = 0.5;
   smax  = 2.0;
   norm_domega = 0.1;
   omega_norm_max = 2.0;
   flowfile = "";
}

// Constructor
IMSRGSolver::IMSRGSolver( Operator &H_in)
//   : H_0(&H_in), H_s(H_in), Eta(H_in), Omega(H_in)// ,dOmega(H_in)
   : H_0(&H_in), H_s(H_in), Eta(H_in) // ,dOmega(H_in)
#ifndef NO_ODE
    ,ode_monitor(*this)
#endif
{
   method = "BCH";
   istep = 0;
   s = 0;
   ds = 0.1;
   ds_max = 0.5;
   smax  = 2.0;
   norm_domega = 0.1;
   omega_norm_max = 2.0;
   flowfile = "";
   modelspace = H_in.GetModelSpace();
   Eta.Erase();
   Eta.SetAntiHermitian();
   Omega.push_back( Eta);
//   Omega.Erase();
//   Omega.SetAntiHermitian();
}

void IMSRGSolver::SetHin( Operator & H_in)
{
   modelspace = H_in.GetModelSpace();
   H_0 = &H_in;
   H_s = H_in;
   Eta = H_in;
   Eta.Erase();
   Eta.SetAntiHermitian();
//   Omega = Eta;
   if (Omega.back().Norm() > 1e-6)
     Omega.push_back(Eta);
}

void IMSRGSolver::Reset()
{
   s=0;
   Eta.Erase();
//   Omega.Erase();
   Omega.resize(0);
   Omega.push_back(Eta);
}


void IMSRGSolver::SetFlowFile(string s)
{
   flowfile = s;
   ofstream flowf;
   if (flowfile != "")
   {
      flowf.open(flowfile,ofstream::out);
      flowf.close();
   }
   WriteFlowStatusHeader(cout);
}

void IMSRGSolver::Solve()
{
   istep = 0;
   generator.Update(&H_s,&Eta);

    // Write details of the flow
   WriteFlowStatus(flowfile);
   WriteFlowStatus(cout);

   for (istep=1;s<smax;++istep)
   {

      double norm_eta = Eta.Norm();
//      double norm_omega = Omega.Norm();
      double norm_omega = Omega.back().Norm();
      if (norm_omega > omega_norm_max)
      {
        H_saved = H_s;
        Omega.push_back(Eta);
        Omega.back().Erase();
        norm_omega = 0;
        cout << "pushing back another Omega." << endl;
      }
      // ds should never be more than 1, as this is over-rotating
      //ds = min(norm_domega / norm_eta / (norm_omega+1.0e-9), ds_max); 
      ds = min(min(norm_domega/norm_eta, norm_domega / norm_eta / (norm_omega+1.0e-9)), ds_max); 
//      if (ds == ds_max) norm_domega /=2;
      if (s+ds > smax) ds = smax-s;
      s += ds;
      Eta *= ds; // Here's the Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
//      Omega = Eta.BCH_Product( Omega ); 
      Omega.back() = Eta.BCH_Product( Omega.back() ); 

      // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
//      H_s = H_0->BCH_Transform( Omega );
      if (Omega.size()<2)
      {
        H_s = H_0->BCH_Transform( Omega.back() );
      }
      else
      {
        H_s = H_saved.BCH_Transform( Omega.back() );
      }
//      for (auto omega : Omega)
//      {
//        H_s = H_s.BCH_Transform( omega );
//      }
        
      generator.Update(&H_s,&Eta);

      // Write details of the flow
      WriteFlowStatus(flowfile);
      WriteFlowStatus(cout);

   }

}


#ifndef NO_ODE
// Implement element-wise division and abs and reduce for Operators.
// This is required for adaptive steppers


/*
Operator operator/ (const Operator& num, const Operator& denom)
{
   Operator quotient = num;
   quotient.ZeroBody /= denom.ZeroBody;
   quotient.OneBody /= denom.OneBody;
   for (int ch=0;ch<quotient.nChannels;++ch)
   {
      for (auto &twobody : quotient.TwoBody.MatEl[ch] )
      {
        int chket = twobody.first;
        twobody.second /= denom.TwoBody.GetMatrix(ch,chket);
      }
   }

}

Operator abs(const Operator& opin)
{
   Operator opout = opin;
   opout.ZeroBody = abs(opout.ZeroBody);
   opout.OneBody = arma::abs(opout.OneBody);
   for (int ch=0;ch<opout.nChannels;++ch)
   {
      for (auto &twobody : opout.TwoBody.MatEl[ch] )
      {
        twobody.second = arma::abs(twobody.second);
      }
   }

}
// Black magic...
namespace boost {namespace numeric {namespace odeint{
template<>
struct vector_space_reduce<Operator&>
{
   template<class Value, class Op>
   Value operator()(const Operator& X, Op op, Value init)
   {
      init = op(init, X.ZeroBody());
      for ( double& v : X.OneBody )
      {
         init = op(init, v) ;
      }
      for (int ch=0;ch<X.nChannels;++ch)
      {
        for (auto &twobody : X.TwoBody[ch] )
        {
          for (double& v : twobody )
          {
             init = op(init, v) ;
          }
        }
      }
      return init;
   }
};
}}}
*/

void IMSRGSolver::Solve_ode()
{

   WriteFlowStatusHeader(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
   runge_kutta4<Operator, double, Operator, double, vector_space_algebra> stepper;
   auto system = std::bind( &IMSRGSolver::ODE_systemH, *this, pl::_1, pl::_2, pl::_3);
   auto monitor = ode_monitor;
   size_t steps = integrate_const(stepper, system, H_s, s, smax, ds, monitor);
   monitor.report();
}

void IMSRGSolver::Solve_ode_adaptive()
{
/*
   WriteFlowStatusHeader(cout);
     WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
//   auto stepper = make_controlled<runge_kutta_dopri5<Operator, double, Operator, double, vector_space_algebra> >(1e-8,1e-8);
   auto system = std::bind( &IMSRGSolver::ODE_systemH, *this, pl::_1, pl::_2, pl::_3);
   auto monitor = ode_monitor;
   auto stepper = make_controlled(1e-8,1e-8, runge_kutta_dopri5<Operator, double, Operator, double, vector_space_algebra>() );
   size_t steps = integrate_adaptive(stepper, system, H_s, s, smax, ds, monitor);
   monitor.report();
*/
}


void IMSRGSolver::ODE_systemH( Operator& x, Operator& dxdt, const double t)
{
   H_s = x;
   s = t;
//   UpdateEta();
   generator.Update(&H_s,&Eta);
//   dxdt = Eta.Commutator(x);
   dxdt = Commutator(Eta,x);
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);
}


void IMSRGSolver::Solve_ode_magnus()
{
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
   runge_kutta4<Operator, double, Operator, double, vector_space_algebra> stepper;
   auto system = std::bind( &IMSRGSolver::ODE_systemOmega, *this, pl::_1, pl::_2, pl::_3);
   auto monitor = ode_monitor;
//   size_t steps = integrate_const(stepper, system, Omega, s, smax, ds, monitor);
   size_t steps = integrate_const(stepper, system, Omega.front(), s, smax, ds, monitor);
   monitor.report();
}

void IMSRGSolver::ODE_systemOmega( Operator& x, Operator& dxdt, const double t)
{
   s = t;
//   Omega = x;
   Omega.front() = x;
//   H_s = H_0->BCH_Transform(Omega);
   H_s = H_0->BCH_Transform(Omega.front());
   generator.Update(&H_s,&Eta);
//   dxdt = Eta - 0.5*Omega.Commutator(Eta);
//   dxdt = Eta - 0.5*Commutator(Omega,Eta);
   dxdt = Eta - 0.5*Commutator(Omega.front(),Eta);
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);

}

#endif


/// Returns \f$ e^{Omega} \mathcal{O} e^{-Omega} \f$
Operator IMSRGSolver::Transform(Operator& OpIn)
{
  Operator OpOut = OpIn;
  for (auto omega : Omega )
  {
    OpOut = OpOut.BCH_Transform( omega );
  }
  return OpOut;
//   return OpIn.BCH_Transform( Omega );
}

Operator IMSRGSolver::Transform(Operator&& OpIn)
{
  Operator OpOut = OpIn;
  for (auto omega : Omega )
  {
    OpOut = OpOut.BCH_Transform( omega );
  }
  return OpOut;
//   return OpIn.BCH_Transform( Omega );
}


/// Returns \f$ e^{-Omega} \mathcal{O} e^{Omega} \f$
Operator IMSRGSolver::InverseTransform(Operator& OpIn)
{

  Operator OpOut = OpIn;
  for (auto omega=Omega.rbegin(); omega !=Omega.rend(); ++omega )
  {
    Operator negomega = -(*omega);
    OpOut = OpOut.BCH_Transform( negomega );
  }
  return OpOut;
//   Operator negomega = -Omega;
//   return OpIn.BCH_Transform( negomega );
}



// count number of equations to be solved
int IMSRGSolver::GetSystemDimension()
{
   int dim = 1; // zero-body part

   int N = H_0->OneBody.n_cols;
   dim += N*(N+1)/2;
   dim += H_0->TwoBody.Dimension();
   return dim;
}



void IMSRGSolver::WriteFlowStatus(string fname)
{
   if (fname !="")
   {
     ofstream ff(fname,ios::app);
     WriteFlowStatus(ff);
   }
}
void IMSRGSolver::WriteFlowStatus(ostream& f)
{
   if ( f.good() )
   {
      int fwidth = 16;
      int fprecision = 9;
      f.setf(ios::fixed);
      f << setw(5) << istep
        << setw(fwidth) << setprecision(3) << s
        << setw(fwidth) << setprecision(fprecision) << H_s.ZeroBody 
        << setw(fwidth) << setprecision(fprecision) << H_s.OneBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << H_s.TwoBodyNorm()
//        << setw(fwidth) << setprecision(fprecision) << Omega.Norm()
        << setw(fwidth) << setprecision(fprecision) << Omega.back().Norm()
        << setw(fwidth) << setprecision(fprecision) << Eta.OneBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << Eta.TwoBodyNorm()
        << setw(fwidth) << setprecision(0)          << H_s.timer["N_Commutators"]
        << setprecision(fprecision)
        << endl;
   }

}

void IMSRGSolver::WriteFlowStatusHeader(string fname)
{
   ofstream ff;
   if (fname !="") ff.open(fname,ios::app);
   WriteFlowStatusHeader(ff);
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
        << endl;
      f << "-----------------------------------------------------------------------------------------------------------------------" << endl;
   }

}
