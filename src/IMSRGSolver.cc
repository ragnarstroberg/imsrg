
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
    : s(0),ds(0.1),ds_max(0.5),
     norm_domega(0.1), omega_norm_max(2.0),method("magnus_euler"),flowfile("")
#ifndef NO_ODE
    ,ode_monitor(*this),ode_mode("H"),ode_e_abs(1e-6),ode_e_rel(1e-6)
#endif
{}

// Constructor
IMSRGSolver::IMSRGSolver( Operator &H_in)
   : modelspace(H_in.GetModelSpace()),H_0(&H_in), FlowingOps(1,H_in), Eta(H_in), // ,dOmega(H_in)
    istep(0), s(0),ds(0.1),ds_max(0.5),
    smax(2.0), norm_domega(0.1), omega_norm_max(2.0),method("magnus_euler"),flowfile("")
#ifndef NO_ODE
    ,ode_monitor(*this),ode_mode("H"),ode_e_abs(1e-6),ode_e_rel(1e-6)
#endif
{
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
//   H_s = H_in;
   FlowingOps[0] = H_in;
   Eta = H_in;
   Eta.Erase();
   Eta.SetAntiHermitian();
   if (Omega.back().Norm() > 1e-6)
   {
     Omega.push_back(Eta);
//        H_saved = H_s;
        H_saved = FlowingOps[0];
    cout << "pushing back another Omega. Omega.size = " << Omega.size() <<" , operator size = " << Omega.front().Size()*sizeof(double)/1024./1024./1024. << " GB" << endl;
   }
}

void IMSRGSolver::Reset()
{
   s=0;
   Eta.Erase();
//   Omega.Erase();
   Omega.resize(0);
   Omega.push_back(Eta);
    cout << "pushing back another Omega. Omega.size = " << Omega.size() <<" , operator size = " << Omega.front().Size()*sizeof(double)/1024./1024./1024. << " GB" << endl;
}

void IMSRGSolver::SetGenerator(string g)
{
  generator.SetType(g);
  if (Omega.back().Norm() > 1e-6)
  {
    Eta.Erase();
    Omega.push_back(Eta);
//        H_saved = H_s;
        H_saved = FlowingOps[0];
    cout << "pushing back another Omega. Omega.size = " << Omega.size() <<" , operator size = " << Omega.front().Size()*sizeof(double)/1024./1024./1024. << " GB" << endl;
  }
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
  if (method == "magnus_euler" or method =="magnus")
    Solve_magnus_euler();
  else if (method == "magnus_modified_euler")
    Solve_magnus_modified_euler();
  else if (method == "flow_adaptive" or method == "flow")
    Solve_ode_adaptive();
  else if (method == "magnus_adaptive")
    Solve_ode_magnus();
  else if (method == "flow_euler")
    Solve_ode();
  else
    cout << "IMSRGSolver: I don't know method " << method << endl;
}

void IMSRGSolver::Solve_magnus_euler()
{
   istep = 0;
//   generator.Update(&H_s,&Eta);
   generator.Update(&FlowingOps[0],&Eta);

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
        H_saved = FlowingOps[0];
        Omega.push_back(Eta);
        Omega.back().Erase();
        norm_omega = 0;
        cout << "pushing back another Omega. Omega.size = " << Omega.size() <<" , operator size = " << Omega.front().Size()*sizeof(double)/1024./1024./1024. << " GB" << endl;
      }
      // ds should never be more than 1, as this is over-rotating
      //ds = min(norm_domega / norm_eta / (norm_omega+1.0e-9), ds_max); 
      ds = min(min(norm_domega/norm_eta, norm_domega / norm_eta / (norm_omega+1.0e-9)), ds_max); 
//      if (ds == ds_max) norm_domega /=2;
      if (s+ds > smax) ds = smax-s;
      s += ds;
      Eta *= ds; // Here's the Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
      Omega.back() = Eta.BCH_Product( Omega.back() ); 

      // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
//      H_s = H_0->BCH_Transform( Omega );
      if (Omega.size()<2)
      {
        FlowingOps[0] = H_0->BCH_Transform( Omega.back() );
      }
      else
      {
        FlowingOps[0] = H_saved.BCH_Transform( Omega.back() );
      }
        
      generator.Update(&FlowingOps[0],&Eta);

      // Write details of the flow
      WriteFlowStatus(flowfile);
      WriteFlowStatus(cout);

   }

}


void IMSRGSolver::Solve_magnus_modified_euler()
{
   istep = 0;
//   generator.Update(&H_s,&Eta);
   generator.Update(&FlowingOps[0],&Eta);

   Operator H_temp;
    // Write details of the flow
   WriteFlowStatus(flowfile);
   WriteFlowStatus(cout);

   for (istep=1;s<smax;++istep)
   {
      double norm_eta = Eta.Norm();
      double norm_omega = Omega.back().Norm();
      if (norm_omega > omega_norm_max)
      {
        H_saved = FlowingOps[0];
        Omega.push_back(Eta);
        Omega.back().Erase();
        norm_omega = 0;
        cout << "pushing back another Omega. Omega.size = " << Omega.size() <<" , operator size = " << Omega.front().Size()*sizeof(double)/1024./1024./1024. << " GB" << endl;
      }
      // ds should never be more than 1, as this is over-rotating
      ds = min(min(norm_domega/norm_eta, norm_domega / norm_eta / (norm_omega+1.0e-9)), ds_max); 
      if (s+ds > smax) ds = smax-s;
      s += ds;

      H_temp = FlowingOps[0] + ds * Commutator(Eta,FlowingOps[0]);
      generator.AddToEta(&H_temp,&Eta);

      Eta *= ds*0.5; // Here's the modified Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
      Omega.back() = Eta.BCH_Product( Omega.back() ); 

      if (Omega.size()<2)
      {
        FlowingOps[0] = H_0->BCH_Transform( Omega.back() );
      }
      else
      {
        FlowingOps[0] = H_saved.BCH_Transform( Omega.back() );
      }
        
      generator.Update(&FlowingOps[0],&Eta);

      // Write details of the flow
      WriteFlowStatus(flowfile);
      WriteFlowStatus(cout);

   }

}


#ifndef NO_ODE

// Implement element-wise division and abs and reduce for Operators.
// This is required for adaptive steppers
vector<Operator> operator/ (const vector<Operator>& num, const vector<Operator>& denom)
{
   vector<Operator> quotient = num;
   for ( size_t i=0;i<num.size();++i )
   {
     quotient[i].ZeroBody /= denom[i].ZeroBody;
     quotient[i].OneBody /= denom[i].OneBody;
     for ( auto& itmat: quotient[i].TwoBody.MatEl )    itmat.second /= denom[i].TwoBody.GetMatrix(itmat.first[0],itmat.first[1]);
   }
   return quotient;
}

vector<Operator> operator* (const double a, const vector<Operator>& X)
{
  vector<Operator> Y = X;
  for ( auto& y : Y )  y *= a;
  return Y;
}

vector<Operator> operator+ ( const vector<Operator>& X, const vector<Operator>& Y)
{
  vector<Operator> Z = X;
  for ( size_t i=0;i<Z.size();++i )  Z[i] += Y[i];
  return Z;
}

// Also need the dubious operation of adding a double to an operator.
vector<Operator> operator+ (const double a, const vector<Operator>& X)
{
   vector<Operator> Y = X;
   for ( auto& y : Y )
   {
     y.ZeroBody += a;
     for( auto& v : y.OneBody ) v += a;
     for ( auto& itmat: y.TwoBody.MatEl )
       for ( auto& v : itmat.second ) v += a;
   }
   return Y;
}

// Return the element-wise absolute value of an operator
// this is needed for ODE adaptive solver
vector<Operator> abs(const vector<Operator>& OpIn)
{
   vector<Operator> OpOut = OpIn;
   for (auto& opout : OpOut )
   {
     opout.ZeroBody = abs(opout.ZeroBody);
     opout.OneBody = arma::abs(opout.OneBody);
     for ( auto& itmat : opout.TwoBody.MatEl )    itmat.second = arma::abs(itmat.second);
   }
   return OpOut;
}

// Apply operation to each element of X and return the result
// this is needed for ODE adaptive solver
// USE THIS FOR BOOST VERSION < 1.56
/*
namespace boost {namespace numeric {namespace odeint{
template<>
struct vector_space_reduce< Operator >
{
   template<class Op>
   double operator()(const Operator& X, Op op, double init)
   {
      init = op(init,X.ZeroBody);
      for ( auto& v : X.OneBody )    init = op(init,v);
      for ( auto& itmat : X.TwoBody.MatEl )
      {
          for (auto& v : itmat.second )    init = op(init,v);
      }
      return init;
   }
};
}}}
*/

// Apply operation to each element of X and return the result
// this is needed for ODE adaptive solver
// USE THIS FOR BOOST VERSION >= 1.56
namespace boost {namespace numeric {namespace odeint{
template<>
struct vector_space_norm_inf< vector<Operator> >
{
   typedef double result_type;
   double operator()(const vector<Operator>& X)
   {
     double norm = 0;
     for ( auto& x : X )
       norm += x.Norm();
     return norm;
   }
};
}}}


void IMSRGSolver::Solve_ode()
{

   ode_mode = "H";
   WriteFlowStatusHeader(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   runge_kutta4< vector<Operator>, double, vector<Operator>, double, vector_space_algebra> stepper;
   auto system = *this;
   auto monitor = ode_monitor;
//   size_t steps = integrate_const(stepper, system, FlowingOps, s, smax, ds, monitor);
   integrate_const(stepper, system, FlowingOps, s, smax, ds, monitor);
   monitor.report();
}

void IMSRGSolver::Solve_ode_adaptive()
{
   ode_mode = "H";
   WriteFlowStatusHeader(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   auto system = *this;
   typedef runge_kutta_dopri5< vector<Operator> , double , vector<Operator> ,double , vector_space_algebra > stepper;
   auto monitor = ode_monitor;
//   size_t steps = integrate_adaptive(make_controlled<stepper>(ode_e_abs,ode_e_rel), system, FlowingOps, s, smax, ds, monitor);
   integrate_adaptive(make_controlled<stepper>(ode_e_abs,ode_e_rel), system, FlowingOps, s, smax, ds, monitor);
   monitor.report();

}

// Evaluate dx/dt for boost ode
//void IMSRGSolver::operator()( const Operator& x, Operator& dxdt, const double t)
void IMSRGSolver::operator()( const vector<Operator>& x, vector<Operator>& dxdt, const double t)
{
   s = t;
   if (ode_mode == "H")
   {
     FlowingOps[0] = x[0];
     auto& H_s = FlowingOps[0];
     generator.Update(&H_s,&Eta);
     for (size_t i=0;i<x.size();++i)
     {
       dxdt[i] = Commutator(Eta,x[i]);
     }
   }
   else if (ode_mode == "Omega")
   {

     double norm_omega = Omega.back().Norm();
     if (norm_omega > omega_norm_max) // This doesn't seem to works so well just yet...
     {
       H_saved = FlowingOps[0];
       Omega.push_back(Eta);
       Omega.back().Erase();
       norm_omega = 0;
       cout << "pushing back another Omega. Omega.size = " << Omega.size() << endl;
     }
     Omega.back() = x.back();
     auto& Omega_s = x.back();
     Operator& H_s = FlowingOps[0];
     if (Omega.size() > 1)
       H_s = H_saved.BCH_Transform(Omega_s);
     else
       H_s = H_0->BCH_Transform(Omega_s);
     generator.Update(&H_s,&Eta);
     dxdt.back() = Eta - 0.5*Commutator(Omega_s,Eta);
   }
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);
}



void IMSRGSolver::Solve_ode_magnus()
{
   ode_mode = "Omega";
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
   runge_kutta4<vector<Operator>, double, vector<Operator>, double, vector_space_algebra> stepper;
   auto system = *this;
   auto monitor = ode_monitor;
//   size_t steps = integrate_const(stepper, system, Omega, s, smax, ds, monitor);
   integrate_const(stepper, system, Omega, s, smax, ds, monitor);
   monitor.report();
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
      auto& H_s = FlowingOps[0];
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
        << setw(7)      << setprecision(0)          << H_s.timer["N_Commutators"]
        << setw(fwidth) << setprecision(fprecision) << H_s.GetMP2_Energy()
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
        << setw(9)      << setprecision(fprecision) << "Ncomm" 
        << setw(12)     << setprecision(fprecision) << "E(MP2)" 
        << endl;
      f << "---------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   }

}
