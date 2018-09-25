
#include "IMSRGSolver.hh"
#include <iomanip>

#ifndef NO_ODE
#include <boost/numeric/odeint.hpp>
#endif

IMSRGSolver::~IMSRGSolver()
{
  CleanupScratch();
}

IMSRGSolver::IMSRGSolver()
    : rw(NULL),s(0),ds(0.1),ds_max(0.5),
     norm_domega(0.1), omega_norm_max(2.0),eta_criterion(1e-6),method("magnus_euler"),
     flowfile(""), n_omega_written(0),max_omega_written(50),magnus_adaptive(true)
     ,ode_monitor(*this),ode_mode("H"),ode_e_abs(1e-6),ode_e_rel(1e-6)
{}

// Constructor
IMSRGSolver::IMSRGSolver( Operator &H_in)
   : modelspace(H_in.GetModelSpace()),rw(NULL), H_0(&H_in), FlowingOps(1,H_in), Eta(H_in), 
    istep(0), s(0),ds(0.1),ds_max(0.5),
    smax(2.0), norm_domega(0.1), omega_norm_max(2.0),eta_criterion(1e-6),method("magnus_euler"),
    flowfile(""), n_omega_written(0),max_omega_written(50),magnus_adaptive(true)
    ,ode_monitor(*this),ode_mode("H"),ode_e_abs(1e-6),ode_e_rel(1e-6)
{
   Eta.Erase();
   Eta.SetAntiHermitian();
//   Omega.push_back( Eta);
   Omega.emplace_back( Eta);
}


void IMSRGSolver::NewOmega()
{
  H_saved = FlowingOps[0];
  cout << "pushing back another Omega. Omega.size = " << Omega.size()
       << " , operator size = " << Omega.front().Size()/1024./1024. << " MB"
       << ",  memory usage = " << profiler.CheckMem()["RSS"]/1024./1024. << " GB"
       << endl;
  if ((rw != NULL) and (rw->GetScratchDir() !=""))
  {
    
    char tmp[512];
    sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), n_omega_written);
    string fname(tmp);
    ofstream ofs(fname, ios::binary);
    Omega.back().WriteBinary(ofs);
    if (Omega.back().GetModelSpace() != Eta.GetModelSpace()) Omega.back() = Eta;
    n_omega_written++;
    cout << "Omega written to file " << fname << "  written " << n_omega_written << " so far." << endl;
    if (n_omega_written > max_omega_written)
    {
      cout << "n_omega_written > max_omega_written.  (" << n_omega_written << " > " << max_omega_written
           << " ) deleting OMEGA files and calling terminate." << endl;
      CleanupScratch();
      terminate();
    }
  }
  else
  {
    Omega.emplace_back(Eta);
  }
  Omega.back().Erase();

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
    NewOmega();
   }
   else
   {
     Omega.back() = Eta;
   }
}



void IMSRGSolver::SetOmega(size_t i, Operator& om)
{
 if ((i+1)> Omega.size())
 {
  Omega.resize(i+1);
 }
 Omega[i] = om;
}

void IMSRGSolver::Reset()
{
   s=0;
   Eta.Erase();
   Omega.resize(0);
   NewOmega();
}

void IMSRGSolver::SetGenerator(string gen)
{
  generator.SetType(gen);
  if (Omega.back().Norm() > 1e-6)
  {
    Eta.Erase();
    NewOmega();
  }
}

void IMSRGSolver::SetFlowFile(string str)
{
   flowfile = str;
   ofstream flowf;
   if (flowfile != "")
   {
      flowf.open(flowfile,ofstream::out);
      flowf.close();
   }
}


void IMSRGSolver::Solve()
{
  if (s<1e-4)
   WriteFlowStatusHeader(cout);


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
  else if (method == "restore_4th_order")
  {
    FlowingOps.emplace_back( Operator( *(FlowingOps[0].GetModelSpace()), 0,0,0,1));
    Solve_ode_adaptive();
  }
  else
    cout << "IMSRGSolver: I don't know method " << method << endl;
}

void IMSRGSolver::UpdateEta()
{
   generator.Update(&FlowingOps[0],&Eta);
}

void IMSRGSolver::Solve_magnus_euler()
{
   istep = 0;
   generator.Update(&FlowingOps[0],&Eta);

   if (generator.GetType() == "shell-model-atan")
   {
     generator.SetDenominatorCutoff(1.0); // do we need this?
   }

    // Write details of the flow
   WriteFlowStatus(flowfile);
   WriteFlowStatus(cout);

   for (istep=1;s<smax;++istep)
   {

      double norm_eta = Eta.Norm();
      if (norm_eta < eta_criterion )
      {
        break;
      }
      double norm_omega = Omega.back().Norm();
      if (norm_omega > omega_norm_max)
      {
        NewOmega();
        norm_omega = 0;
      }
      // ds should never be more than 1, as this is over-rotating
      if (magnus_adaptive)
         ds = min( min( min(norm_domega/norm_eta, norm_domega / norm_eta / (norm_omega+1.0e-9)), omega_norm_max/norm_eta), ds_max); 
      ds = min(ds,smax-s);
//      if (s+ds > smax) ds = smax-s;
      s += ds;
      Eta *= ds; // Here's the Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
      Omega.back() = Eta.BCH_Product( Omega.back() ); 

      // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
      if ((Omega.size()+n_omega_written)<2)
      {
        FlowingOps[0] = H_0->BCH_Transform( Omega.back() );
      }
      else
      {
        FlowingOps[0] = H_saved.BCH_Transform( Omega.back() );
      }

      if (norm_eta<1.0 and generator.GetType() == "shell-model-atan")
      {
        generator.SetDenominatorCutoff(1e-6);
      }
        
      generator.Update(&FlowingOps[0],&Eta);

      // Write details of the flow
      WriteFlowStatus(flowfile);
      WriteFlowStatus(cout);
//      profiler.PrintMemory();

   }

}


void IMSRGSolver::Solve_magnus_modified_euler()
{
   istep = 0;
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
        NewOmega();
        norm_omega = 0;
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

      if ((Omega.size()+n_omega_written)<2)
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
//vector<Operator> operator/ (const vector<Operator>& num, const vector<Operator>& denom)
deque<Operator> operator/ (const deque<Operator>& num, const deque<Operator>& denom)
{
//   vector<Operator> quotient = num;
   deque<Operator> quotient = num;
   for ( size_t i=0;i<num.size();++i )
   {
     quotient[i].ZeroBody /= denom[i].ZeroBody;
     quotient[i].OneBody /= denom[i].OneBody;
     for ( auto& itmat: quotient[i].TwoBody.MatEl )    itmat.second /= denom[i].TwoBody.GetMatrix(itmat.first[0],itmat.first[1]);
   }
   return quotient;
}

//vector<Operator> operator* (const double a, const vector<Operator>& X)
deque<Operator> operator* (const double a, const deque<Operator>& X)
{
//  vector<Operator> Y = X;
  deque<Operator> Y = X;
  for ( auto& y : Y )  y *= a;
  return Y;
}

//vector<Operator> operator+ ( const vector<Operator>& X, const vector<Operator>& Y)
deque<Operator> operator+ ( const deque<Operator>& X, const deque<Operator>& Y)
{
//  vector<Operator> Z = X;
  deque<Operator> Z = X;
  for ( size_t i=0;i<Z.size();++i )  Z[i] += Y[i];
  return Z;
}

// Also need the dubious operation of adding a double to an operator.
//vector<Operator> operator+ (const double a, const vector<Operator>& X)
deque<Operator> operator+ (const double a, const deque<Operator>& X)
{
//   vector<Operator> Y = X;
   deque<Operator> Y = X;
   for ( auto& y : Y )
   {
     y.ZeroBody += a;
     y.OneBody += a;
//     for( auto& v : y.OneBody ) v += a;
     for ( auto& itmat: y.TwoBody.MatEl )
      itmat.second += a;
//       for ( auto& v : itmat.second ) v += a;
   }
   return Y;
}

// Return the element-wise absolute value of an operator
// this is needed for ODE adaptive solver
//vector<Operator> abs(const vector<Operator>& OpIn)
deque<Operator> abs(const deque<Operator>& OpIn)
{
//   vector<Operator> OpOut = OpIn;
   deque<Operator> OpOut = OpIn;
   for (auto& opout : OpOut )
   {
     opout.ZeroBody = std::abs(opout.ZeroBody);
     opout.OneBody = arma::abs(opout.OneBody);
     for ( auto& itmat : opout.TwoBody.MatEl )    itmat.second = arma::abs(itmat.second);
   }
   return OpOut;
}

// Apply operation to each element of X and return the result
// this is needed for ODE adaptive solver
// USE THIS FOR BOOST VERSION < 1.56
#ifdef OLD_BOOST
namespace boost {namespace numeric {namespace odeint{
template<>
struct vector_space_reduce< deque<Operator> >
{
   template<class Op>
   double operator()(const deque<Operator>& X, Op op, double init)
   {
      for (auto& x : X)
      {
        init = op(init,x.ZeroBody);
        for ( auto& v : x.OneBody )    init = op(init,v);
        for ( auto& itmat : x.TwoBody.MatEl )
        {
            for (auto& v : itmat.second )    init = op(init,v);
        }
      }
      return init;
   }
};
}}}
#endif

// Apply operation to each element of X and return the result
// this is needed for ODE adaptive solver
// USE THIS FOR BOOST VERSION >= 1.56
//struct vector_space_norm_inf< vector<Operator> >
#ifndef OLD_BOOST
namespace boost {namespace numeric {namespace odeint{
template<>
struct vector_space_norm_inf< deque<Operator> >
{
   typedef double result_type;
//   double operator()(const vector<Operator>& X)
   double operator()(const deque<Operator>& X)
   {
     double norm = 0;
     for ( auto& x : X )
       norm += x.Norm();
     return norm;
   }
};
}}}
#endif

void IMSRGSolver::Solve_ode()
{

   ode_mode = "H";
   WriteFlowStatusHeader(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
//   runge_kutta4< vector<Operator>, double, vector<Operator>, double, vector_space_algebra> stepper;
   runge_kutta4< deque<Operator>, double, deque<Operator>, double, vector_space_algebra> stepper;
   auto system = *this;
   auto monitor = ode_monitor;
//   size_t steps = integrate_const(stepper, system, FlowingOps, s, smax, ds, monitor);
   integrate_const(stepper, system, FlowingOps, s, smax, ds, monitor);
   monitor.report();
}

void IMSRGSolver::Solve_ode_adaptive()
{
   ode_mode = "H";
   if (method == "restore_4th_order") ode_mode = "Restored";
   WriteFlowStatusHeader(cout);
   WriteFlowStatus(flowfile);
   cout << "done writing header and status" << endl;
   using namespace boost::numeric::odeint;
   auto system = *this;
//   typedef runge_kutta_dopri5< vector<Operator> , double , vector<Operator> ,double , vector_space_algebra > stepper;
   typedef runge_kutta_dopri5< deque<Operator> , double , deque<Operator> ,double , vector_space_algebra > stepper;
//   typedef adams_bashforth_moulton< 4, vector<Operator> , double , vector<Operator> ,double , vector_space_algebra > stepper;
   auto monitor = ode_monitor;
//   size_t steps = integrate_adaptive(make_controlled<stepper>(ode_e_abs,ode_e_rel), system, FlowingOps, s, smax, ds, monitor);
   integrate_adaptive(make_controlled<stepper>(ode_e_abs,ode_e_rel), system, FlowingOps, s, smax, ds, monitor);
   monitor.report();

}

// Evaluate dx/dt for boost ode
//void IMSRGSolver::operator()( const Operator& x, Operator& dxdt, const double t)
//void IMSRGSolver::operator()( const vector<Operator>& x, vector<Operator>& dxdt, const double t)
void IMSRGSolver::operator()( const deque<Operator>& x, deque<Operator>& dxdt, const double t)
{
   s = t;
   if (ode_mode == "H")
   {
     FlowingOps[0] = x[0];
     if (dxdt.size() < x.size()) dxdt.resize(x.size());
     auto& H_s = FlowingOps[0];
     generator.Update(&H_s,&Eta);
     if (Eta.Norm() < eta_criterion)
     {
       for (size_t i=0;i<x.size();++i)
       {
         dxdt[i] = 0*x[i];
       }
     }
     else
     {
       for (size_t i=0;i<x.size();++i)
       {
         dxdt[i] = Commutator(Eta,x[i]);
       }
     }
   }
   else if (ode_mode == "Omega")
   {

     double norm_omega = Omega.back().Norm();
     if (norm_omega > omega_norm_max) // This doesn't seem to works so well just yet...
     {
       NewOmega();
       norm_omega = 0;
     }
     Omega.back() = x.back();
     auto& Omega_s = x.back();
     Operator& H_s = FlowingOps[0];
     if ((Omega.size()+n_omega_written) > 1)
       H_s = H_saved.BCH_Transform(Omega_s);
     else
       H_s = H_0->BCH_Transform(Omega_s);
     generator.Update(&H_s,&Eta);
     if (dxdt.size() < x.size()) dxdt.resize(x.size());
     dxdt.back() = Eta - 0.5*Commutator(Omega_s,Eta);
   }
   else if (ode_mode == "Restored" )
   {
     FlowingOps[0] = x[0];
     FlowingOps[1] = x[1];
     if (dxdt.size() < x.size()) dxdt.resize(x.size());
     dxdt[1] = Operator(x[1]);
     auto& H_s = FlowingOps[0];
     generator.Update(&H_s,&Eta);
     if (Eta.Norm() < eta_criterion)
     {
       for (size_t i=0;i<x.size();++i)
       {
         dxdt[i] = 0*x[i];
       }
     }
     else
     {
       dxdt[0] = Commutator(Eta,x[0]+x[1]);
       dxdt[1].Erase();
       dxdt[1].comm221ss(Eta,x[0]);
       // keep only pp and hh parts of d chi/ ds
       for (auto& a : modelspace->holes)
       {
         for (auto& i : modelspace->particles)
         {
           dxdt[1].OneBody(a,i) = 0;
           dxdt[1].OneBody(i,a) = 0;
         }
       }
       for (size_t i=2;i<x.size();++i)
       {
         dxdt[i] = Commutator(Eta,x[i]);
       }
     }

   }
   WriteFlowStatus(flowfile);
   WriteFlowStatus(cout);
}



void IMSRGSolver::Solve_ode_magnus()
{
   ode_mode = "Omega";
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
//   runge_kutta4<vector<Operator>, double, vector<Operator>, double, vector_space_algebra> stepper;
   runge_kutta4<deque<Operator>, double, deque<Operator>, double, vector_space_algebra> stepper;
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
  return Transform_Partial(OpIn, 0);
}

Operator IMSRGSolver::Transform(Operator&& OpIn)
{
  return Transform_Partial(OpIn, 0);
}




/// Returns \f$ e^{-Omega} \mathcal{O} e^{Omega} \f$
Operator IMSRGSolver::InverseTransform(Operator& OpIn)
{
//  if (OpIn.GetJRank()+OpIn.GetTRank()+OpIn.GetParity()>0)
//  {
//    OpIn.ResetTensorTransformFirstPass();
//  }
  Operator OpOut = OpIn;
  for (auto omega=Omega.rbegin(); omega !=Omega.rend(); ++omega )
  {
    Operator negomega = -(*omega);
    OpOut = OpOut.BCH_Transform( negomega );
  }
  return OpOut;
}

/// Returns \f$ e^{\Omega} \mathcal{O} e^{-\Omega} \f$
/// for the \f$\Omega_i\f$s with index greater than or equal to n.
Operator IMSRGSolver::Transform_Partial(Operator& OpIn, int n)
{
//  cout << "Begin Transform_Partial" << endl;
  Operator OpOut = OpIn;
  if ((rw != NULL) and rw->GetScratchDir() != "")
  {
    Operator omega(OpIn);
    char tmp[512];
    for (int i=n;i<n_omega_written;i++)
    {
     sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
     string fname(tmp);
     ifstream ifs(fname,ios::binary);
     omega.ReadBinary(ifs);
//     if (OpIn.GetJRank()>0) cout << "step " << i << endl;
     OpOut = OpOut.BCH_Transform( omega );
//     if (OpIn.GetJRank()>0)cout << "done" << endl;
    }
  }

  for (size_t i=max(n-n_omega_written,0); i<Omega.size();++i)
  {
//     if (OpIn.GetJRank()>0) cout << "step " << i << endl;
    OpOut = OpOut.BCH_Transform( Omega[i] );
//     if (OpIn.GetJRank()>0)cout << "done" << endl;
  }

  return OpOut;
}


Operator IMSRGSolver::Transform_Partial(Operator&& OpIn, int n)
{
//  cout << "Calling r-value version of Transform_Partial, n = " << n << endl;
  Operator OpOut = OpIn;
  if ((rw != NULL) and rw->GetScratchDir() != "")
  {
    Operator omega(OpIn);
    char tmp[512];
    for (int i=n;i<n_omega_written;i++)
    {
     sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
     string fname(tmp);
     ifstream ifs(fname,ios::binary);
     omega.ReadBinary(ifs);
     OpOut = OpOut.BCH_Transform( omega );
    }
  }

  for (size_t i=max(n-n_omega_written,0); i<Omega.size();++i)
  {
    OpOut = OpOut.BCH_Transform( Omega[i] );
  }
  return OpOut;
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



void IMSRGSolver::CleanupScratch()
{
  if (n_omega_written<=0) return;
  cout << "Cleaning up files written to scratch space" << endl;
  char tmp[512];
  for (int i=0;i<n_omega_written;i++)
  {
    sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
    string fname(tmp);
    if ( remove(tmp) !=0 )
    {
      cout << "Error when attempting to delete " << fname << endl;
    }
  }
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
      f << fixed << setw(5) << istep
        << setw(10) << setprecision(3) << s
        << setw(fwidth) << setprecision(fprecision) << H_s.ZeroBody 
        << setw(fwidth) << setprecision(fprecision) << H_s.Norm()
        << setw(fwidth) << setprecision(fprecision) << H_s.Trace( modelspace->GetAref(), modelspace->GetZref() )
//        << setw(fwidth) << setprecision(fprecision) << H_s.OneBodyNorm()
//        << setw(fwidth) << setprecision(fprecision) << H_s.TwoBodyNorm()
//        << setw(fwidth) << setprecision(fprecision) << Omega.Norm()
        << setw(fwidth) << setprecision(fprecision) << Omega.back().OneBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << Omega.back().TwoBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << Eta.OneBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << Eta.TwoBodyNorm()
        << setw(7)      << setprecision(0)          << profiler.counter["N_ScalarCommutators"] + profiler.counter["N_TensorCommutators"]
        << setw(fwidth) << setprecision(fprecision) << H_s.GetMP2_Energy()
        << setw(7)      << setprecision(0)          << profiler.counter["N_Operators"]
        << setprecision(fprecision)
        << setw(12) << setprecision(3) << profiler.GetTimes()["real"]
        << setw(12) << setprecision(3) << profiler.CheckMem()["RSS"]/1024. << " / " << skipws << profiler.MaxMemUsage()/1024. << fixed
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
      f << fixed << setw(5) << "i"
        << setw(10) << setprecision(3) << "s"
        << setw(fwidth) << setprecision(fprecision) << "E0"
//        << setw(fwidth) << setprecision(fprecision) << "||H_1||" 
//        << setw(fwidth) << setprecision(fprecision) << "||H_2||" 
        << setw(fwidth) << setprecision(fprecision) << "||H||" 
        << setw(fwidth) << setprecision(fprecision) << "Tr(H)/Tr(1)" 
        << setw(fwidth) << setprecision(fprecision) << "||Omega_1||" 
        << setw(fwidth) << setprecision(fprecision) << "||Omega_2||" 
        << setw(fwidth) << setprecision(fprecision) << "||Eta_1||" 
        << setw(fwidth) << setprecision(fprecision) << "||Eta_2||" 
        << setw(7)      << setprecision(fprecision) << "Ncomm" 
        << setw(16)     << setprecision(fprecision) << "E(MP2)" 
        << setw(7)      << setprecision(fprecision) << "N_Ops"
        << setw(16) << setprecision(fprecision) << "Walltime (s)"
        << setw(19) << setprecision(fprecision) << "Memory (MB)"
        << endl;
        for (int x=0;x<175;x++) f << "-";
        f << endl;
   }

}
