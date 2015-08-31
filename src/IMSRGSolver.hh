
#ifndef IMSRGSolver_h
#define IMSRGSolver_h 1

#include <fstream>
#include <string>
#include "Operator.hh"
#include "Generator.hh"

using namespace std;


class IMSRGSolver
{

  public:

//  private:
  ModelSpace* modelspace;
  Operator* H_0; 
  vector<Operator> FlowingOps;
//  Operator H_s;
  Operator H_saved;
  Operator Eta;
//  Operator Omega;
  vector<Operator> Omega;
  Generator generator;
  int istep;
  double s;
  double ds;
  double ds_max;
  double smax;
  double norm_domega;
  double omega_norm_max;
  string method;
  string flowfile;


  ~IMSRGSolver();
  IMSRGSolver();
  IMSRGSolver( Operator& H_in);
  void SetHin( Operator& H_in);
  void Reset();
  void AddOperator(Operator& Op){FlowingOps.push_back(Op);};

  void SetMethod(string m){method=m;};
  void Solve();
  void Solve_magnus_euler();
  void Solve_magnus_modified_euler();

  Operator Transform(Operator& OpIn);
  Operator Transform(Operator&& OpIn);
  Operator InverseTransform(Operator& OpIn);
  Operator GetOmega(int i){return Omega[i];};

  void SetFlowFile(string s);
  void SetDs(double d){ds = d;};
  void SetDsmax(double d){ds_max = d;};
  void SetdOmega(double d){norm_domega = d;};
  void SetSmax(double d){smax = d;};
  void SetGenerator(string g);
  void SetOmegaNormMax(double x){omega_norm_max = x;};
  void SetODETolerance(float x){ode_e_abs=x;ode_e_rel=x;};

  int GetSystemDimension();
  Operator& GetH_s(){return FlowingOps[0];};

  void UpdateOmega();
  void UpdateH();

  void WriteFlowStatus(ostream&);
  void WriteFlowStatusHeader(ostream&);
  void WriteFlowStatus(string);
  void WriteFlowStatusHeader(string);

  void SetDenominatorCutoff(double c){generator.SetDenominatorCutoff(c);};
  void SetDenominatorDelta(double d){generator.SetDenominatorDelta(d);};
  void SetDenominatorDeltaIndex(int i){generator.SetDenominatorDeltaIndex(i);};
  void SetDenominatorDeltaOrbit(string o){generator.SetDenominatorDeltaOrbit(o);};



#ifndef NO_ODE

  // This is used to get flow info from odeint
  class ODE_Monitor
  {
    public:
     ODE_Monitor(IMSRGSolver& solver)
           : imsrgsolver(solver), times(solver.times), E0(solver.E0),
              eta1(solver.eta1),eta2(solver.eta2) {};
     IMSRGSolver& imsrgsolver;
     vector<double>& times;
     vector<double>& E0;
     vector<double>& eta1;
     vector<double>& eta2;
     void operator() (const vector<Operator>& x, double t)
     {
        times.push_back(t);
        E0.push_back(x.front().ZeroBody);
        eta1.push_back(imsrgsolver.Eta.OneBodyNorm());
        eta2.push_back(imsrgsolver.Eta.TwoBodyNorm());
     }
     void report()
     {
        for (size_t i=0; i<times.size(); ++i)
        {
           cout << times[i] << "  " << E0[i] << "  "
                << eta1[i]  << "  " << eta2[i] << endl;
        }
     }
  };


  ODE_Monitor ode_monitor;

  vector<double> times;
  vector<double> E0;
  vector<double> eta1;
  vector<double> eta2;
  string ode_mode;
  float ode_e_abs;
  float ode_e_rel;


  void operator()( const vector<Operator>& x, vector<Operator>& dxdt, const double t);
  void Solve_ode();
  void Solve_ode_adaptive();
  void Solve_ode_magnus();

#endif


};





#endif

