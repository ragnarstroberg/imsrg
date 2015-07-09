
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
  Operator H_s;
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

  void Solve();

  Operator Transform(Operator& OpIn);
  Operator Transform(Operator&& OpIn);
  Operator InverseTransform(Operator& OpIn);
  Operator GetOmega(int i){return Omega[i];};

  void SetFlowFile(string s);
  void SetDs(double d){ds = d;};
  void SetDsmax(double d){ds_max = d;};
  void SetdOmega(double d){norm_domega = d;};
  void SetSmax(double d){smax = d;};
  void SetGenerator(string g){generator.SetType(g);};
  void SetOmegaNormMax(double x){omega_norm_max = x;};

  int GetSystemDimension();

  void UpdateOmega();
  void UpdateH();

  void WriteFlowStatus(ostream&);
  void WriteFlowStatusHeader(ostream&);
  void WriteFlowStatus(string);
  void WriteFlowStatusHeader(string);

  void SetDenominatorCutoff(double c){generator.SetDenominatorCutoff(c);};



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
     void operator() (const Operator& x, double t)
     {
        times.push_back(t);
        E0.push_back(x.ZeroBody);
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


  void Solve_ode();
  void Solve_ode_adaptive();
  void Solve_ode_magnus();
  void ODE_systemH( Operator& x, Operator& dxdt, const double t);
  void ODE_systemOmega( Operator& x, Operator& dxdt, const double t);

#endif


};





#endif

