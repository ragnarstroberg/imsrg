
#ifndef IMSRGSolver_h
#define IMSRGSolver_h 1

#include "Operator.hh"
#include <fstream>
#include <string>

using namespace std;


class IMSRGSolver
{

#ifndef NO_ODE
  // This is used to get flow info from odeint
  class ODE_Monitor
  {
    public:
     ODE_Monitor(IMSRGSolver& solver)
           : imsrgsolver(solver), times(solver.times), E0(solver.E0),
              eta1(solver.eta1),eta2(solver.eta2) {};
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
    private:
     IMSRGSolver& imsrgsolver;
  };
#endif

  public:
  ~IMSRGSolver();
  IMSRGSolver();
  IMSRGSolver( Operator& H_in);
  void SetHin( Operator& H_in);
  void Reset();

  void Solve();
#ifndef NO_ODE
  void Solve_ode();
  void Solve_ode_adaptive();
  void Solve_ode_magnus();
  void ODE_systemH( Operator& x, Operator& dxdt, const double t);
  void ODE_systemOmega( Operator& x, Operator& dxdt, const double t);
#endif
  Operator Transform(Operator& OpIn);
  Operator InverseTransform(Operator& OpIn);
  void SetFlowFile(string s){flowfile = s;};
  void SetDs(double d){ds = d;};
  void SetDsmax(double d){ds_max = d;};
  void SetdOmega(double d){norm_domega = d;};
  void SetSmax(double d){smax = d;};
  void SetGenerator(string g){generator = g;};
  int GetSystemDimension();

  void UpdateEta();
  void UpdateOmega();
  void UpdateH();

  void WriteFlowStatus(ostream&);
  void WriteFlowStatusHeader(ostream&);
  void WriteFlowStatus(string);
  void WriteFlowStatusHeader(string);

  void ConstructGenerator_Wegner();
  void ConstructGenerator_White();
  void ConstructGenerator_Atan();
  void ConstructGenerator_ShellModel();
  void ConstructGenerator_ShellModel_Atan();
  void ConstructGenerator_HartreeFock();
  double Get1bDenominator_ph(int i, int j);
  double Get1bDenominator_pp(int i, int j);
  double Get2bDenominator_pphh(int ch, int ibra, int iket);
  double Get2bDenominator_pppp(int ch, int ibra, int iket);
  double Get2bDenominator_pphp(int ch, int ibra, int iket);


//  private:
  ModelSpace* modelspace;
<<<<<<< HEAD
  Operator Omega;
  Operator H_0;
=======
//  Operator H_0; //Having this copy is dumb.
  Operator* H_0; //Having this copy is dumb.
>>>>>>> master
  Operator H_s;
  Operator Eta;
  Operator Omega;
  int istep;
  double s;
  double ds;
  double ds_max;
  double smax;
  double norm_domega;
  string method;
  string generator;
  string flowfile;

#ifndef NO_ODE
  ODE_Monitor ode_monitor;

  vector<double> times;
  vector<double> E0;
  vector<double> eta1;
  vector<double> eta2;
#endif


};





#endif

