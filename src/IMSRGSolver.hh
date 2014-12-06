
#ifndef IMSRGSolver_h
#define IMSRGSolver_h 1

#include "Operator.hh"
#include <fstream>
#include <string>

using namespace std;


class IMSRGSolver
{

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


  public:
  IMSRGSolver(const Operator& H_in);
  void Solve();
  void Solve_ode();
  void ODE_systemH(const Operator& x, Operator& dxdt, const double t);
  Operator Transform(Operator& OpIn);
  void SetFlowFile(string s){flowfile = s;};
  void SetDs(double d){ds = d;};
  void SetdOmega(double d){norm_domega = d;};
  void SetSmax(double d){smax = d;};
  void SetGenerator(string g){generator = g;};

//  private:
  ModelSpace* modelspace;
  Operator Omega;
  Operator dOmega;
  Operator H_0;
  Operator H_s;
  Operator H_od;
  Operator H_diag;
  Operator Eta;
  int istep;
  int i_full_BCH;
  double s;
  double ds;
  double ds_max;
  double smax;
  double norm_domega;
  string method;
  string generator;
  string flowfile;
  ODE_Monitor ode_monitor;


  vector<double> times;
  vector<double> E0;
  vector<double> eta1;
  vector<double> eta2;


  void UpdateEta();
  void UpdateOmega();
  void UpdateH();

  void WriteFlowStatus(ostream&);

  void ConstructGenerator_Wegner();
  void ConstructGenerator_White();
  void ConstructGenerator_Atan();
  void ConstructGenerator_ShellModel();
  void ConstructGenerator_ShellModel_Atan();
  double GetEpsteinNesbet1bDenominator(int i, int j);
  double GetEpsteinNesbet2bDenominator(int ch, int ibra, int iket);





};





#endif

