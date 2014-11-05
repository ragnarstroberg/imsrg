
#ifndef IMSRGSolver_h
#define IMSRGSolver_h 1

#include "Operator.hh"
#include <fstream>
#include <string>

using namespace std;

class IMSRGSolver
{

  public:
  IMSRGSolver(Operator H_in);
  void Solve();
  Operator Transform(Operator& OpIn);
  void SetFlowFile(string s){flowfile = s;};
  void SetDs(double d){ds = d;};
  void SetSmax(double d){smax = d;};

  private:
  ModelSpace* modelspace;
  Operator Omega;
  Operator dOmega;
  Operator H_0;
  Operator H_s;
  Operator H_od;
  Operator H_diag;
  Operator Eta;
  double s;
  double ds;
  double smax;
  string method;
  string generator;
  string flowfile;

  void UpdateEta();
  void UpdateOmega();
  void UpdateH();

};


#endif

