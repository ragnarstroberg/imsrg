
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

  void UpdateEta();
  void UpdateOmega();
  void UpdateH();

  void WriteFlowStatus(ostream&);

  void ConstructGenerator_Wegner();
  void ConstructGenerator_White();
  void ConstructGenerator_ShellModel();
  void ConstructGenerator_ShellModel1hw(); // Doesn't work yet
  double GetEpsteinNesbet1bDenominator(int i, int j);
  double GetEpsteinNesbet2bDenominator(int ch, int ibra, int iket);

};


#endif

