
#ifndef IMSRGSolver_h
#define IMSRGSolver_h 1

#include "Operator.hh"

class IMSRGSolver
{

  public:
  IMSRGSolver(Operator H_in);
  void Solve();
  Operator Transform(Operator& OpIn);

  private:
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

  void UpdateEta();
  void UpdateOmega();
  void UpdateH();

};


#endif

