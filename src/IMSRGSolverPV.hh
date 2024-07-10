///This is the inherited IMSRGSolverPV class from IMSRGSolver for parity violating interactions

#ifndef IMSRGSolverPV_h
#define IMSRGSolverPV_h 1

#include "IMSRGSolver.hh"
#include "GeneratorPV.hh"


class IMSRGSolverPV: public IMSRGSolver
{
 public:
  Operator* VPT_0;
  Operator VPT_saved;
  Operator* Schiff_0;
  Operator Schiff_saved;
  Operator* Schiffpp_0;
  Operator Schiffpp_saved;
  Operator Etapv;
  std::deque<Operator> FlowingOpsH;
  std::deque<Operator> FlowingOpsV;
  std::deque<Operator> FlowingOpsSchiff;
  std::deque<Operator> FlowingOpsSchiffpp;
  GeneratorPV generatorPV;
  Operator& GetH_s(){return FlowingOpsH[0];};
  Operator& GetVPT_s(){return FlowingOpsV[0];};
  Operator& GetSchiff_s(){return FlowingOpsSchiff[0];};
  Operator& GetSchiffpp_s(){return FlowingOpsSchiffpp[0];};


//  ~IMSRGSolverPV();
  IMSRGSolverPV();
  IMSRGSolverPV(Operator& H_in, Operator& VPT_in, Operator& Schiff_in, Operator& Schifpp_in);
  void Solve_flow_RK4_PV();
  void SetGeneratorPV(std::string g);
};

#endif  

