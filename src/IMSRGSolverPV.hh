/// This is the inherited IMSRGSolverPV class from IMSRGSolver for parity violating interactions

#ifndef IMSRGSolverPV_h
#define IMSRGSolverPV_h 1

#include "IMSRGSolver.hh"
#include "GeneratorPV.hh"

class IMSRGSolverPV : public IMSRGSolver
{
public:
  Operator *VPT_0;
  Operator VPT_saved;
  Operator Etapv;
  std::string method;
  std::deque<Operator> FlowingOpsPV;
  std::deque<Operator> OmegaPV;
  GeneratorPV generatorPV;
  void AddOperatorPV(Operator& Op){FlowingOpsPV.push_back(Op);};
  Operator GetOperatorPV(size_t i){return FlowingOpsPV.at(i);};
  Operator &GetVPT_s() { return FlowingOpsPV[0]; };

  //  ~IMSRGSolverPV();
  IMSRGSolverPV();
  IMSRGSolverPV(Operator &H_in, Operator &VPT_in);
  void Solve_flow_RK4_PV();
  void Solve_magnus_euler_PV();
  void NewOmega_PV();
  void SetMethod(std::string m) { method = m; };
  void SetGeneratorPV(std::string g);
  void SetFlowFilePV(std::string s);
  void WriteFlowStatusPV(std::ostream &);
  void WriteFlowStatusHeaderPV(std::ostream &);
  void WriteFlowStatusPV(std::string);
  void WriteFlowStatusHeaderPV(std::string);
  Operator &GetEtaPV() { return Etapv; };
  GeneratorPV &GetGeneratorPV() { return generatorPV; };

  std::tuple<Operator, Operator> Transform(Operator &OpIn, Operator &OpInPV);
  std::tuple<Operator, Operator> Transform(Operator &&OpIn, Operator &&OpInPV);
  std::tuple<Operator, Operator> Transform_Partial(Operator &OpIn, Operator &OpInPV, int n);
  std::tuple<Operator, Operator> Transform_Partial(Operator &&OpIn, Operator &&OpInPV, int n);
};

#endif
