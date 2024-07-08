///This is the inherited GeneratorPV class from Generator for parity violating interactions


#ifndef GeneratorPV_hh
#define GeneratorPV_hh 1

#include "Generator.hh"


class GeneratorPV : public Generator
{ 
 private:
 Operator * V;
 Operator * Etapv;
  
 public:
 void Update(Operator& H_s, Operator& HPV_s, Operator& Eta_s, Operator& EtaPV_s);
 void AddToEtaPV(Operator& H_s,Operator& HPV_s,Operator& Eta_s,Operator& EtaPV_s); 
 void ConstructGeneratorPV_SingleRef(std::function<double (double,double)>& etafunc );
 void ConstructGeneratorPV_ShellModel(std::function<double (double,double)>& eta_func);
// void ConstructGeneratorPV_ShellModel_NpNh(std::function<double(double,double)>& eta_func);
 //void ConstructGeneratorPV_HartreeFock();
// void ConstructGeneratorPV_1PA(std::function<double(double,double)>& eta_func);
// double Get1bDenominator(int i, int j);

};

#endif

