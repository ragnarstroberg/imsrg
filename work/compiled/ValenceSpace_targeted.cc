#include <stdlib.h>
#include <sstream>
#include "IMSRG.hh"

using namespace imsrg_util;

int main(int argc, char** argv)
{
  // Default parameters
  string inputtbme = "/itch/exch/BlockGen/me2j/chi2b_srg0800_eMax12_lMax10_hwHO020.me2j.gz";
//  string input3bme = "/itch/exch/BlockGen/me3j/chi2b3b400cD-02cE0098_srg0800ho40C_eMax12_EMax12_hwHO020.me3j.gz";
  string input3bme = "none";
  double hw = 20.0;
  int targetMass = -1;
  double smax = 20.0;
  double dsmax = 0.5;
  double ds_0 = 0.5;
  double domega = 0.5;
  double omega_norm_max = 0.25;
//  double omega_norm_max = 2.5;
  string core_generator = "atan";
  string valence_generator = "shell-model-atan";
  string reference_generator = valence_generator;
  int E3max = 12;
  int eMax = 6;
  string flowfile = "default";
  string intfile = "default";
  string fmt2 = "me2j";
  string valence_space = "sd-shell";
  string basis = "HF";
  string method = "magnus";
  int nsteps = 2;
  double ode_tolerance = 1e-6;
  double R2p = 0.770;
  double R2n = -0.1149;
  double DF = 0.033;
  int file2e1max=12;
  int file2e2max=24;
  int file2lmax=10;
  int file3e1max=12;
  int file3e2max=24;
  int file3e3max=12;
  int denominator_delta=0;
  string reference = "O16";


  // Parse some command line options
  for (int iarg=1; iarg<argc; ++iarg)
  {
    string arg = argv[iarg];
    size_t pos = arg.find("=");
    string var = arg.substr(0,pos);
    string val = arg.substr(pos+1);
    cout << var << " => " << val << endl;

    if (var == "eMax" or var == "emax")
      istringstream(val) >> eMax;
    else if (var == "e3Max" or var == "E3max" or var == "E3Max" or var == "e3max")
      istringstream(val) >> E3max;
    else if (var == "A" or var == "targetMass")
      istringstream(val) >> targetMass;
    else if (var == "2bme")
      istringstream(val) >> inputtbme;
    else if (var == "3bme")
      istringstream(val) >> input3bme;
    else if (var == "hw")
      istringstream(val) >> hw;
    else if (var == "flowfile")
      istringstream(val) >> flowfile;
    else if (var == "intfile")
      istringstream(val) >> intfile;
    else if (var == "valence_space")
      istringstream(val) >> valence_space;
    else if (var == "smax")
      istringstream(val) >> smax;
    else if (var == "fmt2")
      istringstream(val) >> fmt2;
    else if (var == "domega")
      istringstream(val) >> domega;
    else if (var == "omega_norm_max")
      istringstream(val) >> omega_norm_max;
    else if (var == "file2e1max")
      istringstream(val) >> file2e1max;
    else if (var == "file2e2max")
      istringstream(val) >> file2e2max;
    else if (var == "file2lmax")
      istringstream(val) >> file2lmax;
    else if (var == "file3e1max")
      istringstream(val) >> file3e1max;
    else if (var == "file3e2max")
      istringstream(val) >> file3e2max;
    else if (var == "file3e3max")
      istringstream(val) >> file3e3max;
    else if (var == "basis")
      istringstream(val) >> basis;
    else if (var == "method")
      istringstream(val) >> method;
    else if (var == "ode_tolerance")
      istringstream(val) >> ode_tolerance;
    else if (var == "nsteps")
      istringstream(val) >> nsteps;
    else if (var == "core_generator")
      istringstream(val) >> core_generator;
    else if (var == "valence_generator")
      istringstream(val) >> valence_generator;
//    else if (var == "reference_generator")
//      istringstream(val) >> reference_generator;
    else if (var == "denominator_delta")
      istringstream(val) >> denominator_delta;
    else if (var == "reference")
      istringstream(val) >> reference;
    else
      cout << "Unknown parameter: " << var << endl;
  }


  char strbuf[200];
  if (flowfile == "default")
  {
    sprintf(strbuf, "output/BCH_flow_SM_%s_A%d_e%d.dat",valence_space.c_str(),targetMass,eMax);
    flowfile = strbuf;
  }
  if ( intfile == "default")
  {
    sprintf(strbuf, "output/SM_%s_%d_e%d",valence_space.c_str(),targetMass,eMax);
    intfile = strbuf;
  }

  ifstream test(inputtbme);
  if( not test.good() )
  {
    cout << "trouble reading " << inputtbme << " exiting. " << endl;
    return 1;
  }
  test.close();
  if (input3bme != "none")
  {
    test.open(input3bme);
    if( not test.good() )
    {
      cout << "trouble reading " << input3bme << " exiting. " << endl;
      return 1;
    }
    test.close();
  }

  ReadWrite rw;

  vector<string> val,core,ref;
  if (valence_space == "sd-shell")
  {
    val = {"p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
    core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
  }
  else if (valence_space =="fp-shell")
  {
    val = {"p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"};
    core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
  }
  else if (valence_space == "d5-shell")
  {
    val = {"p0d5","n0d5"};
    core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
  }
  else if (valence_space == "d5s1-shell")
  {
    val = {"p0d5","n0d5","p1s1","n1s1"};
    core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
  }
  else if (valence_space == "d5d3-shell")
  {
    val = {"p0d5","n0d5","p0d3","n0d3"};
    core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
  }

  if (reference == "O16")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
  else if (reference == "O22")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","n0d5"};
  else if (reference == "O24")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","n0d5","n1s1"};
  else if (reference == "O28")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","n0d5","n1s1","n0d3"};
  else if (reference == "Si28")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5"};
  else if (reference == "Si34")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","n0d3","n1s1"};
  else if (reference == "S32")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p1s1","n1s1"};
  else if (reference == "C12")
    ref = {"p0s1","n0s1","p0p3","n0p3"};
  else if (reference == "Ca40")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
  else if (reference == "Ca48")
    ref = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1","n0f7"};


  ModelSpace modelspace_target(eMax,ref,val);
//  ModelSpace modelspace_core(eMax,valence_space);
  ModelSpace modelspace_core(eMax,core,val);

  modelspace_target.SetHbarOmega(hw);
  modelspace_core.SetHbarOmega(hw);
  if (targetMass > 0)
  {
    cout << "Setting mass to " << targetMass << endl;
    modelspace_target.SetTargetMass(targetMass);
    modelspace_core.SetTargetMass(targetMass);
  }
  modelspace_target.SetN3max(E3max);
  modelspace_core.SetN3max(E3max);
  
  cout << "Making the operator..." << endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace_target,0,0,0,particle_rank);
  Hbare.SetHermitian();

  cout << "Reading interaction..." << endl;
  if (fmt2 == "me2j")
    rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
  else if (fmt2 == "navratil" or fmt2 == "Navratil")
    rw.ReadBareTBME_Navratil(inputtbme, Hbare);

  if (Hbare.particle_rank >=3)
  {
    rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
  }  

  Hbare.CalculateKineticEnergy();
  
//  Operator Tcm = TCM_Op(modelspace_target);
//  Hbare -= Tcm;
  Hbare -= TCM_Op(modelspace_target);
//  Hbare += Trel_Op(modelspace_target);
  
  HartreeFock hf(Hbare);
  hf.Solve();
  cout << "EHF = " << hf.EHF << endl;
  
  if (basis == "HF")
    Hbare = hf.GetNormalOrderedH();
  else if (basis == "oscillator")
    Hbare = Hbare.DoNormalOrdering();

   if (method == "NSmagnus")
  {
    omega_norm_max=500;
    method = "magnus";
  } 

  IMSRGSolver imsrgsolver(Hbare);
  
  imsrgsolver.SetMethod(method);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDsmax(dsmax);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(min(domega,omega_norm_max+1e-6));
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);

  if (nsteps > 1)
  {
     imsrgsolver.SetGenerator(core_generator);
     imsrgsolver.Solve();
     if (method == "magnus")
       smax *= 2;
  }

  imsrgsolver.SetGenerator(valence_generator);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.Solve();

  Operator Hprime = imsrgsolver.GetH_s().UndoNormalOrdering();
  Hprime.SetModelSpace(modelspace_core);
  Hprime = Hprime.DoNormalOrdering();

  IMSRGSolver imsrgsolver2(Hprime);

  imsrgsolver2.SetMethod(method);
  imsrgsolver2.SetGenerator(reference_generator);
  imsrgsolver2.SetSmax(smax);
  imsrgsolver2.SetFlowFile(flowfile+"_2");
  imsrgsolver2.SetDsmax(dsmax);
  imsrgsolver2.SetDenominatorDelta(denominator_delta);
  imsrgsolver2.SetDs(ds_0);
  imsrgsolver2.SetdOmega(min(domega,omega_norm_max+1e-6));
  imsrgsolver2.SetOmegaNormMax(omega_norm_max);
  imsrgsolver2.SetODETolerance(ode_tolerance);

  imsrgsolver2.Solve();

  rw.WriteNuShellX_int(imsrgsolver2.GetH_s(),intfile+".int");
  rw.WriteNuShellX_sps(imsrgsolver2.GetH_s(),intfile+".sp");


  if (method == "magnus")
  {
     vector<Operator> oplist;
     oplist.push_back(R2_p1_Op(modelspace_target));
     oplist.push_back(R2_p2_Op(modelspace_target));
     oplist.push_back(R2CM_Op(modelspace_target));

     for (Operator& op : oplist )
     {
       if (basis == "HF")
          op = hf.TransformToHFBasis(op);
       op = op.DoNormalOrdering();
       op = imsrgsolver.Transform(op);
       op = op.UndoNormalOrdering();
       op.SetModelSpace(modelspace_core);
       op = op.DoNormalOrdering();
       op = imsrgsolver2.Transform(op);
     }

     rw.WriteNuShellX_op(oplist[0],intfile+"_R2p1.int");
     rw.WriteNuShellX_op(oplist[1],intfile+"_R2p2.int");
     rw.WriteNuShellX_op(oplist[2],intfile+"_R2cm.int");
  }
  


  Hbare.PrintTimes();
 
  return 0;
}

