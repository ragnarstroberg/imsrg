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
  int targetMass = 18;
  double smax = 20.0;
  double dsmax = 0.5;
  double ds_0 = 0.5;
  double domega = 0.5;
  double omega_norm_max = 0.25;
//  double omega_norm_max = 2.5;
  string core_generator = "atan";
  string valence_generator = "shell-model-atan";
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
    else if (var == "denominator_delta")
      istringstream(val) >> denominator_delta;
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
  test.open(input3bme);
  if( not test.good() )
  {
    cout << "trouble reading " << input3bme << " exiting. " << endl;
    return 1;
  }
  test.close();

  ReadWrite rw;



  ModelSpace modelspace(eMax,valence_space);

  modelspace.SetHbarOmega(hw);
//  modelspace.SetTargetMass(targetMass);
  modelspace.SetN3max(E3max);
  
  cout << "Making the operator..." << endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
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
  
  Operator Tcm = TCM_Op(modelspace);
  Hbare -= Tcm;
  
  HartreeFock hf(Hbare);
  hf.Solve();
  cout << "EHF = " << hf.EHF << endl;
  
  if (basis == "HF")
    Hbare = hf.GetNormalOrderedH();
  else if (basis == "oscillator")
    Hbare = Hbare.DoNormalOrdering();

//  Operator Rp2 = Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ());
  Operator R2_p1 = R2_p1_Op(modelspace);
  Operator R2_p2 = R2_p2_Op(modelspace);
  Operator R2_cm  = R2CM_Op(modelspace);
  if (basis == "HF")
  {
    R2_p1 = hf.TransformToHFBasis(R2_p1);
    R2_p2 = hf.TransformToHFBasis(R2_p2);
    R2_cm = hf.TransformToHFBasis(R2_cm);
  }
  R2_p1 = R2_p1.DoNormalOrdering();
  R2_p2 = R2_p2.DoNormalOrdering();
  R2_cm = R2_cm.DoNormalOrdering();
//  cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << endl; 
//  cout << " HF charge radius = " << sqrt( Rp2.ZeroBody + R2p + R2n + DF) << endl; 
  
  IMSRGSolver imsrgsolver(Hbare);
  
  imsrgsolver.SetHin(Hbare);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDenominatorDelta(denominator_delta);

  if (method == "magnus")
  {
     imsrgsolver.SetdOmega(domega);
     imsrgsolver.SetOmegaNormMax(omega_norm_max);
     if (nsteps > 1)
     {
       imsrgsolver.SetGenerator(core_generator);
       imsrgsolver.Solve();
       smax *= 2;
     }
     imsrgsolver.SetGenerator(valence_generator);
     imsrgsolver.SetSmax(smax);
     imsrgsolver.Solve();
     rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
     rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
     R2_p1 = imsrgsolver.Transform(R2_p1);
     R2_p2 = imsrgsolver.Transform(R2_p2);
     R2_cm = imsrgsolver.Transform(R2_cm);
     rw.WriteNuShellX_op(R2_p1,intfile+"_R2p1.int");
     rw.WriteNuShellX_op(R2_p2,intfile+"_R2p2.int");
     rw.WriteNuShellX_op(R2_cm,intfile+"_R2cm.int");
  }
  else if (method == "flow")
  {
     imsrgsolver.SetODETolerance(ode_tolerance);
     if (nsteps > 1)
     {
       imsrgsolver.SetGenerator(core_generator);
       imsrgsolver.Solve_ode_adaptive();
     }
     imsrgsolver.SetGenerator(valence_generator);
     imsrgsolver.SetSmax(smax);
     imsrgsolver.Solve_ode_adaptive();
     rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
     rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
  }


  Hbare.PrintTimes();
 
  return 0;
}

