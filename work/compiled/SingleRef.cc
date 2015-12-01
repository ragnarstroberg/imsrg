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
  string generator = "atan";
  int E3max = 12;
  int eMax = 6;
  string flowfile = "default";
  string intfile = "default";
  string fmt2 = "me2j";
  string nucleus = "O16";
  string basis = "HF";
  string method = "magnus";
  double ode_tolerance = 1e-6;
  double R2p = 0.770; // = 0.8775^2
  double R2n = -0.1149;
  double DF = 0.033;
  int file2e1max=12;
  int file2e2max=24;
  int file2lmax=10;
  int file3e1max=12;
  int file3e2max=24;
  int file3e3max=12;


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
    else if (var == "nucleus")
      istringstream(val) >> nucleus;
    else if (var == "smax")
      istringstream(val) >> smax;
    else if (var == "fmt2")
      istringstream(val) >> fmt2;
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
    else if (var == "generator")
      istringstream(val) >> generator;
    else
      cout << "unknown parameter : " << var << endl;
  }


  char strbuf[200];
  if (flowfile == "default")
  {
    sprintf(strbuf, "output/BCH_flow_SR_A%d_emax%d.dat",targetMass,eMax);
    flowfile = strbuf;
  }
  if ( intfile == "default")
  {
    sprintf(strbuf, "output/SDA%d_2stepe%d",targetMass,eMax);
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



  ModelSpace modelspace(eMax,nucleus);

  modelspace.SetHbarOmega(hw);
  if (targetMass > 0)
    modelspace.SetTargetMass(targetMass);
  modelspace.SetN3max(E3max);
  
  cout << "Making the operator..." << endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();

  cout << "Reading interaction..." << endl;
  if (fmt2 == "me2j")
    rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
  else if (fmt2 == "navratil" or fmt2 == "Navratil")
    rw.ReadBareTBME_Navratil(inputtbme, Hbare);

  if (Hbare.particle_rank >=3)
  {
//    rw.Read_Darmstadt_3body(input3bme, Hbare, 12,24,12);
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

  Operator Rp2 = Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ());
  if (basis == "HF")
    Rp2 = hf.TransformToHFBasis(Rp2);
  Rp2 = Rp2.DoNormalOrdering();
  cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << endl; 
  cout << " HF charge radius = " << sqrt( Rp2.ZeroBody + R2p + R2n + DF) << endl; 
  
  IMSRGSolver imsrgsolver(Hbare);

  imsrgsolver.SetMethod(method);
  imsrgsolver.SetGenerator(generator);
  imsrgsolver.SetHin(Hbare);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
//  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetODETolerance(ode_tolerance);
  imsrgsolver.SetdOmega(min(domega,omega_norm_max));
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.Solve();

  if (method == "magnus" or method == "flow-omega")
  {
   Rp2 = imsrgsolver.Transform(Rp2);
   cout << " point proton radius = " << sqrt( Rp2.ZeroBody ) << endl; 
   cout << "charge radius = " << sqrt( Rp2.ZeroBody + R2p + R2n + DF) << endl; 
  }

//  char buf[500];
//  sprintf(buf,"Omega_e%d_%s.op",eMax,nucleus.c_str());
//  string omegafile(buf);
//  rw.WriteOperator(imsrgsolver.GetOmega(0),omegafile);
  Hbare.PrintTimes();

  return 0;
}

