#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "IMSRG.hh"
#include "Parameters.hh"

using namespace imsrg_util;

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
  Parameters PAR(argc,argv);

  string inputtbme = PAR.s("2bme");
  string input3bme = PAR.s("3bme");
  string reference = PAR.s("reference");
  string valence_space = PAR.s("valence_space");
  string basis = PAR.s("basis");
  string method = PAR.s("method");
  string flowfile = PAR.s("flowfile");
  string intfile = PAR.s("intfile");
  string core_generator = PAR.s("core_generator");
  string valence_generator = PAR.s("valence_generator");
  string fmt2 = PAR.s("fmt2");
  string denominator_delta_orbit = PAR.s("denominator_delta_orbit");

  int eMax = PAR.i("emax");
  int E3max = PAR.i("e3max");
  int targetMass = PAR.i("A");
  int nsteps = PAR.i("nsteps");
  int file2e1max = PAR.i("file2e1max");
  int file2e2max = PAR.i("file2e2max");
  int file2lmax = PAR.i("file2lmax");
  int file3e1max = PAR.i("file3e1max");
  int file3e2max = PAR.i("file3e2max");
  int file3e3max = PAR.i("file3e3max");

  double hw = PAR.d("hw");
  double smax = PAR.d("smax");
  double ode_tolerance = PAR.d("ode_tolerance");
  double ds_max = PAR.d("ds_max");
  double ds_0 = PAR.d("ds_0");
  double domega = PAR.d("domega");
  double omega_norm_max = PAR.d("omega_norm_max"); 
  double denominator_delta = PAR.d("denominator_delta");

  // test 2bme file
  ifstream test(inputtbme);
  if( not test.good() )
  {
    cout << "trouble reading " << inputtbme << " exiting. " << endl;
    return 1;
  }
  test.close();
  // test 3bme file
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
  ModelSpace modelspace = reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space);
  ModelSpace ms2;

  if (reference != "default" and reference != "valece_space")
       ms2 = ModelSpace(eMax,valence_space); // Targeted normal ordering. We'll need this for the last step.

  if (modelspace.valence.size() > 0) // we've got a valence space decoupling
  {
    if (nsteps < 0 ) nsteps = 2;
  }
  else // just doing a single reference decoupling
  {
    if (nsteps < 0 ) nsteps = 1;
  }



  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
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
  else if (fmt2 == "oslo" )
    rw.ReadTBME_Oslo(inputtbme, Hbare);

  if (Hbare.particle_rank >=3)
  {
    rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
  }  

//  cout << "Aeff = " << modelspace.GetTargetMass() << endl;
  Hbare += Trel_Op(modelspace);

//  cout << "Just before HF, hole orbits:  ";
//  for (auto& h : Hbare.GetModelSpace()->holes) cout << h << " ";
//  cout << endl;
  HartreeFock hf(Hbare);
  hf.Solve();
  cout << "EHF = " << hf.EHF << endl;
  
  if (basis == "HF")
    Hbare = hf.GetNormalOrderedH();
  else if (basis == "oscillator")
    Hbare = Hbare.DoNormalOrdering();

//  Operator Rp2   = Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ());
  Operator R2_p1 = R2_p1_Op(modelspace);
  Operator R2_p2 = R2_p2_Op(modelspace);
  Operator R2_cm  = R2CM_Op(modelspace);
  if (basis == "HF")
  {
//    Rp2 = hf.TransformToHFBasis(Rp2);
    R2_p1 = hf.TransformToHFBasis(R2_p1);
    R2_p2 = hf.TransformToHFBasis(R2_p2);
    R2_cm = hf.TransformToHFBasis(R2_cm);
  }
//  Rp2 = Rp2.DoNormalOrdering();
  R2_p1 = R2_p1.DoNormalOrdering();
  R2_p2 = R2_p2.DoNormalOrdering();
  R2_cm = R2_cm.DoNormalOrdering();
//  cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << endl; 
//  cout << " HF charge radius = " << sqrt( Rp2.ZeroBody + R2p + R2n + DF) << endl; 
  
  IMSRGSolver imsrgsolver(Hbare);
  
  if (method == "NSmagnus") // "No split" magnus
  {
    omega_norm_max=500;
    method = "magnus";
  }

  imsrgsolver.SetMethod(method);
  imsrgsolver.SetHin(Hbare);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);
  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

  if (nsteps > 1) // two-step decoupling, do core first
  {
    imsrgsolver.SetGenerator(core_generator);
    imsrgsolver.Solve();
    if (method == "magnus") smax *= 2;
  }

  imsrgsolver.SetGenerator(valence_generator);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.Solve();

  // If we're doing targeted normal ordering for the shell model
  // we now re-normal order wrt to the shell model core
  // and do any remaining flow.
  if (reference != "default"  and reference != valence_space)
  {
    if (method == "magnus")
    {
      smax *= 1.5;
      R2_p1 = imsrgsolver.Transform(R2_p1);
      R2_p2 = imsrgsolver.Transform(R2_p2);
      R2_cm = imsrgsolver.Transform(R2_cm);
      imsrgsolver.SetSmax(smax);
    }
    Hbare = imsrgsolver.GetH_s();
    cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << endl;
    Hbare = Hbare.UndoNormalOrdering();
    Hbare.SetModelSpace(ms2);
    cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << endl;
    Hbare = Hbare.DoNormalOrdering();
    imsrgsolver.SetHin(Hbare);
    imsrgsolver.Solve();
  }
  else if (method=="magnus" and modelspace.valence.size()==0) // run-of-the-mill single-ref calculation
  {
//     Rp2 = imsrgsolver.Transform(Rp2);
     R2_p1 = imsrgsolver.Transform(R2_p1);
     R2_p2 = imsrgsolver.Transform(R2_p2);
     R2_cm = imsrgsolver.Transform(R2_cm);
     int Z = modelspace.GetTargetZ();
     int A = modelspace.GetTargetMass();
     double rpp = R2_cm.ZeroBody + (A-2.0)/(A*Z)*R2_p1.ZeroBody -2.0/(A*Z)*R2_p2.ZeroBody;
     double rch = rpp + R2p + R2n + DF;
     cout << "point proton radius: " << sqrt(rpp) << endl;
     cout << "charge radius: " << sqrt(rch) << endl;
//     cout << "point proton radius: " << sqrt(Rp2.ZeroBody) << endl;
  }

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
  if (modelspace.valence.size() > 0)
  {
    rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
    rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");

    if (method == "magnus")
    {
       R2_p1 = imsrgsolver.Transform(R2_p1);
       R2_p2 = imsrgsolver.Transform(R2_p2);
       R2_cm = imsrgsolver.Transform(R2_cm);
       rw.WriteNuShellX_op(R2_p1,intfile+"_R2p1.int");
       rw.WriteNuShellX_op(R2_p2,intfile+"_R2p2.int");
       rw.WriteNuShellX_op(R2_cm,intfile+"_R2cm.int");
    }
  }
  cout << "E0 = " << setprecision(6) << imsrgsolver.GetH_s().ZeroBody << endl;

  Hbare.PrintTimes();
 
  return 0;
}

