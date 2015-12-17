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
  string LECs = PAR.s("LECs");

  int eMax = PAR.i("emax");
  int E3max = PAR.i("e3max");
  int lmax3 = PAR.i("lmax3");
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

  vector<string> opnames = PAR.v("Operators");

  vector<Operator> ops;



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
  rw.SetLECs_preset(LECs);
  ModelSpace modelspace = reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space);
//  ModelSpace ms2;

//  if (reference != "default" and reference != valence_space)
//       ms2 = ModelSpace(eMax,valence_space); // Targeted normal ordering. We'll need this for the last step.

  if (nsteps < 0)
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  modelspace.SetE3max(E3max);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);
  
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

  // Calculate all the desired operators
  for (auto& opname : opnames)
  {
           if (opname == "R2_p1")        ops.emplace_back( R2_1body_Op(modelspace,"proton") );
      else if (opname == "R2_p2")        ops.emplace_back( R2_2body_Op(modelspace,"proton") );
      else if (opname == "R2_n1")        ops.emplace_back( R2_1body_Op(modelspace,"neutron") );
      else if (opname == "R2_n2")        ops.emplace_back( R2_2body_Op(modelspace,"neutron") );
      else if (opname == "Rp2")          ops.emplace_back( Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "Rn2")          ops.emplace_back( Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "Rm2")          ops.emplace_back( Rm2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "E2")           ops.emplace_back( ElectricMultipoleOp(modelspace,2) );
      else if (opname == "M1")           ops.emplace_back( MagneticMultipoleOp(modelspace,1) );
      else if (opname == "Fermi")        ops.emplace_back( AllowedFermi_Op(modelspace) );
      else if (opname == "GamowTeller")  ops.emplace_back( AllowedGamowTeller_Op(modelspace) );
      else if (opname == "R2CM")         ops.emplace_back( R2CM_Op(modelspace) );
      else if (opname == "HCM")          ops.emplace_back( HCM_Op(modelspace) );
      else if (opname.substr(0,4) == "HCM_") // GetHCM with a different frequency, ie HCM_24 for hw=24
      {
         double hw_HCM;
         double hw_save = modelspace.GetHbarOmega();
         istringstream(opname.substr(4,opname.size())) >> hw_HCM;
         modelspace.SetHbarOmega(hw_HCM);
         ops.emplace_back( HCM_Op(modelspace) );
         modelspace.SetHbarOmega(hw_save);
      }
      else //need to remove from the list
      {
         cout << "Unknown operator: " << opname << endl;
      }
  }


  

  for (auto& op : ops)
  {
     if (basis == "HF") op = hf.TransformToHFBasis(op);
     op = op.DoNormalOrdering();
  }
  auto itR2p = find(opnames.begin(),opnames.end(),"Rp2");
  if (itR2p != opnames.end())
  {
    Operator& Rp2 = ops[itR2p-opnames.begin()];
    int Z = modelspace.GetTargetZ();
    int A = modelspace.GetTargetMass();
    cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << endl; 
    cout << " HF charge radius = " << sqrt( Rp2.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << endl; 
  }
  
  if ( method == "HF" )
  {
   Hbare.PrintTimes();
   return 0;
  }

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



  // Transform all the operators
  if (method == "magnus")
  {
    if (ops.size()>0) cout << "transforming operators" << endl;
    for (size_t i=0;i<ops.size();++i)
    {
      cout << opnames[i] << " " << flush;
      ops[i] = imsrgsolver.Transform(ops[i]);
    }
    cout << endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  // If we're doing targeted normal ordering 
  // we now re-normal order wrt to the core
  // and do any remaining flow.
//  if (reference != "default"  and reference != valence_space)
  if ( modelspace.core != modelspace.holes )
  {

    Hbare = imsrgsolver.GetH_s();

    int nOmega = imsrgsolver.GetOmegaSize();
    cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << endl;
    Hbare = Hbare.UndoNormalOrdering();

    ModelSpace ms2(modelspace); // copy the current model space
    ms2.SetReference(ms2.core); // chage the reference determinant
    Hbare.SetModelSpace(ms2);

    cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << endl;
    Hbare = Hbare.DoNormalOrdering();

    imsrgsolver.SetHin(Hbare);
    imsrgsolver.Solve();
    // Change to the new basis, then apply the rest of the transformation to the operators
    for (auto& op : ops)
    {
      op = op.UndoNormalOrdering();
      op.SetModelSpace(ms2);
      op = op.DoNormalOrdering();
      // transform using the remaining omegas
      op = imsrgsolver.Transform_Partial(op,nOmega);
    }
  }



  // Write the output

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
  if (modelspace.valence.size() > 0)
  {
    rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
    rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");

    if (method == "magnus")
    {
       for (int i=0;i<ops.size();++i)
       {
          ops[i] = imsrgsolver.Transform(ops[i]);
          rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
       }
    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    cout << "Core Energy = " << setprecision(6) << imsrgsolver.GetH_s().ZeroBody << endl;
    for (int i=0;i<ops.size();++i)
    {
      Operator& op = ops[i];
      cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
      if ( opnames[i] == "Rp2" )
      {
         int Z = modelspace.GetTargetZ();
         int A = modelspace.GetTargetMass();
         cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << endl; 
         cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << endl; 
      }
    }
  }




  Hbare.PrintTimes();
 
  return 0;
}

