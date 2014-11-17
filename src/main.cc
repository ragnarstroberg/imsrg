#include "ModelSpace.hh"
#include "ReadWrite.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include "IMSRGSolver.hh"
#include "imsrg_util.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

int main(int argc, char**argv)
{

   ReadWrite rw = ReadWrite();
   // If a settings file isn't given as a command line arg, look for the default file
   string settings_file = argc>1 ? string(argv[1]) : "settings.inp";
   rw.ReadSettingsFile(settings_file);

   // These are parameters that may be used in the settings file
   string inputsps	= rw.InputParameters["inputsps"];
   string inputtbme	= rw.InputParameters["inputtbme"];
   string darmstadttbme	= rw.InputParameters["darmstadttbme"];
   string darmstadtEmax	= rw.InputParameters["darmstadtEmax"];
   string flowfile	= rw.InputParameters["flowfile"];
   string ds_str	= rw.InputParameters["ds"];
   string smax_str	= rw.InputParameters["smax"];
   string generator	= rw.InputParameters["generator"];

   if (generator == "") generator = "white";


   cout << "Reading in the modelspace from " << inputsps << endl;
   ModelSpace modelspace = rw.ReadModelSpace(inputsps);

   Operator H_bare =  Operator(&modelspace);
   H_bare.SetHermitian(); // just to be sure

   H_bare.CalculateKineticEnergy();

   if (inputtbme != "")
   {
      cout << "Reading Oslo-style TBME from " << inputtbme << endl;
      rw.ReadBareTBME(inputtbme, H_bare);
      rw.WriteTwoBody(H_bare,"../output/Oslo_H_bare.out");
   }
   else if (darmstadttbme != "")
   {
      int Emax = darmstadtEmax != "" ? atoi(darmstadtEmax.c_str()) : -1;
      cout << "Reading Darmstadt-style TBME from " << darmstadttbme << " with Emax " << Emax << endl;
      rw.ReadBareTBME_Darmstadt(darmstadttbme, H_bare, Emax);
      rw.WriteTwoBody(H_bare,"../output/Darmstadt_H_bare.out");
   }


   cout << "Norm of H_bare = " << H_bare.Norm() << endl;

   HartreeFock  hf = HartreeFock(H_bare);
   hf.Solve();

   Operator H_hf = hf.TransformToHFBasis(H_bare);
   HartreeFock hf2 = HartreeFock(H_hf);
   hf2.Solve();



   cout << "EHF = " << hf.EHF << endl;
   cout << "EHF2 = " << hf2.EHF << endl;
   Operator HFNO = H_hf.DoNormalOrdering();
   Operator HbareNO = H_bare.DoNormalOrdering();

   cout << "Norm of HFNO = " << HFNO.Norm() << endl;

   IMSRGSolver imsrgsolver = IMSRGSolver(HFNO);
//   IMSRGSolver imsrgsolver = IMSRGSolver(HbareNO);
   imsrgsolver.SetFlowFile(flowfile);
   imsrgsolver.SetGenerator(generator);

   if (ds_str != "")
   {
      double ds = strtod(ds_str.c_str(),NULL);
      imsrgsolver.SetDs(ds);
      imsrgsolver.SetdOmega(ds);
   }
   if (smax_str != "")
   {
      double smax = strtod(smax_str.c_str(),NULL);
      imsrgsolver.SetSmax(smax);
   }

/////// THIS IS THE TIME CONSUMING PART //////////
   imsrgsolver.Solve();
//////////////////////////////////////////////////



 /* 
  cout << endl << " density: " << endl;

  int nr_steps = 100;
  vector<double> R(nr_steps,0);
  double dr = 0.1;
  for (int i=0;i<nr_steps;++i) R[i] = i*dr;

  cout << "Calculating occupation numbers..." << endl;

  //vector<double> occ = imsrg_util::GetOccupations(hf,imsrgsolver);
  vector<double> occ = imsrg_util::GetOccupations(hf);

  vector<double>dens_P = imsrg_util::GetDensity(occ,R,modelspace.proton_orbits,modelspace);
  vector<double>dens_N = imsrg_util::GetDensity(occ,R,modelspace.neutron_orbits,modelspace);

  for (int i=0;i<nr_steps;++i)
  {
     cout << R[i] << " " << dens_P[i] << " " << dens_N[i] << endl;
  }
*/

//  Operator Np0s1 = StandardOperators::NumberOp(modelspace,0,0,1,-1); // proton 0s1/2
//  Operator Np0s1_hf = hf.TransformToHFBasis(Np0s1);
//  Operator Np0s1_NO = Np0s1_hf.DoNormalOrdering();
//  Operator Np0s1_final = imsrgsolver.Transform(Np0s1_NO);
//  cout << "proton 0s1/2 occupation = " << Np0s1_final.ZeroBody << endl;

//   rw.WriteValenceOneBody(HFNO,"../output/O16_lmax6_SM_1b_bare.int");
//   rw.WriteValenceTwoBody(HFNO,"../output/O16_lmax6_SM_2b_bare.int");
//   rw.WriteValenceOneBody(imsrgsolver.H_s,"../output/O16_lmax6_SM_1b.int");
//   rw.WriteValenceTwoBody(imsrgsolver.H_s,"../output/O16_lmax6_SM_2b.int");

   rw.WriteNuShellX_int(imsrgsolver.H_s,"../output/Ca40srg.int");
   rw.WriteNuShellX_sps(imsrgsolver.H_s,"../output/Ca40srg.sp");

//   rw.WriteNuShellX_int(HFNO,"../output/He4_bare.int");

  return 0;
}
