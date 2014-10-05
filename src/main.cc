#include "ModelSpace.hh"
#include "ReadWrite.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include <iostream>

using namespace std;

int main()
{


//   char inputsps[100] = "../input/s.sps";
   char inputsps[100] = "../input/sd.sps";
   //char inputsps[100] = "../input/sd4.sps";

//   char inputtbme[100] = "../input/vspsd.int";
//   char inputtbme[100] = "../input/vs.int";
//   char inputtbme[100] = "../input/vjason2.0.int";
//   char inputtbme[100] = "../input/vjason1.8.int";
//   char inputtbme[100] = "../input/tbme1.int";
   //char inputtbme[100] = "../input/tbme2.int";
   char inputtbme[100] = "../input/vsrg1.8_n3lo500_w_coulomb_emax6_hw24.int";

   //ModelSpace * modelspace = new ModelSpace();
   ReadWrite rw = ReadWrite();
   ModelSpace modelspace = rw.ReadModelSpace(inputsps);
   modelspace.SetupKets();
   ModelSpace ms2 = modelspace;
   cout << "Read in the modelspace. Norbits = " << modelspace.GetNumberOrbits() << endl;
   //rw.ReadModelSpace(inputsps,modelspace);
   cout << "Initializing the bare operator..." << endl;
   Operator H_bare =  Operator(&modelspace);
//   Operator H_hf   =  Operator(&modelspace);
   cout << "Calculating the kinetic energy..." << endl;
   rw.CalculateKineticEnergy(&H_bare);
   cout << "Reading in the TBME's..." << endl;
   rw.ReadBareTBME(inputtbme, H_bare);
   cout << "Done reading in the interaction" << endl;

   HartreeFock  hf = HartreeFock(&H_bare);
   cout << "Solving HF equations..." << endl;
   hf.Solve();
   cout << "Transforming to HF basis..." << endl;
//   hf.PrintOrbits();
//   hf.TransformToHFBasis(&H_bare,&H_hf);
   Operator H_hf = hf.TransformToHFBasis(H_bare);
/*
   H_bare.WriteOneBody("../output/Hbare1b.out");
   H_bare.WriteTwoBody("../output/Hbare2b.out");
   H_hf.WriteOneBody("../output/HF1b.out");
   H_hf.WriteTwoBody("../output/HF2b.out");
*/
   cout << "Writing to file..." << endl;
   rw.WriteOneBody(H_bare,"../output/Hbare1b.out");
   rw.WriteTwoBody(H_bare,"../output/Hbare2b.out");
   rw.WriteOneBody(H_hf,"../output/HF1b.out");
   rw.WriteTwoBody(H_hf,"../output/HF2b.out");
//   H_bare->OneBody.print();
//   cout << endl << endl;
//   H_hf->OneBody.print();


   cout << "Repeating the HF procedure on the HF-basis Hamiltonian" << endl;
   cout << " ==============================================================================" << endl;
   HartreeFock  hf2 = HartreeFock(&H_hf);
   hf2.Solve();

   cout << "EHF = " << hf.EHF << endl;
   cout << "EHF2 = " << hf2.EHF << endl;
   cout << "Start normal ordering" << endl;
   Operator HFNO = H_hf.DoNormalOrdering();
   cout << "Normal ordered zero-body part = " << HFNO.ZeroBody << endl;
   cout << "Normal ordered one-body part: " << endl;
//   HFNO.OneBody.print();

  return 0;
}
