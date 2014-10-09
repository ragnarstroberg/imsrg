#include "ModelSpace.hh"
#include "ReadWrite.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include <iostream>

using namespace std;

int main()
{


   char inputsps[100] = "../input/s.sps";
//   char inputsps[100] = "../input/sd.sps";
   //char inputsps[100] = "../input/sd4.sps";

   char inputtbme[100] = "../input/vspsd.int";
//   char inputtbme[100] = "../input/vs.int";
//   char inputtbme[100] = "../input/vjason2.0.int";
//   char inputtbme[100] = "../input/vjason1.8.int";
//   char inputtbme[100] = "../input/tbme1.int";
//   char inputtbme[100] = "../input/tbme2.int";
//   char inputtbme[100] = "../input/vsrg1.8_n3lo500_w_coulomb_emax6_hw24.int";

   //ModelSpace * modelspace = new ModelSpace();
   ReadWrite rw = ReadWrite();
   cout << "Reading in the modelspace" << endl;
   ModelSpace modelspace = rw.ReadModelSpace(inputsps);
   cout << "Setting up the kets" << endl;
   modelspace.SetupKets();
   cout << "Done setting up the kets" << endl;
   cout << "Creating H_bare" << endl;
   Operator H_bare =  Operator(&modelspace);
   cout << "Calculating the kinetic energy" << endl;
   rw.CalculateKineticEnergy(&H_bare);
   cout << "Reading in the TBME " << endl;
   rw.ReadBareTBME(inputtbme, H_bare);

   cout << "setting up HF" << endl;
   HartreeFock  hf = HartreeFock(&H_bare);
   cout << "solving HF" << endl;
   hf.Solve();

   Operator H_hf = hf.TransformToHFBasis(H_bare);
   cout << "After transformation, one body piece is"<< endl;
   H_hf.OneBody.print();
/*
   H_bare.WriteOneBody("../output/Hbare1b.out");
   H_bare.WriteTwoBody("../output/Hbare2b.out");
   H_hf.WriteOneBody("../output/HF1b.out");
   H_hf.WriteTwoBody("../output/HF2b.out");
*/
/*   cout << "Writing to file..." << endl;
   rw.WriteOneBody(H_bare,"../output/Hbare1b.out");
   rw.WriteTwoBody(H_bare,"../output/Hbare2b.out");
   rw.WriteOneBody(H_hf,"../output/HF1b.out");
   rw.WriteTwoBody(H_hf,"../output/HF2b.out");
*/
//   H_bare->OneBody.print();
//   cout << endl << endl;
//   H_hf->OneBody.print();


   cout << "Repeating the HF procedure on the HF-basis Hamiltonian" << endl;
   cout << " ==============================================================================" << endl;
   HartreeFock  hf2 = HartreeFock(&H_hf);
   hf2.Solve();
   cout << "Done solving" << endl;

   cout << "EHF = " << hf.EHF << endl;
   cout << "EHF2 = " << hf2.EHF << endl;
   cout << "Start normal ordering" << endl;
   Operator HFNO = H_hf.DoNormalOrdering();
   Operator HbareNO = H_bare.DoNormalOrdering();
   cout << "Normal ordered zero-body part = " << HFNO.ZeroBody << endl;
   cout << "Normal ordered one-body part: " << endl;
   HFNO.OneBody.print();
   HbareNO.OneBody.print();
   Operator Hcomm = HFNO.Commutator(HbareNO);
//   Operator Hcomm = HbareNO.Commutator(H_bare);
   cout << "Commutator zerobody: " << Hcomm.ZeroBody << endl;
   cout << "Commutator one body:" << endl;;
   Hcomm.OneBody.print();

  return 0;
}
