#include "ModelSpace.hh"
#include "ReadWrite.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include "IMSRGSolver.hh"
#include <iostream>

using namespace std;

int main(int argc, char**argv)
{

   ReadWrite rw = ReadWrite();

   char default_settings_file[200] = "settings.inp";

   if (argc>1)
      rw.ReadSettingsFile(argv[1]);
   else
      rw.ReadSettingsFile(default_settings_file);

   string inputsps = rw.InputParameters["inputsps"];
   string inputtbme = rw.InputParameters["inputtbme"];

   cout << "Reading in the modelspace" << endl;
   ModelSpace modelspace = rw.ReadModelSpace(inputsps.c_str());

   cout << "Setting up the kets" << endl;
   modelspace.SetupKets();

   cout << "Creating H_bare" << endl;
   Operator H_bare =  Operator(&modelspace);

   cout << "Calculating the kinetic energy" << endl;
   rw.CalculateKineticEnergy(&H_bare);

   cout << "Reading in the TBME " << endl;
   rw.ReadBareTBME(inputtbme.c_str(), H_bare);

   cout << "Norm of H_bare = " << H_bare.Norm() << endl;

   cout << "setting up HF" << endl;
   HartreeFock  hf = HartreeFock(&H_bare);

   cout << "solving HF" << endl;
   hf.Solve();

   Operator H_hf = hf.TransformToHFBasis(H_bare);

   cout << "After transformation, one body piece is"<< endl;
   H_hf.OneBody.print();

   cout << "***** Bare Two Body Part *****" << endl;
   H_bare.PrintTwoBody();

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
//   HbareNO.OneBody.print();
//   Operator Hcomm = HFNO.Commutator(HbareNO);
//   cout << "Commutator zerobody: " << Hcomm.ZeroBody << endl;
//   cout << "Commutator one body:" << endl;;
//   Hcomm.OneBody.print();

   //IMSRGSolver imsrgsolver = IMSRGSolver(HFNO);
   IMSRGSolver imsrgsolver = IMSRGSolver(HbareNO);
   imsrgsolver.Solve();

  return 0;
}
