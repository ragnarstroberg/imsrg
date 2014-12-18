#include "ModelSpace.hh"
#include "ReadWrite.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include "IMSRGSolver.hh"
#include "imsrg_util.hh"
#include "AngMom.hh"

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python.hpp>


using namespace boost::python;

BOOST_PYTHON_MODULE(pyIMSRG)
{

   class_<ModelSpace>("ModelSpace",init<>())
      .def("SetHbarOmega", &ModelSpace::SetHbarOmega)
      .def("SetTargetMass", &ModelSpace::SetTargetMass)
      .def("GetHbarOmega", &ModelSpace::GetHbarOmega)
      .def("GetTargetMass", &ModelSpace::GetTargetMass)
   ;

   class_<ReadWrite>("ReadWrite",init<>())
      .def("ReadSettingsFile", &ReadWrite::ReadSettingsFile)
      .def("ReadModelSpace", &ReadWrite::ReadModelSpace)
      .def("ReadBareTBME", &ReadWrite::ReadBareTBME)
      .def("ReadBareTBME_Jason", &ReadWrite::ReadBareTBME_Jason)
      .def("ReadBareTBME_Navratil", &ReadWrite::ReadBareTBME_Navratil)
      .def("ReadBareTBME_Darmstadt", &ReadWrite::ReadBareTBME_Darmstadt)
      .def("WriteOneBody", &ReadWrite::WriteOneBody)
      .def("WriteTwoBody", &ReadWrite::WriteTwoBody)
      .def("WriteValenceOneBody", &ReadWrite::WriteValenceOneBody)
      .def("WriteValenceTwoBody", &ReadWrite::WriteValenceTwoBody)
      .def("WriteNuShellX_sps", &ReadWrite::WriteNuShellX_sps)
      .def("WriteNuShellX_int", &ReadWrite::WriteNuShellX_int)
      .def("WriteAntoine_int", &ReadWrite::WriteAntoine_int)
      .def("SetCoMCorr", &ReadWrite::SetCoMCorr)
      .def_readwrite("InputParameters", &ReadWrite::InputParameters)
   ;

   class_<Operator>("Operator",init<>())
      .def(init< ModelSpace&>())
      .def(self += Operator())
      .def(self + Operator())
      .def(self -= Operator())
      .def(self - Operator())
      .def(self *= double())
      .def(self * double())
      .def(self /= double())
      .def(self / double())
      .def_readwrite("ZeroBody", &Operator::ZeroBody)
      .def_readwrite("OneBody", &Operator::OneBody)
      .def("ScaleOneBody", &Operator::ScaleOneBody)
      .def("ScaleTwoBody", &Operator::ScaleTwoBody)
      .def("DoNormalOrdering", &Operator::DoNormalOrdering)
      .def("CalculateKineticEnergy", &Operator::CalculateKineticEnergy)
      .def("Norm", &Operator::Norm)
      .def("OneBodyNorm", &Operator::OneBodyNorm)
      .def("TwoBodyNorm", &Operator::TwoBodyNorm)
      .def("SetHermitian", &Operator::SetHermitian)
      .def("SetAntiHermitian", &Operator::SetAntiHermitian)
      .def("SetNonHermitian", &Operator::SetNonHermitian)
      .def("Set_BCH_Transform_Threshold", &Operator::Set_BCH_Transform_Threshold)
      .def("Set_BCH_Product_Threshold", &Operator::Set_BCH_Product_Threshold)
      .def("PrintOneBody", &Operator::PrintOneBody)
      .def("PrintTwoBody", &Operator::PrintTwoBody)
   ;



   class_<HartreeFock>("HartreeFock",init<Operator&>())
      .def("Solve",&HartreeFock::Solve)
      .def("TransformToHFBasis",&HartreeFock::TransformToHFBasis)
      .def("GetHbare",&HartreeFock::GetHbare)
      .def_readonly("EHF",&HartreeFock::EHF)
   ;



   class_<IMSRGSolver>("IMSRGSolver",init<Operator>())
      .def("Solve",&IMSRGSolver::Solve)
      .def("Solve_ode",&IMSRGSolver::Solve_ode)
      .def("Solve_ode_magnus",&IMSRGSolver::Solve_ode_magnus)
      .def("Transform",&IMSRGSolver::Transform)
      .def("SetFlowFile",&IMSRGSolver::SetFlowFile)
      .def("SetDs",&IMSRGSolver::SetDs)
      .def("SetdOmega",&IMSRGSolver::SetdOmega)
      .def("SetSmax",&IMSRGSolver::SetSmax)
      .def("SetDsmax",&IMSRGSolver::SetDsmax)
      .def("SetGenerator",&IMSRGSolver::SetGenerator)
      .def_readwrite("H_s", &IMSRGSolver::H_s)
      .def_readwrite("Omega", &IMSRGSolver::Omega)

   ;

   def("TCM_Op",         imsrg_util::TCM_Op);
   def("VCM_Op",         imsrg_util::VCM_Op);
   def("NumberOp",       imsrg_util::NumberOp);
   def("TCM_Op",         imsrg_util::TCM_Op);
   def("HO_density",     imsrg_util::HO_density);
   def("GetOccupationsHF", imsrg_util::GetOccupationsHF);
   def("GetOccupations", imsrg_util::GetOccupations);
   def("GetDensity",     imsrg_util::GetDensity);


}
