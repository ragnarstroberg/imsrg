#include "ModelSpace.hh"
#include "ReadWrite.hh"
#include "Operator.hh"
//#include "Operator3.hh"
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
//      .def("ReadBareTBME_Navratil", &ReadWrite::ReadBareTBME_Navratil<Operator&>)
//      .def("ReadBareTBME_Navratil", &ReadWrite::ReadBareTBME_Navratil<Operator3&>)
      .def("ReadBareTBME_Darmstadt", &ReadWrite::ReadBareTBME_Darmstadt)
      .def("Read_Darmstadt_3body", &ReadWrite::Read_Darmstadt_3body)
      .def("WriteOneBody", &ReadWrite::WriteOneBody)
      .def("WriteTwoBody", &ReadWrite::WriteTwoBody)
      .def("WriteValenceOneBody", &ReadWrite::WriteValenceOneBody)
      .def("WriteValenceTwoBody", &ReadWrite::WriteValenceTwoBody)
      .def("WriteNuShellX_sps", &ReadWrite::WriteNuShellX_sps)
      .def("WriteNuShellX_int", &ReadWrite::WriteNuShellX_int)
      .def("WriteAntoine_int", &ReadWrite::WriteAntoine_int)
      .def("WriteOperator", &ReadWrite::WriteOperator)
      .def("ReadOperator", &ReadWrite::ReadOperator)
      .def("SetCoMCorr", &ReadWrite::SetCoMCorr)
      .def_readwrite("InputParameters", &ReadWrite::InputParameters)
   ;

   class_<Operator>("Operator",init<>())
      .def(init< ModelSpace&>())
      .def(init< ModelSpace&,int,int,int,int>())
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

/*
   class_<Operator3>("Operator3",init<>())
      .def(init< ModelSpace&>())
      .def(init< Operator3&>())
      .def(init< Operator&>())
      .def(self += Operator())
      .def(self + Operator())
      .def(self -= Operator())
      .def(self - Operator())
      .def(self *= double())
      .def(self * double())
      .def(self /= double())
      .def(self / double())
      .def_readwrite("ZeroBody", &Operator3::ZeroBody)
      .def_readwrite("OneBody", &Operator3::OneBody)
      .def("ScaleOneBody", &Operator3::ScaleOneBody)
      .def("ScaleTwoBody", &Operator3::ScaleTwoBody)
      .def("DoNormalOrdering", &Operator3::DoNormalOrdering)
      .def("DoNormalOrdering3", &Operator3::DoNormalOrdering3)
      .def("CalculateKineticEnergy", &Operator3::CalculateKineticEnergy)
      .def("Norm", &Operator3::Norm)
      .def("OneBodyNorm", &Operator3::OneBodyNorm)
      .def("TwoBodyNorm", &Operator3::TwoBodyNorm)
      .def("SetHermitian", &Operator3::SetHermitian)
      .def("SetAntiHermitian", &Operator3::SetAntiHermitian)
      .def("SetNonHermitian", &Operator3::SetNonHermitian)
   ;
*/


   class_<HartreeFock>("HartreeFock",init<Operator&>())
      .def("Solve",&HartreeFock::Solve)
      .def("TransformToHFBasis",&HartreeFock::TransformToHFBasis)
      .def("GetHbare",&HartreeFock::GetHbare)
      .def_readonly("EHF",&HartreeFock::EHF)
   ;

/*
   class_<HartreeFock<Operator>>("HartreeFock",init<Operator&>())
      .def("Solve",&HartreeFock<Operator>::Solve)
      .def("TransformToHFBasis",&HartreeFock<Operator>::TransformToHFBasis)
      .def("GetHbare",&HartreeFock<Operator>::GetHbare)
      .def_readonly("EHF",&HartreeFock<Operator>::EHF)
   ;


   class_<HartreeFock<Operator3>>("HartreeFock3",init<Operator3&>())
      .def("Solve",&HartreeFock<Operator3>::Solve)
      .def("TransformToHFBasis",&HartreeFock<Operator3>::TransformToHFBasis)
      .def("GetHbare",&HartreeFock<Operator3>::GetHbare)
      .def_readonly("EHF",&HartreeFock<Operator3>::EHF)
   ;
*/

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
      .def("GetSystemDimension",&IMSRGSolver::GetSystemDimension)
      .def_readwrite("H_s", &IMSRGSolver::H_s)
      .def_readwrite("Omega", &IMSRGSolver::Omega)
   ;

//   def("TCM_Op",           imsrg_util::TCM_Op<Operator>);
//   def("TCM_Op3",          imsrg_util::TCM_Op<Operator3>);
   def("VCM_Op",           imsrg_util::VCM_Op);
   def("NumberOp",         imsrg_util::NumberOp);
//   def("NumberOp",         imsrg_util::NumberOp<Operator>);
//   def("NumberOp",         imsrg_util::NumberOp<Operator3>);
   def("TCM_Op",           imsrg_util::TCM_Op);
   def("HO_density",       imsrg_util::HO_density);
   def("GetOccupationsHF", imsrg_util::GetOccupationsHF);
//   def("GetOccupationsHF", imsrg_util::GetOccupationsHF<Operator>);
//   def("GetOccupationsHF3", imsrg_util::GetOccupationsHF<Operator3>);
   def("GetOccupations",   imsrg_util::GetOccupations);
//   def("GetOccupations",   imsrg_util::GetOccupations<Operator>);
//   def("GetOccupations3",   imsrg_util::GetOccupations<Operator3>);
   def("GetDensity",       imsrg_util::GetDensity);


}
