#include "ModelSpace.hh"
#include "ReadWrite.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include "IMSRGSolver.hh"
#include "imsrg_util.hh"
#include "AngMom.hh"
#include <string>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


using namespace boost::python;


BOOST_PYTHON_MODULE(pyIMSRG)
{

   class_<vector<string> > ("vector_string")
      .def (vector_indexing_suite< vector<string> >())
   ;

   class_<ModelSpace>("ModelSpace",init<>())
      .def(init<const ModelSpace&>())
      .def(init< int, const std::string&>())
      .def(init< int,vector<string>,vector<string> >())
      .def("Init", &ModelSpace::Init)
      .def("SetHbarOmega", &ModelSpace::SetHbarOmega)
      .def("SetTargetMass", &ModelSpace::SetTargetMass)
      .def("SetN3max", &ModelSpace::SetN3max)
      .def("GetHbarOmega", &ModelSpace::GetHbarOmega)
      .def("GetTargetMass", &ModelSpace::GetTargetMass)
      .def("GetNumberOrbits", &ModelSpace::GetNumberOrbits)
//      .def("PreComputeSixJs", &ModelSpace::PreComputeSixJs)
   ;

   class_<ReadWrite>("ReadWrite",init<>())
      .def("ReadSettingsFile", &ReadWrite::ReadSettingsFile)
      .def("ReadModelSpace", &ReadWrite::ReadModelSpace)
      .def("ReadBareTBME", &ReadWrite::ReadBareTBME)
      .def("ReadBareTBME_Jason", &ReadWrite::ReadBareTBME_Jason)
      .def("ReadBareTBME_Navratil", &ReadWrite::ReadBareTBME_Navratil)
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
      .def("GetOneBody", &Operator::GetOneBody)
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
      .def("GetParticleRank", &Operator::GetParticleRank)
      .def("GetE3max", &Operator::GetE3max)
      .def("SetE3max", &Operator::SetE3max)
      .def("PrintTimes", &Operator::PrintTimes)
      .def("BCH_Transform", &Operator::BCH_Transform)
   ;



   class_<HartreeFock>("HartreeFock",init<Operator&>())
      .def("Solve",&HartreeFock::Solve)
      .def("TransformToHFBasis",&HartreeFock::TransformToHFBasis)
      .def("GetHbare",&HartreeFock::GetHbare)
      .def("GetNormalOrderedH",&HartreeFock::GetNormalOrderedH)
      .def_readonly("EHF",&HartreeFock::EHF)
   ;


   class_<IMSRGSolver>("IMSRGSolver",init<Operator&>())
      .def("Solve",&IMSRGSolver::Solve)
#ifndef NO_ODE
      .def("Solve_ode",&IMSRGSolver::Solve_ode)
      .def("Solve_ode_adaptive",&IMSRGSolver::Solve_ode_adaptive)
      .def("Solve_ode_magnus",&IMSRGSolver::Solve_ode_magnus)
#endif
      .def("Transform",&IMSRGSolver::Transform)
      .def("InverseTransform",&IMSRGSolver::InverseTransform)
      .def("SetFlowFile",&IMSRGSolver::SetFlowFile)
      .def("SetDs",&IMSRGSolver::SetDs)
      .def("SetdOmega",&IMSRGSolver::SetdOmega)
      .def("SetSmax",&IMSRGSolver::SetSmax)
      .def("SetDsmax",&IMSRGSolver::SetDsmax)
      .def("SetHin",&IMSRGSolver::SetHin)
      .def("Reset",&IMSRGSolver::Reset)
      .def("SetGenerator",&IMSRGSolver::SetGenerator)
      .def("GetSystemDimension",&IMSRGSolver::GetSystemDimension)
      .def_readwrite("H_s", &IMSRGSolver::H_s)
      .def_readwrite("Omega", &IMSRGSolver::Omega)
   ;

   def("TCM_Op",           imsrg_util::TCM_Op);
   def("VCM_Op",           imsrg_util::VCM_Op);
   def("R2CM_Op",          imsrg_util::R2CM_Op);
   def("HCM_Op",           imsrg_util::HCM_Op);
   def("NumberOp",         imsrg_util::NumberOp);
   def("RSquaredOp",       imsrg_util::RSquaredOp);
   def("E0Op",             imsrg_util::E0Op);
   def("Isospin2_Op",      imsrg_util::Isospin2_Op);
   def("HO_density",       imsrg_util::HO_density);
   def("GetOccupationsHF", imsrg_util::GetOccupationsHF);
   def("GetOccupations",   imsrg_util::GetOccupations);
   def("GetDensity",       imsrg_util::GetDensity);


}
