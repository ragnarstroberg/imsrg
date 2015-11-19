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
      .def(init< int, int ,int>())
//      .def("Init", &ModelSpace::Init)
      .def("SetHbarOmega", &ModelSpace::SetHbarOmega)
      .def("SetTargetMass", &ModelSpace::SetTargetMass)
      .def("SetN3max", &ModelSpace::SetN3max)
      .def("GetHbarOmega", &ModelSpace::GetHbarOmega)
      .def("GetTargetMass", &ModelSpace::GetTargetMass)
      .def("GetNumberOrbits", &ModelSpace::GetNumberOrbits)
      .def("GetNumberKets", &ModelSpace::GetNumberKets)
//      .def("PreComputeSixJs", &ModelSpace::PreComputeSixJs)
   ;


   class_<Operator>("Operator",init<>())
      .def(init< ModelSpace&>())
      .def(init< ModelSpace&,int,int,int,int>())
      .def(self += Operator())
      .def(self + Operator())
      .def(self -= Operator())
      .def(self - Operator())
      .def( - self)
      .def(self *= double())
      .def(self * double())
      .def(self /= double())
      .def(self / double())
      .def_readwrite("ZeroBody", &Operator::ZeroBody)
      .def_readwrite("OneBody", &Operator::OneBody)
      .def("GetOneBody", &Operator::GetOneBody)
      .def("GetTwoBody", &Operator::GetTwoBody)
      .def("GetTwoBodyDimension", &Operator::GetTwoBodyDimension)
      .def("ScaleOneBody", &Operator::ScaleOneBody)
      .def("ScaleTwoBody", &Operator::ScaleTwoBody)
      .def("DoNormalOrdering", &Operator::DoNormalOrdering)
      .def("UndoNormalOrdering", &Operator::UndoNormalOrdering)
      .def("SetModelSpace", &Operator::SetModelSpace)
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
//      .def("EraseThreeBody", &Operator::EraseThreeBody)
      .def("Size", &Operator::Size)
   ;



   class_<ReadWrite>("ReadWrite",init<>())
      .def("ReadSettingsFile", &ReadWrite::ReadSettingsFile)
      .def("ReadTBME_Oslo", &ReadWrite::ReadTBME_Oslo)
      .def("ReadBareTBME_Jason", &ReadWrite::ReadBareTBME_Jason)
      .def("ReadBareTBME_Navratil", &ReadWrite::ReadBareTBME_Navratil)
      .def("ReadBareTBME_Darmstadt", &ReadWrite::ReadBareTBME_Darmstadt)
      .def("Read_Darmstadt_3body", &ReadWrite::Read_Darmstadt_3body)
      .def("Read3bodyHDF5", &ReadWrite::Read3bodyHDF5)
      .def("Write_me2j", &ReadWrite::Write_me2j)
      .def("Write_me3j", &ReadWrite::Write_me3j)
      .def("WriteTBME_Navratil", &ReadWrite::WriteTBME_Navratil)
      .def("WriteNuShellX_sps", &ReadWrite::WriteNuShellX_sps)
      .def("WriteNuShellX_int", &ReadWrite::WriteNuShellX_int)
      .def("WriteNuShellX_op", &ReadWrite::WriteNuShellX_op)
      .def("WriteAntoine_int", &ReadWrite::WriteAntoine_int)
      .def("WriteOperator", &ReadWrite::WriteOperator)
      .def("WriteOperatorHuman", &ReadWrite::WriteOperatorHuman)
      .def("ReadOperator", &ReadWrite::ReadOperator)
      .def("CompareOperators", &ReadWrite::CompareOperators)
      .def("ReadOneBody_Takayuki", &ReadWrite::ReadOneBody_Takayuki)
      .def("ReadTwoBody_Takayuki", &ReadWrite::ReadTwoBody_Takayuki)
      .def("WriteOneBody_Takayuki", &ReadWrite::WriteOneBody_Takayuki)
      .def("WriteTwoBody_Takayuki", &ReadWrite::WriteTwoBody_Takayuki)
      .def("WriteTensorOneBody", &ReadWrite::WriteTensorOneBody)
      .def("WriteTensorTwoBody", &ReadWrite::WriteTensorTwoBody)
      .def("WriteOneBody_Oslo", &ReadWrite::WriteOneBody_Oslo)
      .def("WriteTwoBody_Oslo", &ReadWrite::WriteTwoBody_Oslo)
      .def("SetCoMCorr", &ReadWrite::SetCoMCorr)
      .def("ReadTwoBodyEngel", &ReadWrite::ReadTwoBodyEngel)
      .def_readwrite("InputParameters", &ReadWrite::InputParameters)
   ;





   class_<HartreeFock>("HartreeFock",init<Operator&>())
      .def("Solve",&HartreeFock::Solve)
      .def("TransformToHFBasis",&HartreeFock::TransformToHFBasis)
      .def("GetHbare",&HartreeFock::GetHbare)
      .def("GetNormalOrderedH",&HartreeFock::GetNormalOrderedH)
      .def("GetOmega",&HartreeFock::GetOmega)
      .def("PrintSPE",&HartreeFock::PrintSPE)
      .def_readonly("EHF",&HartreeFock::EHF)
   ;

   // Define which overloaded version of IMSRGSolver::Transform I want to expose
   Operator (IMSRGSolver::*Transform_ref)(Operator&) = &IMSRGSolver::Transform;

   class_<IMSRGSolver>("IMSRGSolver",init<Operator&>())
      .def("Solve",&IMSRGSolver::Solve)
      .def("Transform",Transform_ref)
      .def("InverseTransform",&IMSRGSolver::InverseTransform)
      .def("SetFlowFile",&IMSRGSolver::SetFlowFile)
      .def("SetMethod",&IMSRGSolver::SetMethod)
      .def("SetEtaCriterion",&IMSRGSolver::SetEtaCriterion)
      .def("SetDs",&IMSRGSolver::SetDs)
      .def("SetdOmega",&IMSRGSolver::SetdOmega)
      .def("SetOmegaNormMax",&IMSRGSolver::SetOmegaNormMax)
      .def("SetSmax",&IMSRGSolver::SetSmax)
      .def("SetDsmax",&IMSRGSolver::SetDsmax)
      .def("SetHin",&IMSRGSolver::SetHin)
      .def("SetODETolerance",&IMSRGSolver::SetODETolerance)
      .def("Reset",&IMSRGSolver::Reset)
      .def("SetGenerator",&IMSRGSolver::SetGenerator)
      .def("SetDenominatorCutoff",&IMSRGSolver::SetDenominatorCutoff)
      .def("SetDenominatorDelta",&IMSRGSolver::SetDenominatorDelta)
      .def("SetDenominatorDeltaOrbit",&IMSRGSolver::SetDenominatorDeltaOrbit)
      .def("GetSystemDimension",&IMSRGSolver::GetSystemDimension)
      .def("GetOmega",&IMSRGSolver::GetOmega)
      .def("GetH_s",&IMSRGSolver::GetH_s,return_value_policy<reference_existing_object>())
      .def_readwrite("Eta", &IMSRGSolver::Eta)
   ;



   class_<IMSRGProfiler>("IMSRGProfiler",init<>())
       .def("PrintTimes",&IMSRGProfiler::PrintTimes)
       .def("PrintCounters",&IMSRGProfiler::PrintCounters)
       .def("PrintAll",&IMSRGProfiler::PrintAll)
       .def("PrintMemory",&IMSRGProfiler::PrintMemory)
   ;

   def("TCM_Op",           imsrg_util::TCM_Op);
   def("Trel_Op",           imsrg_util::Trel_Op);
   def("R2CM_Op",          imsrg_util::R2CM_Op);
   def("HCM_Op",           imsrg_util::HCM_Op);
   def("NumberOp",         imsrg_util::NumberOp);
   def("RSquaredOp",       imsrg_util::RSquaredOp);
   def("E0Op",             imsrg_util::E0Op);
   def("AllowedFermi_Op",             imsrg_util::AllowedFermi_Op);
   def("AllowedGamowTeller_Op",             imsrg_util::AllowedGamowTeller_Op);
   def("ElectricMultipoleOp",             imsrg_util::ElectricMultipoleOp);
   def("MagneticMultipoleOp",             imsrg_util::MagneticMultipoleOp);
   def("Isospin2_Op",      imsrg_util::Isospin2_Op);
   def("HO_density",       imsrg_util::HO_density);
   def("GetOccupationsHF", imsrg_util::GetOccupationsHF);
   def("GetOccupations",   imsrg_util::GetOccupations);
   def("GetDensity",       imsrg_util::GetDensity);
   def("CommutatorTest",   imsrg_util::CommutatorTest);
   def("Calculate_p1p2_all",   imsrg_util::Calculate_p1p2_all);
   def("Single_Ref_1B_Density_Matrix", imsrg_util::Single_Ref_1B_Density_Matrix);
   def("Get_Charge_Density", imsrg_util::Get_Charge_Density);

   def("CG",AngMom::CG);
   def("ThreeJ",AngMom::ThreeJ);
   def("SixJ",AngMom::SixJ);
   def("NineJ",AngMom::NineJ);
   def("NormNineJ",AngMom::NormNineJ);
   def("Moshinsky",AngMom::Moshinsky);


}
