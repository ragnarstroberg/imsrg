#include <Python.h>

#include "IMSRG.hh"
#include <string>

//#include <boost/python/module.hpp>
//#include <boost/python/def.hpp>
//#include <boost/python.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

namespace py = pybind11;

//using namespace boost::python;

  Orbit MS_GetOrbit(ModelSpace& self, int i){ return self.GetOrbit(i);};
  int MS_GetOrbitIndex_Str(ModelSpace& self, string s){ return self.GetOrbitIndex(s);};
  TwoBodyChannel MS_GetTwoBodyChannel(ModelSpace& self, int ch){return self.GetTwoBodyChannel(ch);};

  double TB_GetTBME_J(TwoBodyME& self,int j_bra, int j_ket, int a, int b, int c, int d){return self.GetTBME_J(j_bra,j_ket,a,b,c,d);};
  double TB_GetTBME_J_norm(TwoBodyME& self,int j_bra, int j_ket, int a, int b, int c, int d){return self.GetTBME_J_norm(j_bra,j_ket,a,b,c,d);};

  int TBCGetLocalIndex(TwoBodyChannel& self, int p, int q){ return self.GetLocalIndex( p, q);};

  void ArmaMatPrint( arma::mat& self){ self.print();};
  void OpSetOneBodyME( Operator& self, int i, int j, double v){self.OneBody(i,j) = v;};

  void MS_SetRef(ModelSpace& self, string str){ self.SetReference( str);};

  Operator HF_GetNormalOrderedH(HartreeFock& self){ return self.GetNormalOrderedH();};

//BOOST_PYTHON_MODULE(pyIMSRG)
//PYBIND11_PLUGIN(pyIMSRG)
PYBIND11_MODULE(pyIMSRG, m)
{
//  py::module m("pyIMSRG", "python bindings for IMSRG code");
  m.doc() = "python bindings for IMSRG code";

//   py::class<vector<string> > vector_string("vector_string")
//      .def (vector_indexing_suite< vector<string> >())
//   ;

//   class_<Orbit>("Orbit",init<>())
   py::class_<Orbit>(m,"Orbit")
      .def(py::init<>())
      .def_readwrite("n", &Orbit::n)
      .def_readwrite("l", &Orbit::l)
      .def_readwrite("j2", &Orbit::j2)
      .def_readwrite("tz2", &Orbit::tz2)
      .def_readwrite("occ", &Orbit::occ)
      .def_readwrite("cvq", &Orbit::cvq)
      .def_readwrite("index", &Orbit::index)
   ;

//   class_<TwoBodyChannel>("TwoBodyChannel",init<>())
   py::class_<TwoBodyChannel>(m,"TwoBodyChannel")
      .def(py::init<>())
      .def("GetNumberKets",&TwoBodyChannel::GetNumberKets)
      .def("GetLocalIndex",&TBCGetLocalIndex)
      .def("GetKetIndex",&TwoBodyChannel::GetKetIndex)
   ;

//   class_<ModelSpace>("ModelSpace",init<>())
   py::class_<ModelSpace>(m,"ModelSpace")
      .def(py::init<>())
      .def(py::init<const ModelSpace&>())
      .def(py::init< int, const std::string&>())
      .def(py::init< int, const std::string&, const std::string&>())
      .def(py::init< int,vector<string>,vector<string> >())
      .def(py::init< int,vector<string>,vector<string>,vector<string> >())
      .def("SetHbarOmega", &ModelSpace::SetHbarOmega)
      .def("SetTargetMass", &ModelSpace::SetTargetMass)
      .def("SetE3max", &ModelSpace::SetE3max)
      .def("GetHbarOmega", &ModelSpace::GetHbarOmega)
      .def("GetTargetMass", &ModelSpace::GetTargetMass)
      .def("GetNumberOrbits", &ModelSpace::GetNumberOrbits)
      .def("GetNumberKets", &ModelSpace::GetNumberKets)
      .def("GetOrbit", &MS_GetOrbit)
      .def("GetTwoBodyChannelIndex", &ModelSpace::GetTwoBodyChannelIndex)
      .def("GetTwoBodyChannel", &MS_GetTwoBodyChannel)
      .def("Index2String", &ModelSpace::Index2String)
      .def("ResetFirstPass", &ModelSpace::ResetFirstPass)
      .def("SetReference", &MS_SetRef)
      .def("Init_occ_from_file", &ModelSpace::Init_occ_from_file)
      .def("GetOrbitIndex_fromString", &MS_GetOrbitIndex_Str)
      .def("PreCalculateSixJ", &ModelSpace::PreCalculateSixJ)
      .def_readwrite("core", &ModelSpace::core)
   ;


//   class_<Operator>("Operator",init<>())
   py::class_<Operator>(m,"Operator")
      .def(py::init<>())
      .def(py::init< ModelSpace&>())
      .def(py::init< ModelSpace&,int,int,int,int>())
      .def(py::self += py::self)
      .def(py::self + py::self)
      .def(py::self -= Operator())
      .def(py::self - Operator())
      .def( - py::self)
      .def(py::self *= double())
      .def(py::self * double())
      .def(py::self /= double())
      .def(py::self / double())
      .def(py::self += double())
      .def(py::self + double())
      .def(py::self -= double())
      .def(py::self - double())
      .def_readwrite("ZeroBody", &Operator::ZeroBody)
      .def_readwrite("OneBody", &Operator::OneBody)
      .def_readwrite("TwoBody", &Operator::TwoBody)
      .def("GetOneBody", &Operator::GetOneBody)
      .def("SetOneBody", &Operator::SetOneBody)
      .def("GetTwoBody", &Operator::GetTwoBody)
      .def("SetTwoBody", &Operator::SetTwoBody)
      .def("GetTwoBodyDimension", &Operator::GetTwoBodyDimension)
      .def("ScaleOneBody", &Operator::ScaleOneBody)
      .def("ScaleTwoBody", &Operator::ScaleTwoBody)
      .def("DoNormalOrdering", &Operator::DoNormalOrdering)
      .def("UndoNormalOrdering", &Operator::UndoNormalOrdering)
      .def("SetModelSpace", &Operator::SetModelSpace)
//      .def("CalculateKineticEnergy", &Operator::CalculateKineticEnergy)
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
      .def("MakeReduced", &Operator::MakeReduced)
      .def("MakeNotReduced", &Operator::MakeNotReduced)
      .def("MakeNormalized", &Operator::MakeNormalized)
      .def("MakeUnNormalized", &Operator::MakeUnNormalized)
      .def("GetParticleRank", &Operator::GetParticleRank)
      .def("GetJRank", &Operator::GetJRank)
      .def("GetTRank", &Operator::GetTRank)
      .def("GetParity", &Operator::GetParity)
      .def("GetE3max", &Operator::GetE3max)
      .def("SetE3max", &Operator::SetE3max)
      .def("PrintTimes", &Operator::PrintTimes)
      .def("BCH_Transform", &Operator::BCH_Transform)
      .def("Size", &Operator::Size)
      .def("SetToCommutator", &Operator::SetToCommutator)
      .def("comm110ss", &Operator::comm110ss)
      .def("comm220ss", &Operator::comm220ss)
      .def("comm111ss", &Operator::comm111ss)
      .def("comm121ss", &Operator::comm121ss)
      .def("comm221ss", &Operator::comm221ss)
      .def("comm122ss", &Operator::comm122ss)
      .def("comm222_pp_hh_221ss", &Operator::comm222_pp_hh_221ss)
      .def("comm222_phss", &Operator::comm222_phss)
      .def("comm111st", &Operator::comm111st)
      .def("comm121st", &Operator::comm121st)
      .def("comm122st", &Operator::comm122st)
      .def("comm222_pp_hh_221st", &Operator::comm222_pp_hh_221st)
      .def("comm222_phst", &Operator::comm222_phst)
      .def("MakeNormalized", &Operator::MakeNormalized)
      .def("MakeUnNormalized", &Operator::MakeUnNormalized)
      .def("SetOneBodyME", &OpSetOneBodyME)
   ;

//   class_<arma::mat>("ArmaMat",init<>())
   py::class_<arma::mat>(m,"ArmaMat")
      .def(py::init<>())
      .def("Print",&ArmaMatPrint)
      .def(py::self *= double())
      .def(py::self * double())
      .def(py::self /= double())
      .def(py::self / double())
      .def(py::self += double())
      .def(py::self + double())
      .def(py::self -= double())
      .def(py::self - double())
   ;

//   class_<TwoBodyME>("TwoBodyME",init<>())
   py::class_<TwoBodyME>(m,"TwoBodyME")
      .def(py::init<>())
      .def("GetTBME_J", TB_GetTBME_J)
      .def("GetTBME_J_norm", TB_GetTBME_J_norm)
   ;

//   class_<ReadWrite>("ReadWrite",init<>())
   py::class_<ReadWrite>(m,"ReadWrite")
      .def(py::init<>())
      .def("ReadSettingsFile", &ReadWrite::ReadSettingsFile)
      .def("ReadTBME_Oslo", &ReadWrite::ReadTBME_Oslo)
      .def("ReadBareTBME_Jason", &ReadWrite::ReadBareTBME_Jason)
      .def("ReadBareTBME_Navratil", &ReadWrite::ReadBareTBME_Navratil)
      .def("ReadBareTBME_Darmstadt", &ReadWrite::ReadBareTBME_Darmstadt)
      .def("Read_Darmstadt_3body", &ReadWrite::Read_Darmstadt_3body)
#ifndef NO_HDF5
      .def("Read3bodyHDF5", &ReadWrite::Read3bodyHDF5)
#endif
      .def("Write_me2j", &ReadWrite::Write_me2j)
      .def("Write_me3j", &ReadWrite::Write_me3j)
      .def("WriteTBME_Navratil", &ReadWrite::WriteTBME_Navratil)
      .def("WriteNuShellX_sps", &ReadWrite::WriteNuShellX_sps)
      .def("WriteNuShellX_int", &ReadWrite::WriteNuShellX_int)
      .def("WriteNuShellX_op", &ReadWrite::WriteNuShellX_op)
      .def("ReadNuShellX_int", &ReadWrite::ReadNuShellX_int)
      .def("WriteAntoine_int", &ReadWrite::WriteAntoine_int)
      .def("WriteAntoine_input", &ReadWrite::WriteAntoine_input)
      .def("WriteOperator", &ReadWrite::WriteOperator)
      .def("WriteOperatorHuman", &ReadWrite::WriteOperatorHuman)
      .def("ReadOperator", &ReadWrite::ReadOperator)
      .def("ReadOperatorHuman", &ReadWrite::ReadOperatorHuman)
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
      .def("ReadOperator_Nathan",&ReadWrite::ReadOperator_Nathan)
      .def("ReadTensorOperator_Nathan",&ReadWrite::ReadTensorOperator_Nathan)
      .def("ReadRelCMOpFromJavier",&ReadWrite::ReadRelCMOpFromJavier)
      .def("Set3NFormat",&ReadWrite::Set3NFormat)
      .def_readwrite("InputParameters", &ReadWrite::InputParameters)
   ;





//   class_<HartreeFock>("HartreeFock",init<Operator&>())
   py::class_<HartreeFock>(m,"HartreeFock")
      .def(py::init<Operator&>())
      .def("Solve",&HartreeFock::Solve)
      .def("TransformToHFBasis",&HartreeFock::TransformToHFBasis)
      .def("GetHbare",&HartreeFock::GetHbare)
      .def("GetNormalOrderedH",&HF_GetNormalOrderedH)
      .def("GetOmega",&HartreeFock::GetOmega)
      .def("PrintSPE",&HartreeFock::PrintSPE)
      .def("GetRadialWF_r",&HartreeFock::GetRadialWF_r)
      .def_readonly("EHF",&HartreeFock::EHF)
      .def_readonly("C",&HartreeFock::C)
   ;

   // Define which overloaded version of IMSRGSolver::Transform I want to expose
   Operator (IMSRGSolver::*Transform_ref)(Operator&) = &IMSRGSolver::Transform;

//   class_<IMSRGSolver>("IMSRGSolver",init<Operator&>())
   py::class_<IMSRGSolver>(m,"IMSRGSolver")
      .def(py::init<Operator&>())
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
      .def("SetOmega",&IMSRGSolver::SetOmega)
//      .def("GetH_s",&IMSRGSolver::GetH_s,return_value_policy<reference_existing_object>())
      .def("GetH_s",&IMSRGSolver::GetH_s)
      .def("SetMagnusAdaptive",&IMSRGSolver::SetMagnusAdaptive)
      .def("SetReadWrite", &IMSRGSolver::SetReadWrite)
      .def_readwrite("Eta", &IMSRGSolver::Eta)
   ;



//   class_<IMSRGProfiler>("IMSRGProfiler",init<>())
   py::class_<IMSRGProfiler>(m,"IMSRGProfiler")
      .def(py::init<>())
      .def("PrintTimes",&IMSRGProfiler::PrintTimes)
      .def("PrintCounters",&IMSRGProfiler::PrintCounters)
      .def("PrintAll",&IMSRGProfiler::PrintAll)
      .def("PrintMemory",&IMSRGProfiler::PrintMemory)
   ;





   m.def("TCM_Op",           imsrg_util::TCM_Op);
   m.def("Trel_Op",           imsrg_util::Trel_Op);
   m.def("R2CM_Op",          imsrg_util::R2CM_Op);
   m.def("HCM_Op",           imsrg_util::HCM_Op);
   m.def("NumberOp",         imsrg_util::NumberOp);
   m.def("RSquaredOp",       imsrg_util::RSquaredOp);
   m.def("RpSpinOrbitCorrection", imsrg_util::RpSpinOrbitCorrection);
   m.def("E0Op",             imsrg_util::E0Op);
   m.def("AllowedFermi_Op",             imsrg_util::AllowedFermi_Op);
   m.def("AllowedGamowTeller_Op",             imsrg_util::AllowedGamowTeller_Op);
   m.def("ElectricMultipoleOp",             imsrg_util::ElectricMultipoleOp);
   m.def("MagneticMultipoleOp",             imsrg_util::MagneticMultipoleOp);
   m.def("Sigma_Op", imsrg_util::Sigma_Op);
   m.def("Isospin2_Op",      imsrg_util::Isospin2_Op);
   m.def("LdotS_Op",         imsrg_util::LdotS_Op);
   m.def("HO_density",       imsrg_util::HO_density);
   m.def("HO_Radial_psi",    imsrg_util::HO_Radial_psi);
   m.def("GetOccupationsHF", imsrg_util::GetOccupationsHF);
   m.def("GetOccupations",   imsrg_util::GetOccupations);
   m.def("GetDensity",       imsrg_util::GetDensity);
   m.def("CommutatorTest",   imsrg_util::CommutatorTest);
   m.def("Calculate_p1p2_all",   imsrg_util::Calculate_p1p2_all);
   m.def("Single_Ref_1B_Density_Matrix", imsrg_util::Single_Ref_1B_Density_Matrix);
   m.def("Get_Charge_Density", imsrg_util::Get_Charge_Density);
   m.def("Embed1BodyIn2Body",  imsrg_util::Embed1BodyIn2Body);
   m.def("RadialIntegral",     imsrg_util::RadialIntegral);
   m.def("RadialIntegral_RpowK",     imsrg_util::RadialIntegral_RpowK);
   m.def("FrequencyConversionCoeff", imsrg_util::FrequencyConversionCoeff);
   m.def("OperatorFromString", imsrg_util::OperatorFromString);

   m.def("CG",AngMom::CG);
   m.def("ThreeJ",AngMom::ThreeJ);
   m.def("SixJ",AngMom::SixJ);
   m.def("NineJ",AngMom::NineJ);
   m.def("NormNineJ",AngMom::NormNineJ);
   m.def("Moshinsky",AngMom::Moshinsky);


//  return m.ptr();

}
