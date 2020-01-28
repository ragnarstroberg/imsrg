#include <Python.h>

#include "IMSRG.hh"
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

namespace py = pybind11;


  Orbit MS_GetOrbit(ModelSpace& self, int i){ return self.GetOrbit(i);};
  size_t MS_GetOrbitIndex_Str(ModelSpace& self, std::string s){ return self.GetOrbitIndex(s);};
  TwoBodyChannel MS_GetTwoBodyChannel(ModelSpace& self, int ch){return self.GetTwoBodyChannel(ch);};

  double TB_GetTBME_J(TwoBodyME& self,int j_bra, int j_ket, int a, int b, int c, int d){return self.GetTBME_J(j_bra,j_ket,a,b,c,d);};
  double TB_GetTBME_J_norm(TwoBodyME& self,int j_bra, int j_ket, int a, int b, int c, int d){return self.GetTBME_J_norm(j_bra,j_ket,a,b,c,d);};

  size_t TBCGetLocalIndex(TwoBodyChannel& self, int p, int q){ return self.GetLocalIndex( p, q);};

  void ArmaMatPrint( arma::mat& self){ self.print();};
  void OpSetOneBodyME( Operator& self, int i, int j, double v){self.OneBody(i,j) = v;};

  void MS_SetRef(ModelSpace& self, std::string str){ self.SetReference( str);};

  Operator HF_GetNormalOrderedH(HartreeFock& self){ return self.GetNormalOrderedH();};

//BOOST_PYTHON_MODULE(pyIMSRG)
//PYBIND11_PLUGIN(pyIMSRG)
PYBIND11_MODULE(pyIMSRG, m)
{
  m.doc() = "python bindings for IMSRG code";

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

   py::class_<TwoBodyChannel>(m,"TwoBodyChannel")
      .def(py::init<>())
      .def("GetNumberKets",&TwoBodyChannel::GetNumberKets)
      .def("GetLocalIndex",&TBCGetLocalIndex)
      .def("GetKetIndex",&TwoBodyChannel::GetKetIndex)
   ;

   py::class_<Ket>(m,"Ket")
      .def(py::init<Orbit&,Orbit&>())
      .def_readwrite("p",&Ket::p)
      .def_readwrite("q",&Ket::q)
   ;
   py::class_<Ket3>(m,"Ket3")
      .def(py::init<Orbit&,Orbit&,Orbit&>())
      .def_readwrite("p",&Ket3::p)
      .def_readwrite("q",&Ket3::q)
      .def_readwrite("r",&Ket3::r)
   ;

   py::class_<ModelSpace>(m,"ModelSpace")
      .def(py::init<>())
      .def(py::init<const ModelSpace&>())
      .def(py::init< int, const std::string&>())
      .def(py::init< int, const std::string&, const std::string&>())
      .def(py::init< int,std::vector<std::string>,std::vector<std::string> >())
      .def(py::init< int,std::vector<std::string>,std::vector<std::string>,std::vector<std::string> >())
      .def("SetHbarOmega", &ModelSpace::SetHbarOmega)
      .def("SetTargetMass", &ModelSpace::SetTargetMass)
      .def("SetE3max", &ModelSpace::SetE3max)
      .def("SetEmaxUnocc", &ModelSpace::SetEmaxUnocc)
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
      .def_readwrite("holes",&ModelSpace::holes)
      .def_readwrite("particles",&ModelSpace::particles)
      .def_readwrite("core", &ModelSpace::core)
      .def_readwrite("valence",&ModelSpace::valence)
      .def_readwrite("qspace",&ModelSpace::qspace)
      .def_readwrite("all_orbits",&ModelSpace::all_orbits)
   ;


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
      .def_readwrite("ThreeBody", &Operator::ThreeBody)
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
      .def("Norm", &Operator::Norm)
      .def("OneBodyNorm", &Operator::OneBodyNorm)
      .def("TwoBodyNorm", &Operator::TwoBodyNorm)
      .def("ThreeBodyNorm", &Operator::ThreeBodyNorm)
      .def("SetHermitian", &Operator::SetHermitian)
      .def("SetAntiHermitian", &Operator::SetAntiHermitian)
      .def("SetNonHermitian", &Operator::SetNonHermitian)
      .def("PrintOneBody", &Operator::PrintOneBody)
      .def("PrintTwoBody", &Operator::PrintTwoBody)
      .def("MakeReduced", &Operator::MakeReduced)
      .def("MakeNotReduced", &Operator::MakeNotReduced)
      .def("MakeNormalized", &Operator::MakeNormalized)
      .def("MakeUnNormalized", &Operator::MakeUnNormalized)
      .def("GetParticleRank", &Operator::GetParticleRank)
      .def("SetParticleRank", &Operator::SetParticleRank)
      .def("GetJRank", &Operator::GetJRank)
      .def("GetTRank", &Operator::GetTRank)
      .def("GetParity", &Operator::GetParity)
      .def("GetNumberLegs", &Operator::GetNumberLegs)
      .def("GetE3max", &Operator::GetE3max)
      .def("SetE3max", &Operator::SetE3max)
      .def("PrintTimes", &Operator::PrintTimes)
      .def("Size", &Operator::Size)
      .def("MakeNormalized", &Operator::MakeNormalized)
      .def("MakeUnNormalized", &Operator::MakeUnNormalized)
      .def("SetOneBodyME", &OpSetOneBodyME)
      .def("GetMP2_Energy", &Operator::GetMP2_Energy)
      .def("GetMP3_Energy", &Operator::GetMP3_Energy)
   ;

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

   py::class_<TwoBodyME>(m,"TwoBodyME")
      .def(py::init<>())
      .def("GetTBME_J", TB_GetTBME_J)
      .def("GetTBME_J_norm", TB_GetTBME_J_norm)
   ;

//   py::class_<ThreeBodyME>(m,"ThreeBodyME")
//      .def(py::init<>())
//      .def("SetME", &ThreeBodyME::SetME)
//      .def("GetME", &ThreeBodyME::GetME)
//      .def("GetME_pn", &ThreeBodyME::GetME_pn)
//      .def("RecouplingCoefficient",&ThreeBodyME::RecouplingCoefficient)
//      .def_readonly_static("ABC",&ThreeBodyME::ABC)
//      .def_readonly_static("BCA",&ThreeBodyME::BCA)
//      .def_readonly_static("CAB",&ThreeBodyME::CAB)
//      .def_readonly_static("ACB",&ThreeBodyME::ACB)
//      .def_readonly_static("CBA",&ThreeBodyME::CBA)
//      .def_readonly_static("BAC",&ThreeBodyME::BAC)
//   ;

   py::class_<ThreeBodyMEpn>(m,"ThreeBodyMEpn")
      .def(py::init<>())
      .def("SetME", &ThreeBodyMEpn::SetME)
      .def("GetME", &ThreeBodyMEpn::GetME)
      .def("SetME_pn", &ThreeBodyMEpn::SetME_pn)
      .def("GetME_pn", &ThreeBodyMEpn::GetME_pn)
      .def("RecouplingCoefficient",&ThreeBodyMEpn::RecouplingCoefficient)
      .def("TransformToPN",&ThreeBodyMEpn::TransformToPN)
      .def("SwitchToPN_and_discard",&ThreeBodyMEpn::SwitchToPN_and_discard)
      .def("Print",&ThreeBodyMEpn::Print)
      .def("PrintAll",&ThreeBodyMEpn::PrintAll)
      .def("Erase",&ThreeBodyMEpn::Erase)
      .def_readonly_static("ABC",&ThreeBodyMEpn::ABC)
      .def_readonly_static("BCA",&ThreeBodyMEpn::BCA)
      .def_readonly_static("CAB",&ThreeBodyMEpn::CAB)
      .def_readonly_static("ACB",&ThreeBodyMEpn::ACB)
      .def_readonly_static("CBA",&ThreeBodyMEpn::CBA)
      .def_readonly_static("BAC",&ThreeBodyMEpn::BAC)
   ;

   py::class_<ReadWrite>(m,"ReadWrite")
      .def(py::init<>())
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
      .def("ReadNuShellX_int_iso", &ReadWrite::ReadNuShellX_int_iso)
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
      .def("WriteDaggerOperator",&ReadWrite::WriteDaggerOperator)
      .def("ReadJacobi3NFiles",&ReadWrite::ReadJacobi3NFiles)
      .def("WriteValence3body",&ReadWrite::WriteValence3body)
   ;





   py::class_<HartreeFock>(m,"HartreeFock")
      .def(py::init<Operator&>())
      .def("Solve",&HartreeFock::Solve)
      .def("TransformToHFBasis",&HartreeFock::TransformToHFBasis)
      .def("GetHbare",&HartreeFock::GetHbare)
      .def("GetNormalOrderedH",&HF_GetNormalOrderedH)
      .def("GetOmega",&HartreeFock::GetOmega)
      .def("PrintSPE",&HartreeFock::PrintSPE)
      .def("PrintSPEandWF",&HartreeFock::PrintSPEandWF)
      .def("GetRadialWF_r",&HartreeFock::GetRadialWF_r)
      .def("GetHFPotential",&HartreeFock::GetHFPotential)
      .def("GetAverageHFPotential",&HartreeFock::GetAverageHFPotential)
      .def("GetValence3B",&HartreeFock::GetValence3B)
      .def("FreeVmon",&HartreeFock::FreeVmon)
      .def_static("Vmon3Hash",&HartreeFock::Vmon3Hash)
 // Modifying arguments which were passed by reference causes trouble in python, so instead we bind a lambda function and return a tuple
      .def_static("Vmon3UnHash", [](uint64_t key) { int a,b,c,d,e,f; HartreeFock::Vmon3UnHash(key,a,b,c,d,e,f); return std::make_tuple(a,b,c,d,e,f);}  )
      .def_readonly("EHF",&HartreeFock::EHF)
      .def_readonly("F",&HartreeFock::F) // Fock matrix
      .def_readonly("C",&HartreeFock::C) // Unitary transformation
      .def_readwrite("Vmon3_keys",&HartreeFock::Vmon3_keys)
      .def_readwrite("Vmon3",&HartreeFock::Vmon3)
   ;

   py::class_<HFMBPT,HartreeFock>(m,"HFMBPT")
      .def(py::init<Operator&>())
      .def("GetNaturalOrbitals",&HFMBPT::GetNaturalOrbitals)
      .def("GetNormalOrderedHNAT",&HFMBPT::GetNormalOrderedHNAT)
      .def("PrintSPEandWF",&HFMBPT::PrintSPEandWF)
   ;



   // Define which overloaded version of IMSRGSolver::Transform I want to expose
   Operator (IMSRGSolver::*Transform_ref)(Operator&) = &IMSRGSolver::Transform;

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
      .def("SetHunterGatherer", &IMSRGSolver::SetHunterGatherer)
      .def("AddOperator", &IMSRGSolver::AddOperator)
      .def("GetOperator", &IMSRGSolver::GetOperator)
      .def_readwrite("Eta", &IMSRGSolver::Eta)
   ;



   py::class_<IMSRGProfiler>(m,"IMSRGProfiler")
      .def(py::init<>())
      .def("PrintTimes",&IMSRGProfiler::PrintTimes)
      .def("PrintCounters",&IMSRGProfiler::PrintCounters)
      .def("PrintAll",&IMSRGProfiler::PrintAll)
      .def("PrintMemory",&IMSRGProfiler::PrintMemory)
   ;


   py::class_<Jacobi3BME>(m,"Jacobi3BME")
      .def(py::init<>())
      .def(py::init<int,int,int,int,int>())
      .def("GetDimensionAS",&Jacobi3BME::GetDimensionAS)
      .def("GetDimensionNAS",&Jacobi3BME::GetDimensionNAS)
      .def("GetMatElAS",&Jacobi3BME::GetMatElAS)
      .def("GetMatElNAS",&Jacobi3BME::GetMatElNAS)
      .def("SetEmax",&Jacobi3BME::SetEmax)
      .def("SetE2max",&Jacobi3BME::SetE2max)
      .def("SetE3max",&Jacobi3BME::SetE3max)
      .def("ComputeNAS_MatrixElements",&Jacobi3BME::ComputeNAS_MatrixElements)
      .def("GetLabMatEl",&Jacobi3BME::GetLabMatEl)
      .def("TestReadTcoeffNavratil",&Jacobi3BME::TestReadTcoeffNavratil)
      .def("GetV3mon_all",&Jacobi3BME::GetV3mon_all)
   ;



   py::module Commutator = m.def_submodule("Commutator", "Commutator namespace");
      Commutator.def("Set_BCH_Transform_Threshold", &Commutator::Set_BCH_Transform_Threshold);
      Commutator.def("Set_BCH_Product_Threshold", &Commutator::Set_BCH_Product_Threshold);
      Commutator.def("BCH_Transform", &Commutator::BCH_Transform);
      Commutator.def("BCH_Product", &Commutator::BCH_Product);
      Commutator.def("comm110ss", &Commutator::comm110ss);
      Commutator.def("comm220ss", &Commutator::comm220ss);
      Commutator.def("comm111ss", &Commutator::comm111ss);
      Commutator.def("comm121ss", &Commutator::comm121ss);
      Commutator.def("comm221ss", &Commutator::comm221ss);
      Commutator.def("comm122ss", &Commutator::comm122ss);
      Commutator.def("comm222_pp_hh_221ss", &Commutator::comm222_pp_hh_221ss);
      Commutator.def("comm222_phss", &Commutator::comm222_phss);
      Commutator.def("comm111st", &Commutator::comm111st);
      Commutator.def("comm121st", &Commutator::comm121st);
      Commutator.def("comm122st", &Commutator::comm122st);
      Commutator.def("comm222_pp_hh_221st", &Commutator::comm222_pp_hh_221st);
      Commutator.def("comm222_phst", &Commutator::comm222_phst);


   py::class_<UnitTest>(m,"UnitTest")
//      .def(py::init<>())
      .def(py::init< ModelSpace&>())
      .def("SetRandomSeed",&UnitTest::SetRandomSeed)
      .def("RandomOp",&UnitTest::RandomOp)
      .def("TestCommutators",&UnitTest::TestCommutators)
      .def("TestCommutators3",&UnitTest::TestCommutators3)
      .def("TestDaggerCommutators",&UnitTest::TestDaggerCommutators)
      .def("TestDaggerCommutatorsAlln",&UnitTest::TestDaggerCommutatorsAlln)
      .def("Test3BodyAntisymmetry",&UnitTest::Test3BodyAntisymmetry)
      .def("Test3BodyHermiticity",&UnitTest::Test3BodyHermiticity)
      .def("Test3BodySetGet",&UnitTest::Test3BodySetGet)
   ;


//  py::class_<SymmMatrix<double>>(m,"SymmMatrix")
//     .def(py::init<size_t>())
//     .def(py::init<size_t,int>())
//     .def("Get",&SymmMatrix<double>::Get)
//     .def("Put",&SymmMatrix<double>::Put)
//     .def("FullMatrix",&SymmMatrix<double>::FullMatrix)
//  ;


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
   m.def("HO_Radial_psi",  imsrg_util::HO_Radial_psi);
   m.def("MBPT2_SpectroscopicFactor",  imsrg_util::MBPT2_SpectroscopicFactor);
   m.def("CG",AngMom::CG);
   m.def("ThreeJ",AngMom::ThreeJ);
   m.def("SixJ",AngMom::SixJ);
   m.def("NineJ",AngMom::NineJ);
   m.def("NormNineJ",AngMom::NormNineJ);
   m.def("Moshinsky",AngMom::Moshinsky);
   m.def("TalmiB",AngMom::TalmiB);
   m.def("TalmiI",imsrg_util::TalmiI);
   m.def("Tcoeff",AngMom::Tcoeff);
   m.def("SetUseGooseTank",Commutator::SetUseGooseTank);
   m.def("SetUseIMSRG3",Commutator::SetUseIMSRG3);


//  return m.ptr();

}
