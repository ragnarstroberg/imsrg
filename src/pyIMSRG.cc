#include <Python.h>

#include "IMSRG.hh"
#include <string>
#include <sstream>
#include <vector>
#include "version.hh"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

namespace py = pybind11;

//  Orbit MS_GetOrbit(ModelSpace& self, int i){ return self.GetOrbit(i);};
//  size_t MS_GetOrbitIndex_Str(ModelSpace& self, std::string s){ return self.GetOrbitIndex(s);};
//  TwoBodyChannel MS_GetTwoBodyChannel(ModelSpace& self, int ch){return self.GetTwoBodyChannel(ch);};

//  double TB_GetTBME_J(TwoBodyME& self,int j_bra, int j_ket, int a, int b, int c, int d){return self.GetTBME_J(j_bra,j_ket,a,b,c,d);};
//  double TB_GetTBME_J_norm(TwoBodyME& self,int j_bra, int j_ket, int a, int b, int c, int d){return self.GetTBME_J_norm(j_bra,j_ket,a,b,c,d);};

//  size_t TBCGetLocalIndex(TwoBodyChannel& self, int p, int q){ return self.GetLocalIndex( p, q);};

//  void ArmaMatPrint( arma::mat& self){ self.print();};
//  void OpSetOneBodyME( Operator& self, int i, int j, double v){self.OneBody(i,j) = v;};

//  void MS_SetRef(ModelSpace& self, std::string str){ self.SetReference( str);};
//  void MS_SetRef(ModelSpace& self, const std::set<index_t>& ref){ self.SetReference( ref);};

//  Operator HF_GetNormalOrderedH(HartreeFock& self){ return self.GetNormalOrderedH();};
//  Operator HF_GetNormalOrderedH(HartreeFock& self, int particle_rank=2){ return self.GetNormalOrderedH(particle_rank);};

// BOOST_PYTHON_MODULE(pyIMSRG)
// PYBIND11_PLUGIN(pyIMSRG)
PYBIND11_MODULE(pyIMSRG, m)
{
      m.doc() = "python bindings for IMSRG code";

      py::class_<Orbit>(m, "Orbit")
          .def(py::init<>())
          .def_readwrite("n", &Orbit::n)
          .def_readwrite("l", &Orbit::l)
          .def_readwrite("j2", &Orbit::j2)
          .def_readwrite("tz2", &Orbit::tz2)
          .def_readwrite("occ", &Orbit::occ)
          .def_readwrite("cvq", &Orbit::cvq)
          .def_readwrite("index", &Orbit::index);

      py::class_<TwoBodyChannel>(m, "TwoBodyChannel")
          .def(py::init<>())
          .def("GetNumberKets", &TwoBodyChannel::GetNumberKets)
          //      .def("GetLocalIndex",&TBCGetLocalIndex)
          .def("GetLocalIndex", [](TwoBodyChannel &self, int p, int q)
               { return self.GetLocalIndex(p, q); })
          .def("GetKetIndex", &TwoBodyChannel::GetKetIndex)
          .def("GetKet", [](TwoBodyChannel &self, int i)
               { return self.GetKet(i); })
          .def("GetKetIndex_pp", [](TwoBodyChannel &self)
               { auto& x=self.GetKetIndex_pp(); std::vector<size_t> v(x.begin(),x.end()); return v; })
          .def("GetKetIndex_hh", [](TwoBodyChannel &self)
               { auto& x=self.GetKetIndex_hh(); std::vector<size_t> v(x.begin(),x.end()); return v; })
          .def("GetKetIndex_ph", [](TwoBodyChannel &self)
               { auto& x=self.GetKetIndex_ph(); std::vector<size_t> v(x.begin(),x.end()); return v; })
          .def_readwrite("J", &TwoBodyChannel::J)
          .def_readwrite("parity", &TwoBodyChannel::parity)
          .def_readwrite("Tz", &TwoBodyChannel::Tz);

      py::class_<ThreeBodyChannel>(m, "ThreeBodyChannel")
          .def(py::init<>())
          .def("GetNumber3bKets", &ThreeBodyChannel::GetNumber3bKets)
          .def("GetLocalIndex", &ThreeBodyChannel::GetLocalIndex, py::arg("p"), py::arg("q"), py::arg("r"), py::arg("Jpq"))
          .def("GetKet", [](ThreeBodyChannel &self, int i)
               { return self.GetKet(i); })
          .def_readwrite("twoJ", &ThreeBodyChannel::twoJ)
          .def_readwrite("parity", &ThreeBodyChannel::parity)
          .def_readwrite("twoTz", &ThreeBodyChannel::twoTz);

      py::class_<Ket>(m, "Ket")
          .def(py::init<Orbit &, Orbit &>())
          .def_readwrite("p", &Ket::p)
          .def_readwrite("q", &Ket::q);
      py::class_<Ket3>(m, "Ket3")
          .def(py::init<Orbit &, Orbit &, Orbit &>())
          .def_readwrite("p", &Ket3::p)
          .def_readwrite("q", &Ket3::q)
          .def_readwrite("r", &Ket3::r)
          .def_readwrite("Jpq", &Ket3::Jpq);

      py::class_<ModelSpace>(m, "ModelSpace")
          .def(py::init<>())
          .def(py::init<const ModelSpace &>())
          .def(py::init<int, const std::string &>(), py::arg("emax"), py::arg("reference"))
          .def(py::init<int, const std::string &, const std::string &>(), py::arg("emax"), py::arg("reference"), py::arg("valence"))
          .def(py::init<int, std::vector<std::string>, std::vector<std::string>>(), py::arg("emax"), py::arg("hole_list"), py::arg("valence_list"))
          .def(py::init<int, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>(), py::arg("emax"), py::arg("hole_list"), py::arg("core_list"), py::arg("valence_list"))
          .def("SetHbarOmega", &ModelSpace::SetHbarOmega)
          .def("SetTargetMass", &ModelSpace::SetTargetMass)
          .def("SetTargetZ", &ModelSpace::SetTargetZ)
          .def(
              "AddOrbit", [](ModelSpace &self, int n, int l, int j2, int tz2, double occ, int cvq)
              { self.AddOrbit(n, l, j2, tz2, occ, cvq); },
              py::arg("n"), py::arg("l"), py::arg("j2"), py::arg("tz2"), py::arg("occ"), py::arg("cvq"))
          .def("SetupKets", &ModelSpace::SetupKets)
          .def("Setup3bKets", &ModelSpace::Setup3bKets)
          .def("SetOcc", &ModelSpace::SetOcc, py::arg("n"), py::arg("l"), py::arg("j2"), py::arg("tz2"), py::arg("occ"))
          .def("SetOccNAT", &ModelSpace::SetOccNAT, py::arg("n"), py::arg("l"), py::arg("j2"), py::arg("tz2"), py::arg("occ_nat"))
          .def("SetEmax", &ModelSpace::SetEmax)
          .def("SetE2max", &ModelSpace::SetE2max)
          .def("SetE3max", &ModelSpace::SetE3max)
          .def("SetdE3max", &ModelSpace::SetdE3max)
          .def("SetLmax", &ModelSpace::SetLmax)
          .def("SetEmaxUnocc", &ModelSpace::SetEmaxUnocc)
          .def("SetEmax3Body", &ModelSpace::SetEmax3Body)
          .def("FindEFermi", &ModelSpace::FindEFermi)
          .def("GetHbarOmega", &ModelSpace::GetHbarOmega)
          .def("GetTargetMass", &ModelSpace::GetTargetMass)
          .def("GetTargetZ", &ModelSpace::GetTargetZ)
          .def("GetAref", &ModelSpace::GetAref)
          .def("GetZref", &ModelSpace::GetZref)
          .def("GetNumberOrbits", &ModelSpace::GetNumberOrbits)
          .def("GetNumberKets", &ModelSpace::GetNumberKets)
          .def("GetNumberTwoBodyChannels", &ModelSpace::GetNumberTwoBodyChannels)
          .def("GetNumberThreeBodyChannels", &ModelSpace::GetNumberThreeBodyChannels)
          //      .def("GetOrbit", &MS_GetOrbit)
          .def("GetOrbit", [](ModelSpace &self, int i)
               { return self.GetOrbit(i); })
          .def("GetTwoBodyChannelIndex", &ModelSpace::GetTwoBodyChannelIndex)
          .def("GetTwoBodyChannel", [](ModelSpace &self, int ch)
               { return self.GetTwoBodyChannel(ch); })
          .def("GetThreeBodyChannel", &ModelSpace::GetThreeBodyChannel)
          .def("GetThreeBodyChannelIndex", &ModelSpace::GetThreeBodyChannelIndex, py::arg("twoJ"), py::arg("parity"), py::arg("twoTz"))
          .def("Index2String", &ModelSpace::Index2String)
          .def("ResetFirstPass", &ModelSpace::ResetFirstPass)
          //      .def("SetReference", &MS_SetRef)
          .def("SetReference", [](ModelSpace &self, const std::set<index_t> &ref)
               { self.SetReference(ref); })
          .def("SetReferenceStr", [](ModelSpace &self, std::string s)
               { self.SetReference(s); })
          .def("Init_occ_from_file", &ModelSpace::Init_occ_from_file)
          .def("InitSingleSpecies", &ModelSpace::InitSingleSpecies)
          .def(
              "GetOrbitIndex", [](ModelSpace &self, int n, int l, int j, int tz)
              { return self.GetOrbitIndex(n, l, j, tz); },
              py::arg("n"), py::arg("l"), py::arg("j2"), py::arg("tz2"))
          .def(
              "GetOrbitIndex_fromString", [](ModelSpace &self, std::string s)
              { return self.GetOrbitIndex(s); },
              py::arg("orbstring"))
          .def(
              "GetOneBodyChannels", [](ModelSpace &self, int l, int j, int tz)
              { return self.OneBodyChannels.at({l, j, tz}); },
              py::arg("l"), py::arg("j2"), py::arg("tz2"))
          //      .def("GetOrbitIndex_fromString", &MS_GetOrbitIndex_Str)
          .def("PreCalculateSixJ", &ModelSpace::PreCalculateSixJ)
          .def("SetScalarFirstPass", &ModelSpace::SetScalarFirstPass)
          .def("SetScalar3bFirstPass", &ModelSpace::SetScalar3bFirstPass)
          .def("ClearVectors", &ModelSpace::ClearVectors)
          .def("Print", &ModelSpace::Print)
          .def_readwrite("holes", &ModelSpace::holes)
          .def_readwrite("particles", &ModelSpace::particles)
          .def_readwrite("core", &ModelSpace::core)
          .def_readwrite("valence", &ModelSpace::valence)
          .def_readwrite("qspace", &ModelSpace::qspace)
          .def_readwrite("all_orbits", &ModelSpace::all_orbits);

      py::class_<Operator>(m, "Operator")
          .def(py::init<>())
          .def(py::init<ModelSpace &>())
          .def(py::init<Operator &>())
          .def(py::init<ModelSpace &, int, int, int, int>(), py::arg("modelspace"), py::arg("j_rank"), py::arg("t_rank"), py::arg("parity"), py::arg("particle_rank"))
          .def(py::self += py::self)
          .def(py::self + py::self)
          .def(py::self -= Operator())
          .def(py::self - Operator())
          .def(-py::self)
          .def(py::self *= double())
          .def(py::self * double())
          .def(double() * py::self)
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
          .def("EraseOneBody", &Operator::EraseOneBody)
          .def("EraseTwoBody", &Operator::EraseTwoBody)
          .def("DoNormalOrdering", &Operator::DoNormalOrdering)
          .def("DoNormalOrderingCore", &Operator::DoNormalOrderingCore)
          .def("UndoNormalOrdering", &Operator::UndoNormalOrdering)
          .def("SetModelSpace", &Operator::SetModelSpace)
          .def("Truncate", &Operator::Truncate)
          .def("Norm", &Operator::Norm)
          .def("OneBodyNorm", &Operator::OneBodyNorm)
          .def("TwoBodyNorm", &Operator::TwoBodyNorm)
          .def("ThreeBodyNorm", &Operator::ThreeBodyNorm)
          .def("SetHermitian", &Operator::SetHermitian)
          .def("SetAntiHermitian", &Operator::SetAntiHermitian)
          .def("SetNonHermitian", &Operator::SetNonHermitian)
          .def("IsHermitian", &Operator::IsHermitian)
          .def("IsAntiHermitian", &Operator::IsAntiHermitian)
          .def("PrintOneBody", &Operator::PrintOneBody)
          .def("PrintTwoBody", [](Operator &self)
               { self.PrintTwoBody(); })
          .def("PrintTwoBody_ch", [](Operator &self, int ch)
               { self.PrintTwoBody(ch); })
          .def("PrintTwoBody_chch", [](Operator &self, int ch_bra, int ch_ket)
               { self.PrintTwoBody(ch_bra, ch_ket); })
          //      .def("PrintTwoBody_ch", &Operator::PrintTwoBody)
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
          .def("GetOneBodyChannel", &Operator::GetOneBodyChannel, py::arg("l"), py::arg("j2"), py::arg("tz2"))
          //      .def("SetOneBodyME", &OpSetOneBodyME)
          .def("SetOneBodyME", [](Operator &self, int i, int j, double v)
               { self.OneBody(i, j) = v; })
          .def("GetMP2_Energy", &Operator::GetMP2_Energy)
          .def("GetMP2_3BEnergy", &Operator::GetMP2_Energy)
          .def("GetMP3_Energy", &Operator::GetMP3_Energy)
          .def("GetPPHH_Ladders", &Operator::GetPPHH_Ladders)
          .def(
              "ReadBinary", [](Operator &self, std::string fname)
              { std::ifstream ifs(fname,std::ios::binary);  self.ReadBinary(ifs); },
              py::arg("filename"))
          .def(
              "WriteBinary", [](Operator &self, std::string fname)
              { std::ofstream ofs(fname,std::ios::binary);  self.WriteBinary(ofs); },
              py::arg("filename"))
          //      .def("IsospinProject", &Operator::IsospinProject)
          ;

      py::class_<arma::mat>(m, "ArmaMat")
          .def(py::init<>())
          .def(
              "zeros", [](arma::mat &self, int nrows, int ncols)
              { self.zeros(nrows, ncols); },
              py::arg("nrows"), py::arg("ncols"))
          .def("Print", [](arma::mat &self)
               { self.print(); }) //   &ArmaMatPrint)
          .def("__str__", [](arma::mat &self)
               { std::ostringstream oss; oss << self; return oss.str(); }) //   &ArmaMatPrint)
          .def(
              "save", [](arma::mat &self, std::string fname)
              { self.save(fname); },
              py::arg("filename"))
          .def(
              "load", [](arma::mat &self, std::string fname)
              { self.load(fname); },
              py::arg("filename"))
          //      .def("t", &arma::mat::t) // transpose
          .def("t", [](arma::mat &self)
               {arma::mat x = self.t(); return x; }) // transpose
          .def(py::self *= double())
          //      .def(py::self * double())
          //      .def(double() * py::self)
          //      .def(double() * py::self, [](double x, arma::mat& self){arma::mat out = x * self; return out;} )
          .def(py::self /= double())
          .def(py::self / double())
          //      .def(py::self += ArmaMat())
          //      .def(py::self + ArmaMat())
          //      .def(py::self -= ArmaMat())
          //      .def(py::self - ArmaMat())
          .def(
              "__mul__", [](const arma::mat &A, const arma::mat &B)
              {arma::mat C = A * B; return C; },
              py::is_operator())
          .def(
              "__mul__", [](const arma::mat &B, float A)
              {arma::mat C = A * B; return C; },
              py::is_operator())
          //      .def("__mul__", [](float A, const arma::mat& B){arma::mat C = A * B; return C;}, py::is_operator() )
          .def(
              "__add__", [](const arma::mat &A, const arma::mat &B)
              {arma::mat C = A + B; return C; },
              py::is_operator())
          .def(
              "__sub__", [](const arma::mat &A, const arma::mat &B)
              {arma::mat C = A - B; return C; },
              py::is_operator())
          .def(
              "__call__", [](arma::mat &self, const int i, const int j)
              { return &self(i, j); },
              py::is_operator())
          .def(
              "Set", [](arma::mat &self, const int i, const int j, double x)
              { self(i, j) = x; },
              py::arg("i"), py::arg("j"), py::arg("matel"))
          .def("Getn_rows", [](arma::mat &self)
               { return self.n_rows; })
          .def("Getn_cols", [](arma::mat &self)
               { return self.n_cols; })
          .def("Schur_Prod", [](arma::mat &self, arma::mat &other)
               { arma::mat out = self % other;return out; })
          .def("Norm", [](arma::mat &self)
               { return arma::norm(self, "fro"); })
          .def("trace", [](arma::mat &self)
               { double t =arma::trace(self); return t; })
          .def("sum", [](arma::mat &self)
               {double s= arma::accu(self); return s; });

      py::class_<TwoBodyME>(m, "TwoBodyME")
          .def(py::init<>())
          //      .def("GetTBME_J", TB_GetTBME_J)
          //      .def("GetTBME_J_norm", TB_GetTBME_J_norm)
          .def("GetTBME_J", [](TwoBodyME &self, int Jbra, int Jket, int a, int b, int c, int d)
               { return self.GetTBME_J(Jbra, Jket, a, b, c, d); })
          .def("GetTBME_J_norm", [](TwoBodyME &self, int Jbra, int Jket, int a, int b, int c, int d)
               { return self.GetTBME_J_norm(Jbra, Jket, a, b, c, d); })
          .def(
              "GetTBMEmonopole", [](TwoBodyME &self, int a, int b, int c, int d)
              { return self.GetTBMEmonopole(a, b, c, d); },
              py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"))
          .def("GetTBME_norm", [](TwoBodyME &self, int ch_bra, int ch_ket, int a, int b, int c, int d)
               { return self.GetTBME_norm(ch_bra, ch_ket, a, b, c, d); })
          .def(
              "GetTBMEmonopole_norm", [](TwoBodyME &self, int a, int b, int c, int d)
              { return self.GetTBMEmonopole_norm(a, b, c, d); },
              py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"))
          .def(
              "GetChannelMatrix", [](TwoBodyME &self, int J, int p, int Tz)
              { size_t ch = self.modelspace->GetTwoBodyChannelIndex(J,p,Tz); return self.GetMatrix(ch,ch); },
              py::arg("J"), py::arg("parity"), py::arg("Tz"))
          .def("PrintAll", [](TwoBodyME &self)
               { for (auto& it : self.MatEl){ if (it.second.n_rows>0) { std::cout << it.first[0] << " " << it.first[1] << std::endl << it.second << std::endl;};  } ; })
          .def("PrintMatrix", &TwoBodyME::PrintMatrix, py::arg("ch_bra"), py::arg("ch_ket"))
          .def("Erase", &TwoBodyME::Erase)
          .def("GetTBMEnorm_chij", [](TwoBodyME &self, int ch_bra, int ch_ket, size_t ibra, size_t iket)
               { return self.GetTBME_norm(ch_bra, ch_ket, ibra, iket); })
          .def("SetTBME_chij", [](TwoBodyME &self, int ch_bra, int ch_ket, size_t ibra, size_t iket, double tbme)
               { self.SetTBME(ch_bra, ch_ket, ibra, iket, tbme); })
          .def("Norm", &TwoBodyME::Norm)
          .def(py::self *= double())
          .def(double() * py::self)
          .def(py::self * double())
          .def(py::self + TwoBodyME())
          .def(py::self += TwoBodyME())
          .def(py::self - TwoBodyME())
          .def(py::self -= TwoBodyME());

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

      //   py::class_<ThreeBodyMEpn>(m,"ThreeBodyMEpn")
      py::class_<ThreeBodyME>(m, "ThreeBodyME")
          .def(py::init<>())
          //      .def("SetME", &ThreeBodyMEpn::SetME)
          //      .def("GetME", &ThreeBodyME::GetME)
          .def(
              "GetME_iso", [](ThreeBodyME &self, int Jab, int Jde, int twoJ, int tab, int tde, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f)
              { return self.GetME_iso(Jab, Jde, twoJ, tab, tde, twoTabc, twoTdef, a, b, c, d, e, f); },
              py::arg("Jab"), py::arg("Jde"), py::arg("twoJ"), py::arg("tab"), py::arg("tde"), py::arg("twoTabc"), py::arg("twoTdef"), py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"), py::arg("e"), py::arg("f"))
          //      .def("SetME_pn", &ThreeBodyME::SetME_pn)
          .def("GetME_pn", &ThreeBodyME::GetME_pn)
          .def("SetME_pn_ch", &ThreeBodyME::SetME_pn_ch) // Hopefully not a bad idea to expose this...
          .def("GetME_pn_no2b", &ThreeBodyME::GetME_pn_no2b)
          .def("RecouplingCoefficient", &ThreeBodyME::RecouplingCoefficient)
          .def("TransformToPN", &ThreeBodyME::TransformToPN)
          .def("SwitchToPN_and_discard", &ThreeBodyME::SwitchToPN_and_discard)
          //      .def("Print",&ThreeBodyME::Print)
          //      .def("PrintAll",&ThreeBodyME::PrintAll)
          .def("Erase", &ThreeBodyME::Erase)
          .def("SetMode", &ThreeBodyME::SetMode)
          .def("ReadFile", &ThreeBodyME::ReadFile, py::arg("string_inputs"), py::arg("int_inputs"))
          .def(py::self += ThreeBodyME(), py::is_operator())
          .def(py::self *= double())
          //      .def_readonly_static("ABC",&ThreeBodyME::ABC)
          //      .def_readonly_static("BCA",&ThreeBodyME::BCA)
          //      .def_readonly_static("CAB",&ThreeBodyME::CAB)
          //      .def_readonly_static("ACB",&ThreeBodyME::ACB)
          //      .def_readonly_static("CBA",&ThreeBodyME::CBA)
          //      .def_readonly_static("BAC",&ThreeBodyME::BAC)
          ;

      py::class_<ReadWrite>(m, "ReadWrite")
          .def(py::init<>())
          .def("ReadTBME_Oslo", &ReadWrite::ReadTBME_Oslo)
          .def("ReadTBME_OakRidge", &ReadWrite::ReadTBME_OakRidge, py::arg("spname"), py::arg("tbmename"), py::arg("H"), py::arg("tbme_format") = "ascii")
          .def("ReadBareTBME_Jason", &ReadWrite::ReadBareTBME_Jason)
          .def("ReadBareTBME_Navratil", &ReadWrite::ReadBareTBME_Navratil)
          .def("ReadBareTBME_Darmstadt", &ReadWrite::ReadBareTBME_Darmstadt, py::arg("filename"), py::arg("H"), py::arg("e1max"), py::arg("e2max"), py::arg("lmax"))
          .def("Read_Darmstadt_3body", &ReadWrite::Read_Darmstadt_3body, py::arg("filename"), py::arg("H"), py::arg("e1max"), py::arg("e2max"), py::arg("e3max"))
#ifndef NO_HDF5
          .def("Read3bodyHDF5", &ReadWrite::Read3bodyHDF5)
#endif
          .def("Write_me2j", &ReadWrite::Write_me2j)
          .def("Write_me3j", &ReadWrite::Write_me3j)
          .def("WriteTBME_Navratil", &ReadWrite::WriteTBME_Navratil)
          .def("WriteNuShellX_sps", &ReadWrite::WriteNuShellX_sps, py::arg("op"), py::arg("filename"))
          .def("WriteNuShellX_int", &ReadWrite::WriteNuShellX_int, py::arg("op"), py::arg("filename"))
          .def("WriteNuShellX_op", &ReadWrite::WriteNuShellX_op, py::arg("op"), py::arg("filename"))
          .def("ReadNuShellX_int", &ReadWrite::ReadNuShellX_int, py::arg("op"), py::arg("filename"))
          .def("ReadNuShellX_int_iso", &ReadWrite::ReadNuShellX_int_iso, py::arg("op"), py::arg("filename"))
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
          .def("WriteTokyo", &ReadWrite::WriteTokyo, py::arg("op"), py::arg("filename"), py::arg("mode"))
          .def("WriteTensorTokyo", &ReadWrite::WriteTensorTokyo, py::arg("filename"), py::arg("op"))
          .def(
              "ReadTokyo", [](ReadWrite &self, std::string s, Operator &op)
              { self.ReadTokyo(s, op); },
              py::arg("file_in"), py::arg("op"))
          .def("WriteOneBody_Oslo", &ReadWrite::WriteOneBody_Oslo)
          .def("WriteTwoBody_Oslo", &ReadWrite::WriteTwoBody_Oslo)
          .def("SetCoMCorr", &ReadWrite::SetCoMCorr)
          .def("ReadTwoBodyEngel", &ReadWrite::ReadTwoBodyEngel)
          .def("ReadOperator_Nathan", &ReadWrite::ReadOperator_Nathan)
          .def("ReadTensorOperator_Nathan", &ReadWrite::ReadTensorOperator_Nathan)
          .def("ReadRelCMOpFromJavier", &ReadWrite::ReadRelCMOpFromJavier)
          .def("Set3NFormat", &ReadWrite::Set3NFormat)
          .def("WriteDaggerOperator", &ReadWrite::WriteDaggerOperator)
          .def("ReadJacobi3NFiles", &ReadWrite::ReadJacobi3NFiles)
          .def("WriteValence3body", &ReadWrite::WriteValence3body)
          .def("SetScratchDir", &ReadWrite::SetScratchDir)
          .def("GetScratchDir", &ReadWrite::GetScratchDir)
          .def("CopyFile", &ReadWrite::CopyFile, py::arg("filein"), py::arg("fileout"))
          //      .def("WriteOmega",&ReadWrite::WriteOmega, py::arg("basename"),py::arg("scratch_dir"),py::arg("nOmegas"))
          ;

      py::class_<HartreeFock>(m, "HartreeFock")
          .def(py::init<Operator &>())
          .def("Solve", &HartreeFock::Solve)
          .def("TransformToHFBasis", &HartreeFock::TransformToHFBasis)
          .def("GetHbare", &HartreeFock::GetHbare)
          //      .def("GetNormalOrderedH",&HF_GetNormalOrderedH)
          //      .def("GetNormalOrderedH",&HF_GetNormalOrderedH, py::arg("particle_rank")=2 )
          .def(
              "GetNormalOrderedH", [](HartreeFock &self, int pRank)
              { return self.GetNormalOrderedH(pRank); },
              py::arg("particle_rank") = 2)
          .def(
              "GetNormalOrderedH_Cin", [](HartreeFock &self, arma::mat &C, int pRank)
              { return self.GetNormalOrderedH(C, pRank); },
              py::arg("C"), py::arg("particle_rank") = 2)
          .def("GetOmega", &HartreeFock::GetOmega)
          .def("PrintSPE", &HartreeFock::PrintSPE)
          .def("PrintSPEandWF", &HartreeFock::PrintSPEandWF)
          .def("GetRadialWF_r", &HartreeFock::GetRadialWF_r)
          .def("GetHFPotential", &HartreeFock::GetHFPotential)
          .def("GetAverageHFPotential", &HartreeFock::GetAverageHFPotential)
          .def("GetValence3B", &HartreeFock::GetValence3B)
          .def("FreeVmon", &HartreeFock::FreeVmon)
          .def("UpdateDensityMatrix", &HartreeFock::UpdateDensityMatrix)
          .def("UpdateF", &HartreeFock::UpdateF)
          .def("BuildMonopoleV", &HartreeFock::BuildMonopoleV)
          .def("CalcEHF", &HartreeFock::CalcEHF)
          .def("PrintEHF", &HartreeFock::PrintEHF)
          .def("FillLowestOrbits", &HartreeFock::FillLowestOrbits)
          .def("DiscardNO2Bfrom3N", &HartreeFock::DiscardNO2Bfrom3N)
          .def("FreezeOccupations", &HartreeFock::FreezeOccupations)
          .def("UnFreezeOccupations", &HartreeFock::UnFreezeOccupations)
          .def_static("Vmon3Hash", &HartreeFock::Vmon3Hash)
          // Modifying arguments which were passed by reference causes trouble in python, so instead we bind a lambda function and return a tuple
          .def_static("Vmon3UnHash", [](uint64_t key)
                      { int a,b,c,d,e,f; HartreeFock::Vmon3UnHash(key,a,b,c,d,e,f); return std::make_tuple(a,b,c,d,e,f); })
          .def_readonly("EHF", &HartreeFock::EHF)
          .def_readonly("F", &HartreeFock::F)     // Fock matrix
          .def_readonly("rho", &HartreeFock::rho) // density matrix
                                                  //      .def_readonly("C",&HartreeFock::C) // Unitary transformation
          .def_readwrite("C", &HartreeFock::C)    // Unitary transformation
          .def_readwrite("Vmon3_keys", &HartreeFock::Vmon3_keys)
          .def_readwrite("Vmon3", &HartreeFock::Vmon3);

      py::class_<HFMBPT, HartreeFock>(m, "HFMBPT")
          .def(py::init<Operator &>())
          .def("UseNATOccupations", &HFMBPT::UseNATOccupations)
          .def("GetNaturalOrbitals", &HFMBPT::GetNaturalOrbitals)
          .def("TransformHOToNATBasis", &HFMBPT::TransformHOToNATBasis)
          .def("TransformHFToNATBasis", &HFMBPT::TransformHFToNATBasis)
          .def("GetNormalOrderedHNAT", &HFMBPT::GetNormalOrderedHNAT)
          .def("PrintSPEandWF", &HFMBPT::PrintSPEandWF)
          .def_readwrite("C_HO2NAT", &HFMBPT::C_HO2NAT) // Unitary transformation
          .def_readwrite("C_HF2NAT", &HFMBPT::C_HF2NAT) // Unitary transformation
          ;

      // Define which overloaded version of IMSRGSolver::Transform I want to expose
      //   Operator (IMSRGSolver::*Transform_ref)(Operator&) = &IMSRGSolver::Transform;

      py::class_<IMSRGSolver>(m, "IMSRGSolver")
          .def(py::init<Operator &>())
          .def("Solve", &IMSRGSolver::Solve)
          //      .def("Transform",Transform_ref)
          .def("Transform", [](IMSRGSolver &self, Operator &op)
               { return self.Transform(op); })
          .def("InverseTransform", &IMSRGSolver::InverseTransform)
          .def("SetFlowFile", &IMSRGSolver::SetFlowFile)
          .def("SetMethod", &IMSRGSolver::SetMethod)
          .def("SetEtaCriterion", &IMSRGSolver::SetEtaCriterion)
          .def("SetDs", &IMSRGSolver::SetDs)
          .def("SetdOmega", &IMSRGSolver::SetdOmega)
          .def("SetOmegaNormMax", &IMSRGSolver::SetOmegaNormMax)
          .def("SetSmax", &IMSRGSolver::SetSmax)
          .def("SetDsmax", &IMSRGSolver::SetDsmax)
          .def("SetHin", &IMSRGSolver::SetHin)
          .def("SetODETolerance", &IMSRGSolver::SetODETolerance)
          .def("Reset", &IMSRGSolver::Reset)
          .def("SetGenerator", &IMSRGSolver::SetGenerator)
          .def("SetOnly2bEta", [](IMSRGSolver &self, bool tf)
               { self.GetGenerator().SetOnly2bEta(tf); })
          .def("SetDenominatorCutoff", &IMSRGSolver::SetDenominatorCutoff)
          .def("SetDenominatorDelta", &IMSRGSolver::SetDenominatorDelta)
          .def("SetDenominatorDeltaOrbit", &IMSRGSolver::SetDenominatorDeltaOrbit)
          .def("SetDenominatorPartitioning", &IMSRGSolver::SetDenominatorPartitioning)
          .def("GetSystemDimension", &IMSRGSolver::GetSystemDimension)
          .def("GetOmega", &IMSRGSolver::GetOmega)
          .def("SetOmega", &IMSRGSolver::SetOmega)
          //      .def("GetH_s",&IMSRGSolver::GetH_s,return_value_policy<reference_existing_object>())
          .def("GetH_s", &IMSRGSolver::GetH_s)
          .def("SetH_s", &IMSRGSolver::SetH_s)
          .def("GetS", &IMSRGSolver::GetS)
          .def("SetMagnusAdaptive", &IMSRGSolver::SetMagnusAdaptive)
          .def("SetReadWrite", &IMSRGSolver::SetReadWrite)
          .def("SetHunterGatherer", &IMSRGSolver::SetHunterGatherer)
//          .def("SetPerturbativeTriples", &IMSRGSolver::SetPerturbativeTriples)
//          .def("GetPerturbativeTriples", &IMSRGSolver::GetPerturbativeTriples)
          .def("CalculatePerturbativeTriples", &IMSRGSolver::CalculatePerturbativeTriples)
          .def("AddOperator", &IMSRGSolver::AddOperator)
          .def("GetOperator", &IMSRGSolver::GetOperator)
          .def("EstimateBCHError", &IMSRGSolver::EstimateBCHError)
          .def("UpdateEta", &IMSRGSolver::UpdateEta)
          .def("GetNOmegaWritten", &IMSRGSolver::GetNOmegaWritten)
          .def("GetOmegaSize", &IMSRGSolver::GetOmegaSize)
          //      .def("GetScratchDir",[](IMSRGSolver& self){ return self.rw->GetScratchDir();} )
          .def("GetScratchDir", [](IMSRGSolver &self)
               { return self.scratchdir; })
          .def("FlushOmegaToScratch", &IMSRGSolver::FlushOmegaToScratch)
          .def_readwrite("generator", &IMSRGSolver::generator)
          .def_readwrite("Eta", &IMSRGSolver::Eta)
          .def_readwrite("n_omega_written", &IMSRGSolver::n_omega_written) // I'm not sure I like just directly exposing this...
          ;

      py::class_<Generator>(m, "Generator")
          .def(py::init<>())
          .def("SetType", &Generator::SetType, py::arg("gen_type"))
          .def("SetDenominatorPartitioning", &Generator::SetDenominatorPartitioning, py::arg("Moller_Plessett or Epstein_Nesbet"))
          .def("Update", &Generator::Update, py::arg("H"), py::arg("Eta"))
          .def("GetHod_SingleRef", &Generator::GetHod_SingleRef, py::arg("H"))
          .def("GetHod", &Generator::GetHod, py::arg("H"));

      py::class_<IMSRGProfiler>(m, "IMSRGProfiler")
          .def(py::init<>())
          .def("PrintTimes", &IMSRGProfiler::PrintTimes)
          .def("PrintCounters", &IMSRGProfiler::PrintCounters)
          .def("PrintAll", &IMSRGProfiler::PrintAll)
          .def("PrintMemory", &IMSRGProfiler::PrintMemory)
          .def("Clear", &IMSRGProfiler::Clear);

      py::class_<Jacobi3BME>(m, "Jacobi3BME")
          .def(py::init<>())
          .def(py::init<int, int, int, int, int>())
          .def("GetDimensionAS", &Jacobi3BME::GetDimensionAS)
          .def("GetDimensionNAS", &Jacobi3BME::GetDimensionNAS)
          .def("GetMatElAS", &Jacobi3BME::GetMatElAS)
          .def("GetMatElNAS", &Jacobi3BME::GetMatElNAS)
          .def("SetEmax", &Jacobi3BME::SetEmax)
          .def("SetE2max", &Jacobi3BME::SetE2max)
          .def("SetE3max", &Jacobi3BME::SetE3max)
          .def("ComputeNAS_MatrixElements", &Jacobi3BME::ComputeNAS_MatrixElements)
          .def("GetLabMatEl", &Jacobi3BME::GetLabMatEl)
          .def("TestReadTcoeffNavratil", &Jacobi3BME::TestReadTcoeffNavratil)
          .def("GetV3mon_all", &Jacobi3BME::GetV3mon_all);

      py::module Commutator = m.def_submodule("Commutator", "Commutator namespace");
      Commutator.def("Set_BCH_Transform_Threshold", &Commutator::Set_BCH_Transform_Threshold);
      Commutator.def("Set_BCH_Product_Threshold", &Commutator::Set_BCH_Product_Threshold);
      Commutator.def("Commutator", &Commutator::Commutator);
      Commutator.def("BCH_Transform", &Commutator::BCH_Transform);
      Commutator.def("BCH_Product", &Commutator::BCH_Product);
      Commutator.def("EstimateBCHError", &Commutator::EstimateBCHError);
      Commutator.def("SetUseIMSRG3", &Commutator::SetUseIMSRG3);
      Commutator.def("SetUseIMSRG3N7", &Commutator::SetUseIMSRG3N7);
      Commutator.def("TurnOnTerm", &Commutator::TurnOnTerm);
      Commutator.def("TurnOffTerm", &Commutator::TurnOffTerm);
      Commutator.def("SetThreebodyThreshold", &Commutator::SetThreebodyThreshold);
      Commutator.def("SetIMSRG3Verbose", &Commutator::SetIMSRG3Verbose, py::arg("tf"));
      Commutator.def("SetSingleThread", &Commutator::SetSingleThread, py::arg("tf"));
      Commutator.def("SetOnly2bOmega", &Commutator::SetOnly2bOmega);

      Commutator.def("SetUseFactorizedCorrection", &Commutator::SetUseFactorizedCorrection);
      Commutator.def("SetUseFactorized_GooseTank_Correction_1b", &Commutator::SetUseFactorized_GooseTank_Correction_1b);
      Commutator.def("SetUseFactorized_GooseTank_Correction_2b", &Commutator::SetUseFactorized_GooseTank_Correction_2b);
      Commutator.def("SetUseFactorized_GooseTank_Correction_only_1b", &Commutator::SetUseFactorized_GooseTank_Correction_only_1b);
      Commutator.def("SetUseFactorized_GooseTank_Correction_only_2b", &Commutator::SetUseFactorized_GooseTank_Correction_only_2b);
      Commutator.def("SetUseFactorized_Correct_ZBTerm", &Commutator::SetUseFactorized_Correct_ZBTerm);
      Commutator.def("SetUseFactorized_TypeII_Correction_1b", &Commutator::SetUseFactorized_TypeII_Correction_1b);
      Commutator.def("SetUseFactorized_TypeIII_Correction_1b", &Commutator::SetUseFactorized_TypeIII_Correction_1b);
      Commutator.def("SetUseFactorized_TypeII_Correction_2b", &Commutator::SetUseFactorized_TypeII_Correction_2b);
      Commutator.def("SetUseFactorized_TypeIII_Correction_2b", &Commutator::SetUseFactorized_TypeIII_Correction_2b);
      Commutator.def("UseSlowVersionOfDoubleCommutator", &Commutator::UseSlowVersionOfDoubleCommutator);

      Commutator.def("SetUseFactorized_GT_TypeI_Correction_2b", &Commutator::SetUseFactorized_GT_TypeI_Correction_2b);
      Commutator.def("SetUseFactorized_GT_TypeIV_Correction_2b", &Commutator::SetUseFactorized_GT_TypeIV_Correction_2b);

      // IMSRG(2) commutators
      Commutator.def("comm110ss", &Commutator::comm110ss);
      Commutator.def("comm220ss", &Commutator::comm220ss);
      Commutator.def("comm111ss", &Commutator::comm111ss);
      Commutator.def("comm121ss", &Commutator::comm121ss);
      Commutator.def("comm221ss", &Commutator::comm221ss);
      Commutator.def("comm122ss", &Commutator::comm122ss);
      Commutator.def("comm222_pp_hh_221ss", &Commutator::comm222_pp_hh_221ss);
      Commutator.def("comm222_pp_hhss", &Commutator::comm222_pp_hhss);
      Commutator.def("comm222_phss", &Commutator::comm222_phss);
      // IMSRG(3) commutators
      Commutator.def("comm330ss", &Commutator::comm330ss);
      Commutator.def("comm331ss", &Commutator::comm331ss);
      Commutator.def("comm231ss", &Commutator::comm231ss);
      Commutator.def("comm132ss", &Commutator::comm132ss);
      Commutator.def("comm232ss", &Commutator::comm232ss);
      Commutator.def("comm332_ppph_hhhpss", &Commutator::comm332_ppph_hhhpss);
      Commutator.def("comm332_pphhss", &Commutator::comm332_pphhss);
      Commutator.def("comm223ss", &Commutator::comm223ss);
      Commutator.def("comm133ss", &Commutator::comm133ss);
      Commutator.def("comm233_pp_hhss", &Commutator::comm233_pp_hhss);
      Commutator.def("comm233_phss", &Commutator::comm233_phss);
      Commutator.def("comm333_ppp_hhhss", &Commutator::comm333_ppp_hhhss);
      Commutator.def("comm333_pph_hhpss", &Commutator::comm333_pph_hhpss);
      // scalar-tensor commutators
      Commutator.def("comm111st", &Commutator::comm111st);
      Commutator.def("comm121st", &Commutator::comm121st);
      Commutator.def("comm122st", &Commutator::comm122st);
      Commutator.def("comm222_pp_hh_221st", &Commutator::comm222_pp_hh_221st);
      Commutator.def("comm222_phst", &Commutator::comm222_phst);
      Commutator.def("SetIMSRG3Noqqq", &Commutator::SetIMSRG3Noqqq);
      Commutator.def("SetIMSRG3valence2b", &Commutator::SetIMSRG3valence2b);
      Commutator.def("Discard0bFrom3b", &Commutator::Discard0bFrom3b);
      Commutator.def("Discard1bFrom3b", &Commutator::Discard1bFrom3b);
      Commutator.def("Discard2bFrom3b", &Commutator::Discard2bFrom3b);

      Commutator.def("comm223_231_Factorization", &Commutator::comm223_231_Factorization);
      Commutator.def("comm223_232_Factorization", &Commutator::comm223_232_Factorization);

      Commutator.def("comm223_231_Factorization_slow", &Commutator::comm223_231_Factorization_slow);
      Commutator.def("comm223_232_Factorization_slow", &Commutator::comm223_232_Factorization_slow);
      Commutator.def("comm223_232_Test", &Commutator::comm223_232_Test);
      Commutator.def("comm223_232_Test_2", &Commutator::comm223_232_Test_2);


      py::module ReferenceImplementations = m.def_submodule("ReferenceImplementations", "ReferenceImplementations namespace");
      ReferenceImplementations.def("comm110ss", &ReferenceImplementations::comm110ss);
      ReferenceImplementations.def("comm220ss", &ReferenceImplementations::comm220ss);
      ReferenceImplementations.def("comm111ss", &ReferenceImplementations::comm111ss);
      ReferenceImplementations.def("comm121ss", &ReferenceImplementations::comm121ss);
      ReferenceImplementations.def("comm221ss", &ReferenceImplementations::comm221ss);
      ReferenceImplementations.def("comm122ss", &ReferenceImplementations::comm122ss);
      ReferenceImplementations.def("comm222_pp_hh_221ss", &ReferenceImplementations::comm222_pp_hh_221ss);
      ReferenceImplementations.def("comm222_phss", &ReferenceImplementations::comm222_phss);
      //
      ReferenceImplementations.def("comm223ss", &ReferenceImplementations::comm223ss);
      ReferenceImplementations.def("comm232ss", &ReferenceImplementations::comm232ss);
      ReferenceImplementations.def("comm231ss", &ReferenceImplementations::comm231ss);
      //
      ReferenceImplementations.def("diagram_CIa", &ReferenceImplementations::diagram_CIa);
      ReferenceImplementations.def("diagram_CIb", &ReferenceImplementations::diagram_CIb);
      ReferenceImplementations.def("diagram_CIIa", &ReferenceImplementations::diagram_CIIa);
      ReferenceImplementations.def("diagram_CIIb", &ReferenceImplementations::diagram_CIIb);
      ReferenceImplementations.def("diagram_CIIc", &ReferenceImplementations::diagram_CIIc);
      ReferenceImplementations.def("diagram_CIId", &ReferenceImplementations::diagram_CIId);
      ReferenceImplementations.def("diagram_CIIIa", &ReferenceImplementations::diagram_CIIIa);
      ReferenceImplementations.def("diagram_CIIIb", &ReferenceImplementations::diagram_CIIIb);
      ReferenceImplementations.def("diagram_DIa", &ReferenceImplementations::diagram_DIa);
      ReferenceImplementations.def("diagram_DIb", &ReferenceImplementations::diagram_DIb);
      ReferenceImplementations.def("diagram_DIVa", &ReferenceImplementations::diagram_DIVa);
      ReferenceImplementations.def("diagram_DIVb", &ReferenceImplementations::diagram_DIVb);
      ReferenceImplementations.def("diagram_DIVb_intermediate", &ReferenceImplementations::diagram_DIVb_intermediate);
      ReferenceImplementations.def("comm223_231_BruteForce", &ReferenceImplementations::comm223_231_BruteForce);
      ReferenceImplementations.def("comm223_232_BruteForce", &ReferenceImplementations::comm223_232_BruteForce);




      py::class_<RPA>(m, "RPA")
          .def(py::init<Operator &>())
          .def("ConstructAMatrix", &RPA::ConstructAMatrix, py::arg("J"), py::arg("parity"), py::arg("Tz"), py::arg("Isovector"))
          .def("ConstructBMatrix", &RPA::ConstructBMatrix, py::arg("J"), py::arg("parity"), py::arg("Tz"), py::arg("Isovector"))
          .def("SolveCP", &RPA::SolveCP)
          .def("SolveTDA", &RPA::SolveTDA)
          .def("SolveRPA", &RPA::SolveRPA)
          .def("TransitionToGroundState", &RPA::TransitionToGroundState, py::arg("OpIn"), py::arg("mu"))
          .def("PVCouplingEffectiveCharge", &RPA::PVCouplingEffectiveCharge, py::arg("OpIn"), py::arg("k"), py::arg("l"))
          .def("GetEnergies", [](RPA &self)
               {arma::vec vals = self.GetEnergies(); std::vector<double> vvec; for (auto & v : vals) {vvec.push_back(v);};  return vvec; })
          .def("GetX", [](RPA &self, size_t i)
               {arma::vec vals = self.GetX(i); std::vector<double> vvec; for (auto & v : vals) {vvec.push_back(v);};  return vvec; })
          .def("GetY", [](RPA &self, size_t i)
               {arma::vec vals = self.GetY(i); std::vector<double> vvec; for (auto & v : vals) {vvec.push_back(v);};  return vvec; })
          .def("PrintA", [](RPA &self)
               { std::cout << self.A << std::endl; })
          .def("PrintB", [](RPA &self)
               { std::cout << self.B << std::endl; })
          .def("GetEgs", &RPA::GetEgs);

      py::class_<UnitTest>(m, "UnitTest")
          //      .def(py::init<>())
          .def(py::init<ModelSpace &>())
          .def("SetRandomSeed", &UnitTest::SetRandomSeed)
          .def("RandomOp", &UnitTest::RandomOp)
          .def("TestCommutators", &UnitTest::TestCommutators)
          .def("TestCommutators3", &UnitTest::TestCommutators3)
          .def("TestDaggerCommutators", &UnitTest::TestDaggerCommutators)
          .def("TestDaggerCommutatorsAlln", &UnitTest::TestDaggerCommutatorsAlln)
          .def("Test3BodyAntisymmetry", &UnitTest::Test3BodyAntisymmetry)
          .def("Test3BodyHermiticity", &UnitTest::Test3BodyHermiticity)
          .def("TestRPAEffectiveCharge", &UnitTest::TestRPAEffectiveCharge, py::arg("H"), py::arg("OpIn"), py::arg("k"), py::arg("l"))
          .def("SanityCheck", &UnitTest::SanityCheck)
          .def("Test_comm110ss", &UnitTest::Test_comm110ss)
          .def("Test_comm220ss", &UnitTest::Test_comm220ss)
          .def("Test_comm111ss", &UnitTest::Test_comm111ss)
          .def("Test_comm121ss", &UnitTest::Test_comm121ss)
          .def("Test_comm221ss", &UnitTest::Test_comm221ss)
          .def("Test_comm122ss", &UnitTest::Test_comm122ss)
          .def("Test_comm222_pp_hhss", &UnitTest::Test_comm222_pp_hhss)
          .def("Test_comm222_phss", &UnitTest::Test_comm222_phss)
          .def("Test_comm222_pp_hh_221ss", &UnitTest::Test_comm222_pp_hh_221ss)

          .def("Test_comm222_phst", &UnitTest::Test_comm222_phst)
          ///
          .def("Test_comm330ss", &UnitTest::Test_comm330ss)
          .def("Test_comm331ss", &UnitTest::Test_comm331ss)
          .def("Test_comm231ss", &UnitTest::Test_comm231ss)
          .def("Test_comm132ss", &UnitTest::Test_comm132ss)
          .def("Test_comm232ss", &UnitTest::Test_comm232ss)
          .def("Test_comm223ss", &UnitTest::Test_comm223ss)
          .def("Test_comm133ss", &UnitTest::Test_comm133ss)
          .def("Test_comm332_ppph_hhhpss", &UnitTest::Test_comm332_ppph_hhhpss)
          .def("Test_comm332_pphhss", &UnitTest::Test_comm332_pphhss)
          .def("Test_comm233_pp_hhss", &UnitTest::Test_comm233_pp_hhss)
          .def("Test_comm233_phss", &UnitTest::Test_comm233_phss)
          .def("Test_comm333_ppp_hhhss", &UnitTest::Test_comm333_ppp_hhhss)
          .def("Test_comm333_pph_hhpss", &UnitTest::Test_comm333_pph_hhpss)

          .def("Mscheme_Test_comm110ss", &UnitTest::Mscheme_Test_comm110ss)
          .def("Mscheme_Test_comm220ss", &UnitTest::Mscheme_Test_comm220ss)
          .def("Mscheme_Test_comm111ss", &UnitTest::Mscheme_Test_comm111ss)
          .def("Mscheme_Test_comm121ss", &UnitTest::Mscheme_Test_comm121ss)
          .def("Mscheme_Test_comm221ss", &UnitTest::Mscheme_Test_comm221ss)
          .def("Mscheme_Test_comm122ss", &UnitTest::Mscheme_Test_comm122ss)
          .def("Mscheme_Test_comm222_pp_hhss", &UnitTest::Mscheme_Test_comm222_pp_hhss)
          .def("Mscheme_Test_comm222_phss", &UnitTest::Mscheme_Test_comm222_phss)
          //
          //      .def("Mscheme_Test_comm222_pp_hh_221ss", &UnitTest::Mscheme_Test_comm222_pp_hh_221ss)
          ///
          .def("Mscheme_Test_comm330ss", &UnitTest::Mscheme_Test_comm330ss)
          .def("Mscheme_Test_comm331ss", &UnitTest::Mscheme_Test_comm331ss)
          .def("Mscheme_Test_comm231ss", &UnitTest::Mscheme_Test_comm231ss)
          .def("Mscheme_Test_comm132ss", &UnitTest::Mscheme_Test_comm132ss)
          .def("Mscheme_Test_comm232ss", &UnitTest::Mscheme_Test_comm232ss)
          .def("Mscheme_Test_comm223ss", &UnitTest::Mscheme_Test_comm223ss)
          .def("Mscheme_Test_comm133ss", &UnitTest::Mscheme_Test_comm133ss)
          .def("Mscheme_Test_comm332_ppph_hhhpss", &UnitTest::Mscheme_Test_comm332_ppph_hhhpss)
          .def("Mscheme_Test_comm332_pphhss", &UnitTest::Mscheme_Test_comm332_pphhss)
          .def("Mscheme_Test_comm233_pp_hhss", &UnitTest::Mscheme_Test_comm233_pp_hhss)
          .def("Mscheme_Test_comm233_phss", &UnitTest::Mscheme_Test_comm233_phss)
          .def("Mscheme_Test_comm333_ppp_hhhss", &UnitTest::Mscheme_Test_comm333_ppp_hhhss)
          .def("Mscheme_Test_comm333_pph_hhpss", &UnitTest::Mscheme_Test_comm333_pph_hhpss)
          //      .def("Test3BodySetGet",&UnitTest::Test3BodySetGet)
          .def("GetMschemeMatrixElement_1b", &UnitTest::GetMschemeMatrixElement_1b, py::arg("Op"), py::arg("a"), py::arg("ma"), py::arg("b"), py::arg("mb")) // Op, a,ma, b,mb...
          .def("GetMschemeMatrixElement_2b", &UnitTest::GetMschemeMatrixElement_2b)                                                                          // Op, a,ma, b,mb...
          .def("GetMschemeMatrixElement_3b", &UnitTest::GetMschemeMatrixElement_3b)                                                                          // Op, a,ma, b,mb...

          ;

      //  py::class_<SymmMatrix<double>>(m,"SymmMatrix")
      //     .def(py::init<size_t>())
      //     .def(py::init<size_t,int>())
      //     .def("Get",&SymmMatrix<double>::Get)
      //     .def("Put",&SymmMatrix<double>::Put)
      //     .def("FullMatrix",&SymmMatrix<double>::FullMatrix)
      //  ;

      m.def("BuildVersion", version::BuildVersion);

      m.def("TCM_Op", imsrg_util::TCM_Op);
      m.def("Trel_Op", imsrg_util::Trel_Op);
      m.def("R2CM_Op", imsrg_util::R2CM_Op);
      m.def("HCM_Op", imsrg_util::HCM_Op);
      m.def("NumberOp", imsrg_util::NumberOp);
      m.def("RSquaredOp", imsrg_util::RSquaredOp);
      m.def("RpSpinOrbitCorrection", imsrg_util::RpSpinOrbitCorrection);
      m.def("E0Op", imsrg_util::E0Op);
      m.def("AllowedFermi_Op", imsrg_util::AllowedFermi_Op);
      m.def("AllowedGamowTeller_Op", imsrg_util::AllowedGamowTeller_Op);
      m.def("ElectricMultipoleOp", imsrg_util::ElectricMultipoleOp);
      m.def("MagneticMultipoleOp", imsrg_util::MagneticMultipoleOp);
      m.def("Sigma_Op", imsrg_util::Sigma_Op);
      m.def("Isospin2_Op", imsrg_util::Isospin2_Op);
      m.def("LdotS_Op", imsrg_util::LdotS_Op);
      m.def("HO_density", imsrg_util::HO_density);
      m.def("GetOccupationsHF", imsrg_util::GetOccupationsHF);
      m.def("GetOccupations", imsrg_util::GetOccupations);
      m.def("GetDensity", imsrg_util::GetDensity);
      m.def("CommutatorTest", imsrg_util::CommutatorTest);
      m.def("Calculate_p1p2_all", imsrg_util::Calculate_p1p2_all);
      m.def("Single_Ref_1B_Density_Matrix", imsrg_util::Single_Ref_1B_Density_Matrix);
      m.def("Get_Charge_Density", imsrg_util::Get_Charge_Density);
      m.def("Embed1BodyIn2Body", imsrg_util::Embed1BodyIn2Body);
      m.def("RadialIntegral", imsrg_util::RadialIntegral);
      m.def("RadialIntegral_RpowK", imsrg_util::RadialIntegral_RpowK);
      m.def("RadialIntegral_Gauss", imsrg_util::RadialIntegral_Gauss, py::arg("na"), py::arg("la"), py::arg("nb"), py::arg("lb"), py::arg("sig"));
      m.def("RPA_resummed_1b", imsrg_util::RPA_resummed_1b, py::arg("OpIn"), py::arg("H"), py::arg("mode"));
      m.def("FirstOrderCorr_1b", imsrg_util::FirstOrderCorr_1b, py::arg("OpIn"), py::arg("H"));
      m.def("FrequencyConversionCoeff", imsrg_util::FrequencyConversionCoeff);
      m.def("OperatorFromString", imsrg_util::OperatorFromString);
      m.def("HO_Radial_psi", imsrg_util::HO_Radial_psi, py::arg("n"), py::arg("l"), py::arg("hw"), py::arg("r"));
      m.def("MBPT2_SpectroscopicFactor", imsrg_util::MBPT2_SpectroscopicFactor);
      m.def("SerberTypePotential", imsrg_util::SerberTypePotential, py::arg("modelspace"), py::arg("V0"), py::arg("mu"), py::arg("A"), py::arg("B"), py::arg("C"), py::arg("D"));

      m.def("CG", AngMom::CG);
      m.def("ThreeJ", AngMom::ThreeJ);
      m.def("SixJ", AngMom::SixJ);
      m.def("NineJ", AngMom::NineJ);
      m.def("NormNineJ", AngMom::NormNineJ);
      m.def("Moshinsky", AngMom::Moshinsky, py::arg("N"), py::arg("LAM"), py::arg("n"), py::arg("lam"), py::arg("n1"), py::arg("l1"), py::arg("n2"), py::arg("l2"), py::arg("L"), py::arg("BETA"));
      m.def("TalmiB", AngMom::TalmiB);
      m.def("TalmiI", imsrg_util::TalmiI);
      m.def("Tcoeff", AngMom::Tcoeff);
      m.def("SetUseGooseTank", Commutator::SetUseGooseTank);
      m.def("SetUseIMSRG3", Commutator::SetUseIMSRG3);
      m.def("SetUseIMSRG3N7", Commutator::SetUseIMSRG3N7);
      m.def("FillFactorialLists", AngMom::FillFactorialLists);
      m.def("factorial", AngMom::factorial);
      m.def("double_fact", AngMom::double_fact);

      m.attr("HBARC") = py::float_(PhysConst::HBARC);
      m.attr("M_PROTON") = py::float_(PhysConst::M_PROTON);
      m.attr("M_NEUTRON") = py::float_(PhysConst::M_NEUTRON);
      m.attr("M_NUCLEON") = py::float_(PhysConst::M_NUCLEON);
      m.attr("M_ELECTRON") = py::float_(PhysConst::M_ELECTRON);
      m.attr("M_PION_CHARGED") = py::float_(PhysConst::M_PION_CHARGED);
      m.attr("M_PION_NEUTRAL") = py::float_(PhysConst::M_PION_NEUTRAL);
      m.attr("NUCLEON_VECTOR_G") = py::float_(PhysConst::NUCLEON_VECTOR_G);
      m.attr("NUCLEON_AXIAL_G") = py::float_(PhysConst::NUCLEON_AXIAL_G);
      m.attr("PROTON_SPIN_G") = py::float_(PhysConst::PROTON_SPIN_G);
      m.attr("NEUTRON_SPIN_G") = py::float_(PhysConst::NEUTRON_SPIN_G);
      m.attr("ELECTRON_SPIN_G") = py::float_(PhysConst::ELECTRON_SPIN_G);
      m.attr("ALPHA_FS") = py::float_(PhysConst::ALPHA_FS);
      m.attr("F_PI") = py::float_(PhysConst::F_PI);
      m.attr("HARTREE") = py::float_(PhysConst::HARTREE);
      m.attr("PROTON_RCH2") = py::float_(PhysConst::PROTON_RCH2);
      m.attr("NEUTRON_RCH2") = py::float_(PhysConst::NEUTRON_RCH2);
      m.attr("DARWIN_FOLDY") = py::float_(PhysConst::DARWIN_FOLDY);
}
