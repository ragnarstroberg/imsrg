#ifndef M0nu_hh
#define M0nu_hh

#include "ModelSpace.hh"
#include "Operator.hh"
//#include "imsrg_util.hh"
#include "omp.h"
#include "Pwd.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "PhysicalConstants.hh"
#include <iostream>

//Code to compute the operators for neutrinoless double beta decay
//Mainly implementation of operators found in Charlie Payne master
//thesis found at "file:///Users/antoinebelley/Downloads/ubc_2018_february_payne_charlie.pdf'
//Implemented by Antoine Belley
//Date : 05/2020
namespace M0nu
{
  const double MAGMOM = 4.706;  ///< Magnetic moment
  const double CUTV = 850.0; ///< vector cut off in MeV
  const double CUTA = 1086.0; ///< axial cut off in MeV
  const double R0 = 1.2; ///< Radius prefactor for R = R0A^1/3


  inline double asNorm(int i, int j) {return i==j ? PhysConst::INVSQRT2 : 1.0;}; //for anti-symetrization
  // inline double asNorm(int i, int j) {return i==j ? 0.5 : PhysConst::INVSQRT2;};

  int decimalgen(int a, int b, int c, int d, int maxb, int maxc, int maxd);
  int phase(int x);  
  double HO_Radial_psi_mom(int n, int l, double hw, double p);


  double gv_func(double qsq);
  double ga_func(double qsq);
  double gm_func(double qsq);
  double gp_func(double qsq);
  double hF_VV(double qsq);
  double hGT_AA(double qsq);
  double hGT_AP(double qsq);
  double hGT_PP(double qsq);
  double hGT_MM(double qsq);
  double hT_AA(double qsq);
  double hT_AP(double qsq);
  double hT_PP(double qsq);
  double hT_MM(double qsq);
  double GTFormFactor(double qsq);
  double FermiFormFactor(double qsq); 
  double TensorFormFactor(double qsq);

  double integrate_dq(int n, int l, int np, int lp, int S, int J, double hw, PWD &pwd);
  uint64_t IntHash(int n, int l, int np, int lp, int S, int J);
  void IntUnHash(uint64_t key, uint64_t &n, uint64_t &l, uint64_t &np, uint64_t &lp,uint64_t &S, uint64_t &J);
  std::unordered_map<uint64_t, double> PreCalculateMomentumSpaceIntegrals(int e2max, int Smin, int Smax, double hw, int l_diff, PWD &pwd);
  double GetM0nuIntegral(int e2max, int n, int l, int np, int lp, int J, int S, std::unordered_map<uint64_t, double> &IntList);

  void TalmiMoshinkyTransform(ModelSpace &modelspace, double &sumTMT, double &sumTMTas, int na, int la, int nb, int lb, int nc, int lc, int nd, int ld,
                                                 int Li, int Lf, int S, int J, int rank, int e2max, std::unordered_map<uint64_t, double>& RelativeFrameOPList);

  Operator GamowTeller(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor); ///< Gamow-Teller part of M0nu operator
  Operator Fermi(ModelSpace& modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor); ///< Fermi part of M0nu operator
  Operator Tensor(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor); ///< Tensor part of M0nu operator
  Operator Contact(ModelSpace& modelspace, double regulator_cutoff, int regualtor_power); ///<Contact part of M0nu operator
  Operator DGT_Op(ModelSpace& modelspace); //< Double Gamow-Teller operator

  Operator GamowTellerHeavy(ModelSpace &modelspace,  std::string src, std::function<double(double)> formfactor); ///< Gamow-Teller part of M0nu operator
  Operator FermiHeavy(ModelSpace &modelspace, std::string src, std::function<double(double)> formfactor);       ///< Fermi part of M0nu operator
  Operator TensorHeavy(ModelSpace &modelspace, std::string src, std::function<double(double)> formfactor);      ///< Tensor part of M0nu operator

  double HO_Radial_psi(int n, int l, double hw, double r);
  double fq_radial_GT(double q, double Eclosure, double r12);
  double integrate_dq_radial_GT(double Eclosure, double r12,  int npoints, gsl_integration_glfixed_table * t);
  std::unordered_map<uint64_t,double> PreCalculateM0nuIntegrals_R(int e2max, double hw, double Eclosure, double r12);
  double GetM0nuIntegral_R(int e2max, int n, int l, int np, int lp,int J, double hw, double Eclosure, double r12, std::unordered_map<uint64_t,double> &IntList);
  Operator GamowTeller_R(ModelSpace& modelspace, double Eclosure, double r12);
  Operator DGT_R(ModelSpace& modelspace, double r12);

  long double TalmiB(int na, int la, int nb, int lb, int p);
  double RadialIntegral_Gauss( int na, int la, int nb, int lb, double sigma );
  double integrateRcom(int Ncom, int  Lambda, double  hw, double Rnucl);//,gsl_integration_glfixed_table * t);
  std::unordered_map<uint64_t,double> PreCalculateDGTRComIntegrals(int e2max, double hw, double Rnucl);
  double GetDGTRcomIntegral(int Ncom, int Lam, double hw, std::unordered_map<uint64_t,double> &IntList);
  Operator DGT_R_SurfaceLocalization(ModelSpace& modelspace, double r12);

  double dq_pp(int n, int l, int np, int lp, int S, int J, double hw, PWD &pwd, double p, double pp);
  std::unordered_map<uint64_t, double> PreCalculateMomentumSpaceIntegrands(int e2max, int Smin, int Smax, double hw, int Lrank, PWD &pwd, int index_p, int index_pp);
  // Operator TwoBody_Scalar_integrand(ModelSpace &modelspace, PWD &pwd, int Jrank, int Lrank, double p, double pp, int Smin = 0, int Smax = 1);
  Operator GamowTeller_integrand(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor, int index_p, int index_pp);
  Operator Fermi_integrand(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor, int index_p, int index_pp);
  Operator Tensor_integrand(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor, int index_p, int index_pp);
  Operator Contact_integrand(ModelSpace &modelspace, double regulator_cutoff, int regualtor_power, int index_p, int index_pp);
}
  


#endif
