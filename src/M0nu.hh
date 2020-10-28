#ifndef M0nu_hh
#define M0nu_hh

#include "ModelSpace.hh"
#include "Operator.hh"
#include "imsrg_util.hh"
#include "omp.h"
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
  const double CUTV = 850.0; ///< vector cut off
  const double CUTA = 1086.0; ///< axial cut off
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
  double GTFormFactor(double q); ///< Gamow Teller Form Factor
  double FermiFormFactor(double q); ///< Fermi Form Factor
  double TensorFormFactor(double q); ///< Tensor Form Factor

  double A(double p, double pp, int J, double Eclosure, std::string transition, gsl_integration_glfixed_table * t, int norm, int size);
  uint64_t AHash(int i, int j, int J, int norm);
  void AUnHash(uint64_t key, uint64_t& i, uint64_t& j, uint64_t& J, uint64_t& norm);
  std::unordered_map<uint64_t,double> PrecalculateA(int e2max,double Eclosure, std::string transition, int npoints);
  double GetA(int i, int j, int J, int norm,int l, std::unordered_map<uint64_t,double> &AList);
  double W_fermi_gt(double p, double pp,int index_p,int index_pp, int l, int lp, int J, std::unordered_map<uint64_t,double>& AList);
  double W_tensor(double p, double pp, int index_p, int index_pp, int l, int lp, int J, std::unordered_map<uint64_t,double>& AList);
  
  double fq(double p, double pp,int i, int j, int n, int l, int np, int lp, int J, double hw, std::string transition, double Eclosure, std::string src, std::unordered_map<uint64_t,double>& AList);
  double integrate_dq(int n, int l, int np, int lp, int J, double hw, std::string transition, double Eclosure, std::string src,int npoints, gsl_integration_glfixed_table * t, std::unordered_map<uint64_t,double>& AList);

  uint64_t IntHash(int n, int l, int np, int lp, int J);
  void IntUnHash(uint64_t key, uint64_t& n, uint64_t& l, uint64_t& np, uint64_t& lp, uint64_t& J);
  std::unordered_map<uint64_t,double> PreCalculateM0nuIntegrals(int e2max, double hw, std::string transition, double Eclosure, std::string src);
  double GetM0nuIntegral(int e2max, int n, int l, int np, int lp,int J, double hw, std::string transition,  double Eclosure, std::string src, std::unordered_map<uint64_t,double> &IntList);
  
  Operator GamowTeller(ModelSpace& modelspace, double Eclosure, std::string src); ///< Gamow-Teller part of M0nu operator
  Operator Fermi(ModelSpace& modelspace, double Eclosure, std::string src); ///< Fermi part of M0nu operator
  Operator Tensor(ModelSpace& modelspace, double Eclosure, std::string src); ///< Tensor part of M0nu operator
  Operator DGT_Op(ModelSpace& modelspace); //< Double Gamow-Teller operator

}


#endif
