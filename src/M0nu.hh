#ifndef M0nu_hh
#define M0nu_hh

#include "ModelSpace.hh"
#include "Operator.hh"
#include "imsrg_util.hh"
#include "omp.h"
#include <gsl/gsl_math.h>
#include "PhysicalConstants.hh"
#include <iostream>

//Code to compute the operators for neutrinoless double beta decay
//Mainly implementation of operators found in Charlie Payne master
//thesis found at "file:///Users/antoinebelley/Downloads/ubc_2018_february_payne_charlie.pdf'
//Implemented by Antoine Belley
//Date : 05/2020
namespace M0nu
{
  const double MAGMOM = 4.706;  //Magnetic moment
  const double CUTV = 850.0; //vector cut off
  const double CUTA = 1086.0; //axial cut off
  const double R0 = 1.2;

  inline int intPhase(int i) {return (i%2)==0 ? 1 : -1;}; // calculates the phase of an integer
  inline double asNorm(int i, int j) {return i==j ? PhysConst::INVSQRT2 : 1.0;}; //for anti-symetrization
  inline int PairFN(int a, int b) {return ((a + b)*(a + b + 1)/2) + b;}; // a Cantor pairing function used for cacheing CG in M0nuT

  struct frparams{int n; int l; int np; int lp; double hw; int rho; double q;};
  struct fqparams{int n; int l; int np; int lp; double hw; std::string transition; int rho; double Eclosure; std::string src;};
  
  //double HO_Radial_psi(int n, int l, double hw, double r);
  double frRBME(double r, void *params);
  void   GSLerror(std::string errstr, int status, size_t limit, double epsabs, double epsrel, double abserr, double result, int n, int l, int np, int lp, double q);
  double rbme(int n, int l, int np, int lp, double hw, int rho, double q, std::string src);
  double formfactor(double q, int transition);
  double fq(double q, void *params);
  double integrate_dq(int n, int l, int np, int lp, double hw, std::string transition, int rho, double Eclosure, std::string src);

  uint64_t IntHash(uint64_t n, uint64_t l, uint64_t np, uint64_t lp);
  void IntUnHash(uint64_t key, uint64_t &n, uint64_t &l, uint64_t &np, uint64_t &lp);
  std::unordered_map<uint64_t,double> PreCalculateM0nuIntegrals(int e2max, double hw, std::string transition, int rho, double Eclosure, std::string src);
  double GetM0nuIntegral(int e2max, int n, int l, int np, int lp, double hw, std::string transition, int rho, double Eclosure, std::string src, std::unordered_map<uint64_t,double>& Intlist);

  uint64_t T6jHash(uint64_t l1, uint64_t L1, uint64_t R, uint64_t L2, uint64_t l2);
  void T6jUnHash(uint64_t key, uint64_t &l1, uint64_t &L1, uint64_t &R, uint64_t &L2, uint64_t &l2);
  std::unordered_map<uint64_t,double> PreCalculateM0nuT6j(int e2max);
  double GetM0nuT6j(int l1, int L1, int R, int L2, int l2, std::unordered_map<uint64_t,double>& T6jList);

  Operator GamowTeller(ModelSpace& modelspace, double Eclosure, std::string src);
  Operator Fermi(ModelSpace& modelspace, double Eclosure, std::string src);
  Operator Tensor(ModelSpace& modelspace, double Eclosure, std::string src);

}


#endif
