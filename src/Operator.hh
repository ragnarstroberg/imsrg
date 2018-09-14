///////////////////////////////////////////////////////////////////////////////////
//    Operator.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include "TwoBodyME.hh"
#include "ThreeBodyME.hh"
#include "IMSRGProfiler.hh"
#include <armadillo>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>

//using namespace std;

/// The Operator class provides a generic operator up to three-body, scalar or tensor.
/// The class contains lots of methods and overloaded operators so that the resulting
/// code that uses the operators can look as close as possible to the math that is
/// written down.
class Operator
{
 public:
  //Fields
  ModelSpace * modelspace; ///< Pointer to the associated modelspace
  double ZeroBody; ///< The zero body piece of the operator.
  arma::mat OneBody; ///< The one body piece of the operator, stored in a single NxN armadillo matrix, where N is the number of single-particle orbits.
  TwoBodyME TwoBody; ///< The two body piece of the operator.
  ThreeBodyME ThreeBody; ///< The three body piece of the operator.

  int rank_J; ///< Spherical tensor rank of the operator
  int rank_T; ///< Isotensor rank of the operator
  int parity; ///< Parity of the operator, 0=even 1=odd
  int particle_rank; ///< Maximum particle rank. Should be 2 or 3.

  int E2max; ///< For two-body matrix elements, \f$ e_i + e_j \leq \f$ E2max
  int E3max; ///< For three-body matrix elements, \f$ e_i + e_j + e_k \leq \f$ E3max

  bool hermitian;
  bool antihermitian;
  int nChannels; ///< Number of two-body channels \f$ J,\pi,T_z \f$ associated with the model space


  std::map<std::array<int,3>,std::vector<index_t> > OneBodyChannels;
  IMSRGProfiler profiler;

  static double bch_transform_threshold;
  static double bch_product_threshold;
  static bool use_brueckner_bch;
  static bool use_goose_tank_correction;
  static bool use_goose_tank_correction_titus;



  //Constructors
  ~Operator();
  Operator(); ///< Default constructor
  Operator(ModelSpace&); ///< Construct a 2-body scalar operator
  Operator(ModelSpace&, int Jrank, int Trank, int Parity, int part_rank);
  Operator( const Operator& rhs); ///< Copy constructor
  Operator( Operator&&);

  //Overloaded operators
  Operator& operator=( const Operator& rhs);
  Operator& operator+=( const Operator& rhs);
  Operator operator+( const Operator& rhs) const;
  Operator& operator+=( const double& rhs);
  Operator operator+( const double& rhs) const;
  Operator& operator-=( const Operator& rhs);
  Operator operator-( const Operator& rhs) const;
  Operator operator-( ) const;
  Operator& operator-=( const double& rhs);
  Operator operator-( const double& rhs) const;
  Operator& operator*=( const double rhs);
  Operator operator*( const double rhs) const;
  Operator& operator/=( const double rhs);
  Operator operator/( const double rhs) const;

  Operator& operator=(Operator&& rhs);

  //Methods
  Operator& TempOp(size_t n); ///< Static scratch space for calculations

  // One body setter/getters
  double GetOneBody(int i,int j) {return OneBody(i,j);};
//  void SetOneBody(int i, int j, double val) { OneBody(i,j) = val;};
  void SetOneBody(int i, int j, double val) ;
  int GetTwoBodyDimension(int ch_bra, int ch_ket){ return TwoBody.GetMatrix(ch_bra, ch_ket).n_cols;};
//  double GetTwoBody(int ch_bra, int ch_ket, int i, int j){ return TwoBody.GetMatrix(ch_bra, ch_ket)(i,j);};
  double GetTwoBody(int ch_bra, int ch_ket, int i, int j);
  void SetTwoBody(int J1, int p1, int T1, int J2, int p2, int T2, int i, int j, int k, int l, double v);

  void SetE3max(int e){E3max = e;};
  int GetE3max(){return E3max;};

  // Other setter-getters
  ModelSpace * GetModelSpace();
  void SetModelSpace(ModelSpace &ms){modelspace = &ms;};

  void Erase(); ///< Set all matrix elements to zero.
  void EraseZeroBody(){ZeroBody = 0;}; ///< set zero-body term to zero
  void EraseOneBody(); ///< set all one-body terms to zero
  void EraseTwoBody(); ///< set all two-body terms to zero
  void EraseThreeBody(); ///< set all two-body terms to zero

  void SetHermitian() ;
  void SetAntiHermitian() ;
  void SetNonHermitian() ;
  bool IsHermitian()const {return hermitian;};
  bool IsAntiHermitian()const {return antihermitian;};
  bool IsNonHermitian()const {return not (hermitian or antihermitian);};
  int GetParticleRank()const {return particle_rank;};
  int GetJRank()const {return rank_J;};
  int GetTRank()const {return rank_T;};
  int GetParity()const {return parity;};
  void SetParticleRank(int pr) {particle_rank = pr;};
//  void ResetTensorTransformFirstPass(){tensor_transform_first_pass=true;};

  void MakeReduced();
  void MakeNotReduced();
  void ChangeNormalization(double coeff);
  void MakeNormalized();
  void MakeUnNormalized();

  void ScaleZeroBody(double x);
  void ScaleOneBody(double x);
  void ScaleTwoBody(double x);
  void Symmetrize(); ///< Copy the upper-half triangle to the lower-half triangle for each matrix
  void AntiSymmetrize(); ///< Copy the upper-half triangle to the lower-half triangle with a minus sign.
  void SetUpOneBodyChannels();
  size_t Size();

  void WriteBinary(std::ofstream& ofs);
  void ReadBinary(std::ifstream& ifs);


  // The actually interesting methods
  Operator DoNormalOrdering(); ///< Calls DoNormalOrdering2() or DoNormalOrdering3(), depending on the rank of the operator.
  Operator DoNormalOrdering2(); ///< Returns the normal ordered two-body operator
  Operator DoNormalOrdering3(); ///< Returns the normal ordered three-body operator
  Operator UndoNormalOrdering() const; ///< Returns the operator normal-ordered wrt the vacuum
  Operator Truncate(ModelSpace& ms_new); ///< Returns the operator trunacted to the new model space

  void SetToCommutator(const Operator& X, const Operator& Y);
  void CommutatorScalarScalar( const Operator& X, const Operator& Y) ;
  void CommutatorScalarTensor( const Operator& X, const Operator& Y) ;
  friend Operator Commutator(const Operator& X, const Operator& Y) ; 
//  friend Operator CommutatorScalarScalar( const Operator& X, const Operator& Y) ;
//  friend Operator CommutatorScalarTensor( const Operator& X, const Operator& Y) ;

  Operator BCH_Product(  Operator& )  ; 
  Operator BCH_Transform( const Operator& ) ; 
  Operator Standard_BCH_Transform( const Operator& ) ; 
  Operator Brueckner_BCH_Transform( const Operator& ) ; 

  void CalculateKineticEnergy(); // Deprecated
  void Eye(); ///< set to identity operator -- unused

  double GetMP2_Energy();
  double GetMP3_Energy();
  double MP1_Eval(Operator& );

  void PrintTimes(){profiler.PrintAll();};


  double Norm() const;
  double OneBodyNorm() const;
  double TwoBodyNorm() const;


  double Trace(int Atrace, int Ztrace) const;

  void ScaleFermiDirac(Operator& H, double T, double Efermi);


  void PrintOneBody() const {OneBody.print();};
  void PrintTwoBody(int ch) const {TwoBody.PrintMatrix(ch,ch);};


  static void Set_BCH_Transform_Threshold(double x){bch_transform_threshold=x;};
  static void Set_BCH_Product_Threshold(double x){bch_product_threshold=x;};
  static void SetUseBruecknerBCH(bool tf){use_brueckner_bch = tf;};
  static void SetUseGooseTank(bool tf){use_goose_tank_correction = tf;};

  std::deque<arma::mat> InitializePandya(size_t nch, std::string orientation);
//  void DoPandyaTransformation(std::deque<arma::mat>&, std::deque<arma::mat>&, std::string orientation) const ;
  void DoPandyaTransformation(std::deque<arma::mat>&, std::string orientation) const ;
  void DoPandyaTransformation_SingleChannel(arma::mat& X, int ch_cc, std::string orientation) const ;
  void AddInversePandyaTransformation(const std::deque<arma::mat>&);
  void AddInversePandyaTransformation_SingleChannel(arma::mat& Z, int ch_cc);


  void comm110ss( const Operator& X, const Operator& Y) ; 
  void comm220ss( const Operator& X, const Operator& Y) ;
  void comm111ss( const Operator& X, const Operator& Y) ;
  void comm121ss( const Operator& X, const Operator& Y) ;
  void comm221ss( const Operator& X, const Operator& Y) ;
  void comm122ss( const Operator& X, const Operator& Y) ;
  void comm222_pp_hhss( const Operator& X, const Operator& Y) ;
  void comm222_phss( const Operator& X, const Operator& Y) ;
  void comm222_pp_hh_221ss( const Operator& X, const Operator& Y) ;

//  void GooseTankUpdate( const Operator& Omega, Operator& Nested, Operator& chi);
  void GooseTankUpdate( const Operator& Omega, const Operator& Nested);
//  void goose_tank_ss( const Operator& X, const Operator& Y);

// scalar-tensor commutators

  void ConstructScalarMpp_Mhh(const Operator& X, const Operator& Y, TwoBodyME& Mpp, TwoBodyME& Mhh) const;
  void ConstructScalarMpp_Mhh_GooseTank(const Operator& X, const Operator& Y, TwoBodyME& Mpp, TwoBodyME& Mhh) const;
//  void DoTensorPandyaTransformation(std::map<std::array<int,2>,arma::mat>&, std::map<std::array<int,2>,arma::mat>&) const;
  void DoTensorPandyaTransformation(std::map<std::array<index_t,2>,arma::mat>&) const;
  void DoTensorPandyaTransformation_SingleChannel(arma::mat& X, int ch_bra_cc, int ch_ket_cc) const;
  void AddInverseTensorPandyaTransformation(const std::map<std::array<index_t,2>,arma::mat>&);
  void AddInverseTensorPandyaTransformation_SingleChannel(arma::mat& Zbar, int ch_bra_cc, int ch_ket_cc);

  void comm111st( const Operator& X, const Operator& Y) ;
  void comm121st( const Operator& X, const Operator& Y) ;
  void comm122st( const Operator& X, const Operator& Y) ;
  void comm222_pp_hh_221st( const Operator& X, const Operator& Y) ;
  void comm222_phst( const Operator& X, const Operator& Y) ;

};

/// Non member function, multiply by scalar from left side
Operator operator*(const double lhs, const Operator& rhs);
Operator operator*(const double lhs, const Operator&& rhs);



#endif

