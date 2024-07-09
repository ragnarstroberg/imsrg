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
//#include "ThreeBodyMENO2B.hh"
//#include "ThreeBodyMEpn.hh"
#include "ThreeLegME.hh"
#include "IMSRGProfiler.hh"
#include <armadillo>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <deque>
#include <set>
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
//  ThreeBodyMEpn ThreeBody; ///< The three body piece of the operator.
  ThreeLegME ThreeLeg;  ///< Three-legged operators, used if this is a particle-number-changing operator, i.e. if legs is odd
//  ThreeBodyMENO2B ThreeBodyNO2B; ///< The three body piece of the operator.

  int rank_J; ///< Spherical tensor rank of the operator
  int rank_T; ///< Isotensor rank of the operator
  int parity; ///< Parity of the operator, 0=even 1=odd
  int particle_rank; ///< Maximum particle rank. Should be 2 or 3.
  int legs; ///< The maximum number of particle legs in a diagrammatic representation, e.g. a 2-body operator has 4 legs.

  int E2max; ///< For two-body matrix elements, \f$ e_i + e_j \leq \f$ E2max
  int E3max; ///< For three-body matrix elements, \f$ e_i + e_j + e_k \leq \f$ E3max

  bool hermitian;
  bool antihermitian;
  int nChannels; ///< Number of two-body channels \f$ J,\pi,T_z \f$ associated with the model space

  bool is_reduced; ///< Are the matrix elements stored as reduced matrix elements?



  std::map<std::array<int,3>,std::set<index_t> > OneBodyChannels;  // a set makes more sense for this, because it only contains unique entries
  std::vector< std::set<index_t> > OneBodyChannels_vec; // hopefully this will speed things up.
//  std::set<index_t>& GetOneBodyChannel(int l, int twoj, int twotz) ;
  const std::set<index_t>& GetOneBodyChannel(int l, int twoj, int twotz) const ;
//  std::map<std::array<int,3>,std::vector<index_t> > OneBodyChannels;
  index_t Q_space_orbit; // Orbit with the same quantum numbers as this dagger operator. -1 if it's not a dagger operator.

//  static IMSRGProfiler profiler;
  IMSRGProfiler profiler;


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

  // One body setter/getters
  double GetOneBody(int i,int j) {return OneBody(i,j);};
  void   SetOneBody(int i, int j, double val) ;
  size_t GetTwoBodyDimension(size_t ch_bra, size_t ch_ket){ return TwoBody.GetMatrix(ch_bra, ch_ket).n_cols;};
  double GetTwoBody(size_t ch_bra, size_t ch_ket, size_t i, size_t j);
  void   SetTwoBody(int J1, int p1, int T1, int J2, int p2, int T2, int i, int j, int k, int l, double v);

  void SetE3max(int e){E3max = e;};
  int  GetE3max(){return E3max;};

  // Other setter-getters
//  ModelSpace * GetModelSpace();
  ModelSpace * GetModelSpace() const;  // making this const isn't strictly kosher, but hopefully we shouldn't be in the business of tweaking the modelspace...
  void SetModelSpace(ModelSpace &ms){modelspace = &ms;};

  void Erase(); ///< Set all matrix elements to zero.
  void EraseZeroBody(){ZeroBody = 0;}; ///< set zero-body term to zero
  void EraseOneBody(); ///< set all one-body terms to zero
  void EraseTwoBody(); ///< set all two-body terms to zero
  void EraseThreeBody(); ///< set all two-body terms to zero
  void EraseThreeLeg();

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
  int GetNumberLegs()const {return legs;};
  void SetParticleRank(int pr) { SetNumberLegs(2*pr);};
//  void SetNumberLegs( int l) {legs = l;};
  void SetNumberLegs( int l);
  void SetQSpaceOrbit( index_t q ) {Q_space_orbit = q;};
  index_t GetQSpaceOrbit( ) const {return Q_space_orbit;};

  void MakeReduced();
  void MakeNotReduced();
  bool IsReduced()const { return is_reduced;};
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


  // Undoing normal ordering is equivalent to doing normal ordering with negative occupations.
  // So the occupations na,nb etc are all multiplied by the sign passed to the methods.
  Operator DoNormalOrdering() const; ///< Calls DoNormalOrdering2() or DoNormalOrdering3(), depending on the rank of the operator.
  Operator DoNormalOrdering2(int sign, std::set<index_t> occupied) const; ///< Returns the normal ordered two-body operator
  Operator DoNormalOrdering3(int sign, std::set<index_t> occupied) const; ///< Returns the normal ordered three-body operator
  Operator DoNormalOrderingCore() const; ///< Normal order with respect to core
  Operator DoNormalOrderingFilledValence() const; ///< Normal order with respect to a filled valence space
//  Operator DoNormalOrdering2(int sign=+1) const; ///< Returns the normal ordered two-body operator
//  Operator DoNormalOrdering3(int sign=+1) const; ///< Returns the normal ordered three-body operator
  Operator DoNormalOrderingDagger(int sign, std::set<index_t> occupied) const; ///< Returns the normal ordered dagger operator
  Operator UndoNormalOrdering() const; ///< Returns the operator normal-ordered wrt the vacuum
//  Operator UndoNormalOrdering2() const; ///< Returns the operator normal-ordered wrt the vacuum
//  Operator UndoNormalOrdering2() const;  ///< Returns the operator normal-ordered wrt the vacuum
//  Operator UndoNormalOrdering3() const;  ///< Returns the operator normal-ordered wrt the vacuum
//  Operator UndoNormalOrdering2() const {return this->DoNormalOrdering2(-1);}; ///< Returns the operator normal-ordered wrt the vacuum
//  Operator UndoNormalOrdering3() const {return this->DoNormalOrdering3(-1);}; ///< Returns the operator normal-ordered wrt the vacuum
//  Operator UndoNormalOrderingDagger() const; ///< Returns the operator normal-ordered wrt the vacuum
//  Operator UndoNormalOrderingDagger() const {return this->DoNormalOrderingDagger(-1);}; ///< Returns the operator normal-ordered wrt the vacuum

  Operator Truncate(ModelSpace& ms_new); ///< Returns the operator trunacted to the new model space

  Operator DoIsospinAveraging() const;

  // In principle, these methods should probably be factorized out, but I don't know where to put them...
  double GetMP2_Energy();
//  double GetMP3_Energy();
  std::array<double,3> GetMP3_Energy();
  double MP1_Eval(Operator& );
  double GetMP2_3BEnergy();
  std::array<double,2> GetPPHH_Ladders();

  void PrintTimes(){profiler.PrintAll();};


  double Norm() const;
  double OneBodyNorm() const;
  double TwoBodyNorm() const;
  double ThreeBodyNorm() const;
  double OneLegNorm() const;
  double ThreeLegNorm() const;


  double Trace(int Atrace, int Ztrace) const;

  void ScaleFermiDirac(Operator& H, double T, double Efermi);

  void PrintOneBody() const {OneBody.print();};
  void PrintTwoBody() const {TwoBody.PrintAllMatrices() ;};
  void PrintTwoBody(int ch) const {TwoBody.PrintMatrix(ch,ch);};
  void PrintTwoBody(int ch_bra, int ch_ket) const {TwoBody.PrintMatrix(ch_bra,ch_ket);};

//  arma::vec GetMP2_Impacts() const;
};

/// Non member function, multiply by scalar from left side
Operator operator*(const double lhs, const Operator& rhs);
Operator operator*(const double lhs, const Operator&& rhs);



#endif

