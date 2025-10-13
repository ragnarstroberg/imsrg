///////////////////////////////////////////////////////////////////////////////////
//    ThreeBodyME.hh, part of  imsrg++
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

#ifndef ThreeBodyME_h
#define ThreeBodyME_h 1

#include "ModelSpace.hh"
#include "ThreeBodyStorage.hh"
#include "ThreeBodyStorage_iso.hh"
#include "ThreeBodyStorage_pn.hh"
#include "ThreeBodyStorage_no2b.hh"
#include "ThreeBodyStorage_mono.hh"
#include <fstream>
#include <unordered_map>
#include <memory> // for shared_ptr

typedef double ThreeBME_type;
//typedef float ThreeBME_type;

/// The three-body piece of an operator, stored in nested vectors.
/// The 3BMEs are stored in unnormalized JT coupled form
/// \f$ \langle (abJ_{ab}t_{ab})c | V | (deJ_{de}t_{de})f \rangle_{JT} \f$.
/// To minimize the number of stored matrix elements, only elements with
/// \f$ a\geq b \geq c, a\geq d\geq e \geq f \f$ are stored.
/// The other combinations are obtained on the fly by GetME().
/// The storage format is MatEl[{a,b,c,d,e,f,J,Jab,Jde}][T_index] =
/// \f$ \langle (abJ_{ab}t_{ab})c | V | (deJ_{de}t_{de})f  \rangle_{JT} \f$.
class ThreeBodyME
{
 private:

//  std::shared_ptr<ThreeBodyStorage> threebody_storage;
  std::unique_ptr<ThreeBodyStorage> threebody_storage;

  ModelSpace * modelspace;
//  enum Storage_Mode {isospin,pn};
//  Storage_Mode storage_mode = isospin; // default to isospin

 public:

  // It seems like it should be possible to make these two private, but why bother?
//  std::map<std::array<size_t,2>,size_t> ch_start; // relates {ch_bra,ch_ket} to location in MatEl_pn
//  std::vector<size_t> ch_dim;  // number of 3-body kets in each 3-body channel

  int E3max=0;
  int emax=0; // usually, this should be the emax of the modelspace, but we might want something smaller.
  int herm=1; // +1 for hermitian, -1 for anti-hermitian
  size_t total_dimension;

  int rank_J=0;
  int rank_T=0;
  int parity=0;
  int ISOSPIN_BLOCK_DIMENSION=5; //default is 5, but can be different for charge-changing operators

//  bool isospin_is_allocated = false;
//  bool pn_is_allocated = false;
//  static bool none_allocated;

//  enum Permutation {ABC,BCA,CAB,ACB,BAC,CBA};
  

  // Constructors
//  ~ThreeBodyME(); // dont think we explicitly need this?
  ThreeBodyME(); //
  ThreeBodyME(ModelSpace*); //
  ThreeBodyME(ModelSpace* ms, int e3max);
  ThreeBodyME(const ThreeBodyME& Tbme);
  ThreeBodyME(ModelSpace*, int rank_J, int rank_T, int parity);
  ThreeBodyME(ModelSpace* ms, int e3max, int rank_J, int rank_T, int parity);

  // Overloaded operators
  ThreeBodyME& operator=(const ThreeBodyME&);
  ThreeBodyME& operator*=(const double);
  ThreeBodyME& operator+=(const ThreeBodyME&);
  ThreeBodyME& operator-=(const ThreeBodyME&);
  // TODO Maybe also implmeent * + and - ?


  void Allocate();

  void TransformToPN();
  void TransformToIsospin(); // not implemented yet
  void SwitchToPN_and_discard();// Maybe just make this the default behavior when setting PN-mode ?
  void SetMode(std::string mode);
//  std::string GetMode() const;


  // Various ways to access the matrix elements.
  // The getters should work regardless of how the matrix elements are stored internally
  ThreeBME_type GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const;
  ThreeBME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const;
  ThreeBME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f) const;

  // tensor geteters
  ThreeBME_type GetME_pn(  int Jab_in, int j0, int Jde_in, int j1, int a, int b, int c, int d, int e, int f) const;

  // The setters are only safe when setting in the appropriate formalism.
  // If we're storing in isospin and want to set a pn matrix element, things get messy. Likewise for adding.
  void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ThreeBME_type V);
  void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ThreeBME_type V);
  void SetME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBME_type V);

  void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ThreeBME_type V);
  void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ThreeBME_type V);
  void AddToME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBME_type V);


  // In some cases, for efficiency we may want to set by channel and ket index number, rather than abcdef.
  void AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBME_type V);
  void SetME_pn_ch(  size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBME_type V);

  ThreeBME_type GetME_pn_ch(size_t chbra, size_t chket, size_t ibra, size_t iket) const;
//  ThreeBME_type GetME_pn_ch(size_t chbra, size_t chket, Ket3& bra, Ket3& ket) const; // do we use this??

  // In commutator expressions it often comes up that we want the same matrix element from two operators. This saves an extra recoupling.
  std::vector<double> GetME_pn_TwoOps(int Jab, int Jde, int twoJ, int a, int b, int c, int d, int e, int f, const ThreeBodyME& X, const ThreeBodyME& Y) const;

  ThreeBME_type GetME_pn_no2b(int a, int b, int c, int d, int e, int f,  int J2b) const;
  ThreeBME_type GetME_pn_mono(int a, int b, int c, int d, int e, int f) const;


///// Some other three body methods

  double RecouplingCoefficient(ThreeBodyStorage::Permutation recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const;
  int PermutationPhase( ThreeBodyStorage::Permutation recoupling_case ) const;
  void Permute( ThreeBodyStorage::Permutation perm, size_t a_in, size_t b_in, size_t c_in, size_t& a_out, size_t& b_out, size_t& c_out );
  std::vector<ThreeBodyStorage::Permutation> UniquePermutations( size_t a, size_t b, size_t c ) const;

  // Check that a, b, c fulfill certain truncations and restrictions.
  bool IsKetValid(int Jab, int twoJ, size_t a, size_t b, size_t c) const;
  // Check that a, b, c fulfill the emax truncations.
  bool IsKetInEMaxTruncations(size_t a, size_t b, size_t c) const;
  // Check that orbit a fulfills the 3-body emax truncation.
  bool IsOrbitIn3BodyEMaxTruncation(size_t a) const;
  // Check that orbit a fulfills the 3-body emax truncation.
  bool IsOrbitIn3BodyEMaxTruncation(const Orbit& oa) const;
  size_t GetKetIndex_withRecoupling( int Jab, int twoJ, size_t a, size_t b, size_t c, std::vector<size_t>& ibra, std::vector<double>& recouple) const ;

  // setter-getters
  void SetModelSpace(ModelSpace *ms){modelspace = ms;};
  ModelSpace* GetModelSpace(){return modelspace;};

  int GetE3max(){return E3max;};
  int Getemax(){return emax;};
  void SetE3max(int e);
  void Setemax(int e);
  void SetHermitian();
  void SetAntiHermitian();
  bool IsHermitian() const { return herm==1;};

//  bool Is_PN_Mode() const {return (storage_mode == pn);};
//  bool Is_Isospin_Mode()const {return (storage_mode == isospin);};
  bool Is_PN_Mode() const {return (threebody_storage->GetStorageMode() == "pn");};
  bool Is_Isospin_Mode() const {return (threebody_storage->GetStorageMode() == "isospin");};
  bool IsAllocated() const {return threebody_storage->IsAllocated();};

  std::unordered_map<ThreeBodyStorageChannel,size_t,ThreeBodyStorageChannelHash>& Get_ch_start() const;
  std::vector<size_t>& Get_ch_dim();

  double Norm() const;
  void Erase(); // set all three-body terms to zero
  void Deallocate();
  size_t size();
  void Print() const;
  std::string GetStorageMode() const ;


  void WriteBinary(std::ofstream&);
  void ReadBinary(std::ifstream&);

  void WriteFile(std::vector<std::string> StringInputs, std::vector<int> IntInputs ) const;
  void ReadFile( std::vector<std::string> StringInputs, std::vector<int> IntInputs );

};


#endif
