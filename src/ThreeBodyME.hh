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

  std::shared_ptr<ThreeBodyStorage> threebody_storage;

//  typedef float ThreeBME_type;
//  typedef double ThreeBME_type;
  ModelSpace * modelspace;
//  std::vector<ThreeBME_type> MatEl_iso; // vector for holding isospin matrix elements
//  std::vector<ThreeBME_type> MatEl_pn; // vector for holding proton-neutron matrix elements

  enum Storage_Mode {isospin,pn};
  Storage_Mode storage_mode = isospin; // default to isospin


//  std::unordered_map<size_t, size_t> OrbitIndexHash; // rolls {a,b,c,d,e,f} into a single index for isospin mat el access.


 public:

  // It seems like it should be possible to make these two private, but why bother?
//  std::map<std::array<size_t,2>,size_t> ch_start; // relates {ch_bra,ch_ket} to location in MatEl_pn
//  std::vector<size_t> ch_dim;  // number of 3-body kets in each 3-body channel

  int E3max;
  int emax; // usually, this should be the emax of the modelspace, but we might want something smaller.
  int herm=1; // +1 for hermitian, -1 for anti-hermitian
  size_t total_dimension;

  int rank_J=0;
  int rank_T=0;
  int parity=0;
  int ISOSPIN_BLOCK_DIMENSION=5;

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
  ThreeBodyME& operator*=(const double);
  ThreeBodyME& operator+=(const ThreeBodyME&);
  ThreeBodyME& operator-=(const ThreeBodyME&);
  // TODO Maybe also implmeent * + and - ?


  void Allocate();
//  void Allocate_Isospin(); // These should maybe be private
//  void Allocate_PN(); // These should maybe be private

  void TransformToPN();
  void TransformToIsospin(); // not implemented yet
  void SwitchToPN_and_discard();// Maybe just make this the default behavior when setting PN-mode ?


  // Various ways to access the matrix elements.
  // The getters should work regardless of how the matrix elements are stored internally
  ThreeBME_type GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const;
  ThreeBME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const;
  ThreeBME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f) const;

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



///// Some other three body methods

  double RecouplingCoefficient(ThreeBodyStorage::Permutation recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const;
  int PermutationPhase( ThreeBodyStorage::Permutation recoupling_case ) const;
  void Permute( ThreeBodyStorage::Permutation perm, size_t a_in, size_t b_in, size_t c_in, size_t& a_out, size_t& b_out, size_t& c_out );
  std::vector<ThreeBodyStorage::Permutation> UniquePermutations( size_t a, size_t b, size_t c ) const;

  size_t GetKetIndex_withRecoupling( int Jab, int twoJ, size_t a, size_t b, size_t c, std::vector<size_t>& ibra, std::vector<double>& recouple) const ;

  // setter-getters
  void SetModelSpace(ModelSpace *ms){modelspace = ms;};
  ModelSpace* GetModelSpace(){return modelspace;};

  int GetE3max(){return E3max;};
  int Getemax(){return emax;};
//  void SetE3max(int e){E3max = e;};
//  void Setemax(int e){emax=  e;};
  void SetE3max(int e);
  void Setemax(int e);
  void SetHermitian();
  void SetAntiHermitian();
//  void SetHermitian(){herm = +1;};
//  void SetAntiHermitian(){herm = -1;};

  bool Is_PN_Mode() const {return (storage_mode == pn);};
  bool Is_Isospin_Mode()const {return (storage_mode == isospin);};
  bool IsAllocated() const {return threebody_storage->IsAllocated();};

  std::map<std::array<size_t,2>,size_t>& Get_ch_start();
  std::vector<size_t>& Get_ch_dim();

  double Norm() const;
  void Erase(); // set all three-body terms to zero
  void Deallocate();
  size_t size();
  void Print();
  std::string GetStorageMode();


  void WriteBinary(std::ofstream&);
  void ReadBinary(std::ifstream&);

};


#endif
