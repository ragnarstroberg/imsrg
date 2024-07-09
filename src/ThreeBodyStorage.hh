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

#ifndef ThreeBodyStorage_h
#define ThreeBodyStorage_h 1

#include "ModelSpace.hh"
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <memory> // for shared_ptr

namespace {
  inline size_t hash_combine( size_t lhs, size_t rhs ) {
    lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
    return lhs;
  }
}

struct ThreeBodyStorageChannel {
  size_t ch_bra = 0;
  size_t ch_ket = 0;

  ThreeBodyStorageChannel(size_t cbra, size_t cket) : ch_bra(cbra), ch_ket(cket) {}
};

inline bool operator==(const ThreeBodyStorageChannel& a, const ThreeBodyStorageChannel& b) {
  return (a.ch_bra == b.ch_bra) && (a.ch_ket == b.ch_ket);
}

struct ThreeBodyStorageChannelHash {
  size_t operator()(const ThreeBodyStorageChannel& ch) const {
    return hash_combine(ch.ch_bra, ch.ch_ket);
  }
};


// This is used by mutliple implementations.
class OrbitIsospin
{
  public:
    int idx, n, l, j, e;
    OrbitIsospin(int idx, int n, int l, int j) : idx(idx), n(n), l(l), j(j), e(2*n+l) {};
};

// Abstract base class. Concrete implementations inherit from this.
class ThreeBodyStorage
{
 public:

  typedef double ME_type;
  enum Permutation {ABC,BCA,CAB,ACB,BAC,CBA};

  ModelSpace * modelspace = NULL;
  int emax; // usually, this should be the emax of the modelspace, but we might want something smaller.
  int E2max;
  int E3max;
  int lmax;
  int herm=1; // +1 for hermitian, -1 for anti-hermitian
//  size_t total_dimension; // Maybe this should be private?

  int rank_J=0;
  int rank_T=0;
  int parity=0;
  int ISOSPIN_BLOCK_DIMENSION=5;

  bool is_allocated=false;

  std::unordered_map<ThreeBodyStorageChannel,size_t, ThreeBodyStorageChannelHash> ch_start; // map {ch_bra,ch_ket} => index in the array storage array
  std::vector<size_t> ch_dim;


  // Constructors
  ThreeBodyStorage(); //
  ThreeBodyStorage(ModelSpace*); //
  ThreeBodyStorage(ModelSpace* ms, int e3max);
  ThreeBodyStorage(const ThreeBodyStorage& Tbme);
  ThreeBodyStorage(ModelSpace*, int rank_J, int rank_T, int parity);
  ThreeBodyStorage(ModelSpace* ms, int e3max, int rank_J, int rank_T, int parity);

  virtual ~ThreeBodyStorage() {}; // virtual destructor so we can delete instances of derived classes

//  virtual std::shared_ptr<ThreeBodyStorage> Clone() const =0;
  virtual std::unique_ptr<ThreeBodyStorage> Clone() const =0;

  virtual std::string GetStorageMode() const =0;

  // Overloaded operators
//  virtual ThreeBodyStorage& operator*=(const double) =0;
//  virtual ThreeBodyStorage& operator+=(const ThreeBodyStorage&) =0;
//  virtual ThreeBodyStorage& operator-=(const ThreeBodyStorage&) =0;
  virtual void Multiply(const double) =0;
  virtual void Add(const ThreeBodyStorage&) =0;
  virtual void Subtract(const ThreeBodyStorage&) =0;


  virtual void Allocate() =0;

  ME_type NotImplemented(std::string funcname) const;

  virtual ME_type GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const {return NotImplemented(__func__);};
  virtual ME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const {return NotImplemented(__func__);};
//  virtual ME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f) const =0;

  /// getter for tensor operators
  virtual ME_type GetME_pn(  int Jab_in, int j0, int Jde_in, int j1, int a, int b, int c, int d, int e, int f) const {return NotImplemented(__func__);};

  // The setters are only safe when setting in the appropriate formalism.
  // If we're storing in isospin and want to set a pn matrix element, things get messy. Likewise for adding.
  virtual void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V) {NotImplemented(__func__);};
//  virtual void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ME_type V) =0;
  virtual void SetME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ME_type V) {NotImplemented(__func__);};

  virtual void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V) {NotImplemented(__func__);};
//  virtual void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ME_type V) =0;
  virtual void AddToME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ME_type V) {NotImplemented(__func__);};

  // In some cases, for efficiency we may want to set by channel and ket index number, rather than abcdef.
  virtual void AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type V) {NotImplemented(__func__);};
  virtual void SetME_pn_ch(  size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type V) {NotImplemented(__func__);};

  virtual ME_type GetME_pn_ch(size_t chbra, size_t chket, size_t ibra, size_t iket) const {return NotImplemented(__func__);};
//  virtual ME_type GetME_pn_ch(size_t chbra, size_t chket, Ket3& bra, Ket3& ket) const =0; // do we use this??

  // In commutator expressions it often comes up that we want the same matrix element from two operators. This saves an extra recoupling.
  // We implement this here in the base class because it should(?) be the same for all implementations
  virtual std::vector<ME_type> GetME_pn_TwoOps(int Jab, int Jde, int twoJ, int a, int b, int c, int d, int e, int f, const ThreeBodyStorage& X, const ThreeBodyStorage& Y) const ;


  virtual void SetME_iso_no2b(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int T3, ME_type V){NotImplemented(__func__);};
  virtual ME_type GetME_iso_no2b(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int T3) const ;
  virtual ME_type GetME_pn_no2b(int a, int b, int c, int d, int e, int f,  int J2b) const ; // Get ME with Jab=Jde=J2b, summed over total J with weight 2J+1

  virtual void SetME_iso_mono(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int T3, ME_type V){NotImplemented(__func__);};
  virtual ME_type GetME_iso_mono(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int T3) const ;
  virtual ME_type GetME_pn_mono(int a, int b, int c, int d, int e, int f) const ; // Get ME with Jab=Jde=J2b, summed over total J with weight 2J+1


  virtual double Norm() const =0;
  virtual void Erase() =0; // set all three-body terms to zero
  virtual void Deallocate() =0;
  virtual size_t size() const =0;
  virtual void WriteBinary(std::ofstream& f){NotImplemented(__func__);};
  virtual void ReadBinary(std::ifstream& f){NotImplemented(__func__);};
  virtual void WriteFile(std::vector<std::string>& StringInputs, std::vector<int>& IntInputs ){NotImplemented(__func__);};
  virtual void ReadFile( std::vector<std::string>& StringInputs, std::vector<int>& IntInputs ){NotImplemented(__func__);};
  virtual void Print() {NotImplemented(__func__);};
//  virtual int CountAllocations() const =0;
  virtual int CountAllocations() const {NotImplemented(__func__); return 0;};
  bool IsAllocated() const {return is_allocated;};
  void SetHerm(int h) { herm = h; };
  void SetEmax(int e) { emax = e; };
  void SetE3max(int e) { E3max = e; };


  // We define SortOrbits here, but we use a modified definition for isospin, so we make this virtual
  virtual Permutation SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const;
  // These aren't virtual and can be implemented just once
  double RecouplingCoefficient(Permutation recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const;
  int PermutationPhase( Permutation recoupling_case ) const;
  void Permute( Permutation perm, size_t a_in, size_t b_in, size_t c_in, size_t& a_out, size_t& b_out, size_t& c_out );
  std::vector<Permutation> UniquePermutations( size_t a, size_t b, size_t c ) const;

  // Check that a, b, c fulfill certain truncations and restrictions.
  bool IsKetValid(int Jab, int twoJ, size_t a, size_t b, size_t c) const;
  bool IsKetInEMaxTruncations(size_t a, size_t b, size_t c) const;
  bool IsOrbitIn3BodyEMaxTruncation(size_t a) const;
  bool IsOrbitIn3BodyEMaxTruncation(const Orbit& oa) const;
  size_t GetKetIndex_withRecoupling( int Jab, int twoJ, size_t a, size_t b, size_t c, std::vector<size_t>& ibra, std::vector<double>& recouple) const ;



};

#endif
