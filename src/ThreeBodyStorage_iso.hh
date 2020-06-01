
#ifndef ThreeBodyStorage_iso_h
#define ThreeBodyStorage_iso_h 1
#include "ThreeBodyStorage.hh"



class ThreeBodyStorage_iso : public ThreeBodyStorage
{
 private:
  typedef float isoME_type;
  std::vector<isoME_type> MatEl; // vector for holding proton-neutron matrix elements

  std::unordered_map<size_t, size_t> OrbitIndexHash; // rolls {a,b,c,d,e,f} into a single index for isospin mat el access.
  static bool none_allocated;
  size_t total_dimension;

 public:

  using ThreeBodyStorage::ThreeBodyStorage;
  ThreeBodyStorage_iso( const ThreeBodyStorage_iso& );

  std::shared_ptr<ThreeBodyStorage> Clone() const;

  std::string GetStorageMode() const {return "isospin";};

//  ThreeBodyStorage_iso& operator*=(const double);
//  ThreeBodyStorage_iso& operator+=(const ThreeBodyStorage_iso&);
//  ThreeBodyStorage_iso& operator-=(const ThreeBodyStorage_iso&);

  void Multiply(const double) ;
  void Add(const ThreeBodyStorage&) ;
  void Subtract(const ThreeBodyStorage&) ;

  void Allocate();

  ME_type GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const;
  ME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const;
//  ME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f) const;

  // The setters are only safe when setting in the appropriate formalism.
  // If we're storing in isospin and want to set a pn matrix element, things get messy. Likewise for adding.
  void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V);
//  void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ME_type V);
  void SetME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ME_type V); // not implemented

  void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V);
//  void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ME_type V);
  void AddToME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ME_type V); // not implemented

  // In some cases, for efficiency we may want to set by channel and ket index number, rather than abcdef.
  void AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type V); // not implemented
  void SetME_pn_ch(  size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type V); // not implemented

  ME_type GetME_pn_ch(size_t chbra, size_t chket, size_t ibra, size_t iket) const; // not implemented
//  pnME_type GetME_pn_PN_ch(size_t chbra, size_t chket, Ket3& bra, Ket3& ket) const; // do we use this??


  double Norm() const;
  void Erase(); // Set all elements to zero
  void Deallocate();
  size_t size();
  void WriteBinary(std::ofstream& f);
  void ReadBinary(std::ifstream& f);
  void Print() {};

  void SetHerm(int h) { herm = h; std::cout << __FILE__ << "  setting herm to " << herm; };

  // Internal implementation methods
 private:
  std::vector<std::pair<size_t,double>> AccessME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in) const;

  // We overload the SortOrbits method because we only use even indices in isospin mode.
  ThreeBodyStorage::Permutation SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const;
    // Hashing and un-hashing for orbit index to lookup index in isospin mode
  size_t KeyHash(size_t,size_t,size_t,size_t,size_t,size_t) const;
  void KeyUnhash(size_t& key, size_t& a, size_t& b, size_t& c, size_t& d, size_t& e, size_t& f) const;


};

#endif
