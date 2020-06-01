#ifndef ThreeBodyStorage_pn_h
#define ThreeBodyStorage_pn_h 1

#include "ThreeBodyStorage.hh"



class ThreeBodyStorage_pn : public ThreeBodyStorage
{
 private:
  typedef double pnME_type;
  std::vector<pnME_type> MatEl; // vector for holding proton-neutron matrix elements

  static bool none_allocated;
  size_t total_dimension;

 public:


  using ThreeBodyStorage::ThreeBodyStorage; // use the base class constructors
  ThreeBodyStorage_pn( const ThreeBodyStorage_pn& ); // copy constructor

  std::shared_ptr<ThreeBodyStorage> Clone() const ;

  std::string GetStorageMode() const {return "pn";};

//  ThreeBodyStorage_pn& operator*=(const double);
//  ThreeBodyStorage_pn& operator+=(const ThreeBodyStorage_pn&);
//  ThreeBodyStorage_pn& operator-=(const ThreeBodyStorage_pn&);
  void Multiply(const double) ;
  void Add(const ThreeBodyStorage&) ;
  void Subtract(const ThreeBodyStorage&) ;

  void Allocate();

  ME_type GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const;
  ME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const;
//  ME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f) const;

  // The setters are only safe when setting in the appropriate formalism.
  // If we're storing in isospin and want to set a pn matrix element, things get messy. Likewise for adding.
  void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V); // not implemented
//  void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ME_type V);
  void SetME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, pnME_type V);

  void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V); // not implemented
//  void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ME_type V);
  void AddToME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ME_type V);

  // In some cases, for efficiency we may want to set by channel and ket index number, rather than abcdef.
  void AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type V);
  void SetME_pn_ch(  size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type V);

  ME_type GetME_pn_ch(size_t chbra, size_t chket, size_t ibra, size_t iket) const;
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
  void AccessME(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, size_t& index, int& herm_flip) const;

};



#endif
