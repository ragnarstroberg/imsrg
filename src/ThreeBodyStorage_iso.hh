
#ifndef ThreeBodyStorage_iso_h
#define ThreeBodyStorage_iso_h 1
#include "ThreeBodyStorage.hh"



class ThreeBodyStorage_iso : public ThreeBodyStorage
{
 private:
  typedef float isoME_type;
//  typedef double isoME_type;
  std::vector<isoME_type> MatEl; // vector for holding proton-neutron matrix elements

  std::unordered_map<size_t, size_t> OrbitIndexHash; // rolls {a,b,c,d,e,f} into a single index for isospin mat el access.
  static bool none_allocated;
  size_t total_dimension;

 public:

  using ThreeBodyStorage::ThreeBodyStorage;
  ThreeBodyStorage_iso( const ThreeBodyStorage_iso& );
//  ~ThreeBodyStorage_iso();

//  std::shared_ptr<ThreeBodyStorage> Clone() const override;
  std::unique_ptr<ThreeBodyStorage> Clone() const override;

  std::string GetStorageMode() const override {return "isospin";};

  void Multiply(const double) override;
  void Add(const ThreeBodyStorage&)  override;
  void Subtract(const ThreeBodyStorage&)  override;

  void Allocate() override;

  ME_type GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const override;
  ME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const override;

  // The setters are only safe when setting in the appropriate formalism.
  // If we're storing in isospin and want to set a pn matrix element, things get messy. Likewise for adding.
  void SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V) override;

  void AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V) override;

//  ME_type GetME_pn_ch(size_t chbra, size_t chket, size_t ibra, size_t iket) const override; // not implemented


  double Norm() const override;
  void Erase() override; // Set all elements to zero
  void Deallocate() override;
  size_t size() const override;
  void WriteBinary(std::ofstream& f) override;
  void ReadBinary(std::ifstream& f) override;
//  void Print() override {};

//  void SetHerm(int h) { herm = h; std::cout << __FILE__ << "  setting herm to " << herm; };

  // We override the SortOrbits method because we only use even indices in isospin mode (and a different ordering).
  ThreeBodyStorage::Permutation SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const override;

  // Internal implementation methods
 private:
  std::vector<std::pair<size_t,double>> AccessME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in) const;

    // Hashing and un-hashing for orbit index to lookup index in isospin mode
  size_t KeyHash(size_t,size_t,size_t,size_t,size_t,size_t) const;
  void KeyUnhash(size_t& key, size_t& a, size_t& b, size_t& c, size_t& d, size_t& e, size_t& f) const;


};

#endif
