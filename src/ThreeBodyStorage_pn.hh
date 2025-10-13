#ifndef ThreeBodyStorage_pn_h
#define ThreeBodyStorage_pn_h 1

#include "ThreeBodyStorage.hh"



class ThreeBodyStorage_pn : public ThreeBodyStorage
{
 private:
  typedef double pnME_type;
  std::vector<pnME_type> MatEl; // vector for holding proton-neutron matrix elements

  static bool none_allocated;
  static int number_allocated;
  size_t total_dimension;

 public:


  using ThreeBodyStorage::ThreeBodyStorage; // use the base class constructors
  ThreeBodyStorage_pn( const ThreeBodyStorage_pn& ); // copy constructor
  ~ThreeBodyStorage_pn(); // explicit destructor to count 3 allocations

//  std::shared_ptr<ThreeBodyStorage> Clone() const override ;
  std::unique_ptr<ThreeBodyStorage> Clone() const override ;

  std::string GetStorageMode() const override {return "pn";};

  void Multiply(const double)  override;
  void Add(const ThreeBodyStorage&) override ;
  void Subtract(const ThreeBodyStorage&) override ;

  void Allocate() override;

  ME_type GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const override;
  ME_type GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const override;

  // tensor getters
  ME_type GetME_pn(  int Jab, int j0, int Jde, int j1, int a, int b, int c, int d, int e, int f) const override;



//  void SetME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, pnME_type V) override;
//  void AddToME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ME_type V) override;

  // In some cases, for efficiency we may want to set by channel and ket index number, rather than abcdef.
  void AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type V) override;
  void SetME_pn_ch(  size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type V) override;

  ME_type GetME_pn_ch(size_t chbra, size_t chket, size_t ibra, size_t iket) const override;


  double Norm() const override;
  void Erase() override; // Set all elements to zero
  void Deallocate() override;
  size_t size() const override;
  void WriteBinary(std::ofstream& f) override;
  void ReadBinary(std::ifstream& f) override;

  int CountAllocations() const override {return number_allocated;};
  void Print() override ;

//  void SetHerm(int h) { herm = h; std::cout << __FILE__ << "  setting herm to " << herm; };


  // Internal implementation methods
 private:
  void AccessME(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, size_t& index, int& herm_flip) const;

};



#endif
