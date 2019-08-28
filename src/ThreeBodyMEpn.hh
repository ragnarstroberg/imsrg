#ifndef ThreeBodyMEpn_h
#define ThreeBodyMEpn_h



class ThreeBodyMEpn
{
  typedef float ME_type;

 public:
  ModelSpace * modelspace;
  std::map< array<size_t,2>, SymmMatrix<ME_type>> MatEl;


  void AddToME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n, ME_type) ;
  void SetME_pn(  int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n, ME_type) ;


  ME_type GetME_pn(size_t chbra, size_t chket, size_t ibra, size_t iket) const;
  ME_type GetME_pn(size_t chbra, size_t chket, Ket3& bra, Ket3& ket) const;
  ME_type GetME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n) const;



  size_t GetKetIndex_withRecoupling( int Jab, size_t a, size_t b, size_t c, size_t& ch, double& recouple);

  int ThreeBodyMEpn::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const;


  std::unordered_map<size_t, size_t> OrbitIndexHash; // TODO: reorganize so that we store the pn matrix elements, rather than isospin
  int E3max;
  int emax; // usually, this should be the emax of the modelspace, but we might want something smaller.
  int herm; // +1 for hermitian, -1 for anti-hermitian
  size_t total_dimension;
  const static int ABC;
  const static int BCA;
  const static int CAB;
  const static int ACB;
  const static int BAC;
  const static int CBA;
  
  ThreeBodyMEpn();
  ThreeBodyMEpn(ModelSpace*);
  ThreeBodyMEpn(ModelSpace* ms, int e3max);
  ThreeBodyMEpn(const ThreeBodyME& tbme);

  void Allocate();


  ThreeBodyMEpn& operator*=(const double);
  ThreeBodyMEpn& operator+=(const ThreeBodyMEpn&);
  ThreeBodyMEpn& operator-=(const ThreeBodyMEpn&);


}



#endif
