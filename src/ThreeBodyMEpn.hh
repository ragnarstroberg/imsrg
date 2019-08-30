#ifndef ThreeBodyMEpn_h
#define ThreeBodyMEpn_h

#include "ModelSpace.hh"
#include "SymmMatrix.hh"
#include <vector>
#include <map>
#include <unordered_map>
#include <array>

class ThreeBodyMEpn
{
  typedef float ME_type;

 public:
  ModelSpace * modelspace;
  std::map< std::array<size_t,2>, SymmMatrix<ME_type>> MatEl;


//  std::unordered_map<size_t, size_t> OrbitIndexHash; // TODO: reorganize so that we store the pn matrix elements, rather than isospin
  int E3max;
  int emax; // usually, this should be the emax of the modelspace, but we might want something smaller.
  int herm; // +1 for hermitian, -1 for anti-hermitian
  size_t total_dimension;




  void AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBodyMEpn::ME_type matel);
  void AddToME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n, ME_type) ;
  void SetME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBodyMEpn::ME_type matel);
  void SetME_pn(  int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n, ME_type) ;


  ME_type GetME_pn_ch(size_t chbra, size_t chket, size_t ibra, size_t iket) const;
  ME_type GetME_pn_ch(size_t chbra, size_t chket, Ket3& bra, Ket3& ket) const;
  ME_type GetME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n) const;

  ME_type GetME(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n) ;
  void SetME(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ME_type) ;
  void AddToME(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ME_type) ;


  size_t GetKetIndex_withRecoupling( int twoJ, int Jab, size_t a, size_t b, size_t c, std::vector<size_t>& ibra, std::vector<double>& recouple) const ;

  int SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const;

  int PermutationPhase( int recoupling_case ) const;
  double RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const;


  const static int ABC;
  const static int BCA;
  const static int CAB;
  const static int ACB;
  const static int BAC;
  const static int CBA;
  
  ThreeBodyMEpn();
  ThreeBodyMEpn(ModelSpace*);
  ThreeBodyMEpn(ModelSpace* ms, int e3max);
  ThreeBodyMEpn(const ThreeBodyMEpn& tbme);

  void Allocate();
  size_t size();
  void Erase();
  void SetHermitian(){herm = 1;};
  void SetAntiHermitian(){herm = -1;};
  double Norm() const;


  ThreeBodyMEpn& operator*=(const double);
  ThreeBodyMEpn& operator+=(const ThreeBodyMEpn&);
  ThreeBodyMEpn& operator-=(const ThreeBodyMEpn&);

  void WriteBinary(std::ofstream&){}; // Not implemented yet...
  void ReadBinary(std::ifstream&){};

};



#endif
