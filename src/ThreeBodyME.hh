#ifndef ThreeBodyME_h
#define ThreeBodyME_h 1

#include "ModelSpace.hh"


class ThreeBodyME
{
 public:
  ModelSpace * modelspace;

  map< array<int,9>,array<double,5> >MatEl; //

  int E3max;
  
  ThreeBodyME();
  ThreeBodyME(ModelSpace*);
  ThreeBodyME(ModelSpace* ms, int e3max);

  void Allocate();

  void SetModelSpace(ModelSpace *ms){modelspace = ms;};

//// Three body setter getters
  double AddToME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, double V);
  void   SetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, double V);
  double GetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n);
  double GetME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n);

///// Some other three body methods

  int SortOrbits(int a_in, int b_in, int c_in, int& a,int& b,int& c);
  double RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int J);
  void SetE3max(int e){E3max = e;};
  int GetE3max(){return E3max;};

  void Erase(); // set all two-body terms to zero




};


#endif
