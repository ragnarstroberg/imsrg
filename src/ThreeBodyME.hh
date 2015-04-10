#ifndef ThreeBodyME_h
#define ThreeBodyME_h 1

#include "ModelSpace.hh"

/// The three-body piece of an operator, stored in a map of array<int,9> to array<double,5>.
/// The 3BMEs are stored in unnormalized JT coupled form
/// \f$ \langle (abJ_{ab}t_{ab})c | V | (deJ_{de}t_{de})f \rangle_{JT} \f$.
/// To minimize the number of stored matrix elements, only elements with
/// \f$ a\geq b \geq c, a\geq d\geq e \geq f \f$ are stored.
/// The other combinations are obtained on the fly by GetME().
/// The storage format is MatEl[{a,b,c,d,e,f,J,Jab,Jde}][T_index] =
/// \f$ \langle (abJ_{ab}t_{ab})c | V | (deJ_{de}t_{de})f  \rangle_{JT} \f$.
class ThreeBodyME
{
 public:
  ModelSpace * modelspace;

//  map< array<int,9>,array<double,5> >MatEl; //
  vector<vector<vector<vector<vector<vector<vector<double>>>>>>> MatEl; //
//  double******* MatEl; //

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

  void Erase(); // set all three-body terms to zero




};


#endif
