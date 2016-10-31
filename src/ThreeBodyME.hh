#ifndef ThreeBodyME_h
#define ThreeBodyME_h 1

#include "ModelSpace.hh"
#include <fstream>

//typedef double ThreeBME_type;
typedef float ThreeBME_type;

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
 public:
  ModelSpace * modelspace;
//  vector<vector<vector<vector<vector<vector<vector<ThreeBME_type>>>>>>> MatEl; //
  vector<ThreeBME_type> MatEl;
  vector<vector<vector<vector<vector<vector<size_t>>>>>> OrbitIndex; //
  int E3max;
  size_t total_dimension;
  
  ~ThreeBodyME();
  ThreeBodyME();
  ThreeBodyME(ModelSpace*);
  ThreeBodyME(ModelSpace* ms, int e3max);

  void Allocate();

  void SetModelSpace(ModelSpace *ms){modelspace = ms;};

//// Three body setter getters
  vector<pair<size_t,double>> AccessME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n) const;
  ThreeBME_type AddToME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, ThreeBME_type V);
  void   SetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, ThreeBME_type V);
  ThreeBME_type GetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n) const;
  ThreeBME_type GetME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n) const;

///// Some other three body methods

  int SortOrbits(int a_in, int b_in, int c_in, int& a,int& b,int& c) const;
  double RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int J) const;
  void SetE3max(int e){E3max = e;};
  int GetE3max(){return E3max;};

  void Erase(); // set all three-body terms to zero
  void Deallocate();
  size_t size(){return total_dimension * sizeof(ThreeBME_type);};


  void WriteBinary(ofstream&);
  void ReadBinary(ifstream&);

};


#endif
