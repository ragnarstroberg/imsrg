
#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include <armadillo>
#include <vector>
#include <map>

//#define JMAX 30

using namespace std;

typedef uint64_t orbindx3_t;

class Operator
{
 public:
  //Fields
  double ZeroBody;
  arma::mat OneBody;
//  vector<arma::mat> TwoBody;
//  vector<vector<arma::mat> > TwoBody;
  vector<map<int,arma::mat> > TwoBody;

  map< orbindx3_t, vector< vector<double> > > ThreeBody;

  int rank_J;
  int rank_T;
  int parity;
  int particle_rank;

  vector<vector<int> > TwoBodyTensorChannels;

  int E2max; // I don't do anything with this yet...
  int E3max;



  //Constructors
  // In the future, consider using C++11 rvalues / move constructor to avoid copies in certain cases
  Operator();
  Operator(ModelSpace&);
  Operator(ModelSpace&, int Jrank, int Trank, int Parity, int part_rank);
  Operator( const Operator& rhs);

  void AllocateTwoBody();
  void AllocateThreeBody();

  //Overloaded operators
  Operator& operator=( const Operator& rhs);
  Operator& operator+=( const Operator& rhs);
  Operator operator+( const Operator& rhs) const;
  Operator& operator-=( const Operator& rhs);
  Operator operator-( const Operator& rhs) const;
  Operator& operator*=( const double rhs);
  Operator operator*( const double rhs) const;
  Operator& operator/=( const double rhs);
  Operator operator/( const double rhs) const;

  //Methods
  // One body setter/getters
  double GetOneBody(int i,int j) {return OneBody(i,j);};
  void SetOneBody(int i, int j, double val) { OneBody(i,j) = val;};

  //TwoBody setter/getters
  double GetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d) const;
  void   SetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme);
  double GetTBME(int ch_bra, int ch_ket, Ket &bra, Ket &ket) const;
  void   SetTBME(int ch_bra, int ch_ket, Ket &bra, Ket& ket, double tbme);
  void   AddToTBME(int ch_bra, int ch_ket, Ket &bra, Ket& ket, double tbme);
  double GetTBME(int ch_bra, int ch_ket, int ibra, int iket) const;
  void   SetTBME(int ch_bra, int ch_ket, int ibra, int iket, double tbme);
  void   AddToTBME(int ch_bra, int ch_ket, int ibra, int iket, double tbme);
  double GetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket) const;
  void   SetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket, double tbme);
  void   AddToTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket, double tbme);
  double GetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d) const;
  void   SetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d, double tbme);

  // Scalar setters/getters for backwards compatibility
  double GetTBME(int ch, int a, int b, int c, int d) const;
  void   SetTBME(int ch, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int ch, int a, int b, int c, int d, double tbme);
  double GetTBME(int ch, Ket &bra, Ket &ket) const;
  void   SetTBME(int ch, Ket &bra, Ket& ket, double tbme);
  void   AddToTBME(int ch, Ket &bra, Ket& ket, double tbme);
  double GetTBME(int ch, int ibra, int iket) const;
  void   SetTBME(int ch, int ibra, int iket, double tbme);
  void   AddToTBME(int ch, int ibra, int iket, double tbme);
  double GetTBME(int j, int p, int t, Ket& bra, Ket& ket) const;
  void   SetTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme);
  void   AddToTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme);
  double GetTBME(int j, int p, int t, int a, int b, int c, int d) const;
  void   SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int j, int p, int t, int a, int b, int c, int d, double tbme);


  double GetTBMEmonopole(int a, int b, int c, int d) const;
  double GetTBMEmonopole(Ket & bra, Ket & ket) const;


//// Three body setter getters
  double GetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n);
  double GetThreeBodyME_pn(int Jab_in, int Jde_in, int J2, int Tz, int i, int j, int k, int l, int m, int n);
  void   SetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, double V);


///// Some other three body methods
  orbindx3_t GetThreeBodyOrbitIndex(int a, int b, int c, int d, int e, int f);

  void GetOrbitsFromThreeBodyIndex(orbindx3_t orb_indx, int& a, int& b, int& c, int& d, int& e, int& f);

  void SortThreeBodyOrbits(int& a,int& b,int& c);
  double ThreeBodyRecouplingCoefficient(int a_in, int b_in, int c_in, int a, int b, int c, int Jab_in, int Jab, int J ,char j_or_t);
  void SetE3max(int e){E3max = e;};




  // Other setter-getters
  ModelSpace * GetModelSpace() const {return modelspace;};

  void Erase();
  void EraseZeroBody(){ZeroBody = 0;}; // set zero-body term to zero
  void EraseOneBody(){OneBody.zeros();}; // set all one-body terms to zero
  void EraseTwoBody(); // set all two-body terms to zero

  void SetHermitian() {hermitian=true;antihermitian=false;};
  void SetAntiHermitian() {antihermitian=true;hermitian=false;};
  void SetNonHermitian() {antihermitian=false;hermitian=false;};
  bool IsHermitian()const {return hermitian;};
  bool IsAntiHermitian()const {return antihermitian;};
  bool IsNonHermitian()const {return not (hermitian or antihermitian);};
  int GetParticleRank()const {return particle_rank;};
  int GetJRank()const {return rank_J;};
  int GetTRank()const {return rank_T;};
  int GetParity()const {return parity;};



  void ScaleZeroBody(double x);
  void ScaleOneBody(double x);
  void ScaleTwoBody(double x);

  // The actually interesting methods
  Operator DoNormalOrdering(); 
  Operator DoNormalOrdering2();
  Operator DoNormalOrdering3();

  Operator Commutator(const Operator& opright) const;
  Operator BCH_Product( const Operator& ) const ; 
  Operator BCH_Transform( const Operator& ) const; 

  void CalculateKineticEnergy();
  void Eye(); // set to identity operator

  Operator CommutatorScalarScalar(const Operator& opright) const;
  Operator CommutatorScalarTensor(const Operator& opright) const;
  Operator CommutatorTensorTensor(const Operator& opright) const;

  double Norm() const;
  double OneBodyNorm() const;
  double TwoBodyNorm() const;



  void PrintOneBody() const {OneBody.print();};
  void PrintTwoBody(int ch) const {TwoBody.at(ch).at(ch).print();};


 //private:
  //Fields
  ModelSpace * modelspace;
  bool hermitian;
  bool antihermitian;
  int nChannels;
  static double bch_transform_threshold;
  static double bch_product_threshold;

  //Methods
  void Copy(const Operator& rhs);

  static void Set_BCH_Transform_Threshold(double x){bch_transform_threshold=x;};
  static void Set_BCH_Product_Threshold(double x){bch_product_threshold=x;};

  
  void DoPandyaTransformation(Operator&) const;
  void CalculateCrossCoupled(vector<arma::mat>&, vector<arma::mat>&) const; 

  void comm110ss(const Operator& opright, Operator& opout) const;
  void comm220ss(const Operator& opright, Operator& opout) const;
  void comm111ss(const Operator& opright, Operator& opout) const;
  void comm121ss(const Operator& opright, Operator& opout) const;
  void comm221ss(const Operator& opright, Operator& opout) const;
  void comm122ss(const Operator& opright, Operator& opout) const;
  void comm222_pp_hhss(const Operator& opright, Operator& opout) const;
  void comm222_phss(const Operator& opright, Operator& opout) const;
  void comm222_pp_hh_221ss(const Operator& opright, Operator& opout) const;

// make st and tt commutators

  void comm111st(const Operator& opright, Operator& opout) const;
  void comm121st(const Operator& opright, Operator& opout) const;
  void comm221st(const Operator& opright, Operator& opout) const;
  void comm122st(const Operator& opright, Operator& opout) const;
  void comm222_pp_hhst(const Operator& opright, Operator& opout) const;
  void comm222_pp_hh_221st(const Operator& opright, Operator& opout) const;
  void comm222_phst(const Operator& opright, Operator& opout) const;
  void comm222_phst_pandya(const Operator& opright, Operator& opout) const;

/*
  void comm111tt(const Operator& opright, Operator& opout) const;
  void comm121tt(const Operator& opright, Operator& opout) const;
  void comm221tt(const Operator& opright, Operator& opout) const;
  void comm122tt(const Operator& opright, Operator& opout) const;
  void comm222_pp_hhtt(const Operator& opright, Operator& opout) const;
  void comm222_phtt(const Operator& opright, Operator& opout) const;
  void comm222_pp_hh_221tt(const Operator& opright, Operator& opout) const;
*/

};

// Non member function, multiply by scalar from left side
Operator operator*(const double lhs, const Operator& rhs);

#endif

