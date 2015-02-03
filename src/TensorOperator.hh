
#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include <armadillo>

//#define JMAX 30

using namespace std;

class Operator
{
 public:
  //Fields
  float ZeroBody;
  arma::mat OneBody;
//  vector<arma::mat> TwoBody;
  vector<vector<arma::mat> > TwoBody;

  int rank_J;
  int rank_T;
  int parity;

  vector<vector<int> > TwoBodyTensorChannels;

  //Constructors
  // In the future, consider using C++11 rvalues / move constructor to avoid copies in certain cases
  Operator();
//  Operator(ModelSpace&);
  Operator(ModelSpace&, int Jrank, int Trank, int Parity);
  Operator( const Operator& rhs);

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
  float GetOneBody(int i,int j) {return OneBody(i,j);};
  void SetOneBody(int i, int j, float val) { OneBody(i,j) = val;};

  //TwoBody setter/getters
  double GetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d) const;
  void   SetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme);
  double GetTBME(int ch_bra, int ch_ket, Ket &bra, Ket &ket) const;
  void   SetTBME(int ch_bra, int ch_ket, Ket &bra, Ket& ket, double tbme);
  void   AddToTBME(int ch_bra, int ch_ket, Ket &bra, Ket& ket, double tbme);
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
  double GetTBME(int j, int p, int t, Ket& bra, Ket& ket) const;
  void   SetTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme);
  void   AddToTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme);
  double GetTBME(int j, int p, int t, int a, int b, int c, int d) const;
  void   SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int j, int p, int t, int a, int b, int c, int d, double tbme);


  double GetTBMEmonopole(int a, int b, int c, int d) const;
  double GetTBMEmonopole(Ket & bra, Ket & ket) const;


  // Other setter-getters
  ModelSpace * GetModelSpace() const {return modelspace;};

  void EraseZeroBody(){ZeroBody = 0;}; // set zero-body term to zero
  void EraseOneBody(){OneBody.zeros();}; // set all one-body terms to zero
  void EraseTwoBody(); // set all two-body terms to zero

  void SetHermitian() {hermitian=true;antihermitian=false;};
  void SetAntiHermitian() {antihermitian=true;hermitian=false;};
  void SetNonHermitian() {antihermitian=false;hermitian=false;};
  bool IsHermitian()const {return hermitian;};
  bool IsAntiHermitian()const {return antihermitian;};
  bool IsNonHermitian()const {return not (hermitian or antihermitian);};



  void ScaleZeroBody(double x);
  void ScaleOneBody(double x);
  void ScaleTwoBody(double x);

  // The actually interesting methods
  Operator Commutator(const Operator& opright) const;
  Operator BCH_Product( const Operator& ) const ; 
  Operator BCH_Transform( const Operator& ) const; 
  Operator DoNormalOrdering(); 
  void Eye(); // set to identity operator
  void CalculateKineticEnergy();

  Operator CommutatorScalarScalar(const Operator& opright) const;
  Operator CommutatorScalarTensor(const Operator& opright) const;
  Operator CommutatorTensorTensor(const Operator& opright) const;

  double Norm() const;
  double OneBodyNorm() const;
  double TwoBodyNorm() const;



  void PrintOneBody() const {OneBody.print();};
  void PrintTwoBody(int ch) const {TwoBody[ch].print();};


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
  void comm222_phst(const Operator& opright, Operator& opout) const;
  void comm222_pp_hh_221st(const Operator& opright, Operator& opout) const;


  void comm111tt(const Operator& opright, Operator& opout) const;
  void comm121tt(const Operator& opright, Operator& opout) const;
  void comm221tt(const Operator& opright, Operator& opout) const;
  void comm122tt(const Operator& opright, Operator& opout) const;
  void comm222_pp_hhtt(const Operator& opright, Operator& opout) const;
  void comm222_phtt(const Operator& opright, Operator& opout) const;
  void comm222_pp_hh_221tt(const Operator& opright, Operator& opout) const;

};

// Non member function, multiply by scalar from left side
Operator operator*(const double lhs, const Operator& rhs);

#endif

