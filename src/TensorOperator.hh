
#ifndef TensorOperator_h
#define TensorOperator_h 1

#include "ModelSpace.hh"
#include <armadillo>

//#define JMAX 30

using namespace std;

class TensorOperator
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
  TensorOperator();
//  TensorOperator(ModelSpace&);
  TensorOperator(ModelSpace&, int Jrank, int Trank, int Parity);
  TensorOperator( const TensorOperator& rhs);

  //Overloaded operators
  TensorOperator& operator=( const TensorOperator& rhs);
  TensorOperator& operator+=( const TensorOperator& rhs);
  TensorOperator operator+( const TensorOperator& rhs) const;
  TensorOperator& operator-=( const TensorOperator& rhs);
  TensorOperator operator-( const TensorOperator& rhs) const;
  TensorOperator& operator*=( const double rhs);
  TensorOperator operator*( const double rhs) const;
  TensorOperator& operator/=( const double rhs);
  TensorOperator operator/( const double rhs) const;

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
  TensorOperator Commutator(const TensorOperator& opright) const;
  TensorOperator BCH_Product( const TensorOperator& ) const ; 
  TensorOperator BCH_Transform( const TensorOperator& ) const; 
  TensorOperator DoNormalOrdering(); 
  void Eye(); // set to identity operator
  void CalculateKineticEnergy();

  TensorOperator CommutatorScalarScalar(const TensorOperator& opright) const;
  TensorOperator CommutatorScalarTensor(const TensorOperator& opright) const;
  TensorOperator CommutatorTensorTensor(const TensorOperator& opright) const;

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
  void Copy(const TensorOperator& rhs);

  static void Set_BCH_Transform_Threshold(double x){bch_transform_threshold=x;};
  static void Set_BCH_Product_Threshold(double x){bch_product_threshold=x;};

  
  void CalculateCrossCoupled(vector<arma::mat>&, vector<arma::mat>&) const; 

  void comm110ss(const TensorOperator& opright, TensorOperator& opout) const;
  void comm220ss(const TensorOperator& opright, TensorOperator& opout) const;
  void comm111ss(const TensorOperator& opright, TensorOperator& opout) const;
  void comm121ss(const TensorOperator& opright, TensorOperator& opout) const;
  void comm221ss(const TensorOperator& opright, TensorOperator& opout) const;
  void comm122ss(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_pp_hhss(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_phss(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_pp_hh_221ss(const TensorOperator& opright, TensorOperator& opout) const;

// make st and tt commutators

  void comm111st(const TensorOperator& opright, TensorOperator& opout) const;
  void comm121st(const TensorOperator& opright, TensorOperator& opout) const;
  void comm221st(const TensorOperator& opright, TensorOperator& opout) const;
  void comm122st(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_pp_hhst(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_phst(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_pp_hh_221st(const TensorOperator& opright, TensorOperator& opout) const;


  void comm111tt(const TensorOperator& opright, TensorOperator& opout) const;
  void comm121tt(const TensorOperator& opright, TensorOperator& opout) const;
  void comm221tt(const TensorOperator& opright, TensorOperator& opout) const;
  void comm122tt(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_pp_hhtt(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_phtt(const TensorOperator& opright, TensorOperator& opout) const;
  void comm222_pp_hh_221tt(const TensorOperator& opright, TensorOperator& opout) const;

};

// Non member function, multiply by scalar from left side
TensorOperator operator*(const double lhs, const TensorOperator& rhs);

#endif

