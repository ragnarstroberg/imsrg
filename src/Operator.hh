
#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include <armadillo>

//#define JMAX 30

//class ModelSpace;

using namespace std;

class Operator
{
 public:
  //Fields
  float ZeroBody;
  arma::mat OneBody;
  vector<arma::mat> TwoBody;

  //Constructors
  // In the future, consider using C++11 rvalues / move constructor to avoid copies in certain cases
  Operator();
  Operator(ModelSpace&);
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

//  double GetTBME_NoPhase(int ch, int a, int b, int c, int d) const;

  double GetTBMEmonopole(int a, int b, int c, int d) const;
  double GetTBMEmonopole(Ket & bra, Ket & ket) const;


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

  double Norm() const;
  double OneBodyNorm() const;
  double TwoBodyNorm() const;

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

  void comm110(const Operator& opright, Operator& opout) const;
  void comm220(const Operator& opright, Operator& opout) const;
  void comm111(const Operator& opright, Operator& opout) const;
  void comm121(const Operator& opright, Operator& opout) const;
  void comm221(const Operator& opright, Operator& opout) const;
  void comm122(const Operator& opright, Operator& opout) const;
  void comm222_pp_hh(const Operator& opright, Operator& opout) const;
  void comm222_ph(const Operator& opright, Operator& opout) const;
  void comm222_pp_hh_221(const Operator& opright, Operator& opout) const;


};

// Non member function, multiply by scalar from left side
Operator operator*(const double lhs, const Operator& rhs);

#endif

