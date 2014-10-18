
#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include <armadillo>

#define JMAX 30

using namespace std;

class Operator
{
 public:
  //Fields
  float ZeroBody;
  arma::mat OneBody;
  arma::mat TwoBody[JMAX*2*3];

  //Constructors
  Operator();
  Operator(ModelSpace*);
  Operator( const Operator& rhs){Copy(rhs);};

  //Overloaded operators
  Operator& operator=( const Operator& rhs) {Copy(rhs); return *this;};
  Operator& operator+=( const Operator& rhs);
  Operator operator+( const Operator& rhs);
  Operator& operator*=( const double rhs);
  Operator operator*( const double rhs);

  //Methods
  // One body setter/getters
  float GetOneBody(int i,int j) {return OneBody[i,j];};
  void SetOneBody(int i, int j, float val) { OneBody[i,j] = val;};

  //TwoBody setter/getters
  double GetTBME(int ch, int a, int b, int c, int d) const;
  void   SetTBME(int ch, int a, int b, int c, int d, double tbme);
  double GetTBME(int ch, Ket *bra, Ket *ket) const;
  void   SetTBME(int ch, Ket *bra, Ket* ket, double tbme);
  double GetTBME(int j, int p, int t, Ket* bra, Ket* ket) const;
  void   SetTBME(int j, int p, int t, Ket* bra, Ket* ket, double tbme);
  double GetTBME(int j, int p, int t, int a, int b, int c, int d) const;
  void   SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme);

  double GetTBMEmonopole(int a, int b, int c, int d);
  double GetTBMEmonopole(Ket * bra, Ket * ket);

  // The actually interesting methods
  Operator Commutator(Operator& opright);
  Operator BCH_Product( Operator& ); // not yet implemented
  Operator BCH_Transform( Operator& ); // not yet implemented
  Operator DoNormalOrdering(); // Do normal ordering -- not implemented yet

  double Norm();

  ModelSpace * GetModelSpace() const {return modelspace;};

  void EraseZeroBody(){ZeroBody = 0;}; // set zero-body term to zero
  void EraseOneBody(){OneBody.zeros();}; // set all one-body terms to zero
  void EraseTwoBody(); // set all two-body terms to zero

  void SetHermitian() {hermitian=true;antihermitian=false;};
  void SetAntiHermitian() {antihermitian=true;hermitian=false;};
  void SetNonHermitian() {antihermitian=false;hermitian=false;};
  bool IsHermitian()const {return hermitian;};
  bool IsAntiHermitian()const {return antihermitian;};

  void PrintTwoBody() ;
  void PrintOut() ;

 //private:
  //Fields
  ModelSpace * modelspace;
  bool hermitian;
  bool antihermitian;
  int TwoBodyJmax;
  int nChannels;
  //Methods
  void Copy(const Operator& rhs);
  
  void UpdateCrossCoupled(); // Not implemented, because I don't fully understand it.
  double comm110(Operator& opright);
  double comm220(Operator& opright);
  arma::mat comm111(Operator& opright);
  arma::mat comm121(Operator& opright);
  arma::mat comm211(Operator& opright);
  arma::mat comm221(Operator& opright);
  arma::mat comm122(Operator& opright, int ch);
  arma::mat comm212(Operator& opright, int ch);
  arma::mat comm222_ph(Operator& opright, int ch);
  arma::mat comm222_pphh(Operator& opright, int ch);

  arma::mat P_hole_onebody; // Projector onto one-body hole states, ie the number operator n
  arma::mat P_particle_onebody; // Projector onto one-body particle states, ie nbar
  


};



#endif

