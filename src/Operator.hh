
#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include <armadillo>

#define JMAX 30

using namespace std;

class Operator;
/*
class TwoBodyChannel
{
 public:
   //Fields
   int J;
   int parity;
   int Tz;
   arma::mat TBME;  // matrix of the two body matrix elements, where the indices label 2-body kets.
   arma::mat Proj_pp; // Projector onto pp kets
   arma::mat Proj_hh; // Projector onto hh kets

   // Constructors
   TwoBodyChannel();
   TwoBodyChannel(int j, int p, int t, Operator *op);
   TwoBodyChannel(int N, Operator *op);
   TwoBodyChannel(const TwoBodyChannel& rhs) {Copy(rhs);};
   void Initialize(int N, Operator *op);

   //Overloaded operators
   TwoBodyChannel& operator=(const TwoBodyChannel& rhs) {Copy(rhs); return *this;};
   TwoBodyChannel& operator+=(const TwoBodyChannel&);
   TwoBodyChannel operator+(const TwoBodyChannel& rhs){return TwoBodyChannel(*this) += rhs;};
   TwoBodyChannel& operator-=(const TwoBodyChannel&);
   TwoBodyChannel operator-(const TwoBodyChannel& rhs){return TwoBodyChannel(*this) -= rhs;};

   //Methods
   double GetTBME(int bra, int ket) const ; // Use the modelspace index of the bra and ket
   double GetTBME(Ket *bra, Ket *ket) const ; // Use the bra and ket objects directly
   double GetTBME(int a, int b, int c, int d);
   int GetNumberKets() const {return NumberKets;};
   void SetTBME(int bra, int ket, float tbme);
   void SetTBME(Ket *bra, Ket *ket, float tbme);
   void SetTBME(const arma::mat& tbme){ TBME = tbme;};

   int GetLocalIndex(int ketindex) const { return KetMap[ketindex];}; // modelspace ket index => local ket index
   int GetLocalIndex(int p, int q) const { return KetMap[modelspace->GetKetIndex(p,q)];};
   int GetKetIndex(int i) const { return KetList[i];}; // local ket index => modelspace ket index
   Ket * GetKet(int i) const { return modelspace->GetKet(KetList[i]);}; // get pointer to ket using local index


 private:
   //Fields
   ModelSpace * modelspace;
   vector<int> KetMap;  // eg [ -1, -1, 0, -1, 1, -1, -1, 2 ...] Used for asking what is the local (channel) index of this ket
   vector<int> KetList; // eg [2, 4, 7, ...] Used for looping over all the kets in the channel
   int NumberKets;  // Number of pq configs that participate in this channel
   bool hermitian;
   bool antihermitian;
   //Methods
   bool CheckChannel_ket(int p, int q) const;  // check if |pq> participates in this channel
   bool CheckChannel_ket(Ket *ket) const {return CheckChannel_ket(ket->p,ket->q);};  // check if |pq> participates in this channel
   void Copy(const TwoBodyChannel &);
   
};
*/


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
  Operator operator=( const Operator& rhs) {Copy(rhs); return *this;};
  Operator& operator+=( const Operator& rhs);
  Operator& operator+( const Operator& rhs);
  Operator& operator*=( const double rhs);
  Operator& operator*( const double rhs);

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

  // Access the whole Two body matrix for a given channel
//  arma::mat& GetTwoBody(int ch) const {return TwoBody[ch];};  // not so useful if TwoBody is public...
//  void       SetTwoBody(int ch, const arma::mat& tbme){ TwoBody[ch] = tbme;};
//  arma::mat& GetTwoBody(int j, int p, int t) const {return TwoBody[(t+1)*2*JMAX + p*JMAX + j];};
//  void       SetTwoBody(int j, int p, int t, const arma::mat& tb) {TwoBody[(t+1)*2*JMAX + p*JMAX + j] = tb;};

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

  void PrintOut() ;

 private:
  //Fields
  ModelSpace * modelspace;
  bool hermitian;
  bool antihermitian;
  int TwoBodyJmax;
  int nChannels;
  //Methods
  void Copy(const Operator& rhs);
/*
  double comm110(const arma::mat&, const arma::mat&);
  double comm220(const TwoBodyChannel&,const TwoBodyChannel&);
  arma::mat comm111(const arma::mat&, const arma::mat&);
  arma::mat comm211(const TwoBodyChannel&,const arma::mat&);
  arma::mat comm221(const TwoBodyChannel&, const TwoBodyChannel&);
  TwoBodyChannel comm112(const arma::mat&, const arma::mat&);
  TwoBodyChannel comm212(const TwoBodyChannel&, const arma::mat&);
  TwoBodyChannel comm222(const TwoBodyChannel&,const TwoBodyChannel&);
*/  
  void UpdateCrossCoupled(); // Not implemented, because I don't fully understand it.
  double comm110(Operator& opright);
  double comm220(Operator& opright);
  arma::mat comm111(Operator& opright);
  arma::mat comm211(Operator& opright);
  arma::mat comm221(Operator& opright);
  arma::mat comm212(Operator& opright, int ch);
  arma::mat comm222_ph(Operator& opright, int ch);
  arma::mat comm222_pphh(Operator& opright, int ch);

  arma::mat P_hole_onebody; // Projector onto one-body hole states, ie the number operator n
  arma::mat P_particle_onebody; // Projector onto one-body particle states, ie nbar
  


};



#endif

