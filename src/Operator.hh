
#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include <armadillo>

#define JMAX 30

using namespace std;

class Operator;

class TwoBodyChannel
{
 public:
   //Fields
   int J;
   int parity;
   int Tz;
   arma::mat TBME;  // matrix of the two body matrix elements, where the indices label 2-body kets.

   // Constructors
   TwoBodyChannel();
   TwoBodyChannel(int j, int p, int t, Operator *op);
   TwoBodyChannel(int N, Operator *op);
   TwoBodyChannel(const TwoBodyChannel& rhs) {Copy(rhs);};

   //Overloaded operators
   TwoBodyChannel& operator=(const TwoBodyChannel& rhs) {Copy(rhs); return *this;};
   TwoBodyChannel& operator+=(const TwoBodyChannel&);
   TwoBodyChannel operator+(const TwoBodyChannel& rhs){return TwoBodyChannel(*this) += rhs;};
   TwoBodyChannel& operator-=(const TwoBodyChannel&);
   TwoBodyChannel operator-(const TwoBodyChannel& rhs){return TwoBodyChannel(*this) -= rhs;};

   //Methods
   float GetTBME(int bra, int ket) const ; // Use the modelspace index of the bra and ket
   float GetTBME(Ket *bra, Ket *ket) const ; // Use the bra and ket objects directly
   int GetNumberKets() const {return NumberKets;};
   void SetTBME(int bra, int ket, float tbme);
   void SetTBME(Ket *bra, Ket *ket, float tbme);
   void SetTBME(const arma::mat& tbme){ TBME = tbme;};

   int GetLocalIndex(int ketindex) const { return KetMap[ketindex];}; // modelspace ket index => local ket index
   int GetLocalIndex(int p, int q) const { return KetMap[modelspace->GetKetIndex(p,q)];};
   int GetKetIndex(int i) const { return KetList[i];}; // local ket index => modelspace ket index
   Ket * GetKet(int i) const { return modelspace->GetKet(KetList[i]);};

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



class Operator
{
 public:
  //Fields
  float ZeroBody;
  arma::mat OneBody;
  //TwoBodyChannel TwoBody[JMAX][2][3]; // Might be advantageous to flatten this down to a 1d array to simplify looping
  TwoBodyChannel TwoBody[JMAX*2*3]; // Might be advantageous to flatten this down to a 1d array to simplify looping

  //Constructors
  Operator();
  Operator(ModelSpace*);
  Operator( const Operator& rhs){Copy(rhs);};

  //Overloaded operators
  Operator operator=( const Operator& rhs) {Copy(rhs); return *this;};

  //Methods
  float GetOneBody(int i,int j) {return OneBody[i,j];};
  void SetOneBody(int i, int j, float val) { OneBody[i,j] = val;};
  Operator Commutator(const Operator& opright);
  Operator BCH_Product(const Operator&, const Operator&); // not yet implemented
  Operator BCH_Transform(const Operator&, const Operator&); // not yet implemented
  TwoBodyChannel* GetTwoBodyChannel(int j, int p, int t) {return &(TwoBody[(t+1)*2*JMAX + p*JMAX + j]);};
  TwoBodyChannel* GetTwoBodyChannel(int N) {return &(TwoBody[N]);};
  arma::mat GetTBME(int N) const {return TwoBody[N].TBME;};
  arma::mat GetTBME(int j, int p, int t) const {return TwoBody[(t+1)*2*JMAX + p*JMAX + j].TBME;};
  void SetTwoBodyChannel(int N, const TwoBodyChannel& tbc){ TwoBody[N] = tbc;};
  void SetTwoBodyChannel(int j, int p, int t, const TwoBodyChannel& tbc) {TwoBody[(t+1)*2*JMAX + p*JMAX + j] = tbc;};
  int GetNumberTwoBodyChannels()const {return nChannels;};
  void PrintOut() ;
  void UpdateCrossCoupled(); // Not implemented, because I don't fully understand it.
  ModelSpace * GetModelSpace() const {return modelspace;};
  void SetHermitian() {hermitian=true;antihermitian=false;};
  void SetAntiHermitian() {antihermitian=true;hermitian=false;};
  void SetNonHermitian() {antihermitian=false;hermitian=false;};
  bool IsHermitian()const {return hermitian;};
  bool IsAntiHermitian()const {return antihermitian;};
  void EraseZeroBody(){ZeroBody = 0;}; // set zero-body term to zero
  void EraseOneBody(){OneBody.zeros();}; // set all one-body terms to zero
  void EraseTwoBody(); // set all two-body terms to zero
  Operator DoNormalOrdering(); // Do normal ordering -- not implemented yet

 private:
  //Fields
  ModelSpace * modelspace;
  bool hermitian;
  bool antihermitian;
  int TwoBodyJmax;
  int nChannels;
  //Methods
  void Copy(const Operator& rhs);
  float comm110(const arma::mat&, const arma::mat&);
  float comm220(const TwoBodyChannel&,const TwoBodyChannel&);
  arma::mat comm111(const arma::mat&, const arma::mat&);
  arma::mat comm211(const TwoBodyChannel&,const arma::mat&);
  arma::mat comm221(const TwoBodyChannel&, const TwoBodyChannel&);
  TwoBodyChannel comm112(const arma::mat&, const arma::mat&);
  TwoBodyChannel comm212(const TwoBodyChannel&, const arma::mat&);
  TwoBodyChannel comm222(const TwoBodyChannel&,const TwoBodyChannel&);
  arma::mat P_hole_onebody; // Projector onto one-body hole states, ie the number operator n
  arma::mat P_particle_onebody; // Projector onto one-body particle states, ie nbar
  


};



#endif

