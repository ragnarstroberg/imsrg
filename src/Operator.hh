
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
   ModelSpace * modelspace;
   arma::mat TBME;  // matrix of the two body matrix elements, where the indices label 2-body kets.
   int Number_pq_configs;  // Number of pq configs that participate in this channel

   // Constructors
   TwoBodyChannel();
   TwoBodyChannel(int j, int p, int t, Operator *op);
   TwoBodyChannel(const TwoBodyChannel&);

   //Overloaded operators
   TwoBodyChannel& operator=(const TwoBodyChannel&);
   TwoBodyChannel& operator+=(const TwoBodyChannel&);
   TwoBodyChannel& operator-=(const TwoBodyChannel&);
   TwoBodyChannel& operator+(const TwoBodyChannel&);
   TwoBodyChannel& operator-(const TwoBodyChannel&);

   //Methods
   int Getpqconfig(int a, int b){return pqindex[modelspace->GetKetIndex(a,b)];};  // do a lookup in pqindex
   int Getpqconfig(int ketindex){return pqindex[ketindex];};  // do a lookup in pqindex
   // This pqindex crap needs to be better written. It's confusing me, and I wrote it.
   // Rewrite it so GetTBME(bra,ket) does all this automatically.
   // Think of maybe using a std::map 
   float GetTBME(int bra, int ket); // Use the modelspace index of the bra and ket
   float GetTBME(Ket *bra, Ket *ket); // Use the bra and ket objects directly
   void SetTBME(int bra, int ket, float tbme);
   void SetTBME(Ket *bra, Ket *ket, float tbme);
   int GetKetIndex(int i){ return KetList[i];}; // i goes from 1 to Number_pq_configs. the output is the ket index in the modelspace numbering scheme
   Ket * GetKet(int i){ return modelspace->GetKet(KetList[i]);};

 private:
   //Fields
   vector<int> pqindex;  // lists index of a pq config in the TBME matrix, or else -1
   vector<int> KetList; // A list of kets in this channel, suitable for looping over using GetKetIndex
   bool hermitian;
   bool antihermitian;
   //Methods
   bool CheckChannel_pq(int,int);  // check if |pq> participates in this channel
   
};



class Operator
{
 public:
  //Fields
  float ZeroBody;
  arma::mat OneBody;
  TwoBodyChannel TwoBody[JMAX][2][3];

  //Constructors
  Operator();
  Operator(ModelSpace*);
  //Overloaded operators
  Operator(const Operator&);
  Operator operator=(const Operator&);

  //Methods
  float GetOneBody(int i,int j){return OneBody[i,j];};
  void SetOneBody(int i, int j, float val){ OneBody[i,j] = val;};
  Operator Commutator(const Operator& opleft, const Operator& opright);
  Operator BCH_Product(Operator&, Operator&); // not yet implemented
  Operator BCH_Transform(Operator&, Operator&); // not yet implemented
  TwoBodyChannel* GetTwoBodyChannel(int j, int p, int t){return &(TwoBody[j][p][(t+1)]);};
  void PrintOut();
  void UpdateCrossCoupled(); // Not implemented, because I don't fully understand it.
  ModelSpace * GetModelSpace(){return modelspace;};
  void GetSPEFromModelSpace();
  void SetHermitian() {hermitian=true;antihermitian=false;};
  void SetAntiHermitian() {antihermitian=true;hermitian=false;};
  void SetNonHermitian() {antihermitian=false;hermitian=false;};
  bool IsHermitian(){return hermitian;};
  bool IsAntiHermitian(){return antihermitian;};
  void EraseZeroBody(){ZeroBody = 0;}; // set zero-body term to zero
  void EraseOneBody(){OneBody.zeros();}; // set all one-body terms to zero
  void EraseTwoBody(); // set all two-body terms to zero
  Operator DoNormalOrdering(); // Do normal ordering -- not implemented yet
  void WriteOneBody(const char* filename);
  void WriteTwoBody(const char* filename);

 private:
  //Fields
  ModelSpace * modelspace;
  bool hermitian;
  bool antihermitian;
  int TwoBodyJmax;
  //Methods
  float comm110(arma::mat, arma::mat);
//  float comm210(TwoBodyChannel,arma::mat);
  float comm220(TwoBodyChannel,TwoBodyChannel);
  arma::mat comm111(arma::mat, arma::mat);
  arma::mat comm211(TwoBodyChannel,arma::mat);
  arma::mat comm221(TwoBodyChannel,TwoBodyChannel);
  TwoBodyChannel comm112(arma::mat, arma::mat);
  TwoBodyChannel comm212(TwoBodyChannel,arma::mat);
  TwoBodyChannel comm222(TwoBodyChannel,TwoBodyChannel);
  arma::mat P_hole_onebody; // Projector onto one-body hole states, ie n
  arma::mat P_particle_onebody; // Projector onto one-body particle states, ie nbar
  


};



#endif

