#ifndef ModelSpace_h
#define ModelSpace_h 1

#include <vector>
#include <map>
#include <armadillo>

#define JMAX 30

class Orbit;
class Ket;
class ModelSpace;

using namespace std;



class Orbit
{
 public:
   // Fields
   int n;
   int l;
   int j2;
   int tz2;
   int hvq; // hole=0, valence=1, qspace=2
   float spe;

   //Constructors
   Orbit();
   Orbit(int n ,int l, int j, int t, int hvq ,float spe);
   Orbit(const Orbit&);
   // Methods
   void Set(int n, int l, int j2, int tz2, int hvq, float spe);
};




class Ket  //  | pq >
{
 public:
   // Fields
   int p;
   int q;
   int parity;
   int Tz;
   int Jmin;
   int Jmax;
   int Jstep;
   // Constructors
   Ket(){};
   Ket(ModelSpace * ms, int p, int q);
   // Methods
   int Phase(int J);
   int delta_pq(){return p==q ? 1 : 0;};

 private:
   // Fields
   ModelSpace * ms;
   int phase_prefactor;

};




class TwoBodyChannel
{
 public:
   //Fields
   int J;
   int parity;
   int Tz;
   arma::mat Proj_pp; // Projector onto pp kets
   arma::mat Proj_hh; // Projector onto hh kets
   arma::mat Proj_ph_cc; // Projector onto ph kets, for use with cross-coupled TBMEs

   // Constructors
   TwoBodyChannel();
   TwoBodyChannel(int j, int p, int t, ModelSpace* ms);
   TwoBodyChannel(int N, ModelSpace* ms);
   TwoBodyChannel(const TwoBodyChannel& rhs) {Copy(rhs);};
   void Initialize(int N, ModelSpace* ms);

   //Overloaded operators
   TwoBodyChannel& operator=(const TwoBodyChannel& rhs) {Copy(rhs); return *this;};
   TwoBodyChannel& operator+=(const TwoBodyChannel&);
   TwoBodyChannel operator+(const TwoBodyChannel& rhs){return TwoBodyChannel(*this) += rhs;};
   TwoBodyChannel& operator-=(const TwoBodyChannel&);
   TwoBodyChannel operator-(const TwoBodyChannel& rhs){return TwoBodyChannel(*this) -= rhs;};

   //Methods
   int GetNumberKets() const {return NumberKets;};
   int GetLocalIndex(int ketindex) const { return KetMap[ketindex];}; // modelspace ket index => local ket index
   int GetLocalIndex(int p, int q) const ;
   int GetKetIndex(int i) const { return KetList[i];}; // local ket index => modelspace ket index
   Ket * GetKet(int i) const ; // get pointer to ket using local index
   vector<int> KetIndex_pp;
   vector<int> KetIndex_ph;
   vector<int> KetIndex_hh;


// private:
   //Fields
   ModelSpace * modelspace;
   int NumberKets;  // Number of pq configs that participate in this channel
   vector<int> KetList; // eg [2, 4, 7, ...] Used for looping over all the kets in the channel
   vector<int> KetMap;  // eg [ -1, -1, 0, -1, 1, -1, -1, 2 ...] Used for asking what is the local index of this ket. -1 means the ket doesn't participate in this channel
   //Methods
   virtual bool CheckChannel_ket(int p, int q) const;  // check if |pq> participates in this channel
   bool CheckChannel_ket(Ket *ket) const {return CheckChannel_ket(ket->p,ket->q);};  // check if |pq> participates in this channel
   void Copy(const TwoBodyChannel &);
   
};

class TwoBodyChannel_CC : public TwoBodyChannel
{
  public:
   TwoBodyChannel_CC();
   TwoBodyChannel_CC(int j, int p, int t, ModelSpace* ms);
   TwoBodyChannel_CC(int N, ModelSpace* ms);
   TwoBodyChannel_CC(const TwoBodyChannel& rhs) {Copy(rhs);};
//   void Initialize(int N, ModelSpace* ms);
   bool CheckChannel_ket(int p, int q) const;  // check if |pq> participates in this channel
};



class ModelSpace
{

 public:
   // Fields
   vector<int> holes;
   vector<int> particles;
   vector<int> valence;
   vector<int> qspace;

   // Constructors
   ModelSpace();
   ModelSpace(const ModelSpace&);

   // Overloaded operators
   ModelSpace operator=(const ModelSpace&); 

   // Methods
   void SetupKets();
   void AddOrbit(Orbit orb);
   // Setter/Getters
   Orbit* GetOrbit(int i) const {return Orbits[i];}; 
   Ket* GetKet(int i) const {return Kets[i];};
   Ket* GetKet(int p, int q) const {return Kets[Index2(p,q)];};
   int GetKetIndex(int p, int q) const {return Index2(p,q);}; // convention is p<=q
   int GetKetIndex(Ket * ket) const {return Index2(ket->p,ket->q);}; // convention is p<=q
   int GetNumberOrbits() const {return norbits;};
   int GetNumberKets() const {return Kets.size();};
   void SetHbarOmega(float hw) {hbar_omega = hw;};
   void SetTargetMass(int A) {target_mass = A;};
   float GetHbarOmega() const {return hbar_omega;};
   int GetTargetMass() const {return target_mass;};
   int GetNumberTwoBodyChannels() const {return TwoBodyChannels.size();};
   //TwoBodyChannel& GetTwoBodyChannel(int ch) {return TwoBodyChannels[ch];};
   TwoBodyChannel& GetTwoBodyChannel(int ch) const {return (TwoBodyChannel&) TwoBodyChannels[ch];};
   TwoBodyChannel_CC& GetTwoBodyChannel_CC(int ch) const {return (TwoBodyChannel_CC&) TwoBodyChannels_CC[ch];};
   TwoBodyChannel* GetTwoBodyChannel_ptr(int ch) {return (&TwoBodyChannels[ch]);};
   int GetTwoBodyJmax() const {return TwoBodyJmax;};
   //double GetSixJ(int j1, int j2, int j3, int J1, int J2, int J3);
   double GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3);
   double GetNineJ(double j1, double j2, double j3, double j4, double j5, double j6, double j7, double j8, double j9);
   int GetTwoBodyChannelIndex(int j, int p, int t);
   int phase(double x) {return 1-2*(int(x)%2);};

   int Index1(int n, int l, int j2, int tz2) const {return(2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2 ;};
   //int Index2(int p, int q) const {return q*norbits + p;};
   inline int Index2(int p, int q) const {return q*(q+1)/2 + p;};


 private:
   // Fields
   int norbits;
   int maxj;
   int TwoBodyJmax;
   float hbar_omega;
   int target_mass;
   std::vector<Orbit*> Orbits;
   std::vector<Ket*> Kets;
   std::vector<TwoBodyChannel> TwoBodyChannels;
   std::vector<TwoBodyChannel_CC> TwoBodyChannels_CC;
   int nTwoBodyChannels;
   std::map<long int,double> SixJList;

};




#endif

