#ifndef ModelSpace_h
#define ModelSpace_h 1

//#define ARMA_NO_DEBUG
#include <vector>
#include <unordered_map>
#include <map>
#include <armadillo>
#include "IMSRGProfiler.hh"
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif
#define OCC_CUT 1e-6


using namespace std;

typedef unsigned long long int index_t;

class ModelSpace; //forward declaration so Ket can use ModelSpace

//struct Orbit
class Orbit
{
 public:
   // Fields
   int n;
   int l;
   int j2;
   int tz2;
   double occ; // particle=0, hole=1
//   int ph; // particle=0, hole=1
//   int io; // inside=0, outside=1
   int cvq; // core=0, valence=1, qspace=2
   int index;

   //Constructors
   ~Orbit();
   Orbit();
   Orbit(int n ,int l, int j, int t, double occ, int cvq, int index);
//   Orbit(int n ,int l, int j, int t, int ph, int cvq, int index);
   Orbit(const Orbit&);
//   void swap(Orbit&) throw();
//   Orbit& operator=( const Orbit& );
   // Methods
};




class Ket  //  | pq >
{

 public:
   // Constructors
   ~Ket();
   Ket();
   Ket(Orbit& op, Orbit& oq);
   // Methods
   int Phase(int J);
   int delta_pq(){return dpq;};

   // Fields
   Orbit* op;
   Orbit* oq;
   int p;
   int q;
 private:
   // Fields
   int dpq;
   int phase_prefactor;

};




class TwoBodyChannel
{
 public:
   //Fields
   int J;
   int parity;
   int Tz;

   // Constructors
   ~TwoBodyChannel();
   TwoBodyChannel();
   TwoBodyChannel(int j, int p, int t, ModelSpace* ms);
   TwoBodyChannel(int N, ModelSpace* ms);
   void Initialize(int N, ModelSpace* ms);


   //Methods
   int GetNumberKets() const {return NumberKets;};
   int GetLocalIndex(int ketindex) const { return KetMap[ketindex];}; // modelspace ket index => local ket index
   int GetLocalIndex(int p, int q) const ;
   int GetKetIndex(int i) const { return KetList[i];}; // local ket index => modelspace ket index
   Ket& GetKet(int i)  ; // get pointer to ket using local index

   arma::uvec KetIndex_pp ;
   arma::uvec KetIndex_hh ;
   arma::uvec KetIndex_ph ;
   arma::uvec KetIndex_cc ;
   arma::uvec KetIndex_vc ;
   arma::uvec KetIndex_qc ;
   arma::uvec KetIndex_vv ;
   arma::uvec KetIndex_qv ;
   arma::uvec KetIndex_qq ;
   arma::vec  Ket_occ_hh;
   arma::vec  Ket_unocc_hh;
   arma::vec  Ket_occ_ph;
   arma::vec  Ket_unocc_ph;



   arma::uvec GetKetIndexFromList(vector<index_t>& vec_in);
   arma::uvec& GetKetIndex_pp();
   arma::uvec& GetKetIndex_hh();
   arma::uvec& GetKetIndex_ph();
   arma::uvec& GetKetIndex_cc();
   arma::uvec& GetKetIndex_vc();
   arma::uvec& GetKetIndex_qc();
   arma::uvec& GetKetIndex_vv();
   arma::uvec& GetKetIndex_qv();
   arma::uvec& GetKetIndex_qq();


// private:
   //Fields
   ModelSpace * modelspace;
   int NumberKets;  // Number of pq configs that participate in this channel
   vector<int> KetList; // eg [2, 4, 7, ...] Used for looping over all the kets in the channel
   vector<int> KetMap;  // eg [ -1, -1, 0, -1, 1, -1, -1, 2 ...] Used for asking what is the local index of this ket. -1 means the ket doesn't participate in this channel
   //Methods
   virtual bool CheckChannel_ket(Orbit* op, Orbit* oq) const;  // check if |pq> participates in this channel
   bool CheckChannel_ket(Ket &ket) const {return CheckChannel_ket(ket.op,ket.oq);};  // check if |pq> participates in this channel
   
};



class TwoBodyChannel_CC : public TwoBodyChannel
{
  public:
   ~TwoBodyChannel_CC();
   TwoBodyChannel_CC();
   TwoBodyChannel_CC(int j, int p, int t, ModelSpace* ms);
   TwoBodyChannel_CC(int N, ModelSpace* ms);
   bool CheckChannel_ket(Orbit* op, Orbit* oq) const;  // check if |pq> participates in this channel
//   bool CheckChannel_ket(int p, int q) const;  // check if |pq> participates in this channel
};




class ModelSpace
{

 public:


   // Constructors
   ~ModelSpace();
   ModelSpace();
   ModelSpace(const ModelSpace&); // copy constructor
   ModelSpace( ModelSpace&&); // move constructor
   ModelSpace(int emax, vector<string> hole_list, vector<string> valence_list);
   ModelSpace(int emax, vector<string> hole_list, vector<string> core_list, vector<string> valence_list);
   ModelSpace(int emax, string reference, string valence);
   ModelSpace(int emax, string reference);

   // Overloaded operators
   ModelSpace operator=(const ModelSpace&); 
   ModelSpace operator=(ModelSpace&&); 

   // Methods

   void Init(int emax, string reference, string valence);
   void Init(int emax, map<index_t,double> hole_list, string valence);
   void Init(int emax, map<index_t,double> hole_list, vector<index_t> core_list, vector<index_t> valence_list);
   void Init(int emax, vector<string> hole_list, vector<string> core_list, vector<string> valence_list);
   void Init_occ_from_file(int emax, string valence, string occ_file);

//   vector<index_t> GetOrbitsAZ(int A, int Z);
   map<index_t,double> GetOrbitsAZ(int A, int Z);
   void GetAZfromString(string str, int& A, int& Z);
   vector<index_t> String2Index( vector<string> vs );
   string Index2String(index_t ind);
   void Get0hwSpace(int Aref, int Zref, vector<index_t>& core_list, vector<index_t>& valence_list);
   void ParseCommaSeparatedValenceSpace(string valence, vector<index_t>& core_list, vector<index_t>& valence_list);

   void SetupKets();
   void AddOrbit(Orbit orb);
   void AddOrbit(int n, int l, int j2, int tz2, double occ, int io);
   // Setter/Getters
   Orbit& GetOrbit(int i) {return (Orbit&) Orbits[i];}; 
//   Orbit& GetOrbit(int i) const {return (Orbit&) Orbits[i];}; 
   Ket& GetKet(int i) const {return (Ket&) Kets[i];};
   Ket& GetKet(int p, int q) const {return (Ket&) Kets[Index2(p,q)];};
   int GetOrbitIndex(int n, int l, int j2, int tz2) const {return Index1(n,l,j2,tz2);};
   int GetKetIndex(int p, int q) const {return Index2(p,q);}; // convention is p<=q
   int GetKetIndex(Ket * ket) const {return Index2(ket->p,ket->q);}; // convention is p<=q
   int GetNumberOrbits() const {return norbits;};
   int GetNumberKets() const {return Kets.size();};
   void SetHbarOmega(double hw) {hbar_omega = hw;};
   void SetTargetMass(int A) {target_mass = A;};
   void SetTargetZ(int Z) {target_Z = Z;};
   double GetHbarOmega() const {return hbar_omega;};
   int GetTargetMass() const {return target_mass;};
   int GetTargetZ() const {return target_Z;};
   int GetAref() const {return Aref;};
   int GetZref() const {return Zref;};
   int GetNumberTwoBodyChannels() const {return TwoBodyChannels.size();};
   TwoBodyChannel& GetTwoBodyChannel(int ch) const {return (TwoBodyChannel&) TwoBodyChannels[ch];};
   TwoBodyChannel_CC& GetTwoBodyChannel_CC(int ch) const {return (TwoBodyChannel_CC&) TwoBodyChannels_CC[ch];};
   inline int GetTwoBodyJmax() const {return TwoBodyJmax;};
   inline int GetThreeBodyJmax() const {return ThreeBodyJmax;};
   void SetReference(vector<index_t>);
   void SetReference(map<index_t,double>);
   void SetReference(string);

   int GetEmax(){return Emax;};
   int GetE2max(){return E2max;};
   int GetE3max(){return E3max;};
   int GetLmax2(){return Lmax2;};
   int GetLmax3(){return Lmax3;};
   void SetEmax(int e){Emax=e;};
   void SetE2max(int e){E2max=e;};
   void SetE3max(int e){E3max=e;};
   void SetLmax2(int l){Lmax2=l;};
   void SetLmax3(int l){Lmax3=l;};

   double GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3);
   double GetNineJ(double j1, double j2, double j3, double j4, double j5, double j6, double j7, double j8, double j9);
   double GetMoshinsky( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L); // Inconsistent notation. Not ideal.
   bool SixJ_is_empty(){ return SixJList.empty(); };

   int GetOrbitIndex(string);
   int GetTwoBodyChannelIndex(int j, int p, int t);
   inline int phase(int x) {return (x%2)==0 ? 1 : -1;};
   inline int phase(double x) {return phase(int(x));};

   inline int Index1(int n, int l, int j2, int tz2) const {return(2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2 ;};
//   inline int Index2(int p, int q) const {return q*(q+1)/2 + p;};
   inline int Index2(int p, int q) const {return p*(2*norbits-1-p)/2 + q;};

   void PreCalculateMoshinsky();
   void ClearVectors();
   void CalculatePandyaLookup(int rank_J, int rank_T, int parity); // construct a lookup table for more efficient pandya transformation
//   map<array<int,2>,vector<array<int,2>>>& GetPandyaLookup(int rank_J, int rank_T, int parity);
   map<array<int,2>,array<vector<int>,2>>& GetPandyaLookup(int rank_J, int rank_T, int parity);


   // Data members
   vector<index_t> holes;           // in the reference Slater determinant
//   map<index_t,double> holes;           // in the reference Slater determinant
   vector<index_t> particles;       // above the reference Slater determinant
   vector<index_t> core;            // core for decoupling
   vector<index_t> valence;         // valence space for decoupling
   vector<index_t> qspace;          // above the valence space for decoupling
   vector<index_t> proton_orbits;
   vector<index_t> neutron_orbits;

   vector<index_t> KetIndex_pp; 
   vector<index_t> KetIndex_ph;
   vector<index_t> KetIndex_hh;
   vector<index_t> KetIndex_cc; 
   vector<index_t> KetIndex_vc;
   vector<index_t> KetIndex_qc;
   vector<index_t> KetIndex_vv;
   vector<index_t> KetIndex_qv;
   vector<index_t> KetIndex_qq;
   vector<double> Ket_occ_hh;
   vector<double> Ket_unocc_hh;
   vector<double> Ket_occ_ph;
   vector<double> Ket_unocc_ph;

//   array< array< vector<index_t>, 2>,3> MonopoleKets; //List of kets of a given Tz,parity
   array< array< unordered_map<index_t,index_t>, 2>,3> MonopoleKets; //List of kets of a given Tz,parity

   int Emax;
   int E2max;
   int E3max;
   int Lmax2;
   int Lmax3;
   int OneBodyJmax;
   int TwoBodyJmax;
   int ThreeBodyJmax;
   map<array<int,3>,vector<index_t> > OneBodyChannels;

   vector<unsigned int> SortedTwoBodyChannels;
   vector<unsigned int> SortedTwoBodyChannels_CC;

   static map<string,vector<string>> ValenceSpaces;
//   map< array<int,3>, map< array<int,2>,vector<array<int,2>> > > PandyaLookup;


// private:
   // Fields
   int norbits;
   double hbar_omega;
   int target_mass;
   int target_Z;
   int Aref;
   int Zref;
   int nTwoBodyChannels;
   vector<Orbit> Orbits;
   vector<Ket> Kets;
   vector<TwoBodyChannel> TwoBodyChannels;
   vector<TwoBodyChannel_CC> TwoBodyChannels_CC;
   map< array<int,3>, map< array<int,2>,array<vector<int>,2> > > PandyaLookup;
   bool moshinsky_has_been_precalculated;
   IMSRGProfiler profiler;
//   map<long int,double> SixJList;

   static unordered_map<unsigned long int,double> SixJList;
   static unordered_map<long long unsigned int,double> NineJList;
   static unordered_map<long long unsigned int,double> MoshList;

};




#endif

