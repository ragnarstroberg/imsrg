#ifndef ModelSpace_h
#define ModelSpace_h 1

//#define ARMA_NO_DEBUG
#include <vector>
#include <unordered_map>
#include <map>
#include <armadillo>
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif


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
   int ph; // particle=0, hole=1
   int io; // inside=0, outside=1
   int index;

   //Constructors
   ~Orbit();
   Orbit();
   Orbit(int n ,int l, int j, int t, int ph, int io, int index);
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
   Ket& GetKet(int i) const ; // get pointer to ket using local index

   arma::uvec KetIndex_pp ;
   arma::uvec KetIndex_hh ;
   arma::uvec KetIndex_ph ;
   arma::uvec KetIndex_vv ;
   arma::uvec KetIndex_c_c;
   arma::uvec KetIndex_q_q;
   arma::uvec KetIndex_q_c ;
   arma::uvec KetIndex_v_c ;
   arma::uvec KetIndex_v_q;


   arma::uvec GetKetIndexFromList(vector<index_t>& vec_in);
   arma::uvec& GetKetIndex_pp();
   arma::uvec& GetKetIndex_hh();
   arma::uvec& GetKetIndex_ph();
   arma::uvec& GetKetIndex_vv();
   arma::uvec& GetKetIndex_c_c();  // cc
   arma::uvec& GetKetIndex_q_q();  //qq
   arma::uvec& GetKetIndex_q_c(); // qc
   arma::uvec& GetKetIndex_v_c(); // vc
   arma::uvec& GetKetIndex_v_q(); // qv

// private:
   //Fields
   ModelSpace * modelspace;
   int NumberKets;  // Number of pq configs that participate in this channel
   vector<int> KetList; // eg [2, 4, 7, ...] Used for looping over all the kets in the channel
   vector<int> KetMap;  // eg [ -1, -1, 0, -1, 1, -1, -1, 2 ...] Used for asking what is the local index of this ket. -1 means the ket doesn't participate in this channel
   //Methods
//   virtual bool CheckChannel_ket(int p, int q) const;  // check if |pq> participates in this channel
//   bool CheckChannel_ket(Ket &ket) const {return CheckChannel_ket(ket.p,ket.q);};  // check if |pq> participates in this channel
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
   ModelSpace(int Nmax, vector<string> hole_list, vector<string> inside_list);
   ModelSpace(int Nmax, string);
   ModelSpace(int Nmax, int A, int Z);

   // Overloaded operators
   ModelSpace operator=(const ModelSpace&); 
   ModelSpace operator=(ModelSpace&&); 

   void Init(int Nmax, vector<string> hole_list, vector<string> inside_list);
   void Init_AZ(int nmax, int A, int Z);

   // Common model spaces
//   void Init_He4(int nmax);
//   void Init_O16(int nmax);
//   void Init_Ca40(int nmax);
   void Init_PShell(int nmax);
   void Init_SDShell(int nmax);
   void Init_PSDShell(int nmax);
   void Init_O16PSDShell(int nmax);
   void Init_FPShell(int nmax);
   void Init_SDFPShell(int nmax);
   void Init_SD3F7P3Shell(int nmax);
   void Init_FPG9Shell(int nmax);


   // Methods
   void SetupKets();
   void AddOrbit(Orbit orb);
   void AddOrbit(int n, int l, int j2, int tz2, int ph, int io);
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
   int GetNumberTwoBodyChannels() const {return TwoBodyChannels.size();};
   TwoBodyChannel& GetTwoBodyChannel(int ch) const {return (TwoBodyChannel&) TwoBodyChannels[ch];};
   TwoBodyChannel_CC& GetTwoBodyChannel_CC(int ch) const {return (TwoBodyChannel_CC&) TwoBodyChannels_CC[ch];};
   inline int GetTwoBodyJmax() const {return TwoBodyJmax;};
   inline int GetThreeBodyJmax() const {return ThreeBodyJmax;};

   int GetNmax(){return Nmax;};
   int GetN2max(){return N2max;};
   int GetN3max(){return N3max;};
   void SetNmax(int n){Nmax=n;};
   void SetN2max(int n){N2max=n;};
   void SetN3max(int n){N3max=n;};

   double GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3);
   double GetNineJ(double j1, double j2, double j3, double j4, double j5, double j6, double j7, double j8, double j9);
   double GetMoshinsky( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L); // Inconsistent notation. Not ideal.
   bool SixJ_is_empty(){ return SixJList.empty(); };

   int GetOrbitIndex(string);
   int GetTwoBodyChannelIndex(int j, int p, int t);
   inline int phase(int x) {return (x%2)==0 ? 1 : -1;};
   inline int phase(double x) {return phase(int(x));};

   inline int Index1(int n, int l, int j2, int tz2) const {return(2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2 ;};
   inline int Index2(int p, int q) const {return q*(q+1)/2 + p;};

   void PreCalculateMoshinsky();


   // Data members
   vector<index_t> holes;
   vector<index_t> particles;
   vector<index_t> valence;
   vector<index_t> qspace;
   vector<index_t> hole_qspace;
   vector<index_t> particle_qspace;
   vector<index_t> proton_orbits;
   vector<index_t> neutron_orbits;

   vector<index_t> KetIndex_pp; 
   vector<index_t> KetIndex_ph;
   vector<index_t> KetIndex_hh;
   vector<index_t> KetIndex_vv;
   vector<index_t> KetIndex_c_c; 
   vector<index_t> KetIndex_q_q;
   vector<index_t> KetIndex_q_c;
   vector<index_t> KetIndex_v_c;
   vector<index_t> KetIndex_v_q;

//   array< array< vector<index_t>, 2>,3> MonopoleKets; //List of kets of a given Tz,parity
   array< array< unordered_map<index_t,index_t>, 2>,3> MonopoleKets; //List of kets of a given Tz,parity

   int Nmax;
   int N2max;
   int N3max;
   int OneBodyJmax;
   int TwoBodyJmax;
   int ThreeBodyJmax;
   map<array<int,3>,vector<index_t> > OneBodyChannels;

   vector<unsigned int> SortedTwoBodyChannels;
   vector<unsigned int> SortedTwoBodyChannels_CC;



// private:
   // Fields
   int norbits;
   double hbar_omega;
   int target_mass;
   int target_Z;
   int nTwoBodyChannels;
   vector<Orbit> Orbits;
   vector<Ket> Kets;
   vector<TwoBodyChannel> TwoBodyChannels;
   vector<TwoBodyChannel_CC> TwoBodyChannels_CC;
//   map<long int,double> SixJList;

   static unordered_map<unsigned long int,double> SixJList;
   static unordered_map<long long unsigned int,double> NineJList;
   static unordered_map<long long unsigned int,double> MoshList;

};




#endif

