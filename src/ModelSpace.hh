///////////////////////////////////////////////////////////////////////////////////
//    ModelSpace.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef ModelSpace_h
#define ModelSpace_h 1

//#define ARMA_NO_DEBUG
#include <vector>
#include <unordered_map>
#include <map>
#include <array>
#include <set>
#include <armadillo>
#include "IMSRGProfiler.hh"
//#include "TwoBodyChannel.hh"
//#ifndef SQRT2
//  #define SQRT2 1.4142135623730950488
//#endif
//#define OCC_CUT 1e-6


//using namespace std;

typedef unsigned long long int index_t;

class ModelSpace; //forward declaration so Ket can use ModelSpace

struct Orbit
{
   int n;
   int l;
   int j2;
   int tz2;
   double occ; // particle=0, hole=1
   double occ_nat; // occupation in the natural orbital basis (if we use natural orbitals) particle=0, hole=1
   int cvq; // core=0, valence=1, qspace=2
   int index;

   Orbit(){};
   Orbit(int n, int l, int j2, int tz2, double occ, int cvq, int index)
           : n(n), l(l), j2(j2), tz2(tz2),occ(occ),occ_nat(occ),cvq(cvq),index(index) {};

   Orbit(const Orbit& orb) : n(orb.n), l(orb.l), j2(orb.j2), tz2(orb.tz2),occ(orb.occ),occ_nat(orb.occ),cvq(orb.cvq),index(orb.index) {}
   bool operator==( const Orbit& rhs ) const { return  ( n==rhs.n and l==rhs.l and j2==rhs.j2 and tz2==rhs.tz2 ); };
};



struct Ket  //  | pq >
{
   // Fields
   Orbit* op;
   Orbit* oq;
   size_t p;
   size_t q;
   int dpq;
   int phase_prefactor;

   // Methods
   Ket(){};
    Ket(Orbit& op_in, Orbit& oq_in) : op(&op_in), oq(&oq_in), p(op_in.index), q(oq_in.index), dpq( op_in.index==oq_in.index ? 1 : 0)
    {
     phase_prefactor = ((op->j2+oq->j2)/2 + 1) % 2==0 ? 1 : -1;
    };

   int Phase(int J)  {  return phase_prefactor * (J%2==0 ? 1 : -1); };
   int delta_pq() const {return dpq;};

};


struct Ket3 // | p q r >
{
  // Fields
  size_t p;
  size_t q;
  size_t r;
  Orbit* op;
  Orbit* oq;
  Orbit* oR; // unfortunate that this convention collides with the "or" keyword.
  int Jpq;

  // Methods
  Ket3(){};
  Ket3(Orbit& op_in, Orbit& oq_in, Orbit& oR_in)
     : p(op_in.index), q(oq_in.index), r(oR_in.index), op(&op_in), oq(&oq_in), oR(&oR_in), Jpq(0) {};

  Ket3(Orbit& op_in, Orbit& oq_in, Orbit& oR_in, int jpq)
       : p(op_in.index), q(oq_in.index), r(oR_in.index), op(&op_in), oq(&oq_in), oR(&oR_in), Jpq(jpq) {};

  bool operator==(Ket3& rhs) { return (p==rhs.p and q==rhs.q and r==rhs.r and Jpq==rhs.Jpq);};
  bool operator!=(Ket3& rhs) { return not( (*this)== rhs ) ;};
};





// This is now given its own implementation file, TwoBodyChannel.cc
// but the header info stays here, because of the inextricable dependencies
// between TwoBodyChannel, ModelSpace, Orbit, and Ket. What a mess...
struct TwoBodyChannel
{
// public:
   //Fields
   int J;
   int parity;
   int Tz;

   ModelSpace * modelspace;
   int NumberKets;  // Number of pq configs that participate in this channel
   std::vector<int> KetList; // eg [2, 4, 7, ...] Used for looping over all the kets in the channel
   std::vector<int> KetMap;  // eg [ -1, -1, 0, -1, 1, -1, -1, 2 ...] Used for asking what is the local index of this ket. -1 means the ket doesn't participate in this channel

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


   // Constructors
   virtual ~TwoBodyChannel();
   TwoBodyChannel();
   TwoBodyChannel(int j, int p, int t, ModelSpace* ms);
   TwoBodyChannel(int ch, ModelSpace* ms);
   void Initialize(int ch, ModelSpace* ms);


   //Methods
   size_t GetNumberKets() const {return NumberKets;};
   size_t GetLocalIndex(int ketindex) const { return KetMap[ketindex];}; // modelspace ket index => local ket index
   size_t GetLocalIndex(int p, int q) const ;
   size_t GetKetIndex(int i) const { return KetList[i];}; // local ket index => modelspace ket index
   const Ket& GetKet(int i) const  ; // get pointer to ket using local index
   Ket& GetKet(int i)  ; // get pointer to ket using local index



   arma::uvec GetKetIndexFromList(std::vector<index_t>& vec_in);
   const arma::uvec& GetKetIndex_pp() const;
   const arma::uvec& GetKetIndex_hh() const;
   const arma::uvec& GetKetIndex_ph() const;
   const arma::uvec& GetKetIndex_cc() const;
   const arma::uvec& GetKetIndex_vc() const;
   const arma::uvec& GetKetIndex_qc() const;
   const arma::uvec& GetKetIndex_vv() const;
   const arma::uvec& GetKetIndex_qv() const;
   const arma::uvec& GetKetIndex_qq() const;


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
};







struct ThreeBodyChannel
{

   int twoJ;
   int parity;
   int twoTz;

   ModelSpace * modelspace;
   int NumberKets;  // Number of pqr configs that participate in this channel
   std::vector<int> KetList; // eg [2, 4, 7, ...] Used for looping over all the kets in the channel
   std::vector<int> KetMap;  // eg [ -1, -1, 0, -1, 1, -1, -1, 2 ...] Used for asking what is the local index of this ket. -1 means the ket doesn't participate in this channel

   ThreeBodyChannel(){};
   ThreeBodyChannel(int twoj, int p, int twotz, ModelSpace* ms);
   ThreeBodyChannel(int ch, ModelSpace* ms);
   void Initialize();

   int GetNumber3bKets();

   size_t GetLocalIndex( int p, int q, int r, int Jpq );
   bool CheckChannel_ket( Ket3& ket3);

   Ket3& GetKet(size_t iket) ;
   size_t GetNumberKets();



};






class ModelSpace
{
 public:

   static double OCC_CUT;
   static int NOT_AN_ORBIT;
   static Orbit NULL_ORBIT;

   // Data members

   int Emax; // maximum 2*n + l
   int E2max; // cut on 2-particle eneriges
   int E3max; // cut on 3-particle energies
   int Lmax; // maximum L for a single orbit
   int Lmax2; // maximum L included in a 2-body state (not used)?
   int Lmax3; // maximum L included in a 3-body state (not used)?
   int OneBodyJmax; // maximum J for a single orbit
   int TwoBodyJmax; // maximum J for a 2-body state
   int ThreeBodyJmax; // maximum J for a 3-body state (not used)?
   int EmaxUnocc; // Separate emax cut for orbits with l,j,tz that aren't present in the HF reference

   double dE3max; //  cut on three-body configurations which are considered in the IMSRG(3) commutators, taken relative to the fermi energy.
   double occnat3cut; //  cut on three-body configurations which are considered in the IMSRG(3) commutators, taken relative to the fermi energy.
//   double e_fermi;    // The fermi energy, probably in oscillator units
   std::map<int,double> e_fermi;    // The fermi energy, probably in oscillator units. It's different for protons and neutrons, so index by tz

   int norbits;       // number of single-particle orbits (not counting the degenerate m-projections)
   double hbar_omega; // oscillator frequency in MeV
   int target_mass;   // Particle number of the system we're trying to compute
   int target_Z;      // Proton number of the system we're trying to compute
   int Aref;          // Particle number of the normal-ordering reference
   int Zref;          // Proton number of the normal-ordering reference


   std::vector<Orbit> Orbits; // vector of one-body Orbit structs
   std::vector<Ket> Kets;     // vector of two-body Ket structs
   std::vector<Ket3> Kets3;   // vector of three-body Ket3 structs

   std::map<std::array<int,3>,std::set<index_t> > OneBodyChannels;  // map that takes l,j,tz and gives a list of indicies of orbits with those quantum numbers
//   std::map<std::array<int,3>,std::vector<index_t> > OneBodyChannels;

   int nTwoBodyChannels;  // number of two body channels J,parity,Tz
   std::vector<TwoBodyChannel> TwoBodyChannels;   // vector of TwoBodyChannel structs
   std::vector<TwoBodyChannel_CC> TwoBodyChannels_CC; // Pandya transformed channels (the CC stands for cross-coupled)

   int nThreeBodyChannels;  // number of three body channels J,parity,Tz
   std::vector<ThreeBodyChannel> ThreeBodyChannels; // vector of ThreeBodyChannel structs
   std::map<size_t,size_t> Ket3IndexLookup;
   std::map<size_t,size_t> ThreeBodyChannelLookup;

   std::map< std::array<int,3>, std::map< std::array<size_t,2>,std::array<std::vector<size_t>,2> > > PandyaLookup; // I'm not sure how necessary this is...
   std::unordered_map<size_t,index_t> OrbitLookup; // Map to help look up the orbit index based on its quantum numbers


   std::set<index_t> holes;           // in the reference Slater determinant
   std::set<index_t> particles;       // above the reference Slater determinant
   std::set<index_t> core;            // core for decoupling
   std::set<index_t> valence;         // valence space for decoupling
   std::set<index_t> qspace;          // above the valence space for decoupling
   std::set<index_t> proton_orbits;   // all of the orbits with tz<0
   std::set<index_t> neutron_orbits;  // all the orbits with tz>0
   std::set<index_t> all_orbits;      // all of the orbits, for convenient looping

//   std::vector<index_t> holes;           // in the reference Slater determinant
//   std::vector<index_t> particles;       // above the reference Slater determinant
//   std::vector<index_t> core;            // core for decoupling
//   std::vector<index_t> valence;         // valence space for decoupling
//   std::vector<index_t> qspace;          // above the valence space for decoupling
//   std::vector<index_t> proton_orbits;
//   std::vector<index_t> neutron_orbits;
//   std::vector<index_t> all_orbits;

//   std::set<std::array<int,3>> hole_quantum_numbers; // For checking if an orbit could mix with the hole orbits
   std::set<std::array<int,2>> hole_quantum_numbers; // For checking if an orbit could mix with the hole orbits

   std::vector<index_t> KetIndex_pp; 
   std::vector<index_t> KetIndex_ph;
   std::vector<index_t> KetIndex_hh;
   std::vector<index_t> KetIndex_cc; 
   std::vector<index_t> KetIndex_vc;
   std::vector<index_t> KetIndex_qc;
   std::vector<index_t> KetIndex_vv;
   std::vector<index_t> KetIndex_qv;
   std::vector<index_t> KetIndex_qq;
   std::vector<double> Ket_occ_hh;
   std::vector<double> Ket_unocc_hh;
   std::vector<double> Ket_occ_ph;
   std::vector<double> Ket_unocc_ph;

   std::array< std::array< std::unordered_map<index_t,index_t>, 3>,3> MonopoleKets; //List of kets of a given Tz,parity


   std::vector<unsigned int> SortedTwoBodyChannels;
   std::vector<unsigned int> SortedTwoBodyChannels_CC;

   static std::map< std::string, std::vector<std::string> > ValenceSpaces;


   bool sixj_has_been_precalculated;
   bool moshinsky_has_been_precalculated;
   bool scalar_transform_first_pass;
   bool scalar3b_transform_first_pass;
   std::vector<bool> tensor_transform_first_pass;
   bool single_species; // Is there only one kind of fermion?
   IMSRGProfiler profiler;

   static std::unordered_map<uint64_t,double> SixJList;
   static std::unordered_map<uint64_t,double> NineJList;
   static std::unordered_map<uint64_t,double> MoshList;
   static std::unordered_map<uint64_t,double> T3bList; // Tcoefficients for transforming Jacobi 3b to lab-frame 3b




   // Constructors
//   ~ModelSpace();
   ModelSpace();
   ModelSpace(const ModelSpace&); // copy constructor
   ModelSpace( ModelSpace&&); // move constructor
   ModelSpace(int emax, std::vector<std::string> hole_list, std::vector<std::string> valence_list);
   ModelSpace(int emax, std::vector<std::string> hole_list, std::vector<std::string> core_list, std::vector<std::string> valence_list);
   ModelSpace(int emax, std::string reference, std::string valence);
   ModelSpace(int emax, std::string reference);


   // Overloaded operators
   ModelSpace operator=(const ModelSpace&); 
   ModelSpace operator=(ModelSpace&&); 

   // Methods

   void SetUpOrbits( );

   void Init(int emax, std::string reference, std::string valence);
   void Init(int emax, std::map<index_t,double> hole_list, std::string valence);
   void Init(int emax, std::map<index_t,double> hole_list, std::vector<index_t> core_list, std::vector<index_t> valence_list);// keep this for backward-compatibility
   void Init(int emax, std::map<index_t,double> hole_list, std::set<index_t> core_list, std::set<index_t> valence_list);
   void Init(int emax, std::vector<std::string> hole_list, std::vector<std::string> core_list, std::vector<std::string> valence_list);
   void Init_occ_from_file(int emax, std::string valence, std::string occ_file);
   void InitSingleSpecies( int emax, std::string reference, std::string valence); // Work with just one type of fermion
//   void InitSingleSpecies(int emax, std::string reference, std::string valence); // Work with just one type of fermion

//   std::vector<index_t> GetOrbitsAZ(int A, int Z);
   std::map<index_t,double> GetOrbitsAZ(int A, int Z);
   void GetAZfromString(std::string str, int& A, int& Z);
   std::vector<index_t> String2Index( std::vector<std::string> vs );
   std::string Index2String(index_t ind);
   void Get0hwSpace(int Aref, int Zref, std::set<index_t>& core_list, std::set<index_t>& valence_list);
//   void Get0hwSpace(int Aref, int Zref, std::vector<index_t>& core_list, std::vector<index_t>& valence_list);
   void ParseCommaSeparatedValenceSpace(std::string valence, std::set<index_t>& core_list, std::set<index_t>& valence_list);
//   void ParseCommaSeparatedValenceSpace(std::string valence, std::vector<index_t>& core_list, std::vector<index_t>& valence_list);

   void SetupKets();
   void Setup3bKets();
   void AddOrbit(Orbit orb);
   void AddOrbit(int n, int l, int j2, int tz2, double occ, int io);
   void FindEFermi();
   // Setter/Getters
//   Orbit& GetOrbit(int i) {return (Orbit&) Orbits[i];}; 
   Orbit& GetOrbit(int i); 
   Ket& GetKet(int i) const {return (Ket&) Kets[i];};
   Ket& GetKet(int p, int q) const {return (Ket&) Kets[Index2(p,q)];};
   Ket3& GetKet3(int i) const {return (Ket3&) Kets3[i];};
   size_t GetOrbitIndex(int n, int l, int j2, int tz2) const {return Index1(n,l,j2,tz2);};
   size_t GetKetIndex(int p, int q) const {return Index2(p,q);}; // convention is p<=q
   size_t GetKetIndex(Ket * ket) const {return Index2(ket->p,ket->q);}; // convention is p<=q
   size_t GetKet3Index(int p, int q, int r, int Jpq){return Ket3IndexLookup.at(Ket3IndexHash(p,q,r,Jpq));};
   size_t Ket3IndexHash(size_t p, size_t q, size_t r, size_t Jpq);
   size_t ThreeBodyChannelHash( int twoJ, int parity, int twoTz);

   void SetEmaxUnocc(int e);

   size_t GetNumberOrbits() const {return norbits;};
   size_t GetNumberKets() const {return Kets.size();};

   void SetHbarOmega(double hw) {hbar_omega = hw;};
   void SetTargetMass(int A) {target_mass = A;};
   void SetTargetZ(int Z) {target_Z = Z;};
   double GetHbarOmega() const {return hbar_omega;};
   int GetTargetMass() const {return target_mass;};
   int GetTargetZ() const {return target_Z;};
   int GetAref() const {return Aref;};
   int GetZref() const {return Zref;};
   size_t GetNumberTwoBodyChannels() const {return TwoBodyChannels.size();};
   size_t GetNumberThreeBodyChannels() const {return ThreeBodyChannels.size();};

   TwoBodyChannel& GetTwoBodyChannel(int ch) const {return (TwoBodyChannel&) TwoBodyChannels[ch];};
   TwoBodyChannel_CC& GetTwoBodyChannel_CC(int ch) const {return (TwoBodyChannel_CC&) TwoBodyChannels_CC[ch];};

   ThreeBodyChannel& GetThreeBodyChannel(int ch) const {return (ThreeBodyChannel&) ThreeBodyChannels[ch];};

   int GetTwoBodyJmax() const {return TwoBodyJmax;};
   int GetThreeBodyJmax() const {return ThreeBodyJmax;};
   void SetReference(std::vector<index_t>); // For backwards compatibility
   void SetReference(std::set<index_t>);
   void SetReference(std::map<index_t,double>);
   void SetReference(std::string);

   int GetEmax(){return Emax;};
   int GetE2max(){return E2max;};
   int GetE3max(){return E3max;};
   int GetLmax2(){return Lmax2;};
   int GetLmax(){return Lmax;};
   int GetLmax3(){return Lmax3;};
   void SetEmax(int e){Emax=e;};
   void SetE2max(int e){E2max=e;};
//   void SetE3max(int e){E3max=e;};
   void SetE3max(int e);
//   void SetLmax(int l){Lmax=l;};
   void SetLmax(int l);
   void SetLmax2(int l){Lmax2=l;};
   void SetLmax3(int l){Lmax3=l;};
   void SetdE3max(double e){dE3max = e;};
   void SetOccNat3Cut(double o){occnat3cut = o;};
   double GetdE3max(){return dE3max;};
   double GetOccNat3Cut(){return occnat3cut;}; // setting this to zero or less makes no cut, setting to 0.25 cuts out everything
   void SetEFermi(double ef){e_fermi[-1]=ef; e_fermi[+1]=ef;};
   std::map<int,double> GetEFermi(){ return e_fermi ;};
   void SetEFermi(double ef_proton, double ef_neutron){e_fermi[-1] = ef_proton; e_fermi[1]=ef_neutron;};

   double GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3);
   double GetNineJ(double j1, double j2, double j3, double j4, double j5, double j6, double j7, double j8, double j9);
   double GetMoshinsky( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L); // Inconsistent notation. Not ideal.
   bool SixJ_is_empty(){ return SixJList.empty(); };
   double GetT3b(index_t a,index_t b, index_t c, int Jab, int J, int N1, int L1, int S1, int J1, int N2, int L2,int J2,int J12,int Ncm,int Lcm); // 15 index object. yikes.

   size_t GetOrbitIndex(std::string);
   size_t GetTwoBodyChannelIndex(int j, int p, int t);
   void UnpackTwoBodyChannelIndex( size_t ch, int& j, int& p, int& tz);
   int phase(int x) {return (x%2)==0 ? 1 : -1;};
//   int phase(double x) {return phase(int(x));};

   
   size_t GetThreeBodyChannelIndex( int twoJ, int parity, int twoTz );
   size_t CountThreeBodyStatesInsideCut();


//   size_t Index1(int n, int l, int j2, int tz2) const {return(2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2 ;};
   size_t Index1(int n, int l, int j2, int tz2) const ;
   size_t Index1_hash(int n, int l, int j2, int tz2) const ;
//   size_t Index1(int n, int l, int j2, int tz2) const { size_t indx = (2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2; return (single_species ? indx/2 : indx) ;};
//   inline int Index2(int p, int q) const {return q*(q+1)/2 + p;};
//   size_t Index2(size_t p, size_t q) const {return p*(2*norbits-1-p)/2 + q;};
   size_t Index2(size_t p, size_t q) const ;
//   size_t Index2(size_t p, size_t q) const {return p*(2*all_orbits.size()-1-p)/2 + q;};

   void PreCalculateMoshinsky();
   void PreCalculateSixJ();
   void ClearVectors();
   void ResetFirstPass();
   void CalculatePandyaLookup(int rank_J, int rank_T, int parity); // construct a lookup table for more efficient pandya transformation
//   map<array<int,2>,vector<array<int,2>>>& GetPandyaLookup(int rank_J, int rank_T, int parity);
//   std::map<std::array<int,2>,std::array<std::vector<int>,2>>& GetPandyaLookup(int rank_J, int rank_T, int parity);
   std::map<std::array<size_t,2>,std::array<std::vector<size_t>,2>>& GetPandyaLookup(int rank_J, int rank_T, int parity);
   uint64_t SixJHash(double j1, double j2, double j3, double J1, double J2, double J3);
   void SixJUnHash(uint64_t key, uint64_t& j1, uint64_t& j2, uint64_t& j3, uint64_t& J1, uint64_t& J2, uint64_t& J3);
   uint64_t MoshinskyHash(uint64_t N,uint64_t Lam,uint64_t n,uint64_t lam,uint64_t n1,uint64_t l1,uint64_t n2,uint64_t l2,uint64_t L);
   void MoshinskyUnHash(uint64_t key,uint64_t& N,uint64_t& Lam,uint64_t& n,uint64_t& lam,uint64_t& n1,uint64_t& l1,uint64_t& n2,uint64_t& l2,uint64_t& L);




};




#endif

