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
#include <armadillo>
#include "IMSRGProfiler.hh"
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif
#define OCC_CUT 1e-6


//using namespace std;

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
   int delta_pq() const {return dpq;};

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
   virtual ~TwoBodyChannel();
   TwoBodyChannel();
   TwoBodyChannel(int j, int p, int t, ModelSpace* ms);
   TwoBodyChannel(int N, ModelSpace* ms);
   void Initialize(int N, ModelSpace* ms);


   //Methods
   int GetNumberKets() const {return NumberKets;};
   int GetLocalIndex(int ketindex) const { return KetMap[ketindex];}; // modelspace ket index => local ket index
   int GetLocalIndex(int p, int q) const ;
   int GetKetIndex(int i) const { return KetList[i];}; // local ket index => modelspace ket index
   const Ket& GetKet(int i) const  ; // get pointer to ket using local index
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


// private:
   //Fields
   ModelSpace * modelspace;
   int NumberKets;  // Number of pq configs that participate in this channel
   std::vector<int> KetList; // eg [2, 4, 7, ...] Used for looping over all the kets in the channel
   std::vector<int> KetMap;  // eg [ -1, -1, 0, -1, 1, -1, -1, 2 ...] Used for asking what is the local index of this ket. -1 means the ket doesn't participate in this channel
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
   ModelSpace(int emax, std::vector<std::string> hole_list, std::vector<std::string> valence_list);
   ModelSpace(int emax, std::vector<std::string> hole_list, std::vector<std::string> core_list, std::vector<std::string> valence_list);
   ModelSpace(int emax, std::string reference, std::string valence);
   ModelSpace(int emax, std::string reference);

   // Overloaded operators
   ModelSpace operator=(const ModelSpace&); 
   ModelSpace operator=(ModelSpace&&); 

   // Methods

   void Init(int emax, std::string reference, std::string valence);
   void Init(int emax, std::map<index_t,double> hole_list, std::string valence);
   void Init(int emax, std::map<index_t,double> hole_list, std::vector<index_t> core_list, std::vector<index_t> valence_list);
   void Init(int emax, std::vector<std::string> hole_list, std::vector<std::string> core_list, std::vector<std::string> valence_list);
   void Init_occ_from_file(int emax, std::string valence, std::string occ_file);

//   std::vector<index_t> GetOrbitsAZ(int A, int Z);
   std::map<index_t,double> GetOrbitsAZ(int A, int Z);
   void GetAZfromString(std::string str, int& A, int& Z);
   std::vector<index_t> String2Index( std::vector<std::string> vs );
   std::string Index2String(index_t ind);
   void Get0hwSpace(int Aref, int Zref, std::vector<index_t>& core_list, std::vector<index_t>& valence_list);
   void ParseCommaSeparatedValenceSpace(std::string valence, std::vector<index_t>& core_list, std::vector<index_t>& valence_list);

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
   void SetReference(std::vector<index_t>);
   void SetReference(std::map<index_t,double>);
   void SetReference(std::string);

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

   int GetOrbitIndex(std::string);
   int GetTwoBodyChannelIndex(int j, int p, int t);
   inline int phase(int x) {return (x%2)==0 ? 1 : -1;};
   inline int phase(double x) {return phase(int(x));};

   inline int Index1(int n, int l, int j2, int tz2) const {return(2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2 ;};
//   inline int Index2(int p, int q) const {return q*(q+1)/2 + p;};
   inline int Index2(int p, int q) const {return p*(2*norbits-1-p)/2 + q;};

   void PreCalculateMoshinsky();
   void PreCalculateSixJ();
   void ClearVectors();
   void ResetFirstPass();
   void CalculatePandyaLookup(int rank_J, int rank_T, int parity); // construct a lookup table for more efficient pandya transformation
//   map<array<int,2>,vector<array<int,2>>>& GetPandyaLookup(int rank_J, int rank_T, int parity);
   std::map<std::array<int,2>,std::array<std::vector<int>,2>>& GetPandyaLookup(int rank_J, int rank_T, int parity);
   uint64_t SixJHash(double j1, double j2, double j3, double J1, double J2, double J3);
   void SixJUnHash(uint64_t key, uint64_t& j1, uint64_t& j2, uint64_t& j3, uint64_t& J1, uint64_t& J2, uint64_t& J3);
   uint64_t MoshinskyHash(uint64_t N,uint64_t Lam,uint64_t n,uint64_t lam,uint64_t n1,uint64_t l1,uint64_t n2,uint64_t l2,uint64_t L);
   void MoshinskyUnHash(uint64_t key,uint64_t& N,uint64_t& Lam,uint64_t& n,uint64_t& lam,uint64_t& n1,uint64_t& l1,uint64_t& n2,uint64_t& l2,uint64_t& L);


   // Data members
   std::vector<index_t> holes;           // in the reference Slater determinant
//   map<index_t,double> holes;           // in the reference Slater determinant
   std::vector<index_t> particles;       // above the reference Slater determinant
   std::vector<index_t> core;            // core for decoupling
   std::vector<index_t> valence;         // valence space for decoupling
   std::vector<index_t> qspace;          // above the valence space for decoupling
   std::vector<index_t> proton_orbits;
   std::vector<index_t> neutron_orbits;

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

//   array< array< std::vector<index_t>, 2>,3> MonopoleKets; //List of kets of a given Tz,parity
   std::array< std::array< std::unordered_map<index_t,index_t>, 3>,3> MonopoleKets; //List of kets of a given Tz,parity

   int Emax;
   int E2max;
   int E3max;
   int Lmax2;
   int Lmax3;
   int OneBodyJmax;
   int TwoBodyJmax;
   int ThreeBodyJmax;
   std::map<std::array<int,3>,std::vector<index_t> > OneBodyChannels;

   std::vector<unsigned int> SortedTwoBodyChannels;
   std::vector<unsigned int> SortedTwoBodyChannels_CC;

   static std::map< std::string, std::vector<std::string> > ValenceSpaces;
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
   std::vector<Orbit> Orbits;
   std::vector<Ket> Kets;
   std::vector<TwoBodyChannel> TwoBodyChannels;
   std::vector<TwoBodyChannel_CC> TwoBodyChannels_CC;
   std::map< std::array<int,3>, std::map< std::array<int,2>,std::array<std::vector<int>,2> > > PandyaLookup;
   bool sixj_has_been_precalculated;
   bool moshinsky_has_been_precalculated;
   bool scalar_transform_first_pass;
   std::vector<bool> tensor_transform_first_pass;
   IMSRGProfiler profiler;
//   map<long int,double> SixJList;

   static std::unordered_map<uint64_t,double> SixJList;
   static std::unordered_map<uint64_t,double> NineJList;
   static std::unordered_map<uint64_t,double> MoshList;

};




#endif

