#include "ModelSpace.hh"
#include "AngMom.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

Orbit::~Orbit()
{
//  cout << "In Orbit destructor" << endl;
}

Orbit::Orbit()
: n(-1), l(-1), j2(-1), tz2(-1),ph(-1),cvq(-1),index(-1)
//: n(-1), l(-1), j2(-1), tz2(-1),ph(-1),io(-1),index(-1)
{}

Orbit::Orbit(int n, int l, int j2, int tz2, int ph, int cvq, int index)
: n(n), l(l), j2(j2), tz2(tz2),ph(ph),cvq(cvq),index(index)
//: n(n), l(l), j2(j2), tz2(tz2),ph(ph),io(io),index(index)
{}

Orbit::Orbit(const Orbit& orb)
: n(orb.n), l(orb.l), j2(orb.j2), tz2(orb.tz2),ph(orb.ph),cvq(orb.cvq),index(orb.index)
//: n(orb.n), l(orb.l), j2(orb.j2), tz2(orb.tz2),ph(orb.ph),io(orb.io),index(orb.index)
{}


//************************************************************************
//************************************************************************
//************************************************************************
Ket::~Ket()
{
//  cout << "In Ket destructor" << endl;
}

Ket::Ket()
{}

Ket::Ket(Orbit& op_in, Orbit& oq_in)
: op(&op_in), oq(&oq_in), p(op_in.index), q(oq_in.index)
{
   phase_prefactor = ((op->j2+oq->j2)/2 + 1) % 2==0 ? 1 : -1;
   dpq = p==q ? 1 : 0;
}

int Ket::Phase(int J)
{
   return phase_prefactor * (J%2==0 ? 1 : -1);
}

//************************************************************************
//************************************************************************
//************************************************************************

TwoBodyChannel::~TwoBodyChannel()
{
//  cout << "In TwoBodyChannel destructor" << endl;
}

TwoBodyChannel::TwoBodyChannel()
{}

TwoBodyChannel::TwoBodyChannel(int j, int p, int t, ModelSpace *ms)
{
  Initialize(ms->GetTwoBodyChannelIndex(j,p,t), ms);
}

TwoBodyChannel::TwoBodyChannel(int N, ModelSpace *ms)
{
   Initialize(N,ms);
}

void TwoBodyChannel::Initialize(int N, ModelSpace *ms)
{
   int tbjmax = ms->TwoBodyJmax;
//   J = N%tbjmax;
//   parity = (N/tbjmax)%2;
//   Tz = (N/(2*tbjmax)-1);
   J = N%(tbjmax+1);
   parity = (N/(tbjmax+1))%2;
   Tz = (N/(2*(tbjmax+1))-1);
   modelspace = ms;
   NumberKets = 0;
   int nk = modelspace->GetNumberKets();
   KetMap.resize(nk,-1); // set all values to -1
   for (int i=0;i<nk;i++)
   {
      Ket &ket = modelspace->GetKet(i);
      if ( CheckChannel_ket(ket) )
      {
         KetMap[i] = NumberKets;
         KetList.push_back(i);
         NumberKets++;
      }
   }
   KetIndex_pp = GetKetIndexFromList(modelspace->KetIndex_pp);
   KetIndex_hh = GetKetIndexFromList(modelspace->KetIndex_hh);
   KetIndex_ph = GetKetIndexFromList(modelspace->KetIndex_ph);
   KetIndex_cc = GetKetIndexFromList(modelspace->KetIndex_cc);
   KetIndex_vc = GetKetIndexFromList(modelspace->KetIndex_vc);
   KetIndex_qc = GetKetIndexFromList(modelspace->KetIndex_qc);
   KetIndex_vv = GetKetIndexFromList(modelspace->KetIndex_vv);
   KetIndex_qv = GetKetIndexFromList(modelspace->KetIndex_qv);
   KetIndex_qq = GetKetIndexFromList(modelspace->KetIndex_qq);
//   KetIndex_c_c = GetKetIndexFromList(modelspace->KetIndex_c_c);
//   KetIndex_q_q = GetKetIndexFromList(modelspace->KetIndex_q_q);
//   KetIndex_q_c = GetKetIndexFromList(modelspace->KetIndex_q_c);
//   KetIndex_v_c = GetKetIndexFromList(modelspace->KetIndex_v_c);
//   KetIndex_v_q= GetKetIndexFromList(modelspace->KetIndex_v_q);
}


int TwoBodyChannel::GetLocalIndex(int p, int q) const { return KetMap[modelspace->GetKetIndex(p,q)];}; 

// get pointer to ket using local index
Ket & TwoBodyChannel::GetKet(int i) const { return modelspace->GetKet(KetList[i]);}; 


//bool TwoBodyChannel::CheckChannel_ket(int p, int q) const
bool TwoBodyChannel::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   if ((op->index==oq->index) and (J%2 != 0)) return false; // Pauli principle
   if ((op->l + oq->l)%2 != parity) return false;
   if ((op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}

arma::uvec& TwoBodyChannel::GetKetIndex_pp() { return KetIndex_pp;};
arma::uvec& TwoBodyChannel::GetKetIndex_hh() { return KetIndex_hh;};
arma::uvec& TwoBodyChannel::GetKetIndex_ph() { return KetIndex_ph;};
arma::uvec& TwoBodyChannel::GetKetIndex_cc() { return KetIndex_cc;};
arma::uvec& TwoBodyChannel::GetKetIndex_vc() { return KetIndex_vc;};
arma::uvec& TwoBodyChannel::GetKetIndex_qc() { return KetIndex_qc;};
arma::uvec& TwoBodyChannel::GetKetIndex_vv() { return KetIndex_vv;};
arma::uvec& TwoBodyChannel::GetKetIndex_qv() { return KetIndex_qv;};
arma::uvec& TwoBodyChannel::GetKetIndex_qq() { return KetIndex_qq;};
//arma::uvec& TwoBodyChannel::GetKetIndex_c_c() { return KetIndex_c_c;};
//arma::uvec& TwoBodyChannel::GetKetIndex_q_q() { return KetIndex_q_q;};
//arma::uvec& TwoBodyChannel::GetKetIndex_q_c() { return KetIndex_q_c;};
//arma::uvec& TwoBodyChannel::GetKetIndex_v_c() { return KetIndex_v_c;};
//arma::uvec& TwoBodyChannel::GetKetIndex_v_q(){ return KetIndex_v_q;};


arma::uvec TwoBodyChannel::GetKetIndexFromList(vector<index_t>& vec_in)
{
   vector<index_t> index_list (min(vec_in.size(),KetList.size()));
   auto it = set_intersection(KetList.begin(),KetList.end(),vec_in.begin(),vec_in.end(),index_list.begin());
   index_list.resize(it-index_list.begin());
   for (auto& x : index_list)
   {
     x = KetMap[x];
   }
   return arma::uvec(index_list);
}

//************************************************************************
//************************************************************************
//************************************************************************

TwoBodyChannel_CC::~TwoBodyChannel_CC()
{
//   cout << "In TwoBodyChannel_CC destructor" << endl;
}

TwoBodyChannel_CC::TwoBodyChannel_CC()
{}

TwoBodyChannel_CC::TwoBodyChannel_CC(int j, int p, int t, ModelSpace *ms)
{
  Initialize(ms->GetTwoBodyChannelIndex(j,p,t), ms);
}

TwoBodyChannel_CC::TwoBodyChannel_CC(int N, ModelSpace *ms)
{
   Initialize(N,ms);
}

// Check if orbits pq participate in this cross-coupled two-body channel
// Difference from regular channels:
// no Pauli rule, <pp||nn> is allowed.
bool TwoBodyChannel_CC::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   if ((op->l + oq->l)%2 != parity)    return false;
   if (abs(op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)          return false;
   if (abs(op->j2 - oq->j2) > 2*J)     return false;

   return true;
}


//************************************************************************
//************************************************************************

// Static members

unordered_map<unsigned long int,double> ModelSpace::SixJList;
unordered_map<unsigned long long int,double> ModelSpace::NineJList;
unordered_map<unsigned long long int,double> ModelSpace::MoshList;
map<string,vector<string>> ModelSpace::ValenceSpaces  {
{ "p-shell"  ,         {"p0p3","n0p3","p0p1","n0p1"}},
{ "sd-shell"  ,        {"p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"}},
{ "psd-shell"  ,       {"p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"}},
{ "fp-shell"  ,        {"p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"}},
{ "sdfp-shell"  ,      {"p0d5","n0d5","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"}},
{ "fpg9-shell"  ,      {"p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1","p0g9","n0g9"}},
{ "sd3f7p3-shell"  ,   {"p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p1p3","n1p3"}},
};



ModelSpace::~ModelSpace()
{
//  cout << "In ModelSpace destructor. emax = " << Nmax << endl;
}

ModelSpace::ModelSpace()
{
  cout << "In default constructor" << endl;
   norbits = 0;
   OneBodyJmax = 0;
   TwoBodyJmax = 0;
   ThreeBodyJmax = 0;
   hbar_omega=20;
   target_mass = 16;
   Nmax = 0;
   N2max = 0;
   N3max = 0;
}


ModelSpace::ModelSpace(const ModelSpace& ms)
 :
   holes( ms.holes), particles( ms.particles),
   core(ms.core), valence(ms.valence), qspace( ms.qspace), 
   proton_orbits( ms.proton_orbits),neutron_orbits( ms.neutron_orbits),
   KetIndex_pp( ms.KetIndex_pp), KetIndex_ph( ms.KetIndex_ph), KetIndex_hh( ms.KetIndex_hh),
   KetIndex_cc( ms.KetIndex_cc),
   KetIndex_vc( ms.KetIndex_vc),
   KetIndex_qc( ms.KetIndex_qc),
   KetIndex_vv( ms.KetIndex_vv),
   KetIndex_qv( ms.KetIndex_qv),
   KetIndex_qq( ms.KetIndex_qq),
   Nmax(ms.Nmax), N2max(ms.N2max), N3max(ms.N3max),
   OneBodyJmax(ms.OneBodyJmax), TwoBodyJmax(ms.TwoBodyJmax), ThreeBodyJmax(ms.ThreeBodyJmax),
   OneBodyChannels(ms.OneBodyChannels),
   SortedTwoBodyChannels(ms.SortedTwoBodyChannels),
   SortedTwoBodyChannels_CC(ms.SortedTwoBodyChannels_CC),
   norbits(ms.norbits), hbar_omega(ms.hbar_omega),
   target_mass(ms.target_mass), target_Z(ms.target_Z), Aref(ms.Aref), Zref(ms.Zref),
   Orbits(ms.Orbits), Kets(ms.Kets),
   TwoBodyChannels(ms.TwoBodyChannels), TwoBodyChannels_CC(ms.TwoBodyChannels_CC)
{
//   cout << "In ModelSpace copy constructor" << endl;
//   cout << "Orbits size:  " << ms.Orbits.size() << " -> " << Orbits.size() << endl;
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
}

ModelSpace::ModelSpace(ModelSpace&& ms)
 :
   holes( move(ms.holes)), particles( move(ms.particles)),
   core(move(ms.core)), valence(move(ms.valence)),  qspace( move(ms.qspace)),  
   proton_orbits( move(ms.proton_orbits)),
   neutron_orbits( move(ms.neutron_orbits)),
   KetIndex_pp( move(ms.KetIndex_pp)), KetIndex_ph( move(ms.KetIndex_ph)), KetIndex_hh( move(ms.KetIndex_hh)),
   KetIndex_cc( ms.KetIndex_cc),
   KetIndex_vc( ms.KetIndex_vc),
   KetIndex_qc( ms.KetIndex_qc),
   KetIndex_vv( ms.KetIndex_vv),
   KetIndex_qv( ms.KetIndex_qv),
   KetIndex_qq( ms.KetIndex_qq),
//   KetIndex_vv( move(ms.KetIndex_vv)), KetIndex_c_c( move(ms.KetIndex_c_c)),
//   KetIndex_q_q( move(ms.KetIndex_q_q)),
//   KetIndex_q_c( move(ms.KetIndex_q_c)),
//   KetIndex_v_c( move(ms.KetIndex_v_c)), KetIndex_v_q( move(ms. KetIndex_v_q)),
   Nmax(ms.Nmax), N2max(ms.N2max), N3max(ms.N3max),
   OneBodyJmax(ms.OneBodyJmax), TwoBodyJmax(ms.TwoBodyJmax), ThreeBodyJmax(ms.ThreeBodyJmax),
   OneBodyChannels(move(ms.OneBodyChannels)),
   SortedTwoBodyChannels(move(ms.SortedTwoBodyChannels)),
   SortedTwoBodyChannels_CC(move(ms.SortedTwoBodyChannels_CC)),
   norbits(ms.norbits), hbar_omega(ms.hbar_omega),
   target_mass(ms.target_mass), target_Z(ms.target_Z), Aref(ms.Aref), Zref(ms.Zref),
   Orbits(move(ms.Orbits)), Kets(move(ms.Kets)),
   TwoBodyChannels(move(ms.TwoBodyChannels)), TwoBodyChannels_CC(move(ms.TwoBodyChannels_CC))
{
//   cout << "In ModelSpace move constructor" << endl;
//   cout << "Orbits size:  " << ms.Orbits.size() << " -> " << Orbits.size() << endl;
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
   for (TwoBodyChannel& tbc : ms.TwoBodyChannels)   tbc.modelspace = NULL;
   for (TwoBodyChannel_CC& tbc_cc : ms.TwoBodyChannels_CC)   tbc_cc.modelspace = NULL;
}


// orbit string representation is e.g. p0f7
// Assumes that the core is hole states that aren't in the valence space.
ModelSpace::ModelSpace(int nmax, vector<string> hole_list, vector<string> valence_list)
: Nmax(nmax), N2max(2*nmax), N3max(3*nmax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0),hbar_omega(20),target_mass(16)
{
   Init(nmax, hole_list, valence_list); 
}

// If we don't want the reference to be the core
ModelSpace::ModelSpace(int nmax, vector<string> hole_list, vector<string> core_list, vector<string> valence_list)
: Nmax(nmax), N2max(2*nmax), N3max(3*nmax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0),hbar_omega(20),target_mass(16)
{
   Init(nmax, hole_list, core_list, valence_list); 
}

// DO I NEED THIS? I guess it doesn't hurt too much...
ModelSpace::ModelSpace(int nmax, int A, int Z)
: Nmax(nmax), N2max(2*nmax), N3max(3*nmax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0),hbar_omega(20),target_mass(A),target_Z(Z)
{
   Init(nmax,GetOrbitsAZ(A,Z),GetOrbitsAZ(A,Z),{});
}

// Shortcuts for common modelspaces
ModelSpace::ModelSpace(int nmax, string reference, string valence)
: Nmax(nmax), N2max(2*nmax), N3max(3*nmax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0),hbar_omega(20)
{
  Init(nmax,reference,valence);
}

ModelSpace::ModelSpace(int nmax, string valence)
: Nmax(nmax), N2max(2*nmax), N3max(3*nmax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0),hbar_omega(20)
{
  Init(nmax,valence);
}


// Simple initialization which makes some assumptions
void ModelSpace::Init(int nmax, string valence)
{
  int Ac,Zc=-1;
  vector<index_t> hole_list, valence_list, core_list;

  auto itval = ValenceSpaces.find( valence );

  // if we have a valence space, like sd-shell, then we assume
  // that holes=core and the core consists of all orbits below
  // the lowest valence orbit
  if ( itval != ValenceSpaces.end() )
  {
     valence_list = String2Index(itval->second);
     int z=0;
     while ( find(valence_list.begin(),valence_list.end(),z) == valence_list.end() )
     {
       hole_list.push_back(z);
       z+=2;
     }
     target_Z = hole_list.size();
     int n=1;
     while ( find(valence_list.begin(),valence_list.end(),n) == valence_list.end() )
     {
       hole_list.push_back(n);
       n+=2;
     }
     target_mass = hole_list.size();
     core_list = hole_list;
     for (auto& v : valence_list) core_list.resize( remove(core_list.begin(),core_list.end(),v) - core_list.begin() );
//     cout << "Valence space: " << valence << endl;
  }
  // otherwise, there's no valence space and we have a single reference
  // with holes = core.
  else
  {
     GetAZfromString(valence,Ac,Zc);
//     cout << "Core: " << valence << "  => A = " << Ac << " Z = " << Zc << endl;
     core_list = GetOrbitsAZ(Ac,Zc);
     hole_list = core_list;
     target_mass = Ac;
     target_Z = Zc; 
  }

//  cout << "holes: ";
//  for ( auto &h : hole_list ) cout << h << " ";
//  cout << endl;
//  cout << "core: ";
//  for ( auto &h : core_list ) cout << h << " ";
//  cout << endl;
//  cout << "valence: ";
//  for ( auto &h : valence_list ) cout << h << " ";
//  cout << endl;
  Init(nmax,hole_list,core_list,valence_list);
}

// specify the reference and either the core or valence
void ModelSpace::Init(int nmax, string reference, string valence)
{
  int Aref,Zref;
  int Ac,Zc=-1;
  vector<index_t> hole_list, valence_list, core_list;
  GetAZfromString(reference,Aref,Zref);
//  cout << "Reference: " << reference << "  => A = " << Aref << " Z = " << Zref << endl;
  hole_list = GetOrbitsAZ(Aref,Zref);

  auto itval = ValenceSpaces.find(valence);

  if ( itval != ValenceSpaces.end() ) // we've got a valence space
  {
     valence_list = String2Index(itval->second);
     core_list = hole_list;
     for (auto& v : valence_list) core_list.resize( remove(core_list.begin(),core_list.end(),v) - core_list.begin() );
//     cout << "Valence space: " << valence << endl;
     target_mass = Aref;
     target_Z = Zref;
  }
  else  // no valence space. we've got a single-reference.
  {
     GetAZfromString(valence,Ac,Zc);
//     cout << "Core: " << valence << "  => A = " << Ac << " Z = " << Zc << endl;
     core_list = GetOrbitsAZ(Ac,Zc);
     target_mass = Ac;
     target_Z = Zc; 
  }

  Init(nmax,hole_list,core_list,valence_list);
  
}


void ModelSpace::Init(int nmax, vector<string> hole_list, vector<string> valence_list)
{
   // Assume core is hole states that aren't in the valence space
   vector<string> core_list( hole_list );
   for (auto& v : valence_list) core_list.resize( remove(core_list.begin(),core_list.end(),v) - core_list.begin() );
   Init(nmax,hole_list,core_list,valence_list);
}

void ModelSpace::Init(int nmax, vector<string> hole_list, vector<string> core_list, vector<string> valence_list)
{
   cout << "Creating a model space with Nmax = " << Nmax << "  and hole orbits [";
   for (auto& h : hole_list)  cout << h << " ";
   cout << "]   and core orbits [";
   for (auto& c : core_list)    cout << c << " ";
   cout << "]   and valence orbits [";
   for (auto& v : valence_list)   cout << v << " ";
   cout << "]" << endl;
  Init(nmax, String2Index(hole_list), String2Index(core_list), String2Index(valence_list) );
}

// This is the Init which should inevitably be called
void ModelSpace::Init(int nmax, vector<index_t> hole_list, vector<index_t> core_list, vector<index_t> valence_list)
{
   Orbits.clear();
   particles.clear();
   holes.clear();
   core.clear();
   valence.clear();
//   holes = hole_list;
//   core = core_list;
//   valence = valence_list;
   qspace.clear();
   proton_orbits.clear();
   neutron_orbits.clear();
   OneBodyChannels.clear();
//   cout << "Init" << endl;



   // Make sure no orbits are both core and valence
   for (auto& c : core_list)
   {
     if ( find(valence_list.begin(), valence_list.end(), c) != valence_list.end() )
       cout << "!!!!!!!!!!!!! ModelSpace::Init : Conflicting definition. Orbit " << c << " is in core and valence spaces." << endl;
   }
   

   norbits = (Nmax+1)*(Nmax+2);
   Orbits.resize(norbits);
   for (int N=0; N<=Nmax; ++N)
   {
     for (int l=N; l>=0; l-=2)
     {
       int n = (N-l)/2;
       for (int j2=2*l+1; j2>=2*l-1 and j2>0; j2-=2)
       {
         for (int tz : {-1, 1} )
         {
            int ph = 0;
            int cvq = 2;
            int indx = Index1(n,l,j2,tz);
            if ( find(hole_list.begin(), hole_list.end(), indx) != hole_list.end() ) ph=1; // hole orbit
            if ( find(core_list.begin(), core_list.end(), indx) != core_list.end() ) cvq=0; // core orbit
            if ( find(valence_list.begin(), valence_list.end(), indx) != valence_list.end() ) cvq=1; // valence orbit
            AddOrbit(n,l,j2,tz,ph,cvq);


         }
       }
     }
   }
   SetupKets();
   target_mass = valence_list.size() > 0 ? 2 : 0;
   target_Z = 0;
   for (auto& c : core_list )
   {
     target_mass += (GetOrbit(c).j2+1);
     if (c%2==0) target_Z += (GetOrbit(c).j2+1);
   }
   Aref = 0;
   Zref = 0;
   for (auto& c : hole_list )
   {
     Aref += (GetOrbit(c).j2+1);
     if (c%2==0) Zref += (GetOrbit(c).j2+1);
   }


//   cout << "ModelSpace main Init:  target_mass = " << target_mass << "  target_Z = " << target_Z << endl;
//   cout << "ModelSpace main Init:  Aref = " << Aref << "  Zref = " << Zref << endl;
}


vector<index_t> ModelSpace::String2Index( vector<string> vs )
{
  vector<index_t> vi;
  vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};

  for ( auto& s : vs )
  {
    int tz2 = s[0]=='p' ? -1 : 1;
    int n,j2;
    istringstream( s.substr(1,2) ) >> n;
    int l = find(l_list.begin(),l_list.end(), s[2]) - l_list.begin();
    istringstream( s.substr(3,s.size()) ) >> j2;
    vi.push_back( Index1(n,l,j2,tz2) );
//    cout << s << " -> " << n << " " << l << " " << j2 << " " << tz2 << endl;
  }
  return vi;
}



void ModelSpace::GetAZfromString(string str,int& A, int& Z) // TODO: accept different formats, e.g. 22Na vs Na22
{
  vector<string> periodic_table = {"n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
                        "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
                        "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
                        "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
                        "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb"};
  int i=0;
  while (! isdigit(str[i])) i++;
  string elem = str.substr(0,i);
  stringstream( str.substr(i,str.size()-i)) >> A;
  auto it_elem = find(periodic_table.begin(),periodic_table.end(),elem);
  if (it_elem != periodic_table.end())
  {
    Z = it_elem - periodic_table.begin();
  }
  else
  {
    Z =-1;
   cout << "ModelSpace::GetAZfromString :  Trouble parsing " << str << endl;
  }
//  cout << "GetAZfromString:  " << str << "  -> " << A << " " << Z << endl;
}

// Fill A orbits with Z protons and A-Z neutrons
// assuming a standard shell-model level ordering
vector<index_t> ModelSpace::GetOrbitsAZ(int A, int Z)
{
  int zz = 0;
  int nn = 0; // unfortunate there are so many n's here...
  vector<index_t> orbitsAZ;
  for (int N=0; N<=Nmax; ++N)
  {
    for (int l=N; l>=0; l-=2)
    {
      int n = (N-l)/2;
      int j2 = 2*l+1;
      if (zz < Z)
      {
        orbitsAZ.push_back(Index1(n,l,j2,-1));
        zz += j2+1;
      }
      if (nn < A-Z)
      {
        orbitsAZ.push_back(Index1(n,l,j2,1));
        nn += j2+1;
      }
      if (zz==Z and nn==A-Z) return orbitsAZ; // We're all done here.
      if (zz>Z or nn>A-Z) // Oops. We partially filled a shell.
      {
        cout << "Trouble! No support yet for partially-filled shells!! (A = " <<A << ", Z = " << Z << ") zz = " << zz << "  nn = " << nn <<  endl;
        return orbitsAZ;
      }
    }
    for (int l=N%2; l<=N; l+=2)
    {
      if (l==0) continue;
      int n = (N-l)/2;
      int j2 = 2*l-1;
      if (zz < Z)
      {
        orbitsAZ.push_back(Index1(n,l,j2,-1));
        zz += j2+1;
      }
      if (nn < A-Z)
      {
        orbitsAZ.push_back(Index1(n,l,j2,1));
        nn += j2+1;
      }
//      cout << "A,Z = " << A << "," << Z << "  nn,zz = " << nn << "," << zz << endl;
//      for ( auto& o : orbitsAZ ) cout << o << " ";
//      cout << endl;
      if (zz==Z and nn==A-Z) return orbitsAZ; // We're all done here.
      if (zz>Z or nn>A-Z) // Oops. We partially filled a shell.
      {
        cout << "Trouble! No support yet for partially-filled shells!! (A = " <<A << ", Z = " << Z << ") zz = " << zz << "  nn = " << nn <<  endl;
        return orbitsAZ;
      }
    }
  }
  cout << "Trouble! Model space not big enough to fill A=" << A << " Z="<< Z << "  emax = " << Nmax << endl;
  return orbitsAZ;

}


void ModelSpace::SetReference(vector<index_t> new_reference)
{
  holes = new_reference;
//  cout << "Setting new reference : ";
//  for (auto& h : holes) cout << h << " ";
//  cout << endl;
  Init(Nmax, holes,core,valence);
}


/*
void ModelSpace::Init_AZ(int nmax, int A, int Z)
{
   target_mass = A;
   target_Z = Z;
   Orbits.clear();
   particles.clear();
   holes.clear();
   valence.clear();
   qspace.clear();
   core.clear();
//   particle_qspace.clear();
//   hole_qspace.clear();
   proton_orbits.clear();
   neutron_orbits.clear();
   OneBodyChannels.clear();
   norbits = (nmax+1)*(nmax+2);
   Orbits.resize(norbits);
   int ncore = 0;
   int zcore = 0;
   for (int N=0; N<=nmax; ++N)
   {
     for (int l=N; l>=0; l-=2)
     {
       int n = (N-l)/2;
       int j2 = 2*l + 1;
       if (zcore < Z)
       {
         AddOrbit(n,l,j2,-1,1,1);
         zcore += j2+1;
       }
       else if (zcore == Z)
       {
         AddOrbit(n,l,j2,-1,0,1);
       }
       else 
       {
         cout << "!!!!!!!!!!!!!!!!!!!\nError: Bad input Z.\n!!!!!!!!!!!!!!!!!! " << endl;
         return;
       }
       if (ncore < A-Z)
       {
         AddOrbit(n,l,j2,1,1,1);
         ncore += j2+1;
       }
       else if (ncore == A-Z)
       {
         AddOrbit(n,l,j2,1,0,1);
       }
       else 
       {
         cout << "!!!!!!!!!!!!!!!!!!!\nError: Bad input N.\n!!!!!!!!!!!!!!!!!! " << endl;
         return;
       }
     }
     for (int l=N%2; l<=N; l+=2)
     {
       int n = (N-l)/2;
       int j2 = 2*l - 1;
       if (j2<0) continue;
       if (zcore < Z)
       {
         AddOrbit(n,l,j2,-1,1,1);
         zcore += j2+1;
       }
       else if (zcore == Z)
       {
         AddOrbit(n,l,j2,-1,0,1);
       }
       else 
       {
         cout << "!!!!!!!!!!!!!!!!!!!\nError: Bad input Z.\n!!!!!!!!!!!!!!!!!! " << endl;
         return;
       }
       if (ncore < A-Z)
       {
         AddOrbit(n,l,j2,1,1,1);
         ncore += j2+1;
       }
       else if (ncore == A-Z)
       {
         AddOrbit(n,l,j2,1,0,1);
       }
       else 
       {
         cout << "Error: Bad input N." << endl;
         cout << "!!!!!!!!!!!!!!!!!!!\nError: Bad input N.\n!!!!!!!!!!!!!!!!!! " << endl;
         return;
       }
     }


   }
   SetupKets();

}
*/


// Some of the more common model spaces, for convenience.
/*
void ModelSpace::Init_He4(int nmax)
{
   vector<string> core = {"p0s1","n0s1"};
   vector<string> valence = {};
   target_mass = 4;
   target_Z = 2;
   Init(nmax,core,valence);
}

void ModelSpace::Init_O16(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
   vector<string> valence = {};
   target_mass = 16;
   target_Z = 8;
   Init(nmax,core,valence);
}

void ModelSpace::Init_Ca40(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   vector<string> valence = {};
   target_mass = 40;
   target_Z = 20;
   Init(nmax,core,valence);
}
*/

/*
void ModelSpace::Init_PShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1"};
   vector<string> valence = {"p0p3","n0p3","p0p1","n0p1"};
   target_mass = 6;
   target_Z = 2;
   Init(nmax,core,valence);
}

void ModelSpace::Init_SDShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
   vector<string> valence = {"p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   target_mass = 18;
   target_Z = 8;
   Init(nmax,core,valence);
}

void ModelSpace::Init_PSDShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1"};
   vector<string> valence = {"p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   target_mass = 6;
   target_Z = 3;
   Init(nmax,core,valence);
}

void ModelSpace::Init_O16PSDShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
   vector<string> valence = {"p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   target_mass = 16;
   target_Z = 8;
   Init(nmax,core,valence);
}

void ModelSpace::Init_FPShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   vector<string> valence = {"p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"};
   target_mass = 42;
   target_Z = 20;
   Init(nmax,core,valence);
}

void ModelSpace::Init_SDFPShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
   vector<string> valence = {"p0d5","n0d5","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"};
   target_mass = 18;
   target_Z = 8;
   Init(nmax,core,valence);
}

void ModelSpace::Init_SD3F7P3Shell(int nmax) // Si28 core - for 40Ca region
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5"};
   vector<string> valence = {"p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p1p3","n1p3"};
   target_mass = 40;
   target_Z = 20;
   Init(nmax,core,valence);
}

void ModelSpace::Init_FPG9Shell(int nmax) // Ni56 core, with g9/2
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7"};
   vector<string> valence = {"p0f5","n0f5","p1p3","n1p3","p1p1","n1p1","p0g9","n0g9"};
   target_mass = 58;
   target_Z = 28;
   Init(nmax,core,valence);
}
*/



ModelSpace ModelSpace::operator=(const ModelSpace& ms)
{
   holes =  ms.holes;
   particles =  ms.particles;
   valence = ms.valence;
   qspace =  ms.qspace;
   core = ms.core;
   proton_orbits =  ms.proton_orbits;
   neutron_orbits =  ms.neutron_orbits;
   KetIndex_pp =  ms.KetIndex_pp;
   KetIndex_ph =  ms.KetIndex_ph;
   KetIndex_hh =  ms.KetIndex_hh;
   KetIndex_cc =  ms.KetIndex_cc;
   KetIndex_vc =  ms.KetIndex_vc;
   KetIndex_qc =  ms.KetIndex_qc;
   KetIndex_vv =  ms.KetIndex_vv;
   KetIndex_qv =  ms.KetIndex_qv;
   KetIndex_qq =  ms.KetIndex_qq;
//   KetIndex_vv =  ms.KetIndex_vv;
//   KetIndex_c_c =  ms.KetIndex_c_c;
//   KetIndex_q_q =  ms.KetIndex_q_q;
//   KetIndex_q_c =  ms.KetIndex_q_c;
//   KetIndex_v_c =  ms.KetIndex_v_c;
//   KetIndex_v_q =  ms. KetIndex_v_q;
   Nmax = ms.Nmax;
   N2max = ms.N2max;
   N3max = ms.N3max;
   OneBodyJmax = ms.OneBodyJmax;
   TwoBodyJmax = ms.TwoBodyJmax;
   ThreeBodyJmax = ms.ThreeBodyJmax;
   OneBodyChannels = ms.OneBodyChannels;
   SortedTwoBodyChannels = ms.SortedTwoBodyChannels;
   SortedTwoBodyChannels_CC = ms.SortedTwoBodyChannels_CC;
   norbits = ms.norbits;
   hbar_omega = ms.hbar_omega;
   target_mass = ms.target_mass;
   target_mass = ms.target_Z;
   Aref = ms.Aref;
   Zref = ms.Zref;
   Orbits = ms.Orbits;
   Kets = ms.Kets;
   TwoBodyChannels = ms.TwoBodyChannels;
   TwoBodyChannels_CC = ms.TwoBodyChannels_CC;
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;

//   cout << "In copy assignment for ModelSpace" << endl;
   return ModelSpace(*this);
}



ModelSpace ModelSpace::operator=(ModelSpace&& ms)
{
   holes =  move(ms.holes);
   particles =  move(ms.particles);
   valence = move(ms.valence);
   qspace =  move(ms.qspace);
   core = move(ms.core);
   proton_orbits =  move(ms.proton_orbits);
   neutron_orbits =  move(ms.neutron_orbits);
   KetIndex_pp =  move(ms.KetIndex_pp);
   KetIndex_ph =  move(ms.KetIndex_ph);
   KetIndex_hh =  move(ms.KetIndex_hh);
   KetIndex_cc =  move(ms.KetIndex_cc);
   KetIndex_vc =  move(ms.KetIndex_vc);
   KetIndex_qc =  move(ms.KetIndex_qc);
   KetIndex_vv =  move(ms.KetIndex_vv);
   KetIndex_qv =  move(ms.KetIndex_qv);
   KetIndex_qq =  move(ms.KetIndex_qq);
//   KetIndex_vv =  move(ms.KetIndex_vv);
//   KetIndex_c_c =  move(ms.KetIndex_c_c);
//   KetIndex_q_q =  move(ms.KetIndex_q_q);
//   KetIndex_q_c =  move(ms.KetIndex_q_c);
//   KetIndex_v_c =  move(ms.KetIndex_v_c);
//   KetIndex_v_q =  move(ms. KetIndex_v_q);
   Nmax = move(ms.Nmax);
   N2max = move(ms.N2max);
   N3max = move(ms.N3max);
   OneBodyJmax = move(ms.OneBodyJmax);
   TwoBodyJmax = move(ms.TwoBodyJmax);
   ThreeBodyJmax = move(ms.ThreeBodyJmax);
   OneBodyChannels = move(ms.OneBodyChannels);
   SortedTwoBodyChannels = move(ms.SortedTwoBodyChannels);
   SortedTwoBodyChannels_CC = move(ms.SortedTwoBodyChannels_CC);
   norbits = move(ms.norbits);
   hbar_omega = move(ms.hbar_omega);
   target_mass = move(ms.target_mass);
   target_Z = move(ms.target_Z);
   Aref = move(ms.Aref);
   Zref = move(ms.Zref);
   Orbits = move(ms.Orbits);
   Kets = move(ms.Kets);
   TwoBodyChannels = move(ms.TwoBodyChannels);
   TwoBodyChannels_CC = move(ms.TwoBodyChannels_CC);
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
   for (TwoBodyChannel& tbc : ms.TwoBodyChannels)   tbc.modelspace = NULL;
   for (TwoBodyChannel_CC& tbc_cc : ms.TwoBodyChannels_CC)   tbc_cc.modelspace = NULL;
//   cout << "In move assingment for ModelSpace" << endl;
   return ModelSpace(*this);
}



void ModelSpace::AddOrbit(Orbit orb)
{
//  AddOrbit(orb.n, orb.l, orb.j2, orb.tz2, orb.ph, orb.io);
  AddOrbit(orb.n, orb.l, orb.j2, orb.tz2, orb.ph, orb.cvq);
}

//void ModelSpace::AddOrbit(int n, int l, int j2, int tz2, int ph, int io)
void ModelSpace::AddOrbit(int n, int l, int j2, int tz2, int ph, int cvq)
{
   index_t ind = Index1(n, l, j2, tz2);
//   Orbits[ind] = Orbit(n,l,j2,tz2,ph,io,ind);
   Orbits[ind] = Orbit(n,l,j2,tz2,ph,cvq,ind);

   if (j2 > OneBodyJmax)
   {
      OneBodyJmax = j2;
      TwoBodyJmax = OneBodyJmax;
      ThreeBodyJmax = OneBodyJmax*3-1;
      nTwoBodyChannels = 2*3*(TwoBodyJmax+1);
   }

   if (ph == 0) particles.push_back(ind);
   if (ph == 1) holes.push_back(ind);
//   if (io == 0) valence.push_back(ind);
   if (cvq == 0) core.push_back(ind);
   if (cvq == 1) valence.push_back(ind);
   if (cvq == 2) qspace.push_back(ind);
//   if (io == 1)
//   {
//     qspace.push_back(ind);
//     if (ph == 0) particle_qspace.push_back(ind);
//     if (ph == 1) hole_qspace.push_back(ind);
//   }
   if (tz2<0) proton_orbits.push_back(ind);
   if (tz2>0) neutron_orbits.push_back(ind);

   OneBodyChannels[{l, j2, tz2}].push_back(ind);
}



int ModelSpace::GetOrbitIndex(string orb)
{
  vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};
  int n=-1,l=-1,j2=-1;
  int tz2 = orb[0]=='p' ? -1 : 1;
  stringstream(orb.substr(1,1)) >> n;
  auto it_l = find(l_list.begin(), l_list.end(), orb[2]);
  if ( it_l != l_list.end() )
    l = it_l - l_list.begin();
  else
    cout << "Bad orbit label " << orb << endl;
  stringstream(orb.substr(3)) >> j2;
  return Index1(n,l,j2,tz2);
}

int ModelSpace::GetTwoBodyChannelIndex(int j, int p, int t)
{
   return (t+1)*2*(TwoBodyJmax+1) + p*(TwoBodyJmax+1) + j;
}


void ModelSpace::SetupKets()
{
   int index = 0;

   Kets.resize(Index2(norbits-1,norbits-1)+1);
   for (int p=0;p<norbits;p++)
   {
     for (int q=p;q<norbits;q++)
     {
        index = Index2(p,q);
        Kets[index] = Ket(GetOrbit(p),GetOrbit(q));
     }
   }

  for (index_t index=0;index<Kets.size();++index)
  {
    Ket& ket = Kets[index];
    int Tz = (ket.op->tz2 + ket.oq->tz2)/2;
    int parity = (ket.op->l + ket.oq->l)%2;
//    MonopoleKets[Tz+1][parity].push_back(index);
    MonopoleKets[Tz+1][parity][index]=MonopoleKets[Tz+1][parity].size()-1;
    int php = ket.op->ph;
    int phq = ket.oq->ph;
//    int iop = ket.op->io;
//    int ioq = ket.oq->io;
    int cvq_p = ket.op->cvq;
    int cvq_q = ket.oq->cvq;
    if (php+phq==0) KetIndex_pp.push_back(index);
    if (php+phq==1) KetIndex_ph.push_back(index);
    if (php+phq==2) KetIndex_hh.push_back(index);
    if (cvq_p+cvq_q==0)      KetIndex_cc.push_back(index); // 00
    if (cvq_p+cvq_q==1)      KetIndex_vc.push_back(index); // 01
    if (abs(cvq_p-cvq_q)==2) KetIndex_qc.push_back(index); // 02
    if (cvq_p*cvq_q==1)      KetIndex_vv.push_back(index); // 11
    if (cvq_p+cvq_q==3)      KetIndex_qv.push_back(index); // 12
    if (cvq_p+cvq_q==4)      KetIndex_qq.push_back(index); // 22
//     if (( php + phq)==2) // hh
//     {
//        KetIndex_hh.push_back(index);
//        if ((iop+ioq)==2) // qq
//        {
//           KetIndex_c_c.push_back(index);
//        }
//     }
//     else if ((php + phq) == 0) // pp
//     {
//        KetIndex_pp.push_back(index);
//        if ((iop+ioq)==2) // qq
//        {
//           KetIndex_q_q.push_back(index);
//        }
//     }
//     else //ph
//     {
//        KetIndex_ph.push_back(index);
//        if ((iop+ioq)==2) // qq
//        {
//           KetIndex_q_c.push_back(index);
//        }
//     }
//     if ((iop + ioq) == 0) // vv
//     {
//        KetIndex_vv.push_back(index);
//     }
//
//     if ((iop + ioq) == 1) // vq
//     {
//       if ((iop + php == 2) or (ioq+phq==2) ) // the qspace orbit is a hole
//          KetIndex_v_c.push_back(index);
//       else // v particle_q
//          KetIndex_v_q.push_back(index);
//     }

   }
//   cout << "kets_cc: ";
//   for (auto& k : KetIndex_cc ) cout << k << "->(" << GetKet(k).p << "," << GetKet(k).q << ") ";
//   cout << endl;
//   cout << "kets_vc: ";
//   for (auto& k : KetIndex_vc ) cout << k << "->(" << GetKet(k).p << "," << GetKet(k).q << ") ";
//   cout << endl;
//   cout << "kets_qc: ";
//   for (auto& k : KetIndex_qc ) cout << k << "->(" << GetKet(k).p << "," << GetKet(k).q << ") ";
//   cout << endl;
//   cout << "kets_vv: ";
//   for (auto& k : KetIndex_vv ) cout << k << "->(" << GetKet(k).p << "," << GetKet(k).q << ") ";
//   cout << endl;
//   cout << "kets_qv: ";
//   for (auto& k : KetIndex_qv ) cout << k << "->(" << GetKet(k).p << "," << GetKet(k).q << ") ";
//   cout << endl;
//   cout << "kets_qq: ";
//   for (auto& k : KetIndex_qq ) cout << k << "->(" << GetKet(k).p << "," << GetKet(k).q << ") ";
//   cout << endl;

   SortedTwoBodyChannels.resize(nTwoBodyChannels);
   SortedTwoBodyChannels_CC.resize(nTwoBodyChannels);
   for (int ch=0;ch<nTwoBodyChannels;++ch)
   {
      TwoBodyChannels.push_back(move(TwoBodyChannel(ch,this)));
      TwoBodyChannels_CC.push_back(move(TwoBodyChannel_CC(ch,this)));
      SortedTwoBodyChannels[ch] = ch;
      SortedTwoBodyChannels_CC[ch] = ch;
   }
   // Sort the two body channels in descending order of matrix dimension and discard the size-0 ones.
   // Hopefully this can help with load balancing.
   sort(SortedTwoBodyChannels.begin(),SortedTwoBodyChannels.end(),[this](int i, int j){ return TwoBodyChannels[i].GetNumberKets() > TwoBodyChannels[j].GetNumberKets(); }  );
   sort(SortedTwoBodyChannels_CC.begin(),SortedTwoBodyChannels_CC.end(),[this](int i, int j){ return TwoBodyChannels_CC[i].GetNumberKets() > TwoBodyChannels_CC[j].GetNumberKets(); }  );
   while (  TwoBodyChannels[ SortedTwoBodyChannels.back() ].GetNumberKets() <1 ) SortedTwoBodyChannels.pop_back();
   while (  TwoBodyChannels_CC[ SortedTwoBodyChannels_CC.back() ].GetNumberKets() <1 ) SortedTwoBodyChannels_CC.pop_back();
}




double ModelSpace::GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3)
{
// { j1 j2 j3 }
// { J1 J2 J3 }

//   unsigned long long int key = 20000000000*j1 + 200000000*j2 + 2000000*j3 + 20000*J1 + 200*J2 + 2*J3;
   unsigned long int key = (((unsigned long int) (2*j1)) << 30) +
                           (((unsigned long int) (2*j2)) << 24) +
                           (((unsigned long int) (2*j3)) << 18) +
                           (((unsigned long int) (2*J1)) << 12) +
                           (((unsigned long int) (2*J2)) <<  6) +
                            ((unsigned long int) (2*J3));

   auto it = SixJList.find(key);
   if (it != SixJList.end() ) return it->second;
   double sixj = AngMom::SixJ(j1,j2,j3,J1,J2,J3);
   #pragma omp critical
   SixJList[key] = sixj;
   return sixj;
}



void ModelSpace::PreCalculateMoshinsky()
{
//  if ( not MoshList.empty() ) return; // Already done calculated it...
  #pragma omp parallel for schedule(dynamic,1)
  for (int N=0; N<=N2max/2; ++N)
  {
   unordered_map<unsigned long long int,double> local_MoshList;
   for (int n=0; n<=min(N,N2max/2-N); ++n)
   {
    for (int Lam=0; Lam<=N2max-2*N-2*n; ++Lam)
    {
     int lam_max = (N==n ? min(Lam,N2max-2*N-2*n-Lam) : N2max-2*N-2*n-Lam);
     for (int lam=0; lam<=lam_max; ++lam)
     {
      int e2 = 2*N+Lam + 2*n+lam;
      for (int L=abs(Lam-lam); L<=Lam+lam; ++L)
      {
       for (int n1=0; n1<=N; ++n1)
       {
        for (int n2=0; n2<=min(n1,e2/2-n1); ++n2)
        {
         int l1max = n1==N? min(Lam,e2-2*n1-2*n2) : e2-2*n1-2*n2;
         for (int l1=0; l1<=l1max; ++l1 )
         {
          int l2 = e2-2*n1-2*n2-l1;
          if ( (l1+l2+lam+Lam)%2 >0 ) continue;
          if ( l2<abs(L-l1) or l2>L+l1 ) continue;
          // nmax = 16, lmax = 32 -> good up to emax=32, which I'm nowhere near.
          unsigned long long int key =   ((unsigned long long int) N   << 40)
                                       + ((unsigned long long int) Lam << 34)
                                       + ((unsigned long long int) n   << 30)
                                       + ((unsigned long long int) lam << 26)
                                       + ((unsigned long long int) n1  << 22)
                                       + ((unsigned long long int) l1  << 16)
                                       + ((unsigned long long int) n2  << 12)
                                       + ((unsigned long long int) l2  << 6 )
                                       +  L;

          double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
          local_MoshList[key] = mosh;
         }
        }
       }
      }
     }
    }
   }
   #pragma omp critical
   MoshList.insert( local_MoshList.begin(), local_MoshList.end() );
  }
}



double ModelSpace::GetMoshinsky( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L)
{
   int phase_mosh = 1;
  int switches = 10;

  while (switches > 0)
  {
   switches = 0;
   if (n2>n1 or (n2==n1 and l2>l1))
   {
      swap(n1,n2);
      swap(l1,l2);
      phase_mosh *= phase(Lam+L);
      ++switches;
   }
   if (n>N or (n==N and lam>Lam))
   {
      swap(n,N);
      swap(lam,Lam);
      phase_mosh *= phase(l1 +L);
      ++switches;
   }

   if (n1>N or (n1==N and l1>Lam) or (n1==N and l1==Lam and n2>n) or (n1==N and l1==Lam and n2==n and l2>lam) )
   {
      swap(n1,N);
      swap(l1,Lam);
      swap(n2,n);
      swap(l2,lam);
      ++switches;
//      phase_mosh *= phase(l2+lam); // This phase is given in Moshinsky and Brody, but with the current algorithm, it appears not to be required.
   }
  }

          unsigned long long int key =   ((unsigned long long int) N   << 40)
                                       + ((unsigned long long int) Lam << 34)
                                       + ((unsigned long long int) n   << 30)
                                       + ((unsigned long long int) lam << 26)
                                       + ((unsigned long long int) n1  << 22)
                                       + ((unsigned long long int) l1  << 16)
                                       + ((unsigned long long int) n2  << 12)
                                       + ((unsigned long long int) l2  << 6 )
                                       +  L;

//   unsigned long long int key =  1000000000000 * N
//                                + 100000000000 * Lam
//                                +   1000000000 * n
//                                +    100000000 * lam
//                                +      1000000 * n1
//                                +       100000 * l1
//                                +         1000 * n2
//                                +          100 * l2
//                                +                 L;
   auto it = MoshList.find(key);
   if ( it != MoshList.end() )  return it->second * phase_mosh;

   // if we didn't find it, we need to calculate it.
   double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
//   cout << "Shouldn't be here..." << N << " " << Lam << " " <<  n << " " << lam << " " << n1 << " " << l1 << " " << n2 << " " << l2 << " " << L << endl;
//   #pragma omp atomic
   MoshList[key] = mosh;
   return mosh * phase_mosh;

}




double ModelSpace::GetNineJ(double j1, double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
{
   int k1 = 2*j1;
   int k2 = 2*j2;
   int K12 = 2*J12;
   int k3 = 2*j3;
   int k4 = 2*j4;
   int K34 = 2*J34;
   int K13 = 2*J13;
   int K24 = 2*J24;
   int K = 2*J;

   array<int,9> klist = {k1,k2,K12,k3,k4,K34,K13,K24,K};
   array<double,9> jlist = {j1,j2,J12,j3,j4,J34,J13,J24,J};
   int imin = min_element(klist.begin(),klist.end()) - klist.begin();
   switch (imin)
   {
      case 0:
       klist = {k4,K34,k3,K24,K,K13,k2,K12,k1};
       jlist = {j4,J34,j3,J24,J,J13,j2,J12,j1};
       break;
      case 1:
       klist = {K13,K,K24,k3,K34,k4,k1,K12,k2};
       jlist = {J13,J,J24,j3,J34,j4,j1,J12,j2};
       break;
      case 2:
       klist = {k3,k4,K34,K13,K24,K,k1,k2,K12};
       jlist = {j3,j4,J34,J13,J24,J,j1,j2,J12};
       break;
      case 3:
       klist = {K12,k2,k1,K,K24,K13,K34,k4,k3};
       jlist = {J12,j2,j1,J,J24,J13,J34,j4,j3};
       break;
      case 4:
       klist = {k1,K12,k2,K13,K,K24,k3,K34,k4};
       jlist = {j1,J12,j2,J13,J,J24,j3,J34,j4};
       break;
      case 5:
       klist = {K13,K24,K,k1,k2,K12,k3,k4,K34};
       jlist = {J13,J24,J,j1,j2,J12,j3,j4,J34};
       break;
      case 6:
       klist = {k2,K12,k1,k4,K34,k3,K24,K,K13};
       jlist = {j2,J12,j1,j4,J34,j3,J24,J,J13};
       break;
      case 7:
       klist = {K12,k1,k2,K34,k3,k4,K,K13,K24};
       jlist = {J12,j1,j2,J34,j3,j4,J,J13,J24};
       break;
      case 8:
       break;
   }

   unsigned long long int key =   klist[0];
   unsigned long long int factor = 100;
   for (int i=1; i<9; ++i)
   {
      key += klist[i]*factor;
      factor *=100;
   }
   auto it = NineJList.find(key);
   if (it != NineJList.end() ) return it->second;
   double ninej = AngMom::NineJ(jlist[0],jlist[1],jlist[2],jlist[3],jlist[4],jlist[5],jlist[6],jlist[7],jlist[8]);
   #pragma omp critical
   NineJList[key] = ninej;
   return ninej;

}


