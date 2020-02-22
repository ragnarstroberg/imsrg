#include "ModelSpace.hh"
#include "AngMom.hh"
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include "omp.h"
#include <cstdlib> // for EXIT_FAILURE

double ModelSpace::OCC_CUT = 1e-6;
int ModelSpace::NOT_AN_ORBIT = 99999;
Orbit ModelSpace::NULL_ORBIT = Orbit();



// Static members

std::unordered_map<uint64_t,double> ModelSpace::SixJList;
std::unordered_map<uint64_t,double> ModelSpace::NineJList;
std::unordered_map<uint64_t,double> ModelSpace::MoshList;
std::map< std::string, std::vector<std::string> > ModelSpace::ValenceSpaces  {
{ "s-shell"  ,         {"vacuum", "p0s1","n0s1"}},
{ "p-shell"  ,         {"He4", "p0p3","n0p3","p0p1","n0p1"}},
{ "sp-shell"  ,        {"vacuum", "p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"}},
{ "spsd-shell"  ,      {"vacuum", "p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"}},
{ "spsdpf-shell"  ,    {"vacuum", "p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"}},
{ "sd-shell"  ,        {"O16", "p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"}},
{ "psd-shell"  ,       {"He4", "p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"}},
{ "psdNR-shell"  ,     {"He10","p0p3","p0p1","n0d5","n0d3","n1s1"}}, // protons in p shell, neutrons in sd shell (NR is for neutron-rich)
{ "psdPR-shell"  ,     {"O10","n0p3","n0p1","p0d5","p0d3","p1s1"}}, // neutrons in p shell, protnons in sd shell (PR is for proton-rich)
{ "fp-shell"  ,        {"Ca40","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"}},
{ "sdfp-shell"  ,      {"O16", "p0d5","n0d5","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"}},
{ "sdfpNR-shell"  ,    {"O28", "p0d5","p0d3","p1s1","n0f7","n0f5","n1p3","n1p1"}}, // protons in sd shell, neutrons in fp shell, a la SDPFU from Nowacki/Poves (NR is for neutron-rich)
{ "sdfpPR-shell"  ,    {"Ca28", "n0d5","n0d3","n1s1","p0f7","p0f5","p1p3","p1p1"}}, // neutrons in sd shell, protons in fp shell, (PR is for proton-rich)
{ "fpg9-shell"  ,      {"Ca40","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1","p0g9","n0g9"}},
{ "fpg9NR-shell"  ,    {"Ca40","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1","n0g9"}}, // just add g9/2 for neutrons
{ "fpgdsNR-shell"  ,   {"Ca60","p0f7","p0f5","p1p3","p1p1","n0g9","n0g7","n1d5","n1d3","n2s1"}}, // protons in the fp shell, neutrons in the gds shell
{ "sd3f7p3-shell"  ,   {"Si28","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p1p3","n1p3"}},
{ "gds-shell" ,        {"Zr80","p0g9","n0g9","p0g7","n0g7","p1d5","n1d5","p1d3","n1d3","p2s1","n2s1"}}, // This is a big valence space, more than a few particles will be a serious shell model diagonalization
};



//ModelSpace::~ModelSpace()
//{
////  std::cout << "In ModelSpace destructor. emax = " << Emax << std::endl;
//}

ModelSpace::ModelSpace()
:  Emax(0), E2max(0), E3max(0), Lmax(0), Lmax2(0), Lmax3(0), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0), EmaxUnocc(0), dE3max(999), norbits(0),
  hbar_omega(20), target_mass(16), sixj_has_been_precalculated(false), moshinsky_has_been_precalculated(false),
  scalar_transform_first_pass(true),scalar3b_transform_first_pass(true), tensor_transform_first_pass(40,true), single_species(false)
{
   SetUpOrbits();
//  std::cout << "In default constructor" << std::endl;
}






ModelSpace::ModelSpace(const ModelSpace& ms)
 :
   Emax(ms.Emax), E2max(ms.E2max), E3max(ms.E3max), Lmax(ms.Lmax), Lmax2(ms.Lmax2), Lmax3(ms.Lmax3),
   OneBodyJmax(ms.OneBodyJmax), TwoBodyJmax(ms.TwoBodyJmax), ThreeBodyJmax(ms.ThreeBodyJmax), EmaxUnocc(ms.EmaxUnocc),
   dE3max(ms.dE3max), e_fermi(ms.e_fermi), norbits(ms.norbits), hbar_omega(ms.hbar_omega),
   target_mass(ms.target_mass), target_Z(ms.target_Z), Aref(ms.Aref), Zref(ms.Zref),
   Orbits(ms.Orbits), Kets(ms.Kets), Kets3(ms.Kets3),
   OneBodyChannels(ms.OneBodyChannels), nTwoBodyChannels(ms.nTwoBodyChannels),
   TwoBodyChannels(ms.TwoBodyChannels), TwoBodyChannels_CC(ms.TwoBodyChannels_CC),
   nThreeBodyChannels(ms.nThreeBodyChannels), ThreeBodyChannels(ms.ThreeBodyChannels),
   Ket3IndexLookup(ms.Ket3IndexLookup), ThreeBodyChannelLookup(ms.ThreeBodyChannelLookup),
   PandyaLookup(ms.PandyaLookup), OrbitLookup(ms.OrbitLookup),
   holes( ms.holes), particles( ms.particles),
   core(ms.core), valence(ms.valence), qspace( ms.qspace), 
   proton_orbits( ms.proton_orbits),neutron_orbits( ms.neutron_orbits), all_orbits(ms.all_orbits),
   KetIndex_pp( ms.KetIndex_pp), KetIndex_ph( ms.KetIndex_ph), KetIndex_hh( ms.KetIndex_hh),
   KetIndex_cc( ms.KetIndex_cc),
   KetIndex_vc( ms.KetIndex_vc),
   KetIndex_qc( ms.KetIndex_qc),
   KetIndex_vv( ms.KetIndex_vv),
   KetIndex_qv( ms.KetIndex_qv),
   KetIndex_qq( ms.KetIndex_qq),
   Ket_occ_hh( ms.Ket_occ_hh),
   Ket_unocc_hh( ms.Ket_unocc_hh),
   SortedTwoBodyChannels(ms.SortedTwoBodyChannels),
   SortedTwoBodyChannels_CC(ms.SortedTwoBodyChannels_CC),
   sixj_has_been_precalculated(ms.sixj_has_been_precalculated),
   moshinsky_has_been_precalculated(ms.moshinsky_has_been_precalculated),
   scalar_transform_first_pass(true), scalar3b_transform_first_pass(true), tensor_transform_first_pass(40,true), single_species(ms.single_species)
{
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
}

ModelSpace::ModelSpace(ModelSpace&& ms)
 :
   Emax(ms.Emax), E2max(ms.E2max), E3max(ms.E3max), Lmax(ms.Lmax), Lmax2(ms.Lmax2), Lmax3(ms.Lmax3),
   OneBodyJmax(ms.OneBodyJmax), TwoBodyJmax(ms.TwoBodyJmax), ThreeBodyJmax(ms.ThreeBodyJmax), EmaxUnocc(ms.EmaxUnocc),
   dE3max(ms.dE3max), e_fermi(ms.e_fermi), norbits(ms.norbits), hbar_omega(ms.hbar_omega),
   target_mass(ms.target_mass), target_Z(ms.target_Z), Aref(ms.Aref), Zref(ms.Zref),
   Orbits(std::move(ms.Orbits)), Kets(std::move(ms.Kets)),
   OneBodyChannels(std::move(ms.OneBodyChannels)), nTwoBodyChannels(std::move(ms.nTwoBodyChannels)),
   TwoBodyChannels(std::move(ms.TwoBodyChannels)), TwoBodyChannels_CC(std::move(ms.TwoBodyChannels_CC)),
   nThreeBodyChannels(std::move(ms.nThreeBodyChannels)), ThreeBodyChannels(std::move(ms.ThreeBodyChannels)),
   Ket3IndexLookup(std::move(ms.Ket3IndexLookup)), ThreeBodyChannelLookup(std::move(ThreeBodyChannelLookup)),
   PandyaLookup(std::move(ms.PandyaLookup)), OrbitLookup(ms.OrbitLookup),
   holes( std::move(ms.holes)), particles( std::move(ms.particles)),
   core(std::move(ms.core)), valence(std::move(ms.valence)),  qspace( std::move(ms.qspace)),  
   proton_orbits( std::move(ms.proton_orbits)), neutron_orbits( std::move(ms.neutron_orbits)), all_orbits( std::move(ms.all_orbits)),
   KetIndex_pp( std::move(ms.KetIndex_pp)), KetIndex_ph( std::move(ms.KetIndex_ph)), KetIndex_hh( std::move(ms.KetIndex_hh)),
   KetIndex_cc( ms.KetIndex_cc),
   KetIndex_vc( ms.KetIndex_vc),
   KetIndex_qc( ms.KetIndex_qc),
   KetIndex_vv( ms.KetIndex_vv),
   KetIndex_qv( ms.KetIndex_qv),
   KetIndex_qq( ms.KetIndex_qq),
   Ket_occ_hh( ms.Ket_occ_hh),
   Ket_unocc_hh( ms.Ket_unocc_hh),
   SortedTwoBodyChannels(std::move(ms.SortedTwoBodyChannels)),
   SortedTwoBodyChannels_CC(std::move(ms.SortedTwoBodyChannels_CC)),
   sixj_has_been_precalculated(ms.sixj_has_been_precalculated),
   moshinsky_has_been_precalculated(ms.moshinsky_has_been_precalculated),
   scalar_transform_first_pass(true),scalar3b_transform_first_pass(true), tensor_transform_first_pass(40,true), single_species(ms.single_species)
{
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
   for (TwoBodyChannel& tbc : ms.TwoBodyChannels)   tbc.modelspace = NULL;
   for (TwoBodyChannel_CC& tbc_cc : ms.TwoBodyChannels_CC)   tbc_cc.modelspace = NULL;
}


// orbit std::string representation is e.g. p0f7
// Assumes that the core is hole states that aren't in the valence space.
ModelSpace::ModelSpace(int emax, std::vector<std::string> hole_list, std::vector<std::string> valence_list)
:  Emax(emax), E2max(2*emax), E3max(std::min(14,3*emax)), Lmax(emax), Lmax2(emax), Lmax3(emax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0), EmaxUnocc(emax), dE3max(999), norbits(0), hbar_omega(20), target_mass(16),
     moshinsky_has_been_precalculated(false), scalar_transform_first_pass(true), scalar3b_transform_first_pass(true), tensor_transform_first_pass(40,true), single_species(false)
{
   SetUpOrbits();
   Init(emax, hole_list, hole_list, valence_list); 
}

// If we don't want the reference to be the core
ModelSpace::ModelSpace(int emax, std::vector<std::string> hole_list, std::vector<std::string> core_list, std::vector<std::string> valence_list)
: Emax(emax), E2max(2*emax), E3max(std::min(14,3*emax)), Lmax(emax), Lmax2(emax), Lmax3(emax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0), EmaxUnocc(emax), dE3max(999), norbits(0), hbar_omega(20), target_mass(16),
     sixj_has_been_precalculated(false),moshinsky_has_been_precalculated(false), scalar_transform_first_pass(true), scalar3b_transform_first_pass(true),tensor_transform_first_pass(40,true),single_species(false)
{
   SetUpOrbits();
   Init(emax, hole_list, core_list, valence_list); 
}

// Most conventient interface
ModelSpace::ModelSpace(int emax, std::string reference, std::string valence)
: Emax(emax), E2max(2*emax), E3max(std::min(14,3*emax)), Lmax(emax), Lmax2(emax), Lmax3(emax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0), EmaxUnocc(emax), dE3max(999), hbar_omega(20),
     sixj_has_been_precalculated(false),moshinsky_has_been_precalculated(false), scalar_transform_first_pass(true), scalar3b_transform_first_pass(true),tensor_transform_first_pass(40,true),single_species(false)
{
  SetUpOrbits();
  Init(emax,reference,valence);
}

ModelSpace::ModelSpace(int emax, std::string valence)
: Emax(emax), E2max(2*emax), E3max(std::min(14,3*emax)), Lmax(emax), Lmax2(emax), Lmax3(emax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0), EmaxUnocc(emax), dE3max(999), hbar_omega(20),
     sixj_has_been_precalculated(false),moshinsky_has_been_precalculated(false), scalar_transform_first_pass(true), scalar3b_transform_first_pass(true), tensor_transform_first_pass(40,true),single_species(false)
{
  auto itval = ValenceSpaces.find(valence);
  if ( itval != ValenceSpaces.end() ) // we've got a valence space
     Init(emax,itval->second[0],valence);
  else  // no valence space. we've got a single-reference.
     Init(emax,valence,valence);
}



void ModelSpace::SetUpOrbits( )
{
   std::vector<int> IsospinList = {-1, 1};
   if (single_species) IsospinList.pop_back(); // just use tz = -1 in this case

   norbits = (Emax+1)*(Emax+2);
   Orbits.resize(0);
//   std::cout << "Starting loop in SetUpOrbits " << std::endl;
   for (int N=0; N<=Emax; ++N)
   {
     for (int l=N; l>=0; l-=2)
     {
       if (l>Lmax) continue;
       int n = (N-l)/2;
       for (int j2=2*l+1; j2>=2*l-1 and j2>0; j2-=2)
       {
//         for (int tz : {-1, 1} )
         for (int tz : IsospinList )
         {
            double occ = 0;
            int cvq = 2;
//            int indx = Index1(n,l,j2,tz);
//            if (hole_list.find(indx) != hole_list.end()) occ = hole_list[indx];
//            if ( find(core_list.begin(), core_list.end(), indx) != core_list.end() ) cvq=0; // core orbit
//            if ( find(valence_list.begin(), valence_list.end(), indx) != valence_list.end() ) cvq=1; // valence orbit
//            std::cout << "calling AddOrbit( " << n << " " << l << " " << tz << " " << occ << " " << cvq << " ) " << std::endl;
            AddOrbit(n,l,j2,tz,occ,cvq);
         }
       }
     }
   }
   norbits = all_orbits.size();
//   Orbits.resize(norbits);
//   std::cout << "done with SetUpOrbits" << std::endl;
}


// Specify the reference and either the core or valence
// This is the most convenient interface
void ModelSpace::Init(int emax, std::string reference, std::string valence)
{
//  int Aref,Zref;
  GetAZfromString(reference,Aref,Zref);
  std::map<index_t,double> hole_list = GetOrbitsAZ(Aref,Zref);
  Init(emax,hole_list,valence);
}

void ModelSpace::Init(int emax, std::map<index_t,double> hole_list, std::string valence)
{
  int Ac,Zc;
//  std::vector<index_t> valence_list, core_list;
  std::set<index_t> valence_list, core_list;

  if (valence == "0hw-shell")
  {
    Get0hwSpace(Aref,Zref,core_list,valence_list);
  }
  else if ( valence.find(",")!=std::string::npos ) // interpet as a comma-separated list of core followed by valence orbits
  {
    ParseCommaSeparatedValenceSpace(valence,core_list,valence_list);
  }
  else if ( valence.find("FCI")!=std::string::npos ) // FCI space, so no core, all orbits are valence.
  {
    index_t num_orbits = (emax+1)*(emax+2);
    for (index_t i=0;i<num_orbits;++i) valence_list.insert( i );
//    for (index_t i=0;i<num_orbits;++i) valence_list.push_back( i );
  }
  else // check if it's one of the pre-defined spaces
  {
    auto itval = ValenceSpaces.find(valence);
    std::string core_str;
  
    if ( itval != ValenceSpaces.end() ) // we've got a valence space
    {
       core_str = itval->second[0];
       std::vector<index_t> vvec = String2Index(std::vector<std::string>(itval->second.begin()+1,itval->second.end()));
       for (auto v : vvec ) valence_list.insert(v);
//       valence_list = String2Index(std::vector<std::string>(itval->second.begin()+1,itval->second.end()));
    }
    else  // no valence space. we've got a single-reference.
    {
       core_str = valence;
    }
  
    GetAZfromString(core_str,Ac,Zc);
    for (auto& it_core : GetOrbitsAZ(Ac,Zc) ) core_list.insert(it_core.first);
//    for (auto& it_core : GetOrbitsAZ(Ac,Zc) ) core_list.push_back(it_core.first);
  }

  // This should be done in Init, but it doesn't hurt to do it here as well, I think.
  SetTargetMass(Aref);
  SetTargetZ(Zref);
//  target_mass = Aref;
//  target_Z = Zref;
  Init(emax,hole_list,core_list,valence_list);
  
}


// Specify the model space with std::strings of orbit lists.
// Less convenient, but more flexible
void ModelSpace::Init(int emax, std::vector<std::string> hole_list, std::vector<std::string> core_list, std::vector<std::string> valence_list)
{
   std::cout << "Creating a model space with Emax = " << Emax << "  and hole orbits [";
   for (auto& h : hole_list)  std::cout << h << " ";
   std::cout << "]   and core orbits [";
   for (auto& c : core_list)    std::cout << c << " ";
   std::cout << "]   and valence orbits [";
   for (auto& v : valence_list)   std::cout << v << " ";
   std::cout << "]" << std::endl;
   std::map<index_t,double> hole_map;
   for (auto& h : String2Index(hole_list)) hole_map[h] = 1.0;
  Init(emax, hole_map, String2Index(core_list), String2Index(valence_list) );
}


void ModelSpace::Init_occ_from_file(int emax, std::string valence, std::string occ_file)
{
  index_t orb;
  double occ;
  std::map<index_t,double> hole_list;

  std::ifstream infile(occ_file);
  if (!infile.good())
  {
    std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "Trouble reading file: " << occ_file << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl << std::endl;
  }

  while( infile >> orb >> occ )
  {
    if ( hole_list.find(orb) != hole_list.end() and  std::abs( hole_list[orb] -occ) > 1e-6) 
    {
        std::cout << "Warning: in file " << occ_file << ", redefinition of occupation of orbit "
             << orb << "  " << hole_list[orb] << " => " << occ << std::endl;
    }

    std::cout << "from occ file: " << std::endl;
    hole_list[orb] = occ;
    std::cout << orb << " " << occ << std::endl;
  }

  Init(emax,hole_list,valence);
}

// An attempt at backward compatibility
void ModelSpace::Init(int emax, std::map<index_t,double> hole_list, std::vector<index_t> core_list, std::vector<index_t> valence_list)
{
  std::set<index_t> clist(core_list.begin(),core_list.end());
  std::set<index_t> vlist(valence_list.begin(),valence_list.end());
  Init(emax,hole_list,clist,vlist);
}


// This is the Init which should inevitably be called
void ModelSpace::Init(int emax, std::map<index_t,double> hole_list, std::set<index_t> core_list, std::set<index_t> valence_list)
{
   ClearVectors();
   Emax = emax;
   std::cout << "core list: ";
   for (auto& c : core_list) std::cout << c << " ";
   std::cout << std::endl;
   std::cout << "valence list: ";
   for (auto& v : valence_list) std::cout << v << " ";
   std::cout << std::endl;
   std::cout << "hole list: ";
   for (auto& h : hole_list) std::cout << h.first << " ( " << h.second << " ) ";
   std::cout << std::endl;

   // Make sure no orbits are both core and valence
   for (auto& c : core_list)
   {
     if ( find(valence_list.begin(), valence_list.end(), c) != valence_list.end() )
       std::cout << "!!!!!!!!!!!!! ModelSpace::Init : Conflicting definition. Orbit " << c << " is in core and valence spaces." << std::endl;
   }

   std::vector<int> IsospinList = {-1, 1};
   if (single_species) IsospinList.pop_back(); // just use tz = -1 in this case

   Orbits.resize(0);
   OrbitLookup.clear(); // if we've resized Orbits to zero, there's nothing to point to or look up...

   for (int N=0; N<=Emax; ++N)
   {
     for (int l=N; l>=0; l-=2)
     {
       if (l>Lmax) continue;
       int n = (N-l)/2;
       for (int j2=2*l+1; j2>=2*l-1 and j2>0; j2-=2)
       {
         for (int tz : IsospinList )
         {
            // check if this is a non-hole-type orbit and if we're above the EmaxUnocc cut.

            if ( ((2*n+l)>EmaxUnocc) and (hole_quantum_numbers.find({l,j2})==hole_quantum_numbers.end()) ) continue;
            double occ = 0;
            int cvq = 2;
            AddOrbit(n,l,j2,tz,occ,cvq);  // First, add the orbit to make sure it's in the lookup tables
            int indx = Index1(n,l,j2,tz);
            if (hole_list.find(indx) != hole_list.end()) occ = hole_list[indx];
            if ( find(core_list.begin(), core_list.end(), indx) != core_list.end() ) cvq=0; // core orbit
            if ( find(valence_list.begin(), valence_list.end(), indx) != valence_list.end() ) cvq=1; // valence orbit
            AddOrbit(n,l,j2,tz,occ,cvq); // Add it again with the right occupation and classification. (There's probably a cleaner way to do this.)
         }
       }
     }
   }
   norbits = all_orbits.size();
   double atmp=0;
   double ztmp=0;
   for (auto& h : holes)
   {
     Orbit& oh = GetOrbit(h);
     atmp += (oh.j2+1.)*oh.occ;
     if (oh.tz2 < 0) ztmp += (oh.j2+1.)*oh.occ;
   }
   Aref = round(atmp);
   Zref = round(ztmp);
   if (std::abs(Aref-atmp)>1e-5 or std::abs(Zref-ztmp)>1e-5)
   {
     std::cout << std::endl << "!!!! WARNING  " << __func__ << " recomputed A,Z and got " << atmp << " " << ztmp << std::endl;
   }
   Aref = round(atmp);
   Zref = round(ztmp);
   FindEFermi();
   SetTargetMass(Aref);
   SetTargetZ(Zref);
   SetupKets();
   Setup3bKets();
}


void ModelSpace::InitSingleSpecies(int emax, std::string reference, std::string valence)
{
  single_species = true;
  ClearVectors();
  GetAZfromString(reference,Aref,Zref);
  std::map<index_t,double> hole_list = GetOrbitsAZ(Aref,Zref);
  Init(emax,hole_list,valence);


//  std::vector<index_t> c = core;
//  std::vector<index_t> v = valence;
//  ClearVectors();
//  GetAZfromString(reference,Aref,Zref);
//  std::map<index_t,double> h = GetOrbitsAZ(Aref,Zref);
//  Init(emax, h,c,v);

}

void ModelSpace::SetLmax( int l)
{
  Lmax = l;
  int e = Emax;
  std::map<index_t,double> the_hole_list;
  std::set<index_t> the_core_list;
  std::set<index_t> the_valence_list;
  for ( auto h : holes )
  {
    Orbit& oh = GetOrbit(h);
    if ( oh.l > Lmax )
    {
      std::cout << "!!! NOT GOOD. Hole state " << h << " has l = " << oh.l
                << "  which is bigger than Lmax = " << Lmax << " ... dying now. " << std::endl;
      exit(0);
    }
    the_hole_list[h] = oh.occ;
  }
  for (auto c : core )
  {
    Orbit& oc = GetOrbit(c);
    if ( oc.l > Lmax )
    {
      std::cout << "!!! NOT GOOD. Core state " << c << " has l = " << oc.l
                << "  which is bigger than Lmax = " << Lmax << " ... dying now. " << std::endl;
      exit(0);
    }
    the_core_list.insert(c);
  }
  for (auto v : valence )
  {
    Orbit& ov = GetOrbit(v);
    if ( ov.l > Lmax )
    {
      continue;
//      std::cout << "!!! NOT GOOD. Valence state " << v << " has l = " << ov.l
//                << "  which is bigger than Lmax = " << Lmax << " ... dying now. " << std::endl;
//      exit(0);
    }
    the_valence_list.insert(v);
  }
  Init( e, the_hole_list, the_core_list, the_valence_list );
}


// Get std::vector of orbit indices from std::vector of std::strings
// e.g. "p0f7" gives the index of the proton 0f7/2 orbit.
std::vector<index_t> ModelSpace::String2Index( std::vector<std::string> vs )
{
  std::vector<index_t> vi;
  std::vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};

  for ( auto& s : vs )
  {
    int n,l,j2,tz2;
    tz2 = s[0]=='p' ? -1 : 1;
    std::istringstream( s.substr(1,2) ) >> n;
    l = find(l_list.begin(),l_list.end(), s[2]) - l_list.begin();
    std::istringstream( s.substr(3,s.size()) ) >> j2;
    vi.push_back( Index1(n,l,j2,tz2) );
  }
  return vi;
}


std::string ModelSpace::Index2String( index_t ind)
{
  std::vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};
  Orbit& oi = GetOrbit(ind);
  char c[10];
  char pn = oi.tz2 < 0 ? 'p' : 'n';
  char lstr = l_list[oi.l];
  sprintf(c, "%c%d%c%d",pn,oi.n,lstr,oi.j2);
  return std::string(c) ;
}


void ModelSpace::GetAZfromString(std::string str,int& A, int& Z) // TODO: accept different formats, e.g. 22Na vs Na22
{
  std::vector<std::string> periodic_table = {"n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
                        "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
                        "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
                        "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
                        "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
                        "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md"};
  if (str == "vacuum") str="n0";
  int i=0;
  while (! isdigit(str[i])) i++;
  std::string elem = str.substr(0,i);
  std::stringstream( str.substr(i,str.size()-i)) >> A;
  auto it_elem = find(periodic_table.begin(),periodic_table.end(),elem);
  if (it_elem != periodic_table.end())
  {
    Z = it_elem - periodic_table.begin();
  }
  else
  {
    Z =-1;
   std::cout << "ModelSpace::GetAZfromString :  Trouble parsing " << str << std::endl;
  }
}

// Fill A orbits with Z protons and A-Z neutrons
// assuming a standard shell-model level ordering
std::map<index_t,double> ModelSpace::GetOrbitsAZ(int A, int Z)
{
//  std::cout << "In GetOrbitsAZ  " << A << " " <<Z << std::endl;
//  std::cout << "Size of OrbitLookup = " << OrbitLookup.size() << std::endl;
  int zz = 0;
  int nn = 0; // unfortunate there are so many n's here...
  std::map<index_t,double> holesAZ;
  for (int N=0; N<=Emax; ++N)
  {
    for (int g=2*N+1;g>=-2*N;g-=4)
    {
      int j2 = std::abs(g);
      int l = g<0 ? (j2+1)/2 : (j2-1)/2;
      int n = (N-l)/2;

      if (zz < Z)
      {
        int dz = std::min(Z-zz,j2+1);
        holesAZ[Index1(n,l,j2,-1)] = dz/(j2+1.0);
        zz += dz;
      }
      if (nn < A-Z)
      {
        int dn = std::min(A-Z-nn,j2+1);
        holesAZ[Index1(n,l,j2,1)] = dn/(j2+1.0);
        nn += dn;
      }
      if (zz==Z and nn==A-Z)
      {
         return holesAZ; // We're all done here.
      }
    }
  }
  std::cout << "Trouble! Model space not big enough to fill A=" << A << " Z="<< Z << "  emax = " << Emax << std::endl;
  return holesAZ;

}


/// Find the valence space of one single major oscillator shell each for protons and neutrons (not necessarily
/// the same shell for both) which contains the naive shell-model ground state of the reference.
/// Here, the naive shell model ordering is harmonic oscillator shells split by a spin orbit.
/// For example, if we want to treat C20, with 6 protons and 14 neutrons, we take the 0p shell for protons
/// and 1s0d shell for neutrons.
//void ModelSpace::Get0hwSpace(int Aref, int Zref, std::vector<index_t>& core_list, std::vector<index_t>& valence_list)
void ModelSpace::Get0hwSpace(int Aref, int Zref, std::set<index_t>& core_list, std::set<index_t>& valence_list)
{
  int Nref = Aref-Zref;
  int OSC_protons=0,OSC_neutrons=0;
  while ( (OSC_protons +1)*(OSC_protons +2)*(OSC_protons +3)/3 <= Zref ) OSC_protons++;
  while ( (OSC_neutrons+1)*(OSC_neutrons+2)*(OSC_neutrons+3)/3 <= Nref ) OSC_neutrons++;

  int Zcore = (OSC_protons )*(OSC_protons +1)*(OSC_protons +2)/3;
  int Ncore = (OSC_neutrons)*(OSC_neutrons+1)*(OSC_neutrons+2)/3;

//  for (auto& it_core : GetOrbitsAZ(Zcore+Ncore,Zcore)) core_list.push_back(it_core.first);
  for (auto& it_core : GetOrbitsAZ(Zcore+Ncore,Zcore)) core_list.insert(it_core.first);

  if (Zref>Zcore) // if we have a closed major HO shell of protons, then don't decouple any valence proton orbits
  {
    for (int L=OSC_protons; L>=0; L-=2)
    {
      for (int j2=2*L+1;j2>std::max(2*L-2,0);j2-=2)
      {
//        valence_list.push_back( GetOrbitIndex( (OSC_protons-L)/2, L, j2, -1) );
        valence_list.insert( GetOrbitIndex( (OSC_protons-L)/2, L, j2, -1) );
      }
    }
  }
  if (Nref>Ncore)  // likewise for neutrons
  {
    for (int L=OSC_neutrons; L>=0; L-=2)
    {
      for (int j2=2*L+1;j2>std::max(2*L-2,0);j2-=2)
      {
//        valence_list.push_back( GetOrbitIndex( (OSC_neutrons-L)/2, L, j2, 1) );
        valence_list.insert( GetOrbitIndex( (OSC_neutrons-L)/2, L, j2, 1) );
      }
    }
  }

}


// Parse a std::string containing a comma-separated list of core + valence orbits
// eg, the usual sd shell would look like "O16,p0d5,n0d5,p0d3,n0d3,p1s1,n1s1".
// The number of ways to specify a model space is getting a bit out of hand...
//void ModelSpace::ParseCommaSeparatedValenceSpace(std::string valence, std::vector<index_t>& core_list, std::vector<index_t>& valence_list)
void ModelSpace::ParseCommaSeparatedValenceSpace(std::string valence, std::set<index_t>& core_list, std::set<index_t>& valence_list)
{
  std::istringstream ss(valence);
  std::string orbit_str,core_str;
  getline(ss, core_str, ',');

  int Ac,Zc;
  GetAZfromString(core_str,Ac,Zc);
  for (auto& it_core : GetOrbitsAZ(Ac,Zc) )
  {
    core_list.insert(it_core.first);
//    core_list.push_back(it_core.first);
  }

  while(getline(ss, orbit_str, ','))
  {
    valence_list.insert( String2Index({orbit_str})[0]);
//    valence_list.push_back( String2Index({orbit_str})[0]);
  }
}

// For backwards compatibility
void ModelSpace::SetReference(std::vector<index_t> new_reference)
{
  std::set<index_t> ref( new_reference.begin(),new_reference.end());
  SetReference( ref );
}

void ModelSpace::SetReference(std::set<index_t> new_reference)
{
//  std::cout << "I'm in this SetReference" << std::endl;
//  std::vector<index_t> c = core;
//  std::vector<index_t> v = valence;
  std::set<index_t> c = core;
  std::set<index_t> v = valence;
  std::map<index_t,double> h;
  for (auto r : new_reference) h[r] = 1.0;
  ClearVectors();
  Init(Emax, h,c,v);
}

void ModelSpace::SetReference(std::map<index_t,double> new_reference)
{
//  std::vector<index_t> c = core;
//  std::vector<index_t> v = valence;
  std::set<index_t> c = core;
  std::set<index_t> v = valence;
  if (valence.size()<1) // If we have no valece space, assume it's a single ref and core should equal the reference.
  {
    c.clear();
    for ( auto iter : new_reference )  c.insert(iter.first);
//    c.resize(0);
//    for ( auto iter : new_reference )  c.push_back(iter.first);
  }
  ClearVectors();
  Init(Emax, new_reference,c,v);
}

void ModelSpace::SetReference(std::string new_reference)
{
//  std::vector<index_t> c = core;
//  std::vector<index_t> v = valence;
  std::set<index_t> c = core;
  std::set<index_t> v = valence;
  ClearVectors();
  GetAZfromString(new_reference,Aref,Zref);
  std::map<index_t,double> h = GetOrbitsAZ(Aref,Zref);
  Init(Emax, h,c,v);
}

ModelSpace ModelSpace::operator=(const ModelSpace& ms)
{
   holes =  ms.holes;
   particles =  ms.particles;
   valence = ms.valence;
   qspace =  ms.qspace;
   core = ms.core;
   proton_orbits =  ms.proton_orbits;
   neutron_orbits =  ms.neutron_orbits;
   all_orbits  =  ms.all_orbits;
   KetIndex_pp =  ms.KetIndex_pp;
   KetIndex_ph =  ms.KetIndex_ph;
   KetIndex_hh =  ms.KetIndex_hh;
   KetIndex_cc =  ms.KetIndex_cc;
   KetIndex_vc =  ms.KetIndex_vc;
   KetIndex_qc =  ms.KetIndex_qc;
   KetIndex_vv =  ms.KetIndex_vv;
   KetIndex_qv =  ms.KetIndex_qv;
   KetIndex_qq =  ms.KetIndex_qq;
   Ket_occ_hh  =  ms.Ket_occ_hh;
   Ket_unocc_hh  =  ms.Ket_unocc_hh;
   Emax = ms.Emax;
   E2max = ms.E2max;
   E3max = ms.E3max;
   Lmax2 = ms.Lmax2;
   Lmax3 = ms.Lmax3;
   OneBodyJmax = ms.OneBodyJmax;
   TwoBodyJmax = ms.TwoBodyJmax;
   ThreeBodyJmax = ms.ThreeBodyJmax;
   OneBodyChannels = ms.OneBodyChannels;
   SortedTwoBodyChannels = ms.SortedTwoBodyChannels;
   SortedTwoBodyChannels_CC = ms.SortedTwoBodyChannels_CC;
   norbits = ms.norbits;
   hbar_omega = ms.hbar_omega;
   target_mass = ms.target_mass;
   target_Z = ms.target_Z;
   Aref = ms.Aref;
   Zref = ms.Zref;
   Orbits = ms.Orbits;
   Kets = ms.Kets;
   TwoBodyChannels = ms.TwoBodyChannels;
   TwoBodyChannels_CC = ms.TwoBodyChannels_CC;
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;

//   std::cout << "In copy assignment for ModelSpace" << std::endl;
   return ModelSpace(*this);
}



ModelSpace ModelSpace::operator=(ModelSpace&& ms)
{
   holes =  std::move(ms.holes);
   particles =  std::move(ms.particles);
   valence = std::move(ms.valence);
   qspace =  std::move(ms.qspace);
   core = std::move(ms.core);
   proton_orbits =  std::move(ms.proton_orbits);
   neutron_orbits =  std::move(ms.neutron_orbits);
   all_orbits =  std::move(ms.all_orbits);
   KetIndex_pp =  std::move(ms.KetIndex_pp);
   KetIndex_ph =  std::move(ms.KetIndex_ph);
   KetIndex_hh =  std::move(ms.KetIndex_hh);
   KetIndex_cc =  std::move(ms.KetIndex_cc);
   KetIndex_vc =  std::move(ms.KetIndex_vc);
   KetIndex_qc =  std::move(ms.KetIndex_qc);
   KetIndex_vv =  std::move(ms.KetIndex_vv);
   KetIndex_qv =  std::move(ms.KetIndex_qv);
   KetIndex_qq =  std::move(ms.KetIndex_qq);
   Ket_unocc_hh =  std::move(ms.Ket_unocc_hh);
   Ket_occ_hh =  std::move(ms.Ket_occ_hh);
   Emax = std::move(ms.Emax);
   E2max = std::move(ms.E2max);
   E3max = std::move(ms.E3max);
   Lmax2 = std::move(ms.Lmax2);
   Lmax3 = std::move(ms.Lmax3);
   OneBodyJmax = std::move(ms.OneBodyJmax);
   TwoBodyJmax = std::move(ms.TwoBodyJmax);
   ThreeBodyJmax = std::move(ms.ThreeBodyJmax);
   OneBodyChannels = std::move(ms.OneBodyChannels);
   SortedTwoBodyChannels = std::move(ms.SortedTwoBodyChannels);
   SortedTwoBodyChannels_CC = std::move(ms.SortedTwoBodyChannels_CC);
   norbits = std::move(ms.norbits);
   hbar_omega = std::move(ms.hbar_omega);
   target_mass = std::move(ms.target_mass);
   target_Z = std::move(ms.target_Z);
   Aref = std::move(ms.Aref);
   Zref = std::move(ms.Zref);
   Orbits = std::move(ms.Orbits);
   Kets = std::move(ms.Kets);
   TwoBodyChannels = std::move(ms.TwoBodyChannels);
   TwoBodyChannels_CC = std::move(ms.TwoBodyChannels_CC);
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
   for (TwoBodyChannel& tbc : ms.TwoBodyChannels)   tbc.modelspace = NULL;
   for (TwoBodyChannel_CC& tbc_cc : ms.TwoBodyChannels_CC)   tbc_cc.modelspace = NULL;
   return ModelSpace(*this);
}


Orbit& ModelSpace::GetOrbit(int i)
{
  if (i==NOT_AN_ORBIT) return NULL_ORBIT;
  return (Orbit&) Orbits[i];
} 

void ModelSpace::AddOrbit(Orbit orb)
{
  AddOrbit(orb.n, orb.l, orb.j2, orb.tz2, orb.occ, orb.cvq);
}

void ModelSpace::AddOrbit(int n, int l, int j2, int tz2, double occ, int cvq)
{
//   index_t ind = Index1(n, l, j2, tz2);
//   Orbits[ind] = Orbit(n,l,j2,tz2,occ,cvq,ind);
   size_t hash = Index1_hash(n,l,j2,tz2);
   index_t ind;
   if ( OrbitLookup.find(hash) == OrbitLookup.end() ) // we don't already have this orbit in the list.
   {
     ind = Orbits.size();
     Orbits.emplace_back( Orbit(n,l,j2,tz2,occ,cvq,ind ) );
     OrbitLookup[hash] = ind;
   }
   else  // we already have that one, but we'll replace it with the new info.
   {
     ind = OrbitLookup[hash];
//     std::cout << "WARNING!!!  Adding an orbit that we already had: " << ind << " :  " << n << " " << l << " " << j2 << " " << tz2 << "    this could cause trouble." << std::endl;
//     std::cout << "size of Orbits = " << Orbits.size() << std::endl;
     Orbits[ind] = Orbit(n,l,j2,tz2,occ,cvq,ind);
//     std::cout << "Orbits[" << ind << "]  now points to an orbit who's ind is " << Orbits[ind].index << std::endl;

     // Since we probably sorted into lists last time we added it, we need to remove the index ind from
     // the lists so we can sort it properly according to the new cvq info.
//..     for ( auto& vec : {particles,holes,core,valence,qspace,proton_orbits,neutron_orbits,all_orbits}  )

// now that we use a set for these, we don't need to worry about this.
//     for ( auto vec : {&particles,&holes,&core,&valence,&qspace,&proton_orbits,&neutron_orbits,&all_orbits}  )
//     {
//       vec->erase( std::remove( vec->begin(), vec->end(), ind ), vec->end() );
//     }


   }

//   Orbit& o = Orbits[ind];
//   std::cout << "Now Orbit " << hash << " -> " << ind << "   is " << o.n << " " << o.l << " " << o.j2 << " " << o.tz2 << "    " << o.occ << " " << o.cvq << std::endl;

   if (j2 > OneBodyJmax)
   {
      OneBodyJmax = j2;
      TwoBodyJmax = OneBodyJmax;
      ThreeBodyJmax = OneBodyJmax*3-1; // It doesn't seem like this is actually used anywhere.
      nTwoBodyChannels = 2*3*(TwoBodyJmax+1);
      if (single_species) nTwoBodyChannels = 2*(TwoBodyJmax+1);
   }

   for ( auto orbitlist : {&particles,&holes,&core,&valence,&qspace,&proton_orbits,&neutron_orbits,&all_orbits}  ) orbitlist->erase(ind); // 

   if ( occ < OCC_CUT) particles.insert(ind);
   else holes.insert(ind);
   if (cvq == 0) core.insert(ind);
   if (cvq == 1) valence.insert(ind);
   if (cvq == 2) qspace.insert(ind);
   if (tz2 < 0 ) proton_orbits.insert(ind);
   if (tz2 > 0 ) neutron_orbits.insert(ind);
   all_orbits.insert(ind);


//   if ( occ < OCC_CUT) particles.push_back(ind);
//   else holes.push_back(ind);
//   if (cvq == 0) core.push_back(ind);
//   if (cvq == 1) valence.push_back(ind);
//   if (cvq == 2) qspace.push_back(ind);
//   if (tz2 < 0 ) proton_orbits.push_back(ind);
//   if (tz2 > 0 ) neutron_orbits.push_back(ind);
//   all_orbits.push_back(ind);

//   OneBodyChannels[{l, j2, tz2}].push_back(ind); // (Evidently, we mean one-body channels for an operator with the same symmetries as the Hamiltonian).
   OneBodyChannels[{l, j2, tz2}].insert(ind); // (Evidently, we mean one-body channels for an operator with the same symmetries as the Hamiltonian).
}


void ModelSpace::FindEFermi()
{
  e_fermi = {{-1,0},{+1,0}};
  std::map<int,double> occmax = {{-1,0},{+1,0}};
  for ( auto i : holes )
  {
    Orbit& oi = GetOrbit(i);
    int ei = 2*oi.n+oi.l; // TODO: Generalize this in case we're doing atoms and want to use a different unperturbed energy
    if (ei>e_fermi[oi.tz2])
    {
       e_fermi[oi.tz2] = ei;
       occmax[oi.tz2] = oi.occ;
//       std::cout << __func__ << " updating orbit i= " << i << " ei = " << ei << "  ->  e_fermi = " << e_fermi[-1] << " " << e_fermi[1] << std::endl;
    }
  }
//  std::map<int,double> particle_e_min = {{-1,1e6},{1,1e6}};
  std::map<int,double> particle_e_min = e_fermi;
  for (auto i : particles )
  {
    Orbit& oi = GetOrbit(i);
    int ei = 2*oi.n+oi.l; // TODO: Generalize this in case we're doing atoms and want to use a different unperturbed energy
    if (ei<particle_e_min[oi.tz2])
    {
       particle_e_min[oi.tz2] = ei;
    }
  }
  // If the last level is completely filled, then we defined the fermi surface to 
  // be half way between the highest filled orbit and the lowest unfilled orbit
  if ( occmax[-1] > 0.99)  e_fermi[-1] = 0.5 * ( e_fermi[-1] + particle_e_min[-1] );
  if ( occmax[+1] > 0.99)  e_fermi[+1] = 0.5 * ( e_fermi[+1] + particle_e_min[+1] );

  std::cout << "Fermi energies:  " << e_fermi[-1] << "   " << e_fermi[+1] << std::endl;

}


size_t ModelSpace::GetOrbitIndex(std::string orb)
{
  std::vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};
  int n=-1,l=-1,j2=-1;
  int tz2 = orb[0]=='p' ? -1 : 1;
  std::stringstream(orb.substr(1,1)) >> n;
  auto it_l = find(l_list.begin(), l_list.end(), orb[2]);
  if ( it_l != l_list.end() )
    l = it_l - l_list.begin();
  else
    std::cout << "Bad orbit label " << orb << std::endl;
  std::stringstream(orb.substr(3)) >> j2;
  return Index1(n,l,j2,tz2);
}

size_t ModelSpace::GetTwoBodyChannelIndex(int j, int p, int t)
{
//   return (t+1)*2*(TwoBodyJmax+1) + p*(TwoBodyJmax+1) + j;
   if (single_species)
   {
     return 2*j + p;
   }
   return 6*j + 2*(t+1) + p;
}

void ModelSpace::UnpackTwoBodyChannelIndex( size_t ch, int& j, int& p, int& tz)
{
//     int tbjmax = TwoBodyJmax;
//     j = ch%(tbjmax+1);
//     p = (ch/(tbjmax+1))%2;
//     tz = (ch/(2*(tbjmax+1))-1);
     if (single_species)
     {
      j = ch/2;
      p = ch%2;
      tz = -1;
     }
     else
     {
      j = ch/6;
      tz = (ch%6)/2 -1;
      p = ch%2;
     }
}


size_t ModelSpace::Index1(int n, int l, int j2, int tz2) const 
{
  auto iter = OrbitLookup.find( Index1_hash(n,l,j2,tz2) ) ;
  if ( iter == OrbitLookup.end() ) return NOT_AN_ORBIT;
  else return iter->second;
//  return OrbitLookup.at( Index1_hash( n,l,j2,tz2) );
}


size_t ModelSpace::Index1_hash(int n, int l, int j2, int tz2) const
{
  if ( single_species )
  {
    return ((2*n+l)*(2*n+l+3) + 1-j2 )/2 ;
  }
  else
  {
    return (2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2 ;
  }
}


size_t ModelSpace::Index2(size_t p, size_t q) const 
{
  if ( single_species )
  {
    return p*(2*all_orbits.size()-1-p)/2 + q;
  }
  else
  {
    return p*(2*all_orbits.size()-1-p)/2 + q;
  }
}




void ModelSpace::SetupKets()
{
   Kets.resize(Index2(all_orbits.size()-1,all_orbits.size()-1)+1);
   for (auto p : all_orbits )
   {
     for (auto q : all_orbits )
     {
        if (q<p) continue;
        index_t index = Index2(p,q);
        Kets[index] = Ket(GetOrbit(p),GetOrbit(q));
     }
   }
    for (auto p : all_orbits)
    {
    for (auto q : all_orbits)
    {
     if (q<p) continue;
    index_t index = Index2(p,q);
    Ket& ket = Kets[index];
    int Tz = (ket.op->tz2 + ket.oq->tz2)/2;
    int parity = (ket.op->l + ket.oq->l)%2;
//   The old way this was written led to undefined behavior, depending on when the structure was expanded.
//    MonopoleKets[Tz+1][parity][index] = MonopoleKets[Tz+1][parity].size()-1;
    index_t size = MonopoleKets[Tz+1][parity].size();
    MonopoleKets[Tz+1][parity][index] = size;
    double occp = ket.op->occ;
    double occq = ket.oq->occ;
    int cvq_p = ket.op->cvq;
    int cvq_q = ket.oq->cvq;
    if (cvq_p+cvq_q==0)      KetIndex_cc.push_back(index); // 00
    if (cvq_p+cvq_q==1)      KetIndex_vc.push_back(index); // 01
    if (std::abs(cvq_p-cvq_q)==2) KetIndex_qc.push_back(index); // 02
    if (cvq_p*cvq_q==1)      KetIndex_vv.push_back(index); // 11
    if (cvq_p+cvq_q==3)      KetIndex_qv.push_back(index); // 12
    if (cvq_p+cvq_q==4)      KetIndex_qq.push_back(index); // 22
    if (occp<OCC_CUT and occq<OCC_CUT) KetIndex_pp.push_back(index);
    if ( (occp>OCC_CUT) xor (occq>OCC_CUT) )
    {
       KetIndex_ph.push_back(index);
       Ket_occ_ph.push_back(occp*occq);
       Ket_unocc_ph.push_back((1-occp)*(1-occq));
    }
    if (occp>OCC_CUT and occq>OCC_CUT)
    {
       KetIndex_hh.push_back(index);
       Ket_occ_hh.push_back(occp*occq);
       Ket_unocc_hh.push_back((1-occp)*(1-occq));
    }
   }
   }

   SortedTwoBodyChannels.resize(nTwoBodyChannels);
   SortedTwoBodyChannels_CC.resize(nTwoBodyChannels);
   for (int ch=0;ch<nTwoBodyChannels;++ch)
   {
      TwoBodyChannels.emplace_back(TwoBodyChannel(ch,this));
      TwoBodyChannels_CC.emplace_back(TwoBodyChannel_CC(ch,this));
      SortedTwoBodyChannels[ch] = ch;
      SortedTwoBodyChannels_CC[ch] = ch;
   }
   // Hopefully this can help with load balancing.
   sort(SortedTwoBodyChannels.begin(),SortedTwoBodyChannels.end(),[this](int i, int j){ return TwoBodyChannels[i].GetNumberKets() > TwoBodyChannels[j].GetNumberKets(); }  );
   sort(SortedTwoBodyChannels_CC.begin(),SortedTwoBodyChannels_CC.end(),[this](int i, int j){ return TwoBodyChannels_CC[i].GetNumberKets() > TwoBodyChannels_CC[j].GetNumberKets(); }  );
   while (  TwoBodyChannels[ SortedTwoBodyChannels.back() ].GetNumberKets() <1 ) SortedTwoBodyChannels.pop_back();
   while (  TwoBodyChannels_CC[ SortedTwoBodyChannels_CC.back() ].GetNumberKets() <1 ) SortedTwoBodyChannels_CC.pop_back();
}


void ModelSpace::SetE3max(int e3)
{
  E3max = e3;
  Setup3bKets();
}

/// We keep things relatively simple (?) for now.
/// Just make a vector of all the possible 3b kets
void ModelSpace::Setup3bKets()
{
//  Kets3.resize(0);
  Kets3.clear();
  Ket3IndexLookup.clear(); 
  ThreeBodyChannels.clear();
  ThreeBodyChannelLookup.clear();

  // I'm using a set here because it only stores unique
  // elements, so we don't need to worry about that
  // in the loop.
  std::set<std::array<int,3>> channels_found;

  for ( auto p : all_orbits )
  {
    Orbit& op = GetOrbit(p);
    for ( auto q : all_orbits )
    {
//      if (q>p) continue;
      if (q<p) continue; // this ordering matches the two body storage
      Orbit& oq = GetOrbit(q);
      int Jpq_min = std::abs(op.j2-oq.j2)/2;
      int Jpq_max = (op.j2+oq.j2)/2;
      for ( auto r : all_orbits )
      {
//        if (r>q) continue;
        if (r<q) continue;
        Orbit& oR = GetOrbit(r);
        if ( ( 2*(op.n+oq.n+oR.n)+op.l+oq.l+oR.l) > E3max ) continue;
//        int parity = (op.l+oq.l+oR.l)%2;
//        int twoTz = op.tz2+oq.tz2+oR.tz2;
        for (int Jpq=Jpq_min; Jpq<=Jpq_max; Jpq++)
        {
          if (p==q and Jpq%2>0) continue;
          if (p==q and p==r and  std::abs(op.tz2+oq.tz2+oR.tz2)==3 and op.j2==1) continue;
          Kets3.push_back( Ket3(op,oq,oR,Jpq) );
          Ket3IndexLookup[ Ket3IndexHash(p,q,r,Jpq)] = Kets3.size()-1; // for reverse lookup
//          int twoJ_min = std::abs( 2*Jpq - oR.j2 );
//          int twoJ_max = 2*Jpq + oR.j2;
//          for (int twoJ=twoJ_min; twoJ_max; twoJ+=2) channels_found.insert({twoJ,parity,twoTz});
        }
        
      }
    }
  }
//  std::cout << "Done with loop over orbits. size of Kets3 = " << Kets3.size() << std::endl;

  int twoJ_min=1;
  int twoJ_max= 6*Emax + 1; // 2 * (3*(Emax+1/2)-1)
  for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
  {
   for (int parity=0; parity<=1; parity++)
   {
    for (int twoTz=-3; twoTz<=3; twoTz+=2)
    {
      channels_found.insert({twoJ,parity,twoTz});
    }
   }
  }
//  std::cout << "done with loop over channels found" << std::endl;
  // Now we store all the 3-body channels
  for ( auto JPT : channels_found )
  {
//    std::cout << "JPT " << JPT[0] << " " << JPT[1] << " " << JPT[2] << std::endl;
    ThreeBodyChannels.push_back( ThreeBodyChannel( JPT[0],JPT[1],JPT[2], this )  );
//    std::cout << "pushed it back." << std::endl;
    ThreeBodyChannelLookup[ ThreeBodyChannelHash( JPT[0],JPT[1],JPT[2]) ] = ThreeBodyChannels.size()-1;
//    std::cout << "3Bchan: " << ThreeBodyChannels.size() << " "  << JPT[0] << " " << JPT[1] << " " << JPT[2] << "   hash " << ThreeBodyChannelHash( JPT[0],JPT[1],JPT[2]) << std::endl;
  }
  nThreeBodyChannels = ThreeBodyChannels.size();

  

}


// Turn p,q,r,Jpq into a single index
// For now at least, we are storing p>=q>=r
size_t ModelSpace::Ket3IndexHash(size_t p, size_t q, size_t r, size_t Jpq)
{
  size_t hash = ( p << 24 )
              + ( q << 16 ) 
              + ( r << 8 ) 
              + Jpq;
  return hash;
}



size_t ModelSpace::GetThreeBodyChannelIndex(int twoJ, int parity, int twoTz )
{
  auto iter = ThreeBodyChannelLookup.find( ThreeBodyChannelHash( twoJ, parity, twoTz));
  if ( iter == ThreeBodyChannelLookup.end() ) return -1;
  return iter->second;
//  std::cout << "IN " << __func__ << std::endl;
//  std::cout << "  with Jpt = " << twoJ << " " << parity << " " << twoTz << std::endl;
//  std::cout << " the hash is " << ThreeBodyChannelHash( twoJ , parity, twoTz) << std::endl;
//  std::cout << "I would have thought it was " << (  4*(twoJ-1) + twoTz+3 + parity) << std::endl;
//  std::cout << " An that should point to " <<ThreeBodyChannelLookup.at( ThreeBodyChannelHash( twoJ, parity, twoTz)) << std::endl;
//   return ThreeBodyChannelLookup.at( ThreeBodyChannelHash( twoJ, parity, twoTz)) ;
}

size_t ModelSpace::ThreeBodyChannelHash( int twoJ, int parity, int twoTz)
{
  size_t hash = (4*(twoJ-1) + (twoTz+3) + parity);
  return hash;
}




void ModelSpace::SetEmaxUnocc(int e)
{
  std::cout << __func__ << " " << e << std::endl;
  EmaxUnocc = std::min(Emax,e);
  if (e > Emax)
  {
    std::cout << "WARNING: " << __func__ << "  tried setting EmaxUnocc to " << e
              << "  which is > Emax = " << Emax << ". Setting EmaxUnocc to Emax" << std::endl;
  }

//  std::set<std::array<int,3>> hole_quantum_numbers; // For checking if an orbit could mix with the hole orbits

  hole_quantum_numbers.clear();
  std::map<index_t,double> holemap;
  int max_l = -1;
  for (auto h : holes )
  {
    Orbit& oh = GetOrbit(h);
    holemap[h] = oh.occ;
//    hole_quantum_numbers.insert( { oh.l, oh.j2, oh.tz2} );
//    hole_quantum_numbers.insert( { oh.l, oh.j2} );
    max_l = std::max(max_l, oh.l);
  }

  for (int l=0; l<=std::min(Emax,max_l+2); l++)
  {
    hole_quantum_numbers.insert( {l,2*l+1});
    if (l>0) hole_quantum_numbers.insert( {l,2*l-1} );
  }

//  std::cout << "hole quantum numbers: " << std::endl;
//  for (auto& hqn : hole_quantum_numbers )
//  {
//    std::cout << hqn[0] << " " << hqn[1] << " " << hqn[2] << std::endl;
//    std::cout << hqn[0] << " " << hqn[1]  << std::endl;
//  }

//  for ( auto h : holes ) holemap[h] = GetOrbit(h).occ;
//  Init(Emax, holemap,c,v);
  Init(Emax, holemap,core,valence);

}



void ModelSpace::ClearVectors()
{
   holes.clear();         
   particles.clear();     
   core.clear();          
   valence.clear();       
   qspace.clear();        
   proton_orbits.clear();  
   neutron_orbits.clear();
   all_orbits.clear();
   
   KetIndex_pp.clear();
   KetIndex_ph.clear();
   KetIndex_hh.clear();
   KetIndex_cc.clear();
   KetIndex_vc.clear();
   KetIndex_qc.clear();
   KetIndex_vv.clear();
   KetIndex_qv.clear();
   KetIndex_qq.clear();
   Ket_occ_hh.clear();
   Ket_occ_ph.clear();
   Ket_unocc_hh.clear();
   Ket_unocc_ph.clear();
   for (index_t Tz=0; Tz<3; ++Tz)
   {
     for (index_t parity=0;parity<2; ++parity)
     {
        MonopoleKets[Tz][parity].clear();
     }
   }

   Orbits.clear();
   Kets.clear();
   OneBodyChannels.clear();
   TwoBodyChannels.clear();
   TwoBodyChannels_CC.clear();
   SortedTwoBodyChannels.clear();
   SortedTwoBodyChannels_CC.clear();
   PandyaLookup.clear();
}


// On the first pass through a commutator calculation, we use a single thread so that
// we only calculate the sixJ's we need, but we don't have to worry about race conditions.
// This resets the first pass flag.
void ModelSpace::ResetFirstPass()
{
  scalar_transform_first_pass = true;
  scalar3b_transform_first_pass = true;
  for (size_t i=0;i<tensor_transform_first_pass.size();i++) tensor_transform_first_pass[i] = true;
}

uint64_t ModelSpace::SixJHash(double j1, double j2, double j3, double J1, double J2, double J3)
{
   return   (((uint64_t)(2*j1)) << 50)
          + (((uint64_t)(2*j2)) << 40)
          + (((uint64_t)(2*j3)) << 30)
          + (((uint64_t)(2*J1)) << 20)
          + (((uint64_t)(2*J2)) << 10)
          +  ((uint64_t)(2*J3));

}


void ModelSpace::SixJUnHash(uint64_t key, uint64_t& j1, uint64_t& j2, uint64_t& j3, uint64_t& J1, uint64_t& J2, uint64_t& J3)
{
   j1 = (key >> 50) & 0x3FFL;
   j2 = (key >> 40) & 0x3FFL;
   j3 = (key >> 30) & 0x3FFL;
   J1 = (key >> 20) & 0x3FFL;
   J2 = (key >> 10) & 0x3FFL;
   J3 = (key      ) & 0x3FFL;
}

uint64_t ModelSpace::MoshinskyHash(uint64_t N, uint64_t Lam, uint64_t n, uint64_t lam, uint64_t n1, uint64_t l1, uint64_t n2, uint64_t l2, uint64_t L)
{
   return   (N   << 54)
          + (Lam << 47)
          + (n   << 41)
          + (lam << 34)
          + (n1  << 28)
          + (l1  << 21)
          + (n2  << 15)
          + (l2  << 8 )
          +  L;
}

void ModelSpace::MoshinskyUnHash(uint64_t key,uint64_t& N,uint64_t& Lam,uint64_t& n,uint64_t& lam,uint64_t& n1,uint64_t& l1,uint64_t& n2,uint64_t& l2,uint64_t& L)
{
   N   = (key >> 54) & 0x3FL;
   Lam = (key >> 47) & 0x7FL;
   n   = (key >> 41) & 0x3FL;
   lam = (key >> 34) & 0x7FL;
   n1  = (key >> 28) & 0x3FL;
   l1  = (key >> 21) & 0x7FL;
   n2  = (key >> 15) & 0x3FL;
   l2  = (key >> 8 ) & 0x7FL;
   L   = (key      ) & 0xFFL;
}

double ModelSpace::GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3)
{
// { j1 j2 j3 }
// { J1 J2 J3 }
   uint64_t key = SixJHash(j1,j2,j3,J1,J2,J3);

   const auto it = SixJList.find(key);
   double sixj=0.0;
   if (it != SixJList.end() )
   {
     sixj = it->second;
   }
   else
   {
    sixj = AngMom::SixJ(j1,j2,j3,J1,J2,J3);
    if (omp_get_num_threads()<2)
    {
      #pragma omp critical
      {
        SixJList[key] = sixj;
      }
    }
    else
    {
//      printf("DANGER!!!!!!!  Updating SixJList inside a parellel loop breaks thread safety!\n");
//      printf(" I shouldn't be here in GetSixJ(%.1f %.1f %.1f %.1f %.1f %.1f).  key =%" PRIx64 "   sixj=%f\n",j1,j2,j3,J1,J2,J3,key,sixj); //PRIx64 is portable uint64_t format
      std::cout << "DANGER!!!!!!!  Updating SixJList inside a parellel loop breaks thread safety!" << std::endl;
      std::cout << "  I shouldn't be here in GetSixJ("
                << std::setprecision(1) << std::fixed << j1 << " " << std::setprecision(1) << std::fixed << j2 << " "
                << std::setprecision(1) << std::fixed << j3 << " " << std::setprecision(1) << std::fixed << J1 << " "
                << std::setprecision(1) << std::fixed << J2 << " " << std::setprecision(1) << std::fixed << J3 << "). key = "
                << std::hex << key << "   sixj = " << std::dec << sixj << std::endl;
      profiler.counter["N_CalcSixJ_in_Parallel_loop"] +=1;
//      quick_exit(EXIT_FAILURE);
      exit(EXIT_FAILURE);
    }
   }
   return sixj;
}



/// Loop over all the 6j symbols that we expect to encounter, and 
/// store them in a hash table.
/// Calculate all symbols
/// \f[ \begin{Bmatrix}  ja & jb & J1
///                      jc & jd & J2 \end{Bmatrix}
/// \f]
/// and 
/// \f[ \begin{Bmatrix}  J1 & J2 & J3
///                      ja & jb & jc \end{Bmatrix}
/// \f]
/// where ja,jb,jc are half-integer and J1,J2,J3 are integer.
/// ja,jb,jc run from 1/2 to emax+1/2, while jd runs higher
/// since the 3N recoupling requires it to go up to e(emax+1/2).
/// I haven't yet bothered using the symmetry properties of the
/// 6j symbol.
///
void ModelSpace::PreCalculateSixJ()
{
  if (sixj_has_been_precalculated) return;
  std::cout << "Precalculating SixJ's" << std::endl;
  double t_start = omp_get_wtime();
  std::vector<uint64_t> KEYS;
  for (int j2a=1; j2a<=(2*Emax+1); j2a+=2)
  {
   for (int j2b=1; j2b<=(2*Emax+1); j2b+=2)
   {
    for (int j2c=1; j2c<=(2*Emax+1); j2c+=2)
    {
     // four half-integer j's,  two integer J's
     for (int j2d=1; j2d<=3*(2*Emax+1); j2d+=2)
     {
      // J1 couples a,b, and c,d;  J2 couples a,d and b,c
      int J1_min = std::max( std::abs(j2a-j2b), std::abs(j2c-j2d) );
      int J1_max = std::min( j2a+j2b, j2c+j2d );
      int J2_min = std::max( std::abs(j2a-j2d), std::abs(j2b-j2c) );
      int J2_max = std::min( j2a+j2d, j2b+j2c );
      for (int J1=J1_min; J1<=J1_max; J1+=2)
      {
       for (int J2=J2_min; J2<=J2_max; J2+=2)
       {
         uint64_t key = SixJHash(0.5*j2a,0.5*j2b,0.5*J1,0.5*j2c,0.5*j2d,0.5*J2);
         if ( SixJList.count(key) == 0 ) 
         {
           KEYS.push_back(key);
           SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
         }
       } // for J2
      } // for J1
     } // for j2d

     // three half-integer j's, three integer J's
     // J1 couples a,b, and c,d;  J2 couples a,d and b,c
     int J1_min = std::abs(j2a-j2c) ;
     int J1_max = j2a+j2c;
     int J2_min = std::abs(j2b-j2c) ;
     int J2_max = j2b+j2c;
     for (int J1=J1_min; J1<=J1_max; J1+=2)
     {
      for (int J2=J2_min; J2<=J2_max; J2+=2)
      {
       int J3_min = std::max( std::abs(J1-J2), std::abs(j2a-j2b) );
       int J3_max = std::min( J1+J2, j2a+j2b );
       for (int J3=J3_min; J3<=J3_max; J3+=2)
       {
         uint64_t key = SixJHash(0.5*J1,0.5*J2,0.5*J3,0.5*j2a,0.5*j2b,0.5*j2c);
         if ( SixJList.count(key) == 0 ) 
         {
           KEYS.push_back(key);
           SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
         }
       }// for J3
      }// for J2
     }// for J1
    }// for j2c
   }// for j2b
  }// for j2a

  #pragma omp parallel for schedule(dynamic,1)
  for (size_t i=0;i< KEYS.size(); ++i)
  {
    uint64_t j1,j2,j3,J1,J2,J3;
    uint64_t key = KEYS[i];
    SixJUnHash(key, j1,j2,j3,J1,J2,J3);
    SixJList[key] = AngMom::SixJ(0.5*j1,0.5*j2,0.5*j3,0.5*J1,0.5*J2,0.5*J3);
  }
  sixj_has_been_precalculated = true;
  std::cout << "done calculating sixJs (" << KEYS.size() << " of them)" << std::endl;
  std::cout << "Hash table has " << SixJList.bucket_count() << " buckets and a load factor " << SixJList.load_factor() 
       << "  estimated storage ~ " << ((SixJList.bucket_count()+SixJList.size()) * (sizeof(size_t)+sizeof(void*))) / (1024.*1024.*1024.) << " GB" << std::endl;
  profiler.timer["PreCalculateSixJ"] += omp_get_wtime() - t_start;
}




void ModelSpace::PreCalculateMoshinsky()
{
  if (moshinsky_has_been_precalculated) return;
  double t_start = omp_get_wtime();
  std::cout <<"Calculating moshinsky with Lmax = " << Lmax << std::endl;

  // generating all the keys is fast, so we do this first without parallelization
//  std::vector<unsigned long long int> KEYS;
  std::vector<uint64_t> KEYS;
  for (int N=0; N<=E2max/2; ++N)
  {
   for (int n=0; n<=std::min(N,E2max/2-N); ++n)
   {
//    int Lam_max = std::min( E2max-2*n-2*n, 2*Lmax ); // Here Lmax is the max L of the s.p. basis
    for (int Lam=0; Lam<=E2max-2*N-2*n; ++Lam)
    {
//     int lam_max = std::min( 2*Lmax, (N==n ? std::min(Lam,E2max-2*N-2*n-Lam) : E2max-2*N-2*n-Lam)) ; // Here Lmax is the max L of the s.p. basis
     int lam_max = N==n ? std::min(Lam,E2max-2*N-2*n-Lam) : E2max-2*N-2*n-Lam ; 
     for (int lam=0; lam<=lam_max; ++lam)
     {
      int e2 = 2*N+Lam + 2*n+lam;
      for (int L=std::abs(Lam-lam); L<=Lam+lam; ++L)
      {
       if (L>2*Lmax) continue;
       for (int n1=0; n1<=N; ++n1)
       {
        for (int n2=0; n2<=std::min(n1,e2/2-n1); ++n2)
        {
         int l1max = n1==N? std::min(Lam,e2-2*n1-2*n2) : e2-2*n1-2*n2;
         for (int l1=0; l1<=l1max; ++l1 )
         {
          int l2 = e2-2*n1-2*n2-l1;
          if ( (l1+l2+lam+Lam)%2 >0 ) continue;
          if ( l2<std::abs(L-l1) or l2>L+l1 ) continue;
//          if (l1>Lmax or l2>Lmax) continue;
          if ( (l1>Lmax and l2>Lmax) and (lam>Lmax or Lam>Lmax)) continue;
          if ( (l1>Lmax or l2>Lmax) and (lam>Lmax and Lam>Lmax)) continue;

          uint64_t key = MoshinskyHash(N,Lam,n,lam,n1,l1,n2,l2,L);
          KEYS.push_back(key);
          MoshList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
         }
        }
       }
      }
     }
    }
   }
  }
  // Now we calculate the Moshinsky brackets in parallel
//  std::vector<double> mosh_vals( KEYS.size() );
  #pragma omp parallel for schedule(dynamic,1)
  for (size_t i=0;i< KEYS.size(); ++i)
  {
    uint64_t key = KEYS[i];
    uint64_t N,Lam,n,lam,n1,l1,n2,l2,L;
    MoshinskyUnHash(key,N,Lam,n,lam,n1,l1,n2,l2,L);

    MoshList[key] = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
  }

  moshinsky_has_been_precalculated = true;
  std::cout << "done calculating moshinsky (" << KEYS.size() << " elements)" << std::endl;
  std::cout << "Hash table has " << MoshList.bucket_count() << " buckets and a load factor " << MoshList.load_factor() 
            << "  estimated storage ~ " << ((MoshList.bucket_count()+MoshList.size()) * (sizeof(size_t)+sizeof(void*))) / (1024.*1024.*1024.) << " GB" << std::endl;
  profiler.timer["PreCalculateMoshinsky"] += omp_get_wtime() - t_start;
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
      std::swap(n1,n2);
      std::swap(l1,l2);
      phase_mosh *= phase(Lam+L);
      ++switches;
   }
   if (n>N or (n==N and lam>Lam))
   {
      std::swap(n,N);
      std::swap(lam,Lam);
      phase_mosh *= phase(l1 +L);
      ++switches;
   }

//   if (l1>Lam or (l1==Lam and n1>N) or (l1==Lam and n1==N and l2>lam) or (l1==Lam and n1==N and l2==lam and n2>n) )
   if (n1>N or (n1==N and l1>Lam) or (n1==N and l1==Lam and n2>n) or (n1==N and l1==Lam and n2==n and l2>lam) )
   {
      std::swap(n1,N);
      std::swap(l1,Lam);
      std::swap(n2,n);
      std::swap(l2,lam);
      ++switches;
//      phase_mosh *= phase(l2+lam); // This phase is given in Moshinsky and Brody, but with the current algorithm, it appears not to be required.
   }
  }

   uint64_t key = MoshinskyHash(N,Lam,n,lam,n1,l1,n2,l2,L);


   auto it = MoshList.find(key);
   if ( it != MoshList.end() )  return it->second * phase_mosh;
   if (omp_get_num_threads()>1)
   {
     std::cout << "TROUBLE IN MOSHINSKY LAND!!!!!    <" << N << " " << Lam << " " << n << " " << lam << " | " << n1 << " " << l1 << " " << n2 << " " << l2 << ">_" << L  << std::endl;
   }

   // if we didn't find it, we need to calculate it.
   double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
//   #pragma omp atomic
   MoshList[key] = mosh;
   return mosh * phase_mosh;

}




double ModelSpace::GetNineJ(double j1, double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
{
//   std::cout << "Calling GetNineJ" << std::endl;
   int k1 = 2*j1;
   int k2 = 2*j2;
   int K12 = 2*J12;
   int k3 = 2*j3;
   int k4 = 2*j4;
   int K34 = 2*J34;
   int K13 = 2*J13;
   int K24 = 2*J24;
   int K = 2*J;

   std::array<int,9> klist = {k1,k2,K12,k3,k4,K34,K13,K24,K};
   std::array<double,9> jlist = {j1,j2,J12,j3,j4,J34,J13,J24,J};
   int imin = std::min_element(klist.begin(),klist.end()) - klist.begin();
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
//   unsigned long long int factor = 100;
   unsigned long long int factor = 91;
   for (int i=1; i<9; ++i)
   {
      key += klist[i]*factor;
      factor *=91;
//      factor *=100;
   }
   auto it = NineJList.find(key);
   if (it != NineJList.end() )
   {
     return it->second;
   }
   double ninej = AngMom::NineJ(jlist[0],jlist[1],jlist[2],jlist[3],jlist[4],jlist[5],jlist[6],jlist[7],jlist[8]);
   #pragma omp critical
   NineJList[key] = ninej;
   return ninej;

}


//std::map<std::array<int,2>,std::vector<std::array<int,2>>>& ModelSpace::GetPandyaLookup(int rank_J, int rank_T, int parity)
//std::map<std::array<int,2>,std::array<std::vector<int>,2>>& ModelSpace::GetPandyaLookup(int rank_J, int rank_T, int parity)
std::map<std::array<size_t,2>,std::array<std::vector<size_t>,2>>& ModelSpace::GetPandyaLookup(int rank_J, int rank_T, int parity)
{
   CalculatePandyaLookup(rank_J,rank_T,parity);
   return PandyaLookup[{rank_J,rank_T,parity}];

}


// Generate a lookup table of all the channels that depend on a given set of Pandya-transformed channels
// this is used in the 222ph commutators to avoid calculating things that won't be used.
void ModelSpace::CalculatePandyaLookup(int rank_J, int rank_T, int parity)
{
   if (PandyaLookup.find({rank_J, rank_T, parity})!=PandyaLookup.end()) return; 
//   std::cout << "CalculatePandyaLookup( " << rank_J << ", " << rank_T << ", " << parity << ") " << std::endl;
   double t_start = omp_get_wtime();
//   PandyaLookup[{rank_J,rank_T,parity}] = std::map<std::array<int,2>,std::vector<std::array<int,2>>>();
//   PandyaLookup[{rank_J,rank_T,parity}] = std::map<std::array<int,2>,std::array<std::vector<int>,2>>();
   PandyaLookup[{rank_J,rank_T,parity}] = std::map<std::array<size_t,2>,std::array<std::vector<size_t>,2>>();
   auto& lookup = PandyaLookup[{rank_J,rank_T,parity}];

   size_t ntbc    = TwoBodyChannels.size();
   size_t ntbc_cc = TwoBodyChannels_CC.size();
   for (size_t ch_bra_cc = 0; ch_bra_cc<ntbc_cc; ++ch_bra_cc)
   {
     for (size_t ch_ket_cc = ch_bra_cc; ch_ket_cc<ntbc_cc; ++ch_ket_cc)
     {
//       lookup[{ch_bra_cc,ch_ket_cc}] = std::vector<std::array<int,2>>();
//       lookup[{ch_bra_cc,ch_ket_cc}] = std::array<std::vector<int>,2>(); 
       lookup[{ch_bra_cc,ch_ket_cc}] = std::array<std::vector<size_t>,2>(); 
//       lookup[{ch_bra_cc,ch_ket_cc}] = { <std::vector<int>(), std::vector<int>() }; 
       lookup.at({ch_bra_cc,ch_ket_cc})[0].reserve(ntbc_cc)  ; 
       lookup.at({ch_bra_cc,ch_ket_cc})[1].reserve(ntbc_cc)  ; 
     }
   }

   #pragma omp parallel for schedule(dynamic,1)
   for (size_t ch_bra_cc = 0; ch_bra_cc<ntbc_cc; ++ch_bra_cc)
   {
     TwoBodyChannel_CC& tbc_bra_cc = TwoBodyChannels_CC[ch_bra_cc];
     int twoJ_bra_cc = 2*tbc_bra_cc.J;
     for (size_t ch_ket_cc = ch_bra_cc; ch_ket_cc<ntbc_cc; ++ch_ket_cc)
     {
       TwoBodyChannel_CC& tbc_ket_cc = TwoBodyChannels_CC[ch_ket_cc];
//       lookup[{ch_bra_cc,ch_ket_cc}] = std::vector<std::array<int,2>>();
//       std::vector<int>& bra_list = lookup.at({ch_bra_cc,ch_ket_cc})[0];
//       std::vector<int>& ket_list = lookup.at({ch_bra_cc,ch_ket_cc})[1];
       std::vector<int> bra_list,ket_list;
       int twoJ_ket_cc = 2*tbc_ket_cc.J;
       for (size_t ch_bra=0; ch_bra<ntbc; ++ch_bra)
       {
         TwoBodyChannel& tbc_bra = TwoBodyChannels[ch_bra];
         for (size_t ch_ket=ch_bra; ch_ket<ntbc; ++ch_ket)
         {
           TwoBodyChannel& tbc_ket = TwoBodyChannels[ch_ket];
           if ( std::abs(tbc_bra.J-tbc_ket.J)>rank_J ) continue;
           if ( (tbc_bra.J+tbc_ket.J)<rank_J ) continue;
           if ( std::abs(tbc_bra.Tz-tbc_ket.Tz)>rank_T ) continue;
           if ( (tbc_bra.parity + tbc_ket.parity + parity)%2>0 ) continue;

           bool need_it = false;
           for (size_t ibra=0; ibra<tbc_bra.GetNumberKets(); ++ibra)
           {
             if (need_it) break;
             const Ket& bra = tbc_bra.GetKet(ibra);
             Orbit& oi = *(bra.op);
             Orbit& oj = *(bra.oq);
             for (size_t iket=0; iket<tbc_ket.GetNumberKets(); ++iket)
             {
               const Ket& ket = tbc_ket.GetKet(iket);
               Orbit& ok = *(ket.op);
               Orbit& ol = *(ket.oq);
               bool parity_good = ( (oi.l+ol.l)%2==tbc_bra_cc.parity      and (ok.l+oj.l)%2==tbc_ket_cc.parity );
               bool Tz_good = ( std::abs(oi.tz2+ol.tz2)==2*tbc_bra_cc.Tz   and std::abs(ok.tz2+oj.tz2)==2*tbc_ket_cc.Tz );
               if (single_species) Tz_good = (  -std::abs(oi.tz2+ol.tz2)==2*tbc_bra_cc.Tz   and -std::abs(ok.tz2+oj.tz2)==2*tbc_ket_cc.Tz );
               int j3min = std::abs(oi.j2-ol.j2);
               int j3max = oi.j2+ol.j2;
               int j4min = std::abs(ok.j2-oj.j2);
               int j4max = ok.j2+oj.j2;
               if (  parity_good   and Tz_good
//               if (   (oi.l+ol.l)%2==tbc_bra_cc.parity         and (ok.l+oj.l)%2==tbc_ket_cc.parity
//                         and std::abs(oi.tz2+ol.tz2)==2*tbc_bra_cc.Tz   and std::abs(ok.tz2+oj.tz2)==2*tbc_ket_cc.Tz
                         and j3min<=twoJ_bra_cc and twoJ_bra_cc<=j3max           and j4min<=twoJ_ket_cc and twoJ_ket_cc<=j4max )
               {
                 need_it=true;
                 break;
               }
               if (  parity_good  and Tz_good
//               if (   (oi.l+ol.l)%2==tbc_ket_cc.parity         and (ok.l+oj.l)%2==tbc_bra_cc.parity
//                         and std::abs(oi.tz2+ol.tz2)==2*tbc_ket_cc.Tz   and std::abs(ok.tz2+oj.tz2)==2*tbc_bra_cc.Tz
                         and j3min<=twoJ_ket_cc and twoJ_ket_cc<=j3max           and j4min<=twoJ_bra_cc and twoJ_bra_cc<=j4max )
               {
                 need_it=true;
                 break;
               }


               parity_good = ( (oj.l+ol.l)%2==tbc_bra_cc.parity         and (ok.l+oi.l)%2==tbc_ket_cc.parity );
               Tz_good = ( std::abs(oj.tz2+ol.tz2)==2*tbc_bra_cc.Tz   and std::abs(ok.tz2+oi.tz2)==2*tbc_ket_cc.Tz );
               if (single_species) Tz_good = (  -std::abs(oj.tz2+ol.tz2)==2*tbc_bra_cc.Tz   and -std::abs(ok.tz2+oi.tz2)==2*tbc_ket_cc.Tz );
               j3min = std::abs(oj.j2-ol.j2);
               j3max = oj.j2+ol.j2;
               j4min = std::abs(ok.j2-oi.j2);
               j4max = ok.j2+oi.j2;
               if (  parity_good and Tz_good
//               if (   (oj.l+ol.l)%2==tbc_bra_cc.parity         and (ok.l+oi.l)%2==tbc_ket_cc.parity
//                         and std::abs(oj.tz2+ol.tz2)==2*tbc_bra_cc.Tz   and std::abs(ok.tz2+oi.tz2)==2*tbc_ket_cc.Tz
                         and j3min<=twoJ_bra_cc and twoJ_bra_cc<=j3max           and j4min<=twoJ_ket_cc and twoJ_ket_cc<=j4max )
               {
                 need_it=true;
                 break;
               }
               if (  parity_good  and Tz_good
//               if (   (oj.l+ol.l)%2==tbc_ket_cc.parity         and (ok.l+oi.l)%2==tbc_bra_cc.parity
//                         and std::abs(oj.tz2+ol.tz2)==2*tbc_ket_cc.Tz   and std::abs(ok.tz2+oi.tz2)==2*tbc_bra_cc.Tz
                         and j3min<=twoJ_ket_cc and twoJ_ket_cc<=j3max           and j4min<=twoJ_bra_cc and twoJ_bra_cc<=j4max )
               {
                 need_it=true;
                 break;
               }


             }
           }
           if (need_it)
           {
//             lookup[{ch_bra_cc,ch_ket_cc}].push_back({ch_bra,ch_ket});
             bra_list.push_back(ch_bra);
             ket_list.push_back(ch_ket);
           }
         }
       }
       lookup.at({ch_bra_cc,ch_ket_cc})[0].assign( begin(bra_list),end(bra_list) );
       lookup.at({ch_bra_cc,ch_ket_cc})[1].assign( begin(ket_list),end(ket_list) );
     }
   }
   profiler.timer["CalculatePandyaLookup"] += omp_get_wtime() - t_start;
//   std::cout << "done." << std::endl;
}



