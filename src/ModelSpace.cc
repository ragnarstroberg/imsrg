#include "ModelSpace.hh"
#include "AngMom.hh"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

Orbit::Orbit()
{
}

Orbit::Orbit(int nn, int ll, int jj2, int ttz2, int pph, int iio, double e=0.0)
{
   n = nn;
   l = ll;
   j2 = jj2;
   tz2 = ttz2;
   ph = pph;
   io = iio;
   spe = e;
}

Orbit::Orbit(const Orbit& orb)
{
  n = orb.n;
  l = orb.l;
  j2 = orb.j2;
  tz2 = orb.tz2;
  ph = orb.ph;
  io = orb.io;
  spe = orb.spe;
}

void Orbit::Set(int nn, int ll, int jj2, int ttz2, int pph, int iio, double e)
{
   n = nn;
   l = ll;
   j2 = jj2;
   tz2 = ttz2;
   ph = pph;
   io = iio;
   spe = e;
}


//************************************************************************

Ket::Ket(ModelSpace * modelspace)
{
   ms = modelspace;
}

Ket::Ket(ModelSpace * modelspace, int pp, int qq)
{
   ms = modelspace;
   Setpq(pp,qq);
}

void Ket::Setpq(int pp, int qq)
{
   p = pp;
   q = qq;
   Orbit & op = ms->GetOrbit(p);
   Orbit & oq = ms->GetOrbit(q);
   parity = (op.l + oq.l)%2;
   Tz = (op.tz2 + oq.tz2)/2;
   Jmin = abs(op.j2 - oq.j2)/2;
   Jmax = (op.j2 + oq.j2)/2;
   Jstep = 1;
   E2 = 2*(op.n+oq.n)+op.l+oq.l;
   if (p==q) // Pauli principle
   { 
      Jmax = Jmax-1;
      Jstep=2;
   }
   phase_prefactor = ms->phase((op.j2+oq.j2)/2 + 1);
   dpq = p==q ? 1 : 0;
}

int Ket::Phase(int J)
{
   return phase_prefactor * ms->phase(J);
}

//************************************************************************

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
   J = N%tbjmax;
   parity = (N/tbjmax)%2;
   Tz = (N/(2*tbjmax)-1);
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
}


void TwoBodyChannel::Copy( const TwoBodyChannel& rhs)
{
   J                 = rhs.J;
   parity            = rhs.parity;
   Tz                = rhs.Tz;
   modelspace        = rhs.modelspace;
   NumberKets        = rhs.NumberKets;
   Proj_hh           = rhs.Proj_hh;
   Proj_pp           = rhs.Proj_pp;
   KetMap            = rhs.KetMap;
   KetList           = rhs.KetList;

}


int TwoBodyChannel::GetLocalIndex(int p, int q) const { return KetMap[modelspace->GetKetIndex(p,q)];}; 

// get pointer to ket using local index
Ket & TwoBodyChannel::GetKet(int i) const { return modelspace->GetKet(KetList[i]);}; 


bool TwoBodyChannel::CheckChannel_ket(int p, int q) const
{
   if ((p==q) and (J%2 != 0)) return false; // Pauli principle
   Orbit & op = modelspace->GetOrbit(p);
   Orbit & oq = modelspace->GetOrbit(q);
   if ((op.l + oq.l)%2 != parity) return false;
   if ((op.tz2 + oq.tz2) != 2*Tz) return false;
   if (op.j2 + oq.j2 < 2*J)       return false;
   if (abs(op.j2 - oq.j2) > 2*J)  return false;

   return true;
}

arma::uvec TwoBodyChannel::GetKetIndex_pp() { return GetKetIndexFromList(modelspace->KetIndex_pp);};
arma::uvec TwoBodyChannel::GetKetIndex_hh() { return GetKetIndexFromList(modelspace->KetIndex_hh);};
arma::uvec TwoBodyChannel::GetKetIndex_ph() { return GetKetIndexFromList(modelspace->KetIndex_ph);};
arma::uvec TwoBodyChannel::GetKetIndex_vv() { return GetKetIndexFromList(modelspace->KetIndex_vv);};
arma::uvec TwoBodyChannel::GetKetIndex_holeq_holeq() { return GetKetIndexFromList(modelspace->KetIndex_holeq_holeq);};
arma::uvec TwoBodyChannel::GetKetIndex_particleq_particleq() { return GetKetIndexFromList(modelspace->KetIndex_particleq_particleq);};
arma::uvec TwoBodyChannel::GetKetIndex_particleq_holeq() { return GetKetIndexFromList(modelspace->KetIndex_particleq_holeq);};
arma::uvec TwoBodyChannel::GetKetIndex_v_holeq() { return GetKetIndexFromList(modelspace->KetIndex_v_holeq);};
arma::uvec TwoBodyChannel::GetKetIndex_v_particleq(){ return GetKetIndexFromList(modelspace->KetIndex_v_particleq);};


//arma::uvec TwoBodyChannel::GetKetIndexFromList(vector<unsigned int>& vec_in)
arma::uvec TwoBodyChannel::GetKetIndexFromList(vector<long long unsigned int>& vec_in)
{
//   vector<unsigned int> index_list (min(vec_in.size(),KetList.size()));
   vector<long long unsigned int> index_list (min(vec_in.size(),KetList.size()));
   auto it = set_intersection(KetList.begin(),KetList.end(),vec_in.begin(),vec_in.end(),index_list.begin());
   index_list.resize(it-index_list.begin());
   for (auto& x : index_list)
   {
     x = KetMap[x];
   }
   return arma::uvec(index_list);
}

//************************************************************************

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
bool TwoBodyChannel_CC::CheckChannel_ket(int p, int q) const
{
   Orbit & op = modelspace->GetOrbit(p);
   Orbit & oq = modelspace->GetOrbit(q);
   if ((op.l + oq.l)%2 != parity)    return false;
   if (abs(op.tz2 + oq.tz2) != 2*Tz) return false;
   if (op.j2 + oq.j2 < 2*J)          return false;
   if (abs(op.j2 - oq.j2) > 2*J)     return false;

   return true;
}


//************************************************************************
//************************************************************************

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
{
   cout << "In copy constructor" << endl;
   norbits = ms.norbits;
   hbar_omega = ms.hbar_omega;
   target_mass = ms.target_mass;
   OneBodyJmax = ms.OneBodyJmax;
   TwoBodyJmax = ms.TwoBodyJmax;
   ThreeBodyJmax = ms.ThreeBodyJmax;
   Nmax = ms.Nmax;
   N2max = ms.N2max;
   N3max = ms.N3max;

   holes = ms.holes;
   particles = ms.particles;
   valence = ms.valence;
   qspace = ms.qspace;
   hole_qspace = ms.hole_qspace;
   particle_qspace = ms.particle_qspace;
   proton_orbits = ms.proton_orbits;
   neutron_orbits = ms.neutron_orbits;
   Orbits = ms.Orbits;
   Kets = ms.Kets;
   TwoBodyChannels = ms.TwoBodyChannels;
   TwoBodyChannels_CC = ms.TwoBodyChannels_CC;

   for (Ket& k : Kets)   k.ms = this;
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
}

ModelSpace::ModelSpace(ModelSpace&& ms)
 : norbits(ms.norbits), hbar_omega(ms.hbar_omega), target_mass(ms.target_mass),
   OneBodyJmax(ms.OneBodyJmax), TwoBodyJmax(ms.TwoBodyJmax), ThreeBodyJmax(ms.ThreeBodyJmax),
   Nmax(ms.Nmax), N2max(ms.N2max), N3max(ms.N3max),
   holes( move(ms.holes)), particles( move(ms.particles)), valence(move(ms.valence)),
   qspace( move(ms.qspace)), hole_qspace(move(ms.hole_qspace)), proton_orbits( move(ms.proton_orbits)),
   neutron_orbits( move(ms.neutron_orbits)), Orbits(move(ms.Orbits)), Kets(move(ms.Kets)),
   TwoBodyChannels(move(ms.TwoBodyChannels)), TwoBodyChannels_CC(move(ms.TwoBodyChannels_CC))
{

   for (Ket& k : Kets)   k.ms = this;
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
}


// orbit string representation is e.g. p0f7
ModelSpace::ModelSpace(int nmax, vector<string> hole_list, vector<string> inside_list)
: Nmax(nmax), N2max(2*nmax), N3max(3*nmax)
{
   Init(nmax,hole_list,inside_list);
}

// Shortcuts for common modelspaces
ModelSpace::ModelSpace(int nmax, string str)
: Nmax(nmax), N2max(2*nmax), N3max(3*nmax)
{
  if (str == "skeleton") Init_Skeleton(nmax);
  else if (str == "He4") Init_He4(nmax);
  else if (str == "O16") Init_O16(nmax);
  else if (str == "Ca40") Init_Ca40(nmax);
  else if (str == "p-shell") Init_PShell(nmax);
  else if (str == "sd-shell") Init_SDShell(nmax);
  else if (str == "psd-shell") Init_PSDShell(nmax);
  else if (str == "o16-psd-shell") Init_O16PSDShell(nmax);
  else if (str == "fp-shell") Init_FPShell(nmax);
  else if (str == "sdfp-shell") Init_SDFPShell(nmax);
  else cout << "No such pre-configured model space: " << str << endl;
}

//ModelSpace::ModelSpace(int nmax, vector<string> hole_list, vector<string> inside_list)
void ModelSpace::Init(int nmax, vector<string> hole_list, vector<string> inside_list)
{
   OneBodyJmax = 0;
   TwoBodyJmax = 0;
   ThreeBodyJmax = 0;
   hbar_omega=20;
   target_mass = 16;
   Nmax = nmax;
   N2max = 2*Nmax;
   N3max = 3*Nmax;

   cout << "Creating a model space with Nmax = " << Nmax << "  and hole orbits [";
   for (string& h : hole_list)
   {
       cout << h << " ";
   }
   cout << "]   and valence space [";
   for (string& h : inside_list)
   {
       cout << h << " ";
   }
   cout << "]" << endl;
   
   vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t'};
   vector<char> pn_list = { 'p', 'n' };

   cout << "Generating orbits, etc." << endl;
   norbits = (Nmax+1)*(Nmax+2);
   for (int N=0; N<=Nmax; ++N)
   {
     for (int l=N; l>=0; l-=2)
     {
       int n = (N-l)/2;
       for (int j2=2*l+1; j2>=2*l-1 and j2>0; j2-=2)
       {
         for (int tz=-1; tz<=1; tz+=2)
         {
            int ph = 0;
            int io = 1;
            double spe = 0;
            char orb_string[6];
            sprintf(orb_string, "%c%i%c%i", pn_list[(tz+1)/2], n, l_list[l], j2);
            string orb_str = orb_string;
            auto it_hole = find(hole_list.begin(), hole_list.end(), orb_string);
            if ( it_hole != hole_list.end() )
            {
               ph=1;
               hole_list.erase(it_hole);
            }
            auto it_inside = find(inside_list.begin(), inside_list.end(), orb_string);
            if ( it_inside != inside_list.end() )
            {
               io=0;
               inside_list.erase(it_inside);
            }
            AddOrbit(Orbit(n,l,j2,tz,ph,io,spe));
         }
       }
     }
   }
   cout << "Setting up kets" << endl;
   SetupKets();
   cout << "Done with Init()" << endl;
}



void ModelSpace::Init_Skeleton(int nmax)
{

   OneBodyJmax = 0;
   TwoBodyJmax = 0;
   ThreeBodyJmax = 0;
   hbar_omega=20;
   target_mass = 16;
   Nmax = nmax;
   N2max = 2*Nmax;
   N3max = 3*Nmax;

//   cout << "Creating a skeleton model space with Nmax = " << Nmax << endl;
   
   norbits = (Nmax+1)*(Nmax+2);
   for (int N=0; N<=Nmax; ++N)
   {
     for (int l=N; l>=0; l-=2)
     {
       int n = (N-l)/2;
       for (int j2=2*l+1; j2>=2*l-1 and j2>0; j2-=2)
       {
         for (int tz=-1; tz<=1; tz+=2)
         {
            int ph = 0;
            int io = 1;
            double spe = 0;
            AddOrbit(Orbit(n,l,j2,tz,ph,io,spe));
         }
       }
     }
   }
}

// Some of the more common model spaces, for convenience.
void ModelSpace::Init_He4(int nmax)
{
   vector<string> core = {"p0s1","n0s1"};
   vector<string> valence = {};
   Init(nmax,core,valence);
}

void ModelSpace::Init_O16(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
   vector<string> valence = {};
   Init(nmax,core,valence);
}

void ModelSpace::Init_Ca40(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   vector<string> valence = {};
   Init(nmax,core,valence);
}

void ModelSpace::Init_PShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1"};
   vector<string> valence = {"p0p3","n0p3","p0p1","n0p1"};
   Init(nmax,core,valence);
}

void ModelSpace::Init_SDShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
   vector<string> valence = {"p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   Init(nmax,core,valence);
}

void ModelSpace::Init_PSDShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1"};
   vector<string> valence = {"p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   Init(nmax,core,valence);
}

void ModelSpace::Init_O16PSDShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
   vector<string> valence = {"p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   Init(nmax,core,valence);
}

void ModelSpace::Init_FPShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"};
   vector<string> valence = {"p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"};
   Init(nmax,core,valence);
}

void ModelSpace::Init_SDFPShell(int nmax)
{
   vector<string> core = {"p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"};
   vector<string> valence = {"p0d5","n0d5","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"};
   Init(nmax,core,valence);
}

ModelSpace ModelSpace::operator=(const ModelSpace& ms)
{
   return ModelSpace(ms);
}
ModelSpace ModelSpace::operator=(ModelSpace&& ms)
{
   return ModelSpace(ms);
}

void ModelSpace::AddOrbit(Orbit orb)
{
   int ind = Index1(orb.n, orb.l, orb.j2, orb.tz2);
   if (Orbits.size() <= ind) Orbits.resize(ind+1,Orbit());
   Orbits[ind] = Orbit(orb);
   norbits = Orbits.size();
   if (orb.j2 > OneBodyJmax)
   {
      OneBodyJmax = orb.j2;
      TwoBodyJmax = OneBodyJmax*2;
      ThreeBodyJmax = OneBodyJmax*3-1;
      nTwoBodyChannels = 2*3*(TwoBodyJmax+1);
   }
   if ( 2*orb.n+orb.l > Nmax )
   {
      Nmax = 2*orb.n+orb.l;
      N2max = 2*Nmax;
      N3max = 3*Nmax;
   }

   if (orb.ph == 0) particles.push_back(ind);
   if (orb.ph == 1) holes.push_back(ind);
   if (orb.io == 0) valence.push_back(ind);
   if (orb.io == 1)
   {
     qspace.push_back(ind);
     if (orb.ph == 0) particle_qspace.push_back(ind);
     if (orb.ph == 1) hole_qspace.push_back(ind);
   }
   if (orb.tz2<0) proton_orbits.push_back(ind);
   if (orb.tz2>0) neutron_orbits.push_back(ind);

   OneBodyChannels[{orb.l, orb.j2, orb.tz2}].push_back(ind);

}

int ModelSpace::GetTwoBodyChannelIndex(int j, int p, int t)
{
   return (t+1)*2*(TwoBodyJmax) + p*(TwoBodyJmax) + j;
}


void ModelSpace::SetupKets()
{
   int index = 0;

   for (int p=0;p<norbits;p++)
   {
     for (int q=p;q<norbits;q++)
     {
        index = Index2(p,q);
        if (index >= Kets.size()) Kets.resize(index+1,Ket(this));
        Kets[index].Setpq(p,q);
     }
   }

  for (int index=0;index<Kets.size();++index)
  {
    Ket& ket = Kets[index];
    int php = GetOrbit(ket.p).ph;
    int phq = GetOrbit(ket.q).ph;
    int iop = GetOrbit(ket.p).io;
    int ioq = GetOrbit(ket.q).io;
     if (( php + phq)==2) // hh
     {
        KetIndex_hh.push_back(index);
        if ((iop+ioq)==2) // qq
        {
           KetIndex_holeq_holeq.push_back(index);
        }
     }
     else if ((php + phq) == 0) // pp
     {
        KetIndex_pp.push_back(index);
        if ((iop+ioq)==2) // qq
        {
           KetIndex_particleq_particleq.push_back(index);
        }
     }
     else //ph
     {
        KetIndex_ph.push_back(index);
        if ((iop+ioq)==2) // qq
        {
           KetIndex_particleq_holeq.push_back(index);
        }
     }
     if ((iop + ioq) == 0) // vv
     {
        KetIndex_vv.push_back(index);
     }

     if ((iop + ioq) == 1) // vq
     {
       if ((iop + php == 2) or (ioq+phq==2) ) // the qspace orbit is a hole
          KetIndex_v_holeq.push_back(index);
       else // v particle_q
          KetIndex_v_particleq.push_back(index);
     }

   }

   cout << "Set up TwoBodyChannels" << endl;
   for (int ch=0;ch<nTwoBodyChannels;++ch)
   {
      TwoBodyChannels.push_back(TwoBodyChannel(ch,this));
      TwoBodyChannels_CC.push_back(TwoBodyChannel_CC(ch,this));
   }
}




double ModelSpace::GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3)
{
// { j1 j2 j3 }
// { J1 J2 J3 }

   unsigned long int key = 20000000000*j1 + 200000000*j2 + 2000000*j3 + 20000*J1 + 200*J2 + 2*J3;

   auto it = SixJList.find(key);
   if (it != SixJList.end() ) return it->second;
   double sixj = AngMom::SixJ(j1,j2,j3,J1,J2,J3);
   #pragma omp critical
   SixJList[key] = sixj;
   return sixj;
}



void ModelSpace::PreCalculateMoshinsky()
{
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


//          unsigned long long int key =  1000000000000 * N
//                                       + 100000000000 * Lam
//                                       +   1000000000 * n
//                                       +    100000000 * lam
//                                       +      1000000 * n1
//                                       +       100000 * l1
//                                       +         1000 * n2
//                                       +          100 * l2
//                                       +                 L;

//          if (N==3 and n==3 and Lam==0 and lam==0)
//          cout << N << " " << Lam << " " << n << " " << lam << " " << n1 << " " << l1 << " " << n2 << " " << l2 << " " << L << endl;
          double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
//          MoshList[key] = mosh;
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
   cout << "Shouldn't be here..." << N << " " << Lam << " " <<  n << " " << lam << " " << n1 << " " << l1 << " " << n2 << " " << l2 << " " << L << endl;
//   #pragma omp atomic
   MoshList[key] = mosh;
   return mosh * phase_mosh;

}


/*

double ModelSpace::GetMoshinsky( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L)
{
   unsigned long long int key =  1000000000000 * N
                                + 100000000000 * Lam
                                +   1000000000 * n
                                +    100000000 * lam
                                +      1000000 * n1
                                +       100000 * l1
                                +         1000 * n2
                                +          100 * l2
                                +                 L;
   map<long int,double>::iterator it = MoshList.find(key);
   if ( it != MoshList.end() )  return it->second;

   // if we didn't find it, we need to calculate it.
   double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
   #pragma omp critical
   MoshList[key] = mosh;
   return mosh;

}
*/


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


