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
         int php = ms->GetOrbit(ket.p).ph;
         int phq = ms->GetOrbit(ket.q).ph;
         int iop = ms->GetOrbit(ket.p).io;
         int ioq = ms->GetOrbit(ket.q).io;

         if (( php + phq)==2) // hh
         {
            KetIndex_hh.push_back(NumberKets);
            if ((iop+ioq)==2) // qq
            {
               KetIndex_holeq_holeq.push_back(NumberKets);
            }
         }
         else if ((php + phq) == 0) // pp
         {
            KetIndex_pp.push_back(NumberKets);
            if ((iop+ioq)==2) // qq
            {
               KetIndex_particleq_particleq.push_back(NumberKets);
            }
         }
         else //ph
         {
            KetIndex_ph.push_back(NumberKets);
            if ((iop+ioq)==2) // qq
            {
               KetIndex_particleq_holeq.push_back(NumberKets);
            }
         }


         if ((iop + ioq) == 0) // vv
         {
            KetIndex_vv.push_back(NumberKets);
         }

         if ((iop + ioq) == 1) // vq
         {
           if ((iop + php == 2) or (ioq+phq==2) ) // the qspace orbit is a hole
              KetIndex_v_holeq.push_back(NumberKets);
           else // v particle_q
              KetIndex_v_particleq.push_back(NumberKets);
         }

         NumberKets++;
      }
   }
  
   // Set up projectors which are used in the commutators
   Proj_pp = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   Proj_hh = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   Proj_ph_cc = arma::mat(2*NumberKets, 2*NumberKets, arma::fill::zeros);
   for (int i=0;i<NumberKets;i++)
   {
      Ket &ket = GetKet(i);
      int pha = modelspace->GetOrbit(ket.p).ph;
      int phb = modelspace->GetOrbit(ket.q).ph;
      int j2a = modelspace->GetOrbit(ket.p).j2;
      int j2b = modelspace->GetOrbit(ket.q).j2;

      switch (pha+phb)
      {
         case 0:
           Proj_pp(i,i) = 1;
           break;
         case 2:
           Proj_hh(i,i) = 1;
           break;
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

   KetIndex_pp       = rhs.KetIndex_pp;
   KetIndex_ph       = rhs.KetIndex_ph;
   KetIndex_hh       = rhs.KetIndex_hh;
   KetIndex_vv       = rhs.KetIndex_vv;
   KetIndex_holeq_holeq         = rhs.KetIndex_holeq_holeq;
   KetIndex_particleq_particleq = rhs.KetIndex_particleq_particleq;
   KetIndex_v_holeq     = rhs.KetIndex_v_holeq;
   KetIndex_v_particleq     = rhs.KetIndex_v_particleq;
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

// orbit string representation is e.g. p0f7
ModelSpace::ModelSpace(int nmax, vector<string> hole_list, vector<string> inside_list)
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
   
   vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};
//   map<int,char> pn_list = { (-1,'p'), (1,'n') };
   vector<char> pn_list = { 'p', 'n' };

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
//            sprintf(orb_string, "%c%i%c%i", pn_list[tz], n, l_list[l], j2);
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
   SetupKets();
}



ModelSpace ModelSpace::operator=(const ModelSpace& ms)
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
   for (int ch=0;ch<nTwoBodyChannels;++ch)
   {
      TwoBodyChannels.push_back(TwoBodyChannel(ch,this));
      TwoBodyChannels_CC.push_back(TwoBodyChannel_CC(ch,this));
   }
}

// This needs to be generalized
double ModelSpace::GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3)
{
// { j1 j2 j3 }
// { J1 J2 J3 }
//
// Use 2J in the key so we don't have to worry about half-integers
// Don't really need to store all of them, only need to store
// unique combinations. Since j1,j2,J1,J2 are half-integer
// and j3,J3 are integer, I only swap around the half-integer ones.
   double jlist[4] = {j1,j2,J1,J2};
   int imin = 0;
   double jmin = 9999;
   int k1,k2,k3,K1,K2,K3;
   for (int i=0;i<4;++i)
   {
      if (jlist[i] < jmin)
      {
         imin = i;
         jmin = jlist[i];
      }
   }
   switch (imin)
   {
      case 0:
      k1 = int(2*j1);
      k2 = int(2*j2);
      K1 = int(2*J1);
      K2 = int(2*J2);
      case 1:
      k1 = int(2*j2);
      k2 = int(2*j1);
      K1 = int(2*J2);
      K2 = int(2*J1);
      case 2:
      k1 = int(2*J1);
      k2 = int(2*J2);
      K1 = int(2*j1);
      K2 = int(2*j2);
      case 3:
      k1 = int(2*J2);
      k2 = int(2*J1);
      K1 = int(2*j2);
      K2 = int(2*j1);
   }

   k3 = int(2*j3);
   K3 = int(2*J3);
   long int key = 10000000000*k1 + 100000000*k2 + 1000000*k3 + 10000*K1 + 100*K2 + K3;
   map<long int,double>::iterator it = SixJList.find(key);
   if ( it != SixJList.end() )  return it->second;

   // if we didn't find it, we need to calculate it.
   double sixj = AngMom::SixJ(j1,j2,j3,J1,J2,J3);
   #pragma omp critical
   SixJList[key] = sixj;
   return sixj;

}


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





