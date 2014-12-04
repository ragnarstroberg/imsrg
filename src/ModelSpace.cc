#include "ModelSpace.hh"
#include "AngMom.hh"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

Orbit::Orbit()
{
}

Orbit::Orbit(int nn, int ll, int jj2, int ttz2, int pph, int iio, float e=0.0)
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

void Orbit::Set(int nn, int ll, int jj2, int ttz2, int pph, int iio, float e)
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

Ket::Ket(ModelSpace * modelspace, int pp, int qq)
{
   ms = modelspace;
   p = pp;
   q = qq;
   Orbit * op = ms->GetOrbit(p);
   Orbit * oq = ms->GetOrbit(q);
   parity = (op->l + oq->l)%2;
   Tz = (op->tz2 + oq->tz2)/2;
   Jmin = abs(op->j2 - oq->j2)/2;
   Jmax = (op->j2 + oq->j2)/2;
   Jstep = 1;
   if (p==q) // Pauli principle
   { 
      Jmax--;
      Jstep++;
   }
   phase_prefactor = - ms->phase((op->j2+oq->j2)/2);
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
  Initialize(JMAX*2*t + JMAX*p + j, ms);
}

TwoBodyChannel::TwoBodyChannel(int N, ModelSpace *ms)
{
   Initialize(N,ms);
}

void TwoBodyChannel::Initialize(int N, ModelSpace *ms)
{
   J = N%JMAX;
   parity = (N/JMAX)%2;
   Tz = (N/(2*JMAX)-1);
   modelspace = ms;
   NumberKets = 0;
   int nk = modelspace->GetNumberKets();
   KetMap.resize(nk,-1); // set all values to -1
   for (int i=0;i<nk;i++)
   {
      Ket *ket = modelspace->GetKet(i);
      if ( CheckChannel_ket(ket) )
      {
         KetMap[i] = NumberKets;
         KetList.push_back(i);
         int php = ms->GetOrbit(ket->p)->ph;
         int phq = ms->GetOrbit(ket->q)->ph;
         int iop = ms->GetOrbit(ket->p)->io;
         int ioq = ms->GetOrbit(ket->q)->io;

         if (( php + phq) == 2) // hh
         {
            KetIndex_hh.push_back(NumberKets);
            if ((iop+ioq)==2) // qq
               KetIndex_holeq_holeq.push_back(NumberKets);
         }
         else if ((php + phq) == 0) // pp
         {
            KetIndex_pp.push_back(NumberKets);
            if ((iop+ioq)==2) // qq
               KetIndex_particleq_particleq.push_back(NumberKets);
         }
         else //ph
         {
            KetIndex_ph.push_back(NumberKets);
//            if ((iop+ioq)==2) // qq
//               KetIndex_particleq_holeq.push_back(NumberKets);
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
//            KetIndex_vq.push_back(NumberKets);
         }
//         if ((iop + ioq) == 2) // qq
//         {
//            KetIndex_qq.push_back(NumberKets);
//         }

         NumberKets++;
      }
   }
  
   // Set up projectors which are used in the commutators
   Proj_pp = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   Proj_hh = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   Proj_ph_cc = arma::mat(2*NumberKets, 2*NumberKets, arma::fill::zeros);
   for (int i=0;i<NumberKets;i++)
   {
      Ket *ket = GetKet(i);
      int pha = modelspace->GetOrbit(ket->p)->ph;
      int phb = modelspace->GetOrbit(ket->q)->ph;
      int j2a = modelspace->GetOrbit(ket->p)->j2;
      int j2b = modelspace->GetOrbit(ket->q)->j2;

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
//   KetIndex_qq       = rhs.KetIndex_qq;
//   KetIndex_vq       = rhs.KetIndex_vq;
   KetIndex_vv       = rhs.KetIndex_vv;
   KetIndex_holeq_holeq         = rhs.KetIndex_holeq_holeq;
   KetIndex_particleq_particleq = rhs.KetIndex_particleq_particleq;
   KetIndex_v_holeq     = rhs.KetIndex_v_holeq;
   KetIndex_v_particleq     = rhs.KetIndex_v_particleq;
//   KetIndex_particleq_holeq     = rhs.KetIndex_particleq_holeq;
//   KetIndex_hq       = rhs.KetIndex_hq;
//   KetIndex_vh       = rhs.KetIndex_vh;
}


int TwoBodyChannel::GetLocalIndex(int p, int q) const { return KetMap[modelspace->GetKetIndex(p,q)];}; 

Ket * TwoBodyChannel::GetKet(int i) const { return modelspace->GetKet(KetList[i]);}; // get pointer to ket using local index


bool TwoBodyChannel::CheckChannel_ket(int p, int q) const
{
   if ((p==q) and (J%2 != 0)) return false; // Pauli principle
   Orbit * op = modelspace->GetOrbit(p);
   Orbit * oq = modelspace->GetOrbit(q);
   if ((op->l + oq->l)%2 != parity) return false;
   if ((op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}

//************************************************************************

TwoBodyChannel_CC::TwoBodyChannel_CC()
{}

TwoBodyChannel_CC::TwoBodyChannel_CC(int j, int p, int t, ModelSpace *ms)
{
  Initialize(JMAX*2*t + JMAX*p + j, ms);
}

TwoBodyChannel_CC::TwoBodyChannel_CC(int N, ModelSpace *ms)
{
   Initialize(N,ms);
}



bool TwoBodyChannel_CC::CheckChannel_ket(int p, int q) const
{
   Orbit * op = modelspace->GetOrbit(p);
   Orbit * oq = modelspace->GetOrbit(q);
   if ((op->l + oq->l)%2 != parity) return false;
   if (abs(op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}


//************************************************************************

ModelSpace::ModelSpace()
{
   norbits = 0;
   maxj = 0;
   hbar_omega=20;
   target_mass = 16;
}


ModelSpace::ModelSpace(const ModelSpace& ms)
{
   norbits = 0;
   maxj = 0;
   hbar_omega = ms.hbar_omega;
   target_mass = ms.target_mass;
   int norbits = ms.GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      AddOrbit( Orbit(*ms.GetOrbit(i)) );
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
   if (Orbits.size() <= ind) Orbits.resize(ind+1,NULL);
   Orbits[ind] = (new Orbit(orb));
   norbits = Orbits.size();
   if (orb.j2 > maxj)
   {
      maxj = orb.j2;
      nTwoBodyChannels = 2*3*(JMAX);
   }
   if (orb.ph == 0) particles.push_back(ind);
   if (orb.ph == 1) holes.push_back(ind);
   if (orb.io == 0) valence.push_back(ind);
//   if (orb.io == 1) qspace.push_back(ind);
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
   return (t+1)*2*JMAX + p*JMAX + j;
}


void ModelSpace::SetupKets()
{
   int index = 0;

   for (int p=0;p<norbits;p++)
   {
     for (int q=p;q<norbits;q++)
     {
        index = Index2(p,q);
        if (index >= Kets.size()) Kets.resize(index+1,NULL);
        Kets[index] = new Ket(this,p,q);
     }
   }
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





