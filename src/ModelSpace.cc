#include "ModelSpace.hh"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

Orbit::Orbit()
{
}

Orbit::Orbit(int nn, int ll, int jj2, int ttz2, int hhvq, float e=0.0)
{
   n = nn;
   l = ll;
   j2 = jj2;
   tz2 = ttz2;
   hvq = hhvq;
   spe = e;
}

Orbit::Orbit(const Orbit& orb)
{
  n = orb.n;
  l = orb.l;
  j2 = orb.j2;
  tz2 = orb.tz2;
  hvq = orb.hvq;
  spe = orb.spe;
}

void Orbit::Set(int nn, int ll, int jj2, int ttz2, int hhvq, float e)
{
   n = nn;
   l = ll;
   j2 = jj2;
   this->tz2 = ttz2;
   hvq = hhvq;
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
}

int Ket::Phase(int J)
{
   int exponent = (ms->GetOrbit(p)->j2 + ms->GetOrbit(q)->j2)/2 + J + 1;
   return 1- 2*(exponent%2);
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
   KetMap.resize(nk,-1);
   for (int i=0;i<nk;i++)
   {
      Ket *ket = modelspace->GetKet(i);
      if ( CheckChannel_ket(ket) )
      {
         KetMap[i] = NumberKets;
         KetList.push_back(i);
         NumberKets++;
      }
   }
   Proj_pp = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   Proj_hh = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   //for (int &i: KetList) // C++11 syntax
   for (int i=0;i<NumberKets;i++)
   {
      Ket *ket = GetKet(i);
      if ( modelspace->GetOrbit(ket->p)->hvq ==0 and modelspace->GetOrbit(ket->q)->hvq==0)
      {
         Proj_hh(i,i) = 1;
      }
      if ( modelspace->GetOrbit(ket->p)->hvq >0 and modelspace->GetOrbit(ket->q)->hvq>0)
      {
         Proj_pp(i,i) = 1;
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

ModelSpace::ModelSpace()
{
//   nCore = 0;
   norbits = 0;
   maxj = 0;
   hbar_omega=20;
   target_mass = 16;
}


ModelSpace::ModelSpace(const ModelSpace& ms)
{
//   nCore = 0;
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
//   hole = hole;
//   particle=particle;
//   valence = valence;
//   qspace = qspace;
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
      //nTwoBodyChannels = 2*3*(maxj+1);
      nTwoBodyChannels = 2*3*(JMAX);
   }
//   if (orb.hvq == 0) nCore+=orb.j2+1; // 0 means hole (ie core), 1 means valence, 2 means outside the model space
   if (orb.hvq == 0) hole.push_back(ind);
   if (orb.hvq == 1) valence.push_back(ind);
   if (orb.hvq == 2) qspace.push_back(ind);
   if (orb.hvq >0) particles.push_back(ind);
}



void ModelSpace::SetupKets()
{
   int index = 0;

   for (int p=0;p<norbits;p++)
   {
     for (int q=p;q<norbits;q++)
     {
        //index = q*(q+1)/2 + p;
        //index = q*norbits + p;
        index = Index2(p,q);
        if (index >= Kets.size()) Kets.resize(index+1,NULL);
        Kets[index] = new Ket(this,p,q);
     }
   }
   for (int ch=0;ch<nTwoBodyChannels;++ch)
   {
      TwoBodyChannels.push_back(TwoBodyChannel(ch,this));
   }
}
