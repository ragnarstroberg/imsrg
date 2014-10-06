#include "ModelSpace.hh"
#include <iostream>
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




ModelSpace::ModelSpace()
{
   nCore = 0;
   norbits = 0;
   maxj = 0;
   hbar_omega=20;
   target_mass = 16;
}


ModelSpace::ModelSpace(const ModelSpace& ms)
{
   nCore = 0;
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
   if (orb.hvq == 0) nCore+=orb.j2+1; // 0 means hole (ie core), 1 means valence, 2 means outside the model space
}



void ModelSpace::SetupKets()
{
   maxj = 0;
   int index = 0;
   for (int i=0;i<norbits;i++)
   {
      if (Orbits[i] == NULL) continue;
      int j2 = GetOrbit(i)->j2;
      if ( j2 > maxj) maxj = j2;
   }

   for (int q=0;q<norbits;q++)
   {
//     int tzb = GetOrbit(q)->tz2;
//     int parb = (GetOrbit(q)->l)%2;
     for (int p=0;p<=q;p++)
     {
        index = q*(q+1)/2 + p;
        if (index >= Kets.size()) Kets.resize(index+1,NULL);
        Kets[index] = new Ket(this,p,q);
     }
   }
}
