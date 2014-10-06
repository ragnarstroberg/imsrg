#ifndef ModelSpace_h
#define ModelSpace_h 1

#include <vector>

class Orbit;
class Ket;
class ModelSpace;

class Orbit
{
 public:
   // Fields
   int n;
   int l;
   int j2;
   int tz2;
   int hvq; // hole=0, valence=1, qspace=2
   float spe;

   //Constructors
   Orbit();
   Orbit(int n ,int l, int j, int t, int hvq ,float spe);
   Orbit(const Orbit&);
   // Methods
   void Set(int n, int l, int j2, int tz2, int hvq, float spe);
};


class Ket  //  | pq >
{
 public:
   // Fields
   int p;
   int q;
   int parity;
   int Tz;
   int Jmin;
   int Jmax;
   int Jstep;
   // Constructors
   Ket(){};
   Ket(ModelSpace * ms, int p, int q);
   // Methods
   int Phase(int J);
   int delta_pq(){return p==q ? 1 : 0;};

 private:
   // Fields
   ModelSpace * ms;

};

class ModelSpace
{

 public:

   // Constructors
   ModelSpace();
   ModelSpace(const ModelSpace&);
   // Overloaded operators
   ModelSpace operator=(const ModelSpace&); 
   // Methods
   void SetupKets();
   Orbit* GetOrbit(int i) const {return Orbits[i];}; 
   Ket* GetKet(int i) const {return Kets[i];};
   Ket* GetKet(int p, int q) const {return Kets[GetKetIndex(p,q)];};
   int GetKetIndex(int p, int q) const {return q*(q+1)/2+p;}; // convention is p<=q
   int GetKetIndex(Ket * ket) const {return ket->q*(ket->q+1)/2+ket->p;}; // convention is p<=q
   int GetNumberOrbits() const {return norbits;};
   int GetNumberKets() const {return Kets.size();};
   void AddOrbit(Orbit orb);
   void SetHbarOmega(float hw) {hbar_omega = hw;};
   void SetTargetMass(int A) {target_mass = A;};
   float GetHbarOmega() const {return hbar_omega;};
   int GetTargetMass() const {return target_mass;};
   int Index1(int n, int l, int j2, int tz2){return(2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2 ;};
   int Index2(int p, int q){return p*norbits + q;};

 private:
   // Fields
   int norbits;
   int nCore;
   int maxj;
   float hbar_omega;
   int target_mass;
   std::vector<Orbit*> Orbits;
   std::vector<Ket*> Kets;

};




#endif

