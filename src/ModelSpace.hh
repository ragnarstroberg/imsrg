
#include <vector>

#ifndef ModelSpace_h
#define ModelSpace_h 1


class Orbit;
class Ket;
class ModelSpace;

class Orbit
{
 public:
   int n;
   int l;
   int j2;
   int tz2;
   int hvq; // h=0, v=1, q=-1
   float spe;

   Orbit();
   Orbit(int,int,int,int,int);
   Orbit(int,int,int,int,int,float);
   Orbit(const Orbit&);
   void Set(int n, int l, int j2, int tz2, int hvq, float spe);
};

class Ket  //  | pq >
{
 public:
   int p;
   int q;
   ModelSpace * ms;
   int parity;
   int Tz;
   int Jmin;
   int Jmax;
   int Jstep;
   int Phase(int J);
   int delta_pq(){return p==q?1:0;};
   Ket();
   Ket(ModelSpace * ms, int p, int q);
};

class ModelSpace
{

 public:
   // Fields
   int norbits;
   int nCore;
   int maxj;
   float hbar_omega;
   int target_mass;
   std::vector<Orbit*> Orbits;
   std::vector<Ket*> Kets;
   // Constructors
   ModelSpace();
   ModelSpace(ModelSpace*);
   // Overloaded operators
//   ModelSpace operator=(const ModelSpace&); 

   Orbit* GetOrbit(int i) {return Orbits[i];}; 
   Ket* GetKet(int i) {return Kets[i];};
   Ket* GetKet(int p, int q) {return Kets[GetKetIndex(p,q)];};
   int GetKetIndex(int p, int q) {return q*(q+1)/2+p;}; // convention is a<=b
   int GetKetIndex(Ket * ket) {return ket->q*(ket->q+1)/2+ket->p;}; // convention is a<=b
   int GetNumberOrbits() {return norbits;};
   int GetNumberKets() {return Kets.size();};
   void AddOrbit(Orbit orb);
   void SetupKets();

   int Index1(int n, int l, int j2, int tz2){return(2*n+l)*(2*n+l+3) + 1-j2 + (tz2+1)/2 ;};
//   int Index1(int n, int l, int j2, int tz2){return(2*n+l)*(2*n+l+1) + j2-1 + (tz2+1)/2 ;};
//   int Index2(int a, int b, int J){return J*norbits*norbits + a*norbits + b;};
   int Index2(int p, int q){return p*norbits + q;};

 private:
   int maxn;
   int maxl;
   int maxjw;
   int maxtz;


};




#endif

