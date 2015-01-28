
#ifndef Operator3_h
#define Operator3_h 1

#include "Operator.hh"
#include <unordered_map>
#include <map>

using namespace std;

//typedef unordered_map<long,vector<double> > orbit_map;

class Operator3 : public Operator
{
 public:
   Operator3();
   Operator3(ModelSpace&);
   Operator3(const Operator3&);
   Operator3(const Operator&);

  //Overloaded operators
  Operator3& operator=( const Operator3& rhs);
  Operator3& operator=( const Operator& rhs);
  Operator3& operator+=( const Operator& rhs);
  Operator3 operator+( const Operator& rhs) const;
  Operator3& operator-=( const Operator& rhs);
  Operator3 operator-( const Operator& rhs) const;
  Operator3& operator*=( const double rhs);
  Operator3 operator*( const double rhs) const;
  Operator3& operator/=( const double rhs);
  Operator3 operator/( const double rhs) const;

   void Copy(const Operator&);
   void Copy(const Operator3&);


   void AllocateThreeBody();

   double GetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n);
   double GetThreeBodyME_pn(int Jab_in, int Jde_in, int J2, int Tz, int i, int j, int k, int l, int m, int n);
   void   SetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, double V);

   Operator DoNormalOrdering3();
   long int GetThreeBodyOrbitIndex(int a, int b, int c, int d, int e, int f);

   void SortOrbits(int& a,int& b,int& c);
   double RecouplingCoefficient(int a_in, int b_in, int c_in, int a, int b, int c, int Jab_in, int Jab, int J ,char j_or_t);
   void SetE3max(int e){E3max = e;};


   // Fields
   //unordered_map< long int, vector< vector<double> > > ThreeBody;
   map< long int, vector< vector<double> > > ThreeBody;
   int E3max;

};


#endif
