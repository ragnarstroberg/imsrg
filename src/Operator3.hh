
#ifndef Operator3_h
#define Operator3_h 1

#include "Operator.hh"
#include <unordered_map>

using namespace std;

class Operator3 : public Operator
{
 public:
   Operator3();
   Operator3(ModelSpace&);
   void Copy(const Operator3&);

   void AllocateThreeBody();

   double GetThreeBodyME(int Jab, int Jde, int J2, int tab, int tde, int T2, int i, int j, int k, int l, int m, int n);
   void SetThreeBodyME(int Jab, int Jde, int J2, int tab, int tde, int T2, int i, int j, int k, int l, int m, int n, double V);
   Operator DoNormalOrdering3();
   int GetThreeBodyChannelIndex(int J, int parity, int tab, int tde, int T2);
   long int GetThreeBodyOrbitIndex(int a, int b, int c, int d, int e, int f);

   void SortWithPhase(int& a,int& b,int& c,int& d,int& e, int& f, int& phase);

   // the first vector runs over K2, Tz2 and parity, the vector runs over J
//   vector<vector<arma::mat>> ThreeBody;

   // first index is orbit_index, second index is J_index, 3rd is isospin_index
//   vector< vector< vector<double> > > ThreeBody;
   // first index is channel index, containing (J,T,tab,tde,parity)
   // second index is orbit index, containting (a,b,c,d,e,f)
   // third index is J2 index, containing (Jab,Jde)
//   vector< map<long, vector<double> > > ThreeBodyMatrixElements;
   // (Note, probably better to use unordered_map, i.e. hash table)
//   vector< unordered_map<long, vector<double> > > ThreeBodyMatrixElements;
//   unordered_map<int,unordered_map<long,vector<double> > > ThreeBody;
   vector<unordered_map<long,vector<double> > > ThreeBody;
////   unordered_map<vector<int>, double> ThreeBodyMatrixElements;

};


#endif
