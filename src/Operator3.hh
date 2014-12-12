
#ifndef Operator3_h
#define Operator3_h 1

#include "Operator.hh"

using namespace std;

class Operator3 : Operator
{
 public:
   Operator3();
   Operator3(ModelSpace&);


   double GetThreeBodyME(int J, int K2, int Tz2, int parity, int i, int j, int k, int l, int m, int n);
   Operator DoNormalOrdering3();

   // the first vector runs over K2, Tz2 and parity, the vector runs over J
   vector<vector<arma::mat>> ThreeBody;


}


#endif
