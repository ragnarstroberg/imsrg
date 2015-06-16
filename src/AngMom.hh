
#ifndef AngMom_hh
#define AngMom_hh 1

using namespace std;

//class AngMom
namespace AngMom
{
   int phase(int x);
   double Tri(double j1, double j2, double j3);
   bool Triangle(double j1, double j2, double j3);
   double CG(double ja, double ma, double jb, double mb, double J, double M);
   double ThreeJ(double j1, double j2, double j3, double m1, double m2, double m3);
   double SixJ(double j1, double j2, double j3, double J1, double J2,double J3);
   double NineJ(double j1, double j2, double j3, double j4, double j5, double j6, double j7, double j8, double j9);
   double NormNineJ(double j1, double j2, double j3, double j4, double j5, double j6, double j7, double j8, double j9);
   double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam);
};

#endif

