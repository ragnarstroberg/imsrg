///////////////////////////////////////////////////////////////////////////////////
//    AngMom.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef AngMom_hh
#define AngMom_hh 1

//#include <cmath>
#include <vector>
#include <array>

//#define MOSH_BETA_1 M_PI_4   // Moshinksy beta parameter for mass ratio 1:1
//#define MOSH_BETA_2 atan(sqrt(2)) // Moshinsky beta parameter for  mass ratio m1/m2 = 2
//#define MOSH_BETA_half atan(sqrt(0.5)) // Moshinsky beta parameter for mass ration m1/m2 = 0.5


//class AngMom
namespace AngMom
{

   extern double mosh_beta_1;
   extern double mosh_beta_2;
   extern double mosh_beta_half;

   static std::vector<double> FactorialList;
   static std::vector<double> DoubleFactorialList;
   void FillFactorialLists(int nmax);
   double factorial(int x);
   double double_fact(int x);

   int phase(int x);
   double Tri(double j1, double j2, double j3);
   bool Triangle(double j1, double j2, double j3);
   int Jmin( std::vector< std::array<int,2>> jpairs);
   int Jmax( std::vector< std::array<int,2>> jpairs);
   double CG(double ja, double ma, double jb, double mb, double J, double M);
   double ThreeJ(double j1, double j2, double j3, double m1, double m2, double m3);
   double SixJ(double j1, double j2, double j3, double J1, double J2,double J3);
   double NineJ(double j1, double j2, double j3, double j4, double j5, double j6, double j7, double j8, double j9);
   double NormNineJ(double j1, double j2, double j3, double j4, double j5, double j6, double j7, double j8, double j9);
   double TwelveJ_1_ints(int a1, int a2, int a3, int a4, int b12, int b23, int b34, int b41, int c1, int c2, int c3, int c4);
   double TwelveJ_1(double a1, double a2, double a3, double a4, double b12, double b23, double b34, double b41, double c1, double c2, double c3, double c4);
   double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam, double B=AngMom::mosh_beta_1);
//   double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam, double B=3.1415926/2.0);
//   double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam, double B=MOSH_BETA_1);
   double TalmiB(int n, int l, int nn, int ll, int p);
//   double Tcoeff( LabKet& labket, int Jab, int twoJ, jacobi1_state& jac1, jacobi2_state& jac2, int twoJ12, int Ncm, int Lcm);
   double Tcoeff( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm);
   double Tcoeff_bruteforce( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm);
};

#endif

