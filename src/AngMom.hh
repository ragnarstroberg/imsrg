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

#include <cmath>

#define MOSH_BETA_1 M_PI_4   // Moshinksy beta parameter for mass ratio 1:1
#define MOSH_BETA_2 atan(sqrt(2)) // Moshinsky beta parameter for a mass ratio 2:1

//using namespace std;

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
   double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam, double B=MOSH_BETA_1);
//   double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam);
};

#endif

