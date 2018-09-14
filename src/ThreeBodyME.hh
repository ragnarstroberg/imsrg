///////////////////////////////////////////////////////////////////////////////////
//    ThreeBodyME.hh, part of  imsrg++
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

#ifndef ThreeBodyME_h
#define ThreeBodyME_h 1

#include "ModelSpace.hh"
#include <fstream>
#include <unordered_map>

//typedef double ThreeBME_type;
typedef float ThreeBME_type;

/// The three-body piece of an operator, stored in nested vectors.
/// The 3BMEs are stored in unnormalized JT coupled form
/// \f$ \langle (abJ_{ab}t_{ab})c | V | (deJ_{de}t_{de})f \rangle_{JT} \f$.
/// To minimize the number of stored matrix elements, only elements with
/// \f$ a\geq b \geq c, a\geq d\geq e \geq f \f$ are stored.
/// The other combinations are obtained on the fly by GetME().
/// The storage format is MatEl[{a,b,c,d,e,f,J,Jab,Jde}][T_index] =
/// \f$ \langle (abJ_{ab}t_{ab})c | V | (deJ_{de}t_{de})f  \rangle_{JT} \f$.
class ThreeBodyME
{
 public:
  ModelSpace * modelspace;
  std::vector<ThreeBME_type> MatEl;
  std::unordered_map<size_t, size_t> OrbitIndexHash; //
  int E3max;
  size_t total_dimension;
  const static int ABC;
  const static int BCA;
  const static int CAB;
  const static int ACB;
  const static int BAC;
  const static int CBA;
  
  ~ThreeBodyME();
  ThreeBodyME();
  ThreeBodyME(ModelSpace*);
  ThreeBodyME(ModelSpace* ms, int e3max);

  size_t KeyHash(size_t,size_t,size_t,size_t,size_t,size_t) const;
  void Allocate();

  void SetModelSpace(ModelSpace *ms){modelspace = ms;};

//// Three body setter getters
  std::vector<std::pair<size_t,double>> AccessME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n) const;
  ThreeBME_type AddToME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, ThreeBME_type V);
  void   SetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, ThreeBME_type V);
  ThreeBME_type GetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n) const;
  ThreeBME_type GetME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n) const;

///// Some other three body methods

  int SortOrbits(int a_in, int b_in, int c_in, int& a,int& b,int& c) const;
  inline double RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int J) const;
  void SetE3max(int e){E3max = e;};
  int GetE3max(){return E3max;};

  void Erase(); // set all three-body terms to zero
  void Deallocate();
  size_t size(){return total_dimension * sizeof(ThreeBME_type);};


  void WriteBinary(std::ofstream&);
  void ReadBinary(std::ifstream&);

};


#endif
