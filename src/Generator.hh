///////////////////////////////////////////////////////////////////////////////////
//    Generator.hh, part of  imsrg++
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

#ifndef Generator_hh
#define Generator_hh 1

#include "ModelSpace.hh"
#include "Operator.hh"

#include <string>


class Generator
{
 protected: //DK: Was private. Changed to inherit into GeneratorPV 
  Operator * H;
  Operator * Eta;

 public:

  std::string generator_type;
//  ModelSpace* modelspace;
  double denominator_cutoff;
  double denominator_delta;
  int denominator_delta_index;

  enum denominator_partitioning_t{Epstein_Nesbet,Moller_Plesset,MP_isospin};
  denominator_partitioning_t denominator_partitioning;

  double regulator_length;

  bool only_2b_eta; // even if we're doing IMSRG(3), keep eta to 2b
  bool use_isospin_averaging;

  // The functional form dictating what to do with Hod and a denominator
  static std::function<double(double,double)> wegner_func;
  static std::function<double(double,double)> white_func;
  static std::function<double(double,double)> atan_func;
  static std::function<double(double,double)> imaginarytime_func;
  static std::function<double(double,double)> qtransferatan1_func;


  Generator();
  void SetType(std::string g){generator_type = g;};
  void SetDenominatorPartitioning(std::string dp); 
  std::string GetType(){return generator_type;};
  void Update(Operator& H_s, Operator& Eta_s);
  void AddToEta(Operator& H_s, Operator& Eta_s);
  void SetDenominatorCutoff(double c){denominator_cutoff=c;};
  void SetDenominatorDelta(double d){denominator_delta=d;};         // call SetDenominatorDeltaIndex(-12345) to use it
  void SetDenominatorDeltaIndex(int i){denominator_delta_index=i;};  
  void SetDenominatorDeltaOrbit(std::string orb);
  void SetUseIsospinAveraging( bool tf ){use_isospin_averaging=tf;};

  Operator GetHod(Operator& H);



  void ConstructGenerator_SingleRef(std::function<double (double,double)>& etafunc );
  void ConstructGenerator_SingleRef_3body(std::function<double (double,double)>& etafunc );
  void ConstructGenerator_ShellModel(std::function<double (double,double)>& eta_func);
  void ConstructGenerator_ShellModel_3body(std::function<double (double,double)>& eta_func);
  void ConstructGenerator_ShellModel_NpNh(std::function<double(double,double)>& eta_func);
  void ConstructGenerator_HartreeFock();
  void ConstructGenerator_1PA(std::function<double(double,double)>& eta_func);
  void SetOnly2bEta(bool tf){only_2b_eta = tf;};
  double Get1bDenominator(int i, int j);
  double Get2bDenominator(int ch, int ibra, int iket) { return Get2bDenominator(ch,ch,ibra,iket);};
  double Get2bDenominator(int ch_bra, int ch_ket, int ibra, int iket);
  double Get2bDenominator_Jdep(int ch, int ibra, int iket);
  double Get3bDenominator(int i, int j, int k, int l, int m, int n);

  Operator GetHod_SingleRef( Operator& H );
  Operator GetHod_ShellModel( Operator& H );
  

};


/// Additional helper functions which allow us to conveniently combine sets of orbits.
/// For example, if we want a list of orbits which are either valence or qspace, we do VectorUnion( valence, qspace )
// Templated functions need to be defined in the header file (or else explicitly declared in the .cc file).

/// Base case with one argument
 template <typename T>
 T VectorUnion(const T& v1)
 {
   return v1;
 }
 
// Variadic template to accept an arbitrary number of arguments
 template <typename T, typename... Args>
 T VectorUnion(const T& v1, const T& v2, Args... args)
 {
   T vec(v1.size()+v2.size());
   std::copy(v1.begin(),v1.end(),vec.begin());
   std::copy(v2.begin(),v2.end(),vec.begin()+v1.size());
   return VectorUnion(vec, args...);
 }

// If we use a std::set as the containter, we need to use insert() rather than copy().
 template <typename T, typename... Args>
 std::set<T> VectorUnion(const std::set<T>& s1, const std::set<T>& s2, Args... args)
 {
   std::set<T> s3 = s1;
   s3.insert(s2.begin(),s2.end());
   return VectorUnion( s3, args...);
 }

#endif
