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
 public:

  std::string generator_type;
  Operator * H;
  Operator * Eta;
  ModelSpace* modelspace;
  double denominator_cutoff;
  double denominator_delta;
  int denominator_delta_index;

  enum denominator_partitioning_t{Epstein_Nesbet,Moller_Plesset};
  denominator_partitioning_t denominator_partitioning;

  Operator RspaceRegulator;
  double regulator_length;

  bool only_2b_eta; // even if we're doing IMSRG(3), keep eta to 2b

  // The functional form dictating what to do with Hod and a denominator
  static std::function<double(double,double)> wegner_func;
  static std::function<double(double,double)> white_func;
  static std::function<double(double,double)> atan_func;
  static std::function<double(double,double)> imaginarytime_func;
  static std::function<double(double,double)> qtransferatan1_func;


  Generator();
  void SetType(std::string g){generator_type = g;};
  void SetDenominatorPartitioning(std::string dp){if (dp=="Moller_Plesset"){denominator_partitioning=Moller_Plesset;}else {denominator_partitioning=Epstein_Nesbet;};};
  std::string GetType(){return generator_type;};
  void Update(Operator* H, Operator* Eta);
  void AddToEta(Operator* H, Operator* Eta);
  void SetDenominatorCutoff(double c){denominator_cutoff=c;};
  void SetDenominatorDelta(double d){denominator_delta=d;};
  void SetDenominatorDeltaIndex(int i){denominator_delta_index=i;};
  void SetDenominatorDeltaOrbit(std::string orb);

  void SetRegulatorLength(double r);

// private:
//  void ConstructGenerator_Wegner();
//  void ConstructGenerator_White();

  void ConstructGenerator_SingleRef(std::function<double (double,double)>& etafunc );
  void ConstructGenerator_SingleRef_3body(std::function<double (double,double)>& etafunc );
//  void ConstructGenerator_Atan();
//  void ConstructGenerator_Atan_3body();
//  void ConstructGenerator_ImaginaryTime();
//  void ConstructGenerator_ImaginaryTime_3body();
//  void ConstructGenerator_QTransferAtan(int n);
//  void ConstructGenerator_ShellModel();
  void ConstructGenerator_ShellModel(std::function<double (double,double)>& eta_func);
  void ConstructGenerator_ShellModel_3body(std::function<double (double,double)>& eta_func);
//  void ConstructGenerator_ShellModel_Atan();
//  void ConstructGenerator_ShellModel_Atan_3body();
//  void ConstructGenerator_ShellModel_Wegner();
//  void ConstructGenerator_ShellModel_ImaginaryTime();
//  void ConstructGenerator_ShellModel_ImaginaryTime_3body();
//  void ConstructGenerator_ShellModel_Atan_NpNh();
  void ConstructGenerator_ShellModel_NpNh(std::function<double(double,double)>& eta_func);
  void ConstructGenerator_HartreeFock();
//  void ConstructGenerator_1PA();
  void ConstructGenerator_1PA(std::function<double(double,double)>& eta_func);
  void ConstructGenerator_Rspace();
  void SetOnly2bEta(bool tf){only_2b_eta = tf;};
  double Get1bDenominator(int i, int j);
  double Get2bDenominator(int ch, int ibra, int iket);
  double Get2bDenominator_Jdep(int ch, int ibra, int iket);
  double Get3bDenominator(int i, int j, int k, int l, int m, int n);

  

};

#endif
