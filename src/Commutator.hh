///////////////////////////////////////////////////////////////////////////////////
//    Commutator.hh, part of  imsrg++
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

#ifndef Commutator_h
#define Commutator_h 1

#include "Operator.hh"
#include "TwoBodyME.hh"
#include "ThreeLegME.hh"
#include "TensorCommutators.hh"
#include "IMSRG3Commutators.hh"
#include "DaggerCommutators.hh"
#include "FactorizedDoubleCommutator.hh"
#include "armadillo"
#include <map>
#include <unordered_set>
#include <deque>
#include <array>
#include <string>

namespace Commutator
{

    extern bool single_thread;
    extern bool verbose;

    extern std::map<std::string, bool> comm_term_on; // This allows turning on/off individual IMSRG(3) commutator terms for testing.

    void SetThreebodyThreshold(double x);
    void SetUseIMSRG3(bool tf);
    void SetUseIMSRG3N7(bool tf);
    void SetUseIMSRG3_MP4(bool tf);
    void SetVerbose(bool tf); //Controls if we print out more timing info, useful for optimization work.

    void SetUseIMSRG3N7_Tensor(bool tf);

    void TurnOffTerm(std::string term);
    void TurnOnTerm(std::string term);
    void PrintSettings();

    Operator Commutator(const Operator &X, const Operator &Y);
    Operator CommutatorScalarScalar(const Operator &X, const Operator &Y);
    Operator CommutatorScalarTensor(const Operator &X, const Operator &Y);
    Operator CommutatorScalarDagger(const Operator &X, const Operator &Y);

    void DoPandyaTransformation(const Operator &Z, std::deque<arma::mat> &, std::string orientation);
    void DoPandyaTransformation_SingleChannel(const Operator &Z, arma::mat &X, int ch_cc, std::string orientation);
    void DoPandyaTransformation_SingleChannel_XandY(const Operator &X, const Operator &Y, arma::mat &X2_CC_ph, arma::mat &Y2_CC_ph, int ch_cc);
    void AddInversePandyaTransformation(const std::deque<arma::mat> &Zbar, Operator &Z); // Changed from the above declaration. Not sure how this was compiling...
    void AddInversePandyaTransformation_SingleChannel(Operator &Z, arma::mat &Zbar, int ch_cc);

    void comm110ss(const Operator &X, const Operator &Y, Operator &Z);
    void comm220ss(const Operator &X, const Operator &Y, Operator &Z);
    void comm111ss(const Operator &X, const Operator &Y, Operator &Z);
    void comm121ss(const Operator &X, const Operator &Y, Operator &Z);
    void comm221ss(const Operator &X, const Operator &Y, Operator &Z);
    void comm122ss(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_pp_hhss(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_phss(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_pp_hh_221ss(const Operator &X, const Operator &Y, Operator &Z);

    void comm122ss_slower(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_phss_slower(const Operator &X, const Operator &Y, Operator &Z);

    void ConstructScalarMpp_Mhh(const Operator &X, const Operator &Y, const Operator &Z, TwoBodyME &Mpp, TwoBodyME &Mhh);

    bool check_2b_channel_Tz_parity(const Operator &Op, Orbit &o1, Orbit &o2, Orbit &o3, Orbit &o4);

    void prod110ss(const Operator &X, const Operator &Y, Operator &Z);
    void prod111ss(const Operator &X, const Operator &Y, Operator &Z);
    void prod112ss(const Operator &X, const Operator &Y, Operator &Z);


}

#endif
