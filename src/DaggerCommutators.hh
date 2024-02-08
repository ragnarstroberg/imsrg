///////////////////////////////////////////////////////////////////////////////////
//    DaggerCommutators.hh, part of  imsrg++
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

#ifndef DaggerCommutators_hh
#define DaggerCommutators_hh 1

#include "Operator.hh"
#include "ThreeLegME.hh"
#include <armadillo>

namespace Commutator
{

    // commutator terms involving a dagger operator. 211 means [two legs, one leg] => one leg
    // sd means scalar-dagger
    void comm211sd(const Operator &X, const Operator &Y, Operator &Z);
    void comm231sd(const Operator &X, const Operator &Y, Operator &Z);
    void comm431sd(const Operator &X, const Operator &Y, Operator &Z);
    void comm413_233sd(const Operator &X, const Operator &Y, Operator &Z);
    void comm433sd_pphh(const Operator &X, const Operator &Y, Operator &Z);
    void comm433sd_ph(const Operator &X, const Operator &Y, Operator &Z);
    void comm433sd_ph_dumbway(const Operator &X, const Operator &Y, Operator &Z); // Do it with loops, not matmult. Easier to check, but much slower.

    void comm433_pp_hh_431sd(const Operator &X, const Operator &Y, Operator &Z);
    //  void ConstructDaggerMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, TwoBodyME& Mpp, TwoBodyME& Mhh);
    void ConstructDaggerMpp_Mhh(const Operator &X, const Operator &Y, const Operator &Z, ThreeLegME &Mpp, ThreeLegME &Mhh);
    void DoPandyaTransformation_SingleChannel_Dagger(const Operator &Z, arma::mat &X, int ch_cc);
    void AddInversePandyaTransformation_Dagger(const std::deque<arma::mat> &Zbar, Operator &Z);


}// namespace Commutator


#endif
