///////////////////////////////////////////////////////////////////////////////////
//    TensorCommutators.hh, part of  imsrg++
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

#ifndef TensorCommutators_hh
#define TensorCommutators_hh 1

#include "Operator.hh"
#include <map>
#include <array>
#include <armadillo>

namespace Commutator
{

    void DoTensorPandyaTransformation(const Operator &Z, std::map<std::array<index_t, 2>, arma::mat> &);
    void DoTensorPandyaTransformation_SingleChannel(const Operator &Z, arma::mat &X, int ch_bra_cc, int ch_ket_cc);
    void AddInverseTensorPandyaTransformation(Operator &Z, const std::map<std::array<index_t, 2>, arma::mat> &);
    std::deque<arma::mat> InitializePandya(Operator &Z, size_t nch, std::string orientation);

    void comm111st(const Operator &X, const Operator &Y, Operator &Z);
    void comm121st(const Operator &X, const Operator &Y, Operator &Z);
    void comm122st(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_pp_hh_221st(const Operator &X, const Operator &Y, Operator &Z);
    void comm222_phst(const Operator &X, const Operator &Y, Operator &Z);

    // scalar-tensor with a 3b operator
    void comm331st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm223st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm231st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm232st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm133st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm132st(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test

}// namespace Commutator

#endif
