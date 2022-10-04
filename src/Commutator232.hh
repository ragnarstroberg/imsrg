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

#ifndef COMMUTATOR232_H_
#define COMMUTATOR232_H_

#include <cstdlib>
#include <vector>

#include "ModelSpace.hh"
#include "Operator.hh"

namespace comm232 {

void comm232ss_expand_impl(const Operator &X, const Operator &Y,
                      Operator &Z); // implemented and tested.

namespace internal {

int GetJJ1Max(const Operator &Y);
std::vector<std::size_t>
Get2BChannelsValidIn3BChannel(int jj1max, const ThreeBodyChannel &ch_3b,
                              const Operator &Y);

class Basis3BLookupSingleChannel {
public:
  Basis3BLookupSingleChannel(std::size_t i_ch_3b, std::size_t i_ch_2b,
                             const Operator &Y);

  std::size_t Get3BChannelIndex() const { return i_ch_3b_; }
  std::size_t Get2BChannelIndex() const { return i_ch_2b_; }

  std::size_t GetNum2BStates() const { return states_2b_pq.size(); }
  std::size_t Get2BStateP(std::size_t state_2b_i) const {
    return states_2b_p[state_2b_i];
  }
  std::size_t Get2BStateQ(std::size_t state_2b_i) const {
    return states_2b_q[state_2b_i];
  }

  std::size_t GetNum3BStates() const { return states_3b_pqr.size(); }
  std::size_t Get3BStateP(std::size_t state_3b_i) const {
    return states_3b_p[state_3b_i];
  }
  std::size_t Get3BStateQ(std::size_t state_3b_i) const {
    return states_3b_q[state_3b_i];
  }
  std::size_t Get3BStateR(std::size_t state_3b_i) const {
    return states_3b_r[state_3b_i];
  }

  const std::vector<std::size_t> &
  GetRecouplingIndices(std::size_t state_3b_i) const {
    return states_3b_indices[state_3b_i];
  }
  const std::vector<double> &
  GetRecouplingFactors(std::size_t state_3b_i) const {
    return states_3b_recoupling[state_3b_i];
  }

  const std::vector<std::size_t> Get3rdStateBasis() const { return basis_r; }

  std::size_t GetStateLocalIndex(std::size_t p, std::size_t q,
                                 std::size_t r) const {
    return pqr_lookup[p * wrap_factor * wrap_factor + q * wrap_factor + r];
  }

  int IsStateValid(std::size_t p, std::size_t q, std::size_t r) const {
    return pqr_lookup_validity[p * wrap_factor * wrap_factor + q * wrap_factor +
                               r];
  }

private:
  std::size_t i_ch_3b_ = 0;
  std::size_t i_ch_2b_ = 0;
  std::size_t wrap_factor = 0;
  std::vector<std::size_t> states_2b_pq;
  std::vector<std::size_t> states_2b_p;
  std::vector<std::size_t> states_2b_q;
  std::vector<std::size_t> basis_r;
  std::vector<std::size_t> states_3b_pqr;
  std::vector<std::size_t> states_3b_p;
  std::vector<std::size_t> states_3b_q;
  std::vector<std::size_t> states_3b_r;
  std::vector<std::vector<std::size_t>> states_3b_indices;
  std::vector<std::vector<double>> states_3b_recoupling;
  std::vector<std::size_t> pqr_lookup =
      std::vector<std::size_t>(wrap_factor * wrap_factor * wrap_factor, 0);
  std::vector<int> pqr_lookup_validity =
      std::vector<int>(wrap_factor * wrap_factor * wrap_factor, 0);
};

void UnpackMatrix(std::size_t i_ch_3b,
                  const Basis3BLookupSingleChannel &lookup_bra,
                  const Basis3BLookupSingleChannel &lookup_ket,
                  const Operator &op, arma::mat &mat);

void Comm232Diagram1(std::size_t i_ch_2b_bra, std::size_t i_ch_2b_ket,
                     const ThreeBodyChannel &ch_3b,
                     const Basis3BLookupSingleChannel &lookup_bra,
                     const Basis3BLookupSingleChannel &lookup_ket,
                     const Operator &X, const arma::mat &Y_mat,
                     double overall_factor, Operator &Z);

void Comm232Diagram2(std::size_t i_ch_2b_bra, std::size_t i_ch_2b_ket,
                     const ThreeBodyChannel &ch_3b,
                     const Basis3BLookupSingleChannel &lookup_bra,
                     const Basis3BLookupSingleChannel &lookup_ket,
                     const Operator &X, const arma::mat &Y_mat,
                     double overall_factor, Operator &Z);

void Comm232Diagram3(std::size_t i_ch_2b_bra, std::size_t i_ch_2b_ket,
                     const ThreeBodyChannel &ch_3b,
                     const Basis3BLookupSingleChannel &lookup_bra,
                     const Basis3BLookupSingleChannel &lookup_ket,
                     const Operator &X, const arma::mat &Y_mat,
                     double overall_factor, Operator &Z);

void Comm232Diagram4(std::size_t i_ch_2b_bra, std::size_t i_ch_2b_ket,
                     const ThreeBodyChannel &ch_3b,
                     const Basis3BLookupSingleChannel &lookup_bra,
                     const Basis3BLookupSingleChannel &lookup_ket,
                     const Operator &X, const arma::mat &Y_mat,
                     double overall_factor, Operator &Z);

inline double GetSPStateJ(std::size_t orb_i, const Operator &Y) {
  return Y.modelspace->GetOrbit(orb_i).j2 * 0.5;
}

inline int GetPhase(std::size_t orb_i, std::size_t orb_j, int J_ij,
                    const Operator &Y) {
  int jj_i = Y.modelspace->GetOrbit(orb_i).j2;
  int jj_j = Y.modelspace->GetOrbit(orb_j).j2;
  return Y.modelspace->phase((jj_i + jj_j) / 2 - J_ij);
}

inline double GetSPStateOcc(std::size_t orb_i, const Operator &Y) {
  return Y.modelspace->GetOrbit(orb_i).occ;
}

inline double GetSPStateOccBar(std::size_t orb_i, const Operator &Y) {
  return 1.0 - Y.modelspace->GetOrbit(orb_i).occ;
}

} // namespace internal
} // namespace comm232

#endif // COMMUTATOR232_H_
