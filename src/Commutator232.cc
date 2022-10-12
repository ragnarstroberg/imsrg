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

#include "Commutator232.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "ModelSpace.hh"

#include <omp.h>

template <typename T> void Print(std::string prefix, const T &val) {
  std::cout << prefix << ": " << val << "\n";
}

template <typename T>
void Print(std::string prefix, const T &val, const T &val2) {
  std::cout << prefix << ": " << val << "," << val2 << "\n";
}

struct ChannelPair {
  std::size_t i_ch_3b = 0;
  std::size_t i_ch_2b = 0;

  bool operator==(const ChannelPair &other) const {
    return (i_ch_3b == other.i_ch_3b) && (i_ch_2b == other.i_ch_2b);
  }
};

namespace std {
template <> struct hash<ChannelPair> {
  std::size_t operator()(const ChannelPair &x) const {
    return ((hash<std::size_t>()(x.i_ch_3b) ^
             (hash<std::size_t>()(x.i_ch_2b) << 1)) >>
            1);
  }
};
} // namespace std

namespace comm232 {

void comm232ss_expand_impl_new(const Operator &X, const Operator &Y,
                               Operator &Z) {
  std::cout << "In comm232ss_expand_new\n";
  double tstart = omp_get_wtime();
  Z.modelspace->PreCalculateSixJ();

  int hX = 1;
  if (X.IsAntiHermitian())
    hX = -1;
  int hY = 1;
  if (Y.IsAntiHermitian())
    hY = -1;

  const int emax_3body = Z.modelspace->GetEMax3Body();
  const int e3max = Z.modelspace->GetE3max();
  const int jj1max = internal::ExtractJJ1Max(Z);

  std::size_t num_chans = 0;
  std::size_t num_bytes_3b_basis = 0;

  for (std::size_t i_ch_3b = 0;
       i_ch_3b < Y.modelspace->GetNumberThreeBodyChannels(); i_ch_3b += 1) {
    const ThreeBodyChannel &ch_3b = Z.modelspace->GetThreeBodyChannel(i_ch_3b);
    const std::vector<std::size_t> chans_2b =
        internal::Extract2BChannelsValidIn3BChannel(jj1max, i_ch_3b, Z);

    const auto bases_store = internal::PrestoreBases(i_ch_3b, jj1max, Z, e3max);
    for (const auto &basis : bases_store) {
      num_bytes_3b_basis += basis.second.BasisPQR().NumBytes();
    }

    std::size_t num_2b_blocks = 0;
    std::vector<std::pair<std::size_t, std::size_t>> block_ch_2b_indices;
    std::vector<internal::OneBodyBasis> block_beta_bases;

    for (const auto &basis_ijc_full : bases_store) {
      const std::size_t i_ch_2b_ij = basis_ijc_full.first;
      const TwoBodyChannel &ch_2b_ij =
          Z.modelspace->GetTwoBodyChannel(i_ch_2b_ij);

      for (const auto &basis_abalpha_full : bases_store) {
        const std::size_t i_ch_2b_ab = basis_abalpha_full.first;
        const TwoBodyChannel &ch_2b_ab =
            Z.modelspace->GetTwoBodyChannel(i_ch_2b_ab);

        // 3rd external index alpha constrained by being in state | (ab) J_ab
        // alpha >
        const int tz2_alpha = ch_3b.twoTz - 2 * ch_2b_ab.Tz;
        const int parity_alpha = (ch_3b.parity + ch_2b_ab.parity) % 2;
        const int jj_min_alpha = std::abs(ch_3b.twoJ - ch_2b_ab.J * 2);
        const int jj_max_alpha = std::min(ch_3b.twoJ + ch_2b_ab.J * 2, jj1max);

        // Remaining external index beta constrained by being in state | (alpha
        // beta) J_ij >
        const int tz2_beta = ch_2b_ij.Tz * 2 - tz2_alpha;
        const int parity_beta = (ch_2b_ij.parity + parity_alpha) % 2;
        const int jj_min_beta = std::max(ch_2b_ij.J * 2 - jj_max_alpha, 0);
        const int jj_max_beta = std::min(ch_2b_ij.J * 2 + jj_max_alpha, jj1max);

        // Only basis not externally prestored
        internal::OneBodyBasis basis_beta =
            internal::OneBodyBasis::FromQuantumNumbers(
                Z, jj_min_beta, jj_max_beta, parity_beta, tz2_beta);

        if (basis_beta.BasisSize() == 0) {
          continue;
        }

        block_ch_2b_indices.push_back(std::make_pair(i_ch_2b_ij, i_ch_2b_ab));
        block_beta_bases.push_back(std::move(basis_beta));
      }
    }

    std::vector<std::vector<double>> Z_mats(block_ch_2b_indices.size());
#pragma omp parallel for schedule(static)
    for (std::size_t block_index = 0; block_index < block_ch_2b_indices.size();
         block_index += 1) {
      const std::size_t i_ch_2b_ij = block_ch_2b_indices[block_index].first;
      const std::size_t dim_ij =
          bases_store.at(i_ch_2b_ij).BasisPQ().BasisSize();

      Z_mats[block_index] = std::vector<double>(dim_ij * dim_ij, 0.0);
    }

#pragma omp parallel for schedule(guided, 1)
    for (std::size_t block_index = 0; block_index < block_ch_2b_indices.size();
         block_index += 1) {
      auto &Z_mat = Z_mats[block_index];
      const std::size_t i_ch_2b_ij = block_ch_2b_indices[block_index].first;
      const auto &bases_ijc = bases_store.at(i_ch_2b_ij);
      const TwoBodyChannel &ch_2b_ij =
          Z.modelspace->GetTwoBodyChannel(i_ch_2b_ij);

      const internal::TwoBodyBasis &basis_ij = bases_ijc.BasisPQ();
      const internal::TwoBodyBasis basis_ij_e3max = bases_ijc.BasisPQE3Max();
      const internal::OneBodyBasis &basis_c = bases_ijc.BasisR();
      const internal::ThreeBodyBasis &basis_ijc = bases_ijc.BasisPQR();

      const std::size_t i_ch_2b_ab = block_ch_2b_indices[block_index].second;
      const auto &bases_abalpha = bases_store.at(i_ch_2b_ab);
      const TwoBodyChannel &ch_2b_ab =
          Z.modelspace->GetTwoBodyChannel(i_ch_2b_ab);

      double factor = (sqrt(2 * ch_2b_ab.J + 1) / sqrt(2 * ch_2b_ij.J + 1)) *
                      (ch_3b.twoJ + 1);

      const internal::TwoBodyBasis &basis_ab = bases_abalpha.BasisPQ();
      const internal::TwoBodyBasis &basis_ab_e3max =
          bases_abalpha.BasisPQE3Max();
      const internal::OneBodyBasis &basis_alpha = bases_abalpha.BasisR();
      const internal::ThreeBodyBasis &basis_abalpha = bases_abalpha.BasisPQR();

      // Only basis not externally prestored
      const internal::OneBodyBasis basis_beta = block_beta_bases[block_index];

      std::vector<double> six_js_ij = internal::GenerateSixJMatrixIJ(
          Z, basis_ij, basis_ab_e3max, basis_c, ch_2b_ij.J * 2, ch_3b.twoJ,
          ch_2b_ab.J * 2);
      std::vector<double> six_js_ji = internal::GenerateSixJMatrixJI(
          Z, basis_ij, basis_ab_e3max, basis_c, ch_2b_ij.J * 2, ch_3b.twoJ,
          ch_2b_ab.J * 2);
      std::vector<double> phases =
          internal::GeneratePhases(Z, basis_ij, ch_2b_ij.J * 2);
      std::vector<double> occs =
          internal::GenerateOccsMatrix(Z, basis_ab_e3max, basis_c);

      // This block evaluates [X^(3), Y^(2)] = -1 [Y^(2), X^(3)].
      {
        int comm_factor = -1;
        std::vector<double> X_mat_3b = internal::Generate3BMatrix(
            X, i_ch_3b, basis_ijc, basis_abalpha, basis_ij_e3max, basis_alpha,
            basis_ab_e3max, basis_c);
        std::vector<double> Y_mat_2b = internal::Generate2BMatrix(
            Y, i_ch_2b_ab, basis_ab_e3max, basis_beta, basis_c);

        internal::EvaluateComm232Diagram1(
            comm_factor * hX * hY * factor, i_ch_2b_ij, basis_ab_e3max,
            basis_ij_e3max, basis_ij, basis_alpha, basis_beta, basis_c,
            X_mat_3b, Y_mat_2b, occs, six_js_ij, Z_mat);
        internal::EvaluateComm232Diagram2(
            comm_factor * hX * hY * factor, i_ch_2b_ij, basis_ab_e3max,
            basis_ij_e3max, basis_ij, basis_alpha, basis_beta, basis_c,
            X_mat_3b, Y_mat_2b, occs, six_js_ji, phases, Z_mat);
        internal::EvaluateComm232Diagram3(
            comm_factor * factor * -1, i_ch_2b_ij, basis_ab_e3max,
            basis_ij_e3max, basis_ij, basis_alpha, basis_beta, basis_c,
            X_mat_3b, Y_mat_2b, occs, six_js_ij, Z_mat);
        internal::EvaluateComm232Diagram4(
            comm_factor * factor * -1, i_ch_2b_ij, basis_ab_e3max,
            basis_ij_e3max, basis_ij, basis_alpha, basis_beta, basis_c,
            X_mat_3b, Y_mat_2b, occs, six_js_ji, phases, Z_mat);
      }

      // This block evaluates [X^(2), Y^(3)].
      {
        int comm_factor = 1;
        std::vector<double> Y_mat_3b = internal::Generate3BMatrix(
            Y, i_ch_3b, basis_ijc, basis_abalpha, basis_ij_e3max, basis_alpha,
            basis_ab_e3max, basis_c);
        std::vector<double> X_mat_2b = internal::Generate2BMatrix(
            X, i_ch_2b_ab, basis_ab_e3max, basis_beta, basis_c);

        internal::EvaluateComm232Diagram1(
            comm_factor * hX * hY * factor, i_ch_2b_ij, basis_ab_e3max,
            basis_ij_e3max, basis_ij, basis_alpha, basis_beta, basis_c,
            Y_mat_3b, X_mat_2b, occs, six_js_ij, Z_mat);
        internal::EvaluateComm232Diagram2(
            comm_factor * hX * hY * factor, i_ch_2b_ij, basis_ab_e3max,
            basis_ij_e3max, basis_ij, basis_alpha, basis_beta, basis_c,
            Y_mat_3b, X_mat_2b, occs, six_js_ji, phases, Z_mat);
        internal::EvaluateComm232Diagram3(
            comm_factor * factor * -1, i_ch_2b_ij, basis_ab_e3max,
            basis_ij_e3max, basis_ij, basis_alpha, basis_beta, basis_c,
            Y_mat_3b, X_mat_2b, occs, six_js_ij, Z_mat);
        internal::EvaluateComm232Diagram4(
            comm_factor * factor * -1, i_ch_2b_ij, basis_ab_e3max,
            basis_ij_e3max, basis_ij, basis_alpha, basis_beta, basis_c,
            Y_mat_3b, X_mat_2b, occs, six_js_ji, phases, Z_mat);
      }
    }

    for (std::size_t block_index = 0; block_index < block_ch_2b_indices.size();
         block_index += 1) {
      const auto &Z_mat = Z_mats[block_index];
      const std::size_t i_ch_2b_ij = block_ch_2b_indices[block_index].first;
      const auto basis_ij = bases_store.at(i_ch_2b_ij).BasisPQ();
      const std::size_t dim_ij = basis_ij.BasisSize();
      const auto &i_k_vals = basis_ij.GetPVals();
      const auto &j_l_vals = basis_ij.GetQVals();

#pragma omp parallel for schedule(static) collapse(2)
      for (std::size_t i_ij = 0; i_ij < dim_ij; i_ij += 1) {
        for (std::size_t i_kl = 0; i_kl < dim_ij; i_kl += 1) {
          const auto i = i_k_vals[i_ij];
          const auto j = j_l_vals[i_ij];
          const auto k = i_k_vals[i_kl];
          const auto l = j_l_vals[i_kl];

          Z.TwoBody.AddToTBMENonHermNonNormalized(
              i_ch_2b_ij, i_ch_2b_ij, i, j, k, l, Z_mat[i_ij * dim_ij + i_kl]);
        }
      }
    }
    num_2b_blocks += block_ch_2b_indices.size();
    // Print("NUM_2B_BLOCKS", num_2b_blocks);
    num_chans += num_2b_blocks;
  }

  // Print("NUM_CHANS", num_chans);
  // Print("NUM_BYTES_3B_BASIS", num_bytes_3b_basis);
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}

void comm232ss_expand_impl(const Operator &X, const Operator &Y, Operator &Z) {
  // using namespace internal;
  std::cout << "In comm232ss_expand\n";
  Z.modelspace->PreCalculateSixJ();

  int hX = 1;
  if (X.IsAntiHermitian())
    hX = -1;
  int hY = 1;
  if (Y.IsAntiHermitian())
    hY = -1;

  int emax_3body = Y.modelspace->GetEMax3Body();
  int e3max = Y.modelspace->GetE3max();
  int jj1max = internal::GetJJ1Max(Y);

  for (std::size_t i_ch_3b = 0;
       i_ch_3b < Y.modelspace->GetNumberThreeBodyChannels(); i_ch_3b += 1) {
    // if (i_ch_3b != 0) {
    //   break;
    // }
    const ThreeBodyChannel &ch_3b = Y.modelspace->GetThreeBodyChannel(i_ch_3b);
    const std::vector<std::size_t> chans_2b =
        internal::Get2BChannelsValidIn3BChannel(jj1max, ch_3b, Y);

    std::vector<internal::Basis3BLookupSingleChannel> lookup_vec;
    for (const std::size_t i_ch_2b : chans_2b) {
      lookup_vec.emplace_back(i_ch_3b, i_ch_2b, Y);
      if (lookup_vec.back().GetNum3BStates() == 0) {
        lookup_vec.pop_back();
      }
    }

    for (const auto &lookup_bra : lookup_vec) {
      const std::size_t i_ch_2b_bra = lookup_bra.Get2BChannelIndex();
      // if (i_ch_2b_bra != 0) continue;
      const TwoBodyChannel ch_2b_bra =
          Y.modelspace->GetTwoBodyChannel(i_ch_2b_bra);
      for (const auto &lookup_ket : lookup_vec) {
        const std::size_t i_ch_2b_ket = lookup_ket.Get2BChannelIndex();
        // if (i_ch_2b_ket != 6) continue;
        // std::cout << "BRA" << i_ch_2b_bra << "\n";
        // std::cout << "KET" << i_ch_2b_ket << "\n";
        const TwoBodyChannel ch_2b_ket =
            Y.modelspace->GetTwoBodyChannel(i_ch_2b_ket);

        arma::mat X_ch(lookup_bra.GetNum3BStates(),
                       lookup_ket.GetNum3BStates());
        arma::mat Y_ch(lookup_bra.GetNum3BStates(),
                       lookup_ket.GetNum3BStates());

        UnpackMatrix(i_ch_3b, lookup_bra, lookup_ket, X, X_ch);
        UnpackMatrix(i_ch_3b, lookup_bra, lookup_ket, Y, Y_ch);

        // TODO:
        // Prestore occupations

        Comm232Diagram1(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket,
                        X, Y_ch, hX * hY * 1.0, Z);
        Comm232Diagram1(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket,
                        Y, X_ch, hX * hY * -1.0, Z);

        Comm232Diagram2(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket,
                        X, Y_ch, hX * hY * 1.0, Z);
        Comm232Diagram2(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket,
                        Y, X_ch, hX * hY * -1.0, Z);

        Comm232Diagram3(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket,
                        X, Y_ch, 1.0, Z);
        Comm232Diagram3(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket,
                        Y, X_ch, -1.0, Z);

        Comm232Diagram4(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket,
                        X, Y_ch, 1.0, Z);
        Comm232Diagram4(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket,
                        Y, X_ch, -1.0, Z);

        // Print("ZNORM", Z.TwoBodyNorm());
      }
    }
  }

  // comm232ss_expand_impl_new(X, Y, Z);
}

namespace internal {

std::size_t ExtractWrapFactor(const Operator &Z) {
  const auto &sp_indices = Z.modelspace->orbits_3body_space_;
  std::size_t wrap_factor = 0;
  for (const auto &p : sp_indices) {
    wrap_factor = std::max(wrap_factor, static_cast<std::size_t>(p));
  }
  return wrap_factor + 1;
}

int ExtractJJ1Max(const Operator &Z) {
  int jj1max = 1;
  for (const std::size_t &p : Z.modelspace->orbits_3body_space_) {
    const Orbit op = Z.modelspace->GetOrbit(p);
    jj1max = std::max(op.j2, jj1max);
  }
  return jj1max;
}

std::vector<std::size_t> Extract2BChannelsValidIn3BChannel(int jj1max,
                                                           std::size_t i_ch_3b,
                                                           const Operator &Z) {
  const ThreeBodyChannel &ch_3b = Z.modelspace->GetThreeBodyChannel(i_ch_3b);
  std::vector<std::size_t> valid_ch_2;
  for (std::size_t i_ch_2b = 0;
       i_ch_2b < Z.modelspace->GetNumberTwoBodyChannels(); i_ch_2b += 1) {
    const TwoBodyChannel &ch_2b = Z.modelspace->GetTwoBodyChannel(i_ch_2b);

    if ((ch_2b.J * 2 - jj1max <= ch_3b.twoJ) &&
        (ch_2b.J * 2 + jj1max >= ch_3b.twoJ)) {
      if (std::abs(ch_3b.twoTz - ch_2b.Tz * 2) == 1) {
        valid_ch_2.push_back((i_ch_2b));
      }
    }
  }

  return valid_ch_2;
}

std::vector<std::size_t>
GetLookupIndices(const std::vector<std::size_t> &states,
                 std::size_t lookup_size) {
  std::vector<std::size_t> indices(lookup_size, 0);
  for (std::size_t i_p = 0; i_p < states.size(); i_p += 1) {
    const std::size_t p = states[i_p];
    indices[p] = i_p;
  }
  return indices;
}

std::vector<int> GetLookupValidities(const std::vector<std::size_t> &states,
                                     std::size_t lookup_size) {
  std::vector<int> validities(lookup_size, 0);
  for (std::size_t i_p = 0; i_p < states.size(); i_p += 1) {
    const std::size_t p = states[i_p];
    validities[p] = 1;
  }
  return validities;
}

std::vector<std::size_t> Get1BPIndices(const std::vector<std::size_t> &p_states,
                                       std::size_t wrap_factor) {
  return GetLookupIndices(p_states, wrap_factor);
}

std::vector<int> Get1BPValidities(const std::vector<std::size_t> &p_states,
                                  std::size_t wrap_factor) {
  return GetLookupValidities(p_states, wrap_factor);
}

OneBodyBasis OneBodyBasis::FromQuantumNumbers(const Operator &Z, int j2min,
                                              int j2max, int parity, int tz2) {
  const std::size_t wrap_factor = ExtractWrapFactor(Z);

  std::vector<std::size_t> states_1b(Z.modelspace->orbits_3body_space_.begin(),
                                     Z.modelspace->orbits_3body_space_.end());
  std::sort(states_1b.begin(), states_1b.end());

  std::vector<std::size_t> p_states;
  for (const auto &p : states_1b) {
    const Orbit &op = Z.modelspace->GetOrbit(p);
    if ((op.tz2 == tz2) && (op.l % 2 == parity) && (op.j2 <= j2max) &&
        (op.j2 >= j2min)) {
      p_states.push_back(p);
    }
  }

  return OneBodyBasis(p_states, wrap_factor);
}

std::vector<std::size_t> Get2BPStates(const std::vector<std::size_t> &pq_states,
                                      std::size_t wrap_factor) {
  std::vector<std::size_t> p_states(pq_states.size(), 0);
  std::transform(
      pq_states.begin(), pq_states.end(), p_states.begin(),
      [&wrap_factor](const std::size_t &pq) { return pq / wrap_factor; });
  return p_states;
}

std::vector<std::size_t> Get2BQStates(const std::vector<std::size_t> &pq_states,
                                      std::size_t wrap_factor) {
  std::vector<std::size_t> q_states(pq_states.size(), 0);
  std::transform(
      pq_states.begin(), pq_states.end(), q_states.begin(),
      [&wrap_factor](const std::size_t &pq) { return pq % wrap_factor; });
  return q_states;
}

std::vector<std::size_t>
Get2BPQIndices(const std::vector<std::size_t> &pq_states,
               std::size_t wrap_factor) {
  return GetLookupIndices(pq_states, wrap_factor * wrap_factor);
}

std::vector<int> Get2BPQValidities(const std::vector<std::size_t> &pq_states,
                                   std::size_t wrap_factor) {
  return GetLookupValidities(pq_states, wrap_factor * wrap_factor);
}

TwoBodyBasis TwoBodyBasis::PQInTwoBodyChannel(std::size_t i_ch_2b,
                                              const Operator &Z) {
  const std::size_t wrap_factor = ExtractWrapFactor(Z);
  const TwoBodyChannel &ch_2b = Z.modelspace->GetTwoBodyChannel(i_ch_2b);

  std::vector<std::size_t> states_1b(Z.modelspace->orbits_3body_space_.begin(),
                                     Z.modelspace->orbits_3body_space_.end());
  std::sort(states_1b.begin(), states_1b.end());

  std::vector<std::size_t> pq_states;
  for (const auto &p : states_1b) {
    const Orbit &op = Z.modelspace->GetOrbit(p);
    for (const auto &q : states_1b) {
      const Orbit &oq = Z.modelspace->GetOrbit(q);
      if ((p <= q) && ((p != q) || (ch_2b.J % 2 == 0)) && // Pauli principle
          (op.tz2 + oq.tz2 == ch_2b.Tz * 2) &&
          ((op.l + oq.l) % 2 == ch_2b.parity) &&
          (std::abs(op.j2 - oq.j2) <= ch_2b.J * 2) &&
          (std::abs(op.j2 + oq.j2) >= ch_2b.J * 2)) {
        pq_states.push_back(p * wrap_factor + q);
      }
    }
  }

  return TwoBodyBasis(pq_states, wrap_factor);
}

TwoBodyBasis TwoBodyBasis::PQInTwoBodyChannelWithE3Max(std::size_t i_ch_2b,
                                                       const Operator &Z,
                                                       int e3max) {
  const std::size_t wrap_factor = ExtractWrapFactor(Z);
  const TwoBodyChannel &ch_2b = Z.modelspace->GetTwoBodyChannel(i_ch_2b);

  std::vector<std::size_t> states_1b(Z.modelspace->orbits_3body_space_.begin(),
                                     Z.modelspace->orbits_3body_space_.end());
  std::sort(states_1b.begin(), states_1b.end());

  std::vector<std::size_t> pq_states;
  for (const auto &p : states_1b) {
    const Orbit &op = Z.modelspace->GetOrbit(p);
    const int ep = op.n * 2 + op.l;
    for (const auto &q : states_1b) {
      const Orbit &oq = Z.modelspace->GetOrbit(q);
      const int eq = oq.n * 2 + oq.l;
      if ((p <= q) && ((p != q) || (ch_2b.J % 2 == 0)) && // Pauli principle
          (op.tz2 + oq.tz2 == ch_2b.Tz * 2) &&
          ((op.l + oq.l) % 2 == ch_2b.parity) &&
          (std::abs(op.j2 - oq.j2) <= ch_2b.J * 2) &&
          (std::abs(op.j2 + oq.j2) >= ch_2b.J * 2) && (ep + eq <= e3max)) {
        pq_states.push_back(p * wrap_factor + q);
      }
    }
  }

  return TwoBodyBasis(pq_states, wrap_factor);
}

std::vector<std::size_t>
Get3BPStates(const std::vector<std::size_t> &pqr_states,
             std::size_t wrap_factor) {
  std::vector<std::size_t> p_states(pqr_states.size(), 0);
  std::transform(pqr_states.begin(), pqr_states.end(), p_states.begin(),
                 [&wrap_factor](const std::size_t &pqr) {
                   return (pqr / wrap_factor) / wrap_factor;
                 });
  return p_states;
}

std::vector<std::size_t>
Get3BQStates(const std::vector<std::size_t> &pqr_states,
             std::size_t wrap_factor) {
  std::vector<std::size_t> q_states(pqr_states.size(), 0);
  std::transform(pqr_states.begin(), pqr_states.end(), q_states.begin(),
                 [&wrap_factor](const std::size_t &pqr) {
                   return (pqr / wrap_factor) % wrap_factor;
                 });
  return q_states;
}

std::vector<std::size_t>
Get3BRStates(const std::vector<std::size_t> &pqr_states,
             std::size_t wrap_factor) {
  std::vector<std::size_t> r_states(pqr_states.size(), 0);
  std::transform(
      pqr_states.begin(), pqr_states.end(), r_states.begin(),
      [&wrap_factor](const std::size_t &pqr) { return pqr % wrap_factor; });
  return r_states;
}

std::vector<std::size_t>
Get3BPQRIndices(const std::vector<std::size_t> &pqr_states,
                std::size_t wrap_factor) {
  return GetLookupIndices(pqr_states, wrap_factor * wrap_factor * wrap_factor);
}

std::vector<int> Get3BPQRValidities(const std::vector<std::size_t> &pqr_states,
                                    std::size_t wrap_factor) {
  return GetLookupValidities(pqr_states,
                             wrap_factor * wrap_factor * wrap_factor);
}

void GetRecoupling(std::size_t i_ch_3b, std::size_t i_ch_2b, const Operator &Z,
                   const std::vector<std::size_t> &states_3b_p,
                   const std::vector<std::size_t> &states_3b_q,
                   const std::vector<std::size_t> &states_3b_r,
                   std::vector<std::vector<std::size_t>> &indices,
                   std::vector<std::vector<double>> &recoupling) {
  int twoJ = Z.modelspace->GetThreeBodyChannel(i_ch_3b).twoJ;
  int Jab = Z.modelspace->GetTwoBodyChannel(i_ch_2b).J;

  for (std::size_t i = 0; i < states_3b_p.size(); i += 1) {
    Z.ThreeBody.GetKetIndex_withRecoupling(Jab, twoJ, states_3b_p[i],
                                           states_3b_q[i], states_3b_r[i],
                                           indices[i], recoupling[i]);
  }
}

ThreeBodyBasis ThreeBodyBasis::From2BAnd1BBasis(
    std::size_t i_ch_3b, std::size_t i_ch_2b, const Operator &Z,
    const TwoBodyBasis &basis_pq, const OneBodyBasis &basis_r, int e3max) {
  const std::size_t wrap_factor = ExtractWrapFactor(Z);

  std::size_t num_states = 0;
  for (const std::size_t &pq : basis_pq.GetPQVals()) {
    const std::size_t p = pq / wrap_factor;
    const std::size_t q = pq % wrap_factor;
    const Orbit &op = Z.modelspace->GetOrbit(p);
    const Orbit &oq = Z.modelspace->GetOrbit(q);
    const int ep = op.n * 2 + op.l;
    const int eq = oq.n * 2 + oq.l;
    for (const std::size_t &r : basis_r.GetPVals()) {
      const Orbit &oR = Z.modelspace->GetOrbit(r);
      const int er = oR.n * 2 + oR.l;

      if (ep + eq + er <= e3max) {
        num_states += 1;
      }
    }
  }

  std::vector<std::size_t> pqr_states;
  pqr_states.reserve(num_states);
  for (const std::size_t &pq : basis_pq.GetPQVals()) {
    const std::size_t p = pq / wrap_factor;
    const std::size_t q = pq % wrap_factor;
    const Orbit &op = Z.modelspace->GetOrbit(p);
    const Orbit &oq = Z.modelspace->GetOrbit(q);
    const int ep = op.n * 2 + op.l;
    const int eq = oq.n * 2 + oq.l;
    for (const std::size_t &r : basis_r.GetPVals()) {
      const Orbit &oR = Z.modelspace->GetOrbit(r);
      const int er = oR.n * 2 + oR.l;

      if (ep + eq + er <= e3max) {
        pqr_states.push_back(pq * wrap_factor + r);
      }
    }
  }

  return ThreeBodyBasis(i_ch_3b, i_ch_2b, Z, pqr_states, wrap_factor);
}

std::size_t ThreeBodyBasis::NumBytes() const {
  std::size_t size =
      pqr_states.size() * sizeof(std::size_t) +
      p_states.size() * sizeof(std::size_t) +
      q_states.size() * sizeof(std::size_t) +
      r_states.size() * sizeof(std::size_t) +
      pqr_indices.size() * sizeof(std::size_t) +
      pqr_validities.size() * sizeof(int) +
      pqr_me_indices.size() * sizeof(std::vector<std::size_t>) +
      pqr_me_recoupling_factors.size() * sizeof(std::vector<double>);
  for (const auto &indices : pqr_me_indices) {
    size += indices.size() * sizeof(std::size_t);
  }
  for (const auto &recouplings : pqr_me_recoupling_factors) {
    size += recouplings.size() * sizeof(double);
  }
  return size;
}

std::vector<double> Generate3BMatrix(const Operator &Z, std::size_t i_ch_3b,
                                     const ThreeBodyBasis &basis_ijc,
                                     const ThreeBodyBasis &basis_abalpha,
                                     const TwoBodyBasis &basis_ij,
                                     const OneBodyBasis &basis_alpha,
                                     const TwoBodyBasis &basis_ab,
                                     const OneBodyBasis &basis_c) {
  const std::size_t dim_ij = basis_ij.BasisSize();
  const std::size_t dim_c = basis_c.BasisSize();
  const std::size_t dim_ab = basis_ab.BasisSize();
  const std::size_t dim_alpha = basis_alpha.BasisSize();

  std::vector<double> mat_3b(dim_ij * dim_alpha * dim_ab * dim_c, 0.0);

  const std::size_t dim_ijc = basis_ijc.BasisSize();
  const auto &i_vals = basis_ijc.GetPVals();
  const auto &j_vals = basis_ijc.GetQVals();
  const auto &c_vals = basis_ijc.GetRVals();
  const std::size_t dim_abalpha = basis_abalpha.BasisSize();
  const auto &a_vals = basis_abalpha.GetPVals();
  const auto &b_vals = basis_abalpha.GetQVals();
  const auto &alpha_vals = basis_abalpha.GetRVals();

  for (std::size_t i_ijc = 0; i_ijc < dim_ijc; i_ijc += 1) {
    const std::size_t i = i_vals[i_ijc];
    const std::size_t j = j_vals[i_ijc];
    const std::size_t c = c_vals[i_ijc];
    const std::size_t i_ij = basis_ij.GetLocalIndexForPQ(i, j);
    // Print("i", i);
    // Print("j", j);
    // Print("valid?", basis_ij.GetLocalValidityForPQ(i, j));
    const std::size_t i_c = basis_c.GetLocalIndexForP(c);
    // Print("c", c);
    // Print("valid?", basis_c.GetLocalValidityForP(c));
    for (std::size_t i_abalpha = 0; i_abalpha < dim_abalpha; i_abalpha += 1) {
      const std::size_t a = a_vals[i_abalpha];
      const std::size_t b = b_vals[i_abalpha];
      const std::size_t alpha = alpha_vals[i_abalpha];
      const std::size_t i_ab = basis_ab.GetLocalIndexForPQ(a, b);
      // Print("a", a);
      // Print("b", b);
      // Print("valid?", basis_ab.GetLocalValidityForPQ(a, b));
      const std::size_t i_alpha = basis_alpha.GetLocalIndexForP(alpha);
      // Print("alpha", alpha);
      // Print("valid?", basis_alpha.GetLocalValidityForP(alpha));

      const std::size_t i_ijalphabc =
          Index4B(i_ij, i_alpha, i_ab, i_c, dim_alpha, dim_ab, dim_c);
      double me = 0.0;

      const auto &indices_bra = basis_ijc.GetRecouplingIndices(i_ijc);
      const auto &recoupling_bra = basis_ijc.GetRecouplingFactors(i_ijc);
      const auto &indices_ket = basis_abalpha.GetRecouplingIndices(i_abalpha);
      const auto &recoupling_ket =
          basis_abalpha.GetRecouplingFactors(i_abalpha);
      for (std::size_t i_bra_rec = 0; i_bra_rec < indices_bra.size();
           i_bra_rec += 1) {
        const auto bra_index = indices_bra[i_bra_rec];
        const auto bra_factor = recoupling_bra[i_bra_rec];
        for (std::size_t i_ket_rec = 0; i_ket_rec < indices_ket.size();
             i_ket_rec += 1) {
          const auto ket_index = indices_ket[i_ket_rec];
          const auto ket_factor = recoupling_ket[i_ket_rec];
          me +=
              Z.ThreeBody.GetME_pn_ch(i_ch_3b, i_ch_3b, bra_index, ket_index) *
              bra_factor * ket_factor;
        }
      }

      mat_3b[i_ijalphabc] = me;
    }
  }

  return mat_3b;
}

std::vector<double> Generate2BMatrix(const Operator &Z, std::size_t i_ch_2b_ab,
                                     const TwoBodyBasis &basis_ab,
                                     const OneBodyBasis &basis_beta,
                                     const OneBodyBasis &basis_c) {
  const std::size_t dim_c = basis_c.BasisSize();
  const std::size_t dim_ab = basis_ab.BasisSize();
  const std::size_t dim_beta = basis_beta.BasisSize();

  std::vector<double> mat_2b(dim_beta * dim_ab * dim_c, 0.0);

  const auto &c_vals = basis_c.GetPVals();
  const auto &beta_vals = basis_beta.GetPVals();
  const auto &a_vals = basis_ab.GetPVals();
  const auto &b_vals = basis_ab.GetQVals();

  for (std::size_t i_beta = 0; i_beta < dim_beta; i_beta += 1) {
    const std::size_t beta = beta_vals[i_beta];
    for (std::size_t i_ab = 0; i_ab < dim_ab; i_ab += 1) {
      const std::size_t a = a_vals[i_ab];
      const std::size_t b = b_vals[i_ab];
      for (std::size_t i_c = 0; i_c < dim_c; i_c += 1) {
        const std::size_t c = c_vals[i_c];
        const std::size_t i_betaabc = Index3B(i_beta, i_ab, i_c, dim_ab, dim_c);

        mat_2b[i_betaabc] =
            Z.TwoBody.GetTBME(i_ch_2b_ab, i_ch_2b_ab, a, b, c, beta);
      }
    }
  }

  return mat_2b;
}

std::vector<double> GenerateOccsMatrix(const Operator &Z,
                                       const TwoBodyBasis &basis_ab,
                                       const OneBodyBasis &basis_c) {
  const std::size_t dim_c = basis_c.BasisSize();
  const std::size_t dim_ab = basis_ab.BasisSize();

  std::vector<double> mat(dim_ab * dim_c, 0.0);

  const auto &c_vals = basis_c.GetPVals();
  const auto &a_vals = basis_ab.GetPVals();
  const auto &b_vals = basis_ab.GetQVals();

  for (std::size_t i_ab = 0; i_ab < dim_ab; i_ab += 1) {
    const std::size_t a = a_vals[i_ab];
    const std::size_t b = b_vals[i_ab];

    double norm_factor = 1.0;
    if (a == b) {
      norm_factor = 0.5;
    }
    double n_a = Z.modelspace->GetOrbit(a).occ;
    double n_b = Z.modelspace->GetOrbit(b).occ;
    for (std::size_t i_c = 0; i_c < dim_c; i_c += 1) {
      const std::size_t c = c_vals[i_c];
      double n_c = Z.modelspace->GetOrbit(c).occ;
      const std::size_t i_abc = Index2B(i_ab, i_c, dim_c);

      mat[i_abc] =
          norm_factor * ((n_a * n_b) * (1 - n_c) + (1 - n_a) * (1 - n_b) * n_c);
    }
  }

  return mat;
}

std::vector<double> GenerateSixJMatrixIJ(const Operator &Z,
                                         const TwoBodyBasis &basis_ij,
                                         const TwoBodyBasis &basis_ab,
                                         const OneBodyBasis &basis_c, int JJ_ij,
                                         int JJ_3, int JJ_ab) {
  const std::size_t dim_c = basis_c.BasisSize();
  const std::size_t dim_ab = basis_ab.BasisSize();
  const std::size_t dim_ij = basis_ij.BasisSize();

  std::vector<double> mat(dim_ij * dim_ab * dim_c, 0.0);

  const auto &i_vals = basis_ij.GetPVals();
  const auto &j_vals = basis_ij.GetQVals();
  const auto &c_vals = basis_c.GetPVals();

  for (std::size_t i_ij = 0; i_ij < dim_ij; i_ij += 1) {
    const std::size_t i = i_vals[i_ij];
    const std::size_t j = j_vals[i_ij];
    const int jj_i = Z.modelspace->GetOrbit(i).j2;
    const int jj_j = Z.modelspace->GetOrbit(j).j2;

    std::vector<double> local_sixjs(dim_c, 0.0);
    for (std::size_t i_c = 0; i_c < dim_c; i_c += 1) {
      const std::size_t c = c_vals[i_c];
      const int jj_c = Z.modelspace->GetOrbit(c).j2;

      const double sixj =
          Z.modelspace->GetSixJ(jj_i / 2.0, jj_j / 2.0, JJ_ij / 2.0, jj_c / 2.0,
                                JJ_3 / 2.0, JJ_ab / 2.0);

      local_sixjs[i_c] = sixj;
    }

    for (std::size_t i_ab = 0; i_ab < dim_ab; i_ab += 1) {
      for (std::size_t i_c = 0; i_c < dim_c; i_c += 1) {

        const std::size_t i_ijabc = Index3B(i_ij, i_ab, i_c, dim_ab, dim_c);
        mat[i_ijabc] = local_sixjs[i_c];
      }
    }
  }

  return mat;
}

std::vector<double> GenerateSixJMatrixJI(const Operator &Z,
                                         const TwoBodyBasis &basis_ij,
                                         const TwoBodyBasis &basis_ab,
                                         const OneBodyBasis &basis_c, int JJ_ij,
                                         int JJ_3, int JJ_ab) {
  const std::size_t dim_c = basis_c.BasisSize();
  const std::size_t dim_ab = basis_ab.BasisSize();
  const std::size_t dim_ij = basis_ij.BasisSize();

  std::vector<double> mat(dim_ij * dim_ab * dim_c, 0.0);

  const auto &i_vals = basis_ij.GetPVals();
  const auto &j_vals = basis_ij.GetQVals();
  const auto &c_vals = basis_c.GetPVals();

  for (std::size_t i_ij = 0; i_ij < dim_ij; i_ij += 1) {
    const std::size_t i = i_vals[i_ij];
    const std::size_t j = j_vals[i_ij];
    const int jj_i = Z.modelspace->GetOrbit(i).j2;
    const int jj_j = Z.modelspace->GetOrbit(j).j2;

    std::vector<double> local_sixjs(dim_c, 0.0);
    for (std::size_t i_c = 0; i_c < dim_c; i_c += 1) {
      const std::size_t c = c_vals[i_c];
      const int jj_c = Z.modelspace->GetOrbit(c).j2;

      const double sixj =
          Z.modelspace->GetSixJ(jj_j / 2.0, jj_i / 2.0, JJ_ij / 2.0, jj_c / 2.0,
                                JJ_3 / 2.0, JJ_ab / 2.0);

      local_sixjs[i_c] = sixj;
    }

    for (std::size_t i_ab = 0; i_ab < dim_ab; i_ab += 1) {
      for (std::size_t i_c = 0; i_c < dim_c; i_c += 1) {

        const std::size_t i_ijabc = Index3B(i_ij, i_ab, i_c, dim_ab, dim_c);
        mat[i_ijabc] = local_sixjs[i_c];
      }
    }
  }

  return mat;
}

std::vector<double> GeneratePhases(const Operator &Z,
                                   const TwoBodyBasis &basis_ij, int JJ_ij) {

  const std::size_t dim_ij = basis_ij.BasisSize();

  std::vector<double> mat(dim_ij, 0.0);

  const auto &i_vals = basis_ij.GetPVals();
  const auto &j_vals = basis_ij.GetQVals();

  for (std::size_t i_ij = 0; i_ij < dim_ij; i_ij += 1) {
    const std::size_t i = i_vals[i_ij];
    const std::size_t j = j_vals[i_ij];

    const int phase_factor = -1 * GetPhase(i, j, JJ_ij / 2, Z);
    mat[i_ij] = phase_factor;
  }

  return mat;
}

void EvaluateComm232Diagram1(
    double factor, std::size_t i_ch_2b_ij, const TwoBodyBasis &basis_ab_e3max,
    const TwoBodyBasis &basis_ij_e3max, const TwoBodyBasis &basis_ij,
    const OneBodyBasis &basis_alpha, const OneBodyBasis &basis_beta,
    const OneBodyBasis &basis_c, const std::vector<double> &mat_3b,
    const std::vector<double> &mat_2b, const std::vector<double> &occs,
    const std::vector<double> &six_js_ij, std::vector<double> &Z_mat) {
  const auto dim_abc = basis_ab_e3max.BasisSize() * basis_c.BasisSize();
  const auto dim_kl_e3 = basis_ij_e3max.BasisSize();
  const auto dim_ij = basis_ij.BasisSize();
  const auto dim_alpha_abc = basis_alpha.BasisSize() * dim_abc;

  // i == alpha
  const auto &i_vals = basis_ij.GetPVals();
  // j == beta
  const auto &j_vals = basis_ij.GetQVals();
  const auto &k_vals = basis_ij_e3max.GetPVals();
  const auto &l_vals = basis_ij_e3max.GetQVals();

  const double *occs_slice = occs.data();

  for (std::size_t i_kl = 0; i_kl < dim_kl_e3; i_kl += 1) {
    const auto k = k_vals[i_kl];
    const auto l = l_vals[i_kl];
    std::size_t i_kl_full = basis_ij.GetLocalIndexForPQ(k, l);
    for (std::size_t i_ij = 0; i_ij < dim_ij; i_ij += 1) {
      const auto i = i_vals[i_ij];
      const auto j = j_vals[i_ij];

      if (!basis_alpha.GetLocalValidityForP(i) ||
          !basis_beta.GetLocalValidityForP(j)) {
        continue;
      }
      const auto i_i = basis_alpha.GetLocalIndexForP(i);
      const auto i_j = basis_beta.GetLocalIndexForP(j);

      const double *sixj_slice = &(six_js_ij.data()[i_ij * dim_abc]);
      const double *mat2_slice = &(mat_2b.data()[i_j * dim_abc]);
      const double *mat3_slice =
          &(mat_3b.data()[i_kl * dim_alpha_abc + i_i * dim_abc]);

      const double me =
          Comm232Core(sixj_slice, occs_slice, mat2_slice, mat3_slice, dim_abc);
      Z_mat[i_ij * dim_ij + i_kl_full] += factor * me;
    }
  }
}

void EvaluateComm232Diagram2(
    double factor, std::size_t i_ch_2b_ij, const TwoBodyBasis &basis_ab_e3max,
    const TwoBodyBasis &basis_ij_e3max, const TwoBodyBasis &basis_ij,
    const OneBodyBasis &basis_alpha, const OneBodyBasis &basis_beta,
    const OneBodyBasis &basis_c, const std::vector<double> &mat_3b,
    const std::vector<double> &mat_2b, const std::vector<double> &occs,
    const std::vector<double> &six_js_ji, const std::vector<double> &phases,
    std::vector<double> &Z_mat) {
  const auto dim_abc = basis_ab_e3max.BasisSize() * basis_c.BasisSize();
  const auto dim_kl_e3 = basis_ij_e3max.BasisSize();
  const auto dim_ij = basis_ij.BasisSize();
  const auto dim_alpha_abc = basis_alpha.BasisSize() * dim_abc;

  // i == beta
  const auto &i_vals = basis_ij.GetPVals();
  // j == alpha
  const auto &j_vals = basis_ij.GetQVals();
  const auto &k_vals = basis_ij_e3max.GetPVals();
  const auto &l_vals = basis_ij_e3max.GetQVals();

  const double *occs_slice = occs.data();

  for (std::size_t i_kl = 0; i_kl < dim_kl_e3; i_kl += 1) {
    const auto k = k_vals[i_kl];
    const auto l = l_vals[i_kl];
    std::size_t i_kl_full = basis_ij.GetLocalIndexForPQ(k, l);
    for (std::size_t i_ij = 0; i_ij < dim_ij; i_ij += 1) {
      const double phase = phases[i_ij];
      const auto i = i_vals[i_ij];
      const auto j = j_vals[i_ij];

      if (!basis_beta.GetLocalValidityForP(i) ||
          !basis_alpha.GetLocalValidityForP(j)) {
        continue;
      }
      const auto i_i = basis_beta.GetLocalIndexForP(i);
      const auto i_j = basis_alpha.GetLocalIndexForP(j);

      const double *sixj_slice = &(six_js_ji.data()[i_ij * dim_abc]);
      const double *mat2_slice = &(mat_2b.data()[i_i * dim_abc]);
      const double *mat3_slice =
          &(mat_3b.data()[i_kl * dim_alpha_abc + i_j * dim_abc]);

      const double me =
          Comm232Core(sixj_slice, occs_slice, mat2_slice, mat3_slice, dim_abc);
      Z_mat[i_ij * dim_ij + i_kl_full] += phase * factor * me;
    }
  }
}

void EvaluateComm232Diagram3(
    double factor, std::size_t i_ch_2b_ij, const TwoBodyBasis &basis_ab_e3max,
    const TwoBodyBasis &basis_ij_e3max, const TwoBodyBasis &basis_ij,
    const OneBodyBasis &basis_alpha, const OneBodyBasis &basis_beta,
    const OneBodyBasis &basis_c, const std::vector<double> &mat_3b,
    const std::vector<double> &mat_2b, const std::vector<double> &occs,
    const std::vector<double> &six_js_ij, std::vector<double> &Z_mat) {
  const auto dim_abc = basis_ab_e3max.BasisSize() * basis_c.BasisSize();
  const auto dim_ij_e3 = basis_ij_e3max.BasisSize();
  const auto dim_kl = basis_ij.BasisSize();
  const auto dim_alpha_abc = basis_alpha.BasisSize() * dim_abc;

  // k == alpha
  const auto &k_vals = basis_ij.GetPVals();
  // l == beta
  const auto &l_vals = basis_ij.GetQVals();
  const auto &i_vals = basis_ij_e3max.GetPVals();
  const auto &j_vals = basis_ij_e3max.GetQVals();

  const double *occs_slice = occs.data();

  for (std::size_t i_ij = 0; i_ij < dim_ij_e3; i_ij += 1) {
    const auto i = i_vals[i_ij];
    const auto j = j_vals[i_ij];
    std::size_t i_ij_full = basis_ij.GetLocalIndexForPQ(i, j);
    for (std::size_t i_kl = 0; i_kl < dim_kl; i_kl += 1) {
      const auto k = k_vals[i_kl];
      const auto l = l_vals[i_kl];

      if (!basis_alpha.GetLocalValidityForP(k) ||
          !basis_beta.GetLocalValidityForP(l)) {
        continue;
      }
      const auto i_k = basis_alpha.GetLocalIndexForP(k);
      const auto i_l = basis_beta.GetLocalIndexForP(l);

      const double *sixj_slice = &(six_js_ij.data()[i_kl * dim_abc]);
      const double *mat2_slice = &(mat_2b.data()[i_l * dim_abc]);
      const double *mat3_slice =
          &(mat_3b.data()[i_ij * dim_alpha_abc + i_k * dim_abc]);

      const double me =
          Comm232Core(sixj_slice, occs_slice, mat2_slice, mat3_slice, dim_abc);
      Z_mat[i_ij_full * dim_kl + i_kl] += factor * me;
    }
  }
}

void EvaluateComm232Diagram4(
    double factor, std::size_t i_ch_2b_ij, const TwoBodyBasis &basis_ab_e3max,
    const TwoBodyBasis &basis_ij_e3max, const TwoBodyBasis &basis_ij,
    const OneBodyBasis &basis_alpha, const OneBodyBasis &basis_beta,
    const OneBodyBasis &basis_c, const std::vector<double> &mat_3b,
    const std::vector<double> &mat_2b, const std::vector<double> &occs,
    const std::vector<double> &six_js_ji, const std::vector<double> &phases,
    std::vector<double> &Z_mat) {
  const auto dim_abc = basis_ab_e3max.BasisSize() * basis_c.BasisSize();
  const auto dim_ij_e3 = basis_ij_e3max.BasisSize();
  const auto dim_kl = basis_ij.BasisSize();
  const auto dim_alpha_abc = basis_alpha.BasisSize() * dim_abc;

  // k == beta
  const auto &k_vals = basis_ij.GetPVals();
  // l == alpha
  const auto &l_vals = basis_ij.GetQVals();
  const auto &i_vals = basis_ij_e3max.GetPVals();
  const auto &j_vals = basis_ij_e3max.GetQVals();

  const double *occs_slice = occs.data();

  for (std::size_t i_ij = 0; i_ij < dim_ij_e3; i_ij += 1) {
    const auto i = i_vals[i_ij];
    const auto j = j_vals[i_ij];
    std::size_t i_ij_full = basis_ij.GetLocalIndexForPQ(i, j);
    for (std::size_t i_kl = 0; i_kl < dim_kl; i_kl += 1) {
      const auto k = k_vals[i_kl];
      const auto l = l_vals[i_kl];
      const double phase = phases[i_kl];

      if (!basis_beta.GetLocalValidityForP(k) ||
          !basis_alpha.GetLocalValidityForP(l)) {
        continue;
      }
      const auto i_k = basis_beta.GetLocalIndexForP(k);
      const auto i_l = basis_alpha.GetLocalIndexForP(l);

      const double *sixj_slice = &(six_js_ji.data()[i_kl * dim_abc]);
      const double *mat2_slice = &(mat_2b.data()[i_k * dim_abc]);
      const double *mat3_slice =
          &(mat_3b.data()[i_ij * dim_alpha_abc + i_l * dim_abc]);

      const double me =
          Comm232Core(sixj_slice, occs_slice, mat2_slice, mat3_slice, dim_abc);
      Z_mat[i_ij_full * dim_kl + i_kl] += phase * factor * me;
    }
  }
}

double Comm232Core(const double *slice_sixj, const double *slice_occs,
                   const double *slice_mat_2b, const double *slice_mat_3b,
                   std::size_t slice_size) {
  double val = 0.0;

  for (std::size_t abc = 0; abc < slice_size; abc += 1) {
    // if (slice_mat_3b[abc] != 0.0)
    //   Print("ME_3B", slice_mat_3b[abc], slice_sixj[abc]);
    val += slice_sixj[abc] * slice_occs[abc] * slice_mat_2b[abc] *
           slice_mat_3b[abc];
  }

  return val;
}

std::unordered_map<std::size_t, CollectedBases>
PrestoreBases(std::size_t i_ch_3b, int jj1max, const Operator &Z, int e3max) {
  std::unordered_map<std::size_t, CollectedBases> bases;

  const ThreeBodyChannel &ch_3b = Z.modelspace->GetThreeBodyChannel(i_ch_3b);
  const std::vector<std::size_t> chans_2b =
      internal::Extract2BChannelsValidIn3BChannel(jj1max, i_ch_3b, Z);

  for (const std::size_t &i_ch_2b_ij : chans_2b) {
    const TwoBodyChannel &ch_2b_ij =
        Z.modelspace->GetTwoBodyChannel(i_ch_2b_ij);

    internal::TwoBodyBasis basis_ij =
        internal::TwoBodyBasis::PQInTwoBodyChannel(i_ch_2b_ij, Z);
    if (basis_ij.BasisSize() == 0)
      continue;
    internal::TwoBodyBasis basis_ij_e3max =
        internal::TwoBodyBasis::PQInTwoBodyChannelWithE3Max(i_ch_2b_ij, Z,
                                                            e3max);
    if (basis_ij_e3max.BasisSize() == 0)
      continue;

    // 3rd contracted index c constrained by being in state | (ij) J_ij c >
    const int tz2_c = ch_3b.twoTz - 2 * ch_2b_ij.Tz;
    const int parity_c = (ch_3b.parity + ch_2b_ij.parity) % 2;
    const int jj_min_c = std::abs(ch_3b.twoJ - ch_2b_ij.J * 2);
    const int jj_max_c = std::min(ch_3b.twoJ + ch_2b_ij.J * 2, jj1max);

    internal::OneBodyBasis basis_c = internal::OneBodyBasis::FromQuantumNumbers(
        Z, jj_min_c, jj_max_c, parity_c, tz2_c);
    if (basis_c.BasisSize() == 0)
      continue;

    internal::ThreeBodyBasis basis_ijc =
        internal::ThreeBodyBasis::From2BAnd1BBasis(
            i_ch_3b, i_ch_2b_ij, Z, basis_ij_e3max, basis_c, e3max);
    if (basis_ijc.BasisSize() == 0)
      continue;

    bases.emplace(std::make_pair(
        i_ch_2b_ij,
        CollectedBases(std::move(basis_ij), std::move(basis_ij_e3max),
                       std::move(basis_c), std::move(basis_ijc))));
  }

  return bases;
}

// Previous methods

void UnpackMatrix(std::size_t i_ch_3b,
                  const Basis3BLookupSingleChannel &lookup_bra,
                  const Basis3BLookupSingleChannel &lookup_ket,
                  const Operator &op, arma::mat &mat) {
  // #pragma omp parallel for collapse(2)
  for (std::size_t i_bra = 0; i_bra < lookup_bra.GetNum3BStates(); i_bra += 1) {
    for (std::size_t i_ket = 0; i_ket < lookup_ket.GetNum3BStates();
         i_ket += 1) {
      double me = 0.0;
      const auto &indices_bra = lookup_bra.GetRecouplingIndices(i_bra);
      const auto &recoupling_bra = lookup_bra.GetRecouplingFactors(i_bra);
      const auto &indices_ket = lookup_ket.GetRecouplingIndices(i_ket);
      const auto &recoupling_ket = lookup_ket.GetRecouplingFactors(i_ket);
      for (std::size_t i_bra_rec = 0; i_bra_rec < indices_bra.size();
           i_bra_rec += 1) {
        const auto bra_index = indices_bra[i_bra_rec];
        const auto bra_factor = recoupling_bra[i_bra_rec];
        for (std::size_t i_ket_rec = 0; i_ket_rec < indices_ket.size();
             i_ket_rec += 1) {
          const auto ket_index = indices_ket[i_ket_rec];
          const auto ket_factor = recoupling_ket[i_ket_rec];
          me +=
              op.ThreeBody.GetME_pn_ch(i_ch_3b, i_ch_3b, bra_index, ket_index) *
              bra_factor * ket_factor;
        }
      }
      mat(i_bra, i_ket) = me;
    }
  }
}

// Term 1:
// Externally fixed:
// - J_ij=ch_2b_bra.J,
// - J_ab=ch_2b_ket.J,
// - JJ_3=ch_3b.twoJ
//
// Z^(J_ij)_ijkl =
// 0.5 * hat(J_ab) / hat(J_ij) * hat(JJ_3)^2                // hat_factor
// sum_abc (n_a n_b (1 - n_c) + (1 - n_a) (1 - n_b) n_c)    // occ_factor
// six_j(                                                   // six_j
//  j_i j_j  J_ij
//  j_c JJ_3 J_ab
// )
// (
//   X^(J_ab)_abcj Y^(J_ij,J_ab,JJ_3)_klcabi
// )
//
// This diagram requires a factor of hX * hY in front.
void Comm232Diagram1(std::size_t i_ch_2b_bra, std::size_t i_ch_2b_ket,
                     const ThreeBodyChannel &ch_3b,
                     const Basis3BLookupSingleChannel &lookup_bra,
                     const Basis3BLookupSingleChannel &lookup_ket,
                     const Operator &X, const arma::mat &Y_mat,
                     double overall_factor, Operator &Z) {
  const TwoBodyChannel &ch_2b_bra =
      X.modelspace->GetTwoBodyChannel(i_ch_2b_bra);
  const TwoBodyChannel &ch_2b_ket =
      X.modelspace->GetTwoBodyChannel(i_ch_2b_ket);
  const double J_ij = ch_2b_bra.J;
  const double J_ab = ch_2b_ket.J;
  const double J_3 = ch_3b.twoJ * 0.5;
  // No factor of 1/2 because we use a >= b.
  const double hat_factor =
      sqrt(2.0 * J_ab + 1) / sqrt(2.0 * J_ij + 1) * (2.0 * J_3 + 1);
  // #pragma omp parallel for collapse(2)
  for (std::size_t ij_index = 0; ij_index < lookup_bra.GetNum2BStates();
       ij_index += 1) {
    for (std::size_t kl_index = 0; kl_index < lookup_bra.GetNum2BStates();
         kl_index += 1) {
      const auto i = lookup_bra.Get2BStateP(ij_index);
      const auto j = lookup_bra.Get2BStateQ(ij_index);
      const auto j_i = GetSPStateJ(i, X);
      const auto j_j = GetSPStateJ(j, X);

      const auto k = lookup_bra.Get2BStateP(kl_index);
      const auto l = lookup_bra.Get2BStateQ(kl_index);

      double me = 0.0;

      for (std::size_t ab_index = 0; ab_index < lookup_ket.GetNum2BStates();
           ab_index += 1) {
        const auto a = lookup_ket.Get2BStateP(ab_index);
        const auto b = lookup_ket.Get2BStateQ(ab_index);
        if (!lookup_ket.IsStateValid(a, b, i))
          continue;
        double norm_factor = 1.0;
        if (a == b)
          norm_factor = 0.5;
        const double n_a_n_b = GetSPStateOcc(a, X) * GetSPStateOcc(b, X);
        const double nbar_a_nbar_b =
            GetSPStateOccBar(a, X) * GetSPStateOccBar(b, X);
        const std::size_t abi_index = lookup_ket.GetStateLocalIndex(a, b, i);

        for (const std::size_t c : lookup_bra.Get3rdStateBasis()) {
          if (!lookup_bra.IsStateValid(k, l, c))
            continue;
          const auto j_c = GetSPStateJ(c, X);
          const double n_c = GetSPStateOcc(c, X);
          const double occ_factor =
              norm_factor * (n_a_n_b * (1 - n_c) + nbar_a_nbar_b * n_c);
          const double six_j =
              Z.modelspace->GetSixJ(j_i, j_j, J_ij, j_c, J_3, J_ab);
          const std::size_t klc_index = lookup_bra.GetStateLocalIndex(k, l, c);

          const double me_X_abcj =
              X.TwoBody.GetTBME(i_ch_2b_ket, i_ch_2b_ket, a, b, c, j);
          const double me_Y_klcabi = Y_mat(klc_index, abi_index);
          // if (me_Y_klcabi != 0.0)
          //   Print("ME_3B_SRC", me_Y_klcabi, six_j);

          me += occ_factor * six_j * me_X_abcj * me_Y_klcabi;
        }
      }
      Z.TwoBody.AddToTBMENonHermNonNormalized(i_ch_2b_bra, i_ch_2b_bra, i, j, k,
                                              l,
                                              overall_factor * hat_factor * me);
    }
  }
}

// Term 2:
// Externally fixed:
// - J_ij=ch_2b_bra.J,
// - J_ab=ch_2b_ket.J,
// - JJ_3=ch_3b.twoJ
//
// Z^(J_ij)_ijkl =
// 0.5 * hat(J_ab) / hat(J_ij) * hat(JJ_3)^2               // hat_factor
// -1 (-1)^(j_i + j_j - J_ij)                             // phase_factor
// sum_abc (n_a n_b (1 - n_c) + (1 - n_a) (1 - n_b) n_c)    // occ_factor
// six_j(                                                   // six_j
//  j_j j_i  J_ij
//  j_c JJ_3 J_ab
// )
// (
//   X^(J_ab)_abci Y^(J_ab,J_ij,JJ_3)_klcabj
// )
//
// This diagram requires a factor of hX * hY in front.
void Comm232Diagram2(std::size_t i_ch_2b_bra, std::size_t i_ch_2b_ket,
                     const ThreeBodyChannel &ch_3b,
                     const Basis3BLookupSingleChannel &lookup_bra,
                     const Basis3BLookupSingleChannel &lookup_ket,
                     const Operator &X, const arma::mat &Y_mat,
                     double overall_factor, Operator &Z) {
  const TwoBodyChannel &ch_2b_bra =
      X.modelspace->GetTwoBodyChannel(i_ch_2b_bra);
  const TwoBodyChannel &ch_2b_ket =
      X.modelspace->GetTwoBodyChannel(i_ch_2b_ket);
  const double J_ij = ch_2b_bra.J * 1.0;
  const int J_ij_int = ch_2b_bra.J;
  const double J_ab = ch_2b_ket.J * 1.0;
  const double J_3 = ch_3b.twoJ / 2.0;
  // No factor of 1/2 because we use a >= b.
  const double hat_factor =
      sqrt(2.0 * J_ab + 1) / sqrt(2.0 * J_ij + 1) * (2.0 * J_3 + 1);
  // #pragma omp parallel for collapse(2)
  for (std::size_t ij_index = 0; ij_index < lookup_bra.GetNum2BStates();
       ij_index += 1) {
    for (std::size_t kl_index = 0; kl_index < lookup_bra.GetNum2BStates();
         kl_index += 1) {
      const auto i = lookup_bra.Get2BStateP(ij_index);
      const auto j = lookup_bra.Get2BStateQ(ij_index);
      const auto j_i = GetSPStateJ(i, X);
      const auto j_j = GetSPStateJ(j, X);
      const int phase_factor = -1 * GetPhase(i, j, J_ij_int, X);

      const auto k = lookup_bra.Get2BStateP(kl_index);
      const auto l = lookup_bra.Get2BStateQ(kl_index);

      double me = 0.0;

      for (std::size_t ab_index = 0; ab_index < lookup_ket.GetNum2BStates();
           ab_index += 1) {
        const auto a = lookup_ket.Get2BStateP(ab_index);
        const auto b = lookup_ket.Get2BStateQ(ab_index);
        if (!lookup_ket.IsStateValid(a, b, j))
          continue;
        double norm_factor = 1.0;
        if (a == b)
          norm_factor = 0.5;
        const double n_a_n_b = GetSPStateOcc(a, X) * GetSPStateOcc(b, X);
        const double nbar_a_nbar_b =
            GetSPStateOccBar(a, X) * GetSPStateOccBar(b, X);
        const std::size_t abj_index = lookup_ket.GetStateLocalIndex(a, b, j);

        for (const std::size_t c : lookup_bra.Get3rdStateBasis()) {
          if (!lookup_bra.IsStateValid(k, l, c))
            continue;
          const auto j_c = GetSPStateJ(c, X);
          const double n_c = GetSPStateOcc(c, X);
          const double occ_factor =
              norm_factor * (n_a_n_b * (1 - n_c) + nbar_a_nbar_b * n_c);
          const double six_j =
              Z.modelspace->GetSixJ(j_j, j_i, J_ij, j_c, J_3, J_ab);
          const std::size_t klc_index = lookup_bra.GetStateLocalIndex(k, l, c);

          const double me_X_abci =
              X.TwoBody.GetTBME(i_ch_2b_ket, i_ch_2b_ket, a, b, c, i);
          const double me_Y_klcabj = Y_mat(klc_index, abj_index);
          // if (me_Y_klcabj != 0.0)
          //   Print("ME_3B_SRC", me_Y_klcabj, six_j);

          me += occ_factor * six_j * me_X_abci * me_Y_klcabj;
        }
      }
      Z.TwoBody.AddToTBMENonHermNonNormalized(
          i_ch_2b_bra, i_ch_2b_bra, i, j, k, l,
          overall_factor * hat_factor * phase_factor * me);
    }
  }
}

// Term 3:
// Externally fixed:
// - J_kl=ch_2b_bra.J,
// - J_ab=ch_2b_ket.J,
// - JJ_3=ch_3b.twoJ
//
// Z^(J_kl)_ijkl =
// -0.5 * hat(J_ab) / hat(J_kl) * hat(JJ_3)^2                // hat_factor
// sum_abc (n_a n_b (1 - n_c) + (1 - n_a) (1 - n_b) n_c)    // occ_factor
// six_j(                                                   // six_j
//  j_k j_l  J_kl
//  j_c JJ_3 J_ab
// )
// (
//   X^(J_ab)_abcl Y^(J_kl,J_ab,JJ_3)_ijcabk
// )
void Comm232Diagram3(std::size_t i_ch_2b_bra, std::size_t i_ch_2b_ket,
                     const ThreeBodyChannel &ch_3b,
                     const Basis3BLookupSingleChannel &lookup_bra,
                     const Basis3BLookupSingleChannel &lookup_ket,
                     const Operator &X, const arma::mat &Y_mat,
                     double overall_factor, Operator &Z) {
  const TwoBodyChannel &ch_2b_bra =
      X.modelspace->GetTwoBodyChannel(i_ch_2b_bra);
  const TwoBodyChannel &ch_2b_ket =
      X.modelspace->GetTwoBodyChannel(i_ch_2b_ket);
  const double J_kl = ch_2b_bra.J * 1.0;
  const double J_ab = ch_2b_ket.J * 1.0;
  const double J_3 = ch_3b.twoJ / 2.0;
  // No factor of 1/2 because we use a >= b.
  const double hat_factor =
      -1 * sqrt(2.0 * J_ab + 1) / sqrt(2.0 * J_kl + 1) * (2.0 * J_3 + 1);
  // #pragma omp parallel for collapse(2)
  for (std::size_t ij_index = 0; ij_index < lookup_bra.GetNum2BStates();
       ij_index += 1) {
    for (std::size_t kl_index = 0; kl_index < lookup_bra.GetNum2BStates();
         kl_index += 1) {
      const auto i = lookup_bra.Get2BStateP(ij_index);
      const auto j = lookup_bra.Get2BStateQ(ij_index);

      const auto k = lookup_bra.Get2BStateP(kl_index);
      const auto l = lookup_bra.Get2BStateQ(kl_index);
      const auto j_k = GetSPStateJ(k, X);
      const auto j_l = GetSPStateJ(l, X);

      double me = 0.0;

      for (std::size_t ab_index = 0; ab_index < lookup_ket.GetNum2BStates();
           ab_index += 1) {
        const auto a = lookup_ket.Get2BStateP(ab_index);
        const auto b = lookup_ket.Get2BStateQ(ab_index);
        if (!lookup_ket.IsStateValid(a, b, k))
          continue;
        double norm_factor = 1.0;
        if (a == b)
          norm_factor = 0.5;
        const double n_a_n_b = GetSPStateOcc(a, X) * GetSPStateOcc(b, X);
        const double nbar_a_nbar_b =
            GetSPStateOccBar(a, X) * GetSPStateOccBar(b, X);
        const std::size_t abk_index = lookup_ket.GetStateLocalIndex(a, b, k);

        for (const std::size_t c : lookup_bra.Get3rdStateBasis()) {
          if (!lookup_bra.IsStateValid(i, j, c))
            continue;
          const auto j_c = GetSPStateJ(c, X);
          const double n_c = GetSPStateOcc(c, X);
          const double occ_factor =
              norm_factor * (n_a_n_b * (1 - n_c) + nbar_a_nbar_b * n_c);
          const double six_j =
              Z.modelspace->GetSixJ(j_k, j_l, J_kl, j_c, J_3, J_ab);
          const std::size_t ijc_index = lookup_bra.GetStateLocalIndex(i, j, c);

          const double me_X_abcl =
              X.TwoBody.GetTBME(i_ch_2b_ket, i_ch_2b_ket, a, b, c, l);
          const double me_Y_ijcabk = Y_mat(ijc_index, abk_index);

          me += occ_factor * six_j * me_X_abcl * me_Y_ijcabk;
        }
      }
      Z.TwoBody.AddToTBMENonHermNonNormalized(i_ch_2b_bra, i_ch_2b_bra, i, j, k,
                                              l,
                                              overall_factor * hat_factor * me);
    }
  }
}

// Term 4:
// Externally fixed:
// - J_kl=ch_2b_bra.J,
// - J_ab=ch_2b_ket.J,
// - JJ_3=ch_3b.twoJ
//
// Z^(J_kl)_ijkl =
// -0.5 * hat(J_ab) / hat(J_kl) * hat(JJ_3)^2                // hat_factor
// -1 * (-1)^(j_k + j_l - J_kl)                             // phase_factor
// sum_abc (n_a n_b (1 - n_c) + (1 - n_a) (1 - n_b) n_c)    // occ_factor
// six_j(                                                   // six_j
//  j_k j_l  J_kl
//  j_c JJ_3 J_ab
// )
// (
//   X^(J_ab)_abck Y^(J_kl,J_ab,JJ_3)_ijcabl
// )
void Comm232Diagram4(std::size_t i_ch_2b_bra, std::size_t i_ch_2b_ket,
                     const ThreeBodyChannel &ch_3b,
                     const Basis3BLookupSingleChannel &lookup_bra,
                     const Basis3BLookupSingleChannel &lookup_ket,
                     const Operator &X, const arma::mat &Y_mat,
                     double overall_factor, Operator &Z) {
  const TwoBodyChannel &ch_2b_bra =
      X.modelspace->GetTwoBodyChannel(i_ch_2b_bra);
  const TwoBodyChannel &ch_2b_ket =
      X.modelspace->GetTwoBodyChannel(i_ch_2b_ket);
  const double J_kl = ch_2b_bra.J * 1.0;
  const int J_kl_int = ch_2b_bra.J;
  const double J_ab = ch_2b_ket.J * 1.0;
  const double J_3 = ch_3b.twoJ / 2.0;
  // No factor of 1/2 because we use a >= b.
  const double hat_factor =
      -1 * sqrt(2.0 * J_ab + 1) / sqrt(2.0 * J_kl + 1) * (2.0 * J_3 + 1);
  // #pragma omp parallel for collapse(2)
  for (std::size_t ij_index = 0; ij_index < lookup_bra.GetNum2BStates();
       ij_index += 1) {
    for (std::size_t kl_index = 0; kl_index < lookup_bra.GetNum2BStates();
         kl_index += 1) {
      const auto i = lookup_bra.Get2BStateP(ij_index);
      const auto j = lookup_bra.Get2BStateQ(ij_index);

      const auto k = lookup_bra.Get2BStateP(kl_index);
      const auto l = lookup_bra.Get2BStateQ(kl_index);
      const auto j_k = GetSPStateJ(k, X);
      const auto j_l = GetSPStateJ(l, X);
      const int phase_factor = -1 * GetPhase(k, l, J_kl_int, X);

      double me = 0.0;

      for (std::size_t ab_index = 0; ab_index < lookup_ket.GetNum2BStates();
           ab_index += 1) {
        const auto a = lookup_ket.Get2BStateP(ab_index);
        const auto b = lookup_ket.Get2BStateQ(ab_index);
        if (!lookup_ket.IsStateValid(a, b, l))
          continue;
        double norm_factor = 1.0;
        if (a == b)
          norm_factor = 0.5;
        const double n_a_n_b = GetSPStateOcc(a, X) * GetSPStateOcc(b, X);
        const double nbar_a_nbar_b =
            GetSPStateOccBar(a, X) * GetSPStateOccBar(b, X);
        const std::size_t abl_index = lookup_ket.GetStateLocalIndex(a, b, l);

        for (const std::size_t c : lookup_bra.Get3rdStateBasis()) {
          if (!lookup_bra.IsStateValid(i, j, c))
            continue;
          const auto j_c = GetSPStateJ(c, X);
          const double n_c = GetSPStateOcc(c, X);
          const double occ_factor =
              norm_factor * (n_a_n_b * (1 - n_c) + nbar_a_nbar_b * n_c);
          const double six_j =
              Z.modelspace->GetSixJ(j_l, j_k, J_kl, j_c, J_3, J_ab);
          const std::size_t ijc_index = lookup_bra.GetStateLocalIndex(i, j, c);

          const double me_X_abck =
              X.TwoBody.GetTBME(i_ch_2b_ket, i_ch_2b_ket, a, b, c, k);
          const double me_Y_ijcabl = Y_mat(ijc_index, abl_index);

          me += occ_factor * six_j * me_X_abck * me_Y_ijcabl;
        }
      }
      Z.TwoBody.AddToTBMENonHermNonNormalized(
          i_ch_2b_bra, i_ch_2b_bra, i, j, k, l,
          overall_factor * phase_factor * hat_factor * me);
    }
  }
}

int GetJJ1Max(const Operator &Y) {
  int jj1max = 1;
  for (const std::size_t &p : Y.modelspace->orbits_3body_space_) {
    const Orbit op = Y.modelspace->GetOrbit(p);
    jj1max = std::max(op.j2, jj1max);
  }
  return jj1max;
}

std::vector<std::size_t>
Get2BChannelsValidIn3BChannel(int jj1max, const ThreeBodyChannel &ch_3b,
                              const Operator &Y) {
  std::vector<std::size_t> valid_ch_2;
  for (std::size_t i_ch_2b = 0;
       i_ch_2b < Y.modelspace->GetNumberTwoBodyChannels(); i_ch_2b += 1) {
    const TwoBodyChannel &ch_2b = Y.modelspace->GetTwoBodyChannel(i_ch_2b);

    if ((ch_2b.J * 2 - jj1max <= ch_3b.twoJ) &&
        (ch_2b.J * 2 + jj1max >= ch_3b.twoJ)) {
      if (std::abs(ch_3b.twoTz - ch_2b.Tz * 2) == 1) {
        valid_ch_2.push_back((i_ch_2b));
      }
    }
  }

  return valid_ch_2;
}

std::size_t GetWrapFactor(const Operator &Y) {
  std::size_t wrap_factor = 0;
  for (const std::size_t p : Y.modelspace->orbits_3body_space_) {
    wrap_factor = std::max(p, wrap_factor);
  }
  return wrap_factor + 1;
}

std::vector<std::size_t> Get2BStatesPQ(std::size_t i_ch_2b, const Operator &Y,
                                       std::size_t wrap_factor) {
  const TwoBodyChannel ch_2b = Y.modelspace->GetTwoBodyChannel(i_ch_2b);
  const int e3max = Y.modelspace->GetE3max();
  std::vector<std::size_t> pqs;
  for (const std::size_t &p : Y.modelspace->orbits_3body_space_) {
    const Orbit &op = Y.modelspace->GetOrbit(p);
    for (const std::size_t &q : Y.modelspace->orbits_3body_space_) {
      const Orbit &oq = Y.modelspace->GetOrbit(q);
      if ((p <= q) && ((p != q) || (ch_2b.J % 2 == 0)) && // Pauli principle
          (op.tz2 + oq.tz2 == ch_2b.Tz * 2) &&
          ((op.l + oq.l) % 2 == ch_2b.parity) &&
          (std::abs(op.j2 - oq.j2) <= ch_2b.J * 2) &&
          (op.j2 + oq.j2 >= ch_2b.J * 2)) {
        pqs.push_back(wrap_factor * p + q);
      }
    }
  }

  return pqs;
}

std::vector<std::size_t> Get2BStatesP(const std::vector<std::size_t> &pqs,
                                      std::size_t wrap_factor) {
  std::vector<std::size_t> ps(pqs.size());
  std::transform(pqs.begin(), pqs.end(), ps.begin(),
                 [&](std::size_t pq) { return pq / wrap_factor; });
  return ps;
}

std::vector<std::size_t> Get2BStatesQ(const std::vector<std::size_t> &pqs,
                                      std::size_t wrap_factor) {
  std::vector<std::size_t> ps(pqs.size());
  std::transform(pqs.begin(), pqs.end(), ps.begin(),
                 [&](std::size_t pq) { return pq % wrap_factor; });
  return ps;
}

std::vector<std::size_t> GetBasisR(std::size_t i_ch_2b, std::size_t i_ch_3b,
                                   const Operator &Y) {
  const TwoBodyChannel &ch_2b = Y.modelspace->GetTwoBodyChannel(i_ch_2b);
  const ThreeBodyChannel &ch_3b = Y.modelspace->GetThreeBodyChannel(i_ch_3b);
  const int e3max = Y.modelspace->GetE3max();
  std::vector<std::size_t> rs;
  for (const std::size_t &r : Y.modelspace->orbits_3body_space_) {
    const Orbit &oR = Y.modelspace->GetOrbit(r);
    if (((ch_2b.parity + oR.l) % 2 == ch_3b.parity) &&
        (ch_2b.Tz * 2 + oR.tz2 == ch_3b.twoTz) &&
        (std::abs(ch_2b.J * 2 - oR.j2) <= ch_3b.twoJ) &&
        (std::abs(ch_2b.J * 2 + oR.j2) >= ch_3b.twoJ) &&
        (oR.n * 2 + oR.l <= e3max)) {
      rs.push_back(r);
    }
  }

  return rs;
}

std::vector<std::size_t>
Get3BStatesPQR(const std::vector<std::size_t> &states_2b_p,
               const std::vector<std::size_t> &states_2b_q,
               const std::vector<std::size_t> &basis_r, const Operator &Y,
               std::size_t wrap_factor) {
  std::vector<std::size_t> pqrs;

  for (std::size_t i_pq = 0; i_pq < states_2b_p.size(); i_pq += 1) {
    const std::size_t p = states_2b_p[i_pq];
    const Orbit &op = Y.modelspace->GetOrbit(p);
    int ep = 2 * op.n + op.l;
    const std::size_t q = states_2b_q[i_pq];
    const Orbit &oq = Y.modelspace->GetOrbit(q);
    int eq = 2 * oq.n + oq.l;
    for (const std::size_t &r : basis_r) {
      const Orbit &oR = Y.modelspace->GetOrbit(r);
      int er = 2 * oR.n + oR.l;
      if (ep + eq + er <= Y.modelspace->GetE3max()) {
        pqrs.push_back(wrap_factor * wrap_factor * p + wrap_factor * q + r);
      }
    }
  }
  return pqrs;
}

std::vector<std::size_t> Get3BStatesP(const std::vector<std::size_t> &pqrs,
                                      std::size_t wrap_factor) {
  std::vector<std::size_t> ps(pqrs.size());
  std::transform(pqrs.begin(), pqrs.end(), ps.begin(), [&](std::size_t pqr) {
    return (pqr / wrap_factor) / wrap_factor;
  });
  return ps;
}

std::vector<std::size_t> Get3BStatesQ(const std::vector<std::size_t> &pqrs,
                                      std::size_t wrap_factor) {
  std::vector<std::size_t> ps(pqrs.size());
  std::transform(pqrs.begin(), pqrs.end(), ps.begin(), [&](std::size_t pqr) {
    return (pqr / wrap_factor) % wrap_factor;
  });
  return ps;
}

std::vector<std::size_t> Get3BStatesR(const std::vector<std::size_t> &pqrs,
                                      std::size_t wrap_factor) {
  std::vector<std::size_t> ps(pqrs.size());
  std::transform(pqrs.begin(), pqrs.end(), ps.begin(),
                 [&](std::size_t pqr) { return pqr % wrap_factor; });
  return ps;
}

std::vector<std::size_t> GetPQRLookup(const std::vector<std::size_t> &pqrs,
                                      std::size_t wrap_factor) {
  std::vector<std::size_t> pqr_indices(wrap_factor * wrap_factor * wrap_factor,
                                       0);
  for (std::size_t i_pqr = 0; i_pqr < pqrs.size(); i_pqr += 1) {
    const std::size_t pqr = pqrs[i_pqr];
    pqr_indices[pqr] = i_pqr;
  }
  return pqr_indices;
}

std::vector<int> GetPQRLookupValidity(const std::vector<std::size_t> &pqrs,
                                      std::size_t wrap_factor) {
  std::vector<int> pqr_indices_valid(wrap_factor * wrap_factor * wrap_factor,
                                     0);
  for (const std::size_t &pqr : pqrs) {
    pqr_indices_valid[pqr] = 1;
  }
  return pqr_indices_valid;
}

Basis3BLookupSingleChannel::Basis3BLookupSingleChannel(std::size_t i_ch_3b,
                                                       std::size_t i_ch_2b,
                                                       const Operator &Y)
    : i_ch_3b_(i_ch_3b), i_ch_2b_(i_ch_2b), wrap_factor(GetWrapFactor(Y)),
      states_2b_pq(Get2BStatesPQ(i_ch_2b, Y, wrap_factor)),
      states_2b_p(Get2BStatesP(states_2b_pq, wrap_factor)),
      states_2b_q(Get2BStatesQ(states_2b_pq, wrap_factor)),
      basis_r(GetBasisR(i_ch_2b, i_ch_3b, Y)),
      states_3b_pqr(
          Get3BStatesPQR(states_2b_p, states_2b_q, basis_r, Y, wrap_factor)),
      states_3b_p(Get3BStatesP(states_3b_pqr, wrap_factor)),
      states_3b_q(Get3BStatesQ(states_3b_pqr, wrap_factor)),
      states_3b_r(Get3BStatesR(states_3b_pqr, wrap_factor)),
      // Initialized in contructor
      states_3b_indices(states_3b_pqr.size()),
      states_3b_recoupling(states_3b_pqr.size()),
      pqr_lookup(GetPQRLookup(states_3b_pqr, wrap_factor)),
      pqr_lookup_validity(GetPQRLookupValidity(states_3b_pqr, wrap_factor)) {

  GetRecoupling(i_ch_3b_, i_ch_2b_, Y, states_3b_p, states_3b_q, states_3b_r,
                states_3b_indices, states_3b_recoupling);
}
} // namespace internal
} // namespace comm232