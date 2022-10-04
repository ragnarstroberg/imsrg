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


namespace comm232 {

void comm232ss_expand_impl( const Operator& X, const Operator& Y, Operator& Z )
{
  // using namespace internal;
  std::cout << "In comm232ss_expand\n";
  Z.modelspace->PreCalculateSixJ();

  int emax_3body = Y.modelspace->GetEMax3Body();
  int e3max = Y.modelspace->GetE3max();
  int jj1max = internal::GetJJ1Max(Y);

  for (std::size_t i_ch_3b = 0; i_ch_3b < Y.modelspace->GetNumberThreeBodyChannels(); i_ch_3b += 1) {
    const ThreeBodyChannel& ch_3b = Y.modelspace->GetThreeBodyChannel(i_ch_3b);
    const std::vector<std::size_t> chans_2b = internal::Get2BChannelsValidIn3BChannel(jj1max, ch_3b, Y);

    std::vector<internal::Basis3BLookupSingleChannel> lookup_vec;
    for (const std::size_t i_ch_2b : chans_2b) {
      lookup_vec.emplace_back(i_ch_3b, i_ch_2b, Y);
      if (lookup_vec.back().GetNum3BStates() == 0) {
        lookup_vec.pop_back();
      }
    }

    for (const auto& lookup_bra : lookup_vec) {
      const std::size_t i_ch_2b_bra = lookup_bra.Get2BChannelIndex();
      const TwoBodyChannel ch_2b_bra = Y.modelspace->GetTwoBodyChannel(i_ch_2b_bra);
      for (const auto& lookup_ket : lookup_vec) {
        const std::size_t i_ch_2b_ket = lookup_ket.Get2BChannelIndex();
        const TwoBodyChannel ch_2b_ket = Y.modelspace->GetTwoBodyChannel(i_ch_2b_ket);

        arma::mat X_ch(lookup_bra.GetNum3BStates(), lookup_ket.GetNum3BStates());
        arma::mat Y_ch(lookup_bra.GetNum3BStates(), lookup_ket.GetNum3BStates());

        UnpackMatrix(i_ch_3b, lookup_bra, lookup_ket, X, X_ch);
        UnpackMatrix(i_ch_3b, lookup_bra, lookup_ket, Y, Y_ch);

        // TODO:
        // Prestore occupations

        Comm232Diagram1(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket, X, Y_ch, 1.0, Z);
        Comm232Diagram1(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket, Y, X_ch, -1.0, Z);

        Comm232Diagram2(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket, X, Y_ch, 1.0, Z);
        Comm232Diagram2(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket, Y, X_ch, -1.0, Z);

        Comm232Diagram3(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket, X, Y_ch, 1.0, Z);
        Comm232Diagram3(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket, Y, X_ch, -1.0, Z);

        Comm232Diagram4(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket, X, Y_ch, 1.0, Z);
        Comm232Diagram4(i_ch_2b_bra, i_ch_2b_ket, ch_3b, lookup_bra, lookup_ket, Y, X_ch, -1.0, Z);
      }
    }
  }
}

namespace internal {

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
// - J_ij=ch_2b_ket.J,
// - J_ab=ch_2b_bra.J,
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
//   X^(J_ab)_cjab Y^(J_ab,J_ij,JJ_3)_abiklc
// )
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
  const double J_ij = ch_2b_ket.J;
  const double J_ab = ch_2b_bra.J;
  const double J_3 = ch_3b.twoJ * 0.5;
  // No factor of 1/2 because we use a >= b.
  const double hat_factor =
      sqrt(2.0 * J_ab + 1) / sqrt(2.0 * J_ij + 1) * (2.0 * J_3 + 1);
// #pragma omp parallel for collapse(2)
  for (std::size_t ij_index = 0; ij_index < lookup_ket.GetNum2BStates();
       ij_index += 1) {
    for (std::size_t kl_index = 0; kl_index < lookup_ket.GetNum2BStates();
         kl_index += 1) {
      const auto i = lookup_ket.Get2BStateP(ij_index);
      const auto j = lookup_ket.Get2BStateQ(ij_index);
      const auto j_i = GetSPStateJ(i, X);
      const auto j_j = GetSPStateJ(j, X);

      const auto k = lookup_ket.Get2BStateP(kl_index);
      const auto l = lookup_ket.Get2BStateQ(kl_index);

      double me = 0.0;

      for (std::size_t ab_index = 0; ab_index < lookup_bra.GetNum2BStates();
           ab_index += 1) {
        const auto a = lookup_bra.Get2BStateP(ab_index);
        const auto b = lookup_bra.Get2BStateQ(ab_index);
        if (!lookup_bra.IsStateValid(a, b, i))
          continue;
        double norm_factor = 1.0;
        if (a == b)
          norm_factor = 0.5;
        const double n_a_n_b = GetSPStateOcc(a, X) * GetSPStateOcc(b, X);
        const double nbar_a_nbar_b =
            GetSPStateOccBar(a, X) * GetSPStateOccBar(b, X);
        const std::size_t abi_index = lookup_bra.GetStateLocalIndex(a, b, i);

        for (const std::size_t c : lookup_ket.Get3rdStateBasis()) {
          if (!lookup_ket.IsStateValid(k, l, c))
            continue;
          const auto j_c = GetSPStateJ(c, X);
          const double n_c = GetSPStateOcc(c, X);
          const double occ_factor =
              norm_factor * (n_a_n_b * (1 - n_c) + nbar_a_nbar_b * n_c);
          const double six_j =
              Z.modelspace->GetSixJ(j_i, j_j, J_ij, j_c, J_3, J_ab);
          const std::size_t klc_index = lookup_ket.GetStateLocalIndex(k, l, c);

          const double me_X_cjab =
              X.TwoBody.GetTBME(i_ch_2b_bra, i_ch_2b_bra, c, j, a, b);
          const double me_Y_abiklc = Y_mat(abi_index, klc_index);

          me += occ_factor * six_j * me_X_cjab * me_Y_abiklc;
        }
      }
      Z.TwoBody.AddToTBMENonHermNonNormalized(i_ch_2b_ket, i_ch_2b_ket, i, j, k,
                                              l,
                                              overall_factor * hat_factor * me);
    }
  }
}

// Term 2:
// Externally fixed:
// - J_ij=ch_2b_ket.J,
// - J_ab=ch_2b_bra.J,
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
//   X^(J_ab)_ciab Y^(J_ab,J_ij,JJ_3)_abjklc
// )
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
  const double J_ij = ch_2b_ket.J * 1.0;
  const int J_ij_int = ch_2b_ket.J;
  const double J_ab = ch_2b_bra.J * 1.0;
  const double J_3 = ch_3b.twoJ / 2.0;
  // No factor of 1/2 because we use a >= b.
  const double hat_factor =
      sqrt(2.0 * J_ab + 1) / sqrt(2.0 * J_ij + 1) * (2.0 * J_3 + 1);
// #pragma omp parallel for collapse(2)
  for (std::size_t ij_index = 0; ij_index < lookup_ket.GetNum2BStates();
       ij_index += 1) {
    for (std::size_t kl_index = 0; kl_index < lookup_ket.GetNum2BStates();
         kl_index += 1) {
      const auto i = lookup_ket.Get2BStateP(ij_index);
      const auto j = lookup_ket.Get2BStateQ(ij_index);
      const auto j_i = GetSPStateJ(i, X);
      const auto j_j = GetSPStateJ(j, X);
      const int phase_factor = -1 * GetPhase(i, j, J_ij_int, X);

      const auto k = lookup_ket.Get2BStateP(kl_index);
      const auto l = lookup_ket.Get2BStateQ(kl_index);

      double me = 0.0;

      for (std::size_t ab_index = 0; ab_index < lookup_bra.GetNum2BStates();
           ab_index += 1) {
        const auto a = lookup_bra.Get2BStateP(ab_index);
        const auto b = lookup_bra.Get2BStateQ(ab_index);
        if (!lookup_bra.IsStateValid(a, b, j))
          continue;
        double norm_factor = 1.0;
        if (a == b)
          norm_factor = 0.5;
        const double n_a_n_b = GetSPStateOcc(a, X) * GetSPStateOcc(b, X);
        const double nbar_a_nbar_b =
            GetSPStateOccBar(a, X) * GetSPStateOccBar(b, X);
        const std::size_t abj_index = lookup_bra.GetStateLocalIndex(a, b, j);

        for (const std::size_t c : lookup_ket.Get3rdStateBasis()) {
          if (!lookup_ket.IsStateValid(k, l, c))
            continue;
          const auto j_c = GetSPStateJ(c, X);
          const double n_c = GetSPStateOcc(c, X);
          const double occ_factor =
              norm_factor * (n_a_n_b * (1 - n_c) + nbar_a_nbar_b * n_c);
          const double six_j =
              Z.modelspace->GetSixJ(j_j, j_i, J_ij, j_c, J_3, J_ab);
          const std::size_t klc_index = lookup_ket.GetStateLocalIndex(k, l, c);

          const double me_X_ciab =
              X.TwoBody.GetTBME(i_ch_2b_bra, i_ch_2b_bra, c, i, a, b);
          const double me_Y_abjklc = Y_mat(abj_index, klc_index);

          me += occ_factor * six_j * me_X_ciab * me_Y_abjklc;
        }
      }
      Z.TwoBody.AddToTBMENonHermNonNormalized(
          i_ch_2b_ket, i_ch_2b_ket, i, j, k, l,
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

  int GetJJ1Max(const Operator& Y) {
    int jj1max = 1;
    for (const std::size_t& p : Y.modelspace->orbits_3body_space_) {
      const Orbit op = Y.modelspace->GetOrbit(p);
      jj1max = std::max(op.j2, jj1max);
    }
    return jj1max;
  }

  std::vector<std::size_t> Get2BChannelsValidIn3BChannel(int jj1max, const ThreeBodyChannel& ch_3b, const Operator& Y) {
    std::vector<std::size_t> valid_ch_2;
    for (std::size_t i_ch_2b = 0; i_ch_2b < Y.modelspace->GetNumberTwoBodyChannels(); i_ch_2b += 1) {
      const TwoBodyChannel& ch_2b = Y.modelspace->GetTwoBodyChannel(i_ch_2b);

      if (
        (ch_2b.J * 2 - jj1max <= ch_3b.twoJ) &&
        (ch_2b.J * 2 + jj1max >= ch_3b.twoJ)
      ) {
        if (std::abs(ch_3b.twoTz - ch_2b.Tz * 2) == 1) {
          valid_ch_2.push_back((i_ch_2b));
        }
      }
    }

    return valid_ch_2;
  }

  std::size_t GetWrapFactor(const Operator& Y) {
    std::size_t wrap_factor = 0;
    for (const std::size_t p : Y.modelspace->orbits_3body_space_) {
      wrap_factor = std::max(p, wrap_factor);
    }
    return wrap_factor + 1;
  }

  std::vector<std::size_t> Get2BStatesPQ(std::size_t i_ch_2b, const Operator& Y, std::size_t wrap_factor) {
    const TwoBodyChannel ch_2b = Y.modelspace->GetTwoBodyChannel(i_ch_2b);
    const int e3max = Y.modelspace->GetE3max();
    std::vector<std::size_t> pqs;
    for (const std::size_t& p : Y.modelspace->orbits_3body_space_) {
      const Orbit& op = Y.modelspace->GetOrbit(p);
      for (const std::size_t& q : Y.modelspace->orbits_3body_space_) {
        const Orbit& oq = Y.modelspace->GetOrbit(q);
        if (
          (p <= q) &&
          ((p != q) || (ch_2b.J % 2 == 0)) &&  // Pauli principle
          (op.tz2 + oq.tz2 == ch_2b.Tz * 2) &&
          ((op.l + oq.l) % 2 == ch_2b.parity) &&
          (std::abs(op.j2 - oq.j2) <= ch_2b.J * 2) &&
          (op.j2 + oq.j2 >= ch_2b.J * 2)
        ) {
          pqs.push_back(wrap_factor * p + q);
        } 
      }
    }

    return pqs;
  }

  std::vector<std::size_t> Get2BStatesP(const std::vector<std::size_t>& pqs, std::size_t wrap_factor) {
    std::vector<std::size_t> ps(pqs.size());
    std::transform(pqs.begin(), pqs.end(), ps.begin(), [&](std::size_t pq) { return pq / wrap_factor; });
    return ps;
  }

  std::vector<std::size_t> Get2BStatesQ(const std::vector<std::size_t>& pqs, std::size_t wrap_factor) {
    std::vector<std::size_t> ps(pqs.size());
    std::transform(pqs.begin(), pqs.end(), ps.begin(), [&](std::size_t pq) { return pq % wrap_factor; });
    return ps;
  }

  std::vector<std::size_t> GetBasisR(std::size_t i_ch_2b, std::size_t i_ch_3b, const Operator& Y) {
    const TwoBodyChannel& ch_2b = Y.modelspace->GetTwoBodyChannel(i_ch_2b);
    const ThreeBodyChannel& ch_3b = Y.modelspace->GetThreeBodyChannel(i_ch_3b);
    const int e3max = Y.modelspace->GetE3max();
    std::vector<std::size_t> rs;
    for (const std::size_t& r : Y.modelspace->orbits_3body_space_) {
      const Orbit& oR = Y.modelspace->GetOrbit(r);
      if (
        ((ch_2b.parity + oR.l) % 2 == ch_3b.parity) &&
        (ch_2b.Tz * 2 + oR.tz2 == ch_3b.twoTz) &&
        (std::abs(ch_2b.J * 2 - oR.j2) <= ch_3b.twoJ) &&
        (std::abs(ch_2b.J * 2 + oR.j2) >= ch_3b.twoJ) &&
        (oR.n * 2 + oR.l <= e3max)
      ) {
        rs.push_back(r);
      }
    }

    return rs;
  }

  std::vector<std::size_t> Get3BStatesPQR(const std::vector<std::size_t>& states_2b_p, const std::vector<std::size_t>& states_2b_q, const std::vector<std::size_t>& basis_r, const Operator& Y, std::size_t wrap_factor) {
    std::vector<std::size_t> pqrs;

    for (std::size_t i_pq = 0; i_pq < states_2b_p.size(); i_pq += 1) {
      const std::size_t p = states_2b_p[i_pq];
      const Orbit& op = Y.modelspace->GetOrbit(p);
      int ep = 2 * op.n + op.l;
      const std::size_t q = states_2b_q[i_pq];
      const Orbit& oq = Y.modelspace->GetOrbit(q);
      int eq = 2 * oq.n + oq.l;
      for (const std::size_t& r : basis_r) {
        const Orbit& oR = Y.modelspace->GetOrbit(r);
        int er = 2 * oR.n + oR.l;
        if (ep + eq + er <= Y.modelspace->GetE3max()) {
          pqrs.push_back(wrap_factor * wrap_factor * p + wrap_factor * q + r);
        }
      }
    }
    return pqrs;
  }

  std::vector<std::size_t> Get3BStatesP(const std::vector<std::size_t>& pqrs, std::size_t wrap_factor) {
    std::vector<std::size_t> ps(pqrs.size());
    std::transform(pqrs.begin(), pqrs.end(), ps.begin(), [&](std::size_t pqr) { return (pqr / wrap_factor) / wrap_factor; });
    return ps;
  }

  std::vector<std::size_t> Get3BStatesQ(const std::vector<std::size_t>& pqrs, std::size_t wrap_factor) {
    std::vector<std::size_t> ps(pqrs.size());
    std::transform(pqrs.begin(), pqrs.end(), ps.begin(), [&](std::size_t pqr) { return (pqr / wrap_factor) % wrap_factor; });
    return ps;
  }

  std::vector<std::size_t> Get3BStatesR(const std::vector<std::size_t>& pqrs, std::size_t wrap_factor) {
    std::vector<std::size_t> ps(pqrs.size());
    std::transform(pqrs.begin(), pqrs.end(), ps.begin(), [&](std::size_t pqr) { return pqr % wrap_factor; });
    return ps;
  }

  std::vector<std::size_t> GetPQRLookup(const std::vector<std::size_t>& pqrs, std::size_t wrap_factor) {
    std::vector<std::size_t> pqr_indices(wrap_factor * wrap_factor * wrap_factor, 0);
    for (std::size_t i_pqr = 0; i_pqr < pqrs.size(); i_pqr += 1) {
      const std::size_t pqr = pqrs[i_pqr];
      pqr_indices[pqr] = i_pqr;
    }
    return pqr_indices;
  }

  std::vector<int> GetPQRLookupValidity(const std::vector<std::size_t>& pqrs, std::size_t wrap_factor) {
    std::vector<int> pqr_indices_valid(wrap_factor * wrap_factor * wrap_factor, 0);
    for (const std::size_t& pqr : pqrs) {
      pqr_indices_valid[pqr] = 1;
    }
    return pqr_indices_valid;
  }

  void GetRecoupling(
    std::size_t i_ch_3b, std::size_t i_ch_2b,
    const Operator& Y,
    const std::vector<std::size_t>& states_3b_p,
    const std::vector<std::size_t>& states_3b_q,
    const std::vector<std::size_t>& states_3b_r,
    std::vector<std::vector<std::size_t>>& indices, std::vector<std::vector<double>>& recoupling) {
      int twoJ = Y.modelspace->GetThreeBodyChannel(i_ch_3b).twoJ;
      int Jab = Y.modelspace->GetTwoBodyChannel(i_ch_2b).J;

      for (std::size_t i = 0; i < states_3b_p.size(); i += 1) {
        Y.ThreeBody.GetKetIndex_withRecoupling(Jab, twoJ, states_3b_p[i], states_3b_q[i], states_3b_r[i], indices[i], recoupling[i]);
      }
    }


  Basis3BLookupSingleChannel::Basis3BLookupSingleChannel(std::size_t i_ch_3b, std::size_t i_ch_2b, const Operator& Y) :
      i_ch_3b_(i_ch_3b),
      i_ch_2b_(i_ch_2b),
      wrap_factor(GetWrapFactor(Y)),
      states_2b_pq(Get2BStatesPQ(i_ch_2b, Y, wrap_factor)),
      states_2b_p(Get2BStatesP(states_2b_pq, wrap_factor)),
      states_2b_q(Get2BStatesQ(states_2b_pq, wrap_factor)),
      basis_r(GetBasisR(i_ch_2b, i_ch_3b, Y)),
      states_3b_pqr(Get3BStatesPQR(states_2b_p, states_2b_q, basis_r, Y, wrap_factor)),
      states_3b_p(Get3BStatesP(states_3b_pqr, wrap_factor)),
      states_3b_q(Get3BStatesQ(states_3b_pqr, wrap_factor)),
      states_3b_r(Get3BStatesR(states_3b_pqr, wrap_factor)),
      // Initialized in contructor
      states_3b_indices(states_3b_pqr.size()),
      states_3b_recoupling(states_3b_pqr.size()),
      pqr_lookup(GetPQRLookup(states_3b_pqr, wrap_factor)),
      pqr_lookup_validity(GetPQRLookupValidity(states_3b_pqr, wrap_factor))
      {

        GetRecoupling(i_ch_3b_, i_ch_2b_, Y, states_3b_p, states_3b_q, states_3b_r, states_3b_indices, states_3b_recoupling);
      }
}
}