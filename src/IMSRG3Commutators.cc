
#include "IMSRG3Commutators.hh"
#include "Commutator.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh"


namespace Commutator
{

  bool imsrg3_no_qqq = false;
  bool imsrg3_valence_2b = false;
  bool discard_0b_from_3b = false;
  bool discard_1b_from_3b = false;
  bool discard_2b_from_3b = false;
//  bool Commutator::verbose = false;
  bool perturbative_triples = false;

  double imsrg3_dE6max = 1e20;
  double threebody_threshold = 0;

//  void SetIMSRG3Verbose(bool tf) { Commutator::verbose = tf; };

  void Discard0bFrom3b(bool tf) { discard_0b_from_3b = tf; };
  void Discard1bFrom3b(bool tf) { discard_1b_from_3b = tf; };
  void Discard2bFrom3b(bool tf) { discard_2b_from_3b = tf; };


  void SetThreebodyThreshold(double x)
  {
    threebody_threshold = x;
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  ////////////   BEGIN SCALAR-SCALAR COMMUTATORS WITH 3-body ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  //  For the three-body operators, we often encounter formulas involving the J-coupled
  //  permutation operator P(ij/k)^{J1,J}, which is defined as
  //
  //  P(ij/k)^{J1,J} = 1 - sum_J2 sqrt{(2J1+1)(2J2+1)} (-1)^2(ji+jj+jk)   {ji jk J1} P_ik
  //                                                                      {jk J  J2}
  //
  //                     - sum_J2 sqrt{(2J1+1)(2J2+1)} (-1)^{jj+jk+J1+J2) {jk ji J1} P_jk
  //                                                                      {jk J  J2}
  //

  //*****************************************************************************************
  //
  //    ~~~~~~~~~~~~~~~~        Uncoupled expression:
  //   /\      /\      /\         Z_0 = 1/36 sum_abcdef (nanbnc n`dn`en`f) (X_abcdef Y_defabc - Y_abcdef X_defabc)
  // a(  )d  b(  )e  c(  )f
  //   \/      \/      \/      Coupled expression:
  //    ~~~~~~~~~~~~~~~~        Z_0 = 1/36 sum_abcdef sum_J1,J2,J  (nanbnc n`dn`en`f)  (2J+1)
  //                                    (X^{J1J2J}_abcdef Y^{J1J2J}_defabc - Y^{J1J2J}_abcdef X^{J1J2J}_defabc)
  //    Verified with UnitTest
  //
  /// Three body commutator expression
  /// \f{equation}{
  /// Z_0 = \frac{1}{36} \sum_{abcdef} \sum_{J_{1}J_{2}J} n_a n_b n_c \bar{n}_d\bar{n}_e\bar{n}_f (2J+1)
  ///        ( X^{J_{1}J_{2}J}_{abcdef} Y^{J_1J_2J}_{defabc} - Y^{J_{1}J_{2}J}_{abcdef} X^{J_1J_2J}_{defabc})
  /// \f}
  ///
  void comm330ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    double z0 = 0;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    if (X3.Norm() < 1e-8 or Y3.Norm() < 1e-8)
      return;
    if ((X.GetJRank() != Y.GetJRank()) or (X.GetParity() != Y.GetParity()) or (X.GetTRank() != Y.GetTRank()))
      return;

    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    // roll ch3 and ibra into a single index to improve load balancing
    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : X.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = X.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        //       bra_ket_channels.push_back( { it.first.ch_bra,it.first.ch_ket, static_cast<size_t>(ibra) } ); // (ch_bra, ch_ket,ibra)
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, ibra}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();

#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z0)
    for (size_t ibra_ket = 0; ibra_ket < n_bra_ket_ch; ibra_ket++)
    {

      size_t ch3bra = bra_ket_channels[ibra_ket][0];
      size_t ch3ket = bra_ket_channels[ibra_ket][1];
      size_t ibra = bra_ket_channels[ibra_ket][2];

      //    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      auto &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch3bra);
      auto &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch3ket);
      int twoJ = Tbc_bra.twoJ;
      size_t nkets = Tbc_ket.GetNumberKets();
      //    for (size_t ibra=0; ibra<nkets; ibra++)
      //    {
      Ket3 &bra = Tbc_bra.GetKet(ibra);
      int ea = 2 * bra.op->n + bra.op->l;
      int eb = 2 * bra.oq->n + bra.oq->l;
      int ec = 2 * bra.oR->n + bra.oR->l;
      int tza = bra.op->tz2;
      int tzb = bra.oq->tz2;
      int tzc = bra.oR->tz2;
      double na = bra.op->occ;
      double nb = bra.oq->occ;
      double nc = bra.oR->occ;
      double occnat_a = bra.op->occ_nat;
      double occnat_b = bra.oq->occ_nat;
      double occnat_c = bra.oR->occ_nat;
      double abc_symm = 6;
      if (bra.p == bra.q and bra.q == bra.r)
        abc_symm = 1;
      else if (bra.p == bra.q or bra.q == bra.r or bra.p==bra.r)
        abc_symm = 3;
      if ((std::abs(ea - e_fermi[tza]) + std::abs(eb - e_fermi[tzb]) + std::abs(ec - e_fermi[tzc])) > Z.modelspace->GetdE3max())
        continue;
      if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
        continue;
      size_t iket_max = nkets;
      if (ch3bra == ch3ket)
        iket_max = ibra; // dont need iket=ibra because the commutator will be zero

      for (size_t iket = 0; iket < iket_max; iket++)
      {
        Ket3 &ket = Tbc_ket.GetKet(iket);
        double nd = ket.op->occ;
        double ne = ket.oq->occ;
        double nf = ket.oR->occ;
        double occnat_d = ket.op->occ_nat;
        double occnat_e = ket.oq->occ_nat;
        double occnat_f = ket.oR->occ_nat;
        int ed = 2 * ket.op->n + ket.op->l;
        int ee = 2 * ket.oq->n + ket.oq->l;
        int ef = 2 * ket.oR->n + ket.oR->l;
        int tzd = ket.op->tz2;
        int tze = ket.oq->tz2;
        int tzf = ket.oR->tz2;
        if ((std::abs(ed - e_fermi[tzd]) + std::abs(ee - e_fermi[tze]) + std::abs(ef - e_fermi[tzf])) > Z.modelspace->GetdE3max())
          continue;
        if ((occnat_d * (1 - occnat_d) * occnat_e * (1 - occnat_e) * occnat_f * (1 - occnat_f)) < Z.modelspace->GetOccNat3Cut())
          continue;
        // account for the iket>ibra case, which we don't do explicitly
        double occfactor = na * nb * nc * (1 - nd) * (1 - ne) * (1 - nf) - (1 - na) * (1 - nb) * (1 - nc) * nd * ne * nf;
        if (std::abs(occfactor) < 1e-6)
          continue;
        double def_symm = 6;
        if (ket.p == ket.q and ket.q == ket.r)
          def_symm = 1;
        else if (ket.p == ket.q or ket.q == ket.r or ket.p==ket.r)
          def_symm = 3;

        double xabcdef = X3.GetME_pn_ch(ch3bra, ch3ket, ibra, iket);
        double yabcdef = Y3.GetME_pn_ch(ch3bra, ch3ket, ibra, iket);
        double xdefabc = X3.GetME_pn_ch(ch3ket, ch3bra, iket, ibra);
        double ydefabc = Y3.GetME_pn_ch(ch3ket, ch3bra, iket, ibra);

        z0 += 1. / 36 * occfactor * abc_symm * def_symm * (twoJ + 1) * (xabcdef * ydefabc - yabcdef * xdefabc);

      } // for iket
    }   // for ch3  ibra_ket

    Z.ZeroBody += z0;
    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  //*****************************************************************************************
  //                   |i
  //    *~~~~~~~[X]~~~~*        Uncoupled expression:
  //   / \      / \    |          Z_ij = 1/12 sum_abcde (nanb n`c n`dn`e) (X_abicde Y_cdeabj - Y_abicde X_cdeabj)
  // a(   )c  b(   )d  |e
  //   \ /      \ /    |       Coupled expression:
  //    *~~~~~~~[Y]~~~~*        Z_ij = 1/12 sum_abcde sum_J1,J2,J  (nanb n`c n`dn`e)  (2J+1)/(2ji+1)
  //                   |j               (X^{J1J2J}_abicde Y^{J1J2J}_cdeabj - Y^{J1J2J}_abicde X^{J1J2J}_cdeabj)
  //
  //  Verfied with UnitTest
  //
  void comm331ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z1 = Z.OneBody;
    int herm = Z.IsAntiHermitian() ? -1 : 1;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();
    int emax_3body = Z.modelspace->GetEMax3Body();
    int e3max = Z.modelspace->GetE3max();

    if (X.GetParticleRank() < 3)
      return;
    if (Y.GetParticleRank() < 3)
      return;

    bool x_channel_diag = X.GetParity() == 0 and X.GetTRank() == 0;
    bool y_channel_diag = Y.GetParity() == 0 and Y.GetTRank() == 0;

    size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
    size_t norb = Z.modelspace->GetNumberOrbits();
    int num_threads = omp_get_max_threads();
    std::vector<arma::mat> Z_mats;
    Z_mats.reserve(num_threads);
    for (int i = 0; i < num_threads; i += 1)
    {
      Z_mats.push_back(arma::zeros(norb, norb));
    }
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (size_t i = 0; i < norb; i++)
    {
      for (size_t ch_ab = 0; ch_ab < nch2; ch_ab++)
      {
        int tid = omp_get_thread_num();
        Orbit &oi = Z.modelspace->GetOrbit(i);
        int ei = 2 * oi.n + oi.l;
        if (ei > emax_3body)
          continue;
        double occnat_i = oi.occ_nat;
        int tzi = oi.tz2;
        auto &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
        int Jab = tbc_ab.J;
        size_t nkets_ab = tbc_ab.GetNumberKets();
        int twoJ_min = std::abs(2 * Jab - oi.j2);
        int twoJ_max = 2 * Jab + oi.j2;

        for (auto j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
        {
          if (j > i)
            continue;
          double zij = 0;
          Orbit &oj = Z.modelspace->GetOrbit(j);
          double occnat_j = oj.occ_nat;
          int ej = 2 * oj.n + oj.l;
          if (ej > emax_3body)
            continue;
          int tzj = oj.tz2;

          for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
          {
            Ket &ket_ab = tbc_ab.GetKet(iket_ab);
            size_t a = ket_ab.p;
            size_t b = ket_ab.q;
            //           if (std::abs( ket_ab.op->occ * ket_ab.oq->occ )<1e-6 ) continue;
            double na = ket_ab.op->occ;
            double nb = ket_ab.oq->occ;
            if ((std::abs(na * nb) + std::abs((1 - na) * (1 - nb))) < 1e-8)
              continue;
            int ea = 2 * ket_ab.op->n + ket_ab.op->l;
            int eb = 2 * ket_ab.oq->n + ket_ab.oq->l;
            if (ea > emax_3body)
              continue;
            if (eb > emax_3body)
              continue;
            if (ea + eb + ei > e3max)
              continue;
            if (ea + eb + ej > e3max)
              continue;
            int tza = ket_ab.op->tz2;
            int tzb = ket_ab.oq->tz2;
            double occnat_a = ket_ab.op->occ_nat;
            double occnat_b = ket_ab.oq->occ_nat;
            //      double occnat_c = bra.oR->occ_nat;
            if ((std::abs(ea - e_fermi[tza]) + std::abs(eb - e_fermi[tzb]) + std::abs(ej - e_fermi[tzj])) > Z.modelspace->GetdE3max())
              continue;
            if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) < Z.modelspace->GetOccNat3Cut())
              continue;

            int parity_abi = (tbc_ab.parity + oi.l) % 2;
            int twoTz_abi = tbc_ab.Tz * 2 + oi.tz2;
            int parity_abj = (tbc_ab.parity + oj.l) % 2;
            int twoTz_abj = tbc_ab.Tz * 2 + oj.tz2;

            for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
            {
              double Jfactor = (twoJ + 1.0) / (oi.j2 + 1);
              double ab_symmetry_factor = (a == b) ? 1.0 : 2.0;

              for (int parity_cde : {0, 1})
              {
                for (int twoTz_cde : {-3, -1, 1, 3})
                {
                  bool x_abi_good = (parity_abi + parity_cde) % 2 == X.GetParity() and std::abs(twoTz_abi - twoTz_cde) / 2 == X.GetTRank();
                  bool y_abi_good = (parity_abi + parity_cde) % 2 == Y.GetParity() and std::abs(twoTz_abi - twoTz_cde) / 2 == Y.GetTRank();
                  bool x_abj_good = (parity_abj + parity_cde) % 2 == X.GetParity() and std::abs(twoTz_abj - twoTz_cde) / 2 == X.GetTRank();
                  bool y_abj_good = (parity_abj + parity_cde) % 2 == Y.GetParity() and std::abs(twoTz_abj - twoTz_cde) / 2 == Y.GetTRank();
                  if (not((x_abi_good and y_abj_good) or (x_abj_good and y_abi_good)))
                    continue;
                  size_t ch_cde = Z.modelspace->GetThreeBodyChannelIndex(twoJ, parity_cde, twoTz_cde);
                  if (ch_cde == size_t(-1))
                    continue; // maybe that channel doesn't exist
                  auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch_cde);
                  size_t nkets3 = Tbc.GetNumberKets();

                  for (size_t iket_cde = 0; iket_cde < nkets3; iket_cde++)
                  {
                    Ket3 &ket_cde = Tbc.GetKet(iket_cde);
                    double nc = ket_cde.op->occ;
                    double nd = ket_cde.oq->occ;
                    double ne = ket_cde.oR->occ;
                    double occfactor = na * nb * (1 - nc) * (1 - nd) * (1 - ne) + (1 - na) * (1 - nb) * nc * nd * ne; // fixes mistake found by Matthias Heinz Oct 2022
                    if (std::abs(occfactor) < 1e-8)
                      continue;
                    size_t c = ket_cde.p;
                    size_t d = ket_cde.q;
                    size_t e = ket_cde.r;
                    int ec = 2 * ket_cde.op->n + ket_cde.op->l;
                    int ed = 2 * ket_cde.oq->n + ket_cde.oq->l;
                    int ee = 2 * ket_cde.oR->n + ket_cde.oR->l;
                    if (ec > emax_3body)
                      continue;
                    if (ed > emax_3body)
                      continue;
                    if (ee > emax_3body)
                      continue;
                    if (ec + ed + ee > e3max)
                      continue;
                    double occnat_c = ket_cde.op->occ_nat;
                    double occnat_d = ket_cde.oq->occ_nat;
                    double occnat_e = ket_cde.oR->occ_nat;
                    int tzc = ket_cde.op->tz2;
                    int tzd = ket_cde.oq->tz2;
                    int tze = ket_cde.oR->tz2;
                    if ((std::abs(ec - e_fermi[tzc]) + std::abs(ed - e_fermi[tzd]) + std::abs(ee - e_fermi[tze])) > Z.modelspace->GetdE3max())
                      continue;
                    if ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_e * (1 - occnat_e)) < Z.modelspace->GetOccNat3Cut())
                      continue;
                    double cde_symmetry_factor = 6;
                    if (c == d and d == e)
                      cde_symmetry_factor = 1;
                    else if (c == d or d == e)
                      cde_symmetry_factor = 3;
                    int Jcd = ket_cde.Jpq;

                    double x_abicde = 0, y_abicde = 0, x_cdeabj = 0, y_cdeabj = 0;
                    // If needed we should use GetME_pn_TwoOps here.
                    if (x_channel_diag and y_channel_diag)
                    {
                      auto xy_abicde = X3.GetME_pn_TwoOps(Jab, Jcd, twoJ, a, b, i, c, d, e, X3, Y3);
                      auto xy_cdeabj = X3.GetME_pn_TwoOps(Jcd, Jab, twoJ, c, d, e, a, b, j, X3, Y3);
                      x_abicde = xy_abicde[0];
                      y_abicde = xy_abicde[1];
                      x_cdeabj = xy_cdeabj[0];
                      y_cdeabj = xy_cdeabj[1];
                    }
                    else
                    {
                      if (x_abi_good and y_abj_good)
                      {
                        x_abicde = X3.GetME_pn(Jab, Jcd, twoJ, a, b, i, c, d, e);
                        y_cdeabj = Y3.GetME_pn(Jcd, Jab, twoJ, c, d, e, a, b, j);
                      }
                      if (y_abi_good and x_abj_good)
                      {
                        y_abicde = Y3.GetME_pn(Jab, Jcd, twoJ, a, b, i, c, d, e);
                        x_cdeabj = X3.GetME_pn(Jcd, Jab, twoJ, c, d, e, a, b, j);
                      }
                    }
                    zij += 1. / 12 * ab_symmetry_factor * cde_symmetry_factor * occfactor * Jfactor * (x_abicde * y_cdeabj - y_abicde * x_cdeabj);

                  } // for iket_cde
                }   // for twoTz_cde
              }     // for parity_cde
            }       // for twoJ
          }         // for iket_ab

          Z_mats[tid](i, j) += zij;
          if (i != j)
            Z_mats[tid](j, i) += herm * zij;
        } // for j
      }   // for ch_ab
    }     // for i

    for (int tid = 0; tid < num_threads; tid += 1)
    {
      for (std::size_t i = 0; i < norb; i += 1)
      {
        for (std::size_t j = 0; j < norb; j += 1)
        {
          Z1(i, j) += Z_mats[tid](i, j);
        }
      }
    }

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  //*****************************************************************************************
  //                   |i
  //    *~~[X]~~*      |        Uncoupled expression:
  //   / \      / \    |          Z_ij = 1/4 sum_abcd (nanb n`cn`d) (X_abcd Y_cdiabj - Y_abicdj X_cdab)
  // a(   )c  b(   )d  |
  //   \ /      \ /    |       Coupled expression:
  //   *~~~~~~~[Y]~~~~~*        Z_ij = 1/4sum_abcd sum_J1J  (nanb n`c n`d)  (2J+1)/(2ji+1)
  //                   |j               (X^{J1}_abcd Y^{J1J1J}_cdiabj - Y^{J1J1J}_abicdj X^{J1}_cdab)
  //                   |
  //                              We only sum a<=b and c<=d, so we do not explicitly include the factor 1/4,
  //                              except for the a==b, or c==d case, where there is no double counting.
  //  Verfied with UnitTest
  //
  void comm231ss(const Operator &X, const Operator &Y, Operator &Z)
  {

    double tstart = omp_get_wtime();
    auto &X2 = X.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y2 = Y.TwoBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z1 = Z.OneBody;
    int hZ = Z.IsHermitian() ? +1 : -1;
    int x_particle_rank = X.GetParticleRank();
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    int herm = Z.IsHermitian() ? 1 : -1;

    size_t norb = Z.modelspace->GetNumberOrbits();
    std::vector<std::array<size_t, 2>> ij_pairs;
    for (size_t i = 0; i < norb; i++)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      if (!Z.ThreeBody.IsOrbitIn3BodyEMaxTruncation(oi))
        continue;
      for (auto j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        if (j > i)
          continue;
        if (!Z.ThreeBody.IsOrbitIn3BodyEMaxTruncation(j))
          continue;
        ij_pairs.push_back({i, j});
      }
    }
    size_t nij = ij_pairs.size();

    #pragma omp parallel for schedule(dynamic, 1) 
    for (size_t index_ij = 0; index_ij < nij; index_ij++)
    {
      double t_local = omp_get_wtime();
      size_t i = ij_pairs[index_ij][0];
      size_t j = ij_pairs[index_ij][1];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      int ei = 2 * oi.n + oi.l;
      double d_ei = std::abs(ei - e_fermi[oi.tz2]);
      double occnat_i = oi.occ_nat;
      Orbit &oj = Z.modelspace->GetOrbit(j);
      int ej = 2 * oj.n + oj.l;
      double d_ej = std::abs(ej - e_fermi[oj.tz2]);
      double occnat_j = oj.occ_nat;
      double zij = 0;

      // We do the same thing for X2 Y3 and Y2 X3, so instead of duplicating code,
      // put it inside a loop where xy_iter=0 handles X2 Y3 and xy_iter handles Y2 X3
      for (int xy_iter = 0; xy_iter <= 1; xy_iter++)
      {
        auto &OP2 = xy_iter == 0 ? X : Y;
        auto &OP3 = xy_iter == 0 ? Y : X;
        if (OP3.ThreeBodyNorm()<1e-8 ) continue; // no use if there isn't a 3body piece
        int itersign = xy_iter == 0 ? +1 : -1;

        // Do the loop over 2b channel blocks
        for (auto &it_ch : OP2.TwoBody.MatEl)
        {
          size_t ch_bra = it_ch.first[0];
          size_t ch_ket = it_ch.first[1];
          TwoBodyChannel &tbc_bra = OP2.modelspace->GetTwoBodyChannel(ch_bra);
          TwoBodyChannel &tbc_ket = OP2.modelspace->GetTwoBodyChannel(ch_ket);
          int J = tbc_bra.J;
          size_t nbras = tbc_bra.GetNumberKets();
          size_t nkets = tbc_ket.GetNumberKets();

          for (size_t ibra = 0; ibra < nbras; ibra++)
          {
            Ket &bra = tbc_bra.GetKet(ibra);
            int a = bra.p;
            int b = bra.q;
            if (!OP3.ThreeBody.IsOrbitIn3BodyEMaxTruncation(a))
              continue;
            if (!OP3.ThreeBody.IsOrbitIn3BodyEMaxTruncation(b))
              continue;
            double na = bra.op->occ;
            double nb = bra.oq->occ;
            int ea = 2 * bra.op->n + bra.op->l;
            int eb = 2 * bra.oq->n + bra.oq->l;
            double occnat_a = bra.op->occ_nat;
            double occnat_b = bra.oq->occ_nat;
            double d_ea = std::abs(ea - e_fermi[bra.op->tz2]);
            double d_eb = std::abs(eb - e_fermi[bra.oq->tz2]);
            double occnat_abi = occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_i * (1 - occnat_i);
            double occnat_abj = occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j);
            bool abi_ok = (occnat_abi > Z.modelspace->GetOccNat3Cut()) and ((ea + eb + ei) <= Z.modelspace->E3max) and ((d_ea + d_eb + d_ei) <= Z.modelspace->dE3max);
            bool abj_ok = (occnat_abj > Z.modelspace->GetOccNat3Cut()) and ((ea + eb + ej) <= Z.modelspace->E3max) and ((d_ea + d_eb + d_ej) <= Z.modelspace->dE3max);
            if (not(abi_ok or abj_ok))
              continue;

            // check stuff

            for (size_t iket = 0; iket < nkets; iket++)
            {
              if ((ch_bra == ch_ket) and (iket > ibra))
                continue;
              Ket &ket = tbc_ket.GetKet(iket);
              int c = ket.p;
              int d = ket.q;
              if (!OP3.ThreeBody.IsOrbitIn3BodyEMaxTruncation(c))
                continue;
              if (!OP3.ThreeBody.IsOrbitIn3BodyEMaxTruncation(d))
                continue;
              double nc = ket.op->occ;
              double nd = ket.oq->occ;
              int ec = 2 * ket.op->n + ket.op->l;
              int ed = 2 * ket.oq->n + ket.oq->l;
              double d_ec = std::abs(ec - e_fermi[ket.op->tz2]);
              double d_ed = std::abs(ec - e_fermi[ket.oq->tz2]);
              double occnat_c = bra.op->occ_nat;
              double occnat_d = bra.oq->occ_nat;

              double occnat_cdj = occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_j * (1 - occnat_j);
              double occnat_cdi = occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_i * (1 - occnat_i);

              bool cdi_ok = (occnat_cdi > Z.modelspace->GetOccNat3Cut()) and ((ec + ed + ei) <= Z.modelspace->E3max) and ((d_ec + d_ed + d_ei) <= Z.modelspace->dE3max);
              bool cdj_ok = (occnat_cdj > Z.modelspace->GetOccNat3Cut()) and ((ec + ed + ej) <= Z.modelspace->E3max) and ((d_ec + d_ed + d_ej) <= Z.modelspace->dE3max);

              double prefactor = na * nb * (1 - nc) * (1 - nd) - (1 - na) * (1 - nb) * nc * nd;
              if (std::abs(prefactor) < 1e-8)
                continue;
              double Xabcd = OP2.TwoBody.GetTBME(ch_bra, ch_ket, bra, ket);
              double Xcdab = OP2.TwoBody.GetTBME(ch_ket, ch_bra, ket, bra);

              if (a == b)
                prefactor /= 2;
              if (c == d)
                prefactor /= 2;
              int twoJ_min = std::abs(2 * J - oi.j2);
              int twoJ_max = 2 * J + oi.j2;
              for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
              {
                double yabicdj = 0;
                double ycdiabj = 0;
                if (cdi_ok and abj_ok)
                {
                  ycdiabj = OP3.ThreeBody.GetME_pn(J, J, twoJ, c, d, i, a, b, j);
                }
                if (abi_ok and cdj_ok)
                {
                  yabicdj = OP3.ThreeBody.GetME_pn(J, J, twoJ, a, b, i, c, d, j);
                }
                zij += prefactor * itersign * (twoJ + 1) * ((Xabcd * ycdiabj - yabicdj * Xcdab));
              } // for twoJ

            } // for iket
          }   // for ibra

        } // for it_ch, X 2b channel blocks

      } // for xy_iter

      Z1(i, j) += zij / (oi.j2 + 1.0);
      if (i != j)
      {
        Z1(j, i) += hZ * zij / (oi.j2 + 1.0);
      }
    } // for i

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }


  //*****************************************************************************************
  //
  // i|  j|   *~~[X]  Uncoupled expression:
  //  |   |  / \          Z_ijkl = sum_ab (nan`b-n`anb) (X_ab * Y_ijbkla)
  //  |   | (a  )b
  //  |   |  \ /
  //  *~~[Y]~~*       Coupled expression:
  //  |   |              Z_{ijkl}^{J} = sum_ab (nan`b-n`anb) sum_J' (2J'+1)/(2J+1) ( X_ab * Y_{ijbkla}^{J,J,J'} )
  // k|  l|
  //
  //
  //                 When we call Symmetrize, we copy the upper triangle <ibra|iket> with ibra<=iket onto the lower triangle.
  //                 Of course, calling AddToTBME updates both according to the symmetry, so we don't need to worry so much...
  //
  //  Verfied with UnitTest
  //
  void comm132ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    auto &X1 = X.OneBody;
    auto &X3 = X.ThreeBody;
    auto &Y1 = Y.OneBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    int hZ = Z.IsHermitian() ? +1 : -1;

    bool x_channel_diag = X.GetParity() == 0 and X.GetTRank() == 0;
    bool y_channel_diag = Y.GetParity() == 0 and Y.GetTRank() == 0;

    //  int x_particle_rank = X.GetParticleRank();
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    Z.modelspace->PreCalculateSixJ(); // Presumably this has already been called. If so, this does nothing. But we want to make sure.

    int norb = Z.modelspace->GetNumberOrbits();
    size_t nch = Z.modelspace->GetNumberTwoBodyChannels();
    // The only potentially thread-unsafe part of this loop is the access to 3b matrix elements which might require
    // recoupling, leading to a 6j. If these are precomputed, there is no thread safety issue, so no need to check first_pass
    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
    //  for (size_t ch=0; ch<nch; ch++)
    //  {
    std::vector<size_t> ch_bra_list;
    std::vector<size_t> ch_ket_list;
    for (auto &it_ch : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(it_ch.first[0]);
      ch_ket_list.push_back(it_ch.first[1]);
    }
    size_t nbraket = ch_bra_list.size();
//    for (auto &it_ch : Z.TwoBody.MatEl)
//#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ibraket=0; ibraket<nbraket; ibraket++)
    {
//      size_t ch_bra = it_ch.first[0];
//      size_t ch_ket = it_ch.first[1];
      size_t ch_bra = ch_bra_list[ibraket];
      size_t ch_ket = ch_ket_list[ibraket];
      auto &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      auto &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
      for (int ibra = 0; ibra < nbras; ibra++) // <ij| states
      {
        for (int iket = 0; iket < nkets; iket++) // |kl> states
        {
          if ((ch_bra == ch_ket) and iket > ibra)
            continue; // Use Hermiticity to avoid doing twice the work
          Ket &bra = tbc_bra.GetKet(ibra);
          int i = bra.p;
          int j = bra.q;
          if (!Y3.IsOrbitIn3BodyEMaxTruncation(i))
            continue;
          if (!Y3.IsOrbitIn3BodyEMaxTruncation(j))
            continue;
          int ei = 2 * bra.op->n + bra.op->l;
          int ej = 2 * bra.oq->n + bra.oq->l;
          double d_ei = std::abs(ei - e_fermi[bra.op->tz2]);
          double d_ej = std::abs(ej - e_fermi[bra.oq->tz2]);
          double occnat_i = bra.op->occ_nat;
          double occnat_j = bra.oq->occ_nat;
          Ket &ket = tbc_ket.GetKet(iket);
          int k = ket.p;
          int l = ket.q;
          if (!Y3.IsOrbitIn3BodyEMaxTruncation(k))
            continue;
          if (!Y3.IsOrbitIn3BodyEMaxTruncation(l))
            continue;
          int ek = 2 * ket.op->n + ket.op->l;
          int el = 2 * ket.oq->n + ket.oq->l;
          double d_ek = std::abs(ek - e_fermi[ket.op->tz2]);
          double d_el = std::abs(el - e_fermi[ket.oq->tz2]);
          double occnat_k = ket.op->occ_nat;
          double occnat_l = ket.oq->occ_nat;
          double zijkl = 0;
          // TODO loop over 1-body channels to do the 3-body recoupling in the outer loop
          // Then loop over a and b within that channel.
          for (int a = 0; a < norb; a++)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            if (!Y3.IsOrbitIn3BodyEMaxTruncation(oa))
              continue;

            int ea = 2 * oa.n + oa.l;
            double d_ea = std::abs(ea - e_fermi.at(oa.tz2));
            double occnat_a = oa.occ_nat;
            if ((ek + el + ea) > Z.modelspace->E3max)
              continue;
            if ((d_ek + d_el + d_ea) > Z.modelspace->dE3max)
              continue;
            if ((occnat_k * (1 - occnat_k) * occnat_l * (1 - occnat_l) * occnat_a * (1 - occnat_a)) < Z.modelspace->GetOccNat3Cut())
              continue;

            int twoJ_min = std::abs(oa.j2 - 2 * J);
            int twoJ_max = oa.j2 + 2 * J;
            // this is less efficient than doing the loop twice, but less code
            // and easily lets us treat the case of Y having nonzero Tz or odd parity
            for (int b_loop = 0; b_loop <= 1; b_loop++)
            {
              if (x_channel_diag and y_channel_diag and b_loop > 0)
                continue;
              if ( b_loop ==0 and Y.ThreeBodyNorm()<1e-8)
                 continue;
              if ( b_loop ==1 and X.ThreeBodyNorm()<1e-8)
                 continue;
              std::set<size_t> blist;
              //             std::cout << "HERE AT LINE " << __LINE__ <<  " and a is " << oa.l << " " << oa.j2 << " " << oa.tz2 << std::endl;
              //             for ( auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2})) blist.insert(b);
              if (b_loop == 0)
                for (auto b : X.GetOneBodyChannel(oa.l, oa.j2, oa.tz2))
                  blist.insert(b);
              //             for ( auto b : Y.OneBodyChannels.at({oa.l,oa.j2,oa.tz2})) blist.insert(b);
              if (b_loop == 1)
                for (auto b : Y.GetOneBodyChannel(oa.l, oa.j2, oa.tz2))
                  blist.insert(b);
              //             for ( auto b : Z.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) // TODO: We can make this a<=b or a>=b, I think. Just need to mind some factors of 2
              for (auto b : blist)
              {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                if (!Y3.IsOrbitIn3BodyEMaxTruncation(ob))
                  continue;

                int eb = 2 * ob.n + ob.l;
                double d_eb = std::abs(eb - e_fermi[ob.tz2]);
                double occnat_b = ob.occ_nat;
                if ((ei + ej + eb) > Z.modelspace->E3max)
                  continue;
                if ((d_ei + d_ej + d_eb) > Z.modelspace->dE3max)
                  continue;
                if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_b * (1 - occnat_b)) < Z.modelspace->GetOccNat3Cut())
                  continue;
                double occfactor = oa.occ - ob.occ;
                if (std::abs(occfactor) < 1e-6)
                  continue;
                for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                {
                  //                 double xijbkla = X3.GetME_pn(J,J,twoJ,i,j,b,k,l,a);

                  double xijbkla = 0, yijbkla = 0;

                  if (x_channel_diag and y_channel_diag)
                  {
                    auto xandy = Y3.GetME_pn_TwoOps(J, J, twoJ, i, j, b, k, l, a, X3, Y3);
                    xijbkla = xandy[0];
                    yijbkla = xandy[1];
                  }
                  else
                  {
                    if (b_loop == 0)
                      yijbkla = Y3.GetME_pn(J, J, twoJ, i, j, b, k, l, a);
                    if (b_loop == 1)
                      xijbkla = X3.GetME_pn(J, J, twoJ, i, j, b, k, l, a);
                  }
                  //                 double yijbkla = Y3.GetME_pn(J,J,twoJ,i,j,b,k,l,a);

                  zijkl += occfactor * (twoJ + 1.) / (2 * J + 1) * (X1(a, b) * yijbkla - Y1(a, b) * xijbkla);
                }
              }
            } // for b_loop
          }
          // normalize the tbme
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          //        Z2.AddToTBMENonHerm(ch, ch, ibra, iket, zijkl );
          Z2.AddToTBMENonHerm(ch_bra, ch_ket, ibra, iket, zijkl);
          if ((ch_bra == ch_ket) and ibra != iket)
          {
            Z2.AddToTBMENonHerm(ch_ket, ch_bra, iket, ibra, zijkl * hZ);
          }
        }
      }
    }

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  void comm232ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();

    //  comm232ss_expand_full(X, Y, Z);
    comm232ss_srs_optimized(X, Y, Z);

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  //*****************************************************************************************
  //
  //  |         |
  // i|        j|     Uncoupled expression:
  //  *~~[X]~*  |         Z_ijkl = -1/2 sum_abc (nanbn`c+n`an`bnc) ( (1-Pij) X_icab * Y_abjklc - (1-Pkl) Yijcabl * Xabkc )
  // a|    b/c\ |
  //  |    /   \|
  //  *~~[Y]~~~~*      Coupled expression:
  //  |   |              Z_{ijkl}^{J} = -1/2 sum_abc (nanbn`c+n`an`bnc) sum_J'J" (2J'+1)(2J"+1)/sqrt(2J+1) (-1)^{2J"+J'-J}
  // k|  l|                           *  [   (1 - (-1)^{i+j-J}Pij) (-1)^{j-c} { j  J" J' } X_icab^{J'} * Y_{abjklc}^{J'JJ"}
  //                                                                          { c  i  J  }
  //
  //                                       -(1 - (-1)^{k+l-J}Pkl)  (-1)^{l-c} { l  J" J' } Y_{ijcabl}^{JJ'J"} * X_{abkc}^{J'}  ]
  //                                                                          { c  k  J  }
  //
  //   The factor 1/2 out front is absorbed by the fact that we only do the ordering a<=b  (need to deal carefully with a==b)
  //
  //  For efficiency, we unfold this diagram to look like
  //  |
  // i|
  //  *~~[X]~*
  //  |a  b / \c
  //  |    /   \
//  *~~[Y]~~~*
  //  |   |    |
  //  |k  |l   |j
  //
  //  Verified with UnitTest
  //
  // This is the time hog of the n^7 scaling terms   (seems to be doing better...)
  // Now this is fine. The 223 commutator is taking all the time.
  //

  // This optimized implementation is originally by Ragnar.
  // Strategy: Enumerate states |i> and states |klj`> for the output operator Z
  //            and form a dense matrix <i|Z|klj`> which will be filled via matrix multiplication:
  //             <i|Z|klj`> = <i|X|abc`> <abc`|Y|klj`> - <i|Y|abc`><abc`|X|klj`>
  //           This means we need to identify the appropriate intermediate pph states |abc`>
  //           If X,Y,Z are all channel-diagonal, then |i>, |klj`>, |abc`> all have the same JpTz quantum numbers.
  //           If X or Y have odd parity, then |abc`> will have the same or opposite parity to |i>, depending XY or YX.
  //           If X or Y change Tz, then this is the most complicated case, because multiple Tz could be possible in |abc`>
  //           Generally, the relevant |abc`> will be different for the XY and YX cases, unless X and Y are channel-diagonal.
  //           As a last step, we need to transform <i|Z|klj`> back to <ij|Z|kl>.
  void comm232ss_srs_optimized(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    //  auto& X2 = X.TwoBody;
    auto &X3 = X.ThreeBody;
    //  auto& Y2 = Y.TwoBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    bool x_has_3 = X3.IsAllocated();
    bool y_has_3 = Y3.IsAllocated();
    bool x_channel_diag = X.GetTRank() == 0 and X.GetParity() == 0;
    bool y_channel_diag = Y.GetTRank() == 0 and Y.GetParity() == 0;
    bool z_channel_diag = Z.GetTRank() == 0 and Z.GetParity() == 0;

    // Eventually, this version should be able to handle both cases (channel diagonal and not)
    // without too much of a performance hit. But that's not the case for now. This is a quick
    // and dirty way of doing things.
    if (not(x_channel_diag and y_channel_diag and z_channel_diag))
    {
      comm232ss_slow(X, Y, Z);
      return;
    }
    else
    {
      comm232ss_srs_optimized_old(X, Y, Z);
      return;
    }

    int hermX = X.IsHermitian() ? 1 : -1;
    int hermY = Y.IsHermitian() ? 1 : -1;

    Z.modelspace->PreCalculateSixJ(); // If we already did this, this does nothing.

    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    // Bookkeeping strategy:
    // We need to enumerate all JpTz channels for |i>, and all JpTz channels for |klj`>
    // I believe that for a given operator X,Y,Z, a state |klj`> belongs with a unique channel |i>
    // For Tz changing operators, tzi = tzk+tzl-tzj +- Tz_op.
    // So each sub-block can be indexed by the single-particle channels.
    // However, the channel to which |klj`> belongs depends on the operator in question.
    // Things we need to do
    //     - Given <i| and |klj`> look up the matrix element <i|Z|klj`>
    //     - Allocate and fill the dense matrices <i|X|abc`> <abc`|X|klj`>, and same for Y
    //       but we only need to do this for one given <i| channel at a time.
    //
    //  map   i    -> index in dense matrix of X,Y,Z
    //  map klj`J  -> index in dense matrix of X,Y,Z
    //  map abc`J  -> index in dense matrix of X,Y  (we can reuse the map for klj`J
    //  map j_i,p_i,tz_i -> dense matrix <i|Z|klj`>
    //  map j_i,p_i,tz_1, index_i -> orbit index for i
    //  map j_i,p_i,tz_i, index_kljJ -> orbit indices k,l,j and Jkl
    //
    // We need additional bookkeeping machinery to store recoupling coefficients
    // This may be something to do in a second round of optimization.

    auto Hash_obc_index = [](int j2, int parity, int tz2)
    { return 2 * j2 - 1 + tz2 + parity; };

    int j2max_1b = 2 * Z.modelspace->GetEmax() + 1;
    size_t max_obc = Hash_obc_index(j2max_1b, 1, 1);

    size_t norb = Z.modelspace->GetNumberOrbits();
    std::vector<size_t> dense_index_i(norb, -1); // orbit index -> index in its corresponding dense matrix. Should correspond to radial qnumber n.

    std::vector<std::map<std::array<size_t, 4>, size_t>> dense_index_kljJ(max_obc + 1); //  3 orbit indices k,l,j and Jkl  -> index in corresponding dense matrix
    // and we want to do reverse lookup as well.
    std::vector<std::vector<size_t>> dense_i_lookup(max_obc + 1);                   // reverse lookup
    std::vector<std::vector<std::array<size_t, 4>>> dense_kljJ_lookup(max_obc + 1); // reverse lookup by one-body channel and dense index -> {k,l,j,Jkl}
    std::vector<std::vector<std::array<size_t, 4>>> dense_abcJ_lookup(max_obc + 1); // reverse lookup by one-body channel and dense index -> {k,l,j,Jkl}

    std::vector<arma::mat> ZMAT_list(max_obc + 1); // use a vector because we can use a perfect hash 2j2-1 + tz + parity

    for (auto &iter_obc_i : Z.modelspace->OneBodyChannels)
    {
      int j2_i = iter_obc_i.first[1];
      int parity_i = iter_obc_i.first[0] % 2; // it's stored as l for some reason.
      int tz2_i = iter_obc_i.first[2];

      // get the index of this one-body channel in our storage structures
      size_t obc_hash_i = Hash_obc_index(j2_i, parity_i, tz2_i);

      // loop over one-body states in this channel, count them and record their position in dense_index_i
      size_t index_i = 0;
      for (size_t i : iter_obc_i.second)
      {

        dense_index_i[i] = index_i;
        dense_i_lookup[obc_hash_i].push_back(i);
        index_i++;
      }

      // do the same for the pph states klj`J
      size_t index_klj = 0;
      size_t index_abc = 0;
      for (auto &iter_obc_j : Z.modelspace->OneBodyChannels)
      {
        int j2_j = iter_obc_j.first[1];
        int parity_j = iter_obc_j.first[0] % 2; // it's stored as l for some reason.
        int tz2_j = iter_obc_j.first[2];

        int parity_kl = (parity_i + parity_j + Z.GetParity()) % 2;
        int Tz_kl = (tz2_i + tz2_j) / 2 - Z.GetTRank();
        if (std::abs(Tz_kl) > 1)
          Tz_kl += 2 * Z.GetTRank(); // ambiguity because Trank=1 can change dTz by +- 1.

        size_t Jkl_min = std::abs(j2_i - j2_j) / 2;
        size_t Jkl_max = (j2_i + j2_j) / 2;

        for (size_t Jkl = Jkl_min; Jkl <= Jkl_max; Jkl++)
        {
          size_t ch_kl = Z.modelspace->GetTwoBodyChannelIndex(Jkl, parity_kl, Tz_kl);
          TwoBodyChannel &tbc_kl = Z.modelspace->GetTwoBodyChannel(ch_kl);
          size_t nkets_kl = tbc_kl.GetNumberKets();
          for (size_t iket = 0; iket < nkets_kl; iket++)
          {
            Ket &ket = tbc_kl.GetKet(iket);
            size_t k = ket.p;
            size_t l = ket.q;
            double na = ket.op->occ;
            double nb = ket.oq->occ;
            for (size_t j : iter_obc_j.second)
            {
              double nc = Z.modelspace->GetOrbit(j).occ;
              dense_index_kljJ[obc_hash_i][{k, l, j, Jkl}] = index_klj;
              dense_kljJ_lookup[obc_hash_i].push_back({k, l, j, Jkl});
              index_klj++;
              if ((na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc) > 1e-8) // Check that this has the right occupations to be an internal state
              {
                dense_abcJ_lookup[obc_hash_i].push_back({k, l, j, Jkl});
                index_abc++;
              }
            }
          }
        } // for Jkl
      }   // for iter_obc_j
          //     size_t zmat_index = 2*j2_i -1  + tz2_i + parity_i;
      ZMAT_list[obc_hash_i] = arma::mat(index_i - 0, index_klj - 0);
    } // for iter_obc_i  <= OneBodyChannels

    // Here we are outside any loops.

    // Now we loop through each dense matrix sub-block  <i|Z|klj> and fill the
    // necessary dense matrices for X and Y to be used in the matmult.
    for (auto &iter_obc_i : Z.modelspace->OneBodyChannels)
    {
      int j2_i = iter_obc_i.first[1];
      int parity_i = iter_obc_i.first[0] % 2; // it's stored as l for some reason.
      int tz2_i = iter_obc_i.first[2];
      size_t obc_hash_i = Hash_obc_index(j2_i, parity_i, tz2_i);
      double ji = 0.5 * j2_i;

      auto &ZMAT = ZMAT_list[obc_hash_i];
      //          |i
      //          |___,,     X
      //         a|  b| \c
      //          |___|__\   Y
      //          |   |   \
     //          |k  |l   \j

      // If Z is not channel diagonal, but X is, then there is a different relationship between
      // the one-body channel i and the pph channel klj. We can account for this by finding the
      // appropriate one-body channel to use in the lookup.
      // The j,p,tz of state klj is:   j_klj = j_i
      //                               parity_klj = (parity_i + Z.parity)%2
      //                    tk+tk-tj = tz_klj = tz_i +- Z.TRank.
      // Now the channel abc_X has      j_abc_X = j_i
      //                           parity_abc_X = (parity_i + X.parity)%2
      //                               tz_abc_X = tz_i +- X.Trank
      //
      // We want to find the one-body channel i_abc_X that, based on the symmetries of Z
      // would give us the channel abc_X. Equating the channels abc and klj above we have
      // j_iX = j_i
      // (parity_iX+Z.parity)%2 = (parity_i + X.parity)%2 => parity_iX = (parity_i + Z.parity+X.parirty)%2
      // tz_iX +- Z.Trank = tz_i +- X.Trank => tz_iX = tz_i +- (X.Trank - Z.Trank)
      //
      //        int parity_kl = ( parity_i+parity_j + Z.GetParity() )%2;
      //        int Tz_kl = (tz2_i + tz2_j)/2 - Z.GetTRank();
      //        if ( std::abs(Tz_kl) > 1 ) Tz_kl += 2*Z.GetTRank();
      int j2_iX = j2_i;
      int parity_iX = (parity_i + X.GetParity() + Z.GetParity()) % 2;
      int tz2_iX = tz2_i - (X.GetTRank() - Z.GetTRank()) * 2;
      if (std::abs(tz2_iX) != 1)
        tz2_iX = tz2_i + (X.GetTRank() - Z.GetTRank()) * 2;

      int j2_iY = j2_i;
      int parity_iY = (parity_i + Y.GetParity() + Z.GetParity()) % 2;
      int tz2_iY = tz2_i - (Y.GetTRank() - Z.GetTRank()) * 2;
      if (std::abs(tz2_iY) != 1)
        tz2_iY = tz2_i + (Y.GetTRank() - Z.GetTRank()) * 2;

      size_t obc_hash_iX = Hash_obc_index(j2_iX, parity_iX, tz2_iX);
      size_t obc_hash_iY = Hash_obc_index(j2_iY, parity_iY, tz2_iY);

      size_t dim_i = ZMAT.n_rows;
      size_t dim_klj = ZMAT.n_cols;
      //     size_t dim_abc_X = dense_abcJ_lookup[obc_hash_i].size();
      //     size_t dim_abc_Y = dim_abc_X;
      size_t dim_abc_X = dense_abcJ_lookup[obc_hash_iX].size();
      size_t dim_abc_Y = dense_abcJ_lookup[obc_hash_iY].size();
      arma::mat X2MAT(dim_i, dim_abc_X, arma::fill::zeros);
      arma::mat Y2MAT(dim_i, dim_abc_Y, arma::fill::zeros);
      arma::mat X3MAT(dim_abc_Y, dim_klj, arma::fill::zeros);
      arma::mat Y3MAT(dim_abc_X, dim_klj, arma::fill::zeros);

      // TODO: if case Z is not channel-diagonal, we need to be more careful with the abc states.
      // we will do that once the diagonal case is working...
      // XY_case=0 means we're considering the X2 * Y3 term, and XY_case=1 means the Y2 * X3 term
      // These can be different if X or Y is not channel diagonal.
      for (int XY_case = 0; XY_case < 2; XY_case++)
      {
        size_t dim_abc = XY_case == 0 ? dim_abc_X : dim_abc_Y;
        size_t obc_hash = XY_case == 0 ? obc_hash_iX : obc_hash_iY;
        if ( XY_case==0 and Y.ThreeBodyNorm()<1e-8)
           continue;
        if ( XY_case==1 and X.ThreeBodyNorm()<1e-8)
           continue;
        // fill <i|X|abc>
        for (size_t index_abc = 0; index_abc < dim_abc; index_abc++)
        {
          //        auto& abcJ = dense_abcJ_lookup[obc_hash_i].at(index_abc);
          auto &abcJ = dense_abcJ_lookup[obc_hash].at(index_abc);
          size_t a = abcJ[0];
          size_t b = abcJ[1];
          size_t c = abcJ[2];
          int Jab = abcJ[3];
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          Orbit &oc = Z.modelspace->GetOrbit(c);

          double ja = 0.5 * oa.j2;
          double jb = 0.5 * ob.j2;
          double jc = 0.5 * oc.j2;
          double occ_abc = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
          if (a == b)
            occ_abc *= 0.5; // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b

          // fill <i|X|abc>  The occupation factor goes with the 2-body matrix.
          for (size_t index_i = 0; index_i < dim_i; index_i++)
          {

            size_t i = dense_i_lookup[obc_hash_i].at(index_i);
            if (XY_case == 0)
            {
              double xiabc = -sqrt((2 * Jab + 1.)) * occ_abc * X.TwoBody.GetTBME_J(Jab, c, i, a, b);
              X2MAT(index_i, index_abc) = xiabc;
            }
            // THIS ONLY WORKS FOR CHANNEL DIAGONAL OPERATORS
            if (XY_case == 1)
            {
              double yiabc = -sqrt((2 * Jab + 1.)) * occ_abc * Y.TwoBody.GetTBME_J(Jab, c, i, a, b);
              Y2MAT(index_i, index_abc) = yiabc;
            }
          } // for index_i

          size_t twoJ = j2_i;
          for (size_t index_klj = 0; index_klj < dim_klj; index_klj++)
          {
            auto &kljJ = dense_kljJ_lookup[obc_hash_i].at(index_klj);
            size_t k = kljJ[0];
            size_t l = kljJ[1];
            size_t j = kljJ[2];
            int Jkl = kljJ[3];
            Orbit &ok = Z.modelspace->GetOrbit(k);
            Orbit &ol = Z.modelspace->GetOrbit(l);
            Orbit &oj = Z.modelspace->GetOrbit(j);
            double jk = 0.5 * ok.j2;
            double jl = 0.5 * ol.j2;
            double jj = 0.5 * oj.j2;

            double yabcklj = 0;
            double xabcklj = 0;

            size_t twoJp_min = std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * Jkl - oc.j2));
            size_t twoJp_max = std::min(2 * Jab + oj.j2, 2 * Jkl + oc.j2);

            for (size_t twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
            {
              double sixj = Z.modelspace->GetSixJ(jj, ji, Jkl, jc, 0.5 * twoJp, Jab);

              if (XY_case == 0)
              {
                bool y3_good = true;
                double yabjklc = 0;
                if (y3_good)
                  yabjklc = Y3.GetME_pn(Jab, Jkl, twoJp, a, b, j, k, l, c);
                yabcklj += (twoJp + 1.) / sqrt(2 * Jkl + 1.) * sixj * yabjklc;
              }
              if (XY_case == 1)
              {
                bool x3_good = true;
                double xabjklc = 0;
                if (x3_good)
                  xabjklc = X3.GetME_pn(Jab, Jkl, twoJp, a, b, j, k, l, c);
                xabcklj += (twoJp + 1.) / sqrt(2 * Jkl + 1.) * sixj * xabjklc;
              }
            } // for twoJp

            if (XY_case == 0)
              Y3MAT(index_abc, index_klj) = yabcklj;
            if (XY_case == 1)
              X3MAT(index_abc, index_klj) = xabcklj;

          } // for index_klj

        } // for index_abc_X
      }   // for XY_case

      // Now we do the mat mult.
      ZMAT = (X2MAT * Y3MAT - Y2MAT * X3MAT);

    } // for iter_obc_i  <= OneBodyChannels

    // Here we are outside of any loops.

    // Convert <i|Z|klj`> => <ij|Z|kl>  (with J coupling implied).
    // We need to access

    // It looks like nothing dangerous happens inside the loop, so we don't need to check for first_pass
    // Roll ch and ibra into a single loop for better load balancing
    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.TwoBody.MatEl)
    {
      size_t chbra = it.first[0];
      size_t chket = it.first[1];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(chbra);
      size_t nbras = tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras; ibra++)
      {
        bra_ket_channels.push_back({chbra, chket, ibra}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t nbra_ket_channels = bra_ket_channels.size();
    //  for ( size_t ch=0; ch<nch; ch++)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nbra_ket_channels; ich++)
    {
      size_t chbra = bra_ket_channels[ich][0];
      size_t chket = bra_ket_channels[ich][1];
      size_t ibra = bra_ket_channels[ich][2];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(chbra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(chket);
      //    size_t nkets = tbc.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      size_t J = tbc_bra.J;
      //    for (size_t ibra=0; ibra<nkets; ibra++)
      //    {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      if (imsrg3_valence_2b and (oi.cvq != 1 or oj.cvq != 1))
        continue;

      size_t iket_min = 0;
      if (chbra == chket)
        iket_min = ibra;
      for (size_t iket = iket_min; iket < nkets; iket++)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        size_t k = ket.p;
        size_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int phase_ij = X.modelspace->phase((oi.j2 + oj.j2) / 2 - J);
        int phase_kl = X.modelspace->phase((ok.j2 + ol.j2) / 2 - J);

        if (imsrg3_valence_2b and (ok.cvq != 1 or ol.cvq != 1))
          continue;

        // TODO probably can condense this into a for loop over i,j,k,l

        size_t obc_hash_i = Hash_obc_index(oi.j2, oi.l % 2, oi.tz2);
        size_t obc_hash_j = Hash_obc_index(oj.j2, oj.l % 2, oj.tz2);
        size_t obc_hash_k = Hash_obc_index(ok.j2, ok.l % 2, ok.tz2);
        size_t obc_hash_l = Hash_obc_index(ol.j2, ol.l % 2, ol.tz2);

        auto &ZMAT_i = ZMAT_list[obc_hash_i];
        auto &ZMAT_j = ZMAT_list[obc_hash_j];
        auto &ZMAT_k = ZMAT_list[obc_hash_k];
        auto &ZMAT_l = ZMAT_list[obc_hash_l];

        size_t ind_i = dense_index_i[i];
        size_t ind_j = dense_index_i[j];
        size_t ind_k = dense_index_i[k];
        size_t ind_l = dense_index_i[l];

        size_t ind_klj = -1;
        size_t ind_kli = -1;
        size_t ind_ijl = -1;
        size_t ind_ijk = -1;

        ind_klj = dense_index_kljJ[obc_hash_i][{k, l, j, J}];
        ind_kli = dense_index_kljJ[obc_hash_j][{k, l, i, J}]; // {3,5,3, 1}
        ind_ijl = dense_index_kljJ[obc_hash_k][{i, j, l, J}];
        ind_ijk = dense_index_kljJ[obc_hash_l][{i, j, k, J}];

        double Z_jkli = 0;
        double Z_iklj = 0;
        double Z_lijk = 0;
        double Z_kijl = 0;

        if (ind_klj < ZMAT_i.n_cols)
          Z_iklj = ZMAT_i(ind_i, ind_klj);
        if (ind_kli < ZMAT_j.n_cols)
          Z_jkli = ZMAT_j(ind_j, ind_kli);
        if (ind_ijl < ZMAT_k.n_cols)
          Z_kijl = ZMAT_k(ind_k, ind_ijl);
        if (ind_ijk < ZMAT_l.n_cols)
          Z_lijk = ZMAT_l(ind_l, ind_ijk);

        double zijkl = phase_ij * Z_iklj - Z_jkli + hermX * hermY * (Z_lijk - phase_kl * Z_kijl);

        // normalize the tbme
        zijkl /= sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));

        Z2.AddToTBME(chbra, chket, ibra, iket, zijkl);
      } // for iket
      //    }//for ibra
    } // for ch

    //  if (Z.modelspace->scalar3b_transform_first_pass)   Z.profiler.timer["comm232_first_pass"] += omp_get_wtime() - tstart;
//    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  // This optimized implementation is originally by Ragnar.
  void comm232ss_srs_optimized_old(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    //  auto& X2 = X.TwoBody;
    auto &X3 = X.ThreeBody;
    //  auto& Y2 = Y.TwoBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    bool x_has_3 = X3.IsAllocated();
    bool y_has_3 = Y3.IsAllocated();
    bool x_channel_diag = X.GetTRank() == 0 and X.GetParity() == 0;
    bool y_channel_diag = Y.GetTRank() == 0 and Y.GetParity() == 0;
    bool z_channel_diag = Z.GetTRank() == 0 and Z.GetParity() == 0;

    int hermX = X.IsHermitian() ? 1 : -1;
    int hermY = Y.IsHermitian() ? 1 : -1;

    Z.modelspace->PreCalculateSixJ(); // If we already did this, this does nothing.

    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    // Hash for lookup
    auto Hash_comm232_key = [](std::array<size_t, 5> &kljJJ)
    { return kljJJ[0] + (kljJJ[1] << 8) + (kljJJ[2] << 16) + (kljJJ[3] << 24) + (kljJJ[4] << 32); };

    size_t nch = Z.modelspace->GetNumberTwoBodyChannels();

    // first, enumerate the one-body channels |i> => ji,parityi,tzi
    std::map<std::array<int, 3>, std::vector<size_t>> local_one_body_channels;          //  maps {j,parity,tz} => vector of index_j
    std::map<std::array<int, 3>, std::vector<size_t>> external_local_one_body_channels; //  maps {j,parity,tz} => vector of index_j
    for (auto j : Z.modelspace->all_orbits)
    {
      Orbit &oj = Z.modelspace->GetOrbit(j);
      //    if (!Y3.IsOrbitIn3BodyEMaxTruncation(oj)) continue;
      std::array<int, 3> obc = {oj.j2, oj.l % 2, oj.tz2};
      if (local_one_body_channels.find(obc) == local_one_body_channels.end())
        local_one_body_channels[obc] = {j};
      else
        local_one_body_channels[obc].push_back(j);
      if (imsrg3_valence_2b and (oj.cvq != 1))
        continue;
      if (external_local_one_body_channels.find(obc) == external_local_one_body_channels.end())
        external_local_one_body_channels[obc] = {j};
      else
        external_local_one_body_channels[obc].push_back(j);
    }

    // TODO: This orgainization needs to be rethought to better accommodate isospin changing operators.
    //
    // next, figure out which three-body states |klj`> and |abc`> exist, make a list, and give them an
    // index for where they'll sit in the matrix
    // For an isospin-changing operator, there can be multiple possible channels for the three-body state
    // and the isospin projection can be +- 3/2, which is not possible for a single-particle state
    std::map<std::array<int, 3>, std::vector<std::array<size_t, 4>>> klj_list; //  maps {j,parity,tz} => {kljJ}
                                                                               //  std::map<std::array<int,3>,arma::mat> ZMAT_list; //  maps {j,parity,tz} => matrix <i|Z|klj`>
    std::map<std::array<int, 4>, arma::mat> ZMAT_list;                         //  maps {j,parity,tz,t2zklj} => matrix <i|Z|klj`>

    std::vector<std::array<int, 3>> obc_keys;
    for (auto &iter_i : external_local_one_body_channels)
      obc_keys.push_back(iter_i.first);
    size_t nkeys = obc_keys.size();

    for (size_t ikey = 0; ikey < nkeys; ikey++)
    {
      auto &obc_key = obc_keys[ikey]; // obc_key is {j2, parity, tz2} for a single orbit
      std::vector<size_t> &obc_orbits = external_local_one_body_channels[obc_key];
      int j2i = obc_key[0];
      int parityi = obc_key[1];
      int tz2i = obc_key[2];
      klj_list[obc_key] = {};
      auto &klj_list_i = klj_list[obc_key];

      for (auto &iter_j : local_one_body_channels)
      {
        int j2j = iter_j.first[0];
        int parityj = iter_j.first[1];
        int tz2j = iter_j.first[2];

        int Jkl_min = std::abs(j2i - j2j) / 2;
        int Jkl_max = (j2i + j2j) / 2;
        //      int Tzkl = (tz2i + tz2j)/2;
        //      int paritykl = (parityi+parityj)%2;
        int paritykl = (parityi + parityj + Z.GetParity()) % 2;
        for (int Tzkl = -1; Tzkl <= 1; Tzkl++)
        {
          if (std::abs((tz2i + tz2j) / 2 - Tzkl) != Z.GetTRank())
            continue;
          for (int Jkl = Jkl_min; Jkl <= Jkl_max; Jkl++)
          {

            size_t ch_kl = Z.modelspace->GetTwoBodyChannelIndex(Jkl, paritykl, Tzkl);
            TwoBodyChannel &tbc_kl = Z.modelspace->GetTwoBodyChannel(ch_kl);
            size_t nkets_kl = tbc_kl.GetNumberKets();
            for (size_t iket_kl = 0; iket_kl < nkets_kl; iket_kl++)
            {
              Ket &ket_kl = tbc_kl.GetKet(iket_kl);
              if (!Y3.IsOrbitIn3BodyEMaxTruncation(ket_kl.p))
                continue;
              if (!Y3.IsOrbitIn3BodyEMaxTruncation(ket_kl.q))
                continue;

              int e_k = 2 * ket_kl.op->n + ket_kl.op->l;
              int e_l = 2 * ket_kl.oq->n + ket_kl.oq->l;
              if ((e_k + e_l) > Z.modelspace->GetE3max())
                continue;

              double de_k = std::abs(e_k - e_fermi[ket_kl.op->tz2]);
              double de_l = std::abs(e_l - e_fermi[ket_kl.oq->tz2]);
              if ((de_k + de_l) > Z.modelspace->GetdE3max())
                continue;
              //          double occnat_k = ket_kl.op->occ_nat;
              //          double occnat_l = ket_kl.oq->occ_nat;
              for (size_t j : iter_j.second)
              {
                //            if (!Y3.IsOrbitIn3BodyEMaxTruncation(j)) continue;
                // if (!Z.ThreeBody.IsKetInEMaxTruncations(ket_kl.p, ket_kl.q, j)) continue;
                Orbit &oj = Z.modelspace->GetOrbit(j);
                if (!Y3.IsOrbitIn3BodyEMaxTruncation(oj))
                  continue;
                //            double occnat_j = oj.occ_nat;
                //            if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_j*(1-occnat_j) ) < Z.modelspace->GetOccNat3Cut() ) continue;

                klj_list_i.push_back({ket_kl.p, ket_kl.q, j, (size_t)tbc_kl.J});

              } // for j
            }   // for iket_kl
          }     // for Jkl
        }       // for Tzkl
      }         // for iter_j in local one body channels
      if (z_channel_diag)
      {
        size_t dim_i = obc_orbits.size();          // how many 1b states are in this jpt channel
        size_t dim_klj = klj_list[obc_key].size(); // how many 2p1h states are in this jpt channel
                                                   //       ZMAT_list[obc_key] =  arma::mat(dim_i,dim_klj) ;  //allocate the matrix <i|Z|klj`>
        std::array<int, 4> zmat_key = {j2i, parityi, tz2i, tz2i};
        //       ZMAT_list[obc_key] =  arma::mat(dim_i,dim_klj) ;  //allocate the matrix <i|Z|klj`>
        ZMAT_list[zmat_key] = arma::mat(dim_i, dim_klj); // allocate the matrix <i|Z|klj`>
      }

    } // for ikey

    if (not z_channel_diag)
    {
      for (size_t ikey = 0; ikey < nkeys; ikey++)
      {
        auto &obc_key_i = obc_keys[ikey];
        int j2i = obc_key_i[0];
        int parityi = obc_key_i[1];
        int tz2i = obc_key_i[2];
        size_t dim_i = external_local_one_body_channels[obc_key_i].size();
        for (size_t kljkey = 0; kljkey < nkeys; kljkey++)
        {
          auto &obc_key_klj = obc_keys[kljkey];
          int j2_klj = obc_key_klj[0];
          int parity_klj = obc_key_klj[1];
          int tz2_klj = obc_key_klj[2];

          if (j2i != j2_klj)
            continue; // we're only handling scalar operators
          if ((parityi + parity_klj % 2) != Z.GetParity())
            continue;
          if (std::abs(tz2i - tz2_klj) / 2 != Z.GetTRank())
            continue;
          // for an isospin changing operator, it's possible that two klj channels can connect with channel i.
          //          auto& klj_list_i = klj_list[obc_key_i];
          size_t dim_klj = klj_list[obc_key_klj].size(); // how many 2p1h states are in this jpt channel
          std::array<int, 4> zmat_key = {j2i, parityi, tz2i, tz2_klj};
          //          ZMAT_list[obc_key_i] =  arma::mat(dim_i,dim_klj) ;  //allocate the matrix <i|Z|klj`>
          ZMAT_list[zmat_key] = arma::mat(dim_i, dim_klj); // allocate the matrix <i|Z|klj`>
        }                                                  // for kljkey
      }                                                    // for ikey
    }

// This loop is what takes all the time.
// Outer loop over one-body quantum numbers.
// For the usual situation, these same quantum numbers correspond to |klj`> and |abc`>.
// But if X or Y changes parity or isospin, then the three-body states will correspond to a different set of one-body quantum numbers.
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ikey = 0; ikey < nkeys; ikey++)
    {

      auto &obc_key_i = obc_keys[ikey]; // this is an array {j2,parity,tz2}
      int j2i = obc_key_i[0];
      int parityi = obc_key_i[1];
      int tz2i = obc_key_i[2];
      std::vector<size_t> &obc_orbits = external_local_one_body_channels[obc_key_i]; // the orbits that have quantum numbers {j2,parity,tz2}
      auto &klj_list_i = klj_list[obc_key_i];                                        // list of 3-body pph states |klj`> with quantum numbers JPTz such that <i|Z|klj`> is nonzero

      // identify the possible quantum numbers for the abc channel.
      // they are determined by the pph channels connected to i by X and Y, i.e. <i|X|abc> and <i|Y|abc>

      for (int parity_klj : {0, 1})
      {
        for (int tz2_klj : {-3, -1, 1, 3})
        {
          if ((parity_klj + parityi) % 2 != Z.GetParity())
            continue;
          if (std::abs(tz2i - tz2_klj) / 2 != Z.GetTRank())
            continue;

          std::array<int, 3> obc_key_klj = {j2i, parity_klj, tz2_klj};
          //       auto& klj_list_i = klj_list[obc_key]; // list of 3-body pph states |klj`> with quantum numbers JPTz such that <i|Z|klj`> is nonzero
          auto &klj_list_klj = klj_list[obc_key_klj]; // list of 3-body pph states |klj`> with quantum numbers JPTz such that <i|Z|klj`> is nonzero

          std::vector<size_t> abc_list_i; // keep track of which klj states should go in the abc loop
          std::vector<double> abc_occ_list_i;
          std::vector<size_t> abc_list_klj; // keep track of which klj states should go in the abc loop
          std::vector<double> abc_occ_list_klj;
          for (int abc_loop = 0; abc_loop <= 1; abc_loop++)
          {
            auto &key_this_loop = abc_loop == 0 ? obc_key_i : obc_key_klj;
            if (z_channel_diag and abc_loop > 0)
              continue;
            for (size_t i_kljJ = 0; i_kljJ < klj_list[key_this_loop].size(); i_kljJ++)
            {
              auto &kljJ = klj_list[key_this_loop][i_kljJ];
              size_t a = kljJ[0];
              size_t b = kljJ[1];
              size_t c = kljJ[2];
              Orbit &oa = Z.modelspace->GetOrbit(a);
              Orbit &ob = Z.modelspace->GetOrbit(b);
              Orbit &oc = Z.modelspace->GetOrbit(c);
              //         if (!Y3.IsOrbitIn3BodyEMaxTruncation(oa)) continue;
              //         if (!Y3.IsOrbitIn3BodyEMaxTruncation(ob)) continue;
              //         if (!Y3.IsOrbitIn3BodyEMaxTruncation(oc)) continue;
              double e_a = 2 * oa.n + oa.l;
              double e_b = 2 * ob.n + ob.l;
              double de_a = std::abs(e_a - e_fermi[oa.tz2]);
              double de_b = std::abs(e_b - e_fermi[ob.tz2]);
              if ((e_a + e_b) > Z.modelspace->GetE3max())
                continue;
              if ((de_a + de_b) > Z.modelspace->GetdE3max())
                continue;
              double na = oa.occ;
              double nb = ob.occ;
              double nc = oc.occ;
              double occupation_factor = na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc;
              if (std::abs(occupation_factor) < 1e-6)
                continue;
              if (a == b)
                occupation_factor *= 0.5; // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b
              if (abc_loop == 0)
              {
                abc_list_i.push_back(i_kljJ);
                abc_occ_list_i.push_back(occupation_factor);
              }
              else
              {
                abc_list_klj.push_back(i_kljJ);
                abc_occ_list_klj.push_back(occupation_factor);
              }
            }
          } // for abc_loop

          size_t dim_i = obc_orbits.size();              // how many sp states are in this jpt channel
                                                         //       size_t dim_klj = klj_list[obc_key].size(); // how many 3-body pph states in this jpt channel
          size_t dim_klj = klj_list[obc_key_klj].size(); // how many 3-body pph states in this jpt channel
                                                         //       size_t dim_abc = abc_list.size(); // how many 3-body states which contribute to the |abc`> sum
          size_t dim_abc_i = abc_list_i.size();          // how many 3-body states which contribute to the |abc`> sum
          size_t dim_abc_klj = abc_list_klj.size();      // how many 3-body states which contribute to the |abc`> sum

          if (x_channel_diag and y_channel_diag)
            dim_abc_klj = dim_abc_i;

          // Now allocate the matrices in this channel
          arma::mat X2MAT, Y2MAT, X3MAT, Y3MAT;
          if (x_channel_diag)
          {
            X2MAT.zeros(dim_i, dim_abc_i);
            X3MAT.zeros(dim_abc_klj, dim_klj);
          }
          else
          {
            X2MAT.zeros(dim_i, dim_abc_klj);
            X3MAT.zeros(dim_abc_i, dim_klj);
          }
          if (y_channel_diag)
          {
            Y2MAT.zeros(dim_i, dim_abc_i);
            Y3MAT.zeros(dim_abc_klj, dim_klj);
          }
          else
          {
            Y2MAT.zeros(dim_i, dim_abc_klj);
            Y3MAT.zeros(dim_abc_i, dim_klj);
          }

          //       // Now allocate the matrices in this channel
          //       arma::mat X2MAT(dim_i,   dim_abc, arma::fill::zeros);
          //       arma::mat Y2MAT(dim_i,   dim_abc, arma::fill::zeros);
          //       arma::mat X3MAT(dim_abc, dim_klj, arma::fill::zeros);
          //       arma::mat Y3MAT(dim_abc, dim_klj, arma::fill::zeros);

          // figure out which recouplings we'll need when filling the matrices
          std::set<std::array<size_t, 5>> kljJJ_needed; // we use a set to avoid repeats

          for (int abc_loop = 0; abc_loop <= 1; abc_loop++)
          {
            if (z_channel_diag and abc_loop > 0)
              continue;
            int dim_abc = abc_loop == 0 ? dim_abc_i : dim_abc_klj;
            auto &abc_list = abc_loop == 0 ? abc_list_i : abc_list_klj;
            auto &klj_list_thisloop = abc_loop == 0 ? klj_list_i : klj_list_klj;

            for (size_t ind_abc = 0; ind_abc < dim_abc; ind_abc++)
            {
              //            auto& abcJ = klj_list_i[ abc_list[ind_abc] ];
              auto &abcJ = klj_list_thisloop[abc_list[ind_abc]];
              size_t a = abcJ[0];
              size_t b = abcJ[1];
              size_t c = abcJ[2];
              int Jab = abcJ[3];
              Orbit &oa = Z.modelspace->GetOrbit(a);
              Orbit &ob = Z.modelspace->GetOrbit(b);
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double e_a = 2 * oa.n + oa.l;
              double e_b = 2 * ob.n + ob.l;
              double e_c = 2 * oc.n + oc.l;
              double de_a = std::abs(e_a - e_fermi[oa.tz2]);
              double de_b = std::abs(e_b - e_fermi[ob.tz2]);
              double de_c = std::abs(e_c - e_fermi[oc.tz2]);
              double occnat_a = oa.occ_nat;
              double occnat_b = ob.occ_nat;
              double occnat_c = oc.occ_nat;
              int j2c = oc.j2;
              for (size_t ind_klj = 0; ind_klj < dim_klj; ind_klj++)
              {
                //              auto& kljJ = klj_list_i[ ind_klj ];
                auto &kljJ = klj_list_klj[ind_klj];
                size_t k = kljJ[0];
                size_t l = kljJ[1];
                size_t j = kljJ[2];
                int Jkl = kljJ[3];
                Orbit &oj = Z.modelspace->GetOrbit(j);
                Orbit &ok = Z.modelspace->GetOrbit(k);
                Orbit &ol = Z.modelspace->GetOrbit(l);

                double e_j = 2 * oj.n + oj.l;
                double e_k = 2 * ok.n + ok.l;
                double e_l = 2 * ol.n + ol.l;
                double de_j = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
                double de_k = std::abs(2 * ok.n + ok.l - e_fermi[ok.tz2]);
                double de_l = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);

                double occnat_j = oj.occ_nat;
                double occnat_k = ok.occ_nat;
                double occnat_l = ol.occ_nat;
                if ((e_k + e_l + e_c) > Z.modelspace->GetE3max())
                  continue;
                if ((e_a + e_b + e_j) > Z.modelspace->GetE3max())
                  continue;
                if ((de_k + de_l + de_c) > Z.modelspace->GetdE3max())
                  continue;
                if ((de_a + de_b + de_j) > Z.modelspace->GetdE3max())
                  continue;
                if ((occnat_k * (1 - occnat_k) * occnat_l * (1 - occnat_l) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
                  continue;
                if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) < Z.modelspace->GetOccNat3Cut())
                  continue;
                if (imsrg3_no_qqq and ((ok.cvq + ol.cvq + oc.cvq) > 5 or (oa.cvq + ob.cvq + oj.cvq) > 5))
                  continue;
                if (imsrg3_valence_2b and (ok.cvq != 1 or ol.cvq != 1 or oj.cvq != 1))
                  continue;
                int twoJp_min = std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * Jkl - j2c));
                int twoJp_max = std::min(2 * Jab + oj.j2, 2 * Jkl + j2c);
                for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                {
                  kljJJ_needed.insert({k, l, c, (size_t)Jkl, (size_t)twoJp});
                  kljJJ_needed.insert({a, b, j, (size_t)Jab, (size_t)twoJp});
                } // for twoJp
              }   // for ind_klj
            }     // for ind_abc
          }       // for abc_loop

          // now compute the couplings once and store them
          // this is beneficial because the same recoupling happens on the bra and ket sides, and this way we avoid doing the same recoupling multiple times

          std::vector<size_t> recouple_info;                          // channel, start_pointer, n_terms
          std::vector<size_t> recouple_kets;                          // concatenated list of ket indices
          std::vector<double> recouple_coeff;                         // concatenated list of recoupling coefficients
          std::unordered_map<size_t, size_t> recoupling_cache_lookup; // a lookup table so we can find the info easily

          for (auto kljJJ : kljJJ_needed)
          {
            size_t k = kljJJ[0];
            size_t l = kljJJ[1];
            size_t j = kljJJ[2];
            size_t Jkl = kljJJ[3];
            size_t twoJp = kljJJ[4];
            size_t hash = Hash_comm232_key(kljJJ);

            if (!Y.ThreeBody.IsKetValid(Jkl, twoJp, k, l, j))
              continue;
            std::vector<size_t> iketlist;
            std::vector<double> recouplelist;
            size_t ch_check = Y.ThreeBody.GetKetIndex_withRecoupling(Jkl, twoJp, k, l, j, iketlist, recouplelist);

            recoupling_cache_lookup[hash] = recouple_info.size();
            recouple_info.push_back(ch_check);
            recouple_info.push_back(recouple_coeff.size());
            recouple_info.push_back(iketlist.size());
            //     recouple_info.push_back( recoupling_cache.size() );

            for (size_t ik : iketlist)
              recouple_kets.push_back(ik); // indices of kets to loop over
            for (double rec : recouplelist)
              recouple_coeff.push_back(rec); // recoupling coefficents (sixJ etc)
          }

          // Now fill the matrices

          for (int abc_loop = 0; abc_loop <= 1; abc_loop++)
          {
            if (z_channel_diag and abc_loop > 0)
              continue;
            int dim_abc = abc_loop == 0 ? dim_abc_i : dim_abc_klj;
            auto &abc_list = abc_loop == 0 ? abc_list_i : abc_list_klj;
            auto &abc_occ_list = abc_loop == 0 ? abc_occ_list_i : abc_occ_list_klj;
            auto &klj_list_thisloop = abc_loop == 0 ? klj_list_i : klj_list_klj;

            for (size_t ind_abc = 0; ind_abc < dim_abc; ind_abc++) // rows of X2MAT
            {

              //            auto& abcJ = klj_list_i[ abc_list[ind_abc] ];
              auto &abcJ = klj_list_thisloop[abc_list[ind_abc]];
              size_t a = abcJ[0];
              size_t b = abcJ[1];
              size_t c = abcJ[2];
              int Jab = abcJ[3];
              Orbit &oa = Z.modelspace->GetOrbit(a);
              Orbit &ob = Z.modelspace->GetOrbit(b);
              Orbit &oc = Z.modelspace->GetOrbit(c);
              int j2c = oc.j2;
              double jc = 0.5 * j2c;
              double occ_abc = abc_occ_list[ind_abc]; // TODO: we can probably incorporate that hat factor with the occupation factor

              double de_a = std::abs(2 * oa.n + oa.l - e_fermi[oa.tz2]);
              double de_b = std::abs(2 * ob.n + ob.l - e_fermi[ob.tz2]);
              double de_c = std::abs(2 * oc.n + oc.l - e_fermi[oc.tz2]);
              double occnat_a = oa.occ_nat;
              double occnat_b = ob.occ_nat;
              double occnat_c = oc.occ_nat;

              // fill the 2 body part
              for (size_t ind_i = 0; ind_i < dim_i; ind_i++) // ind_i indexes the list of sp states in this one body channel. columns of X2MAT
              {
                size_t i = obc_orbits[ind_i]; // i is the orbit index as per ModelSpace

                if ((x_channel_diag and abc_loop == 0) or ((not x_channel_diag) and abc_loop == 1))
                {
                  X2MAT(ind_i, ind_abc) = -sqrt((2 * Jab + 1.)) * occ_abc * X.TwoBody.GetTBME_J(Jab, c, i, a, b);
                }
                if ((y_channel_diag and abc_loop == 0) or ((not y_channel_diag) and abc_loop == 1))
                {
                  Y2MAT(ind_i, ind_abc) = -sqrt((2 * Jab + 1.)) * occ_abc * Y.TwoBody.GetTBME_J(Jab, c, i, a, b);
                }
              } // for ind_i

              // now the 3 body part
              for (size_t ind_klj = 0; ind_klj < dim_klj; ind_klj++)
              {
                // <abi Jab twoJ | X | klc Jkl twoJ >
                //              auto& kljJ = klj_list_i[ ind_klj ];
                auto &kljJ = klj_list_klj[ind_klj];
                size_t k = kljJ[0];
                size_t l = kljJ[1];
                size_t j = kljJ[2];
                int Jkl = kljJ[3];
                Orbit &oj = Z.modelspace->GetOrbit(j);
                Orbit &ok = Z.modelspace->GetOrbit(k);
                Orbit &ol = Z.modelspace->GetOrbit(l);
                double ji = 0.5 * j2i;
                double jj = 0.5 * oj.j2;
                //           double jk = 0.5 * ok.j2;
                //           double jl = 0.5 * ol.j2;

                double de_j = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
                double de_k = std::abs(2 * ok.n + ok.l - e_fermi[ok.tz2]);
                double de_l = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);
                double occnat_j = oj.occ_nat;
                double occnat_k = ok.occ_nat;
                double occnat_l = ol.occ_nat;
                // Are these checks redundant?
                if ((de_a + de_b + de_j) > Z.modelspace->GetdE3max())
                  continue;
                if ((de_k + de_l + de_c) > Z.modelspace->GetdE3max())
                  continue;
                if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) < Z.modelspace->GetOccNat3Cut())
                  continue;
                if ((occnat_k * (1 - occnat_k) * occnat_l * (1 - occnat_l) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
                  continue;
                if (imsrg3_no_qqq and ((ok.cvq + ol.cvq + oc.cvq) > 5 or (oa.cvq + ob.cvq + oj.cvq) > 5))
                  continue;
                if (imsrg3_valence_2b and (ok.cvq != 1 or ol.cvq != 1 or oj.cvq != 1))
                  continue;

                int dTz = (oa.tz2 + ob.tz2 + oj.tz2 - ok.tz2 - ol.tz2 - oc.tz2);
                int dparity = (oa.l + ob.l + oj.l + ok.l + ol.l + oc.l) % 2;
                bool x3_good = x_has_3 and (std::abs(dTz) == X.GetTRank()) and (dparity == X.GetParity());
                bool y3_good = y_has_3 and (std::abs(dTz) == Y.GetTRank()) and (dparity == Y.GetParity());

                int twoJp_min = std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * Jkl - j2c));
                int twoJp_max = std::min(2 * Jab + oj.j2, 2 * Jkl + j2c);

                for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                {
                  if (!Y.ThreeBody.IsKetValid(Jkl, twoJp, k, l, c))
                    continue;
                  if (!Y.ThreeBody.IsKetValid(Jab, twoJp, a, b, j))
                    continue;
                  //              double sixj = X.modelspace->GetSixJ(Jkl,jj,ji,  Jab, jc, 0.5*twoJp );
                  double sixj = X.modelspace->GetSixJ(jj, ji, Jkl, jc, 0.5 * twoJp, Jab);

                  std::array<size_t, 5> abjJJ = {a, b, j, (size_t)Jab, (size_t)twoJp};
                  std::array<size_t, 5> klcJJ = {k, l, c, (size_t)Jkl, (size_t)twoJp};
                  size_t hash_abj = Hash_comm232_key(abjJJ);
                  size_t hash_klc = Hash_comm232_key(klcJJ);

                  size_t pointer_abj = recoupling_cache_lookup.at(hash_abj);
                  size_t pointer_klc = recoupling_cache_lookup.at(hash_klc);
                  //              size_t ch_check = recoupling_cache[pointer_abj];

                  size_t ch_abj = recouple_info[pointer_abj];
                  size_t ch_klc = recouple_info[pointer_klc];

                  size_t rec_pointer_abj = recouple_info[pointer_abj + 1];
                  size_t rec_pointer_klc = recouple_info[pointer_klc + 1];

                  size_t listsize_abj = recouple_info[pointer_abj + 2];
                  size_t listsize_klc = recouple_info[pointer_klc + 2];

                  //              size_t ch_abj = recoupling_cache[pointer_abj];
                  //              size_t ch_klc = recoupling_cache[pointer_klc];
                  //              size_t listsize_abj = recoupling_cache[pointer_abj+1];
                  //              size_t listsize_klc = recoupling_cache[pointer_klc+1];

                  double xabjklc = 0;
                  double yabjklc = 0;
                  for (size_t ilist_abj = 0; ilist_abj < listsize_abj; ilist_abj++)
                  {
                    size_t ibra_abj = recouple_kets[rec_pointer_abj + ilist_abj];
                    double recouple_bra = recouple_coeff[rec_pointer_abj + ilist_abj];
                    //                size_t ibra_abj = (size_t) recoupling_cache[pointer_abj+2+ilist_abj];
                    //                double recouple_bra =      recoupling_cache[pointer_abj+2+listsize_abj+ilist_abj];
                    for (size_t ilist_klc = 0; ilist_klc < listsize_klc; ilist_klc++)
                    {

                      size_t iket_klc = recouple_kets[rec_pointer_klc + ilist_klc];
                      double recouple_ket = recouple_coeff[rec_pointer_klc + ilist_klc];
                      //                  size_t iket_klc = (size_t) recoupling_cache[pointer_klc+2+ilist_klc];
                      //                  double recouple_ket =      recoupling_cache[pointer_klc+2+listsize_klc+ilist_klc];

                      //                  xabjklc += recouple_bra*recouple_ket * X3.GetME_pn_ch(ch_check,ch_check, ibra_abj, iket_klc );
                      if (x3_good)
                        xabjklc += recouple_bra * recouple_ket * X3.GetME_pn_ch(ch_abj, ch_klc, ibra_abj, iket_klc);
                      if (y3_good)
                        yabjklc += recouple_bra * recouple_ket * Y3.GetME_pn_ch(ch_abj, ch_klc, ibra_abj, iket_klc);
                    }
                  }

                  // if x is channel diagonal, the channel of abc should match klj. This happens if abc_loop==1.
                  // otherwise, we want abc_loop=1. This omits the case where neither X nor Y are channel diagonal.
                  // Ideally, we should handle that case as well, but I'm hoping I rewrite this routine before that comes up.
                  if ((x_channel_diag and abc_loop == 1) or (y_channel_diag and abc_loop == 0))
                  {
                    X3MAT(ind_abc, ind_klj) += (twoJp + 1.) / sqrt(2 * Jkl + 1.) * sixj * xabjklc;
                  }
                  if ((y_channel_diag and abc_loop == 1) or (x_channel_diag and abc_loop == 0))
                  {
                    Y3MAT(ind_abc, ind_klj) += (twoJp + 1.) / sqrt(2 * Jkl + 1.) * sixj * yabjklc;
                  }

                } // for twoJp
              }   // for ind_klj
            }     // for ind_abc
          }       // for abc_loop

          // now we do the mat mult
          //    ZMAT_list[obc_key] =  (  X2MAT * Y3MAT - Y2MAT * X3MAT  ) ;
          std::array<int, 4> zmat_key = {j2i, parityi, tz2i, tz2_klj};
          ZMAT_list[zmat_key] = (X2MAT * Y3MAT - Y2MAT * X3MAT);

        } // for tz2_klj;
      }   // for parity_klj

    } // for iter_i in local one body channels

    // now we need to unpack all that mess and store it in Z
    // It looks like nothing dangerous happens inside the loop, so we don't need to check for first_pass
    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
    // Roll ch and ibra into a single loop for better load balancing
    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.TwoBody.MatEl)
    {
      size_t chbra = it.first[0];
      size_t chket = it.first[1];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(chbra);
      size_t nbras = tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras; ibra++)
      {
        bra_ket_channels.push_back({chbra, chket, ibra}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t nbra_ket_channels = bra_ket_channels.size();
    //  for ( size_t ch=0; ch<nch; ch++)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nbra_ket_channels; ich++)
    {
      size_t chbra = bra_ket_channels[ich][0];
      size_t chket = bra_ket_channels[ich][1];
      size_t ibra = bra_ket_channels[ich][2];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(chbra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(chket);
      //    size_t nkets = tbc.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      size_t J = tbc_bra.J;
      //    for (size_t ibra=0; ibra<nkets; ibra++)
      //    {
      Ket &bra = tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      if (imsrg3_valence_2b and (oi.cvq != 1 or oj.cvq != 1))
        continue;

      //      auto& lobc_i = local_one_body_channels.at({oi.j2,oi.l%2,oi.tz2});
      //      auto& lobc_j = local_one_body_channels.at({oj.j2,oj.l%2,oj.tz2});
      auto &lobc_i = external_local_one_body_channels.at({oi.j2, oi.l % 2, oi.tz2});
      auto &lobc_j = external_local_one_body_channels.at({oj.j2, oj.l % 2, oj.tz2});
      size_t ind_i = std::distance(lobc_i.begin(), std::find(lobc_i.begin(), lobc_i.end(), i));
      size_t ind_j = std::distance(lobc_j.begin(), std::find(lobc_j.begin(), lobc_j.end(), j));

      //      auto& list_i = klj_list.at({oi.j2,oi.l%2,oi.tz2});
      //      auto& list_j = klj_list.at({oj.j2,oj.l%2,oj.tz2});
      //      auto& ZMat_i = ZMAT_list.at({oi.j2,oi.l%2,oi.tz2});
      //      auto& ZMat_j = ZMAT_list.at({oj.j2,oj.l%2,oj.tz2});
      //      int tz2_klj = ok.tz2+ol.tz2-oj.tz2;
      //      int tz2_kli = ok.tz2+ol.tz2-oi.tz2;
      //      auto& ZMat_i = ZMAT_list.at({oi.j2,oi.l%2,oi.tz2, tz2_klj});
      //      auto& ZMat_j = ZMAT_list.at({oj.j2,oj.l%2,oj.tz2, tz2_kli});

      size_t iket_min = 0;
      if (chbra == chket)
        iket_min = ibra;
      for (size_t iket = iket_min; iket < nkets; iket++)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        size_t k = ket.p;
        size_t l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        int phase_ij = X.modelspace->phase((oi.j2 + oj.j2) / 2 - J);
        int phase_kl = X.modelspace->phase((ok.j2 + ol.j2) / 2 - J);

        if (imsrg3_valence_2b and (ok.cvq != 1 or ol.cvq != 1))
          continue;

        //        auto& lobc_k = local_one_body_channels.at({ok.j2,ok.l%2,ok.tz2});
        //        auto& lobc_l = local_one_body_channels.at({ol.j2,ol.l%2,ol.tz2});
        auto &lobc_k = external_local_one_body_channels.at({ok.j2, ok.l % 2, ok.tz2});
        auto &lobc_l = external_local_one_body_channels.at({ol.j2, ol.l % 2, ol.tz2});
        size_t ind_k = std::distance(lobc_k.begin(), std::find(lobc_k.begin(), lobc_k.end(), k));
        size_t ind_l = std::distance(lobc_l.begin(), std::find(lobc_l.begin(), lobc_l.end(), l));

        //        auto& list_k = klj_list.at({ok.j2,ok.l%2,ok.tz2});
        //        auto& list_l = klj_list.at({ol.j2,ol.l%2,ol.tz2});
        //        auto& ZMat_k = ZMAT_list.at({ok.j2,ok.l%2,ok.tz2});
        //        auto& ZMat_l = ZMAT_list.at({ol.j2,ol.l%2,ol.tz2});

        int tz2_klj = ok.tz2 + ol.tz2 - oj.tz2;
        int tz2_kli = ok.tz2 + ol.tz2 - oi.tz2;
        int tz2_ijl = oi.tz2 + oj.tz2 - ol.tz2;
        int tz2_ijk = oi.tz2 + oj.tz2 - ok.tz2;
        int parity_klj = (ok.l + ol.l + oj.l) % 2;
        int parity_kli = (ok.l + ol.l + oi.l) % 2;
        int parity_ijk = (oi.l + oj.l + ok.l) % 2;
        int parity_ijl = (oi.l + oj.l + ol.l) % 2;

        auto &list_i = klj_list.at({oi.j2, parity_klj, tz2_klj});
        auto &list_j = klj_list.at({oj.j2, parity_kli, tz2_kli});
        auto &list_k = klj_list.at({ok.j2, parity_ijl, tz2_ijl});
        auto &list_l = klj_list.at({ol.j2, parity_ijk, tz2_ijk});

        auto &ZMat_i = ZMAT_list.at({oi.j2, oi.l % 2, oi.tz2, tz2_klj});
        auto &ZMat_j = ZMAT_list.at({oj.j2, oj.l % 2, oj.tz2, tz2_kli});

        auto &ZMat_k = ZMAT_list.at({ok.j2, ok.l % 2, ok.tz2, tz2_ijl});
        auto &ZMat_l = ZMAT_list.at({ol.j2, ol.l % 2, ol.tz2, tz2_ijk});

        std::array<size_t, 4> key_klj = {k, l, j, J};
        std::array<size_t, 4> key_kli = {k, l, i, J};
        std::array<size_t, 4> key_ijk = {i, j, k, J};
        std::array<size_t, 4> key_ijl = {i, j, l, J};

        size_t ind_klj = std::distance(list_i.begin(), std::find(list_i.begin(), list_i.end(), key_klj));
        size_t ind_kli = std::distance(list_j.begin(), std::find(list_j.begin(), list_j.end(), key_kli));
        size_t ind_ijl = std::distance(list_k.begin(), std::find(list_k.begin(), list_k.end(), key_ijl));
        size_t ind_ijk = std::distance(list_l.begin(), std::find(list_l.begin(), list_l.end(), key_ijk));
        double Z_jkli = 0;
        double Z_iklj = 0;
        double Z_lijk = 0;
        double Z_kijl = 0;
        //        std::cout << "ijkl " << i << " " << j << " " << k << " " << l << std::endl;
        //        std::cout << "tz2_xxx" << tz2_klj << " " << tz2_kli << " " << tz2_ijl << " " << tz2_ijk << std::endl;
        //        std::cout << "dimensions " << ZMat_i.n_rows << " " << ZMat_i.n_cols << " , "
        //                                   << ZMat_j.n_rows << " " << ZMat_j.n_cols << " , "
        //                                   << ZMat_k.n_rows << " " << ZMat_k.n_cols << " , "
        //                                   << ZMat_l.n_rows << " " << ZMat_l.n_cols << std::endl;
        //        std::cout << "indices  " << ind_i << " " << ind_klj << " , "
        //                                 << ind_j << " " << ind_kli << " , "
        //                                 << ind_k << " " << ind_ijl << " , "
        //                                 << ind_l << " " << ind_ijk << std::endl;
        if (ind_klj < list_i.size())
          Z_iklj = ZMat_i(ind_i, ind_klj);
        if (ind_kli < list_j.size())
          Z_jkli = ZMat_j(ind_j, ind_kli);
        if (ind_ijl < list_k.size())
          Z_kijl = ZMat_k(ind_k, ind_ijl);
        if (ind_ijk < list_l.size())
          Z_lijk = ZMat_l(ind_l, ind_ijk);

        double zijkl = phase_ij * Z_iklj - Z_jkli + hermX * hermY * (Z_lijk - phase_kl * Z_kijl);

        // normalize the tbme
        zijkl /= sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));

        //        zijkl *= -1.0 / sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
        Z2.AddToTBME(chbra, chket, ibra, iket, zijkl);
      } // for iket
      //    }//for ibra
    } // for ch
      //  if (Z.modelspace->scalar3b_transform_first_pass)   Z.profiler.timer["comm232_first_pass"] += omp_get_wtime() - tstart;
    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  // This is the comm232ss implementation by Matthias.
  // "Expand" refers to the fact that the 3BMEs are expanded into a full set
  // channel by channel to optimize the contraction.
  //
  // In the 232 commutator, if emax != emax_3body (the truncation for the 3B operator)
  // one of the external indices is in the emax model space
  // instead of the emax_3body model space. This implementation accounts for that properly
  // rather than simply truncating all indices to the emax_3body level.
  // This comes at a substantial computational cost, especially for large emax.
//  void comm232ss_expand_full(const Operator &X, const Operator &Y, Operator &Z)
//  {
//    comm232::comm232ss_expand_impl_full(X, Y, Z);
//  }

  // This is the comm232ss implementation by Matthias.
  // "Expand" refers to the fact that the 3BMEs are expanded into a full set
  // channel by channel to optimize the contraction.
  //
  // In the 232 commutator, if emax != emax_3body (the truncation for the 3B operator)
  // one of the external indices is in the emax model space
  // instead of the emax_3body model space. This implementation truncates
  // all indices to the emax_3body level, effectively performing a truncation
  // on the commutator level rather than on the matrix element level.
  // This implementation is fairly cheap and pretty efficient on big nodes with many threads.
//  void comm232ss_expand_reduced(const Operator &X, const Operator &Y, Operator &Z)
//  {
//    comm232::comm232ss_expand_impl_red(X, Y, Z);
//  }

//  //// Begin modifications to comm232ss by Matthias.
//  namespace
//  {
//    static inline size_t Hash_comm232_key2(const std::array<size_t, 5> &kljJJ)
//    {
//      return kljJJ[0] + (kljJJ[1] << 8) + (kljJJ[2] << 16) + (kljJJ[3] << 24) +
//             (kljJJ[4] << 32);
//    };
//
//    static inline std::array<size_t, 5> Unhash_comm232_key2(size_t hash)
//    {
//      std::array<size_t, 5> kljJJ;
//      size_t mask = 0xFF;
//      for (size_t i = 0; i < 5; i += 1)
//      {
//        kljJJ[i] = (hash >> i * 8) & mask;
//      }
//      return kljJJ;
//    };
//  } // namespace
//
//  /// These should probably be moved to some other place...
//  template <typename V>
//  size_t GetVectorSize(const std::vector<V> &v)
//  {
//    return v.size() * sizeof(V);
//  }
//
//  template <typename K, typename V>
//  size_t GetMapSize(const std::map<K, std::vector<V>> &m)
//  {
//    size_t x = 0;
//    x += m.size() * (sizeof(K) + sizeof(std::vector<V>));
//    for (const auto &y : m)
//    {
//      x += GetVectorSize(y.second);
//    }
//    return x;
//  }
//
//  template <typename K, typename V>
//  size_t GetMapSizeFlat(const std::unordered_map<K, V> &m)
//  {
//    size_t x = 0;
//    x += m.size() * (sizeof(K) + sizeof(V));
//    return x;
//  }
//
//  template <typename K>
//  size_t GetSetSizeFlat(const std::unordered_set<K> &m)
//  {
//    size_t x = 0;
//    x += m.size() * (sizeof(K));
//    return x;
//  }

  //*****************************************************************************************
  //
  //  |         |
  // i|        j|     Uncoupled expression:
  //  *~~[X]~*  |         Z_ijkl = -1/2 sum_abc (nanbn`c+n`an`bnc) ( (1-Pij) X_icab * Y_abjklc - (1-Pkl) Yijcabl * Xabkc )
  // a|    b/c\ |
  //  |    /   \|
  //  *~~[Y]~~~~*      Coupled expression:
  //  |   |              Z_{ijkl}^{J} = -1/2 sum_abc (nanbn`c+n`an`bnc) sum_J'J" (2J'+1)(2J"+1)/sqrt(2J+1) (-1)^{2J"+J'-J}
  // k|  l|                           *  [   (1 - (-1)^{i+j-J}Pij) (-1)^{j-c} { j  J" J' } X_icab^{J'} * Y_{abjklc}^{J'JJ"}
  //                                                                          { c  i  J  }
  //
  //                                       -(1 - (-1)^{k+l-J}Pkl)  (-1)^{l-c} { l  J" J' } Y_{ijcabl}^{JJ'J"} * X_{abkc}^{J'}  ]
  //                                                                          { c  k  J  }
  //
  //   The factor 1/2 out front is absorbed by the fact that we only do the ordering a<=b  (need to deal carefully with a==b)
  //
  //  For efficiency, we unfold this diagram to look like
  //  |
  // i|
  //  *~~[X]~*
  //  |a  b / \c
  //  |    /   \
//  *~~[Y]~~~*
  //  |   |    |
  //  |k  |l   |j
  //
  //  Verified with UnitTest
  //
  // This is the time hog of the n^7 scaling terms   (seems to be doing better...)
  // Now this is fine. The 223 commutator is taking all the time.
  //
  // This is a variant of comm232ss_srs_optimized that was broken up by Matthias
  // and tweaked to get improved performance when using many threads in large model spaces.
//  void comm232ss_mh_optimized(const Operator &X, const Operator &Y, Operator &Z)
//  {
//    size_t lookups_size = 0;
//    double tstart = omp_get_wtime();
//    //  auto& X2 = X.TwoBody;
//    auto &X3 = X.ThreeBody;
//    //  auto& Y2 = Y.TwoBody;
//    auto &Y3 = Y.ThreeBody;
//    auto &Z2 = Z.TwoBody;
//
//    bool x_has_3 = X3.IsAllocated();
//    bool y_has_3 = Y3.IsAllocated();
//
//    int hermX = X.IsHermitian() ? 1 : -1;
//    int hermY = Y.IsHermitian() ? 1 : -1;
//
//    Z.modelspace->PreCalculateSixJ(); // if this has already been done, this does nothing.
//
//    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();
//
//    // Hash for lookup
//    //  auto Hash_comm232_key = [] ( std::array<size_t,5>& kljJJ ){  return kljJJ[0] + (kljJJ[1] << 8) + (kljJJ[2] << 16) + (kljJJ[3] << 24) + (kljJJ[4]<<32);};
//
//    size_t nch = Z.modelspace->GetNumberTwoBodyChannels();
//
//    // first, enumerate the one-body channels |i> => ji,parityi,tzi
//    std::map<std::array<int, 3>, std::vector<size_t>> local_one_body_channels;          //  maps {j,parity,tz} => vector of index_j
//    std::map<std::array<int, 3>, std::vector<size_t>> external_local_one_body_channels; //  maps {j,parity,tz} => vector of index_j
//    comm232_new_Determine1BChannels(Z, Y, local_one_body_channels, external_local_one_body_channels);
//    lookups_size += GetMapSize(local_one_body_channels);
//    lookups_size += GetMapSize(external_local_one_body_channels);
//
//    // next, figure out which three-body states |klj`> and |abc`> exist, make a list, and give them an
//    // index for where they'll sit in the matrix
//    std::map<std::array<int, 3>, std::vector<std::array<size_t, 4>>> klj_list; //  maps {j,parity,tz} => {kljJ}
//    std::map<std::array<int, 3>, arma::mat> ZMAT_list;                         //  maps {j,parity,tz} => matrix <i|Z|klj`>
//
//    std::vector<std::array<int, 3>> obc_keys;
//    for (auto &iter_i : external_local_one_body_channels)
//      obc_keys.push_back(iter_i.first);
//    size_t nkeys = obc_keys.size();
//    lookups_size += GetVectorSize(obc_keys);
//
//    comm232_new_Populate1BChannel(Z, Y, e_fermi, local_one_body_channels, external_local_one_body_channels, obc_keys, klj_list, ZMAT_list);
//
//    lookups_size += GetMapSize(klj_list);
//
//    // This loop is what takes all the time.
//    // Outer loop over one-body quantum numbers, which also label the three-body pph states |klj`> and |abc`>
//    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
//    // #pragma omp parallel for schedule(dynamic,1)
//    for (size_t ikey = 0; ikey < nkeys; ikey++)
//    {
//      double tloopbody_start = omp_get_wtime();
//
//      const auto &obc_key = obc_keys[ikey]; // this is an array {j2,parity,tz2}
//      int j2i = obc_key[0];
//      const std::vector<size_t> &obc_orbits = external_local_one_body_channels[obc_key]; // the orbits that have quantum numbers {j2,parity,tz2}
//      const auto &klj_list_i = klj_list[obc_key];                                        // list of 3-body pph stats |klj`> with quantum numbers {j2,parity,tz2}
//
//      std::vector<size_t> abc_list; // keep track of which klj states should go in the abc loop
//      std::vector<double> abc_occ_list;
//      comm232_new_Determine3BStatesIn1BChannel(Y, Z, klj_list, obc_key, e_fermi, abc_list, abc_occ_list);
//      lookups_size += GetVectorSize(abc_list);
//      lookups_size += GetVectorSize(abc_occ_list);
//
//      size_t dim_i = obc_orbits.size();          // how many sp states are in this jpt channel
//      size_t dim_klj = klj_list[obc_key].size(); // how many 3-body pph states in this jpt channel
//      size_t dim_abc = abc_list.size();          // how many 3-body states which contribute to the |abc`> sum
//
//      double tmatalloc_start = omp_get_wtime();
//      // Now allocate the matrices in this channel
//      arma::mat X2MAT(dim_i, dim_abc, arma::fill::zeros);
//      arma::mat Y2MAT(dim_i, dim_abc, arma::fill::zeros);
//      arma::mat X3MAT(dim_abc, dim_klj, arma::fill::zeros);
//      arma::mat Y3MAT(dim_abc, dim_klj, arma::fill::zeros);
//      Z.profiler.timer["_" + std::string(__func__) + ", mat alloc"] += omp_get_wtime() - tmatalloc_start;
//
//      // figure out which recouplings we'll need when filling the matrices
//      std::unordered_set<size_t> kljJJ_needed; // we use a set to avoid repeats
//      comm232_new_GenerateRequiredRecouplings(Z, Y, abc_list, klj_list_i, e_fermi, dim_abc, dim_klj, kljJJ_needed);
//      lookups_size += GetSetSizeFlat(kljJJ_needed);
//
//      // now compute the couplings once and store them in a hash table
//      std::vector<double> recoupling_cache;                       // a vector storing all the relevant recoupling info
//      std::unordered_map<size_t, size_t> recoupling_cache_lookup; // a lookup table so we can find the info easily
//      comm232_new_ComputeRequiredRecouplings(Y, kljJJ_needed, recoupling_cache, recoupling_cache_lookup);
//      lookups_size += GetVectorSize(recoupling_cache);
//      lookups_size += GetMapSizeFlat(recoupling_cache_lookup);
//
//      // Now fill the matrices
//      comm232_new_FillMatrices(Z, X, Y, dim_abc, dim_i, dim_klj, j2i, x_has_3, y_has_3, abc_list, klj_list_i, e_fermi, abc_occ_list, obc_orbits, recoupling_cache, recoupling_cache_lookup, X2MAT, Y2MAT, X3MAT, Y3MAT);
//
//      double tmatmul_start = omp_get_wtime();
//      // now we do the mat mult
//      ZMAT_list[obc_key] = (X2MAT * Y3MAT - Y2MAT * X3MAT);
//      Z.profiler.timer["_" + std::string(__func__) + ", mat mul"] += omp_get_wtime() - tmatmul_start;
//      Z.profiler.timer["_" + std::string(__func__) + ", loop body"] += omp_get_wtime() - tloopbody_start;
//
//    } // for iter_i in local one body channels
//
//    std::cout << "Prestoring all lookups would require " << lookups_size / (1024.0 * 1024.0) << " MB\n";
//
//    // now we need to unpack all that mess and store it in Z
//    comm232_new_Unpack2BResult(X, Y, nch, external_local_one_body_channels, klj_list, ZMAT_list, Z);
//    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
//  }

//  void comm232_new_Determine1BChannels(
//      const Operator &Z,
//      const Operator &Y,
//      std::map<std::array<int, 3>, std::vector<size_t>> &local_one_body_channels,
//      std::map<std::array<int, 3>, std::vector<size_t>> &external_local_one_body_channels)
//  {
//    const auto tstart = omp_get_wtime();
//
//    const auto &Y3 = Y.ThreeBody;
//    for (auto j : Z.modelspace->all_orbits)
//    {
//      Orbit &oj = Z.modelspace->GetOrbit(j);
//      // if (!Y3.IsOrbitIn3BodyEMaxTruncation(oj)) continue;
//      std::array<int, 3> obc = {oj.j2, oj.l % 2, oj.tz2};
//      if (local_one_body_channels.find(obc) == local_one_body_channels.end())
//        local_one_body_channels[obc] = {j};
//      else
//        local_one_body_channels[obc].push_back(j);
//      if (imsrg3_valence_2b and (oj.cvq != 1))
//        continue;
//      if (external_local_one_body_channels.find(obc) == external_local_one_body_channels.end())
//        external_local_one_body_channels[obc] = {j};
//      else
//        external_local_one_body_channels[obc].push_back(j);
//    }
//    Y.profiler.timer[__func__] += omp_get_wtime() - tstart;
//  }
//
//  void comm232_new_Populate1BChannel(
//      const Operator &Z,
//      const Operator &Y,
//      const std::map<int, double> &e_fermi,
//      const std::map<std::array<int, 3>, std::vector<size_t>> &local_one_body_channels,
//      const std::map<std::array<int, 3>, std::vector<size_t>> &external_local_one_body_channels,
//      const std::vector<std::array<int, 3>> &obc_keys,
//      std::map<std::array<int, 3>, std::vector<std::array<size_t, 4>>> &klj_list,
//      std::map<std::array<int, 3>, arma::mat> &ZMAT_list)
//  {
//    const auto tstart = omp_get_wtime();
//
//    const auto &Y3 = Y.ThreeBody;
//    size_t nkeys = obc_keys.size();
//    for (size_t ikey = 0; ikey < nkeys; ikey++)
//    {
//      auto &obc_key = obc_keys[ikey];
//      const std::vector<size_t> &obc_orbits = external_local_one_body_channels.at(obc_key);
//      int j2i = obc_key[0];
//      int parityi = obc_key[1];
//      int tz2i = obc_key[2];
//      klj_list[obc_key] = {};
//      auto &klj_list_i = klj_list[obc_key];
//
//      for (auto &iter_j : local_one_body_channels)
//      {
//        int j2j = iter_j.first[0];
//        int parityj = iter_j.first[1];
//        int tz2j = iter_j.first[2];
//
//        int Jkl_min = std::abs(j2i - j2j) / 2;
//        int Jkl_max = (j2i + j2j) / 2;
//        int Tzkl = (tz2i + tz2j) / 2;
//        int paritykl = (parityi + parityj) % 2;
//        for (int Jkl = Jkl_min; Jkl <= Jkl_max; Jkl++)
//        {
//
//          size_t ch_kl = Z.modelspace->GetTwoBodyChannelIndex(Jkl, paritykl, Tzkl);
//          TwoBodyChannel &tbc_kl = Z.modelspace->GetTwoBodyChannel(ch_kl);
//          size_t nkets_kl = tbc_kl.GetNumberKets();
//          for (size_t iket_kl = 0; iket_kl < nkets_kl; iket_kl++)
//          {
//            Ket &ket_kl = tbc_kl.GetKet(iket_kl);
//            if (!Y3.IsOrbitIn3BodyEMaxTruncation(ket_kl.p))
//              continue;
//            if (!Y3.IsOrbitIn3BodyEMaxTruncation(ket_kl.q))
//              continue;
//            int e_k = 2 * ket_kl.op->n + ket_kl.op->l;
//            int e_l = 2 * ket_kl.oq->n + ket_kl.oq->l;
//            // if ( (e_k + e_l) > Z.modelspace->GetE3max()) continue;
//            double de_k = std::abs(e_k - e_fermi.at(ket_kl.op->tz2));
//            double de_l = std::abs(e_l - e_fermi.at(ket_kl.oq->tz2));
//            if ((de_k + de_l) > Z.modelspace->GetdE3max())
//              continue;
//            // double occnat_k = ket_kl.op->occ_nat;
//            // double occnat_l = ket_kl.oq->occ_nat;
//            for (size_t j : iter_j.second)
//            {
//              // if (!Y3.IsOrbitIn3BodyEMaxTruncation(j)) continue;
//              // if (!Z.ThreeBody.IsKetInEMaxTruncations(ket_kl.p, ket_kl.q, j)) continue;
//              Orbit &oj = Z.modelspace->GetOrbit(j);
//              if (!Y3.IsOrbitIn3BodyEMaxTruncation(oj))
//                continue;
//              // double occnat_j = oj.occ_nat;
//              // if ( (occnat_k*(1-occnat_k) * occnat_l*(1-occnat_l) * occnat_j*(1-occnat_j) ) < Z.modelspace->GetOccNat3Cut() ) continue;
//
//              klj_list_i.push_back({ket_kl.p, ket_kl.q, j, (size_t)tbc_kl.J});
//
//            } // for j
//          }   // for iket_kl
//        }     // for Jkl
//      }       // for iter_j in local one body channels
//
//      size_t dim_i = obc_orbits.size(); // how many sp states are in this jpt channel
//      size_t dim_klj = klj_list[obc_key].size();
//      ZMAT_list[obc_key] = arma::mat(dim_i, dim_klj);
//    } // for ikey
//    Y.profiler.timer[__func__] += omp_get_wtime() - tstart;
//  }
//
//  void comm232_new_Determine3BStatesIn1BChannel(
//      const Operator &Y,
//      const Operator &Z,
//      const std::map<std::array<int, 3>, std::vector<std::array<size_t, 4>>> &klj_list,
//      std::array<int, 3> obc_key,
//      const std::map<int, double> &e_fermi,
//      std::vector<size_t> &abc_list,
//      std::vector<double> &abc_occ_list)
//  {
//    const auto tstart = omp_get_wtime();
//
//    const auto &Y3 = Y.ThreeBody;
//    for (size_t i_kljJ = 0; i_kljJ < klj_list.at(obc_key).size(); i_kljJ++)
//    {
//      auto &kljJ = klj_list.at(obc_key)[i_kljJ];
//      size_t a = kljJ[0];
//      size_t b = kljJ[1];
//      size_t c = kljJ[2];
//      Orbit &oa = Z.modelspace->GetOrbit(a);
//      Orbit &ob = Z.modelspace->GetOrbit(b);
//      Orbit &oc = Z.modelspace->GetOrbit(c);
//      // if (!Y3.IsOrbitIn3BodyEMaxTruncation(oa)) continue;
//      // if (!Y3.IsOrbitIn3BodyEMaxTruncation(ob)) continue;
//      // if (!Y3.IsOrbitIn3BodyEMaxTruncation(oc)) continue;
//      double e_a = 2 * oa.n + oa.l;
//      double e_b = 2 * ob.n + ob.l;
//      double de_a = std::abs(e_a - e_fermi.at(oa.tz2));
//      double de_b = std::abs(e_b - e_fermi.at(ob.tz2));
//      if ((e_a + e_b) > Z.modelspace->GetE3max())
//        continue;
//      if ((de_a + de_b) > Z.modelspace->GetdE3max())
//        continue;
//      double na = oa.occ;
//      double nb = ob.occ;
//      double nc = oc.occ;
//      double occupation_factor = na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc;
//      if (std::abs(occupation_factor) < 1e-6)
//        continue;
//      if (a == b)
//        occupation_factor *= 0.5; // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b
//      abc_list.push_back(i_kljJ);
//      abc_occ_list.push_back(occupation_factor);
//    }
//    Y.profiler.timer[__func__] += omp_get_wtime() - tstart;
//  }
//
//  void comm232_new_GenerateRequiredRecouplings(
//      const Operator &Z,
//      const Operator &Y,
//      const std::vector<size_t> &abc_list,
//      const std::vector<std::array<size_t, 4>> &klj_list_i,
//      const std::map<int, double> &e_fermi,
//      size_t dim_abc,
//      size_t dim_klj,
//      std::unordered_set<size_t> &kljJJ_needed)
//  {
//    const auto tstart_parallel = omp_get_wtime();
//
//    const auto &Y3 = Y.ThreeBody;
//    // Set up one set per thread
//    int num_threads = omp_get_max_threads();
//    std::vector<std::unordered_set<size_t>> kljJJ_needed_threadsafe_vec(num_threads);
//
//// Loop in parallel over many indices
//#pragma omp parallel for collapse(2)
//    for (size_t ind_abc = 0; ind_abc < dim_abc; ind_abc++)
//    {
//      for (size_t ind_klj = 0; ind_klj < dim_klj; ind_klj++)
//      {
//        // Get thread id for current thread
//        int thread_id = omp_get_thread_num();
//
//        // for abc loop
//        auto &abcJ = klj_list_i[abc_list[ind_abc]];
//        size_t a = abcJ[0];
//        size_t b = abcJ[1];
//        size_t c = abcJ[2];
//        int Jab = abcJ[3];
//        Orbit &oa = Z.modelspace->GetOrbit(a);
//        Orbit &ob = Z.modelspace->GetOrbit(b);
//        Orbit &oc = Z.modelspace->GetOrbit(c);
//        double e_a = 2 * oa.n + oa.l;
//        double e_b = 2 * ob.n + ob.l;
//        double e_c = 2 * oc.n + oc.l;
//        double de_a = std::abs(e_a - e_fermi.at(oa.tz2));
//        double de_b = std::abs(e_b - e_fermi.at(ob.tz2));
//        double de_c = std::abs(e_c - e_fermi.at(oc.tz2));
//        double occnat_a = oa.occ_nat;
//        double occnat_b = ob.occ_nat;
//        double occnat_c = oc.occ_nat;
//        int j2c = oc.j2;
//
//        // for klj loop
//        auto &kljJ = klj_list_i[ind_klj];
//        size_t k = kljJ[0];
//        size_t l = kljJ[1];
//        size_t j = kljJ[2];
//        int Jkl = kljJ[3];
//        Orbit &oj = Z.modelspace->GetOrbit(j);
//        Orbit &ok = Z.modelspace->GetOrbit(k);
//        Orbit &ol = Z.modelspace->GetOrbit(l);
//
//        double e_j = 2 * oj.n + oj.l;
//        double e_k = 2 * ok.n + ok.l;
//        double e_l = 2 * ol.n + ol.l;
//        double de_j = std::abs(e_j - e_fermi.at(oj.tz2));
//        double de_k = std::abs(e_k - e_fermi.at(ok.tz2));
//        double de_l = std::abs(e_l - e_fermi.at(ol.tz2));
//
//        double occnat_j = oj.occ_nat;
//        double occnat_k = ok.occ_nat;
//        double occnat_l = ol.occ_nat;
//        if ((e_k + e_l + e_c) > Z.modelspace->GetE3max())
//          continue;
//        if ((e_a + e_b + e_j) > Z.modelspace->GetE3max())
//          continue;
//        if ((de_k + de_l + de_c) > Z.modelspace->GetdE3max())
//          continue;
//        if ((de_a + de_b + de_j) > Z.modelspace->GetdE3max())
//          continue;
//        if ((occnat_k * (1 - occnat_k) * occnat_l * (1 - occnat_l) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
//          continue;
//        if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) < Z.modelspace->GetOccNat3Cut())
//          continue;
//        if (imsrg3_no_qqq and ((ok.cvq + ol.cvq + oc.cvq) > 5 or (oa.cvq + ob.cvq + oj.cvq) > 5))
//          continue;
//        if (imsrg3_valence_2b and (ok.cvq != 1 or ol.cvq != 1 or oj.cvq != 1))
//          continue;
//        int twoJp_min = std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * Jkl - j2c));
//        int twoJp_max = std::min(2 * Jab + oj.j2, 2 * Jkl + j2c);
//        for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
//        {
//          // Write to set corresponding to local thread
//          kljJJ_needed_threadsafe_vec[thread_id].insert(Hash_comm232_key2({k, l, c, (size_t)Jkl, (size_t)twoJp}));
//          kljJJ_needed_threadsafe_vec[thread_id].insert(Hash_comm232_key2({a, b, j, (size_t)Jab, (size_t)twoJp}));
//        } // for twoJp
//      }   // for ind_klj
//    }     // for ind_abc
//
//    Y.profiler.timer[std::string(__func__) + " parallel"] += omp_get_wtime() - tstart_parallel;
//
//    const auto tstart_merge = omp_get_wtime();
//
//    // Merge sets
//    for (const auto &thread_safe_kljJJ : kljJJ_needed_threadsafe_vec)
//    {
//      for (const auto &el : thread_safe_kljJJ)
//      {
//        kljJJ_needed.insert(el);
//      }
//    }
//
//    Y.profiler.timer[std::string(__func__) + " merge"] += omp_get_wtime() - tstart_merge;
//  }
//
//  void comm232_new_ComputeRequiredRecouplings(
//      const Operator &Y,
//      const std::unordered_set<size_t> &kljJJ_needed,
//      std::vector<double> &recoupling_cache,
//      std::unordered_map<size_t, size_t> &recoupling_cache_lookup)
//  {
//    const auto tstart_parallel = omp_get_wtime();
//
//    // Read kljJJ_needed into vector so we can use OpenMP pragma on for loop
//    std::vector<size_t> kljJJ_needed_vec(kljJJ_needed.cbegin(), kljJJ_needed.cend());
//
//    // Create thread-safe recoupling_cache and recoupling_cache_lookup, which we need to reduce later
//    int num_threads = omp_get_max_threads();
//
//    // Define chunk_size for dynamic scheduling
//    // Giving each thread 0.1% of the problem per call to the OpenMP runtime seems reasonable.
//    // On most clusters, this means threads will make somewhere between 5 and 20 calls to the runtime.
//    int chunk_size = std::max(kljJJ_needed_vec.size() / 1000, 1UL);
//
//    std::vector<std::vector<double>> recoupling_cache_threadsafe(num_threads);
//    std::vector<std::unordered_map<size_t, size_t>> recoupling_cache_lookup_threadsafe(num_threads);
//
//#pragma omp parallel for schedule(dynamic, chunk_size)
//    for (size_t index = 0; index < kljJJ_needed_vec.size(); index += 1)
//    {
//      size_t hash = kljJJ_needed_vec[index];
//      int thread_id = omp_get_thread_num();
//      const auto kljJJ = Unhash_comm232_key2(hash);
//      size_t k = kljJJ[0];
//      size_t l = kljJJ[1];
//      size_t j = kljJJ[2];
//      size_t Jkl = kljJJ[3];
//      size_t twoJp = kljJJ[4];
//
//      if (!Y.ThreeBody.IsKetValid(Jkl, twoJp, k, l, j))
//        continue;
//      std::vector<size_t> iketlist;
//      std::vector<double> recouplelist;
//      size_t ch_check = Y.ThreeBody.GetKetIndex_withRecoupling(Jkl, twoJp, k, l, j, iketlist, recouplelist);
//      size_t listsize = iketlist.size();
//      recoupling_cache_lookup_threadsafe[thread_id][hash] = recoupling_cache_threadsafe[thread_id].size();
//      recoupling_cache_threadsafe[thread_id].push_back(ch_check);
//      recoupling_cache_threadsafe[thread_id].push_back(listsize);
//      for (size_t ik : iketlist)
//        recoupling_cache_threadsafe[thread_id].push_back(ik); // Store the size_t (an unsigned int) as a double. I hope this doesn't cause problems...
//      for (double rec : recouplelist)
//        recoupling_cache_threadsafe[thread_id].push_back(rec); // Point to where this element lives in the vector
//    }
//    Y.profiler.timer[std::string(__func__) + " parallel"] += omp_get_wtime() - tstart_parallel;
//
//    const auto tstart_merge = omp_get_wtime();
//    // Copying once we've reserved the right size should be cheap relative to recoupling logic,
//    // so this can be done in serial
//    size_t total_cache_size = 0;
//    size_t total_cache_lookup_size = 0;
//    for (int thread_id = 0; thread_id < num_threads; thread_id += 1)
//    {
//      total_cache_size += recoupling_cache_threadsafe[thread_id].size();
//      total_cache_lookup_size += recoupling_cache_lookup_threadsafe[thread_id].size();
//    }
//
//    recoupling_cache.reserve(total_cache_size);
//    recoupling_cache_lookup.reserve(total_cache_lookup_size);
//    for (int thread_id = 0; thread_id < num_threads; thread_id += 1)
//    {
//      for (const auto &kv_pair : recoupling_cache_lookup_threadsafe[thread_id])
//      {
//        const size_t hash = kv_pair.first;
//        const size_t index = kv_pair.second;
//
//        const auto ch_check = static_cast<size_t>(recoupling_cache_threadsafe[thread_id][index]);
//        const auto list_size = static_cast<size_t>(recoupling_cache_threadsafe[thread_id][index + 1]);
//        recoupling_cache_lookup[hash] = recoupling_cache.size();
//        recoupling_cache.push_back(ch_check);
//        recoupling_cache.push_back(list_size);
//
//        // Copy iketlist
//        for (size_t offset = 0; offset < list_size; offset += 1)
//        {
//          recoupling_cache.push_back(recoupling_cache_threadsafe[thread_id][index + 2 + offset]);
//        }
//
//        // Copy recouple list
//        for (size_t offset = 0; offset < list_size; offset += 1)
//        {
//          recoupling_cache.push_back(recoupling_cache_threadsafe[thread_id][index + 2 + offset + list_size]);
//        }
//      }
//    }
//    Y.profiler.timer[std::string(__func__) + " merge"] += omp_get_wtime() - tstart_merge;
//  }
//
//  void comm232_new_FillMatrices(const Operator &Z, const Operator &X, const Operator Y, size_t dim_abc, size_t dim_i, size_t dim_klj,
//                                int j2i,
//                                bool x_has_3,
//                                bool y_has_3,
//                                const std::vector<size_t> &abc_list,
//                                const std::vector<std::array<size_t, 4>> &klj_list_i,
//                                const std::map<int, double> &e_fermi,
//                                const std::vector<double> &abc_occ_list,
//                                const std::vector<size_t> &obc_orbits,
//                                const std::vector<double> &recoupling_cache,
//                                const std::unordered_map<size_t, size_t> &recoupling_cache_lookup,
//                                arma::mat &X2MAT,
//                                arma::mat &Y2MAT,
//                                arma::mat &X3MAT,
//                                arma::mat &Y3MAT)
//  {
//    const auto tstart = omp_get_wtime();
//
//    const auto &Y3 = Y.ThreeBody;
//    const auto &X3 = X.ThreeBody;
//
//// With a bit more work, we can 2 loops for the 3-body part
//#pragma omp parallel for
//    for (size_t ind_abc = 0; ind_abc < dim_abc; ind_abc++) // rows of X2MAT
//    {
//
//      auto &abcJ = klj_list_i[abc_list[ind_abc]];
//      size_t a = abcJ[0];
//      size_t b = abcJ[1];
//      size_t c = abcJ[2];
//      int Jab = abcJ[3];
//      Orbit &oa = Z.modelspace->GetOrbit(a);
//      Orbit &ob = Z.modelspace->GetOrbit(b);
//      Orbit &oc = Z.modelspace->GetOrbit(c);
//      // if (!Y3.IsOrbitIn3BodyEMaxTruncation(oa)) continue;
//      // if (!Y3.IsOrbitIn3BodyEMaxTruncation(ob)) continue;
//      // if (!Y3.IsOrbitIn3BodyEMaxTruncation(oc)) continue;
//      int j2c = oc.j2;
//      double jc = 0.5 * j2c;
//      double occ_abc = abc_occ_list[ind_abc]; // TODO: we can probably incorporate that hat factor with the occupation factor
//
//      double de_a = std::abs(2 * oa.n + oa.l - e_fermi.at(oa.tz2));
//      double de_b = std::abs(2 * ob.n + ob.l - e_fermi.at(ob.tz2));
//      double de_c = std::abs(2 * oc.n + oc.l - e_fermi.at(oc.tz2));
//      double occnat_a = oa.occ_nat;
//      double occnat_b = ob.occ_nat;
//      double occnat_c = oc.occ_nat;
//
//      // fill the 2 body part
//      for (size_t ind_i = 0; ind_i < dim_i; ind_i++) // ind_i indexes the list of sp states in this one body channel. columns of X2MAT
//      {
//        size_t i = obc_orbits[ind_i]; // i is the orbit index as per ModelSpace
//
//        X2MAT(ind_i, ind_abc) = -sqrt((2 * Jab + 1.)) * occ_abc * X.TwoBody.GetTBME_J(Jab, c, i, a, b);
//        Y2MAT(ind_i, ind_abc) = -sqrt((2 * Jab + 1.)) * occ_abc * Y.TwoBody.GetTBME_J(Jab, c, i, a, b);
//      } // for ind_i
//    }
//
//    size_t max_num_threads = omp_get_max_threads();
//    // Each thread makes on average 10 calls to the runtime per commutator.
//    // This is fairly low while allowing the guided scheduling to do its job.
//    size_t chunk_size = std::max((dim_abc * dim_klj) / max_num_threads / 10, 1UL);
//// now the 3 body part
//#pragma omp parallel for schedule(guided, chunk_size) collapse(2)
//    for (size_t ind_abc = 0; ind_abc < dim_abc; ind_abc++) // rows of X2MAT
//    {
//      for (size_t ind_klj = 0; ind_klj < dim_klj; ind_klj++)
//      {
//        // For abc
//        auto &abcJ = klj_list_i[abc_list[ind_abc]];
//        size_t a = abcJ[0];
//        size_t b = abcJ[1];
//        size_t c = abcJ[2];
//        int Jab = abcJ[3];
//        Orbit &oa = Z.modelspace->GetOrbit(a);
//        Orbit &ob = Z.modelspace->GetOrbit(b);
//        Orbit &oc = Z.modelspace->GetOrbit(c);
//        // MH: Would be great if these checks were redundant.
//        // if (!Y3.IsOrbitIn3BodyEMaxTruncation(oa)) continue;
//        // if (!Y3.IsOrbitIn3BodyEMaxTruncation(ob)) continue;
//        // if (!Y3.IsOrbitIn3BodyEMaxTruncation(oc)) continue;
//        int j2c = oc.j2;
//        double jc = 0.5 * j2c;
//        double occ_abc = abc_occ_list[ind_abc]; // TODO: we can probably incorporate that hat factor with the occupation factor
//
//        double de_a = std::abs(2 * oa.n + oa.l - e_fermi.at(oa.tz2));
//        double de_b = std::abs(2 * ob.n + ob.l - e_fermi.at(ob.tz2));
//        double de_c = std::abs(2 * oc.n + oc.l - e_fermi.at(oc.tz2));
//        double occnat_a = oa.occ_nat;
//        double occnat_b = ob.occ_nat;
//        double occnat_c = oc.occ_nat;
//
//        // For ilj
//        // <abi Jab twoJ | X | klc Jkl twoJ >
//        auto &kljJ = klj_list_i[ind_klj];
//        size_t k = kljJ[0];
//        size_t l = kljJ[1];
//        size_t j = kljJ[2];
//        int Jkl = kljJ[3];
//        Orbit &oj = Z.modelspace->GetOrbit(j);
//        Orbit &ok = Z.modelspace->GetOrbit(k);
//        Orbit &ol = Z.modelspace->GetOrbit(l);
//        // MH: Would be great if these checks were redundant.
//        // if (!Y3.IsOrbitIn3BodyEMaxTruncation(oj)) continue;
//        // if (!Y3.IsOrbitIn3BodyEMaxTruncation(ok)) continue;
//        // if (!Y3.IsOrbitIn3BodyEMaxTruncation(ol)) continue;
//        double ji = 0.5 * j2i;
//        double jj = 0.5 * oj.j2;
//        //        double jk = 0.5 * ok.j2;
//        //        double jl = 0.5 * ol.j2;
//
//        double de_j = std::abs(2 * oj.n + oj.l - e_fermi.at(oj.tz2));
//        double de_k = std::abs(2 * ok.n + ok.l - e_fermi.at(ok.tz2));
//        double de_l = std::abs(2 * ol.n + ol.l - e_fermi.at(ol.tz2));
//        double occnat_j = oj.occ_nat;
//        double occnat_k = ok.occ_nat;
//        double occnat_l = ol.occ_nat;
//        // Are these checks redundant?
//        // MH: Would be great if these checks were redundant.
//        if ((de_a + de_b + de_j) > Z.modelspace->GetdE3max())
//          continue;
//        if ((de_k + de_l + de_c) > Z.modelspace->GetdE3max())
//          continue;
//        if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) < Z.modelspace->GetOccNat3Cut())
//          continue;
//        if ((occnat_k * (1 - occnat_k) * occnat_l * (1 - occnat_l) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
//          continue;
//        if (imsrg3_no_qqq and ((ok.cvq + ol.cvq + oc.cvq) > 5 or (oa.cvq + ob.cvq + oj.cvq) > 5))
//          continue;
//        if (imsrg3_valence_2b and (ok.cvq != 1 or ol.cvq != 1 or oj.cvq != 1))
//          continue;
//
//        int twoJp_min = std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * Jkl - j2c));
//        int twoJp_max = std::min(2 * Jab + oj.j2, 2 * Jkl + j2c);
//
//        for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
//        {
//          if (!Y.ThreeBody.IsKetValid(Jkl, twoJp, k, l, c))
//            continue;
//          if (!Y.ThreeBody.IsKetValid(Jab, twoJp, a, b, j))
//            continue;
//          //           double sixj = X.modelspace->GetSixJ(Jkl,jj,ji,  Jab, jc, 0.5*twoJp );
//          double sixj = X.modelspace->GetSixJ(jj, ji, Jkl, jc, 0.5 * twoJp, Jab);
//
//          std::array<size_t, 5> abjJJ = {a, b, j, (size_t)Jab, (size_t)twoJp};
//          std::array<size_t, 5> klcJJ = {k, l, c, (size_t)Jkl, (size_t)twoJp};
//          size_t hash_abj = Hash_comm232_key2(abjJJ);
//          size_t hash_klc = Hash_comm232_key2(klcJJ);
//          size_t pointer_abj = recoupling_cache_lookup.at(hash_abj);
//          size_t pointer_klc = recoupling_cache_lookup.at(hash_klc);
//          size_t ch_check = recoupling_cache[pointer_abj];
//          size_t listsize_abj = recoupling_cache[pointer_abj + 1];
//          size_t listsize_klc = recoupling_cache[pointer_klc + 1];
//
//          double xabjklc = 0;
//          double yabjklc = 0;
//          for (size_t ilist_abj = 0; ilist_abj < listsize_abj; ilist_abj++)
//          {
//            size_t ibra_abj = (size_t)recoupling_cache[pointer_abj + 2 + ilist_abj];
//            double recouple_bra = recoupling_cache[pointer_abj + 2 + listsize_abj + ilist_abj];
//            for (size_t ilist_klc = 0; ilist_klc < listsize_klc; ilist_klc++)
//            {
//              size_t iket_klc = (size_t)recoupling_cache[pointer_klc + 2 + ilist_klc];
//              double recouple_ket = recoupling_cache[pointer_klc + 2 + listsize_klc + ilist_klc];
//
//              //               xabjklc += recouple_bra*recouple_ket * X3.GetME_pn_ch(ch_check,ch_check, ibra_abj, iket_klc );
//              if (x_has_3)
//                xabjklc += recouple_bra * recouple_ket * X3.GetME_pn_ch(ch_check, ch_check, ibra_abj, iket_klc);
//              if (y_has_3)
//                yabjklc += recouple_bra * recouple_ket * Y3.GetME_pn_ch(ch_check, ch_check, ibra_abj, iket_klc);
//            }
//          }
//
//          X3MAT(ind_abc, ind_klj) += (twoJp + 1.) / sqrt(2 * Jkl + 1.) * sixj * xabjklc;
//          Y3MAT(ind_abc, ind_klj) += (twoJp + 1.) / sqrt(2 * Jkl + 1.) * sixj * yabjklc;
//
//        } // for twoJp
//      }   // for ind_klj
//    }     // for ind_abc
//    Y.profiler.timer[__func__] += omp_get_wtime() - tstart;
//  }
//
//  void comm232_new_Unpack2BResult(
//      const Operator &X,
//      const Operator &Y,
//      size_t nch,
//      const std::map<std::array<int, 3>, std::vector<size_t>> &external_local_one_body_channels,
//      const std::map<std::array<int, 3>, std::vector<std::array<size_t, 4>>> &klj_list,
//      const std::map<std::array<int, 3>, arma::mat> &ZMAT_list,
//      Operator &Z)
//  {
//    const auto tstart = omp_get_wtime();
//
//    const auto &Y3 = Y.ThreeBody;
//    int hermX = X.IsHermitian() ? 1 : -1;
//    int hermY = Y.IsHermitian() ? 1 : -1;
//    auto &Z2 = Z.TwoBody;
//    for (size_t ch = 0; ch < nch; ch++)
//    {
//      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
//      size_t nkets = tbc.GetNumberKets();
//      size_t J = tbc.J;
//#pragma omp parallel for schedule(static) collapse(2)
//      for (size_t ibra = 0; ibra < nkets; ibra++)
//      {
//        for (size_t iket = 0; iket < nkets; iket++)
//        {
//          Ket &bra = tbc.GetKet(ibra);
//          size_t i = bra.p;
//          size_t j = bra.q;
//          Orbit &oi = Z.modelspace->GetOrbit(i);
//          Orbit &oj = Z.modelspace->GetOrbit(j);
//          if (imsrg3_valence_2b and (oi.cvq != 1 or oj.cvq != 1))
//            continue;
//
//          //      auto& lobc_i = local_one_body_channels.at({oi.j2,oi.l%2,oi.tz2});
//          //      auto& lobc_j = local_one_body_channels.at({oj.j2,oj.l%2,oj.tz2});
//          auto &lobc_i = external_local_one_body_channels.at({oi.j2, oi.l % 2, oi.tz2});
//          auto &lobc_j = external_local_one_body_channels.at({oj.j2, oj.l % 2, oj.tz2});
//          size_t ind_i = std::distance(lobc_i.begin(), std::find(lobc_i.begin(), lobc_i.end(), i));
//          size_t ind_j = std::distance(lobc_j.begin(), std::find(lobc_j.begin(), lobc_j.end(), j));
//
//          auto &list_i = klj_list.at({oi.j2, oi.l % 2, oi.tz2});
//          auto &list_j = klj_list.at({oj.j2, oj.l % 2, oj.tz2});
//          auto &ZMat_i = ZMAT_list.at({oi.j2, oi.l % 2, oi.tz2});
//          auto &ZMat_j = ZMAT_list.at({oj.j2, oj.l % 2, oj.tz2});
//
//          Ket &ket = tbc.GetKet(iket);
//          size_t k = ket.p;
//          size_t l = ket.q;
//          Orbit &ok = Z.modelspace->GetOrbit(k);
//          Orbit &ol = Z.modelspace->GetOrbit(l);
//          int phase_ij = X.modelspace->phase((oi.j2 + oj.j2) / 2 - J);
//          int phase_kl = X.modelspace->phase((ok.j2 + ol.j2) / 2 - J);
//
//          if (imsrg3_valence_2b and (ok.cvq != 1 or ol.cvq != 1))
//            continue;
//
//          //        auto& lobc_k = local_one_body_channels.at({ok.j2,ok.l%2,ok.tz2});
//          //        auto& lobc_l = local_one_body_channels.at({ol.j2,ol.l%2,ol.tz2});
//          auto &lobc_k = external_local_one_body_channels.at({ok.j2, ok.l % 2, ok.tz2});
//          auto &lobc_l = external_local_one_body_channels.at({ol.j2, ol.l % 2, ol.tz2});
//          size_t ind_k = std::distance(lobc_k.begin(), std::find(lobc_k.begin(), lobc_k.end(), k));
//          size_t ind_l = std::distance(lobc_l.begin(), std::find(lobc_l.begin(), lobc_l.end(), l));
//
//          auto &list_k = klj_list.at({ok.j2, ok.l % 2, ok.tz2});
//          auto &list_l = klj_list.at({ol.j2, ol.l % 2, ol.tz2});
//          auto &ZMat_k = ZMAT_list.at({ok.j2, ok.l % 2, ok.tz2});
//          auto &ZMat_l = ZMAT_list.at({ol.j2, ol.l % 2, ol.tz2});
//
//          std::array<size_t, 4> key_klj = {k, l, j, J};
//          std::array<size_t, 4> key_kli = {k, l, i, J};
//          std::array<size_t, 4> key_ijk = {i, j, k, J};
//          std::array<size_t, 4> key_ijl = {i, j, l, J};
//
//          size_t ind_klj = std::distance(list_i.begin(), std::find(list_i.begin(), list_i.end(), key_klj));
//          size_t ind_kli = std::distance(list_j.begin(), std::find(list_j.begin(), list_j.end(), key_kli));
//          size_t ind_ijl = std::distance(list_k.begin(), std::find(list_k.begin(), list_k.end(), key_ijl));
//          size_t ind_ijk = std::distance(list_l.begin(), std::find(list_l.begin(), list_l.end(), key_ijk));
//          double Z_jkli = 0;
//          double Z_iklj = 0;
//          double Z_lijk = 0;
//          double Z_kijl = 0;
//          if (ind_klj < list_i.size())
//            Z_iklj = ZMat_i(ind_i, ind_klj);
//          if (ind_kli < list_j.size())
//            Z_jkli = ZMat_j(ind_j, ind_kli);
//          if (ind_ijl < list_k.size())
//            Z_kijl = ZMat_k(ind_k, ind_ijl);
//          if (ind_ijk < list_l.size())
//            Z_lijk = ZMat_l(ind_l, ind_ijk);
//
//          double zijkl = phase_ij * Z_iklj - Z_jkli + hermX * hermY * (Z_lijk - phase_kl * Z_kijl);
//
//          // normalize the tbme
//          zijkl /= sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
//          //        zijkl *= -1.0 / sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
//          Z2.AddToTBMENonHerm(ch, ch, ibra, iket, zijkl);
//        } // for iket
//      }   // for ibra
//    }     // for ch
//    Y.profiler.timer[__func__] += omp_get_wtime() - tstart;
//  }
//
//  //// END additions to comm232ss by Matthias

  void comm232ss_debug(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    //  auto& X2 = X.TwoBody;
    auto &X3 = X.ThreeBody;
    //  auto& Y2 = Y.TwoBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    bool x_has_3 = X3.IsAllocated();
    bool y_has_3 = Y3.IsAllocated();
    //  bool x_has_3 = X3.is_allocated;
    //  bool y_has_3 = Y3.is_allocated;

    int hermX = X.IsHermitian() ? 1 : -1;
    int hermY = Y.IsHermitian() ? 1 : -1;

    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    int nch = Z.modelspace->GetNumberTwoBodyChannels();

    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (int ch = 0; ch < nch; ch++)
    {
      auto &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nkets = tbc.GetNumberKets();

      //    double tstart = omp_get_wtime();
      // The strategy used here is the following. We reorganize <ic|X|ab> as <i|X|abc'>  and  <abj|Y|klc> as <abc'|Y|klj'>
      // where the prime indicates time reversal. Then we can cast things as a matrix multiplication
      //  <i|Z|klj'> = <i|X|abc'><abc'|Y|klj'>   and then we need to transform Z back to <ij|Z|kl>
      //
      //   i|   c|      i|                a| b| j|       a| b|   /c'
      //    |____|  ==>  |____             |__|__|  ==>   |__|__/
      //    |    |       |    |\           |  |  |        |  |  \
    //   a|   b|      a|   b| \c'       k| l| c|       k| l|   \j'
      //
      // The matrices for X and Z will clearly not be square. The left side is just a single particle orbit, while the right has 3 orbits.
      // Here, I determine which single orbits will be needed in this two-body channel.
      std::set<size_t> ij_orbits_set;
      std::set<int> j_jvals_set;
      for (int ibra = 0; ibra < nkets; ibra++)
      {
        Ket &bra = tbc.GetKet(ibra);
        ij_orbits_set.insert(bra.p);
        ij_orbits_set.insert(bra.q);
        j_jvals_set.insert(bra.op->j2);
        j_jvals_set.insert(bra.oq->j2);
      }
      std::vector<size_t> ij_orbits;
      std::vector<int> j_jvals;
      std::map<size_t, size_t> ij_orbits_lookup;
      std::map<int, size_t> jvals_lookup;
      for (auto i : ij_orbits_set)
      {
        ij_orbits.push_back(i); // list of all the indices an orbit (either i or j ) in bra could have
        ij_orbits_lookup[i] = ij_orbits.size() - 1;
      }
      for (auto j : j_jvals_set)
      {
        j_jvals.push_back(j); // list of possible j values (angular momentum) an orbit in bra could have
        jvals_lookup[j] = j_jvals.size() - 1;
      }
      size_t nij = ij_orbits.size();
      size_t njvals = j_jvals.size();
      //    std::cout << "   ch = " << ch << "   njvals = " << njvals << std::endl;

      Z.profiler.timer["_comm232_block1"] += omp_get_wtime() - tstart;
      tstart = omp_get_wtime();

      // next, we make a list of the abc' combinations that will inter into the sum
      std::vector<size_t> ch_ab_list;
      std::vector<size_t> iket_ab_list;
      std::vector<size_t> a_list;
      std::vector<size_t> b_list;
      std::vector<size_t> c_list;
      std::vector<double> occ_abc_list;

      for (int ch_ab = 0; ch_ab < nch; ch_ab++)
      {
        auto &tbc_ab = X.modelspace->GetTwoBodyChannel(ch_ab);
        //      int Jab = tbc_ab.J;
        size_t nkets_ab = tbc_ab.GetNumberKets();
        for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
        {
          Ket &ket_ab = tbc_ab.GetKet(iket_ab);
          int a = ket_ab.p;
          int b = ket_ab.q;
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          //        int ea = 2*oa.n + oa.l;
          //        int eb = 2*ob.n + ob.l;

          for (auto c : Z.modelspace->all_orbits)
          {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            //          double jc = 0.5*oc.j2;
            //          int ec = 2*oc.n + oc.l;
            //          if ( (std::abs( ea-e_fermi[oa.tz2]) + std::abs(eb-e_fermi[ob.tz2]) + std::abs(ec-e_fermi[oc.tz2])) > Z.modelspace->GetdE3max() ) continue;
            double occfactor = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
            if (std::abs(occfactor) < 1e-6)
              continue;
            if (std::abs(oa.tz2 + ob.tz2 - oc.tz2) > 2)
              continue;
            if (a == b)
              occfactor *= 0.5; // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b
                                //          if (  (std::abs( 2*Jab -oc.j2)>oj.j2)  or  ((2*Jab+oc.j2)<oj.j2) ) continue;
                                //          if (  (std::abs( 2*Jab -oc.j2)>std::max(oi.j2,oj.j2))  or  ((2*Jab+oc.j2)<std::min(oi.j2,oj.j2)) ) continue;
                                //          if ( tbc_ab.parity + oj.l + oc.l
            ch_ab_list.push_back(ch_ab);
            iket_ab_list.push_back(iket_ab);
            a_list.push_back(a);
            b_list.push_back(b);
            c_list.push_back(c);
            occ_abc_list.push_back(occfactor);
          } // for c
        }   // for iket_ab
      }     // for ch_ab

      Z.profiler.timer["_comm232_block2"] += omp_get_wtime() - tstart;
      tstart = omp_get_wtime();

      // allocate the matrices
      size_t n_abc = ch_ab_list.size();

      arma::mat X2MAT(nij, n_abc, arma::fill::zeros);
      arma::mat Y3MAT(n_abc, nij * njvals * nkets, arma::fill::zeros);
      arma::mat Y2MAT(nij, n_abc, arma::fill::zeros);
      arma::mat X3MAT(n_abc, nij * njvals * nkets, arma::fill::zeros);

      // Now fill the X2 mat and Y3 mat
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
      //      #pragma omp parallel for schedule(dynamic,1)
      for (size_t ind_abc = 0; ind_abc < n_abc; ind_abc++)
      {
        size_t a = a_list[ind_abc];
        size_t b = b_list[ind_abc];
        size_t c = c_list[ind_abc];
        int Jab = X.modelspace->GetTwoBodyChannel(ch_ab_list[ind_abc]).J;
        Orbit &oa = X.modelspace->GetOrbit(a);
        Orbit &ob = X.modelspace->GetOrbit(b);
        Orbit &oc = X.modelspace->GetOrbit(c);

        double jc = 0.5 * oc.j2;
        double d_ea = std::abs(2 * oa.n + oa.l - e_fermi[oa.tz2]);
        double d_eb = std::abs(2 * ob.n + ob.l - e_fermi[ob.tz2]);
        double d_ec = std::abs(2 * oc.n + oc.l - e_fermi[oc.tz2]);
        //        double occnat_a = oa.occ_nat;
        //        double occnat_b = ob.occ_nat;
        double occnat_c = oc.occ_nat;
        for (size_t ind_i = 0; ind_i < nij; ind_i++)
        {
          size_t i = ij_orbits[ind_i];
          //      if (ch==0) std::cout << "  ind_i = " << ind_i << "   i = " << i << std::endl;
          Orbit &oi = X.modelspace->GetOrbit(i);
          double ji = 0.5 * oi.j2;
          double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
          //      double occnat_i = oi.occ_nat;

          if ((d_ea + d_eb + d_ei) > Z.modelspace->GetdE3max())
            continue;

          X2MAT(ind_i, ind_abc) = -sqrt((2 * Jab + 1.)) * occ_abc_list[ind_abc] * X.TwoBody.GetTBME_J(Jab, c, i, a, b);
          Y2MAT(ind_i, ind_abc) = -sqrt((2 * Jab + 1.)) * occ_abc_list[ind_abc] * Y.TwoBody.GetTBME_J(Jab, c, i, a, b);

          if ((oa.l + ob.l + oi.l + tbc.parity + oc.l) % 2 != Y.parity)
            continue; // TODO: this will cause problems if we try something with Y.parity =1
          if (std::abs((oa.tz2 + ob.tz2 + oi.tz2) - (tbc.Tz * 2 + oc.tz2)) > 2 * Y.rank_T)
            continue;

          int twoJp_min = std::max(std::abs(2 * Jab - oi.j2), std::abs(2 * J - oc.j2));
          int twoJp_max = std::min(2 * Jab + oi.j2, 2 * J + oc.j2);

          // Why is this part structured like this?
          // In tests, the time for this entire commutator routine was dominated by the time to access 3-body matrix elements.
          // That's because I'm asking for them in an order different from the one in which they're stored, which means recoupling.
          // Fortunately, both the X and Y block need exactly the same recoupling, and we can pull the recoupling for the bra
          // side a few loops out.
          for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
          {
            if (!Y.ThreeBody.IsKetValid(Jab, twoJp, a, b, i))
              continue;
            std::vector<size_t> ibra_list;
            std::vector<double> recouple_bra_list;
            size_t ch_check = Y.ThreeBody.GetKetIndex_withRecoupling(Jab, twoJp, a, b, i, ibra_list, recouple_bra_list);

            for (size_t ind_jj = 0; ind_jj < njvals; ind_jj++)
            {
              int jj2 = j_jvals[ind_jj];
              double jj = 0.5 * jj2;
              double sixj = X.modelspace->GetSixJ(J, ji, jj, Jab, jc, 0.5 * twoJp);
              if (std::abs(sixj) < 1e-7)
                continue;
              // now loop over |kl> states
              for (int iket = 0; iket < nkets; iket++)
              {
                size_t index_kli = (ind_jj + ind_i * njvals) * nkets + iket;
                Ket &ket = tbc.GetKet(iket);
                int k = ket.p;
                int l = ket.q;
                double d_ek = std::abs(2 * ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
                double d_el = std::abs(2 * ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
                double occnat_k = ket.op->occ_nat;
                double occnat_l = ket.oq->occ_nat;
                if ((d_ek + d_el + d_ec) > Z.modelspace->GetdE3max())
                  continue;
                if ((occnat_k * (1 - occnat_k) * occnat_l * (1 - occnat_l) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
                  continue;
                if (!Y.ThreeBody.IsKetValid(J, twoJp, k, l, c))
                  continue;
                std::vector<size_t> iket_list;
                std::vector<double> recouple_ket_list;
                ch_check = Y.ThreeBody.GetKetIndex_withRecoupling(J, twoJp, k, l, c, iket_list, recouple_ket_list);

                double xabiklc = 0;
                double yabiklc = 0;
                for (size_t I = 0; I < ibra_list.size(); I++)
                {
                  for (size_t J = 0; J < iket_list.size(); J++)
                  {
                    // I explicitly check if x and y have 3-body components because the call GetME_pn_ch goes straight to the data array
                    // without a safety net. If the 3-body structure isn't allocated, then bad things will happen.
                    if (x_has_3)
                      xabiklc += recouple_bra_list[I] * recouple_ket_list[J] * X3.GetME_pn_ch(ch_check, ch_check, ibra_list[I], iket_list[J]);
                    if (y_has_3)
                      yabiklc += recouple_bra_list[I] * recouple_ket_list[J] * Y3.GetME_pn_ch(ch_check, ch_check, ibra_list[I], iket_list[J]);
                  }
                }

                X3MAT(ind_abc, index_kli) += (twoJp + 1) / sqrt(2 * J + 1) * sixj * xabiklc;
                Y3MAT(ind_abc, index_kli) += (twoJp + 1) / sqrt(2 * J + 1) * sixj * yabiklc;

              } // for iket
            }
          } // for ind_jj
        }   // for ind_abc

      } // for ind_i

      Z.profiler.timer["comm232_block3"] += omp_get_wtime() - tstart;
      tstart = omp_get_wtime();
      /// now we're back out to the ch loop level.

      // finally do the matrix multiplication
      arma::mat ZMat = (X2MAT * Y3MAT - Y2MAT * X3MAT);

      Z.profiler.timer["comm232_block4"] += omp_get_wtime() - tstart;
      tstart = omp_get_wtime();

      // now convert back from <i|Z|klj'> to <ij|Z|kl>
      //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
#pragma omp parallel for schedule(dynamic, 1)
      for (int ibra = 0; ibra < nkets; ibra++)
      {
        Ket &bra = tbc.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        size_t ind_i = ij_orbits_lookup[i];
        size_t ind_j = ij_orbits_lookup[j];
        size_t ind_ji = jvals_lookup[oi.j2];
        size_t ind_jj = jvals_lookup[oj.j2];

        for (int iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc.GetKet(iket);
          int k = ket.p;
          int l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          //        double jk = 0.5*ok.j2;
          //        double jl = 0.5*ol.j2;
          size_t ind_k = ij_orbits_lookup[k];
          size_t ind_l = ij_orbits_lookup[l];
          size_t ind_jk = jvals_lookup[ok.j2];
          size_t ind_jl = jvals_lookup[ol.j2];

          int phase_ij = X.modelspace->phase((oi.j2 + oj.j2) / 2 - J);
          int phase_kl = X.modelspace->phase((ok.j2 + ol.j2) / 2 - J);
          double zijkl = ZMat(ind_j, (ind_jj + ind_i * njvals) * nkets + iket);
          zijkl -= phase_ij * ZMat(ind_i, (ind_ji + ind_j * njvals) * nkets + iket);

          zijkl -= hermX * hermY * ZMat(ind_l, (ind_jl + ind_k * njvals) * nkets + ibra);
          zijkl += hermX * hermY * phase_kl * ZMat(ind_k, (ind_jk + ind_l * njvals) * nkets + ibra);

          // normalize the tbme
          zijkl *= -1.0 / sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
          Z2.AddToTBME(ch, ch, ibra, iket, zijkl);

        } // for iket
      }   // for ibra

      Z.profiler.timer["_comm232_block5"] += omp_get_wtime() - tstart;

    } // for ch

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  // the old way that also works. It's slower but easier to read.
  // For now, this is used if we want to treat an operator that changes Tz or parity.
  //  The fast implementation should be updated to do that...
  void comm232ss_slow(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X2 = X.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y2 = Y.TwoBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    int nch = Z.modelspace->GetNumberTwoBodyChannels();

    std::vector<std::array<size_t, 2>> channels;
    for (auto &iter : Z.TwoBody.MatEl)
      channels.push_back(iter.first);
    size_t nchans = channels.size();
    //  for (int ch=0; ch<nch; ch++)
    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nchans; ich++)
    {
      size_t ch_bra = channels[ich][0];
      size_t ch_ket = channels[ich][1];
      auto &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      auto &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        int ket_min = (ch_bra == ch_ket) ? ibra : 0;
        for (int iket = ket_min; iket < nkets; iket++)
        {
          double zijkl = 0;
          Ket &ket = tbc_ket.GetKet(iket);
          int k = ket.p;
          int l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          double jk = 0.5 * ok.j2;
          double jl = 0.5 * ol.j2;
          for (auto c : Z.modelspace->all_orbits)
          {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double jc = 0.5 * oc.j2;

            for (int ch_ab = 0; ch_ab < nch; ch_ab++)
            {
              auto &tbc_ab = X.modelspace->GetTwoBodyChannel(ch_ab);
              int Jab = tbc_ab.J;
              size_t nkets_ab = tbc_ab.GetNumberKets();
              for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
              {
                Ket &ket_ab = tbc_ab.GetKet(iket_ab);
                int a = ket_ab.p;
                int b = ket_ab.q;
                Orbit &oa = Z.modelspace->GetOrbit(a);
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double occfactor = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
                if (std::abs(occfactor) < 1e-6)
                  continue;
                if (a == b)
                  occfactor *= 0.5; // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b

                bool xicab_good = ((oi.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(oi.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                bool yicab_good = ((oi.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(oi.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                bool xjcab_good = ((oj.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(oj.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                bool yjcab_good = ((oj.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(oj.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                bool xabck_good = ((ok.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(ok.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                bool yabck_good = ((ok.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(ok.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                bool xabcl_good = ((ol.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(ol.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                bool yabcl_good = ((ol.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(ol.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);

                // Xicab term
                if ((xicab_good or yicab_good) and (std::abs(oi.j2 - oc.j2) <= 2 * Jab) and (oi.j2 + oc.j2 >= 2 * Jab))
                {

                  int twoJ_min = std::max(std::abs(oc.j2 - 2 * J), std::abs(oj.j2 - 2 * Jab));
                  int twoJ_max = std::min(oc.j2 + 2 * J, oj.j2 + 2 * Jab);
                  double xciab = X2.GetTBME_J(Jab, c, i, a, b);
                  double yciab = Y2.GetTBME_J(Jab, c, i, a, b);
                  int phasefactor = Z.modelspace->phase((oi.j2 + oj.j2) / 2 - J);
                  for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                  {

                    double Jtot = 0.5 * twoJ;
                    double sixj = Z.modelspace->GetSixJ(jj, ji, J, jc, Jtot, Jab);
                    double hatfactor = (twoJ + 1) * sqrt((2 * Jab + 1.) / (2 * J + 1));
                    double xabjklc = yicab_good ? X3.GetME_pn(Jab, J, twoJ, a, b, j, k, l, c) : 0;
                    double yabjklc = xicab_good ? Y3.GetME_pn(Jab, J, twoJ, a, b, j, k, l, c) : 0;
                    zijkl += occfactor * hatfactor * phasefactor * sixj * (xciab * yabjklc - yciab * xabjklc);
                  }
                }

                // Xjcab term
                if ((xjcab_good or yjcab_good) and (std::abs(oj.j2 - oc.j2) <= 2 * Jab) and (oj.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoJ_min = std::max(std::abs(oc.j2 - 2 * J), std::abs(oi.j2 - 2 * Jab));
                  int twoJ_max = std::min(oc.j2 + 2 * J, oi.j2 + 2 * Jab);
                  double xcjab = X2.GetTBME_J(Jab, c, j, a, b);
                  double ycjab = Y2.GetTBME_J(Jab, c, j, a, b);
                  int phasefactor = 1;
                  for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                  {
                    double Jtot = 0.5 * twoJ;
                    double sixj = Z.modelspace->GetSixJ(ji, jj, J, jc, Jtot, Jab);
                    double hatfactor = (twoJ + 1) * sqrt((2 * Jab + 1.) / (2 * J + 1));
                    double xabiklc = yjcab_good ? X3.GetME_pn(Jab, J, twoJ, a, b, i, k, l, c) : 0;
                    double yabiklc = xjcab_good ? Y3.GetME_pn(Jab, J, twoJ, a, b, i, k, l, c) : 0;
                    zijkl -= occfactor * hatfactor * phasefactor * sixj * (xcjab * yabiklc - ycjab * xabiklc);
                  }
                }

                // Xabck term
                if ((xabck_good or yabck_good) and (std::abs(ok.j2 - oc.j2) <= 2 * Jab) and (ok.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoJ_min = std::max(std::abs(oc.j2 - 2 * J), std::abs(ol.j2 - 2 * Jab));
                  int twoJ_max = std::min(oc.j2 + 2 * J, ol.j2 + 2 * Jab);
                  double xabck = X2.GetTBME_J(Jab, a, b, c, k);
                  double yabck = Y2.GetTBME_J(Jab, a, b, c, k);
                  int phasefactor = Z.modelspace->phase((ok.j2 + ol.j2) / 2 - J);
                  for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                  {
                    double Jtot = 0.5 * twoJ;
                    double sixj = Z.modelspace->GetSixJ(jl, jk, J, jc, Jtot, Jab);
                    double hatfactor = (twoJ + 1) * sqrt((2 * Jab + 1.) / (2 * J + 1));
                    double xijcabl = yabck_good ? X3.GetME_pn(J, Jab, twoJ, i, j, c, a, b, l) : 0;
                    double yijcabl = xabck_good ? Y3.GetME_pn(J, Jab, twoJ, i, j, c, a, b, l) : 0;
                    zijkl -= occfactor * hatfactor * phasefactor * sixj * (yijcabl * xabck - xijcabl * yabck);
                  }
                }

                // Xabcl term
                if ((xabcl_good or yabcl_good) and (std::abs(ol.j2 - oc.j2) <= 2 * Jab) and (ol.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoJ_min = std::max(std::abs(oc.j2 - 2 * J), std::abs(ok.j2 - 2 * Jab));
                  int twoJ_max = std::min(oc.j2 + 2 * J, ok.j2 + 2 * Jab);
                  double xabcl = X2.GetTBME_J(Jab, a, b, c, l);
                  double yabcl = Y2.GetTBME_J(Jab, a, b, c, l);
                  int phasefactor = 1;
                  for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                  {
                    double Jtot = 0.5 * twoJ;
                    double sixj = Z.modelspace->GetSixJ(jk, jl, J, jc, Jtot, Jab);
                    double hatfactor = (twoJ + 1) * sqrt((2 * Jab + 1.) / (2 * J + 1));
                    double xijcabk = yabcl_good ? X3.GetME_pn(J, Jab, twoJ, i, j, c, a, b, k) : 0;
                    double yijcabk = xabcl_good ? Y3.GetME_pn(J, Jab, twoJ, i, j, c, a, b, k) : 0;
                    zijkl += occfactor * hatfactor * phasefactor * sixj * (yijcabk * xabcl - xijcabk * yabcl);
                  }
                }

              } // for iket_ab
            }   // for ch2
          }     // for c

          // normalize the tbme
          zijkl *= -1.0 / sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
        } // for iket
      }   // for ibra
    }     // for ch

  } // comm232ss

  //*****************************************************************************************
  //
  // i|  j|
  //  |   |           Uncoupled expression:
  //  *~~[X]~~~~*        Z_ijkl = 1/6 sum_abcd (nanbncn`d -n`an`bn`cnd)(X_ijdabc Y_abckld - Y_ijdabc Xabckld)
  //  |   |     /\      .
  // a|  b|   c(  )d    .
  //  |   |     \/      .
  //  *~~[Y]~~~~*      Coupled expression:
  //  |   |              Z_{ijkl}^{J} = 1/6 sum_abcd (nanbncn`d-n`an`bn`cnd) 1/(2J+1) sum_J1J' (2J'+1)(X_{ijdabc}^{J J1 J'} Y_{abckld}^{J1 J J'} - X<->Y )
  // k|  l|
  //
  //
  //  Checked with UnitTest and passed.
  //
  void comm332_ppph_hhhpss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    bool x3_allocated = X3.IsAllocated();
    bool y3_allocated = Y3.IsAllocated();
    if (not(x3_allocated and y3_allocated))
      return;

    int hX = X.IsHermitian() ? 1 : -1;
    int hY = Y.IsHermitian() ? 1 : -1;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    int nch = Z.modelspace->GetNumberTwoBodyChannels();
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
    Z.modelspace->PreCalculateSixJ(); // if this is already done, then this does nothing.
                                      //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch = 0; ch < nch; ch++)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nkets = tbc.GetNumberKets();

      std::vector<size_t> ch3_abc_list;
      std::vector<size_t> iket_abc_list;
      std::vector<size_t> d_list;
      std::vector<double> factor_list;

      // figure out how big the abcd side of the matrices should be, and store some stuff for later lookup
      for (size_t ch_abc = 0; ch_abc < nch3; ch_abc++)
      {
        ThreeBodyChannel &Tbc_abc = Z.modelspace->GetThreeBodyChannel(ch_abc);
        if (std::abs(Tbc_abc.twoTz - 2 * tbc.Tz) > 1)
          continue;
        int twoJ = Tbc_abc.twoJ;
        size_t nkets_abc = Tbc_abc.GetNumberKets();
        for (size_t iket_abc = 0; iket_abc < nkets_abc; iket_abc++)
        {
          Ket3 &ket_abc = Tbc_abc.GetKet(iket_abc);
          size_t a = ket_abc.p;
          size_t b = ket_abc.q;
          size_t c = ket_abc.r;
          double occ_abc = ket_abc.op->occ * ket_abc.oq->occ * ket_abc.oR->occ;
          double occ_abc_bar = (1 - ket_abc.op->occ) * (1 - ket_abc.oq->occ) * (1 - ket_abc.oR->occ);
          double d_ea = std::abs(2 * ket_abc.op->n + ket_abc.op->l - e_fermi[ket_abc.op->tz2]);
          double d_eb = std::abs(2 * ket_abc.oq->n + ket_abc.oq->l - e_fermi[ket_abc.oq->tz2]);
          double d_ec = std::abs(2 * ket_abc.oR->n + ket_abc.oR->l - e_fermi[ket_abc.oR->tz2]);
          double occnat_a = ket_abc.op->occ_nat;
          double occnat_b = ket_abc.oq->occ_nat;
          double occnat_c = ket_abc.oR->occ_nat;
          if ((d_ea + d_eb + d_ec) > Z.modelspace->dE3max)
            continue;
          if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
            continue;
          if ((std::abs(occ_abc) == 0) and (std::abs(occ_abc_bar) == 0))
            continue;
          //        int Jab = ket_abc.Jpq;

          double symm_factor = 6; // 6 possible orderings of abc. If a==b, then only 3 orderings, and if a==b==c, then only 1 ordering.
          if (a == b and b == c)
            symm_factor = 1;
          else if (a == b or a == c or b == c)
            symm_factor = 3;

          for (auto d : Z.modelspace->all_orbits)
          {
            Orbit &od = Z.modelspace->GetOrbit(d);
            double nd = od.occ;
            //          double d_ed = std::abs( 2*od.n + od.l - e_fermi[od.tz2]);
            //          double occnat_d = od.occ_nat;
            double occfactor = occ_abc * (1 - nd) - occ_abc_bar * nd;
            if (std::abs(occfactor) < 1e-6)
              continue;
            if ((std::abs(2 * J - od.j2) > twoJ) or (2 * J + od.j2) < twoJ)
              continue;
            if ((od.l + tbc.parity + Tbc_abc.parity) % 2 > 0)
              continue;
            if ((Tbc_abc.twoTz) != (od.tz2 + 2 * tbc.Tz))
              continue;

            ch3_abc_list.push_back(ch_abc);
            iket_abc_list.push_back(iket_abc);
            d_list.push_back(d);
            factor_list.push_back((twoJ + 1.) / (2 * J + 1) * occfactor * symm_factor / 6.);
          } // for d
        }   // for iket_abc
      }     // for ch_abc

      // initialize the matrices
      size_t dim_abcd = ch3_abc_list.size();
      arma::mat XMAT(nkets, dim_abcd, arma::fill::zeros);
      arma::mat YMAT(nkets, dim_abcd, arma::fill::zeros);

      // fill the matrices  TODO: I think we can make this faster by first doing the recoupling lookup for ijd
      for (size_t index_abcd = 0; index_abcd < dim_abcd; index_abcd++)
      {
        size_t ch_abc = ch3_abc_list[index_abcd];
        size_t iket_abc = iket_abc_list[index_abcd];
        double factor = factor_list[index_abcd];
        ThreeBodyChannel &Tbc_abc = Z.modelspace->GetThreeBodyChannel(ch_abc);
        Ket3 &ket_abc = Tbc_abc.GetKet(iket_abc);
        int Jab = ket_abc.Jpq;
        int twoJ = Tbc_abc.twoJ;
        size_t a = ket_abc.p;
        size_t b = ket_abc.q;
        size_t c = ket_abc.r;
        size_t d = d_list[index_abcd];
        Orbit &od = Z.modelspace->GetOrbit(d);
        double d_ed = std::abs(2 * od.n + od.l - e_fermi[od.tz2]);
        double occnat_d = od.occ_nat;

        for (int ibra = 0; ibra < nkets; ibra++)
        {
          Ket &bra = tbc.GetKet(ibra);
          int i = bra.p;
          int j = bra.q;
          double d_ei = std::abs(2 * bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
          double d_ej = std::abs(2 * bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
          double occnat_i = bra.op->occ_nat;
          double occnat_j = bra.oq->occ_nat;
          if ((d_ei + d_ej + d_ed) > Z.modelspace->dE3max)
            continue;
          if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_d * (1 - occnat_d)) < Z.modelspace->GetOccNat3Cut())
            continue;
          double norm2b = (i == j) ? 1. / sqrt(2) : 1;

          XMAT(ibra, index_abcd) = norm2b * X3.GetME_pn(J, Jab, twoJ, i, j, d, a, b, c);
          YMAT(ibra, index_abcd) = norm2b * factor * Y3.GetME_pn(J, Jab, twoJ, i, j, d, a, b, c);

        } // for ibra
      }   // for index_abcd

      // do the mat mult
      arma::mat ZMAT = hY * XMAT * YMAT.t() - hX * YMAT * XMAT.t();

      // now unpack
      for (int ibra = 0; ibra < nkets; ibra++)
      {
        for (int iket = ibra; iket < nkets; iket++)
        {
          Z2.AddToTBME(ch, ch, ibra, iket, ZMAT(ibra, iket));
        }
      }
    } // for ch
    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  // the old slow way
  /*
  void comm332_ppph_hhhpss( const Operator& X, const Operator& Y, Operator& Z )
  {
    double tstart = omp_get_wtime();
    auto& X3 = X.ThreeBody;
    auto& Y3 = Y.ThreeBody;
    auto& Z2 = Z.TwoBody;

    int nch = Z.modelspace->GetNumberTwoBodyChannels();
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
    #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (int ch=0; ch<nch; ch++)
    {
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nkets = tbc.GetNumberKets();
      for (int ibra=0; ibra<nkets; ibra++)
      {
        Ket& bra = tbc.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        for (int iket=ibra; iket<nkets; iket++)
        {
          Ket& ket = tbc.GetKet(iket);
          int k = ket.p;
          int l = ket.q;

          double zijkl = 0;

          for (size_t ch_abc=0; ch_abc<nch3; ch_abc++)
          {
            ThreeBodyChannel& Tbc_abc = Z.modelspace->GetThreeBodyChannel(ch_abc);
            if ( std::abs( Tbc_abc.twoTz - 2*tbc.Tz) > 1) continue;
            int twoJ = Tbc_abc.twoJ;
            size_t nkets_abc = Tbc_abc.GetNumberKets();
            for (size_t iket_abc=0; iket_abc<nkets_abc; iket_abc++)
            {
              Ket3& ket_abc = Tbc_abc.GetKet(iket_abc);
              size_t a = ket_abc.p;
              size_t b = ket_abc.q;
              size_t c = ket_abc.r;
              double occ_abc = ket_abc.op->occ * ket_abc.oq->occ * ket_abc.oR->occ;
              double occ_abc_bar = (1-ket_abc.op->occ) * (1-ket_abc.oq->occ) * (1-ket_abc.oR->occ);
              if ( (std::abs(occ_abc)==0) and (std::abs(occ_abc_bar)==0) ) continue;
              int Jab = ket_abc.Jpq;

                double symm_factor = 6;  // 6 possible orderings of abc. If a==b, then only 3 orderings, and if a==b==c, then only 1 ordering.
                if ( a==b and b==c )
                   symm_factor = 1;
                else if (a==b or a==c or b==c )
                   symm_factor = 3;

                for ( auto d : Z.modelspace->all_orbits )
                {
                  Orbit& od = Z.modelspace->GetOrbit(d);
                  double nd = od.occ;
                  double occfactor = occ_abc*(1-nd) - occ_abc_bar*nd ;
                  if (std::abs(occfactor)<1e-6) continue;
                  if ( (std::abs( 2*J-od.j2) > twoJ) or (2*J+od.j2)<twoJ ) continue;
                  if ( (od.l+tbc.parity + Tbc_abc.parity)%2>0) continue;
                  if ( (Tbc_abc.twoTz)!=(od.tz2 + 2*tbc.Tz) ) continue;
                  zijkl += symm_factor * 1./6  * (twoJ+1.)/(2*J+1) * occfactor
                           *(  X3.GetME_pn(J,Jab,twoJ, i,j,d,a,b,c) * Y3.GetME_pn(Jab, J, twoJ, a,b,c,k,l,d)
                             - Y3.GetME_pn(J,Jab,twoJ, i,j,d,a,b,c) * X3.GetME_pn(Jab, J, twoJ, a,b,c,k,l,d) );

                }// for d
            }// for iket_abc
          }// for ch_abc
          zijkl /=  sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
          Z2.AddToTBME(ch,ch,ibra,iket,zijkl);
        }// for iket
      }// for ibra
    }// for ch
    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }
  */

  //*****************************************************************************************
  //
  //           i| j|
  //            |  |  Uncoupled expression:
  //   *~~~[X]~~*  |            Z_ijkl = 1/6 sum_abcd (nanbncn`d -n`an`bn`cnd)(X_ijdabc Y_abckld - Y_ijdabc Xabckld)
  //  / \  / \  |  |                  + 1/4 (1-Pij)(1-Pkl) sum_{abcd} nanbn`cn`d ( X_abicdk Y_cdjabl - Y_abicdk X_cdjabl )
  // (a c)(b d) |  |
  //  \ /  \ /  |  |
  //   *~~~[Y]~~|~~*
  //            |  |   Coupled expression:
  //           k| l|     Z_{ijkl}^{J} = 1/6 sum_abcd (nanbncn`d-n`an`bn`cnd) 1/(2J+1) sum_J1J' (2J'+1)(X_{ijdabc}^{J J1 J'} Y_{abckld}^{J1 J J'} - X<->Y )
  //                                    + 1/4  (1-Pij^J)(1-Pkl^J) sum_{abcd} nanbn`cn`d sum_{J1 J2 J' J"} (2J'+1)(2J"+1)  (-1)^{2J+J2+J'+J"+2j-i-k}
  //                                      { p  q  J }
  //                                    * { J1 J" s } ( X_{abicdk}^{J1 J2 J'} Y_{cdjabl}^{J2 J1 J"} - X<->Y )
  //                                      { J' J2 r }
  //
  //  We recouple the expression so that it factorizes and can be cast as matrix multiplication
  //
  //            i| /k`
  //             |/
  //    *~~~[X]~~*
  //   / \  / \
//  (a c)(b d)
  //   \ /  \ /
  //    *~~~[Y]~~~*
  //              | \
//             l|  \j'
  //
  //        Tested with UnitTest and passed.
  //
  void comm332_pphhss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    //  double t_internal = omp_get_wtime();
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();
    //  std::cout << "  fermi levels " << e_fermi[-1] << " " << e_fermi[+1] << std::endl;

    bool x3_allocated = X3.IsAllocated();
    bool y3_allocated = Y3.IsAllocated();
    if (not(x3_allocated and y3_allocated))
      return;

    size_t nch = Z.modelspace->GetNumberTwoBodyChannels();
    size_t nch_CC = Z.modelspace->TwoBodyChannels_CC.size();

    std::deque<arma::mat> Z_bar(nch_CC);
    std::deque<arma::mat> Z_bar_flipbra(nch_CC);
    std::deque<arma::mat> Z_bar_flipket(nch_CC);
    std::deque<arma::mat> Z_bar_flipboth(nch_CC);

    // before entering the parallel loop, lets allocate the Pandya-transformed matrices Zbar.
    // we keep 4 because we want i<j and i>j and k<l and k>l
    for (size_t ch = 0; ch < nch_CC; ++ch)
    {
      TwoBodyChannel_CC &tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_CC.GetNumberKets();

      Z_bar[ch].zeros(nKets_cc, nKets_cc);
      Z_bar_flipbra[ch].zeros(nKets_cc, nKets_cc);
      Z_bar_flipket[ch].zeros(nKets_cc, nKets_cc);
      Z_bar_flipboth[ch].zeros(nKets_cc, nKets_cc);
    }

    double occnat_factor_max = 0;
    for (auto i : Z.modelspace->all_orbits)
    {
      double occnat_i = Z.modelspace->GetOrbit(i).occ_nat;
      occnat_factor_max = std::max(occnat_factor_max, occnat_i * (1 - occnat_i));
    }

    Z.modelspace->PreCalculateSixJ(); // if it's been done, this won't do it again. But make sure for thread safety.
    // But we still need some SixJs which aren't computed in the current implementation of PreCalculateSixJ().
    // we need symbols of the form { J2b J2b J2b }
    //                             { j1b j1b J3b }
    if (Z.modelspace->scalar3b_transform_first_pass)
    {
      double t_internal = omp_get_wtime();
      int j_1b_max = Z.modelspace->OneBodyJmax;
      int j_2b_max = Z.modelspace->TwoBodyJmax;
      for (int j2i = 1; j2i <= j_1b_max; j2i += 2)
      {
        for (int j2j = 1; j2j <= j_1b_max; j2j += 2)
        {
          for (int J1 = 0; J1 <= j_2b_max; J1++)
          {
            for (int J2 = 0; J2 <= j_2b_max; J2++)
            {
              int J3min = std::max(std::abs(j2i - j2j) / 2, std::abs(J1 - J2));
              int J3max = std::min(std::abs(j2i + j2j) / 2, std::abs(J1 + J2));
              int twoJ_min = std::max(std::abs(2 * J1 - j2j), std::abs(2 * J2 - j2i));
              int twoJ_max = std::min(2 * J1 + j2j, 2 * J2 + j2i);
              for (int J3 = J3min; J3 <= J3max; J3++)
              {
                for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                {
                  Z.modelspace->GetSixJ(J1, J2, J3, 0.5 * j2i, 0.5 * j2j, 0.5 * twoJ);
                }
              }
            }
          }
        }
      }
      Z.profiler.timer["__precalc6J_comm332_pphhss"] += omp_get_wtime() - t_internal;
    }

    // Loop through Pandya-transformed channels and compute the matrix Zbar by mat mult
    //  #pragma omp parallel for schedule(dynamic,1)
    //  There are some SixJs in this loop which are not automatically computed with PreCalculateSixJ, so
    //  we need to be careful the first time we call this.
#pragma omp parallel for schedule(dynamic, 1) // if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch = 0; ch < nch_CC; ++ch)
    {
      TwoBodyChannel_CC &tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_CC.GetNumberKets();

      int J_ph = tbc_CC.J;
      //      int parity_ph = tbc_CC.parity;
      //      int Tz_ph = tbc_CC.Tz;

      size_t n_abcd_states = 0;

      // count the abcd states that can act in this channel
      for (size_t ch_ab = 0; ch_ab < nch; ch_ab++)
      {
        TwoBodyChannel &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
        int Jab = tbc_ab.J;
        int nkets_ab = tbc_ab.GetNumberKets();

        for (size_t ch_cd = 0; ch_cd < nch; ch_cd++)
        {
          TwoBodyChannel &tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
          int Jcd = tbc_cd.J;
          int nkets_cd = tbc_cd.GetNumberKets();
          if (std::abs(Jab - Jcd) > J_ph or (Jab + Jcd) < J_ph)
            continue;
          if ((tbc_ab.parity + tbc_cd.parity + tbc_CC.parity) % 2 > 0)
            continue;
          if (std::abs(tbc_cd.Tz - tbc_ab.Tz) != tbc_CC.Tz)
            continue;
          for (int iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
          {
            Ket &ket_ab = tbc_ab.GetKet(iket_ab);
            //            int a = ket_ab.p;
            //            int b = ket_ab.q;
            double na = ket_ab.op->occ;
            double nb = ket_ab.oq->occ;
            double occnat_a = ket_ab.op->occ_nat;
            double occnat_b = ket_ab.oq->occ_nat;
            double d_ea = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
            double d_eb = std::abs(2 * ket_ab.oq->n + ket_ab.oq->l - e_fermi[ket_ab.oq->tz2]);
            if (na * nb < 1e-8 and (1 - na) * (1 - nb) < 1e-8)
              continue;
            if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_factor_max) < Z.modelspace->GetOccNat3Cut())
              continue;
            if ((d_ea + d_eb) > Z.modelspace->GetdE3max())
              continue;

            for (int iket_cd = 0; iket_cd < nkets_cd; iket_cd++)
            {
              Ket &ket_cd = tbc_cd.GetKet(iket_cd);
              //              int c = ket_cd.p;
              //              int d = ket_cd.q;
              double nc = ket_cd.op->occ;
              double nd = ket_cd.oq->occ;
              double occnat_c = ket_cd.op->occ_nat;
              double occnat_d = ket_cd.oq->occ_nat;
              double d_ec = std::abs(2 * ket_cd.op->n + ket_cd.op->l - e_fermi[ket_cd.op->tz2]);
              double d_ed = std::abs(2 * ket_cd.oq->n + ket_cd.oq->l - e_fermi[ket_cd.oq->tz2]);
              double occupation_factor = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);
              if (std::abs(occupation_factor) < 1e-6)
                continue;
              if ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_factor_max) < Z.modelspace->GetOccNat3Cut())
                continue;
              if ((d_ec + d_ed) > Z.modelspace->GetdE3max())
                continue;

              n_abcd_states++;

            } // for iket_cd
          }   // for iket_ab
        }     // for ch_cd
      }       // for ch_ab

      // Now that we know how big they are, allocate the matrices
      arma::mat X3_ij(nKets_cc, n_abcd_states, arma::fill::zeros);
      arma::mat X3_ji(nKets_cc, n_abcd_states, arma::fill::zeros);
      arma::mat Y3_ij(n_abcd_states, nKets_cc, arma::fill::zeros);
      arma::mat Y3_ji(n_abcd_states, nKets_cc, arma::fill::zeros);

      // Fill the matrices, which are basically <ij`| X | cd a`b`> type things (the ticks ` mean time-reversed state)
      for (size_t ibra_CC = 0; ibra_CC < nKets_cc; ibra_CC++)
      {
        Ket &bra_CC = tbc_CC.GetKet(ibra_CC);
        int i = bra_CC.p;
        int j = bra_CC.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);

        int j2i = oi.j2;
        int j2j = oj.j2;
        double ji = 0.5 * j2i;
        double jj = 0.5 * j2j;
        //        int phase_ij = Z.modelspace->phase( (j2i+j2j)/2);

        double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;

        size_t ind_abcd = -1;
        for (size_t ch_ab = 0; ch_ab < nch; ch_ab++)
        {
          TwoBodyChannel &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
          int Jab = tbc_ab.J;
          int nkets_ab = tbc_ab.GetNumberKets();

          for (size_t ch_cd = 0; ch_cd < nch; ch_cd++)
          {
            TwoBodyChannel &tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
            int Jcd = tbc_cd.J;
            int nkets_cd = tbc_cd.GetNumberKets();
            if (std::abs(Jab - Jcd) > J_ph or (Jab + Jcd) < J_ph)
              continue;
            if ((tbc_ab.parity + tbc_cd.parity + tbc_CC.parity) % 2 > 0)
              continue;
            if (std::abs(tbc_cd.Tz - tbc_ab.Tz) != tbc_CC.Tz)
              continue;

            // I erroneously thought I'd fixed a bug here in commit 7ffbe... but it actually gave wrong answers. This seems to have corrected the oopsie. --SRS
            int twoJp_min = std::min(std::max(std::abs(2 * Jab - j2i), std::abs(2 * Jcd - j2j)), std::max(std::abs(2 * Jab - j2j), std::abs(2 * Jcd - j2i)));
            int twoJp_max = std::max(std::min((2 * Jab + j2i), (2 * Jcd + j2j)), std::min((2 * Jab + j2j), (2 * Jcd + j2i)));
            //            int twoJp_min = std::max( {std::abs(2*Jab-j2i), std::abs(2*Jcd-j2j), std::abs(2*Jab-j2j), std::abs(2*Jcd-j2i) });
            //            int twoJp_max = std::min( {(2*Jab+j2i), (2*Jcd+j2j), (2*Jab+j2j), (2*Jcd+j2i) });
            std::vector<double> sixj_ij_list;
            std::vector<double> sixj_ji_list;
            for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
            {
              //              sixj_ij_list.push_back( (twoJp+1) * Z.modelspace->GetSixJ(ji,jj,J_ph, Jcd,Jab,0.5*twoJp) * Z.modelspace->phase((j2i+twoJp)/2) );
              //              sixj_ji_list.push_back( (twoJp+1) * Z.modelspace->GetSixJ(jj,ji,J_ph, Jcd,Jab,0.5*twoJp) * Z.modelspace->phase((j2j+twoJp)/2) );
              double sixj_ij = 0;
              double sixj_ji = 0;
              if (AngMom::Triangle(Jcd, jj, 0.5 * twoJp) and AngMom::Triangle(Jab, ji, 0.5 * twoJp))
              {
                sixj_ij = Z.modelspace->GetSixJ(Jcd, Jab, J_ph, ji, jj, 0.5 * twoJp);
              }
              if (AngMom::Triangle(Jcd, ji, 0.5 * twoJp) and AngMom::Triangle(Jab, jj, 0.5 * twoJp))
              {
                sixj_ji = Z.modelspace->GetSixJ(Jcd, Jab, J_ph, jj, ji, 0.5 * twoJp);
              }

              //              sixj_ij_list.push_back( (twoJp+1) * Z.modelspace->GetSixJ(Jcd,Jab,J_ph, ji,jj,0.5*twoJp) * Z.modelspace->phase((j2i+twoJp)/2) );
              //              sixj_ji_list.push_back( (twoJp+1) * Z.modelspace->GetSixJ(Jcd,Jab,J_ph, jj,ji,0.5*twoJp) * Z.modelspace->phase((j2j+twoJp)/2) );
              sixj_ij_list.push_back((twoJp + 1) * sixj_ij * Z.modelspace->phase((j2i + twoJp) / 2));
              sixj_ji_list.push_back((twoJp + 1) * sixj_ji * Z.modelspace->phase((j2j + twoJp) / 2));
            }

            bool abi_cdj_Tz_ok = ((2 * tbc_ab.Tz + oi.tz2) == (2 * tbc_cd.Tz + oj.tz2));
            bool abj_cdi_Tz_ok = ((2 * tbc_ab.Tz + oj.tz2) == (2 * tbc_cd.Tz + oi.tz2));

            for (int iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
            {
              Ket &ket_ab = tbc_ab.GetKet(iket_ab);
              int a = ket_ab.p;
              int b = ket_ab.q;
              double na = ket_ab.op->occ;
              double nb = ket_ab.oq->occ;
              if (na * nb < 1e-8 and (1 - na) * (1 - nb) < 1e-8)
                continue;
              double d_ea = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
              double d_eb = std::abs(2 * ket_ab.oq->n + ket_ab.oq->l - e_fermi[ket_ab.oq->tz2]);
              double occnat_a = ket_ab.op->occ_nat;
              double occnat_b = ket_ab.oq->occ_nat;
              bool keep_abi = ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_i * (1 - occnat_i)) >= Z.modelspace->GetOccNat3Cut());
              bool keep_abj = ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) >= Z.modelspace->GetOccNat3Cut());
              keep_abi = keep_abi and (d_ea + d_eb + d_ei <= Z.modelspace->dE3max);
              keep_abj = keep_abj and (d_ea + d_eb + d_ej <= Z.modelspace->dE3max);

              if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_factor_max) < Z.modelspace->GetOccNat3Cut())
                continue;
              if ((d_ea + d_eb) > Z.modelspace->GetdE3max())
                continue;

              std::vector<size_t> ch_abi_list;
              std::vector<size_t> ch_abj_list;
              std::vector<std::vector<size_t>> kets_abi_list;
              std::vector<std::vector<size_t>> kets_abj_list;
              std::vector<std::vector<double>> recouple_abi_list;
              std::vector<std::vector<double>> recouple_abj_list;

              if (keep_abi or keep_abj)
              {
                for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                {
                  std::vector<size_t> ketlist;
                  std::vector<double> reclist;
                  size_t ch_check = -1;
                  if (keep_abi and abi_cdj_Tz_ok and (std::abs(2 * Jab - oi.j2) <= twoJp) and (2 * Jab + oi.j2 >= twoJp) and Y.ThreeBody.IsKetValid(Jab, twoJp, a, b, i))
                  {
                    ch_check = Y.ThreeBody.GetKetIndex_withRecoupling(Jab, twoJp, a, b, i, ketlist, reclist);
                  }
                  ch_abi_list.push_back(ch_check);
                  kets_abi_list.push_back(ketlist);
                  recouple_abi_list.push_back(reclist);
                  ketlist.clear();
                  reclist.clear();
                  ch_check = -1;
                  if (keep_abj and abj_cdi_Tz_ok and (std::abs(2 * Jab - oj.j2) <= twoJp) and (2 * Jab + oj.j2 >= twoJp) and Y.ThreeBody.IsKetValid(Jab, twoJp, a, b, j))
                  {
                    ch_check = Y.ThreeBody.GetKetIndex_withRecoupling(Jab, twoJp, a, b, j, ketlist, reclist);
                  }
                  ch_abj_list.push_back(ch_check);
                  kets_abj_list.push_back(ketlist);
                  recouple_abj_list.push_back(reclist);
                }
              }

              for (int iket_cd = 0; iket_cd < nkets_cd; iket_cd++)
              {
                Ket &ket_cd = tbc_cd.GetKet(iket_cd);
                int c = ket_cd.p;
                int d = ket_cd.q;
                double nc = ket_cd.op->occ;
                double nd = ket_cd.oq->occ;
                double d_ec = std::abs(2 * ket_cd.op->n + ket_cd.op->l - e_fermi[ket_cd.op->tz2]);
                double d_ed = std::abs(2 * ket_cd.oq->n + ket_cd.oq->l - e_fermi[ket_cd.oq->tz2]);
                double occnat_c = ket_cd.op->occ_nat;
                double occnat_d = ket_cd.oq->occ_nat;
                double occupation_factor = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);
                if (std::abs(occupation_factor) < 1e-6)
                  continue;

                if ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_factor_max) < Z.modelspace->GetOccNat3Cut())
                  continue;
                if ((d_ec + d_ed) > Z.modelspace->GetdE3max())
                  continue;

                double symmetry_factor = 1; // we only sum a<=b and c<=d, so we undercount by a factor of 4, canceling the 1/4 in the formula
                if (a == b)
                  symmetry_factor *= 0.5; // if a==b or c==d, then the permutation doesn't give a new state, so there's less undercounting
                if (c == d)
                  symmetry_factor *= 0.5;

                ind_abcd++; // increment the matrix index.

                if (not(keep_abi or keep_abj))
                  continue; // need this after the ind_abcd++... maybe be more clever?

                bool keep_cdi = ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_i * (1 - occnat_i)) >= Z.modelspace->GetOccNat3Cut());
                bool keep_cdj = ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_j * (1 - occnat_j)) >= Z.modelspace->GetOccNat3Cut());
                keep_cdi = keep_cdi and (d_ec + d_ed + d_ei <= Z.modelspace->dE3max);
                keep_cdj = keep_cdj and (d_ec + d_ed + d_ej <= Z.modelspace->dE3max);

                if (not(keep_cdi or keep_cdj))
                  continue;

                std::vector<size_t> ch_cdi_list;
                std::vector<size_t> ch_cdj_list;
                std::vector<std::vector<size_t>> kets_cdi_list;
                std::vector<std::vector<size_t>> kets_cdj_list;
                std::vector<std::vector<double>> recouple_cdi_list;
                std::vector<std::vector<double>> recouple_cdj_list;

                for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                {
                  std::vector<size_t> ketlist;
                  std::vector<double> reclist;
                  size_t ch_check = -1;
                  if (keep_abj and keep_cdi and abj_cdi_Tz_ok and (std::abs(2 * Jcd - oi.j2) <= twoJp) and (2 * Jcd + oi.j2 >= twoJp) and Y.ThreeBody.IsKetValid(Jcd, twoJp, c, d, i))
                  {
                    ch_check = Y.ThreeBody.GetKetIndex_withRecoupling(Jcd, twoJp, c, d, i, ketlist, reclist);
                  }
                  ch_cdi_list.push_back(ch_check);
                  kets_cdi_list.push_back(ketlist);
                  recouple_cdi_list.push_back(reclist);
                  ketlist.clear();
                  reclist.clear();
                  ch_check = -1;
                  if (keep_abi and keep_cdj and abi_cdj_Tz_ok and (std::abs(2 * Jcd - oj.j2) <= twoJp) and (2 * Jcd + oj.j2 >= twoJp) and Y.ThreeBody.IsKetValid(Jcd, twoJp, c, d, j))
                  {
                    ch_check = Y.ThreeBody.GetKetIndex_withRecoupling(Jcd, twoJp, c, d, j, ketlist, reclist);
                  }
                  ch_cdj_list.push_back(ch_check);
                  kets_cdj_list.push_back(ketlist);
                  recouple_cdj_list.push_back(reclist);
                } // for twoJp

                for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                {
                  double sixj_ij = sixj_ij_list[(twoJp - twoJp_min) / 2]; // <-- these also include a phase and a (twoJp+1)
                  double sixj_ji = sixj_ji_list[(twoJp - twoJp_min) / 2];

                  if (keep_abi and keep_cdj and abi_cdj_Tz_ok)
                  {
                    size_t ch_abi = ch_abi_list[(twoJp - twoJp_min) / 2];
                    auto &kets_abi = kets_abi_list[(twoJp - twoJp_min) / 2];
                    auto &kets_cdj = kets_cdj_list[(twoJp - twoJp_min) / 2];
                    auto &recouple_abi = recouple_abi_list[(twoJp - twoJp_min) / 2];
                    auto &recouple_cdj = recouple_cdj_list[(twoJp - twoJp_min) / 2];
                    double xabicdj = 0;
                    double ycdjabi = 0;
                    for (size_t Iabi = 0; Iabi < kets_abi.size(); Iabi++)
                    {
                      for (size_t Icdj = 0; Icdj < kets_cdj.size(); Icdj++)
                      {
                        xabicdj += recouple_abi[Iabi] * recouple_cdj[Icdj] * X3.GetME_pn_ch(ch_abi, ch_abi, kets_abi[Iabi], kets_cdj[Icdj]);
                        ycdjabi += recouple_abi[Iabi] * recouple_cdj[Icdj] * Y3.GetME_pn_ch(ch_abi, ch_abi, kets_cdj[Icdj], kets_abi[Iabi]);
                      }
                    }
                    X3_ij(ibra_CC, ind_abcd) += sixj_ij * xabicdj;
                    Y3_ij(ind_abcd, ibra_CC) += sixj_ij * ycdjabi * occupation_factor * symmetry_factor;
                  }

                  if (keep_abj and keep_cdi and abj_cdi_Tz_ok)
                  {
                    size_t ch_abj = ch_abj_list[(twoJp - twoJp_min) / 2];
                    auto &kets_abj = kets_abj_list[(twoJp - twoJp_min) / 2];
                    auto &kets_cdi = kets_cdi_list[(twoJp - twoJp_min) / 2];
                    auto &recouple_abj = recouple_abj_list[(twoJp - twoJp_min) / 2];
                    auto &recouple_cdi = recouple_cdi_list[(twoJp - twoJp_min) / 2];
                    double xabjcdi = 0;
                    double ycdiabj = 0;
                    for (size_t Iabj = 0; Iabj < kets_abj.size(); Iabj++)
                    {
                      for (size_t Icdi = 0; Icdi < kets_cdi.size(); Icdi++)
                      {
                        xabjcdi += recouple_abj[Iabj] * recouple_cdi[Icdi] * X3.GetME_pn_ch(ch_abj, ch_abj, kets_abj[Iabj], kets_cdi[Icdi]);
                        ycdiabj += recouple_abj[Iabj] * recouple_cdi[Icdi] * Y3.GetME_pn_ch(ch_abj, ch_abj, kets_cdi[Icdi], kets_abj[Iabj]);
                      }
                    }
                    X3_ji(ibra_CC, ind_abcd) += sixj_ji * xabjcdi;
                    Y3_ji(ind_abcd, ibra_CC) += sixj_ji * ycdiabj * occupation_factor * symmetry_factor;
                  }

                } // for twoJp

              } // for iket_cd
            }   // for iket_ab
          }     // for ch_cd
        }       // for ch_ab

      } // for ibra_CC

      Z_bar[ch] = X3_ij * Y3_ij;
      Z_bar_flipbra[ch] = X3_ji * Y3_ij;
      Z_bar_flipket[ch] = X3_ij * Y3_ji;
      Z_bar_flipboth[ch] = X3_ji * Y3_ji;

    } // for ch  (ph-coupled 2-body channels)

    //  Z.profiler.timer["comm332_CC_loops"] += omp_get_wtime() - t_internal;
    //  t_internal = omp_get_wtime();

    // Now transform Zbar to Z

    //   #pragma omp parallel for schedule(dynamic,1)
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch = 0; ch < nch; ch++)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      size_t nkets = tbc.GetNumberKets();
      for (size_t ibra = 0; ibra < nkets; ibra++)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        int phase_ij = Z.modelspace->phase((oi.j2 + oj.j2) / 2 - J);
        for (size_t iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          double jk = 0.5 * ok.j2;
          double jl = 0.5 * ol.j2;
          int phase_kl = Z.modelspace->phase((ok.j2 + ol.j2) / 2 - J);

          double zijkl = 0;

          int parity_cc = (oi.l + ol.l) % 2;
          int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
          int J_ph_min = std::max(std::abs(oi.j2 - ol.j2), std::abs(ok.j2 - oj.j2)) / 2;
          int J_ph_max = std::min((oi.j2 + ol.j2), (ok.j2 + oj.j2)) / 2;
          for (int J_ph = J_ph_min; J_ph <= J_ph_max; J_ph++)
          {
            int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(J_ph, parity_cc, Tz_cc);
            TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
            int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
            int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
            double sixj = Z.modelspace->GetSixJ(ji, jj, J, jk, jl, J_ph);
            double zbar_ilkj;
            double zbar_jkli;
            if (i <= l and k <= j) // il kj
            {
              zbar_ilkj = Z_bar[ch_cc](indx_il, indx_kj);
              zbar_jkli = Z_bar_flipboth[ch_cc](indx_kj, indx_il);
            }
            else if (i > l and k <= j) /// li kj
            {
              zbar_ilkj = Z_bar_flipbra[ch_cc](indx_il, indx_kj);
              zbar_jkli = Z_bar_flipbra[ch_cc](indx_kj, indx_il);
            }
            else if (i <= l and k > j) // il jk
            {
              zbar_ilkj = Z_bar_flipket[ch_cc](indx_il, indx_kj);
              zbar_jkli = Z_bar_flipket[ch_cc](indx_kj, indx_il);
            }
            else // li jk
            {
              zbar_ilkj = Z_bar_flipboth[ch_cc](indx_il, indx_kj);
              zbar_jkli = Z_bar[ch_cc](indx_kj, indx_il);
            }

            zijkl += (2 * J_ph + 1) * sixj * (zbar_ilkj + phase_ij * phase_kl * zbar_jkli);
          }

          parity_cc = (oj.l + ol.l) % 2;
          Tz_cc = std::abs(oj.tz2 - ol.tz2) / 2;
          J_ph_min = std::max(std::abs(oj.j2 - ol.j2), std::abs(ok.j2 - oi.j2)) / 2;
          J_ph_max = std::min((oj.j2 + ol.j2), (ok.j2 + oi.j2)) / 2;
          for (int J_ph = J_ph_min; J_ph <= J_ph_max; J_ph++)
          {
            int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(J_ph, parity_cc, Tz_cc);
            TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
            int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
            int indx_ki = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
            double sixj = Z.modelspace->GetSixJ(jj, ji, J, jk, jl, J_ph);
            double zbar_jlki;
            double zbar_iklj;

            if (j <= l and k <= i)
            {
              zbar_jlki = Z_bar[ch_cc](indx_jl, indx_ki); // <-- this one
              zbar_iklj = Z_bar_flipboth[ch_cc](indx_ki, indx_jl);
            }
            else if (j <= l and k > i)
            {
              zbar_jlki = Z_bar_flipket[ch_cc](indx_jl, indx_ki);
              zbar_iklj = Z_bar_flipket[ch_cc](indx_ki, indx_jl);
            }
            else if (j > l and k <= i)
            {
              zbar_jlki = Z_bar_flipbra[ch_cc](indx_jl, indx_ki);
              zbar_iklj = Z_bar_flipbra[ch_cc](indx_ki, indx_jl);
            }
            else if (j > l and k > i)
            {
              zbar_jlki = Z_bar_flipboth[ch_cc](indx_jl, indx_ki);
              zbar_iklj = Z_bar[ch_cc](indx_ki, indx_jl);
            }

            zijkl -= (2 * J_ph + 1) * sixj * (phase_ij * zbar_jlki + phase_kl * zbar_iklj);
          }

          // make it a normalized TBME
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          // the AddToTBME routine automatically takes care of the hermitian conjugate as well
          Z2.AddToTBME(ch, ch, ibra, iket, zijkl);

        } // for iket
      }   // for ibra
    }     // for ch

    //  Z.profiler.timer["comm332_main_loop"] += omp_get_wtime() - t_internal;

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  void comm332_pphhss_debug(const Operator &X, const Operator &Y, Operator &Z)
  // void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z )
  {
    //  std::cout << "================== BEGIN " << __func__ << " ===================" << std::endl;
    double tstart = omp_get_wtime();
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    bool x3_allocated = X3.IsAllocated();
    bool y3_allocated = Y3.IsAllocated();
    if (not(x3_allocated and y3_allocated))
      return;

    int nch = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (int ch = 0; ch < nch; ch++)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nkets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nkets; ibra++)
      {
        Ket &bra = tbc.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        int ji2 = bra.op->j2;
        int jj2 = bra.oq->j2;
        double ji = 0.5 * ji2;
        double jj = 0.5 * jj2;
        double d_ei = std::abs(2 * bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
        double d_ej = std::abs(2 * bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
        double occnat_i = bra.op->occ_nat;
        double occnat_j = bra.oq->occ_nat;
        for (int iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc.GetKet(iket);
          int k = ket.p;
          int l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          int jk2 = ket.op->j2;
          int jl2 = ket.oq->j2;
          double jk = 0.5 * jk2;
          double jl = 0.5 * jl2;
          double d_ek = std::abs(2 * ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
          double d_el = std::abs(2 * ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
          double occnat_k = ket.op->occ_nat;
          double occnat_l = ket.oq->occ_nat;

          int phase_ij = Z.modelspace->phase((ji2 + jj2) / 2 - J);
          int phase_kl = Z.modelspace->phase((jk2 + jl2) / 2 - J);
          double zijkl = 0;

          int J_ph_min = std::min({std::abs(ji2 - jl2), std::abs(jj2 - jk2), std::abs(jj2 - jl2), std::abs(ji2 - jk2)}) / 2;
          int J_ph_max = std::max({ji2 + jl2, jj2 + jk2, jj2 + jl2, ji2 + jk2}) / 2;
          //        int J_ph_min = std::max( std::abs(ji2 - jl2), std::abs(jj2-jk2) )/2;
          //        int J_ph_max = std::min( ji2+jl2, jj2+jk2 )/2;

          for (int J_ph = J_ph_min; J_ph <= J_ph_max; J_ph++)
          {
            double zbar_ilkj = 0;
            double zbar_jlki = 0;
            double zbar_iklj = 0;
            double zbar_jkli = 0;

            for (int ch_ab = 0; ch_ab < nch; ch_ab++)
            {
              TwoBodyChannel &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
              int Jab = tbc_ab.J;
              int nkets_ab = tbc_ab.GetNumberKets();

              for (int ch_cd = 0; ch_cd < nch; ch_cd++)
              {
                TwoBodyChannel &tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
                int Jcd = tbc_cd.J;
                int nkets_cd = tbc_cd.GetNumberKets();
                if (std::abs(Jab - Jcd) > J_ph or (Jab + Jcd) < J_ph)
                  continue;

                for (int iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
                {
                  Ket &ket_ab = tbc_ab.GetKet(iket_ab);
                  int a = ket_ab.p;
                  int b = ket_ab.q;
                  double na = ket_ab.op->occ;
                  double nb = ket_ab.oq->occ;
                  if (na * nb < 1e-8 and (1 - na) * (1 - nb) < 1e-8)
                    continue;
                  double d_ea = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
                  double d_eb = std::abs(2 * ket_ab.oq->n + ket_ab.oq->l - e_fermi[ket_ab.oq->tz2]);
                  double occnat_a = ket_ab.op->occ_nat;
                  double occnat_b = ket_ab.oq->occ_nat;
                  bool keep_abi = ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_i * (1 - occnat_i)) >= Z.modelspace->GetOccNat3Cut());
                  bool keep_abj = ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) >= Z.modelspace->GetOccNat3Cut());
                  bool keep_abk = ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_k * (1 - occnat_k)) >= Z.modelspace->GetOccNat3Cut());
                  bool keep_abl = ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_l * (1 - occnat_l)) >= Z.modelspace->GetOccNat3Cut());
                  keep_abi = keep_abi and (d_ea + d_eb + d_ei <= Z.modelspace->dE3max);
                  keep_abj = keep_abj and (d_ea + d_eb + d_ej <= Z.modelspace->dE3max);
                  keep_abk = keep_abk and (d_ea + d_eb + d_ek <= Z.modelspace->dE3max);
                  keep_abl = keep_abl and (d_ea + d_eb + d_el <= Z.modelspace->dE3max);
                  if (not(keep_abi or keep_abj))
                    continue;
                  if (not(keep_abk or keep_abl))
                    continue;

                  for (int iket_cd = 0; iket_cd < nkets_cd; iket_cd++)
                  {
                    Ket &ket_cd = tbc_cd.GetKet(iket_cd);
                    int c = ket_cd.p;
                    int d = ket_cd.q;
                    double nc = ket_cd.op->occ;
                    double nd = ket_cd.oq->occ;
                    double d_ec = std::abs(2 * ket_cd.op->n + ket_cd.op->l - e_fermi[ket_cd.op->tz2]);
                    double d_ed = std::abs(2 * ket_cd.oq->n + ket_cd.oq->l - e_fermi[ket_cd.oq->tz2]);
                    double occnat_c = ket_cd.op->occ_nat;
                    double occnat_d = ket_cd.oq->occ_nat;
                    //                  double jc = 0.5*ket_cd.op->j2;
                    //                  double jd = 0.5*ket_cd.oq->j2;
                    double occupation_factor = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);
                    if (std::abs(occupation_factor) < 1e-6)
                      continue;

                    bool keep_cdi = ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_i * (1 - occnat_i)) >= Z.modelspace->GetOccNat3Cut());
                    bool keep_cdj = ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_j * (1 - occnat_j)) >= Z.modelspace->GetOccNat3Cut());
                    bool keep_cdk = ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_k * (1 - occnat_k)) >= Z.modelspace->GetOccNat3Cut());
                    bool keep_cdl = ((occnat_c * (1 - occnat_c) * occnat_d * (1 - occnat_d) * occnat_l * (1 - occnat_l)) >= Z.modelspace->GetOccNat3Cut());
                    keep_cdi = keep_cdi and (d_ec + d_ed + d_ei <= Z.modelspace->dE3max);
                    keep_cdj = keep_cdj and (d_ec + d_ed + d_ej <= Z.modelspace->dE3max);
                    keep_cdk = keep_cdk and (d_ec + d_ed + d_ek <= Z.modelspace->dE3max);
                    keep_cdl = keep_cdl and (d_ec + d_ed + d_el <= Z.modelspace->dE3max);
                    if (not(keep_cdk or keep_cdl))
                      continue;
                    if (not(keep_cdi or keep_cdj))
                      continue;

                    double symmetry_factor = 1; // we only sum a<=b and c<=d, so we undercount by a factor of 4, canceling the 1/4 in the formula
                    if (a == b)
                      symmetry_factor *= 0.5; // if a==b or c==d, then the permutation doesn't give a new state, so there's less undercounting
                    if (c == d)
                      symmetry_factor *= 0.5;

                    double X3_il = 0;
                    double X3_jl = 0;
                    double X3_ik = 0;
                    double X3_jk = 0;
                    double Y3_il = 0;
                    double Y3_jl = 0;
                    double Y3_ik = 0;
                    double Y3_jk = 0;

                    if (keep_abi and keep_cdl and keep_cdj and keep_abk and ((tbc_ab.parity + tbc_cd.parity + oi.l + ol.l) % 2 == X.parity) and (std::abs(2 * tbc_ab.Tz + oi.tz2 - 2 * tbc_cd.Tz - ol.tz2) == 2 * X.rank_T) and (std::abs(ji2 - jl2) <= 2 * J_ph) and ((ji2 + jl2) >= 2 * J_ph) and (std::abs(jj2 - jk2) <= 2 * J_ph) and ((jj2 + jk2) >= 2 * J_ph))
                    {
                      int twoJp_min = std::max(std::abs(2 * Jab - ji2), std::abs(2 * Jcd - jl2));
                      int twoJp_max = std::min((2 * Jab + ji2), (2 * Jcd + jl2));
                      //                    if (i==0 and j==1 and k==0 and l==1)
                      //                    {
                      ////                      std::cout << "Getting a contribution from Xabicdl = " << a << " " << b << " " << i << " " << c << " " << d << " " << l << std::endl;
                      //                      std::cout << "Getting a contribution from X = " << a << " " << b << " " << i << " " << c << " " << d << " " << l << "   " << na << " " << nb << " " << nc << " " << nd << " -> " << occupation_factor << std::endl;
                      //                    }
                      for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                      {
                        double sixj_x = Z.modelspace->GetSixJ(ji, jl, J_ph, Jcd, Jab, 0.5 * twoJp);
                        X3_il += (twoJp + 1) * sixj_x * Z.modelspace->phase((ji2 + twoJp) / 2) * X3.GetME_pn(Jab, Jcd, twoJp, a, b, i, c, d, l);
                      }

                      twoJp_min = std::max(std::abs(2 * Jcd - jj2), std::abs(2 * Jab - jk2));
                      twoJp_max = std::min((2 * Jcd + jj2), (2 * Jab + jk2));
                      for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                      {
                        double sixj_x = Z.modelspace->GetSixJ(jj, jk, J_ph, Jab, Jcd, 0.5 * twoJp);
                        Y3_jk += (twoJp + 1) * sixj_x * Z.modelspace->phase((jk2 + twoJp) / 2) * Y3.GetME_pn(Jcd, Jab, twoJp, c, d, j, a, b, k);
                      }
                    }

                    if (keep_abj and keep_cdl and keep_cdi and keep_abk and ((tbc_ab.parity + tbc_cd.parity + oj.l + ol.l) % 2 == X.parity) and (std::abs(2 * tbc_ab.Tz + oj.tz2 - 2 * tbc_cd.Tz - ol.tz2) == 2 * X.rank_T) and (std::abs(jj2 - jl2) <= 2 * J_ph) and ((jj2 + jl2) >= 2 * J_ph) and (std::abs(ji2 - jk2) <= 2 * J_ph) and ((ji2 + jk2) >= 2 * J_ph))
                    {
                      int twoJp_min = std::max(std::abs(2 * Jab - jj2), std::abs(2 * Jcd - jl2));
                      int twoJp_max = std::min((2 * Jab + jj2), (2 * Jcd + jl2));

                      for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                      {
                        double sixj_x = Z.modelspace->GetSixJ(jj, jl, J_ph, Jcd, Jab, 0.5 * twoJp);
                        X3_jl += (twoJp + 1) * sixj_x * Z.modelspace->phase((jj2 + twoJp) / 2) * X3.GetME_pn(Jab, Jcd, twoJp, a, b, j, c, d, l);
                      }
                      twoJp_min = std::max(std::abs(2 * Jcd - ji2), std::abs(2 * Jab - jk2));
                      twoJp_max = std::min((2 * Jcd + ji2), (2 * Jab + jk2));
                      for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                      {
                        double sixj_x = Z.modelspace->GetSixJ(ji, jk, J_ph, Jab, Jcd, 0.5 * twoJp);
                        Y3_ik += (twoJp + 1) * sixj_x * Z.modelspace->phase((jk2 + twoJp) / 2) * Y3.GetME_pn(Jcd, Jab, twoJp, c, d, i, a, b, k);
                      }
                      //                    if (i==0 and j==1 and k==0 and l==1)
                      //                    {
                      ////                      std::cout << "Getting a contribution from Xabjcdl = " << a << " " << b << " " << j << " " << c << " " << d << " " << l << std::endl;
                      //                      std::cout << "Getting a contribution from X = " << a << " " << b << " " << j << " " << c << " " << d << " " << l  << "   "  << na << " " << nb << " " << nc << " " << nd << " -> "<< occupation_factor << " X3_jl = " << X3_jl << "  Y3_ik = " << Y3_ik << std::endl;
                      //                    }
                    }

                    if (keep_abi and keep_cdk and keep_cdj and keep_abl and ((tbc_ab.parity + tbc_cd.parity + oi.l + ok.l) % 2 == X.parity) and (std::abs(2 * tbc_ab.Tz + oi.tz2 - 2 * tbc_cd.Tz - ok.tz2) == 2 * X.rank_T) and (std::abs(ji2 - jk2) <= 2 * J_ph) and ((ji2 + jk2) >= 2 * J_ph) and (std::abs(jj2 - jl2) <= 2 * J_ph) and ((jj2 + jl2) >= 2 * J_ph))
                    {
                      int twoJp_min = std::max(std::abs(2 * Jab - ji2), std::abs(2 * Jcd - jk2));
                      int twoJp_max = std::min((2 * Jab + ji2), (2 * Jcd + jk2));
                      //                    if (i==0 and j==1 and k==0 and l==1)
                      //                    {
                      ////                      std::cout << "Getting a contribution from Xabicdk = " << a << " " << b << " " << i << " " << c << " " << d << " " << k << std::endl;
                      //                      std::cout << "Getting a contribution from X = " << a << " " << b << " " << i << " " << c << " " << d << " " << k  << "   "  << na << " " << nb << " " << nc << " " << nd << " -> "<< occupation_factor << std::endl;
                      //                    }
                      for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                      {
                        double sixj_x = Z.modelspace->GetSixJ(ji, jk, J_ph, Jcd, Jab, 0.5 * twoJp);
                        X3_ik += (twoJp + 1) * sixj_x * Z.modelspace->phase((ji2 + twoJp) / 2) * X3.GetME_pn(Jab, Jcd, twoJp, a, b, i, c, d, k);
                      }
                      twoJp_min = std::max(std::abs(2 * Jcd - jj2), std::abs(2 * Jab - jl2));
                      twoJp_max = std::min((2 * Jcd + jj2), (2 * Jab + jl2));
                      for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                      {
                        double sixj_x = Z.modelspace->GetSixJ(jj, jl, J_ph, Jab, Jcd, 0.5 * twoJp);
                        Y3_jl += (twoJp + 1) * sixj_x * Z.modelspace->phase((jl2 + twoJp) / 2) * Y3.GetME_pn(Jcd, Jab, twoJp, c, d, j, a, b, l);
                      }
                    }

                    if (keep_abj and keep_cdk and keep_cdi and keep_abl and ((tbc_ab.parity + tbc_cd.parity + oj.l + ok.l) % 2 == X.parity) and (std::abs(2 * tbc_ab.Tz + oj.tz2 - 2 * tbc_cd.Tz - ok.tz2) == 2 * X.rank_T) and (std::abs(jj2 - jk2) <= 2 * J_ph) and ((jj2 + jk2) >= 2 * J_ph) and (std::abs(ji2 - jl2) <= 2 * J_ph) and ((ji2 + jl2) >= 2 * J_ph))
                    {
                      int twoJp_min = std::max(std::abs(2 * Jab - jj2), std::abs(2 * Jcd - jk2));
                      int twoJp_max = std::min((2 * Jab + jj2), (2 * Jcd + jk2));
                      //                    if (i==0 and j==1 and k==0 and l==1)
                      //                    {
                      ////                      std::cout << "Getting a contribution from Xabjcdk = " << a << " " << b << " " << j << " " << c << " " << d << " " << k << std::endl;
                      //                      std::cout << "Getting a contribution from X = " << a << " " << b << " " << j << " " << c << " " << d << " " << k  << "   "  << na << " " << nb << " " << nc << " " << nd << " -> "<< occupation_factor << std::endl;
                      //                    }
                      for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                      {
                        double sixj_x = Z.modelspace->GetSixJ(jj, jk, J_ph, Jcd, Jab, 0.5 * twoJp);
                        X3_jk += (twoJp + 1) * sixj_x * Z.modelspace->phase((jj2 + twoJp) / 2) * X3.GetME_pn(Jab, Jcd, twoJp, a, b, j, c, d, k);
                      }
                      twoJp_min = std::max(std::abs(2 * Jcd - ji2), std::abs(2 * Jab - jl2));
                      twoJp_max = std::min((2 * Jcd + ji2), (2 * Jab + jl2));
                      for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                      {
                        double sixj_x = Z.modelspace->GetSixJ(ji, jl, J_ph, Jab, Jcd, 0.5 * twoJp);
                        Y3_il += (twoJp + 1) * sixj_x * Z.modelspace->phase((jl2 + twoJp) / 2) * Y3.GetME_pn(Jcd, Jab, twoJp, c, d, i, a, b, l);
                      }
                    }

                    zbar_ilkj += occupation_factor * symmetry_factor * X3_il * Y3_jk;
                    zbar_jlki += occupation_factor * symmetry_factor * X3_jl * Y3_ik;
                    zbar_iklj += occupation_factor * symmetry_factor * X3_ik * Y3_jl;
                    zbar_jkli += occupation_factor * symmetry_factor * X3_jk * Y3_il;
                    //                  zbar_ilkj += occupation_factor * symmetry_factor * ( X3_il * Y3_jk + phase_ij*phase_kl*X3_jk * Y3_il);
                    //                  zbar_jlki += occupation_factor * symmetry_factor * ( X3_jl * Y3_ik + phase_ij*phase_kl*X3_ik * Y3_jl);

                  } // for iket_cd
                }   // for iket_ab
              }     // for ch_cd
            }       // for ch_ab

            double sixj1 = Z.modelspace->GetSixJ(ji, jj, J, jk, jl, J_ph);
            double sixj1_flip = Z.modelspace->GetSixJ(jj, ji, J, jk, jl, J_ph);
            //        int phase_ij = Z.modelspace->phase( (ji2+jj2)/2-J);
            //        int phase_kl = Z.modelspace->phase( (jk2+jl2)/2-J);
            zijkl += (2 * J_ph + 1) * (sixj1 * (zbar_ilkj + phase_ij * phase_kl * zbar_jkli) - sixj1_flip * (phase_ij * zbar_jlki + phase_kl * zbar_iklj));
            //        zijkl += (2*J_ph+1) * (  sixj1 * zbar_ilkj   - sixj1_flip* phase_ij*zbar_jlki );
            //        if ( i==0 and j==1 and k==0 and l==1)
            //        {
            //           std::cout << " Jph = " << J_ph << " zbar_ilkj " << zbar_ilkj << "   zbar_jkli " << zbar_jkli << "   zbar_jlki " << zbar_jlki << "   zbar_iklj " << zbar_iklj << "  ->  " << zijkl << std::endl;
            //        }

          } // for J_ph
          // make it a normalized TBME
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          // the AddToTBME routine automatically takes care of the hermitian conjugate as well
          Z2.AddToTBME(ch, ch, ibra, iket, zijkl);
        } // for iket
      }   // for ibra
    }     // for ch
    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }



  //*****************************************************************************************
  //
  // i|  j|  k|       Uncoupled expression:
  //  |   |   *--[X]     Z_ijklmn = sum_a P(ij/k) (X_ka * Y_ijalmn) - P(lm/n) (Y_ijklma * X_an)
  //  |   |   |a
  //  *~~[Y]~~*       Coupled expression:
  //  |   |   |          Z_{ijklmn}^{J1,J2,J} = sum_a P(ij/k)^{J1,J} (X_ka * Y_{ijalmn}^{J1,J2,J})
  // l|  m|  n|                                      -P(lm/n)^{J2,J} (Y_{ijklma}^{J1,J2,J} * X_an)
  //
  //  Checked with UnitTest and passed
  //  We still have disagreement for an isospin changing operator (as of Feb 8 2024).
  //

  void comm133ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    if (Commutator::verbose)
      std::cout << __func__ << std::endl;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;

    int hermX = X.IsHermitian() ? 1 : -1;
    int hermY = Y.IsHermitian() ? 1 : -1;

    bool x_channel_diag = (X.GetTRank() == 0 and X.GetParity() == 0);
    bool y_channel_diag = (Y.GetTRank() == 0 and Y.GetParity() == 0);
    bool z_channel_diag = (Z.GetTRank() == 0 and Z.GetParity() == 0);

    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    bool x3_allocated = X3.IsAllocated();
    bool y3_allocated = Y3.IsAllocated();

    Z.modelspace->PreCalculateSixJ(); // If this has already been done, this does nothing.
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();

    std::vector<std::array<size_t, 4>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      // By default, assume X is channell diagonal, so the internal channel matches the outer channel for X.
      size_t ch_internal_xy = it.first.ch_bra;
      size_t ch_internal_yx = it.first.ch_ket;
      if (not x_channel_diag) // if X is not channel-diagonal, check if Y is.
      {
        if (y_channel_diag)
        {
          size_t ch_internal_xy = it.first.ch_ket;
          size_t ch_internal_yx = it.first.ch_bra;
        }
        else
        {
          std::cout << " Uh Oh. " << __func__ << "  called with both X and Y not channel diagonal. This is not yet implemented." << std::endl;
          return;
        }
      }
      bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, ch_internal_xy, ch_internal_yx}); // (ch_bra, ch_ket,ibra)
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();

    // TODO Identify when it helps to parallelize in the outer loop, and when it's better to give the threads to the matmult step.
    // Memory-wise, it's better to let the threads do mat mult
    // On my MacBook with 8 threads, linking against OpenBLAS, letting the matmult have the threads is better.
    // On the CRC machines with up to 24 threads, linking against MKL, parallelizing at the channel level is better by a factor 10.
    //  -SRS
    //  #pragma omp parallel for schedule(dynamic,1)
    //  for (size_t ch3=0; ch3<nch3; ch3++)
    for (size_t ich_bra_ket = 0; ich_bra_ket < n_bra_ket_ch; ich_bra_ket++)
    {
      double t_internal = omp_get_wtime();

      size_t ch_bra = bra_ket_channels[ich_bra_ket][0];
      size_t ch_ket = bra_ket_channels[ich_bra_ket][1];
      size_t ch_internal_xy = bra_ket_channels[ich_bra_ket][2];
      size_t ch_internal_yx = bra_ket_channels[ich_bra_ket][3];

      auto Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch_bra);
      auto Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch_ket);
      auto Tbc_xy = Z.modelspace->GetThreeBodyChannel(ch_internal_xy);
      auto Tbc_yx = Z.modelspace->GetThreeBodyChannel(ch_internal_yx);

      int twoJ = Tbc_bra.twoJ;

      int nbras = Tbc_bra.GetNumberKets();
      int nkets = Tbc_ket.GetNumberKets();
      int nxy = Tbc_xy.GetNumberKets();
      int nyx = Tbc_yx.GetNumberKets();

      std::vector<size_t> bras_kept;
      std::vector<size_t> kets_kept;
      std::vector<size_t> xy_kept;
      std::vector<size_t> yx_kept;

      std::map<size_t, size_t> ket_kept_lookup;
      std::map<size_t, size_t> bra_kept_lookup;
      std::map<size_t, size_t> xy_kept_lookup;
      std::map<size_t, size_t> yx_kept_lookup;


      /// Use a lambda to define a function object for a repeated snippet of code.
      /// Something similar is used in other commutator expressions and eventually this should
      /// probably be moved to something like ThreeBodyStorage::IsKetValid().
      auto keep = [&](ThreeBodyChannel &Tbc, size_t iket)
      {
        Ket3 &ket = Tbc.GetKet(iket);
        double d_ei = std::abs(2 * ket.op->n + ket.op->l - e_fermi.at(ket.op->tz2));
        double d_ej = std::abs(2 * ket.oq->n + ket.oq->l - e_fermi.at(ket.oq->tz2));
        double d_ek = std::abs(2 * ket.oR->n + ket.oR->l - e_fermi.at(ket.oR->tz2));
        double occnat_i = ket.op->occ_nat;
        double occnat_j = ket.oq->occ_nat;
        double occnat_k = ket.oR->occ_nat;
        if ((d_ei + d_ej + d_ek) > Z.modelspace->GetdE3max())
          return false;
        if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
          return false;
        return true;
      };

      for (size_t iket = 0; iket < nbras; iket++)
      {
        if (not keep(Tbc_bra, iket))
          continue;
        bras_kept.push_back(iket);
        bra_kept_lookup[iket] = bras_kept.size() - 1;
      }
      for (size_t iket = 0; iket < nkets; iket++)
      {
        if (not keep(Tbc_ket, iket))
          continue;
        kets_kept.push_back(iket);
        ket_kept_lookup[iket] = kets_kept.size() - 1;
      }
      for (size_t iket = 0; iket < nxy; iket++)
      {
        if (not keep(Tbc_xy, iket))
          continue;
        xy_kept.push_back(iket);
        xy_kept_lookup[iket] = xy_kept.size() - 1;
      }
      for (size_t iket = 0; iket < nyx; iket++)
      {
        if (not keep(Tbc_yx, iket))
          continue;
        yx_kept.push_back(iket);
        yx_kept_lookup[iket] = yx_kept.size() - 1;
      }

      size_t nbras_kept = bras_kept.size();
      size_t nkets_kept = kets_kept.size();
      size_t nxy_kept = xy_kept.size();
      size_t nyx_kept = yx_kept.size();

      Z.profiler.timer["_" + std::string(__func__) + "_check_kept"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      arma::mat X1MAT_bra(nbras_kept, nxy_kept, arma::fill::zeros);
      arma::mat X1MAT_ket(nyx_kept, nkets_kept, arma::fill::zeros);
      arma::mat Y1MAT_bra(nbras_kept, nyx_kept, arma::fill::zeros);
      arma::mat Y1MAT_ket(nxy_kept, nkets_kept, arma::fill::zeros);
      arma::mat X3MAT_bra(nbras_kept, nxy_kept, arma::fill::zeros);
      arma::mat X3MAT_ket(nyx_kept, nkets_kept, arma::fill::zeros);
      arma::mat Y3MAT_bra(nbras_kept, nyx_kept, arma::fill::zeros);
      arma::mat Y3MAT_ket(nxy_kept, nkets_kept, arma::fill::zeros);
      arma::mat Z3MAT(nbras_kept, nkets_kept, arma::fill::zeros);

      Z.profiler.timer["_" + std::string(__func__) + "_allocate_matrices"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      // kept_lookup is a map   Full index => Kept index, so iter_bra.first gives the full index, and iter_bra.second is the
      // index for the 3-body state we keep in this commutator

#pragma omp parallel for schedule(guided) collapse(2)
      for (std::size_t local_ket_index = 0; local_ket_index < kets_kept.size(); local_ket_index += 1)
      {
        for (size_t local_yx_index = 0; local_yx_index < yx_kept.size(); local_yx_index++)
        {
          std::size_t iket = kets_kept[local_ket_index];
          size_t iyx = yx_kept[local_yx_index];
          if (x3_allocated)
             X3MAT_ket(local_yx_index, local_ket_index) = X3.GetME_pn_ch(ch_internal_yx, ch_ket, iyx, iket);
          if (x_channel_diag and y_channel_diag and y3_allocated)
            Y3MAT_ket(local_yx_index, local_ket_index) = Y3.GetME_pn_ch(ch_internal_yx, ch_ket, iyx, iket);
        }
      }

      // if we're not channel diagonal, we need to deal with various cases explicitly.
      if (not(x_channel_diag and y_channel_diag))
      {
      if (y3_allocated)
      {
#pragma omp parallel for schedule(guided) collapse(2)
        for (std::size_t local_ket_index = 0; local_ket_index < kets_kept.size(); local_ket_index += 1)
        {
          for (size_t local_xy_index = 0; local_xy_index < xy_kept.size(); local_xy_index++)
          {
            std::size_t iket = kets_kept[local_ket_index];
            size_t ixy = xy_kept[local_xy_index];
            Y3MAT_ket(local_xy_index, local_ket_index) = Y3.GetME_pn_ch(ch_internal_xy, ch_ket, ixy, iket);
          }
        }
      }

      if (x3_allocated)
      {
#pragma omp parallel for schedule(guided) collapse(2)
        for (std::size_t local_bra_index = 0; local_bra_index < bras_kept.size(); local_bra_index += 1)
        {
          for (size_t local_xy_index = 0; local_xy_index < xy_kept.size(); local_xy_index++)
          {
            std::size_t ibra = bras_kept[local_bra_index];
            size_t ixy = xy_kept[local_xy_index];
            X3MAT_bra(local_bra_index, local_xy_index) = X3.GetME_pn_ch(ch_bra, ch_internal_xy, ibra, ixy);
          }
        }
      }

      if (y3_allocated)
      {
#pragma omp parallel for schedule(guided) collapse(2)
        for (std::size_t local_bra_index = 0; local_bra_index < bras_kept.size(); local_bra_index += 1)
        {
          for (size_t local_yx_index = 0; local_yx_index < yx_kept.size(); local_yx_index++)
          {
            std::size_t ibra = bras_kept[local_bra_index];
            size_t iyx = yx_kept[local_yx_index];
            Y3MAT_bra(local_bra_index, local_yx_index) = Y3.GetME_pn_ch(ch_bra, ch_internal_yx, ibra, iyx);
          }
        }
      }
      }

      Z.profiler.timer["_" + std::string(__func__) + "_fill_3Bmatrices"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      //    arma::mat X1MAT_bra( nbras_kept, nxy_kept,   arma::fill::zeros);
      //    arma::mat X1MAT_ket( nyx_kept,   nkets_kept, arma::fill::zeros);
      //    arma::mat Y1MAT_bra( nbras_kept, nyx_kept,   arma::fill::zeros);
      //    arma::mat Y1MAT_ket( nxy_kept,   nkets_kept, arma::fill::zeros);

      auto check_ijk = [&](int I, int J, int K, int Jij)
      {
        if (!Z3.IsKetValid(Jij, twoJ, I, J, K))
          return false; // this checks triangle conditions and emax cuts.
        Orbit &oI = Z.modelspace->GetOrbit(I);
        Orbit &oJ = Z.modelspace->GetOrbit(J);
        Orbit &oK = Z.modelspace->GetOrbit(K);

        // Check that this passes various cuts on 3-body states
        double d_eI = std::abs(2 * oI.n + oI.l - e_fermi.at(oI.tz2));
        double d_eJ = std::abs(2 * oJ.n + oJ.l - e_fermi.at(oJ.tz2));
        double d_eK = std::abs(2 * oK.n + oK.l - e_fermi.at(oK.tz2));
        if ((d_eI + d_eJ + d_eK) > Z.modelspace->GetdE3max())
          return false;
        double nnbar_I = oI.occ_nat * (1 - oI.occ_nat);
        double nnbar_J = oJ.occ_nat * (1 - oJ.occ_nat);
        double nnbar_K = oK.occ_nat * (1 - oK.occ_nat);
        if ((nnbar_I * nnbar_J * nnbar_K) < Z.modelspace->GetOccNat3Cut())
          return false;
        return true;
      };

#pragma omp parallel for schedule(dynamic, 1)
      for (size_t local_index_bra = 0; local_index_bra < nbras_kept; local_index_bra++)
      {
        size_t ibra = bras_kept[local_index_bra];
        Ket3 &bra = Tbc_bra.GetKet(ibra);
        int Jij = bra.Jpq;

        // Below, we do the same thing for 3 different permutations: ajk  iak  ija
        // where a is summed over all states in the same one-body channel as the orbit it replaces.
        // So
        std::vector<size_t> ijk = {bra.p, bra.q, bra.r};
        for (size_t iperm = 0; iperm < 3; iperm++)
        {

          std::vector<size_t> IJK = ijk;
          Orbit &oreplace = Z.modelspace->GetOrbit(IJK[iperm]);
          for (auto a : X.GetOneBodyChannel(oreplace.l, oreplace.j2, oreplace.tz2))
          {
            IJK[iperm] = a; // so now IJK is {a,j,k}  or {i,a,k} or {i,j,a}, depending on the permutation.
            size_t I = IJK[0];
            size_t J = IJK[1];
            size_t K = IJK[2];
            if (not check_ijk(I, J, K, Jij))
              continue;

//            std::cout << " " << __FILE__ << " line " << __LINE__ << "  ijk = " << ijk[0] << " " << ijk[1] << " " << ijk[2] <<  " -> " << IJK[0] << " " << IJK[1] << " " << IJK[2] << std::endl;

            /// If IJK is not in the ordering that we use to store the matrix elements, we need to do some recoupling.
            std::vector<size_t> ket_list;
            std::vector<double> recouple_list;
            //Z3.GetKetIndex_withRecoupling(Jij, twoJ, I, J, K, ket_list, recouple_list);
            size_t ch_return = Z3.GetKetIndex_withRecoupling(Jij, twoJ, I, J, K, ket_list, recouple_list);
            if (ch_return != ch_internal_xy) ket_list.clear(); // if this ket is in the wrong channel, don't add it. (This can happen for isospin changing operators).

            for (size_t ilist = 0; ilist < ket_list.size(); ilist++)
            {
              auto iter_find = xy_kept_lookup.find(ket_list[ilist]);
              if (iter_find == xy_kept_lookup.end())
                continue;
              size_t local_index_xy = iter_find->second;
              double recouple = recouple_list[ilist];
              X1MAT_bra(local_index_bra, local_index_xy) += X1(ijk[iperm], a) * recouple;
              if (x_channel_diag and y_channel_diag)
                Y1MAT_bra(local_index_bra, local_index_xy) += Y1(ijk[iperm], a) * recouple;
            }
          } // a

          if (not(x_channel_diag and y_channel_diag))
          {
          
            for (auto a : Y.GetOneBodyChannel(oreplace.l, oreplace.j2, oreplace.tz2))
            {
              IJK[iperm] = a; // so now IJK is {a,j,k}  or {i,a,k} or {i,j,a}, depending on the permutation.
              size_t I = IJK[0];
              size_t J = IJK[1];
              size_t K = IJK[2];
              if (not check_ijk(I, J, K, Jij))
                continue;

//            std::cout << " " << __FILE__ << " line " << __LINE__ << "  ijk = " << ijk[0] << " " << ijk[1] << " " << ijk[2] <<  " -> " << IJK[0] << " " << IJK[1] << " " << IJK[2] << std::endl;

              /// If IJK is not in the ordering that we use to store the matrix elements, we need to do some recoupling.
              std::vector<size_t> ket_list;
              std::vector<double> recouple_list;
              size_t ch_return = Z3.GetKetIndex_withRecoupling(Jij, twoJ, I, J, K, ket_list, recouple_list);
              if (ch_return != ch_internal_yx) ket_list.clear(); // if this ket is in the wrong channel, don't add it. (This can happen for isospin changing operators).

              for (size_t ilist = 0; ilist < ket_list.size(); ilist++)
              {
//                std::cout << "     ilist = " << ilist << std::endl;
                auto iter_find = yx_kept_lookup.find(ket_list[ilist]);
                if (iter_find == yx_kept_lookup.end())
                  continue;
//                std::cout << "        found it. " << std::endl;
                size_t local_index_yx = iter_find->second;
                double recouple = recouple_list[ilist];
                Y1MAT_bra(local_index_bra, local_index_yx) += Y1(ijk[iperm], a) * recouple;
//                std::cout << "        adding  " << Y1(ijk[iperm],a) << " * " << recouple << "  ->  " << Y1MAT_bra(local_index_bra, local_index_yx) << std::endl;
              }
            } // a
          }

        } // iperm
      }   // for ibra

      // now do the same thing for X1MAT_ket and Y1MAT_ket
      if (not(x_channel_diag and y_channel_diag))
      {
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t local_index_ket = 0; local_index_ket < nkets_kept; local_index_ket++)
        {
          size_t iket = kets_kept[local_index_ket];
          Ket3 &ket = Tbc_ket.GetKet(iket);
          int Jij = ket.Jpq;

          // Below, we do the same thing for 3 different permutations: ajk  iak  ija
          // where a is summed over all states in the same one-body channel as the orbit it replaces.
          // So
          std::vector<size_t> ijk = {ket.p, ket.q, ket.r};
          for (size_t iperm = 0; iperm < 3; iperm++)
          {

            std::vector<size_t> IJK = ijk;
            Orbit &oreplace = Z.modelspace->GetOrbit(IJK[iperm]);
            for (auto a : Y.GetOneBodyChannel(oreplace.l, oreplace.j2, oreplace.tz2))
            {
              IJK[iperm] = a; // so now IJK is {a,j,k}  or {i,a,k} or {i,j,a}, depending on the permutation.
              size_t I = IJK[0];
              size_t J = IJK[1];
              size_t K = IJK[2];
              if (not check_ijk(I, J, K, Jij))
                continue;

              /// If IJK is not in the ordering that we use to store the matrix elements, we need to do some recoupling.
              std::vector<size_t> ket_list;
              std::vector<double> recouple_list;
              size_t ch_return = Z3.GetKetIndex_withRecoupling(Jij, twoJ, I, J, K, ket_list, recouple_list);
              if (ch_return != ch_internal_xy) ket_list.clear(); // if this ket is in the wrong channel, don't add it. (This can happen for isospin changing operators).

              for (size_t ilist = 0; ilist < ket_list.size(); ilist++)
              {
                auto iter_find = xy_kept_lookup.find(ket_list[ilist]);
                if (iter_find == xy_kept_lookup.end())
                  continue;
                size_t local_index_xy = iter_find->second;
                double recouple = recouple_list[ilist];
                Y1MAT_ket(local_index_xy, local_index_ket) += Y1(a, ijk[iperm]) * recouple;
              }
            } // a
            for (auto a : X.GetOneBodyChannel(oreplace.l, oreplace.j2, oreplace.tz2))
            {
              IJK[iperm] = a; // so now IJK is {a,j,k}  or {i,a,k} or {i,j,a}, depending on the permutation.
              size_t I = IJK[0];
              size_t J = IJK[1];
              size_t K = IJK[2];
              if (not check_ijk(I, J, K, Jij))
                continue;

              /// If IJK is not in the ordering that we use to store the matrix elements, we need to do some recoupling.
              std::vector<size_t> ket_list;
              std::vector<double> recouple_list;
              size_t ch_return = Z3.GetKetIndex_withRecoupling(Jij, twoJ, I, J, K, ket_list, recouple_list);
              if (ch_return != ch_internal_yx) ket_list.clear(); // if this ket is in the wrong channel, don't add it. (This can happen for isospin changing operators).

              for (size_t ilist = 0; ilist < ket_list.size(); ilist++)
              {
                auto iter_find = yx_kept_lookup.find(ket_list[ilist]);
                if (iter_find == yx_kept_lookup.end())
                  continue;
                size_t local_index_yx = iter_find->second;
                double recouple = recouple_list[ilist];
                X1MAT_ket(local_index_yx, local_index_ket) += X1(a, ijk[iperm]) * recouple;
              }
            } // a

          } // iperm
        }   // for iket
      }     // if channel diag

      Z.profiler.timer["_" + std::string(__func__) + "_fill_1Bmatrices"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      // Do the matrix multiplication
      //    Z3MAT = X1MAT*Y3MAT - Y1MAT*X3MAT;
      Z3MAT = X1MAT_bra * Y3MAT_ket;
      Z3MAT -= Y1MAT_bra * X3MAT_ket;

      if (z_channel_diag)
      {
        Z3MAT -= hermX * hermY * Z3MAT.t();
      }
      else
      {
        Z3MAT -= Y3MAT_bra * X1MAT_ket;
        Z3MAT += X3MAT_bra * Y1MAT_ket;
      }


      Z.profiler.timer["_" + std::string(__func__) + "_matmult"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

// unpack the result
#pragma omp parallel for schedule(dynamic, 1000) collapse(2)
      for (std::size_t local_bra_index = 0; local_bra_index < bras_kept.size(); local_bra_index += 1)
      {
        for (std::size_t local_ket_index = 0; local_ket_index < kets_kept.size(); local_ket_index += 1)
        {
          std::size_t ibra = bras_kept[local_bra_index];
          std::size_t iket = kets_kept[local_ket_index];
          if ((ch_bra == ch_ket) and (iket < ibra))
            continue;
          Z3.AddToME_pn_ch(ch_bra, ch_ket, ibra, iket, Z3MAT(local_bra_index, local_ket_index));
        }
      }

      Z.profiler.timer["_" + std::string(__func__) + "_unpack"] += omp_get_wtime() - t_internal;

    } // for ch3

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }



  //*****************************************************************************************
  //
  //  i|  j|  k|   Uncoupled expression:
  //   |~X~|   |      Z_ijklmn  =  P(ij/k)P(lm/n) sum_a  (X_ijla * Y_akmn - Y_ijla * X_akmn)
  //   |   |a  |
  //   |   |~Y~|   Coupled expression:
  //  l|  m|  n|     Z_{ijklmn}^{J1,J2,J}  = Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9
  //                 where each of those terms corresponds to a permutation from the uncoupled expression
  //
  // We cast this as mat-mult by flipping things around and recoupling:
  // i|  j|  k|           l\ |i j|
  //  |~X~|   |             \|~X~|
  //  |   |a  |     -->          |a
  //  |   |~Y~|                  |~Y~|\
  // l|  m|  n|                 m|  n| \k
  //  checked with UnitTest, and it looks good.
  //
  //// Trying the more straighforward way again
  /////////////////////////////////////////////////
  void comm223ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double tstart = omp_get_wtime();
    if (Commutator::verbose)
      std::cout << __func__ << std::endl;
    auto &Z3 = Z.ThreeBody;
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    if ((std::abs(X2.Norm() * Y2.Norm()) < 1e-6) and not Z.modelspace->scalar3b_transform_first_pass)
      return;

    Z.modelspace->PreCalculateSixJ(); // If we already did this, this does nothing.

    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    const std::array<ThreeBodyStorage::Permutation, 3> index_perms = {ThreeBodyStorage::ABC, ThreeBodyStorage::CBA, ThreeBodyStorage::ACB};

    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, static_cast<size_t>(ibra)}); // (ch_bra, ch_ket,ibra)
      }
    }

    // If we're doing perturbative triples, we can avoid allocating any 3-body structures.
    // In this case, Z3 isn't allocated, so we need to manually specify which channels to loop over.
    if (perturbative_triples)
    {
      bra_ket_channels.clear(); // just in case we somehow had things allocated, we don't want to double count.
      size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
      for (size_t ch3 = 0; ch3 < nch3; ch3++)
      {
        ThreeBodyChannel &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
        size_t nkets = Tbc.GetNumberKets();
        for (size_t ibra = 0; ibra < nkets; ibra++)
        {
          auto &bra = Tbc.GetKet(ibra); // Will we need this bra? Check here to help with load balancing in parallel loop below
          // if ( not((bra.op->cvq + bra.oq->cvq + bra.oR->cvq) == 0 or (bra.op->cvq > 0 and bra.oq->cvq > 0 and bra.oR->cvq > 0))) continue;
          double occ_ijk = (bra.op->occ) * (bra.oq->occ) * (bra.oR->occ);
          double unocc_ijk = (1 - bra.op->occ) * (1 - bra.oq->occ) * (1 - bra.oR->occ);

          if ( (occ_ijk<1e-8) and (unocc_ijk<1e-8 ) ) continue;

          bra_ket_channels.push_back({ch3, ch3, ibra}); // (ch_bra, ch_ket,ibra)
        }
      }
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();

    double Emp2 = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Emp2)
    for (size_t ibra_ket = 0; ibra_ket < n_bra_ket_ch; ibra_ket++)
    {
      size_t ch3bra = bra_ket_channels[ibra_ket][0];
      size_t ch3ket = bra_ket_channels[ibra_ket][1];
      size_t ibra = bra_ket_channels[ibra_ket][2];
      auto &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch3bra);
      auto &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch3ket);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      size_t nkets3 = Tbc_ket.GetNumberKets();
      int twoJ = Tbc_bra.twoJ; // Scalar commutator so J is the same in bra and ket channel
      double Jtot = 0.5 * twoJ;

      if (nbras3 == 0 or nkets3 == 0)
        continue; // nothing to be done here...

      /// set up permutation stuff here so we can reuse it
      std::vector<std::array<size_t, 3>> ijk;
      std::vector<int> J1p_min;
      std::vector<int> J1p_max;
      std::vector<std::vector<double>> recouple_ijk;

      auto &bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      Orbit &ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;
      double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
      double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
      double d_ek = std::abs(2 * ok.n + ok.l - e_fermi[ok.tz2]);
      double occnat_i = oi.occ_nat;
      double occnat_j = oj.occ_nat;
      double occnat_k = ok.occ_nat;

      double occ_ijk = (bra.op->occ) * (bra.oq->occ) * (bra.oR->occ);
      double unocc_ijk = (1 - bra.op->occ) * (1 - bra.oq->occ) * (1 - bra.oR->occ);

      if ((d_ei + d_ej + d_ek) > Z.modelspace->dE3max)
        continue;
      if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
        continue;

      if (perturbative_triples and  (occ_ijk<1e-8) and (unocc_ijk<1e-8) )
        continue;
      if (imsrg3_no_qqq and (oi.cvq + oj.cvq + ok.cvq) > 5)
        continue; // Need at least one core or valence particle

      int J1 = bra.Jpq;

      size_t iket_max = nkets3 - 1;
      if (ch3bra == ch3ket)
        iket_max = ibra;
      for (size_t iket = 0; iket <= iket_max; iket++)
      {
        auto &ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit &ol = Z.modelspace->GetOrbit(l);
        Orbit &om = Z.modelspace->GetOrbit(m);
        Orbit &on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        double d_el = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs(2 * om.n + om.l - e_fermi[om.tz2]);
        double d_en = std::abs(2 * on.n + on.l - e_fermi[on.tz2]);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_n = on.occ_nat;
        double occ_lmn = (ket.op->occ) * (ket.oq->occ) * (ket.oR->occ);
        double unocc_lmn = (1 - ket.op->occ) * (1 - ket.oq->occ) * (1 - ket.oR->occ);
        if ((d_el + d_em + d_en) > Z.modelspace->dE3max)
          continue;
        if ((occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_n * (1 - occnat_n)) < Z.modelspace->GetOccNat3Cut())
          continue;

        if (perturbative_triples and (std::abs(occ_ijk*unocc_lmn - unocc_ijk*occ_lmn)<1e-8) )
          continue;
        if (imsrg3_no_qqq and (ol.cvq + om.cvq + on.cvq) > 5)
          continue;
        int J2 = ket.Jpq;

        std::vector<std::array<size_t, 3>> lmn;
        std::vector<int> J2p_min;
        std::vector<int> J2p_max;
        std::vector<std::vector<double>> recouple_lmn;

        double zijklmn = 0;
        /// BEGIN THE SLOW BIT...

        for (int perm_ijk = 0; perm_ijk < 3; perm_ijk++)
        {

          size_t I1, I2, I3;
          Z3.Permute(index_perms[perm_ijk], i, j, k, I1, I2, I3);
          Orbit &o1 = Z.modelspace->GetOrbit(I1);
          Orbit &o2 = Z.modelspace->GetOrbit(I2);
          Orbit &o3 = Z.modelspace->GetOrbit(I3);
          int j1pmin = J1;
          int j1pmax = J1;
          if (index_perms[perm_ijk] != ThreeBodyStorage::ABC)
          {
            j1pmin = std::max(std::abs(o1.j2 - o2.j2), std::abs(twoJ - o3.j2)) / 2;
            j1pmax = std::min(o1.j2 + o2.j2, twoJ + o3.j2) / 2;
          }

          double j3 = 0.5 * o3.j2;

          for (int J1p = j1pmin; J1p <= j1pmax; J1p++)
          {

            double rec_ijk = Z3.RecouplingCoefficient(index_perms[perm_ijk], ji, jj, jk, J1p, J1, twoJ);
            rec_ijk *= Z3.PermutationPhase(index_perms[perm_ijk]); // do we get a fermionic minus sign?
            if (I1 == I2)
              rec_ijk *= PhysConst::SQRT2;
            int ch12 = Z.modelspace->GetTwoBodyChannelIndex(J1p, (o1.l + o2.l) % 2, (o1.tz2 + o2.tz2) / 2);

            TwoBodyChannel &tbc12 = Z.modelspace->GetTwoBodyChannel(ch12);
            size_t nkets_12 = tbc12.GetNumberKets();
            size_t ind_12 = tbc12.GetLocalIndex(std::min(I1, I2), std::max(I1, I2));
            if (ind_12 > nkets_12)
              continue;
            double phase_12 = I1 > I2 ? -Z.modelspace->phase((o1.j2 + o2.j2 - 2 * J1p) / 2) : 1;

            for (int perm_lmn = 0; perm_lmn < 3; perm_lmn++)
            {
              size_t I4, I5, I6;
              Z3.Permute(index_perms[perm_lmn], l, m, n, I4, I5, I6);
              Orbit &o4 = Z.modelspace->GetOrbit(I4);
              Orbit &o5 = Z.modelspace->GetOrbit(I5);
              Orbit &o6 = Z.modelspace->GetOrbit(I6);
              double j6 = 0.5 * o6.j2;

              for (int parity_a : {0, 1})
              {
                //              int parity_a = (o1.l+o2.l+o6.l)%2; // Need to fix this if we want to treat parity changing operators.  ... Fixed it? SRS Jan 10.
                //              if (  (o1.l+o2.l+o6.l+parity_a) != X.GetParity() and (o4.l+o6.l+o3.l+parity_a) != X.GetParity() ) continue;
                //              if (  (o1.l+o2.l+o6.l+parity_a) != Y.GetParity() and (o4.l+o6.l+o3.l+parity_a) != Y.GetParity() ) continue;
                for (int tz2a : {-1, 1})
                {
                  int dTz126a = (o1.tz2 + o2.tz2 - o6.tz2 - tz2a) / 2;
                  int dTz3a45 = (o3.tz2 + tz2a - o4.tz2 - o5.tz2) / 2;

                  bool X_126a_good = ((o1.l + o2.l + o6.l + parity_a) % 2 == X.GetParity()) and (std::abs(dTz126a) == X.GetTRank());
                  bool Y_126a_good = ((o1.l + o2.l + o6.l + parity_a) % 2 == Y.GetParity()) and (std::abs(dTz126a) == Y.GetTRank());
                  bool X_3a45_good = ((o4.l + o5.l + o3.l + parity_a) % 2 == X.GetParity()) and (std::abs(dTz3a45) == X.GetTRank());
                  bool Y_3a45_good = ((o4.l + o5.l + o3.l + parity_a) % 2 == Y.GetParity()) and (std::abs(dTz3a45) == Y.GetTRank());

                  if (not((X_126a_good and Y_3a45_good) or (Y_126a_good and X_3a45_good)))
                    continue;
                  //                if (  std::abs(dTz126a) != X.GetTRank() and std::abs( dTz3a45) != X.GetTRank() ) continue;
                  //                if (  std::abs(dTz126a) != Y.GetTRank() and std::abs( dTz3a45) != Y.GetTRank() ) continue;

                  int ch6a = Z.modelspace->GetTwoBodyChannelIndex(J1p, (o6.l + parity_a) % 2, (o6.tz2 + tz2a) / 2);
                  TwoBodyChannel &tbc6a = Z.modelspace->GetTwoBodyChannel(ch6a);
                  size_t nkets_6a = tbc6a.GetNumberKets();

                  int j2pmin = J2;
                  int j2pmax = J2;
                  if (index_perms[perm_lmn] != ThreeBodyStorage::ABC)
                  {
                    j2pmin = std::max(std::abs(o4.j2 - o5.j2), std::abs(twoJ - o6.j2)) / 2;
                    j2pmax = std::min(o4.j2 + o5.j2, twoJ + o6.j2) / 2;
                  }

                  for (int J2p = j2pmin; J2p <= j2pmax; J2p++)
                  {

                    double rec_lmn = Z3.RecouplingCoefficient(index_perms[perm_lmn], jl, jm, jn, J2p, J2, twoJ);
                    rec_lmn *= Z3.PermutationPhase(index_perms[perm_lmn]); // do we get a fermionic minus sign?
                    if (I4 == I5)
                      rec_lmn *= PhysConst::SQRT2;

                    //                 int ch2 = Z.modelspace->GetTwoBodyChannelIndex(J2p, (o4.l+o5.l)%2, (o4.tz2+o5.tz2)/2 );
                    int ch45 = Z.modelspace->GetTwoBodyChannelIndex(J2p, (o4.l + o5.l) % 2, (o4.tz2 + o5.tz2) / 2);
                    int ch3a = Z.modelspace->GetTwoBodyChannelIndex(J2p, (o3.l + parity_a) % 2, (o3.tz2 + tz2a) / 2);

                    //                 TwoBodyChannel& tbc2 = Z.modelspace->GetTwoBodyChannel(ch2);
                    TwoBodyChannel &tbc45 = Z.modelspace->GetTwoBodyChannel(ch45);
                    TwoBodyChannel &tbc3a = Z.modelspace->GetTwoBodyChannel(ch3a);
                    //                 size_t nkets_2 = tbc2.GetNumberKets();
                    size_t nkets_45 = tbc45.GetNumberKets();
                    size_t nkets_3a = tbc3a.GetNumberKets();
                    //                 size_t ind_45 = tbc2.GetLocalIndex(std::min(I4,I5),std::max(I4,I5));
                    size_t ind_45 = tbc45.GetLocalIndex(std::min(I4, I5), std::max(I4, I5));
                    //                 if (ind_45>nkets_2) continue;
                    if (ind_45 > nkets_45)
                      continue;
                    double phase_45 = I4 > I5 ? -Z.modelspace->phase((o4.j2 + o5.j2 - 2 * J2p) / 2) : 1;

                    //                 const auto XMAT2 = X.TwoBody.GetMatrix(ch2,ch2).col(ind_45);
                    //                 const auto YMAT2 = Y.TwoBody.GetMatrix(ch2,ch2).col(ind_45);

                    double hat_factor = sqrt((2 * J1p + 1) * (2 * J2p + 1));

                    int j2a_min = std::max(std::abs(o6.j2 - 2 * J1p), std::abs(o3.j2 - 2 * J2p));
                    int j2a_max = std::min(o6.j2 + 2 * J1p, o3.j2 + 2 * J2p);

                    for (auto &it_obc : Z.modelspace->OneBodyChannels)
                    {
                      int la = it_obc.first[0];
                      int j2a = it_obc.first[1];
                      if ((j2a < j2a_min) or (j2a > j2a_max))
                        continue;
                      if (la % 2 != parity_a)
                        continue;
                      if (it_obc.first[2] != tz2a)
                        continue;
                      double ja = 0.5 * j2a;

                      double sixj;
                      if (twoJ <= 2 * Z.modelspace->GetEmax() + 1)
                      {
                        sixj = Z.modelspace->GetCachedSixJ(o3.j2, twoJ, J1p, o6.j2, j2a, J2p);
                      }
                      else
                      {
                        sixj = j2a < twoJ ? Z.modelspace->GetSixJ(j6, ja, J1p, j3, 0.5 * twoJ, J2p)
                                          : Z.modelspace->GetSixJ(j6, 0.5 * twoJ, J2p, j3, ja, J1p);
                      }
                      double prefactor = rec_ijk * rec_lmn * sixj * hat_factor * phase_12 * phase_45;

                      for (size_t a : it_obc.second)
                      {

                        // This should be the straighforward but maybe less efficient way to do it.
                        double phase_6a = phase_12;
                        double phase_3a = phase_45;
                        if (I3 == a)
                          phase_3a *= PhysConst::SQRT2;
                        if (I6 == a)
                          phase_6a *= PhysConst::SQRT2;

                        double x_126a = X_126a_good ? X.TwoBody.GetTBME_norm(ch12, ch6a, I1, I2, I6, a) : 0;
                        double y_126a = Y_126a_good ? Y.TwoBody.GetTBME_norm(ch12, ch6a, I1, I2, I6, a) : 0;
                        double x_3a45 = X_3a45_good ? X.TwoBody.GetTBME_norm(ch3a, ch45, I3, a, I4, I5) : 0;
                        double y_3a45 = Y_3a45_good ? Y.TwoBody.GetTBME_norm(ch3a, ch45, I3, a, I4, I5) : 0;

                        zijklmn += prefactor * phase_6a * phase_3a * (x_126a * y_3a45 - y_126a * x_3a45);


                      } // for a
                        //                }// for j2a
                    }   // for it_obc
                  }     // for J2p
                }       // for tz2a
              }         // for parity_a

            }           // for perm_lmn
           }             // for J1p

        }               // for perm_ijk

        // If we're doing perturbative triples, we don't need to store the full 3N, we just want the energy contribution.
        // The relevant one-body matrix elements should have been put in Z.
        if (perturbative_triples)
        {

          // double occ_ijk = (bra.op->occ) * (bra.oq->occ) * (bra.oR->occ);
          // double occ_lmn = (ket.op->occ) * (ket.oq->occ) * (ket.oR->occ);
          // double unocc_ijk = (1 - bra.op->occ) * (1 - bra.oq->occ) * (1 - bra.oR->occ);
          // double unocc_lmn = (1 - ket.op->occ) * (1 - ket.oq->occ) * (1 - ket.oR->occ);

          double symm_ijk = 6;
          if (i == j and i == k)
            symm_ijk = 1;
          else if (i == j or i == k or j == k)
            symm_ijk = 3;
          double symm_lmn = 6;
          if (l == m and l == n)
            symm_lmn = 1;
          else if (l == m or l == n or m == n)
            symm_lmn = 3;

          double Eijk = Z.OneBody(i, i) + Z.OneBody(j, j) + Z.OneBody(k, k);
          double Elmn = Z.OneBody(l, l) + Z.OneBody(m, m) + Z.OneBody(n, n);
          Emp2 += 1. / 36 * symm_ijk * symm_lmn * (twoJ + 1) * zijklmn * zijklmn * (occ_ijk * unocc_lmn - occ_lmn * unocc_ijk) / (Eijk - Elmn);
        }
        else
        {
          Z3.AddToME_pn_ch(ch3bra, ch3ket, ibra, iket, zijklmn);
        }

      } // for iket
    }   // for i_bra_ket

    // If we're doing perturbative triples, store the result in the zero-body part of Z since this function returns void.
    if (perturbative_triples)
    {
      Z.ZeroBody = Emp2;
    }

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }






  bool check_2b_channel_Tz_parity(const Operator &Op, Orbit &o1, Orbit &o2, Orbit &o3, Orbit &o4)
  {
    return (((o1.l + o2.l + o3.l + o4.l) % 2 == Op.parity) and (std::abs(o1.tz2 + o2.tz2 - o3.tz2 - o4.tz2) == 2 * Op.rank_T));
  }

  //*****************************************************************************************
  //
  //  |     |    |     Uncoupled expression:
  // i|    j|   k|       Z_ijklmn = 1/2 sum_{ab} (n`an`b - nanb) P(ij/k) X_{ijab} Y_{abklmn} - P(lm/n) Y_{ijkabn} X_{ablm}
  //  *~~X~~*    |
  // a|    b|    |
  //  *~~~~[Y]~~~*     Coupled expression:
  // l|    m|   n|     Z_{ijklmn}^{J1,J2,J} = 1/2 sum_{ab} (n`an`b - nanb) ( P(ij/k)^{J1,J} X_{ijab}^{J1} Y_{abklmn}^{J1,J2,J}
  //  |     |    |                                                         - P(lm/n)^{J2,J} Y_{ijkabn}^{J1,J2,J} X_{ablm}^{J2} )
  //  |     |    |
  //
  //
  //      Checked with UnitTest and passed
  //
  void comm233_pp_hhss(const Operator &X, const Operator &Y, Operator &Z)
  {
    //  std::cout << "ENTER " <<__func__ << std::endl;
    double tstart = omp_get_wtime();
    if (Commutator::verbose)
      std::cout << __func__ << std::endl;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;

    int hermX = X.IsHermitian() ? 1 : -1;
    int hermY = Y.IsHermitian() ? 1 : -1;

    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    //  int norbs = Z.modelspace->GetNumberOrbits();
    //  double X3NORM = X3.Norm();
    //  double Y3NORM = Y3.Norm();
    //  bool x3_allocated = X3.is_allocated;
    //  bool y3_allocated = Y3.is_allocated;
    bool x3_allocated = X3.IsAllocated();
    bool y3_allocated = Y3.IsAllocated();
    //  std::cout << "Begin the ch3 loop " << std::endl;

    Z.modelspace->PreCalculateSixJ();

    size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch3 = 0; ch3 < nch3; ch3++)
    {
      auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      int twoJ = Tbc.twoJ;
      double Jtot = 0.5 * twoJ;
      size_t nkets = Tbc.GetNumberKets();
      std::vector<size_t> kets_kept;
      std::map<size_t, size_t> kept_lookup;
      for (size_t iket = 0; iket < nkets; iket++)
      {
        Ket3 &ket = Tbc.GetKet(iket);
        double d_ei = std::abs(2 * ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
        double d_ej = std::abs(2 * ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
        double d_ek = std::abs(2 * ket.oR->n + ket.oR->l - e_fermi[ket.oR->tz2]);
        double occnat_i = ket.op->occ_nat;
        double occnat_j = ket.oq->occ_nat;
        double occnat_k = ket.oR->occ_nat;
        if ((d_ei + d_ej + d_ek) > Z.modelspace->GetdE3max())
          continue;
        if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
          continue;
        kets_kept.push_back(iket);
        kept_lookup[iket] = kets_kept.size() - 1;
      }
      size_t nkets_kept = kets_kept.size();

      arma::mat X2MAT(nkets_kept, nkets_kept, arma::fill::zeros);
      arma::mat Y2MAT(nkets_kept, nkets_kept, arma::fill::zeros);
      arma::mat X3MAT(nkets_kept, nkets_kept, arma::fill::zeros);
      arma::mat Y3MAT(nkets_kept, nkets_kept, arma::fill::zeros);
      arma::mat Z3MAT(nkets_kept, nkets_kept, arma::fill::zeros);

      //    std::cout << "   fill the 2b matrices ch = " << ch3 << std::endl;

      for (size_t index_bra = 0; index_bra < nkets_kept; index_bra++)
      {
        size_t ibra = kets_kept[index_bra];
        Ket3 &bra = Tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        Orbit &ok = Z.modelspace->GetOrbit(k);
        double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
        double d_ek = std::abs(2 * ok.n + ok.l - e_fermi[ok.tz2]);
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;
        double occnat_k = ok.occ_nat;
        int Jij = bra.Jpq;
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        double jk = 0.5 * ok.j2;

        for (size_t ch_ab = 0; ch_ab < nch2; ch_ab++)
        {
          auto &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
          size_t nkets_ab = tbc_ab.GetNumberKets();
          int Jab = tbc_ab.J;

          //        std::cout << "   direct " << std::endl;
          if (((2 * tbc_ab.Tz + ok.tz2) == Tbc.twoTz) and ((tbc_ab.parity + ok.l + Tbc.parity) % 2 == 0) and Jab == Jij)
          {
            for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
            {
              auto &ket_ab = tbc_ab.GetKet(iket_ab);
              size_t a = ket_ab.p;
              size_t b = ket_ab.q;
              double d_ea = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
              double d_eb = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
              double occnat_a = ket_ab.op->occ_nat;
              double occnat_b = ket_ab.oq->occ_nat;
              if ((d_ea + d_eb + d_ek) > Z.modelspace->dE3max)
                continue;
              if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
                continue;
              double occfactor = 1 - ket_ab.op->occ - ket_ab.oq->occ;
              if (std::abs(occfactor) < 1e-6)
                continue;
              double symm_factor = (a == b) ? 1.0 : 2.0;
              double xijab = X2.GetTBME_J(Jab, Jab, i, j, a, b);
              double yijab = Y2.GetTBME_J(Jab, Jab, i, j, a, b);

              if (!Z3.IsKetValid(Jab, twoJ, a, b, k))
                continue;

              std::vector<size_t> ket_list;
              std::vector<double> recouple_list;

              //              size_t ch_check = Z3.GetKetIndex_withRecoupling( Jab, twoJ, a, b, k,  ket_list,  recouple_list );
              Z3.GetKetIndex_withRecoupling(Jab, twoJ, a, b, k, ket_list, recouple_list);
              for (size_t ilist = 0; ilist < ket_list.size(); ilist++)
              {
                auto iter_find = kept_lookup.find(ket_list[ilist]);
                if (iter_find == kept_lookup.end())
                  continue;
                size_t index_ket = iter_find->second;
                double recouple = recouple_list[ilist];
                X2MAT(index_bra, index_ket) += 0.5 * symm_factor * occfactor * recouple * xijab;
                Y2MAT(index_bra, index_ket) += 0.5 * symm_factor * occfactor * recouple * yijab;
              } // for ilist

            } // for iket_ab
          }   // if J, Tz and parity check for term 1

          //        std::cout << "   Pik " << std::endl;
          // Permute  Pik
          if (((2 * tbc_ab.Tz + oi.tz2) == Tbc.twoTz) and ((tbc_ab.parity + oi.l + Tbc.parity) % 2 == 0) and (std::abs(2 * Jab - oi.j2) <= twoJ) and ((2 * Jab + oi.j2) >= twoJ) and (std::abs(2 * Jab - oj.j2) <= ok.j2) and ((2 * Jab + oj.j2) >= ok.j2))
          {
            double sixj = Z.modelspace->GetSixJ(ji, jj, Jij, jk, Jtot, Jab);
            if (std::abs(sixj) > 1e-6)
            {
              double hats = sqrt((2 * Jij + 1) * (2 * Jab + 1));
              for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
              {
                auto &ket_ab = tbc_ab.GetKet(iket_ab);
                size_t a = ket_ab.p;
                size_t b = ket_ab.q;
                double d_ea = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
                double d_eb = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
                double occnat_a = ket_ab.op->occ_nat;
                double occnat_b = ket_ab.oq->occ_nat;
                if ((d_ea + d_eb + d_ei) > Z.modelspace->dE3max)
                  continue;
                if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_i * (1 - occnat_i)) < Z.modelspace->GetOccNat3Cut())
                  continue;
                double occfactor = 1 - ket_ab.op->occ - ket_ab.oq->occ;
                if (std::abs(occfactor) < 1e-6)
                  continue;
                double symm_factor = (a == b) ? 1.0 : 2.0;
                double xkjab = X2.GetTBME_J(Jab, Jab, k, j, a, b);
                double ykjab = Y2.GetTBME_J(Jab, Jab, k, j, a, b);

                if (!Z3.IsKetValid(Jab, twoJ, a, b, i))
                  continue;
                std::vector<size_t> ket_list;
                std::vector<double> recouple_list;

                //              std::cout << " line " << __LINE__ << "  calling GetKetIndex_withRecoupling " << Jab << " " << twoJ << " " << a << " " << b << " " << i << std::endl;
                //              size_t ch_check = Z3.GetKetIndex_withRecoupling( Jab, twoJ, a, b, i,  ket_list,  recouple_list );
                Z3.GetKetIndex_withRecoupling(Jab, twoJ, a, b, i, ket_list, recouple_list);
                //              std::cout << "    size of lists " << ket_list.size() << " " << recouple_list.size() << std::endl;
                for (size_t ilist = 0; ilist < ket_list.size(); ilist++)
                {
                  auto iter_find = kept_lookup.find(ket_list[ilist]);
                  if (iter_find == kept_lookup.end())
                    continue;
                  size_t index_ket = iter_find->second;
                  double recouple = recouple_list[ilist];
                  X2MAT(index_bra, index_ket) += 0.5 * hats * sixj * symm_factor * occfactor * recouple * xkjab;
                  Y2MAT(index_bra, index_ket) += 0.5 * hats * sixj * symm_factor * occfactor * recouple * ykjab;
                } // for ilist
                  //              std::cout << "  ok. line " << __LINE__ << std::endl;

              } // for iket_ab
            }   // if sixj nonzero
          }     // if J, Tz and parity check for term 2

          //        std::cout << "   Pjk " << std::endl;
          // Permute  Pjk
          if (((2 * tbc_ab.Tz + oj.tz2) == Tbc.twoTz) and ((tbc_ab.parity + oj.l + Tbc.parity) % 2 == 0) and (std::abs(2 * Jab - oj.j2) <= twoJ) and ((2 * Jab + oj.j2) >= twoJ) and (std::abs(2 * Jab - oi.j2) <= ok.j2) and ((2 * Jab + oi.j2) >= ok.j2))
          {
            double sixj = Z.modelspace->GetSixJ(jj, ji, Jij, jk, Jtot, Jab);
            if (std::abs(sixj) > 1e-6)
            {
              double hats = sqrt((2 * Jij + 1) * (2 * Jab + 1));
              int phase = Z.modelspace->phase((oj.j2 + ok.j2) / 2 - Jij - Jab);
              for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
              {
                auto &ket_ab = tbc_ab.GetKet(iket_ab);
                size_t a = ket_ab.p;
                size_t b = ket_ab.q;
                double d_ea = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
                double d_eb = std::abs(2 * ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2]);
                double occnat_a = ket_ab.op->occ_nat;
                double occnat_b = ket_ab.oq->occ_nat;
                if ((d_ea + d_eb + d_ej) > Z.modelspace->dE3max)
                  continue;
                if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) < Z.modelspace->GetOccNat3Cut())
                  continue;
                double occfactor = 1 - ket_ab.op->occ - ket_ab.oq->occ;
                if (std::abs(occfactor) < 1e-6)
                  continue;
                double symm_factor = (a == b) ? 1.0 : 2.0;
                double xikab = X2.GetTBME_J(Jab, Jab, i, k, a, b);
                double yikab = Y2.GetTBME_J(Jab, Jab, i, k, a, b);

                if (!Z3.IsKetValid(Jab, twoJ, a, b, j))
                  continue;

                std::vector<size_t> ket_list;
                std::vector<double> recouple_list;

                //              size_t ch_check = Z3.GetKetIndex_withRecoupling( Jab, twoJ, a, b, j,  ket_list,  recouple_list );
                Z3.GetKetIndex_withRecoupling(Jab, twoJ, a, b, j, ket_list, recouple_list);
                for (size_t ilist = 0; ilist < ket_list.size(); ilist++)
                {
                  auto iter_find = kept_lookup.find(ket_list[ilist]);
                  if (iter_find == kept_lookup.end())
                    continue;
                  size_t index_ket = iter_find->second;
                  double recouple = recouple_list[ilist];
                  X2MAT(index_bra, index_ket) -= 0.5 * phase * hats * sixj * symm_factor * occfactor * recouple * xikab;
                  Y2MAT(index_bra, index_ket) -= 0.5 * phase * hats * sixj * symm_factor * occfactor * recouple * yikab;

                } // for ilist

              } // for iket_ab
            }   // if sixj nonzero
          }     // if  J, Tz and parity check for term 3

        } // for ch_ab

      } // for index_bra

      //     std::cout << "Fill the 3b matrices " << std::endl;

      // Fill X3 and Y3
      // kept_lookup is a map   Full index => Kept index, so iter_bra.first gives the full index, and iter_bra.second is the
      // index for the 3-body state we keep in this commutator
      for (auto &iter_bra : kept_lookup)
      {
        for (auto &iter_ket : kept_lookup)
        {
          if (x3_allocated)
            X3MAT(iter_bra.second, iter_ket.second) = X3.GetME_pn_ch(ch3, ch3, iter_bra.first, iter_ket.first);

          if (y3_allocated)
            Y3MAT(iter_bra.second, iter_ket.second) = Y3.GetME_pn_ch(ch3, ch3, iter_bra.first, iter_ket.first);
        }
      }

      // Do the matrix multiplication
      //    Z3MAT = X2MAT*Y3MAT - Y2MAT*X3MAT  +  hermX*hermY * ( X3MAT.t()*Y2MAT.t() - Y3MAT.t()*X2MAT.t() );
      Z3MAT = X2MAT * Y3MAT - Y2MAT * X3MAT;
      Z3MAT -= hermX * hermY * Z3MAT.t();

      //    std::cout << "Unpack..." << std::endl;

      // unpack the result
      for (auto &iter_bra : kept_lookup)
      {
        for (auto &iter_ket : kept_lookup)
        {
          if (iter_ket.first > iter_bra.first)
            continue;
          Z3.AddToME_pn_ch(ch3, ch3, iter_bra.first, iter_ket.first, Z3MAT(iter_bra.second, iter_ket.second));
        }
      }

    } // for ch3
    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  // old much slower way
  /*

  void comm233_pp_hhss_debug( const Operator& X, const Operator& Y, Operator& Z )
  {

    double tstart = omp_get_wtime();
    auto& X2 = X.TwoBody;
    auto& Y2 = Y.TwoBody;
    auto& X3 = X.ThreeBody;
    auto& Y3 = Y.ThreeBody;
    auto& Z3 = Z.ThreeBody;
    std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
    size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
    #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch=0; ch<=nch3; ch++)
    {
      auto Tbc = Z.modelspace->GetThreeBodyChannel(ch);
      size_t nket3 = Tbc.GetNumberKets();
      int twoJ = Tbc.twoJ;
      double Jtot = 0.5 * twoJ;
      for (size_t ibra=0; ibra<nket3; ibra++)
      {
        auto& bra = Tbc.GetKet(ibra);
        int J1 = bra.Jpq;
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit& oi = Z.modelspace->GetOrbit(i);
        Orbit& oj = Z.modelspace->GetOrbit(j);
        Orbit& ok = Z.modelspace->GetOrbit(k);
        double ji = 0.5*oi.j2;
        double jj = 0.5*oj.j2;
        double jk = 0.5*ok.j2;
        double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
        double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
        if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
        for (size_t iket=0; iket<=ibra; iket++)
        {
          auto& ket = Tbc.GetKet(iket);
          int J2 = ket.Jpq;
          int l = ket.p;
          int m = ket.q;
          int n = ket.r;
          Orbit& ol = Z.modelspace->GetOrbit(l);
          Orbit& om = Z.modelspace->GetOrbit(m);
          Orbit& on = Z.modelspace->GetOrbit(n);
          double jl = 0.5*ol.j2;
          double jm = 0.5*om.j2;
          double jn = 0.5*on.j2;
          double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
          double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
          double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
          if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;

          double z_ijklmn = 0;

          for (size_t ch_ab=0; ch_ab<nch2; ch_ab++)
          {
            auto& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
            size_t nkets_ab = tbc_ab.GetNumberKets();
            int Jab = tbc_ab.J;

            if ( ((2*tbc_ab.Tz + ok.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + ok.l + Tbc.parity)%2==0) and Jab == J1)
            {
              for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
              {
                auto& ket_ab = tbc_ab.GetKet(iket_ab);
                size_t a = ket_ab.p;
                size_t b = ket_ab.q;
                double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                if ( (d_ea + d_eb + d_ek) > Z.modelspace->dE3max ) continue;
                double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                if (std::abs(occfactor)<1e-6) continue;
                double symm_factor =  (a==b) ? 1.0 : 2.0;
                double xijab = X2.GetTBME_J(Jab,Jab,i,j,a,b);
                double yijab = Y2.GetTBME_J(Jab,Jab,i,j,a,b);
                double xabklmn = X3.GetME_pn(J1,J2,twoJ,a,b,k,l,m,n);
                double yabklmn = Y3.GetME_pn(J1,J2,twoJ,a,b,k,l,m,n);
                z_ijklmn += 0.5 * symm_factor * occfactor * (xijab * yabklmn - yijab * xabklmn);
  //              if (ch==0 and ibra==4 and iket==4)
  //              {
  //                 std::cout << "  iket_ab = " << iket_ab << "   adding to z   0.5 * " << symm_factor << " * " << occfactor << " * " << xijab << " " << yabklmn << "  =>  " << z_ijklmn << std::endl;;
  //              }
              }// for iket_ab
            } // if J, Tz and parity check for term 1


            // Permute  Pik
            if ( ((2*tbc_ab.Tz + oi.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + oi.l + Tbc.parity)%2==0)
                       and (std::abs(2*Jab-oi.j2)<=twoJ) and ((2*Jab+oi.j2)>=twoJ) )
            {
               double sixj = Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,Jab);
               if (std::abs(sixj)>1e-6)
               {
                double hats = sqrt( (2*J1+1)*(2*Jab+1));
                for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
                {
                  auto& ket_ab = tbc_ab.GetKet(iket_ab);
                  size_t a = ket_ab.p;
                  size_t b = ket_ab.q;
                  double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                  double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                  if ( (d_ea + d_eb + d_ei) > Z.modelspace->dE3max ) continue;
                  double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                  if (std::abs(occfactor)<1e-6) continue;
                  double symm_factor =  (a==b) ? 1.0 : 2.0;
                  double xkjab = X2.GetTBME_J(Jab,Jab,k,j,a,b);
                  double ykjab = Y2.GetTBME_J(Jab,Jab,k,j,a,b);
                  double xabilmn = X3.GetME_pn(Jab,J2,twoJ,a,b,i,l,m,n);
                  double yabilmn = Y3.GetME_pn(Jab,J2,twoJ,a,b,i,l,m,n);
                  z_ijklmn += 0.5  * symm_factor * occfactor * hats * sixj * (xkjab * yabilmn - ykjab * xabilmn);
  //              if (ch==0 and ibra==4 and iket==4)
  //              {
  //                 std::cout << "Pik  iket_ab = " << iket_ab << "   adding to z   0.5 * " << symm_factor << " * " << occfactor << " * " << hats
  //                           << " * " << sixj << " * " << xkjab << " " << yabilmn << "  =>  " << z_ijklmn << std::endl;;
  //              }
                }// for iket_ab
              }// if sixj nonzero
            } // if J, Tz and parity check for term 2

            // Permute  Pjk
            if ( ((2*tbc_ab.Tz + oj.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + oj.l + Tbc.parity)%2==0)
                       and (std::abs(2*Jab-oj.j2)<=twoJ) and ((2*Jab+oj.j2)>=twoJ) )
            {
               double sixj = Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,Jab);
               if (std::abs(sixj)>1e-6)
               {
                 double hats = sqrt( (2*J1+1)*(2*Jab+1));
                 int phase = Z.modelspace->phase( (oj.j2+ok.j2)/2-J1-Jab);
                 for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
                 {
                   auto& ket_ab = tbc_ab.GetKet(iket_ab);
                   size_t a = ket_ab.p;
                   size_t b = ket_ab.q;
                   double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                   double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                   if ( (d_ea + d_eb + d_ej) > Z.modelspace->dE3max ) continue;
                   double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                   if (std::abs(occfactor)<1e-6) continue;
                   double symm_factor =  (a==b) ? 1.0 : 2.0;
                   double xikab = X2.GetTBME_J(Jab,Jab,i,k,a,b);
                   double yikab = Y2.GetTBME_J(Jab,Jab,i,k,a,b);
                   double xabjlmn = X3.GetME_pn(Jab,J2,twoJ,a,b,j,l,m,n);
                   double yabjlmn = Y3.GetME_pn(Jab,J2,twoJ,a,b,j,l,m,n);
                   z_ijklmn -= 0.5 * symm_factor  * occfactor * phase * hats * sixj * (xikab * yabjlmn - yikab * xabjlmn);
  //              if (ch==0 and ibra==4 and iket==4)
  //              {
  //                 std::cout << "Pjk  iket_ab = " << iket_ab << "   adding to z   0.5 * " << symm_factor << " * " << occfactor << " * "
  //                           << phase << " * " << hats << " * " << sixj << " * " << xikab << " " << yabjlmn << "  =>  " << z_ijklmn << std::endl;
  //              }
                 }// for iket_ab
              }// if sixj nonzero
            } // if  J, Tz and parity check for term 3

            // direct Y3 X2 term
            if ( ((2*tbc_ab.Tz + on.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + on.l + Tbc.parity)%2==0) and Jab == J2 )
            {
              for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
              {
                auto& ket_ab = tbc_ab.GetKet(iket_ab);
                size_t a = ket_ab.p;
                size_t b = ket_ab.q;
                double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                if ( (d_ea + d_eb + d_en) > Z.modelspace->dE3max ) continue;
                double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                if (std::abs(occfactor)<1e-6) continue;
                double symm_factor =  (a==b) ? 1.0 : 2.0;
                double xablm = X2.GetTBME_J(Jab,Jab,a,b,l,m);
                double yablm = Y2.GetTBME_J(Jab,Jab,a,b,l,m);
                double xijkabn = X3.GetME_pn(J1,J2,twoJ,i,j,k,a,b,n);
                double yijkabn = Y3.GetME_pn(J1,J2,twoJ,i,j,k,a,b,n);
                z_ijklmn -= 0.5 * symm_factor  * occfactor * (xablm * yijkabn - yablm * xijkabn);
              }// for iket_ab
            } // if  Tz and parity check for term 4


            // Permute  Pln
            if ( ((2*tbc_ab.Tz + ol.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + ol.l + Tbc.parity)%2==0)
                       and (std::abs(2*Jab-ol.j2)<=twoJ) and ((2*Jab+ol.j2)>=twoJ) )
            {
               double sixj = Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,Jab);
               if (std::abs(sixj)>1e-6)
               {
                 double hats = sqrt( (2*J2+1)*(2*Jab+1));
                 for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
                 {
                   auto& ket_ab = tbc_ab.GetKet(iket_ab);
                   size_t a = ket_ab.p;
                   size_t b = ket_ab.q;
                   double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                   double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                   if ( (d_ea + d_eb + d_el) > Z.modelspace->dE3max ) continue;
                   double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                   if (std::abs(occfactor)<1e-6) continue;
                   double symm_factor =  (a==b) ? 1.0 : 2.0;
                   double xabnm = X2.GetTBME_J(Jab,Jab,a,b,n,m);
                   double yabnm = Y2.GetTBME_J(Jab,Jab,a,b,n,m);
                   double xijkabl = X3.GetME_pn(J1,Jab,twoJ,i,j,k,a,b,l);
                   double yijkabl = Y3.GetME_pn(J1,Jab,twoJ,i,j,k,a,b,l);
                   z_ijklmn -= 0.5 * symm_factor  * occfactor * hats * sixj * (xabnm * yijkabl - yabnm * xijkabl);
                 }// for iket_ab
               }// if sixj nonzero
            } // if  Tz and parity check for term 5


            // Permute  Pmn
            if ( ((2*tbc_ab.Tz + om.tz2) == Tbc.twoTz)  and ( (tbc_ab.parity + om.l + Tbc.parity)%2==0)
                       and (std::abs(2*Jab-om.j2)<=twoJ) and ((2*Jab+om.j2)>=twoJ) )
            {
                double sixj = Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,Jab);
               if (std::abs(sixj)>1e-6)
               {
                double hats = sqrt( (2*J2+1)*(2*Jab+1));
                int phase = Z.modelspace->phase( (om.j2+on.j2)/2-J2-Jab);
                for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
                {
                  auto& ket_ab = tbc_ab.GetKet(iket_ab);
                  size_t a = ket_ab.p;
                  size_t b = ket_ab.q;
                  double d_ea = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                  double d_eb = std::abs( 2*ket_ab.op->n + ket_ab.op->l - e_fermi[ket_ab.op->tz2] );
                  if ( (d_ea + d_eb + d_em) > Z.modelspace->dE3max ) continue;
                  double occfactor = 1-ket_ab.op->occ-ket_ab.oq->occ;
                  if (std::abs(occfactor)<1e-6) continue;
                  double symm_factor =  (a==b) ? 1.0 : 2.0;
                  double xabln = X2.GetTBME_J(Jab,Jab,a,b,l,n);
                  double yabln = Y2.GetTBME_J(Jab,Jab,a,b,l,n);
                  double xijkabm = X3.GetME_pn(J1,Jab,twoJ,i,j,k,a,b,m);
                  double yijkabm = Y3.GetME_pn(J1,Jab,twoJ,i,j,k,a,b,m);
                  z_ijklmn += 0.5 * symm_factor * phase  * occfactor * hats * sixj * (xabln * yijkabm - yabln * xijkabm);
                }// for iket_ab
              }// if sixj nonzero
            } // if  Tz and parity check for term 6

          }// for ch_ab

  //        if (ch==0 and ibra<5 and iket<5)
  //        {
  //          std::cout << __func__ << "   ch = " << ch << "  ibra =  " <<  ibra << "  iket = " << iket
  //                    << "   ijklmn " << bra.p << " " << bra.q << " " << bra.r << " " << ket.p << " " << ket.q << " " << ket.r
  //                    << "  Jij Jlm twoJ " << bra.Jpq << " " << ket.Jpq << "  " << twoJ
  //                    << "  zijklmn = " << z_ijklmn << std::endl;
  //        }
          // no normalization to be done. we store un-normalized 3body matrix elements.
          if (std::abs(z_ijklmn)>1e-6)
          {
  //            Z3.AddToME_pn( J1, J2, twoJ, i,j,k,l,m,n, z_ijklmn );
              Z3.AddToME_pn_ch( ch,ch,ibra,iket, z_ijklmn );
          }

        }// for iket
      }// for ibra
    }// for ch

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }
  */

  //*****************************************************************************************
  //
  // i|   j|        k|  Uncoupled expression:
  //  |    |         |    C_ijklmn = sum_{ab} (n`anb -nan`b) P(ij/k)P(lm/n) Y_{ijalmb} X_{bkan}
  //  |    |   *~~X~~|
  //  |    | a/ \b   |
  //  |    |  \ /    |
  //  |~~~[Y]~~*     |  Coupled expression:
  //  |    |         |  C_{ijklmn}^{J1,J2,J} =  P^{J1,J}(ij/k)P^{J2,j}_{lm/n) sum_{ab} (n`anb -nan`b)
  // l|   m|        n|                          * sum_{J',J3} (2J'+1)(2J3+1) (-1)^{2J+k-n+J1-J2}
  //  |    |         |                          { J   J2  n  }
  //                                          * { J1  J'  a  } * Y_{ijalmb}^{J1,J2,J'} X_{bkan}^{J3}
  //                                            { k   b   J3 }
  //
  //      Checked with UnitTest and passed
  //
  //  This gets flipped around to look like
  //
  //                     k|  /n
  //                      | /
  //                *~~X~~*
  //              a/ \b
  //               \ /
  //    ,*~~~~[Y]~~*
  //   / |   / |
  // i/ l|  /j |m
  //
  void comm233_phss(const Operator &X, const Operator &Y, Operator &Z)
  {

    double tstart = omp_get_wtime();
    if (Commutator::verbose)
      std::cout << __func__ << std::endl;
    double t_internal = omp_get_wtime();
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();
    int hX = X.IsHermitian() ? 1 : -1;
    int hY = Y.IsHermitian() ? 1 : -1;

    bool x3_allocated = X3.IsAllocated();
    bool y3_allocated = Y3.IsAllocated();

    size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();          // number of regular 2b channels
    size_t nch2_CC = Z.modelspace->TwoBodyChannels_CC.size();        // number of particle hole channels
    std::deque<std::vector<size_t>> ph_kets_cc(nch2_CC);             // list of ph 2b kets for each CC channel
    std::deque<std::map<size_t, size_t>> ph_kets_cc_lookup(nch2_CC); // map i j1 to  matrix index
    std::deque<std::vector<double>> occfactors(nch2_CC);             // occupation factors for ab`

    auto hash_key_lmij = [](size_t l, size_t m, size_t Jlm, size_t i, size_t j, size_t Jij)
    { return (l + (m << 12) + (i << 24) + (j << 36) + (Jlm << 48) + (Jij << 56)); };
    auto unhash_key_lmij = [](size_t &l, size_t &m, size_t &Jlm, size_t &i, size_t &j, size_t &Jij, size_t key)
    { l=(key & 0xFFFL);m=((key>>12)&0xFFFL);i=((key>>24)&0xFFFL);j=((key>>36)&0xFFFL); Jlm=((key>>48)&0xFFL); Jij=((key>>56)&0xFFL); };                                                                     //  0xF = 15 = 1111 (4bits), so 0xFFF is 12 bits of 1's. 0xFFFL makes it a long
    std::deque<std::unordered_map<size_t, size_t>> ph_kets3_lookup(nch2_CC); // map  lm Jlm  (ij Jij)` to matrix index
                                                                             //  std::deque<std::map<std::array<int,6>,size_t>> ph_kets3_lookup(nch2_CC); // map  lm Jlm  (ij Jij)` to matrix index

    std::deque<arma::mat> Z3_ph(nch2_CC); // matrices for ph transformed Z3. This will be the big one.

    double occnat_factor_max = 0;
    for (auto i : Z.modelspace->all_orbits)
    {
      double occnat_i = Z.modelspace->GetOrbit(i).occ_nat;
      occnat_factor_max = std::max(occnat_factor_max, occnat_i * (1 - occnat_i));
    }

    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch_cc = 0; ch_cc < nch2_CC; ch_cc++)
    {
      TwoBodyChannel_CC &tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      //    int Jph = tbc_CC.J;
      size_t nkets_CC = tbc_CC.GetNumberKets();
      // Figure out which 2b ph kets will contribute to ab`
      for (size_t iket = 0; iket < nkets_CC; iket++)
      {
        Ket &ket = tbc_CC.GetKet(iket);
        double occfactor = ket.op->occ - ket.oq->occ;
        if (std::abs(occfactor) > 1e-7)
        {
          ph_kets_cc[ch_cc].push_back(iket);
          ph_kets_cc_lookup[ch_cc][iket] = ph_kets_cc[ch_cc].size() - 1;
          occfactors[ch_cc].push_back(occfactor);
        }
      } // for iket

      //    size_t nkets_ph = ph_kets_cc[ch_cc].size();

      auto &ph_kets3 = ph_kets3_lookup[ch_cc];
      size_t nkets3 = 0;
      // next, figure out how many states we'll need for the 3b matrices
      for (size_t ch_lm = 0; ch_lm < nch2; ch_lm++)
      {
        TwoBodyChannel &tbc_lm = Z.modelspace->GetTwoBodyChannel(ch_lm);
        size_t nkets_lm = tbc_lm.GetNumberKets();
        for (size_t ch_ij = 0; ch_ij < nch2; ch_ij++)
        {
          TwoBodyChannel &tbc_ij = Z.modelspace->GetTwoBodyChannel(ch_ij);
          if ((tbc_lm.parity + tbc_ij.parity) % 2 != tbc_CC.parity)
            continue;
          if (std::abs(tbc_lm.Tz - tbc_ij.Tz) != tbc_CC.Tz)
            continue; // This abs bit may make things expensive...
          if ((std::abs(tbc_lm.J - tbc_ij.J) > tbc_CC.J) or (tbc_lm.J + tbc_ij.J) < tbc_CC.J)
            continue;

          size_t nkets_ij = tbc_ij.GetNumberKets();
          for (size_t iket_lm = 0; iket_lm < nkets_lm; iket_lm++)
          {
            Ket &ket_lm = tbc_lm.GetKet(iket_lm);
            size_t l = ket_lm.p;
            size_t m = ket_lm.q;
            int Jlm = tbc_lm.J;
            Orbit &ol = Z.modelspace->GetOrbit(l);
            Orbit &om = Z.modelspace->GetOrbit(m);
            double d_el = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);
            double d_em = std::abs(2 * om.n + om.l - e_fermi[om.tz2]);
            double occnat_l = ol.occ_nat;
            double occnat_m = om.occ_nat;
            if ((occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_factor_max) < Z.modelspace->GetOccNat3Cut())
              continue;
            if ((d_el + d_em) > Z.modelspace->dE3max)
              continue;

            for (size_t iket_ij = 0; iket_ij < nkets_ij; iket_ij++)
            {
              Ket &ket_ij = tbc_ij.GetKet(iket_ij);
              size_t i = ket_ij.p;
              size_t j = ket_ij.q;
              int Jij = tbc_ij.J;
              // Only compute one ordering of   lmij  vs ijlm  since we can get the other ordering from symmetry
              if ((std::min(i, j) > std::min(l, m)) or ((std::min(i, j) == std::min(l, m)) and (std::max(i, j) > std::max(l, m))) or ((std::min(i, j) == std::min(l, m)) and (std::max(i, j) == std::max(l, m)) and (Jij > Jlm)))
                continue;

              Orbit &oi = Z.modelspace->GetOrbit(i);
              Orbit &oj = Z.modelspace->GetOrbit(j);
              double occnat_i = oi.occ_nat;
              double occnat_j = oj.occ_nat;
              double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
              double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
              if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_factor_max) < Z.modelspace->GetOccNat3Cut())
                continue;
              if ((d_ei + d_ej) > Z.modelspace->dE3max)
                continue;

              size_t key = hash_key_lmij(l, m, (size_t)Jlm, i, j, (size_t)Jij);
              ph_kets3[key] = nkets3;
              //           ph_kets3[{l,m,Jlm,i,j,Jij}] = nkets3;
              nkets3++;
            }
          }

        } // for ch_ij
      }   // for ch_lm

      // so now we can allocate the ph 3-body matrices

      Z3_ph[ch_cc].zeros(2 * nkets_CC, nkets3); // Zbar_knlmij
    }                                           // for ch_cc

    Z.profiler.timer["_comm233_setup_lists"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    // Compute the particle-hole recoupled matrix elements of Z
    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ch_cc = 0; ch_cc < nch2_CC; ch_cc++)
    {
      TwoBodyChannel_CC &tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int Jph = tbc_CC.J;
      size_t nkets_CC = tbc_CC.GetNumberKets();
      size_t nkets_ph = ph_kets_cc[ch_cc].size();
      size_t nkets3 = Z3_ph[ch_cc].n_cols;

      arma::mat X2_ph(2 * nkets_CC, 2 * nkets_ph, arma::fill::zeros); // we need to store both permutations ij ji because the ph transformed 2bmes don't have as nice a symmetry
      arma::mat Y2_ph(2 * nkets_CC, 2 * nkets_ph, arma::fill::zeros);

      arma::mat X3_ph(2 * nkets_ph, nkets3, arma::fill::zeros); // Xbar_ablmij
      arma::mat Y3_ph(2 * nkets_ph, nkets3, arma::fill::zeros);

      // now we fill the matrices
      for (size_t ibra_ij = 0; ibra_ij < nkets_CC; ibra_ij++)
      {
        Ket &bra_ij = tbc_CC.GetKet(ibra_ij);
        size_t i = bra_ij.p;
        size_t j = bra_ij.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        for (size_t iket_ab = 0; iket_ab < nkets_ph; iket_ab++)
        {
          Ket &ket_ab = tbc_CC.GetKet(ph_kets_cc[ch_cc].at(iket_ab));
          size_t a = ket_ab.p;
          size_t b = ket_ab.q;
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double ja = 0.5 * oa.j2;
          double jb = 0.5 * ob.j2;
          //        double occ_factor = occfactors[ch_cc][iket_ab];

          double occ_factor = oa.occ - ob.occ;

          // XJph_ij`ab` = -sum J'  (2J'+1) { i j J_ph } XJ'_ibaj
          //                                { a b J'   }

          double xijab = 0;
          double xjiab = 0;
          double yijab = 0;
          double yjiab = 0;
          int Jprime_min = std::max(std::abs(oi.j2 - ob.j2), std::abs(oa.j2 - oj.j2)) / 2;
          int Jprime_max = std::min(oi.j2 + ob.j2, oa.j2 + oj.j2) / 2;
          for (int Jprime = Jprime_min; Jprime <= Jprime_max; Jprime++)
          {
            double sixj_ij = Z.modelspace->GetSixJ(ji, jj, Jph, ja, jb, Jprime);
            xijab -= (2 * Jprime + 1) * sixj_ij * X2.GetTBME_J(Jprime, i, b, a, j);
            yijab -= (2 * Jprime + 1) * sixj_ij * Y2.GetTBME_J(Jprime, i, b, a, j);
          }
          Jprime_min = std::max(std::abs(oj.j2 - ob.j2), std::abs(oa.j2 - oi.j2)) / 2;
          Jprime_max = std::min(oj.j2 + ob.j2, oa.j2 + oi.j2) / 2;
          for (int Jprime = Jprime_min; Jprime <= Jprime_max; Jprime++)
          {
            double sixj_ji = Z.modelspace->GetSixJ(jj, ji, Jph, ja, jb, Jprime);
            xjiab -= (2 * Jprime + 1) * sixj_ji * X2.GetTBME_J(Jprime, j, b, a, i);
            yjiab -= (2 * Jprime + 1) * sixj_ji * Y2.GetTBME_J(Jprime, j, b, a, i);
          }
          int phase_X = hX * Z.modelspace->phase((oi.j2 + oj.j2 + oa.j2 + ob.j2) / 2);
          int phase_Y = hY * Z.modelspace->phase((oi.j2 + oj.j2 + oa.j2 + ob.j2) / 2);

          X2_ph(ibra_ij, iket_ab) = xijab * occ_factor;
          Y2_ph(ibra_ij, iket_ab) = yijab * occ_factor;
          X2_ph(ibra_ij + nkets_CC, iket_ab) = xjiab * occ_factor;
          Y2_ph(ibra_ij + nkets_CC, iket_ab) = yjiab * occ_factor;
          X2_ph(ibra_ij, iket_ab + nkets_ph) = phase_X * xjiab * (-occ_factor);
          Y2_ph(ibra_ij, iket_ab + nkets_ph) = phase_Y * yjiab * (-occ_factor);
          X2_ph(ibra_ij + nkets_CC, iket_ab + nkets_ph) = phase_X * xijab * (-occ_factor);
          Y2_ph(ibra_ij + nkets_CC, iket_ab + nkets_ph) = phase_Y * yijab * (-occ_factor);

        } // for iket_ab
      }   // for ibra_ij

      // now fill the 3-body matrices
      auto &ph_kets3 = ph_kets3_lookup[ch_cc];
      for (auto iter_lmij : ph_kets3)
      {
        size_t index_lmij = iter_lmij.second;

        size_t key = iter_lmij.first;
        size_t l, m, i, j, Jlm_tmp, Jij_tmp;
        unhash_key_lmij(l, m, Jlm_tmp, i, j, Jij_tmp, key);
        int Jlm = (int)Jlm_tmp;
        int Jij = (int)Jij_tmp;

        //      size_t index_ijlm = ph_kets3.at({i,j,Jij,l,m,Jlm});
        Orbit &ol = Z.modelspace->GetOrbit(l);
        Orbit &om = Z.modelspace->GetOrbit(m);
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double occnat_l = ol.occ_nat;
        double occnat_m = om.occ_nat;
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;
        double d_el = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);
        double d_em = std::abs(2 * om.n + om.l - e_fermi[om.tz2]);
        double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);

        //      for ( size_t iket_ab=0; iket_ab<nkets_ph; iket_ab++)
        for (size_t iket_ab = 0; iket_ab < 2 * nkets_ph; iket_ab++)
        {
          Ket &ket_ab = tbc_CC.GetKet(ph_kets_cc[ch_cc].at(iket_ab % nkets_ph));

          size_t a = (iket_ab < nkets_ph) ? ket_ab.p : ket_ab.q;
          size_t b = (iket_ab < nkets_ph) ? ket_ab.q : ket_ab.p;
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double ja = 0.5 * oa.j2;
          double jb = 0.5 * ob.j2;

          double occnat_a = oa.occ_nat;
          double occnat_b = ob.occ_nat;
          double d_ea = std::abs(2 * oa.n + oa.l - e_fermi[oa.tz2]);
          double d_eb = std::abs(2 * ob.n + ob.l - e_fermi[ob.tz2]);
          double occ_factor = oa.occ - ob.occ;
          if (std::abs(occ_factor) < 1e-6)
            continue;

          if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_a * (1 - occnat_a)) < Z.modelspace->GetOccNat3Cut())
            continue;
          if ((occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_b * (1 - occnat_b)) < Z.modelspace->GetOccNat3Cut())
            continue;
          if ((d_ei + d_ej + d_ea) > Z.modelspace->dE3max)
            continue;
          if ((d_el + d_em + d_eb) > Z.modelspace->dE3max)
            continue;

          double X3bar_ablmij = 0;
          double Y3bar_ablmij = 0;

          if ((ob.l + ol.l + om.l + oa.l + oi.l + oj.l) % 2 > 0)
            continue;
          if ((ob.tz2 + ol.tz2 + om.tz2) != (oa.tz2 + oi.tz2 + oj.tz2))
            continue;

          int twoJp_min = std::max(std::abs(oa.j2 - 2 * Jij), std::abs(ob.j2 - 2 * Jlm));
          int twoJp_max = std::min((oa.j2 + 2 * Jij), (ob.j2 + 2 * Jlm));
          for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
          {
            int phase_y = Z.modelspace->phase((oa.j2 + twoJp) / 2);

            double xablmij = 0;
            double yablmij = 0;
            // get the X3 and Y3 elements in one go (to avoid doing the recoupling twice)
            if (x3_allocated and y3_allocated)
            {
              auto xandy = Y3.GetME_pn_TwoOps(Jij, Jlm, twoJp, i, j, a, l, m, b, X3, Y3);
              xablmij = xandy[0];
              yablmij = xandy[1];
            }
            else if (y3_allocated)
            {
              yablmij = Y3.GetME_pn(Jij, Jlm, twoJp, i, j, a, l, m, b);
            }
            else if (x3_allocated)
            {
              xablmij = X3.GetME_pn(Jij, Jlm, twoJp, i, j, a, l, m, b);
            }

            double sixjy = Z.modelspace->GetSixJ(Jij, Jlm, Jph, jb, ja, 0.5 * twoJp);
            X3bar_ablmij -= phase_y * (twoJp + 1) * sixjy * xablmij;
            Y3bar_ablmij -= phase_y * (twoJp + 1) * sixjy * yablmij;
          } // for twoJp

          size_t index_ab = iket_ab;
          //        size_t index_ba =  (iket_ab + nkets_ph) % (2*nkets_ph);
          X3_ph(index_ab, index_lmij) = X3bar_ablmij;
          Y3_ph(index_ab, index_lmij) = Y3bar_ablmij;

        } // for iket_ab
      }   // for iter_lmij

      Z3_ph[ch_cc] = X2_ph * Y3_ph - Y2_ph * X3_ph;

    } // for ch_cc
      //  std::cout << "DONE WITH THE CC BIT" << std::endl;

    Z.profiler.timer["_comm233_fill_matrices"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    // Now transform it back into standard coupling
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ch3 = 0; ch3 < nch3; ch3++)
    {
      auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      int twoJ = Tbc.twoJ;
      double Jtot = 0.5 * twoJ;
      for (size_t ibra = 0; ibra < nkets3; ibra++)
      {
        auto &bra = Tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        Orbit &ok = Z.modelspace->GetOrbit(k);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        double jk = 0.5 * ok.j2;
        double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
        double d_ek = std::abs(2 * ok.n + ok.l - e_fermi[ok.tz2]);
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;
        double occnat_k = ok.occ_nat;
        if ((d_ei + d_ej + d_ek) > Z.modelspace->dE3max)
          continue;
        if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
          continue;
        int J1 = bra.Jpq;

        std::vector<std::array<size_t, 3>> ijk = {{i, j, k}, {k, j, i}, {i, k, j}};

        std::vector<int> J1p_min = {J1, std::max(std::abs(ok.j2 - oj.j2), std::abs(twoJ - oi.j2)) / 2, std::max(std::abs(oi.j2 - ok.j2), std::abs(twoJ - oj.j2)) / 2};
        std::vector<int> J1p_max = {J1, std::min(ok.j2 + oj.j2, twoJ + oi.j2) / 2, std::min(oi.j2 + ok.j2, twoJ + oj.j2) / 2};
        std::vector<std::vector<double>> recouple_ijk = {{1}, {}, {}};
        for (int J1p = J1p_min[1]; J1p <= J1p_max[1]; J1p++)
          recouple_ijk[1].push_back(sqrt((2 * J1 + 1) * (2 * J1p + 1)) * Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p));

        for (int J1p = J1p_min[2]; J1p <= J1p_max[2]; J1p++)
          recouple_ijk[2].push_back(-Z.modelspace->phase((oj.j2 + ok.j2) / 2 + J1 + J1p) * sqrt((2 * J1 + 1) * (2 * J1p + 1)) * Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p));

        for (size_t iket = ibra; iket < nkets3; iket++)
        {
          auto &ket = Tbc.GetKet(iket);
          size_t l = ket.p;
          size_t m = ket.q;
          size_t n = ket.r;
          Orbit &ol = Z.modelspace->GetOrbit(l);
          Orbit &om = Z.modelspace->GetOrbit(m);
          Orbit &on = Z.modelspace->GetOrbit(n);
          double jl = 0.5 * ol.j2;
          double jm = 0.5 * om.j2;
          double jn = 0.5 * on.j2;
          double d_el = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);
          double d_em = std::abs(2 * om.n + om.l - e_fermi[om.tz2]);
          double d_en = std::abs(2 * on.n + on.l - e_fermi[on.tz2]);
          double occnat_l = ol.occ_nat;
          double occnat_m = om.occ_nat;
          double occnat_n = on.occ_nat;
          if ((d_el + d_em + d_en) > Z.modelspace->dE3max)
            continue;
          if ((occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_n * (1 - occnat_n)) < Z.modelspace->GetOccNat3Cut())
            continue;
          int J2 = ket.Jpq;

          std::vector<std::array<size_t, 3>> lmn = {{l, m, n}, {n, m, l}, {l, n, m}};
          std::vector<int> J2p_min = {J2, std::max(std::abs(on.j2 - om.j2), std::abs(twoJ - ol.j2)) / 2, std::max(std::abs(ol.j2 - on.j2), std::abs(twoJ - om.j2)) / 2};
          std::vector<int> J2p_max = {J2, std::min(on.j2 + om.j2, twoJ + ol.j2) / 2, std::min(ol.j2 + on.j2, twoJ + om.j2) / 2};
          std::vector<std::vector<double>> recouple_lmn = {{1}, {}, {}};

          for (int J2p = J2p_min[1]; J2p <= J2p_max[1]; J2p++)
            recouple_lmn[1].push_back(sqrt((2 * J2 + 1) * (2 * J2p + 1)) * Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p));

          for (int J2p = J2p_min[2]; J2p <= J2p_max[2]; J2p++)
            recouple_lmn[2].push_back(-Z.modelspace->phase((om.j2 + on.j2) / 2 + J2 + J2p) * sqrt((2 * J2 + 1) * (2 * J2p + 1)) * Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p));

          double z_ijklmn = 0;

          for (int perm_ijk = 0; perm_ijk < 3; perm_ijk++)
          {
            size_t I1 = ijk[perm_ijk][0];
            size_t I2 = ijk[perm_ijk][1];
            size_t I3 = ijk[perm_ijk][2];
            Orbit &o1 = Z.modelspace->GetOrbit(I1);
            Orbit &o2 = Z.modelspace->GetOrbit(I2);
            Orbit &o3 = Z.modelspace->GetOrbit(I3);
            //          double occnat_1 = o1.occ_nat;
            //          double occnat_2 = o2.occ_nat;
            //          double d_e1 = std::abs( 2*o1.n + o1.l - e_fermi[o1.tz2]);
            //          double d_e2 = std::abs( 2*o2.n + o2.l - e_fermi[o2.tz2]);
            double j3 = 0.5 * o3.j2;
            for (int J1p = J1p_min[perm_ijk]; J1p <= J1p_max[perm_ijk]; J1p++)
            {
              if ((I1 == I2) and J1p % 2 != 0)
                continue;
              double rec_ijk = recouple_ijk[perm_ijk].at(J1p - J1p_min[perm_ijk]);
              for (int perm_lmn = 0; perm_lmn < 3; perm_lmn++)
              {
                size_t I4 = lmn[perm_lmn][0];
                size_t I5 = lmn[perm_lmn][1];
                size_t I6 = lmn[perm_lmn][2];
                Orbit &o4 = Z.modelspace->GetOrbit(I4);
                Orbit &o5 = Z.modelspace->GetOrbit(I5);
                Orbit &o6 = Z.modelspace->GetOrbit(I6);
                //              double occnat_4 = o4.occ_nat;
                //              double occnat_5 = o5.occ_nat;
                //              double d_e4 = std::abs( 2*o4.n + o4.l - e_fermi[o4.tz2]);
                //              double d_e5 = std::abs( 2*o5.n + o5.l - e_fermi[o5.tz2]);
                double j6 = 0.5 * o6.j2;
                for (int J2p = J2p_min[perm_lmn]; J2p <= J2p_max[perm_lmn]; J2p++)
                {
                  if ((I4 == I5) and J2p % 2 != 0)
                    continue;
                  double rec_lmn = recouple_lmn[perm_lmn].at(J2p - J2p_min[perm_lmn]);

                  double z_123456 = 0;
                  int Jph_min = std::max(std::abs(2 * J1p - 2 * J2p), std::abs(o3.j2 - o6.j2)) / 2;
                  int Jph_max = std::min(2 * J1p + 2 * J2p, o3.j2 + o6.j2) / 2;
                  int Tz_ph = std::abs(o3.tz2 - o6.tz2) / 2;
                  int parity_ph = (o3.l + o6.l) % 2;

                  for (int Jph = Jph_min; Jph <= Jph_max; Jph++)
                  {
                    size_t ch_ph = Z.modelspace->GetTwoBodyChannelIndex(Jph, parity_ph, Tz_ph);
                    TwoBodyChannel_CC &tbc_CC = Z.modelspace->GetTwoBodyChannel_CC(ch_ph);
                    //                   size_t nkets_CC = tbc_CC.GetNumberKets();

                    // we only compute one ordering of the right side of Zbar, but we can get the other ordering by symmetry
                    // so we check here which ordering we need, and pick up the corresponding phase factor if needed.
                    size_t index_36 = tbc_CC.GetLocalIndex(I3, I6);
                    size_t key_4512 = hash_key_lmij(std::min(I4, I5), std::max(I4, I5), (size_t)J2p, std::min(I1, I2), std::max(I1, I2), (size_t)J1p);
                    int phase_45 = (I4 > I5) ? -Z.modelspace->phase((o4.j2 + o5.j2) / 2 - J2p) : 1;
                    int phase_12 = (I1 > I2) ? -Z.modelspace->phase((o1.j2 + o2.j2) / 2 - J1p) : 1;
                    int flip_phase = 1;

                    // if any of these are true, we're looking for the not-stored ordering
                    if ((std::min(I1, I2) > std::min(I4, I5)) or ((std::min(I1, I2) == std::min(I4, I5)) and (std::max(I1, I2) > std::max(I4, I5))) or ((std::min(I1, I2) == std::min(I4, I5)) and (std::max(I1, I2) == std::max(I4, I5)) and (J1p > J2p)))
                    {
                      index_36 = tbc_CC.GetLocalIndex(I6, I3);
                      key_4512 = hash_key_lmij(std::min(I1, I2), std::max(I1, I2), (size_t)J1p, std::min(I4, I5), std::max(I4, I5), (size_t)J2p);
                      flip_phase = Z.modelspace->phase((o3.j2 - o6.j2) / 2);
                    }
                    size_t index_4512 = ph_kets3_lookup[ch_ph].at(key_4512);

                    // < 36` Jph | Z | 45 J2p  (12 J1p)` ; Jph >
                    double z_364512 = flip_phase * phase_45 * phase_12 * Z3_ph[ch_ph](index_36, index_4512);

                    // Swapped the ordering here to match with what's done in PreCalculateSixJ()
                    //                   double sixj_ph = Z.modelspace->GetSixJ(J1p,j3,Jtot,  j6,J2p,Jph);
                    double sixj_ph = Z.modelspace->GetSixJ(J1p, J2p, Jph, j6, j3, Jtot);
                    int phase_ph = -Z.modelspace->phase((twoJ + o3.j2) / 2);
                    z_123456 -= phase_ph * (2 * Jph + 1) * sixj_ph * z_364512;
                  } // for Jph

                  // include the recoupling coefficients needed for the permutations
                  z_ijklmn += rec_ijk * rec_lmn * z_123456;

                } // for J2p
              }   // for perm_lmn
            }     // for J1p
          }       // for perm_ijk

          Z3.AddToME_pn_ch(ch3, ch3, ibra, iket, z_ijklmn);
        } // for iket
      }   // for ibra
    }     // for ch3
    Z.profiler.timer["_comm233_pph_recouple"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  void comm233_phss_debug(const Operator &X, const Operator &Y, Operator &Z)
  {

    double tstart = omp_get_wtime();
    if (Commutator::verbose)
      std::cout << __func__ << std::endl;
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch3 = 0; ch3 < nch3; ch3++)
    {
      auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      int twoJ = Tbc.twoJ;
      double Jtot = 0.5 * twoJ;
      for (size_t ibra = 0; ibra < nkets3; ibra++)
      {
        auto &bra = Tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        Orbit &ok = Z.modelspace->GetOrbit(k);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        double jk = 0.5 * ok.j2;
        double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
        double d_ek = std::abs(2 * ok.n + ok.l - e_fermi[ok.tz2]);
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;
        double occnat_k = ok.occ_nat;
        if ((d_ei + d_ej + d_ek) > Z.modelspace->dE3max)
          continue;
        if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
          continue;
        int J1 = bra.Jpq;
        for (size_t iket = ibra; iket < nkets3; iket++)
        {
          auto &ket = Tbc.GetKet(iket);
          size_t l = ket.p;
          size_t m = ket.q;
          size_t n = ket.r;
          Orbit &ol = Z.modelspace->GetOrbit(l);
          Orbit &om = Z.modelspace->GetOrbit(m);
          Orbit &on = Z.modelspace->GetOrbit(n);
          double jl = 0.5 * ol.j2;
          double jm = 0.5 * om.j2;
          double jn = 0.5 * on.j2;
          double d_el = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);
          double d_em = std::abs(2 * om.n + om.l - e_fermi[om.tz2]);
          double d_en = std::abs(2 * on.n + on.l - e_fermi[on.tz2]);
          double occnat_l = ol.occ_nat;
          double occnat_m = om.occ_nat;
          double occnat_n = on.occ_nat;
          if ((d_el + d_em + d_en) > Z.modelspace->dE3max)
            continue;
          if ((occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_n * (1 - occnat_n)) < Z.modelspace->GetOccNat3Cut())
            continue;
          int J2 = ket.Jpq;

          double z_ijklmn = 0;

          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double ja = 0.5 * oa.j2;
            double occnat_a = oa.occ_nat;
            double d_ea = std::abs(2 * oa.n + oa.l - e_fermi[oa.tz2]);
            bool keep_ija = ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_a * (1 - occnat_a)) >= Z.modelspace->GetOccNat3Cut());
            bool keep_kja = ((occnat_k * (1 - occnat_k) * occnat_j * (1 - occnat_j) * occnat_a * (1 - occnat_a)) >= Z.modelspace->GetOccNat3Cut());
            bool keep_ika = ((occnat_i * (1 - occnat_i) * occnat_k * (1 - occnat_k) * occnat_a * (1 - occnat_a)) >= Z.modelspace->GetOccNat3Cut());
            keep_ija = keep_ija and (d_ei + d_ej + d_ea <= Z.modelspace->dE3max);
            keep_kja = keep_kja and (d_ek + d_ej + d_ea <= Z.modelspace->dE3max);
            keep_ika = keep_ika and (d_ei + d_ek + d_ea <= Z.modelspace->dE3max);
            if (not(keep_ija or keep_kja or keep_ika))
              continue;

            for (auto &b : Z.modelspace->all_orbits)
            {
              //           if (b<a) continue;
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double jb = 0.5 * ob.j2;
              double occnat_b = ob.occ_nat;
              double d_eb = std::abs(2 * ob.n + ob.l - e_fermi[ob.tz2]);
              double occ_factor = oa.occ - ob.occ;
              if (std::abs(occ_factor) < 1e-6)
                continue;

              bool keep_lmb = ((occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_b * (1 - occnat_b)) >= Z.modelspace->GetOccNat3Cut());
              bool keep_lnb = ((occnat_l * (1 - occnat_l) * occnat_n * (1 - occnat_n) * occnat_b * (1 - occnat_b)) >= Z.modelspace->GetOccNat3Cut());
              bool keep_nmb = ((occnat_n * (1 - occnat_n) * occnat_m * (1 - occnat_m) * occnat_b * (1 - occnat_b)) >= Z.modelspace->GetOccNat3Cut());
              keep_lmb = keep_lmb and (d_el + d_em + d_eb <= Z.modelspace->dE3max);
              keep_lnb = keep_lnb and (d_el + d_en + d_eb <= Z.modelspace->dE3max);
              keep_nmb = keep_nmb and (d_en + d_em + d_eb <= Z.modelspace->dE3max);
              if (not(keep_lmb or keep_lnb or keep_nmb))
                continue;

              // Direct term i.e. Z1
              if ((ob.l + ok.l + oa.l + on.l) % 2 == 0 and ((ob.tz2 + ok.tz2) == (oa.tz2 + on.tz2)) and keep_ija and keep_lmb)
              {
                int Jx_min = std::max(std::abs(ob.j2 - ok.j2), std::abs(oa.j2 - on.j2)) / 2;
                int Jx_max = std::min(ob.j2 + ok.j2, oa.j2 + on.j2) / 2;
                int twoJJ_min = std::max(std::abs(oa.j2 - 2 * J1), std::abs(ob.j2 - 2 * J2));
                int twoJJ_max = std::min(oa.j2 + 2 * J1, ob.j2 + 2 * J2);
                int phase = -Z.modelspace->phase((ok.j2 + on.j2) / 2 + J1 + J2);
                for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                {
                  double xijalmb = X3.GetME_pn(J1, J2, twoJJ, i, j, a, l, m, b);
                  double yijalmb = Y3.GetME_pn(J1, J2, twoJJ, i, j, a, l, m, b);
                  double JJtot = 0.5 * twoJJ;
                  for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                  {
                    double xbkan = X2.GetTBME_J(Jx, b, k, a, n);
                    double ybkan = Y2.GetTBME_J(Jx, b, k, a, n);
                    double hats = (twoJJ + 1) * (2 * Jx + 1);
                    double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2, jk, J1, Jtot, Jx, ja, jn);
                    z_ijklmn -= occ_factor * hats * phase * ninej * (xbkan * yijalmb - ybkan * xijalmb);
                  } // for Jx
                }   // for twoJ
              }     // Z1 block

              // Z2  ~ Pik Z1     Xbian Ykjalmb
              if ((ob.l + oi.l + oa.l + on.l) % 2 == 0 and ((ob.tz2 + oi.tz2) == (oa.tz2 + on.tz2)) and keep_kja and keep_lmb)
              {
                int J1p_min = std::abs(ok.j2 - oj.j2) / 2;
                int J1p_max = (ok.j2 + oj.j2) / 2;
                int Jx_min = std::max(std::abs(ob.j2 - oi.j2), std::abs(oa.j2 - on.j2)) / 2;
                int Jx_max = std::min(ob.j2 + oi.j2, oa.j2 + on.j2) / 2;
                for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                {
                  double sixj = Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p);
                  if (std::abs(sixj) < 1e-6)
                    continue;
                  int twoJJ_min = std::max(std::abs(ob.j2 - 2 * J2), std::abs(oa.j2 - 2 * J1p));
                  int twoJJ_max = std::min(ob.j2 + 2 * J2, oa.j2 + 2 * J1p);
                  int phase = Z.modelspace->phase((oi.j2 + on.j2) / 2 + J1p + J2);
                  for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                  {
                    double xkjalmb = X3.GetME_pn(J1p, J2, twoJJ, k, j, a, l, m, b);
                    double ykjalmb = Y3.GetME_pn(J1p, J2, twoJJ, k, j, a, l, m, b);
                    double JJtot = 0.5 * twoJJ;
                    for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                    {
                      double xbian = X2.GetTBME_J(Jx, b, i, a, n);
                      double ybian = Y2.GetTBME_J(Jx, b, i, a, n);
                      double hats = (twoJJ + 1) * (2 * Jx + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1));
                      double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2, ji, J1p, Jtot, Jx, ja, jn);
                      z_ijklmn += occ_factor * hats * phase * sixj * ninej * (xbian * ykjalmb - ybian * xkjalmb);
                    } // for Jx
                  }   // for twoJ
                }     // for J1p
              }       // Z2 block

              // Z3  ~ Pjk Z1      Xbjan Yikalmb
              if ((ob.l + oj.l + oa.l + on.l) % 2 == 0 and ((ob.tz2 + oj.tz2) == (oa.tz2 + on.tz2)) and keep_ika and keep_lmb)
              {
                int J1p_min = std::abs(oi.j2 - ok.j2) / 2;
                int J1p_max = (oi.j2 + ok.j2) / 2;
                int Jx_min = std::max(std::abs(ob.j2 - oj.j2), std::abs(oa.j2 - on.j2)) / 2;
                int Jx_max = std::min(ob.j2 + oj.j2, oa.j2 + on.j2) / 2;
                for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                {
                  double sixj = Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p);
                  if (std::abs(sixj) < 1e-6)
                    continue;
                  int twoJJ_min = std::max(std::abs(ob.j2 - 2 * J2), std::abs(oa.j2 - 2 * J1p));
                  int twoJJ_max = std::min(ob.j2 + 2 * J2, oa.j2 + 2 * J1p);
                  int phase = Z.modelspace->phase((ok.j2 + on.j2) / 2 + J1 + J2);
                  for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                  {
                    double xikalmb = X3.GetME_pn(J1p, J2, twoJJ, i, k, a, l, m, b);
                    double yikalmb = Y3.GetME_pn(J1p, J2, twoJJ, i, k, a, l, m, b);
                    double JJtot = 0.5 * twoJJ;
                    for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                    {
                      double xbjan = X2.GetTBME_J(Jx, b, j, a, n);
                      double ybjan = Y2.GetTBME_J(Jx, b, j, a, n);
                      double hats = (twoJJ + 1) * (2 * Jx + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1));
                      double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2, jj, J1p, Jtot, Jx, ja, jn);
                      z_ijklmn += occ_factor * hats * phase * sixj * ninej * (xbjan * yikalmb - ybjan * xikalmb);
                    } // for Jx
                  }   // for twoJ
                }     // for J1p
              }       // Z3 block

              // Z4  ~ Pln Z1    Xbkal Y ijanmb
              if ((ob.l + ok.l + oa.l + ol.l) % 2 == 0 and ((ob.tz2 + ok.tz2) == (oa.tz2 + ol.tz2)) and keep_ija and keep_nmb)
              {
                int J2p_min = std::abs(on.j2 - om.j2) / 2;
                int J2p_max = (on.j2 + om.j2) / 2;
                int Jx_min = std::max(std::abs(ob.j2 - ok.j2), std::abs(oa.j2 - ol.j2)) / 2;
                int Jx_max = std::min(ob.j2 + ok.j2, oa.j2 + ol.j2) / 2;
                for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                {
                  double sixj = Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p);
                  if (std::abs(sixj) < 1e-6)
                    continue;
                  int twoJJ_min = std::max(std::abs(ob.j2 - 2 * J2p), std::abs(oa.j2 - 2 * J1));
                  int twoJJ_max = std::min(ob.j2 + 2 * J2p, oa.j2 + 2 * J1);
                  int phase = Z.modelspace->phase((ok.j2 + ol.j2) / 2 + J1 + J2p);
                  for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                  {
                    double xijanmb = X3.GetME_pn(J1, J2p, twoJJ, i, j, a, n, m, b);
                    double yijanmb = Y3.GetME_pn(J1, J2p, twoJJ, i, j, a, n, m, b);
                    double JJtot = 0.5 * twoJJ;
                    for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                    {
                      double xbkal = X2.GetTBME_J(Jx, b, k, a, l);
                      double ybkal = Y2.GetTBME_J(Jx, b, k, a, l);
                      double hats = (twoJJ + 1) * (2 * Jx + 1) * sqrt((2 * J2 + 1) * (2 * J2p + 1));
                      double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2p, jk, J1, Jtot, Jx, ja, jl);
                      z_ijklmn += occ_factor * hats * phase * sixj * ninej * (xbkal * yijanmb - ybkal * xijanmb);
                    } // for Jx
                  }   // for twoJ
                }     // for J1p
              }       // Z4 block

              // Z5  ~ Pmn Z1    Xbkam Yijalnb
              if ((ob.l + ok.l + oa.l + om.l) % 2 == 0 and ((ob.tz2 + ok.tz2) == (oa.tz2 + om.tz2)) and keep_ija and keep_lnb)
              {
                int J2p_min = std::abs(ol.j2 - on.j2) / 2;
                int J2p_max = (ol.j2 + on.j2) / 2;
                int Jx_min = std::max(std::abs(ob.j2 - ok.j2), std::abs(oa.j2 - om.j2)) / 2;
                int Jx_max = std::min(ob.j2 + ok.j2, oa.j2 + om.j2) / 2;
                for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                {
                  double sixj = Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p);
                  if (std::abs(sixj) < 1e-6)
                    continue;
                  int twoJJ_min = std::max(std::abs(ob.j2 - 2 * J2p), std::abs(oa.j2 - 2 * J1));
                  int twoJJ_max = std::min(ob.j2 + 2 * J2p, oa.j2 + 2 * J1);
                  int phase = Z.modelspace->phase((ok.j2 + on.j2) / 2 + J1 + J2);
                  for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                  {
                    double xijalnb = X3.GetME_pn(J1, J2p, twoJJ, i, j, a, l, n, b);
                    double yijalnb = Y3.GetME_pn(J1, J2p, twoJJ, i, j, a, l, n, b);
                    double JJtot = 0.5 * twoJJ;
                    for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                    {
                      double xbkam = X2.GetTBME_J(Jx, b, k, a, m);
                      double ybkam = Y2.GetTBME_J(Jx, b, k, a, m);
                      double hats = (twoJJ + 1) * (2 * Jx + 1) * sqrt((2 * J2 + 1) * (2 * J2p + 1));
                      double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2p, jk, J1, Jtot, Jx, ja, jm);
                      z_ijklmn += occ_factor * hats * phase * sixj * ninej * (xbkam * yijalnb - ybkam * xijalnb);
                    } // for Jx
                  }   // for twoJ
                }     // for J1p
              }       // Z5 block

              // Z6  ~ Pik Pln Z1     Xbial Ykjanmb
              if ((ob.l + oi.l + oa.l + ol.l) % 2 == 0 and ((ob.tz2 + oi.tz2) == (oa.tz2 + ol.tz2)) and keep_kja and keep_nmb)
              {
                int J1p_min = std::abs(ok.j2 - oj.j2) / 2;
                int J1p_max = (ok.j2 + oj.j2) / 2;
                int J2p_min = std::abs(on.j2 - om.j2) / 2;
                int J2p_max = (on.j2 + om.j2) / 2;
                int Jx_min = std::max(std::abs(ob.j2 - oi.j2), std::abs(oa.j2 - ol.j2)) / 2;
                int Jx_max = std::min(ob.j2 + oi.j2, oa.j2 + ol.j2) / 2;
                for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                {
                  double sixj1 = Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p);
                  if (std::abs(sixj1) < 1e-6)
                    continue;
                  for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                  {
                    double sixj2 = Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p);
                    if (std::abs(sixj2) < 1e-6)
                      continue;
                    int twoJJ_min = std::max(std::abs(ob.j2 - 2 * J2p), std::abs(oa.j2 - 2 * J1p));
                    int twoJJ_max = std::min(ob.j2 + 2 * J2p, oa.j2 + 2 * J1p);
                    int phase = -Z.modelspace->phase((oi.j2 + ol.j2) / 2 + J1p + J2p);
                    for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                    {
                      double xkjanmb = X3.GetME_pn(J1p, J2p, twoJJ, k, j, a, n, m, b);
                      double ykjanmb = Y3.GetME_pn(J1p, J2p, twoJJ, k, j, a, n, m, b);
                      double JJtot = 0.5 * twoJJ;
                      for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                      {
                        double xbial = X2.GetTBME_J(Jx, b, i, a, l);
                        double ybial = Y2.GetTBME_J(Jx, b, i, a, l);
                        double hats = (twoJJ + 1) * (2 * Jx + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1) * (2 * J2 + 1) * (2 * J2p + 1));
                        double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2p, ji, J1p, Jtot, Jx, ja, jl);
                        z_ijklmn -= occ_factor * hats * phase * sixj1 * sixj2 * ninej * (xbial * ykjanmb - ybial * xkjanmb);
                      } // for Jx
                    }   // for twoJ
                  }     // for J2p
                }       // for J1p
              }         // Z6 block

              // Z7  ~ Pjk Pln Z1     Xbjal Yikanmb
              if ((ob.l + oj.l + oa.l + ol.l) % 2 == 0 and ((ob.tz2 + oj.tz2) == (oa.tz2 + ol.tz2)) and keep_ika and keep_nmb)
              {
                int J1p_min = std::abs(oi.j2 - ok.j2) / 2;
                int J1p_max = (oi.j2 + ok.j2) / 2;
                int J2p_min = std::abs(on.j2 - om.j2) / 2;
                int J2p_max = (on.j2 + om.j2) / 2;
                int Jx_min = std::max(std::abs(ob.j2 - oj.j2), std::abs(oa.j2 - ol.j2)) / 2;
                int Jx_max = std::min(ob.j2 + oj.j2, oa.j2 + ol.j2) / 2;
                for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                {
                  double sixj1 = Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p);
                  if (std::abs(sixj1) < 1e-6)
                    continue;
                  for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                  {
                    double sixj2 = Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p);
                    if (std::abs(sixj2) < 1e-6)
                      continue;
                    int twoJJ_min = std::max(std::abs(ob.j2 - 2 * J2p), std::abs(oa.j2 - 2 * J1p));
                    int twoJJ_max = std::min(ob.j2 + 2 * J2p, oa.j2 + 2 * J1p);
                    int phase = -Z.modelspace->phase((ok.j2 + ol.j2) / 2 + J1 + J2p);
                    for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                    {
                      double xikanmb = X3.GetME_pn(J1p, J2p, twoJJ, i, k, a, n, m, b);
                      double yikanmb = Y3.GetME_pn(J1p, J2p, twoJJ, i, k, a, n, m, b);
                      double JJtot = 0.5 * twoJJ;
                      for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                      {
                        double xbjal = X2.GetTBME_J(Jx, b, j, a, l);
                        double ybjal = Y2.GetTBME_J(Jx, b, j, a, l);
                        double hats = (twoJJ + 1) * (2 * Jx + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1) * (2 * J2 + 1) * (2 * J2p + 1));
                        double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2p, jj, J1p, Jtot, Jx, ja, jl);
                        z_ijklmn -= occ_factor * hats * phase * sixj1 * sixj2 * ninej * (xbjal * yikanmb - ybjal * xikanmb);
                      } // for Jx
                    }   // for twoJ
                  }     // for J2p
                }       // for J1p
              }         // Z7 block

              // Z8  ~ Pik mln Z1     Xbiam Ykjalnb
              if ((ob.l + oi.l + oa.l + om.l) % 2 == 0 and ((ob.tz2 + oi.tz2) == (oa.tz2 + om.tz2)) and keep_kja and keep_lnb)
              {
                int J1p_min = std::abs(ok.j2 - oj.j2) / 2;
                int J1p_max = (ok.j2 + oj.j2) / 2;
                int J2p_min = std::abs(ol.j2 - on.j2) / 2;
                int J2p_max = (ol.j2 + on.j2) / 2;
                int Jx_min = std::max(std::abs(ob.j2 - oi.j2), std::abs(oa.j2 - om.j2)) / 2;
                int Jx_max = std::min(ob.j2 + oi.j2, oa.j2 + om.j2) / 2;
                for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                {
                  double sixj1 = Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p);
                  if (std::abs(sixj1) < 1e-6)
                    continue;
                  for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                  {
                    double sixj2 = Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p);
                    if (std::abs(sixj2) < 1e-6)
                      continue;
                    int twoJJ_min = std::max(std::abs(ob.j2 - 2 * J2p), std::abs(oa.j2 - 2 * J1p));
                    int twoJJ_max = std::min(ob.j2 + 2 * J2p, oa.j2 + 2 * J1p);
                    int phase = -Z.modelspace->phase((oi.j2 + on.j2) / 2 + J1p + J2);
                    for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                    {
                      double xkjalnb = X3.GetME_pn(J1p, J2p, twoJJ, k, j, a, l, n, b);
                      double ykjalnb = Y3.GetME_pn(J1p, J2p, twoJJ, k, j, a, l, n, b);
                      double JJtot = 0.5 * twoJJ;
                      for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                      {
                        double xbiam = X2.GetTBME_J(Jx, b, i, a, m);
                        double ybiam = Y2.GetTBME_J(Jx, b, i, a, m);
                        double hats = (twoJJ + 1) * (2 * Jx + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1) * (2 * J2 + 1) * (2 * J2p + 1));
                        double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2p, ji, J1p, Jtot, Jx, ja, jm);
                        z_ijklmn -= occ_factor * hats * phase * sixj1 * sixj2 * ninej * (xbiam * ykjalnb - ybiam * xkjalnb);
                      } // for Jx
                    }   // for twoJ
                  }     // for J2p
                }       // for J1p
              }         // Z8 block

              // Z9  ~ Pjk mln Z1     Xbjam Yikalnb
              if ((ob.l + oj.l + oa.l + om.l) % 2 == 0 and ((ob.tz2 + oj.tz2) == (oa.tz2 + om.tz2)) and keep_ika and keep_lnb)
              {
                int J1p_min = std::abs(oi.j2 - ok.j2) / 2;
                int J1p_max = (oi.j2 + ok.j2) / 2;
                int J2p_min = std::abs(ol.j2 - on.j2) / 2;
                int J2p_max = (ol.j2 + on.j2) / 2;
                int Jx_min = std::max(std::abs(ob.j2 - oj.j2), std::abs(oa.j2 - om.j2)) / 2;
                int Jx_max = std::min(ob.j2 + oj.j2, oa.j2 + om.j2) / 2;
                for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                {
                  double sixj1 = Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p);
                  if (std::abs(sixj1) < 1e-6)
                    continue;
                  for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                  {
                    double sixj2 = Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p);
                    if (std::abs(sixj2) < 1e-6)
                      continue;
                    int twoJJ_min = std::max(std::abs(ob.j2 - 2 * J2p), std::abs(oa.j2 - 2 * J1p));
                    int twoJJ_max = std::min(ob.j2 + 2 * J2p, oa.j2 + 2 * J1p);
                    int phase = -Z.modelspace->phase((ok.j2 + on.j2) / 2 + J1 + J2);
                    for (int twoJJ = twoJJ_min; twoJJ <= twoJJ_max; twoJJ += 2)
                    {
                      double xikalnb = X3.GetME_pn(J1p, J2p, twoJJ, i, k, a, l, n, b);
                      double yikalnb = Y3.GetME_pn(J1p, J2p, twoJJ, i, k, a, l, n, b);
                      double JJtot = 0.5 * twoJJ;
                      for (int Jx = Jx_min; Jx <= Jx_max; Jx++)
                      {
                        double xbjam = X2.GetTBME_J(Jx, b, j, a, m);
                        double ybjam = Y2.GetTBME_J(Jx, b, j, a, m);
                        double hats = (twoJJ + 1) * (2 * Jx + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1) * (2 * J2 + 1) * (2 * J2p + 1));
                        double ninej = Z.modelspace->GetNineJ(jb, JJtot, J2p, jj, J1p, Jtot, Jx, ja, jm);
                        z_ijklmn -= occ_factor * hats * phase * sixj1 * sixj2 * ninej * (xbjam * yikalnb - ybjam * xikalnb);
                      } // for Jx
                    }   // for twoJ
                  }     // for J2p
                }       // for J1p
              }         // Z9 block

            } // for b
          }   // for a
          Z3.AddToME_pn_ch(ch3, ch3, ibra, iket, z_ijklmn);
        } // for iket
      }   // for ibra
    }     // for ch3

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  //*****************************************************************************************
  //
  //  |    |    |     Uncoupled expression:
  // i|   j|   k|       Z_ijklmn = 1/6 sum_{abc} (nanbnc + n`an`bn`c) ( X_{ijkabc} Y_{abclmn} - Y_{ijkabc} X_{abclmn}
  //  *~~~[X]~~~*
  // a|   b|   c|
  //  *~~~[Y]~~~*     Coupled expression:
  // l|   m|   n|     Z_{ijklmn}^{J1,J2,J} = 1/6 sum_{abc} sum_{J'} (nanbnc + n`an`bn`c) ( X_{ijkabc}^{J1,J',J} Y_{abclmn}^{J',J2,J} - X<->Y )
  //  |    |    |
  //
  //  Checked with UnitTest and passed
  //
  void comm333_ppp_hhhss(const Operator &X, const Operator &Y, Operator &Z)
  {

    double tstart = omp_get_wtime();
    if (Commutator::verbose)
      std::cout << __func__ << std::endl;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    int hX = X.IsHermitian() ? 1 : -1;
    int hY = Y.IsHermitian() ? 1 : -1;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    bool x3_allocated = X3.IsAllocated();
    bool y3_allocated = Y3.IsAllocated();
    if (not(x3_allocated and y3_allocated))
      return;

    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ch3 = 0; ch3 < nch3; ch3++)
    {
      auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      std::vector<size_t> abc_kets;
      std::vector<double> abc_factors;

      for (size_t iket_abc = 0; iket_abc < nkets3; iket_abc++)
      {
        auto &ket_abc = Tbc.GetKet(iket_abc);
        double na = ket_abc.op->occ;
        double nb = ket_abc.oq->occ;
        double nc = ket_abc.oR->occ;
        //       double occ_factor = na*nb*nc - (1-na)*(1-nb)*(1-nc); // typo??? Should be plus, not minus
        double occ_factor = na * nb * nc + (1 - na) * (1 - nb) * (1 - nc); // fixed by SRS July 19 2022
        if (std::abs(occ_factor) < 1e-6)
          continue;
        double d_ea = std::abs(2 * ket_abc.op->n + ket_abc.op->l - e_fermi[ket_abc.op->tz2]);
        double d_eb = std::abs(2 * ket_abc.oq->n + ket_abc.oq->l - e_fermi[ket_abc.oq->tz2]);
        double d_ec = std::abs(2 * ket_abc.oR->n + ket_abc.oR->l - e_fermi[ket_abc.oR->tz2]);
        double occnat_a = ket_abc.op->occ_nat;
        double occnat_b = ket_abc.oq->occ_nat;
        double occnat_c = ket_abc.oR->occ_nat;
        if ((d_ea + d_eb + d_ec) > Z.modelspace->GetdE3max())
          continue;
        if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
          continue;
        double symm_factor = 1;
        if ((ket_abc.p == ket_abc.q) and (ket_abc.p == ket_abc.r))
          symm_factor = 1. / 6;
        else if ((ket_abc.p == ket_abc.q) or (ket_abc.q == ket_abc.r))
          symm_factor = 3. / 6;
        abc_kets.push_back(iket_abc);
        abc_factors.push_back(occ_factor * symm_factor);
      } // for iket_abc

      std::vector<size_t> bras_kept;
      for (size_t ibra = 0; ibra < nkets3; ibra++)
      {
        Ket3 &bra = Tbc.GetKet(ibra);
        double d_ei = std::abs(2 * bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
        double d_ej = std::abs(2 * bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
        double d_ek = std::abs(2 * bra.oR->n + bra.oR->l - e_fermi[bra.oR->tz2]);
        double occnat_i = bra.op->occ_nat;
        double occnat_j = bra.oq->occ_nat;
        double occnat_k = bra.oR->occ_nat;
        if ((d_ei + d_ej + d_ek) > Z.modelspace->GetdE3max())
          continue;
        if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
          continue;
        bras_kept.push_back(ibra);
      }

      size_t dim_abc = abc_kets.size();
      arma::mat XMAT(nkets3, dim_abc, arma::fill::zeros);
      arma::mat YMAT(nkets3, dim_abc, arma::fill::zeros);
      for (size_t index_abc = 0; index_abc < dim_abc; index_abc++)
      {
        size_t iket_abc = abc_kets[index_abc];
        double factor = abc_factors[index_abc];
        for (size_t ibra : bras_kept)
        {
          double xijkabc = X3.GetME_pn_ch(ch3, ch3, ibra, iket_abc);
          double yijkabc = Y3.GetME_pn_ch(ch3, ch3, ibra, iket_abc);
          XMAT(ibra, index_abc) = factor * xijkabc;
          YMAT(ibra, index_abc) = yijkabc;
        } // for ibra
      }   // for index_abc

      // Now do the mat mult
      arma::mat ZMAT = hY * XMAT * YMAT.t();
      ZMAT -= hX * hY * ZMAT.t();

      // Store the results
      for (size_t ibra : bras_kept)
      {
        for (size_t iket : bras_kept)
        {
          if (iket < ibra)
            continue;
          Z3.AddToME_pn_ch(ch3, ch3, ibra, iket, ZMAT(ibra, iket));

        } // for iket
      }   // for ibra
    }     // for ch3

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  //*****************************************************************************************
  //
  //  |i  |j      k/     Uncoupled expression:
  //  |   |       /      Z_ijklmn =-1/2 P(ij/k)P(lm/n) sum_{abc} (nanbn`c + na`n`bnc) ( X_{ijcabn} Y_{abklmc} - Y_{ijcabn} X_{abklmc} )
  //  *~~[X]~*   /
  //  |   |  |\ /        Coupled expression:
  // a|  b| c| /         Z_{ijklmn}^{J1,J2,J} = +1/2 P^{J1,J}(ij/k)P^{J2,j}_{lm/n) sum_{abc} (nanbn`c+n`an`bnc) sum_{J',J",J3} (2J'+1)(2J"+1)
  //  |   |  |/ \                                  { k   J1  J  }
  //  *~~[Y]~*   \                               * { J3  J'  n  } * ( X_{ijcabn}^{J1,J3,J'} Y_{abklmc}^{J3,J2,J"} - X<->Y )
  //  |   |       \                                { J"  c   J2 }
  //  |l  |m       \ n
  //
  //
  //  Checked with UnitTest and passed
  // TODO: Go back and implement the dE3max cut on the internal terms as well
  //
  void comm333_pph_hhpss(const Operator &X, const Operator &Y, Operator &Z)
  {

    double tstart = omp_get_wtime();
    if (Commutator::verbose)
      std::cout << __func__ << std::endl;
    double t_internal = omp_get_wtime();

    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    int hX = X.IsHermitian() ? 1 : -1;
    int hY = Y.IsHermitian() ? 1 : -1;

    bool x3_allocated = X3.IsAllocated();
    bool y3_allocated = Y3.IsAllocated();
    if (not(x3_allocated and y3_allocated))
      return;

    size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();

    std::map<std::array<int, 3>, size_t> channels_pph; // map twoJ_ph, parity_ph, twoTz_ph  -> channel index
    std::vector<std::array<int, 3>> channel_list_pph;  // vector  channel_index -> twoJ-ph, parity_ph, twoTz_ph.

    //  std::vector< std::map<std::array<int,4>, size_t>> ket_lookup_pph; // in a given channel, map  i,j,k,Jij -> matrix index

    auto hash_key_ijnJ = [](size_t i, size_t j, size_t n, size_t Jij)
    { return (i + (j << 12) + (n << 24) + (Jij << 36)); };
    auto unhash_key_ijnJ = [](size_t &i, size_t &j, size_t &n, size_t &Jij, size_t key)
    { i=(key & 0xFFFL);j=((key>>12)&0xFFFL);n=((key>>24)&0xFFFL);Jij=((key>>36)&0xFFFL); };                                                            //  0xF = 15 = 1111 (4bits), so 0xFFF is 12 bits of 1's. 0xFFFL makes it a long
    std::vector<std::unordered_map<size_t, size_t>> ket_lookup_pph; // in a given channel, map  i,j,n,Jij -> matrix index

    std::vector<std::vector<std::array<int, 4>>> ket_lookup_abc;
    std::vector<std::vector<double>> occs_abc_list;

    std::vector<arma::mat> Zbar;

    // we're looking at pph type states. We don't make a cut on the occupations in this channel
    // but if the first two orbits ij can't possibly make it past the occnat cut later on, we don't bother including them
    // To figure out if they'll make it, we want to know the biggest contribution that can possibly be made by a third orbit.
    double occnat_factor_max = 0;
    for (auto i : Z.modelspace->all_orbits)
    {
      double occnat_i = Z.modelspace->GetOrbit(i).occ_nat;
      occnat_factor_max = std::max(occnat_factor_max, occnat_i * (1 - occnat_i));
    }

    // Generate all the pph type channels J,parity,Tz. Note that these will be
    // not the same as the standard 3-body channels since we aren't enforcing antisymmetry with the 3rd particle
    int twoJph_min = 1;
    int twoJph_max = 6 * Z.modelspace->GetEmax() + 3;
    size_t nch_pph = 0;
    for (int twoJph = twoJph_min; twoJph <= twoJph_max; twoJph += 2)
    {
      for (int parity_ph = 0; parity_ph <= 1; parity_ph++)
      {
        for (int twoTz_ph = -3; twoTz_ph <= 3; twoTz_ph += 2)
        {
          //        std::vector<std::array<int,4>> good_kets_ijn;
          size_t ngood_ijn = 0;

          std::unordered_map<size_t, size_t> good_kets_ijn;
          //        std::map<std::array<int,4>, size_t> good_kets_ijn;

          std::vector<std::array<int, 4>> good_kets_abc;
          //        std::vector<size_t> kets_abc;
          std::vector<double> occs_abc;
          for (size_t chij = 0; chij < nch2; chij++)
          {
            TwoBodyChannel &tbc_ij = Z.modelspace->GetTwoBodyChannel(chij);
            int Jij = tbc_ij.J;
            int parity_ij = tbc_ij.parity;
            int Tz_ij = tbc_ij.Tz;
            //          for (size_t n : Z.modelspace->all_orbits )
            for (int n : Z.modelspace->all_orbits)
            {
              Orbit &on = Z.modelspace->GetOrbit(n);
              if ((on.l + parity_ij) % 2 != parity_ph)
                continue;
              if ((2 * Tz_ij - on.tz2) != twoTz_ph)
                continue; // note the minus sign because n is a hole
              if ((std::abs(2 * Jij - on.j2) > twoJph) or ((2 * Jij + on.j2) < twoJph))
                continue;

              size_t nkets_ij = tbc_ij.GetNumberKets();
              for (size_t iket_ij = 0; iket_ij < nkets_ij; iket_ij++)
              {

                // check that we pass our cuts on E3max, occupations, etc.
                Ket &ket_ij = tbc_ij.GetKet(iket_ij);
                //              size_t i = ket_ij.p;
                //              size_t j = ket_ij.q;
                int i = ket_ij.p;
                int j = ket_ij.q;
                double occnat_i = ket_ij.op->occ_nat;
                double occnat_j = ket_ij.oq->occ_nat;
                double d_ei = std::abs(2 * ket_ij.op->n + ket_ij.op->l - e_fermi[ket_ij.op->tz2]);
                double d_ej = std::abs(2 * ket_ij.oq->n + ket_ij.oq->l - e_fermi[ket_ij.oq->tz2]);
                // if i and j cant make it past the OccNat and dE3max cuts, don't bother including it
                if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_factor_max) < Z.modelspace->GetOccNat3Cut())
                  continue;
                if ((d_ei + d_ej) > Z.modelspace->dE3max)
                  continue;

                size_t key = hash_key_ijnJ(i, j, n, size_t(Jij));
                good_kets_ijn[key] = ngood_ijn;
                //              good_kets_ijn[ {i,j,n,Jij} ] = ngood_ijn;
                ngood_ijn++;
                double occ_factor = ket_ij.op->occ * ket_ij.oq->occ * (1 - on.occ) + (1 - ket_ij.op->occ) * (1 - ket_ij.oq->occ) * on.occ;

                if (i == j)
                  occ_factor *= 0.5; // because we only sum b<a
                if (std::abs(occ_factor) > 1e-8)
                {
                  good_kets_abc.push_back({i, j, n, Jij});
                  occs_abc.push_back(occ_factor);
                }
              }
            }
          }
          size_t ngood_abc = good_kets_abc.size();
          if ((ngood_ijn < 1) or (ngood_abc < 1))
            continue;

          channels_pph[{twoJph, parity_ph, twoTz_ph}] = nch_pph;
          channel_list_pph.push_back({twoJph, parity_ph, twoTz_ph});
          ket_lookup_pph.push_back(good_kets_ijn);
          ket_lookup_abc.push_back(good_kets_abc);
          occs_abc_list.push_back(occs_abc);
          Zbar.push_back(arma::mat(ngood_ijn, ngood_ijn, arma::fill::zeros));
          nch_pph++;
        } // for twoTz_ph
      }   // for parity_ph
    }     // for twoJph
          //  std::cout << "Nominally " << nch_pph << "  pph channels " << std::endl;

    Z.profiler.timer["_comm333_pph_setup_lists"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    // Now that we've made our lists and allocated the matrices, we fill the matrices inside the parallel block.
    // #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ch_pph = 0; ch_pph < nch_pph; ch_pph++)
    {
      int twoJph = channel_list_pph[ch_pph][0];
      double Jph = 0.5 * twoJph;
      auto &good_kets_ijn = ket_lookup_pph[ch_pph];
      auto &good_kets_abc = ket_lookup_abc[ch_pph];
      auto &occs_abc = occs_abc_list[ch_pph];
      size_t ngood_ijn = good_kets_ijn.size();
      size_t ngood_abc = good_kets_abc.size();
      arma::mat Xbar_ijnabc(ngood_ijn, ngood_abc, arma::fill::zeros);
      arma::mat Ybar_abcijn(ngood_abc, ngood_ijn, arma::fill::zeros);

      // now we fill Xbar an Ybar. This is where the heavy lifting is.
      for (auto &iter_ijn : good_kets_ijn)
      {
        size_t Iijn = iter_ijn.second;

        size_t key = iter_ijn.first;
        size_t i, j, n, Jij_tmp;
        unhash_key_ijnJ(i, j, n, Jij_tmp, key);
        int Jij = (int)Jij_tmp;
        //     size_t i = iter_ijn.first[0];
        //     size_t j = iter_ijn.first[1];
        //     size_t n = iter_ijn.first[2];
        //     int  Jij = iter_ijn.first[3];

        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        Orbit &on = Z.modelspace->GetOrbit(n);

        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;
        double occnat_n = on.occ_nat;
        double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
        double d_en = std::abs(2 * on.n + on.l - e_fermi[on.tz2]);
        double jn = 0.5 * on.j2;
        for (size_t Iabc = 0; Iabc < ngood_abc; Iabc++)
        {
          auto &abc_info = good_kets_abc[Iabc];
          size_t a = abc_info[0];
          size_t b = abc_info[1];
          size_t c = abc_info[2];
          double occ_factor = occs_abc[Iabc];
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          Orbit &oc = Z.modelspace->GetOrbit(c);
          double jc = 0.5 * oc.j2;
          int Jab = abc_info[3];

          double occnat_a = oa.occ_nat;
          double occnat_b = ob.occ_nat;
          double occnat_c = oc.occ_nat;
          double d_ea = std::abs(2 * oa.n + oa.l - e_fermi[oa.tz2]);
          double d_eb = std::abs(2 * ob.n + ob.l - e_fermi[ob.tz2]);
          double d_ec = std::abs(2 * oc.n + oc.l - e_fermi[oc.tz2]);

          if ((occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_n * (1 - occnat_n)) < Z.modelspace->GetOccNat3Cut())
            continue;
          if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_c * (1 - occnat_c)) < Z.modelspace->GetOccNat3Cut())
            continue;
          if ((d_ea + d_eb + d_en) > Z.modelspace->dE3max)
            continue;
          if ((d_ei + d_ej + d_ec) > Z.modelspace->dE3max)
            continue;

          double xbar = 0;
          double ybar = 0;
          int twoJprime_min = std::max(std::abs(2 * Jij - oc.j2), std::abs(2 * Jab - on.j2));
          int twoJprime_max = std::min((2 * Jij + oc.j2), (2 * Jab + on.j2));
          for (int twoJprime = twoJprime_min; twoJprime <= twoJprime_max; twoJprime += 2)
          {
            // We use this if/else construct because the precomputed SixJ symbols have
            // the largest half-integer on the bottom row.
            //         double sixj_ijn = Z.modelspace->GetSixJ( Jij, jn, Jph,  Jab, jc, 0.5*twoJprime);
            double sixj_ijn = (2 * Jph < twoJprime) ? Z.modelspace->GetSixJ(jn, Jph, Jij, jc, 0.5 * twoJprime, Jab)
                                                    : Z.modelspace->GetSixJ(jn, 0.5 * twoJprime, Jab, jc, Jph, Jij);

            // TODO This comes up often enough that we should just make ThreeBodyMEpn do this.
            std::vector<double> xandy = Y3.GetME_pn_TwoOps(Jij, Jab, twoJprime, i, j, c, a, b, n, X3, Y3);
            double xijcabn = xandy[0];
            double yijcabn = xandy[1];

            xbar += (twoJprime + 1) * sixj_ijn * xijcabn;
            ybar += (twoJprime + 1) * sixj_ijn * yijcabn;

          } // for twoJprime

          Xbar_ijnabc(Iijn, Iabc) = xbar;
          Ybar_abcijn(Iabc, Iijn) = hY * ybar * occ_factor;
        } // for Iabc
      }   // for iter_ijn

      auto &zbar = Zbar[ch_pph];
      zbar = Xbar_ijnabc * Ybar_abcijn;
      zbar -= hX * hY * zbar.t();

    } // for ch_pph , end of parallel block

    Z.profiler.timer["_comm333_pph_fill_matrices"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    // Now that we've computed Zbar (i.e. the pph transformed commutator), we need to
    // transform that back to Z, and do all the permutations to ensure antisymmetry
    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ch3 = 0; ch3 < nch3; ch3++)
    {
      auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      int twoJ = Tbc.twoJ;
      double Jtot = 0.5 * twoJ;
      for (size_t ibra = 0; ibra < nkets3; ibra++)
      {
        auto &bra = Tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        Orbit &ok = Z.modelspace->GetOrbit(k);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        double jk = 0.5 * ok.j2;
        double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
        double d_ek = std::abs(2 * ok.n + ok.l - e_fermi[ok.tz2]);
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;
        double occnat_k = ok.occ_nat;
        if ((d_ei + d_ej + d_ek) > Z.modelspace->dE3max)
          continue;
        if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
          continue;
        int J1 = bra.Jpq;

        // Here are the permutations for ijk
        std::vector<std::array<size_t, 3>> ijk = {{i, j, k}, {k, j, i}, {i, k, j}};

        std::vector<int> J1p_min = {J1, std::max(std::abs(ok.j2 - oj.j2), std::abs(twoJ - oi.j2)) / 2, std::max(std::abs(oi.j2 - ok.j2), std::abs(twoJ - oj.j2)) / 2};
        std::vector<int> J1p_max = {J1, std::min(ok.j2 + oj.j2, twoJ + oi.j2) / 2, std::min(oi.j2 + ok.j2, twoJ + oj.j2) / 2};
        //      std::vector<int> J1p_min = {J1,  std::abs(ok.j2-oj.j2)/2,   std::abs(oi.j2-ok.j2)/2 };
        //      std::vector<int> J1p_max = {J1,  (ok.j2+oj.j2)/2 , (oi.j2+ok.j2)/2 };
        std::vector<std::vector<double>> recouple_ijk = {{1}, {}, {}};
        for (int J1p = J1p_min[1]; J1p <= J1p_max[1]; J1p++)
          recouple_ijk[1].push_back(sqrt((2 * J1 + 1) * (2 * J1p + 1)) * Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p));

        for (int J1p = J1p_min[2]; J1p <= J1p_max[2]; J1p++)
          recouple_ijk[2].push_back(-Z.modelspace->phase((oj.j2 + ok.j2) / 2 + J1 + J1p) * sqrt((2 * J1 + 1) * (2 * J1p + 1)) * Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p));

        for (size_t iket = 0; iket <= ibra; iket++)
        {
          auto &ket = Tbc.GetKet(iket);
          size_t l = ket.p;
          size_t m = ket.q;
          size_t n = ket.r;
          Orbit &ol = Z.modelspace->GetOrbit(l);
          Orbit &om = Z.modelspace->GetOrbit(m);
          Orbit &on = Z.modelspace->GetOrbit(n);
          double jl = 0.5 * ol.j2;
          double jm = 0.5 * om.j2;
          double jn = 0.5 * on.j2;
          double d_el = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);
          double d_em = std::abs(2 * om.n + om.l - e_fermi[om.tz2]);
          double d_en = std::abs(2 * on.n + on.l - e_fermi[on.tz2]);
          double occnat_l = ol.occ_nat;
          double occnat_m = om.occ_nat;
          double occnat_n = on.occ_nat;
          if ((d_el + d_em + d_en) > Z.modelspace->dE3max)
            continue;
          if ((occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_n * (1 - occnat_n)) < Z.modelspace->GetOccNat3Cut())
            continue;
          int J2 = ket.Jpq;

          // permutations for lmn
          std::vector<std::array<size_t, 3>> lmn = {{l, m, n}, {n, m, l}, {l, n, m}};

          std::vector<int> J2p_min = {J2, std::max(std::abs(on.j2 - om.j2), std::abs(twoJ - ol.j2)) / 2, std::max(std::abs(ol.j2 - on.j2), std::abs(twoJ - om.j2)) / 2};
          std::vector<int> J2p_max = {J2, std::min(on.j2 + om.j2, twoJ + ol.j2) / 2, std::min(ol.j2 + on.j2, twoJ + om.j2) / 2};
          //        std::vector<int> J2p_min = {J2,  std::abs(on.j2-om.j2)/2,   std::abs(ol.j2-on.j2)/2 };
          //        std::vector<int> J2p_max = {J2,  (on.j2+om.j2)/2 , (ol.j2+on.j2)/2 };
          std::vector<std::vector<double>> recouple_lmn = {{1}, {}, {}};

          for (int J2p = J2p_min[1]; J2p <= J2p_max[1]; J2p++)
            recouple_lmn[1].push_back(sqrt((2 * J2 + 1) * (2 * J2p + 1)) * Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p));

          for (int J2p = J2p_min[2]; J2p <= J2p_max[2]; J2p++)
            recouple_lmn[2].push_back(-Z.modelspace->phase((om.j2 + on.j2) / 2 + J2 + J2p) * sqrt((2 * J2 + 1) * (2 * J2p + 1)) * Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p));

          double z_ijklmn = 0;

          for (int perm_ijk = 0; perm_ijk < 3; perm_ijk++)
          {
            size_t I1 = ijk[perm_ijk][0];
            size_t I2 = ijk[perm_ijk][1];
            size_t I3 = ijk[perm_ijk][2];
            Orbit &o1 = Z.modelspace->GetOrbit(I1);
            Orbit &o2 = Z.modelspace->GetOrbit(I2);
            Orbit &o3 = Z.modelspace->GetOrbit(I3);
            double j3 = 0.5 * o3.j2;
            for (int J1p = J1p_min[perm_ijk]; J1p <= J1p_max[perm_ijk]; J1p++)
            {
              if ((I1 == I2) and J1p % 2 > 0)
                continue;
              double rec_ijk = recouple_ijk[perm_ijk].at(J1p - J1p_min[perm_ijk]);
              for (int perm_lmn = 0; perm_lmn < 3; perm_lmn++)
              {
                size_t I4 = lmn[perm_lmn][0];
                size_t I5 = lmn[perm_lmn][1];
                size_t I6 = lmn[perm_lmn][2];
                Orbit &o4 = Z.modelspace->GetOrbit(I4);
                Orbit &o5 = Z.modelspace->GetOrbit(I5);
                Orbit &o6 = Z.modelspace->GetOrbit(I6);
                double j6 = 0.5 * o6.j2;

                int twoTz_ph = (o1.tz2 + o2.tz2 - o6.tz2);
                int parity_ph = (o1.l + o2.l + o6.l) % 2;

                for (int J2p = J2p_min[perm_lmn]; J2p <= J2p_max[perm_lmn]; J2p++)
                {
                  double rec_lmn = recouple_lmn[perm_lmn].at(J2p - J2p_min[perm_lmn]);
                  if ((I4 == I5) and J2p % 2 > 0)
                    continue;

                  size_t key_126 = hash_key_ijnJ(std::min(I1, I2), std::max(I1, I2), I6, J1p);
                  size_t key_453 = hash_key_ijnJ(std::min(I4, I5), std::max(I4, I5), I3, J2p);

                  int twoJph_min = std::max(std::abs(2 * J1p - o6.j2), std::abs(2 * J2p - o3.j2));
                  int twoJph_max = std::min((2 * J1p + o6.j2), (2 * J2p + o3.j2));
                  if (twoJph_min > twoJph_max)
                    continue;

                  //                std::vector<double> Zbar_126453 ( (twoJph_max-twoJph_min)/2+1, 0. );
                  int phase_12 = I1 > I2 ? -Z.modelspace->phase((o1.j2 + o2.j2) / 2 - J1p) : 1;
                  int phase_45 = I4 > I5 ? -Z.modelspace->phase((o4.j2 + o5.j2) / 2 - J2p) : 1;

                  double z_123456 = 0;
                  for (int twoJph = twoJph_min; twoJph <= twoJph_max; twoJph += 2)
                  {
                    auto iter_ch = channels_pph.find({twoJph, parity_ph, twoTz_ph});
                    if (iter_ch == channels_pph.end())
                      continue;
                    size_t ch_pph = iter_ch->second;

                    size_t index_126 = ket_lookup_pph[ch_pph].at(key_126);
                    size_t index_453 = ket_lookup_pph[ch_pph].at(key_453);
                    //                  size_t index_126 = ket_lookup_pph[ch_pph].at({ std::min(I1,I2),std::max(I1,I2),I6,J1p} );
                    //                  size_t index_453 = ket_lookup_pph[ch_pph].at({ std::min(I4,I5),std::max(I4,I5),I3,J2p} );
                    double zbar = phase_12 * phase_45 * Zbar[ch_pph](index_126, index_453);

                    double Jph = 0.5 * twoJph;
                    double sixj_ph = (Jph > Jtot) ? Z.modelspace->GetSixJ(j3, Jtot, J1p, j6, Jph, J2p) : Z.modelspace->GetSixJ(j3, Jph, J2p, j6, Jtot, J1p);
                    //                  double sixj_ph = Z.modelspace->GetSixJ( J1p, j3, Jtot,  J2p, j6, Jph );
                    //                  double sixj_ph;
                    //                  if (Jph>Jtot)   sixj_ph =  Z.modelspace->GetSixJ(  j3, Jtot,J1p , j6, Jph,J2p );
                    //                  else    sixj_ph = Z.modelspace->GetSixJ(  j3,Jph,J2p,  j6  , Jtot,J1p );
                    //                  double sixj_ph = Z.modelspace->GetSixJ(  j3, Jtot,J1p , j6, Jph,J2p );
                    z_123456 += (twoJph + 1) * sixj_ph * zbar;
                  } // for Jph
                  z_ijklmn += rec_ijk * rec_lmn * z_123456;

                } // for J2p
              }   // for perm_lmn
            }     // for J1p
          }       // for perm_ijk

          Z3.AddToME_pn_ch(ch3, ch3, ibra, iket, z_ijklmn);
        } // for iket
      }   // for ibra
    }     // for ch3

    Z.profiler.timer["_comm333_pph_recouple"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

  /*
  //
  // The old slow way, just to make sure we didn't screw anything up by going faster
  //
  void comm333_pph_hhpss_debug( const Operator& X, const Operator& Y, Operator& Z )
  {

    double tstart = omp_get_wtime();

    auto& X3 = X.ThreeBody;
    auto& Y3 = Y.ThreeBody;
    auto& Z3 = Z.ThreeBody;
    std::map<int,double> e_fermi = Z.modelspace->GetEFermi();

    size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();

    #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch3=0; ch3<nch3; ch3++)
    {
      auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      int twoJ = Tbc.twoJ;
      double Jtot = 0.5 * twoJ;
      for (size_t ibra=0; ibra<nkets3; ibra++)
      {
        auto& bra = Tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit& oi = Z.modelspace->GetOrbit(i);
        Orbit& oj = Z.modelspace->GetOrbit(j);
        Orbit& ok = Z.modelspace->GetOrbit(k);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        double jk = 0.5 * ok.j2;
        double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
        double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;
        double occnat_k = ok.occ_nat;
        if ( (d_ei + d_ej + d_ek) > Z.modelspace->dE3max ) continue;
        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < Z.modelspace->GetOccNat3Cut() ) continue;
        int J1 = bra.Jpq;



        std::vector<std::array<size_t,3>> ijk = { {i,j,k}, {k,j,i}, {i,k,j} };
        std::vector<int> J1p_min = {J1,  std::abs(ok.j2-oj.j2)/2,   std::abs(oi.j2-ok.j2)/2 };
        std::vector<int> J1p_max = {J1,  (ok.j2+oj.j2)/2 , (oi.j2+ok.j2)/2 };
        std::vector<std::vector<double>> recouple_ijk = {{1},{},{} };
        for (int J1p=J1p_min[1]; J1p<=J1p_max[1]; J1p++)
             recouple_ijk[1].push_back( sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(ji,jj,J1,jk,Jtot,J1p) );

        for (int J1p=J1p_min[2]; J1p<=J1p_max[2]; J1p++)
             recouple_ijk[2].push_back( -Z.modelspace->phase((oj.j2+ok.j2)/2+J1+J1p)*sqrt( (2*J1+1)*(2*J1p+1)) * Z.modelspace->GetSixJ(jj,ji,J1,jk,Jtot,J1p) );


  //      for (size_t iket=ibra; iket<nkets3; iket++)
        for (size_t iket=0; iket<=ibra; iket++)
        {
          auto& ket = Tbc.GetKet(iket);
          size_t l = ket.p;
          size_t m = ket.q;
          size_t n = ket.r;
          Orbit& ol = Z.modelspace->GetOrbit(l);
          Orbit& om = Z.modelspace->GetOrbit(m);
          Orbit& on = Z.modelspace->GetOrbit(n);
          double jl = 0.5 * ol.j2;
          double jm = 0.5 * om.j2;
          double jn = 0.5 * on.j2;
          double d_el = std::abs( 2*ol.n + ol.l - e_fermi[ol.tz2]);
          double d_em = std::abs( 2*om.n + om.l - e_fermi[om.tz2]);
          double d_en = std::abs( 2*on.n + on.l - e_fermi[on.tz2]);
          double occnat_l = ol.occ_nat;
          double occnat_m = om.occ_nat;
          double occnat_n = on.occ_nat;
          if ( (d_el + d_em + d_en) > Z.modelspace->dE3max ) continue;
          if ( (occnat_l*(1-occnat_l) * occnat_m*(1-occnat_m) * occnat_n*(1-occnat_n) ) < Z.modelspace->GetOccNat3Cut() ) continue;
          int J2 = ket.Jpq;


          std::vector<std::array<size_t,3>> lmn = { {l,m,n}, {n,m,l}, {l,n,m} };
          std::vector<int> J2p_min = {J2,  std::abs(on.j2-om.j2)/2,   std::abs(ol.j2-on.j2)/2 };
          std::vector<int> J2p_max = {J2,  (on.j2+om.j2)/2 , (ol.j2+on.j2)/2 };
          std::vector<std::vector<double>> recouple_lmn = {{1},{},{} };

          for (int J2p=J2p_min[1]; J2p<=J2p_max[1]; J2p++)
             recouple_lmn[1].push_back( sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jl,jm,J2,jn,Jtot,J2p) );

          for (int J2p=J2p_min[2]; J2p<=J2p_max[2]; J2p++)
               recouple_lmn[2].push_back( -Z.modelspace->phase((om.j2+on.j2)/2+J2+J2p)*sqrt( (2*J2+1)*(2*J2p+1)) * Z.modelspace->GetSixJ(jm,jl,J2,jn,Jtot,J2p) );


          double z_ijklmn = 0;


          for ( int perm_ijk=0; perm_ijk<3; perm_ijk++ )
          {
            size_t I1 = ijk[perm_ijk][0];
            size_t I2 = ijk[perm_ijk][1];
            size_t I3 = ijk[perm_ijk][2];
            Orbit& o1 = Z.modelspace->GetOrbit( I1 );
            Orbit& o2 = Z.modelspace->GetOrbit( I2 );
            Orbit& o3 = Z.modelspace->GetOrbit( I3 );
            double occnat_1 = o1.occ_nat;
            double occnat_2 = o2.occ_nat;
            double occnat_3 = o3.occ_nat;
            double d_e1 = std::abs( 2*o1.n + o1.l - e_fermi[o1.tz2]);
            double d_e2 = std::abs( 2*o2.n + o2.l - e_fermi[o2.tz2]);
            double d_e3 = std::abs( 2*o3.n + o3.l - e_fermi[o3.tz2]);
            double j3 = 0.5*o3.j2;
            for (int J1p=J1p_min[perm_ijk]; J1p<=J1p_max[perm_ijk]; J1p++)
            {
              double rec_ijk = recouple_ijk[perm_ijk].at(J1p-J1p_min[perm_ijk]);
              for ( int perm_lmn=0; perm_lmn<3; perm_lmn++ )
              {
                size_t I4 = lmn[perm_lmn][0];
                size_t I5 = lmn[perm_lmn][1];
                size_t I6 = lmn[perm_lmn][2];
                Orbit& o4 = Z.modelspace->GetOrbit( I4 );
                Orbit& o5 = Z.modelspace->GetOrbit( I5 );
                Orbit& o6 = Z.modelspace->GetOrbit( I6 );
                double occnat_4 = o4.occ_nat;
                double occnat_5 = o5.occ_nat;
                double occnat_6 = o6.occ_nat;
                double d_e4 = std::abs( 2*o4.n + o4.l - e_fermi[o4.tz2]);
                double d_e5 = std::abs( 2*o5.n + o5.l - e_fermi[o5.tz2]);
                double d_e6 = std::abs( 2*o6.n + o6.l - e_fermi[o6.tz2]);
                double j6 = 0.5*o6.j2;
                for (int J2p=J2p_min[perm_lmn]; J2p<=J2p_max[perm_lmn]; J2p++)
                {
                  double rec_lmn = recouple_lmn[perm_lmn].at(J2p-J2p_min[perm_lmn]);

                  for (size_t ch2=0; ch2<nch2; ch2++)
                  {
                    auto& tbc_ab = Z.modelspace->GetTwoBodyChannel(ch2);
                    if ( std::abs(Tbc.twoTz-2*tbc_ab.Tz)==5) continue; // TODO there are probably other checks at the channel level...
                    size_t nkets_ab = tbc_ab.GetNumberKets();
                    int Jab = tbc_ab.J;

                    for (size_t iket_ab=0; iket_ab<nkets_ab; iket_ab++)
                    {
                      Ket& ket_ab = tbc_ab.GetKet(iket_ab);
                      size_t a = ket_ab.p;
                      size_t b = ket_ab.q;
                      Orbit& oa = Z.modelspace->GetOrbit(a);
                      Orbit& ob = Z.modelspace->GetOrbit(b);
                      if (std::abs(oa.occ * ob.occ)<1e-6 and std::abs( (1-oa.occ)*(1-ob.occ))<1e-6) continue;

                      double occnat_a = oa.occ_nat;
                      double occnat_b = ob.occ_nat;
                      double d_ea = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
                      double d_eb = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
                      bool keep_ab3 = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_3*(1-occnat_3) ) >= Z.modelspace->GetOccNat3Cut();
                      bool keep_ab6 = (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_6*(1-occnat_6) ) >= Z.modelspace->GetOccNat3Cut();
                      keep_ab3 = keep_ab3 and (d_ea+d_eb+d_e3 <= Z.modelspace->dE3max );
                      keep_ab6 = keep_ab6 and (d_ea+d_eb+d_e6 <= Z.modelspace->dE3max );

                      if ( not ( keep_ab3 or keep_ab6) ) continue;

                      for (auto c : Z.modelspace->all_orbits)
                      {
                        Orbit& oc = Z.modelspace->GetOrbit(c);
                        double occ_factor = oa.occ * ob.occ * (1-oc.occ) + (1-oa.occ)*(1-ob.occ)*oc.occ;
                        if (std::abs(occ_factor)<1e-6) continue;
                        if (a==b) occ_factor *=0.5; // because we only sum b<a
                        double occnat_c = oc.occ_nat;
                        double d_ec = std::abs( 2*oc.n + oc.l - e_fermi[oc.tz2]);
                        double jc = 0.5 * oc.j2;
                        bool keep_12c = (occnat_1*(1-occnat_1) * occnat_2*(1-occnat_2) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut();
                        bool keep_45c = (occnat_4*(1-occnat_4) * occnat_5*(1-occnat_5) * occnat_c*(1-occnat_c) ) >= Z.modelspace->GetOccNat3Cut();
                        keep_12c = keep_12c and (d_e1+d_e2+d_ec <= Z.modelspace->dE3max );
                        keep_45c = keep_45c and (d_e4+d_e5+d_ec <= Z.modelspace->dE3max );


                        if ( ((oa.l+ob.l+o3.l+o4.l+o5.l+oc.l)%2==0) and ((o1.l+o2.l+oc.l+oa.l+ob.l+o6.l)%2==0)
                         and  ( (oa.tz2+ob.tz2+o3.tz2)==(o4.tz2+o5.tz2+oc.tz2) ) and ( (o1.tz2+o2.tz2+oc.tz2)==(oa.tz2+ob.tz2+o6.tz2))
                          and keep_ab3 and keep_45c and keep_12c and keep_ab6)
  //                        )//and keep_ab3 and keep_45c and keep_12c and keep_ab6)
                        {

                          int twoJx_min = std::max( std::abs(2*Jab - o3.j2), std::abs(2*J2p - oc.j2) );
                          int twoJx_max = std::min( 2*Jab+o3.j2 , 2*J2p + oc.j2 );
                          int twoJy_min = std::max( std::abs(2*J1p - oc.j2), std::abs(2*Jab - o6.j2) );
                          int twoJy_max = std::min( 2*J1p+oc.j2 , 2*Jab + o6.j2 );
                          if (twoJx_min <= twoJx_max and twoJy_min<=twoJy_max)
                          {
                            std::vector<double> xab345c( (twoJx_max-twoJx_min)/2+1, 0);
                            std::vector<double> yab345c( (twoJx_max-twoJx_min)/2+1, 0);
                            std::vector<double> x12cab6( (twoJy_max-twoJy_min)/2+1, 0);
                            std::vector<double> y12cab6( (twoJy_max-twoJy_min)/2+1, 0);
                            for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                            {
                              size_t iJx = (twoJx-twoJx_min)/2;
                              xab345c[iJx] = X3.GetME_pn(Jab,J2p,twoJx, a,b,I3,I4,I5,c);
                              yab345c[iJx] = Y3.GetME_pn(Jab,J2p,twoJx, a,b,I3,I4,I5,c);
                            }
                            for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                            {
                               size_t iJy = (twoJy-twoJy_min)/2;
                               x12cab6[iJy] = X3.GetME_pn(J1p,Jab,twoJy, I1,I2,c,a,b,I6);
                               y12cab6[iJy] = Y3.GetME_pn(J1p,Jab,twoJy, I1,I2,c,a,b,I6);
                            }
                            for (int twoJx=twoJx_min; twoJx<=twoJx_max; twoJx+=2)
                            {
                              double JJx = 0.5 * twoJx;
                              size_t iJx = (twoJx-twoJx_min)/2;
                              for (int twoJy=twoJy_min; twoJy<=twoJy_max; twoJy+=2)
                              {
                                 double JJy = 0.5 * twoJy;
                                 size_t iJy = (twoJy-twoJy_min)/2;
                                 double hats = (twoJx+1)*(twoJy+1);
                                 double ninej = Z.modelspace->GetNineJ( j3,Jab,JJx, J1p,JJy,jc, Jtot,j6,J2p);
                                  z_ijklmn +=  rec_ijk * rec_lmn * occ_factor * hats * ninej * ( xab345c[iJx]*y12cab6[iJy] - yab345c[iJx]*x12cab6[iJy] );
                              }// for twoJy
                            }// for twoJx
                          }
                        }// Z1 block

                       }// for c
                     }// for iket_ab
                   }// for ch2


                      }// for J2p
                    }// for perm_lmn
                  }// for J1p
                }// for perm_ijk

          Z3.AddToME_pn_ch( ch3, ch3, ibra, iket, z_ijklmn);
        }// for iket
      }// for ibra
    }//for ch3

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }
  */

  // void comm333_pph_hhpss( const Operator& X, const Operator& Y, Operator& Z )
  void comm333_pph_hhpss_debug(const Operator &X, const Operator &Y, Operator &Z)
  {

    double tstart = omp_get_wtime();
    if (Commutator::verbose)
      std::cout << __func__ << std::endl;

    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();

#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch3 = 0; ch3 < nch3; ch3++)
    {
      auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      int twoJ = Tbc.twoJ;
      double Jtot = 0.5 * twoJ;
      for (size_t ibra = 0; ibra < nkets3; ibra++)
      {
        auto &bra = Tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        Orbit &ok = Z.modelspace->GetOrbit(k);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        double jk = 0.5 * ok.j2;
        double d_ei = std::abs(2 * oi.n + oi.l - e_fermi[oi.tz2]);
        double d_ej = std::abs(2 * oj.n + oj.l - e_fermi[oj.tz2]);
        double d_ek = std::abs(2 * ok.n + ok.l - e_fermi[ok.tz2]);
        double occnat_i = oi.occ_nat;
        double occnat_j = oj.occ_nat;
        double occnat_k = ok.occ_nat;
        if ((d_ei + d_ej + d_ek) > Z.modelspace->dE3max)
          continue;
        if ((occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_k * (1 - occnat_k)) < Z.modelspace->GetOccNat3Cut())
          continue;
        int J1 = bra.Jpq;
        //      for (size_t iket=ibra; iket<nkets3; iket++)
        for (size_t iket = 0; iket <= ibra; iket++)
        {
          auto &ket = Tbc.GetKet(iket);
          size_t l = ket.p;
          size_t m = ket.q;
          size_t n = ket.r;
          Orbit &ol = Z.modelspace->GetOrbit(l);
          Orbit &om = Z.modelspace->GetOrbit(m);
          Orbit &on = Z.modelspace->GetOrbit(n);
          double jl = 0.5 * ol.j2;
          double jm = 0.5 * om.j2;
          double jn = 0.5 * on.j2;
          double d_el = std::abs(2 * ol.n + ol.l - e_fermi[ol.tz2]);
          double d_em = std::abs(2 * om.n + om.l - e_fermi[om.tz2]);
          double d_en = std::abs(2 * on.n + on.l - e_fermi[on.tz2]);
          double occnat_l = ol.occ_nat;
          double occnat_m = om.occ_nat;
          double occnat_n = on.occ_nat;
          if ((d_el + d_em + d_en) > Z.modelspace->dE3max)
            continue;
          if ((occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_n * (1 - occnat_n)) < Z.modelspace->GetOccNat3Cut())
            continue;
          int J2 = ket.Jpq;

          //              if ( not ((i==2 and j==4 and k==5 and l==4 and m==4 and n==5)
          //              if ( not ((i==4 and j==4 and k==5 and l==2 and m==4 and n==5)
          //               or (i==5 and j==4 and k==2 and l==5 and m==4 and n==4)) ) continue;
          double z_ijklmn = 0;

          for (size_t ch2 = 0; ch2 < nch2; ch2++)
          {
            auto &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch2);
            if (std::abs(Tbc.twoTz - 2 * tbc_ab.Tz) == 5)
              continue; // TODO there are probably other checks at the channel level...
            size_t nkets_ab = tbc_ab.GetNumberKets();
            int Jab = tbc_ab.J;

            for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
            {
              Ket &ket_ab = tbc_ab.GetKet(iket_ab);
              size_t a = ket_ab.p;
              size_t b = ket_ab.q;
              Orbit &oa = Z.modelspace->GetOrbit(a);
              Orbit &ob = Z.modelspace->GetOrbit(b);
              if (std::abs(oa.occ * ob.occ) < 1e-6 and std::abs((1 - oa.occ) * (1 - ob.occ)) < 1e-6)
                continue;

              double occnat_a = oa.occ_nat;
              double occnat_b = ob.occ_nat;
              double d_ea = std::abs(2 * oa.n + oa.l - e_fermi[oa.tz2]);
              double d_eb = std::abs(2 * ob.n + ob.l - e_fermi[ob.tz2]);
              bool keep_abi = (occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_i * (1 - occnat_i)) >= Z.modelspace->GetOccNat3Cut();
              bool keep_abj = (occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_j * (1 - occnat_j)) >= Z.modelspace->GetOccNat3Cut();
              bool keep_abk = (occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_k * (1 - occnat_k)) >= Z.modelspace->GetOccNat3Cut();
              bool keep_abl = (occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_l * (1 - occnat_l)) >= Z.modelspace->GetOccNat3Cut();
              bool keep_abm = (occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_m * (1 - occnat_m)) >= Z.modelspace->GetOccNat3Cut();
              bool keep_abn = (occnat_a * (1 - occnat_a) * occnat_b * (1 - occnat_b) * occnat_n * (1 - occnat_n)) >= Z.modelspace->GetOccNat3Cut();
              keep_abi = keep_abi and (d_ea + d_eb + d_ei <= Z.modelspace->dE3max);
              keep_abj = keep_abj and (d_ea + d_eb + d_ej <= Z.modelspace->dE3max);
              keep_abk = keep_abk and (d_ea + d_eb + d_ek <= Z.modelspace->dE3max);
              keep_abl = keep_abl and (d_ea + d_eb + d_el <= Z.modelspace->dE3max);
              keep_abm = keep_abm and (d_ea + d_eb + d_em <= Z.modelspace->dE3max);
              keep_abn = keep_abn and (d_ea + d_eb + d_en <= Z.modelspace->dE3max);

              for (auto c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double occ_factor = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
                if (std::abs(occ_factor) < 1e-6)
                  continue;
                if (a == b)
                  occ_factor *= 0.5; // because we only sum b<a
                double occnat_c = oc.occ_nat;
                double d_ec = std::abs(2 * oc.n + oc.l - e_fermi[oc.tz2]);
                double jc = 0.5 * oc.j2;
                bool keep_ijc = (occnat_i * (1 - occnat_i) * occnat_j * (1 - occnat_j) * occnat_c * (1 - occnat_c)) >= Z.modelspace->GetOccNat3Cut();
                bool keep_kjc = (occnat_k * (1 - occnat_k) * occnat_j * (1 - occnat_j) * occnat_c * (1 - occnat_c)) >= Z.modelspace->GetOccNat3Cut();
                bool keep_ikc = (occnat_i * (1 - occnat_i) * occnat_k * (1 - occnat_k) * occnat_c * (1 - occnat_c)) >= Z.modelspace->GetOccNat3Cut();
                bool keep_lmc = (occnat_l * (1 - occnat_l) * occnat_m * (1 - occnat_m) * occnat_c * (1 - occnat_c)) >= Z.modelspace->GetOccNat3Cut();
                bool keep_nmc = (occnat_n * (1 - occnat_n) * occnat_m * (1 - occnat_m) * occnat_c * (1 - occnat_c)) >= Z.modelspace->GetOccNat3Cut();
                bool keep_lnc = (occnat_l * (1 - occnat_l) * occnat_n * (1 - occnat_n) * occnat_c * (1 - occnat_c)) >= Z.modelspace->GetOccNat3Cut();
                keep_ijc = keep_ijc and (d_ei + d_ej + d_ec <= Z.modelspace->dE3max);
                keep_kjc = keep_kjc and (d_ek + d_ej + d_ec <= Z.modelspace->dE3max);
                keep_ikc = keep_ikc and (d_ei + d_ek + d_ec <= Z.modelspace->dE3max);
                keep_lmc = keep_lmc and (d_el + d_em + d_ec <= Z.modelspace->dE3max);
                keep_nmc = keep_nmc and (d_en + d_em + d_ec <= Z.modelspace->dE3max);
                keep_lnc = keep_lnc and (d_el + d_en + d_ec <= Z.modelspace->dE3max);

                // Direct Z1 term
                if (((oa.l + ob.l + ok.l + ol.l + om.l + oc.l) % 2 == 0) and ((oi.l + oj.l + oc.l + oa.l + ob.l + on.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + ok.tz2) == (ol.tz2 + om.tz2 + oc.tz2)) and ((oi.tz2 + oj.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + on.tz2)) and keep_abk and keep_lmc and keep_ijc and keep_abn)
                {
                  int twoJx_min = std::max(std::abs(2 * Jab - ok.j2), std::abs(2 * J2 - oc.j2));
                  int twoJx_max = std::min(2 * Jab + ok.j2, 2 * J2 + oc.j2);
                  int twoJy_min = std::max(std::abs(2 * J1 - oc.j2), std::abs(2 * Jab - on.j2));
                  int twoJy_max = std::min(2 * J1 + oc.j2, 2 * Jab + on.j2);
                  if (twoJx_min <= twoJx_max and twoJy_min <= twoJy_max)
                  {
                    std::vector<double> xabklmc((twoJx_max - twoJx_min) / 2 + 1, 0);
                    std::vector<double> yabklmc((twoJx_max - twoJx_min) / 2 + 1, 0);
                    std::vector<double> xijcabn((twoJy_max - twoJy_min) / 2 + 1, 0);
                    std::vector<double> yijcabn((twoJy_max - twoJy_min) / 2 + 1, 0);
                    for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                    {
                      //                  double JJx = 0.5 * twoJx;
                      size_t iJx = (twoJx - twoJx_min) / 2;
                      xabklmc[iJx] = X3.GetME_pn(Jab, J2, twoJx, a, b, k, l, m, c);
                      yabklmc[iJx] = Y3.GetME_pn(Jab, J2, twoJx, a, b, k, l, m, c);
                      //                  double xabklmc = X3.GetME_pn(Jab,J2,twoJx, a,b,k,l,m,c);
                      //                  double yabklmc = Y3.GetME_pn(Jab,J2,twoJx, a,b,k,l,m,c);
                    }
                    for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                    {
                      size_t iJy = (twoJy - twoJy_min) / 2;
                      xijcabn[iJy] = X3.GetME_pn(J1, Jab, twoJy, i, j, c, a, b, n);
                      yijcabn[iJy] = Y3.GetME_pn(J1, Jab, twoJy, i, j, c, a, b, n);
                    }
                    for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                    {
                      double JJx = 0.5 * twoJx;
                      size_t iJx = (twoJx - twoJx_min) / 2;
                      for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                      {
                        double JJy = 0.5 * twoJy;
                        size_t iJy = (twoJy - twoJy_min) / 2;
                        double hats = (twoJx + 1) * (twoJy + 1);
                        double ninej = Z.modelspace->GetNineJ(jk, Jab, JJx, J1, JJy, jc, Jtot, jn, J2);
                        //                     double xijcabn = X3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,n);
                        //                     double yijcabn = Y3.GetME_pn(J1,Jab,twoJy, i,j,c,a,b,n);
                        //                     z_ijklmn +=  occ_factor * hats * ninej * ( xabklmc*yijcabn - yabklmc*xijcabn );
                        z_ijklmn += occ_factor * hats * ninej * (xabklmc[iJx] * yijcabn[iJy] - yabklmc[iJx] * xijcabn[iJy]);
                      } // for twoJy
                    }   // for twoJx
                  }
                } // Z1 block
                  //              std::cout << "  after Z1  " << z_ijklmn << std::endl;

                // Z2 term  Pik    Xabilmc Ykjcabn
                if (((oa.l + ob.l + oi.l + ol.l + om.l + oc.l) % 2 == 0) and ((ok.l + oj.l + oc.l + oa.l + ob.l + on.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + oi.tz2) == (ol.tz2 + om.tz2 + oc.tz2)) and ((ok.tz2 + oj.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + on.tz2)) and keep_abi and keep_lmc and keep_kjc and keep_abn)
                {
                  int twoJx_min = std::max(std::abs(2 * Jab - oi.j2), std::abs(2 * J2 - oc.j2));
                  int twoJx_max = std::min(2 * Jab + oi.j2, 2 * J2 + oc.j2);
                  int J1p_min = std::abs(ok.j2 - oj.j2) / 2;
                  int J1p_max = (ok.j2 + oj.j2) / 2;
                  for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                  {
                    double JJx = 0.5 * twoJx;
                    double xabilmc = X3.GetME_pn(Jab, J2, twoJx, a, b, i, l, m, c);
                    double yabilmc = Y3.GetME_pn(Jab, J2, twoJx, a, b, i, l, m, c);

                    for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                    {
                      int twoJy_min = std::max(std::abs(2 * J1p - oc.j2), std::abs(2 * Jab - on.j2));
                      int twoJy_max = std::min(2 * J1p + oc.j2, 2 * Jab + on.j2);
                      double sixj = Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p);
                      if (std::abs(sixj) < 1e-6)
                        continue;
                      int phase = -1;
                      for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                      {
                        double JJy = 0.5 * twoJy;
                        double hats = (twoJx + 1) * (twoJy + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1));
                        double ninej = Z.modelspace->GetNineJ(ji, Jab, JJx, J1p, JJy, jc, Jtot, jn, J2);
                        double xkjcabn = X3.GetME_pn(J1p, Jab, twoJy, k, j, c, a, b, n);
                        double ykjcabn = Y3.GetME_pn(J1p, Jab, twoJy, k, j, c, a, b, n);
                        z_ijklmn -= occ_factor * phase * hats * sixj * ninej * (xabilmc * ykjcabn - yabilmc * xkjcabn);
                      } // for twoJy
                    }   // for J1p
                  }     // for twoJx
                }       // Z2 block
                        //              std::cout << "  after Z2  " << z_ijklmn << std::endl;

                // Z3 term  Pjk    Xabjlmc Yikcabn
                if (((oa.l + ob.l + oj.l + ol.l + om.l + oc.l) % 2 == 0) and ((oi.l + ok.l + oc.l + oa.l + ob.l + on.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + oj.tz2) == (ol.tz2 + om.tz2 + oc.tz2)) and ((oi.tz2 + ok.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + on.tz2)) and keep_abj and keep_lmc and keep_ikc and keep_abn)
                {
                  int twoJx_min = std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * J2 - oc.j2));
                  int twoJx_max = std::min(2 * Jab + oj.j2, 2 * J2 + oc.j2);
                  int J1p_min = std::abs(oi.j2 - ok.j2) / 2;
                  int J1p_max = (oi.j2 + ok.j2) / 2;
                  for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                  {
                    double JJx = 0.5 * twoJx;
                    double xabjlmc = X3.GetME_pn(Jab, J2, twoJx, a, b, j, l, m, c);
                    double yabjlmc = Y3.GetME_pn(Jab, J2, twoJx, a, b, j, l, m, c);

                    for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                    {
                      int twoJy_min = std::max(std::abs(2 * J1p - oc.j2), std::abs(2 * Jab - on.j2));
                      int twoJy_max = std::min(2 * J1p + oc.j2, 2 * Jab + on.j2);
                      double sixj = Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p);
                      if (std::abs(sixj) < 1e-6)
                        continue;
                      int phase = Z.modelspace->phase((oj.j2 + ok.j2) / 2 + J1 + J1p);
                      for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                      {
                        double JJy = 0.5 * twoJy;
                        double hats = (twoJx + 1) * (twoJy + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1));
                        double ninej = Z.modelspace->GetNineJ(jj, Jab, JJx, J1p, JJy, jc, Jtot, jn, J2);
                        double xikcabn = X3.GetME_pn(J1p, Jab, twoJy, i, k, c, a, b, n);
                        double yikcabn = Y3.GetME_pn(J1p, Jab, twoJy, i, k, c, a, b, n);
                        z_ijklmn -= occ_factor * phase * hats * sixj * ninej * (xabjlmc * yikcabn - yabjlmc * xikcabn);
                      } // for twoJy
                    }   // for J1p
                  }     // for twoJx
                }       // Z3 block
                        //              std::cout << "  after Z3  " << z_ijklmn << std::endl;

                // Z4 term  Pln    Xabknmc Yijcabl
                if (((oa.l + ob.l + ok.l + on.l + om.l + oc.l) % 2 == 0) and ((oi.l + oj.l + oc.l + oa.l + ob.l + ol.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + ok.tz2) == (on.tz2 + om.tz2 + oc.tz2)) and ((oi.tz2 + oj.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + ol.tz2)) and keep_abk and keep_nmc and keep_ijc and keep_abl)
                {
                  int twoJy_min = std::max(std::abs(2 * J1 - oc.j2), std::abs(2 * Jab - ol.j2));
                  int twoJy_max = std::min(2 * J1 + oc.j2, 2 * Jab + ol.j2);
                  int J2p_min = std::abs(on.j2 - om.j2) / 2;
                  int J2p_max = (on.j2 + om.j2) / 2;
                  for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                  {
                    double xijcabl = X3.GetME_pn(J1, Jab, twoJy, i, j, c, a, b, l);
                    double yijcabl = Y3.GetME_pn(J1, Jab, twoJy, i, j, c, a, b, l);
                    double JJy = 0.5 * twoJy;
                    for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                    {
                      double sixj = Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p);
                      int phase = -1;
                      if (std::abs(sixj) < 1e-6)
                        continue;
                      int twoJx_min = std::max(std::abs(2 * Jab - ok.j2), std::abs(2 * J2p - oc.j2));
                      int twoJx_max = std::min(2 * Jab + ok.j2, 2 * J2p + oc.j2);
                      for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                      {
                        double JJx = 0.5 * twoJx;
                        double xabknmc = X3.GetME_pn(Jab, J2p, twoJx, a, b, k, n, m, c);
                        double yabknmc = Y3.GetME_pn(Jab, J2p, twoJx, a, b, k, n, m, c);
                        double hats = (twoJx + 1) * (twoJy + 1) * sqrt((2 * J2 + 1) * (2 * J2p + 1));
                        double ninej = Z.modelspace->GetNineJ(jk, Jab, JJx, J1, JJy, jc, Jtot, jl, J2p);
                        z_ijklmn -= occ_factor * phase * hats * sixj * ninej * (xabknmc * yijcabl - yabknmc * xijcabl);
                      } // for twoJx
                    }   // for J2p
                  }     // for twoJy
                }       // Z4 block
                        //              std::cout << "  after Z4  " << z_ijklmn << std::endl;

                // Z5 term  Pmn    Xabklnc Yijcabm
                if (((oa.l + ob.l + ok.l + ol.l + on.l + oc.l) % 2 == 0) and ((oi.l + oj.l + oc.l + oa.l + ob.l + om.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + ok.tz2) == (ol.tz2 + on.tz2 + oc.tz2)) and ((oi.tz2 + oj.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + om.tz2)) and keep_abk and keep_lnc and keep_ijc and keep_abm)
                {
                  int twoJy_min = std::max(std::abs(2 * J1 - oc.j2), std::abs(2 * Jab - om.j2));
                  int twoJy_max = std::min(2 * J1 + oc.j2, 2 * Jab + om.j2);
                  int J2p_min = std::abs(ol.j2 - on.j2) / 2;
                  int J2p_max = (ol.j2 + on.j2) / 2;
                  for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                  {
                    double xijcabm = X3.GetME_pn(J1, Jab, twoJy, i, j, c, a, b, m);
                    double yijcabm = Y3.GetME_pn(J1, Jab, twoJy, i, j, c, a, b, m);
                    double JJy = 0.5 * twoJy;
                    for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                    {
                      double sixj = Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p);
                      int phase = Z.modelspace->phase((om.j2 + on.j2) / 2 + J2 + J2p);
                      if (std::abs(sixj) < 1e-6)
                        continue;
                      int twoJx_min = std::max(std::abs(2 * Jab - ok.j2), std::abs(2 * J2p - oc.j2));
                      int twoJx_max = std::min(2 * Jab + ok.j2, 2 * J2p + oc.j2);
                      for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                      {
                        double JJx = 0.5 * twoJx;
                        double xabklnc = X3.GetME_pn(Jab, J2p, twoJx, a, b, k, l, n, c);
                        double yabklnc = Y3.GetME_pn(Jab, J2p, twoJx, a, b, k, l, n, c);
                        double hats = (twoJx + 1) * (twoJy + 1) * sqrt((2 * J2 + 1) * (2 * J2p + 1));
                        double ninej = Z.modelspace->GetNineJ(jk, Jab, JJx, J1, JJy, jc, Jtot, jm, J2p);
                        z_ijklmn -= occ_factor * phase * hats * sixj * ninej * (xabklnc * yijcabm - yabklnc * xijcabm);
                      } // for twoJx
                    }   // for J2p
                  }     // for twoJy
                }       // Z5 block
                        //              std::cout << "  after Z5  " << z_ijklmn << std::endl;

                // Z6 term  Pik Pln    Xabinmc Ykjcabl
                if (((oa.l + ob.l + oi.l + on.l + om.l + oc.l) % 2 == 0) and ((ok.l + oj.l + oc.l + oa.l + ob.l + ol.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + oi.tz2) == (on.tz2 + om.tz2 + oc.tz2)) and ((ok.tz2 + oj.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + ol.tz2)) and keep_abi and keep_nmc and keep_kjc and keep_abl)
                {
                  int J1p_min = std::abs(ok.j2 - oj.j2) / 2;
                  int J1p_max = (ok.j2 + oj.j2) / 2;
                  int J2p_min = std::abs(on.j2 - om.j2) / 2;
                  int J2p_max = (on.j2 + om.j2) / 2;
                  for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                  {
                    double sixj2 = Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p);
                    if (std::abs(sixj2) < 1e-6)
                      continue;
                    int twoJx_min = std::max(std::abs(2 * Jab - oi.j2), std::abs(2 * J2p - oc.j2));
                    int twoJx_max = std::min(2 * Jab + oi.j2, 2 * J2p + oc.j2);
                    for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                    {
                      double JJx = 0.5 * twoJx;
                      double xabinmc = X3.GetME_pn(Jab, J2p, twoJx, a, b, i, n, m, c);
                      double yabinmc = Y3.GetME_pn(Jab, J2p, twoJx, a, b, i, n, m, c);

                      for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                      {
                        int twoJy_min = std::max(std::abs(2 * J1p - oc.j2), std::abs(2 * Jab - ol.j2));
                        int twoJy_max = std::min(2 * J1p + oc.j2, 2 * Jab + ol.j2);
                        double sixj1 = Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p);
                        if (std::abs(sixj1) < 1e-6)
                          continue;
                        int phase = 1;
                        for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                        {
                          //                         if (twoJy==5) continue;
                          double JJy = 0.5 * twoJy;
                          double hats = (twoJx + 1) * (twoJy + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1) * (2 * J2 + 1) * (2 * J2p + 1));
                          double ninej = Z.modelspace->GetNineJ(ji, Jab, JJx, J1p, JJy, jc, Jtot, jl, J2p);
                          double xkjcabl = X3.GetME_pn(J1p, Jab, twoJy, k, j, c, a, b, l);
                          double ykjcabl = Y3.GetME_pn(J1p, Jab, twoJy, k, j, c, a, b, l);
                          z_ijklmn += occ_factor * phase * hats * sixj1 * sixj2 * ninej * (xabinmc * ykjcabl - yabinmc * xkjcabl);
                        } // for twoJy
                      }   // for J1p
                    }     // for twoJx
                  }       // for J2p
                }         // Z6 block
                          //              std::cout << "  after Z6  " << z_ijklmn << std::endl;

                // Z7 term  Pjk Pln    Xabjnmc Yikcabl
                if (((oa.l + ob.l + oj.l + on.l + om.l + oc.l) % 2 == 0) and ((oi.l + ok.l + oc.l + oa.l + ob.l + ol.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + oj.tz2) == (on.tz2 + om.tz2 + oc.tz2)) and ((oi.tz2 + ok.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + ol.tz2)) and keep_abj and keep_nmc and keep_ikc and keep_abl)
                {
                  int J1p_min = std::abs(oi.j2 - ok.j2) / 2;
                  int J1p_max = (oi.j2 + ok.j2) / 2;
                  int J2p_min = std::abs(on.j2 - om.j2) / 2;
                  int J2p_max = (on.j2 + om.j2) / 2;
                  for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                  {
                    double sixj2 = Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p);
                    if (std::abs(sixj2) < 1e-6)
                      continue;
                    int twoJx_min = std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * J2p - oc.j2));
                    int twoJx_max = std::min(2 * Jab + oj.j2, 2 * J2p + oc.j2);
                    for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                    {
                      double JJx = 0.5 * twoJx;
                      double xabjnmc = X3.GetME_pn(Jab, J2p, twoJx, a, b, j, n, m, c);
                      double yabjnmc = Y3.GetME_pn(Jab, J2p, twoJx, a, b, j, n, m, c);

                      for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                      {
                        int twoJy_min = std::max(std::abs(2 * J1p - oc.j2), std::abs(2 * Jab - ol.j2));
                        int twoJy_max = std::min(2 * J1p + oc.j2, 2 * Jab + ol.j2);
                        double sixj1 = Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p);
                        if (std::abs(sixj1) < 1e-6)
                          continue;
                        int phase = -Z.modelspace->phase((oj.j2 + ok.j2) / 2 + J1 + J1p);
                        for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                        {
                          double JJy = 0.5 * twoJy;
                          double hats = (twoJx + 1) * (twoJy + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1) * (2 * J2 + 1) * (2 * J2p + 1));
                          double ninej = Z.modelspace->GetNineJ(jj, Jab, JJx, J1p, JJy, jc, Jtot, jl, J2p);
                          double xikcabl = X3.GetME_pn(J1p, Jab, twoJy, i, k, c, a, b, l);
                          double yikcabl = Y3.GetME_pn(J1p, Jab, twoJy, i, k, c, a, b, l);
                          z_ijklmn += occ_factor * phase * hats * sixj1 * sixj2 * ninej * (xabjnmc * yikcabl - yabjnmc * xikcabl);
                        } // for twoJy
                      }   // for J1p
                    }     // for twoJx
                  }       // for J2p
                }         // Z7 block
                          //              std::cout << "  after Z7  " << z_ijklmn << std::endl;

                // Z8 term  Pik Pmn    Xabilnc Ykjcabm
                if (((oa.l + ob.l + oi.l + ol.l + on.l + oc.l) % 2 == 0) and ((ok.l + oj.l + oc.l + oa.l + ob.l + om.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + oi.tz2) == (ol.tz2 + on.tz2 + oc.tz2)) and ((ok.tz2 + oj.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + om.tz2)) and keep_abi and keep_lnc and keep_kjc and keep_abm)
                {
                  int J1p_min = std::abs(ok.j2 - oj.j2) / 2;
                  int J1p_max = (ok.j2 + oj.j2) / 2;
                  int J2p_min = std::abs(ol.j2 - on.j2) / 2;
                  int J2p_max = (ol.j2 + on.j2) / 2;
                  for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                  {
                    double sixj2 = Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p);
                    if (std::abs(sixj2) < 1e-6)
                      continue;
                    int twoJx_min = std::max(std::abs(2 * Jab - oi.j2), std::abs(2 * J2p - oc.j2));
                    int twoJx_max = std::min(2 * Jab + oj.j2, 2 * J2p + oc.j2);
                    for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                    {
                      double JJx = 0.5 * twoJx;
                      double xabilnc = X3.GetME_pn(Jab, J2p, twoJx, a, b, i, l, n, c);
                      double yabilnc = Y3.GetME_pn(Jab, J2p, twoJx, a, b, i, l, n, c);

                      for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                      {
                        int twoJy_min = std::max(std::abs(2 * J1p - oc.j2), std::abs(2 * Jab - om.j2));
                        int twoJy_max = std::min(2 * J1p + oc.j2, 2 * Jab + om.j2);
                        double sixj1 = Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p);
                        if (std::abs(sixj1) < 1e-6)
                          continue;
                        int phase = -Z.modelspace->phase((om.j2 + on.j2) / 2 + J2 + J2p);
                        for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                        {
                          double JJy = 0.5 * twoJy;
                          double hats = (twoJx + 1) * (twoJy + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1) * (2 * J2 + 1) * (2 * J2p + 1));
                          double ninej = Z.modelspace->GetNineJ(ji, Jab, JJx, J1p, JJy, jc, Jtot, jm, J2p);
                          double xkjcabm = X3.GetME_pn(J1p, Jab, twoJy, k, j, c, a, b, m);
                          double ykjcabm = Y3.GetME_pn(J1p, Jab, twoJy, k, j, c, a, b, m);
                          z_ijklmn += occ_factor * phase * hats * sixj1 * sixj2 * ninej * (xabilnc * ykjcabm - yabilnc * xkjcabm);
                        } // for twoJy
                      }   // for J1p
                    }     // for twoJx
                  }       // for J2p
                }         // Z8 block
                          //              std::cout << "  after Z8  " << z_ijklmn << std::endl;

                // Z9 term  Pjk Pmn    Xabjlnc Yikcabm
                if (((oa.l + ob.l + oj.l + ol.l + on.l + oc.l) % 2 == 0) and ((oi.l + ok.l + oc.l + oa.l + ob.l + om.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + oj.tz2) == (ol.tz2 + on.tz2 + oc.tz2)) and ((oi.tz2 + ok.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + om.tz2)) and keep_abj and keep_lnc and keep_ikc and keep_abm)
                {
                  int J1p_min = std::abs(oi.j2 - ok.j2) / 2;
                  int J1p_max = (oi.j2 + ok.j2) / 2;
                  int J2p_min = std::abs(ol.j2 - on.j2) / 2;
                  int J2p_max = (ol.j2 + on.j2) / 2;
                  for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                  {
                    double sixj2 = Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p);
                    if (std::abs(sixj2) < 1e-6)
                      continue;
                    int twoJx_min = std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * J2p - oc.j2));
                    int twoJx_max = std::min(2 * Jab + oj.j2, 2 * J2p + oc.j2);
                    for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                    {
                      double JJx = 0.5 * twoJx;
                      double xabjlnc = X3.GetME_pn(Jab, J2p, twoJx, a, b, j, l, n, c);
                      double yabjlnc = Y3.GetME_pn(Jab, J2p, twoJx, a, b, j, l, n, c);

                      for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
                      {
                        int twoJy_min = std::max(std::abs(2 * J1p - oc.j2), std::abs(2 * Jab - om.j2));
                        int twoJy_max = std::min(2 * J1p + oc.j2, 2 * Jab + om.j2);
                        double sixj1 = Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p);
                        if (std::abs(sixj1) < 1e-6)
                          continue;
                        int phase = Z.modelspace->phase((oj.j2 + ok.j2 + om.j2 + on.j2) / 2 + J1 + J1p + J2 + J2p);
                        for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                        {
                          double JJy = 0.5 * twoJy;
                          double hats = (twoJx + 1) * (twoJy + 1) * sqrt((2 * J1 + 1) * (2 * J1p + 1) * (2 * J2 + 1) * (2 * J2p + 1));
                          double ninej = Z.modelspace->GetNineJ(jj, Jab, JJx, J1p, JJy, jc, Jtot, jm, J2p);
                          double xikcabm = X3.GetME_pn(J1p, Jab, twoJy, i, k, c, a, b, m);
                          double yikcabm = Y3.GetME_pn(J1p, Jab, twoJy, i, k, c, a, b, m);
                          z_ijklmn += occ_factor * phase * hats * sixj1 * sixj2 * ninej * (xabjlnc * yikcabm - yabjlnc * xikcabm);
                        } // for twoJy
                      }   // for J1p
                    }     // for twoJx
                  }       // for J2p
                }         // Z9 block
                          //              std::cout << "  after Z9  " << z_ijklmn << std::endl;

              } // for c
            }   // for iket_ab
          }     // for ch2
                //        std::cout << "  " << i << " " << j << " " << k << " " << l << " " << m << " " << n << "  J1 J2 two J " << J1 << " " << J2 << " " << twoJ << "   Z = " << z_ijklmn << std::endl;
          Z3.AddToME_pn_ch(ch3, ch3, ibra, iket, z_ijklmn);
        } // for iket
      }   // for ibra
    }     // for ch3

    Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
  }

}// namespace Commutator
