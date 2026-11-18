//=====================================================================================================================
//  comm332_ppph_hhhpss with precompiled bookkeeping (scalar X, Y). Notation: see the reader's guide in PerfPlan.inc.
//=====================================================================================================================
//  Formula (IMSRG3Commutators.cc):
//    Z_ijkl = 1/6 sum_abcd (n'a n'b n'c nd - na nb nc n'd) (X_{ijd,abc} Y_{abc,kld} - Y_{ijd,abc} X_{abc,kld})
//    Z^J_ijkl = 1/6 sum_abcd (...) 1/(2J+1) sum_{J1 J'} (2J'+1) (X^{J J1 J'}_{ijd,abc} Y^{J1 J J'}_{abc,kld} - X <-> Y)
//  Per 2B channel (J) the sum over intermediates q = (3B ket abc in channel J', orbit d) is a GEMM:
//    Step 1  list the q with nonzero occupation, weight f_q = (2J'+1)/(2J+1) occ_q symm_abc/6
//    Step 2  XM(ij, q) = X_{(ij)J d, abc} = (R X)(P(gij, d), abc)      (half-unpacked R X from the cache: the bra (ij; d)
//            needs recoupling to stored order, the ket abc is a stored ket), YM(ij, q) = f_q (R Y)(P(gij, d), abc)
//    Step 3  ZM = hY XM YM^T - hX YM XM^T, added to Z2 (upper triangle).
void comm332_ppph_hhhpss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  auto &X3 = X.ThreeBody; auto &Y3 = Y.ThreeBody; auto &Z2 = Z.TwoBody;
  if (not(X3.IsAllocated() and Y3.IsAllocated())) return;
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0)
  {
    comm332_ppph_hhhpss(X, Y, Z);
    return;
  }
  int hX = X.IsHermitian() ? 1 : -1, hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(X, Z);
  const std::vector<arma::mat> &RX = perfplan::unpack_rx(X, P);
  const std::vector<arma::mat> &RY = perfplan::unpack_rx(Y, P);
  ModelSpace &ms = *Z.modelspace;
  int nch = ms.GetNumberTwoBodyChannels();
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  int norb = P.norb;
  // Z2 is a [3B, 3B] contraction: (ij d) and (kl d) are 3B states, so only kets with both orbits and orbits d with
  // e <= emax_3body contribute (rows: those kets of the channel)
  auto in3 = [&](size_t o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l <= ms.GetEMax3Body(); };
#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ch++)
  {
    TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
    int J = tbc.J;
    std::vector<int> rows;
    for (size_t ib = 0; ib < tbc.GetNumberKets(); ib++) if (in3(tbc.GetKet(ib).p) and in3(tbc.GetKet(ib).q)) rows.push_back(ib);
    int nkets = rows.size();
    std::vector<size_t> ch3_list, iket_list, d_list; std::vector<double> fac;   // Step 1
    for (size_t ch_abc = 0; ch_abc < nch3; ch_abc++)
    {
      ThreeBodyChannel &T = ms.GetThreeBodyChannel(ch_abc);
      if (std::abs(T.twoTz - 2 * tbc.Tz) > 1) continue;
      int twoJ = T.twoJ;
      for (size_t ik = 0; ik < T.GetNumberKets(); ik++)
      {
        Ket3 &k3 = T.GetKet(ik);
        size_t a = k3.p, b = k3.q, c = k3.r;
        double occ = k3.op->occ * k3.oq->occ * k3.oR->occ, occb = (1 - k3.op->occ) * (1 - k3.oq->occ) * (1 - k3.oR->occ);
        if (occ == 0 and occb == 0) continue;
        double symm = 6;
        if (a == b and b == c) symm = 1; else if (a == b or a == c or b == c) symm = 3;
        for (auto d : ms.all_orbits)
        {
          if (not in3(d)) continue;
          Orbit &od = ms.GetOrbit(d);
          double occf = occb * od.occ - occ * (1 - od.occ);
          if (std::abs(occf) < 1e-6) continue;
          if (std::abs(2 * J - od.j2) > twoJ or 2 * J + od.j2 < twoJ) continue;
          if ((od.l + tbc.parity + T.parity) % 2 > 0) continue;
          if (T.twoTz != od.tz2 + 2 * tbc.Tz) continue;
          ch3_list.push_back(ch_abc); iket_list.push_back(ik); d_list.push_back(d);
          fac.push_back((twoJ + 1.) / (2 * J + 1) * occf * symm / 6.);
        }
      }
    }
    size_t dim = ch3_list.size();   // Step 2
    arma::mat XM(nkets, dim, arma::fill::zeros), YM(nkets, dim, arma::fill::zeros);
    for (size_t q = 0; q < dim; q++)
    {
      size_t c3 = ch3_list[q], ik = iket_list[q], d = d_list[q];
      if (RX[c3].n_elem == 0) continue;
      for (int ib = 0; ib < nkets; ib++)
      {
        Ket &bra = tbc.GetKet(rows[ib]);
        int g = P.g2idx[((size_t)bra.p * norb + bra.q) * (P.J2max + 1) + J];
        if (g < 0) continue;
        int pidx = P.pidx(c3, g, d);
        if (pidx < 0) continue;
        double norm2b = (bra.p == bra.q) ? 1. / sqrt(2) : 1;
        XM(ib, q) = norm2b * RX[c3](pidx, ik);
        YM(ib, q) = norm2b * fac[q] * RY[c3](pidx, ik);
      }
    }
    arma::mat ZM = hY * XM * YM.t() - hX * YM * XM.t();   // Step 3
    for (int ib = 0; ib < nkets; ib++)
      for (int ik = ib; ik < nkets; ik++)
        Z2.AddToTBME(ch, ch, rows[ib], rows[ik], ZM(ib, ik));
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
