//=====================================================================================================================
//  comm332_ppph_hhhpss, scalar X and isospin-changing Y. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
// comm332_ppph_hhhpss for a scalar X and a channel-changing Y (J rank 0, parity 0, T rank > 0), precompiled bookkeeping.
// Same maths as comm332_ppph_hhhpss_tz: per 2B channel c, intermediates L(c) = {(abc, d)}, X_c(ij, q) = X_{ijd,abc};
//   Z(cb, ck) = hY X_cb F_cb Y_{ck<-L(cb)}^T - hX Y_{cb<-L(ck)} F_ck X_ck^T,   Y_{c'<-L(c)}(kl, q) = Y_{kld,abc}
// with X_{ijd,abc} = (R X)(P(ij, d), abc) and Y_{kld,abc} = (R Y)(P(kl, d), abc) from the unpack caches.
// Included after PerfPlan332p.inc and PerfPlan333tz.inc (inside namespace Commutator).
void comm332_ppph_hhhpss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced())
  {
    comm332_ppph_hhhpss_tz(X, Y, Z);
    return;
  }
  auto &Z2 = Z.TwoBody;
  if (not(X.ThreeBody.IsAllocated() and Y.ThreeBody.IsAllocated())) return;
  int hX = X.IsHermitian() ? 1 : -1, hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(Z, Z);
  const std::vector<arma::mat> &RX = perfplan::unpack_rx(X, P);
  const perfplan::WOffCache &WY = perfplan::unpack_off(Y, P);
  const std::vector<arma::mat> WYrk = perfplan::unpack_off_rk(WY, P);   // R(cb) Y, held while this term runs
  int nch = ms.GetNumberTwoBodyChannels();
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  int norb = P.norb;
  struct L { std::vector<int> ch3, ik, d; std::vector<double> fac; arma::mat XM; std::vector<int> rows; };
  std::vector<L> Ls(nch);
  // Z2 is a [3B, 3B] contraction: (ij d) and (kl d) are 3B states, so only kets with both orbits and orbits d with
  // e <= emax_3body contribute (rows: those kets of the channel)
  auto in3 = [&](size_t o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l <= ms.GetEMax3Body(); };

#pragma omp parallel for schedule(dynamic, 1)
  for (int ch = 0; ch < nch; ch++)
  {
    TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
    auto &Lc = Ls[ch];
    for (size_t ib = 0; ib < tbc.GetNumberKets(); ib++) if (in3(tbc.GetKet(ib).p) and in3(tbc.GetKet(ib).q)) Lc.rows.push_back(ib);
    int J = tbc.J, nkets = Lc.rows.size();
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
          Lc.ch3.push_back(ch_abc); Lc.ik.push_back(ik); Lc.d.push_back(d);
          Lc.fac.push_back((twoJ + 1.) / (2 * J + 1) * occf * symm / 6.);
        }
      }
    }
    size_t dim = Lc.ch3.size();
    Lc.XM.zeros(nkets, dim);
    for (size_t q = 0; q < dim; q++)
    {
      size_t c3 = Lc.ch3[q], ik = Lc.ik[q], d = Lc.d[q];
      if (RX[c3].n_elem == 0) continue;
      for (int ib = 0; ib < nkets; ib++)
      {
        Ket &bra = tbc.GetKet(Lc.rows[ib]);
        int g = P.g2idx[((size_t)bra.p * norb + bra.q) * (P.J2max + 1) + J];
        if (g < 0) continue;
        int pidx = P.pidx(c3, g, d);
        if (pidx < 0) continue;
        Lc.XM(ib, q) = ((bra.p == bra.q) ? 1. / sqrt(2) : 1.) * RX[c3](pidx, ik);
      }
    }
  }
  // Y_{c' <- L(c)}(kl, q) * fac_q
  auto Yfrom = [&](int cp, int c) -> arma::mat
  {
    TwoBodyChannel &tbc = ms.GetTwoBodyChannel(cp);
    const std::vector<int> &rows = Ls[cp].rows;
    int J = tbc.J, nkets = rows.size();
    auto &Lc = Ls[c];
    size_t dim = Lc.ch3.size();
    arma::mat YM(nkets, dim, arma::fill::zeros);
    for (size_t q = 0; q < dim; q++)
    {
      int c3 = Lc.ch3[q], ik = Lc.ik[q], d = Lc.d[q];
      Orbit &od = ms.GetOrbit(d);
      int twoJ = ms.GetThreeBodyChannel(c3).twoJ;
      for (int ib = 0; ib < nkets; ib++)
      {
        Ket &bra = tbc.GetKet(rows[ib]);
        int g = P.g2idx[((size_t)bra.p * norb + bra.q) * (P.J2max + 1) + J];
        if (g < 0) continue;
        Orbit &o1 = ms.GetOrbit(bra.p); Orbit &o2 = ms.GetOrbit(bra.q);
        int c3p = P.ch3tab[P.i3(twoJ, (o1.l + o2.l + od.l) % 2, o1.tz2 + o2.tz2 + od.tz2)];
        if (c3p < 0) continue;
        int k = WY.idx[(size_t)c3p * WY.nch + c3];
        if (k < 0) continue;
        int pidx = P.pidx(c3p, g, d);
        if (pidx < 0) continue;
        YM(ib, q) = ((bra.p == bra.q) ? 1. / sqrt(2) : 1.) * Lc.fac[q] * WYrk[k](pidx, ik);
      }
    }
    return YM;
  };
  std::vector<std::array<size_t, 2>> blocks;
  for (auto &it : Z2.MatEl) blocks.push_back({it.first[0], it.first[1]});
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < blocks.size(); ib++)
  {
    size_t cb = blocks[ib][0], ck = blocks[ib][1];
    if (cb == ck) continue;
    if (ms.GetTwoBodyChannel(cb).J != ms.GetTwoBodyChannel(ck).J) continue;
    size_t nb = Ls[cb].rows.size(), nk = Ls[ck].rows.size();
    if (nb == 0 or nk == 0) continue;
    arma::mat ZM(nb, nk, arma::fill::zeros);
    if (Ls[cb].ch3.size() > 0) ZM += hY * Ls[cb].XM * Yfrom(ck, cb).t();
    if (Ls[ck].ch3.size() > 0) ZM -= hX * Yfrom(cb, ck) * Ls[ck].XM.t();
    for (size_t i = 0; i < nb; i++)
      for (size_t k = 0; k < nk; k++) Z2.AddToTBME(cb, ck, Ls[cb].rows[i], Ls[ck].rows[k], ZM(i, k));
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
