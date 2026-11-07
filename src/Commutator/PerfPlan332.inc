//=====================================================================================================================
//  comm332_pphhss with precompiled bookkeeping (scalar X, Y, default cuts). Notation: see the reader's guide in PerfPlan.inc.
//=====================================================================================================================
//  Formula (the pphh part of the [3,3] -> 2 commutator, IMSRG3Commutators.cc):
//    Z_ijkl = 1/4 (1 - Pij)(1 - Pkl) sum_{abcd} na nb n'c n'd (X_{abi,cdk} Y_{cdj,abl} - Y_{abi,cdk} X_{cdj,abl})
//  evaluated, as in the original, in cross-coupled channels (Jph, parity, |dTz|):
//    Step 1  unpack W_X, W_Y (cached).
//    Step 2  Xbar(pq, abcd) = sum_J' S(p, q, J') W_X^{J'}(P(gab, p), P(gcd, q)),   S = (2J'+1) {Jcd Jab Jph; jp jq J'} (-1)^{(jp + J')}
//            Ybar(abcd, pq) = occ(abcd) sum_J' S(p, q, J') W_Y^{J'}(P(gcd, q), P(gab, p))   (rows p q in both orders)
//    Step 3  Zbar = Xbar Ybar for the four row/column orderings (GEMMs).
//    Step 4  back-transformation into Z2 with {ji jj J; jk jl Jph} 6j symbols: unchanged from comm332_pphhss.
//  Only step 2 differs from the original: the 3B Pandya sums read W instead of GetME_pn. W_Y(q, p) is taken as
//  hY W_Y(p, q) (same column as the W_X read) for Hermitian or anti-Hermitian Y.
//  Included after PerfPlan233.inc (uses Plan233's 6j table).
namespace perfplan
{
  struct Plan332
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1;
    std::vector<std::vector<std::array<int, 2>>> abcd;   // ch_CC -> (g_ab, g_cd)
    std::vector<std::vector<double>> fac;                 // occupation_factor * symmetry_factor
    // Z2 here is a [3B, 3B] contraction, so all four external orbits and the pairs (ab), (cd) lie in the 3B space
    // (e <= emax_3body). rows3: the cross-coupled kets with both orbits in the 3B space; rowOf: CC ket -> row or -1.
    std::vector<std::vector<int>> rows3, rowOf;
  };
  static Plan332 &get332(const Operator &Z, Plan &P)
  {
    static Plan332 Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ) return Q;
    Q = Plan332(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max();
    size_t nch_CC = ms.TwoBodyChannels_CC.size(), nch = ms.GetNumberTwoBodyChannels();
    int norb = P.norb;
    Q.abcd.resize(nch_CC); Q.fac.resize(nch_CC); Q.rows3.resize(nch_CC); Q.rowOf.resize(nch_CC);
    auto in3 = [&](int o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l <= ms.GetEMax3Body(); };
    for (size_t ch = 0; ch < nch_CC; ch++)
    {
      auto &tbc_CC = ms.GetTwoBodyChannel_CC(ch);
      int J_ph = tbc_CC.J;
      Q.rowOf[ch].assign(tbc_CC.GetNumberKets(), -1);
      for (size_t ib = 0; ib < tbc_CC.GetNumberKets(); ib++)
      {
        Ket &kt = tbc_CC.GetKet(ib);
        if (in3(kt.p) and in3(kt.q)) { Q.rowOf[ch][ib] = Q.rows3[ch].size(); Q.rows3[ch].push_back(ib); }
      }
      for (size_t ch_ab = 0; ch_ab < nch; ch_ab++)        // same order as the original
      {
        auto &tbc_ab = ms.GetTwoBodyChannel(ch_ab);
        for (size_t ch_cd = 0; ch_cd < nch; ch_cd++)
        {
          auto &tbc_cd = ms.GetTwoBodyChannel(ch_cd);
          if (std::abs(tbc_ab.J - tbc_cd.J) > J_ph or (tbc_ab.J + tbc_cd.J) < J_ph) continue;
          if ((tbc_ab.parity + tbc_cd.parity + tbc_CC.parity) % 2 > 0) continue;
          if (std::abs(tbc_cd.Tz - tbc_ab.Tz) != tbc_CC.Tz) continue;
          for (size_t iab = 0; iab < tbc_ab.GetNumberKets(); iab++)
          {
            Ket &kab = tbc_ab.GetKet(iab);
            double na = kab.op->occ, nb = kab.oq->occ;
            if (na * nb < 1e-8 and (1 - na) * (1 - nb) < 1e-8) continue;
            int gab = P.g2idx[((size_t)kab.p * norb + kab.q) * (P.J2max + 1) + tbc_ab.J];
            if (gab >= P.ng3) continue;   // (ab) outside the 3B space
            for (size_t icd = 0; icd < tbc_cd.GetNumberKets(); icd++)
            {
              Ket &kcd = tbc_cd.GetKet(icd);
              double nc = kcd.op->occ, nd = kcd.oq->occ;
              double occf = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);
              if (std::abs(occf) < 1e-6) continue;
              double sym = 1;
              if (kab.p == kab.q) sym *= 0.5;
              if (kcd.p == kcd.q) sym *= 0.5;
              int gcd = P.g2idx[((size_t)kcd.p * norb + kcd.q) * (P.J2max + 1) + tbc_cd.J];
              if (gcd >= P.ng3) continue;   // (cd) outside the 3B space
              Q.abcd[ch].push_back({gab, gcd}); Q.fac[ch].push_back(occf * sym);
            }
          }
        }
      }
    }
    return Q;
  }
} // namespace perfplan

void comm332_pphhss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0)
  {
    comm332_pphhss(X, Y, Z);
    return;
  }
  auto &X3 = X.ThreeBody; auto &Y3 = Y.ThreeBody; auto &Z2 = Z.TwoBody;
  if (not(X3.IsAllocated() and Y3.IsAllocated())) return;
  perfplan::Plan &P = perfplan::get(X, Z);
  perfplan::Plan233 &S = perfplan::get233(Z, P);   // 6j table {a b c; d e f} of shape (2J, 2J, 2J, j2, j2, tJ)
  perfplan::Plan332 &Q = perfplan::get332(Z, P);
  ModelSpace &ms = *Z.modelspace;
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  int norb = P.norb;
  double t = omp_get_wtime();
  const std::vector<perfplan::wmat> &WX = perfplan::unpack(X, P);
  const std::vector<perfplan::wmat> &WY = perfplan::unpack(Y, P);
  Z.profiler.timer["_comm332_pphh_plan_unpack"] += omp_get_wtime() - t; t = omp_get_wtime();
  // W_Y = R Y R^T is (anti)symmetric like Y: W_Y(p2, p1) = hY W_Y(p1, p2), read in the same column as W_X(p1, p2)
  const bool symY = Y.IsHermitian() or Y.IsAntiHermitian();
  const double hY = Y.IsHermitian() ? 1.0 : -1.0;
  auto eo = [&](int o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l; };
  int emx = 0;
  for (int o = 0; o < norb; o++) emx = std::max(emx, eo(o));
  const bool ecut = P.E3max < 3 * emx;
  size_t nch_CC = ms.TwoBodyChannels_CC.size();
  std::deque<arma::mat> Z_bar(nch_CC), Z_bar_flipbra(nch_CC), Z_bar_flipket(nch_CC), Z_bar_flipboth(nch_CC);
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ch = 0; ch < nch_CC; ++ch)
  {
    TwoBodyChannel_CC &tbc_CC = ms.GetTwoBodyChannel_CC(ch);
    const std::vector<int> &rows = Q.rows3[ch];
    size_t nKets_cc = rows.size();   // rows: CC kets with both orbits in the 3B space
    int J_ph = tbc_CC.J;
    auto &abcd = Q.abcd[ch]; auto &fac = Q.fac[ch];
    size_t nab = abcd.size();
    arma::mat X3_ij(nKets_cc, nab, arma::fill::zeros), X3_ji(nKets_cc, nab, arma::fill::zeros);
    arma::mat Y3_ij(nab, nKets_cc, arma::fill::zeros), Y3_ji(nab, nKets_cc, arma::fill::zeros);
    std::vector<int> eab(nab), ecd(nab);   // pair energies, for the E3max pre-check
    for (size_t s = 0; s < nab; s++) { eab[s] = eo(P.g2[abcd[s][0]][0]) + eo(P.g2[abcd[s][0]][1]); ecd[s] = eo(P.g2[abcd[s][1]][0]) + eo(P.g2[abcd[s][1]][1]); }
    for (size_t ib = 0; ib < nKets_cc; ib++)   // Step 2: Pandya fill from W
    {
      Ket &bra = tbc_CC.GetKet(rows[ib]);
      int i = bra.p, j = bra.q;
      Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
      const int ei = eo(i), ej = eo(j);
      for (size_t s = 0; s < nab; s++)
      {
        // E3max < 3 emax: the ij term needs P(gab, i) and P(gcd, j), the ji term P(gab, j) and P(gcd, i)
        if (ecut and (eab[s] + ei > P.E3max or ecd[s] + ej > P.E3max) and (eab[s] + ej > P.E3max or ecd[s] + ei > P.E3max)) continue;
        int gab = abcd[s][0], gcd = abcd[s][1];
        int a = P.g2[gab][0], b = P.g2[gab][1], Jab = P.g2[gab][2], Jcd = P.g2[gcd][2];
        Orbit &oa = ms.GetOrbit(a); Orbit &ob = ms.GetOrbit(b);
        int twoJp_min = std::min(std::max(std::abs(2 * Jab - oi.j2), std::abs(2 * Jcd - oj.j2)), std::max(std::abs(2 * Jab - oj.j2), std::abs(2 * Jcd - oi.j2)));
        int twoJp_max = std::max(std::min(2 * Jab + oi.j2, 2 * Jcd + oj.j2), std::min(2 * Jab + oj.j2, 2 * Jcd + oi.j2));
        double xij = 0, yij = 0, xji = 0, yji = 0;
        int par_abi = (oa.l + ob.l + oi.l) % 2, tz_abi = oa.tz2 + ob.tz2 + oi.tz2;
        int par_abj = (oa.l + ob.l + oj.l) % 2, tz_abj = oa.tz2 + ob.tz2 + oj.tz2;
        for (int tJp = twoJp_min; tJp <= twoJp_max; tJp += 2)
        {
          double s_ij = (tJp + 1) * S.six6(2 * Jcd, 2 * Jab, 2 * J_ph, oi.j2, oj.j2, tJp) * ms.phase((oi.j2 + tJp) / 2);
          double s_ji = (tJp + 1) * S.six6(2 * Jcd, 2 * Jab, 2 * J_ph, oj.j2, oi.j2, tJp) * ms.phase((oj.j2 + tJp) / 2);
          if (s_ij != 0)
          {
            int c3 = P.ch3tab[P.i3(tJp, par_abi, tz_abi)];
            if (c3 >= 0)
            {
              int p1 = P.pidx(c3, gab, i), p2 = P.pidx(c3, gcd, j);
              if (p1 >= 0 and p2 >= 0) { xij += s_ij * WX[c3](p1, p2); yij += s_ij * (symY ? hY * WY[c3](p1, p2) : WY[c3](p2, p1)); }
            }
          }
          if (s_ji != 0)
          {
            int c3 = P.ch3tab[P.i3(tJp, par_abj, tz_abj)];
            if (c3 >= 0)
            {
              int p1 = P.pidx(c3, gab, j), p2 = P.pidx(c3, gcd, i);
              if (p1 >= 0 and p2 >= 0) { xji += s_ji * WX[c3](p1, p2); yji += s_ji * (symY ? hY * WY[c3](p1, p2) : WY[c3](p2, p1)); }
            }
          }
        }
        X3_ij(ib, s) = xij; Y3_ij(s, ib) = yij * fac[s];
        X3_ji(ib, s) = xji; Y3_ji(s, ib) = yji * fac[s];
      }
    }
    // Step 3
    Z_bar[ch] = X3_ij * Y3_ij; Z_bar_flipbra[ch] = X3_ji * Y3_ij; Z_bar_flipket[ch] = X3_ij * Y3_ji; Z_bar_flipboth[ch] = X3_ji * Y3_ji;
  }
  Z.profiler.timer["_comm332_pphh_plan_fill_gemm"] += omp_get_wtime() - t; t = omp_get_wtime();
  size_t nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ch = 0; ch < nch; ch++)      // Step 4: back-transformation, unchanged from comm332_pphhss
  {
    TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
    int J = tbc.J;
    size_t nkets = tbc.GetNumberKets();
    for (size_t ibra = 0; ibra < nkets; ibra++)
    {
      Ket &bra = tbc.GetKet(ibra);
      size_t i = bra.p, j = bra.q;
      Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
      if (2 * oi.n + oi.l > ms.GetEMax3Body() or 2 * oj.n + oj.l > ms.GetEMax3Body()) continue;   // zero: i, j outside the 3B space
      double ji = 0.5 * oi.j2, jj = 0.5 * oj.j2;
      int phase_ij = ms.phase((oi.j2 + oj.j2) / 2 - J);
      for (size_t iket = ibra; iket < nkets; iket++)
      {
        Ket &ket = tbc.GetKet(iket);
        size_t k = ket.p, l = ket.q;
        Orbit &ok = ms.GetOrbit(k); Orbit &ol = ms.GetOrbit(l);
        if (2 * ok.n + ok.l > ms.GetEMax3Body() or 2 * ol.n + ol.l > ms.GetEMax3Body()) continue;   // zero: k, l outside the 3B space
        double jk = 0.5 * ok.j2, jl = 0.5 * ol.j2;
        int phase_kl = ms.phase((ok.j2 + ol.j2) / 2 - J);
        double zijkl = 0;
        int parity_cc = (oi.l + ol.l) % 2, Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        int J_ph_min = std::max(std::abs(oi.j2 - ol.j2), std::abs(ok.j2 - oj.j2)) / 2, J_ph_max = std::min(oi.j2 + ol.j2, ok.j2 + oj.j2) / 2;
        for (int J_ph = J_ph_min; J_ph <= J_ph_max; J_ph++)
        {
          int ch_cc = ms.GetTwoBodyChannelIndex(J_ph, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = ms.GetTwoBodyChannel_CC(ch_cc);
          int indx_il = Q.rowOf[ch_cc][tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l))], indx_kj = Q.rowOf[ch_cc][tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j))];
          double sixj = ms.GetSixJ(ji, jj, J, jk, jl, J_ph);
          double zbar_ilkj, zbar_jkli;
          if (i <= l and k <= j) { zbar_ilkj = Z_bar[ch_cc](indx_il, indx_kj); zbar_jkli = Z_bar_flipboth[ch_cc](indx_kj, indx_il); }
          else if (i > l and k <= j) { zbar_ilkj = Z_bar_flipbra[ch_cc](indx_il, indx_kj); zbar_jkli = Z_bar_flipbra[ch_cc](indx_kj, indx_il); }
          else if (i <= l and k > j) { zbar_ilkj = Z_bar_flipket[ch_cc](indx_il, indx_kj); zbar_jkli = Z_bar_flipket[ch_cc](indx_kj, indx_il); }
          else { zbar_ilkj = Z_bar_flipboth[ch_cc](indx_il, indx_kj); zbar_jkli = Z_bar[ch_cc](indx_kj, indx_il); }
          zijkl += (2 * J_ph + 1) * sixj * (zbar_ilkj + phase_ij * phase_kl * zbar_jkli);
        }
        parity_cc = (oj.l + ol.l) % 2; Tz_cc = std::abs(oj.tz2 - ol.tz2) / 2;
        J_ph_min = std::max(std::abs(oj.j2 - ol.j2), std::abs(ok.j2 - oi.j2)) / 2; J_ph_max = std::min(oj.j2 + ol.j2, ok.j2 + oi.j2) / 2;
        for (int J_ph = J_ph_min; J_ph <= J_ph_max; J_ph++)
        {
          int ch_cc = ms.GetTwoBodyChannelIndex(J_ph, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = ms.GetTwoBodyChannel_CC(ch_cc);
          int indx_jl = Q.rowOf[ch_cc][tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l))], indx_ki = Q.rowOf[ch_cc][tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i))];
          double sixj = ms.GetSixJ(jj, ji, J, jk, jl, J_ph);
          double zbar_jlki = 0, zbar_iklj = 0;
          if (j <= l and k <= i) { zbar_jlki = Z_bar[ch_cc](indx_jl, indx_ki); zbar_iklj = Z_bar_flipboth[ch_cc](indx_ki, indx_jl); }
          else if (j <= l and k > i) { zbar_jlki = Z_bar_flipket[ch_cc](indx_jl, indx_ki); zbar_iklj = Z_bar_flipket[ch_cc](indx_ki, indx_jl); }
          else if (j > l and k <= i) { zbar_jlki = Z_bar_flipbra[ch_cc](indx_jl, indx_ki); zbar_iklj = Z_bar_flipbra[ch_cc](indx_ki, indx_jl); }
          else { zbar_jlki = Z_bar_flipboth[ch_cc](indx_jl, indx_ki); zbar_iklj = Z_bar[ch_cc](indx_ki, indx_jl); }
          zijkl -= (2 * J_ph + 1) * sixj * (phase_ij * zbar_jlki + phase_kl * zbar_iklj);
        }
        zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
        Z2.AddToTBME(ch, ch, ibra, iket, zijkl);
      }
    }
  }
  Z.profiler.timer["_comm332_pphh_plan_back"] += omp_get_wtime() - t;
  perfplan::release_memory();
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
