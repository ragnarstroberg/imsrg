//=====================================================================================================================
//  comm332_pphhss, scalar X and isospin-changing Y. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
// comm332_pphhss for a scalar X and a channel-changing Y (J rank 0, parity 0, T rank > 0), precompiled bookkeeping.
// Same maths as comm332_pphhss_tz (signed rows, both orderings, Zbar(row channel, d, Delta) = Xbar Ybar), with the
// 3B Pandya sums read from the unpacked W (X diagonal, Y off-diagonal blocks), the bookkeeping cached, and table
// lookups in the back-transformation. Included after PerfPlan223.inc, PerfPlan233.inc and PerfPlan333tz.inc.
namespace perfplan
{
  struct Plan332tz
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1, TY = -1, pY = -1;
    std::vector<std::array<std::vector<std::array<int, 2>>, 3>> rows;   // ch_cc -> [d+1] -> ordered (p, q)
    std::vector<std::vector<int>> rowpos;                               // ch_cc -> CC local index (0..2n-1) -> position
    struct Abcd { int gab, gcd; double factor; };
    std::vector<std::array<std::vector<Abcd>, 3>> mids;                  // ch_cc -> [d+1]
    struct Item { int ch_row, d, idelta, ch_col, dp; };
    std::vector<Item> items;
    std::vector<std::array<long, 9>> item_of;
    std::vector<std::array<int, 2>> tasks;                               // (ch_row, d+1)
  };
  static Plan332tz &get332tz(const Operator &Z, Plan &P, Plan233 &S, int TY, int pY)
  {
    static Plan332tz Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ and Q.TY == TY and Q.pY == pY) return Q;
    Q = Plan332tz(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max(); Q.TY = TY; Q.pY = pY;
    size_t nch = ms.GetNumberTwoBodyChannels(), nch_CC = ms.TwoBodyChannels_CC.size();
    int norb = P.norb;
    Q.rows.resize(nch_CC); Q.rowpos.resize(nch_CC); Q.mids.resize(nch_CC); Q.item_of.resize(nch_CC);
    for (size_t ch = 0; ch < nch_CC; ch++)
    {
      auto &tbc = ms.GetTwoBodyChannel_CC(ch);
      size_t n = tbc.GetNumberKets();
      Q.rowpos[ch].assign(2 * n, -1);
      for (size_t r = 0; r < 2 * n; r++)
      {
        Ket &ket = tbc.GetKet(r % n);
        int p = (r < n) ? ket.p : ket.q, q = (r < n) ? ket.q : ket.p;
        // [3B, 3B] -> 2B: all external orbits lie in the 3B space
        if (2 * ms.GetOrbit(p).n + ms.GetOrbit(p).l > ms.GetEMax3Body() or 2 * ms.GetOrbit(q).n + ms.GetOrbit(q).l > ms.GetEMax3Body()) continue;
        int d = (ms.GetOrbit(p).tz2 - ms.GetOrbit(q).tz2) / 2;
        Q.rowpos[ch][r] = Q.rows[ch][d + 1].size();
        Q.rows[ch][d + 1].push_back({p, q});
      }
      int J_ph = tbc.J;
      for (size_t ch_ab = 0; ch_ab < nch; ch_ab++)
      {
        auto &tab = ms.GetTwoBodyChannel(ch_ab);
        for (size_t ch_cd = 0; ch_cd < nch; ch_cd++)
        {
          auto &tcd = ms.GetTwoBodyChannel(ch_cd);
          if (std::abs(tab.J - tcd.J) > J_ph or tab.J + tcd.J < J_ph) continue;
          if ((tab.parity + tcd.parity + tbc.parity) % 2 > 0) continue;
          int dsign = tcd.Tz - tab.Tz;
          if (std::abs(dsign) != tbc.Tz) continue;
          for (size_t iab = 0; iab < tab.GetNumberKets(); iab++)
          {
            Ket &kab = tab.GetKet(iab);
            double na = kab.op->occ, nb = kab.oq->occ;
            if (na * nb < 1e-8 and (1 - na) * (1 - nb) < 1e-8) continue;
            int gab = P.g2idx[((size_t)kab.p * norb + kab.q) * (P.J2max + 1) + tab.J];
            if (gab >= P.ng3) continue;   // (ab) outside the 3B space
            for (size_t icd = 0; icd < tcd.GetNumberKets(); icd++)
            {
              Ket &kcd = tcd.GetKet(icd);
              double nc = kcd.op->occ, nd = kcd.oq->occ;
              double occf = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);
              if (std::abs(occf) < 1e-6) continue;
              double sym = 1;
              if (kab.p == kab.q) sym *= 0.5;
              if (kcd.p == kcd.q) sym *= 0.5;
              int gcd = P.g2idx[((size_t)kcd.p * norb + kcd.q) * (P.J2max + 1) + tcd.J];
              if (gcd >= P.ng3) continue;   // (cd) outside the 3B space
              Q.mids[ch][dsign + 1].push_back({gab, gcd, occf * sym});
            }
          }
        }
      }
    }
    for (size_t ch = 0; ch < nch_CC; ch++)
    {
      Q.item_of[ch].fill(-1);
      auto &tbc = ms.GetTwoBodyChannel_CC(ch);
      for (int d = -1; d <= 1; d++)
      {
        if (Q.rows[ch][d + 1].empty() or Q.mids[ch][d + 1].empty()) continue;
        bool any = false;
        for (int idelta = 0; idelta < (TY == 0 ? 1 : 2); idelta++)
        {
          int dp = d - (idelta == 0 ? TY : -TY);
          if (std::abs(dp) > 1) continue;
          int cc = S.chph[(size_t)tbc.J * 4 + ((tbc.parity + pY) % 2) * 2 + std::abs(dp)];
          if (cc < 0 or Q.rows[cc][dp + 1].empty()) continue;
          Q.item_of[ch][(d + 1) * 3 + idelta] = Q.items.size();
          Q.items.push_back({(int)ch, d, idelta, cc, dp});
          any = true;
        }
        if (any) Q.tasks.push_back({(int)ch, d + 1});
      }
    }
    return Q;
  }
} // namespace perfplan

void comm332_pphhss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced())
  {
    comm332_pphhss_tz(X, Y, Z);
    return;
  }
  if (not(X.ThreeBody.IsAllocated() and Y.ThreeBody.IsAllocated())) return;
  int TY = Y.GetTRank();
  perfplan::Plan &P = perfplan::get(Z, Z);
  perfplan::Plan233 &S = perfplan::get233(Z, P);
  perfplan::Plan223 &S6 = perfplan::get223(Z, P);   // 6j {j3 J J12; j6 ja J45} table, used as {ji jj J; jk jl Jph}
  perfplan::Plan332tz &Q = perfplan::get332tz(Z, P, S, TY, Y.GetParity());
  int norb = P.norb;
  double t = omp_get_wtime();
  const std::vector<perfplan::wmat> &WX = perfplan::unpack(X, P);
  const perfplan::WOffCache &WY = perfplan::unpack_off(Y, P);
  // sum_J' (2J'+1) {Jcd Jab Jph; jp jq J'} (-1)^((jp+J')) O_{(ab)p, (cd)q}   (isY: O_{(cd)q, (ab)p})
  auto pandya = [&](bool isY, int J_ph, int p, int q, const perfplan::Plan332tz::Abcd &m) -> double
  {
    Orbit &op = ms.GetOrbit(p); Orbit &oq = ms.GetOrbit(q);
    int Jab = P.g2[m.gab][2], Jcd = P.g2[m.gcd][2];
    Orbit &oa = ms.GetOrbit(P.g2[m.gab][0]); Orbit &ob = ms.GetOrbit(P.g2[m.gab][1]);
    Orbit &oc = ms.GetOrbit(P.g2[m.gcd][0]); Orbit &od = ms.GetOrbit(P.g2[m.gcd][1]);
    int par_abp = (oa.l + ob.l + op.l) % 2, tz_abp = oa.tz2 + ob.tz2 + op.tz2;
    int par_cdq = (oc.l + od.l + oq.l) % 2, tz_cdq = oc.tz2 + od.tz2 + oq.tz2;
    int lo = std::max(std::abs(2 * Jab - op.j2), std::abs(2 * Jcd - oq.j2)), hi = std::min(2 * Jab + op.j2, 2 * Jcd + oq.j2);
    double sum = 0;
    for (int tJp = lo; tJp <= hi; tJp += 2)
    {
      double s6 = S.six6(2 * Jcd, 2 * Jab, 2 * J_ph, op.j2, oq.j2, tJp);
      if (s6 == 0) continue;
      int c1 = P.ch3tab[P.i3(tJp, par_abp, tz_abp)], c2 = P.ch3tab[P.i3(tJp, par_cdq, tz_cdq)];
      if (c1 < 0 or c2 < 0) continue;
      int p1 = P.pidx(c1, m.gab, p), p2 = P.pidx(c2, m.gcd, q);
      if (p1 < 0 or p2 < 0) continue;
      double me;
      if (isY) { int k = WY.idx[(size_t)c2 * WY.nch + c1]; if (k < 0) continue; me = WY.mats[k](p2, p1); }
      else { if (c1 != c2) continue; me = WX[c1](p1, p2); }
      sum += (tJp + 1) * s6 * ms.phase((op.j2 + tJp) / 2) * me;
    }
    return sum;
  };
  std::vector<arma::mat> Zbar(Q.items.size());
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t it = 0; it < Q.tasks.size(); it++)
  {
    int ch_row = Q.tasks[it][0], d = Q.tasks[it][1] - 1;
    int J_ph = ms.GetTwoBodyChannel_CC(ch_row).J;
    auto &rows_r = Q.rows[ch_row][d + 1];
    auto &mid = Q.mids[ch_row][d + 1];
    arma::mat Xbar(rows_r.size(), mid.size());
    for (size_t im = 0; im < mid.size(); im++)
      for (size_t ir = 0; ir < rows_r.size(); ir++) Xbar(ir, im) = pandya(false, J_ph, rows_r[ir][0], rows_r[ir][1], mid[im]);
    for (int idelta = 0; idelta < 2; idelta++)
    {
      long itm = Q.item_of[ch_row][(d + 1) * 3 + idelta];
      if (itm < 0) continue;
      auto &cols_c = Q.rows[Q.items[itm].ch_col][Q.items[itm].dp + 1];
      arma::mat Ybar(mid.size(), cols_c.size());
      for (size_t ic = 0; ic < cols_c.size(); ic++)
        for (size_t im = 0; im < mid.size(); im++) Ybar(im, ic) = mid[im].factor * pandya(true, J_ph, cols_c[ic][0], cols_c[ic][1], mid[im]);
      Zbar[itm] = Xbar * Ybar;
    }
  }
  Z.profiler.timer["_comm332_pphh_plantz_fill_gemm"] += omp_get_wtime() - t; t = omp_get_wtime();
  auto zbar = [&](int J_ph, int p, int q, int r, int s) -> double
  {
    Orbit &op = ms.GetOrbit(p); Orbit &oq = ms.GetOrbit(q); Orbit &orr = ms.GetOrbit(r); Orbit &os = ms.GetOrbit(s);
    int d = (op.tz2 - oq.tz2) / 2, dp = (orr.tz2 - os.tz2) / 2;
    int Delta = d - dp, idelta = (Delta == TY) ? 0 : ((Delta == -TY) ? 1 : -1);
    if (idelta < 0) return 0.0;
    int ch_row = S.chph[(size_t)J_ph * 4 + ((op.l + oq.l) % 2) * 2 + std::abs(d)];
    if (ch_row < 0) return 0.0;
    long itm = Q.item_of[ch_row][(d + 1) * 3 + idelta];
    if (itm < 0) return 0.0;
    int ch_col = Q.items[itm].ch_col;
    int kr = S.ccLocal[ch_row][(size_t)p * norb + q], kc = S.ccLocal[ch_col][(size_t)r * norb + s];
    if (kr < 0 or kc < 0) return 0.0;
    long ir = Q.rowpos[ch_row][kr], ic = Q.rowpos[ch_col][kc];
    if (ir < 0 or ic < 0) return 0.0;
    return Zbar[itm](ir, ic);
  };
  auto six = [&](int ji, int jj, int J, int jk, int jl, int Jph) { return S6.six6(ji, jj, 2 * J, jk, jl, 2 * Jph); };
  auto &Z2 = Z.TwoBody;
  std::vector<std::array<size_t, 3>> bk;
  for (auto &it : Z2.MatEl)
    for (size_t ibra = 0; ibra < ms.GetTwoBodyChannel(it.first[0]).GetNumberKets(); ibra++) bk.push_back({it.first[0], it.first[1], ibra});
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < bk.size(); ib++)
  {
    size_t cb = bk[ib][0], ck = bk[ib][1], ibra = bk[ib][2];
    TwoBodyChannel &tb = ms.GetTwoBodyChannel(cb); TwoBodyChannel &tk = ms.GetTwoBodyChannel(ck);
    if (tb.J != tk.J) continue;
    int J = tb.J;
    Ket &bra = tb.GetKet(ibra);
    int i = bra.p, j = bra.q;
    Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
    if (2 * oi.n + oi.l > ms.GetEMax3Body() or 2 * oj.n + oj.l > ms.GetEMax3Body()) continue;   // zero outside the 3B space
    int phase_ij = ms.phase((oi.j2 + oj.j2) / 2 - J);
    for (size_t iket = (cb == ck ? ibra : 0); iket < tk.GetNumberKets(); iket++)
    {
      Ket &ket = tk.GetKet(iket);
      int k = ket.p, l = ket.q;
      Orbit &ok = ms.GetOrbit(k); Orbit &ol = ms.GetOrbit(l);
      if (2 * ok.n + ok.l > ms.GetEMax3Body() or 2 * ol.n + ol.l > ms.GetEMax3Body()) continue;   // zero outside the 3B space
      int phase_kl = ms.phase((ok.j2 + ol.j2) / 2 - J);
      double zijkl = 0;
      for (int J_ph = std::max(std::abs(oi.j2 - ol.j2), std::abs(ok.j2 - oj.j2)) / 2; J_ph <= std::min(oi.j2 + ol.j2, ok.j2 + oj.j2) / 2; J_ph++)
        zijkl += (2 * J_ph + 1) * six(oi.j2, oj.j2, J, ok.j2, ol.j2, J_ph) * (zbar(J_ph, i, l, k, j) + phase_ij * phase_kl * zbar(J_ph, j, k, l, i));
      for (int J_ph = std::max(std::abs(oj.j2 - ol.j2), std::abs(ok.j2 - oi.j2)) / 2; J_ph <= std::min(oj.j2 + ol.j2, ok.j2 + oi.j2) / 2; J_ph++)
        zijkl -= (2 * J_ph + 1) * six(oj.j2, oi.j2, J, ok.j2, ol.j2, J_ph) * (phase_ij * zbar(J_ph, j, l, k, i) + phase_kl * zbar(J_ph, i, k, l, j));
      zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
      Z2.AddToTBME(cb, ck, ibra, iket, zijkl);
    }
  }
  Z.profiler.timer["_comm332_pphh_plantz_back"] += omp_get_wtime() - t;
  perfplan::release_memory();
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
