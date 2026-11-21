//=====================================================================================================================
//  comm233_phss, scalar X and isospin-changing Y. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
// comm233_phss for a scalar X and a channel-changing Y (J rank 0, parity 0, T rank > 0), precompiled bookkeeping.
// Same maths as comm233_phss_tz:
//  - rows: CC channel (J, pi, |d|) kets in both orderings, split by the signed d = (tz_p - tz_q)/2;
//  - an item (row channel, d, Delta) pairs with the column channel (J, pi, d' = d - Delta); columns (lm Jlm)(ij Jij), both
//    orderings, with Tz_lm - Tz_ij = d';
//      Zbar = X2ph(rows, mids(J,pi,d)) Y3ph(mids, cols) - Y2ph(rows, mids(J,pi,d')) X3ph(mids, cols)
//  - 3B Pandya from the unpacked W (X: diagonal blocks, Y: off-diagonal blocks), 2B Pandya from the channel matrices;
//  - back-transformation into W_Z(cb, ck) and Z(cb, ck) = A(cb)^T W_Z A(ck), as in comm233_phss_plan.
//  E3max < 3 emax_3body ("blocked", as in comm233_phss_plan): Zbar(row (p q), col (glm, gij)) is only read when P(gij, p)
//  and P(glm, q) exist, i.e. e(gij) <= E3max - e_p and e(glm) <= E3max - e_q. Rows are grouped in classes of equal
//  (e_p, e_q). Per column channel the columns of each glm are sorted by e(gij); per class the stored rows are, for each
//  glm with e(glm) <= E3max - e_q, a zero row followed by the prefix of its columns with e(gij) <= E3max - e_p. Zbar is
//  computed in column chunks (which also bounds the 3B Pandya matrices) and scattered into these blocks.
// Included after PerfPlan233.inc, PerfPlan223.inc and PerfPlan333tz.inc (inside namespace Commutator).
namespace perfplan
{
  struct Col233tz
  {
    int J, par, d;
    std::vector<std::array<int, 2>> cols;   // (glm, gij)
    std::vector<long> base;                 // [class of glm] -> offset into tab, or -1
    std::vector<int> nIof;                  // [class of glm] -> size of the gij class
    std::vector<int> tab;                   // column index, or cols.size() (zero sentinel)
    // blocked storage (see the header): glm positions, 1-based position of each column in the e(gij)-sorted segment of
    // its glm (segpos[cols.size()] = 0: the zero row), and per row class the first stored row of each glm (-1: none)
    std::vector<int> glmPos, glms, segpos;
    std::vector<std::vector<int>> off;
    std::vector<size_t> nrows;
    // per column: pair energies, pair J and pair classes of (glm, gij) (all < 256); ord = columns sorted by (e(gij), e(glm))
    std::vector<unsigned char> cEij, cElm, cJij, cJlm, cCij, cClm;
    std::vector<int> ord;
  };
  struct Plan233tz
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1, TY = -1, pY = -1;
    std::vector<int> cls, clsLoc, clsSize;  // pair class (Tz+1)*2 + parity, index within the class, class sizes
    std::vector<std::array<std::vector<int>, 3>> rowsD;               // ch_cc -> [d+1] -> CC local indices r (0..2n-1)
    std::vector<std::vector<int>> rowpos;                             // ch_cc -> [r] -> position in its d list
    std::vector<std::array<std::vector<std::array<int, 2>>, 3>> midsD;   // ch_cc -> [d+1] -> ordered (a, b), n_a != n_b
    std::vector<Col233tz> colch;
    struct Item { int ch_row, d, ic; };
    std::vector<Item> items;
    std::vector<std::array<long, 9>> item_of;   // ch_cc -> [(d+1)*3 + idelta]
    std::vector<std::vector<int>> items_of_col;
    // blocked storage: row class (e_p (E3b+1) + e_q) and position within the class, per (ch_cc, d) row; class sizes
    bool blocked = false; int E3b = 0, ncls = 0;
    std::vector<std::array<std::vector<int>, 3>> rowCls, rowKpos, clsCount;
  };

  static Plan233tz &get233tz(const Operator &Z, Plan &P, int TY, int pY)
  {
    static Plan233tz Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ and Q.TY == TY and Q.pY == pY) return Q;
    Q = Plan233tz(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max(); Q.TY = TY; Q.pY = pY;
    int ng = P.g2.size();
    // classes only for pairs of the 3B space: the (g_lm, g_ij) columns are pairs of 3B states
    Q.cls.assign(ng, -1); Q.clsLoc.assign(ng, -1); Q.clsSize.assign(6, 0);
    for (int g = 0; g < P.ng3; g++)
    {
      Orbit &op = ms.GetOrbit(P.g2[g][0]); Orbit &oq = ms.GetOrbit(P.g2[g][1]);
      int c = ((op.tz2 + oq.tz2) / 2 + 1) * 2 + (op.l + oq.l) % 2;
      Q.cls[g] = c; Q.clsLoc[g] = Q.clsSize[c]++;
    }
    std::vector<std::vector<int>> byCls(6);
    for (int g = 0; g < P.ng3; g++) byCls[Q.cls[g]].push_back(g);
    size_t nch2_CC = ms.TwoBodyChannels_CC.size();
    Q.rowsD.resize(nch2_CC); Q.rowpos.resize(nch2_CC); Q.midsD.resize(nch2_CC); Q.item_of.resize(nch2_CC);
    for (size_t ch = 0; ch < nch2_CC; ch++)
    {
      auto &tbc = ms.GetTwoBodyChannel_CC(ch);
      size_t n = tbc.GetNumberKets();
      Q.rowpos[ch].assign(2 * n, -1);
      for (size_t r = 0; r < 2 * n; r++)
      {
        Ket &ket = tbc.GetKet(r % n);
        int p = (r < n) ? ket.p : ket.q, q = (r < n) ? ket.q : ket.p;
        Orbit &op = ms.GetOrbit(p); Orbit &oq = ms.GetOrbit(q);
        // rows (I3 I6) are read only for orbits of the 3B space, and the summed pairs (a b) need P(., a), P(., b)
        if (2 * op.n + op.l > ms.GetEMax3Body() or 2 * oq.n + oq.l > ms.GetEMax3Body()) continue;
        int d = (op.tz2 - oq.tz2) / 2;
        Q.rowpos[ch][r] = Q.rowsD[ch][d + 1].size();
        Q.rowsD[ch][d + 1].push_back(r);
        if (std::abs(op.occ - oq.occ) > 1e-7) Q.midsD[ch][d + 1].push_back({p, q});
      }
    }
    std::map<std::array<int, 3>, int> colIndex;
    for (size_t ch = 0; ch < nch2_CC; ch++)
    {
      Q.item_of[ch].fill(-1);
      auto &tbc = ms.GetTwoBodyChannel_CC(ch);
      for (int d = -1; d <= 1; d++)
      {
        if (Q.rowsD[ch][d + 1].empty()) continue;
        for (int idelta = 0; idelta < (TY == 0 ? 1 : 2); idelta++)
        {
          int Delta = idelta == 0 ? TY : -TY;
          std::array<int, 3> key = {tbc.J, (tbc.parity + pY) % 2, d - Delta};
          auto it = colIndex.find(key);
          int ic;
          if (it == colIndex.end()) { ic = Q.colch.size(); colIndex[key] = ic; Q.colch.push_back(Col233tz()); auto &C = Q.colch.back(); C.J = key[0]; C.par = key[1]; C.d = key[2]; }
          else ic = it->second;
          Q.item_of[ch][(d + 1) * 3 + idelta] = Q.items.size();
          Q.items.push_back({(int)ch, d, ic});
        }
      }
    }
    Q.items_of_col.resize(Q.colch.size());
    for (size_t it = 0; it < Q.items.size(); it++) Q.items_of_col[Q.items[it].ic].push_back(it);
    for (auto &C : Q.colch)
    {
      C.base.assign(6, -1); C.nIof.assign(6, 0);
      size_t off = 0;
      for (int cL = 0; cL < 6; cL++)
      {
        int TzL = cL / 2 - 1, parL = cL % 2, TzI = TzL - C.d, parI = (C.par + 2 - parL) % 2;
        if (TzI < -1 or TzI > 1) continue;
        int cI = (TzI + 1) * 2 + parI;
        C.base[cL] = off; C.nIof[cL] = Q.clsSize[cI];
        off += (size_t)Q.clsSize[cL] * Q.clsSize[cI];
      }
      C.tab.assign(off, -1);
      for (int cL = 0; cL < 6; cL++)
      {
        if (C.base[cL] < 0) continue;
        int TzL = cL / 2 - 1, parL = cL % 2, TzI = TzL - C.d, parI = (C.par + 2 - parL) % 2;
        int cI = (TzI + 1) * 2 + parI;
        for (int glm : byCls[cL])
          for (int gij : byCls[cI])
          {
            int Jlm = P.g2[glm][2], Jij = P.g2[gij][2];
            if (std::abs(Jlm - Jij) > C.J or Jlm + Jij < C.J) continue;
            C.tab[C.base[cL] + (size_t)Q.clsLoc[glm] * C.nIof[cL] + Q.clsLoc[gij]] = C.cols.size();
            C.cols.push_back({glm, gij});
          }
      }
      for (auto &v : C.tab) if (v < 0) v = C.cols.size();
    }
    // ---- blocked storage
    auto eo = [&](int o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l; };
    auto ep = [&](int g) { return eo(P.g2[g][0]) + eo(P.g2[g][1]); };
    for (auto &C : Q.colch)
    {
      size_t ncol = C.cols.size();
      C.cEij.resize(ncol); C.cElm.resize(ncol); C.cJij.resize(ncol); C.cJlm.resize(ncol); C.cCij.resize(ncol); C.cClm.resize(ncol);
      C.ord.resize(ncol);
      for (size_t col = 0; col < ncol; col++)
      {
        int glm = C.cols[col][0], gij = C.cols[col][1];
        C.cEij[col] = ep(gij); C.cElm[col] = ep(glm); C.cJij[col] = P.g2[gij][2]; C.cJlm[col] = P.g2[glm][2];
        C.cCij[col] = Q.cls[gij]; C.cClm[col] = Q.cls[glm]; C.ord[col] = col;
      }
      std::stable_sort(C.ord.begin(), C.ord.end(), [&](int a, int b) {
        return C.cEij[a] != C.cEij[b] ? C.cEij[a] < C.cEij[b] : C.cElm[a] < C.cElm[b]; });
    }
    const int E = ms.GetEMax3Body(), E3 = Q.E3max;
    Q.E3b = E; Q.ncls = (E + 1) * (E + 1);
    Q.rowCls.resize(nch2_CC); Q.rowKpos.resize(nch2_CC); Q.clsCount.resize(nch2_CC);
    for (size_t ch = 0; ch < nch2_CC; ch++)
    {
      auto &tbc = ms.GetTwoBodyChannel_CC(ch);
      size_t n = tbc.GetNumberKets();
      for (int d = 0; d < 3; d++)
      {
        Q.clsCount[ch][d].assign(Q.ncls, 0);
        for (int r : Q.rowsD[ch][d])
        {
          Ket &ket = tbc.GetKet(r % n);
          int pp = (r < (int)n) ? ket.p : ket.q, qq = (r < (int)n) ? ket.q : ket.p;
          int k = eo(pp) * (E + 1) + eo(qq);
          Q.rowCls[ch][d].push_back(k); Q.rowKpos[ch][d].push_back(Q.clsCount[ch][d][k]++);
        }
      }
    }
    for (auto &C : Q.colch)
    {
      size_t ncol = C.cols.size();
      C.glmPos.assign(P.ng3, -1);
      for (auto &c : C.cols) if (C.glmPos[c[0]] < 0) { C.glmPos[c[0]] = C.glms.size(); C.glms.push_back(c[0]); }
      std::vector<std::vector<int>> seg(C.glms.size());
      for (size_t col = 0; col < ncol; col++) seg[C.glmPos[C.cols[col][0]]].push_back(col);
      C.segpos.assign(ncol + 1, 0);
      std::vector<std::vector<int>> segE(C.glms.size());   // e(gij) along each sorted segment
      for (size_t gp = 0; gp < seg.size(); gp++)
      {
        std::stable_sort(seg[gp].begin(), seg[gp].end(), [&](int a, int b) { return ep(C.cols[a][1]) < ep(C.cols[b][1]); });
        for (size_t t = 0; t < seg[gp].size(); t++) { C.segpos[seg[gp][t]] = t + 1; segE[gp].push_back(ep(C.cols[seg[gp][t]][1])); }
      }
      C.off.assign(Q.ncls, {}); C.nrows.assign(Q.ncls, 0);
      for (int k = 0; k < Q.ncls; k++)
      {
        int ec = k / (E + 1), ed = k % (E + 1);
        C.off[k].assign(C.glms.size(), -1);
        size_t nr = 0;
        for (size_t gp = 0; gp < C.glms.size(); gp++)
        {
          if (ep(C.glms[gp]) > E3 - ed) continue;
          C.off[k][gp] = nr;
          nr += 1 + (std::upper_bound(segE[gp].begin(), segE[gp].end(), E3 - ec) - segE[gp].begin());
        }
        C.nrows[k] = nr;
      }
    }
    {
      double used = 0, full = 0;
      for (auto &it : Q.items)
      {
        auto &C = Q.colch[it.ic];
        full += (double)(C.cols.size() + 1) * Q.rowsD[it.ch_row][it.d + 1].size();
        for (int k = 0; k < Q.ncls; k++) used += (double)C.nrows[k] * Q.clsCount[it.ch_row][it.d + 1][k];
      }
      Q.blocked = used < 0.9 * full;
      if (std::getenv("IMSRG_PLAN_DEBUG") != nullptr)
        std::cout << "PLANDEBUG comm233_ph tz: blocked " << Q.blocked << ", stored " << used << " of " << full << " entries" << std::endl;
      if (std::getenv("IMSRG_PLAN_DEBUG") != nullptr)
      {
        double gc = 0, gt = 0, gb = 0;   // GB: columns + per-column data, column tables, blocked offsets
        for (auto &C : Q.colch)
        {
          gc += (8.0 + 6 + 4 + 4) * 1e-9 * C.cols.size();   // cols, 6 bytes, ord, segpos
          gt += 4e-9 * C.tab.size() + 8e-9 * (C.base.size() + C.nIof.size());
          for (auto &o : C.off) gb += 4e-9 * o.size();
          gb += 4e-9 * (C.glmPos.size() + C.glms.size());
        }
        std::cout << "PLANDEBUG comm233_ph tz plan memory (GB): columns " << gc << ", column tables " << gt << ", blocked offsets " << gb << std::endl;
      }
    }
    return Q;
  }
} // namespace perfplan

void comm233_phss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced())
  {
    comm233_phss_tz(X, Y, Z);
    return;
  }
  auto &Z3 = Z.ThreeBody;
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  int TY = Y.GetTRank(), hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(Z, Z);
  const perfplan::Cut3 &cutZ = perfplan::cut3(*Z.modelspace);
  perfplan::Plan233 &Q = perfplan::get233(Z, P);
  perfplan::Plan223 &Q2 = perfplan::get223(Z, P);
  perfplan::Plan233tz &T = perfplan::get233tz(Z, P, TY, Y.GetParity());
  int norb = P.norb, nJ = P.J2max + 1;
  size_t nch2 = ms.GetNumberTwoBodyChannels(), nch2_CC = ms.TwoBodyChannels_CC.size();
  double t = omp_get_wtime();
  static const std::vector<arma::mat> empty;
  static const std::vector<perfplan::wmat> emptyW;
  const std::vector<perfplan::wmat> &WX = x3 ? perfplan::unpack(X, P) : emptyW;
  static perfplan::WOffCache emptyoff;
  const perfplan::WOffCache &WY = y3 ? perfplan::unpack_off(Y, P) : emptyoff;
  Z.profiler.timer["_comm233_ph_plantz_unpack"] += omp_get_wtime() - t; t = omp_get_wtime();

  // 2B matrices: X diagonal, Y both orientations
  std::vector<const arma::mat *> XM(nch2);
  for (size_t ch = 0; ch < nch2; ch++) XM[ch] = &X.TwoBody.GetMatrix(ch, ch);
  std::vector<arma::mat> YF(nch2 * nch2);
  for (auto &kv : Y.TwoBody.MatEl)
  {
    size_t cb = kv.first[0], ck = kv.first[1];
    YF[cb * nch2 + ck] = kv.second;
    if (cb != ck) YF[ck * nch2 + cb] = hY * kv.second.t();
  }
  // un-normalised TBME <pq; J| O |rs; J> (= GetTBME_J), isY selects the operator
  auto tbme = [&](bool isY, int J, int p, int q, int r, int s) -> double
  {
    int gb = P.g2idx[((size_t)std::min(p, q) * norb + std::max(p, q)) * (P.J2max + 1) + J];
    int gk = P.g2idx[((size_t)std::min(r, s) * norb + std::max(r, s)) * (P.J2max + 1) + J];
    if (gb < 0 or gk < 0) return 0.0;
    int cb = Q2.chOf[gb], ck = Q2.chOf[gk];
    double v;
    if (isY) { const arma::mat &M = YF[(size_t)cb * nch2 + ck]; if (M.n_elem == 0) return 0.0; v = M(Q2.locOf[gb], Q2.locOf[gk]); }
    else { if (cb != ck) return 0.0; v = (*XM[cb])(Q2.locOf[gb], Q2.locOf[gk]); }
    if (p > q) v *= Q2.kph[gb];
    if (r > s) v *= Q2.kph[gk];
    if (p == q) v *= PhysConst::SQRT2;
    if (r == s) v *= PhysConst::SQRT2;
    return v;
  };
  // 2B Pandya: -sum_J' (2J'+1) {jp jq J; ja jb J'} O_{pb,aq}^{J'} (n_a - n_b)
  auto pandya2 = [&](bool isY, int Jph, int p, int q, int a, int b) -> double
  {
    Orbit &op = ms.GetOrbit(p); Orbit &oq = ms.GetOrbit(q); Orbit &oa = ms.GetOrbit(a); Orbit &ob = ms.GetOrbit(b);
    int Jmin = std::max(std::abs(op.j2 - ob.j2), std::abs(oa.j2 - oq.j2)) / 2, Jmax = std::min(op.j2 + ob.j2, oa.j2 + oq.j2) / 2;
    double sum = 0;
    for (int Jp = Jmin; Jp <= Jmax; Jp++)
      sum -= (2 * Jp + 1) * ms.GetSixJ(0.5 * op.j2, 0.5 * oq.j2, Jph, 0.5 * oa.j2, 0.5 * ob.j2, Jp) * tbme(isY, Jp, p, b, a, q);
    return sum * (oa.occ - ob.occ);
  };
  // 3B Pandya: -sum_J' (-1)^((ja+J')) (2J'+1) {Jij Jlm Jph; jb ja J'} <(ij) a; J'| O |(lm) b; J'>
  // Columns cl[0 .. ncl). P(gij, a) and P(glm, b) only exist for e(gij) + e_a <= E3max and e(glm) + e_b <= E3max, so other
  // columns are skipped (exact zeros). Per (a b), the 3B channels, the W blocks and the 6j weights depend only on the pair
  // classes, pair J values and J' of a column, so they are small tables.
  const int nT3 = (P.tJmax3 + 1) / 2, nJ3 = P.J2max3 + 1;
  auto fill3 = [&](arma::mat &O3ph, bool isY, const std::vector<std::array<int, 2>> &mids, const perfplan::Col233tz &C,
                   const int *cl, size_t ncl)
  {
    size_t nm = mids.size();
    O3ph.zeros(nm, ncl);
    std::vector<int> chbT(6 * nT3), chkT(6 * nT3);
    std::vector<const perfplan::wmat *> Wt(36 * nT3);
    std::vector<double> w6((size_t)nJ3 * nJ3 * nT3);
    for (size_t im = 0; im < nm; im++)
    {
      int a = mids[im][0], b = mids[im][1];
      Orbit &oa = ms.GetOrbit(a); Orbit &ob = ms.GetOrbit(b);
      if (std::abs(oa.occ - ob.occ) < 1e-6) continue;
      const int lim_ij = T.E3max - (2 * oa.n + oa.l), lim_lm = T.E3max - (2 * ob.n + ob.l);
      for (int c = 0; c < 6; c++)   // pair class c = (Tz + 1) * 2 + parity
        for (int t = 0; t < nT3; t++)
        {
          int tzp = 2 * (c / 2 - 1), parp = c % 2;
          chbT[c * nT3 + t] = P.ch3tab[P.i3(2 * t + 1, (parp + oa.l) % 2, tzp + oa.tz2)];
          chkT[c * nT3 + t] = P.ch3tab[P.i3(2 * t + 1, (parp + ob.l) % 2, tzp + ob.tz2)];
        }
      for (int ci = 0; ci < 6; ci++)
        for (int ck = 0; ck < 6; ck++)
          for (int t = 0; t < nT3; t++)
          {
            int chb = chbT[ci * nT3 + t], chk = chkT[ck * nT3 + t];
            const perfplan::wmat *W = nullptr;
            if (chb >= 0 and chk >= 0)
            {
              if (isY) { int k = WY.idx[(size_t)chb * WY.nch + chk]; if (k >= 0) W = &WY.mats[k]; }
              else if (chb == chk) W = &WX[chb];
              if (W != nullptr and W->n_elem == 0) W = nullptr;
            }
            Wt[(ci * 6 + ck) * nT3 + t] = W;
          }
      std::fill(w6.begin(), w6.end(), NAN);   // NaN = not computed yet
      for (size_t j = 0; j < ncl; j++)
      {
        int col = cl != nullptr ? cl[j] : (int)j;   // nullptr: columns 0 .. ncl
        if (C.cEij[col] > lim_ij or C.cElm[col] > lim_lm) continue;
        int glm = C.cols[col][0], gij = C.cols[col][1];
        int Jlm = C.cJlm[col], Jij = C.cJij[col], ci = C.cCij[col], ck = C.cClm[col];
        int lo = std::max(std::abs(oa.j2 - 2 * Jij), std::abs(ob.j2 - 2 * Jlm)), hi = std::min(oa.j2 + 2 * Jij, ob.j2 + 2 * Jlm);
        double s = 0;
        for (int tJp = lo; tJp <= hi; tJp += 2)
        {
          int t = tJp / 2;
          const perfplan::wmat *W = Wt[(ci * 6 + ck) * nT3 + t];
          if (W == nullptr) continue;
          int p1 = P.pidx(chbT[ci * nT3 + t], gij, a), p2 = P.pidx(chkT[ck * nT3 + t], glm, b);
          if (p1 < 0 or p2 < 0) continue;
          double &w = w6[((size_t)Jij * nJ3 + Jlm) * nT3 + t];
          if (std::isnan(w)) w = ms.phase((oa.j2 + tJp) / 2) * (tJp + 1) * Q.six6(2 * Jij, 2 * Jlm, 2 * C.J, ob.j2, oa.j2, tJp);
          s -= w * (*W)(p1, p2);
        }
        O3ph(im, j) = s;
      }
    }
  };
  auto midsOf = [&](int J, int par, int d) -> const std::vector<std::array<int, 2>> *
  {
    if (std::abs(d) > 1) return nullptr;
    int ch = Q.chph[(size_t)J * 4 + par * 2 + std::abs(d)];
    if (ch < 0) return nullptr;
    return &T.midsD[ch][d + 1];
  };

  // ---- particle-hole Zbar per item, stored transposed (cols + sentinel row) x rows
  std::vector<arma::mat> ZbT(T.items.size());
  std::vector<std::vector<arma::mat>> Zst(T.blocked ? T.items.size() : 0);   // blocked: [item][row class]
  size_t ncolch = T.colch.size();
  double tt_p2 = 0, tt_f3 = 0, tt_mm = 0, tt_sc = 0;   // thread time (blocked branch): 2B Pandya, 3B fill, GEMM, scatter
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : tt_p2, tt_f3, tt_mm, tt_sc)
  for (size_t ic = 0; ic < ncolch; ic++)
  {
    auto &C = T.colch[ic];
    size_t ncol = C.cols.size();
    const std::vector<std::array<int, 2>> *mid_c = midsOf(C.J, C.par, C.d);
    if (T.blocked)
    {
      const int E = T.E3b, E3 = T.E3max;
      auto &its = T.items_of_col[ic];
      std::vector<arma::mat> X2s(its.size()), Y2s(its.size());
      // per item: its row classes, the rows of each class (in row order, so the position in the list is the column of
      // Zst), and the smallest e_p, e_q over its rows (columns beyond E3max - those are never read)
      std::vector<std::vector<int>> ucls(its.size());
      std::vector<std::vector<std::vector<int>>> urows(its.size());
      std::vector<int> minEc(its.size(), 1 << 20), minEd(its.size(), 1 << 20);
      std::vector<std::vector<arma::mat>> X2k(its.size()), Y2k(its.size());   // [item][class]: rows of that class
      double tq = omp_get_wtime();
      for (size_t u = 0; u < its.size(); u++)   // 2B Pandya matrices (rows x mids), once per item
      {
        auto &item = T.items[its[u]];
        auto &tbc = ms.GetTwoBodyChannel_CC(item.ch_row);
        size_t n = tbc.GetNumberKets();
        auto &rows = T.rowsD[item.ch_row][item.d + 1];
        auto &mid_r = T.midsD[item.ch_row][item.d + 1];
        auto rowpq = [&](size_t irow, int &p, int &q) { int r = rows[irow]; Ket &k = tbc.GetKet(r % n); p = (r < (int)n) ? k.p : k.q; q = (r < (int)n) ? k.q : k.p; };
        if (y3 and mid_r.size() > 0)
        {
          X2s[u].set_size(rows.size(), mid_r.size());
          for (size_t irow = 0; irow < rows.size(); irow++)
          { int p, q; rowpq(irow, p, q); for (size_t im = 0; im < mid_r.size(); im++) X2s[u](irow, im) = pandya2(false, tbc.J, p, q, mid_r[im][0], mid_r[im][1]); }
        }
        if (x3 and mid_c != nullptr and mid_c->size() > 0)
        {
          Y2s[u].set_size(rows.size(), mid_c->size());
          for (size_t irow = 0; irow < rows.size(); irow++)
          { int p, q; rowpq(irow, p, q); for (size_t im = 0; im < mid_c->size(); im++) Y2s[u](irow, im) = pandya2(true, tbc.J, p, q, (*mid_c)[im][0], (*mid_c)[im][1]); }
        }
        Zst[its[u]].resize(T.ncls);
        urows[u].resize(T.ncls);
        auto &rc = T.rowCls[item.ch_row][item.d + 1];
        for (size_t irow = 0; irow < rows.size(); irow++) urows[u][rc[irow]].push_back(irow);
        for (int k = 0; k < T.ncls; k++)
        {
          if (urows[u][k].empty()) continue;
          Zst[its[u]][k].zeros(C.nrows[k], urows[u][k].size());
          ucls[u].push_back(k);
          minEc[u] = std::min(minEc[u], k / (E + 1)); minEd[u] = std::min(minEd[u], k % (E + 1));
        }
        X2k[u].resize(T.ncls); Y2k[u].resize(T.ncls);
        for (int k : ucls[u])
        {
          arma::uvec rk(urows[u][k].size());
          for (size_t t = 0; t < rk.n_elem; t++) rk[t] = urows[u][k][t];
          if (X2s[u].n_elem > 0) X2k[u][k] = X2s[u].rows(rk);
          if (Y2s[u].n_elem > 0) Y2k[u][k] = Y2s[u].rows(rk);
        }
      }
      tt_p2 += omp_get_wtime() - tq;
      const size_t CH = 4096;   // column chunk, in the order C.ord (sorted by e(gij), then e(glm))
      for (size_t c0 = 0; c0 < ncol; c0 += CH)
      {
        size_t c1 = std::min(ncol, c0 + CH), nch = c1 - c0;
        const int *cl = &C.ord[c0];
        int minij = C.cEij[cl[0]], minlm = 1 << 20;   // sorted: the first column has the smallest e(gij)
        for (size_t j = 0; j < nch; j++) minlm = std::min(minlm, (int)C.cElm[cl[j]]);
        arma::mat X3c; bool haveX3 = false;
        for (size_t u = 0; u < its.size(); u++)
        {
          auto &item = T.items[its[u]];
          auto &mid_r = T.midsD[item.ch_row][item.d + 1];
          if (ucls[u].empty() or minij > E3 - minEc[u] or minlm > E3 - minEd[u]) continue;   // nothing in this chunk is read
          tq = omp_get_wtime();
          if (not haveX3 and x3 and mid_c != nullptr and mid_c->size() > 0) fill3(X3c, false, *mid_c, C, cl, nch);
          haveX3 = true;
          arma::mat Y3c;
          if (X2s[u].n_elem > 0) fill3(Y3c, true, mid_r, C, cl, nch);
          tt_f3 += omp_get_wtime() - tq;
          // per row class: the chunk columns it reads (e(gij) <= E3max - e_p, e(glm) <= E3max - e_q), then
          // Zk = X2(rows of k) Y3(:, cols) - Y2(rows of k) X3(:, cols), written into its block
          std::vector<arma::uword> js;
          for (int k : ucls[u])
          {
            const int lij = E3 - k / (E + 1), llm = E3 - k % (E + 1);
            js.clear();
            for (size_t j = 0; j < nch; j++) if (C.cEij[cl[j]] <= lij and C.cElm[cl[j]] <= llm) js.push_back(j);
            if (js.empty()) continue;
            double tq2 = omp_get_wtime();
            arma::uvec jv(js);
            const size_t nr = urows[u][k].size();
            arma::mat Zk(nr, js.size(), arma::fill::zeros);
            if (X2s[u].n_elem > 0) Zk += X2k[u][k] * Y3c.cols(jv);
            if (Y2s[u].n_elem > 0 and X3c.n_elem > 0) Zk -= Y2k[u][k] * X3c.cols(jv);
            double tq3 = omp_get_wtime(); tt_mm += tq3 - tq2;
            arma::mat &M = Zst[its[u]][k];
            const size_t ld = M.n_rows;
            for (size_t t = 0; t < js.size(); t++)
            {
              int col = cl[js[t]];
              double *dst = M.memptr() + C.off[k][C.glmPos[C.cols[col][0]]] + C.segpos[col];
              const double *zc = Zk.colptr(t);
              for (size_t r = 0; r < nr; r++) dst[r * ld] = zc[r];
            }
            tt_sc += omp_get_wtime() - tq3;
          }
        }
      }
      continue;
    }
    arma::mat X3ph;
    if (x3 and mid_c != nullptr and mid_c->size() > 0 and ncol > 0) fill3(X3ph, false, *mid_c, C, nullptr, ncol);
    for (int it : T.items_of_col[ic])
    {
      auto &item = T.items[it];
      auto &tbc = ms.GetTwoBodyChannel_CC(item.ch_row);
      size_t n = tbc.GetNumberKets();
      auto &rows = T.rowsD[item.ch_row][item.d + 1];
      auto &mid_r = T.midsD[item.ch_row][item.d + 1];
      size_t nrow = rows.size();
      arma::mat Zb(nrow, ncol, arma::fill::zeros);
      auto rowpq = [&](size_t irow, int &p, int &q) { int r = rows[irow]; Ket &k = tbc.GetKet(r % n); p = (r < (int)n) ? k.p : k.q; q = (r < (int)n) ? k.q : k.p; };
      if (y3 and mid_r.size() > 0 and ncol > 0)
      {
        arma::mat X2ph(nrow, mid_r.size());
        for (size_t irow = 0; irow < nrow; irow++)
        {
          int p, q; rowpq(irow, p, q);
          for (size_t im = 0; im < mid_r.size(); im++) X2ph(irow, im) = pandya2(false, tbc.J, p, q, mid_r[im][0], mid_r[im][1]);
        }
        arma::mat Y3ph; fill3(Y3ph, true, mid_r, C, nullptr, ncol);
        Zb += X2ph * Y3ph;
      }
      if (X3ph.n_elem > 0)
      {
        arma::mat Y2ph(nrow, mid_c->size());
        for (size_t irow = 0; irow < nrow; irow++)
        {
          int p, q; rowpq(irow, p, q);
          for (size_t im = 0; im < mid_c->size(); im++) Y2ph(irow, im) = pandya2(true, tbc.J, p, q, (*mid_c)[im][0], (*mid_c)[im][1]);
        }
        Zb -= Y2ph * X3ph;
      }
      ZbT[it] = Zb.t();
      ZbT[it].resize(ncol + 1, nrow);   // zero sentinel row
    }
  }
  Z.profiler.timer["_comm233_ph_plantz_fill_gemm"] += omp_get_wtime() - t; t = omp_get_wtime();
  if (T.blocked)
  {
    Z.profiler.timer["_comm233_ph_plantz_tt_pandya2"] += tt_p2; Z.profiler.timer["_comm233_ph_plantz_tt_fill3"] += tt_f3;
    Z.profiler.timer["_comm233_ph_plantz_tt_gemm"] += tt_mm; Z.profiler.timer["_comm233_ph_plantz_tt_scatter"] += tt_sc;
  }

  // ---- back-transformation into the stored blocks (cb, ck), cb != ck
  std::vector<std::array<size_t, 2>> blocks;
  for (auto &it : Z3.Get_ch_start()) blocks.push_back({it.first.ch_bra, it.first.ch_ket});
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < blocks.size(); ib++)
  {
    size_t cb = blocks[ib][0], ck = blocks[ib][1];
    if (cb == ck) continue;
    auto &Tb = ms.GetThreeBodyChannel(cb); auto &Tk = ms.GetThreeBodyChannel(ck);
    if (Tb.twoJ != Tk.twoJ) continue;
    int twoJ = Tb.twoJ;
    size_t nPb = P.Pst[cb].size(), nPk = P.Pst[ck].size(), nb = Tb.GetNumberKets(), nk = Tk.GetNumberKets();
    if (nPb == 0 or nPk == 0 or nb == 0 or nk == 0) continue;
    arma::mat WZ(nPb, nPk, arma::fill::zeros);
    std::vector<std::vector<int>> rowsBy(norb), colsBy(norb);
    for (size_t p = 0; p < nPb; p++) rowsBy[P.Pst[cb][p][1]].push_back(p);
    for (size_t p = 0; p < nPk; p++) colsBy[P.Pst[ck][p][1]].push_back(p);
    std::vector<double> wt(nJ * nJ);
    for (int I3 = 0; I3 < norb; I3++)
    {
      auto &Rw = rowsBy[I3];
      if (Rw.empty()) continue;
      Orbit &o3 = ms.GetOrbit(I3);
      size_t nr = Rw.size();
      std::vector<int> rg(nr), rJ(nr), rloc(nr); std::vector<char> hasR(nJ, 0);
      for (size_t r = 0; r < nr; r++) { rg[r] = P.Pst[cb][Rw[r]][0]; rJ[r] = P.g2[rg[r]][2]; rloc[r] = T.clsLoc[rg[r]]; hasR[rJ[r]] = 1; }
      Orbit &o1 = ms.GetOrbit(P.g2[rg[0]][0]); Orbit &o2 = ms.GetOrbit(P.g2[rg[0]][1]);
      int Tz12 = (o1.tz2 + o2.tz2) / 2;
      double phase_ph = -ms.phase((twoJ + o3.j2) / 2);
      for (int I6 = 0; I6 < norb; I6++)
      {
        auto &Cl = colsBy[I6];
        if (Cl.empty()) continue;
        Orbit &o6 = ms.GetOrbit(I6);
        size_t nc = Cl.size();
        int h0 = P.Pst[ck][Cl[0]][0];
        Orbit &o4 = ms.GetOrbit(P.g2[h0][0]); Orbit &o5 = ms.GetOrbit(P.g2[h0][1]);
        int Tz45 = (o4.tz2 + o5.tz2) / 2;
        int d36 = (o3.tz2 - o6.tz2) / 2, dp = Tz45 - Tz12, Delta = d36 - dp;
        int idelta = (Delta == TY) ? 0 : ((Delta == -TY) ? 1 : -1);
        if (idelta < 0) continue;
        int parity_ph = (o3.l + o6.l) % 2;
        int cL = T.cls[h0];
        for (int Jph = std::abs(o3.j2 - o6.j2) / 2; Jph <= (o3.j2 + o6.j2) / 2; Jph++)
        {
          int ch_ph = Q.chph[(size_t)Jph * 4 + parity_ph * 2 + std::abs(d36)];
          if (ch_ph < 0) continue;
          long itm = T.item_of[ch_ph][(d36 + 1) * 3 + idelta];
          if (itm < 0) continue;
          auto &C = T.colch[T.items[itm].ic];
          if (C.base[cL] < 0) continue;
          int k = Q.ccLocal[ch_ph][(size_t)I3 * norb + I6];
          if (k < 0) continue;
          int irow = T.rowpos[ch_ph][k];
          double pref = -phase_ph * (2 * Jph + 1);
          for (int J12 = 0; J12 < nJ; J12++)
            if (hasR[J12])
              for (int J45 = 0; J45 < nJ; J45++) wt[J12 * nJ + J45] = pref * Q.six6(2 * J12, 2 * J45, 2 * Jph, o6.j2, o3.j2, twoJ);
          const int *tabL = &C.tab[C.base[cL]];
          int nI = C.nIof[cL];
          if (T.blocked)
          {
            // P(g12, I3) and P(g45, I6) exist, so every read below lies in the stored prefix of its glm (see the header)
            int kc = T.rowCls[ch_ph][d36 + 1][irow];
            const double *zr = Zst[itm][kc].colptr(T.rowKpos[ch_ph][d36 + 1][irow]);
            const std::vector<int> &offk = C.off[kc];
            for (size_t q = 0; q < nc; q++)
            {
              int g45 = P.Pst[ck][Cl[q]][0], J45 = P.g2[g45][2];
              int gp = C.glmPos[g45];
              if (gp < 0 or offk[gp] < 0) continue;   // no column with this glm: all entries are zero
              const double *zg = zr + offk[gp];
              const int *trow = tabL + (size_t)T.clsLoc[g45] * nI;
              double *wc = WZ.colptr(Cl[q]);
              for (size_t r = 0; r < nr; r++) wc[Rw[r]] += wt[rJ[r] * nJ + J45] * zg[C.segpos[trow[rloc[r]]]];
            }
            continue;
          }
          const double *zr = ZbT[itm].colptr(irow);
          for (size_t q = 0; q < nc; q++)
          {
            int g45 = P.Pst[ck][Cl[q]][0], J45 = P.g2[g45][2];
            const int *trow = tabL + (size_t)T.clsLoc[g45] * nI;
            double *wc = WZ.colptr(Cl[q]);
            for (size_t r = 0; r < nr; r++) wc[Rw[r]] += wt[rJ[r] * nJ + J45] * zr[trow[rloc[r]]];
          }
        }
      }
    }
    arma::mat Zk = perfplan::AtWA(P.A[cb], WZ, P.A[ck]);
    for (size_t ibra = 0; ibra < nb; ibra++)
      for (size_t iket = 0; iket < nk; iket++)
        if (perfplan::pass3(cutZ, cb, ibra) and perfplan::pass3(cutZ, ck, iket)) Z3.AddToME_pn_ch(cb, ck, ibra, iket, Zk(ibra, iket));
  }
  Z.profiler.timer["_comm233_ph_plantz_back"] += omp_get_wtime() - t;
  perfplan::release_memory();
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
