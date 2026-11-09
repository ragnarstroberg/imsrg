//=====================================================================================================================
//  comm333_pph_hhpss, scalar X and isospin-changing Y (+ the off-diagonal unpack used by all *_tz kernels)
//=====================================================================================================================
//  Isospin-changing operators (all *_plan_tz kernels): an operator Y with J rank 0, parity 0 and isospin rank T_Y > 0
//  (e.g. the Fermi operator T+- in the transformation for delta_C) changes Tz by Delta = +-T_Y. Its 3B part therefore
//  connects different 3B channels: the stored blocks are Y(cb, ck), bra channel cb != ket channel ck, with
//  2Tz(cb) = 2Tz(ck) + 2 Delta. In the product basis each block unpacks to
//      W_Y(cb, ck) = R(cb) Y(cb, ck) R(ck)^T          (WOffCache below; both orientations (cb, ck) and (ck, cb)),
//  and a result built in the product basis is antisymmetrised per block as  Z(cb, ck) = A(cb)^T W_Z(cb, ck) A(ck).
//  X stays scalar and uses the ordinary unpack of PerfPlan.inc. Every tz kernel computes the same matrix elements as
//  the corresponding *_tz term in IMSRG3Commutators.cc (validated against ReferenceImplementations and the m-scheme
//  tests, see PLAN_KERNELS.md).
//
// comm333_pph_hhpss for a scalar X and a channel-changing Y (J rank 0, parity 0, T rank > 0: e.g. the Fermi operator
// in the transform), with precompiled bookkeeping. Same maths as comm333_pph_hhpss_tz:
//   A(r,c)    = hY Xbar(r) D(r) Ybar(r->c)^T            r, c pph channels with twoTz_c = twoTz_r - 2 Delta, Delta = +-T_Y
//   Zbar(r,c) = A(r,c) - hX hY A(c,r)^T
// with Xbar from the unpacked W_X (diagonal blocks) and Ybar from the unpacked off-diagonal blocks
// W_Y(cb, ck) = R(cb) Y(cb, ck) R(ck)^T. The back-transformation goes into each stored block (cb, ck) of Z:
//   W_Z(cb, ck) += (2Jph+1) 6j Zbar(r, c),   Z(cb, ck) = A(cb)^T W_Z A(ck).
// Weights use ModelSpace::GetSixJ (the plan's sparse 6j table only holds the scalar combinations).
// E3max < 3 emax_3body ("blocked", as in comm333_pph_hhpss_plan): the back-transformation reads Zbar(r, c)(e1, e2), with
// e1 = (g1, n1) and e2 = (g2, n2), only if e(g1) + e(n2) <= E3max and e(g2) + e(n1) <= E3max. The ext states of a channel
// are sorted by pair energy within each single orbit, so for the orbit block (n1, n2) the used part is the corner
// Cr.pre[n1 norb + n2] x Cc.pre[n2 norb + n1]. For such channel pairs only these corners of A and Zbar are computed and
// stored (A(c, r) has the transposed corners, so Zbar = A(r, c) - hX hY A(c, r)^T works corner by corner).
// Included after PerfPlan.inc (inside namespace Commutator).
namespace perfplan
{
  // Unpacked off-diagonal blocks of a channel-changing operator, both orientations, cached per Commutator call.
  // mats = R(cb) Y R(ck)^T, kk = Y(cb, ck) (stored kets), blk = (cb, ck) per block; R(cb) Y is computed on demand (unpack_off_rk)
  struct WOffCache { long epoch = -1; const void *key = nullptr; size_t nch = 0; std::vector<wmat> mats; std::vector<arma::mat> kk; std::vector<int> idx;
                     std::vector<std::array<size_t, 2>> blk; };
  static WOffCache woff;
  static const WOffCache &unpack_off(const Operator &op, Plan &P)
  {
    std::lock_guard<std::mutex> lock(wmtx);
    if (scope_active and woff.epoch == epoch and woff.key == &op.ThreeBody) return woff;   // reuse only inside a commutator scope
    ModelSpace &ms = *op.modelspace;
    size_t nch = ms.GetNumberThreeBodyChannels();
    woff = WOffCache(); woff.nch = nch; woff.idx.assign(nch * nch, -1);
    std::vector<std::array<size_t, 2>> bl;
    for (auto &it : op.ThreeBody.Get_ch_start())
    {
      size_t cb = it.first.ch_bra, ck = it.first.ch_ket;
      if (cb == ck or P.Pst[cb].empty() or P.Pst[ck].empty()) continue;
      woff.idx[cb * nch + ck] = bl.size(); bl.push_back({cb, ck});
      woff.idx[ck * nch + cb] = bl.size(); bl.push_back({ck, cb});
    }
    woff.mats.resize(bl.size()); woff.kk.resize(bl.size()); woff.blk = bl;
    auto &T3 = op.ThreeBody;
    const Cut3 &cut = cut3(ms);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < bl.size(); b++)
    {
      size_t cb = bl[b][0], ck = bl[b][1];
      size_t nb = ms.GetThreeBodyChannel(cb).GetNumberKets(), nk = ms.GetThreeBodyChannel(ck).GetNumberKets();
      arma::mat Yk(nb, nk);
      for (size_t i = 0; i < nb; i++)
        for (size_t k = 0; k < nk; k++) Yk(i, k) = T3.GetME_pn_ch(cb, ck, i, k);
      mask3(cut, cb, ck, Yk);   // dE3max / occnat3 cut (no-op without cuts)
      arma::mat rk = AtM(P.Rt[cb], Yk);          // R(cb) Y
      woff.mats[b] = to_wmat(DS(rk, P.Rt[ck]));  // (R(cb) Y) R(ck)^T
      woff.kk[b] = std::move(Yk);
    }
    woff.epoch = epoch; woff.key = &op.ThreeBody;
    return woff;
  }
  // R(cb) Y(cb, ck) per block of a WOffCache, computed on demand (held only by the terms that read it)
  static std::vector<arma::mat> unpack_off_rk(const WOffCache &w, Plan &P)
  {
    std::vector<arma::mat> rk(w.kk.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < w.kk.size(); b++) rk[b] = AtM(P.Rt[w.blk[b][0]], w.kk[b]);
    return rk;
  }
} // namespace perfplan

void comm333_pph_hhpss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced())
  {
    comm333_pph_hhpss_tz(X, Y, Z);
    return;
  }
  auto &Z3 = Z.ThreeBody;
  if (not(X.ThreeBody.IsAllocated() and Y.ThreeBody.IsAllocated())) return;
  int hX = X.IsHermitian() ? 1 : -1, hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(X, Z);
  const perfplan::Cut3 &cutZ = perfplan::cut3(*Z.modelspace);
  int norb = P.norb, nJ = P.J2max + 1;
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  double t = omp_get_wtime();
  const std::vector<perfplan::wmat> &WX = perfplan::unpack(X, P);
  const perfplan::WOffCache &WY = perfplan::unpack_off(Y, P);
  Z.profiler.timer["_comm333_pph_plantz_unpack"] += omp_get_wtime() - t; t = omp_get_wtime();
  auto six = [&](int a, int b, int c, int d, int e, int f) { return ms.GetSixJ(0.5 * a, 0.5 * b, 0.5 * c, 0.5 * d, 0.5 * e, 0.5 * f); };

  // Out(e, a) = sum_J' (2J'+1) 6j(jn Jph J12; jc J' Jab) W(P(g12, c), P(gab, n)),  e = (g12, n) in Ce.ext, a = (gab, c) in Ca.intl,
  // W = getW(channel of (g12, c), channel of (gab, n))
  auto fillbar = [&](arma::mat &Out, const perfplan::Ch &Ce, const perfplan::Ch &Ca, auto getW)
  {
    size_t ne = Ce.ext.size(), ni = Ca.intl.size();
    Out.zeros(ne, ni);
    if (ne == 0 or ni == 0) return;
    std::map<int, std::vector<int>> eByN, aByC;
    for (size_t e = 0; e < ne; e++) eByN[Ce.ext[e][1]].push_back(e);
    for (size_t a = 0; a < ni; a++) aByC[Ca.intl[a][1]].push_back(a);
    std::vector<int> er, pr, rJ, ac, pc, cJ; std::vector<double> wt(nJ * nJ);
    for (auto &en : eByN)
    {
      int n = en.first; Orbit &on = ms.GetOrbit(n);
      int g0 = Ce.ext[en.second[0]][0];
      Orbit &oi = ms.GetOrbit(P.g2[g0][0]); Orbit &oj = ms.GetOrbit(P.g2[g0][1]);
      for (auto &acn : aByC)
      {
        int cc = acn.first; Orbit &oc = ms.GetOrbit(cc);
        int h0 = Ca.intl[acn.second[0]][0];
        Orbit &oa = ms.GetOrbit(P.g2[h0][0]); Orbit &ob = ms.GetOrbit(P.g2[h0][1]);
        int par_b = (oi.l + oj.l + oc.l) % 2, tz_b = oi.tz2 + oj.tz2 + oc.tz2;
        int par_k = (oa.l + ob.l + on.l) % 2, tz_k = oa.tz2 + ob.tz2 + on.tz2;
        for (int tJp = 1; tJp <= P.tJmax; tJp += 2)
        {
          int chb = P.ch3tab[P.i3(tJp, par_b, tz_b)], chk = P.ch3tab[P.i3(tJp, par_k, tz_k)];
          if (chb < 0 or chk < 0) continue;
          const perfplan::wmat *W = getW(chb, chk);
          if (W == nullptr or W->n_elem == 0) continue;
          er.clear(); pr.clear(); rJ.clear(); ac.clear(); pc.clear(); cJ.clear();
          for (int e : en.second)
          {
            int g12 = Ce.ext[e][0], J12 = P.g2[g12][2];
            if (std::abs(2 * J12 - oc.j2) > tJp or 2 * J12 + oc.j2 < tJp) continue;
            int p1 = P.pidx(chb, g12, cc);
            if (p1 >= 0) { er.push_back(e); pr.push_back(p1); rJ.push_back(J12); }
          }
          if (er.empty()) continue;
          for (int a : acn.second)
          {
            int gab = Ca.intl[a][0], Jab = P.g2[gab][2];
            if (std::abs(2 * Jab - on.j2) > tJp or 2 * Jab + on.j2 < tJp) continue;
            int p2 = P.pidx(chk, gab, n);
            if (p2 >= 0) { ac.push_back(a); pc.push_back(p2); cJ.push_back(Jab); }
          }
          if (ac.empty()) continue;
          std::fill(wt.begin(), wt.end(), NAN);   // NaN = not computed yet
          for (size_t q = 0; q < ac.size(); q++)
          {
            const auto *wc = W->colptr(pc[q]);
            double *oc_ = Out.colptr(ac[q]);
            for (size_t r = 0; r < er.size(); r++)
            {
              double &w = wt[rJ[r] * nJ + cJ[q]];
              if (std::isnan(w)) w = (tJp + 1) * six(on.j2, Ce.twoJph, 2 * rJ[r], oc.j2, tJp, 2 * cJ[q]);
              oc_[er[r]] += w * wc[pr[r]];
            }
          }
        }
      }
    }
  };

  // channel pairs (r, c)
  size_t npph = P.pph.size();
  int TY = Y.GetTRank(), pY = Y.GetParity();   // Y changes 2Tz by +-2 TY and the parity by pY
  std::vector<std::array<size_t, 2>> pairs;
  std::vector<long> pairIdx(npph * npph, -1);
  for (size_t r = 0; r < npph; r++)
    for (int id = 0; id < (TY == 0 ? 1 : 2); id++)   // Delta = +TY, -TY; parity-changing only (TY = 0): one partner
    {
      int d = id == 0 ? TY : -TY;
      auto &C = P.pph[r];
      int tz = C.twoTz - 2 * d;
      if (tz < -9 or tz > 9) continue;
      int c = P.pphtab[P.i3(C.twoJph, (C.par + pY) % 2, tz)];
      if (c < 0) continue;
      pairIdx[r * npph + c] = pairs.size(); pairs.push_back({r, (size_t)c});
    }
  // blocked layout of each pair: corner (n1, n2) at boff[n1 norb + n2], leading dimension Cr.pre[n1 norb + n2]
  std::vector<char> blk(pairs.size(), 0);
  std::vector<std::vector<size_t>> boff(pairs.size());
  for (size_t ip = 0; ip < pairs.size(); ip++)
  {
    auto &Cr = P.pph[pairs[ip][0]]; auto &Cc = P.pph[pairs[ip][1]];
    if (Cr.pre.empty() or Cc.pre.empty()) continue;
    boff[ip].assign((size_t)norb * norb, 0);
    size_t used = 0;
    for (int n1 = 0; n1 < norb; n1++)
      for (int n2 = 0; n2 < norb; n2++)
      {
        boff[ip][(size_t)n1 * norb + n2] = used;
        used += (size_t)Cr.pre[(size_t)n1 * norb + n2] * Cc.pre[(size_t)n2 * norb + n1];
      }
    blk[ip] = used < 0.9 * (double)Cr.ext.size() * Cc.ext.size();
    if (not blk[ip]) boff[ip].clear();
  }
  if (std::getenv("IMSRG_PLAN_DEBUG") != nullptr)
  {
    size_t nb = 0; for (char v : blk) nb += v;
    std::cout << "PLANDEBUG comm333_pph tz: " << nb << " of " << pairs.size() << " channel pairs blocked" << std::endl;
  }
  auto bsize = [&](size_t ip) { auto &Cr = P.pph[pairs[ip][0]]; auto &Cc = P.pph[pairs[ip][1]]; size_t z = (size_t)norb * norb - 1;
                                return boff[ip][z] + (size_t)Cr.pre[z] * Cc.pre[z]; };
  std::vector<arma::mat> A(pairs.size());
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ip = 0; ip < pairs.size(); ip++)
  {
    size_t r = pairs[ip][0], c = pairs[ip][1];
    auto &Cr = P.pph[r]; auto &Cc = P.pph[c];
    if (blk[ip]) A[ip].zeros(bsize(ip), 1);
    else if (Cr.intl.empty()) { A[ip].zeros(Cr.ext.size(), Cc.ext.size()); continue; }
    if (Cr.intl.empty()) continue;
    arma::mat Xbar, Ybar;
    fillbar(Xbar, Cr, Cr, [&](int chb, int chk) -> const perfplan::wmat * { return chb == chk ? &WX[chb] : nullptr; });
    fillbar(Ybar, Cc, Cr, [&](int chb, int chk) -> const perfplan::wmat * { int k = WY.idx[(size_t)chb * WY.nch + chk]; return k < 0 ? nullptr : &WY.mats[k]; });
    for (size_t a = 0; a < Cr.intl.size(); a++) Xbar.col(a) *= hY * Cr.occ_int[a];
    if (not blk[ip]) { A[ip] = Xbar * Ybar.t(); continue; }
    for (int n1 = 0; n1 < norb; n1++)   // corners only
      for (int n2 = 0; n2 < norb; n2++)
      {
        int m = Cr.pre[(size_t)n1 * norb + n2], k = Cc.pre[(size_t)n2 * norb + n1];
        if (m == 0 or k == 0) continue;
        // intl (gab, c) contribute only if P(gab, n1) and P(gab, n2) exist: the prefix e(gab) <= E3max - max(e_n1, e_n2)
        Orbit &on1 = ms.GetOrbit(n1); Orbit &on2 = ms.GetOrbit(n2);
        int ki = Cr.ipre[std::max(2 * on1.n + on1.l, 2 * on2.n + on2.l)];
        if (ki == 0) continue;
        arma::mat Ab(A[ip].memptr() + boff[ip][(size_t)n1 * norb + n2], m, k, false, true);
        Ab = Xbar.submat(Cr.noff[n1], 0, Cr.noff[n1] + m - 1, ki - 1) * Ybar.submat(Cc.noff[n2], 0, Cc.noff[n2] + k - 1, ki - 1).t();
      }
  }
  // Zbar(r, c) = A(r, c) - hX hY A(c, r)^T, in place: each pair (r, c) is processed together with its partner (c, r),
  // from the unmodified A of both (one temporary block), so A and Zbar are never held at the same time
  const double hh = hX * hY;
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ip = 0; ip < pairs.size(); ip++)
  {
    size_t r = pairs[ip][0], c = pairs[ip][1];
    long iq = pairIdx[c * npph + r];
    if (iq < 0 or (size_t)iq <= ip) continue;   // no partner (Zbar = A), or handled with the partner
    if (not blk[ip])
    {
      arma::mat tmp = A[ip];
      A[ip] -= hh * A[iq].t();
      A[iq] -= hh * tmp.t();
      continue;
    }
    auto &Cr = P.pph[r]; auto &Cc = P.pph[c];
    for (int n1 = 0; n1 < norb; n1++)   // corner (n1, n2) of (r, c) is the transpose of corner (n2, n1) of (c, r)
      for (int n2 = 0; n2 < norb; n2++)
      {
        int m = Cr.pre[(size_t)n1 * norb + n2], k = Cc.pre[(size_t)n2 * norb + n1];
        if (m == 0 or k == 0) continue;
        arma::mat Ap(A[ip].memptr() + boff[ip][(size_t)n1 * norb + n2], m, k, false, true);
        arma::mat Aq(A[iq].memptr() + boff[iq][(size_t)n2 * norb + n1], k, m, false, true);
        arma::mat tmp = Ap;
        Ap -= hh * Aq.t();
        Aq -= hh * tmp.t();
      }
  }
  std::vector<arma::mat> Zbar = std::move(A);
  Z.profiler.timer["_comm333_pph_plantz_fill_gemm"] += omp_get_wtime() - t; t = omp_get_wtime();

  // back-transformation into the stored blocks (cb, ck) of Z
  std::vector<std::array<size_t, 2>> blocks;
  for (auto &it : Z3.Get_ch_start()) blocks.push_back({it.first.ch_bra, it.first.ch_ket});
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < blocks.size(); ib++)
  {
    size_t cb = blocks[ib][0], ck = blocks[ib][1];
    if (cb == ck) continue;   // T rank > 0: no diagonal blocks
    auto &Tb = ms.GetThreeBodyChannel(cb); auto &Tk = ms.GetThreeBodyChannel(ck);
    if (Tb.twoJ != Tk.twoJ) continue;
    int twoJ = Tb.twoJ;
    size_t nPb = P.Pst[cb].size(), nPk = P.Pst[ck].size(), nb = Tb.GetNumberKets(), nk = Tk.GetNumberKets();
    if (nPb == 0 or nPk == 0 or nb == 0 or nk == 0) continue;
    arma::mat WZ(nPb, nPk, arma::fill::zeros);
    std::vector<std::vector<int>> rowsBy(norb), colsBy(norb);
    for (size_t p = 0; p < nPb; p++) rowsBy[P.Pst[cb][p][1]].push_back(p);
    for (size_t p = 0; p < nPk; p++) colsBy[P.Pst[ck][p][1]].push_back(p);
    std::vector<int> rp, re, rJ, cp, ce, cJ; std::vector<double> wt(nJ * nJ);
    for (int I3 = 0; I3 < norb; I3++)
    {
      auto &Rw = rowsBy[I3];
      if (Rw.empty()) continue;
      Orbit &o3 = ms.GetOrbit(I3);
      int g0 = P.Pst[cb][Rw[0]][0];
      Orbit &o1 = ms.GetOrbit(P.g2[g0][0]); Orbit &o2 = ms.GetOrbit(P.g2[g0][1]);
      for (int I6 = 0; I6 < norb; I6++)
      {
        auto &Cl = colsBy[I6];
        if (Cl.empty()) continue;
        Orbit &o6 = ms.GetOrbit(I6);
        int h0 = P.Pst[ck][Cl[0]][0];
        Orbit &o4 = ms.GetOrbit(P.g2[h0][0]); Orbit &o5 = ms.GetOrbit(P.g2[h0][1]);
        int par_r = (o1.l + o2.l + o6.l) % 2, tz_r = o1.tz2 + o2.tz2 - o6.tz2;
        int par_c = (o4.l + o5.l + o3.l) % 2, tz_c = o4.tz2 + o5.tz2 - o3.tz2;
        for (int tJph = 1; tJph <= P.tJmax; tJph += 2)
        {
          int r = P.pphtab[P.i3(tJph, par_r, tz_r)], c = P.pphtab[P.i3(tJph, par_c, tz_c)];
          if (r < 0 or c < 0) continue;
          long ip = pairIdx[(size_t)r * npph + c];
          if (ip < 0) continue;
          auto &Cr = P.pph[r]; auto &Cc = P.pph[c];
          rp.clear(); re.clear(); rJ.clear(); cp.clear(); ce.clear(); cJ.clear();
          for (int p1 : Rw)
          {
            int g12 = P.Pst[cb][p1][0], J12 = P.g2[g12][2];
            if (std::abs(2 * J12 - o6.j2) > tJph or 2 * J12 + o6.j2 < tJph) continue;
            int e1 = g12 < P.ng3 ? Cr.extIdx[(size_t)g12 * norb + I6] : -1;
            if (e1 >= 0) { rp.push_back(p1); re.push_back(e1); rJ.push_back(J12); }
          }
          if (rp.empty()) continue;
          for (int p2 : Cl)
          {
            int g45 = P.Pst[ck][p2][0], J45 = P.g2[g45][2];
            if (std::abs(2 * J45 - o3.j2) > tJph or 2 * J45 + o3.j2 < tJph) continue;
            int e2 = g45 < P.ng3 ? Cc.extIdx[(size_t)g45 * norb + I3] : -1;
            if (e2 >= 0) { cp.push_back(p2); ce.push_back(e2); cJ.push_back(J45); }
          }
          if (cp.empty()) continue;
          std::fill(wt.begin(), wt.end(), NAN);
          const arma::mat &Zb = Zbar[ip];
          if (blk[ip])   // corner (I6, I3): every read lies inside it (P(g12, I3) and P(g45, I6) exist)
          {
            const double *z0 = Zb.memptr() + boff[ip][(size_t)I6 * norb + I3];
            size_t ld = Cr.pre[(size_t)I6 * norb + I3];
            for (size_t k = 0; k < rp.size(); k++) re[k] -= Cr.noff[I6];
            for (size_t q = 0; q < cp.size(); q++)
            {
              const double *zc = z0 + (size_t)(ce[q] - Cc.noff[I3]) * ld;
              double *wc = WZ.colptr(cp[q]);
              for (size_t k = 0; k < rp.size(); k++)
              {
                double &w = wt[rJ[k] * nJ + cJ[q]];
                if (std::isnan(w)) w = (tJph + 1) * six(o3.j2, twoJ, 2 * rJ[k], o6.j2, tJph, 2 * cJ[q]);
                wc[rp[k]] += w * zc[re[k]];
              }
            }
            continue;
          }
          for (size_t q = 0; q < cp.size(); q++)
          {
            const double *zc = Zb.colptr(ce[q]);
            double *wc = WZ.colptr(cp[q]);
            for (size_t k = 0; k < rp.size(); k++)
            {
              double &w = wt[rJ[k] * nJ + cJ[q]];
              if (std::isnan(w)) w = (tJph + 1) * six(o3.j2, twoJ, 2 * rJ[k], o6.j2, tJph, 2 * cJ[q]);
              wc[rp[k]] += w * zc[re[k]];
            }
          }
        }
      }
    }
    arma::mat Zk = perfplan::AtWA(P.A[cb], WZ, P.A[ck]);
    for (size_t ibra = 0; ibra < nb; ibra++)
      for (size_t iket = 0; iket < nk; iket++)
        if (perfplan::pass3(cutZ, cb, ibra) and perfplan::pass3(cutZ, ck, iket)) Z3.AddToME_pn_ch(cb, ck, ibra, iket, Zk(ibra, iket));
  }
  Z.profiler.timer["_comm333_pph_plantz_back"] += omp_get_wtime() - t;
  perfplan::release_memory();
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
