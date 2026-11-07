//=====================================================================================================================
//  comm132ss and comm231ss with precompiled bookkeeping (scalar X, Y, default cuts). Notation: PerfPlan.inc.
//=====================================================================================================================
//  The 3B matrix elements GetME_pn(J, J, J', ...) of the original terms are read as gathered blocks of the unpacked
//  W = R X R^T (W^{J'}(P((ij)J, b), P((kl)J, a)) = X^{J,J,J'}_{ijb,kla}).
//
//  comm132ss (IMSRG3Commutators.cc):
//    Z^J_ijkl = sum_ab (na n'b - n'a nb) sum_J' (2J'+1)/(2J+1) (X_ab Y^{J,J,J'}_{ijb,kla} - Y_ab X^{J,J,J'}_{ijb,kla})
//    Per 2B channel: Step 1, for each (a, b, J') the block W(P(ij, b), P(kl, a)) over all kets (ij), (kl);
//                    Step 2, M += (occupations, X_ab, (2J'+1)/(2J+1)) * block;  Z2 += M.
//
//  comm231ss (further down):
//    Z_ij = 1/4 sum_abcd (na nb n'c n'd) (2J+1)/(2ji+1) (X^{J1}_{ab,cd} Y^{J1 J1 J}_{cdi,abj} - Y^{J1 J1 J}_{abi,cdj} X^{J1}_{cd,ab})
//    Per (i, j) and 2B channel J1: G(cd, ab) = W^J(P(cd, i), P(ab, j)), then Z_ij += sum (F o (X2 o G^T - G o X2^T)),
//    F the occupation weights of the ket pairs (o: elementwise product).
//  Included after PerfPlan.inc (inside namespace Commutator).

void comm132ss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0)
  {
    comm132ss(X, Y, Z);
    return;
  }
  if (X.GetParticleRank() < 3 and Y.GetParticleRank() < 3) return;
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  if (not(x3 or y3)) return;
  int hZ = Z.IsHermitian() ? +1 : -1;
  perfplan::Plan &P = perfplan::get(Z, Z);
  static const std::vector<arma::mat> empty;
  static const std::vector<perfplan::wmat> emptyW;
  const std::vector<perfplan::wmat> &WX = x3 ? perfplan::unpack(X, P) : emptyW;
  const std::vector<perfplan::wmat> &WY = y3 ? perfplan::unpack(Y, P) : emptyW;
  int norb = P.norb;
  auto &Z2 = Z.TwoBody;
  size_t nch = ms.GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ch = 0; ch < nch; ch++)
  {
    auto &tbc = ms.GetTwoBodyChannel(ch);
    int J = tbc.J;
    size_t nk = tbc.GetNumberKets();
    if (nk == 0) continue;
    std::vector<int> gk(nk);
    for (size_t k = 0; k < nk; k++) gk[k] = P.g2idx[((size_t)tbc.GetKet(k).p * norb + tbc.GetKet(k).q) * (P.J2max + 1) + J];
    Orbit &o1 = ms.GetOrbit(tbc.GetKet(0).p); Orbit &o2 = ms.GetOrbit(tbc.GetKet(0).q);
    int par12 = (o1.l + o2.l) % 2, tz12 = o1.tz2 + o2.tz2;
    arma::mat M(nk, nk, arma::fill::zeros);
    std::vector<arma::uword> rb, rk, ca, ck;
    for (int a = 0; a < norb; a++)
    {
      Orbit &oa = ms.GetOrbit(a);
      for (auto b : X.GetOneBodyChannel(oa.l, oa.j2, oa.tz2))
      {
        Orbit &ob = ms.GetOrbit(b);
        double occfactor = oa.occ - ob.occ;
        if (std::abs(occfactor) < 1e-6) continue;
        double xab = X.OneBody(a, b), yab = Y.OneBody(a, b);
        for (int twoJ = std::abs(oa.j2 - 2 * J); twoJ <= oa.j2 + 2 * J; twoJ += 2)
        {
          int ch3 = P.ch3tab[P.i3(twoJ, (par12 + oa.l) % 2, tz12 + oa.tz2)];
          if (ch3 < 0) continue;
          rb.clear(); rk.clear(); ca.clear(); ck.clear();
          for (size_t k = 0; k < nk; k++)
          {
            if (gk[k] < 0) continue;
            int pr = P.pidx(ch3, gk[k], b), pc = P.pidx(ch3, gk[k], a);
            if (pr >= 0) { rb.push_back(pr); rk.push_back(k); }
            if (pc >= 0) { ca.push_back(pc); ck.push_back(k); }
          }
          if (rb.empty() or ca.empty()) continue;
          arma::uvec urb(rb), urk(rk), uca(ca), uck(ck);
          double w = occfactor * (twoJ + 1.) / (2 * J + 1);
          arma::mat blk(rb.size(), ca.size(), arma::fill::zeros);
          if (y3 and xab != 0 and WY[ch3].n_elem > 0) blk += (w * xab) * PLAN_W_D(WY[ch3].submat(urb, uca));
          if (x3 and yab != 0 and WX[ch3].n_elem > 0) blk -= (w * yab) * PLAN_W_D(WX[ch3].submat(urb, uca));   // empty W: no stored kets
          M.submat(urk, uck) += blk;
        }
      }
    }
    for (size_t ibra = 0; ibra < nk; ibra++)
    {
      Ket &bra = tbc.GetKet(ibra);
      for (size_t iket = 0; iket <= ibra; iket++)
      {
        Ket &ket = tbc.GetKet(iket);
        double zijkl = M(ibra, iket) / sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
        Z2.AddToTBMENonHerm(ch, ch, ibra, iket, zijkl);
        if (ibra != iket) Z2.AddToTBMENonHerm(ch, ch, iket, ibra, zijkl * hZ);
      }
    }
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}

//   Z_ij = sum_{ab >= cd (stored order)} f_abcd sum_J' (2J'+1)/(2j_i+1) ( X_abcd Y_{(cd)i,(ab)j}^{J'} - Y_{(ab)i,(cd)j}^{J'} X_cdab )  - (X <-> Y)
//   f_abcd = (n_a n_b nbar_c nbar_d - nbar_a nbar_b n_c n_d), halved for a == b and for c == d
void comm231ss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  // emax_3body < emax: the existing comm231ss is faster here (this kernel loops over the full 2B space), so use it
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0
      or ms.GetEMax3Body() < ms.GetEmax())
  {
    comm231ss(X, Y, Z);
    return;
  }
  int hZ = Z.IsHermitian() ? +1 : -1;
  perfplan::Plan &P = perfplan::get(Z, Z);
  bool use[2] = {Y.ThreeBody.IsAllocated() and Y.ThreeBodyNorm() >= 1e-8, X.ThreeBody.IsAllocated() and X.ThreeBodyNorm() >= 1e-8};
  static const std::vector<arma::mat> empty;
  static const std::vector<perfplan::wmat> emptyW;
  const std::vector<perfplan::wmat> &WY = use[0] ? perfplan::unpack(Y, P) : emptyW;
  const std::vector<perfplan::wmat> &WX = use[1] ? perfplan::unpack(X, P) : emptyW;
  int norb = P.norb;
  size_t nch = ms.GetNumberTwoBodyChannels();
  // per 2B channel: the 2B matrices (as GetTBME), the weight F (lower triangle of the ket pairs) and the ket pairs
  struct ChData { std::vector<int> g; arma::mat X2, Y2, F; int par = 0, tz = 0; };
  std::vector<ChData> CD(nch);
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ch = 0; ch < nch; ch++)
  {
    auto &tbc = ms.GetTwoBodyChannel(ch);
    size_t nk = tbc.GetNumberKets();
    if (nk == 0) continue;
    auto &D = CD[ch];
    D.g.resize(nk); D.X2.set_size(nk, nk); D.Y2.set_size(nk, nk); D.F.zeros(nk, nk);
    Orbit &o1 = ms.GetOrbit(tbc.GetKet(0).p); Orbit &o2 = ms.GetOrbit(tbc.GetKet(0).q);
    D.par = (o1.l + o2.l) % 2; D.tz = o1.tz2 + o2.tz2;
    for (size_t k = 0; k < nk; k++) D.g[k] = P.g2idx[((size_t)tbc.GetKet(k).p * norb + tbc.GetKet(k).q) * (P.J2max + 1) + tbc.J];
    for (size_t ib = 0; ib < nk; ib++)
    {
      Ket &bra = tbc.GetKet(ib);
      double na = bra.op->occ, nb = bra.oq->occ;
      for (size_t ik = 0; ik < nk; ik++)
      {
        Ket &ket = tbc.GetKet(ik);
        D.X2(ib, ik) = X.TwoBody.GetTBME(ch, ch, bra, ket);
        D.Y2(ib, ik) = Y.TwoBody.GetTBME(ch, ch, bra, ket);
        if (ik > ib) continue;
        double nc = ket.op->occ, nd = ket.oq->occ;
        double pf = na * nb * (1 - nc) * (1 - nd) - (1 - na) * (1 - nb) * nc * nd;
        if (std::abs(pf) < 1e-8) continue;
        if (bra.p == bra.q) pf /= 2;
        if (ket.p == ket.q) pf /= 2;
        D.F(ib, ik) = pf;
      }
    }
  }
  std::vector<std::array<int, 2>> ij;
  for (int i = 0; i < norb; i++)
  {
    Orbit &oi = ms.GetOrbit(i);
    for (auto j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      if ((int)j <= i) ij.push_back({i, (int)j});
  }
  std::vector<double> zres(ij.size(), 0.0);
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t iij = 0; iij < ij.size(); iij++)
  {
    int i = ij[iij][0], j = ij[iij][1];
    Orbit &oi = ms.GetOrbit(i);
    double zij = 0;
    std::vector<arma::uword> ri, rki, cj, ckj;
    for (size_t ch = 0; ch < nch; ch++)
    {
      auto &D = CD[ch];
      size_t nk = D.g.size();
      if (nk == 0) continue;
      int J = ms.GetTwoBodyChannel(ch).J;
      for (int twoJ = std::abs(2 * J - oi.j2); twoJ <= 2 * J + oi.j2; twoJ += 2)
      {
        int ch3 = P.ch3tab[P.i3(twoJ, (D.par + oi.l) % 2, D.tz + oi.tz2)];
        if (ch3 < 0) continue;
        ri.clear(); rki.clear(); cj.clear(); ckj.clear();
        for (size_t k = 0; k < nk; k++)
        {
          if (D.g[k] < 0) continue;
          int p1 = P.pidx(ch3, D.g[k], i), p2 = P.pidx(ch3, D.g[k], j);
          if (p1 >= 0) { ri.push_back(p1); rki.push_back(k); }
          if (p2 >= 0) { cj.push_back(p2); ckj.push_back(k); }
        }
        if (ri.empty() or cj.empty()) continue;
        arma::uvec uri(ri), urki(rki), ucj(cj), uckj(ckj);
        for (int pass = 0; pass < 2; pass++)   // pass 0: X2 Y3, pass 1: Y2 X3 (with a minus sign)
        {
          if (not use[pass]) continue;
          const perfplan::wmat &W3 = (pass == 0) ? WY[ch3] : WX[ch3];
          if (W3.n_elem == 0) continue;   // channel without stored kets
          const arma::mat &O2 = (pass == 0) ? D.X2 : D.Y2;
          // G(cd, ab) = W(P(cd, i), P(ab, j)), embedded in the full (nk x nk) ket-pair space
          arma::mat G(nk, nk, arma::fill::zeros);
          G.submat(urki, uckj) = PLAN_W_D(W3.submat(uri, ucj));
          double s = arma::accu(D.F % (O2 % G.t() - G % O2.t()));
          zij += ((pass == 0) ? 1.0 : -1.0) * (twoJ + 1) * s;
        }
      }
    }
    zres[iij] = zij;
  }
  for (size_t iij = 0; iij < ij.size(); iij++)
  {
    int i = ij[iij][0], j = ij[iij][1];
    double zij = zres[iij] / (ms.GetOrbit(i).j2 + 1.0);
    Z.OneBody(i, j) += zij;
    if (i != j) Z.OneBody(j, i) += hZ * zij;
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
