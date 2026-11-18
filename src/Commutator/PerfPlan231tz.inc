//=====================================================================================================================
//  comm231ss, scalar X and Y changing isospin and/or parity. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
//  Formula as comm231ss (PerfPlan132.inc). Following the structure of the existing comm231ss, for each pair (i, j) of
//  Z and each stored 2B block (ch_bra, ch_ket) of the 2B operator O2 (ab in ch_bra, cd in ch_ket):
//    z_ij += sign sum_J' (2J'+1) sum_{ab,cd} F(ab,cd) [ O2_{ab,cd} O3^{J'}_{(cd)i,(ab)j} - O3^{J'}_{(ab)i,(cd)j} O2_{cd,ab} ]
//  for (O2, O3, sign) = (X, Y, +1) and (Y, X, -1); F = na nb n'c n'd - n'a n'b nc nd (halved for a == b and for c == d,
//  only cd <= ab when ch_bra == ch_ket). O3 is read from the unpacked W blocks (diagonal for X, off-diagonal for Y) as
//  two gathered submatrices G1 = W(P(cd, i), P(ab, j)) and G2 = W(P(ab, i), P(cd, j)). Then Z_ij = z_ij/(2ji+1) and
//  Z_ji = hZ Z_ij.
//  Included after PerfPlan333tz.inc (inside namespace Commutator).
void comm231ss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced() or ms.GetEMax3Body() < ms.GetEmax())   // see comm231ss_plan
  {
    comm231ss(X, Y, Z);
    return;
  }
  int hZ = Z.IsHermitian() ? +1 : -1;
  perfplan::Plan &P = perfplan::get(Z, Z);
  bool x3 = X.ThreeBody.IsAllocated() and X.ThreeBodyNorm() >= 1e-8, y3 = Y.ThreeBody.IsAllocated() and Y.ThreeBodyNorm() >= 1e-8;
  static const std::vector<arma::mat> empty;
  static const std::vector<perfplan::wmat> emptyW;
  static perfplan::WOffCache emptyoff;
  const std::vector<perfplan::wmat> &WX = x3 ? perfplan::unpack(X, P) : emptyW;
  const perfplan::WOffCache &WY = y3 ? perfplan::unpack_off(Y, P) : emptyoff;
  int norb = P.norb;
  // the 3B block W(ch3b, ch3k) of O3 (nullptr if absent)
  auto getW = [&](bool isY, int ch3b, int ch3k) -> const perfplan::wmat *
  {
    if (ch3b < 0 or ch3k < 0) return nullptr;
    if (not isY) return (ch3b == ch3k and WX[ch3b].n_elem > 0) ? &WX[ch3b] : nullptr;
    int k = WY.idx[(size_t)ch3b * WY.nch + ch3k];
    return (k < 0 or WY.mats[k].n_elem == 0) ? nullptr : &WY.mats[k];
  };
  // per 2B block of each 2B operator: kets' pair indices, matrices, weights
  struct Blk { size_t chb, chk; int J, parb, tzb, park, tzk; std::vector<int> gb, gk; arma::mat O_bk, O_kb, F; };
  auto blocks_of = [&](const Operator &O2)
  {
    std::vector<Blk> v;
    for (auto &it : O2.TwoBody.MatEl)
    {
      Blk B; B.chb = it.first[0]; B.chk = it.first[1];
      auto &tb = ms.GetTwoBodyChannel(B.chb); auto &tk = ms.GetTwoBodyChannel(B.chk);
      size_t nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
      if (nb == 0 or nk == 0 or tb.J != tk.J) continue;
      B.J = tb.J;
      Orbit &b1 = ms.GetOrbit(tb.GetKet(0).p); Orbit &b2 = ms.GetOrbit(tb.GetKet(0).q);
      Orbit &k1 = ms.GetOrbit(tk.GetKet(0).p); Orbit &k2 = ms.GetOrbit(tk.GetKet(0).q);
      B.parb = (b1.l + b2.l) % 2; B.tzb = b1.tz2 + b2.tz2; B.park = (k1.l + k2.l) % 2; B.tzk = k1.tz2 + k2.tz2;
      B.gb.resize(nb); B.gk.resize(nk);
      for (size_t k = 0; k < nb; k++) B.gb[k] = P.g2idx[((size_t)tb.GetKet(k).p * norb + tb.GetKet(k).q) * (P.J2max + 1) + B.J];
      for (size_t k = 0; k < nk; k++) B.gk[k] = P.g2idx[((size_t)tk.GetKet(k).p * norb + tk.GetKet(k).q) * (P.J2max + 1) + B.J];
      B.O_bk.set_size(nb, nk); B.O_kb.set_size(nb, nk); B.F.zeros(nb, nk);
      for (size_t ib = 0; ib < nb; ib++)
      {
        Ket &bra = tb.GetKet(ib);
        double na = bra.op->occ, nb_ = bra.oq->occ;
        for (size_t ik = 0; ik < nk; ik++)
        {
          Ket &ket = tk.GetKet(ik);
          B.O_bk(ib, ik) = O2.TwoBody.GetTBME(B.chb, B.chk, bra, ket);
          B.O_kb(ib, ik) = O2.TwoBody.GetTBME(B.chk, B.chb, ket, bra);
          if (B.chb == B.chk and ik > ib) continue;
          double nc = ket.op->occ, nd = ket.oq->occ;
          double pf = na * nb_ * (1 - nc) * (1 - nd) - (1 - na) * (1 - nb_) * nc * nd;
          if (std::abs(pf) < 1e-8) continue;
          if (bra.p == bra.q) pf /= 2;
          if (ket.p == ket.q) pf /= 2;
          B.F(ib, ik) = pf;
        }
      }
      v.push_back(std::move(B));
    }
    return v;
  };
  std::vector<Blk> BX = y3 ? blocks_of(X) : std::vector<Blk>(), BY = x3 ? blocks_of(Y) : std::vector<Blk>();
  // pairs (i, j) of Z, j <= i as in comm231ss
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
    Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
    double zij = 0;
    for (int pass = 0; pass < 2; pass++)   // pass 0: X2 with Y3, pass 1: Y2 with X3 (sign -1)
    {
      const std::vector<Blk> &BB = (pass == 0) ? BX : BY;
      bool isY = (pass == 0);
      for (auto &B : BB)
      {
        size_t nb = B.gb.size(), nk = B.gk.size();
        for (int twoJ = std::abs(2 * B.J - oi.j2); twoJ <= 2 * B.J + oi.j2; twoJ += 2)
        {
          // G1(cd, ab) = W(P(cd, i), P(ab, j)),  G2(ab, cd) = W(P(ab, i), P(cd, j))
          const perfplan::wmat *W1 = getW(isY, P.ch3tab[P.i3(twoJ, (B.park + oi.l) % 2, B.tzk + oi.tz2)], P.ch3tab[P.i3(twoJ, (B.parb + oj.l) % 2, B.tzb + oj.tz2)]);
          const perfplan::wmat *W2 = getW(isY, P.ch3tab[P.i3(twoJ, (B.parb + oi.l) % 2, B.tzb + oi.tz2)], P.ch3tab[P.i3(twoJ, (B.park + oj.l) % 2, B.tzk + oj.tz2)]);
          if (W1 == nullptr and W2 == nullptr) continue;
          int c1b = P.ch3tab[P.i3(twoJ, (B.park + oi.l) % 2, B.tzk + oi.tz2)], c1k = P.ch3tab[P.i3(twoJ, (B.parb + oj.l) % 2, B.tzb + oj.tz2)];
          int c2b = P.ch3tab[P.i3(twoJ, (B.parb + oi.l) % 2, B.tzb + oi.tz2)], c2k = P.ch3tab[P.i3(twoJ, (B.park + oj.l) % 2, B.tzk + oj.tz2)];
          double s = 0;
          for (size_t ib = 0; ib < nb; ib++)
          {
            if (B.gb[ib] < 0) continue;
            for (size_t ik = 0; ik < nk; ik++)
            {
              double f = B.F(ib, ik);
              if (f == 0 or B.gk[ik] < 0) continue;
              double g1 = 0, g2 = 0;
              if (W1)
              {
                int p = P.pidx(c1b, B.gk[ik], i), q = P.pidx(c1k, B.gb[ib], j);
                if (p >= 0 and q >= 0) g1 = (*W1)(p, q);
              }
              if (W2)
              {
                int p = P.pidx(c2b, B.gb[ib], i), q = P.pidx(c2k, B.gk[ik], j);
                if (p >= 0 and q >= 0) g2 = (*W2)(p, q);
              }
              s += f * (B.O_bk(ib, ik) * g1 - g2 * B.O_kb(ib, ik));
            }
          }
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
