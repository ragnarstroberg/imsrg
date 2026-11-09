//=====================================================================================================================
//  comm132ss, scalar X and Y changing isospin and/or parity. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
//  Formula as comm132ss (PerfPlan132.inc). Per stored block (chb, chk) of Z2 (bra pairs ij in chb, ket pairs kl in chk):
//    Z^J_{ij,kl} = sum_ab (na - nb) sum_J' (2J'+1)/(2J+1) [ X_ab W_Y^{J'}(P(ij, b), P(kl, a)) - Y_ab W_X^{J'}(P(ij, b), P(kl, a)) ]
//  In the first term a, b are in the same orbit channel (X scalar) and W_Y is an off-diagonal block of the unpacked
//  Y (3B channels of (ij)b and (kl)a differ); in the second a, b are connected by Y and W_X is diagonal (the two
//  product states are in the same 3B channel).
//  Included after PerfPlan333tz.inc (inside namespace Commutator).
void comm132ss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced())
  {
    comm132ss(X, Y, Z);
    return;
  }
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  if (not(x3 or y3)) return;
  perfplan::Plan &P = perfplan::get(Z, Z);
  static const std::vector<arma::mat> empty;
  static const std::vector<perfplan::wmat> emptyW;
  static perfplan::WOffCache emptyoff;
  const std::vector<perfplan::wmat> &WX = x3 ? perfplan::unpack(X, P) : emptyW;
  const perfplan::WOffCache &WY = y3 ? perfplan::unpack_off(Y, P) : emptyoff;
  int norb = P.norb;
  auto &Z2 = Z.TwoBody;
  std::vector<std::array<size_t, 2>> blocks;
  for (auto &it : Z2.MatEl) if (it.first[0] != it.first[1]) blocks.push_back({it.first[0], it.first[1]});
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ibk = 0; ibk < blocks.size(); ibk++)
  {
    size_t chb = blocks[ibk][0], chk = blocks[ibk][1];
    auto &tb = ms.GetTwoBodyChannel(chb); auto &tk = ms.GetTwoBodyChannel(chk);
    if (tb.J != tk.J) continue;
    int J = tb.J;
    size_t nb = tb.GetNumberKets(), nk = tk.GetNumberKets();
    if (nb == 0 or nk == 0) continue;
    std::vector<int> gb(nb), gk(nk);
    for (size_t k = 0; k < nb; k++) gb[k] = P.g2idx[((size_t)tb.GetKet(k).p * norb + tb.GetKet(k).q) * (P.J2max + 1) + J];
    for (size_t k = 0; k < nk; k++) gk[k] = P.g2idx[((size_t)tk.GetKet(k).p * norb + tk.GetKet(k).q) * (P.J2max + 1) + J];
    Orbit &b1 = ms.GetOrbit(tb.GetKet(0).p); Orbit &b2 = ms.GetOrbit(tb.GetKet(0).q);
    Orbit &k1 = ms.GetOrbit(tk.GetKet(0).p); Orbit &k2 = ms.GetOrbit(tk.GetKet(0).q);
    int parb = (b1.l + b2.l) % 2, tzb = b1.tz2 + b2.tz2, park = (k1.l + k2.l) % 2, tzk = k1.tz2 + k2.tz2;
    arma::mat M(nb, nk, arma::fill::zeros);
    // one (a, b, twoJ') contribution: M += w * W(P(ij, b), P(kl, a)) over the ket pairs of the two channels
    auto add = [&](const perfplan::wmat &W, int ch3b, int ch3k, int a, int b, double w)
    {
      std::vector<arma::uword> rb, rkk, ca, ckk;
      for (size_t k = 0; k < nb; k++)
      {
        if (gb[k] < 0) continue;
        int pr = P.pidx(ch3b, gb[k], b);
        if (pr >= 0) { rb.push_back(pr); rkk.push_back(k); }
      }
      for (size_t k = 0; k < nk; k++)
      {
        if (gk[k] < 0) continue;
        int pc = P.pidx(ch3k, gk[k], a);
        if (pc >= 0) { ca.push_back(pc); ckk.push_back(k); }
      }
      if (rb.empty() or ca.empty()) return;
      M.submat(arma::uvec(rkk), arma::uvec(ckk)) += w * PLAN_W_D(W.submat(arma::uvec(rb), arma::uvec(ca)));
    };
    for (int a = 0; a < norb; a++)
    {
      Orbit &oa = ms.GetOrbit(a);
      // X_ab Y_{(ij)b,(kl)a}: a, b in the same orbit channel
      if (y3)
        for (auto b : X.GetOneBodyChannel(oa.l, oa.j2, oa.tz2))
        {
          Orbit &ob = ms.GetOrbit(b);
          double occ = oa.occ - ob.occ, xab = X.OneBody(a, b);
          if (std::abs(occ) < 1e-6 or xab == 0) continue;
          for (int twoJ = std::abs(oa.j2 - 2 * J); twoJ <= oa.j2 + 2 * J; twoJ += 2)
          {
            int ch3b = P.ch3tab[P.i3(twoJ, (parb + ob.l) % 2, tzb + ob.tz2)], ch3k = P.ch3tab[P.i3(twoJ, (park + oa.l) % 2, tzk + oa.tz2)];
            if (ch3b < 0 or ch3k < 0) continue;
            int iw = WY.idx[(size_t)ch3b * WY.nch + ch3k];
            if (iw < 0 or WY.mats[iw].n_elem == 0) continue;
            add(WY.mats[iw], ch3b, ch3k, a, b, occ * (twoJ + 1.) / (2 * J + 1) * xab);
          }
        }
      // - Y_ab X_{(ij)b,(kl)a}: a, b connected by Y, both product states in the same 3B channel
      if (x3)
        for (auto b : Y.GetOneBodyChannel(oa.l, oa.j2, oa.tz2))
        {
          Orbit &ob = ms.GetOrbit(b);
          double occ = oa.occ - ob.occ, yab = Y.OneBody(a, b);
          if (std::abs(occ) < 1e-6 or yab == 0) continue;
          for (int twoJ = std::abs(oa.j2 - 2 * J); twoJ <= oa.j2 + 2 * J; twoJ += 2)
          {
            int ch3b = P.ch3tab[P.i3(twoJ, (parb + ob.l) % 2, tzb + ob.tz2)], ch3k = P.ch3tab[P.i3(twoJ, (park + oa.l) % 2, tzk + oa.tz2)];
            if (ch3b < 0 or ch3b != ch3k or WX[ch3b].n_elem == 0) continue;
            add(WX[ch3b], ch3b, ch3b, a, b, -occ * (twoJ + 1.) / (2 * J + 1) * yab);
          }
        }
    }
    for (size_t ibra = 0; ibra < nb; ibra++)
    {
      Ket &bra = tb.GetKet(ibra);
      for (size_t iket = 0; iket < nk; iket++)
      {
        Ket &ket = tk.GetKet(iket);
        Z2.AddToTBMENonHerm(chb, chk, ibra, iket, M(ibra, iket) / sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq())));
      }
    }
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
