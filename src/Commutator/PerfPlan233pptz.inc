//=====================================================================================================================
//  comm233_pp_hhss, scalar X and isospin-changing Y. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
// comm233_pp_hhss for a scalar X and a channel-changing Y (J rank 0, parity 0, T rank > 0), precompiled bookkeeping.
// Generalises comm233_pp_hhss_plan to the stored blocks (cb, ck) of Z (cb != ck):
//   M(x, y) = A(x)^T [ D_X(x) (R Y)(x, y) - D_Y(x <- y) (R X)(y) ],     Z(cb, ck) = M(cb, ck) - hX hY M(ck, cb)^T
// D is the 2B operator embedded in the P basis (same third orbit c and pair J; 0.5 symm (1 - n_a - n_b) on the
// columns, sqrt2 for identical pairs), D_Y couples pair channels that differ in Tz.
// Included after PerfPlan233pp.inc and PerfPlan333tz.inc (inside namespace Commutator).
void comm233_pp_hhss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced())
  {
    comm233_pp_hhss_tz(X, Y, Z);
    return;
  }
  auto &Z3 = Z.ThreeBody;
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  if (not(x3 or y3)) return;
  int hXY = (X.IsHermitian() ? 1 : -1) * (Y.IsHermitian() ? 1 : -1), hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(Z, Z);
  const perfplan::Cut3 &cutZ = perfplan::cut3(*Z.modelspace);
  perfplan::Plan233pp &Q = perfplan::get233pp(Z, P);
  size_t nch3 = ms.GetNumberThreeBodyChannels(), nch2 = ms.GetNumberTwoBodyChannels();
  static const std::vector<arma::mat> empty;
  static perfplan::WOffCache emptyoff;
  const std::vector<arma::mat> &RX = x3 ? perfplan::unpack_rx(X, P) : empty;
  const perfplan::WOffCache &WY = y3 ? perfplan::unpack_off(Y, P) : emptyoff;
  const std::vector<arma::mat> WYrk = perfplan::unpack_off_rk(WY, P);   // R(cb) Y, held while this term runs
  std::vector<arma::mat> YF(nch2 * nch2);
  for (auto &kv : Y.TwoBody.MatEl)
  {
    size_t cb = kv.first[0], ck = kv.first[1];
    YF[cb * nch2 + ck] = kv.second;
    if (cb != ck) YF[ck * nch2 + cb] = hY * kv.second.t();
  }
  // group lookup per channel: (c, J) -> group index
  auto groupKey = [&](size_t ch, const perfplan::Grp233pp &g) { int ip = g.rows(0); return std::make_pair(P.Pst[ch][ip][1], P.g2[P.Pst[ch][ip][0]][2]); };
  auto Mxy = [&](size_t x, size_t y) -> arma::mat
  {
    size_t nPx = P.Pst[x].size(), nky = ms.GetThreeBodyChannel(y).GetNumberKets();
    arma::mat G(nPx, nky, arma::fill::zeros);
    int kxy = WY.idx.empty() ? -1 : WY.idx[x * WY.nch + y];
    std::map<std::pair<int, int>, int> ymap;
    for (size_t ig = 0; ig < Q.grp[y].size(); ig++) ymap[groupKey(y, Q.grp[y][ig])] = ig;
    for (auto &g : Q.grp[x])
    {
      arma::mat acc(g.rows.n_elem, nky, arma::fill::zeros);
      if (y3 and kxy >= 0)
      {
        arma::mat Dx = X.TwoBody.GetMatrix(g.ch2, g.ch2).submat(g.loc, g.loc);
        Dx.each_row() %= g.scale; Dx.each_col() %= g.rnorm;
        acc += Dx * WYrk[kxy].rows(g.rows);
      }
      if (x3 and RX[y].n_elem > 0)
      {
        auto it = ymap.find(groupKey(x, g));
        if (it != ymap.end())
        {
          auto &h = Q.grp[y][it->second];
          const arma::mat &M = YF[(size_t)g.ch2 * nch2 + h.ch2];
          if (M.n_elem > 0)
          {
            arma::mat Dy = M.submat(g.loc, h.loc);
            Dy.each_row() %= h.scale; Dy.each_col() %= g.rnorm;
            acc -= Dy * RX[y].rows(h.rows);
          }
        }
      }
      G.rows(g.rows) = acc;
    }
    return perfplan::AtM(P.A[x], G);
  };
  std::vector<std::array<size_t, 2>> blocks;
  for (auto &it : Z3.Get_ch_start()) blocks.push_back({it.first.ch_bra, it.first.ch_ket});
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < blocks.size(); ib++)
  {
    size_t cb = blocks[ib][0], ck = blocks[ib][1];
    if (cb == ck) continue;
    size_t nb = ms.GetThreeBodyChannel(cb).GetNumberKets(), nk = ms.GetThreeBodyChannel(ck).GetNumberKets();
    if (P.Pst[cb].empty() or P.Pst[ck].empty() or nb == 0 or nk == 0) continue;
    arma::mat Zk = Mxy(cb, ck) - hXY * Mxy(ck, cb).t();
    for (size_t ibra = 0; ibra < nb; ibra++)
      for (size_t iket = 0; iket < nk; iket++)
        if (perfplan::pass3(cutZ, cb, ibra) and perfplan::pass3(cutZ, ck, iket)) Z3.AddToME_pn_ch(cb, ck, ibra, iket, Zk(ibra, iket));
  }
  (void)nch3;
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
