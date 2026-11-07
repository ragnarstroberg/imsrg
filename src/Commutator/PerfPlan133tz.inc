//=====================================================================================================================
//  comm133ss, scalar X and isospin-changing Y. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
// comm133ss for a scalar X and a channel-changing Y (J rank 0, parity 0, T rank > 0), precompiled bookkeeping.
// Same maths as the non-channel-diagonal branch of comm133ss, per stored block (cb, ck), cb != ck:
//   Z = X1bra(cb) Y3(cb,ck) - Y1bra(cb<-ck) X3(ck) - Y3(cb,ck) X1ket(ck) + X3(cb) Y1ket(cb<-ck)
// where the 1B operators embedded in the 3B space are built exactly as in the original (bra side: O1(orb, a) with the
// bra orbit replaced; ket side: O1(a, orb) with the ket orbit replaced), as sparse matrices with a cached pattern.
// Included after PerfPlan133.inc and PerfPlan333tz.inc (inside namespace Commutator).
namespace perfplan
{
  struct Emb133 { size_t nr = 0, nc = 0; arma::uvec rowind, colptr; std::vector<arma::uword> slot; std::vector<std::array<int, 2>> o12; std::vector<double> coef; };
  // bra side (rows = kets of chR, recoupled into chC) or ket side (cols = kets of chC, recoupled into chR)
  static Emb133 build133(const Operator &O, ThreeBodyME &Z3, ModelSpace &ms, size_t chR, size_t chC, bool braSide)
  {
    struct C { arma::uword col, row; int o1, o2; double coef; };
    std::vector<C> cs;
    std::vector<size_t> klist; std::vector<double> rlist;
    auto &TR = ms.GetThreeBodyChannel(chR); auto &TC = ms.GetThreeBodyChannel(chC);
    auto &Tsrc = braSide ? TR : TC;
    size_t chDst = braSide ? chC : chR, nDst = braSide ? TC.GetNumberKets() : TR.GetNumberKets();
    int twoJ = Tsrc.twoJ;
    for (size_t is = 0; is < Tsrc.GetNumberKets(); is++)
    {
      Ket3 &k3 = Tsrc.GetKet(is);
      int Jij = k3.Jpq;
      std::array<size_t, 3> ijk = {k3.p, k3.q, k3.r};
      for (int iperm = 0; iperm < 3; iperm++)
      {
        std::array<size_t, 3> IJK = ijk;
        Orbit &orep = ms.GetOrbit(ijk[iperm]);
        for (auto a : O.GetOneBodyChannel(orep.l, orep.j2, orep.tz2))
        {
          IJK[iperm] = a;
          if (not Z3.IsKetValid(Jij, twoJ, IJK[0], IJK[1], IJK[2])) continue;
          klist.clear(); rlist.clear();
          size_t chr = Z3.GetKetIndex_withRecoupling(Jij, twoJ, IJK[0], IJK[1], IJK[2], klist, rlist);
          if (chr != chDst) continue;
          for (size_t t = 0; t < klist.size(); t++)
          {
            if (klist[t] >= nDst) continue;
            if (braSide) cs.push_back({(arma::uword)klist[t], (arma::uword)is, (int)ijk[iperm], (int)a, rlist[t]});
            else cs.push_back({(arma::uword)is, (arma::uword)klist[t], (int)a, (int)ijk[iperm], rlist[t]});
          }
        }
      }
    }
    Emb133 E; E.nr = TR.GetNumberKets(); E.nc = TC.GetNumberKets();
    std::vector<size_t> ord(cs.size());
    for (size_t t = 0; t < ord.size(); t++) ord[t] = t;
    std::stable_sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return cs[x].col != cs[y].col ? cs[x].col < cs[y].col : cs[x].row < cs[y].row; });
    std::vector<arma::uword> rows, colptr(E.nc + 1, 0);
    E.slot.resize(cs.size()); E.o12.resize(cs.size()); E.coef.resize(cs.size());
    for (size_t u = 0; u < ord.size(); u++)
    {
      const C &c = cs[ord[u]];
      if (u == 0 or cs[ord[u - 1]].col != c.col or cs[ord[u - 1]].row != c.row) { rows.push_back(c.row); colptr[c.col + 1]++; }
      E.slot[u] = rows.size() - 1; E.o12[u] = {c.o1, c.o2}; E.coef[u] = c.coef;
    }
    for (size_t k = 0; k < E.nc; k++) colptr[k + 1] += colptr[k];
    E.rowind = arma::uvec(rows); E.colptr = arma::uvec(colptr);
    return E;
  }
  static arma::sp_mat value133(const Emb133 &E, const arma::mat &O1)
  {
    arma::vec v(E.rowind.n_elem, arma::fill::zeros);
    for (size_t t = 0; t < E.slot.size(); t++) v(E.slot[t]) += O1(E.o12[t][0], E.o12[t][1]) * E.coef[t];
    return arma::sp_mat(E.rowind, E.colptr, v, E.nr, E.nc, false);
  }
  struct Plan133tz
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1, TY = -1, pY = -1;
    std::vector<std::array<size_t, 2>> blocks;
    std::vector<std::array<Emb133, 4>> emb;   // per block: X1bra(cb), Y1bra(cb<-ck), X1ket(ck), Y1ket(cb<-ck)
  };
  static Plan133tz &get133tz(const Operator &X, const Operator &Y, Operator &Z, Plan &P)
  {
    static Plan133tz Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    int TY = Y.GetTRank(), pY = Y.GetParity();
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ and Q.TY == TY and Q.pY == pY) return Q;
    Q = Plan133tz(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max(); Q.TY = TY; Q.pY = pY;
    for (auto &it : Z.ThreeBody.Get_ch_start())
      if (it.first.ch_bra != it.first.ch_ket) Q.blocks.push_back({it.first.ch_bra, it.first.ch_ket});
    Q.emb.resize(Q.blocks.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ib = 0; ib < Q.blocks.size(); ib++)
    {
      size_t cb = Q.blocks[ib][0], ck = Q.blocks[ib][1];
      Q.emb[ib][0] = build133(X, Z.ThreeBody, ms, cb, cb, true);
      Q.emb[ib][1] = build133(Y, Z.ThreeBody, ms, cb, ck, true);
      Q.emb[ib][2] = build133(X, Z.ThreeBody, ms, ck, ck, false);
      Q.emb[ib][3] = build133(Y, Z.ThreeBody, ms, cb, ck, false);
    }
    return Q;
  }
} // namespace perfplan

void comm133ss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced())
  {
    comm133ss(X, Y, Z);
    return;
  }
  auto &Z3 = Z.ThreeBody;
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  perfplan::Plan &P = perfplan::get(Z, Z);
  const perfplan::Cut3 &cutZ = perfplan::cut3(*Z.modelspace);
  perfplan::Plan133tz &Q = perfplan::get133tz(X, Y, Z, P);
  static const std::vector<arma::mat> empty;
  static perfplan::WOffCache emptyoff;
  const std::vector<arma::mat> &XK = x3 ? perfplan::unpack_k(X, P) : empty;
  const perfplan::WOffCache &WY = y3 ? perfplan::unpack_off(Y, P) : emptyoff;
  const arma::mat &X1 = X.OneBody, &Y1 = Y.OneBody;
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < Q.blocks.size(); ib++)
  {
    size_t cb = Q.blocks[ib][0], ck = Q.blocks[ib][1];
    size_t nb = ms.GetThreeBodyChannel(cb).GetNumberKets(), nk = ms.GetThreeBodyChannel(ck).GetNumberKets();
    if (nb == 0 or nk == 0) continue;
    arma::mat Zk(nb, nk, arma::fill::zeros);
    int k = y3 ? WY.idx[cb * WY.nch + ck] : -1;
    if (k >= 0 and WY.kk[k].n_elem > 0)
    {
      const arma::mat &Y3 = WY.kk[k];
      Zk += perfplan::SD(perfplan::value133(Q.emb[ib][0], X1), Y3);
      Zk -= Y3 * perfplan::value133(Q.emb[ib][2], X1);
    }
    if (x3)
    {
      if (XK[ck].n_elem > 0) Zk -= perfplan::SD(perfplan::value133(Q.emb[ib][1], Y1), XK[ck]);
      if (XK[cb].n_elem > 0) Zk += XK[cb] * perfplan::value133(Q.emb[ib][3], Y1);
    }
    for (size_t ibra = 0; ibra < nb; ibra++)
      for (size_t iket = 0; iket < nk; iket++)
        if (perfplan::pass3(cutZ, cb, ibra) and perfplan::pass3(cutZ, ck, iket)) Z3.AddToME_pn_ch(cb, ck, ibra, iket, Zk(ibra, iket));
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
