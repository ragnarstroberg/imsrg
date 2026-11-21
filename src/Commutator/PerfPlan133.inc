//=====================================================================================================================
//  comm133ss with precompiled bookkeeping (scalar X, Y, default dE3max / occnat3 cuts). Notation: PerfPlan.inc.
//=====================================================================================================================
//  Formula (IMSRG3Commutators.cc):
//    Z_ijklmn = sum_a P(ij/k) X_{ka} Y_{ija,lmn} - P(lm/n) Y_{ijk,lma} X_{an}     - (X <-> Y)
//  Per 3B channel this is a matrix product in the stored kets, as in the original:
//    Z3 = X1 Y3 - Y1 X3,     Z3 -= hX hY Z3^T,
//  where X1 (n_k x n_k) is the 1B operator embedded in the 3B space: X1(k', k) = sum over the three positions of the
//  orbit replaced in |k>, with the recoupling coefficients of GetKetIndex_withRecoupling (~10 nonzeros per row), and
//  Y3, X3 are the stored channel blocks (the K of the unpack cache).
//  Step 1  (plan, once per model space) the sparsity pattern of X1 and, per nonzero, which 1B element and coefficient.
//  Step 2  (per call) fill the values of X1 and Y1 from X.OneBody, Y.OneBody.
//  Step 3  Z3 = X1 Y3 - Y1 X3 column by column (contiguous reads; armadillo's sparse x dense walks dense rows), then
//          Z3 - hX hY Z3^T into the stored block.
//  Included after PerfPlan.inc (inside namespace Commutator).
namespace perfplan
{
  struct Ch133
  {
    arma::uvec rowind, colptr;            // CSC structure of X1MAT (nk x nk)
    std::vector<arma::uword> slot;        // contribution -> position in the CSC values
    std::vector<std::array<int, 2>> o12;  // contribution -> (o1, o2): adds X1(o1, o2) * coef
    std::vector<double> coef;
  };
  struct Plan133 { const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1; std::vector<Ch133> ch; };

  static Plan133 &get133(const Operator &X, const Operator &Z, Plan &P)
  {
    static Plan133 Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ) return Q;
    Q = Plan133(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max();
    size_t nch3 = ms.GetNumberThreeBodyChannels();
    Q.ch.resize(nch3);
    auto &Z3 = Z.ThreeBody;
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nch3; ich++)
    {
      size_t ch = P.ord3[ich];
      auto &T = ms.GetThreeBodyChannel(ch);
      size_t nk = T.GetNumberKets();
      int twoJ = T.twoJ;
      struct C { arma::uword col, row; int o1, o2; double coef; };
      std::vector<C> cs;
      std::vector<size_t> klist; std::vector<double> rlist;
      for (size_t ibra = 0; ibra < nk; ibra++)
      {
        Ket3 &bra = T.GetKet(ibra);
        int Jij = bra.Jpq;
        std::array<size_t, 3> ijk = {bra.p, bra.q, bra.r};
        for (int iperm = 0; iperm < 3; iperm++)
        {
          std::array<size_t, 3> IJK = ijk;
          Orbit &orep = ms.GetOrbit(ijk[iperm]);
          for (auto a : X.GetOneBodyChannel(orep.l, orep.j2, orep.tz2))
          {
            IJK[iperm] = a;
            if (not Z3.IsKetValid(Jij, twoJ, IJK[0], IJK[1], IJK[2])) continue;
            klist.clear(); rlist.clear();
            size_t chr = Z3.GetKetIndex_withRecoupling(Jij, twoJ, IJK[0], IJK[1], IJK[2], klist, rlist);
            if (chr != ch) continue;
            for (size_t t = 0; t < klist.size(); t++)
              if (klist[t] < nk) cs.push_back({(arma::uword)klist[t], (arma::uword)ibra, (int)ijk[iperm], (int)a, rlist[t]});
          }
        }
      }
      std::vector<size_t> ord(cs.size());
      for (size_t t = 0; t < ord.size(); t++) ord[t] = t;
      std::stable_sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return cs[x].col != cs[y].col ? cs[x].col < cs[y].col : cs[x].row < cs[y].row; });
      Ch133 &H = Q.ch[ch];
      std::vector<arma::uword> rows; std::vector<arma::uword> colptr(nk + 1, 0);
      H.slot.resize(cs.size()); H.o12.resize(cs.size()); H.coef.resize(cs.size());
      for (size_t u = 0; u < ord.size(); u++)
      {
        const C &c = cs[ord[u]];
        if (rows.empty() or u == 0 or cs[ord[u - 1]].col != c.col or cs[ord[u - 1]].row != c.row)
        {
          rows.push_back(c.row);
          colptr[c.col + 1]++;
        }
        H.slot[u] = rows.size() - 1; H.o12[u] = {c.o1, c.o2}; H.coef[u] = c.coef;
      }
      for (size_t k = 0; k < nk; k++) colptr[k + 1] += colptr[k];
      H.rowind = arma::uvec(rows); H.colptr = arma::uvec(colptr);
    }
    return Q;
  }
} // namespace perfplan

void comm133ss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0
      )
  {
    comm133ss(X, Y, Z);
    return;
  }
  auto &Z3 = Z.ThreeBody;
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  if (not(x3 or y3)) return;
  int hXY = (X.IsHermitian() ? 1 : -1) * (Y.IsHermitian() ? 1 : -1);
  perfplan::Plan &P = perfplan::get(X, Z);
  perfplan::Plan133 &Q = perfplan::get133(X, Z, P);
  const perfplan::Cut3 &cutZ = perfplan::cut3(*Z.modelspace);
  static const std::vector<arma::mat> empty;
  const std::vector<arma::mat> &XK = x3 ? perfplan::unpack_k(X, P) : empty;
  const std::vector<arma::mat> &YK = y3 ? perfplan::unpack_k(Y, P) : empty;
  auto &X1 = X.OneBody; auto &Y1 = Y.OneBody;
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < nch3; ich++)
  {
    size_t ch = P.ord3[ich];
    size_t nk = Z.modelspace->GetThreeBodyChannel(ch).GetNumberKets();
    auto &H = Q.ch[ch];
    if (nk == 0 or H.rowind.n_elem == 0) continue;
    if ((x3 and XK[ch].n_elem == 0) or (y3 and YK[ch].n_elem == 0)) continue;
    // Steps 2 + 3: Zk = SX YK - SY XK with SX, SY the embedded 1B operators (shared CSC pattern), one column at a time:
    // contiguous reads of YK, XK and writes within one column of Zk (armadillo's sparse x dense walks dense rows)
    arma::mat Zk(nk, nk, arma::fill::zeros);
    size_t nnz = H.rowind.n_elem;
    std::vector<double> vx(nnz, 0.0), vy(nnz, 0.0);
    if (y3) for (size_t t = 0; t < H.slot.size(); t++) vx[H.slot[t]] += X1(H.o12[t][0], H.o12[t][1]) * H.coef[t];
    if (x3) for (size_t t = 0; t < H.slot.size(); t++) vy[H.slot[t]] += Y1(H.o12[t][0], H.o12[t][1]) * H.coef[t];
    const arma::uword *cp = H.colptr.memptr(), *ri = H.rowind.memptr();
    for (size_t j = 0; j < nk; j++)
    {
      double *zc = Zk.colptr(j);
      const double *yc = y3 ? YK[ch].colptr(j) : nullptr, *xc = x3 ? XK[ch].colptr(j) : nullptr;
      for (size_t c = 0; c < nk; c++)
      {
        double dy = y3 ? yc[c] : 0.0, dx = x3 ? xc[c] : 0.0;
        if (dy == 0 and dx == 0) continue;
        for (arma::uword t = cp[c]; t < cp[c + 1]; t++) zc[ri[t]] += vx[t] * dy - vy[t] * dx;
      }
    }
    arma::mat Zt = Zk - hXY * Zk.t();   // Z3 - hX hY Z3^T
    perfplan::mask3(cutZ, ch, ch, Zt);   // dE3max / occnat3 cut on the output
    Z3.AddChannelMatrix_pn(ch, Zt, -1);   // = AddToME_pn_ch for iket >= ibra
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
