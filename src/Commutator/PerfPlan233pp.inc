//=====================================================================================================================
//  comm233_pp_hhss with precompiled bookkeeping (scalar X, Y, default dE3max / occnat3 cuts). Notation: PerfPlan.inc.
//=====================================================================================================================
//  Formula (IMSRG3Commutators.cc):
//    Z_ijklmn = 1/2 sum_{ab} (n'a n'b - na nb) [ P(ij/k) X_{ij,ab} Y_{abk,lmn} - P(lm/n) Y_{ijk,abn} X_{ab,lm} ]  - (X <-> Y)
//  The original builds, per 3B channel, X2MAT = A^T D_X R (the 2B operator embedded in the 3B space, including P(ik),
//  P(jk)) element by element and the dense X3MAT from GetME_pn_ch, then does two dense n_k x n_k x n_k GEMMs. Here
//    Z3 = A^T (D_X (R Y3) - D_Y (R X3)),     Z3 -= hX hY Z3^T,
//  with
//    R Y3, R X3   from the unpack cache (the 3B operator with its bra in the product basis);
//    D            the 2B operator acting on the pair of a product state: block diagonal in (third orbit c, pair J),
//                 D((ij)J c, (ab)J c) = X^J_{ij,ab} * 0.5 symm (1 - na - nb), sqrt2 for identical pair orbits
//                 (GetTBME is un-normalised), as in the original;
//    A^T          the antisymmetriser P(ij/k) of the bra (sparse, see PerfPlan.inc).
//  Step 1  per (c, J) group of P rows: Y-part block = D_X block * (R Y3)(group rows, :), dense and small.
//  Step 2  Z3 = A^T (sum of the blocks), Z3 - hX hY Z3^T into the stored block.
//  Included after PerfPlan.inc (inside namespace Commutator).
namespace perfplan
{
  struct Grp233pp { int ch2; arma::uvec rows, loc; arma::rowvec scale; arma::vec rnorm; };   // rnorm: sqrt2 for a==b (GetTBME is unnormalised)
  struct Plan233pp
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1;
    std::vector<std::vector<Grp233pp>> grp;   // ch3 -> groups of P rows with the same (c, J_pair)
  };

  static Plan233pp &get233pp(const Operator &Z, Plan &P)
  {
    static Plan233pp Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ) return Q;
    Q = Plan233pp(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max();
    int norb = P.norb;
    size_t nch3 = ms.GetNumberThreeBodyChannels();
    Q.grp.resize(nch3);
    for (size_t ch = 0; ch < nch3; ch++)
    {
      std::map<std::pair<int, int>, std::vector<int>> byCJ;   // (c, J) -> P rows
      for (size_t ip = 0; ip < P.Pst[ch].size(); ip++)
      {
        int g = P.Pst[ch][ip][0], c = P.Pst[ch][ip][1];
        byCJ[{c, P.g2[g][2]}].push_back(ip);
      }
      for (auto &kv : byCJ)
      {
        Grp233pp G;
        auto &rows = kv.second;
        int g0 = P.Pst[ch][rows[0]][0];
        Orbit &op0 = ms.GetOrbit(P.g2[g0][0]); Orbit &oq0 = ms.GetOrbit(P.g2[g0][1]);
        G.ch2 = ms.GetTwoBodyChannelIndex(P.g2[g0][2], (op0.l + oq0.l) % 2, (op0.tz2 + oq0.tz2) / 2);
        auto &tbc = ms.GetTwoBodyChannel(G.ch2);
        G.rows.set_size(rows.size()); G.loc.set_size(rows.size()); G.scale.set_size(rows.size()); G.rnorm.set_size(rows.size());
        for (size_t r = 0; r < rows.size(); r++)
        {
          int g = P.Pst[ch][rows[r]][0];
          int a = P.g2[g][0], b = P.g2[g][1];
          G.rows(r) = rows[r];
          G.loc(r) = tbc.GetLocalIndex(a, b);
          double occf = 1 - ms.GetOrbit(a).occ - ms.GetOrbit(b).occ;
          G.rnorm(r) = (a == b) ? PhysConst::SQRT2 : 1.0;
          G.scale(r) = (std::abs(occf) < 1e-6) ? 0.0 : 0.5 * ((a == b) ? 1.0 : 2.0) * occf * G.rnorm(r);
        }
        Q.grp[ch].push_back(std::move(G));
      }
    }
    (void)norb;
    return Q;
  }
} // namespace perfplan

void comm233_pp_hhss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0
      )
  {
    comm233_pp_hhss(X, Y, Z);
    return;
  }
  auto &Z3 = Z.ThreeBody;
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  if (not(x3 or y3)) return;
  int hXY = (X.IsHermitian() ? 1 : -1) * (Y.IsHermitian() ? 1 : -1);
  perfplan::Plan &P = perfplan::get(X, Z);
  perfplan::Plan233pp &Q = perfplan::get233pp(Z, P);
  const perfplan::Cut3 &cutZ = perfplan::cut3(*Z.modelspace);
  static const std::vector<arma::mat> empty;
  const std::vector<arma::mat> &RX = x3 ? perfplan::unpack_rx(X, P) : empty;
  const std::vector<arma::mat> &RY = y3 ? perfplan::unpack_rx(Y, P) : empty;
  auto &X2 = X.TwoBody; auto &Y2 = Y.TwoBody;
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < nch3; ich++)
  {
    size_t ch = P.ord3[ich];
    size_t nP = P.Pst[ch].size(), nk = Z.modelspace->GetThreeBodyChannel(ch).GetNumberKets();
    if (nP == 0 or nk == 0) continue;
    arma::mat G(nP, nk, arma::fill::zeros);
    for (auto &g : Q.grp[ch])
    {
      arma::mat acc(g.rows.n_elem, nk, arma::fill::zeros);   // Step 1: D_X (R Y3) - D_Y (R X3) on this (c, J) group
      if (y3)
      {
        arma::mat Dx = X2.GetMatrix(g.ch2, g.ch2).submat(g.loc, g.loc);
        Dx.each_row() %= g.scale; Dx.each_col() %= g.rnorm;
        acc += Dx * RY[ch].rows(g.rows);
      }
      if (x3)
      {
        arma::mat Dy = Y2.GetMatrix(g.ch2, g.ch2).submat(g.loc, g.loc);
        Dy.each_row() %= g.scale; Dy.each_col() %= g.rnorm;
        acc -= Dy * RX[ch].rows(g.rows);
      }
      G.rows(g.rows) = acc;
    }
    arma::mat M = perfplan::AtM(P.A[ch], G);   // Step 2: A^T (...), then M - hX hY M^T
    arma::mat Zk = M - hXY * M.t();
    perfplan::mask3(cutZ, ch, ch, Zk);   // dE3max / occnat3 cut on the output
    Z3.AddChannelMatrix_pn(ch, Zk, +1);   // = AddToME_pn_ch for iket <= ibra
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
