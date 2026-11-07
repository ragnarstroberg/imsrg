//=====================================================================================================================
//  comm331ss with precompiled bookkeeping (scalar X, Y, default cuts). Notation: see the reader's guide in PerfPlan.inc.
//=====================================================================================================================
//  Formula (IMSRG3Commutators.cc):
//    Z_ij = 1/12 sum_abcde (na nb n'c n'd n'e + n'a n'b nc nd ne) (2J+1)/(2ji+1) (X^J_{abi,cde} Y^J_{cde,abj} - Y <-> X)
//  (summed over the couplings J1 of ab, J2 of cde and J). In the plan, cde is a stored ket and (ab)i a product state:
//    X_{(ab)i,cde} = (R X)(P(ab, i), cde)  and  X_{cde,(ab)j} = hX (R X)(P(ab, j), cde)   (half-unpacked R X, cached).
//  The occupation factor is a sum of two (ab) x (cde) products, so for each pair ab one small GEMM over the 3B kets
//  gives all (i, j):
//    Step 1  rows S = {P(ab, i)} of R X and R Y for this pair, the cde columns weighted by w (occupations, 1/12 s_ab s_cde)
//    Step 2  M = (RX_S w) RY_S^T  (and Y <-> X), added to Z_ij with (2J+1)/(2ji+1).
//  Included after PerfPlan.inc (inside namespace Commutator).
void comm331ss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0)
  {
    comm331ss(X, Y, Z);
    return;
  }
  if (X.GetParticleRank() < 3 or Y.GetParticleRank() < 3) return;
  if (not(X.ThreeBody.IsAllocated() and Y.ThreeBody.IsAllocated())) return;
  int herm = Z.IsAntiHermitian() ? -1 : 1;
  int hX = X.IsHermitian() ? 1 : -1, hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(X, Z);
  const std::vector<arma::mat> &RX = perfplan::unpack_rx(X, P);
  const std::vector<arma::mat> &RY = perfplan::unpack_rx(Y, P);
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  int norb = P.norb;
  int nthr = omp_get_max_threads();
  std::vector<arma::mat> Zt(nthr, arma::mat(norb, norb, arma::fill::zeros));
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < nch3; ich++)
  {
    size_t ch = P.ord3[ich];
    auto &T = ms.GetThreeBodyChannel(ch);
    size_t nP = P.Pst[ch].size(), nk = T.GetNumberKets();
    if (nP == 0 or nk == 0 or RX[ch].n_elem == 0 or RY[ch].n_elem == 0) continue;
    arma::mat &Zl = Zt[omp_get_thread_num()];
    // ket weights: s_cde (1-nc)(1-nd)(1-ne) and s_cde nc nd ne
    arma::rowvec g1(nk), g2(nk);
    for (size_t k = 0; k < nk; k++)
    {
      Ket3 &ket = T.GetKet(k);
      double nc = ket.op->occ, nd = ket.oq->occ, ne = ket.oR->occ;
      double s = 6;
      if (ket.p == ket.q and ket.q == ket.r) s = 1;
      else if (ket.p == ket.q or ket.q == ket.r) s = 3;
      g1(k) = s * (1 - nc) * (1 - nd) * (1 - ne);
      g2(k) = s * nc * nd * ne;
    }
    // P rows grouped by pair g
    std::vector<std::array<int, 2>> gp(nP);
    for (size_t ip = 0; ip < nP; ip++) gp[ip] = {P.Pst[ch][ip][0], (int)ip};
    std::sort(gp.begin(), gp.end());
    for (size_t s0 = 0; s0 < nP;)
    {
      size_t s1 = s0;
      while (s1 < nP and gp[s1][0] == gp[s0][0]) s1++;
      int g = gp[s0][0];
      int a = P.g2[g][0], b = P.g2[g][1];
      double na = ms.GetOrbit(a).occ, nb = ms.GetOrbit(b).occ;
      double f1 = na * nb, f2 = (1 - na) * (1 - nb);
      size_t ns = s1 - s0;
      if (std::abs(f1) + std::abs(f2) >= 1e-8)
      {
        arma::uvec S(ns); std::vector<int> ov(ns);
        for (size_t u = 0; u < ns; u++) { S(u) = gp[s0 + u][1]; ov[u] = P.Pst[ch][S(u)][1]; }
        arma::rowvec w = f1 * g1 + f2 * g2;
        arma::mat A = RX[ch].rows(S);
        A.each_row() %= w;
        arma::mat M = A * RY[ch].rows(S).t();
        double sab = (a == b) ? 1.0 : 2.0;
        for (size_t u = 0; u < ns; u++)
        {
          Orbit &oi = ms.GetOrbit(ov[u]);
          double fac = sab / 12.0 * (T.twoJ + 1.0) / (oi.j2 + 1);
          for (size_t v = 0; v < ns; v++)
          {
            int i = ov[u], j = ov[v];
            if (j > i) continue;
            Orbit &oj = ms.GetOrbit(j);
            if (oj.l != oi.l or oj.j2 != oi.j2 or oj.tz2 != oi.tz2) continue;
            double zij = fac * (hY * M(u, v) - hX * M(v, u));
            Zl(i, j) += zij;
            if (i != j) Zl(j, i) += herm * zij;
          }
        }
      }
      s0 = s1;
    }
  }
  for (auto &m : Zt) Z.OneBody += m;
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
