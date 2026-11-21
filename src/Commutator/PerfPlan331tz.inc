//=====================================================================================================================
//  comm331ss, scalar X and Y changing isospin and/or parity. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
//  Formula as comm331ss (PerfPlan331.inc). With X diagonal in the 3B channels and Y connecting channel cb to ck,
//  for i in channel cb (as P(ab, i)) and j in ck (as P(ab, j)):
//    Z_ij = f_ab,i sum_cde w(cde) [ hY (R X)_cb(P(ab,i), cde) (R Y)_{ck<-cb}(P(ab,j), cde)
//                                  - hX (R Y)_{cb<-ck}(P(ab,i), cde) (R X)_ck(P(ab,j), cde) ]
//  with (R Y)_{x<-y} = R(x) Y(x, y) from the off-diagonal unpack (both orientations), w the occupation weights of the
//  intermediate kets cde (as in comm331ss_plan) and f_ab,i = s_ab / 12 (2J+1)/(2j_i+1). For cb = ck it reduces to the
//  scalar kernel. One small GEMM per (pair ab, ordered channel pair).
//  Included after PerfPlan333tz.inc (inside namespace Commutator).
void comm331ss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced())
  {
    comm331ss(X, Y, Z);
    return;
  }
  if (X.GetParticleRank() < 3 or Y.GetParticleRank() < 3) return;
  if (not(X.ThreeBody.IsAllocated() and Y.ThreeBody.IsAllocated())) return;
  int hX = X.IsHermitian() ? 1 : -1, hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(X, Z);
  const std::vector<arma::mat> &RX = perfplan::unpack_rx(X, P);
  const perfplan::WOffCache &WY = perfplan::unpack_off(Y, P);
  const std::vector<arma::mat> WYrk = perfplan::unpack_off_rk(WY, P);   // R(cb) Y, held while this term runs
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  int norb = P.norb;
  // occupation weights of the stored kets (g1: all particles, g2: all holes), with the ket symmetry factor
  std::vector<arma::rowvec> G1(nch3), G2(nch3);
  for (size_t ch = 0; ch < nch3; ch++)
  {
    auto &T = ms.GetThreeBodyChannel(ch);
    size_t nk = T.GetNumberKets();
    G1[ch].set_size(nk); G2[ch].set_size(nk);
    for (size_t k = 0; k < nk; k++)
    {
      Ket3 &ket = T.GetKet(k);
      double nc = ket.op->occ, nd = ket.oq->occ, ne = ket.oR->occ;
      double s = 6;
      if (ket.p == ket.q and ket.q == ket.r) s = 1;
      else if (ket.p == ket.q or ket.q == ket.r) s = 3;
      G1[ch](k) = s * (1 - nc) * (1 - nd) * (1 - ne);
      G2[ch](k) = s * nc * nd * ne;
    }
  }
  // ordered channel pairs (cb, ck) with Y(cb, ck) and Y(ck, cb) available
  std::vector<std::array<size_t, 2>> pairs;
  for (size_t cb = 0; cb < nch3; cb++)
    for (size_t ck = 0; ck < nch3; ck++)
      if (cb != ck and WY.idx[cb * WY.nch + ck] >= 0 and WY.idx[ck * WY.nch + cb] >= 0) pairs.push_back({cb, ck});
  int nthr = omp_get_max_threads();
  std::vector<arma::mat> Zt(nthr, arma::mat(norb, norb, arma::fill::zeros));
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ip = 0; ip < pairs.size(); ip++)
  {
    size_t cb = pairs[ip][0], ck = pairs[ip][1];
    auto &Tb = ms.GetThreeBodyChannel(cb);
    if (RX[cb].n_elem == 0 or RX[ck].n_elem == 0) continue;
    const arma::mat &RYkb = WYrk[WY.idx[ck * WY.nch + cb]];   // R(ck) Y(ck, cb): rows P in ck, cols stored kets of cb
    const arma::mat &RYbk = WYrk[WY.idx[cb * WY.nch + ck]];   // R(cb) Y(cb, ck): rows P in cb, cols stored kets of ck
    if (RYkb.n_elem == 0 or RYbk.n_elem == 0) continue;
    arma::mat &Zl = Zt[omp_get_thread_num()];
    // P rows of both channels grouped by pair g
    std::map<int, std::vector<arma::uword>> rb, rk;
    for (size_t p = 0; p < P.Pst[cb].size(); p++) rb[P.Pst[cb][p][0]].push_back(p);
    for (size_t p = 0; p < P.Pst[ck].size(); p++) rk[P.Pst[ck][p][0]].push_back(p);
    for (auto &gb : rb)
    {
      auto it = rk.find(gb.first);
      if (it == rk.end()) continue;
      int g = gb.first, a = P.g2[g][0], b = P.g2[g][1];
      double na = ms.GetOrbit(a).occ, nb = ms.GetOrbit(b).occ;
      double f1 = na * nb, f2 = (1 - na) * (1 - nb);
      if (std::abs(f1) + std::abs(f2) < 1e-8) continue;
      arma::uvec Sb(gb.second), Sk(it->second);
      arma::rowvec wb = f1 * G1[cb] + f2 * G2[cb], wk = f1 * G1[ck] + f2 * G2[ck];
      arma::mat A1 = RX[cb].rows(Sb); A1.each_row() %= wb;
      arma::mat M1 = A1 * RYkb.rows(Sk).t();                     // sum over cde in cb
      arma::mat A2 = RYbk.rows(Sb); A2.each_row() %= wk;
      arma::mat M2 = A2 * RX[ck].rows(Sk).t();                   // sum over cde in ck
      double sab = (a == b) ? 1.0 : 2.0;
      for (size_t u = 0; u < Sb.n_elem; u++)
      {
        int i = P.Pst[cb][Sb(u)][1];
        double fac = sab / 12.0 * (Tb.twoJ + 1.0) / (ms.GetOrbit(i).j2 + 1);
        for (size_t v = 0; v < Sk.n_elem; v++)
        {
          int j = P.Pst[ck][Sk(v)][1];
          if (ms.GetOrbit(j).j2 != ms.GetOrbit(i).j2) continue;   // J = 0 operator: ji = jj
          Zl(i, j) += fac * (hY * M1(u, v) - hX * M2(u, v));
        }
      }
    }
  }
  for (auto &m : Zt) Z.OneBody += m;
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
