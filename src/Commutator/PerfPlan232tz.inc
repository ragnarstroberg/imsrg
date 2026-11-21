//=====================================================================================================================
//  comm232ss, scalar X and isospin-changing Y. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
// comm232ss for a scalar X and a channel-changing Y (J rank 0, parity 0, T rank > 0), precompiled bookkeeping.
// Generalises comm232ss_plan: the pph states |klj`; j_i> are grouped in channels K = (j2, parity, tz2) with tz2 in
// {-3,-1,1,3} (Y shifts the charge, so K need not be the channel of a real orbit). For an orbit channel o of i and
// K' = (j2_o, par_o, tz2_o + 2 Delta):
//   ZMAT(o, K')(i, klj) = X2(i, abc in K(o)) Y3(abc in K(o), klj in K') - Y2(i, abc in K') X3(abc in K', klj in K')
// and the same 4-term assembly into every stored block of Z2. Included after PerfPlan232.inc and PerfPlan333tz.inc.
namespace perfplan
{
  struct KC232
  {
    int j2, par, tz2;
    std::vector<std::array<int, 2>> klj;   // (g_kl, j)
    std::vector<int> abc;                  // indices into klj with nonzero occupation weight
    std::vector<double> occ;
  };
  struct Plan232tz
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1;
    std::vector<KC232> kc;
    std::vector<int> kcTab;                    // [(j2 * 2 + par) * 4 + (tz2 + 3) / 2] -> kc or -1
    std::vector<int> kcOfObc;                  // orbit channel -> kc
    std::vector<std::vector<int>> colIdx;      // kc -> [g * norb + j] -> klj index or -1
    int key(int j2, int par, int tz2) const { return (j2 * 2 + par) * 4 + (tz2 + 3) / 2; }
  };
  static Plan232tz &get232tz(const Operator &Z, Plan &P, Plan232 &Q2)
  {
    static Plan232tz Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ) return Q;
    Q = Plan232tz(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max();
    int norb = P.norb, ng = P.ng3;   // klj columns index product states P(kl, j): pairs of the 3B space only
    Q.kcTab.assign((P.j2max + 1) * 2 * 4, -1);
    for (int j2 = 1; j2 <= P.j2max; j2 += 2)
      for (int par = 0; par <= 1; par++)
        for (int tz2 = -3; tz2 <= 3; tz2 += 2)
        {
          KC232 C; C.j2 = j2; C.par = par; C.tz2 = tz2;
          std::vector<int> ci((size_t)ng * norb, -1);
          for (int g = 0; g < ng; g++)
          {
            int k = P.g2[g][0], l = P.g2[g][1], Jkl = P.g2[g][2];
            Orbit &ok = ms.GetOrbit(k); Orbit &ol = ms.GetOrbit(l);
            for (int j = 0; j < norb; j++)
            {
              Orbit &oj = ms.GetOrbit(j);
              if (2 * oj.n + oj.l > ms.GetEMax3Body()) continue;   // P(kl, j) needs j in the 3B space
              if ((ok.l + ol.l + oj.l) % 2 != par) continue;
              if (ok.tz2 + ol.tz2 - oj.tz2 != tz2) continue;
              if (std::abs(2 * Jkl - oj.j2) > j2 or 2 * Jkl + oj.j2 < j2) continue;
              ci[(size_t)g * norb + j] = C.klj.size();
              double na = ok.occ, nb = ol.occ, nc = oj.occ;
              double w = na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc;
              if (k == l) w *= 0.5;
              if (std::abs(na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc) > 1e-8) { C.abc.push_back(C.klj.size()); C.occ.push_back(w); }
              C.klj.push_back({g, j});
            }
          }
          if (C.klj.empty()) continue;
          Q.kcTab[Q.key(j2, par, tz2)] = Q.kc.size();
          Q.kc.push_back(std::move(C)); Q.colIdx.push_back(std::move(ci));
        }
    Q.kcOfObc.resize(Q2.obc.size());
    for (size_t o = 0; o < Q2.obc.size(); o++) Q.kcOfObc[o] = Q.kcTab[Q.key(Q2.obc[o].j2, Q2.obc[o].par, Q2.obc[o].tz2)];
    return Q;
  }
} // namespace perfplan

void comm232ss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or X.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0 or (Y.GetParity() == 0 and Y.GetTRank() == 0)
      or Y.IsReduced() or Z.IsReduced() or imsrg3_valence_2b)
  {
    comm232ss(X, Y, Z);
    return;
  }
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  if (not(x3 or y3)) return;
  int hX = X.IsHermitian() ? 1 : -1, hY = Y.IsHermitian() ? 1 : -1, TY = Y.GetTRank(), pY = Y.GetParity();
  perfplan::Plan &P = perfplan::get(Z, Z);
  perfplan::Plan232 &Q2 = perfplan::get232(Z, P);
  perfplan::Plan232tz &Q = perfplan::get232tz(Z, P, Q2);
  int norb = P.norb, nJ = P.J2max + 1;
  static const std::vector<arma::mat> empty;
  static const std::vector<perfplan::wmat> emptyW;
  static perfplan::WOffCache emptyoff;
  const std::vector<perfplan::wmat> &WX = x3 ? perfplan::unpack(X, P) : emptyW;
  const perfplan::WOffCache &WY = y3 ? perfplan::unpack_off(Y, P) : emptyoff;
  auto six = [&](int a, int b, int c, int d, int e, int f) { return Q2.six6(a, b, c, d, e, f); };

  // Out(ia, ik) += sum_J' (2J'+1)/sqrt(2Jkl+1) {jj ji Jkl; jc J' Jab} W(P(gab, j), P(gkl, c)),  rows abc of A, cols klj of B
  auto fill3 = [&](arma::mat &Out, const perfplan::KC232 &A, const perfplan::KC232 &B, int j2i, auto getW)
  {
    size_t da = A.abc.size(), dk = B.klj.size();
    Out.zeros(da, dk);
    if (da == 0 or dk == 0) return;
    std::vector<std::vector<int>> rowsByC(norb), colsByJ(norb);
    for (size_t ia = 0; ia < da; ia++) rowsByC[A.klj[A.abc[ia]][1]].push_back(ia);
    for (size_t ik = 0; ik < dk; ik++) colsByJ[B.klj[ik][1]].push_back(ik);
    std::vector<int> r_ia, r_p, r_J, q_ik, q_p, q_J;
    std::vector<double> wt(nJ * nJ);
    for (int c = 0; c < norb; c++)
    {
      auto &Rc = rowsByC[c];
      if (Rc.empty()) continue;
      Orbit &oc = ms.GetOrbit(c);
      int g0 = A.klj[A.abc[Rc[0]]][0];
      Orbit &oa0 = ms.GetOrbit(P.g2[g0][0]); Orbit &ob0 = ms.GetOrbit(P.g2[g0][1]);
      for (int j = 0; j < norb; j++)
      {
        auto &Qj = colsByJ[j];
        if (Qj.empty()) continue;
        Orbit &oj = ms.GetOrbit(j);
        int h0 = B.klj[Qj[0]][0];
        Orbit &ok0 = ms.GetOrbit(P.g2[h0][0]); Orbit &ol0 = ms.GetOrbit(P.g2[h0][1]);
        int par_b = (oa0.l + ob0.l + oj.l) % 2, tz_b = oa0.tz2 + ob0.tz2 + oj.tz2;
        int par_k = (ok0.l + ol0.l + oc.l) % 2, tz_k = ok0.tz2 + ol0.tz2 + oc.tz2;
        for (int tJp = 1; tJp <= P.tJmax; tJp += 2)
        {
          int chb = P.ch3tab[P.i3(tJp, par_b, tz_b)], chk = P.ch3tab[P.i3(tJp, par_k, tz_k)];
          if (chb < 0 or chk < 0) continue;
          const perfplan::wmat *W = getW(chb, chk);
          if (W == nullptr or W->n_elem == 0) continue;
          r_ia.clear(); r_p.clear(); r_J.clear(); q_ik.clear(); q_p.clear(); q_J.clear();
          for (int ia : Rc)
          {
            int gab = A.klj[A.abc[ia]][0], Jab = P.g2[gab][2];
            if (std::abs(2 * Jab - oj.j2) > tJp or 2 * Jab + oj.j2 < tJp) continue;
            int p1 = P.pidx(chb, gab, j);
            if (p1 >= 0) { r_ia.push_back(ia); r_p.push_back(p1); r_J.push_back(Jab); }
          }
          if (r_ia.empty()) continue;
          for (int ik : Qj)
          {
            int gkl = B.klj[ik][0], Jkl = P.g2[gkl][2];
            if (std::abs(2 * Jkl - oc.j2) > tJp or 2 * Jkl + oc.j2 < tJp) continue;
            int p2 = P.pidx(chk, gkl, c);
            if (p2 >= 0) { q_ik.push_back(ik); q_p.push_back(p2); q_J.push_back(Jkl); }
          }
          if (q_ik.empty()) continue;
          for (int Jab = 0; Jab < nJ; Jab++)
            for (int Jkl = 0; Jkl < nJ; Jkl++)
              wt[Jab * nJ + Jkl] = (tJp + 1.) / sqrt(2 * Jkl + 1.) * six(oj.j2, j2i, 2 * Jkl, oc.j2, tJp, 2 * Jab);
          for (size_t q = 0; q < q_ik.size(); q++)
          {
            const double *wq = &wt[q_J[q]];
            const auto *wc = W->colptr(q_p[q]); double *oc_ = Out.colptr(q_ik[q]);
            for (size_t r = 0; r < r_ia.size(); r++) oc_[r_ia[r]] += wq[r_J[r] * nJ] * wc[r_p[r]];
          }
        }
      }
    }
  };
  // f2 X^{Jab}_{ci,ab} for i in orbit channel o, abc of channel K
  auto fill2 = [&](arma::mat &Out, const Operator &O, const perfplan::OBC232 &Co, const perfplan::KC232 &K)
  {
    size_t di = Co.orbs.size(), da = K.abc.size();
    Out.zeros(di, da);
    for (size_t ia = 0; ia < da; ia++)
    {
      int gab = K.klj[K.abc[ia]][0], c = K.klj[K.abc[ia]][1];
      int a = P.g2[gab][0], b = P.g2[gab][1], Jab = P.g2[gab][2];
      double f2 = -sqrt(2 * Jab + 1.) * K.occ[ia];
      for (size_t ii = 0; ii < di; ii++) Out(ii, ia) = f2 * O.TwoBody.GetTBME_J(Jab, c, Co.orbs[ii], a, b);
    }
  };
  auto getWX = [&](int chb, int chk) -> const perfplan::wmat * { return chb == chk ? &WX[chb] : nullptr; };
  auto getWY = [&](int chb, int chk) -> const perfplan::wmat * { int k = WY.idx[(size_t)chb * WY.nch + chk]; return k < 0 ? nullptr : &WY.mats[k]; };

  // tasks (o, K')
  size_t nobc = Q2.obc.size(), nkc = Q.kc.size();
  std::vector<std::array<int, 2>> tasks;
  std::vector<int> taskOf(nobc * nkc, -1);
  for (size_t o = 0; o < nobc; o++)
    for (int id = 0; id < (TY == 0 ? 1 : 2); id++)   // Delta = +TY, -TY; parity-changing only (TY = 0): one partner
    {
      int d = id == 0 ? TY : -TY;
      auto &Co = Q2.obc[o];
      int tz2 = Co.tz2 + 2 * d;
      if (tz2 < -3 or tz2 > 3) continue;
      int k = Q.kcTab[Q.key(Co.j2, (Co.par + pY) % 2, tz2)];
      if (k < 0 or Q.kcOfObc[o] < 0) continue;
      taskOf[o * nkc + k] = tasks.size(); tasks.push_back({(int)o, k});
    }
  std::vector<arma::mat> ZMAT(tasks.size());
  double t = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t it = 0; it < tasks.size(); it++)
  {
    int o = tasks[it][0], k = tasks[it][1];
    auto &Co = Q2.obc[o];
    auto &A = Q.kc[Q.kcOfObc[o]];
    auto &B = Q.kc[k];
    ZMAT[it].zeros(Co.orbs.size(), B.klj.size());
    if (y3 and A.abc.size() > 0)
    {
      arma::mat X2, Y3m;
      fill2(X2, X, Co, A);
      fill3(Y3m, A, B, Co.j2, getWY);
      ZMAT[it] += X2 * Y3m;
    }
    if (x3 and B.abc.size() > 0)
    {
      arma::mat Y2, X3m;
      fill2(Y2, Y, Co, B);
      fill3(X3m, B, B, Co.j2, getWX);
      ZMAT[it] -= Y2 * X3m;
    }
  }
  Z.profiler.timer["_comm232_plantz_fill_gemm"] += omp_get_wtime() - t; t = omp_get_wtime();

  // assembly into every stored block of Z2 (all bra/ket pairs for off-diagonal blocks)
  auto &Z2 = Z.TwoBody;
  std::vector<std::array<size_t, 3>> bk;
  for (auto &it : Z2.MatEl)
    for (size_t ibra = 0; ibra < ms.GetTwoBodyChannel(it.first[0]).GetNumberKets(); ibra++) bk.push_back({it.first[0], it.first[1], ibra});
  auto get = [&](int r, int g, int s) -> double
  {
    int o = Q2.obcOf[r];
    Orbit &or_ = ms.GetOrbit(r); Orbit &os = ms.GetOrbit(s);
    Orbit &og1 = ms.GetOrbit(P.g2[g][0]); Orbit &og2 = ms.GetOrbit(P.g2[g][1]);
    int tz2 = og1.tz2 + og2.tz2 - os.tz2;
    if (tz2 < -3 or tz2 > 3) return 0.0;
    int k = Q.kcTab[Q.key(or_.j2, (og1.l + og2.l + os.l) % 2, tz2)];   // channel of the (g s) states
    if (k < 0) return 0.0;
    int itk = taskOf[(size_t)o * nkc + k];
    if (itk < 0) return 0.0;
    int c = (g >= 0 and g < P.ng3) ? Q.colIdx[k][(size_t)g * norb + s] : -1;
    return c < 0 ? 0.0 : ZMAT[itk](Q2.posOf[r], c);
  };
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < bk.size(); ib++)
  {
    size_t chb = bk[ib][0], chk = bk[ib][1], ibra = bk[ib][2];
    TwoBodyChannel &tbcb = ms.GetTwoBodyChannel(chb); TwoBodyChannel &tbck = ms.GetTwoBodyChannel(chk);
    if (tbcb.J != tbck.J) continue;
    int J = tbcb.J;
    Ket &bra = tbcb.GetKet(ibra);
    int i = bra.p, j = bra.q;
    Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
    int gij = P.g2idx[((size_t)i * norb + j) * (P.J2max + 1) + J];
    for (size_t iket = (chb == chk ? ibra : 0); iket < tbck.GetNumberKets(); iket++)
    {
      Ket &ket = tbck.GetKet(iket);
      int k = ket.p, l = ket.q;
      Orbit &ok = ms.GetOrbit(k); Orbit &ol = ms.GetOrbit(l);
      int gkl = P.g2idx[((size_t)k * norb + l) * (P.J2max + 1) + J];
      int phase_ij = ms.phase((oi.j2 + oj.j2) / 2 - J), phase_kl = ms.phase((ok.j2 + ol.j2) / 2 - J);
      double zijkl = phase_ij * get(i, gkl, j) - get(j, gkl, i) + hX * hY * (get(l, gij, k) - phase_kl * get(k, gij, l));
      zijkl /= sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
      Z2.AddToTBME(chb, chk, ibra, iket, zijkl);
    }
  }
  Z.profiler.timer["_comm232_plantz_back"] += omp_get_wtime() - t;
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
