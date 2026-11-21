//=====================================================================================================================
//  comm232ss with precompiled bookkeeping (scalar X, Y). Notation: see the reader's guide in PerfPlan.inc.
//=====================================================================================================================
//  Formula (IMSRG3Commutators.cc):
//    Z_ijkl = -1/2 sum_abc (na nb n'c + n'a n'b nc) ( (1 - Pij) X_{ic,ab} Y_{abj,klc} - (1 - Pkl) Y_{ijc,abl} X_{ab,kc} )
//  Follows the structure of comm232ss_srs_optimized: per one-body channel of i, with intermediate states abc = (gab, c)
//  (nonzero occupation factor occ_abc) and pph states klj = (gkl, j),
//    ZMAT(i, klj) = X2(i, abc) Y3(abc, klj) - Y2(i, abc) X3(abc, klj)
//    Step 1  X2(i, abc) = -sqrt(2Jab+1) occ_abc X^{Jab}_{ci,ab}                 (un-normalised 2B matrix element)
//    Step 2  Y3(abc, klj) = sum_J' (2J'+1)/sqrt(2Jkl+1) {jj ji Jkl; jc J' Jab} W_Y^{J'}(P(gab, j), P(gkl, c))
//            gathered from the unpacked W in blocks of fixed (c, j), with the 6j in a small (Jab, Jkl) table
//    Step 3  ZMAT = X2 Y3 - Y2 X3 (GEMM), then the same 4-term assembly into Z2 as the original ((1 - Pij), (1 - Pkl)).
//  Included after PerfPlan332.inc (inside namespace Commutator).
namespace perfplan
{
  struct OBC232
  {
    int j2, par, tz2;
    std::vector<int> orbs;                         // i in the channel
    std::vector<std::array<int, 2>> klj;            // (g_kl, j)
    std::vector<int> abc;                           // indices into klj with nonzero occupation weight
    std::vector<double> occ;                        // occupation weight (incl. 1/2 for a == b)
  };
  struct Plan232
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1;
    std::vector<OBC232> obc;
    std::vector<int> obcOf, posOf;                  // orbit -> one-body channel, position in it
    std::vector<std::vector<int>> colIdx;           // obc -> [g_kl * norb + j] -> klj index or -1
    std::vector<double> six; int d6[6];
    size_t i6(int a, int b, int c, int d, int e, int f) const
    { return ((((((size_t)a * d6[1] + b) * d6[2] + c) * d6[3] + d) * d6[4] + e) * d6[5] + f); }
    // bounds-checked read: arguments outside the table (orbits or pairs outside the 3B space) give 0
    double six6(int a, int b, int c, int d, int e, int f) const
    {
      if (a < 0 or b < 0 or c < 0 or d < 0 or e < 0 or f < 0 or a >= d6[0] or b >= d6[1] or c >= d6[2] or d >= d6[3] or e >= d6[4] or f >= d6[5]) return 0.0;
      return six[i6(a, b, c, d, e, f)];
    }
  };
  static Plan232 &get232(const Operator &Z, Plan &P)
  {
    static Plan232 Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ) return Q;
    Q = Plan232(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max();
    int norb = P.norb, ng = P.ng3;   // klj columns index product states P(kl, j): pairs of the 3B space only
    Q.obcOf.assign(norb, -1); Q.posOf.assign(norb, -1);
    for (auto &it : ms.OneBodyChannels)
    {
      OBC232 C; C.j2 = it.first[1]; C.par = it.first[0] % 2; C.tz2 = it.first[2];
      for (auto i : it.second) { Q.obcOf[i] = Q.obc.size(); Q.posOf[i] = C.orbs.size(); C.orbs.push_back(i); }
      Q.obc.push_back(C);
    }
    Q.colIdx.resize(Q.obc.size());
    for (size_t o = 0; o < Q.obc.size(); o++)
    {
      auto &C = Q.obc[o];
      Q.colIdx[o].assign((size_t)ng * norb, -1);
      for (int g = 0; g < ng; g++)
      {
        int k = P.g2[g][0], l = P.g2[g][1], Jkl = P.g2[g][2];
        Orbit &ok = ms.GetOrbit(k); Orbit &ol = ms.GetOrbit(l);
        for (int j = 0; j < norb; j++)
        {
          Orbit &oj = ms.GetOrbit(j);
          if (2 * oj.n + oj.l > ms.GetEMax3Body()) continue;   // P(kl, j) needs j in the 3B space
          // |klj`; j_i> : parity and Tz of (k l) minus j must match the channel of i
          if ((ok.l + ol.l + oj.l) % 2 != C.par) continue;
          if (ok.tz2 + ol.tz2 - oj.tz2 != C.tz2) continue;
          if (std::abs(2 * Jkl - oj.j2) > C.j2 or 2 * Jkl + oj.j2 < C.j2) continue;
          Q.colIdx[o][(size_t)g * norb + j] = C.klj.size();
          double na = ok.occ, nb = ol.occ, nc = oj.occ;
          double w = na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc;
          if (k == l) w *= 0.5;
          if (std::abs(na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc) > 1e-8) { C.abc.push_back(C.klj.size()); C.occ.push_back(w); }
          C.klj.push_back({g, j});
        }
      }
    }
    // 6j {jj ji Jkl; jc J' Jab}: twice-integer arguments (j2, j2, 2J, j2, tJ, 2J)
    // jj, Jkl, jc, J', Jab belong to 3B states; ji (the row orbit i) runs over all orbits
    Q.d6[0] = P.j2max3 + 1; Q.d6[1] = P.j2max + 1; Q.d6[2] = 2 * P.J2max3 + 1; Q.d6[3] = P.j2max3 + 1; Q.d6[4] = P.tJmax3 + 2; Q.d6[5] = 2 * P.J2max3 + 1;
    Q.six.assign((size_t)Q.d6[0] * Q.d6[1] * Q.d6[2] * Q.d6[3] * Q.d6[4] * Q.d6[5], 0.0);
    for (int a = 1; a < Q.d6[0]; a += 2)
      for (int b = 1; b < Q.d6[1]; b += 2)
        for (int c = 0; c < Q.d6[2]; c += 2)
          for (int d = 1; d < Q.d6[3]; d += 2)
            for (int e = 1; e < Q.d6[4]; e += 2)
              for (int f = 0; f < Q.d6[5]; f += 2)
              {
                if (std::abs(a - b) > c or a + b < c or std::abs(a - e) > f or a + e < f) continue;
                if (std::abs(d - b) > f or d + b < f or std::abs(d - e) > c or d + e < c) continue;
                Q.six[Q.i6(a, b, c, d, e, f)] = AngMom::SixJ(0.5 * a, 0.5 * b, 0.5 * c, 0.5 * d, 0.5 * e, 0.5 * f);
              }
    return Q;
  }
} // namespace perfplan

void comm232ss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  bool x3 = X.ThreeBody.IsAllocated(), y3 = Y.ThreeBody.IsAllocated();
  if (not(x3 or y3)) return;
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0)
  {
    comm232ss(X, Y, Z);
    return;
  }
  int hermX = X.IsHermitian() ? 1 : -1, hermY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(x3 ? X : Y, Z);
  perfplan::Plan232 &Q = perfplan::get232(Z, P);
  ModelSpace &ms = *Z.modelspace;
  int norb = P.norb;
  static const std::vector<arma::mat> empty;
  static const std::vector<perfplan::wmat> emptyW;
  bool use_y3 = y3 and Y.ThreeBodyNorm() >= 1e-8, use_x3 = x3 and X.ThreeBodyNorm() >= 1e-8;
  const std::vector<perfplan::wmat> &WX = use_x3 ? perfplan::unpack(X, P) : emptyW;
  const std::vector<perfplan::wmat> &WY = use_y3 ? perfplan::unpack(Y, P) : emptyW;
  size_t nobc = Q.obc.size();
  std::vector<arma::mat> ZMAT(nobc);
  double t = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t o = 0; o < nobc; o++)
  {
    auto &C = Q.obc[o];
    size_t di = C.orbs.size(), dk = C.klj.size(), da = C.abc.size();
    ZMAT[o].zeros(di, dk);
    if (di == 0 or dk == 0 or da == 0) continue;
    arma::mat X2(di, da, arma::fill::zeros), Y2(di, da, arma::fill::zeros), X3m(da, dk, arma::fill::zeros), Y3m(da, dk, arma::fill::zeros);
    int j2i = C.j2;
    for (size_t ia = 0; ia < da; ia++)   // Step 1: 2B part
    {
      int col = C.abc[ia];
      int gab = C.klj[col][0], c = C.klj[col][1];
      int a = P.g2[gab][0], b = P.g2[gab][1], Jab = P.g2[gab][2];
      double f2 = -sqrt(2 * Jab + 1.) * C.occ[ia];
      for (size_t ii = 0; ii < di; ii++)
      {
        int i = C.orbs[ii];
        if (use_y3) X2(ii, ia) = f2 * X.TwoBody.GetTBME_J(Jab, c, i, a, b);
        if (use_x3) Y2(ii, ia) = f2 * Y.TwoBody.GetTBME_J(Jab, c, i, a, b);
      }
    }
    // Step 2: 3B part, blocked by (c, j): rows ia = (gab, c) and columns ik = (gkl, j) with fixed c, j share parity / Tz, so
    // for each J' the block is a gathered submatrix W(P(gab, j), P(gkl, c)) times a (Jab, Jkl) weight table.
    {
      int nJ = P.J2max + 1;
      std::vector<std::vector<int>> rowsByC(norb), colsByJ(norb);
      for (size_t ia = 0; ia < da; ia++) rowsByC[C.klj[C.abc[ia]][1]].push_back(ia);
      for (size_t ik = 0; ik < dk; ik++) colsByJ[C.klj[ik][1]].push_back(ik);
      std::vector<int> r_ia, r_p, r_J, q_ik, q_p, q_J;
      std::vector<double> wt(nJ * nJ);
      for (int c = 0; c < norb; c++)
      {
        auto &Rc = rowsByC[c];
        if (Rc.empty()) continue;
        Orbit &oc = ms.GetOrbit(c);
        int g0 = C.klj[C.abc[Rc[0]]][0];
        Orbit &oa0 = ms.GetOrbit(P.g2[g0][0]); Orbit &ob0 = ms.GetOrbit(P.g2[g0][1]);
        for (int j = 0; j < norb; j++)
        {
          auto &Qj = colsByJ[j];
          if (Qj.empty()) continue;
          Orbit &oj = ms.GetOrbit(j);
          int par = (oa0.l + ob0.l + oj.l) % 2, tz = oa0.tz2 + ob0.tz2 + oj.tz2;
          for (int tJp = 1; tJp <= P.tJmax; tJp += 2)
          {
            int ch3 = P.ch3tab[P.i3(tJp, par, tz)];
            if (ch3 < 0) continue;
            r_ia.clear(); r_p.clear(); r_J.clear(); q_ik.clear(); q_p.clear(); q_J.clear();
            for (int ia : Rc)
            {
              // imsrg3_no_qqq (as in comm232ss_srs_optimized_old): skip 3B states (a b j) made of three q orbits
              if (imsrg3_no_qqq and ms.GetOrbit(P.g2[C.klj[C.abc[ia]][0]][0]).cvq + ms.GetOrbit(P.g2[C.klj[C.abc[ia]][0]][1]).cvq + oj.cvq > 5) continue;
              int gab = C.klj[C.abc[ia]][0], Jab = P.g2[gab][2];
              if (std::abs(2 * Jab - oj.j2) > tJp or 2 * Jab + oj.j2 < tJp) continue;
              int p1 = P.pidx(ch3, gab, j);
              if (p1 >= 0) { r_ia.push_back(ia); r_p.push_back(p1); r_J.push_back(Jab); }
            }
            if (r_ia.empty()) continue;
            for (int ik : Qj)
            {
              if (imsrg3_no_qqq and ms.GetOrbit(P.g2[C.klj[ik][0]][0]).cvq + ms.GetOrbit(P.g2[C.klj[ik][0]][1]).cvq + oc.cvq > 5) continue;   // (k l c)
              int gkl = C.klj[ik][0], Jkl = P.g2[gkl][2];
              if (std::abs(2 * Jkl - oc.j2) > tJp or 2 * Jkl + oc.j2 < tJp) continue;
              int p2 = P.pidx(ch3, gkl, c);
              if (p2 >= 0) { q_ik.push_back(ik); q_p.push_back(p2); q_J.push_back(Jkl); }
            }
            if (q_ik.empty()) continue;
            for (int Jab = 0; Jab < nJ; Jab++)
              for (int Jkl = 0; Jkl < nJ; Jkl++)
                wt[Jab * nJ + Jkl] = (tJp + 1.) / sqrt(2 * Jkl + 1.) * Q.six6(oj.j2, j2i, 2 * Jkl, oc.j2, tJp, 2 * Jab);
            size_t nr = r_ia.size();
            for (size_t q = 0; q < q_ik.size(); q++)
            {
              const double *wq = &wt[q_J[q]];
              if (use_y3)
              {
                const auto *wc = WY[ch3].colptr(q_p[q]); double *yc = Y3m.colptr(q_ik[q]);
                for (size_t r = 0; r < nr; r++) yc[r_ia[r]] += wq[r_J[r] * nJ] * wc[r_p[r]];
              }
              if (use_x3)
              {
                const auto *wc = WX[ch3].colptr(q_p[q]); double *xc = X3m.colptr(q_ik[q]);
                for (size_t r = 0; r < nr; r++) xc[r_ia[r]] += wq[r_J[r] * nJ] * wc[r_p[r]];
              }
            }
          }
        }
      }
    }
    ZMAT[o] = X2 * Y3m - Y2 * X3m;   // Step 3
  }
  Z.profiler.timer["_comm232_plan_fill_gemm"] += omp_get_wtime() - t; t = omp_get_wtime();
  auto &Z2 = Z.TwoBody;
  std::vector<std::array<size_t, 2>> bk;
  for (auto &it : Z2.MatEl)
    for (size_t ibra = 0; ibra < ms.GetTwoBodyChannel(it.first[0]).GetNumberKets(); ibra++) bk.push_back({it.first[0], ibra});
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < bk.size(); ib++)
  {
    size_t ch = bk[ib][0], ibra = bk[ib][1];
    TwoBodyChannel &tbc = ms.GetTwoBodyChannel(ch);
    int J = tbc.J;
    Ket &bra = tbc.GetKet(ibra);
    int i = bra.p, j = bra.q;
    Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
    int gij = P.g2idx[((size_t)i * norb + j) * (P.J2max + 1) + J];
    if (imsrg3_valence_2b and (oi.cvq != 1 or oj.cvq != 1)) continue;   // valence 2B output only
    for (size_t iket = ibra; iket < tbc.GetNumberKets(); iket++)
    {
      Ket &ket = tbc.GetKet(iket);
      int k = ket.p, l = ket.q;
      Orbit &ok = ms.GetOrbit(k); Orbit &ol = ms.GetOrbit(l);
      if (imsrg3_valence_2b and (ok.cvq != 1 or ol.cvq != 1)) continue;
      int gkl = P.g2idx[((size_t)k * norb + l) * (P.J2max + 1) + J];
      int phase_ij = ms.phase((oi.j2 + oj.j2) / 2 - J), phase_kl = ms.phase((ok.j2 + ol.j2) / 2 - J);
      auto get = [&](int r, int g, int s) -> double
      {
        int o = Q.obcOf[r];
        int c = (g >= 0 and g < P.ng3) ? Q.colIdx[o][(size_t)g * norb + s] : -1;
        return c < 0 ? 0.0 : ZMAT[o](Q.posOf[r], c);
      };
      double Z_iklj = get(i, gkl, j), Z_jkli = get(j, gkl, i), Z_kijl = get(k, gij, l), Z_lijk = get(l, gij, k);
      double zijkl = phase_ij * Z_iklj - Z_jkli + hermX * hermY * (Z_lijk - phase_kl * Z_kijl);
      zijkl /= sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
      Z2.AddToTBME(ch, ch, ibra, iket, zijkl);
    }
  }
  Z.profiler.timer["_comm232_plan_back"] += omp_get_wtime() - t;
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
