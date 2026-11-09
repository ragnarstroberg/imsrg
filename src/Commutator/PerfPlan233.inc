//=====================================================================================================================
//  comm233_phss with precompiled bookkeeping (scalar X, Y). Notation: see the reader's guide in PerfPlan.inc.
//=====================================================================================================================
//  Formula (IMSRG3Commutators.cc):
//    Z_ijklmn = sum_{ab} (n'a nb - na n'b) P(ij/k) P(lm/n) Y_{ija,lmb} X_{bk,an}  -  (X <-> Y),
//  evaluated, as in the original, as a matrix product in particle-hole (cross-coupled, "CC") channels
//  (Jph, parity, |dTz|) after a Pandya transformation of the 2B and the 3B operator:
//    Step 1 (unpack)  W_X = R X R^T, W_Y = R Y R^T (cached).
//    Step 2 (2B)      X2_ph(kn, ab) = -(na - nb) sum_J' (2J'+1) {jk jn Jph; ja jb J'} X^{J'}_{kb,an}, rows k n in both
//                     orders, columns the ph kets (ab) in both orders (unchanged from comm233_phss).
//    Step 3 (3B)      X3bar(col, ab) = - sum_J' (-1)^{(ja + J')} (2J'+1) {Jij Jlm Jph; jb ja J'} W_X^{J'}(P(gij, a), P(glm, b)),
//                     col = (glm, gij), the pairs in stored order (ij <= lm); stored transposed as X3T(col, ab).
//    Step 4 (GEMM)    Z3_ph = Y3T X2_ph^T - X3T Y2_ph^T          (rows: columns col, columns: CC kets (I3 I6), both orders)
//    Step 5 (back)    W_Z(P(g12, I3), P(g45, I6)) = sum_Jph (2Jph+1) {J12 J45 Jph; j6 j3 J} (phase) Z3_ph(col(g45, g12), (I3 I6))
//                                                 + (flip phase) Z3_ph(col(g12, g45), (I6 I3))   [the "patch-P" flip]
//    Step 6           Z += A^T W_Z A  per 3B channel.
//  E3max < 3 emax ("blocked" CC channels): Z3_ph(col, (p q)) is only read when P(gij, p) and P(glm, q) exist, i.e.
//  e(gij) + e_p <= E3max and e(glm) + e_q <= E3max. CC kets are grouped by (e_p, e_q); per group only those rows are
//  computed (step 4 on contiguous row ranges, see CC233::perm) and stored, per glm a zero row plus the needed prefix of
//  its columns sorted by e(gij); step 5 reads them through CC233::ctabL / Cls::off. Fill and back are unchanged.
//
//  Notes for imsrg++ readers: steps 2, 5 and the phases are those of comm233_phss (including the patch-P fix of the
//  flip phase for two anti-Hermitian operators); only the 3B Pandya sum (step 3) reads W instead of GetME_pn, and the
//  back-transformation goes through W_Z and A. Scalar X, Y only (isospin-changing Y: PerfPlan233tz.inc).
//  Included after PerfPlan.inc (inside namespace Commutator).
namespace perfplan
{
  struct CC233
  {
    std::vector<std::array<int, 2>> cols;   // column -> (g_lm, g_ij), stored ordering (ij <= lm)
    // compact column table: for g_lm and a partner class k (pair Tz, parity) compatible with the channel,
    // ctab[cbase[g_lm * 6 + k] + pos(g_ij)] -> column, or cols.size() (a zero sentinel row of Z3_ph); pos = rank within class
    std::vector<int> ctab;
    std::vector<long> cbase;
    // E3max < 3 emax ("blocked"): Z3_ph(col, k) is only read for e(g_ij) + e(first orbit of k) <= E3max and
    // e(g_lm) + e(second orbit of k) <= E3max. The CC kets k (both orders) are grouped in classes of equal orbit
    // energies (c, d). Per class the stored rows are, for each g_lm with e(g_lm) <= E3max - d, a zero row followed by
    // the prefix of that g_lm's columns (sorted by e(g_ij)) with e(g_ij) <= E3max - c. ctabL is ctab with the column
    // replaced by its 1-based position in the energy-sorted segment of its g_lm (0 = missing -> the zero row).
    bool blocked = false;
    std::vector<int> ctabL;
    std::vector<int> kcls, kpos;                       // CC ket k -> class, position in the class
    struct Cls { int c = 0, d = 0; std::vector<int> kets; std::vector<long> off; size_t nrows = 0; std::vector<int> rows;
                 std::vector<std::array<int, 2>> rng; };   // [start, start + len) of the perm positions this class needs
    std::vector<Cls> cls;                              // rows: column index per stored row (-1 = zero row)
    // X3T, Y3T rows in perm order: columns sorted by (e(g_ij), e(g_lm)), so the rows of a class are <= 2 emax + 1 ranges
    std::vector<int> perm, iperm;                      // perm position -> column, column -> perm position
    std::vector<int> pdst;                             // [class * ncols + perm position] -> stored row in that class
  };
  struct Plan233
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1;
    std::vector<CC233> cc;
    std::vector<int> rank;                  // g -> rank of (p, q, J) in lexicographic order (stored-ordering test)
    std::vector<int> cls, posC;             // g -> pair class (Tz+1)*2 + parity; position within its class in rank order
    std::vector<int> chph;                  // [Jph * 2 * 2 + parity * 2 + Tz] -> CC channel (Tz = |dTz|/2 in {0,1})
    std::vector<std::vector<int>> ccLocal;  // ch_cc -> [p * norb + q] -> local index in the CC channel or -1
    std::vector<double> six;                // 6j {a b c; d e f}, twice-integer arguments
    std::vector<double> six_back;           // back-transformation 6j {J12 J45 Jph; j6 j3 J}, contiguous in (J12, J45)
    int bnJ = 0, bnj = 0, bnPh = 0;
    size_t ib6(int tJ, int j6, int j3, int Jph) const { return ((((size_t)tJ * bnj + j6) * bnj + j3) * bnPh + Jph) * bnJ * bnJ; }
    int d6[6];
    size_t i6(int a, int b, int c, int d, int e, int f) const
    { return ((((((size_t)a * d6[1] + b) * d6[2] + c) * d6[3] + d) * d6[4] + e) * d6[5] + f); }
    // bounds-checked read: arguments outside the table (orbits or pairs outside the 3B space) give 0
    double six6(int a, int b, int c, int d, int e, int f) const
    {
      if (a < 0 or b < 0 or c < 0 or d < 0 or e < 0 or f < 0 or a >= d6[0] or b >= d6[1] or c >= d6[2] or d >= d6[3] or e >= d6[4] or f >= d6[5]) return 0.0;
      return six[i6(a, b, c, d, e, f)];
    }
  };

  static Plan233 &get233(const Operator &Z, Plan &P)
  {
    static Plan233 Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ) return Q;
    double t0 = omp_get_wtime();
    Q = Plan233(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max();
    int ng = P.g2.size(), ng3 = P.ng3;
    size_t nch2_CC = ms.TwoBodyChannels_CC.size();
    Q.cc.resize(nch2_CC);
    auto before = [&](int ga, int gb)   // ga <= gb in (min, max, J) lexicographic order
    {
      auto &A = P.g2[ga]; auto &B = P.g2[gb];
      if (A[0] != B[0]) return A[0] < B[0];
      if (A[1] != B[1]) return A[1] < B[1];
      return A[2] <= B[2];
    };
    {
      std::vector<int> ord(ng);
      for (int g = 0; g < ng; g++) ord[g] = g;
      std::sort(ord.begin(), ord.end(), [&](int a, int b) { return P.g2[a] < P.g2[b]; });
      Q.rank.assign(ng, 0);
      for (int r = 0; r < ng; r++) Q.rank[ord[r]] = r;
    }
    Q.cls.resize(ng); Q.posC.resize(ng);
    std::vector<std::vector<int>> byCls(6);   // pairs of each class in rank order
    {
      std::vector<int> ord(ng);
      for (int g = 0; g < ng; g++) ord[g] = g;
      std::sort(ord.begin(), ord.end(), [&](int a, int b) { return Q.rank[a] < Q.rank[b]; });
      for (int g : ord)
      {
        if (g >= P.ng3) continue;   // columns are pairs of 3B states: pairs of the 3B space only
        Orbit &o1 = ms.GetOrbit(P.g2[g][0]); Orbit &o2 = ms.GetOrbit(P.g2[g][1]);
        int k = ((o1.tz2 + o2.tz2) / 2 + 1) * 2 + (o1.l + o2.l) % 2;
        Q.cls[g] = k; Q.posC[g] = byCls[k].size(); byCls[k].push_back(g);
      }
    }
    auto eo = [&](int o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l; };
    auto ep = [&](int g) { return eo(P.g2[g][0]) + eo(P.g2[g][1]); };
    for (size_t ch = 0; ch < nch2_CC; ch++)
    {
      auto &tbc_CC = ms.GetTwoBodyChannel_CC(ch);
      auto &C = Q.cc[ch];
      std::vector<int> colIdx((size_t)ng3 * ng3, -1);
      for (int glm = 0; glm < ng3; glm++)
        for (int k = 0; k < 6; k++)
        for (int gij : byCls[k])   // columns of one glm ordered by (class, rank): nearly contiguous reads in the back
        {
          int Jlm = P.g2[glm][2], Jij = P.g2[gij][2];
          Orbit &ol = ms.GetOrbit(P.g2[glm][0]); Orbit &om = ms.GetOrbit(P.g2[glm][1]);
          Orbit &oi = ms.GetOrbit(P.g2[gij][0]); Orbit &oj = ms.GetOrbit(P.g2[gij][1]);
          if ((ol.l + om.l + oi.l + oj.l) % 2 != tbc_CC.parity) continue;
          if (std::abs((ol.tz2 + om.tz2) / 2 - (oi.tz2 + oj.tz2) / 2) != tbc_CC.Tz) continue;
          if (std::abs(Jlm - Jij) > tbc_CC.J or Jlm + Jij < tbc_CC.J) continue;
          if (not before(gij, glm)) continue;
          colIdx[(size_t)glm * ng3 + gij] = C.cols.size();
          C.cols.push_back({glm, gij});
        }
      int sentinel = C.cols.size();
      C.cbase.assign((size_t)ng3 * 6, -1);
      std::vector<int> segpos(C.cols.size(), 0);   // 1-based position of a column in the energy-sorted segment of its g_lm
      std::vector<std::vector<int>> seg(ng3);      // g_lm -> its columns sorted by e(g_ij)
      for (size_t col = 0; col < C.cols.size(); col++) seg[C.cols[col][0]].push_back(col);
      for (int glm = 0; glm < ng3; glm++)
      {
        std::stable_sort(seg[glm].begin(), seg[glm].end(), [&](int a, int b) { return ep(C.cols[a][1]) < ep(C.cols[b][1]); });
        for (size_t t = 0; t < seg[glm].size(); t++) segpos[seg[glm][t]] = t + 1;
      }
      for (int glm = 0; glm < ng3; glm++)
      {
        Orbit &ol = ms.GetOrbit(P.g2[glm][0]); Orbit &om = ms.GetOrbit(P.g2[glm][1]);
        int Tzlm = (ol.tz2 + om.tz2) / 2, parlm = (ol.l + om.l) % 2;
        for (int k = 0; k < 6; k++)
        {
          int Tzk = k / 2 - 1, park = k % 2;
          if ((parlm + park) % 2 != tbc_CC.parity or std::abs(Tzlm - Tzk) != tbc_CC.Tz or byCls[k].empty()) continue;
          C.cbase[(size_t)glm * 6 + k] = C.ctab.size();
          for (int gij : byCls[k]) { int c = colIdx[(size_t)glm * ng3 + gij]; C.ctab.push_back(c < 0 ? sentinel : c); C.ctabL.push_back(c < 0 ? 0 : segpos[c]); }
        }
      }
      // classes of CC kets and their stored rows
      {
        size_t nkc = tbc_CC.GetNumberKets();
        C.kcls.assign(2 * nkc, -1); C.kpos.assign(2 * nkc, -1);
        std::map<std::pair<int, int>, int> cid;
        for (size_t k = 0; k < 2 * nkc; k++)
        {
          Ket &kt = tbc_CC.GetKet(k % nkc);
          int first = (k < nkc) ? kt.p : kt.q, second = (k < nkc) ? kt.q : kt.p;
          std::pair<int, int> key = {eo(first), eo(second)};
          auto it = cid.find(key);
          if (it == cid.end()) { it = cid.insert({key, (int)C.cls.size()}).first; C.cls.emplace_back(); C.cls.back().c = key.first; C.cls.back().d = key.second; }
          C.kcls[k] = it->second; C.kpos[k] = C.cls[it->second].kets.size(); C.cls[it->second].kets.push_back(k);
        }
        size_t used = 0, full = (C.cols.size() + 1) * 2 * nkc;
        const int e3b = ms.GetEMax3Body();
        for (auto &K : C.cls)
        {
          K.off.assign(ng3, -1);
          // Z3_ph(col, (p q)) is read only when P(g_ij, p) and P(g_lm, q) exist: classes with an orbit outside the
          // 3B space (e > emax_3body) need no rows
          if (K.c > e3b or K.d > e3b) { K.nrows = 0; continue; }
          for (int glm = 0; glm < ng3; glm++)
          {
            if (seg[glm].empty() or ep(glm) > Q.E3max - K.d) continue;
            K.off[glm] = K.rows.size();
            K.rows.push_back(-1);
            for (int col : seg[glm]) { if (ep(C.cols[col][1]) > Q.E3max - K.c) break; K.rows.push_back(col); }
          }
          K.nrows = K.rows.size();
          used += K.nrows * K.kets.size();
        }
        C.blocked = used < 0.9 * full;
        if (C.blocked)
        {
          size_t ncol = C.cols.size();
          C.perm.resize(ncol);
          for (size_t col = 0; col < ncol; col++) C.perm[col] = col;
          std::stable_sort(C.perm.begin(), C.perm.end(), [&](int a, int b) {
            int ba = ep(C.cols[a][1]), bb = ep(C.cols[b][1]);
            return ba != bb ? ba < bb : ep(C.cols[a][0]) < ep(C.cols[b][0]); });
          C.iperm.resize(ncol);
          for (size_t u = 0; u < ncol; u++) C.iperm[C.perm[u]] = u;
          C.pdst.assign(C.cls.size() * ncol, -1);
          for (size_t t = 0; t < C.cls.size(); t++)
          {
            auto &K = C.cls[t];
            if (K.c > e3b or K.d > e3b) continue;   // no rows (see above)
            size_t j = 0;
            while (j < ncol)
            {
              int b = ep(C.cols[C.perm[j]][1]);
              size_t j0 = j;
              while (j < ncol and ep(C.cols[C.perm[j]][1]) == b) j++;
              if (b > Q.E3max - K.c) continue;
              size_t len = 0;
              while (j0 + len < j and ep(C.cols[C.perm[j0 + len]][0]) <= Q.E3max - K.d) len++;
              if (len == 0) continue;
              K.rng.push_back({(int)j0, (int)len});
              for (size_t u = j0; u < j0 + len; u++)
              {
                int col = C.perm[u];
                C.pdst[t * ncol + u] = K.off[C.cols[col][0]] + segpos[col];
              }
            }
          }
        }
        if (not C.blocked) { C.cls.clear(); C.kcls.clear(); C.kpos.clear(); C.ctabL.clear(); C.ctabL.shrink_to_fit(); }
      }
    }
    int norb = P.norb;
    Q.chph.assign((size_t)(2 * P.J2max + 2) * 4, -1);
    for (int J = 0; J <= 2 * P.J2max + 1; J++)          // exactly what the original does: GetTwoBodyChannelIndex(Jph, parity, |dTz|)
      for (int par = 0; par <= 1; par++)
        for (int Tz = 0; Tz <= 1; Tz++)
        {
          int ch = ms.GetTwoBodyChannelIndex(J, par, Tz);
          if (ch >= 0 and ch < (int)nch2_CC) Q.chph[(size_t)J * 4 + par * 2 + Tz] = ch;
        }
    Q.ccLocal.resize(nch2_CC);
    for (size_t ch = 0; ch < nch2_CC; ch++)
    {
      auto &tbc_CC = ms.GetTwoBodyChannel_CC(ch);
      Q.ccLocal[ch].assign((size_t)norb * norb, -1);
      for (int p = 0; p < norb; p++)                       // exactly what the original does: GetLocalIndex(p, q)
        for (int q = 0; q < norb; q++)
        {
          size_t k = tbc_CC.GetLocalIndex(p, q);
          if (k < 2 * tbc_CC.GetNumberKets()) Q.ccLocal[ch][(size_t)p * norb + q] = (int)k;   // reversed order maps into the second half
        }
    }
    // 6j table of shape (2J, 2J, 2J, j2, j2, tJ): fill uses {Jij Jlm Jph; jb ja J'}, back uses {J1p J2p Jph; j6 j3 Jtot}
    // all arguments belong to 3B states (pairs and orbits of the 3B space)
    Q.d6[0] = 2 * P.J2max3 + 1; Q.d6[1] = 2 * P.J2max3 + 1; Q.d6[2] = 2 * P.J2max3 + P.j2max3 + 2;
    Q.d6[3] = P.j2max3 + 1; Q.d6[4] = P.j2max3 + 1; Q.d6[5] = P.tJmax3 + 2;
    Q.six.assign((size_t)Q.d6[0] * Q.d6[1] * Q.d6[2] * Q.d6[3] * Q.d6[4] * Q.d6[5], 0.0);
    for (int a = 0; a < Q.d6[0]; a += 2)
      for (int b = 0; b < Q.d6[1]; b += 2)
        for (int c = 0; c < Q.d6[2]; c += 2)
          for (int d = 1; d < Q.d6[3]; d += 2)
            for (int e = 1; e < Q.d6[4]; e += 2)
              for (int f = 1; f < Q.d6[5]; f += 2)
              {
                if (std::abs(a - b) > c or a + b < c or std::abs(a - e) > f or a + e < f) continue;
                if (std::abs(d - b) > f or d + b < f or std::abs(d - e) > c or d + e < c) continue;
                Q.six[Q.i6(a, b, c, d, e, f)] = AngMom::SixJ(0.5 * a, 0.5 * b, 0.5 * c, 0.5 * d, 0.5 * e, 0.5 * f);
              }
    Q.bnJ = P.J2max + 1; Q.bnj = P.j2max3 + 1; Q.bnPh = P.j2max3 + 1;   // bnJ: the kernels' nJ stride
    Q.six_back.assign((size_t)(P.tJmax3 + 2) * Q.bnj * Q.bnj * Q.bnPh * Q.bnJ * Q.bnJ, 0.0);
    for (int tJ = 1; tJ <= P.tJmax3 + 1; tJ += 2)
      for (int j6 = 1; j6 <= P.j2max3; j6 += 2)
        for (int j3 = 1; j3 <= P.j2max3; j3 += 2)
          for (int Jph = 0; Jph < Q.bnPh; Jph++)
          {
            if (2 * Jph >= Q.d6[2]) continue;
            double *b = &Q.six_back[Q.ib6(tJ, j6, j3, Jph)];
            for (int J12 = 0; J12 < Q.bnJ; J12++)
              for (int J45 = 0; J45 < Q.bnJ; J45++)
                if (tJ < Q.d6[5]) b[J12 * Q.bnJ + J45] = Q.six6(2 * J12, 2 * J45, 2 * Jph, j6, j3, tJ);
          }
    Z.profiler.timer["_comm233_ph_plan_build"] += omp_get_wtime() - t0;
    if (std::getenv("IMSRG_PLAN_DEBUG") != nullptr)
    {
      double gc = 0, gt = 0, gp = 0, gr = 0, go = 0;   // GB: cols, ctab/ctabL/cbase, perm/iperm/pdst, class rows, class off + kcls/kpos
      for (auto &C : Q.cc)
      {
        gc += 8e-9 * C.cols.size(); gt += 4e-9 * (C.ctab.size() + C.ctabL.size()) + 8e-9 * C.cbase.size();
        gp += 4e-9 * (C.perm.size() + C.iperm.size() + C.pdst.size());
        for (auto &K : C.cls) { gr += 4e-9 * K.rows.size() + 4e-9 * K.kets.size(); go += 8e-9 * K.off.size(); }
        go += 4e-9 * (C.kcls.size() + C.kpos.size());
      }
      std::cout << "PLANDEBUG comm233_ph plan memory (GB): cols " << gc << ", ctab " << gt << ", perm/pdst " << gp << ", class rows " << gr
                << ", class off " << go << ", 6j " << 8e-9 * (Q.six.size() + Q.six_back.size()) << std::endl;
    }
    return Q;
  }
} // namespace perfplan

void comm233_phss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0)
  {
    comm233_phss(X, Y, Z);
    return;
  }
  auto &X2 = X.TwoBody; auto &Y2 = Y.TwoBody; auto &X3 = X.ThreeBody; auto &Y3 = Y.ThreeBody; auto &Z3 = Z.ThreeBody;
  bool x3 = X3.IsAllocated(), y3 = Y3.IsAllocated();
  if (not(x3 or y3)) return;
  int hX = X.IsHermitian() ? 1 : -1, hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(x3 ? X : Y, Z);
  perfplan::Plan233 &Q = perfplan::get233(Z, P);
  ModelSpace &ms = *Z.modelspace;
  const perfplan::Cut3 &cutZ = perfplan::cut3(ms);
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  int norb = P.norb;
  double t = omp_get_wtime();
  static const std::vector<arma::mat> empty;
  static const std::vector<perfplan::wmat> emptyW;
  // Step 1: unpack (cached)
  const std::vector<perfplan::wmat> &WX = x3 ? perfplan::unpack(X, P) : emptyW;
  const std::vector<perfplan::wmat> &WY = y3 ? perfplan::unpack(Y, P) : emptyW;
  Z.profiler.timer["_comm233_ph_plan_unpack"] += omp_get_wtime() - t; t = omp_get_wtime();

  size_t nch2_CC = ms.TwoBodyChannels_CC.size();
  std::vector<arma::mat> Z3_ph(nch2_CC);
  std::vector<std::vector<arma::mat>> Z3b(nch2_CC);   // blocked: per class (rows x kets of the class)

  double tt_f2 = 0, tt_f3 = 0, tt_mm = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : tt_f2, tt_f3, tt_mm)
  for (size_t ch_cc = 0; ch_cc < nch2_CC; ch_cc++)
  {
    double tq0 = omp_get_wtime();
    TwoBodyChannel_CC &tbc_CC = ms.GetTwoBodyChannel_CC(ch_cc);
    int Jph = tbc_CC.J;
    size_t nkets_CC = tbc_CC.GetNumberKets();
    auto &C = Q.cc[ch_cc];
    size_t nkets3 = C.cols.size();
    std::vector<size_t> ph_kets;   // summed ph pairs (a b): X3T(col, ab) = 0 unless both orbits are in the 3B space
    auto in3 = [&](const Orbit *o) { return 2 * o->n + o->l <= ms.GetEMax3Body(); };
    for (size_t iket = 0; iket < nkets_CC; iket++)
    {
      Ket &ket = tbc_CC.GetKet(iket);
      if (std::abs(ket.op->occ - ket.oq->occ) > 1e-7 and in3(ket.op) and in3(ket.oq)) ph_kets.push_back(iket);
    }
    size_t nkets_ph = ph_kets.size();
    const bool blk = C.blocked;
    if (blk)
    {
      Z3b[ch_cc].resize(C.cls.size());
      for (size_t t = 0; t < C.cls.size(); t++) Z3b[ch_cc][t].zeros(C.cls[t].nrows, C.cls[t].kets.size());
    }
    else Z3_ph[ch_cc].zeros(nkets3 + 1, 2 * nkets_CC);   // transposed layout, last row = zero sentinel
    if (nkets_ph == 0 or nkets3 == 0) continue;
    arma::mat X2_ph(2 * nkets_CC, 2 * nkets_ph, arma::fill::zeros), Y2_ph(2 * nkets_CC, 2 * nkets_ph, arma::fill::zeros);
    // Step 2: 2B Pandya transform X2_ph, Y2_ph (unchanged from comm233_phss)
    for (size_t ibra_ij = 0; ibra_ij < nkets_CC; ibra_ij++)
    {
      Ket &bra_ij = tbc_CC.GetKet(ibra_ij);
      size_t i = bra_ij.p, j = bra_ij.q;
      Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
      if (not(in3(&oi) and in3(&oj))) continue;   // Z3_ph columns (I3 I6) are only read for orbits of the 3B space
      double ji = 0.5 * oi.j2, jj = 0.5 * oj.j2;
      for (size_t iket_ab = 0; iket_ab < nkets_ph; iket_ab++)
      {
        Ket &ket_ab = tbc_CC.GetKet(ph_kets[iket_ab]);
        size_t a = ket_ab.p, b = ket_ab.q;
        Orbit &oa = ms.GetOrbit(a); Orbit &ob = ms.GetOrbit(b);
        double ja = 0.5 * oa.j2, jb = 0.5 * ob.j2;
        double occ_factor = oa.occ - ob.occ;
        double xijab = 0, xjiab = 0, yijab = 0, yjiab = 0;
        int Jmin = std::max(std::abs(oi.j2 - ob.j2), std::abs(oa.j2 - oj.j2)) / 2, Jmax = std::min(oi.j2 + ob.j2, oa.j2 + oj.j2) / 2;
        for (int Jp = Jmin; Jp <= Jmax; Jp++)
        {
          double s6 = ms.GetSixJ(ji, jj, Jph, ja, jb, Jp);
          xijab -= (2 * Jp + 1) * s6 * X2.GetTBME_J(Jp, i, b, a, j);
          yijab -= (2 * Jp + 1) * s6 * Y2.GetTBME_J(Jp, i, b, a, j);
        }
        Jmin = std::max(std::abs(oj.j2 - ob.j2), std::abs(oa.j2 - oi.j2)) / 2; Jmax = std::min(oj.j2 + ob.j2, oa.j2 + oi.j2) / 2;
        for (int Jp = Jmin; Jp <= Jmax; Jp++)
        {
          double s6 = ms.GetSixJ(jj, ji, Jph, ja, jb, Jp);
          xjiab -= (2 * Jp + 1) * s6 * X2.GetTBME_J(Jp, j, b, a, i);
          yjiab -= (2 * Jp + 1) * s6 * Y2.GetTBME_J(Jp, j, b, a, i);
        }
        int phase_X = hX * ms.phase((oi.j2 + oj.j2 + oa.j2 + ob.j2) / 2), phase_Y = hY * ms.phase((oi.j2 + oj.j2 + oa.j2 + ob.j2) / 2);
        X2_ph(ibra_ij, iket_ab) = xijab * occ_factor; Y2_ph(ibra_ij, iket_ab) = yijab * occ_factor;
        X2_ph(ibra_ij + nkets_CC, iket_ab) = xjiab * occ_factor; Y2_ph(ibra_ij + nkets_CC, iket_ab) = yjiab * occ_factor;
        X2_ph(ibra_ij, iket_ab + nkets_ph) = phase_X * xjiab * (-occ_factor); Y2_ph(ibra_ij, iket_ab + nkets_ph) = phase_Y * yjiab * (-occ_factor);
        X2_ph(ibra_ij + nkets_CC, iket_ab + nkets_ph) = phase_X * xijab * (-occ_factor); Y2_ph(ibra_ij + nkets_CC, iket_ab + nkets_ph) = phase_Y * yijab * (-occ_factor);
      }
    }
    double tq1 = omp_get_wtime();
    // Step 3: 3B part from W lookups:  X3bar_(ab),(lm ij) = - sum_J' (-1)^{(j2a + 2J')/2} (2J'+1) {Jij Jlm Jph; jb ja J'} X^{J'}_{(ij)a,(lm)b}
    // Stored transposed (cols x 2 nkets_ph); (a b) outer so that, with columns ordered (lm outer, ij inner) and the
    // third-orbit-major P basis, each W column P(lm, b) is read nearly sequentially. 6j weights: small per-(a b) table.
    arma::mat X3T(nkets3, 2 * nkets_ph, arma::fill::zeros), Y3T(nkets3, 2 * nkets_ph, arma::fill::zeros);
    {
      int nJ = P.J2max + 1, nT = P.tJmax / 2 + 2;
      std::vector<double> wt((size_t)nJ * nJ * nT);
      std::vector<int> cJij(nkets3), cJlm(nkets3), cpar(nkets3), ctz(nkets3), cgij(nkets3), cglm(nkets3);
      std::vector<int> cEij(nkets3), cElm(nkets3);   // pair energies, for the E3max pre-check
      auto eo3 = [&](int o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l; };
      for (size_t col = 0; col < nkets3; col++)
      {
        int glm = C.cols[col][0], gij = C.cols[col][1];
        cEij[col] = eo3(P.g2[gij][0]) + eo3(P.g2[gij][1]); cElm[col] = eo3(P.g2[glm][0]) + eo3(P.g2[glm][1]);
        Orbit &oi = ms.GetOrbit(P.g2[gij][0]); Orbit &oj = ms.GetOrbit(P.g2[gij][1]);
        cJij[col] = P.g2[gij][2]; cJlm[col] = P.g2[glm][2]; cpar[col] = (oi.l + oj.l) % 2; ctz[col] = oi.tz2 + oj.tz2;
        cgij[col] = gij; cglm[col] = glm;
      }
      for (size_t iab = 0; iab < 2 * nkets_ph; iab++)
      {
        Ket &ket_ab = tbc_CC.GetKet(ph_kets[iab % nkets_ph]);
        size_t a = (iab < nkets_ph) ? ket_ab.p : ket_ab.q, b = (iab < nkets_ph) ? ket_ab.q : ket_ab.p;
        Orbit &oa = ms.GetOrbit(a); Orbit &ob = ms.GetOrbit(b);
        if (std::abs(oa.occ - ob.occ) < 1e-6) continue;
        std::fill(wt.begin(), wt.end(), NAN);
        double *xc = X3T.colptr(iab), *yc = Y3T.colptr(iab);
        // E3max < 3 emax: P(gij, a) and P(glm, b) only exist below these pair energies; checking them first avoids
        // most of the (random) index-table lookups
        const int lim_ij = P.E3max - eo3(a), lim_lm = P.E3max - eo3(b);
        for (size_t col = 0; col < nkets3; col++)
        {
          if (C.blocked and (cEij[col] > lim_ij or cElm[col] > lim_lm)) continue;
          int Jij = cJij[col], Jlm = cJlm[col];
          int par = (cpar[col] + oa.l) % 2, tz = ctz[col] + oa.tz2;
          int lo = std::max(std::abs(oa.j2 - 2 * Jij), std::abs(ob.j2 - 2 * Jlm)), hi = std::min(oa.j2 + 2 * Jij, ob.j2 + 2 * Jlm);
          double xs = 0, ys = 0;
          double *wrow = &wt[((size_t)Jij * nJ + Jlm) * nT];
          for (int tJp = lo; tJp <= hi; tJp += 2)
          {
            int ch3 = P.ch3tab[P.i3(tJp, par, tz)];
            if (ch3 < 0) continue;
            const int *pt = P.PIdxT[ch3].data();   // orbit-major: consecutive columns read nearby entries
            int p1 = pt[(size_t)a * P.ng3 + cgij[col]], p2 = pt[(size_t)b * P.ng3 + cglm[col]];
            if (p1 < 0 or p2 < 0) continue;
            double &w = wrow[tJp / 2];
            if (std::isnan(w)) w = ms.phase((oa.j2 + tJp) / 2) * (tJp + 1) * Q.six6(2 * Jij, 2 * Jlm, 2 * Jph, ob.j2, oa.j2, tJp);
            if (x3) xs -= w * WX[ch3](p1, p2);
            if (y3) ys -= w * WY[ch3](p1, p2);
          }
          size_t wr = blk ? C.iperm[col] : col;   // blocked: rows in perm order
          xc[wr] = xs; yc[wr] = ys;
        }
      }
    }
    double tq2 = omp_get_wtime();
    if (blk)
    {
      // Step 4 (blocked): Z3_ph = Y3T X2_ph^T - X3T Y2_ph^T on the needed entries only. The rows of X3T, Y3T are in perm order, so each
      // class needs a few contiguous row ranges: BLAS reads them in place (leading dimension nkets3), no gathers.
      size_t nd = 2 * nkets_ph, ncol = nkets3;
      const double one = 1.0, zero = 0.0, mone = -1.0;
      arma::blas_int lda = ncol, bk = nd;
      arma::mat XK, YK, Zs;
      for (size_t t = 0; t < C.cls.size(); t++)
      {
        auto &K = C.cls[t];
        if (K.rng.empty()) continue;
        size_t nkt = K.kets.size();
        arma::uvec kk(nkt);
        for (size_t u = 0; u < nkt; u++) kk(u) = K.kets[u];
        XK = X2_ph.rows(kk); YK = Y2_ph.rows(kk);
        arma::mat &B = Z3b[ch_cc][t];
        arma::blas_int n = nkt, ldb = nkt;
        for (auto &rg : K.rng)
        {
          arma::blas_int m = rg[1];
          Zs.set_size(rg[1], nkt);
          arma::blas::gemm<double>("N", "T", &m, &n, &bk, &one, Y3T.memptr() + rg[0], &lda, XK.memptr(), &ldb, &zero, Zs.memptr(), &m);
          arma::blas::gemm<double>("N", "T", &m, &n, &bk, &mone, X3T.memptr() + rg[0], &lda, YK.memptr(), &ldb, &one, Zs.memptr(), &m);
          const int *dst = &C.pdst[t * ncol + rg[0]];
          for (size_t q = 0; q < nkt; q++)
          {
            double *bc = B.colptr(q); const double *sc = Zs.colptr(q);
            for (int u = 0; u < rg[1]; u++) bc[dst[u]] = sc[u];
          }
        }
      }
      tt_f2 += tq1 - tq0; tt_f3 += tq2 - tq1; tt_mm += omp_get_wtime() - tq2;
      continue;
    }
    Z3_ph[ch_cc] = Y3T * X2_ph.t() - X3T * Y2_ph.t();   // Step 4: stored transposed: (cols) x (2 nkets_CC)
    tt_f2 += tq1 - tq0; tt_f3 += tq2 - tq1; tt_mm += omp_get_wtime() - tq2;
    Z3_ph[ch_cc].resize(nkets3 + 1, 2 * nkets_CC);                   // zero sentinel row
  }
  Z.profiler.timer["_comm233_ph_plan_tt_fill2"] += tt_f2; Z.profiler.timer["_comm233_ph_plan_tt_fill3"] += tt_f3; Z.profiler.timer["_comm233_ph_plan_tt_gemm"] += tt_mm;
  Z.profiler.timer["_comm233_ph_plan_fill_gemm"] += omp_get_wtime() - t; t = omp_get_wtime();

  double tt_wz = 0, tt_a = 0, n_it = 0, n_wz = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : tt_wz, tt_a, n_it, n_wz)
  for (size_t ich = 0; ich < nch3; ich++)
  {
    size_t ch = P.ord3[ich];
    auto &T = ms.GetThreeBodyChannel(ch);
    size_t nP = P.Pst[ch].size(), nk = T.GetNumberKets();
    if (nP == 0 or nk == 0) continue;
    int twoJ = T.twoJ;
    double tb0 = omp_get_wtime();
    arma::mat WZ(nP, nP, arma::fill::zeros); n_wz += (double)nP * nP;
    // Step 5 (back-transformation into W_Z), blocked by the third orbits (I3 of the row, I6 of the column): for fixed (I3, I6, Jph) the CC row of Z3_ph is
    // fixed, so it is copied once into a contiguous vector and gathered with the (nearly sequential) column index table.
    std::vector<std::vector<int>> byOrb(norb), gOf(norb), JOf(norb), rkOf(norb), pcOf(norb);
    std::vector<std::vector<char>> hasJ(norb, std::vector<char>(P.J2max + 1, 0));
    for (size_t p = 0; p < nP; p++) byOrb[P.Pst[ch][p][1]].push_back(p);
    for (int o = 0; o < norb; o++)   // sorted by pair rank, so the rank tests below are prefixes (no branches)
    {
      auto &L = byOrb[o];
      std::sort(L.begin(), L.end(), [&](int a, int b) { return Q.rank[P.Pst[ch][a][0]] < Q.rank[P.Pst[ch][b][0]]; });
      for (int p : L)
      {
        int g = P.Pst[ch][p][0];
        gOf[o].push_back(g); JOf[o].push_back(P.g2[g][2]); rkOf[o].push_back(Q.rank[g]); pcOf[o].push_back(Q.posC[g]); hasJ[o][P.g2[g][2]] = 1;
      }
    }
    int nJ = P.J2max + 1;
    std::vector<double> wt(nJ * nJ);
    std::vector<double> Bnv, Bfv;   // block buffers, grown once and reused (no allocation per block)
    for (int I3 = 0; I3 < norb; I3++)
    {
      auto &Rw = byOrb[I3];
      if (Rw.empty()) continue;
      const int *rg = gOf[I3].data(), *rJ = JOf[I3].data(), *rrk = rkOf[I3].data(), *rpc = pcOf[I3].data();
      int clsR = Q.cls[rg[0]];   // all rows (g12, I3) of the channel share the pair class
      size_t nr = Rw.size();
      Orbit &o3 = ms.GetOrbit(I3);
      int phase_ph = -ms.phase((twoJ + o3.j2) / 2);
      for (int I6 = 0; I6 < norb; I6++)
      {
        auto &Cl = byOrb[I6];
        if (Cl.empty()) continue;
        Orbit &o6 = ms.GetOrbit(I6);
        int Tz_ph = std::abs(o3.tz2 - o6.tz2) / 2, parity_ph = (o3.l + o6.l) % 2;
        double fph = -hX * hY * ms.phase((o3.j2 - o6.j2) / 2);
        int Jph_lo = std::abs(o3.j2 - o6.j2) / 2, Jph_hi = (o3.j2 + o6.j2) / 2;
        const int *cg = gOf[I6].data(), *cJ = JOf[I6].data(), *crk = rkOf[I6].data(), *cpc = pcOf[I6].data();
        int clsC = Q.cls[cg[0]];
        size_t nc = Cl.size();
        // (I3, I6) block accumulated in cache-sized buffers: Bn(r, q) for rank(g12) <= rank(g45) (walks colIdx[g45][g12]
        // and the Z3_ph row with g12), Bf(q, r) for the flipped case (walks colIdx[g12][g45] with g45), then scattered.
        if (Bnv.size() < nr * nc) { Bnv.resize(nr * nc); Bfv.resize(nr * nc); }
        std::fill(Bnv.begin(), Bnv.begin() + nr * nc, 0.0); std::fill(Bfv.begin(), Bfv.begin() + nr * nc, 0.0);
        double *Bn = Bnv.data(), *Bf = Bfv.data();   // Bn(r, q) = Bn[r + q * nr], Bf(q, r) = Bf[q + r * nc]
        bool any = false;
        for (int Jph = Jph_lo; Jph <= Jph_hi; Jph++)
        {
          int ch_ph = Q.chph[(size_t)Jph * 4 + parity_ph * 2 + Tz_ph];
          if (ch_ph < 0) continue;
          auto &C = Q.cc[ch_ph];
          auto &Zp = Z3_ph[ch_ph];
          int row_nf = Q.ccLocal[ch_ph][(size_t)I3 * norb + I6], row_f = Q.ccLocal[ch_ph][(size_t)I6 * norb + I3];
          const bool blk = C.blocked;
          int ncolZ = blk ? (int)C.kcls.size() : (int)Zp.n_cols;
          const double *znf = nullptr, *zf = nullptr;
          const long *offnf = nullptr, *offf = nullptr;
          if (blk)
          {
            if (row_nf >= 0 and row_nf < ncolZ) { int t = C.kcls[row_nf]; znf = Z3b[ch_ph][t].colptr(C.kpos[row_nf]); offnf = C.cls[t].off.data(); }
            if (row_f >= 0 and row_f < ncolZ) { int t = C.kcls[row_f]; zf = Z3b[ch_ph][t].colptr(C.kpos[row_f]); offf = C.cls[t].off.data(); }
          }
          else
          {
            znf = (row_nf >= 0 and row_nf < ncolZ) ? Zp.colptr(row_nf) : nullptr;   // contiguous
            zf = (row_f >= 0 and row_f < ncolZ) ? Zp.colptr(row_f) : nullptr;
          }
          if (znf == nullptr and zf == nullptr) continue;
          const std::vector<int> &ct = blk ? C.ctabL : C.ctab;
          any = true; n_it += (double)nr * nc;
          double pref = -phase_ph * (2 * Jph + 1);
          // weight pref * 6j(J12 J45 Jph; j6 j3 J) depends only on (J12, J45): small table
          {
            if (Jph >= Q.bnPh or twoJ >= P.tJmax3 + 2 or o6.j2 >= Q.bnj or o3.j2 >= Q.bnj) continue;   // 6j = 0 (triangle)
            const double *b6 = &Q.six_back[Q.ib6(twoJ, o6.j2, o3.j2, Jph)];
            for (int k = 0; k < nJ * nJ; k++) wt[k] = pref * b6[k];
          }
          if (znf != nullptr)
            for (size_t q = 0; q < nc; q++)
            {
              int rk45 = crk[q];
              long base = C.cbase[(size_t)cg[q] * 6 + clsR];
              if (base < 0) continue;
              const int *colrow = &ct[base];   // column of (g45, g12) = colrow[pos(g12)], pos increasing with r
              const double *wq = &wt[cJ[q]];
              double *bn = Bn + q * nr;
              const double *zq = znf;
              if (blk) { long o = offnf[cg[q]]; if (o < 0) continue; zq = znf + o; }   // segment of g_lm = g45
              size_t split = std::upper_bound(rrk, rrk + nr, rk45) - rrk;   // rank(g12) <= rank(g45)
              for (size_t r = 0; r < split; r++) bn[r] += wq[rJ[r] * nJ] * zq[colrow[rpc[r]]];
            }
          if (zf != nullptr)
            for (size_t r = 0; r < nr; r++)
            {
              int rk12 = rrk[r];
              long base = C.cbase[(size_t)rg[r] * 6 + clsC];
              if (base < 0) continue;
              const int *colrow = &ct[base];   // column of (g12, g45) = colrow[pos(g45)]
              const double *wr = &wt[rJ[r] * nJ];
              double *bf = Bf + r * nc;
              const double *zr = zf;
              if (blk) { long o = offf[rg[r]]; if (o < 0) continue; zr = zf + o; }   // segment of g_lm = g12
              size_t split = std::lower_bound(crk, crk + nc, rk12) - crk;   // rank(g45) < rank(g12)
              for (size_t q = 0; q < split; q++) bf[q] += fph * wr[cJ[q]] * zr[colrow[cpc[q]]];
            }
        }
        if (not any) continue;
        for (size_t q = 0; q < nc; q++)
        {
          double *wc = WZ.colptr(Cl[q]);
          const double *bn = Bn + q * nr;
          for (size_t r = 0; r < nr; r++) wc[Rw[r]] += bn[r] + Bf[q + r * nc];
        }
      }
    }
    double tb1 = omp_get_wtime();
    arma::mat Zk = perfplan::AtWA(P.A[ch], WZ, P.A[ch]);   // Step 6: Z = A^T W_Z A
    perfplan::mask3(cutZ, ch, ch, Zk);   // dE3max / occnat3 cut on the output
    Z3.AddChannelMatrix_pn(ch, Zk, -1);   // = AddToME_pn_ch for iket >= ibra
    tt_wz += tb1 - tb0; tt_a += omp_get_wtime() - tb1;
  }
  Z.profiler.timer["_comm233_ph_plan_tt_back_wz"] += tt_wz;
  Z.profiler.timer["_comm233_ph_plan_tt_back_A"] += tt_a;
  Z.profiler.timer["_comm233_ph_plan_n_it"] += n_it;
  Z.profiler.timer["_comm233_ph_plan_n_wz"] += n_wz;
  Z.profiler.timer["_comm233_ph_plan_back"] += omp_get_wtime() - t;
  perfplan::release_memory();
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
