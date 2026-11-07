//=====================================================================================================================
//  Precompiled-bookkeeping ("perf plan") kernels for the IMSRG(3) commutators: shared machinery + comm333_pph_hhpss
//=====================================================================================================================
//  Enabled with IMSRG_PERF_PLAN=1 (see the dispatch at the top of each term in IMSRG3Commutators.cc). Every kernel
//  computes the same matrix elements as the term it replaces; only the bookkeeping differs. This file holds the pieces
//  that all kernels share; the other PerfPlan*.inc files hold one term each and are included after this one.
//
//  Reader's guide (notation used in all PerfPlan*.inc files)
//  ---------------------------------------------------------
//  * 3B channel  T = (J, parity, 2Tz).  Stored kets |k> = |(ab)Jab c; J>, a <= b <= c, antisymmetric
//    (the storage of ThreeBodyStorage_pn); n_k = number of stored kets of the channel. X(k, l) = <k|X|l>.
//
//  * Product basis P ("pair-antisymmetric product states"): |p> = |(ab)Jab c; J> with a <= b and ANY third orbit c,
//    ordered by the third orbit c (so all p with the same c form a contiguous range).
//    Pst[ch][p] = (g, c), where g indexes the 2B ket (a, b, Jab) in g2[]. PIdx gives p from (g, c), or -1.
//
//  * R (n_P x n_k): expresses a product state in stored kets,  |p> = sum_k R(p, k) |k>
//    (the coefficients of ThreeBodyME::GetKetIndex_withRecoupling). Hence the product-basis matrix of X is
//        W = R X R^T,     W(p, q) = <p|X|q> = sum_{k,l} R(p, k) X(k, l) R(q, l).
//    "Unpacking" an operator means computing, per channel, K = X (stored block), R X, and W.
//
//  * A (n_P x n_k): the antisymmetriser. A result first built in the product basis, W_Z(p, q), is turned into stored
//    matrix elements by
//        Z(k, l) = sum_{p,q} A(p, k) W_Z(p, q) A(q, l),   i.e.  Z = A^T W_Z A,
//    where A(:, k) holds the direct, P(ik) and P(jk) terms of the J-coupled permutation operator P(ij/k)^{J1,J}
//    (6j symbols and phases exactly as in the original back-transformations). A^T ... A implements P(ij/k) P(lm/n).
//
//  * Particle-hole (Pandya) channels, for comm333_pph (below): |e> = |(ij)J12 , n^-1 ; Jph> with pair g = (i <= j, J12)
//    and "hole" orbit n (ext states), and the internal states a = (gab, c) with a nonzero occupation factor (intl).
//
//  * Unpack cache (WCache below): within one Commutator() call X and Y do not change, so they are unpacked once and
//    shared by all terms. Across calls (nested commutators of a BCH series share Omega) each channel's stored matrix
//    is compared bitwise with the cached one and only changed channels are recomputed.
//
//  * E3max < 3 emax: |p> exists only if e_a + e_b + e_c <= E3max (e = 2n + l). Intermediates that are labelled by
//    states without that cut (the pph ext states here, the column pairs of comm233_ph) are computed and stored only
//    where both product states they connect exist ("blocked"); for E3max = 3 emax nothing is cut.
//
//  comm333_pph_hhpss  (formula as in IMSRG3Commutators.cc)
//  -----------------
//    Z_ijklmn = -1/2 P(ij/k) P(lm/n) sum_{abc} (na nb n'c + n'a n'b nc) (X_{ijc,abn} Y_{abk,lmc} - Y_{ijc,abn} X_{abk,lmc})
//  The plan evaluates it as a product in particle-hole channels Jph:
//    Step 1 (unpack)   W_X = R X R^T, W_Y = R Y R^T per 3B channel (cached).
//    Step 2 (fill)     Xbar(e, a) = sum_J' (2J'+1) {jn Jph J12; jc J' Jab} W_X^{J'}(P(g12, c), P(gab, n)),
//                      e = (g12, n) ext, a = (gab, c) intl, and the same for Ybar from W_Y; Ybar(:, a) *= occ(a)
//                      (occupation factor of the intl state, including 1/2 for identical pair orbits).
//    Step 3 (GEMM)     Zbar = hY Xbar Ybar^T,  Zbar -= hX hY Zbar^T          (h = +1 Hermitian, -1 anti-Hermitian)
//    Step 4 (back)     W_Z(P(g12, I3), P(g45, I6)) = sum_Jph (2Jph+1) {j3 J J12; j6 Jph J45} Zbar^{Jph}(e(g12, I6), e(g45, I3))
//    Step 5 (antisym.) Z += A^T W_Z A  per 3B channel.
//  The 6j weights of steps 2 and 4 depend only on the pair J values inside a block of fixed single orbits, so they are
//  small tables, and the gathers run over contiguous P ranges. Steps 2-3 run over pph channels, steps 4-5 over 3B
//  channels, both largest first (dynamic scheduling).
//  For E3max < 3 emax, step 3 only computes the corners of Zbar that step 4 reads (Ch::blocked below).
//
//  Notes for imsrg++ readers: the plan structures depend only on the ModelSpace (and the occupations) and are rebuilt
//  when either changes. Scalar X, Y only; everything else goes to the existing comm333_pph_hhpss (the
//  isospin-changing case is comm333_pph_hhpss_plan_tz in PerfPlan333tz.inc).
//  Included at the end of IMSRG3Commutators.cc (inside namespace Commutator).

namespace perfplan
{
  struct Ch
  {
    int twoJph, par, twoTz;
    std::vector<std::array<int, 2>> ext;   // (g, n): pair g = (i<=j, Jij), "hole" n
    std::vector<std::array<int, 2>> intl;  // (g, c) with nonzero occupation factor
    std::vector<double> occ_int;
    std::vector<int> extIdx;               // [g*norb + n] -> ext index or -1
    // Within each single orbit n the ext states are sorted by the pair energy e(g). The back-transformation only reads
    // Zbar(e1, e2) with e(g1) + e(n2) <= E3max and e(g2) + e(n1) <= E3max (both P states must exist), so for the block of
    // single orbits (n1, n2) the usable part is the corner pre[n1 n2] x pre[n2 n1]. For E3max < 3 emax only these corners
    // are computed and stored ("blocked"), contiguous at boff[n1 norb + n2] with leading dimension pre[n1 norb + n2].
    std::vector<int> noff;                 // ext states with single orbit n are [noff[n], noff[n+1])
    std::vector<int> pre;
    std::vector<size_t> boff;
    size_t bsize = 0;
    bool blocked = false;
    std::vector<int> ipre;                 // [e] -> number of intl states (sorted by pair energy) with e(gab) <= E3max - e
  };

  struct Plan
  {
    const ModelSpace *ms = nullptr;
    std::vector<double> occ;
    int E3max = -1, norb = 0, J2max = 0, j2max = 0, tJmax = 0;
    int j2max3 = 0, J2max3 = 0, tJmax3 = 0;   // the same over the orbits of the 3B space (e <= emax_3body): 6j table sizes
    std::vector<std::array<int, 3>> g2;  // g -> (p, q, J), all 2B kets; those with both orbits in the 3B space first
    int ng3 = 0;                         // g < ng3: both orbits have e <= emax_3body (only these pairs occur in 3B states)
    std::vector<int> g2idx;              // (p*norb+q)*(J2max+1)+J -> g
    std::vector<int> ch3tab;             // [twoJ][par][twoTz+9] -> ch3 or -1
    std::vector<std::vector<std::array<int, 2>>> Pst;  // ch3 -> (g, c)
    std::vector<std::vector<int>> PIdx;                 // ch3 -> [g*norb + c] -> P index or -1, for g < ng3
    std::vector<std::vector<int>> PIdxT;                // ch3 -> [c*ng3 + g] -> P index or -1 (the same table, orbit-major)
    // P index of the product state (pair g, third orbit c) in 3B channel ch, or -1 (also for pairs outside the 3B space)
    int pidx(size_t ch, int g, int c) const { return (g >= 0 and g < ng3) ? PIdx[ch][(size_t)g * norb + c] : -1; }
    std::vector<arma::sp_mat> R, A;                     // ch3 -> n_P x n_k
    std::vector<arma::sp_mat> Rt;                       // ch3 -> R^T (n_k x n_P), for column-contiguous products
    std::vector<Ch> pph;
    std::vector<int> pphtab;                            // [twoJph][par][twoTz+9] -> pph channel or -1
    std::vector<double> six;                            // dense 6j table, twice-integer arguments
    std::vector<size_t> ord3, ordpph;                   // channels by decreasing cost (largest first: less tail imbalance)
    int d6[6];
    size_t i6(int a, int b, int c, int d, int e, int f) const
    {
      return ((((((size_t)a * d6[1] + b) * d6[2] + c) * d6[3] + d) * d6[4] + e) * d6[5] + f);
    }
    // bounds-checked read: arguments outside the table (orbits or pairs outside the 3B space) give 0
    double six6(int a, int b, int c, int d, int e, int f) const
    {
      if (a < 0 or b < 0 or c < 0 or d < 0 or e < 0 or f < 0 or a >= d6[0] or b >= d6[1] or c >= d6[2] or d >= d6[3] or e >= d6[4] or f >= d6[5]) return 0.0;
      return six[i6(a, b, c, d, e, f)];
    }
    double sixj(int a, int b, int c, int d, int e, int f) const { return six6(a, b, c, d, e, f); }
    size_t i3(int twoJ, int par, int twoTz) const { return ((size_t)twoJ * 2 + par) * 19 + (twoTz + 9); }
  };

  // Give freed memory back to the operating system (glibc keeps it otherwise), called after the terms with large
  // intermediates, so that the next term does not start from a fragmented, inflated resident size.
  static double t_release = 0, t_rx = 0;   // thread-safe enough for timing (called outside parallel regions)
  static void release_memory()
  {
#if defined(__GLIBC__)
    double t0 = omp_get_wtime();
    malloc_trim(0);
    t_release += omp_get_wtime() - t0;
#endif
  }
  // W = R X R^T: double precision by default; compiled with -DIMSRG_PLAN_W_FLOAT it is stored in single precision (half
  // the memory of the largest plan intermediate; results then agree with double precision to about 1e-7 relative).
#ifdef IMSRG_PLAN_W_FLOAT
  typedef arma::fmat wmat;
  static inline wmat to_wmat(arma::mat &&m) { return arma::conv_to<arma::fmat>::from(m); }
#else
  typedef arma::mat wmat;
  static inline wmat to_wmat(arma::mat &&m) { return std::move(m); }
#endif
#ifdef IMSRG_PLAN_W_FLOAT
#define PLAN_W_D(x) arma::conv_to<arma::mat>::from(x)   // a W expression used with double matrices
#else
#define PLAN_W_D(x) (x)
#endif
  static long plan_gen = 0;   // incremented whenever the plan is rebuilt (invalidates unpacked data)
  static Plan &get(const Operator &Xop, const Operator &Z)
  {
    static Plan P;
    static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    std::vector<double> occ;
    for (auto i : ms.all_orbits) occ.push_back(ms.GetOrbit(i).occ);
    if (P.ms == &ms and P.E3max == ms.GetE3max() and P.occ == occ)
      return P;
    double t0 = omp_get_wtime();
    plan_gen++;
    P = Plan();
    P.ms = &ms; P.occ = occ; P.E3max = ms.GetE3max();
    int norb = ms.GetNumberOrbits(); P.norb = norb;
    for (int i = 0; i < norb; i++) P.j2max = std::max(P.j2max, ms.GetOrbit(i).j2);
    P.J2max = P.j2max;   // max pair J = j2max (two orbits with j = j2max/2)
    P.tJmax = 2 * P.J2max + P.j2max;
    int e3b = ms.GetEMax3Body();
    for (int i = 0; i < norb; i++) { Orbit &oi = ms.GetOrbit(i); if (2 * oi.n + oi.l <= e3b) P.j2max3 = std::max(P.j2max3, oi.j2); }
    P.J2max3 = P.j2max3; P.tJmax3 = 2 * P.J2max3 + P.j2max3;
    // 2B kets. Pairs whose orbits both lie in the 3B space (e <= emax_3body) come first (g < ng3): only they occur in
    // 3B states, so the per-channel tables below are sized by ng3. For emax_3body = emax, ng3 = ng and the order is
    // the same as before.
    auto in3 = [&](size_t o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l <= e3b; };
    P.g2idx.assign((size_t)norb * norb * (P.J2max + 1), -1);
    for (int pass = 0; pass < 2; pass++)
    {
      for (size_t ch = 0; ch < ms.GetNumberTwoBodyChannels(); ch++)
      {
        auto &tbc = ms.GetTwoBodyChannel(ch);
        for (size_t k = 0; k < tbc.GetNumberKets(); k++)
        {
          auto &ket = tbc.GetKet(k);
          if ((in3(ket.p) and in3(ket.q)) != (pass == 0)) continue;
          P.g2idx[((size_t)ket.p * norb + ket.q) * (P.J2max + 1) + tbc.J] = P.g2.size();
          P.g2.push_back({(int)ket.p, (int)ket.q, tbc.J});
        }
      }
      if (pass == 0) P.ng3 = P.g2.size();
    }
    int ng3 = P.ng3;
    // 3B channel table
    size_t nch3 = ms.GetNumberThreeBodyChannels();
    P.ch3tab.assign((size_t)(P.tJmax + 2) * 2 * 19, -1);
    for (size_t ch = 0; ch < nch3; ch++)
    {
      auto &T = ms.GetThreeBodyChannel(ch);
      if (T.twoJ <= P.tJmax + 1) P.ch3tab[P.i3(T.twoJ, T.parity, T.twoTz)] = ch;
    }
    // 6j table dims (twice-integer arguments): (j2, tJ, 2J, j2, tJ, 2J)
    // all arguments come from 3B states (pph channels without ext states are dropped), so 3B-space sizes suffice
    P.d6[0] = P.j2max3 + 1; P.d6[1] = P.tJmax3 + 2; P.d6[2] = 2 * P.J2max3 + 1;
    P.d6[3] = P.j2max3 + 1; P.d6[4] = P.tJmax3 + 2; P.d6[5] = 2 * P.J2max3 + 1;
    P.six.assign((size_t)P.d6[0] * P.d6[1] * P.d6[2] * P.d6[3] * P.d6[4] * P.d6[5], 0.0);
    std::vector<char> have(P.six.size(), 0);
    auto need6 = [&](int a, int b, int c, int d, int e, int f)
    {
      size_t k = P.i6(a, b, c, d, e, f);
      if (!have[k]) { have[k] = 1; P.six[k] = AngMom::SixJ(0.5 * a, 0.5 * b, 0.5 * c, 0.5 * d, 0.5 * e, 0.5 * f); }
    };
    // P basis, R
    P.Pst.resize(nch3); P.PIdx.resize(nch3); P.PIdxT.resize(nch3); P.R.resize(nch3); P.A.resize(nch3); P.Rt.resize(nch3);
    std::vector<size_t> idx; std::vector<double> rec;
    for (size_t ch = 0; ch < nch3; ch++)
    {
      auto &T = ms.GetThreeBodyChannel(ch);
      P.PIdx[ch].assign((size_t)ng3 * norb, -1);
      for (int c = 0; c < norb; c++)   // third-orbit-major order: blocks with fixed c are contiguous (cache locality)
      {
        Orbit &oc = ms.GetOrbit(c);
        if (not in3(c)) continue;
        for (int g = 0; g < ng3; g++)
        {
          int p = P.g2[g][0], q = P.g2[g][1], J = P.g2[g][2];
          Orbit &op = ms.GetOrbit(p); Orbit &oq = ms.GetOrbit(q);
          if ((op.l + oq.l + oc.l) % 2 != T.parity or op.tz2 + oq.tz2 + oc.tz2 != T.twoTz) continue;
          if (std::abs(2 * J - oc.j2) > T.twoJ or 2 * J + oc.j2 < T.twoJ) continue;
          if (not Xop.ThreeBody.IsKetValid(J, T.twoJ, p, q, c)) continue;
          P.PIdx[ch][(size_t)g * norb + c] = P.Pst[ch].size();
          P.Pst[ch].push_back({g, c});
        }
      }
      size_t nP = P.Pst[ch].size(), nk = T.GetNumberKets();
      P.PIdxT[ch].assign((size_t)ng3 * norb, -1);
      for (int g = 0; g < ng3; g++)
        for (int c = 0; c < norb; c++) P.PIdxT[ch][(size_t)c * ng3 + g] = P.PIdx[ch][(size_t)g * norb + c];
      std::vector<arma::uword> rr, cc; std::vector<double> vv;
      for (size_t ip = 0; ip < nP; ip++)
      {
        int g = P.Pst[ch][ip][0], c = P.Pst[ch][ip][1];
        idx.clear(); rec.clear();
        size_t chr = Xop.ThreeBody.GetKetIndex_withRecoupling(P.g2[g][2], T.twoJ, P.g2[g][0], P.g2[g][1], c, idx, rec);
        if (chr != ch) continue;
        for (size_t t = 0; t < idx.size(); t++) { rr.push_back(ip); cc.push_back(idx[t]); vv.push_back(rec[t]); }
      }
      arma::umat loc(2, rr.size());
      for (size_t t = 0; t < rr.size(); t++) { loc(0, t) = rr[t]; loc(1, t) = cc[t]; }
      P.R[ch] = arma::sp_mat(true, loc, arma::vec(vv), nP, nk);
      P.Rt[ch] = P.R[ch].t();
      // A: ket (l m n; J2) -> sum over the 3 permutations (lmn),(nml),(lnm) and recoupled J2p of P states, as in
      // the back-transformation of comm333_pph_hhpss (rec_lmn * phase_45)
      std::vector<arma::uword> ar, ac; std::vector<double> av;
      double Jtot = 0.5 * T.twoJ;
      for (size_t k = 0; k < nk; k++)
      {
        auto &ket = T.GetKet(k);
        size_t l = ket.p, m = ket.q, n = ket.r;
        Orbit &ol = ms.GetOrbit(l); Orbit &om = ms.GetOrbit(m); Orbit &on = ms.GetOrbit(n);
        int J2 = ket.Jpq;
        std::array<std::array<size_t, 3>, 3> lmn = {{{l, m, n}, {n, m, l}, {l, n, m}}};
        int mn[3] = {J2, std::max(std::abs(on.j2 - om.j2), std::abs(T.twoJ - ol.j2)) / 2, std::max(std::abs(ol.j2 - on.j2), std::abs(T.twoJ - om.j2)) / 2};
        int mx[3] = {J2, std::min(on.j2 + om.j2, T.twoJ + ol.j2) / 2, std::min(ol.j2 + on.j2, T.twoJ + om.j2) / 2};
        for (int pm = 0; pm < 3; pm++)
        {
          size_t I4 = lmn[pm][0], I5 = lmn[pm][1], I6 = lmn[pm][2];
          Orbit &o4 = ms.GetOrbit(I4); Orbit &o5 = ms.GetOrbit(I5);
          for (int J2p = mn[pm]; J2p <= mx[pm]; J2p++)
          {
            if (I4 == I5 and J2p % 2 > 0) continue;
            double r = 1.0;
            if (pm == 1) r = sqrt((2 * J2 + 1) * (2 * J2p + 1)) * ms.GetSixJ(0.5 * ol.j2, 0.5 * om.j2, J2, 0.5 * on.j2, Jtot, J2p);
            if (pm == 2) r = -ms.phase((om.j2 + on.j2) / 2 + J2 + J2p) * sqrt((2 * J2 + 1) * (2 * J2p + 1)) * ms.GetSixJ(0.5 * om.j2, 0.5 * ol.j2, J2, 0.5 * on.j2, Jtot, J2p);
            int ph45 = I4 > I5 ? -ms.phase((o4.j2 + o5.j2) / 2 - J2p) : 1;
            int g = P.g2idx[((size_t)std::min(I4, I5) * norb + std::max(I4, I5)) * (P.J2max + 1) + J2p];
            if (g < 0) continue;
            int ip = P.pidx(ch, g, I6);
            if (ip < 0) continue;
            ar.push_back(ip); ac.push_back(k); av.push_back(r * ph45);
          }
        }
      }
      arma::umat aloc(2, ar.size());
      for (size_t t = 0; t < ar.size(); t++) { aloc(0, t) = ar[t]; aloc(1, t) = ac[t]; }
      P.A[ch] = arma::sp_mat(true, aloc, arma::vec(av), nP, nk);   // duplicates are summed
    }
    auto eorb = [&](int o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l; };
    auto epair = [&](int g) { return eorb(P.g2[g][0]) + eorb(P.g2[g][1]); };
    // pph channels (same construction as comm333_pph_hhpss; default cuts only)
    P.pphtab.assign((size_t)(P.tJmax + 2) * 2 * 19, -1);
    for (int twoJph = 1; twoJph <= P.tJmax; twoJph += 2)
      for (int par = 0; par <= 1; par++)
        for (int twoTz = -3; twoTz <= 3; twoTz += 2)
        {
          Ch C; C.twoJph = twoJph; C.par = par; C.twoTz = twoTz;
          // ext and intl states (pair g, orbit n) are only read through product states P(g, .) and P(., n), so only pairs
          // and orbits of the 3B space are needed
          C.extIdx.assign((size_t)ng3 * norb, -1);
          for (int n = 0; n < norb; n++)   // single-orbit-major order (cache locality in the gathers), pair energy within n
          {
            Orbit &on = ms.GetOrbit(n);
            if (not in3(n)) continue;
            std::vector<int> gs;
            for (int g = 0; g < ng3; g++)
            {
              int i = P.g2[g][0], j = P.g2[g][1], Jij = P.g2[g][2];
              Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
              if ((on.l + oi.l + oj.l) % 2 != par) continue;
              if ((oi.tz2 + oj.tz2 - on.tz2) != twoTz) continue;
              if (std::abs(2 * Jij - on.j2) > twoJph or 2 * Jij + on.j2 < twoJph) continue;
              gs.push_back(g);
            }
            std::stable_sort(gs.begin(), gs.end(), [&](int a, int b) { return epair(a) < epair(b); });
            for (int g : gs)
            {
              int i = P.g2[g][0], j = P.g2[g][1];
              Orbit &oi = ms.GetOrbit(i); Orbit &oj = ms.GetOrbit(j);
              C.extIdx[(size_t)g * norb + n] = C.ext.size();
              C.ext.push_back({g, n});
              double occf = oi.occ * oj.occ * (1 - on.occ) + (1 - oi.occ) * (1 - oj.occ) * on.occ;
              if (i == j) occf *= 0.5;
              if (std::abs(occf) > 1e-8) { C.intl.push_back({g, n}); C.occ_int.push_back(occf); }
            }
          }
          if (C.ext.empty()) continue;   // channels without intl states are kept: column channels for channel-changing Y
          C.noff.assign(norb + 1, 0);
          for (auto &e : C.ext) C.noff[e[1] + 1]++;
          for (int n = 0; n < norb; n++) C.noff[n + 1] += C.noff[n];
          C.pre.assign((size_t)norb * norb, 0); C.boff.assign((size_t)norb * norb, 0);
          for (int n1 = 0; n1 < norb; n1++)
            for (int n2 = 0; n2 < norb; n2++)
            {
              int lim = P.E3max - eorb(n2), k = C.noff[n1];
              while (k < C.noff[n1 + 1] and epair(C.ext[k][0]) <= lim) k++;
              C.pre[(size_t)n1 * norb + n2] = k - C.noff[n1];
            }
          size_t used = 0;
          for (int n1 = 0; n1 < norb; n1++)
            for (int n2 = 0; n2 < norb; n2++)
            {
              C.boff[(size_t)n1 * norb + n2] = used;
              used += (size_t)C.pre[(size_t)n1 * norb + n2] * C.pre[(size_t)n2 * norb + n1];
            }
          C.bsize = used;
          C.blocked = used < 0.9 * (double)C.ext.size() * C.ext.size();
          {
            // intl sorted by pair energy: for the block (n1, n2) only intl (gab, c) with e(gab) <= E3max - max(e(n1), e(n2))
            // contribute (P(gab, n) must exist), a prefix of length ipre[max(e(n1), e(n2))]
            std::vector<size_t> io(C.intl.size());
            for (size_t t = 0; t < io.size(); t++) io[t] = t;
            std::stable_sort(io.begin(), io.end(), [&](size_t a, size_t b) { return epair(C.intl[a][0]) < epair(C.intl[b][0]); });
            std::vector<std::array<int, 2>> in2; std::vector<double> oc2;
            for (size_t t : io) { in2.push_back(C.intl[t]); oc2.push_back(C.occ_int[t]); }
            C.intl = std::move(in2); C.occ_int = std::move(oc2);
            int emx = 0;
            for (int o = 0; o < norb; o++) emx = std::max(emx, eorb(o));
            C.ipre.assign(emx + 1, 0);
            for (int b = 0; b <= emx; b++)
            {
              int k = 0;
              while (k < (int)C.intl.size() and epair(C.intl[k][0]) <= P.E3max - b) k++;
              C.ipre[b] = k;
            }
          }
          P.pphtab[P.i3(twoJph, par, twoTz)] = P.pph.size();
          P.pph.push_back(std::move(C));
        }
    // the 6j values that occur
    for (auto &C : P.pph)
      for (auto &e : C.ext)
        for (auto &a : C.intl)
        {
          int jn = ms.GetOrbit(e[1]).j2, Jij = P.g2[e[0]][2], jc = ms.GetOrbit(a[1]).j2, Jab = P.g2[a[0]][2];
          for (int tJp = std::max(std::abs(2 * Jij - jc), std::abs(2 * Jab - jn)); tJp <= std::min(2 * Jij + jc, 2 * Jab + jn); tJp += 2)
            need6(jn, C.twoJph, 2 * Jij, jc, tJp, 2 * Jab);
        }
    for (size_t ch = 0; ch < nch3; ch++)
    {
      int tJ = ms.GetThreeBodyChannel(ch).twoJ;
      if (tJ > P.tJmax) continue;
      std::vector<std::array<int, 2>> types;   // distinct (J12, j3)
      std::set<std::array<int, 2>> seen;
      for (auto &s : P.Pst[ch]) { std::array<int, 2> t = {P.g2[s[0]][2], ms.GetOrbit(s[1]).j2}; if (seen.insert(t).second) types.push_back(t); }
      for (auto &t1 : types)
        for (auto &t2 : types)
          for (int tJph = 1; tJph <= P.tJmax; tJph += 2)
            if (std::abs(2 * t1[0] - t2[1]) <= tJph and tJph <= 2 * t1[0] + t2[1] and std::abs(2 * t2[0] - t1[1]) <= tJph and tJph <= 2 * t2[0] + t1[1])
              need6(t1[1], tJ, 2 * t1[0], t2[1], tJph, 2 * t2[0]);
    }
    Z.profiler.timer["_comm333_pph_plan_build"] += omp_get_wtime() - t0;
    P.ord3.resize(nch3); P.ordpph.resize(P.pph.size());
    for (size_t ch = 0; ch < nch3; ch++) P.ord3[ch] = ch;
    for (size_t c = 0; c < P.pph.size(); c++) P.ordpph[c] = c;
    std::stable_sort(P.ord3.begin(), P.ord3.end(), [&](size_t a, size_t b) { return P.Pst[a].size() > P.Pst[b].size(); });
    auto cpph = [&](size_t c) { double ne = P.pph[c].ext.size(), ni = P.pph[c].intl.size(); return ne * ne * ni; };
    std::stable_sort(P.ordpph.begin(), P.ordpph.end(), [&](size_t a, size_t b) { return cpph(a) > cpph(b); });
    if (std::getenv("IMSRG_PLAN_DEBUG") != nullptr)
    {
      double gi = 0, gR = 0, gx = 0;   // GB: index tables, sparse R and A, pph tables
      for (size_t ch = 0; ch < nch3; ch++)
      {
        gi += 4e-9 * (P.PIdx[ch].size() + P.PIdxT[ch].size()) + 8e-9 * P.Pst[ch].size();
        gR += 12e-9 * (2.0 * P.R[ch].n_nonzero + P.A[ch].n_nonzero);
      }
      for (auto &C : P.pph) gx += 4e-9 * C.extIdx.size() + 8e-9 * (C.ext.size() + C.intl.size() + C.pre.size() + C.boff.size()) + 8e-9 * C.occ_int.size();
      std::cout << "PLANDEBUG plan memory (GB): index tables " << gi << ", R/Rt/A " << gR << ", pph " << gx << ", 6j " << 8e-9 * P.six.size() << std::endl;
      double tp = 0, t3 = 0;
      for (size_t c = 0; c < P.pph.size(); c++) tp += cpph(c);
      for (size_t ch = 0; ch < nch3; ch++) t3 += (double)P.Pst[ch].size() * P.Pst[ch].size();
      std::cout << "PLANDEBUG pph channels " << P.pph.size() << " cost share of top 5:";
      for (size_t r = 0; r < std::min<size_t>(5, P.pph.size()); r++) std::cout << " " << cpph(P.ordpph[r]) / tp << " (ne " << P.pph[P.ordpph[r]].ext.size() << " ni " << P.pph[P.ordpph[r]].intl.size() << ")";
      std::cout << std::endl << "PLANDEBUG ch3 " << nch3 << " nP^2 share of top 5:";
      for (size_t r = 0; r < std::min<size_t>(5, nch3); r++) std::cout << " " << (double)P.Pst[P.ord3[r]].size() * P.Pst[P.ord3[r]].size() / t3 << " (nP " << P.Pst[P.ord3[r]].size() << ")";
      std::cout << std::endl;
    }
    return P;
  }

  // Ab^T G with column-contiguous access: out(k, l) = Ab.col(k) . G.col(l)
  static arma::mat AtM(const arma::sp_mat &Ab, const arma::mat &G)
  {
    size_t nb = Ab.n_cols, nl = G.n_cols;
    arma::mat out(nb, nl);
    for (size_t l = 0; l < nl; l++)
    {
      const double *gc = G.colptr(l);
      double *oc = out.colptr(l);
      for (size_t k = 0; k < nb; k++)
      {
        double s = 0;
        for (arma::uword t = Ab.col_ptrs[k]; t < Ab.col_ptrs[k + 1]; t++) s += Ab.values[t] * gc[Ab.row_indices[t]];
        oc[k] = s;
      }
    }
    return out;
  }

  // M S for dense M and CSC sparse S, as column axpys: out.col(j) = sum_{(i, v) in S.col(j)} v M.col(i)
  static arma::mat DS(const arma::mat &M, const arma::sp_mat &S)
  {
    size_t nr = M.n_rows, nc = S.n_cols;
    arma::mat out(nr, nc, arma::fill::zeros);
    for (size_t j = 0; j < nc; j++)
    {
      double *oc = out.colptr(j);
      for (arma::uword t = S.col_ptrs[j]; t < S.col_ptrs[j + 1]; t++)
      {
        const double v = S.values[t];
        const double *mc = M.colptr(S.row_indices[t]);
        for (size_t i = 0; i < nr; i++) oc[i] += v * mc[i];
      }
    }
    return out;
  }

  // S D for CSC sparse S and dense D, column by column (armadillo's sparse x dense walks the rows of D)
  static arma::mat SD(const arma::sp_mat &S, const arma::mat &D)
  {
    size_t nr = S.n_rows, nc = D.n_cols, nm = S.n_cols;
    arma::mat out(nr, nc, arma::fill::zeros);
    for (size_t j = 0; j < nc; j++)
    {
      const double *dc = D.colptr(j);
      double *oc = out.colptr(j);
      for (size_t c = 0; c < nm; c++)
      {
        const double d = dc[c];
        if (d == 0) continue;
        for (arma::uword t = S.col_ptrs[c]; t < S.col_ptrs[c + 1]; t++) oc[S.row_indices[t]] += S.values[t] * d;
      }
    }
    return out;
  }

  // Zk = Ab^T WZ Ak with column-contiguous access only (armadillo's sparse^T x dense walks dense rows)
  static arma::mat AtWA(const arma::sp_mat &Ab, const arma::mat &WZ, const arma::sp_mat &Ak) { return AtM(Ab, DS(WZ, Ak)); }


  // dE3max / occnat3 cut, as in the IMSRG(3) terms of IMSRG3Commutators.cc: a 3B state |abc> takes part (as a bra, a ket
  // or an intermediate) only if
  //     |e_a - eF| + |e_b - eF| + |e_c - eF| <= dE3max      and      prod_{x = a,b,c} occnat_x (1 - occnat_x) >= occnat3cut,
  // with e = 2n + l and the Fermi energy eF of the orbit's tz. The cut depends only on the three orbits, so it is the same
  // for a product state |(ab)J c> and for a stored ket with the same orbits, and it commutes with the antisymmetriser.
  // The plan therefore masks each unpacked operator (stored kets failing the cut -> zero rows and columns of X, hence of
  // R X and W) and each 3B output block. keep[ch][k] = 1 if stored ket k of 3B channel ch passes.
  struct Cut3
  {
    bool active = false;
    const ModelSpace *ms = nullptr; double dE3max = 0, occcut = 0; std::vector<double> occnat; std::map<int, double> ef;
    std::vector<std::vector<char>> keep;
  };
  static Cut3 &cut3(ModelSpace &ms)
  {
    static Cut3 C; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    std::vector<double> occnat;
    for (auto i : ms.all_orbits) occnat.push_back(ms.GetOrbit(i).occ_nat);
    double de3 = ms.GetdE3max(), oc = ms.GetOccNat3Cut();
    std::map<int, double> ef = ms.GetEFermi();
    if (C.ms == &ms and C.dE3max == de3 and C.occcut == oc and C.occnat == occnat and C.ef == ef and C.keep.size() == ms.GetNumberThreeBodyChannels())
      return C;
    C = Cut3(); C.ms = &ms; C.dE3max = de3; C.occcut = oc; C.occnat = occnat; C.ef = ef;
    C.active = (de3 < 900 or oc > 0);
    size_t nch3 = ms.GetNumberThreeBodyChannels();
    C.keep.resize(nch3);
    if (not C.active) return C;
    auto de = [&](Orbit *o) { return std::abs(2 * o->n + o->l - ef[o->tz2]); };
    auto oo = [&](Orbit *o) { return o->occ_nat * (1 - o->occ_nat); };
    for (size_t ch = 0; ch < nch3; ch++)
    {
      auto &T = ms.GetThreeBodyChannel(ch);
      size_t nk = T.GetNumberKets();
      C.keep[ch].assign(nk, 1);
      for (size_t k = 0; k < nk; k++)
      {
        Ket3 &kt = T.GetKet(k);
        if (de(kt.op) + de(kt.oq) + de(kt.oR) > de3 or oo(kt.op) * oo(kt.oq) * oo(kt.oR) < oc) C.keep[ch][k] = 0;
      }
    }
    return C;
  }
  // zero the rows (bra channel cb) and columns (ket channel ck) of kets that fail the cut
  static void mask3(const Cut3 &C, size_t cb, size_t ck, arma::mat &M)
  {
    if (not C.active) return;
    const auto &kb = C.keep[cb], &kk = C.keep[ck];
    for (size_t j = 0; j < M.n_cols; j++)
    {
      double *mc = M.colptr(j);
      if (not kk[j]) { std::fill(mc, mc + M.n_rows, 0.0); continue; }
      for (size_t i = 0; i < M.n_rows; i++) if (not kb[i]) mc[i] = 0;
    }
  }
  static bool pass3(const Cut3 &C, size_t ch, size_t k) { return not C.active or C.keep[ch][k]; }

  // Unpacked W = R X R^T, R X and X per 3B channel. Shared without any check between the terms of one CommutatorScalarScalar
  // call: inside the scope perfplan_begin_commutator() ... perfplan_end_commutator() X and Y do not change. Across calls
  // (and for direct calls of single terms) the stored matrix of each channel is compared bitwise with the cached one and
  // only changed channels are recomputed. Two slots, least-recently-used replacement, so the two operators of one kernel
  // call never evict each other.
  static long epoch = 0;
  static bool scope_active = false;
  static long use_counter = 0;
  struct WCache
  {
    long epoch = -1; const void *key = nullptr; long last_use = -1;
    const Plan *plan = nullptr; long plan_gen = -1;   // the plan the data were unpacked with
    std::vector<wmat> W;                              // R X R^T
    std::vector<arma::mat> K;                         // stored-ket matrix (R X is computed on demand, see unpack_rx)
  };
  static WCache wslot[2];
  static std::mutex wmtx;
  static WCache &unpack_slot(const Operator &op, Plan &P);
  static const std::vector<wmat> &unpack(const Operator &op, Plan &P) { return unpack_slot(op, P).W; }
  static std::vector<arma::mat> unpack_rx(const Operator &op, Plan &P);   // R X (n_P x n_k), computed on demand
  static const std::vector<arma::mat> &unpack_k(const Operator &op, Plan &P) { return unpack_slot(op, P).K; }     // X (n_k x n_k)
  static WCache &unpack_slot(const Operator &op, Plan &P)
  {
    std::lock_guard<std::mutex> lock(wmtx);
    const void *key = &op.ThreeBody;
    if (scope_active)
      for (auto &c : wslot)
        if (c.epoch == epoch and c.key == key) { c.last_use = ++use_counter; return c; }
    // Not unpacked yet in this commutator. Take the slot last used by the same operator object if there is one (e.g. Omega
    // in the nested commutators of a BCH series), else the least recently used one. Channels whose stored matrix is
    // bitwise equal to the cached one keep their R X and W, which depend on it alone, so the reuse never changes results.
    double t0 = omp_get_wtime();
    WCache *same = nullptr;
    for (auto &c : wslot) if (c.key == key) same = &c;
    WCache &c = same ? *same : ((wslot[0].last_use <= wslot[1].last_use) ? wslot[0] : wslot[1]);
    ModelSpace &ms = *op.modelspace;
    size_t nch3 = ms.GetNumberThreeBodyChannels();
    const Cut3 &cut = cut3(ms);
    bool reuse = c.plan == &P and c.plan_gen == plan_gen and c.K.size() == nch3;
    if (not reuse) { c.W.assign(nch3, wmat()); c.K.assign(nch3, arma::mat()); }
    auto &T3 = op.ThreeBody;
    double nreused = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : nreused)
    for (size_t ich = 0; ich < nch3; ich++)
    {
      size_t ch = P.ord3[ich];
      size_t nk = ms.GetThreeBodyChannel(ch).GetNumberKets();
      if (P.Pst[ch].empty()) continue;
      if (nk == 0)   // P states that all vanish on antisymmetrisation (possible for E3max < 3 emax): W = 0
      {
        size_t nP = P.Pst[ch].size();
        if (c.W[ch].n_rows != nP) { c.W[ch].zeros(nP, nP); c.K[ch].reset(); }
        continue;
      }
      arma::mat Xk(nk, nk);
      T3.GetChannelMatrix_pn(ch, Xk);   // = GetME_pn_ch(ch, ch, i, k), one channel lookup
      mask3(cut, ch, ch, Xk);           // dE3max / occnat3 cut (no-op without cuts)
      if (reuse and c.K[ch].n_rows == nk and c.K[ch].n_cols == nk
          and std::memcmp(c.K[ch].memptr(), Xk.memptr(), nk * nk * sizeof(double)) == 0)
      {
        nreused += 1;
        continue;
      }
      arma::mat RXc = AtM(P.Rt[ch], Xk);   // R X
      c.W[ch] = to_wmat(DS(RXc, P.Rt[ch]));   // (R X) R^T
      c.K[ch] = std::move(Xk);
    }
    c.epoch = epoch; c.key = key; c.last_use = ++use_counter; c.plan = &P; c.plan_gen = plan_gen;
    if (std::getenv("IMSRG_PLAN_DEBUG") != nullptr and not reuse)
    {
      double gW = 0, gK = 0;
      for (size_t ch = 0; ch < nch3; ch++) { gW += 8e-9 * c.W[ch].n_elem; gK += 8e-9 * c.K[ch].n_elem; }
      std::cout << "PLANDEBUG unpack cache (GB) of one operator: W " << gW << ", K " << gK << std::endl;
    }
    op.profiler.timer["_perfplan_unpack"] += omp_get_wtime() - t0;
    op.profiler.timer["_perfplan_unpack_reused_channels"] += nreused;
    return c;
  }
  // R X per channel (n_P x n_k), from the cached stored-ket blocks K. Not cached: the few terms that read it hold it only
  // while they run, so it is not resident during the memory-heavy terms. Cheap (a sparse x dense product per channel).
  static std::vector<arma::mat> unpack_rx(const Operator &op, Plan &P)
  {
    WCache &c = unpack_slot(op, P);
    double t0 = omp_get_wtime();
    size_t nch3 = c.K.size();
    std::vector<arma::mat> RX(nch3);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nch3; ich++)
    {
      size_t ch = P.ord3[ich];
      if (P.Pst[ch].empty()) continue;
      if (c.K[ch].n_elem == 0) { RX[ch].zeros(P.Pst[ch].size(), 0); continue; }   // no stored kets
      RX[ch] = AtM(P.Rt[ch], c.K[ch]);
    }
    t_rx += omp_get_wtime() - t0;
    op.profiler.timer["_perfplan_rx_on_demand"] = t_rx;
    return RX;
  }
} // namespace perfplan

void perfplan_begin_commutator() { std::lock_guard<std::mutex> lock(perfplan::wmtx); perfplan::epoch++; perfplan::scope_active = true; }
void perfplan_end_commutator() { std::lock_guard<std::mutex> lock(perfplan::wmtx); perfplan::epoch++; perfplan::scope_active = false; }

void comm333_pph_hhpss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  auto &X3 = X.ThreeBody; auto &Y3 = Y.ThreeBody; auto &Z3 = Z.ThreeBody;
  if (not(X3.IsAllocated() and Y3.IsAllocated())) return;
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0)
  {
    comm333_pph_hhpss(X, Y, Z);
    return;
  }
  int hX = X.IsHermitian() ? 1 : -1;
  int hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(X, Z);
  ModelSpace &ms = *Z.modelspace;
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  int norb = P.norb;
  const perfplan::Cut3 &cutZ = perfplan::cut3(ms);
  double t = omp_get_wtime();
  // Step 1 (unpack): W_X = R X R^T, W_Y = R Y R^T per 3B channel (cached, see the reader's guide at the top)
  const std::vector<perfplan::wmat> &WX = perfplan::unpack(X, P);
  const std::vector<perfplan::wmat> &WY = perfplan::unpack(Y, P);
  Z.profiler.timer["_comm333_pph_plan_unpack"] += omp_get_wtime() - t; t = omp_get_wtime();
  // Steps 2 + 3, per pph channel: fill Xbar, Ybar from W_X, W_Y (6j-weighted gathers), then Zbar = hY Xbar Ybar^T - (transpose)
  size_t npph = P.pph.size();
  std::vector<arma::mat> Zbar(npph);
  double tt_fill = 0, tt_gemm = 0;   // summed thread times
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : tt_fill, tt_gemm)
  for (size_t ic = 0; ic < npph; ic++)
  {
    size_t c = P.ordpph[ic];
    double tc0 = omp_get_wtime();
    auto &C = P.pph[c];
    size_t ne = C.ext.size(), ni = C.intl.size();
    if (ni == 0) continue;   // Zbar = 0 (left empty)
    // Xbar from X and Ybar from Y, both (ext x int); Zbar = hY Xbar D Ybar^T with D = occ_int, antisymmetrised below
    arma::mat Xbar(ne, ni, arma::fill::zeros), Ybar(ne, ni, arma::fill::zeros);
    int nJ = P.J2max + 1;
    std::vector<int> rJ, cJ; std::vector<double> wt(nJ * nJ); std::vector<char> haveR(nJ), haveC(nJ);
    std::vector<arma::uword> er, pr, ac2, pc;   // reused per block (capacity kept)
    // Blocked by the single orbits: ext states (g12, n) with the same n and int states (gab, c) with the same c share
    // parity/Tz, so for each J' the block is a gathered submatrix WX(P(g12, c), P(gab, n)) times block-constant weights.
    std::map<int, std::vector<arma::uword>> eByN, aByC;
    std::vector<int> epx(ne), epa(ni);   // pair energies of the ext and intl states
    auto eo = [&](int o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l; };
    for (size_t e = 0; e < ne; e++) epx[e] = eo(P.g2[C.ext[e][0]][0]) + eo(P.g2[C.ext[e][0]][1]);
    for (size_t a = 0; a < ni; a++) epa[a] = eo(P.g2[C.intl[a][0]][0]) + eo(P.g2[C.intl[a][0]][1]);
    for (size_t e = 0; e < ne; e++) eByN[C.ext[e][1]].push_back(e);
    for (size_t a = 0; a < ni; a++) aByC[C.intl[a][1]].push_back(a);
    for (auto &en : eByN)
    {
      int n = en.first; auto &E = en.second;
      Orbit &on = ms.GetOrbit(n);
      int g0 = C.ext[E[0]][0];
      Orbit &oi = ms.GetOrbit(P.g2[g0][0]); Orbit &oj = ms.GetOrbit(P.g2[g0][1]);
      for (auto &ac : aByC)
      {
        int cc = ac.first; auto &Av = ac.second;
        Orbit &oc = ms.GetOrbit(cc);
        int par = (oi.l + oj.l + oc.l) % 2, tz = oi.tz2 + oj.tz2 + oc.tz2;
        for (int tJp = 1; tJp <= P.tJmax; tJp += 2)
        {
          int ch3 = P.ch3tab[P.i3(tJp, par, tz)];
          if (ch3 < 0) continue;
          er.clear(); pr.clear(); ac2.clear(); pc.clear();
          // ext rows of orbit n and intl columns of orbit cc are sorted by pair energy: P(g12, cc) and P(gab, n) exist only
          // below E3max - e(cc) and E3max - e(n), so the scans stop there
          const int limE = P.E3max - (2 * oc.n + oc.l), limA = P.E3max - (2 * on.n + on.l);
          for (auto e : E)
          {
            if (C.blocked and epx[e] > limE) break;
            int g12 = C.ext[e][0], Jij = P.g2[g12][2];
            if (std::abs(2 * Jij - oc.j2) > tJp or 2 * Jij + oc.j2 < tJp) continue;
            int p1 = P.pidx(ch3, g12, cc);
            if (p1 >= 0) { er.push_back(e); pr.push_back(p1); }
          }
          if (er.empty()) continue;
          for (auto a : Av)
          {
            if (C.blocked and epa[a] > limA) break;
            int gab = C.intl[a][0], Jab = P.g2[gab][2];
            if (std::abs(2 * Jab - on.j2) > tJp or 2 * Jab + on.j2 < tJp) continue;
            int p2 = P.pidx(ch3, gab, n);
            if (p2 >= 0) { ac2.push_back(a); pc.push_back(p2); }
          }
          if (ac2.empty()) continue;
          // weight (2J'+1) 6j(jn Jph J12; jc J' Jab) depends only on (J12, Jab): small table, fused gather-scatter
          rJ.resize(er.size()); cJ.resize(ac2.size());
          std::fill(haveR.begin(), haveR.end(), 0); std::fill(haveC.begin(), haveC.end(), 0);
          for (size_t r = 0; r < er.size(); r++) { rJ[r] = P.g2[C.ext[er[r]][0]][2]; haveR[rJ[r]] = 1; }
          for (size_t q = 0; q < ac2.size(); q++) { cJ[q] = P.g2[C.intl[ac2[q]][0]][2]; haveC[cJ[q]] = 1; }
          for (int J12 = 0; J12 < nJ; J12++)
            if (haveR[J12])
              for (int Jab = 0; Jab < nJ; Jab++)
                if (haveC[Jab]) wt[J12 * nJ + Jab] = (tJp + 1) * P.sixj(on.j2, C.twoJph, 2 * J12, oc.j2, tJp, 2 * Jab);
          const perfplan::wmat &WXc = WX[ch3], &WYc = WY[ch3];
          for (size_t q = 0; q < ac2.size(); q++)
          {
            const auto *wxc = WXc.colptr(pc[q]); const auto *wyc = WYc.colptr(pc[q]);
            double *xb = Xbar.colptr(ac2[q]), *yb = Ybar.colptr(ac2[q]);
            const double *wq = &wt[cJ[q]];
            for (size_t r = 0; r < er.size(); r++)
            {
              double w = wq[rJ[r] * nJ];
              xb[er[r]] += w * wxc[pr[r]];
              yb[er[r]] += w * wyc[pr[r]];
            }
          }
        }
      }
    }
    for (size_t a = 0; a < ni; a++) Ybar.col(a) *= C.occ_int[a];
    double tc1 = omp_get_wtime();
    if (not C.blocked)
    {
      Zbar[c] = hY * Xbar * Ybar.t();   // (dsyr2k for the symmetric hX hY = -1 case is no faster in OpenBLAS: same flops)
      Zbar[c] -= hX * hY * Zbar[c].t();
    }
    else
    {
      // only the usable corners: with M = Xbar Ybar^T, Zbar_{12} = hY M_{12} - hX M_{21}^T, Zbar_{21} = hY M_{21} - hX M_{12}^T.
      // M_{12} and M_{21} are written straight into the two stored blocks (BLAS on the row ranges, inner dimension cut to
      // the contributing intl prefix), then combined pairwise in place.
      Zbar[c].set_size(C.bsize, 1);
      double *zb = Zbar[c].memptr();
      const double one = 1.0, zero = 0.0;
      const arma::blas_int ldx = ne;
      auto eo = [&](int o) { Orbit &oo = ms.GetOrbit(o); return 2 * oo.n + oo.l; };
      for (int n1 = 0; n1 < norb; n1++)
        for (int n2 = n1; n2 < norb; n2++)
        {
          int m = C.pre[(size_t)n1 * norb + n2], k = C.pre[(size_t)n2 * norb + n1];
          if (m == 0 or k == 0) continue;
          int ka = C.ipre[std::max(eo(n1), eo(n2))];
          double *b12 = zb + C.boff[(size_t)n1 * norb + n2], *b21 = zb + C.boff[(size_t)n2 * norb + n1];
          int r1 = C.noff[n1], r2 = C.noff[n2];
          if (ka == 0)
          {
            std::fill(b12, b12 + (size_t)m * k, 0.0);
            if (n1 != n2) std::fill(b21, b21 + (size_t)m * k, 0.0);
            continue;
          }
          arma::blas_int bm = m, bk = k, bka = ka;
          arma::blas::gemm<double>("N", "T", &bm, &bk, &bka, &one, Xbar.memptr() + r1, &ldx, Ybar.memptr() + r2, &ldx, &zero, b12, &bm);
          if (n1 == n2)
          {
            for (int q = 0; q < k; q++)
            {
              b12[q + (size_t)q * m] *= (hY - hX);
              for (int r = q + 1; r < m; r++)
              {
                double a = b12[r + (size_t)q * m], b = b12[q + (size_t)r * m];
                b12[r + (size_t)q * m] = hY * a - hX * b; b12[q + (size_t)r * m] = hY * b - hX * a;
              }
            }
            continue;
          }
          arma::blas::gemm<double>("N", "T", &bk, &bm, &bka, &one, Xbar.memptr() + r2, &ldx, Ybar.memptr() + r1, &ldx, &zero, b21, &bk);
          for (int q = 0; q < k; q++)
            for (int r = 0; r < m; r++)
            {
              double a = b12[r + (size_t)q * m], b = b21[q + (size_t)r * k];
              b12[r + (size_t)q * m] = hY * a - hX * b; b21[q + (size_t)r * k] = hY * b - hX * a;
            }
        }
    }
    tt_fill += tc1 - tc0; tt_gemm += omp_get_wtime() - tc1;
  }
  Z.profiler.timer["_comm333_pph_plan_threadtime_fill"] += tt_fill;
  Z.profiler.timer["_comm333_pph_plan_threadtime_gemm"] += tt_gemm;
  Z.profiler.timer["_comm333_pph_plan_fill_gemm"] += omp_get_wtime() - t; t = omp_get_wtime();
  // Step 4 (back): W_Z(P(g12, I3), P(g45, I6)) = sum_Jph (2Jph+1) {j3 J J12; j6 Jph J45} Zbar^{Jph}(e(g12, I6), e(g45, I3)),
  // Step 5: Z += A^T W_Z A  (both per 3B channel)
  double tt_wz = 0, tt_a = 0, tt_add = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : tt_wz, tt_a, tt_add)
  for (size_t ich = 0; ich < nch3; ich++)
  {
    size_t ch = P.ord3[ich];
    auto &T = ms.GetThreeBodyChannel(ch);
    size_t nP = P.Pst[ch].size(), nk = T.GetNumberKets();
    if (nP == 0 or nk == 0) continue;
    double tb0 = omp_get_wtime();
    arma::mat WZ(nP, nP, arma::fill::zeros);
    // Blocked by the third orbit: for fixed (I3, I6) all rows (g12, I3) share parity/Tz with I6, so for each Jph the
    // block is a gathered submatrix Zbar(ext(g12, I6), ext(g45, I3)) times (2Jph+1) 6j(j3 J J12; j6 Jph J45).
    std::vector<std::vector<arma::uword>> byOrb(norb);
    for (size_t p = 0; p < nP; p++) byOrb[P.Pst[ch][p][1]].push_back(p);
    std::vector<int> rJ, cJ;
    std::vector<arma::uword> rp, re, cpos, ce;   // reused per block (capacity kept)
    std::vector<double> wt((P.J2max + 1) * (P.J2max + 1));
    std::vector<char> haveR(P.J2max + 1), haveC(P.J2max + 1);
    for (int I3 = 0; I3 < norb; I3++)
    {
      auto &Rw = byOrb[I3];
      if (Rw.empty()) continue;
      Orbit &o3 = ms.GetOrbit(I3);
      int g0 = P.Pst[ch][Rw[0]][0];
      Orbit &o1 = ms.GetOrbit(P.g2[g0][0]); Orbit &o2 = ms.GetOrbit(P.g2[g0][1]);
      for (int I6 = 0; I6 < norb; I6++)
      {
        auto &Cl = byOrb[I6];
        if (Cl.empty()) continue;
        Orbit &o6 = ms.GetOrbit(I6);
        int par = (o1.l + o2.l + o6.l) % 2, tz = o1.tz2 + o2.tz2 - o6.tz2;
        for (int tJph = 1; tJph <= P.tJmax; tJph += 2)
        {
          int cp = P.pphtab[P.i3(tJph, par, tz)];
          if (cp < 0 or Zbar[cp].n_elem == 0) continue;
          auto &C = P.pph[cp];
          rp.clear(); re.clear(); cpos.clear(); ce.clear();
          for (auto p1 : Rw)
          {
            int g12 = P.Pst[ch][p1][0];
            if (std::abs(2 * P.g2[g12][2] - o6.j2) > tJph or 2 * P.g2[g12][2] + o6.j2 < tJph) continue;
            int e1 = g12 < P.ng3 ? C.extIdx[(size_t)g12 * norb + I6] : -1;
            if (e1 >= 0) { rp.push_back(p1); re.push_back(e1); }
          }
          if (rp.empty()) continue;
          for (auto p2 : Cl)
          {
            int g45 = P.Pst[ch][p2][0];
            if (std::abs(2 * P.g2[g45][2] - o3.j2) > tJph or 2 * P.g2[g45][2] + o3.j2 < tJph) continue;
            int e2 = g45 < P.ng3 ? C.extIdx[(size_t)g45 * norb + I3] : -1;
            if (e2 >= 0) { cpos.push_back(p2); ce.push_back(e2); }
          }
          if (cpos.empty()) continue;
          // the weight (2Jph+1) 6j(j3 J J12; j6 Jph J45) depends only on (J12, J45): small table, fused gather-scatter
          int nJ = P.J2max + 1;
          rJ.resize(rp.size()); cJ.resize(cpos.size());
          std::fill(wt.begin(), wt.end(), 0.0);
          std::fill(haveR.begin(), haveR.end(), 0); std::fill(haveC.begin(), haveC.end(), 0);
          for (size_t r = 0; r < rp.size(); r++) { rJ[r] = P.g2[P.Pst[ch][rp[r]][0]][2]; haveR[rJ[r]] = 1; }
          for (size_t q = 0; q < cpos.size(); q++) { cJ[q] = P.g2[P.Pst[ch][cpos[q]][0]][2]; haveC[cJ[q]] = 1; }
          for (int J12 = 0; J12 < nJ; J12++)
            if (haveR[J12])
              for (int J45 = 0; J45 < nJ; J45++)
                if (haveC[J45]) wt[J12 * nJ + J45] = (tJph + 1) * P.sixj(o3.j2, T.twoJ, 2 * J12, o6.j2, tJph, 2 * J45);
          const arma::mat &Zb = Zbar[cp];
          const double *zbase = Zb.memptr(); size_t ldz = Zb.n_rows;
          if (C.blocked)
          {
            // block (I6, I3), local positions; used entries lie inside the stored corner (both P states exist)
            size_t b = (size_t)I6 * norb + I3;
            ldz = C.pre[b];
            zbase = Zb.memptr() + C.boff[b] - C.noff[I6] - (std::ptrdiff_t)C.noff[I3] * (std::ptrdiff_t)ldz;
          }
          for (size_t q = 0; q < cpos.size(); q++)
          {
            const double *zc = zbase + (std::ptrdiff_t)ce[q] * (std::ptrdiff_t)ldz;
            double *wc = WZ.colptr(cpos[q]);
            const double *wq = &wt[cJ[q]];
            for (size_t r = 0; r < rp.size(); r++)
              wc[rp[r]] += wq[rJ[r] * nJ] * zc[re[r]];
          }
        }
      }
    }
    double tb1 = omp_get_wtime();
    arma::mat Zk = perfplan::AtWA(P.A[ch], WZ, P.A[ch]);
    double tb2 = omp_get_wtime();
    perfplan::mask3(cutZ, ch, ch, Zk);   // dE3max / occnat3 cut on the output
    Z3.AddChannelMatrix_pn(ch, Zk, +1);   // = AddToME_pn_ch for iket <= ibra
    tt_wz += tb1 - tb0; tt_a += tb2 - tb1; tt_add += omp_get_wtime() - tb2;
  }
  Z.profiler.timer["_comm333_pph_plan_tt_back_wz"] += tt_wz;
  Z.profiler.timer["_comm333_pph_plan_tt_back_A"] += tt_a;
  Z.profiler.timer["_comm333_pph_plan_tt_back_add"] += tt_add;
  Z.profiler.timer["_comm333_pph_plan_back"] += omp_get_wtime() - t;
  perfplan::release_memory();
  Z.profiler.timer["_perfplan_release_memory"] = perfplan::t_release;   // cumulative, all plan kernels
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
