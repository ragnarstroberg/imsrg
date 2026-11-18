//=====================================================================================================================
//  comm223ss with precompiled bookkeeping (scalar X, Y, default cuts). Notation: see the reader's guide in PerfPlan.inc.
//=====================================================================================================================
//  Formula (IMSRG3Commutators.cc):
//    Z_ijklmn = P(ij/k) P(lm/n) sum_a (X_{ij,la} Y_{ak,mn} - Y_{ij,la} X_{ak,mn})
//  (2B x 2B -> 3B, one contracted orbit a). In the product basis, with rows (12 J12; 3) and columns (45 J45; 6),
//    W_Z[(12 J12; 3), (45 J45; 6)] = sum_a sqrt((2J12+1)(2J45+1)) {j3 J J12; j6 ja J45} (X_{12,6a} Y_{3a,45} - Y_{12,6a} X_{3a,45})
//  with un-normalised (6a), (3a) matrix elements, and then  Z = A223^T W_Z A223.  A223 is this term's own
//  antisymmetriser: it contains exactly the permutation recoupling of the original (RecouplingCoefficient *
//  PermutationPhase * sqrt2 for identical pair orbits * ket phase), which differs in convention from P.A.
//  Steps, per 3B channel and per block of fixed third orbits (I3 for the rows, I6 for the columns):
//    Step 1  row factors    UX(r, a) = X_{12, I6 a},  UY(r, a) = Y_{12, I6 a}    (rows r = P(g12, I3), a over the orbits
//                           whose pair (I6 a) has the parity and Tz of g12)
//    Step 2  column factors VY(q, a) = Y_{I3 a, 45},  VX(q, a) = X_{I3 a, 45},   grouped by J45
//    Step 3  W_Z(r, q) += sum_a w(J12, a; J45) (UX(r, a) VY(q, a) - UY(r, a) VX(q, a)),  w = the hat * 6j factor above,
//            written as column axpys (the blocks, ~40 x 10 x 10, are too small for BLAS)
//    Step 4  Z += A223^T W_Z A223.
//  Included after PerfPlan.inc (inside namespace Commutator).
namespace perfplan
{
  struct Plan223
  {
    const ModelSpace *ms = nullptr; std::vector<double> occ; int E3max = -1;
    std::vector<arma::sp_mat> A;       // ch3 -> n_P x n_k
    std::vector<int> chOf, locOf;      // global 2B ket g -> channel, local index
    std::vector<double> kph;           // g -> Ket::Phase(J) (phase for the reversed order)
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

  static Plan223 &get223(Operator &Z, Plan &P)
  {
    static Plan223 Q; static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    ModelSpace &ms = *Z.modelspace;
    if (Q.ms == &ms and Q.E3max == ms.GetE3max() and Q.occ == P.occ) return Q;
    double t0 = omp_get_wtime();
    Q = Plan223(); Q.ms = &ms; Q.occ = P.occ; Q.E3max = ms.GetE3max();
    int ng = P.g2.size(), norb = P.norb;
    Q.chOf.resize(ng); Q.locOf.resize(ng); Q.kph.resize(ng);
    for (int g = 0; g < ng; g++)
    {
      int p = P.g2[g][0], q = P.g2[g][1], J = P.g2[g][2];
      Orbit &op = ms.GetOrbit(p); Orbit &oq = ms.GetOrbit(q);
      int ch = ms.GetTwoBodyChannelIndex(J, (op.l + oq.l) % 2, (op.tz2 + oq.tz2) / 2);
      auto &tbc = ms.GetTwoBodyChannel(ch);
      size_t loc = tbc.GetLocalIndex(p, q);
      Q.chOf[g] = ch; Q.locOf[g] = loc; Q.kph[g] = tbc.GetKet(loc).Phase(J);
    }
    auto &Z3 = Z.ThreeBody;
    const std::array<ThreeBodyStorage::Permutation, 3> perms = {ThreeBodyStorage::ABC, ThreeBodyStorage::CBA, ThreeBodyStorage::ACB};
    size_t nch3 = ms.GetNumberThreeBodyChannels();
    Q.A.resize(nch3);
    for (size_t ch = 0; ch < nch3; ch++)
    {
      auto &T = ms.GetThreeBodyChannel(ch);
      size_t nP = P.Pst[ch].size(), nk = T.GetNumberKets();
      int twoJ = T.twoJ;
      std::vector<arma::uword> ar, ac; std::vector<double> av;
      for (size_t k = 0; k < nk; k++)
      {
        auto &ket = T.GetKet(k);
        size_t l = ket.p, m = ket.q, n = ket.r;
        double jl = 0.5 * ms.GetOrbit(l).j2, jm = 0.5 * ms.GetOrbit(m).j2, jn = 0.5 * ms.GetOrbit(n).j2;
        int J2 = ket.Jpq;
        for (int ip = 0; ip < 3; ip++)
        {
          size_t I4, I5, I6;
          Z3.Permute(perms[ip], l, m, n, I4, I5, I6);
          Orbit &o4 = ms.GetOrbit(I4); Orbit &o5 = ms.GetOrbit(I5); Orbit &o6 = ms.GetOrbit(I6);
          int jmin = J2, jmax = J2;
          if (perms[ip] != ThreeBodyStorage::ABC)
          {
            jmin = std::max(std::abs(o4.j2 - o5.j2), std::abs(twoJ - o6.j2)) / 2;
            jmax = std::min(o4.j2 + o5.j2, twoJ + o6.j2) / 2;
          }
          for (int J2p = jmin; J2p <= jmax; J2p++)
          {
            int g = P.g2idx[((size_t)std::min(I4, I5) * norb + std::max(I4, I5)) * (P.J2max + 1) + J2p];
            if (g < 0) continue;
            int pidx = P.pidx(ch, g, I6);
            if (pidx < 0) continue;
            double rec = Z3.RecouplingCoefficient(perms[ip], jl, jm, jn, J2p, J2, twoJ) * Z3.PermutationPhase(perms[ip]);
            if (I4 == I5) rec *= PhysConst::SQRT2;
            if (I4 > I5) rec *= Q.kph[g];
            ar.push_back(pidx); ac.push_back(k); av.push_back(rec);
          }
        }
      }
      arma::umat loc(2, ar.size());
      for (size_t t = 0; t < ar.size(); t++) { loc(0, t) = ar[t]; loc(1, t) = ac[t]; }
      Q.A[ch] = arma::sp_mat(true, loc, arma::vec(av), nP, nk);
    }
    // 6j {j3 J J12; j6 ja J45}, twice-integer arguments (j2, tJ, 2J, j2, j2, 2J)
    // j3, J, J12, j6, J45 belong to 3B states; ja (summed) runs over all orbits
    Q.d6[0] = P.j2max3 + 1; Q.d6[1] = P.tJmax3 + 2; Q.d6[2] = 2 * P.J2max3 + 1; Q.d6[3] = P.j2max3 + 1; Q.d6[4] = P.j2max + 1; Q.d6[5] = 2 * P.J2max3 + 1;
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
    Z.profiler.timer["_comm223_plan_build"] += omp_get_wtime() - t0;
    return Q;
  }
} // namespace perfplan

void comm223ss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0)
  {
    comm223ss(X, Y, Z);
    return;
  }
  perfplan::Plan &P = perfplan::get(Z, Z);
  perfplan::Plan223 &Q = perfplan::get223(Z, P);
  const perfplan::Cut3 &cutZ = perfplan::cut3(*Z.modelspace);
  ModelSpace &ms = *Z.modelspace;
  auto &Z3 = Z.ThreeBody;
  size_t nch3 = ms.GetNumberThreeBodyChannels();
  size_t nch2 = ms.GetNumberTwoBodyChannels();
  int norb = P.norb;
  std::vector<const arma::mat *> XM(nch2), YM(nch2);
  for (size_t ch = 0; ch < nch2; ch++) { XM[ch] = &X.TwoBody.GetMatrix(ch, ch); YM[ch] = &Y.TwoBody.GetMatrix(ch, ch); }
  // pair (p, a) with J -> (g, phase incl. sqrt2 for p == a), as the original's un-normalised MEs
  auto pairref = [&](int p, int a, int J, int &g, double &ph)
  {
    g = P.g2idx[((size_t)std::min(p, a) * norb + std::max(p, a)) * (P.J2max + 1) + J];
    if (g < 0) return;
    ph = (p > a) ? Q.kph[g] : 1.0;
    if (p == a) ph *= PhysConst::SQRT2;
  };
  double t = omp_get_wtime();
  const bool flags = perturbative_triples or imsrg3_no_qqq or imsrg3_only_vvv;
  double emp2_tot = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : emp2_tot)
  for (size_t ich = 0; ich < nch3; ich++)
  {
    size_t ch = P.ord3[ich];
    auto &T = ms.GetThreeBodyChannel(ch);
    size_t nP = P.Pst[ch].size(), nk = T.GetNumberKets();
    if (nP == 0 or nk == 0) continue;
    int twoJ = T.twoJ;
    arma::mat WZ(nP, nP, arma::fill::zeros);
    // Blocked by the third orbit: rows R = P states (g12, I3), cols C = P states (g45, I6). For fixed (I3, I6),
    //   WZ(R, C) = sum_a hat * 6j(j3 J J12; j6 ja J45) [X(g12; I6 a) Y(I3 a; g45) - Y(g12; I6 a) X(I3 a; g45)]
    // computed as small GEMMs over a, one per distinct J45 among the columns.
    std::vector<std::vector<arma::uword>> byOrb(norb);
    for (size_t p = 0; p < nP; p++) byOrb[P.Pst[ch][p][1]].push_back(p);
    int nJ = P.J2max + 1;
    std::vector<double> wt;
    for (int I3 = 0; I3 < norb; I3++)
    {
      auto &R = byOrb[I3];
      if (R.empty()) continue;
      Orbit &o3 = ms.GetOrbit(I3);
      size_t nr = R.size();
      size_t r0 = R[0];
      if (R.back() - r0 + 1 != nr) throw std::runtime_error("comm223ss_plan: P states with one third orbit are not contiguous");
      std::vector<int> rJ(nr); std::vector<char> hasJ(nJ, 0);
      for (size_t r = 0; r < nr; r++) { rJ[r] = P.g2[P.Pst[ch][R[r]][0]][2]; hasJ[rJ[r]] = 1; }
      // all pairs g12 of these rows share parity and Tz (fixed by the 3B channel and I3)
      int g0 = P.Pst[ch][R[0]][0];
      Orbit &o1 = ms.GetOrbit(P.g2[g0][0]); Orbit &o2 = ms.GetOrbit(P.g2[g0][1]);
      int par12 = (o1.l + o2.l) % 2, tz12 = o1.tz2 + o2.tz2;
      for (int I6 = 0; I6 < norb; I6++)
      {
        auto &Cs = byOrb[I6];
        if (Cs.empty()) continue;
        Orbit &o6 = ms.GetOrbit(I6);
        size_t nc = Cs.size();
        // only orbits a with (I6, a) in the parity / Tz of g12 contribute
        std::vector<int> al;
        for (int a = 0; a < norb; a++)
        {
          Orbit &oa = ms.GetOrbit(a);
          if ((o6.l + oa.l) % 2 == par12 and o6.tz2 + oa.tz2 == tz12) al.push_back(a);
        }
        size_t na = al.size();
        if (na == 0) continue;
        // Step 1: row factors over a: X(g12; I6 a), Y(g12; I6 a) (un-normalised, with ordering phase)
        arma::mat UX(nr, na, arma::fill::zeros), UY(nr, na, arma::fill::zeros);
        for (size_t r = 0; r < nr; r++)
        {
          int g12 = P.Pst[ch][R[r]][0], J12 = rJ[r], ch12 = Q.chOf[g12], loc12 = Q.locOf[g12];
          for (size_t ia = 0; ia < na; ia++)
          {
            int g6a; double ph6a = 1;
            pairref(I6, al[ia], J12, g6a, ph6a);
            if (g6a < 0 or Q.chOf[g6a] != ch12) continue;
            UX(r, ia) = (*XM[ch12])(loc12, Q.locOf[g6a]) * ph6a; UY(r, ia) = (*YM[ch12])(loc12, Q.locOf[g6a]) * ph6a;
          }
        }
        // Step 2: column factors over a: Y(I3 a; g45), X(I3 a; g45); grouped by J45
        std::map<int, std::vector<size_t>> colsByJ;
        for (size_t c = 0; c < nc; c++) colsByJ[P.g2[P.Pst[ch][Cs[c]][0]][2]].push_back(c);
        for (auto &cj : colsByJ)
        {
          int J45 = cj.first; auto &cl = cj.second; size_t ncj = cl.size();
          arma::mat VY(ncj, na, arma::fill::zeros), VX(ncj, na, arma::fill::zeros);
          arma::uvec cidx(ncj);
          for (size_t q = 0; q < ncj; q++)
          {
            int g45 = P.Pst[ch][Cs[cl[q]]][0], ch45 = Q.chOf[g45], loc45 = Q.locOf[g45];
            cidx(q) = Cs[cl[q]];
            for (size_t ia = 0; ia < na; ia++)
            {
              int g3a; double ph3a = 1;
              pairref(I3, al[ia], J45, g3a, ph3a);
              if (g3a < 0 or Q.chOf[g3a] != ch45) continue;
              VY(q, ia) = (*YM[ch45])(Q.locOf[g3a], loc45) * ph3a; VX(q, ia) = (*XM[ch45])(Q.locOf[g3a], loc45) * ph3a;
            }
          }
          // weights hat * 6j(j3 J J12; j6 ja J45) depend only on (J12, a): small table
          wt.assign((size_t)nJ * na, 0.0);
          for (int J12 = 0; J12 < nJ; J12++)
            if (hasJ[J12])
            {
              double hat = sqrt((2 * J12 + 1) * (2 * J45 + 1));
              for (size_t ia = 0; ia < na; ia++) wt[J12 * na + ia] = hat * Q.six6(o3.j2, twoJ, 2 * J12, o6.j2, ms.GetOrbit(al[ia]).j2, 2 * J45);
            }
          arma::mat WXr(nr, na), WYr(nr, na);
          for (size_t ia = 0; ia < na; ia++)
            for (size_t r = 0; r < nr; r++)
            {
              double w = wt[rJ[r] * na + ia];
              WXr(r, ia) = w * UX(r, ia); WYr(r, ia) = w * UY(r, ia);
            }
          // Step 3: WZ(R, C) += WXr VY^T - WYr VX^T. The rows R are the contiguous P range [r0, r0 + nr) (third-orbit-major P
          // basis), and the GEMMs are too small for BLAS to pay off: contiguous column axpys instead.
          for (size_t q = 0; q < ncj; q++)
          {
            double *wc = WZ.colptr(cidx(q)) + r0;
            for (size_t ia = 0; ia < na; ia++)
            {
              double vy = VY(q, ia), vx = VX(q, ia);
              if (vy == 0 and vx == 0) continue;
              const double *ux = WXr.colptr(ia), *uy = WYr.colptr(ia);
              for (size_t r = 0; r < nr; r++) wc[r] += ux[r] * vy - uy[r] * vx;
            }
          }
        }
      }
    }
    arma::mat Zk = perfplan::AtWA(Q.A[ch], WZ, Q.A[ch]);   // Step 4: Z = A223^T W_Z A223
    perfplan::mask3(cutZ, ch, ch, Zk);   // dE3max / occnat3 cut on the output
    if (flags)
    {
      // imsrg3_no_qqq / imsrg3_only_vvv / perturbative_triples filters on the output states, as in comm223ss
      auto &Tc = ms.GetThreeBodyChannel(ch);
      std::vector<double> occ(nk), unocc(nk); std::vector<char> ok(nk, 1), vvv(nk, 0); std::vector<int> symm(nk, 6);
      for (size_t k = 0; k < nk; k++)
      {
        Ket3 &kt = Tc.GetKet(k);
        occ[k] = kt.op->occ * kt.oq->occ * kt.oR->occ; unocc[k] = (1 - kt.op->occ) * (1 - kt.oq->occ) * (1 - kt.oR->occ);
        int cvq = kt.op->cvq + kt.oq->cvq + kt.oR->cvq;
        vvv[k] = (kt.op->cvq == 1 and kt.oq->cvq == 1 and kt.oR->cvq == 1);
        if (perturbative_triples and occ[k] < 1e-8 and unocc[k] < 1e-8) ok[k] = 0;
        if (imsrg3_no_qqq and cvq > 5) ok[k] = 0;
        if (imsrg3_only_vvv and not vvv[k]) ok[k] = 0;
        if (kt.p == kt.q and kt.p == kt.r) symm[k] = 1; else if (kt.p == kt.q or kt.p == kt.r or kt.q == kt.r) symm[k] = 3;
      }
      double emp2 = 0;
      for (size_t ib = 0; ib < nk; ib++)
        for (size_t ik = 0; ik <= ib; ik++)   // iket <= ibra, as the original loops
        {
          bool use = ok[ib] and ok[ik];
          if (perturbative_triples and std::abs(occ[ib] * unocc[ik] - unocc[ib] * occ[ik]) < 1e-8) use = false;
          if (perturbative_triples and pert_trip_novvv and vvv[ib] and vvv[ik]) use = false;
          if (not use) { Zk(ib, ik) = 0; continue; }
          if (perturbative_triples)
          {
            Ket3 &b = Tc.GetKet(ib); Ket3 &k = Tc.GetKet(ik);
            double Eb = Z.OneBody(b.p, b.p) + Z.OneBody(b.q, b.q) + Z.OneBody(b.r, b.r);
            double Ek = Z.OneBody(k.p, k.p) + Z.OneBody(k.q, k.q) + Z.OneBody(k.r, k.r);
            double z = Zk(ib, ik);
            emp2 += 1. / 36 * symm[ib] * symm[ik] * (twoJ + 1) * z * z * (occ[ib] * unocc[ik] - occ[ik] * unocc[ib]) / (Eb - Ek);
          }
        }
      if (perturbative_triples) { emp2_tot += emp2; continue; }   // energy only, Z3 is not stored
    }
    Z3.AddChannelMatrix_pn(ch, Zk, +1);   // = AddToME_pn_ch for iket <= ibra
  }
  if (perturbative_triples) Z.ZeroBody = emp2_tot;   // as comm223ss: the energy goes into the zero-body part
  Z.profiler.timer["_comm223_plan_main"] += omp_get_wtime() - t;
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
