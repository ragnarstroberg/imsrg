//=====================================================================================================================
//  comm223ss, scalar X and isospin-changing Y. Notation: PerfPlan.inc, PerfPlan333tz.inc
//=====================================================================================================================
// comm223ss for a scalar X and a channel-changing Y (J rank 0, parity 0, T rank > 0), precompiled bookkeeping.
// Same maths as comm223ss (which handles the channel pairs itself), per stored block (cb, ck) of Z, cb != ck:
//   W_Z[(12 J12; 3), (45 J45; 6)] = sum_a hat 6j(j3 J J12; j6 ja J45) (X_{12,6a} Y_{3a,45} - Y_{12,6a} X_{3a,45})
// then Z(cb, ck) = A223(cb)^T W_Z A223(ck). The two terms sum over different orbit sets a (X conserves Tz, Y shifts it).
// Included after PerfPlan223.inc (inside namespace Commutator).
void comm223ss_plan_tz(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  perfplan::Plan &P = perfplan::get(Z, Z);
  const perfplan::Cut3 &cutZ = perfplan::cut3(*Z.modelspace);
  perfplan::Plan223 &Q = perfplan::get223(Z, P);
  auto &Z3 = Z.ThreeBody;
  int norb = P.norb, nJ = P.J2max + 1;
  size_t nch2 = ms.GetNumberTwoBodyChannels();
  int hY = Y.IsHermitian() ? 1 : -1;
  std::vector<const arma::mat *> XM(nch2);
  for (size_t ch = 0; ch < nch2; ch++) XM[ch] = &X.TwoBody.GetMatrix(ch, ch);
  // Y blocks in both orientations: (chb, chk) stored, (chk, chb) = hY M^T (J equal, so no extra phase)
  std::vector<arma::mat> YF(nch2 * nch2);
  for (auto &kv : Y.TwoBody.MatEl)
  {
    size_t cb = kv.first[0], ck = kv.first[1];
    YF[cb * nch2 + ck] = kv.second;
    if (cb != ck) YF[ck * nch2 + cb] = hY * kv.second.t();
  }
  auto pairref = [&](int p, int a, int J, int &g, double &ph)
  {
    g = P.g2idx[((size_t)std::min(p, a) * norb + std::max(p, a)) * (P.J2max + 1) + J];
    if (g < 0) return;
    ph = (p > a) ? Q.kph[g] : 1.0;
    if (p == a) ph *= PhysConst::SQRT2;
  };
  auto Yel = [&](int gb, int gk) -> double
  {
    const arma::mat &M = YF[(size_t)Q.chOf[gb] * nch2 + Q.chOf[gk]];
    return M.n_elem == 0 ? 0.0 : M(Q.locOf[gb], Q.locOf[gk]);
  };
  std::vector<std::array<size_t, 2>> blocks;
  for (auto &it : Z3.Get_ch_start()) blocks.push_back({it.first.ch_bra, it.first.ch_ket});
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ib = 0; ib < blocks.size(); ib++)
  {
    size_t cb = blocks[ib][0], ck = blocks[ib][1];
    if (cb == ck) continue;
    auto &Tb = ms.GetThreeBodyChannel(cb); auto &Tk = ms.GetThreeBodyChannel(ck);
    if (Tb.twoJ != Tk.twoJ) continue;
    int twoJ = Tb.twoJ;
    size_t nPb = P.Pst[cb].size(), nPk = P.Pst[ck].size(), nb = Tb.GetNumberKets(), nk = Tk.GetNumberKets();
    if (nPb == 0 or nPk == 0 or nb == 0 or nk == 0) continue;
    arma::mat WZ(nPb, nPk, arma::fill::zeros);
    std::vector<std::vector<arma::uword>> rowsBy(norb), colsBy(norb);
    for (size_t p = 0; p < nPb; p++) rowsBy[P.Pst[cb][p][1]].push_back(p);
    for (size_t p = 0; p < nPk; p++) colsBy[P.Pst[ck][p][1]].push_back(p);
    std::vector<double> wt;
    for (int I3 = 0; I3 < norb; I3++)
    {
      auto &R = rowsBy[I3];
      if (R.empty()) continue;
      Orbit &o3 = ms.GetOrbit(I3);
      size_t nr = R.size();
      arma::uvec ridx(R);
      std::vector<int> rg(nr), rJ(nr); std::vector<char> hasJ(nJ, 0);
      for (size_t r = 0; r < nr; r++) { rg[r] = P.Pst[cb][R[r]][0]; rJ[r] = P.g2[rg[r]][2]; hasJ[rJ[r]] = 1; }
      Orbit &o1 = ms.GetOrbit(P.g2[rg[0]][0]); Orbit &o2 = ms.GetOrbit(P.g2[rg[0]][1]);
      int par12 = (o1.l + o2.l) % 2, tz12 = o1.tz2 + o2.tz2;
      for (int I6 = 0; I6 < norb; I6++)
      {
        auto &Cs = colsBy[I6];
        if (Cs.empty()) continue;
        Orbit &o6 = ms.GetOrbit(I6);
        size_t nc = Cs.size();
        int h0 = P.Pst[ck][Cs[0]][0];
        Orbit &o4 = ms.GetOrbit(P.g2[h0][0]); Orbit &o5 = ms.GetOrbit(P.g2[h0][1]);
        int par45 = (o4.l + o5.l) % 2, tz45 = o4.tz2 + o5.tz2;
        // term 1: X(12; 6a) with (6a) in the channel of 12; term 2: X(3a; 45) with (3a) in the channel of 45
        std::vector<int> al1, al2;
        for (int a = 0; a < norb; a++)
        {
          Orbit &oa = ms.GetOrbit(a);
          if ((o6.l + oa.l) % 2 == par12 and o6.tz2 + oa.tz2 == tz12) al1.push_back(a);
          if ((o3.l + oa.l) % 2 == par45 and o3.tz2 + oa.tz2 == tz45) al2.push_back(a);
        }
        size_t na1 = al1.size(), na2 = al2.size();
        if (na1 == 0 and na2 == 0) continue;
        arma::mat U1(nr, na1, arma::fill::zeros), U2(nr, na2, arma::fill::zeros);
        for (size_t r = 0; r < nr; r++)
        {
          int g12 = rg[r], J12 = rJ[r], ch12 = Q.chOf[g12], loc12 = Q.locOf[g12];
          for (size_t ia = 0; ia < na1; ia++)
          {
            int g6a; double ph = 1; pairref(I6, al1[ia], J12, g6a, ph);
            if (g6a < 0 or Q.chOf[g6a] != ch12) continue;
            U1(r, ia) = (*XM[ch12])(loc12, Q.locOf[g6a]) * ph;
          }
          for (size_t ia = 0; ia < na2; ia++)
          {
            int g6a; double ph = 1; pairref(I6, al2[ia], J12, g6a, ph);
            if (g6a < 0) continue;
            U2(r, ia) = Yel(g12, g6a) * ph;
          }
        }
        std::map<int, std::vector<size_t>> colsByJ;
        for (size_t c = 0; c < nc; c++) colsByJ[P.g2[P.Pst[ck][Cs[c]][0]][2]].push_back(c);
        for (auto &cj : colsByJ)
        {
          int J45 = cj.first; auto &cl = cj.second; size_t ncj = cl.size();
          arma::mat V1(ncj, na1, arma::fill::zeros), V2(ncj, na2, arma::fill::zeros);
          arma::uvec cidx(ncj);
          for (size_t q = 0; q < ncj; q++)
          {
            int g45 = P.Pst[ck][Cs[cl[q]]][0], ch45 = Q.chOf[g45], loc45 = Q.locOf[g45];
            cidx(q) = Cs[cl[q]];
            for (size_t ia = 0; ia < na1; ia++)
            {
              int g3a; double ph = 1; pairref(I3, al1[ia], J45, g3a, ph);
              if (g3a < 0) continue;
              V1(q, ia) = Yel(g3a, g45) * ph;
            }
            for (size_t ia = 0; ia < na2; ia++)
            {
              int g3a; double ph = 1; pairref(I3, al2[ia], J45, g3a, ph);
              if (g3a < 0 or Q.chOf[g3a] != ch45) continue;
              V2(q, ia) = (*XM[ch45])(Q.locOf[g3a], loc45) * ph;
            }
          }
          arma::mat acc(nr, ncj, arma::fill::zeros);
          for (int term = 0; term < 2; term++)
          {
            auto &al = term == 0 ? al1 : al2;
            size_t na = al.size();
            if (na == 0) continue;
            arma::mat &U = term == 0 ? U1 : U2;
            arma::mat &V = term == 0 ? V1 : V2;
            wt.assign((size_t)nJ * na, 0.0);
            for (int J12 = 0; J12 < nJ; J12++)
              if (hasJ[J12])
              {
                double hat = sqrt((2 * J12 + 1) * (2 * J45 + 1));
                for (size_t ia = 0; ia < na; ia++) wt[J12 * na + ia] = hat * Q.six6(o3.j2, twoJ, 2 * J12, o6.j2, ms.GetOrbit(al[ia]).j2, 2 * J45);
              }
            arma::mat Wr(nr, na);
            for (size_t ia = 0; ia < na; ia++)
              for (size_t r = 0; r < nr; r++) Wr(r, ia) = wt[rJ[r] * na + ia] * U(r, ia);
            if (term == 0) acc += Wr * V.t();
            else acc -= Wr * V.t();
          }
          WZ.submat(ridx, cidx) += acc;
        }
      }
    }
    arma::mat Zk = perfplan::AtWA(Q.A[cb], WZ, Q.A[ck]);
    for (size_t ibra = 0; ibra < nb; ibra++)
      for (size_t iket = 0; iket < nk; iket++)
        if (perfplan::pass3(cutZ, cb, ibra) and perfplan::pass3(cutZ, ck, iket)) Z3.AddToME_pn_ch(cb, ck, ibra, iket, Zk(ibra, iket));
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
