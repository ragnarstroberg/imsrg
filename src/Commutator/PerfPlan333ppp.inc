//=====================================================================================================================
//  comm333_ppp_hhhss with precompiled bookkeeping (scalar X, Y, default cuts). Notation: see PerfPlan.inc.
//=====================================================================================================================
//  Formula (IMSRG3Commutators.cc):
//    Z_ijklmn = 1/6 sum_{abc} (na nb nc + n'a n'b n'c) (X_{ijk,abc} Y_{abc,lmn} - Y_{ijk,abc} X_{abc,lmn})
//  Per 3B channel a GEMM over the intermediate stored kets abc, as in the original:
//    Z3 = hY X(:, abc) F Y(:, abc)^T,     Z3 -= hX hY Z3^T,     F = occupation * symmetry factor of the ket abc.
//  The plan only takes the stored-ket blocks X, Y from the unpack cache (K) and writes with AddChannelMatrix_pn; on
//  its own this term is cheap, the gain is in sharing the unpacked blocks with the other terms.
//  Included after PerfPlan.inc (inside namespace Commutator).
void comm333_ppp_hhhss_plan(const Operator &X, const Operator &Y, Operator &Z)
{
  double tstart = omp_get_wtime();
  ModelSpace &ms = *Z.modelspace;
  if (X.GetTRank() != 0 or Y.GetTRank() != 0 or X.GetParity() != 0 or Y.GetParity() != 0 or X.GetJRank() != 0 or Y.GetJRank() != 0
      )
  {
    comm333_ppp_hhhss(X, Y, Z);
    return;
  }
  if (not(X.ThreeBody.IsAllocated() and Y.ThreeBody.IsAllocated())) return;
  int hX = X.IsHermitian() ? 1 : -1, hY = Y.IsHermitian() ? 1 : -1;
  perfplan::Plan &P = perfplan::get(Z, Z);
  const perfplan::Cut3 &cutZ = perfplan::cut3(ms);
  const std::vector<arma::mat> &XK = perfplan::unpack_k(X, P);
  const std::vector<arma::mat> &YK = perfplan::unpack_k(Y, P);
  auto &Z3 = Z.ThreeBody;
  size_t nch3 = ms.GetNumberThreeBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t ich = 0; ich < nch3; ich++)
  {
    size_t ch = P.ord3[ich];
    auto &T = ms.GetThreeBodyChannel(ch);
    size_t nk = T.GetNumberKets();
    if (nk == 0 or XK[ch].n_elem == 0 or YK[ch].n_elem == 0) continue;
    std::vector<arma::uword> abc; std::vector<double> fac;
    for (size_t k = 0; k < nk; k++)
    {
      Ket3 &ket = T.GetKet(k);
      double na = ket.op->occ, nb = ket.oq->occ, nc = ket.oR->occ;
      double occ_factor = na * nb * nc + (1 - na) * (1 - nb) * (1 - nc);
      if (std::abs(occ_factor) < 1e-6) continue;
      double symm = 1;
      if (ket.p == ket.q and ket.p == ket.r) symm = 1. / 6;
      else if (ket.p == ket.q or ket.q == ket.r) symm = 3. / 6;
      abc.push_back(k); fac.push_back(occ_factor * symm);
    }
    if (abc.empty()) continue;
    arma::uvec uabc(abc);
    arma::mat XM = XK[ch].cols(uabc);
    XM.each_row() %= arma::rowvec(fac);
    arma::mat ZM = hY * XM * YK[ch].cols(uabc).t();
    arma::mat Zt = ZM - hX * hY * ZM.t();
    perfplan::mask3(cutZ, ch, ch, Zt);   // dE3max / occnat3 cut on the output
    Z3.AddChannelMatrix_pn(ch, Zt, -1);   // = AddToME_pn_ch for iket >= ibra
  }
  Z.profiler.timer[__func__] += omp_get_wtime() - tstart;
}
