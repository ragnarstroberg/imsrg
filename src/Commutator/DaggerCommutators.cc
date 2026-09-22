
#include "DaggerCommutators.hh"
#include "Commutator.hh"
#include "PhysicalConstants.hh"


namespace Commutator
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Begin Scalar-Dagger Commutator Terms////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///  In the below formulas, the number in parentheses indicates the number of external legs of a diagram, rather
  ///  than the usual particle rank. In this notation, a one-body operator has 2 external legs and so would be indicated
  ///  with a (2). The (1) then means a single creation/annihilation operator, and (3) means two creators and one annihilator.

  //*****************************************************************************************
  // [X^(2), Y^(1)]^(1)
  //
  //      |i                 i|
  //      |                   |
  //     (X)         ---      |  (Y)
  //       a\                 |  /a
  //         (Y)             (X)/
  //
  //
  // Sum_ia X_ia Y_a^{lambda}
  // Adapted from
  //    void Operator::comm111ss( const Operator & X, const Operator& Y)
  // This could be modified to be more efficent since we only need one column from Y to go to one column in Z.
  // But this is so far from being the bottleneck that it makes more sense to be clear.
  void comm211sd(const Operator &X, const Operator &Y, Operator &Z)
  {
    Z.OneBody += X.OneBody * Y.OneBody;
  }

  //*****************************************************************************************
  // [X^(2),Y^(3)]^(1)
  //         |i             |i
  //         |              |
  //  (X)_ b |         a _(Y)
  //    \_\  |   ---    /_/
  //    a  (Y)         (X) b
  //
  // Sum_ab Sum_J (na n`b - n`a nb) X_ab Y_bia  hat(J)^2/hat(lambda)^2
  // or
  // Sum_ab Sum_J (2J+1)/(2lambda+1)  (na n`b) ( X_ab Y_bia - X_ba Y_aib )
  //
  // Adapted from
  //    void Operator::comm121ss( const Operator& X, const Operator& Y)
  void comm231sd(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    if (Y.legs < 3)
      return;
    //   index_t norbits = Z.modelspace->GetNumberOrbits();
    index_t Q = Z.GetQSpaceOrbit();
    Orbit &oQ = Z.modelspace->GetOrbit(Q);
    //   for (index_t i=0;i<norbits;++i)
    //   #pragma omp parallel for
    //   for (auto i : Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
    for (auto i : Z.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2))
    {
      for (auto &a : Z.modelspace->holes) // C++11 syntax
      {
        Orbit &oa = Z.modelspace->GetOrbit(a);
        //             for (index_t b=0; b<norbits; ++b)
        //             for (auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2} ))
        for (auto b : X.GetOneBodyChannel(oa.l, oa.j2, oa.tz2))
        {
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double nanb = oa.occ * (1 - ob.occ); // despite what the name suggests, nanb is n_a * (1-n_b).
          if (std::abs(nanb) < ModelSpace::OCC_CUT)
            continue;
          int Jmin = std::abs(oa.j2 - oQ.j2) / 2;
          int Jmax = (oa.j2 + oQ.j2) / 2;
          double Ymon_bia = 0;
          double Ymon_aib = 0;
          for (int J = Jmin; J <= Jmax; J++)
          {
            Ymon_bia += (2 * J + 1.0) / (oQ.j2 + 1.0) * Y.ThreeLeg.GetME_J(J, b, i, a);
            Ymon_aib += (2 * J + 1.0) / (oQ.j2 + 1.0) * Y.ThreeLeg.GetME_J(J, a, i, b);
          }
          //                Z.OneBody(i,Q) += nanb * X.OneBody(a,b) * Ymon_bia - X.OneBody(b,a) * Ymon_aib;
          Z.OneBody(i, 0) += nanb * X.OneBody(a, b) * Ymon_bia - X.OneBody(b, a) * Ymon_aib;
          //                  Z.OneBody(i,Q) += (ob.j2+1) * nanb *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,Q) ;  // Is this still the right way to do this? Do we need to worry about normalization? (It looks ok).
          //                  Z.OneBody(i,Q) -= (oa.j2+1) * nanb *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,Q) ;  // GetTBMEmonopole returns unnormalized TBME summed over J times (2J+1)/((2ji+1)*(2jj+1))
        }
      }
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //*****************************************************************************************
  // [X^(2),Y^(3)]^(3)  and [X^(4),Y^(1)]^(3), combined
  //
  //     |   |  /       |  | /             \ /         \ |  (Y)
  //    (X)  | /   __   | (Y)              (X)    __    \| /
  //      \  |/         | /        and     / \          (X)
  //       (Y)         (X)                /  (Y)         |
  //
  //  ZJ_ijkQ = sum_a  X_ia YJ_ajkQ  + X_ja YJ_iakQ + X_ak YJ_ijaQ + XJ_ijka Y_aQ
  //
  //  I adapted this, but then gave up and brute forced it because the intermediate matrix
  //  stuff was too obfuscated, and I'm the one who wrote it... -SRS 14/04/2018
  // Adapted from
  //  void Operator::comm122ss( const Operator& X, const Operator& Y )
  void comm413_233sd(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;

    index_t Qorbit = Z.GetQSpaceOrbit();
    Orbit &oQ = Z.modelspace->GetOrbit(Qorbit);
    //   int norb = Z.modelspace->GetNumberOrbits();

    int n_nonzero = Z.modelspace->SortedTwoBodyChannels.size();
    //   #pragma omp parallel for schedule(dynamic,1)
    for (int ich = 0; ich < n_nonzero; ++ich)
    {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);

      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0; indx_ij < npq; ++indx_ij)
      {
        Ket &bra = tbc.GetKet(indx_ij);
        int i = bra.p;
        int j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double norm_ij = (i == j) ? PhysConst::INVSQRT2 : 1.0;
        //         for ( int k=0; k<norb; ++k )
        for (auto k : Z.modelspace->all_orbits)
        {
          Orbit &ok = Z.modelspace->GetOrbit(k);
          //           if (not tbc.CheckChannel_ket(&ok,&oQ)) continue;
          //           if ( ( (ok.tz2+oQ.tz2)!=2*tbc.Tz ) or ( (ok.l+oQ.l)%2!=tbc.parity ) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue;
          if (((ok.tz2 - oQ.tz2) != 2 * tbc.Tz) or ((ok.l + oQ.l) % 2 != tbc.parity) or ((ok.j2 + oQ.j2) < 2 * tbc.J) or (std::abs(ok.j2 - oQ.j2) > 2 * tbc.J))
            continue;

          //           double norm_kQ = (k==Qorbit) ? PhysConst::INVSQRT2 : 1.0;
          double norm_kQ = 1.0;
          double cijk = 0;

          //           std::cout << "INNER LOOPs ijk = " << i << " " << j << " " << k << " Z has " << Z.GetNumberLegs() << " legs " << " and there are " << Z.ThreeLeg.MatEl.size() << " 3 leg channels " << std::endl;
          // The first three loops are the [2,3]->3 bit
          //           for ( int a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
          for (int a : X.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
          {
            //             cijk += X1(i,a) * Y.TwoBody.GetTBME_norm(ch,ch,a,j,k,Qorbit);
            //             cijk += X1(i,a) * Y.TwoBody.GetTBME(ch,ch,a,j,k,Qorbit);
            cijk += X1(i, a) * Y.ThreeLeg.GetME(ch, a, j, k);
          }
          //           std::cout << "check 1. size of X.OneBodyChannels is " << X.OneBodyChannels.Size() << std::endl;
          //           for ( int a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
          for (int a : X.GetOneBodyChannel(oj.l, oj.j2, oj.tz2))
          {
            //             cijk += X1(j,a) * Y.TwoBody.GetTBME_norm(ch,ch,i,a,k,Qorbit);
            //             cijk += X1(j,a) * Y.TwoBody.GetTBME(ch,ch,i,a,k,Qorbit);
            cijk += X1(j, a) * Y.ThreeLeg.GetME(ch, i, a, k);
          }
          //           std::cout << "check 2" << std::endl;
          //           for ( int a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
          for (int a : X.GetOneBodyChannel(ok.l, ok.j2, ok.tz2))
          {
            //             cijk += X1(a,k) * Y.TwoBody.GetTBME_norm(ch,ch,i,j,a,Qorbit); // should this have a minus sign?
            //             cijk -= X1(a,k) * Y.TwoBody.GetTBME_norm(ch,ch,i,j,a,Qorbit);
            //             cijk -= X1(a,k) * Y.TwoBody.GetTBME(ch,ch,i,j,a,Qorbit);
            cijk -= Y.ThreeLeg.GetME(ch, i, j, a) * X1(a, k);
          }

          //           std::cout << " DOING 413 bit "  << std::endl;
          // This here loop is the [4,1]->3 bit
          //           for ( int a : Y.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
          for (int a : Y.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2))
          {
            //             cijk += Y1(a,Qorbit) * X.TwoBody.GetTBME_norm(ch,ch,i,j,k,a);   // This determines the normalization for the (adagger adagger a) term.
            //             cijk +=  X.TwoBody.GetTBME(ch,ch,i,j,k,a) * Y1(a,Qorbit) ;   // This determines the normalization for the (adagger adagger a) term.
            cijk += X.TwoBody.GetTBME(ch, ch, i, j, k, a) * Y1(a, 0); // This determines the normalization for the (adagger adagger a) term.
                                                                      //             if (i==0 and j==8 and k==0)
                                                                      //             {
                                                                      //                std::cout << "   " << __func__ << "  a = " << a << "  Xijka, Ya = " << X.TwoBody.GetTBME(ch,ch,i,j,k,a) << " , " << Y1(a,Qorbit) << "  cijk = " << cijk << std::endl;
                                                                      //             }
          }

          cijk *= (norm_ij * norm_kQ); // We normalize like a scalar TBME so that the setter/getters make sense. We'll fix the wrong |kQ> normalization when writing to file.

          Z.ThreeLeg.AddToME(ch, i, j, k, cijk); // AddToME  assumes a normalized matrix element.
        }
      }
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void ConstructDaggerMpp_Mhh(const Operator &X, const Operator &Y, const Operator &Z, ThreeLegME &Mpp, ThreeLegME &Mhh)
  {
    int nch = Z.modelspace->SortedTwoBodyChannels.size();
#ifndef OPENBLAS_NOUSEOMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int ich = 0; ich < nch; ++ich)
    {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);

      auto &LHS = X.TwoBody.GetMatrix(ch, ch);
      auto &RHS = Y.ThreeLeg.GetMatrix(ch);

      auto &Matrixpp = Mpp.GetMatrix(ch);
      auto &Matrixhh = Mhh.GetMatrix(ch);

      auto &kets_pp = tbc.GetKetIndex_pp();
      auto &kets_hh = tbc.GetKetIndex_hh();
      auto &kets_ph = tbc.GetKetIndex_ph();
      auto &nanb = tbc.Ket_occ_hh;
      auto &nbarnbar_hh = tbc.Ket_unocc_hh;
      auto &nbarnbar_ph = tbc.Ket_unocc_ph;

      Matrixpp = LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh = LHS.cols(kets_hh) * arma::diagmat(nanb) * RHS.rows(kets_hh);
      if (kets_hh.size() > 0)
        Matrixpp += LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) * RHS.rows(kets_hh);
      if (kets_ph.size() > 0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) * RHS.rows(kets_ph);

    } // for ch
  }

  /// Since comm443_pp_hhsd() and comm431sd() both require the construction of
  /// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
  /// only calculate the intermediates once.
  /// We do the same thing in the standard scalar-scalar commutators.
  void comm433_pp_hh_431sd(const Operator &X, const Operator &Y, Operator &Z)
  {

    double t = omp_get_wtime();

    //   static TwoBodyME Mpp = Z.TwoBody;
    //   static TwoBodyME Mhh = Z.TwoBody;
    ThreeLegME Mpp = Z.ThreeLeg;
    ThreeLegME Mhh = Z.ThreeLeg;
    Mpp *= 0;
    Mhh *= 0;

    ConstructDaggerMpp_Mhh(X, Y, Z, Mpp, Mhh);

    Z.ThreeLeg += Mpp;
    Z.ThreeLeg -= Mhh;

    Z.profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;

    // Now, the one body part
    t = omp_get_wtime();
    int Q = Z.GetQSpaceOrbit();
    Orbit &oQ = Z.modelspace->GetOrbit(Q);
    //   auto ilist = Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2});
    //   std::vector<index_t> ilist(Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}).begin(), Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}).end());
    std::vector<index_t> ilist(Z.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2).begin(), Z.GetOneBodyChannel(oQ.l, oQ.j2, oQ.tz2).end());
    int ni = ilist.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (int i_ind = 0; i_ind < ni; ++i_ind)
    {
      int i = ilist[i_ind];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double cijJ = 0;
      for (int ch = 0; ch < Z.nChannels; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        double Jfactor = (2 * tbc.J + 1.0);
        // Sum c over holes and include the nbar_a * nbar_b terms
        for (auto &c : Z.modelspace->holes)
        {
          Orbit &oc = Z.modelspace->GetOrbit(c);
          //            cijJ += Jfactor * oc.occ * Mpp.GetTBME(ch,c,i,c,Q);     // We use the GetTBME, which returns an unnormalized matrix element, as required.
          //            cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,Q);
          cijJ += Jfactor * oc.occ * Mpp.GetME(ch, c, i, c); // We use the GetTBME, which returns an unnormalized matrix element, as required.
          cijJ += Jfactor * (1 - oc.occ) * Mhh.GetME(ch, c, i, c);
        }
        // Sum c over particles and include the n_a * n_b terms
        for (auto &c : Z.modelspace->particles)
        {
          //            cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,Q);
          cijJ += Jfactor * Mhh.GetME(ch, c, i, c);
        }
      }
      Z.OneBody(i, 0) += cijJ / (oi.j2 + 1.0); // The factor of 1/2 in the formula is absorbed by the fact that the mat-mult only sums a<=b.
                                               //      Z.OneBody(i,Q) += cijJ /(oi.j2+1.0);   // The factor of 1/2 in the formula is absorbed by the fact that the mat-mult only sums a<=b.
    }                                          // for i
    Z.profiler.timer["pphh 413sd"] += omp_get_wtime() - t;
  }

  //*****************************************************************************************
  // [X^(4),Y^(3)]^(3)]  ph piece, the slow way
  //
  //   |           |       |           |
  //   |      _(Y)_|       |_(X)_      |
  //   |     /\            |    /\     |
  //   |    (  )      __   |   (  )    |
  //   |_(X)_\/            |    \/_(Y)_|
  //   |                   |
  //
  //  Straightfoward and very slow implementation, only used for unit testing the
  //  faster mat-mult implementation.  The two implementations agree as of Nov 2, 2018 - SRS
  //  Benchmarked against the m-scheme version in UnitTest.cc, and it agrees (although
  //  some changes were required compared to the impoementation tested in Nov 2018).
  //  So this again agrees with the more efficient way, and with the m-scheme verion as of Sep 30 2019 - SRS
  void comm433sd_ph_dumbway(const Operator &X, const Operator &Y, Operator &Z)
  {

    double t_start = omp_get_wtime();
    //   int norb = Z.modelspace->GetNumberOrbits();
    int nch = Z.modelspace->SortedTwoBodyChannels.size();
    index_t Q = Y.GetQSpaceOrbit();
    Orbit &oQ = Z.modelspace->GetOrbit(Q);
    double jQ = 0.5 * oQ.j2;
    auto kets_ph = Z.modelspace->KetIndex_ph;
    for (auto khh : Z.modelspace->KetIndex_hh) // ph kets needs to include hh list in case there are fractionally occupied states
    {
      kets_ph.push_back(khh);
    }

    for (int ich = 0; ich < nch; ++ich)
    {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      size_t nkets = tbc.GetNumberKets();
      for (size_t ibra = 0; ibra < nkets; ibra++)
      {
        auto &bra = tbc.GetKet(ibra);
        std::vector<index_t> i_cases = {bra.p, bra.q};
        std::vector<index_t> j_cases = {bra.q, bra.p};
        std::vector<int> ijsign_cases = {+1, bra.Phase(J)}; // This accounts for the 1 - (-1)^(ji +jj-J) Pij  factor.
        double norm_ij = (bra.p == bra.q) ? PhysConst::INVSQRT2 : 1.0;

        for (auto k : Z.modelspace->all_orbits)
        {
          Orbit &ok = Z.modelspace->GetOrbit(k);
          //            if ( not tbc.CheckChannel_ket(&ok,&oQ) ) continue;
          //            if ( ( (ok.tz2+oQ.tz2)!=2*tbc.Tz) or ( (ok.l+oQ.l)%2!=tbc.parity) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue; // if |kQ> doesn't live in this channel, move along.
          if (((ok.tz2 - oQ.tz2) != 2 * tbc.Tz) or ((ok.l + oQ.l) % 2 != tbc.parity) or ((ok.j2 + oQ.j2) < 2 * tbc.J) or (std::abs(ok.j2 - oQ.j2) > 2 * tbc.J))
            continue; // if |kQ> doesn't live in this channel, move along.
          double zijk = 0.;
          //            double norm_kQ = (k==Q) ? PhysConst::INVSQRT2 : 1.0;
          double norm_kQ = 1.0;
          double jk = 0.5 * ok.j2;

          for (int ijcase = 0; ijcase <= 1; ijcase++)
          {
            //         if (ijcase==1) continue;
            index_t i = i_cases[ijcase];
            index_t j = j_cases[ijcase];
            //          std::cout << "~~~~ i,j = " << i << " " << j << std::endl;
            int ijsign = ijsign_cases[ijcase];
            Orbit &oi = Z.modelspace->GetOrbit(i);
            Orbit &oj = Z.modelspace->GetOrbit(j);
            double ji = 0.5 * oi.j2;
            double jj = 0.5 * oj.j2;
            //          for ( size_t k=0; k<norb; k++)

            int Jprime_min = std::max(std::abs(ji - jQ), std::abs(jj - jk));
            int Jprime_max = std::min(ji + jQ, jj + jk);
            for (int Jprime = Jprime_min; Jprime <= Jprime_max; ++Jprime)
            {
              double sixjprime = Z.modelspace->GetSixJ(ji, jj, J, jk, jQ, Jprime);
              if (std::abs(sixjprime) < 1e-8)
                continue;

              double XYprod = 0;

              for (auto &iketab : kets_ph)
              {
                auto &ketab = Z.modelspace->GetKet(iketab);
                if ((ketab.op->j2 + ketab.oq->j2) < 2 * Jprime or std::abs(ketab.op->j2 - ketab.oq->j2) > 2 * Jprime)
                  continue;
                std::vector<index_t> ab_cases = {ketab.p, ketab.q};
                for (int abcase = 0; abcase <= 1; abcase++)
                {
                  index_t a = ab_cases[abcase];
                  index_t b = ab_cases[1 - abcase];

                  Orbit &oa = Z.modelspace->GetOrbit(a);
                  Orbit &ob = Z.modelspace->GetOrbit(b);
                  double ja = 0.5 * oa.j2;
                  double jb = 0.5 * ob.j2;
                  double nanb = oa.occ - ob.occ;

                  //                  if ( not (a==0 and b==10) or (a==10 and b==0) ) continue;
                  //                  if ( not ((a==1 and b==8) or (a==8 and b==1)) ) continue;

                  // we can probably also check triangle for (a,b,Jprime)

                  int JA_min = std::max(std::abs(ja - jj), std::abs(jk - jb));
                  int JA_max = std::min(ja + jj, jk + jb);

                  int JB_min = std::max(std::abs(ja - jQ), std::abs(ji - jb));
                  int JB_max = std::min(ja + jQ, ji + jb);
                  double matelX = 0;
                  double matelY = 0;
                  for (int JA = JA_min; JA <= JA_max; ++JA)
                  {
                    double sixjA = Z.modelspace->GetSixJ(ja, jb, Jprime, jk, jj, JA);
                    matelX -= (2 * JA + 1) * sixjA * X.TwoBody.GetTBME_J(JA, a, j, k, b); // GetTBME_J returns an un-normalized matrix element
                  }
                  for (int JB = JB_min; JB <= JB_max; ++JB)
                  {
                    double sixjB = Z.modelspace->GetSixJ(ja, jb, Jprime, ji, jQ, JB);
                    //                    matelY -= (2*JB+1) * sixjB * Y.TwoBody.GetTBME_J(JB,i,b,a,Q);   // GetTBME_J returns an un-normalized matrix element
                    matelY -= (2 * JB + 1) * sixjB * Y.ThreeLeg.GetME_J(JB, i, b, a); // GetTBME_J returns an un-normalized matrix element
                  }
                  XYprod += nanb * matelX * matelY;
                  //                  if ( ( (i==0 and j==1) or (i==1 and j==0)) and k==1  and std::abs(nanb)>0)
                  //                  if ( ( (i==0 and j==1) or (i==1 and j==0)) and k==1 and ((a==0 and b==10) or (b==0 and a==10)) )
                  //                  if ( ( (i==0 and j==1) or (i==1 and j==0)) and k==1 and ((a==1 and b==8) or (b==1 and a==8)) )
                  //                  {
                  //                    std::cout << "   J = " << J << "  a,b = " << a << " " << b << "  Jprime = " << Jprime << "  nanb " << nanb
                  //                              << "  Xbar_abkj Ybar_iQab " << matelX << " " << matelY << " => XYprod = " << XYprod << std::endl;
                  //                  }
                }
              } // loop over ab cases
              zijk -= ijsign * (2 * Jprime + 1) * sixjprime * XYprod;
              //                  if (( (i==0 and j==1) or (i==1 and j==0)) and k==1)
              //                  {
              //                    std::cout << "i,j = " << i << " " << j << " ijsign(J=" << J << ") = " << ijsign << "  2Jprime+1 " << 2*Jprime+1 << "  sixjprime " << sixjprime << "  XYprod " << XYprod << "  => zijk " << zijk
              //                               <<  std::endl;
              //                  }
            } // loop over ab kets

            //                  if (( (i==0 and j==1) or (i==1 and j==0)) and k==1)
            //                  {
            //                    std::cout << " zijk = " << zijk <<  " norm_ij norm_kQ = " << norm_ij << " " << norm_kQ << "    before setting, the ME is " << Z.ThreeLeg.GetME_J(J,i,j,k) << std::endl << std::endl;
            //                  }

          } // loop over ij cases

          //            Z.TwoBody.AddToTBME_J(J,i,j,k,Q,  zijk * norm_ij*norm_kQ);   // AddToTBME_J assumes a normalized matrix element.
          Z.ThreeLeg.AddToME_J(J, bra.p, bra.q, k, zijk * norm_ij * norm_kQ); // AddToTBME_J assumes a normalized matrix element.
        }                                                                     // loop over k
      }                                                                       // for ibra
    }                                                                         // for ich
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //*****************************************************************************************
  // [X^(4),Y^(3)]^(3)]  ph piece
  //
  //   |           |       |           |
  //   |      _(Y)_|       |_(X)_      |
  //   |     /\            |    /\     |
  //   |    (  )      __   |   (  )    |
  //   |_(X)_\/            |    \/_(Y)_|
  //   |                   |
  //
  // For formula and details of implementation, see Operator::comm222_phss. The only change
  // made here to accommodate a dagger operator is to replace XY - YX  with -YX. This
  // ensures that we don't include contributions from the Qspace orbit to X.
  // Adapted from
  //    void Operator::comm222_phss( const Operator& X, const Operator& Y )

  void comm433sd_ph(const Operator &X, const Operator &Y, Operator &Z)
  {

    double t_start = omp_get_wtime();

    int hx = X.IsHermitian() ? 1 : -1; // assuming X is either hermitian or antihermitian. I'm sure this will come back to bite me some day.

    // Construct the intermediate matrix Z_bar and fill it with zeros.
    //   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup( Z.GetJRank(), Z.GetTRank(), Z.GetParity() );
    size_t nch = Z.modelspace->SortedTwoBodyChannels_CC.size();
    size_t norb = Z.modelspace->GetNumberOrbits();
    std::deque<arma::mat> Z_bar(Z.nChannels);
    std::vector<bool> lookup_empty(Z.nChannels, true);
    for (size_t ich = 0; ich < nch; ++ich)
    {
      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      index_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros(norb, 2 * nKets_cc); // Z_iQ`kj`   we want both orderings of k and j, but Q is fixed.
                                           //      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
      lookup_empty[ich] = false;
    }

// loop over channels for Z_bar
#ifndef OPENBLAS_NOUSEOMP
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar_transform_first_pass)
#endif
    for (size_t ich = 0; ich < nch; ++ich)
    {
      if (lookup_empty.at(ich))
        continue;
      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC.at(ich);
      const TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      size_t nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

      arma::mat Y_bar_ph;
      arma::mat X_bar_ph;

      DoPandyaTransformation_SingleChannel_Dagger(Y, Y_bar_ph, ch);    // Generate  YbarJ'_iQ`ab` * (na-nb)      for a<=b and a>b
      DoPandyaTransformation_SingleChannel(X, X_bar_ph, ch, "normal"); // Generate XbarJ'_ab`kj`                  for a<=b and a>b
      auto &Zbar_ch = Z_bar.at(ch);                                    // Marginally useful aliasing to avoid lookups...

      // Zbar will be
      // ZbarJ'_iQ`kj` = Sum_ab   (na-nb) <iQ`|YbarJ'|ab`> <ab`|XbarJ'|kj`>
      // So then we can perform an inverse Pandya transform (which is really just another Pandya transform) to get
      // ZJ_ijkQ = - Sum_J' (2J'+1) { ji  jj  J  }  ZbarJ'_iQ`kj`
      //                            { jk  jQ  J' }

      // The shape of Xt_bar_ph is  ( rows, columns) = ( 2*nph_kets,  nKets_cc)  =>    [     Xbar    ]  |  ab`, element of nph_kets, which has ph` and hp`
      //                                                                               [             ]  v
      //                                                                                    ->
      //                                                                                    kj` is an element of kets_cc
      //
      // The shape of Y_bar_ph is   (rows, columns) = (norb, 2*nph_kets)         =>    [     Ybar    ]  | iQ`  runs over all orbits i, with Q fixed.
      //                                                                               [             ]  v
      //                                                                                    ->
      //                                                                                    ab` is an element of nph_kets, here we have ph` and hp`
      //
      // So we should multiply Ybar * Xbar to sum over ab` and get out Zbar      =>      [     Zbar    ]  |  iQ`
      //                                                                                 [             ]  v
      //                                                                                      ->
      //                                                                                      kj`
      //
      // Really, we need  <iQ`|Zbar|kj`> and <iQ`|Zbar|jk`>

      if (Y_bar_ph.size() < 1 or X_bar_ph.size() < 1) // for an armadillo matrix, .size()  gives the total number of elements, i.e. n_rows * n_cols
      {
        //        Zbar_ch = arma::zeros( norb, 2*nKets_cc );  // This seems unnecessary since we initialized things earlier... try getting rid of it and see if things break.
        continue;
      }

      // get the phases for taking the transpose
      arma::mat PhaseMat(nKets_cc, nKets_cc, arma::fill::ones);
      for (index_t iket = 0; iket < nKets_cc; iket++)
      {
        const Ket &ket = tbc_cc.GetKet(iket);
        if (Z.modelspace->phase((ket.op->j2 + ket.oq->j2) / 2) > 0)
          continue;
        PhaseMat.col(iket) *= -1;
        PhaseMat.row(iket) *= -1;
      }
      arma::uvec phkets = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph());
      auto PhaseMatX = PhaseMat.rows(phkets) * hx;

      //                                           [      |     ]
      //     create full Y matrix from the half:   [  Yhp | Y'ph]   where the prime indicates multiplication by (-1)^(i+j+k+l) h_x
      //                                           [      |     ]   Flipping hp <-> ph and multiplying by the phase is equivalent to
      //                                           [  Yph | Y'hp]   having kets |kj> with k>j.
      Zbar_ch = Y_bar_ph * join_horiz(X_bar_ph, join_vert(X_bar_ph.tail_rows(nph_kets) % PhaseMatX,
                                                          X_bar_ph.head_rows(nph_kets) % PhaseMatX));
    }

    Z.profiler.timer["Build Z_bar_Dagger"] += omp_get_wtime() - t_start;
    //   std::cout << "Done building Z_bar_dagger" << std::endl;

    // Perform inverse Pandya transform on Z_bar to get Z
    t_start = omp_get_wtime();
    AddInversePandyaTransformation_Dagger(Z_bar, Z);

    //   Z.modelspace->scalar_transform_first_pass = false;
    Z.profiler.timer["InversePandyaTransformation_Dagger"] += omp_get_wtime() - t_start;

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //**************************************************************************
  // DAGGER VARIETY
  //
  //  X^J_iQ`ab` = - sum_J' { i Q J } (2J'+1) X^J'_ibaQ
  //                        { a b J'}
  /// A modification of the Pandya transform for the case of a dagger operator.
  /// We embed the dagger operator as a scalar with a ficticious fourth index Q,
  /// which is outside the Hilbert space. So our relevant transformation is
  /// \f[
  ///  \bar{X}^{J}_{i\bar{Q}a\bar{b}} = - \sum_{J'} (2J'+1)
  ///  \left\{ \begin{array}{lll}
  ///  j_i  &  j_Q  &  J \\
///  j_a  &  j_b  &  J' \\
///  \end{array} \right\}
  ///  X^{J}_{ibaQ}
  /// \f]
  /// where the overbar indicates time-reversed orbits.
  /// Since we use this for the commutator expression, the k and j orbits will be a and b,
  /// i.e. the orbits we sum over with a factor (na-nb). We include the (na-nb) factor
  /// in the output matrix.
  void DoPandyaTransformation_SingleChannel_Dagger(const Operator &Z, arma::mat &TwoBody_CC_ph, int ch_cc)
  {
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    //   int nKets_cc = tbc_cc.GetNumberKets();
    size_t norb = Z.modelspace->GetNumberOrbits();
    arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph());
    int nph_kets = kets_ph.n_rows;
    int J_cc = tbc_cc.J;

    //   TwoBody_CC_ph.zeros( nKets_cc, 2*nph_kets);
    TwoBody_CC_ph.zeros(norb, 2 * nph_kets); // factor 2 because we want both orderings ab and ba.

    index_t Q = Z.GetQSpaceOrbit();
    Orbit &oQ = Z.modelspace->GetOrbit(Q);
    double jQ = oQ.j2 * 0.5;

    // loop over cross-coupled ph kets |ab> in this channel
    // (this is the side that gets summed over in the matrix multiplication)
    for (int iket = 0; iket < nph_kets; ++iket)
    {
      Ket &ket_cc = tbc_cc.GetKet(kets_ph[iket]);
      int a = ket_cc.p;
      int b = ket_cc.q;
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Orbit &ob = Z.modelspace->GetOrbit(b);
      double ja = oa.j2 * 0.5;
      double jb = ob.j2 * 0.5;

      // Here we make some lists so that we don't need to repeat code when switching (a,b) -> (b,a).
      // Instead we can just iterate over the two cases.
      std::vector<int> ab = {a, b};
      std::vector<double> jab = {ja, jb};
      std::vector<double> nanb = {(oa.occ - ob.occ), (ob.occ - oa.occ)};

      // we loop over both orderings, a<b and a>b. Here, this is done by exchanging a<->b and taking a minus sign due to the (na-nb) factor.
      for (int ab_case = 0; ab_case <= 1; ab_case++)
      {
        size_t ab1 = ab[ab_case];
        size_t ab2 = ab[1 - ab_case];
        double jab1 = jab[ab_case];
        double jab2 = jab[1 - ab_case];
        size_t indx_ab = iket + ab_case * nph_kets;

        // loop over orbits i, and check if <iQ| is in the desired channel.
        for (size_t i = 0; i < norb; i++)
        {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          if (not tbc_cc.CheckChannel_ket(&oi, &oQ))
            continue;         //  <iQ|  isn't in this channel, so move along. (Note, for this check, the ordering of i,Q doesn't matter).
          size_t indx_iQ = i; // since Q is fixed, we can label the bra <iQ| by the index i.
          double ji = oi.j2 * 0.5;

          int Jmin = std::max(std::abs(jab1 - jQ), std::abs(ji - jab2));
          int Jmax = std::min(jab1 + jQ, ji + jab2);
          double Xbar = 0;
          for (int J_std = Jmin; J_std <= Jmax; ++J_std)
          {
            double sixj = Z.modelspace->GetSixJ(ji, jQ, J_cc, jab1, jab2, J_std);
            if (std::abs(sixj) < 1e-8)
              continue;
            //             double tbme = Z.TwoBody.GetTBME_J(J_std,i,ab2,ab1,Q);
            double tbme = Z.ThreeLeg.GetME_J(J_std, i, ab2, ab1);
            Xbar -= (2 * J_std + 1) * sixj * tbme;
          }
          TwoBody_CC_ph(indx_iQ, indx_ab) = Xbar * nanb[ab_case];
        }
      }
    }
  }

  // Same idea as the scalar commutator version, with enough modifications to make it not wrong (I hope).
  // One big difference is that the Pandya-transformed operator Zbar has its bra indices numbered
  // by orbit index rather than ket index, because it's   <iQ|Zbar|kj>  and Q is fixed.
  // Otherwise, it should look pretty similar.
  void AddInversePandyaTransformation_Dagger(const std::deque<arma::mat> &Zbar, Operator &Z)
  {
    // Do the inverse Pandya transform
    int n_nonzeroChannels = Z.modelspace->SortedTwoBodyChannels.size();
    size_t norb = Z.modelspace->GetNumberOrbits();
    size_t Q = Z.GetQSpaceOrbit();
    Orbit &oQ = Z.modelspace->GetOrbit(Q);
    double jQ = 0.5 * oQ.j2;

    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
    //   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
    for (int ich = 0; ich < n_nonzeroChannels; ++ich)
    {
      size_t ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      size_t nKets = tbc.GetNumberKets();

      // the bra is the <ij| part of Z_ijkQ
      for (size_t ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);

        // Two-element lists for use below, because we need to antisymmetrize with 1-Pij, including a phase factor (-1)^(ji+jj-J)
        // The code would look like the exact same thing copy-pasted twice, with i and j exchanged in the second copy.
        std::vector<size_t> ij_switcheroo = {bra.p, bra.q};
        std::vector<int> phaseij = {+1, bra.Phase(J)};

        double norm_ij = (bra.p == bra.q) ? PhysConst::INVSQRT2 : 1.0;

        for (size_t k = 0; k < norb; k++)
        {
          Orbit &ok = Z.modelspace->GetOrbit(k);
          //          if ( not tbc.CheckChannel_ket(&ok, &oQ) ) continue; // if |kQ> doesn't live in this channel, move along.
          //          if ( ( (ok.tz2+oQ.tz2)!=2*tbc.Tz) or ( (ok.l+oQ.l)%2!=tbc.parity) or ( (ok.j2+oQ.j2)<2*tbc.J) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue; // if |kQ> doesn't live in this channel, move along.
          if (((ok.tz2 - oQ.tz2) != 2 * tbc.Tz) or ((ok.l + oQ.l) % 2 != tbc.parity) or ((ok.j2 + oQ.j2) < 2 * tbc.J) or (std::abs(ok.j2 - oQ.j2) > 2 * tbc.J))
            continue; // if |kQ> doesn't live in this channel, move along.

          double jk = ok.j2 / 2.;
          //          double norm_kQ = (k==Q) ? PhysConst::INVSQRT2  : 1.0;
          double norm_kQ = 1.0;

          double commijk = 0; // collect the sum in commijk

          for (int ij_case = 0; ij_case <= 1; ij_case++) // loop over the two cases corresponding to 1-Pij (with the appropriate phase on the exchange term)
          {
            size_t i = ij_switcheroo[ij_case];
            size_t j = ij_switcheroo[1 - ij_case];
            int Pij = phaseij[ij_case];

            Orbit &oi = Z.modelspace->GetOrbit(i);
            Orbit &oj = Z.modelspace->GetOrbit(j);
            double ji = oi.j2 / 2.;
            double jj = oj.j2 / 2.;

            // now loop over the Jprime entering in the Pandya transformation
            int parity_cc = (oi.l + oQ.l) % 2;
            //            int Tz_cc = std::abs(oi.tz2+oQ.tz2)/2;
            int Tz_cc = std::abs(oi.tz2 - oQ.tz2) / 2;
            int Jmin = std::max(std::abs(int(ji - jQ)), std::abs(int(jk - jj)));
            int Jmax = std::min(ji + jQ, jk + jj);
            for (int Jprime = Jmin; Jprime <= Jmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(ji, jj, J, jk, jQ, Jprime);
              if (std::abs(sixj) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              size_t indx_iQ = i;
              size_t indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k)) + (k > j ? nkets_cc : 0);
              double me1 = Zbar.at(ch_cc)(indx_iQ, indx_kj);
              commijk -= (2 * Jprime + 1.) * sixj * me1 * Pij;
            }

          } // for ij_case
          commijk *= norm_ij * norm_kQ;
          //            Z.TwoBody.AddToTBME_J(J,i,j,k,Q,  commijk );   // AddToTBME_J assumes a normalized matrix element.
          //            Z.ThreeLeg.AddToME_J(J,i,j,k,  commijk );   // AddToTBME_J assumes a normalized matrix element.
          Z.ThreeLeg.AddToME_J(J, bra.p, bra.q, k, commijk); // AddToTBME_J assumes a normalized matrix element.
        }                                                    // for k
      }                                                      // for ibra, which is <ij|
    }                                                        // for ich, which is J etc
  }



}
