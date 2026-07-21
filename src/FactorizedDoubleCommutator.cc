
#include "FactorizedDoubleCommutator.hh"
#include "Commutator.hh"
#include "PhysicalConstants.hh"

namespace Commutator
{

  namespace FactorizedDoubleCommutator
  {

    bool use_goose_tank_1b = true; // always include
    bool use_goose_tank_2b = true; // always include

    bool use_1b_intermediates = true;
    bool use_2b_intermediates = true;
    bool use_goose_tank_only_1b = false; // only calculate Goose Tanks
    bool use_goose_tank_only_2b = false; // only calculate Goose Tanks
    bool use_TypeII_1b = true;
    bool use_TypeIII_1b = true;
    bool use_TypeII_2b = true;
    bool use_TypeIII_2b = true;
    bool use_GT_TypeI_2b = true;
    bool use_GT_TypeIV_2b = true;
    //  bool SlowVersion = false;

    void SetUse_GooseTank_1b(bool tf)
    {
      use_goose_tank_1b = tf;
    }

    void SetUse_GooseTank_2b(bool tf)
    {
      use_goose_tank_2b = tf;
    }

    void SetUse_1b_Intermediates(bool tf)
    {
      use_1b_intermediates = tf;
    }
    void SetUse_2b_Intermediates(bool tf)
    {
      use_2b_intermediates = tf;
    }

    void SetUse_GooseTank_only_1b(bool tf)
    {
      use_goose_tank_only_1b = tf;
    }

    void SetUse_GooseTank_only_2b(bool tf)
    {
      use_goose_tank_only_2b = tf;
    }

    void SetUse_TypeII_1b(bool tf)
    {
      use_TypeII_1b = tf;
    }

    void SetUse_TypeIII_1b(bool tf)
    {
      use_TypeIII_1b = tf;
    }

    void SetUse_TypeII_2b(bool tf)
    {
      use_TypeII_2b = tf;
    }

    void SetUse_TypeIII_2b(bool tf)
    {
      use_TypeIII_2b = tf;
    }

    void SetUse_GT_TypeI_2b(bool tf)
    {
      use_GT_TypeI_2b = tf;
    }

    void SetUse_GT_TypeIV_2b(bool tf)
    {
      use_GT_TypeIV_2b = tf;
    }

    //  void UseSlowVersion(bool tf)
    //  {
    //    SlowVersion = tf;
    //  }

    // factorize double commutator [Eta, [Eta, Gamma]_3b ]_1b
    void comm223_231(const Operator &Eta, const Operator &Gamma, Operator &Z)
    {

      // Do some extra work to check if the operators being passed in are reduced,
      // because the comm223_231 routines assume not-reduced matrix elements
      const Operator * Etanred   = &Eta;  // Pointer to the non-reduced version of the operator
      const Operator * Gammanred = &Gamma;
      Operator Etatmp,Gammatmp;  // Declared, but not yet allocated, for reasons of scope.
      if ( Eta.IsReduced() ) // comm223_231 doesn't expect reduced operators. Need to make it not reduced.
      {
         Etatmp = Eta; // Actual copy, not a reference
         Etatmp.MakeNotReduced();
         Etanred = &Etatmp; // Now Etanred points to the not-reduced copy Etatmp
      }
      if ( Gamma.IsReduced() )
      {
         Gammatmp = Gamma;
         Gammatmp.MakeNotReduced();
         Gammanred = &Gammatmp;  // Now Gammanred points to the not-reduced copy Gammatmp
      }

      bool z_was_reduced = Z.IsReduced();
      if (z_was_reduced) Z.MakeNotReduced();


      if (use_1b_intermediates)
      {
        comm223_231_chi1b(*Etanred, *Gammanred, Z); // topology with 1-body intermediate (fast)
      }
      if (use_2b_intermediates)
      {
        comm223_231_chi2b(*Etanred, *Gammanred, Z); // topology with 2-body intermediate (slow)
      }
      if (z_was_reduced) Z.MakeReduced();

      return;
    } // comm223_231

    ////////////////////////////////////////////////////////////////////////////
    /// factorized 223_231 double commutator with 1b intermediate
    ////////////////////////////////////////////////////////////////////////////
    void comm223_231_chi1b(const Operator &Eta, const Operator &Gamma, Operator &Z)
    {

      double t_internal = omp_get_wtime(); // timer
      double t_start = omp_get_wtime();    // timer

      Z.modelspace->PreCalculateSixJ();

      // determine symmetry
      int hEta = Eta.IsHermitian() ? 1 : -1;
      int hGamma = Gamma.IsHermitian() ? 1 : -1;
      // int hZ = Z.IsHermitian() ? 1 : -1;
      int hZ = hGamma;
      // ###########################################################
      //  diagram I
      // The intermediate one body operator
      //  Chi_221_a :
      //          eta | d
      //         _____|
      //       /\     |
      //   a  (  ) b  | c
      //       \/_____|
      //          eta |
      //              | e
      // Chi_221_a = \sum \hat(J_0) ( nnnn - ... ) eta eta

      auto Chi_221_a = Z.OneBody;
      Chi_221_a.zeros(); // Set all elements to zero

      int nch = Z.modelspace->GetNumberTwoBodyChannels();
      int norbits = Z.modelspace->all_orbits.size();
      std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());

//      TwoBodyME intermediateTB = Z.TwoBody;
      TwoBodyME intermediateTB = Eta.TwoBody;
      intermediateTB.Erase();

// full matrix
#pragma omp parallel for schedule(dynamic, 1)
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();

        arma::mat Eta_matrix = Eta.TwoBody.GetMatrix(ch, ch);
        arma::mat Eta_matrix_nnnn = Eta_matrix;

        for (int ibra = 0; ibra < nKets; ++ibra)
        {
          Ket &bra = tbc.GetKet(ibra);
          double n_i = bra.op->occ;
          double n_j = bra.oq->occ;

          for (int iket = ibra; iket < nKets; ++iket)
          {
            Ket &ket = tbc.GetKet(iket);
            double n_k = ket.op->occ;
            double n_l = ket.oq->occ;
            double occfactor = n_i * n_j * (1 - n_k) * (1 - n_l) - (1 - n_i) * (1 - n_j) * n_k * n_l;

            Eta_matrix_nnnn(ibra, iket) *= occfactor;
            Eta_matrix_nnnn(iket, ibra) *= -occfactor;

          } // for iket
        } // for ibra

        arma::mat tmp = 2 * (2 * J0 + 1) * Eta_matrix_nnnn * Eta_matrix;
        //        tmp += tmp.t();

        intermediateTB.GetMatrix(ch, ch) = tmp + tmp.t();
      } // for ch

      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_intermediateTB"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

#pragma omp parallel for schedule(dynamic, 1)
      for (int indexd = 0; indexd < norbits; ++indexd)
      {
        auto d = allorb_vec[indexd];
        Orbit &od = Z.modelspace->GetOrbit(d);
        double dj2 = od.j2 + 1.0;
//        for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
        for (auto &e : Eta.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
//        for (auto e : Z.modelspace->all_orbits) // delta_jd je
        {
          if (e > d)
            continue;
          Orbit &oe = Z.modelspace->GetOrbit(e);
//          if (oe.j2 != od.j2) continue;
          double eta_de = 0;

          for (auto &b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);

            int Jmin = std::abs(od.j2 - ob.j2) / 2;
            int Jmax = (od.j2 + ob.j2) / 2;
            int dJ = 1;
            if (d == b or e == b)
            {
              dJ = 2;
              Jmin += Jmin % 2;
            }
            for (int J0 = Jmin; J0 <= Jmax; J0++)
            {
              eta_de += intermediateTB.GetTBME_J(J0, J0, b, d, b, e);
            }
          }
          Chi_221_a(d, e) += eta_de / dj2;
          if (d != e)
            Chi_221_a(e, d) += eta_de / dj2;
        } // e
      } // d
//      std::cout << "Chi_221_a: " << std::endl << Chi_221_a << std::endl;

      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_indexd"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

#pragma omp parallel for schedule(dynamic, 1)
      for (int indexp = 0; indexp < norbits; ++indexp)
      {
        auto p = allorb_vec[indexp];
        Orbit &op = Z.modelspace->GetOrbit(p);
        for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
        {
          if (q > p)
            continue;
          Orbit &oq = Z.modelspace->GetOrbit(q);
          double zij = 0;

          for (auto &d : Z.modelspace->all_orbits)
          {
            Orbit &od = Z.modelspace->GetOrbit(d);
//            for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
//            for (auto e : Z.modelspace->all_orbits) // delta_jd je
            for (auto &e : Eta.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
            {
              Orbit &oe = Z.modelspace->GetOrbit(e);
              if ( oe.j2 != od.j2 ) continue;

              int J1min = std::abs(od.j2 - oq.j2) / 2;
              int J1max = (od.j2 + oq.j2) / 2;
              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                zij += (2 * J1 + 1) * Chi_221_a(d, e) * Gamma.TwoBody.GetTBME_J(J1, J1, e, p, d, q);
              }
            }
          }

          Z.OneBody(p, q) += 0.5 * zij / (op.j2 + 1.0);
          if (p != q)
            Z.OneBody(q, p) += 0.5 * hZ * zij / (op.j2 + 1.0);
          //--------------------------------------------------
        } // for q
      } // for p
//      std::cout << "fI : " << std::endl << Z.OneBody << std::endl;

      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_indexp"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      if (use_goose_tank_only_1b)
      {
        Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
        return;
      }

      // *********************************************************************************** //
      //                                  Diagram III                                        //
      // *********************************************************************************** //

      // ###########################################################
      //  diagram III_a and diagram III_b
      // ###########################################################WW
      // The one body operator
      //  Chi_221_b :
      //          eta | d
      //         _____|
      //       /\     |
      //   a  (  ) b  | c
      //       \/_____|
      //        Gamma |
      //              | e
      // Chi_221_b = \sum \hat(J_0) ( nnnn - ... ) eta Gamma
      // non-Hermit
      auto Chi_221_b = Z.OneBody;
      Chi_221_b.zeros(); // Set all elements to zero

#pragma omp parallel for
      for (int indexd = 0; indexd < norbits; ++indexd)
      {
        auto d = allorb_vec[indexd];
        Orbit &od = Z.modelspace->GetOrbit(d);
        double n_d = od.occ;
        double nbar_d = 1.0 - n_d;

//        for (auto e : Z.modelspace->all_orbits) // delta_jd je
        for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
        {
          // if (e > d)
          //   continue;
          Orbit &oe = Z.modelspace->GetOrbit(e);
//          if ( oe.j2 != od.j2 ) continue;
          double n_e = oe.occ;
          double nbar_e = 1.0 - n_e;

          double eta_de = 0;

          for (int ch = 0; ch < nch; ++ch)
          {
            TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
            int J0 = tbc.J;
            int nKets = tbc.GetNumberKets();
            for (int ibra = 0; ibra < nKets; ++ibra)
            {
              Ket &bra = tbc.GetKet(ibra);
              int b = bra.p;
              int c = bra.q;

              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;

              Orbit &oc = Z.modelspace->GetOrbit(c);
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &a : Z.modelspace->all_orbits)
              {
                Orbit &oa = Z.modelspace->GetOrbit(a);
                double n_a = oa.occ;
                double nbar_a = 1.0 - n_a;

                double occfactor = (nbar_a * nbar_e * n_b * n_c - nbar_b * nbar_c * n_a * n_e);
                if (std::abs(occfactor) < 1e-6)
                  continue;
                double MEs = (2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, c, a, e) * Gamma.TwoBody.GetTBME_J(J0, J0, a, d, b, c);
                eta_de += MEs;
                if (b != c)
                  eta_de += MEs;
              }
            }
          }
          Chi_221_b(d, e) += eta_de / (od.j2 + 1.0);
        } // e
      } // d
//      std::cout << "Chi_221_b " << std::endl << Chi_221_b << std::endl;

      //  diagram III_a and diagram III_b together
#pragma omp parallel for
      for (int indexd = 0; indexd < norbits; ++indexd)
      {
        auto p = allorb_vec[indexd];
        Orbit &op = Z.modelspace->GetOrbit(p);
        for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
        {
          if (q > p)
            continue;
          Orbit &oq = Z.modelspace->GetOrbit(q);
          double zij_a = 0;
          double zij_b = 0;
          // loop abcde
          for (auto &d : Z.modelspace->all_orbits)
          {
            Orbit &od = Z.modelspace->GetOrbit(d);

            for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
            {
              Orbit &oe = Z.modelspace->GetOrbit(e);

              int J1min = std::abs(oe.j2 - op.j2) / 2;
              int J1max = (oe.j2 + op.j2) / 2;

              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                double etaME = (2 * J1 + 1) * Eta.TwoBody.GetTBME_J(J1, J1, e, p, d, q);
                zij_a += Chi_221_b(d, e) * etaME;
                zij_b += hZ * Chi_221_b(e, d) * etaME;
              }
            }
          }

          Z.OneBody(p, q) += 0.5 * (zij_a - zij_b) / (op.j2 + 1.0);
          if (p != q)
            Z.OneBody(q, p) += 0.5 * hZ * (zij_a - zij_b) / (op.j2 + 1.0);
          //--------------------------------------------------
        } // for q
      } // for p
//      std::cout << "fII : " << std::endl << Z.OneBody << std::endl;

      Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
      return;
    } // comm223_231_chi1b



// PREVIOUS IMPLEMENTATION. SAVING THIS IN CASE I SCREW THINGS UP WHILE MODIFYING IT.
/*
    ////////////////////////////////////////////////////////////////////////////
    /// factorized 223_231 double commutator with 2b intermediate
    ////////////////////////////////////////////////////////////////////////////
    void comm223_231_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z)
    {

      double t_internal = omp_get_wtime(); // timer
      double t_start = omp_get_wtime();    // timer

      Z.modelspace->PreCalculateSixJ();

      // determine symmetry
      int hEta = Eta.IsHermitian() ? 1 : -1;
      int hGamma = Gamma.IsHermitian() ? 1 : -1;
      // int hZ = Z.IsHermitian() ? 1 : -1;
      int hZ = hGamma;

      // *********************************************************************************** //
      //                               Diagram II_b and II_d                                 //
      // *********************************************************************************** //
      // ###########################################################
      //  diagram II_b and II_d
      //
      // The two body operator
      //  Chi_222_b :
      //        c  | eta |  p
      //           |_____|
      //        b  |_____|  e
      //           | eta |
      //        a  |     |  d
      // ###########################################################

      int nch = Z.modelspace->GetNumberTwoBodyChannels();
      int norbits = Z.modelspace->all_orbits.size();

      TwoBodyME intermediateTB = Z.TwoBody;
      intermediateTB.Erase();

      std::vector<int> bra_channels;
      std::vector<int> ket_channels;

      for (auto &itmat : Z.TwoBody.MatEl)
      {
        bra_channels.push_back(itmat.first[0]);
        ket_channels.push_back(itmat.first[1]);
      }
      int nbra_ket_ch = bra_channels.size();

// fill  Gamma_matrix Eta_matrix, Eta_matrix_nnnn
#pragma omp parallel for schedule(dynamic, 1)
//      for (int ch = 0; ch < nch; ++ch)
      for (int ich=0; ich<nbra_ket_ch; ich++)
      {
        size_t ch_bra = bra_channels[ich];
        size_t ch_ket = ket_channels[ich];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int Jbra = tbc_bra.J;
        int nBras = tbc_bra.GetNumberKets();
        int nKets = tbc_ket.GetNumberKets();

        const arma::mat& Eta_mat_bra = Eta.TwoBody.GetMatrix(ch_bra, ch_bra);
        const arma::mat& Eta_mat_ket = Eta.TwoBody.GetMatrix(ch_ket, ch_ket);
        arma::mat Eta_mat_nnnn_bra = Eta_mat_bra;
        arma::mat Eta_mat_nnnn_ket = Eta_mat_ket;
        const arma::mat& Gamma_mat = Gamma.TwoBody.GetMatrix(ch_bra, ch_ket);

//        arma::mat Eta_matrix = Eta.TwoBody.GetMatrix(ch, ch);
//        arma::mat Eta_matrix_nnnn = Eta_matrix;
//        arma::mat Gamma_matrix = Gamma.TwoBody.GetMatrix(ch, ch);

        for (int ibra = 0; ibra < nBras; ++ibra)
        {
          Ket &bra = tbc_bra.GetKet(ibra);
          double n_i = bra.op->occ;
          double n_j = bra.oq->occ;

          for (int iket = 0; iket < nBras; ++iket)
          {
            Ket &ket = tbc_bra.GetKet(iket);
            double n_k = ket.op->occ;
            double n_l = ket.oq->occ;
            double occfactor = n_i * n_j * (1 - n_k) * (1 - n_l) - (1 - n_i) * (1 - n_j) * n_k * n_l;

            Eta_mat_nnnn_bra(ibra, iket) *= occfactor;
          } // for iket
        } // for ibra
        for (int ibra = 0; ibra < nKets; ++ibra)
        {
          Ket &bra = tbc_ket.GetKet(ibra);
          double n_i = bra.op->occ;
          double n_j = bra.oq->occ;

          for (int iket = 0; iket < nKets; ++iket)
          {
            Ket &ket = tbc_ket.GetKet(iket);
            double n_k = ket.op->occ;
            double n_l = ket.oq->occ;
            double occfactor = n_i * n_j * (1 - n_k) * (1 - n_l) - (1 - n_i) * (1 - n_j) * n_k * n_l;

            Eta_mat_nnnn_ket(ibra, iket) *= occfactor;
          } // for iket
        } // for ibra

//        arma::mat Chi_222_b = 4 * (2 * J0 + 1) * Eta_matrix * Eta_matrix_nnnn * Gamma_matrix;
        arma::mat Chi_222_b =  Eta_mat_bra * Eta_mat_nnnn_bra * Gamma_mat;
        if ( ch_bra==ch_ket)
        {
          Chi_222_b += Chi_222_b.t();
        }
        else
        {
           Chi_222_b -= Gamma_mat * Eta_mat_nnnn_ket * Eta_mat_ket;
        }
        Chi_222_b *= 4 * (2 * Jbra + 1) ;

//        Chi_222_b += Chi_222_b.t();

//        intermediateTB.GetMatrix(ch, ch) = Chi_222_b;
        intermediateTB.GetMatrix(ch_bra, ch_ket) = Chi_222_b;

      } // for ch

      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_chi2_pp_fill_chi"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      //  diagram II_b and II_d
      //
      // IIb_pq = 1/4 1/(2 jp + 1) \sum_acdJ0 Chi_222_b_cpad * Gamma_bar_adcq
      // IIb_pq = - 1/4 1/(2 jp + 1) \sum_abeJ0  Chi_222_a_bqae * Gamma_bar_bqae
      // ###########################################################
      std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
      for (int indexp = 0; indexp < norbits; ++indexp)
      {
        auto p = allorb_vec[indexp];
        Orbit &op = Z.modelspace->GetOrbit(p);
        double jp = op.j2 / 2.;
        for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
        {
          if (q > p)
            continue;
          Orbit &oq = Z.modelspace->GetOrbit(q);
          double jq = oq.j2 / 2.;
          double zpq = 0;

          // loop abcde
          for (auto &c : Z.modelspace->all_orbits)
          {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double jc = oc.j2 / 2.;

            int J0min = std::abs(oc.j2 - op.j2) / 2;
            int J0max = (oc.j2 + op.j2) / 2;
            for (int J0 = J0min; J0 <= J0max; J0++)
            {
              zpq += intermediateTB.GetTBME_J(J0, J0, c, p, c, q);
            }
          }
          Z.OneBody(p, q) += 0.25 * zpq / (op.j2 + 1.0);
          if (p != q)
            Z.OneBody(q, p) += 0.25 * hZ * zpq / (op.j2 + 1.0);
          //--------------------------------------------------
        } // for q
      } // for p

      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_chi2_pp_fill1b"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      // *********************************************************************************** //
      //                               Diagram II_a and II_c                                 //
      // *********************************************************************************** //

      // ###########################################################
      //  diagram II_a
      //
      //  Pandya transform
      //  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
      //                        { k l J'}
      int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();

      // The two body operator
      //  Chi_222_a :
      //            eta |
      //           _____|
      //          /\    |
      //   |     (  )
      //   |_____ \/
      //   | eta
      //
      //  Chi_222_a = \sum_pq (nbar_e * nbar_d * n_f * n_c - nbar_f * nbar_c * n_e * n_d )

      // TODO: Modify this to handle Tz-changing operators.
      // To do this, the intermediate object ~ Eta * Eta * Gamma should be generalized to
      // allow different CC channels from the bra and the ket. The Eta * Eta part will conserve
      // the channel still. Then in the loop over p,q,e below, need to find the CC channels
      // for the bra and the ket in the look up. Shouldn't be too terrible.
      std::deque<arma::mat> IntermediateTwobody(n_nonzero);
/// Pandya transformation
#pragma omp parallel for schedule(dynamic, 1)
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        int J_cc = tbc_cc.J;

        // We don't need to double the rows of the Eta_bar matrix, since these lead to rows in
        // the IntermediateTwobody matrix which can be obtained by symmetry from the top rows.
        // So this will save us some time, and when filling the resuliting one-body matrix we just
        // ensure we always access with the ordering that is stored. -SRS
        arma::mat Eta_bar = arma::mat(nKets_cc, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
        // arma::mat Eta_bar = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
        arma::mat Eta_bar_nnnn = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
        arma::mat Gamma_bar = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);    // SRS ADDED
                                                                                           //      Gamma_bar[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
        arma::mat Chi_222_a;
        //      IntermediateTwobody[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        // transform operator
        // loop over cross-coupled ph bras <ab| in this channel
        //      for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc)
        for (int ibra_cc = 0; ibra_cc < nKets_cc; ++ibra_cc)
        {
          int a, b;
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
          a = bra_cc.p;
          b = bra_cc.q;
          // if (ibra_cc >= nKets_cc and a == b)
          //   continue;

          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          double nbar_a = 1 - n_a;
          double ja = oa.j2 * 0.5;

          Orbit &ob = Z.modelspace->GetOrbit(b);
          double n_b = ob.occ;
          double nbar_b = 1 - n_b;
          double jb = ob.j2 * 0.5;

          // loop over cross-coupled kets |cd> in this channel
          for (int iket_cc = ibra_cc; iket_cc < nKets_cc * 2; ++iket_cc)
          {
            if ((iket_cc % nKets_cc) < ibra_cc)
              continue; // We'll get these from symmetry
            int c, d;
            if (iket_cc < nKets_cc)
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
              c = ket_cc_cd.p;
              d = ket_cc_cd.q;
            }
            else
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
              d = ket_cc_cd.p;
              c = ket_cc_cd.q;
            }

            Orbit &oc = Z.modelspace->GetOrbit(c);
            double n_c = oc.occ;
            double nbar_c = 1 - n_c;
            double jc = oc.j2 * 0.5;

            Orbit &od = Z.modelspace->GetOrbit(d);
            double n_d = od.occ;
            double nbar_d = 1 - n_d;
            double jd = od.j2 * 0.5;

            double occ_factor = nbar_c * nbar_b * n_a * n_d - n_c * n_b * nbar_a * nbar_d;

            // Check the isospin projection. If this isn't conserved in the usual channel,
            // then all the xcbad and yadcb will be zero and we don't need to bother computing SixJs.
            if (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != Gamma.GetTRank() and std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != Eta.GetTRank())
              continue;

            int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
            int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
            double Xbar = 0;
            double Ybar = 0;
            int dJ_std = 1;
            if ((a == d or b == c))
            {
              dJ_std = 2;
              jmin += jmin % 2;
            }
            for (int J_std = jmin; J_std <= jmax; J_std += dJ_std)
            {

              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                Xbar -= (2 * J_std + 1) * sixj1 * Eta.TwoBody.GetTBME_J(J_std, a, d, c, b);
                Ybar -= (2 * J_std + 1) * sixj1 * Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
              }
            }
            double flip_phase = Z.modelspace->phase((oa.j2 + ob.j2 + oc.j2 + od.j2) / 2);

            if (iket_cc < nKets_cc or (iket_cc >= nKets_cc and c != d))
            {
              // direct term
              Gamma_bar(ibra_cc, iket_cc) = Ybar;
              Eta_bar(ibra_cc, iket_cc) = Xbar;
              Eta_bar_nnnn(ibra_cc, iket_cc) = Xbar * occ_factor;

              if (iket_cc != ibra_cc)
              {
                // Hermiticity: Xbar_cdab = hX * Xbar_abcd.  We get a minus sign on the occupation factor
                Gamma_bar(iket_cc, ibra_cc) = hGamma * Ybar;
                Eta_bar_nnnn(iket_cc, ibra_cc) = hEta * Xbar * (-occ_factor);
                if (iket_cc < nKets_cc)
                {
                  Eta_bar(iket_cc, ibra_cc) = hEta * Xbar;
                }
              }
            }

            if (a != b)
            {
              // By exchange symmetry Xbar_badc = phase * hX * Xbar_abcd.  For Eta_bar_nnnn, this also requires swapping labels in the occupations -> minus sign.
              Gamma_bar(ibra_cc + nKets_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Ybar * flip_phase * hGamma;
              Eta_bar_nnnn(ibra_cc + nKets_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Xbar * flip_phase * hEta * (-occ_factor);
              // Eta_bar( ibra_cc+nKets_cc, (iket_cc + nKets_cc)%(2*nKets_cc)) = Xbar * flip_phase * hEta;
            }

            if (iket_cc >= nKets_cc or (iket_cc < nKets_cc and c != d))
            {
              // Combined exchange symmetry and hermiticity
              // Xbar_dcba = phase * Xbar_abcd
              Gamma_bar((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Ybar * flip_phase;
              Eta_bar_nnnn((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Xbar * flip_phase * (occ_factor);
              if (iket_cc >= nKets_cc)
                Eta_bar((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Xbar * flip_phase;
            }

          } // for iket-cc

          //-------------------
        } // for ibra_cc

        IntermediateTwobody[ch_cc] = (2 * J_cc + 1) * Eta_bar * Eta_bar_nnnn * Gamma_bar;
        // IntermediateTwobody[ch_cc] = (2 * J_cc + 1) * Eta_bar * Gamma_bar;
      } // for ch_cc

      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_chi2_ph_fill_chi"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      // ###########################################################
      // diagram II_a
      //
      //  IIa_pq = 1/ (2 jp + 1) \sum_abeJ3 Chi_222_a_peab * Gamma_bar_abqe
      //
      // diagram II_c
      //
      //  IIc_pq = - 1/ (2 jp + 1) \sum_abe J3 Chi_222_a_eqab * Gamma_bar_abep
      // ###########################################################

// may be worth rolling p and q loops together for better load balancing
#pragma omp parallel for schedule(dynamic, 1)
      for (int indexd = 0; indexd < norbits; ++indexd)
      {
        auto p = allorb_vec[indexd];
        Orbit &op = Z.modelspace->GetOrbit(p);
        double jp = op.j2 / 2.;
        double j2hat2 = (op.j2 + 1.0);
        for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
        {
          if (q > p)
            continue;
          double zij = 0;
          for (auto &e : Z.modelspace->all_orbits) // delta_jp jq
          {
            Orbit &oe = Z.modelspace->GetOrbit(e);

            int Jtmin = std::abs(op.j2 - oe.j2) / 2;
            int Jtmax = (op.j2 + oe.j2) / 2;
            int parity_cc = (op.l + oe.l) % 2;
            int Tz_cc = std::abs(op.tz2 - oe.tz2) / 2;
            double zij = 0;
            for (int Jt = Jtmin; Jt <= Jtmax; Jt++)
            {

              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jt, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              // Make sure we access the element <ab|X|cd> with a<=b. If we want the other ordering, we get a minus sign.
              int ind_pe = tbc_cc.GetLocalIndex(p, e);
              int ind_qe = tbc_cc.GetLocalIndex(q, e);
              int ind_eq = tbc_cc.GetLocalIndex(e, q);
              int ind_ep = tbc_cc.GetLocalIndex(e, p);

              if (p <= e)
              {
                zij += IntermediateTwobody[ch_cc](ind_pe, ind_qe);
              }
              else
              {
                zij -= IntermediateTwobody[ch_cc](ind_ep, ind_eq);
              }

              if (e <= q)
              {
                zij -= IntermediateTwobody[ch_cc](ind_eq, ind_ep);
              }
              else
              {
                zij += IntermediateTwobody[ch_cc](ind_qe, ind_pe);
              }
            }
            Z.OneBody(p, q) += zij / j2hat2;
            if (p != q)
              Z.OneBody(q, p) += hZ * zij / j2hat2;
          }
        }
      }

      Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
      return;

    } // comm223_231_chi2b




*/



    ////////////////////////////////////////////////////////////////////////////
    /// factorized 223_231 double commutator with 2b intermediate
    ////////////////////////////////////////////////////////////////////////////
    void comm223_231_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z)
    {

      double t_internal = omp_get_wtime(); // timer
      double t_start = omp_get_wtime();    // timer

      Z.modelspace->PreCalculateSixJ();

      // determine symmetry
      int hEta = Eta.IsHermitian() ? 1 : -1;
      int hGamma = Gamma.IsHermitian() ? 1 : -1;
      // int hZ = Z.IsHermitian() ? 1 : -1;
      int hZ = hGamma;

      // *********************************************************************************** //
      //                               Diagram II_b and II_d                                 //
      // *********************************************************************************** //
      // ###########################################################
      //  diagram II_b and II_d
      //
      // The two body operator
      //  Chi_222_b :
      //        c  | eta |  p
      //           |_____|
      //        b  |_____|  e
      //           | eta |
      //        a  |     |  d
      // ###########################################################

      int nch = Z.modelspace->GetNumberTwoBodyChannels();
      int norbits = Z.modelspace->all_orbits.size();

      TwoBodyME intermediateTB = Z.TwoBody;
      intermediateTB.Erase();

      std::vector<int> bra_channels;
      std::vector<int> ket_channels;

      for (auto &itmat : Z.TwoBody.MatEl)
      {
        bra_channels.push_back(itmat.first[0]);
        ket_channels.push_back(itmat.first[1]);
      }
      int nbra_ket_ch = bra_channels.size();

// fill  Gamma_matrix Eta_matrix, Eta_matrix_nnnn
#pragma omp parallel for schedule(dynamic, 1)
//      for (int ch = 0; ch < nch; ++ch)
      for (int ich=0; ich<nbra_ket_ch; ich++)
      {
        size_t ch_bra = bra_channels[ich];
        size_t ch_ket = ket_channels[ich];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int Jbra = tbc_bra.J;
        int nBras = tbc_bra.GetNumberKets();
        int nKets = tbc_ket.GetNumberKets();

        const arma::mat& Eta_mat_bra = Eta.TwoBody.GetMatrix(ch_bra, ch_bra);
        const arma::mat& Eta_mat_ket = Eta.TwoBody.GetMatrix(ch_ket, ch_ket);
        arma::mat Eta_mat_nnnn_bra = Eta_mat_bra;
        arma::mat Eta_mat_nnnn_ket = Eta_mat_ket;
        const arma::mat& Gamma_mat = Gamma.TwoBody.GetMatrix(ch_bra, ch_ket);

//        arma::mat Eta_matrix = Eta.TwoBody.GetMatrix(ch, ch);
//        arma::mat Eta_matrix_nnnn = Eta_matrix;
//        arma::mat Gamma_matrix = Gamma.TwoBody.GetMatrix(ch, ch);

        for (int ibra = 0; ibra < nBras; ++ibra)
        {
          Ket &bra = tbc_bra.GetKet(ibra);
          double n_i = bra.op->occ;
          double n_j = bra.oq->occ;

          for (int iket = 0; iket < nBras; ++iket)
          {
            Ket &ket = tbc_bra.GetKet(iket);
            double n_k = ket.op->occ;
            double n_l = ket.oq->occ;
            double occfactor = n_i * n_j * (1 - n_k) * (1 - n_l) - (1 - n_i) * (1 - n_j) * n_k * n_l;

            Eta_mat_nnnn_bra(ibra, iket) *= occfactor;
          } // for iket
        } // for ibra
        for (int ibra = 0; ibra < nKets; ++ibra)
        {
          Ket &bra = tbc_ket.GetKet(ibra);
          double n_i = bra.op->occ;
          double n_j = bra.oq->occ;

          for (int iket = 0; iket < nKets; ++iket)
          {
            Ket &ket = tbc_ket.GetKet(iket);
            double n_k = ket.op->occ;
            double n_l = ket.oq->occ;
            double occfactor = n_i * n_j * (1 - n_k) * (1 - n_l) - (1 - n_i) * (1 - n_j) * n_k * n_l;

            Eta_mat_nnnn_ket(ibra, iket) *= occfactor;
          } // for iket
        } // for ibra

//        arma::mat Chi_222_b = 4 * (2 * J0 + 1) * Eta_matrix * Eta_matrix_nnnn * Gamma_matrix;
        arma::mat Chi_222_b =  Eta_mat_bra * Eta_mat_nnnn_bra * Gamma_mat;
        if ( ch_bra==ch_ket)
        {
          Chi_222_b += Chi_222_b.t();
        }
        else
        {
           Chi_222_b -= Gamma_mat * Eta_mat_nnnn_ket * Eta_mat_ket;
        }
        Chi_222_b *= 4 * (2 * Jbra + 1) ;

//        Chi_222_b += Chi_222_b.t();

//        intermediateTB.GetMatrix(ch, ch) = Chi_222_b;
        intermediateTB.GetMatrix(ch_bra, ch_ket) = Chi_222_b;

      } // for ch

      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_chi2_pp_fill_chi"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      //  diagram II_b and II_d
      //
      // IIb_pq = 1/4 1/(2 jp + 1) \sum_acdJ0 Chi_222_b_cpad * Gamma_bar_adcq
      // IIb_pq = - 1/4 1/(2 jp + 1) \sum_abeJ0  Chi_222_a_bqae * Gamma_bar_bqae
      // ###########################################################
      std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
      for (int indexp = 0; indexp < norbits; ++indexp)
      {
        auto p = allorb_vec[indexp];
        Orbit &op = Z.modelspace->GetOrbit(p);
        double jp = op.j2 / 2.;
        for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
        {
          if (q > p)
            continue;
          Orbit &oq = Z.modelspace->GetOrbit(q);
          double jq = oq.j2 / 2.;
          double zpq = 0;

          // loop abcde
          for (auto &c : Z.modelspace->all_orbits)
          {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double jc = oc.j2 / 2.;

            int J0min = std::abs(oc.j2 - op.j2) / 2;
            int J0max = (oc.j2 + op.j2) / 2;
            for (int J0 = J0min; J0 <= J0max; J0++)
            {
              zpq += intermediateTB.GetTBME_J(J0, J0, c, p, c, q);
            }
          }
          Z.OneBody(p, q) += 0.25 * zpq / (op.j2 + 1.0);
          if (p != q)
            Z.OneBody(q, p) += 0.25 * hZ * zpq / (op.j2 + 1.0);
          //--------------------------------------------------
        } // for q
      } // for p

      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_chi2_pp_fill1b"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }


      // *********************************************************************************** //
      //                               Diagram II_a and II_c                                 //
      // *********************************************************************************** //

      // ###########################################################
      //  diagram II_a
      //
      //  Pandya transform
      //  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
      //                        { k l J'}
      size_t n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();

      // The two body operator
      //  Chi_222_a :
      //            eta |
      //           _____|
      //          /\    |
      //   |     (  )
      //   |_____ \/
      //   | eta
      //
      //  Chi_222_a = \sum_pq (nbar_e * nbar_d * n_f * n_c - nbar_f * nbar_c * n_e * n_d )

      // TODO: Modify this to handle Tz-changing operators.
      // To do this, the intermediate object ~ Eta * Eta * Gamma should be generalized to
      // allow different CC channels from the bra and the ket. The Eta * Eta part will conserve
      // the channel still. Then in the loop over p,q,e below, need to find the CC channels
      // for the bra and the ket in the look up. Shouldn't be too terrible.
//      std::deque<arma::mat> IntermediateTwobody(n_nonzero);
//      std::cout << " " << __func__ << " line " << __LINE__ << std::endl;

      std::deque<arma::mat> EtaEta(n_nonzero);
      std::map<std::array<size_t,2>,arma::mat> EtaEtaGamma;
      std::vector<size_t> ch_bra_list;
      std::vector<size_t> ch_ket_list;

      for (size_t ch_cc_bra = 0; ch_cc_bra < n_nonzero; ++ch_cc_bra)
      {
        TwoBodyChannel_CC &tbc_cc_bra = Z.modelspace->GetTwoBodyChannel_CC(ch_cc_bra);
        size_t nBras_cc = tbc_cc_bra.GetNumberKets();
        EtaEta[ch_cc_bra] = arma::zeros(nBras_cc, nBras_cc*2);

        for (size_t ch_cc_ket = 0; ch_cc_ket < n_nonzero; ++ch_cc_ket)
        {
          TwoBodyChannel_CC &tbc_cc_ket = Z.modelspace->GetTwoBodyChannel_CC(ch_cc_ket);
          int nKets_cc = tbc_cc_ket.GetNumberKets();
          
          if ( tbc_cc_ket.J != tbc_cc_bra.J ) continue;  // We're still working with rotational scalar operators.
          if ( (tbc_cc_ket.parity + tbc_cc_bra.parity + Gamma.GetParity() )%2 != 0) continue; // conserve parity
          if ( std::abs( tbc_cc_ket.Tz - tbc_cc_bra.Tz) > Gamma.GetTRank() ) continue;
          EtaEtaGamma[{ch_cc_bra,ch_cc_ket}] = arma::zeros( nBras_cc, nKets_cc );
          ch_bra_list.push_back( ch_cc_bra );
          ch_ket_list.push_back( ch_cc_ket );
        }
      }
//      std::cout << " " << __func__ << " line " << __LINE__ << std::endl;

/// Pandya transformation
#pragma omp parallel for schedule(dynamic, 1)
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        size_t nKets_cc = tbc_cc.GetNumberKets();
        int J_cc = tbc_cc.J;

        // We don't need to double the rows of the Eta_bar matrix, since these lead to rows in
        // the IntermediateTwobody matrix which can be obtained by symmetry from the top rows.
        // So this will save us some time, and when filling the resuliting one-body matrix we just
        // ensure we always access with the ordering that is stored. -SRS
        arma::mat Eta_bar      = arma::mat(nKets_cc,     nKets_cc * 2, arma::fill::zeros); // SRS ADDED
        arma::mat Eta_bar_nnnn = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
//        arma::mat Gamma_bar    = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);    // SRS ADDED
                                                                                           //      Gamma_bar[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
//        arma::mat Chi_222_a;  // We don't do anything with this?
        //      IntermediateTwobody[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        // transform operator
        // loop over cross-coupled ph bras <ab| in this channel
        //      for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc)
        for (size_t ibra_cc = 0; ibra_cc < nKets_cc; ++ibra_cc)
        {
          int a, b;
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
          a = bra_cc.p;
          b = bra_cc.q;
          // if (ibra_cc >= nKets_cc and a == b)
          //   continue;

          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          double nbar_a = 1 - n_a;
          double ja = oa.j2 * 0.5;

          Orbit &ob = Z.modelspace->GetOrbit(b);
          double n_b = ob.occ;
          double nbar_b = 1 - n_b;
          double jb = ob.j2 * 0.5;

          // loop over cross-coupled kets |cd> in this channel
          for (int iket_cc = ibra_cc; iket_cc < nKets_cc * 2; ++iket_cc)
          {
            if ((iket_cc % nKets_cc) < ibra_cc)
              continue; // We'll get these from symmetry
            int c, d;
            if (iket_cc < nKets_cc)
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
              c = ket_cc_cd.p;
              d = ket_cc_cd.q;
            }
            else
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
              d = ket_cc_cd.p;
              c = ket_cc_cd.q;
            }

            Orbit &oc = Z.modelspace->GetOrbit(c);
            double n_c = oc.occ;
            double nbar_c = 1 - n_c;
            double jc = oc.j2 * 0.5;

            Orbit &od = Z.modelspace->GetOrbit(d);
            double n_d = od.occ;
            double nbar_d = 1 - n_d;
            double jd = od.j2 * 0.5;

            double occ_factor = nbar_c * nbar_b * n_a * n_d - n_c * n_b * nbar_a * nbar_d;

            // Check the isospin projection. If this isn't conserved in the usual channel,
            // then all the xcbad and yadcb will be zero and we don't need to bother computing SixJs.
//            if (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != Gamma.GetTRank() and std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != Eta.GetTRank())
//              continue;
            if ( std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != Eta.GetTRank())
              continue;

            int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
            int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
            double Xbar = 0;
            double Ybar = 0;
            int dJ_std = 1;
            if ((a == d or b == c))
            {
              dJ_std = 2;
              jmin += jmin % 2;
            }
            for (int J_std = jmin; J_std <= jmax; J_std += dJ_std)
            {

              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                Xbar -= (2 * J_std + 1) * sixj1 * Eta.TwoBody.GetTBME_J(J_std, a, d, c, b);
//                Ybar -= (2 * J_std + 1) * sixj1 * Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
              }
            }
            double flip_phase = Z.modelspace->phase((oa.j2 + ob.j2 + oc.j2 + od.j2) / 2);

            if (iket_cc < nKets_cc or (iket_cc >= nKets_cc and c != d))
            {
              // direct term
//              Gamma_bar(ibra_cc, iket_cc) = Ybar;
              Eta_bar(ibra_cc, iket_cc) = Xbar;
              Eta_bar_nnnn(ibra_cc, iket_cc) = Xbar * occ_factor;

              if (iket_cc != ibra_cc)
              {
                // Hermiticity: Xbar_cdab = hX * Xbar_abcd.  We get a minus sign on the occupation factor
//                Gamma_bar(iket_cc, ibra_cc) = hGamma * Ybar;
                Eta_bar_nnnn(iket_cc, ibra_cc) = hEta * Xbar * (-occ_factor);
                if (iket_cc < nKets_cc)
                {
                  Eta_bar(iket_cc, ibra_cc) = hEta * Xbar;
                }
              }
            }

            if (a != b)
            {
              // By exchange symmetry Xbar_badc = phase * hX * Xbar_abcd.  For Eta_bar_nnnn, this also requires swapping labels in the occupations -> minus sign.
//              Gamma_bar(ibra_cc + nKets_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Ybar * flip_phase * hGamma;
              Eta_bar_nnnn(ibra_cc + nKets_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Xbar * flip_phase * hEta * (-occ_factor);
              // Eta_bar( ibra_cc+nKets_cc, (iket_cc + nKets_cc)%(2*nKets_cc)) = Xbar * flip_phase * hEta;
            }

            if (iket_cc >= nKets_cc or (iket_cc < nKets_cc and c != d))
            {
              // Combined exchange symmetry and hermiticity
              // Xbar_dcba = phase * Xbar_abcd
//              Gamma_bar((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Ybar * flip_phase;
              Eta_bar_nnnn((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Xbar * flip_phase * (occ_factor);
              if (iket_cc >= nKets_cc)
                Eta_bar((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Xbar * flip_phase;
            }

          } // for iket-cc

          //-------------------
        } // for ibra_cc

        EtaEta[ch_cc] = (2 * J_cc + 1) * Eta_bar * Eta_bar_nnnn ;
//        std::cout << " line " << __LINE__ << "  ch_cc = " << ch_cc
//                  << " -> " << tbc_cc.J << " " << tbc_cc.Tz << " " << tbc_cc.parity << " nkets = " << tbc_cc.GetNumberKets()
//                  << "  dimension " << EtaEta[ch_cc].n_rows << " x " << EtaEta[ch_cc].n_cols << std::endl;
       
//        IntermediateTwobody[ch_cc] = EtaEta[ch_cc] * Gamma_bar;
//        EtaEtaGamma[{ch_cc,ch_cc}] = EtaEta[ch_cc] * Gamma_bar;
//        IntermediateTwobody[ch_cc] = (2 * J_cc + 1) * Eta_bar * Eta_bar_nnnn * Gamma_bar;
      } // for ch_cc


//      std::cout << " " << __func__ << " line " << __LINE__ << std::endl;
      if (Commutator::verbose)
      {
        Z.profiler.timer["_231_F_chi2_ph_fill_chi"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }


      int nch_braket = ch_bra_list.size();
      for (int ich=0; ich<nch_braket; ich++)
      {
         size_t ch_bra_cc = ch_bra_list[ich];
         size_t ch_ket_cc = ch_ket_list[ich];
         TwoBodyChannel_CC &tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
         TwoBodyChannel_CC &tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
         size_t nBras_cc = tbc_bra_cc.GetNumberKets();
         size_t nKets_cc = tbc_ket_cc.GetNumberKets();
         int J_cc = tbc_bra_cc.J;

         arma::mat Gamma_bar    = arma::mat(nBras_cc * 2, nKets_cc * 2, arma::fill::zeros);    // SRS ADDED
//         std::cout << " begin ibra,iket loops. ch_bra,ch_ket = " << ch_bra_cc << " ,  " << ch_ket_cc
//                   << " Nbras,Nkets = " << nBras_cc << " , " << nKets_cc
//                   << std::endl;

         for (size_t ibra_cc = 0; ibra_cc < nBras_cc; ++ibra_cc)
         {
           int a, b;
           Ket &bra_cc = tbc_bra_cc.GetKet(ibra_cc);
           a = bra_cc.p;
           b = bra_cc.q;
           // if (ibra_cc >= nKets_cc and a == b)
           //   continue;

           Orbit &oa = Z.modelspace->GetOrbit(a);
           double n_a = oa.occ;
           double nbar_a = 1 - n_a;
           double ja = oa.j2 * 0.5;

           Orbit &ob = Z.modelspace->GetOrbit(b);
           double n_b = ob.occ;
           double nbar_b = 1 - n_b;
           double jb = ob.j2 * 0.5;

           // loop over cross-coupled kets |cd> in this channel
           for (int iket_cc = ibra_cc; iket_cc < nKets_cc * 2; ++iket_cc)
           {
             if ((iket_cc % nKets_cc) < ibra_cc)
               continue; // We'll get these from symmetry
             int c, d;
             if (iket_cc < nKets_cc)
             {
               Ket &ket_cc_cd = tbc_ket_cc.GetKet(iket_cc);
               c = ket_cc_cd.p;
               d = ket_cc_cd.q;
             }
             else
             {
               Ket &ket_cc_cd = tbc_ket_cc.GetKet(iket_cc - nKets_cc);
               d = ket_cc_cd.p;
               c = ket_cc_cd.q;
             }

             Orbit &oc = Z.modelspace->GetOrbit(c);
             double n_c = oc.occ;
             double nbar_c = 1 - n_c;
             double jc = oc.j2 * 0.5;

             Orbit &od = Z.modelspace->GetOrbit(d);
             double n_d = od.occ;
             double nbar_d = 1 - n_d;
             double jd = od.j2 * 0.5;


             int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
              int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
              double Xbar = 0;
              double Ybar = 0;
              int dJ_std = 1;
              if ((a == d or b == c))
              {
                dJ_std = 2;
                jmin += jmin % 2;
              }
              for (int J_std = jmin; J_std <= jmax; J_std += dJ_std)
              {
  
                double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
                if (std::abs(sixj1) > 1e-8)
                {
                  Ybar -= (2 * J_std + 1) * sixj1 * Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
                }
              }
              double flip_phase = Z.modelspace->phase((oa.j2 + ob.j2 + oc.j2 + od.j2) / 2);

              if (iket_cc < nKets_cc or (iket_cc >= nKets_cc and c != d))
              {
                // direct term
                Gamma_bar(ibra_cc, iket_cc) = Ybar;
  
                if (ch_bra_cc==ch_ket_cc and iket_cc != ibra_cc)
                {
                  // Hermiticity: Xbar_cdab = hX * Xbar_abcd.  We get a minus sign on the occupation factor
                  Gamma_bar(iket_cc, ibra_cc) = hGamma * Ybar;
                }
              }
              if (a != b)
              {
                // By exchange symmetry Xbar_badc = phase * hX * Xbar_abcd.  For Eta_bar_nnnn, this also requires swapping labels in the occupations -> minus sign.
                Gamma_bar(ibra_cc + nBras_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Ybar * flip_phase * hGamma;
                // Eta_bar( ibra_cc+nKets_cc, (iket_cc + nKets_cc)%(2*nKets_cc)) = Xbar * flip_phase * hEta;
              }
  
              if (ch_bra_cc==ch_ket_cc  and ( iket_cc >= nKets_cc or (iket_cc < nKets_cc and c != d)) )
              {
                // Combined exchange symmetry and hermiticity
                // Xbar_dcba = phase * Xbar_abcd
                Gamma_bar((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nBras_cc) = Ybar * flip_phase;
              }
           }// for iket_cc
         }// for ibra_cc

//        std::cout << "  about to do mat mult. ch_bra, ch_ket " << ch_bra_cc << " " << ch_ket_cc
//                  << " : " << tbc_bra_cc.Tz << " " << tbc_bra_cc.parity << " , " << tbc_ket_cc.Tz << " " << tbc_ket_cc.parity
//                  << "  dimensions " << EtaEta[ch_bra_cc].n_rows << " x " << EtaEta[ch_bra_cc].n_cols
//                  << "  " << Gamma_bar.n_rows << " x " << Gamma_bar.n_cols
//                  << std::endl;

        EtaEtaGamma[{ch_bra_cc,ch_ket_cc}] = EtaEta[ch_bra_cc] * Gamma_bar;
      }

//      std::cout << " " << __func__ << " line " << __LINE__ << std::endl;


      // ###########################################################
      // diagram II_a
      //
      //  IIa_pq = 1/ (2 jp + 1) \sum_abeJ3 Chi_222_a_peab * Gamma_bar_abqe
      //
      // diagram II_c
      //
      //  IIc_pq = - 1/ (2 jp + 1) \sum_abe J3 Chi_222_a_eqab * Gamma_bar_abep
      // ###########################################################

// may be worth rolling p and q loops together for better load balancing
#pragma omp parallel for schedule(dynamic, 1)
      for (int indexd = 0; indexd < norbits; ++indexd)
      {
        auto p = allorb_vec[indexd];
        Orbit &op = Z.modelspace->GetOrbit(p);
        double jp = op.j2 / 2.;
        double j2hat2 = (op.j2 + 1.0);
        for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
        {
          if (q > p)
            continue;
          Orbit &oq = Z.modelspace->GetOrbit(q);
          double zij = 0;
          for (auto &e : Z.modelspace->all_orbits) // delta_jp jq
          {
            Orbit &oe = Z.modelspace->GetOrbit(e);

            int Jtmin = std::abs(op.j2 - oe.j2) / 2;
            int Jtmax = (op.j2 + oe.j2) / 2;
            int parity_pe = (op.l + oe.l) % 2;
            int parity_qe = (oq.l + oe.l) % 2;
            int Tz_pe = std::abs(op.tz2 - oe.tz2) / 2;
            int Tz_qe = std::abs(oq.tz2 - oe.tz2) / 2;
            double zij = 0;
            for (int Jt = Jtmin; Jt <= Jtmax; Jt++)
            {

              size_t ch_cc_pe = Z.modelspace->GetTwoBodyChannelIndex(Jt, parity_pe, Tz_pe);
              size_t ch_cc_qe = Z.modelspace->GetTwoBodyChannelIndex(Jt, parity_qe, Tz_qe);
              TwoBodyChannel_CC &tbc_cc_pe = Z.modelspace->GetTwoBodyChannel_CC(ch_cc_pe);
              TwoBodyChannel_CC &tbc_cc_qe = Z.modelspace->GetTwoBodyChannel_CC(ch_cc_qe);
              // Make sure we access the element <ab|X|cd> with a<=b. If we want the other ordering, we get a minus sign.
              int ind_pe = tbc_cc_pe.GetLocalIndex(p, e);
              int ind_ep = tbc_cc_pe.GetLocalIndex(e, p);
              int ind_qe = tbc_cc_qe.GetLocalIndex(q, e);
              int ind_eq = tbc_cc_qe.GetLocalIndex(e, q);

              if (p <= e)
              {
//                zij += IntermediateTwobody[ch_cc_pe,ch_cc_qe](ind_pe, ind_qe);
                zij += EtaEtaGamma[{ch_cc_pe,ch_cc_qe}](ind_pe, ind_qe);
              }
              else
              {
                zij -= EtaEtaGamma[{ch_cc_pe,ch_cc_qe}](ind_ep, ind_eq);
              }

              if (e <= q)
              {
                zij -= EtaEtaGamma[{ch_cc_qe,ch_cc_pe}](ind_eq, ind_ep);
              }
              else
              {
                zij += EtaEtaGamma[{ch_cc_qe,ch_cc_pe}](ind_qe, ind_pe);
              }
            }
            Z.OneBody(p, q) += zij / j2hat2;
            if (p != q)
              Z.OneBody(q, p) += hZ * zij / j2hat2;
          }
        }
      }

      Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
      return;

    } // comm223_231_chi2b





    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////

    void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z)
    {

      // Do some extra work to check if the operators being passed in are reduced,
      // because the comm223_232 routines assume not-reduced matrix elements
      const Operator * Etanred   = &Eta;  // Pointer to the non-reduced version of the operator
      const Operator * Gammanred = &Gamma;
      Operator Etatmp,Gammatmp;  // Declared, but not yet allocated, for reasons of scope.
      if ( Eta.IsReduced() ) // comm223_231 doesn't expect reduced operators. Need to make it not reduced.
      {
         Etatmp = Eta; // Actual copy, not a reference
         Etatmp.MakeNotReduced();
         Etanred = &Etatmp; // Now Etanred points to the not-reduced copy Etatmp
      }
      if ( Gamma.IsReduced() )
      {
         Gammatmp = Gamma;
         Gammatmp.MakeNotReduced();
         Gammanred = &Gammatmp;  // Now Gammanred points to the not-reduced copy Gammatmp
      }

      bool z_was_reduced = Z.IsReduced();
      if (z_was_reduced) Z.MakeNotReduced();



      if (use_1b_intermediates)
      {
        comm223_232_chi1b(*Etanred, *Gammanred, Z); // topology with 1-body intermediate (fast)
      }
      if (use_2b_intermediates)
      {
        comm223_232_chi2b(*Etanred, *Gammanred, Z); // topology with 2-body intermediate (slow)
      }

      if (z_was_reduced) Z.MakeReduced(); // put it back the way we found it...

      return;
    }


    ////////////////////////////////////////////////////////////////////////////
    /// factorized 223_232 double commutator with 1b intermediate
    ////////////////////////////////////////////////////////////////////////////
    void comm223_232_chi1b(const Operator &Eta, const Operator &Gamma, Operator &Z)
    {
      // global variables
      double t_start = omp_get_wtime();
      double t_internal = omp_get_wtime();
      Z.modelspace->PreCalculateSixJ();
      int norbits = Z.modelspace->all_orbits.size();
      // Two Body channels
      std::vector<size_t> ch_bra_list, ch_ket_list;
      for (auto &iter : Z.TwoBody.MatEl)
      {
        ch_bra_list.push_back(iter.first[0]);
        ch_ket_list.push_back(iter.first[1]);
      }
      size_t nch = ch_bra_list.size();
      // int nch = Z.modelspace->GetNumberTwoBodyChannels(); // number of TB channels
      // int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC(); // number of CC channels
      auto &Z2 = Z.TwoBody;

      bool Z_is_scalar = true;
      if (Z.TwoBody.rank_T != 0)
      {
        Z_is_scalar = false;
      }
      // determine symmetry
      int hEta = Eta.IsHermitian() ? 1 : -1;
      int hGamma = Gamma.IsHermitian() ? 1 : -1;
      int hZ = hGamma;
      // ####################################################################################
      //                      Factorization of Ia, Ib, IVa and IVb
      // ####################################################################################

      arma::mat CHI_I = Gamma.OneBody * 0;
      arma::mat CHI_II = Gamma.OneBody * 0;

      // The intermidate one body operator
      //  CHI_I :                            //  CHI_II :
      //          eta | p                    //          eta | p
      //         _____|                      //         _____|
      //       /\     |                      //       /\     |
      //   a  (  ) b  | c                    //   a  (  ) b  | c
      //       \/_____|                      //       \/~~~~~|
      //          eta |                      //        gamma |
      //              | q                    //              | q
      //-------------------------------------------------------------------------------------
      // CHI_I_pq  = 1/2 \sum_abcJ2 \hat(J_2) ( \bar{n}_a \bar{n}_c n_b - \bar{n}_c n_a n_c )
      //             eta^J2_bpac eta^J2_acbq
      //
      // CHI_II_pq = 1/2 \sum_abcJ2 \hat(J_2) ( \bar{n}_b \bar{n}_c n_a - \bar{n}_a n_b n_c )
      //             eta^J2_bcaq gamma^J2_apbc
      //-------------------------------------------------------------------------------------

      // Rolling the two loops into one helps with load balancing
      std::vector<index_t> p_list, q_list;
      for (auto p : Z.modelspace->all_orbits)
      {
        Orbit &op = Z.modelspace->GetOrbit(p);
        for (auto q : Eta.OneBodyChannels.at({op.l, op.j2, op.tz2}))
        {
          p_list.push_back(p);
          q_list.push_back(q);
        }
      }
      size_t ipq_max = p_list.size();
/// Build the intermediate one-body operators
#pragma omp parallel for schedule(dynamic)
      for (size_t ipq = 0; ipq < ipq_max; ipq++)
      {
        index_t p = p_list[ipq];
        index_t q = q_list[ipq];
        //      for (size_t p = 0; p < norbits; p++)
        //      {
        Orbit &op = Z.modelspace->GetOrbit(p);
        //        for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
        //        {
        Orbit &oq = Z.modelspace->GetOrbit(q);

        double chi_pq = 0;
        double chiY_pq = 0;

        for (auto a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;
          if (nbar_a < 1e-6)
            continue;

          for (auto i : Z.modelspace->holes)
          {
            Orbit &oi = Z.modelspace->GetOrbit(i);
            double n_i = oi.occ;

            for (auto j : Z.modelspace->holes)
            {
              Orbit &oj = Z.modelspace->GetOrbit(j);
              double n_j = oj.occ;

              double occfactor = nbar_a * n_i * n_j;
              if (occfactor < 1.e-7)
              {
                continue;
              }

              int J2min = std::max(std::abs(oa.j2 - oq.j2), std::abs(oi.j2 - oj.j2)) / 2;
              int J2max = std::min(oa.j2 + oq.j2, oi.j2 + oj.j2) / 2;

              for (int J2 = J2min; J2 <= J2max; J2++)
              {
                double xijaq = Eta.TwoBody.GetTBME_J(J2, J2, i, j, a, q);
                double xapij, yapij;
                if (Z_is_scalar)
                {
                  Eta.TwoBody.GetTBME_J_twoOps(Gamma.TwoBody, J2, J2, a, p, i, j, xapij, yapij);
                }
                else
                {
                  xapij = Eta.TwoBody.GetTBME_J(J2, J2, a, p, i, j);
                  yapij = Gamma.TwoBody.GetTBME_J(J2, J2, a, p, i, j);
                }
                chi_pq += 0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * xapij * xijaq;
                chiY_pq += 0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * yapij * xijaq;
              }
            } // for j

            for (auto b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double occfactor = nbar_a * nbar_b * n_i;

              if (std::abs(occfactor) < 1e-7)
                continue;

              int J2min = std::max({std::abs(oa.j2 - ob.j2), std::abs(oi.j2 - oq.j2), std::abs(oi.j2 - op.j2)}) / 2;
              int J2max = std::min({oa.j2 + ob.j2, oi.j2 + oq.j2, oi.j2 + op.j2}) / 2;

              for (int J2 = J2min; J2 <= J2max; J2++)
              {
                double xabiq = Eta.TwoBody.GetTBME_J(J2, J2, a, b, i, q);
                double xipab, yipab;

                if (Z_is_scalar)
                {
                  Eta.TwoBody.GetTBME_J_twoOps(Gamma.TwoBody, J2, J2, i, p, a, b, xipab, yipab);
                }
                else
                {
                  xipab = Eta.TwoBody.GetTBME_J(J2, J2, i, p, a, b);
                  yipab = Gamma.TwoBody.GetTBME_J(J2, J2, i, p, a, b);
                }

                chi_pq += 0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * xipab * xabiq;
                chiY_pq += 0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * yipab * xabiq;
              }
            } // for b

          } // for i
        } // for a
        CHI_I(p, q) = chi_pq;
        CHI_II(p, q) = chiY_pq;
        //        } // for q
      } // for p

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

/// Now use the intermediate to form the double commutator
#pragma omp parallel for schedule(dynamic, 1)
      for (int ich = 0; ich < nch; ich++)
      {
        size_t ch_bra = ch_bra_list[ich];
        size_t ch_ket = ch_ket_list[ich];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int J = tbc_bra.J;
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();
        for (size_t ibra = 0; ibra < nbras; ibra++)
        {
          Ket &bra = tbc_bra.GetKet(ibra);
          index_t p = bra.p;
          index_t q = bra.q;
          Orbit &op = Z.modelspace->GetOrbit(p);
          Orbit &oq = Z.modelspace->GetOrbit(q);
          int phasepq = bra.Phase(J);

          int ketmin = 0;
          if (ch_bra == ch_ket)
            ketmin = ibra;
          for (size_t iket = ketmin; iket < nkets; iket++)
          {
            Ket &ket = tbc_ket.GetKet(iket);
            index_t r = ket.p;
            index_t s = ket.q;
            Orbit &oR = Z.modelspace->GetOrbit(r);
            Orbit &os = Z.modelspace->GetOrbit(s);
            double zpqrs = 0;
            double gamma2b, eta2b;

            // SRS Modified so that we don't have to re-look up the channel and bra/ket indices.
            // This requires a bit more work to get the normalization and phases right
            // but at emax=6 is speeds this up by a factor of 2.
            for (auto b : Eta.OneBodyChannels.at({op.l, op.j2, op.tz2}))
            {
              auto ibra_bq = tbc_bra.GetLocalIndex(std::min(b, q), std::max(b, q));
              if (ibra_bq < 0 or ibra_bq > nbras)
                continue;
              double norm = (b == q ? PhysConst::SQRT2 : 1) * (p == q ? 1 / PhysConst::SQRT2 : 1);
              if (b > q)
                norm *= bra.Phase(tbc_bra.J);
              zpqrs += norm * CHI_I(p, b) * Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_bq, iket);
              if (Z_is_scalar)
                zpqrs += norm * hZ * CHI_II(b, p) * Eta.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_bq, iket);
              // zpqrs += CHI_I(p, b) * Gamma.TwoBody.GetTBME_J(J, J, b, q, r, s);
              // zpqrs += hZ * CHI_II(b, p) * Eta.TwoBody.GetTBME_J(J, J, b, q, r, s);
            }
            for (auto b : Eta.OneBodyChannels.at({oq.l, oq.j2, oq.tz2}))
            {
              auto ibra_pb = tbc_bra.GetLocalIndex(std::min(p, b), std::max(p, b));
              if (ibra_pb < 0 or ibra_pb > nbras)
                continue;
              double norm = (b == p ? PhysConst::SQRT2 : 1) * (p == q ? 1 / PhysConst::SQRT2 : 1);
              if (p > b)
                norm *= bra.Phase(tbc_bra.J);
              zpqrs += norm * CHI_I(q, b) * Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_pb, iket);
              if (Z_is_scalar)
                zpqrs += norm * hZ * CHI_II(b, q) * Eta.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra_pb, iket);
              // zpqrs += CHI_I(q, b) *     Gamma.TwoBody.GetTBME_J(J, J, p, b, r, s);
              // zpqrs += hZ * CHI_II(b, q) * Eta.TwoBody.GetTBME_J(J, J, p, b, r, s);
            }
            for (auto b : Eta.OneBodyChannels.at({oR.l, oR.j2, oR.tz2}))
            {
              auto iket_bs = tbc_ket.GetLocalIndex(std::min(b, s), std::max(b, s));
              if (iket_bs < 0 or iket_bs > nkets)
                continue;
              double norm = (b == s ? PhysConst::SQRT2 : 1) * (r == s ? 1 / PhysConst::SQRT2 : 1);
              if (b > s)
                norm *= ket.Phase(tbc_ket.J);
              zpqrs += norm * Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_bs) * CHI_I(b, r);
              if (Z_is_scalar)
                zpqrs -= norm * Eta.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_bs) * CHI_II(b, r);
              // zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_I(b, r);
              // zpqrs -=   Eta.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_II(b, r);
            }
            for (auto b : Eta.OneBodyChannels.at({os.l, os.j2, os.tz2}))
            {
              auto iket_rb = tbc_ket.GetLocalIndex(std::min(r, b), std::max(r, b));
              if (iket_rb < 0 or iket_rb > nkets)
                continue;
              double norm = (b == r ? PhysConst::SQRT2 : 1) * (r == s ? 1 / PhysConst::SQRT2 : 1);
              if (r > b)
                norm *= ket.Phase(tbc_ket.J);
              zpqrs += norm * Gamma.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_rb) * CHI_I(b, s);
              if (Z_is_scalar)
                zpqrs -= norm * Eta.TwoBody.GetTBME_norm(ch_bra, ch_ket, ibra, iket_rb) * CHI_II(b, s);
              // zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_I(b, s);
              // zpqrs -=   Eta.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_II(b, s);
            }

            //            if (p == q)
            //              zpqrs /= PhysConst::SQRT2;
            //            if (r == s)
            //              zpqrs /= PhysConst::SQRT2;
            Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpqrs);
          } // for iket
        } // for ibra

      } // for itmat

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      CHI_I.clear();
      CHI_II.clear();

      // Timer
      Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
      return;
    } // comm223_232_chi1b

    ////////////////////////////////////////////////////////////////////////////
    /// factorized 223_232 double commutator with 2b intermediate
    ////////////////////////////////////////////////////////////////////////////
    void comm223_232_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z)
    {
      // global variables
      double t_start = omp_get_wtime();
      double t_internal = omp_get_wtime();
      Z.modelspace->PreCalculateSixJ();
      int norbits = Z.modelspace->all_orbits.size();
      // Two Body channels
      std::vector<size_t> ch_bra_list, ch_ket_list;
      for (auto &iter : Z.TwoBody.MatEl)
      {
        ch_bra_list.push_back(iter.first[0]);
        ch_ket_list.push_back(iter.first[1]);
      }
      size_t nch = ch_bra_list.size();
      int nch_eta = Eta.modelspace->GetNumberTwoBodyChannels();
      // int nch = Z.modelspace->GetNumberTwoBodyChannels(); // number of TB channels
      int n_nonzero = Eta.modelspace->GetNumberTwoBodyChannels_CC(); // number of CC channels
      auto &Z2 = Z.TwoBody;

      bool Z_is_scalar = true;
      if (Z.TwoBody.rank_T != 0)
      {
        Z_is_scalar = false;
      }
      // determine symmetry
      int hEta = Eta.IsHermitian() ? 1 : -1;
      int hGamma = Gamma.IsHermitian() ? 1 : -1;
      // int hZ = Z.IsHermitian() ? 1 : -1;
      int hZ = hGamma;

      // *********************************************************************************** //
      //                             Diagram II and III                                      //
      // *********************************************************************************** //

      //______________________________________________________________________
      // global array
      std::deque<arma::mat> bar_Eta(n_nonzero);   // released
      std::deque<arma::mat> bar_Gamma(n_nonzero); // released
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        bar_Eta[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        bar_Gamma[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      std::deque<arma::mat> barCHI_III(n_nonzero);    //  released
      std::deque<arma::mat> bar_CHI_V(n_nonzero);     // released
      std::deque<arma::mat> bar_CHI_VI(n_nonzero);    //  released
      std::deque<arma::mat> bar_CHI_VI_II(n_nonzero); //  released

      /// Pandya transformation
      /// construct bar_Gamma, bar_Eta, nnnbar_Eta, nnnbar_Eta_d, bar_CHI_VI
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        int J_cc = tbc_cc.J;
        if (nKets_cc < 1)
        {
          continue;
        }

        arma::mat nnnbar_Eta = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        arma::mat nnnbar_Eta_d = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        for (int ibra_cc = 0; ibra_cc < nKets_cc; ++ibra_cc)
        {
          int a, b;
          if (ibra_cc < nKets_cc)
          {
            Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
            a = bra_cc.p;
            b = bra_cc.q;
          }
          else
          {
            Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nKets_cc);
            b = bra_cc.p;
            a = bra_cc.q;
          }
          if (ibra_cc >= nKets_cc and a == b)
            continue;

          Orbit &oa = Z.modelspace->GetOrbit(a);
          double ja = oa.j2 * 0.5;
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          Orbit &ob = Z.modelspace->GetOrbit(b);
          double jb = ob.j2 * 0.5;
          double n_b = ob.occ;
          double nbar_b = 1.0 - n_b;

          // loop over cross-coupled kets |cd> in this channel
          for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc)
          {
            if ((iket_cc % nKets_cc) < ibra_cc)
              continue; // We'll get these from symmetry
            int c, d;
            if (iket_cc < nKets_cc)
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
              c = ket_cc_cd.p;
              d = ket_cc_cd.q;
            }
            else
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
              d = ket_cc_cd.p;
              c = ket_cc_cd.q;
            }
            // if (iket_cc >= nKets_cc and c == d)
            //  continue;

            Orbit &oc = Z.modelspace->GetOrbit(c);
            double jc = oc.j2 * 0.5;
            double n_c = oc.occ;
            double nbar_c = 1.0 - n_c;

            Orbit &od = Z.modelspace->GetOrbit(d);
            double jd = od.j2 * 0.5;
            double n_d = od.occ;
            double nbar_d = 1.0 - n_d;

            double occ_AbarBC = (nbar_a * n_b * n_c + n_a * nbar_b * nbar_c);
            double occ_ABbarD = (n_a * nbar_b * n_d + nbar_a * n_b * nbar_d);

            double occ_BCDbar = (n_b * n_c * nbar_d + nbar_b * nbar_c * n_d);
            double occ_ACbarD = (n_a * nbar_c * n_d + nbar_a * n_c * nbar_d);

            int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
            int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
            double Etabar = 0;
            double Gammabar = 0;
            int dJ_std = 1;
            if ((a == d or b == c))
            {
              dJ_std = 2;
              jmin += jmin % 2;
            }
            for (int J_std = jmin; J_std <= jmax; J_std += dJ_std)
            {
              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                Etabar -= (2 * J_std + 1) * sixj1 * Eta.TwoBody.GetTBME_J(J_std, a, d, c, b);
                Gammabar -= (2 * J_std + 1) * sixj1 * Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
              }
            }

            double flip_phase = Z.modelspace->phase((oa.j2 + ob.j2 + oc.j2 + od.j2) / 2);
            if (iket_cc < nKets_cc or (iket_cc >= nKets_cc and c != d))
            {
              // direct term
              bar_Gamma[ch_cc](ibra_cc, iket_cc) = Gammabar;
              bar_Eta[ch_cc](ibra_cc, iket_cc) = Etabar;
              nnnbar_Eta(ibra_cc, iket_cc) = Etabar * occ_AbarBC;
              nnnbar_Eta_d(ibra_cc, iket_cc) = Etabar * occ_ABbarD;

              if (iket_cc != ibra_cc)
              {
                // Hermiticity: Xbar_cdab = hX * Xbar_abcd.
                bar_Gamma[ch_cc](iket_cc, ibra_cc) = hGamma * Gammabar;
                bar_Eta[ch_cc](iket_cc, ibra_cc) = hEta * Etabar;
                nnnbar_Eta(iket_cc, ibra_cc) = hEta * Etabar * occ_ACbarD;
                nnnbar_Eta_d(iket_cc, ibra_cc) = hEta * Etabar * occ_BCDbar;
              }
            }

            if (a != b)
            {
              // By exchange symmetry Xbar_badc = phase * hX * Xbar_abcd.
              bar_Gamma[ch_cc](ibra_cc + nKets_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Gammabar * flip_phase * hGamma;
              bar_Eta[ch_cc](ibra_cc + nKets_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Etabar * flip_phase * hEta;
              nnnbar_Eta(ibra_cc + nKets_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Etabar * flip_phase * hEta * occ_ABbarD;
              nnnbar_Eta_d(ibra_cc + nKets_cc, (iket_cc + nKets_cc) % (2 * nKets_cc)) = Etabar * flip_phase * hEta * occ_AbarBC;
            }

            if (iket_cc >= nKets_cc or (iket_cc < nKets_cc and c != d))
            {
              // Combined exchange symmetry and hermiticity
              // Xbar_dcba = phase * Xbar_abcd
              bar_Gamma[ch_cc]((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Gammabar * flip_phase;
              bar_Eta[ch_cc]((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Etabar * flip_phase;
              nnnbar_Eta((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Etabar * flip_phase * occ_BCDbar;
              nnnbar_Eta_d((iket_cc + nKets_cc) % (2 * nKets_cc), ibra_cc + nKets_cc) = Etabar * flip_phase * occ_ACbarD;
            }
          }
          //-------------------
        }

        barCHI_III[ch_cc] = bar_Eta[ch_cc] * nnnbar_Eta;
        bar_CHI_V[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta;
        bar_CHI_VI[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta_d;
        bar_CHI_VI_II[ch_cc] = hEta * (nnnbar_Eta_d).t() * bar_Gamma[ch_cc];
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      //-------------------------------------------------------------------------------
      // intermediate operator for diagram IIa  and IIc
      // Theintermediate two body operator
      //  Chi_III :
      //            eta |
      //           _____|
      //          /\    |
      //   |     (  )
      //   |_____ \/
      //   | eta
      TwoBodyME Chi_III_Op = Eta.TwoBody;
      Chi_III_Op.Erase();
      // Inverse Pandya transformation
      //  X^J_ijkl  = - ( 1- P_ij )  sum_J' (2J'+1)  { i j J }  \bar{X}^J'_il`kj`
      //                                             { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch_eta; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();
        for (int ibra = 0; ibra < nKets; ++ibra)
        {
          Ket &bra = tbc.GetKet(ibra);
          size_t i = bra.p;
          size_t j = bra.q;
          Orbit &oi = *(bra.op);
          Orbit &oj = *(bra.oq);
          int ji = oi.j2;
          int jj = oj.j2;

          for (int iket = 0; iket < nKets * 2; ++iket)
          {
            size_t k, l;
            if (iket < nKets)
            {
              Ket &ket = tbc.GetKet(iket);
              k = ket.p;
              l = ket.q;
            }
            else
            {
              Ket &ket = tbc.GetKet(iket - nKets);
              l = ket.p;
              k = ket.q;
            }

            Orbit &ok = Z.modelspace->GetOrbit(k);
            Orbit &ol = Z.modelspace->GetOrbit(l);
            int jk = ok.j2;
            int jl = ol.j2;
            double commij = 0;
            double commji = 0;

            // ijkl
            int parity_cc = (oi.l + ol.l) % 2;
            int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
            int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
            int Jpmax = std::min(ji + jl, jj + jk) / 2;

            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {

              double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              if (std::abs(sixj) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
              int indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k));
              if (indx_il < 0 or indx_kj < 0)
                continue;
              indx_il += (i > l ? nkets_cc : 0);
              indx_kj += (k > j ? nkets_cc : 0);

              double me1 = barCHI_III[ch_cc](indx_il, indx_kj);
              commij -= (2 * Jprime + 1) * sixj * me1;
            }

            // jikl, exchange i and j
            parity_cc = (oi.l + ok.l) % 2;
            Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
            Jpmin = std::max(std::abs(int(jj - jl)), std::abs(int(jk - ji))) / 2;
            Jpmax = std::min(int(jj + jl), int(jk + ji)) / 2;

            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);

              if (std::abs(sixj) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              int indx_ik = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
              int indx_lj = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));

              if (indx_ik < 0 or indx_lj < 0)
                continue;
              indx_ik += (k > i ? nkets_cc : 0);
              indx_lj += (j > l ? nkets_cc : 0);
              double me1 = barCHI_III[ch_cc](indx_lj, indx_ik);
              commji -= (2 * Jprime + 1) * sixj * me1;
            }

            double zijkl = (commij - Z.modelspace->phase((ji + jj) / 2 - J0) * commji);
            //          Chi_III[ch](ibra, iket) += zijkl;
            if (i == j)
              zijkl /= PhysConst::SQRT2;
            if (k == l)
              zijkl /= PhysConst::SQRT2;
            if (iket < nKets)
              Chi_III_Op.GetMatrix(ch, ch)(ibra, iket) += zijkl;
            if (iket >= nKets)
              Chi_III_Op.GetMatrix(ch, ch)(ibra, iket % nKets) -= zijkl * Z.modelspace->phase((jk + jl) / 2 - J0);
          }
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      //------------------------------------------------------------------------------
      // intermediate operator for diagram IIIc and IIId
      //------------------------------------------------------------------------------
      //  The intermediate two body operator
      //  Chi_VI :
      //            eta |
      //           _____|
      //          /\    |
      //         (  )     |
      //          \/~~~~~~|
      //            gamma |
      //
      //  Chi_VI_cdqh = \sum_{ab} (nbar_a * n_b * n_c - nbar_a * nbar_b * n_c )
      //                 ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                { J3 J4 J0 } { J3 J4 J0 }
      //                { jc jd ja } { jq jh jb }
      //
      //                \bar{Eta}_bq`ac`  Gamma_dahb
      //-------------------------------------------------------------------------------
      TwoBodyME Chi_VI_Op = Gamma.TwoBody;
      Chi_VI_Op.Erase();

      TwoBodyME Chi_VI_II_Op = Gamma.TwoBody;
      Chi_VI_II_Op.Erase();
      // BUILD CHI_VI
      // Inverse Pandya transformation
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        size_t ch_bra = ch_bra_list[ch];
        size_t ch_ket = ch_ket_list[ch];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();
        if (nbras == 0 or nkets == 0)
          continue;
        int J0 = tbc_bra.J;
        for (int ibra = 0; ibra < nbras * 2; ++ibra)
        {
          size_t i, j;
          if (ibra < nbras)
          {
            Ket &bra = tbc_bra.GetKet(ibra);
            i = bra.p;
            j = bra.q;
          }
          else
          {
            Ket &bra = tbc_bra.GetKet(ibra - nbras);
            i = bra.q;
            j = bra.p;
          }
          // if (ibra >= nbras and i == j)
          //   continue;

          Orbit &oi = Z.modelspace->GetOrbit(i);
          int ji = oi.j2;
          Orbit &oj = Z.modelspace->GetOrbit(j);
          int jj = oj.j2;

          for (int iket = 0; iket < nkets * 2; ++iket)
          {
            size_t k, l;
            if (iket < nkets)
            {
              Ket &ket = tbc_ket.GetKet(iket);
              k = ket.p;
              l = ket.q;
            }
            else
            {
              Ket &ket = tbc_ket.GetKet(iket - nkets);
              k = ket.q;
              l = ket.p;
            }
            // if (iket >= nkets and k == l)
            //   continue;

            Orbit &ok = Z.modelspace->GetOrbit(k);
            Orbit &ol = Z.modelspace->GetOrbit(l);
            int jk = ok.j2;
            int jl = ol.j2;
            double commijkl = 0;
            double commijlk = 0;

            double commijkld = 0;
            double commjikld = 0;

            // ijkl, direct term        -->  il kj
            int parity_cc = (oi.l + ol.l) % 2;
            int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
            int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
            int Jpmax = std::min(ji + jl, jj + jk) / 2;
            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              if (std::abs(sixj) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 1)
                continue;

              int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
              int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
              if (indx_il < 0 or indx_kj < 0)
                continue;

              indx_il += (i > l ? nkets_cc : 0);
              indx_kj += (k > j ? nkets_cc : 0);
              double me1 = bar_CHI_VI[ch_cc](indx_il, indx_kj);
              commijkl -= (2 * Jprime + 1) * sixj * me1;
              double me12 = bar_CHI_VI_II[ch_cc](indx_il, indx_kj);
              commijkld -= (2 * Jprime + 1) * sixj * me12;
            }

            // ijlk,  exchange k and l -->  ik lj
            parity_cc = (oi.l + ok.l) % 2;
            Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
            Jpmin = std::max(std::abs(ji - jk), std::abs(jj - jl)) / 2;
            Jpmax = std::min(ji + jk, jj + jl) / 2;
            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jl * 0.5, jk * 0.5, Jprime);
              if (std::abs(sixj) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              int indx_ik = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
              int indx_lj = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
              if (indx_ik < 0 or indx_lj < 0)
                continue;

              int indx_ki = indx_ik;
              int indx_jl = indx_lj;

              // exchange k and l
              indx_ik += (i > k ? nkets_cc : 0);
              indx_lj += (l > j ? nkets_cc : 0);
              double me2 = bar_CHI_VI[ch_cc](indx_ik, indx_lj);
              commijlk -= (2 * Jprime + 1) * sixj * me2;

              indx_ki += (k > i ? nkets_cc : 0);
              indx_jl += (j > l ? nkets_cc : 0);
              double me21 = bar_CHI_VI_II[ch_cc](indx_jl, indx_ki);
              commjikld -= (2 * Jprime + 1) * sixj * me21;
            }

            double zijkl = (commijkl - Z.modelspace->phase((jk + jl) / 2 - J0) * commijlk);
            double zijkl_II = (commijkld - Z.modelspace->phase((ji + jj) / 2 - J0) * commjikld);

            if (i == j)
            {
              zijkl /= PhysConst::SQRT2;
              zijkl_II /= PhysConst::SQRT2;
            }
            if (k == l)
            {
              zijkl /= PhysConst::SQRT2;
              zijkl_II /= PhysConst::SQRT2;
            }

            if (iket < nkets)
            {
              if (ibra < nbras)
              {
                Chi_VI_Op.GetMatrix(ch_bra, ch_ket)(ibra, iket) += zijkl;
              }
              if (ibra >= nbras)
              {
                Chi_VI_Op.GetMatrix(ch_bra, ch_ket)(ibra % nbras, iket) -= zijkl * Z.modelspace->phase((ji + jj) / 2 - J0);
              }
            }

            if (ibra < nbras)
            {
              if (iket < nkets)
              {
                Chi_VI_II_Op.GetMatrix(ch_bra, ch_ket)(ibra, iket) += zijkl_II;
              }
              if (iket >= nkets)
              {
                Chi_VI_II_Op.GetMatrix(ch_bra, ch_ket)(ibra, iket % nkets) -= zijkl_II * Z.modelspace->phase((jk + jl) / 2 - J0);
              }
            }
          }
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      // release memory
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        bar_CHI_VI[ch_cc].clear();
        bar_CHI_VI_II[ch_cc].clear();
      }
      bar_CHI_VI.clear();
      bar_CHI_VI_II.clear();

// Diagram IIa and Diagram IIc
// Diagram IIIc and Diagram IIId
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        size_t ch_bra = ch_bra_list[ch];
        size_t ch_ket = ch_ket_list[ch];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();

        if (nbras < 1 or nkets < 1)
          continue;
        // Diagram IIa and IIc
        arma::mat Multi_matirx = Chi_III_Op.GetMatrix(ch_bra, ch_bra) * Gamma.TwoBody.GetMatrix(ch_bra, ch_ket);
        Multi_matirx += hZ * Gamma.TwoBody.GetMatrix(ch_bra, ch_ket) * (Chi_III_Op.GetMatrix(ch_ket, ch_ket).t());
        // Diagram IIIc and Diagram IIId
        Multi_matirx += -Eta.TwoBody.GetMatrix(ch_bra) * Chi_VI_Op.GetMatrix(ch_bra, ch_ket) - (Chi_VI_II_Op.GetMatrix(ch_bra, ch_ket) * Eta.TwoBody.GetMatrix(ch_ket));
        Z2.GetMatrix(ch_bra, ch_ket) += Multi_matirx;
      } // J0 channel

      // release memory
      Chi_III_Op.Deallocate();
      Chi_VI_Op.Deallocate();
      Chi_VI_II_Op.Deallocate();

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      // Diagram IIb and Diagram IId
      std::deque<arma::mat> barCHI_III_RC(n_nonzero); // released Recoupled bar CHI_III
      std::deque<arma::mat> bar_CHI_V_RC(n_nonzero);  // released
      /// build intermediate bar operator
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        barCHI_III_RC[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        bar_CHI_V_RC[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      //------------------------------------------------------------------------------
      // intermediate operator for diagram IIb and IId    IIIa and IIIb
      //------------------------------------------------------------------------------
      /// Pandya transformation only recouple the angula momentum
      /// IIb and IId
      /// diagram IIIa - diagram IIIb
      //  \bar{X}^J_ijkl  = sum_J' (2J'+1)  { i j J }  (-)^(j+k+J')  \bar{X}^J'_il`jk`
      //                                    { k l J'}
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        int J_cc = tbc_cc.J;
        for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc)
        {
          int a, b;
          if (ibra_cc < nKets_cc)
          {
            Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
            a = bra_cc.p;
            b = bra_cc.q;
          }
          else
          {
            Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nKets_cc);
            b = bra_cc.p;
            a = bra_cc.q;
          }
          if (ibra_cc >= nKets_cc and a == b)
            continue;
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double ja = oa.j2 * 0.5;
          double jb = ob.j2 * 0.5;

          // loop over cross-coupled kets |cd> in this channel
          for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc)
          {
            int c, d;
            if (iket_cc < nKets_cc)
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
              c = ket_cc_cd.p;
              d = ket_cc_cd.q;
            }
            else
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
              d = ket_cc_cd.p;
              c = ket_cc_cd.q;
            }
            if (iket_cc >= nKets_cc and c == d)
              continue;
            Orbit &oc = Z.modelspace->GetOrbit(c);
            Orbit &od = Z.modelspace->GetOrbit(d);
            double jc = oc.j2 * 0.5;
            double jd = od.j2 * 0.5;

            int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
            int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
            double XbarIIbd = 0;
            double XbarIIIab = 0;

            for (int J_std = jmin; J_std <= jmax; J_std++)
            {
              // int phaseFactor = Z.modelspace->phase(J_std + (oc.j2 + ob.j2) / 2);
              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                int parity_cc = (oa.l + od.l) % 2;
                int Tz_cc = std::abs(oa.tz2 - od.tz2) / 2;
                int ch_cc_old = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_cc, Tz_cc);

                TwoBodyChannel_CC &tbc_cc_old = Z.modelspace->GetTwoBodyChannel_CC(ch_cc_old);
                int nkets = tbc_cc_old.GetNumberKets();
                int indx_ad = tbc_cc_old.GetLocalIndex(std::min(int(a), int(d)), std::max(int(a), int(d)));
                int indx_bc = tbc_cc_old.GetLocalIndex(std::min(int(b), int(c)), std::max(int(b), int(c)));
                if (indx_ad >= 0 and indx_bc >= 0)
                {
                  int indx_cb = indx_bc;
                  if (a > d)
                    indx_ad += nkets;
                  if (b > c)
                    indx_bc += nkets;
                  if (c > b)
                    indx_cb += nkets;
                  XbarIIbd -= Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J_std) * (2 * J_std + 1) * sixj1 * (barCHI_III[ch_cc_old](indx_bc, indx_ad) + barCHI_III[ch_cc_old](indx_ad, indx_bc));
                  XbarIIIab += Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J_std) * (2 * J_std + 1) * sixj1 * (bar_CHI_V[ch_cc_old](indx_ad, indx_bc) - hZ * bar_CHI_V[ch_cc_old](indx_bc, indx_ad));
                }
              }
            }
            barCHI_III_RC[ch_cc](ibra_cc, iket_cc) = XbarIIbd;
            bar_CHI_V_RC[ch_cc](ibra_cc, iket_cc) = XbarIIIab;
          }
          //-------------------
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      /// release memory
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        barCHI_III[ch_cc].clear();
        bar_CHI_V[ch_cc].clear();
      }
      barCHI_III.clear();
      bar_CHI_V.clear();

      // ##########################################
      //      diagram IIIa - diagram IIIb
      // ##########################################
      std::deque<arma::mat> CHI_V_final(n_nonzero);
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        CHI_V_final[ch_cc] = bar_Eta[ch_cc] * bar_CHI_V_RC[ch_cc];
      }
      /// release memory
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        bar_CHI_V_RC[ch_cc].clear();
      }
      bar_CHI_V_RC.clear();

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      //  Inverse Pandya transformation
      //  diagram IIIa - diagram IIIb
      //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
      //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
      //                                    { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        int ch_bra = ch_bra_list[ch];
        int ch_ket = ch_ket_list[ch];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();
        int J0 = tbc_bra.J;

        if (nbras == 0 or nkets == 0)
          continue;

        for (int ibra = 0; ibra < nbras; ++ibra)
        {
          Ket &bra = tbc_bra.GetKet(ibra);
          size_t i = bra.p;
          size_t j = bra.q;
          Orbit &oi = *(bra.op);
          Orbit &oj = *(bra.oq);
          int ji = oi.j2;
          int jj = oj.j2;
          int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

          int ketmin = 0;
          if (ch_bra == ch_ket)
            ketmin = ibra;
          for (int iket = ketmin; iket < nkets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc_ket.GetKet(iket);
            k = ket.p;
            l = ket.q;

            Orbit &ok = Z.modelspace->GetOrbit(k);
            Orbit &ol = Z.modelspace->GetOrbit(l);
            int jk = ok.j2;
            int jl = ol.j2;
            double commijkl = 0;
            double commjikl = 0;
            double commijlk = 0;
            double commjilk = 0;

            // jikl, direct term        -->  jl  ki
            // ijlk, exchange ij and kl -->  lj  ik
            int parity_cc = (oi.l + ok.l) % 2;
            int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
            int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
            int Jpmax = std::min(jj + jl, ji + jk) / 2;
            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj1 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              if (std::abs(sixj1) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);

              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 1)
                continue;

              int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
              int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
              if (indx_jl < 0 or indx_ik < 0)
                continue;

              int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
              // direct term
              indx_jl += (j > l ? nkets_cc : 0);
              indx_ik += (i > k ? nkets_cc : 0);
              double me1 = CHI_V_final[ch_cc](indx_jl, indx_ik);
              commjikl -= phase1 * (2 * Jprime + 1) * sixj1 * me1;

              int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);
              // exchange ij and kl
              double me2 = CHI_V_final[ch_cc](indx_ik, indx_jl);
              commijlk -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
            }

            // ijkl,  exchange i and j -->  il  kj
            // jilk,  exchange k and l -->  jk li
            parity_cc = (oi.l + ol.l) % 2;
            Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
            Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
            Jpmax = std::min(ji + jl, jj + jk) / 2;
            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj1 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              if (std::abs(sixj1) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 1)
                continue;

              int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
              int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
              if (indx_il < 0 or indx_jk < 0)
                continue;

              int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
              // exchange k and l
              indx_il += (i > l ? nkets_cc : 0);
              indx_jk += (j > k ? nkets_cc : 0);
              double me1 = CHI_V_final[ch_cc](indx_jk, indx_il);
              commjilk -= phase1 * (2 * Jprime + 1) * sixj1 * me1;

              int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);
              // exchange i and j
              double me2 = CHI_V_final[ch_cc](indx_il, indx_jk);
              commijkl -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
            }

            double zijkl = (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
            zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

            if (i == j)
              zijkl /= PhysConst::SQRT2;
            if (k == l)
              zijkl /= PhysConst::SQRT2;

            Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, phaseFactor * zijkl);
          }
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        CHI_V_final[ch_cc].clear();
      }
      CHI_V_final.clear();

      // ######################################
      // declare CHI_IV
      //------------------------------------------------------------------------------
      // intermediate operator for diagram IIe and IIf
      //------------------------------------------------------------------------------
      // The intermediate two body operator
      //  Chi_IV :
      //        q  | eta |  b
      //           |_____|
      //        a  |_____|  c
      //           | eta |
      //        p  |     |  d
      /////////////////////////////////////////////////////////////////////////////////
      std::deque<arma::mat> CHI_IV(nch_eta); // released
      for (int ch = 0; ch < nch_eta; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int nKets = tbc.GetNumberKets();
        // Not symmetric
        CHI_IV[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      }

      //------------------------------------------------------------------------------
      //                      Factorization of IIIe and IIIf
      //------------------------------------------------------------------------------
      //
      // The intermediate two body operator
      //  CHI_VII :
      //        g  |     |  c
      //           |~~~~~|
      //        a  |     |  b
      //           |_____|
      //        h  |     |  d
      //------------------------------------------------------------------------------
      std::deque<arma::mat> CHI_VII(nch); // released
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int nKets = tbc.GetNumberKets();
        // Not symmetric
        CHI_VII[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      // this loop appears to be broken.
      // full matrix
#pragma omp parallel for
      for (int ch = 0; ch < nch_eta; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();
        arma::mat Eta_matrix(2 * nKets, 2 * nKets);
        arma::mat Eta_matrix_c(2 * nKets, 2 * nKets);
        arma::mat Eta_matrix_d(2 * nKets, 2 * nKets);
        arma::mat Gamma_matrix(2 * nKets, 2 * nKets);

        for (int ibra = 0; ibra < nKets; ++ibra)
        {
          Ket &bra = tbc.GetKet(ibra);
          size_t i = bra.p;
          size_t j = bra.q;
          Orbit &oi = *(bra.op);
          Orbit &oj = *(bra.oq);
          int ji = oi.j2;
          int jj = oj.j2;
          double n_i = oi.occ;
          double bar_n_i = 1. - n_i;
          double n_j = oj.occ;
          double bar_n_j = 1. - n_j;

          for (int iket = 0; iket < nKets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc.GetKet(iket);
            k = ket.p;
            l = ket.q;
            Orbit &ok = Z.modelspace->GetOrbit(k);
            Orbit &ol = Z.modelspace->GetOrbit(l);
            int jk = ok.j2;
            int jl = ol.j2;
            double n_k = ok.occ;
            double bar_n_k = 1. - n_k;
            double n_l = ol.occ;
            double bar_n_l = 1. - n_l;

            double occfactor_k = (bar_n_i * bar_n_j * n_k + n_i * n_j * bar_n_k);
            double occfactor_l = (bar_n_i * bar_n_j * n_l + n_i * n_j * bar_n_l);

            double EtaME = Eta.TwoBody.GetTBME_J(J0, i, j, k, l);
            double GammaME = Gamma.TwoBody.GetTBME_J(J0, i, j, k, l);

            Eta_matrix(ibra, iket) = EtaME;
            Eta_matrix_c(ibra, iket) = occfactor_k * EtaME;
            Eta_matrix_d(ibra, iket) = occfactor_l * EtaME;
            if (Z_is_scalar)
              Gamma_matrix(ibra, iket) = GammaME;
            if (i != j)
            {
              int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
              Eta_matrix(ibra + nKets, iket) = phase * EtaME;
              Eta_matrix_c(ibra + nKets, iket) = occfactor_k * phase * EtaME;
              Eta_matrix_d(ibra + nKets, iket) = occfactor_l * phase * EtaME;
              if (Z_is_scalar)
                Gamma_matrix(ibra + nKets, iket) = phase * GammaME;
              if (k != l)
              {
                phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
                Eta_matrix(ibra + nKets, iket + nKets) = phase * EtaME;
                Eta_matrix_c(ibra + nKets, iket + nKets) = occfactor_l * phase * EtaME;
                Eta_matrix_d(ibra + nKets, iket + nKets) = occfactor_k * phase * EtaME;
                if (Z_is_scalar)
                  Gamma_matrix(ibra + nKets, iket + nKets) = phase * GammaME;

                phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                Eta_matrix(ibra, iket + nKets) = phase * EtaME;
                Eta_matrix_c(ibra, iket + nKets) = occfactor_l * phase * EtaME;
                Eta_matrix_d(ibra, iket + nKets) = occfactor_k * phase * EtaME;
                if (Z_is_scalar)
                  Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
              }
            }
            else
            {
              if (k != l)
              {
                int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                Eta_matrix(ibra, iket + nKets) = phase * EtaME;
                Eta_matrix_c(ibra, iket + nKets) = occfactor_l * phase * EtaME;
                Eta_matrix_d(ibra, iket + nKets) = occfactor_k * phase * EtaME;
                if (Z_is_scalar)
                  Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
              }
            }
          }
        }
        // TODO: We can use symmetry here so that we don't have to use the full matrix.
        CHI_IV[ch] = Eta_matrix * Eta_matrix_c;
        CHI_IV[ch] += (Eta_matrix * Eta_matrix_d).t();
        if (Z_is_scalar)
          CHI_VII[ch] = Gamma_matrix * Eta_matrix_d + hEta * Eta_matrix_d.t() * Gamma_matrix;
      }

      //// Full Isospin tensor CHI_VII
      //// diagram IIIe and IIIf
      if (not Z_is_scalar)
      {
#pragma omp parallel for
        for (int ch = 0; ch < nch; ++ch)
        {
          int ch_bra = ch_bra_list[ch];
          int ch_ket = ch_ket_list[ch];
          TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
          TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
          size_t nbras = tbc_bra.GetNumberKets();
          size_t nKets = tbc_ket.GetNumberKets();
          if (nbras == 0 or nKets == 0)
            continue;
          int J0 = tbc_bra.J;

          arma::mat Eta_matrix_bra(2 * nbras, 2 * nbras);
          arma::mat Eta_matrix_ket(2 * nKets, 2 * nKets);
          arma::mat Gamma_matrix(2 * nbras, 2 * nKets);

          // full Gamma_matrix Eta_matrix_bra
          for (int ibra = 0; ibra < nbras; ++ibra)
          {
            Ket &bra = tbc_bra.GetKet(ibra);
            size_t i = bra.p;
            size_t j = bra.q;
            Orbit &oi = *(bra.op);
            Orbit &oj = *(bra.oq);
            int ji = oi.j2;
            int jj = oj.j2;
            double n_i = oi.occ;
            double bar_n_i = 1. - n_i;
            double n_j = oj.occ;
            double bar_n_j = 1. - n_j;

            // full gamma
            for (int iket = 0; iket < nKets; ++iket)
            {
              size_t k, l;
              Ket &ket = tbc_ket.GetKet(iket);
              k = ket.p;
              l = ket.q;
              Orbit &ok = Z.modelspace->GetOrbit(k);
              Orbit &ol = Z.modelspace->GetOrbit(l);
              int jk = ok.j2;
              int jl = ol.j2;
              double n_k = ok.occ;
              double bar_n_k = 1. - n_k;
              double n_l = ol.occ;
              double bar_n_l = 1. - n_l;

              double GammaME = Gamma.TwoBody.GetTBME_J(J0, i, j, k, l);
              Gamma_matrix(ibra, iket) = GammaME;
              if (i != j)
              {
                int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
                Gamma_matrix(ibra + nbras, iket) = phase * GammaME;
                if (k != l)
                {
                  phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
                  Gamma_matrix(ibra + nbras, iket + nKets) = phase * GammaME;

                  phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                  Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
                }
              }
              else
              {
                if (k != l)
                {
                  int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                  Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
                }
              }
            }

            // full Eta_matrix_bra
            for (int iket = 0; iket < nbras; ++iket)
            {
              size_t k, l;
              Ket &ket = tbc_bra.GetKet(iket);
              k = ket.p;
              l = ket.q;
              Orbit &ok = Z.modelspace->GetOrbit(k);
              Orbit &ol = Z.modelspace->GetOrbit(l);
              int jk = ok.j2;
              int jl = ol.j2;
              double n_k = ok.occ;
              double bar_n_k = 1. - n_k;
              double n_l = ol.occ;
              double bar_n_l = 1. - n_l;

              double occfactor_k = (bar_n_i * bar_n_j * n_k + n_i * n_j * bar_n_k);
              double occfactor_l = (bar_n_i * bar_n_j * n_l + n_i * n_j * bar_n_l);

              double EtaME = Eta.TwoBody.GetTBME_J(J0, i, j, k, l);

              Eta_matrix_bra(ibra, iket) = occfactor_l * EtaME;
              if (i != j)
              {
                int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
                Eta_matrix_bra(ibra + nbras, iket) = occfactor_l * phase * EtaME;
                if (k != l)
                {
                  phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
                  Eta_matrix_bra(ibra + nbras, iket + nbras) = occfactor_k * phase * EtaME;

                  phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                  Eta_matrix_bra(ibra, iket + nbras) = occfactor_k * phase * EtaME;
                }
              }
              else
              {
                if (k != l)
                {
                  int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                  Eta_matrix_bra(ibra, iket + nbras) = occfactor_k * phase * EtaME;
                }
              }
            }
          }

          // full Eta_matrix_ket
          for (int ibra = 0; ibra < nKets; ++ibra)
          {
            Ket &bra = tbc_ket.GetKet(ibra);
            size_t i = bra.p;
            size_t j = bra.q;
            Orbit &oi = *(bra.op);
            Orbit &oj = *(bra.oq);
            int ji = oi.j2;
            int jj = oj.j2;
            double n_i = oi.occ;
            double bar_n_i = 1. - n_i;
            double n_j = oj.occ;
            double bar_n_j = 1. - n_j;

            for (int iket = 0; iket < nKets; ++iket)
            {
              size_t k, l;
              Ket &ket = tbc_ket.GetKet(iket);
              k = ket.p;
              l = ket.q;
              Orbit &ok = Z.modelspace->GetOrbit(k);
              Orbit &ol = Z.modelspace->GetOrbit(l);
              int jk = ok.j2;
              int jl = ol.j2;
              double n_k = ok.occ;
              double bar_n_k = 1. - n_k;
              double n_l = ol.occ;
              double bar_n_l = 1. - n_l;
              double occfactor_k = (bar_n_i * bar_n_j * n_k + n_i * n_j * bar_n_k);
              double occfactor_l = (bar_n_i * bar_n_j * n_l + n_i * n_j * bar_n_l);

              double EtaME = Eta.TwoBody.GetTBME_J(J0, i, j, k, l);

              Eta_matrix_ket(ibra, iket) = occfactor_l * EtaME;
              if (i != j)
              {
                int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
                Eta_matrix_ket(ibra + nKets, iket) = occfactor_l * phase * EtaME;
                if (k != l)
                {
                  phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
                  Eta_matrix_ket(ibra + nKets, iket + nKets) = occfactor_k * phase * EtaME;

                  phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                  Eta_matrix_ket(ibra, iket + nKets) = occfactor_k * phase * EtaME;
                }
              }
              else
              {
                if (k != l)
                {
                  int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                  Eta_matrix_ket(ibra, iket + nKets) = occfactor_k * phase * EtaME;
                }
              }
            }
          }

          CHI_VII[ch] = Gamma_matrix * Eta_matrix_ket + hEta * Eta_matrix_bra.t() * Gamma_matrix;
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      //_______________________________________________________________________________
      std::deque<arma::mat> bar_CHI_IV(n_nonzero);     // released
      std::deque<arma::mat> bar_CHI_VII_CC(n_nonzero); // released
      /// build intermediate bar operator
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        if (nKets_cc < 1)
          continue;

        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        bar_CHI_IV[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        bar_CHI_VII_CC[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      /// Pandya transformation only recouple the angula momentum
      /// IIe and IIf                 barCHI_III_RC   bar_CHI_IV
      /// diagram IIIe and IIIf       bar_CHI_VII_CC
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        if (nKets_cc < 1)
        {
          continue;
        }

        int J_cc = tbc_cc.J;
        for (int ibra_cc = 0; ibra_cc < nKets_cc * 2; ++ibra_cc)
        {
          int a, b;
          if (ibra_cc < nKets_cc)
          {
            Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
            a = bra_cc.p;
            b = bra_cc.q;
          }
          else
          {
            Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nKets_cc);
            b = bra_cc.p;
            a = bra_cc.q;
          }
          if (ibra_cc >= nKets_cc and a == b)
            continue;
          Orbit &oa = Z.modelspace->GetOrbit(a);
          Orbit &ob = Z.modelspace->GetOrbit(b);
          double ja = oa.j2 * 0.5;
          double jb = ob.j2 * 0.5;

          // loop over cross-coupled kets |cd> in this channel
          for (int iket_cc = 0; iket_cc < nKets_cc * 2; ++iket_cc)
          {
            int c, d;
            if (iket_cc < nKets_cc)
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
              c = ket_cc_cd.p;
              d = ket_cc_cd.q;
            }
            else
            {
              Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nKets_cc);
              d = ket_cc_cd.p;
              c = ket_cc_cd.q;
            }
            if (iket_cc >= nKets_cc and c == d)
              continue;
            Orbit &oc = Z.modelspace->GetOrbit(c);
            Orbit &od = Z.modelspace->GetOrbit(d);
            double jc = oc.j2 * 0.5;
            double jd = od.j2 * 0.5;

            int Tz_J2_bc = (ob.tz2 + oc.tz2) / 2;
            int Tz_J2_ad = (oa.tz2 + od.tz2) / 2;
            int parity_J2 = (ob.l + oc.l) % 2;

            int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
            int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
            double XbarIIef = 0;
            double XbarIIIef = 0;
            for (int J_std = jmin; J_std <= jmax; J_std++)
            {
              int phaseFactor = Z.modelspace->phase(J_std + (oc.j2 + ob.j2) / 2);
              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                int ch_J2_bc = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_J2, Tz_J2_bc);
                int ch_J2_ad = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_J2, Tz_J2_ad);

                TwoBodyChannel &tbc_J2_bc = Z.modelspace->GetTwoBodyChannel(ch_J2_bc);
                TwoBodyChannel &tbc_J2_ad = Z.modelspace->GetTwoBodyChannel(ch_J2_ad);
                int nkets_bc = tbc_J2_bc.GetNumberKets();
                int nkets_ad = tbc_J2_ad.GetNumberKets();
                if (nkets_bc < 1 or nkets_ad < 1)
                  continue;

                int indx_bc = tbc_J2_bc.GetLocalIndex(std::min(int(b), int(c)), std::max(int(b), int(c)));
                int indx_ad = tbc_J2_ad.GetLocalIndex(std::min(int(a), int(d)), std::max(int(a), int(d)));

                if (indx_ad < 0 or indx_bc < 0)
                  continue;
                if (a > d)
                  indx_ad += nkets_ad;

                int indx_cb = indx_bc;
                if (b > c)
                  indx_bc += nkets_bc;
                if (c > b)
                  indx_cb += nkets_bc;

                int index_ch = ch_J2_ad;
                if (not Z_is_scalar)
                {
                  if (ch_J2_ad > ch_J2_bc)
                  {
                    index_ch = ch_J2_ad;
                    ch_J2_ad = ch_J2_bc;
                    ch_J2_bc = index_ch;
                  }
                  index_ch = -1;
                  for (size_t i = 0; i < nch; i++)
                  {
                    if (ch_bra_list[i] == ch_J2_ad and ch_ket_list[i] == ch_J2_bc)
                    {
                      index_ch = i;
                    }
                  }
                }

                XbarIIef -= (2 * J_std + 1) * sixj1 * CHI_IV[ch_J2_bc](indx_ad, indx_cb);
                if (std::abs(Tz_J2_ad - Tz_J2_bc) == Z.TwoBody.rank_T)
                  XbarIIIef += phaseFactor * (2 * J_std + 1) * sixj1 * CHI_VII[index_ch](indx_ad, indx_bc);
              }
            }
            bar_CHI_IV[ch_cc](ibra_cc, iket_cc) = XbarIIef;
            bar_CHI_VII_CC[ch_cc](ibra_cc, iket_cc) = XbarIIIef;
          }
          //-------------------
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      /// release memory
      for (int ch = 0; ch < nch_eta; ++ch)
      {
        CHI_IV[ch].clear();
      }
      for (int ch = 0; ch < nch; ++ch)
      {
        CHI_VII[ch].clear();
      }
      CHI_IV.clear();
      CHI_VII.clear();

      /////////////////////////////////////////////
      //     diagram    IIe and IIf
      /////////////////////////////////////////////
      std::deque<arma::mat> bar_CHI_gamma(n_nonzero); // released
      /// initial bar_CHI_V
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        bar_CHI_gamma[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      // calculate bat_chi_IV * bar_gamma
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        bar_CHI_gamma[ch_cc] = bar_CHI_IV[ch_cc] * bar_Gamma[ch_cc];
      }
      // release memroy
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        bar_CHI_IV[ch_cc].clear();
        //bar_Gamma[ch_cc].clear();
      }
      bar_CHI_IV.clear();
      //bar_Gamma.clear();

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      //  Inverse Pandya transformation
      //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
      //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
      //                                    { k l J'}
      //  Diagram e
      //  II(e)^J_ijkl  = - 1/2 ( 1- P_ij ) ( 1- P_kl ) sum_J' (2J'+1)  { i j J }  \bar{bar_CHI_gamma}^J'_il`kj`
      //                                                                { k l J'}
      //  Diagram f
      //  II(f)^J_ijkl  = - 1/2 ( 1- P_ij ) ( 1- P_kl ) sum_J' (2J'+1)  { i j J }  \bar{bar_CHI_gamma_II}^J'_il`kj`
      //
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        size_t ch_bra = ch_bra_list[ch];
        size_t ch_ket = ch_ket_list[ch];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nKets = tbc_ket.GetNumberKets();
        if (nbras == 0 or nKets == 0)
          continue;

        int J0 = tbc_bra.J;
        for (int ibra = 0; ibra < nbras; ++ibra)
        {
          Ket &bra = tbc_bra.GetKet(ibra);
          size_t i = bra.p;
          size_t j = bra.q;
          Orbit &oi = *(bra.op);
          Orbit &oj = *(bra.oq);
          int ji = oi.j2;
          int jj = oj.j2;

          int ketmin = 0;
          if (ch_bra == ch_ket)
            ketmin = ibra;
          for (int iket = ketmin; iket < nKets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc_ket.GetKet(iket);
            k = ket.p;
            l = ket.q;

            Orbit &ok = Z.modelspace->GetOrbit(k);
            Orbit &ol = Z.modelspace->GetOrbit(l);
            int jk = ok.j2;
            int jl = ol.j2;
            double commijkl = 0;
            double commjikl = 0;
            double commijlk = 0;
            double commjilk = 0;

            // ijkl direct term
            int parity_cc = (oi.l + ol.l) % 2;
            int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
            int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
            int Jpmax = std::min(ji + jl, jj + jk) / 2;
            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              if (std::abs(sixj) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 1)
                continue;

              int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
              int indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k));
              if (indx_il < 0 or indx_kj < 0)
                continue;
              // jilk, exchange i and j, k and l     ->  jk   li
              int indx_jk = indx_kj + (j > k ? nkets_cc : 0);
              int indx_li = indx_il + (l > i ? nkets_cc : 0);
              double me1 = bar_CHI_gamma[ch_cc](indx_jk, indx_li);
              commjilk -= (2 * Jprime + 1) * sixj * me1;

              // ijkl direct term
              indx_il += (i > l ? nkets_cc : 0);
              indx_kj += (k > j ? nkets_cc : 0);
              me1 = bar_CHI_gamma[ch_cc](indx_il, indx_kj);
              commijkl -= (2 * Jprime + 1) * sixj * me1;
            }

            // jikl, exchange i and j    ->  jl ki
            parity_cc = (oi.l + ok.l) % 2;
            Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
            Jpmin = std::max(std::abs(int(jj - jl)), std::abs(int(jk - ji))) / 2;
            Jpmax = std::min(int(jj + jl), int(jk + ji)) / 2;
            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              if (std::abs(sixj) < 1e-8)
                continue;

              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 1)
                continue;

              int indx_ki = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
              int indx_jl = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
              if (indx_ki < 0 or indx_jl < 0)
                continue;

              // ijlk, exchange k and l     ->  ik lj
              int indx_ik = indx_ki + (i > k ? nkets_cc : 0);
              int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
              double me1 = bar_CHI_gamma[ch_cc](indx_ik, indx_lj);
              commijlk -= (2 * Jprime + 1) * sixj * me1;

              // jikl, exchange i and j    ->  jl ki
              indx_ki += (k > i ? nkets_cc : 0);
              indx_jl += (j > l ? nkets_cc : 0);
              me1 = bar_CHI_gamma[ch_cc](indx_jl, indx_ki);
              commjikl -= (2 * Jprime + 1) * sixj * me1;
            }

            double zijkl = (commijkl - Z.modelspace->phase((ji + jj) / 2 - J0) * commjikl);
            zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commijlk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commjilk);

            if (i == j)
              zijkl /= PhysConst::SQRT2;
            if (k == l)
              zijkl /= PhysConst::SQRT2;

            Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zijkl);
          }
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)

      {
        bar_CHI_gamma[ch_cc].clear();
      }
      bar_CHI_gamma.clear();

      // ######################################################
      //                 Diagram IIb and IId
      // ######################################################
      std::deque<arma::mat> CHI_III_final(n_nonzero);
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc_bra = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nbras = tbc_cc_bra.GetNumberKets();
        if (nbras < 1)
          continue;
        // Not symmetric
        CHI_III_final[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
      }

#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        CHI_III_final[ch_cc] = bar_Gamma[ch_cc] * barCHI_III_RC[ch_cc];
      }
      /// release memory
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        barCHI_III_RC[ch_cc].clear();
        bar_Gamma[ch_cc].clear();
      }
      barCHI_III_RC.clear();
      bar_Gamma.clear();


      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }
      //  Inverse Pandya transformation
      //  diagram IIb and IId
      //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
      //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
      //                                    { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        size_t ch_bra = ch_bra_list[ch];
        size_t ch_ket = ch_ket_list[ch];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nKets = tbc_ket.GetNumberKets();
        if (nbras == 0 or nKets == 0)
          continue;

        int J0 = tbc_bra.J;
        for (int ibra = 0; ibra < nbras; ++ibra)
        {
          Ket &bra = tbc_bra.GetKet(ibra);
          size_t i = bra.p;
          size_t j = bra.q;
          Orbit &oi = *(bra.op);
          Orbit &oj = *(bra.oq);
          int ji = oi.j2;
          int jj = oj.j2;
          int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

          int ketmin = 0;
          if (ch_bra == ch_ket)
            ketmin = ibra;
          for (int iket = ketmin; iket < nKets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc_ket.GetKet(iket);
            k = ket.p;
            l = ket.q;

            Orbit &ok = Z.modelspace->GetOrbit(k);
            Orbit &ol = Z.modelspace->GetOrbit(l);
            int jk = ok.j2;
            int jl = ol.j2;
            double commijkl = 0;
            double commjikl = 0;
            double commijlk = 0;
            double commjilk = 0;

            // jikl, direct term        -->  jl  ki
            // ijlk, exchange ij and kl -->  lj  ik
            int parity_cc = (oi.l + ok.l) % 2;
            int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
            int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
            int Jpmax = std::min(jj + jl, ji + jk) / 2;

            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj1 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              if (std::abs(sixj1) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);

              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 1)
                continue;

              int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
              int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
              if (indx_jl < 0 or indx_ik < 0)
                continue;

              int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
              // direct term
              indx_jl += (j > l ? nkets_cc : 0);
              indx_ik += (i > k ? nkets_cc : 0);
              double me1 = CHI_III_final[ch_cc](indx_jl, indx_ik);
              commjikl -= phase1 * (2 * Jprime + 1) * sixj1 * me1;

              int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);
              // exchange ij and kl
              double me2 = CHI_III_final[ch_cc](indx_ik, indx_jl);
              commijlk -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
            }

            // ijkl,  exchange i and j -->  il  kj
            // jilk,  exchange k and l -->  jk li
            parity_cc = (oi.l + ol.l) % 2;
            Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
            Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
            Jpmax = std::min(ji + jl, jj + jk) / 2;

            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj1 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              if (std::abs(sixj1) < 1e-8)
                continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 1)
                continue;

              int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
              int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
              if (indx_il < 0 or indx_jk < 0)
                continue;

              int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
              // exchange k and l
              indx_il += (i > l ? nkets_cc : 0);
              indx_jk += (j > k ? nkets_cc : 0);
              double me1 = CHI_III_final[ch_cc](indx_jk, indx_il);
              commjilk -= phase1 * (2 * Jprime + 1) * sixj1 * me1;

              int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);
              // exchange i and j
              double me2 = CHI_III_final[ch_cc](indx_il, indx_jk);
              commijkl -= phase2 * (2 * Jprime + 1) * sixj1 * me2;
            }

            double zijkl = (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
            zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

            if (i == j)
              zijkl /= PhysConst::SQRT2;
            if (k == l)
              zijkl /= PhysConst::SQRT2;

            Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, phaseFactor * zijkl);
          }
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      for (int ch = 0; ch < n_nonzero; ++ch)
      {
        CHI_III_final[ch].clear();
      }
      CHI_III_final.clear();


      /////////////////////////////////////////////
      //     diagram    IIIe and IIIf
      /////////////////////////////////////////////
      std::deque<arma::mat> bar_CHI_VII_CC_ef(n_nonzero); // released
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        bar_CHI_VII_CC_ef[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      // bar_CHI_VII_CC_ef = bar_CHI_VII_CC * bar_Eta
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets = tbc_cc.GetNumberKets();
        if (nKets < 1)
          continue;
        bar_CHI_VII_CC_ef[ch_cc] = bar_CHI_VII_CC[ch_cc] * bar_Eta[ch_cc];
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      // release bar_Eta
      // release bar_CHI_VII_CC
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        bar_Eta[ch_cc].clear();
        bar_CHI_VII_CC[ch_cc].clear();
      }
      bar_Eta.clear();
      bar_CHI_VII_CC.clear();

      //  Inverse Pandya transformation
      //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
      //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
      //                                    { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        int ch_bra = ch_bra_list[ch];
        int ch_ket = ch_ket_list[ch];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();
        int J0 = tbc_bra.J;
        if (nbras == 0 or nkets == 0)
          continue;

        for (int ibra = 0; ibra < nbras; ++ibra)
        {
          Ket &bra = tbc_bra.GetKet(ibra);
          size_t i = bra.p;
          size_t j = bra.q;
          Orbit &oi = *(bra.op);
          Orbit &oj = *(bra.oq);
          int ji = oi.j2;
          int jj = oj.j2;
          int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

          int ketmin = 0;
          if (ch_bra == ch_ket)
            ketmin = ibra;
          for (int iket = ketmin; iket < nkets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc_ket.GetKet(iket);
            k = ket.p;
            l = ket.q;

            Orbit &ok = Z.modelspace->GetOrbit(k);
            Orbit &ol = Z.modelspace->GetOrbit(l);
            int jk = ok.j2;
            int jl = ol.j2;
            double commijkl = 0;
            double commjikl = 0;
            double commijlk = 0;
            double commjilk = 0;

            // jikl, direct term        -->  jl  ki
            // ijlk, exchange ij and kl -->  lj  ik
            int parity_cc = (oi.l + ok.l) % 2;
            int Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
            int Jpmin = std::max(std::abs(jj - jl), std::abs(ji - jk)) / 2;
            int Jpmax = std::min(jj + jl, ji + jk) / 2;
            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              double sixj2 = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jl * 0.5, jk * 0.5, Jprime);
              if (std::abs(sixj) < 1e-8)
                continue;

              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 1)
                continue;

              int indx_lj = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
              int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
              if (indx_lj < 0 or indx_ik < 0)
                continue;

              // direct term
              int indx_jl = indx_lj + (j > l ? nkets_cc : 0);
              int indx_ki = indx_ik + (k > i ? nkets_cc : 0);
              double me1 = bar_CHI_VII_CC_ef[ch_cc](indx_jl, indx_ki);
              commjikl -= (2 * Jprime + 1) * sixj * me1;

              // exchange ij and kl
              indx_ik += (i > k ? nkets_cc : 0);
              indx_lj += (l > j ? nkets_cc : 0);
              double me2 = bar_CHI_VII_CC_ef[ch_cc](indx_ik, indx_lj);
              commijlk -= (2 * Jprime + 1) * sixj2 * me2;
            }

            // ijkl,  exchange i and j -->  il  kj
            // jilk,  exchange k and l -->  jk li
            parity_cc = (oi.l + ol.l) % 2;
            Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
            Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
            Jpmax = std::min(ji + jl, jj + jk) / 2;
            for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
              double sixj2 = Z.modelspace->GetSixJ(jj * 0.5, ji * 0.5, J0, jl * 0.5, jk * 0.5, Jprime);
              if (std::abs(sixj) < 1e-8)
                continue;

              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);

              int nkets_cc = tbc_cc.GetNumberKets();
              if (nkets_cc < 0)
                continue;

              int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
              int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
              if (indx_il < 0 or indx_kj < 0)
                continue;

              // exchange k and l
              int indx_jk = indx_kj + (j > k ? nkets_cc : 0);
              int indx_li = indx_il + (l > i ? nkets_cc : 0);
              double me2 = bar_CHI_VII_CC_ef[ch_cc](indx_jk, indx_li);
              commjilk -= (2 * Jprime + 1) * sixj2 * me2;

              // exchange i and j
              indx_il += (i > l ? nkets_cc : 0);
              indx_kj += (k > j ? nkets_cc : 0);
              double me1 = bar_CHI_VII_CC_ef[ch_cc](indx_il, indx_kj);
              commijkl -= (2 * Jprime + 1) * sixj * me1;
            }

            double zijkl = (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
            zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

            if (i == j)
              zijkl /= PhysConst::SQRT2;
            if (k == l)
              zijkl /= PhysConst::SQRT2;

            Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, phaseFactor * 0.5 * zijkl);
          }
        }
      }

      if (Commutator::verbose)
      {
        Z.profiler.timer["_" + std::string(__func__) + "_" + std::to_string(__LINE__)] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
      }

      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        bar_CHI_VII_CC_ef[ch_cc].clear();
      }
      bar_CHI_VII_CC_ef.clear();

      // Timer
      Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
      return;
    } //  comm223_232_chi2b

  } // namespace FactorizedDoubleCommutator
} // namespace Commutator
