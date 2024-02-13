
#include "FactorizedDoubleCommutator.hh"
#include "PhysicalConstants.hh"



namespace Commutator
{


 namespace FactorizedDoubleCommutator
 {


  bool use_goose_tank_1b = true;       // always include
  bool use_goose_tank_2b = true;       // always include

  bool use_1b_intermediates;
  bool use_2b_intermediates;
  bool use_goose_tank_only_1b = false; // only calculate Goose Tanks
  bool use_goose_tank_only_2b = false; // only calculate Goose Tanks
  bool use_TypeII_1b = true;
  bool use_TypeIII_1b = true;
  bool use_TypeII_2b = true;
  bool use_TypeIII_2b = true;
  bool use_GT_TypeI_2b = true;
  bool use_GT_TypeIV_2b = true;
  bool SlowVersion = false;


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



  void UseSlowVersion(bool tf)
  {
    SlowVersion = tf;
  }

  // factorize double commutator [Eta, [Eta, Gamma]_3b ]_1b
  void comm223_231(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {

    if( use_1b_intermediates)
    {
      comm223_231_chi1b( Eta, Gamma, Z);  // topology with 1-body intermediate (fast)
    }
    if( use_2b_intermediates)
    {
      comm223_231_chi2b( Eta, Gamma, Z);  // topology with 2-body intermediate (slow)
    }
    return;
  } //comm223_231




////////////////////////////////////////////////////////////////////////////
/// factorized 223_231 double commutator with 1b intermediate
////////////////////////////////////////////////////////////////////////////
  void comm223_231_chi1b(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {

    double t_internal = omp_get_wtime(); // timer
    double t_start = omp_get_wtime(); // timer

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

    std::deque<arma::mat> Eta_matrix(nch);      // used
    std::deque<arma::mat> Eta_matrix_nnnn(nch); // used
    std::deque<arma::mat> Gamma_matrix(nch);
    std::deque<arma::mat> intermediateTB(nch);
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      Eta_matrix[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      Eta_matrix_nnnn[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      Gamma_matrix[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      intermediateTB[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

// full matrix
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
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

          double occfactor_k = (bar_n_k * bar_n_l * n_i * n_j - bar_n_i * bar_n_j * n_k * n_l);
          double EtaME = Eta.TwoBody.GetTBME_J(J0, i, j, k, l);
          double GammaME = Gamma.TwoBody.GetTBME_J(J0, i, j, k, l);
          Eta_matrix[ch](ibra, iket) = EtaME;
          Gamma_matrix[ch](ibra, iket) = GammaME;
          Eta_matrix_nnnn[ch](ibra, iket) = occfactor_k * EtaME;
          if (i != j)
          {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Eta_matrix[ch](ibra + nKets, iket) = phase * EtaME;
            Gamma_matrix[ch](ibra + nKets, iket) = phase * GammaME;
            Eta_matrix_nnnn[ch](ibra + nKets, iket) = occfactor_k * phase * EtaME;
            if (k != l)
            {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Eta_matrix[ch](ibra + nKets, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra + nKets, iket + nKets) = phase * GammaME;
              Eta_matrix_nnnn[ch](ibra + nKets, iket + nKets) = occfactor_k * phase * EtaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix[ch](ibra, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra, iket + nKets) = phase * GammaME;
              Eta_matrix_nnnn[ch](ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
          else
          {
            if (k != l)
            {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix[ch](ibra, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra, iket + nKets) = phase * GammaME;
              Eta_matrix_nnnn[ch](ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
        }
      }
    }

    Z.profiler.timer["_231_F_loop1"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

//    if (use_goose_tank)
    if (use_goose_tank_1b)
    {
#pragma omp parallel for 
      for (int ch = 0; ch < nch; ++ch)
      {
//        std::cout << " ch = " << ch << " dimensions " << Eta_matrix_nnnn[ch].n_rows << " x " << Eta_matrix_nnnn[ch].n_cols
//                  << "    " << Eta_matrix[ch].n_rows << " x " << Eta_matrix[ch].n_cols << std::endl;
        intermediateTB[ch] = Eta_matrix_nnnn[ch] * Eta_matrix[ch];
        intermediateTB[ch] -= Eta_matrix[ch] * Eta_matrix_nnnn[ch].t();  // Can this be obtained from the previous line without another matmult?
      }

    Z.profiler.timer["_231_F_intermediateTB"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

#pragma omp parallel for
      for (int indexd = 0; indexd < norbits; ++indexd)
      {
        auto d = allorb_vec[indexd];
        Orbit &od = Z.modelspace->GetOrbit(d);
        double dj2 = od.j2 + 1.0;
        for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
        {
          if (e > d)
            continue;
          double eta_de = 0;
          for (int ch = 0; ch < nch; ++ch)
          {
            TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
            int J0 = tbc.J;
            int nKets = tbc.GetNumberKets();

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);

              int indx_bd = tbc.GetLocalIndex(std::min(int(b), int(d)), std::max(int(b), int(d)));
              int indx_be = tbc.GetLocalIndex(std::min(int(b), int(e)), std::max(int(b), int(e)));
              if (indx_bd < 0 or indx_be < 0)
              {
                continue;
              }
              if (b > d)
                indx_bd += nKets;
              if (b > e)
                indx_be += nKets;
              double doubleEta = (2 * J0 + 1) * intermediateTB[ch](indx_bd, indx_be);
              eta_de += doubleEta;
            }
          }
          Chi_221_a(d, e) += eta_de / dj2;
          if (d != e)
            Chi_221_a(e, d) += eta_de / dj2;
        } // e
      }   // d

    Z.profiler.timer["_231_F_indexd"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

#pragma omp parallel for
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
            for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
            {
              Orbit &oe = Z.modelspace->GetOrbit(e);

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
      }   // for p
      // std::cout << "diagram I  " << Z.OneBodyNorm() << std::endl;
      // Z.EraseOneBody();
    }
    Z.profiler.timer["_231_F_indexp"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

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

      for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
      {
        // if (e > d)
        //   continue;
        Orbit &oe = Z.modelspace->GetOrbit(e);
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
    }   // d

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
    }   // for p

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return;
  } //comm223_231_chi1b



////////////////////////////////////////////////////////////////////////////
/// factorized 223_231 double commutator with 2b intermediate
////////////////////////////////////////////////////////////////////////////
  void comm223_231_chi2b(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {

    double t_internal = omp_get_wtime(); // timer
    double t_start = omp_get_wtime(); // timer

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


    int nch = Z.modelspace->GetNumberTwoBodyChannels();
    int norbits = Z.modelspace->all_orbits.size();
    std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());

//    std::deque<arma::mat> Eta_matrix(nch);      // used
//    std::deque<arma::mat> Eta_matrix_nnnn(nch); // used
//    std::deque<arma::mat> Gamma_matrix(nch);
    std::deque<arma::mat> intermediateTB(nch);
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
//      Eta_matrix[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
//      Eta_matrix_nnnn[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
//      Gamma_matrix[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      intermediateTB[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

// fill  Gamma_matrix Eta_matrix, Eta_matrix_nnnn
#pragma omp parallel for schedule(dynamic,1)
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();

      arma::mat Eta_matrix = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      arma::mat Eta_matrix_nnnn = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      arma::mat Gamma_matrix = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);

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

          double occfactor_k = (bar_n_k * bar_n_l * n_i * n_j - bar_n_i * bar_n_j * n_k * n_l);
//          double EtaME = Eta.TwoBody.GetTBME_J(J0, i, j, k, l);
//          double GammaME = Gamma.TwoBody.GetTBME_J(J0, i, j, k, l);
          double EtaME = Eta.TwoBody.GetTBME(ch,ch,bra,ket);
          double GammaME = Gamma.TwoBody.GetTBME(ch,ch,bra,ket);
          Eta_matrix(ibra, iket) = EtaME;  // Eta_ijkl
          Gamma_matrix(ibra, iket) = GammaME;
          Eta_matrix_nnnn(ibra, iket) = occfactor_k * EtaME;
          if (i != j)
          {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Eta_matrix(ibra + nKets, iket) = phase * EtaME;  // Eta_jikl
            Gamma_matrix(ibra + nKets, iket) = phase * GammaME;
            Eta_matrix_nnnn(ibra + nKets, iket) = occfactor_k * phase * EtaME;
            if (k != l)
            {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Eta_matrix(ibra + nKets, iket + nKets) = phase * EtaME; // Eta_jilk
              Gamma_matrix(ibra + nKets, iket + nKets) = phase * GammaME;
              Eta_matrix_nnnn(ibra + nKets, iket + nKets) = occfactor_k * phase * EtaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix(ibra, iket + nKets) = phase * EtaME;  // Eta_ijlk
              Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
              Eta_matrix_nnnn(ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
          else
          {
            if (k != l)
            {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix(ibra, iket + nKets) = phase * EtaME; // Eta_ijlk
              Gamma_matrix(ibra, iket + nKets) = phase * GammaME;
              Eta_matrix_nnnn(ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
        }// for iket
      }// for ibra

      arma::mat Chi_222_b = (2 * J0 + 1) * Eta_matrix * Eta_matrix_nnnn;

      intermediateTB[ch] = Chi_222_b * Gamma_matrix;
      intermediateTB[ch] += Gamma_matrix * Chi_222_b.t();
    }// for ch

    Z.profiler.timer["_231_F_chi2_pp_fill_chi"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();


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
    //  diagram II_b and II_d
    //
    // IIb_pq = 1/4 1/(2 jp + 1) \sum_acdJ0 Chi_222_b_cpad * Gamma_bar_adcq
    // IIb_pq = - 1/4 1/(2 jp + 1) \sum_abeJ0  Chi_222_a_bqae * Gamma_bar_bqae
    // ###########################################################
#pragma omp parallel for schedule(dynamic,1)
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
        double zij = 0;

        // loop abcde
        for (auto &c : Z.modelspace->all_orbits)
        {
          Orbit &oc = Z.modelspace->GetOrbit(c);
          double jc = oc.j2 / 2.;

          int J0min = std::abs(oc.j2 - op.j2) / 2;
          int J0max = (oc.j2 + op.j2) / 2;
          int parity = (op.l + oc.l) % 2;
          int Tz_J0 = (op.tz2 + oc.tz2) / 2;
          for (int J0 = J0min; J0 <= J0max; J0++)
          {
            int ch_J0 = Z.modelspace->GetTwoBodyChannelIndex(J0, parity, Tz_J0);
            TwoBodyChannel &tbc_J0 = Z.modelspace->GetTwoBodyChannel(ch_J0);
            int nkets = tbc_J0.GetNumberKets();
            int indx_cp = tbc_J0.GetLocalIndex(std::min(int(c), int(p)), std::max(int(c), int(p)));
            int indx_cq = tbc_J0.GetLocalIndex(std::min(int(c), int(q)), std::max(int(c), int(q)));
            if (indx_cp < 0 or indx_cq < 0)
            {
              continue;
            }
            if (c > p)
              indx_cp += nkets;
            if (c > q)
              indx_cq += nkets;

            zij += intermediateTB[ch_J0](indx_cp, indx_cq);
          }
        }
        Z.OneBody(p, q) += 0.25 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += 0.25 * hZ * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p


    Z.profiler.timer["_231_F_chi2_pp_fill1b"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();



    // *********************************************************************************** //
    //                                  Diagram II                                         //
    // *********************************************************************************** //

    // ###########################################################
    //  diagram II_a
    //
    //  Pandya transform
    //  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
    //                        { k l J'}
    int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();

    std::deque<arma::mat> IntermediateTwobody(n_nonzero);
/// Pandya transformation
#pragma omp parallel for schedule(dynamic,1)
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      int J_cc = tbc_cc.J;

      arma::mat Eta_bar = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
      arma::mat Eta_bar_nnnn = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED 
      arma::mat Gamma_bar = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
//      Gamma_bar[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros); // SRS ADDED
      arma::mat Chi_222_a;
      IntermediateTwobody[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      // transform operator
      // loop over cross-coupled ph bras <ab| in this channel
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
        double n_a = oa.occ;
        double nbar_a = 1 - n_a;
        double ja = oa.j2 * 0.5;

        Orbit &ob = Z.modelspace->GetOrbit(b);
        double n_b = ob.occ;
        double nbar_b = 1 - n_b;
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
          Orbit &oc = Z.modelspace->GetOrbit(c);
          double n_c = oc.occ;
          double nbar_c = 1 - n_c;
          double jc = oc.j2 * 0.5;

          Orbit &od = Z.modelspace->GetOrbit(d);
          double n_d = od.occ;
          double nbar_d = 1 - n_d;
          double jd = od.j2 * 0.5;

          double occ_factor = nbar_c * nbar_b * n_a * n_d - n_c * n_b * nbar_a * nbar_d;
          if (iket_cc >= nKets_cc and c == d)
            continue;

          // Check the isospin projection. If this isn't conserved in the usual channel,
          // then all the xcbad and yadcb will be zero and we don't need to bother computing SixJs.
          // if (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != 0)
          //   continue;

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
              Xbar -= (2 * J_std + 1) * sixj1 *   Eta.TwoBody.GetTBME_J(J_std, a, d, c, b);
              Ybar -= (2 * J_std + 1) * sixj1 * Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
            }
          }
          Eta_bar(ibra_cc, iket_cc) = Xbar;
          Eta_bar_nnnn(ibra_cc, iket_cc) = Xbar * occ_factor;
          Gamma_bar(ibra_cc, iket_cc) = Ybar;
        }// for iket-cc

        //-------------------
      }// for ibra_cc
      // SRS EDIT
      IntermediateTwobody[ch_cc] = (2*J_cc+1) * Eta_bar * Eta_bar_nnnn * Gamma_bar;
    }// for ch_cc

    Z.profiler.timer["_231_F_chi2_ph_fill_chi"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
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
#pragma omp parallel for  schedule(dynamic,1)
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
            int iket_cc = tbc_cc.GetLocalIndex(p, e);
            int jket_cc = tbc_cc.GetLocalIndex(q, e);
            int iket_cc2 = tbc_cc.GetLocalIndex(e, q);
            int jket_cc2 = tbc_cc.GetLocalIndex(e, p);
            zij += IntermediateTwobody[ch_cc](iket_cc, jket_cc);
            zij -= IntermediateTwobody[ch_cc](iket_cc2, jket_cc2);
          }
          Z.OneBody(p, q) += zij / j2hat2;
          if (p != q)
            Z.OneBody(q, p) += hZ * zij / j2hat2;
        }
      }
    }

    Z.profiler.timer["_231_F_chi2_ph_fill1b"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();


    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return;
  } //comm223_231_chi2b















////////////////////////////////////////////////////////////

void comm223_231_slow(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {
    double t_internal = omp_get_wtime(); // timer
    double t_start = omp_get_wtime(); // timer

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

    if (use_goose_tank_1b)
    {
#pragma omp parallel for
      for (int indexd = 0; indexd < norbits; ++indexd)
      {
        auto d = allorb_vec[indexd];
        Orbit &od = Z.modelspace->GetOrbit(d);
        double n_d = od.occ;
        double nbar_d = 1.0 - n_d;

        for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
        {
          if (e > d)
            continue;
          Orbit &oe = Z.modelspace->GetOrbit(e);
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
              int a = bra.p;
              int c = bra.q;

              Orbit &oa = Z.modelspace->GetOrbit(a);
              double n_a = oa.occ;
              double nbar_a = 1.0 - n_a;

              Orbit &oc = Z.modelspace->GetOrbit(c);
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &b : Z.modelspace->all_orbits)
              {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double n_b = ob.occ;
                double nbar_b = 1.0 - n_b;

                double occfactor = (nbar_a * nbar_c * n_b * n_d - nbar_b * nbar_d * n_a * n_c - nbar_b * nbar_e * n_a * n_c + nbar_a * nbar_c * n_b * n_e);
                if (std::abs(occfactor) < 1e-6)
                  continue;
                double doubleEta = (2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, d, a, c) * Eta.TwoBody.GetTBME_J(J0, J0, a, c, b, e);
                eta_de += doubleEta;
                if (a != c)
                  eta_de += doubleEta;
              }
            }
          }
          Chi_221_a(d, e) += eta_de / (od.j2 + 1.0);
          if (d != e)
            Chi_221_a(e, d) += eta_de / (od.j2 + 1.0);
        } // e
      }   // d

#pragma omp parallel for
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
            for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
            {
              Orbit &oe = Z.modelspace->GetOrbit(e);

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
      }   // for p
      // std::cout << "diagram I  " << Z.OneBodyNorm() << std::endl;
      // Z.EraseOneBody();
    }
    if (use_goose_tank_only_1b)
    {
      Z.profiler.timer[__func__] += omp_get_wtime() - t_internal;
      return;
    }

    // *********************************************************************************** //
    //                                  Diagram II                                         //
    // *********************************************************************************** //


    Z.profiler.timer["231_diagram_I"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime(); // timer
    if (use_TypeII_1b)
    {
      // ###########################################################
      //  diagram II_a
      //
      //  Pandya transform
      //  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
      //                        { k l J'}
      int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC();
      std::deque<arma::mat> Eta_bar(n_nonzero);
      std::deque<arma::mat> Gamma_bar(n_nonzero);
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        Eta_bar[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        Gamma_bar[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        int J_cc = tbc_cc.J;
        // transform operator
        // loop over cross-coupled ph bras <ab| in this channel
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
          double ja = oa.j2 * 0.5;

          Orbit &ob = Z.modelspace->GetOrbit(b);
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
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double jc = oc.j2 * 0.5;

            Orbit &od = Z.modelspace->GetOrbit(d);
            double jd = od.j2 * 0.5;

            if (iket_cc >= nKets_cc and c == d)
              continue;

            // Check the isospin projection. If this isn't conserved in the usual channel,
            // then all the xcbad and yadcb will be zero and we don't need to bother computing SixJs.
            // if (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != 0)
            //   continue;

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
            Eta_bar[ch_cc](ibra_cc, iket_cc) = Xbar;
            Gamma_bar[ch_cc](ibra_cc, iket_cc) = Ybar;
          }

          //-------------------
        }
      }

      Z.profiler.timer["231_diagram_IIa"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime(); // timer
      
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
      //              \bar{eta}_pedc * \bar{eta}_cdab
      std::deque<arma::mat> Chi_222_a(n_nonzero);
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        Chi_222_a[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      //  (nbar_e * nbar_d * n_f * n_c - nbar_f * nbar_c * n_e * n_d ) <ab|cd> <cd|ef> in cross-coupled
#pragma omp parallel for
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        int J3 = tbc_cc.J;
        //----------------------------------
        // transform operator
        // loop over cross-coupled ph bras <ab| in this channel
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
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          Orbit &ob = Z.modelspace->GetOrbit(b);
          double n_b = ob.occ;
          double nbar_b = 1.0 - n_b;

          // loop over cross-coupled kets |cd> in this channel
          for (int jket_cc = 0; jket_cc < nKets_cc * 2; ++jket_cc)
          {
            int e, f;
            if (jket_cc < nKets_cc)
            {
              Ket &ket_cc_ef = tbc_cc.GetKet(jket_cc);
              e = ket_cc_ef.p;
              f = ket_cc_ef.q;
            }
            else
            {
              Ket &ket_cc_ef = tbc_cc.GetKet(jket_cc - nKets_cc);
              f = ket_cc_ef.p;
              e = ket_cc_ef.q;
            }
            if (jket_cc >= nKets_cc and e == f)
              continue;

            Orbit &oe = Z.modelspace->GetOrbit(e);
            double n_e = oe.occ;
            double nbar_e = 1.0 - n_e;

            Orbit &of = Z.modelspace->GetOrbit(f);
            double n_f = of.occ;
            double nbar_f = 1.0 - n_f;

            double chi_ME = 0.;
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
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              Orbit &od = Z.modelspace->GetOrbit(d);
              double n_d = od.occ;
              double nbar_d = 1.0 - n_d;

              double occfactor = (nbar_e * nbar_d * n_f * n_c - nbar_f * nbar_c * n_e * n_d);
              chi_ME += occfactor * (2 * J3 + 1) * Eta_bar[ch_cc](ibra_cc, iket_cc) * Eta_bar[ch_cc](iket_cc, jket_cc);
            }
            Chi_222_a[ch_cc](ibra_cc, jket_cc) = chi_ME;
          }
        }
      }

      std::deque<arma::mat> IntermediateTwobody(n_nonzero);
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        IntermediateTwobody[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

#pragma omp parallel for
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        IntermediateTwobody[ch_cc] = Chi_222_a[ch_cc] * Gamma_bar[ch_cc];
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

#pragma omp parallel for
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
              int iket_cc = tbc_cc.GetLocalIndex(p, e);
              int jket_cc = tbc_cc.GetLocalIndex(q, e);
              int iket_cc2 = tbc_cc.GetLocalIndex(e, q);
              int jket_cc2 = tbc_cc.GetLocalIndex(e, p);
              zij += IntermediateTwobody[ch_cc](iket_cc, jket_cc);   // II_a
              zij -= IntermediateTwobody[ch_cc](iket_cc2, jket_cc2); // II_c
            }
            Z.OneBody(p, q) += zij / j2hat2;
            if (p != q)
              Z.OneBody(q, p) += hZ * zij / j2hat2;
          }
        }
      }
      // std::cout << "diagram IIa and IIc " << Z.OneBodyNorm() << std::endl;
      // Z.EraseOneBody();

      // ###########################################################
      //  diagram II_b
      //
      // The two body operator
      //  Chi_222_b :
      //        c  | eta |  p
      //           |_____|
      //        b  |_____|  e
      //           | eta |
      //        a  |     |  d
      // ###########################################################

      //************************************************
      //  THIS PART CAN BE FURTHER POLISHED
      //************************************************

      TwoBodyME Chi_222_b = Eta.TwoBody;
      Chi_222_b.Erase();
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();

        for (int ibra = 0; ibra < nKets; ++ibra)
        {
          Ket bra = tbc.GetKet(ibra);
          int c = bra.p;
          int p = bra.q;

          Orbit op = Z.modelspace->GetOrbit(p);
          double jp = op.j2 / 2.;
          double n_p = op.occ;
          double nbar_p = 1.0 - n_p;

          Orbit oc = Z.modelspace->GetOrbit(c);
          double jc = oc.j2 / 2.;
          double n_c = oc.occ;
          double nbar_c = 1.0 - n_c;

          for (int jbra = 0; jbra < nKets; ++jbra)
          {
            Ket bra_j = tbc.GetKet(jbra);
            int a = bra_j.p;
            int d = bra_j.q;

            Orbit oa = Z.modelspace->GetOrbit(a);
            double ja = oa.j2 / 2.;
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;

            Orbit od = Z.modelspace->GetOrbit(d);
            double jd = od.j2 / 2.;
            double n_d = od.occ;
            double nbar_d = 1.0 - n_d;

            //---------------------------------
            double zij = 0;
            double zijtest = 0.;

            for (int kbra = 0; kbra < nKets; ++kbra)
            {
              Ket bra_k = tbc.GetKet(kbra);
              int b = bra_k.p;
              int e = bra_k.q;

              Orbit ob = Z.modelspace->GetOrbit(b);
              double jb = ob.j2 / 2.;
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;

              Orbit oe = Z.modelspace->GetOrbit(e);
              double je = oe.j2 / 2.;
              double n_e = oe.occ;
              double nbar_e = 1.0 - n_e;
              double occfactor = (nbar_a * nbar_d * n_b * n_e - nbar_b * nbar_e * n_a * n_d);
              if (std::abs(occfactor) < 1e-6)
                continue;

              double eta1 = Eta.TwoBody.GetTBME_J(J0, J0, c, p, b, e);
              double eta2 = Eta.TwoBody.GetTBME_J(J0, J0, b, e, a, d);

              zij += (2 * J0 + 1) * occfactor * eta1 * eta2;
              if (e != b)
              {
                zij += (2 * J0 + 1) * occfactor * eta1 * eta2;
              }
            }

            Chi_222_b.GetMatrix(ch, ch)(ibra, jbra) += zij;
          }
        }
        //--------------------------------------------------
      } // for p


      // ###########################################################
      //  diagram II_b
      // IIb_pq = 1/4 1/(2 jp + 1) \sum_acdJ0 Chi_222_b_cpad * Gamma_bar_adcq

#pragma omp parallel for
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
          double zij = 0;

          // loop abcde
          for (auto &c : Z.modelspace->all_orbits)
          {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double jc = oc.j2 / 2.;

            int J0min = std::abs(oc.j2 - op.j2) / 2;
            int J0max = (oc.j2 + op.j2) / 2;

            for (int J0 = J0min; J0 <= J0max; J0++)
            {
              for (auto &a : Z.modelspace->all_orbits)
              {
                Orbit &oa = Z.modelspace->GetOrbit(a);
                double ja = oa.j2 / 2.;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double jd = od.j2 / 2.;

                  zij += Chi_222_b.GetTBME_J_norm(J0, J0, c, p, a, d) * Gamma.TwoBody.GetTBME_J(J0, J0, a, d, c, q);
                } // J0
              }
            }
          }
          Z.OneBody(p, q) += 0.25 * zij / (op.j2 + 1.0);
          if (p != q)
            Z.OneBody(q, p) += 0.25 * hZ * zij / (op.j2 + 1.0);
          //--------------------------------------------------
        } // for q
      }   // for p
          // std::cout<< "diagram IIb " << Z.OneBodyNorm() << std::endl;
          // Z.EraseOneBody();

      // ###########################################################
      //  diagram II_d
      //
      // IId_pq = - 1/4 1/(2 jp + 1) \sum_abeJ0  Chi_222_a_bqae * Gamma_bar_bqae
      // ###########################################################
#pragma omp parallel for
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


          for (int ch = 0; ch < nch; ++ch)
          {
            TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
            int J0 = tbc.J;
            int nKets = tbc.GetNumberKets();
            for (int jbra = 0; jbra < nKets; ++jbra)
            {
              Ket &bra_j = tbc.GetKet(jbra);
              int a = bra_j.p;
              int e = bra_j.q;


              for (auto &b : Z.modelspace->all_orbits)
              {
                double MEs = Chi_222_b.GetTBME_J_norm(J0, J0, b, q, a, e) * Gamma.TwoBody.GetTBME_J(J0, J0, b, p, a, e);
                zij -= MEs;
                if (a != e)
                {
                  zij -= MEs;
                }
              }
            }
          }


          Z.OneBody(p, q) -= 0.25 * zij / (op.j2 + 1.0);
          if (p != q)
            Z.OneBody(q, p) -= 0.25 * hZ * zij / (op.j2 + 1.0);
          //--------------------------------------------------

        } // for q
      }   // for p
          // std::cout<< "diagram IId " << Z.OneBodyNorm() << std::endl;
          // Z.EraseOneBody();
    }

    // *********************************************************************************** //
    //                                  Diagram III                                        //
    // *********************************************************************************** //
    if (use_TypeIII_1b)
    {
      // ###########################################################
      //  diagram III_a and diagram III_b
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

        for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
        {
          Orbit &oe = Z.modelspace->GetOrbit(e);
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
      }   // d

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
      }   // for p
          // std::cout<< "diagram IIIa and IIIb " << Z.OneBodyNorm() << std::endl;
          // Z.EraseOneBody();
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_internal;
    return;
  }





////////////////////////////////////////////////////////////////////////////////////

  void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {

    if( use_1b_intermediates)
    {
      comm223_232_chi1b( Eta, Gamma, Z);  // topology with 1-body intermediate (fast)
    }
    if( use_2b_intermediates)
    {
      comm223_232_chi2b( Eta, Gamma, Z);  // topology with 2-body intermediate (slow)
    }

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
    int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC(); // number of CC channels
    auto &Z2 = Z.TwoBody;

    // determine symmetry
    int hEta = Eta.IsHermitian() ? 1 : -1;
    int hGamma = Gamma.IsHermitian() ? 1 : -1;
    // int hZ = Z.IsHermitian() ? 1 : -1;
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

/// Build the intermediate one-body operators
#pragma omp parallel for schedule(dynamic)
      for (size_t p = 0; p < norbits; p++)
      {
        Orbit &op = Z.modelspace->GetOrbit(p);
        for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
        {
          Orbit &oq = Z.modelspace->GetOrbit(q);

          double chi_pq = 0;
          double chiY_pq = 0;

          for (auto a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            if ( nbar_a<1e-6) continue;

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
                  double xapij = Eta.TwoBody.GetTBME_J(J2, J2, a, p, i, j);
                  double yapij = Gamma.TwoBody.GetTBME_J(J2, J2, a, p, i, j);

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

                if ( std::abs(occfactor)<1e-7) continue;

                int J2min = std::max({std::abs(oa.j2 - ob.j2), std::abs(oi.j2 - oq.j2), std::abs(oi.j2 - op.j2)}) / 2;
                int J2max = std::min({oa.j2 + ob.j2, oi.j2 + oq.j2, oi.j2 + op.j2}) / 2;

                for (int J2 = J2min; J2 <= J2max; J2++)
                {
                  double xipab = Eta.TwoBody.GetTBME_J(J2, J2, i, p, a, b);
                  double xabiq = Eta.TwoBody.GetTBME_J(J2, J2, a, b, i, q);
                  double yipab = Gamma.TwoBody.GetTBME_J(J2, J2, i, p, a, b);

                  chi_pq += 0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * xipab * xabiq;
                  chiY_pq += 0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * yipab * xabiq;
                }
              } // for b

            } // for i
          }   // for a
          CHI_I(p, q) = chi_pq ;
          CHI_II(p, q) = chiY_pq  ;
        } // for q
      }   // for p

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
          size_t p = bra.p;
          size_t q = bra.q;
          Orbit &op = Z.modelspace->GetOrbit(p);
          Orbit &oq = Z.modelspace->GetOrbit(q);
          int phasepq = bra.Phase(J);
          for (size_t iket = ibra; iket < nkets; iket++)
          {
            Ket &ket = tbc_ket.GetKet(iket);
            size_t r = ket.p;
            size_t s = ket.q;
            Orbit &oR = Z.modelspace->GetOrbit(r);
            Orbit &os = Z.modelspace->GetOrbit(s);
            double zpqrs = 0;

            for (auto b : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
            {
              zpqrs += CHI_I(p, b) * Gamma.TwoBody.GetTBME_J(J, J, b, q, r, s);
              zpqrs += hZ * CHI_II(b, p) * Eta.TwoBody.GetTBME_J(J, J, b, q, r, s); // tricky minus sign.
            }                                                                       // for a
            for (auto b : Z.OneBodyChannels.at({oq.l, oq.j2, oq.tz2}))
            {
              zpqrs += CHI_I(q, b) * Gamma.TwoBody.GetTBME_J(J, J, p, b, r, s);
              zpqrs += hZ * CHI_II(b, q) * Eta.TwoBody.GetTBME_J(J, J, p, b, r, s); // tricky minus sign.
            }                                                                       // for a
            for (auto b : Z.OneBodyChannels.at({oR.l, oR.j2, oR.tz2}))
            {
              zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_I(b, r);
              zpqrs -= Eta.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_II(b, r);
            } // for a
            for (auto b : Z.OneBodyChannels.at({os.l, os.j2, os.tz2}))
            {
              zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_I(b, s);
              zpqrs -= Eta.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_II(b, s);
            } // for a

            // normalize
            if (p == q)
              zpqrs /= PhysConst::SQRT2;
            if (r == s)
              zpqrs /= PhysConst::SQRT2;
            Z2.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
          } // for iket
        }   // for ibra

      } // for itmat

      CHI_I.clear();
      CHI_II.clear();
    

    // Timer
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return;
  }// comm223_232_chi1b





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
    // int nch = Z.modelspace->GetNumberTwoBodyChannels(); // number of TB channels
    int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC(); // number of CC channels
    auto &Z2 = Z.TwoBody;

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
    std::deque<arma::mat> bar_Eta(n_nonzero);      // released
    std::deque<arma::mat> nnnbar_Eta(n_nonzero);   // released
    std::deque<arma::mat> nnnbar_Eta_d(n_nonzero); // released
    std::deque<arma::mat> bar_Gamma(n_nonzero);    // released
    /// initial bar_Eta and nnnbar_Eta
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_Eta[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      nnnbar_Eta[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      nnnbar_Eta_d[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      bar_Gamma[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    }

    /// Pandya transformation
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
          double jc = oc.j2 * 0.5;
          double n_c = oc.occ;
          double nbar_c = 1.0 - n_c;

          Orbit &od = Z.modelspace->GetOrbit(d);
          double jd = od.j2 * 0.5;
          double n_d = od.occ;
          double nbar_d = 1.0 - n_d;

          double occfactor_c = (nbar_c * nbar_b * n_a + nbar_a * n_c * n_b);
          double occfactor_d = (nbar_d * nbar_a * n_b + nbar_b * n_a * n_d);
          int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
          int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
          double Xbar = 0;
          double Ybar = 0;
          double Zbar = 0;
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
              double temp_eta = Eta.TwoBody.GetTBME_J(J_std, a, d, c, b);
              Xbar -=               (2 * J_std + 1) * sixj1 * temp_eta;
              Ybar -= occfactor_c * (2 * J_std + 1) * sixj1 * temp_eta;
              Zbar -= occfactor_d * (2 * J_std + 1) * sixj1 * temp_eta;
              Gammabar -= (2 * J_std + 1) * sixj1 * Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
            }
          }
          bar_Eta[ch_cc](ibra_cc, iket_cc) = Xbar;
          nnnbar_Eta[ch_cc](ibra_cc, iket_cc) = Ybar;
          nnnbar_Eta_d[ch_cc](ibra_cc, iket_cc) = Zbar;
          bar_Gamma[ch_cc](ibra_cc, iket_cc) = Gammabar;
        }
        //-------------------
      }
    }

    Z.profiler.timer[std::string(__func__) + " Global array"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    //-------------------------------------------------------------------------------
    std::deque<arma::mat> bar_CHI_VI(n_nonzero); //  released
    /// build intermediate bar operator
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_CHI_VI[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    }

    // BUILD bar_CHI_VI
    // bar_CHI_VI = (bar{na} nb nc + bar{nb} bar{nc} na) \bar{Gamma}_dhba \bar{Eta}_baqc
#pragma omp parallel for
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      bar_CHI_VI[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta_d[ch_cc];
    }

    // release memory
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      nnnbar_Eta_d[ch_cc].clear();
    }
    nnnbar_Eta_d.clear();

    // BUILD bar_CHI_III
    std::deque<arma::mat> barCHI_III(n_nonzero); //  released
    /// build intermediate bar operator
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      barCHI_III[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    }

// matrix multiplication
#pragma omp parallel for
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      barCHI_III[ch_cc] = bar_Eta[ch_cc] * nnnbar_Eta[ch_cc];
    }

    // build Chi_III
    std::deque<arma::mat> Chi_III(nch); // released
    std::deque<arma::mat> CHI_VI(nch);  // released
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      Chi_III[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      CHI_VI[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

    // Todo : Combine the following two blocks
    //
    // Inverse Pandya transformation
    //  X^J_ijkl  = - ( 1- P_ij )  sum_J' (2J'+1)  { i j J }  \bar{X}^J'_il`kj`
    //                                             { k l J'}
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
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
          Chi_III[ch](ibra, iket) += zijkl;
        }
      }
    }

    // BUILD CHI_VI
    // Inverse Pandya transformation
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets * 2; ++ibra)
      {
        size_t i, j;
        if (ibra < nKets)
        {
          Ket &bra = tbc.GetKet(ibra);
          i = bra.p;
          j = bra.q;
        }
        else
        {
          Ket &bra = tbc.GetKet(ibra - nKets);
          i = bra.q;
          j = bra.p;
        }
        if (ibra >= nKets and i == j)
          continue;

        Orbit &oi = Z.modelspace->GetOrbit(i);
        int ji = oi.j2;
        Orbit &oj = Z.modelspace->GetOrbit(j);
        int jj = oj.j2;

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
          double commijkl = 0;
          double commijlk = 0;

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
            int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
            int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
            if (indx_il < 0 or indx_kj < 0)
              continue;

            indx_il += (i > l ? nkets_cc : 0);
            indx_kj += (k > j ? nkets_cc : 0);
            double me1 = bar_CHI_VI[ch_cc](indx_il, indx_kj);
            commijkl -= (2 * Jprime + 1) * sixj * me1;
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

            // exchange k and l
            indx_ik += (i > k ? nkets_cc : 0);
            indx_lj += (l > j ? nkets_cc : 0);
            double me2 = bar_CHI_VI[ch_cc](indx_ik, indx_lj);
            commijlk -= (2 * Jprime + 1) * sixj * me2;
          }

          double zijkl = (commijkl - Z.modelspace->phase((jk + jl) / 2 - J0) * commijlk);

          CHI_VI[ch](ibra, iket) += zijkl;
        }
      }
    }
    // Todo End ###########################################


    // release memory
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      bar_CHI_VI[ch_cc].clear();
    }
    bar_CHI_VI.clear();

    std::deque<arma::mat> Eta_matrix(nch);   // released
    std::deque<arma::mat> Gamma_matrix(nch); // released
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      Eta_matrix[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      Gamma_matrix[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

    // full matrix
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
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
          Eta_matrix[ch](ibra, iket) = EtaME;
          Gamma_matrix[ch](ibra, iket) = GammaME;
          if (i != j)
          {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Eta_matrix[ch](ibra + nKets, iket) = phase * EtaME;
            Gamma_matrix[ch](ibra + nKets, iket) = phase * GammaME;
            if (k != l)
            {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Eta_matrix[ch](ibra + nKets, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra + nKets, iket + nKets) = phase * GammaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix[ch](ibra, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra, iket + nKets) = phase * GammaME;
            }
          }
          else
          {
            if (k != l)
            {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix[ch](ibra, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra, iket + nKets) = phase * GammaME;
            }
          }
        }
      }
    }

// Diagram IIIc and Diagram IIId
// IIa_pgqh = \sum_ad Chi_III^J0_pgda * Gamma^J0_daqh
// IIc_pgqh = \sum_ad Gamma^J0_pgad * Chi_III^J0_adqh
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      arma::Mat<double> Multi_matirx(nKets * 2, nKets * 2);
      Multi_matirx = Chi_III[ch] * Gamma_matrix[ch];
      Multi_matirx -= (Eta_matrix[ch] * CHI_VI[ch]);
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          double z_IIa = Multi_matirx(ibra, iket);
          double z_IIc = hZ * Multi_matirx(iket, ibra);

          if (p == g)
          {
            z_IIa /= PhysConst::SQRT2;
            z_IIc /= PhysConst::SQRT2;
          }
          if (q == h)
          {
            z_IIa /= PhysConst::SQRT2;
            z_IIc /= PhysConst::SQRT2;
          }
          Z2.AddToTBME(ch, ch, ibra, iket, z_IIa + z_IIc);
        } // iket
      }   // ibra
    }     // J0 channel


    // release memory
    for (int ch = 0; ch < nch; ++ch)
    {
      Chi_III[ch].clear();
      CHI_VI[ch].clear();
    }
    Chi_III.clear();
    CHI_VI.clear();

    //____________________ bar_CHI_V
    std::deque<arma::mat> bar_CHI_V(n_nonzero); // released
    /// initial bar_CHI_V
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_CHI_V[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    }

    // build bar_CHI_V
#pragma omp parallel for
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      bar_CHI_V[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta[ch_cc];
    }
    // release memory
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      nnnbar_Eta[ch_cc].clear();
    }
    nnnbar_Eta.clear();

    // ****************************************
    // eta_c and eta_d
    std::deque<arma::mat> Eta_matrix_c(nch); // released
    std::deque<arma::mat> Eta_matrix_d(nch); // released
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      Eta_matrix_c[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      Eta_matrix_d[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

    // full matrix
  #pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
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
          Eta_matrix_c[ch](ibra, iket) = occfactor_k * EtaME;
          Eta_matrix_d[ch](ibra, iket) = occfactor_l * EtaME;
          if (i != j)
          {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Eta_matrix_c[ch](ibra + nKets, iket) = occfactor_k * phase * EtaME;
            Eta_matrix_d[ch](ibra + nKets, iket) = occfactor_l * phase * EtaME;
            if (k != l)
            {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Eta_matrix_c[ch](ibra + nKets, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra + nKets, iket + nKets) = occfactor_k * phase * EtaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix_c[ch](ibra, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
          else
          {
            if (k != l)
            {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix_c[ch](ibra, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
        }
      }
    }
    // eta_c and eta_d
    // ****************************************

    std::deque<arma::mat> CHI_VII(nch); // released
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      // Not symmetric
      CHI_VII[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      CHI_VII[ch] = Gamma_matrix[ch] * Eta_matrix_d[ch];
    }
    // release memory
    for (int ch = 0; ch < nch; ++ch)
    {
      Gamma_matrix[ch].clear();
    }
    Gamma_matrix.clear();

    // declare CHI_IV
    std::deque<arma::mat> CHI_IV(nch);  // released
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      // Not symmetric
      CHI_IV[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      CHI_IV[ch] = Eta_matrix[ch] * Eta_matrix_c[ch];
      CHI_IV[ch] += (Eta_matrix[ch] * Eta_matrix_d[ch]).t();
    }

    // release memory
    for (int ch = 0; ch < nch; ++ch)
    {
      Eta_matrix[ch].clear();
      Eta_matrix_c[ch].clear();
      Eta_matrix_d[ch].clear();
    }
    Eta_matrix.clear();
    Eta_matrix_c.clear();
    Eta_matrix_d.clear();


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

    /// Pandya transformation only recouple the angula momentum
    /// IIb and IId
    /// diagram IIIa - diagram IIIb
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
          double XbarIIef = 0;
          double XbarIIIab = 0;
          double XbarIIIef = 0;
          for (int J_std = jmin; J_std <= jmax; J_std++)
          {
            int phaseFactor = Z.modelspace->phase(J_std + (oc.j2 + ob.j2) / 2);
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
                if (a > d)
                  indx_ad += nkets;
                if (b > c)
                  indx_bc += nkets;
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

    /// release memory
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      barCHI_III[ch_cc].clear();
      bar_CHI_V[ch_cc].clear();
    }
    barCHI_III.clear();
    bar_CHI_V.clear();

    //_______________________________________________________________________________
    std::deque<arma::mat> bar_CHI_IV(n_nonzero);    // released
    std::deque<arma::mat> bar_CHI_VII_CC(n_nonzero);// released
    /// build intermediate bar operator
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_CHI_IV[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      bar_CHI_VII_CC[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    }

    /// Pandya transformation only recouple the angula momentum
    /// IIe and IIf
    /// diagram IIIe and IIIf
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
          double XbarIIef = 0;
          double XbarIIIab = 0;
          double XbarIIIef = 0;
          for (int J_std = jmin; J_std <= jmax; J_std++)
          {
            int phaseFactor = Z.modelspace->phase(J_std + (oc.j2 + ob.j2) / 2);
            double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
            if (std::abs(sixj1) > 1e-8)
            {
              /// IIIef
              int parity_cc = (oa.l + od.l) % 2;
              int Tz_J2 = (ob.tz2 + oc.tz2) / 2;
              int ch_J2 = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_cc, Tz_J2);
              TwoBodyChannel &tbc_J2 = Z.modelspace->GetTwoBodyChannel(ch_J2);
              int nkets = tbc_J2.GetNumberKets();
              int indx_ad = tbc_J2.GetLocalIndex(std::min(int(a), int(d)), std::max(int(a), int(d)));
              int indx_bc = tbc_J2.GetLocalIndex(std::min(int(b), int(c)), std::max(int(b), int(c)));
              if (indx_ad >= 0 and indx_bc >= 0)
              {
                int indx_cb = indx_bc;
                if (a > d)
                  indx_ad += nkets;
                if (b > c)
                  indx_bc += nkets;
                if (c > b)
                  indx_cb += nkets;

                XbarIIef -= (2 * J_std + 1) * sixj1 * CHI_IV[ch_J2](indx_ad, indx_cb);
                XbarIIIef += phaseFactor * (2 * J_std + 1) * sixj1 * CHI_VII[ch_J2](indx_ad, indx_bc);
              }
            }
          }
          bar_CHI_IV[ch_cc](ibra_cc, iket_cc) = XbarIIef;
          bar_CHI_VII_CC[ch_cc](ibra_cc, iket_cc) = XbarIIIef;
        }

        //-------------------
      }
    }

    /// release memory    
    for (int ch = 0; ch < nch; ++ch)
    {
      CHI_IV[ch].clear();
      CHI_VII[ch].clear();
    }
    CHI_IV.clear();
    CHI_VII.clear();

    std::deque<arma::mat> CHI_III_final(nch);  // released
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      // Not symmetric
      CHI_III_final[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

#pragma omp parallel for
    for (int ch = 0; ch < n_nonzero; ++ch)
    {
      CHI_III_final[ch] = bar_Gamma[ch] * barCHI_III_RC[ch];
    }
    /// release memory
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      barCHI_III_RC[ch_cc].clear();
    }
    barCHI_III_RC.clear();

    //  Inverse Pandya transformation
    //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
    //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
    //                                    { k l J'}
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
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
        int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          size_t k, l;
          Ket &ket = tbc.GetKet(iket);
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
            int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
            int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
            if (indx_jl < 0 or indx_ik < 0)
              continue;

            int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
            indx_jl += (j > l ? nkets_cc : 0);
            int indx_ki = indx_ik + (k > i ? nkets_cc : 0);
            indx_ik += (i > k ? nkets_cc : 0);
            int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
            int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);

            // IIb and IId
            // direct term
            double meIIb = CHI_III_final[ch_cc](indx_jl, indx_ik);
            commjikl -= phase1 * (2 * Jprime + 1) * sixj1 * meIIb;
            // exchange ij and kl
            double meIId = CHI_III_final[ch_cc](indx_ik, indx_jl);
            commijlk -= phase2 * (2 * Jprime + 1) * sixj1 * meIId;


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
            int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
            int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
            if (indx_il < 0 or indx_jk < 0)
              continue;

            int indx_kj = indx_jk + (k > j ? nkets_cc : 0);
            indx_jk += (j > k ? nkets_cc : 0);
            int indx_li = indx_il + (l > i ? nkets_cc : 0);
            indx_il += (i > l ? nkets_cc : 0);
            int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
            int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);

            // IIb and IId
            // exchange k and l
            double meIIb = CHI_III_final[ch_cc](indx_jk, indx_il);
            commjilk -= phase1 * (2 * Jprime + 1) * sixj1 * meIIb;
            // exchange i and j
            double meIId = CHI_III_final[ch_cc](indx_il, indx_jk);
            commijkl -= phase2 * (2 * Jprime + 1) * sixj1 * meIId;

          }

          double zijkl = (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
          zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, phaseFactor * zijkl);
        }
      }
    }
    
    for (int ch = 0; ch < nch; ++ch)
    {
      CHI_III_final[ch].clear();
    }
    CHI_III_final.clear();

    // ##########################################
    std::deque<arma::mat> CHI_V_final(nch);    // released
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      // Not symmetric
      CHI_V_final[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

#pragma omp parallel for
    for (int ch = 0; ch < n_nonzero; ++ch)
    {
      CHI_V_final[ch] = bar_Eta[ch] * bar_CHI_V_RC[ch];
    }
    /// release memory
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      bar_CHI_V_RC[ch_cc].clear();
    }
    bar_CHI_V_RC.clear();

    //  Inverse Pandya transformation
    //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
    //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
    //                                    { k l J'}
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
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
        int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          size_t k, l;
          Ket &ket = tbc.GetKet(iket);
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
            int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
            int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
            if (indx_jl < 0 or indx_ik < 0)
              continue;

            int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
            indx_jl += (j > l ? nkets_cc : 0);
            int indx_ki = indx_ik + (k > i ? nkets_cc : 0);
            indx_ik += (i > k ? nkets_cc : 0);
            int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
            int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);

            // IIIa and IIIb
            // direct term
            double meIIIa = CHI_V_final[ch_cc](indx_jl, indx_ik);
            commjikl -= phase1 * (2 * Jprime + 1) * sixj1 * meIIIa;
            // exchange ij and kl
            double meIIIb = CHI_V_final[ch_cc](indx_ik, indx_jl);
            commijlk -= phase2 * (2 * Jprime + 1) * sixj1 * meIIIb;
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
            int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
            int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
            if (indx_il < 0 or indx_jk < 0)
              continue;

            int indx_kj = indx_jk + (k > j ? nkets_cc : 0);
            indx_jk += (j > k ? nkets_cc : 0);
            int indx_li = indx_il + (l > i ? nkets_cc : 0);
            indx_il += (i > l ? nkets_cc : 0);
            int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
            int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);

            /// IIIa and IIIb
            // exchange k and l
            double meIIIa = CHI_V_final[ch_cc](indx_jk, indx_il);
            commjilk -= phase1 * (2 * Jprime + 1) * sixj1 * meIIIa;
            // exchange i and j
            double meIIIb = CHI_V_final[ch_cc](indx_il, indx_jk);
            commijkl -= phase2 * (2 * Jprime + 1) * sixj1 * meIIIb;
          }

          double zijkl = (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
          zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, phaseFactor * zijkl);
        }
      }
    }

    for (int ch = 0; ch < nch; ++ch)
    {
      CHI_V_final[ch].clear();
    }
    CHI_V_final.clear();

    // ######################################

    std::deque<arma::mat> bar_CHI_gamma(n_nonzero);   // released
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
      bar_Gamma[ch_cc].clear();
    }
    bar_CHI_IV.clear();
    bar_Gamma.clear();


    //  Inverse Pandya transformation
    //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
    //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
    //                                    { k l J'}
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
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
        int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          size_t k, l;
          Ket &ket = tbc.GetKet(iket);
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
            int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
            int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
            if (indx_jl < 0 or indx_ik < 0)
              continue;


            int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
            indx_jl += (j > l ? nkets_cc : 0);
            int indx_ki = indx_ik + (k > i ? nkets_cc : 0);
            indx_ik += (i > k ? nkets_cc : 0);
            int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
            int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);


            // IIe and IIf
            // jilk, exchange i and j, k and l     ->  jk   li
            double meIIef = bar_CHI_gamma[ch_cc](indx_ik, indx_lj);
            commjilk -= 0.5 * phaseFactor * (2 * Jprime + 1) * sixj1 * meIIef;
            // jikl, exchange i and j    ->  jl ki
            meIIef = bar_CHI_gamma[ch_cc](indx_jl, indx_ki);
            commijkl -= 0.5 * phaseFactor * (2 * Jprime + 1) * sixj1 * meIIef;
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
            int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
            int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
            if (indx_il < 0 or indx_jk < 0)
              continue;

            int indx_kj = indx_jk + (k > j ? nkets_cc : 0);
            indx_jk += (j > k ? nkets_cc : 0);
            int indx_li = indx_il + (l > i ? nkets_cc : 0);
            indx_il += (i > l ? nkets_cc : 0);
            int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
            int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);

            // IIe and IIf
            // ijlk, exchange k and l     ->  ik lj
            double meIIef = bar_CHI_gamma[ch_cc](indx_jk, indx_li);
            commijlk -= 0.5 * phaseFactor * (2 * Jprime + 1) * sixj1 * meIIef;
            // jikl, exchange i and j    ->  jl ki
            meIIef = bar_CHI_gamma[ch_cc](indx_il, indx_kj);
            commjikl -= 0.5 * phaseFactor * (2 * Jprime + 1) * sixj1 * meIIef;
          }

          double zijkl = (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
          zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, phaseFactor * zijkl);
        }
      }
    }

    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)

    {
      bar_CHI_gamma[ch_cc].clear();
    }
    bar_CHI_gamma.clear();
    // ########################################################################



    std::deque<arma::mat> bar_CHI_VII_CC_ef(n_nonzero);    // released
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
      arma::Mat<double> Multi_matirx(nKets * 2, nKets * 2);
      Multi_matirx = bar_CHI_VII_CC[ch_cc] * bar_Eta[ch_cc];
      bar_CHI_VII_CC_ef[ch_cc] = Multi_matirx + hZ * Multi_matirx.t();
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
        int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          size_t k, l;
          Ket &ket = tbc.GetKet(iket);
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
            int indx_jl = tbc_cc.GetLocalIndex(std::min(j, l), std::max(j, l));
            int indx_ik = tbc_cc.GetLocalIndex(std::min(k, i), std::max(k, i));
            if (indx_jl < 0 or indx_ik < 0)
              continue;

            int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
            indx_jl += (j > l ? nkets_cc : 0);
            int indx_ki = indx_ik + (k > i ? nkets_cc : 0);
            indx_ik += (i > k ? nkets_cc : 0);
            int phase1 = Z.modelspace->phase(Jprime + (ji + jk) / 2);
            int phase2 = Z.modelspace->phase(Jprime + (jj + jl) / 2);

            /// IIIe and IIIf
            // direct term
            double meIIIef = bar_CHI_VII_CC_ef[ch_cc](indx_jl, indx_ki);
            commjikl -= 0.5 * (2 * Jprime + 1) * sixj1 * meIIIef;
            // exchange ij and kl
            double meIIIef2 = bar_CHI_VII_CC_ef[ch_cc](indx_ik, indx_lj);
            commijlk -= 0.5 * (2 * Jprime + 1) * sixj1 * meIIIef2;
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
            int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
            int indx_jk = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
            if (indx_il < 0 or indx_jk < 0)
              continue;

            int indx_kj = indx_jk + (k > j ? nkets_cc : 0);
            indx_jk += (j > k ? nkets_cc : 0);
            int indx_li = indx_il + (l > i ? nkets_cc : 0);
            indx_il += (i > l ? nkets_cc : 0);
            int phase1 = Z.modelspace->phase(Jprime + (ji + jl) / 2);
            int phase2 = Z.modelspace->phase(Jprime + (jj + jk) / 2);

            /// IIIe and IIIf
            // exchange k and l
            double meIIIef = bar_CHI_VII_CC_ef[ch_cc](indx_jk, indx_li);
            commjilk -= 0.5 * (2 * Jprime + 1) * sixj1 * meIIIef;
            // exchange i and j
            double meIIIef2 = bar_CHI_VII_CC_ef[ch_cc](indx_il, indx_kj);
            commijkl -= 0.5 * (2 * Jprime + 1) * sixj1 * meIIIef2;


          }

          double zijkl = (commjikl - Z.modelspace->phase((ji + jj) / 2 - J0) * commijkl);
          zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commjilk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commijlk);

          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, phaseFactor * zijkl);
        }
      }
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










  void comm223_232_slow(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {
    // global variables
    double t_start = omp_get_wtime();
    double t_internal = omp_get_wtime();
    double t_type = omp_get_wtime();
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
    int n_nonzero = Z.modelspace->GetNumberTwoBodyChannels_CC(); // number of CC channels
    auto &Z2 = Z.TwoBody;

    // determine symmetry
    int hEta = Eta.IsHermitian() ? 1 : -1;
    int hGamma = Gamma.IsHermitian() ? 1 : -1;
    int hZ = hGamma;
    // ####################################################################################
    //                      Factorization of Ia, Ib, IVa and IVb
    // ####################################################################################

    if (use_goose_tank_2b)
    {
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
#pragma omp parallel for schedule(dynamic)
      for (size_t p = 0; p < norbits; p++)
      {
        Orbit &op = Z.modelspace->GetOrbit(p);
        for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
        {
          Orbit &oq = Z.modelspace->GetOrbit(q);

          double chi_pq = 0;
          double chiY_pq = 0;

          for (auto a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;

            for (auto i : Z.modelspace->holes)
            {
              Orbit &oi = Z.modelspace->GetOrbit(i);
              double n_i = oi.occ;

              for (auto j : Z.modelspace->holes)
              {
                Orbit &oj = Z.modelspace->GetOrbit(j);
                double n_j = oj.occ;

                double occfactor = nbar_a * n_i * n_j;
                int J2min = std::max(std::abs(oa.j2 - oq.j2), std::abs(oi.j2 - oj.j2)) / 2;
                int J2max = std::min(oa.j2 + oq.j2, oi.j2 + oj.j2) / 2;

                for (int J2 = J2min; J2 <= J2max; J2++)
                {
                  double xijaq = Eta.TwoBody.GetTBME_J(J2, J2, i, j, a, q);
                  double xapij = Eta.TwoBody.GetTBME_J(J2, J2, a, p, i, j);
                  double yapij = Gamma.TwoBody.GetTBME_J(J2, J2, a, p, i, j);

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

                int J2min = std::max({std::abs(oa.j2 - ob.j2), std::abs(oi.j2 - oq.j2), std::abs(oi.j2 - op.j2)}) / 2;
                int J2max = std::min({oa.j2 + ob.j2, oi.j2 + oq.j2, oi.j2 + op.j2}) / 2;

                for (int J2 = J2min; J2 <= J2max; J2++)
                {
                  double xipab = Eta.TwoBody.GetTBME_J(J2, J2, i, p, a, b);
                  double xabiq = Eta.TwoBody.GetTBME_J(J2, J2, a, b, i, q);
                  double yipab = Gamma.TwoBody.GetTBME_J(J2, J2, i, p, a, b);

                  chi_pq += 0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * xipab * xabiq;
                  chiY_pq += 0.5 * occfactor * (2 * J2 + 1) / (oq.j2 + 1) * yipab * xabiq;
                }
              } // for b

            } // for i
          }   // for a
          CHI_I(p, q) = chi_pq;
          CHI_II(p, q) = chiY_pq;
        } // for q
      }   // for p

#pragma omp parallel for schedule(dynamic, 1)
      for (int ich = 0; ich < nch; ich++)
      {
        size_t ch_bra = ch_bra_list[ich];
        size_t ch_ket = ch_ket_list[ich];
        //     size_t ch_bra = itmat.first[0];
        //     size_t ch_ket = itmat.first[1];
        TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int J = tbc_bra.J;
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();
        for (size_t ibra = 0; ibra < nbras; ibra++)
        {
          Ket &bra = tbc_bra.GetKet(ibra);
          size_t p = bra.p;
          size_t q = bra.q;
          Orbit &op = Z.modelspace->GetOrbit(p);
          Orbit &oq = Z.modelspace->GetOrbit(q);
          int phasepq = bra.Phase(J);
          //       for ( size_t iket=0; iket<nkets; iket++)
          for (size_t iket = ibra; iket < nkets; iket++)
          {
            Ket &ket = tbc_ket.GetKet(iket);
            size_t r = ket.p;
            size_t s = ket.q;
            Orbit &oR = Z.modelspace->GetOrbit(r);
            Orbit &os = Z.modelspace->GetOrbit(s);
            double zpqrs = 0;

            for (auto b : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
            {
              if (use_GT_TypeI_2b)
                zpqrs += CHI_I(p, b) * Gamma.TwoBody.GetTBME_J(J, J, b, q, r, s);
              if (use_GT_TypeIV_2b)
                zpqrs += hZ * CHI_II(b, p) * Eta.TwoBody.GetTBME_J(J, J, b, q, r, s); // tricky minus sign.
            }                                                                         // for a
            for (auto b : Z.OneBodyChannels.at({oq.l, oq.j2, oq.tz2}))
            {
              if (use_GT_TypeI_2b)
                zpqrs += CHI_I(q, b) * Gamma.TwoBody.GetTBME_J(J, J, p, b, r, s);
              if (use_GT_TypeIV_2b)
                zpqrs += hZ * CHI_II(b, q) * Eta.TwoBody.GetTBME_J(J, J, p, b, r, s); // tricky minus sign.
            }                                                                         // for a
            for (auto b : Z.OneBodyChannels.at({oR.l, oR.j2, oR.tz2}))
            {
              if (use_GT_TypeI_2b)
                zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_I(b, r);
              if (use_GT_TypeIV_2b)
                zpqrs -= Eta.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_II(b, r);
            } // for a
            for (auto b : Z.OneBodyChannels.at({os.l, os.j2, os.tz2}))
            {
              if (use_GT_TypeI_2b)
                zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_I(b, s);
              if (use_GT_TypeIV_2b)
                zpqrs -= Eta.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_II(b, s);
            } // for a

            // normalize
            if (p == q)
              zpqrs /= PhysConst::SQRT2;
            if (r == s)
              zpqrs /= PhysConst::SQRT2;
            Z2.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
          } // for iket
        }   // for ibra


      } // for itmat

      Z.profiler.timer[std::string(__func__) + " Diagram I and Diagram IV"] += omp_get_wtime() - t_type;
      t_type = omp_get_wtime();

      Z.profiler.timer[std::string(__func__) + " Diagram I and Diagram IV"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();
    }
    if (use_goose_tank_only_2b)
    {
      return;
    }

    // std::cout << "diagram I and IV " << Z.TwoBodyNorm() << std::endl;

    //______________________________________________________________________
    // global array
    std::deque<arma::mat> bar_Eta(n_nonzero);
    std::deque<arma::mat> nnnbar_Eta(n_nonzero);
    std::deque<arma::mat> nnnbar_Eta_d(n_nonzero);
    std::deque<arma::mat> bar_Gamma(n_nonzero);
    /// initial bar_Eta and nnnbar_Eta
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_Eta[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      nnnbar_Eta[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      nnnbar_Eta_d[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      bar_Gamma[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
    }

    /// Pandya transformation
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
          double jc = oc.j2 * 0.5;
          double n_c = oc.occ;
          double nbar_c = 1.0 - n_c;

          Orbit &od = Z.modelspace->GetOrbit(d);
          double jd = od.j2 * 0.5;
          double n_d = od.occ;
          double nbar_d = 1.0 - n_d;

          double occfactor_c = (nbar_c * nbar_b * n_a + nbar_a * n_c * n_b);
          double occfactor_d = (nbar_d * nbar_a * n_b + nbar_b * n_a * n_d);
          int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
          int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
          double Xbar = 0;
          double Ybar = 0;
          double Zbar = 0;
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
              double temp_eta = Eta.TwoBody.GetTBME_J(J_std, a, d, c, b);
              Xbar -= (2 * J_std + 1) * sixj1 * temp_eta;
              Ybar -= occfactor_c * (2 * J_std + 1) * sixj1 * temp_eta;
              Zbar -= occfactor_d * (2 * J_std + 1) * sixj1 * temp_eta;
              Gammabar -= (2 * J_std + 1) * sixj1 * Gamma.TwoBody.GetTBME_J(J_std, a, d, c, b);
            }
          }
          bar_Eta[ch_cc](ibra_cc, iket_cc) = Xbar;
          nnnbar_Eta[ch_cc](ibra_cc, iket_cc) = Ybar;
          nnnbar_Eta_d[ch_cc](ibra_cc, iket_cc) = Zbar;
          bar_Gamma[ch_cc](ibra_cc, iket_cc) = Gammabar;
        }
        //-------------------
      }
    }

    std::deque<arma::mat> Eta_matrix(nch);
    std::deque<arma::mat> Eta_matrix_c(nch);
    std::deque<arma::mat> Eta_matrix_d(nch);
    std::deque<arma::mat> Gamma_matrix(nch);
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      Eta_matrix[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      Eta_matrix_c[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      Eta_matrix_d[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      Gamma_matrix[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

    // full matrix
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
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
          Eta_matrix[ch](ibra, iket) = EtaME;
          Gamma_matrix[ch](ibra, iket) = GammaME;
          Eta_matrix_c[ch](ibra, iket) = occfactor_k * EtaME;
          Eta_matrix_d[ch](ibra, iket) = occfactor_l * EtaME;
          if (i != j)
          {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Eta_matrix[ch](ibra + nKets, iket) = phase * EtaME;
            Gamma_matrix[ch](ibra + nKets, iket) = phase * GammaME;
            Eta_matrix_c[ch](ibra + nKets, iket) = occfactor_k * phase * EtaME;
            Eta_matrix_d[ch](ibra + nKets, iket) = occfactor_l * phase * EtaME;
            if (k != l)
            {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Eta_matrix[ch](ibra + nKets, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra + nKets, iket + nKets) = phase * GammaME;
              Eta_matrix_c[ch](ibra + nKets, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra + nKets, iket + nKets) = occfactor_k * phase * EtaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix[ch](ibra, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra, iket + nKets) = phase * GammaME;
              Eta_matrix_c[ch](ibra, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
          else
          {
            if (k != l)
            {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix[ch](ibra, iket + nKets) = phase * EtaME;
              Gamma_matrix[ch](ibra, iket + nKets) = phase * GammaME;
              Eta_matrix_c[ch](ibra, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra, iket + nKets) = occfactor_k * phase * EtaME;
            }
          }
        }
      }
    }

    Z.profiler.timer[std::string(__func__) + " Global array"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    // *********************************************************************************** //
    //                                  Diagram II                                         //
    // *********************************************************************************** //
    if (use_TypeII_2b)
    {
      // ####################################################################################
      //                      Factorization of IIa and IIc
      // ####################################################################################
      // Theintermediate two body operator
      //  Chi_III :
      //            eta |
      //           _____|
      //          /\    |
      //   |     (  )
      //   |_____ \/
      //   | eta
      //
      //  Pandya transform
      //  \bar{X}^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
      //                              { k l J'}
      //
      //  \bar{Chi}_III = \sum_{bc J2 J3} (nbar_b * nbar_d * n_c - nbar_c * n_b * n_d )
      //                  \bar{eta}_pacb * \bar{eta}_cbdg
      //-------------------------------------------------------------------------------

      std::deque<arma::mat> barCHI_III(n_nonzero);
      std::deque<arma::mat> barCHI_III_RC(n_nonzero); // Recoupled bar CHI_III
      /// build intermediate bar operator
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        barCHI_III[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        barCHI_III_RC[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

#pragma omp parallel for
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        barCHI_III[ch_cc] = bar_Eta[ch_cc] * nnnbar_Eta[ch_cc];
      }

      // build Chi_III
      std::deque<arma::mat> Chi_III(nch);
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int nKets = tbc.GetNumberKets();
        Chi_III[ch] = arma::mat(nKets, nKets * 2, arma::fill::zeros);
      }

      // Inverse Pandya transformation
      //  X^J_ijkl  = - ( 1- P_ij )  sum_J' (2J'+1)  { i j J }  \bar{X}^J'_il`kj`
      //                                             { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
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
            Chi_III[ch](ibra, iket) += zijkl;
          }
        }
      }

      // IIa_pgqh = \sum_ad Chi_III^J0_pgda * Gamma^J0_daqh
      // IIc_pgqh = \sum_ad Gamma^J0_pgad * Chi_III^J0_adqh
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();
        arma::Mat<double> Multi_matirx(nKets, nKets * 2);
        Multi_matirx = Chi_III[ch] * Gamma_matrix[ch];
        for (int ibra = 0; ibra < nKets; ++ibra)
        {
          Ket &bra = tbc.GetKet(ibra);
          size_t p = bra.p;
          size_t g = bra.q;

          for (int iket = ibra; iket < nKets; ++iket)
          {
            Ket &ket = tbc.GetKet(iket);
            size_t q = ket.p;
            size_t h = ket.q;

            double z_IIa = Multi_matirx(ibra, iket);
            double z_IIc = hZ * Multi_matirx(iket, ibra);

            if (p == g)
            {
              z_IIa /= PhysConst::SQRT2;
              z_IIc /= PhysConst::SQRT2;
            }
            if (q == h)
            {
              z_IIa /= PhysConst::SQRT2;
              z_IIc /= PhysConst::SQRT2;
            }
            Z2.AddToTBME(ch, ch, ibra, iket, z_IIa + z_IIc);
          } // iket
        }   // ibra
      }     // J0 channel

      Z.profiler.timer[std::string(__func__) + " Diagram IIa and IIc"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      // ####################################################################################
      //                      Factorization of IIb and IId
      // ####################################################################################
      //   diagram IIb
      //
      //   II(b)^J0_pghq = - P(p/g) * P(q/h) * \sum_{ad J4 J5} ( 2 * J4 + 1 ) ( 2 * J5 + 1 )
      //
      //                   { J0 J4 J5 } { J0 J4 J5 }
      //                   { jd jq jh } { ja jp jg }
      //
      //                   (\bar(Chi_III)^J5_dqap)^T Gamma^J4_gahd
      //
      //-------------------------------------------------------------------------------------
      //   diagram IId
      //
      //   II(d)^J0_pgqh = - P(p/g) * P(q/h) * \sum_{abcd J4 J5} ( 2 * J4 + 1 ) ( 2 * J5 + 1 )
      //
      //                   { J0 J5 J4 } { J5 J4 J0 }
      //                   { jd jp jg } { jq jh ja }
      //
      //                   (-)^(jd + jg + jh + ja) * (\bar(Chi_III)^J5_dgha)^T Gamma^J4_dpaq
      //-------------------------------------------------------------------------------------

      std::deque<arma::mat> CHI_III_final(nch);
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int nKets = tbc.GetNumberKets();
        // Not symmetric
        CHI_III_final[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      }

      /// Pandya transformation only recouple the angula momentum
      /// IIb and IId
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
          double ja = oa.j2 * 0.5;

          Orbit &ob = Z.modelspace->GetOrbit(b);
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
            double jc = oc.j2 * 0.5;
            Orbit &od = Z.modelspace->GetOrbit(d);
            double jd = od.j2 * 0.5;

            int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
            int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
            double Xbar = 0;

            for (int J_std = jmin; J_std <= jmax; J_std++)
            {
              int parity_cc = (oa.l + od.l) % 2;
              int Tz_cc = std::abs(oa.tz2 - od.tz2) / 2;
              int ch_cc_old = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_cc, Tz_cc);

              TwoBodyChannel_CC &tbc_cc_old = Z.modelspace->GetTwoBodyChannel_CC(ch_cc_old);
              int nkets = tbc_cc_old.GetNumberKets();
              int indx_ad = tbc_cc_old.GetLocalIndex(std::min(int(a), int(d)), std::max(int(a), int(d)));
              int indx_bc = tbc_cc_old.GetLocalIndex(std::min(int(b), int(c)), std::max(int(b), int(c)));
              if (indx_ad < 0 or indx_bc < 0)
                continue;

              if (a > d)
                indx_ad += nkets;
              if (b > c)
                indx_bc += nkets;

              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                Xbar -= Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J_std) * (2 * J_std + 1) * sixj1 * (barCHI_III[ch_cc_old](indx_bc, indx_ad) + barCHI_III[ch_cc_old](indx_ad, indx_bc));
              }
            }

            barCHI_III_RC[ch_cc](ibra_cc, iket_cc) = Xbar;
          }

          //-------------------
        }
      }

#pragma omp parallel for
      for (int ch = 0; ch < n_nonzero; ++ch)
      {
        CHI_III_final[ch] = bar_Gamma[ch] * barCHI_III_RC[ch];
      }

      //  Inverse Pandya transformation
      //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
      //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
      //                                    { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
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
          int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

          for (int iket = ibra; iket < nKets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc.GetKet(iket);
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

            Z2.AddToTBME(ch, ch, ibra, iket, phaseFactor * zijkl);
          }
        }
      }

      Z.profiler.timer[std::string(__func__) + " Diagram IIb and IId"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      // ####################################################################################
      //                      Factorization of IIe and IIf
      // ####################################################################################

      // ###########################################################
      //
      // The intermediate two body operator
      //  Chi_IV :
      //        q  | eta |  b
      //           |_____|
      //        a  |_____|  c
      //           | eta |
      //        p  |     |  d
      // ###########################################################

      std::deque<arma::mat> CHI_IV(nch);
      std::deque<arma::mat> CHI_IV_II(nch);
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int nKets = tbc.GetNumberKets();
        // Not symmetric
        CHI_IV[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
        CHI_IV_II[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      }
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();
        CHI_IV[ch] = Eta_matrix[ch] * Eta_matrix_c[ch];
        CHI_IV_II[ch] = (Eta_matrix[ch] * Eta_matrix_d[ch]).t();
      }

      // build bar_CHI_IV
      std::deque<arma::mat> bar_CHI_IV(n_nonzero);
      std::deque<arma::mat> bar_CHI_IV_II(n_nonzero);
      std::deque<arma::mat> bar_CHI_gamma(n_nonzero);
      std::deque<arma::mat> bar_CHI_gamma_II(n_nonzero);
      /// initial bar_CHI_IV
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        bar_CHI_IV[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        bar_CHI_IV_II[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        bar_CHI_gamma[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        bar_CHI_gamma_II[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      /// Pandya transformation
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
            double jc = oc.j2 * 0.5;
            double n_c = oc.occ;
            double nbar_c = 1.0 - n_c;

            Orbit &od = Z.modelspace->GetOrbit(d);
            double jd = od.j2 * 0.5;
            double n_d = od.occ;
            double nbar_d = 1.0 - n_d;

            int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
            int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
            double Xbar = 0;
            double Xbar1 = 0;
            int dJ_std = 1;
            if ((a == d or b == c))
            {
              dJ_std = 2;
              jmin += jmin % 2;
            }
            for (int J_std = jmin; J_std <= jmax; J_std += dJ_std)
            {

              int parity = (oa.l + od.l) % 2;
              int Tz = (oa.tz2 + od.tz2) / 2;
              int ch = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity, Tz);
              TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
              int nkets = tbc.GetNumberKets();
              int indx_ad = tbc.GetLocalIndex(std::min(int(a), int(d)), std::max(int(a), int(d)));
              int indx_cb = tbc.GetLocalIndex(std::min(int(c), int(b)), std::max(int(c), int(b)));
              if (indx_ad < 0 or indx_cb < 0)
                continue;
              if (a > d)
                indx_ad += nkets;
              if (c > b)
                indx_cb += nkets;

              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                Xbar -= (2 * J_std + 1) * sixj1 * CHI_IV[ch](indx_ad, indx_cb);
                Xbar1 -= (2 * J_std + 1) * sixj1 * CHI_IV_II[ch](indx_ad, indx_cb);
              }
            }

            bar_CHI_IV[ch_cc](ibra_cc, iket_cc) = Xbar;
            bar_CHI_IV_II[ch_cc](ibra_cc, iket_cc) = Xbar1;
          }

          //-------------------
        }
      }

      // calculate bat_chi_IV * bar_gamma
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        bar_CHI_gamma[ch_cc] = bar_CHI_IV[ch_cc] * bar_Gamma[ch_cc];
        bar_CHI_gamma_II[ch_cc] = bar_CHI_IV_II[ch_cc] * bar_Gamma[ch_cc];
      }

      //  Diagram e
      //  II(e)^J_ijkl  = - 1/2 ( 1- P_ij ) ( 1- P_kl ) sum_J' (2J'+1)  { i j J }  \bar{bar_CHI_gamma}^J'_il`kj`
      //                                                                { k l J'}
      //  Diagram f
      //  II(f)^J_ijkl  = - 1/2 ( 1- P_ij ) ( 1- P_kl ) sum_J' (2J'+1)  { i j J }  \bar{bar_CHI_gamma_II}^J'_il`kj`
      //
      //  Inverse Pandya transformation                                                              { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
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

          for (int iket = ibra; iket < nKets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc.GetKet(iket);
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
              int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
              int indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k));
              if (indx_il < 0 or indx_kj < 0)
                continue;
              // jilk, exchange i and j, k and l     ->  jk   li
              int indx_jk = indx_kj + (j > k ? nkets_cc : 0);
              int indx_li = indx_il + (l > i ? nkets_cc : 0);
              double me1 = bar_CHI_gamma[ch_cc](indx_jk, indx_li) + bar_CHI_gamma_II[ch_cc](indx_jk, indx_li);
              commjilk -= (2 * Jprime + 1) * sixj * me1;

              // ijkl direct term
              indx_il += (i > l ? nkets_cc : 0);
              indx_kj += (k > j ? nkets_cc : 0);
              me1 = bar_CHI_gamma[ch_cc](indx_il, indx_kj) + bar_CHI_gamma_II[ch_cc](indx_il, indx_kj);
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
              int indx_ki = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k));
              int indx_jl = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j));
              if (indx_ki < 0 or indx_jl < 0)
                continue;

              // ijlk, exchange k and l     ->  ik lj
              int indx_ik = indx_ki + (i > k ? nkets_cc : 0);
              int indx_lj = indx_jl + (l > j ? nkets_cc : 0);
              double me1 = bar_CHI_gamma[ch_cc](indx_ik, indx_lj) + bar_CHI_gamma_II[ch_cc](indx_ik, indx_lj);
              commijlk -= (2 * Jprime + 1) * sixj * me1;

              // jikl, exchange i and j    ->  jl ki
              indx_ki += (k > i ? nkets_cc : 0);
              indx_jl += (j > l ? nkets_cc : 0);
              me1 = bar_CHI_gamma[ch_cc](indx_jl, indx_ki) + +bar_CHI_gamma_II[ch_cc](indx_jl, indx_ki);
              commjikl -= (2 * Jprime + 1) * sixj * me1;
            }

            double zijkl = (commijkl - Z.modelspace->phase((ji + jj) / 2 - J0) * commjikl);
            zijkl += (-Z.modelspace->phase((jl + jk) / 2 - J0) * commijlk + Z.modelspace->phase((jk + jl + ji + jj) / 2) * commjilk);

            if (i == j)
              zijkl /= PhysConst::SQRT2;
            if (k == l)
              zijkl /= PhysConst::SQRT2;

            Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zijkl);
          }
        }
      }

      Z.profiler.timer[std::string(__func__) + " Diagram IIe and IIf"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      Z.profiler.timer[std::string(__func__) + " Diagram II"] += omp_get_wtime() - t_type;
      t_type = omp_get_wtime();
    }
    // *********************************************************************************** //
    //                                Diagram III                                          //
    // *********************************************************************************** //
    if (use_TypeIII_2b)
    {
      //------------------------------------------------------------------------------
      //  The intermediate two body operator
      //  Chi_V :
      //            eta |
      //           _____|
      //          /\    |
      //         (  )     |
      //          \/~~~~~~|
      //            gamma |
      //
      //  Pandya transform
      //  \bar{X}^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
      //                              { k l J'}
      //
      //  \bar{Chi}_V = \sum_{ab} (nbar_a * nbar_c * n_b - nbar_b * n_a * n_c )
      //                \bar{Gamma}_dq`ba` * \bar{eta}_ba`cp`
      //-------------------------------------------------------------------------------
      // The operator nnnbar_Eta and bar_Gamma are used to construct Chi_V
      std::deque<arma::mat> bar_CHI_V(n_nonzero);
      std::deque<arma::mat> bar_CHI_V_RC(n_nonzero);
      /// initial bar_CHI_V
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        bar_CHI_V[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        bar_CHI_V_RC[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      // build bar_CHI_V
#pragma omp parallel for
      for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
      {
        bar_CHI_V[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta[ch_cc];
      }

      // ################################################################################## //
      //                      Factorization of IIIa and IIIb                                //
      // ################################################################################## //
      //   diagram IIIa
      //
      //   III(a)^J0_pgqh = P(p/g) * P(q/h) * \sum_{cd J3 J5}
      //                   ( 2 * J3 + 1 ) ( 2 * J5 + 1 )
      //
      //                   { J3 J0 J5 } { J3 J0 J5 }
      //                   { jp jc jg } { jq jd jh }
      //
      //                   bar_CHI_V^J5_dq`cp` eta^J3_gchd
      //
      //   diagram IIIb
      //
      //   III(b)^J0_pgqh = - P(p/g) * P(q/h) * \sum_{cd J3 J5}
      //                    ( 2 * J3 + 1 ) ( 2 * J5 + 1 )
      //
      //                    { J0 J3 J5 } { J0 J3 J5 }
      //                    { jc jq jh } { jd jp jg }
      //
      //                    bar_CHI_V^J5_cq`dp` eta^J3_gdhc
      // ####################################################################################

      /// Pandya transformation only recouple the angula momentum
      /// diagram IIIa - diagram IIIb
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
            double jc = oc.j2 * 0.5;
            double n_c = oc.occ;
            double nbar_c = 1.0 - n_c;

            Orbit &od = Z.modelspace->GetOrbit(d);
            double jd = od.j2 * 0.5;
            double n_d = od.occ;
            double nbar_d = 1.0 - n_d;

            int jmin = std::max(std::abs(oa.j2 - od.j2), std::abs(oc.j2 - ob.j2)) / 2;
            int jmax = std::min(oa.j2 + od.j2, oc.j2 + ob.j2) / 2;
            double Xbar = 0;

            for (int J_std = jmin; J_std <= jmax; J_std++)
            {
              int parity_cc = (oa.l + od.l) % 2;
              int Tz_cc = std::abs(oa.tz2 - od.tz2) / 2;
              int ch_cc_old = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_cc, Tz_cc);

              TwoBodyChannel_CC &tbc_cc_old = Z.modelspace->GetTwoBodyChannel_CC(ch_cc_old);
              int nkets = tbc_cc_old.GetNumberKets();
              int indx_ad = tbc_cc_old.GetLocalIndex(std::min(int(a), int(d)), std::max(int(a), int(d)));
              int indx_bc = tbc_cc_old.GetLocalIndex(std::min(int(b), int(c)), std::max(int(b), int(c)));
              if (indx_ad < 0 or indx_bc < 0)
                continue;

              if (a > d)
                indx_ad += nkets;
              if (b > c)
                indx_bc += nkets;

              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                Xbar += Z.modelspace->phase((ob.j2 + oc.j2) / 2 + J_std) * (2 * J_std + 1) * sixj1 * (bar_CHI_V[ch_cc_old](indx_ad, indx_bc) - hZ * bar_CHI_V[ch_cc_old](indx_bc, indx_ad));
              }
            }

            bar_CHI_V_RC[ch_cc](ibra_cc, iket_cc) = Xbar;
          }

          //-------------------
        }
      }

      std::deque<arma::mat> CHI_V_final(nch);
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int nKets = tbc.GetNumberKets();
        // Not symmetric
        CHI_V_final[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      }

#pragma omp parallel for
      for (int ch = 0; ch < n_nonzero; ++ch)
      {
        CHI_V_final[ch] = bar_Eta[ch] * bar_CHI_V_RC[ch];
      }

      //  Inverse Pandya transformation
      //  X^J_ijkl  = - ( 1- P_ij ) ( 1- P_kl ) (-)^{J + ji + jj}  sum_J' (2J'+1)
      //                (-)^{J' + ji + jk}  { j i J }  \bar{X}^J'_jl`ki`
      //                                    { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
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
          int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

          for (int iket = ibra; iket < nKets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc.GetKet(iket);
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

            Z2.AddToTBME(ch, ch, ibra, iket, phaseFactor * zijkl);
          }
        }
      }

      Z.profiler.timer[std::string(__func__) + " Diagram IIIa and IIIb"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      // ####################################################################################
      //                      Factorization of IIIc and IIId
      // ####################################################################################

      //------------------------------------------------------------------------------
      //  The intermediate two body operator
      //  Chi_V :
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

      std::deque<arma::mat> bar_CHI_VI(n_nonzero);
      /// initial bar_Eta and nnnbar_Eta
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        bar_CHI_VI[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      // BUILD bar_CHI_VI
      // bar_CHI_VI = (bar{na} nb nc + bar{nb} bar{nc} na) \bar{Gamma}_dhba \bar{Eta}_baqc
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        bar_CHI_VI[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta_d[ch_cc];
      }

      std::deque<arma::mat> CHI_VI(nch);
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int nKets = tbc.GetNumberKets();
        // Not symmetric
        CHI_VI[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      }

      // BUILD CHI_VI
      // Inverse Pandya transformation
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();
        for (int ibra = 0; ibra < nKets * 2; ++ibra)
        {
          size_t i, j;
          if (ibra < nKets)
          {
            Ket &bra = tbc.GetKet(ibra);
            i = bra.p;
            j = bra.q;
          }
          else
          {
            Ket &bra = tbc.GetKet(ibra - nKets);
            i = bra.q;
            j = bra.p;
          }
          if (ibra >= nKets and i == j)
            continue;

          Orbit &oi = Z.modelspace->GetOrbit(i);
          int ji = oi.j2;
          Orbit &oj = Z.modelspace->GetOrbit(j);
          int jj = oj.j2;

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
            double commijkl = 0;
            double commijlk = 0;

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
              int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l));
              int indx_kj = tbc_cc.GetLocalIndex(std::min(k, j), std::max(k, j));
              if (indx_il < 0 or indx_kj < 0)
                continue;

              indx_il += (i > l ? nkets_cc : 0);
              indx_kj += (k > j ? nkets_cc : 0);
              double me1 = bar_CHI_VI[ch_cc](indx_il, indx_kj);
              commijkl -= (2 * Jprime + 1) * sixj * me1;
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

              // exchange k and l
              indx_ik += (i > k ? nkets_cc : 0);
              indx_lj += (l > j ? nkets_cc : 0);
              double me2 = bar_CHI_VI[ch_cc](indx_ik, indx_lj);
              commijlk -= (2 * Jprime + 1) * sixj * me2;
            }

            double zijkl = (commijkl - Z.modelspace->phase((jk + jl) / 2 - J0) * commijlk);

            CHI_VI[ch](ibra, iket) += zijkl;
          }
        }
      }

// Diagram IIIc and Diagram IIId
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();
        if (nKets < 1)
          continue;
        arma::Mat<double> Multi_matirx(nKets * 2, nKets * 2);
        Multi_matirx = Eta_matrix[ch] * CHI_VI[ch];
        for (int ibra = 0; ibra < nKets; ++ibra)
        {
          Ket &bra = tbc.GetKet(ibra);
          size_t p = bra.p;
          size_t g = bra.q;

          for (int iket = ibra; iket < nKets; ++iket)
          {
            Ket &ket = tbc.GetKet(iket);
            size_t q = ket.p;
            size_t h = ket.q;

            double zpgqhIIIc = Multi_matirx(ibra, iket);
            double zpgqhIIId = hZ * Multi_matirx(iket, ibra);

            if (p == g)
            {
              zpgqhIIIc /= PhysConst::SQRT2;
              zpgqhIIId /= PhysConst::SQRT2;
            }
            if (q == h)
            {
              zpgqhIIIc /= PhysConst::SQRT2;
              zpgqhIIId /= PhysConst::SQRT2;
            }
            Z2.AddToTBME(ch, ch, ibra, iket, -zpgqhIIIc - zpgqhIIId);
          } // iket
        }   // ibra

      } // J0 channel

      Z.profiler.timer[std::string(__func__) + " Diagram IIIc and IIId"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      // ####################################################################################
      //                      Factorization of IIIe and IIIf
      // ####################################################################################

      // ###########################################################
      //
      // The intermediate two body operator
      //  CHI_VII :
      //        g  |     |  c
      //           |~~~~~|
      //        a  |     |  b
      //           |_____|
      //        h  |     |  d
      // ###########################################################

      std::deque<arma::mat> CHI_VII(nch);
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int nKets = tbc.GetNumberKets();
        // Not symmetric
        CHI_VII[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
      }

#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J0 = tbc.J;
        int nKets = tbc.GetNumberKets();
        CHI_VII[ch] = Gamma_matrix[ch] * Eta_matrix_d[ch];
      }

      std::deque<arma::mat> bar_CHI_VII_CC(n_nonzero);
      std::deque<arma::mat> bar_CHI_VII_CC_ef(n_nonzero);
      /// initial bar_Eta and nnnbar_Eta
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
        int nKets_cc = tbc_cc.GetNumberKets();
        // because the restriction a<b in the bar and ket vector, if we want to store the full
        // Pandya transformed matrix, we twice the size of matrix
        bar_CHI_VII_CC[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
        bar_CHI_VII_CC_ef[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      }

      /// Modified Pandya transformation
      //  \bar{X}^J_ij`kl` =   sum_J' (-)^{jk + jj + J'} { i j J } (2J'+1) X^J'_ilkj
      //                                                 { k l J'}
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
              int phaseFactor = Z.modelspace->phase(J_std + (oc.j2 + ob.j2) / 2);
              double sixj1 = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
              if (std::abs(sixj1) > 1e-8)
              {
                int parity_J2 = (ob.l + oc.l) % 2;
                int Tz_J2 = (ob.tz2 + oc.tz2) / 2;
                int ch_J2 = Z.modelspace->GetTwoBodyChannelIndex(J_std, parity_J2, Tz_J2);
                TwoBodyChannel &tbc_J2 = Z.modelspace->GetTwoBodyChannel(ch_J2);
                int nkets = tbc_J2.GetNumberKets();
                int indx_ad = tbc_J2.GetLocalIndex(std::min(int(a), int(d)), std::max(int(a), int(d)));
                int indx_bc = tbc_J2.GetLocalIndex(std::min(int(b), int(c)), std::max(int(b), int(c)));
                if (indx_ad < 0 or indx_bc < 0)
                  continue;
                if (a > d)
                  indx_ad += nkets;
                if (b > c)
                  indx_bc += nkets;
                Ybar += phaseFactor * (2 * J_std + 1) * sixj1 * CHI_VII[ch_J2](indx_ad, indx_bc);
              }
            }
            bar_CHI_VII_CC[ch_cc](ibra_cc, iket_cc) = Ybar;
          }
          //-------------------
        }
      }

      // Diagram IIIe and Diagram IIIf
      Z.profiler.timer[std::string(__func__) + " Diagram IIIe and IIIf  -a"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      // bar_CHI_VII_CC_ef = bar_CHI_VII_CC * bar_Eta
#pragma omp parallel for
      for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
      {
        bar_CHI_VII_CC_ef[ch_cc] = bar_CHI_VII_CC[ch_cc] * bar_Eta[ch_cc];
        bar_CHI_VII_CC_ef[ch_cc] += hZ * bar_Eta[ch_cc].t() * bar_CHI_VII_CC[ch_cc].t();
      }

      Z.profiler.timer[std::string(__func__) + " Diagram IIIe and IIIf  -b"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      //  Inverse Pandya transformation
      //  X^J_ijkl  = - 1/2 ( 1- P_ij ) ( 1- P_kl ) sum_J' (2J'+1) (-)^{J0 + jp + jg}
      //                { j i J }  \bar{X}^J'_jl`ki`
      //                { k l J'}
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
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
          int phaseFactor = Z.modelspace->phase(J0 + (ji + jj) / 2);

          for (int iket = ibra; iket < nKets; ++iket)
          {
            size_t k, l;
            Ket &ket = tbc.GetKet(iket);
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
              // if (std::abs(sixj) < 1e-8)
              //   continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
              TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
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

            Z2.AddToTBME(ch, ch, ibra, iket, phaseFactor * 0.5 * zijkl);
          }
        }
      }

      Z.profiler.timer[std::string(__func__) + " Diagram IIIe and IIIf  -c"] += omp_get_wtime() - t_internal;
      t_internal = omp_get_wtime();

      Z.profiler.timer[std::string(__func__) + " Diagram III"] += omp_get_wtime() - t_type;
    }
    // Timer
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return;
  }




 } // namespace FactorizedDoubleCommutator
} // namespace Commutator
