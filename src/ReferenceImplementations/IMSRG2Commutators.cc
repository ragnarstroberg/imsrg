
#include "ReferenceImplementations.hh"
#include "ModelSpace.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"


/// Straightforward implementation of J-coupled commutator expressions
/// without optimizations. This should be benchmarked against the
/// mscheme implementation and then left untouched.
namespace ReferenceImplementations
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Z0 =  sum_ab (2ja+1) (na-nb) (Xab Yba )
  ///                   =  sum_ab (2ja+1) na (1-nb) (Xab Yba - Yab Xba)
  ///
  ///
  void comm110ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    double z0 = 0;
    for (size_t a : Z.modelspace->all_orbits)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (size_t b : Z.modelspace->all_orbits)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        z0 += (oa.j2 + 1) * oa.occ * (1 - ob.occ) * (X1(a, b) * Y1(b, a) - Y1(a, b) * X1(b, a));
      }
    }
    Z.ZeroBody += z0;
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression : Z0 = 1/4 sum_abcd sum_J (2J+1) ( n_a n_b nbar_c nbar_d ) ( x_abcd y_cdab - y_abcd x_cdab )
  ///
  ///
  void comm220ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    double z0 = 0;

    for (size_t a : Z.modelspace->all_orbits)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (size_t b : Z.modelspace->all_orbits)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);

        for (size_t c : Z.modelspace->all_orbits)
        {
          Orbit &oc = Z.modelspace->GetOrbit(c);
          for (size_t d : Z.modelspace->all_orbits)
          {
            Orbit &od = Z.modelspace->GetOrbit(d);

            int Jmin = std::max(std::abs(oa.j2 - ob.j2), std::abs(oc.j2 - od.j2)) / 2;
            int Jmax = std::min(oa.j2 + ob.j2, oc.j2 + od.j2) / 2;
            for (int J = Jmin; J <= Jmax; J++)
            {
              double xabcd = X2.GetTBME_J(J, J, a, b, c, d);
              double xcdab = X2.GetTBME_J(J, J, c, d, a, b);
              double yabcd = Y2.GetTBME_J(J, J, a, b, c, d);
              double ycdab = Y2.GetTBME_J(J, J, c, d, a, b);
              z0 += 1. / 4 * (2 * J + 1) * (oa.occ * ob.occ * (1 - oc.occ) * (1 - od.occ)) * (xabcd * ycdab - yabcd * xcdab);
            }
          }
        }
      }
    }
    Z.ZeroBody += z0;
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij = sum_a  (Xia Yaj - Yia Xaj)
  ///
  ///
  void comm111ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    auto &Z1 = Z.OneBody;
    for (size_t i : Z.modelspace->all_orbits)
    {
      for (size_t j : Z.modelspace->all_orbits)
      {
        for (size_t a : Z.modelspace->all_orbits)
        {
          Z1(i, j) += X1(i, a) * Y1(a, j) - Y1(i, a) * X1(a, j);
        }
      }
    }
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij = sum_ab sum_J (2J+1)/(2ji+1) (na-nb) ( Xab YJbiaj - Yab XJbiaj )
  ///
  ///
  void comm121ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &Z1 = Z.OneBody;

    for (size_t i : Z.modelspace->all_orbits)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        Orbit &oj = Z.modelspace->GetOrbit(j);

        for (size_t a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          for (size_t b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            int Jmin = std::max(std::abs(ob.j2 - oi.j2), std::abs(oa.j2 - oj.j2)) / 2;
            int Jmax = std::min(ob.j2 + oi.j2, oa.j2 + oj.j2) / 2;
            for (int J = Jmin; J <= Jmax; J++)
            {
              double xbiaj = X2.GetTBME_J(J, J, b, i, a, j);
              double ybiaj = Y2.GetTBME_J(J, J, b, i, a, j);
              Z1(i, j) += (2 * J + 1) / (oi.j2 + 1.) * (oa.occ - ob.occ) * (X1(a, b) * ybiaj - Y1(a, b) * xbiaj);
            }
          }
        }
      }
    }
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij = 1/2 sum_abc (na nb (1-nc) + (1-na)(1-nb)nc) (2J+1)/(2ji+1) (XJciab YJabcj - YJciab XJabcj)
  ///
  ///
  void comm221ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &Z1 = Z.OneBody;

    for (size_t i : Z.modelspace->all_orbits)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        Orbit &oj = Z.modelspace->GetOrbit(j);

            for (size_t c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);


        for (size_t a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          for (size_t b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
//            for (size_t c : Z.modelspace->all_orbits)
//            {
//              Orbit &oc = Z.modelspace->GetOrbit(c);

              //               if ( (oi.l+oc.l + oa.l+ob.l)%2>0) continue;
              //               if ( (oi.tz2+oc.tz2) != (oa.tz2+ob.tz2) ) continue;

              int Jmin = std::max(std::abs(ob.j2 - oa.j2), std::abs(oc.j2 - oi.j2)) / 2;
              int Jmax = std::min(ob.j2 + oa.j2, oc.j2 + oi.j2) / 2;
              for (int J = Jmin; J <= Jmax; J++)
              {
                double xciab = X2.GetTBME_J(J, J, c, i, a, b);
                double yciab = Y2.GetTBME_J(J, J, c, i, a, b);
                double xabcj = X2.GetTBME_J(J, J, a, b, c, j);
                double yabcj = Y2.GetTBME_J(J, J, a, b, c, j);
                double zij = 1. / 2 * (2 * J + 1) / (oi.j2 + 1.) * (oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ) * (xciab * yabcj - yciab * xabcj);
                Z1(i, j) += zij;
              } // J
            } // c
          } // b
        } // a
      } // j
    } // i

//    Z.PrintOneBody();
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJijkl = sum_a (Xia YJajkl + Xja YJiakl - YJijal Xak - YJijka Xal) -   {X<->Y}
  ///
  ///
  /// SRS Modified on Nov 3 2022 to accommodate operators that change parity or isospin.
  void comm122ss(const Operator &X, const Operator &Y, Operator &Z)
  {

    double t_start = omp_get_wtime();
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &Z2 = Z.TwoBody;

    std::vector<size_t> ch_bra_list, ch_ket_list;
    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();

    //   int nch = Z.modelspace->GetNumberTwoBodyChannels();
    #pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;

          double zijkl = 0;
          for (size_t a : Z.modelspace->all_orbits)
          {
            //           Orbit& oa = Z.modelspace->GetOrbit(a);
            double xajkl = X2.GetTBME_J(J, J, a, j, k, l);
            double yajkl = Y2.GetTBME_J(J, J, a, j, k, l);
            double xiakl = X2.GetTBME_J(J, J, i, a, k, l);
            double yiakl = Y2.GetTBME_J(J, J, i, a, k, l);
            double xijal = X2.GetTBME_J(J, J, i, j, a, l);
            double yijal = Y2.GetTBME_J(J, J, i, j, a, l);
            double xijka = X2.GetTBME_J(J, J, i, j, k, a);
            double yijka = Y2.GetTBME_J(J, J, i, j, k, a);

            zijkl += X1(i, a) * yajkl + X1(j, a) * yiakl - yijal * X1(a, k) - yijka * X1(a, l);
            zijkl -= Y1(i, a) * xajkl + Y1(j, a) * xiakl - xijal * Y1(a, k) - xijka * Y1(a, l);
          } // a
          // Need to normalize here, because AddToTBME expects a normalized TBME.
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);

        } // iket
      } // ibra
    } // ch
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJijkl = 1/2 sum_ab ((1-na)(1-nb) - na nb)  (XJijab YJabkl - YJijab XJabkl)
  ///
  ///

  void comm222_pp_hhss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &Z2 = Z.TwoBody;

    std::vector<size_t> ch_bra_list, ch_ket_list;

    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();
    size_t norb = Z.modelspace->GetNumberOrbits();

    #pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;

        size_t ket_min = (ch_bra == ch_ket) ? ibra : 0;
        for (int iket = ket_min; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;

          double zijkl = 0;
          for (size_t a : Z.modelspace->all_orbits)
          {
            for (size_t b : Z.modelspace->all_orbits)
            {
              if (b < a)
                continue;
              Orbit &oa = Z.modelspace->GetOrbit(a);
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double xijab = X2.GetTBME_J(J, J, i, j, a, b);
              double yijab = Y2.GetTBME_J(J, J, i, j, a, b);
              double xabkl = X2.GetTBME_J(J, J, a, b, k, l);
              double yabkl = Y2.GetTBME_J(J, J, a, b, k, l);
              double flip_factor = (a == b) ? 1 : 2; // looping over kets only gets a<=b. So we need a factor of 2 for the other ordering.
              zijkl += 1. / 2 * flip_factor * ((1 - oa.occ) * (1 - ob.occ) - oa.occ * ob.occ) * (xijab * yabkl - yijab * xabkl);
            }
          }
          // Need to normalize here, because AddToTBME expects a normalized TBME.
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);

        } // iket
      } // ibra
    } // ch
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  } // comm222_pp_hhss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJijkl =   sum_ab (na-nb) sum_J'  (2J'+1) (1- (-1)^{i+j-J} Pij) { i j J  } ( ~XJ'ilab ~YJ'abkj - ~YJ'ilab ~XJ'abkj)
  ///                                                                                 { k l J' }
  ///           in this expression, the ~X and ~Y are the Pandya-transformed matrix elements
  ///           ~XJabcd = - sum_J' (2J'+1) { a b J  } XJ'_adcb
  ///                                      { c d J' }
  ///
  void comm222_phss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &Z2 = Z.TwoBody;
    if (X.IsReduced() or Y.IsReduced())
    {
      std::cout << "  " << __FILE__ << "  line " << __LINE__ << "  calling " << __func__ << "  with operators that are in reduced form. Things will be wrong." << std::endl;
    }

    std::vector<size_t> ch_bra_list, ch_ket_list;
    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();

    //   int nch = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);

      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      //     for (int ibra=0; ibra<nkets; ibra++)
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);

          double zijkl = 0;
          int Jpmin = std::min(std::max(std::abs(oi.j2 - ol.j2), std::abs(oj.j2 - ok.j2)),
                               std::max(std::abs(oj.j2 - ol.j2), std::abs(oi.j2 - ok.j2))) /
                      2;
          int Jpmax = std::max(std::min(oi.j2 + ol.j2, oj.j2 + ok.j2),
                               std::min(oj.j2 + ol.j2, oi.j2 + ok.j2)) /
                      2;

          for (int Jp = Jpmin; Jp <= Jpmax; Jp++)
          {

            double zbar_ilkj = 0;
            double zbar_jlki = 0;

            for (size_t a : Z.modelspace->all_orbits)
            {
              Orbit &oa = Z.modelspace->GetOrbit(a);
              for (size_t b : Z.modelspace->all_orbits)
              {
                Orbit &ob = Z.modelspace->GetOrbit(b);

                // "Direct" term
                double xbar_ilab = 0;
                double ybar_ilab = 0;

                int Jppmin = std::max(std::abs(oi.j2 - ob.j2), std::abs(oa.j2 - ol.j2)) / 2;
                int Jppmax = std::min(oi.j2 + ob.j2, oa.j2 + ol.j2) / 2;
                for (int Jpp = Jppmin; Jpp <= Jppmax; Jpp++)
                {
                  double sixj_ilab = AngMom::SixJ(oi.j2 * 0.5, ol.j2 * 0.5, Jp, oa.j2 * 0.5, ob.j2 * 0.5, Jpp);
                  double xibal = X2.GetTBME_J(Jpp, i, b, a, l);
                  double yibal = Y2.GetTBME_J(Jpp, i, b, a, l);
                  xbar_ilab -= (2 * Jpp + 1) * sixj_ilab * xibal;
                  ybar_ilab -= (2 * Jpp + 1) * sixj_ilab * yibal;
                }

                double xbar_abkj = 0;
                double ybar_abkj = 0;
                Jppmin = std::max(std::abs(oa.j2 - oj.j2), std::abs(ok.j2 - ob.j2)) / 2;
                Jppmax = std::min(oa.j2 + oj.j2, ok.j2 + ob.j2) / 2;
                for (int Jpp = Jppmin; Jpp <= Jppmax; Jpp++)
                {
                  double sixj_abkj = AngMom::SixJ(oa.j2 * 0.5, ob.j2 * 0.5, Jp, ok.j2 * 0.5, oj.j2 * 0.5, Jpp);
                  double xajkb = X2.GetTBME_J(Jpp, a, j, k, b);
                  double yajkb = Y2.GetTBME_J(Jpp, a, j, k, b);
                  xbar_abkj -= (2 * Jpp + 1) * sixj_abkj * xajkb;
                  ybar_abkj -= (2 * Jpp + 1) * sixj_abkj * yajkb;
                }

                //                zbar_ilkj += (oa.occ - ob.occ) * (xbar_ilab * ybar_abkj - ybar_ilab * xbar_abkj);
                zbar_ilkj += (oa.occ - ob.occ) * (xbar_ilab * ybar_abkj - ybar_ilab * xbar_abkj);

                // "Exchange" term obtained by exchanging i<->j.
                double xbar_jlab = 0;
                double ybar_jlab = 0;
                Jppmin = std::max(std::abs(oj.j2 - ob.j2), std::abs(oa.j2 - ol.j2)) / 2;
                Jppmax = std::min(oj.j2 + ob.j2, oa.j2 + ol.j2) / 2;
                for (int Jpp = Jppmin; Jpp <= Jppmax; Jpp++)
                {
                  double sixj_jlab = AngMom::SixJ(oj.j2 * 0.5, ol.j2 * 0.5, Jp, oa.j2 * 0.5, ob.j2 * 0.5, Jpp);
                  double xjbal = X2.GetTBME_J(Jpp, j, b, a, l);
                  double yjbal = Y2.GetTBME_J(Jpp, j, b, a, l);
                  xbar_jlab -= (2 * Jpp + 1) * sixj_jlab * xjbal;
                  ybar_jlab -= (2 * Jpp + 1) * sixj_jlab * yjbal;
                }

                double xbar_abki = 0;
                double ybar_abki = 0;
                Jppmin = std::max(std::abs(oa.j2 - oi.j2), std::abs(ok.j2 - ob.j2)) / 2;
                Jppmax = std::min(oa.j2 + oi.j2, ok.j2 + ob.j2) / 2;
                for (int Jpp = Jppmin; Jpp <= Jppmax; Jpp++)
                {
                  double sixj_abki = AngMom::SixJ(oa.j2 * 0.5, ob.j2 * 0.5, Jp, ok.j2 * 0.5, oi.j2 * 0.5, Jpp);
                  double xaikb = X2.GetTBME_J(Jpp, a, i, k, b);
                  double yaikb = Y2.GetTBME_J(Jpp, a, i, k, b);
                  xbar_abki -= (2 * Jpp + 1) * sixj_abki * xaikb;
                  ybar_abki -= (2 * Jpp + 1) * sixj_abki * yaikb;
                }

                zbar_jlki += (oa.occ - ob.occ) * (xbar_jlab * ybar_abki - ybar_jlab * xbar_abki);

              } // b
            } // a

            double sixj_ijkl = AngMom::SixJ(oi.j2 * 0.5, oj.j2 * 0.5, J, ok.j2 * 0.5, ol.j2 * 0.5, Jp);
            double sixj_jikl = AngMom::SixJ(oj.j2 * 0.5, oi.j2 * 0.5, J, ok.j2 * 0.5, ol.j2 * 0.5, Jp);
            int phase_ij = AngMom::phase((oi.j2 + oj.j2 - 2 * J) / 2);
            zijkl += (2 * Jp + 1) * sixj_ijkl * zbar_ilkj;            // Direct
            zijkl -= (2 * Jp + 1) * sixj_jikl * zbar_jlki * phase_ij; // Exchange, with phase

          } // Jp

          // Need to normalize here, because AddToTBME expects a normalized TBME.
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);

        } // iket
      } // ibra
    } // ch

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  } // comm222_phss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// In the optimized version, an intermediate is generated
  ///  and used for both, but here we just call both sequentially.
  ///
  void comm222_pp_hh_221ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    comm222_pp_hhss(X, Y, Z);
    comm221ss(X, Y, Z);
  }


} // namespace ReferenceImplementations
////////////////////////////////////////////////////////////////////
