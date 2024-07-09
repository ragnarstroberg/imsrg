
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
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression : Z0 = 1/4 sum_abcd sum_J (2J+1) ( n_a n_b nbar_c nbar_d ) ( x_abcd y_cdab - y_abcd x_cdab )
  ///
  ///
  void comm220ss(const Operator &X, const Operator &Y, Operator &Z)
  {
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
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij = sum_a  (Xia Yaj - Yia Xaj)
  ///
  ///
  void comm111ss(const Operator &X, const Operator &Y, Operator &Z)
  {
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
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij = sum_ab sum_J (2J+1)/(2ji+1) (na-nb) ( Xab YJbiaj - Yab XJbiaj )
  ///
  ///
  void comm121ss(const Operator &X, const Operator &Y, Operator &Z)
  {
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
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij = 1/2 sum_abc (na nb (1-nc) + (1-na)(1-nb)nc) (2J+1)/(2ji+1) (XJciab YJabcj - YJciab XJabcj)
  ///
  ///
  void comm221ss(const Operator &X, const Operator &Y, Operator &Z)
  {
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
            for (size_t c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);

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
                Z1(i, j) += 1. / 2 * (2 * J + 1) / (oi.j2 + 1.) * (oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ) * (xciab * yabcj - yciab * xabcj);

              } // J
            } // c
          } // b
        } // a

      } // j
    } // i
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
    X.profiler.timer["ref_" + std::string(__func__)] += omp_get_wtime() - t_start;
  }

  /*
  void comm122ss( const Operator& X, const Operator& Y, Operator& Z )
  {
     auto& X1 = X.OneBody;
     auto& Y1 = Y.OneBody;
     auto& X2 = X.TwoBody;
     auto& Y2 = Y.TwoBody;
     auto& Z2 = Z.TwoBody;

     int nch = Z.modelspace->GetNumberTwoBodyChannels();
     #pragma omp parallel for schedule(dynamic,1)
     for (int ch=0; ch<nch; ch++)
     {
       TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
       int J = tbc.J;
       int nkets = tbc.GetNumberKets();
       for (int ibra=0; ibra<nkets; ibra++)
       {
         Ket& bra = tbc.GetKet(ibra);
         size_t i = bra.p;
         size_t j = bra.q;

         for (int iket=ibra; iket<nkets; iket++)
         {
           Ket& ket = tbc.GetKet(iket);
           size_t k = ket.p;
           size_t l = ket.q;

           double zijkl = 0;
           for ( size_t a : Z.modelspace->all_orbits )
           {
  //           Orbit& oa = Z.modelspace->GetOrbit(a);
             double xajkl = X2.GetTBME_J(J,J, a,j,k,l);
             double yajkl = Y2.GetTBME_J(J,J, a,j,k,l);
             double xiakl = X2.GetTBME_J(J,J, i,a,k,l);
             double yiakl = Y2.GetTBME_J(J,J, i,a,k,l);
             double xijal = X2.GetTBME_J(J,J, i,j,a,l);
             double yijal = Y2.GetTBME_J(J,J, i,j,a,l);
             double xijka = X2.GetTBME_J(J,J, i,j,k,a);
             double yijka = Y2.GetTBME_J(J,J, i,j,k,a);

             zijkl += X1(i,a) * yajkl + X1(j,a) * yiakl - yijal * X1(a,k) - yijka * X1(a,l);
             zijkl -= Y1(i,a) * xajkl + Y1(j,a) * xiakl - xijal * Y1(a,k) - xijka * Y1(a,l);
           }//a
           // Need to normalize here, because AddToTBME expects a normalized TBME.
           if (i==j) zijkl /= PhysConst::SQRT2;
           if (k==l) zijkl /= PhysConst::SQRT2;

           Z2.AddToTBME(ch,ch, ibra,iket, zijkl);

         }//iket
       }//ibra
     }//ch

  }*/

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJijkl = 1/2 sum_ab ((1-na)(1-nb) - na nb)  (XJijab YJabkl - YJijab XJabkl)
  ///
  ///

  void comm222_pp_hhss(const Operator &X, const Operator &Y, Operator &Z)
  {
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

      //   int nch = Z.modelspace->GetNumberTwoBodyChannels();
      //   #pragma omp parallel for schedule(dynamic,1)
      //   for (int ch=0; ch<nch; ch++)
      //   {
      //     TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
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
          //         int Jpmin = std::max( std::abs(oi.j2-ol.j2) , std::abs(oj.j2-ok.j2) )/2;
          //         int Jpmax = std::min( oi.j2+ol.j2 , oj.j2+ok.j2 )/2;
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

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Z0 = 1/36 sum_abcdef sum_J1J2J na nb nc (1-nd)(1-ne)(1-nf) (2J+1) (XJ1J2J_abcdef YJ1J2J_defabc - YJ1J2J_abcdef XJ1J2J_defabc)
  ///
  ///
  void comm330ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    double z0 = 0;
    if (Z.GetTRank() != 0)
      return; // zero body term has Tz=0
    if (Z.GetParity() != 0)
      return; // zero body term has even parity

    size_t norb = Z.modelspace->GetNumberOrbits();
    //   for ( size_t a : Z.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z0)
    for (size_t a = 0; a < norb; a++)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (size_t b : Z.modelspace->all_orbits)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);

        for (size_t c : Z.modelspace->all_orbits)
        {
          Orbit &oc = Z.modelspace->GetOrbit(c);
          if (std::abs(oa.occ * ob.occ * oc.occ) < 1e-7)
            continue;
          for (size_t d : Z.modelspace->all_orbits)
          {
            Orbit &od = Z.modelspace->GetOrbit(d);

            for (size_t e : Z.modelspace->all_orbits)
            {
              Orbit &oe = Z.modelspace->GetOrbit(e);
              for (size_t f : Z.modelspace->all_orbits)
              {
                Orbit &of = Z.modelspace->GetOrbit(f);

                if ((oa.l + ob.l + oc.l + od.l + oe.l + of.l) % 2 != X.GetParity())
                  continue;
                if (std::abs(oa.tz2 + ob.tz2 + oc.tz2 - od.tz2 - oe.tz2 - of.tz2) != X.GetTRank())
                  continue;
                if (std::abs(oa.occ * ob.occ * oc.occ * (1 - od.occ) * (1 - oe.occ) * (1 - of.occ)) < 1e-9)
                  continue;

                int J1min = std::abs(oa.j2 - ob.j2) / 2;
                int J2min = std::abs(od.j2 - oe.j2) / 2;
                int J1max = (oa.j2 + ob.j2) / 2;
                int J2max = (od.j2 + oe.j2) / 2;

                for (int J1 = J1min; J1 <= J1max; J1++)
                {
                  for (int J2 = J2min; J2 <= J2max; J2++)
                  {
                    int twoJmin = std::max(std::abs(oc.j2 - 2 * J1), std::abs(of.j2 - 2 * J2));
                    int twoJmax = std::min(oc.j2 + 2 * J1, of.j2 + 2 * J2);
                    for (int twoJ = twoJmin; twoJ <= twoJmax; twoJ += 2)
                    {
                      double xabcdef = X3.GetME_pn(J1, J2, twoJ, a, b, c, d, e, f);
                      double yabcdef = Y3.GetME_pn(J1, J2, twoJ, a, b, c, d, e, f);
                      double xdefabc = X3.GetME_pn(J2, J1, twoJ, d, e, f, a, b, c);
                      double ydefabc = Y3.GetME_pn(J2, J1, twoJ, d, e, f, a, b, c);
                      z0 += 1. / 36 * (twoJ + 1) * (oa.occ * ob.occ * oc.occ * (1 - od.occ) * (1 - oe.occ) * (1 - of.occ)) * (xabcdef * ydefabc - yabcdef * xdefabc);

                    } // twoJ
                  } // J2
                } // J1
              } // f
            } // e
          } // d
        } // c
      } // b
    } // a
    Z.ZeroBody += z0;
  } // comm330ss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij = 1/12 sum_abcde sum_J1J2J na nb (1-nc) (1-nd)(1-ne) (2J+1)/(2ji+1) (XJ1J2J_abicde YJ1J2J_cdeabj - YJ1J2J_abicde XJ1J2J_cdeabj)
  ///
  ///
  void comm331ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z1 = Z.OneBody;

    size_t norb = Z.modelspace->GetNumberOrbits();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < norb; i++)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double zij = 0;
        for (size_t a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          for (size_t b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);

            int J1min = std::abs(oa.j2 - ob.j2) / 2;
            int J1max = (oa.j2 + ob.j2) / 2;

            for (size_t c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              for (size_t d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                int J2min = std::abs(oc.j2 - od.j2) / 2;
                int J2max = (oc.j2 + od.j2) / 2;

                for (size_t e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);

                  double occupation_factor = oa.occ * ob.occ * (1 - oc.occ) * (1 - od.occ) * (1 - oe.occ) // fixed mistake found by Matthias Heinz Oct 2022
                                             + (1 - oa.occ) * (1 - ob.occ) * oc.occ * od.occ * oe.occ;

                  if (std::abs(occupation_factor) < 1e-8)
                    continue;

                  for (int J1 = J1min; J1 <= J1max; J1++)
                  {
                    if (a == b and J1 % 2 > 0)
                      continue;

                    for (int J2 = J2min; J2 <= J2max; J2++)
                    {
                      if (c == d and J2 % 2 > 0)
                        continue;
                      int twoJmin = std::max(std::abs(2 * J1 - oi.j2), std::abs(2 * J2 - oe.j2));
                      int twoJmax = std::min(2 * J1 + oi.j2, 2 * J2 + oe.j2);
                      for (int twoJ = twoJmin; twoJ <= twoJmax; twoJ += 2)
                      {
                        double xabicde = X3.GetME_pn(J1, J2, twoJ, a, b, i, c, d, e);
                        double yabicde = Y3.GetME_pn(J1, J2, twoJ, a, b, i, c, d, e);
                        double xcdeabj = X3.GetME_pn(J2, J1, twoJ, c, d, e, a, b, j);
                        double ycdeabj = Y3.GetME_pn(J2, J1, twoJ, c, d, e, a, b, j);
                        zij += 1. / 12 * occupation_factor * (twoJ + 1.) / (oi.j2 + 1) * (xabicde * ycdeabj - yabicde * xcdeabj);
                      }
                    } // Jj
                  } // J1
                } // e
              } // d
            } // c
          } // b
        } // a
        Z1(i, j) += zij;
      } // j
    } // i
  } // comm331ss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij = 1/4 sum_abcd sum_J2J na nb (1-nc) (1-nd) (2J+1)/(2ji+1) (XJ1_abcd YJ1J1J_cdiabj + XJ1J2J_abicdj YJ1_cdab    - YJ1_abcd XJ1J1J_cdiabj - YJ1J1J_abicdj XJ1_cdab)
  ///
  ///
  void comm231ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z1 = Z.OneBody;

    size_t norb = Z.modelspace->GetNumberOrbits();
    //   for ( size_t i : Z.modelspace->all_orbits )
    //   #pragma omp parallel for schedule(dynamic,1)
    for (size_t i = 0; i < norb; i++)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      //     for ( size_t j : Z.modelspace->all_orbits )
      for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double zij = 0;
        for (size_t a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          for (size_t b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);

            int J1min = std::abs(oa.j2 - ob.j2) / 2;
            int J1max = (oa.j2 + ob.j2) / 2;

            for (size_t c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              for (size_t d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                int J2min = std::abs(oc.j2 - od.j2) / 2;
                int J2max = (oc.j2 + od.j2) / 2;

                //                if ( (oa.l+ob.l+oi.l + oc.l+od.l+oj.l + Z.parity)%2 >0) continue;
                //                if ( std::abs(oa.tz2+ob.tz2+oi.tz2 - oc.tz2-od.tz2-oe.tz2) != Z.GetTRank()) continue;
                if (std::abs(oa.occ * ob.occ * (1 - oc.occ) * (1 - od.occ)) < 1e-8)
                  continue;

                for (int J1 = J1min; J1 <= J1max; J1++)
                {
                  if (a == b and J1 % 2 > 0)
                    continue;
                  if (c == d and J1 % 2 > 0)
                    continue;

                  double xabcd = X2.GetTBME_J(J1, J1, a, b, c, d);
                  double xcdab = X2.GetTBME_J(J1, J1, c, d, a, b);
                  double yabcd = Y2.GetTBME_J(J1, J1, a, b, c, d);
                  double ycdab = Y2.GetTBME_J(J1, J1, c, d, a, b);

                  int twoJmin = std::max(std::abs(2 * J1 - oi.j2), std::abs(2 * J1 - oj.j2));
                  int twoJmax = std::min(2 * J1 + oi.j2, 2 * J1 + oj.j2);
                  for (int twoJ = twoJmin; twoJ <= twoJmax; twoJ += 2)
                  {
                    double xabicdj = X3.GetME_pn(J1, J1, twoJ, a, b, i, c, d, j);
                    double yabicdj = Y3.GetME_pn(J1, J1, twoJ, a, b, i, c, d, j);
                    double xcdiabj = X3.GetME_pn(J1, J1, twoJ, c, d, i, a, b, j);
                    double ycdiabj = Y3.GetME_pn(J1, J1, twoJ, c, d, i, a, b, j);
                    zij += 1. / 4 * oa.occ * ob.occ * (1 - oc.occ) * (1 - od.occ) * (twoJ + 1.) / (oi.j2 + 1) * (xabcd * ycdiabj + xabicdj * ycdab - yabcd * xcdiabj - yabicdj * xcdab);
                  } // J
                } // J1
              } // d
            } // c
          } // b
        } // a
        Z1(i, j) += zij;
      } // j
    } // i
  } // comm231ss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1ijkl =  sum_ab sum_J (na-nb) sum_J  (2J+1)/(2J1+1) ( Xab YJ1J1J_ijbkla - Yab XJ1J1J_ijbkla )
  ///
  ///
  void comm132ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    std::vector<size_t> ch_bra_list, ch_ket_list;

    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++)
    //    int nch2 = Z.modelspace->GetNumberTwoBodyChannels();
    // #pragma omp parallel for schedule(dynamic, 1)
    //    for (int ch2 = 0; ch2 < nch2; ch2++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      //      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch2);
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
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
            Orbit &oa = Z.modelspace->GetOrbit(a);
            for (size_t b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double xab = X1(a, b);
              double yab = Y1(a, b);
              if (std::abs(xab) < 1e-8 and std::abs(yab) < 1e-8)
                continue;
              int twoJmin = std::max(std::abs(2 * J1 - oa.j2), std::abs(2 * J1 - ob.j2));
              int twoJmax = std::min(2 * J1 + oa.j2, 2 * J1 + ob.j2);
              for (int twoJ = twoJmin; twoJ <= twoJmax; twoJ += 2)
              {
                double xijbkla = X3.GetME_pn(J1, J1, twoJ, i, j, b, k, l, a);
                double yijbkla = Y3.GetME_pn(J1, J1, twoJ, i, j, b, k, l, a);
                zijkl += (oa.occ - ob.occ) * (twoJ + 1.) / (2 * J1 + 1) * (X1(a, b) * yijbkla - Y1(a, b) * xijbkla);

              } // twoJ

            } // b
          } // a
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
        } // iket
      } // ibra
    } // ch2

  } // comm132ss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1ijkl = -1/2 sum_abc sum_J2J3 ( nanb(1-nc) + (1-na)(1-nb)nc ) sqrt{ (2J2+1)/(2J1+1)} (2J3+1)
  ///                                 [ (1-PJ1_ij) { i j  J1 } ( XJ2_cjab YJ2J1J3_abiklc - YJ2_cjab XJ2J1J3_abiklc )
  ///                                              { c J3 J2 }
  ///                                  -(1-PJ1_kl) { k l  J1 } ( YJ1J23_ijcabk XJ2_abcl - XJ1J2J3_ijcabk YJ2_abcl  ]
  ///                                              { c J3 J2 }
  ///
  void comm232ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X2 = X.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y2 = Y.TwoBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    int nch = Z.modelspace->GetNumberTwoBodyChannels();

    std::vector<std::array<size_t, 2>> channels;
    for (auto &iter : Z.TwoBody.MatEl)
      channels.push_back(iter.first);
    size_t nchans = channels.size();
    //  for (int ch=0; ch<nch; ch++)
    //  #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar3b_transform_first_pass)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nchans; ich++)
    {
      size_t ch_bra = channels[ich][0];
      size_t ch_ket = channels[ich][1];
      auto &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      auto &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        int ket_min = (ch_bra == ch_ket) ? ibra : 0;
        for (int iket = ket_min; iket < nkets; iket++)
        {
          double zijkl = 0;
          Ket &ket = tbc_ket.GetKet(iket);
          int k = ket.p;
          int l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          double jk = 0.5 * ok.j2;
          double jl = 0.5 * ol.j2;
          for (auto c : Z.modelspace->all_orbits)
          {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double jc = 0.5 * oc.j2;

            for (int ch_ab = 0; ch_ab < nch; ch_ab++)
            {
              auto &tbc_ab = X.modelspace->GetTwoBodyChannel(ch_ab);
              int Jab = tbc_ab.J;
              size_t nkets_ab = tbc_ab.GetNumberKets();
              for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
              {
                Ket &ket_ab = tbc_ab.GetKet(iket_ab);
                int a = ket_ab.p;
                int b = ket_ab.q;
                Orbit &oa = Z.modelspace->GetOrbit(a);
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double occfactor = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
                if (std::abs(occfactor) < 1e-6)
                  continue;
                if (a == b)
                  occfactor *= 0.5; // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b

                bool xicab_good = ((oi.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(oi.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                bool yicab_good = ((oi.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(oi.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                bool xjcab_good = ((oj.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(oj.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                bool yjcab_good = ((oj.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(oj.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                bool xabck_good = ((ok.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(ok.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                bool yabck_good = ((ok.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(ok.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                bool xabcl_good = ((ol.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(ol.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                bool yabcl_good = ((ol.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(ol.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);

                // Xicab term
                if ((xicab_good or yicab_good) and (std::abs(oi.j2 - oc.j2) <= 2 * Jab) and (oi.j2 + oc.j2 >= 2 * Jab))
                {

                  int twoJ_min = std::max(std::abs(oc.j2 - 2 * J), std::abs(oj.j2 - 2 * Jab));
                  int twoJ_max = std::min(oc.j2 + 2 * J, oj.j2 + 2 * Jab);
                  double xciab = X2.GetTBME_J(Jab, c, i, a, b);
                  double yciab = Y2.GetTBME_J(Jab, c, i, a, b);
                  int phasefactor = Z.modelspace->phase((oi.j2 + oj.j2) / 2 - J);
                  for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                  {

                    double Jtot = 0.5 * twoJ;
                    double sixj = Z.modelspace->GetSixJ(jj, ji, J, jc, Jtot, Jab);
                    double hatfactor = (twoJ + 1) * sqrt((2 * Jab + 1.) / (2 * J + 1));
                    double xabjklc = yicab_good ? X3.GetME_pn(Jab, J, twoJ, a, b, j, k, l, c) : 0;
                    double yabjklc = xicab_good ? Y3.GetME_pn(Jab, J, twoJ, a, b, j, k, l, c) : 0;
                    zijkl += occfactor * hatfactor * phasefactor * sixj * (xciab * yabjklc - yciab * xabjklc);
                  }
                }

                // Xjcab term
                if ((xjcab_good or yjcab_good) and (std::abs(oj.j2 - oc.j2) <= 2 * Jab) and (oj.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoJ_min = std::max(std::abs(oc.j2 - 2 * J), std::abs(oi.j2 - 2 * Jab));
                  int twoJ_max = std::min(oc.j2 + 2 * J, oi.j2 + 2 * Jab);
                  double xcjab = X2.GetTBME_J(Jab, c, j, a, b);
                  double ycjab = Y2.GetTBME_J(Jab, c, j, a, b);
                  int phasefactor = 1;
                  for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                  {
                    double Jtot = 0.5 * twoJ;
                    double sixj = Z.modelspace->GetSixJ(ji, jj, J, jc, Jtot, Jab);
                    double hatfactor = (twoJ + 1) * sqrt((2 * Jab + 1.) / (2 * J + 1));
                    double xabiklc = yjcab_good ? X3.GetME_pn(Jab, J, twoJ, a, b, i, k, l, c) : 0;
                    double yabiklc = xjcab_good ? Y3.GetME_pn(Jab, J, twoJ, a, b, i, k, l, c) : 0;
                    zijkl -= occfactor * hatfactor * phasefactor * sixj * (xcjab * yabiklc - ycjab * xabiklc);
                  }
                }

                // Xabck term
                if ((xabck_good or yabck_good) and (std::abs(ok.j2 - oc.j2) <= 2 * Jab) and (ok.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoJ_min = std::max(std::abs(oc.j2 - 2 * J), std::abs(ol.j2 - 2 * Jab));
                  int twoJ_max = std::min(oc.j2 + 2 * J, ol.j2 + 2 * Jab);
                  double xabck = X2.GetTBME_J(Jab, a, b, c, k);
                  double yabck = Y2.GetTBME_J(Jab, a, b, c, k);
                  int phasefactor = Z.modelspace->phase((ok.j2 + ol.j2) / 2 - J);
                  for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                  {
                    double Jtot = 0.5 * twoJ;
                    double sixj = Z.modelspace->GetSixJ(jl, jk, J, jc, Jtot, Jab);
                    double hatfactor = (twoJ + 1) * sqrt((2 * Jab + 1.) / (2 * J + 1));
                    double xijcabl = yabck_good ? X3.GetME_pn(J, Jab, twoJ, i, j, c, a, b, l) : 0;
                    double yijcabl = xabck_good ? Y3.GetME_pn(J, Jab, twoJ, i, j, c, a, b, l) : 0;
                    zijkl -= occfactor * hatfactor * phasefactor * sixj * (yijcabl * xabck - xijcabl * yabck);
                  }
                }

                // Xabcl term
                if ((xabcl_good or yabcl_good) and (std::abs(ol.j2 - oc.j2) <= 2 * Jab) and (ol.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoJ_min = std::max(std::abs(oc.j2 - 2 * J), std::abs(ok.j2 - 2 * Jab));
                  int twoJ_max = std::min(oc.j2 + 2 * J, ok.j2 + 2 * Jab);
                  double xabcl = X2.GetTBME_J(Jab, a, b, c, l);
                  double yabcl = Y2.GetTBME_J(Jab, a, b, c, l);
                  int phasefactor = 1;
                  for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
                  {
                    double Jtot = 0.5 * twoJ;
                    double sixj = Z.modelspace->GetSixJ(jk, jl, J, jc, Jtot, Jab);
                    double hatfactor = (twoJ + 1) * sqrt((2 * Jab + 1.) / (2 * J + 1));
                    double xijcabk = yabcl_good ? X3.GetME_pn(J, Jab, twoJ, i, j, c, a, b, k) : 0;
                    double yijcabk = xabcl_good ? Y3.GetME_pn(J, Jab, twoJ, i, j, c, a, b, k) : 0;
                    zijkl += occfactor * hatfactor * phasefactor * sixj * (yijcabk * xabcl - xijcabk * yabcl);
                  }
                }

              } // for iket_ab
            } // for ch2
          } // for c

          // normalize the tbme
          zijkl *= -1.0 / sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
        } // for iket
      } // for ibra
    } // for ch

  } // comm232ss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1ijkl = 1/6 sum_abcd sum_J2J3  (na nb nc(1-nd) - (1-na)(1-nb)(1-nc)nd) (2J3+1)/(2J1+1)
  ///                                       *  ( XJ1J2J3_ijdabc YJ2J1J3_abckld  - YJ1J2J3_ijdabc XJ2J1J3_abckld )
  ///
  void comm332_ppph_hhhpss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;

    int nch2 = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch2 = 0; ch2 < nch2; ch2++)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch2);
      int J1 = tbc.J;
      int nkets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nkets; ibra++)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        for (int iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;
          double zijkl = 0;
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

                  if (std::abs(oa.occ * ob.occ * oc.occ * (1 - od.occ) - (1 - oa.occ) * (1 - ob.occ) * (1 - oc.occ) * od.occ) < 1e-7)
                    continue;
                  if ((oi.l + oj.l + od.l + oa.l + ob.l + oc.l) % 2 > 0)
                    continue;
                  if ((oi.tz2 + oj.tz2 + od.tz2) - (oa.tz2 + ob.tz2 + oc.tz2) != 0)
                    continue;

                  int J2min = std::abs(oa.j2 - ob.j2) / 2;
                  int J2max = (oa.j2 + ob.j2) / 2;
                  for (int J2 = J2min; J2 <= J2max; J2++)
                  {

                    int twoJmin = std::max(std::abs(2 * J1 - od.j2), std::abs(2 * J2 - oc.j2));
                    int twoJmax = std::min(2 * J1 + od.j2, 2 * J2 + oc.j2);
                    for (int twoJ = twoJmin; twoJ <= twoJmax; twoJ += 2)
                    {
                      double xijdabc = X3.GetME_pn(J1, J2, twoJ, i, j, d, a, b, c);
                      double yijdabc = Y3.GetME_pn(J1, J2, twoJ, i, j, d, a, b, c);
                      double xabckld = X3.GetME_pn(J2, J1, twoJ, a, b, c, k, l, d);
                      double yabckld = Y3.GetME_pn(J2, J1, twoJ, a, b, c, k, l, d);
                      zijkl += 1. / 6 * (oa.occ * ob.occ * oc.occ * (1 - od.occ) - (1 - oa.occ) * (1 - ob.occ) * (1 - oc.occ) * od.occ) * (twoJ + 1.) / (2 * J1 + 1) * (xijdabc * yabckld - yijdabc * xabckld);
                    } // twoJ
                  } // J2
                } // d
              } // c
            } // b
          } // a
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          Z2.AddToTBME(ch2, ch2, ibra, iket, zijkl);
        } // iket
      } // ibra
    } // ch2

  } // comm332_ppph_hhhpss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ_ijkn = (1-PJij)(1-PJkl) -sum_Jph (2Jph+1) { i j J   } ZbarJph_ilkj
  ///                                                             { k l Jph }
  ///
  ///  where         ZbarJph_ilkj = 1/4 sum_abcd sum Jab Jcd  ( na nb(1-nc)(1-nd) -  (1-na)(1-nb)nc nd ) XJph_il;(abJab)(cdJcd) YJph_(cdJcd)(abJab);kj
  ///
  ///  and           XJph_il;(abJab)(cdJcd) = -sum_J (-1)^{i+J} {i   l   Jph } XJabJcdJ_abicdl
  ///                                                           {Jcd Jab J   }
  ///
  void comm332_pphhss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    int nch = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (int ch = 0; ch < nch; ch++)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nkets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nkets; ibra++)
      {
        Ket &bra = tbc.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        int ji2 = bra.op->j2;
        int jj2 = bra.oq->j2;
        double ji = 0.5 * ji2;
        double jj = 0.5 * jj2;
        double d_ei = std::abs(2 * bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
        double d_ej = std::abs(2 * bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
        double occnat_i = bra.op->occ_nat;
        double occnat_j = bra.oq->occ_nat;
        for (int iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc.GetKet(iket);
          int k = ket.p;
          int l = ket.q;
          int jk2 = ket.op->j2;
          int jl2 = ket.oq->j2;
          double jk = 0.5 * jk2;
          double jl = 0.5 * jl2;
          double d_ek = std::abs(2 * ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
          double d_el = std::abs(2 * ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
          double occnat_k = ket.op->occ_nat;
          double occnat_l = ket.oq->occ_nat;

          double zijkl = 0;
          // Now the loops on the right hand side
          for (int ch_ab = 0; ch_ab < nch; ch_ab++)
          {
            TwoBodyChannel &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
            int Jab = tbc_ab.J;
            int nkets = tbc_ab.GetNumberKets();
            for (int iket_ab = 0; iket_ab < nkets; iket_ab++)
            {
              Ket &ket_ab = tbc_ab.GetKet(iket_ab);
              int a = ket_ab.p;
              int b = ket_ab.q;
              double na = ket_ab.op->occ;
              double nb = ket_ab.oq->occ;

              for (int ch_cd = 0; ch_cd < nch; ch_cd++)
              {
                TwoBodyChannel &tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
                int Jcd = tbc_cd.J;
                int nkets = tbc_cd.GetNumberKets();
                for (int iket_cd = 0; iket_cd < nkets; iket_cd++)
                {
                  Ket &ket_cd = tbc_cd.GetKet(iket_cd);
                  int c = ket_cd.p;
                  int d = ket_cd.q;
                  double nc = ket_cd.op->occ;
                  double nd = ket_cd.oq->occ;
                  //                double jc = 0.5*ket_cd.op->j2;
                  //                double jd = 0.5*ket_cd.oq->j2;
                  double occupation_factor = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);
                  if (std::abs(occupation_factor) < 1e-6)
                    continue;

                  double symmetry_factor = 1; // we only sum a<=b and c<=d, so we undercount by a factor of 4, canceling the 1/4 in the formula
                  if (a == b)
                    symmetry_factor *= 0.5; // if a==b or c==d, then the permutation doesn't give a new state, so there's less undercounting
                  if (c == d)
                    symmetry_factor *= 0.5;

                  // Figure out which range of twoJp and twoJpp we will need
                  int twoJp_min = std::max(std::min(std::abs(2 * Jab - ji2), std::abs(2 * Jab - jj2)), std::min(std::abs(2 * Jcd - jl2), std::abs(2 * Jcd - jk2)));
                  int twoJpp_min = std::max(std::min(std::abs(2 * Jcd - jj2), std::abs(2 * Jcd - ji2)), std::min(std::abs(2 * Jab - jk2), std::abs(2 * Jab - jl2)));
                  int twoJp_max = std::min(2 * Jab + std::max(ji2, jj2), 2 * Jcd + std::max(jk2, jl2));
                  int twoJpp_max = std::min(2 * Jcd + std::max(ji2, jj2), 2 * Jab + std::max(jk2, jl2));

                  if (twoJpp_max < twoJpp_min)
                    continue;
                  for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                  {

                    double xabicdl = X3.GetME_pn(Jab, Jcd, twoJp, a, b, i, c, d, l);
                    double xabicdk = X3.GetME_pn(Jab, Jcd, twoJp, a, b, i, c, d, k);
                    double xabjcdl = X3.GetME_pn(Jab, Jcd, twoJp, a, b, j, c, d, l);
                    double xabjcdk = X3.GetME_pn(Jab, Jcd, twoJp, a, b, j, c, d, k);

                    for (int twoJpp = twoJpp_min; twoJpp <= twoJpp_max; twoJpp += 2)
                    {
                      double Jp = 0.5 * twoJp;
                      double Jpp = 0.5 * twoJpp;
                      double hatfactor = (twoJp + 1) * (twoJpp + 1);

                      // I think having these in the inner loop may be disastrous for performance
                      double ninej1 = Z.modelspace->GetNineJ(Jab, jk, Jpp, ji, J, jj, Jp, jl, Jcd);
                      double ninej2 = Z.modelspace->GetNineJ(Jab, jk, Jpp, jj, J, ji, Jp, jl, Jcd); // permute i<->j
                      double ninej3 = Z.modelspace->GetNineJ(Jab, jl, Jpp, ji, J, jj, Jp, jk, Jcd); // permute k<->l
                      double ninej4 = Z.modelspace->GetNineJ(Jab, jl, Jpp, jj, J, ji, Jp, jk, Jcd); // permute i<->j and k<->l

                      // These phase factors account for the minus signs associated with the permutations
                      // so that all the permuted terms should just be added with their phase (no extra minus sign).
                      int phase1 = Z.modelspace->phase((ji2 + jk2 + twoJp + twoJpp) / 2);
                      int phase2 = Z.modelspace->phase((ji2 + jk2 + twoJp + twoJpp) / 2 - J); // from permuting i<->j
                      int phase3 = Z.modelspace->phase((ji2 + jk2 + twoJp + twoJpp) / 2 - J); // from permuting k<->l
                      int phase4 = Z.modelspace->phase((ji2 + jk2 + twoJp + twoJpp) / 2);     // from permuting i<->j and k<->l

                      double ycdjabk = Y3.GetME_pn(Jcd, Jab, twoJpp, c, d, j, a, b, k);
                      double ycdjabl = Y3.GetME_pn(Jcd, Jab, twoJpp, c, d, j, a, b, l);
                      double ycdiabk = Y3.GetME_pn(Jcd, Jab, twoJpp, c, d, i, a, b, k);
                      double ycdiabl = Y3.GetME_pn(Jcd, Jab, twoJpp, c, d, i, a, b, l);

                      zijkl += symmetry_factor * hatfactor * occupation_factor * (1 * phase1 * ninej1 * xabicdl * ycdjabk + 1 * phase2 * ninej2 * xabjcdl * ycdiabk // i<->j
                                                                                  + 1 * phase3 * ninej3 * xabicdk * ycdjabl                                         // k<->l
                                                                                  + 1 * phase4 * ninej4 * xabjcdk * ycdiabl                                         // i<->j and k<->l
                                                                                 );

                    } // for twoJpp
                  } // for twoJp

                } // for iket_cd
              } // for ch_cd
            } // for iket_ab
          } // for ch_ab
          // make it a normalized TBME
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          // the AddToTBME routine automatically takes care of the hermitian conjugate as well
          Z2.AddToTBME(ch, ch, ibra, iket, zijkl);
        } // for iket
      } // for ibra
    } // for ch
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1J2J3_ijklmn =  sum_a  PJ1J3(ij/k) ( X_ka YJ1J2J3_ijalmn - Y_ka XJ1J2J3_ijalmn)
  ///                                       - PJ2J3(lm/n) ( YJ1J2J3_ijklma X_an - XJ1J2J3_ijklma Y_an )
  ///
  ///
  void comm133ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    int norbs = Z.modelspace->GetNumberOrbits();
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
    Z.modelspace->PreCalculateSixJ();

    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, static_cast<size_t>(ibra)}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < n_bra_ket_ch; ich++)
    {
      size_t ch_bra = bra_ket_channels[ich][0];
      size_t ch_ket = bra_ket_channels[ich][1];
      size_t ibra = bra_ket_channels[ich][2];
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch_bra);
      ThreeBodyChannel &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch_ket);

      int twoJ = Tbc_bra.twoJ;
      size_t nbras = Tbc_bra.GetNumberKets();
      size_t nkets = Tbc_ket.GetNumberKets();
      Ket3 &bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      Orbit &ok = Z.modelspace->GetOrbit(k);
      int Jij = bra.Jpq;
      size_t ket_min = (ch_bra == ch_ket) ? ibra : 0;
      for (size_t iket = ket_min; iket < nkets; iket++)
      {
        Ket3 &ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit &ol = Z.modelspace->GetOrbit(l);
        Orbit &om = Z.modelspace->GetOrbit(m);
        Orbit &on = Z.modelspace->GetOrbit(n);
        int Jlm = ket.Jpq;

        double zsum = 0;
        // First, connect on the bra side
        
        for (auto a : X.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
        {
          zsum += X1(i, a) * Y3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
        }
        for (auto a : Y.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
        {
          zsum -= Y1(i, a) * X3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
        }
        for (auto a : X.GetOneBodyChannel(oj.l, oj.j2, oj.tz2))
        {
          zsum += X1(j, a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
        }
        for (auto a : Y.GetOneBodyChannel(oj.l, oj.j2, oj.tz2))
        {
          zsum -= Y1(j, a) * X3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
        }
        for (auto a : X.GetOneBodyChannel(ok.l, ok.j2, ok.tz2))
        {
          zsum += X1(k, a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
        }
        for (auto a : Y.GetOneBodyChannel(ok.l, ok.j2, ok.tz2))
        {
          zsum -= Y1(k, a) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
        }

        // Now connect on the ket side
        for (auto a : X.GetOneBodyChannel(ol.l, ol.j2, ol.tz2))
        {
          zsum -= X1(a, l) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
        }
        for (auto a : Y.GetOneBodyChannel(ol.l, ol.j2, ol.tz2))
        {
          zsum += Y1(a, l) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
        }
        for (auto a : X.GetOneBodyChannel(om.l, om.j2, om.tz2))
        {
          zsum -= X1(a, m) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
        }
        for (auto a : Y.GetOneBodyChannel(om.l, om.j2, om.tz2))
        {
          zsum += Y1(a, m) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
        }       
        for (auto a : X.GetOneBodyChannel(on.l, on.j2, on.tz2))
        {
          zsum -= X1(a, n) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
        }
        for (auto a : Y.GetOneBodyChannel(on.l, on.j2, on.tz2))
        {
          zsum += Y1(a, n) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
        }

        Z3.AddToME_pn_ch(ch_bra, ch_ket, ibra, iket, zsum);

      } // for iket
    } // for ich  -> {ch_bra,ch_ket,ibra}

  } // comm133ss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1J2J3_ijklmn =  sum_a  PJ1J3(ij/k) PJ1J2(lm/n) sqrt( (2J1+1)(2J2+1)) { n  a  J1 } ( XJ1_ijna YJ2_kalm - YJ1_ijna XJ2_kalm )
  ///                                                                                       { k  J3 J2 }
  ///
  ///
  void comm223ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &Z3 = Z.ThreeBody;
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    Z.modelspace->PreCalculateSixJ();
    // Permutations of indices which are needed to produce antisymmetrized matrix elements  P(ij/k) |ijk> = |ijk> - |kji> - |ikj>
    const std::array<ThreeBodyStorage::Permutation, 3> index_perms = {ThreeBodyStorage::ABC, ThreeBodyStorage::CBA, ThreeBodyStorage::ACB};

    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, static_cast<size_t>(ibra)}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();


#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ibra_ket = 0; ibra_ket < n_bra_ket_ch; ibra_ket++)
    {
      size_t ch3bra = bra_ket_channels[ibra_ket][0];
      size_t ch3ket = bra_ket_channels[ibra_ket][1];
      size_t ibra = bra_ket_channels[ibra_ket][2];
      auto &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch3bra);
      auto &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch3ket);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      size_t nkets3 = Tbc_ket.GetNumberKets();
      int twoJ = Tbc_bra.twoJ; // Scalar commutator so J is the same in bra and ket channel
      double Jtot = 0.5 * twoJ;

      auto &bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      Orbit &ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;

      int J1 = bra.Jpq;

      size_t iket_max = nkets3-1;
      if (ch3bra == ch3ket)
        iket_max = ibra;
      for (size_t iket = 0; iket <= iket_max; iket++)
      {
        auto &ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit &ol = Z.modelspace->GetOrbit(l);
        Orbit &om = Z.modelspace->GetOrbit(m);
        Orbit &on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        int J2 = ket.Jpq;

        double zijklmn = 0;
        /// BEGIN THE SLOW BIT...

        // Now we need to loop over the permutations in ijk and then lmn
        for (auto perm_ijk : index_perms)
        {
          size_t I1, I2, I3;
          Z3.Permute(perm_ijk, i, j, k, I1, I2, I3);
          Orbit &o1 = Z.modelspace->GetOrbit(I1);
          Orbit &o2 = Z.modelspace->GetOrbit(I2);
          Orbit &o3 = Z.modelspace->GetOrbit(I3);


          int J1p_min = J1;
          int J1p_max = J1;
          if (perm_ijk != ThreeBodyStorage::ABC)
          {
            J1p_min = std::max(std::abs(o1.j2 - o2.j2), std::abs(twoJ - o3.j2)) / 2;
            J1p_max = std::min(o1.j2 + o2.j2, twoJ + o3.j2) / 2;
          }

          int parity_12 = (o1.l + o2.l) % 2;
          int Tz_12 = (o1.tz2 + o2.tz2) / 2;

          double j3 = 0.5 * o3.j2;

          for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
          {

            double rec_ijk = Z3.RecouplingCoefficient(perm_ijk, ji, jj, jk, J1p, J1, twoJ);
            rec_ijk *= Z3.PermutationPhase(perm_ijk); // do we get a fermionic minus sign?


            for (auto perm_lmn : index_perms)
            {
              size_t I4, I5, I6;
              Z3.Permute(perm_lmn, l, m, n, I4, I5, I6);
              Orbit &o4 = Z.modelspace->GetOrbit(I4);
              Orbit &o5 = Z.modelspace->GetOrbit(I5);
              Orbit &o6 = Z.modelspace->GetOrbit(I6);

              double j6 = 0.5 * o6.j2;

              int J2p_min = J2;
              int J2p_max = J2;
              if (perm_lmn != ThreeBodyStorage::ABC)
              {
                J2p_min = std::max(std::abs(o4.j2 - o5.j2), std::abs(twoJ - o6.j2)) / 2;
                J2p_max = std::min(o4.j2 + o5.j2, twoJ + o6.j2) / 2;
              }


              for (int tz2a : {-1, 1})
              {
                int dTz_126a = (o1.tz2 + o2.tz2 - o6.tz2 - tz2a)/2;
                int dTz_3a45 = (o3.tz2 + tz2a - o4.tz2 - o5.tz2 )/2;
                if (std::abs(dTz_126a) != X.GetTRank() and std::abs(dTz_126a) != Y.GetTRank())  continue;
                if (std::abs(dTz_3a45) != X.GetTRank() and std::abs(dTz_3a45) != Y.GetTRank())  continue;

                for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                {

                  double rec_lmn = Z3.RecouplingCoefficient(perm_lmn, jl, jm, jn, J2p, J2, twoJ);
                  rec_lmn *= Z3.PermutationPhase(perm_lmn); // do we get a fermionic minus sign?

                  int j2a_min = std::max(std::abs(o6.j2 - 2 * J1p), std::abs(o3.j2 - 2 * J2p));
                  int j2a_max = std::min(o6.j2 + 2 * J1p, o3.j2 + 2 * J2p);

                  for (int j2a = j2a_min; j2a <= j2a_max; j2a += 2)
                  {
                    double ja = 0.5 * j2a;

                    double sixj;
                    if (twoJ <= 2 * Z.modelspace->GetEmax() + 1)
                    {
                      sixj = Z.modelspace->GetCachedSixJ(o3.j2, twoJ, J1p, o6.j2, j2a, J2p);
                    }
                    else
                    {
                      sixj = j2a < twoJ ? Z.modelspace->GetSixJ(j6, ja, J1p, j3, 0.5 * twoJ, J2p)
                                        : Z.modelspace->GetSixJ(j6, 0.5 * twoJ, J2p, j3, ja, J1p);
                    }

                    for (size_t a : Z.modelspace->all_orbits)
                    {
                      Orbit& oa = Z.modelspace->GetOrbit(a);
                      if ( oa.j2 != j2a ) continue;
                      if ( oa.tz2 != tz2a ) continue;
                      double x_126a = X2.GetTBME_J(J1p, J1p, I1, I2, I6, a);
                      double y_126a = Y2.GetTBME_J(J1p, J1p, I1, I2, I6, a);
                      double x_3a45 = X2.GetTBME_J(J2p, J2p, I3, a, I4, I5);
                      double y_3a45 = Y2.GetTBME_J(J2p, J2p, I3, a, I4, I5);

                      zijklmn += rec_ijk * rec_lmn * sixj * sqrt((2 * J1p + 1) * (2 * J2p + 1)) * (x_126a * y_3a45 - y_126a * x_3a45);
                     
                    } // for a
                  }   // for j2a
                }     // for J2p
              }       // for tz2a
            }         // for perm_lmn
          }           // for J1p
        }             // for perm_ijk

        Z3.AddToME_pn_ch(ch3bra, ch3ket, ibra, iket, zijklmn); // this needs to be modified for beta decay
      }                                                        // for iket
    }                                                          // for ch3

  } // comm233ss

  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1J2J3_ijklmn =  1/2 sum_ab ((1-na)(1-nb)-na nb) [ PJ1J3(ij/k) ( XJ1_ijab YJ1J2J3_abklmn - YJ1_ijab XJ1J2J3_abklmn )
  ///                                                                   -PJ1J2(lm/n) ( YJ1J2J3_ijkabn XJ2_ablm - XJ1J2J3_ijkabn YJ2_ablm )
  ///
  ///
  /// THIS VERSION IS STILL TOO SLOW FOR GOING BEYOND EMAX=2...
  void comm233_pp_hhss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;

    // Permutations of indices which are needed to produce antisymmetrized matrix elements  P(ij/k) |ijk> = |ijk> - |kji> - |ikj>
    const std::array<ThreeBodyStorage::Permutation, 3> index_perms = {ThreeBodyStorage::ABC, ThreeBodyStorage::CBA, ThreeBodyStorage::ACB};

    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, static_cast<size_t>(ibra)}); // (ch_bra, ch_ket,ibra)
      }
    }

    size_t n_bra_ket_ch = bra_ket_channels.size();

    //  std::vector< std::array<size_t,2> > bra_ket_channels;
    //  for ( auto& it : Z.ThreeBody.Get_ch_start() )
    //  {
    //     bra_ket_channels.push_back( { it.first[0],it.first[1] } ); // (ch_bra, ch_ket)
    //  }

    //  size_t n_bra_ket_ch = bra_ket_channels.size();
    //  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
    //  for (size_t ch3=0; ch3<nch3; ch3++)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ibra_ket = 0; ibra_ket < n_bra_ket_ch; ibra_ket++)
    {
      size_t ch3bra = bra_ket_channels[ibra_ket][0];
      size_t ch3ket = bra_ket_channels[ibra_ket][1];
      size_t ibra = bra_ket_channels[ibra_ket][2];
      auto &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch3bra);
      auto &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch3ket);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      size_t nkets3 = Tbc_ket.GetNumberKets();
      int twoJ = Tbc_bra.twoJ; // Scalar commutator so J is the same in bra and ket channel
      double Jtot = 0.5 * twoJ;

      //    for (size_t ibra=0; ibra<nbras3; ibra++)
      //    {
      auto &bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      Orbit &ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;

      int J1 = bra.Jpq;

      size_t iket_max = nkets3;
      if (ch3bra == ch3ket)
        iket_max = ibra;
      for (size_t iket = 0; iket <= iket_max; iket++)
      {
        auto &ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit &ol = Z.modelspace->GetOrbit(l);
        Orbit &om = Z.modelspace->GetOrbit(m);
        Orbit &on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        int J2 = ket.Jpq;

        double zijklmn = 0;

        // Now we need to loop over the permutations in ijk and then lmn
        for (auto perm_ijk : index_perms) // {ijk} -> {123}
        {
          size_t I1, I2, I3;
          Z3.Permute(perm_ijk, i, j, k, I1, I2, I3);
          Orbit &o1 = Z.modelspace->GetOrbit(I1);
          Orbit &o2 = Z.modelspace->GetOrbit(I2);
          Orbit &o3 = Z.modelspace->GetOrbit(I3);

          int J1p_min = J1;
          int J1p_max = J1;
          if (perm_ijk != ThreeBodyStorage::ABC)
          {
            J1p_min = std::max(std::abs(o1.j2 - o2.j2), std::abs(twoJ - o3.j2)) / 2;
            J1p_max = std::min(o1.j2 + o2.j2, twoJ + o3.j2) / 2;
          }

          for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
          {
            double Pijk = Z3.PermutationPhase(perm_ijk) * Z3.RecouplingCoefficient(perm_ijk, ji, jj, jk, J1p, J1, twoJ);

            for (size_t a : Z.modelspace->all_orbits)
            {
              Orbit &oa = Z.modelspace->GetOrbit(a);
              for (size_t b : Z.modelspace->all_orbits)
              {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                if (J1p < std::abs(oa.j2 - ob.j2) / 2 or J1p > (oa.j2 + ob.j2) / 2)
                  continue;
                if ((o1.l + o2.l + oa.l + ob.l + X.parity) % 2 > 0 and (o1.l + o2.l + oa.l + ob.l + Y.parity) % 2 > 0)
                  continue;
                if ((std::abs(o1.tz2 + o2.tz2 - oa.tz2 - ob.tz2) != 2 * X.GetTRank()) and (std::abs(o1.tz2 + o2.tz2 - oa.tz2 - ob.tz2) != 2 * Y.GetTRank()))
                  continue;
                if (std::abs(((1 - oa.occ) * (1 - ob.occ) - oa.occ * ob.occ)) < 1e-7)
                  continue;

                double x12ab = X2.GetTBME_J(J1p, J1p, I1, I2, a, b);
                double y12ab = Y2.GetTBME_J(J1p, J1p, I1, I2, a, b);
                auto x_y = X3.GetME_pn_TwoOps(J1p, J2, twoJ, a, b, I3, l, m, n, X3, Y3);
                double xab3lmn = x_y[0];
                double yab3lmn = x_y[1];
                //                    double xab3lmn = X3.GetME_pn( J1p, J2, twoJ, a,b,I3,l,m,n );
                //                    double yab3lmn = Y3.GetME_pn( J1p, J2, twoJ, a,b,I3,l,m,n );

                zijklmn += 1. / 2 * ((1 - oa.occ) * (1 - ob.occ) - oa.occ * ob.occ) * Pijk * (x12ab * yab3lmn - y12ab * xab3lmn);
              } // b
            } // a
          }
        } // for perm_ijk

        for (auto perm_lmn : index_perms) // {lmn} -> {123}
        {
          size_t I1, I2, I3;
          Z3.Permute(perm_lmn, l, m, n, I1, I2, I3);
          Orbit &o1 = Z.modelspace->GetOrbit(I1);
          Orbit &o2 = Z.modelspace->GetOrbit(I2);
          Orbit &o3 = Z.modelspace->GetOrbit(I3);

          int J2p_min = J2;
          int J2p_max = J2;
          if (perm_lmn != ThreeBodyStorage::ABC)
          {
            J2p_min = std::max(std::abs(o1.j2 - o2.j2), std::abs(twoJ - o3.j2)) / 2;
            J2p_max = std::min(o1.j2 + o2.j2, twoJ + o3.j2) / 2;
          }
          for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
          {
            double Plmn = Z3.PermutationPhase(perm_lmn) * Z3.RecouplingCoefficient(perm_lmn, jl, jm, jn, J2p, J2, twoJ);

            for (size_t a : Z.modelspace->all_orbits)
            {
              Orbit &oa = Z.modelspace->GetOrbit(a);
              for (size_t b : Z.modelspace->all_orbits)
              {
                Orbit &ob = Z.modelspace->GetOrbit(b);

                if (J2p < std::abs(oa.j2 - ob.j2) / 2 or J2p > (oa.j2 + ob.j2) / 2)
                  continue;
                if ((o1.l + o2.l + oa.l + ob.l + X.parity) % 2 > 0 and (o1.l + o2.l + oa.l + ob.l + Y.parity) % 2 > 0)
                  continue;
                if ((std::abs(o1.tz2 + o2.tz2 - oa.tz2 - ob.tz2) != 2 * X.GetTRank()) and (std::abs(o1.tz2 + o2.tz2 - oa.tz2 - ob.tz2) != 2 * Y.GetTRank()))
                  continue;
                if (std::abs(((1 - oa.occ) * (1 - ob.occ) - oa.occ * ob.occ)) < 1e-7)
                  continue;

                double xab12 = X2.GetTBME_J(J2p, J2p, a, b, I1, I2);
                double yab12 = Y2.GetTBME_J(J2p, J2p, a, b, I1, I2);

                auto x_y = X3.GetME_pn_TwoOps(J1, J2p, twoJ, i, j, k, a, b, I3, X3, Y3);
                double xijkab3 = x_y[0];
                double yijkab3 = x_y[1];
                //                    double xijkab3 = X3.GetME_pn( J1, J2p, twoJ, i,j,k, a,b,I3 );
                //                    double yijkab3 = Y3.GetME_pn( J1, J2p, twoJ, i,j,k, a,b,I3 );
                zijklmn -= 1. / 2 * ((1 - oa.occ) * (1 - ob.occ) - oa.occ * ob.occ) * Plmn * (yijkab3 * xab12 - xijkab3 * yab12);
              } // b
            } // a
          }
        } // for perm_lmn

        //          }//b
        //        }//a

        Z3.AddToME_pn_ch(ch3bra, ch3ket, ibra, iket, zijklmn); // this needs to be modified for beta decay
      } // for iket
        //    }// for ibra
    } // for ch3

  } // comm233_pp_hhss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1J2J3_ijklmn =  - sum_ab sum_J4J5 (na-nb) PJ1J3(ij/k) PJ2J3(lm/n) (-1)^{k+n+J1+J2} { b  J5 J2 } ( XJ4_bkan YJ1J2J5_ijalmb - YJ4_bkan XJ1J2J5_ijalmb )
  ///                                                                                                     { k  J1 J3 }
  ///                                                                                                     { J4 a  n  }
  ///
  void comm233_phss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;

    // Permutations of indices which are needed to produce antisymmetrized matrix elements  P(ij/k) |ijk> = |ijk> - |kji> - |ikj>
    const std::array<ThreeBodyStorage::Permutation, 3> index_perms = {ThreeBodyStorage::ABC, ThreeBodyStorage::CBA, ThreeBodyStorage::ACB};

    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, static_cast<size_t>(ibra)}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();

    //  std::vector< std::array<size_t,2> > bra_ket_channels;
    //  for ( auto& it : Z.ThreeBody.Get_ch_start() )
    //  {
    //     bra_ket_channels.push_back( { it.first[0],it.first[1] } ); // (ch_bra, ch_ket)
    //  }

    //  Z.modelspace->PreCalculateNineJ();

    //  size_t n_bra_ket_ch = bra_ket_channels.size();
    //  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
    //  for (size_t ch3=0; ch3<nch3; ch3++)
    //  #pragma omp parallel for schedule(dynamic,1)  /// The 9js used here aren't generated by PreCalculateNineJ, so the parallel loop causes trouble.
    for (size_t ibra_ket = 0; ibra_ket < n_bra_ket_ch; ibra_ket++)
    {
      size_t ch3bra = bra_ket_channels[ibra_ket][0];
      size_t ch3ket = bra_ket_channels[ibra_ket][1];
      size_t ibra = bra_ket_channels[ibra_ket][2];
      auto &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch3bra);
      auto &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch3ket);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      size_t nkets3 = Tbc_ket.GetNumberKets();
      int twoJ = Tbc_bra.twoJ; // Scalar commutator so J is the same in bra and ket channel
      double Jtot = 0.5 * twoJ;

      //    for (size_t ibra=0; ibra<nbras3; ibra++)
      //    {
      auto &bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      Orbit &ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;

      int J1 = bra.Jpq;
      size_t iket_max = nkets3;
      if (ch3bra == ch3ket)
        iket_max = ibra;
      for (size_t iket = 0; iket <= iket_max; iket++)
      {
        auto &ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit &ol = Z.modelspace->GetOrbit(l);
        Orbit &om = Z.modelspace->GetOrbit(m);
        Orbit &on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        int J2 = ket.Jpq;

        double zijklmn = 0;

        // Now we need to loop over the permutations in ijk and then lmn
        for (auto perm_ijk : index_perms) // {ijk} -> {123}
        {
          //          if (perm_ijk != index_perms[0]) continue;
          size_t I1, I2, I3;
          Z3.Permute(perm_ijk, i, j, k, I1, I2, I3);
          Orbit &o1 = Z.modelspace->GetOrbit(I1);
          Orbit &o2 = Z.modelspace->GetOrbit(I2);
          Orbit &o3 = Z.modelspace->GetOrbit(I3);

          int J1p_min = J1;
          int J1p_max = J1;
          if (perm_ijk != ThreeBodyStorage::ABC)
          {
            J1p_min = std::max(std::abs(o1.j2 - o2.j2), std::abs(twoJ - o3.j2)) / 2;
            J1p_max = std::min(o1.j2 + o2.j2, twoJ + o3.j2) / 2;
          }

          for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
          {
            double Pijk = Z3.PermutationPhase(perm_ijk) * Z3.RecouplingCoefficient(perm_ijk, ji, jj, jk, J1p, J1, twoJ);

            for (auto perm_lmn : index_perms) // {lmn} -> {456}
            {
              //               if (perm_lmn != index_perms[0]) continue;
              size_t I4, I5, I6;
              Z3.Permute(perm_lmn, l, m, n, I4, I5, I6);
              Orbit &o4 = Z.modelspace->GetOrbit(I4);
              Orbit &o5 = Z.modelspace->GetOrbit(I5);
              Orbit &o6 = Z.modelspace->GetOrbit(I6);

              int J2p_min = J2;
              int J2p_max = J2;
              if (perm_lmn != ThreeBodyStorage::ABC)
              {
                J2p_min = std::max(std::abs(o4.j2 - o5.j2), std::abs(twoJ - o6.j2)) / 2;
                J2p_max = std::min(o4.j2 + o5.j2, twoJ + o6.j2) / 2;
              }
              for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
              {
                double Plmn = Z3.PermutationPhase(perm_lmn) * Z3.RecouplingCoefficient(perm_lmn, jl, jm, jn, J2p, J2, twoJ);

                for (size_t a : Z.modelspace->all_orbits)
                {
                  Orbit &oa = Z.modelspace->GetOrbit(a);
                  for (size_t b : Z.modelspace->all_orbits)
                  {
                    Orbit &ob = Z.modelspace->GetOrbit(b);

                    if (std::abs(oa.occ - ob.occ) < 1e-7)
                      continue;

                    if ((o3.l + ob.l + o6.l + oa.l + X.parity) % 2 > 0 and (o3.l + ob.l + o6.l + oa.l + Y.parity) % 2 > 0)
                      continue;
                    if (std::abs(o3.tz2 + ob.tz2 - o6.tz2 - oa.tz2) != X.GetTRank() and std::abs(o3.tz2 + ob.tz2 - o6.tz2 - oa.tz2) != Y.GetTRank())
                      continue;

                    int Ja6_min = std::max(std::abs(oa.j2 - o6.j2), std::abs(ob.j2 - o3.j2)) / 2;
                    int Ja6_max = std::min(oa.j2 + o6.j2, ob.j2 + o3.j2) / 2;

                    int twoJp_min = std::max(std::abs(oa.j2 - 2 * J1p), std::abs(ob.j2 - 2 * J2p));
                    int twoJp_max = std::min(oa.j2 + 2 * J1p, ob.j2 + 2 * J2p);

                    //          if ( (std::abs(o1.tz2+o2.tz2-oa.tz2-ob.tz2)!=2*X.GetTRank()) and (std::abs(o1.tz2+o2.tz2-oa.tz2-ob.tz2)!=2*Y.GetTRank()) ) continue;

                    for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                    {
                      double x12a45b = X3.GetME_pn(J1p, J2p, twoJp, I1, I2, a, I4, I5, b);
                      double y12a45b = Y3.GetME_pn(J1p, J2p, twoJp, I1, I2, a, I4, I5, b);

                      for (int Ja6 = Ja6_min; Ja6 <= Ja6_max; Ja6++)
                      {
                        double xb3a6 = X2.GetTBME_J(Ja6, b, I3, a, I6);
                        double yb3a6 = Y2.GetTBME_J(Ja6, b, I3, a, I6);

                        double ninej = Z.modelspace->GetNineJ(ob.j2 * 0.5, twoJp * 0.5, J2p, o3.j2 * 0.5, J1p, twoJ * 0.5, Ja6, oa.j2 * 0.5, o6.j2 * 0.5);
                        //                          double ninej = AngMom::NineJ( ob.j2*0.5, twoJp*0.5, J2p,  o3.j2*0.5, J1p, twoJ*0.5, Ja6, oa.j2*0.5, o6.j2*0.5);
                        int phase = AngMom::phase((o3.j2 + o6.j2) / 2 + J1p + J2p + twoJ);

                        zijklmn -= (oa.occ - ob.occ) * Pijk * Plmn * (2 * Ja6 + 1) * (twoJp + 1) * phase * ninej * (xb3a6 * y12a45b - yb3a6 * x12a45b);

                      } // Ja6
                    } // twoJp
                  } // b
                } // a
              } // J2p
            } // perm_lmn
          } // J1p
        } // perm_ijk
        Z3.AddToME_pn_ch(ch3bra, ch3ket, ibra, iket, zijklmn); // this needs to be modified for beta decay
      } // for iket
        //    }//ibra
    } // ch
    std::cout << "Ref " << __func__ << " Done" << std::endl;
  } // comm233_phss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1J2J3_ijklmn =  1/6 sum_abc sum_J4 (na nb nc + (1-na)(1-nb)(1-nc)) ( XJ1J4J3_ijkabc YJ4J2J3_abclmn - YJ1J4J3_ijkabc XJ4J2J3_abclmn )
  ///
  void comm333_ppp_hhhss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;

    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, static_cast<size_t>(ibra)}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();

    //  std::vector< std::array<size_t,2> > bra_ket_channels;
    //  for ( auto& it : Z.ThreeBody.Get_ch_start() )
    //  {
    //     bra_ket_channels.push_back( { it.first[0],it.first[1] } ); // (ch_bra, ch_ket)
    //  }
    //
    //
    //  size_t n_bra_ket_ch = bra_ket_channels.size();
    //  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
    //  for (size_t ch3=0; ch3<nch3; ch3++)
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ibra_ket = 0; ibra_ket < n_bra_ket_ch; ibra_ket++)
    {
      size_t ch3bra = bra_ket_channels[ibra_ket][0];
      size_t ch3ket = bra_ket_channels[ibra_ket][1];
      size_t ibra = bra_ket_channels[ibra_ket][2];
      auto &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch3bra);
      auto &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch3ket);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      size_t nkets3 = Tbc_ket.GetNumberKets();

      //    for (size_t ibra=0; ibra<nbras3; ibra++)
      //    {

      size_t iket_max = nkets3;
      if (ch3bra == ch3ket)
        iket_max = ibra;
      for (size_t iket = 0; iket <= iket_max; iket++)
      {

        double zijklmn = 0;

        for (size_t iket_abc = 0; iket_abc < nkets3; iket_abc++)
        {
          Ket3 &ket_abc = Tbc_ket.GetKet(iket_abc);
          double counting_factor = 6.0; // in general, 6 different permutations of abc
          if (ket_abc.p == ket_abc.r)
            counting_factor = 1.0; // if a=c, we must have a=b=c because we store a<=b<=c, and there is only one permutation
          else if ((ket_abc.p == ket_abc.q) or (ket_abc.q == ket_abc.r))
            counting_factor = 3; // if two orbits match, only 3 different permutations

          double nanbnc = ket_abc.op->occ * ket_abc.oq->occ * ket_abc.oR->occ;
          double nanbnc_bar = (1 - ket_abc.op->occ) * (1 - ket_abc.oq->occ) * (1 - ket_abc.oR->occ);

          double xijkabc = X3.GetME_pn_ch(ch3bra, ch3ket, ibra, iket_abc);
          double yijkabc = Y3.GetME_pn_ch(ch3bra, ch3ket, ibra, iket_abc);
          double xabclmn = X3.GetME_pn_ch(ch3ket, ch3ket, iket_abc, iket);
          double yabclmn = Y3.GetME_pn_ch(ch3ket, ch3ket, iket_abc, iket);
          zijklmn += 1. / 6 * counting_factor * (nanbnc + nanbnc_bar) * (xijkabc * yabclmn - yijkabc * xabclmn);
        } // for abc

        Z3.AddToME_pn_ch(ch3bra, ch3ket, ibra, iket, zijklmn); // this needs to be modified for beta decay

      } // iket : lmn
      //    }//ibra : ijk
    } // chbra, chket

  } // comm333_ppp_hhhss

  void comm333_pph_hhpss(const Operator &X, const Operator &Y, Operator &Z)
  {

    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z3 = Z.ThreeBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    size_t nch2 = Z.modelspace->GetNumberTwoBodyChannels();
    size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();

#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (size_t ch3 = 0; ch3 < nch3; ch3++)
    {
      auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      int twoJ = Tbc.twoJ;
      double Jtot = 0.5 * twoJ;
      for (size_t ibra = 0; ibra < nkets3; ibra++)
      {
        auto &bra = Tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        Orbit &ok = Z.modelspace->GetOrbit(k);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        double jk = 0.5 * ok.j2;
        int J1 = bra.Jpq;

        std::vector<std::array<size_t, 3>> ijk = {{i, j, k}, {k, j, i}, {i, k, j}};
        std::vector<int> J1p_min = {J1, std::abs(ok.j2 - oj.j2) / 2, std::abs(oi.j2 - ok.j2) / 2};
        std::vector<int> J1p_max = {J1, (ok.j2 + oj.j2) / 2, (oi.j2 + ok.j2) / 2};
        std::vector<std::vector<double>> recouple_ijk = {{1}, {}, {}};
        for (int J1p = J1p_min[1]; J1p <= J1p_max[1]; J1p++)
          recouple_ijk[1].push_back(sqrt((2 * J1 + 1) * (2 * J1p + 1)) * Z.modelspace->GetSixJ(ji, jj, J1, jk, Jtot, J1p));

        for (int J1p = J1p_min[2]; J1p <= J1p_max[2]; J1p++)
          recouple_ijk[2].push_back(-Z.modelspace->phase((oj.j2 + ok.j2) / 2 + J1 + J1p) * sqrt((2 * J1 + 1) * (2 * J1p + 1)) * Z.modelspace->GetSixJ(jj, ji, J1, jk, Jtot, J1p));

        //      for (size_t iket=ibra; iket<nkets3; iket++)
        for (size_t iket = 0; iket <= ibra; iket++)
        {
          auto &ket = Tbc.GetKet(iket);
          size_t l = ket.p;
          size_t m = ket.q;
          size_t n = ket.r;
          Orbit &ol = Z.modelspace->GetOrbit(l);
          Orbit &om = Z.modelspace->GetOrbit(m);
          Orbit &on = Z.modelspace->GetOrbit(n);
          double jl = 0.5 * ol.j2;
          double jm = 0.5 * om.j2;
          double jn = 0.5 * on.j2;
          int J2 = ket.Jpq;

          std::vector<std::array<size_t, 3>> lmn = {{l, m, n}, {n, m, l}, {l, n, m}};
          std::vector<int> J2p_min = {J2, std::abs(on.j2 - om.j2) / 2, std::abs(ol.j2 - on.j2) / 2};
          std::vector<int> J2p_max = {J2, (on.j2 + om.j2) / 2, (ol.j2 + on.j2) / 2};
          std::vector<std::vector<double>> recouple_lmn = {{1}, {}, {}};

          for (int J2p = J2p_min[1]; J2p <= J2p_max[1]; J2p++)
            recouple_lmn[1].push_back(sqrt((2 * J2 + 1) * (2 * J2p + 1)) * Z.modelspace->GetSixJ(jl, jm, J2, jn, Jtot, J2p));

          for (int J2p = J2p_min[2]; J2p <= J2p_max[2]; J2p++)
            recouple_lmn[2].push_back(-Z.modelspace->phase((om.j2 + on.j2) / 2 + J2 + J2p) * sqrt((2 * J2 + 1) * (2 * J2p + 1)) * Z.modelspace->GetSixJ(jm, jl, J2, jn, Jtot, J2p));

          double z_ijklmn = 0;

          for (int perm_ijk = 0; perm_ijk < 3; perm_ijk++)
          {
            size_t I1 = ijk[perm_ijk][0];
            size_t I2 = ijk[perm_ijk][1];
            size_t I3 = ijk[perm_ijk][2];
            Orbit &o1 = Z.modelspace->GetOrbit(I1);
            Orbit &o2 = Z.modelspace->GetOrbit(I2);
            Orbit &o3 = Z.modelspace->GetOrbit(I3);

            double j3 = 0.5 * o3.j2;
            for (int J1p = J1p_min[perm_ijk]; J1p <= J1p_max[perm_ijk]; J1p++)
            {
              double rec_ijk = recouple_ijk[perm_ijk].at(J1p - J1p_min[perm_ijk]);
              for (int perm_lmn = 0; perm_lmn < 3; perm_lmn++)
              {
                size_t I4 = lmn[perm_lmn][0];
                size_t I5 = lmn[perm_lmn][1];
                size_t I6 = lmn[perm_lmn][2];
                Orbit &o4 = Z.modelspace->GetOrbit(I4);
                Orbit &o5 = Z.modelspace->GetOrbit(I5);
                Orbit &o6 = Z.modelspace->GetOrbit(I6);
                double j6 = 0.5 * o6.j2;
                for (int J2p = J2p_min[perm_lmn]; J2p <= J2p_max[perm_lmn]; J2p++)
                {
                  double rec_lmn = recouple_lmn[perm_lmn].at(J2p - J2p_min[perm_lmn]);

                  for (size_t ch2 = 0; ch2 < nch2; ch2++)
                  {
                    auto &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch2);
                    if (std::abs(Tbc.twoTz - 2 * tbc_ab.Tz) == 5)
                      continue; // TODO there are probably other checks at the channel level...
                    size_t nkets_ab = tbc_ab.GetNumberKets();
                    int Jab = tbc_ab.J;

                    for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
                    {
                      Ket &ket_ab = tbc_ab.GetKet(iket_ab);
                      size_t a = ket_ab.p;
                      size_t b = ket_ab.q;
                      Orbit &oa = Z.modelspace->GetOrbit(a);
                      Orbit &ob = Z.modelspace->GetOrbit(b);
                      if (std::abs(oa.occ * ob.occ) < 1e-6 and std::abs((1 - oa.occ) * (1 - ob.occ)) < 1e-6)
                        continue;

                      for (auto c : Z.modelspace->all_orbits)
                      {
                        Orbit &oc = Z.modelspace->GetOrbit(c);
                        double occ_factor = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
                        if (std::abs(occ_factor) < 1e-6)
                          continue;
                        if (a == b)
                          occ_factor *= 0.5; // because we only sum b<a
                        double jc = 0.5 * oc.j2;

                        if (((oa.l + ob.l + o3.l + o4.l + o5.l + oc.l) % 2 == 0) and ((o1.l + o2.l + oc.l + oa.l + ob.l + o6.l) % 2 == 0) and ((oa.tz2 + ob.tz2 + o3.tz2) == (o4.tz2 + o5.tz2 + oc.tz2)) and ((o1.tz2 + o2.tz2 + oc.tz2) == (oa.tz2 + ob.tz2 + o6.tz2)))
                        {

                          int twoJx_min = std::max(std::abs(2 * Jab - o3.j2), std::abs(2 * J2p - oc.j2));
                          int twoJx_max = std::min(2 * Jab + o3.j2, 2 * J2p + oc.j2);
                          int twoJy_min = std::max(std::abs(2 * J1p - oc.j2), std::abs(2 * Jab - o6.j2));
                          int twoJy_max = std::min(2 * J1p + oc.j2, 2 * Jab + o6.j2);
                          if (twoJx_min > twoJx_max or twoJy_min > twoJy_max)
                            continue;
                          //                        if (twoJx_min <= twoJx_max and twoJy_min<=twoJy_max)
                          //                        {
                          std::vector<double> xab345c((twoJx_max - twoJx_min) / 2 + 1, 0);
                          std::vector<double> yab345c((twoJx_max - twoJx_min) / 2 + 1, 0);
                          std::vector<double> x12cab6((twoJy_max - twoJy_min) / 2 + 1, 0);
                          std::vector<double> y12cab6((twoJy_max - twoJy_min) / 2 + 1, 0);
                          for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                          {
                            size_t iJx = (twoJx - twoJx_min) / 2;
                            xab345c[iJx] = X3.GetME_pn(Jab, J2p, twoJx, a, b, I3, I4, I5, c);
                            yab345c[iJx] = Y3.GetME_pn(Jab, J2p, twoJx, a, b, I3, I4, I5, c);
                          }
                          for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                          {
                            size_t iJy = (twoJy - twoJy_min) / 2;
                            x12cab6[iJy] = X3.GetME_pn(J1p, Jab, twoJy, I1, I2, c, a, b, I6);
                            y12cab6[iJy] = Y3.GetME_pn(J1p, Jab, twoJy, I1, I2, c, a, b, I6);
                          }
                          for (int twoJx = twoJx_min; twoJx <= twoJx_max; twoJx += 2)
                          {
                            double JJx = 0.5 * twoJx;
                            size_t iJx = (twoJx - twoJx_min) / 2;
                            for (int twoJy = twoJy_min; twoJy <= twoJy_max; twoJy += 2)
                            {
                              double JJy = 0.5 * twoJy;
                              size_t iJy = (twoJy - twoJy_min) / 2;
                              double hats = (twoJx + 1) * (twoJy + 1);
                              double ninej = Z.modelspace->GetNineJ(j3, Jab, JJx, J1p, JJy, jc, Jtot, j6, J2p);
                              z_ijklmn += rec_ijk * rec_lmn * occ_factor * hats * ninej * (xab345c[iJx] * y12cab6[iJy] - yab345c[iJx] * x12cab6[iJy]);
                            } // for twoJy
                          } // for twoJx
                            //                        }
                        } // Z1 block

                      } // for c
                    } // for iket_ab
                  } // for ch2

                } // for J2p
              } // for perm_lmn
            } // for J1p
          } // for perm_ijk

          Z3.AddToME_pn_ch(ch3, ch3, ibra, iket, z_ijklmn);
        } // for iket
      } // for ibra
    } // for ch3

  } // comm333_pph_hhpss

  /// Start of scalar-tensor commutators

  // This has not yet been validated, and is almost certainly wrong.
  void comm222_phst(const Operator &X, const Operator &Y, Operator &Z)
  {
    int lambda = Y.GetJRank();
    for (auto &iter : Z.TwoBody.MatEl)
    {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      auto &Z2MAT = iter.second;
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      for (size_t ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t q = bra.q;
        Orbit &op = Z.modelspace->GetOrbit(p);
        Orbit &oq = Z.modelspace->GetOrbit(q);
        for (size_t iket = 0; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          Orbit &oR = Z.modelspace->GetOrbit(r);
          Orbit &os = Z.modelspace->GetOrbit(s);
          double zpqrs = 0;

          std::vector<std::vector<size_t>> perms = {{p, q, r, s}, {q, p, r, s}, {p, q, s, r}, {q, p, s, r}};
          std::vector<int> phases = {1,
                                     -AngMom::phase((op.j2 + oq.j2) / 2 - J1),
                                     -AngMom::phase((oR.j2 + os.j2) / 2 - J2),
                                     AngMom::phase((op.j2 + oq.j2) / 2 - J1) * AngMom::phase((oR.j2 + os.j2) / 2 - J2)};

          for (int iperm = 0; iperm < 4; iperm++)
          {

            int I1 = perms[iperm][0];
            int I2 = perms[iperm][1];
            int I3 = perms[iperm][2];
            int I4 = perms[iperm][3];
            int phaseperm = phases[iperm];
            Orbit &o1 = Z.modelspace->GetOrbit(I1);
            Orbit &o2 = Z.modelspace->GetOrbit(I2);
            Orbit &o3 = Z.modelspace->GetOrbit(I3);
            Orbit &o4 = Z.modelspace->GetOrbit(I4);

            int J3min = std::abs(o1.j2 - o4.j2) / 2;
            int J3max = (o1.j2 + o4.j2) / 2;
            int J4min = std::abs(o3.j2 - o2.j2) / 2;
            int J4max = (o3.j2 + o2.j2) / 2;
            for (int J3 = J3min; J3 <= J3max; J3++)
            {
              for (int J4 = J4min; J4 <= J4max; J4++)
              {
                if (std::abs(J3 - J4) > lambda)
                  continue;
                if ((J3 + J4) < lambda)
                  continue;
                double Zbar_1432 = 0;

                for (size_t a : Z.modelspace->all_orbits)
                {
                  Orbit &oa = Z.modelspace->GetOrbit(a);
                  for (size_t b : Z.modelspace->all_orbits)
                  {
                    Orbit &ob = Z.modelspace->GetOrbit(b);
                    double nanb = oa.occ - ob.occ;
                    if (std::abs(nanb) < 1e-9)
                      continue;

                    double Xbar_psab = 0;
                    double Ybar_abrq = 0;

                    int J5min = std::max(std::abs(o1.j2 - o4.j2), std::abs(oa.j2 - ob.j2)) / 2.;
                    int J5max = std::min(o1.j2 + o4.j2, oa.j2 + ob.j2) / 2.;
                    for (int J5 = J5min; J5 <= J5max; J5++)
                    {
                      double sixj = AngMom::SixJ(o1.j2 / 2., o4.j2 / 2., J3, oa.j2 / 2., ob.j2 / 2., J5);
                      Xbar_psab -= (2 * J5 + 1) * sixj * X.TwoBody.GetTBME_J(J5, J5, I1, b, a, I4);
                    }
                    int J6min = std::abs(oa.j2 - o2.j2) / 2.;
                    int J6max = (oa.j2 + o2.j2) / 2.;
                    int J7min = std::abs(ob.j2 - o3.j2) / 2.;
                    int J7max = (ob.j2 + o3.j2) / 2.;
                    for (int J6 = J6min; J6 <= J6max; J6++)
                    {
                      for (int J7 = J7min; J7 <= J7max; J7++)
                      {
                        if (std::abs(J6 - J7) > lambda or (J6 + J7) < lambda)
                          continue;
                        double ninejY = AngMom::NineJ(oa.j2 / 2., o2.j2 / 2., J6, ob.j2 / 2., o3.j2 / 2., J7, J3, J4, lambda);
                        double hats = sqrt((2 * J3 + 1) * (2 * J4 + 1) * (2 * J6 + 1) * (2 * J7 + 1));
                        int phaseY = AngMom::phase((ob.j2 + o2.j2) / 2 + J4 + J7);
                        Ybar_abrq -= hats * phaseY * ninejY * Y.TwoBody.GetTBME_J(J6, J7, a, I2, I3, b);
                        if (ch_bra == 1)
                        {
                          std::cout << "  ========= J6,J7= " << J6 << " " << J7 << "  " << hats << " * " << phaseY << " * " << ninejY << " * " << Y.TwoBody.GetTBME_J(J6, J7, a, I2, I3, b) << " I got Ya23b from Y2.GetTBME_J(" << J6 << " " << J7 << " " << a << " " << I2 << " " << I3 << " " << b << " )" << std::endl;
                        }
                      }
                    }
                    Zbar_1432 += nanb * Xbar_psab * Ybar_abrq;

                    if (ch_bra == 1)
                    {
                      std::cout << "     a b " << a << " " << b << " " << nanb << " *  ( " << Xbar_psab << " * " << Ybar_abrq << " ) -> zbar_1432 = " << Zbar_1432 << std::endl;
                    }

                  } // for b
                } // for a

                double ninej = AngMom::NineJ(o1.j2 / 2., o4.j2 / 2., J3, o2.j2 / 2., o3.j2 / 2., J4, J1, J2, lambda);
                //                         zpqrs += phaseperm* sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) ) * nanb * ninej * AngMom::phase((o2.j2+o4.j2)/2 +J2+J4) * Xbar_psab * Ybar_abrq;
                zpqrs += phaseperm * sqrt((2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1)) * ninej * AngMom::phase((o2.j2 + o4.j2) / 2 + J2 + J4) * Zbar_1432;

                if (ch_bra == 1)
                {
                  std::cout << __func__ << " " << __LINE__ << " channel J = " << tbc_bra.J << "  J3,J4 " << J3 << " " << J4 << " iperm " << iperm << " : "
                            << phaseperm << " * " << sqrt((2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1)) << " * " << ninej << " * " << AngMom::phase((o2.j2 + o4.j2) / 2 + J2 + J4) << " * " << Zbar_1432 << " (<-Zbar_1432, zpqrs ->)  " << zpqrs << std::endl;
                }
              } // for J4
            } // for J3
          } // for iperm

          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;

          Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zpqrs);
        } // for iket
      } // for ibra
        //    if ( ch_bra==3) std::cout << __func__ << "   ch=3 Jpt " << tbc_bra.J << " " << tbc_bra.parity << " " << tbc_bra.Tz << "   " << tbc_bra.GetKet(0).p << " " << tbc_bra.GetKet(0).q << " Z2 = " << std::endl << iter.second << std::endl;
    } // for iter
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// scalar-tensor commutators with 3-body

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij^\lamda = 1/12 sum_abcde sum_J1J2J (na nb n¯c n¯d n¯e + ¯na n¯b nc nd ne) (2J+1)/(2ji+1)
  ///                             * (-)^(j + lamda + J1 + j0)
  ///                             * (X^J1,j0,J2,j0;0_abicde Y^J2,j0,J1,j1;lamda_cdeabj
  ///                             -  Y^J1,j0,J2,j1;lamda_abicde  X^J2,j1,J1,j0;0_cdeabj )
  ///
  ///
  void comm331st(const Operator &X, const Operator &Y, Operator &Z)
  {
    Z.modelspace->PreCalculateSixJ();
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z1 = Z.OneBody;

    double tstart = omp_get_wtime();
    int Lambda = Z.GetJRank();
    int hZ = Z.IsHermitian() ? +1 : -1;

    size_t norb = Z.modelspace->GetNumberOrbits();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < norb; i++)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        if (j < i)
          continue;

        Orbit &oj = Z.modelspace->GetOrbit(j);
        if (std::abs(oj.j2 - oi.j2) > Lambda * 2 or (oi.j2 + oj.j2) < Lambda * 2)
          continue;

        double zij = 0;
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

                for (size_t e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);

                  // Tensor operator may change parity. So orbit i and j may not have the same parity
                  if (
                      ((oa.l + ob.l + oi.l + oc.l + od.l + oe.l + X.GetParity()) % 2 != 0) and ((oa.l + ob.l + oj.l + oc.l + od.l + oe.l + Y.GetParity()) % 2 != 0) and
                      ((oa.l + ob.l + oi.l + oc.l + od.l + oe.l + Y.GetParity()) % 2 != 0) and ((oa.l + ob.l + oj.l + oc.l + od.l + oe.l + X.GetParity()) % 2 != 0))
                    continue;

                  if (
                      (std::abs(oa.tz2 + ob.tz2 + oi.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * X.GetTRank()) and (std::abs(oa.tz2 + ob.tz2 + oj.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * Y.GetTRank()) and
                      (std::abs(oa.tz2 + ob.tz2 + oi.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * Y.GetTRank()) and (std::abs(oa.tz2 + ob.tz2 + oj.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * X.GetTRank()))
                    continue;

                  double occupation_factor = oa.occ * ob.occ * (1 - oc.occ) * (1 - od.occ) * (1 - oe.occ) // fixed mistake found by Matthias Heinz Oct 2022
                                             + (1 - oa.occ) * (1 - ob.occ) * oc.occ * od.occ * oe.occ;
                  // occupation_factor = 1.;
                  if (std::abs(occupation_factor) < 1e-7)
                    continue;

                  int J1min = std::abs(oa.j2 - ob.j2) / 2;
                  int J1max = (oa.j2 + ob.j2) / 2;

                  int J2min = std::abs(oc.j2 - od.j2) / 2;
                  int J2max = (oc.j2 + od.j2) / 2;

                  for (int J1 = J1min; J1 <= J1max; J1++)
                  {
                    if (a == b and J1 % 2 > 0)
                      continue;

                    for (int J2 = J2min; J2 <= J2max; J2++)
                    {
                      if (c == d and J2 % 2 > 0)
                        continue;

                      int j0min = std::max(std::abs(2 * J1 - oi.j2), std::abs(2 * J2 - oe.j2));
                      int j0max = std::min(2 * J1 + oi.j2, 2 * J2 + oe.j2);

                      int j1min = std::abs(2 * J1 - oj.j2);
                      int j1max = 2 * J1 + oj.j2;

                      for (int j0 = j0min; j0 <= j0max; j0 += 2)
                      {
                        for (int j1 = j1min; j1 <= j1max; j1 += 2)
                        {
                          if (std::abs(j0 - j1) > Lambda * 2 or (j0 + j1) < Lambda * 2)
                            continue;

                          double sixj1 = AngMom::SixJ(oj.j2 / 2., oi.j2 / 2., Lambda, j0 / 2., j1 / 2., J1);
                          if (std::abs(sixj1) < 1.e-6)
                            continue;
                          double xabicde = X3.GetME_pn(J1, J2, j0, a, b, i, c, d, e);     // scalar
                          double ycdeabj = Y3.GetME_pn(J2, j0, J1, j1, c, d, e, a, b, j); // tensor

                          int phase = AngMom::phase((oj.j2 + j0) / 2 + J1 + Lambda);
                          zij += 1. / 12 * phase * occupation_factor * sqrt(j0 + 1) * sqrt(j1 + 1) * sixj1 * (xabicde * ycdeabj); //
                        } // j1
                      } // j0

                      j1min = std::max(std::abs(2 * J1 - oj.j2), std::abs(2 * J2 - oe.j2));
                      j1max = std::min(2 * J1 + oj.j2, 2 * J2 + oe.j2);

                      j0min = std::abs(2 * J1 - oi.j2);
                      j0max = 2 * J1 + oi.j2;

                      for (int j0 = j0min; j0 <= j0max; j0 += 2)
                      {
                        for (int j1 = j1min; j1 <= j1max; j1 += 2)
                        {
                          if (std::abs(j0 - j1) > Lambda * 2 or (j0 + j1) < Lambda * 2)
                            continue;

                          double sixj1 = AngMom::SixJ(oj.j2 / 2., oi.j2 / 2., Lambda, j0 / 2., j1 / 2., J1);
                          if (std::abs(sixj1) < 1.e-6)
                            continue;

                          double yabicde = Y3.GetME_pn(J1, j0, J2, j1, a, b, i, c, d, e); // tensor
                          double xcdeabj = X3.GetME_pn(J2, J1, j1, c, d, e, a, b, j);     // sclar

                          int phase = AngMom::phase((oj.j2 + j0) / 2 + J1 + Lambda);
                          zij -= 1. / 12 * phase * occupation_factor * sqrt(j0 + 1) * sqrt(j1 + 1) * sixj1 * (yabicde * xcdeabj); //
                        } // j1
                      } // j0

                    } // J2
                  } // J1
                } // e
              } // d
            } // c
          } // b
        } // a
        Z1(i, j) += zij;
        if ( i != j)
          Z1(j, i) += AngMom::phase((oi.j2 - oj.j2)/2 ) * hZ * zij;
      } // j
    } // i

    X.profiler.timer[__func__] += omp_get_wtime() - tstart;
  } // comm331st

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:   Z^(J1j1,J2j2)Lamda_ijklmn = sum_{a, J3} PJ1j1(ij/k) PJ1j2(lm/n) sqrt( (2J1+1) (2j1+1)) (2j2+1) (2J3+1)) 
  ///                                           (  (-)^j2+jn+Lamda+ J3 + 1
  ///                                           { jk ja J3 } { J3 J2 Lamda }  XJ1_ijna YJ3J2_Lamda_kalm
  ///                                           { jn j1 J1 } { j2 j1 jn    }
  ///                                            - (-)^j2+jk+Lamda+ J1
  ///                                           { jn ja J3 } { J3 J1 Lamda }  YJ1J3_Lamda_ijna XJ2_kalm
  ///                                           { jk j2 J2 } { j1 j2 jk    }
  ///
  void comm223st(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &Z3 = Z.ThreeBody;
    int Lambda = Z.GetJRank();
    Z.modelspace->PreCalculateSixJ();
    // Permutations of indices which are needed to produce antisymmetrized matrix elements  P(ij/k) |ijk> = |ijk> - |kji> - |ikj>
    const std::array<ThreeBodyStorage::Permutation, 3> index_perms = {ThreeBodyStorage::ABC, ThreeBodyStorage::CBA, ThreeBodyStorage::ACB};

    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, static_cast<size_t>(ibra)}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ibra_ket = 0; ibra_ket < n_bra_ket_ch; ibra_ket++)
    {
      size_t ch3bra = bra_ket_channels[ibra_ket][0];
      size_t ch3ket = bra_ket_channels[ibra_ket][1];
      size_t ibra = bra_ket_channels[ibra_ket][2];
      auto &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch3bra);
      auto &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch3ket);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      size_t nkets3 = Tbc_ket.GetNumberKets();


      // int twoJ = Tbc_bra.twoJ; // Scalar commutator so J is the same in bra and ket channel
      int twoj1 = Tbc_bra.twoJ; 
      int twoj2 = Tbc_ket.twoJ; 

      auto &bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      Orbit &ok = Z.modelspace->GetOrbit(k);
      double ji = 0.5 * oi.j2;
      double jj = 0.5 * oj.j2;
      double jk = 0.5 * ok.j2;

      int J1 = bra.Jpq;

      size_t iket_max = nkets3 - 1;
      if (ch3bra == ch3ket)
        iket_max = ibra;
      for (size_t iket = 0; iket <= iket_max; iket++)
      {
        auto &ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit &ol = Z.modelspace->GetOrbit(l);
        Orbit &om = Z.modelspace->GetOrbit(m);
        Orbit &on = Z.modelspace->GetOrbit(n);
        double jl = 0.5 * ol.j2;
        double jm = 0.5 * om.j2;
        double jn = 0.5 * on.j2;
        int J2 = ket.Jpq;

        double zijklmn = 0;
        /// BEGIN THE SLOW BIT...

        // Now we need to loop over the permutations in ijk and then lmn
        for (auto perm_ijk : index_perms)
        {
          size_t I1, I2, I3;
          Z3.Permute(perm_ijk, i, j, k, I1, I2, I3);
          Orbit &o1 = Z.modelspace->GetOrbit(I1);
          Orbit &o2 = Z.modelspace->GetOrbit(I2);
          Orbit &o3 = Z.modelspace->GetOrbit(I3);

          int J1p_min = J1;
          int J1p_max = J1;
          if (perm_ijk != ThreeBodyStorage::ABC)
          {
            J1p_min = std::max(std::abs(o1.j2 - o2.j2), std::abs(twoj1 - o3.j2)) / 2;
            J1p_max = std::min(o1.j2 + o2.j2, twoj1 + o3.j2) / 2;
          }

          int parity_12 = (o1.l + o2.l) % 2;
          int Tz_12 = (o1.tz2 + o2.tz2) / 2;

          double j3 = 0.5 * o3.j2;

          for (int J1p = J1p_min; J1p <= J1p_max; J1p++)
          {

            double rec_ijk = Z3.RecouplingCoefficient(perm_ijk, ji, jj, jk, J1p, J1, twoj1);
            rec_ijk *= Z3.PermutationPhase(perm_ijk); // do we get a fermionic minus sign?

            for (auto perm_lmn : index_perms)
            {
              size_t I4, I5, I6;
              Z3.Permute(perm_lmn, l, m, n, I4, I5, I6);
              Orbit &o4 = Z.modelspace->GetOrbit(I4);
              Orbit &o5 = Z.modelspace->GetOrbit(I5);
              Orbit &o6 = Z.modelspace->GetOrbit(I6);

              double j6 = 0.5 * o6.j2;

              int J2p_min = J2;
              int J2p_max = J2;
              if (perm_lmn != ThreeBodyStorage::ABC)
              {
                J2p_min = std::max(std::abs(o4.j2 - o5.j2), std::abs(twoj2 - o6.j2)) / 2;
                J2p_max = std::min(o4.j2 + o5.j2, twoj2 + o6.j2) / 2;
              }

              for (size_t a : Z.modelspace->all_orbits)
              {
                Orbit &oa = Z.modelspace->GetOrbit(a);
                int j2a = oa.j2;

                int dTz_126a = (o1.tz2 + o2.tz2 - o6.tz2 - oa.tz2) / 2;
                int dTz_3a45 = (o3.tz2 + oa.tz2 - o4.tz2 - o5.tz2) / 2;
                if (std::abs(dTz_126a) != X.GetTRank() and std::abs(dTz_126a) != Y.GetTRank())
                  continue;
                if (std::abs(dTz_3a45) != X.GetTRank() and std::abs(dTz_3a45) != Y.GetTRank())
                  continue;

                for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                {
                  double rec_lmn = Z3.RecouplingCoefficient(perm_lmn, jl, jm, jn, J2p, J2, twoj2);
                  rec_lmn *= Z3.PermutationPhase(perm_lmn); // do we get a fermionic minus sign?

                  // direct term
                  int J3_min = std::abs(o3.j2 - j2a) / 2;
                  int J3_max = ( o3.j2 + j2a ) / 2;
                  for (int J3 = J3_min ; J3 <= J3_max; J3++)
                  {
                    if ( J3 +  J2p < Lambda or  std::abs(J3 - J2p) > Lambda )
                      continue;
                    double sixj_1, sixj_2;
                    sixj_1 = AngMom::SixJ(o3.j2 / 2., j2a / 2.,      J3,
                                          o6.j2 / 2., twoj1 / 2.,    J1p);
                    sixj_2 = AngMom::SixJ(J3,         J2p,           Lambda,
                                          twoj2 / 2., twoj1 / 2.,    o6.j2 / 2.);

                    int phase = AngMom::phase((o6.j2 + twoj2) / 2 + J3 + Lambda);
                    double facotrs = sqrt((2 * J1p + 1) * (2 * J3 + 1) * (twoj1 + 1) * (twoj2 + 1)); 
                    double x_126a = X2.GetTBME_J(J1p, J1p, I1, I2, I6, a);
                    double y_3a45 = Y2.GetTBME_J(J3, J2p, I3, a, I4, I5);
                    zijklmn += rec_ijk * rec_lmn * sixj_1 * sixj_2 * phase * facotrs * (x_126a * y_3a45); 
                  } // J3

                  // exchange term
                  J3_min = std::abs(o6.j2 - j2a) / 2;
                  J3_max = ( o6.j2 + j2a ) / 2;
                  for (int J3 = J3_min ; J3 <= J3_max; J3++)
                  {
                    if ( J3 +  J1p < Lambda or  std::abs(J3 - J1p) > Lambda )
                      continue;
                    double sixj_1, sixj_2;
                    sixj_1 = AngMom::SixJ(o6.j2 / 2., j2a / 2.,      J3,
                                          o3.j2 / 2., twoj2 / 2.,    J2p);
                    sixj_2 = AngMom::SixJ(J3,         J1p,           Lambda,
                                          twoj1 / 2., twoj2 / 2.,    o3.j2 / 2.);

                    int phase = AngMom::phase((o3.j2 + twoj2) / 2 + J1p + Lambda);
                    double facotrs = sqrt((2 * J2p + 1) * (2 * J3 + 1) * (twoj1 + 1) * (twoj2 + 1)); 

                    double y_126a = Y2.GetTBME_J(J1p, J3, I1, I2, I6, a);
                    double x_3a45 = X2.GetTBME_J(J2p, J2p, I3, a, I4, I5);
                    zijklmn -= rec_ijk * rec_lmn * sixj_1 * sixj_2 * phase * facotrs * ( y_126a * x_3a45 );
                  } // J3
                } // for J2p
              } // for orbit a 
            } // for perm_lmn
          } // for J1p
        } // for perm_ijk

        Z3.AddToME_pn_ch(ch3bra, ch3ket, ibra, iket, zijklmn); // this needs to be modified for beta decay
      } // for iket
    } // for ch3

  } // comm233st

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Zij^Lamda = 1/4 sum_abcd sum_J1j1j2 na nb (1-nc) (1-nd) sqrt((2j1+1)(2j2+1)) 
  ///                            (-)^(j_j + Lamda + J1 + j1) { jj  ji  Lamda }     
  ///                                                        { j1  j2  J1    } 
  ///                            ( X^J1_abcd Y^(J1j1, J1j2)Lamda_cdiabj - Y^(J1j1, J1j2)Lamda_abicdj XJ1_cdab )
  ///
  ///                            + 1/4 sum_abcd sum_J1J2j1 na nb (1-nc) (1-nd) (2j1+1) 
  ///                            [(-)^(j_i + J1 + j1) { Lamda  J2  J1  }  X(J1j1,J2j1)0_abicdj Y(J1,J2)Lamda_cdab   
  ///                                                 { j1     ji  jj  } 
  ///                            -(-)^(j_i + J2 + j1) { Lamda  J2  J1  }  Y(J1,J2)Lamda_abcd X(J2j1,J1j1)0_cdiabj ]
  ///                                                 { j1     jj  ji  } 
  ///
  void comm231st(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z1 = Z.OneBody;
    int Lambda = Z.GetJRank();
    Z.modelspace->PreCalculateSixJ();
    size_t norb = Z.modelspace->GetNumberOrbits();

  #pragma omp parallel for schedule(dynamic,1)
    for (size_t i = 0; i < norb; i++)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double zij = 0;

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

                double occ_factor =  oa.occ * ob.occ * (1 - oc.occ) * (1 - od.occ);
                if (std::abs(occ_factor) < 1e-8)
                  continue;

                int J1min = std::max(std::abs(oa.j2 - ob.j2), std::abs(oc.j2 - od.j2)) / 2;
                int J1max = std::min(oa.j2 + ob.j2, oc.j2 + od.j2) / 2;
                for (int J1 = J1min; J1 <= J1max; J1++)
                {
                  if (a == b and J1 % 2 > 0)
                    continue;
                  if (c == d and J1 % 2 > 0)
                    continue;

                  double xabcd = X2.GetTBME_J(J1, J1, a, b, c, d);
                  double xcdab = X2.GetTBME_J(J1, J1, c, d, a, b);

                  int twoj1min = std::abs(2 * J1 - oi.j2);
                  int twoj1max = 2 * J1 + oi.j2;
                  int twoj2min = std::abs(2 * J1 - oj.j2);
                  int twoj2max = 2 * J1 + oj.j2;
                  for (int twoj1 = twoj1min; twoj1 <= twoj1max; twoj1 += 2)
                  {
                    for (int twoj2 = twoj2min; twoj2 <= twoj2max; twoj2 += 2)
                    {
                      double sixj = AngMom::SixJ(oi.j2 / 2., oj.j2 / 2.,  Lambda,
                                                 twoj2 / 2., twoj1 / 2.,  J1);
                      int phase = AngMom::phase((oj.j2 + twoj1) / 2 + J1 + Lambda);
                      double facotrs = sqrt((twoj1 + 1) * (twoj2 + 1)); 
                      double yabicdj = Y3.GetME_pn(J1, twoj1, J1, twoj2, a, b, i, c, d, j);
                      double ycdiabj = Y3.GetME_pn(J1, twoj1, J1, twoj2, c, d, i, a, b, j);
                      zij += 1. / 4 * occ_factor * phase * facotrs * sixj * (xabcd * ycdiabj- yabicdj * xcdab);
                    } // j2
                  } // j1
                } // J1

                J1min = std::abs(oa.j2 - ob.j2) / 2;
                J1max = (oa.j2 + ob.j2) / 2;
                int J2min = std::abs(oc.j2 - od.j2) / 2;
                int J2max = (oc.j2 + od.j2) / 2;
                for (int J1 = J1min; J1 <= J1max; J1++)
                {
                  if (a == b and J1 % 2 > 0)
                    continue;
                  for (int J2 = J2min; J2 <= J2max; J2++)
                  {
                    if (c == d and J2 % 2 > 0)
                      continue;
                    double yabcd = Y2.GetTBME_J(J1, J2, a, b, c, d);
                    double ycdab = Y2.GetTBME_J(J2, J1, c, d, a, b);
                    
                    int twoj1min = std::max(std::abs(oi.j2 - 2 * J1), std::abs(oj.j2 - 2 * J2)) ;
                    int twoj1max = std::min(oi.j2 + 2 * J1, oj.j2 + 2 * J2) ;
                    for (int twoj1 = twoj1min; twoj1 <= twoj1max; twoj1 += 2)
                    {
                      double sixj = AngMom::SixJ(Lambda,     J2,          J1,
                                                 twoj1 / 2., oi.j2 / 2.,  oj.j2 / 2.);
                      int phase = AngMom::phase((oi.j2 + twoj1) / 2 + J1);
                      double facotrs = (twoj1 + 1); 
                      double xabicdj = X3.GetME_pn(J1, J2, twoj1, a, b, i, c, d, j);
                      zij += 1. / 4 * occ_factor * phase * sixj * facotrs * (xabicdj * ycdab);
                    } // j1

                    twoj1min = std::max(std::abs(oi.j2 - 2 * J2), std::abs(oj.j2 - 2 * J1)) ;
                    twoj1max = std::min(oi.j2 + 2 * J2, oj.j2 + 2 * J1) ;
                    for (int twoj1 = twoj1min; twoj1 <= twoj1max; twoj1 += 2)
                    {
                      double sixj = AngMom::SixJ(Lambda,     J2,          J1,
                                                 twoj1 / 2., oj.j2 / 2.,  oi.j2 / 2.);
                      int phase = AngMom::phase((oi.j2 + twoj1) / 2 + J2);
                      double facotrs = (twoj1 + 1); 
                      double xcdiabj = X3.GetME_pn(J2, J1, twoj1, c, d, i, a, b, j);
                      zij -= 1. / 4 * occ_factor * phase * sixj * facotrs * (yabcd * xcdiabj);
                    } // j1

                  } // J2
                } // J1
              } // d
            } // c
          } // b
        } // a
        Z1(i, j) += zij;
      } // j
    } // i
  } // comm231st

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Z^(J1, J2)Lamda_ijkl = -1/2 sum_abc sum_{J3,j1,j2} ( nanb(1-nc) + (1-na)(1-nb)nc ) 
  ///                                           (1-PJ1_ij) sqrt{ (2J1+1)(2J3+1)(2j1+1)(2j2+1) } 
  ///                                            (-)^(J1 + Lamda + jc + j2)  { ji jj J1 } { Lamda  J1  J2 } 
  ///                                                                        { jc j1 J3 } { jc     j2  j1 }
  ///                                                       X^J2_cjab Y^(J3j1,J2j2)Lamda_abiklc
  ///                                       -1/2 sum_abc sum_{J3,J4,j1,j2} ( nanb(1-nc) + (1-na)(1-nb)nc ) 
  ///                                           (1-PJ1_ij) sqrt{ (2J1+1)(2J3+1) }    (2j1+1)(2j2+1)
  ///                                            (-)^(J1 + ji+ jc + j3) { Lamda J4 J3 } { jc j1 J2 } { Lamda  J1  J2 } 
  ///                                                                   { jc    jj j2 } { ji j2 J4 } { ji     j2  jj }
  ///                                                        Y^(J3,J4)Lamda_cjab X^(J4j1,J2j1)0_abiklc
  ///                                        1/2 sum_abc sum_{J3,j1,j2} ( nanb(1-nc) + (1-na)(1-nb)nc ) 
  ///                                           (1-PJ2_kl) sqrt{ (2J2+1)(2J3+1)(2j1+1)(2j2+1) } 
  ///                                            (-)^(J1 + Lamda + jc + j2)  { J2 Lamda J1 } { jk  jl  J2 } 
  ///                                                                        { j1 jc    j2 } { jc  j2  J3 }
  ///                                                        Y^(J1j1,J3j2)Lamda_ijcabk XJ3_abcl
  ///                                        1/2 sum_abc sum_{J3,J4j1,j2} ( nanb(1-nc) + (1-na)(1-nb)nc ) 
  ///                                           (1-PJ2_kl) sqrt{ (2J2+1)(2J4+1)}    (2j1+1)(2j2+1) 
  ///                                            (-)^(J1 + jk + jc + J3)  { jc j1 J1 } { Lamda  J3  J4 } { J2 J1 Lamda } 
  ///                                                                     { jk j2 j3 } { jc     jl  j2 } { j2 jl jk    }
  ///                                                        X^(J1j1,J3j1)_ijcabk Y^(J3,J4)Lamda_abcl
  ///
  void comm232st(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X2 = X.TwoBody;
    auto &X3 = X.ThreeBody;
    auto &Y2 = Y.TwoBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;
    int Lambda = Z.GetJRank();
    Z.modelspace->PreCalculateSixJ();
    int nch = Z.modelspace->GetNumberTwoBodyChannels();

    std::vector<std::array<size_t, 2>> channels;
    for (auto &iter : Z.TwoBody.MatEl)
      channels.push_back(iter.first);
    size_t nchans = channels.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < nchans; ich++)
    {
      size_t ch_bra = channels[ich][0];
      size_t ch_ket = channels[ich][1];
      auto &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      auto &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double ji = 0.5 * oi.j2;
        double jj = 0.5 * oj.j2;
        int ket_min = (ch_bra == ch_ket) ? ibra : 0;
        for (int iket = ket_min; iket < nkets; iket++)
        {
          double zijkl = 0;
          Ket &ket = tbc_ket.GetKet(iket);
          int k = ket.p;
          int l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          double jk = 0.5 * ok.j2;
          double jl = 0.5 * ol.j2;
          for (auto c : Z.modelspace->all_orbits)
          {
            Orbit &oc = Z.modelspace->GetOrbit(c);
            double jc = 0.5 * oc.j2;

            for (int ch_ab = 0; ch_ab < nch; ch_ab++)
            {
              auto &tbc_ab = X.modelspace->GetTwoBodyChannel(ch_ab);
              int Jab = tbc_ab.J;
              size_t nkets_ab = tbc_ab.GetNumberKets();
              for (size_t iket_ab = 0; iket_ab < nkets_ab; iket_ab++)
              {
                Ket &ket_ab = tbc_ab.GetKet(iket_ab);
                int a = ket_ab.p;
                int b = ket_ab.q;
                Orbit &oa = Z.modelspace->GetOrbit(a);
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double occfactor = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
                if (std::abs(occfactor) < 1e-6)
                  continue;
                if (a == b)
                  occfactor *= 0.5; // we sum a<=b, and drop the 1/2, but we still need the 1/2 for a==b

                // P_ij
                bool xcjab_good = ((oj.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(oj.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                // Xcjab term  // J3 = Jab
                if ((xcjab_good) and (std::abs(oj.j2 - oc.j2) <= 2 * Jab) and (oj.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoj1_min = std::abs(oi.j2 - 2 * Jab);
                  int twoj1_max = oi.j2 + 2 * Jab;
                  int twoj2_min = std::abs(oc.j2 - 2 * J2);
                  int twoj2_max = oc.j2 + 2 * J2;
                  double xcjab = X2.GetTBME_J(Jab, c, j, a, b);

                  for (int twoj1 = twoj1_min; twoj1 <= twoj1_max; twoj1 += 2)
                  {
                    for (int twoj2 = twoj2_min; twoj2 <= twoj2_max; twoj2 += 2)
                    {
                      int phasefactor = Z.modelspace->phase((oc.j2 + twoj2) / 2 + J1 + Lambda);
                      double hatfactor = sqrt((2 * Jab + 1.) * (2 * J1 + 1) * (twoj1 + 1) * ( twoj2 + 1));
                      double yabiklc = Y3.GetME_pn(Jab, twoj1, J2, twoj2, a, b, i, k, l, c);
                      double sixj = AngMom::SixJ(oi.j2 / 2.,     oj.j2 / 2.,      J1,
                                                 oc.j2 / 2.,     twoj1 / 2.,      Jab);
                      sixj *=       AngMom::SixJ(Lambda,         J1,              J2,
                                                 oc.j2 / 2.,     twoj2 / 2.,      twoj1 / 2.);
                      zijkl -= occfactor * hatfactor * phasefactor * sixj * xcjab * yabiklc;
                    } // j2
                  } // j1
                }

                bool xciab_good = ((oi.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(oi.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                // xciab   // J3 = Jab
                if ((xciab_good) and (std::abs(oi.j2 - oc.j2) <= 2 * Jab) and (oi.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoj1_min = std::abs(oj.j2 - 2 * Jab);
                  int twoj1_max = oj.j2 + 2 * Jab;
                  int twoj2_min = std::abs(oc.j2 - 2 * J2);
                  int twoj2_max = oc.j2 + 2 * J2;
                  double xciab = X2.GetTBME_J(Jab, c, i, a, b);
                  for (int twoj1 = twoj1_min; twoj1 <= twoj1_max; twoj1 += 2)
                  {
                    for (int twoj2 = twoj2_min; twoj2 <= twoj2_max; twoj2 += 2)
                    {
                      int phasefactor = Z.modelspace->phase((oi.j2 + oj.j2 + oc.j2 + twoj2) / 2  + Lambda);
                      double hatfactor = sqrt((2 * Jab + 1.) * (2 * J1 + 1) * (twoj1 + 1) * ( twoj2 + 1));
                      double yabjklc = Y3.GetME_pn(Jab, twoj1, J2, twoj2, a, b, j, k, l, c);
                      double sixj = AngMom::SixJ(oj.j2 / 2.,     oi.j2 / 2.,      J1,
                                                 oc.j2 / 2.,     twoj1 / 2.,      Jab);
                      sixj *=       AngMom::SixJ(Lambda,         J1,              J2,
                                                 oc.j2 / 2.,     twoj2 / 2.,      twoj1 / 2.);
                      zijkl += occfactor * hatfactor * phasefactor * sixj * xciab * yabjklc;
                    } // j2
                  } // j1
                }

                bool ycjab_good = ((oj.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(oj.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                // ycjab   // J4 = Jab
                if ((ycjab_good) )
                {
                  int twoj1_min = std::max({std::abs(oi.j2 - 2 * Jab), std::abs(oc.j2 - 2 * J2)});
                  int twoj1_max = std::min({oi.j2 + 2 * Jab, oc.j2 + 2 * J2});

                  for (int twoj1 = twoj1_min; twoj1 <= twoj1_max; twoj1 += 2)
                  {
                    double xabiklc = X3.GetME_pn(Jab, J2, twoj1, a, b, i, k, l, c);
                    int twoj2_min = std::max({std::abs(oc.j2 - 2 * Jab), std::abs(oi.j2 - 2 * J2), std::abs(oj.j2 - 2 * Lambda)});
                    int twoj2_max = std::min({oc.j2 + 2 * Jab, oi.j2 + 2 * J2, oj.j2 + 2 * Lambda});
                    for (int twoj2 = twoj2_min; twoj2 <= twoj2_max; twoj2 += 2)
                    {
                      int J3min = std::abs(oj.j2 - oc.j2) / 2;
                      int J3max = (oj.j2 + oc.j2) / 2;
                      for (int J3 = J3min; J3 <= J3max; J3++ )
                      {
                        int phasefactor = Z.modelspace->phase((oi.j2 + oc.j2) / 2 + J1 + J3);
                        double hatfactor = sqrt((2 * J3 + 1.) * (2 * J1 + 1) ) * (twoj1 + 1) * ( twoj2 + 1);
                        double ycjab = Y2.GetTBME_J(J3, Jab, c, j, a, b);
              
                        double sixj = AngMom::SixJ(Lambda,         Jab,            J3,
                                                  oc.j2 / 2.,     oj.j2 / 2.,      twoj2 / 2.);
                        sixj *=       AngMom::SixJ(oc.j2 / 2.,    twoj1 / 2.,      J2,
                                                   oi.j2 / 2.,    twoj2 / 2.,      Jab);
                        sixj *=       AngMom::SixJ(J1,            Lambda,          J2,
                                                   twoj2 / 2.,    oi.j2 / 2.,      oj.j2 / 2.);
                        zijkl -= occfactor * hatfactor * phasefactor * sixj * ycjab * xabiklc;
                      } // J3
                    } // j2
                  } // j1
                }

                bool yciab_good = ((oi.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(oi.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                // yciab   // J4 = Jab  (ci)^J3
                if ((yciab_good) )
                {
                  int twoj1_min = std::max({std::abs(oj.j2 - 2 * Jab), std::abs(oc.j2 - 2 * J2)});
                  int twoj1_max = std::min({oj.j2 + 2 * Jab, oc.j2 + 2 * J2});
                  for (int twoj1 = twoj1_min; twoj1 <= twoj1_max; twoj1 += 2)
                  {
                    double xabjklc = X3.GetME_pn(Jab, J2, twoj1, a, b, j, k, l, c);
                    int twoj2_min = std::max({std::abs(oc.j2 - 2 * Jab), std::abs(oj.j2 - 2 * J2), std::abs(oi.j2 - 2 * Lambda)});
                    int twoj2_max = std::min({oc.j2 + 2 * Jab, oj.j2 + 2 * J2, oi.j2 + 2 * Lambda});
                    for (int twoj2 = twoj2_min; twoj2 <= twoj2_max; twoj2 += 2)
                    {
                      int J3min = std::abs(oi.j2 - oc.j2) / 2;
                      int J3max = (oi.j2 + oc.j2) / 2;
                      for (int J3 = J3min; J3 <= J3max; J3++ )
                      {
                        if (i == c and J3 % 2 == 1 )
                          continue;
                      
                        int phasefactor = Z.modelspace->phase((oi.j2 + oc.j2) / 2 + J3);
                        double hatfactor = sqrt((2 * J3 + 1.) * (2 * J1 + 1) ) * (twoj1 + 1) * ( twoj2 + 1);
                        double yciab = Y2.GetTBME_J(J3, Jab, c, i, a, b);
              
                        double sixj = AngMom::SixJ(Lambda,        Jab,             J3,
                                                  oc.j2 / 2.,     oi.j2 / 2.,      twoj2 / 2.);

                               sixj *=AngMom::SixJ(oc.j2 / 2.,    twoj1 / 2.,      J2,
                                                   oj.j2 / 2.,    twoj2 / 2.,      Jab);

                               sixj *=AngMom::SixJ(J1,            Lambda,          J2,
                                                   twoj2 / 2.,    oj.j2 / 2.,      oi.j2 / 2.);
                        zijkl -= occfactor * hatfactor * phasefactor * sixj * yciab * xabjklc;
                      } // J3
                    } // j2
                  } // j1
                }

                // P_kl
                bool xabcl_good = ((ol.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(ol.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                // Xabcl term  // J3 = Jab
                if ((xabcl_good) and (std::abs(ol.j2 - oc.j2) <= 2 * Jab) and (ol.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoj1_min = std::abs(oc.j2 - 2 * J1);
                  int twoj1_max = oc.j2 + 2 * J1;
                  int twoj2_min = std::abs(ok.j2 - 2 * Jab);
                  int twoj2_max = ok.j2 + 2 * Jab;
                  double xabcl = X2.GetTBME_J(Jab, a, b, c, l);
                  for (int twoj1 = twoj1_min; twoj1 <= twoj1_max; twoj1 += 2)
                  {
                    for (int twoj2 = twoj2_min; twoj2 <= twoj2_max; twoj2 += 2)
                    {
                      int phasefactor = Z.modelspace->phase((oc.j2 + twoj2) / 2 + J1 + Lambda);
                      double hatfactor = sqrt((2 * Jab + 1.) * (2 * J2 + 1) * (twoj1 + 1) * ( twoj2 + 1));
                      double yijcabk = Y3.GetME_pn(J1, twoj1, Jab, twoj2, i, j, c, a, b, k);

                      double sixj = AngMom::SixJ(J2,             Lambda,          J1,
                                                 twoj1 / 2.,     oc.j2 / 2.,      twoj2 / 2.);
                      sixj *=       AngMom::SixJ(ok.j2 / 2.,     ol.j2 / 2.,      J2,
                                                 oc.j2 / 2.,     twoj2 / 2.,      Jab);
                      zijkl += occfactor * hatfactor * phasefactor * sixj * xabcl * yijcabk;
                    } // j2
                  } // j1
                }

                bool xabck_good = ((ok.l + oc.l + tbc_ab.parity) % 2 == X.parity) and (std::abs(ok.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * X.rank_T);
                // Xabck term  // J3 = Jab
                if ((xabck_good) and (std::abs(ok.j2 - oc.j2) <= 2 * Jab) and (ok.j2 + oc.j2 >= 2 * Jab))
                {
                  int twoj1_min = std::abs(oc.j2 - 2 * J1);
                  int twoj1_max = oc.j2 + 2 * J1;
                  int twoj2_min = std::abs(ol.j2 - 2 * Jab);
                  int twoj2_max = ol.j2 + 2 * Jab;
                  double xabck = X2.GetTBME_J(Jab, a, b, c, k);
                  for (int twoj1 = twoj1_min; twoj1 <= twoj1_max; twoj1 += 2)
                  {
                    for (int twoj2 = twoj2_min; twoj2 <= twoj2_max; twoj2 += 2)
                    {
                      int phasefactor = Z.modelspace->phase((oc.j2 + twoj2 + ok.j2 + ol.j2) / 2 + J1 + J2 + Lambda);
                      double hatfactor = sqrt((2 * Jab + 1.) * (2 * J2 + 1) * (twoj1 + 1) * ( twoj2 + 1));
                      double yijcabkl = Y3.GetME_pn(J1, twoj1, Jab, twoj2, i, j, c, a, b, l);

                      double sixj = AngMom::SixJ(J2,             Lambda,          J1,
                                                 twoj1 / 2.,     oc.j2 / 2.,      twoj2 / 2.);
                      sixj *=       AngMom::SixJ(ol.j2 / 2.,     ok.j2 / 2.,      J2,
                                                 oc.j2 / 2.,     twoj2 / 2.,      Jab);
                      zijkl -= occfactor * hatfactor * phasefactor * sixj * xabck * yijcabkl;
                    } // j2
                  } // j1
                }

                bool yabcl_good = ((ol.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(ol.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                // yabcl   // J3 = Jab   (cl)^J4
                if ((yabcl_good) )
                {
                  int twoj1_min = std::max(std::abs(oc.j2 - 2 * J1), std::abs(ok.j2 - 2 * Jab));
                  int twoj1_max = std::min(oc.j2 + 2 * J1, ok.j2 + 2 * Jab);
                  for (int twoj1 = twoj1_min; twoj1 <= twoj1_max; twoj1 += 2)
                  {
                    double xijcabk = X3.GetME_pn(J1, Jab, twoj1, i, j, c, a, b, k);
                    int twoj2_min = std::max({std::abs(oc.j2 - 2 * Jab), std::abs(ok.j2 - 2 * J1), std::abs(ol.j2 - 2 * Lambda)});
                    int twoj2_max = std::min({oc.j2 + 2 * Jab, ok.j2 + 2 * J1, ol.j2 + 2 * Lambda});
                    for (int twoj2 = twoj2_min; twoj2 <= twoj2_max; twoj2 += 2)
                    {
                      int J4min = std::abs(ol.j2 - oc.j2) / 2;
                      int J4max = (ol.j2 + oc.j2) / 2;
                      for (int J4 = J4min; J4 <= J4max; J4 += 1)
                      {
                        int phasefactor = Z.modelspace->phase((ok.j2 + oc.j2) / 2 + J1 + Jab);
                        double hatfactor = sqrt((2 * J2 + 1.) * (2 * J4 + 1) ) * (twoj1 + 1) * ( twoj2 + 1);
                        double yabcl = Y2.GetTBME_J(Jab, J4, a, b, c, l);
                        double sixj = AngMom::SixJ(oc.j2 / 2.,    twoj1 / 2.,      J1,
                                                   ok.j2 / 2.,    twoj2 / 2.,      Jab);
                        sixj *=       AngMom::SixJ(Lambda,        Jab,             J4,
                                                   oc.j2 / 2.,    ol.j2 / 2.,      twoj2 / 2.);
                        sixj *=       AngMom::SixJ(J2,            J1,              Lambda,
                                                   twoj2 / 2.,    ol.j2 / 2.,      ok.j2 / 2.);
                        zijkl += occfactor * hatfactor * phasefactor * sixj * yabcl * xijcabk;
                      } // J3
                    } // j2
                  } // j1
                }

                bool yabck_good = ((ok.l + oc.l + tbc_ab.parity) % 2 == Y.parity) and (std::abs(ok.tz2 + oc.tz2 - 2 * tbc_ab.Tz) == 2 * Y.rank_T);
                // yabck   // J3 = Jab   (ck)^J4
                if ((yabck_good) )
                {
                  int twoj1_min = std::max(std::abs(oc.j2 - 2 * J1), std::abs(ol.j2 - 2 * Jab));
                  int twoj1_max = std::min(oc.j2 + 2 * J1, ol.j2 + 2 * Jab);
                  for (int twoj1 = twoj1_min; twoj1 <= twoj1_max; twoj1 += 2)
                  {
                    double xijcabl = X3.GetME_pn(J1, Jab, twoj1, i, j, c, a, b, l);
                    int twoj2_min = std::max({std::abs(oc.j2 - 2 * Jab), std::abs(ol.j2 - 2 * J1), std::abs(ok.j2 - 2 * Lambda)});
                    int twoj2_max = std::min({oc.j2 + 2 * Jab, ol.j2 + 2 * J1, ok.j2 + 2 * Lambda});
                    for (int twoj2 = twoj2_min; twoj2 <= twoj2_max; twoj2 += 2)
                    {
                      int J4min = std::abs(ok.j2 - oc.j2) / 2;
                      int J4max = (ok.j2 + oc.j2) / 2;
                      for (int J4 = J4min; J4 <= J4max; J4 += 1)
                      {
                        int phasefactor = Z.modelspace->phase((ok.j2 + oc.j2) / 2 + J1 + J2 + Jab);
                        double hatfactor = sqrt((2 * J2 + 1.) * (2 * J4 + 1) ) * (twoj1 + 1) * ( twoj2 + 1);
                        double yabck = Y2.GetTBME_J(Jab, J4, a, b, c, k);
                        double sixj = AngMom::SixJ(oc.j2 / 2.,    twoj1 / 2.,      J1,
                                                   ol.j2 / 2.,    twoj2 / 2.,      Jab);
                        sixj *=       AngMom::SixJ(Lambda,        Jab,             J4,
                                                   oc.j2 / 2.,    ok.j2 / 2.,      twoj2 / 2.);
                        sixj *=       AngMom::SixJ(J2,            J1,              Lambda,
                                                   twoj2 / 2.,    ok.j2 / 2.,      ol.j2 / 2.);
                        zijkl += occfactor * hatfactor * phasefactor * sixj * yabck * xijcabl;
                      } // J3
                    } // j2
                  } // j1
                }
              } // for iket_ab
            } // for ch2 J3
          } // for c

          // normalize the tbme
          zijkl *= -1.0 / sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
        } // for iket
      } // for ibra
    } // for ch

  } // comm232st

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Z^(J1j1, J2j2)Lamda_ijklmn =   PJ3(ij/k) sum_a  ( X_ka Y^(J1j1, J2j2)Lamda_ijalmn 
  ///                                                      -(-)^(J1+Lambda+j1+ja)sqrt(2 j1 + 1) sqrt(2 j2 + 1) 
  ///                                                           { jk ja Lamda} Y^Lamda_ka X^(J1j2, J2j2)0_ijalmn)
  ///                                                           { j2 j1 J1   }
  ///                                                                                 
  ///                                             - PJ3(lm/n) sum_a ( Y^(J1j1, J2j2)Lamda_ijklma X_an 
  ///                                                       -(-)^(J2+Lambda+j1+jn) sqrt(2 j1 + 1) sqrt(2 j2 + 1) 
  ///                                                           { jn ja Lamda} X^(J1j1,J2j1)0_ijklma Y^Lamda_an )
  ///                                                           { j1 j2 J2   }
  ///
  ///
  void comm133st(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    auto &Z3 = Z.ThreeBody;
    int Lambda = Z.GetJRank();
    Z.modelspace->PreCalculateSixJ();

    std::vector<std::array<size_t, 3>> bra_ket_channels;
    for (auto &it : Z.ThreeBody.Get_ch_start())
    {
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(it.first.ch_bra);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        bra_ket_channels.push_back({it.first.ch_bra, it.first.ch_ket, static_cast<size_t>(ibra)}); // (ch_bra, ch_ket,ibra)
      }
    }
    size_t n_bra_ket_ch = bra_ket_channels.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ich = 0; ich < n_bra_ket_ch; ich++)
    {
      size_t ch_bra = bra_ket_channels[ich][0];
      size_t ch_ket = bra_ket_channels[ich][1];
      size_t ibra = bra_ket_channels[ich][2];
      ThreeBodyChannel &Tbc_bra = Z.modelspace->GetThreeBodyChannel(ch_bra);
      ThreeBodyChannel &Tbc_ket = Z.modelspace->GetThreeBodyChannel(ch_ket);

      int twoj1 = Tbc_bra.twoJ;
      int twoj2 = Tbc_ket.twoJ;
      // std::cout<< ch_bra << "   " << ch_ket << "   "<< twoj1 << "   "<< twoj2 <<"   "<< Tbc_bra.parity <<"   "<< Tbc_ket.parity<<"   "<< Tbc_bra.twoTz <<"   "<< Tbc_ket.twoTz << std::endl;

      size_t nbras = Tbc_bra.GetNumberKets();
      size_t nkets = Tbc_ket.GetNumberKets();
      Ket3 &bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      Orbit &ok = Z.modelspace->GetOrbit(k);
      int Jij = bra.Jpq;
      size_t ket_min = (ch_bra == ch_ket) ? ibra : 0;
      for (size_t iket = ket_min; iket < nkets; iket++)
      {
        Ket3 &ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit &ol = Z.modelspace->GetOrbit(l);
        Orbit &om = Z.modelspace->GetOrbit(m);
        Orbit &on = Z.modelspace->GetOrbit(n);
        int Jlm = ket.Jpq;

        double zsum = 0;
        // First, connect on the bra side
        for (auto a : X.GetOneBodyChannel(ok.l, ok.j2, ok.tz2))
        {
          // eq.(9.1)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != ok.j2)
            continue;
          zsum += X1(k, a) * Y3.GetME_pn(Jij, twoj1, Jlm, twoj2, i, j, a, l, m, n);
        }
        for (auto a : Y.GetOneBodyChannel(ok.l, ok.j2, ok.tz2))
        {          
          // eq.(9.2)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != ok.j2)
            continue;

          int phasefactor = Z.modelspace->phase((oa.j2 + twoj1) / 2 + Jij + Lambda);
          double sixj = AngMom::SixJ(ok.j2 / 2.,     oa.j2 / 2.,      Lambda,
                                     twoj2 / 2.,     twoj1 / 2.,      Jij);
          double hatfactor = sqrt( (twoj1 + 1) * ( twoj2 + 1) );
          zsum -= phasefactor * hatfactor * sixj * Y1(k, a) * X3.GetME_pn(Jij, Jlm, twoj2, i, j, a, l, m, n);
        }
        for (auto a : X.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
        {
          // eq.(9.3)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != oi.j2)
            continue;
          int J3min = std::abs(ok.j2 - oj.j2) / 2;
          int J3max = (ok.j2 + oj.j2) / 2;
          for (int J3 = J3min; J3 <= J3max; J3++ )
          {
            double hatfactor = sqrt( ( 2 * Jij + 1) * ( 2 * J3 + 1) );
            double sixj = AngMom::SixJ(oi.j2 / 2.,     oj.j2 / 2.,      Jij,
                                       ok.j2 / 2.,     twoj1 / 2.,      J3);
            zsum += sixj * hatfactor * X1(i, a) * Y3.GetME_pn(J3, twoj1, Jlm, twoj2, k, j, a, l, m, n);
          }
        }
        
        // for( auto a : X.modelspace->all_orbits )
        for (auto a : Y.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
        {
          // eq.(9.4)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != oi.j2)
            continue;

          int J3min = std::abs(ok.j2 - oj.j2) / 2;
          int J3max = (ok.j2 + oj.j2) / 2;
          for (int J3 = J3min; J3 <= J3max; J3++ )
          {
            int phasefactor = Z.modelspace->phase((oa.j2 + twoj1) / 2 + J3 + Lambda);
            double hatfactor = sqrt( ( 2 * Jij + 1) * ( 2 * J3 + 1) * (twoj1 + 1) * ( twoj2 + 1) );
            double sixj = AngMom::SixJ(twoj2 / 2.,     oa.j2 / 2.,      J3,
                                       oi.j2 / 2.,     twoj1 / 2.,      Lambda);
            sixj       *= AngMom::SixJ(oi.j2 / 2.,     oj.j2 / 2.,      Jij,
                                       ok.j2 / 2.,     twoj1 / 2.,      J3);
            zsum -= phasefactor * hatfactor * sixj * Y1(i, a) * X3.GetME_pn(J3, Jlm, twoj2, k, j, a, l, m, n);
          }
        }
        for (auto a : X.GetOneBodyChannel(oj.l, oj.j2, oj.tz2))
        {
          // eq.(9.5)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != oj.j2)
            continue;
          int J3min = std::abs(ok.j2 - oi.j2) / 2;
          int J3max = (ok.j2 + oi.j2) / 2;
          for (int J3 = J3min; J3 <= J3max; J3++ )
          {
            int phasefactor = Z.modelspace->phase((oj.j2 + ok.j2) / 2 + Jij + J3);
            double hatfactor = sqrt( ( 2 * Jij + 1) * ( 2 * J3 + 1) );

            double sixj = AngMom::SixJ(oj.j2 / 2.,     oi.j2 / 2.,      Jij,
                                       ok.j2 / 2.,     twoj1 / 2.,      J3);
            
            zsum -= phasefactor * hatfactor * sixj * X1(j, a) * Y3.GetME_pn(J3, twoj1, Jlm, twoj2, i, k, a, l, m, n);
          }
        }
        for (auto a : Y.GetOneBodyChannel(oj.l, oj.j2, oj.tz2))
        {
          // eq.(9.6)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != oj.j2)
            continue;

          int J3min = std::abs(ok.j2 - oi.j2) / 2;
          int J3max = (ok.j2 + oi.j2) / 2;
          int phasefactor = Z.modelspace->phase((oa.j2 + ok.j2 + oj.j2 + twoj1) / 2 + Jij + Lambda);
          for (int J3 = J3min; J3 <= J3max; J3++ )
          {
            double hatfactor = sqrt( ( 2 * Jij + 1) * ( 2 * J3 + 1) * (twoj1 + 1) * ( twoj2 + 1) );
            double sixj = AngMom::SixJ(twoj2 / 2.,     oa.j2 / 2.,      J3,
                                       oj.j2 / 2.,     twoj1 / 2.,      Lambda);
            sixj       *= AngMom::SixJ(oj.j2 / 2.,     oi.j2 / 2.,      Jij,
                                       ok.j2 / 2.,     twoj1 / 2.,      J3);
            zsum += phasefactor * hatfactor * sixj * Y1(j, a) * X3.GetME_pn(J3, Jlm, twoj2, i, k, a, l, m, n);
          }
        }

        // Now connect on the ket side
        for (auto a : X.GetOneBodyChannel(on.l, on.j2, on.tz2))
        {
          // eq.(10.1)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != on.j2)
            continue;
          zsum -= Y3.GetME_pn(Jij, twoj1, Jlm, twoj2, i, j, k, l, m, a) * X1(a, n);
        }
        for (auto a : Y.GetOneBodyChannel(on.l, on.j2, on.tz2))
        {
          // eq.(10.2)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != on.j2)
            continue;
          int phasefactor = Z.modelspace->phase(( on.j2 + twoj1 ) / 2 + Jlm + Lambda);
          double hatfactor = sqrt( (twoj1 + 1) * ( twoj2 + 1) );
          double sixj = AngMom::SixJ(on.j2 / 2.,     oa.j2 / 2.,      Lambda,
                                     twoj1 / 2.,     twoj2 / 2.,      Jlm);
          zsum += phasefactor * hatfactor * sixj * X3.GetME_pn(Jij, Jlm, twoj1, i, j, k, l, m, a) * Y1(a, n);
        }
        for (auto a : X.GetOneBodyChannel(ol.l, ol.j2, ol.tz2))
        {
          // eq.(10.3)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != ol.j2)
            continue;
          int J3min = std::abs(on.j2 - om.j2) / 2;
          int J3max = (on.j2 + om.j2) / 2;
          for (int J3 = J3min; J3 <= J3max; J3++ )
          {
            double hatfactor = sqrt( (2 * Jlm + 1) * ( 2 * J3 + 1) );
            double sixj = AngMom::SixJ(ol.j2 / 2.,     om.j2 / 2.,      Jlm,
                                       on.j2 / 2.,     twoj2 / 2.,      J3);
            zsum -= hatfactor * sixj * Y3.GetME_pn(Jij, twoj1, J3, twoj2, i, j, k, n, m, a) * X1(a, l);
          }
        }
        for (auto a : Y.GetOneBodyChannel(ol.l, ol.j2, ol.tz2))
        {
          // eq.(10.4)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != ol.j2)
            continue;
          int J3min = std::abs(on.j2 - om.j2) / 2;
          int J3max = (on.j2 + om.j2) / 2;
          for (int J3 = J3min; J3 <= J3max; J3++ )
          {
            if ( oa.j2 + ol.j2 < Lambda * 2 and std::abs( oa.j2 - ol.j2 ) > Lambda * 2 )
              continue;
            int phasefactor = Z.modelspace->phase(( ol.j2 + twoj1 ) / 2 + J3 + Lambda);
            double hatfactor = sqrt( (twoj1 + 1) * ( twoj2 + 1) * ( 2 * J3 + 1) * ( 2 * Jlm + 1) );
            double sixj = AngMom::SixJ(on.j2 / 2.,     om.j2 / 2.,      J3,
                                       ol.j2 / 2.,     twoj2 / 2.,      Jlm);
            sixj       *= AngMom::SixJ(ol.j2 / 2.,     oa.j2 / 2.,      Lambda,
                                       twoj1 / 2.,     twoj2 / 2.,      J3);
            zsum += phasefactor * hatfactor * sixj * X3.GetME_pn(Jij, J3, twoj1, i, j, k, n, m, a) * Y1(a, l);
          }
        }
        for (auto a : X.GetOneBodyChannel(om.l, om.j2, om.tz2))
        {
          // eq.(10.5)
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != om.j2)
            continue;
          int J3min = std::abs(ol.j2 - on.j2) / 2;
          int J3max = (ol.j2 + on.j2) / 2;
          for (int J3 = J3min; J3 <= J3max; J3++ )
          {
            int phasefactor = Z.modelspace->phase(( om.j2 + on.j2 ) / 2 + J3 + Jlm);
            double hatfactor = sqrt( ( 2 * J3 + 1) * ( 2 * Jlm + 1) );
            double sixj = AngMom::SixJ(om.j2 / 2.,     ol.j2 / 2.,      Jlm,
                                       on.j2 / 2.,     twoj2 / 2.,      J3);
            zsum += phasefactor * hatfactor * sixj * Y3.GetME_pn(Jij, twoj1, J3, twoj2, i, j, k, l, n, a) * X1(a, m);
          }
        }
        for (auto a : Y.GetOneBodyChannel(om.l, om.j2, om.tz2))
        {
          // eq.(10.6)    
          Orbit &oa = Z.modelspace->GetOrbit(a);
          if (oa.j2 != om.j2)
            continue;      
          int J3min = std::abs(ol.j2 - on.j2) / 2;
          int J3max = (ol.j2 + on.j2) / 2;
          for (int J3 = J3min; J3 <= J3max; J3++ )
          {
            if ( oa.j2 + om.j2 < Lambda * 2 and std::abs( oa.j2 - om.j2 ) > Lambda * 2 )
              continue;
            int phasefactor = Z.modelspace->phase(( on.j2 + twoj1 ) / 2 + Jlm + Lambda);
            double hatfactor = sqrt( ( 2 * J3 + 1) * ( 2 * Jlm + 1) *  (twoj1 + 1) * ( twoj2 + 1));
            double sixj = AngMom::SixJ(on.j2 / 2.,     ol.j2 / 2.,      J3,
                                       om.j2 / 2.,     twoj2 / 2.,      Jlm);
            sixj       *= AngMom::SixJ(om.j2 / 2.,     oa.j2 / 2.,      Lambda,
                                       twoj1 / 2.,     twoj2 / 2.,      J3);                           
            zsum += phasefactor * hatfactor * sixj * X3.GetME_pn(Jij, J3, twoj1, i, j, k, l, n, a) * Y1(a, m);
          }
        }

        Z3.AddToME_pn_ch(ch_bra, ch_ket, ibra, iket, zsum);

      } // for iket
    } // for ich  -> {ch_bra,ch_ket,ibra}

  } // comm133st


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Z^(J1,J2)Lamda_ijkl =  sum_{a,b,j1,j2} (na-nb) (-)^(J1+ja+j2+Lamda) Sqrt[(2j1+1) (2j2+1)]
  ///                                                       { Lamda  J1     J2 } A_ab B^(J1j1,J2j2)Lamda_ijbkla
  ///                                                       { ja     j2     j1 }
  //                                        sum_{a,b,j1}    (na-nb) (-)^(J2+ja+j1) (2j1+1)
  ///                                                       { J2     Lamda  J1 } B_ab A^(J1j1,J2j1)0_ijbkla
  ///                                                       { jb     j1     ja }
  ///
  void comm132st(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;
    int Lambda = Z.GetJRank();
    Z.modelspace->PreCalculateSixJ();

    std::vector<size_t> ch_bra_list, ch_ket_list;

    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      //      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch2);
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
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
            Orbit &oa = Z.modelspace->GetOrbit(a);
            for (size_t b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              if ( oa.j2 != ob.j2 )
                continue;
              
              int twoj1min = std::abs(2 * J1 - ob.j2);
              int twoj1max = 2 * J1 + ob.j2;
              int twoj2min = std::abs(2 * J2 - oa.j2);
              int twoj2max = 2 * J2 + oa.j2;
              double xab = X1(a, b);
              if (std::abs(xab) > 1e-8)
                for (int twoj1 = twoj1min; twoj1 <= twoj1max; twoj1 += 2)
                {
                  for (int twoj2 = twoj2min; twoj2 <= twoj2max; twoj2 += 2)
                  {
                    int phasefactor = Z.modelspace->phase((oa.j2 + twoj2) / 2 + J1 + Lambda);
                    double hatfactor = sqrt( (twoj1 + 1) * ( twoj2 + 1) );
                    double sixj = AngMom::SixJ( Lambda    ,     J1,              J2,
                                                oa.j2 / 2.,     twoj2 / 2.,      twoj1 / 2.);
                    double yijbkla = Y3.GetME_pn(J1, twoj1, J2, twoj2, i, j, b, k, l, a);
                    zijkl += phasefactor * hatfactor * sixj * (oa.occ - ob.occ) * (xab * yijbkla );
                  } // twoj2
                } // twoj1


              twoj1min = std::max(std::abs(2 * J1 - ob.j2), std::abs(2 * J2 - oa.j2));
              twoj1max = std::min(2 * J1 + ob.j2, 2 * J2 + oa.j2);
              double yab = Y1(a, b);
              if (std::abs(yab) > 1e-8)
                for (int twoj1 = twoj1min; twoj1 <= twoj1max; twoj1 += 2)
                {
                  double xijbkla = X3.GetME_pn(J1, J2, twoj1, i, j, b, k, l, a);
                  if (std::abs(xijbkla) < 1e-8)
                    continue;
                  int phasefactor = Z.modelspace->phase((oa.j2 + twoj1) / 2 + J2);
                  double hatfactor = (twoj1 + 1);
                  double sixj = AngMom::SixJ( J2    ,         Lambda,          J1,
                                              ob.j2 / 2.,     twoj1 / 2.,      oa.j2 / 2.);
                  zijkl -= phasefactor * hatfactor * sixj * (oa.occ - ob.occ) * (yab * xijbkla);
                } // twoj1
            } // b
          } // a
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
        } // iket
      } // ibra
    } // ch2
  } // comm132st


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    Z^(J1,J2)Lamda_ijkl = 1/6 sum_abcd sum_{J3,j1,j2} (-)^( J1+Lamda+jd+j2 ) sqrt( (2j1+1) * (2j2+1) )
  ///                                       *  (na nb nc(1-nd) + (1-na)(1-nb)(1-nc)nd)  { J1 J2 Lamda }
  ///                                                                                   { j2 j1 jd    }
  ///                                       *  ( X^(J1j1,J3j1)0_ijdabc Y^(J3j1,J2j2)Lamda_abckld  
  ///                                          - Y^(J1j1,J3j2)Lamda_ijdabc X^(J3j2,J2j2)0_abckld )
  ///
  void comm332_ppph_hhhpst(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;
    int Lambda = Z.GetJRank();
    Z.modelspace->PreCalculateSixJ();

    std::vector<int> bra_channels;
    std::vector<int> ket_channels;
    for (auto &itmat : Z.TwoBody.MatEl)
    {
      bra_channels.push_back(itmat.first[0]);
      ket_channels.push_back(itmat.first[1]);
    }
    int nmat = bra_channels.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch2 = 0; ch2 < nmat; ch2++)
    {
      int ch_bra = bra_channels[ch2];
      int ch_ket = ket_channels[ch2];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        for (int iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;
          double zijkl = 0;
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
                  double occfactor = oa.occ * ob.occ * oc.occ * (1 - od.occ) - (1 - oa.occ) * (1 - ob.occ) * (1 - oc.occ) * od.occ;
                  if (std::abs(occfactor) < 1e-7)
                    continue;
                  //if ((oi.l + oj.l + od.l + oa.l + ob.l + oc.l) % 2 > 0)
                  //  continue;
                  //if ((oi.tz2 + oj.tz2 + od.tz2) - (oa.tz2 + ob.tz2 + oc.tz2) != 0)
                  //  continue;
                  int J3min = std::abs(oa.j2 - ob.j2) / 2;
                  int J3max = (oa.j2 + ob.j2) / 2;
                  for (int J3 = J3min; J3 <= J3max; J3++)
                  {
                    int twoj1min = std::abs(2 * J1 - od.j2);
                    int twoj1max = 2 * J1 + od.j2;
                    int twoj2min = std::abs(2 * J2 - od.j2);
                    int twoj2max = 2 * J2 + od.j2;
                    for (int twoj1 = twoj1min; twoj1 <= twoj1max; twoj1 += 2)
                    {
                      for (int twoj2 = twoj2min; twoj2 <= twoj2max; twoj2 += 2)
                      {
                        int phasefactor = Z.modelspace->phase((od.j2 + twoj2) / 2 + J1 + Lambda);
                        double hatfactor = sqrt( (twoj1 + 1) * ( twoj2 + 1) );
                        double sixj = AngMom::SixJ(J1,            J2,              Lambda,
                                                   twoj2 / 2.,    twoj1 / 2.,      od.j2 / 2.);
                        // Eq.(7.1)
                        if (  twoj1 >= std::abs(2 * J3 - oc.j2) and twoj1 <= (2 * J3 + oc.j2) )
                        {
                          double xijdabc = X3.GetME_pn(J1, J3, twoj1, i, j, d, a, b, c);
                          double yabckld = Y3.GetME_pn(J3, twoj1, J2, twoj2, a, b, c, k, l, d);
                          zijkl += 1. / 6 * occfactor * phasefactor * hatfactor * sixj * (xijdabc * yabckld);
                        }
                        // Eq.(7.2)
                        if (  twoj2 >= std::abs(2 * J3 - oc.j2) and twoj2 <= (2 * J3 + oc.j2) )
                        {
                          double yijdabc = Y3.GetME_pn(J1, twoj1, J3, twoj2, i, j, d, a, b, c);
                          double xabckld = X3.GetME_pn(J3, J2, twoj2, a, b, c, k, l, d);
                          zijkl -= 1. / 6 * occfactor * phasefactor * hatfactor * sixj * (yijdabc * xabckld);
                        }
                      } // twoj1
                    } // twoj2
                  } // J2
                } // d
              } // c
            } // b
          } // a
          
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
        } // iket
      } // ibra
    } // ch2

  } // comm332_ppph_hhhpst




/**
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ_ijkn = (1-PJij)(1-PJkl) -sum_Jph (2Jph+1) { i j J   } ZbarJph_ilkj
  ///                                                             { k l Jph }
  ///
  ///  where         ZbarJph_ilkj = 1/4 sum_abcd sum Jab Jcd  ( na nb(1-nc)(1-nd) -  (1-na)(1-nb)nc nd ) XJph_il;(abJab)(cdJcd) YJph_(cdJcd)(abJab);kj
  ///
  ///  and           XJph_il;(abJab)(cdJcd) = -sum_J (-1)^{i+J} {i   l   Jph } XJabJcdJ_abicdl
  ///                                                           {Jcd Jab J   }
  ///
  void comm332_pphhss(const Operator &X, const Operator &Y, Operator &Z)
  {
    auto &X3 = X.ThreeBody;
    auto &Y3 = Y.ThreeBody;
    auto &Z2 = Z.TwoBody;
    std::map<int, double> e_fermi = Z.modelspace->GetEFermi();

    int nch = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar3b_transform_first_pass)
    for (int ch = 0; ch < nch; ch++)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nkets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nkets; ibra++)
      {
        Ket &bra = tbc.GetKet(ibra);
        int i = bra.p;
        int j = bra.q;
        int ji2 = bra.op->j2;
        int jj2 = bra.oq->j2;
        double ji = 0.5 * ji2;
        double jj = 0.5 * jj2;
        double d_ei = std::abs(2 * bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
        double d_ej = std::abs(2 * bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
        double occnat_i = bra.op->occ_nat;
        double occnat_j = bra.oq->occ_nat;
        for (int iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc.GetKet(iket);
          int k = ket.p;
          int l = ket.q;
          int jk2 = ket.op->j2;
          int jl2 = ket.oq->j2;
          double jk = 0.5 * jk2;
          double jl = 0.5 * jl2;
          double d_ek = std::abs(2 * ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
          double d_el = std::abs(2 * ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
          double occnat_k = ket.op->occ_nat;
          double occnat_l = ket.oq->occ_nat;

          double zijkl = 0;
          // Now the loops on the right hand side
          for (int ch_ab = 0; ch_ab < nch; ch_ab++)
          {
            TwoBodyChannel &tbc_ab = Z.modelspace->GetTwoBodyChannel(ch_ab);
            int Jab = tbc_ab.J;
            int nkets = tbc_ab.GetNumberKets();
            for (int iket_ab = 0; iket_ab < nkets; iket_ab++)
            {
              Ket &ket_ab = tbc_ab.GetKet(iket_ab);
              int a = ket_ab.p;
              int b = ket_ab.q;
              double na = ket_ab.op->occ;
              double nb = ket_ab.oq->occ;

              for (int ch_cd = 0; ch_cd < nch; ch_cd++)
              {
                TwoBodyChannel &tbc_cd = Z.modelspace->GetTwoBodyChannel(ch_cd);
                int Jcd = tbc_cd.J;
                int nkets = tbc_cd.GetNumberKets();
                for (int iket_cd = 0; iket_cd < nkets; iket_cd++)
                {
                  Ket &ket_cd = tbc_cd.GetKet(iket_cd);
                  int c = ket_cd.p;
                  int d = ket_cd.q;
                  double nc = ket_cd.op->occ;
                  double nd = ket_cd.oq->occ;
                  //                double jc = 0.5*ket_cd.op->j2;
                  //                double jd = 0.5*ket_cd.oq->j2;
                  double occupation_factor = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);
                  if (std::abs(occupation_factor) < 1e-6)
                    continue;

                  double symmetry_factor = 1; // we only sum a<=b and c<=d, so we undercount by a factor of 4, canceling the 1/4 in the formula
                  if (a == b)
                    symmetry_factor *= 0.5; // if a==b or c==d, then the permutation doesn't give a new state, so there's less undercounting
                  if (c == d)
                    symmetry_factor *= 0.5;

                  // Figure out which range of twoJp and twoJpp we will need
                  int twoJp_min = std::max(std::min(std::abs(2 * Jab - ji2), std::abs(2 * Jab - jj2)), std::min(std::abs(2 * Jcd - jl2), std::abs(2 * Jcd - jk2)));
                  int twoJpp_min = std::max(std::min(std::abs(2 * Jcd - jj2), std::abs(2 * Jcd - ji2)), std::min(std::abs(2 * Jab - jk2), std::abs(2 * Jab - jl2)));
                  int twoJp_max = std::min(2 * Jab + std::max(ji2, jj2), 2 * Jcd + std::max(jk2, jl2));
                  int twoJpp_max = std::min(2 * Jcd + std::max(ji2, jj2), 2 * Jab + std::max(jk2, jl2));

                  if (twoJpp_max < twoJpp_min)
                    continue;
                  for (int twoJp = twoJp_min; twoJp <= twoJp_max; twoJp += 2)
                  {

                    double xabicdl = X3.GetME_pn(Jab, Jcd, twoJp, a, b, i, c, d, l);
                    double xabicdk = X3.GetME_pn(Jab, Jcd, twoJp, a, b, i, c, d, k);
                    double xabjcdl = X3.GetME_pn(Jab, Jcd, twoJp, a, b, j, c, d, l);
                    double xabjcdk = X3.GetME_pn(Jab, Jcd, twoJp, a, b, j, c, d, k);

                    for (int twoJpp = twoJpp_min; twoJpp <= twoJpp_max; twoJpp += 2)
                    {
                      double Jp = 0.5 * twoJp;
                      double Jpp = 0.5 * twoJpp;
                      double hatfactor = (twoJp + 1) * (twoJpp + 1);

                      // I think having these in the inner loop may be disastrous for performance
                      double ninej1 = Z.modelspace->GetNineJ(Jab, jk, Jpp, ji, J, jj, Jp, jl, Jcd);
                      double ninej2 = Z.modelspace->GetNineJ(Jab, jk, Jpp, jj, J, ji, Jp, jl, Jcd); // permute i<->j
                      double ninej3 = Z.modelspace->GetNineJ(Jab, jl, Jpp, ji, J, jj, Jp, jk, Jcd); // permute k<->l
                      double ninej4 = Z.modelspace->GetNineJ(Jab, jl, Jpp, jj, J, ji, Jp, jk, Jcd); // permute i<->j and k<->l

                      // These phase factors account for the minus signs associated with the permutations
                      // so that all the permuted terms should just be added with their phase (no extra minus sign).
                      int phase1 = Z.modelspace->phase((ji2 + jk2 + twoJp + twoJpp) / 2);
                      int phase2 = Z.modelspace->phase((ji2 + jk2 + twoJp + twoJpp) / 2 - J); // from permuting i<->j
                      int phase3 = Z.modelspace->phase((ji2 + jk2 + twoJp + twoJpp) / 2 - J); // from permuting k<->l
                      int phase4 = Z.modelspace->phase((ji2 + jk2 + twoJp + twoJpp) / 2);     // from permuting i<->j and k<->l

                      double ycdjabk = Y3.GetME_pn(Jcd, Jab, twoJpp, c, d, j, a, b, k);
                      double ycdjabl = Y3.GetME_pn(Jcd, Jab, twoJpp, c, d, j, a, b, l);
                      double ycdiabk = Y3.GetME_pn(Jcd, Jab, twoJpp, c, d, i, a, b, k);
                      double ycdiabl = Y3.GetME_pn(Jcd, Jab, twoJpp, c, d, i, a, b, l);

                      zijkl += symmetry_factor * hatfactor * occupation_factor * (1 * phase1 * ninej1 * xabicdl * ycdjabk + 1 * phase2 * ninej2 * xabjcdl * ycdiabk // i<->j
                                                                                  + 1 * phase3 * ninej3 * xabicdk * ycdjabl                                         // k<->l
                                                                                  + 1 * phase4 * ninej4 * xabjcdk * ycdiabl                                         // i<->j and k<->l
                                                                                 );

                    } // for twoJpp
                  } // for twoJp

                } // for iket_cd
              } // for ch_cd
            } // for iket_ab
          } // for ch_ab
          // make it a normalized TBME
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          // the AddToTBME routine automatically takes care of the hermitian conjugate as well
          Z2.AddToTBME(ch, ch, ibra, iket, zijkl);
        } // for iket
      } // for ibra
    } // for ch
  }
*/





  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// Start double nested commutators
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //      *------*     |p
  //     /\     | \    |       diagrams Ia, and Ia* (which is just the Hermitian conjugate, with p<->q).
  //   a(  )i  b|  )   |
  //     \/     |  |j  |
  //      *~~~~~*  )   |
  //           c| /    |
  //             *-----*
  //                   |
  //                   |q
  //
  void diagram_CIa(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto a : Z.modelspace->particles)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (auto b : Z.modelspace->particles)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        for (auto i : Z.modelspace->holes)
        {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          for (auto j : Z.modelspace->holes)
          {
            Orbit &oj = Z.modelspace->GetOrbit(j);
            for (auto c : X.OneBodyChannels.at({oj.l, oj.j2, oj.tz2}))
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              if ((1 - oc.occ) < 1e-4)
                continue;
              int J1min = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
              int J1max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                double xijab = X.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                double yabic = Y.TwoBody.GetTBME_J(J1, J1, a, b, i, c);
                double xabic = X.TwoBody.GetTBME_J(J1, J1, a, b, i, c);
                for (auto p : Z.modelspace->all_orbits)
                {
                  Orbit &op = Z.modelspace->GetOrbit(p);
                  for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
                  {
                    Orbit &oq = Z.modelspace->GetOrbit(q);
                    double zpq = 0;
                    int J2min = std::max(std::abs(oj.j2 - op.j2), std::abs(oj.j2 - oq.j2)) / 2;
                    int J2max = std::min(oj.j2 + op.j2, oj.j2 + oq.j2) / 2;
                    for (int J2 = J2min; J2 <= J2max; J2++)
                    {
                      double xcpjq = X.TwoBody.GetTBME_J(J2, J2, c, p, j, q);
                      double xcqjp = X.TwoBody.GetTBME_J(J2, J2, c, q, j, p);
                      double ycpjq = Y.TwoBody.GetTBME_J(J2, J2, c, p, j, q);
                      double ycqjp = Y.TwoBody.GetTBME_J(J2, J2, c, q, j, p);
                      zpq -= 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oj.j2 + 1) * (1 - oc.occ) * xijab * yabic * xcpjq;
                      zpq -= 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oj.j2 + 1) * (1 - oc.occ) * xijab * yabic * xcqjp * hX * hX * hY;

                      /// Also include the XXY and YXX contributions
                      zpq += 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oj.j2 + 1) * (1 - oc.occ) * xijab * xabic * ycpjq;
                      zpq += 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oj.j2 + 1) * (1 - oc.occ) * xijab * xabic * ycqjp * hX * hX * hY;
                    }
                    Z.OneBody(p, q) += zpq;
                  }
                }
              }
            }
          }
        }
      }
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //
  //      *------*     |p
  //     /\     | \    |       diagrams Ib, and Ib* (which is just the Hermitian conjugate, with p<->q).
  //   a(  )i  j|  )   |       eventually, this can probably be combined with diagram Ia.
  //     \/     |  |b  |
  //      *~~~~~*  )   |
  //           k| /    |
  //             *-----*
  //                   |
  //                   |q
  //
  void diagram_CIb(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto a : Z.modelspace->particles)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (auto b : Z.modelspace->particles)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        for (auto i : Z.modelspace->holes)
        {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          for (auto j : Z.modelspace->holes)
          {
            Orbit &oj = Z.modelspace->GetOrbit(j);
            for (auto k : X.OneBodyChannels.at({ob.l, ob.j2, ob.tz2}))
            {
              Orbit &ok = Z.modelspace->GetOrbit(k);
              if (ok.occ < 1e-4)
                continue;
              int J1min = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
              int J1max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                double xijab = X.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                double yakij = Y.TwoBody.GetTBME_J(J1, J1, a, k, i, j);
                double xakij = X.TwoBody.GetTBME_J(J1, J1, a, k, i, j);
                for (auto p : Z.modelspace->all_orbits)
                {
                  Orbit &op = Z.modelspace->GetOrbit(p);
                  for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
                  {
                    Orbit &oq = Z.modelspace->GetOrbit(q);
                    double zpq = 0;
                    int J2min = std::max(std::abs(ob.j2 - op.j2), std::abs(ok.j2 - oq.j2)) / 2;
                    int J2max = std::min(ob.j2 + op.j2, ok.j2 + oq.j2) / 2;
                    for (int J2 = J2min; J2 <= J2max; J2++)
                    {
                      double xbpkq = X.TwoBody.GetTBME_J(J2, J2, b, p, k, q);
                      double xbqkp = X.TwoBody.GetTBME_J(J2, J2, b, q, k, p);
                      double ybpkq = Y.TwoBody.GetTBME_J(J2, J2, b, p, k, q);
                      double ybqkp = Y.TwoBody.GetTBME_J(J2, J2, b, q, k, p);
                      zpq += 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (ok.j2 + 1) * ok.occ * xijab * yakij * xbpkq;
                      zpq += 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (ok.j2 + 1) * ok.occ * xijab * yakij * xbqkp * hX * hX * hY;

                      /// Also include the XXY and YXX contributions
                      zpq -= 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (ok.j2 + 1) * ok.occ * xijab * xakij * ybpkq;
                      zpq -= 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (ok.j2 + 1) * ok.occ * xijab * xakij * ybqkp * hX * hX * hY;
                    }
                    Z.OneBody(p, q) += zpq;
                  }
                }
              }
            }
          }
        }
      }
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void diagram_CIIa(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto a : Z.modelspace->particles)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (auto b : Z.modelspace->particles)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        for (auto i : Z.modelspace->holes)
        {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          for (auto j : Z.modelspace->holes)
          {
            Orbit &oj = Z.modelspace->GetOrbit(j);
            for (auto c : Z.modelspace->particles)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              int J1min = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
              int J1max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                double xijab = X.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                for (auto p : Z.modelspace->all_orbits)
                {
                  Orbit &op = Z.modelspace->GetOrbit(p);
                  for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
                  {
                    Orbit &oq = Z.modelspace->GetOrbit(q);
                    double zpq = 0;

                    double yabcq = Y.TwoBody.GetTBME_J(J1, J1, a, b, c, q);
                    double yabcp = Y.TwoBody.GetTBME_J(J1, J1, a, b, c, p);
                    double xcpij = X.TwoBody.GetTBME_J(J1, J1, c, p, i, j);
                    double xcqij = X.TwoBody.GetTBME_J(J1, J1, c, q, i, j);

                    double xabcq = X.TwoBody.GetTBME_J(J1, J1, a, b, c, q);
                    double xabcp = X.TwoBody.GetTBME_J(J1, J1, a, b, c, p);
                    double ycpij = Y.TwoBody.GetTBME_J(J1, J1, c, p, i, j);
                    double ycqij = Y.TwoBody.GetTBME_J(J1, J1, c, q, i, j);
                    zpq += 0.25 * (2 * J1 + 1.) / (op.j2 + 1) * xijab * yabcq * xcpij;
                    zpq += 0.25 * (2 * J1 + 1.) / (op.j2 + 1) * xijab * yabcp * xcqij * hX * hX * hY;
                    /// Also include the XXY and YXX contributions
                    zpq -= 0.25 * (2 * J1 + 1.) / (op.j2 + 1) * xijab * xabcq * ycpij;
                    zpq -= 0.25 * (2 * J1 + 1.) / (op.j2 + 1) * xijab * xabcp * ycqij * hX * hX * hY;
                    Z.OneBody(p, q) += zpq;
                  }
                }
              }
            }
          }
        }
      }
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void diagram_CIIb(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    Z.modelspace->PreCalculateSixJ();
    size_t norbits = Z.modelspace->GetNumberOrbits();

    //   for ( auto p : Z.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic)
    for (size_t p = 0; p < norbits; p++)
    {
      Orbit &op = Z.modelspace->GetOrbit(p);
      for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
      {
        Orbit &oq = Z.modelspace->GetOrbit(q);
        double zpq = 0;

        for (auto a : Z.modelspace->particles)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          for (auto b : Z.modelspace->particles)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            for (auto i : Z.modelspace->holes)
            {
              Orbit &oi = Z.modelspace->GetOrbit(i);
              for (auto j : Z.modelspace->holes)
              {
                Orbit &oj = Z.modelspace->GetOrbit(j);
                for (auto c : Z.modelspace->particles)
                {
                  Orbit &oc = Z.modelspace->GetOrbit(c);

                  int J4min = std::max(std::abs(op.j2 - oc.j2), std::abs(oj.j2 - ob.j2)) / 2;
                  int J4max = std::min(op.j2 + oc.j2, oj.j2 + ob.j2) / 2;
                  for (int J4 = J4min; J4 <= J4max; J4++)
                  {

                    double xijab = 0;
                    int J1min = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
                    int J1max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
                    for (int J1 = J1min; J1 <= J1max; J1++)
                    {
                      // double xijab = X.TwoBody.GetTBME_J(J1,J1,i,j,a,b);
                      double sixj1 = AngMom::phase(J1) * (2 * J1 + 1) * Z.modelspace->GetSixJ(oa.j2 * 0.5, oi.j2 * 0.5, J4, oj.j2 * 0.5, ob.j2 * 0.5, J1);
                      xijab += sixj1 * X.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                    }

                    double ybpjc = 0;
                    double ybqjc = 0;
                    double xbpjc = 0;
                    double xbqjc = 0;
                    int J2min = std::max(std::abs(ob.j2 - op.j2), std::abs(oj.j2 - oc.j2)) / 2;
                    int J2max = std::min(ob.j2 + op.j2, oj.j2 + oc.j2) / 2;
                    for (int J2 = J2min; J2 <= J2max; J2++)
                    {
                      //                      double ybpjc = Y.TwoBody.GetTBME_J(J2,J2,b,p,j,c);
                      //                      double ybqjc = Y.TwoBody.GetTBME_J(J2,J2,b,q,j,c);
                      //                      double xbpjc = X.TwoBody.GetTBME_J(J2,J2,b,p,j,c);
                      //                      double xbqjc = X.TwoBody.GetTBME_J(J2,J2,b,q,j,c);
                      double sixj2 = AngMom::phase(J2) * (2 * J2 + 1) * Z.modelspace->GetSixJ(op.j2 * 0.5, oc.j2 * 0.5, J4, oj.j2 * 0.5, ob.j2 * 0.5, J2);
                      ybpjc += sixj2 * Y.TwoBody.GetTBME_J(J2, J2, b, p, j, c);
                      ybqjc += sixj2 * Y.TwoBody.GetTBME_J(J2, J2, b, q, j, c);
                      xbpjc += sixj2 * X.TwoBody.GetTBME_J(J2, J2, b, p, j, c);
                      xbqjc += sixj2 * X.TwoBody.GetTBME_J(J2, J2, b, q, j, c);
                    }

                    double xaciq = 0;
                    double xacip = 0;
                    double yaciq = 0;
                    double yacip = 0;
                    int J3min = std::max(std::abs(oa.j2 - oc.j2), std::abs(oi.j2 - oq.j2)) / 2;
                    int J3max = std::min(oa.j2 + oc.j2, oi.j2 + oq.j2) / 2;
                    for (int J3 = J3min; J3 <= J3max; J3++)
                    {
                      // double xaciq = X.TwoBody.GetTBME_J(J3,J3,a,c,i,q);
                      // double xacip = X.TwoBody.GetTBME_J(J3,J3,a,c,i,p);
                      // double yaciq = Y.TwoBody.GetTBME_J(J3,J3,a,c,i,q);
                      // double yacip = Y.TwoBody.GetTBME_J(J3,J3,a,c,i,p);
                      double sixj3 = AngMom::phase(J3) * (2 * J3 + 1) * Z.modelspace->GetSixJ(oa.j2 * 0.5, oi.j2 * 0.5, J4, op.j2 * 0.5, oc.j2 * 0.5, J3);
                      xaciq += sixj3 * X.TwoBody.GetTBME_J(J3, J3, a, c, i, q);
                      xacip += sixj3 * X.TwoBody.GetTBME_J(J3, J3, a, c, i, p);
                      yaciq += sixj3 * Y.TwoBody.GetTBME_J(J3, J3, a, c, i, q);
                      yacip += sixj3 * Y.TwoBody.GetTBME_J(J3, J3, a, c, i, p);
                    }

                    //    int J4min = std::max( std::abs( op.j2-oc.j2), std::abs( oj.j2-ob.j2) )/2;
                    //    int J4max = std::min(  op.j2+oc.j2,  oj.j2+ob.j2 )/2;
                    //    double sixj_prod = 0;
                    //    for (int J4=J4min; J4<=J4max; J4++)
                    //    {
                    //       double sixj1 = Z.modelspace->GetSixJ( oa.j2*0.5, oi.j2*0.5, J4,   oj.j2*0.5, ob.j2*0.5, J1);
                    //       double sixj2 = Z.modelspace->GetSixJ( op.j2*0.5, oc.j2*0.5, J4,   oj.j2*0.5, ob.j2*0.5, J2);
                    //       double sixj3 = Z.modelspace->GetSixJ( oa.j2*0.5, oi.j2*0.5, J4,   op.j2*0.5, oc.j2*0.5, J3);
                    //       sixj_prod += (2*J4+1) * AngMom::phase(J1+J2+J3+(ob.j2+oj.j2)/2) * sixj1*sixj2*sixj3;
                    //    }

                    //                           zpq -= (2*J1+1.)*(2*J2+1)*(2*J3+1)/(op.j2+1) * sixj_prod *  xijab * ybpjc * xaciq;
                    //                           zpq -= (2*J1+1.)*(2*J2+1)*(2*J3+1)/(op.j2+1) * sixj_prod *  xijab * ybqjc * xacip * hX*hX*hY;
                    zpq -= (2 * J4 + 1.) / (op.j2 + 1) * AngMom::phase((ob.j2 + oj.j2) / 2) * xijab * ybpjc * xaciq;
                    zpq -= (2 * J4 + 1.) / (op.j2 + 1) * AngMom::phase((ob.j2 + oj.j2) / 2) * xijab * ybqjc * xacip * hX * hX * hY;
                    /// Also include the XXY and YXX contributions
                    zpq += (2 * J4 + 1.) / (op.j2 + 1) * AngMom::phase((ob.j2 + oj.j2) / 2) * xijab * xbpjc * yaciq;
                    zpq += (2 * J4 + 1.) / (op.j2 + 1) * AngMom::phase((ob.j2 + oj.j2) / 2) * xijab * xbqjc * yacip * hX * hX * hY;
                    //}
                    //}
                  } // for J4
                }
              } // for j
            } // for i
          } // for b
        } // for a
        Z.OneBody(p, q) += zpq;
      } // for q
    } // for p
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void diagram_CIIc(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto a : Z.modelspace->particles)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (auto b : Z.modelspace->particles)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        for (auto i : Z.modelspace->holes)
        {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          for (auto j : Z.modelspace->holes)
          {
            Orbit &oj = Z.modelspace->GetOrbit(j);
            for (auto k : Z.modelspace->holes)
            {
              Orbit &ok = Z.modelspace->GetOrbit(k);
              int J1min = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
              int J1max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                double xijab = X.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                for (auto p : Z.modelspace->all_orbits)
                {
                  Orbit &op = Z.modelspace->GetOrbit(p);
                  for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
                  {
                    Orbit &oq = Z.modelspace->GetOrbit(q);
                    double zpq = 0;
                    int J2min = std::max(std::abs(ob.j2 - ok.j2), std::abs(oj.j2 - oq.j2)) / 2;
                    int J2max = std::min(ob.j2 + ok.j2, oj.j2 + oq.j2) / 2;
                    for (int J2 = J2min; J2 <= J2max; J2++)
                    {
                      double ybkjq = Y.TwoBody.GetTBME_J(J2, J2, b, k, j, q);
                      double ybkjp = Y.TwoBody.GetTBME_J(J2, J2, b, k, j, p);
                      double xbkjq = X.TwoBody.GetTBME_J(J2, J2, b, k, j, q);
                      double xbkjp = X.TwoBody.GetTBME_J(J2, J2, b, k, j, p);

                      int J3min = std::max(std::abs(oa.j2 - op.j2), std::abs(oi.j2 - ok.j2)) / 2;
                      int J3max = std::min(oa.j2 + op.j2, oi.j2 + ok.j2) / 2;
                      for (int J3 = J3min; J3 <= J3max; J3++)
                      {
                        double xapik = X.TwoBody.GetTBME_J(J3, J3, a, p, i, k);
                        double xaqik = X.TwoBody.GetTBME_J(J3, J3, a, q, i, k);
                        double yapik = Y.TwoBody.GetTBME_J(J3, J3, a, p, i, k);
                        double yaqik = Y.TwoBody.GetTBME_J(J3, J3, a, q, i, k);

                        int J4min = std::max(std::abs(op.j2 - ok.j2), std::abs(oj.j2 - ob.j2)) / 2;
                        int J4max = std::min(op.j2 + ok.j2, oj.j2 + ob.j2) / 2;
                        double sixj_prod = 0;
                        for (int J4 = J4min; J4 <= J4max; J4++)
                        {
                          double sixj1 = Z.modelspace->GetSixJ(oa.j2 * 0.5, oi.j2 * 0.5, J4, oj.j2 * 0.5, ob.j2 * 0.5, J1);
                          double sixj2 = Z.modelspace->GetSixJ(ok.j2 * 0.5, op.j2 * 0.5, J4, oj.j2 * 0.5, ob.j2 * 0.5, J2);
                          double sixj3 = Z.modelspace->GetSixJ(oa.j2 * 0.5, oi.j2 * 0.5, J4, ok.j2 * 0.5, op.j2 * 0.5, J3);
                          sixj_prod += (2 * J4 + 1) * AngMom::phase(J1 + J2 + J3 + (ob.j2 + oj.j2) / 2) * sixj1 * sixj2 * sixj3;
                        }

                        zpq += (2 * J1 + 1.) * (2 * J2 + 1) * (2 * J3 + 1) / (op.j2 + 1) * sixj_prod * xijab * ybkjq * xapik;
                        zpq += (2 * J1 + 1.) * (2 * J2 + 1) * (2 * J3 + 1) / (op.j2 + 1) * sixj_prod * xijab * ybkjp * xaqik * hX * hX * hY;
                        /// Also include the XXY and YXX contributions
                        zpq -= (2 * J1 + 1.) * (2 * J2 + 1) * (2 * J3 + 1) / (op.j2 + 1) * sixj_prod * xijab * xbkjq * yapik;
                        zpq -= (2 * J1 + 1.) * (2 * J2 + 1) * (2 * J3 + 1) / (op.j2 + 1) * sixj_prod * xijab * xbkjp * yaqik * hX * hX * hY;
                      }
                    }
                    Z.OneBody(p, q) += zpq;
                  }
                }
              }
            }
          }
        }
      }
    }

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void diagram_CIId(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto a : Z.modelspace->particles)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (auto b : Z.modelspace->particles)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        for (auto i : Z.modelspace->holes)
        {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          for (auto j : Z.modelspace->holes)
          {
            Orbit &oj = Z.modelspace->GetOrbit(j);
            for (auto k : Z.modelspace->holes)
            {
              Orbit &ok = Z.modelspace->GetOrbit(k);
              int J1min = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
              int J1max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                double xijab = X.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                for (auto p : Z.modelspace->all_orbits)
                {
                  Orbit &op = Z.modelspace->GetOrbit(p);
                  for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
                  {
                    Orbit &oq = Z.modelspace->GetOrbit(q);
                    double zpq = 0;

                    double ykpij = Y.TwoBody.GetTBME_J(J1, J1, k, p, i, j);
                    double ykqij = Y.TwoBody.GetTBME_J(J1, J1, k, q, i, j);
                    double xabkq = X.TwoBody.GetTBME_J(J1, J1, a, b, k, q);
                    double xabkp = X.TwoBody.GetTBME_J(J1, J1, a, b, k, p);
                    double xkpij = X.TwoBody.GetTBME_J(J1, J1, k, p, i, j);
                    double xkqij = X.TwoBody.GetTBME_J(J1, J1, k, q, i, j);
                    double yabkq = Y.TwoBody.GetTBME_J(J1, J1, a, b, k, q);
                    double yabkp = Y.TwoBody.GetTBME_J(J1, J1, a, b, k, p);
                    zpq -= 0.25 * (2 * J1 + 1.) / (op.j2 + 1) * xijab * ykpij * xabkq;
                    zpq -= 0.25 * (2 * J1 + 1.) / (op.j2 + 1) * xijab * ykqij * xabkp * hX * hX * hY;
                    /// Also include the XXY and YXX contributions
                    zpq += 0.25 * (2 * J1 + 1.) / (op.j2 + 1) * xijab * xkpij * yabkq;
                    zpq += 0.25 * (2 * J1 + 1.) / (op.j2 + 1) * xijab * xkqij * yabkp * hX * hX * hY;
                    Z.OneBody(p, q) += zpq;
                  }
                }
              }
            }
          }
        }
      }
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void diagram_CIIIa(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto a : Z.modelspace->particles)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (auto b : Z.modelspace->particles)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        for (auto i : Z.modelspace->holes)
        {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          for (auto j : Z.modelspace->holes)
          {
            Orbit &oj = Z.modelspace->GetOrbit(j);
            for (auto c : X.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}))
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              if ((1 - oc.occ) < 1e-4)
                continue;
              int J1min = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
              int J1max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                double xijab = X.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                double yijab = Y.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                double xcbij = X.TwoBody.GetTBME_J(J1, J1, c, b, i, j);
                double ycbij = Y.TwoBody.GetTBME_J(J1, J1, c, b, i, j);
                for (auto p : Z.modelspace->all_orbits)
                {
                  Orbit &op = Z.modelspace->GetOrbit(p);
                  for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
                  {
                    Orbit &oq = Z.modelspace->GetOrbit(q);
                    double zpq = 0;
                    int J2min = std::max(std::abs(op.j2 - oa.j2), std::abs(oq.j2 - oc.j2)) / 2;
                    int J2max = std::min(op.j2 + oa.j2, oq.j2 + oc.j2) / 2;
                    for (int J2 = J2min; J2 <= J2max; J2++)
                    {
                      // The "time reversed" diagram is identical, so we just get a factor of 2
                      double ypaqc = Y.TwoBody.GetTBME_J(J2, J2, p, a, q, c);
                      double xpaqc = X.TwoBody.GetTBME_J(J2, J2, p, a, q, c);
                      zpq -= 2 * 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oa.j2 + 1) * (1 - oc.occ) * xijab * ypaqc * xcbij;
                      /// Also include the XXY and YXX contributions
                      zpq += 1 * 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oa.j2 + 1) * (1 - oc.occ) * xijab * xpaqc * ycbij;
                      zpq += 1 * 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oa.j2 + 1) * (1 - oc.occ) * yijab * xpaqc * xcbij;
                    }
                    Z.OneBody(p, q) += zpq;
                  }
                }
              }
            }
          }
        }
      }
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void diagram_CIIIb(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto a : Z.modelspace->particles)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      for (auto b : Z.modelspace->particles)
      {
        Orbit &ob = Z.modelspace->GetOrbit(b);
        for (auto i : Z.modelspace->holes)
        {
          Orbit &oi = Z.modelspace->GetOrbit(i);
          for (auto j : Z.modelspace->holes)
          {
            Orbit &oj = Z.modelspace->GetOrbit(j);
            for (auto k : X.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
            {
              Orbit &ok = Z.modelspace->GetOrbit(k);
              if (ok.occ < 1e-4)
                continue;
              int J1min = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
              int J1max = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
              for (int J1 = J1min; J1 <= J1max; J1++)
              {
                double xijab = X.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                double yijab = Y.TwoBody.GetTBME_J(J1, J1, i, j, a, b);
                double xabkj = X.TwoBody.GetTBME_J(J1, J1, a, b, k, j);
                double yabkj = Y.TwoBody.GetTBME_J(J1, J1, a, b, k, j);
                for (auto p : Z.modelspace->all_orbits)
                {
                  Orbit &op = Z.modelspace->GetOrbit(p);
                  for (auto q : Z.OneBodyChannels.at({op.l, op.j2, op.tz2}))
                  {
                    Orbit &oq = Z.modelspace->GetOrbit(q);
                    double zpq = 0;
                    int J2min = std::max(std::abs(op.j2 - ok.j2), std::abs(oq.j2 - oi.j2)) / 2;
                    int J2max = std::min(op.j2 + ok.j2, oq.j2 + oi.j2) / 2;
                    for (int J2 = J2min; J2 <= J2max; J2++)
                    {
                      double ypkqi = Y.TwoBody.GetTBME_J(J2, J2, p, k, q, i);
                      double xpkqi = X.TwoBody.GetTBME_J(J2, J2, p, k, q, i);
                      zpq += 2 * 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oi.j2 + 1) * ok.occ * xijab * ypkqi * xabkj;
                      /// Also include the XXY and YXX contributions
                      zpq -= 1 * 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oi.j2 + 1) * ok.occ * xijab * xpkqi * yabkj;
                      zpq -= 1 * 0.5 * (2 * J1 + 1.) * (2 * J2 + 1) / (op.j2 + 1) / (oi.j2 + 1) * ok.occ * yijab * xpkqi * xabkj;
                    }
                    Z.OneBody(p, q) += zpq;
                  }
                }
              }
            }
          }
        }
      }
    }
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  // Initially, write it as the straightforward N^8 sum.
  // Then worry about factorizing it.
  void diagram_DIa(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto itmat : Z.TwoBody.MatEl)
    {
      size_t ch_bra = itmat.first[0];
      size_t ch_ket = itmat.first[1];
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
        int phasepq = bra.Phase(J);
        //        for ( size_t iket=0; iket<nkets; iket++)
        for (size_t iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          int phasers = ket.Phase(J);
          double zpqrs = 0;
          for (auto a : Z.modelspace->particles)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            for (auto b : Z.modelspace->particles)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              for (auto c : Z.modelspace->particles)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                for (auto i : Z.modelspace->holes)
                {
                  Orbit &oi = Z.modelspace->GetOrbit(i);
                  double xcqrs = X.TwoBody.GetTBME_J(J, J, c, q, r, s);
                  double xpcrs = X.TwoBody.GetTBME_J(J, J, p, c, r, s);
                  double xpqrc = X.TwoBody.GetTBME_J(J, J, p, q, r, c);
                  double xpqcs = X.TwoBody.GetTBME_J(J, J, p, q, c, s);
                  int dcp = oc.j2 == bra.op->j2;
                  int dcq = oc.j2 == bra.oq->j2;
                  int dcr = oc.j2 == ket.op->j2;
                  int dcs = oc.j2 == ket.oq->j2;
                  //                     int J1min = std::max( std::abs( oa.j2-ob.j2), std::abs(oi.j2-oc.j2) )/2;
                  //                     int J1max = std::min(  oa.j2+ob.j2, oi.j2+oc.j2 )/2;
                  int J1min = std::abs(oa.j2 - ob.j2) / 2;
                  int J1max = (oa.j2 + ob.j2) / 2;
                  for (int J1 = J1min; J1 <= J1max; J1++)
                  {
                    double yabic = Y.TwoBody.GetTBME_J(J1, J1, a, b, i, c);
                    double xipab = X.TwoBody.GetTBME_J(J1, J1, i, p, a, b);
                    double xabir = X.TwoBody.GetTBME_J(J1, J1, a, b, i, r);
                    double xabis = X.TwoBody.GetTBME_J(J1, J1, a, b, i, s);
                    double xiqab = X.TwoBody.GetTBME_J(J1, J1, i, p, a, b);
                    // zpqrs += 0.5*(2*J1+1)/(bra.op->j2+1) * ( xipab * yabic * xcqrs * dcp  +phasepq* xiqab * yabic * xcprs * dcq);  // Includes the (1-Ppq)
                    // zpqrs += 0.5*(2*J1+1)/(bra.op->j2+1) * ( xpqcs * yabic * xabir * dcr  +phasers* xpqcr * yabic * xabis * dcs) * hY; // DIa*   includes (1-Prs)
                    zpqrs += 0.5 * (2 * J1 + 1) / (bra.op->j2 + 1) * (xipab * yabic * xcqrs * dcp + xiqab * yabic * xpcrs * dcq);      // Includes the (1-Ppq)
                    zpqrs += 0.5 * (2 * J1 + 1) / (bra.op->j2 + 1) * (xpqcs * yabic * xabir * dcr + xpqrc * yabic * xabis * dcs) * hY; // DIa*   includes (1-Prs)
                  }
                } // for j
              } // for i
            } // for b
          } // for a
          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      } // for ibra

    } // for itmat

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  /*
  void diagram_DIb( const Operator& X, const Operator& Y, Operator& Z )
  {
    double t_start = omp_get_wtime();
     int hX = X.IsHermitian() ? +1 : -1;
     int hY = Y.IsHermitian() ? +1 : -1;
    arma::mat CHI_XX = Y.OneBody * 0;
    arma::mat CHI_XY = Y.OneBody * 0;
    size_t norb = Z.modelspace->GetNumberOrbits();

    #pragma omp parallel for schedule(dynamic)
    for ( size_t p=0; p<norb; p++)
    {
       Orbit& op = Z.modelspace->GetOrbit(p);
       for ( auto q : X.OneBodyChannels.at({op.l,op.j2,op.tz2}) )
       {
           Orbit& oq = Z.modelspace->GetOrbit(q);

           double chi_pq = 0;
           double chiY_pq = 0;

           for ( auto a : Z.modelspace->particles)
            {
              Orbit& oa = Z.modelspace->GetOrbit(a);

              for ( auto i : Z.modelspace->holes)
              {
                Orbit& oi = Z.modelspace->GetOrbit(i);


                for ( auto j : Z.modelspace->holes)
                {
                    Orbit& oj = Z.modelspace->GetOrbit(j);

                    int J1min = std::max( std::abs( oa.j2-oq.j2), std::abs(oi.j2-oj.j2) )/2;
                    int J1max = std::min(  oa.j2+oq.j2, oi.j2+oj.j2 )/2;

                    for ( int J1=J1min; J1<=J1max; J1++)
                    {

                       double xijaq = X.TwoBody.GetTBME_J(J1,J1,i,j,a,q);
                       double xapij = X.TwoBody.GetTBME_J(J1,J1,a,p,i,j);
                       double yapij = Y.TwoBody.GetTBME_J(J1,J1,a,p,i,j);

                       //chi_pq  += 0.5*(2*J1+1)/(oq.j2+1) * (1-oq.occ) * xijap * xaqij ;
                      // chiY_pq += 0.5*(2*J1+1)/(oq.j2+1) * (1-oq.occ) * xijap * yaqij ;
                       chi_pq  += 0.5*(2*J1+1)/(oq.j2+1) * (1) * xijaq * xapij ;
                       chiY_pq += 0.5*(2*J1+1)/(oq.j2+1) * (1) * xijaq * yapij ;

                       ///   chi_rk += 0.5*(2*J1+1)/(bra.op->j2+1) * ( xijra * ykaij  *dkr  );  // Includes the (1-Prs)
                          // chi_pq +=                             * ( xijpa & yqaij * dpq )  where ij holes, a particle.

                    }
                }//for j

                //for ( auto b : Z.modelspace->particles)
                //{
                //    Orbit& ob = Z.modelspace->GetOrbit(b);

                //    int J1min = std::max( std::abs( oa.j2-ob.j2), std::abs(oi.j2-oq.j2) )/2;
                //    int J1max = std::min(  oa.j2+ob.j2, oi.j2+oq.j2 )/2;

                //    for ( int J1=J1min; J1<=J1max; J1++)
                //    {

                //       double xipab = X.TwoBody.GetTBME_J(J1,J1,i,p,a,b);
                //       double xabiq = X.TwoBody.GetTBME_J(J1,J1,a,b,i,q);
                //       double yabiq = Y.TwoBody.GetTBME_J(J1,J1,a,b,i,q);

                //       chi_pq  += 0.5*(2*J1+1)/(oq.j2+1) * 1 * xipab * xabiq ;
                //       chiY_pq += 0.5*(2*J1+1)/(oq.j2+1) * 1 * xipab * yabiq ;
                //    }
                //}//for b

              }//for i
            }//for a
            CHI_XX(p,q) = chi_pq;
            CHI_XY(p,q) = chiY_pq;
            }//for q
        }// for p








     for ( auto itmat : Z.TwoBody.MatEl )
     {
        size_t ch_bra = itmat.first[0];
        size_t ch_ket = itmat.first[1];
        TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int J = tbc_bra.J;
        size_t nbras = tbc_bra.GetNumberKets();
        size_t nkets = tbc_ket.GetNumberKets();
        for ( size_t ibra=0; ibra<nbras; ibra++)
        {
          Ket& bra = tbc_bra.GetKet(ibra);
          size_t p = bra.p;
          size_t q = bra.q;
  //        for ( size_t iket=0; iket<nkets; iket++)
          for ( size_t iket=ibra; iket<nkets; iket++)
          {
             Ket& ket = tbc_ket.GetKet(iket);
             size_t r = ket.p;
             size_t s = ket.q;
             double zpqrs = 0;
             for ( auto k : Z.modelspace->holes)
             {
                 Orbit& ok = Z.modelspace->GetOrbit(k);
                 double xpqks = X.TwoBody.GetTBME_J(J,J,p,q,k,s);
                 double xpqrk = X.TwoBody.GetTBME_J(J,J,p,q,r,k);
                 double xpkrs = X.TwoBody.GetTBME_J(J,J,p,k,r,s);
                 double xkqrs = X.TwoBody.GetTBME_J(J,J,k,q,r,s);

  //             zpqrs += xpqks * chi_rk + xpqrk * chi_sk + xkqrs * chi_pk + xpkrs * chi_qk;
               zpqrs += xpqks * CHI_XY(k,r) + xpqrk * CHI_XY(k,s) + CHI_XY(p,k) * xkqrs  +  CHI_XY(q,k) * xpkrs ;
             }//for k
             // normalize
             if (p==q ) zpqrs /= PhysConst::SQRT2;
             if (r==s ) zpqrs /= PhysConst::SQRT2;
             Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
          }// for iket
        }// for ibra

     }// for itmat

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }
  */

  void diagram_DIb(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto itmat : Z.TwoBody.MatEl)
    {
      size_t ch_bra = itmat.first[0];
      size_t ch_ket = itmat.first[1];
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
        //        for ( size_t iket=0; iket<nkets; iket++)
        for (size_t iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          double zpqrs = 0;
          for (auto k : Z.modelspace->holes)
          {
            Orbit &ok = Z.modelspace->GetOrbit(k);
            double xpqks = X.TwoBody.GetTBME_J(J, J, p, q, k, s);
            double xpqrk = X.TwoBody.GetTBME_J(J, J, p, q, r, k);
            double xpkrs = X.TwoBody.GetTBME_J(J, J, p, k, r, s);
            double xkqrs = X.TwoBody.GetTBME_J(J, J, k, q, r, s);
            int dkp = ok.j2 == bra.op->j2;
            int dkq = ok.j2 == bra.oq->j2;
            int dkr = ok.j2 == ket.op->j2;
            int dks = ok.j2 == ket.oq->j2;

            double chi_rk = 0;
            double chi_sk = 0;
            double chi_pk = 0;
            double chi_qk = 0;

            for (auto a : Z.modelspace->particles)
            {
              Orbit &oa = Z.modelspace->GetOrbit(a);
              for (auto i : Z.modelspace->holes)
              {
                Orbit &oi = Z.modelspace->GetOrbit(i);
                for (auto j : Z.modelspace->holes)
                {
                  Orbit &oj = Z.modelspace->GetOrbit(j);

                  //                     int J1min = std::max( std::abs( oa.j2-ob.j2), std::abs(oi.j2-oc.j2) )/2;
                  //                     int J1max = std::min(  oa.j2+ob.j2, oi.j2+oc.j2 )/2;
                  int J1min = std::abs(oi.j2 - oj.j2) / 2;
                  int J1max = (oi.j2 + oj.j2) / 2;
                  for (int J1 = J1min; J1 <= J1max; J1++)
                  {
                    double ykaij = Y.TwoBody.GetTBME_J(J1, J1, k, a, i, j);
                    double yijka = Y.TwoBody.GetTBME_J(J1, J1, i, j, k, a);
                    double xijra = X.TwoBody.GetTBME_J(J1, J1, i, j, r, a);
                    double xijsa = X.TwoBody.GetTBME_J(J1, J1, i, j, s, a);
                    double xpaij = X.TwoBody.GetTBME_J(J1, J1, p, a, i, j);
                    double xqaij = X.TwoBody.GetTBME_J(J1, J1, q, a, i, j);
                    // zpqrs += 0.5*(2*J1+1)/(bra.op->j2+1) * ( xijra * ykaij * xpqks *dkr  + xijsa * ykaij * xpqrk *dks);  // Includes the (1-Prs)
                    // zpqrs += 0.5*(2*J1+1)/(bra.op->j2+1) * ( xkqrs * ykaij * xpaij *dkp  + xpkrs * ykaij * xqaij *dkq) * hY; // DIa*   includes (1-Ppq)

                    chi_rk += 0.5 * (2 * J1 + 1) / (ket.op->j2 + 1) * (xijra * ykaij * dkr); // Includes the (1-Prs)
                    chi_sk += 0.5 * (2 * J1 + 1) / (ket.op->j2 + 1) * (xijsa * ykaij * dks); // Includes the (1-Prs)
                    chi_pk += 0.5 * (2 * J1 + 1) / (ket.op->j2 + 1) * (yijka * xpaij * dkp); // DIa*   includes (1-Ppq)
                    chi_qk += 0.5 * (2 * J1 + 1) / (ket.op->j2 + 1) * (yijka * xqaij * dkq); // DIa*   includes (1-Ppq)
                  }
                } // for j
              } // for i
            } // for a
            zpqrs += xpqks * chi_rk + xpqrk * chi_sk + xkqrs * chi_pk + xpkrs * chi_qk;
          } // for k
          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      } // for ibra

    } // for itmat

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void diagram_DIVa(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto itmat : Z.TwoBody.MatEl)
    {
      size_t ch_bra = itmat.first[0];
      size_t ch_ket = itmat.first[1];
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
        //       int phasepq = bra.Phase(J);
        //       for ( size_t iket=0; iket<nkets; iket++)
        for (size_t iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          //          int phasers = ket.Phase(J);
          double zpqrs = 0;
          for (auto a : Z.modelspace->particles)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            for (auto b : Z.modelspace->particles)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              for (auto i : Z.modelspace->holes)
              {
                Orbit &oi = Z.modelspace->GetOrbit(i);
                for (auto j : Z.modelspace->holes)
                {
                  Orbit &oj = Z.modelspace->GetOrbit(j);
                  double yjqrs = Y.TwoBody.GetTBME_J(J, J, j, q, r, s);
                  double ypjrs = Y.TwoBody.GetTBME_J(J, J, p, j, r, s);
                  double ypqjs = Y.TwoBody.GetTBME_J(J, J, p, q, j, s);
                  double ypqrj = Y.TwoBody.GetTBME_J(J, J, p, q, r, j);
                  double xjqrs = X.TwoBody.GetTBME_J(J, J, j, q, r, s);
                  double xpjrs = X.TwoBody.GetTBME_J(J, J, p, j, r, s);
                  double xpqjs = X.TwoBody.GetTBME_J(J, J, p, q, j, s);
                  double xpqrj = X.TwoBody.GetTBME_J(J, J, p, q, r, j);
                  int djp = oj.j2 == bra.op->j2;
                  int djq = oj.j2 == bra.oq->j2;
                  int djr = oj.j2 == ket.op->j2;
                  int djs = oj.j2 == ket.oq->j2;
                  //                    int J1min = std::max( std::abs( oa.j2-ob.j2), std::abs(oi.j2-oc.j2) )/2;
                  //                    int J1max = std::min(  oa.j2+ob.j2, oi.j2+oc.j2 )/2;
                  int J1min = std::max(std::abs(oa.j2 - ob.j2), std::abs(oi.j2 - oj.j2)) / 2;
                  int J1max = std::min(oa.j2 + ob.j2, oi.j2 + oj.j2) / 2;
                  //                    int J1min =  std::abs( oa.j2-ob.j2 )/2;
                  //                    int J1max = ( oa.j2+ob.j2 )/2;
                  for (int J1 = J1min; J1 <= J1max; J1++)
                  {
                    // if ( std::abs( oa.occ*ob.occ)<1e-3) continue;
                    double xipab = X.TwoBody.GetTBME_J(J1, J1, i, p, a, b);
                    double xiqab = X.TwoBody.GetTBME_J(J1, J1, i, q, a, b);
                    double xirab = X.TwoBody.GetTBME_J(J1, J1, i, r, a, b);
                    double xisab = X.TwoBody.GetTBME_J(J1, J1, i, s, a, b);
                    double xabij = X.TwoBody.GetTBME_J(J1, J1, a, b, i, j);
                    double yabij = Y.TwoBody.GetTBME_J(J1, J1, a, b, i, j);

                    zpqrs += 0.5 * (2 * J1 + 1) / (oj.j2 + 1) * (xipab * yjqrs * xabij * djp + xiqab * ypjrs * xabij * djq);           // Includes the (1-Ppq)
                    zpqrs += 0.5 * (2 * J1 + 1) / (oj.j2 + 1) * (xirab * ypqjs * xabij * djr + xisab * ypqrj * xabij * djs) * hX * hX; // DIVa*   includes (1-Prs)
                    // Need to include the XXY and YXX term
                    zpqrs -= 0.5 * (2 * J1 + 1) / (oj.j2 + 1) * (xipab * xjqrs * yabij * djp + xiqab * xpjrs * yabij * djq);           // Includes the (1-Ppq)
                    zpqrs -= 0.5 * (2 * J1 + 1) / (oj.j2 + 1) * (xirab * xpqjs * yabij * djr + xisab * xpqrj * yabij * djs) * hX * hY; // DIVa*   includes (1-Prs)
                  }
                } // for j
              } // for i
            } // for b
          } // for a
          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      } // for ibra

    } // for itmat

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  void diagram_DIVb(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;
    for (auto itmat : Z.TwoBody.MatEl)
    {
      size_t ch_bra = itmat.first[0];
      size_t ch_ket = itmat.first[1];
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
        //       for ( size_t iket=0; iket<nkets; iket++)
        for (size_t iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          double zpqrs = 0;
          for (auto a : Z.modelspace->particles)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            for (auto b : Z.modelspace->particles)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              for (auto i : Z.modelspace->holes)
              {
                Orbit &oi = Z.modelspace->GetOrbit(i);
                for (auto j : Z.modelspace->holes)
                {
                  Orbit &oj = Z.modelspace->GetOrbit(j);
                  double ypqbs = Y.TwoBody.GetTBME_J(J, J, p, q, b, s);
                  double ypqrb = Y.TwoBody.GetTBME_J(J, J, p, q, r, b);
                  double ybqrs = Y.TwoBody.GetTBME_J(J, J, b, q, r, s);
                  double ypbrs = Y.TwoBody.GetTBME_J(J, J, p, b, r, s);
                  double xpqbs = X.TwoBody.GetTBME_J(J, J, p, q, b, s);
                  double xpqrb = X.TwoBody.GetTBME_J(J, J, p, q, r, b);
                  double xbqrs = X.TwoBody.GetTBME_J(J, J, b, q, r, s);
                  double xpbrs = X.TwoBody.GetTBME_J(J, J, p, b, r, s);
                  int J1min = std::max(std::abs(oa.j2 - ob.j2), std::abs(oi.j2 - oj.j2)) / 2;
                  int J1max = std::min(oa.j2 + ob.j2, oi.j2 + oj.j2) / 2;

                  int dbp = ob.j2 == bra.op->j2;
                  int dbq = ob.j2 == bra.oq->j2;
                  int dbr = ob.j2 == ket.op->j2;
                  int dbs = ob.j2 == ket.oq->j2;

                  for (int J1 = J1min; J1 <= J1max; J1++)
                  {

                    double xijar = X.TwoBody.GetTBME_J(J1, J1, i, j, a, r);
                    double xijas = X.TwoBody.GetTBME_J(J1, J1, i, j, a, s);
                    double xijap = X.TwoBody.GetTBME_J(J1, J1, i, j, a, p);
                    double xijaq = X.TwoBody.GetTBME_J(J1, J1, i, j, a, q);
                    double xabij = X.TwoBody.GetTBME_J(J1, J1, a, b, i, j);
                    double yabij = Y.TwoBody.GetTBME_J(J1, J1, a, b, i, j);

                    // Here I suspect the * term just amounts to a factor of 2
                    zpqrs += 0.5 * (2 * J1 + 1) / (ob.j2 + 1) * (xijar * ypqbs * xabij * dbr + xijas * ypqrb * xabij * dbs);           // Includes the (1-Prs)
                    zpqrs += 0.5 * (2 * J1 + 1) / (ob.j2 + 1) * (xijap * ybqrs * xabij * dbp + xijaq * ypbrs * xabij * dbq) * hX * hX; // DIVb*   includes (1-Ppq)
                                                                                                                                       //                       // Need to include the XXY and YXX term
                    zpqrs -= 0.5 * (2 * J1 + 1) / (ob.j2 + 1) * (xijar * xpqbs * yabij * dbr + xijas * xpqrb * yabij * dbs);           // Includes the (1-Prs)
                    zpqrs -= 0.5 * (2 * J1 + 1) / (ob.j2 + 1) * (xijap * xbqrs * yabij * dbp + xijaq * xpbrs * yabij * dbq) * hX * hY; // DIVb*   includes (1-Ppq)
                  }
                } // for j
              } // for i
            } // for b
          } // for a
          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      } // for ibra

    } // for itmat

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  /// Now do the DIVa and DIVb together by constructing a one-body intermediate
  /// This reduces the scaling to N^5 + N^5.
  void diagram_DIVb_intermediate(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;

    arma::mat CHI_XX = Y.OneBody * 0;
    arma::mat CHI_XY = Y.OneBody * 0;
    //  size_t nch = Z.modelspace->GetNumberTwoBodyChannels();
    size_t norb = Z.modelspace->GetNumberOrbits();

    //  for ( auto p : Z.modelspace->all_orbits)
#pragma omp parallel for schedule(dynamic)
    for (size_t p = 0; p < norb; p++)
    {
      Orbit &op = Z.modelspace->GetOrbit(p);
      for (auto q : X.OneBodyChannels.at({op.l, op.j2, op.tz2}))
      {
        Orbit &oq = Z.modelspace->GetOrbit(q);

        double chi_pq = 0;
        double chiY_pq = 0;

        for (auto a : Z.modelspace->particles)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);

          for (auto i : Z.modelspace->holes)
          {
            Orbit &oi = Z.modelspace->GetOrbit(i);

            for (auto j : Z.modelspace->holes)
            {
              Orbit &oj = Z.modelspace->GetOrbit(j);

              int J1min = std::max(std::abs(oa.j2 - oq.j2), std::abs(oi.j2 - oj.j2)) / 2;
              int J1max = std::min(oa.j2 + oq.j2, oi.j2 + oj.j2) / 2;

              for (int J1 = J1min; J1 <= J1max; J1++)
              {

                double xijaq = X.TwoBody.GetTBME_J(J1, J1, i, j, a, q);
                double xapij = X.TwoBody.GetTBME_J(J1, J1, a, p, i, j);
                double yapij = Y.TwoBody.GetTBME_J(J1, J1, a, p, i, j);

                chi_pq += 0.5 * (2 * J1 + 1) / (oq.j2 + 1) * xapij * xijaq;
                chiY_pq += 0.5 * (2 * J1 + 1) / (oq.j2 + 1) * yapij * xijaq;
              }
            } // for j

            for (auto b : Z.modelspace->particles)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);

              int J1min = std::max({std::abs(oa.j2 - ob.j2), std::abs(oi.j2 - oq.j2), std::abs(oi.j2 - op.j2)}) / 2;
              int J1max = std::min({oa.j2 + ob.j2, oi.j2 + oq.j2, oi.j2 + op.j2}) / 2;

              for (int J1 = J1min; J1 <= J1max; J1++)
              {

                double xipab = X.TwoBody.GetTBME_J(J1, J1, i, p, a, b);
                double xabiq = X.TwoBody.GetTBME_J(J1, J1, a, b, i, q);
                double yabiq = Y.TwoBody.GetTBME_J(J1, J1, a, b, i, q);
                double yipab = Y.TwoBody.GetTBME_J(J1, J1, i, p, a, b);

                chi_pq += 0.5 * (2 * J1 + 1) / (oq.j2 + 1) * xipab * xabiq;
                //                chiY_pq += 0.5 * (2 * J1 + 1) / (oq.j2 + 1) * xipab * yabiq;
                chiY_pq += 0.5 * (2 * J1 + 1) / (oq.j2 + 1) * yipab * xabiq; // JUST TRYING THIS WITHOUT CONFIRMING...
              }
            } // for b

          } // for i
        } // for a
        CHI_XX(p, q) = chi_pq;
        CHI_XY(p, q) = chiY_pq;
      } // for q
    } // for p

    std::vector<size_t> ch_bra_list, ch_ket_list;
    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    size_t nch = ch_bra_list.size();

    //   int nch = Z.modelspace->GetNumberTwoBodyChannels();
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
        int phasepq = bra.Phase(J);
        //       for ( size_t iket=0; iket<nkets; iket++)
        for (size_t iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          int phasers = ket.Phase(J);
          double zpqrs = 0;

          for (auto b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);

            zpqrs += CHI_XX(p, b) * Y.TwoBody.GetTBME_J(J, J, b, q, r, s) + CHI_XX(q, b) * Y.TwoBody.GetTBME_J(J, J, p, b, r, s);
            zpqrs += Y.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_XX(b, r) + Y.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_XX(b, s);
            ////            zpqrs -= CHI_XY(p, b) * X.TwoBody.GetTBME_J(J, J, b, q, r, s) + CHI_XY(q, b) * X.TwoBody.GetTBME_J(J, J, p, b, r, s);
            zpqrs += CHI_XY(b, p) * X.TwoBody.GetTBME_J(J, J, b, q, r, s) + CHI_XY(b, q) * X.TwoBody.GetTBME_J(J, J, p, b, r, s); // tricky minus sign
            zpqrs -= X.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_XY(b, r) + X.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_XY(b, s);

            // zpqrs += CHI_XY(b, p) * X.TwoBody.GetTBME_J(J, J, b, q, r, s) + CHI_XY(b, q) * X.TwoBody.GetTBME_J(J, J, p, b, r, s);
            // zpqrs -= X.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_XY(b, r) + X.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_XY(b, s);

          } // for a

          //          for ( auto j : Z.modelspace->holes)
          //          {
          //              Orbit& oj = Z.modelspace->GetOrbit(j);
          //
          //
          //              zpqrs += Y.TwoBody.GetTBME_J(J,J,j,q,r,s)* CHI_XX(p,j) +  Y.TwoBody.GetTBME_J(J,J,p,j,r,s)* CHI_XX(q,j);
          //              zpqrs += Y.TwoBody.GetTBME_J(J,J,p,q,j,s)* CHI_XX(r,j) +  Y.TwoBody.GetTBME_J(J,J,p,q,r,j)* CHI_XX(s,j);
          //              zpqrs -= X.TwoBody.GetTBME_J(J,J,j,q,r,s)* CHI_XY(p,j) +  X.TwoBody.GetTBME_J(J,J,p,j,r,s)* CHI_XY(q,j);
          //              zpqrs -= X.TwoBody.GetTBME_J(J,J,p,q,j,s)* CHI_XY(r,j) +  X.TwoBody.GetTBME_J(J,J,p,q,r,j)* CHI_XY(s,j);
          //
          //          }

          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      } // for ibra

    } // for itmat

    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  /*
  void diagram_DIVb_intermediate( const Operator& X, const Operator& Y, Operator& Z )
  {
    double t_start = omp_get_wtime();
    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;


    arma::mat CHI_XX = Y.OneBody * 0;
    arma::mat CHI_XY = Y.OneBody * 0;
  //  size_t nch = Z.modelspace->GetNumberTwoBodyChannels();
    size_t norb = Z.modelspace->GetNumberOrbits();

  //  for ( auto p : Z.modelspace->all_orbits)
    #pragma omp parallel for schedule(dynamic)
    for ( size_t p=0; p<norb; p++)
    {
       Orbit& op = Z.modelspace->GetOrbit(p);
       for ( auto q : X.OneBodyChannels.at({op.l,op.j2,op.tz2}) )
       {
           Orbit& oq = Z.modelspace->GetOrbit(q);

           double chi_pq = 0;
           double chiY_pq = 0;

           for ( auto a : Z.modelspace->particles)
            {
              Orbit& oa = Z.modelspace->GetOrbit(a);

              for ( auto i : Z.modelspace->holes)
              {
                Orbit& oi = Z.modelspace->GetOrbit(i);


                for ( auto j : Z.modelspace->holes)
                {
                    Orbit& oj = Z.modelspace->GetOrbit(j);

                    int J1min = std::max( std::abs( oa.j2-oq.j2), std::abs(oi.j2-oj.j2) )/2;
                    int J1max = std::min(  oa.j2+oq.j2, oi.j2+oj.j2 )/2;

                    for ( int J1=J1min; J1<=J1max; J1++)
                    {

                       double xijap = X.TwoBody.GetTBME_J(J1,J1,i,j,a,p);
                       double xaqij = X.TwoBody.GetTBME_J(J1,J1,a,q,i,j);
                       double yaqij = Y.TwoBody.GetTBME_J(J1,J1,a,q,i,j);

                       chi_pq  += 0.5*(2*J1+1)/(oq.j2+1) * (1-oq.occ) * xijap * xaqij ;
                       chiY_pq += 0.5*(2*J1+1)/(oq.j2+1) * (1-oq.occ) * xijap * yaqij ;

                       ///   chi_rk += 0.5*(2*J1+1)/(bra.op->j2+1) * ( xijra * ykaij  *dkr  );  // Includes the (1-Prs)

                    }
                }//for j

                for ( auto b : Z.modelspace->particles)
                {
                    Orbit& ob = Z.modelspace->GetOrbit(b);

                    int J1min = std::max( std::abs( oa.j2-ob.j2), std::abs(oi.j2-oq.j2) )/2;
                    int J1max = std::min(  oa.j2+ob.j2, oi.j2+oq.j2 )/2;

                    for ( int J1=J1min; J1<=J1max; J1++)
                    {

                       double xipab = X.TwoBody.GetTBME_J(J1,J1,i,p,a,b);
                       double xabiq = X.TwoBody.GetTBME_J(J1,J1,a,b,i,q);
                       double yabiq = Y.TwoBody.GetTBME_J(J1,J1,a,b,i,q);

                       chi_pq  += 0.5*(2*J1+1)/(oq.j2+1) * op.occ * xipab * xabiq ;
                       chiY_pq += 0.5*(2*J1+1)/(oq.j2+1) * op.occ * xipab * yabiq ;
                    }
                }//for b

              }//for i
            }//for a
            CHI_XX(p,q) = chi_pq;
            CHI_XY(p,q) = chiY_pq;
            }//for q
        }// for p





     std::vector<size_t> ch_bra_list,ch_ket_list;
     for ( auto& iter : Z.TwoBody.MatEl )
     {
        ch_bra_list.push_back( iter.first[0] );
        ch_ket_list.push_back( iter.first[1] );
     }
     size_t nch = ch_bra_list.size();

  //   int nch = Z.modelspace->GetNumberTwoBodyChannels();
     #pragma omp parallel for schedule(dynamic,1)
     for (int ich=0; ich<nch; ich++)
     {
       size_t ch_bra = ch_bra_list[ich];
       size_t ch_ket = ch_ket_list[ich];
  //     size_t ch_bra = itmat.first[0];
  //     size_t ch_ket = itmat.first[1];
       TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
       TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
       int J = tbc_bra.J;
       size_t nbras = tbc_bra.GetNumberKets();
       size_t nkets = tbc_ket.GetNumberKets();
       for ( size_t ibra=0; ibra<nbras; ibra++)
       {
         Ket& bra = tbc_bra.GetKet(ibra);
         size_t p = bra.p;
         size_t q = bra.q;
         int phasepq = bra.Phase(J);
  //       for ( size_t iket=0; iket<nkets; iket++)
         for ( size_t iket=ibra; iket<nkets; iket++)
         {
            Ket& ket = tbc_ket.GetKet(iket);
            size_t r = ket.p;
            size_t s = ket.q;
            int phasers = ket.Phase(J);
            double zpqrs = 0;

           for ( auto b : Z.modelspace->particles)
           {
             Orbit& ob = Z.modelspace->GetOrbit(b);


              zpqrs +=  Y.TwoBody.GetTBME_J(J,J,p,q,b,s) * CHI_XX(r,b) +  Y.TwoBody.GetTBME_J(J,J,p,q,r,b)* CHI_XX(s,b);
              zpqrs +=  Y.TwoBody.GetTBME_J(J,J,b,q,r,s) * CHI_XX(p,b) +  Y.TwoBody.GetTBME_J(J,J,p,b,r,s)* CHI_XX(q,b);
              zpqrs -=  X.TwoBody.GetTBME_J(J,J,p,q,b,s) * CHI_XY(r,b) +  X.TwoBody.GetTBME_J(J,J,p,q,r,b)* CHI_XY(s,b);
              zpqrs -=  X.TwoBody.GetTBME_J(J,J,b,q,r,s) * CHI_XY(p,b) +  X.TwoBody.GetTBME_J(J,J,p,b,r,s)* CHI_XY(q,b);

            }//for a

            for ( auto j : Z.modelspace->holes)
            {
                Orbit& oj = Z.modelspace->GetOrbit(j);


                zpqrs += Y.TwoBody.GetTBME_J(J,J,j,q,r,s)* CHI_XX(p,j) +  Y.TwoBody.GetTBME_J(J,J,p,j,r,s)* CHI_XX(q,j);
                zpqrs += Y.TwoBody.GetTBME_J(J,J,p,q,j,s)* CHI_XX(r,j) +  Y.TwoBody.GetTBME_J(J,J,p,q,r,j)* CHI_XX(s,j);
                zpqrs -= X.TwoBody.GetTBME_J(J,J,j,q,r,s)* CHI_XY(p,j) +  X.TwoBody.GetTBME_J(J,J,p,j,r,s)* CHI_XY(q,j);
                zpqrs -= X.TwoBody.GetTBME_J(J,J,p,q,j,s)* CHI_XY(r,j) +  X.TwoBody.GetTBME_J(J,J,p,q,r,j)* CHI_XY(s,j);

            }

            // normalize
            if (p==q ) zpqrs /= PhysConst::SQRT2;
            if (r==s ) zpqrs /= PhysConst::SQRT2;
            Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
         }// for iket
       }// for ibra

    }// for itmat


    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }
  */

  // [Omega, [Omega, H]]
  void comm223_231_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {
    Z.modelspace->PreCalculateSixJ();
    bool EraseOB = false;
     EraseOB = true;

    int orbit_print = 2;

    // determine symmetry
    int hEta = Eta.IsHermitian() ? 1 : -1;
    int hGamma = Gamma.IsHermitian() ? 1 : -1;
    // int hZ = Z.IsHermitian() ? 1 : -1;
    int hZ = hGamma;
    // ####################################################################################
    //   diagram I
    //
    //   I_pq = 1/2 \delta_{jp jq} / (2jp + 1) sum_{abcde J0 J1} \delta_{je jd}
    //          (\barn_a \barn_c nd nb -\barn_b \barn_d na nc - \barn_b \barn_e na nc
    //          + \barn_a \barn_c nb ne ) (2J_0 + 1) (2J_1 + 1) / (2jd + 1)
    //          eta^J0_bdac  eta^J0_acbe   Gamma^J1_epdq
    // ####################################################################################
    for (auto &p : Z.modelspace->all_orbits)
    {
      Orbit &op = Z.modelspace->GetOrbit(p);
      for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
      {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        double zij = 0;

        // loop abcde
        for (auto &a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          for (auto &b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double n_b = ob.occ;
            double nbar_b = 1.0 - n_b;
            for (auto &c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                double n_d = od.occ;
                double nbar_d = 1.0 - n_d;

                for (auto &e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);
                  if (oe.j2 != od.j2)
                    continue;
                  double n_e = oe.occ;
                  double nbar_e = 1.0 - n_e;
                  double occfactor = (nbar_a * nbar_c * n_b * n_d - nbar_b * nbar_d * n_a * n_c - nbar_b * nbar_e * n_a * n_c + nbar_a * nbar_c * n_b * n_e);
                  if (std::abs(occfactor) < 1e-6)
                    continue;
                  // condition

                  if ( (ob.l+od.l+oa.l+oc.l)%2 != Eta.GetParity() ) continue;
                  if ( (oa.l+oc.l+ob.l+oe.l)%2 != Eta.GetParity() ) continue;
                  if ( (oe.l+op.l+od.l+oq.l)%2 != Gamma.GetParity() ) continue;

                  if ( std::abs(ob.tz2+od.tz2-oa.tz2-oc.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(oa.tz2+oc.tz2-ob.tz2-oe.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(oe.tz2+op.tz2-od.tz2-oq.tz2) != Gamma.GetTRank() ) continue;

                  int J0min = std::max({std::abs(oa.j2 - oc.j2), std::abs(ob.j2 - od.j2), std::abs(ob.j2 - oe.j2)}) / 2;
                  int J0max = std::min({oa.j2 + oc.j2, ob.j2 + od.j2, ob.j2 + oe.j2}) / 2;

                  int J1min = std::max({std::abs(oe.j2 - op.j2), std::abs(od.j2 - oq.j2)}) / 2;
                  int J1max = std::min({oe.j2 + op.j2, od.j2 + oq.j2}) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    for (int J1 = J1min; J1 <= J1max; J1++)
                    {
                      zij += (2 * J0 + 1) * (2 * J1 + 1) / (od.j2 + 1.0) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, d, a, c) * Eta.TwoBody.GetTBME_J(J0, J0, a, c, b, e) * Gamma.TwoBody.GetTBME_J(J1, J1, e, p, d, q);
                      if (p==orbit_print and q==orbit_print)
                      {
                         std::cout << " I " << a << " " << b << " " << c << " " << d << "  " << J0 << " " << J1 << "   "
                                   << (2 * J0 + 1) * (2 * J1 + 1) / (od.j2 + 1.0) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, d, a, c) * Eta.TwoBody.GetTBME_J(J0, J0, a, c, b, e) * Gamma.TwoBody.GetTBME_J(J1, J1, e, p, d, q)
                                   << std::endl;
                      }
                    }
                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) += 0.5 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += 0.5 * hZ * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q

    }   // for p
//    std::cout << "diagram I  " << Z.OneBodyNorm() << std::endl;
//    std::cout << Z.OneBody << std::endl;
    if(EraseOB)
      Z.EraseOneBody();

    // ####################################################################################
    //  diagram II_a
    //
    //   IIa_pq = - \delta_{jp jq} / (2jp + 1) sum_{abcde J0 J1 J2 J3}
    //            ( \barn_a \barn_c nd nb - \barn_b \barn_d na nc )
    //            (2J_0 + 1) (2J_1 + 1)  (2J_2 + 1)  (2J_3 + 1)
    //            { ja jb J3 } { jp je J3 } { jb jp J2 }
    //            { jd jc J0 } { jd jc J1 } { je ja J3 }
    //            eta^J0_dbac  eta^J1_pcde   Gamma^J2_aeqb
    // ####################################################################################
    for (auto &p : Z.modelspace->all_orbits)
    {
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
        for (auto &a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double ja = oa.j2 / 2.;
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          for (auto &b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double jb = ob.j2 / 2.;
            double n_b = ob.occ;
            double nbar_b = 1.0 - n_b;
            for (auto &c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2 / 2.;
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                double jd = od.j2 / 2.;
                double n_d = od.occ;
                double nbar_d = 1.0 - n_d;

                for (auto &e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);
                  double je = oe.j2 / 2.;
                  double n_e = oe.occ;
                  double nbar_e = 1.0 - n_e;
                  double occfactor = (nbar_a * nbar_c * n_b * n_d - nbar_b * nbar_d * n_a * n_c);
                  if (std::abs(occfactor) < 1e-6)
                    continue;

                  if ( (ob.l+od.l+oa.l+oc.l)%2 != Eta.GetParity() ) continue;
                  if ( (oc.l+op.l+od.l+oe.l)%2 != Eta.GetParity() ) continue;
                  if ( (oa.l+oe.l+ob.l+oq.l)%2 != Gamma.GetParity() ) continue;

                  if ( std::abs(ob.tz2+od.tz2-oa.tz2-oc.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(oc.tz2+op.tz2-od.tz2-oe.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(oa.tz2+oe.tz2-ob.tz2-oq.tz2) != Gamma.GetTRank() ) continue;

                  int J0min = std::max({std::abs(ob.j2 - od.j2), std::abs(oa.j2 - oc.j2)}) / 2;
                  int J0max = std::min({ob.j2 + od.j2, oa.j2 + oc.j2}) / 2;

                  int J1min = std::max({std::abs(oc.j2 - op.j2), std::abs(od.j2 - oe.j2)}) / 2;
                  int J1max = std::min({oc.j2 + op.j2, od.j2 + oe.j2}) / 2;

                  int J2min = std::max({std::abs(oe.j2 - oa.j2), std::abs(ob.j2 - oq.j2)}) / 2;
                  int J2max = std::min({oe.j2 + oa.j2, ob.j2 + oq.j2}) / 2;

                  int J3min = std::max({std::abs(oa.j2 - ob.j2), std::abs(op.j2 - oe.j2), std::abs(oc.j2 - od.j2)}) / 2;
                  int J3max = std::min({oa.j2 + ob.j2, op.j2 + oe.j2, od.j2 + oc.j2}) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    for (int J1 = J1min; J1 <= J1max; J1++)
                    {
                      for (int J2 = J2min; J2 <= J2max; J2++)
                      {
                        for (int J3 = J3min; J3 <= J3max; J3++)
                        {
                          double phasefactor = Z.modelspace->phase(J0 + J1 + J2 + (oc.j2 + od.j2) / 2);
                          double sixj = Z.modelspace->GetSixJ(ja, jb, J3, jd, jc, J0);
                          sixj *= Z.modelspace->GetSixJ(jp, je, J3, jd, jc, J1);
                          sixj *= Z.modelspace->GetSixJ(jb, jp, J2, je, ja, J3);
                          zij += phasefactor * (2 * J0 + 1) * (2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * sixj * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, d, a, c) * Eta.TwoBody.GetTBME_J(J1, J1, c, p, d, e) * Gamma.TwoBody.GetTBME_J(J2, J2, a, e, b, q);
                      if (p==orbit_print and q==orbit_print)
                      {
                         std::cout << " IIa " << a << " " << b << " " << c << " " << d << "  " << J0 << " " << J1 << "   "
                                   << phasefactor * (2 * J0 + 1) * (2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * sixj * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, d, a, c) * Eta.TwoBody.GetTBME_J(J1, J1, c, p, d, e) * Gamma.TwoBody.GetTBME_J(J2, J2, a, e, b, q)
                                   << std::endl;
                      }
                        }
                      }
                    }
                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) += zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += hZ * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q

    }   // for p
//    std::cout << "diagram IIa " << Z.OneBodyNorm() << std::endl;
//    std::cout << Z.OneBody << std::endl;
    if(EraseOB)
      Z.EraseOneBody();

    // ####################################################################################
    //  diagram II_b
    //
    //   Ib_pq = 1/4 \delta_{jp jq} / (2jp + 1) sum_{abcde J0}
    //          (\barn_a \barn_d nb ne -\barn_b \barn_e na nd ) (2J_0 + 1)
    //          eta^J0_bdac  eta^J0_pcbe   Gamma^J0_adcq
    // ####################################################################################
    for (auto &p : Z.modelspace->all_orbits)
    {
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
        for (auto &a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double ja = oa.j2 / 2.;
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          for (auto &b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double jb = ob.j2 / 2.;
            double n_b = ob.occ;
            double nbar_b = 1.0 - n_b;
            for (auto &c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2 / 2.;
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                double jd = od.j2 / 2.;
                double n_d = od.occ;
                double nbar_d = 1.0 - n_d;

                for (auto &e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);
                  double je = oe.j2 / 2.;
                  double n_e = oe.occ;
                  double nbar_e = 1.0 - n_e;
                  double occfactor = (nbar_a * nbar_d * n_b * n_e - nbar_b * nbar_e * n_a * n_d);
                  if (std::abs(occfactor) < 1e-6)
                    continue;
                  if ( (ob.l+oe.l+oa.l+od.l)%2 != Eta.GetParity() ) continue;
                  if ( (oc.l+op.l+ob.l+oe.l)%2 != Eta.GetParity() ) continue;
                  if ( (oa.l+od.l+oc.l+oq.l)%2 != Gamma.GetParity() ) continue;

                  if ( std::abs(ob.tz2+oe.tz2-oa.tz2-od.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(oc.tz2+op.tz2-ob.tz2-oe.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(oa.tz2+od.tz2-oc.tz2-oq.tz2) != Gamma.GetTRank() ) continue;

                  int J0min = std::max({std::abs(oa.j2 - od.j2), std::abs(ob.j2 - oe.j2), std::abs(oc.j2 - op.j2)}) / 2;
                  int J0max = std::min({oa.j2 + od.j2, ob.j2 + oe.j2, op.j2 + oc.j2}) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    zij += (2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, e, a, d) * Eta.TwoBody.GetTBME_J(J0, J0, c, p, b, e) * Gamma.TwoBody.GetTBME_J(J0, J0, a, d, c, q);
                      if (p==orbit_print and q==orbit_print)
                      {
                         std::cout << " IIb " << a << " " << b << " " << c << " " << d << "  " << e << " " << J0 << "   "
                                   << (2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, e, a, d) * Eta.TwoBody.GetTBME_J(J0, J0, c, p, b, e) * Gamma.TwoBody.GetTBME_J(J0, J0, a, d, c, q)
                                   << std::endl;
                      }

                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) += 0.25 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += hZ * 0.25 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
//    std::cout << "diagram IIb " << Z.OneBodyNorm() << std::endl;
//    std::cout << Z.OneBody << std::endl;
    if(EraseOB)
      Z.EraseOneBody();

    // ####################################################################################
    //  diagram II_c
    //
    //   IIc_pq = \delta_{jp jq} / (2jp + 1) sum_{abcde J0 J1 J2 J3}
    //            ( \barn_a \barn_d nb ne - \barn_b \barn_e na nd )
    //            (2J_0 + 1) (2J_1 + 1)  (2J_2 + 1)  (2J_3 + 1)
    //            { jd je J3 } { jc jp J3 } { je jc J2 }
    //            { jb ja J0 } { jb ja J1 } { jp jd J3 }
    //            eta^J0_beda  eta^J1_cabq  Gamma^J2_dpce
    // ####################################################################################
    for (auto &p : Z.modelspace->all_orbits)
    {
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
        for (auto &a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double ja = oa.j2 / 2.;
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          for (auto &b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double jb = ob.j2 / 2.;
            double n_b = ob.occ;
            double nbar_b = 1.0 - n_b;
            for (auto &c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2 / 2.;
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                double jd = od.j2 / 2.;
                double n_d = od.occ;
                double nbar_d = 1.0 - n_d;

                for (auto &e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);
                  double je = oe.j2 / 2.;
                  double n_e = oe.occ;
                  double nbar_e = 1.0 - n_e;
                  double occfactor = (nbar_a * nbar_d * n_b * n_e - nbar_b * nbar_e * n_a * n_d);
                  if (std::abs(occfactor) < 1e-6)
                    continue;

                  if ( (ob.l+oe.l+od.l+oa.l)%2 != Eta.GetParity() ) continue;
                  if ( (oc.l+oa.l+ob.l+oq.l)%2 != Eta.GetParity() ) continue;
                  if ( (od.l+op.l+oc.l+oe.l)%2 != Gamma.GetParity() ) continue;

                  if ( std::abs(ob.tz2+oe.tz2-od.tz2-oa.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(oc.tz2+oa.tz2-ob.tz2-oq.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(od.tz2+op.tz2-oc.tz2-oe.tz2) != Gamma.GetTRank() ) continue;

                  int J0min = std::abs(oa.j2 - od.j2) / 2;
                  int J0max = (oa.j2 + od.j2) / 2;

                  int J1min = std::abs(oa.j2 - oc.j2) / 2;
                  int J1max = (oa.j2 + oc.j2) / 2;

                  int J2min = std::abs(oc.j2 - oe.j2) / 2;
                  int J2max = (oc.j2 + oe.j2) / 2;

                  int J3min = std::abs(od.j2 - oe.j2) / 2;
                  int J3max = (od.j2 + oe.j2) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    for (int J1 = J1min; J1 <= J1max; J1++)
                    {
                      for (int J2 = J2min; J2 <= J2max; J2++)
                      {
                        for (int J3 = J3min; J3 <= J3max; J3++)
                        {
                          double sixj = Z.modelspace->GetSixJ(jd, je, J3, jb, ja, J0);
                          sixj *= Z.modelspace->GetSixJ(jc, jp, J3, jb, ja, J1);
                          sixj *= Z.modelspace->GetSixJ(je, jc, J2, jp, jd, J3);
                          zij += (2 * J0 + 1) * (2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * sixj * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, e, d, a) * Eta.TwoBody.GetTBME_J(J1, J1, c, a, b, q) * Gamma.TwoBody.GetTBME_J(J2, J2, d, p, c, e);
                      if (p==orbit_print and q==orbit_print)
                      {
                         std::cout << " IIc " << a << " " << b << " " << c << " " << d << "  " << J0 << " " << J1 << "   "
                                   << (2 * J0 + 1) * (2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * sixj * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, e, d, a) * Eta.TwoBody.GetTBME_J(J1, J1, c, a, b, q) * Gamma.TwoBody.GetTBME_J(J2, J2, d, p, c, e)
                                   << std::endl;
                      }
                        }
                      }
                    }
                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) += zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += hZ * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q

    }   // for p
//    std::cout << "diagram IIc " << Z.OneBodyNorm() << std::endl;
//    std::cout << Z.OneBody << std::endl;
    if(EraseOB)
      Z.EraseOneBody();

    // ####################################################################################
    //  diagram II_d
    //
    //   IId_pq = - 1/4 \delta_{jp jq} / (2jp + 1) sum_{abcde J0}
    //          (\barn_c \barn_d na ne -\barn_a \barn_e nc nd ) (2J_0 + 1)
    //          eta^J0_aecd  eta^J0_cdbq  Gamma^J0_bqae
    // ####################################################################################
    for (auto &p : Z.modelspace->all_orbits)
    {
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
        for (auto &a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double ja = oa.j2 / 2.;
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          for (auto &b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double jb = ob.j2 / 2.;
            double n_b = ob.occ;
            double nbar_b = 1.0 - n_b;
            for (auto &c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2 / 2.;
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                double jd = od.j2 / 2.;
                double n_d = od.occ;
                double nbar_d = 1.0 - n_d;

                for (auto &e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);
                  double je = oe.j2 / 2.;
                  double n_e = oe.occ;
                  double nbar_e = 1.0 - n_e;
                  double occfactor = (nbar_c * nbar_d * n_a * n_e - nbar_a * nbar_e * n_c * n_d);
                  if (std::abs(occfactor) < 1e-6)
                    continue;

                  int J0min = std::abs(oa.j2 - oe.j2) / 2;
                  int J0max = (oa.j2 + oe.j2) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    zij += (2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, a, e, c, d) * Eta.TwoBody.GetTBME_J(J0, J0, c, d, b, q) * Gamma.TwoBody.GetTBME_J(J0, J0, b, p, a, e);
                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) -= 0.25 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) -= hZ * 0.25 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
//    std::cout << "diagram IId " << Z.OneBodyNorm() << std::endl;
//    std::cout << Z.OneBody << std::endl;
    if(EraseOB)
      Z.EraseOneBody();

    // ####################################################################################
    //  diagram III_a
    //
    //   IIIa_pq = 1/2 \delta_{jp jq} / (2jp + 1) sum_{abcde J0 J1}
    //          (\barn_a \barn_e nb nc -\barn_b \barn_c na ne ) (2J_0 + 1) (2J_1 + 1) / (2J_d + 1)
    //          eta^J0_bcae  eta^J1_epdq  Gamma^J0_adbc
    // ####################################################################################
    for (auto &p : Z.modelspace->all_orbits)
    {
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
        for (auto &a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double ja = oa.j2 / 2.;
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          for (auto &b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double jb = ob.j2 / 2.;
            double n_b = ob.occ;
            double nbar_b = 1.0 - n_b;
            for (auto &c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2 / 2.;
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                double jd = od.j2 / 2.;
                double n_d = od.occ;
                double nbar_d = 1.0 - n_d;

                for (auto &e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);
                  if (od.j2 != oe.j2)
                    continue;

                  double je = oe.j2 / 2.;
                  double n_e = oe.occ;
                  double nbar_e = 1.0 - n_e;
                  double occfactor = (nbar_a * nbar_e * n_b * n_c - nbar_b * nbar_c * n_a * n_e);
                  if (std::abs(occfactor) < 1e-6)
                    continue;

                  int J0min = std::abs(ob.j2 - oc.j2) / 2;
                  int J0max = (ob.j2 + oc.j2) / 2;

                  int J1min = std::abs(oe.j2 - op.j2) / 2;
                  int J1max = (oe.j2 + op.j2) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    for (int J1 = J1min; J1 <= J1max; J1++)
                    {
                      zij += (2 * J0 + 1) * (2 * J1 + 1) / (od.j2 + 1.) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, c, a, e) * Eta.TwoBody.GetTBME_J(J1, J1, e, p, d, q) * Gamma.TwoBody.GetTBME_J(J0, J0, a, d, b, c);
                    }
                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) += 0.5 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += hZ * 0.5 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
//    std::cout << "diagram IIIa " << Z.OneBodyNorm() << std::endl;
//    std::cout << Z.OneBody << std::endl;
    if(EraseOB)
      Z.EraseOneBody();

    // ###########################################################
    //  diagram III_b
    //
    //   IIIb_pq = - 1/2 \delta_{jp jq} / (2jp + 1) sum_{abcde J0 J1} \delta_{je jd}
    //          (\barn_a \barn_c nb nd -\barn_b \barn_d na nc ) (2J_0 + 1) (2J_1 + 1) / (2J_d + 1)
    //          eta^J0_bdac  eta^J1_epdq  Gamma^J0_acbe
    // ####################################################################################
    for (auto &p : Z.modelspace->all_orbits)
    {
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
        for (auto &a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double ja = oa.j2 / 2.;
          double n_a = oa.occ;
          double nbar_a = 1.0 - n_a;

          for (auto &b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double jb = ob.j2 / 2.;
            double n_b = ob.occ;
            double nbar_b = 1.0 - n_b;
            for (auto &c : Z.modelspace->all_orbits)
            {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2 / 2.;
              double n_c = oc.occ;
              double nbar_c = 1.0 - n_c;

              for (auto &d : Z.modelspace->all_orbits)
              {
                Orbit &od = Z.modelspace->GetOrbit(d);
                double jd = od.j2 / 2.;
                double n_d = od.occ;
                double nbar_d = 1.0 - n_d;

                for (auto &e : Z.modelspace->all_orbits)
                {
                  Orbit &oe = Z.modelspace->GetOrbit(e);
                  if (od.j2 != oe.j2)
                    continue;

                  double je = oe.j2 / 2.;
                  double n_e = oe.occ;
                  double nbar_e = 1.0 - n_e;
                  double occfactor = (nbar_a * nbar_c * n_b * n_d - nbar_b * nbar_d * n_a * n_c);
                  if (std::abs(occfactor) < 1e-6)
                    continue;

                  int J0min = std::abs(ob.j2 - od.j2) / 2;
                  int J0max = (ob.j2 + od.j2) / 2;

                  int J1min = std::abs(oe.j2 - op.j2) / 2;
                  int J1max = (oe.j2 + op.j2) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    for (int J1 = J1min; J1 <= J1max; J1++)
                    {
                      zij += (2 * J0 + 1) * (2 * J1 + 1) / (od.j2 + 1.) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, d, a, c) * Eta.TwoBody.GetTBME_J(J1, J1, e, p, d, q) * Gamma.TwoBody.GetTBME_J(J0, J0, a, c, b, e);
                    }
                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) -= 0.5 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) -= hZ * 0.5 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
//    std::cout << "diagram IIIb " << Z.OneBodyNorm() << std::endl;
//    std::cout << Z.OneBody << std::endl;
    if(EraseOB)
      Z.EraseOneBody();

    return;
  }

  void comm223_232_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {
    // global variables
    Z.modelspace->PreCalculateSixJ();
    int nch = Z.modelspace->GetNumberTwoBodyChannels(); // number of TB channels
    int norbits = Z.modelspace->all_orbits.size();
    std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
    auto &Z2 = Z.TwoBody;
    bool EraseTB = false;
     EraseTB = true;

    // determine symmetry
    int hEta = Eta.IsHermitian() ? 1 : -1;
    int hGamma = Gamma.IsHermitian() ? 1 : -1;
    // int hZ = Z.IsHermitian() ? 1 : -1;
    int hZ = hGamma;
    // ####################################################################################
    //   diagram Ia
    //
    //   I(a)^J0_pgqh = 1/2 * P_pg  1/ (2 jp + 1) \sum_abcd J2 \delta_{jd, jp}
    //                   ( \bar{n_a} \bar{n_c} n_b + \bar{n_b} n_a n_c )
    //                   ( 2 * J2 + 1 ) * eta^J2_bpac eta^J2_acbd Gamma^J0_dgqh
    // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);

        int phase_pg = bra.Phase(J0);
        double denominator_p = (op.j2 + 1.0);
        double denominator_g = (og.j2 + 1.0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  bool delta_jpjd = (od.j2 == op.j2);
                  bool delta_jgjd = (od.j2 == og.j2);

                  int j2min = std::max({std::abs(oa.j2 - oc.j2), std::abs(ob.j2 - od.j2), std::abs(ob.j2 - op.j2)}) / 2;
                  int j2max = std::min({oa.j2 + oc.j2, ob.j2 + od.j2, ob.j2 + op.j2}) / 2;
                  double occfactor = (nbar_a * nbar_c * n_b + nbar_b * n_a * n_c);

                  if (fabs(occfactor) < 1.e-7)
                    continue;

                  if (delta_jpjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
                      zpgqh += occfactor * (2 * J2 + 1) / denominator_p * Eta.TwoBody.GetTBME_J(J2, b, p, a, c) * Eta.TwoBody.GetTBME_J(J2, a, c, b, d) * Gamma.TwoBody.GetTBME_J(J0, d, g, q, h);
                      //                    zpgqh -= occfactor * (2 * J2 + 1) / denominator_p * Eta.TwoBody.GetTBME_J(J2, b, p, a, c) * Eta.TwoBody.GetTBME_J(J0, d, g, q, h) * Gamma.TwoBody.GetTBME_J(J2, a, c, b, d);
                    }

                  // exchanging  p <-> g
                  j2min = std::max({std::abs(oa.j2 - oc.j2), std::abs(ob.j2 - od.j2), std::abs(ob.j2 - og.j2)}) / 2;
                  j2max = std::min({oa.j2 + oc.j2, ob.j2 + od.j2, ob.j2 + og.j2}) / 2;
                  if (delta_jgjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
                      // zpgqh += phase_pg * occfactor * (2 * J2 + 1) / denominator_g * Eta.TwoBody.GetTBME_J(J2, b, g, a, c) * Eta.TwoBody.GetTBME_J(J2, a, c, b, d) * Gamma.TwoBody.GetTBME_J(J0, d, p, q, h);
                      zpgqh += occfactor * (2 * J2 + 1) / denominator_g * Eta.TwoBody.GetTBME_J(J2, b, g, a, c) * Eta.TwoBody.GetTBME_J(J2, a, c, b, d) * Gamma.TwoBody.GetTBME_J(J0, p, d, q, h);
                    }
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram Ia   " << Z.TwoBodyNorm() << std::endl;
    //Z.TwoBody.PrintMatrix(0,0);
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram Ib
      //
      //   I(b)^J0_pgqh = 1/2 *  P_qh  *  1/ (2 jq + 1) \sum_abcd J2  \delta_{jd, jq}
      //                   ( \bar{n_a} n_b n_c + \bar{n_b} \bar{n_c} n_a )
      //                   ( 2 * J2 + 1 ) * eta^J2_adbc eta^J2_bcaq Gamma^J0_pgdh
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();

      for (int iket = 0; iket < nKets; ++iket)
      {
        Ket &ket = tbc.GetKet(iket);
        size_t q = ket.p;
        size_t h = ket.q;
        Orbit &oq = *(ket.op);
        Orbit &oh = *(ket.oq);

        for (int ibra = 0; ibra <= iket; ++ibra)
        {
          Ket &bra = tbc.GetKet(ibra);
          size_t p = bra.p;
          size_t g = bra.q;
          Orbit &op = *(bra.op);
          Orbit &og = *(bra.oq);

          int phase_qh = ket.Phase(J0);
          double denominator_q = (oq.j2 + 1.0);
          double denominator_h = (oh.j2 + 1.0);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);

                  double occfactor = (nbar_a * n_b * n_c + nbar_b * nbar_c * n_a);
                  if (fabs(occfactor) < 1.e-7)
                    continue;

                  bool delta_jqjd = (od.j2 == oq.j2);
                  bool delta_jhjd = (od.j2 == oh.j2);

                  int j2min = std::max({std::abs(oa.j2 - od.j2) / 2, std::abs(ob.j2 - oc.j2) / 2, std::abs(oa.j2 - oq.j2) / 2});
                  int j2max = std::min({(oa.j2 + od.j2) / 2, (ob.j2 + oc.j2) / 2, (oa.j2 + oq.j2) / 2});

                  if (delta_jqjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
                      zpgqh += occfactor * (2 * J2 + 1) / denominator_q * Eta.TwoBody.GetTBME_J(J2, a, d, b, c) * Eta.TwoBody.GetTBME_J(J2, b, c, a, q) * Gamma.TwoBody.GetTBME_J(J0, p, g, d, h);
                    }

                  j2min = std::max({std::abs(oa.j2 - od.j2) / 2, std::abs(ob.j2 - oc.j2) / 2, std::abs(oa.j2 - oh.j2) / 2});
                  j2max = std::min({(oa.j2 + od.j2) / 2, (ob.j2 + oc.j2) / 2, (oa.j2 + oh.j2) / 2});
                  if (delta_jhjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
                      // zpgqh += phase_qh * occfactor * (2 * J2 + 1) / denominator_h * Eta.TwoBody.GetTBME_J(J2, a, d, b, c) * Eta.TwoBody.GetTBME_J(J2, b, c, a, h) * Gamma.TwoBody.GetTBME_J(J0, p, g, d, q);
                      zpgqh += occfactor * (2 * J2 + 1) / denominator_h * Eta.TwoBody.GetTBME_J(J2, a, d, b, c) * Eta.TwoBody.GetTBME_J(J2, b, c, a, h) * Gamma.TwoBody.GetTBME_J(J0, p, g, q, d);
                       if (ch==0 and ibra==2 and iket==4)
                       {
                          std::cout << " Ib " << occfactor * (2 * J2 + 1) / denominator_h * Eta.TwoBody.GetTBME_J(J2, a, d, b, c) * Eta.TwoBody.GetTBME_J(J2, b, c, a, h) * Gamma.TwoBody.GetTBME_J(J0, p, g, q, d)
                                    << std::endl;
                       }
                    }
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram Ib   " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IVa
      //
      //   IV(a)^J0_pgqh = - P(q/h) * 1/2 1/ (2jq + 1) \sum_abcd J2 \delta_{jd, jq}
      //                   ( \bar{n_a} n_b n_c + \bar{n_b} \bar{n_c} n_a )
      //                   ( 2 * J2 + 1 ) * eta^J2_bcaq eta^J0_pgdh Gamma^J2_adbc
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();

      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          int phase_qh = ket.Phase(J0);
          double denominator_q = (oq.j2 + 1.0);
          double denominator_h = (oh.j2 + 1.0);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);

                  double occfactor = (nbar_a * n_b * n_c + nbar_b * nbar_c * n_a);
                  if (fabs(occfactor) < 1.e-7)
                    continue;

                  int delta_jqjd = (od.j2 == oq.j2);
                  int delta_jhjd = (od.j2 == oh.j2);

                  int j2min = std::max({std::abs(oa.j2 - od.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(oa.j2 - oq.j2) / 2});
                  int j2max = std::min({(oa.j2 + od.j2) / 2, (oc.j2 + ob.j2) / 2, (oa.j2 + oq.j2) / 2});

                  if (delta_jqjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
                      zpgqh -= occfactor * (2 * J2 + 1) / denominator_q * Eta.TwoBody.GetTBME_J(J2, b, c, a, q) * Eta.TwoBody.GetTBME_J(J0, p, g, d, h) * Gamma.TwoBody.GetTBME_J(J2, a, d, b, c);
                    }

                  j2min = std::max({std::abs(oa.j2 - od.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(oa.j2 - oh.j2) / 2});
                  j2max = std::min({(oa.j2 + od.j2) / 2, (oc.j2 + ob.j2) / 2, (oa.j2 + oh.j2) / 2});

                  if (delta_jhjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
                      // zpgqh -= phase_qh * occfactor * (2 * J2 + 1) / denominator_h * Eta.TwoBody.GetTBME_J(J2, b, c, a, h) * Eta.TwoBody.GetTBME_J(J0, p, g, d, q) * Gamma.TwoBody.GetTBME_J(J2, a, d, b, c);
                      zpgqh -= occfactor * (2 * J2 + 1) / denominator_h * Eta.TwoBody.GetTBME_J(J2, b, c, a, h) * Eta.TwoBody.GetTBME_J(J0, p, g, q, d) * Gamma.TwoBody.GetTBME_J(J2, a, d, b, c);
                       if (ch==0 and ibra==2 and iket==4)
                       {
                          std::cout << " IVa " << -occfactor * (2 * J2 + 1) / denominator_h * Eta.TwoBody.GetTBME_J(J2, b, c, a, h) * Eta.TwoBody.GetTBME_J(J0, p, g, q, d) * Gamma.TwoBody.GetTBME_J(J2, a, d, b, c)
                                    << std::endl;
                       }
                    }

                  // *****************************************
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram IVa  " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IVb
      //
      //   IV(b)^J0_pgqh = -  P(p/g) * 1/2 1/ (2jp + 1) \sum_abcd J2 \delta_{jd, jp}
      //                   ( \bar{n_a} \bar{n_c} n_b + \bar{n_b} n_a n_c )
      //                   ( 2 * J2 + 1 ) * eta^J2_bpac eta^J0_dgqh Gamma^J2_acbd
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        int phase_pg = bra.Phase(J0);
        double denominator_p = (op.j2 + 1.0);
        double denominator_g = (og.j2 + 1.0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);

                  double occfactor = (nbar_a * nbar_c * n_b + nbar_b * n_a * n_c);
                  if (fabs(occfactor) < 1.e-7)
                    continue;

                  int delta_jpjd = (od.j2 == op.j2);
                  int delta_jgjd = (od.j2 == og.j2);

                  int j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(od.j2 - ob.j2) / 2, std::abs(ob.j2 - op.j2) / 2});
                  int j2max = std::min({(oa.j2 + oc.j2) / 2, (od.j2 + ob.j2) / 2, (ob.j2 + op.j2) / 2});

                  if (delta_jpjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
                      zpgqh -= occfactor * (2 * J2 + 1) / denominator_p * Eta.TwoBody.GetTBME_J(J2, b, p, a, c) * Eta.TwoBody.GetTBME_J(J0, d, g, q, h) * Gamma.TwoBody.GetTBME_J(J2, a, c, b, d);
                    }

                  j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(od.j2 - ob.j2) / 2, std::abs(ob.j2 - og.j2) / 2});
                  j2max = std::min({(oa.j2 + oc.j2) / 2, (od.j2 + ob.j2) / 2, (ob.j2 + og.j2) / 2});

                  if (delta_jgjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
                      // zpgqh -= phase_pg * occfactor * (2 * J2 + 1) / denominator_g * Eta.TwoBody.GetTBME_J(J2, b, g, a, c) * Eta.TwoBody.GetTBME_J(J0, d, p, q, h) * Gamma.TwoBody.GetTBME_J(J2, a, c, b, d);
                      zpgqh -= occfactor * (2 * J2 + 1) / denominator_g * Eta.TwoBody.GetTBME_J(J2, b, g, a, c) * Eta.TwoBody.GetTBME_J(J0, p, d, q, h) * Gamma.TwoBody.GetTBME_J(J2, a, c, b, d);
                    }
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zpgqh);

        } // iket

      }   // ibra
    }     // J0 channel
    //std::cout << "diagram IVb  " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIa
      //
      //   II(a)^J0_pgqh = - P_pg \sum_{abcd J2 J3 J4} ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                   { jd jg J4 } { jp ja J4 } { jg jp J0 }
      //                   { jc jb J2 } { jc jb J3 } { ja jd J4 }
      //
      //                   ( \bar{n_b} \bar{n_d} n_c + \bar{n_c} n_b n_d )
      //                   eta^J2_cgdb eta^J3_pbca Gamma^J0_daqh
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;

        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;

          double zpgqh = 0.;

          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_b * nbar_d * n_c + nbar_c * n_b * n_d);
                  if (fabs(occfactor) < 1.e-7)
                    continue;

                  if ( (oc.l+og.l+od.l+ob.l)%2 != Eta.GetParity() ) continue;
                  if ( (op.l+ob.l+oc.l+oa.l)%2 != Eta.GetParity() ) continue;
                  if ( (od.l+oa.l+oq.l+oh.l)%2 != Gamma.GetParity() ) continue;

                  if ( std::abs(oc.tz2+og.tz2-od.tz2-ob.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(op.tz2+ob.tz2-oc.tz2-oa.tz2) != Eta.GetTRank() ) continue;
                  if ( std::abs(od.tz2+oa.tz2-oq.tz2-oh.tz2) != Gamma.GetTRank() ) continue;

                  /// direct term
                  int j2min = std::max(std::abs(oc.j2 - og.j2), std::abs(ob.j2 - od.j2)) / 2;
                  int j2max = std::min(oc.j2 + og.j2, ob.j2 + od.j2) / 2;

                  int j3min = std::max(std::abs(oa.j2 - oc.j2), std::abs(ob.j2 - op.j2)) / 2;
                  int j3max = std::min(oa.j2 + oc.j2, ob.j2 + op.j2) / 2;

                  int j4min = std::max({std::abs(od.j2 - og.j2), std::abs(oa.j2 - op.j2), std::abs(oc.j2 - ob.j2)}) / 2;
                  int j4max = std::min({od.j2 + og.j2, oa.j2 + op.j2, oc.j2 + ob.j2}) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jd, jg, J4, jc, jb, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jp, ja, J4, jc, jb, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jg, jp, J0, ja, jd, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jd, jg, J4, jc, jb, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jp, ja, J4, jc, jb, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jg, jp, J0, ja, jd, J4);

                        zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, c, g, d, b) * Eta.TwoBody.GetTBME_J(J3, p, b, c, a) * Gamma.TwoBody.GetTBME_J(J0, d, a, q, h);
                      }
                    }
                  }

                  /// exchange term, exchange pg
                  j2min = std::max(std::abs(oc.j2 - op.j2), std::abs(ob.j2 - od.j2)) / 2;
                  j2max = std::min(oc.j2 + op.j2, ob.j2 + od.j2) / 2;

                  j3min = std::max(std::abs(oa.j2 - oc.j2), std::abs(ob.j2 - og.j2)) / 2;
                  j3max = std::min(oa.j2 + oc.j2, ob.j2 + og.j2) / 2;

                  j4min = std::max({std::abs(od.j2 - op.j2), std::abs(oa.j2 - og.j2), std::abs(oc.j2 - ob.j2)}) / 2;
                  j4max = std::min({od.j2 + op.j2, oa.j2 + og.j2, oc.j2 + ob.j2}) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jd, jp, J4, jc, jb, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jg, ja, J4, jc, jb, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jp, jg, J0, ja, jd, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jd, jp, J4, jc, jb, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jg, ja, J4, jc, jb, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jp, jg, J0, ja, jd, J4);

                        zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, c, p, d, b) * Eta.TwoBody.GetTBME_J(J3, g, b, c, a) * Gamma.TwoBody.GetTBME_J(J0, d, a, q, h);
                       if (ch==0 and ibra==2 and iket==4)
                       {
                          std::cout << " IIa " << -phase_pg * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, c, p, d, b) * Eta.TwoBody.GetTBME_J(J3, g, b, c, a) * Gamma.TwoBody.GetTBME_J(J0, d, a, q, h)
                                    << std::endl;
                       }
                      }
                    }
                  }
                  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
                } // d
              } // c
            } // b
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram IIa  " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIb
      //
      //   II(b)^J0_pghq = - P(p/g) * P(q/h) * \sum_{abcd J2 J3 J4 J5} ( 2 * J2 + 1 ) ( 2 * J3 + 1 )
      //                   ( 2 * J4 + 1 ) ( 2 * J5 + 1 )
      //
      //                   { jq jd J5 } { jp ja J5 } { J0 J5 J4 } { J0 J4 J5 }
      //                   { jc jb J2 } { jc jb J3 } { jd jh jq } { ja jp jg }
      //
      //                   ( \bar{n_b} n_c n_d + \bar{n_c} \bar{n_d} n_b )
      //                   eta^J2_dcbq eta^J3_bpac Gamma^J4_gahd
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_b * n_c * n_d + nbar_c * nbar_d * n_b);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  // direct term
                  int j2min = std::max(std::abs(oc.j2 - od.j2), std::abs(ob.j2 - oq.j2)) / 2;
                  int j2max = std::min(oc.j2 + od.j2, ob.j2 + oq.j2) / 2;

                  int j3min = std::max(std::abs(ob.j2 - op.j2), std::abs(oa.j2 - oc.j2)) / 2;
                  int j3max = std::min(ob.j2 + op.j2, oa.j2 + oc.j2) / 2;

                  int j4min = std::max(std::abs(og.j2 - oa.j2), std::abs(oh.j2 - od.j2)) / 2;
                  int j4max = std::min(og.j2 + oa.j2, oh.j2 + od.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oq.j2 - od.j2) / 2, std::abs(op.j2 - oa.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(J0 - J4)});
                        int j5max = std::min({(oq.j2 + od.j2) / 2, (op.j2 + oa.j2) / 2, (oc.j2 + ob.j2) / 2, J0 + J4});
                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(jq, jd, J5, jc, jb, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(jp, ja, J5, jc, jb, J3);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(J0, J5, J4, jd, jh, jq);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J4, J5, ja, jp, jg);

                          double sixj1 = Z.modelspace->GetSixJ(jq, jd, J5, jc, jb, J2);
                          double sixj2 = Z.modelspace->GetSixJ(jp, ja, J5, jc, jb, J3);
                          double sixj3 = Z.modelspace->GetSixJ(J0, J5, J4, jd, jh, jq);
                          double sixj4 = Z.modelspace->GetSixJ(J0, J4, J5, ja, jp, jg);

                          zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, d, c, b, q) * Eta.TwoBody.GetTBME_J(J3, b, p, a, c) * Gamma.TwoBody.GetTBME_J(J4, g, a, h, d);
                        }
                      }
                    }
                  }

                  // exchange q <-> h
                  j2min = std::max(std::abs(oc.j2 - od.j2), std::abs(ob.j2 - oh.j2)) / 2;
                  j2max = std::min(oc.j2 + od.j2, ob.j2 + oh.j2) / 2;

                  j3min = std::max(std::abs(ob.j2 - op.j2), std::abs(oa.j2 - oc.j2)) / 2;
                  j3max = std::min(ob.j2 + op.j2, oa.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(og.j2 - oa.j2), std::abs(oq.j2 - od.j2)) / 2;
                  j4max = std::min(og.j2 + oa.j2, oq.j2 + od.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oh.j2 - od.j2) / 2, std::abs(op.j2 - oa.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(J0 - J4)});
                        int j5max = std::min({(oh.j2 + od.j2) / 2, (op.j2 + oa.j2) / 2, (oc.j2 + ob.j2) / 2, J0 + J4});
                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(jh, jd, J5, jc, jb, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(jp, ja, J5, jc, jb, J3);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(J0, J5, J4, jd, jq, jh);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J4, J5, ja, jp, jg);

                          double sixj1 = Z.modelspace->GetSixJ(jh, jd, J5, jc, jb, J2);
                          double sixj2 = Z.modelspace->GetSixJ(jp, ja, J5, jc, jb, J3);
                          double sixj3 = Z.modelspace->GetSixJ(J0, J5, J4, jd, jq, jh);
                          double sixj4 = Z.modelspace->GetSixJ(J0, J4, J5, ja, jp, jg);

                          zpgqh -= phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, d, c, b, h) * Eta.TwoBody.GetTBME_J(J3, b, p, a, c) * Gamma.TwoBody.GetTBME_J(J4, g, a, q, d);
                        }
                      }
                    }
                  }

                  // exchange p <-> g
                  j2min = std::max(std::abs(oc.j2 - od.j2), std::abs(ob.j2 - oq.j2)) / 2;
                  j2max = std::min(oc.j2 + od.j2, ob.j2 + oq.j2) / 2;

                  j3min = std::max(std::abs(ob.j2 - og.j2), std::abs(oa.j2 - oc.j2)) / 2;
                  j3max = std::min(ob.j2 + og.j2, oa.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(op.j2 - oa.j2), std::abs(oh.j2 - od.j2)) / 2;
                  j4max = std::min(op.j2 + oa.j2, oh.j2 + od.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oq.j2 - od.j2) / 2, std::abs(og.j2 - oa.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(J0 - J4)});
                        int j5max = std::min({(oq.j2 + od.j2) / 2, (og.j2 + oa.j2) / 2, (oc.j2 + ob.j2) / 2, J0 + J4});
                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(jq, jd, J5, jc, jb, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(jg, ja, J5, jc, jb, J3);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(J0, J5, J4, jd, jh, jq);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J4, J5, ja, jg, jp);

                          double sixj1 = Z.modelspace->GetSixJ(jq, jd, J5, jc, jb, J2);
                          double sixj2 = Z.modelspace->GetSixJ(jg, ja, J5, jc, jb, J3);
                          double sixj3 = Z.modelspace->GetSixJ(J0, J5, J4, jd, jh, jq);
                          double sixj4 = Z.modelspace->GetSixJ(J0, J4, J5, ja, jg, jp);

                          zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, d, c, b, q) * Eta.TwoBody.GetTBME_J(J3, b, g, a, c) * Gamma.TwoBody.GetTBME_J(J4, p, a, h, d);
                        }
                      }
                    }
                  }

                  // exchange p <-> g  and q <-> h
                  j2min = std::max(std::abs(oc.j2 - od.j2), std::abs(ob.j2 - oh.j2)) / 2;
                  j2max = std::min(oc.j2 + od.j2, ob.j2 + oh.j2) / 2;

                  j3min = std::max(std::abs(ob.j2 - og.j2), std::abs(oa.j2 - oc.j2)) / 2;
                  j3max = std::min(ob.j2 + og.j2, oa.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(op.j2 - oa.j2), std::abs(oq.j2 - od.j2)) / 2;
                  j4max = std::min(op.j2 + oa.j2, oq.j2 + od.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oh.j2 - od.j2) / 2, std::abs(og.j2 - oa.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(J0 - J4)});
                        int j5max = std::min({(oh.j2 + od.j2) / 2, (og.j2 + oa.j2) / 2, (oc.j2 + ob.j2) / 2, J0 + J4});
                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(jh, jd, J5, jc, jb, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(jg, ja, J5, jc, jb, J3);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(J0, J5, J4, jd, jq, jh);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J4, J5, ja, jg, jp);

                          double sixj1 = Z.modelspace->GetSixJ(jh, jd, J5, jc, jb, J2);
                          double sixj2 = Z.modelspace->GetSixJ(jg, ja, J5, jc, jb, J3);
                          double sixj3 = Z.modelspace->GetSixJ(J0, J5, J4, jd, jq, jh);
                          double sixj4 = Z.modelspace->GetSixJ(J0, J4, J5, ja, jg, jp);

                          zpgqh -= phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, d, c, b, h) * Eta.TwoBody.GetTBME_J(J3, b, g, a, c) * Gamma.TwoBody.GetTBME_J(J4, p, a, q, d);
                        }
                      }
                    }
                  }

                  // ******************************************
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram IIb  " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIc
      //
      //   II(c)^J0_pgqh = - P_qh * \sum_{abcd J2 J3 J4} ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                   { jq jd J4 } { ja jh J4 } { jq jh J0 }
      //                   { jc jb J2 } { jc jb J3 } { ja jd J4 }
      //
      //                   ( \bar{n_b} n_c n_d + \bar{n_c} \bar{n_d} n_b )
      //                   eta^J2_cdqb eta^J3_abch Gamma^J0_pgad
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();

      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;

          int phase_qh = ket.Phase(J0);
          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_b * n_c * n_d + nbar_c * nbar_d * n_b);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  int j2min = std::max(std::abs(oc.j2 - od.j2), std::abs(ob.j2 - oq.j2)) / 2;
                  int j2max = std::min(oc.j2 + od.j2, ob.j2 + oq.j2) / 2;

                  int j3min = std::max(std::abs(oa.j2 - ob.j2), std::abs(oh.j2 - oc.j2)) / 2;
                  int j3max = std::min(oa.j2 + ob.j2, oh.j2 + oc.j2) / 2;

                  int j4min = std::max({std::abs(oq.j2 - od.j2), std::abs(oa.j2 - oh.j2), std::abs(oc.j2 - ob.j2)}) / 2;
                  int j4max = std::min({oq.j2 + od.j2, oa.j2 + oh.j2, oc.j2 + ob.j2}) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jq, jd, J4, jc, jb, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(ja, jh, J4, jc, jb, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, ja, jd, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jq, jd, J4, jc, jb, J2);
                        double sixj2 = Z.modelspace->GetSixJ(ja, jh, J4, jc, jb, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, ja, jd, J4);

                        zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, c, d, q, b) * Eta.TwoBody.GetTBME_J(J3, a, b, c, h) * Gamma.TwoBody.GetTBME_J(J0, p, g, a, d);
                      }
                    }
                  }

                  // exchanging q and h
                  j2min = std::max(std::abs(oc.j2 - od.j2), std::abs(ob.j2 - oh.j2)) / 2;
                  j2max = std::min(oc.j2 + od.j2, ob.j2 + oh.j2) / 2;

                  j3min = std::max(std::abs(oa.j2 - ob.j2), std::abs(oq.j2 - oc.j2)) / 2;
                  j3max = std::min(oa.j2 + ob.j2, oq.j2 + oc.j2) / 2;

                  j4min = std::max({std::abs(oh.j2 - od.j2), std::abs(oa.j2 - oq.j2), std::abs(oc.j2 - ob.j2)}) / 2;
                  j4max = std::min({oh.j2 + od.j2, oa.j2 + oq.j2, oc.j2 + ob.j2}) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jh, jd, J4, jc, jb, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(ja, jq, J4, jc, jb, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, ja, jd, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jh, jd, J4, jc, jb, J2);
                        double sixj2 = Z.modelspace->GetSixJ(ja, jq, J4, jc, jb, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, ja, jd, J4);

                        zpgqh -= phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, c, d, h, b) * Eta.TwoBody.GetTBME_J(J3, a, b, c, q) * Gamma.TwoBody.GetTBME_J(J0, p, g, a, d);
                      }
                    }
                  }

                  //  8888888888888888888888888888888888888888888888888888888888888
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram IIc  " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IId
      //
      //   II(d)^J0_pgqh = - P(p/g) * P(q/h) * \sum_{abcd J2 J3 J4 J5} ( 2 * J2 + 1 )
      //                   ( 2 * J3 + 1 ) ( 2 * J4 + 1 ) ( 2 * J5 + 1 )
      //
      //                   { jd jg J5 } { ja jh J5 } { J0 J5 J4 } { J5 J4 J0 }
      //                   { jc jb J2 } { jc jb J3 } { jd jp jg } { jq jh ja }
      //
      //                   ( \bar{n_c} n_b n_d - \bar{n_b} \bar{n_d} n_c )
      //                   eta^J2_gcbd eta^J3_bahc Gamma^J4_dpaq
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_c * n_b * n_d + nbar_b * nbar_d * n_c);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  // direct term
                  int j2min = std::max(std::abs(oc.j2 - og.j2), std::abs(ob.j2 - od.j2)) / 2;
                  int j2max = std::min(oc.j2 + og.j2, ob.j2 + od.j2) / 2;

                  int j3min = std::max(std::abs(ob.j2 - oa.j2), std::abs(oh.j2 - oc.j2)) / 2;
                  int j3max = std::min(ob.j2 + oa.j2, oh.j2 + oc.j2) / 2;

                  int j4min = std::max(std::abs(oq.j2 - oa.j2), std::abs(op.j2 - od.j2)) / 2;
                  int j4max = std::min(oq.j2 + oa.j2, op.j2 + od.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(og.j2 - od.j2) / 2, std::abs(oh.j2 - oa.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(J0 - J4)});
                        int j5max = std::min({(og.j2 + od.j2) / 2, (oh.j2 + oa.j2) / 2, (oc.j2 + ob.j2) / 2, J0 + J4});

                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(jq, jd, J5, jc, jb, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(jp, ja, J5, jc, jb, J3);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(J0, J5, J4, jd, jh, jq);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J4, J5, ja, jp, jg);

                          double sixj1 = Z.modelspace->GetSixJ(jd, jg, J5, jc, jb, J2);
                          double sixj2 = Z.modelspace->GetSixJ(ja, jh, J5, jc, jb, J3);
                          double sixj3 = Z.modelspace->GetSixJ(J0, J5, J4, jd, jp, jg);
                          double sixj4 = Z.modelspace->GetSixJ(J5, J4, J0, jq, jh, ja);

                          zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, g, c, b, d) * Eta.TwoBody.GetTBME_J(J3, b, a, h, c) * Gamma.TwoBody.GetTBME_J(J4, d, p, a, q);
                        }
                      }
                    }
                  }

                  // exchange q and h
                  j2min = std::max(std::abs(oc.j2 - og.j2), std::abs(ob.j2 - od.j2)) / 2;
                  j2max = std::min(oc.j2 + og.j2, ob.j2 + od.j2) / 2;

                  j3min = std::max(std::abs(ob.j2 - oa.j2), std::abs(oq.j2 - oc.j2)) / 2;
                  j3max = std::min(ob.j2 + oa.j2, oq.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(oh.j2 - oa.j2), std::abs(op.j2 - od.j2)) / 2;
                  j4max = std::min(oh.j2 + oa.j2, op.j2 + od.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(og.j2 - od.j2) / 2, std::abs(oq.j2 - oa.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(J0 - J4)});
                        int j5max = std::min({(og.j2 + od.j2) / 2, (oq.j2 + oa.j2) / 2, (oc.j2 + ob.j2) / 2, J0 + J4});

                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(jh, jd, J5, jc, jb, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(jp, ja, J5, jc, jb, J3);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(J0, J5, J4, jd, jq, jh);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J4, J5, ja, jp, jg);

                          double sixj1 = Z.modelspace->GetSixJ(jd, jg, J5, jc, jb, J2);
                          double sixj2 = Z.modelspace->GetSixJ(ja, jq, J5, jc, jb, J3);
                          double sixj3 = Z.modelspace->GetSixJ(J0, J5, J4, jd, jp, jg);
                          double sixj4 = Z.modelspace->GetSixJ(J5, J4, J0, jh, jq, ja);

                          zpgqh -= phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, g, c, b, d) * Eta.TwoBody.GetTBME_J(J3, b, a, q, c) * Gamma.TwoBody.GetTBME_J(J4, d, p, a, h);
                        }
                      }
                    }
                  }

                  // exchange p and g
                  j2min = std::max(std::abs(oc.j2 - op.j2), std::abs(ob.j2 - od.j2)) / 2;
                  j2max = std::min(oc.j2 + op.j2, ob.j2 + od.j2) / 2;

                  j3min = std::max(std::abs(ob.j2 - oa.j2), std::abs(oh.j2 - oc.j2)) / 2;
                  j3max = std::min(ob.j2 + oa.j2, oh.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(oq.j2 - oa.j2), std::abs(og.j2 - od.j2)) / 2;
                  j4max = std::min(oq.j2 + oa.j2, og.j2 + od.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(op.j2 - od.j2) / 2, std::abs(oh.j2 - oa.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(J0 - J4)});
                        int j5max = std::min({(op.j2 + od.j2) / 2, (oh.j2 + oa.j2) / 2, (oc.j2 + ob.j2) / 2, J0 + J4});

                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(jq, jd, J5, jc, jb, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(jg, ja, J5, jc, jb, J3);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(J0, J5, J4, jd, jh, jq);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J4, J5, ja, jg, jp);

                          double sixj1 = Z.modelspace->GetSixJ(jd, jp, J5, jc, jb, J2);
                          double sixj2 = Z.modelspace->GetSixJ(ja, jh, J5, jc, jb, J3);
                          double sixj3 = Z.modelspace->GetSixJ(J0, J5, J4, jd, jg, jp);
                          double sixj4 = Z.modelspace->GetSixJ(J5, J4, J0, jq, jh, ja);

                          zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, p, c, b, d) * Eta.TwoBody.GetTBME_J(J3, b, a, h, c) * Gamma.TwoBody.GetTBME_J(J4, d, g, a, q);
                        }
                      }
                    }
                  }

                  // exchange p g and q h
                  j2min = std::max(std::abs(oc.j2 - op.j2), std::abs(ob.j2 - od.j2)) / 2;
                  j2max = std::min(oc.j2 + op.j2, ob.j2 + od.j2) / 2;

                  j3min = std::max(std::abs(ob.j2 - oa.j2), std::abs(oq.j2 - oc.j2)) / 2;
                  j3max = std::min(ob.j2 + oa.j2, oq.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(oh.j2 - oa.j2), std::abs(og.j2 - od.j2)) / 2;
                  j4max = std::min(oh.j2 + oa.j2, og.j2 + od.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(op.j2 - od.j2) / 2, std::abs(oq.j2 - oa.j2) / 2, std::abs(oc.j2 - ob.j2) / 2, std::abs(J0 - J4)});
                        int j5max = std::min({(op.j2 + od.j2) / 2, (oq.j2 + oa.j2) / 2, (oc.j2 + ob.j2) / 2, J0 + J4});

                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(jh, jd, J5, jc, jb, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(jg, ja, J5, jc, jb, J3);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(J0, J5, J4, jd, jq, jh);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J4, J5, ja, jg, jp);

                          double sixj1 = Z.modelspace->GetSixJ(jd, jp, J5, jc, jb, J2);
                          double sixj2 = Z.modelspace->GetSixJ(ja, jq, J5, jc, jb, J3);
                          double sixj3 = Z.modelspace->GetSixJ(J0, J5, J4, jd, jg, jp);
                          double sixj4 = Z.modelspace->GetSixJ(J5, J4, J0, jh, jq, ja);

                          zpgqh -= phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, p, c, b, d) * Eta.TwoBody.GetTBME_J(J3, b, a, q, c) * Gamma.TwoBody.GetTBME_J(J4, d, g, a, h);
                        }
                      }
                    }
                  }
                }
              }
            }
          } // a
          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram IId  " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIe
      //
      //   II(e)^J0_pgqh = - 1/2 * P(p/g) * P(q/h) *  \sum_{abcd J2 J3 J4}
      //                     ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                   { jp jh J4 } { jq jg J4 } { jh jq J0 }
      //                   { jb jd J2 } { jb jd J3 } { jg jp J4 }
      //
      //                   ( \bar{n_b} n_a n_c + \bar{n_a} \bar{n_c} n_b )
      //                   eta^J2_acbh eta^J2_pdac Gamma^J3_bgqd
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;

          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_b * n_a * n_c + nbar_a * nbar_c * n_b);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  // direct term
                  int j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(ob.j2 - oh.j2) / 2, std::abs(op.j2 - od.j2) / 2});
                  int j2max = std::min({(oa.j2 + oc.j2) / 2, (ob.j2 + oh.j2) / 2, (op.j2 + od.j2) / 2});

                  int j3min = std::max(std::abs(ob.j2 - og.j2), std::abs(oq.j2 - od.j2)) / 2;
                  int j3max = std::min(ob.j2 + og.j2, oq.j2 + od.j2) / 2;

                  int j4min = std::max({std::abs(op.j2 - oh.j2) / 2, std::abs(oq.j2 - og.j2) / 2, std::abs(ob.j2 - od.j2) / 2});
                  int j4max = std::min({(op.j2 + oh.j2) / 2, (oq.j2 + og.j2) / 2, (ob.j2 + od.j2) / 2});

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jp, jh, J4, jb, jd, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jq, jg, J4, jb, jd, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, jg, jp, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jp, jh, J4, jb, jd, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jq, jg, J4, jb, jd, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, jg, jp, J4);

                        zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, a, c, b, h) * Eta.TwoBody.GetTBME_J(J2, p, d, a, c) * Gamma.TwoBody.GetTBME_J(J3, b, g, q, d);
                      }
                    }
                  }

                  // exchange q <-> h
                  j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(ob.j2 - oq.j2) / 2, std::abs(op.j2 - od.j2) / 2});
                  j2max = std::min({(oa.j2 + oc.j2) / 2, (ob.j2 + oq.j2) / 2, (op.j2 + od.j2) / 2});

                  j3min = std::max(std::abs(ob.j2 - og.j2), std::abs(oh.j2 - od.j2)) / 2;
                  j3max = std::min(ob.j2 + og.j2, oh.j2 + od.j2) / 2;

                  j4min = std::max({std::abs(op.j2 - oq.j2) / 2, std::abs(oh.j2 - og.j2) / 2, std::abs(ob.j2 - od.j2) / 2});
                  j4max = std::min({(op.j2 + oq.j2) / 2, (oh.j2 + og.j2) / 2, (ob.j2 + od.j2) / 2});

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jp, jq, J4, jb, jd, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jh, jg, J4, jb, jd, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, jg, jp, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jp, jq, J4, jb, jd, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jh, jg, J4, jb, jd, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, jg, jp, J4);

                        zpgqh -= phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, a, c, b, q) * Eta.TwoBody.GetTBME_J(J2, p, d, a, c) * Gamma.TwoBody.GetTBME_J(J3, b, g, h, d);
                      }
                    }
                  }

                  // exchange p <-> g
                  j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(ob.j2 - oh.j2) / 2, std::abs(og.j2 - od.j2) / 2});
                  j2max = std::min({(oa.j2 + oc.j2) / 2, (ob.j2 + oh.j2) / 2, (og.j2 + od.j2) / 2});

                  j3min = std::max(std::abs(ob.j2 - op.j2), std::abs(oq.j2 - od.j2)) / 2;
                  j3max = std::min(ob.j2 + op.j2, oq.j2 + od.j2) / 2;

                  j4min = std::max({std::abs(og.j2 - oh.j2) / 2, std::abs(oq.j2 - op.j2) / 2, std::abs(ob.j2 - od.j2) / 2});
                  j4max = std::min({(og.j2 + oh.j2) / 2, (oq.j2 + op.j2) / 2, (ob.j2 + od.j2) / 2});

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jg, jh, J4, jb, jd, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jq, jp, J4, jb, jd, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, jp, jg, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jg, jh, J4, jb, jd, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jq, jp, J4, jb, jd, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, jp, jg, J4);

                        zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, a, c, b, h) * Eta.TwoBody.GetTBME_J(J2, g, d, a, c) * Gamma.TwoBody.GetTBME_J(J3, b, p, q, d);
                      }
                    }
                  }

                  // exchange p <-> g  and q <-> h
                  j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(ob.j2 - oq.j2) / 2, std::abs(og.j2 - od.j2) / 2});
                  j2max = std::min({(oa.j2 + oc.j2) / 2, (ob.j2 + oq.j2) / 2, (og.j2 + od.j2) / 2});

                  j3min = std::max(std::abs(ob.j2 - op.j2), std::abs(oh.j2 - od.j2)) / 2;
                  j3max = std::min(ob.j2 + op.j2, oh.j2 + od.j2) / 2;

                  j4min = std::max({std::abs(og.j2 - oq.j2) / 2, std::abs(oh.j2 - op.j2) / 2, std::abs(ob.j2 - od.j2) / 2});
                  j4max = std::min({(og.j2 + oq.j2) / 2, (oh.j2 + op.j2) / 2, (ob.j2 + od.j2) / 2});

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jg, jq, J4, jb, jd, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jh, jp, J4, jb, jd, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, jp, jg, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jg, jq, J4, jb, jd, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jh, jp, J4, jb, jd, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, jp, jg, J4);

                        zpgqh -= phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, a, c, b, q) * Eta.TwoBody.GetTBME_J(J2, g, d, a, c) * Gamma.TwoBody.GetTBME_J(J3, b, p, h, d);
                      }
                    }
                  }
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram IIe  " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIf
      //
      //   II(f)^J0_pgqh = - 1/2 * P(p/g) * P(q/h) * \sum_{abcd J2 J3 J4}
      //                     ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                   { jh jp J4 } { jg jq J4 } { jh jq J0 }
      //                   { jb jd J2 } { jb jd J3 } { jg jp J4 }
      //
      //                   ( \bar{n_a} \bar{n_c} n_b + \bar{n_b} n_a n_c )
      //                   eta^J2_pbac eta^J2_acdh Gamma^J3_dgqb
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_a * nbar_c * n_b + nbar_b * n_a * n_c);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  // direct term
                  int j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(od.j2 - oh.j2) / 2, std::abs(op.j2 - ob.j2) / 2});
                  int j2max = std::min({(oa.j2 + oc.j2) / 2, (od.j2 + oh.j2) / 2, (op.j2 + ob.j2) / 2});

                  int j3min = std::max(std::abs(od.j2 - og.j2), std::abs(oq.j2 - ob.j2)) / 2;
                  int j3max = std::min(od.j2 + og.j2, oq.j2 + ob.j2) / 2;

                  int j4min = std::max({std::abs(op.j2 - oh.j2) / 2, std::abs(oq.j2 - og.j2) / 2, std::abs(ob.j2 - od.j2) / 2});
                  int j4max = std::min({(op.j2 + oh.j2) / 2, (oq.j2 + og.j2) / 2, (ob.j2 + od.j2) / 2});

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jh, jp, J4, jb, jd, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jg, jq, J4, jb, jd, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, jg, jp, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jh, jp, J4, jb, jd, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jg, jq, J4, jb, jd, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, jg, jp, J4);

                        zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, p, b, a, c) * Eta.TwoBody.GetTBME_J(J2, a, c, d, h) * Gamma.TwoBody.GetTBME_J(J3, d, g, q, b);
                      }
                    }
                  }

                  // exchange h <-> q
                  j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(od.j2 - oq.j2) / 2, std::abs(op.j2 - ob.j2) / 2});
                  j2max = std::min({(oa.j2 + oc.j2) / 2, (od.j2 + oq.j2) / 2, (op.j2 + ob.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - og.j2), std::abs(oh.j2 - ob.j2)) / 2;
                  j3max = std::min(od.j2 + og.j2, oh.j2 + ob.j2) / 2;

                  j4min = std::max({std::abs(op.j2 - oq.j2) / 2, std::abs(oh.j2 - og.j2) / 2, std::abs(ob.j2 - od.j2) / 2});
                  j4max = std::min({(op.j2 + oq.j2) / 2, (oh.j2 + og.j2) / 2, (ob.j2 + od.j2) / 2});

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jq, jp, J4, jb, jd, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jg, jh, J4, jb, jd, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, jg, jp, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jq, jp, J4, jb, jd, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jg, jh, J4, jb, jd, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, jg, jp, J4);

                        zpgqh -= phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, p, b, a, c) * Eta.TwoBody.GetTBME_J(J2, a, c, d, q) * Gamma.TwoBody.GetTBME_J(J3, d, g, h, b);
                      }
                    }
                  }

                  // exchange p <-> g
                  j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(od.j2 - oh.j2) / 2, std::abs(og.j2 - ob.j2) / 2});
                  j2max = std::min({(oa.j2 + oc.j2) / 2, (od.j2 + oh.j2) / 2, (og.j2 + ob.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - op.j2), std::abs(oq.j2 - ob.j2)) / 2;
                  j3max = std::min(od.j2 + op.j2, oq.j2 + ob.j2) / 2;

                  j4min = std::max({std::abs(og.j2 - oh.j2) / 2, std::abs(oq.j2 - op.j2) / 2, std::abs(ob.j2 - od.j2) / 2});
                  j4max = std::min({(og.j2 + oh.j2) / 2, (oq.j2 + op.j2) / 2, (ob.j2 + od.j2) / 2});

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jh, jg, J4, jb, jd, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jp, jq, J4, jb, jd, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, jp, jg, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jh, jg, J4, jb, jd, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jp, jq, J4, jb, jd, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, jp, jg, J4);

                        zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, g, b, a, c) * Eta.TwoBody.GetTBME_J(J2, a, c, d, h) * Gamma.TwoBody.GetTBME_J(J3, d, p, q, b);
                      }
                    }
                  }

                  // exchange p <-> g and q h
                  j2min = std::max({std::abs(oa.j2 - oc.j2) / 2, std::abs(od.j2 - oq.j2) / 2, std::abs(og.j2 - ob.j2) / 2});
                  j2max = std::min({(oa.j2 + oc.j2) / 2, (od.j2 + oq.j2) / 2, (og.j2 + ob.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - op.j2), std::abs(oh.j2 - ob.j2)) / 2;
                  j3max = std::min(od.j2 + op.j2, oh.j2 + ob.j2) / 2;

                  j4min = std::max({std::abs(og.j2 - oq.j2) / 2, std::abs(oh.j2 - op.j2) / 2, std::abs(ob.j2 - od.j2) / 2});
                  j4max = std::min({(og.j2 + oq.j2) / 2, (oh.j2 + op.j2) / 2, (ob.j2 + od.j2) / 2});

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jq, jg, J4, jb, jd, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(jp, jh, J4, jb, jd, J3);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, jp, jg, J4);

                        double sixj1 = Z.modelspace->GetSixJ(jq, jg, J4, jb, jd, J2);
                        double sixj2 = Z.modelspace->GetSixJ(jp, jh, J4, jb, jd, J3);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, jp, jg, J4);

                        zpgqh -= phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, g, b, a, c) * Eta.TwoBody.GetTBME_J(J2, a, c, d, q) * Gamma.TwoBody.GetTBME_J(J3, d, p, h, b);
                      }
                    }
                  }
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zpgqh);

        } // iket
      }   // ibra
    }     // J0 channel
    //std::cout << "diagram IIf  " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIIa
      //
      //   III(a)^J0_pgqh = P(p/g) * P(q/h) * \sum_{abcd J2 J3 J4 J5}
      //                   ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 ) ( 2 * J5 + 1 )
      //
      //                   { ja jb J5 } { J3 J0 J5 } { jq jd J5 } { J3 J0 J5 }
      //                   { jp jc J2 } { jp jc jg } { ja jb J4 } { jq jd jh }
      //
      //                   ( \bar{n_a} \bar{n_c} n_b + \bar{n_b} n_a n_c )
      //                   eta^J2_bpca eta^J3_gchd Gamma^J4_dabq
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_a * nbar_c * n_b + nbar_b * n_a * n_c);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  // direct term
                  int j2min = std::max(std::abs(ob.j2 - op.j2), std::abs(oc.j2 - oa.j2)) / 2;
                  int j2max = std::min(ob.j2 + op.j2, oc.j2 + oa.j2) / 2;

                  int j3min = std::max(std::abs(og.j2 - oc.j2), std::abs(oh.j2 - od.j2)) / 2;
                  int j3max = std::min(og.j2 + oc.j2, oh.j2 + od.j2) / 2;

                  int j4min = std::max(std::abs(od.j2 - oa.j2), std::abs(ob.j2 - oq.j2)) / 2;
                  int j4max = std::min(od.j2 + oa.j2, ob.j2 + oq.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oq.j2 - od.j2) / 2, std::abs(op.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                        int j5max = std::min({(oa.j2 + ob.j2) / 2, (oq.j2 + od.j2) / 2, (op.j2 + oc.j2) / 2, J0 + J3});

                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jb, J5, jp, jc, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(J3, J0, J5, jp, jc, jg);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jd, J5, ja, jb, J4);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J3, J0, J5, jq, jd, jh);

                          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J5, jp, jc, J2);
                          double sixj2 = Z.modelspace->GetSixJ(J3, J0, J5, jp, jc, jg);
                          double sixj3 = Z.modelspace->GetSixJ(jq, jd, J5, ja, jb, J4);
                          double sixj4 = Z.modelspace->GetSixJ(J3, J0, J5, jq, jd, jh);

                          zpgqh += occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, p, c, a) * Eta.TwoBody.GetTBME_J(J3, g, c, h, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, q);
                       if (ch==0 and ibra==2 and iket==4)
                       {
                          std::cout << " IIIa " << occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, p, c, a) * Eta.TwoBody.GetTBME_J(J3, g, c, h, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, q)
                                    << std::endl;
                       }
                        }
                      }
                    }
                  }

                  // exchange q <-> h
                  j2min = std::max(std::abs(ob.j2 - op.j2), std::abs(oc.j2 - oa.j2)) / 2;
                  j2max = std::min(ob.j2 + op.j2, oc.j2 + oa.j2) / 2;

                  j3min = std::max(std::abs(og.j2 - oc.j2), std::abs(oq.j2 - od.j2)) / 2;
                  j3max = std::min(og.j2 + oc.j2, oq.j2 + od.j2) / 2;

                  j4min = std::max(std::abs(od.j2 - oa.j2), std::abs(ob.j2 - oh.j2)) / 2;
                  j4max = std::min(od.j2 + oa.j2, ob.j2 + oh.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oh.j2 - od.j2) / 2, std::abs(op.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                        int j5max = std::min({(oa.j2 + ob.j2) / 2, (oh.j2 + od.j2) / 2, (op.j2 + oc.j2) / 2, J0 + J3});

                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jb, J5, jp, jc, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(J3, J0, J5, jp, jc, jg);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jd, J5, ja, jb, J4);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J3, J0, J5, jh, jd, jq);

                          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J5, jp, jc, J2);
                          double sixj2 = Z.modelspace->GetSixJ(J3, J0, J5, jp, jc, jg);
                          double sixj3 = Z.modelspace->GetSixJ(jh, jd, J5, ja, jb, J4);
                          double sixj4 = Z.modelspace->GetSixJ(J3, J0, J5, jh, jd, jq);

                          zpgqh += phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, p, c, a) * Eta.TwoBody.GetTBME_J(J3, g, c, q, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, h);
                       if (ch==0 and ibra==2 and iket==4)
                       {
                          std::cout << " IIIa " << phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, p, c, a) * Eta.TwoBody.GetTBME_J(J3, g, c, q, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, h)
                                    << std::endl;
                       }
                        }
                      }
                    }
                  }

                  // exchange p <-> g
                  j2min = std::max(std::abs(ob.j2 - og.j2), std::abs(oc.j2 - oa.j2)) / 2;
                  j2max = std::min(ob.j2 + og.j2, oc.j2 + oa.j2) / 2;

                  j3min = std::max(std::abs(op.j2 - oc.j2), std::abs(oh.j2 - od.j2)) / 2;
                  j3max = std::min(op.j2 + oc.j2, oh.j2 + od.j2) / 2;

                  j4min = std::max(std::abs(od.j2 - oa.j2), std::abs(ob.j2 - oq.j2)) / 2;
                  j4max = std::min(od.j2 + oa.j2, ob.j2 + oq.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oq.j2 - od.j2) / 2, std::abs(og.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                        int j5max = std::min({(oa.j2 + ob.j2) / 2, (oq.j2 + od.j2) / 2, (og.j2 + oc.j2) / 2, J0 + J3});

                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jb, J5, jg, jc, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(J3, J0, J5, jg, jc, jp);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jd, J5, ja, jb, J4);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J3, J0, J5, jq, jd, jh);

                          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J5, jg, jc, J2);
                          double sixj2 = Z.modelspace->GetSixJ(J3, J0, J5, jg, jc, jp);
                          double sixj3 = Z.modelspace->GetSixJ(jq, jd, J5, ja, jb, J4);
                          double sixj4 = Z.modelspace->GetSixJ(J3, J0, J5, jq, jd, jh);

                          zpgqh += phase_pg * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, g, c, a) * Eta.TwoBody.GetTBME_J(J3, p, c, h, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, q);
                       if (ch==0 and ibra==2 and iket==4)
                       {
                          std::cout << " IIIa " << phase_pg * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, g, c, a) * Eta.TwoBody.GetTBME_J(J3, p, c, h, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, q)
                                    << std::endl;
                       }
                        }
                      }
                    }
                  }

                  // exchange p <-> g and q <-> h
                  j2min = std::max(std::abs(ob.j2 - og.j2), std::abs(oc.j2 - oa.j2)) / 2;
                  j2max = std::min(ob.j2 + og.j2, oc.j2 + oa.j2) / 2;

                  j3min = std::max(std::abs(op.j2 - oc.j2), std::abs(oq.j2 - od.j2)) / 2;
                  j3max = std::min(op.j2 + oc.j2, oq.j2 + od.j2) / 2;

                  j4min = std::max(std::abs(od.j2 - oa.j2), std::abs(ob.j2 - oh.j2)) / 2;
                  j4max = std::min(od.j2 + oa.j2, ob.j2 + oh.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oh.j2 - od.j2) / 2, std::abs(og.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                        int j5max = std::min({(oa.j2 + ob.j2) / 2, (oh.j2 + od.j2) / 2, (og.j2 + oc.j2) / 2, J0 + J3});

                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jb, J5, jg, jc, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(J3, J0, J5, jg, jc, jp);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jd, J5, ja, jb, J4);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J3, J0, J5, jh, jd, jq);

                          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J5, jg, jc, J2);
                          double sixj2 = Z.modelspace->GetSixJ(J3, J0, J5, jg, jc, jp);
                          double sixj3 = Z.modelspace->GetSixJ(jh, jd, J5, ja, jb, J4);
                          double sixj4 = Z.modelspace->GetSixJ(J3, J0, J5, jh, jd, jq);

                          zpgqh += phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, g, c, a) * Eta.TwoBody.GetTBME_J(J3, p, c, q, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, h);
                       if (ch==0 and ibra==2 and iket==4)
                       {
                          std::cout << " IIIa " << phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, g, c, a) * Eta.TwoBody.GetTBME_J(J3, p, c, q, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, h)
                                    << std::endl;
                       }
                        }
                      }
                    }
                  }
                }
              }
            }
          } // a
          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIIa " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIIb
      //
      //   III(b)^J0_pgqh = P(p/g) * P(q/h) * \sum_{abcd J2 J3 J4 J5}
      //                   ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 ) ( 2 * J5 + 1 )
      //
      //                   { ja jb J5 }  { jd jp J5 } { J0 J3 J5 } { J0 J3 J5 }
      //                   { jc jq J2 }  { ja jb J4 } { jc jq jh } { jd jp jg }
      //
      //                   ( \bar{n_a} n_b n_c + \bar{n_b} \bar{n_c} n_a )
      //                   eta^J2_cbaq eta^J3_gdhc Gamma^J4_apdb
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_a * n_b * n_c + nbar_b * nbar_c * n_a);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  // direct term
                  int j2min = std::max(std::abs(ob.j2 - oc.j2), std::abs(oa.j2 - oq.j2)) / 2;
                  int j2max = std::min(ob.j2 + oc.j2, oa.j2 + oq.j2) / 2;

                  int j3min = std::max(std::abs(og.j2 - od.j2), std::abs(oh.j2 - oc.j2)) / 2;
                  int j3max = std::min(og.j2 + od.j2, oh.j2 + oc.j2) / 2;

                  int j4min = std::max(std::abs(oa.j2 - op.j2), std::abs(od.j2 - ob.j2)) / 2;
                  int j4max = std::min(oa.j2 + op.j2, od.j2 + ob.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(op.j2 - od.j2) / 2, std::abs(oq.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                        int j5max = std::min({(oa.j2 + ob.j2) / 2, (op.j2 + od.j2) / 2, (oq.j2 + oc.j2) / 2, J0 + J3});
                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jb, J5, jc, jq, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(J0, J3, J5, jc, jq, jh);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(jd, jp, J5, ja, jb, J4);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J3, J5, jd, jp, jg);

                          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J5, jc, jq, J2);
                          double sixj2 = Z.modelspace->GetSixJ(J0, J3, J5, jc, jq, jh);
                          double sixj3 = Z.modelspace->GetSixJ(jd, jp, J5, ja, jb, J4);
                          double sixj4 = Z.modelspace->GetSixJ(J0, J3, J5, jd, jp, jg);

                          zpgqh += occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, c, b, a, q) * Eta.TwoBody.GetTBME_J(J3, g, d, h, c) * Gamma.TwoBody.GetTBME_J(J4, a, p, d, b);
                        }
                      }
                    }
                  }

                  // exchange q <-> h
                  j2min = std::max(std::abs(ob.j2 - oc.j2), std::abs(oa.j2 - oh.j2)) / 2;
                  j2max = std::min(ob.j2 + oc.j2, oa.j2 + oh.j2) / 2;

                  j3min = std::max(std::abs(og.j2 - od.j2), std::abs(oq.j2 - oc.j2)) / 2;
                  j3max = std::min(og.j2 + od.j2, oq.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(oa.j2 - op.j2), std::abs(od.j2 - ob.j2)) / 2;
                  j4max = std::min(oa.j2 + op.j2, od.j2 + ob.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(op.j2 - od.j2) / 2, std::abs(oh.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                        int j5max = std::min({(oa.j2 + ob.j2) / 2, (op.j2 + od.j2) / 2, (oh.j2 + oc.j2) / 2, J0 + J3});
                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jb, J5, jc, jh, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(J0, J3, J5, jc, jh, jq);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(jd, jp, J5, ja, jb, J4);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J3, J5, jd, jp, jg);

                          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J5, jc, jh, J2);
                          double sixj2 = Z.modelspace->GetSixJ(J0, J3, J5, jc, jh, jq);
                          double sixj3 = Z.modelspace->GetSixJ(jd, jp, J5, ja, jb, J4);
                          double sixj4 = Z.modelspace->GetSixJ(J0, J3, J5, jd, jp, jg);

                          zpgqh += phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, c, b, a, h) * Eta.TwoBody.GetTBME_J(J3, g, d, q, c) * Gamma.TwoBody.GetTBME_J(J4, a, p, d, b);
                        }
                      }
                    }
                  }

                  // exchange p <-> g
                  j2min = std::max(std::abs(ob.j2 - oc.j2), std::abs(oa.j2 - oq.j2)) / 2;
                  j2max = std::min(ob.j2 + oc.j2, oa.j2 + oq.j2) / 2;

                  j3min = std::max(std::abs(op.j2 - od.j2), std::abs(oh.j2 - oc.j2)) / 2;
                  j3max = std::min(op.j2 + od.j2, oh.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(oa.j2 - og.j2), std::abs(od.j2 - ob.j2)) / 2;
                  j4max = std::min(oa.j2 + og.j2, od.j2 + ob.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(og.j2 - od.j2) / 2, std::abs(oq.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                        int j5max = std::min({(oa.j2 + ob.j2) / 2, (og.j2 + od.j2) / 2, (oq.j2 + oc.j2) / 2, J0 + J3});
                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jb, J5, jc, jq, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(J0, J3, J5, jc, jq, jh);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(jd, jg, J5, ja, jb, J4);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J3, J5, jd, jg, jp);

                          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J5, jc, jq, J2);
                          double sixj2 = Z.modelspace->GetSixJ(J0, J3, J5, jc, jq, jh);
                          double sixj3 = Z.modelspace->GetSixJ(jd, jg, J5, ja, jb, J4);
                          double sixj4 = Z.modelspace->GetSixJ(J0, J3, J5, jd, jg, jp);

                          zpgqh += phase_pg * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, c, b, a, q) * Eta.TwoBody.GetTBME_J(J3, p, d, h, c) * Gamma.TwoBody.GetTBME_J(J4, a, g, d, b);
                        }
                      }
                    }
                  }

                  // exchange p <-> g and q <-> h
                  j2min = std::max(std::abs(ob.j2 - oc.j2), std::abs(oa.j2 - oh.j2)) / 2;
                  j2max = std::min(ob.j2 + oc.j2, oa.j2 + oh.j2) / 2;

                  j3min = std::max(std::abs(op.j2 - od.j2), std::abs(oq.j2 - oc.j2)) / 2;
                  j3max = std::min(op.j2 + od.j2, oq.j2 + oc.j2) / 2;

                  j4min = std::max(std::abs(oa.j2 - og.j2), std::abs(od.j2 - ob.j2)) / 2;
                  j4max = std::min(oa.j2 + og.j2, od.j2 + ob.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        int j5min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(og.j2 - od.j2) / 2, std::abs(oh.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                        int j5max = std::min({(oa.j2 + ob.j2) / 2, (og.j2 + od.j2) / 2, (oh.j2 + oc.j2) / 2, J0 + J3});
                        for (int J5 = j5min; J5 <= j5max; J5++)
                        {
                          // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jb, J5, jc, jh, J2);
                          // double sixj2 = Z.modelspace->GetCachedSixJ(J0, J3, J5, jc, jh, jq);
                          // double sixj3 = Z.modelspace->GetCachedSixJ(jd, jg, J5, ja, jb, J4);
                          // double sixj4 = Z.modelspace->GetCachedSixJ(J0, J3, J5, jd, jg, jp);

                          double sixj1 = Z.modelspace->GetSixJ(ja, jb, J5, jc, jh, J2);
                          double sixj2 = Z.modelspace->GetSixJ(J0, J3, J5, jc, jh, jq);
                          double sixj3 = Z.modelspace->GetSixJ(jd, jg, J5, ja, jb, J4);
                          double sixj4 = Z.modelspace->GetSixJ(J0, J3, J5, jd, jg, jp);

                          zpgqh += phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, c, b, a, h) * Eta.TwoBody.GetTBME_J(J3, p, d, q, c) * Gamma.TwoBody.GetTBME_J(J4, a, g, d, b);
                        }
                      }
                    }
                  }
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIIb " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIIc
      //
      //   III(c)^J0_pgqh = - P(q/h) * \sum_{abcd J2 J3 J4} ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                   { jq jb J4 } { J3 J0 J4 } { J4 J3 J0 }
      //                   { jc ja J2 } { jc ja jd } { jh jq Jb }
      //
      //                   ( \bar{n_a} n_b n_c + \bar{n_b} \bar{n_c} n_a )
      //                   eta^J2_bcaq eta^J0_pgcd Gamma^J3_dahb
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();

      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;

          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_a * n_b * n_c + nbar_b * nbar_c * n_a);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  int j2min = std::max(std::abs(ob.j2 - oc.j2), std::abs(oa.j2 - oq.j2)) / 2;
                  int j2max = std::min(ob.j2 + oc.j2, oa.j2 + oq.j2) / 2;

                  int j3min = std::max(std::abs(od.j2 - oa.j2), std::abs(oh.j2 - ob.j2)) / 2;
                  int j3max = std::min(od.j2 + oa.j2, oh.j2 + ob.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(oq.j2 - ob.j2) / 2, std::abs(oc.j2 - oa.j2) / 2, std::abs(J0 - J3)});
                      int j4max = std::min({(oq.j2 + ob.j2) / 2, (oc.j2 + oa.j2) / 2, J0 + J3});
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jq, jb, J4, jc, ja, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J3, J0, J4, jc, ja, jd);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(J4, J3, J0, jh, jq, jb);

                        double sixj1 = Z.modelspace->GetSixJ(jq, jb, J4, jc, ja, J2);
                        double sixj2 = AngMom::SixJ(J3, J0, J4, jc, ja, jd);
                        double sixj3 = AngMom::SixJ(J4, J3, J0, jh, jq, jb);

                        zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, b, c, a, q) * Eta.TwoBody.GetTBME_J(J0, p, g, c, d) * Gamma.TwoBody.GetTBME_J(J3, d, a, h, b);
                      }
                    }
                  }

                  // exchanging q and h
                  j2min = std::max(std::abs(ob.j2 - oc.j2), std::abs(oa.j2 - oh.j2)) / 2;
                  j2max = std::min(ob.j2 + oc.j2, oa.j2 + oh.j2) / 2;

                  j3min = std::max(std::abs(od.j2 - oa.j2), std::abs(oq.j2 - ob.j2)) / 2;
                  j3max = std::min(od.j2 + oa.j2, oq.j2 + ob.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(oh.j2 - ob.j2) / 2, std::abs(oc.j2 - oa.j2) / 2, std::abs(J0 - J3)});
                      int j4max = std::min({(oh.j2 + ob.j2) / 2, (oc.j2 + oa.j2) / 2, J0 + J3});
                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(jh, jb, J4, jc, ja, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J3, J0, J4, jc, ja, jd);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(J4, J3, J0, jq, jh, jb);

                        double sixj1 = Z.modelspace->GetSixJ(jh, jb, J4, jc, ja, J2);
                        double sixj2 = AngMom::SixJ(J3, J0, J4, jc, ja, jd);
                        double sixj3 = AngMom::SixJ(J4, J3, J0, jq, jh, jb);

                        zpgqh -= phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, b, c, a, h) * Eta.TwoBody.GetTBME_J(J0, p, g, c, d) * Gamma.TwoBody.GetTBME_J(J3, d, a, q, b);
                      }
                    }
                  }

                  // ************************
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIIc " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIId
      //
      //   III(d)^J0_pgqh = - P(p/g) *  \sum_{abcd J2 J3 J4}
      //                      ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                   { ja jp J4 } { J0 J3 J4 } { J4 J3 J0 }
      //                   { jb jc J2 } { jb jc jd } { jg jp ja }
      //
      //                   ( \bar{n_a} \bar{n_c} n_b + \bar{n_b} n_a n_c )
      //                   eta^J2_bpac eta^J0_cdqh Gamma^J3_gadb
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;

          double zpgqh = 0.;

          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_a * nbar_c * n_b + nbar_b * n_a * n_c);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  int j2min = std::max(std::abs(ob.j2 - op.j2), std::abs(oa.j2 - oc.j2)) / 2;
                  int j2max = std::min(ob.j2 + op.j2, oa.j2 + oc.j2) / 2;

                  int j3min = std::max(std::abs(og.j2 - oa.j2), std::abs(od.j2 - ob.j2)) / 2;
                  int j3max = std::min(og.j2 + oa.j2, od.j2 + ob.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(oa.j2 - op.j2) / 2, std::abs(ob.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                      int j4max = std::min({(oa.j2 + op.j2) / 2, (ob.j2 + oc.j2) / 2, J0 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jp, J4, jb, jc, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J0, J3, J4, jb, jc, jd);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(J4, J3, J0, jg, jp, ja);

                        double sixj1 = Z.modelspace->GetSixJ(ja, jp, J4, jb, jc, J2);
                        double sixj2 = AngMom::SixJ(J0, J3, J4, jb, jc, jd);
                        double sixj3 = AngMom::SixJ(J4, J3, J0, jg, jp, ja);

                        zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, b, p, a, c) * Eta.TwoBody.GetTBME_J(J0, c, d, q, h) * Gamma.TwoBody.GetTBME_J(J3, g, a, d, b);
                      }
                    }
                  }

                  // exchange p and g
                  j2min = std::max(std::abs(ob.j2 - og.j2), std::abs(oa.j2 - oc.j2)) / 2;
                  j2max = std::min(ob.j2 + og.j2, oa.j2 + oc.j2) / 2;

                  j3min = std::max(std::abs(op.j2 - oa.j2), std::abs(od.j2 - ob.j2)) / 2;
                  j3max = std::min(op.j2 + oa.j2, od.j2 + ob.j2) / 2;

                  // int j4min = std::max(std::abs(oa.j2 - op.j2), std::abs(ob.j2 - oc.j2)) / 2;
                  // int j4max = std::min(oa.j2 + op.j2, ob.j2 + oc.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(oa.j2 - og.j2) / 2, std::abs(ob.j2 - oc.j2) / 2, std::abs(J0 - J3)});
                      int j4max = std::min({(oa.j2 + og.j2) / 2, (ob.j2 + oc.j2) / 2, J0 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(ja, jg, J4, jb, jc, J2);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J0, J3, J4, jb, jc, jd);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(J4, J3, J0, jp, jg, ja);

                        double sixj1 = Z.modelspace->GetSixJ(ja, jg, J4, jb, jc, J2);
                        double sixj2 = AngMom::SixJ(J0, J3, J4, jb, jc, jd);
                        double sixj3 = AngMom::SixJ(J4, J3, J0, jp, jg, ja);

                        zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, b, g, a, c) * Eta.TwoBody.GetTBME_J(J0, c, d, q, h) * Gamma.TwoBody.GetTBME_J(J3, p, a, d, b);
                      }
                    }
                  }

                  // ****************************************
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIId " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIIe
      //
      //   III(e)^J0_pgqh = - 1/2 * P(p/g) * P(q/h) * \sum_{abcd J2 J3 J4}
      //                   ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                   { jh jq J0 } { J2 J3 J4 } { J2 J3 J4 }
      //                   { jg jp J4 } { jp jh jd } { jq jg Jc }
      //
      //                   ( \bar{n_d} n_a n_b + \bar{n_a} \bar{n_b} n_d )
      //                   eta^J2_abhd eta^J3_dpcq Gamma^J2_gcab
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;

          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_d * n_a * n_b + nbar_a * nbar_b * n_d);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  // direct term
                  int j2min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oh.j2 - od.j2) / 2, std::abs(og.j2 - oc.j2) / 2});
                  int j2max = std::min({(oa.j2 + ob.j2) / 2, (oh.j2 + od.j2) / 2, (og.j2 + oc.j2) / 2});

                  int j3min = std::max(std::abs(od.j2 - op.j2), std::abs(oc.j2 - oq.j2)) / 2;
                  int j3max = std::min(od.j2 + op.j2, oc.j2 + oq.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(op.j2 - oh.j2) / 2, std::abs(oq.j2 - og.j2) / 2, std::abs(J2 - J3)});
                      int j4max = std::min({(op.j2 + oh.j2) / 2, (oq.j2 + og.j2) / 2, J2 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jp, jh, jd);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jq, jg, jc);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, jg, jp, J4);

                        double sixj1 = Z.modelspace->GetSixJ(J2, J3, J4, jp, jh, jd);
                        double sixj2 = Z.modelspace->GetSixJ(J2, J3, J4, jq, jg, jc);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, jg, jp, J4);
                        zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, a, b, h, d) * Eta.TwoBody.GetTBME_J(J3, d, p, c, q) * Gamma.TwoBody.GetTBME_J(J2, g, c, a, b);
                      }
                    }
                  }

                  // exchange q <-> h
                  j2min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oq.j2 - od.j2) / 2, std::abs(og.j2 - oc.j2) / 2});
                  j2max = std::min({(oa.j2 + ob.j2) / 2, (oq.j2 + od.j2) / 2, (og.j2 + oc.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - op.j2), std::abs(oc.j2 - oh.j2)) / 2;
                  j3max = std::min(od.j2 + op.j2, oc.j2 + oh.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(op.j2 - oq.j2) / 2, std::abs(oh.j2 - og.j2) / 2, std::abs(J2 - J3)});
                      int j4max = std::min({(op.j2 + oq.j2) / 2, (oh.j2 + og.j2) / 2, J2 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jp, jq, jd);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jh, jg, jc);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, jg, jp, J4);

                        double sixj1 = Z.modelspace->GetSixJ(J2, J3, J4, jp, jq, jd);
                        double sixj2 = Z.modelspace->GetSixJ(J2, J3, J4, jh, jg, jc);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, jg, jp, J4);

                        zpgqh -= phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, a, b, q, d) * Eta.TwoBody.GetTBME_J(J3, d, p, c, h) * Gamma.TwoBody.GetTBME_J(J2, g, c, a, b);
                      }
                    }
                  }

                  // exchange p <-> g
                  j2min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oh.j2 - od.j2) / 2, std::abs(op.j2 - oc.j2) / 2});
                  j2max = std::min({(oa.j2 + ob.j2) / 2, (oh.j2 + od.j2) / 2, (op.j2 + oc.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - og.j2), std::abs(oc.j2 - oq.j2)) / 2;
                  j3max = std::min(od.j2 + og.j2, oc.j2 + oq.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(og.j2 - oh.j2) / 2, std::abs(oq.j2 - op.j2) / 2, std::abs(J2 - J3)});
                      int j4max = std::min({(og.j2 + oh.j2) / 2, (oq.j2 + op.j2) / 2, J2 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jg, jh, jd);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jq, jp, jc);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, jp, jg, J4);

                        double sixj1 = Z.modelspace->GetSixJ(J2, J3, J4, jg, jh, jd);
                        double sixj2 = Z.modelspace->GetSixJ(J2, J3, J4, jq, jp, jc);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, jp, jg, J4);

                        zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, a, b, h, d) * Eta.TwoBody.GetTBME_J(J3, d, g, c, q) * Gamma.TwoBody.GetTBME_J(J2, p, c, a, b);
                      }
                    }
                  }

                  // exchange p <-> g and q <-> h
                  j2min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oq.j2 - od.j2) / 2, std::abs(op.j2 - oc.j2) / 2});
                  j2max = std::min({(oa.j2 + ob.j2) / 2, (oq.j2 + od.j2) / 2, (op.j2 + oc.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - og.j2), std::abs(oc.j2 - oh.j2)) / 2;
                  j3max = std::min(od.j2 + og.j2, oc.j2 + oh.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(og.j2 - oq.j2) / 2, std::abs(oh.j2 - op.j2) / 2, std::abs(J2 - J3)});
                      int j4max = std::min({(og.j2 + oq.j2) / 2, (oh.j2 + op.j2) / 2, J2 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jg, jq, jd);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jh, jp, jc);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, jp, jg, J4);

                        double sixj1 = Z.modelspace->GetSixJ(J2, J3, J4, jg, jq, jd);
                        double sixj2 = Z.modelspace->GetSixJ(J2, J3, J4, jh, jp, jc);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, jp, jg, J4);

                        zpgqh -= phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, a, b, q, d) * Eta.TwoBody.GetTBME_J(J3, d, g, c, h) * Gamma.TwoBody.GetTBME_J(J2, p, c, a, b);
                      }
                    }
                  }
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIIe " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

      // ####################################################################################
      //   diagram IIIf
      //
      //   III(f)^J0_pgqh = - 1/2 * P(p/g) * P(q/h) * \sum_{abcd J2 J3 J4}
      //                    ( 2 * J2 + 1 ) ( 2 * J3 + 1 ) ( 2 * J4 + 1 )
      //
      //                   { jh jq J0 } { J2 J3 J4 } { J2 J3 J4 }
      //                   { jg jp J4 } { jp jh jd } { jq jg Jc }
      //
      //                   ( \bar{n_c} n_a n_b + \bar{n_a} \bar{n_b} n_c )
      //                   eta^J2_gcab eta^J3_dpcq Gamma^J2_abhd
      // ####################################################################################
#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        for (int iket = ibra; iket < nKets; ++iket)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;
          Orbit &oq = *(ket.op);
          Orbit &oh = *(ket.oq);
          double jq = oq.j2 * 0.5;
          double jh = oh.j2 * 0.5;
          int phase_qh = ket.Phase(J0);

          double zpgqh = 0.;
          for (auto &a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double n_a = oa.occ;
            double nbar_a = 1.0 - n_a;
            double ja = oa.j2 * 0.5;

            for (auto &b : Z.modelspace->all_orbits)
            {
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double n_b = ob.occ;
              double nbar_b = 1.0 - n_b;
              double jb = ob.j2 * 0.5;

              for (auto &c : Z.modelspace->all_orbits)
              {
                Orbit &oc = Z.modelspace->GetOrbit(c);
                double n_c = oc.occ;
                double nbar_c = 1.0 - n_c;
                double jc = oc.j2 * 0.5;

                for (auto &d : Z.modelspace->all_orbits)
                {
                  Orbit &od = Z.modelspace->GetOrbit(d);
                  double n_d = od.occ;
                  double nbar_d = 1.0 - n_d;
                  double jd = od.j2 * 0.5;

                  double occfactor = (nbar_c * n_a * n_b + nbar_a * nbar_b * n_c);
                  if (fabs(occfactor) < 1.e-7)
                    continue;
                  // direct term
                  int j2min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oh.j2 - od.j2) / 2, std::abs(og.j2 - oc.j2) / 2});
                  int j2max = std::min({(oa.j2 + ob.j2) / 2, (oh.j2 + od.j2) / 2, (og.j2 + oc.j2) / 2});

                  int j3min = std::max(std::abs(od.j2 - op.j2), std::abs(oc.j2 - oq.j2)) / 2;
                  int j3max = std::min(od.j2 + op.j2, oc.j2 + oq.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(op.j2 - oh.j2) / 2, std::abs(oq.j2 - og.j2) / 2, std::abs(J2 - J3)});
                      int j4max = std::min({(op.j2 + oh.j2) / 2, (oq.j2 + og.j2) / 2, J2 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jq, jg, jc);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jp, jh, jd);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, jg, jp, J4);

                        double sixj1 = Z.modelspace->GetSixJ(J2, J3, J4, jq, jg, jc);
                        double sixj2 = Z.modelspace->GetSixJ(J2, J3, J4, jp, jh, jd);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, jg, jp, J4);

                        zpgqh -= occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, g, c, a, b) * Eta.TwoBody.GetTBME_J(J3, d, p, c, q) * Gamma.TwoBody.GetTBME_J(J2, a, b, h, d);
                      }
                    }
                  }

                  // exchange q <-> h
                  j2min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oq.j2 - od.j2) / 2, std::abs(og.j2 - oc.j2) / 2});
                  j2max = std::min({(oa.j2 + ob.j2) / 2, (oq.j2 + od.j2) / 2, (og.j2 + oc.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - op.j2), std::abs(oc.j2 - oh.j2)) / 2;
                  j3max = std::min(od.j2 + op.j2, oc.j2 + oh.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(op.j2 - oq.j2) / 2, std::abs(oh.j2 - og.j2) / 2, std::abs(J2 - J3)});
                      int j4max = std::min({(op.j2 + oq.j2) / 2, (oh.j2 + og.j2) / 2, J2 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jh, jg, jc);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jp, jq, jd);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, jg, jp, J4);

                        double sixj1 = Z.modelspace->GetSixJ(J2, J3, J4, jh, jg, jc);
                        double sixj2 = Z.modelspace->GetSixJ(J2, J3, J4, jp, jq, jd);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, jg, jp, J4);

                        zpgqh -= phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, g, c, a, b) * Eta.TwoBody.GetTBME_J(J3, d, p, c, h) * Gamma.TwoBody.GetTBME_J(J2, a, b, q, d);
                      }
                    }
                  }

                  // exchange p <-> g
                  j2min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oh.j2 - od.j2) / 2, std::abs(op.j2 - oc.j2) / 2});
                  j2max = std::min({(oa.j2 + ob.j2) / 2, (oh.j2 + od.j2) / 2, (op.j2 + oc.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - og.j2), std::abs(oc.j2 - oq.j2)) / 2;
                  j3max = std::min(od.j2 + og.j2, oc.j2 + oq.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(op.j2 - oq.j2) / 2, std::abs(oh.j2 - og.j2) / 2, std::abs(J2 - J3)});
                      int j4max = std::min({(op.j2 + oq.j2) / 2, (oh.j2 + og.j2) / 2, J2 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jq, jp, jc);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jg, jh, jd);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jh, jq, J0, jp, jg, J4);

                        double sixj1 = Z.modelspace->GetSixJ(J2, J3, J4, jq, jp, jc);
                        double sixj2 = Z.modelspace->GetSixJ(J2, J3, J4, jg, jh, jd);
                        double sixj3 = Z.modelspace->GetSixJ(jh, jq, J0, jp, jg, J4);

                        zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, p, c, a, b) * Eta.TwoBody.GetTBME_J(J3, d, g, c, q) * Gamma.TwoBody.GetTBME_J(J2, a, b, h, d);
                      }
                    }
                  }

                  // exchange p <-> g and q <-> h
                  j2min = std::max({std::abs(oa.j2 - ob.j2) / 2, std::abs(oq.j2 - od.j2) / 2, std::abs(op.j2 - oc.j2) / 2});
                  j2max = std::min({(oa.j2 + ob.j2) / 2, (oq.j2 + od.j2) / 2, (op.j2 + oc.j2) / 2});

                  j3min = std::max(std::abs(od.j2 - og.j2), std::abs(oc.j2 - oh.j2)) / 2;
                  j3max = std::min(od.j2 + og.j2, oc.j2 + oh.j2) / 2;

                  for (int J2 = j2min; J2 <= j2max; J2++)
                  {
                    for (int J3 = j3min; J3 <= j3max; J3++)
                    {
                      int j4min = std::max({std::abs(op.j2 - oh.j2) / 2, std::abs(oq.j2 - og.j2) / 2, std::abs(J2 - J3)});
                      int j4max = std::min({(op.j2 + oh.j2) / 2, (oq.j2 + og.j2) / 2, J2 + J3});

                      for (int J4 = j4min; J4 <= j4max; J4++)
                      {
                        // double sixj1 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jh, jp, jc);
                        // double sixj2 = Z.modelspace->GetCachedSixJ(J2, J3, J4, jg, jq, jd);
                        // double sixj3 = Z.modelspace->GetCachedSixJ(jq, jh, J0, jp, jg, J4);

                        double sixj1 = Z.modelspace->GetSixJ(J2, J3, J4, jh, jp, jc);
                        double sixj2 = Z.modelspace->GetSixJ(J2, J3, J4, jg, jq, jd);
                        double sixj3 = Z.modelspace->GetSixJ(jq, jh, J0, jp, jg, J4);

                        zpgqh -= phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, p, c, a, b) * Eta.TwoBody.GetTBME_J(J3, d, g, c, h) * Gamma.TwoBody.GetTBME_J(J2, a, b, q, d);
                      }
                    }
                  }
                  // **************************************************************
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIIf " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

    return;
  }

  void comm223_231(const Operator &Eta, const Operator &Gamma, Operator &Z)
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
    } // d

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
    } // for p
    // std::cout << "diagram I  " << Z.OneBodyNorm() << std::endl;
    // Z.EraseOneBody();

    // *********************************************************************************** //
    //                                  Diagram II                                         //
    // *********************************************************************************** //

    Z.profiler.timer["231_diagram_I"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime(); // timer

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
    } // for p
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
    } // for p
      // std::cout<< "diagram IId " << Z.OneBodyNorm() << std::endl;
      // Z.EraseOneBody();

    // *********************************************************************************** //
    //                                  Diagram III                                        //
    // *********************************************************************************** //

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
    } // d

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
      // std::cout<< "diagram IIIa and IIIb " << Z.OneBodyNorm() << std::endl;
      // Z.EraseOneBody();

    Z.profiler.timer[__func__] += omp_get_wtime() - t_internal;
    return;
  }

  void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z)
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
        } // for a
        CHI_I(p, q) = chi_pq;
        CHI_II(p, q) = chiY_pq;
      } // for q
    } // for p

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
            zpqrs += CHI_I(p, b) * Gamma.TwoBody.GetTBME_J(J, J, b, q, r, s);
            zpqrs += hZ * CHI_II(b, p) * Eta.TwoBody.GetTBME_J(J, J, b, q, r, s); // tricky minus sign.
          } // for a
          for (auto b : Z.OneBodyChannels.at({oq.l, oq.j2, oq.tz2}))
          {
            zpqrs += CHI_I(q, b) * Gamma.TwoBody.GetTBME_J(J, J, p, b, r, s);
            zpqrs += hZ * CHI_II(b, q) * Eta.TwoBody.GetTBME_J(J, J, p, b, r, s); // tricky minus sign.
          } // for a
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
      } // for ibra

    } // for itmat

    Z.profiler.timer[std::string(__func__) + " Diagram I and Diagram IV"] += omp_get_wtime() - t_type;
    t_type = omp_get_wtime();

    Z.profiler.timer[std::string(__func__) + " Diagram I and Diagram IV"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

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
      } // ibra
    } // J0 channel

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

    // *********************************************************************************** //
    //                                Diagram III                                          //
    // *********************************************************************************** //

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
      } // ibra

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

    // Timer
    Z.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return;
  }

} // namespace ReferenceImplementations
////////////////////////////////////////////////////////////////////
