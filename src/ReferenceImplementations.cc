
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
            }   // c
          }     // b
        }       // a

      } // j
    }   // i
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
      }   // ibra
    }     // ch
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

    int nch = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ch = 0; ch < nch; ch++)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nkets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nkets; ibra++)
      {
        Ket &bra = tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;

        for (int iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;

          double zijkl = 0;
          for (int iket_ab = 0; iket_ab < nkets; iket_ab++)
          {
            Ket &ket_ab = tbc.GetKet(iket_ab);
            size_t a = ket_ab.p;
            size_t b = ket_ab.q;
            Orbit &oa = Z.modelspace->GetOrbit(a);
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double xijab = X2.GetTBME(ch, ch, bra, ket_ab);
            double yijab = Y2.GetTBME(ch, ch, bra, ket_ab);
            double xabkl = X2.GetTBME(ch, ch, ket_ab, ket);
            double yabkl = Y2.GetTBME(ch, ch, ket_ab, ket);
            double flip_factor = (a == b) ? 1 : 2; // looping over kets only gets a<=b. So we need a factor of 2 for the other ordering.
            zijkl += 1. / 2 * flip_factor * ((1 - oa.occ) * (1 - ob.occ) - oa.occ * ob.occ) * (xijab * yabkl - yijab * xabkl);
          } // iket_ab

          // Need to normalize here, because AddToTBME expects a normalized TBME.
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch, ch, ibra, iket, zijkl);

        } // iket
      }   // ibra
    }     // ch
  }       // comm222_pp_hhss

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
            }   // a

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
      }   // ibra
    }     // ch

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

                if ((oa.l + ob.l + oc.l + od.l + oe.l + of.l) % 2 > 0)
                  continue;
                if ((oa.tz2 + ob.tz2 + oc.tz2) != (od.tz2 + oe.tz2 + of.tz2))
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
                  }   // J2
                }     // J1
              }       // f
            }         // e
          }           // d
        }             // c
      }               // b
    }                 // a
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
                  if ((oa.l + ob.l + oi.l + oc.l + od.l + oe.l + Z.parity) % 2 > 0)
                    continue;
                  if (std::abs(oa.tz2 + ob.tz2 + oi.tz2 - oc.tz2 - od.tz2 - oe.tz2) != Z.GetTRank())
                    continue;
                  //                if ( std::abs( oa.occ * ob.occ *(1-oc.occ)*(1-od.occ)*(1-oe.occ) )<1e-8 ) continue;

                  double cde_symmetry_factor = 1;                                                         // why this???
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
                        zij += 1. / 12 * cde_symmetry_factor * occupation_factor * (twoJ + 1.) / (oi.j2 + 1) * (xabicde * ycdeabj - yabicde * xcdeabj);
                      }
                    } // Jj
                  }   // J1
                }     // e
              }       // d
            }         // c
          }           // b
        }             // a
        Z1(i, j) += zij;
      } // j
    }   // i
  }     // comm331ss

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
                }   // J1
              }     // d
            }       // c
          }         // b
        }           // a
        Z1(i, j) += zij;
      } // j
    }   // i
  }     // comm231ss

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
        for (int iket = ibra; iket < nkets; iket++)
        {
          Ket &ket = tbc.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;
          double zijkl = 0;
          for (size_t a : Z.modelspace->all_orbits)
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            //           if (a != 0) continue; // TODO: REMOVE THIS!!!!
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

                //                if (i==0 and j==10 and k==10 and l==10)
                //                {
                //                    std::cout << "a b " << a << " " << b << "  yab = " << yab << " " << "xijbkla = " << xijbkla << " J1 twoJ = " << J1 << " " << twoJ << "   z = " <<zijkl << std::endl;
                //                }
              } // twoJ

            } // b
          }   // a
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          Z2.AddToTBME(ch2, ch2, ibra, iket, zijkl);
        } // iket
      }   // ibra
    }     // ch2

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
            }   // for ch2
          }     // for c

          // normalize the tbme
          zijkl *= -1.0 / sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
        } // for iket
      }   // for ibra
    }     // for ch

  } // comm232ss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  /// Expression:    ZJ1ijkl = 1/6 sum_abcd sum_J2J3  (na nb nc(1-nd) + (1-na)(1-nb)(1-nc)nd) (2J3+1)/(2J1+1)
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

                  if (std::abs(oa.occ * ob.occ * oc.occ * (1 - od.occ) + (1 - oa.occ) * (1 - ob.occ) * (1 - oc.occ) * od.occ) < 1e-7)
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
                  }   // J2
                }     // d
              }       // c
            }         // b
          }           // a
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          Z2.AddToTBME(ch2, ch2, ibra, iket, zijkl);
        } // iket
      }   // ibra
    }     // ch2

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
                  }   // for twoJp

                } // for iket_cd
              }   // for ch_cd
            }     // for iket_ab
          }       // for ch_ab
          // make it a normalized TBME
          zijkl /= sqrt((1. + bra.delta_pq()) * (1. + ket.delta_pq()));
          // the AddToTBME routine automatically takes care of the hermitian conjugate as well
          Z2.AddToTBME(ch, ch, ibra, iket, zijkl);
        } // for iket
      }   // for ibra
    }     // for ch
  }

  /*

  void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z )
  {
     auto& X3 = X.ThreeBody;
     auto& Y3 = Y.ThreeBody;
     auto& Z2 = Z.TwoBody;

     int nch2 = Z.modelspace->GetNumberTwoBodyChannels();
     #pragma omp parallel for schedule(dynamic,1)
     for (int ch2=0; ch2<nch2; ch2++)
     {
       TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch2);
       int J1 = tbc.J;
       int nkets = tbc.GetNumberKets();
       for (int ibra=0; ibra<nkets; ibra++)
       {
         Ket& bra = tbc.GetKet(ibra);
         size_t i=bra.p;
         size_t j=bra.q;
         Orbit& oi = Z.modelspace->GetOrbit(i);
         Orbit& oj = Z.modelspace->GetOrbit(j);
         for (int iket=ibra; iket<nkets; iket++)
         {
           Ket& ket = tbc.GetKet(iket);
           size_t k=ket.p;
           size_t l=ket.q;

           double zijkl = 0;

           // loop over permutations  ijkl -> 1234

           double zbar_1432 = 0;





           zijkl /= sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
           Z2.AddToTBME(ch2,ch2, ibra,iket, zijkl);
         }//iket
       }//ibra
     }//ch2

  }
  */

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
    //  #pragma omp parallel for schedule(dynamic,1)
    for (size_t ch3 = 0; ch3 < nch3; ch3++)
    {
      auto &Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
      int twoJ = Tbc.twoJ;
      size_t nkets = Tbc.GetNumberKets();
      for (size_t ibra = 0; ibra < nkets; ibra++)
      {
        Ket3 &bra = Tbc.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        Orbit &ok = Z.modelspace->GetOrbit(k);
        int Jij = bra.Jpq;
        for (size_t iket = ibra; iket < nkets; iket++)
        {
          Ket3 &ket = Tbc.GetKet(iket);
          size_t l = ket.p;
          size_t m = ket.q;
          size_t n = ket.r;
          Orbit &ol = Z.modelspace->GetOrbit(l);
          Orbit &om = Z.modelspace->GetOrbit(m);
          Orbit &on = Z.modelspace->GetOrbit(n);
          int Jlm = ket.Jpq;

          double zsum = 0;
          // First, connect on the bra side
          //        for (auto a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
          for (auto a : X.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
          {
            zsum += X1(i, a) * Y3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
            zsum -= Y1(i, a) * X3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
          }
          //        for (auto a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
          for (auto a : X.GetOneBodyChannel(oj.l, oj.j2, oj.tz2))
          {
            zsum += X1(j, a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
            zsum -= Y1(j, a) * X3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
          }
          for (auto a : X.GetOneBodyChannel(ok.l, ok.j2, ok.tz2))
          {
            zsum += X1(k, a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
            zsum -= Y1(k, a) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
          }

          // Now connect on the ket side
          for (auto a : X.GetOneBodyChannel(ol.l, ol.j2, ol.tz2))
          {
            zsum -= X1(a, l) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
            zsum += Y1(a, l) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
          }
          for (auto a : X.GetOneBodyChannel(om.l, om.j2, om.tz2))
          {
            zsum -= X1(a, m) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
            zsum += Y1(a, m) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
          }
          for (auto a : X.GetOneBodyChannel(on.l, on.j2, on.tz2))
          {
            zsum -= X1(a, n) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
            zsum += Y1(a, n) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
          }

          //        Z3.AddToME_pn(Jij, Jlm, twoJ, i,j,k,l,m,n, zsum );
          Z3.AddToME_pn_ch(ch3, ch3, ibra, iket, zsum);

        } // for iket
      }   // for ibra
    }     // for ch3

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

              int parity_a = (o1.l + o2.l + o6.l) % 2; // Need to fix this if we want to treat parity changing operators.
              for (int tz2a : {-1, 1})
              {
                int dTz = o1.tz2 + o2.tz2 - o6.tz2 - tz2a;
                if (std::abs(dTz) != X.GetTRank() and std::abs(dTz) != Y.GetTRank())
                  continue;

                for (int J2p = J2p_min; J2p <= J2p_max; J2p++)
                {

                  double rec_lmn = Z3.RecouplingCoefficient(perm_lmn, jl, jm, jn, J2p, J2, twoJ);
                  rec_lmn *= Z3.PermutationPhase(perm_lmn); // do we get a fermionic minus sign?

                  int j2a_min = std::max(std::abs(o6.j2 - 2 * J1p), std::abs(o3.j2 - 2 * J2p));
                  int j2a_max = std::min(o6.j2 + 2 * J1p, o3.j2 + 2 * J2p);

                  for (int j2a = j2a_min; j2a <= j2a_max; j2a += 2)
                  {
                    double ja = 0.5 * j2a;
                    int la = (j2a - 1) / 2 + ((j2a - 1) / 2 + parity_a) % 2;
                    if (la > Z.modelspace->GetEmax())
                      continue;

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

                    for (size_t a : Z.GetOneBodyChannel(la, j2a, tz2a))
                    {
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
                                                               //    }// for ibra
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
            }   // a
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
            }   // a
          }
        } // for perm_lmn

        //          }//b
        //        }//a

        Z3.AddToME_pn_ch(ch3bra, ch3ket, ibra, iket, zijklmn); // this needs to be modified for beta decay
      }                                                        // for iket
                                                               //    }// for ibra
    }                                                          // for ch3

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

                      }                                        // Ja6
                    }                                          // twoJp
                  }                                            // b
                }                                              // a
              }                                                // J2p
            }                                                  // perm_lmn
          }                                                    // J1p
        }                                                      // perm_ijk
        Z3.AddToME_pn_ch(ch3bra, ch3ket, ibra, iket, zijklmn); // this needs to be modified for beta decay
      }                                                        // for iket
                                                               //    }//ibra
    }                                                          // ch
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
                          }   // for twoJx
                              //                        }
                        }     // Z1 block

                      } // for c
                    }   // for iket_ab
                  }     // for ch2

                } // for J2p
              }   // for perm_lmn
            }     // for J1p
          }       // for perm_ijk

          Z3.AddToME_pn_ch(ch3, ch3, ibra, iket, z_ijklmn);
        } // for iket
      }   // for ibra
    }     // for ch3

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
                }   // for a

                double ninej = AngMom::NineJ(o1.j2 / 2., o4.j2 / 2., J3, o2.j2 / 2., o3.j2 / 2., J4, J1, J2, lambda);
                //                         zpqrs += phaseperm* sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) ) * nanb * ninej * AngMom::phase((o2.j2+o4.j2)/2 +J2+J4) * Xbar_psab * Ybar_abrq;
                zpqrs += phaseperm * sqrt((2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1)) * ninej * AngMom::phase((o2.j2 + o4.j2) / 2 + J2 + J4) * Zbar_1432;

                if (ch_bra == 1)
                {
                  std::cout << __func__ << " " << __LINE__ << " channel J = " << tbc_bra.J << "  J3,J4 " << J3 << " " << J4 << " iperm " << iperm << " : "
                            << phaseperm << " * " << sqrt((2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1)) << " * " << ninej << " * " << AngMom::phase((o2.j2 + o4.j2) / 2 + J2 + J4) << " * " << Zbar_1432 << " (<-Zbar_1432, zpqrs ->)  " << zpqrs << std::endl;
                }
              } // for J4
            }   // for J3
          }     // for iperm

          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;

          Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zpqrs);
        } // for iket
      }   // for ibra
          //    if ( ch_bra==3) std::cout << __func__ << "   ch=3 Jpt " << tbc_bra.J << " " << tbc_bra.parity << " " << tbc_bra.Tz << "   " << tbc_bra.GetKet(0).p << " " << tbc_bra.GetKet(0).q << " Z2 = " << std::endl << iter.second << std::endl;
    }     // for iter
  }

  ////// Start double nested commutators

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
            }   // for i
          }     // for b
        }       // for a
        Z.OneBody(p, q) += zpq;
      } // for q
    }   // for p
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
              }   // for i
            }     // for b
          }       // for a
          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      }   // for ibra

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
              }   // for i
            }     // for a
            zpqrs += xpqks * chi_rk + xpqrk * chi_sk + xkqrs * chi_pk + xpkrs * chi_qk;
          } // for k
          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      }   // for ibra

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
              }   // for i
            }     // for b
          }       // for a
          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      }   // for ibra

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
              }   // for i
            }     // for b
          }       // for a
          // normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;
          Z.TwoBody.AddToTBME(ch_bra, ch_ket, bra, ket, zpqrs);
        } // for iket
      }   // for ibra

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

              int J1min = std::max(std::abs(oa.j2 - ob.j2), std::abs(oi.j2 - oq.j2)) / 2;
              int J1max = std::min(oa.j2 + ob.j2, oi.j2 + oq.j2) / 2;

              for (int J1 = J1min; J1 <= J1max; J1++)
              {

                double xipab = X.TwoBody.GetTBME_J(J1, J1, i, p, a, b);
                double xabiq = X.TwoBody.GetTBME_J(J1, J1, a, b, i, q);
                double yabiq = Y.TwoBody.GetTBME_J(J1, J1, a, b, i, q);

                chi_pq += 0.5 * (2 * J1 + 1) / (oq.j2 + 1) * xipab * xabiq;
                chiY_pq += 0.5 * (2 * J1 + 1) / (oq.j2 + 1) * xipab * yabiq;
              }
            } // for b

          } // for i
        }   // for a
        CHI_XX(p, q) = chi_pq;
        CHI_XY(p, q) = chiY_pq;
      } // for q
    }   // for p

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
            zpqrs -= CHI_XY(p, b) * X.TwoBody.GetTBME_J(J, J, b, q, r, s) + CHI_XY(q, b) * X.TwoBody.GetTBME_J(J, J, p, b, r, s);
            zpqrs -= X.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_XY(b, r) + X.TwoBody.GetTBME_J(J, J, p, q, r, b) * CHI_XY(b, s);

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
      }   // for ibra

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
  void comm223_321_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {
    Z.modelspace->PreCalculateSixJ();
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

                  int J0min = std::max({std::abs(oa.j2 - oc.j2), std::abs(ob.j2 - od.j2), std::abs(ob.j2 - oe.j2)}) / 2;
                  int J0max = std::min({oa.j2 + oc.j2, ob.j2 + od.j2, ob.j2 + oe.j2}) / 2;

                  int J1min = std::max({std::abs(oe.j2 - op.j2), std::abs(od.j2 - oq.j2)}) / 2;
                  int J1max = std::min({oe.j2 + op.j2, od.j2 + oq.j2}) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    for (int J1 = J1min; J1 <= J1max; J1++)
                    {
                      zij += (2 * J0 + 1) * (2 * J1 + 1) / (od.j2 + 1.0) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, d, a, c) * Eta.TwoBody.GetTBME_J(J0, J0, a, c, b, e) * Gamma.TwoBody.GetTBME_J(J1, J1, e, p, d, q);
                    }
                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) += 0.5 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += 0.5 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p

    // std::cout << "diagram I  " << Z.OneBodyNorm() << std::endl;
    // Z.EraseOneBody();

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
          Z.OneBody(q, p) += zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
    // std::cout << "diagram IIa " << Z.OneBodyNorm() << std::endl;
    // Z.EraseOneBody();

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

                  int J0min = std::max({std::abs(oa.j2 - od.j2), std::abs(ob.j2 - oe.j2), std::abs(oc.j2 - op.j2)}) / 2;
                  int J0max = std::min({oa.j2 + od.j2, ob.j2 + oe.j2, op.j2 + oc.j2}) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    zij += (2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, e, a, d) * Eta.TwoBody.GetTBME_J(J0, J0, c, p, b, e) * Gamma.TwoBody.GetTBME_J(J0, J0, a, d, c, q);
                  } // J0
                }
              }
            }
          }
        }
        Z.OneBody(p, q) += 0.25 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += 0.25 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
    // std::cout << "diagram IIb " << Z.OneBodyNorm() << std::endl;
    // Z.EraseOneBody();

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
          Z.OneBody(q, p) += zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
    // std::cout << "diagram IIc " << Z.OneBodyNorm() << std::endl;
    // Z.EraseOneBody();

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
          Z.OneBody(q, p) -= 0.25 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
    // std::cout << "diagram IId " << Z.OneBodyNorm() << std::endl;
    // Z.EraseOneBody();

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
          Z.OneBody(q, p) += 0.5 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
    // std::cout << "diagram IIIa " << Z.OneBodyNorm() << std::endl;
    // Z.EraseOneBody();

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
          Z.OneBody(q, p) -= 0.5 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    }   // for p
    // std::cout << "diagram IIIb " << Z.OneBodyNorm() << std::endl;
    // Z.EraseOneBody();

    return;
  }



} // namespace ReferenceImplementations
////////////////////////////////////////////////////////////////////
