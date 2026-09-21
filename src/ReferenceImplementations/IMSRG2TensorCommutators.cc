
#include "ReferenceImplementations.hh"
#include "ModelSpace.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"


/// Straightforward implementation of J-coupled commutator expressions
/// without optimizations. This should be benchmarked against the
/// mscheme implementation and then left untouched.
namespace ReferenceImplementations
{
  /// Start of scalar-tensor commutators

  void comm111st(const Operator &X, const Operator &Y, Operator &Z)
  {

    double t_start = omp_get_wtime();
    for (auto i : Z.modelspace->all_orbits)
    {
      for (auto j : Z.modelspace->all_orbits)
      {
        double zij = 0;
        for (auto a : Z.modelspace->all_orbits)
        {
          zij += X.OneBody(i, a) * Y.OneBody(a, j) - Y.OneBody(i, a) * X.OneBody(a, j);
        } // for a
        Z.OneBody(i, j) += zij;
      } // for j
    } // for i
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  } // comm111st

  // Streightforward implementation of second term of equation (B2) in Parzuchowski et al PRC 96, 034324 (2017).
  // Note that equation (B4) of that paper has a typo. The Jhat should be Jhat^2.
  //
  void comm121st(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    Z.modelspace->PreCalculateNineJ();
    int lambda = Y.GetJRank();
    for (auto i : Z.modelspace->all_orbits)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double ji = 0.5 * oi.j2;
      for (auto j : Z.modelspace->all_orbits)
      {
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double jj = 0.5 * oj.j2;
        double zij = 0;
        for (auto a : Z.modelspace->all_orbits)
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);
          double jahat = sqrt(oa.j2 + 1);
          double ja = 0.5 * oa.j2;
          double na = oa.occ;
          for (auto b : Z.modelspace->all_orbits)
          {
            Orbit &ob = Z.modelspace->GetOrbit(b);
            double jb = 0.5 * ob.j2;
            double nb = ob.occ;
            double Xbar_ijab = 0;
            int J2min = std::max(std::abs(oi.j2 - ob.j2), std::abs(oa.j2 - oj.j2)) / 2;
            int J2max = std::min(oi.j2 + ob.j2, oj.j2 + oa.j2) / 2;
            for (int J2 = J2min; J2 <= J2max; J2++)
            {
              double sixj = AngMom::SixJ(ji, jj, lambda, ja, jb, J2);
              Xbar_ijab -= (2 * J2 + 1) * sixj * X.TwoBody.GetTBME_J(J2, J2, i, b, a, j);
            }
            double Ybar_ijab = 0;
            int J3min = std::abs(oi.j2 - ob.j2) / 2;
            int J3max = (oi.j2 + ob.j2) / 2;
            for (int J3 = J3min; J3 <= J3max; J3++)
            {
              int J4min = std::max(std::abs(oj.j2 - oa.j2), std::abs(2 * J3 - 2 * lambda)) / 2;
              int J4max = std::min(oj.j2 + oa.j2, 2 * J3 + 2 * lambda) / 2;
              for (int J4 = J4min; J4 <= J4max; J4++)
              {
                //                double ninej = Z.modelspace->GetNineJ( ji, jb, J3,  jj, ja, J4, lambda, 0, lambda);
                double ninej = AngMom::NineJ(ji, jb, J3, jj, ja, J4, lambda, 0, lambda);
                Ybar_ijab -= sqrt((2 * lambda + 1) * (2 * 0 + 1) * (2 * J3 + 1) * (2 * J4 + 1)) * AngMom::phase(jj + jb + 0 + J4) * ninej * Y.TwoBody.GetTBME_J(J3, J4, i, b, a, j);
              }
            }
            zij -= (na - nb) * (Xbar_ijab * Y.OneBody(a, b) - jahat * Ybar_ijab * X.OneBody(a, b));
          } // for b
        } // for a
        Z.OneBody(i, j) += zij;
      } // for j
    } // for i
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  } // comm121st

  // This seems to be working --SRS 12/17/2024.
  void comm122st(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    Z.modelspace->PreCalculateSixJ();
    int lambda = Y.GetJRank();
    int phase_lambda = Z.modelspace->phase(lambda); // (-1)^lambda
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
      double J1J2hat = sqrt((2 * J1 + 1) * (2 * J2 + 1));
      for (size_t ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t q = bra.q;
        Orbit &op = Z.modelspace->GetOrbit(p);
        Orbit &oq = Z.modelspace->GetOrbit(q);
        double jp = op.j2 * 0.5;
        double jq = oq.j2 * 0.5;
        for (size_t iket = 0; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          Orbit &oR = Z.modelspace->GetOrbit(r);
          Orbit &os = Z.modelspace->GetOrbit(s);
          double jr = oR.j2 * 0.5;
          double js = os.j2 * 0.5;
          double zpqrs = 0;
          if ((ch_bra == ch_ket) and (iket > ibra))
            continue;

          for (auto a : X.GetOneBodyChannel(op.l, op.j2, op.tz2))
          {
            zpqrs += X.OneBody(p, a) * Y.TwoBody.GetTBME_J(J1, J2, a, q, r, s);
          }
          for (auto a : X.GetOneBodyChannel(oq.l, oq.j2, oq.tz2))
          {
            zpqrs += X.OneBody(q, a) * Y.TwoBody.GetTBME_J(J1, J2, p, a, r, s);
          }
          for (auto a : X.GetOneBodyChannel(oR.l, oR.j2, oR.tz2))
          {
            zpqrs -= X.OneBody(a, r) * Y.TwoBody.GetTBME_J(J1, J2, p, q, a, s);
          }
          for (auto a : X.GetOneBodyChannel(os.l, os.j2, os.tz2))
          {
            zpqrs -= X.OneBody(a, s) * Y.TwoBody.GetTBME_J(J1, J2, p, q, r, a);
          }

          for (auto a : Y.GetOneBodyChannel(op.l, op.j2, op.tz2))
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double ja = oa.j2 * 0.5;
            double sixj = AngMom::SixJ(J2, J1, lambda, jp, ja, jq);
            zpqrs -= J1J2hat * phase_lambda * Z.modelspace->phase(jp + jq + J2) * sixj * Y.OneBody(p, a) * X.TwoBody.GetTBME_J(J2, J2, a, q, r, s);
          }
          for (auto a : Y.GetOneBodyChannel(oq.l, oq.j2, oq.tz2))
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double ja = oa.j2 * 0.5;
            double sixj = AngMom::SixJ(J2, J1, lambda, jq, ja, jp);
            zpqrs += J1J2hat * phase_lambda * Z.modelspace->phase(J1 - J2) * sixj * Y.OneBody(q, a) * X.TwoBody.GetTBME_J(J2, J2, a, p, r, s);
          }
          for (auto a : Y.GetOneBodyChannel(os.l, os.j2, os.tz2))
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double ja = oa.j2 * 0.5;
            double sixj = AngMom::SixJ(J1, J2, lambda, js, ja, jr);
            zpqrs += J1J2hat * phase_lambda * Z.modelspace->phase(jr + js - J1) * sixj * Y.OneBody(a, s) * X.TwoBody.GetTBME_J(J1, J1, p, q, r, a);
          }
          for (auto a : Y.GetOneBodyChannel(oR.l, oR.j2, oR.tz2))
          {
            Orbit &oa = Z.modelspace->GetOrbit(a);
            double ja = oa.j2 * 0.5;
            double sixj = AngMom::SixJ(J1, J2, lambda, jr, ja, js);
            zpqrs -= J1J2hat * phase_lambda * Z.modelspace->phase(J1 + J2) * sixj * Y.OneBody(a, r) * X.TwoBody.GetTBME_J(J1, J1, p, q, s, a);
          }

          // Normalize
          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;

          Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zpqrs);
        } // for iket
      } // for ibra

    } // for ch_bra/ch_ket
    // Z.PrintTwoBody();
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  }

  // This has not yet been validated, and is almost certainly wrong.
  void comm221st(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    int lambda = Y.GetJRank();
    for (auto i : Z.modelspace->all_orbits)
    {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double ji = 0.5 * oi.j2;
      for (auto j : Z.modelspace->all_orbits)
      {
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double jj = 0.5 * oj.j2;
        double zij = 0;

        for (auto c : Z.modelspace->all_orbits)
        {
          Orbit &oc = Z.modelspace->GetOrbit(c);
          double jc = oc.j2 * 0.5;
          int J1min = std::abs(oi.j2 - oc.j2) / 2;
          int J1max = (oi.j2 + oc.j2) / 2;
          int J2min = std::abs(oj.j2 - oc.j2) / 2;
          int J2max = (oj.j2 + oc.j2) / 2;
          for (int J1 = J1min; J1 <= J1max; J1++)
          {
            for (int J2 = J2min; J2 <= J2max; J2++)
            {
              double sixj = AngMom::SixJ(J1, J2, lambda, jj, ji, jc);
              double hats = sqrt((2 * J1 + 1) * (2 * J2 + 1));
              int phase = AngMom::phase(jj + jc + J1 + lambda);

              for (auto a : Z.modelspace->all_orbits)
              {
                Orbit &oa = Z.modelspace->GetOrbit(a);
                for (auto b : Z.modelspace->all_orbits)
                {
                  Orbit &ob = Z.modelspace->GetOrbit(b);
                  double nanbnc = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;

                  double Xciab = X.TwoBody.GetTBME_J(J1, J1, c, i, a, b);
                  double Xabcj = X.TwoBody.GetTBME_J(J2, J2, a, b, c, j);
                  double Yciab = Y.TwoBody.GetTBME_J(J1, J2, c, i, a, b);
                  double Yabcj = Y.TwoBody.GetTBME_J(J1, J2, a, b, c, j);
                  zij += 0.5 * nanbnc * hats * phase * sixj * (Xciab * Yabcj - Yciab * Xabcj);

                } // for b
              } // for a

            } // for J2
          } // for J1
        } // for c
        Z.OneBody(i, j) += zij;

      } // for j
    } // for i

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  } // comm221st

  // [2,2]->2 ph commutator. Tested against m-scheme expressions and it agrees.
  // Formula:
  //
  //    Zpqrs_J1,J2,lam = sum_ab s sum_J3J4 hat{J1J2J3J4} (na-nb) [ [1 - Ppq(J1)][1-Prs(J2)]  * (-1)^(q+s+J2+J4) { p  s  J3  } Xbar_psab^J3 * Ybar_abrq^J3,J4,lam ]
  void comm222_phst(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
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
          if ((ch_bra == ch_ket) and iket > ibra)
            continue;
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          Orbit &oR = Z.modelspace->GetOrbit(r);
          Orbit &os = Z.modelspace->GetOrbit(s);
          double zpqrs = 0;

          // Loop over the 3 permutations  1, Ppq, Prs, PpqPrs, with the dummy indices 1,2,3,4.
          // This avoids copy-pasing the same code 4 times.
          std::vector<std::vector<size_t>> perms = {{p, q, r, s}, {q, p, r, s}, {p, q, s, r}, {q, p, s, r}};
          std::vector<int> phases = {+1,
                                     -AngMom::phase((op.j2 + oq.j2) / 2 - J1),
                                     -AngMom::phase((oR.j2 + os.j2) / 2 - J2),
                                     +AngMom::phase((op.j2 + oq.j2) / 2 - J1) * AngMom::phase((oR.j2 + os.j2) / 2 - J2)};

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
            double j1 = o1.j2 * 0.5;
            double j2 = o2.j2 * 0.5;
            double j3 = o3.j2 * 0.5;
            double j4 = o4.j2 * 0.5;

            //            int J3min = std::abs(o1.j2 - o4.j2) / 2;
            //            int J3max = (o1.j2 + o4.j2) / 2;
            //            int J4min = std::abs(o3.j2 - o2.j2) / 2;
            //            int J4max = (o3.j2 + o2.j2) / 2;
            int J3min = std::abs(j1 - j4);
            int J3max = j1 + j4;
            int J4min = std::abs(j2 - j3);
            int J4max = j2 + j3;
            for (int J3 = J3min; J3 <= J3max; J3++)
            {
              for (int J4 = J4min; J4 <= J4max; J4++)
              {
                if (not AngMom::Triangle(J3, J4, lambda))
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
                    double ja = oa.j2 * 0.5;
                    double jb = ob.j2 * 0.5;
                    if (not AngMom::Triangle(ja, jb, J3))
                      continue;

                    double Xbar_psab = 0;
                    double Ybar_abrq = 0;

                    int J5min = std::max(std::abs(o1.j2 - ob.j2), std::abs(oa.j2 - o4.j2)) / 2.;
                    int J5max = std::min(o1.j2 + ob.j2, oa.j2 + o4.j2) / 2.;
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
                        if (not AngMom::Triangle(J6, J7, lambda))
                          continue;
                        double ninejY = AngMom::NineJ(oa.j2 / 2., o2.j2 / 2., J6, ob.j2 / 2., o3.j2 / 2., J7, J3, J4, lambda);
                        double hats = sqrt((2 * J3 + 1) * (2 * J4 + 1) * (2 * J6 + 1) * (2 * J7 + 1));
                        int phaseY = AngMom::phase((ob.j2 + o2.j2) / 2 + J4 + J7);
                        Ybar_abrq -= hats * phaseY * ninejY * Y.TwoBody.GetTBME_J(J6, J7, a, I2, I3, b);
                      }
                    }
                    Zbar_1432 += nanb * Xbar_psab * Ybar_abrq;

                  } // for b
                } // for a

                double ninej = AngMom::NineJ(o1.j2 / 2., o4.j2 / 2., J3, o2.j2 / 2., o3.j2 / 2., J4, J1, J2, lambda);
                zpqrs += phaseperm * sqrt((2 * J1 + 1) * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1)) * ninej * AngMom::phase((o2.j2 + o4.j2) / 2 + J2 + J4) * Zbar_1432;

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
    } // for iter
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  }//comm222_phst

  void comm222_pp_hhst(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
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
          if ((ch_bra == ch_ket) and iket > ibra)
            continue;
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          Orbit &oR = Z.modelspace->GetOrbit(r);
          Orbit &os = Z.modelspace->GetOrbit(s);
          double zpqrs = 0;

          for (auto a : Z.modelspace->all_orbits)
          {
            for (auto b : Z.modelspace->all_orbits)
            {
              Orbit &oa = Z.modelspace->GetOrbit(a);
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double nanb = 1 - oa.occ - ob.occ;
              double Xpqab = X.TwoBody.GetTBME_J(J1, J1, p, q, a, b);
              double Xabrs = X.TwoBody.GetTBME_J(J2, J2, a, b, r, s);
              double Ypqab = Y.TwoBody.GetTBME_J(J1, J2, p, q, a, b);
              double Yabrs = Y.TwoBody.GetTBME_J(J1, J2, a, b, r, s);
              zpqrs += 0.5 * nanb * (Xpqab * Yabrs - Ypqab * Xabrs);
            } // for b
          } // for a

          if (p == q)
            zpqrs /= PhysConst::SQRT2;
          if (r == s)
            zpqrs /= PhysConst::SQRT2;

          Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zpqrs);

        } // for iket
      } // for ibra
    } // for ch_bra,ch_ket

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
  } // comm222_pp_hhst

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// scalar-tensor commutators with 3-body


} // namespace ReferenceImplementations
////////////////////////////////////////////////////////////////////
