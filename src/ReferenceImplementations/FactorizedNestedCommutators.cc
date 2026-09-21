
#include "ReferenceImplementations.hh"
#include "ModelSpace.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"
#include <deque>


/// Straightforward implementation of J-coupled commutator expressions
/// without optimizations. This should be benchmarked against the
/// mscheme implementation and then left untouched.
namespace ReferenceImplementations
{

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
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
                      double sixj1 = AngMom::phase(J1) * (2 * J1 + 1) * AngMom::SixJ(oa.j2 * 0.5, oi.j2 * 0.5, J4, oj.j2 * 0.5, ob.j2 * 0.5, J1);
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
                      double sixj2 = AngMom::phase(J2) * (2 * J2 + 1) * AngMom::SixJ(op.j2 * 0.5, oc.j2 * 0.5, J4, oj.j2 * 0.5, ob.j2 * 0.5, J2);
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
                      double sixj3 = AngMom::phase(J3) * (2 * J3 + 1) * AngMom::SixJ(oa.j2 * 0.5, oi.j2 * 0.5, J4, op.j2 * 0.5, oc.j2 * 0.5, J3);
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
                    //       double sixj1 = AngMom::SixJ( oa.j2*0.5, oi.j2*0.5, J4,   oj.j2*0.5, ob.j2*0.5, J1);
                    //       double sixj2 = AngMom::SixJ( op.j2*0.5, oc.j2*0.5, J4,   oj.j2*0.5, ob.j2*0.5, J2);
                    //       double sixj3 = AngMom::SixJ( oa.j2*0.5, oi.j2*0.5, J4,   op.j2*0.5, oc.j2*0.5, J3);
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
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
                          double sixj1 = AngMom::SixJ(oa.j2 * 0.5, oi.j2 * 0.5, J4, oj.j2 * 0.5, ob.j2 * 0.5, J1);
                          double sixj2 = AngMom::SixJ(ok.j2 * 0.5, op.j2 * 0.5, J4, oj.j2 * 0.5, ob.j2 * 0.5, J2);
                          double sixj3 = AngMom::SixJ(oa.j2 * 0.5, oi.j2 * 0.5, J4, ok.j2 * 0.5, op.j2 * 0.5, J3);
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

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
    double t_start = omp_get_wtime();
    Z.modelspace->PreCalculateSixJ();
    bool EraseOB = false;
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
//      for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
      for (auto q : Z.modelspace->all_orbits ) // delta_jp jq
      {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if ( oq.j2 != op.j2) continue;
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

                  if ((ob.l + od.l + oa.l + oc.l) % 2 != Eta.GetParity())
                    continue;
                  if ((oa.l + oc.l + ob.l + oe.l) % 2 != Eta.GetParity())
                    continue;
                  if ((oe.l + op.l + od.l + oq.l) % 2 != Gamma.GetParity())
                    continue;

                  if (std::abs(ob.tz2 + od.tz2 - oa.tz2 - oc.tz2) != Eta.GetTRank())
                    continue;
                  if (std::abs(oa.tz2 + oc.tz2 - ob.tz2 - oe.tz2) != Eta.GetTRank())
                    continue;
                  if (std::abs(oe.tz2 + op.tz2 - od.tz2 - oq.tz2) != Gamma.GetTRank())
                    continue;

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
          Z.OneBody(q, p) += 0.5 * hZ * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q

    } // for p
        std::cout << "diagram I  " << Z.OneBodyNorm() << std::endl;
        std::cout << Z.OneBody << std::endl;
    if (EraseOB)
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
//      for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
      for (auto &q : Z.modelspace->all_orbits)
      {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if ( oq.j2 != op.j2 ) continue;
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

                  if ((ob.l + od.l + oa.l + oc.l) % 2 != Eta.GetParity())
                    continue;
                  if ((oc.l + op.l + od.l + oe.l) % 2 != Eta.GetParity())
                    continue;
                  if ((oa.l + oe.l + ob.l + oq.l) % 2 != Gamma.GetParity())
                    continue;

                  if (std::abs(ob.tz2 + od.tz2 - oa.tz2 - oc.tz2) != Eta.GetTRank())
                    continue;
                  if (std::abs(oc.tz2 + op.tz2 - od.tz2 - oe.tz2) != Eta.GetTRank())
                    continue;
                  if (std::abs(oa.tz2 + oe.tz2 - ob.tz2 - oq.tz2) != Gamma.GetTRank())
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
                          double sixj = AngMom::SixJ(ja, jb, J3, jd, jc, J0);
                          sixj *= AngMom::SixJ(jp, je, J3, jd, jc, J1);
                          sixj *= AngMom::SixJ(jb, jp, J2, je, ja, J3);
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
          Z.OneBody(q, p) += hZ * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q

    } // for p
        std::cout << "diagram IIa " << Z.OneBodyNorm() << std::endl;
        std::cout << Z.OneBody << std::endl;
    if (EraseOB)
      Z.EraseOneBody();

    // ####################################################################################
    //  diagram II_b
    //
    //   Ib_pq = 1/4 \delta_{jp jq} / (2jp + 1) sum_{abcde J0}
    //          (\barn_a \barn_d nb ne -\barn_b \barn_e na nd ) (2J_0 + 1)
    //          eta^J0_bdac  eta^J0_pcbe   Gamma^J0_adcq
    // ####################################################################################
    arma::mat fII_b = Z.OneBody * 0;
    for (auto &p : Z.modelspace->all_orbits)
    {
      Orbit &op = Z.modelspace->GetOrbit(p);
      double jp = op.j2 / 2.;
//      for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
      for (auto q : Z.modelspace->all_orbits)
      {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if ( oq.j2 != op.j2 ) continue;
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
                  if ((ob.l + oe.l + oa.l + od.l) % 2 != Eta.GetParity())
                    continue;
                  if ((oc.l + op.l + ob.l + oe.l) % 2 != Eta.GetParity())
                    continue;
                  if ((oa.l + od.l + oc.l + oq.l) % 2 != Gamma.GetParity())
                    continue;

                  if (std::abs(ob.tz2 + oe.tz2 - oa.tz2 - od.tz2) != Eta.GetTRank())
                    continue;
                  if (std::abs(oc.tz2 + op.tz2 - ob.tz2 - oe.tz2) != Eta.GetTRank())
                    continue;
                  if (std::abs(oa.tz2 + od.tz2 - oc.tz2 - oq.tz2) != Gamma.GetTRank())
                    continue;

                  int J0min = std::max({std::abs(oa.j2 - od.j2), std::abs(ob.j2 - oe.j2), std::abs(oc.j2 - op.j2)}) / 2;
                  int J0max = std::min({oa.j2 + od.j2, ob.j2 + oe.j2, op.j2 + oc.j2}) / 2;

                  for (int J0 = J0min; J0 <= J0max; J0++)
                  {
                    zij += (2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, e, a, d) * Eta.TwoBody.GetTBME_J(J0, J0, c, p, b, e) * Gamma.TwoBody.GetTBME_J(J0, J0, a, d, c, q);
                    zij += (2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, e, a, d) * Eta.TwoBody.GetTBME_J(J0, J0, c, q, b, e) * Gamma.TwoBody.GetTBME_J(J0, J0, a, d, c, p);


//                    if (p==0 and q==0 and a==0 and d==0 and c==0 and J0==0)
//                    {
//                        std::cout << "abcJ = " << a << " " << d << " " << c << " " << J0 << " Gamma = " << Gamma.TwoBody.GetTBME_J(J0, J0, a, d, c, q)
//                                  << "  chi = " << 0.25*(2 * J0 + 1) * occfactor * Eta.TwoBody.GetTBME_J(J0, J0, b, e, a, d) * Eta.TwoBody.GetTBME_J(J0, J0, c, p, b, e)
//                                  << " be = " << b << " " << e
//                                  <<  "    zij = " << zij << "=>  " <<  0.25*zij / (op.j2+1.0)<< std::endl;
//                    }

                  } // J0
                }
              }
            }
          }
        }
        fII_b(p,q) = 0.25 * zij / (op.j2 + 1.0);
        fII_b(q,p) = hZ * fII_b(p,q);
        Z.OneBody(p, q) += 0.25 * zij / (op.j2 + 1.0);
        if (p != q)
          Z.OneBody(q, p) += hZ * 0.25 * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q
    } // for p
        std::cout << "diagram IIb " << std::endl << fII_b << std::endl;
//        std::cout << "diagram IIb " << Z.OneBodyNorm() << std::endl;
        std::cout << Z.OneBody << std::endl;
    if (EraseOB)
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
//      for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
      for (auto &q : Z.modelspace->all_orbits)
      {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if ( oq.j2 != op.j2 ) continue;
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

                  if ((ob.l + oe.l + od.l + oa.l) % 2 != Eta.GetParity())
                    continue;
                  if ((oc.l + oa.l + ob.l + oq.l) % 2 != Eta.GetParity())
                    continue;
                  if ((od.l + op.l + oc.l + oe.l) % 2 != Gamma.GetParity())
                    continue;

                  if (std::abs(ob.tz2 + oe.tz2 - od.tz2 - oa.tz2) != Eta.GetTRank())
                    continue;
                  if (std::abs(oc.tz2 + oa.tz2 - ob.tz2 - oq.tz2) != Eta.GetTRank())
                    continue;
                  if (std::abs(od.tz2 + op.tz2 - oc.tz2 - oe.tz2) != Gamma.GetTRank())
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
                          double sixj = AngMom::SixJ(jd, je, J3, jb, ja, J0);
                          sixj *= AngMom::SixJ(jc, jp, J3, jb, ja, J1);
                          sixj *= AngMom::SixJ(je, jc, J2, jp, jd, J3);
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
          Z.OneBody(q, p) += hZ * zij / (op.j2 + 1.0);
        //--------------------------------------------------
      } // for q

    } // for p
        std::cout << "diagram IIc " << Z.OneBodyNorm() << std::endl;
        std::cout << Z.OneBody << std::endl;
    if (EraseOB)
      Z.EraseOneBody();


    // SRS This diagram is almost exactly the same as II_b, so it's easier to combine them
    // ####################################################################################
    //  diagram II_d
    //
    //   IId_pq = - 1/4 \delta_{jp jq} / (2jp + 1) sum_{abcde J0}
    //          (\barn_c \barn_d na ne -\barn_a \barn_e nc nd ) (2J_0 + 1)
    //          eta^J0_aecd  eta^J0_cdbq  Gamma^J0_bqae
    // ####################################################################################
/*
    for (auto &p : Z.modelspace->all_orbits)
    {
      Orbit &op = Z.modelspace->GetOrbit(p);
      double jp = op.j2 / 2.;
//      for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
      for (auto q : Z.modelspace->all_orbits)
      {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if ( oq.j2 != op.j2 ) continue;
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
    } // for p
        std::cout << "diagram IId " << Z.OneBodyNorm() << std::endl;
        std::cout << Z.OneBody << std::endl;
    if (EraseOB)
      Z.EraseOneBody();

*/

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
//      for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
      for (auto q : Z.modelspace->all_orbits)
      {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if ( oq.j2 != op.j2 ) continue;
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
    } // for p
        std::cout << "diagram IIIa " << Z.OneBodyNorm() << std::endl;
        std::cout << Z.OneBody << std::endl;
    if (EraseOB)
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
//      for (auto &q : Z.GetOneBodyChannel(op.l, op.j2, op.tz2)) // delta_jp jq
      for (auto q : Z.modelspace->all_orbits)
      {
        if (q > p)
          continue;
        Orbit &oq = Z.modelspace->GetOrbit(q);
        if ( oq.j2 != op.j2 ) continue;
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
    } // for p
        std::cout << "diagram IIIb " << Z.OneBodyNorm() << std::endl;
        std::cout << Z.OneBody << std::endl;
    if (EraseOB)
      Z.EraseOneBody();

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
    return;
  }




  void comm223_232_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {
    double t_start = omp_get_wtime();
    // global variables
    Z.modelspace->PreCalculateSixJ();
    int norbits = Z.modelspace->all_orbits.size();
    std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
    auto &Z2 = Z.TwoBody;
    bool EraseTB = false;
    // EraseTB = true;

    // determine symmetry
    int hEta = Eta.IsHermitian() ? 1 : -1;
    int hGamma = Gamma.IsHermitian() ? 1 : -1;
    // int hZ = Z.IsHermitian() ? 1 : -1;
    int hZ = hGamma;

    std::vector<size_t> ch_bra_list, ch_ket_list;
    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();
    // int nch = Z.modelspace->GetNumberTwoBodyChannels(); // number of TB channels

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
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);

      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J0 = tbc_bra.J; // J scalar
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);

        int phase_pg = bra.Phase(J0);
        double denominator_p = (op.j2 + 1.0);
        double denominator_g = (og.j2 + 1.0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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
                    }

                  // exchanging  p <-> g
                  j2min = std::max({std::abs(oa.j2 - oc.j2), std::abs(ob.j2 - od.j2), std::abs(ob.j2 - og.j2)}) / 2;
                  j2max = std::min({oa.j2 + oc.j2, ob.j2 + od.j2, ob.j2 + og.j2}) / 2;
                  if (delta_jgjd)
                    for (int J2 = j2min; J2 <= j2max; J2++)
                    {
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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram Ia " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();

      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;

        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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
                    }
                }
              }
            }
          } // a

          if (p == g)
            zpgqh /= PhysConst::SQRT2;
          if (q == h)
            zpgqh /= PhysConst::SQRT2;

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram Ib " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();

      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IVa " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        int phase_pg = bra.Phase(J0);
        double denominator_p = (op.j2 + 1.0);
        double denominator_g = (og.j2 + 1.0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IVb " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;

        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                        double sixj1 = AngMom::SixJ(jd, jg, J4, jc, jb, J2);
                        double sixj2 = AngMom::SixJ(jp, ja, J4, jc, jb, J3);
                        double sixj3 = AngMom::SixJ(jg, jp, J0, ja, jd, J4);

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

                        double sixj1 = AngMom::SixJ(jd, jp, J4, jc, jb, J2);
                        double sixj2 = AngMom::SixJ(jg, ja, J4, jc, jb, J3);
                        double sixj3 = AngMom::SixJ(jp, jg, J0, ja, jd, J4);

                        zpgqh -= phase_pg * occfactor * sixj1 * sixj2 * sixj3 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * Eta.TwoBody.GetTBME_J(J2, c, p, d, b) * Eta.TwoBody.GetTBME_J(J3, g, b, c, a) * Gamma.TwoBody.GetTBME_J(J0, d, a, q, h);
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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIa " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();

      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                        double sixj1 = AngMom::SixJ(jq, jd, J4, jc, jb, J2);
                        double sixj2 = AngMom::SixJ(ja, jh, J4, jc, jb, J3);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, ja, jd, J4);

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

                        double sixj1 = AngMom::SixJ(jh, jd, J4, jc, jb, J2);
                        double sixj2 = AngMom::SixJ(ja, jq, J4, jc, jb, J3);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, ja, jd, J4);

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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIc " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();

      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                          double sixj1 = AngMom::SixJ(jq, jd, J5, jc, jb, J2);
                          double sixj2 = AngMom::SixJ(jp, ja, J5, jc, jb, J3);
                          double sixj3 = AngMom::SixJ(J0, J5, J4, jd, jh, jq);
                          double sixj4 = AngMom::SixJ(J0, J4, J5, ja, jp, jg);

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

                          double sixj1 = AngMom::SixJ(jh, jd, J5, jc, jb, J2);
                          double sixj2 = AngMom::SixJ(jp, ja, J5, jc, jb, J3);
                          double sixj3 = AngMom::SixJ(J0, J5, J4, jd, jq, jh);
                          double sixj4 = AngMom::SixJ(J0, J4, J5, ja, jp, jg);

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

                          double sixj1 = AngMom::SixJ(jq, jd, J5, jc, jb, J2);
                          double sixj2 = AngMom::SixJ(jg, ja, J5, jc, jb, J3);
                          double sixj3 = AngMom::SixJ(J0, J5, J4, jd, jh, jq);
                          double sixj4 = AngMom::SixJ(J0, J4, J5, ja, jg, jp);

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

                          double sixj1 = AngMom::SixJ(jh, jd, J5, jc, jb, J2);
                          double sixj2 = AngMom::SixJ(jg, ja, J5, jc, jb, J3);
                          double sixj3 = AngMom::SixJ(J0, J5, J4, jd, jq, jh);
                          double sixj4 = AngMom::SixJ(J0, J4, J5, ja, jg, jp);

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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIb " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                          double sixj1 = AngMom::SixJ(jd, jg, J5, jc, jb, J2);
                          double sixj2 = AngMom::SixJ(ja, jh, J5, jc, jb, J3);
                          double sixj3 = AngMom::SixJ(J0, J5, J4, jd, jp, jg);
                          double sixj4 = AngMom::SixJ(J5, J4, J0, jq, jh, ja);

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

                          double sixj1 = AngMom::SixJ(jd, jg, J5, jc, jb, J2);
                          double sixj2 = AngMom::SixJ(ja, jq, J5, jc, jb, J3);
                          double sixj3 = AngMom::SixJ(J0, J5, J4, jd, jp, jg);
                          double sixj4 = AngMom::SixJ(J5, J4, J0, jh, jq, ja);

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

                          double sixj1 = AngMom::SixJ(jd, jp, J5, jc, jb, J2);
                          double sixj2 = AngMom::SixJ(ja, jh, J5, jc, jb, J3);
                          double sixj3 = AngMom::SixJ(J0, J5, J4, jd, jg, jp);
                          double sixj4 = AngMom::SixJ(J5, J4, J0, jq, jh, ja);

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

                          double sixj1 = AngMom::SixJ(jd, jp, J5, jc, jb, J2);
                          double sixj2 = AngMom::SixJ(ja, jq, J5, jc, jb, J3);
                          double sixj3 = AngMom::SixJ(J0, J5, J4, jd, jg, jp);
                          double sixj4 = AngMom::SixJ(J5, J4, J0, jh, jq, ja);

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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IId " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                        double sixj1 = AngMom::SixJ(jp, jh, J4, jb, jd, J2);
                        double sixj2 = AngMom::SixJ(jq, jg, J4, jb, jd, J3);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, jg, jp, J4);

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

                        double sixj1 = AngMom::SixJ(jp, jq, J4, jb, jd, J2);
                        double sixj2 = AngMom::SixJ(jh, jg, J4, jb, jd, J3);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, jg, jp, J4);

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

                        double sixj1 = AngMom::SixJ(jg, jh, J4, jb, jd, J2);
                        double sixj2 = AngMom::SixJ(jq, jp, J4, jb, jd, J3);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, jp, jg, J4);

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

                        double sixj1 = AngMom::SixJ(jg, jq, J4, jb, jd, J2);
                        double sixj2 = AngMom::SixJ(jh, jp, J4, jb, jd, J3);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, jp, jg, J4);

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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIe " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                        double sixj1 = AngMom::SixJ(jh, jp, J4, jb, jd, J2);
                        double sixj2 = AngMom::SixJ(jg, jq, J4, jb, jd, J3);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, jg, jp, J4);

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

                        double sixj1 = AngMom::SixJ(jq, jp, J4, jb, jd, J2);
                        double sixj2 = AngMom::SixJ(jg, jh, J4, jb, jd, J3);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, jg, jp, J4);

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

                        double sixj1 = AngMom::SixJ(jh, jg, J4, jb, jd, J2);
                        double sixj2 = AngMom::SixJ(jp, jq, J4, jb, jd, J3);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, jp, jg, J4);

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

                        double sixj1 = AngMom::SixJ(jq, jg, J4, jb, jd, J2);
                        double sixj2 = AngMom::SixJ(jp, jh, J4, jb, jd, J3);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, jp, jg, J4);

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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIf " << Z.TwoBodyNorm() << std::endl;
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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                          double sixj1 = AngMom::SixJ(ja, jb, J5, jp, jc, J2);
                          double sixj2 = AngMom::SixJ(J3, J0, J5, jp, jc, jg);
                          double sixj3 = AngMom::SixJ(jq, jd, J5, ja, jb, J4);
                          double sixj4 = AngMom::SixJ(J3, J0, J5, jq, jd, jh);

                          zpgqh += occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, p, c, a) * Eta.TwoBody.GetTBME_J(J3, g, c, h, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, q);
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

                          double sixj1 = AngMom::SixJ(ja, jb, J5, jp, jc, J2);
                          double sixj2 = AngMom::SixJ(J3, J0, J5, jp, jc, jg);
                          double sixj3 = AngMom::SixJ(jh, jd, J5, ja, jb, J4);
                          double sixj4 = AngMom::SixJ(J3, J0, J5, jh, jd, jq);

                          zpgqh += phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, p, c, a) * Eta.TwoBody.GetTBME_J(J3, g, c, q, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, h);
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

                          double sixj1 = AngMom::SixJ(ja, jb, J5, jg, jc, J2);
                          double sixj2 = AngMom::SixJ(J3, J0, J5, jg, jc, jp);
                          double sixj3 = AngMom::SixJ(jq, jd, J5, ja, jb, J4);
                          double sixj4 = AngMom::SixJ(J3, J0, J5, jq, jd, jh);

                          zpgqh += phase_pg * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, g, c, a) * Eta.TwoBody.GetTBME_J(J3, p, c, h, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, q);
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

                          double sixj1 = AngMom::SixJ(ja, jb, J5, jg, jc, J2);
                          double sixj2 = AngMom::SixJ(J3, J0, J5, jg, jc, jp);
                          double sixj3 = AngMom::SixJ(jh, jd, J5, ja, jb, J4);
                          double sixj4 = AngMom::SixJ(J3, J0, J5, jh, jd, jq);

                          zpgqh += phase_pg * phase_qh * occfactor * sixj1 * sixj2 * sixj3 * sixj4 * (2 * J2 + 1) * (2 * J3 + 1) * (2 * J4 + 1) * (2 * J5 + 1) * Eta.TwoBody.GetTBME_J(J2, b, g, c, a) * Eta.TwoBody.GetTBME_J(J3, p, c, q, d) * Gamma.TwoBody.GetTBME_J(J4, d, a, b, h);
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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpgqh);

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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                          double sixj1 = AngMom::SixJ(ja, jb, J5, jc, jq, J2);
                          double sixj2 = AngMom::SixJ(J0, J3, J5, jc, jq, jh);
                          double sixj3 = AngMom::SixJ(jd, jp, J5, ja, jb, J4);
                          double sixj4 = AngMom::SixJ(J0, J3, J5, jd, jp, jg);

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

                          double sixj1 = AngMom::SixJ(ja, jb, J5, jc, jh, J2);
                          double sixj2 = AngMom::SixJ(J0, J3, J5, jc, jh, jq);
                          double sixj3 = AngMom::SixJ(jd, jp, J5, ja, jb, J4);
                          double sixj4 = AngMom::SixJ(J0, J3, J5, jd, jp, jg);

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

                          double sixj1 = AngMom::SixJ(ja, jb, J5, jc, jq, J2);
                          double sixj2 = AngMom::SixJ(J0, J3, J5, jc, jq, jh);
                          double sixj3 = AngMom::SixJ(jd, jg, J5, ja, jb, J4);
                          double sixj4 = AngMom::SixJ(J0, J3, J5, jd, jg, jp);

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

                          double sixj1 = AngMom::SixJ(ja, jb, J5, jc, jh, J2);
                          double sixj2 = AngMom::SixJ(J0, J3, J5, jc, jh, jq);
                          double sixj3 = AngMom::SixJ(jd, jg, J5, ja, jb, J4);
                          double sixj4 = AngMom::SixJ(J0, J3, J5, jd, jg, jp);

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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpgqh);

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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();

      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                        double sixj1 = AngMom::SixJ(jq, jb, J4, jc, ja, J2);
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

                        double sixj1 = AngMom::SixJ(jh, jb, J4, jc, ja, J2);
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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpgqh);

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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                        double sixj1 = AngMom::SixJ(ja, jp, J4, jb, jc, J2);
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

                        double sixj1 = AngMom::SixJ(ja, jg, J4, jb, jc, J2);
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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zpgqh);

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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                        double sixj1 = AngMom::SixJ(J2, J3, J4, jp, jh, jd);
                        double sixj2 = AngMom::SixJ(J2, J3, J4, jq, jg, jc);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, jg, jp, J4);
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

                        double sixj1 = AngMom::SixJ(J2, J3, J4, jp, jq, jd);
                        double sixj2 = AngMom::SixJ(J2, J3, J4, jh, jg, jc);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, jg, jp, J4);

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

                        double sixj1 = AngMom::SixJ(J2, J3, J4, jg, jh, jd);
                        double sixj2 = AngMom::SixJ(J2, J3, J4, jq, jp, jc);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, jp, jg, J4);

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

                        double sixj1 = AngMom::SixJ(J2, J3, J4, jg, jq, jd);
                        double sixj2 = AngMom::SixJ(J2, J3, J4, jh, jp, jc);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, jp, jg, J4);

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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zpgqh);

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
      // TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;
        Orbit &op = *(bra.op);
        Orbit &og = *(bra.oq);
        double jp = op.j2 * 0.5;
        double jg = og.j2 * 0.5;
        int phase_pg = bra.Phase(J0);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
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

                        double sixj1 = AngMom::SixJ(J2, J3, J4, jq, jg, jc);
                        double sixj2 = AngMom::SixJ(J2, J3, J4, jp, jh, jd);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, jg, jp, J4);

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

                        double sixj1 = AngMom::SixJ(J2, J3, J4, jh, jg, jc);
                        double sixj2 = AngMom::SixJ(J2, J3, J4, jp, jq, jd);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, jg, jp, J4);

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

                        double sixj1 = AngMom::SixJ(J2, J3, J4, jq, jp, jc);
                        double sixj2 = AngMom::SixJ(J2, J3, J4, jg, jh, jd);
                        double sixj3 = AngMom::SixJ(jh, jq, J0, jp, jg, J4);

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

                        double sixj1 = AngMom::SixJ(J2, J3, J4, jh, jp, jc);
                        double sixj2 = AngMom::SixJ(J2, J3, J4, jg, jq, jd);
                        double sixj3 = AngMom::SixJ(jq, jh, J0, jp, jg, J4);

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

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, 0.5 * zpgqh);

        } // iket
      } // ibra
    } // J0 channel
    std::cout << "diagram IIIf " << Z.TwoBodyNorm() << std::endl;
    if (EraseTB)
      Z.EraseTwoBody();

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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

      for (auto &e : Eta.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
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
//          for (auto &e : Z.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
          for (auto &e : Eta.GetOneBodyChannel(od.l, od.j2, od.tz2)) // delta_jd je
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

    Z.profiler.timer["231_diagram_I"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime(); // timer

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

            double sixj1 = AngMom::SixJ(ja, jb, J_cc, jc, jd, J_std);
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

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
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
      // std::cout << iter.first[0] << "    " << iter.first[1] <<std::endl;
    }
    size_t nch = ch_bra_list.size();

    int nch_eta = Eta.modelspace->GetNumberTwoBodyChannels();
    for (int ch = 0; ch < nch_eta; ++ch)
    {
      TwoBodyChannel &tbc = Eta.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      // std::cout << ch << "    " << nKets << "    " << tbc.J <<std::endl;
    }

    bool Z_is_scalar = true;
    if (Z.TwoBody.rank_T != 0)
    {
      Z_is_scalar = false;
    }
    int n_nonzero = Eta.modelspace->GetNumberTwoBodyChannels_CC(); // number of CC channels
    auto &Z2 = Z.TwoBody;

    // determine symmetry
    int hEta = Eta.IsHermitian() ? 1 : -1;
    int hGamma = Gamma.IsHermitian() ? 1 : -1;
    int hZ = hGamma;
    // ####################################################################################
    //                      Factorization of Ia, Ib, IVa and IVb
    // ####################################################################################

    arma::mat CHI_I = Eta.OneBody * 0;
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
      for (auto q : Eta.OneBodyChannels.at({op.l, op.j2, op.tz2}))
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

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (size_t iket = ketmin; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t r = ket.p;
          size_t s = ket.q;
          Orbit &oR = Z.modelspace->GetOrbit(r);
          Orbit &os = Z.modelspace->GetOrbit(s);
          double zpqrs = 0;

          for (auto b : Eta.OneBodyChannels.at({op.l, op.j2, op.tz2}))
          {
            zpqrs += CHI_I(p, b) * Gamma.TwoBody.GetTBME_J(J, J, b, q, r, s);
            zpqrs += hZ * CHI_II(b, p) * Eta.TwoBody.GetTBME_J(J, J, b, q, r, s); // tricky minus sign.
          } // for a
          for (auto b : Eta.OneBodyChannels.at({oq.l, oq.j2, oq.tz2}))
          {
            zpqrs += CHI_I(q, b) * Gamma.TwoBody.GetTBME_J(J, J, p, b, r, s);
            zpqrs += hZ * CHI_II(b, q) * Eta.TwoBody.GetTBME_J(J, J, p, b, r, s); // tricky minus sign.
          } // for a
          for (auto b : Eta.OneBodyChannels.at({oR.l, oR.j2, oR.tz2}))
          {
            zpqrs += Gamma.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_I(b, r);
            zpqrs -= Eta.TwoBody.GetTBME_J(J, J, p, q, b, s) * CHI_II(b, r);
          } // for a
          for (auto b : Eta.OneBodyChannels.at({os.l, os.j2, os.tz2}))
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
    } // for channels

    Z.profiler.timer[std::string(__func__) + " Diagram I and Diagram IV"] += omp_get_wtime() - t_type;
    t_type = omp_get_wtime();

    Z.profiler.timer[std::string(__func__) + " Diagram I and Diagram IV"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
    // std::cout << "diagram I and IV " << Z.TwoBodyNorm() << std::endl;

    // *********************************************************************************** //
    //                                  Diagram II                                         //
    // *********************************************************************************** //
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
      if (nKets_cc == 0)
        continue;

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
            double sixj1 = AngMom::SixJ(ja, jb, J_cc, jc, jd, J_std);
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

    std::deque<arma::mat> Eta_matrix(nch_eta);
    std::deque<arma::mat> Eta_matrix_c(nch_eta);
    std::deque<arma::mat> Eta_matrix_d(nch_eta);
    std::deque<arma::mat> Gamma_matrix(nch);
    for (int ch = 0; ch < nch_eta; ++ch)
    {
      size_t ch_bra, ch_ket;
      ch_bra = ch;
      ch_ket = ch;
      // find index
      size_t index_ch_Z = ch;
      for (size_t i = 0; i < ch_bra_list.size(); ++i)
      {
        if (ch_bra_list[i] == ch)
        {
          ch_bra = ch_bra_list[i];
          ch_ket = ch_ket_list[i];
          index_ch_Z = i;
          break;
        }
      }
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nKets = tbc_ket.GetNumberKets();

      // std::cout << ch_bra << "  " << ch_ket << "  " << nbras << "  " << nKets << std::endl;

      if (nbras == 0 or nKets == 0)
        continue;
      Eta_matrix[ch] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
      Eta_matrix_c[ch] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
      Eta_matrix_d[ch] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
      if (Z_is_scalar == false and ch_bra != ch_ket)
      {
        Gamma_matrix[index_ch_Z] = arma::mat(nbras * 2, nKets * 2, arma::fill::zeros);
      }
      else if (Z_is_scalar)
        Gamma_matrix[ch] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
    }

    // full matrix
#pragma omp parallel for
    for (int ch = 0; ch < nch_eta; ++ch)
    {
      TwoBodyChannel &tbc = Eta.modelspace->GetTwoBodyChannel(ch);
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
          Eta_matrix_c[ch](ibra, iket) = occfactor_k * EtaME;
          Eta_matrix_d[ch](ibra, iket) = occfactor_l * EtaME;
          if (Z_is_scalar)
            Gamma_matrix[ch](ibra, iket) = GammaME;
          if (i != j)
          {
            int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
            Eta_matrix[ch](ibra + nKets, iket) = phase * EtaME;
            Eta_matrix_c[ch](ibra + nKets, iket) = occfactor_k * phase * EtaME;
            Eta_matrix_d[ch](ibra + nKets, iket) = occfactor_l * phase * EtaME;
            if (Z_is_scalar)
              Gamma_matrix[ch](ibra + nKets, iket) = phase * GammaME;
            if (k != l)
            {
              phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
              Eta_matrix[ch](ibra + nKets, iket + nKets) = phase * EtaME;
              Eta_matrix_c[ch](ibra + nKets, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra + nKets, iket + nKets) = occfactor_k * phase * EtaME;
              if (Z_is_scalar)
                Gamma_matrix[ch](ibra + nKets, iket + nKets) = phase * GammaME;

              phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix[ch](ibra, iket + nKets) = phase * EtaME;
              Eta_matrix_c[ch](ibra, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra, iket + nKets) = occfactor_k * phase * EtaME;
              if (Z_is_scalar)
                Gamma_matrix[ch](ibra, iket + nKets) = phase * GammaME;
            }
          }
          else
          {
            if (k != l)
            {
              int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
              Eta_matrix[ch](ibra, iket + nKets) = phase * EtaME;
              Eta_matrix_c[ch](ibra, iket + nKets) = occfactor_l * phase * EtaME;
              Eta_matrix_d[ch](ibra, iket + nKets) = occfactor_k * phase * EtaME;
              if (Z_is_scalar)
                Gamma_matrix[ch](ibra, iket + nKets) = phase * GammaME;
            }
          }
        }
      }
    }

    if (not Z_is_scalar)
    {
#pragma omp parallel for
      for (int ch = 0; ch < nch; ++ch)
      {
        // find index
        size_t ch_bra = ch_bra_list[ch];
        size_t ch_ket = ch_ket_list[ch];
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
          double n_i = oi.occ;
          double bar_n_i = 1. - n_i;
          double n_j = oj.occ;
          double bar_n_j = 1. - n_j;

          for (int iket = 0; iket < nkets; ++iket)
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

            double GammaME = Gamma.TwoBody.GetTBME_J(J0, i, j, k, l);
            Gamma_matrix[ch](ibra, iket) = GammaME;

            if (i != j)
            {
              int phase = Z.modelspace->phase((ji + jj) / 2 + J0 + 1);
              Gamma_matrix[ch](ibra + nbras, iket) = phase * GammaME;

              if (k != l)
              {
                phase = Z.modelspace->phase((ji + jj + jk + jl) / 2);
                Gamma_matrix[ch](ibra + nbras, iket + nkets) = phase * GammaME;

                phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                Gamma_matrix[ch](ibra, iket + nkets) = phase * GammaME;
              }
            }
            else
            {
              if (k != l)
              {
                int phase = Z.modelspace->phase((jk + jl) / 2 + J0 + 1);
                Gamma_matrix[ch](ibra, iket + nkets) = phase * GammaME;
              }
            }
          }
        }
      }
    }

    Z.profiler.timer[std::string(__func__) + " Global array"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

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
      TwoBodyChannel_CC &tbc_cc = Eta.modelspace->GetTwoBodyChannel_CC(ch_cc);
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
    std::deque<arma::mat> Chi_III(nch_eta);
    for (int ch = 0; ch < nch_eta; ++ch)
    {
      TwoBodyChannel &tbc_ket = Eta.modelspace->GetTwoBodyChannel(ch);
      size_t nkets = tbc_ket.GetNumberKets();
      if (nkets == 0)
        continue;
      Chi_III[ch] = arma::mat(nkets, nkets * 2, arma::fill::zeros);
    }

// Inverse Pandya transformation
//  X^J_ijkl  = - ( 1- P_ij )  sum_J' (2J'+1)  { i j J }  \bar{X}^J'_il`kj`
//                                             { k l J'}
#pragma omp parallel for
    for (int ch = 0; ch < nch_eta; ++ch)
    {
      TwoBodyChannel &tbc_bra = Eta.modelspace->GetTwoBodyChannel(ch);
      size_t nkets = tbc_bra.GetNumberKets();
      if (nkets == 0)
        continue;

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      for (int ibra = 0; ibra < nkets; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = *(bra.op);
        Orbit &oj = *(bra.oq);
        int ji = oi.j2;
        int jj = oj.j2;

        for (int iket = 0; iket < nkets * 2; ++iket)
        {
          size_t k, l;
          if (iket < nkets)
          {
            Ket &ket = tbc_bra.GetKet(iket);
            k = ket.p;
            l = ket.q;
          }
          else
          {
            Ket &ket = tbc_bra.GetKet(iket - nkets);
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

            double sixj = AngMom::SixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
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
            double sixj = AngMom::SixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);

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
      size_t ch_bra = ch_bra_list[ch];
      size_t ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();

      if (nbras == 0 or nkets == 0)
        continue;

      int J0 = tbc_bra.J;
      // int nKets = tbc.GetNumberKets();
      arma::Mat<double> Multi_matirx(nbras, nkets * 2);
      arma::Mat<double> Multi_matirxII(nbras * 2, nkets);
      // Multi_matirx = Chi_III[ch] * Gamma_matrix[ch];

      Multi_matirx = Chi_III[ch_bra] * Gamma_matrix[ch];
      Multi_matirxII = (Gamma_matrix[ch]) * (Chi_III[ch_ket].t());

      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;

          double z_IIa = Multi_matirx(ibra, iket);
          // double z_IIc = hZ * Multi_matirx(iket, ibra);
          double z_IIc = hZ * Multi_matirxII(ibra, iket);
          if (Z_is_scalar)
            z_IIc = hZ * Multi_matirx(iket, ibra);

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
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, z_IIa + z_IIc);
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

            double sixj1 = AngMom::SixJ(ja, jb, J_cc, jc, jd, J_std);
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
      TwoBodyChannel_CC &tbc_cc_bra = Z.modelspace->GetTwoBodyChannel_CC(ch);
      int nbras = tbc_cc_bra.GetNumberKets();
      if (nbras < 1)
        continue;
      CHI_III_final[ch] = bar_Gamma[ch] * barCHI_III_RC[ch];
    }

//  Inverse Pandya transformation
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
            double sixj1 = AngMom::SixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
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
            double sixj1 = AngMom::SixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
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

    std::deque<arma::mat> CHI_IV(nch_eta);
    for (int ch = 0; ch < nch_eta; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      // Not symmetric
      CHI_IV[ch] = arma::mat(nKets * 2, nKets * 2, arma::fill::zeros);
    }

#pragma omp parallel for
    for (int ch = 0; ch < nch_eta; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J0 = tbc.J;
      int nKets = tbc.GetNumberKets();
      CHI_IV[ch] = Eta_matrix[ch] * Eta_matrix_c[ch] + (Eta_matrix[ch] * Eta_matrix_d[ch]).t();
    }

    // build bar_CHI_IV
    std::deque<arma::mat> bar_CHI_IV(n_nonzero);
    std::deque<arma::mat> bar_CHI_gamma(n_nonzero);
    /// initial bar_CHI_IV
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      if (nKets_cc == 0)
        continue;

      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_CHI_IV[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
      bar_CHI_gamma[ch_cc] = arma::mat(nKets_cc * 2, nKets_cc * 2, arma::fill::zeros);
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

            double sixj1 = AngMom::SixJ(ja, jb, J_cc, jc, jd, J_std);
            if (std::abs(sixj1) > 1e-8)
            {
              Xbar -= (2 * J_std + 1) * sixj1 * CHI_IV[ch](indx_ad, indx_cb);
            }
          }

          bar_CHI_IV[ch_cc](ibra_cc, iket_cc) = Xbar;
        }

        //-------------------
      }
    }

// calculate bat_chi_IV * bar_gamma
#pragma omp parallel for
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
      bar_CHI_gamma[ch_cc] = bar_CHI_IV[ch_cc] * bar_Gamma[ch_cc];
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
            double sixj = AngMom::SixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
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
            double sixj = AngMom::SixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
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
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_CHI_V[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
      bar_CHI_V_RC[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
    }

// build bar_CHI_V
#pragma omp parallel for
    for (size_t ch_cc = 0; ch_cc < n_nonzero; ch_cc++)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
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
      int nbras_cc = tbc_cc.GetNumberKets();
      if (nbras_cc < 1)
        continue;

      int J_cc = tbc_cc.J;
      for (int ibra_cc = 0; ibra_cc < nbras_cc * 2; ++ibra_cc)
      {
        int a, b;
        if (ibra_cc < nbras_cc)
        {
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
          a = bra_cc.p;
          b = bra_cc.q;
        }
        else
        {
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nbras_cc);
          b = bra_cc.p;
          a = bra_cc.q;
        }
        if (ibra_cc >= nbras_cc and a == b)
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
        for (int iket_cc = 0; iket_cc < nbras_cc * 2; ++iket_cc)
        {
          int c, d;
          if (iket_cc < nbras_cc)
          {
            Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
            c = ket_cc_cd.p;
            d = ket_cc_cd.q;
          }
          else
          {
            Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nbras_cc);
            d = ket_cc_cd.p;
            c = ket_cc_cd.q;
          }
          if (iket_cc >= nbras_cc and c == d)
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
            int nbras = tbc_cc_old.GetNumberKets();
            if (nbras < 1)
              continue;
            int indx_ad = tbc_cc_old.GetLocalIndex(std::min(int(a), int(d)), std::max(int(a), int(d)));
            int indx_bc = tbc_cc_old.GetLocalIndex(std::min(int(b), int(c)), std::max(int(b), int(c)));
            if (indx_ad < 0 or indx_bc < 0)
              continue;

            if (a > d)
              indx_ad += nbras;
            if (b > c)
              indx_bc += nbras;

            double sixj1 = AngMom::SixJ(ja, jb, J_cc, jc, jd, J_std);
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

    std::deque<arma::mat> CHI_V_final(n_nonzero);
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
      // Not symmetric
      CHI_V_final[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
    }

#pragma omp parallel for
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
      CHI_V_final[ch_cc] = bar_Eta[ch_cc] * bar_CHI_V_RC[ch_cc];
    }

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
            double sixj1 = AngMom::SixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
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
            double sixj1 = AngMom::SixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
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
    std::deque<arma::mat> bar_CHI_VId(n_nonzero);
    /// initial bar_Eta and nnnbar_Eta
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_CHI_VI[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
      bar_CHI_VId[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
    }

// BUILD bar_CHI_VI
// bar_CHI_VI = (bar{na} nb nc + bar{nb} bar{nc} na) \bar{Gamma}_dhba \bar{Eta}_baqc
#pragma omp parallel for
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
      bar_CHI_VI[ch_cc] = bar_Gamma[ch_cc] * nnnbar_Eta_d[ch_cc];
      bar_CHI_VId[ch_cc] = hEta * (nnnbar_Eta_d[ch_cc]).t() * bar_Gamma[ch_cc];
    }

    std::deque<arma::mat> CHI_VI(nch);
    std::deque<arma::mat> CHI_VI_II(nch);
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
      // Not symmetric
      CHI_VI[ch] = arma::mat(nbras * 2, nkets * 2, arma::fill::zeros);
      CHI_VI_II[ch] = arma::mat(nbras * 2, nkets * 2, arma::fill::zeros);
    }

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
        if (ibra >= nbras and i == j)
          continue;

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
          if (iket >= nkets and k == l)
            continue;

          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);
          int jk = ok.j2;
          int jl = ol.j2;

          double commijkl = 0;
          double commjikl = 0;
          double commijlk = 0;
          double commjilk = 0;

          double commijkld = 0;
          double commjikld = 0;
          double commijlkd = 0;
          double commjilkd = 0;

          // ijkl, direct term        -->  il kj
          int parity_cc = (oi.l + ol.l) % 2;
          int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
          int Jpmin = std::max(std::abs(ji - jl), std::abs(jj - jk)) / 2;
          int Jpmax = std::min(ji + jl, jj + jk) / 2;
          for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
          {
            double sixj = AngMom::SixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
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

            commijkld -= (2 * Jprime + 1) * sixj * bar_CHI_VId[ch_cc](indx_il, indx_kj);
          }

          // ijlk,  exchange k and l -->  ik lj
          parity_cc = (oi.l + ok.l) % 2;
          Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
          Jpmin = std::max(std::abs(ji - jk), std::abs(jj - jl)) / 2;
          Jpmax = std::min(ji + jk, jj + jl) / 2;
          for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
          {
            double sixj = AngMom::SixJ(ji * 0.5, jj * 0.5, J0, jl * 0.5, jk * 0.5, Jprime);
            if (std::abs(sixj) < 1e-8)
              continue;
            int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
            TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
            int nkets_cc = tbc_cc.GetNumberKets();
            if (nkets_cc < 1)
              continue;
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
            double me21 = bar_CHI_VId[ch_cc](indx_jl, indx_ki);
            commjikld -= (2 * Jprime + 1) * sixj * me21;
          }

          double zijkl = (commijkl - Z.modelspace->phase((jk + jl) / 2 - J0) * commijlk);
          CHI_VI[ch](ibra, iket) += zijkl;

          zijkl = (commijkld - Z.modelspace->phase((ji + jj) / 2 - J0) * commjikld);
          CHI_VI_II[ch](ibra, iket) += zijkl;
        }
      }
    }

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
      int J0 = tbc_bra.J;

      arma::Mat<double> Multi_matirx(nbras * 2, nkets * 2);
      Multi_matirx = Eta_matrix[ch_bra] * CHI_VI[ch] + (CHI_VI_II[ch] * Eta_matrix[ch_ket]);

      for (int ibra = 0; ibra < nbras; ++ibra)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t p = bra.p;
        size_t g = bra.q;

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; ++iket)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t q = ket.p;
          size_t h = ket.q;

          double zpgqhIIIcd = Multi_matirx(ibra, iket);
          if (p == g)
          {
            zpgqhIIIcd /= PhysConst::SQRT2;
          }
          if (q == h)
          {
            zpgqhIIIcd /= PhysConst::SQRT2;
          }
          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, -zpgqhIIIcd);
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
      int ch_bra = ch_bra_list[ch];
      int ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      if (nbras == 0 or nkets == 0)
        continue;
      // Not symmetric
      CHI_VII[ch] = arma::mat(nbras * 2, nkets * 2, arma::fill::zeros);
    }

#pragma omp parallel for
    for (int ch = 0; ch < nch; ++ch)
    {
      int ch_bra = ch_bra_list[ch];
      int ch_ket = ch_ket_list[ch];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      size_t nbras = tbc_bra.GetNumberKets();
      size_t nkets = tbc_ket.GetNumberKets();
      if (nbras == 0 or nkets == 0)
        continue;
      /// III_e
      CHI_VII[ch] = Gamma_matrix[ch] * Eta_matrix_d[ch_ket];
      /// III_f
      CHI_VII[ch] += hEta * Eta_matrix_d[ch_bra].t() * Gamma_matrix[ch];
    }

    std::deque<arma::mat> bar_CHI_VII_CC(n_nonzero);
    std::deque<arma::mat> bar_CHI_VII_CC_ef(n_nonzero);
    /// initial bar_Eta and nnnbar_Eta
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
      // because the restriction a<b in the bar and ket vector, if we want to store the full
      // Pandya transformed matrix, we twice the size of matrix
      bar_CHI_VII_CC[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
      bar_CHI_VII_CC_ef[ch_cc] = arma::mat(nbras * 2, nbras * 2, arma::fill::zeros);
    }

    /// Modified Pandya transformation
    //  \bar{X}^J_ij`kl` =   sum_J' (-)^{jk + jj + J'} { i j J } (2J'+1) X^J'_ilkj
    //                                                 { k l J'}
#pragma omp parallel for
    for (int ch_cc = 0; ch_cc < n_nonzero; ++ch_cc)
    {
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras_cc = tbc_cc.GetNumberKets();
      if (nbras_cc < 1)
        continue;

      int J_cc = tbc_cc.J;
      for (int ibra_cc = 0; ibra_cc < nbras_cc * 2; ++ibra_cc)
      {
        int a, b;
        if (ibra_cc < nbras_cc)
        {
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc);
          a = bra_cc.p;
          b = bra_cc.q;
        }
        else
        {
          Ket &bra_cc = tbc_cc.GetKet(ibra_cc - nbras_cc);
          b = bra_cc.p;
          a = bra_cc.q;
        }
        if (ibra_cc >= nbras_cc and a == b)
          continue;
        Orbit &oa = Z.modelspace->GetOrbit(a);
        Orbit &ob = Z.modelspace->GetOrbit(b);
        double ja = oa.j2 * 0.5;
        double jb = ob.j2 * 0.5;

        // loop over cross-coupled kets |cd> in this channel
        for (int iket_cc = 0; iket_cc < nbras_cc * 2; ++iket_cc)
        {
          int c, d;
          if (iket_cc < nbras_cc)
          {
            Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc);
            c = ket_cc_cd.p;
            d = ket_cc_cd.q;
          }
          else
          {
            Ket &ket_cc_cd = tbc_cc.GetKet(iket_cc - nbras_cc);
            d = ket_cc_cd.p;
            c = ket_cc_cd.q;
          }
          if (iket_cc >= nbras_cc and c == d)
            continue;
          Orbit &oc = Z.modelspace->GetOrbit(c);
          Orbit &od = Z.modelspace->GetOrbit(d);
          double jc = oc.j2 * 0.5;
          double jd = od.j2 * 0.5;

          int Tz_J2_bc = (ob.tz2 + oc.tz2) / 2;
          int Tz_J2_ad = (oa.tz2 + od.tz2) / 2;
          int parity_J2 = (ob.l + oc.l) % 2;

          if (std::abs(Tz_J2_ad - Tz_J2_bc) != Z.TwoBody.rank_T)
          {
            continue;
          }

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
            double sixj1 = AngMom::SixJ(ja, jb, J_cc, jc, jd, J_std);
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
              if (b > c)
                indx_bc += nkets_bc;

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
              Ybar += phaseFactor * (2 * J_std + 1) * sixj1 * CHI_VII[index_ch](indx_ad, indx_bc);
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
      TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nbras = tbc_cc.GetNumberKets();
      if (nbras < 1)
        continue;
      bar_CHI_VII_CC_ef[ch_cc] = bar_CHI_VII_CC[ch_cc] * bar_Eta[ch_cc];
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
            double sixj = AngMom::SixJ(jj * 0.5, ji * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
            double sixj2 = AngMom::SixJ(ji * 0.5, jj * 0.5, J0, jl * 0.5, jk * 0.5, Jprime);
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
            double sixj = AngMom::SixJ(ji * 0.5, jj * 0.5, J0, jk * 0.5, jl * 0.5, Jprime);
            double sixj2 = AngMom::SixJ(jj * 0.5, ji * 0.5, J0, jl * 0.5, jk * 0.5, Jprime);
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

    Z.profiler.timer[std::string(__func__) + " Diagram IIIe and IIIf  -c"] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();

    Z.profiler.timer[std::string(__func__) + " Diagram III"] += omp_get_wtime() - t_type;

    // Timer
    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
    return;
  }


  

  //////////////////////////////////////////////////
  // Equation B5a from PRC 110 044317
  void comm223_231_fI(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {

    double t_start = omp_get_wtime();
    // Build intermediate chi_ij equation B6a.
    arma::mat chi_alpha = 0* Z.OneBody;
    for (auto i : Z.modelspace->all_orbits )
    {
       Orbit& oi = Z.modelspace->GetOrbit(i);
//       for ( auto j : Z.modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2})  )
       for ( auto j : Z.modelspace->all_orbits  )
       {
         Orbit& oj = Z.modelspace->GetOrbit(j);
         if (oi.j2 != oj.j2) continue;
         double chi_ij = 0;
         for (auto a : Z.modelspace->all_orbits)
         {
            Orbit& oa = Z.modelspace->GetOrbit(a);
            for (auto b : Z.modelspace->all_orbits)
            {
               Orbit& ob = Z.modelspace->GetOrbit(b);
               for (auto c : Z.modelspace->all_orbits)
               {
                  Orbit& oc = Z.modelspace->GetOrbit(c);
                  int Jmin = AngMom::Jmin({ {oa.j2,ob.j2}, {oi.j2,oc.j2}  })/2;
                  int Jmax = AngMom::Jmax({ {oa.j2,ob.j2}, {oi.j2,oc.j2}  })/2;
                  for (int J=Jmin; J<=Jmax; J++)
                  {
                     double Omega_ciab = Eta.TwoBody.GetTBME_J(J,J,c,i,a,b);
                     double Omega_abcj = Eta.TwoBody.GetTBME_J(J,J,a,b,c,j);
                     double occfactor  = (1-oa.occ)*(1-ob.occ)*oc.occ*oi.occ - oa.occ*ob.occ*(1-oc.occ)*(1-oi.occ);
                            occfactor += (1-oa.occ)*(1-ob.occ)*oc.occ*oj.occ - oa.occ*ob.occ*(1-oc.occ)*(1-oj.occ);
                     chi_ij += 0.5*(2*J+1)/(oi.j2+1.) * occfactor * Omega_ciab * Omega_abcj;

                  }// for J
               }// for c
            }// for b
         }// for a
       chi_alpha(i,j) = chi_ij;
       }// for j
    }// for i

    // Done making intermediate chi_ij.
    std::cout << "chi_alpha = " << std::endl << chi_alpha << std::endl;
    
    for (auto i : Z.modelspace->all_orbits )
    {
       Orbit& oi = Z.modelspace->GetOrbit(i);
       for ( auto j : Z.GetOneBodyChannel(oi.l,oi.j2,oi.tz2)  )
       {
         Orbit& oj = Z.modelspace->GetOrbit(j);
         double fI_ij = 0;
         for (auto a : Z.modelspace->all_orbits)
         {
            Orbit& oa = Z.modelspace->GetOrbit(a);
            for (auto b : Z.modelspace->all_orbits)
            {
               Orbit& ob = Z.modelspace->GetOrbit(b);
               int Jmin = AngMom::Jmin({ {oa.j2,oj.j2}, {oi.j2,ob.j2}  })/2;
               int Jmax = AngMom::Jmax({ {oa.j2,oj.j2}, {oi.j2,ob.j2}  })/2;
               for (int J=Jmin; J<=Jmax; J++)
               {
                  double Gamma_biaj = Gamma.TwoBody.GetTBME_J(J,J, b,i,a,j);
                  fI_ij += (2*J+1)/(oi.j2+1.) * chi_alpha(a,b) * Gamma_biaj ;

               }// for J
            }// for b
         }// for a
         Z.OneBody(i,j) += fI_ij;
       }// for j
    }// for i

    std::cout << "fI : " << std::endl << Z.OneBody << std::endl;

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
    return;
  }



  //////////////////////////////////////////////////
  // Equation B5b from PRC 110 044317
  // Note that equation B6b for the intermediate has a typo.
  // It should have i <-> j swapped. It is fixed in this routine.
  void comm223_231_fII(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {

    double t_start = omp_get_wtime();
     // Build intermediate chi_ij equation B6b.
     arma::mat chi_beta = 0* Z.OneBody;
     for (auto i : Z.modelspace->all_orbits )
     {
        Orbit& oi = Z.modelspace->GetOrbit(i);
        for ( auto j : Gamma.GetOneBodyChannel(oi.l,oi.j2,oi.tz2)  )
        {
          Orbit& oj = Z.modelspace->GetOrbit(j);
          double chi_ij = 0;
          for (auto a : Z.modelspace->all_orbits)
          {
             Orbit& oa = Z.modelspace->GetOrbit(a);
             for (auto b : Z.modelspace->all_orbits)
             {
                Orbit& ob = Z.modelspace->GetOrbit(b);
                for (auto c : Z.modelspace->all_orbits)
                {
                   Orbit& oc = Z.modelspace->GetOrbit(c);
                   int Jmin = AngMom::Jmin({ {oa.j2,ob.j2}, {oi.j2,oc.j2}  })/2;
                   int Jmax = AngMom::Jmax({ {oa.j2,ob.j2}, {oi.j2,oc.j2}  })/2;
                   for (int J=Jmin; J<=Jmax; J++)
                   {
                      double Omega_cjab = Eta.TwoBody.GetTBME_J(J,J,c,j,a,b);
                      double Gamma_abci = Gamma.TwoBody.GetTBME_J(J,J,a,b,c,i);
                      double occfactor = (1-oa.occ)*(1-ob.occ)*oc.occ*oj.occ - oa.occ*ob.occ*(1-oc.occ)*(1-oj.occ);
                      chi_ij += (2*J+1)/(oi.j2+1.0) * occfactor * Omega_cjab * Gamma_abci;
                   }// for J
                }// for c
             }// for b
          }// for a
        chi_beta(i,j) = chi_ij;
        }// for j
     }// for i
    std::cout << "chi_beta = " << std::endl << chi_beta << std::endl;

     // Done making intermediate chi_ij.
     arma::mat fII = 0*Z.OneBody;
     
     for (auto i : Z.modelspace->all_orbits )
     {
        Orbit& oi = Z.modelspace->GetOrbit(i);
        for ( auto j : Z.GetOneBodyChannel(oi.l,oi.j2,oi.tz2)  )
        {
          Orbit& oj = Z.modelspace->GetOrbit(j);
          double fII_ij = 0;
          for (auto a : Z.modelspace->all_orbits)
          {
             Orbit& oa = Z.modelspace->GetOrbit(a);
             for (auto b : Z.modelspace->all_orbits)
             {
                Orbit& ob = Z.modelspace->GetOrbit(b);
                int Jmin = AngMom::Jmin({ {oa.j2,oj.j2}, {oi.j2,ob.j2}  })/2;
                int Jmax = AngMom::Jmax({ {oa.j2,oj.j2}, {oi.j2,ob.j2}  })/2;
                for (int J=Jmin; J<=Jmax; J++)
                {
                   double Omega_biaj = Eta.TwoBody.GetTBME_J(J,J, b,i,a,j);
                   fII_ij += 0.5*(2*J+1)/(oi.j2+1.) * ( chi_beta(a,b) - chi_beta(b,a) ) * Omega_biaj;  
                }// for J
             }// for b
          }// for a
          Z.OneBody(i,j) += fII_ij;
          fII(i,j) = fII_ij;
        }// for j
     }// for i

    std::cout << "fII : " << std::endl << fII << std::endl;
//    std::cout << "fII : " << std::endl << Z.OneBody << std::endl;

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
    return;
  }



  ////
  //// Eq B5c, with intermediate defined in B6c.
  void comm223_231_fIIIa(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {
    double t_start = omp_get_wtime();

    
    arma::mat fIIIa = 0*Z.OneBody;

    std::unordered_map<int,arma::mat> chi_g;
    size_t nch = Z.modelspace->GetNumberTwoBodyChannels_CC();
    for (size_t ch=0; ch<nch; ch++)
    {
       TwoBodyChannel_CC& tbc = Z.modelspace->GetTwoBodyChannel_CC(ch);
       int J = tbc.J;
       size_t nkets = tbc.GetNumberKets();

       arma::mat Omega_bar = arma::zeros(2*nkets,2*nkets);
       for (size_t IJ=0; IJ<nkets; IJ++)
       {
          Ket& ketij = tbc.GetKet(IJ);
          size_t i = ketij.p;
          size_t j = ketij.q;
          int j2i = ketij.op->j2 ;
          int j2j = ketij.oq->j2 ;
          double ji = j2i * 0.5;
          double jj = j2j * 0.5;
          for (size_t KL=0; KL<nkets; KL++)
          {
            Ket& ketkl = tbc.GetKet(KL);
            size_t k = ketkl.p;
            size_t l = ketkl.q;
            double nk = ketkl.op->occ;
            double nl = ketkl.oq->occ;
            int j2k = ketkl.op->j2 ;
            int j2l = ketkl.oq->j2 ;
            double jk = j2k * 0.5;
            double jl = j2l * 0.5;

            double omegabar_ijkl = 0;
            double omegabar_jikl = 0;
            double omegabar_ijlk = 0;
            double omegabar_jilk = 0;
            int JJmin = AngMom::Jmin({ {j2i,j2l} , {j2j,j2k}  }) / 2;
            int JJmax = AngMom::Jmax({ {j2i,j2l} , {j2j,j2k}  }) / 2;
            for (int JJ=JJmin; JJ<=JJmax; JJ++)
            {
                double sixj1 = Z.modelspace->GetSixJ( ji,jj,J,  jk,jl,JJ);
                double Omega_ilkj = Eta.TwoBody.GetTBME_J(JJ,JJ,i,l,k,j);
                double Omega_jkli = Eta.TwoBody.GetTBME_J(JJ,JJ,j,k,l,i);
                omegabar_ijkl += -(2*JJ+1) * sixj1 * Omega_ilkj;
                omegabar_jilk += -(2*JJ+1) * sixj1 * Omega_jkli;
            }

            JJmin = AngMom::Jmin({ {j2j,j2l} , {j2i,j2k}  }) / 2;
            JJmax = AngMom::Jmax({ {j2j,j2l} , {j2i,j2k}  }) / 2;
            for (int JJ=JJmin; JJ<=JJmax; JJ++)
            {
                double sixj2 = Z.modelspace->GetSixJ( jj,ji,J,  jk,jl,JJ);
                double Omega_jlki = Eta.TwoBody.GetTBME_J(JJ,JJ,j,l,k,i);
                double Omega_iklj = Eta.TwoBody.GetTBME_J(JJ,JJ,i,k,l,j);
                omegabar_jikl += -(2*JJ+1) * sixj2 * Omega_jlki;
                omegabar_ijlk += -(2*JJ+1) * sixj2 * Omega_iklj;
            }

          Omega_bar(IJ,       KL    )  = omegabar_ijkl;
          Omega_bar(IJ+nkets, KL    )  = omegabar_jikl;
          Omega_bar(IJ,       KL+nkets) = omegabar_ijlk;
          Omega_bar(IJ+nkets, KL+nkets) = omegabar_jilk;
          }// for iket
       }// for ibra

       chi_g[ch] = arma::zeros(2*nkets,2*nkets);
       for (size_t IJ=0; IJ<nkets; IJ++)
       {
          Ket& ketij = tbc.GetKet(IJ);
          size_t JI = IJ + nkets;
          for (size_t KL=0; KL<nkets; KL++)
          {
            Ket& ketkl = tbc.GetKet(KL);
            size_t LK = KL + nkets;
            double nk = ketkl.op->occ;
            double nl = ketkl.oq->occ;

            double chibar_ijkl = 0;
            double chibar_jikl = 0;
            double chibar_ijlk = 0;
            double chibar_jilk = 0;
            // we loop over a<=b, so we need a factor 2 to account for a>b
            // the case a==b is handled automatically because we're using normalized matrix elements
            for (size_t AB=0; AB<nkets; AB++)
            {
               Ket& ketab = tbc.GetKet(AB);
               size_t BA = AB + nkets;
               double na = ketab.op->occ;
               double nb = ketab.oq->occ;
               double occ_abkl = na*(1-nb)*(1-nk)*nl - (1-na)*nb*nk*(1-nl);
               double occ_bakl = nb*(1-na)*(1-nk)*nl - (1-nb)*na*nk*(1-nl);
               double occ_ablk = na*(1-nb)*(1-nl)*nk - (1-na)*nb*nl*(1-nk);
               double occ_balk = nb*(1-na)*(1-nl)*nk - (1-nb)*na*nl*(1-nk);
               double OmegaBar_ijab = Omega_bar( IJ, AB );
               double OmegaBar_ijba = Omega_bar( IJ, BA );
               double OmegaBar_jiab = Omega_bar( JI, AB );
               double OmegaBar_jiba = Omega_bar( JI, BA );
               double OmegaBar_abkl = Omega_bar( AB, KL );
               double OmegaBar_bakl = Omega_bar( BA, KL );
               double OmegaBar_ablk = Omega_bar( AB, LK );
               double OmegaBar_balk = Omega_bar( BA, LK );
               chibar_ijkl += (2*J+1)  * occ_abkl * OmegaBar_ijab * OmegaBar_abkl; 

               chibar_jikl += (2*J+1)  * occ_abkl * OmegaBar_jiab * OmegaBar_abkl; 

               chibar_ijlk += (2*J+1)  * occ_ablk * OmegaBar_ijab * OmegaBar_ablk; 

               chibar_jilk += (2*J+1)  * occ_ablk * OmegaBar_jiab * OmegaBar_ablk; 


               if ( ketab.p != ketab.q )
               {
               chibar_ijkl += (2*J+1)  * occ_bakl * OmegaBar_ijba * OmegaBar_bakl; 
               chibar_jikl += (2*J+1)  * occ_bakl * OmegaBar_jiba * OmegaBar_bakl; 
               chibar_ijlk += (2*J+1)  * occ_balk * OmegaBar_ijba * OmegaBar_balk; 
               chibar_jilk += (2*J+1)  * occ_balk * OmegaBar_jiba * OmegaBar_balk; 
               }

            }
            chi_g[ch](IJ,KL) = chibar_ijkl;
            chi_g[ch](JI,KL) = chibar_jikl;
            chi_g[ch](IJ,LK) = chibar_ijlk;
            chi_g[ch](JI,LK) = chibar_jilk;
          }// for iket
       }// for ibra
    }//for ch

    std::cout << "successfully build chi_g" << std::endl;


    for (auto i : Z.modelspace->all_orbits)
    {
       Orbit& oi = Z.modelspace->GetOrbit(i);
       for (auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
       {
          Orbit& oj = Z.modelspace->GetOrbit(j);
          double fIIIa_ij = 0;
          for (auto a : Z.modelspace->all_orbits)
          {
             Orbit& oa = Z.modelspace->GetOrbit(a);
             for (auto b : Z.modelspace->all_orbits)
             {
                Orbit& ob = Z.modelspace->GetOrbit(b);
                int parity = ( oa.l+ob.l )%2;
                int Tz = std::abs( oa.tz2 - ob.tz2)/2;
                for (auto c : Z.modelspace->all_orbits)
                {
                   Orbit& oc = Z.modelspace->GetOrbit(c);
                   if (  (oi.l+oc.l)%2 != parity ) continue;
                   if ( std::abs(oi.tz2-oc.tz2) != 2*Tz ) continue;
                   int Jmin = AngMom::Jmin( { {oa.j2,ob.j2}, {oc.j2,oi.j2}  }) /2;
                   int Jmax = AngMom::Jmax( { {oa.j2,ob.j2}, {oc.j2,oi.j2}  }) /2;
                   for (int J=Jmin; J<=Jmax; J++)
                   {
                       size_t ch = Z.modelspace->GetTwoBodyChannelIndex(J,parity,Tz);
                       TwoBodyChannel_CC& tbc = Z.modelspace->GetTwoBodyChannel_CC(ch);
                       size_t nkets = tbc.GetNumberKets();
                       size_t ab = tbc.GetLocalIndex( a,b );
                       size_t ba = tbc.GetLocalIndex( b,a );
                       size_t ic = tbc.GetLocalIndex( i,c );
                       size_t cj = tbc.GetLocalIndex( c,j );

//                       std::cout << "I want to access element ic,ab => " << ic << " " << ab
//                                 << "  : " << i << " " << c << " " << a << " " << b << "  with J = " << J << " " << parity << " "<< Tz
//                                 << " and chi_g[ch] has dimension " << chi_g[ch].n_rows << " " << chi_g[ch].n_cols << std::endl;
                       double chibar_icab = ( ic < 2*nkets ) ?  chi_g[ch](ic,ab) : 0;
                       double chibar_cjab = ( cj < 2*nkets ) ?  chi_g[ch](cj,ab) : 0;
//                       double chibar_icab = chi_g[ch](ic,ab);
//                       double chibar_cjab = chi_g[ch](cj,ab);

                       double Gammabar_abjc = 0;
                       double Gammabar_abci = 0;

                       int JJmin = AngMom::Jmin({ {oa.j2,oi.j2} , {ob.j2,oc.j2} })/2;
                       int JJmax = AngMom::Jmin({ {oa.j2,oi.j2} , {ob.j2,oc.j2} })/2;
                       for (int JJ=JJmin; JJ<=JJmax; JJ++)
                       {
                          double sixj2 = Z.modelspace->GetSixJ( oa.j2*0.5, ob.j2*0.5, J,  oc.j2*0.5, oi.j2*0.5, JJ);
                          Gammabar_abci -= (2*JJ+1) * sixj2 * Gamma.TwoBody.GetTBME_J(JJ,JJ,a,i,c,b);
                       }

                       JJmin = AngMom::Jmin({ {oa.j2,oc.j2} , {ob.j2,oj.j2} })/2;
                       JJmax = AngMom::Jmin({ {oa.j2,oc.j2} , {ob.j2,oj.j2} })/2;
                       for (int JJ=JJmin; JJ<=JJmax; JJ++)
                       {
                          double sixj1 = Z.modelspace->GetSixJ( oa.j2*0.5, ob.j2*0.5, J,  oj.j2*0.5, oc.j2*0.5, JJ);
                          Gammabar_abjc -= (2*JJ+1) * sixj1 * Gamma.TwoBody.GetTBME_J(JJ,JJ,a,c,j,b);
                       }

                       
                       fIIIa_ij += 1/(oi.j2+1.0) * ( chibar_icab * Gammabar_abjc - chibar_cjab * Gammabar_abci );
                   }// for J
                }// for c
             }// for b
          }//for a
          Z.OneBody(i,j) += fIIIa_ij;
          fIIIa(i,j) = fIIIa_ij;
       }// for j
    }// for i
//       arma::mat Gamma_bar = arma::zeros(2*nkets,2*nkets);
//       for (size_t IJ=0; IJ<nkets; IJ++)
//       {
//          Ket& ketij = tbc.GetKet(IJ);
//          double ji = ketij.op->j2 * 0.5;
//          double jj = ketij.oq->j2 * 0.5;
//          for (size_t KL=0; KL<nkets; KL++)
//          {
//            Ket& ketKL = tbc.GetKet(KL);
//            double nk = ketKL.op->occ;
//            double nl = ketKL.oq->occ;
//            double jk = ketKL.op->j2 * 0.5;
//            double jl = ketKL.oq->j2 * 0.5;
//
//            double omegabar_ijkl = 0;
//            double omegabar_jikl = 0;
//            double omegabar_ijlk = 0;
//            double omegabar_jilk = 0;
//            int JJmin = AngMom::Jmin({ {ji,jl} , {jj,jk}  });
//            int JJmax = AngMom::Jmax({ {ji,jl} , {jj,jk}  });
//            for (int JJ=JJmin; JJ<=JJmax; JJ++)
//            {
//                double sixj1 = Z.modelspace->GetSixJ( ji,jj,J,  jk,jl,JJ);
//                double sixj2 = Z.modelspace->GetSixJ( jj,ji,J,  jk,jl,JJ);
//                double Omega_ilkj = Eta.TwoBody.GetTBME_J(JJ,JJ,i,l,k,j);
//                double Omega_jlki = Eta.TwoBody.GetTBME_J(JJ,JJ,j,l,k,i);
//                double Omega_iklj = Eta.TwoBody.GetTBME_J(JJ,JJ,i,k,l,j);
//                double Omega_jkli = Eta.TwoBody.GetTBME_J(JJ,JJ,j,k,l,i);
//                omegabar_ijkl += -(2*JJ+1) * sixj1 * Omega_ilkj;
//                omegabar_jikl += -(2*JJ+1) * sixj2 * Omega_jlki;
//                omegabar_ijlk += -(2*JJ+1) * sixj2 * Omega_iklj;
//                omegabar_jilk += -(2*JJ+1) * sixj1 * Omega_jkli;
//            }
//          Omega_bar(IJ,       KL    )  = omegabar_ijkl;
//          Omega_bar(IJ+nkets, KL    )  = omegabar_jikl;
//          Omega_bar(IJ,       KL+nkets) = omegabar_ijlk;
//          Omega_bar(IJ+nkets, KL+nkets) = omegabar_jilk;
//          }// for iket
//       }// for ibra



//    }//for ch
    std::cout << "fIIIa : " << std::endl << fIIIa << std::endl;


    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
    return;
  }



 //// Eq B5d, with intermediate defined in B6d.
  void comm223_231_fIIIb(const Operator &Eta, const Operator &Gamma, Operator &Z)
  {
    double t_start = omp_get_wtime();

    // The intermediate has the symmetry of Eta*Eta, so it should be a scalar, with positive parity
    // The intermediate is non-hermitian because of the occupation factors.

    TwoBodyME chi_d( Z.modelspace, 0,0,0);
    chi_d.SetNonHermitian();
    size_t nch = Z.modelspace->GetNumberTwoBodyChannels();
    for (size_t ch=0; ch<nch; ch++)
    {
       TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
       int J = tbc.J;
       size_t nkets = tbc.GetNumberKets();
       for (size_t ij=0; ij<nkets; ij++)
       {
          Ket& ketij = tbc.GetKet(ij);
            double ni = ketij.op->occ;
            double nj = ketij.oq->occ;
          for (size_t kl=0; kl<nkets; kl++)
          {
            Ket& ketkl = tbc.GetKet(kl);
            double nk = ketkl.op->occ;
            double nl = ketkl.oq->occ;

            double chi_ijkl = 0;
            // we loop over a<=b, so we need a factor 2 to account for a>b
            // the case a==b is handled automatically because we're using normalized matrix elements
            for (size_t ab=0; ab<nkets; ab++)
            {
               Ket& ketab = tbc.GetKet(ab);
               double na = ketab.op->occ;
               double nb = ketab.oq->occ;
//               double occ_factor = na*nb*(1-nk)*(1-nl) - (1-na)*(1-nb)*nk*nl;
//               double occ_factor = na*nb*(1-nk)*(1-nl) - (1-na)*(1-nb)*nk*nl;
               double occ_factor = na*nb*(1-ni)*(1-nj) - (1-na)*(1-nb)*ni*nj;
               double Omega_ijab = Eta.TwoBody.GetTBME_norm(ch,ch,ij,ab);
               double Omega_abkl = Eta.TwoBody.GetTBME_norm(ch,ch,ab,kl);
               chi_ijkl += 2 * (2*J+1) / 4.0 * occ_factor * Omega_ijab * Omega_abkl; // the extra factor of 2 is for a>b.
            }
            double normalization = 1.0;
//            if (ketij.p == ketij.q) normalization /= PhysConst::SQRT2;
//            if (ketkl.p == ketkl.q) normalization /= PhysConst::SQRT2;
            chi_d.SetTBME(ch,ch,ij,kl, chi_ijkl * normalization );
//            if ( J==0 and ketkl.p==0 and ketkl.q==0  )
//            {
//              std::cout << "      building chid   " << ketij.p << " " << ketij.q << " " << ketkl.p << " " << ketkl.q << "   " << chi_ijkl << std::endl;
//            }

          }// for iket
       }// for ibra
    }//for ch


    arma::mat fIIIb = 0*Z.OneBody;
    for ( auto i : Z.modelspace->all_orbits )
    {
      Orbit& oi = Z.modelspace->GetOrbit(i);
      for ( auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
        double fIIIb_ij = 0;
        for ( auto a : Z.modelspace->all_orbits )
        {
           Orbit& oa = Z.modelspace->GetOrbit(a);
           for ( auto b : Z.modelspace->all_orbits )
           {
              Orbit& ob = Z.modelspace->GetOrbit(b);
              for ( auto c : Z.modelspace->all_orbits )
              {
                Orbit& oc = Z.modelspace->GetOrbit(c);
                int Jmin = AngMom::Jmin({ {oa.j2,ob.j2}, {oi.j2,oc.j2}}) /2;
                int Jmax = AngMom::Jmax({ {oa.j2,ob.j2}, {oi.j2,oc.j2}}) /2;
                for (int J=Jmin; J<=Jmax; J++)
                {
//                   double chi_ciab = chi_d.GetTBME_J(J,J,c,i,a,b);
//                   double chi_cjab = chi_d.GetTBME_J(J,J,c,j,a,b);

                   double chi_abci = chi_d.GetTBME_J(J,J,a,b,c,i);
                   double chi_abcj = chi_d.GetTBME_J(J,J,a,b,c,j);
                   double Gamma_cjab = Gamma.TwoBody.GetTBME_J(J,J,c,j,a,b);
                   double Gamma_ciab = Gamma.TwoBody.GetTBME_J(J,J,c,i,a,b);
//                   double chi_abcj = chi_d.GetTBME_J(J,J,a,b,c,j);
//                   double Gamma_ciab = Gamma.TwoBody.GetTBME_J(J,J,c,j,a,b);
//                   double Gamma_abcj = Gamma.TwoBody.GetTBME_J(J,J,a,b,c,j);
//                   double Gamma_abci = Gamma.TwoBody.GetTBME_J(J,J,a,b,c,i);
////                   fIIIb_ij += 1/(oi.j2+1.0) * (chi_ciab * Gamma_abcj - chi_abcj * Gamma_ciab);
                   fIIIb_ij += 1/(oi.j2+1.0) * (chi_abcj * Gamma_ciab + chi_abci * Gamma_cjab);
//                   fIIIb_ij += 1/(oi.j2+1.0) * (chi_ciab * Gamma_abcj );
//                    if (i==0 and j==0 and std::abs(Gamma_abcj)>1e-6 and a==0 and b==0 and c==0 and J==0)
//                    {
//                        std::cout << "** abcJ = " << a << " " << b << " " << c << " " << J << " Gamma = " << Gamma_abcj << "  chi = " << chi_ciab <<  "    zij = " << fIIIb_ij << std::endl;
//                    }
                }// for J
              }// for c
           }// for b
        }// for a
        Z.OneBody(i,j) += fIIIb_ij;
        fIIIb(i,j) = fIIIb_ij;
      }// for j
    }// for i

    std::cout << "fIIIb : " << std::endl << fIIIb << std::endl;
//    std::cout << "fIIIb : " << std::endl << Z.OneBody << std::endl;

    Z.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
    return;
  }



} // namespace ReferenceImplementations
////////////////////////////////////////////////////////////////////
