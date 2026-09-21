
#include "ReferenceImplementations.hh"
#include "ModelSpace.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"
#include <deque>


namespace ReferenceImplementations
{


  double TriplesGuess(const Operator &Omega, const Operator &H)
  {
    double t_start = omp_get_wtime();
    double Etrip = 0;
 
    int norb = H.modelspace->GetNumberOrbits();
    #pragma omp parallel for collapse(2) reduction (+:Etrip)
    for (int a=0; a<norb; a++ )
    {
      for (int b=0; b<norb; b++ )
      {
        Orbit& oa = H.modelspace->GetOrbit(a);
        Orbit& ob = H.modelspace->GetOrbit(b);
        if (  ((1-oa.occ)<1e-6) or ((1-ob.occ)<1e-6) ) continue;

        for (auto c : H.modelspace->particles )
        {
          Orbit& oc = H.modelspace->GetOrbit(c);
          for (auto i : H.modelspace->holes )
          {
            Orbit& oi = H.modelspace->GetOrbit(i);
            for (auto j : H.modelspace->holes )
            {
              Orbit& oj = H.modelspace->GetOrbit(j);
              for (auto k : H.OneBodyChannels.at({oj.l,oj.j2,oj.tz2})  )
              {
                Orbit& ok = H.modelspace->GetOrbit(k);
                if (ok.occ < 1e-6) continue;
                for (auto l : H.modelspace->holes )
                {
                  Orbit& ol = H.modelspace->GetOrbit(l);
                  for (auto m : H.modelspace->holes )
                  {
                    Orbit& om = H.modelspace->GetOrbit(m);
                    double occ_factor = (1-oa.occ)*(1-ob.occ)*(1-oc.occ) * oi.occ * oj.occ * ok.occ * ol.occ * om.occ ;
                    int J1min = AngMom::Jmin({  {oa.j2,ob.j2}, {oi.j2,oj.j2} }) / 2;
                    int J1max = AngMom::Jmax({  {oa.j2,ob.j2}, {oi.j2,oj.j2} }) / 2;
                    int J2min = AngMom::Jmin({  {ol.j2,om.j2}, {oc.j2,oj.j2} }) / 2;
                    int J2max = AngMom::Jmax({  {ol.j2,om.j2}, {oc.j2,oj.j2} }) / 2;
                    double denom = H.OneBody(a,a) + H.OneBody(b,b) + H.OneBody(c,c) - H.OneBody(i,i) - H.OneBody(l,l) - H.OneBody(m,m);
                    double omega_prod = 0;
                    double h_prod = 0;
                    for (int J1=J1min; J1<=J1max; J1++)
                    {
                      omega_prod += (2*J1+1) * Omega.TwoBody.GetTBME_J(J1,J1,i,k,a,b) * Omega.TwoBody.GetTBME_J(J1,J1,a,b,i,j);
                    }
                    for (int J2=J2min; J2<=J2max; J2++)
                    {
                      h_prod += (2*J2+1) * H.TwoBody.GetTBME_J(J2,J2,l,m,k,c) * H.TwoBody.GetTBME_J(J2,J2,j,c,l,m);
                    }
                    Etrip += 0.25 * occ_factor * omega_prod * h_prod / denom / (oj.j2+1);
                  }// for m
                }// for l
              }// for k
            }// for j
          }//for i
        }// for c
      }// for b
    }// for a


    #pragma omp parallel for collapse(2) reduction (+:Etrip)
    for (int a=0; a<norb; a++ )
    {
      for (int b=0; b<norb; b++ )
      {
        Orbit& oa = H.modelspace->GetOrbit(a);
        Orbit& ob = H.modelspace->GetOrbit(b);
        if (  ((1-oa.occ)<1e-6) or ((1-ob.occ)<1e-6) ) continue;

        for (auto c : H.modelspace->particles )
        {
          Orbit& oc = H.modelspace->GetOrbit(c);
          for (auto d : H.OneBodyChannels.at({ob.l,ob.j2,ob.tz2}) )
          {
            Orbit& od = H.modelspace->GetOrbit(d);
            if ( (1-od.occ)<1e-6) continue;
            for (auto e : H.modelspace->particles )
            {
              Orbit& oe = H.modelspace->GetOrbit(e);
              for (auto i : H.modelspace->holes )
              {
                Orbit& oi = H.modelspace->GetOrbit(i);
                for (auto j : H.modelspace->holes )
                {
                  Orbit& oj = H.modelspace->GetOrbit(j);
                  for (auto k : H.modelspace->holes )
                  {
                    Orbit& ok = H.modelspace->GetOrbit(k);
                    double occ_factor = (1-oa.occ)*(1-ob.occ)*(1-oc.occ) *(1-od.occ) *(1-oe.occ) * oi.occ * oj.occ * ok.occ ;
                    int J1min = AngMom::Jmin({  {oa.j2,ob.j2}, {oi.j2,oj.j2} }) / 2;
                    int J1max = AngMom::Jmax({  {oa.j2,ob.j2}, {oi.j2,oj.j2} }) / 2;
                    int J2min = AngMom::Jmin({  {oc.j2,oe.j2}, {ob.j2,ok.j2} }) / 2;
                    int J2max = AngMom::Jmax({  {oc.j2,oe.j2}, {ob.j2,ok.j2} }) / 2;
                    double denom = H.OneBody(a,a) + H.OneBody(c,c) + H.OneBody(e,e) - H.OneBody(i,i) - H.OneBody(j,j) - H.OneBody(k,k);
                    double omega_prod = 0;
                    double h_prod = 0;
                    for (int J1=J1min; J1<=J1max; J1++)
                    {
                      omega_prod += (2*J1+1) * Omega.TwoBody.GetTBME_J(J1,J1,i,j,a,d) * Omega.TwoBody.GetTBME_J(J1,J1,a,b,i,j);
                    }
                    for (int J2=J2min; J2<=J2max; J2++)
                    {
                      h_prod += (2*J2+1) * H.TwoBody.GetTBME_J(J2,J2,d,k,c,e) * H.TwoBody.GetTBME_J(J2,J2,c,e,b,k);
                    }
                    Etrip += 0.25 * occ_factor * omega_prod * h_prod / denom / (ob.j2+1);
                  }// for m
                }// for l
              }// for k
            }// for j
          }//for i
        }// for c
      }// for b
    }// for a




   
    H.profiler.timer[ "ReferenceImplementations::" + std::string(__func__)] += omp_get_wtime() - t_start;
    return Etrip;
  }



} // namespace ReferenceImplementations
////////////////////////////////////////////////////////////////////
