
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


double GetDenom(const Operator& H, const std::vector<index_t>& holes, const std::vector<index_t>& particles)
{
   double denom = 0;
   for ( auto& h : holes )
   {
      denom += H.OneBody(h,h);
   }
   for ( auto& p : particles )
   {
      denom -= H.OneBody(p,p);
   }
   return denom;

}

double GetMP4_term( const Operator& H, int diagram)
{
   double E = 0;
   if      (diagram==1)   E = GetMP4_F1(H);
   else if (diagram==2)   E = GetMP4_F2(H);
   else if (diagram==3)   E = GetMP4_F3(H);
   else if (diagram==4)   E = GetMP4_F4(H);
   else if (diagram==5)   E = GetMP4_F5(H);
   else if (diagram==6)   E = GetMP4_F6(H);
   else if (diagram==7)   E = GetMP4_F7(H);
   else if (diagram==8)   E = GetMP4_F8(H);
   else if (diagram==9)   E = GetMP4_F9(H);
   else if (diagram==10)   E = GetMP4_F10(H);
   else if (diagram==11)   E = GetMP4_F11(H);
   else if (diagram==12)   E = GetMP4_F12(H);
   else if (diagram==13)   E = GetMP4_F13(H);
   else if (diagram==14)   E = GetMP4_F14(H);
   else if (diagram==15)   E = GetMP4_F15(H);
   else if (diagram==16)   E = GetMP4_F16(H);
   else if (diagram==17)   E = GetMP4_F17(H);
   else if (diagram==18)   E = GetMP4_F18(H);
   else if (diagram==19)   E = GetMP4_F19(H);
   else if (diagram==20)   E = GetMP4_F20(H);
   else if (diagram==21)   E = GetMP4_F21(H);
   else if (diagram==22)   E = GetMP4_F22(H);
   else if (diagram==23)   E = GetMP4_F23(H);
   else if (diagram==24)   E = GetMP4_F24(H);
   else if (diagram==25)   E = GetMP4_F25(H);
   else if (diagram==26)   E = GetMP4_F26(H);
   else if (diagram==27)   E = GetMP4_F27(H);
   else if (diagram==28)   E = GetMP4_F28(H);
   else if (diagram==29)   E = GetMP4_F29(H);
   else if (diagram==30)   E = GetMP4_F30(H);
   else if (diagram==31)   E = GetMP4_F31(H);
   else if (diagram==32)   E = GetMP4_F32(H);
   else if (diagram==33)   E = GetMP4_F33(H);
   else if (diagram==34)   E = GetMP4_F34(H);
   else if (diagram==35)   E = GetMP4_F35(H);
   else if (diagram==36)   E = GetMP4_F36(H);
   else if (diagram==37)   E = GetMP4_F37(H);
   else if (diagram==38)   E = GetMP4_F38(H);
   else if (diagram==39)   E = GetMP4_F39(H);
   else
   {
      std::cout << __func__ << " Term " << diagram << " not yet implemented." << std::endl;
   }
   return E;
}

// Diagram F1 (as numbered by ADG)   corresponds to diagram 4 from Shavitt & Bartlett
// mscheme expression: F1 = 1/4 sum_abcijklm (v_abij v_ijak v_kclm v_lmbc) / (eps_abij eps_bk eps_bclm)
// agrees.
double GetMP4_F1( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F1 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               if ( ok.j2 != ob.j2) continue;
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_kb = GetDenom(H,{k},{b});
                   double e_lmbc = GetDenom(H,{l,m},{b,c});
                   double denom = e_ijab * e_kb * e_lmbc;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{ok.j2,oc.j2},{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{ok.j2,oc.j2},{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijak = H.TwoBody.GetTBME_J(J0,J0,i,j,a,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vkclm = H.TwoBody.GetTBME_J(J1,J1,k,c,l,m);
                       double vlmbc = H.TwoBody.GetTBME_J(J1,J1,l,m,b,c);
                       F1 += 1./4 * (2*J0+1) * (2*J1+1) / (ob.j2+1) * vabij * vijak * vkclm * vlmbc / denom;
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F1;
}

// Diagram F2 (as numbered by ADG)   complex conjugate diagram: F3
// mscheme expression: F2 = -1/4 sum_abcdijkl (v_abij v_ijak v_cdbl v_klcd) / (eps_ijab eps_kb eps_klcd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F2( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F2 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 if ( ok.j2 != ob.j2) continue;
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_kb = GetDenom(H,{k},{b});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_kb * e_klcd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{ob.j2,ol.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{ob.j2,ol.j2},{ok.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijak = H.TwoBody.GetTBME_J(J0,J0,i,j,a,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdbl = H.TwoBody.GetTBME_J(J1,J1,c,d,b,l);
                       double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
                       F2 += -1./4 * (2*J0+1) * (2*J1+1) / (ob.j2+1) * vabij * vijak * vcdbl * vklcd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F2;
}



// Diagram F3 (as numbered by ADG)   complex conjugate diagram: F2
// mscheme expression: F3 = -1/4 sum_abcdijkl (v_abij v_icab v_jdkl v_klcd) / (eps_ijab eps_jc eps_klcd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F3( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F3 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               if ( oj.j2 != oc.j2) continue;
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jc = GetDenom(H,{j},{c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_jc * e_klcd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oi.j2,oc.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oi.j2,oc.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oj.j2,od.j2},{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oj.j2,od.j2},{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vicab = H.TwoBody.GetTBME_J(J0,J0,i,c,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vjdkl = H.TwoBody.GetTBME_J(J1,J1,j,d,k,l);
                       double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
                       F3 += -1./4 * (2*J0+1) * (2*J1+1) / (oc.j2+1) * vabij * vicab * vjdkl * vklcd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F3;
}



// Diagram F4 (as numbered by ADG)
// mscheme expression: F4 = 1/4 sum_abcdeijk (v_abij v_icab v_deck v_jkde) / (eps_ijab eps_jc eps_jkde)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F4( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F4 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 if ( oj.j2 != oc.j2) continue;
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jc = GetDenom(H,{j},{c});
                   double e_jkde = GetDenom(H,{j,k},{d,e});
                   double denom = e_ijab * e_jc * e_jkde;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oi.j2,oc.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oi.j2,oc.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{od.j2,oe.j2},{oc.j2,ok.j2},{oj.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{od.j2,oe.j2},{oc.j2,ok.j2},{oj.j2,ok.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vicab = H.TwoBody.GetTBME_J(J0,J0,i,c,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vdeck = H.TwoBody.GetTBME_J(J1,J1,d,e,c,k);
                       double vjkde = H.TwoBody.GetTBME_J(J1,J1,j,k,d,e);
                       F4 += 1./4 * (2*J0+1) * (2*J1+1) / (oc.j2+1) * vabij * vicab * vdeck * vjkde / denom;
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F4;
}


// Diagram F5 (as numbered by ADG)
// mscheme expression: F5 = 1/16 sum_abijklmn (v_abij v_ijkl v_klmn v_mnab) / (eps^ij_ab eps^kl_ab eps^mn_ab)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F5( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F5 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto i: H.modelspace->holes )
       {
         Orbit& oi = H.modelspace->GetOrbit(i);
         for ( auto j: H.modelspace->holes )
         {
           Orbit& oj = H.modelspace->GetOrbit(j);
           for ( auto k: H.modelspace->holes )
           {
             Orbit& ok = H.modelspace->GetOrbit(k);
             for ( auto l: H.modelspace->holes )
             {
               Orbit& ol = H.modelspace->GetOrbit(l);
               for ( auto m: H.modelspace->holes )
               {
                 Orbit& om = H.modelspace->GetOrbit(m);
                 for ( auto n: H.modelspace->holes )
                 {
                   Orbit& on = H.modelspace->GetOrbit(n);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double e_mnab = GetDenom(H,{m,n},{a,b});
                   double denom = e_ijab * e_klab * e_mnab;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ok.j2,ol.j2},{om.j2,on.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ok.j2,ol.j2},{om.j2,on.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijkl = H.TwoBody.GetTBME_J(J0,J0,i,j,k,l);
                     double vklmn = H.TwoBody.GetTBME_J(J0,J0,k,l,m,n);
                     double vmnab = H.TwoBody.GetTBME_J(J0,J0,m,n,a,b);
                     F5 += 1./16 * (2*J0+1) * vabij * vijkl * vklmn * vmnab / denom;
                   }// for J0
                 }// for n
               }// for m
             }// for l
           }// for k
         }// for j
       }// for i
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F5;
}



// Diagram F6 (as numbered by ADG)   complex conjugate diagram: F8
// mscheme expression: F6 = 1/2 sum_abcijklm (v_abij v_ijkl v_kcam v_lmbc) / (eps^ij_ab eps^kl_ab eps^lm_bc)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F6( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F6 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double e_lmbc = GetDenom(H,{l,m},{b,c});
                   double denom = e_ijab * e_klab * e_lmbc;
                   int phase_exponent = (ob.j2+ol.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ok.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{ok.j2,oc.j2},{oa.j2,om.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{ok.j2,oc.j2},{oa.j2,om.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijkl = H.TwoBody.GetTBME_J(J0,J0,i,j,k,l);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vkcam = H.TwoBody.GetTBME_J(J1,J1,k,c,a,m);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vlmbc = H.TwoBody.GetTBME_J(J2,J2,l,m,b,c);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, ok.j2/2., J2, om.j2/2., ol.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F6 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vijkl * vkcam * vlmbc / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F6;
}



// Diagram F7 (as numbered by ADG)   complex conjugate diagram: F14
// mscheme expression: F7 = 1/16 sum_abcdijkl (v_abij v_ijkl v_cdab v_klcd) / (eps^ij_ab eps^kl_ab eps^kl_cd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F7( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F7 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_klab * e_klcd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijkl = H.TwoBody.GetTBME_J(J0,J0,i,j,k,l);
                     double vcdab = H.TwoBody.GetTBME_J(J0,J0,c,d,a,b);
                     double vklcd = H.TwoBody.GetTBME_J(J0,J0,k,l,c,d);
                     F7 += 1./16 * (2*J0+1) * vabij * vijkl * vcdab * vklcd / denom;
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F7;
}



// Diagram F8 (as numbered by ADG)   complex conjugate diagram: F6
// mscheme expression: F8 = 1/2 sum_abcijklm (v_abij v_icak v_jklm v_lmbc) / (eps^ij_ab eps^jk_bc eps^lm_bc)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F8( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F8 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jkbc = GetDenom(H,{j,k},{b,c});
                   double e_lmbc = GetDenom(H,{l,m},{b,c});
                   double denom = e_ijab * e_jkbc * e_lmbc;
                   int phase_exponent = (ob.j2+oj.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,ok.j2},{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,ok.j2},{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicak = H.TwoBody.GetTBME_J(J1,J1,i,c,a,k);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjklm = H.TwoBody.GetTBME_J(J2,J2,j,k,l,m);
                         double vlmbc = H.TwoBody.GetTBME_J(J2,J2,l,m,b,c);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, oi.j2/2., J2, ok.j2/2., oj.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F8 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vicak * vjklm * vlmbc / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F8;
}



// Diagram F9 (as numbered by ADG)
// mscheme expression: F9 = sum_abcdijkl (v_abik v_icaj v_jdcl v_klbd) / (eps^ik_ab eps^jk_cb eps^kl_bd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F9( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F9 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F9)
   for ( auto a : part_vec )
   {
     for ( auto b: part_vec )
     {
       for ( auto c: part_vec )
       {
     Orbit& oa = H.modelspace->GetOrbit(a);
       Orbit& ob = H.modelspace->GetOrbit(b);
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_jkcb = GetDenom(H,{j,k},{c,b});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ikab * e_jkcb * e_klbd;
                   int phase_exponent = (ob.j2+oc.j2+oj.j2+ok.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,oj.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,oj.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{oc.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{oc.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   // J4 does not label any two-body matrix element; its range comes purely from
                   // the triangle conditions of the four 6j symbols it appears in.
                   int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{oc.j2,oj.j2},{od.j2,ol.j2}} ) /2;
                   int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{oc.j2,oj.j2},{od.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicaj = H.TwoBody.GetTBME_J(J1,J1,i,c,a,j);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdcl = H.TwoBody.GetTBME_J(J2,J2,j,d,c,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oc.j2/2., oj.j2/2., J4, oa.j2/2., oi.j2/2., J1 );
                             double sixj3 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, ok.j2/2., ob.j2/2., J3 );
                             double sixj4 = H.modelspace->GetSixJ( ol.j2/2., od.j2/2., J4, oj.j2/2., oc.j2/2., J2 );
                             F9 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                   * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vicaj * vjdcl * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F9;
}



// Diagram F10 (as numbered by ADG)
// mscheme expression: F10 = -sum_abcdijkl (v_abij v_icak v_jdcl v_klbd) / (eps^ij_ab eps^jk_cb eps^kl_bd)
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F10( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F10 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F10)
//   #pragma omp parallel for
   for ( auto a : part_vec )
   {
     for ( auto b: part_vec )
     {
       for ( auto c: part_vec )
       {
     Orbit& oa = H.modelspace->GetOrbit(a);
       Orbit& ob = H.modelspace->GetOrbit(b);
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jkcb = GetDenom(H,{j,k},{c,b});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ijab * e_jkcb * e_klbd;
                   int phase_exponent = (od.j2+ol.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{oc.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{oc.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicak = H.TwoBody.GetTBME_J(J1,J1,i,c,a,k);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdcl = H.TwoBody.GetTBME_J(J2,J2,j,d,c,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           // J4 does not label any two-body matrix element. Two of the 6j symbols
                           // couple it to single-particle j's, but the other two couple it to J2
                           // and J3 directly, so its range can only be fixed once J2,J3 are known.
                           int J4_min = AngMom::Jmin( {{oj.j2,ob.j2},{oc.j2,ok.j2},{2*J2,2*J3}} ) /2;
                           int J4_max = AngMom::Jmax( {{oj.j2,ob.j2},{oc.j2,ok.j2},{2*J2,2*J3}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( oj.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oc.j2/2., ok.j2/2., J4, oa.j2/2., oi.j2/2., J1 );
                             double sixj3 = H.modelspace->GetSixJ( J3, J4, J2, oj.j2/2., od.j2/2., ob.j2/2. );
                             double sixj4 = H.modelspace->GetSixJ( J4, J2, J3, ol.j2/2., ok.j2/2., oc.j2/2. );
                             F10 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabij * vicak * vjdcl * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F10;
}



// Diagram F11 (as numbered by ADG)
// mscheme expression: F11 = -sum_abcdijkl (v_abik v_icaj v_jdbl v_klcd) / (eps^ik_ab eps^jk_bc eps^kl_cd)
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F11( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F11 =0;
//   #pragma omp parallel for
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F11)
   for ( auto a : part_vec )
   {
     for ( auto b: part_vec )
     {
       for ( auto c: part_vec )
       {
     Orbit& oa = H.modelspace->GetOrbit(a);
       Orbit& ob = H.modelspace->GetOrbit(b);
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_jkbc = GetDenom(H,{j,k},{b,c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ikab * e_jkbc * e_klcd;
                   int phase_exponent = (od.j2+ol.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,oj.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,oj.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{ob.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{ob.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicaj = H.TwoBody.GetTBME_J(J1,J1,i,c,a,j);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdbl = H.TwoBody.GetTBME_J(J2,J2,j,d,b,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklcd = H.TwoBody.GetTBME_J(J3,J3,k,l,c,d);
                           // As in F10, J4 is tied to J2,J3 directly through two of the 6j symbols.
                           int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{oc.j2,oj.j2},{2*J2,2*J3}} ) /2;
                           int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{oc.j2,oj.j2},{2*J2,2*J3}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oc.j2/2., oj.j2/2., J4, oa.j2/2., oi.j2/2., J1 );
                             double sixj3 = H.modelspace->GetSixJ( J3, J4, J2, ob.j2/2., ol.j2/2., ok.j2/2. );
                             double sixj4 = H.modelspace->GetSixJ( J3, J2, J4, oj.j2/2., oc.j2/2., od.j2/2. );
                             F11 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vicaj * vjdbl * vklcd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F11;
}



// Diagram F12 (as numbered by ADG)
// mscheme expression: F12 = sum_abcdijkl (v_abij v_icak v_jdbl v_klcd) / (eps^ij_ab eps^jk_bc eps^kl_cd)
// Fixed the triangle condidions. Now it works.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F12( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F12 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jkbc = GetDenom(H,{j,k},{b,c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_jkbc * e_klcd;
                   int phase_exponent = (ob.j2+oc.j2+oj.j2+ok.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{ob.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{ob.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   // Here all four 6j symbols tie J4 to single-particle j's only.
                   // SRS: J4 has triangle cond with (j,b) (a,i) (c,k) (l,d).  It appears Claude messed up the conditions
//                   int J4_min = AngMom::Jmin( {{oj.j2,ob.j2},{oc.j2,ok.j2},{ol.j2,od.j2},{ol.j2,ok.j2}} ) /2;
//                   int J4_max = AngMom::Jmax( {{oj.j2,ob.j2},{oc.j2,ok.j2},{ol.j2,od.j2},{ol.j2,ok.j2}} ) /2;
                   int J4_min = AngMom::Jmin( {{oj.j2,ob.j2},{oc.j2,ok.j2},{ol.j2,od.j2},{oa.j2,oi.j2}} ) /2;
                   int J4_max = AngMom::Jmax( {{oj.j2,ob.j2},{oc.j2,ok.j2},{ol.j2,od.j2},{oa.j2,oi.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicak = H.TwoBody.GetTBME_J(J1,J1,i,c,a,k);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdbl = H.TwoBody.GetTBME_J(J2,J2,j,d,b,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklcd = H.TwoBody.GetTBME_J(J3,J3,k,l,c,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( oj.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oc.j2/2., ok.j2/2., J4, oa.j2/2., oi.j2/2., J1 );
                             double sixj3 = H.modelspace->GetSixJ( ol.j2/2., od.j2/2., J4, oj.j2/2., ob.j2/2., J2 );
                             double sixj4 = H.modelspace->GetSixJ( oc.j2/2., od.j2/2., J3, ol.j2/2., ok.j2/2., J4 );
                             F12 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabij * vicak * vjdbl * vklcd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F12;
}



// Diagram F13 (as numbered by ADG)   complex conjugate diagram: F15
// mscheme expression: F13 = 1/2 sum_abcdeijk (v_abij v_icak v_debc v_jkde) / (eps^ij_ab eps^jk_bc eps^jk_de)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F13( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F13 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jkbc = GetDenom(H,{j,k},{b,c});
                   double e_jkde = GetDenom(H,{j,k},{d,e});
                   double denom = e_ijab * e_jkbc * e_jkde;
                   int phase_exponent = (ob.j2+oj.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{od.j2,oe.j2},{ob.j2,oc.j2},{oj.j2,ok.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{od.j2,oe.j2},{ob.j2,oc.j2},{oj.j2,ok.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicak = H.TwoBody.GetTBME_J(J1,J1,i,c,a,k);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vdebc = H.TwoBody.GetTBME_J(J2,J2,d,e,b,c);
                         double vjkde = H.TwoBody.GetTBME_J(J2,J2,j,k,d,e);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, oi.j2/2., J2, ok.j2/2., oj.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F13 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vicak * vdebc * vjkde / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F13;
}



// Diagram F14 (as numbered by ADG)   complex conjugate diagram: F7
// mscheme expression: F14 = 1/16 sum_abcdijkl (v_abij v_cdab v_ijkl v_klcd) / (eps^ij_ab eps^ij_cd eps^kl_cd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F14( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F14 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijcd = GetDenom(H,{i,j},{c,d});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_ijcd * e_klcd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2},{ok.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2},{ok.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vcdab = H.TwoBody.GetTBME_J(J0,J0,c,d,a,b);
                     double vijkl = H.TwoBody.GetTBME_J(J0,J0,i,j,k,l);
                     double vklcd = H.TwoBody.GetTBME_J(J0,J0,k,l,c,d);
                     F14 += 1./16 * (2*J0+1) * vabij * vcdab * vijkl * vklcd / denom;
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F14;
}



// Diagram F15 (as numbered by ADG)   complex conjugate diagram: F13
// mscheme expression: F15 = 1/2 sum_abcdeijk (v_abij v_cdab v_ieck v_jkde) / (eps^ij_ab eps^ij_cd eps^jk_de)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F15( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F15 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijcd = GetDenom(H,{i,j},{c,d});
                   double e_jkde = GetDenom(H,{j,k},{d,e});
                   double denom = e_ijab * e_ijcd * e_jkde;
                   int phase_exponent = (od.j2+oj.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vcdab = H.TwoBody.GetTBME_J(J0,J0,c,d,a,b);
                     int J1_min = AngMom::Jmin( {{oi.j2,oe.j2},{oc.j2,ok.j2}} ) /2;
                     int J1_max = AngMom::Jmax( {{oi.j2,oe.j2},{oc.j2,ok.j2}} ) /2;
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vieck = H.TwoBody.GetTBME_J(J1,J1,i,e,c,k);
                       int J2_min = AngMom::Jmin( {{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                       int J2_max = AngMom::Jmax( {{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjkde = H.TwoBody.GetTBME_J(J2,J2,j,k,d,e);
                         double ninej = H.modelspace->GetNineJ( oj.j2/2., oi.j2/2., J0, ok.j2/2., J1, oc.j2/2., J2, oe.j2/2., od.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F15 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vcdab * vieck * vjkde / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F15;
}



// Diagram F16 (as numbered by ADG)
// mscheme expression: F16 = 1/16 sum_abcdefij (v_abij v_cdab v_efcd v_ijef) / (eps^ij_ab eps^ij_cd eps^ij_ef)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F16( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F16 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto f: H.modelspace->particles )
             {
               Orbit& of = H.modelspace->GetOrbit(f);
               for ( auto i: H.modelspace->holes )
               {
                 Orbit& oi = H.modelspace->GetOrbit(i);
                 for ( auto j: H.modelspace->holes )
                 {
                   Orbit& oj = H.modelspace->GetOrbit(j);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijcd = GetDenom(H,{i,j},{c,d});
                   double e_ijef = GetDenom(H,{i,j},{e,f});
                   double denom = e_ijab * e_ijcd * e_ijef;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2},{oe.j2,of.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2},{oe.j2,of.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vcdab = H.TwoBody.GetTBME_J(J0,J0,c,d,a,b);
                     double vefcd = H.TwoBody.GetTBME_J(J0,J0,e,f,c,d);
                     double vijef = H.TwoBody.GetTBME_J(J0,J0,i,j,e,f);
                     F16 += 1./16 * (2*J0+1) * vabij * vcdab * vefcd * vijef / denom;
                   }// for J0
                 }// for j
               }// for i
             }//for f
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F16;
}



// Diagram F17 (as numbered by ADG)
// mscheme expression: F17 = 1/4 sum_abcijklm (v_abil v_icjk v_jkcm v_lmab) / (eps^il_ab eps^jkl_cab eps^lm_ab)
// Claude was OFF BY A MINUS SIGN. I fixed it.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F17( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F17 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   if ( om.j2 != oi.j2) continue;
                   double e_ilab = GetDenom(H,{i,l},{a,b});
                   double e_jklcab = GetDenom(H,{j,k,l},{c,a,b});
                   double e_lmab = GetDenom(H,{l,m},{a,b});
                   double denom = e_ilab * e_jklcab * e_lmab;
                   int phase_exponent = (oc.j2+ol.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ol.j2},{ol.j2,om.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ol.j2},{ol.j2,om.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,ok.j2},{oc.j2,om.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,ok.j2},{oc.j2,om.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabil = H.TwoBody.GetTBME_J(J0,J0,a,b,i,l);
                     double vlmab = H.TwoBody.GetTBME_J(J0,J0,l,m,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicjk = H.TwoBody.GetTBME_J(J1,J1,i,c,j,k);
                       double vjkcm = H.TwoBody.GetTBME_J(J1,J1,j,k,c,m);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F17 -= 1./4 * phase * (2*J0+1) * (2*J1+1) / (oi.j2+1) * vabil * vicjk * vjkcm * vlmab / denom;
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F17;
}



// Diagram F18 (as numbered by ADG)
// mscheme expression: F18 = 1/2 sum_abcijklm (v_abij v_ickl v_jkcm v_lmab) / (eps^ij_ab eps^jkl_cab eps^lm_ab)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F18( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F18 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jklcab = GetDenom(H,{j,k,l},{c,a,b});
                   double e_lmab = GetDenom(H,{l,m},{a,b});
                   double denom = e_ijab * e_jklcab * e_lmab;
                   int phase_exponent = (oc.j2+oj.j2+ok.j2+ol.j2)/2;  // + J0, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ol.j2,om.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ol.j2,om.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{ok.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,ok.j2},{oc.j2,om.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,ok.j2},{oc.j2,om.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vlmab = H.TwoBody.GetTBME_J(J0,J0,l,m,a,b);
                     double phase = ( (J0+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vickl = H.TwoBody.GetTBME_J(J1,J1,i,c,k,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjkcm = H.TwoBody.GetTBME_J(J2,J2,j,k,c,m);
                         double ninej = H.modelspace->GetNineJ( oj.j2/2., oi.j2/2., J0, ok.j2/2., J1, ol.j2/2., J2, oc.j2/2., om.j2/2. );
                         F18 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vickl * vjkcm * vlmab / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F18;
}


// Diagram F19 (as numbered by ADG)
// mscheme expression: F19 = 1/2 sum_abcijklm (v_abil v_icjk v_jkam v_lmbc) / (eps^il_ab eps^jkl_abc eps^lm_bc)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F19( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F19 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ilab = GetDenom(H,{i,l},{a,b});
                   double e_jklabc = GetDenom(H,{j,k,l},{a,b,c});
                   double e_lmbc = GetDenom(H,{l,m},{b,c});
                   double denom = e_ilab * e_jklabc * e_lmbc;
                   int phase_exponent = (ob.j2+ol.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ol.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,ok.j2},{oa.j2,om.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,ok.j2},{oa.j2,om.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabil = H.TwoBody.GetTBME_J(J0,J0,a,b,i,l);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicjk = H.TwoBody.GetTBME_J(J1,J1,i,c,j,k);
                       double vjkam = H.TwoBody.GetTBME_J(J1,J1,j,k,a,m);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vlmbc = H.TwoBody.GetTBME_J(J2,J2,l,m,b,c);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, oi.j2/2., J2, om.j2/2., ol.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F19 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabil * vicjk * vjkam * vlmbc / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F19;
}



// Diagram F20 (as numbered by ADG)
// mscheme expression: F20 = sum_abcijklm (v_abij v_ickl v_jkam v_lmbc) / (eps^ij_ab eps^jkl_abc eps^lm_bc)
// There was a bug in Takayuki's Fortran implementation. Now both agree.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F20( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F20 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F20)
   for ( auto a : part_vec )
   {
     for ( auto b: part_vec )
     {
       for ( auto c: part_vec )
       {
     Orbit& oa = H.modelspace->GetOrbit(a);
       Orbit& ob = H.modelspace->GetOrbit(b);
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jklabc = GetDenom(H,{j,k,l},{a,b,c});
                   double e_lmbc = GetDenom(H,{l,m},{b,c});
                   double denom = e_ijab * e_jklabc * e_lmbc;
                   int phase_exponent = (oc.j2+oi.j2+oj.j2+om.j2)/2;  // + J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{ok.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,ok.j2},{oa.j2,om.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,ok.j2},{oa.j2,om.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vickl = H.TwoBody.GetTBME_J(J1,J1,i,c,k,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjkam = H.TwoBody.GetTBME_J(J2,J2,j,k,a,m);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vlmbc = H.TwoBody.GetTBME_J(J3,J3,l,m,b,c);
                           double phase = ( (J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           // J4 does not label any two-body matrix element. Two of the 6j symbols
                           // couple it to single-particle j's; the other two couple it to J1 and J3.
                           int J4_min = AngMom::Jmin( {{oi.j2,ob.j2},{oj.j2,oa.j2},{om.j2,ok.j2},{2*J1,2*J3}} ) /2;
                           int J4_max = AngMom::Jmax( {{oi.j2,ob.j2},{oj.j2,oa.j2},{om.j2,ok.j2},{2*J1,2*J3}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( oj.j2/2., oa.j2/2., J4, ob.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( J1, J3, J4, ob.j2/2., oi.j2/2., oc.j2/2. );
                             double sixj3 = H.modelspace->GetSixJ( om.j2/2., ok.j2/2., J4, oj.j2/2., oa.j2/2., J2 );
                             double sixj4 = H.modelspace->GetSixJ( J1, J3, J4, om.j2/2., ok.j2/2., ol.j2/2. );
                             F20 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabij * vickl * vjkam * vlmbc / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F20;
}



// Diagram F21 (as numbered by ADG)   complex conjugate diagram: F25
// mscheme expression: F21 = -sum_abcdijkl (v_abik v_icjl v_jdac v_klbd) / (eps^ik_ab eps^jkl_acb eps^kl_bd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F21( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F21 =0;
//   #pragma omp parallel for
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F21)
   for ( auto a : part_vec )
   {
     for ( auto b: part_vec )
     {
       for ( auto c: part_vec )
       {
     Orbit& oa = H.modelspace->GetOrbit(a);
       Orbit& ob = H.modelspace->GetOrbit(b);
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_jklacb = GetDenom(H,{j,k,l},{a,c,b});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ikab * e_jklacb * e_klbd;
                   int phase_exponent = (oa.j2+ob.j2+oc.j2+od.j2+oi.j2+oj.j2+ok.j2+ol.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{oa.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{oa.j2,oc.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicjl = H.TwoBody.GetTBME_J(J1,J1,i,c,j,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdac = H.TwoBody.GetTBME_J(J2,J2,j,d,a,c);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           // J4 is tied to two single-particle pairs and to J1,J2 directly.
                           int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{od.j2,ol.j2},{2*J1,2*J2}} ) /2;
                           int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{od.j2,ol.j2},{2*J1,2*J2}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( J1, J2, J4, oa.j2/2., oi.j2/2., oc.j2/2. );
                             double sixj3 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, ok.j2/2., ob.j2/2., J3 );
                             double sixj4 = H.modelspace->GetSixJ( J1, J4, J2, od.j2/2., oj.j2/2., ol.j2/2. );
                             F21 += -phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vicjl * vjdac * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F21;
}



// Diagram F22 (as numbered by ADG)   complex conjugate diagram: F26
// mscheme expression: F22 = 1/2 sum_abcdijkl (v_abij v_ickl v_jdac v_klbd) / (eps^ij_ab eps^jkl_acb eps^kl_bd)
// Missing minus sign? Yes, dropped the minus sign from AMC
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F22( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F22 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jklacb = GetDenom(H,{j,k,l},{a,c,b});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ijab * e_jklacb * e_klbd;
                   int phase_exponent = (ob.j2+oj.j2)/2;  // + J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{oa.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{oa.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vickl = H.TwoBody.GetTBME_J(J1,J1,i,c,k,l);
                       double vklbd = H.TwoBody.GetTBME_J(J1,J1,k,l,b,d);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdac = H.TwoBody.GetTBME_J(J2,J2,j,d,a,c);
                         double ninej = H.modelspace->GetNineJ( oa.j2/2., ob.j2/2., J0, oc.j2/2., J1, oi.j2/2., J2, od.j2/2., oj.j2/2. );
                         double phase = ( (J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F22 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vickl * vjdac * vklbd / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F22;
}



// Diagram F23 (as numbered by ADG)   complex conjugate diagram: F29
// mscheme expression: F23 = 1/2 sum_abcdijkl (v_abik v_icjl v_jdab v_klcd) / (eps^ik_ab eps^jkl_abc eps^kl_cd)
// missing minus sign? yep. fixed it.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F23( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F23 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_jklabc = GetDenom(H,{j,k,l},{a,b,c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ikab * e_jklabc * e_klcd;
                   int phase_exponent = (oc.j2+ok.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2},{oj.j2,od.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2},{oj.j2,od.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     double vjdab = H.TwoBody.GetTBME_J(J0,J0,j,d,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicjl = H.TwoBody.GetTBME_J(J1,J1,i,c,j,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vklcd = H.TwoBody.GetTBME_J(J2,J2,k,l,c,d);
                         double ninej = H.modelspace->GetNineJ( ok.j2/2., oi.j2/2., J0, ol.j2/2., J1, oj.j2/2., J2, oc.j2/2., od.j2/2. );
                         double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F23 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabik * vicjl * vjdab * vklcd / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F23;
}



// Diagram F24 (as numbered by ADG)   complex conjugate diagram: F30
// mscheme expression: F24 = -1/4 sum_abcdijkl (v_abij v_ickl v_jdab v_klcd) / (eps^ij_ab eps^jkl_abc eps^kl_cd)
// missing minus sign. fixed it.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F24( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F24 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             if ( oi.j2 != od.j2) continue;
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jklabc = GetDenom(H,{j,k,l},{a,b,c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_jklabc * e_klcd;
                   int phase_exponent = (oc.j2+oj.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oj.j2,od.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oj.j2,od.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vjdab = H.TwoBody.GetTBME_J(J0,J0,j,d,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vickl = H.TwoBody.GetTBME_J(J1,J1,i,c,k,l);
                       double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F24 += 1./4 * phase * (2*J0+1) * (2*J1+1) / (od.j2+1) * vabij * vickl * vjdab * vklcd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F24;
}



// Diagram F25 (as numbered by ADG)   complex conjugate diagram: F21
// mscheme expression: F25 = -sum_abcdijkl (v_abik v_cdaj v_ijcl v_klbd) / (eps^ik_ab eps^ijk_cbd eps^kl_bd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F25( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F25 =0;
//   #pragma omp parallel for
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F25)
   for ( auto a : part_vec )
   {
     for ( auto b: part_vec )
     {
       for ( auto c: part_vec )
       {
     Orbit& oa = H.modelspace->GetOrbit(a);
       Orbit& ob = H.modelspace->GetOrbit(b);
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_ijkcbd = GetDenom(H,{i,j,k},{c,b,d});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ikab * e_ijkcbd * e_klbd;
                   int phase_exponent = (oa.j2+ob.j2+oc.j2+od.j2+oi.j2+oj.j2+ok.j2+ol.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,oj.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,oj.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oi.j2,oj.j2},{oc.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oi.j2,oj.j2},{oc.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdaj = H.TwoBody.GetTBME_J(J1,J1,c,d,a,j);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vijcl = H.TwoBody.GetTBME_J(J2,J2,i,j,c,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           // J4 is tied to two single-particle pairs and to J1,J2 directly.
                           int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{od.j2,ol.j2},{2*J1,2*J2}} ) /2;
                           int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{od.j2,ol.j2},{2*J1,2*J2}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( J2, J1, J4, oa.j2/2., oi.j2/2., oj.j2/2. );
                             double sixj3 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, ok.j2/2., ob.j2/2., J3 );
                             double sixj4 = H.modelspace->GetSixJ( J1, J2, J4, ol.j2/2., od.j2/2., oc.j2/2. );
                             F25 += -phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vcdaj * vijcl * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F25;
}



// Diagram F26 (as numbered by ADG)   complex conjugate diagram: F22
// mscheme expression: F26 = 1/2 sum_abcdijkl (v_abij v_cdak v_ijcl v_klbd) / (eps^ij_ab eps^ijk_cbd eps^kl_bd)
// Missing minus sign due to bug in AMC
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F26( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F26 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijkcbd = GetDenom(H,{i,j,k},{c,b,d});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ijab * e_ijkcbd * e_klbd;
                   int phase_exponent = (ob.j2+ok.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,ol.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijcl = H.TwoBody.GetTBME_J(J0,J0,i,j,c,l);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdak = H.TwoBody.GetTBME_J(J1,J1,c,d,a,k);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vklbd = H.TwoBody.GetTBME_J(J2,J2,k,l,b,d);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, od.j2/2., J1, oc.j2/2., J2, ok.j2/2., ol.j2/2. );
                         F26 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vcdak * vijcl * vklbd / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F26;
}



// Diagram F27 (as numbered by ADG)
// mscheme expression: F27 = 1/4 sum_abcdeijk (v_abjk v_cdai v_iecd v_jkbe) / (eps^jk_ab eps^ijk_cdb eps^jk_be)
// missing minus sign?  yep. fixed it.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F27( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F27 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             if ( oe.j2 != oa.j2) continue;
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_jkab = GetDenom(H,{j,k},{a,b});
                   double e_ijkcdb = GetDenom(H,{i,j,k},{c,d,b});
                   double e_jkbe = GetDenom(H,{j,k},{b,e});
                   double denom = e_jkab * e_ijkcdb * e_jkbe;
                   int phase_exponent = (ob.j2+oi.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oj.j2,ok.j2},{ob.j2,oe.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oj.j2,ok.j2},{ob.j2,oe.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,oi.j2},{oi.j2,oe.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,oi.j2},{oi.j2,oe.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabjk = H.TwoBody.GetTBME_J(J0,J0,a,b,j,k);
                     double vjkbe = H.TwoBody.GetTBME_J(J0,J0,j,k,b,e);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdai = H.TwoBody.GetTBME_J(J1,J1,c,d,a,i);
                       double viecd = H.TwoBody.GetTBME_J(J1,J1,i,e,c,d);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F27 -= 1./4 * phase * (2*J0+1) * (2*J1+1) / (oa.j2+1) * vabjk * vcdai * viecd * vjkbe / denom;
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F27;
}



// Diagram F28 (as numbered by ADG)
// mscheme expression: F28 = 1/2 sum_abcdeijk (v_abij v_cdak v_iecd v_jkbe) / (eps^ij_ab eps^ijk_cdb eps^jk_be)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F28( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F28 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijkcdb = GetDenom(H,{i,j,k},{c,d,b});
                   double e_jkbe = GetDenom(H,{j,k},{b,e});
                   double denom = e_ijab * e_ijkcdb * e_jkbe;
                   int phase_exponent = (ob.j2+oj.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,ok.j2},{oi.j2,oe.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,ok.j2},{oi.j2,oe.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,ok.j2},{ob.j2,oe.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,ok.j2},{ob.j2,oe.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdak = H.TwoBody.GetTBME_J(J1,J1,c,d,a,k);
                       double viecd = H.TwoBody.GetTBME_J(J1,J1,i,e,c,d);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjkbe = H.TwoBody.GetTBME_J(J2,J2,j,k,b,e);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oe.j2/2., J1, oi.j2/2., J2, ok.j2/2., oj.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F28 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vcdak * viecd * vjkbe / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F28;
}



// Diagram F29 (as numbered by ADG)   complex conjugate diagram: F23
// mscheme expression: F29 = 1/2 sum_abcdijkl (v_abik v_cdaj v_ijbl v_klcd) / (eps^ik_ab eps^ijk_bcd eps^kl_cd)
// minus sign was missing. fixed it.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F29( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F29 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_ijkbcd = GetDenom(H,{i,j,k},{b,c,d});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ikab * e_ijkbcd * e_klcd;
                   int phase_exponent = (ob.j2+ok.j2)/2;  // + J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,oj.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,oj.j2},{ok.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oi.j2,oj.j2},{ob.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oi.j2,oj.j2},{ob.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdaj = H.TwoBody.GetTBME_J(J1,J1,c,d,a,j);
                       double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vijbl = H.TwoBody.GetTBME_J(J2,J2,i,j,b,l);
                         double ninej = H.modelspace->GetNineJ( oa.j2/2., ob.j2/2., J0, oj.j2/2., J2, oi.j2/2., J1, ol.j2/2., ok.j2/2. );
                         double phase = ( (J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F29 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabik * vcdaj * vijbl * vklcd / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F29;
}



// Diagram F30 (as numbered by ADG)   complex conjugate diagram: F24
// mscheme expression: F30 = -1/4 sum_abcdijkl (v_abij v_cdak v_ijbl v_klcd) / (eps^ij_ab eps^ijk_bcd eps^kl_cd)
// missing minus sign? yep. fixed it.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F30( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F30 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   if ( ol.j2 != oa.j2) continue;
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijkbcd = GetDenom(H,{i,j,k},{b,c,d});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_ijkbcd * e_klcd;
                   int phase_exponent = (ob.j2+ok.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ob.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ob.j2,ol.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,ok.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,ok.j2},{ok.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijbl = H.TwoBody.GetTBME_J(J0,J0,i,j,b,l);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdak = H.TwoBody.GetTBME_J(J1,J1,c,d,a,k);
                       double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F30 += 1./4 * phase * (2*J0+1) * (2*J1+1) / (oa.j2+1) * vabij * vcdak * vijbl * vklcd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F30;
}



// Diagram F31 (as numbered by ADG)
// mscheme expression: F31 = 1/2 sum_abcdeijk (v_abjk v_cdai v_iebc v_jkde) / (eps^jk_ab eps^ijk_bcd eps^jk_de)
// missing minus sign. fixed it.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F31( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F31 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_jkab = GetDenom(H,{j,k},{a,b});
                   double e_ijkbcd = GetDenom(H,{i,j,k},{b,c,d});
                   double e_jkde = GetDenom(H,{j,k},{d,e});
                   double denom = e_jkab * e_ijkbcd * e_jkde;
                   int phase_exponent = (ob.j2+oc.j2+od.j2+oi.j2)/2;  // + J0, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,oi.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,oi.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oi.j2,oe.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oi.j2,oe.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabjk = H.TwoBody.GetTBME_J(J0,J0,a,b,j,k);
                     double vjkde = H.TwoBody.GetTBME_J(J0,J0,j,k,d,e);
                     double phase = ( (J0+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdai = H.TwoBody.GetTBME_J(J1,J1,c,d,a,i);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double viebc = H.TwoBody.GetTBME_J(J2,J2,i,e,b,c);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, od.j2/2., J2, oi.j2/2., oe.j2/2. );
                         F31 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabjk * vcdai * viebc * vjkde / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F31;
}



// Diagram F32 (as numbered by ADG)
// mscheme expression: F32 = sum_abcdeijk (v_abij v_cdak v_iebc v_jkde) / (eps^ij_ab eps^ijk_bcd eps^jk_de)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F32( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F32 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F32)
   for ( auto a : part_vec )
   {
     for ( auto b: part_vec )
     {
       for ( auto c: part_vec )
       {
     Orbit& oa = H.modelspace->GetOrbit(a);
       Orbit& ob = H.modelspace->GetOrbit(b);
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijkbcd = GetDenom(H,{i,j,k},{b,c,d});
                   double e_jkde = GetDenom(H,{j,k},{d,e});
                   double denom = e_ijab * e_ijkbcd * e_jkde;
                   int phase_exponent = (oa.j2+ob.j2+oe.j2+ok.j2)/2;  // + J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oi.j2,oe.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oi.j2,oe.j2},{ob.j2,oc.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdak = H.TwoBody.GetTBME_J(J1,J1,c,d,a,k);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double viebc = H.TwoBody.GetTBME_J(J2,J2,i,e,b,c);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vjkde = H.TwoBody.GetTBME_J(J3,J3,j,k,d,e);
                           double phase = ( (J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           // J4 is tied to two single-particle pairs and to J1,J3 directly.
                           int J4_min = AngMom::Jmin( {{oi.j2,ob.j2},{oc.j2,oe.j2},{2*J3,2*J1}} ) /2;
                           int J4_max = AngMom::Jmax( {{oi.j2,ob.j2},{oc.j2,oe.j2},{2*J3,2*J1}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( oi.j2/2., ob.j2/2., J4, oa.j2/2., oj.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( J3, J1, J4, oa.j2/2., oj.j2/2., ok.j2/2. );
                             double sixj3 = H.modelspace->GetSixJ( oc.j2/2., oe.j2/2., J4, oi.j2/2., ob.j2/2., J2 );
                             double sixj4 = H.modelspace->GetSixJ( J3, J4, J1, oc.j2/2., od.j2/2., oe.j2/2. );
                             F32 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabij * vcdak * viebc * vjkde / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F32;
}

// Diagram F33 (as numbered by ADG)
// mscheme expression: F33 = -1/4 sum_abcdijkl (v_abik v_cdjl v_ijcd v_klab) / (eps^ik_ab eps^ijkl_cdab eps^kl_ab)
// minus sign error. Fixed.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F33( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F33 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   if ( ol.j2 != oi.j2) continue;
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_ijklcdab = GetDenom(H,{i,j,k,l},{c,d,a,b});
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double denom = e_ikab * e_ijklcdab * e_klab;
                   int phase_exponent = (oj.j2+ok.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2},{ok.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oj.j2,ol.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oj.j2,ol.j2},{oi.j2,oj.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     double vklab = H.TwoBody.GetTBME_J(J0,J0,k,l,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdjl = H.TwoBody.GetTBME_J(J1,J1,c,d,j,l);
                       double vijcd = H.TwoBody.GetTBME_J(J1,J1,i,j,c,d);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F33 += 1./4 * phase * (2*J0+1) * (2*J1+1) / (oi.j2+1) * vabik * vcdjl * vijcd * vklab / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F33;
}



// Diagram F34 (as numbered by ADG)
// mscheme expression: F34 = 1/16 sum_abcdijkl (v_abij v_cdkl v_ijcd v_klab) / (eps^ij_ab eps^ijkl_cdab eps^kl_ab)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F34( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F34 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijklcdab = GetDenom(H,{i,j,k,l},{c,d,a,b});
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double denom = e_ijab * e_ijklcdab * e_klab;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2},{ok.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2},{ok.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vcdkl = H.TwoBody.GetTBME_J(J0,J0,c,d,k,l);
                     double vijcd = H.TwoBody.GetTBME_J(J0,J0,i,j,c,d);
                     double vklab = H.TwoBody.GetTBME_J(J0,J0,k,l,a,b);
                     F34 += 1./16 * (2*J0+1) * vabij * vcdkl * vijcd * vklab / denom;
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F34;
}



// Diagram F35 (as numbered by ADG)
// mscheme expression: F35 = -1/4 sum_abcdijkl (v_abkl v_cdij v_ijac v_klbd) / (eps^kl_ab eps^ijkl_acbd eps^kl_bd)
// minus sign error. Fixed.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F35( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F35 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           if ( od.j2 != oa.j2) continue;
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double e_ijklacbd = GetDenom(H,{i,j,k,l},{a,c,b,d});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_klab * e_ijklacbd * e_klbd;
                   int phase_exponent = (ob.j2+oc.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabkl = H.TwoBody.GetTBME_J(J0,J0,a,b,k,l);
                     double vklbd = H.TwoBody.GetTBME_J(J0,J0,k,l,b,d);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdij = H.TwoBody.GetTBME_J(J1,J1,c,d,i,j);
                       double vijac = H.TwoBody.GetTBME_J(J1,J1,i,j,a,c);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F35 += 1./4 * phase * (2*J0+1) * (2*J1+1) / (oa.j2+1) * vabkl * vcdij * vijac * vklbd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F35;
}



// Diagram F36 (as numbered by ADG)
// mscheme expression: F36 = sum_abcdijkl (v_abik v_cdjl v_ijac v_klbd) / (eps^ik_ab eps^ijkl_acbd eps^kl_bd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F36( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F36 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_ijklacbd = GetDenom(H,{i,j,k,l},{a,c,b,d});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ikab * e_ijklacbd * e_klbd;
                   int phase_exponent = (ob.j2+oc.j2+oj.j2+ok.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oj.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oj.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   // J4 does not label any two-body matrix element; its range comes purely from
                   // the (single-particle) triangle conditions of the four 6j symbols it appears in.
                   int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{oj.j2,oc.j2},{od.j2,ol.j2}} ) /2;
                   int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{oj.j2,oc.j2},{od.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdjl = H.TwoBody.GetTBME_J(J1,J1,c,d,j,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vijac = H.TwoBody.GetTBME_J(J2,J2,i,j,a,c);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oj.j2/2., oc.j2/2., J4, oa.j2/2., oi.j2/2., J2 );
                             double sixj3 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, ok.j2/2., ob.j2/2., J3 );
                             double sixj4 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, oj.j2/2., oc.j2/2., J1 );
                             F36 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vcdjl * vijac * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F36;
}



// Diagram F37 (as numbered by ADG)
// mscheme expression: F37 = -1/4 sum_abcdijkl (v_abij v_cdkl v_ijac v_klbd) / (eps^ij_ab eps^ijkl_acbd eps^kl_bd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F37( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F37 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         if ( oc.j2 != ob.j2) continue;
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijklacbd = GetDenom(H,{i,j,k,l},{a,c,b,d});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ijab * e_ijklacbd * e_klbd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijac = H.TwoBody.GetTBME_J(J0,J0,i,j,a,c);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdkl = H.TwoBody.GetTBME_J(J1,J1,c,d,k,l);
                       double vklbd = H.TwoBody.GetTBME_J(J1,J1,k,l,b,d);
                       F37 += -1./4 * (2*J0+1) * (2*J1+1) / (ob.j2+1) * vabij * vcdkl * vijac * vklbd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F37;
}



// Diagram F38 (as numbered by ADG)
// mscheme expression: F38 = 1/16 sum_abcdijkl (v_abkl v_cdij v_ijab v_klcd) / (eps^kl_ab eps^ijkl_abcd eps^kl_cd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F38( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F38 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double e_ijklabcd = GetDenom(H,{i,j,k,l},{a,b,c,d});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_klab * e_ijklabcd * e_klcd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{ok.j2,ol.j2},{oc.j2,od.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{ok.j2,ol.j2},{oc.j2,od.j2},{oi.j2,oj.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabkl = H.TwoBody.GetTBME_J(J0,J0,a,b,k,l);
                     double vcdij = H.TwoBody.GetTBME_J(J0,J0,c,d,i,j);
                     double vijab = H.TwoBody.GetTBME_J(J0,J0,i,j,a,b);
                     double vklcd = H.TwoBody.GetTBME_J(J0,J0,k,l,c,d);
                     F38 += 1./16 * (2*J0+1) * vabkl * vcdij * vijab * vklcd / denom;
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F38;
}



// Diagram F39 (as numbered by ADG)
// mscheme expression: F39 = -1/4 sum_abcdijkl (v_abik v_cdjl v_ijab v_klcd) / (eps^ik_ab eps^ijkl_abcd eps^kl_cd)
// agrees.
// Written by AI (Claude sonnet 5) based on expressions and example from Ragnar. Tested and corrected by Ragnar.
double GetMP4_F39( const Operator& H)
{
   double t_start = omp_get_wtime();

   double F39 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 if ( ok.j2 != oj.j2) continue;
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_ijklabcd = GetDenom(H,{i,j,k,l},{a,b,c,d});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ikab * e_ijklabcd * e_klcd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oj.j2,ol.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oj.j2,ol.j2},{ok.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     double vijab = H.TwoBody.GetTBME_J(J0,J0,i,j,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdjl = H.TwoBody.GetTBME_J(J1,J1,c,d,j,l);
                       double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
                       F39 += -1./4 * (2*J0+1) * (2*J1+1) / (oj.j2+1) * vabik * vcdjl * vijab * vklcd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F39;
}



}
