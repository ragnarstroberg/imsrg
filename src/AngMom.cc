
#include "AngMom.hh"
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <iostream>


namespace AngMom
{
 
// double fct(double x)
// {
//    if (x<0) return 0;
//    double f=1;
//    for (double g=2;g<=x;g=g+1) {f *= g;}
//    return f;
// }

 double fct(double x)
 {
   return gsl_sf_fact(int(x));
 }

 int phase(int x)
 {
    return x%2==0 ? 1 : -1;
 }

 
 double Tri(double j1, double j2, double j3)
 {
    return fct(j1+j2-j3) * fct(j1-j2+j3) * fct(-j1+j2+j3)/fct(j1+j2+j3+1);
 }

 double CG(double ja, double ma, double jb, double mb, double J, double M)
 {
    return phase(ja-jb+M) * sqrt(2*J+1) * gsl_sf_coupling_3j(int(2*ja),int(2*jb),int(2*J),int(2*ma),int(2*mb),int(2*-M));
 }

 double SixJ(double j1, double j2, double j3, double J1, double J2,double J3)
 {
   return gsl_sf_coupling_6j(int(2*j1),int(2*j2),int(2*j3),int(2*J1),int(2*J2),int(2*J3));
 }

 double NineJ(double j1,double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
 {
   return gsl_sf_coupling_9j(int(2*j1),int(2*j2),int(2*J12),int(2*j3),int(2*j4),int(2*J34),int(2*J13),int(2*J24),int(2*J));
 }

 double NormNineJ(double j1,double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
 {
   return sqrt( (2*J12+1)*(2*J34+1)*(2*J13+1)*(2*J24+1) ) * NineJ(j1,j2,J12,j3,j4,J34,J13,J24,J);
 }


/* 
 //#####################################################################
 //  Wigner 6-J symbol { j1 j2 j3 }  See definition and Racah formula at
 //                    { J1 J2 J3 }  http://mathworld.wolfram.com/Wigner6j-Symbol.html
 double SixJ(double j1, double j2, double j3, double J1, double J2,double J3)
 {
   
   double triads[4][3] = {{j1,j2,j3},{j1,J2,J3},{J1,j2,J3},{J1,J2,j3}};
   
   for (int i=0;i<4;i++)
   {
      if ( (triads[i][0]+triads[i][1]<triads[i][2])
        || (triads[i][0]-triads[i][1]>triads[i][2])
        || ((int)(2*(triads[i][0]+triads[i][1]+triads[i][2]))%2>0) )
          return 0;
   } 
   double sixj = 0;
   for (double t=j1+j2+j3; t<=(j1+j2+j3+J1+J2+J3); t+=1)
   {
      double ff = fct(t-j1-j2-j3)*fct(t-j1-J2-J3)*
                 fct(t-J1-j2-J3)*fct(t-J1-J2-j3)*
                 fct(j1+j2+J1+J2-t)*fct(j2+j3+J2+J3-t)*
                 fct(j3+j1+J3+J1-t) ;
      if (ff>0)
      {
         sixj += pow(-1,t) * fct(t+1)/ff;
      }
   
   }
   for (int i=0;i<4;i++)
   {
      sixj *= sqrt(Tri( triads[i][0],triads[i][1],triads[i][2]));
   }
   return sixj;
 }
*/

/*
double NineJ(double j1,double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
{
   float ninej = 0;
   int ph = 1-2*(abs(int(J-j1))%2);
   for (float g = fabs(J-j1); g<=J+j1; g+=1)
   {
      ninej +=  ph * (2*g+1)
                * SixJ(j1,j2,J12,J34,J,g)
                * SixJ(j3,j4,J34,j2,g,J24)
                * SixJ(J13,J24,J,g,j1,j3);
   }
   return ninej;
}
*/




//  Talmi-Moshinsky Bracket, using the algorthm of Buck et al. Nuc. Phys. A 600 (1996) 387-402
//  Their phase convention differs from Moshinsky's by a factor (-1)**(L+l+lam).
//  I correct for this in order to stick to Moshinsky's phase convention.
//
double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam)
{
   // check energy conservation
   int f1 = 2*n1 + l1;
   int f2 = 2*n2 + l2;
   int F  = 2*N  + L;
   int f  = 2*n  + l;
   if (f1+f2 != f+F) return 0;

   // cos(2 beta) = (m2w2 - m1w1)/(m2w2 + m1w1)
   // for m1=m2 and w1=w2, take beta = pi/4
   double B = M_PI_4;
   double cB = cos(B);
   double sB = sin(B);

   double mosh1 = phase((l1+l2+L+l)/2) / pow(2.,(l1+l2+L+l)/4.0);
   mosh1 *= sqrt( gsl_sf_fact(n1) * gsl_sf_fact(n2) * gsl_sf_fact(N) * gsl_sf_fact(n) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n1+l1)+1) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n2+l2)+1) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(N+ L)+1)  );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n+ l)+1)  );

   double mosh2 = 0;
   for (int la=0;la<=min(f1,F);++la)
   {
    for (int lb=(la+l1)%2; lb<=l1+la; lb+=2)
    {
     double cg_ab =  CG(la,0,lb,0,l1,0);
     if (cg_ab == 0) continue;
     for (int lc=(la+L)%2; lc<=la+L; lc+=2)
     {
      double cg_ac = CG(la,0,lc,0,L,0);
      if (cg_ac == 0) continue;
      for (int ld=(lb+l)%2; ld<=min(lb+l,l2+lc); ld+=2)
      {
       double cg_bd = CG(lb,0,ld,0,l,0);
       if (cg_bd == 0) continue;
       double cg_cd = CG(lc,0,ld,0,l2,0);
       if (cg_cd == 0) continue;
       for (int a=0;a<=int((f1-la-lb)/2);++a)
       {
          int b = int((f1-2*a-la-lb)/2);
          int c = int((F-2*a-la-lc)/2);
          int d = int((f2-F+2*a+la-ld)/2);
          if (b<0 or c<0 or d<0) continue;
          double mosh3 = phase(la+lb+lc) * pow(2.,(la+lb+lc+ld)/2.);
          mosh3 *= pow(sB,2*a+la+2*d+ld);
          mosh3 *= pow(cB,2*b+lb+2*c+lc);
          mosh3 *= (2*la+1)/gsl_sf_fact(a)/gsl_sf_doublefact(2.*(a+la)+1);
          mosh3 *= (2*lb+1)/gsl_sf_fact(b)/gsl_sf_doublefact(2.*(b+lb)+1);
          mosh3 *= (2*lc+1)/gsl_sf_fact(c)/gsl_sf_doublefact(2.*(c+lc)+1);
          mosh3 *= (2*ld+1)/gsl_sf_fact(d)/gsl_sf_doublefact(2.*(d+ld)+1);
          mosh3 *= NineJ(la,lb,l1,lc,ld,l2,L,l,lam);
          mosh3 *= cg_ab * cg_ac * cg_bd * cg_cd;
          mosh2 += mosh3;

       }
      }
     }
    }
   }
   int extraphase = phase(l+L+lam);
   cout << "mosh1 = " << mosh1 << "   mosh2 = " << mosh2 << endl;
   return extraphase * mosh1 * mosh2;

}



} //end namespace AngMom

