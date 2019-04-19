
#include "AngMom.hh"
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <iostream>



namespace AngMom
{

 int phase(int x)
 {
    return x%2==0 ? 1 : -1;
 }

 
 double Tri(double j1, double j2, double j3)
 {
    return gsl_sf_fact(j1+j2-j3) * gsl_sf_fact(j1-j2+j3) * gsl_sf_fact(-j1+j2+j3)/gsl_sf_fact(j1+j2+j3+1);
 }
 bool Triangle(double j1, double j2, double j3)
 {
   if (round(2*(j1+j2))<round(2*j3)) return false;
   if (std::abs(round(2*(j1-j2)))>round(2*j3)) return false;
   return true;
 }

 double ThreeJ(double j1, double j2, double j3, double m1, double m2, double m3)
 {
   return gsl_sf_coupling_3j(int(2*j1), int(2*j2), int(2*j3), int(2*m1), int(2*m2), int(2*m3));
 }

 double CG(double ja, double ma, double jb, double mb, double J, double M)
 {
    return phase(ja-jb+M) * sqrt(2*J+1) * ThreeJ(ja,jb,J,ma,mb,-M);
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





//  Talmi-Moshinsky Bracket, using the algorthm of Buck et al. Nuc. Phys. A 600 (1996) 387-402
//  Their phase convention differs from Moshinsky's by a factor (-1)**(L+l+lam).
//  I correct for this in order to stick to Moshinsky's phase convention.
//
double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam, double B)
{
   // check energy conservation
   int f1 = 2*n1 + l1;
   int f2 = 2*n2 + l2;
   int F  = 2*N  + L;
   int f  = 2*n  + l;
   if (f1+f2 != f+F) return 0;

   // cos(2 beta) = (m2w2 - m1w1)/(m2w2 + m1w1)
   // for m1=m2 and w1=w2, take beta = pi/4
//   double B = M_PI_4;
//   double B = M_PI_4;
   double cB = cos(B);
   double sB = sin(B);

   double mosh1 = phase((l1+l2+L+l)/2) / pow(2.,(l1+l2+L+l)/4.0);
   mosh1 *= sqrt( gsl_sf_fact(n1) * gsl_sf_fact(n2) * gsl_sf_fact(N) * gsl_sf_fact(n) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n1+l1)+1) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n2+l2)+1) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(N+ L)+1)  );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n+ l)+1)  );

   double mosh2 = 0;
   for (int la=0;la<=std::min(f1,F);++la)
   {
//    for (int lb=(la+l1)%2; lb<=l1+la; lb+=2)
    for (int lb=(la+l1)%2; lb<=std::min(l1+la,f1-la); lb+=2)
    {
     double cg_ab =  CG(la,0,lb,0,l1,0);
     if (cg_ab == 0) continue;
     int amax = int((f1-la-lb)/2);
     if (amax<0) continue;
//     int bmax = int((f1-2*0-la-lb)/2); //TODO: This requirement can probably be incorporated in the limits for lb...
//     if (bmax<0) continue;
//     for (int lc=(la+L)%2; lc<=la+L; lc+=2)
     for (int lc=(la+L)%2; lc<=std::min(la+L,F-la); lc+=2)
     {
      double cg_ac = CG(la,0,lc,0,L,0);
      if (cg_ac == 0) continue;
//      int cmax = int((F-2*0-la-lc)/2);
//      if (cmax<0) continue;
      int ldmax = std::min( std::min(lb+l,l2+lc), f2-F+2*amax+la);
//      for (int ld=(lb+l)%2; ld<=std::min(lb+l,l2+lc); ld+=2)
      for (int ld=(lb+l)%2; ld<=ldmax; ld+=2)
      {
//       int dmax = int((f2-F+2*amax+la-ld)/2);
//       if (dmax<0) continue;
       double cg_bd = CG(lb,0,ld,0,l,0);
       if (cg_bd == 0) continue;
       double cg_cd = CG(lc,0,ld,0,l2,0);
       if (cg_cd == 0) continue;
       double ninej = NineJ(la,lb,l1,lc,ld,l2,L,l,lam);
       double mosh3_pre = phase(la+lb+lc) * pow(2.,(la+lb+lc+ld)/2.) *  ninej
                          * cg_ab * cg_ac * cg_bd * cg_cd
                          * (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1);
       double mosh3 = 0;
//       for (int a=0;a<=int((f1-la-lb)/2);++a)
       for (int a=0;a<=amax;++a)
       {
          int b = int((f1-2*a-la-lb)/2);
          int c = int((F-2*a-la-lc)/2);
          int d = int((f2-F+2*a+la-ld)/2);
          if (b<0 or c<0 or d<0) continue;
//          double mosh3 = phase(la+lb+lc) * pow(2.,(la+lb+lc+ld)/2.);
//          mosh3 *= pow(sB,2*a+la+2*d+ld);
//          mosh3 *= pow(cB,2*b+lb+2*c+lc);
//          mosh3 *= NineJ(la,lb,l1,lc,ld,l2,L,l,lam);
//          mosh3 *= ninej;
//          mosh3 *= cg_ab * cg_ac * cg_bd * cg_cd;
//          mosh3 *= (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1);
//          mosh3 /= ( gsl_sf_fact(a) * gsl_sf_doublefact(2.*(a+la)+1)
           mosh3 += pow(sB,2*a+la+2*d+ld) * pow(cB,2*b+lb+2*c+lc)
                  / ( gsl_sf_fact(a) * gsl_sf_doublefact(2.*(a+la)+1)
                    * gsl_sf_fact(b) * gsl_sf_doublefact(2.*(b+lb)+1)
                    * gsl_sf_fact(c) * gsl_sf_doublefact(2.*(c+lc)+1)
                    * gsl_sf_fact(d) * gsl_sf_doublefact(2.*(d+ld)+1) );
//          mosh2 += mosh3;

       }
       mosh2 += mosh3 * mosh3_pre;
      }
     }
    }
   }
   int extraphase = phase(l+L+lam);
//   cout << "mosh1 = " << mosh1 << "   mosh2 = " << mosh2 << endl;
   return extraphase * mosh1 * mosh2;

}




// Copied directly from Moshinsky & Brody 'Tables of Transformation Brackets for Nuclear Shell-Model Calculations' 1967
// Still need to check these against the tabulated values

double TalmiB(int n, int l, int nn, int ll, int p)
{
  double B = 0;
  int alpha = std::max(0, p-(l+ll)/2-nn);
  int beta  = std::min(n, p-(l+ll)/2);
  for (int k=alpha; k<=beta; k++)
  {
    B  +=   gsl_sf_fact(l+k)*gsl_sf_fact(p-(l-ll)/2-k) /
            ( gsl_sf_fact(k)*gsl_sf_fact(2*l+2*k+1)*gsl_sf_fact(n-k)*gsl_sf_fact(2*p-l+ll-2*k+1)*gsl_sf_fact(nn-p+(l+ll)/2+k)*gsl_sf_fact(p-(l+ll)/2-k)  ); 
  }
  B *= phase(p-(l+ll)/2) * gsl_sf_fact(2*p+1) / ( pow(2,n+nn)*gsl_sf_fact(p) )
         * sqrt(gsl_sf_fact(n)*gsl_sf_fact(nn)*gsl_sf_fact(2*n+2*l+1)*gsl_sf_fact(2*nn+2*ll+1) / (gsl_sf_fact(n+l)*gsl_sf_fact(nn+ll))  );

  return B;
}





//  double Tcoeff( LabKet& labket, int Jab, int twoJ, jacobi1_state& jac1, jacobi2_state& jac2, int twoJ12, int Ncm, int Lcm)
  double Tcoeff( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm)
  {
//    int na = labket.p.n;
//    int la = labket.p.l;
//    double ja = 0.5*labket.p.j2;
    double ja = 0.5*j2a;
    double sa = 0.5;
//    int nb = labket.q.n;
//    int lb = labket.q.l;
//    double jb = 0.5*labket.q.j2;
    double jb = 0.5*j2b;
    double sb = 0.5;
//    int nc = labket.r.n;
//    int lc = labket.r.l;
//    double jc = 0.5*labket.r.j2;
    double jc = 0.5*j2c;
    double sc = 0.5;
//
//    int N1 = jac1.n;
//    int L1 = jac1.l;
//    int S1 = jac1.s;
//    int J1 = jac1.j;
//    int N2 = jac2.n;
//    int L2 = jac2.l;
    double S2 = 0.5; 
//    int twoJ2 = jac2.j2;

    // Limits
    int Lab_min = std::max( std::abs(la-lb), std::abs(Jab-S1) );
    int Lab_max = std::min( la+lb , Jab+S1 );


    double tcoeff=0.;
    for (int Lab=Lab_min; Lab<=Lab_max; Lab++)
    {
      double ninejab = (Lab+1) * NineJ(la, lb, Lab, sa, sb, S1, ja, jb, Jab );
      int L12_min = std::max(twoJ12-S1-1, std::abs(L1-L2) );
      int L12_max = std::min((twoJ12+S1+1), L1+L2 );
      for (int L12=L12_min; L12<=L12_max; L12+=2)
      {

        int L_min = std::max( std::abs(Lcm-L12), std::abs(Lab-lc) );
        int L_max = std::min( Lcm+L12, Lab+lc );
        for (int L=L_min; L<=L_max; L++)
        {
          double t12=0;
          int twoS12_min = std::max( std::abs( 2*S1-1 ), std::max( std::abs(twoJ-2*L ), std::abs(2*L12-twoJ12) ) );
          int twoS12_max = std::min(2*S1 + 1 , std::min( twoJ+2*L, 2*L12+twoJ12) );
          for (int twoS12=twoS12_min; twoS12<=twoS12_max; twoS12+=2)  // use twice S12, to make it an integer
          {
            t12 += (twoS12+1) * NineJ(L1,L2,L12,S1,S2,0.5*twoS12,J1,0.5*twoJ2,0.5*twoJ12)
                              * NineJ(Lab,lc,L12,S1,sc,0.5*twoS12,Jab,jc,0.5*twoJ)
                              * SixJ(Lcm, L12,L,0.5*twoS12,0.5*twoJ,0.5*twoJ12) * phase( lc + Lab + L + L1 + (twoJ + twoS12)/2 );
          }

          double tmosh=0;
          int Lambda_min = std::max( std::abs(Lcm-L2), std::abs(L1-L));
          int Lambda_max = std::min( Lcm+L2, L1+L );
          for (int Lambda=Lambda_min; Lambda<=Lambda_max; Lambda++)
          {
            double sixjLambda = SixJ(Lcm,L2,Lambda,L1,L,L12);
            int curlyL_min = std::max( std::abs(lc-Lambda), std::abs(L1-Lab) );
            int curlyL_max = std::min( lc+Lambda, L1+Lab );
            for (int curlyL=curlyL_min; curlyL<=curlyL_max; curlyL+=2)
            {
              int curlyN = ((Ncm+nc) + (Lcm +lc-curlyL)/2);
              // increment moshinsky loop bit
              tmosh += sixjLambda * SixJ(lc,curlyL,Lambda,L1,L,Lab) * phase( Lambda)
                       * Moshinsky( curlyN,curlyL,  N1,L1,      nb,lb,     na,la,     Lab,  MOSH_BETA_1)
                       * Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  MOSH_BETA_2);

            }
          }
          // combine 12-loop and moshinsky loop bits;
          tcoeff += (2*L12+1) * (2*L+1)  * ninejab * t12 * tmosh;
        }
      }
    }
 
  // multiply by all the j-hat symbols.   
    double jhats = sqrt( (2*ja+1)*(2*jb+1)*(2*jc+1)*(2*Jab+1)*(twoJ+1)*(2*J1+1)*(twoJ2+1)*(2*S1+1) );
    return tcoeff * jhats;
  }








} //end namespace AngMom

