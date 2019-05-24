
#include "AngMom.hh"
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <iostream>


/// Namespace for angular momentum coupling functions
namespace AngMom
{

 /// returns \f$ (-1)^x \f$
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

 /// Wigner 3j symbol
 double ThreeJ(double j1, double j2, double j3, double m1, double m2, double m3)
 {
   return gsl_sf_coupling_3j(int(2*j1), int(2*j2), int(2*j3), int(2*m1), int(2*m2), int(2*m3));
 }

 /// Clebsch-Gordan coefficient
 double CG(double ja, double ma, double jb, double mb, double J, double M)
 {
    return phase(ja-jb+M) * sqrt(2*J+1) * ThreeJ(ja,jb,J,ma,mb,-M);
 }

 /// Wigner 6j symbol
 double SixJ(double j1, double j2, double j3, double J1, double J2,double J3)
 {
//   return gsl_sf_coupling_6j(int(2*j1),int(2*j2),int(2*j3),int(2*J1),int(2*J2),int(2*J3));
   return SixJ_int(int(2*j1),int(2*j2),int(2*j3),int(2*J1),int(2*J2),int(2*J3));
 }

 double SixJ_int(int twoj1, int twoj2, int twoj3, int twoJ1, int twoJ2, int twoJ3)
 {
   return gsl_sf_coupling_6j(twoj1,twoj2,twoj3,twoJ1,twoJ2,twoJ3);
 }

 /// Wigner 9j symbol
 double NineJ(double j1,double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
 {
   return gsl_sf_coupling_9j(int(2*j1),int(2*j2),int(2*J12),int(2*j3),int(2*j4),int(2*J34),int(2*J13),int(2*J24),int(2*J));
 }


 double NineJ_int(int j1,int j2, int J12, int j3, int j4, int J34, int J13, int J24, int J)
 {
   return gsl_sf_coupling_9j(j1,j2,J12,j3,j4,J34,J13,J24,J);
 }

 /// Normalized Wigner 9j symbol
 double NormNineJ(double j1,double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
 {
   return sqrt( (2*J12+1)*(2*J34+1)*(2*J13+1)*(2*J24+1) ) * NineJ(j1,j2,J12,j3,j4,J34,J13,J24,J);
 }

// double General3nJ_FirstKind( std::vector<int>& jvals )
 double General3nJ_FirstKind( std::vector<int>& a, std::vector<int>& b, std::vector<int>& c )
 {
   int n = a.size();
   int S = 0;
   int xmin=0;
   int xmax=2147483647; // maximum 32-bit signed integer
   for ( int i=0; i<n; i++)
   {
     S += a[i] + b[i] + c[i];
     xmin = std::max( xmin, std::abs(a[i]-c[i]) );
     xmax = std::min( xmax, a[i]+c[i] );
   }

   std::cout << "n = " << n << "  xmin,xmax = " << xmin << " " << xmax << std::endl;
   double sum_3nj = 0;
   for (int twox=xmin; twox<=xmax; twox+=2)
   {
     double prodx = (twox+1) * AngMom::phase((S+(n-1)*twox)/2);
     std::cout << "before i, prodx = " << prodx << std::endl;
     for (int i=0; i<n-1; i++)
     {
       prodx *= SixJ_int( a[i], a[i+1], b[i], c[i+1], c[i], twox );
       std::cout << "i = " << i << "   SixJ( " << a[i] << " " << a[i+1] << " " << b[i] << " " << c[i+1] << " " << c[i] << " " << twox << " ) = "
                 << SixJ_int( a[i], a[i+1], b[i], c[i+1], c[i], twox ) << "   prodx = " << prodx << std::endl;
     }
     prodx *= SixJ_int( a[n-1], c[0], b[n-1], a[0], c[n-1], twox );
       std::cout << "i = " << n-1 << "   SixJ( " << a[n-1] << " " << c[0] << " " << b[n-1] << " " << a[0] << " " << c[n-1] << " " << twox << " ) = "
                 << SixJ_int( a[n-1], c[0], b[n-1], a[0], c[n-1], twox ) << "  prodx = " << prodx << std::endl;
     sum_3nj += prodx;
     std::cout << "sum_3nj = " << sum_3nj << std::endl;
   }
   return sum_3nj;
 }


 double FancyNineJ( int twoj1, int twoj2, int twoJ12, int twoj3, int twoj4, int twoJ34, int twoJ13, int twoJ24, int twoJ )
 {
   std::vector<int> a = { twoj1, twoj2, twoJ12};
   std::vector<int> b = { twoj3, twoj4, twoJ34};
   std::vector<int> c = { twoJ13, twoJ24, twoJ};
   double ninej = General3nJ_FirstKind( a,b,c );
//   double ninej=1.0;
   return ninej;
 }

 double FancyTwelveJ(int a1, int a2, int a3, int a4, int b12, int b23, int b34, int b41, int c1, int c2, int c3, int c4)
 {
   std::vector<int> a = { a1,a2,a3,a4};
   std::vector<int> b = { b12,b23,b34,b41};
   std::vector<int> c = { c1,c2,c3,c4};
   double twelvej = General3nJ_FirstKind( a,b,c );
//   double ninej=1.0;
   return twelvej;
 }

 /// 12j symbol of the first kind, Varshalovich pg 362 eq (6).
 /// Accepts twice the j value for each entry
 // TODO: It might help to take advantage of the simplification if one of the arguments is zero
 double TwelveJ_1_ints(int a1, int a2, int a3, int a4, int b12, int b23, int b34, int b41, int c1, int c2, int c3, int c4)
 {
   int S = (a1+a2+a3+a4+b12+b23+b34+b41+c1+c2+c3+c4);
   int twoXmin  = std::max( std::abs(a1-c1), std::max( std::abs(a2-c2), std::max( std::abs(a3-c3), std::abs(a4-c4)) ) );
   int twoXmax  = std::min( a1+c1, std::min( a2+c2, std::min( a3+c3, a4+c4 ) ) );
//   twoXmin=0; twoXmax = 15;
   double twelvej = 0.;
   std::cout << "[  " << a1 << " " << a2 << " " << a3 << " " << a3 << " "
             << b12 << " " << b23 << " " << b34 << " " << b41 << " "
             << c1 << " " << c2 << " " << c3 << " " << c4 << " ]" << std::endl;
//   for (int twoX=twoXmin; twoX<=twoXmax; twoX+=2)
   for (int twoX=twoXmin; twoX<=twoXmax; twoX+=1)
   {
     double sixj1 = gsl_sf_coupling_6j(a1,a2,b12,c2,c1,twoX);
     double sixj2 = gsl_sf_coupling_6j(a2,a3,b23,c3,c2,twoX);
     double sixj3 = gsl_sf_coupling_6j(a3,a4,b34,c4,c3,twoX);
//     double sixj4 = gsl_sf_coupling_6j(a4,a1,b41,c1,c4,twoX);
     double sixj4 = gsl_sf_coupling_6j(a4,c1,b41,a1,c4,twoX);
     twelvej += (twoX+1) * AngMom::phase( (S-twoX)/2 ) * sixj1 * sixj2 * sixj3 * sixj4;
     std::cout << "In TwelveJ. twoX = " << twoX << "  -> "
               << sixj1 << " " << sixj2 << " " << sixj3 << " " << sixj4 << "   " << twelvej << std::endl;
   }
   std::cout << "Done with the loop" << std::endl;
   return twelvej;
 }

 /// 12j symbol of the first kind. Accepts the integer or half-integr value of each entry (i.e. not twice the value)
 double TwelveJ_1(double a1, double a2, double a3, double a4, double b12, double b23, double b34, double b41, double c1, double c2, double c3, double c4)
 {
   return AngMom::TwelveJ_1_ints( a1*2, a2*2, a3*2, a4*2, b12*2, b23*2, b34*2, b41*2, c1*2, c2*2, c3*2, c4*2);
 }

///  Talmi-Moshinsky Bracket, using the algorthm of Buck et al. Nuc. Phys. A 600 (1996) 387-402
///  Their phase convention differs from Moshinsky's by a factor \f$(-1)^{L+l+\Lambda}\f$.
///  I correct for this in order to stick to Moshinsky's phase convention.
///  I also use the mass-ratio parameter d from Kamuntavicius et al, NPA 695,191 (2001)
///  rather than the beta used by Buck and Merchant. They are related by
///  \f$\sin^2(\beta) = \frac{d}{1+d}\$f,  \f$\cos^2(\beta) = \frac{1}{1+d}\f$.
//
//double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam, double B)
double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam, double d)
{
   // check energy conservation
   int f1 = 2*n1 + l1;
   int f2 = 2*n2 + l2;
   int F  = 2*N  + L;
   int f  = 2*n  + l;
   if (f1+f2 != f+F) return 0;

   // cos(2 beta) = (m2*hw2 - m1*hw1)/(m2*hw2 + m1*hw1)
   // tan(beta) = sqrt( m1*hw1 / m2*hw2 )
   // for m1=m2 and w1=w2, take beta = pi/4, defined as MOSH_BETA_1
   // for m1=2*m2 and hw1=hw2, tan(beta) = sqrt(2), defined as MOSH_BETA_2 
   // for m2=2*m1 and hw1=hw2, tan(beta) = sqrt(0.5), defined as MOSH_BETA_half
   // the d parameter in Kamuntavicius et al NPA 695,191 (2001) is related to beta
   // as d/(1+d) = sin^2(beta) , 1/(1+d) = cos^2(beta)
   // so d corresponds to beta = atan( sqrt( d ) )
//   double cB = cos(B);
//   double sB = sin(B);
   double cB = sqrt(1/(1.0+d));
   double sB = sqrt(d/(1.0+d));

   double mosh1 = phase((l1+l2+L+l)/2) / pow(2.,(l1+l2+L+l)/4.0);
   mosh1 *= sqrt( gsl_sf_fact(n1) * gsl_sf_fact(n2) * gsl_sf_fact(N) * gsl_sf_fact(n) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n1+l1)+1) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n2+l2)+1) );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(N+ L)+1)  );
   mosh1 *= sqrt( gsl_sf_doublefact(2*(n+ l)+1)  );

   double mosh2 = 0;
   for (int la=0;la<=std::min(f1,F);++la)
   {
    for (int lb=(la+l1)%2; lb<=std::min(l1+la,f1-la); lb+=2)
    {
     double cg_ab =  CG(la,0,lb,0,l1,0);
     if (cg_ab == 0) continue;
     int amax = int((f1-la-lb)/2);
     if (amax<0) continue;
     for (int lc=(la+L)%2; lc<=std::min(la+L,F-la); lc+=2)
     {
      double cg_ac = CG(la,0,lc,0,L,0);
      if (cg_ac == 0) continue;
      int ldmax = std::min( std::min(lb+l,l2+lc), f2-F+2*amax+la);
      for (int ld=(lb+l)%2; ld<=ldmax; ld+=2)
      {
       double cg_bd = CG(lb,0,ld,0,l,0);
       if (cg_bd == 0) continue;
       double cg_cd = CG(lc,0,ld,0,l2,0);
       if (cg_cd == 0) continue;
       double ninej = NineJ(la,lb,l1,lc,ld,l2,L,l,lam);
       double mosh3_pre = phase(la+lb+lc) * pow(2.,(la+lb+lc+ld)/2.) *  ninej
                          * cg_ab * cg_ac * cg_bd * cg_cd
                          * (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1);
       double mosh3 = 0;
       for (int a=0;a<=amax;++a)
       {
          int b = int((f1-2*a-la-lb)/2);
          int c = int((F-2*a-la-lc)/2);
          int d = int((f2-F+2*a+la-ld)/2);
          if (b<0 or c<0 or d<0) continue;

           mosh3 += pow(sB,2*a+la+2*d+ld) * pow(cB,2*b+lb+2*c+lc)
                  / ( gsl_sf_fact(a) * gsl_sf_doublefact(2.*(a+la)+1)
                    * gsl_sf_fact(b) * gsl_sf_doublefact(2.*(b+lb)+1)
                    * gsl_sf_fact(c) * gsl_sf_doublefact(2.*(c+lc)+1)
                    * gsl_sf_fact(d) * gsl_sf_doublefact(2.*(d+ld)+1) );

       }
       mosh2 += mosh3 * mosh3_pre;
      }
     }
    }
   }
   int extraphase = phase(l+L+lam); // See note at the top about switching phase conventions
   return extraphase * mosh1 * mosh2;

}




/// Copied directly from Moshinsky & Brody 'Tables of Transformation Brackets for Nuclear Shell-Model Calculations' 1967
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





  double Tcoeff_bruteforce( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm)
  {
    double ja = 0.5*j2a;
    double jb = 0.5*j2b;
    double jc = 0.5*j2c;
    double sa = 0.5;
    double sb = 0.5;
    double sc = 0.5;
    double S2 = 0.5;


//    std::cout << std::endl << " " << __func__ << " ( " << na << " " << la << " " << j2a << "  " << nb << " " << lb << " " << j2b << "  " << nc << " " << lc << " " << j2c
//              << "   " << Jab << " " << twoJ << "  " << N1 << " " << L1 << " " << S1 << "   " << N2 << " " << L2 << "  " << twoJ2 << " " << twoJ12 << "  " << Ncm << " " << Lcm << " ) " << std::endl;

    double tcoeff = 0;
    int L12_min = std::abs(L1 - L2);
    int L12_max = L1+L2;
    int Lab_min = std::max( std::abs(la-lb), std::abs(Jab-S1) ) ;
    int Lab_max = std::min( la+lb ,  Jab+S1 ) ;

//    std::cout << " *** la,lb,Jab,S1 = " << la << " " << lb << " " << Jab << " " << S1 << "  Lab min/max = " << Lab_min << " " << Lab_max << std::endl;
//    std::cout << " *** L1,L2,J12,S1 = " << L1 << " " << L2 << " " << twoJ12 << " " << S1 << ",      L12 min/max = " << L12_min << " " << L12_max << std::endl;
    for ( int twoS12=1; twoS12<=(2*S1+1); twoS12+=2)
    {
      double S12 = 0.5 * twoS12;
      for ( int L12=L12_min; L12<=L12_max; L12+=1)
      {
        if (  (std::abs(twoS12-2*L12) > twoJ12)  or (twoS12+2*L12<twoJ12) ) continue; // check triangle condition for (S12 L12, J12)
        for ( int Lab=Lab_min; Lab<=Lab_max; Lab+=1 )
        {
          int L_min = std::max( std::abs(Lab-lc) ,  std::abs(twoJ-twoS12)/2 );
          int L_max = std::min( Lab+lc , (twoJ+twoS12)/2 );
          for (int L=L_min; L<=L_max; L+=1 )
          {
            int curlyL_min = std::abs( L1-Lab );
            int curlyL_max = L1+Lab;
            for (int curlyL=curlyL_min; curlyL<=curlyL_max; curlyL+=1)
            {
              int Lambda_min = std::max( std::abs(lc - curlyL), std::max( std::abs(L1-L), std::abs( L2-Lcm ) )  );
              int Lambda_max = std::min( lc+curlyL, std::min( L1+L, L2+Lcm) );
              for (int Lambda=Lambda_min; Lambda<=Lambda_max; Lambda+=1)
              {

//                // energy cons. 2*curlyN + curlyL + 2*N1+L1 = 2*nb+lb * 2*na+la from first moshinsky bracket
                if( (2*na+la+2*nb+lb - 2*N1-L1-curlyL)%2>0 ) continue;
                int curlyN = (2*na+la+2*nb+lb - 2*N1-L1-curlyL)/2;
                if (curlyN<0) continue;
//                // more energy conservation
                if ( (2*na+la+2*nb+lb+2*nc+lc) != (2*Ncm+Lcm + 2*N1+L1 + 2*N2+L2) ) continue;

//                double mosh1 = AngMom::Moshinsky(curlyN,curlyL, N1,L1, na,la, nb,lb, Lab, MOSH_BETA_1);
//                double mosh2 = AngMom::Moshinsky(Ncm,Lcm, N2,L2, curlyN,curlyL, nc,lc, Lambda, MOSH_BETA_2);
                double mosh1 = AngMom::Moshinsky(curlyN,curlyL, N1,L1, na,la, nb,lb, Lab, 1.0);
                double mosh2 = AngMom::Moshinsky(Ncm,Lcm, N2,L2, curlyN,curlyL, nc,lc, Lambda, 2.0);
                double ninej1 = AngMom::NineJ( la,lb,Lab, sa,sb,S1, ja,jb,Jab );
                double ninej2 = AngMom::NineJ( Lab,lc,L, S1,sc,S12, Jab,jc,0.5*twoJ) ;
                double ninej3 = AngMom::NineJ( L1,L2,L12, S1,S2,S12, J1,0.5*twoJ2,0.5*twoJ12 );
                double sixj1 = AngMom::SixJ( lc,curlyL,Lambda, L1,L,Lab );
                double sixj2 = AngMom::SixJ( Lcm,L2,Lambda, L1,L,L12 );
                double sixj3 = AngMom::SixJ( Lcm,L12,L, S12,0.5*twoJ,0.5*twoJ12 );
              // matching notation to Nogga et al PRC 73 064002 (2006)
              //  Nogga     Roth
              //  L12       Lab
              //  J12       Jab
              //  n12       N1
              //  l12       L1
              //  s12       S1
              //  j12       J1
              //  n3        N2
              //  l3        L2
              //  I3        J2
              //  S3        S12
              //  L3        L12
              //  curlyL    L
              //  curlyL12  curlyL
              //  N12       curlyN
              //  ncm       Ncm
              //  lcm       Lcm
              //  curlyJ    J
              //  J         J12
              double hats_all =  sqrt( (j2a+1)*(j2b+1)*(j2c+1)*(2*Jab+1)*(twoJ12+1)*(2*J1+1)*(twoJ2+1) * (2*S1+1) ) * (twoS12+1)*(2*Lab+1)*(2*L+1)*(2*L12+1)  * (2*Lambda+1);
              int phase_123 = AngMom::phase( lc+Lambda+Lab+L+L1+(twoS12+twoJ)/2 );

                int phasemosh1 = AngMom::phase( curlyL + L1 + Lab + 0*(curlyL+L1-la-lb)/2 );
                int phasemosh2 = AngMom::phase( Lcm + L2 + Lambda + 0*(Lcm+L2-curlyL-lc)/2);
              tcoeff += hats_all * phase_123 * ninej1 *  sixj1 *  ninej2 *  ninej3 *  sixj2 * sixj3 * mosh2 * mosh1 * phasemosh1 * phasemosh2;

              }
            }
          }
        }
      }
    }
//    std::cout << "*#*#*  " << tcoeff << std::endl;
    return tcoeff;
  }



/// Compute the T coefficient for transforming between a lab-frame 3-body state \f$|a,b,c, Jab, J\rangle\f$ where a=> \f$n_a,la,ja\f$   to a non-antisymmetrized Jacobi 3-body state
/// specified by \f$[ \{|N_1 L_1 S_1 J_1> \otimes |N_2 L_2 J_2> \}^{J_{12}} \otimes |N_{cm},L_{cm}\rangle ]^{J}\f$
/// The formula is given in the appendix of Nogga et al, Phys. Rev. C 73, 064002 (2006).
/// However, when writing this, I adopted the notation of  Roth et al. PRC 90, 024325 (2014).
/// Unfortunately, there appear to be mistakes in the latter paper, and they have been corrected here.
/// Specifically, in order to get agreement with the results of Nogga et al, the following corrections to the second line of equation (10) in Roth et al were needed
/// * The \f$\hat{J}\f$ should be \f$\hat{J}_{12}\f$
/// * The \f$\hat{\mathcal{L}}^2\f$ should be removed
/// * The Moshinsky bracket \f$\langle\langle N L, N_1 L_1 ;L_{ab}|n_b l_b, n_a l_a \rangle\rangle \f$ should have the right side flipped, i.e. \f$|n_a l_a, n_b l_b>>\f$
///
/// In Nogga et al, the Moshinsky brackets are referred to  Kamuntavicius et al. NPA 695 191 (2001).
/// This reference uses oscillator wave functions with a phase \f$(-1)^n\f$, and they state that they have a corresponding relative phase of \f$(-1)^{ (L+l-l_1-l_2)/2}\f$ compared to that of Buck and Merchant.
/// Buck and Merchant, in turn have a phase of  \f$(-1)^{L+l+\Lambda}\f$  relative to Moshinsky's phase, which is what the function AngMom::Moshinsky() returns.
/// So to get the Kamuntavicius convention, we must multiply by \f$(-1)^{L+l+\Lambda +(L+l-l_1-l_2)/2 }\f$
/// However, comparing with Petr Navratil's T coefficients, I don't need the additional phase, i.e. I just use the Buck/Merchant convention.
/// After all that mess, the quantity which is computed here is
/// \f{equation}{
/// \begin{aligned}
/// T &= \sum_{\mathcal{N},\mathcal{L},\Lambda}\sum_{L,L_{ab},L_{12},S_{12}} (-1)^{l_c+\Lambda+L_{ab}+L+L_1+S_{12}+J} \\
/// &\times \hat{j}_a \hat{j}_b \hat{j}_c \hat{J}_{ab} \hat{J}_{1} \hat{J}_{2} \hat{J}_{12} \hat{S}_{1}  \hat{S}_{12}^2 \hat{L}_{ab}^2 \hat{L}^2 \hat{L}_{12}^2 \hat{\Lambda}^2
/// \langle\langle \mathcal{N}\mathcal{L},N_1 L_1; L_{ab}| n_a l_a n_b l_b \rangle\rangle_1 \\
/// &\times \langle\langle N_{cm}L_{cm},N_2 L_2 ; \Lambda | \mathcal{N}\mathcal{L},n_c l_c \rangle\rangle_2
/// \left\{\begin{array}{lll} l_a    & l_b & L_{ab} \\ s_a & s_b & S_1    \\ j_a    & j_b & J_{ab} \end{array}\right\}
/// \left\{\begin{array}{lll} L_{ab} & l_c & L      \\ S_1 & s_c & S_{12} \\ J_{ab} & j_c & J \end{array}\right\}
/// \left\{\begin{array}{lll} L_1    & L_2 & L_{12} \\ S_1 & S_2 & S_{12} \\ J_1    & J_2 & J_{12} \end{array}\right\}  \\
/// &\times
/// \left\{\begin{array}{lll} l_c    & \mathcal{L} & \Lambda \\ L_1    & L & L_{ab}  \end{array}\right\}
/// \left\{\begin{array}{lll} L_{cm} & L_2         & \Lambda \\ L_1    & L & L_{12}  \end{array}\right\}
/// \left\{\begin{array}{lll} L_{cm} & L_{12}      &  L      \\ S_{12} & J & J_{12}  \end{array}\right\}
/// \end{aligned}
/// \f}
  double Tcoeff( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm)
  {
    double ja = 0.5*j2a;
    double sa = 0.5;
    double jb = 0.5*j2b;
    double sb = 0.5;
    double jc = 0.5*j2c;
    double sc = 0.5;
//
    double S2 = 0.5; 
//    std::cout << std::endl << " " << __func__ << " ( " << na << " " << la << " " << j2a << "  " << nb << " " << lb << " " << j2b << "  " << nc << " " << lc << " " << j2c
//              << "   " << Jab << " " << twoJ << "  " << N1 << " " << L1 << " " << S1 << "   " << N2 << " " << L2 << "  " << twoJ2 << " " << twoJ12 << "  " << Ncm << " " << Lcm << " ) " << std::endl;

    // Limits
    int Lab_min = std::max( std::abs(la-lb), std::abs(Jab-S1) );
    int Lab_max = std::min( la+lb , Jab+S1 );
//    std::cout << "   Lab min/max = " << Lab_min << " " << Lab_max << std::endl;


    double tcoeff=0.;
    for (int Lab=Lab_min; Lab<=Lab_max; Lab++)
    {
      double ninejab = (2*Lab+1) * NineJ(la, lb, Lab, sa, sb, S1, ja, jb, Jab ); 
//      std::cout << "   Lab , ninej = " << Lab << " " << ninejab << std::endl;
      if (std::abs(ninejab)<1e-8) continue; // <- This may not help very much
      int L12_min = std::max((twoJ12-2*S1-1)/2, std::abs(L1-L2) );
      int L12_max = std::min((twoJ12+2*S1+1)/2, L1+L2 );
      for (int L12=L12_min; L12<=L12_max; L12+=1)
      {

        int L_min = std::max( std::abs(Lcm-L12), std::abs(Lab-lc) );
        int L_max = std::min( Lcm+L12, Lab+lc );
//        std::cout << "   L12, L_min,L_max = " << L12 << " " << L_min << " " << L_max << std::endl;
        for (int L=L_min; L<=L_max; L++)
        {
          double t12=0;
          int twoS12_min = std::max( std::abs( 2*S1-1 ), std::max( std::abs(twoJ-2*L ), std::abs(2*L12-twoJ12) ) );
          int twoS12_max = std::min(2*S1 + 1 , std::min( twoJ+2*L, 2*L12+twoJ12) );
          for (int twoS12=twoS12_min; twoS12<=twoS12_max; twoS12+=2)  // use twice S12, to make it an integer
          {
            t12 += (twoS12+1) * NineJ(L1,L2,L12,S1,S2,0.5*twoS12,J1,0.5*twoJ2,0.5*twoJ12)
                              * NineJ(Lab,lc,L,S1,sc,0.5*twoS12,Jab,jc,0.5*twoJ)
                              * SixJ(Lcm, L12,L,0.5*twoS12,0.5*twoJ,0.5*twoJ12) * phase( lc + Lab + L + L1 + (twoJ + twoS12)/2 );
          }
          if (std::abs(t12)<1e-9) continue;
//          std::cout << "   Lab,L12,L = " << Lab << " " << L12 << " " << L << "   t12 = " << t12 << std::endl;

          double tmosh=0;
          int Lambda_min = std::max( std::abs(Lcm-L2), std::abs(L1-L));
          int Lambda_max = std::min( Lcm+L2, L1+L );
          for (int Lambda=Lambda_min; Lambda<=Lambda_max; Lambda++)
          {
            double sixjLambda = (2*Lambda+1) * SixJ(Lcm,L2,Lambda,L1,L,L12);
            if ( std::abs(sixjLambda)<1e-9) continue;
            
            int curlyL_min = std::max( std::abs(lc-Lambda), std::abs(L1-Lab) );
            int curlyL_max = std::min( lc+Lambda, L1+Lab );
            for (int curlyL=curlyL_min; curlyL<=curlyL_max; curlyL+=1)
            {
              int curlyN = na+nb-N1  +(la+lb-L1-curlyL)/2;
              if (curlyN<0) continue;
              if ( (la+lb+L1+curlyL)%2>0 ) continue;
              if ( 2*curlyN+curlyL +2*nc+lc != 2*Ncm+Lcm + 2*N2+L2) continue;
              double sixjCurly =  SixJ(lc,curlyL,Lambda,L1,L,Lab) ;  // I think the 2*curlyL+1 in eq (10) in the Roth paper isn't supposed to be there
              if (std::abs( sixjCurly)<1e-9) continue;
              // increment moshinsky loop bit
              int moshphase1 = AngMom::phase( curlyL+L1+Lab + 0*(curlyL+L1 - lb-la)/2 );  // converting phase conventions
              int moshphase2 = AngMom::phase( Lcm+L2+Lambda + 0*(Lcm+L2 - curlyL-lc)/2 ); // converting phase conventions (comment at the beginning of the function)
              tmosh += sixjLambda * sixjCurly * phase( Lambda) * moshphase1 * moshphase2
                       * Moshinsky( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab,  1.0)
                       * Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  2.0);
//                       * Moshinsky( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab,  MOSH_BETA_1)
//                       * Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  MOSH_BETA_2);

//              std::cout << "curlyN,curlyL = " << curlyN << " " << curlyL << "  "
//                        << sixjLambda << " " << sixjCurly << " " << moshphase1 << " " << moshphase2 << "   "
//                        << Moshinsky( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab,  1.0) << "   "
//                        << Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  2.0)
//                        << std::endl;
            }
          }
//          std::cout << "     tmosh = " << tmosh << std::endl;
          // combine 12-loop and moshinsky loop bits;
          tcoeff += (2*L12+1) * (2*L+1)  * ninejab * t12 * tmosh;
        }
      }
    }
 
  // multiply by all the j-hat symbols.   
  //  We also have (2S12+1) included in t12,  (2L12+1)(2L+1) included in tcoeff a few lines above,  (2Lab+1) included in ninejab, (2Lambda+1) included in sixjLambda, and (2curlyL+1) included sixjCurly
    double jhats = sqrt( (2*ja+1)*(2*jb+1)*(2*jc+1)*(2*Jab+1)*(twoJ12+1)*(2*J1+1)*(twoJ2+1)*(2*S1+1) );
//    std::cout << "  ->  " << tcoeff << " * " << jhats << "  =  " << tcoeff * jhats << std::endl;
    return tcoeff * jhats;
  }




  double Tcoeff_reorder( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm)
  {
    double ja = 0.5*j2a;
    double sa = 0.5;
    double jb = 0.5*j2b;
    double sb = 0.5;
    double jc = 0.5*j2c;
    double sc = 0.5;
//
    double S2 = 0.5; 
//    std::cout << std::endl << " " << __func__ << " ( " << na << " " << la << " " << j2a << "  " << nb << " " << lb << " " << j2b << "  " << nc << " " << lc << " " << j2c
//              << "   " << Jab << " " << twoJ << "  " << N1 << " " << L1 << " " << S1 << "   " << N2 << " " << L2 << "  " << twoJ2 << " " << twoJ12 << "  " << Ncm << " " << Lcm << " ) " << std::endl;

    // Limits
    int Lab_min = std::max( std::abs(la-lb), std::abs(Jab-S1) );
    int Lab_max = std::min( la+lb , Jab+S1 );
    int L12_min = std::abs(L1-L2);
    int L12_max = L1+L2;
    int twoS12_min = std::abs(2*S1-1);
    int twoS12_max = 2*S1+1;
    if ( Lab_max<Lab_min or L12_max<L12_min or twoS12_max<twoS12_min ) return 0.0;
//    std::cout << "   Lab min/max = " << Lab_min << " " << Lab_max << std::endl;

    std::vector<double> ninej_ab(Lab_max-Lab_min+1, 0.);
    for (int Lab=Lab_min; Lab<=Lab_max; Lab++)
    {
      ninej_ab[Lab-Lab_min] = AngMom::phase(Lab) * (2*Lab+1) * NineJ(la, lb, Lab, sa, sb, S1, ja, jb, Jab ); 
    }
    std::vector<double> ninej_12(2*(L12_max-L12_min+1),0);
    for (int L12=L12_min; L12<=L12_max; L12++)
    {
      for (int twoS12=twoS12_min; twoS12<=twoS12_max; twoS12+=2)
      {
        ninej_12[ 2*(L12-L12_min)+twoS12/2] = AngMom::phase((twoS12+twoJ)/2)*(2*L12+1)*(twoS12+1) * NineJ(L1,L2,L12,S1,S2,0.5*twoS12,J1,0.5*twoJ2,0.5*twoJ12);
      }
    }
//    std::cout << " done with precompute step " << std::endl;

    double tcoeff=0.;
    for (int Lab=Lab_min; Lab<=Lab_max; Lab++)
    {
//      double ninejab = (2*Lab+1) * NineJ(la, lb, Lab, sa, sb, S1, ja, jb, Jab ); 
//      std::cout << "   Lab , ninej = " << Lab << " " << ninejab << std::endl;
      double ninejab = ninej_ab[Lab-Lab_min];
      if (std::abs(ninejab)<1e-8) continue; 

      // Moshinsky's are expensive, so let's do them in the outer loop (maybe even more outer than this??)
      int curlyL_min =  std::abs(L1-Lab) ;
      int curlyL_max =  L1+Lab ;
      if ( (curlyL_min + L1 + la + lb)%2>0 ) curlyL_min++; // parity for the first Moshinsky bracket
      if ( (curlyL_min + lc + Lcm + L2)%2>0 ) continue; // parity for the second Moshinsky bracket
      for (int curlyL=curlyL_min; curlyL<=curlyL_max; curlyL+=2)
      {
        int curlyN = na+nb-N1  +(la+lb-L1-curlyL)/2;
        if (curlyN<0) continue;
//        if ( (la+lb+L1+curlyL)%2>0 ) continue; 
        if ( 2*curlyN+curlyL +2*nc+lc != 2*Ncm+Lcm + 2*N2+L2) continue; // energy conservation in second Moshinsky bracket
        int mosh1phase = AngMom::phase( curlyL+L1+Lab + 0*(curlyL+L1 - lb-la)/2 );
        double mosh1 = mosh1phase  * Moshinsky( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab,  1.0);
        if (std::abs(mosh1)<1e-9) continue;

        int Lambda_min = std::max( std::abs(Lcm-L2), std::abs(curlyL-lc) );
        int Lambda_max = std::min( Lcm+L2, curlyL+lc );
        for (int Lambda=Lambda_min; Lambda<=Lambda_max; Lambda++) 
        {
          int moshphase2 = AngMom::phase( Lcm+L2+Lambda + 0*(Lcm+L2 - curlyL-lc)/2 ); // converting phase conventions (comment at the beginning of the function)
          double mosh2 = moshphase2 * (2*Lambda+1) * Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  2.0);
          if (std::abs(mosh2)<1e-9) continue;

          // This inner loop is all just summing over sixj's and a ninej, and can probably be re-expressed with some fancy recoupling
          double sum_L = 0;
          int L_min = std::max( std::abs(Lambda-L1), std::abs(Lab-lc) );
          int L_max = std::min( Lambda+L1, Lab+lc );
          for (int L=L_min; L<=L_max; L++)
          {
            double sixj1 =  AngMom::phase(L+Lambda) * (2*L+1) * AngMom::SixJ( lc, curlyL, Lambda, L1,L,Lab);
            if (std::abs(sixj1)<1e-9) continue;

            double sum_12 = 0;
            for (int twoS12=twoS12_min; twoS12<=twoS12_max; twoS12+=2)
            {
              if ( std::abs(2*L-twoS12)>twoJ or 2*L+twoS12<twoJ ) continue; // triangle condition for ninejL
              double ninejL =  NineJ( Lab, lc, L,
                                      S1,  sc, 0.5*twoS12,
                                      Jab, jc, 0.5*twoJ);
              if ( std::abs(ninejL)<1e-9) continue;
              for (int L12=L12_min; L12<=L12_max; L12++)
              {
                double ninej12 = ninej_12[2*(L12-L12_min)+twoS12/2];
                if (std::abs(ninej12)<1e-9) continue;
                // Don't bother with triangle checks here, because the sixj functions already do that
                double sixj2 = AngMom::SixJ( Lcm,L12,L,0.5*twoS12,0.5*twoJ,0.5*twoJ12);
                double sixj3 = AngMom::SixJ( Lcm,L2,Lambda,L1,L,L12);

                sum_12 += ninejL * ninej12 * sixj2 * sixj3;
              } // for L12
            } // for twoS12
            sum_L += sixj1 * sum_12;
          } // for L
          tcoeff += ninejab * mosh1 * mosh2 * sum_L;
        } // for Lambda
      } // for curlyL
    } // for Lab

  // multiply by all the j-hat symbols.   
  // accounting for all the hat factors: 
  //  We also have (2S12+1)(2L12+1) included in ninej12,  (2L+1) included in sixj1,  (2Lab+1) included in ninejab, (2Lambda+1) included in mosh2, and (2curlyL+1) included in mosh1
    // hat factors:  sqrt( (2*ja+1)*(2*jb+1)*(2*jc+1)*(2*Jab+1)*(twoJ12+1)*(2*J1+1)*(twoJ2+1)*(2*S1+1) ) *   (2*Lab+1) * (twoS12+1) (2*L12+1) * (2*L+1) * (2*Lambda+1) 
    //                                                 global                                                (ninejab) * (    ninej12       ) * (sixj1)   (   mosh2  )
    // phases:    ( L12 + S12 + Lcm + J + L1 + L2 - L12 + Lcm + L2 + L1 + L  + curlyL + lc - Lambda + lc + Lab - L + lc + L1 + curlyL + L )
    // reduces to (       S12       + J                                               + lc - Lambda      + Lab - L      + L1              )
    //            (L1 + lc) * ( S12+J ) * (Lab) * (Lambda +L)
    //              global     ninej12   ninejab     sixj1
    //
    // phase( lc + Lab + L + L1 + (twoJ + twoS12)/2 ) * phase( Lambda)
    double jhats = sqrt( (2*ja+1)*(2*jb+1)*(2*jc+1)*(2*Jab+1)*(twoJ12+1)*(2*J1+1)*(twoJ2+1)*(2*S1+1) ) ;
    int globalphase = AngMom::phase(L1+lc);
    return tcoeff * jhats * globalphase;
  }






  double Tcoeff_fancy( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm)
  {
    double ja = 0.5*j2a;
    double sa = 0.5;
    double jb = 0.5*j2b;
    double sb = 0.5;
//    double jc = 0.5*j2c;
//    double sc = 0.5;
    int Eab = 2*(na+nb)+la+lb;

//    double S2 = 0.5;

    int Lab_min = std::max( std::abs(la-lb), std::abs(Jab-S1) );
    int Lab_max = std::min( la+lb , Jab+S1 );

    double tcoeff=0.;
    for (int Lab=Lab_min; Lab<=Lab_max; Lab++)
    {
      double ninejab = (2*Lab+1) * NineJ(la, lb, Lab, sa, sb, S1, ja, jb, Jab ); 
      if (std::abs(ninejab)<1e-8) continue; // <- This may not help very much

      double curlyL_sum = 0;
      int curlyL_min = std::max( std::abs(L1-Lab), std::abs(J1-Jab) );
      int curlyL_max = std::min( L1+Lab, J1+Jab );
      if ( (L1+Lab+curlyL_min)%2 > 0 ) curlyL_min++; // conserve parity in the Moshinsky bracket
      for (int curlyL=curlyL_min; curlyL<=curlyL_max; curlyL+=2)
      {
        int curlyN = ( Eab - 2*N1-L1 - curlyL )/2;
        if (curlyN<0) continue;
        double mosh1 = AngMom::Moshinsky(curlyN,curlyL,N1,L1, na,la, nb,lb, Lab, 1.0);
//        int mosh1phase = AngMom::phase( curlyN+N1+na+nb ); // need to confirm which convention is being used here...
        int mosh1phase = AngMom::phase( curlyL+L1+Lab + (curlyN+N1+na+nb)*1 ); // need to confirm which convention is being used here...

        double Lambda_sum = 0;
        int Lambda_min = std::max( std::abs(Lcm-L2), std::abs(curlyL-lc) );
        int Lambda_max = std::min( Lcm+L2, curlyL+lc );
        for (int Lambda=Lambda_min; Lambda<=Lambda_max; Lambda++)
        {
          double mosh2 = AngMom::Moshinsky( Ncm,Lcm, N2,L2,  curlyN,curlyL, nc,lc, Lambda, 2.0);
//          int mosh2phase = AngMom::phase( Ncm+N2+curlyN+nc);
          int mosh2phase = AngMom::phase( Lcm+L2+Lambda + (Ncm+N2+curlyN+nc)*1 );

          double sixj =  AngMom::SixJ( J1,  curlyL, Jab,
                                       Lab, S1,     L1  );

          if (std::abs(mosh2)<1e-8 or std::abs(sixj)<1e-8) continue;
          double twelvej = AngMom::TwelveJ_1_ints(  twoJ,   j2c,      1,        twoJ2,
                                                    2*Jab,  2*lc,     2*L2,     twoJ12,
                                                    2*J1,   2*curlyL, 2*Lambda, 2*Lcm );

          Lambda_sum +=  AngMom::phase(Lambda)*(2*Lambda+1) * mosh2 * mosh2phase * sixj * twelvej ;
          std::cout << "    inner loop: " << AngMom::phase(Lambda) << " " << 2*Lambda+1 << " " << mosh2 << " " << mosh2phase
                    << " " << sixj << " " << twelvej << "   ->  " << Lambda_sum << std::endl;
        } // for Lambda
        curlyL_sum += AngMom::phase(curlyL) * mosh1 * mosh1phase * Lambda_sum;
        std::cout << "  curlyL_sum: " << AngMom::phase(curlyL) << " " << mosh1 << " " <<  mosh1phase << " " << Lambda_sum << std::endl;
      } // for curlyL
      tcoeff += (2*Lab+1) * ninejab * curlyL_sum;
    } // for Lab
    double hats = sqrt( (2*S1+1)*(j2a+1)*(j2b+1)*(j2c+1)*(2*Jab+1)*(2*J1+1)*(twoJ2+1)*(twoJ12+1) );
    int overallphase = AngMom::phase(L1 + Lcm + S1 + Jab + (twoJ2-1)/2 );
    return tcoeff * hats * overallphase;
  }


} //end namespace AngMom

