
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
   return gsl_sf_coupling_6j(int(2*j1),int(2*j2),int(2*j3),int(2*J1),int(2*J2),int(2*J3));
 }

 /// Wigner 9j symbol
 double NineJ(double j1,double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
 {
   return gsl_sf_coupling_9j(int(2*j1),int(2*j2),int(2*J12),int(2*j3),int(2*j4),int(2*J34),int(2*J13),int(2*J24),int(2*J));
 }

 /// Normalized Wigner 9j symbol
 double NormNineJ(double j1,double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
 {
   return sqrt( (2*J12+1)*(2*J34+1)*(2*J13+1)*(2*J24+1) ) * NineJ(j1,j2,J12,j3,j4,J34,J13,J24,J);
 }





///  Talmi-Moshinsky Bracket, using the algorthm of Buck et al. Nuc. Phys. A 600 (1996) 387-402
///  Their phase convention differs from Moshinsky's by a factor \f$(-1)^{L+l+\Lambda}\f$.
///  I correct for this in order to stick to Moshinsky's phase convention.
//
double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lam, double B)
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

                double mosh1 = AngMom::Moshinsky(curlyN,curlyL, N1,L1, na,la, nb,lb, Lab, MOSH_BETA_1);
                double mosh2 = AngMom::Moshinsky(Ncm,Lcm, N2,L2, curlyN,curlyL, nc,lc, Lambda, MOSH_BETA_2);
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



//                std::cout << "      ##" << L12 << "## " <<phase_123 << " " << hats_all << " " << phasemosh1 << " " << mosh1 << " " << phasemosh2 << " " << mosh2 << " " << ninej1 << " " << ninej2 << " " << ninej3 << " " << sixj1 << " " << sixj2 << " " << sixj3 << "   :  "
//                          << tcoeff << "  "
//                          << std::endl;
//                std::cout << "       hat-breakdown: sqrt( " <<  (j2a+1)<< " * " <<(j2b+1)<< " * " <<(j2c+1)<< " * " <<(2*Jab+1)<< " * " <<(twoJ+1)<< " * " <<(2*J1+1)<< " * " <<(twoJ2+1)<< " * " << (2*S1+1) << " )  * "
//                          << (twoS12+1) << " * " << (2*Lab+1) << " * " << (2*L+1) << " * " << (2*curlyL+1) << " * " << (2*Lambda+1) << std::endl;
//                }
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


    double tcoeff=0.;
    for (int Lab=Lab_min; Lab<=Lab_max; Lab++)
    {
      double ninejab = (2*Lab+1) * NineJ(la, lb, Lab, sa, sb, S1, ja, jb, Jab );
      int L12_min = std::max((twoJ12-2*S1-1)/2, std::abs(L1-L2) );
      int L12_max = std::min((twoJ12+2*S1+1)/2, L1+L2 );
      for (int L12=L12_min; L12<=L12_max; L12+=1)
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
                              * NineJ(Lab,lc,L,S1,sc,0.5*twoS12,Jab,jc,0.5*twoJ)
                              * SixJ(Lcm, L12,L,0.5*twoS12,0.5*twoJ,0.5*twoJ12) * phase( lc + Lab + L + L1 + (twoJ + twoS12)/2 );
          }
          if (std::abs(t12)<1e-9) continue;

          double tmosh=0;
          int Lambda_min = std::max( std::abs(Lcm-L2), std::abs(L1-L));
          int Lambda_max = std::min( Lcm+L2, L1+L );
          for (int Lambda=Lambda_min; Lambda<=Lambda_max; Lambda++)
          {
            double sixjLambda = (2*Lambda+1) * SixJ(Lcm,L2,Lambda,L1,L,L12);
            
            int curlyL_min = std::max( std::abs(lc-Lambda), std::abs(L1-Lab) );
            int curlyL_max = std::min( lc+Lambda, L1+Lab );
            for (int curlyL=curlyL_min; curlyL<=curlyL_max; curlyL+=1)
            {
              int curlyN = na+nb-N1  +(la+lb-L1-curlyL)/2;
              if (curlyN<0) continue;
              if ( (la+lb+L1+curlyL)%2>0 ) continue;
              if ( 2*curlyN+curlyL +2*nc+lc != 2*Ncm+Lcm + 2*N2+L2) continue;
              double sixjCurly =  SixJ(lc,curlyL,Lambda,L1,L,Lab) ;  // I think the 2*curlyL+1 in eq (10) in the Roth paper isn't supposed to be there
              // increment moshinsky loop bit
              int moshphase1 = AngMom::phase( curlyL+L1+Lab + 0*(curlyL+L1 - lb-la)/2 );  // converting phase conventions
              int moshphase2 = AngMom::phase( Lcm+L2+Lambda + 0*(Lcm+L2 - curlyL-lc)/2 ); // converting phase conventions (comment at the beginning of the function)
              tmosh += sixjLambda * sixjCurly * phase( Lambda) * moshphase1 * moshphase2
                       * Moshinsky( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab,  MOSH_BETA_1)
                       * Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  MOSH_BETA_2);

            }
          }
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








} //end namespace AngMom

