

#include "DarkMatterNREFT.hh"
#include "PhysicalConstants.hh"
#include <gsl/gsl_sf_gamma.h>  // for gsl_sf_gamma,  gsl_sf_fact, gsl_sf_doublefact
#include <gsl/gsl_sf_hyperg.h> // for gsl_sf_hyperg_1F1
#include <gsl/gsl_sf_coupling.h> // for gsl_sf_coupling_3j, gsl_sf_coupling_6j

//#ifndef HBARC
// #define HBARC 197.3269718 // hc in MeV * fm
//#endif
//#ifndef M_NUCLEON
// #define M_NUCLEON 938.9185 // average nucleon mass in MeV
//#endif

namespace DM_NREFT
{
  using PhysConst::HBARC;
  using PhysConst::M_NUCLEON;

// From Daniel Gazda's fortran implementation:
//
//  double precision function jho(na0 , la, nb0, lb, L, y)
//    use mod_gsl
//    implicit none
//    integer, intent(in) :: na0, la, nb0, lb, L
//    double precision, intent(in) :: y
//    integer :: k1, k, na, nb
//    double precision :: f1, f2, f3
//
//    ! shift the radial numbers
//    na = na0 + 1
//    nb = nb0 + 1
//    ! since the code below expects n >= 1 and I use n >= 0 in the rest of the code
//
//    jho = 0.0d0
//
//    f1 = 2.d0**L/doublefactorial(2*L+1)* &
//         y**(dble(L)/2.d0)*exp(-y) &
//         *sqrt(factorial(nb-1)*factorial(na-1) &
//         *dgamma(dble(na+la)+0.5d0)*dgamma(dble(nb+lb)+0.5d0))
//
//    do k = 0, nb - 1
//       f2 = dble((-1)**k)/factorial(k)/factorial(nb - 1 - k) &
//            /dgamma(dble(lb+k)+3.d0/2.d0)
//       f3 = 0.0d0
//       do k1 = 0,na-1
//          f3 = f3 + dble((-1)**k1)/factorial(k1)/factorial(na - 1 - k1) &
//               *dgamma((dble(L + lb + la + 2*k + 2*k1 + 3))/2.d0) &
//               /dgamma(dble(la + k1) + 3.0d0/2.0d0) &
//               *hg1f1((dble(L - lb - la - 2*k - 2*k1))/2.d0, dble(L) + 3.0d0/2.0d0, y)
//       end do
//       jho = jho + f2 * f3
//    end do
//    jho = f1 * jho
//  end function jho
//


  double jho( int na0, int la, int nb0, int lb, int L, double y)
  {
     int na = na0 + 1;
     int nb = nb0 + 1;

     double f1 = pow(2.0,L) / gsl_sf_doublefact(2*L+1)
                * pow(y, 0.5*L ) * exp(-y)
                * sqrt( gsl_sf_fact(nb-1)*gsl_sf_fact(na-1) * gsl_sf_gamma(na+la+0.5) * gsl_sf_gamma(nb+lb+0.5) );

     double j_ho = 0.0;
     for (int k=0; k<nb; k++)
     {
       double f2 = pow(-1,k) / gsl_sf_fact(k) / gsl_sf_fact(nb-1-k) / gsl_sf_gamma(lb + k + 1.5);
       double f3 = 0.0;
       for (int k1=0; k1<na; k1++)
       {
         f3 += pow(-1,k1) / gsl_sf_fact(k1) / gsl_sf_fact(na-1-k1)
             * gsl_sf_gamma( 0.5*(L+lb+la+2*k+2*k1+3.0) )
             / gsl_sf_gamma(la + k1 + 1.5)
             * gsl_sf_hyperg_1F1( 0.5*(L - lb - la-2*k-2*k1), L+1.5, y   ) ;
       }
       j_ho += f2*f3;
     }
     j_ho *= f1;
     return j_ho;
  }




//  double precision function jdmho(na0, la, nb0, lb, L, y)
//    use mod_gsl
//    implicit none
//    integer, intent(in) :: na0, la, nb0, lb, L
//    double precision, intent(in) :: y
//    integer :: k1, k, na, nb
//    double precision :: f1, f2, f3
//
//    na = na0 + 1
//    nb = nb0 + 1
//    ! since the code below expects n >= 1
//
//    jdmho = 0.0d0
//
//    f1 = 2.d0**(L-1)/doublefactorial(2*L+1) &
//         *y**(dble(L-1)/2.d0)*exp(-y) &
//         *sqrt(factorial(nb-1)*factorial(na-1) &
//         *dgamma(dble(na+la)+1.d0/2.d0)*dgamma(dble(nb+lb)+0.5d0))
//
//    do k = 0, nb - 1
//       f2 = dble((-1)**k)/factorial(k)/factorial(nb-1-k) / dgamma(dble(lb+k)+3.d0/2.d0)
//       f3 = 0.0d0
//       do k1=0,na-1
//          f3 = f3 + dble((-1)**k1)/factorial(k1)/factorial(na-1-k1) &
//               *dgamma((dble(L+lb+la+2*k+2*k1+2))/2.d0) &
//               /dgamma(dble(la+k1)+3.d0/2.d0) &
//               *(-(dble(L+lb+la+2*k+2*k1+2))/2.d0 &
//               *hg1f1((dble(L-lb-la-2*k-2*k1-1))/2.d0 &
//               ,dble(L)+3.d0/2.d0,y) &
//               +dble(2*k) &
//               *hg1f1((dble(L-lb-la-2*k-2*k1+1))/2.d0, dble(L)+3.d0/2.d0, y))
//       end do
//       jdmho = jdmho + f2 * f3
//    end do
//    jdmho = f1 * jdmho
//  end function jdmho

  double jdmho( int na0, int la, int nb0, int lb, int L, double y)
  {
     int na = na0 + 1;
     int nb = nb0 + 1;

     double f1 = pow(2,L-1) / gsl_sf_doublefact(2*L+1)
               * pow(y, 0.5*(L-1) ) * exp(-y)
                 * sqrt( gsl_sf_fact(nb-1)*gsl_sf_fact(na-1) * gsl_sf_gamma(na+la+0.5) * gsl_sf_gamma(nb+lb+0.5) );

     double j_dmho = 0.0;
     for (int k=0; k<nb; k++)
     {
       double f2 = pow(-1,k) / gsl_sf_fact(k) / gsl_sf_fact(nb-1-k) / gsl_sf_gamma(lb + k + 1.5);
       double f3 = 0.0;
       for (int k1=0; k1<na; k1++)
       {
         f3 += pow(-1,k1) / gsl_sf_fact(k1) / gsl_sf_fact(na-1-k1)
             * gsl_sf_gamma( 0.5*(L+lb+la+2*k+2*k1+2.0) )
             / gsl_sf_gamma(la + k1 + 1.5)
             * ( -(L+lb+la+2*k+2*k1+2) * 0.5 
                * gsl_sf_hyperg_1F1( 0.5*(L - lb - la - 2*k - 2*k1 - 1), L+1.5, y )
                + (2*k)
                * gsl_sf_hyperg_1F1( 0.5*(L - lb - la - 2*k - 2*k1 + 1), L+1.5, y ) // different from two lines above by -1 -> +1.
               );
       }
       j_dmho += f2*f3;
     }
     j_dmho *= f1;
     return j_dmho;

  }



//  double precision function jdpho(na0,la,nb0,lb,L,y)
//    use mod_gsl
//    implicit none
//    integer, intent(in) :: na0, la, nb0, lb, L
//    double precision, intent(in) :: y
//    integer :: k1, k, na, nb
//    double precision :: f1, f2, f3
//
//    na = na0 + 1
//    nb = nb0 + 1
//    ! since the code below expects n_rad >= 1
//
//    jdpho = 0.0d0
//
//    f1 = 2.d0**(L-1)/doublefactorial(2*L+1) &
//         *y**(dble(L-1)/2.d0)*exp(-y) &
//         *sqrt(factorial(nb-1)*factorial(na-1) &
//         *dgamma(dble(na+la)+0.5d0)*dgamma(dble(nb+lb)+0.5d0))
//    do k = 0, nb - 1
//       f2 = dble((-1)**k)/factorial(k)/factorial(nb-1-k) &
//            /dgamma(dble(lb+k)+3.d0/2.d0)
//       f3 = 0.0d0
//       do k1 = 0, na-1
//          f3 = f3+dble((-1)**k1)/factorial(k1)/factorial(na-1-k1) &
//               *dgamma((dble(L+lb+la+2*k+2*k1+2))/2.d0) &
//               /dgamma(dble(la+k1)+3.d0/2.d0) &
//               *(-(dble(L+lb+la+2*k+2*k1+2))/2.d0 &
//               *hg1f1((dble(L-lb-la-2*k-2*k1-1))/2.d0 &
//               ,dble(L)+3.d0/2.d0,y) &
//               +dble(2*lb+2*k+1) &
//               *hg1f1((dble(L-lb-la-2*k-2*k1+1))/2.d0 ,dble(L)+3.d0/2.d0, y))
//       end do
//       jdpho = jdpho+f2*f3
//    end do
//    jdpho = f1*jdpho
//
//  end function jdpho

  double jdpho( int na0, int la, int nb0, int lb, int L, double y)
  {
     int na = na0 + 1;
     int nb = nb0 + 1;

     double f1 = pow(2,L-1) / gsl_sf_doublefact(2*L+1)
               * pow(y, 0.5*(L-1) ) * exp(-y)
                 * sqrt( gsl_sf_fact(nb-1)*gsl_sf_fact(na-1)
                 * gsl_sf_gamma(na+la+0.5) * gsl_sf_gamma(nb+lb+0.5) );

     double j_dpho = 0.0;
     for (int k=0; k<nb; k++)
     {
       double f2 = pow(-1,k) / gsl_sf_fact(k) / gsl_sf_fact(nb-1-k)
                   / gsl_sf_gamma(lb + k + 1.5);
       double f3 = 0.0;
       for (int k1=0; k1<na; k1++)
       {
         f3 += pow(-1,k1) / gsl_sf_fact(k1) / gsl_sf_fact(na-1-k1)
             * gsl_sf_gamma( 0.5*(L+lb+la+2*k+2*k1+2.0) )
             / gsl_sf_gamma(la + k1 + 1.5)
             * ( -(L+lb+la+2*k+2*k1+2) * 0.5 
                * gsl_sf_hyperg_1F1( 0.5*(L - lb - la - 2*k - 2*k1 - 1), L+1.5, y )
                + (2*lb+2*k+1)
                * gsl_sf_hyperg_1F1( 0.5*(L - lb - la - 2*k - 2*k1 + 1), L+1.5, y ) // different from two lines above by -1 -> +1.
               );
       }
       j_dpho += f2*f3;
     }
     j_dpho *= f1;
     return j_dpho;

  }




// From Daniel Gazda's fortran implementation:
//   br( j ) = sqrt( 2J+1)
//   br2( j ) = J+1
//   wnj is the Wigner n-j symbol
//
//  double precision function M(na, la, ja2, nb, lb, jb2, J, y)
//    ! < na la ja2/2 | M_J (q r) | nb lb jb2/2 >
//    use mod_parameters, only: pi
//    use mod_jho, only: jho
//    use mod_gsl, only: w3j, w6j
//    implicit none
//    integer, intent(in) :: na, la, ja2, nb, lb, jb2, J
//    double precision, intent(in) :: y
//
//    M = 1.0d0/sqrt(4.d0*pi)*(-1.d0)**(dble(2*J+jb2+1)/2.d0) &
//         *br(la) * br(lb) * br2(ja2) * br2(jb2) * br(J) &
//         *w6j(2*la, ja2, 1, jb2, 2*lb, 2*J) &
//         *w3j(2*la, 2*J, 2*lb, 0, 0, 0) &
//         *jho(na, la, nb, lb, 2*J, y)
//
//  end function M


  Operator M( ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    double b2 =  1.0 / (modelspace.GetHbarOmega() * M_NUCLEON) ; // b^2 in MeV^-2
    double y = q*q*b2/4.0;
    int Tz = 0;
    int parity = J%2; // normal parity operator
//    int isofactor = IsoSV=="+" ? 1 : -1; // "+" labels isoscalar and "-" is isovector
    std::array<int,2> isofactor = {1,1};  // isoscalar, proton + neutron.
    if (IsoSV =="-" ) isofactor[1] = -1;  // isovector, neutrons get a minus sign.
    else if (IsoSV == "p") isofactor[1] = 0; // proton only,  neutrons don't contribute.
    else if (IsoSV == "n") isofactor[0] = 0; // neutron only, protons don't contribute.
    std::cout << "In " << __func__ << " and isofactor is " << isofactor[0] << " " << isofactor[1] << std::endl;
    Operator M_op( modelspace, J, Tz, parity, 2);
    int norb = modelspace.GetNumberOrbits();
    for (int a=0; a<norb; a++)
    {
      Orbit& oa = modelspace.GetOrbit(a);
      int na  = oa.n;
      int la  = oa.l;
      int j2a = oa.j2;
//      for (int b=0; b<norb; b++)
      for (int b : M_op.OneBodyChannels.at({la,j2a,oa.tz2}) )
      {
        Orbit& ob = modelspace.GetOrbit(b);
        int nb  = ob.n;
        int lb  = ob.l;
        int j2b = ob.j2;

        double mab = 1.0 / sqrt(4*M_PI)* pow(-1,(2*J+j2b+1)*0.5)
                     * sqrt( (2*la+1.) * (2*lb+1.) * (j2a+1.) * (j2b+1.) * (2*J+1.) )
                     * gsl_sf_coupling_6j(2*la, j2a, 1, j2b, 2*lb, 2*J )
                     * gsl_sf_coupling_3j(2*la, 2*J, 2*lb, 0,0,0)
//bhu                     * jho( na, la, nb, lb, 2*J, y );
                     * jho( na, la, nb, lb, J, y );
//        M_op.OneBody(a,b) = mab * pow( isofactor, (oa.tz2+1)/2 );
	// This formula assumes we're storing things as reduced matrix elements. If J=0, then we aren't, so we get a factor sqrt(2j+1)
	if (J==0 and parity==0)
	{
	   mab /= sqrt( j2a+1.0);
        }
        M_op.OneBody(a,b) = mab  * isofactor[(oa.tz2+1)/2];
      }
    }

    return M_op;

  }



//  double precision function Ms(na, la, ja2, nb, lb, jb2, J, L, y)
//    ! < na la ja2/2 | M_JL (q r) | nb lb jb2/2 >
//    use mod_parameters, only: pi
//    use mod_jho, only: jho
//    use mod_gsl, only: w3j,w9j
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J,L
//    double precision, intent(in) :: y
//
//    Ms = sqrt(3.d0/(2.d0*pi))*(-1.d0)**la &
//         *br(la)*br(lb)*br2(ja2)*br2(jb2)*br(L)*br(J) &
//         *w9j(2*la,2*lb,2*L,1,1,2,ja2,jb2,2*J) &
//         *w3j(2*la,2*L,2*lb,0,0,0) &
//         *jho(na,la,nb,lb,L,y)
//
//  end function Ms

  Operator Ms( ModelSpace& modelspace, std::string IsoSV, int J, int L, double q )
  {
    double b2 =  1.0 / (modelspace.GetHbarOmega() * M_NUCLEON) ;
    double y = q*q*b2/4.0;
    int Tz = 0;
    int parity = L%2; // normal parity operator
//    int isofactor = IsoSV=="+" ? 1 : -1; // "+" labels isoscalar and "-" is isovector
    std::array<int,2> isofactor = {1,1};  // isoscalar, proton + neutron.
    if (IsoSV =="-" ) isofactor[1] = -1;  // isovector, neutrons get a minus sign.
    else if (IsoSV == "p") isofactor[1] = 0; // proton only,  neutrons don't contribute.
    else if (IsoSV == "n") isofactor[0] = 0; // neutron only, protons don't contribute.
    Operator Ms_op( modelspace, J, Tz, parity, 2);
    int norb = modelspace.GetNumberOrbits();
    for (int a=0; a<norb; a++)
    {
      Orbit& oa = modelspace.GetOrbit(a);
      int na  = oa.n;
      int la  = oa.l;
      int j2a = oa.j2;
//      for (int b=0; b<norb; b++)
      for (int b : Ms_op.OneBodyChannels.at({la,j2a,oa.tz2}) )
      {
        Orbit& ob = modelspace.GetOrbit(b);
        int nb  = ob.n;
        int lb  = ob.l;
        int j2b = ob.j2;

        double mab =  sqrt(1.5/M_PI)* pow(-1,la)
                     * sqrt( (2*la+1.) * (2*lb+1.) * (j2a+1.) * (j2b+1.) * (2*L+1) * (2*J+1.) )
                     * gsl_sf_coupling_9j(2*la, 2*lb, 2*L, 1, 1, 2, j2a, j2b, 2*J )
                     * gsl_sf_coupling_3j(2*la, 2*L, 2*lb, 0,0,0)
                     * jho( na, la, nb, lb, L, y );
	// This formula assumes we're storing things as reduced matrix elements. If J=0, then we aren't, so we get a factor sqrt(2j+1)
	if (J==0 and parity==0)
	{
	   mab /= sqrt( j2a+1.0);
        }
        Ms_op.OneBody(a,b) = mab * isofactor[(oa.tz2+1)/2];
      }
    }

    return Ms_op;
  }



//  double precision function Mg(na, la, ja2, nb, lb, jb2, J, L, y)
//    ! < na la ja2/2 | M_JL (q r)  1/q grad | nb lb jb2/2 >
//    use mod_parameters, only: pi
//    use mod_jho, only: jho,jdmho,jdpho
//    use mod_gsl, only: w3j,w6j,w9j
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J,L
//    double precision, intent(in) :: y
//
//    Mg = 0.0d0
//    if (lb == 0) return ! function gives zero even for undefined lb=0 !!!
//    if (L < 0) return ! and also for nonsense L < 0 !!!
//
//    Mg = 1.d0/sqrt(4.d0*pi)*(-1.d0)**(dble(2*L+jb2+1)/2.d0) &
//         *br(la)*br2(ja2)*br2(jb2)*br(L)*br(J) &
//         *w6j(2*la,ja2,1,jb2,2*lb,2*J) &
//         *(-sqrt(dble(lb+1))*br(lb+1) &
//         *w6j(2*L,2,2*J,2*lb,2*la,2*(lb+1)) &
//         *w3j(2*la,2*L,2*(lb+1),0,0,0) &
//         *jdmho(na,la,nb,lb,L,y) &
//         +sqrt(dble(lb))*br(lb-1) &
//         *w6j(2*L,2,2*J,2*lb,2*la,2*(lb-1)) &
//         *w3j(2*la,2*L,2*(lb-1),0,0,0) &
//         *jdpho(na,la,nb,lb,L,y))
//  end function Mg


  Operator Mg( ModelSpace& modelspace, std::string IsoSV, int J, int L, double q )
  {
    double b2 =  1.0 / (modelspace.GetHbarOmega() * M_NUCLEON) ;
    double y = q*q*b2/4.0;
    int Tz = 0;
    int parity = (L+1)%2; // abnormal parity operator
//    int isofactor = IsoSV=="+" ? 1 : -1; // "+" labels isoscalar and "-" is isovector
    std::array<int,2> isofactor = {1,1};  // isoscalar, proton + neutron.
    if (IsoSV =="-" ) isofactor[1] = -1;  // isovector, neutrons get a minus sign.
    else if (IsoSV == "p") isofactor[1] = 0; // proton only,  neutrons don't contribute.
    else if (IsoSV == "n") isofactor[0] = 0; // neutron only, protons don't contribute.
    Operator Mg_op(modelspace, J, Tz, parity, 2);
    int norb = modelspace.GetNumberOrbits();
    for (int a=0; a<norb; a++)
    {
      Orbit& oa = modelspace.GetOrbit(a);
      int na  = oa.n;
      int la  = oa.l;
      int j2a = oa.j2;
      //for (int b=0; b<norb; b++)
      for (int b : Mg_op.OneBodyChannels.at({la,j2a,oa.tz2}) )
      {
        Orbit& ob = modelspace.GetOrbit(b);
        int nb  = ob.n;
        int lb  = ob.l;
        int j2b = ob.j2;
        if (la==0) continue;
        if (lb==0) continue;

        double mab = 1.0 / sqrt(4.0*M_PI) * pow( -1, (2*L+j2b+1.0)*0.5 )
                     * sqrt( (2*la+1) * (j2a+1) * (j2b+1) * (2*L+1) * (2*J+1) )
                     * gsl_sf_coupling_6j( 2*la, j2a, 1, j2b, 2*lb, 2*J)
                     * (- sqrt( lb+1.0) * sqrt( 2*(lb+1)+1 )
                         * gsl_sf_coupling_6j(2*L,2,2*J,2*lb,2*la,2*(lb+1))
                         * gsl_sf_coupling_3j(2*la,2*L,2*(lb+1),0,0,0)
                         * jdmho(na,la,nb,lb,L,y)
                        + sqrt(lb) * sqrt( 2*(lb-1)+1 )
                         * gsl_sf_coupling_6j(2*L,2,2*J,2*lb,2*la,2*(lb-1))
                         * gsl_sf_coupling_3j(2*la,2*L,2*(lb-1),0,0,0)
                         * jdpho(na,la,nb,lb,L,y)
                        );
	// This formula assumes we're storing things as reduced matrix elements. If J=0, then we aren't, so we get a factor sqrt(2j+1)
	if (J==0 and parity==0)
	{
	   mab /= sqrt( j2a+1.0);
        }
        Mg_op.OneBody(a,b) = mab * isofactor[(oa.tz2+1)/2];
      }
    }

    return Mg_op;
  }


//  double precision function Sigma(na,la,ja2,nb,lb,jb2,J,y)
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    Sigma=Ms(na,la,ja2,nb,lb,jb2,J,J,y)
//  end function Sigma

  Operator Sigma( ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    Operator Sigma_op = Ms( modelspace, IsoSV, J, J, q);
    Sigma_op.SetAntiHermitian();
    return Sigma_op;
  }



//  double precision function Sigmap(na,la,ja2,nb,lb,jb2,J,y)
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    Sigmap = -sqrt(dble(J))/br(J) * Ms(na,la,ja2,nb,lb,jb2,J,J+1,y)
//    if (J == 0) return
//    Sigmap = Sigmap + sqrt(dble(J+1))/br(J) * Ms(na,la,ja2,nb,lb,jb2,J,J-1,y)
//  end function Sigmap

  Operator Sigmap( ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    Operator Sigmap_op = -sqrt(J/(2*J+1.)) * Ms( modelspace, IsoSV, J, J+1, q);
    if ( J > 0 )  Sigmap_op += sqrt( (J+1.)/(2*J+1.) ) * Ms(modelspace, IsoSV, J, J-1, q);
    return Sigmap_op;
  }





//  double precision function Sigmapp(na,la,ja2,nb,lb,jb2,J,y)
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    Sigmapp = sqrt(dble(J+1))/br(J) * Ms(na,la,ja2,nb,lb,jb2,J,J+1,y)
//    if (.not. J == 0) then
//       Sigmapp = Sigmapp + sqrt(dble(J))/br(J) * Ms(na,la,ja2,nb,lb,jb2,J,J-1,y)
//    end if
//  end function Sigmapp

  Operator Sigmapp( ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    Operator Sigmapp_op = sqrt( (J+1.0)/(2*J+1.) ) * Ms( modelspace, IsoSV, J, J+1, q);
    if ( J > 0 )  Sigmapp_op += sqrt(J/(2*J+1.)) * Ms( modelspace, IsoSV, J, J-1, q);
    return Sigmapp_op;
  }



//  double precision function Delta(na,la,ja2,nb,lb,jb2,J,y)
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    Delta = Mg(na,la,ja2,nb,lb,jb2,J,J,y)
//  end function Delta


  Operator Delta( ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    return Mg(modelspace, IsoSV, J, J, q);
  }



//  double precision function Deltap(na,la,ja2,nb,lb,jb2,J,y)
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    Deltap = -sqrt(dble(J))/br(J)*Mg(na,la,ja2,nb,lb,jb2,J,J+1,y) &
//         +sqrt(dble(J+1))/br(J)*Mg(na,la,ja2,nb,lb,jb2,J,J-1,y)
//  end function Deltap


  Operator Deltap(  ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    Operator Deltap_op = -sqrt(J/(2*J+1.)) * Mg(modelspace, IsoSV, J, J+1, q);
    if ( J > 0 )  Deltap_op += sqrt( (J+1.)/(2*J+1.) ) * Mg(modelspace, IsoSV, J, J-1, q);
    Deltap_op.SetAntiHermitian();
    return Deltap_op;
  }



//  double precision function PhiF(la,ja2,jb2)
//    use mod_parameters, only: pi
//    implicit none
//    integer, intent(in) :: la,ja2,jb2
//    PhiF = (-1.d0)**(la+1)*6.d0*br(la)*br2(ja2)*br2(jb2)/sqrt(4.d0*pi)
//  end function PhiF


  double PhiF( int la, int j2a, int j2b )
  {
    return pow(-1,la+1) * 6 * sqrt( (2*la+1) * (j2a+1) * (j2b+1)) / sqrt(4*M_PI);
  }




//  double precision function PhiS1(na,la,ja2,nb,lb,jb2,J,y)
//    use mod_jho, only: jdmho
//    use mod_gsl, only: w3j,w6j,w9j
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    integer :: lam
//
//    PhiS1=0.0d0
//    do lam = J, J + 1
//       PhiS1=PhiS1+(-1.d0)**dble(J-lam+1)*br(lam)**2 &
//            *w6j(2*(J+1),2,2*lam,2,2*J,2) &
//            *w6j(2*(J+1),2,2*lam,2*lb,2*la,2*(lb+1)) &
//            *w9j(2*la,2*lb,2*lam,1,1,2,ja2,jb2,2*J)
//    end do
//
//    PhiS1 = PhiS1*br(lb+1)*sqrt(dble(lb+1)) &
//         *w3j(2*la,2*(J+1),2*(lb+1),0,0,0) &
//         *jdmHO(na,la,nb,lb,J+1,y)
//
//  end function PhiS1



  double PhiS1( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y )
  {
    double phis1 = 0.0;
    for (int lam = J; lam<=J+1; lam++ )
    {
      phis1 += pow(-1, J-lam+1) * sqrt( (2*lam+1)*(2*lam+1) )
               * gsl_sf_coupling_6j(2*(J+1),2,2*lam,2,2*J,2)
               * gsl_sf_coupling_6j(2*(J+1),2,2*lam,2*lb,2*la,2*(lb+1))
               * gsl_sf_coupling_9j(2*la,2*lb,2*lam,1,1,2,j2a,j2b,2*J);
    }

    phis1 *= sqrt(2*(lb+1)+1)*sqrt(lb+1)
             * gsl_sf_coupling_3j(2*la,2*(J+1),2*(lb+1),0,0,0)
             * jdmho(na,la,nb,lb,J+1,y);

    return phis1;
  }




//  double precision function PhiS2(na,la,ja2,nb,lb,jb2,J,y)
//    use mod_jho, only: jdpho
//    use mod_gsl, only: w3j,w6j,w9j
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    integer :: lam
//
//    PhiS2=0.0d0
//    if (lb == 0) return
//
//    do lam=J,J+1
//       PhiS2=PhiS2+(-1.d0)**dble(J-lam)*br(lam)**2 &
//            *w6j(2*(J+1),2,2*lam,2,2*J,2) &
//            *w6j(2*(J+1),2,2*lam,2*lb,2*la,2*(lb-1)) &
//            *w9j(2*la,2*lb,2*lam,1,1,2,ja2,jb2,2*J)
//    end do
//
//    PhiS2=PhiS2*br(lb-1)*sqrt(dble(lb)) &
//         *w3j(2*la,2*(J+1),2*(lb-1),0,0,0) &
//         *jdpHO(na,la,nb,lb,J+1,y)
//
//  end function PhiS2



  double PhiS2( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y )
  {
    double phis2 = 0.0;
    if (lb==0) return phis2;

    for (int lam = J; lam<=J+1; lam++ )
    {
      phis2 += pow(-1, J-lam) * sqrt( (2*lam+1)*(2*lam+1) )
               * gsl_sf_coupling_6j(2*(J+1),2,2*lam,2,2*J,2)
               * gsl_sf_coupling_6j(2*(J+1),2,2*lam,2*lb,2*la,2*(lb-1))
               * gsl_sf_coupling_9j(2*la,2*lb,2*lam,1,1,2,j2a,j2b,2*J);
    }

    phis2 *= sqrt(2*(lb-1)+1)*sqrt(lb)
             * gsl_sf_coupling_3j(2*la,2*(J+1),2*(lb-1),0,0,0)
             * jdpho(na,la,nb,lb,J+1,y);

    return phis2;
  }


//  double precision function PhiS3(na,la,ja2,nb,lb,jb2,J,y)
//    use mod_jho, only: jdmho
//    use mod_gsl, only: w3j,w6j,w9j
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    integer :: lam
//
//    PhiS3 = 0.0d0
//    if (J == 0) return
//
//    do lam=J-1,J
//       PhiS3 = PhiS3+(-1.d0)**dble(J-lam+1)*br(lam)**2 &
//            *w6j(2*(J-1),2,2*lam,2,2*J,2) &
//            *w6j(2*(J-1),2,2*lam,2*lb,2*la,2*(lb+1)) &
//            *w9j(2*la,2*lb,2*lam,1,1,2,ja2,jb2,2*J)
//    end do
//
//    PhiS3 = PhiS3*br(lb+1)*sqrt(dble(lb+1)) &
//         *w3j(2*la,2*(J-1),2*(lb+1),0,0,0) &
//         *jdmHO(na,la,nb,lb,J-1,y)
//
//  end function PhiS3


  double PhiS3( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y )
  {
    double phis3 = 0.0;
    if (J==0) return phis3;

    for (int lam = J-1; lam<=J; lam++ )
    {
      phis3 += pow(-1, J-lam+1) * sqrt( (2*lam+1)*(2*lam+1) )
               * gsl_sf_coupling_6j(2*(J-1),2,2*lam,2,2*J,2)
               * gsl_sf_coupling_6j(2*(J-1),2,2*lam,2*lb,2*la,2*(lb+1))
               * gsl_sf_coupling_9j(2*la,2*lb,2*lam,1,1,2,j2a,j2b,2*J);
    }

    phis3 *= sqrt(2*(lb+1)+1)*sqrt(lb+1)
             * gsl_sf_coupling_3j(2*la,2*(J-1),2*(lb+1),0,0,0)
             * jdmho(na,la,nb,lb,J-1,y);

    return phis3;
  }


//  double precision function PhiS4(na,la,ja2,nb,lb,jb2,J,y)
//    use mod_jho, only: jdpho
//    use mod_gsl, only: w3j,w6j,w9j
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    integer :: lam
//
//    PhiS4=0.0d0
//    if (J == 0 .or. lb == 0) return
//
//    do lam=J-1,J
//       PhiS4 = PhiS4+(-1.d0)**dble(J-lam)*br(lam)**2 &
//            *w6j(2*(J-1),2,2*lam,2,2*J,2) &
//            *w6j(2*(J-1),2,2*lam,2*lb,2*la,2*(lb-1)) &
//            *w9j(2*la,2*lb,2*lam,1,1,2,ja2,jb2,2*J)
//    end do
//
//    PhiS4 = PhiS4*br(lb-1)*sqrt(dble(lb)) &
//         *w3j(2*la,2*(J-1),2*(lb-1),0,0,0) &
//         *jdpHO(na,la,nb,lb,J-1,y)
//  end function PhiS4

  double PhiS4( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y )
  {
    double phis4 = 0.0;
    if (J==0 or lb==0) return phis4;

    for (int lam = J-1; lam<=J; lam++ )
    {
      phis4 += pow(-1, J-lam) * sqrt( (2*lam+1)*(2*lam+1) )
               * gsl_sf_coupling_6j(2*(J-1),2,2*lam,2,2*J,2)
               * gsl_sf_coupling_6j(2*(J-1),2,2*lam,2*lb,2*la,2*(lb-1))
               * gsl_sf_coupling_9j(2*la,2*lb,2*lam,1,1,2,j2a,j2b,2*J);
    }

    phis4 *= sqrt(2*(lb-1)+1)*sqrt(lb)
             * gsl_sf_coupling_3j(2*la,2*(J-1),2*(lb-1),0,0,0)
             * jdpho(na,la,nb,lb,J-1,y);

    return phis4;
  }





//  double precision function Phip(na,la,ja2,nb,lb,jb2,J,y)
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//
//    Phip = PhiF(la,ja2,jb2) * ( &
//         -br(J+1)*sqrt(dble(J)) &
//         *(PhiS1(na,la,ja2,nb,lb,jb2,J,y) &
//         +PhiS2(na,la,ja2,nb,lb,jb2,J,y)))
//
//    if (J == 0) return
//
//    Phip = Phip + PhiF(la,ja2,jb2) * ( &
//         br(J-1)*sqrt(dble(J+1)) &
//         *(PhiS3(na,la,ja2,nb,lb,jb2,J,y) &
//         +PhiS4(na,la,ja2,nb,lb,jb2,J,y)))
//
//  end function Phip


  Operator Phip( ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {

    double b2 =  1.0 / (modelspace.GetHbarOmega() * M_NUCLEON) ;
    double y = q*q*b2/4.0;
    int Tz = 0;
    int parity = J%2; // normal parity operator
//    int isofactor = IsoSV=="+" ? 1 : -1; // "+" labels isoscalar and "-" is isovector
    std::array<int,2> isofactor = {1,1};  // isoscalar, proton + neutron.
    if (IsoSV =="-" ) isofactor[1] = -1;  // isovector, neutrons get a minus sign.
    else if (IsoSV == "p") isofactor[1] = 0; // proton only,  neutrons don't contribute.
    else if (IsoSV == "n") isofactor[0] = 0; // neutron only, protons don't contribute.
    Operator Phip_op(modelspace, J, Tz, parity, 2);
    Phip_op.SetAntiHermitian();
    int norb = modelspace.GetNumberOrbits();
    for (int a=0; a<norb; a++)
    {
      Orbit& oa = modelspace.GetOrbit(a);
      int na  = oa.n;
      int la  = oa.l;
      int j2a = oa.j2;
//      for (int b=0; b<norb; b++)
      for ( int b : Phip_op.OneBodyChannels.at({la,j2a,oa.tz2}) )
      {
        Orbit& ob = modelspace.GetOrbit(b);
        int nb  = ob.n;
        int lb  = ob.l;
        int j2b = ob.j2;

        double phip_ab = PhiF(la,j2a,j2b) * ( -sqrt(2*(J+1)+1) * sqrt(J)
                             * ( PhiS1(na,la,j2a,nb,lb,j2b,J,y) + PhiS2(na,la,j2a,nb,lb,j2b,J,y) ) );

        if (J>0) phip_ab += PhiF(la,j2a,j2b) * ( sqrt(2*(J-1)+1) * sqrt(J+1)
                             * ( PhiS3(na,la,j2a,nb,lb,j2b,J,y) + PhiS4(na,la,j2a,nb,lb,j2b,J,y) ) );

	// This formula assumes we're storing things as reduced matrix elements. If J=0, then we aren't, so we get a factor sqrt(2j+1)
	if (J==0 and parity==0)
	{
	   phip_ab /= sqrt( j2a+1.0);
        }
//        Phip_op.OneBody(a,b) = phip_ab * pow( isofactor, (oa.tz2+1)/2 );
        Phip_op.OneBody(a,b) = phip_ab * isofactor[(oa.tz2+1)/2];
      }
    }
    return Phip_op;
  }


//   if (normal_conditions(la,ja2,lb,jb2,J)) then
//      ME1b%Phitp(ia,ib,J) = Phip(na,la,ja2,nb,lb,jb2,J,y) &
//           +Ms(na,la,ja2,nb,lb,jb2,J,J,y)/2.d0


  Operator Phitp( ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    Operator Phitp_op = Phip(modelspace,IsoSV,J,q) + 0.5 * Ms(modelspace, IsoSV, J, J, q);
    return Phitp_op;
  }





//double precision function Phipp(na,la,ja2,nb,lb,jb2,J,y)
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//
//    Phipp = PhiF(la,ja2,jb2) &
//         *(br(J+1)*sqrt(dble(J+1)) &
//         *(PhiS1(na,la,ja2,nb,lb,jb2,J,y) &
//         +PhiS2(na,la,ja2,nb,lb,jb2,J,y)))
//
//    if (.not. J == 0) then
//
//       Phipp = Phipp + PhiF(la,ja2,jb2) &
//            *(br(J-1)*sqrt(dble(J)) &
//            *(PhiS3(na,la,ja2,nb,lb,jb2,J,y) &
//            + PhiS4(na,la,ja2,nb,lb,jb2,J,y)))
//
//    end if
//
//  end function Phipp


  Operator Phipp(   ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    double b2 =  1.0 / (modelspace.GetHbarOmega() * M_NUCLEON) ;
    double y = q*q*b2/4.0;
    int Tz = 0;
    int parity = J%2; // normal parity operator
//    int isofactor = IsoSV=="+" ? 1 : -1; // "+" labels isoscalar and "-" is isovector
    std::array<int,2> isofactor = {1,1};  // isoscalar, proton + neutron.
    if (IsoSV =="-" ) isofactor[1] = -1;  // isovector, neutrons get a minus sign.
    else if (IsoSV == "p") isofactor[1] = 0; // proton only,  neutrons don't contribute.
    else if (IsoSV == "n") isofactor[0] = 0; // neutron only, protons don't contribute.
    Operator Phipp_op(modelspace, J, Tz, parity, 2);
    int norb = modelspace.GetNumberOrbits();
    for (int a=0; a<norb; a++)
    {
      Orbit& oa = modelspace.GetOrbit(a);
      int na  = oa.n;
      int la  = oa.l;
      int j2a = oa.j2;
      // IsoSV = "p" => proton, "n" => neutron, "+" => proton + neutron,   "-" => proton - neutron
//      for (int b=0; b<norb; b++)
      for ( int b : Phipp_op.OneBodyChannels.at({la,j2a,oa.tz2}) )
      {
        Orbit& ob = modelspace.GetOrbit(b);
        int nb  = ob.n;
        int lb  = ob.l;
        int j2b = ob.j2;

        double phipp_ab = PhiF(la,j2a,j2b) * ( sqrt(2*(J+1)+1) * sqrt(J+1)
                             * ( PhiS1(na,la,j2a,nb,lb,j2b,J,y) + PhiS2(na,la,j2a,nb,lb,j2b,J,y) ) );

        if (J>0) phipp_ab += PhiF(la,j2a,j2b) * ( sqrt(2*(J-1)+1) * sqrt(J)
                             * ( PhiS3(na,la,j2a,nb,lb,j2b,J,y) + PhiS4(na,la,j2a,nb,lb,j2b,J,y) ) );

        // This formula assumes we're storing things as reduced matrix elements. If J=0, then we aren't, so we get a factor sqrt(2j+1)
	if (J==0 and parity==0)
	{
	   phipp_ab /= sqrt( j2a+1.0);
	}
        Phipp_op.OneBody(a,b) = phipp_ab *  isofactor[(oa.tz2+1)/2];
      }
    }
    return Phipp_op;
  }




//  double precision function Omega(na,la,ja2,nb,lb,jb2,J,y)
//    use mod_parameters, only: pi
//    use mod_gsl, only: w3j,w6j
//    use mod_jho, only: jdmho,jdpho
//    implicit none
//    integer, intent(in) :: na,la,ja2,nb,lb,jb2,J
//    double precision, intent(in) :: y
//    double precision :: O1,O2
//    if (jb2 == 2*lb+1) then
//       O1=-jdmho(na,la,nb,lb,J,y)
//    else
//       O1=0.0d0
//    end if
//
//    if (jb2 == 2*lb-1) then
//       O2=jdpho(na,la,nb,lb,J,y)
//    else
//       O2=0.0d0
//    end if
//    Omega = 1.d0/sqrt(4.d0*pi)*(-1.d0)**la &
//         *br(la)*br2(ja2)*br2(jb2)*br(jb2-lb)*br(J) &
//         *w6j(2*la,ja2,1,jb2,2*(jb2-lb),2*J) &
//         *w3j(2*la,2*J,2*(jb2-lb),0,0,0) &
//         *(O1+O2)
//  end function Omega


  Operator Omega(  ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    double b2 =  1.0 / (modelspace.GetHbarOmega() * M_NUCLEON) ;
    double y = q*q*b2/4.0;
    int Tz = 0;
    int parity = (J+1)%2; // abnormal parity operator
//    int isofactor = IsoSV=="+" ? 1 : -1; // "+" labels isoscalar and "-" is isovector
    std::array<int,2> isofactor = {1,1};  // isoscalar, proton + neutron.
    if (IsoSV =="-" ) isofactor[1] = -1;  // isovector, neutrons get a minus sign.
    else if (IsoSV == "p") isofactor[1] = 0; // proton only,  neutrons don't contribute.
    else if (IsoSV == "n") isofactor[0] = 0; // neutron only, protons don't contribute.
    Operator Omega_op(modelspace, J, Tz, parity, 2);
    std::cout << "Warning: Omega is neither Symmetrize nor AntiSymmetrize !!!" << std::endl;
    int norb = modelspace.GetNumberOrbits();
    for (int a=0; a<norb; a++)
    {
      Orbit& oa = modelspace.GetOrbit(a);
      int na  = oa.n;
      int la  = oa.l;
      int j2a = oa.j2;
//      for (int b=0; b<norb; b++)
      for ( int b : Omega_op.OneBodyChannels.at({la,j2a,oa.tz2}) )
      {
        Orbit& ob = modelspace.GetOrbit(b);
        int nb  = ob.n;
        int lb  = ob.l;
        int j2b = ob.j2;

        double O1=0.0;
        double O2=0.0;
        if (j2b==2*lb+1)
        {
          O1 = -jdmho(na,la,nb,lb,J,y);
        }
        else if (j2b==2*lb-1)
        {
          O2 = jdpho(na,la,nb,lb,J,y);
        }
        else
        {
           continue;
        }
        double omega_ab = 1.0/sqrt(4*M_PI) * modelspace.phase(la)
                          * sqrt( (2*la+1) * (j2a+1) * (j2b+1) * (2*(j2b-lb)+1) * (2*J+1) )
                          * gsl_sf_coupling_6j(2*la,j2a,1,j2b,2*(j2b-lb),2*J)
                          * gsl_sf_coupling_3j(2*la,2*J,2*(j2b-lb),0,0,0)
                          * (O1+O2);
//        Omega_op.OneBody(a,b) = omega_ab * pow( isofactor, (oa.tz2+1)/2 );
        // This formula assumes we're storing things as reduced matrix elements. If J=0, then we aren't, so we get a factor sqrt(2j+1)
	if (J==0 and parity==0)
	{
	   omega_ab /= sqrt( j2a+1.0);
	}
        Omega_op.OneBody(a,b) = omega_ab * isofactor[(oa.tz2+1)/2];
      }
    }
    return Omega_op;
  }








  Operator Omegat(  ModelSpace& modelspace, std::string IsoSV, int J, double q )
  {
    Operator Omegat_op = Omega(modelspace,IsoSV,J,q) + 0.5 * Sigmapp(modelspace, IsoSV, J, q);
    Omegat_op.SetAntiHermitian();
    return Omegat_op;
  }






} // namespace DM_NREFT

