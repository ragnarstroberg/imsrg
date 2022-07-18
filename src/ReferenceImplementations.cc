
#include "ReferenceImplementations.hh"
#include "ModelSpace.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"

/// Straightforward implementation of J-coupled commutator expressions
/// without optimizations. This should be benchmarked against the
/// mscheme implementation and then left untouched.
namespace ReferenceImplementations{

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    Z0 =  sum_ab (2ja+1) (na-nb) (Xab Yba )
///                   =  sum_ab (2ja+1) na (1-nb) (Xab Yba - Yab Xba)
///
///
void comm110ss( const Operator& X, const Operator& Y, Operator& Z )  
{
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;
   double z0 = 0;
   for ( size_t a : Z.modelspace->all_orbits )
   {
      Orbit& oa = Z.modelspace->GetOrbit(a);
      for ( size_t b : Z.modelspace->all_orbits )
      {
         Orbit& ob = Z.modelspace->GetOrbit(b);
         z0 += (oa.j2+1) * oa.occ *(1-ob.occ) *  (  X1(a,b) * Y1(b,a) - Y1(a,b) * X1(b,a)  ) ;
      }
   }
   Z.ZeroBody += z0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression : Z0 = 1/4 sum_abcd sum_J (2J+1) ( n_a n_b nbar_c nbar_d ) ( x_abcd y_cdab - y_abcd x_cdab )
///
///
void comm220ss( const Operator& X, const Operator& Y, Operator& Z ) 
{
   auto& X2 = X.TwoBody;
   auto& Y2 = Y.TwoBody;
   double z0 = 0;

   for ( size_t a : Z.modelspace->all_orbits )
   {
      Orbit& oa = Z.modelspace->GetOrbit(a);
      for ( size_t b : Z.modelspace->all_orbits )
      {
         Orbit& ob = Z.modelspace->GetOrbit(b);

         for ( size_t c : Z.modelspace->all_orbits )
         {
            Orbit& oc = Z.modelspace->GetOrbit(c);
            for ( size_t d : Z.modelspace->all_orbits )
            {
               Orbit& od = Z.modelspace->GetOrbit(d);

               int Jmin = std::max(  std::abs( oa.j2-ob.j2), std::abs( oc.j2-od.j2) )/2;
               int Jmax = std::min( oa.j2+ob.j2,  oc.j2+od.j2)/2;
               for (int J=Jmin; J<=Jmax; J++)
               {
                  double xabcd = X2.GetTBME_J( J,J, a,b,c,d);
                  double xcdab = X2.GetTBME_J( J,J, c,d,a,b);
                  double yabcd = Y2.GetTBME_J( J,J, a,b,c,d);
                  double ycdab = Y2.GetTBME_J( J,J, c,d,a,b);
                  z0 += 1./4 * (2*J+1) * ( oa.occ * ob.occ * (1-oc.occ)*(1-od.occ))  * ( xabcd * ycdab - yabcd * xcdab);
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
void comm111ss( const Operator& X, const Operator& Y, Operator& Z ) 
{
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;
   auto& Z1 = Z.OneBody;
   for ( size_t i : Z.modelspace->all_orbits )
   {
      for ( size_t j : Z.modelspace->all_orbits )
      {
         for ( size_t a : Z.modelspace->all_orbits )
         {
            Z1(i,j) += X1(i,a) * Y1(a,j) - Y1(i,a) * X1(a,j);
         }
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    Zij = sum_ab sum_J (2J+1)/(2ji+1) (na-nb) ( Xab YJbiaj - Yab XJbiaj )
///
///
void comm121ss( const Operator& X, const Operator& Y, Operator& Z ) 
{
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;
   auto& X2 = X.TwoBody;
   auto& Y2 = Y.TwoBody;
   auto& Z1 = Z.OneBody;

   for ( size_t i : Z.modelspace->all_orbits )
   {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    for ( size_t j : Z.GetOneBodyChannel( oi.l,oi.j2,oi.tz2)  )
    {
        Orbit& oj = Z.modelspace->GetOrbit(j);

        for ( size_t a : Z.modelspace->all_orbits )
        {
          Orbit& oa = Z.modelspace->GetOrbit(a);
          for ( size_t b : Z.modelspace->all_orbits )
          {
             Orbit& ob = Z.modelspace->GetOrbit(b);
             int Jmin = std::max( std::abs(ob.j2-oi.j2) , std::abs(oa.j2-oj.j2) )/2;
             int Jmax = std::min( ob.j2+oi.j2 , oa.j2+oj.j2 ) / 2;
             for ( int J=Jmin; J<=Jmax; J++)
             {
               double xbiaj = X2.GetTBME_J(J,J, b,i,a,j);
               double ybiaj = Y2.GetTBME_J(J,J, b,i,a,j);
               Z1(i,j) += (2*J+1)/(oi.j2+1.) * (oa.occ-ob.occ) * ( X1(a,b)*ybiaj - Y1(a,b)*xbiaj );
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
void comm221ss( const Operator& X, const Operator& Y, Operator& Z ) 
{
   auto& X2 = X.TwoBody;
   auto& Y2 = Y.TwoBody;
   auto& Z1 = Z.OneBody;

   for ( size_t i : Z.modelspace->all_orbits )
   {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    for ( size_t j : Z.GetOneBodyChannel( oi.l,oi.j2,oi.tz2)  )
    {
        Orbit& oj = Z.modelspace->GetOrbit(j);


        for ( size_t a : Z.modelspace->all_orbits )
        {
          Orbit& oa = Z.modelspace->GetOrbit(a);
          for ( size_t b : Z.modelspace->all_orbits )
          {
             Orbit& ob = Z.modelspace->GetOrbit(b);
             for ( size_t c : Z.modelspace->all_orbits )
             {
               Orbit& oc = Z.modelspace->GetOrbit(c);

//               if ( (oi.l+oc.l + oa.l+ob.l)%2>0) continue;
//               if ( (oi.tz2+oc.tz2) != (oa.tz2+ob.tz2) ) continue;

               int Jmin = std::max( std::abs(ob.j2-oa.j2) , std::abs(oc.j2-oi.j2) )/2;
               int Jmax = std::min( ob.j2+oa.j2 , oc.j2+oi.j2 ) / 2;
               for ( int J=Jmin; J<=Jmax; J++)
               {
                 double xciab = X2.GetTBME_J(J,J, c,i,a,b);
                 double yciab = Y2.GetTBME_J(J,J, c,i,a,b);
                 double xabcj = X2.GetTBME_J(J,J, a,b,c,j);
                 double yabcj = Y2.GetTBME_J(J,J, a,b,c,j);
                 Z1(i,j) += 1./2 * (2*J+1)/(oi.j2+1.) * (oa.occ * ob.occ * (1-oc.occ) + (1-oa.occ)*(1-ob.occ)*oc.occ) * ( xciab*yabcj - yciab*xabcj );

               }//J
             }//c
          }//b
        }//a


    }//j
   }//i

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    ZJijkl = sum_a (Xia YJajkl + Xja YJiakl - YJijal Xak - YJijka Xal) -   {X<->Y}
///
///
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

}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    ZJijkl = 1/2 sum_ab ((1-na)(1-nb) - na nb)  (XJijab YJabkl - YJijab XJabkl)
///
///

void comm222_pp_hhss( const Operator& X, const Operator& Y, Operator& Z)
{
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
         for ( int iket_ab=0; iket_ab<nkets; iket_ab++ )
         {
           Ket& ket_ab = tbc.GetKet(iket_ab);
           size_t a = ket_ab.p;
           size_t b = ket_ab.q;
           Orbit& oa = Z.modelspace->GetOrbit(a);
           Orbit& ob = Z.modelspace->GetOrbit(b);
           double xijab = X2.GetTBME(ch,ch, bra, ket_ab);
           double yijab = Y2.GetTBME(ch,ch, bra, ket_ab);
           double xabkl = X2.GetTBME(ch,ch, ket_ab, ket);
           double yabkl = Y2.GetTBME(ch,ch, ket_ab, ket);
           double flip_factor = (a==b) ? 1 : 2; // looping over kets only gets a<=b. So we need a factor of 2 for the other ordering.
           zijkl += 1./2 * flip_factor * ( (1-oa.occ)*(1-ob.occ) - oa.occ * ob.occ ) * ( xijab * yabkl  - yijab * xabkl );
         }// iket_ab

         // Need to normalize here, because AddToTBME expects a normalized TBME.
         if (i==j) zijkl /= PhysConst::SQRT2;
         if (k==l) zijkl /= PhysConst::SQRT2;

         Z2.AddToTBME(ch,ch, ibra,iket, zijkl);

       }// iket
     }// ibra
   }// ch
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    ZJijkl =   sum_ab (na-nb) sum_J'  (2J'+1) (1- (-1)^{i+j-J} Pij) { i j J  } ( ~XJ'ilab ~YJ'abkj - ~YJ'ilab ~XJ'abkj)
///                                                                                 { k l J' }
///           in this expression, the ~X and ~Y are the Pandya-transformed matrix elements
///           ~XJabcd = - sum_J' (2J'+1) { a b J  } XJ'_adcb
///                                      { c d J' }
///
void comm222_phss( const Operator& X, const Operator& Y, Operator& Z)
{
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
       Orbit& oi = Z.modelspace->GetOrbit(i);
       Orbit& oj = Z.modelspace->GetOrbit(j);

       for (int iket=ibra; iket<nkets; iket++)
       {
         Ket& ket = tbc.GetKet(iket);
         size_t k = ket.p;
         size_t l = ket.q;
         Orbit& ok = Z.modelspace->GetOrbit(k);
         Orbit& ol = Z.modelspace->GetOrbit(l);

         double zijkl = 0;
//         int Jpmin = std::max( std::abs(oi.j2-ol.j2) , std::abs(oj.j2-ok.j2) )/2;
//         int Jpmax = std::min( oi.j2+ol.j2 , oj.j2+ok.j2 )/2;
         int Jpmin = std::min(  std::max( std::abs(oi.j2-ol.j2) , std::abs(oj.j2-ok.j2)   ),
                                std::max( std::abs(oj.j2-ol.j2) , std::abs(oi.j2-ok.j2)   ) )/2;
         int Jpmax = std::max(  std::min( oi.j2+ol.j2 , oj.j2+ok.j2 ),
                                std::min( oj.j2+ol.j2 , oi.j2+ok.j2 )  )/2;

         for (int Jp=Jpmin; Jp<=Jpmax; Jp++ )
         {

           double zbar_ilkj = 0;
           double zbar_jlki = 0;
           double sixj_ijkl = AngMom::SixJ( oi.j2*0.5, oj.j2*0.5, J,  ok.j2*0.5, ol.j2*0.5, Jp);
           double sixj_jikl = AngMom::SixJ( oj.j2*0.5, oi.j2*0.5, J,  ok.j2*0.5, ol.j2*0.5, Jp);

           for (size_t a : Z.modelspace->all_orbits)
           {
             Orbit& oa = Z.modelspace->GetOrbit(a);
             for (size_t b : Z.modelspace->all_orbits)
             {
                Orbit& ob = Z.modelspace->GetOrbit(b);

                // "Direct" term
                double xbar_ilab=0;
                double ybar_ilab=0;

                int Jppmin = std::max( std::abs(oi.j2-ob.j2) , std::abs(oa.j2-ol.j2) )/2;
                int Jppmax = std::min( oi.j2+ob.j2 , oa.j2+ol.j2 )/2;
                for (int Jpp=Jppmin; Jpp<=Jppmax; Jpp++)
                {
                  double sixj_ilab = AngMom::SixJ( oi.j2*0.5, ol.j2*0.5, Jp,  oa.j2*0.5, ob.j2*0.5, Jpp);
                  double xibal = X2.GetTBME_J(Jpp, i,b,a,l);
                  double yibal = Y2.GetTBME_J(Jpp, i,b,a,l);
                  xbar_ilab -= (2*Jpp+1) * sixj_ilab * xibal;
                  ybar_ilab -= (2*Jpp+1) * sixj_ilab * yibal;
                }


                double xbar_abkj=0;
                double ybar_abkj=0;
                Jppmin = std::max( std::abs(oa.j2-oj.j2) , std::abs(ok.j2-ob.j2) )/2;
                Jppmax = std::min( oa.j2+oj.j2 , ok.j2+ob.j2 )/2;
                for (int Jpp=Jppmin; Jpp<=Jppmax; Jpp++)
                {
                  double sixj_abkj = AngMom::SixJ( oa.j2*0.5, ob.j2*0.5, Jp,  ok.j2*0.5, oj.j2*0.5, Jpp);
                  double xajkb = X2.GetTBME_J(Jpp, a,j,k,b);
                  double yajkb = Y2.GetTBME_J(Jpp, a,j,k,b);
                  xbar_abkj -= (2*Jpp+1) * sixj_abkj * xajkb;
                  ybar_abkj -= (2*Jpp+1) * sixj_abkj * yajkb;
                }



                zbar_ilkj += (oa.occ - ob.occ) * (xbar_ilab * ybar_abkj - ybar_ilab * xbar_abkj);


                // "Exchange" term obtained by exchanging i<->j.
                double xbar_jlab=0;
                double ybar_jlab=0;
                Jppmin = std::max( std::abs(oj.j2-ob.j2) , std::abs(oa.j2-ol.j2) )/2;
                Jppmax = std::min( oj.j2+ob.j2 , oa.j2+ol.j2 )/2;
                for (int Jpp=Jppmin; Jpp<=Jppmax; Jpp++)
                {
                  double sixj_jlab = AngMom::SixJ( oj.j2*0.5, ol.j2*0.5, Jp,  oa.j2*0.5, ob.j2*0.5, Jpp);
                  double xjbal = X2.GetTBME_J(Jpp, j,b,a,l);
                  double yjbal = Y2.GetTBME_J(Jpp, j,b,a,l);
                  xbar_jlab -= (2*Jpp+1) * sixj_jlab * xjbal;
                  ybar_jlab -= (2*Jpp+1) * sixj_jlab * yjbal;
                }


                double xbar_abki=0;
                double ybar_abki=0;
                Jppmin = std::max( std::abs(oa.j2-oi.j2) , std::abs(ok.j2-ob.j2) )/2;
                Jppmax = std::min( oa.j2+oi.j2 , ok.j2+ob.j2 )/2;
                for (int Jpp=Jppmin; Jpp<=Jppmax; Jpp++)
                {
                  double sixj_abki = AngMom::SixJ( oa.j2*0.5, ob.j2*0.5, Jp,  ok.j2*0.5, oi.j2*0.5, Jpp);
                  double xaikb = X2.GetTBME_J(Jpp, a,i,k,b);
                  double yaikb = Y2.GetTBME_J(Jpp, a,i,k,b);
                  xbar_abki -= (2*Jpp+1) * sixj_abki * xaikb;
                  ybar_abki -= (2*Jpp+1) * sixj_abki * yaikb;
                }


                zbar_jlki += (oa.occ - ob.occ) * (xbar_jlab * ybar_abki - ybar_jlab * xbar_abki);



             }//b
           }//a
           int phase_ij = AngMom::phase( (oi.j2 + oj.j2 - 2*J)/2 );
           zijkl += (2*Jp+1) * sixj_ijkl * zbar_ilkj;  // Direct
           zijkl -= (2*Jp+1) * sixj_jikl * zbar_jlki * phase_ij;   //Exchange, with phase


         }// Jp


         // Need to normalize here, because AddToTBME expects a normalized TBME.
         if (i==j) zijkl /= PhysConst::SQRT2;
         if (k==l) zijkl /= PhysConst::SQRT2;

         Z2.AddToTBME(ch,ch, ibra,iket, zijkl);

       }// iket
     }// ibra
   }// ch

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// In the optimized version, an intermediate is generated
///  and used for both, but here we just call both sequentially.
///
void comm222_pp_hh_221ss( const Operator& X, const Operator& Y, Operator& Z)
{
  comm222_pp_hhss(X,Y,Z);
  comm221ss(X,Y,Z);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    Z0 = 1/36 sum_abcdef sum_J1J2J na nb nc (1-nd)(1-ne)(1-nf) (2J+1) (XJ1J2J_abcdef YJ1J2J_defabc - YJ1J2J_abcdef XJ1J2J_defabc)
///
///
void comm330ss( const Operator& X, const Operator& Y, Operator& Z )
{
   auto& X3 = X.ThreeBody;
   auto& Y3 = Y.ThreeBody;
   double z0 = 0;

   size_t norb = Z.modelspace->GetNumberOrbits();
//   for ( size_t a : Z.modelspace->all_orbits )
   #pragma omp parallel for schedule(dynamic,1) reduction(+:z0)
   for ( size_t a=0; a<norb; a++ )
   {
    Orbit& oa = Z.modelspace->GetOrbit(a);
    for ( size_t b : Z.modelspace->all_orbits )
    {
       Orbit& ob = Z.modelspace->GetOrbit(b);

       for ( size_t c : Z.modelspace->all_orbits )
       {
         Orbit& oc = Z.modelspace->GetOrbit(c);
         if ( std::abs ( oa.occ * ob.occ * oc.occ) <1e-7) continue;
         for ( size_t d : Z.modelspace->all_orbits )
         {
           Orbit& od = Z.modelspace->GetOrbit(d);

           for ( size_t e : Z.modelspace->all_orbits )
           {
             Orbit& oe = Z.modelspace->GetOrbit(e);
             for ( size_t f : Z.modelspace->all_orbits )
             {
                Orbit& of = Z.modelspace->GetOrbit(f);

                if ( (oa.l+ob.l+oc.l + od.l+oe.l+of.l)%2>0 ) continue;
                if ( (oa.tz2+ob.tz2+oc.tz2) != (od.tz2+oe.tz2+of.tz2) ) continue;
                if ( std::abs ( oa.occ * ob.occ * oc.occ* (1-od.occ)*(1-oe.occ)*(1-of.occ)) <1e-9) continue;

                int J1min = std::abs( oa.j2-ob.j2)/2;
                int J2min = std::abs( od.j2-oe.j2)/2;
                int J1max = (oa.j2+ob.j2)/2;
                int J2max = (od.j2+oe.j2)/2;

                for (int J1=J1min; J1<=J1max; J1++)
                {
                  for (int J2=J2min; J2<=J2max; J2++)
                  {
                    int twoJmin = std::max( std::abs( oc.j2-2*J1) , std::abs(of.j2-2*J2));
                    int twoJmax = std::min( oc.j2+2*J1 , of.j2 + 2*J2);
                    for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
                    {
                       double xabcdef = X3.GetME_pn( J1,J2, twoJ, a,b,c,d,e,f);
                       double yabcdef = Y3.GetME_pn( J1,J2, twoJ, a,b,c,d,e,f);
                       double xdefabc = X3.GetME_pn( J2,J1, twoJ, d,e,f,a,b,c);
                       double ydefabc = Y3.GetME_pn( J2,J1, twoJ, d,e,f,a,b,c);
                      z0 += 1./36 * (twoJ+1) * ( oa.occ * ob.occ * oc.occ* (1-od.occ)*(1-oe.occ)*(1-of.occ))  * ( xabcdef * ydefabc - yabcdef * xdefabc);

                    }// twoJ
                 }//J2
               }//J1
             }//f
           }//e
         }//d
       }//c
    }//b
   }//a
   Z.ZeroBody += z0;
}        




//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    Zij = 1/12 sum_abcde sum_J1J2J na nb (1-nc) (1-nd)(1-ne) (2J+1)/(2ji+1) (XJ1J2J_abicde YJ1J2J_cdeabj - YJ1J2J_abicde XJ1J2J_cdeabj)
///
///
void comm331ss( const Operator& X, const Operator& Y, Operator& Z )
{
   auto& X3 = X.ThreeBody;
   auto& Y3 = Y.ThreeBody;
   auto& Z1 = Z.OneBody;

   size_t norb = Z.modelspace->GetNumberOrbits();
//   for ( size_t i : Z.modelspace->all_orbits )
   #pragma omp parallel for schedule(dynamic,1)
   for ( size_t i=0; i<norb; i++ )
   {
     Orbit& oi = Z.modelspace->GetOrbit(i);
//     for ( size_t j : Z.modelspace->all_orbits )
     for ( size_t j : Z.GetOneBodyChannel(oi.l, oi.j2,oi.tz2) )
     {
       Orbit& oj = Z.modelspace->GetOrbit(j);
       double zij = 0;
       for ( size_t a : Z.modelspace->all_orbits )
       {
         Orbit& oa = Z.modelspace->GetOrbit(a);
         for ( size_t b : Z.modelspace->all_orbits )
         {
          Orbit& ob = Z.modelspace->GetOrbit(b);

          int J1min = std::abs( oa.j2-ob.j2)/2;
          int J1max = (oa.j2+ob.j2)/2;

          for ( size_t c : Z.modelspace->all_orbits )
          {
            Orbit& oc = Z.modelspace->GetOrbit(c);
            for ( size_t d : Z.modelspace->all_orbits )
            {
              Orbit& od = Z.modelspace->GetOrbit(d);
              int J2min = std::abs( oc.j2-od.j2)/2;
              int J2max = (oc.j2+od.j2)/2;

              for ( size_t e : Z.modelspace->all_orbits )
              {
                Orbit& oe = Z.modelspace->GetOrbit(e);
                if ( (oa.l+ob.l+oi.l + oc.l+od.l+oe.l + Z.parity)%2 >0) continue;
                if ( std::abs(oa.tz2+ob.tz2+oi.tz2 - oc.tz2-od.tz2-oe.tz2) != Z.GetTRank()) continue;
                if ( std::abs( oa.occ * ob.occ *(1-oc.occ)*(1-od.occ)*(1-oe.occ) )<1e-8 ) continue;

                for (int J1=J1min; J1<=J1max; J1++)
                {
                  if ( a==b and J1%2>0) continue;

                  for (int J2=J2min; J2<=J2max; J2++)
                  {
                    if ( c==d and J2%2>0) continue;
                    int twoJmin = std::max( std::abs( 2*J1-oi.j2) , std::abs( 2*J2-oe.j2) );
                    int twoJmax = std::min(  2*J1+oi.j2 ,  2*J2+oe.j2 );
                    for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
                    {
                      double xabicde = X3.GetME_pn( J1, J2, twoJ, a,b,i,c,d,e);
                      double yabicde = Y3.GetME_pn( J1, J2, twoJ, a,b,i,c,d,e);
                      double xcdeabj = X3.GetME_pn( J2, J1, twoJ, c,d,e,a,b,j);
                      double ycdeabj = Y3.GetME_pn( J2, J1, twoJ, c,d,e,a,b,j);
                      zij += 1./12 * oa.occ * ob.occ *(1-oc.occ)*(1-od.occ)*(1-oe.occ) * (twoJ+1.)/(oi.j2+1) * ( xabicde * ycdeabj - yabicde * xcdeabj );
                    }               
                  }//Jj
                }//J1
              }//e
            }//d
          }//c
         }//b
       }//a
       Z1(i,j) += zij;
     }//j
   }//i
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    Zij = 1/4 sum_abcd sum_J2J na nb (1-nc) (1-nd) (2J+1)/(2ji+1) (XJ1_abcd YJ1J1J_cdiabj + XJ1J2J_abicdj YJ1_cdab    - YJ1_abcd XJ1J1J_cdiabj - YJ1J1J_abicdj XJ1_cdab)
///
///
void comm231ss( const Operator& X, const Operator& Y, Operator& Z )
{
   auto& X2 = X.TwoBody;
   auto& Y2 = Y.TwoBody;
   auto& X3 = X.ThreeBody;
   auto& Y3 = Y.ThreeBody;
   auto& Z1 = Z.OneBody;

   size_t norb = Z.modelspace->GetNumberOrbits();
//   for ( size_t i : Z.modelspace->all_orbits )
//   #pragma omp parallel for schedule(dynamic,1)
   for ( size_t i=0; i<norb; i++ )
   {
     Orbit& oi = Z.modelspace->GetOrbit(i);
//     for ( size_t j : Z.modelspace->all_orbits )
     for ( size_t j : Z.GetOneBodyChannel(oi.l, oi.j2,oi.tz2) )
     {
       Orbit& oj = Z.modelspace->GetOrbit(j);
       double zij = 0;
       for ( size_t a : Z.modelspace->all_orbits )
       {
         Orbit& oa = Z.modelspace->GetOrbit(a);
         for ( size_t b : Z.modelspace->all_orbits )
         {
          Orbit& ob = Z.modelspace->GetOrbit(b);

          int J1min = std::abs( oa.j2-ob.j2)/2;
          int J1max = (oa.j2+ob.j2)/2;

          for ( size_t c : Z.modelspace->all_orbits )
          {
            Orbit& oc = Z.modelspace->GetOrbit(c);
            for ( size_t d : Z.modelspace->all_orbits )
            {
              Orbit& od = Z.modelspace->GetOrbit(d);
              int J2min = std::abs( oc.j2-od.j2)/2;
              int J2max = (oc.j2+od.j2)/2;

//                if ( (oa.l+ob.l+oi.l + oc.l+od.l+oj.l + Z.parity)%2 >0) continue;
//                if ( std::abs(oa.tz2+ob.tz2+oi.tz2 - oc.tz2-od.tz2-oe.tz2) != Z.GetTRank()) continue;
                if ( std::abs( oa.occ * ob.occ *(1-oc.occ)*(1-od.occ) )<1e-8 ) continue;

                for (int J1=J1min; J1<=J1max; J1++)
                {
                  if ( a==b and J1%2>0) continue;
                  if ( c==d and J1%2>0) continue;

                  double xabcd = X2.GetTBME_J( J1,J1, a,b,c,d);
                  double xcdab = X2.GetTBME_J( J1,J1, c,d,a,b);
                  double yabcd = Y2.GetTBME_J( J1,J1, a,b,c,d);
                  double ycdab = Y2.GetTBME_J( J1,J1, c,d,a,b);

                  int twoJmin = std::max( std::abs( 2*J1-oi.j2) , std::abs( 2*J1-oj.j2) );
                  int twoJmax = std::min(  2*J1+oi.j2 ,  2*J1+oj.j2 );
                  for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
                  {
                    double xabicdj = X3.GetME_pn( J1, J1, twoJ, a,b,i,c,d,j);
                    double yabicdj = Y3.GetME_pn( J1, J1, twoJ, a,b,i,c,d,j);
                    double xcdiabj = X3.GetME_pn( J1, J1, twoJ, c,d,i,a,b,j);
                    double ycdiabj = Y3.GetME_pn( J1, J1, twoJ, c,d,i,a,b,j);
                    zij += 1./4 * oa.occ * ob.occ *(1-oc.occ)*(1-od.occ) * (twoJ+1.)/(oi.j2+1) * ( xabcd * ycdiabj + xabicdj *ycdab - yabcd * xcdiabj - yabicdj * xcdab );
                  }//J
                }//J1
            }//d
          }//c
         }//b
       }//a
       Z1(i,j) += zij;
     }//j
   }//i
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    ZJ1ijkl =  sum_ab sum_J (na-nb) sum_J  (2J+1)/(2J1+1) ( Xab YJ1J1J_ijbkla - Yab XJ1J1J_ijbkla )
///
///
void comm132ss( const Operator& X, const Operator& Y, Operator& Z )
{
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;
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
       for (int iket=ibra; iket<nkets; iket++)
       {
         Ket& ket = tbc.GetKet(iket);
         size_t k=ket.p;
         size_t l=ket.q;
         double zijkl = 0;
         for ( size_t a : Z.modelspace->all_orbits)
         {
           Orbit& oa = Z.modelspace->GetOrbit(a);
           for ( size_t b : Z.modelspace->all_orbits)
           {
              Orbit& ob = Z.modelspace->GetOrbit(b);
              double xab = X1(a,b);
              double yab = Y1(a,b);
              if (  std::abs(xab)<1e-8  and std::abs(yab)<1e-8) continue;
              int twoJmin = std::max( std::abs(2*J1-oa.j2) , std::abs(2*J1-ob.j2) );
              int twoJmax = std::min( 2*J1+oa.j2 , 2*J1+ob.j2 );
              for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
              {
                double xijbkla = X3.GetME_pn( J1,J1,twoJ, i,j,b,k,l,a);
                double yijbkla = Y3.GetME_pn( J1,J1,twoJ, i,j,b,k,l,a);
                zijkl +=  (oa.occ - ob.occ) * (twoJ+1.)/(2*J1+1) * ( X1(a,b) * yijbkla - Y1(a,b)*xijbkla);
              }//twoJ
           }//b
         }//a
         zijkl /= sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
         Z2.AddToTBME(ch2,ch2, ibra,iket, zijkl);
       }//iket
     }//ibra
   }//ch2

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    ZJ1ijkl = -1/2 sum_abc sum_J2J3 ( nanb(1-nc) + (1-na)(1-nb)nc ) sqrt{ (2J2+1)/(2J1+1)} (2J3+1) 
///                                 [ (1-PJ1_ij) { i j  J1 } ( XJ2_cjab YJ2J1J3_abiklc - YJ2_cjab XJ2J1J3_abiklc )
///                                              { c J3 J2 }
///                                  -(1-PJ1_kl) { k l  J1 } ( YJ1J23_ijcabk XJ2_abcl - XJ1J2J3_ijcabk YJ2_abcl  ]
///                                              { c J3 J2 }
///
void comm232ss( const Operator& X, const Operator& Y, Operator& Z )
{}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    ZJ1ijkl = 1/6 sum_abcd sum_J2J3  (na nb nc(1-nd) + (1-na)(1-nb)(1-nc)nd) (2J3+1)/(2J1+1)
///                                       *  ( XJ1J2J3_ijdabc YJ2J1J3_abckld  - YJ1J2J3_ijdabc XJ2J1J3_abckld )
///                                 
void comm332_ppph_hhhpss( const Operator& X, const Operator& Y, Operator& Z )
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
         for ( size_t a : Z.modelspace->all_orbits)
         {
           Orbit& oa = Z.modelspace->GetOrbit(a);
           for ( size_t b : Z.modelspace->all_orbits)
           {
             Orbit& ob = Z.modelspace->GetOrbit(b);

             for ( size_t c : Z.modelspace->all_orbits)
             {
               Orbit& oc = Z.modelspace->GetOrbit(c);
               for ( size_t d : Z.modelspace->all_orbits)
               {
                 Orbit& od = Z.modelspace->GetOrbit(d);

                 if ( std::abs(oa.occ*ob.occ*oc.occ*(1-od.occ) + (1-oa.occ)*(1-ob.occ)*(1-oc.occ)*od.occ )<1e-7 ) continue;
                 if ( (oi.l+oj.l+od.l+oa.l+ob.l+oc.l)%2>0) continue;
                 if ( (oi.tz2+oj.tz2+od.tz2)-(oa.tz2+ob.tz2+oc.tz2) !=0 ) continue;

                 int J2min = std::abs(oa.j2-ob.j2)/2;
                 int J2max = (oa.j2+ob.j2)/2;
                 for (int J2=J2min; J2<=J2max; J2++)
                 {

                   int twoJmin = std::max( std::abs(2*J1-od.j2) , std::abs(2*J2-oc.j2) );
                   int twoJmax = std::min( 2*J1+od.j2 , 2*J2+oc.j2 );
                   for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
                   {
                     double xijdabc = X3.GetME_pn( J1,J2,twoJ, i,j,d,a,b,c);
                     double yijdabc = Y3.GetME_pn( J1,J2,twoJ, i,j,d,a,b,c);
                     double xabckld = X3.GetME_pn( J2,J1,twoJ, a,b,c,k,l,d);
                     double yabckld = Y3.GetME_pn( J2,J1,twoJ, a,b,c,k,l,d);
                     zijkl += 1./6 * ( oa.occ*ob.occ*oc.occ*(1-od.occ) - (1-oa.occ)*(1-ob.occ)*(1-oc.occ)*od.occ ) * (twoJ+1.)/(2*J1+1)
                                 * ( xijdabc *yabckld - yijdabc * xabckld);
                   }//twoJ
                 }//J2
               }//d
             }//c
           }//b
         }//a
         zijkl /= sqrt((1.+bra.delta_pq())*(1.+ket.delta_pq()));
         Z2.AddToTBME(ch2,ch2, ibra,iket, zijkl);
       }//iket
     }//ibra
   }//ch2

}

/// This one is nasty.
void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z )
{}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Expression:    ZJ1J2J3_ijklmn =  sum_a  PJ1J3(ij/k) ( X_ka YJ1J2J3_ijalmn - Y_ka XJ1J2J3_ijalmn)
///                                       - PJ2J3(lm/n) ( YJ1J2J3_ijklma X_an - XJ1J2J3_ijklma Y_an )
///  
///  
void comm133ss( const Operator& X, const Operator& Y, Operator& Z )
{
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  auto& X1 = X.OneBody;
  auto& Y1 = Y.OneBody;
  int norbs = Z.modelspace->GetNumberOrbits();
  size_t nch3 = Z.modelspace->GetNumberThreeBodyChannels();
//  #pragma omp parallel for schedule(dynamic,1) 
  for (size_t ch3=0; ch3<=nch3; ch3++)
  {
    auto& Tbc = Z.modelspace->GetThreeBodyChannel(ch3);
    int twoJ = Tbc.twoJ;
    size_t nkets = Tbc.GetNumberKets();
    for (size_t ibra=0; ibra<nkets; ibra++)
    {
      Ket3& bra = Tbc.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      Orbit& ok = Z.modelspace->GetOrbit(k);
      int Jij = bra.Jpq;
      for (size_t iket=ibra; iket<nkets; iket++)
      {
        Ket3& ket = Tbc.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        Orbit& ol = Z.modelspace->GetOrbit(l);
        Orbit& om = Z.modelspace->GetOrbit(m);
        Orbit& on = Z.modelspace->GetOrbit(n);
        int Jlm = ket.Jpq;
        
        double zsum =0;
        // First, connect on the bra side
        for (auto a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
        {
          zsum += X1(i,a) * Y3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
          zsum -= Y1(i,a) * X3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
        }
        for (auto a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
        {
          zsum += X1(j,a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
          zsum -= Y1(j,a) * X3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
        }
        for (auto a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
        {
          zsum += X1(k,a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
          zsum -= Y1(k,a) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
        }

        // Now connect on the ket side
        for (auto a : X.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
        {
          zsum -= X1(a,l) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
          zsum += Y1(a,l) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
        }
        for (auto a : X.OneBodyChannels.at({om.l,om.j2,om.tz2}) )
        {
          zsum -= X1(a,m) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
          zsum += Y1(a,m) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
        }
        for (auto a : X.OneBodyChannels.at({on.l,on.j2,on.tz2}) )
        {
          zsum -= X1(a,n) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
          zsum += Y1(a,n) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
        }
  
//        Z3.AddToME_pn(Jij, Jlm, twoJ, i,j,k,l,m,n, zsum );
        Z3.AddToME_pn_ch(ch3,ch3, ibra, iket, zsum );

      }// for iket
    }// for ibra
  }// for ch3
    
}


void comm223ss( const Operator& X, const Operator& Y, Operator& Z )
{}



} // namespace ReferenceImplementations
////////////////////////////////////////////////////////////////////

