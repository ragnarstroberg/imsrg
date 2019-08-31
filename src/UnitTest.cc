
#include "UnitTest.hh"
#include "AngMom.hh"
#include <armadillo>
#include <random>
#include <string>


uint64_t UnitTest::random_seed = 1;

UnitTest::UnitTest(ModelSpace& ms)
 : modelspace( & ms )
{
}


Operator UnitTest::RandomOp( ModelSpace& modelspace, int jrank, int tz, int parity, int particle_rank, int hermitian )
{
  Operator Rando(modelspace, jrank, tz, parity, particle_rank);
  if ( hermitian == -1 ) Rando.SetAntiHermitian();
  else if ( hermitian == +1 ) Rando.SetHermitian();
  else  Rando.SetNonHermitian();


  arma::mat symmetry_allowed = arma::zeros( arma::size( Rando.OneBody) );
  for ( auto i : modelspace.all_orbits)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    for ( auto j : Rando.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      symmetry_allowed(i,j) = 1;
    }
  }
  Rando.OneBody.randn() ;
  Rando.OneBody = Rando.OneBody % symmetry_allowed ;
  Rando.OneBody += hermitian * Rando.OneBody.t();

  if ( particle_rank > 1 )
  {
    for (auto& itmat : Rando.TwoBody.MatEl )
    {
      itmat.second.randn();
      itmat.second += hermitian * itmat.second.t();
    }
  }

 // If needed, fill the 3-body piece. Since the storage scheme has
 // places for matrix elements that are zero by symmetry, we need
 // to actually loop through and make sure we don't give a finite value
 // to something that should be zero.
//  if ( particle_rank > 2 )
//  {
//    std::default_random_engine generator(random_seed);
//    double mean = 0;
//    double stddev = 1;
//    std::normal_distribution<double> distribution(mean,stddev);
//
//    for ( auto& iter : Rando.ThreeBody.OrbitIndexHash )
//    {
//      size_t a,b,c,d,e,f;
//      size_t key = iter.first;
//      size_t lookup_abcdef = iter.second;
//      Rando.ThreeBody.KeyUnhash(key,a,b,c,d,e,f);
//
//      size_t key_check = Rando.ThreeBody.KeyHash(a,b,c,d,e,f);
//      // while we're unit testing, we may as well make sure the hash works in both directions
//      if (key_check != key )
//      {
//       std::cout << "WARNING!!!!!!!!   key_check != key : " << key_check << " != " << key << std::endl;
//      }
//
//      Orbit& oa = Rando.modelspace->GetOrbit(a);
//      Orbit& ob = Rando.modelspace->GetOrbit(b);
//      Orbit& oc = Rando.modelspace->GetOrbit(c);
//      Orbit& od = Rando.modelspace->GetOrbit(d);
//      Orbit& oe = Rando.modelspace->GetOrbit(e);
//      Orbit& of = Rando.modelspace->GetOrbit(f);
//
//      int Jab_min = std::abs(oa.j2-ob.j2)/2;
//      int Jab_max = (oa.j2+ob.j2)/2;
//      int Jde_min = std::abs(od.j2-oe.j2)/2;
//      int Jde_max = (od.j2+oe.j2)/2;
//      int J_index=0;
//
//      for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
//      {
//       for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
//       {
//         int J2_min = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jde-of.j2));
//         int J2_max = std::min( 2*Jab+oc.j2, 2*Jde+of.j2);
//         for (int J2=J2_min; J2<=J2_max; J2+=2)
//         {
////           J_index += (J2-J2_min)/2*5;
//           for (int tab=0; tab<=1; tab++)
//           {
//             for (int tde=0; tde<=1; tde++)
//             {
//               for (int twoT=1; twoT<=std::min(2*tab+1,2*tde+1); twoT+=2)
//               {
////                 int Tindex = 2*tab + tde + (twoT-1)/2;
//                 int Tindex =   (J2-J2_min)/2*5 +  2*tab + tde + (twoT-1)/2;
//                 size_t location = lookup_abcdef + J_index + Tindex;
//
//                 if ( ( a==b and (tab+Jab)%2==0 )
//                   or ( d==e and (tde+Jde)%2==0 )
//                   or ( a==b and a==c and twoT==3 and oa.j2<3 )
//                   or ( d==e and d==f and twoT==3 and od.j2<3 ))
//                 {
//                    Rando.ThreeBody.MatEl[location] = 0.0; // zero by symmetry 
////                    if (a==4 and b==4 and c==0 and d==4 and e==4 and f==0)
////                    {
////                       std::cout << "SKIPPING Jab,Jde,J2, tab,tde,twoT  " << Jab << " " << Jde << " " << J2 << "   " << tab << " " << tde << " " << twoT
////                                 << " (location = " << location  << ")   to " << Rando.ThreeBody.MatEl[location] << std::endl;
////                    }
////                       std::cout << " skip.... location = " << lookup_abcdef << " + " << J_index << " + " << Tindex << " = " << location
////                                 << " abcdef  " << a << " " << b << " " << c << " " << d << " " << e << " " << f
////                                 << "  key is " << key   << "   element 874 is " << Rando.ThreeBody.MatEl[874] << std::endl;
//                 }
//                 else
//                 {
//                    Rando.ThreeBody.MatEl[location] = distribution(generator); // make it a random number
//                       if (a==d and b==e and c==f) // if we're storing the hermitian conjugate, make sure they're compatible
//                       {
//                         if ( Jab>Jde or ( Jab==Jde and tab>tde ) )
//                         {
////                            std::cout << "match. location = " << location << "  abcdef = " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
////                            std::cout << " before ... element 874 is " << Rando.ThreeBody.MatEl[874] << std::endl;
//                            double me_flip = Rando.ThreeBody.GetME(Jde, Jab, J2, tde, tab, twoT, d,e,f,a,b,c);
//                            Rando.ThreeBody.MatEl[location] = hermitian * me_flip;
////                       std::cout << "  after ... element 874 is " << Rando.ThreeBody.MatEl[874] << std::endl;
//                         }
//                         
//                       }
////                    if (a==4 and b==4 and c==0 and d==4 and e==4 and f==0)
////                    {
////                       std::cout << "SETTING Jab,Jde,J2, tab,tde,twoT  " << Jab << " " << Jde << " " << J2 << "   " << tab << " " << tde << " " << twoT
////                                 << " (location = " << lookup_abcdef << " + " << J_index << " + " << Tindex << " = " << location  << ")   to " << Rando.ThreeBody.MatEl[location] << std::endl;
////                    }
////                       std::cout << " set ... element 874 is " << Rando.ThreeBody.MatEl[874] << std::endl;
//                 }
//               } //T2
//             } // tde
//           } // tab
//         } //J2
//       J_index += (J2_max-J2_min+2)/2*5;
//
//       } //Jde
//      } //Jab
//      
//    }
//
//
////    for ( size_t i=0; i<Rando.ThreeBody.MatEl.size(); i++)
////    {
////      Rando.ThreeBody.MatEl[i] = distribution(generator);
////    }
//
//  }


  std::cout << "In  " << __func__  << "  norm of 1b : " << Rando.OneBodyNorm();
  if (particle_rank > 1) std::cout << "   norm of 2b " << Rando.TwoBodyNorm();
  if (particle_rank > 2) std::cout << "   norm of 3b " << Rando.ThreeBodyNorm();
  std::cout << std::endl;
  
  return Rando;
}



void UnitTest::Test3BodyAntisymmetry()
{

  arma::arma_rng::set_seed( random_seed );
  Operator Y = RandomOp(*modelspace, 0, 0, 0, 3, 1);



//  int a=0,b=0,c=3,d=2,e=2,f=5;
//  int a=0,b=0,c=5,d=2,e=2,f=3;
//  int a=4,b=4,c=0,d=4,e=4,f=0;
  int a=5,b=0,c=0,d=5,e=0,f=0;
//  int a=5,b=4,c=0,d=5,e=4,f=0;

  Orbit& oa = Y.modelspace->GetOrbit(a);
  Orbit& ob = Y.modelspace->GetOrbit(b);
  Orbit& oc = Y.modelspace->GetOrbit(c);
  Orbit& od = Y.modelspace->GetOrbit(d);
  Orbit& oe = Y.modelspace->GetOrbit(e);
  Orbit& of = Y.modelspace->GetOrbit(f);

  double yJ0_abcdef = Y.ThreeBody.GetME_pn(0,0,1,a,b,c,d,e,f);
  double yJ1_abcdef = Y.ThreeBody.GetME_pn(1,0,1,a,b,c,d,e,f);
  double yJ0_cbadef = Y.ThreeBody.GetME_pn(0,0,1,c,b,a,d,e,f);
  double yJ1_cbadef = Y.ThreeBody.GetME_pn(1,0,1,c,b,a,d,e,f);

  std::cout << "yJ0_abcdef = " << yJ0_abcdef << std::endl;
  std::cout << "yJ1_abcdef = " << yJ1_abcdef << std::endl;
  std::cout << "yJ0_cbadef = " << yJ0_cbadef << std::endl;
  std::cout << "    should be - " << AngMom::phase(oa.j2+ob.j2+oc.j2) << " * 1 * 1 * " << AngMom::SixJ(0.5,0.5,0,0.5,0.5,0) << " * " << yJ0_abcdef
            << " = " << -AngMom::phase(oa.j2+ob.j2+oc.j2) * AngMom::SixJ(0.5,0.5,0,0.5,0.5,0) * yJ0_abcdef << std::endl;
  std::cout << "yJ1_cbadef = " << yJ1_cbadef << std::endl;
  std::cout << "    should be - " << AngMom::phase(oa.j2+ob.j2+oc.j2) << " * sqrt(3) * 1 * " << AngMom::SixJ(0.5,0.5,1,0.5,0.5,0) << " * " << yJ0_abcdef
            << " = " << -AngMom::phase(oa.j2+ob.j2+oc.j2) * sqrt(3) * AngMom::SixJ(0.5,0.5,1,0.5,0.5,0) * yJ0_abcdef << std::endl;

//  return;


  for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
  {
  for (int mb=-ob.j2; mb<=ob.j2; mb+=2)
  {
  for (int mc=-oc.j2; mc<=oc.j2; mc+=2)
  {
  for (int md=-od.j2; md<=od.j2; md+=2)
  {
  for (int me=-oe.j2; me<=oe.j2; me+=2)
  {
   int mf=ma+mb+mc-md-me;
   if ( std::abs(mf)>of.j2 ) continue;
  std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << " @mvals: " << ma << " " << mb << " " << mc << " " << md << " " << me << " " << mf << std::endl;
  std::cout << std::endl << "@abc" << std::endl;
  double yabcdef = GetMschemeMatrixElement_3b(Y, a,ma, b,mb ,c,mc, d,md, e,me, f,mf);
  if (std::abs(yabcdef)<1e-8) continue;
//  std::cout << "@bca" << std::endl;
//  double ybcadef = GetMschemeMatrixElement_3b(Y, b,mb, c,mc, a,ma, d,md, e,me, f,mf);
//  std::cout << "@cab" << std::endl;
//  double ycabdef = GetMschemeMatrixElement_3b(Y, c,mc, a,ma, b,mb, d,md, e,me, f,mf);
//
//  std::cout << "@acb" << std::endl;
//  double yacbdef = GetMschemeMatrixElement_3b(Y, a,ma, c,mc, b,mb, d,md, e,me, f,mf);
//  std::cout << "@bac" << std::endl;
//  double ybacdef = GetMschemeMatrixElement_3b(Y, b,mb, a,ma, c,mc, d,md, e,me, f,mf);
  std::cout << std::endl << "@cba" << std::endl;
  double ycbadef = GetMschemeMatrixElement_3b(Y, c,mc, b,mb, a,ma, d,md, e,me, f,mf);

  std::cout << " @@mvals: " << ma << " " << mb << " " << mc << " " << md << " " << me << " " << mf << std::endl;
//  std::cout << "ANTISYMMETRY:  abc: " << yabcdef << " bca: " << ybcadef << " cab: " << ycabdef
//            << "    acb: " << yacbdef << " bac: " << ybacdef << " cba: " << ycbadef << std::endl;
  std::cout << "ANTISYMMETRY:  abc: " << yabcdef  << " cba: " << ycbadef << std::endl;
  return;

  }
  }
  }
  }
  }


//  for (int a : Y.modelspace->all_orbits )
//  {
//    for (int b : Y.modelspace->all_orbits )
//    {
//      for (int c : Y.modelspace->all_orbits )
//      {
//      }
//    }
//  }


}



void UnitTest::TestCommutators()
{
  arma::arma_rng::set_seed( random_seed );
  Operator X = RandomOp(*modelspace, 0, 0, 0, 2, -1);
  Operator Y = RandomOp(*modelspace, 0, 0, 0, 2, +1);

  bool all_good = true;

  all_good &= Test_comm110ss( X, Y );
  all_good &= Test_comm220ss( X, Y );

  all_good &= Test_comm111ss( X, Y );
  all_good &= Test_comm121ss( X, Y );
  all_good &= Test_comm221ss( X, Y );

  all_good &= Test_comm122ss( X, Y );
  all_good &= Test_comm222_pp_hhss(  X, Y );
  all_good &= Test_comm222_phss(  X, Y );

  if ( all_good )
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }

}


void UnitTest::TestCommutators3()
{
  std::cout << " random_seed = " << random_seed << std::endl;
  arma::arma_rng::set_seed( random_seed );
  modelspace->PreCalculateSixJ();
  Operator X = RandomOp(*modelspace, 0, 0, 0, 3, -1);
  random_seed ++;
  Operator Y = RandomOp(*modelspace, 0, 0, 0, 3, +1);
  random_seed --;

  bool all_good = true;

  all_good &= Test_comm330ss( X, Y );

  if ( all_good )
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }
}



double UnitTest::GetMschemeMatrixElement_1b( const Operator& Op, int a, int ma, int b, int mb )
{
  double matel = 0;
  int Jop = Op.GetJRank();
  if ( Jop == 0) // scalar operator
  {
    if ( ma==mb )
    {
      matel = Op.OneBody(a,b);
    }
  }
  else // reduced matrix element using Wigner-Eckart thm:   <a ma | Op JM | b mb > = (-1)^2J * <J M b mb | a ma > / sqrt(2ja+1) < a || Op J || b >
  {
    Orbit& oa = Op.modelspace->GetOrbit(a);
    Orbit& ob = Op.modelspace->GetOrbit(b);
    matel = AngMom::phase(2*Jop) * AngMom::CG(Jop, 0.5*(ma-mb), 0.5*ob.j2, 0.5*mb, 0.5*oa.j2, 0.5*ma) / sqrt(oa.j2+1.0) * Op.OneBody(a,b);
  }

  return matel;
}



double UnitTest::GetMschemeMatrixElement_2b( const Operator& Op, int a, int ma, int b, int mb, int c, int mc, int d, int md )
{
  double matel = 0;
  int Jop = Op.GetJRank();
  Orbit& oa = Op.modelspace->GetOrbit(a);
  Orbit& ob = Op.modelspace->GetOrbit(b);
  Orbit& oc = Op.modelspace->GetOrbit(c);
  Orbit& od = Op.modelspace->GetOrbit(d);
  if (a==b and ma==mb) return 0;
  if (c==d and mc==md) return 0;
  if ( Jop == 0) // scalar operator
  {
    if ( (ma+mb) == (mc+md) )
    {
      int Jmin = std::max(  std::abs( oa.j2-ob.j2) , std::abs(oc.j2-od.j2)  ) / 2;
      int Jmax = std::min( oa.j2+ob.j2  ,  oc.j2+od.j2 ) / 2;
      int M = (ma+mb)/2;
      for (int J=Jmin; J<=Jmax; J++)
      {
        // We take the un-normalized TBME (the one with the tilde) so we don't
        // need to worry about normalization factors.
        double clebsch_ab = AngMom::CG( 0.5*oa.j2, 0.5*ma, 0.5*ob.j2, 0.5*mb,  J, M );
        double clebsch_cd = AngMom::CG( 0.5*oc.j2, 0.5*mc, 0.5*od.j2, 0.5*md,  J, M );
        matel += clebsch_ab * clebsch_cd * Op.TwoBody.GetTBME_J(J,a,b,c,d);
      }
    }
  }
  else
  {
    std::cout << " WARNING!!! " << __func__ << "   not yet implemented for tensor operator " << std::endl;
    exit(0);
  }

  return matel;
}


double UnitTest::GetMschemeMatrixElement_3b( const Operator& Op, int a, int ma, int b, int mb, int c, int mc, int d, int md, int e, int me, int f, int mf )
{
  double matel = 0;
  int Jop = Op.GetJRank();
  Orbit& oa = Op.modelspace->GetOrbit(a);
  Orbit& ob = Op.modelspace->GetOrbit(b);
  Orbit& oc = Op.modelspace->GetOrbit(c);
  Orbit& od = Op.modelspace->GetOrbit(d);
  Orbit& oe = Op.modelspace->GetOrbit(e);
  Orbit& of = Op.modelspace->GetOrbit(f);
  if (a==b and ma==mb) return 0;
  if (a==c and ma==mc) return 0;
  if (b==c and mb==mc) return 0;
  if (d==e and md==me) return 0;
  if (d==f and md==mf) return 0;
  if (e==f and me==mf) return 0;

//  if (a==b and a==c and oa.j2<2) return 0;
//  if (d==e and d==f and od.j2<2) return 0;

  if ( Jop != 0)  return 0;// Haven't even thought about tensor 3N stuff yet...

  if ( (ma+mb+mc) != (md+me+mf) ) return 0;

  int Mab = (ma+mb)/2;
  int Mde = (md+me)/2;
  int twoM = ma+mb+mc;
  int Jab_min = std::max( std::abs(Mab), std::abs( oa.j2-ob.j2) / 2);
  int Jab_max =  (oa.j2+ob.j2 ) / 2;
  int Jde_min = std::max( std::abs(Mde), std::abs( od.j2-oe.j2) / 2);
  int Jde_max =  (od.j2+oe.j2 ) / 2;

  for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
  {
    double clebsch_ab = AngMom::CG(0.5*oa.j2, 0.5*ma,  0.5*ob.j2, 0.5*mb,  Jab, Mab );
    if ( std::abs(clebsch_ab)<1e-6 ) continue;
    for (int Jde=Jde_min; Jde<=Jde_max; Jde++)
    {
      double clebsch_de = AngMom::CG(0.5*od.j2, 0.5*md,  0.5*oe.j2, 0.5*me,  Jde, Mde );
      if ( std::abs(clebsch_de)<1e-6 ) continue;
      int twoJ_min = std::max( {std::abs(twoM), std::abs(2*Jab-oc.j2),  std::abs(2*Jde-of.j2) } );
      int twoJ_max = std::min( 2*Jab+oc.j2,  2*Jde+of.j2 );
      for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
      {
        double clebsch_abc = AngMom::CG( Jab, Mab, 0.5*oc.j2, 0.5*mc,  0.5*twoJ, 0.5*twoM );
        double clebsch_def = AngMom::CG( Jde, Mde, 0.5*of.j2, 0.5*mf,  0.5*twoJ, 0.5*twoM );
        if (std::abs(clebsch_abc)<1e-6 or std::abs(clebsch_def)<1e-6) continue; // avoid the look up, with possible 6js...

        double meJ = Op.ThreeBody.GetME_pn(Jab, Jde, twoJ, a,b,c,d,e,f);
//        matel += clebsch_ab * clebsch_abc * clebsch_de * clebsch_def * Op.ThreeBody.GetME_pn(Jab, Jde, twoJ, a,b,c,d,e,f);
        matel += clebsch_ab * clebsch_abc * clebsch_de * clebsch_def * meJ;

//        if ((a==0 and b==0 and c==3 and d==2 and e==5 and f==2) or (d==0 and e==0 and f==3 and a==2 and b==5 and c==2)) 
//        if ((a==0 and b==0 and c==3 and d==2 and e==2 and f==5 and Jde==1) ) 
//        if ((a==0 and b==0 and c==3 and d==2 and e==2 and f==5 ) or  (a==0 and b==3 and c==0 and d==2 and e==2 and f==5 ) or (a==3 and b==0 and c==0 and d==2 and e==2 and f==5 )) 
//        if ((a==0 and b==0 and c==5 and d==0 and e==0 and f==5 ) or  (a==0 and b==5 and c==0 and d==0 and e==0 and f==5 ) or (a==5 and b==0 and c==0 and d==0 and e==0 and f==5 )) 
//        {
//        std::cout << "$abc: " << a << " " << b << " " << c << " def: " << d << " " << e << " " << f << std::endl;
//        std::cout << "$m vals: " << ma << " " << mb << " " << mc << "  " << md << " " << me << " " << mf << std::endl;
//        std::cout << "        Jab Jde twoJ " << Jab << " " << Jde << " " << twoJ
//                  << " clebsch: " << clebsch_ab << " " << clebsch_de << " " << clebsch_abc << " " << clebsch_def
//                  << "   matel_J " << meJ 
//                  << "  matel = " << matel << std::endl;
//        }
//        if (a==0 and b==0 and c==3 and d==2 and e==5 and f==2) std::cout << "    003252: Jab,Jde,J " << Jab << " " << Jde << " " << twoJ << "   -> " << Op.ThreeBody.GetME_pn(Jab, Jde, twoJ, a,b,c,d,e,f) << " -> " << matel << std::endl;
//        if (a==0 and b==0 and c==3 and d==2 and e==2 and f==5) std::cout << "    003225: Jab,Jde,J " << Jab << " " << Jde << " " << twoJ << "   -> " << Op.ThreeBody.GetME_pn(Jab, Jde, twoJ, a,b,c,d,e,f) << " -> " << matel << std::endl;
//        if (a==2 and b==2 and c==5 and d==0 and e==0 and f==3) std::cout << "   225003: Jab,Jde,J " << Jab << " " << Jde << " " << twoJ << "   -> " << Op.ThreeBody.GetME_pn(Jab, Jde, twoJ, a,b,c,d,e,f) << " -> " << matel << std::endl;
      }
    }
  }
//  if (a==0 and b==0 and c==3 and d==2 and e==5 and f==2) std::cout << "  ---done--> " << matel << std::endl;
//  if (a==0 and b==0 and c==3 and d==2 and e==2 and f==5) std::cout << "  ---done--> " << matel << std::endl;
//  if (a==2 and b==2 and c==5 and d==0 and e==0 and f==3) std::cout << "  ---done--> " << matel << std::endl;
  return matel;

}



/// M-Scheme Formula:
///
/// Z0 = 1/4 * sum_ab na(1-nb) (Xab * Yba - Yab * Xba)
///
bool UnitTest::Test_comm110ss( const Operator& X, const Operator& Y )
{
   
  Operator Z_J( Y );
  Z_J.Erase();

  Commutator::comm110ss( X, Y, Z_J);

  double Z0_m = 0;
  for (auto a : X.modelspace->all_orbits )
  {
   Orbit& oa = X.modelspace->GetOrbit(a);
   for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
   {
     double na = oa.occ;
     for (auto b : X.modelspace->all_orbits )
     {
      Orbit& ob = X.modelspace->GetOrbit(b);
      double nb = ob.occ;
      for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
      {
        double Xab = GetMschemeMatrixElement_1b( X, a, ma, b, mb );
        double Xba = GetMschemeMatrixElement_1b( X, b, mb, a, ma );
        double Yba = GetMschemeMatrixElement_1b( Y, b, mb, a, ma );
        double Yab = GetMschemeMatrixElement_1b( Y, a, ma, b, mb );

        Z0_m += na * (1-nb) * ( Xab*Yba - Yab*Xba );
      }
     }
   }
  }

  double summed_error =  Z0_m - Z_J.ZeroBody;
  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << Z0_m << " " << Z_J.ZeroBody
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}


/// M-Scheme Formula:
///
/// Z0 = 1/4 sum_abcd nanb (1-nc)(1-nd) ( Xabcd * Ycdab - Yabcd * Xcdab )
///
bool UnitTest::Test_comm220ss( const Operator& X, const Operator& Y) 
{
  Operator Z_J( Y );
  Z_J.Erase();

  Commutator::comm220ss( X, Y, Z_J);

  double Z0_m = 0;
  for (auto a : X.modelspace->all_orbits )
  {
   Orbit& oa = X.modelspace->GetOrbit(a);
   double na = oa.occ;
   for (auto b : X.modelspace->all_orbits )
   {
    Orbit& ob = X.modelspace->GetOrbit(b);
    double nb = ob.occ;
    for (auto c : X.modelspace->all_orbits )
    {
     Orbit& oc = X.modelspace->GetOrbit(c);
     double nc = oc.occ;
     for (auto d : X.modelspace->all_orbits )
     {
      Orbit& od = X.modelspace->GetOrbit(d);
      if ( (oa.l+ob.l+oc.l+od.l)%2 >0) continue;
      if ( (oa.tz2+ob.tz2) != (oc.tz2+od.tz2) ) continue;
      double nd = od.occ;


      for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
      {
        for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
        {
          for (int mc=-oc.j2; mc<=oc.j2; mc+=2 )
          {
            int md = ma+mb-mc;

            double Xabcd = GetMschemeMatrixElement_2b( X, a,ma, b,mb, c,mc, d,md );
            double Yabcd = GetMschemeMatrixElement_2b( Y, a,ma, b,mb, c,mc, d,md );
            double Xcdab = GetMschemeMatrixElement_2b( X, c,mc, d,md, a,ma, b,mb );
            double Ycdab = GetMschemeMatrixElement_2b( Y, c,mc, d,md, a,ma, b,mb );

            Z0_m += (1./4) * na * nb * (1.0-nc) * (1.0-nd) * ( Xabcd * Ycdab - Yabcd * Xcdab );

          } // for mc
        } // for mb
      } // for ma
     } // for d
    } // for c
   } // for b
  } // for a

  double summed_error =  Z0_m - Z_J.ZeroBody;
  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << Z0_m << " " << Z_J.ZeroBody
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}




/// M-Scheme Formula:
///
///  Zij = sum_a (Xia * Yaj - Yia * Xaj)
///
bool UnitTest::Test_comm111ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();


  Commutator::comm111ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      Orbit& oj = X.modelspace->GetOrbit(j);
      int mj = oj.j2;
      double Zm_ij = 0;
      if ( (oi.j2==oj.j2) and (oi.l==oj.l) and (oi.tz2==oj.tz2) )
      {
        for (auto a : X.modelspace->all_orbits)
        {
         Orbit& oa = X.modelspace->GetOrbit(a);
         for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
         {
           double Xia = GetMschemeMatrixElement_1b( X, i, mi, a, ma );
           double Yia = GetMschemeMatrixElement_1b( Y, i, mi, a, ma );
           double Xaj = GetMschemeMatrixElement_1b( X, a, ma, j, mj );
           double Yaj = GetMschemeMatrixElement_1b( Y, a, ma, j, mj );

           Zm_ij += Xia * Yaj - Yia * Xaj;

         }
        }
        double ZJ_ij = Z_J.OneBody(i,j);
        double err = Zm_ij - ZJ_ij;
        if (std::abs(err)>1e-6)
        {
          std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij
                    << "   ZJ_ij = " << Z_J.OneBody(i,j) << "   err = " << err << std::endl; 
        }
        summed_error += err*err;
        sum_m += Zm_ij*Zm_ij;
        sum_J += ZJ_ij*ZJ_ij;
      } 
    }
  }

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.OneBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 1b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}



/// M-Scheme Formula:
/// Zij = [X1b,Y2b] - [Y1b,X2b]
//
/// Zij = sum_ab na*(1-nb) ( Xab * Ybiaj - Yaibj * Xba
///                        - Yab * Xbiaj + Xaibj * Yba )
///
bool UnitTest::Test_comm121ss( const Operator& X, const Operator& Y)
{

  Operator Z_J( Y );
  Z_J.Erase();


  Commutator::comm121ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      Orbit& oj = X.modelspace->GetOrbit(j);
      int mj = oj.j2;
      double Zm_ij = 0;
      if ( (oi.j2==oj.j2) and (oi.l==oj.l) and (oi.tz2==oj.tz2) )
      {
        for (auto a : X.modelspace->all_orbits)
        {
          Orbit& oa = X.modelspace->GetOrbit(a);
          double na = oa.occ;
          for (auto b : X.modelspace->all_orbits)
          {
            Orbit& ob = X.modelspace->GetOrbit(b);
            double nb = ob.occ;
            for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
            {
             for (int mb=-ob.j2; mb<=ob.j2; mb+=2)
             {
               double Xab = GetMschemeMatrixElement_1b( X, a, ma, b, mb );
               double Yab = GetMschemeMatrixElement_1b( Y, a, ma, b, mb );
               double Xba = GetMschemeMatrixElement_1b( X, b, mb, a, ma );
               double Yba = GetMschemeMatrixElement_1b( Y, b, mb, a, ma );
               double Xbiaj = GetMschemeMatrixElement_2b( X, b,mb, i,mi, a,ma, j,mj );
               double Ybiaj = GetMschemeMatrixElement_2b( Y, b,mb, i,mi, a,ma, j,mj );
               double Xaibj = GetMschemeMatrixElement_2b( X, a,ma, i,mi, b,mb, j,mj );
               double Yaibj = GetMschemeMatrixElement_2b( Y, a,ma, i,mi, b,mb, j,mj );

               Zm_ij += na * (1-nb) * ( Xab * Ybiaj - Yaibj * Xba
                                      - Yab * Xbiaj + Xaibj * Yba );
              
             }// for mb
            }// for ma
          }// for b 
        }// for a
      } //if ji=jj etc. 
    double ZJ_ij = Z_J.OneBody(i,j);
    double err = Zm_ij - ZJ_ij;
    if (std::abs(err)>1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij
                << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl; 
    }
    summed_error += err*err;
    sum_m += Zm_ij*Zm_ij;
    sum_J += ZJ_ij*ZJ_ij;
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.OneBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 1b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}




/// M-Scheme Formula:
//
// Zij = 1/2 sum_abc [ na*nb*(1-nc) +(1-na)*(1-nb)*nc ] * (Xciab * Yabcj - Yciab * Xabcj)
//
// THIS ONE DOESN'T WORK YET!!!!
//
bool UnitTest::Test_comm221ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();


//  Commutator::comm222_pp_hh_221ss( X, Y, Z_J ) ; 
  Commutator::comm221ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      Orbit& oj = X.modelspace->GetOrbit(j);
      int mj = oj.j2;
      double Zm_ij = 0;
      if ( (oi.j2==oj.j2) and (oi.l==oj.l) and (oi.tz2==oj.tz2) )
      {
        for (auto a : X.modelspace->all_orbits )
        {
          Orbit& oa = X.modelspace->GetOrbit(a);
          double na = oa.occ;
          for (auto b : X.modelspace->all_orbits )
          {

            Orbit& ob = X.modelspace->GetOrbit(b);
            double nb = ob.occ;
            for (auto c : X.modelspace->all_orbits )
            {
               Orbit& oc = X.modelspace->GetOrbit(c);
               double nc = oc.occ;
               if ( std::abs(na*nb*(1-nc)+(1-na)*(1-nb)*nc )<1e-6) continue;
               if ( (oi.l+oc.l + oa.l+ob.l)%2>0) continue;
               if ( (oi.tz2+oc.tz2) != (oa.tz2+ob.tz2) ) continue;

               for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
               {
                 for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
                 {

                   int mc = ma+mb-mi;
                   if ( std::abs(mc)>oc.j2 ) continue;
                   double Xciab = GetMschemeMatrixElement_2b( X, c,mc, i,mi, a,ma, b,mb );
                   double Yciab = GetMschemeMatrixElement_2b( Y, c,mc, i,mi, a,ma, b,mb );
                   double Xabcj = GetMschemeMatrixElement_2b( X, a,ma, b,mb, c,mc, j,mj );
                   double Yabcj = GetMschemeMatrixElement_2b( Y, a,ma, b,mb, c,mc, j,mj );

                   Zm_ij += (1./2) * ( na*nb*(1-nc)+(1-na)*(1-nb)*nc ) * (Xciab * Yabcj - Yciab * Xabcj );

                 }// for ma
               }// for mb
            }// for c
          }// for b
        }// for a

      } //if ji=jj etc. 

      double ZJ_ij = Z_J.OneBody(i,j);
      double err = Zm_ij - ZJ_ij;
      if (std::abs(err)>1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij
                  << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl; 
      }
      summed_error += err*err;
      sum_m += Zm_ij*Zm_ij;
      sum_J += ZJ_ij*ZJ_ij;
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.OneBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 1b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}



/// M-Scheme Formula:
//
// Zijkl = sum_a   (Xia * Yajkl + Xja * Yiakl - Yijal * Xak - Yijka * Xal)
//               - (Yia * Xajkl + Yja * Xiakl - Xijal * Yak - Xijka * Yal)
//
bool UnitTest::Test_comm122ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();


  Commutator::comm122ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      if (j<i) continue;
      Orbit& oj = X.modelspace->GetOrbit(j);
      if ( not  (oi.l==oj.l and oi.j2==oj.j2 and oi.tz2==oj.tz2)  ) continue;
      for (auto k : X.modelspace->all_orbits )
      {
        Orbit& ok = X.modelspace->GetOrbit(k);
        if (k<i) continue;
        for (auto l : X.modelspace->all_orbits )
        {
          if (l<k) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          if ( not  (ok.l==ol.l and ok.j2==ol.j2 and ok.tz2==ol.tz2) ) continue;
          if ( (oi.l+oj.l+ok.l+ol.l)%2>0 ) continue;
          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
            for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
            {

              int ml = mi + mj - mk;
              if ( std::abs(ml)>ol.j2 ) continue;

              double Zm_ijkl = 0;

              for (auto a : X.modelspace->all_orbits )
              {
               Orbit& oa = X.modelspace->GetOrbit(a);

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
                   double Xia = GetMschemeMatrixElement_1b( X, i, mi, a, ma );
                   double Xja = GetMschemeMatrixElement_1b( X, j, mj, a, ma );
                   double Xak = GetMschemeMatrixElement_1b( X, a, ma, k, mk );
                   double Xal = GetMschemeMatrixElement_1b( X, a, ma, l, ml );
                   double Yia = GetMschemeMatrixElement_1b( Y, i, mi, a, ma );
                   double Yja = GetMschemeMatrixElement_1b( Y, j, mj, a, ma );
                   double Yak = GetMschemeMatrixElement_1b( Y, a, ma, k, mk );
                   double Yal = GetMschemeMatrixElement_1b( Y, a, ma, l, ml );
                   double Xajkl = GetMschemeMatrixElement_2b( X, a,ma, j,mj, k,mk, l,ml );
                   double Xiakl = GetMschemeMatrixElement_2b( X, i,mi, a,ma, k,mk, l,ml );
                   double Xijal = GetMschemeMatrixElement_2b( X, i,mi, j,mj, a,ma, l,ml );
                   double Xijka = GetMschemeMatrixElement_2b( X, i,mi, j,mj, k,mk, a,ma );
                   double Yajkl = GetMschemeMatrixElement_2b( Y, a,ma, j,mj, k,mk, l,ml );
                   double Yiakl = GetMschemeMatrixElement_2b( Y, i,mi, a,ma, k,mk, l,ml );
                   double Yijal = GetMschemeMatrixElement_2b( Y, i,mi, j,mj, a,ma, l,ml );
                   double Yijka = GetMschemeMatrixElement_2b( Y, i,mi, j,mj, k,mk, a,ma );

                   Zm_ijkl +=  ( Xia*Yajkl + Xja*Yiakl - Yijal*Xak - Yijka*Xal)
                              -( Yia*Xajkl + Yja*Xiakl - Xijal*Yak - Xijka*Yal) ;
                 
                }// for ma
              }// for a

              double ZJ_ijkl = GetMschemeMatrixElement_2b( Z_J, i,mi, j,mj, k,mk, l,ml ) ;
              double err = Zm_ijkl - ZJ_ijkl;
              if (std::abs(err)>1e-6)
              {
                std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                          << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl; 
              }
              summed_error += err*err;
              sum_m += Zm_ijkl*Zm_ijkl;
              sum_J += ZJ_ijkl*ZJ_ijkl;

            }// for mk
          }// for mj
        }// for l
      }// for k
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}




/// M-Scheme Formula:
//
// Zijkl = 1/2 sum_ab [(1-na)*(1-nb) - na*nb ] * ( Xijab * Yabkl - Yijab * Xabkl )
//
bool UnitTest::Test_comm222_pp_hhss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();


  Commutator::comm222_pp_hhss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      if (j<i) continue;
      Orbit& oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits )
      {
        Orbit& ok = X.modelspace->GetOrbit(k);
        if (k<i) continue;
        for (auto l : X.modelspace->all_orbits )
        {
          if (l<k) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {
             if ( std::abs(mk)>ol.j2 ) continue;
             int ml = mi + mj - mk;

             double Zm_ijkl = 0;
             for (auto a : X.modelspace->all_orbits )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
              for (auto b : X.modelspace->all_orbits )
              {
                Orbit& ob = X.modelspace->GetOrbit(b);
                double nb = ob.occ;
                if ( std::abs( ( (1-na)*(1-nb) - na*nb ) ) < 1e-6) continue;

                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
                   int mb = mi + mj - ma;
                   if (std::abs(mb)>ob.j2) continue;
                
                  double Xijab = GetMschemeMatrixElement_2b( X, i,mi, j,mj, a,ma, b,mb );
                  double Xabkl = GetMschemeMatrixElement_2b( X, a,ma, b,mb, k,mk, l,ml );
                  double Yijab = GetMschemeMatrixElement_2b( Y, i,mi, j,mj, a,ma, b,mb );
                  double Yabkl = GetMschemeMatrixElement_2b( Y, a,ma, b,mb, k,mk, l,ml );
                  Zm_ijkl += (1./2) * ( (1-na)*(1-nb) - na*nb ) * (Xijab * Yabkl - Yijab * Xabkl );
                }// for ma
              }// for b
             }// for a


             double ZJ_ijkl = GetMschemeMatrixElement_2b( Z_J, i,mi, j,mj, k,mk, l,ml ) ;
             double err = Zm_ijkl - ZJ_ijkl;
             if (std::abs(err)>1e-6)
             {
               std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                         << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl; 
             }
             summed_error += err*err;
             sum_m += Zm_ijkl*Zm_ijkl;
             sum_J += ZJ_ijkl*ZJ_ijkl;

           }// for mk
          }// for mj
        }// for l
      }// for k
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}



/// M-Scheme Formula:
//
// Zijkl = sum_ab [na*(1-nb) - nb*(1-na)] * (1-Pij)(1-Pkl) ( Xaibk Ybjal  )
//
bool UnitTest::Test_comm222_phss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();

  Commutator::comm222_phss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      if (j<i) continue;
      Orbit& oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits )
      {
        Orbit& ok = X.modelspace->GetOrbit(k);
        if (k<i) continue;
        for (auto l : X.modelspace->all_orbits )
        {
          if (l<k) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;


         for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
         {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {

             int ml = mi + mj - mk;
             if ( std::abs(ml)>ol.j2 ) continue;

             double Zm_ijkl = 0;
             for (auto a : X.modelspace->all_orbits )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
              for (auto b : X.modelspace->all_orbits )
              {
                Orbit& ob = X.modelspace->GetOrbit(b);
                double nb = ob.occ;
                if ( std::abs(  na*(1-nb) - nb*(1-na) ) < 1e-6) continue;

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
                 for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
                 {

                   double Xaibk = GetMschemeMatrixElement_2b( X, a,ma, i,mi, b,mb, k,mk );
                   double Ybjal = GetMschemeMatrixElement_2b( Y, b,mb, j,mj, a,ma, l,ml );
    
                   // Pij
                   double Xajbk = GetMschemeMatrixElement_2b( X, a,ma, j,mj, b,mb, k,mk );
                   double Ybial = GetMschemeMatrixElement_2b( Y, b,mb, i,mi, a,ma, l,ml );
    
                   // Pkl
                   double Xaibl = GetMschemeMatrixElement_2b( X, a,ma, i,mi, b,mb, l,ml );
                   double Ybjak = GetMschemeMatrixElement_2b( Y, b,mb, j,mj, a,ma, k,mk );
    
                   // PijPkl
                   double Xajbl = GetMschemeMatrixElement_2b( X, a,ma, j,mj, b,mb, l,ml );
                   double Ybiak = GetMschemeMatrixElement_2b( Y, b,mb, i,mi, a,ma, k,mk );
    
                   Zm_ijkl += ( na*(1-nb) - nb*(1-na) ) * ( Xaibk*Ybjal - Xajbk*Ybial - Xaibl*Ybjak + Xajbl*Ybiak );


                 }// for mb
                }// for ma
              }// for b
             }// for a

             double ZJ_ijkl = GetMschemeMatrixElement_2b( Z_J, i,mi, j,mj, k,mk, l,ml ) ;
             double err = Zm_ijkl - ZJ_ijkl;
             if (std::abs(err)>1e-6)
             {
               std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                         << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl; 
             }
             summed_error += err*err;
             sum_m += Zm_ijkl*Zm_ijkl;
             sum_J += ZJ_ijkl*ZJ_ijkl;

           }// for mk
         }// for mj
        }// for l
      }// for k
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}








/// M-Scheme Formula:
//
// Z0b = 1/36 sum_abcdef  na*nb*nc*(1-nd)(1-ne)(1-nf) ( X_abcdef * Y_defabc - Yabcdef * Xdefabc )
//
bool UnitTest::Test_comm330ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();

  Commutator::comm330ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double Z0_m = 0;
//  int norbits = X.modelspace->GetNumberOrbits();
//  #pragma omp parallel for schedule(dynamic,1) reduction(+:Z0_m)
//  for (int a=0; a<norbits; a++ )
  for (auto a : X.modelspace->all_orbits )
  {
   Orbit& oa = X.modelspace->GetOrbit(a);
   double na = oa.occ;
   for (auto b : X.modelspace->all_orbits )
   {
    Orbit& ob = X.modelspace->GetOrbit(b);
    double nb = ob.occ;
    for (auto c : X.modelspace->all_orbits )
    {
     Orbit& oc = X.modelspace->GetOrbit(c);
     double nc = oc.occ;
     if ( std::abs(na*nb*nc)<1e-6 ) continue;

     for (auto d : X.modelspace->all_orbits )
     {
      Orbit& od = X.modelspace->GetOrbit(d);
      double nd = od.occ;

      for (auto e : X.modelspace->all_orbits )
      {
       Orbit& oe = X.modelspace->GetOrbit(e);
       double ne = oe.occ;
       for (auto f : X.modelspace->all_orbits )
       {
        Orbit& of = X.modelspace->GetOrbit(f);
        double nf = of.occ;
        if ( std::abs((1.-nd)*(1.-ne)*(1.-nf))<1e-6 ) continue;

        // check parity and isospin projection
        if ( (oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0 ) continue;
        if ( (oa.tz2+ob.tz2+oc.tz2) !=(od.tz2+oe.tz2+of.tz2) ) continue;
         double dz=0;

//        if (a==0 and b==0 and c==3 and d==2 and e==5 and f==2)
        if (a==0 and b==0 and c==3 and d==2 and e==2 and f==5)
        {
         int Jab = 0;
         int twoJ = 3;
         for (int Jde=1; Jde<=2; Jde++)
         {
         double xabcdef = 0;
         for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
         {
          for (int mb=-ob.j2; mb<=ob.j2; mb+=2)
          {
           if (a==b and ma==mb) continue;
           for (int mc=-oc.j2; mc<=oc.j2; mc+=2)
           {
             if (a==c and ma==mc) continue;
             if (b==c and mb==mc) continue;
             for (int md=-od.j2; md<=od.j2; md+=2)
             {
              for (int me=-oe.j2; me<=oe.j2; me+=2)
              {
                int mf = ma+mb+mc - me-md;
                if ( std::abs(mf)>of.j2) continue;
//                if (d==e and md==me) continue;
//                if (d==f and md==mf) continue;
//                if (e==f and me==mf) continue;
                int Mab = (ma+mb)/2;
                int Mde = (md+me)/2;
                int twoM = ma+mb+mc;
                if (twoM != twoJ) continue;
                double cg1 = AngMom::CG(0.5*oa.j2,0.5*ma, 0.5*ob.j2,0.5*mb,   Jab, Mab);
                double cg2 = AngMom::CG(Jab,Mab,          0.5*oc.j2,0.5*mc,   0.5*twoJ, 0.5*twoM);
                double cg3 = AngMom::CG(0.5*od.j2,0.5*md, 0.5*oe.j2,0.5*me,   Jde, Mde);
                double cg4 = AngMom::CG(Jde,Mde ,         0.5*of.j2,0.5*mf,   0.5*twoJ, 0.5*twoM);
                xabcdef += cg1*cg2*cg3*cg4* GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, d,md, e,me, f,mf );
                double xabcfed = GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, f,mf, e,me, d,md );
                double xabcdfe = GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, d,md, f,mf, e,me );
                double xabcedf = GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, e,me, d,md, f,mf );
                double xabcefd = GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, e,me, f,mf, d,md );
                double xabcfde = GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, f,mf, d,md, e,me );
                double x_plain = GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, d,md, e,me, f,mf );
                std::cout << "m values: " << ma << " " << mb << " " << mc << "   " << md << " " <<me << " " << mf << std::endl;
                std::cout << "ANTISYMMETRY CHECK:  def :" << x_plain << "  efd: " << xabcefd  << "  fde: " << xabcfde
                          << "  fed: " << xabcfed << "  dfe: " << xabcdfe << "  edf: " << xabcedf  << std::endl;
              }
              }}}}// for ma ... md
           std::cout << "  Jde = " << Jde << "  xabcdef = " << xabcdef << std::endl;
           std::cout << "  should be " << X.ThreeBody.GetME_pn(Jab,Jde,twoJ,a,b,c,d,e,f) << std::endl;
            }

        for (int Jde=1; Jde<=2; Jde++)
        {
        double xdede = 0;
        for (int ma=-od.j2; ma<=od.j2; ma+=2)
         {
          for (int mb=-oe.j2; mb<=oe.j2; mb+=2)
          {
             for (int md=-od.j2; md<=od.j2; md+=2)
             {
                int me = ma+mb - md;
                int Mab = (ma+mb)/2;
                if (Mab != Jde) continue;
                double cg1 = AngMom::CG(0.5*od.j2,0.5*ma, 0.5*oe.j2,0.5*mb, Jde, Mab);
                double cg2 = AngMom::CG(0.5*od.j2,0.5*md, 0.5*oe.j2,0.5*me, Jde, Mab);
                xdede += cg1 * cg2 * GetMschemeMatrixElement_2b( Y, d,ma, e,mb, d,md, e,me);
             }
           }
          }
          std::cout << "trying 2b.  Jde = " << Jde << "  xdede = " << xdede << std::endl << "   should be " << Y.TwoBody.GetTBME_J(Jde,d,e,d,e) << std::endl;
          }


         }


         for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
         {
          for (int mb=-ob.j2; mb<=ob.j2; mb+=2)
          {
           if ((a==b) and (ma==mb)) continue;
           for (int mc=-oc.j2; mc<=oc.j2; mc+=2)
           {
           if ((a==c) and (ma==mc)) continue;
           if ((b==c) and (mb==mc)) continue;

//             if ( (ma+mb+mc) != 3) continue; // <--- This is just for testing REMOVE!!!

             for (int md=-od.j2; md<=od.j2; md+=2)
             {
              for (int me=-oe.j2; me<=oe.j2; me+=2)
              {
                if ((d==e) and (md==me)) continue;
                int mf = ma+mb+mc - me-md;
                if ( std::abs(mf)>of.j2 ) continue;
                if ((d==f) and (md==mf)) continue;
                if ((e==f) and (me==mf)) continue;


                double Xabcdef = GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, d,md, e,me, f,mf );
                double Yabcdef = GetMschemeMatrixElement_3b( Y, a,ma, b,mb, c,mc, d,md, e,me, f,mf );
                double Xdefabc = GetMschemeMatrixElement_3b( X, d,md, e,me, f,mf, a,ma, b,mb, c,mc );
                double Ydefabc = GetMschemeMatrixElement_3b( Y, d,md, e,me, f,mf, a,ma, b,mb, c,mc );

//                dz += (1./36) * na*nb*nc*(1-nd)*(1-ne)*(1-nf) * (Xabcdef * Ydefabc - Yabcdef * Xdefabc );
                dz += (2./36) * na*nb*nc*(1-nd)*(1-ne)*(1-nf) * (Xabcdef * Ydefabc  );
//                Z0_m += (1./36) * na*nb*nc*(1-nd)*(1-ne)*(1-nf) * (Xabcdef * Ydefabc - Yabcdef * Xdefabc );
//                std::cout << "abcdef" << a << " " << b << " " << c << " " << d << " " << e << " " << f
//                          << " Xabcdef Xdefabc  Yabcdef  Ydefabc   " << Xabcdef << " " << Xdefabc << " " << Yabcdef << " " << Ydefabc
//                          << "   Z0_m = " << Z0_m << std::endl;

//            if (a==0 and b==0 and c==3 and d==2 and e==2 and f==5)
            if (a==0 and b==0 and c==3 and d==2 and e==5 and f==2)
            {
//              std::cout << "     ma mb mc md me mf = " << ma << " " << mb << " " << mc << " " << md << " " << me << " " << mf << "   dz = " << dz << std::endl;
              std::cout << "     ma mb mc md me mf = " << ma << " " << mb << " " << mc << " " << md << " " << me << " " << mf << std::endl;
                std::cout << " X,Y X', Y' " << Xabcdef << " " << Yabcdef << " " << Xdefabc << " " << Ydefabc
                          << " this dz = " << (2./36.)*   Xabcdef * Ydefabc << "  summed dz = " << dz << std::endl;
            }

              }// for me
             }// for md
           }// for mc
          }// for mb
         }// for ma
         Z0_m += dz;
         std::cout << " xxDEBUG: abcdef " << a << " " << b << " " << c << " " << d << " " << e << " " << f
                   << "    dz = " << dz << "   z0 " << Z0_m << std::endl;
       }// for f
      }// for e
     }// for d
    }// for c
   }// for b
  }// for a
  std::cout << "By the way, hermitian? X: " << X.IsHermitian() << "   Y: " << Y.IsHermitian()
            << "  anti?  X: " << X.IsAntiHermitian() << "  Y: " << Y.IsAntiHermitian() << std::endl;
  std::cout << " 3body: X  " << X.ThreeBody.herm << "   Y  " << Y.ThreeBody.herm << std::endl;

  std::cout << "prove it X...  " << X.ThreeBody.GetME_pn(0,2,3,0,0,3,2,5,2)  << "   " << X.ThreeBody.GetME_pn(2,0,3,2,2,5,0,0,3) << std::endl;
  std::cout << "prove it Y...  " << Y.ThreeBody.GetME_pn(0,2,3,0,0,3,2,5,2)  << "   " << Y.ThreeBody.GetME_pn(2,0,3,2,5,2,0,0,3) << std::endl;
  std::cout << "Norms : " << X.ThreeBodyNorm() << "   " << Y.ThreeBodyNorm() << std::endl;

  std::cout << "hm...: " << X.ThreeBody.GetME_pn(1,2,3, 0,0,3,2,2,5 ) << "    " << Y.ThreeBody.GetME_pn(1,2,3, 0,0,3,2,2,5 ) << std::endl;

  double summed_error =  Z0_m - Z_J.ZeroBody;
  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << Z0_m << " " << Z_J.ZeroBody
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}









