
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
  if ( particle_rank > 20 )
  {
    Rando.ThreeBody.SwitchToPN_and_discard();
    std::default_random_engine generator(random_seed);
    double mean = 0;
    double stddev = 1;
    std::normal_distribution<double> distribution(mean,stddev);

    size_t nch = modelspace.GetNumberThreeBodyChannels();
    for (size_t ch=0; ch<nch; ch++)
    {
      ThreeBodyChannel& Tbc = modelspace.GetThreeBodyChannel(ch);
      size_t nkets = Tbc.GetNumberKets();
      for (size_t ibra=0; ibra<nkets; ibra++)
      {
        for (size_t iket=0; iket<=ibra; iket++)
        {
          double random_me = distribution(generator);
          Rando.ThreeBody.SetME_pn_PN_ch( ch,ch, ibra,iket, random_me );
        }
      }
    }

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
  }


  std::cout << "In  " << __func__  << "  norm of 1b : " << Rando.OneBodyNorm();
  if (particle_rank > 1) std::cout << "   norm of 2b " << Rando.TwoBodyNorm();
  if (particle_rank > 2) std::cout << "   norm of 3b " << Rando.ThreeBodyNorm();
  std::cout << std::endl;
  
  return Rando;
}


// Generate a random dagger operator, i.e. one that changes particle number by +1
// for use in testing the dagger commutator routines.
//DaggerOperator UnitTest::RandomDaggerOp(ModelSpace& modelspace, index_t Q)
Operator UnitTest::RandomDaggerOp(ModelSpace& modelspace, index_t Q)
{
//   Operator dag(modelspace,Q);
   Operator dag(modelspace);
   dag.SetNumberLegs(3);
   dag.SetQSpaceOrbit(Q);
   dag.SetNonHermitian();

   std::default_random_engine generator(random_seed);
   double mean = 0;
   double stddev = 1;
   std::normal_distribution<double> distribution(mean,stddev);

   // Only orbits with the same l,j,tz quantum numbers as the Q orbit
   // can be non-zero, so we fill a single column of the one-body matrix.
   Orbit& oQ = modelspace.GetOrbit(Q);
   for (auto i : dag.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
   {
     double random_me = distribution(generator);
     dag.OneBody(i,Q)= random_me;
   }

   size_t nch = modelspace.GetNumberTwoBodyChannels();
   for (size_t ch=0; ch<nch; ch++)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
     size_t nkets = tbc.GetNumberKets();
     for (size_t ibra=0; ibra<nkets; ibra++)
     {
       Ket& bra = tbc.GetKet(ibra);
       for ( auto k : modelspace.all_orbits )
       {
         Orbit& ok = modelspace.GetOrbit(k);
         // check whether |kQ> lives in this channel
         if ( (ok.tz2+oQ.tz2 != 2*tbc.Tz) or ((ok.l+oQ.l)%2 != tbc.parity) or ( (ok.j2+oQ.j2)<2*tbc.J ) or ( std::abs(ok.j2-oQ.j2)>2*tbc.J) ) continue;
//         if ( not tbc.CheckChannel_ket( &ok, &oQ) ) continue;
         
         double random_me = distribution(generator);
//         dag.TwoBody.SetTBME(ch,ch,bra.p,bra.q,k,Q, random_me);
         dag.ThreeLeg.SetME(ch,bra.p,bra.q,k, random_me);
       }
     }
//     std::cout << "ch = " << ch << "  matrix looks like " << std::endl << dag.ThreeLeg.GetMatrix(ch) << std::endl << std::endl;
   }

   return dag;
}

/*
Operator UnitTest::RandomDaggerOp(ModelSpace& modelspace, index_t Q)
{
   Operator dag(modelspace);
   dag.SetNumberLegs(3);
   dag.SetQSpaceOrbit(Q);
   dag.SetNonHermitian();

   std::default_random_engine generator(random_seed);
   double mean = 0;
   double stddev = 1;
   std::normal_distribution<double> distribution(mean,stddev);

   // Only orbits with the same l,j,tz quantum numbers as the Q orbit
   // can be non-zero, so we fill a single column of the one-body matrix.
   Orbit& oQ = modelspace.GetOrbit(Q);
   for (auto i : dag.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
   {
     double random_me = distribution(generator);
     dag.OneBody(i,Q)= random_me;
   }

   size_t nch = modelspace.GetNumberTwoBodyChannels();
   for (size_t ch=0; ch<nch; ch++)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
     size_t nkets = tbc.GetNumberKets();
     for (size_t ibra=0; ibra<nkets; ibra++)
     {
       Ket& bra = tbc.GetKet(ibra);
       for ( auto k : modelspace.all_orbits )
       {
         Orbit& ok = modelspace.GetOrbit(k);
         // check whether |kQ> lives in this channel
         if ( not tbc.CheckChannel_ket( &ok, &oQ) ) continue;
         
         double random_me = distribution(generator);
         dag.TwoBody.SetTBME(ch,ch,bra.p,bra.q,k,Q, random_me);
       }
     }
     std::cout << "ch = " << ch << "  matrix looks like " << std::endl << dag.TwoBody.GetMatrix(ch,ch) << std::endl << std::endl;
   }

   return dag;
}

*/



void UnitTest::Test3BodyHermiticity(Operator& Y)
{

  std::cout << " Hermitian? " << Y.IsHermitian() << std::endl;

  int a=0,b=0,c=10,d=0,e=10,f=0;
  double yabcdef = Y.ThreeBody.GetME_pn(0,0,1,a,b,c,d,e,f);
  double ydefabc = Y.ThreeBody.GetME_pn(0,0,1,d,e,f,a,b,c);
  std::cout << " yabcdef = " << yabcdef << "   ydefabc = " << ydefabc << std::endl;
//  int ma=1,mb=-1,mc=1,md=1,me=-1,mf=1;
//  double yabcdef = GetMschemeMatrixElement_3b(Y, a,ma, b,mb, c,mc, d,md, e,me, f,mf)

}

//void UnitTest::Test3BodyAntisymmetry()
void UnitTest::Test3BodyAntisymmetry(Operator& Y)
{

  arma::arma_rng::set_seed( random_seed );
//  Operator Y = RandomOp(*modelspace, 0, 0, 0, 3, 1);



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

//  double yJ0_abcdef = Y.ThreeBody.GetME_pn(0,0,1,a,b,c,d,e,f);
//  double yJ1_abcdef = Y.ThreeBody.GetME_pn(1,0,1,a,b,c,d,e,f);
//  double yJ0_cbadef = Y.ThreeBody.GetME_pn(0,0,1,c,b,a,d,e,f);
//  double yJ1_cbadef = Y.ThreeBody.GetME_pn(1,0,1,c,b,a,d,e,f);

//  std::cout << "yJ0_abcdef = " << yJ0_abcdef << std::endl;
//  std::cout << "yJ1_abcdef = " << yJ1_abcdef << std::endl;
//  std::cout << "yJ0_cbadef = " << yJ0_cbadef << std::endl;
//  std::cout << "    should be - " << AngMom::phase(oa.j2+ob.j2+oc.j2) << " * 1 * 1 * " << AngMom::SixJ(0.5,0.5,0,0.5,0.5,0) << " * " << yJ0_abcdef
//            << " = " << -AngMom::phase(oa.j2+ob.j2+oc.j2) * AngMom::SixJ(0.5,0.5,0,0.5,0.5,0) * yJ0_abcdef << std::endl;
//  std::cout << "yJ1_cbadef = " << yJ1_cbadef << std::endl;
//  std::cout << "    should be - " << AngMom::phase(oa.j2+ob.j2+oc.j2) << " * sqrt(3) * 1 * " << AngMom::SixJ(0.5,0.5,1,0.5,0.5,0) << " * " << yJ0_abcdef
//            << " = " << -AngMom::phase(oa.j2+ob.j2+oc.j2) * sqrt(3) * AngMom::SixJ(0.5,0.5,1,0.5,0.5,0) * yJ0_abcdef << std::endl;

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
  std::cout << "@bca" << std::endl;
  double ybcadef = GetMschemeMatrixElement_3b(Y, b,mb, c,mc, a,ma, d,md, e,me, f,mf);
  std::cout << "@cab" << std::endl;
  double ycabdef = GetMschemeMatrixElement_3b(Y, c,mc, a,ma, b,mb, d,md, e,me, f,mf);

  std::cout << "@acb" << std::endl;
  double yacbdef = GetMschemeMatrixElement_3b(Y, a,ma, c,mc, b,mb, d,md, e,me, f,mf);
  std::cout << "@bac" << std::endl;
  double ybacdef = GetMschemeMatrixElement_3b(Y, b,mb, a,ma, c,mc, d,md, e,me, f,mf);
  std::cout << std::endl << "@cba" << std::endl;
  double ycbadef = GetMschemeMatrixElement_3b(Y, c,mc, b,mb, a,ma, d,md, e,me, f,mf);

  std::cout << " @@mvals: " << ma << " " << mb << " " << mc << " " << md << " " << me << " " << mf << std::endl;
  std::cout << "ANTISYMMETRY:  abc: " << yabcdef << " bca: " << ybcadef << " cab: " << ycabdef
            << "    acb: " << yacbdef << " bac: " << ybacdef << " cba: " << ycbadef << std::endl;
//  std::cout << "ANTISYMMETRY:  abc: " << yabcdef  << " cba: " << ycbadef << std::endl;
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
  double t_start = omp_get_wtime();
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

  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
}


//void UnitTest::TestCommutators3()
void UnitTest::TestCommutators3(Operator& X, Operator& Y)
{
  double t_start = omp_get_wtime();
//  std::cout << " random_seed = " << random_seed << std::endl;
//  arma::arma_rng::set_seed( random_seed );
//  modelspace->PreCalculateSixJ();
//  Operator X = RandomOp(*modelspace, 0, 0, 0, 3, -1);
//  random_seed ++;
//  Operator Y = RandomOp(*modelspace, 0, 0, 0, 3, +1);
//  random_seed --;

  X += RandomOp(*modelspace, 0, 0, 0, 2, -1);
  Y += RandomOp(*modelspace, 0, 0, 0, 2, +1);
  bool all_good = true;

//  all_good &= Test_comm330ss( X, Y );
//  all_good &= Test_comm331ss( X, Y );
//  all_good &= Test_comm231ss( X, Y );
//  all_good &= Test_comm132ss( X, Y );
  all_good &= Test_comm223ss( X, Y );

  if ( all_good )
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }
  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
}





void UnitTest::TestDaggerCommutators(index_t Q)
{
  double t_start = omp_get_wtime();
  arma::arma_rng::set_seed( random_seed );
  Operator X = RandomOp(*modelspace, 0, 0, 0, 2, -1);
  Operator Y = RandomDaggerOp(*modelspace, Q);

  bool all_good = true;

  all_good &= Test_comm211sd( X, Y );
  all_good &= Test_comm231sd( X, Y );
  all_good &= Test_comm431sd( X, Y );
  all_good &= Test_comm413sd( X, Y );
  all_good &= Test_comm233sd( X, Y );
  all_good &= Test_comm433_pp_hh_sd( X, Y );
  all_good &= Test_comm433sd_ph( X, Y );

  if ( all_good )
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }
  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
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

  if (a==b) Jab_max = std::min( Jab_max,  oa.j2-1);
  if (d==e) Jde_max = std::min( Jde_max,  od.j2-1);

  for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
  {
    if (a==b and (Jab%2)>0 ) continue;
    double clebsch_ab = AngMom::CG(0.5*oa.j2, 0.5*ma,  0.5*ob.j2, 0.5*mb,  Jab, Mab );
    if ( std::abs(clebsch_ab)<1e-6 ) continue;
    for (int Jde=Jde_min; Jde<=Jde_max; Jde++)
    {
      if (d==e and (Jde%2)>0 ) continue;
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
        if ((a==0 and b==0 and c==10 and d==0 and e==10 and f==0 ) ) 
        {
        std::cout << "$abc: " << a << " " << b << " " << c << " def: " << d << " " << e << " " << f << std::endl;
        std::cout << "$m vals: " << ma << " " << mb << " " << mc << "  " << md << " " << me << " " << mf << std::endl;
        std::cout << "        Jab Jde twoJ " << Jab << " " << Jde << " " << twoJ
                  << " clebsch: " << clebsch_ab << " " << clebsch_de << " " << clebsch_abc << " " << clebsch_def
                  << "   matel_J " << meJ 
                  << "  matel = " << matel << std::endl;
        }
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




// The following two routines implicitly assume that we've chosen the apropriate projection quantum numbers.
// This may or may not be the most straightforward way to do things.
double UnitTest::GetMschemeMatrixElement_1leg( const Operator& Op, int a, int ma )
{
  
  return GetMschemeMatrixElement_1b( Op, a, ma, Op.GetQSpaceOrbit(), ma) ;
}


double UnitTest::GetMschemeMatrixElement_3leg( const Operator& Op, int a, int ma, int b, int mb, int c, int mc )
{

  double matel = 0;
  int Jop = Op.GetJRank();
  int Q = Op.GetQSpaceOrbit();
  Orbit& oa = Op.modelspace->GetOrbit(a);
  Orbit& ob = Op.modelspace->GetOrbit(b);
  Orbit& oc = Op.modelspace->GetOrbit(c);
  Orbit& oQ = Op.modelspace->GetOrbit(Q);
  int mQ = ma+mb-mc;
  if (mQ > oQ.j2) return 0;
  if (a==b and ma==mb) return 0;
//  if (c==d and mc==md) return 0;
  if ( Jop == 0) // scalar operator
  {
    if ( (ma+mb) == (mc+mQ) )
    {
      int Jmin = std::max(  std::abs( oa.j2-ob.j2) , std::abs(oc.j2-oQ.j2)  ) / 2;
      int Jmax = std::min( oa.j2+ob.j2  ,  oc.j2+oQ.j2 ) / 2;
      int M = (ma+mb)/2;
      for (int J=Jmin; J<=Jmax; J++)
      {
        // We take the un-normalized TBME (the one with the tilde) so we don't
        // need to worry about normalization factors.
        double clebsch_ab = AngMom::CG( 0.5*oa.j2, 0.5*ma, 0.5*ob.j2, 0.5*mb,  J, M );
        double clebsch_cd = AngMom::CG( 0.5*oc.j2, 0.5*mc, 0.5*oQ.j2, 0.5*mQ,  J, M );

        matel += clebsch_ab * clebsch_cd * Op.ThreeLeg.GetME_J(J,a,b,c);


//        if (a==0 and b==1 and c==1)
//        {
//
//         std::cout << "*******  JM=" << J << " " << M << "  clebsch = " << clebsch_ab << " " << clebsch_cd <<  " * " << Op.ThreeLeg.GetME_J(J,a,b,c) << "  ->  " << clebsch_ab * clebsch_cd * Op.ThreeLeg.GetME_J(J,a,b,c) << "   matel = " << matel << std::endl;
//        }
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
// THIS ONE DOESN'T WORK YET!!!!  (OK, it works now. All is well...)
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
  size_t norbits = X.modelspace->GetNumberOrbits();
  #pragma omp parallel for schedule(dynamic,1) reduction(+:Z0_m)
  for (size_t a=0; a<norbits; a++ )
//  for (auto a : X.modelspace->all_orbits )
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

                dz += (1./36) * na*nb*nc*(1-nd)*(1-ne)*(1-nf) * (Xabcdef * Ydefabc - Yabcdef * Xdefabc );
//                dz += (2./36) * na*nb*nc*(1-nd)*(1-ne)*(1-nf) * (Xabcdef * Ydefabc  );

              }// for me
             }// for md
           }// for mc
          }// for mb
         }// for ma
         Z0_m += dz;
       }// for f
      }// for e
     }// for d
    }// for c
   }// for b
  }// for a

  double summed_error =  Z0_m - Z_J.ZeroBody;
  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << Z0_m << " " << Z_J.ZeroBody
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}



/// M-Scheme Formula:
//
// Z_ij = 1/4 sum_abcde (nanb n`c n`dn`e) (X_abicde Y_cdeabj - Y_abicde X_cdeabj)
//
// this is very slow...
bool UnitTest::Test_comm331ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();

  Commutator::comm331ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize( );

  std::cout << __func__ <<  "   Done with J-scheme commutator. Begin m-scheme version..." << std::endl;

  double summed_error =  0;
  double sum_m = 0;
  double sum_J = 0;

  for ( auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for ( auto j : Z_J.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
//      Orbit& oj = X.modelspace->GetOrbit(j);
      int mj = mi;
      double Zm_ij = 0;
      size_t norb = X.modelspace->GetNumberOrbits();
//      for ( auto a : X.modelspace->all_orbits )
      #pragma omp parallel for schedule(dynamic,1) reduction(+:Zm_ij)
      for (size_t a=0; a<norb; a++)
      {
        Orbit& oa = X.modelspace->GetOrbit(a);
        double na = oa.occ;
        for ( auto b : X.modelspace->all_orbits )
        {
          Orbit& ob = X.modelspace->GetOrbit(b);
          double nb = ob.occ;
          for ( auto c : X.modelspace->all_orbits )
          {
            Orbit& oc = X.modelspace->GetOrbit(c);
            double nc = oc.occ;
            for ( auto d : X.modelspace->all_orbits )
            {
              Orbit& od = X.modelspace->GetOrbit(d);
              double nd = od.occ;
              for ( auto e : X.modelspace->all_orbits )
              {
                Orbit& oe = X.modelspace->GetOrbit(e);
                double ne = oe.occ;
                for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
                {
                 for (int mb=-ob.j2; mb<=ob.j2; mb+=2)
                 {
                  for (int mc=-oc.j2; mc<=oc.j2; mc+=2)
                  {
                   for (int md=-od.j2; md<=od.j2; md+=2)
                   {
                    int me = ma+mb+mi - mc-md;
// Z_ij = 1/4 sum_abcde (nanb n`c n`dn`e) (X_abicde Y_cdeabj - Y_abicde X_cdeabj)
                    double xabicde = GetMschemeMatrixElement_3b( X, a,ma, b,mb, i,mi, c,mc, d,md, e,me );
                    double ycdeabj = GetMschemeMatrixElement_3b( Y, c,mc, d,md, e,me, a,ma, b,mb, j,mj );
                    double xcdeabj = GetMschemeMatrixElement_3b( X, c,mc, d,md, e,me, a,ma, b,mb, j,mj );
                    double yabicde = GetMschemeMatrixElement_3b( Y, a,ma, b,mb, i,mi, c,mc, d,md, e,me );
                    Zm_ij += 1./4 * na*nb*(1-nc)*(1-nd)*(1-ne) * ( xabicde * ycdeabj - yabicde * xcdeabj );
                   }// for md
                  }// for mc
                 }// for mb
                }// for ma
              }// for e
            } // for d
          } // for c
        } // for b
      }// for a

      double ZJ_ij = Z_J.OneBody(i,j);
      double err = Zm_ij - ZJ_ij;
      if (std::abs(err)>1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j 
                  << "   zij = " << Zm_ij << "   ZJ_iij = " << ZJ_ij << "   err = " << err << std::endl; 
      }
      summed_error += err*err;
      sum_m += Zm_ij*Zm_ij;
      sum_J += ZJ_ij*ZJ_ij;

    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}





/// M-Scheme Formula:
//
//  Z_ij = 1/4 sum_abcd (nanb n`cn`d) (X_abcd Y_cdiabj - Y_abicdj X_cdab)
//
bool UnitTest::Test_comm231ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();

  Commutator::comm231ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize( );

  std::cout << __func__ <<  "   Done with J-scheme commutator. Begin m-scheme version..." << std::endl;

  double summed_error =  0;
  double sum_m = 0;
  double sum_J = 0;

  for ( auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for ( auto j : Z_J.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
//      Orbit& oj = X.modelspace->GetOrbit(j);
      int mj = mi;
      double Zm_ij = 0;

      for ( auto a : X.modelspace->all_orbits )
      {
        Orbit& oa = X.modelspace->GetOrbit(a);
        double na = oa.occ;
        for ( auto b : X.modelspace->all_orbits )
        {
          Orbit& ob = X.modelspace->GetOrbit(b);
          double nb = ob.occ;
          for ( auto c : X.modelspace->all_orbits )
          {
            Orbit& oc = X.modelspace->GetOrbit(c);
            double nc = oc.occ;
            for ( auto d : X.modelspace->all_orbits )
            {
              Orbit& od = X.modelspace->GetOrbit(d);
              double nd = od.occ;
              for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
              {
               for (int mb=-ob.j2; mb<=ob.j2; mb+=2)
               {
                for (int mc=-oc.j2; mc<=oc.j2; mc+=2)
                {
                 for (int md=-od.j2; md<=od.j2; md+=2)
                 {
                      if ( (ma+mb) != (mc+md) ) continue;
                      double xabcd = GetMschemeMatrixElement_2b( X, a,ma, b,mb, c,mc, d,md );
                      double xcdab = GetMschemeMatrixElement_2b( X, c,mc, d,md, a,ma, b,mb );
                      double yabicdj = GetMschemeMatrixElement_3b( Y, a,ma, b,mb, i,mi, c,mc, d,md, j,mj );
                      double ycdiabj = GetMschemeMatrixElement_3b( Y, c,mc, d,md, i,mi, a,ma, b,mb, j,mj );

                      double yabcd = GetMschemeMatrixElement_2b( Y, a,ma, b,mb, c,mc, d,md );
                      double ycdab = GetMschemeMatrixElement_2b( Y, c,mc, d,md, a,ma, b,mb );
                      double xabicdj = GetMschemeMatrixElement_3b( X, a,ma, b,mb, i,mi, c,mc, d,md, j,mj );
                      double xcdiabj = GetMschemeMatrixElement_3b( X, c,mc, d,md, i,mi, a,ma, b,mb, j,mj );

                      Zm_ij += (1./4) * na*nb*(1-nc)*(1-nd) * ( ( xabcd*ycdiabj - yabicdj*xcdab )
                                                              - ( yabcd*xcdiabj - xabicdj*ycdab ) );

                 }// for md
                }// for mc
               }// for mb
              }// for ma

            } // for d
          } // for c
        } // for b
      }// for a


      double ZJ_ij = Z_J.OneBody(i,j);
      double err = Zm_ij - ZJ_ij;
      if (std::abs(err)>1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j 
                  << "   zij = " << Zm_ij << "   ZJ_iij = " << ZJ_ij << "   err = " << err << std::endl; 
      }
      summed_error += err*err;
      sum_m += Zm_ij*Zm_ij;
      sum_J += ZJ_ij*ZJ_ij;

    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}









/// M-Scheme Formula:
//
// Z_ijkl = sum_ab (na-nb) (X_ab * Y_ijbkla)
//
bool UnitTest::Test_comm132ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();


  Commutator::comm132ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  std::cout << "X one body : " << std::endl << X.OneBody << std::endl;
  std::cout << "Y one body : " << std::endl << Y.OneBody << std::endl;


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
             if ( std::abs(mk)>ol.j2 ) continue; // this just eliminates some redundant matrix elements (I think?)
             int ml = mi + mj - mk;
             if ( std::abs(ml)>ol.j2) continue;

             double Zm_ijkl = 0;
             for (auto a : X.modelspace->all_orbits )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
//              for (auto b : X.modelspace->all_orbits )
              for (auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
              {
                Orbit& ob = X.modelspace->GetOrbit(b);
                double nb = ob.occ;

                if ( (oi.l+oj.l+ob.l+ok.l+ol.l+oa.l)%2>0) continue;
//                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
//                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;
                if ( (oi.tz2+oj.tz2+ob.tz2) != (ok.tz2+ol.tz2+oa.tz2) ) continue;

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
//                   int mb = mi + mj - ma;
                   int mb = ma;
                
                  double xab = GetMschemeMatrixElement_1b( X, a,ma, b,mb );
                  double yab = GetMschemeMatrixElement_1b( Y, a,ma, b,mb );
                  double xijbkla = GetMschemeMatrixElement_3b( X, i,mi, j,mj, b,mb, k,mk, l,ml, a,ma );
                  double yijbkla = GetMschemeMatrixElement_3b( Y, i,mi, j,mj, b,mb, k,mk, l,ml, a,ma );

//                  Zm_ijkl +=  ( na*(1-nb) - nb*(1-na) ) * (xab * yijbkla - yab * xijbkla  );
//                  Zm_ijkl +=  ( na - nb ) * (xab * yijbkla - yab * xijbkla  );
//                 Zm_ijkl +=  ( na - nb ) * (xab * yijbkla   );
                  Zm_ijkl +=  ( na - nb ) * (-yab * xijbkla   );
                  if ( i==0 and j==0 and k==0 and l==10 and mi==mk and std::abs(na-nb)>1e-5 )
                  {
                   std::cout << " X  a,b = " << a << " " << b << " {m} = " << ma << " " << mb << "  yab = " << yab << " xijbkla = " << xijbkla << "   z= " << Zm_ijkl << std::endl;
                   std::cout << " Y  a,b = " << a << " " << b << " {m} = " << ma << " " << mb << "  xab = " << xab << " yijbkla = " << yijbkla << "   z= " << Zm_ijkl << std::endl;
                   std::cout << "    herm. conj. of X : " << GetMschemeMatrixElement_3b( X, k,mk, l,ml, a,ma, i,mi, j,mj, b,mb ) << std::endl;
                   std::cout << "    herm. conj. of Y : " << GetMschemeMatrixElement_3b( Y, k,mk, l,ml, a,ma, i,mi, j,mj, b,mb ) << std::endl;
                  }
                }// for ma
              }// for b
             }// for a


             double ZJ_ijkl = GetMschemeMatrixElement_2b( Z_J, i,mi, j,mj, k,mk, l,ml ) ;
             double err = Zm_ijkl - ZJ_ijkl;
             if (std::abs(err)>1e-6)
             {
               std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                         << " {m} = " << mi << " " << mj << " " << mk << " " << ml 
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


//
// M-scheme formula
//
// Z_ijklmn = (1-Pik-Pjk)(1-Plm-Pln) sum_a (X_ijla Y_akmn - Y_ijla X_akmn)
//
bool UnitTest::Test_comm223ss( const Operator& X, const Operator& Y )
{
  Operator Z_J( *(Y.modelspace), 0,0,0,3 );
//  Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();


  Commutator::comm223ss( X, Y, Z_J);

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
        if (k<j) continue;


        for (auto l : X.modelspace->all_orbits )
        {
          if (l<i) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          for ( auto m : X.modelspace->all_orbits )
          {
            if (m<l) continue;
            Orbit& om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits )
            {
              if (n<m) continue;
              Orbit& on = X.modelspace->GetOrbit(n);
              if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2 !=0 ) continue;
              if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;

              // loop over projections
              for (int m_i=-oi.j2; m_i<=oi.j2; m_i+=2)
              {
               for (int m_j=-oj.j2; m_j<=oj.j2; m_j+=2)
               {
                for (int m_k=-ok.j2; m_k<=ok.j2; m_k+=2)
                {
                 for (int m_l=-ol.j2; m_l<=ol.j2; m_l+=2)
                 {
                  for (int m_m=-om.j2; m_m<=om.j2; m_m+=2)
                  {
                   for (int m_n=-on.j2; m_n<=on.j2; m_n+=2)
                   {
                     if ( (m_i+m_j+m_k) != (m_l+m_m+m_n) ) continue;
                     double z_ijklmn = 0;
                     for ( auto a : X.modelspace->all_orbits )
                     {
                       Orbit& oa = X.modelspace->GetOrbit(a);
                       for ( int m_a=-oa.j2; m_a<=oa.j2; m_a++)
                       {
                         // "direct" term
                         double x_ijla = GetMschemeMatrixElement_2b( X, i,m_i, j,m_j, l,m_l, a,m_a );
                         double y_akmn = GetMschemeMatrixElement_2b( Y, a,m_a, k,m_k, m,m_m, n,m_n );
                         double y_ijla = GetMschemeMatrixElement_2b( Y, i,m_i, j,m_j, l,m_l, a,m_a );
                         double x_akmn = GetMschemeMatrixElement_2b( X, a,m_a, k,m_k, m,m_m, n,m_n );
                         // Pik
                         double x_kjla = GetMschemeMatrixElement_2b( X, k,m_k, j,m_j, l,m_l, a,m_a );
                         double y_aimn = GetMschemeMatrixElement_2b( Y, a,m_a, i,m_i, m,m_m, n,m_n );
                         double y_kjla = GetMschemeMatrixElement_2b( Y, k,m_k, j,m_j, l,m_l, a,m_a );
                         double x_aimn = GetMschemeMatrixElement_2b( X, a,m_a, i,m_i, m,m_m, n,m_n );
                         // Pjk
                         double x_ikla = GetMschemeMatrixElement_2b( X, i,m_i, k,m_k, l,m_l, a,m_a );
                         double y_ajmn = GetMschemeMatrixElement_2b( Y, a,m_a, j,m_j, m,m_m, n,m_n );
                         double y_ikla = GetMschemeMatrixElement_2b( Y, i,m_i, k,m_k, l,m_l, a,m_a );
                         double x_ajmn = GetMschemeMatrixElement_2b( X, a,m_a, j,m_j, m,m_m, n,m_n );
                         // Plm
                         double x_ijma = GetMschemeMatrixElement_2b( X, i,m_i, j,m_j, m,m_m, a,m_a );
                         double y_akln = GetMschemeMatrixElement_2b( Y, a,m_a, k,m_k, l,m_l, n,m_n );
                         double y_ijma = GetMschemeMatrixElement_2b( Y, i,m_i, j,m_j, m,m_m, a,m_a );
                         double x_akln = GetMschemeMatrixElement_2b( X, a,m_a, k,m_k, l,m_l, n,m_n );
                         // Pln
                         double x_ijna = GetMschemeMatrixElement_2b( X, i,m_i, j,m_j, n,m_n, a,m_a );
                         double y_akml = GetMschemeMatrixElement_2b( Y, a,m_a, k,m_k, m,m_m, l,m_l );
                         double y_ijna = GetMschemeMatrixElement_2b( Y, i,m_i, j,m_j, n,m_n, a,m_a );
                         double x_akml = GetMschemeMatrixElement_2b( X, a,m_a, k,m_k, m,m_m, l,m_l );
                         // Pik Plm
                         double x_kjma = GetMschemeMatrixElement_2b( X, k,m_k, j,m_j, m,m_m, a,m_a );
                         double y_ailn = GetMschemeMatrixElement_2b( Y, a,m_a, i,m_i, l,m_l, n,m_n );
                         double y_kjma = GetMschemeMatrixElement_2b( Y, k,m_k, j,m_j, m,m_m, a,m_a );
                         double x_ailn = GetMschemeMatrixElement_2b( X, a,m_a, i,m_i, l,m_l, n,m_n );
                         // Pik Pln
                         double x_kjna = GetMschemeMatrixElement_2b( X, k,m_k, j,m_j, n,m_n, a,m_a );
                         double y_aiml = GetMschemeMatrixElement_2b( Y, a,m_a, i,m_i, m,m_m, l,m_l );
                         double y_kjna = GetMschemeMatrixElement_2b( Y, k,m_k, j,m_j, n,m_n, a,m_a );
                         double x_aiml = GetMschemeMatrixElement_2b( X, a,m_a, i,m_i, m,m_m, l,m_l );
                         // Pjk Plm
                         double x_ikma = GetMschemeMatrixElement_2b( X, i,m_i, k,m_k, m,m_m, a,m_a );
                         double y_ajln = GetMschemeMatrixElement_2b( Y, a,m_a, j,m_j, l,m_l, n,m_n );
                         double y_ikma = GetMschemeMatrixElement_2b( Y, i,m_i, k,m_k, m,m_m, a,m_a );
                         double x_ajln = GetMschemeMatrixElement_2b( X, a,m_a, j,m_j, l,m_l, n,m_n );
                         // Pjk Pln
                         double x_ikna = GetMschemeMatrixElement_2b( X, i,m_i, k,m_k, n,m_n, a,m_a );
                         double y_ajml = GetMschemeMatrixElement_2b( Y, a,m_a, j,m_j, m,m_m, l,m_l );
                         double y_ikna = GetMschemeMatrixElement_2b( Y, i,m_i, k,m_k, n,m_n, a,m_a );
                         double x_ajml = GetMschemeMatrixElement_2b( X, a,m_a, j,m_j, m,m_m, l,m_l );
                         z_ijklmn += x_ijla * y_akmn - y_ijla * x_akmn;
                         z_ijklmn -= x_kjla * y_aimn - y_kjla * x_aimn;
                         z_ijklmn -= x_ikla * y_ajmn - y_ikla * x_ajmn;
                         z_ijklmn -= x_ijma * y_akln - y_ijma * x_akln;
                         z_ijklmn -= x_ijna * y_akml - y_ijna * x_akml;
                         z_ijklmn += x_kjma * y_ailn - y_kjma * x_ailn;
                         z_ijklmn += x_kjna * y_aiml - y_kjna * x_aiml;
                         z_ijklmn += x_ikma * y_ajln - y_ikma * x_ajln;
                         z_ijklmn += x_ikna * y_ajml - y_ikna * x_ajml;
                       }
                     }
                    
                    double ZJ_ijklmn = GetMschemeMatrixElement_3b( Z_J, i,m_i, j,m_j, k,m_k, l,m_l, m,m_m, n,m_n );
                    double err = z_ijklmn - ZJ_ijklmn;
                    if (std::abs(err)>1e-6)
                    {
                      std::cout << "Trouble in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl; 
                    }
                    summed_error += err*err;
                    sum_m += z_ijklmn*z_ijklmn;
                    sum_J += ZJ_ijklmn*ZJ_ijklmn;
                   } // for m_n
                  } // for m_m
                 } // for m_l
                } // for m_k
               } // for m_j
              } // for m_i
              
            } // for n
          } // for m
        } // for l
      } // for k
    } // for j
  } // for i



  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
//  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

 
}


//***************************************************************
//// Now we move on to the scalar-dagger commutator expressions.
//***************************************************************


/// M-Scheme Formula:
//
//  Z_i = sum_a X_ia Y_a
//
bool UnitTest::Test_comm211sd( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.SetNonHermitian();
  Z_J.Erase();


  Commutator::comm211sd( X, Y, Z_J);

//  std::cout << "X one body : " << std::endl << X.OneBody << std::endl;
//  std::cout << "Y one body : " << std::endl << Y.OneBody << std::endl;


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  size_t Q = Z_J.GetQSpaceOrbit();
//  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);


  for (auto i : Z_J.modelspace->all_orbits )
  {
    Orbit& oi = Z_J.modelspace->GetOrbit(i);
    int mi = oi.j2;
    double Zm_i = 0;
    for (auto a : X.modelspace->all_orbits)
    {
     Orbit& oa = X.modelspace->GetOrbit(a);
     for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
     {
       double Xia = GetMschemeMatrixElement_1b( X, i, mi, a, ma );
       double Ya = GetMschemeMatrixElement_1leg( Y, a, ma );

       Zm_i += Xia * Ya ;
     }
    }
    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err)>1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i <<  "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl; 
    }
    summed_error += err*err;
    sum_m += Zm_i*Zm_i;
    sum_J += ZJ_i*ZJ_i;
  }





  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}



/// M-Scheme Formula:
//
//  Z_i = sum_ab (na-nb)X_ab Y_bia
//
bool UnitTest::Test_comm231sd( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.SetNonHermitian();
  Z_J.Erase();


  Commutator::comm231sd( X, Y, Z_J);


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  size_t Q = Z_J.GetQSpaceOrbit();
//  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);



  for (auto i : Z_J.modelspace->all_orbits )
  {
    Orbit& oi = Z_J.modelspace->GetOrbit(i);
    int mi = oi.j2;
    double Zm_i = 0;
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
         for (int mb=-ob.j2; mb<=ob.j2;mb+=2)
         {
           double Xab = GetMschemeMatrixElement_1b( X, a, ma, b,mb );
           double Ybia = GetMschemeMatrixElement_3leg( Y, b,mb, i,mi, a,ma );

           Zm_i += (na-nb) * Xab * Ybia ;
         }
       }
     }
    }
    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err)>1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i <<  "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl; 
    }
    summed_error += err*err;
    sum_m += Zm_i*Zm_i;
    sum_J += ZJ_i*ZJ_i;
  }



  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}




/// M-Scheme Formula:
//
//  Z_i = 1/2 sum_abc (na nb`nc`-na`nb nc) X_aibc Y_bca
//
bool UnitTest::Test_comm431sd( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.SetNonHermitian();
  Z_J.Erase();


  Commutator::comm433_pp_hh_431sd( X, Y, Z_J);


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  size_t Q = Z_J.GetQSpaceOrbit();
//  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);


  for (auto i : Z_J.modelspace->all_orbits )
  {
    Orbit& oi = Z_J.modelspace->GetOrbit(i);
    int mi = oi.j2;
    double Zm_i = 0;
    for (auto a : X.modelspace->all_orbits)
    {
     Orbit& oa = X.modelspace->GetOrbit(a);
     double na = oa.occ;
     for (auto b : X.modelspace->all_orbits)
     {
       Orbit& ob = X.modelspace->GetOrbit(b);
       double nb = ob.occ;
       for (auto c : X.modelspace->all_orbits)
       {
         Orbit& oc = X.modelspace->GetOrbit(c);
         double nc = oc.occ;
         for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
         {
           for (int mb=-ob.j2; mb<=ob.j2;mb+=2)
           {
             for (int mc=-oc.j2; mc<=oc.j2;mc+=2)
             {
               double Xaibc = GetMschemeMatrixElement_2b( X, a, ma, i,mi, b,mb, c,mc );
               double Ybca = GetMschemeMatrixElement_3leg( Y, b,mb, c,mc, a,ma );

               Zm_i += 1./2 * ( na*(1-nb)*(1-nc) + (1-na)*nb*nc ) * Xaibc * Ybca ;
             }
           }
         }
       }
     }
    }
    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err)>1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i <<  "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl; 
    }
    summed_error += err*err;
    sum_m += Zm_i*Zm_i;
    sum_J += ZJ_i*ZJ_i;
  }



  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}





/// M-Scheme Formula:
//
//
// 413:
// Zijk = sum_a Xijka*Ya
//
bool UnitTest::Test_comm413sd( const Operator& Xin, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.Erase();

  Operator X = Xin;
  X.EraseOneBody(); // we do this so that the 233 term does not contribute.

  Commutator::comm413_233sd( X, Y, Z_J);

  size_t Q = Z_J.GetQSpaceOrbit();
  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);

  int mQ = oQ.j2;

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
//          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
//          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {
             if ( (mi+mj)!=(mk+mQ) ) continue;
//             if ( std::abs(mk)>ol.j2 ) continue;
//             int ml = mi + mj - mk;

             double Zm_ijk = 0;
             for (auto a : X.modelspace->all_orbits )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
//              double na = oa.occ;

//                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
//                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
                  if (ma != mQ) continue;

                  double Xijka =  GetMschemeMatrixElement_2b( X, i,mi, j,mj, k,mk, a,ma ) ;
                  double Ya = GetMschemeMatrixElement_1leg(Y, a,ma) ;

                  Zm_ijk += Xijka * Ya; // 413 term
//                  if (i==0 and j==8 and k==0)
//                  {
//                    std::cout << "  #### a = " << a << "  m vals " << mi << " " << mj << " " << mk << " " << ma << "  " << mQ << "   Xijka , Ya = " << Xijka << " " << Ya << "   Zm_ijk = " << Zm_ijk << std::endl;
//                  }

                }// for ma
             }// for a


             double ZJ_ijk = GetMschemeMatrixElement_3leg( Z_J, i,mi, j,mj, k,mk ) ;
             double err = Zm_ijk - ZJ_ijk;
             if (std::abs(err)>1e-6)
             {
               std::cout << "Trouble in " << __func__ << "  i,j,k = " << i << " " << j << " " << k 
                         << "  mvals " << mi << " " << mj << " " << mk << "   " << mQ
                         << "   Zm_ijk = " << Zm_ijk << "   ZJ_ijk = " << ZJ_ijk << "   err = " << err << std::endl; 
             }
             summed_error += err*err;
             sum_m += Zm_ijk*Zm_ijk;
             sum_J += ZJ_ijk*ZJ_ijk;

           }// for mk
          }// for mj
      }// for k
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.ThreeLegNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}


/// M-Scheme Formula:
//
// 233:
// Zijk =  sum_a ( Xia*Yajk + Xja*Yjak - Yija * Xak )
//
bool UnitTest::Test_comm233sd( const Operator& X, const Operator& Yin )
{

  Operator Y = Yin;

  Operator Z_J( Y );
  Z_J.Erase();

  Y.EraseOneBody(); // we do this so the 413 term does not contribute

  Commutator::comm413_233sd( X, Y, Z_J);

//  if ( Z_J.IsHermitian() )
//     Z_J.Symmetrize();
//  else if (Z_J.IsAntiHermitian() )
//     Z_J.AntiSymmetrize();

  size_t Q = Z_J.GetQSpaceOrbit();
  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);

  int mQ = oQ.j2;

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
//          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
//          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {
             if ( (mi+mj)!=(mk+mQ) ) continue;
//             if ( std::abs(mk)>ol.j2 ) continue;
//             int ml = mi + mj - mk;

             double Zm_ijk = 0;
             for (auto a : X.modelspace->all_orbits )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
//              double na = oa.occ;

//                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
//                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
                  double Xia = GetMschemeMatrixElement_1b(X, i,mi, a,ma);
                  double Xja = GetMschemeMatrixElement_1b(X, j,mj, a,ma);
                  double Xak = GetMschemeMatrixElement_1b(X, a,ma, k,mk);

                  double Yajk = GetMschemeMatrixElement_3leg(Y, a,ma, j,mj, k,mk);
                  double Yiak = GetMschemeMatrixElement_3leg(Y, i,mi, a,ma, k,mk);
                  double Yija = GetMschemeMatrixElement_3leg(Y, i,mi, j,mj, a,ma);

                  Zm_ijk += Xia*Yajk + Xja*Yiak - Yija*Xak; // 223 term


                }// for ma
             }// for a


             double ZJ_ijk = GetMschemeMatrixElement_3leg( Z_J, i,mi, j,mj, k,mk ) ;
             double err = Zm_ijk - ZJ_ijk;
             if (std::abs(err)>1e-6)
             {
               std::cout << "Trouble in " << __func__ << "  i,j,k = " << i << " " << j << " " << k 
                         << "  mvals " << mi << " " << mj << " " << mk
                         << "   Zm_ijk = " << Zm_ijk << "   ZJ_ijk = " << ZJ_ijk << "   err = " << err << std::endl; 
             }
             summed_error += err*err;
             sum_m += Zm_ijk*Zm_ijk;
             sum_J += ZJ_ijk*ZJ_ijk;

           }// for mk
          }// for mj
      }// for k
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.ThreeLegNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}




/// M-Scheme Formula:
//
// 433 pp/hh:
// Zijk =  1/2 sum_ab (nanb - na`nb`)  Xijab*Yabk 
//
bool UnitTest::Test_comm433_pp_hh_sd( const Operator& X, const Operator& Yin )
{

  Operator Y = Yin;

  Operator Z_J( Y );
  Z_J.Erase();

  Commutator::comm433_pp_hh_431sd( X, Y, Z_J);

  size_t Q = Z_J.GetQSpaceOrbit();
  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);

  int mQ = oQ.j2;

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
//          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
//          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {
             if ( (mi+mj)!=(mk+mQ) ) continue;
//             if ( std::abs(mk)>ol.j2 ) continue;
//             int ml = mi + mj - mk;

             double Zm_ijk = 0;
             for (auto a : X.modelspace->all_orbits )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
              for (auto b : X.modelspace->all_orbits )
              {
               Orbit& ob = X.modelspace->GetOrbit(b);
               double nb = ob.occ;

//                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
//                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
                  for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
                  {
                    double Xijab = GetMschemeMatrixElement_2b(X, i,mi, j,mj, a,ma, b,mb);
                    double Yabk = GetMschemeMatrixElement_3leg(Y, a,ma, b,mb, k,mk);

                    Zm_ijk += 1./2 * ((1-na)*(1-nb) - na*nb) * Xijab * Yabk ;

                  }// for mb
                }// for ma
              }// for b
             }// for a


             double ZJ_ijk = GetMschemeMatrixElement_3leg( Z_J, i,mi, j,mj, k,mk ) ;
             double err = Zm_ijk - ZJ_ijk;
             if (std::abs(err)>1e-6)
             {
               std::cout << "Trouble in " << __func__ << "  i,j,k = " << i << " " << j << " " << k 
                         << "  mvals " << mi << " " << mj << " " << mk
                         << "   Zm_ijk = " << Zm_ijk << "   ZJ_ijk = " << ZJ_ijk << "   err = " << err << std::endl; 
             }
             summed_error += err*err;
             sum_m += Zm_ijk*Zm_ijk;
             sum_J += ZJ_ijk*ZJ_ijk;

           }// for mk
          }// for mj
      }// for k
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.ThreeLegNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}




/// M-Scheme Formula:
//
// 433 pp/hh:
// Zijk =  sum_ab (na - nb)  Xiakb Ybja
//
bool UnitTest::Test_comm433sd_ph( const Operator& X, const Operator& Yin )
{

  Operator Y = Yin;

  Operator Z_J( Y );
  Z_J.Erase();

  Commutator::comm433sd_ph( X, Y, Z_J);
//  Commutator::comm433sd_ph_dumbway( X, Y, Z_J);

  size_t Q = Z_J.GetQSpaceOrbit();
  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);

  int mQ = oQ.j2;

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
//          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
//          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;


//          if (not (i==0 and j==1 and k==1) ) continue;

          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {
             if ( (mi+mj)!=(mk+mQ) ) continue;
//             if ( std::abs(mk)>ol.j2 ) continue;
//             int ml = mi + mj - mk;

             double Zm_ijk = 0;
             for (auto a : X.modelspace->all_orbits )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
              for (auto b : X.modelspace->all_orbits )
              {
               Orbit& ob = X.modelspace->GetOrbit(b);
               double nb = ob.occ;

//                  if ( not (a==0 and b==10) or (a==10 and b==0) ) continue;
//                  if ( not ((a==1 and b==8) or (a==8 and b==1)) ) continue;
//                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
//                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
                  for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
                  {
                    double Xiakb = GetMschemeMatrixElement_2b(X, i,mi, a,ma, k,mk, b,mb);
                    double Ybja = GetMschemeMatrixElement_3leg(Y, b,mb, j,mj, a,ma);
                    double Xjakb = GetMschemeMatrixElement_2b(X, j,mj, a,ma, k,mk, b,mb);
                    double Ybia = GetMschemeMatrixElement_3leg(Y, b,mb, i,mi, a,ma);

                    Zm_ijk += (na-nb) * ( Xiakb * Ybja  -  Xjakb * Ybia ) ;
//                    if (i==0 and j==1 and k==1 and mi==1 and mj==-1 and mk==-1 and ( (a==0 and b==10) or (b==0 and a==10) ) )
//                    if (i==0 and j==1 and k==1 and mi==1 and mj==-1 and mk==-1 and std::abs(na-nb)>0 )
//                    if (i==0 and j==1 and k==1 and mi==1 and mj==-1 and mk==-1 and ( (a==1 and b==8) or (b==1 and a==8) ) )
//                    {
//                      std::cout << "  a,b " << a << " " << b << "   ma,mb " << ma << " " << mb << "   na nb  " << na << " " << nb
//                                << "  Xiakb Ybja " << Xiakb << " " << Ybja << "   Xjakb Ybia " << Xjakb << " " << Ybia << "  => Zm_ijk = " << Zm_ijk << std::endl;
//                    }
//                    Zm_ijk -= (na-nb) * (     Xjakb * Ybia ) ;

                  }// for mb
                }// for ma
              }// for b
             }// for a


             double ZJ_ijk = GetMschemeMatrixElement_3leg( Z_J, i,mi, j,mj, k,mk ) ;
             double err = Zm_ijk - ZJ_ijk;
             if (std::abs(err)>1e-6)
             {
               std::cout << "Trouble in " << __func__ << "  i,j,k = " << i << " " << j << " " << k 
                         << "  mvals " << mi << " " << mj << " " << mk
                         << "   Zm_ijk = " << Zm_ijk << "   ZJ_ijk = " << ZJ_ijk << "   err = " << err << std::endl; 
             }
             summed_error += err*err;
             sum_m += Zm_ijk*Zm_ijk;
             sum_J += ZJ_ijk*ZJ_ijk;

           }// for mk
          }// for mj
      }// for k
    }// for j
  }// for i

  bool passed = std::abs( summed_error ) <1e-6 ;
  std::string passfail = passed ? "PASS " : "FAIL";
  if ( Z_J.ThreeLegNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ <<   "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;

}






