
#include "UnitTest.hh"
#include "AngMom.hh"
#include <armadillo>
#include <random>
#include <string>
#include "imsrg_util.hh"


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
     dag.OneBody(i,0)= random_me;
//     dag.OneBody(i,Q)= random_me;
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



void UnitTest::Test3BodySetGet(Operator& Y)
{
  int i=2,j=0,k=0,l=2,m=0,n=0;
  int Jij=1;
  int Jlm=1;
  int twoJ=1;

  double meread = Y.ThreeBody.GetME_pn(Jij,Jlm,twoJ,i,j,k,l,m,n);
  double setval = 10.0;
  Y.ThreeBody.SetME_pn(Jij,Jlm,twoJ,i,j,k,l,m,n,setval);
  double getval = Y.ThreeBody.GetME_pn(Jij,Jlm,twoJ,i,j,k,l,m,n);

  std::cout << "Before setting, ME = " << meread << std::endl;
  std::cout << "Set to " << setval << std::endl;
  std::cout << "After setting, get " << getval << std::endl;

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

  Operator Xherm = Y;
  Xherm.TwoBody.Erase();
  Xherm.ThreeBody.Erase();
  X += RandomOp(*modelspace, 0, 0, 0, 2, -1);
  Y += RandomOp(*modelspace, 0, 0, 0, 2, +1);
  Xherm += RandomOp(*modelspace, 0, 0, 0, 2, +1);
  X.ThreeBody.Erase();
  Commutator::comm223ss( Xherm, Y, X); // Make the 3-body part of X equal to the commutator of 2 hermitian 2b operators
  bool all_good = true;

//  all_good &= Test_comm330ss( X, Y );
//  all_good &= Test_comm331ss( X, Y );
//  all_good &= Test_comm231ss( X, Y );
//  all_good &= Test_comm132ss( X, Y );
//  all_good &= Test_comm232ss( X, Y );
//  all_good &= Test_comm223ss( X, Y );
//  all_good &= Test_comm133ss( X, Y );

//  all_good &= Test_comm332_ppph_hhhpss( X, Y ); 
//  all_good &= Test_comm332_pphhss( X, Y );  

//  all_good &= Test_comm233_pp_hhss( X, Y );   
//  all_good &= Test_comm233_ph_ss( X, Y );  
//  all_good &= Test_comm333_ppp_hhh_ss( X, Y );  
  all_good &= Test_comm333_pph_hhp_ss( X, Y );  



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

  Orbit& oQ = Y.modelspace->GetOrbit(Q);
  
  Y *= 0;
  std::vector<Operator> Yn;
  std::vector<Operator> Zn;
  for (size_t n : Y.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
  {
    Y.OneBody(n,0) = 1.0;
    Yn.push_back( 0*Y );
    Yn.back().OneBody(n,0) = 1.0;
    Zn.push_back(  Commutator::Commutator( X, Yn.back() ) );
    Zn.back() = Commutator::Commutator( X, Zn.back() );
  }
  
  Operator Z = Commutator::Commutator( X, Y);
  Z = Commutator::Commutator( X, Z);

  std::cout << "Y 1b: " << std::endl << Y.OneBody << std::endl;
  std::cout << "Yn 1b: " << std::endl;
  for (auto y : Yn ) std::cout << y.OneBody << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "Z 1b: " << std::endl << Z.OneBody << std::endl;
  std::cout << "Zn 1b: " << std::endl;
  for (auto z : Zn ) std::cout << z.OneBody << std::endl;
  std::cout << std::endl;


  Operator Zdiff = Z;
  std::cout << "legs for Z and Zdiff : " << Z.legs << " " << Zdiff.legs << std::endl;
  for ( auto z : Zn ) Zdiff -= z;
  std::cout << "norm of Zdiff = " << Zdiff.Norm() << "  norm Z = " << Z.Norm()
            << "  norm of Zn = ";
  for (auto z : Zn ) std::cout << z.Norm() << " ";
  std::cout << std::endl;

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



void UnitTest::TestDaggerCommutatorsAlln(index_t Q)
{
  double t_start = omp_get_wtime();
  arma::arma_rng::set_seed( random_seed );
  Operator X = RandomOp(*modelspace, 0, 0, 0, 2, -1);
//  Operator Y = RandomDaggerOp(*modelspace, Q);

  Operator Y = imsrg_util::DaggerAlln_Op(*modelspace, Q);
  Operator Z = Commutator::Commutator(X, Commutator::Commutator(X,Y));
  std::vector<Operator> Yn;
  std::vector<Operator> Zn;
  Orbit& oQ = X.modelspace->GetOrbit(Q);
  for ( auto nQ : X.modelspace->OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
  {
    Yn.push_back( imsrg_util::Dagger_Op(*modelspace, nQ) );
//    Yn.back().OneBody(nQ,nQ) = 0;
//    Yn.back().OneBody(nQ,Q) = 1.0;
//    Yn.back().OneBody(nQ,0) = 0;
    Yn.back().OneBody(nQ,0) = 1.0;
    Zn.push_back( Commutator::Commutator( X, Commutator::Commutator(X,Yn.back()) ) );
  }
  std::cout << "Y one body: "<< std::endl << Y.OneBody << std::endl;
  std::cout << "Z one body: "<< std::endl << Z.OneBody << std::endl;
  for ( size_t n=0; n< Yn.size(); n++)
  {
    std::cout << " Y " << n << std::endl << Yn[n].OneBody << std::endl;
    std::cout << " Z " << n << std::endl << Zn[n].OneBody << std::endl;
  }
  Operator Zsum = Zn.front();
  for (size_t n=1; n<Zn.size(); n++)   Zsum += Zn[n];
  Operator Zdiff = Z - Zsum;
  double Znorm = Zdiff.Norm();
  std::cout << "Norms: " << Z.Norm() << "   vs  ";
  for (auto z : Zn) std::cout << z.Norm() << " ";
  std::cout << std::endl << "Norm of Zdiff = " << Znorm << std::endl;

  

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
//        std::cout << " line " << __LINE__ << "   J1,J2,J = " << Jab << " " << Jde << " " << twoJ << "  cgs: " << clebsch_ab << " " << clebsch_abc << " " << clebsch_de << " " << clebsch_def << "   =  " << clebsch_ab * clebsch_abc * clebsch_de * clebsch_def << " * " << meJ << std::endl;

//        if ((a==0 and b==0 and c==3 and d==2 and e==5 and f==2) or (d==0 and e==0 and f==3 and a==2 and b==5 and c==2)) 
//        if ((a==0 and b==0 and c==3 and d==2 and e==2 and f==5 and Jde==1) ) 
//        if ((a==0 and b==0 and c==3 and d==2 and e==2 and f==5 ) or  (a==0 and b==3 and c==0 and d==2 and e==2 and f==5 ) or (a==3 and b==0 and c==0 and d==2 and e==2 and f==5 )) 
//        if ((a==0 and b==0 and c==5 and d==0 and e==0 and f==5 ) or  (a==0 and b==5 and c==0 and d==0 and e==0 and f==5 ) or (a==5 and b==0 and c==0 and d==0 and e==0 and f==5 )) 
//        if ((a==0 and b==0 and c==10 and d==0 and e==10 and f==0 ) ) 
//        if ((a==2 and b==0 and c==0 and d==2 and e==0 and f==0 ) ) 
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




// The following two routines implicitly assume that we've chosen the apropriate projection quantum numbers.
// This may or may not be the most straightforward way to do things.
double UnitTest::GetMschemeMatrixElement_1leg( const Operator& Op, int a, int ma )
{
  
  return Op.OneBody(a,0) ;
//  return GetMschemeMatrixElement_1b( Op, a, ma, Op.GetQSpaceOrbit(), ma) ;
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
    
//                   Zm_ijkl += ( na*(1-nb) - nb*(1-na) ) * ( Xaibk*Ybjal - Xajbk*Ybial - Xaibl*Ybjak + Xajbl*Ybiak );
                   Zm_ijkl += (na - nb) * ( Xaibk*Ybjal - Xajbk*Ybial - Xaibl*Ybjak + Xajbl*Ybiak );


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
      if (j<i) continue;
      int mj = mi;
      double Zm_ij = 0;
      size_t norb = X.modelspace->GetNumberOrbits();
//      std::cout << __func__ << "  ij = " << i << " " << j << std::endl;
//      for ( auto a : X.modelspace->all_orbits )
      #pragma omp parallel for schedule(dynamic,1) reduction(+:Zm_ij)
      for (size_t a=0; a<norb; a++)
      {
        Orbit& oa = X.modelspace->GetOrbit(a);
        double na = oa.occ;
        if ( na<1e-6 ) continue;
        for ( auto b : X.modelspace->all_orbits )
        {
          Orbit& ob = X.modelspace->GetOrbit(b);
          double nb = ob.occ;
          if ( nb<1e-6 ) continue;
          for ( auto c : X.modelspace->all_orbits )
          {
            Orbit& oc = X.modelspace->GetOrbit(c);
            double nc = oc.occ;
            if ( (1-nc)<1e-6 ) continue;
            for ( auto d : X.modelspace->all_orbits )
            {
              Orbit& od = X.modelspace->GetOrbit(d);
              double nd = od.occ;
              if ( (1-nd)<1e-6 ) continue;
              for ( auto e : X.modelspace->all_orbits )
              {
                Orbit& oe = X.modelspace->GetOrbit(e);
                if ( (oa.l+ob.l+oi.l + oc.l+od.l+oe.l)%2 > 0) continue;
                if ( (oa.tz2+ob.tz2+oi.tz2) != (oc.tz2+od.tz2+oe.tz2)) continue;
                double ne = oe.occ;
                if ( (1-ne)<1e-6 ) continue;
                double occfactor = na*nb*(1-nc)*(1-nd)*(1-ne);
                if (std::abs(occfactor)<1e-7) continue;
//                std::cout << "abcde " << a << " " << b << " " << c << " " << d << " " << e << " " << e << std::endl;

                for (int ma=-oa.j2; ma<=oa.j2; ma+=2)
                {
                 for (int mb=-ob.j2; mb<=ob.j2; mb+=2)
                 {
                  for (int mc=-oc.j2; mc<=oc.j2; mc+=2)
                  {
                   for (int md=-od.j2; md<=od.j2; md+=2)
                   {
                    int me = ma+mb+mi - mc-md;
                    if (std::abs(me) > oe.j2) continue;
// Z_ij = 1/4 sum_abcde (nanb n`c n`dn`e) (X_abicde Y_cdeabj - Y_abicde X_cdeabj)
                    double xabicde = GetMschemeMatrixElement_3b( X, a,ma, b,mb, i,mi, c,mc, d,md, e,me );
                    double ycdeabj = GetMschemeMatrixElement_3b( Y, c,mc, d,md, e,me, a,ma, b,mb, j,mj );
                    double xcdeabj = GetMschemeMatrixElement_3b( X, c,mc, d,md, e,me, a,ma, b,mb, j,mj );
                    double yabicde = GetMschemeMatrixElement_3b( Y, a,ma, b,mb, i,mi, c,mc, d,md, e,me );
                    Zm_ij += 1./12 * occfactor * ( xabicde * ycdeabj - yabicde * xcdeabj );
//                    Zm_ij += 1./4 * occfactor * ( xabicde * ycdeabj - yabicde * xcdeabj );
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
                  << "   zij = " << Zm_ij << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl; 
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

      size_t norb = X.modelspace->GetNumberOrbits();
//      for ( auto a : X.modelspace->all_orbits )
      #pragma omp parallel for schedule(dynamic,1) reduction(+:Zm_ij)
      for ( size_t a=0; a<norb; a++ )
      {
        Orbit& oa = X.modelspace->GetOrbit(a);
        double na = oa.occ;
        if (na<1e-6) continue;
        for ( auto b : X.modelspace->all_orbits )
        {
          Orbit& ob = X.modelspace->GetOrbit(b);
          double nb = ob.occ;
          if (nb<1e-6) continue;
          for ( auto c : X.modelspace->all_orbits )
          {
            Orbit& oc = X.modelspace->GetOrbit(c);
            double nc = oc.occ;
            if ((1-nc)<1e-6) continue;
            for ( auto d : X.modelspace->all_orbits )
            {
              Orbit& od = X.modelspace->GetOrbit(d);
              double nd = od.occ;
              if ((1-nd)<1e-6) continue;
              if ( (oa.l+ob.l+oc.l+od.l)%2 > 0 ) continue;
              if ( (oa.tz2+ob.tz2) !=(oc.tz2+od.tz2) ) continue;
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
                  << "   zij = " << Zm_ij << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl; 
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

//                  if ( mi != mk) continue;
//             if (not (i==0 and j==0 and k==0 and l==10)) continue;
             double Zm_ijkl = 0;
             for (auto a : X.modelspace->all_orbits )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
              for (auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
              {
                Orbit& ob = X.modelspace->GetOrbit(b);
                double nb = ob.occ;


                for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                {
                   int mb = ma;
                
                  double xab = GetMschemeMatrixElement_1b( X, a,ma, b,mb );
                  double yab = GetMschemeMatrixElement_1b( Y, a,ma, b,mb );
                  double xijbkla = GetMschemeMatrixElement_3b( X, i,mi, j,mj, b,mb, k,mk, l,ml, a,ma );
                  double yijbkla = GetMschemeMatrixElement_3b( Y, i,mi, j,mj, b,mb, k,mk, l,ml, a,ma );

                  Zm_ijkl +=  ( na - nb ) * (1*xab * yijbkla - 1*yab * xijbkla   );
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



/// M-Scheme Formula:
//
// Z_ijkl = -1/2 * sum_abc (n_a*n_b*nbar_c + nbar_a*nbar_b*n_c) * [  Xicab*Yabjklc - Xjcab*Yabiklc - Yijcabl*Xabkc + Yijcabk*Xablc 
//                                                                 - Yicab*Xabjklc + Yjcab*Xabiklc + Xijcabl*Yabkc - Xijcabk*Yablc ]
//
bool UnitTest::Test_comm232ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();


  Commutator::comm232ss( X, Y, Z_J);
//  Z_J.Erase();
//  Commutator::comm232ss_debug( X, Y, Z_J);
//  Commutator::comm232ss_slow( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;
  size_t norb = X.modelspace->GetNumberOrbits();
  std::cout << "Begin m-scheme loops" << std::endl;

  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
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
          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue; // check parity
          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue; // check isospin projection

//          for (int mi=-oi.j2; mi<=oi.j2; mi+=2)
          for (int mi= oi.j2; mi<=oi.j2; mi+=2)
          {
          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {
             if ( std::abs(mk)>ol.j2 ) continue; // this just eliminates some redundant matrix elements (I think?)
             int ml = mi + mj - mk;
             if ( std::abs(ml)>ol.j2) continue;

//             if (i>0 or j>0 or k>0 or l>0) continue;
//             if (not (i==0 and j==1 and k==0 and l==1)) continue;
//             if (not (i==0 and j==1 and k==2 and l==5)) continue;
//             if (not (i==0 and j==0 and k==2 and l==2)) continue;
//             if (not (mi==-1 and mj==1 and mk==1 and ml==-1)) continue;
//             if (not (i==2 and j==2 and k==0 and l==0)) continue;
//             if (not (mi==1 and mj==-1 and mk==-1 and ml==1)) continue;

             double Zm_ijkl = 0;

//             for (auto a : X.modelspace->all_orbits )
             #pragma omp parallel for schedule(dynamic,1) reduction(+:Zm_ijkl)
             for (size_t a=0; a<norb; a++ )
             {
              Orbit& oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
              for (auto b : X.modelspace->all_orbits )
              {
                Orbit& ob = X.modelspace->GetOrbit(b);
                double nb = ob.occ;

//                if ( (oi.l+oj.l+ob.l+ok.l+ol.l+oa.l)%2>0) continue;
//                if ( (oi.tz2+oj.tz2+ob.tz2) != (ok.tz2+ol.tz2+oa.tz2) ) continue;

//                     if (b>a) continue ;

                for (auto c : X.modelspace->all_orbits )
                {
                  Orbit& oc = X.modelspace->GetOrbit(c);
                  double nc = oc.occ;

//                  if ( not ((a==0 and b==3) or  (a==3 and b==0)) ) continue;
//                  if ( not ( (a==2 and b==2) or (a==0 and b==0) ) ) continue;
//                  if ( not (a==5 and b==4 and c==1) ) continue;
                  for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                  {
                    for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
                    {
                      for (int mc=-oc.j2; mc<=oc.j2; mc+=2 )
                      {
                  
                        double xicab = GetMschemeMatrixElement_2b( X, i,mi, c,mc, a,ma, b,mb );
                        double yicab = GetMschemeMatrixElement_2b( Y, i,mi, c,mc, a,ma, b,mb );
                        double xjcab = GetMschemeMatrixElement_2b( X, j,mj, c,mc, a,ma, b,mb );
                        double yjcab = GetMschemeMatrixElement_2b( Y, j,mj, c,mc, a,ma, b,mb );
                        double xabkc = GetMschemeMatrixElement_2b( X, a,ma, b,mb, k,mk, c,mc );
                        double yabkc = GetMschemeMatrixElement_2b( Y, a,ma, b,mb, k,mk, c,mc );
                        double xablc = GetMschemeMatrixElement_2b( X, a,ma, b,mb, l,ml, c,mc );
                        double yablc = GetMschemeMatrixElement_2b( Y, a,ma, b,mb, l,ml, c,mc );

                        double xabjklc = GetMschemeMatrixElement_3b( X, a,ma, b,mb, j,mj, k,mk, l,ml, c,mc );
                        double yabjklc = GetMschemeMatrixElement_3b( Y, a,ma, b,mb, j,mj, k,mk, l,ml, c,mc );
                        double xabiklc = GetMschemeMatrixElement_3b( X, a,ma, b,mb, i,mi, k,mk, l,ml, c,mc );
                        double yabiklc = GetMschemeMatrixElement_3b( Y, a,ma, b,mb, i,mi, k,mk, l,ml, c,mc );
                        double xijcabl = GetMschemeMatrixElement_3b( X, i,mi, j,mj, c,mc, a,ma, b,mb, l,ml );
                        double yijcabl = GetMschemeMatrixElement_3b( Y, i,mi, j,mj, c,mc, a,ma, b,mb, l,ml );
                        double xijcabk = GetMschemeMatrixElement_3b( X, i,mi, j,mj, c,mc, a,ma, b,mb, k,mk );
                        double yijcabk = GetMschemeMatrixElement_3b( Y, i,mi, j,mj, c,mc, a,ma, b,mb, k,mk );

                        double occfactor = na*nb*(1-nc) + (1-na)*(1-nb)*nc;

                        Zm_ijkl += -0.5* occfactor * ( xicab*yabjklc - xjcab*yabiklc - yijcabl*xabkc + yijcabk*xablc
                                                     - yicab*xabjklc + yjcab*xabiklc + xijcabl*yabkc - xijcabk*yablc );
//                        Zm_ijkl += -0.5* occfactor * ( xicab*yabjklc - xjcab*yabiklc - yijcabl*xabkc + yijcabk*xablc
//                                                     - yicab*xabjklc + yjcab*xabiklc + xijcabl*yabkc - xijcabk*yablc );

                      }// for mc
                    }// for mb
                  }// for ma
                }// for c
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
          }// for mi
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
// Z_ijkl = 1/6 * sum_abcd (n_a*n_b*n_c*nbar_d - nbar_a*nbar_b*nbar_c*n_d) * [  Xijdabc*Yabckld - Yijdabc*Xabckld  ]
//
bool UnitTest::Test_comm332_ppph_hhhpss( const Operator& X, const Operator& Y ) // test not yet implemented
{
  Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();


  Commutator::comm332_ppph_hhhpss( X, Y, Z_J);


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
          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue; // check parity
          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue; // check isospin projection

          for (int mi= oi.j2; mi<=oi.j2; mi+=2) // notice the lack of minus sign in the initial mi value
          {
          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {
             if ( std::abs(mk)>ol.j2 ) continue; // this just eliminates some redundant matrix elements (I think?)
             int ml = mi + mj - mk;
             if ( std::abs(ml)>ol.j2) continue;


             double Zm_ijkl = 0;
             size_t norb = X.modelspace->GetNumberOrbits();
             #pragma omp parallel for schedule(dynamic,1)  reduction(+:Zm_ijkl)
             for (size_t a=0; a<norb; a++ )
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
                    double nd = od.occ;
                    double occfactor = na*nb*nc*(1-nd) - (1-na)*(1-nb)*(1-nc)*nd;

                    // These make the contribution trivially zero, so I skip them in the name of efficiently
                    // testing the more complicated part. Commenting them out allows to check that the trivial stuff is right.
                    if ( (oi.l+oj.l+od.l + oa.l+ob.l+oc.l)%2 >0 ) continue;
                    if ( (oi.tz2+oj.tz2+od.tz2) != (oa.tz2+ob.tz2+oc.tz2)) continue;
                    if (std::abs(occfactor)<1e-7) continue;

                    for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                    {
                      for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
                      {
                        for (int mc=-oc.j2; mc<=oc.j2; mc+=2 )
                        {
                          for (int md=-od.j2; md<=od.j2; md+=2 )
                          {
                            if ( (mi+mj+md) != (ma+mb+mc) ) continue;

                            double xijdabc = GetMschemeMatrixElement_3b( X, i,mi, j,mj, d,md, a,ma, b,mb, c,mc );
                            double yijdabc = GetMschemeMatrixElement_3b( Y, i,mi, j,mj, d,md, a,ma, b,mb, c,mc );
                            double xabckld = GetMschemeMatrixElement_3b( X, a,ma, b,mb, c,mc, k,mk, l,ml, d,md );
                            double yabckld = GetMschemeMatrixElement_3b( Y, a,ma, b,mb, c,mc, k,mk, l,ml, d,md );

                           Zm_ijkl += 1./6 * occfactor * (xijdabc * yabckld - yijdabc * xabckld );

                          }// for md
                        }// for mc
                      }// for mb
                    }// for ma
                  }// for d
                }// for c
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
          }// for mi
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
// Z_ijkl = 1/4 (1-Pij)(1-Pkl) sum_abcd (n_a*n_b*nbar_c*nbar_d) * (  Xabicdk*Ycdjabl - Yabicdk*Xcdjabl  )
//
// Z_ijkl = 1/4  sum_abcd (n_a*n_b*nbar_c*nbar_d) * (  Xabicdk*Ycdjabl - Yabicdk*Xcdjabl 
//                                                   - Xabjcdk*Ycdiabl + Yabjcdk*Xcdiabl
//                                                   - Xabicdl*Ycdjabk + Yabicdl*Xcdjabk
//                                                   + Xabjcdl*Ycdiabk - Yabjcdl*Xcdiabk )
//
bool UnitTest::Test_comm332_pphhss( const Operator& X, const Operator& Y ) // test not yet implemented
{
  Operator Z_J( Y );
  Operator Z_J_old( Y );
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J_old.SetHermitian();
  Z_J_old.Erase();

//  Z_J.modelspace->SetOccNat3Cut(1e-5);
//  Z_J.modelspace->SetdE3max(1);

//  Commutator::comm332_pphhss_debug( X, Y, Z_J);
  Commutator::comm332_pphhss_debug( X, Y, Z_J_old);
//  Z_J.Erase();
  Commutator::comm332_pphhss( X, Y, Z_J);

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
          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue; // check parity
          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue; // check isospin projection

//          for (int mi=-oi.j2; mi<=oi.j2; mi+=2)
          for (int mi= oi.j2; mi<=oi.j2; mi+=2)
          {
          for (int mj=-oj.j2; mj<=oj.j2; mj+=2)
          {
           for (int mk=-ok.j2; mk<=ok.j2; mk+=2)
           {
             if ( std::abs(mk)>ol.j2 ) continue; // this just eliminates some redundant matrix elements (I think?)
             int ml = mi + mj - mk;
             if ( std::abs(ml)>ol.j2) continue;

//             if (i>0 or j>0 or k>0 or l>0) continue;
//             if (not (i==0 and j==1 and k==0 and l==1)) continue;
//             if (not (i==0 and j==1 and k==2 and l==5)) continue;
//             if (not (i==0 and j==1 and k==2 and l==3)) continue;
//             if (not (i==0 and j==2 and k==0 and l==2)) continue;
//             if (not (mi==-1 and mj==1 and mk==1 and ml==-1)) continue;
//             if (not (i==2 and j==2 and k==0 and l==0)) continue;
//             if (not (mi==1 and mj==-1 and mk==-1 and ml==1)) continue;
//             if (not (mi==1 and mj==1 and mk==1 and ml==1)) continue;

             double Zm_ijkl = 0;
             size_t norb = X.modelspace->GetNumberOrbits();
//             for (auto a : X.modelspace->all_orbits )
             #pragma omp parallel for schedule(dynamic,1)  reduction(+:Zm_ijkl)
             for ( size_t a=0; a<norb; a++ )
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
                    double nd = od.occ;

                    if (  (oa.l+ob.l+oi.l + oc.l+od.l+ok.l)%2>0
                      and (oa.l+ob.l+oi.l + oc.l+od.l+ol.l)%2>0
                      and (oa.l+ob.l+oj.l + oc.l+od.l+ol.l)%2>0
                      and (oa.l+ob.l+oj.l + oc.l+od.l+ok.l)%2>0 ) continue; // check parity for X and Y

                    if (  (oa.tz2+ob.tz2+oi.tz2) != (oc.tz2+od.tz2+ok.tz2)
                      and (oa.tz2+ob.tz2+oi.tz2) != (oc.tz2+od.tz2+ol.tz2)
                      and (oa.tz2+ob.tz2+oj.tz2) != (oc.tz2+od.tz2+ol.tz2)
                      and (oa.tz2+ob.tz2+oj.tz2) != (oc.tz2+od.tz2+ok.tz2) ) continue; // check isospin projection


                    for (int ma=-oa.j2; ma<=oa.j2; ma+=2 )
                    {
                      for (int mb=-ob.j2; mb<=ob.j2; mb+=2 )
                      {
                        for (int mc=-oc.j2; mc<=oc.j2; mc+=2 )
                        {
                          for (int md=-od.j2; md<=od.j2; md+=2 )
                          {
                    // Z_ijkl = 1/4  sum_abcd (n_a*n_b*nbar_c*nbar_d) * (  Xabicdk*Ycdjabl - Yabicdk*Xcdjabl 
//                                                   - Xabjcdk*Ycdiabl + Yabjcdk*Xcdiabl
//                                                   - Xabicdl*Ycdjabk + Yabicdl*Xcdjabk
//                                                   + Xabjcdl*Ycdiabk - Yabjcdl*Xcdiabk )

                          double xabicdk = GetMschemeMatrixElement_3b( X, a,ma, b,mb, i,mi, c,mc, d,md, k,mk );
                          double xabicdl = GetMschemeMatrixElement_3b( X, a,ma, b,mb, i,mi, c,mc, d,md, l,ml );
                          double xabjcdk = GetMschemeMatrixElement_3b( X, a,ma, b,mb, j,mj, c,mc, d,md, k,mk );
                          double xabjcdl = GetMschemeMatrixElement_3b( X, a,ma, b,mb, j,mj, c,mc, d,md, l,ml );

                          double ycdiabk = GetMschemeMatrixElement_3b( Y, c,mc, d,md, i,mi, a,ma, b,mb, k,mk );
                          double ycdiabl = GetMschemeMatrixElement_3b( Y, c,mc, d,md, i,mi, a,ma, b,mb, l,ml );
                          double ycdjabk = GetMschemeMatrixElement_3b( Y, c,mc, d,md, j,mj, a,ma, b,mb, k,mk );
                          double ycdjabl = GetMschemeMatrixElement_3b( Y, c,mc, d,md, j,mj, a,ma, b,mb, l,ml );

                          double occfactor = (1-na)*(1-nb)*nc*nd - na*nb*(1-nc)*(1-nd);

//                          Zm_ijkl += 1./4 * occfactor * (xabicdl*ycdjabk 
                          double dz = 1./4 * occfactor * (xabicdl*ycdjabk 
                                                       -  xabjcdl*ycdiabk 
                                                       -  xabicdk*ycdjabl 
                                                       +  xabjcdk*ycdiabl   );
                          Zm_ijkl += dz;

//                          if (std::abs(dz)>1.0e-6 or (a==0 and b==1 and c==4 and d==5) )
//                          if ( (a==0 and b==1 and c==4 and d==5) )
//                          if (std::abs(dz)>1.0e-6  )
//                          {
//                          std::cout << "abcd " << a << " " << b << " " << c << " " << d << " {m} " << ma << " " << mb << " " << mc << " " << md
//                                    << "    occfactor " << occfactor
////                                    << "  xabjcdk " << xabjcdk << "  ycdiabk " << ycdiabk 
//                                    << "  xabicdl " << xabicdl << "  ycdjabk " << ycdjabk 
////                                    << "  xabicdk " << xabicdk << "  xabicdl " << xabicdl << "  xabjcdk  " << xabjcdk << "  xabjcdl " << xabjcdl
////                                    << "  ycdiabk " << ycdiabk << "  ycdiabl " << ycdiabl << "  ycdjabk  " << ycdjabk << "  ycdjabl " << ycdjabl
//                                    << "  dz = " << dz << "  => " << Zm_ijkl << std::endl;
//                          }

                          }// for md
                        }// for mc
                      }// for mb
                    }// for ma
                  }// for d
                }// for c
              }// for b
             }// for a


             double ZJ_ijkl = GetMschemeMatrixElement_2b( Z_J, i,mi, j,mj, k,mk, l,ml ) ;
             double ZJ_old_ijkl = GetMschemeMatrixElement_2b( Z_J_old, i,mi, j,mj, k,mk, l,ml ) ;
             double err = Zm_ijkl - ZJ_ijkl;
//             if (std::abs(ZJ_ijkl-ZJ_old_ijkl)>1e-6)
             if (std::abs(err)>1e-6)
             {
               std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                         << " {m} = " << mi << " " << mj << " " << mk << " " << ml 
//                         << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl; 
                         << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << " ZJ_old_ijkl " << ZJ_old_ijkl << "   err = " << err
                         << "  J err = " << ZJ_ijkl - ZJ_old_ijkl << std::endl; 
             }
             summed_error += err*err;
             sum_m += Zm_ijkl*Zm_ijkl;
             sum_J += ZJ_ijkl*ZJ_ijkl;

           }// for mk
          }// for mj
          }// for mi
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
  Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();


  Commutator::comm223ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;
  std::cout << "start mscheme loops " << std::endl;

  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      if (j>i) continue;
      Orbit& oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits )
      {
        Orbit& ok = X.modelspace->GetOrbit(k);
        if (k>j) continue;


        for (auto l : X.modelspace->all_orbits )
        {
          if (l>i) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          for ( auto m : X.modelspace->all_orbits )
          {
            if (m>l) continue;
            Orbit& om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits )
            {
              if (n>m) continue;
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
                 if ( (i==j and m_i==m_j) or (i==k and m_i==m_k) or j==k and m_j==m_k) continue;

                 for (int m_l=-ol.j2; m_l<=ol.j2; m_l+=2)
                 {
                  for (int m_m=-om.j2; m_m<=om.j2; m_m+=2)
                  {
                   for (int m_n=-on.j2; m_n<=on.j2; m_n+=2)
                   {
                     if ( (m_i+m_j+m_k) != (m_l+m_m+m_n) ) continue;
                     if ( (l==m and m_l==m_m) or (l==n and m_l==m_n) or m==n and m_m==m_n) continue;
                     double z_ijklmn = 0;


                     size_t norb = X.modelspace->GetNumberOrbits();
//                     for ( auto a : X.modelspace->all_orbits )
                     #pragma omp parallel for schedule(dynamic,1) reduction(+:z_ijklmn)
                     for ( size_t a=0;a<norb; a++ )
                     {
                       Orbit& oa = X.modelspace->GetOrbit(a);
                       for ( int m_a=-oa.j2; m_a<=oa.j2; m_a+=2)
                       {
                         // "direct" term  Z1
                         double x_ijla = GetMschemeMatrixElement_2b( X, i,m_i, j,m_j, l,m_l, a,m_a );
                         double y_akmn = GetMschemeMatrixElement_2b( Y, a,m_a, k,m_k, m,m_m, n,m_n );
                         double y_ijla = GetMschemeMatrixElement_2b( Y, i,m_i, j,m_j, l,m_l, a,m_a );
                         double x_akmn = GetMschemeMatrixElement_2b( X, a,m_a, k,m_k, m,m_m, n,m_n );
                         // Pik  Z2
                         double x_kjla = GetMschemeMatrixElement_2b( X, k,m_k, j,m_j, l,m_l, a,m_a );
                         double y_aimn = GetMschemeMatrixElement_2b( Y, a,m_a, i,m_i, m,m_m, n,m_n );
                         double y_kjla = GetMschemeMatrixElement_2b( Y, k,m_k, j,m_j, l,m_l, a,m_a );
                         double x_aimn = GetMschemeMatrixElement_2b( X, a,m_a, i,m_i, m,m_m, n,m_n );
                         // Pjk   Z3
                         double x_ikla = GetMschemeMatrixElement_2b( X, i,m_i, k,m_k, l,m_l, a,m_a );
                         double y_ajmn = GetMschemeMatrixElement_2b( Y, a,m_a, j,m_j, m,m_m, n,m_n );
                         double y_ikla = GetMschemeMatrixElement_2b( Y, i,m_i, k,m_k, l,m_l, a,m_a );
                         double x_ajmn = GetMschemeMatrixElement_2b( X, a,m_a, j,m_j, m,m_m, n,m_n );
                         // Plm   Z4
                         double x_ijma = GetMschemeMatrixElement_2b( X, i,m_i, j,m_j, m,m_m, a,m_a );
                         double y_akln = GetMschemeMatrixElement_2b( Y, a,m_a, k,m_k, l,m_l, n,m_n );
                         double y_ijma = GetMschemeMatrixElement_2b( Y, i,m_i, j,m_j, m,m_m, a,m_a );
                         double x_akln = GetMschemeMatrixElement_2b( X, a,m_a, k,m_k, l,m_l, n,m_n );
                         // Pln   Z5
                         double x_ijna = GetMschemeMatrixElement_2b( X, i,m_i, j,m_j, n,m_n, a,m_a );
                         double y_akml = GetMschemeMatrixElement_2b( Y, a,m_a, k,m_k, m,m_m, l,m_l );
                         double y_ijna = GetMschemeMatrixElement_2b( Y, i,m_i, j,m_j, n,m_n, a,m_a );
                         double x_akml = GetMschemeMatrixElement_2b( X, a,m_a, k,m_k, m,m_m, l,m_l );
                         // Pik Plm    Z6
                         double x_kjma = GetMschemeMatrixElement_2b( X, k,m_k, j,m_j, m,m_m, a,m_a );
                         double y_ailn = GetMschemeMatrixElement_2b( Y, a,m_a, i,m_i, l,m_l, n,m_n );
                         double y_kjma = GetMschemeMatrixElement_2b( Y, k,m_k, j,m_j, m,m_m, a,m_a );
                         double x_ailn = GetMschemeMatrixElement_2b( X, a,m_a, i,m_i, l,m_l, n,m_n );
                         // Pik Pln    Z7
                         double x_kjna = GetMschemeMatrixElement_2b( X, k,m_k, j,m_j, n,m_n, a,m_a );
                         double y_aiml = GetMschemeMatrixElement_2b( Y, a,m_a, i,m_i, m,m_m, l,m_l );
                         double y_kjna = GetMschemeMatrixElement_2b( Y, k,m_k, j,m_j, n,m_n, a,m_a );
                         double x_aiml = GetMschemeMatrixElement_2b( X, a,m_a, i,m_i, m,m_m, l,m_l );
                         // Pjk Plm    Z8
                         double x_ikma = GetMschemeMatrixElement_2b( X, i,m_i, k,m_k, m,m_m, a,m_a );
                         double y_ajln = GetMschemeMatrixElement_2b( Y, a,m_a, j,m_j, l,m_l, n,m_n );
                         double y_ikma = GetMschemeMatrixElement_2b( Y, i,m_i, k,m_k, m,m_m, a,m_a );
                         double x_ajln = GetMschemeMatrixElement_2b( X, a,m_a, j,m_j, l,m_l, n,m_n );
                         // Pjk Pln    Z9
                         double x_ikna = GetMschemeMatrixElement_2b( X, i,m_i, k,m_k, n,m_n, a,m_a );
                         double y_ajml = GetMschemeMatrixElement_2b( Y, a,m_a, j,m_j, m,m_m, l,m_l );
                         double y_ikna = GetMschemeMatrixElement_2b( Y, i,m_i, k,m_k, n,m_n, a,m_a );
                         double x_ajml = GetMschemeMatrixElement_2b( X, a,m_a, j,m_j, m,m_m, l,m_l );
//

                         z_ijklmn += x_ijla * y_akmn - y_ijla * x_akmn;
                         z_ijklmn -= x_kjla * y_aimn - y_kjla * x_aimn;
                         z_ijklmn -= x_ikla * y_ajmn - y_ikla * x_ajmn;
                         z_ijklmn -= x_ijma * y_akln - y_ijma * x_akln;
                         z_ijklmn -= x_ijna * y_akml - y_ijna * x_akml;
                         z_ijklmn += x_kjma * y_ailn - y_kjma * x_ailn;
                         z_ijklmn += x_kjna * y_aiml - y_kjna * x_aiml;
                         z_ijklmn += x_ikma * y_ajln - y_ikma * x_ajln;
                         z_ijklmn += x_ikna * y_ajml - y_ikna * x_ajml;

                       }// for m_a
                     }// for a
                    
                    double ZJ_ijklmn = GetMschemeMatrixElement_3b( Z_J, i,m_i, j,m_j, k,m_k, l,m_l, m,m_m, n,m_n );
                    double err = z_ijklmn - ZJ_ijklmn;
                    if (std::abs(err)>1e-6 )
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



bool UnitTest::Test_comm133ss( const Operator& X, const Operator& Y )
{

  Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  Commutator::comm133ss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
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

//              if ( not (i==1 and j==0 and k==0 and l==1 and m==0 and n==0) ) continue;

//              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i= oi.j2; m_i<=oi.j2; m_i+=2)
              {
               for (int m_j=-oj.j2; m_j<=oj.j2; m_j+=2)
               {
                for (int m_k=-ok.j2; m_k<=ok.j2; m_k+=2)
                {
                 if ( (i==j and m_i==m_j) or (i==k and m_i==m_k) or j==k and m_j==m_k) continue;

                 for (int m_l=-ol.j2; m_l<=ol.j2; m_l+=2)
                 {
                  for (int m_m=-om.j2; m_m<=om.j2; m_m+=2)
                  {
                   for (int m_n=-on.j2; m_n<=on.j2; m_n+=2)
                   {
                     if ( (m_i+m_j+m_k) != (m_l+m_m+m_n) ) continue;
                     if ( (l==m and m_l==m_m) or (l==n and m_l==m_n) or m==n and m_m==m_n) continue;

//                     if (not (m_i==1 and m_j==1 and m_k==-1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                     double z_ijklmn = 0;

                     size_t norb = X.modelspace->GetNumberOrbits();
//                     for ( auto a : X.modelspace->all_orbits )
                     #pragma omp parallel for schedule(dynamic,1) reduction(+:z_ijklmn)
                     for ( size_t a=0; a<norb; a++ )
                     {
//                       if (a!=11) continue;
                       Orbit& oa = X.modelspace->GetOrbit(a);
                       for ( int m_a=-oa.j2; m_a<=oa.j2; m_a+=2)
                       {
                         double xia = X.OneBody(i,a);
                         double xja = X.OneBody(j,a);
                         double xka = X.OneBody(k,a);
                         double xal = X.OneBody(a,l);
                         double xam = X.OneBody(a,m);
                         double xan = X.OneBody(a,n);

                         double yia = Y.OneBody(i,a);
                         double yja = Y.OneBody(j,a);
                         double yka = Y.OneBody(k,a);
                         double yal = Y.OneBody(a,l);
                         double yam = Y.OneBody(a,m);
                         double yan = Y.OneBody(a,n);

                         double xijalmn = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, a, m_a, l, m_l, m, m_m, n,m_n );
                         double xiaklmn = GetMschemeMatrixElement_3b( X, i,m_i, a,m_a, k, m_k, l, m_l, m, m_m, n,m_n );
                         double xajklmn = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k, m_k, l, m_l, m, m_m, n,m_n );
                         double xijkamn = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, k, m_k, a, m_a, m, m_m, n,m_n );
                         double xijklan = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, k, m_k, l, m_l, a, m_a, n,m_n );
                         double xijklma = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, k, m_k, l, m_l, m, m_m, a,m_a );

                         double yijalmn = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, a, m_a, l, m_l, m, m_m, n,m_n );
                         double yiaklmn = GetMschemeMatrixElement_3b( Y, i,m_i, a,m_a, k, m_k, l, m_l, m, m_m, n,m_n );
                         double yajklmn = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k, m_k, l, m_l, m, m_m, n,m_n );
                         double yijkamn = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, k, m_k, a, m_a, m, m_m, n,m_n );
                         double yijklan = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, k, m_k, l, m_l, a, m_a, n,m_n );
                         double yijklma = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, k, m_k, l, m_l, m, m_m, a,m_a );

                         z_ijklmn += xia * yajklmn  +  xja * yiaklmn  +  xka * yijalmn - yijkamn * xal  -  yijklan * xam  -  yijklma * xan;
                         z_ijklmn -= yia * xajklmn  +  yja * xiaklmn  +  yka * xijalmn - xijkamn * yal  -  xijklan * yam  -  xijklma * yan;
                       }// for m_a
                     }// for a

                    double ZJ_ijklmn = GetMschemeMatrixElement_3b( Z_J, i,m_i, j,m_j, k,m_k, l,m_l, m,m_m, n,m_n );
                    double err = z_ijklmn - ZJ_ijklmn;
                    if (std::abs(err)>1e-6 )
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



// M-Scheme Formula:
//
// Z_ijkl = 1/2 sum_ab  (1-na-nb) [  P(ij/k) ( Xijab Yabklmn - Yijab Xabklmn )  -  P(lm/n) (Xablm Yijkabn - Yablm Xijkabn) ]
//
//       = 1/2 sum_ab (1-na-nb) [ Xijab*Yabklmn - Yijab*Xabklmn - Xkjab*Yabilmn + Ykjab*Xabilmn - Xikab*Yabjlmn + Yikab*Xabjlmn
//                              - Yijkabn*Xablm + Xijkabn*Yablm + Yijkabl*Xabnm - Xijkabl*Yabnm + Yijkabm*Xabln - Xijkabm*Yabln ]
//
bool UnitTest::Test_comm233_pp_hhss( const Operator& X, const Operator& Y ) // test not yet implemented
{

  Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  Commutator::comm233_pp_hhss( X, Y, Z_J);
//  Z_J.Erase();
//  Commutator::comm233_pp_hhss_debug( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;



  size_t norb = X.modelspace->GetNumberOrbits();
  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
//      if (j<i) continue;
      Orbit& oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits )
      {
        Orbit& ok = X.modelspace->GetOrbit(k);
//        if (k<j) continue;


        for (auto l : X.modelspace->all_orbits )
        {
//          if (l<i) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          for ( auto m : X.modelspace->all_orbits )
          {
//            if (m<l) continue;
            Orbit& om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits )
            {
//              if (n<m) continue;
              Orbit& on = X.modelspace->GetOrbit(n);
              if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2 !=0 ) continue;
              if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;

//              if ( not (i==1 and j==0 and k==0 and l==1 and m==0 and n==0) ) continue;
//              if ( not (i==0 and j==0 and k==2 and l==2 and m==2 and n==4) ) continue;
//              if ( not (i==4 and j==2 and k==2 and l==2 and m==0 and n==0) ) continue;

//              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i= oi.j2; m_i<=oi.j2; m_i+=2)
              {
               for (int m_j=-oj.j2; m_j<=oj.j2; m_j+=2)
               {
                for (int m_k=-ok.j2; m_k<=ok.j2; m_k+=2)
                {
                 if ( (i==j and m_i==m_j) or (i==k and m_i==m_k) or j==k and m_j==m_k) continue;

                 for (int m_l=-ol.j2; m_l<=ol.j2; m_l+=2)
                 {
                  for (int m_m=-om.j2; m_m<=om.j2; m_m+=2)
                  {
                   for (int m_n=-on.j2; m_n<=on.j2; m_n+=2)
                   {
                     if ( (m_i+m_j+m_k) != (m_l+m_m+m_n) ) continue;
                     if ( (l==m and m_l==m_m) or (l==n and m_l==m_n) or m==n and m_m==m_n) continue;

//                     if (not (m_i==1 and m_j==-1 and m_k==1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                     double z_ijklmn = 0;

//                     std::cout << " {m} " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n << std::endl;
//                     for ( auto a : X.modelspace->all_orbits )
                     #pragma omp parallel for schedule(dynamic,1) reduction(+:z_ijklmn)
                     for ( size_t a=0; a<norb; a++ )
                     {
                       Orbit& oa = X.modelspace->GetOrbit(a);
                      for ( size_t b=0; b<norb; b++ )
                      {
//                        if (not (( a==4 and b==5) or (a==5 and b==4))  ) continue;
//                        if (not (( a==2 and b==2) or (a==2 and b==2))  ) continue;
                       Orbit& ob = X.modelspace->GetOrbit(b);
                       double occfactor = 1 - oa.occ - ob.occ;
//                       if (a!=11) continue;
                       for ( int m_a=-oa.j2; m_a<=oa.j2; m_a+=2)
                       {
                        for ( int m_b=-ob.j2; m_b<=ob.j2; m_b+=2)
                        {

///       = 1/2 sum_ab (1-na-nb) [ Xijab*Yabklmn - Yijab*Xabklmn - Xkjab*Yabilmn + Ykjab*Xabilmn - Xikab*Yabjlmn + Yikab*Xabjlmn
///                              - Yijkabn*Xablm + Xijkabn*Yablm + Yijkabl*Xabnm - Xijkabl*Yabnm + Yijkabm*Xabln - Xijkabm*Yabln ]
                         double xabilmn = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, i,m_i, l,m_l, m,m_m, n,m_n );
                         double xabjlmn = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, j,m_j, l,m_l, m,m_m, n,m_n );
                         double xabklmn = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, k,m_k, l,m_l, m,m_m, n,m_n );
                         double xijkabl = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, k,m_k, a,m_a, b,m_b, l,m_l );
                         double xijkabm = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, k,m_k, a,m_a, b,m_b, m,m_m );
                         double xijkabn = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, k,m_k, a,m_a, b,m_b, n,m_n );

                         double yabilmn = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, i,m_i, l,m_l, m,m_m, n,m_n );
                         double yabjlmn = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, j,m_j, l,m_l, m,m_m, n,m_n );
                         double yabklmn = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, k,m_k, l,m_l, m,m_m, n,m_n );
                         double yijkabl = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, k,m_k, a,m_a, b,m_b, l,m_l );
                         double yijkabm = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, k,m_k, a,m_a, b,m_b, m,m_m );
                         double yijkabn = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, k,m_k, a,m_a, b,m_b, n,m_n );

                         double xijab = GetMschemeMatrixElement_2b( X, i,m_i, j,m_j, a,m_a, b,m_b  );
                         double xkjab = GetMschemeMatrixElement_2b( X, k,m_k, j,m_j, a,m_a, b,m_b  );
                         double xikab = GetMschemeMatrixElement_2b( X, i,m_i, k,m_k, a,m_a, b,m_b  );
                         double xablm = GetMschemeMatrixElement_2b( X, a,m_a, b,m_b, l,m_l, m,m_m  );
                         double xabnm = GetMschemeMatrixElement_2b( X, a,m_a, b,m_b, n,m_n, m,m_m  );
                         double xabln = GetMschemeMatrixElement_2b( X, a,m_a, b,m_b, l,m_l, n,m_n  );

                         double yijab = GetMschemeMatrixElement_2b( Y, i,m_i, j,m_j, a,m_a, b,m_b  );
                         double ykjab = GetMschemeMatrixElement_2b( Y, k,m_k, j,m_j, a,m_a, b,m_b  );
                         double yikab = GetMschemeMatrixElement_2b( Y, i,m_i, k,m_k, a,m_a, b,m_b  );
                         double yablm = GetMschemeMatrixElement_2b( Y, a,m_a, b,m_b, l,m_l, m,m_m  );
                         double yabnm = GetMschemeMatrixElement_2b( Y, a,m_a, b,m_b, n,m_n, m,m_m  );
                         double yabln = GetMschemeMatrixElement_2b( Y, a,m_a, b,m_b, l,m_l, n,m_n  );


//                         z_ijklmn += 0.5*occfactor*( (xijab*yabklmn - yijab*xabklmn) - (xkjab*yabilmn - ykjab*xabilmn) - (xikab*yabjlmn - yikab*xabjlmn)
//                                                   - (yijkabn*xablm - xijkabn*yablm) + (yijkabl*xabnm - xijkabl*yabnm) + (yijkabm*xabln + xijkabm*yabln) );
                         double dz = 0.5*occfactor*( (1*xijab*yabklmn - 1*yijab*xabklmn) - (1*xkjab*yabilmn - 1*ykjab*xabilmn) - (1*xikab*yabjlmn - 1*yikab*xabjlmn)
                                                   - (1*yijkabn*xablm - 1*xijkabn*yablm) + (1*yijkabl*xabnm - 1*xijkabl*yabnm) + (1*yijkabm*xabln - 1*xijkabm*yabln) );
                         z_ijklmn += dz;
//                         if (std::abs(dz)>1e-6)
//                         {
//                           std::cout << "a,b = " << a << " " << b << "  ma mb = " << m_a << " " << m_b << "  dz = " << dz << "  => " << z_ijklmn << std::endl;
//                           std::cout << "    " << xijab * yabklmn << "  " << yijab*xabklmn << " " << xkjab*yabilmn << " " << ykjab*xabilmn << " " << xikab*yabjlmn << " " << yikab*xabjlmn << std::endl;
//                           std::cout << "    " << yijkabn*xablm << " " << xijkabn*yablm << " " << yijkabl*xabnm << " " << xijkabl*yabnm << " " << yijkabm*xabln << " " << xijkabm*yabln << std::endl;
//                           std::cout << "   [[  " <<  xijab*yabklmn << " " << xkjab*yabilmn << " " <<  xikab*yabjlmn << " ]]  -->  "
//                                     <<   xkjab << " " << yabilmn  << "   ... " << k << " " << j << " " << m_k << " " << m_j << "   " << a << " " << b << " " << j << "  " << m_a << " " << m_b << " " << m_j     << std::endl;
//                         }
                        }// for m_b
                       }// for m_a
                      }// for b
                     }// for a

//                    double zJ_111 = Z_J.ThreeBody.GetME_pn(1,1,1,1,0,0,1,0,0);
//                    std::cout << " zJ_111 = " << zJ_111 << std::endl;
                    double ZJ_ijklmn = GetMschemeMatrixElement_3b( Z_J, i,m_i, j,m_j, k,m_k, l,m_l, m,m_m, n,m_n );
                    double err = z_ijklmn - ZJ_ijklmn;
                    if (std::abs(err)>1e-6 )
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




// M-Scheme Formula:
//
// Z_ijkl = - sum_ab (na-nb) [  P(i/jk)P(l/mn) ( Xbial Yajkbmn - Ybial Xajkbmn )   ]
//
//       =  - sum_ab (na-nb) [ Xbial Yajkbmn - Ybial Xajkbmn  -  Xbjal Yaikbmn + Ybjal Xaikbmn  -  Xbkal Yajibmn + Ybkal Xajibmn
//                           - Xbiam Yajkbln + Ybiam Xajkbln  +  Xbjam Yaikbln - Ybjam Xaikbln  +  Xbkam Yajibln - Ybkam Xajibln
//                           - Xbian Yajkbml + Ybian Xajkbml  +  Xbjan Yaikbml - Ybjan Xaikbml  +  Xbkan Yajibml - Ybkan Xajibml  ]
//                              
//
bool UnitTest::Test_comm233_ph_ss( const Operator& X, const Operator& Y ) // test not yet implemented
{
   Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  Commutator::comm233_phss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      if (j>i) continue;
      Orbit& oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits )
      {
        Orbit& ok = X.modelspace->GetOrbit(k);
        if (k>j) continue;


        for (auto l : X.modelspace->all_orbits )
        {
          if (l>i) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          for ( auto m : X.modelspace->all_orbits )
          {
            if (m>l) continue;
            Orbit& om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits )
            {
              if (n>m) continue;
              Orbit& on = X.modelspace->GetOrbit(n);
              if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2 !=0 ) continue;
              if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;

//              if ( not (i==5 and j==4 and k==0 and l==5 and m==4 and n==0) ) continue;
//              if ( not (i==3 and j==2 and k==0 and l==3 and m==2 and n==0) ) continue;
//              if ( not (i==4 and j==1 and k==0 and l==4 and m==3 and n==2) ) continue;

//              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i=-oi.j2; m_i<=oi.j2; m_i+=2)
              {
               for (int m_j=-oj.j2; m_j<=oj.j2; m_j+=2)
               {
                for (int m_k=-ok.j2; m_k<=ok.j2; m_k+=2)
                {
                 if ( (i==j and m_i==m_j) or (i==k and m_i==m_k) or j==k and m_j==m_k) continue;

                 for (int m_l=-ol.j2; m_l<=ol.j2; m_l+=2)
                 {
                  for (int m_m=-om.j2; m_m<=om.j2; m_m+=2)
                  {
                   for (int m_n=-on.j2; m_n<=on.j2; m_n+=2)
                   {
                     if ( (m_i+m_j+m_k) != (m_l+m_m+m_n) ) continue;
                     if ( (l==m and m_l==m_m) or (l==n and m_l==m_n) or m==n and m_m==m_n) continue;

                     if (not (m_i==1 and m_j==-1 and m_k==1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                     double z_ijklmn = 0;

                     size_t norb = X.modelspace->GetNumberOrbits();
//                     for ( auto a : X.modelspace->all_orbits )
                     #pragma omp parallel for schedule(dynamic,1) reduction(+:z_ijklmn)
                     for ( size_t a=0; a<norb; a++ )
                     {
                       Orbit& oa = X.modelspace->GetOrbit(a);
                      for ( size_t b=0; b<norb; b++ )
                      {
                       Orbit& ob = X.modelspace->GetOrbit(b);
                       double occfactor = oa.occ - ob.occ;
//                       if (a!=11) continue;
//                       if (a!=3 and a!=5) continue;
                       for ( int m_a=-oa.j2; m_a<=oa.j2; m_a+=2)
                       {
                        for ( int m_b=-ob.j2; m_b<=ob.j2; m_b+=2)
                        {

                         double xajkbmn = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, n,m_n );
                         double xaikbmn = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, n,m_n );
                         double xajibmn = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, n,m_n );
                         double xajkbln = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, l,m_l, n,m_n );
                         double xaikbln = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, l,m_l, n,m_n );
                         double xajibln = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, l,m_l, n,m_n );
                         double xajkbml = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, l,m_l );
                         double xaikbml = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, l,m_l );
                         double xajibml = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, l,m_l );

                         double xbial = GetMschemeMatrixElement_2b( X, b,m_b, i,m_i, a,m_a, l,m_l  );
                         double xbjal = GetMschemeMatrixElement_2b( X, b,m_b, j,m_j, a,m_a, l,m_l  );
                         double xbkal = GetMschemeMatrixElement_2b( X, b,m_b, k,m_k, a,m_a, l,m_l  );
                         double xbiam = GetMschemeMatrixElement_2b( X, b,m_b, i,m_i, a,m_a, m,m_m  );
                         double xbjam = GetMschemeMatrixElement_2b( X, b,m_b, j,m_j, a,m_a, m,m_m  );
                         double xbkam = GetMschemeMatrixElement_2b( X, b,m_b, k,m_k, a,m_a, m,m_m  );
                         double xbian = GetMschemeMatrixElement_2b( X, b,m_b, i,m_i, a,m_a, n,m_n  );
                         double xbjan = GetMschemeMatrixElement_2b( X, b,m_b, j,m_j, a,m_a, n,m_n  );
                         double xbkan = GetMschemeMatrixElement_2b( X, b,m_b, k,m_k, a,m_a, n,m_n  );

                         double yajkbmn = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, n,m_n );
                         double yaikbmn = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, n,m_n );
                         double yajibmn = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, n,m_n );
                         double yajkbln = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, l,m_l, n,m_n );
                         double yaikbln = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, l,m_l, n,m_n );
                         double yajibln = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, l,m_l, n,m_n );
                         double yajkbml = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, l,m_l );
                         double yaikbml = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, l,m_l );
                         double yajibml = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, l,m_l );

                         double ybial = GetMschemeMatrixElement_2b( Y, b,m_b, i,m_i, a,m_a, l,m_l  );
                         double ybjal = GetMschemeMatrixElement_2b( Y, b,m_b, j,m_j, a,m_a, l,m_l  );
                         double ybkal = GetMschemeMatrixElement_2b( Y, b,m_b, k,m_k, a,m_a, l,m_l  );
                         double ybiam = GetMschemeMatrixElement_2b( Y, b,m_b, i,m_i, a,m_a, m,m_m  );
                         double ybjam = GetMschemeMatrixElement_2b( Y, b,m_b, j,m_j, a,m_a, m,m_m  );
                         double ybkam = GetMschemeMatrixElement_2b( Y, b,m_b, k,m_k, a,m_a, m,m_m  );
                         double ybian = GetMschemeMatrixElement_2b( Y, b,m_b, i,m_i, a,m_a, n,m_n  );
                         double ybjan = GetMschemeMatrixElement_2b( Y, b,m_b, j,m_j, a,m_a, n,m_n  );
                         double ybkan = GetMschemeMatrixElement_2b( Y, b,m_b, k,m_k, a,m_a, n,m_n  );
//////////////////----------------------------------------------------
                    
                         double xijalmb = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, a,m_a, l,m_l, m,m_m, b,m_b );
                         double yijalmb = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, a,m_a, l,m_l, m,m_m, b,m_b );
                         double xkjalmb = GetMschemeMatrixElement_3b( X, k,m_k, j,m_j, a,m_a, l,m_l, m,m_m, b,m_b );
                         double ykjalmb = GetMschemeMatrixElement_3b( Y, k,m_k, j,m_j, a,m_a, l,m_l, m,m_m, b,m_b );
                         double xijanmb = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, a,m_a, n,m_n, m,m_m, b,m_b );
                         double yijanmb = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, a,m_a, n,m_n, m,m_m, b,m_b );
                         double xikalmb = GetMschemeMatrixElement_3b( X, i,m_i, k,m_k, a,m_a, l,m_l, m,m_m, b,m_b );
                         double yikalmb = GetMschemeMatrixElement_3b( Y, i,m_i, k,m_k, a,m_a, l,m_l, m,m_m, b,m_b );
                         double xijalnb = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, a,m_a, l,m_l, n,m_n, b,m_b );
                         double yijalnb = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, a,m_a, l,m_l, n,m_n, b,m_b );
                         double xkjanmb = GetMschemeMatrixElement_3b( X, k,m_k, j,m_j, a,m_a, n,m_n, m,m_m, b,m_b );
                         double ykjanmb = GetMschemeMatrixElement_3b( Y, k,m_k, j,m_j, a,m_a, n,m_n, m,m_m, b,m_b );
                         double xikanmb = GetMschemeMatrixElement_3b( X, i,m_i, k,m_k, a,m_a, n,m_n, m,m_m, b,m_b );
                         double yikanmb = GetMschemeMatrixElement_3b( Y, i,m_i, k,m_k, a,m_a, n,m_n, m,m_m, b,m_b );
                         double xkjalnb = GetMschemeMatrixElement_3b( X, k,m_k, j,m_j, a,m_a, l,m_l, n,m_n, b,m_b );
                         double ykjalnb = GetMschemeMatrixElement_3b( Y, k,m_k, j,m_j, a,m_a, l,m_l, n,m_n, b,m_b );
                         double xikalnb = GetMschemeMatrixElement_3b( X, i,m_i, k,m_k, a,m_a, l,m_l, n,m_n, b,m_b );
                         double yikalnb = GetMschemeMatrixElement_3b( Y, i,m_i, k,m_k, a,m_a, l,m_l, n,m_n, b,m_b );

//  Zijkl =  - sum_ab (na-nb) [ Xbial Yajkbmn - Ybial Xajkbmn  -  Xbjal Yaikbmn + Ybjal Xaikbmn  -  Xbkal Yajibmn + Ybkal Xajibmn
//                            - Xbiam Yajkbln + Ybiam Xajkbln  +  Xbjam Yaikbln - Ybjam Xaikbln  +  Xbkam Yajibln - Ybkam Xajibln
//                            - Xbian Yajkbml + Ybian Xajkbml  +  Xbjan Yaikbml - Ybjan Xaikbml  +  Xbkan Yajibml - Ybkan Xajibml  ]
                          // xijalmn = xaijblm = xajibml
//                          z_ijklmn -= occfactor * (xbkan * yajibml - ybkan * xajibml);
//                           z_ijklmn -= occfactor * (xbkan * yijalmb - ybkan*xijalmb);
                           double dz = -occfactor *1* (xbkan * yijalmb - ybkan * xijalmb);
                                  dz += occfactor *1* (xbian * ykjalmb - ybian * xkjalmb);
                                  dz += occfactor *1* (xbjan * yikalmb - ybjan * xikalmb);
                                  dz += occfactor *1* (xbkal * yijanmb - ybkal * xijanmb);
                                  dz += occfactor *1* (xbkam * yijalnb - ybkam * xijalnb);
                                  dz -= occfactor *1* (xbial * ykjanmb - ybial * xkjanmb);
                                  dz -= occfactor *1* (xbjal * yikanmb - ybjal * xikanmb);
                                  dz -= occfactor *1* (xbiam * ykjalnb - ybiam * xkjalnb);
                                  dz -= occfactor *1* (xbjam * yikalnb - ybjam * xikalnb);
                           z_ijklmn += dz;

                        }// for m_b
                       }// for m_a
                      }// for b
                     }// for a

                    double ZJ_ijklmn = GetMschemeMatrixElement_3b( Z_J, i,m_i, j,m_j, k,m_k, l,m_l, m,m_m, n,m_n );
                    double err = z_ijklmn - ZJ_ijklmn;
                    if (std::abs(err)>1e-6 )
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



// M-Scheme Formula:
//
// Z_ijkl = 1/6 sum_abc [na*nb*nc + (1-na)*(1-nb)(1-nc) ] * [  Xijkabc*Yabclmn - Yijkabc*Xabclmn  ]
//                              
//
bool UnitTest::Test_comm333_ppp_hhh_ss( const Operator& X, const Operator& Y ) // test not yet implemented
{
   Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  Commutator::comm333_ppp_hhhss( X, Y, Z_J);

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
//    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
//      if (j<i) continue;
      Orbit& oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits )
      {
        Orbit& ok = X.modelspace->GetOrbit(k);
//        if (k<j) continue;


        for (auto l : X.modelspace->all_orbits )
        {
//          if (l<i) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          for ( auto m : X.modelspace->all_orbits )
          {
//            if (m<l) continue;
            Orbit& om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits )
            {
              if (n<m) continue;
              Orbit& on = X.modelspace->GetOrbit(n);
              if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2 !=0 ) continue;
              if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;

//              if ( not (i==1 and j==0 and k==0 and l==1 and m==0 and n==0) ) continue;

//              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i= oi.j2; m_i<=oi.j2; m_i+=2)
              {
               for (int m_j=-oj.j2; m_j<=oj.j2; m_j+=2)
               {
                for (int m_k=-ok.j2; m_k<=ok.j2; m_k+=2)
                {
                 if ( (i==j and m_i==m_j) or (i==k and m_i==m_k) or j==k and m_j==m_k) continue;

                 for (int m_l=-ol.j2; m_l<=ol.j2; m_l+=2)
                 {
                  for (int m_m=-om.j2; m_m<=om.j2; m_m+=2)
                  {
                   for (int m_n=-on.j2; m_n<=on.j2; m_n+=2)
                   {
                     if ( (m_i+m_j+m_k) != (m_l+m_m+m_n) ) continue;
                     if ( (l==m and m_l==m_m) or (l==n and m_l==m_n) or m==n and m_m==m_n) continue;

//                     if (not (m_i==1 and m_j==1 and m_k==-1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                     double z_ijklmn = 0;

                     size_t norb = X.modelspace->GetNumberOrbits();
//                     for ( auto a : X.modelspace->all_orbits )
                     #pragma omp parallel for schedule(dynamic,1) reduction(+:z_ijklmn)
                     for ( size_t a=0; a<norb; a++ )
                     {
                       Orbit& oa = X.modelspace->GetOrbit(a);
                      for ( size_t b=0; b<norb; b++ )
                      {
                       Orbit& ob = X.modelspace->GetOrbit(b);
                      for ( size_t c=0; c<norb; c++ )
                      {
                       Orbit& oc = X.modelspace->GetOrbit(c);
                       double occfactor = oa.occ*ob.occ*oc.occ - (1-oa.occ)*(1-ob.occ)*(1-oc.occ);
//                       if (a!=11) continue;
                       for ( int m_a=-oa.j2; m_a<=oa.j2; m_a+=2)
                       {
                        for ( int m_b=-ob.j2; m_b<=ob.j2; m_b+=2)
                        {
                         for ( int m_c=-oc.j2; m_c<=oc.j2; m_c+=2)
                         {

// Z_ijkl = 1/6 sum_abc [na*nb*nc + (1-na)*(1-nb)(1-nc) ] * [  Xijkabc*Yabclmn - Yijkabc*Xabclmn  ]
                         double xijkabc = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, k,m_k, a,m_a, b,m_b, c,m_c );
                         double xabclmn = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, c,m_c, l,m_l, m,m_m, n,m_n );
                         double yijkabc = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, k,m_k, a,m_a, b,m_b, c,m_c );
                         double yabclmn = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, c,m_c, l,m_l, m,m_m, n,m_n );

                         z_ijklmn += 1./6 * occfactor * ( xijkabc*yabclmn - yijkabc*xabclmn);

                         }// for m_c
                        }// for m_b
                       }// for m_a
                       }// for c
                      }// for b
                     }// for a

                    double ZJ_ijklmn = GetMschemeMatrixElement_3b( Z_J, i,m_i, j,m_j, k,m_k, l,m_l, m,m_m, n,m_n );
                    double err = z_ijklmn - ZJ_ijklmn;
                    if (std::abs(err)>1e-6 )
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


// M-Scheme Formula:
//
// Z_ijkl = - 1/2 sum_abc [na*nb*(1-nc) + (1-na)*(1-nb)*nc ] * P(ij/k) P(l/mn) * [ Xabkcmn*Ycijabl - Yabkcmn*Xcijabl   ]
//                              
//        = - 1/2 sum_abc [na*nb*(1-nc) + (1-na)*(1-nb)*nc ] *
//                [ (Xabkcmn*Ycijabl - Yabkcmn*Xcijabl) - (Xabicmn*Yckjabl - Yabicmn*Xckjabl)  - (Xabjcmn*Ycikabl - Yabjcmn*Xcikabl) 
//                 -(Xabkcln*Ycijabm - Yabkcln*Xcijabm) + (Xabicln*Yckjabm - Yabicln*Xckjabm)  + (Xabjcln*Ycikabm - Yabjcln*Xcikabm) 
//                 -(Xabkcml*Ycijabn - Yabkcml*Xcijabn) + (Xabicml*Yckjabn - Yabicml*Xckjabn)  + (Xabjcml*Ycikabn - Yabjcml*Xcikabn) ]
//
bool UnitTest::Test_comm333_pph_hhp_ss( const Operator& X, const Operator& Y ) // test not yet implemented
{
   Operator Z_J( Y );
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

//  Commutator::comm333ss( X, Y, Z_J); // TODO: make a sub-call to just the pph hhp part
  Commutator::comm333_pph_hhpss( X, Y, Z_J); // TODO: make a sub-call to just the pph hhp part

  if ( Z_J.IsHermitian() )
     Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian() )
     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits )
  {
    Orbit& oi = X.modelspace->GetOrbit(i);
//    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits )
    {
      if (j<i) continue;
//      if (j>i) continue;
      Orbit& oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits )
      {
        Orbit& ok = X.modelspace->GetOrbit(k);
        if (k<j) continue;
//        if (k>j) continue;


        for (auto l : X.modelspace->all_orbits )
        {
          if (l>i) continue;
          Orbit& ol = X.modelspace->GetOrbit(l);
          for ( auto m : X.modelspace->all_orbits )
          {
//            if (m<l) continue;
//            if (m>l) continue;
            Orbit& om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits )
            {
//              if (n<m) continue; 
//              if (n>m) continue;
              Orbit& on = X.modelspace->GetOrbit(n);
              if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2 !=0 ) continue;
              if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;

//              if ( not (i==4 and j==0 and k==0 and l==4 and m==2 and n==2) ) continue;
//              if ( not (i==4 and j==4 and k==5 and l==2 and m==4 and n==5) ) continue;

              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
//              double Z001 = Z_J.ThreeBody.GetME_pn(0,0,1,i,j,k,l,m,n);
//              double Z011 = Z_J.ThreeBody.GetME_pn(0,1,1,i,j,k,l,m,n);
//              double Z101 = Z_J.ThreeBody.GetME_pn(1,0,1,i,j,k,l,m,n);
//              double Z111 = Z_J.ThreeBody.GetME_pn(1,1,1,i,j,k,l,m,n);
//              double Z113 = Z_J.ThreeBody.GetME_pn(1,1,3,i,j,k,l,m,n);
//              std::cout << "Z001 = " <<Z001 << " Z011 = " << Z011 << "  Z101 = " << Z101 << "  Z111 = " << Z111 << "  Z113 = " << Z113 << std::endl; 
              // loop over projections
              for (int m_i= oi.j2; m_i<=oi.j2; m_i+=2)
              {
               for (int m_j=-oj.j2; m_j<=oj.j2; m_j+=2)
               {
                for (int m_k=-ok.j2; m_k<=ok.j2; m_k+=2)
                {
                 if ( (i==j and m_i==m_j) or (i==k and m_i==m_k) or j==k and m_j==m_k) continue;

                 for (int m_l=-ol.j2; m_l<=ol.j2; m_l+=2)
                 {
                  for (int m_m=-om.j2; m_m<=om.j2; m_m+=2)
                  {
                   for (int m_n=-on.j2; m_n<=on.j2; m_n+=2)
                   {
                     if ( (m_i+m_j+m_k) != (m_l+m_m+m_n) ) continue;
                     if ( (l==m and m_l==m_m) or (l==n and m_l==m_n) or (m==n and m_m==m_n)) continue;

                     if (not (m_i==1 and m_j==-1 and m_k==1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                     double z_ijklmn = 0;

                     size_t norb = X.modelspace->GetNumberOrbits();
//                     for ( auto a : X.modelspace->all_orbits )
                     #pragma omp parallel for schedule(dynamic,1) reduction(+:z_ijklmn)
                     for ( size_t a=0; a<norb; a++ )
                     {
                       Orbit& oa = X.modelspace->GetOrbit(a);
                      for ( size_t b=0; b<norb; b++ )
                      {
                       Orbit& ob = X.modelspace->GetOrbit(b);
                      for ( size_t c=0; c<norb; c++ )
                      {
//                       if (not (a==0 and b==0 and c==2)) continue;
//                       if (not (((a==0 and b==1) or (a==1 and b==0))and c==2 )) continue;
//                       if (not (((a==0 and b==1))and c==2 )) continue;
//                       if (not (((a==0 and b==1))and c==4 )) continue;
//                       if (not (((a==0 and b==1) or a==1 and b==0)and c==2 )) continue;
//                       if (not (((a==0 and b==1) or a==1 and b==0)and c<9992 )) continue;
//                       if (not (((a==0 and b==1))and c<1002 )) continue;
//                       if ( (((a==0 and b==1) or (a==1 and b==0))and c<1002 )) continue;
//                       if ( (((a>1 or b>1) )and c<1002 )) continue;
//              if (a!=b) continue;
                       Orbit& oc = X.modelspace->GetOrbit(c);
                       double occfactor = oa.occ * ob.occ * (1-oc.occ) + (1-oa.occ)*(1-ob.occ)*oc.occ ;
                       for ( int m_a=-oa.j2; m_a<=oa.j2; m_a+=2)
                       {
                        for ( int m_b=-ob.j2; m_b<=ob.j2; m_b+=2)
                        {
                         for ( int m_c=-oc.j2; m_c<=oc.j2; m_c+=2)
                         {


                         double xabklmc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, k,m_k, l,m_l, m,m_m, c,m_c );// direct
                         double xijcabn = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, c,m_c, a,m_a, b,m_b, n,m_n );// direct

                         double xabilmc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, i,m_i, l,m_l, m,m_m, c,m_c );// Pik
                         double xkjcabn = GetMschemeMatrixElement_3b( X, k,m_k, j,m_j, c,m_c, a,m_a, b,m_b, n,m_n );// Pik
                         double xabjlmc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, j,m_j, l,m_l, m,m_m, c,m_c );// Pjk
                         double xikcabn = GetMschemeMatrixElement_3b( X, i,m_i, k,m_k, c,m_c, a,m_a, b,m_b, n,m_n );// Pjk
                         double xabknmc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, k,m_k, n,m_n, m,m_m, c,m_c );// Pln
                         double xijcabl = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, c,m_c, a,m_a, b,m_b, l,m_l );// Pln
                         double xabklnc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, k,m_k, l,m_l, n,m_n, c,m_c );// Pmn
                         double xijcabm = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, c,m_c, a,m_a, b,m_b, m,m_m );// Pmn

                         double xabinmc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, i,m_i, n,m_n, m,m_m, c,m_c );// Pik Pln
                         double xkjcabl = GetMschemeMatrixElement_3b( X, k,m_k, j,m_j, c,m_c, a,m_a, b,m_b, l,m_l );// Pik Pln
                         double xabjnmc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, j,m_j, n,m_n, m,m_m, c,m_c );// Pjk Pln
                         double xikcabl = GetMschemeMatrixElement_3b( X, i,m_i, k,m_k, c,m_c, a,m_a, b,m_b, l,m_l );// Pjk Pln
                         double xabilnc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, i,m_i, l,m_l, n,m_n, c,m_c );// Pik Pmn
                         double xkjcabm = GetMschemeMatrixElement_3b( X, k,m_k, j,m_j, c,m_c, a,m_a, b,m_b, m,m_m );// Pik Pmn
                         double xabjlnc = GetMschemeMatrixElement_3b( X, a,m_a, b,m_b, j,m_j, l,m_l, n,m_n, c,m_c );// Pjk Pmn
                         double xikcabm = GetMschemeMatrixElement_3b( X, i,m_i, k,m_k, c,m_c, a,m_a, b,m_b, m,m_m );// Pjk Pmn

                         double yabklmc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, k,m_k, l,m_l, m,m_m, c,m_c );// direct
                         double yijcabn = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, c,m_c, a,m_a, b,m_b, n,m_n );// direct

                         double yabilmc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, i,m_i, l,m_l, m,m_m, c,m_c );// Pik
                         double ykjcabn = GetMschemeMatrixElement_3b( Y, k,m_k, j,m_j, c,m_c, a,m_a, b,m_b, n,m_n );// Pik
                         double yabjlmc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, j,m_j, l,m_l, m,m_m, c,m_c );// Pjk
                         double yikcabn = GetMschemeMatrixElement_3b( Y, i,m_i, k,m_k, c,m_c, a,m_a, b,m_b, n,m_n );// Pjk
                         double yabknmc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, k,m_k, n,m_n, m,m_m, c,m_c );// Pln
                         double yijcabl = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, c,m_c, a,m_a, b,m_b, l,m_l );// Pln
                         double yabklnc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, k,m_k, l,m_l, n,m_n, c,m_c );// Pmn
                         double yijcabm = GetMschemeMatrixElement_3b( Y, i,m_i, j,m_j, c,m_c, a,m_a, b,m_b, m,m_m );// Pmn

                         double yabinmc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, i,m_i, n,m_n, m,m_m, c,m_c );// Pik Pln
                         double ykjcabl = GetMschemeMatrixElement_3b( Y, k,m_k, j,m_j, c,m_c, a,m_a, b,m_b, l,m_l );// Pik Pln
                         double yabjnmc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, j,m_j, n,m_n, m,m_m, c,m_c );// Pjk Pln
                         double yikcabl = GetMschemeMatrixElement_3b( Y, i,m_i, k,m_k, c,m_c, a,m_a, b,m_b, l,m_l );// Pjk Pln
                         double yabilnc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, i,m_i, l,m_l, n,m_n, c,m_c );// Pik Pmn
                         double ykjcabm = GetMschemeMatrixElement_3b( Y, k,m_k, j,m_j, c,m_c, a,m_a, b,m_b, m,m_m );// Pik Pmn
                         double yabjlnc = GetMschemeMatrixElement_3b( Y, a,m_a, b,m_b, j,m_j, l,m_l, n,m_n, c,m_c );// Pjk Pmn
                         double yikcabm = GetMschemeMatrixElement_3b( Y, i,m_i, k,m_k, c,m_c, a,m_a, b,m_b, m,m_m );// Pjk Pmn


//        =   1/2 sum_abc [na*nb*(1-nc) + (1-na)*(1-nb)*nc ] *
//                [ (Xabkcmn*Ycijabl - Yabkcmn*Xcijabl) - (Xabicmn*Yckjabl - Yabicmn*Xckjabl)  - (Xabjcmn*Ycikabl - Yabjcmn*Xcikabl) 
//                 -(Xabkcln*Ycijabm - Yabkcln*Xcijabm) + (Xabicln*Yckjabm - Yabicln*Xckjabm)  + (Xabjcln*Ycikabm - Yabjcln*Xcikabm) 
//                 -(Xabkcml*Ycijabn - Yabkcml*Xcijabn) + (Xabicml*Yckjabn - Yabicml*Xckjabn)  + (Xabjcml*Ycikabn - Yabjcml*Xcikabn) ]
                           double dz = 0;
                           dz += xabklmc * yijcabn - yabklmc * xijcabn ; // Z1
                           dz -= xabilmc * ykjcabn - yabilmc * xkjcabn ; // Z2
                           dz -= xabjlmc * yikcabn - yabjlmc * xikcabn ; // Z3
                           dz -= xabknmc * yijcabl - yabknmc * xijcabl ; // Z4
                           dz -= xabklnc * yijcabm - yabklnc * xijcabm ; // Z5
                           dz += xabinmc * ykjcabl - yabinmc * xkjcabl ; // Z6
                           dz += xabjnmc * yikcabl - yabjnmc * xikcabl ; // Z7
                           dz += xabilnc * ykjcabm - yabilnc * xkjcabm ; // Z8
                           dz += xabjlnc * yikcabm - yabjlnc * xikcabm ; // Z9
                           z_ijklmn += 0.5 * occfactor * dz;
//                           z_ijklmn += 0.5 * occfactor * dz; // extra for ab<=>ba
//                           if (std::abs(dz)>1e-6) std::cout << "  abc  " << a << " " << b << " "  << c << "  {m}  " << m_a << " " << m_b << " " << m_c
//                              << "    x " << xabinmc << "  y " << ykjcabl << "    " << dz << "   " << 0.5*occfactor*dz <<  "    =>  " << z_ijklmn << std::endl;

                          }// for m_c
                         }// for m_b
                        }// for m_a
                       }// for c
                      }// for b
                     }// for a

//                    double X011 = 1*X.ThreeBody.GetME_pn(0,1,1,0,1,4,5,4,2);
//                    double X111 = 1*X.ThreeBody.GetME_pn(1,1,1,0,1,4,5,4,2);
//                    double X113 = 1*X.ThreeBody.GetME_pn(1,1,3,0,1,4,5,4,2);
//                    double X103 = 1*X.ThreeBody.GetME_pn(1,0,3,0,1,4,5,4,2);
//
//                    double Y003 = 1*Y.ThreeBody.GetME_pn(0,0,3,5,4,2,0,1,2);
//                    double Y103 = 1*Y.ThreeBody.GetME_pn(1,0,3,5,4,2,0,1,2);
//                    double Y013 = 1*Y.ThreeBody.GetME_pn(0,1,3,5,4,2,0,1,2);
//                    double Y111 = 1*Y.ThreeBody.GetME_pn(1,1,1,5,4,2,0,1,2);
//                    double Y113 = 1*Y.ThreeBody.GetME_pn(1,1,3,5,4,2,0,1,2);
//                    double Y115 = 1*Y.ThreeBody.GetME_pn(1,1,5,5,4,2,0,1,2);
//
//                    double XA = -1/sqrt(12)*X011 + 1./6*X111 - 1./(3*sqrt(10))*X113 + sqrt(1./6)*X103;
//                    double XB = +1/sqrt(12)*X011 + 1./6*X111 - 1./(3*sqrt(10))*X113 + sqrt(1./6)*X103;
//                    double XC = +1/3.*X111 + 1./(3*sqrt(10))*X113 + 1./(sqrt(6))*X103;
//                    double XD = -sqrt(3./10)*X113 + 1./sqrt(2)*X103;
//
//                    double YA =  1./2*Y003 -1./(2*sqrt(15))*Y103 - 1./(2*sqrt(15))*Y013 + 1./6*Y111 + 1./30*Y113 + 3./10*Y115;
//                    double YB = -1./2*Y003 +1./(2*sqrt(15))*Y103 - 1./(2*sqrt(15))*Y013 + 1./6*Y111 + 1./30*Y113 + 3./10*Y115;
//                    double YC =  -2./(sqrt(15))*Y013  -1./6*Y111  - 2./15*Y113 + 3./10*Y115;
//                    double YD =  1./(sqrt(5))*Y013 - sqrt(3./25)*Y113 + sqrt(3./25)*Y115;
//
//                    double byhand = XA*YA + XB*YB + XC*YC + XD*YD;
//                    std::cout << "A " << XA << " " << YA << "   " << XA*YA << std::endl;
//                    std::cout << "B " << XB << " " << YB << "   " << XB*YB << std::endl;
//                    std::cout << "C " << XC << " " << YC << "   " << XC*YC << std::endl;
//                    std::cout << "D " << XD << " " << YD << "   " << XD*YD << std::endl;
//                    std::cout << "I think by hand is " << byhand << std::endl;
//
//                    double Z0013 = 0.5 * X011*Y003;
//                    double Z0113 = 0.5 * sqrt(5)/3*X111*Y013;
//                    double Z0133 = 0.5 * sqrt(8.)/3 * X113 * Y013;
//
//                    double Z1013 = -sqrt(3)/2 * sqrt(5)/3 * X011 * Y103;
//                    double Z1111 = -sqrt(3)/2 * -1./9 * X111 * Y111;
//                    double Z1113 = -sqrt(3)/2 * +1./9 * X111 * Y113;
//                    double Z1115 = -sqrt(3)/2 * X111 * Y115;
//                    double Z1131 = -sqrt(3)/2 * -sqrt(10.)/9 * X113 * Y111;
//                    double Z1133 = -sqrt(3)/2 * 7*4/(9*sqrt(10)) * X113 * Y113;
//                    double Z1135 = -sqrt(3)/2 * -2/(sqrt(10)) * X113 * Y115;
//                   
//                    std::cout << "I think X113 = " << X113 << "  and  Y111 = " << Y111 << std::endl;
//                    std::cout << "I think the J1p=0 piece should be " << Z0013 << " + " << Z0113 << " + " << Z0133 << "  = " << Z0013+Z0113+Z0133 << std::endl;
//                    std::cout << "and I think the J1p=1 piece should be "
//                    << Z1013 << " + " << Z1111 << " + " << Z1113 << " + " << Z1115 << " + " << Z1131 << " + " << Z1133 << " + " << Z1135 << " = "
//                    << Z1013  +  Z1111  +  Z1113  + Z1115  +  Z1131  +  Z1133 +  Z1135 << std::endl;
//                    std::cout << " so total Jscheme should be " 
//                    << Z0013+Z0113+Z0133 + Z1013  +  Z1111  +  Z1113  + Z1115  +  Z1131  +  Z1133 +  Z1135
//                    << "   =>  " << -2/sqrt(12) * (Z0013+Z0113+Z0133 + Z1013  +  Z1111  +  Z1113  + Z1115  +  Z1131  +  Z1133 +  Z1135) << std::endl;
                    

                    double ZJ_ijklmn = GetMschemeMatrixElement_3b( Z_J, i,m_i, j,m_j, k,m_k, l,m_l, m,m_m, n,m_n );
                    double ZJ_hermconj = GetMschemeMatrixElement_3b( Z_J, l,m_l, m,m_m, n,m_n , i,m_i, j,m_j, k,m_k); 
                    double err = z_ijklmn - ZJ_ijklmn;
//                    std::cout << z_ijklmn << std::endl;
                    if (std::abs(err)>1e-6 )
                    {
                      std::cout << "Trouble in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err <<  "  hc = " << ZJ_hermconj << std::endl; 
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
    double ZJ_i = Z_J.OneBody(i,0);
//    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err)>1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i <<  "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i,0) << "   err = " << err << std::endl; 
//                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl; 
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
    double ZJ_i = Z_J.OneBody(i,0);
//    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err)>1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i <<  "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i,0) << "   err = " << err << std::endl; 
//                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl; 
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
    double ZJ_i = Z_J.OneBody(i,0);
//    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err)>1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i <<  "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i,0) << "   err = " << err << std::endl; 
//                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl; 
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






