
#include "UnitTest.hh"
#include "AngMom.hh"
#include <armadillo>
#include <random>
#include <string>
#include "Commutator.hh"
#include "IMSRG3Commutators.hh"
#include "BCH.hh"
#include "FactorizedDoubleCommutator.hh"
#include "imsrg_util.hh"
#include "version.hh"
#include "ReferenceImplementations.hh" // for commutators
#include "IMSRGSolver.hh"              // for Perturbative Triples

#include <omp.h>

uint64_t UnitTest::random_seed = 1;

UnitTest::UnitTest(ModelSpace &ms)
    : modelspace(&ms)
{
}

Operator UnitTest::RandomOp(ModelSpace &modelspace, int jrank, int tz, int parity, int particle_rank, int hermitian)
{
  Operator Rando(modelspace, jrank, tz, parity, particle_rank);
  if (hermitian == -1)
    Rando.SetAntiHermitian();
  else if (hermitian == +1)
    Rando.SetHermitian();
  else
    Rando.SetNonHermitian();

  random_seed++;

  arma::arma_rng::set_seed(random_seed);

  arma::mat symmetry_allowed = arma::zeros(arma::size(Rando.OneBody));
  for (auto i : modelspace.all_orbits)
  {
    Orbit &oi = modelspace.GetOrbit(i);
    for (auto j : Rando.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
    {
      symmetry_allowed(i, j) = 1;
      if (Rando.IsReduced() and j > i) // if jrank>0, we're storing reduced matrix elements, and we get a phase from Wigner-Eckart.
      {
        Orbit &oj = modelspace.GetOrbit(j);
        symmetry_allowed(i, j) *= AngMom::phase((oi.j2 - oj.j2) / 2);
      }
    }
  }
  Rando.OneBody.randn();
  Rando.OneBody += hermitian * Rando.OneBody.t();
  Rando.OneBody = Rando.OneBody % symmetry_allowed;

  if (particle_rank > 1)
  {
    for (auto &itmat : Rando.TwoBody.MatEl)
    {
      itmat.second.randn();
      if (itmat.first[0] == itmat.first[1])
      {
        itmat.second += hermitian * itmat.second.t();
      }
    }
  }

  // If needed, fill the 3-body piece. Since the storage scheme has
  // places for matrix elements that are zero by symmetry, we need
  // to actually loop through and make sure we don't give a finite value
  // to something that should be zero.
  // In addition, because we use antisymmetrized matrix elements |abc Jab J>,
  // if a=b, we require Jab is even. This requirement becomes more complicated
  // when b=c. It's easier to randomly fill the matrix elements involving |bba Jbb J>
  // and then obtain the |abb Jab J> by recoupling.
  if (particle_rank > 2)
  {
    Rando.ThreeBody.SwitchToPN_and_discard();
    std::default_random_engine generator(random_seed);
    double mean = 0;
    double stddev = 1;
    std::normal_distribution<double> distribution(mean, stddev);
    int herm = Rando.IsHermitian() ? 1 : -1;

    size_t bufsize = 80;
    std::vector<double> randbuf(bufsize);
    //    for (size_t i=0; i<bufsize; i++)   randbuf[i] = distribution(generator);

    for (auto &iter_braket : Rando.ThreeBody.Get_ch_start())
    {
      size_t ch_bra = iter_braket.first.ch_bra;
      size_t ch_ket = iter_braket.first.ch_ket;
      ThreeBodyChannel &Tbc_bra = modelspace.GetThreeBodyChannel(ch_bra);
      ThreeBodyChannel &Tbc_ket = modelspace.GetThreeBodyChannel(ch_ket);
      int twoj1 = Tbc_bra.twoJ;
      int twoj2 = Tbc_ket.twoJ;
      size_t nbras = Tbc_bra.GetNumberKets();

      for (size_t i = 0; i < bufsize; i++)
        randbuf[i] = distribution(generator);

      for (size_t ibra = 0; ibra < nbras; ibra++)
      {
        Ket3 &bra = Tbc_bra.GetKet(ibra);
        // We store ibra <= iket
        size_t nkets = Tbc_ket.GetNumberKets();
        size_t iket_min = ch_bra == ch_ket ? ibra : 0;
        if (ch_bra == ch_ket and herm == -1)
          iket_min = ibra + 1;
        for (size_t iket = iket_min; iket < nkets; iket++)
        {
          Ket3 &ket = Tbc_ket.GetKet(iket);

          // When we have repeated orbits, different Jpq
          if (herm == -1 and ch_bra == ch_ket)
          {
            if ((bra.p == ket.p) and (bra.q == ket.q) and (bra.r == ket.r)) //  < abc | abc >  so possible zero from antihermiticity
            {
              if (bra.q == bra.r)
                continue; // < abb | abb >  in this case, states with different Jab are not orthogonal, so this makes antihermiticty messy. just skip it
            }
          }

          double random_me = 0;

          // If we have |abb> coupled to Jab, then there are linear dependencies due to antisymmetry. This is more transparent if we recouple and enforce Jbb = even.
          if ((bra.q == bra.r) or (ket.q == ket.r))
          {

            std::vector<ThreeBodyStorage::Permutation> perm_list_abc = {ThreeBodyStorage::ABC};
            std::vector<ThreeBodyStorage::Permutation> perm_list_def = {ThreeBodyStorage::ABC};

            if (bra.p != bra.q and bra.q == bra.r)
              perm_list_abc = {ThreeBodyStorage::CBA};
            if (ket.p != ket.q and ket.q == ket.r)
              perm_list_def = {ThreeBodyStorage::CBA};
            if (bra.p == bra.q and bra.q == bra.r)
              perm_list_abc = {ThreeBodyStorage::ABC, ThreeBodyStorage::CBA, ThreeBodyStorage::ACB};
            if (ket.p == ket.q and ket.q == ket.r)
              perm_list_def = {ThreeBodyStorage::ABC, ThreeBodyStorage::CBA, ThreeBodyStorage::ACB};
            for (auto perm_abc : perm_list_abc)
            {
              size_t a, b, c, d, e, f;
              Rando.ThreeBody.Permute(perm_abc, bra.p, bra.q, bra.r, a, b, c);
              double ja = 0.5 * Rando.modelspace->GetOrbit(a).j2;
              double jb = 0.5 * Rando.modelspace->GetOrbit(b).j2;
              double jc = 0.5 * Rando.modelspace->GetOrbit(c).j2;
              int Jab_min = std::abs(ja - jb);
              int Jab_max = (ja + jb);
              for (int Jab = Jab_min; Jab <= Jab_max; Jab++)
              {
                if ((a == b) and Jab % 2 != 0)
                  continue;
                // double recouple_bra = Rando.ThreeBody.PermutationPhase(perm_abc)  * Rando.ThreeBody.RecouplingCoefficient( perm_abc, ja, jb, jc, Jab, bra.Jpq, twoJ);
                double recouple_bra = Rando.ThreeBody.PermutationPhase(perm_abc) * Rando.ThreeBody.RecouplingCoefficient(perm_abc, ja, jb, jc, bra.Jpq, Jab, twoj1);
                if (std::abs(recouple_bra) < 1e-9)
                  continue;

                for (auto perm_def : perm_list_def)
                {
                  Rando.ThreeBody.Permute(perm_def, ket.p, ket.q, ket.r, d, e, f);
                  double jd = 0.5 * Rando.modelspace->GetOrbit(d).j2;
                  double je = 0.5 * Rando.modelspace->GetOrbit(e).j2;
                  double jf = 0.5 * Rando.modelspace->GetOrbit(f).j2;
                  int Jde_min = std::abs(jd - je);
                  int Jde_max = (jd + je);
                  for (int Jde = Jde_min; Jde <= Jde_max; Jde++)
                  {
                    if ((d == e) and Jde % 2 != 0)
                      continue;
                    // double recouple_ket = Rando.ThreeBody.PermutationPhase(perm_def) * Rando.ThreeBody.RecouplingCoefficient( perm_def, jd, je, jf, Jde, ket.Jpq, twoJ);
                    double recouple_ket = Rando.ThreeBody.PermutationPhase(perm_def) * Rando.ThreeBody.RecouplingCoefficient(perm_def, jd, je, jf, ket.Jpq, Jde, twoj2);
                    double fakematel = randbuf[(Jab + Jde) % bufsize];
                    if (ibra == iket and Jde > Jab)
                      fakematel *= herm;
                    random_me += recouple_bra * recouple_ket * fakematel;
                  }

                } // for perm_def
              }
            } // for perm_abc
          }
          else // otherwise, we have the easy case, all the matrix elements are linearly independent and can be set independently
          {
            double random_me = distribution(generator);
          }

          Rando.ThreeBody.SetME_pn_ch(ch_bra, ch_ket, ibra, iket, random_me);
        }
      }
    }
  }

  std::cout << "In  " << __func__ << "  norm of 1b : " << Rando.OneBodyNorm();
  if (particle_rank > 1)
    std::cout << "   norm of 2b " << Rando.TwoBodyNorm();
  if (particle_rank > 2)
    std::cout << "   norm of 3b " << Rando.ThreeBodyNorm();
  std::cout << std::endl;

  return Rando;
}

// Generate a random dagger operator, i.e. one that changes particle number by +1
// for use in testing the dagger commutator routines.
// DaggerOperator UnitTest::RandomDaggerOp(ModelSpace& modelspace, index_t Q)
Operator UnitTest::RandomDaggerOp(ModelSpace &modelspace, index_t Q)
{
  //   Operator dag(modelspace,Q);
  Operator dag(modelspace);
  dag.SetNumberLegs(3);
  dag.SetQSpaceOrbit(Q);
  dag.SetNonHermitian();

  std::default_random_engine generator(random_seed);
  double mean = 0;
  double stddev = 1;
  std::normal_distribution<double> distribution(mean, stddev);

  // Only orbits with the same l,j,tz quantum numbers as the Q orbit
  // can be non-zero, so we fill a single column of the one-body matrix.
  Orbit &oQ = modelspace.GetOrbit(Q);
  for (auto i : dag.OneBodyChannels.at({oQ.l, oQ.j2, oQ.tz2}))
  {
    double random_me = distribution(generator);
    dag.OneBody(i, 0) = random_me;
    //     dag.OneBody(i,Q)= random_me;
  }

  size_t nch = modelspace.GetNumberTwoBodyChannels();
  for (size_t ch = 0; ch < nch; ch++)
  {
    TwoBodyChannel &tbc = modelspace.GetTwoBodyChannel(ch);
    size_t nkets = tbc.GetNumberKets();
    for (size_t ibra = 0; ibra < nkets; ibra++)
    {
      Ket &bra = tbc.GetKet(ibra);
      for (auto k : modelspace.all_orbits)
      {
        Orbit &ok = modelspace.GetOrbit(k);
        // check whether |kQ> lives in this channel
        if ((ok.tz2 + oQ.tz2 != 2 * tbc.Tz) or ((ok.l + oQ.l) % 2 != tbc.parity) or ((ok.j2 + oQ.j2) < 2 * tbc.J) or (std::abs(ok.j2 - oQ.j2) > 2 * tbc.J))
          continue;
        //         if ( not tbc.CheckChannel_ket( &ok, &oQ) ) continue;

        double random_me = distribution(generator);
        //         dag.TwoBody.SetTBME(ch,ch,bra.p,bra.q,k,Q, random_me);
        dag.ThreeLeg.SetME(ch, bra.p, bra.q, k, random_me);
      }
    }
    //     std::cout << "ch = " << ch << "  matrix looks like " << std::endl << dag.ThreeLeg.GetMatrix(ch) << std::endl << std::endl;
  }

  return dag;
}

void UnitTest::Test3BodyHermiticity(Operator &Y)
{

  std::cout << " Hermitian? " << Y.IsHermitian() << std::endl;

  int a = 0, b = 0, c = 10, d = 0, e = 10, f = 0;
  double yabcdef = Y.ThreeBody.GetME_pn(0, 0, 1, a, b, c, d, e, f);
  double ydefabc = Y.ThreeBody.GetME_pn(0, 0, 1, d, e, f, a, b, c);
  std::cout << " yabcdef = " << yabcdef << "   ydefabc = " << ydefabc << std::endl;
  //  int ma=1,mb=-1,mc=1,md=1,me=-1,mf=1;
  //  double yabcdef = GetMschemeMatrixElement_3b(Y, a,ma, b,mb, c,mc, d,md, e,me, f,mf)
}

// void UnitTest::Test3BodyAntisymmetry()
void UnitTest::Test3BodyAntisymmetry(Operator &Y)
{

  arma::arma_rng::set_seed(random_seed);
  //  Operator Y = RandomOp(*modelspace, 0, 0, 0, 3, 1);

  bool all_good = true;

  std::vector<std::array<int, 6>> indices;
  //  indices.push_back(  {4,4,0,4,4,0} );
  //  indices.push_back(  {0,4,4,0,2,2} );
  indices.push_back({3, 3, 2, 3, 3, 2});
  //  indices.push_back(  {2,0,2,0,2,4} );
  //  indices.push_back(  {3,3,3,1,3,3} );

  //  int a=0,b=0,c=3,d=2,e=2,f=5;
  //  int a=0,b=0,c=5,d=2,e=2,f=3;
  //  int a=4,b=4,c=0,d=4,e=4,f=0;
  //  int a=5,b=0,c=0,d=5,e=0,f=0;
  //  int a=5,b=4,c=0,d=5,e=4,f=0;

  for (auto ind : indices)
  {
    int a, b, c, d, e, f;
    a = ind[0];
    b = ind[1];
    c = ind[2];
    d = ind[3];
    e = ind[4];
    f = ind[5];

    std::cout << " abcdef " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
    Orbit &oa = Y.modelspace->GetOrbit(a);
    Orbit &ob = Y.modelspace->GetOrbit(b);
    Orbit &oc = Y.modelspace->GetOrbit(c);
    Orbit &od = Y.modelspace->GetOrbit(d);
    Orbit &oe = Y.modelspace->GetOrbit(e);
    Orbit &of = Y.modelspace->GetOrbit(f);

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

    for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
    {
      for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
      {
        for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
        {
          for (int md = -od.j2; md <= od.j2; md += 2)
          {
            for (int me = -oe.j2; me <= oe.j2; me += 2)
            {
              int mf = ma + mb + mc - md - me;
              if (std::abs(mf) > of.j2)
                continue;
              double yabcdef = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, c, mc, d, md, e, me, f, mf);
              if (std::abs(yabcdef) < 1e-8)
                continue;
              //  if ( not (  ma==-3 and mb==-1 and mc==1 and md==-3 and me==3 and mf==-3) ) continue;
              //  std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
              //  std::cout << " @mvals: " << ma << " " << mb << " " << mc << " " << md << " " << me << " " << mf << std::endl;
              //  std::cout << std::endl << "@abc" << std::endl;
              //  std::cout << "@bca" << std::endl;
              double ybcadef = GetMschemeMatrixElement_3b(Y, b, mb, c, mc, a, ma, d, md, e, me, f, mf);
              //  std::cout << "@cab" << std::endl;
              double ycabdef = GetMschemeMatrixElement_3b(Y, c, mc, a, ma, b, mb, d, md, e, me, f, mf);

              //  std::cout << "@acb" << std::endl;
              double yacbdef = GetMschemeMatrixElement_3b(Y, a, ma, c, mc, b, mb, d, md, e, me, f, mf);
              //  std::cout << "@bac" << std::endl;
              double ybacdef = GetMschemeMatrixElement_3b(Y, b, mb, a, ma, c, mc, d, md, e, me, f, mf);
              //  std::cout << std::endl << "@cba" << std::endl;
              double ycbadef = GetMschemeMatrixElement_3b(Y, c, mc, b, mb, a, ma, d, md, e, me, f, mf);

              std::cout << " @@mvals: " << ma << " " << mb << " " << mc << " " << md << " " << me << " " << mf << std::endl;
              std::cout << "ANTISYMMETRY:  abc: " << yabcdef << " bca: " << ybcadef << " cab: " << ycabdef
                        << "    acb: " << yacbdef << " bac: " << ybacdef << " cba: " << ycbadef << std::endl;

              all_good &= (std::abs(yabcdef - ybcadef) < 1e-7 and std::abs(yabcdef - ycabdef) < 1e-7 and std::abs(yabcdef + yacbdef) < 1e-7 and std::abs(yabcdef + ybacdef) < 1e-7 and std::abs(yabcdef + ycbadef) < 1e-7);
              //  std::cout << "ANTISYMMETRY:  abc: " << yabcdef  << " cba: " << ycbadef << std::endl;
              //  return;

            } // me
          } // md
        } // mc
      } // mb
    } // ma
  }

  if (all_good)
  {
    std::cout << "PASSED ANTISYMMETRY TEST " << std::endl;
  }
  else
  {
    std::cout << "!!!! FAILED ANTISYMMETRY TEST!!!!!!!!!!!!!!!!!" << std::endl;
  }
}




bool UnitTest::TestNormalOrdering(Operator& Op)
{
   bool passed = true;
   Operator OpNO = Op.DoNormalOrdering();


   // Check the zero body piece
   double opM_zero = Op.ZeroBody;
   // one-body -> zero body
   for (auto a : Op.modelspace->holes)
   {
     Orbit& oa = Op.modelspace->GetOrbit(a);
     for (int m2a= -oa.j2; m2a<= oa.j2; m2a+=2)
     {
       opM_zero += oa.occ * GetMschemeMatrixElement_1b(  Op,  a,  m2a, a, m2a );
     }
   }

   // two-body -> zero body
   for (auto a : Op.modelspace->holes)
   {
     Orbit& oa = Op.modelspace->GetOrbit(a);
     for (auto b : Op.modelspace->holes)
     {
        if (b>a) continue;
        Orbit& ob = Op.modelspace->GetOrbit(b);
        for (int m2a= -oa.j2; m2a<= oa.j2; m2a+=2)
        {
           for (int m2b= -ob.j2; m2b<= ob.j2; m2b+=2)
           {
              if ( a==b and m2b>=m2a) continue; 
              opM_zero += oa.occ * ob.occ * GetMschemeMatrixElement_2b(  Op,  a,m2a, b,m2b, a,m2a, b,m2b );
           }
        }
        
     }
   }

   // three-body -> zero body
   if ( Op.GetParticleRank()>=3 )
   {
     for (auto a : Op.modelspace->holes)
     {
       Orbit& oa = Op.modelspace->GetOrbit(a);
       for (auto b : Op.modelspace->holes)
       {
          if (b>a) continue;
          Orbit& ob = Op.modelspace->GetOrbit(b);
          for (auto c : Op.modelspace->holes)
          {
             if (c>b) continue;
             Orbit& oc = Op.modelspace->GetOrbit(c);
             for (int m2a= -oa.j2; m2a<= oa.j2; m2a+=2)
             {
                for (int m2b= -ob.j2; m2b<= ob.j2; m2b+=2)
                {
                   if ( a==b and m2b>=m2a) continue; 
                   for (int m2c= -oc.j2; m2c<= oc.j2; m2c+=2)
                   {
                      if ( b==c and m2c>=m2b) continue; 
                      opM_zero += oa.occ * ob.occ * oc.occ * GetMschemeMatrixElement_3b(  Op,  a,m2a, b,m2b, c,m2c, a,m2a, b,m2b, c,m2c );
                   }
                }
             }
          }
       }
     }
   }



   double diff0b = OpNO.ZeroBody - opM_zero;
   std::cout << "0b:   m-scheme = " << opM_zero << "   J-scheme = " << OpNO.ZeroBody << "    diff = " << diff0b << std::endl;
   passed &= ( std::abs(diff0b) < 1e-6);


   // next check the 1b piece
   arma::mat opM_one =  Op.OneBody;
      // two-body -> one-body
     for (auto i : Op.modelspace->all_orbits)
     {
       Orbit& oi = Op.modelspace->GetOrbit(i);
       for (auto j : Op.modelspace->all_orbits)
       {
          Orbit& oj = Op.modelspace->GetOrbit(j);
          double Oij = 0;
          int m2i = 1;
          int m2j = 1;
          for (auto a : Op.modelspace->holes)
          {
             Orbit& oa = Op.modelspace->GetOrbit(a);
             for (int m2a= -oa.j2; m2a<= oa.j2; m2a+=2)
             {
                 Oij += oa.occ * GetMschemeMatrixElement_2b(  Op,  i,m2i, a,m2a, j,m2j, a,m2a );
             }
          }
          if (not Op.IsReduced() )
          {
             opM_one(i,j) += Oij;
          }
          else
          {
             int lambda = Op.GetJRank();
             int mu = (m2i -m2j)/2;
             double cg = AngMom::CG(oj.j2*0.5, m2j*0.5,  lambda,mu,  oi.j2*0.5, m2i*0.5);

             if (std::abs(cg)>1e-7)
             {
               double opMij = opM_one(i,j);
               opM_one(i,j) += sqrt(oi.j2+1) / cg * Oij;
             }
             else if ( std::abs(Oij)>1e-6)
             {
                std::cout << "Clebsch was bad for ij = " << i << " " << j << std::endl;
             }
          }
       }
     }

      // three-body -> one-body
     for (auto i : Op.modelspace->all_orbits)
     {
       Orbit& oi = Op.modelspace->GetOrbit(i);
       for (auto j : Op.modelspace->all_orbits)
       {
          Orbit& oj = Op.modelspace->GetOrbit(j);
          double Oij = 0;
          int m2i = 1;
          int m2j = 1;
          for (auto a : Op.modelspace->holes)
          {
             Orbit& oa = Op.modelspace->GetOrbit(a);
             for (int m2a= -oa.j2; m2a<= oa.j2; m2a+=2)
             {
               for (auto b : Op.modelspace->holes)
               {
                 if ( b>a) continue;
                 Orbit& ob = Op.modelspace->GetOrbit(b);
                 for (int m2b= -ob.j2; m2b<= ob.j2; m2b+=2)
                 {
                    if ( a==b and m2b>=m2a) continue;
                    Oij += oa.occ * ob.occ * GetMschemeMatrixElement_3b(  Op,  i,m2i, a,m2a, b,m2b, j,m2j, a,m2a, b,m2b );
                 }
               }
            }
          }
          if (not Op.IsReduced() )
          {
             opM_one(i,j) += Oij;
          }
          else
          {
             int lambda = Op.GetJRank();
             int mu = (m2i -m2j)/2;
             double cg = AngMom::CG(oj.j2*0.5, m2j*0.5,  lambda,mu,  oi.j2*0.5, m2i*0.5);

             if (std::abs(cg)>1e-7)
             {
               double opMij = opM_one(i,j);
               opM_one(i,j) += sqrt(oi.j2+1) / cg * Oij;
             }
             else if ( std::abs(Oij)>1e-6)
             {
                std::cout << "Clebsch was bad for ij = " << i << " " << j << std::endl;
             }
          }
       }
     }

     arma::mat diff1b = OpNO.OneBody - opM_one;
     std::cout << "1b:   m-scheme = " << arma::norm( opM_one ) << "   J-scheme = " << arma::norm(OpNO.OneBody)  << "    diff = " << arma::norm(diff1b) << std::endl;
     passed &= ( arma::norm(diff1b) < 1e-6);


    if ( Op.GetParticleRank()>=3)
    {
      double summed_diff_2b = 0;
      double summed_Oijkl = 0;
      // finally, check the 2b piece
      for (auto i : Op.modelspace->all_orbits)
      {
        Orbit& oi = Op.modelspace->GetOrbit(i);
        for (auto j : Op.modelspace->all_orbits)
        {
           Orbit& oj = Op.modelspace->GetOrbit(j);
           for (auto k : Op.modelspace->all_orbits)
           {
             Orbit& ok = Op.modelspace->GetOrbit(k);
             for (auto l : Op.modelspace->all_orbits)
             {
                Orbit& ol = Op.modelspace->GetOrbit(l);
                double Oijkl = 0;
                int m2i = 1;
                int m2j = 1;
                int m2k = 1;
                int m2l = 1;
                for (auto a : Op.modelspace->holes)
                {
                   Orbit& oa = Op.modelspace->GetOrbit(a);
                   for (int m2a= -oa.j2; m2a<= oa.j2; m2a+=2)
                   {
                       Oijkl += oa.occ * GetMschemeMatrixElement_3b(  Op,  i,m2i, j,m2j, a,m2a, k,m2k, l,m2l, a,m2a );
                   }
                }
                double oijkl_before = GetMschemeMatrixElement_2b( Op, i,m2i,j,m2j,k,m2k,l,m2l);
                double oijkl_after  = GetMschemeMatrixElement_2b( OpNO, i,m2i,j,m2j,k,m2k,l,m2l);
                double diff = ( oijkl_after - oijkl_before - Oijkl);
                if (std::abs(diff)>1e-6)
                {
                  std::cout << "Trouble !  " << Oijkl << " != " << oijkl_after - oijkl_before << std::endl;
                }
                summed_diff_2b += std::abs(diff);
                summed_Oijkl += std::abs(Oijkl);
             }// for l
           }//for k
        }//for j
      }//for i
      std::cout << "2b:   summed_diff_2b = " << summed_diff_2b << "  sum |Oijkl| = " << summed_Oijkl << std::endl;
      passed &= std::abs(summed_diff_2b) < 1e-6;
     }

   return passed;
}




bool UnitTest::TestCommutators()
{
  double t_start = omp_get_wtime();
  arma::arma_rng::set_seed(random_seed);
  Operator X = RandomOp(*modelspace, 0, 0, 0, 3, -1);
  Operator Y = RandomOp(*modelspace, 0, 0, 0, 3, +1);
  modelspace->PreCalculateSixJ();

  bool all_good = true;

  if (Commutator::comm_term_on["comm110ss"])
    all_good &= Test_comm110ss(X, Y);
  if (Commutator::comm_term_on["comm220ss"])
    all_good &= Test_comm220ss(X, Y);

  if (Commutator::comm_term_on["comm111ss"])
    all_good &= Test_comm111ss(X, Y);
  if (Commutator::comm_term_on["comm121ss"])
    all_good &= Test_comm121ss(X, Y);
  if (Commutator::comm_term_on["comm221ss"])
    all_good &= Test_comm221ss(X, Y);

  if (Commutator::comm_term_on["comm122ss"])
    all_good &= Test_comm122ss(X, Y);
  if (Commutator::comm_term_on["comm222_pp_hhss"])
    all_good &= Test_comm222_pp_hhss(X, Y);
  if (Commutator::comm_term_on["comm222_phss"])
    all_good &= Test_comm222_phss(X, Y);
  //  if ( Commutator::comm_term_on["comm222_pp_hh_221ss"] )   all_good &= Test_comm222_pp_hh_221ss( X, Y );
  if (Commutator::comm_term_on["comm222_pp_hh"] and Commutator::comm_term_on["comm221ss"])
    all_good &= Test_comm222_pp_hh_221ss(X, Y);

  if (Commutator::comm_term_on["comm330ss"])
    all_good &= Test_comm330ss(X, Y);
  if (Commutator::comm_term_on["comm331ss"])
    all_good &= Test_comm331ss(X, Y);
  if (Commutator::comm_term_on["comm231ss"])
    all_good &= Test_comm231ss(X, Y);
  if (Commutator::comm_term_on["comm132ss"])
    all_good &= Test_comm132ss(X, Y);
  if (Commutator::comm_term_on["comm232ss"])
    all_good &= Test_comm232ss(X, Y);
  if (Commutator::comm_term_on["comm223ss"])
    all_good &= Test_comm223ss(X, Y);
  if (Commutator::comm_term_on["comm133ss"])
    all_good &= Test_comm133ss(X, Y);

  if (Commutator::comm_term_on["comm332_ppph_hhhpss"])
    all_good &= Test_comm332_ppph_hhhpss(X, Y);
  if (Commutator::comm_term_on["comm332_pphhss"])
    all_good &= Test_comm332_pphhss(X, Y);
  if (Commutator::comm_term_on["comm233_pp_hhss"])
    all_good &= Test_comm233_pp_hhss(X, Y);
  if (Commutator::comm_term_on["comm233_phss"])
    all_good &= Test_comm233_phss(X, Y);

  if (Commutator::comm_term_on["comm333_ppp_hhhss"])
    all_good &= Test_comm333_ppp_hhhss(X, Y);
  if (Commutator::comm_term_on["comm333_pph_hhpss"])
    all_good &= Test_comm333_pph_hhpss(X, Y);

  if (all_good)
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }

  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return all_good;
}


//bool UnitTest::TestCommutators_Tensor()
bool UnitTest::TestCommutators_Tensor(Operator& X, Operator& Y)
{
  double t_start = omp_get_wtime();
//  arma::arma_rng::set_seed(random_seed);
//  Operator X = RandomOp(*modelspace, 0, 0, 0, 3, -1); // generator-like. Jrank=0, even parity, Trank=0, antihermitian
////  Operator Y = RandomOp(*modelspace, 0, 0, 0, 3, +1); // Jrank = 1,  even parity, Trank =0, hermitian
//  Operator Y = RandomOp(*modelspace, 1, 0, 0, 3, +1); // Jrank = 1,  even parity, Trank =0, hermitian
  X.modelspace->PreCalculateSixJ();

  bool all_good = true;

//  Operator Ynred(Y);
//  Ynred.MakeNotReduced();
//  all_good &= Mscheme_Test_comm121ss(X,Ynred);
//  all_good &= Mscheme_Test_comm121st(X,Y);
//  all_good &= Mscheme_Test_comm221st(X,Y);
//  all_good &= Mscheme_Test_comm222_pp_hhst(X,Y);
//  all_good &= Mscheme_Test_comm222_phst(X,Y);
//  X.OneBody *= 0;

  if (Commutator::comm_term_on["comm111st"])
               all_good &= Test_comm111st(X, Y);
  if (Commutator::comm_term_on["comm121st"])
               all_good &= Test_comm121st(X, Y);
  if (Commutator::comm_term_on["comm122st"])
               all_good &= Test_comm122st(X, Y);
  if (Commutator::comm_term_on["comm221st"])
               all_good &= Test_comm221st(X, Y);
  if (Commutator::comm_term_on["comm222_pp_hhst"])
               all_good &= Test_comm222_pp_hhst(X, Y);
  if (Commutator::comm_term_on["comm222_phst"])
               all_good &= Test_comm222_phst(X, Y);

  // 3n7 tensor commutators
  if (Commutator::comm_term_on["comm331st"])
               all_good &= Test_comm331st(X, Y);
  if (Commutator::comm_term_on["comm231st"])
               all_good &= Test_comm231st(X, Y);
  if (Commutator::comm_term_on["comm232st"])
               all_good &= Test_comm232st(X, Y);
  if (Commutator::comm_term_on["comm132st"])
               all_good &= Test_comm132st(X, Y);
  if (Commutator::comm_term_on["comm223st"])
               all_good &= Test_comm223st(X, Y);
  if (Commutator::comm_term_on["comm133st"])
               all_good &= Test_comm133st(X, Y);

  // 3n8 and 3n9 commutators
  if (Commutator::comm_term_on["comm332_pphhst"])
               all_good &= Test_comm332_pphhst(X, Y);
  if (Commutator::comm_term_on["comm332_ppph_hhhpst"])
               all_good &= Test_comm332_ppph_hhhpst(X, Y);
  if (Commutator::comm_term_on["comm233_pp_hhst"])
               all_good &= Test_comm233_pp_hhst(X, Y);
  if (Commutator::comm_term_on["comm233_phst"])
               all_good &= Test_comm233_phst(X, Y);
  if (Commutator::comm_term_on["comm333_ppp_hhhst"])
               all_good &= Test_comm333_ppp_hhhst(X, Y);
  if (Commutator::comm_term_on["comm333_pph_hhpst"])
               all_good &= Test_comm333_pph_hhpst(X, Y);








  if (all_good)
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }

  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return all_good;
}




bool UnitTest::TestCommutators_IsospinChanging()
{
  double t_start = omp_get_wtime();
  arma::arma_rng::set_seed(random_seed);
  Operator X = RandomOp(*modelspace, 0, 0, 0, 3, -1);
  Operator Y = RandomOp(*modelspace, 0, 1, 0, 3, +1);
  modelspace->PreCalculateSixJ();

  bool all_good = true;

  if (Commutator::comm_term_on["comm110ss"])
    all_good &= Test_comm110ss(X, Y);
  if (Commutator::comm_term_on["comm220ss"])
    all_good &= Test_comm220ss(X, Y);

  if (Commutator::comm_term_on["comm111ss"])
    all_good &= Test_comm111ss(X, Y);
  if (Commutator::comm_term_on["comm121ss"])
    all_good &= Test_comm121ss(X, Y);
  if (Commutator::comm_term_on["comm221ss"])
    all_good &= Test_comm221ss(X, Y);

  if (Commutator::comm_term_on["comm122ss"])
    all_good &= Test_comm122ss(X, Y);
  if (Commutator::comm_term_on["comm222_pp_hhss"])
    all_good &= Test_comm222_pp_hhss(X, Y);
  if (Commutator::comm_term_on["comm222_phss"])
    all_good &= Test_comm222_phss(X, Y);
  //  if ( Commutator::comm_term_on["comm222_pp_hh_221ss"] )   all_good &= Test_comm222_pp_hh_221ss( X, Y );
  if (Commutator::comm_term_on["comm222_pp_hh"] and Commutator::comm_term_on["comm221ss"])
    all_good &= Test_comm222_pp_hh_221ss(X, Y);

  if (Commutator::comm_term_on["comm330ss"])
    all_good &= Test_comm330ss(X, Y);
  if (Commutator::comm_term_on["comm331ss"])
    all_good &= Test_comm331ss(X, Y);
  if (Commutator::comm_term_on["comm231ss"])
    all_good &= Test_comm231ss(X, Y);
  if (Commutator::comm_term_on["comm132ss"])
    all_good &= Test_comm132ss(X, Y);
  if (Commutator::comm_term_on["comm232ss"])
    all_good &= Test_comm232ss(X, Y);
  if (Commutator::comm_term_on["comm223ss"])
    all_good &= Test_comm223ss(X, Y);
  if (Commutator::comm_term_on["comm133ss"])
    all_good &= Test_comm133ss(X, Y);

  if (Commutator::comm_term_on["comm332_ppph_hhhpss"])
    all_good &= Test_comm332_ppph_hhhpss(X, Y);
  if (Commutator::comm_term_on["comm332_pphhss"])
    all_good &= Test_comm332_pphhss(X, Y);
  if (Commutator::comm_term_on["comm233_pp_hhss"])
    all_good &= Test_comm233_pp_hhss(X, Y);
  if (Commutator::comm_term_on["comm233_phss"])
    all_good &= Test_comm233_phss(X, Y);

  if (Commutator::comm_term_on["comm333_ppp_hhhss"])
    all_good &= Test_comm333_ppp_hhhss(X, Y);
  if (Commutator::comm_term_on["comm333_pph_hhpss"])
    all_good &= Test_comm333_pph_hhpss(X, Y);

  if (all_good)
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }

  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return all_good;
}

bool UnitTest::TestCommutators_ParityChanging()
{
  double t_start = omp_get_wtime();
  arma::arma_rng::set_seed(random_seed);
  Operator X = RandomOp(*modelspace, 0, 0, 0, 2, -1);
  Operator Y = RandomOp(*modelspace, 0, 0, 1, 2, +1);
  modelspace->PreCalculateSixJ();
  Y.MakeNotReduced();

  bool all_good = true;

  if (Commutator::comm_term_on["comm110ss"])
    all_good &= Test_comm110ss(X, Y);
  if (Commutator::comm_term_on["comm220ss"])
    all_good &= Test_comm220ss(X, Y);

  if (Commutator::comm_term_on["comm111ss"])
    all_good &= Test_comm111ss(X, Y);
  if (Commutator::comm_term_on["comm121ss"])
    all_good &= Test_comm121ss(X, Y);
  if (Commutator::comm_term_on["comm221ss"])
    all_good &= Test_comm221ss(X, Y);

  if (Commutator::comm_term_on["comm122ss"])
    all_good &= Test_comm122ss(X, Y);
  if (Commutator::comm_term_on["comm222_pp_hhss"])
    all_good &= Test_comm222_pp_hhss(X, Y);
  if (Commutator::comm_term_on["comm222_phss"])
    all_good &= Test_comm222_phss(X, Y);
  //  if ( Commutator::comm_term_on["comm222_pp_hh_221ss"] )   all_good &= Test_comm222_pp_hh_221ss( X, Y );
  if (Commutator::comm_term_on["comm222_pp_hh"] and Commutator::comm_term_on["comm221ss"])
    all_good &= Test_comm222_pp_hh_221ss(X, Y);

  if (Commutator::comm_term_on["comm330ss"])
    all_good &= Test_comm330ss(X, Y);
  if (Commutator::comm_term_on["comm331ss"])
    all_good &= Test_comm331ss(X, Y);
  if (Commutator::comm_term_on["comm231ss"])
    all_good &= Test_comm231ss(X, Y);
  if (Commutator::comm_term_on["comm132ss"])
    all_good &= Test_comm132ss(X, Y);
  if (Commutator::comm_term_on["comm232ss"])
    all_good &= Test_comm232ss(X, Y);
  if (Commutator::comm_term_on["comm223ss"])
    all_good &= Test_comm223ss(X, Y);
  if (Commutator::comm_term_on["comm133ss"])
    all_good &= Test_comm133ss(X, Y);

  if (Commutator::comm_term_on["comm332_ppph_hhhpss"])
    all_good &= Test_comm332_ppph_hhhpss(X, Y);
  if (Commutator::comm_term_on["comm332_pphhss"])
    all_good &= Test_comm332_pphhss(X, Y);
  if (Commutator::comm_term_on["comm233_pp_hhss"])
    all_good &= Test_comm233_pp_hhss(X, Y);
  if (Commutator::comm_term_on["comm233_phss"])
    all_good &= Test_comm233_phss(X, Y);

  if (Commutator::comm_term_on["comm333_ppp_hhhss"])
    all_good &= Test_comm333_ppp_hhhss(X, Y);
  if (Commutator::comm_term_on["comm333_pph_hhpss"])
    all_good &= Test_comm333_pph_hhpss(X, Y);

  if (all_good)
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }

  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  return all_good;
}

// void UnitTest::TestCommutators3()
//  TODO: THIS IS DEPRECATED. I SHOULD PROBABLY JUST REMOVE IT... -SRS
void UnitTest::TestCommutators3(Operator &X, Operator &Y, std::vector<std::string> &skiplist)
{
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  std::cout << __func__ << "  is deprecated. Use UnitTest::TestCommutators() instead." << std::endl;
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  //  double t_start = omp_get_wtime();
  //  std::cout << " random_seed = " << random_seed << std::endl;
  //  arma::arma_rng::set_seed( random_seed );
  //  modelspace->PreCalculateSixJ();
  //  Operator X = RandomOp(*modelspace, 0, 0, 0, 3, -1);
  //  random_seed ++;
  //  Operator Y = RandomOp(*modelspace, 0, 0, 0, 3, +1);
  //  random_seed --;
  //  std::cout << "BEGIN " << __func__ << std::endl;

  //  Operator Xherm = Y;
  //  Xherm.TwoBody.Erase();
  //  Xherm.ThreeBody.Erase();
  //  X += RandomOp(*modelspace, 0, 0, 0, 2, -1);
  //  Y += RandomOp(*modelspace, 0, 0, 0, 2, +1);
  //  Xherm += RandomOp(*modelspace, 0, 0, 0, 2, +1);
  //  X.ThreeBody.Erase();
  //  Commutator::SetSingleThread(false);
  //  Commutator::comm223ss( Xherm, Y, X); // Make the 3-body part of X equal to the commutator of 2 hermitian 2b operators
  //  Commutator::comm133ss( Xherm, X, Y);
  //
  ////  X.modelspace->scalar3b_transform_first_pass = false;
  ////  Commutator::comm223ss( Xherm, Y, X);
  //  bool all_good = true;
  //
  //  std::cout << "3body norms of X and Y: " << X.ThreeBodyNorm() << " " << Y.ThreeBodyNorm() << std::endl;
  //
  //  all_good &= Test_comm330ss(X,Y);
  //
  //  if ( std::find(skiplist.begin(),skiplist.end(), "allN6") == skiplist.end() )
  //  {
  //     if ( std::find(skiplist.begin(),skiplist.end(), "comm330ss") == skiplist.end())   all_good &= Test_comm330ss( X, Y );
  //     if ( std::find(skiplist.begin(),skiplist.end(), "comm231ss") == skiplist.end())   all_good &= Test_comm231ss( X, Y );
  //     if ( std::find(skiplist.begin(),skiplist.end(), "comm132ss") == skiplist.end())   all_good &= Test_comm132ss( X, Y );
  //  }
  //
  //  if ( std::find(skiplist.begin(),skiplist.end(), "allN7") == skiplist.end() )
  //  {
  //     if ( std::find(skiplist.begin(),skiplist.end(), "comm331ss") == skiplist.end())   all_good &= Test_comm331ss( X, Y );
  //     if ( std::find(skiplist.begin(),skiplist.end(), "comm232ss") == skiplist.end())   all_good &= Test_comm232ss( X, Y );
  //     if ( std::find(skiplist.begin(),skiplist.end(), "comm223ss") == skiplist.end())   all_good &= Test_comm223ss( X, Y );
  //     if ( std::find(skiplist.begin(),skiplist.end(), "comm133ss") == skiplist.end())   all_good &= Test_comm133ss( X, Y );
  //  }
  //
  //  if ( std::find(skiplist.begin(),skiplist.end(), "allN8") == skiplist.end() )
  //  {
  //    if ( std::find(skiplist.begin(),skiplist.end(), "comm332_ppph_hhhpss") == skiplist.end())   all_good &= Test_comm332_ppph_hhhpss( X, Y );
  //    if ( std::find(skiplist.begin(),skiplist.end(), "comm332_pphhss") == skiplist.end())   all_good &= Test_comm332_pphhss( X, Y );
  //    if ( std::find(skiplist.begin(),skiplist.end(), "comm233_pp_hhss") == skiplist.end())   all_good &= Test_comm233_pp_hhss( X, Y );
  //    if ( std::find(skiplist.begin(),skiplist.end(), "comm233_phss") == skiplist.end())   all_good &= Test_comm233_phss( X, Y );
  //  }
  //  if ( std::find(skiplist.begin(),skiplist.end(), "allN9") == skiplist.end() )
  //  {
  //    if ( std::find(skiplist.begin(),skiplist.end(), "comm333_ppp_hhhss") == skiplist.end())   all_good &= Test_comm333_ppp_hhhss( X, Y );
  //    if ( std::find(skiplist.begin(),skiplist.end(), "comm333_pph_hhpss") == skiplist.end())   all_good &= Test_comm333_pph_hhpss( X, Y );
  //  }
  //
  //
  //
  //  if ( all_good )
  //  {
  //    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  //  }
  //  else
  //  {
  //    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  //  }
  //  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

void UnitTest::TestDaggerCommutators(index_t Q)
{
  double t_start = omp_get_wtime();
  arma::arma_rng::set_seed(random_seed);
  Operator X = RandomOp(*modelspace, 0, 0, 0, 2, -1);
  Operator Y = RandomDaggerOp(*modelspace, Q);

  Orbit &oQ = Y.modelspace->GetOrbit(Q);

  Y *= 0;
  std::vector<Operator> Yn;
  std::vector<Operator> Zn;
  for (size_t n : Y.OneBodyChannels.at({oQ.l, oQ.j2, oQ.tz2}))
  {
    Y.OneBody(n, 0) = 1.0;
    Yn.push_back(0 * Y);
    Yn.back().OneBody(n, 0) = 1.0;
    Zn.push_back(Commutator::Commutator(X, Yn.back()));
    Zn.back() = Commutator::Commutator(X, Zn.back());
  }

  Operator Z = Commutator::Commutator(X, Y);
  Z = Commutator::Commutator(X, Z);

  std::cout << "Y 1b: " << std::endl
            << Y.OneBody << std::endl;
  std::cout << "Yn 1b: " << std::endl;
  for (auto y : Yn)
    std::cout << y.OneBody << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "Z 1b: " << std::endl
            << Z.OneBody << std::endl;
  std::cout << "Zn 1b: " << std::endl;
  for (auto z : Zn)
    std::cout << z.OneBody << std::endl;
  std::cout << std::endl;

  Operator Zdiff = Z;
  std::cout << "legs for Z and Zdiff : " << Z.legs << " " << Zdiff.legs << std::endl;
  for (auto z : Zn)
    Zdiff -= z;
  std::cout << "norm of Zdiff = " << Zdiff.Norm() << "  norm Z = " << Z.Norm()
            << "  norm of Zn = ";
  for (auto z : Zn)
    std::cout << z.Norm() << " ";
  std::cout << std::endl;

  bool all_good = true;

  all_good &= Test_comm211sd(X, Y);
  all_good &= Test_comm231sd(X, Y);
  all_good &= Test_comm431sd(X, Y);
  all_good &= Test_comm413sd(X, Y);
  all_good &= Test_comm233sd(X, Y);
  all_good &= Test_comm433_pp_hh_sd(X, Y);
  all_good &= Test_comm433sd_ph(X, Y);

  if (all_good)
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
  arma::arma_rng::set_seed(random_seed);
  Operator X = RandomOp(*modelspace, 0, 0, 0, 2, -1);
  //  Operator Y = RandomDaggerOp(*modelspace, Q);

  Operator Y = imsrg_util::DaggerAlln_Op(*modelspace, Q);
  Operator Z = Commutator::Commutator(X, Commutator::Commutator(X, Y));
  std::vector<Operator> Yn;
  std::vector<Operator> Zn;
  Orbit &oQ = X.modelspace->GetOrbit(Q);
  for (auto nQ : X.modelspace->OneBodyChannels.at({oQ.l, oQ.j2, oQ.tz2}))
  {
    Yn.push_back(imsrg_util::Dagger_Op(*modelspace, nQ));
    //    Yn.back().OneBody(nQ,nQ) = 0;
    //    Yn.back().OneBody(nQ,Q) = 1.0;
    //    Yn.back().OneBody(nQ,0) = 0;
    Yn.back().OneBody(nQ, 0) = 1.0;
    Zn.push_back(Commutator::Commutator(X, Commutator::Commutator(X, Yn.back())));
  }
  std::cout << "Y one body: " << std::endl
            << Y.OneBody << std::endl;
  std::cout << "Z one body: " << std::endl
            << Z.OneBody << std::endl;
  for (size_t n = 0; n < Yn.size(); n++)
  {
    std::cout << " Y " << n << std::endl
              << Yn[n].OneBody << std::endl;
    std::cout << " Z " << n << std::endl
              << Zn[n].OneBody << std::endl;
  }
  Operator Zsum = Zn.front();
  for (size_t n = 1; n < Zn.size(); n++)
    Zsum += Zn[n];
  Operator Zdiff = Z - Zsum;
  double Znorm = Zdiff.Norm();
  std::cout << "Norms: " << Z.Norm() << "   vs  ";
  for (auto z : Zn)
    std::cout << z.Norm() << " ";
  std::cout << std::endl
            << "Norm of Zdiff = " << Znorm << std::endl;

  bool all_good = true;

  all_good &= Test_comm211sd(X, Y);
  all_good &= Test_comm231sd(X, Y);
  all_good &= Test_comm431sd(X, Y);
  all_good &= Test_comm413sd(X, Y);
  all_good &= Test_comm233sd(X, Y);
  all_good &= Test_comm433_pp_hh_sd(X, Y);
  all_good &= Test_comm433sd_ph(X, Y);

  if (all_good)
  {
    std::cout << " Done with " << __func__ << " and all is well" << std::endl;
  }
  else
  {
    std::cout << " Done with " << __func__ << " and at least one test failed" << std::endl;
  }
  X.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

double UnitTest::GetMschemeMatrixElement_1b(const Operator &Op, int a, int ma, int b, int mb)
{
  double matel = 0;
  int Jop = Op.GetJRank();
  int Top = Op.GetTRank();
//  if (Jop == 0 and Top == 0) // scalar operator
  if (not Op.IsReduced() ) // scalar operator
  {
    if (ma == mb)
    {
      matel = Op.OneBody(a, b);
    }
  }
  else // reduced matrix element using Wigner-Eckart thm:   <a ma | Op JM | b mb > = (-1)^2J * <J M b mb | a ma > / sqrt(2ja+1) < a || Op J || b >
  {
    Orbit &oa = Op.modelspace->GetOrbit(a);
    Orbit &ob = Op.modelspace->GetOrbit(b);
    if (oa.j2 + ob.j2 >= 2 * Jop and std::abs(oa.j2 - ob.j2) <= 2 * Jop and  2 * Jop >= std::abs(ma - mb))
    {
      matel = AngMom::CG(0.5 * ob.j2, 0.5 * mb, Jop, 0.5 * (ma - mb), 0.5 * oa.j2, 0.5 * ma) / sqrt(oa.j2 + 1.0) * Op.OneBody(a, b);
//      std::cout << __func__ << "  " << AngMom::CG(0.5 * ob.j2, 0.5 * mb, Jop, 0.5 * (ma - mb), 0.5 * oa.j2, 0.5 * ma) << " / " << sqrt(oa.j2 + 1.0) << " * " << Op.OneBody(a, b) << "  = " << matel << std::endl;
    }
  }
  return matel;
}

double UnitTest::GetMschemeMatrixElement_2b(const Operator &Op, int a, int ma, int b, int mb, int c, int mc, int d, int md)
{
  double matel = 0;
  int Jop = Op.GetJRank();
  int Top = Op.GetTRank();
  Orbit &oa = Op.modelspace->GetOrbit(a);
  Orbit &ob = Op.modelspace->GetOrbit(b);
  Orbit &oc = Op.modelspace->GetOrbit(c);
  Orbit &od = Op.modelspace->GetOrbit(d);
  if (a == b and ma == mb)
    return 0;
  if (c == d and mc == md)
    return 0;

//  if (Jop == 0) // scalar operator
  if (not Op.IsReduced() ) // scalar operator

  {
    if ((ma + mb) == (mc + md))
    {
      int Jmin = std::max(std::abs(oa.j2 - ob.j2), std::abs(oc.j2 - od.j2)) / 2;
      int Jmax = std::min(oa.j2 + ob.j2, oc.j2 + od.j2) / 2;
      int M = (ma + mb) / 2;
      for (int J = Jmin; J <= Jmax; J++)
      {
        // We take the un-normalized TBME (the one with the tilde) so we don't
        // need to worry about normalization factors.
        double clebsch_ab = AngMom::CG(0.5 * oa.j2, 0.5 * ma, 0.5 * ob.j2, 0.5 * mb, J, M);
        double clebsch_cd = AngMom::CG(0.5 * oc.j2, 0.5 * mc, 0.5 * od.j2, 0.5 * md, J, M);
        matel += clebsch_ab * clebsch_cd * Op.TwoBody.GetTBME_J(J, a, b, c, d);
      }
    }
  }
  else
  {
    if (abs(ma + mb - mc - md) <= 2 * Jop)
    {
      int J1_min = std::abs(oa.j2 - ob.j2) / 2;
      int J1_max = (oa.j2 + ob.j2) / 2;

      int J2_min = std::abs(oc.j2 - od.j2) / 2;
      int J2_max = (oc.j2 + od.j2) / 2;

      int M1 = (ma + mb) / 2;
      int M2 = (mc + md) / 2;
      for (int J1 = J1_min; J1 <= J1_max; J1++)
      {
        if (std::abs(M1) > J1)
          continue;
        if ( a == b and J1 % 2 > 0 )
          continue;
        
        for (int J2 = J2_min; J2 <= J2_max; J2++)
        {
          if (std::abs(M2) > J2)
            continue;
          if ( c == d and J2 % 2 > 0 )
            continue;

          if (J1 + J2 < Jop or std::abs(J1 - J2) > Jop)
            continue;

          // We take the un-normalized TBME (the one with the tilde) so we don't
          // need to worry about normalization factors.
          double clebsch_ab = AngMom::CG(0.5 * oa.j2, 0.5 * ma, 0.5 * ob.j2, 0.5 * mb, J1, M1);
          double clebsch_cd = AngMom::CG(0.5 * oc.j2, 0.5 * mc, 0.5 * od.j2, 0.5 * md, J2, M2);
          double clebsch_J1J2J = AngMom::CG(J2, M2, Jop, (M1 - M2), J1, M1);
          matel += clebsch_J1J2J * clebsch_ab * clebsch_cd / sqrt(2 * J1 + 1.0) * Op.TwoBody.GetTBME_J(J1, J2, a, b, c, d);
//          std::cout << __func__ << "  J1 J2 " << J1 << " " << J2 << " " << clebsch_J1J2J << " * "<< clebsch_ab << " * " << clebsch_cd << " / " << sqrt(2 * J1 + 1.0) << " * " << Op.TwoBody.GetTBME_J(J1, J2, a, b, c, d)
//                    << "  =  " << clebsch_J1J2J * clebsch_ab * clebsch_cd / sqrt(2 * J1 + 1.0) * Op.TwoBody.GetTBME_J(J1, J2, a, b, c, d) << "   =>  " << matel << std::endl;
        }
      }
    }
  }

  return matel;
}

double UnitTest::GetMschemeMatrixElement_3b(const Operator &Op, int a, int ma, int b, int mb, int c, int mc, int d, int md, int e, int me, int f, int mf)
{
  double matel = 0;
  int Jop = Op.GetJRank();
  int Tzop = Op.GetTRank();
  int Pop = Op.GetParity();
  Orbit &oa = Op.modelspace->GetOrbit(a);
  Orbit &ob = Op.modelspace->GetOrbit(b);
  Orbit &oc = Op.modelspace->GetOrbit(c);
  Orbit &od = Op.modelspace->GetOrbit(d);
  Orbit &oe = Op.modelspace->GetOrbit(e);
  Orbit &of = Op.modelspace->GetOrbit(f);
  if (a == b and ma == mb)
    return 0;
  if (a == c and ma == mc)
    return 0;
  if (b == c and mb == mc)
    return 0;
  if (d == e and md == me)
    return 0;
  if (d == f and md == mf)
    return 0;
  if (e == f and me == mf)
    return 0;

  if ((oa.l + ob.l + oc.l + od.l + oe.l + of.l + Pop) % 2 != 0)
    return 0;
  if (std::abs(oa.tz2 + ob.tz2 + oc.tz2 - od.tz2 - oe.tz2 - of.tz2) != 2 * Tzop)
    return 0;

  if (Jop == 0 and Tzop == 0 and Pop == 0)
  {
    if ((ma + mb + mc) != (md + me + mf))
      return 0;

    int Mab = (ma + mb) / 2;
    int Mde = (md + me) / 2;
    int twoM = ma + mb + mc;
    int Jab_min = std::max(std::abs(Mab), std::abs(oa.j2 - ob.j2) / 2);
    int Jab_max = (oa.j2 + ob.j2) / 2;
    int Jde_min = std::max(std::abs(Mde), std::abs(od.j2 - oe.j2) / 2);
    int Jde_max = (od.j2 + oe.j2) / 2;

    if (a == b)
      Jab_max = std::min(Jab_max, oa.j2 - 1);
    if (d == e)
      Jde_max = std::min(Jde_max, od.j2 - 1);

    for (int Jab = Jab_min; Jab <= Jab_max; Jab++)
    {
      if (a == b and (Jab % 2) > 0)
        continue;

      double clebsch_ab = AngMom::CG(0.5 * oa.j2, 0.5 * ma, 0.5 * ob.j2, 0.5 * mb, Jab, Mab);
      if (std::abs(clebsch_ab) < 1e-8)
        continue;

      for (int Jde = Jde_min; Jde <= Jde_max; Jde++)
      {
        if (d == e and (Jde % 2) > 0)
          continue;

        double clebsch_de = AngMom::CG(0.5 * od.j2, 0.5 * md, 0.5 * oe.j2, 0.5 * me, Jde, Mde);
        if (std::abs(clebsch_de) < 1e-8)
          continue;

        int twoJ_min = std::max({std::abs(twoM), std::abs(2 * Jab - oc.j2), std::abs(2 * Jde - of.j2)});
        int twoJ_max = std::min(2 * Jab + oc.j2, 2 * Jde + of.j2);
        for (int twoJ = twoJ_min; twoJ <= twoJ_max; twoJ += 2)
        {
          double clebsch_abc = AngMom::CG(Jab, Mab, 0.5 * oc.j2, 0.5 * mc, 0.5 * twoJ, 0.5 * twoM);
          double clebsch_def = AngMom::CG(Jde, Mde, 0.5 * of.j2, 0.5 * mf, 0.5 * twoJ, 0.5 * twoM);
          if (std::abs(clebsch_abc) < 1e-8 or std::abs(clebsch_def) < 1e-8)
            continue; // avoid the look up, with possible 6js...
          double meJ = Op.ThreeBody.GetME_pn(Jab, Jde, twoJ, a, b, c, d, e, f);
          matel += clebsch_ab * clebsch_abc * clebsch_de * clebsch_def * meJ;
          //        std::cout << __func__ << "  Jab, Jde, twoJ  " << Jab << " " << Jde << " " << twoJ
          //                  << "  cab cabc cde cdef " << clebsch_ab << " " << clebsch_abc << " " << clebsch_de << " " << clebsch_def
          //                  << "   meJ = " << meJ << "   ->  matel = " << matel << std::endl;
        }
      }
    }
  }
  else
  {
    if (std::abs(ma + mb + mc - md - me - mf) > 2 * Jop)
      return 0;

    int Mab = (ma + mb) / 2;
    int Mde = (md + me) / 2;
    int twoMabc = ma + mb + mc;
    int twoMdef = md + me + mf;
    int Jab_min = std::max(std::abs(Mab), std::abs(oa.j2 - ob.j2) / 2);
    int Jab_max = (oa.j2 + ob.j2) / 2;
    int Jde_min = std::max(std::abs(Mde), std::abs(od.j2 - oe.j2) / 2);
    int Jde_max = (od.j2 + oe.j2) / 2;
    int twoTm = ma + mb + mc - md - me - mf;

    if (a == b)
      Jab_max = std::min(Jab_max, oa.j2 - 1);
    if (d == e)
      Jde_max = std::min(Jde_max, od.j2 - 1);

    for (int Jab = Jab_min; Jab <= Jab_max; Jab++)
    {
      if (a == b and (Jab % 2) > 0)
        continue;
      double clebsch_ab = AngMom::CG(0.5 * oa.j2, 0.5 * ma, 0.5 * ob.j2, 0.5 * mb, Jab, Mab);
      if (std::abs(clebsch_ab) < 1e-8)
        continue;
      for (int Jde = Jde_min; Jde <= Jde_max; Jde++)
      {
        if (d == e and (Jde % 2) > 0)
          continue;
        double clebsch_de = AngMom::CG(0.5 * od.j2, 0.5 * md, 0.5 * oe.j2, 0.5 * me, Jde, Mde);
        if (std::abs(clebsch_de) < 1e-8)
          continue;

        int j1_min = std::max(std::abs(twoMabc), std::abs(2 * Jab - oc.j2));
        int j1_max = 2 * Jab + oc.j2;
        if (a == b and b == c)
          j1_max = std::min(j1_max, 3 * oa.j2 - 4);
        else if ((c == a or c == b) and a != b)
        {
          j1_max = std::min(j1_max, oa.j2 + ob.j2 + oc.j2 - 1);
        }

        int j2_min = std::max(std::abs(twoMdef), std::abs(2 * Jde - of.j2));
        int j2_max = 2 * Jde + of.j2;
        if (d == e and e == f)
          j2_max = std::min(j2_max, 3 * od.j2 - 4);
        else if ((f == d or f == e) and d != e)
        {
          j2_max = std::min(j2_max, od.j2 + oe.j2 + of.j2 - 1);
        }

        for (int twoj1 = j1_min; twoj1 <= j1_max; twoj1 += 2)
        {

          for (int twoj2 = j2_min; twoj2 <= j2_max; twoj2 += 2)
          {
            if (twoj1 + twoj2 < 2 * Jop or std::abs(twoj1 - twoj2) > 2 * Jop)
            {
              continue;
            }
            double clebsch_abc = AngMom::CG(Jab, Mab, 0.5 * oc.j2, 0.5 * mc, 0.5 * twoj1, 0.5 * twoMabc);
            double clebsch_def = AngMom::CG(Jde, Mde, 0.5 * of.j2, 0.5 * mf, 0.5 * twoj2, 0.5 * twoMdef);
            double CG_WE = AngMom::CG(0.5 * twoj2, 0.5 * twoMdef, Jop, 0.5 * twoTm, 0.5 * twoj1, 0.5 * twoMabc);

            if (std::abs(clebsch_abc) < 1e-8 or std::abs(clebsch_def) < 1e-8 or std::abs(CG_WE) < 1e-8)
              continue; // avoid the look up, with possible 6js...

            double meJ = Op.ThreeBody.GetME_pn(Jab, twoj1, Jde, twoj2, a, b, c, d, e, f);
            matel += CG_WE / sqrt(twoj1 + 1.) * clebsch_ab * clebsch_abc * clebsch_de * clebsch_def * meJ;
            //        std::cout << __func__ << "  Jab, Jde, twoJ  " << Jab << " " << Jde << " " << twoJ
            //        << "  cab cabc cde cdef " << clebsch_ab << " " << clebsch_abc << " " << clebsch_de << " " << clebsch_def
            //        << "   meJ = " << meJ << "   ->  matel = " << matel << std::endl;
          }
        }
      }
    }
  }

  return matel;
}

// The following two routines implicitly assume that we've chosen the apropriate projection quantum numbers.
// This may or may not be the most straightforward way to do things.
double UnitTest::GetMschemeMatrixElement_1leg(const Operator &Op, int a, int ma)
{

  return Op.OneBody(a, 0);
  //  return GetMschemeMatrixElement_1b( Op, a, ma, Op.GetQSpaceOrbit(), ma) ;
}

double UnitTest::GetMschemeMatrixElement_3leg(const Operator &Op, int a, int ma, int b, int mb, int c, int mc)
{

  double matel = 0;
  int Jop = Op.GetJRank();
  int Q = Op.GetQSpaceOrbit();
  Orbit &oa = Op.modelspace->GetOrbit(a);
  Orbit &ob = Op.modelspace->GetOrbit(b);
  Orbit &oc = Op.modelspace->GetOrbit(c);
  Orbit &oQ = Op.modelspace->GetOrbit(Q);
  int mQ = ma + mb - mc;
  if (mQ > oQ.j2)
    return 0;
  if (a == b and ma == mb)
    return 0;
  //  if (c==d and mc==md) return 0;
  if (Jop == 0) // scalar operator
  {
    if ((ma + mb) == (mc + mQ))
    {
      int Jmin = std::max(std::abs(oa.j2 - ob.j2), std::abs(oc.j2 - oQ.j2)) / 2;
      int Jmax = std::min(oa.j2 + ob.j2, oc.j2 + oQ.j2) / 2;
      int M = (ma + mb) / 2;
      for (int J = Jmin; J <= Jmax; J++)
      {
        // We take the un-normalized TBME (the one with the tilde) so we don't
        // need to worry about normalization factors.
        double clebsch_ab = AngMom::CG(0.5 * oa.j2, 0.5 * ma, 0.5 * ob.j2, 0.5 * mb, J, M);
        double clebsch_cd = AngMom::CG(0.5 * oc.j2, 0.5 * mc, 0.5 * oQ.j2, 0.5 * mQ, J, M);

        matel += clebsch_ab * clebsch_cd * Op.ThreeLeg.GetME_J(J, a, b, c);

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

bool UnitTest::Test_against_ref_impl(const Operator &X, const Operator &Y, commutator_func ComOpt, commutator_func ComRef, std::string output_tag)
{
  int z_Jrank = X.GetJRank() + Y.GetJRank(); // I sure hope this is zero.
  int z_Trank = X.GetTRank() + Y.GetTRank();
  int z_parity = (X.GetParity() + Y.GetParity()) % 2;
  int z_particlerank = Commutator::use_imsrg3 ? 3: 2;
  int hx = X.IsHermitian() ? +1 : -1;
  int hy = Y.IsHermitian() ? +1 : -1;

  Operator Xtmp, Ytmp;        // Declared, but not yet allocated, for reasons of scope.
  const Operator *Xnred = &X; // Pointer to the non-reduced version of the operator
  const Operator *Ynred = &Y;
  if (X.IsReduced() and z_Jrank ==0) // CommutatorScalarScalar doesn't expect reduced operators. Need to make it not reduced.
  {
    Xtmp = X;
    Xtmp.MakeNotReduced();
    Xnred = &Xtmp; // Now Xnred points to the not-reduced copy Xtmp
  }
  if (Y.IsReduced() and z_Jrank ==0)
  {
    Ytmp = Y;
    Ytmp.MakeNotReduced();
    Ynred = &Ytmp; // Now Ynred points to the not-reduced copy Ytmp
  }

  
  ModelSpace &ms = *(Y.GetModelSpace());
  Operator Z(ms, z_Jrank, z_Trank, z_parity, z_particlerank);
  if (hx*hy > 0)
  {
    Z.SetAntiHermitian();
  }
  if ( z_particlerank > 2 )
  {
     Z.ThreeBody.SetMode("pn");
  } 
//  Operator Z(Y);
//  Z.Erase();
  Operator Zref(Z);


  if (Z.IsReduced() and z_Jrank==0)
    Z.MakeNotReduced();

  if ((X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()))
  {
    Z.SetAntiHermitian();
    Zref.SetAntiHermitian();
  }
    
  else if ((X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()))
  {
    Z.SetHermitian();
    Zref.SetHermitian();
  }
  else
  {
    Z.SetNonHermitian();
    Zref.SetNonHermitian();
  }

  if (Zref.IsReduced() and z_Jrank==0)
    Zref.MakeNotReduced();


  ComOpt(*Xnred, *Ynred, Z);
  if ( not Z.IsReduced() and ((Z.GetParity() != 0) or (Z.GetTRank() != 0) and z_Jrank==0) )
  {  
    Z.MakeReduced(); // If Z changes parity or Tz, we by default store it as reduced. So make it as expected. Is that a good idea? Not sure....
  }

//  double tstart = omp_get_wtime();
  ComRef(*Xnred, *Ynred, Zref);
  if ( not Zref.IsReduced() and ((Zref.GetParity() != 0) or (Zref.GetTRank() != 0) and z_Jrank==0) )
  {
    Zref.MakeReduced(); // If Z changes parity or Tz, we by default store it as reduced. So make it as expected. Is that a good idea? Not sure....
  }
//  Z.profiler.timer["_ref_" + output_tag] += omp_get_wtime() - tstart;
  // std::cout<<Z.Norm()<<" "<<Z.ZeroBody<<std::endl;
  // std::cout << Zref.Norm() << " " << Zref.ZeroBody << std::endl;
  double normOpt = Z.Norm() + Z.ZeroBody;
  double normRef = Zref.Norm() + Zref.ZeroBody;
  Z -= Zref;
  double summed_error = Z.Norm() + Z.ZeroBody;

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << output_tag << "  sum_ref, sum_opt = " << normRef << " " << normOpt
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;

  return passed;
}

bool UnitTest::Test_comm110ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm110ss, ReferenceImplementations::comm110ss, "comm110ss");
}

bool UnitTest::Test_comm220ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm220ss, ReferenceImplementations::comm220ss, "comm220ss");
}

bool UnitTest::Test_comm111ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm111ss, ReferenceImplementations::comm111ss, "comm111ss");
}

bool UnitTest::Test_comm121ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm121ss, ReferenceImplementations::comm121ss, "comm121ss");
}

bool UnitTest::Test_comm221ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm221ss, ReferenceImplementations::comm221ss, "comm221ss");
}

bool UnitTest::Test_comm122ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm122ss, ReferenceImplementations::comm122ss, "comm122ss");
}

bool UnitTest::Test_comm222_pp_hhss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm222_pp_hhss, ReferenceImplementations::comm222_pp_hhss, "comm222_pp_hhss");
}

bool UnitTest::Test_comm222_phss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm222_phss, ReferenceImplementations::comm222_phss, "comm222_phss");
}
bool UnitTest::Test_comm222_pp_hh_221ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm222_pp_hh_221ss, ReferenceImplementations::comm222_pp_hh_221ss, "comm222_pp_hh_221ss");
}

/// scalar-tensor commutators
bool UnitTest::Test_comm111st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm111st, ReferenceImplementations::comm111st, "comm111st");
}
bool UnitTest::Test_comm121st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm121st, ReferenceImplementations::comm121st, "comm121st");
}
bool UnitTest::Test_comm122st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm122st, ReferenceImplementations::comm122st, "comm122st");
}
bool UnitTest::Test_comm221st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm221st, ReferenceImplementations::comm221st, "comm221st");
}
bool UnitTest::Test_comm222_pp_hhst(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm222_pp_hhst, ReferenceImplementations::comm222_pp_hhst, "comm222_pp_hhst");
}
bool UnitTest::Test_comm222_phst(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm222_phst, ReferenceImplementations::comm222_phst, "comm222_phst");
}



/// IMSRG(3) commutators
bool UnitTest::Test_comm330ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm330ss, ReferenceImplementations::comm330ss, "comm330ss");
}

bool UnitTest::Test_comm331ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm331ss, ReferenceImplementations::comm331ss, "comm331ss");
}

bool UnitTest::Test_comm231ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm231ss, ReferenceImplementations::comm231ss, "comm231ss");
}

bool UnitTest::Test_comm132ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm132ss, ReferenceImplementations::comm132ss, "comm132ss");
}

bool UnitTest::Test_comm332_ppph_hhhpss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm332_ppph_hhhpss, ReferenceImplementations::comm332_ppph_hhhpss, "comm332_ppph_hhhpss");
}

bool UnitTest::Test_comm332_pphhss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm332_pphhss, ReferenceImplementations::comm332_pphhss, "comm332_pphhss");
}

bool UnitTest::Test_comm133ss(const Operator &X, const Operator &Y)
{
  //Operator Xmod = X;
  ////  Xmod.ThreeBody.Erase();
  //Operator Ymod = Y;
  //Ymod.ThreeBody.Erase();
  //return Test_against_ref_impl(Xmod, Ymod, Commutator::comm133ss, ReferenceImplementations::comm133ss, "comm133ss");
  return Test_against_ref_impl(X,Y,  Commutator::comm133ss,  ReferenceImplementations::comm133ss,  "comm133ss");
  //  return Test_against_ref_impl(X,Y,  ReferenceImplementations::comm133ss,  ReferenceImplementations::comm133ss,  "comm133ss");
}

bool UnitTest::Test_comm223ss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm223ss, ReferenceImplementations::comm223ss, "comm223ss");
  //  return Test_against_ref_impl(X,Y,  ReferenceImplementations::comm223ss,  ReferenceImplementations::comm223ss,  "comm223ss");
}

bool UnitTest::Test_comm233_pp_hhss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm233_pp_hhss, ReferenceImplementations::comm233_pp_hhss, "comm233_pp_hhss");
}

bool UnitTest::Test_comm333_ppp_hhhss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm333_ppp_hhhss, ReferenceImplementations::comm333_ppp_hhhss, "comm333_ppp_hhhss");
}

/// THIS IS SO SLOW....
bool UnitTest::Test_comm233_phss(const Operator &X, const Operator &Y)
{
  //  return Test_against_ref_impl(X,Y,  Commutator::comm233_phss,  ReferenceImplementations::comm233_phss,  "comm233_phss");
  return Test_against_ref_impl(X, Y, Commutator::comm233_phss, Commutator::comm233_phss, "comm233_phss");
}

bool UnitTest::Test_comm232ss(const Operator &X, const Operator &Y)
{
  //  return Test_against_ref_impl(X,Y,  Commutator::comm232ss,  Commutator::comm232ss_slow,  "comm232ss");
  //  return Test_against_ref_impl(X,Y,  Commutator::comm232ss,  Commutator::comm232ss_debug,  "comm232ss");
  //  return Test_against_ref_impl(X,Y,  Commutator::comm232ss,  ReferenceImplementations::comm232ss,  "comm232ss");
  return Test_against_ref_impl(X, Y, Commutator::comm232ss, ReferenceImplementations::comm232ss, "comm232ss");
  //  return Test_against_ref_impl(X,Y,  Commutator::comm232ss,  Commutator::comm232ss_srs_optimized_old,  "comm232ss");
  // return Test_against_ref_impl(X,Y,  Commutator::comm232ss,  Commutator::comm232ss_expand_new, "comm232ss");
}

/// INCREDIBLY SLOW...
bool UnitTest::Test_comm333_pph_hhpss(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm333_pph_hhpss, ReferenceImplementations::comm333_pph_hhpss, "comm333_pph_hhpss");
}

bool UnitTest::Test_comm331st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm331st, ReferenceImplementations::comm331st, "comm331st");
}
bool UnitTest::Test_comm231st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm231st, ReferenceImplementations::comm231st, "comm231st");
}
bool UnitTest::Test_comm232st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm232st, ReferenceImplementations::comm232st, "comm232st");
}
bool UnitTest::Test_comm132st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm132st, ReferenceImplementations::comm132st, "comm132st");
}
bool UnitTest::Test_comm223st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm223st, ReferenceImplementations::comm223st, "comm223st");
}
bool UnitTest::Test_comm133st(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm133st, ReferenceImplementations::comm133st, "comm133st");
}

bool UnitTest::Test_comm332_pphhst(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm332_pphhst, ReferenceImplementations::comm332_pphhst, "comm332_pphhst");
}
bool UnitTest::Test_comm332_ppph_hhhpst(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm332_ppph_hhhpst, ReferenceImplementations::comm332_ppph_hhhpst, "comm332_ppph_hhhpst");
}
bool UnitTest::Test_comm233_pp_hhst(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm233_pp_hhst, ReferenceImplementations::comm233_pp_hhst, "comm233_pp_hhst");
}
bool UnitTest::Test_comm233_phst(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm233_phst, ReferenceImplementations::comm233_phst, "comm233_phst");
}
bool UnitTest::Test_comm333_ppp_hhhst(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm333_ppp_hhhst, ReferenceImplementations::comm333_ppp_hhhst, "comm333_ppp_hhhst");
}
bool UnitTest::Test_comm333_pph_hhpst(const Operator &X, const Operator &Y)
{
  return Test_against_ref_impl(X, Y, Commutator::comm333_pph_hhpst, ReferenceImplementations::comm333_pph_hhpst, "comm333_pph_hhpst");
}




    void comm332_pphhst(const Operator &X, const Operator &Y, Operator &Z);       // PASS the unit test
    void comm332_ppph_hhhpst(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test
    void comm233_pp_hhst(const Operator &X, const Operator &Y, Operator &Z);      // PASS the unit test
    void comm233_phst(const Operator &X, const Operator &Y, Operator &Z);         // PASS the unit test
    void comm333_ppp_hhhst(const Operator &X, const Operator &Y, Operator &Z);    // PASS the unit test
    void comm333_pph_hhpst(const Operator &X, const Operator &Y, Operator &Z);    // PASS the unit test


/// M-Scheme Formula:
///
/// Z0 =  sum_ab na(1-nb) (Xab * Yba - Yab * Xba)
///
bool UnitTest::Mscheme_Test_comm110ss(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

  Commutator::comm110ss(X, Y, Z_J);

  double Z0_m = 0;
  for (auto a : X.modelspace->all_orbits)
  {
    Orbit &oa = X.modelspace->GetOrbit(a);
    for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
    {
      double na = oa.occ;
      for (auto b : X.modelspace->all_orbits)
      {
        Orbit &ob = X.modelspace->GetOrbit(b);
        double nb = ob.occ;
        for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
        {
          double Xab = GetMschemeMatrixElement_1b(X, a, ma, b, mb);
          double Xba = GetMschemeMatrixElement_1b(X, b, mb, a, ma);
          double Yba = GetMschemeMatrixElement_1b(Y, b, mb, a, ma);
          double Yab = GetMschemeMatrixElement_1b(Y, a, ma, b, mb);

          Z0_m += na * (1 - nb) * (Xab * Yba - Yab * Xba);
        }
      }
    }
  }

  double summed_error = Z0_m - Z_J.ZeroBody;
  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << Z0_m << " " << Z_J.ZeroBody
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
///
/// Z0 = 1/4 sum_abcd nanb (1-nc)(1-nd) ( Xabcd * Ycdab - Yabcd * Xcdab )
///
bool UnitTest::Mscheme_Test_comm220ss(const Operator &X, const Operator &Y)
{
  Operator Z_J(Y);
  Z_J.Erase();

  Commutator::comm220ss(X, Y, Z_J);

  double Z0_m = 0;
  for (auto a : X.modelspace->all_orbits)
  {
    Orbit &oa = X.modelspace->GetOrbit(a);
    double na = oa.occ;
    for (auto b : X.modelspace->all_orbits)
    {
      Orbit &ob = X.modelspace->GetOrbit(b);
      double nb = ob.occ;
      for (auto c : X.modelspace->all_orbits)
      {
        Orbit &oc = X.modelspace->GetOrbit(c);
        double nc = oc.occ;
        for (auto d : X.modelspace->all_orbits)
        {
          Orbit &od = X.modelspace->GetOrbit(d);
          if ((oa.l + ob.l + oc.l + od.l) % 2 > 0)
            continue;
          if ((oa.tz2 + ob.tz2) != (oc.tz2 + od.tz2))
            continue;
          double nd = od.occ;

          for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
          {
            for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
            {
              for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
              {
                int md = ma + mb - mc;

                double Xabcd = GetMschemeMatrixElement_2b(X, a, ma, b, mb, c, mc, d, md);
                double Yabcd = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, c, mc, d, md);
                double Xcdab = GetMschemeMatrixElement_2b(X, c, mc, d, md, a, ma, b, mb);
                double Ycdab = GetMschemeMatrixElement_2b(Y, c, mc, d, md, a, ma, b, mb);

                Z0_m += (1. / 4) * na * nb * (1.0 - nc) * (1.0 - nd) * (Xabcd * Ycdab - Yabcd * Xcdab);

              } // for mc
            } // for mb
          } // for ma
        } // for d
      } // for c
    } // for b
  } // for a

  double summed_error = Z0_m - Z_J.ZeroBody;
  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << Z0_m << " " << Z_J.ZeroBody
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
///
///  Zij = sum_a (Xia * Yaj - Yia * Xaj)
///
// bool UnitTest::Test_comm111ss( const Operator& X, const Operator& Y )
bool UnitTest::Mscheme_Test_comm111ss(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

  Operator Xcpy(X);
  Operator Ycpy(Y);
  Ycpy.MakeReduced();
  Commutator::comm111st(Xcpy, Ycpy, Z_J);
  Z_J.MakeNotReduced();

  //  Commutator::comm111ss( X, Y, Z_J);
  //  Commutator::comm111st( X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      Orbit &oj = X.modelspace->GetOrbit(j);
      int mj = oj.j2;
      double Zm_ij = 0;
      if ((oi.j2 == oj.j2) and (oi.l == oj.l) and (oi.tz2 == oj.tz2))
      {
        for (auto a : X.modelspace->all_orbits)
        {
          Orbit &oa = X.modelspace->GetOrbit(a);
          for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
          {
            double Xia = GetMschemeMatrixElement_1b(X, i, mi, a, ma);
            double Yia = GetMschemeMatrixElement_1b(Y, i, mi, a, ma);
            double Xaj = GetMschemeMatrixElement_1b(X, a, ma, j, mj);
            double Yaj = GetMschemeMatrixElement_1b(Y, a, ma, j, mj);

            Zm_ij += Xia * Yaj - Yia * Xaj;
          }
        }
        double ZJ_ij = Z_J.OneBody(i, j);
        double err = Zm_ij - ZJ_ij;
        if (std::abs(err) > 1e-6)
        {
          std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij
                    << "   ZJ_ij = " << Z_J.OneBody(i, j) << "   err = " << err << std::endl;
        }
        summed_error += err * err;
        sum_m += Zm_ij * Zm_ij;
        sum_J += ZJ_ij * ZJ_ij;
      }
    }
  }

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.OneBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 1b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
/// Zij = [X1b,Y2b] - [Y1b,X2b]
//
/// Zij = sum_ab na*(1-nb) ( Xab * Ybiaj - Yaibj * Xba
///                        - Yab * Xbiaj + Xaibj * Yba )
///
// bool UnitTest::Test_comm121ss( const Operator& X, const Operator& Y)
bool UnitTest::Mscheme_Test_comm121ss(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

//  Operator Xcpy(X);
//  Operator Ycpy(Y);
//  Ycpy.MakeReduced();
//  Commutator::comm121ss(Xcpy, Ycpy, Z_J);
  Commutator::comm121ss(X,Y, Z_J);
//  Z_J.MakeNotReduced();


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      Orbit &oj = X.modelspace->GetOrbit(j);
      int mj = oj.j2;
      double Zm_ij = 0;
      if ((oi.j2 == oj.j2) and (oi.l == oj.l) and (oi.tz2 == oj.tz2))
      {
        for (auto a : X.modelspace->all_orbits)
        {
          Orbit &oa = X.modelspace->GetOrbit(a);
          double na = oa.occ;
          for (auto b : X.modelspace->all_orbits)
          {
            Orbit &ob = X.modelspace->GetOrbit(b);
            double nb = ob.occ;
            for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
            {
              for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
              {
                double Xab = GetMschemeMatrixElement_1b(X, a, ma, b, mb);
                double Yab = GetMschemeMatrixElement_1b(Y, a, ma, b, mb);
                double Xba = GetMschemeMatrixElement_1b(X, b, mb, a, ma);
                double Yba = GetMschemeMatrixElement_1b(Y, b, mb, a, ma);
                double Xbiaj = GetMschemeMatrixElement_2b(X, b, mb, i, mi, a, ma, j, mj);
                double Ybiaj = GetMschemeMatrixElement_2b(Y, b, mb, i, mi, a, ma, j, mj);
                double Xaibj = GetMschemeMatrixElement_2b(X, a, ma, i, mi, b, mb, j, mj);
                double Yaibj = GetMschemeMatrixElement_2b(Y, a, ma, i, mi, b, mb, j, mj);

                Zm_ij += na * (1 - nb) * (Xab * Ybiaj - Yaibj * Xba - Yab * Xbiaj + Xaibj * Yba);

              } // for mb
            } // for ma
          } // for b
        } // for a
      } // if ji=jj etc.
      double ZJ_ij = Z_J.OneBody(i, j);
      double err = Zm_ij - ZJ_ij;
      if (std::abs(err) > 1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij
                  << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl;
      }
      summed_error += err * err;
      sum_m += Zm_ij * Zm_ij;
      sum_J += ZJ_ij * ZJ_ij;
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.OneBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 1b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}


bool UnitTest::Mscheme_Test_comm121st(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();
  Operator Z_ref(Z_J);

  Commutator::comm121st(X, Y, Z_J);
  ReferenceImplementations::comm121st(X, Y, Z_ref);
  double lambda = Y.GetJRank();
//  std::cout << "lambda = " << lambda << std::endl;
//  std::cout << "Are they reduced ? " << X.IsReduced() << " " << Y.IsReduced() << " " << Z_J.IsReduced() << std::endl;
//  std::cout << "Hermitian? " << X.IsHermitian() << " " << Y.IsHermitian() << " " << Z_J.IsHermitian() << std::endl;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      Orbit &oj = X.modelspace->GetOrbit(j);
      int mj = oj.j2;
      int mu = (mi-mj)/2;
      double Zm_ij = 0;
//      if ((oi.j2 == oj.j2) and (oi.l == oj.l) and (oi.tz2 == oj.tz2))
      if (true)
      {
        for (auto a : X.modelspace->all_orbits)
        {
          Orbit &oa = X.modelspace->GetOrbit(a);
          double na = oa.occ;
          for (auto b : X.modelspace->all_orbits)
          {
            Orbit &ob = X.modelspace->GetOrbit(b);
            double nb = ob.occ;
            for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
            {
              for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
              {
                double Xab = GetMschemeMatrixElement_1b(X, a, ma, b, mb);
                double Yab = GetMschemeMatrixElement_1b(Y, a, ma, b, mb);
                double Xba = GetMschemeMatrixElement_1b(X, b, mb, a, ma);
                double Yba = GetMschemeMatrixElement_1b(Y, b, mb, a, ma);
                double Xbiaj = GetMschemeMatrixElement_2b(X, b, mb, i, mi, a, ma, j, mj);
                double Ybiaj = GetMschemeMatrixElement_2b(Y, b, mb, i, mi, a, ma, j, mj);
                double Xaibj = GetMschemeMatrixElement_2b(X, a, ma, i, mi, b, mb, j, mj);
                double Yaibj = GetMschemeMatrixElement_2b(Y, a, ma, i, mi, b, mb, j, mj);

                Zm_ij += na * (1 - nb) * (Xab * Ybiaj - Yaibj * Xba - Yab * Xbiaj + Xaibj * Yba);

//                if ( (  (i==2 and j==4) or (i==4 and j==2)) and (a==0 and b==8))
//                {
//                   std::cout << "## ij " << i << " " << j << "   ab " << a << " " << b << " : " << Xab << " " << Xba << "  " << Yab << " "<< Yba
//                             << "   " << Xbiaj << " " << Xaibj << "   " << Ybiaj << " " << Yaibj << "      =  "
//                             << na * (1 - nb) * (Xab * Ybiaj - Yaibj * Xba - Yab * Xbiaj + Xaibj * Yba)  << "  =>  " << Zm_ij
//                             << "   projections " << mi << " " << mj << " " << ma << " " << mb << "  mu = " << mu  << std::endl;
//                }

              } // for mb
            } // for ma
          } // for b
        } // for a
      } // if ji=jj etc.
//      double ZJ_ij = Z_J.OneBody(i, j);
      double ZJ_ij = GetMschemeMatrixElement_1b(Z_J, i, mi, j, mj) ;
      double Zr_ij = GetMschemeMatrixElement_1b(Z_ref, i, mi, j, mj) ;
      // The J coupled matrix element should be reduced.
//      double cg = AngMom::CG( oj.j2/2. , mj/2., lambda,mu,  oi.j2/2., mi/2.);
//      std::cout << "clebsch gordan < " << oj.j2/2. << " "<<  mj/2. << " " << lambda << " " << mu << " | " << oi.j2/2. << mi/2. << " > = " << cg << std::endl;
//      ZJ_ij *= cg/sqrt(oi.j2+1);
//      double err = Zm_ij - ZJ_ij;
      double err = Zm_ij - Zr_ij;
      if (std::abs(err) > 1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij
                  << "   ZJ_ij = " << Zr_ij << "   err = " << err <<  "  optimized: " << ZJ_ij << std::endl;
      }
//      else
//      {
//        std::cout << "     OK for i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij<< "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl;
//      }
      summed_error += err * err;
      sum_m += Zm_ij * Zm_ij;
      sum_J += ZJ_ij * ZJ_ij;
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.OneBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 1b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Zij = 1/2 sum_abc [ na*nb*(1-nc) +(1-na)*(1-nb)*nc ] * (Xciab * Yabcj - Yciab * Xabcj)
//
//
bool UnitTest::Mscheme_Test_comm221ss(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

  Operator Xcpy(X);
  Operator Ycpy(Y);
  Ycpy.MakeReduced();
  Commutator::comm222_pp_hh_221st(Xcpy, Ycpy, Z_J);
  Z_J.MakeNotReduced();
  //  Commutator::comm222_pp_hh_221ss( X, Y, Z_J ) ;
  //  Commutator::comm221ss( X, Y, Z_J);
  //  Commutator::comm222_pp_hh_221st( X, Y, Z_J);
  Z_J.EraseTwoBody();
  //  ReferenceImplementations::comm221ss( X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      Orbit &oj = X.modelspace->GetOrbit(j);
      int mj = oj.j2;
      double Zm_ij = 0;
      if ((oi.j2 == oj.j2) and (oi.l == oj.l) and (oi.tz2 == oj.tz2))
      {
        for (auto a : X.modelspace->all_orbits)
        {
          Orbit &oa = X.modelspace->GetOrbit(a);
          double na = oa.occ;
          for (auto b : X.modelspace->all_orbits)
          {

            Orbit &ob = X.modelspace->GetOrbit(b);
            double nb = ob.occ;
            for (auto c : X.modelspace->all_orbits)
            {
              Orbit &oc = X.modelspace->GetOrbit(c);
              double nc = oc.occ;
              if (std::abs(na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc) < 1e-6)
                continue;
              if ((oi.l + oc.l + oa.l + ob.l) % 2 > 0)
                continue;
              if ((oi.tz2 + oc.tz2) != (oa.tz2 + ob.tz2))
                continue;

              for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
              {
                for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                {

                  int mc = ma + mb - mi;
                  if (std::abs(mc) > oc.j2)
                    continue;
                  double Xciab = GetMschemeMatrixElement_2b(X, c, mc, i, mi, a, ma, b, mb);
                  double Yciab = GetMschemeMatrixElement_2b(Y, c, mc, i, mi, a, ma, b, mb);
                  double Xabcj = GetMschemeMatrixElement_2b(X, a, ma, b, mb, c, mc, j, mj);
                  double Yabcj = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, c, mc, j, mj);

                  Zm_ij += (1. / 2) * (na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc) * (Xciab * Yabcj - Yciab * Xabcj);

                } // for ma
              } // for mb
            } // for c
          } // for b
        } // for a

      } // if ji=jj etc.

      double ZJ_ij = Z_J.OneBody(i, j);
      double err = Zm_ij - ZJ_ij;
      if (std::abs(err) > 1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij
                  << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl;
      }
      summed_error += err * err;
      sum_m += Zm_ij * Zm_ij;
      sum_J += ZJ_ij * ZJ_ij;
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.OneBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 1b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Zij = 1/2 sum_abc [ na*nb*(1-nc) +(1-na)*(1-nb)*nc ] * (Xciab * Yabcj - Yciab * Xabcj)
//
//
bool UnitTest::Mscheme_Test_comm221st(const Operator &X, const Operator &Y)
{


  Operator Z_J(Y);
  Z_J.Erase();
  Operator Z_ref(Z_J);

  Commutator::comm222_pp_hh_221st(X, Y, Z_J);
  ReferenceImplementations::comm221st(X, Y, Z_ref);
  double lambda = Y.GetJRank();




  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      Orbit &oj = X.modelspace->GetOrbit(j);
      int mj = oj.j2;
      int mu = mi - mj;
      double Zm_ij = 0;
//      if ((oi.j2 == oj.j2) and (oi.l == oj.l) and (oi.tz2 == oj.tz2))
      if (true)
      {
        for (auto a : X.modelspace->all_orbits)
        {
          Orbit &oa = X.modelspace->GetOrbit(a);
          double na = oa.occ;
          for (auto b : X.modelspace->all_orbits)
          {

            Orbit &ob = X.modelspace->GetOrbit(b);
            double nb = ob.occ;
            for (auto c : X.modelspace->all_orbits)
            {
              Orbit &oc = X.modelspace->GetOrbit(c);
              double nc = oc.occ;
              if (std::abs(na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc) < 1e-6)
                continue;
//              if ((oi.l + oc.l + oa.l + ob.l) % 2 > 0)
//                continue;
//              if ((oi.tz2 + oc.tz2) - (oa.tz2 + ob.tz2))
//                continue;

              for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
              {
                for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                {

                  double Xciab=0;
                  double Yciab=0;
                  double Xabcj=0;
                  double Yabcj=0;
//                  int mc = ma + mb - mi;
                  int mc = ma + mb - mi;
                  if (std::abs(mc) <= oc.j2)
                  {
                     Xciab = GetMschemeMatrixElement_2b(X, c, mc, i, mi, a, ma, b, mb);
                     Yabcj = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, c, mc, j, mj);
                  }
                  mc = ma + mb - mj;
                  if (std::abs(mc) <= oc.j2)
                  {
                     Xabcj = GetMschemeMatrixElement_2b(X, a, ma, b, mb, c, mc, j, mj);
                     Yciab = GetMschemeMatrixElement_2b(Y, c, mc, i, mi, a, ma, b, mb);
                  }

                  Zm_ij += (1. / 2) * (na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc) * (Xciab * Yabcj - Yciab * Xabcj);

                } // for ma
              } // for mb
            } // for c
          } // for b
        } // for a

      } // if ji=jj etc.

//      double ZJ_ij = Z_J.OneBody(i, j);

      double ZJ_ij = GetMschemeMatrixElement_1b(Z_J, i, mi, j, mj) ;
      double Zr_ij = GetMschemeMatrixElement_1b(Z_ref, i, mi, j, mj) ;

//      double err = Zm_ij - ZJ_ij;
      double err = Zm_ij - Zr_ij;
      if (std::abs(err) > 1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j << "   Zm_ij = " << Zm_ij
                  << "   ZJ_ij = " << ZJ_ij << "  Zr_ij = " << Zr_ij  << "   err = " << err << std::endl;
      }
      summed_error += err * err;
      sum_m += Zm_ij * Zm_ij;
      sum_J += ZJ_ij * ZJ_ij;
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.OneBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 1b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}




/// M-Scheme Formula:
//
// Zijkl = sum_a   (Xia * Yajkl + Xja * Yiakl - Yijal * Xak - Yijka * Xal)
//               - (Yia * Xajkl + Yja * Xiakl - Xijal * Yak - Xijka * Yal)
//
bool UnitTest::Mscheme_Test_comm122ss(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

  Operator Xcpy(X);
  Operator Ycpy(Y);
  Ycpy.MakeReduced();
  Commutator::comm122st(Xcpy, Ycpy, Z_J);
  Z_J.MakeNotReduced();
  //  Commutator::comm122ss( X, Y, Z_J);
  //  Commutator::comm122st( X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      if (not(oi.l == oj.l and oi.j2 == oj.j2 and oi.tz2 == oj.tz2))
        continue;
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if (not(ok.l == ol.l and ok.j2 == ol.j2 and ok.tz2 == ol.tz2))
            continue;
          if ((oi.l + oj.l + ok.l + ol.l) % 2 > 0)
            continue;
          if ((oi.tz2 + oj.tz2) != (ok.tz2 + ol.tz2))
            continue;

          for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
          {
            for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
            {

              int ml = mi + mj - mk;
              if (std::abs(ml) > ol.j2)
                continue;

              double Zm_ijkl = 0;

              for (auto a : X.modelspace->all_orbits)
              {
                Orbit &oa = X.modelspace->GetOrbit(a);

                for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                {
                  double Xia = GetMschemeMatrixElement_1b(X, i, mi, a, ma);
                  double Xja = GetMschemeMatrixElement_1b(X, j, mj, a, ma);
                  double Xak = GetMschemeMatrixElement_1b(X, a, ma, k, mk);
                  double Xal = GetMschemeMatrixElement_1b(X, a, ma, l, ml);
                  double Yia = GetMschemeMatrixElement_1b(Y, i, mi, a, ma);
                  double Yja = GetMschemeMatrixElement_1b(Y, j, mj, a, ma);
                  double Yak = GetMschemeMatrixElement_1b(Y, a, ma, k, mk);
                  double Yal = GetMschemeMatrixElement_1b(Y, a, ma, l, ml);
                  double Xajkl = GetMschemeMatrixElement_2b(X, a, ma, j, mj, k, mk, l, ml);
                  double Xiakl = GetMschemeMatrixElement_2b(X, i, mi, a, ma, k, mk, l, ml);
                  double Xijal = GetMschemeMatrixElement_2b(X, i, mi, j, mj, a, ma, l, ml);
                  double Xijka = GetMschemeMatrixElement_2b(X, i, mi, j, mj, k, mk, a, ma);
                  double Yajkl = GetMschemeMatrixElement_2b(Y, a, ma, j, mj, k, mk, l, ml);
                  double Yiakl = GetMschemeMatrixElement_2b(Y, i, mi, a, ma, k, mk, l, ml);
                  double Yijal = GetMschemeMatrixElement_2b(Y, i, mi, j, mj, a, ma, l, ml);
                  double Yijka = GetMschemeMatrixElement_2b(Y, i, mi, j, mj, k, mk, a, ma);

                  Zm_ijkl += (Xia * Yajkl + Xja * Yiakl - Yijal * Xak - Yijka * Xal) - (Yia * Xajkl + Yja * Xiakl - Xijal * Yak - Xijka * Yal);

                } // for ma
              } // for a

              double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
              double err = Zm_ijkl - ZJ_ijkl;
              if (std::abs(err) > 1e-6)
              {
                std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                          << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
              }
              summed_error += err * err;
              sum_m += Zm_ijkl * Zm_ijkl;
              sum_J += ZJ_ijkl * ZJ_ijkl;

            } // for mk
          } // for mj
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Zijkl = 1/2 sum_ab [(1-na)*(1-nb) - na*nb ] * ( Xijab * Yabkl - Yijab * Xabkl )
//
bool UnitTest::Mscheme_Test_comm222_pp_hhss(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

    Commutator::comm222_pp_hhss( X, Y, Z_J);

//  Operator Xcpy(X);
//  Operator Ycpy(Y);
//  Ycpy.MakeReduced();
//  Commutator::comm222_pp_hh_221ss(Xcpy, Ycpy, Z_J);
//  Z_J.MakeNotReduced();
  //  Commutator::comm222_pp_hh_221st( X, Y, Z_J);
  Z_J.EraseOneBody();
  //  Commutator::comm222_pp_hh_221ss( X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if ((oi.l + oj.l + ok.l + ol.l) % 2 > 0)
            continue;
          if ((oi.tz2 + oj.tz2) != (ok.tz2 + ol.tz2))
            continue;

          for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
          {
            for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
            {
              if (std::abs(mk) > ol.j2)
                continue;
              int ml = mi + mj - mk;

              double Zm_ijkl = 0;
              for (auto a : X.modelspace->all_orbits)
              {
                Orbit &oa = X.modelspace->GetOrbit(a);
                double na = oa.occ;
                for (auto b : X.modelspace->all_orbits)
                {
                  Orbit &ob = X.modelspace->GetOrbit(b);
                  double nb = ob.occ;
                  if (std::abs(((1 - na) * (1 - nb) - na * nb)) < 1e-6)
                    continue;

                  if ((oi.l + oj.l + oa.l + ob.l) % 2 > 0)
                    continue;
                  if ((oi.tz2 + oj.tz2) != (oa.tz2 + ob.tz2))
                    continue;

                  for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                  {
                    int mb = mi + mj - ma;
                    if (std::abs(mb) > ob.j2)
                      continue;

                    double Xijab = GetMschemeMatrixElement_2b(X, i, mi, j, mj, a, ma, b, mb);
                    double Xabkl = GetMschemeMatrixElement_2b(X, a, ma, b, mb, k, mk, l, ml);
                    double Yijab = GetMschemeMatrixElement_2b(Y, i, mi, j, mj, a, ma, b, mb);
                    double Yabkl = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, k, mk, l, ml);
                    Zm_ijkl += (1. / 2) * ((1 - na) * (1 - nb) - na * nb) * (Xijab * Yabkl - Yijab * Xabkl);
                  } // for ma
                } // for b
              } // for a

              double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
              double err = Zm_ijkl - ZJ_ijkl;
              if (std::abs(err) > 1e-6)
              {
                std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                          << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
              }
              summed_error += err * err;
              sum_m += Zm_ijkl * Zm_ijkl;
              sum_J += ZJ_ijkl * ZJ_ijkl;

            } // for mk
          } // for mj
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}


/// M-Scheme Formula:
//
// Zijkl = 1/2 sum_ab [(1-na)*(1-nb) - na*nb ] * ( Xijab * Yabkl - Yijab * Xabkl )
//
bool UnitTest::Mscheme_Test_comm222_pp_hhst(const Operator &X, const Operator &Y)
{


  double lambda = Y.GetJRank();

  int parityZ = (X.GetParity() + Y.GetParity()) % 2;
  int TzZ = X.GetTRank() + Y.GetTRank();
  int JrankZ = Y.GetJRank();

  Operator Z_J(*(Y.modelspace), JrankZ, TzZ, parityZ, 2);
  Operator Z_ref(Z_J);


  Commutator::comm222_phst(X, Y, Z_J);
  ReferenceImplementations::comm222_phst(X, Y, Z_ref);


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if ((oi.l + oj.l + ok.l + ol.l + Z_J.GetParity()) % 2 > 0)
            continue;
          if (std::abs((oi.tz2 + oj.tz2) - (ok.tz2 + ol.tz2)) != Z_J.GetTRank() * 2)
            continue;

          for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
          {
            for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
            {
              if (std::abs(mk) > ol.j2)
                continue;
              int ml = mi + mj - mk;

              double Zm_ijkl = 0;
              for (auto a : X.modelspace->all_orbits)
              {
                Orbit &oa = X.modelspace->GetOrbit(a);
                double na = oa.occ;
                for (auto b : X.modelspace->all_orbits)
                {
                  Orbit &ob = X.modelspace->GetOrbit(b);
                  double nb = ob.occ;
                  if (std::abs(((1 - na) * (1 - nb) - na * nb)) < 1e-6)
                    continue;

                  if ((oi.l + oj.l + oa.l + ob.l) % 2 > 0)
                    continue;
                  if ((oi.tz2 + oj.tz2) != (oa.tz2 + ob.tz2))
                    continue;

                  for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                  {
                    int mb = mi + mj - ma;
                    if (std::abs(mb) > ob.j2)
                      continue;

                    double Xijab = GetMschemeMatrixElement_2b(X, i, mi, j, mj, a, ma, b, mb);
                    double Xabkl = GetMschemeMatrixElement_2b(X, a, ma, b, mb, k, mk, l, ml);
                    double Yijab = GetMschemeMatrixElement_2b(Y, i, mi, j, mj, a, ma, b, mb);
                    double Yabkl = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, k, mk, l, ml);
                    Zm_ijkl += (1. / 2) * ((1 - na) * (1 - nb) - na * nb) * (Xijab * Yabkl - Yijab * Xabkl);
                  } // for ma
                } // for b
              } // for a

              double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
              double err = Zm_ijkl - ZJ_ijkl;
              if (std::abs(err) > 1e-6)
              {
                std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                          << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
              }
              summed_error += err * err;
              sum_m += Zm_ijkl * Zm_ijkl;
              sum_J += ZJ_ijkl * ZJ_ijkl;

            } // for mk
          } // for mj
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}



/// M-Scheme Formula:
//
// Zijkl = sum_ab [na*(1-nb) - nb*(1-na)] * (1-Pij)(1-Pkl) ( Xaibk Ybjal  )
//
bool UnitTest::Mscheme_Test_comm222_phss(const Operator &X, const Operator &Y)
{

  int parityZ = (X.GetParity() + Y.GetParity()) % 2;
  int TzZ = X.GetTRank() + Y.GetTRank();
  int JrankZ = 0;

  Operator Z_J(*(Y.modelspace), JrankZ, TzZ, parityZ, 2);
  Operator Z_ref(*(Y.modelspace), JrankZ, TzZ, parityZ, 2);

  Commutator::comm222_phss(X, Y, Z_J);
  ReferenceImplementations::comm222_phss( X, Y, Z_ref);


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if ((oi.l + oj.l + ok.l + ol.l + Z_J.GetParity()) % 2 > 0)
            continue;
          if (std::abs((oi.tz2 + oj.tz2) - (ok.tz2 + ol.tz2)) != Z_J.GetTRank() * 2)
            continue;

          for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
          {
            for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
            {

              int ml = mi + mj - mk;
              if (std::abs(ml) > ol.j2)
                continue;

              double Zm_ijkl = 0;
              for (auto a : X.modelspace->all_orbits)
              {
                Orbit &oa = X.modelspace->GetOrbit(a);
                double na = oa.occ;
                for (auto b : X.modelspace->all_orbits)
                {
                  Orbit &ob = X.modelspace->GetOrbit(b);
                  double nb = ob.occ;
                  double Zsave = Zm_ijkl;
                  //                if ( std::abs(  na*(1-nb) - nb*(1-na) ) < 1e-6) continue;
                  if (std::abs(na - nb) < 1e-6)
                    continue;

                  for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                  {
                    for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                    {

                      double Yaibk = GetMschemeMatrixElement_2b(Y, a, ma, i, mi, b, mb, k, mk);
                      double Xbjal = GetMschemeMatrixElement_2b(X, b, mb, j, mj, a, ma, l, ml);

                      // Pij
                      double Yajbk = GetMschemeMatrixElement_2b(Y, a, ma, j, mj, b, mb, k, mk);
                      double Xbial = GetMschemeMatrixElement_2b(X, b, mb, i, mi, a, ma, l, ml);

                      // Pkl
                      double Yaibl = GetMschemeMatrixElement_2b(Y, a, ma, i, mi, b, mb, l, ml);
                      double Xbjak = GetMschemeMatrixElement_2b(X, b, mb, j, mj, a, ma, k, mk);

                      // PijPkl
                      double Yajbl = GetMschemeMatrixElement_2b(Y, a, ma, j, mj, b, mb, l, ml);
                      double Xbiak = GetMschemeMatrixElement_2b(X, b, mb, i, mi, a, ma, k, mk);

                      Zm_ijkl -= (na - nb) * (Yaibk * Xbjal - Yajbk * Xbial - Yaibl * Xbjak + Yajbl * Xbiak);



                    } // for mb
                  } // for ma
                } // for b
              } // for a

              //             double ZJ_ijkl = GetMschemeMatrixElement_2b( Z_J, i,mi, j,mj, k,mk, l,ml ) ;
              double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
              double Zr_ijkl = GetMschemeMatrixElement_2b(Z_ref, i, mi, j, mj, k, mk, l, ml);
              double err = Zm_ijkl - Zr_ijkl;
//              double err = Zm_ijkl - ZJ_ijkl;
              if (std::abs(err) > 1e-6)
              {
                std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                          << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "  Zr_ijkl = " << Zr_ijkl <<  "   err = " << err << std::endl;
              }
              summed_error += err * err;
              sum_m += Zm_ijkl * Zm_ijkl;
              sum_J += ZJ_ijkl * ZJ_ijkl;

            } // for mk
          } // for mj
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}


bool UnitTest::Mscheme_Test_comm222_phst(const Operator &X, const Operator &Y)
{

  int parityZ = (X.GetParity() + Y.GetParity()) % 2;
  int TzZ = X.GetTRank() + Y.GetTRank();
  int JrankZ =  Y.GetJRank();

  Operator Z_J(*(Y.modelspace), JrankZ, TzZ, parityZ, 2);


//  Commutator::comm222_phst(X, Y, Z_J);
  ReferenceImplementations::comm222_phst( X, Y, Z_J);



  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if ((oi.l + oj.l + ok.l + ol.l + Z_J.GetParity()) % 2 > 0)
            continue;
          if (std::abs((oi.tz2 + oj.tz2) - (ok.tz2 + ol.tz2)) != Z_J.GetTRank() * 2)
            continue;

          for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
          {
            for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
            {

              int ml = mi + mj - mk;
              if (std::abs(ml) > ol.j2)
                continue;

              double Zm_ijkl = 0;
              for (auto a : X.modelspace->all_orbits)
              {
                Orbit &oa = X.modelspace->GetOrbit(a);
                double na = oa.occ;
                for (auto b : X.modelspace->all_orbits)
                {
                  Orbit &ob = X.modelspace->GetOrbit(b);
                  double nb = ob.occ;
                  double Zsave = Zm_ijkl;
                  //                if ( std::abs(  na*(1-nb) - nb*(1-na) ) < 1e-6) continue;
                  if (std::abs(na - nb) < 1e-6)
                    continue;

                  for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                  {
                    for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                    {

                      double Xaibk = GetMschemeMatrixElement_2b(X, a, ma, i, mi, b, mb, k, mk);
                      double Ybjal = GetMschemeMatrixElement_2b(Y, b, mb, j, mj, a, ma, l, ml);

                      // Pij
                      double Xajbk = GetMschemeMatrixElement_2b(X, a, ma, j, mj, b, mb, k, mk);
                      double Ybial = GetMschemeMatrixElement_2b(Y, b, mb, i, mi, a, ma, l, ml);

                      // Pkl
                      double Xaibl = GetMschemeMatrixElement_2b(X, a, ma, i, mi, b, mb, l, ml);
                      double Ybjak = GetMschemeMatrixElement_2b(Y, b, mb, j, mj, a, ma, k, mk);

                      // PijPkl
                      double Xajbl = GetMschemeMatrixElement_2b(X, a, ma, j, mj, b, mb, l, ml);
                      double Ybiak = GetMschemeMatrixElement_2b(Y, b, mb, i, mi, a, ma, k, mk);

                      Zm_ijkl += (na - nb) * (Xaibk * Ybjal - Xajbk * Ybial - Xaibl * Ybjak + Xajbl * Ybiak);

                    } // for mb
                  } // for ma
                } // for b
              } // for a

              //             double ZJ_ijkl = GetMschemeMatrixElement_2b( Z_J, i,mi, j,mj, k,mk, l,ml ) ;
              double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
              double err = Zm_ijkl - ZJ_ijkl;
              if (std::abs(err) > 1e-6)
              {
                std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                          << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
              }
              summed_error += err * err;
              sum_m += Zm_ijkl * Zm_ijkl;
              sum_J += ZJ_ijkl * ZJ_ijkl;

            } // for mk
          } // for mj
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z0b = 1/36 sum_abcdef  na*nb*nc*(1-nd)(1-ne)(1-nf) ( X_abcdef * Y_defabc - Yabcdef * Xdefabc )
//
// bool UnitTest::Test_comm330ss( const Operator& X, const Operator& Y )
bool UnitTest::Mscheme_Test_comm330ss(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;
  Operator Z_J(Y);
  Z_J.Erase();

  Commutator::comm330ss(X, Y, Z_J);
  //  ReferenceImplementations::comm330ss( X, Y, Z_J);

  //  if ( Z_J.IsHermitian() )
  //     Z_J.Symmetrize();
  //  else if (Z_J.IsAntiHermitian() )
  //     Z_J.AntiSymmetrize();

  double Z0_m = 0;
  size_t norbits = X.modelspace->GetNumberOrbits();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Z0_m)
  for (size_t a = 0; a < norbits; a++)
  //  for (auto a : X.modelspace->all_orbits )
  {
    Orbit &oa = X.modelspace->GetOrbit(a);
    double na = oa.occ;
    for (auto b : X.modelspace->all_orbits)
    {
      Orbit &ob = X.modelspace->GetOrbit(b);
      double nb = ob.occ;
      for (auto c : X.modelspace->all_orbits)
      {
        Orbit &oc = X.modelspace->GetOrbit(c);
        double nc = oc.occ;
        if (std::abs(na * nb * nc) < 1e-9)
          continue;

        for (auto d : X.modelspace->all_orbits)
        {
          Orbit &od = X.modelspace->GetOrbit(d);
          double nd = od.occ;

          for (auto e : X.modelspace->all_orbits)
          {
            Orbit &oe = X.modelspace->GetOrbit(e);
            double ne = oe.occ;
            for (auto f : X.modelspace->all_orbits)
            {
              Orbit &of = X.modelspace->GetOrbit(f);
              double nf = of.occ;
              if (std::abs((1. - nd) * (1. - ne) * (1. - nf)) < 1e-9)
                continue;

              // check parity and isospin projection
              //        if ( (oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0 ) continue;
              //        if ( (oa.tz2+ob.tz2+oc.tz2) !=(od.tz2+oe.tz2+of.tz2) ) continue;
              if (((oa.l + ob.l + oc.l + od.l + oe.l + of.l + X.GetParity()) % 2 != 0) or ((oa.l + ob.l + oc.l + od.l + oe.l + of.l + Y.GetParity()) % 2 != 0))
                continue;
              if ((std::abs(oa.tz2 + ob.tz2 + oc.tz2 - od.tz2 - oe.tz2 - of.tz2) != 2 * X.GetTRank()) or (std::abs(oa.tz2 + ob.tz2 + oc.tz2 - od.tz2 - oe.tz2 - of.tz2) != 2 * Y.GetTRank()))
                continue;
              double dz = 0;

              std::cout << " abcdef " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;

              for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
              {
                for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                {
                  if ((a == b) and (ma == mb))
                    continue;
                  for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                  {
                    if ((a == c) and (ma == mc))
                      continue;
                    if ((b == c) and (mb == mc))
                      continue;

                    for (int md = -od.j2; md <= od.j2; md += 2)
                    {
                      for (int me = -oe.j2; me <= oe.j2; me += 2)
                      {
                        if ((d == e) and (md == me))
                          continue;
                        int mf = ma + mb + mc - me - md;
                        if (std::abs(mf) > of.j2)
                          continue;
                        if ((d == f) and (md == mf))
                          continue;
                        if ((e == f) and (me == mf))
                          continue;

                        double Xabcdef = GetMschemeMatrixElement_3b(X, a, ma, b, mb, c, mc, d, md, e, me, f, mf);
                        double Yabcdef = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, c, mc, d, md, e, me, f, mf);
                        double Xdefabc = GetMschemeMatrixElement_3b(X, d, md, e, me, f, mf, a, ma, b, mb, c, mc);
                        double Ydefabc = GetMschemeMatrixElement_3b(Y, d, md, e, me, f, mf, a, ma, b, mb, c, mc);

                        dz += (1. / 36) * na * nb * nc * (1 - nd) * (1 - ne) * (1 - nf) * (Xabcdef * Ydefabc - Yabcdef * Xdefabc);

                      } // for me
                    } // for md
                  } // for mc
                } // for mb
              } // for ma
              Z0_m += dz;
            } // for f
          } // for e
        } // for d
      } // for c
    } // for b
  } // for a

  double summed_error = Z0_m - Z_J.ZeroBody;
  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << Z0_m << " " << Z_J.ZeroBody
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ij = 1/12 sum_abcde (nanb n`c n`dn`e + n`an`b ncndne) (X_abicde Y_cdeabj - Y_abicde X_cdeabj)
//
// this is very slow...
bool UnitTest::Mscheme_Test_comm331ss(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

  Commutator::comm331ss(X, Y, Z_J);
  //  ReferenceImplementations::comm331ss( X, Y, Z_J);

  //  if ( Z_J.IsHermitian() )
  //     Z_J.Symmetrize();
  //  else if (Z_J.IsAntiHermitian() )
  //     Z_J.AntiSymmetrize( );

  std::cout << __func__ << "   Done with J-scheme commutator. Begin m-scheme version..." << std::endl;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : Z_J.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
    {
      //      Orbit& oj = X.modelspace->GetOrbit(j);
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      int mj = mi;
      double Zm_ij = 0;
      size_t norb = X.modelspace->GetNumberOrbits();
      //      for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ij)
      for (size_t a = 0; a < norb; a++)
      {
        Orbit &oa = X.modelspace->GetOrbit(a);
        double na = oa.occ;
        //        if ( na<1e-6 ) continue;
        for (auto b : X.modelspace->all_orbits)
        {
          Orbit &ob = X.modelspace->GetOrbit(b);
          double nb = ob.occ;
          //          if ( nb<1e-6 ) continue;
          for (auto c : X.modelspace->all_orbits)
          {
            Orbit &oc = X.modelspace->GetOrbit(c);
            double nc = oc.occ;
            //            if ( (1-nc)<1e-6 ) continue;
            for (auto d : X.modelspace->all_orbits)
            {
              Orbit &od = X.modelspace->GetOrbit(d);
              double nd = od.occ;
              //              if ( (1-nd)<1e-6 ) continue;
              for (auto e : X.modelspace->all_orbits)
              {
                Orbit &oe = X.modelspace->GetOrbit(e);
                //                if ( (oa.l+ob.l+oi.l + oc.l+od.l+oe.l)%2 > 0) continue;
                //                if ( (oa.tz2+ob.tz2+oi.tz2) != (oc.tz2+od.tz2+oe.tz2)) continue;

                if (((oa.l + ob.l + oi.l + oc.l + od.l + oe.l + X.GetParity()) % 2 != 0) and ((oa.l + ob.l + oi.l + oc.l + od.l + oe.l + Y.GetParity()) % 2 != 0))
                  continue;
                if ((std::abs(oa.tz2 + ob.tz2 + oi.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * X.GetTRank()) and (std::abs(oa.tz2 + ob.tz2 + oi.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * Y.GetTRank()))
                  continue;
                if (((oa.l + ob.l + oj.l + oc.l + od.l + oe.l + X.GetParity()) % 2 != 0) and ((oa.l + ob.l + oj.l + oc.l + od.l + oe.l + Y.GetParity()) % 2 != 0))
                  continue;
                if ((std::abs(oa.tz2 + ob.tz2 + oj.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * X.GetTRank()) and (std::abs(oa.tz2 + ob.tz2 + oj.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * Y.GetTRank()))
                  continue;

                double ne = oe.occ;
                //                if ( (1-ne)<1e-6 ) continue;
                // Mistake found by Matthias Heinz Oct 14 2022                double occfactor = na*nb*(1-nc)*(1-nd)*(1-ne);
                double occfactor = na * nb * (1 - nc) * (1 - nd) * (1 - ne) + (1 - na) * (1 - nb) * nc * nd * ne;
                if (std::abs(occfactor) < 1e-7)
                  continue;

                for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                {
                  for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                  {
                    for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                    {
                      for (int md = -od.j2; md <= od.j2; md += 2)
                      {
                        int me = ma + mb + mi - mc - md;
                        if (std::abs(me) > oe.j2)
                          continue;
                        double xabicde = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, c, mc, d, md, e, me);
                        double ycdeabj = GetMschemeMatrixElement_3b(Y, c, mc, d, md, e, me, a, ma, b, mb, j, mj);
                        double xcdeabj = GetMschemeMatrixElement_3b(X, c, mc, d, md, e, me, a, ma, b, mb, j, mj);
                        double yabicde = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, i, mi, c, mc, d, md, e, me);
                        Zm_ij += 1. / 12 * occfactor * (xabicde * ycdeabj - yabicde * xcdeabj);
                      } // for md
                    } // for mc
                  } // for mb
                } // for ma
              } // for e
            } // for d
          } // for c
        } // for b
      } // for a

      double ZJ_ij = Z_J.OneBody(i, j);
      double err = Zm_ij - ZJ_ij;
      if (std::abs(err) > 1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j
                  << "   zij = " << Zm_ij << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl;
      }
      summed_error += err * err;
      sum_m += Zm_ij * Zm_ij;
      sum_J += ZJ_ij * ZJ_ij;

    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
//  Z_ij = 1/4 sum_abcd (nanb n`cn`d) (X_abcd Y_cdiabj - Y_abicdj X_cdab)
//
bool UnitTest::Mscheme_Test_comm231ss(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;

  Operator Z_J(Y);
  Z_J.Erase();

  Commutator::comm231ss(X, Y, Z_J);
  // ReferenceImplementations::comm231ss( X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  std::cout << __func__ << "   Done with J-scheme commutator. Begin m-scheme version..." << std::endl;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;
  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : Z_J.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
    {
      //      Orbit& oj = X.modelspace->GetOrbit(j);
      int mj = mi;
      double Zm_ij = 0;

      size_t norb = X.modelspace->GetNumberOrbits();
      //      for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ij)
      for (size_t a = 0; a < norb; a++)
      {
        Orbit &oa = X.modelspace->GetOrbit(a);
        double na = oa.occ;
        if (na < 1e-6)
          continue;
        for (auto b : X.modelspace->all_orbits)
        {
          Orbit &ob = X.modelspace->GetOrbit(b);
          double nb = ob.occ;
          if (nb < 1e-6)
            continue;
          for (auto c : X.modelspace->all_orbits)
          {
            Orbit &oc = X.modelspace->GetOrbit(c);
            double nc = oc.occ;
            if ((1 - nc) < 1e-6)
              continue;
            for (auto d : X.modelspace->all_orbits)
            {
              Orbit &od = X.modelspace->GetOrbit(d);
              double nd = od.occ;
              if ((1 - nd) < 1e-6)
                continue;
              //              if ( (oa.l+ob.l+oc.l+od.l)%2 > 0 ) continue;
              //              if ( (oa.tz2+ob.tz2) !=(oc.tz2+od.tz2) ) continue;
              if (((oa.l + ob.l + oc.l + od.l + parityX) % 2 > 0) and ((oa.l + ob.l + oc.l + od.l + parityY) % 2 > 0))
                continue;
              if ((std::abs(oa.tz2 + ob.tz2 - oc.tz2 - od.tz2) != 2 * TzX) and (std::abs(oa.tz2 + ob.tz2 - oc.tz2 - od.tz2) != 2 * TzY))
                continue;
              for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
              {
                for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                {
                  for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                  {
                    for (int md = -od.j2; md <= od.j2; md += 2)
                    {
                      if ((ma + mb) != (mc + md))
                        continue;
                      double xabcd = GetMschemeMatrixElement_2b(X, a, ma, b, mb, c, mc, d, md);
                      double xcdab = GetMschemeMatrixElement_2b(X, c, mc, d, md, a, ma, b, mb);
                      double yabicdj = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, i, mi, c, mc, d, md, j, mj);
                      double ycdiabj = GetMschemeMatrixElement_3b(Y, c, mc, d, md, i, mi, a, ma, b, mb, j, mj);

                      double yabcd = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, c, mc, d, md);
                      double ycdab = GetMschemeMatrixElement_2b(Y, c, mc, d, md, a, ma, b, mb);
                      double xabicdj = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, c, mc, d, md, j, mj);
                      double xcdiabj = GetMschemeMatrixElement_3b(X, c, mc, d, md, i, mi, a, ma, b, mb, j, mj);

                      Zm_ij += (1. / 4) * na * nb * (1 - nc) * (1 - nd) * ((xabcd * ycdiabj - yabicdj * xcdab) - (yabcd * xcdiabj - xabicdj * ycdab));

                    } // for md
                  } // for mc
                } // for mb
              } // for ma

            } // for d
          } // for c
        } // for b
      } // for a

      double ZJ_ij = Z_J.OneBody(i, j);
      double err = Zm_ij - ZJ_ij;
      if (std::abs(err) > 1e-6)
      {
        std::cout << "Trouble in " << __func__ << "  i,j = " << i << " " << j
                  << "   zij = " << Zm_ij << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl;
      }
      summed_error += err * err;
      sum_m += Zm_ij * Zm_ij;
      sum_J += ZJ_ij * ZJ_ij;

    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ijkl = sum_ab (na-nb) (X_ab * Y_ijbkla)
//
bool UnitTest::Mscheme_Test_comm132ss(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;

  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();

  Commutator::comm132ss(X, Y, Z_J);
  //  ReferenceImplementations::comm132ss( X, Y, Z_J);

  //  if ( Z_J.IsHermitian() )
  //     Z_J.Symmetrize();
  //  else if (Z_J.IsAntiHermitian() )
  //     Z_J.AntiSymmetrize();

  std::cout << "X one body : " << std::endl
            << X.OneBody << std::endl;
  std::cout << "Y one body : " << std::endl
            << Y.OneBody << std::endl;

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = (parityX + parityY) % 2;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          //          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
          //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;
          if ((oi.l + oj.l + ok.l + ol.l + parityZ) % 2 > 0)
            continue;
          if (std::abs(oi.tz2 + oj.tz2 - ok.tz2 - ol.tz2) != TzZ * 2)
            continue;

          for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
          {
            for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
            {
              if (std::abs(mk) > ol.j2)
                continue; // this just eliminates some redundant matrix elements (I think?)
              int ml = mi + mj - mk;
              if (std::abs(ml) > ol.j2)
                continue;

              //                  if ( mi != mk) continue;
              //             if (not (i==0 and j==0 and k==0 and l==10)) continue;
              double Zm_ijkl = 0;
              for (auto a : X.modelspace->all_orbits)
              {
                Orbit &oa = X.modelspace->GetOrbit(a);
                double na = oa.occ;
                for (auto b : X.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}))
                {
                  Orbit &ob = X.modelspace->GetOrbit(b);
                  double nb = ob.occ;

                  for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                  {
                    int mb = ma;

                    double xab = GetMschemeMatrixElement_1b(X, a, ma, b, mb);
                    double yab = GetMschemeMatrixElement_1b(Y, a, ma, b, mb);
                    double xijbkla = GetMschemeMatrixElement_3b(X, i, mi, j, mj, b, mb, k, mk, l, ml, a, ma);
                    double yijbkla = GetMschemeMatrixElement_3b(Y, i, mi, j, mj, b, mb, k, mk, l, ml, a, ma);

                    Zm_ijkl += (na - nb) * (xab * yijbkla - yab * xijbkla);

                  } // for ma
                } // for b
              } // for a

              double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
              double err = Zm_ijkl - ZJ_ijkl;
              if (std::abs(err) > 1e-6)
              {
                std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                          << " {m} = " << mi << " " << mj << " " << mk << " " << ml
                          << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
              }
              summed_error += err * err;
              sum_m += Zm_ijkl * Zm_ijkl;
              sum_J += ZJ_ijkl * ZJ_ijkl;

            } // for mk
          } // for mj
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ijkl = -1/2 * sum_abc (n_a*n_b*nbar_c + nbar_a*nbar_b*n_c) * [  Xicab*Yabjklc - Xjcab*Yabiklc - Yijcabl*Xabkc + Yijcabk*Xablc
//                                                                 - Yicab*Xabjklc + Yjcab*Xabiklc + Xijcabl*Yabkc - Xijcabk*Yablc ]
//
bool UnitTest::Mscheme_Test_comm232ss(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;

  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();

  Commutator::comm232ss(X, Y, Z_J);
  //  ReferenceImplementations::comm232ss( X, Y, Z_J);
  //  Z_J.Erase();
  //  Commutator::comm232ss_debug( X, Y, Z_J);
  //  Commutator::comm232ss_slow( X, Y, Z_J);

  //  if ( Z_J.IsHermitian() )
  //     Z_J.Symmetrize();
  //  else if (Z_J.IsAntiHermitian() )
  //     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;
  size_t norb = X.modelspace->GetNumberOrbits();
  std::cout << "Begin m-scheme loops" << std::endl;

  for (auto i : X.modelspace->all_orbits)
  {
    std::cout << " orbit " << i << " of " << X.modelspace->all_orbits.size() << std::endl;
    Orbit &oi = X.modelspace->GetOrbit(i);
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k or (k == i and l < j))
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          //          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue; // check parity
          //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue; // check isospin projection
          //         if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue; // check parity
          //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue; // check isospin projection

          int m_max = std::min(oi.j2 + oj.j2, ok.j2 + ol.j2);
          //          for (int mi=-oi.j2; mi<=oi.j2; mi+=2)
          for (int mi = oi.j2; mi <= oi.j2; mi += 2)
          {
            for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
            {
              if (mi + mj < 0)
                continue;
              if (mi + mj != m_max)
                continue;
              for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
              {
                if (std::abs(mk) > ol.j2)
                  continue; // this just eliminates some redundant matrix elements (I think?)
                int ml = mi + mj - mk;
                if (std::abs(ml) > ol.j2)
                  continue;

                double Zm_ijkl = 0;

                //             for (auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ijkl)
                for (size_t a = 0; a < norb; a++)
                {
                  Orbit &oa = X.modelspace->GetOrbit(a);
                  double na = oa.occ;
                  for (auto b : X.modelspace->all_orbits)
                  {
                    Orbit &ob = X.modelspace->GetOrbit(b);
                    double nb = ob.occ;

                    for (auto c : X.modelspace->all_orbits)
                    {
                      Orbit &oc = X.modelspace->GetOrbit(c);
                      double nc = oc.occ;
                      double occfactor = na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc;
                      if (std::abs(occfactor) < 1e-8)
                        continue;

                      for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                      {
                        for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                        {
                          for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                          {

                            double xicab = GetMschemeMatrixElement_2b(X, i, mi, c, mc, a, ma, b, mb);
                            double yicab = GetMschemeMatrixElement_2b(Y, i, mi, c, mc, a, ma, b, mb);
                            double xjcab = GetMschemeMatrixElement_2b(X, j, mj, c, mc, a, ma, b, mb);
                            double yjcab = GetMschemeMatrixElement_2b(Y, j, mj, c, mc, a, ma, b, mb);
                            double xabkc = GetMschemeMatrixElement_2b(X, a, ma, b, mb, k, mk, c, mc);
                            double yabkc = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, k, mk, c, mc);
                            double xablc = GetMschemeMatrixElement_2b(X, a, ma, b, mb, l, ml, c, mc);
                            double yablc = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, l, ml, c, mc);

                            double xabjklc = GetMschemeMatrixElement_3b(X, a, ma, b, mb, j, mj, k, mk, l, ml, c, mc);
                            double yabjklc = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, j, mj, k, mk, l, ml, c, mc);
                            double xabiklc = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, k, mk, l, ml, c, mc);
                            double yabiklc = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, i, mi, k, mk, l, ml, c, mc);
                            double xijcabl = GetMschemeMatrixElement_3b(X, i, mi, j, mj, c, mc, a, ma, b, mb, l, ml);
                            double yijcabl = GetMschemeMatrixElement_3b(Y, i, mi, j, mj, c, mc, a, ma, b, mb, l, ml);
                            double xijcabk = GetMschemeMatrixElement_3b(X, i, mi, j, mj, c, mc, a, ma, b, mb, k, mk);
                            double yijcabk = GetMschemeMatrixElement_3b(Y, i, mi, j, mj, c, mc, a, ma, b, mb, k, mk);

                            //                        double occfactor = na*nb*(1-nc) + (1-na)*(1-nb)*nc;

                            Zm_ijkl += -0.5 * occfactor * (xicab * yabjklc - xjcab * yabiklc - yijcabl * xabkc + yijcabk * xablc - yicab * xabjklc + yjcab * xabiklc + xijcabl * yabkc - xijcabk * yablc);
                            //                        Zm_ijkl += -0.5* occfactor * ( xicab*yabjklc - xjcab*yabiklc - yijcabl*xabkc + yijcabk*xablc
                            //                                                     - yicab*xabjklc + yjcab*xabiklc + xijcabl*yabkc - xijcabk*yablc );

                          } // for mc
                        } // for mb
                      } // for ma
                    } // for c
                  } // for b
                } // for a

                double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
                double err = Zm_ijkl - ZJ_ijkl;
                if (std::abs(err) > 1e-6)
                {
                  std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                            << " {m} = " << mi << " " << mj << " " << mk << " " << ml
                            << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
                }
                summed_error += err * err;
                sum_m += Zm_ijkl * Zm_ijkl;
                sum_J += ZJ_ijkl * ZJ_ijkl;

              } // for mk
            } // for mj
          } // for mi
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ijkl = 1/6 * sum_abcd (n_a*n_b*n_c*nbar_d - nbar_a*nbar_b*nbar_c*n_d) * [  Xijdabc*Yabckld - Yijdabc*Xabckld  ]
// THIS HAS AN OVERALL MINUS SIGN ERROR! (Thanks to Victor Vaida for pointing this out)
// Corrected Sep 3 2025. SRS.
//
bool UnitTest::Mscheme_Test_comm332_ppph_hhhpss(const Operator &X, const Operator &Y) // test not yet implemented
{
  std::cout << __func__ << std::endl;
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();

  ///Commutator::comm332_ppph_hhhpss(X, Y, Z_J);
  ReferenceImplementations::comm332_ppph_hhhpss(X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if ((oi.l + oj.l + ok.l + ol.l) % 2 > 0)
            continue; // check parity
          // if ((oi.tz2 + oj.tz2) != (ok.tz2 + ol.tz2))
          //  continue; // check isospin projection

          for (int mi = oi.j2; mi <= oi.j2; mi += 2) // notice the lack of minus sign in the initial mi value
          {
            for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
            {
              for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
              {
                if (std::abs(mk) > ol.j2)
                  continue; // this just eliminates some redundant matrix elements (I think?)
                int ml = mi + mj - mk;
                if (std::abs(ml) > ol.j2)
                  continue;

                double Zm_ijkl = 0;
                size_t norb = X.modelspace->GetNumberOrbits();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ijkl)
                for (size_t a = 0; a < norb; a++)
                {
                  Orbit &oa = X.modelspace->GetOrbit(a);
                  double na = oa.occ;
                  for (auto b : X.modelspace->all_orbits)
                  {
                    Orbit &ob = X.modelspace->GetOrbit(b);
                    double nb = ob.occ;

                    for (auto c : X.modelspace->all_orbits)
                    {
                      Orbit &oc = X.modelspace->GetOrbit(c);
                      double nc = oc.occ;

                      for (auto d : X.modelspace->all_orbits)
                      {
                        Orbit &od = X.modelspace->GetOrbit(d);
                        double nd = od.occ;
//                        double occfactor = na * nb * nc * (1 - nd) - (1 - na) * (1 - nb) * (1 - nc) * nd;
                        double occfactor = (1 - na) * (1 - nb) * (1 - nc) * nd -  na * nb * nc * (1 - nd); // Corrected Sep 3 2025 (SRS) 

                        // These make the contribution trivially zero, so I skip them in the name of efficiently
                        // testing the more complicated part. Commenting them out allows to check that the trivial stuff is right.
                        if ((oi.l + oj.l + od.l + oa.l + ob.l + oc.l) % 2 > 0)
                          continue;
                        if ((oi.tz2 + oj.tz2 + od.tz2) != (oa.tz2 + ob.tz2 + oc.tz2))
                          continue;
                        if (std::abs(occfactor) < 1e-7)
                          continue;

                        for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                        {
                          for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                          {
                            for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                            {
                              for (int md = -od.j2; md <= od.j2; md += 2)
                              {
                                if ((mi + mj + md) != (ma + mb + mc))
                                  continue;

                                double xijdabc = GetMschemeMatrixElement_3b(X, i, mi, j, mj, d, md, a, ma, b, mb, c, mc);
                                double yijdabc = GetMschemeMatrixElement_3b(Y, i, mi, j, mj, d, md, a, ma, b, mb, c, mc);
                                double xabckld = GetMschemeMatrixElement_3b(X, a, ma, b, mb, c, mc, k, mk, l, ml, d, md);
                                double yabckld = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, c, mc, k, mk, l, ml, d, md);

                                Zm_ijkl += 1. / 6 * occfactor * (xijdabc * yabckld - yijdabc * xabckld);

                              } // for md
                            } // for mc
                          } // for mb
                        } // for ma
                      } // for d
                    } // for c
                  } // for b
                } // for a

                double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
                double err = Zm_ijkl - ZJ_ijkl;
                if (std::abs(err) > 1e-6)
                {
                  std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                            << " {m} = " << mi << " " << mj << " " << mk << " " << ml
                            << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
                }
                summed_error += err * err;
                sum_m += Zm_ijkl * Zm_ijkl;
                sum_J += ZJ_ijkl * ZJ_ijkl;

              } // for mk
            } // for mj
          } // for mi
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
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
bool UnitTest::Mscheme_Test_comm332_pphhss(const Operator &X, const Operator &Y) // test not yet implemented
{
  std::cout << __func__ << std::endl;
  Operator Z_J(Y);
  Operator Z_J_old(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J_old.SetHermitian();
  Z_J_old.Erase();

  //  Z_J.modelspace->SetOccNat3Cut(1e-5);
  //  Z_J.modelspace->SetdE3max(1);

  //  Commutator::comm332_pphhss_debug( X, Y, Z_J);
  Commutator::comm332_pphhss_debug(X, Y, Z_J_old);
  //  Z_J.Erase();
  Commutator::comm332_pphhss(X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if ((oi.l + oj.l + ok.l + ol.l) % 2 > 0)
            continue; // check parity
          if ((oi.tz2 + oj.tz2) != (ok.tz2 + ol.tz2))
            continue; // check isospin projection

          //          for (int mi=-oi.j2; mi<=oi.j2; mi+=2)
          for (int mi = oi.j2; mi <= oi.j2; mi += 2)
          {
            for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
            {
              for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
              {
                if (std::abs(mk) > ol.j2)
                  continue; // this just eliminates some redundant matrix elements (I think?)
                int ml = mi + mj - mk;
                if (std::abs(ml) > ol.j2)
                  continue;

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
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ijkl)
                for (size_t a = 0; a < norb; a++)
                {
                  Orbit &oa = X.modelspace->GetOrbit(a);
                  double na = oa.occ;
                  for (auto b : X.modelspace->all_orbits)
                  {
                    Orbit &ob = X.modelspace->GetOrbit(b);
                    double nb = ob.occ;

                    for (auto c : X.modelspace->all_orbits)
                    {
                      Orbit &oc = X.modelspace->GetOrbit(c);
                      double nc = oc.occ;

                      for (auto d : X.modelspace->all_orbits)
                      {
                        Orbit &od = X.modelspace->GetOrbit(d);
                        double nd = od.occ;

                        if ((oa.l + ob.l + oi.l + oc.l + od.l + ok.l) % 2 > 0 and (oa.l + ob.l + oi.l + oc.l + od.l + ol.l) % 2 > 0 and (oa.l + ob.l + oj.l + oc.l + od.l + ol.l) % 2 > 0 and (oa.l + ob.l + oj.l + oc.l + od.l + ok.l) % 2 > 0)
                          continue; // check parity for X and Y

                        if ((oa.tz2 + ob.tz2 + oi.tz2) != (oc.tz2 + od.tz2 + ok.tz2) and (oa.tz2 + ob.tz2 + oi.tz2) != (oc.tz2 + od.tz2 + ol.tz2) and (oa.tz2 + ob.tz2 + oj.tz2) != (oc.tz2 + od.tz2 + ol.tz2) and (oa.tz2 + ob.tz2 + oj.tz2) != (oc.tz2 + od.tz2 + ok.tz2))
                          continue; // check isospin projection

                        for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                        {
                          for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                          {
                            for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                            {
                              for (int md = -od.j2; md <= od.j2; md += 2)
                              {
                                // Z_ijkl = 1/4  sum_abcd (n_a*n_b*nbar_c*nbar_d) * (  Xabicdk*Ycdjabl - Yabicdk*Xcdjabl
                                //                                                   - Xabjcdk*Ycdiabl + Yabjcdk*Xcdiabl
                                //                                                   - Xabicdl*Ycdjabk + Yabicdl*Xcdjabk
                                //                                                   + Xabjcdl*Ycdiabk - Yabjcdl*Xcdiabk )

                                double xabicdk = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, c, mc, d, md, k, mk);
                                double xabicdl = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, c, mc, d, md, l, ml);
                                double xabjcdk = GetMschemeMatrixElement_3b(X, a, ma, b, mb, j, mj, c, mc, d, md, k, mk);
                                double xabjcdl = GetMschemeMatrixElement_3b(X, a, ma, b, mb, j, mj, c, mc, d, md, l, ml);

                                double ycdiabk = GetMschemeMatrixElement_3b(Y, c, mc, d, md, i, mi, a, ma, b, mb, k, mk);
                                double ycdiabl = GetMschemeMatrixElement_3b(Y, c, mc, d, md, i, mi, a, ma, b, mb, l, ml);
                                double ycdjabk = GetMschemeMatrixElement_3b(Y, c, mc, d, md, j, mj, a, ma, b, mb, k, mk);
                                double ycdjabl = GetMschemeMatrixElement_3b(Y, c, mc, d, md, j, mj, a, ma, b, mb, l, ml);

                                double occfactor = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);

                                //                          Zm_ijkl += 1./4 * occfactor * (xabicdl*ycdjabk
                                double dz = 1. / 4 * occfactor * (xabicdl * ycdjabk - xabjcdl * ycdiabk - xabicdk * ycdjabl + xabjcdk * ycdiabl);
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

                              } // for md
                            } // for mc
                          } // for mb
                        } // for ma
                      } // for d
                    } // for c
                  } // for b
                } // for a

                double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
                double ZJ_old_ijkl = GetMschemeMatrixElement_2b(Z_J_old, i, mi, j, mj, k, mk, l, ml);
                double err = Zm_ijkl - ZJ_ijkl;
                //             if (std::abs(ZJ_ijkl-ZJ_old_ijkl)>1e-6)
                if (std::abs(err) > 1e-6)
                {
                  std::cout << "Trouble in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                            << " {m} = " << mi << " " << mj << " " << mk << " " << ml
                            //                         << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
                            << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << " ZJ_old_ijkl " << ZJ_old_ijkl << "   err = " << err
                            << "  J err = " << ZJ_ijkl - ZJ_old_ijkl << std::endl;
                }
                summed_error += err * err;
                sum_m += Zm_ijkl * Zm_ijkl;
                sum_J += ZJ_ijkl * ZJ_ijkl;

              } // for mk
            } // for mj
          } // for mi
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

//
// M-scheme formula
//
// Z_ijklmn = (1-Pik-Pjk)(1-Plm-Pln) sum_a (X_ijla Y_akmn - Y_ijla X_akmn)
//
bool UnitTest::Mscheme_Test_comm223ss(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = (parityX + parityY) % 2;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;
  int JrankZ = X.GetJRank() + Y.GetJRank();

  Operator Z_J(*(X.modelspace), JrankZ, TzZ, parityZ, 3);
  //  Operator Z_J( Y );
  Z_J.SetHermitian();
  //  Z_J.Erase();
  Z_J.ThreeBody.SetMode("pn");
  Z_J.ThreeBody.Allocate();

  std::cout << "Calling J-scheme commutator" << std::endl;
  Commutator::comm223ss(X, Y, Z_J);
  //  ReferenceImplementations::comm223ss( X, Y, Z_J);

  //  Commutator::comm223ss_new( X, Y, Z_J);
  //  if (X.modelspace->GetEmax() > 1) return true;

  //  if ( Z_J.IsHermitian() )
  //     Z_J.Symmetrize();
  //  else if (Z_J.IsAntiHermitian() )
  //     Z_J.AntiSymmetrize();
  //  std::cout << "  done." << std::endl;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;
  std::cout << "start mscheme loops " << std::endl;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j > i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k > j)
          continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l > i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            if (m > l)
              continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n > m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              //              if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2 !=0 ) continue;
              //              if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l + parityZ) % 2 != 0)
                continue;
              if (std::abs(oi.tz2 + oj.tz2 + ok.tz2 - ol.tz2 - om.tz2 - on.tz2) != 2 * TzZ)
                continue;

              // loop over projections
              for (int m_i = -oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ((m_i + m_j + m_k) != (m_l + m_m + m_n))
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;
                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                            {
                              // "direct" term  Z1
                              double x_ijla = GetMschemeMatrixElement_2b(X, i, m_i, j, m_j, l, m_l, a, m_a);
                              double y_akmn = GetMschemeMatrixElement_2b(Y, a, m_a, k, m_k, m, m_m, n, m_n);
                              double y_ijla = GetMschemeMatrixElement_2b(Y, i, m_i, j, m_j, l, m_l, a, m_a);
                              double x_akmn = GetMschemeMatrixElement_2b(X, a, m_a, k, m_k, m, m_m, n, m_n);
                              // Pik  Z2
                              double x_kjla = GetMschemeMatrixElement_2b(X, k, m_k, j, m_j, l, m_l, a, m_a);
                              double y_aimn = GetMschemeMatrixElement_2b(Y, a, m_a, i, m_i, m, m_m, n, m_n);
                              double y_kjla = GetMschemeMatrixElement_2b(Y, k, m_k, j, m_j, l, m_l, a, m_a);
                              double x_aimn = GetMschemeMatrixElement_2b(X, a, m_a, i, m_i, m, m_m, n, m_n);
                              // Pjk   Z3
                              double x_ikla = GetMschemeMatrixElement_2b(X, i, m_i, k, m_k, l, m_l, a, m_a);
                              double y_ajmn = GetMschemeMatrixElement_2b(Y, a, m_a, j, m_j, m, m_m, n, m_n);
                              double y_ikla = GetMschemeMatrixElement_2b(Y, i, m_i, k, m_k, l, m_l, a, m_a);
                              double x_ajmn = GetMschemeMatrixElement_2b(X, a, m_a, j, m_j, m, m_m, n, m_n);
                              // Plm   Z4
                              double x_ijma = GetMschemeMatrixElement_2b(X, i, m_i, j, m_j, m, m_m, a, m_a);
                              double y_akln = GetMschemeMatrixElement_2b(Y, a, m_a, k, m_k, l, m_l, n, m_n);
                              double y_ijma = GetMschemeMatrixElement_2b(Y, i, m_i, j, m_j, m, m_m, a, m_a);
                              double x_akln = GetMschemeMatrixElement_2b(X, a, m_a, k, m_k, l, m_l, n, m_n);
                              // Pln   Z5
                              double x_ijna = GetMschemeMatrixElement_2b(X, i, m_i, j, m_j, n, m_n, a, m_a);
                              double y_akml = GetMschemeMatrixElement_2b(Y, a, m_a, k, m_k, m, m_m, l, m_l);
                              double y_ijna = GetMschemeMatrixElement_2b(Y, i, m_i, j, m_j, n, m_n, a, m_a);
                              double x_akml = GetMschemeMatrixElement_2b(X, a, m_a, k, m_k, m, m_m, l, m_l);
                              // Pik Plm    Z6
                              double x_kjma = GetMschemeMatrixElement_2b(X, k, m_k, j, m_j, m, m_m, a, m_a);
                              double y_ailn = GetMschemeMatrixElement_2b(Y, a, m_a, i, m_i, l, m_l, n, m_n);
                              double y_kjma = GetMschemeMatrixElement_2b(Y, k, m_k, j, m_j, m, m_m, a, m_a);
                              double x_ailn = GetMschemeMatrixElement_2b(X, a, m_a, i, m_i, l, m_l, n, m_n);
                              // Pik Pln    Z7
                              double x_kjna = GetMschemeMatrixElement_2b(X, k, m_k, j, m_j, n, m_n, a, m_a);
                              double y_aiml = GetMschemeMatrixElement_2b(Y, a, m_a, i, m_i, m, m_m, l, m_l);
                              double y_kjna = GetMschemeMatrixElement_2b(Y, k, m_k, j, m_j, n, m_n, a, m_a);
                              double x_aiml = GetMschemeMatrixElement_2b(X, a, m_a, i, m_i, m, m_m, l, m_l);
                              // Pjk Plm    Z8
                              double x_ikma = GetMschemeMatrixElement_2b(X, i, m_i, k, m_k, m, m_m, a, m_a);
                              double y_ajln = GetMschemeMatrixElement_2b(Y, a, m_a, j, m_j, l, m_l, n, m_n);
                              double y_ikma = GetMschemeMatrixElement_2b(Y, i, m_i, k, m_k, m, m_m, a, m_a);
                              double x_ajln = GetMschemeMatrixElement_2b(X, a, m_a, j, m_j, l, m_l, n, m_n);
                              // Pjk Pln    Z9
                              double x_ikna = GetMschemeMatrixElement_2b(X, i, m_i, k, m_k, n, m_n, a, m_a);
                              double y_ajml = GetMschemeMatrixElement_2b(Y, a, m_a, j, m_j, m, m_m, l, m_l);
                              double y_ikna = GetMschemeMatrixElement_2b(Y, i, m_i, k, m_k, n, m_n, a, m_a);
                              double x_ajml = GetMschemeMatrixElement_2b(X, a, m_a, j, m_j, m, m_m, l, m_l);
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

                            } // for m_a
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "Trouble in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  //  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

bool UnitTest::Mscheme_Test_comm133ss(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;

  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  Commutator::comm133ss(X, Y, Z_J);
  //  ReferenceImplementations::comm133ss( X, Y, Z_J);

  /*
    if ( Z_J.IsHermitian() )
       Z_J.Symmetrize();
    else if (Z_J.IsAntiHermitian() )
       Z_J.AntiSymmetrize();
  */

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = (parityX + parityY) % 2;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < j)
          continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l < i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            if (m < l)
              continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n < m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              //              if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2 !=0 ) continue;
              //              if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l + parityZ) % 2 != 0)
                continue;
              if (std::abs(oi.tz2 + oj.tz2 + ok.tz2 - ol.tz2 - om.tz2 - on.tz2) != TzZ)
                continue;

              //              if ( not (i==1 and j==0 and k==0 and l==1 and m==0 and n==0) ) continue;

              // loop over projections
              for (int m_i = oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if (m_i + m_j + m_k < 0)
                      continue;
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ((m_i + m_j + m_k) != (m_l + m_m + m_n))
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          //                     if (not (m_i==1 and m_j==1 and m_k==-1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            //                       if (a!=11) continue;
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            double xia = X.OneBody(i, a);
                            double xja = X.OneBody(j, a);
                            double xka = X.OneBody(k, a);
                            double xal = X.OneBody(a, l);
                            double xam = X.OneBody(a, m);
                            double xan = X.OneBody(a, n);

                            double yia = Y.OneBody(i, a);
                            double yja = Y.OneBody(j, a);
                            double yka = Y.OneBody(k, a);
                            double yal = Y.OneBody(a, l);
                            double yam = Y.OneBody(a, m);
                            double yan = Y.OneBody(a, n);
                            double xajklmn = 0;
                            double xiaklmn = 0;
                            double xijalmn = 0;
                            double xijkamn = 0;
                            double xijklan = 0;
                            double xijklma = 0;
                            double yajklmn = 0;
                            double yiaklmn = 0;
                            double yijalmn = 0;
                            double yijkamn = 0;
                            double yijklan = 0;
                            double yijklma = 0;
                            for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                            {

                              //                         double xijalmn = GetMschemeMatrixElement_3b( X, i,m_i, j,m_j, a, m_a, l, m_l, m, m_m, n,m_n );

                              //                         if ( m_a==m_k and oa.j2==ok.j2 and oa.l==ok.l and oa.tz2==ok.tz2)
                              if (m_a == m_k and (std::abs(xka) > 1e-8 or std::abs(yka) > 1e-8))
                              {
                                xijalmn += GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, a, m_a, l, m_l, m, m_m, n, m_n);
                                yijalmn += GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, a, m_a, l, m_l, m, m_m, n, m_n);
                              }

                              //                         if ( m_a==m_j and oa.j2==oj.j2 and oa.l==oj.l and oa.tz2==oj.tz2)
                              if (m_a == m_j and (std::abs(xja) > 1e-8 or std::abs(yja) > 1e-8))
                              {
                                xiaklmn += GetMschemeMatrixElement_3b(X, i, m_i, a, m_a, k, m_k, l, m_l, m, m_m, n, m_n);
                                yiaklmn += GetMschemeMatrixElement_3b(Y, i, m_i, a, m_a, k, m_k, l, m_l, m, m_m, n, m_n);
                              }

                              //                         if ( m_a==m_i and oa.j2==oi.j2 and oa.l==oi.l and oa.tz2==oi.tz2)
                              if (m_a == m_i and (std::abs(xia) > 1e-8 or std::abs(yia) > 1e-8))
                              {
                                xajklmn += GetMschemeMatrixElement_3b(X, a, m_a, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                                yajklmn += GetMschemeMatrixElement_3b(Y, a, m_a, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                              }

                              //                         if ( m_a==m_l and oa.j2==ol.j2 and oa.l==ol.l and oa.tz2==ol.tz2)
                              if (m_a == m_l and (std::abs(xal) > 1e-8 or std::abs(yal) > 1e-8))
                              {
                                xijkamn += GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, m, m_m, n, m_n);
                                yijkamn += GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, m, m_m, n, m_n);
                              }

                              //                         if ( m_a==m_m and oa.j2==om.j2 and oa.l==om.l and oa.tz2==om.tz2)
                              if (m_a == m_m and (std::abs(xam) > 1e-8 or std::abs(yam) > 1e-8))
                              {
                                xijklan += GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, l, m_l, a, m_a, n, m_n);
                                yijklan += GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, l, m_l, a, m_a, n, m_n);
                              }

                              //                         if ( m_a==m_n and oa.j2==on.j2 and oa.l==on.l and oa.tz2==on.tz2)
                              if (m_a == m_n and (std::abs(xan) > 1e-8 or std::abs(yan) > 1e-8))
                              {
                                xijklma += GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, a, m_a);
                                yijklma += GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, a, m_a);
                              }

                            } // for m_a
                            z_ijklmn += xia * yajklmn + xja * yiaklmn + xka * yijalmn - yijkamn * xal - yijklan * xam - yijklma * xan;
                            z_ijklmn -= yia * xajklmn + yja * xiaklmn + yka * xijalmn - xijkamn * yal - xijklan * yam - xijklma * yan;
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "Trouble in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  //  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
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
bool UnitTest::Mscheme_Test_comm233_pp_hhss(const Operator &X, const Operator &Y) // test not yet implemented
{
  std::cout << __func__ << std::endl;

  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  Commutator::comm233_pp_hhss(X, Y, Z_J);
  //  Z_J.Erase();
  //  Commutator::comm233_pp_hhss_debug( X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;

  size_t norb = X.modelspace->GetNumberOrbits();
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < j)
          continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l < i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            if (m < l)
              continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n < m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l) % 2 != 0)
                continue;
              if ((oi.tz2 + oj.tz2 + ok.tz2) != (ol.tz2 + om.tz2 + on.tz2))
                continue;

              //              if ( not (i==0 and j==0 and k==1 and l==0 and m==2 and n==3) ) continue;
              //              if ( not (i==1 and j==0 and k==0 and l==1 and m==0 and n==0) ) continue;
              //              if ( not (i==0 and j==0 and k==2 and l==2 and m==2 and n==4) ) continue;
              //              if ( not (i==4 and j==2 and k==2 and l==2 and m==0 and n==0) ) continue;

              //              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i = oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ((m_i + m_j + m_k) != (m_l + m_m + m_n))
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          //                     if (not (m_i==1 and m_j==-1 and m_k==1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                          double z_ijklmn = 0;

                          //                     std::cout << " {m} " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n << std::endl;
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (size_t b = 0; b < norb; b++)
                            {
                              //                        if (not (( a==4 and b==5) or (a==5 and b==4))  ) continue;
                              //                        if (not (( a==2 and b==2) or (a==2 and b==2))  ) continue;
                              Orbit &ob = X.modelspace->GetOrbit(b);
                              double occfactor = 1 - oa.occ - ob.occ;
                              //                       if (a!=11) continue;
                              for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                              {
                                for (int m_b = -ob.j2; m_b <= ob.j2; m_b += 2)
                                {

                                  ///       = 1/2 sum_ab (1-na-nb) [ Xijab*Yabklmn - Yijab*Xabklmn - Xkjab*Yabilmn + Ykjab*Xabilmn - Xikab*Yabjlmn + Yikab*Xabjlmn
                                  ///                              - Yijkabn*Xablm + Xijkabn*Yablm + Yijkabl*Xabnm - Xijkabl*Yabnm + Yijkabm*Xabln - Xijkabm*Yabln ]
                                  double xabilmn = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, i, m_i, l, m_l, m, m_m, n, m_n);
                                  double xabjlmn = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, j, m_j, l, m_l, m, m_m, n, m_n);
                                  double xabklmn = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, k, m_k, l, m_l, m, m_m, n, m_n);
                                  double xijkabl = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, l, m_l);
                                  double xijkabm = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, m, m_m);
                                  double xijkabn = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, n, m_n);

                                  double yabilmn = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, i, m_i, l, m_l, m, m_m, n, m_n);
                                  double yabjlmn = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, j, m_j, l, m_l, m, m_m, n, m_n);
                                  double yabklmn = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, k, m_k, l, m_l, m, m_m, n, m_n);
                                  double yijkabl = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, l, m_l);
                                  double yijkabm = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, m, m_m);
                                  double yijkabn = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, n, m_n);

                                  double xijab = GetMschemeMatrixElement_2b(X, i, m_i, j, m_j, a, m_a, b, m_b);
                                  double xkjab = GetMschemeMatrixElement_2b(X, k, m_k, j, m_j, a, m_a, b, m_b);
                                  double xikab = GetMschemeMatrixElement_2b(X, i, m_i, k, m_k, a, m_a, b, m_b);
                                  double xablm = GetMschemeMatrixElement_2b(X, a, m_a, b, m_b, l, m_l, m, m_m);
                                  double xabnm = GetMschemeMatrixElement_2b(X, a, m_a, b, m_b, n, m_n, m, m_m);
                                  double xabln = GetMschemeMatrixElement_2b(X, a, m_a, b, m_b, l, m_l, n, m_n);

                                  double yijab = GetMschemeMatrixElement_2b(Y, i, m_i, j, m_j, a, m_a, b, m_b);
                                  double ykjab = GetMschemeMatrixElement_2b(Y, k, m_k, j, m_j, a, m_a, b, m_b);
                                  double yikab = GetMschemeMatrixElement_2b(Y, i, m_i, k, m_k, a, m_a, b, m_b);
                                  double yablm = GetMschemeMatrixElement_2b(Y, a, m_a, b, m_b, l, m_l, m, m_m);
                                  double yabnm = GetMschemeMatrixElement_2b(Y, a, m_a, b, m_b, n, m_n, m, m_m);
                                  double yabln = GetMschemeMatrixElement_2b(Y, a, m_a, b, m_b, l, m_l, n, m_n);

                                  //                         z_ijklmn += 0.5*occfactor*( (xijab*yabklmn - yijab*xabklmn) - (xkjab*yabilmn - ykjab*xabilmn) - (xikab*yabjlmn - yikab*xabjlmn)
                                  //                                                   - (yijkabn*xablm - xijkabn*yablm) + (yijkabl*xabnm - xijkabl*yabnm) + (yijkabm*xabln + xijkabm*yabln) );
                                  double dz = 0.5 * occfactor * ((1 * xijab * yabklmn - 1 * yijab * xabklmn) - (1 * xkjab * yabilmn - 1 * ykjab * xabilmn) - (1 * xikab * yabjlmn - 1 * yikab * xabjlmn) - (1 * yijkabn * xablm - 1 * xijkabn * yablm) + (1 * yijkabl * xabnm - 1 * xijkabl * yabnm) + (1 * yijkabm * xabln - 1 * xijkabm * yabln));
                                  z_ijklmn += dz;
                                  //                         if (std::abs(dz)>1e-6)
                                  //                         {
                                  //                           std::cout << "a,b = " << a << " " << b << "  ma mb = " << m_a << " " << m_b << "  dz = " << dz << "  => " << z_ijklmn << std::endl;
                                  //                           std::cout << "    " << xijab * yabklmn << "  " << yijab*xabklmn << " " << xkjab*yabilmn << " " << ykjab*xabilmn << " " << xikab*yabjlmn << " " << yikab*xabjlmn << std::endl;
                                  //                           std::cout << "    " << yijkabn*xablm << " " << xijkabn*yablm << " " << yijkabl*xabnm << " " << xijkabl*yabnm << " " << yijkabm*xabln << " " << xijkabm*yabln << std::endl;
                                  //                           std::cout << "   [[  " <<  xijab*yabklmn << " " << xkjab*yabilmn << " " <<  xikab*yabjlmn << " ]]  -->  "
                                  //                                     <<   xkjab << " " << yabilmn  << "   ... " << k << " " << j << " " << m_k << " " << m_j << "   " << a << " " << b << " " << j << "  " << m_a << " " << m_b << " " << m_j     << std::endl;
                                  //                         }
                                } // for m_b
                              } // for m_a
                            } // for b
                          } // for a

                          //                    double zJ_111 = Z_J.ThreeBody.GetME_pn(1,1,1,1,0,0,1,0,0);
                          //                    std::cout << " zJ_111 = " << zJ_111 << std::endl;
                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "Trouble in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  //  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
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
bool UnitTest::Mscheme_Test_comm233_phss(const Operator &X, const Operator &Y) // test not yet implemented
{
  std::cout << __func__ << std::endl;
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  Commutator::comm233_phss(X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j > i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k > j)
          continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l > i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            if (m > l)
              continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n > m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l) % 2 != 0)
                continue;
              if ((oi.tz2 + oj.tz2 + ok.tz2) != (ol.tz2 + om.tz2 + on.tz2))
                continue;

              //              if ( not (i==5 and j==4 and k==0 and l==5 and m==4 and n==0) ) continue;
              //              if ( not (i==3 and j==2 and k==0 and l==3 and m==2 and n==0) ) continue;
              //              if ( not (i==4 and j==1 and k==0 and l==4 and m==3 and n==2) ) continue;

              //              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i = -oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ((m_i + m_j + m_k) != (m_l + m_m + m_n))
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          if (not(m_i == 1 and m_j == -1 and m_k == 1 and m_l == 1 and m_m == -1 and m_n == 1))
                            continue;
                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (size_t b = 0; b < norb; b++)
                            {
                              Orbit &ob = X.modelspace->GetOrbit(b);
                              double occfactor = oa.occ - ob.occ;
                              //                       if (a!=11) continue;
                              //                       if (a!=3 and a!=5) continue;
                              for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                              {
                                for (int m_b = -ob.j2; m_b <= ob.j2; m_b += 2)
                                {

                                  //                         double xajkbmn = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, n,m_n );
                                  //                         double xaikbmn = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, n,m_n );
                                  //                         double xajibmn = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, n,m_n );
                                  //                         double xajkbln = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, l,m_l, n,m_n );
                                  //                         double xaikbln = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, l,m_l, n,m_n );
                                  //                         double xajibln = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, l,m_l, n,m_n );
                                  //                         double xajkbml = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, l,m_l );
                                  //                         double xaikbml = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, l,m_l );
                                  //                         double xajibml = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, l,m_l );

                                  double xbial = GetMschemeMatrixElement_2b(X, b, m_b, i, m_i, a, m_a, l, m_l);
                                  double xbjal = GetMschemeMatrixElement_2b(X, b, m_b, j, m_j, a, m_a, l, m_l);
                                  double xbkal = GetMschemeMatrixElement_2b(X, b, m_b, k, m_k, a, m_a, l, m_l);
                                  double xbiam = GetMschemeMatrixElement_2b(X, b, m_b, i, m_i, a, m_a, m, m_m);
                                  double xbjam = GetMschemeMatrixElement_2b(X, b, m_b, j, m_j, a, m_a, m, m_m);
                                  double xbkam = GetMschemeMatrixElement_2b(X, b, m_b, k, m_k, a, m_a, m, m_m);
                                  double xbian = GetMschemeMatrixElement_2b(X, b, m_b, i, m_i, a, m_a, n, m_n);
                                  double xbjan = GetMschemeMatrixElement_2b(X, b, m_b, j, m_j, a, m_a, n, m_n);
                                  double xbkan = GetMschemeMatrixElement_2b(X, b, m_b, k, m_k, a, m_a, n, m_n);

                                  //                         double yajkbmn = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, n,m_n );
                                  //                         double yaikbmn = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, n,m_n );
                                  //                         double yajibmn = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, n,m_n );
                                  //                         double yajkbln = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, l,m_l, n,m_n );
                                  //                         double yaikbln = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, l,m_l, n,m_n );
                                  //                         double yajibln = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, l,m_l, n,m_n );
                                  //                         double yajkbml = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, l,m_l );
                                  //                         double yaikbml = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, l,m_l );
                                  //                         double yajibml = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, l,m_l );

                                  double ybial = GetMschemeMatrixElement_2b(Y, b, m_b, i, m_i, a, m_a, l, m_l);
                                  double ybjal = GetMschemeMatrixElement_2b(Y, b, m_b, j, m_j, a, m_a, l, m_l);
                                  double ybkal = GetMschemeMatrixElement_2b(Y, b, m_b, k, m_k, a, m_a, l, m_l);
                                  double ybiam = GetMschemeMatrixElement_2b(Y, b, m_b, i, m_i, a, m_a, m, m_m);
                                  double ybjam = GetMschemeMatrixElement_2b(Y, b, m_b, j, m_j, a, m_a, m, m_m);
                                  double ybkam = GetMschemeMatrixElement_2b(Y, b, m_b, k, m_k, a, m_a, m, m_m);
                                  double ybian = GetMschemeMatrixElement_2b(Y, b, m_b, i, m_i, a, m_a, n, m_n);
                                  double ybjan = GetMschemeMatrixElement_2b(Y, b, m_b, j, m_j, a, m_a, n, m_n);
                                  double ybkan = GetMschemeMatrixElement_2b(Y, b, m_b, k, m_k, a, m_a, n, m_n);
                                  //////////////////----------------------------------------------------

                                  double xijalmb = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double yijalmb = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double xkjalmb = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double ykjalmb = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double xijanmb = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double yijanmb = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double xikalmb = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double yikalmb = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double xijalnb = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double yijalnb = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double xkjanmb = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double ykjanmb = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double xikanmb = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double yikanmb = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double xkjalnb = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double ykjalnb = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double xikalnb = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double yikalnb = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, a, m_a, l, m_l, n, m_n, b, m_b);

                                  //  Zijkl =  - sum_ab (na-nb) [ Xbial Yajkbmn - Ybial Xajkbmn  -  Xbjal Yaikbmn + Ybjal Xaikbmn  -  Xbkal Yajibmn + Ybkal Xajibmn
                                  //                            - Xbiam Yajkbln + Ybiam Xajkbln  +  Xbjam Yaikbln - Ybjam Xaikbln  +  Xbkam Yajibln - Ybkam Xajibln
                                  //                            - Xbian Yajkbml + Ybian Xajkbml  +  Xbjan Yaikbml - Ybjan Xaikbml  +  Xbkan Yajibml - Ybkan Xajibml  ]
                                  // xijalmn = xaijblm = xajibml
                                  //                          z_ijklmn -= occfactor * (xbkan * yajibml - ybkan * xajibml);
                                  //                           z_ijklmn -= occfactor * (xbkan * yijalmb - ybkan*xijalmb);
                                  double dz = -occfactor * 1 * (xbkan * yijalmb - ybkan * xijalmb);
                                  dz += occfactor * 1 * (xbian * ykjalmb - ybian * xkjalmb);
                                  dz += occfactor * 1 * (xbjan * yikalmb - ybjan * xikalmb);
                                  dz += occfactor * 1 * (xbkal * yijanmb - ybkal * xijanmb);
                                  dz += occfactor * 1 * (xbkam * yijalnb - ybkam * xijalnb);
                                  dz -= occfactor * 1 * (xbial * ykjanmb - ybial * xkjanmb);
                                  dz -= occfactor * 1 * (xbjal * yikanmb - ybjal * xikanmb);
                                  dz -= occfactor * 1 * (xbiam * ykjalnb - ybiam * xkjalnb);
                                  dz -= occfactor * 1 * (xbjam * yikalnb - ybjam * xikalnb);
                                  z_ijklmn += dz;

                                } // for m_b
                              } // for m_a
                            } // for b
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "Trouble in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  //  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

// M-Scheme Formula:
//
// Z_ijkl = 1/6 sum_abc [na*nb*nc + (1-na)*(1-nb)(1-nc) ] * [  Xijkabc*Yabclmn - Yijkabc*Xabclmn  ]
//
//
bool UnitTest::Mscheme_Test_comm333_ppp_hhhss(const Operator &X, const Operator &Y) // test not yet implemented
{
  std::cout << __func__ << std::endl;
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  Commutator::comm333_ppp_hhhss(X, Y, Z_J);

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      //      if (j<i) continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        //        if (k<j) continue;

        for (auto l : X.modelspace->all_orbits)
        {
          //          if (l<i) continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            //            if (m<l) continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n < m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l) % 2 != 0)
                continue;
              if ((oi.tz2 + oj.tz2 + ok.tz2) != (ol.tz2 + om.tz2 + on.tz2))
                continue;

              //              if ( not (i==1 and j==0 and k==0 and l==1 and m==0 and n==0) ) continue;

              //              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i = oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ((m_i + m_j + m_k) != (m_l + m_m + m_n))
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          //                     if (not (m_i==1 and m_j==1 and m_k==-1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (size_t b = 0; b < norb; b++)
                            {
                              Orbit &ob = X.modelspace->GetOrbit(b);
                              for (size_t c = 0; c < norb; c++)
                              {
                                Orbit &oc = X.modelspace->GetOrbit(c);
                                //                       double occfactor = oa.occ*ob.occ*oc.occ - (1-oa.occ)*(1-ob.occ)*(1-oc.occ); // why a minus sign??? That's wrong.
                                double occfactor = oa.occ * ob.occ * oc.occ + (1 - oa.occ) * (1 - ob.occ) * (1 - oc.occ); // Fixed by SRS July 19 2022
                                                                                                                          //                       if (a!=11) continue;
                                for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                                {
                                  for (int m_b = -ob.j2; m_b <= ob.j2; m_b += 2)
                                  {
                                    for (int m_c = -oc.j2; m_c <= oc.j2; m_c += 2)
                                    {

                                      // Z_ijkl = 1/6 sum_abc [na*nb*nc + (1-na)*(1-nb)(1-nc) ] * [  Xijkabc*Yabclmn - Yijkabc*Xabclmn  ]
                                      double xijkabc = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, c, m_c);
                                      double xabclmn = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, c, m_c, l, m_l, m, m_m, n, m_n);
                                      double yijkabc = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, c, m_c);
                                      double yabclmn = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, c, m_c, l, m_l, m, m_m, n, m_n);

                                      z_ijklmn += 1. / 6 * occfactor * (xijkabc * yabclmn - yijkabc * xabclmn);

                                    } // for m_c
                                  } // for m_b
                                } // for m_a
                              } // for c
                            } // for b
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "Trouble in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  //  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
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
bool UnitTest::Mscheme_Test_comm333_pph_hhpss(const Operator &X, const Operator &Y) // test not yet implemented
{
  std::cout << __func__ << std::endl;
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.Allocate();

  //  Commutator::comm333ss( X, Y, Z_J); // TODO: make a sub-call to just the pph hhp part
  Commutator::comm333_pph_hhpss(X, Y, Z_J); // TODO: make a sub-call to just the pph hhp part

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      //      if (j>i) continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < j)
          continue;
        //        if (k>j) continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l > i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            //            if (m<l) continue;
            //            if (m>l) continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              //              if (n<m) continue;
              //              if (n>m) continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l) % 2 != 0)
                continue;
              if ((oi.tz2 + oj.tz2 + ok.tz2) != (ol.tz2 + om.tz2 + on.tz2))
                continue;

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
              for (int m_i = oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ((m_i + m_j + m_k) != (m_l + m_m + m_n))
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          if (not(m_i == 1 and m_j == -1 and m_k == 1 and m_l == 1 and m_m == -1 and m_n == 1))
                            continue;
                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (size_t b = 0; b < norb; b++)
                            {
                              Orbit &ob = X.modelspace->GetOrbit(b);
                              for (size_t c = 0; c < norb; c++)
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
                                Orbit &oc = X.modelspace->GetOrbit(c);
                                double occfactor = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
                                for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                                {
                                  for (int m_b = -ob.j2; m_b <= ob.j2; m_b += 2)
                                  {
                                    for (int m_c = -oc.j2; m_c <= oc.j2; m_c += 2)
                                    {

                                      double xabklmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, k, m_k, l, m_l, m, m_m, c, m_c); // direct
                                      double xijcabn = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, n, m_n); // direct

                                      double xabilmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, i, m_i, l, m_l, m, m_m, c, m_c); // Pik
                                      double xkjcabn = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, n, m_n); // Pik
                                      double xabjlmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, j, m_j, l, m_l, m, m_m, c, m_c); // Pjk
                                      double xikcabn = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, n, m_n); // Pjk
                                      double xabknmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, k, m_k, n, m_n, m, m_m, c, m_c); // Pln
                                      double xijcabl = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, l, m_l); // Pln
                                      double xabklnc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, k, m_k, l, m_l, n, m_n, c, m_c); // Pmn
                                      double xijcabm = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, m, m_m); // Pmn

                                      double xabinmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, i, m_i, n, m_n, m, m_m, c, m_c); // Pik Pln
                                      double xkjcabl = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, l, m_l); // Pik Pln
                                      double xabjnmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, j, m_j, n, m_n, m, m_m, c, m_c); // Pjk Pln
                                      double xikcabl = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, l, m_l); // Pjk Pln
                                      double xabilnc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, i, m_i, l, m_l, n, m_n, c, m_c); // Pik Pmn
                                      double xkjcabm = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, m, m_m); // Pik Pmn
                                      double xabjlnc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, j, m_j, l, m_l, n, m_n, c, m_c); // Pjk Pmn
                                      double xikcabm = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, m, m_m); // Pjk Pmn

                                      double yabklmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, k, m_k, l, m_l, m, m_m, c, m_c); // direct
                                      double yijcabn = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, n, m_n); // direct

                                      double yabilmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, i, m_i, l, m_l, m, m_m, c, m_c); // Pik
                                      double ykjcabn = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, n, m_n); // Pik
                                      double yabjlmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, j, m_j, l, m_l, m, m_m, c, m_c); // Pjk
                                      double yikcabn = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, n, m_n); // Pjk
                                      double yabknmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, k, m_k, n, m_n, m, m_m, c, m_c); // Pln
                                      double yijcabl = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, l, m_l); // Pln
                                      double yabklnc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, k, m_k, l, m_l, n, m_n, c, m_c); // Pmn
                                      double yijcabm = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, m, m_m); // Pmn

                                      double yabinmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, i, m_i, n, m_n, m, m_m, c, m_c); // Pik Pln
                                      double ykjcabl = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, l, m_l); // Pik Pln
                                      double yabjnmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, j, m_j, n, m_n, m, m_m, c, m_c); // Pjk Pln
                                      double yikcabl = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, l, m_l); // Pjk Pln
                                      double yabilnc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, i, m_i, l, m_l, n, m_n, c, m_c); // Pik Pmn
                                      double ykjcabm = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, m, m_m); // Pik Pmn
                                      double yabjlnc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, j, m_j, l, m_l, n, m_n, c, m_c); // Pjk Pmn
                                      double yikcabm = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, m, m_m); // Pjk Pmn

                                      //        =   1/2 sum_abc [na*nb*(1-nc) + (1-na)*(1-nb)*nc ] *
                                      //                [ (Xabkcmn*Ycijabl - Yabkcmn*Xcijabl) - (Xabicmn*Yckjabl - Yabicmn*Xckjabl)  - (Xabjcmn*Ycikabl - Yabjcmn*Xcikabl)
                                      //                 -(Xabkcln*Ycijabm - Yabkcln*Xcijabm) + (Xabicln*Yckjabm - Yabicln*Xckjabm)  + (Xabjcln*Ycikabm - Yabjcln*Xcikabm)
                                      //                 -(Xabkcml*Ycijabn - Yabkcml*Xcijabn) + (Xabicml*Yckjabn - Yabicml*Xckjabn)  + (Xabjcml*Ycikabn - Yabjcml*Xcikabn) ]
                                      double dz = 0;
                                      dz += xabklmc * yijcabn - yabklmc * xijcabn; // Z1
                                      dz -= xabilmc * ykjcabn - yabilmc * xkjcabn; // Z2
                                      dz -= xabjlmc * yikcabn - yabjlmc * xikcabn; // Z3
                                      dz -= xabknmc * yijcabl - yabknmc * xijcabl; // Z4
                                      dz -= xabklnc * yijcabm - yabklnc * xijcabm; // Z5
                                      dz += xabinmc * ykjcabl - yabinmc * xkjcabl; // Z6
                                      dz += xabjnmc * yikcabl - yabjnmc * xikcabl; // Z7
                                      dz += xabilnc * ykjcabm - yabilnc * xkjcabm; // Z8
                                      dz += xabjlnc * yikcabm - yabjlnc * xikcabm; // Z9
                                      z_ijklmn += 0.5 * occfactor * dz;
                                      //                           z_ijklmn += 0.5 * occfactor * dz; // extra for ab<=>ba
                                      //                           if (std::abs(dz)>1e-6) std::cout << "  abc  " << a << " " << b << " "  << c << "  {m}  " << m_a << " " << m_b << " " << m_c
                                      //                              << "    x " << xabinmc << "  y " << ykjcabl << "    " << dz << "   " << 0.5*occfactor*dz <<  "    =>  " << z_ijklmn << std::endl;

                                    } // for m_c
                                  } // for m_b
                                } // for m_a
                              } // for c
                            } // for b
                          } // for a

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

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double ZJ_hermconj = GetMschemeMatrixElement_3b(Z_J, l, m_l, m, m_m, n, m_n, i, m_i, j, m_j, k, m_k);
                          double err = z_ijklmn - ZJ_ijklmn;
                          //                    std::cout << z_ijklmn << std::endl;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "Trouble in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << "  hc = " << ZJ_hermconj << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  //  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
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
bool UnitTest::Test_comm211sd(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.SetNonHermitian();
  Z_J.Erase();

  Commutator::comm211sd(X, Y, Z_J);

  //  std::cout << "X one body : " << std::endl << X.OneBody << std::endl;
  //  std::cout << "Y one body : " << std::endl << Y.OneBody << std::endl;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  //  size_t Q = Z_J.GetQSpaceOrbit();
  //  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);

  for (auto i : Z_J.modelspace->all_orbits)
  {
    Orbit &oi = Z_J.modelspace->GetOrbit(i);
    int mi = oi.j2;
    double Zm_i = 0;
    for (auto a : X.modelspace->all_orbits)
    {
      Orbit &oa = X.modelspace->GetOrbit(a);
      for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
      {
        double Xia = GetMschemeMatrixElement_1b(X, i, mi, a, ma);
        double Ya = GetMschemeMatrixElement_1leg(Y, a, ma);

        Zm_i += Xia * Ya;
      }
    }
    double ZJ_i = Z_J.OneBody(i, 0);
    //    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err) > 1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i << "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i, 0) << "   err = " << err << std::endl;
      //                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl;
    }
    summed_error += err * err;
    sum_m += Zm_i * Zm_i;
    sum_J += ZJ_i * ZJ_i;
  }

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
//  Z_i = sum_ab (na-nb)X_ab Y_bia
//
bool UnitTest::Test_comm231sd(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.SetNonHermitian();
  Z_J.Erase();

  Commutator::comm231sd(X, Y, Z_J);

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  //  size_t Q = Z_J.GetQSpaceOrbit();
  //  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);

  for (auto i : Z_J.modelspace->all_orbits)
  {
    Orbit &oi = Z_J.modelspace->GetOrbit(i);
    int mi = oi.j2;
    double Zm_i = 0;
    for (auto a : X.modelspace->all_orbits)
    {
      Orbit &oa = X.modelspace->GetOrbit(a);
      double na = oa.occ;
      for (auto b : X.modelspace->all_orbits)
      {
        Orbit &ob = X.modelspace->GetOrbit(b);
        double nb = ob.occ;
        for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
        {
          for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
          {
            double Xab = GetMschemeMatrixElement_1b(X, a, ma, b, mb);
            double Ybia = GetMschemeMatrixElement_3leg(Y, b, mb, i, mi, a, ma);

            Zm_i += (na - nb) * Xab * Ybia;
          }
        }
      }
    }
    double ZJ_i = Z_J.OneBody(i, 0);
    //    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err) > 1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i << "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i, 0) << "   err = " << err << std::endl;
      //                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl;
    }
    summed_error += err * err;
    sum_m += Zm_i * Zm_i;
    sum_J += ZJ_i * ZJ_i;
  }

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
//  Z_i = 1/2 sum_abc (na nb`nc`-na`nb nc) X_aibc Y_bca
//
bool UnitTest::Test_comm431sd(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.SetNonHermitian();
  Z_J.Erase();

  Commutator::comm433_pp_hh_431sd(X, Y, Z_J);

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  //  size_t Q = Z_J.GetQSpaceOrbit();
  //  Orbit& oQ = Z_J.modelspace->GetOrbit(Q);

  for (auto i : Z_J.modelspace->all_orbits)
  {
    Orbit &oi = Z_J.modelspace->GetOrbit(i);
    int mi = oi.j2;
    double Zm_i = 0;
    for (auto a : X.modelspace->all_orbits)
    {
      Orbit &oa = X.modelspace->GetOrbit(a);
      double na = oa.occ;
      for (auto b : X.modelspace->all_orbits)
      {
        Orbit &ob = X.modelspace->GetOrbit(b);
        double nb = ob.occ;
        for (auto c : X.modelspace->all_orbits)
        {
          Orbit &oc = X.modelspace->GetOrbit(c);
          double nc = oc.occ;
          for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
          {
            for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
            {
              for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
              {
                double Xaibc = GetMschemeMatrixElement_2b(X, a, ma, i, mi, b, mb, c, mc);
                double Ybca = GetMschemeMatrixElement_3leg(Y, b, mb, c, mc, a, ma);

                Zm_i += 1. / 2 * (na * (1 - nb) * (1 - nc) + (1 - na) * nb * nc) * Xaibc * Ybca;
              }
            }
          }
        }
      }
    }
    double ZJ_i = Z_J.OneBody(i, 0);
    //    double ZJ_i = Z_J.OneBody(i,Q);
    double err = Zm_i - ZJ_i;
    if (std::abs(err) > 1e-6)
    {
      std::cout << "Trouble in " << __func__ << "  i = " << i << "   Zm_i = " << Zm_i
                << "   ZJ_i = " << Z_J.OneBody(i, 0) << "   err = " << err << std::endl;
      //                << "   ZJ_i = " << Z_J.OneBody(i,Q) << "   err = " << err << std::endl;
    }
    summed_error += err * err;
    sum_m += Zm_i * Zm_i;
    sum_J += ZJ_i * ZJ_i;
  }

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
//
// 413:
// Zijk = sum_a Xijka*Ya
//
bool UnitTest::Test_comm413sd(const Operator &Xin, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

  Operator X = Xin;
  X.EraseOneBody(); // we do this so that the 233 term does not contribute.

  Commutator::comm413_233sd(X, Y, Z_J);

  size_t Q = Z_J.GetQSpaceOrbit();
  Orbit &oQ = Z_J.modelspace->GetOrbit(Q);

  int mQ = oQ.j2;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        //          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
        //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

        for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
        {
          for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
          {
            if ((mi + mj) != (mk + mQ))
              continue;
            //             if ( std::abs(mk)>ol.j2 ) continue;
            //             int ml = mi + mj - mk;

            double Zm_ijk = 0;
            for (auto a : X.modelspace->all_orbits)
            {
              Orbit &oa = X.modelspace->GetOrbit(a);
              //              double na = oa.occ;

              //                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
              //                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

              for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
              {
                if (ma != mQ)
                  continue;

                double Xijka = GetMschemeMatrixElement_2b(X, i, mi, j, mj, k, mk, a, ma);
                double Ya = GetMschemeMatrixElement_1leg(Y, a, ma);

                Zm_ijk += Xijka * Ya; // 413 term
                                      //                  if (i==0 and j==8 and k==0)
                                      //                  {
                                      //                    std::cout << "  #### a = " << a << "  m vals " << mi << " " << mj << " " << mk << " " << ma << "  " << mQ << "   Xijka , Ya = " << Xijka << " " << Ya << "   Zm_ijk = " << Zm_ijk << std::endl;
                                      //                  }

              } // for ma
            } // for a

            double ZJ_ijk = GetMschemeMatrixElement_3leg(Z_J, i, mi, j, mj, k, mk);
            double err = Zm_ijk - ZJ_ijk;
            if (std::abs(err) > 1e-6)
            {
              std::cout << "Trouble in " << __func__ << "  i,j,k = " << i << " " << j << " " << k
                        << "  mvals " << mi << " " << mj << " " << mk << "   " << mQ
                        << "   Zm_ijk = " << Zm_ijk << "   ZJ_ijk = " << ZJ_ijk << "   err = " << err << std::endl;
            }
            summed_error += err * err;
            sum_m += Zm_ijk * Zm_ijk;
            sum_J += ZJ_ijk * ZJ_ijk;

          } // for mk
        } // for mj
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.ThreeLegNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// 233:
// Zijk =  sum_a ( Xia*Yajk + Xja*Yjak - Yija * Xak )
//
bool UnitTest::Test_comm233sd(const Operator &X, const Operator &Yin)
{

  Operator Y = Yin;

  Operator Z_J(Y);
  Z_J.Erase();

  Y.EraseOneBody(); // we do this so the 413 term does not contribute

  Commutator::comm413_233sd(X, Y, Z_J);

  //  if ( Z_J.IsHermitian() )
  //     Z_J.Symmetrize();
  //  else if (Z_J.IsAntiHermitian() )
  //     Z_J.AntiSymmetrize();

  size_t Q = Z_J.GetQSpaceOrbit();
  Orbit &oQ = Z_J.modelspace->GetOrbit(Q);

  int mQ = oQ.j2;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        //          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
        //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

        for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
        {
          for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
          {
            if ((mi + mj) != (mk + mQ))
              continue;
            //             if ( std::abs(mk)>ol.j2 ) continue;
            //             int ml = mi + mj - mk;

            double Zm_ijk = 0;
            for (auto a : X.modelspace->all_orbits)
            {
              Orbit &oa = X.modelspace->GetOrbit(a);
              //              double na = oa.occ;

              //                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
              //                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

              for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
              {
                double Xia = GetMschemeMatrixElement_1b(X, i, mi, a, ma);
                double Xja = GetMschemeMatrixElement_1b(X, j, mj, a, ma);
                double Xak = GetMschemeMatrixElement_1b(X, a, ma, k, mk);

                double Yajk = GetMschemeMatrixElement_3leg(Y, a, ma, j, mj, k, mk);
                double Yiak = GetMschemeMatrixElement_3leg(Y, i, mi, a, ma, k, mk);
                double Yija = GetMschemeMatrixElement_3leg(Y, i, mi, j, mj, a, ma);

                Zm_ijk += Xia * Yajk + Xja * Yiak - Yija * Xak; // 223 term

              } // for ma
            } // for a

            double ZJ_ijk = GetMschemeMatrixElement_3leg(Z_J, i, mi, j, mj, k, mk);
            double err = Zm_ijk - ZJ_ijk;
            if (std::abs(err) > 1e-6)
            {
              std::cout << "Trouble in " << __func__ << "  i,j,k = " << i << " " << j << " " << k
                        << "  mvals " << mi << " " << mj << " " << mk
                        << "   Zm_ijk = " << Zm_ijk << "   ZJ_ijk = " << ZJ_ijk << "   err = " << err << std::endl;
            }
            summed_error += err * err;
            sum_m += Zm_ijk * Zm_ijk;
            sum_J += ZJ_ijk * ZJ_ijk;

          } // for mk
        } // for mj
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.ThreeLegNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// 433 pp/hh:
// Zijk =  1/2 sum_ab (nanb - na`nb`)  Xijab*Yabk
//
bool UnitTest::Test_comm433_pp_hh_sd(const Operator &X, const Operator &Yin)
{

  Operator Y = Yin;

  Operator Z_J(Y);
  Z_J.Erase();

  Commutator::comm433_pp_hh_431sd(X, Y, Z_J);

  size_t Q = Z_J.GetQSpaceOrbit();
  Orbit &oQ = Z_J.modelspace->GetOrbit(Q);

  int mQ = oQ.j2;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        //          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
        //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

        for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
        {
          for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
          {
            if ((mi + mj) != (mk + mQ))
              continue;
            //             if ( std::abs(mk)>ol.j2 ) continue;
            //             int ml = mi + mj - mk;

            double Zm_ijk = 0;
            for (auto a : X.modelspace->all_orbits)
            {
              Orbit &oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
              for (auto b : X.modelspace->all_orbits)
              {
                Orbit &ob = X.modelspace->GetOrbit(b);
                double nb = ob.occ;

                //                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
                //                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

                for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                {
                  for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                  {
                    double Xijab = GetMschemeMatrixElement_2b(X, i, mi, j, mj, a, ma, b, mb);
                    double Yabk = GetMschemeMatrixElement_3leg(Y, a, ma, b, mb, k, mk);

                    Zm_ijk += 1. / 2 * ((1 - na) * (1 - nb) - na * nb) * Xijab * Yabk;

                  } // for mb
                } // for ma
              } // for b
            } // for a

            double ZJ_ijk = GetMschemeMatrixElement_3leg(Z_J, i, mi, j, mj, k, mk);
            double err = Zm_ijk - ZJ_ijk;
            if (std::abs(err) > 1e-6)
            {
              std::cout << "Trouble in " << __func__ << "  i,j,k = " << i << " " << j << " " << k
                        << "  mvals " << mi << " " << mj << " " << mk
                        << "   Zm_ijk = " << Zm_ijk << "   ZJ_ijk = " << ZJ_ijk << "   err = " << err << std::endl;
            }
            summed_error += err * err;
            sum_m += Zm_ijk * Zm_ijk;
            sum_J += ZJ_ijk * ZJ_ijk;

          } // for mk
        } // for mj
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.ThreeLegNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// 433 pp/hh:
// Zijk =  sum_ab (na - nb)  Xiakb Ybja
//
bool UnitTest::Test_comm433sd_ph(const Operator &X, const Operator &Yin)
{

  Operator Y = Yin;

  Operator Z_J(Y);
  Z_J.Erase();

  Commutator::comm433sd_ph(X, Y, Z_J);
  //  Commutator::comm433sd_ph_dumbway( X, Y, Z_J);

  size_t Q = Z_J.GetQSpaceOrbit();
  Orbit &oQ = Z_J.modelspace->GetOrbit(Q);

  int mQ = oQ.j2;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        //          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
        //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;

        //          if (not (i==0 and j==1 and k==1) ) continue;

        for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
        {
          for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
          {
            if ((mi + mj) != (mk + mQ))
              continue;
            //             if ( std::abs(mk)>ol.j2 ) continue;
            //             int ml = mi + mj - mk;

            double Zm_ijk = 0;
            for (auto a : X.modelspace->all_orbits)
            {
              Orbit &oa = X.modelspace->GetOrbit(a);
              double na = oa.occ;
              for (auto b : X.modelspace->all_orbits)
              {
                Orbit &ob = X.modelspace->GetOrbit(b);
                double nb = ob.occ;

                //                  if ( not (a==0 and b==10) or (a==10 and b==0) ) continue;
                //                  if ( not ((a==1 and b==8) or (a==8 and b==1)) ) continue;
                //                if ( (oi.l+oj.l+oa.l+ob.l)%2>0) continue;
                //                if ( (oi.tz2+oj.tz2) != (oa.tz2+ob.tz2) ) continue;

                for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                {
                  for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                  {
                    double Xiakb = GetMschemeMatrixElement_2b(X, i, mi, a, ma, k, mk, b, mb);
                    double Ybja = GetMschemeMatrixElement_3leg(Y, b, mb, j, mj, a, ma);
                    double Xjakb = GetMschemeMatrixElement_2b(X, j, mj, a, ma, k, mk, b, mb);
                    double Ybia = GetMschemeMatrixElement_3leg(Y, b, mb, i, mi, a, ma);

                    Zm_ijk += (na - nb) * (Xiakb * Ybja - Xjakb * Ybia);
                    //                    if (i==0 and j==1 and k==1 and mi==1 and mj==-1 and mk==-1 and ( (a==0 and b==10) or (b==0 and a==10) ) )
                    //                    if (i==0 and j==1 and k==1 and mi==1 and mj==-1 and mk==-1 and std::abs(na-nb)>0 )
                    //                    if (i==0 and j==1 and k==1 and mi==1 and mj==-1 and mk==-1 and ( (a==1 and b==8) or (b==1 and a==8) ) )
                    //                    {
                    //                      std::cout << "  a,b " << a << " " << b << "   ma,mb " << ma << " " << mb << "   na nb  " << na << " " << nb
                    //                                << "  Xiakb Ybja " << Xiakb << " " << Ybja << "   Xjakb Ybia " << Xjakb << " " << Ybia << "  => Zm_ijk = " << Zm_ijk << std::endl;
                    //                    }
                    //                    Zm_ijk -= (na-nb) * (     Xjakb * Ybia ) ;

                  } // for mb
                } // for ma
              } // for b
            } // for a

            double ZJ_ijk = GetMschemeMatrixElement_3leg(Z_J, i, mi, j, mj, k, mk);
            double err = Zm_ijk - ZJ_ijk;
            if (std::abs(err) > 1e-6)
            {
              std::cout << "Trouble in " << __func__ << "  i,j,k = " << i << " " << j << " " << k
                        << "  mvals " << mi << " " << mj << " " << mk
                        << "   Zm_ijk = " << Zm_ijk << "   ZJ_ijk = " << ZJ_ijk << "   err = " << err << std::endl;
            }
            summed_error += err * err;
            sum_m += Zm_ijk * Zm_ijk;
            sum_J += ZJ_ijk * ZJ_ijk;

          } // for mk
        } // for mj
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "PASS " : "FAIL";
  if (Z_J.ThreeLegNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

bool UnitTest::TestRPAEffectiveCharge(const Operator &H, const Operator &OpIn, size_t a, size_t b)
{
  bool passed = true;

  int L = OpIn.GetJRank();
  int M = 0;
  double omega = H.OneBody(a, a) - H.OneBody(b, b);

  // number of m-scheme orbits
  int n_mscheme = 0;
  for (auto k : H.modelspace->all_orbits)
    n_mscheme += H.modelspace->GetOrbit(k).j2 + 1;

  std::cout << "n_mscheme = " << n_mscheme << std::endl;
  // Make column vector Tkl and Tkleff
  // and Matrix Mklpq
  arma::vec Tkl(n_mscheme * n_mscheme, arma::fill::zeros);
  arma::vec Tkl_eff(n_mscheme * n_mscheme, arma::fill::zeros);
  arma::mat Mklpq(n_mscheme * n_mscheme, n_mscheme * n_mscheme, arma::fill::zeros);

  size_t index_ab = 0;
  size_t index_kl = 0;
  for (auto k : H.modelspace->all_orbits)
  {
    Orbit &ok = H.modelspace->GetOrbit(k);
    double jk = 0.5 * ok.j2;
    //   std::cout << "  k = " << k << std::endl;
    for (auto l : H.modelspace->all_orbits)
    {
      Orbit &ol = H.modelspace->GetOrbit(l);
      double jl = 0.5 * ol.j2;

      //    std::cout << "   l = " << l << std::endl;
      for (int twomk = -ok.j2; twomk <= ok.j2; twomk += 2)
      {
        for (int twoml = -ol.j2; twoml <= ol.j2; twoml += 2)
        {
          //        std::cout << "index_kl = " << index_kl << std::endl;
          double mk = 0.5 * twomk;
          double ml = 0.5 * twoml;
          // Operator OpIn stores reduced matrix elements, so use Wigner Eckart
          Tkl(index_kl) = OpIn.OneBody(k, l) * AngMom::CG(jl, ml, L, M, jk, mk) / sqrt(2 * jk + 1.);
          //        std::cout << "Tkl is " << Tkl(index_kl) << std::endl;
          if ((k == a) and (l == b) and (twomk == 1) and (twoml == 1))
            index_ab = index_kl; // So we can go back and read what we want

          size_t index_pq = 0;
          // construct the matrix Mklpq = vkqlp (nq-np)/(omega-ep+eq)
          for (auto p : H.modelspace->all_orbits)
          {
            Orbit &op = H.modelspace->GetOrbit(p);
            double ep = H.OneBody(p, p);
            double np = op.occ;
            for (auto q : H.modelspace->all_orbits)
            {
              //             std::cout << "    p,q = " << p << " , " << q << std::endl;
              Orbit &oq = H.modelspace->GetOrbit(q);
              double eq = H.OneBody(q, q);
              double nq = oq.occ;
              for (int twomp = -op.j2; twomp <= op.j2; twomp += 2)
              {
                for (int twomq = -oq.j2; twomq <= oq.j2; twomq += 2)
                {
                  //                  double mp=0.5*twomp;
                  //                  double mq=0.5*twomq;
                  //                  double vkqlp = GetMschemeMatrixElement_2b( H, k,mk, q,mq, l,ml, p,mp );
                  //                  std::cout << "      about to get vkqlp " << "   index_pq = " << index_pq << std::endl;
                  if (std::abs(np - nq) > 1e-2)
                  {
                    double vkqlp = GetMschemeMatrixElement_2b(H, k, twomk, q, twomq, l, twoml, p, twomp);
                    Mklpq(index_kl, index_pq) = vkqlp * (nq - np) / (omega - ep + eq);
                  }
                  //                  std::cout <<"  and set it." << std::endl;
                  index_pq++;
                } // for twomq
              } // for twomp
            } // for q
          } // for p

          index_kl++;
        } // for twoml
      } // for twomk
    } // for l
  } // for k

  for (int iter = 0; iter < 12; iter++)
  {
    Tkl_eff = Tkl + Mklpq * Tkl_eff;
    std::cout << "iter = " << iter << "   Tkl = " << Tkl(index_ab) << "   Tkleff = " << Tkl_eff(index_ab) << "  => e = " << std::setprecision(8) << Tkl_eff(index_ab) / Tkl(index_ab) << std::endl;
  }

  return passed;
}

bool UnitTest::TestFactorizedDoubleCommutators()
{
  bool passed = true;

  int jrank = 0;
  int tz = 0;
  int parity = 0;
  int particle_rank = 2;

  Operator eta = RandomOp(*modelspace, jrank, tz, parity, particle_rank, -1);
  Operator H = RandomOp(*modelspace, jrank, tz, parity, particle_rank, +1);
  Operator OpOut_direct(*modelspace, jrank, tz, parity, 3);
  Operator OpOut_factorized(*modelspace, jrank, tz, parity, 2);
  OpOut_direct.ThreeBody.SetMode("pn");


  Commutator::comm223ss(eta, H, OpOut_direct);
  Commutator::comm231ss(eta, OpOut_direct, OpOut_direct);
  Commutator::comm232ss(eta, OpOut_direct, OpOut_direct);

  OpOut_direct.ThreeBody.Erase();
  // OpOut_factorized.EraseOneBody();  
  // OpOut_factorized.TwoBody.Erase();

  Commutator::FactorizedDoubleCommutator::SetUse_1b_Intermediates(true);
  Commutator::FactorizedDoubleCommutator::SetUse_2b_Intermediates(true);

  Commutator::FactorizedDoubleCommutator::comm223_231(eta, H, OpOut_factorized);
  Commutator::FactorizedDoubleCommutator::comm223_232(eta, H, OpOut_factorized);

  std::cout << "Norm of OpOut_direct:     " << OpOut_direct.Norm()     << "  1b : " << OpOut_direct.OneBodyNorm()     << "  2b : " << OpOut_direct.TwoBodyNorm() << std::endl;
  std::cout << "Norm of OpOut_factorized: " << OpOut_factorized.Norm() << "  1b : " << OpOut_factorized.OneBodyNorm() << "  2b : " << OpOut_factorized.TwoBodyNorm() << std::endl;

  OpOut_direct -= OpOut_factorized;
  passed = OpOut_direct.Norm() < 1e-6;
  if (not passed)
  {
    std::cout << __func__ << "  Uh Oh. Norm of difference is " << OpOut_direct.Norm()
              << "  1b part: " << OpOut_direct.OneBodyNorm()   << "  2b part: " << OpOut_direct.TwoBodyNorm()
              << std::endl;
  }

  return passed;
}



bool UnitTest::TestPerturbativeTriples()
{
  bool passed = true;

  Operator H = RandomOp( *modelspace, 0,0,0, 2, +1);
  Operator Omega = RandomOp( *modelspace, 0,0,0, 2, -1);
  H *= 1e3/H.Norm();
  Omega *= 1.0/Omega.Norm();

  IMSRGSolver imsrgsolver( H );
  imsrgsolver.SetGenerator("white");
  imsrgsolver.SetDenominatorPartitioning("Moller_Plesset");

  Omega += 10* imsrgsolver.generator.GetHod(Omega);
  imsrgsolver.SetOmega(0, Omega);

  double triples = imsrgsolver.CalculatePerturbativeTriples();

  /// Now we do it out explicitly as a cross check.
  /// As described in PRC 110, 044316 (2024), section IIIB
  /// The triples correction is
  /// Delta E[3] = 1/2 [Omegabar, Wbar]_0
  /// with
  /// Wbar = [Omega_2, Htilde_2]_3  (eq 26)
  /// Htilde = sum_k 1/(k+1)! [Omega,H](k)  (eq 27)
  /// Omegabar = Wbar^od / Deltabar (eq 21)
  BCH::SetBCHSkipiEq1(true);
  Operator Htilde = imsrgsolver.Transform(H);
  BCH::SetBCHSkipiEq1(false);

  Operator W = H; // We will use the 1b part of H in the denominators below
  W.ZeroBody =0;
  W.ThreeBody.SetMode("pn");
  W.SetParticleRank(3);
  Commutator::comm223ss( Omega, Htilde, W);


  Operator Omegabar = W*0;
  Omegabar.SetAntiHermitian();

  imsrgsolver.generator.Update(W,Omegabar); // Omegabar = Wbar_od / Deltabar

  Commutator::comm330ss(Omegabar,W,W);
  W.ZeroBody *=0.5;
  double diff = W.ZeroBody - triples;
  passed = std::abs( diff ) < 1e-6;
  if (not passed)
  {
    std::cout << __func__ << " Uh oh. Using IMSRGSolver.CalculatePerturbativeTriples() I get " << triples
              << " , but calculating directly with commutators I get " << W.ZeroBody << "     diff = " << diff
              << std::endl;
  }
  else
  {
    std::cout << "W0 = " << W.ZeroBody << "  triples = " << triples << " difference = " << diff << std::endl;
  }

  return passed;
}

bool UnitTest::SanityCheck()
{

  std::cout << "BUILD VERSION = " << version::BuildVersion() << std::endl;
  std::cout << "Test simple Clebsch-Gordan coeff..." << std::endl;
  double cg1 = AngMom::CG(0.5, 0.5, 0.5, -0.5, 0, 0);
  if (std::abs(cg1 - sqrt(0.5)) > 1e-6)
  {
    std::cout << __FILE__ << "  " << __func__ << " failed on line " << __LINE__ << std::endl;
    return false;
  }

  std::cout << "Construct a model space..." << std::endl;
  int emax = 2;
  std::string ref = "He4";
  auto ms = ModelSpace(emax, ref, ref);
  //  int A,Z;
  double A, Z;
  ms.GetAZfromString("Pb208", A, Z);
  //  if ( not (A==208 and Z==82) )
  if ((std::abs(A - 208.) > 1e-10) or (std::abs(Z - 82.) > 1e-10))
  {
    std::cout << __FILE__ << "  " << __func__ << " failed on line " << __LINE__ << std::endl;
    return false;
  }

  // try making a mixed reference
  double amix, zmix;
  ms.GetAZfromString("mix_A4.5_Z2.5", amix, zmix);
  if ((std::abs(amix - 4.5) > 1e-10) or (std::abs(zmix - 2.5) > 1e-10))
  {
    std::cout << __FILE__ << "  " << __func__ << " failed on line " << __LINE__ << std::endl;
    std::cout << "amix = " << amix << " (should be 4.5)  zmix = " << zmix << " (should be 2.5) " << std::endl;
    return false;
  }

  std::cout << "Construct the kinetic energy operator..." << std::endl;
  Operator trel = imsrg_util::Trel_Op(ms);
  double normT = trel.Norm();
  std::cout << "...it should have a non-zero norm...  norm = " << normT << std::endl;
  if (std::abs(normT) < 1e-8)
  {
    std::cout << __FILE__ << "  " << __func__ << " failed on line " << __LINE__ << std::endl;
    return false;
  }
  Operator comTT = Commutator::Commutator(trel, trel);
  double normcomTT = comTT.Norm();
  std::cout << "...and check that it commutes with itself...   || [T,T] || = " << normcomTT << std::endl;
  if (std::abs(normcomTT) > 1e-8)
  {
    std::cout << __FILE__ << "  " << __func__ << " failed on line " << __LINE__ << std::endl;
    return false;
  }

  std::cout << __func__ << " :  Things look ok! " << std::endl;
  return true;
}

/// M-Scheme Formula:
//
// Z_ij = 1/12 sum_abcde (nanb n`c n`dn`e + n`an`b ncndne) (X_abicde Y_cdeabj - Y_abicde X_cdeabj)
//
bool UnitTest::Mscheme_Test_comm331st(const Operator &X, const Operator &Y)
{

  Operator Z_J(Y);
  Z_J.Erase();

  int Lambda = Y.GetJRank();
  int TRank = Y.GetTRank();

  Operator Y_copy(Y);
  if (Lambda == 0 and TRank == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }

  // Commutator::comm331st( X, Y, Z_J);
  ReferenceImplementations::comm331st(X, Y_copy, Z_J);

  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  std::cout << __func__ << "   Done with J-scheme commutator. Begin m-scheme version..." << std::endl;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    mi = 1;
    for (auto j : Z_J.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);

      if (std::abs(oj.j2 - oi.j2) > Lambda * 2 or (oi.j2 + oj.j2) < Lambda * 2)
        continue;
      for (int mi = -oi.j2; mi <= oi.j2; mi += 2)
      {
        for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
        {
          int twoTm = mi - mj;
          if (std::abs(twoTm) > 2 * Lambda)
            continue;

          double Zm_ij = 0;
          size_t norb = X.modelspace->GetNumberOrbits();

          // #pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ij)
          for (size_t a = 0; a < norb; a++)
          {
            Orbit &oa = X.modelspace->GetOrbit(a);
            double na = oa.occ;

            for (auto b : X.modelspace->all_orbits)
            {
              Orbit &ob = X.modelspace->GetOrbit(b);
              double nb = ob.occ;

              for (auto c : X.modelspace->all_orbits)
              {
                Orbit &oc = X.modelspace->GetOrbit(c);
                double nc = oc.occ;

                for (auto d : X.modelspace->all_orbits)
                {
                  Orbit &od = X.modelspace->GetOrbit(d);
                  double nd = od.occ;

                  for (auto e : X.modelspace->all_orbits)
                  {
                    Orbit &oe = X.modelspace->GetOrbit(e);
                    double ne = oe.occ;

                    // Tensor operator may change parity. So orbit i and j may not have the same parity
                    if (
                        ((oa.l + ob.l + oi.l + oc.l + od.l + oe.l + X.GetParity()) % 2 != 0) and ((oa.l + ob.l + oj.l + oc.l + od.l + oe.l + Y.GetParity()) % 2 != 0) and
                        ((oa.l + ob.l + oi.l + oc.l + od.l + oe.l + Y.GetParity()) % 2 != 0) and ((oa.l + ob.l + oj.l + oc.l + od.l + oe.l + X.GetParity()) % 2 != 0))
                      continue;

                    if (
                        (std::abs(oa.tz2 + ob.tz2 + oi.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * X.GetTRank()) and (std::abs(oa.tz2 + ob.tz2 + oj.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * Y.GetTRank()) and
                        (std::abs(oa.tz2 + ob.tz2 + oi.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * Y.GetTRank()) and (std::abs(oa.tz2 + ob.tz2 + oj.tz2 - oc.tz2 - od.tz2 - oe.tz2) != 2 * X.GetTRank()))
                      continue;

                    double occfactor = na * nb * (1 - nc) * (1 - nd) * (1 - ne) + (1 - na) * (1 - nb) * nc * nd * ne;
                    if (std::abs(occfactor) < 1e-7)
                      continue;

                    for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                    {
                      for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                      {
                        for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                        {
                          for (int md = -od.j2; md <= od.j2; md += 2)
                          {
                            for (int me = -oe.j2; me <= oe.j2; me += 2)
                            {
                              if (std::abs(ma + mb + mi - mc - md - me) <= 2 * X.GetJRank() and std::abs(ma + mb + mj - mc - md - me) <= 2 * Y.GetJRank() and (mc + md + me - ma - mb - mj) == twoTm)
                              {
                                double xabicde = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, c, mc, d, md, e, me);
                                double ycdeabj = GetMschemeMatrixElement_3b(Y, c, mc, d, md, e, me, a, ma, b, mb, j, mj);
                                Zm_ij += 1. / 12 * occfactor * (xabicde * ycdeabj);
                              }

                              if (std::abs(ma + mb + mi - mc - md - me) <= 2 * Y.GetJRank() and std::abs(ma + mb + mj - mc - md - me) <= 2 * X.GetJRank() and ma + mb + mi - mc - md - me == twoTm)
                              {
                                double xcdeabj = GetMschemeMatrixElement_3b(X, c, mc, d, md, e, me, a, ma, b, mb, j, mj);
                                double yabicde = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, i, mi, c, mc, d, md, e, me);
                                Zm_ij -= 1. / 12 * occfactor * (yabicde * xcdeabj);
                              }

                            } // for me
                          } // for md
                        } // for mc
                      } // for mb
                    } // for ma
                  } // for e
                } // for d
              } // for c
            } // for b
          } // for a

          double ZJ_ij = GetMschemeMatrixElement_1b(Z_J, i, mi, j, mj);
          
          double err = Zm_ij - ZJ_ij;

          if (std::abs(err) > 1e-6)
          {
            std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j = " << i << " " << j
                      << "   zij = " << Zm_ij << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl;
          }
          summed_error += err * err;
          sum_m += Zm_ij * Zm_ij;
          sum_J += ZJ_ij * ZJ_ij;
        }
      }
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-scheme formula
//
// Z_ijklmn = (1-Pik-Pjk)(1-Plm-Pln) sum_a (X_ijla Y_akmn - Y_ijla X_akmn)
//
bool UnitTest::Mscheme_Test_comm223st(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = (parityX + parityY) % 2;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;
  int JrankZ = X.GetJRank() + Y.GetJRank();
  int Lambda = Y.GetJRank();

  Operator Z_J(*(X.modelspace), JrankZ, TzZ, parityZ, 3);
  //  Operator Z_J( Y );
  Z_J.SetHermitian();
  //  Z_J.Erase();
  Z_J.ThreeBody.SetMode("pn");
  Z_J.ThreeBody.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and TzY == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }

  std::cout << "Calling J-scheme commutator" << std::endl;
  ReferenceImplementations::comm223st(X, Y_copy, Z_J);

  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;
  std::cout << "start mscheme loops " << std::endl;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j > i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k > j)
          continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l > i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            if (m > l)
              continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n > m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              //              if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2 !=0 ) continue;
              //              if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l + parityZ) % 2 != 0)
                continue;
              if (std::abs(oi.tz2 + oj.tz2 + ok.tz2 - ol.tz2 - om.tz2 - on.tz2) != 2 * TzZ)
                continue;

              // loop over projections
              for (int m_i = -oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          int twoTm = m_i + m_j + m_k - m_l - m_m - m_n;
                          if (std::abs(twoTm) > 2 * Lambda)
                            continue;

                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;
                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                            {
                              // "direct" term  Z1
                              double x_ijla = GetMschemeMatrixElement_2b(X, i, m_i, j, m_j, l, m_l, a, m_a);
                              double y_akmn = GetMschemeMatrixElement_2b(Y, a, m_a, k, m_k, m, m_m, n, m_n);
                              double y_ijla = GetMschemeMatrixElement_2b(Y, i, m_i, j, m_j, l, m_l, a, m_a);
                              double x_akmn = GetMschemeMatrixElement_2b(X, a, m_a, k, m_k, m, m_m, n, m_n);
                              // Pik  Z2
                              double x_kjla = GetMschemeMatrixElement_2b(X, k, m_k, j, m_j, l, m_l, a, m_a);
                              double y_aimn = GetMschemeMatrixElement_2b(Y, a, m_a, i, m_i, m, m_m, n, m_n);
                              double y_kjla = GetMschemeMatrixElement_2b(Y, k, m_k, j, m_j, l, m_l, a, m_a);
                              double x_aimn = GetMschemeMatrixElement_2b(X, a, m_a, i, m_i, m, m_m, n, m_n);
                              // Pjk   Z3
                              double x_ikla = GetMschemeMatrixElement_2b(X, i, m_i, k, m_k, l, m_l, a, m_a);
                              double y_ajmn = GetMschemeMatrixElement_2b(Y, a, m_a, j, m_j, m, m_m, n, m_n);
                              double y_ikla = GetMschemeMatrixElement_2b(Y, i, m_i, k, m_k, l, m_l, a, m_a);
                              double x_ajmn = GetMschemeMatrixElement_2b(X, a, m_a, j, m_j, m, m_m, n, m_n);
                              // Plm   Z4
                              double x_ijma = GetMschemeMatrixElement_2b(X, i, m_i, j, m_j, m, m_m, a, m_a);
                              double y_akln = GetMschemeMatrixElement_2b(Y, a, m_a, k, m_k, l, m_l, n, m_n);
                              double y_ijma = GetMschemeMatrixElement_2b(Y, i, m_i, j, m_j, m, m_m, a, m_a);
                              double x_akln = GetMschemeMatrixElement_2b(X, a, m_a, k, m_k, l, m_l, n, m_n);
                              // Pln   Z5
                              double x_ijna = GetMschemeMatrixElement_2b(X, i, m_i, j, m_j, n, m_n, a, m_a);
                              double y_akml = GetMschemeMatrixElement_2b(Y, a, m_a, k, m_k, m, m_m, l, m_l);
                              double y_ijna = GetMschemeMatrixElement_2b(Y, i, m_i, j, m_j, n, m_n, a, m_a);
                              double x_akml = GetMschemeMatrixElement_2b(X, a, m_a, k, m_k, m, m_m, l, m_l);
                              // Pik Plm    Z6
                              double x_kjma = GetMschemeMatrixElement_2b(X, k, m_k, j, m_j, m, m_m, a, m_a);
                              double y_ailn = GetMschemeMatrixElement_2b(Y, a, m_a, i, m_i, l, m_l, n, m_n);
                              double y_kjma = GetMschemeMatrixElement_2b(Y, k, m_k, j, m_j, m, m_m, a, m_a);
                              double x_ailn = GetMschemeMatrixElement_2b(X, a, m_a, i, m_i, l, m_l, n, m_n);
                              // Pik Pln    Z7
                              double x_kjna = GetMschemeMatrixElement_2b(X, k, m_k, j, m_j, n, m_n, a, m_a);
                              double y_aiml = GetMschemeMatrixElement_2b(Y, a, m_a, i, m_i, m, m_m, l, m_l);
                              double y_kjna = GetMschemeMatrixElement_2b(Y, k, m_k, j, m_j, n, m_n, a, m_a);
                              double x_aiml = GetMschemeMatrixElement_2b(X, a, m_a, i, m_i, m, m_m, l, m_l);
                              // Pjk Plm    Z8
                              double x_ikma = GetMschemeMatrixElement_2b(X, i, m_i, k, m_k, m, m_m, a, m_a);
                              double y_ajln = GetMschemeMatrixElement_2b(Y, a, m_a, j, m_j, l, m_l, n, m_n);
                              double y_ikma = GetMschemeMatrixElement_2b(Y, i, m_i, k, m_k, m, m_m, a, m_a);
                              double x_ajln = GetMschemeMatrixElement_2b(X, a, m_a, j, m_j, l, m_l, n, m_n);
                              // Pjk Pln    Z9
                              double x_ikna = GetMschemeMatrixElement_2b(X, i, m_i, k, m_k, n, m_n, a, m_a);
                              double y_ajml = GetMschemeMatrixElement_2b(Y, a, m_a, j, m_j, m, m_m, l, m_l);
                              double y_ikna = GetMschemeMatrixElement_2b(Y, i, m_i, k, m_k, n, m_n, a, m_a);
                              double x_ajml = GetMschemeMatrixElement_2b(X, a, m_a, j, m_j, m, m_m, l, m_l);
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

                            } // for m_a
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;

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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  //  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
//  Z_ij = 1/4 sum_abcd (nanb n`cn`d) (X_abcd Y_cdiabj - Y_abicdj X_cdab)
//
bool UnitTest::Mscheme_Test_comm231st(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;
  int Lambda = Y.GetJRank();
  int Top = Y.GetJRank();

  Operator Z_J(Y);
  Z_J.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and Top == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }

  std::cout << "Calling J-scheme commutator" << std::endl;
  ReferenceImplementations::comm231st(X, Y_copy, Z_J);

  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  std::cout << __func__ << "   Done with J-scheme commutator. Begin m-scheme version..." << std::endl;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;
  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    int mi = oi.j2;
    for (auto j : Z_J.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
    {
      int mj = mi;
      double Zm_ij = 0;

      size_t norb = X.modelspace->GetNumberOrbits();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ij)
      for (size_t a = 0; a < norb; a++)
      {
        Orbit &oa = X.modelspace->GetOrbit(a);
        double na = oa.occ;
        if (na < 1e-6)
          continue;
        for (auto b : X.modelspace->all_orbits)
        {
          Orbit &ob = X.modelspace->GetOrbit(b);
          double nb = ob.occ;
          if (nb < 1e-6)
            continue;
          for (auto c : X.modelspace->all_orbits)
          {
            Orbit &oc = X.modelspace->GetOrbit(c);
            double nc = oc.occ;
            if ((1 - nc) < 1e-6)
              continue;
            for (auto d : X.modelspace->all_orbits)
            {
              Orbit &od = X.modelspace->GetOrbit(d);
              double nd = od.occ;
              if ((1 - nd) < 1e-6)
                continue;
              if (((oa.l + ob.l + oc.l + od.l + parityX) % 2 > 0) and ((oa.l + ob.l + oc.l + od.l + parityY) % 2 > 0))
                continue;
              if ((std::abs(oa.tz2 + ob.tz2 - oc.tz2 - od.tz2) != 2 * TzX) and (std::abs(oa.tz2 + ob.tz2 - oc.tz2 - od.tz2) != 2 * TzY))
                continue;
              for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
              {
                for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                {
                  for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                  {
                    for (int md = -od.j2; md <= od.j2; md += 2)
                    {
                      if ((ma + mb) != (mc + md))
                        continue;
                      double xabcd = GetMschemeMatrixElement_2b(X, a, ma, b, mb, c, mc, d, md);
                      double xcdab = GetMschemeMatrixElement_2b(X, c, mc, d, md, a, ma, b, mb);
                      double yabicdj = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, i, mi, c, mc, d, md, j, mj);
                      double ycdiabj = GetMschemeMatrixElement_3b(Y, c, mc, d, md, i, mi, a, ma, b, mb, j, mj);

                      double yabcd = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, c, mc, d, md);
                      double ycdab = GetMschemeMatrixElement_2b(Y, c, mc, d, md, a, ma, b, mb);
                      double xabicdj = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, c, mc, d, md, j, mj);
                      double xcdiabj = GetMschemeMatrixElement_3b(X, c, mc, d, md, i, mi, a, ma, b, mb, j, mj);

                      Zm_ij += (1. / 4) * na * nb * (1 - nc) * (1 - nd) * ((xabcd * ycdiabj - yabicdj * xcdab) - (yabcd * xcdiabj - xabicdj * ycdab));

                    } // for md
                  } // for mc
                } // for mb
              } // for ma

            } // for d
          } // for c
        } // for b
      } // for a

      double ZJ_ij = GetMschemeMatrixElement_1b(Z_J, i, mi, j, mj);
      double err = Zm_ij - ZJ_ij;
      if (std::abs(err) > 1e-6)
      {
        std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j = " << i << " " << j
                  << "   zij = " << Zm_ij << "   ZJ_ij = " << ZJ_ij << "   err = " << err << std::endl;
      }
      summed_error += err * err;
      sum_m += Zm_ij * Zm_ij;
      sum_J += ZJ_ij * ZJ_ij;

    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ijkl = -1/2 * sum_abc (n_a*n_b*nbar_c + nbar_a*nbar_b*n_c) * [  Xicab*Yabjklc - Xjcab*Yabiklc - Yijcabl*Xabkc + Yijcabk*Xablc
//                                                                 - Yicab*Xabjklc + Yjcab*Xabiklc + Xijcabl*Yabkc - Xijcabk*Yablc ]
//
bool UnitTest::Mscheme_Test_comm232st(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;
  int Lambda = Y.GetJRank();
  int Top = Y.GetJRank();
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and Top == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }
  // Commutator::comm232ss(X, Y, Z_J);
  ReferenceImplementations::comm232st(X, Y_copy, Z_J);
  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  //  if ( Z_J.IsHermitian() )
  //     Z_J.Symmetrize();
  //  else if (Z_J.IsAntiHermitian() )
  //     Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;
  size_t norb = X.modelspace->GetNumberOrbits();
  std::cout << "Begin m-scheme loops" << std::endl;

  for (auto i : X.modelspace->all_orbits)
  {
    std::cout << " orbit " << i << " of " << X.modelspace->all_orbits.size() << std::endl;
    Orbit &oi = X.modelspace->GetOrbit(i);
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k or (k == i and l < j))
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          //          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue; // check parity
          //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue; // check isospin projection
          //         if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue; // check parity
          //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue; // check isospin projection

          int m_max = std::min(oi.j2 + oj.j2, ok.j2 + ol.j2);
          //          for (int mi=-oi.j2; mi<=oi.j2; mi+=2)
          for (int mi = oi.j2; mi <= oi.j2; mi += 2)
          {
            for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
            {
              if (mi + mj < 0)
                continue;
              if (mi + mj != m_max)
                continue;
              for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
              {
                for (int ml = -ol.j2; ml <= ol.j2; ml += 2)
                {

                  double Zm_ijkl = 0;

                  //             for (auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ijkl)
                  for (size_t a = 0; a < norb; a++)
                  {
                    Orbit &oa = X.modelspace->GetOrbit(a);
                    double na = oa.occ;
                    for (auto b : X.modelspace->all_orbits)
                    {
                      Orbit &ob = X.modelspace->GetOrbit(b);
                      double nb = ob.occ;

                      for (auto c : X.modelspace->all_orbits)
                      {
                        Orbit &oc = X.modelspace->GetOrbit(c);
                        double nc = oc.occ;
                        double occfactor = na * nb * (1 - nc) + (1 - na) * (1 - nb) * nc;
                        if (std::abs(occfactor) < 1e-8)
                          continue;

                        for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                        {
                          for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                          {
                            for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                            {

                              double xicab = GetMschemeMatrixElement_2b(X, i, mi, c, mc, a, ma, b, mb);
                              double yicab = GetMschemeMatrixElement_2b(Y, i, mi, c, mc, a, ma, b, mb);
                              double xjcab = GetMschemeMatrixElement_2b(X, j, mj, c, mc, a, ma, b, mb);
                              double yjcab = GetMschemeMatrixElement_2b(Y, j, mj, c, mc, a, ma, b, mb);
                              double xabkc = GetMschemeMatrixElement_2b(X, a, ma, b, mb, k, mk, c, mc);
                              double yabkc = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, k, mk, c, mc);
                              double xablc = GetMschemeMatrixElement_2b(X, a, ma, b, mb, l, ml, c, mc);
                              double yablc = GetMschemeMatrixElement_2b(Y, a, ma, b, mb, l, ml, c, mc);

                              double xabjklc = GetMschemeMatrixElement_3b(X, a, ma, b, mb, j, mj, k, mk, l, ml, c, mc);
                              double yabjklc = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, j, mj, k, mk, l, ml, c, mc);
                              double xabiklc = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, k, mk, l, ml, c, mc);
                              double yabiklc = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, i, mi, k, mk, l, ml, c, mc);
                              double xijcabl = GetMschemeMatrixElement_3b(X, i, mi, j, mj, c, mc, a, ma, b, mb, l, ml);
                              double yijcabl = GetMschemeMatrixElement_3b(Y, i, mi, j, mj, c, mc, a, ma, b, mb, l, ml);
                              double xijcabk = GetMschemeMatrixElement_3b(X, i, mi, j, mj, c, mc, a, ma, b, mb, k, mk);
                              double yijcabk = GetMschemeMatrixElement_3b(Y, i, mi, j, mj, c, mc, a, ma, b, mb, k, mk);

                              Zm_ijkl += -0.5 * occfactor * (xicab * yabjklc - xjcab * yabiklc - yijcabl * xabkc + yijcabk * xablc - yicab * xabjklc + yjcab * xabiklc + xijcabl * yabkc - xijcabk * yablc);

                            } // for mc
                          } // for mb
                        } // for ma
                      } // for c
                    } // for b
                  } // for a

                  double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
                  double err = Zm_ijkl - ZJ_ijkl;
                  if (std::abs(err) > 1e-6)
                  {
                    std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                              << " {m} = " << mi << " " << mj << " " << mk << " " << ml
                              << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
                  }
                  summed_error += err * err;
                  sum_m += Zm_ijkl * Zm_ijkl;
                  sum_J += ZJ_ijkl * ZJ_ijkl;

                } // for ml
              } // for mk
            } // for mj
          } // for mi
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ijkl = P(ij/k) * sum_a  [ X_ka * Y_ijalmn - X_ijalmn * Y_ka ]
//          -P(lm/n) * sum_a [ X_an * Y_ijklma + X_ijklma * Y_an ]
//
bool UnitTest::Mscheme_Test_comm133st(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;
  int Lambda = Y.GetJRank();
  int Top = Y.GetJRank();
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.ThreeBody.SetMode("pn");
  Z_J.ThreeBody.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and Top == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }
//  Commutator::comm133st(X, Y_copy, Z_J);
  ReferenceImplementations::comm133st(X, Y_copy, Z_J);
  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

/*
  std::cout << "X one body : " << std::endl
            << X.OneBody << std::endl;
  std::cout << "Y one body : " << std::endl
            << Y.OneBody << std::endl;
*/

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzY;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < j)
          continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l < i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            if (m < l)
              continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n < m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);

              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l + parityY) % 2 != 0)
                continue;
              if (std::abs(oi.tz2 + oj.tz2 + ok.tz2 - ol.tz2 - om.tz2 - on.tz2) != 2 * TzZ)
                continue;
              // loop over projections
              for (int m_i = oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if (m_i + m_j + m_k < 0)
                      continue;
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          //if ((m_i + m_j + m_k) != (m_l + m_m + m_n))
                          //  continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          double z_ijklmn = 0;
                          int Tm = m_i + m_j + m_k - m_l - m_m - m_n;
                          if (std::abs(Tm) > Lambda * 2)
                            continue;

                          size_t norb = X.modelspace->GetNumberOrbits();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                            {
                              // Bra
                              double xka = 0;
                              double yijalmn = 0;
//                              if (oa.j2 == ok.j2)
                              {  
                                xka = GetMschemeMatrixElement_1b(X, k, m_k, a, m_a);
                                yijalmn = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, a, m_a, l, m_l, m, m_m, n, m_n);
                                z_ijklmn += xka * yijalmn;
                              }

                              double yka = 0;
                              double xijalmn = 0;
//                              if (oa.j2 == ok.j2 and  m_i + m_j + m_a - m_l - m_m - m_n == 0)
//                              if (  m_i + m_j + m_a - m_l - m_m - m_n == 0)
                              {
                                yka = GetMschemeMatrixElement_1b(Y, k, m_k, a, m_a);
                                xijalmn = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, a, m_a, l, m_l, m, m_m, n, m_n);
                                z_ijklmn -= yka * xijalmn;
                              }

                              double xja = 0;
                              double yiaklmn = 0;
//                              if (oa.j2 == oj.j2 and m_a == m_j and m_i + m_k + m_a - m_l - m_m - m_n == Tm)
                              {
                                xja = GetMschemeMatrixElement_1b(X, j, m_j, a, m_a);
                                yiaklmn = GetMschemeMatrixElement_3b(Y, i, m_i, a, m_a, k, m_k, l, m_l, m, m_m, n, m_n);
                                z_ijklmn += xja * yiaklmn;
                              }

                              double xiaklmn = 0;
                              double yja = 0;
//                              if (oa.j2 == oj.j2 and m_i + m_k + m_a - m_l - m_m - m_n == 0)
                              {
                                xiaklmn = GetMschemeMatrixElement_3b(X, i, m_i, a, m_a, k, m_k, l, m_l, m, m_m, n, m_n);
                                yja = GetMschemeMatrixElement_1b(Y, j, m_j, a, m_a);
                                z_ijklmn -= yja * xiaklmn;
                              }

                              double xia = 0;
                              double yajklmn = 0;
//                              if (oa.j2 == oi.j2 and m_a == m_i and m_j + m_k + m_a - m_l - m_m - m_n == Tm)
                              {  
                                xia = GetMschemeMatrixElement_1b(X, i, m_i, a, m_a);
                                yajklmn = GetMschemeMatrixElement_3b(Y, a, m_a, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                                z_ijklmn += xia * yajklmn;
                              }

                              double yia = 0;
                              double xajklmn = 0;
//                              if (oa.j2 == oi.j2 and m_j + m_k + m_a - m_l - m_m - m_n == 0)
                              {
                                xajklmn = GetMschemeMatrixElement_3b(X, a, m_a, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                                yia = GetMschemeMatrixElement_1b(Y, i, m_i, a, m_a);
                                z_ijklmn -= yia * xajklmn;
                              }

                              /// Ket
                              double xal = 0;
                              double yijkamn = 0;
//                              if (oa.j2 == ol.j2 and m_a == m_l and m_i + m_j + m_k - m_a - m_m - m_n == Tm)
                              { 
                                xal = GetMschemeMatrixElement_1b(X, a, m_a, l, m_l);
                                yijkamn = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, m, m_m, n, m_n);
                                z_ijklmn -= yijkamn * xal;
                              }

                              double yal = 0;
                              double xijkamn = 0;
//                              if (oa.j2 == ol.j2 and m_i + m_j + m_k - m_a - m_m - m_n == 0)
                              {
                                xijkamn = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, m, m_m, n, m_n);
                                yal = GetMschemeMatrixElement_1b(Y, a, m_a, l, m_l);
                                z_ijklmn += xijkamn * yal;
                              }

                              double xam = 0;
                              double yijklan = 0;
//                              if (oa.j2 == om.j2 and m_a == m_m and m_i + m_j + m_k - m_a - m_l - m_n == Tm)
                              { 
                                xam = GetMschemeMatrixElement_1b(X, a, m_a, m, m_m);
                                yijklan = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, l, m_l, a, m_a, n, m_n);
                                z_ijklmn -= yijklan * xam;
                              }

                              double yam = 0;
                              double xijklan = 0;
//                              if (oa.j2 == om.j2 and m_i + m_j + m_k - m_a - m_l - m_n == 0)
//                              if (oa.j2 == om.j2)
                              {
                                yam = GetMschemeMatrixElement_1b(Y, a, m_a, m, m_m);
                                xijklan = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, l, m_l, a, m_a, n, m_n);
                                z_ijklmn += xijklan * yam;
                              }

                              double xan = 0;
                              double yijklma = 0;
//                              if (oa.j2 == on.j2 and m_a == m_n and m_i + m_j + m_k - m_a - m_l - m_m == Tm)
                              {
                                xan = GetMschemeMatrixElement_1b(X, a, m_a, n, m_n);
                                yijklma = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, a, m_a);
                                z_ijklmn -= yijklma * xan;
                              }


                              double yan = 0;
                              double xijklma = 0;
//                              if (oa.j2 == on.j2 and m_i + m_j + m_k - m_a - m_l - m_m == 0)
                              {
                                xijklma = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, a, m_a);
                                yan = GetMschemeMatrixElement_1b(Y, a, m_a, n, m_n);
                                z_ijklmn += xijklma * yan;
                              }


                            } // for m_a
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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
  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  //  if ( Z_J.TwoBodyNorm() < 1e-6 ) std::cout << "WARNING " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ijkl = sum_ab (na-nb) (Xab * Yijbkla - Yab * Xijbkla)
//
bool UnitTest::Mscheme_Test_comm132st(const Operator &X, const Operator &Y)
{
  std::cout << __func__ << std::endl;
  int Lambda = Y.GetJRank();
  int Top = Y.GetJRank();
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and Top == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }
  // Commutator::comm132st(X, Y_copy, Z_J);
  ReferenceImplementations::comm132st(X, Y_copy, Z_J);
  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  //  if ( Z_J.IsHermitian() )
  //     Z_J.Symmetrize();
  //  else if (Z_J.IsAntiHermitian() )
  //     Z_J.AntiSymmetrize();

  //std::cout << "X one body : " << std::endl
  //          << X.OneBody << std::endl;
  //std::cout << "Y one body : " << std::endl
  //          << Y.OneBody << std::endl;

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = (parityX + parityY) % 2;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          //          if ( (oi.l+oj.l+ok.l+ol.l)%2>0) continue;
          //          if ( (oi.tz2+oj.tz2) != (ok.tz2+ol.tz2) ) continue;
          if ((oi.l + oj.l + ok.l + ol.l + parityY) % 2 > 0)
            continue;
          if (std::abs(oi.tz2 + oj.tz2 - ok.tz2 - ol.tz2) != TzZ * 2)
            continue;
          for (int mi = -oi.j2; mi <= oi.j2; mi += 2)
          {
            for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
            {
              for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
              {
                for (int ml = -ol.j2; ml <= ol.j2; ml += 2)
                {
                  int Tm = mi + mj - mk - ml;
                  if (std::abs(Tm) > 2 * Lambda or Tm < 0)
                    continue;

                  double Zm_ijkl = 0;
                  for (auto a : X.modelspace->all_orbits)
                  {
                    Orbit &oa = X.modelspace->GetOrbit(a);
                    double na = oa.occ;

                    for (auto b : X.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}))
                    {
//                      if (a == b)
//                        continue;
                      
                      Orbit &ob = X.modelspace->GetOrbit(b);
                      double nb = ob.occ;

                      for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                      {
                        for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                        {
                          double xab = GetMschemeMatrixElement_1b(X, a, ma, b, mb);
                          double yijbkla = 0;
                          if ( std::abs(xab) > 1.e-8 )
                          {
                            yijbkla = GetMschemeMatrixElement_3b(Y, i, mi, j, mj, b, mb, k, mk, l, ml, a, ma);
                          }
                          Zm_ijkl += (na - nb) * ( xab * yijbkla );
                        } // for mb
                      } // for ma
                    } // for b

                    for (auto b : Y.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}))
                    {
//                      if (a == b)
//                        continue;
                      Orbit &ob = Y.modelspace->GetOrbit(b);
                      double nb = ob.occ;
//                      if ( oa.j2 != ob.j2 )
//                        continue;
                      for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                      {
                        for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                        {                          
                          double xijbkla = 0;
                          double yab = GetMschemeMatrixElement_1b(Y, a, ma, b, mb);
                          if ( std::abs(yab) > 1.e-8 )
                          {
                            xijbkla = GetMschemeMatrixElement_3b(X, i, mi, j, mj, b, mb, k, mk, l, ml, a, ma);
                          }
                          Zm_ijkl -= (na - nb) * ( yab * xijbkla);
                        } // for mb
                      } // for ma
                    } // for b

                  } // for a

                  double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
                  double err = Zm_ijkl - ZJ_ijkl;
                  if (std::abs(err) > 1e-6)
                  {
                    std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                              << " {m} = " << mi << " " << mj << " " << mk << " " << ml
                              << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
                  }
                  summed_error += err * err;
                  sum_m += Zm_ijkl * Zm_ijkl;
                  sum_J += ZJ_ijkl * ZJ_ijkl;
                } // for ml
              } // for mk
            } // for mj
          } // for mi
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "\033[31mWARNING\033[0m " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ijkl = 1/6 * sum_abcd (n_a*n_b*n_c*nbar_d - nbar_a*nbar_b*nbar_c*n_d) * [  Xijdabc*Yabckld - Yijdabc*Xabckld  ]
//
bool UnitTest::Mscheme_Test_comm332_ppph_hhhpst(const Operator &X, const Operator &Y) 
{
  std::cout << __func__ << std::endl;
  int Lambda = Y.GetJRank();
  int Top = Y.GetJRank();
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.SetParticleRank(2);
  Z_J.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and Top == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }

  ReferenceImplementations::comm332_ppph_hhhpst(X, Y_copy, Z_J);
  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }
  std::cout << __func__ << "   Done with J-scheme commutator. Begin m-scheme version..." << std::endl;
  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = (parityX + parityY) % 2;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if ((oi.l + oj.l + ok.l + ol.l + parityY) % 2 > 0)
            continue; // check parity
          if ((oi.tz2 + oj.tz2 - ok.tz2 - ol.tz2) != TzZ * 2)
            continue; // check isospin projection

          for (int mi = -oi.j2; mi <= oi.j2; mi += 2) 
          {
            for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
            {
              for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
              {
                for (int ml = -ol.j2; ml <= ol.j2; ml += 2)
                {

                  if (  mi + mk < 0 or mk +ml < 0) // save some cpu time
                    continue;
                  
                  double Zm_ijkl = 0;
                  size_t norb = X.modelspace->GetNumberOrbits();
  #pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ijkl)
                  for (size_t a = 0; a < norb; a++)
                  {
                    Orbit &oa = X.modelspace->GetOrbit(a);
                    double na = oa.occ;
                    for (auto b : X.modelspace->all_orbits)
                    {
                      Orbit &ob = X.modelspace->GetOrbit(b);
                      double nb = ob.occ;

                      for (auto c : X.modelspace->all_orbits)
                      {
                        Orbit &oc = X.modelspace->GetOrbit(c);
                        double nc = oc.occ;

                        for (auto d : X.modelspace->all_orbits)
                        {
                          Orbit &od = X.modelspace->GetOrbit(d);
                          double nd = od.occ;
                          double occfactor = na * nb * nc * (1 - nd) - (1 - na) * (1 - nb) * (1 - nc) * nd;

                          // These make the contribution trivially zero, so I skip them in the name of efficiently
                          // testing the more complicated part. Commenting them out allows to check that the trivial stuff is right.
                          if ((oi.l + oj.l + od.l + oa.l + ob.l + oc.l) % 2 > 0)
                             continue;
                          // if ((oi.tz2 + oj.tz2 + od.tz2) != (oa.tz2 + ob.tz2 + oc.tz2))
                          //  continue;
                          if (std::abs(occfactor) < 1e-7)
                            continue;

                          for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                          {
                            for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                            {
                              for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                              {
                                for (int md = -od.j2; md <= od.j2; md += 2)
                                {
                                  double xijdabc = GetMschemeMatrixElement_3b(X, i, mi, j, mj, d, md, a, ma, b, mb, c, mc);
                                  double yijdabc = GetMschemeMatrixElement_3b(Y, i, mi, j, mj, d, md, a, ma, b, mb, c, mc);
                                  double xabckld = GetMschemeMatrixElement_3b(X, a, ma, b, mb, c, mc, k, mk, l, ml, d, md);
                                  double yabckld = GetMschemeMatrixElement_3b(Y, a, ma, b, mb, c, mc, k, mk, l, ml, d, md);
                                  Zm_ijkl += 1. / 6 * occfactor * (xijdabc * yabckld - yijdabc * xabckld);
                                } // for md
                              } // for mc
                            } // for mb
                          } // for ma
                        } // for d
                      } // for c
                    } // for b
                  } // for a

                  double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
                  double err = Zm_ijkl - ZJ_ijkl;
                  if (std::abs(err) > 1e-6)
                  {
                    std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                              << " {m} = " << mi << " " << mj << " " << mk << " " << ml
                              << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
                  }
                  summed_error += err * err;
                  sum_m += Zm_ijkl * Zm_ijkl;
                  sum_J += ZJ_ijkl * ZJ_ijkl;
                } // for ml
              } // for mk
            } // for mj
          } // for mi
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "\033[31mWARNING\033[0m " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

/// M-Scheme Formula:
//
// Z_ijkl = 1/4 (1-Pij)(1-Pkl) sum_abcd (n_a*n_b*nbar_c*nbar_d) * (  Xabicdk*Ycdjabl  )
//
//
bool UnitTest::Mscheme_Test_comm332_pphhst(const Operator &X, const Operator &Y) 
{
  std::cout << __func__ << std::endl;
  int Lambda = Y.GetJRank();
  int Top = Y.GetTRank();
  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and Top == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }
  ReferenceImplementations::comm332_pphhst(X, Y_copy, Z_J);
  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = (parityX + parityY) % 2;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < i)
          continue;
        for (auto l : X.modelspace->all_orbits)
        {
          if (l < k)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          if ((oi.l + oj.l + ok.l + ol.l + parityY) % 2 > 0)
            continue; // check parity
          if ((oi.tz2 + oj.tz2 - ok.tz2 - ol.tz2) != TzZ * 2)
            continue; // check isospin projection

          for (int mi = oi.j2; mi <= oi.j2; mi += 2)
          {
            for (int mj = -oj.j2; mj <= oj.j2; mj += 2)
            {
              for (int mk = -ok.j2; mk <= ok.j2; mk += 2)
              {
                for (int ml = -ol.j2; ml <= ol.j2; ml += 2)
                {
                  int twoJz = mi + mj - mk - ml;
                  if (std::abs(twoJz) > 2 * Lambda )
                    continue;                  
                  
                  double Zm_ijkl = 0;
                  size_t norb = X.modelspace->GetNumberOrbits();
                  //             for (auto a : X.modelspace->all_orbits )
  #pragma omp parallel for schedule(dynamic, 1) reduction(+ : Zm_ijkl)
                  for (size_t a = 0; a < norb; a++)
                  {
                    Orbit &oa = X.modelspace->GetOrbit(a);
                    double na = oa.occ;
                    for (auto b : X.modelspace->all_orbits)
                    {
                      Orbit &ob = X.modelspace->GetOrbit(b);
                      double nb = ob.occ;

                      for (auto c : X.modelspace->all_orbits)
                      {
                        Orbit &oc = X.modelspace->GetOrbit(c);
                        double nc = oc.occ;

                        for (auto d : X.modelspace->all_orbits)
                        {
                          Orbit &od = X.modelspace->GetOrbit(d);
                          double nd = od.occ;

                          //if ((oa.l + ob.l + oi.l + oc.l + od.l + ok.l) % 2 > 0 and (oa.l + ob.l + oi.l + oc.l + od.l + ol.l) % 2 > 0 and (oa.l + ob.l + oj.l + oc.l + od.l + ol.l) % 2 > 0 and (oa.l + ob.l + oj.l + oc.l + od.l + ok.l) % 2 > 0)
                          //  continue; // check parity for X and Y

                          //if ((oa.tz2 + ob.tz2 + oi.tz2) != (oc.tz2 + od.tz2 + ok.tz2) and (oa.tz2 + ob.tz2 + oi.tz2) != (oc.tz2 + od.tz2 + ol.tz2) and (oa.tz2 + ob.tz2 + oj.tz2) != (oc.tz2 + od.tz2 + ol.tz2) and (oa.tz2 + ob.tz2 + oj.tz2) != (oc.tz2 + od.tz2 + ok.tz2))
                          //  continue; // check isospin projection

                          for (int ma = -oa.j2; ma <= oa.j2; ma += 2)
                          {
                            for (int mb = -ob.j2; mb <= ob.j2; mb += 2)
                            {
                              for (int mc = -oc.j2; mc <= oc.j2; mc += 2)
                              {
                                for (int md = -od.j2; md <= od.j2; md += 2)
                                {

                                  double xabicdk = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, c, mc, d, md, k, mk);
                                  double xabicdl = GetMschemeMatrixElement_3b(X, a, ma, b, mb, i, mi, c, mc, d, md, l, ml);
                                  double xabjcdk = GetMschemeMatrixElement_3b(X, a, ma, b, mb, j, mj, c, mc, d, md, k, mk);
                                  double xabjcdl = GetMschemeMatrixElement_3b(X, a, ma, b, mb, j, mj, c, mc, d, md, l, ml);

                                  double ycdiabk = GetMschemeMatrixElement_3b(Y, c, mc, d, md, i, mi, a, ma, b, mb, k, mk);
                                  double ycdiabl = GetMschemeMatrixElement_3b(Y, c, mc, d, md, i, mi, a, ma, b, mb, l, ml);
                                  double ycdjabk = GetMschemeMatrixElement_3b(Y, c, mc, d, md, j, mj, a, ma, b, mb, k, mk);
                                  double ycdjabl = GetMschemeMatrixElement_3b(Y, c, mc, d, md, j, mj, a, ma, b, mb, l, ml);

                                  double occfactor = (1 - na) * (1 - nb) * nc * nd - na * nb * (1 - nc) * (1 - nd);
                                  double dz = 1. / 4 * occfactor * (xabicdl * ycdjabk - xabjcdl * ycdiabk - xabicdk * ycdjabl + xabjcdk * ycdiabl);
                                  Zm_ijkl += dz;

                                } // for md
                              } // for mc
                            } // for mb
                          } // for ma
                        } // for d
                      } // for c
                    } // for b
                  } // for a

                  double ZJ_ijkl = GetMschemeMatrixElement_2b(Z_J, i, mi, j, mj, k, mk, l, ml);
                  double err = Zm_ijkl - ZJ_ijkl;

                  if (std::abs(err) > 1e-6)
                  {
                    std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l = " << i << " " << j << " " << k << " " << l
                              << " {m} = " << mi << " " << mj << " " << mk << " " << ml
                              << "   Zm_ijkl = " << Zm_ijkl << "   ZJ_ijkl = " << ZJ_ijkl << "   err = " << err << std::endl;
                  }
                  summed_error += err * err;
                  sum_m += Zm_ijkl * Zm_ijkl;
                  sum_J += ZJ_ijkl * ZJ_ijkl;

                } // for ml
              } // for mk
            } // for mj
          } // for mi
        } // for l
      } // for k
    } // for j
  } // for i

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  if (Z_J.TwoBodyNorm() < 1e-6)
    std::cout << "\033[31mWARNING\033[0m " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
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
bool UnitTest::Mscheme_Test_comm233_pp_hhst(const Operator &X, const Operator &Y) // test not yet implemented
{
  std::cout << __func__ << std::endl;

  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = parityY;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;
  int JrankZ = Y.GetJRank();
  int Lambda = Y.GetJRank();

  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.SetParticleRank(3);
  Z_J.ThreeBody.SetMode("pn");
  Z_J.ThreeBody.Erase();


  Operator Y_copy(Y);
  if (Lambda == 0 and TzY == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }
  ReferenceImplementations::comm233_pp_hhst(X, Y_copy, Z_J);
  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();


  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;

  size_t norb = X.modelspace->GetNumberOrbits();
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < j)
          continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l < i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            if (m < l)
              continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n < m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l + parityY) % 2 != 0)
                continue;
              if ((oi.tz2 + oj.tz2 + ok.tz2 - ol.tz2 - om.tz2 - on.tz2) != 2 * TzY )
                continue;
              //              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i = oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;
                    if ( m_i + m_j + m_k < 0  )
                      continue;
                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ( std::abs(m_i + m_j + m_k - m_l - m_m - m_n) > 2 * Lambda )
                            continue;

                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          //                     if (not (m_i==1 and m_j==-1 and m_k==1 and m_l==1 and m_m==-1 and m_n==1) ) continue;
                          double z_ijklmn = 0;

                          //                     std::cout << " {m} " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n << std::endl;
                          //                     for ( auto a : X.modelspace->all_orbits )
//#pragma omp parallel for schedule(dynamic, 1) 
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (size_t b = 0; b < norb; b++)
                            {
                              //                        if (not (( a==4 and b==5) or (a==5 and b==4))  ) continue;
                              //                        if (not (( a==2 and b==2) or (a==2 and b==2))  ) continue;
                              Orbit &ob = X.modelspace->GetOrbit(b);
                              double occfactor = (1 - oa.occ)*(1 - ob.occ) - oa.occ * ob.occ;
                              //                       if (a!=11) continue;
                              for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                              {
                                for (int m_b = -ob.j2; m_b <= ob.j2; m_b += 2)
                                {

                                  ///       = 1/2 sum_ab (1-na-nb) [ Xijab*Yabklmn - Yijab*Xabklmn - Xkjab*Yabilmn + Ykjab*Xabilmn - Xikab*Yabjlmn + Yikab*Xabjlmn
                                  ///                              - Yijkabn*Xablm + Xijkabn*Yablm + Yijkabl*Xabnm - Xijkabl*Yabnm + Yijkabm*Xabln - Xijkabm*Yabln ]
                                  double xabilmn = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, i, m_i, l, m_l, m, m_m, n, m_n);
                                  double xabjlmn = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, j, m_j, l, m_l, m, m_m, n, m_n);
                                  double xabklmn = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, k, m_k, l, m_l, m, m_m, n, m_n);
                                  double xijkabl = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, l, m_l);
                                  double xijkabm = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, m, m_m);
                                  double xijkabn = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, n, m_n);

                                  double yabilmn = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, i, m_i, l, m_l, m, m_m, n, m_n);
                                  double yabjlmn = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, j, m_j, l, m_l, m, m_m, n, m_n);
                                  double yabklmn = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, k, m_k, l, m_l, m, m_m, n, m_n);
                                  double yijkabl = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, l, m_l);
                                  double yijkabm = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, m, m_m);
                                  double yijkabn = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, n, m_n);

                                  double xijab = GetMschemeMatrixElement_2b(X, i, m_i, j, m_j, a, m_a, b, m_b);
                                  double xkjab = GetMschemeMatrixElement_2b(X, k, m_k, j, m_j, a, m_a, b, m_b);
                                  double xikab = GetMschemeMatrixElement_2b(X, i, m_i, k, m_k, a, m_a, b, m_b);
                                  double xablm = GetMschemeMatrixElement_2b(X, a, m_a, b, m_b, l, m_l, m, m_m);
                                  double xabnm = GetMschemeMatrixElement_2b(X, a, m_a, b, m_b, n, m_n, m, m_m);
                                  double xabln = GetMschemeMatrixElement_2b(X, a, m_a, b, m_b, l, m_l, n, m_n);

                                  double yijab = GetMschemeMatrixElement_2b(Y, i, m_i, j, m_j, a, m_a, b, m_b);
                                  double ykjab = GetMschemeMatrixElement_2b(Y, k, m_k, j, m_j, a, m_a, b, m_b);
                                  double yikab = GetMschemeMatrixElement_2b(Y, i, m_i, k, m_k, a, m_a, b, m_b);
                                  double yablm = GetMschemeMatrixElement_2b(Y, a, m_a, b, m_b, l, m_l, m, m_m);
                                  double yabnm = GetMschemeMatrixElement_2b(Y, a, m_a, b, m_b, n, m_n, m, m_m);
                                  double yabln = GetMschemeMatrixElement_2b(Y, a, m_a, b, m_b, l, m_l, n, m_n);

                                  //                         z_ijklmn += 0.5*occfactor*( (xijab*yabklmn - yijab*xabklmn) - (xkjab*yabilmn - ykjab*xabilmn) - (xikab*yabjlmn - yikab*xabjlmn)
                                  //                                                   - (yijkabn*xablm - xijkabn*yablm) + (yijkabl*xabnm - xijkabl*yabnm) + (yijkabm*xabln + xijkabm*yabln) );
                                  double dz = 0.5 * occfactor * ((1 * xijab * yabklmn - 1 * yijab * xabklmn) 
                                                               - (1 * xkjab * yabilmn - 1 * ykjab * xabilmn) 
                                                               - (1 * xikab * yabjlmn - 1 * yikab * xabjlmn) 
                                                               - (1 * yijkabn * xablm - 1 * xijkabn * yablm) 
                                                               + (1 * yijkabl * xabnm - 1 * xijkabl * yabnm) 
                                                               + (1 * yijkabm * xabln - 1 * xijkabm * yabln));
                                  z_ijklmn += dz;
                                } // for m_b
                              } // for m_a
                            } // for b
                          } // for a

                          //                    double zJ_111 = Z_J.ThreeBody.GetME_pn(1,1,1,1,0,0,1,0,0);
                          //                    std::cout << " zJ_111 = " << zJ_111 << std::endl;
                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  if ( Z_J.ThreeBodyNorm() < 1e-6 ) 
    std::cout << "\033[31mWARNING\033[0m "  << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
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
bool UnitTest::Mscheme_Test_comm233_phst(const Operator &X, const Operator &Y) 
{
  std::cout << __func__ << std::endl;
  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = parityY;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;
  int JrankZ = Y.GetJRank();
  int Lambda = Y.GetJRank();

  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.SetParticleRank(3);
  Z_J.ThreeBody.SetMode("pn");
  Z_J.ThreeBody.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }
  ReferenceImplementations::comm233_phst(X, Y_copy, Z_J);
  if (Lambda == 0)
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j > i)
        continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k > j)
          continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l > i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            if (m > l)
              continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n > m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l + parityZ) % 2 != 0)
                continue;
              if ((oi.tz2 + oj.tz2 + ok.tz2 - ol.tz2 - om.tz2 - on.tz2) != 2 * TzZ)
                continue;

              //              if ( not (i==5 and j==4 and k==0 and l==5 and m==4 and n==0) ) continue;
              //              if ( not (i==3 and j==2 and k==0 and l==3 and m==2 and n==0) ) continue;
              //              if ( not (i==4 and j==1 and k==0 and l==4 and m==3 and n==2) ) continue;

              //              std::cout << " ijklmn " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
              // loop over projections
              for (int m_i = -oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if (std::abs(m_i + m_j + m_k - m_l - m_m - m_n) > Lambda * 2)
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          // if (not(m_i == 1 and m_j == -1 and m_k == 1 and m_l == 1 and m_m == -1 and m_n == 1))
                          //   continue;
                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (size_t b = 0; b < norb; b++)
                            {
                              Orbit &ob = X.modelspace->GetOrbit(b);
                              double occfactor = oa.occ - ob.occ;
                              for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                              {
                                for (int m_b = -ob.j2; m_b <= ob.j2; m_b += 2)
                                {

                                  //                         double xajkbmn = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, n,m_n );
                                  //                         double xaikbmn = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, n,m_n );
                                  //                         double xajibmn = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, n,m_n );
                                  //                         double xajkbln = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, l,m_l, n,m_n );
                                  //                         double xaikbln = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, l,m_l, n,m_n );
                                  //                         double xajibln = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, l,m_l, n,m_n );
                                  //                         double xajkbml = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, l,m_l );
                                  //                         double xaikbml = GetMschemeMatrixElement_3b( X, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, l,m_l );
                                  //                         double xajibml = GetMschemeMatrixElement_3b( X, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, l,m_l );

                                  double xbial = GetMschemeMatrixElement_2b(X, b, m_b, i, m_i, a, m_a, l, m_l);
                                  double xbjal = GetMschemeMatrixElement_2b(X, b, m_b, j, m_j, a, m_a, l, m_l);
                                  double xbkal = GetMschemeMatrixElement_2b(X, b, m_b, k, m_k, a, m_a, l, m_l);
                                  double xbiam = GetMschemeMatrixElement_2b(X, b, m_b, i, m_i, a, m_a, m, m_m);
                                  double xbjam = GetMschemeMatrixElement_2b(X, b, m_b, j, m_j, a, m_a, m, m_m);
                                  double xbkam = GetMschemeMatrixElement_2b(X, b, m_b, k, m_k, a, m_a, m, m_m);
                                  double xbian = GetMschemeMatrixElement_2b(X, b, m_b, i, m_i, a, m_a, n, m_n);
                                  double xbjan = GetMschemeMatrixElement_2b(X, b, m_b, j, m_j, a, m_a, n, m_n);
                                  double xbkan = GetMschemeMatrixElement_2b(X, b, m_b, k, m_k, a, m_a, n, m_n);

                                  //                         double yajkbmn = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, n,m_n );
                                  //                         double yaikbmn = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, n,m_n );
                                  //                         double yajibmn = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, n,m_n );
                                  //                         double yajkbln = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, l,m_l, n,m_n );
                                  //                         double yaikbln = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, l,m_l, n,m_n );
                                  //                         double yajibln = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, l,m_l, n,m_n );
                                  //                         double yajkbml = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, k,m_k, b,m_b, m,m_m, l,m_l );
                                  //                         double yaikbml = GetMschemeMatrixElement_3b( Y, a,m_a, i,m_i, k,m_k, b,m_b, m,m_m, l,m_l );
                                  //                         double yajibml = GetMschemeMatrixElement_3b( Y, a,m_a, j,m_j, i,m_i, b,m_b, m,m_m, l,m_l );

                                  double ybial = GetMschemeMatrixElement_2b(Y, b, m_b, i, m_i, a, m_a, l, m_l);
                                  double ybjal = GetMschemeMatrixElement_2b(Y, b, m_b, j, m_j, a, m_a, l, m_l);
                                  double ybkal = GetMschemeMatrixElement_2b(Y, b, m_b, k, m_k, a, m_a, l, m_l);
                                  double ybiam = GetMschemeMatrixElement_2b(Y, b, m_b, i, m_i, a, m_a, m, m_m);
                                  double ybjam = GetMschemeMatrixElement_2b(Y, b, m_b, j, m_j, a, m_a, m, m_m);
                                  double ybkam = GetMschemeMatrixElement_2b(Y, b, m_b, k, m_k, a, m_a, m, m_m);
                                  double ybian = GetMschemeMatrixElement_2b(Y, b, m_b, i, m_i, a, m_a, n, m_n);
                                  double ybjan = GetMschemeMatrixElement_2b(Y, b, m_b, j, m_j, a, m_a, n, m_n);
                                  double ybkan = GetMschemeMatrixElement_2b(Y, b, m_b, k, m_k, a, m_a, n, m_n);
                                  //////////////////----------------------------------------------------

                                  double xijalmb = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double yijalmb = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double xkjalmb = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double ykjalmb = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double xijanmb = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double yijanmb = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double xikalmb = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double yikalmb = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, a, m_a, l, m_l, m, m_m, b, m_b);
                                  double xijalnb = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double yijalnb = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double xkjanmb = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double ykjanmb = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double xikanmb = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double yikanmb = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, a, m_a, n, m_n, m, m_m, b, m_b);
                                  double xkjalnb = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double ykjalnb = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double xikalnb = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, a, m_a, l, m_l, n, m_n, b, m_b);
                                  double yikalnb = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, a, m_a, l, m_l, n, m_n, b, m_b);

                                  //  Zijkl =  - sum_ab (na-nb) [ Xbial Yajkbmn - Ybial Xajkbmn  -  Xbjal Yaikbmn + Ybjal Xaikbmn  -  Xbkal Yajibmn + Ybkal Xajibmn
                                  //                            - Xbiam Yajkbln + Ybiam Xajkbln  +  Xbjam Yaikbln - Ybjam Xaikbln  +  Xbkam Yajibln - Ybkam Xajibln
                                  //                            - Xbian Yajkbml + Ybian Xajkbml  +  Xbjan Yaikbml - Ybjan Xaikbml  +  Xbkan Yajibml - Ybkan Xajibml  ]
                                  // xijalmn = xaijblm = xajibml
                                  //                          z_ijklmn -= occfactor * (xbkan * yajibml - ybkan * xajibml);
                                  //                           z_ijklmn -= occfactor * (xbkan * yijalmb - ybkan*xijalmb);
                                  double dz = -occfactor * 1 * (xbkan * yijalmb - ybkan * xijalmb);
                                  dz += occfactor * 1 * (xbian * ykjalmb - ybian * xkjalmb);
                                  dz += occfactor * 1 * (xbjan * yikalmb - ybjan * xikalmb);
                                  dz += occfactor * 1 * (xbkal * yijanmb - ybkal * xijanmb);
                                  dz += occfactor * 1 * (xbkam * yijalnb - ybkam * xijalnb);
                                  dz -= occfactor * 1 * (xbial * ykjanmb - ybial * xkjanmb);
                                  dz -= occfactor * 1 * (xbjal * yikanmb - ybjal * xikanmb);
                                  dz -= occfactor * 1 * (xbiam * ykjalnb - ybiam * xkjalnb);
                                  dz -= occfactor * 1 * (xbjam * yikalnb - ybjam * xikalnb);
                                  z_ijklmn += dz;

                                } // for m_b
                              } // for m_a
                            } // for b
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-5)
                          {
                            std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  if ( Z_J.ThreeBodyNorm() < 1e-6 ) 
    std::cout << "\033[31mWARNING\033[0m " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}

// M-Scheme Formula:
//
// Z_ijkl = 1/6 sum_abc [na*nb*nc + (1-na)*(1-nb)(1-nc) ] * [  Xijkabc*Yabclmn - Yijkabc*Xabclmn  ]
//
//
bool UnitTest::Mscheme_Test_comm333_ppp_hhhst(const Operator &X, const Operator &Y) 
{
  std::cout << __func__ << std::endl;
  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = parityY;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;
  int JrankZ = Y.GetJRank();
  int Lambda = Y.GetJRank();
  Y.modelspace->PreCalculateSixJ();

  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.SetParticleRank(3);
  Z_J.ThreeBody.SetMode("pn");
  Z_J.ThreeBody.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and TzY == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }
  ReferenceImplementations::comm333_ppp_hhhst(X, Y_copy, Z_J);
  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      //      if (j<i) continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        //        if (k<j) continue;

        for (auto l : X.modelspace->all_orbits)
        {
          //          if (l<i) continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            //            if (m<l) continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              if (n < m)
                continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l + parityZ) % 2 != 0)
                continue;
              if (oi.tz2 + oj.tz2 + ok.tz2 - ol.tz2 - om.tz2 - on.tz2 != 2 * TzZ)
                continue;

              // loop over projections
              for (int m_i = oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ((m_i + m_j + m_k) != (m_l + m_m + m_n))
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;
                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (size_t b = 0; b < norb; b++)
                            {
                              Orbit &ob = X.modelspace->GetOrbit(b);
                              for (size_t c = 0; c < norb; c++)
                              {
                                Orbit &oc = X.modelspace->GetOrbit(c);
                                //                       double occfactor = oa.occ*ob.occ*oc.occ - (1-oa.occ)*(1-ob.occ)*(1-oc.occ); // why a minus sign??? That's wrong.
                                double occfactor = oa.occ * ob.occ * oc.occ + (1 - oa.occ) * (1 - ob.occ) * (1 - oc.occ); // Fixed by SRS July 19 2022
                                                                                                                          //                       if (a!=11) continue;
                                for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                                {
                                  for (int m_b = -ob.j2; m_b <= ob.j2; m_b += 2)
                                  {
                                    for (int m_c = -oc.j2; m_c <= oc.j2; m_c += 2)
                                    {
                                      double xijkabc = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, c, m_c);
                                      double xabclmn = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, c, m_c, l, m_l, m, m_m, n, m_n);
                                      double yijkabc = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, k, m_k, a, m_a, b, m_b, c, m_c);
                                      double yabclmn = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, c, m_c, l, m_l, m, m_m, n, m_n);

                                      z_ijklmn += 1. / 6 * occfactor * (xijkabc * yabclmn - yijkabc * xabclmn);
                                    } // for m_c
                                  } // for m_b
                                } // for m_a
                              } // for c
                            } // for b
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double err = z_ijklmn - ZJ_ijklmn;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ?  "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  if ( Z_J.ThreeBodyNorm() < 1e-4 ) 
    std::cout << "\033[31mWARNING\033[0m " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
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
bool UnitTest::Mscheme_Test_comm333_pph_hhpst(const Operator &X, const Operator &Y) 
{
  std::cout << __func__ << std::endl;
  int parityX = X.GetParity();
  int parityY = Y.GetParity();
  int parityZ = parityY;
  int TzX = X.GetTRank();
  int TzY = Y.GetTRank();
  int TzZ = TzX + TzY;
  int JrankZ = Y.GetJRank();
  int Lambda = Y.GetJRank();
  Y.modelspace->PreCalculateSixJ();

  Operator Z_J(Y);
  Z_J.SetHermitian();
  Z_J.Erase();
  Z_J.SetParticleRank(3);
  Z_J.ThreeBody.SetMode("pn");
  Z_J.ThreeBody.Erase();

  Operator Y_copy(Y);
  if (Lambda == 0 and TzY == 0)
  {
    Y_copy.MakeReduced();
    Z_J.MakeReduced();
  }
  ReferenceImplementations::comm333_pph_hhpst(X, Y_copy, Z_J);
  if (Lambda == 0 and !Y.IsReduced())
  {
    // Y_copy.MakeNotReduced();
    Z_J.MakeNotReduced();
  }

  if (Z_J.IsHermitian())
    Z_J.Symmetrize();
  else if (Z_J.IsAntiHermitian())
    Z_J.AntiSymmetrize();

  double summed_error = 0;
  double sum_m = 0;
  double sum_J = 0;

  std::cout << __func__ << "    Start M loops" << std::endl;
  for (auto i : X.modelspace->all_orbits)
  {
    Orbit &oi = X.modelspace->GetOrbit(i);
    //    int mi = oi.j2;
    for (auto j : X.modelspace->all_orbits)
    {
      if (j < i)
        continue;
      //      if (j>i) continue;
      Orbit &oj = X.modelspace->GetOrbit(j);
      for (auto k : X.modelspace->all_orbits)
      {
        Orbit &ok = X.modelspace->GetOrbit(k);
        if (k < j)
          continue;
        //        if (k>j) continue;

        for (auto l : X.modelspace->all_orbits)
        {
          if (l > i)
            continue;
          Orbit &ol = X.modelspace->GetOrbit(l);
          for (auto m : X.modelspace->all_orbits)
          {
            //            if (m<l) continue;
            //            if (m>l) continue;
            Orbit &om = X.modelspace->GetOrbit(m);
            for (auto n : X.modelspace->all_orbits)
            {
              //              if (n<m) continue;
              //              if (n>m) continue;
              Orbit &on = X.modelspace->GetOrbit(n);
              if ((oi.l + oj.l + ok.l + ol.l + om.l + on.l + parityZ) % 2 != 0)
                continue;
              // if ((oi.tz2 + oj.tz2 + ok.tz2 - ol.tz2 - om.tz2 - on.tz2) != 2 * TzZ)
              //   continue;

              // loop over projections
              for (int m_i = oi.j2; m_i <= oi.j2; m_i += 2)
              {
                for (int m_j = -oj.j2; m_j <= oj.j2; m_j += 2)
                {
                  for (int m_k = -ok.j2; m_k <= ok.j2; m_k += 2)
                  {
                    if ((i == j and m_i == m_j) or (i == k and m_i == m_k) or (j == k and m_j == m_k))
                      continue;

                    for (int m_l = -ol.j2; m_l <= ol.j2; m_l += 2)
                    {
                      for (int m_m = -om.j2; m_m <= om.j2; m_m += 2)
                      {
                        for (int m_n = -on.j2; m_n <= on.j2; m_n += 2)
                        {
                          if ( std::abs(m_i + m_j + m_k - m_l - m_m - m_n) > 2 * Lambda)
                            continue;
                          if ((l == m and m_l == m_m) or (l == n and m_l == m_n) or (m == n and m_m == m_n))
                            continue;

                          // if (not(m_i == 1 and m_j == -1 and m_k == 1 and m_l == 1 and m_m == -1 and m_n == 1))
                          //   continue;

                          if ( (m_i + m_j + m_k < 0 and m_l + m_m + m_n < 0))
                            continue;

                          double z_ijklmn = 0;

                          size_t norb = X.modelspace->GetNumberOrbits();
                          //                     for ( auto a : X.modelspace->all_orbits )
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : z_ijklmn)
                          for (size_t a = 0; a < norb; a++)
                          {
                            Orbit &oa = X.modelspace->GetOrbit(a);
                            for (size_t b = 0; b < norb; b++)
                            {
                              Orbit &ob = X.modelspace->GetOrbit(b);
                              for (size_t c = 0; c < norb; c++)
                              {
                                Orbit &oc = X.modelspace->GetOrbit(c);
                                double occfactor = oa.occ * ob.occ * (1 - oc.occ) + (1 - oa.occ) * (1 - ob.occ) * oc.occ;
                                for (int m_a = -oa.j2; m_a <= oa.j2; m_a += 2)
                                {
                                  for (int m_b = -ob.j2; m_b <= ob.j2; m_b += 2)
                                  {
                                    for (int m_c = -oc.j2; m_c <= oc.j2; m_c += 2)
                                    {

                                      double xabklmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, k, m_k, l, m_l, m, m_m, c, m_c); // direct
                                      double xijcabn = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, n, m_n); // direct

                                      double xabilmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, i, m_i, l, m_l, m, m_m, c, m_c); // Pik
                                      double xkjcabn = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, n, m_n); // Pik
                                      double xabjlmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, j, m_j, l, m_l, m, m_m, c, m_c); // Pjk
                                      double xikcabn = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, n, m_n); // Pjk
                                      double xabknmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, k, m_k, n, m_n, m, m_m, c, m_c); // Pln
                                      double xijcabl = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, l, m_l); // Pln
                                      double xabklnc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, k, m_k, l, m_l, n, m_n, c, m_c); // Pmn
                                      double xijcabm = GetMschemeMatrixElement_3b(X, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, m, m_m); // Pmn

                                      double xabinmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, i, m_i, n, m_n, m, m_m, c, m_c); // Pik Pln
                                      double xkjcabl = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, l, m_l); // Pik Pln
                                      double xabjnmc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, j, m_j, n, m_n, m, m_m, c, m_c); // Pjk Pln
                                      double xikcabl = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, l, m_l); // Pjk Pln
                                      double xabilnc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, i, m_i, l, m_l, n, m_n, c, m_c); // Pik Pmn
                                      double xkjcabm = GetMschemeMatrixElement_3b(X, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, m, m_m); // Pik Pmn
                                      double xabjlnc = GetMschemeMatrixElement_3b(X, a, m_a, b, m_b, j, m_j, l, m_l, n, m_n, c, m_c); // Pjk Pmn
                                      double xikcabm = GetMschemeMatrixElement_3b(X, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, m, m_m); // Pjk Pmn

                                      double yabklmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, k, m_k, l, m_l, m, m_m, c, m_c); // direct
                                      double yijcabn = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, n, m_n); // direct

                                      double yabilmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, i, m_i, l, m_l, m, m_m, c, m_c); // Pik
                                      double ykjcabn = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, n, m_n); // Pik
                                      double yabjlmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, j, m_j, l, m_l, m, m_m, c, m_c); // Pjk
                                      double yikcabn = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, n, m_n); // Pjk
                                      double yabknmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, k, m_k, n, m_n, m, m_m, c, m_c); // Pln
                                      double yijcabl = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, l, m_l); // Pln
                                      double yabklnc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, k, m_k, l, m_l, n, m_n, c, m_c); // Pmn
                                      double yijcabm = GetMschemeMatrixElement_3b(Y, i, m_i, j, m_j, c, m_c, a, m_a, b, m_b, m, m_m); // Pmn

                                      double yabinmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, i, m_i, n, m_n, m, m_m, c, m_c); // Pik Pln
                                      double ykjcabl = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, l, m_l); // Pik Pln
                                      double yabjnmc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, j, m_j, n, m_n, m, m_m, c, m_c); // Pjk Pln
                                      double yikcabl = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, l, m_l); // Pjk Pln
                                      double yabilnc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, i, m_i, l, m_l, n, m_n, c, m_c); // Pik Pmn
                                      double ykjcabm = GetMschemeMatrixElement_3b(Y, k, m_k, j, m_j, c, m_c, a, m_a, b, m_b, m, m_m); // Pik Pmn
                                      double yabjlnc = GetMschemeMatrixElement_3b(Y, a, m_a, b, m_b, j, m_j, l, m_l, n, m_n, c, m_c); // Pjk Pmn
                                      double yikcabm = GetMschemeMatrixElement_3b(Y, i, m_i, k, m_k, c, m_c, a, m_a, b, m_b, m, m_m); // Pjk Pmn

                                      //        =   1/2 sum_abc [na*nb*(1-nc) + (1-na)*(1-nb)*nc ] *
                                      //                [ (Xabkcmn*Ycijabl - Yabkcmn*Xcijabl) - (Xabicmn*Yckjabl - Yabicmn*Xckjabl)  - (Xabjcmn*Ycikabl - Yabjcmn*Xcikabl)
                                      //                 -(Xabkcln*Ycijabm - Yabkcln*Xcijabm) + (Xabicln*Yckjabm - Yabicln*Xckjabm)  + (Xabjcln*Ycikabm - Yabjcln*Xcikabm)
                                      //                 -(Xabkcml*Ycijabn - Yabkcml*Xcijabn) + (Xabicml*Yckjabn - Yabicml*Xckjabn)  + (Xabjcml*Ycikabn - Yabjcml*Xcikabn) ]
                                      double dz = 0;
                                      dz += xabklmc * yijcabn - yabklmc * xijcabn; // Z1
                                      dz -= xabilmc * ykjcabn - yabilmc * xkjcabn; // Z2
                                      dz -= xabjlmc * yikcabn - yabjlmc * xikcabn; // Z3
                                      dz -= xabknmc * yijcabl - yabknmc * xijcabl; // Z4
                                      dz -= xabklnc * yijcabm - yabklnc * xijcabm; // Z5
                                      dz += xabinmc * ykjcabl - yabinmc * xkjcabl; // Z6
                                      dz += xabjnmc * yikcabl - yabjnmc * xikcabl; // Z7
                                      dz += xabilnc * ykjcabm - yabilnc * xkjcabm; // Z8
                                      dz += xabjlnc * yikcabm - yabjlnc * xikcabm; // Z9
                                      z_ijklmn += 0.5 * occfactor * dz;
                                    } // for m_c
                                  } // for m_b
                                } // for m_a
                              } // for c
                            } // for b
                          } // for a

                          double ZJ_ijklmn = GetMschemeMatrixElement_3b(Z_J, i, m_i, j, m_j, k, m_k, l, m_l, m, m_m, n, m_n);
                          double ZJ_hermconj = GetMschemeMatrixElement_3b(Z_J, l, m_l, m, m_m, n, m_n, i, m_i, j, m_j, k, m_k);
                          double err = z_ijklmn - ZJ_ijklmn;
                          //                    std::cout << z_ijklmn << std::endl;
                          if (std::abs(err) > 1e-6)
                          {
                            std::cout << "\033[31mTrouble\033[0m in " << __func__ << "  i,j,k,l,m,n = " << i << " " << j << " " << k << " " << l << " " << m << " " << n
                                      << " {m} = " << m_i << " " << m_j << " " << m_k << " " << m_l << " " << m_m << " " << m_n
                                      << "   Zm_ijklmn = " << z_ijklmn << "   ZJ_ijklmn = " << ZJ_ijklmn << "   err = " << err << "  hc = " << ZJ_hermconj << std::endl;
                          }
                          summed_error += err * err;
                          sum_m += z_ijklmn * z_ijklmn;
                          sum_J += ZJ_ijklmn * ZJ_ijklmn;
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

  bool passed = std::abs(summed_error) < 1e-6;
  std::string passfail = passed ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m";
  if ( Z_J.ThreeBodyNorm() < 1e-6 ) 
    std::cout << "\033[31mWARNING\033[0m " << __func__ << "||Z_J 2b|| = 0. Trivial test?" << std::endl;
  std::cout << "   " << __func__ << "  sum_m, sum_J = " << sum_m << " " << sum_J
            << "    summed error = " << summed_error << "  => " << passfail << std::endl;
  return passed;
}













