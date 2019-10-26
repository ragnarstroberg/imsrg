
#include "Generator.hh"
#include "Commutator.hh"
#include "Operator.hh"
#include "imsrg_util.hh" // for VectorUnion
#include "AngMom.hh" // for TalmiB and Moshinsky
#include "PhysicalConstants.hh" // for HBARC and M_NUCLEON

#include "omp.h"
#include <string>
#include <gsl/gsl_sf_gamma.h> // for gsl_sf_gamma_inc_P  used in the RspaceRegulator

//using namespace imsrg_util;

Generator::Generator()
  : generator_type("white"), denominator_cutoff(1e-6)  , denominator_delta(0), denominator_delta_index(-1)
{}


void Generator::Update(Operator * H_s, Operator * Eta_s)
{
   Eta_s->EraseOneBody();
   Eta_s->EraseTwoBody();
   AddToEta(H_s,Eta_s);
}


void Generator::AddToEta(Operator * H_s, Operator * Eta_s)
{
   double start_time = omp_get_wtime();
   H = H_s;
   Eta = Eta_s;
   modelspace = H->GetModelSpace();

        if (generator_type == "wegner")                       ConstructGenerator_Wegner(); // never tested, probably doesn't work.
   else if (generator_type == "shell-model-wegner")           ConstructGenerator_ShellModel_Wegner(); // never tested, probably doesn't work.
   else if (generator_type == "white")                        ConstructGenerator_White();
   else if (generator_type == "atan")                         ConstructGenerator_Atan();
   else if (generator_type == "imaginary-time")               ConstructGenerator_ImaginaryTime();
   else if (generator_type == "qtransfer-atan")               ConstructGenerator_QTransferAtan(1);
   else if (generator_type == "shell-model")                  ConstructGenerator_ShellModel();
   else if (generator_type == "shell-model-atan")             ConstructGenerator_ShellModel_Atan();
   else if (generator_type == "shell-model-atan-npnh")        ConstructGenerator_ShellModel_Atan_NpNh();
   else if (generator_type == "shell-model-imaginary-time")   ConstructGenerator_ShellModel_ImaginaryTime();
   else if (generator_type == "hartree-fock")                 ConstructGenerator_HartreeFock();
   else if (generator_type == "1PA")                          ConstructGenerator_1PA();
   else if (generator_type.find("qtransfer-atan") != std::string::npos )
   {
     int n;
     std::istringstream( generator_type.substr( generator_type.find("_")+1) ) >> n;
     ConstructGenerator_QTransferAtan(n);
   }
   else if (generator_type == "rspace")                       ConstructGenerator_Rspace();
   else
   {
      std::cout << "Error. Unkown generator_type: " << generator_type << std::endl;
   }
   Eta->profiler.timer["UpdateEta"] += omp_get_wtime() - start_time;

}


// Old method used to test some things out. Not typically used.
void Generator::SetDenominatorDeltaOrbit(std::string orb)
{
  if (orb == "all")
     SetDenominatorDeltaIndex(-12345);
  else
  {
     SetDenominatorDeltaIndex( modelspace->GetOrbitIndex(orb) );
     std::cout << "Setting denominator delta orbit " << orb << " => " << modelspace->GetOrbitIndex(orb) << std::endl;
  }
}


// Epstein-Nesbet energy denominators for White-type generator_types
double Generator::Get1bDenominator(int i, int j) 
{
   double ni = modelspace->GetOrbit(i).occ;
   double nj = modelspace->GetOrbit(j).occ;
   
   double denominator = H->OneBody(i,i) - H->OneBody(j,j);
   denominator += ( ni-nj ) * H->TwoBody.GetTBMEmonopole(i,j,i,j);

   if (denominator_delta_index==-12345 or i == denominator_delta_index or j==denominator_delta_index)
     denominator += denominator_delta;

   if (std::abs(denominator)<denominator_cutoff)
     denominator = denominator_cutoff;
//     denominator *= denominator_cutoff/(std::abs(denominator)+1e-6);

   return denominator;
}


double Generator::Get2bDenominator(int ch, int ibra, int iket) 
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int i = bra.p;
   int j = bra.q;
   int k = ket.p;
   int l = ket.q;
   double denominator = H->OneBody(i,i)+ H->OneBody(j,j) - H->OneBody(k,k) - H->OneBody(l,l);
   if (denominator_delta_index == -12345) denominator += denominator_delta;
   double ni = bra.op->occ;
   double nj = bra.oq->occ;
   double nk = ket.op->occ;
   double nl = ket.oq->occ;

   denominator       += ( 1-ni-nj ) * H->TwoBody.GetTBMEmonopole(i,j,i,j); // pp'pp'
   denominator       -= ( 1-nk-nl ) * H->TwoBody.GetTBMEmonopole(k,l,k,l); // hh'hh'
   denominator       += ( ni-nk ) * H->TwoBody.GetTBMEmonopole(i,k,i,k); // phph
   denominator       += ( ni-nl ) * H->TwoBody.GetTBMEmonopole(i,l,i,l); // ph'ph'
   denominator       += ( nj-nk ) * H->TwoBody.GetTBMEmonopole(j,k,j,k); // p'hp'h
   denominator       += ( nj-nl ) * H->TwoBody.GetTBMEmonopole(j,l,j,l); // p'h'p'h'

   if (std::abs(denominator)<denominator_cutoff)
     denominator = denominator_cutoff;
//     denominator *= denominator_cutoff/(std::abs(denominator)+1e-6);
   return denominator;
}

double Generator::Get3bDenominator( int i, int j, int k, int l, int m, int n )
{
  auto h1 = H->OneBody;
  double denominator = h1(i,i) + h1(j,j) + h1(k,k) - h1(l,l) - h1(m,m) - h1(n,n);
  return denominator;
}

// Keep the Jdependence for the Gamma_ijij and Gamma_klkl terms, because it's
// relatively unambiguous to work out
double Generator::Get2bDenominator_Jdep(int ch, int ibra, int iket) 
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int i = bra.p;
   int j = bra.q;
   int k = ket.p;
   int l = ket.q;
   double denominator = H->OneBody(i,i)+ H->OneBody(j,j) - H->OneBody(k,k) - H->OneBody(l,l);
   if (denominator_delta_index == -12345) denominator += denominator_delta;
   double ni = bra.op->occ;
   double nj = bra.oq->occ;
   double nk = ket.op->occ;
   double nl = ket.oq->occ;

   denominator       += ( 1-ni-nj ) * H->TwoBody.GetTBME(tbc.J,i,j,i,j); // pp'pp'
   denominator       -= ( 1-nk-nl ) * H->TwoBody.GetTBME(tbc.J,k,l,k,l); // hh'hh'
   denominator       += ( ni-nk ) * H->TwoBody.GetTBMEmonopole(i,k,i,k); // phph
   denominator       += ( ni-nl ) * H->TwoBody.GetTBMEmonopole(i,l,i,l); // ph'ph'
   denominator       += ( nj-nk ) * H->TwoBody.GetTBMEmonopole(j,k,j,k); // p'hp'h
   denominator       += ( nj-nl ) * H->TwoBody.GetTBMEmonopole(j,l,j,l); // p'h'p'h'

   if (std::abs(denominator)<denominator_cutoff)
     denominator = denominator_cutoff;
//     denominator *= denominator_cutoff/(std::abs(denominator)+1e-6);
   return denominator;
}


//// I haven't used this, so I don't know if it's right.
//void Generator::ConstructGenerator_Wegner()
//{
//   Operator H_diag = *H;
//   H_diag.ZeroBody = 0;
//   for (auto& a : modelspace->holes)
//   {
////      index_t a = it_a.first;
//      for (auto& b : modelspace->valence)
//      {
//         H_diag.OneBody(a,b) =0;
//         H_diag.OneBody(b,a) =0;
//      }
//   }
//
//   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
//   {  // Note, should also decouple the v and q spaces
//      // This is wrong. The projection operator should be different.
//      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
//      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_pp(), tbc.GetKetIndex_ph() ).zeros();
//      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_hh(), tbc.GetKetIndex_ph() ).zeros();
//      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_ph(), tbc.GetKetIndex_pp() ).zeros();
//      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_ph(), tbc.GetKetIndex_hh() ).zeros();
//      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_pp(), tbc.GetKetIndex_hh() ).zeros();
//      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_hh(), tbc.GetKetIndex_pp() ).zeros();
//   }
//   Eta->SetToCommutator(H_diag,*H);
////   *Eta = Commutator(H_diag,*H);
//}




void Generator::ConstructGenerator_Wegner()
{
   Operator H_od = 0 * (*H);
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         H_od.OneBody(i,a) = H->OneBody(i,a);
         H_od.OneBody(a,i) = H_od.OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<H_od.nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& HOD2 =  H_od.TwoBody.GetMatrix(ch);
      arma::mat& H2   =    H->TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            HOD2(ibra,iket) = H2(ibra,iket) ;
            HOD2(iket,ibra) = HOD2(ibra,iket) ; 
         }
      }
    }
   *Eta = Commutator::Commutator(*H,H_od);
}








void Generator::ConstructGenerator_White()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = H->OneBody(i,a)/denominator;
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}






void Generator::ConstructGenerator_Atan()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
//            double denominator = Get2bDenominator_Jdep(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
    // Three body piece
    int E3cut = 4;
    double t_start = omp_get_wtime();
//    std::cout << "checking ranks in generator. " << Eta->GetParticleRank() << "  and " << H->GetParticleRank() << "   Norms: " << Eta->ThreeBodyNorm() << "  " << H->ThreeBodyNorm()  << std::endl;
    if ( Eta->GetParticleRank()>2 and H->GetParticleRank()>2 )
    {
     std::cout << "looping in generator 3-body part " << std::endl;
    for (auto a : modelspace->core )
    {
     Orbit& oa = modelspace->GetOrbit(a);
     for (auto b : modelspace->core )
     {
      Orbit& ob = modelspace->GetOrbit(b);
      int Jab_min = std::abs(oa.j2-ob.j2)/2;
      int Jab_max = (oa.j2+ob.j2)/2;
      for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
      {
       for (auto c : modelspace->core )
       {
        Orbit& oc = modelspace->GetOrbit(c);
        if ( (2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l) > E3cut ) continue;

        for ( auto i : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
        {
         Orbit& oi = modelspace->GetOrbit(i);
         for ( auto j : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
         {
          Orbit& oj = modelspace->GetOrbit(j);
          int Jij_min = std::abs(oi.j2-oj.j2)/2;
          int Jij_max = (oi.j2+oj.j2)/2;
          for (int Jij=Jij_min; Jij<=Jij_max; Jij++)
          {
           for ( auto k : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
           {
            Orbit& ok = modelspace->GetOrbit(k);
            if ( (2*(oi.n+oj.n+ok.n)+oi.l+oj.l+ok.l) > 2*E3cut ) continue;
            if ( (oa.l+ob.l+oc.l+oi.l+oj.l+ok.l)%2>0 ) continue;
            if ( (oa.tz2+ob.tz2+oc.tz2) != (oi.tz2+oj.tz2+ok.tz2) ) continue;
            double denominator = Get3bDenominator( a,b,c, i,j,k ) ;

            int twoJ_min = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jij-ok.j2) );
            int twoJ_max = std::min( 2*Jab+oc.j2, 2*Jij+ok.j2 );
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ++)
            {
              double ME_od = H->ThreeBody.GetME_pn( Jab, Jij, twoJ, a,b,c,i,j,k);
              double eta = 0.5*atan(2*ME_od / denominator);
//              double eta = (ME_od / denominator);
//              if (std::abs(eta)>1e-6)
//              {
//                std::cout << "Nonzero eta !  <" << a << " " << b << " " << c << ", " << Jab << " | " << i << " " << j << " " << k << ", " << Jij << " > _ " << twoJ << "  => " << ME_od << " / " << denominator << " = " << eta << std::endl;
//              }
              Eta->ThreeBody.AddToME_pn( Jab, Jij, twoJ, a,b,c,i,j,k,  eta);
            }
           }
          }
         }
        }
       }
      } 
     } 
    }
    }
    H->profiler.timer["Update Eta 3body"] += omp_get_wtime() - t_start;
}



/// Imaginary time generator \f[ \eta = sgn(\Delta) h_{od} \]
void Generator::ConstructGenerator_ImaginaryTime()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
         if (denominator==0) denominator = 1;
         Eta->OneBody(i,a) += H->OneBody(i,a) *denominator/std::abs(denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            if (denominator==0) denominator = 1;
            ETA2(ibra,iket) += H2(ibra,iket) * denominator / std::abs(denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}



void Generator::ConstructGenerator_QTransferAtan(int n)
{
   using PhysConst::M_NUCLEON;
   using PhysConst::HBARC;
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = pow(std::abs(denominator)*M_NUCLEON/HBARC/HBARC, 0.5*n ) * 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
//            double denominator = Get2bDenominator_Jdep(ch,ibra,iket);
            ETA2(ibra,iket) = pow(std::abs(denominator)*M_NUCLEON/HBARC/HBARC, 0.5*n) * 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}



void Generator::ConstructGenerator_ShellModel()
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : imsrg_util::VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : imsrg_util::VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = H->OneBody(i,a)/denominator;
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }


   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H->TwoBody.GetMatrix(ch);

      // Decouple the core
      for ( auto& iket : imsrg_util::VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}




void Generator::ConstructGenerator_ShellModel_Atan()
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : imsrg_util::VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : imsrg_util::VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }


   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H->TwoBody.GetMatrix(ch);

      // Decouple the core
      for ( auto& iket : imsrg_util::VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator) ;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}




void Generator::ConstructGenerator_ShellModel_Wegner()
{
   Operator H_od = 0 * (*H);
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : imsrg_util::VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : imsrg_util::VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         H_od.OneBody(i,a) = H->OneBody(i,a);
         H_od.OneBody(a,i) = H_od.OneBody(i,a);
      }
   }


   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& HOD2 =  H_od.TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H->TwoBody.GetMatrix(ch);

      // Decouple the core
      for ( auto& iket : imsrg_util::VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            HOD2(ibra,iket) = H2(ibra,iket) ;
            HOD2(iket,ibra) = HOD2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
         {
            HOD2(ibra,iket) = H2(ibra,iket) ;
            HOD2(iket,ibra) = HOD2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
    *Eta = Commutator::Commutator(*H,H_od);
}




/// Imaginary time generator for a valence space
void Generator::ConstructGenerator_ShellModel_ImaginaryTime()
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : imsrg_util::VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : imsrg_util::VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
         if (denominator==0) denominator = 1;
         Eta->OneBody(i,a) += H->OneBody(i,a) *denominator/std::abs(denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }


   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H->TwoBody.GetMatrix(ch);

      // Decouple the core
      for ( auto& iket : imsrg_util::VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            if (denominator==0) denominator = 1;
            ETA2(ibra,iket) = H2(ibra,iket) * denominator / std::abs(denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            if (denominator==0) denominator = 1;
            ETA2(ibra,iket) = H2(ibra,iket) * denominator / std::abs(denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}



void Generator::ConstructGenerator_ShellModel_Atan_NpNh()
{
  ConstructGenerator_ShellModel_Atan();

//  std::cout << "In ShellModel_Atat_NpNh, adding to Eta" << std::endl;
  // decouple f_cc'
  for ( auto& c : modelspace->core )
  {
   for ( auto& cprime : modelspace->core )
   {
     if (cprime<=c) continue;
     double denominator = Get1bDenominator(c,cprime);
     Eta->OneBody(c,cprime) = 0.5*atan(2*H->OneBody(c,cprime) / denominator );
//     std::cout << "c,cprime = " << c << " " << cprime << "  etacc' = " << Eta->OneBody(c,cprime) << std::endl;
     Eta->OneBody(cprime,c) = - Eta->OneBody(c,cprime);
   }
  }

  int nchan = modelspace->GetNumberTwoBodyChannels();
  for (int ch=0;ch<nchan;++ch)
  {
     TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
//     std::cout << "ch = " << ch << "  vc size = " << tbc.GetKetIndex_vc().size() << "   qc size = " << tbc.GetKetIndex_qc().size() << std::endl;
     arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
     arma::mat& H2 =  H->TwoBody.GetMatrix(ch);
  // decouple Gamma_qcvc'
     for (auto& iket : tbc.GetKetIndex_vc())
     {
       Ket& ket = tbc.GetKet(iket);
       for (auto& ibra : tbc.GetKetIndex_qc())
       {
         Ket& bra = tbc.GetKet(ibra);
//         std::cout << bra.p << " " << bra.q << " " << ket.p << " " << ket.q << std::endl;
         if ((ket.p==bra.p) or (ket.p==bra.q) or (ket.q==bra.p) or (ket.q==bra.q) ) continue;
         double denominator = Get2bDenominator(ch,ibra,iket);
         ETA2(ibra,iket) = 0.5*atan( 2*H2(ibra,iket) / denominator );
         ETA2(iket,ibra) = -ETA2(ibra,iket);
//         std::cout << "   qcvc': " << ket.p << " " << ket.q << " " << bra.p << " " << bra.q << "    " << ETA2(ibra,iket) << std::endl;
       }
     }
  // decouple Gamma_pcc'c''
     for (auto& iket : tbc.GetKetIndex_cc())
     {
       Ket& ket = tbc.GetKet(iket);
       for (auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_vc(), tbc.GetKetIndex_qc() )  )
       {
         Ket& bra = tbc.GetKet(ibra);
         if ((ket.p==bra.p) or (ket.p==bra.q) or (ket.q==bra.p) or (ket.q==bra.q) ) continue;
         double denominator = Get2bDenominator(ch,ibra,iket);
         ETA2(ibra,iket) = 0.5*atan( 2*H2(ibra,iket) / denominator );
         ETA2(iket,ibra) = -ETA2(ibra,iket);
       }
     }
  }
//  std::cout << "all done" << std::endl;
}


void Generator::ConstructGenerator_HartreeFock()
{
   Eta->SetParticleRank(1);
   // One body piece -- eliminate ph bits
//   unsigned int norbits = modelspace->GetNumberOrbits();
//   for (unsigned int i=0;i<norbits;++i)
   for (auto i : modelspace->all_orbits)
   {
//      for (unsigned int j=0; j<i; ++j)
      for (auto j : modelspace->all_orbits)
      {
         if (j>i) continue;
         double denominator = Get1bDenominator(i,j);
         Eta->OneBody(i,j) = H->OneBody(i,j)/denominator;
         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
      }
   } 
}




// So far this is useless
void Generator::ConstructGenerator_1PA()
{

   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : imsrg_util::VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : imsrg_util::VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }


   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H->TwoBody.GetMatrix(ch);

      // Decouple the core
      for ( auto& iket : imsrg_util::VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }
    }

}









// Weight each matrix element in eta by the matrix element of the regulator
// which in the current implementation is a step function in R
void Generator::ConstructGenerator_Rspace()
{

   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) =  RspaceRegulator.OneBody(i,a) * 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      auto& R2 = RspaceRegulator.TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
         for ( auto& ibra : imsrg_util::VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
//            double denominator = Get2bDenominator_Jdep(ch,ibra,iket);
            ETA2(ibra,iket) =  R2(ibra,iket)  *  0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }

}








void Generator::SetRegulatorLength(double r0)
{
  using PhysConst::M_NUCLEON;
  using PhysConst::HBARC;
//  std::cout << "setting Regulator length to " << r0 << std::endl;
  regulator_length = r0;

  RspaceRegulator = Operator(*modelspace);
//  std::cout << "Succeeded in creating RspaceRegulator" << std::endl;
  double oscillator_b_squared = HBARC*HBARC / M_NUCLEON / modelspace->GetHbarOmega();
  double x0 = r0*r0/oscillator_b_squared; // express length in units of oscillator length: r0^2 = x0 * b^2
  
  // First, do the one-body part
//  int norb = modelspace->GetNumberOrbits();
//  for (int a=0; a<norb; a++)
  for (auto a : modelspace->all_orbits )
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for ( int b : RspaceRegulator.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
    {
      Orbit& ob = modelspace->GetOrbit(b);
      double matel = 0;
      for (int p=oa.l; p<= oa.n+ob.n+oa.l; p++)
      {
        double I_p = gsl_sf_gamma_inc_P(p+1.5,x0);  // Talmi integral for a step function Theta(r0 -r );
        double B_p = AngMom::TalmiB(oa.n, oa.l, ob.n, ob.l, p);
        matel += I_p * B_p;
      }
      RspaceRegulator.OneBody(a,b) = std::abs(matel);
      RspaceRegulator.OneBody(b,a) = std::abs(matel);
    }
  }
  // That wasn't so bad... Now for the Moshinsky fun.

  std::cout << "Done with the 1b part." << std::endl;
  std::cout << RspaceRegulator.OneBody << std::endl << std::endl;

  double sa,sb,sc,sd;
  sa=sb=sc=sd=0.5;

  int nchan = modelspace->GetNumberTwoBodyChannels();
//  modelspace.PreCalculateMoshinsky();
//  #pragma omp parallel for schedule(dynamic,1)  // In order to make this parallel, need to precompute the Moshinsky brackets.
  for (int ch=0; ch<nchan; ++ch)
  {
//     std::cout << "ch = " << ch << std::endl;
     TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
     int nkets = tbc.GetNumberKets();
     int J = tbc.J;
     for (int ibra=0;ibra<nkets;++ibra)
     {
        Ket & bra = tbc.GetKet(ibra);
        Orbit & oa = modelspace->GetOrbit(bra.p);
        Orbit & ob = modelspace->GetOrbit(bra.q);
        int na = oa.n;
        int nb = ob.n;
        int la = oa.l;
        int lb = ob.l;
        double ja = oa.j2/2.0;
        double jb = ob.j2/2.0;
        int fab = 2*na + 2*nb + la + lb;
//        if ( fab > E2max) continue;
        for (int iket=ibra;iket<nkets;++iket)
        {
//           std::cout << "   ibra,iket = " << ibra << " " << iket << std::endl;
           
           Ket & ket = tbc.GetKet(iket);
           Orbit & oc = modelspace->GetOrbit(ket.p);
           Orbit & od = modelspace->GetOrbit(ket.q);
           int nc = oc.n;
           int nd = od.n;
           int lc = oc.l;
           int ld = od.l;
           double jc = oc.j2/2.0;
           double jd = od.j2/2.0;
           int fcd = 2*nc + 2*nd + lc + ld;
//           if ( fcd > E2max) continue;

           double matel = 0;

           // First, transform to LS coupling using 9j coefficients
           for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
           {
             for (int Sab=0; Sab<=1; ++Sab)
             {
               if ( std::abs(Lab-Sab)>J or Lab+Sab<J) continue;


               double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
               if (njab == 0) continue;
               int Scd = Sab;
               int Lcd = Lab;
               double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
               if (njcd == 0) continue;

               // Next, transform to rel / com coordinates with Moshinsky tranformation
               for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
               {
                 for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
                 {
                   int Lam_cd = Lam_ab; // the projector depends only on the relative distance r, so center of mass bit must be the same. 
                   int N_cd = N_ab;
                   for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
                   {
                      if (Lab<std::abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
                      // factor to account for antisymmetrization
        
                      int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace->phase( lam_ab + Sab ))/ 2;
                      if ( asymm_factor ==0 ) continue;
        
                      int lam_cd = lam_ab; // conserve lam and Lam
                      int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation
                      int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation

                      if (n_cd < 0) continue;
                      if  (n_ab != n_cd and N_ab != N_cd) continue;
        
                      double mosh_ab = modelspace->GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
        
                      if (std::abs(mosh_ab)<1e-8) continue;
        
                      double mosh_cd = modelspace->GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                      if (std::abs(mosh_cd)<1e-8) continue;
//                      std::cout << "n_ab, lam_ab = " << n_ab << " " << lam_ab << "    n_cd, lam_cd = " << n_cd << " " << lam_cd << std::endl;

                      double me_rel = 0;
                      for (int p=oa.l; p<= oa.n+ob.n+oa.l; p++)
                      {
                        double I_p = gsl_sf_gamma_inc_P(p+1.5,x0);  // Talmi integral for a step function Theta(r0 -r );
                        double B_p = AngMom::TalmiB(oa.n, oa.l, ob.n, ob.l, p);
                        me_rel += I_p * B_p;
                      }


//                      int i_equiv = modelspace->GetOrbitIndex( n_ab, lam_ab, 2*lam_ab+1, 1 );  // the matrix element is independent of spin/isospin, so just pick one.
//                      int j_equiv = modelspace->GetOrbitIndex( n_cd, lam_cd, 2*lam_cd+1, 1 );  // the matrix element is independent of spin/isospin, so just pick one.
        
                      double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                      matel +=  me_rel * prefactor ; // don't need to do the whole Talmi integral dance. We already did that.
//                      matel +=  RspaceRegulator.OneBody(i_equiv,j_equiv) * prefactor ; // don't need to do the whole Talmi integral dance. We already did that.
        
                   } // lam_ab
                 } // Lam_ab
               } // N_ab
        
             } // Sab
           } // Lab

           RspaceRegulator.TwoBody.SetTBME(ch, ibra, iket, std::abs(matel) );
             
        } // iket
     } // ibra
  } // ch
  RspaceRegulator *= 1./r0*r0*4*3.14159;
//  RspaceRegulator *= 100 * RspaceRegulator.Norm();
//  std::cout << "Done with the 2b part. Norm is " << RspaceRegulator.OneBodyNorm() << "  " << RspaceRegulator.TwoBodyNorm() <<  std::endl;

}




