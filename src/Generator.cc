
#include "Generator.hh"
#include "imsrg_util.hh" // for VectorUnion

using namespace imsrg_util;

Generator::Generator()
  : generator_type("white"), denominator_cutoff(0)  , denominator_delta(0), denominator_delta_index(-1)
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

   if (generator_type == "wegner") // never tested, probably doesn't work.
   {
     ConstructGenerator_Wegner();
   } 
   else if (generator_type == "white")
   {
     ConstructGenerator_White();
   } 
   else if (generator_type == "atan")
   {
     ConstructGenerator_Atan();
   } 
   else if (generator_type == "shell-model")
   {
     ConstructGenerator_ShellModel();
   }
   else if (generator_type == "shell-model-atan")
   {
     ConstructGenerator_ShellModel_Atan();
   }
//   else if (generator_type == "shell-model-atan-cut")
//   {
//     ConstructGenerator_ShellModel_Atan_Cut();
//   }
   else if (generator_type == "hartree-fock")
   {
     ConstructGenerator_HartreeFock();
   }
   else
   {
      cout << "Error. Unkown generator_type: " << generator_type << endl;
   }

   Eta->profiler.timer["UpdateEta"] += omp_get_wtime() - start_time;
}


void Generator::SetDenominatorDeltaOrbit(string orb)
{
  if (orb == "all")
     SetDenominatorDeltaIndex(-12345);
  else
     SetDenominatorDeltaIndex( modelspace->GetOrbitIndex(orb) );
  cout << "Setting denominator delta orbit " << orb << " => " << modelspace->GetOrbitIndex(orb) << endl;
}

// Epstein-Nesbet energy denominators for White-type generator_types
double Generator::Get1bDenominator(int i, int j) 
{
   double ni = modelspace->GetOrbit(i).occ;
   double nj = modelspace->GetOrbit(j).occ;
   
   double denominator = H->OneBody(i,i) - H->OneBody(j,j);
   if (denominator_delta_index==-12345 or i == denominator_delta_index or j==denominator_delta_index)
     denominator += denominator_delta;
//   if (ni != nj)
     denominator += ( ni-nj ) * H->TwoBody.GetTBMEmonopole(i,j,i,j);

   if (abs(denominator)<denominator_cutoff)
     denominator *= denominator_cutoff/(abs(denominator)+1e-6);

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

   if (abs(denominator)<denominator_cutoff)
     denominator *= denominator_cutoff/(abs(denominator)+1e-6);
   return denominator;
}





// I haven't used this, so I don't know if it's right.
void Generator::ConstructGenerator_Wegner()
{
   Operator H_diag = *H;
   H_diag.ZeroBody = 0;
//   for (auto& a : modelspace->holes)
   for (auto& it_a : modelspace->holes)
   {
      index_t a = it_a.first;
      for (auto& b : modelspace->valence)
      {
         H_diag.OneBody(a,b) =0;
         H_diag.OneBody(b,a) =0;
      }
   }

   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {  // Note, should also decouple the v and q spaces
      // This is wrong. The projection operator should be different.
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_pp(), tbc.GetKetIndex_ph() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_hh(), tbc.GetKetIndex_ph() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_ph(), tbc.GetKetIndex_pp() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_ph(), tbc.GetKetIndex_hh() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_pp(), tbc.GetKetIndex_hh() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_hh(), tbc.GetKetIndex_pp() ).zeros();
   }
//   H_diag.Symmetrize();

//   *Eta = H_diag.Commutator(*H);
   *Eta = Commutator(H_diag,*H);
}




void Generator::ConstructGenerator_White()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
//   for ( auto& a : modelspace->hole_qspace)
   {
//      for ( auto& i : modelspace->particles)
//      for ( auto& i : modelspace->particle_qspace)
//      for ( auto& i : modelspace->valence)
      for ( auto& i : VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) += H->OneBody(i,a)/denominator;
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
//      for ( auto& i : modelspace->qspace)
//      {
//         double denominator = Get1bDenominator(i,a);
//         Eta->OneBody(i,a) += H->OneBody(i,a)/denominator;
//         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
//      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
//      for ( auto& iket : tbc.GetKetIndex_c_c() )
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
//         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
//         for ( auto& ibra : tbc.GetKetIndex_vv() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
//         for ( auto& ibra : tbc.GetKetIndex_v_q() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
      }
    }
}






void Generator::ConstructGenerator_Atan()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
//   for ( auto& a : modelspace->hole_qspace)
   {
//     cout << "ConstructGenerator_Atan, a = " << a << "  norbits = " << modelspace->GetNumberOrbits() << endl;
//      for ( auto& i : modelspace->particles)
//      for ( auto& i : modelspace->particle_qspace)
//      for ( auto& i : modelspace->valence)
//      cout << "    i = ";
      for ( auto& i : VectorUnion(modelspace->valence,modelspace->qspace) )
      {
//         cout << i << " ";
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) += 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
//      cout << endl;
//      for ( auto& i : modelspace->valence)
//      for ( auto& i : modelspace->qspace)
//      {
//         double denominator = Get1bDenominator(i,a);
//         Eta->OneBody(i,a) += 0.5*atan(2*H->OneBody(i,a)/denominator);
//         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
//      }

   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
//      for ( auto& iket : tbc.GetKetIndex_c_c() )
//      cout << "ch = " << ch << "(" << tbc.J << "," << tbc.parity << "," << tbc.Tz << ")  ket_cc: ";
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
//         for ( auto& ibra : tbc.GetKetIndex_q_q() )
//         for ( auto& ibra : tbc.GetKetIndex_qq() )
//         cout << iket << "->(" << modelspace->GetKet(iket).p <<"," << modelspace->GetKet(iket).q  << ") ";
         for ( auto& ibra : VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
//         for ( auto& ibra : tbc.GetKetIndex_vv() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
////         for ( auto& ibra : tbc.GetKetIndex_v_q() )
//         for ( auto& ibra : tbc.GetKetIndex_vq() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
      }
//      cout << endl;
    }
}







void Generator::ConstructGenerator_ShellModel()
{
   // One body piece -- make sure the valence one-body part is diagonal
   // no excitations out of the core

   Eta->Erase();

   for ( auto& a : VectorUnion(modelspace->core, modelspace->valence))
   {
//      Orbit& oa = modelspace->GetOrbit(a);
      for (auto& i : VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = H->OneBody(i,a)/denominator;
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }


//   index_t norb = modelspace->GetNumberOrbits();
//   for (index_t i=0; i<norb; ++i)
//   {
//     for (index_t j=i;j<norb; ++j)
//     {
//         double denominator = Get1bDenominator(i,j);
//         Eta->OneBody(i,j) += H->OneBody(i,j)/denominator;
//         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
//     }
//   }
//   for ( auto& i : modelspace->particles )
//   {
////      for ( auto& j : modelspace->holes )
////      {
////         double denominator = Get1bDenominator(i,j);
////         Eta->OneBody(i,j) += H->OneBody(i,j)/denominator;
////         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
////      }
////      for ( auto& j : modelspace->particles )
//      for ( auto& j : VectorUnion( modelspace->holes, modelspace->particles ) )
//      {
//         if (i==j) continue;
//         double denominator = Get1bDenominator(i,j);
//         Eta->OneBody(i,j) += H->OneBody(i,j)/denominator;
//         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
//      }
//   }
//   for ( auto& i : modelspace->holes )
//   {
//      for ( auto& j : modelspace->holes )
//      {
//         if (i==j) continue;
//         double denominator = -Get1bDenominator(i,j);
//         Eta->OneBody(i,j) += H->OneBody(i,j)/denominator;
//         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
//      }
//   }


   // Two body piece 
   // we want to drive to zero any terms that take 
   //   |hh> to |pp>
   //   |vh> to |pp>
   //   |vv> to |qq> or |vq>

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& ETA2 = Eta->TwoBody.GetMatrix(ch);
      auto& H2 = H->TwoBody.GetMatrix(ch);


      // Decouple vv states from pq states
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         // < qq' | vv' >
//         for ( auto& ibra : tbc.GetKetIndex_q_q() ) 
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv()) ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

//         // < vq | vv' > 
//         for ( auto& ibra : tbc.GetKetIndex_v_q() ) 
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
      }


      // Decouple hh states

//      for ( auto& iket : tbc.GetKetIndex_c_c() )
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {

         // < qq' | hh' >
//         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_vv() )  )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         // < vq | hh' >
//         for ( auto& ibra : tbc.GetKetIndex_v_q() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
//
//         // < vv | hh' >
//         for ( auto& ibra : tbc.GetKetIndex_vv() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }

      }


      // Decouple vh states

//      for ( auto& iket : tbc.GetKetIndex_v_c() )
      for ( auto& iket : tbc.GetKetIndex_vc() )
      {
         // < qq | vh >
//         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_vv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

//         // < vq | vh >
//         for ( auto& ibra : tbc.GetKetIndex_v_q() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
//
//         // < vv | vh >
//         for ( auto& ibra : tbc.GetKetIndex_vv() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }

      }


    }
}





void Generator::ConstructGenerator_ShellModel_Atan()
{
   // One body piece -- make sure the valence one-body part is diagonal
//   unsigned int norbits = modelspace->GetNumberOrbits();
//   for ( auto& a : modelspace->valence)
   for ( auto& a : VectorUnion(modelspace->core, modelspace->valence))
   {
//      Orbit& oa = modelspace->GetOrbit(a);
//      for (unsigned int j=0; j<norbits; ++j)
//      for (auto j : H->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}))
      for (auto& i : VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
//         Eta->OneBody(i,j) += 0.5*atan(2*H->OneBody(i,j)/denominator);
         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }
   // Also make sure the core is decoupled
//   for ( auto& a : modelspace->hole_qspace)
//   for ( auto& a : modelspace->core)
//   {
////      for ( auto& i : VectorUnion( modelspace->valence, modelspace->qspace ) )
//      for ( auto& i : modelspace->qspace)
//      {
//         double denominator = Get1bDenominator(i,a);
////         Eta->OneBody(i,a) += 0.5*atan(2*H->OneBody(i,a)/denominator);
//         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
//         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
//      }
////      for ( auto& i : modelspace->valence)
////      {
////         double denominator = Get1bDenominator(i,a);
////         Eta->OneBody(i,a) += 0.5*atan(2*H->OneBody(i,a)/denominator);
////         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
////      }
//
//   }


   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H->TwoBody.GetMatrix(ch);

      // Decouple the core
//      for ( auto& iket : tbc.GetKetIndex_c_c() )
//      for ( auto& iket : tbc.GetKetIndex_cc() )
      for ( auto& iket : VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
//         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
//         for ( auto& ibra : tbc.GetKetIndex_vv() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
//         for ( auto& ibra : tbc.GetKetIndex_v_q() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
//         for ( auto& ibra : tbc.GetKetIndex_q_q() ) 
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator) ;
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator) ;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
//         for ( auto& ibra : tbc.GetKetIndex_v_q() ) 
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
 
      }

      // Decouple vh states

//      for ( auto& iket : tbc.GetKetIndex_v_c() )
//      for ( auto& iket : tbc.GetKetIndex_vc() )
//      {
////         for ( auto& ibra : tbc.GetKetIndex_q_q() )
//         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
//
//         for ( auto& ibra : tbc.GetKetIndex_v_q() )
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
//        for ( auto& ibra : tbc.GetKetIndex_vv() ) 
//         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
//            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
//         }
//      }

    }
}


/*
// Obsolete
void Generator::ConstructGenerator_ShellModel_Atan_Cut()
{
   // One body piece -- make sure the valence one-body part is diagonal
   unsigned int norbits = modelspace->GetNumberOrbits();
   for ( auto& i : modelspace->valence)
   {
      for (unsigned int j=0; j<norbits; ++j)
      {
         if (i==j) continue;
         double denominator = Get1bDenominator(i,j);
         if (abs(denominator) > denominator_cutoff);
           Eta->OneBody(i,j) += 0.5*atan(2*H->OneBody(i,j)/denominator);
         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
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
      for ( auto& iket : tbc.GetKetIndex_c_c() )
      {
         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
         if (abs(denominator) > denominator_cutoff);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_vv() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
         if (abs(denominator) > denominator_cutoff);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_v_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
         if (abs(denominator) > denominator_cutoff);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         for ( auto& ibra : tbc.GetKetIndex_q_q() ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
         if (abs(denominator) > denominator_cutoff);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator) ;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_v_q() ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
         if (abs(denominator) > denominator_cutoff);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
 
      }

      // Decouple vh states

      for ( auto& iket : tbc.GetKetIndex_v_c() )
      {
         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
         if (abs(denominator) > denominator_cutoff);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         for ( auto& ibra : tbc.GetKetIndex_v_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
         if (abs(denominator) > denominator_cutoff);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
        for ( auto& ibra : tbc.GetKetIndex_vv() ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
         if (abs(denominator) > denominator_cutoff);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}

*/







void Generator::ConstructGenerator_HartreeFock()
{
   Eta->SetParticleRank(1);
   // One body piece -- eliminate ph bits
   unsigned int norbits = modelspace->GetNumberOrbits();
   for (unsigned int i=0;i<norbits;++i)
   {
      for (unsigned int j=0; j<i; ++j)
      {
         double denominator = Get1bDenominator(i,j);
         Eta->OneBody(i,j) += H->OneBody(i,j)/denominator;
         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
      }
   } 
}












