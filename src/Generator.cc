
#include "Generator.hh"

Generator::Generator()
  : generator_type("white"), denominator_delta(0), denominator_delta_index(-1)
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
   else if (generator_type == "shell-model-atan-cut")
   {
     ConstructGenerator_ShellModel_Atan_Cut();
   }
   else if (generator_type == "hartree-fock")
   {
     ConstructGenerator_HartreeFock();
   }
   else
   {
      cout << "Error. Unkown generator_type: " << generator_type << endl;
   }

   Eta->timer["UpdateEta"] += omp_get_wtime() - start_time;
}


void Generator::SetDenominatorDeltaOrbit(string orb)
{
  if (orb == "all")
     SetDenominatorDeltaIndex(-12345);
  else
     SetDenominatorDeltaIndex( modelspace->GetOrbitIndex(orb) );
}

// Epstein-Nesbet energy denominators for White-type generator_types
double Generator::Get1bDenominator(int i, int j) 
{
   int ni = modelspace->GetOrbit(i).ph;
   int nj = modelspace->GetOrbit(j).ph;
   
   double denominator = H->OneBody(i,i) - H->OneBody(j,j);
   if (denominator_delta_index==-12345 or i == denominator_delta_index or j==denominator_delta_index)
     denominator += denominator_delta;
   if (ni != nj)
     denominator += ( ni-nj ) * H->TwoBody.GetTBMEmonopole(i,j,i,j);
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
   double denominator = H->OneBody(i,i)+ H->OneBody(j,j) - H->OneBody(k,k) - H->OneBody(l,l) + denominator_delta;
   if (denominator_delta_index == -12345) denominator += denominator_delta;
   int ni = bra.op->ph;
   int nj = bra.oq->ph;
   int nk = ket.op->ph;
   int nl = ket.oq->ph;

   denominator       += ( 1-ni-nj ) * H->TwoBody.GetTBMEmonopole(i,j,i,j); // pp'pp'
   denominator       -= ( 1-nk-nl ) * H->TwoBody.GetTBMEmonopole(k,l,k,l); // hh'hh'
   denominator       += ( ni-nk ) * H->TwoBody.GetTBMEmonopole(i,k,i,k); // phph
   denominator       += ( ni-nl ) * H->TwoBody.GetTBMEmonopole(i,l,i,l); // ph'ph'
   denominator       += ( nj-nk ) * H->TwoBody.GetTBMEmonopole(j,k,j,k); // p'hp'h
   denominator       += ( nj-nl ) * H->TwoBody.GetTBMEmonopole(j,l,j,l); // p'h'p'h'

   return denominator;
}





// I haven't used this, so I don't know if it's right.
void Generator::ConstructGenerator_Wegner()
{
   Operator H_diag = *H;
   H_diag.ZeroBody = 0;
   for (auto& a : modelspace->holes)
   {
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




/*
void Generator::ConstructGenerator_White()
{
   // One body piece -- eliminate ph bits
   for ( auto& i : modelspace->particles)
   {
      for ( auto& a : modelspace->holes)
      {
         //double denominator = GetEpsteinNesbet1bDenominator(i,a);
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = H->OneBody(i,a)/denominator;
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits. Note that the hh'hp pieces are accounted
   // for in the normal ordered one-body part.
   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto& ETA2 = Eta->TwoBody.GetMatrix(ch);
      auto& H2 = H->TwoBody.GetMatrix(ch);
      for ( auto& ibra : tbc.GetKetIndex_pp() )
      {
         for ( auto& iket : tbc.GetKetIndex_hh() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) =  H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}
*/

void Generator::ConstructGenerator_White()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->hole_qspace)
   {
//      for ( auto& i : modelspace->particles)
      for ( auto& i : modelspace->particle_qspace)
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) += H->OneBody(i,a)/denominator;
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
      for ( auto& i : modelspace->valence)
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) += H->OneBody(i,a)/denominator;
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_c_c() )
      {
         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_vv() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_v_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}






void Generator::ConstructGenerator_Atan()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->hole_qspace)
   {
//      for ( auto& i : modelspace->particles)
      for ( auto& i : modelspace->particle_qspace)
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) += 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
      for ( auto& i : modelspace->valence)
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) += 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }

   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta->nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_c_c() )
      {
         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_vv() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_v_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}







void Generator::ConstructGenerator_ShellModel()
{
   // One body piece -- make sure the valence one-body part is diagonal
   // no excitations out of the core

   Eta->Erase();

   for ( auto& i : modelspace->particles )
   {
      for ( auto& j : modelspace->holes )
      {
         double denominator = Get1bDenominator(i,j);
         Eta->OneBody(i,j) += H->OneBody(i,j)/denominator;
         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
      }
      for ( auto& j : modelspace->particles )
      {
         if (i==j) continue;
         double denominator = Get1bDenominator(i,j);
         Eta->OneBody(i,j) += H->OneBody(i,j)/denominator;
         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
      }
   }
   for ( auto& i : modelspace->holes )
   {
      for ( auto& j : modelspace->holes )
      {
         if (i==j) continue;
         double denominator = -Get1bDenominator(i,j);
         Eta->OneBody(i,j) += H->OneBody(i,j)/denominator;
         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
      }
   }


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
         for ( auto& ibra : tbc.GetKetIndex_q_q() ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vq | vv' > 
         for ( auto& ibra : tbc.GetKetIndex_v_q() ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }


      // Decouple hh states

      for ( auto& iket : tbc.GetKetIndex_c_c() )
      {

         // < qq' | hh' >
         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         // < vq | hh' >
         for ( auto& ibra : tbc.GetKetIndex_v_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vv | hh' >
         for ( auto& ibra : tbc.GetKetIndex_vv() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }


      // Decouple vh states

      for ( auto& iket : tbc.GetKetIndex_v_c() )
      {
         // < qq | vh >
         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vq | vh >
         for ( auto& ibra : tbc.GetKetIndex_v_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vv | vh >
         for ( auto& ibra : tbc.GetKetIndex_vv() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }


    }
}





void Generator::ConstructGenerator_ShellModel_Atan()
{
   // One body piece -- make sure the valence one-body part is diagonal
   unsigned int norbits = modelspace->GetNumberOrbits();
   for ( auto& i : modelspace->valence)
   {
      for (unsigned int j=0; j<norbits; ++j)
      {
         if (i==j) continue;
         double denominator = Get1bDenominator(i,j);
         Eta->OneBody(i,j) += 0.5*atan(2*H->OneBody(i,j)/denominator);
         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
      }
   }
   // Also make sure the core is decoupled
   for ( auto& a : modelspace->hole_qspace)
   {
      for ( auto& i : modelspace->particle_qspace)
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) += 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
      for ( auto& i : modelspace->valence)
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) += 0.5*atan(2*H->OneBody(i,a)/denominator);
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
      for ( auto& iket : tbc.GetKetIndex_c_c() )
      {
         for ( auto& ibra : tbc.GetKetIndex_q_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_vv() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_v_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
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
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator) ;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         for ( auto& ibra : tbc.GetKetIndex_v_q() ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
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
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         for ( auto& ibra : tbc.GetKetIndex_v_q() )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
        for ( auto& ibra : tbc.GetKetIndex_vv() ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) += 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}


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












