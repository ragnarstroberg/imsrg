
#include "Generator.hh"
#include "imsrg_util.hh" // for VectorUnion

#include "omp.h"

using namespace imsrg_util;

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

        if (generator_type == "wegner")           ConstructGenerator_Wegner(); // never tested, probably doesn't work.
   else if (generator_type == "white")            ConstructGenerator_White();
   else if (generator_type == "atan")             ConstructGenerator_Atan();
   else if (generator_type == "imaginary-time")            ConstructGenerator_ImaginaryTime();
   else if (generator_type == "shell-model")      ConstructGenerator_ShellModel();
   else if (generator_type == "shell-model-atan") ConstructGenerator_ShellModel_Atan();
   else if (generator_type == "shell-model-atan-npnh") ConstructGenerator_ShellModel_Atan_NpNh();
   else if (generator_type == "shell-model-imaginary-time")            ConstructGenerator_ShellModel_ImaginaryTime();
   else if (generator_type == "hartree-fock")     ConstructGenerator_HartreeFock();
   else if (generator_type == "1PA")     ConstructGenerator_1PA();
   else
   {
      cout << "Error. Unkown generator_type: " << generator_type << endl;
   }
   Eta->profiler.timer["UpdateEta"] += omp_get_wtime() - start_time;
}


// Old method used to test some things out. Not typically used.
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



// I haven't used this, so I don't know if it's right.
void Generator::ConstructGenerator_Wegner()
{
   Operator H_diag = *H;
   H_diag.ZeroBody = 0;
   for (auto& a : modelspace->holes)
   {
//      index_t a = it_a.first;
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
   Eta->SetToCommutator(H_diag,*H);
//   *Eta = Commutator(H_diag,*H);
}




void Generator::ConstructGenerator_White()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : VectorUnion(modelspace->valence,modelspace->qspace) )
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
      for ( auto& iket : tbc.GetKetIndex_cc() )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
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
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : VectorUnion(modelspace->valence,modelspace->qspace) )
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
         for ( auto& ibra : VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}



/// Imaginary time generator \f[ \eta = sgn(\Delta) h_{od} \]
void Generator::ConstructGenerator_ImaginaryTime()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : VectorUnion(modelspace->valence,modelspace->qspace) )
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
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            if (denominator==0) denominator = 1;
            ETA2(ibra,iket) += H2(ibra,iket) * denominator / std::abs(denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}



void Generator::ConstructGenerator_ShellModel()
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : VectorUnion( modelspace->valence, modelspace->qspace ) )
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
      for ( auto& iket : VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
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
   for ( auto& a : VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : VectorUnion( modelspace->valence, modelspace->qspace ) )
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
      for ( auto& iket : VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator) ;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}





/// Imaginary time generator for a valence space
void Generator::ConstructGenerator_ShellModel_ImaginaryTime()
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : VectorUnion( modelspace->valence, modelspace->qspace ) )
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
      for ( auto& iket : VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
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
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
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

//  cout << "In ShellModel_Atat_NpNh, adding to Eta" << endl;
  // decouple f_cc'
  for ( auto& c : modelspace->core )
  {
   for ( auto& cprime : modelspace->core )
   {
     if (cprime<=c) continue;
     double denominator = Get1bDenominator(c,cprime);
     Eta->OneBody(c,cprime) = 0.5*atan(2*H->OneBody(c,cprime) / denominator );
//     cout << "c,cprime = " << c << " " << cprime << "  etacc' = " << Eta->OneBody(c,cprime) << endl;
     Eta->OneBody(cprime,c) = - Eta->OneBody(c,cprime);
   }
  }

  int nchan = modelspace->GetNumberTwoBodyChannels();
  for (int ch=0;ch<nchan;++ch)
  {
     TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
//     cout << "ch = " << ch << "  vc size = " << tbc.GetKetIndex_vc().size() << "   qc size = " << tbc.GetKetIndex_qc().size() << endl;
     arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
     arma::mat& H2 =  H->TwoBody.GetMatrix(ch);
  // decouple Gamma_qcvc'
     for (auto& iket : tbc.GetKetIndex_vc())
     {
       Ket& ket = tbc.GetKet(iket);
       for (auto& ibra : tbc.GetKetIndex_qc())
       {
         Ket& bra = tbc.GetKet(ibra);
//         cout << bra.p << " " << bra.q << " " << ket.p << " " << ket.q << endl;
         if ((ket.p==bra.p) or (ket.p==bra.q) or (ket.q==bra.p) or (ket.q==bra.q) ) continue;
         double denominator = Get2bDenominator(ch,ibra,iket);
         ETA2(ibra,iket) = 0.5*atan( 2*H2(ibra,iket) / denominator );
         ETA2(iket,ibra) = -ETA2(ibra,iket);
//         cout << "   qcvc': " << ket.p << " " << ket.q << " " << bra.p << " " << bra.q << "    " << ETA2(ibra,iket) << endl;
       }
     }
  // decouple Gamma_pcc'c''
     for (auto& iket : tbc.GetKetIndex_cc())
     {
       Ket& ket = tbc.GetKet(iket);
       for (auto& ibra : VectorUnion( tbc.GetKetIndex_vc(), tbc.GetKetIndex_qc() )  )
       {
         Ket& bra = tbc.GetKet(ibra);
         if ((ket.p==bra.p) or (ket.p==bra.q) or (ket.q==bra.p) or (ket.q==bra.q) ) continue;
         double denominator = Get2bDenominator(ch,ibra,iket);
         ETA2(ibra,iket) = 0.5*atan( 2*H2(ibra,iket) / denominator );
         ETA2(iket,ibra) = -ETA2(ibra,iket);
       }
     }
  }
//  cout << "all done" << endl;
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




// So far this is useless
void Generator::ConstructGenerator_1PA()
{
  ConstructGenerator_Atan();
  int nchan = modelspace->GetNumberTwoBodyChannels();
  for (int ch=0;ch<nchan;++ch)
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
    arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
    arma::mat& H2 =  H->TwoBody.GetMatrix(ch);
    // decouple Gamma_ppph'
    for (auto& iket : tbc.GetKetIndex_ph())
    {
      Ket& ket = tbc.GetKet(iket);
      if ((2*ket.oq->n + ket.oq->l)<3) continue;
      for (auto& ibra : tbc.GetKetIndex_pp())
      {
        Ket& bra = tbc.GetKet(ibra);
//        cout << bra.p << " " << bra.q << " " << ket.p << " " << ket.q << endl;
        if ((ket.p==bra.p) or (ket.p==bra.q) or (ket.q==bra.p) or (ket.q==bra.q) ) continue;
        double denominator = Get2bDenominator(ch,ibra,iket);
        ETA2(ibra,iket) = 0.5*atan( 2*H2(ibra,iket) / denominator );
        ETA2(iket,ibra) = -ETA2(ibra,iket);
      }
    }
  }
}








