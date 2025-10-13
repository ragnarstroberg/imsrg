#include "GeneratorPV.hh"
// #include "Generator.hh"
#include "PhysicalConstants.hh" // for HBARC and M_NUCLEON
// #include "Commutator.hh"
// #include "Operator.hh"
// GeneratorPV::GeneratorPV():
//              Generator()
//{}

using PhysConst::HBARC;
using PhysConst::M_NUCLEON;

void GeneratorPV::Update(Operator &H_s, Operator &HPV_s, Operator &Eta_s, Operator &EtaPV_s)
{
   Eta_s.Erase();
   EtaPV_s.Erase();
   AddToEtaPV(H_s, HPV_s, Eta_s, EtaPV_s);
}

void GeneratorPV::AddToEtaPV(Operator &H_s, Operator &HPV_s, Operator &Eta_s, Operator &EtaPV_s)
{
   double start_time = omp_get_wtime();
   H = &H_s;
   Eta = &Eta_s;
   V = &HPV_s;
   V -> MakeNotReduced();
   Etapv = &EtaPV_s;
   Etapv->MakeNotReduced();
   //   modelspace = H->GetModelSpace();
   if (generator_type == "wegner")
      ConstructGeneratorPV_SingleRef(wegner_func); // never tested, probably doesn't work.
   else if (generator_type == "white")
      ConstructGeneratorPV_SingleRef(white_func);
   else if (generator_type == "atan")
      ConstructGeneratorPV_SingleRef(atan_func);
   else if (generator_type == "imaginary-time")
      ConstructGeneratorPV_SingleRef(imaginarytime_func);
   else if (generator_type == "qtransfer-atan")
      ConstructGeneratorPV_SingleRef(qtransferatan1_func);
   else if (generator_type == "shell-model-wegner")
      ConstructGeneratorPV_ShellModel(wegner_func); // never tested, probably doesn't work.
   else if (generator_type == "shell-model")
   {
      ConstructGeneratorPV_ShellModel(white_func);
   }
      
   else if (generator_type == "shell-model-atan")
      ConstructGeneratorPV_ShellModel(atan_func);
   //   else if (generator_type == "shell-model-atan-npnh")        ConstructGeneratorPV_ShellModel_NpNh(atan_func);
   else if (generator_type == "shell-model-imaginary-time")
      ConstructGeneratorPV_ShellModel(imaginarytime_func);
   //   else if (generator_type == "hartree-fock")                 ConstructGeneratorPV_HartreeFock();
   //   else if (generator_type == "1PA")                          ConstructGeneratorPV_1PA(atan_func);
   else if (generator_type.find("qtransfer-atan") != std::string::npos)
   {
      int n;
      std::istringstream(generator_type.substr(generator_type.find("_") + 1)) >> n;
      std::function<double(double, double)> qtransferatanN_func = [n](double Hod, double denom)
      { return pow(std::abs(denom) * M_NUCLEON / HBARC / HBARC, 0.5 * n) * atan_func(Hod, denom); };
      //     ConstructGenerator_QTransferAtan(n);
      ConstructGeneratorPV_SingleRef(qtransferatanN_func);
   }
   //   else if (generator_type == "rspace")                       ConstructGenerator_Rspace();
   else
   {
      std::cout << "Error. Unkown generator_type: " << generator_type << std::endl;
   }
   V->MakeReduced();
   Etapv->MakeReduced();
   Eta->profiler.timer["UpdateEta"] += omp_get_wtime() - start_time;
   Etapv->profiler.timer["UpdateEta"] += omp_get_wtime() - start_time;
}

void GeneratorPV::ConstructGeneratorPV_SingleRef(std::function<double(double, double)> &etafunc)
{
   double start_time = omp_get_wtime();
   // If V is hermitian, than it is real and Etapv is antihermitian.
   // If V is not antihermitian, then it is imaginary and Etapv is hermitian.
   int herm = V->IsHermitian() ? -1 : +1;
   // One body piece -- eliminate ph bits
   for (auto &a : H->modelspace->core)
   {
      for (auto &i : VectorUnion(H->modelspace->valence, H->modelspace->qspace))
      {
         
         double denominator = Get1bDenominator(i, a);
         Eta->OneBody(i, a) = etafunc(H->OneBody(i, a), denominator);
         Eta->OneBody(a, i) = -Eta->OneBody(i, a);
         Etapv->OneBody(i, a) = etafunc(V->OneBody(i, a), denominator); // commented Beatriz 27/08/24 to make etapv(1b)=0
         Etapv->OneBody(a, i) = herm*Etapv->OneBody(i, a); // old version
         Eta->profiler.timer["UpdateEta1beta"] += omp_get_wtime() - start_time;
         Etapv->profiler.timer["UpdateEta1betapv"] += omp_get_wtime() - start_time;

         //         std::cout<<"numeratorH="<<H->OneBody(i,a)<<__LINE__<<std::endl;
         //         std::cout<<"i orbit="<< i << " " <<__LINE__<<std::endl;
         //         std::cout<<"a orbit="<< a << " " <<__LINE__<<std::endl;
      }
   }
   if (only_1b_eta)
      return;
   for (auto &iter : Eta->TwoBody.MatEl)
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      //      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      TwoBodyChannel &tbc_bra = H->modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = H->modelspace->GetTwoBodyChannel(ch_ket);
      //      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat &ETA2 = iter.second;
      //      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      arma::mat &H2 = H->TwoBody.GetMatrix(ch_bra, ch_ket);

      for (auto &iket : tbc_ket.GetKetIndex_cc()) // cc means core-core ('holes' refer to the reference state)
      {
         for (auto &ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv()))
         {
            double denominator = Get2bDenominator(ch_bra, ch_ket, ibra, iket);
            ETA2(ibra, iket) = etafunc(H2(ibra, iket), denominator);
            ETA2(iket, ibra) = -ETA2(ibra, iket); // Eta needs to be antisymmetric
            // std::cout << ch_bra << ch_ket << ibra << iket << " " << ETA2(ibra, iket)<<" "<< Get2bDenominator(ch_bra, ch_ket, ibra, iket) << std::endl;
            //   std::cout << __func__ << "  line " << __LINE__ << " bra,ket " << bra.p << " " << bra.q << " , " << ket.p << " " << ket.q  << "  J = " << tbc_bra.J << "   numerator /denom = " << H2(ibra,iket) << " / " << denominator << std::endl;//added Beatriz 11/09/24
         }
      }
      Eta->profiler.timer["UpdateEta2beta"] += omp_get_wtime() - start_time;
   }
   
   for (auto &iter : Etapv->TwoBody.MatEl)
   {
      Etapv->profiler.timer["UpdateEtapv2biter"] += omp_get_wtime() - start_time;
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      //      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      TwoBodyChannel &tbc_bra = V->modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = V->modelspace->GetTwoBodyChannel(ch_ket);

      Etapv->profiler.timer["UpdateEtapv2bodychannelbraket"] += omp_get_wtime() - start_time;
      //      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat &ETAPV2 = iter.second;
      //      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
      arma::mat &V2 = V->TwoBody.GetMatrix(ch_bra, ch_ket);
      Etapv->profiler.timer["UpdateEtapvV"] += omp_get_wtime() - start_time;

      for (auto &iket : tbc_ket.GetKetIndex_cc()) // cc means core-core ('holes' refer to the reference state)
      {
         for (auto &ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv()))
         {
            double denominator = Get2bDenominator(ch_bra, ch_ket, ibra, iket);
            ETAPV2(ibra, iket) =  etafunc(V2(ibra, iket), denominator);
         }
      }
      Etapv->profiler.timer["UpdateEta2betapv"] += omp_get_wtime() - start_time;
   }
}

void GeneratorPV::ConstructGeneratorPV_ShellModel(std::function<double(double, double)> &eta_func)
{
   // If V is hermitian, than it is real and Etapv is antihermitian.
   // If V is not antihermitian, then it is imaginary and Etapv is hermitian.
   int herm = V->IsHermitian() ? -1 : +1;
   // One body piece -- make sure the valence one-body part is diagonal
   for (auto &a : VectorUnion(H->modelspace->core, H->modelspace->valence))
   {
      for (auto &i : VectorUnion(H->modelspace->valence, H->modelspace->qspace))
      {
         if (i == a)
            continue;
         double denominator = Get1bDenominator(i, a);
         Eta->OneBody(i, a) = eta_func(H->OneBody(i, a), denominator);
         Eta->OneBody(a, i) = -Eta->OneBody(i, a);
         Etapv->OneBody(i, a) = eta_func(V->OneBody(i, a), denominator);
         Etapv->OneBody(a, i) = herm*Etapv->OneBody(i, a); // old version
         // std::cout << "  looping in generatorPV,Eta 1b part = " << Eta->OneBody(i,a) << std::endl;
         // std::cout << "  looping in generatorPV,Etapv 1b part = " << Etapv->OneBody(i,a) << std::endl;
      }
   }
   if (only_1b_eta)
      return;
   for (auto &iter : Eta->TwoBody.MatEl)
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel &tbc_bra = H->modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = H->modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat &ETA2 = iter.second;
      arma::mat &H2 = H->TwoBody.GetMatrix(ch_bra, ch_ket);

      // Decouple the core
      for (auto &iket : VectorUnion(tbc_ket.GetKetIndex_cc(), tbc_ket.GetKetIndex_vc())) // cc means core-core, vc means valence-core
      {
         for (auto &ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv()))
         {
            double denominator = Get2bDenominator(ch_bra, ch_ket, ibra, iket);
            ETA2(ibra, iket) = eta_func(H2(ibra, iket), denominator);
            ETA2(iket, ibra) = -ETA2(ibra, iket); // Eta needs to be antisymmetric
            // std::cout << "  looping in generatorPV,Eta 2b part = " << ETA2(ibra,iket) << std::endl;
         }
      }

      // Decouple the valence space
      for (auto &iket : tbc_ket.GetKetIndex_vv())
      {
         //         auto& ket = tbc.GetKet(iket);
         for (auto &ibra : VectorUnion(tbc_bra.GetKetIndex_qv(), tbc_bra.GetKetIndex_qq()))
         {
          
            double denominator = Get2bDenominator(ch_bra, ch_ket, ibra, iket);
            ETA2(ibra, iket) = eta_func(H2(ibra, iket), denominator);
            ETA2(iket, ibra) = -ETA2(ibra, iket); // Eta needs to be antisymmetric
            // if (ETA2(ibra, iket) != 0)
            // {
            //    std::cout << std::setprecision(6) << std::fixed;
            //    std::cout << "  looping in generatorPV,Eta 2b part = " << ETA2(ibra, iket) << std::endl;
            // }
            
         }
      }
   }
   for (auto &iter : Etapv->TwoBody.MatEl)
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel &tbc_bra = V->modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = V->modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat &ETAPV2 = iter.second;
      arma::mat &V2 = V->TwoBody.GetMatrix(ch_bra, ch_ket);
      // Decouple the core
      for (auto &iket : VectorUnion(tbc_ket.GetKetIndex_cc(),tbc_ket.GetKetIndex_vc()))
      {
         for (auto &ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv()))
         {
            // std::cout << "ch_bra=" << ch_bra << " ch_ket=" << ch_ket << " ibra=" << ibra << " iket=" << iket << std::endl;
            // std::cout<< ETAPV2 <<std::endl;
            double denominator = Get2bDenominator(ch_bra, ch_ket, ibra, iket);
            ETAPV2(ibra, iket) = eta_func(V2(ibra, iket), denominator);
         }
      }

      // Decouple the valence space
      for (auto &iket : tbc_ket.GetKetIndex_vv())
      {
         for (auto &ibra : VectorUnion(tbc_bra.GetKetIndex_qv(), tbc_bra.GetKetIndex_qq()))
         {
            // std::cout << "ch_bra=" << ch_bra << " ch_ket=" << ch_ket << " ibra=" << ibra << " iket=" << iket << std::endl;
            // std::cout << ETAPV2 <<std::endl; 
            double denominator = Get2bDenominator(ch_bra, ch_ket, ibra, iket);
            ETAPV2(ibra, iket) = eta_func(V2(ibra, iket), denominator);
         }
      }
   }
}
