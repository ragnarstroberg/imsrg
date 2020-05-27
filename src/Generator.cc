
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
using PhysConst::M_NUCLEON;
using PhysConst::HBARC;

std::function<double(double,double)> Generator::wegner_func = [] (double Hod, double denom){ return Hod * denom;};
std::function<double(double,double)> Generator::white_func = [] (double Hod, double denom){ return Hod / denom;};
std::function<double(double,double)> Generator::atan_func = [] (double Hod, double denom){ return 0.5 * atan(2*Hod / denom);};
std::function<double(double,double)> Generator::imaginarytime_func = [] (double Hod, double denom){ return Hod * ( (denom>0)-(denom<0) ) ;};//signum function

std::function<double(double,double)> Generator::qtransferatan1_func = [](double Hod, double denom){return pow(std::abs(denom)*M_NUCLEON/HBARC/HBARC, 0.5*1) * atan_func(Hod, denom);};

Generator::Generator()
  : generator_type("white"), denominator_cutoff(1e-6)  , denominator_delta(0), denominator_delta_index(-1), only_2b_eta(false)
{}


void Generator::Update(Operator * H_s, Operator * Eta_s)
{
   Eta_s->EraseOneBody();
   Eta_s->EraseTwoBody();
   Eta_s->EraseThreeBody();// do we need this? Maybe?
   AddToEta(H_s,Eta_s);
}


void Generator::AddToEta(Operator * H_s, Operator * Eta_s)
{
   double start_time = omp_get_wtime();
   H = H_s;
   Eta = Eta_s;
   modelspace = H->GetModelSpace();

        if (generator_type == "wegner")                       ConstructGenerator_SingleRef(wegner_func); // never tested, probably doesn't work.
   else if (generator_type == "white")                        ConstructGenerator_SingleRef(white_func);
   else if (generator_type == "atan")                         ConstructGenerator_SingleRef(atan_func) ;
   else if (generator_type == "imaginary-time")               ConstructGenerator_SingleRef(imaginarytime_func);
   else if (generator_type == "qtransfer-atan")               ConstructGenerator_SingleRef(qtransferatan1_func);
   else if (generator_type == "shell-model-wegner")           ConstructGenerator_ShellModel(wegner_func); // never tested, probably doesn't work.
   else if (generator_type == "shell-model")                  ConstructGenerator_ShellModel(white_func);
   else if (generator_type == "shell-model-atan")             ConstructGenerator_ShellModel(atan_func);
   else if (generator_type == "shell-model-atan-npnh")        ConstructGenerator_ShellModel_NpNh(atan_func);
   else if (generator_type == "shell-model-imaginary-time")   ConstructGenerator_ShellModel(imaginarytime_func);
   else if (generator_type == "hartree-fock")                 ConstructGenerator_HartreeFock();
   else if (generator_type == "1PA")                          ConstructGenerator_1PA(atan_func);
   else if (generator_type.find("qtransfer-atan") != std::string::npos )
   {
     int n;
     std::istringstream( generator_type.substr( generator_type.find("_")+1) ) >> n;
     std::function<double(double,double)> qtransferatanN_func = [n](double Hod, double denom){return pow(std::abs(denom)*M_NUCLEON/HBARC/HBARC, 0.5*n) * atan_func(Hod, denom);};
//     ConstructGenerator_QTransferAtan(n);
     ConstructGenerator_SingleRef( qtransferatanN_func );
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





void Generator::ConstructGenerator_SingleRef(std::function<double (double,double)>& etafunc )
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
//         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(i,a) = etafunc( H->OneBody(i,a), denominator);
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
//            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(ibra,iket) = etafunc( H2(ibra,iket), denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }

    if ( Eta->GetParticleRank()>2 and H->GetParticleRank()>2 and not only_2b_eta )
    {
       double t_start = omp_get_wtime();
//       ConstructGenerator_Atan_3body();
       ConstructGenerator_SingleRef_3body( etafunc );
       H->profiler.timer["Update Eta 3body"] += omp_get_wtime() - t_start;
    }// if particle rank >3
}





// Off-diagonal pieces are <abc|ijk> = <ppp|ccc> where c is core and p is either valence or q (that is, not core).
//void Generator::ConstructGenerator_Atan_3body()
void Generator::ConstructGenerator_SingleRef_3body(std::function<double (double,double)>& etafunc )
{
     double t_start = omp_get_wtime();
     size_t ncore = modelspace->core.size();
     std::vector<size_t> corevec;
     for (auto a : modelspace->core) corevec.push_back(a);
     std::map<int,double> e_fermi = modelspace->GetEFermi();
     std::cout << __func__ << "  looping in generator 3-body part .  Size of H3 = " << H->ThreeBodyNorm() << std::endl;
//    for (auto a : modelspace->core )
     size_t nch3 = modelspace->GetNumberThreeBodyChannels();
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t ch3=0; ch3<nch3; ch3++)
    {
      ThreeBodyChannel& Tbc = modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      for (size_t ibra=0; ibra<nkets3; ibra++)
      {
        Ket3& bra = Tbc.GetKet(ibra);
         if ( (  (bra.op->cvq==0) or (bra.oq->cvq==0) or (bra.oR->cvq==0) ) ) continue; //cvq==0 means core orbit


        double d_ea = std::abs( 2*bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
        double d_eb = std::abs( 2*bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
        double d_ec = std::abs( 2*bra.oR->n + bra.oR->l - e_fermi[bra.oR->tz2]);
        double occnat_a = bra.op->occ_nat;
        double occnat_b = bra.oq->occ_nat;
        double occnat_c = bra.oR->occ_nat;
        if ( d_ea + d_eb + d_ec > modelspace->GetdE3max() ) continue;
        if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_c*(1-occnat_c) ) < modelspace->GetOccNat3Cut() ) continue ;
        size_t a = bra.p;
        size_t b = bra.q;
        size_t c = bra.r;

        
        for (size_t iket=0; iket<nkets3; iket++)
        {
           Ket3& ket = Tbc.GetKet(iket);
           if ( not (  (ket.op->cvq==0) and (ket.oq->cvq==0) and (ket.oR->cvq==0) ) ) continue; //cvq==0 means core orbit
           double d_ei = std::abs( 2*ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
           double d_ej = std::abs( 2*ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
           double d_ek = std::abs( 2*ket.oR->n + ket.oR->l - e_fermi[ket.oR->tz2]);
           double occnat_i = ket.op->occ_nat;
           double occnat_j = ket.oq->occ_nat;
           double occnat_k = ket.oR->occ_nat;
           if ( d_ei + d_ej + d_ek > modelspace->GetdE3max() ) continue;
           if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < modelspace->GetOccNat3Cut() ) continue ;
           size_t i = ket.p;
           size_t j = ket.q;
           size_t k = ket.r;

           double denominator = Get3bDenominator( a,b,c, i,j,k ) ;

           double ME_od = H->ThreeBody.GetME_pn_PN_ch(ch3,ch3,ibra,iket );
           double eta =  etafunc( ME_od, denominator);

           Eta->ThreeBody.AddToME_pn_PN_ch( ch3,ch3,ibra,iket,  eta); // hermitian conjugate automatically gets added
           
        }// for iket
      }// for ibra

    }// for ch3

    std::cout << "Norm of Eta3 = " << Eta->ThreeBodyNorm() << std::endl;
    H->profiler.timer[__func__] += omp_get_wtime() - t_start;
}




void Generator::ConstructGenerator_ShellModel(std::function<double (double,double)>& eta_func)
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : imsrg_util::VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : imsrg_util::VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
//         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(i,a) = eta_func(H->OneBody(i,a), denominator);
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
//            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(ibra,iket) = eta_func(H2(ibra,iket) , denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
//         auto& ket = tbc.GetKet(iket);
         for ( auto& ibra : imsrg_util::VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
         {
//            auto& bra = tbc.GetKet(ibra);
            double denominator = Get2bDenominator(ch,ibra,iket);
//            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator) ;
            ETA2(ibra,iket) = eta_func(H2(ibra,iket) , denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }

    if ( Eta->GetParticleRank()>2 and H->GetParticleRank()>2 and not only_2b_eta )
    {
//       ConstructGenerator_ShellModel_Atan_3body();
       ConstructGenerator_ShellModel_3body(eta_func);
    }

}




// off diagonal pieces are  ppp|ccc  ppp|ccv ppp|cvv  qpp|vvv   where p is either v or q
//void Generator::ConstructGenerator_SingleRef_3body(std::function<double (double,double)>& etafunc )
void Generator::ConstructGenerator_ShellModel_3body(std::function<double (double,double)>& etafunc )
{
     double t_start = omp_get_wtime();
     size_t ncore = modelspace->core.size();
     std::vector<size_t> corevec;
     for (auto a : modelspace->core) corevec.push_back(a);
     std::map<int,double> e_fermi = modelspace->GetEFermi();
     std::cout << __func__ << "  looping in generator 3-body part .  Size of H3 = " << H->ThreeBodyNorm() << std::endl;
//    for (auto a : modelspace->core )
     size_t nch3 = modelspace->GetNumberThreeBodyChannels();
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t ch3=0; ch3<nch3; ch3++)
    {
      ThreeBodyChannel& Tbc = modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      for (size_t ibra=0; ibra<nkets3; ibra++)
      {
        Ket3& bra = Tbc.GetKet(ibra);
//        int cvq_i = bra.op->cvq;
//        int cvq_j = bra.oq->cvq;
//        int cvq_k = bra.oR->cvq;
        if (   (bra.op->cvq==0) or (bra.oq->cvq==0) or (bra.oR->cvq==0)  ) continue; //cvq==0 means core, so we want all v or q in the bra. 
        double d_ei = std::abs( 2*bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
        double d_ej = std::abs( 2*bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
        double d_ek = std::abs( 2*bra.oR->n + bra.oR->l - e_fermi[bra.oR->tz2]);
        double occnat_i = bra.op->occ_nat;
        double occnat_j = bra.oq->occ_nat;
        double occnat_k = bra.oR->occ_nat;
        if ( d_ei + d_ej + d_ek > modelspace->GetdE3max() ) continue;
        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < modelspace->GetOccNat3Cut() ) continue ;
        size_t i = bra.p;
        size_t j = bra.q;
        size_t k = bra.r;
        
        for (size_t iket=0; iket<nkets3; iket++)
        {
           Ket3& ket = Tbc.GetKet(iket);
           if (   (ket.op->cvq==2) or (ket.oq->cvq==2) or (ket.oR->cvq==2)  ) continue; //cvq==2 means q, i.e. not core or valence. we want all c or v in ket. 
           if (  (bra.op->cvq==1) and (bra.oq->cvq==1) and (bra.oR->cvq==1) and (ket.op->cvq==1) and (ket.oq->cvq==1) and (ket.oR->cvq==1) ) continue;// no vvvvvv
           double d_ea = std::abs( 2*ket.op->n + ket.op->l - e_fermi[ket.op->tz2]);
           double d_eb = std::abs( 2*ket.oq->n + ket.oq->l - e_fermi[ket.oq->tz2]);
           double d_ec = std::abs( 2*ket.oR->n + ket.oR->l - e_fermi[ket.oR->tz2]);
           double occnat_a = ket.op->occ_nat;
           double occnat_b = ket.oq->occ_nat;
           double occnat_c = ket.oR->occ_nat;
           if ( d_ea + d_eb + d_ec > modelspace->GetdE3max() ) continue;
           if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_c*(1-occnat_c) ) < modelspace->GetOccNat3Cut() ) continue ;
           size_t a = ket.p;
           size_t b = ket.q;
           size_t c = ket.r;

           double denominator = Get3bDenominator( i,j,k, a,b,c ) ;

           double ME_od = H->ThreeBody.GetME_pn_PN_ch(ch3,ch3,ibra,iket );
           double eta =  etafunc( ME_od, denominator);

           Eta->ThreeBody.AddToME_pn_PN_ch( ch3,ch3,ibra,iket,  eta); // hermitian conjugate automatically gets added
           
        }// for iket
      }// for ibra

    }// for ch3

    std::cout << "Norm of Eta3 = " << Eta->ThreeBodyNorm() << std::endl;
    H->profiler.timer[__func__] += omp_get_wtime() - t_start;
}




/*

// off diagonal pieces are  ppp|ccc  ppp|ccv ppp|cvv  qpp|vvv   where p is either v or q
//void Generator::ConstructGenerator_ShellModel_Atan_3body()
void Generator::ConstructGenerator_ShellModel_3body(std::function<double (double,double)>& eta_func)
{
    // Three body piece
//    int E3cut = 99;
    int E3cut = modelspace->E3max;
    double t_start = omp_get_wtime();
    std::map<int,double> e_fermi = modelspace->GetEFermi();

     std::cout << "looping in generator 3-body part .  Size of H3 = " << H->ThreeBodyNorm() << std::endl;
    std::vector<size_t> alist;
    for (auto a : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) ) alist.push_back(a);
    size_t a_size = alist.size();

//    for (auto a : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace) )
    for (size_t ind_a=0; ind_a<a_size; ind_a++ )
    {
     auto a = alist[ind_a];
     Orbit& oa = modelspace->GetOrbit(a);
     double d_ea = std::abs( 2*oa.n + oa.l - e_fermi[oa.tz2]);
     for (auto b : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace ) )
     {
      if (b>a) continue;
      Orbit& ob = modelspace->GetOrbit(b);
      double d_eb = std::abs( 2*ob.n + ob.l - e_fermi[ob.tz2]);
      int Jab_min = std::abs(oa.j2-ob.j2)/2;
      int Jab_max = (oa.j2+ob.j2)/2;
      for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
      {
       for (auto c : imsrg_util::VectorUnion(modelspace->valence,modelspace->qspace ) )
       {
        if (c>b) continue;
        Orbit& oc = modelspace->GetOrbit(c);
        double d_ec = std::abs( 2*oc.n + oc.l - e_fermi[oc.tz2]);
        if ( d_ea + d_eb + d_ec > modelspace->GetdE3max() ) continue;
//        if ( (2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l) > E3cut ) continue;
        if ( (2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l) > E3cut ) continue;

        for ( auto i : imsrg_util::VectorUnion(modelspace->core,modelspace->valence) )
        {
         Orbit& oi = modelspace->GetOrbit(i);
         double d_ei = std::abs( 2*oi.n + oi.l - e_fermi[oi.tz2]);
//         // if any bra indices match any ket indices, then the normal ordering will take care of it
//         if ( i==a or i==b or i==c) continue;
         for ( auto j : imsrg_util::VectorUnion(modelspace->core,modelspace->valence) )
         {
          if (j>i) continue;
//          if ( j==a or j==b or j==c) continue;
          Orbit& oj = modelspace->GetOrbit(j);
          double d_ej = std::abs( 2*oj.n + oj.l - e_fermi[oj.tz2]);
          int Jij_min = std::abs(oi.j2-oj.j2)/2;
          int Jij_max = (oi.j2+oj.j2)/2;
          for (int Jij=Jij_min; Jij<=Jij_max; Jij++)
          {
           for ( auto k : imsrg_util::VectorUnion(modelspace->core,modelspace->valence) )
           {
            if (k>j) continue;
//            if ( k==a or k==b or k==c) continue;
            // check if everything is valence, becase we don't want vvv|vvv to be off-diagonal
            if (   (std::find( modelspace->valence.begin(), modelspace->valence.end(), a) != modelspace->valence.end()) // a is valence
              and  (std::find( modelspace->valence.begin(), modelspace->valence.end(), b) != modelspace->valence.end()) // b is valence
              and  (std::find( modelspace->valence.begin(), modelspace->valence.end(), c) != modelspace->valence.end()) // c is valence
              and  (std::find( modelspace->valence.begin(), modelspace->valence.end(), i) != modelspace->valence.end()) // i is valence
              and  (std::find( modelspace->valence.begin(), modelspace->valence.end(), j) != modelspace->valence.end()) // j is valence
              and  (std::find( modelspace->valence.begin(), modelspace->valence.end(), k) != modelspace->valence.end()) // k is valence
              ) continue;

            Orbit& ok = modelspace->GetOrbit(k);
            double d_ek = std::abs( 2*ok.n + ok.l - e_fermi[ok.tz2]);
            if ( d_ei + d_ej + d_ek > modelspace->GetdE3max() ) continue;

//            if ( (2*(oi.n+oj.n+ok.n)+oi.l+oj.l+ok.l) > 2*E3cut ) continue;
            if ( (2*(oi.n+oj.n+ok.n)+oi.l+oj.l+ok.l) > E3cut ) continue;
            if ( (oa.l+ob.l+oc.l+oi.l+oj.l+ok.l)%2>0 ) continue;
            if ( (oa.tz2+ob.tz2+oc.tz2) != (oi.tz2+oj.tz2+ok.tz2) ) continue;
            double denominator = Get3bDenominator( a,b,c, i,j,k ) ;

            int twoJ_min = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jij-ok.j2) );
            int twoJ_max = std::min( 2*Jab+oc.j2, 2*Jij+ok.j2 );
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
            {
              double ME_od = H->ThreeBody.GetME_pn( Jab, Jij, twoJ, a,b,c,i,j,k);
              double eta = 0.5*atan(2*ME_od / denominator);
//              double eta = (ME_od / denominator);
//              if (std::abs(eta)>1e-6)
//              if (std::abs(eta)>1e-3)
//              {
//                std::cout << "big eta !  <" << a << " " << b << " " << c << ", " << Jab << " | " << i << " " << j << " " << k << ", " << Jij << " > _ " << twoJ << "  => " << ME_od << " / " << denominator << " = " << eta << std::endl;
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
    }// for a


}

*/

//void Generator::ConstructGenerator_ShellModel_Atan_NpNh()
void Generator::ConstructGenerator_ShellModel_NpNh(std::function<double(double,double)>& eta_func)
{
//  ConstructGenerator_ShellModel_Atan();
  ConstructGenerator_ShellModel(atan_func);

//  std::cout << "In ShellModel_Atat_NpNh, adding to Eta" << std::endl;
  // decouple f_cc'
  for ( auto& c : modelspace->core )
  {
   for ( auto& cprime : modelspace->core )
   {
     if (cprime<=c) continue;
     double denominator = Get1bDenominator(c,cprime);
//     Eta->OneBody(c,cprime) = 0.5*atan(2*H->OneBody(c,cprime) / denominator );
     Eta->OneBody(c,cprime) = eta_func(H->OneBody(c,cprime), denominator );
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
//         ETA2(ibra,iket) = 0.5*atan( 2*H2(ibra,iket) / denominator );
         ETA2(ibra,iket) = eta_func( H2(ibra,iket) , denominator );
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
//         ETA2(ibra,iket) = 0.5*atan( 2*H2(ibra,iket) / denominator );
         ETA2(ibra,iket) = eta_func( H2(ibra,iket) , denominator );
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
   for (auto i : modelspace->all_orbits)
   {
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
//void Generator::ConstructGenerator_1PA()
void Generator::ConstructGenerator_1PA(std::function<double(double,double)>& eta_func)
{

   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : imsrg_util::VectorUnion(modelspace->core, modelspace->valence))
   {
      for (auto& i : imsrg_util::VectorUnion( modelspace->valence, modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
//         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
         Eta->OneBody(i,a) = eta_func(H->OneBody(i,a),denominator);
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
//            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(ibra,iket) = eta_func(H2(ibra,iket) , denominator);
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




