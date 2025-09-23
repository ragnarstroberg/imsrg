#include <random>
#include "Generator.hh"
#include "Commutator.hh"
#include "Operator.hh"
#include "PhysicalConstants.hh" // for HBARC and M_NUCLEON

#include "omp.h"
#include <string>

using PhysConst::M_NUCLEON;
using PhysConst::HBARC;

std::function<double(double,double)> Generator::wegner_func = [] (double Hod, double denom){ return Hod * denom;};
std::function<double(double,double)> Generator::white_func = [] (double Hod, double denom){ return Hod / denom;};
std::function<double(double,double)> Generator::atan_func = [] (double Hod, double denom){ return 0.5 * atan(2*Hod / denom);};
std::function<double(double,double)> Generator::imaginarytime_func = [] (double Hod, double denom){ return Hod * ( (denom>0)-(denom<0) ) ;};//signum function

std::function<double(double,double)> Generator::qtransferatan1_func = [](double Hod, double denom){return pow(std::abs(denom)*M_NUCLEON/HBARC/HBARC, 0.5*1) * atan_func(Hod, denom);};

Generator::Generator()
  : generator_type("white"),/* modelspace(NULL),*/ denominator_cutoff(1e-6)  , denominator_delta(0), denominator_delta_index(-1), denominator_partitioning(Epstein_Nesbet),  only_2b_eta(false), use_isospin_averaging(false)
{}



void Generator::SetDenominatorPartitioning(std::string dp)
{
   if (dp=="Moller_Plesset" or dp=="Moller-Plesset" or dp=="MollerPlesset" or dp=="MP") denominator_partitioning=Moller_Plesset;
   else if ( dp=="MP_isospin") denominator_partitioning=MP_isospin;
   else denominator_partitioning=Epstein_Nesbet;
}

void Generator::Update(Operator& H_s, Operator& Eta_s)
{
   Eta_s.Erase();
   AddToEta(H_s,Eta_s);
   if (use_isospin_averaging)
   {
//      Eta_s = Eta_s.DoIsospinAveraging();
      Eta_s = Eta_s.UndoNormalOrdering().DoIsospinAveraging().DoNormalOrdering();
   }
}


void Generator::AddToEta(Operator& H_s, Operator& Eta_s)
{
   double start_time = omp_get_wtime();
   H = &H_s;
   Eta = &Eta_s;
//   modelspace = H->GetModelSpace();

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
     ConstructGenerator_SingleRef( qtransferatanN_func );
   }
   else
   {
      std::cout << "Error. Unkown generator_type: " << generator_type << std::endl;
   }
   Eta->profiler.timer["UpdateEta"] += omp_get_wtime() - start_time;

}



Operator Generator::GetHod(Operator& H)
{
   for (auto sr : {"wegner","white","atan","imaginary-time","qtransfer-atan"})
   {
      if (generator_type == sr )  return GetHod_SingleRef(H);
   }
   for (auto sm : {"shell-model-wegner","shell-model","shell-model-atan","shell-model-imaginary-time"})
   {
      if (generator_type == sm )  return GetHod_ShellModel(H);
   }
   std::cout << "GetHod not implemented for generator type " << generator_type << "   so you get zero." << std::endl;
   return 0*H; 
}



Operator Generator::GetHod_CC(Operator& H, std:: string lr)
{
   if(lr == "left"){
   return GetHod_SingleRef_left (H);
   }
   else{
   return GetHod_SingleRef_right (H);
   }
}


// Old method used to test some things out. Not typically used.
void Generator::SetDenominatorDeltaOrbit(std::string orb)
{
  if (orb == "all")
     SetDenominatorDeltaIndex(-12345);
  else
  {
     SetDenominatorDeltaIndex( H->modelspace->GetOrbitIndex(orb) );
     std::cout << "Setting denominator delta orbit " << orb << " => " << H->modelspace->GetOrbitIndex(orb) << std::endl;
  }
}


// Epstein-Nesbet energy denominators for White-type generator_types
double Generator::Get1bDenominator(int i, int j) 
{
   double ni = H->modelspace->GetOrbit(i).occ;
   double nj = H->modelspace->GetOrbit(j).occ;
   
   double denominator = H->OneBody(i,i) - H->OneBody(j,j);
   if ( denominator_partitioning == Epstein_Nesbet)
   {
      denominator += ( ni-nj ) * H->TwoBody.GetTBMEmonopole(i,j,i,j);
   }
   else if ( denominator_partitioning == MP_isospin )
   {
      // Average the SP energies of protons and neutrons to make the denominator isospin symmetric
      Orbit& oi = H->modelspace->GetOrbit(i);
      Orbit& oj = H->modelspace->GetOrbit(j);
      size_t ii = H->modelspace->GetOrbitIndex( oi.n,oi.l,oi.j2,-oi.tz2);
      size_t jj = H->modelspace->GetOrbitIndex( oj.n,oj.l,oj.j2,-oj.tz2);
      denominator += H->OneBody(ii,ii) - H->OneBody(jj,jj);
      denominator /=2;
   }

   if (denominator_delta_index==-12345 or i == denominator_delta_index or j==denominator_delta_index)
     denominator += denominator_delta;

   // if (std::abs(denominator) < denominator_cutoff)
   // {
   //    if (denominator > 0)
   //    {
   //       denominator = denominator_cutoff;
   //    }
   //    else
   //    {
   //       denominator = - denominator_cutoff;
   //    }    
   // }
   if (std::abs(denominator)<denominator_cutoff)
     denominator = denominator_cutoff;

   return denominator;
}


double Generator::Get2bDenominator(int ch_bra, int ch_ket, int ibra, int iket) 
{
   TwoBodyChannel& tbc_bra = H->modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket = H->modelspace->GetTwoBodyChannel(ch_ket);
   Ket & bra = tbc_bra.GetKet(ibra);
   Ket & ket = tbc_ket.GetKet(iket);
   int i = bra.p;
   int j = bra.q;
   int k = ket.p;
   int l = ket.q;
   double denominator = H->OneBody(i,i)+ H->OneBody(j,j) - H->OneBody(k,k) - H->OneBody(l,l);

   if ( denominator_partitioning == MP_isospin )
   {
      // Average the SP energies of protons and neutrons to make the denominator isospin symmetric
      Orbit& oi = H->modelspace->GetOrbit(i);
      Orbit& oj = H->modelspace->GetOrbit(j);
      Orbit& ok = H->modelspace->GetOrbit(k);
      Orbit& ol = H->modelspace->GetOrbit(l);
      size_t ii = H->modelspace->GetOrbitIndex( oi.n,oi.l,oi.j2,-oi.tz2);
      size_t jj = H->modelspace->GetOrbitIndex( oj.n,oj.l,oj.j2,-oj.tz2);
      size_t kk = H->modelspace->GetOrbitIndex( ok.n,ok.l,ok.j2,-ok.tz2);
      size_t ll = H->modelspace->GetOrbitIndex( ol.n,ol.l,ol.j2,-ol.tz2);
      denominator += H->OneBody(ii,ii) + H->OneBody(jj,jj) - H->OneBody(kk,kk) - H->OneBody(ll,ll);
      denominator /=2;
   }

   if (denominator_delta_index == -12345) denominator += denominator_delta;
   double ni = bra.op->occ;
   double nj = bra.oq->occ;
   double nk = ket.op->occ;
   double nl = ket.oq->occ;

   if ( denominator_partitioning == Epstein_Nesbet)
   {
     denominator       += ( 1-ni-nj ) * H->TwoBody.GetTBMEmonopole(i,j,i,j); // pp'pp'
     denominator       -= ( 1-nk-nl ) * H->TwoBody.GetTBMEmonopole(k,l,k,l); // hh'hh'
     denominator       += ( ni-nk )   * H->TwoBody.GetTBMEmonopole(i,k,i,k); // phph
     denominator       += ( ni-nl )   * H->TwoBody.GetTBMEmonopole(i,l,i,l); // ph'ph'
     denominator       += ( nj-nk )   * H->TwoBody.GetTBMEmonopole(j,k,j,k); // p'hp'h
     denominator       += ( nj-nl )   * H->TwoBody.GetTBMEmonopole(j,l,j,l); // p'h'p'h'
   }

   // if (std::abs(denominator) < denominator_cutoff)
   // {
   //   if (denominator > 0)
   //   {
   //       denominator = denominator_cutoff;
   //   }
   //   else
   //   {
   //       denominator = -denominator_cutoff;
   //   }
   // }
   if (std::abs(denominator)<denominator_cutoff)
     denominator = denominator_cutoff;

   return denominator;
}

double Generator::Get3bDenominator( int i, int j, int k, int l, int m, int n )
{
  auto h1 = H->OneBody;
  double denominator = h1(i,i) + h1(j,j) + h1(k,k) - h1(l,l) - h1(m,m) - h1(n,n);
  if ( denominator_partitioning == Epstein_Nesbet)
  { // This  is unpleasant and almost certainly not necessary, so not implementing it yet...
  }
  return denominator;
}

// Keep the Jdependence for the Gamma_ijij and Gamma_klkl terms, because it's
// relatively unambiguous to work out
double Generator::Get2bDenominator_Jdep(int ch, int ibra, int iket) 
{
   TwoBodyChannel& tbc = H->modelspace->GetTwoBodyChannel(ch);
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
   return denominator;
}







void Generator::ConstructGenerator_SingleRef(std::function<double (double,double)>& etafunc )
{
   // One body piece -- eliminate ph bits
   for ( auto& a : H->modelspace->core)
   {
      for ( auto& i : VectorUnion(H->modelspace->valence, H->modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = etafunc( H->OneBody(i,a), denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for ( auto& iter : Eta->TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H->modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H->modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat& ETA2 =  iter.second;
      arma::mat& H2 = H->TwoBody.GetMatrix(ch_bra,ch_ket);
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
            double denominator = Get2bDenominator(ch_bra,ch_ket,ibra,iket);
            ETA2(ibra,iket) = etafunc( H2(ibra,iket), denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }

    if ( Eta->GetParticleRank()>2 and H->GetParticleRank()>2 and not only_2b_eta )
    {
       double t_start = omp_get_wtime();
       ConstructGenerator_SingleRef_3body( etafunc );
       H->profiler.timer["Update Eta 3body"] += omp_get_wtime() - t_start;
    }// if particle rank >3
}


 


// Off-diagonal pieces are <abc|ijk> = <ppp|ccc> where c is core and p is either valence or q (that is, not core).
//void Generator::ConstructGenerator_Atan_3body()
void Generator::ConstructGenerator_SingleRef_3body(std::function<double (double,double)>& etafunc )
{
     double t_start = omp_get_wtime();
//     size_t ncore = H->modelspace->core.size();
     std::vector<size_t> corevec;
     for (auto a : H->modelspace->core) corevec.push_back(a);
     std::map<int,double> e_fermi = H->modelspace->GetEFermi();
//     std::cout << __func__ << "  looping in generator 3-body part .  Size of H3 = " << H->ThreeBodyNorm() << std::endl;
//    for (auto a : H->modelspace->core )
     size_t nch3 = H->modelspace->GetNumberThreeBodyChannels();
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t ch3=0; ch3<nch3; ch3++)
    {
      ThreeBodyChannel& Tbc = H->modelspace->GetThreeBodyChannel(ch3);
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
        if ( d_ea + d_eb + d_ec > H->modelspace->GetdE3max() ) continue;
        if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_c*(1-occnat_c) ) < H->modelspace->GetOccNat3Cut() ) continue ;
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
           if ( d_ei + d_ej + d_ek > H->modelspace->GetdE3max() ) continue;
           if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < H->modelspace->GetOccNat3Cut() ) continue ;
           size_t i = ket.p;
           size_t j = ket.q;
           size_t k = ket.r;

           double denominator = Get3bDenominator( a,b,c, i,j,k ) ;

           double ME_od = H->ThreeBody.GetME_pn_ch(ch3,ch3,ibra,iket );
           double eta =  etafunc( ME_od, denominator);

           Eta->ThreeBody.AddToME_pn_ch( ch3,ch3,ibra,iket,  eta); // hermitian conjugate automatically gets added
           
        }// for iket
      }// for ibra

    }// for ch3

//    std::cout << "Norm of Eta3 = " << std::setprecision(8) <<  Eta->ThreeBodyNorm() << std::endl;
    H->profiler.timer[__func__] += omp_get_wtime() - t_start;
}



 
void Generator::ConstructGenerator_ShellModel(std::function<double (double,double)>& eta_func)
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : VectorUnion( H->modelspace->core, H->modelspace->valence))
   {
      for (auto& i : VectorUnion( H->modelspace->valence, H->modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = eta_func(H->OneBody(i,a), denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }


   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = H->modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = H->modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H->TwoBody.GetMatrix(ch);

      // Decouple the core
      for ( auto& iket : VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = eta_func(H2(ibra,iket) , denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }

      // Decouple the valence space
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) ) 
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = eta_func(H2(ibra,iket) , denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }

    if ( Eta->GetParticleRank()>2 and H->GetParticleRank()>2 and not only_2b_eta )
    {
       ConstructGenerator_ShellModel_3body(eta_func);
    }

}





// off diagonal pieces are  ppp|ccc  ppp|ccv ppp|cvv  qpp|vvv   where p is either v or q
void Generator::ConstructGenerator_ShellModel_3body(std::function<double (double,double)>& etafunc )
{
     double t_start = omp_get_wtime();
     std::vector<size_t> corevec;
     for (auto a : H->modelspace->core) corevec.push_back(a);
     std::map<int,double> e_fermi = H->modelspace->GetEFermi();
//     std::cout << __func__ << "  looping in generator 3-body part .  Size of H3 = " << H->ThreeBodyNorm() << std::endl;
     size_t nch3 = H->modelspace->GetNumberThreeBodyChannels();
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t ch3=0; ch3<nch3; ch3++)
    {
      ThreeBodyChannel& Tbc = H->modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      for (size_t ibra=0; ibra<nkets3; ibra++)
      {
        Ket3& bra = Tbc.GetKet(ibra);
        if (   (bra.op->cvq==0) or (bra.oq->cvq==0) or (bra.oR->cvq==0)  ) continue; //cvq==0 means core, so we want all v or q in the bra. 
        double d_ei = std::abs( 2*bra.op->n + bra.op->l - e_fermi[bra.op->tz2]);
        double d_ej = std::abs( 2*bra.oq->n + bra.oq->l - e_fermi[bra.oq->tz2]);
        double d_ek = std::abs( 2*bra.oR->n + bra.oR->l - e_fermi[bra.oR->tz2]);
        double occnat_i = bra.op->occ_nat;
        double occnat_j = bra.oq->occ_nat;
        double occnat_k = bra.oR->occ_nat;
        if ( d_ei + d_ej + d_ek > H->modelspace->GetdE3max() ) continue;
        if ( (occnat_i*(1-occnat_i) * occnat_j*(1-occnat_j) * occnat_k*(1-occnat_k) ) < H->modelspace->GetOccNat3Cut() ) continue ;
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
           if ( d_ea + d_eb + d_ec > H->modelspace->GetdE3max() ) continue;
           if ( (occnat_a*(1-occnat_a) * occnat_b*(1-occnat_b) * occnat_c*(1-occnat_c) ) < H->modelspace->GetOccNat3Cut() ) continue ;
           size_t a = ket.p;
           size_t b = ket.q;
           size_t c = ket.r;

           double denominator = Get3bDenominator( i,j,k, a,b,c ) ;

           double ME_od = H->ThreeBody.GetME_pn_ch(ch3,ch3,ibra,iket );
           double eta =  etafunc( ME_od, denominator);

           Eta->ThreeBody.AddToME_pn_ch( ch3,ch3,ibra,iket,  eta); // hermitian conjugate automatically gets added
           
        }// for iket
      }// for ibra

    }// for ch3

//    std::cout << "Norm of Eta3 = " << Eta->ThreeBodyNorm() << std::endl;
    H->profiler.timer[__func__] += omp_get_wtime() - t_start;
}




void Generator::ConstructGenerator_ShellModel_NpNh(std::function<double(double,double)>& eta_func)
{
  ConstructGenerator_ShellModel(atan_func);

  for ( auto& c : H->modelspace->core )
  {
   for ( auto& cprime : H->modelspace->core )
   {
     if (cprime<=c) continue;
     double denominator = Get1bDenominator(c,cprime);
     Eta->OneBody(c,cprime) = eta_func(H->OneBody(c,cprime), denominator );
     Eta->OneBody(cprime,c) = - Eta->OneBody(c,cprime);
   }
  }

  int nchan = H->modelspace->GetNumberTwoBodyChannels();
  for (int ch=0;ch<nchan;++ch)
  {
     TwoBodyChannel& tbc = H->modelspace->GetTwoBodyChannel(ch);
     arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
     arma::mat& H2 =  H->TwoBody.GetMatrix(ch);
  // decouple Gamma_qcvc'
     for (auto& iket : tbc.GetKetIndex_vc())
     {
       Ket& ket = tbc.GetKet(iket);
       for (auto& ibra : tbc.GetKetIndex_qc())
       {
         Ket& bra = tbc.GetKet(ibra);
         if ((ket.p==bra.p) or (ket.p==bra.q) or (ket.q==bra.p) or (ket.q==bra.q) ) continue;
         double denominator = Get2bDenominator(ch,ibra,iket);
         ETA2(ibra,iket) = eta_func( H2(ibra,iket) , denominator );
         ETA2(iket,ibra) = -ETA2(ibra,iket);
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
         ETA2(ibra,iket) = eta_func( H2(ibra,iket) , denominator );
         ETA2(iket,ibra) = -ETA2(ibra,iket);
       }
     }
  }
}


void Generator::ConstructGenerator_HartreeFock()
{
   Eta->SetParticleRank(1);
   // One body piece -- eliminate ph bits
   for (auto i : H->modelspace->all_orbits)
   {
      for (auto j : H->modelspace->all_orbits)
      {
         if (j>i) continue;
         double denominator = Get1bDenominator(i,j);
         Eta->OneBody(i,j) = H->OneBody(i,j)/denominator;
         Eta->OneBody(j,i) = - Eta->OneBody(i,j);
      }
   } 
}




// So far this is useless
void Generator::ConstructGenerator_1PA(std::function<double(double,double)>& eta_func)
{

   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : VectorUnion(H->modelspace->core, H->modelspace->valence))
   {
      for (auto& i : VectorUnion( H->modelspace->valence, H->modelspace->qspace ) )
      {
         if (i==a) continue;
         double denominator = Get1bDenominator(i,a);
         Eta->OneBody(i,a) = eta_func(H->OneBody(i,a),denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
      }
   }


   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = H->modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = H->modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H->TwoBody.GetMatrix(ch);

      // Decouple the core
      for ( auto& iket : VectorUnion( tbc.GetKetIndex_cc(), tbc.GetKetIndex_vc() ) )
      {
         for ( auto& ibra : VectorUnion( tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv(), tbc.GetKetIndex_qq() ) )
         {
            double denominator = Get2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) = eta_func(H2(ibra,iket) , denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }
    }

}






Operator  Generator::GetHod_SingleRef(Operator& H )
{
   Operator Hod = 0.0* H;
   // One body piece -- eliminate ph bits
   for ( auto& a : H.modelspace->core)
   {
      for ( auto& i : VectorUnion(H.modelspace->valence,H.modelspace->qspace) )
      {
         Hod.OneBody(i,a) = H.OneBody(i,a);
         Hod.OneBody(a,i) = H.OneBody(a,i);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat& H2 =  iter.second;
      arma::mat& Hod2 = Hod.TwoBody.GetMatrix(ch_bra,ch_ket);
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
            Hod2(ibra,iket) =  H2(ibra,iket);
            Hod2(iket,ibra) =  H2(iket,ibra) ; // Eta needs to be antisymmetric
         }
      }
    }

   // Three body --- ppp hhh bits.
    size_t nch3 = H.modelspace->GetNumberThreeBodyChannels();
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t ch3=0; ch3<nch3; ch3++)
    {
      ThreeBodyChannel& Tbc = H.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      for (size_t ibra=0; ibra<nkets3; ibra++)
      {
        Ket3& bra = Tbc.GetKet(ibra);
        // bra should be ppp where p is eiher v or q
        if ( (  (bra.op->cvq==0) or (bra.oq->cvq==0) or (bra.oR->cvq==0) ) ) continue; //cvq==0 means core orbit
        
        for (size_t iket=0; iket<nkets3; iket++)
        {
           Ket3& ket = Tbc.GetKet(iket);
           // ket should be ccc
           if ( not (  (ket.op->cvq==0) and (ket.oq->cvq==0) and (ket.oR->cvq==0) ) ) continue; //cvq==0 means core orbit

           double h_abcijk = H.ThreeBody.GetME_pn_ch(ch3,ch3,ibra,iket );

           Hod.ThreeBody.SetME_pn_ch( ch3,ch3,ibra,iket,  h_abcijk); // hermitian conjugate automatically gets added
           
        }// for iket
      }// for ibra

    }// for ch3

   return Hod;
}





 
Operator Generator::GetHod_ShellModel(Operator& H)
{
   Operator Hod = 0.0* H;
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& a : VectorUnion( H.modelspace->core, H.modelspace->valence))
   {
      for (auto& i : VectorUnion( H.modelspace->valence, H.modelspace->qspace ) )
      {
         if (i==a) continue;
         Hod.OneBody(i,a) = H.OneBody(i,a);
         Hod.OneBody(a,i) = H.OneBody(a,i);
      }
   }


   // Two body piece -- eliminate pp'hh' bits
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat& H2 =  iter.second;
      arma::mat& Hod2 = Hod.TwoBody.GetMatrix(ch_bra,ch_ket);

      // Decouple the core
      for ( auto& iket : VectorUnion( tbc_ket.GetKetIndex_cc(), tbc_ket.GetKetIndex_vc() ) ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra :  VectorUnion( tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv(), tbc_bra.GetKetIndex_qq() ) )
         {
            Hod2(ibra,iket) =  H2(ibra,iket);
            Hod2(iket,ibra) =  H2(iket,ibra) ; // Eta needs to be antisymmetric
         }
      }

      // Decouple the valence space
      for ( auto& iket :  tbc_ket.GetKetIndex_vv() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra :  VectorUnion( tbc_bra.GetKetIndex_qv(), tbc_bra.GetKetIndex_qq() ) )
         {
            Hod2(ibra,iket) =  H2(ibra,iket);
            Hod2(iket,ibra) =  H2(iket,ibra) ; // Eta needs to be antisymmetric
         }
      }

    }


 
    // off-diagonal:   <ppp|ccc>, <ppp|vcc>, <ppp|vvc>, <qpp|vvv>  where p is v or q
    //                 
    size_t nch3 = H.modelspace->GetNumberThreeBodyChannels();
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t ch3=0; ch3<nch3; ch3++)
    {
      ThreeBodyChannel& Tbc = H.modelspace->GetThreeBodyChannel(ch3);
      size_t nkets3 = Tbc.GetNumberKets();
      for (size_t ibra=0; ibra<nkets3; ibra++)
      {
        Ket3& bra = Tbc.GetKet(ibra);
        if (   (bra.op->cvq==0) or (bra.oq->cvq==0) or (bra.oR->cvq==0)  ) continue; //cvq==0 means core, so we want all v or q in the bra. 

        
        for (size_t iket=0; iket<nkets3; iket++)
        {
           Ket3& ket = Tbc.GetKet(iket);
           if (   (ket.op->cvq==2) or (ket.oq->cvq==2) or (ket.oR->cvq==2)  ) continue; //cvq==2 means q, i.e. not core or valence. we want all c or v in ket. 
           if (  (bra.op->cvq==1) and (bra.oq->cvq==1) and (bra.oR->cvq==1) and (ket.op->cvq==1) and (ket.oq->cvq==1) and (ket.oR->cvq==1) ) continue;// no vvvvvv


           double ME_od = H.ThreeBody.GetME_pn_ch(ch3,ch3,ibra,iket );

           Hod.ThreeBody.SetME_pn_ch( ch3,ch3,ibra,iket,  ME_od); // hermitian conjugate automatically gets added

           
        }// for iket
      }// for ibra

    }// for ch3



    return Hod;
}
 



Operator  Generator::GetHod_SingleRef_left(Operator& H )
{
   Operator Hod = 0.0* H;
   // One body piece -- eliminate ph bits
   for ( auto& i : H.modelspace->core)
   {
      for ( auto& a : VectorUnion(H.modelspace->valence,H.modelspace->qspace) )
      {
          // left is excitation
         //Hod.OneBody(i,a) = H.OneBody(i,a);
         Hod.OneBody(a,i) = H.OneBody(a,i);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat& H2 =  iter.second;
      arma::mat& Hod2 = Hod.TwoBody.GetMatrix(ch_bra,ch_ket);
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
            Hod2(ibra,iket) =  H2(ibra,iket);
            //Hod2(iket,ibra) =  H2(iket,ibra) ; // Eta needs to be antisymmetric
         }
      }
    }


   return Hod;
}



Operator  Generator::GetHod_SingleRef_right(Operator& H )
{
   Operator Hod = 0.0* H;
   // One body piece -- eliminate ph bits
   for ( auto& a : H.modelspace->core)
   {
      for ( auto& i : VectorUnion(H.modelspace->valence,H.modelspace->qspace) )
      {
         Hod.OneBody(i,a) = H.OneBody(i,a);
        // Hod.OneBody(a,i) = H.OneBody(a,i);
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat& H2 =  iter.second;
      arma::mat& Hod2 = Hod.TwoBody.GetMatrix(ch_bra,ch_ket);
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
       //     Hod2(ibra,iket) =  H2(ibra,iket);
            Hod2(iket,ibra) =  H2(iket,ibra) ; // Eta needs to be antisymmetric
         }
      }
    }


   return Hod;
}


Operator  Generator::Get_adjoint(Operator& H )
{
   Operator Hod = 0.0* H;
   // One body piece -- eliminate ph bits
   for ( auto& i : H.modelspace->core)
   {
	
        Orbit &oi = H.modelspace->GetOrbit(i);
      for ( auto& a : VectorUnion(H.modelspace->valence,H.modelspace->qspace) )
      {
        Orbit &oa = H.modelspace->GetOrbit(a);

	int phase_ia = H.modelspace->phase((oi.j2-oa.j2)/2);

         Hod.OneBody(i,a) = H.OneBody(a,i)*phase_ia;
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat& H2 =  iter.second;
      arma::mat& Hod2 = Hod.TwoBody.GetMatrix(ch_ket,ch_bra);
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
       //     Hod2(ibra,iket) =  H2(ibra,iket);
	int phase_jj=H.modelspace->phase(tbc_bra.J-tbc_ket.J);

            Hod2(iket,ibra) =  H2(ibra,iket)*phase_jj ;
         }
      }
    }


   return Hod;
}








double  Generator::GetOverlap(Operator& H1, Operator& H2 )
{
   double ovlp = 0;
   // One body piece -- eliminate ph bits
   for ( auto& i : H1.modelspace->holes)
   {

      Orbit &oi = H1.modelspace->GetOrbit(i);

      for ( auto& a : VectorUnion(H1.modelspace->valence,H1.modelspace->qspace) )
      {
         ovlp += H1.OneBody(a,i) * H2.OneBody(a,i)* (oi.j2+1);
      }
   }

   for ( auto& iter : H1.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];

      TwoBodyChannel& tbc_bra = H1.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H1.modelspace->GetTwoBodyChannel(ch_ket);

      arma::mat& Hmat1 =  iter.second;
      arma::mat& Hmat2 = H2.TwoBody.GetMatrix(ch_bra,ch_ket);

      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
             ovlp += Hmat1(ibra,iket)*Hmat2(ibra,iket)*(2*tbc_bra.J +1);
         }
      }
    }
   return ovlp;

}







Operator  Generator::GetEOM_ladder(Operator& H , int herm)
{

  // Generate a (anti-)hermit operator, regardless the hermitian of operator of H
  // int herm, 0 for hermit, and 1 for antihermit
    int hZ = H.IsHermitian() ? +1 : -1;

    int herm_phase =0;
    Operator Hod = 0.0* H;

  if(herm == 0){
   Hod.SetHermitian();
     herm_phase =1;
}
  if(herm == 1){
   Hod.SetAntiHermitian();
     herm_phase = -1 ;
}
   

   // One body piece -- eliminate ph bits
   // 
   for ( auto& i : H.modelspace->core)
   {
      Orbit& oi = H.modelspace->GetOrbit(i);
      for ( auto& a : VectorUnion(H.modelspace->valence,H.modelspace->qspace) )
      {
      Orbit& oa = H.modelspace->GetOrbit(a);
         Hod.OneBody(a,i) = H.OneBody(a,i);
	int phase_ia = H.modelspace->phase((oa.j2-oi.j2)/2);
         Hod.OneBody(i,a) = H.OneBody(a,i)*phase_ia*herm_phase;
	
      }
   }


//Hod.EraseTwoBody();

   // Two body piece only stored half channel, no need to change
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);

      arma::mat& H2 =  iter.second;

      arma::mat& Hod2 = Hod.TwoBody.GetMatrix(ch_bra,ch_ket);


// diagonal channel
	// diagnal channel, add <ab|ij> and <ij||ab> will be automatically added with proper factor from (anti)hermitian
if(ch_bra==ch_ket){

      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
	{
	Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));
}
}



}
 
 if(ch_bra != ch_ket){


// off-diagonal channel
   

      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
	{
        //Hod2(ibra,iket)= H2(ibra,iket);
	Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));
}
}


      for ( auto& ibra : tbc_bra.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& iket : VectorUnion(tbc_ket.GetKetIndex_qq(), tbc_ket.GetKetIndex_vv(), tbc_ket.GetKetIndex_qv() ) )
	{
        //Hod2(ibra,iket)= H2(ibra,iket)*hZ*herm_phase;
	Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket)*hZ*herm_phase);
}
}

}


}


   return Hod;
}



double  Generator::GetEOM_Overlap(Operator& H1, Operator& H2 )
{
   double ovlp = 0;
   // One body piece -- eliminate ph bits
   for ( auto& i : H1.modelspace->holes)
   {

      Orbit &oi = H1.modelspace->GetOrbit(i);

      for ( auto& a : VectorUnion(H1.modelspace->valence,H1.modelspace->qspace) )
      {
         ovlp += H1.OneBody(a,i) * H2.OneBody(a,i);
      }
   }

   for ( auto& iter : H1.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];

      TwoBodyChannel& tbc_bra = H1.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H1.modelspace->GetTwoBodyChannel(ch_ket);

      arma::mat& Hmat1 =  iter.second;
      arma::mat& Hmat2 = H2.TwoBody.GetMatrix(ch_bra,ch_ket);

// diagonal channel
     if(ch_bra==ch_ket){
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
             ovlp += Hmat1(ibra,iket)*Hmat2(ibra,iket);
         }
      }
}

// off-diagonal channel
     if(ch_bra!=ch_ket){

      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
             ovlp += Hmat1(ibra,iket)*Hmat2(ibra,iket);
         }
      }

      for ( auto& iket : VectorUnion(tbc_ket.GetKetIndex_qq(), tbc_ket.GetKetIndex_vv(), tbc_ket.GetKetIndex_qv() ) )
      {
         for ( auto& ibra : tbc_bra.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
         {
             ovlp += Hmat1(ibra,iket)*Hmat2(ibra,iket);
         }
      }


}


    }
   return ovlp;

}


void  Generator::force_decouple(Operator& H )
{

   

   for ( auto& i : H.modelspace->valence)
   {
      Orbit& oi = H.modelspace->GetOrbit(i);
      for ( auto& a : H.modelspace->valence )
      {
      Orbit& oa = H.modelspace->GetOrbit(a);
//
      //   H.OneBody(a,i) = 0.;
//
      //   H.OneBody(i,a) = 0.;
      }
   }

   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);



// vvvv
//
      for ( auto& iket : tbc_ket.GetKetIndex_vv() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : tbc_bra.GetKetIndex_vv()  )
	{
      //  H.TwoBody.SetTBME(ch_bra,ch_ket, ibra, iket, 0.);
}
}
}



   // One body piece -- eliminate ph bits
   //  particle hole excitation
   //
   for ( auto& i : H.modelspace->core)
   {
      Orbit& oi = H.modelspace->GetOrbit(i);
      for ( auto& a : VectorUnion(H.modelspace->valence,H.modelspace->qspace) )
      {
      Orbit& oa = H.modelspace->GetOrbit(a);
         H.OneBody(a,i) = 0.;
         H.OneBody(i,a) = 0.;
      }
   }

   for ( auto& i : H.modelspace->valence)
   {
      Orbit& oi = H.modelspace->GetOrbit(i);
      for ( auto& a : H.modelspace->qspace )
      {
      Orbit& oa = H.modelspace->GetOrbit(a);

         H.OneBody(a,i) = 0.;

         H.OneBody(i,a) = 0.;
      }
   }






//Hod.EraseTwoBody();

   // Two body piece only stored half channel, no need to change
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);



// PPHH
//
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
	{
        H.TwoBody.SetTBME(ch_bra,ch_ket, ibra, iket, 0.);
}
}


// PPvq
//
      for ( auto& iket : tbc_ket.GetKetIndex_vc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
	{
	H.TwoBody.SetTBME(ch_bra,ch_ket, ibra, iket, 0.);
}
}


// PPvv
//
      for ( auto& iket : tbc_ket.GetKetIndex_vv() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv() ) )
	{
      H.TwoBody.SetTBME(ch_bra,ch_ket, ibra, iket, 0.);
}
}

}

}





Operator  Generator::GetVSEOM_ladder(Operator& H , int herm)
{

  // Generate a (anti-)hermit operator, regardless the hermitian of operator of H
  // int herm, 0 for hermit, and 1 for antihermit
    Operator Hod = 0.0* H;
    int herm_phase=1;

  if(herm == 0){
   Hod.SetHermitian();
     herm_phase =1;
}
  if(herm == 1){
   Hod.SetAntiHermitian();
     herm_phase = -1 ;
}
   

   // One body piece -- eliminate ph bits
   //  particle hole excitation
   //
   for ( auto& i : H.modelspace->core)
   {
      Orbit& oi = H.modelspace->GetOrbit(i);
      for ( auto& a : VectorUnion(H.modelspace->valence,H.modelspace->qspace) )
      {
      Orbit& oa = H.modelspace->GetOrbit(a);
         Hod.OneBody(a,i) = H.OneBody(a,i);
        Hod.OneBody(i,a) = H.OneBody(a,i)*herm_phase;
      }
   }
//
   for ( auto& i : H.modelspace->valence)
   {
      Orbit& oi = H.modelspace->GetOrbit(i);
      for ( auto& a : H.modelspace->qspace )
      {
      Orbit& oa = H.modelspace->GetOrbit(a);

         Hod.OneBody(a,i) = H.OneBody(a,i);
         Hod.OneBody(i,a) = H.OneBody(a,i)*herm_phase;
      }
   }

//Hod.EraseTwoBody();

   // Two body piece only stored half channel, no need to change
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);

      arma::mat& H2 =  iter.second;

      arma::mat& Hod2 = Hod.TwoBody.GetMatrix(ch_bra,ch_ket);


// PPHH
//
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
	{
Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));
}
}


// PPvc
//
      for ( auto& iket : tbc_ket.GetKetIndex_vc() ) // cc means core-core ('holes' refer to the reference state)
      {
        Ket & dket = tbc_ket.GetKet(iket);
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() )  )
//         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv() ) )
	{
     Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));

}
}


// PPvv
//
      for ( auto& iket : tbc_ket.GetKetIndex_vv() ) // cc means core-core ('holes' refer to the reference state)
      {
        Ket & dket = tbc_ket.GetKet(iket);

         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv() ) )
	{
     Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));
}
}

}


   return Hod;
}




double  Generator::GetVSEOM_Overlap(Operator& H )
{
    // [S+S', S-S']=[S',S] - [S,S']=2[S',S]
  double  ovlp=0;
  double  ovlp1=0;
  double  ovlp2=0;

//  double c2= 0.95542499;
//  double c1= 0.29523395;
// double c1= 1;
//  double c2= 0;
 //
//emax 2
  
 double c1= 0.95161831;
 double c2= 0.30728258;

// double c2= -0.95161831;
// double c1= 0.30728258;
 
// emax 3
//double c1 = 0.9675254 ;
//double c2= 0.25277383 ;
//emax 4
//double c2=  0.24297425 ;
//double c1 = 0.97003274;



    ovlp+=H.ZeroBody;
    

   for ( auto& i : H.modelspace->valence)
   {
      Orbit& oi = H.modelspace->GetOrbit(i);
        if(oi.tz2 != 1)continue;
       // if(std::abs(H.OneBody(i,i))<1e-6)continue;
      //  std::cout << oi.j2 << " " << H.OneBody(i,i) << std::endl;
        if(oi.j2==3)ovlp1+=2*H.OneBody(i,i)*c1*c1;
        if(oi.j2==1)ovlp1+=2*H.OneBody(i,i)*c2*c2;

    //    std::cout << std::setprecision (15) << H.OneBody(i,i) << std::endl;
   }


    //    std::cout << "here one body" << ovlp1 << std::endl;

   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);
      arma::mat& H2 =  iter.second;
      if(tbc_bra.J != 0)continue ;
      if(tbc_bra.Tz != 1)continue ;
      if(tbc_bra.parity != 0)continue ;

      for ( auto& iket : tbc_ket.GetKetIndex_vv() ) // cc means core-core ('holes' refer to the reference state)
      {

        Ket & dket = tbc_bra.GetKet(iket);
         for ( auto& ibra :  tbc_bra.GetKetIndex_vv() )
	{
        Ket & dbra = tbc_bra.GetKet(ibra);
//        std::cout << dbra.p << dbra.q << dket.p<<dket.q <<std::endl;
        if(ibra == iket){
            if(dbra.p == 3 and dbra.q ==3)ovlp2+=c1*c1*H2(ibra,iket);
            if(dbra.p == 5 and dbra.q ==5)ovlp2+=c2*c2*H2(ibra,iket);
  //          std::cout << dbra.p << dbra.q<<dket.p<<dket.q<<" " <<std::setprecision (15)<< H2(ibra,iket) << " " << ovlp2 <<std::endl;
        }
        if(ibra != iket){
            // do not suppose the hermitian of hamiltinian
            if(dbra.p == 3 and dbra.q ==3)ovlp2+=c1*c2*H2(ibra,iket);
            if(dbra.p == 5 and dbra.q ==5)ovlp2+=c1*c2*H2(ibra,iket);
//            std::cout << dbra.p << dbra.q<<dket.p<<dket.q<<" " << std::setprecision (15) << H2(ibra,iket) << " " << ovlp2 <<std::endl;
  //          std::cout << H2(ibra,iket) << "  " << H2(iket,ibra) <<" "<<ovlp2 << " " <<c1*c2*H2(ibra,iket) << " " << std::endl ;
        }
}
}}

//std::cout << " overlap, 0b, 1b,2b : "<<ovlp <<" " << ovlp1 << " " <<ovlp2<< std::endl;
return(ovlp+ovlp1+ovlp2);

}








void  Generator::PrintEOMVS_ladder(Operator& H, int info)
{

    double a[6];
    int idx=0;


   // Two body piece only stored half channel, no need to change
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);

      arma::mat& H2 =  iter.second;

      if(tbc_bra.J !=0 )continue;
      if(tbc_bra.parity !=0 )continue;
      if(tbc_bra.Tz!=1 )continue;



// PPvv
//
      for ( auto& iket : VectorUnion(tbc_ket.GetKetIndex_qq(),tbc_ket.GetKetIndex_qv(),tbc_ket.GetKetIndex_vv())   ) // cc means core-core ('holes' refer to the 2233reference state)
      {
        Ket & dket = tbc_ket.GetKet(iket);
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(),tbc_bra.GetKetIndex_qv(), tbc_bra.GetKetIndex_vv() ) )
	{
        Ket & dbra = tbc_bra.GetKet(ibra);

        if(abs(H2(ibra,iket)) >0.0000001){
            if(info==1) std::cout<< dbra.p<<dbra.q<<dket.p<<dket.q<<" " << H2(ibra,iket) <<std::endl;
        }
  // if(dket.p==5 and dket.q==5){
  //   Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));}
}
}


      for ( auto& iket : tbc_ket.GetKetIndex_vv() ) // cc means core-core ('holes' refer to the reference state)
      {
        Ket & dket = tbc_ket.GetKet(iket);
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv() ) )
	{
        Ket & dbra = tbc_bra.GetKet(ibra);

        if(abs(H2(ibra,iket)) >0.0000001){
            if(info==1) std::cout<< dbra.p<<dbra.q<<dket.p<<dket.q<<std::endl;
        a[idx]=H2(ibra,iket);
        idx+=1;}
  // if(dket.p==5 and dket.q==5){
  //   Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));}
}
}

}


 for (int i = 0; i < 6; ++i) {
        std::cout << a[i] << " "; // Print element and a space
    }
    std::cout << std::endl;

}



Operator  Generator::GetVSEOM_ladder_part(Operator& H , int herm, int part)
{

  // Generate a (anti-)hermit operator, regardless the hermitian of operator of H
  // int herm, 0 for hermit, and 1 for antihermit
    Operator Hod = 0.0* H;
    int herm_phase=1;

  if(herm == 0){
   Hod.SetHermitian();
     herm_phase =1;
}
  if(herm == 1){
   Hod.SetAntiHermitian();
     herm_phase = -1 ;
}
   

   // Two body piece only stored half channel, no need to change
   for ( auto& iter : H.TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = H.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = H.modelspace->GetTwoBodyChannel(ch_ket);

      arma::mat& H2 =  iter.second;

      arma::mat& Hod2 = Hod.TwoBody.GetMatrix(ch_bra,ch_ket);


      if(part == 0){
// PPHH
//
      for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
	{
//Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));
}
}}


    if(part==1){
// PPvq
//
      for ( auto& iket : tbc_ket.GetKetIndex_vc() ) // cc means core-core ('holes' refer to the reference state)
      {
        Ket & dket = tbc_ket.GetKet(iket);
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() )  )
	{
     Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));
}
}
}


if(part==2){
// PPvv
//
      for ( auto& iket : tbc_ket.GetKetIndex_vv() ) // cc means core-core ('holes' refer to the reference state)
      {
        Ket & dket = tbc_ket.GetKet(iket);

         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_qv() ) )
	{
     Hod.TwoBody.AddToTBME(ch_bra,ch_ket, ibra, iket, H2(ibra,iket));
}
}
}

}


   return Hod;
}

