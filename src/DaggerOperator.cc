
#include "DaggerOperator.hh"
#include "IMSRGProfiler.hh"



//////////////////// DESTRUCTOR //////////////////////////////////////////
DaggerOperator::~DaggerOperator()
{
  profiler.counter["N_Operators"] --;
}

/////////////////// CONSTRUCTORS /////////////////////////////////////////
DaggerOperator::DaggerOperator() : Operator() 
{}

DaggerOperator::DaggerOperator(ModelSpace& ms) : Operator( ms) ///< Construct a 2-body scalar operator
{}

DaggerOperator::DaggerOperator(ModelSpace& ms, int Jrank, int Trank, int Parity, int part_rank) : Operator(ms,  Jrank, Trank, Parity, part_rank)
{}

DaggerOperator::DaggerOperator( const DaggerOperator& rhs) : Operator( rhs ) ///< Copy constructor
{}

DaggerOperator::DaggerOperator( DaggerOperator&& rhs) : Operator( rhs )
{}


DaggerOperator::DaggerOperator(ModelSpace& ms, index_t Q) : Operator( ms )
{ 
  SetQSpaceOrbit(Q);
  OneBody(Q,Q)= 1.0;
}




DaggerOperator& DaggerOperator::operator=( const DaggerOperator& rhs) = default;






//Operator Operator::Commutator( Operator& opright)
/// Returns \f$ Z = [X,Y] \f$
/// @relates DaggerOperator
//NOT SURE IF I NEED TO PUT THIS HERE OR MAKE THE FRIEND METHOD IN OPERATOR.CC
//RECOGNIZE A DAGGER OPERATOR.
DaggerOperator Commutator( const Operator& X, const DaggerOperator& Y)
{
  int jrank = max(X.rank_J,Y.rank_J);
  int trank = max(X.rank_T,Y.rank_T);
  int parity = (X.parity+Y.parity)%2;
  int particlerank = max(X.particle_rank,Y.particle_rank);
  DaggerOperator Z(*(Y.modelspace),jrank,trank,parity,particlerank);
  Z.SetToCommutator(X,Y);
  return Z;
}


void DaggerOperator::SetToCommutator( const Operator& X, const DaggerOperator& Y)
{
//   profiler.counter["N_Commutators"] += 1;
   double t_start = omp_get_wtime();
   DaggerOperator& Z = *this;
   modelspace->PreCalculateSixJ();
   int xrank = X.rank_J + X.rank_T + X.parity;
//   int yrank = Y.rank_J + Y.rank_T + Y.parity;
   if (xrank==0)
   {
         Z.CommutatorScalarDagger(X,Y); // [S,S]
   }
   else
   {
      std::cout << "!!!!!!!!!!!!TROUBLE!!!!!!!!!!!!!!!!" << std::endl;
      cout << " Tensor-Dagger commutator not yet implemented." << endl;
   }
   profiler.timer["Commutator"] += omp_get_wtime() - t_start;
}




/// Commutator where \f$ X \f$ is a scalar and \f$Y\f$ is a dagger.
/// Should be called through Commutator()
void DaggerOperator::CommutatorScalarDagger( const Operator& X, const DaggerOperator& Y) 
{
   profiler.counter["N_DaggerCommutators"] += 1;
   double t_css = omp_get_wtime();
   DaggerOperator& Z = *this;
   Z = Y;
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseTwoBody();

   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();

//   double t_start = omp_get_wtime();
   Z.comm211sd( X, Y ) ; 
   Z.comm231sd( X, Y ) ;
//   comm431sd( X, Y ) ;
   Z.comm413_233sd( X, Y ) ; 
//   comm433sd_pphh( X, Y ) ; 
   Z.comm433sd_ph( X, Y ) ; 
   Z.comm433_pp_hh_431sd( X, Y ) ; 


//   // NEED TO FIGURE OUT IF SYMMETRIZE IS A SANE THING TO DO WITH DAGGER OPERATORS
//   if ( Z.IsHermitian() )
//      Z.Symmetrize();
//   else if (Z.IsAntiHermitian() )
//      Z.AntiSymmetrize();

   profiler.timer["CommutatorScalarScalar"] += omp_get_wtime() - t_css;

}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////// BEGIN IMPLEMENTATION OF COMMUTATOR EXPRESSIONS //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


//*****************************************************************************************
// [X^(2), Y^(1)]^(1)
//
//      |i                 i|   
//      |                   |   
//     (X)         ---      |  (Y)
//       a\                 |  /a
//         (Y)             (X)/
//           
//
// Sum_ia X_ia Y_a^{lambda}
// Adapted from
//    void Operator::comm111ss( const Operator & X, const Operator& Y) 
// This could be modified since we only need one column from Y to go to one column in Z.
void DaggerOperator::comm211sd( const Operator& X, const DaggerOperator& Y )
{
   DaggerOperator& Z = *this;
   Z.OneBody += X.OneBody*Y.OneBody;
}



//*****************************************************************************************
// [X^(2),Y^(3)]^(1)
//         |i             |i
//         |              |
//  (X)_ b |         a _(Y)
//    \_\  |   ---    /_/
//    a  (Y)         (X) b
//
// Sum_iab Sum_J X_ab Y_bia  hat(J)/hat(lambda)
//
//
// Adapted from
//    void Operator::comm121ss( const Operator& X, const Operator& Y) 
void DaggerOperator::comm231sd( const Operator& X, const DaggerOperator& Y) 
{
   DaggerOperator& Z = *this;
   index_t norbits = modelspace->GetNumberOrbits();
   index_t j = Z.GetQSpaceOrbit();
   Orbit& oj = modelspace->GetOrbit(j);
//   for (index_t i=0;i<norbits;++i)
//   #pragma omp parallel for 
   for (auto i : OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
   {
//      Orbit &oi = modelspace->GetOrbit(i);
      // j is Q-space orbit
//      index_t jmin = Z.IsNonHermitian() ? 0 : i;
//      for (auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
//      {
//          if (j<jmin) continue; // only calculate upper triangle
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             for (index_t b=0; b<norbits; ++b)
             {
                Orbit &ob = modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ);
                if (abs(nanb)<OCC_CUT) continue;
                if (Y.particle_rank>1)
                {
                  Z.OneBody(i,j) += (ob.j2+1) * nanb *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                  Z.OneBody(i,j) -= (oa.j2+1) * nanb *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,j) ;
                }
                // comm211 part
//                if (X.particle_rank>1)
//                {
//                  Z.OneBody(i,j) -= (ob.j2+1) * nanb * Y.OneBody(a,b) * X.TwoBody.GetTBMEmonopole(b,i,a,j) ;
//                  Z.OneBody(i,j) += (oa.j2+1) * nanb * Y.OneBody(b,a) * X.TwoBody.GetTBMEmonopole(a,i,b,j) ;
//                }
             }
          }
//      }
   }
}



//*****************************************************************************************
// [X^(4),Y^(3)]^(1)
//
//      i |                  i|  
//       (Y) ___              |   _(X) 
//          |\  \             |  / / |
//          | \  \      __    | | /  /                                                                                 
//          \  \ |            | |/__/
//           \_(X)            (Y)      
//                                       
//  Zij = sum 1/4 1/(2j_i +1)  sum_abcJ (2J+1)  (na nb n_c - n_a n_b nc) AJ_ciab BJ_abc
//
//  Adapted from
//    void Operator::comm221ss( const Operator& X, const Operator& Y) 
void DaggerOperator::comm431sd( const Operator& X, const DaggerOperator& Y)
{

   DaggerOperator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();

   static TwoBodyME Mpp = Z.TwoBody;
   static TwoBodyME Mhh = Z.TwoBody;

   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo.
   int nch = modelspace->SortedTwoBodyChannels.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = X.TwoBody.GetMatrix(ch,ch);
      auto& RHS = Y.TwoBody.GetMatrix(ch,ch);

      auto& Matrixpp = Mpp.GetMatrix(ch,ch);
      auto& Matrixhh = Mhh.GetMatrix(ch,ch);

      auto& kets_pp = tbc.GetKetIndex_pp();
      auto& kets_hh = tbc.GetKetIndex_hh();
      auto& kets_ph = tbc.GetKetIndex_ph();
      auto& nanb = tbc.Ket_occ_hh;
      auto& nbarnbar_hh = tbc.Ket_unocc_hh;
      auto& nbarnbar_ph = tbc.Ket_unocc_ph;
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;

// Commented this out because it accounts for the AB + BA structure, which shouldn't be there
// for the dagger operator. -SRS 14/04/2018
//      if (Z.IsHermitian())
//      {
//         Matrixpp +=  Matrixpp.t();
//         Matrixhh +=  Matrixhh.t();
//      }
//      else if (Z.IsAntiHermitian()) // i.e. LHS and RHS are both hermitian or ant-hermitian
//      {
//         Matrixpp -=  Matrixpp.t();
//         Matrixhh -=  Matrixhh.t();
//      }
//      else
//      {
//        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
//        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
//        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
//        if (kets_ph.size()>0)
//          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
//      }


   } //for ch

   index_t j = Z.GetQSpaceOrbit();
   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
//      int jmin = Z.IsNonHermitian() ? 0 : i;
//      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
//      {
//         if (j<jmin) continue;
         double cijJ = 0;
         for (int ch=0;ch<nChannels;++ch)
         {
            TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
            double Jfactor = (2*tbc.J+1.0);
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : modelspace->holes)
            {
               Orbit& oc = modelspace->GetOrbit(c);
               cijJ += Jfactor * oc.occ     * Mpp.GetTBME(ch,c,i,c,j);
               cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,j);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,j);
            }
         }
         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
//      } // for j
   }


}





//*****************************************************************************************
// [X^(2),Y^(3)]^(3)  and [X^(4),Y^(1)]^(3), combined
//
//     |   |  /       |  | /             \ /         \ |  (Y)
//    (X)  | /   __   | (Y)              (X)    __    \| / 
//      \  |/         | /        and     /  \         (X)
//       (Y)         (X)                /   (Y)        | 
//
//  ZJ_ijkQ = sum_a  X_ia YJ_ajkQ  + X_ja YJ_iakQ + X_ak YJ_ijaQ + XJ_ijka Y_aQ
//
//  I adapted this, but then gave up and brute forced it because the intermediate matrix
//  stuff was too obfuscated, and I'm the one who wrote it... -SRS 14/04/2018
// Adapted from
//  void Operator::comm122ss( const Operator& X, const Operator& Y ) 
void DaggerOperator::comm413_233sd( const Operator& X, const DaggerOperator& Y) 
{

   DaggerOperator& Z = *this;
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;
//   int hZ = Z.IsHermitian() ? 1 : -1;

   index_t Qorbit = Z.GetQSpaceOrbit();
   Orbit& oQ = modelspace->GetOrbit(Qorbit);

   int n_nonzero = modelspace->SortedTwoBodyChannels.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
//      auto& X2 = X.TwoBody.GetMatrix(ch,ch);
//      auto& Y2 = Y.TwoBody.GetMatrix(ch,ch);
//      auto& Z2 = Z.TwoBody.GetMatrix(ch,ch);
//      arma::mat W2(size(Z2),arma::fill::zeros); // temporary intermediate matrix

      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = modelspace->GetOrbit(i);
         Orbit& oj = modelspace->GetOrbit(j);
         for ( int k : OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
         {
           Orbit& ok = modelspace->GetOrbit(k);
           double cijk = 0;
           for ( int a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
           {
             cijk += X1(i,a) * Y.TwoBody.GetTBME_norm(ch,ch,a,j,k,Qorbit);
           }
           for ( int a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
           {
             cijk += X1(j,a) * Y.TwoBody.GetTBME_norm(ch,ch,i,a,k,Qorbit);
           }
           for ( int a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
           {
             cijk += X1(a,k) * Y.TwoBody.GetTBME_norm(ch,ch,i,j,a,Qorbit);
           }
           for ( int a : Y.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
           {
             cijk += Y1(a,Qorbit) * X.TwoBody.GetTBME_norm(ch,ch,i,j,k,a);
           }
           Z.TwoBody.SetTBME(ch,ch,i,j,k,Qorbit, cijk);
         }
      }
   }

}





//*****************************************************************************************
// [X^(4),Y^(3)]^(3)  pp - hh contribution
//  
//  ZJ_ijkQ = 1/2 sum_ab XJ_ijab YJ_abkQ  (_na _nb - na nb )
//
//    | |      
//    (X)           \ \   /(Y)
//    | |    __      \ \ / / |
//    (Y)             \(X)/  |
//    |                      |
//
//
// Adapted from
//    void Operator::comm222_pp_hhss( const Operator& X, const Operator& Y ) 
void DaggerOperator::comm433sd_pphh( const Operator& X, const DaggerOperator& Y)  
{
   DaggerOperator& Z = *this;

   static TwoBodyME Mpp = Z.TwoBody;
   static TwoBodyME Mhh = Z.TwoBody;

   double t = omp_get_wtime();
   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo.
   int nch = modelspace->SortedTwoBodyChannels.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = X.TwoBody.GetMatrix(ch,ch);
      auto& RHS = Y.TwoBody.GetMatrix(ch,ch);
      auto& OUT = Z.TwoBody.GetMatrix(ch,ch);

      auto& Matrixpp = Mpp.GetMatrix(ch,ch);
      auto& Matrixhh = Mhh.GetMatrix(ch,ch);

      auto& kets_pp = tbc.GetKetIndex_pp();
      auto& kets_hh = tbc.GetKetIndex_hh();
      auto& kets_ph = tbc.GetKetIndex_ph();
      auto& nanb = tbc.Ket_occ_hh;
//      auto& nabar_nbbar = tbc.Ket_unocc_hh;
      auto& nbarnbar_hh = tbc.Ket_unocc_hh;
      auto& nbarnbar_ph = tbc.Ket_unocc_ph;
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;


//      if (Z.IsHermitian())
//      {
//         Matrixpp +=  Matrixpp.t();
//         Matrixhh +=  Matrixhh.t();
//      }
//      else if (Z.IsAntiHermitian()) // i.e. LHS and RHS are both hermitian or ant-hermitian
//      {
//         Matrixpp -=  Matrixpp.t();
//         Matrixhh -=  Matrixhh.t();
//      }
//      else
//      {
//        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
//        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
//        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
//        if (kets_ph.size()>0)
//          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
//      }


      // The two body part
      OUT += Matrixpp - Matrixhh;
   } //for ch
   profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;
}








void DaggerOperator::ConstructDaggerMpp_Mhh(const Operator& X, const Operator& Y, TwoBodyME& Mpp, TwoBodyME& Mhh) const
{
   int nch = modelspace->SortedTwoBodyChannels.size();
//   Operator& Z = *this;
//   bool z_is_hermitian = IsHermitian();
//   bool z_is_antihermitian = IsAntiHermitian();
//   bool z_is_hermitian = X.IsHermitian() xor Y.IsHermitian();
//   bool z_is_antihermitian = (X.IsHermitian() == Y.IsHermitian()) and (X.IsAntiHermitian() == Y.IsAntiHermitian());
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = X.TwoBody.GetMatrix(ch,ch);
      auto& RHS = Y.TwoBody.GetMatrix(ch,ch);
//      auto& OUT = Z.TwoBody.GetMatrix(ch,ch);

      auto& Matrixpp = Mpp.GetMatrix(ch,ch);
      auto& Matrixhh = Mhh.GetMatrix(ch,ch);

      auto& kets_pp = tbc.GetKetIndex_pp();
      auto& kets_hh = tbc.GetKetIndex_hh();
      auto& kets_ph = tbc.GetKetIndex_ph();
      auto& nanb = tbc.Ket_occ_hh;
      auto& nbarnbar_hh = tbc.Ket_unocc_hh;
      auto& nbarnbar_ph = tbc.Ket_unocc_ph;
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;


//      if (z_is_hermitian)
//      {
//         Matrixpp +=  Matrixpp.t();
//         Matrixhh +=  Matrixhh.t();
//      }
//      else if (z_is_antihermitian) // i.e. LHS and RHS are both hermitian or ant-hermitian
//      {
//         Matrixpp -=  Matrixpp.t();
//         Matrixhh -=  Matrixhh.t();
//      }
//      else
//      {
//        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
//        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
//        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
//        if (kets_ph.size()>0)
//          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
//      }


      // The two body part
//      OUT += Matrixpp - Matrixhh;
   } //for ch

}


/// Since comm443_pp_hhsd() and comm431sd() both require the construction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
/// We do the same thing in the standard scalar-scalar commutators.
void DaggerOperator::comm433_pp_hh_431sd( const Operator& X, const Operator& Y )  
{

//   int herm = Z.IsHermitian() ? 1 : -1;
   DaggerOperator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();

   static TwoBodyME Mpp = Z.TwoBody;
   static TwoBodyME Mhh = Z.TwoBody;

   double t = omp_get_wtime();
   ConstructDaggerMpp_Mhh( X, Y, Mpp, Mhh);

//   Z.TwoBody += (Mpp - Mhh);
   Z.TwoBody += Mpp;
   Z.TwoBody -= Mhh;
//   OUT += Matrixpp - Matrixhh;
   profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;

   t = omp_get_wtime();
   // The one body part
   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      int jmin = Z.IsNonHermitian() ? 0 : i;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         if (j<jmin) continue;
         double cijJ = 0;
         for (int ch=0;ch<nChannels;++ch)
         {
            TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
            double Jfactor = (2*tbc.J+1.0);
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : modelspace->holes)
            {
               Orbit& oc = modelspace->GetOrbit(c);
               cijJ += Jfactor * oc.occ * Mpp.GetTBME(ch,c,i,c,j); 
               cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,j); 
            }
            // Sum c over particles and include the n_a * n_b terms
            for (auto& c : modelspace->particles)
            {
               cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,j);
            }
         }
         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   } // for i
   profiler.timer["pphh One Body bit"] += omp_get_wtime() - t;
}



//*****************************************************************************************
// [X^(4),Y^(3)]^(3)]  ph piece
//                                             
//   |           |       |           |           
//   |      _(Y)_|       |_(X)_      |            
//   |     /\            |    /\     |
//   |    (  )      __   |   (  )    |            
//   |_(X)_\/            |    \/_(Y)_|            
//   |                   |           |
//
// For formula and details of implementation, see Operator::comm222_phss. The only change
// made here to accommodate a dagger operator is to replace XY - YX  with -YX. This
// ensures that we don't include contributions from the Qspace orbit to X.
// Adapted from
//    void Operator::comm222_phss( const Operator& X, const Operator& Y ) 
void DaggerOperator::comm433sd_ph( const Operator& X, const DaggerOperator& Y)  
{

   DaggerOperator& Z = *this;
   int hy = Y.IsHermitian() ? 1 : -1;
   double t_start = omp_get_wtime();


   // Construct the intermediate matrix Z_bar
   const auto& pandya_lookup = modelspace->GetPandyaLookup(rank_J, rank_T, parity);
   int nch = modelspace->SortedTwoBodyChannels_CC.size();
   t_start = omp_get_wtime();
   deque<arma::mat> Z_bar (nChannels );
   vector<bool> lookup_empty(nChannels,true);
   for (int ich=0;ich<nch;++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels_CC[ich];
      index_t nKets_cc = modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros( nKets_cc, 2*nKets_cc );
      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
   }


   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->scalar_transform_first_pass)
   #endif
   for (int ich=0; ich<nch; ++ich )
   {
      if (lookup_empty.at(ich)) continue;
      int ch = modelspace->SortedTwoBodyChannels_CC.at(ich);
      const TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

      arma::mat Y_bar_ph;
      arma::mat Xt_bar_ph;

      Y.DoPandyaTransformation_SingleChannel(Y_bar_ph,ch,"normal");
      X.DoPandyaTransformation_SingleChannel(Xt_bar_ph,ch,"transpose");
      auto& Zbar_ch = Z_bar.at(ch);


      if (Y_bar_ph.size()<1 or Xt_bar_ph.size()<1)
      {
        Zbar_ch = arma::zeros( Xt_bar_ph.n_rows, Y_bar_ph.n_cols*2);
        continue;
      }

      // get the phases for taking the transpose
      arma::mat PhaseMat(nKets_cc, nKets_cc, arma::fill::ones );
      for (index_t iket=0;iket<nKets_cc;iket++)
      {
         const Ket& ket = tbc_cc.GetKet(iket);
         if ( modelspace->phase( (ket.op->j2 + ket.oq->j2)/2 ) > 0) continue;
         PhaseMat.col( iket ) *= -1;
         PhaseMat.row( iket ) *= -1;
      }
      arma::uvec phkets = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
      auto PhaseMatY = PhaseMat.rows(phkets) * hy;

//                                           [      |     ]
//     create full Y matrix from the half:   [  Yhp | Y'ph]   where the prime indicates multiplication by (-1)^(i+j+k+l) h_y
//                                           [      |     ]   Flipping hp <-> ph and multiplying by the phase is equivalent to
//                                           [  Yph | Y'hp]   having kets |kj> with k>j.
      Zbar_ch =  Xt_bar_ph * join_horiz(Y_bar_ph, join_vert(   Y_bar_ph.tail_rows(nph_kets)%PhaseMatY,
                                                               Y_bar_ph.head_rows(nph_kets)%PhaseMatY) );



      // If Z is hermitian, then XY is anti-hermitian, and so XY - YX = XY + (XY)^T
      // For the dagger operator, we're actually only interested in the -YX piece
      // so the only change from the scalar-scalar case is += becomes = and -= becomes = -.
      // I hope that works. -SRS 14/04/2018
      if ( Z.IsHermitian() )
      {
         Zbar_ch.head_cols(nKets_cc) = Zbar_ch.head_cols(nKets_cc).t();
      }
      else
      {
         Zbar_ch.head_cols(nKets_cc) - -Zbar_ch.head_cols(nKets_cc).t();
      }
      Zbar_ch.tail_cols(nKets_cc) += Zbar_ch.tail_cols(nKets_cc).t()%PhaseMat;

   }

   profiler.timer["Build Z_bar"] += omp_get_wtime() - t_start;

   // Perform inverse Pandya transform on Z_bar to get Z
   t_start = omp_get_wtime();
   Z.AddInversePandyaTransformation(Z_bar);

   modelspace->scalar_transform_first_pass = false;
   profiler.timer["InversePandyaTransformation"] += omp_get_wtime() - t_start;

}



