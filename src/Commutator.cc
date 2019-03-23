
#include "Commutator.hh"
#include "Operator.hh"
#include "TwoBodyME.hh"
#include "ThreeBodyME.hh"
#include "armadillo"
#include <map>
#include <deque>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>



namespace Commutator {

bool use_goose_tank_correction = false;
bool use_brueckner_bch = false;
bool use_imsrg3 = false;
double bch_transform_threshold = 1e-9;
double bch_product_threshold = 1e-4;


void Set_BCH_Transform_Threshold(double x)
{bch_transform_threshold=x;}

void Set_BCH_Product_Threshold(double x)
{bch_product_threshold=x;}

void SetUseBruecknerBCH(bool tf)
{use_brueckner_bch = tf;}

void SetUseGooseTank(bool tf)
{use_goose_tank_correction = tf;}

void SetUseIMSRG3(bool tf)
{use_imsrg3 = tf;}


//Operator Operator::Commutator( Operator& opright)
/// Returns \f$ Z = [X,Y] \f$
Operator Commutator( const Operator& X, const Operator& Y)
{
//  int jrank = std::max(X.rank_J,Y.rank_J);
//  int trank = std::max(X.rank_T,Y.rank_T);
//  int parity = (X.parity+Y.parity)%2;
//  int particlerank = std::max(X.particle_rank,Y.particle_rank);
  int xrank = X.rank_J + X.rank_T + X.parity;
  int yrank = Y.rank_J + Y.rank_T + Y.parity;
  int xlegs = X.GetNumberLegs();
  int ylegs = Y.GetNumberLegs();

  X.modelspace->PreCalculateSixJ();

  if (xrank==0)
   {
      if ( (xlegs%2==0) and (ylegs%2==1) )
      {
        return CommutatorScalarDagger(X,Y);
      }
      if ( (xlegs%2==1) and (ylegs%2==0) )
      {
        return -CommutatorScalarDagger(Y,X);
      }
      if (yrank==0)
      {
         return CommutatorScalarScalar(X,Y); // [S,S]
      }
      else
      {
         return CommutatorScalarTensor(X,Y); // [S,T]
      }
   }
   else if(yrank==0)
   {
      return -CommutatorScalarTensor(Y,X); // [T,S]
   }
   else
   {
      std::cout << "In Tensor-Tensor because X.rank_J = " << X.rank_J << "  X.rank_T = " << X.rank_T << "  X.parity = " << X.parity << "   ";
      std::cout <<                        "  Y.rank_J = " << Y.rank_J << "  Y.rank_T = " << Y.rank_T << "  Y.parity = " << Y.parity << std::endl;
      std::cout << " Tensor-Tensor commutator not yet implemented." << std::endl;
   }
   return  0*Y;
}


/// Commutator where \f$ X \f$ and \f$Y\f$ are scalar operators.
/// Should be called through Commutator()
Operator CommutatorScalarScalar( const Operator& X, const Operator& Y) 
{
   std::cout << "Calling CommutatorScalarScalar..." << std::endl;
   X.profiler.counter["N_ScalarCommutators"] += 1;
   double t_css = omp_get_wtime();
   Operator Z( *(Y.GetModelSpace()), std::max(X.GetJRank(),Y.GetJRank()), std::max(X.GetTRank(),Y.GetTRank()), (X.GetParity()+Y.GetParity())%2, std::max(X.GetParticleRank(),Y.GetParticleRank()) );

   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();

   if ( not Z.IsAntiHermitian() )
   {
      comm110ss(X, Y, Z);
      if (X.particle_rank>1 and Y.particle_rank>1)
        comm220ss(X, Y, Z) ;
   }

   double t_start = omp_get_wtime();
//   Z.comm111ss(X, Y);
   comm111ss(X, Y, Z);
   X.profiler.timer["comm111ss"] += omp_get_wtime() - t_start;

    t_start = omp_get_wtime();
//   Z.comm121ss(X,Y);
   comm121ss(X, Y, Z);
   X.profiler.timer["comm121ss"] += omp_get_wtime() - t_start;

    t_start = omp_get_wtime();
//   Z.comm122ss(X,Y); 
   comm122ss(X, Y, Z); 
   X.profiler.timer["comm122ss"] += omp_get_wtime() - t_start;

   if (X.particle_rank>1 and Y.particle_rank>1)
   {
     t_start = omp_get_wtime();
//     Z.comm222_pp_hh_221ss(X, Y);
     comm222_pp_hh_221ss(X, Y, Z);
     X.profiler.timer["comm222_pp_hh_221ss"] += omp_get_wtime() - t_start;
      
     t_start = omp_get_wtime();
//     Z.comm222_phss(X, Y);
     comm222_phss(X, Y, Z);
     X.profiler.timer["comm222_phss"] += omp_get_wtime() - t_start;
   }

   std::cout << "Made it here. About to check if use_imsrg3 is true" << std::endl;

   if (use_imsrg3)
   {
     if (X.GetParticleRank()>2 and Y.GetParticleRank()>2)
     {
       t_start = omp_get_wtime();
       comm330ss(X, Y, Z);
       X.profiler.timer["comm330ss"] += omp_get_wtime() - t_start;
       t_start = omp_get_wtime();
       comm331ss(X, Y, Z);
       X.profiler.timer["comm331ss"] += omp_get_wtime() - t_start;
     }

     t_start = omp_get_wtime();
     comm231ss(X, Y, Z);
     X.profiler.timer["comm231ss"] += omp_get_wtime() - t_start;

     t_start = omp_get_wtime();
     comm132ss(X, Y, Z);
     X.profiler.timer["comm132ss"] += omp_get_wtime() - t_start;

     t_start = omp_get_wtime();
     comm232ss(X, Y, Z);
     X.profiler.timer["comm232ss"] += omp_get_wtime() - t_start;

     t_start = omp_get_wtime();
     comm133ss(X, Y, Z);
     X.profiler.timer["comm133ss"] += omp_get_wtime() - t_start;

//     std::cout << "Begin comm223ss. particle rank of Z = " << Z.GetParticleRank() << std::endl;
     t_start = omp_get_wtime();
     comm223ss(X, Y, Z);
     X.profiler.timer["comm223ss"] += omp_get_wtime() - t_start;
   }


   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();

   X.profiler.timer["CommutatorScalarScalar"] += omp_get_wtime() - t_css;
   return Z;

}


/// Commutator \f$[X,Y]\f$ where \f$ X \f$ is a scalar operator and \f$Y\f$ is a tensor operator.
/// Should be called through Commutator()
Operator CommutatorScalarTensor( const Operator& X, const Operator& Y) 
{
   X.profiler.counter["N_TensorCommutators"] += 1;
   double t_cst = omp_get_wtime();
   Operator Z = Y; // This ensures the commutator has the same tensor rank as Y
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseTwoBody();

   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();

   double t_start = omp_get_wtime();
   comm111st(X, Y, Z);
   X.profiler.timer["comm111st"] += omp_get_wtime() - t_start;
   t_start = omp_get_wtime();
   comm121st(X, Y, Z);
   X.profiler.timer["comm121st"] += omp_get_wtime() - t_start;

   t_start = omp_get_wtime();
   comm122st(X, Y, Z);
   X.profiler.timer["comm122st"] += omp_get_wtime() - t_start;
   t_start = omp_get_wtime();
   comm222_pp_hh_221st(X, Y, Z);
   X.profiler.timer["comm222_pp_hh_221st"] += omp_get_wtime() - t_start;
   t_start = omp_get_wtime();
   comm222_phst(X, Y, Z);
   X.profiler.timer["comm222_phst"] += omp_get_wtime() - t_start;

   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();
   X.profiler.timer["CommutatorScalarTensor"] += omp_get_wtime() - t_cst;
   return Z;
}




/// Commutator where \f$ X \f$ is a scalar and \f$Y\f$ is a dagger (i.e. it creates an additional particle).
/// Should be called through Commutator()
Operator CommutatorScalarDagger( const Operator& X, const Operator& Y) 
{
   X.profiler.counter["N_DaggerCommutators"] += 1;
   double t_css = omp_get_wtime();
   Operator Z = Y;
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseTwoBody();

//   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
//   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
//   else Z.SetNonHermitian();

   comm211sd( X, Y, Z ) ; 
   comm231sd( X, Y, Z ) ;
   comm413_233sd( X, Y, Z ) ; 
   comm433_pp_hh_431sd( X, Y, Z ) ; 
   comm433sd_ph( X, Y, Z ) ; 
//   comm433sd_ph_dumbway( X, Y, Z ) ; 

   X.profiler.timer["CommutatorScalarDagger"] += omp_get_wtime() - t_css;
   return Z;

}







//*****************************************************************************************
/// BCH_Transform(X,Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
/// with all commutators truncated at the two-body level.
Operator BCH_Transform(  const Operator& OpIn, const Operator& Omega)
{
   return use_brueckner_bch ? Brueckner_BCH_Transform( OpIn, Omega ) :  Standard_BCH_Transform( OpIn, Omega );
}

/// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
/// with all commutators truncated at the two-body level.
Operator Standard_BCH_Transform( const Operator& OpIn, const Operator &Omega)
{
   double t_start = omp_get_wtime();
   int max_iter = 40;
   int warn_iter = 12;
   double nx = OpIn.Norm();
   double ny = Omega.Norm();
   Operator OpOut = OpIn;
   double factorial_denom = 1.0;
   Operator goosetank_chi;  // auxiliary one-body operator used to recover 4th-order quadruples.
   if (use_goose_tank_correction)
   {
     goosetank_chi = OpIn;
     goosetank_chi.SetParticleRank(1);
     goosetank_chi.Erase();
   }
   if (nx>bch_transform_threshold)
   {
     Operator OpNested = OpIn;
     double epsilon = nx * exp(-2*ny) * bch_transform_threshold / (2*ny);
     for (int i=1; i<=max_iter; ++i)
     {

        if (use_goose_tank_correction  )
        {
          auto chi_last = goosetank_chi.OneBody;
          goosetank_chi = GooseTankUpdate( Omega, OpNested);
          OpNested.OneBody += chi_last;  // add the chi from the previous step to OpNested.
        }
        

        OpNested = Commutator(Omega,OpNested); // the ith nested commutator
        factorial_denom /= i;
        OpOut += factorial_denom * OpNested;
  
        if (OpOut.rank_J > 0)
        {
            std::cout << "Tensor BCH, i=" << i << "  Norm = " << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.OneBodyNorm() << " " 
                                                              << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.TwoBodyNorm() << " "
                                                              << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.Norm() << std::endl;
        }
        epsilon *= i+1;
        if (OpNested.Norm() < epsilon)  break;
        if (i == warn_iter)  std::cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << std::endl;
        else if (i == max_iter)   std::cout << "Warning: BCH_Transform didn't coverge after "<< max_iter << " nested commutators" << std::endl;
     }
   }
   std::cout << "Done with BCH_Transform, 3-body norm of OpOut = " << OpOut.ThreeBodyNorm() << std::endl;
   OpIn.profiler.timer["BCH_Transform"] += omp_get_wtime() - t_start;
   return OpOut;
}


//  Update the auxiliary one-body operator chi, using Omega and the ith nested commutator
//  This has not been tested for tensor commutators, but it *should* work.
//  Of course, there's not as clean a motivation in terms of perturbation theory for the tensors...
// 
Operator GooseTankUpdate( const Operator& Omega, const Operator& OpNested)
{
   double t_start = omp_get_wtime();
   Operator goosetank_chi = Operator( *(OpNested.modelspace), OpNested.rank_J, OpNested.rank_T, OpNested.parity, 1) ;
   goosetank_chi.EraseOneBody();
   if (goosetank_chi.rank_J==0 )
   {
     comm221ss( Omega, OpNested, goosetank_chi );  // update chi.
   }
   else
   {
     comm222_pp_hh_221st( Omega, OpNested, goosetank_chi );  // update chi.
   }
   goosetank_chi.Symmetrize(); // the commutator call only does half the matrix, so we symmetrize
   int norbits = OpNested.modelspace->GetNumberOrbits();
   for (int i=0;i<norbits;++i)  // enforce n_in_j + nbar_i nbar_j
   {
     Orbit &oi = OpNested.modelspace->GetOrbit(i);
     for (int j=0;j<norbits;++j)
     {
      Orbit &oj = OpNested.modelspace->GetOrbit(j);
      goosetank_chi.OneBody(i,j) *=  oi.occ*oj.occ + (1.0-oi.occ)*(1.0-oj.occ) ;
      }
   }
   OpNested.profiler.timer["GooseTankUpdate"] += omp_get_wtime() - t_start;
   return goosetank_chi;
}




/// Variation of the BCH transformation procedure
/// requested by a one Dr. T.D. Morris
/// \f[ e^{\Omega_1 + \Omega_2} X e^{-\Omega_1 - \Omega_2}
///    \rightarrow 
///  e^{\Omega_2} e^{\Omega_1}  X e^{-\Omega_1} e^{-\Omega_2} \f]
Operator Brueckner_BCH_Transform( const Operator& OpIn, const Operator& Omega)
{
   Operator Omega1 = Omega;
   Operator Omega2 = Omega;
   Omega1.SetParticleRank(1);
   Omega1.EraseTwoBody();
   Omega2.EraseOneBody();
   Operator OpOut = Standard_BCH_Transform(OpIn, Omega1);
   OpOut = Standard_BCH_Transform(OpOut, Omega2);
   return OpOut;
}


//*****************************************************************************************
// Baker-Campbell-Hausdorff formula
//  returns Z, where
//  exp(Z) = exp(X) * exp(Y).
//  Z = X + Y + 1/2[X, Y]
//     + 1/12 [X,[X,Y]] + 1/12 [Y,[Y,X]]
//     - 1/24 [Y,[X,[X,Y]]]
//     - 1/720 [Y,[Y,[Y,[Y,X]]]] - 1/720 [X,[X,[X,[X,Y]]]]
//     + ...
//*****************************************************************************************
/// X.BCH_Product(Y) returns \f$Z\f$ such that \f$ e^{Z} = e^{X}e^{Y}\f$
/// by employing the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + Y + \frac{1}{2}[X,Y] + \frac{1}{12}([X,[X,Y]]+[Y,[Y,X]]) + \ldots \f]
//*****************************************************************************************
Operator BCH_Product(  Operator& X, Operator& Y)
{
   double tstart = omp_get_wtime();
   double nx = X.Norm();
   std::vector<double> bernoulli = {1.0, -0.5, 1./6, 0.0, -1./30,  0.0 ,  1./42,     0,  -1./30};
   std::vector<double> factorial = {1.0,  1.0,  2.0, 6.0,    24.,  120.,   720., 5040.,  40320.};


   Operator Z = X + Y;
//   if (use_goose_tank_correction) return Z; // Not sure why this is here
//   Operator Nested = Y;
//   Nested.SetToCommutator(Y,X)
   Operator Nested = Commutator(Y,X);  // [Y,X]


   double nxy = Nested.Norm();
   // We assume X is small, but just in case, we check if we should include the [X,[X,Y]] term.
   if ( nxy*nx > bch_product_threshold)
   {
     Z += (1./12) * Commutator(Nested,X);
   }
   
   int k = 1;
   // k=1 adds 1/2[X,Y],  k=2 adds 1/12 [Y,[Y,X]], k=4 adds -1/720 [Y,[Y,[Y,[Y,X]]]], and so on.
   while( Nested.Norm() > bch_product_threshold and k<9)
   {
     if ((k<2) or (k%2==0))
        Z += (bernoulli[k]/factorial[k]) * Nested;

     Nested = Commutator(Y,Nested);
     k++;
   }

   X.profiler.timer["BCH_Product"] += omp_get_wtime() - tstart;
   return Z;
}















///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// Below is the implementation of the commutators in the various channels
///////////////////////////////////////////////////////////////////////////////////////////

//*****************************************************************************************
//                ____Y    __         ____X
//          X ___(_)             Y___(_) 
//
//  [X1,Y1](0) = Sum_ab (2j_a+1) x_ab y_ba  (n_a-n_b) 
//             = Sum_a  (2j_a+1)  (xy-yx)_aa n_a
//
// -- AGREES WITH NATHAN'S RESULTS
/// \f[
///  [X_{1)},Y_{(1)}]_{(0)} = \sum_{a} n_a (2j_a+1) \left(X_{(1)}Y_{(1)}-Y_{(1)}X_{(1)}\right)_{aa}
/// \f]
//void Operator::comm110ss( const Operator& X, const Operator& Y) 
void comm110ss( const Operator& X, const Operator& Y, Operator& Z) 
{
  if (X.IsHermitian() and Y.IsHermitian()) return ; // I think this is the case
  if (X.IsAntiHermitian() and Y.IsAntiHermitian()) return ; // I think this is the case

   arma::mat xyyx = X.OneBody*Y.OneBody - Y.OneBody*X.OneBody;
   for ( auto& a : Z.modelspace->holes) 
   {
      Orbit& oa = Z.modelspace->GetOrbit(a);
      Z.ZeroBody += (oa.j2+1) * oa.occ * xyyx(a,a);
   }
}


//*****************************************************************************************
//         __Y__       __X__
//        ()_ _()  -  ()_ _()
//           X           Y
//
//  [ X^(2), Y^(2) ]^(0) = 1/2 Sum_abcd  Sum_J (2J+1) x_abcd y_cdab (n_a n_b nbar_c nbar_d)
//                       = 1/2 Sum_J (2J+1) Sum_abcd x_abcd y_cdab (n_a n_b nbar_c nbar_d)  
//                       = 1/2 Sum_J (2J+1) Sum_ab  (X*P_pp*Y)_abab  P_hh
//
//  -- AGREES WITH NATHAN'S RESULTS (within < 1%)
/// \f[
/// [X_{(2)},Y_{(2)}]_{(0)} = \frac{1}{2} \sum_{J} (2J+1) \sum_{abcd} (n_a n_b \bar{n}_c \bar{n}_d) \tilde{X}_{abcd}^{J} \tilde{Y}_{cdab}^{J}
/// \f]
/// may be rewritten as
/// \f[
/// [X_{(2)},Y_{(2)}]_{(0)} = 2 \sum_{J} (2J+1) Tr(X_{hh'pp'}^{J} Y_{pp'hh'}^{J})
/// \f] where we obtain a factor of four from converting two unrestricted sums to restricted sums, i.e. \f$\sum_{ab} \rightarrow \sum_{a\leq b} \f$,
/// and using the normalized TBME.
//void Operator::comm220ss( const Operator& X, const Operator& Y) 
void comm220ss( const Operator& X, const Operator& Y, Operator& Z) 
{
   if (X.IsHermitian() and Y.IsHermitian()) return; // I think this is the case
   if (X.IsAntiHermitian() and Y.IsAntiHermitian()) return; // I think this is the case

   for (int ch=0;ch<Z.nChannels;++ch)
   {
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      auto hh = tbc.GetKetIndex_hh();
      auto ph = tbc.GetKetIndex_ph();
      auto pp = tbc.GetKetIndex_pp();
      arma::uvec nbar_indices = arma::join_cols(hh,ph);
      nbar_indices = arma::join_cols(nbar_indices,pp);
      if (hh.size()==0 ) continue;
      auto nn = tbc.Ket_occ_hh;
      arma::vec nbarnbar = arma::join_cols(tbc.Ket_unocc_hh, tbc.Ket_unocc_ph);
      auto & X2 = X.TwoBody.GetMatrix(ch).submat(hh,nbar_indices);
      arma::mat Y2 = Y.TwoBody.GetMatrix(ch).submat(nbar_indices,hh);
      Y2.head_rows(nbarnbar.size()).each_col() %= nbarnbar;
      Z.ZeroBody += 2 * (2*tbc.J+1) * arma::sum(arma::diagvec( X2 * Y2 ) % nn); // This could be made more efficient, but who cares?
   }
}

//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |
//
// -- AGREES WITH NATHAN'S RESULTS
/// \f[
/// [X_{(1)},Y_{(1)}]_{(1)} = X_{(1)}Y_{(1)} - Y_{(1)}X_{(1)}
/// \f]
//void Operator::comm111ss( Operator & Y, Operator& Z) 
//void Operator::comm111ss( const Operator & X, const Operator& Y) 
void comm111ss( const Operator & X, const Operator& Y, Operator& Z) 
{
   Z.OneBody += X.OneBody*Y.OneBody - Y.OneBody*X.OneBody;
}

//*****************************************************************************************
//                                       |
//      i |              i |             |
//        |    ___.Y       |__X__        |
//        |___(_)    _     |   (_)__.    |  [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab 
//      j | X            j |        Y    |        * sum_J (2J+1) x_biaj^(J)  
//                                       |      
//---------------------------------------*        = 1/(2j+1) sum_a n_a sum_J (2J+1)
//                                                  * sum_b y_ab x_biaj - yba x_aibj
//
//                     (note: I think this should actually be)
//                                                = sum_ab (n_a nbar_b) sum_J (2J+1)/(2j_i+1)
//                                                      * y_ab xbiag - yba x_aibj
//
// -- AGREES WITH NATHAN'S RESULTS 
/// Returns \f$ [X_{(1)},Y_{(2)}] - [Y_{(1)},X_{(2)}] \f$, where
/// \f[
/// [X_{(1)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1}\sum_{ab} (n_a \bar{n}_b) \sum_{J} (2J+1) (X_{ab} Y^J_{biaj} - X_{ba} Y^J_{aibj})
/// \f]
//void Operator::comm121ss( const Operator& X, const Operator& Y) 
void comm121ss( const Operator& X, const Operator& Y, Operator& Z) 
{
   index_t norbits = Z.modelspace->GetNumberOrbits();
   #pragma omp parallel for 
   for (index_t i=0;i<norbits;++i)
   {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      index_t jmin = Z.IsNonHermitian() ? 0 : i;
      for (auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
      {
          if (j<jmin) continue; // only calculate upper triangle
          for (auto& a : Z.modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = Z.modelspace->GetOrbit(a);
             for (index_t b=0; b<norbits; ++b)
             {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ);
                if (std::abs(nanb)<OCC_CUT) continue;
                if (Y.particle_rank>1)
                {
                  Z.OneBody(i,j) += (ob.j2+1) * nanb *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                  Z.OneBody(i,j) -= (oa.j2+1) * nanb *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,j) ;
                }
                // comm211 part
                if (X.particle_rank>1)
                {
                  Z.OneBody(i,j) -= (ob.j2+1) * nanb * Y.OneBody(a,b) * X.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                  Z.OneBody(i,j) += (oa.j2+1) * nanb * Y.OneBody(b,a) * X.TwoBody.GetTBMEmonopole(a,i,b,j) ;
                }
             }
          }
      }
   }
}



//*****************************************************************************************
//
//      i |              i |            [X2,Y2](1)  =  1/(2(2j_i+1)) sum_J (2J+1) 
//        |__Y__           |__X__           * sum_abc (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
//        |    /\          |    /\          * (x_ciab y_abcj - y_ciab xabcj)
//        |   (  )   _     |   (  )                                                                                      
//        |____\/          |____\/       = 1/(2(2j+1)) sum_J (2J+1)
//      j | X            j |  Y            *  sum_c ( Pp*X*Phh*Y*Pp - Pp*Y*Phh*X*Pp)  - (Ph*X*Ppp*Y*Ph - Ph*Y*Ppp*X*Ph)_cicj
//                                     
//
// -- AGREES WITH NATHAN'S RESULTS 
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
// \f[
// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2(2j_i+1)}\sum_{J}(2J+1)\sum_{c}
// \left( \mathcal{P}_{pp} (X \mathcal{P}_{hh} Y^{J} 
// - Y^{J} \mathcal{P}_{hh} X^{J}) \mathcal{P}_{pp}
//  - \mathcal{P}_{hh} (X^{J} \mathcal{P}_{pp} Y^{J} 
//  -  Y^{J} \mathcal{P}_{pp} X^{J}) \mathcal{P}_{hh} \right)_{cicj}
// \f]
/// \f[
/// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2(2j_i+1)}\sum_{J}(2J+1)\sum_{abc} (\bar{n}_a\bar{n}_bn_c + n_an_b\bar{n}_c)
///  (X^{J}_{ciab} Y^{J}_{abcj} - Y^{J}_{ciab}X^{J}_{abcj})
/// \f]
/// This may be rewritten as
/// \f[
/// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1} \sum_{c} \sum_{J} (2J+1) \left( n_c \mathcal{M}^{J}_{pp,icjc} + \bar{n}_c\mathcal{M}^{J}_{hh,icjc} \right)
/// \f]
/// With the intermediate matrix \f[ \mathcal{M}^{J}_{pp} \equiv \frac{1}{2} (X^{J}\mathcal{P}_{pp} Y^{J} - Y^{J}\mathcal{P}_{pp}X^{J}) \f]
/// and likewise for \f$ \mathcal{M}^{J}_{hh} \f$
//void Operator::comm221ss( const Operator& X, const Operator& Y) 
void comm221ss( const Operator& X, const Operator& Y, Operator& Z) 
{

   double t_start = omp_get_wtime();
   int norbits = Z.modelspace->GetNumberOrbits();

   static TwoBodyME Mpp = Y.TwoBody;
   static TwoBodyME Mhh = Y.TwoBody;

   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo.
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

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


      if (Z.IsHermitian())
      {
         Matrixpp +=  Matrixpp.t();
         Matrixhh +=  Matrixhh.t();
      }
      else if (Z.IsAntiHermitian()) // i.e. LHS and RHS are both hermitian or ant-hermitian
      {
         Matrixpp -=  Matrixpp.t();
         Matrixhh -=  Matrixhh.t();
      }
      else
      {
        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
        if (kets_ph.size()>0)
          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
      }


   } //for ch

   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      int jmin = Z.IsNonHermitian() ? 0 : i;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         if (j<jmin) continue;
         double cijJ = 0;
         for (int ch=0;ch<Z.nChannels;++ch)
         {
            TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
            double Jfactor = (2*tbc.J+1.0);
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : Z.modelspace->holes)
            {
               Orbit& oc = Z.modelspace->GetOrbit(c);
               cijJ += Jfactor * oc.occ     * Mpp.GetTBME(ch,c,i,c,j);
               cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,j);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : Z.modelspace->particles)
            {
               cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,j);
            }
         }
         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   }

   X.profiler.timer["comm221ss"] += omp_get_wtime() - t_start;

}





//*****************************************************************************************
//
//    |     |               |      |           [X2,Y1](2) = sum_a ( Y_ia X_ajkl + Y_ja X_iakl - Y_ak X_ijal - Y_al X_ijka )
//    |     |___.Y          |__X___|         
//    |     |         _     |      |          
//    |_____|               |      |_____.Y        
//    |  X  |               |      |            
//
// -- AGREES WITH NATHAN'S RESULTS
/// Returns \f$ [X_{(1)},Y_{(2)}]_{(2)} - [Y_{(1)},X_{(2)}]_{(2)} \f$, where
/// \f[
/// [X_{(1)},Y_{(2)}]^{J}_{ijkl} = \sum_{a} ( X_{ia}Y^{J}_{ajkl} + X_{ja}Y^{J}_{iakl} - X_{ak} Y^{J}_{ijal} - X_{al} Y^{J}_{ijka} )
/// \f]
/// here, all TBME are unnormalized, i.e. they should have a tilde.
// This is still too slow...
//void Operator::comm122ss( Operator& Y, Operator& Z ) 
//void Operator::comm122ss( const Operator& X, const Operator& Y ) 
void comm122ss( const Operator& X, const Operator& Y, Operator& Z ) 
{
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;
   int hZ = Z.IsHermitian() ? 1 : -1;

   int n_nonzero = Z.modelspace->SortedTwoBodyChannels.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      auto& X2 = X.TwoBody.GetMatrix(ch,ch);
      auto& Y2 = Y.TwoBody.GetMatrix(ch,ch);
      auto& Z2 = Z.TwoBody.GetMatrix(ch,ch);
      arma::mat W2(size(Z2),arma::fill::zeros); // temporary intermediate matrix

      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = Z.modelspace->GetOrbit(i);
         Orbit& oj = Z.modelspace->GetOrbit(j);
         int flipphaseij = - Z.modelspace->phase((oi.j2+oj.j2)/2-tbc.J);

         // make lists of the indices we want, then do matrix multiplication.
         // there may be a more efficient way to find these
         std::vector<index_t> ind1_ia, ind1_ja,ind2_aj,ind2_ai;
         std::vector<double> factor_ia,factor_ja;
         for (int a : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
         {
            size_t ind2 = tbc.GetLocalIndex( std::min(a,j), std::max(a,j) );
            if (ind2<0 or ind2>=tbc.GetNumberKets()) continue;
            ind1_ia.push_back(a);
            ind2_aj.push_back(ind2);
            factor_ia.push_back( a>j ? flipphaseij : (a==j ? SQRT2 : 1));
         }
         if (i!=j)
         {
           for (int a : Z.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
           {
              size_t ind2 = tbc.GetLocalIndex( std::min(a,i), std::max(a,i) );
              if (ind2<0 or ind2>=tbc.GetNumberKets()) continue;
              ind1_ja.push_back(a);
              ind2_ai.push_back(ind2);
              factor_ja.push_back( i>a ? flipphaseij : (i==a ? SQRT2 : 1));
           }

         }

         arma::uvec u_ind1_ia(ind1_ia);
         arma::uvec u_ind1_ja(ind1_ja);
         arma::uvec u_ind2_aj(ind2_aj);
         arma::uvec u_ind2_ai(ind2_ai);
         arma::vec  v_factor_ia(factor_ia);
         arma::vec  v_factor_ja(factor_ja);
         if (i==j)
         {
           v_factor_ia /= SQRT2;
           v_factor_ja /= SQRT2;
         }

         // This is fairly obfuscated, but hopefully faster for bigger calculations
         if (X.particle_rank>1 and Y.particle_rank>1)
         {
            W2.col(indx_ij) = join_horiz(    Y2.cols(join_vert( u_ind2_aj,u_ind2_ai))  , X2.cols(join_vert(u_ind2_aj,u_ind2_ai) ) )
                             * join_vert(  join_vert( X1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, X1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja ),
                                          -join_vert( Y1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, Y1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja ));
         }
         else if (X.particle_rank<2 and Y.particle_rank>1)
         {
            W2.col(indx_ij) =     Y2.cols(join_vert( u_ind2_aj,u_ind2_ai))    
                             *   join_vert( X1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, X1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja );
         }
         else if (X.particle_rank>1 and Y.particle_rank<2)
         {
            W2.col(indx_ij) =      -X2.cols(join_vert(u_ind2_aj,u_ind2_ai) ) 
                                      *    join_vert( Y1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, Y1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja );
         }

      if (i==j) W2.col(indx_ij) *= 2;

      }
      Z2 -= W2 + hZ*W2.t();
   }

}





//*****************************************************************************************
//
//  |     |      |     |   
//  |__Y__|      |__x__|   [X2,Y2](2)_pp(hh) = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
//  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
//  |__X__|      |__Y__|   
//  |     |      |     |   
//
// -- AGREES WITH NATHAN'S RESULTS
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
/// Calculates the part of the commutator \f$ [X_{(2)},Y_{(2)}]_{(2)} \f$ which involves particle-particle
/// or hole-hole intermediate states.
/// \f[
/// [X_{(2)},Y_{(2)}]^{J}_{ijkl} = \frac{1}{2} \sum_{ab} (\bar{n}_a\bar{n}_b - n_an_b) (X^{J}_{ijab}Y^{J}_{ablk} - Y^{J}_{ijab}X^{J}_{abkl})
/// \f]
/// This may be written as
/// \f[
/// [X_{(2)},Y_{(2)}]^{J} = \mathcal{M}^{J}_{pp} - \mathcal{M}^{J}_{hh}
/// \f]
/// With the intermediate matrices
/// \f[
/// \mathcal{M}^{J}_{pp} \equiv \frac{1}{2}(X^{J} \mathcal{P}_{pp} Y^{J} - Y^{J} \mathcal{P}_{pp} X^{J})
/// \f]
/// and likewise for \f$ \mathcal{M}^{J}_{hh} \f$.
//void Operator::comm222_pp_hhss( Operator& opright, Operator& opout ) 
//void Operator::comm222_pp_hhss( const Operator& X, const Operator& Y ) 
void comm222_pp_hhss( const Operator& X, const Operator& Y, Operator& Z ) 
{

   static TwoBodyME Mpp = Z.TwoBody;
   static TwoBodyME Mhh = Z.TwoBody;

   double t = omp_get_wtime();
   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo.
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

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


      if (Z.IsHermitian())
      {
         Matrixpp +=  Matrixpp.t();
         Matrixhh +=  Matrixhh.t();
      }
      else if (Z.IsAntiHermitian()) // i.e. LHS and RHS are both hermitian or ant-hermitian
      {
         Matrixpp -=  Matrixpp.t();
         Matrixhh -=  Matrixhh.t();
      }
      else
      {
        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
        if (kets_ph.size()>0)
          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
      }


      // The two body part
      OUT += Matrixpp - Matrixhh;
   } //for ch
   X.profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;
}








void ConstructScalarMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, TwoBodyME& Mpp, TwoBodyME& Mhh)
{
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
   bool z_is_hermitian = Z.IsHermitian();
   bool z_is_antihermitian = Z.IsAntiHermitian();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

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


      if (z_is_hermitian)
      {
         Matrixpp +=  Matrixpp.t();
         Matrixhh +=  Matrixhh.t();
      }
      else if (z_is_antihermitian) // i.e. LHS and RHS are both hermitian or ant-hermitian
      {
         Matrixpp -=  Matrixpp.t();
         Matrixhh -=  Matrixhh.t();
      }
      else
      {
        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
        if (kets_ph.size()>0)
          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
      }


   } //for ch

}


/// Since comm222_pp_hhss() and comm221ss() both require the construction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
void comm222_pp_hh_221ss( const Operator& X, const Operator& Y, Operator& Z )  
{

//   int herm = Z.IsHermitian() ? 1 : -1;
//   Operator& Z = *this;
   int norbits = Z.modelspace->GetNumberOrbits();

   static TwoBodyME Mpp = Z.TwoBody;
   static TwoBodyME Mhh = Z.TwoBody;

   double t = omp_get_wtime();
   ConstructScalarMpp_Mhh( X, Y, Z, Mpp, Mhh);

//   Z.TwoBody += (Mpp - Mhh);
   Z.TwoBody += Mpp;
   Z.TwoBody -= Mhh;
//   OUT += Matrixpp - Matrixhh;
   X.profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;

   t = omp_get_wtime();
   // The one body part
   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      int jmin = Z.IsNonHermitian() ? 0 : i;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         if (j<jmin) continue;
         double cijJ = 0;
         for (int ch=0;ch<Z.nChannels;++ch)
         {
            TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
            double Jfactor = (2*tbc.J+1.0);
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : Z.modelspace->holes)
            {
               Orbit& oc = Z.modelspace->GetOrbit(c);
               cijJ += Jfactor * oc.occ * Mpp.GetTBME(ch,c,i,c,j); 
               cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,j); 
            }
            // Sum c over particles and include the n_a * n_b terms
            for (auto& c : Z.modelspace->particles)
            {
               cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,j);
            }
         }
         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   } // for i
   X.profiler.timer["pphh One Body bit"] += omp_get_wtime() - t;
}



//**************************************************************************
//
//  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
//                        { k l J'}
// SCALAR VARIETY
/// The scalar Pandya transformation is defined as
/// \f[
///  \bar{X}^{J}_{i\bar{j}k\bar{l}} = - \sum_{J'} (2J'+1)
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  X^{J}_{ilkj}
/// \f]
/// where the overbar indicates time-reversed orbits.
//void Operator::DoPandyaTransformation_SingleChannel(arma::mat& TwoBody_CC_ph, int ch_cc, std::string orientation="normal") const
void DoPandyaTransformation_SingleChannel(const Operator& Z, arma::mat& TwoBody_CC_ph, int ch_cc, std::string orientation="normal")
{
   int herm = Z.IsHermitian() ? 1 : -1;
   TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
   int nKets_cc = tbc_cc.GetNumberKets();
   arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
   int nph_kets = kets_ph.n_rows;
   int J_cc = tbc_cc.J;

   if (orientation=="normal") TwoBody_CC_ph.zeros( 2*nph_kets, nKets_cc);
   else if (orientation=="transpose") TwoBody_CC_ph.zeros( nKets_cc, 2*nph_kets);
   else
   {
     std::cout << __PRETTY_FUNCTION__ << " =>  Unknown orientation input  " << orientation << ". Don't know what to do with this." << std::endl;
     return;
   }

   // loop over cross-coupled ph bras <ab| in this channel
   // (this is the side that gets summed over in the matrix multiplication)
   for (int ibra=0; ibra<nph_kets; ++ibra)
   {
      Ket & bra_cc = tbc_cc.GetKet( kets_ph[ibra] );
      // we want to evaluate a<=b and a>=b, so to avoid code duplication, we turn this into a loop over the two orderings
      std::vector<size_t> ab_switcheroo = { bra_cc.p, bra_cc.q };
      for ( int ab_case=0; ab_case<=1; ab_case++)
      {
        int a = ab_switcheroo[ab_case];   // this little bit gives us a,b if ab_case=0 and b,a if ab_case=1
        int b = ab_switcheroo[1-ab_case];
        size_t bra_shift = ab_case*nph_kets;  // if we switch a<->b, we offset the bra index by nph_kets
//        int a = bra_cc.p;
//        int b = bra_cc.q;

        Orbit & oa = Z.modelspace->GetOrbit(a);
        Orbit & ob = Z.modelspace->GetOrbit(b);
        double ja = oa.j2*0.5;
        double jb = ob.j2*0.5;
        double na_nb_factor = oa.occ - ob.occ;

        // loop over cross-coupled kets |cd> in this channel
        for (int iket_cc=0; iket_cc<nKets_cc; ++iket_cc)
        {
           Ket & ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
           int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
           int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
           Orbit & oc = Z.modelspace->GetOrbit(c);
           Orbit & od = Z.modelspace->GetOrbit(d);
           double jc = oc.j2*0.5;
           double jd = od.j2*0.5;


           int jmin = std::max(std::abs(ja-jd),std::abs(jc-jb));
           int jmax = std::min(ja+jd,jc+jb);
           double Xbar = 0;
           for (int J_std=jmin; J_std<=jmax; ++J_std)
           {
              double sixj = Z.modelspace->GetSixJ(ja,jb,J_cc,jc,jd,J_std);
              if (std::abs(sixj) < 1e-8) continue;
              double tbme = Z.TwoBody.GetTBME_J(J_std,a,d,c,b);
              Xbar -= (2*J_std+1) * sixj * tbme ;
           }
           if (orientation=="normal")
           {
             TwoBody_CC_ph(ibra+bra_shift, iket_cc) = Xbar;
           }
           else // "transpose"
           {
             TwoBody_CC_ph(iket_cc,ibra+bra_shift) = herm * Xbar * na_nb_factor;  // we slap the (na-nb) on the transposed one.
           }

//           // Exchange (a <-> b) to account for the (n_a - n_b) term
//           jmin = std::max(std::abs(jb-jd),std::abs(jc-ja));
//           jmax = std::min(jb+jd,jc+ja);
//           Xbar = 0;
//           for (int J_std=jmin; J_std<=jmax; ++J_std)
//           {
//              double sixj = Z.modelspace->GetSixJ(jb,ja,J_cc,jc,jd,J_std);
//              if (std::abs(sixj) < 1e-8) continue;
//              double tbme = Z.TwoBody.GetTBME_J(J_std,b,d,c,a);
//              Xbar -= (2*J_std+1) * sixj * tbme ;
//           }
//           if (orientation=="normal")
//           {
//             TwoBody_CC_ph(ibra+nph_kets,iket_cc) = Xbar;
//           }
//           else  // "transpose"
//           {
//             TwoBody_CC_ph(iket_cc,ibra+nph_kets) = herm * Xbar * -na_nb_factor;
//           }

        }
      }
   }
}


//void Operator::DoPandyaTransformation(deque<arma::mat>& TwoBody_CC_ph, std::string orientation="normal") const
void DoPandyaTransformation(const Operator& Z, std::deque<arma::mat>& TwoBody_CC_ph, std::string orientation="normal")
{
   // loop over cross-coupled channels
   int n_nonzero = Z.modelspace->SortedTwoBodyChannels_CC.size();
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch_cc = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      DoPandyaTransformation_SingleChannel(Z, TwoBody_CC_ph[ch_cc] , ch_cc, orientation);
   }
}



void AddInversePandyaTransformation_SingleChannel( Operator& Z,  arma::mat& Zbar, int ch_cc)
{
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
   TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
   int Jcc = tbc_cc.J;
   int Jhat2 = 2*Jcc+1;
   int parity_cc = tbc_cc.parity;
   int Tz_cc = tbc_cc.Tz;
   int nkets_cc = tbc_cc.GetNumberKets();
   int n_nonzeroChannels = Z.modelspace->SortedTwoBodyChannels.size();
//   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nKets = tbc.GetNumberKets();
      auto& Zmat = Z.TwoBody.GetMatrix(ch,ch);

      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = Z.modelspace->GetOrbit(i);
         Orbit & oj = Z.modelspace->GetOrbit(j);
         double ji = 0.5*oi.j2;
         double jj = 0.5*oj.j2;
         int ketmin = Z.IsHermitian() ? ibra : ibra+1;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = Z.modelspace->GetOrbit(k);
            Orbit & ol = Z.modelspace->GetOrbit(l);
            double jk = 0.5*ok.j2;
            double jl = 0.5*ol.j2;

            double commij = 0;
            double commji = 0;

            int jmin = std::max(std::abs(int(ji-jl)),std::abs(int(jk-jj)));
            int jmax = std::min(int(ji+jl),int(jk+jj));
            if ( ((oi.l+ol.l)%2==parity_cc)  and  (std::abs(oi.tz2+ol.tz2)==Tz_cc*2) and (Jcc>=jmin) and (Jcc<=jmax) )
            {
               double sixj = Z.modelspace->GetSixJ(ji,jj,J,jk,jl,Jcc);
               int indx_il = tbc_cc.GetLocalIndex(i,l) ;
               int indx_kj = tbc_cc.GetLocalIndex( std::min(j,k), std::max(j,k) ) +(k>j?nkets_cc:0);
               commij -= Jhat2 * sixj * Zbar(indx_il,indx_kj) ;

            }

            if (k==l)
            {
              commji = commij;
            }
            else if (i==j)
            {
              commji = Z.modelspace->phase(ji+jj+jk+jl) * commij;
            }
            else
            {
              // now loop over the cross coupled TBME's
              jmin = std::max(std::abs(int(jj-jl)),std::abs(int(jk-ji)));
              jmax = std::min(int(jj+jl),int(jk+ji));
              if ( (oi.l+ok.l)%2==parity_cc  and  std::abs(oi.tz2+ok.tz2)==Tz_cc*2 and Jcc>=jmin and Jcc<=jmax)
              {
                 double sixj = Z.modelspace->GetSixJ(jj,ji,J,jk,jl,Jcc);
                 int indx_ik = tbc_cc.GetLocalIndex(i,k) ;
                 int indx_lj = tbc_cc.GetLocalIndex(std::min(l,j),std::max(l,j)) +(l>j?nkets_cc:0);
                 // we always have i<=k so we should always flip Z_jlki = (-1)^{i+j+k+l} Z_iklj
                 // the phase we get from that flip combines with the phase from Pij, to give the phase included below
                 commji -= Jhat2 *  sixj *  Zbar(indx_ik, indx_lj) ;


              }
            }
            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            #pragma omp atomic
            Zmat(ibra,iket) -= (commij - Z.modelspace->phase(jk+jl-J ) * commji) / norm;
         }
      }
   }
}



//void Operator::AddInversePandyaTransformation(const deque<arma::mat>& Zbar)
void AddInversePandyaTransformation(const std::deque<arma::mat>& Zbar, Operator& Z)
{
    // Do the inverse Pandya transform
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
   int n_nonzeroChannels = Z.modelspace->SortedTwoBodyChannels.size();

   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nKets = tbc.GetNumberKets();

      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = Z.modelspace->GetOrbit(i);
         Orbit & oj = Z.modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;
         int ketmin = Z.IsHermitian() ? ibra : ibra+1;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = Z.modelspace->GetOrbit(k);
            Orbit & ol = Z.modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;

            double commij = 0;
            double commji = 0;

            int parity_cc = (oi.l+ol.l)%2;
            int Tz_cc = std::abs(oi.tz2+ol.tz2)/2;
            int Jpmin = std::max(std::abs(int(ji-jl)),std::abs(int(jk-jj)));
            int Jpmax = std::min(int(ji+jl),int(jk+jj));
            for (int Jprime=Jpmin; Jprime<=Jpmax; ++Jprime)
            {
               double sixj = Z.modelspace->GetSixJ(ji,jj,J,jk,jl,Jprime);
               if (std::abs(sixj)<1e-8) continue;
               int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
               int nkets_cc = tbc_cc.GetNumberKets();
               int indx_il = tbc_cc.GetLocalIndex(std::min(i,l),std::max(i,l)) +(i>l?nkets_cc:0);
               int indx_kj = tbc_cc.GetLocalIndex(std::min(j,k),std::max(j,k)) +(k>j?nkets_cc:0);
               double me1 = Zbar.at(ch_cc)(indx_il,indx_kj);
               commij -= (2*Jprime+1) * sixj * me1;
            }

            if (k==l)
            {
              commji = commij;
            }
            else if (i==j)
            {
              commji = Z.modelspace->phase(ji+jj+jk+jl) * commij;
            }
            else
            {
              // now loop over the cross coupled TBME's
              parity_cc = (oi.l+ok.l)%2;
              Tz_cc = std::abs(oi.tz2+ok.tz2)/2;
              Jpmin = std::max(std::abs(int(jj-jl)),std::abs(int(jk-ji)));
              Jpmax = std::min(int(jj+jl),int(jk+ji));
              for (int Jprime=Jpmin; Jprime<=Jpmax; ++Jprime)
              {
                 double sixj = Z.modelspace->GetSixJ(jj,ji,J,jk,jl,Jprime);
                 if (std::abs(sixj)<1e-8) continue;
                 int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
                 TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
                 int nkets_cc = tbc_cc.GetNumberKets();
                 int indx_ik = tbc_cc.GetLocalIndex(std::min(i,k),std::max(i,k)) +(i>k?nkets_cc:0);
                 int indx_lj = tbc_cc.GetLocalIndex(std::min(l,j),std::max(l,j)) +(l>j?nkets_cc:0);
                 // we always have i<=k so we should always flip Z_jlki = (-1)^{i+j+k+l} Z_iklj
                 // the phase we get from that flip combines with the phase from Pij, to give the phase included below
                 double me1 = Zbar.at(ch_cc)(indx_ik, indx_lj) ;
                 commji -= (2*Jprime+1) *  sixj * me1;

              }
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            Z.TwoBody.GetMatrix(ch,ch)(ibra,iket) -= (commij - Z.modelspace->phase(jk+jl-J ) * commji) / norm;
         }
      }
   }
}

///*************************************
/// convenience function
/// called by comm222_phss
///*************************************
//deque<arma::mat> Operator::InitializePandya(size_t nch, std::string orientation="normal")
std::deque<arma::mat> InitializePandya(Operator& Z, size_t nch, std::string orientation="normal")
{
   std::deque<arma::mat> X(nch);
   int n_nonzero = Z.modelspace->SortedTwoBodyChannels_CC.size();
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch_cc = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();
      if (orientation=="normal")
         X[ch_cc] = arma::mat(2*nph_kets,   nKets_cc, arma::fill::zeros);
      else if (orientation=="transpose")
         X[ch_cc] = arma::mat(nKets_cc, 2*nph_kets,   arma::fill::zeros);
   }
   return X;
}

//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.     
//                                             
//   |          |      |          |           
//   |     __Y__|      |     __X__|            
//   |    /\    |      |    /\    |
//   |   (  )   |  _   |   (  )   |            
//   |____\/    |      |____\/    |            
//   |  X       |      |  Y       |            
//           
//            
// -- This appears to agree with Nathan's results
//
/// Calculates the part of \f$ [X_{(2)},Y_{(2)}]_{ijkl} \f$ which involves ph intermediate states, here indicated by \f$ Z^{J}_{ijkl} \f$
/// \f[
/// Z^{J}_{ijkl} = \sum_{ab}(n_a\bar{n}_b-\bar{n}_an_b)\sum_{J'} (2J'+1)
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
/// \left( \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} - 
///   \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \right)
///  -(-1)^{j_i+j_j-J}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_i  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
/// \left( \bar{X}^{J'}_{j\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{i}} - 
///   \bar{Y}^{J'}_{j\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{i}} \right)
/// \right]
/// \f]
/// This is implemented by defining an intermediate matrix
/// \f[
/// \bar{Z}^{J}_{i\bar{l}k\bar{j}} \equiv \sum_{ab}(n_a\bar{n}_b)
/// \left[ \left( \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} - 
///   \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \right)
/// -\left( \bar{X}^{J'}_{i\bar{l}b\bar{a}}\bar{Y}^{J'}_{b\bar{a}k\bar{j}} - 
///    \bar{Y}^{J'}_{i\bar{l}b\bar{a}}\bar{X}^{J'}_{b\bar{a}k\bar{j}} \right)\right]
/// \f]
/// The Pandya-transformed matrix elements are obtained with DoPandyaTransformation().
/// The matrices \f$ \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} \f$
/// and \f$ \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \f$
/// are related by a Hermitian conjugation, which saves two matrix multiplications.
/// The commutator is then given by
/// \f[
/// Z^{J}_{ijkl} = \sum_{J'} (2J'+1)
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  \bar{Z}^{J'}_{i\bar{l}k\bar{j}}
///  -(-1)^{j_i+j_j-J}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_i  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  \bar{Z}^{J'}_{j\bar{l}k\bar{i}}
///  \right]
///  \f]
///
//void Operator::comm222_phss( const Operator& X, const Operator& Y ) 
void comm222_phss( const Operator& X, const Operator& Y, Operator& Z ) 
{

   int hy = Y.IsHermitian() ? 1 : -1;
   // Create Pandya-transformed hp and ph matrix elements
   double t_start = omp_get_wtime();

   // Construct the intermediate matrix Z_bar
   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup(Z.GetJRank(), Z.GetTRank(), Z.GetParity());
   size_t nch = Z.modelspace->SortedTwoBodyChannels_CC.size();
   t_start = omp_get_wtime();
   std::deque<arma::mat> Z_bar ( Z.nChannels );
   std::vector<bool> lookup_empty(Z.nChannels,true);
   for (size_t ich=0;ich<nch;++ich)
   {
      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      index_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros( nKets_cc, 2*nKets_cc );
      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
   }


   #ifndef OPENBLAS_NOUSEOMP
//   #pragma omp parallel for schedule(dynamic,1)
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   #endif
   for (size_t ich=0; ich<nch; ++ich )
   {
      if (lookup_empty.at(ich)) continue;
      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC.at(ich);
      const TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      size_t nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

      arma::mat Y_bar_ph;
      arma::mat Xt_bar_ph;

      DoPandyaTransformation_SingleChannel(Y,Y_bar_ph,ch,"normal");
      DoPandyaTransformation_SingleChannel(X,Xt_bar_ph,ch,"transpose");
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
         if ( Z.modelspace->phase( (ket.op->j2 + ket.oq->j2)/2 ) > 0) continue;
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
      if ( Z.IsHermitian() )
      {
         Zbar_ch.head_cols(nKets_cc) += Zbar_ch.head_cols(nKets_cc).t();
      }
      else
      {
         Zbar_ch.head_cols(nKets_cc) -= Zbar_ch.head_cols(nKets_cc).t();
      }
      Zbar_ch.tail_cols(nKets_cc) += Zbar_ch.tail_cols(nKets_cc).t()%PhaseMat;

   }

   X.profiler.timer["Build Z_bar"] += omp_get_wtime() - t_start;

   // Perform inverse Pandya transform on Z_bar to get Z
   t_start = omp_get_wtime();

   // Actually, the Pandya transform has a minus sign in the definition,
   // and the ph commutator has an overall minus sign, so we're technically subtracting
   // the inverse Pandya transformation. Also, the inverse Pandya transformation
   // is just the regular Pandya transformation. The distinction in the code
   // is because some other commutator-specific things are done at the same time.
   AddInversePandyaTransformation(Z_bar, Z);

   Z.modelspace->scalar_transform_first_pass = false;
   X.profiler.timer["InversePandyaTransformation"] += omp_get_wtime() - t_start;

}




//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
////////////   BEGIN SCALAR-SCALAR COMMUTATORS WITH 3-body ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


//  For the three-body operators, we often encounter formulas involving the J-coupled
//  permutation operator P(ij/k)^{J1,J}, which is defined as
//
//  P(ij/k)^{J1,J} = 1 - sum_J2 sqrt{(2J1+1)(2J2+1)} (-1)^2(ji+jj+jk)   {ji jk J1} P_ik
//                                                                      {jk J  J2}
//
//                     - sum_J2 sqrt{(2J1+1)(2J2+1)} (-1)^{jj+jk+J1+J2) {jk ji J1} P_jk
//                                                                      {jk J  J2}
//



//*****************************************************************************************
//
//    ~~~~~~~~~~~~~~~~        Uncoupled expression: 
//   /\      /\      /\         Z_0 = sum_abcdef (nanbnc n`dn`en`f) (X_abcdef Y_defabc - Y_abcdef X_defabc)
// a(  )d  b(  )e  c(  )f
//   \/      \/      \/      Coupled expression:
//    ~~~~~~~~~~~~~~~~        Z_0 = sum_abcdef sum_J1,J2,J  (nanbnc n`dn`en`f)  (2J+1)
//                                    (X^{J1J2J}_abcdef Y^{J1J2J}_defabc - Y^{J1J2J}_abcdef X^{J1J2J}_defabc)
//
void comm330ss( const Operator& X, const Operator& Y, Operator& Z )
{
  std::cout << "start comm330" << std::endl;
  int E3cut = 7;
  double z0 = 0;
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  std::cout << "Norm X3 " << X3.Norm() << "  Norm Y3 " << Y3.Norm() << std::endl;
  if (X3.Norm()<1e-6 or Y3.Norm()<1e-6 ) return;
  #pragma omp parallel for schedule(dynamic,1) reduction(+ : z0)
//  for ( auto a : Z.modelspace->holes )
  for ( size_t indxa=0; indxa<Z.modelspace->holes.size(); indxa++ )
  {
   auto a = Z.modelspace->holes[indxa];
   Orbit& oa = Z.modelspace->GetOrbit(a);
   for( auto b : Z.modelspace->holes )
   {
    Orbit& ob = Z.modelspace->GetOrbit(b);
    int Jab_min = std::abs(oa.j2-ob.j2)/2;
    int Jab_max = std::abs(oa.j2+ob.j2)/2;
    for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
    {
     if (a==b and Jab%2>0) continue;
     for ( auto c : Z.modelspace->holes )
     {
      Orbit& oc = Z.modelspace->GetOrbit(c);
      for ( auto d : Z.modelspace->particles )
      {
       Orbit& od = Z.modelspace->GetOrbit(d);
       for ( auto e : Z.modelspace->particles )
       {
        Orbit& oe = Z.modelspace->GetOrbit(e);
        int Jde_min = std::abs(od.j2-oe.j2)/2;
        int Jde_max = std::abs(od.j2+oe.j2)/2;
        for ( int Jde=Jde_min; Jde<=Jde_max; Jde++)
        {
         if (d==e and Jde%2>0) continue;
         for ( auto f : Z.modelspace->particles )
         {
          Orbit& of = Z.modelspace->GetOrbit(f);
//          std::cout << "abcdef = " << a << " " << b << " " << c << " " << d << " " << e << " " << f
//                    << "   " << Jab << " " << Jde << "  " << "  " << 2*(od.n+oe.n+of.n)+od.l+oe.l+of.l  << std::endl;
          if ( ( 2*(od.n+oe.n+of.n)+od.l+oe.l+of.l) > E3cut ) continue;
          if ( (oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0 ) continue;
          if ( (oa.tz2+ob.tz2+oc.tz2) != (od.tz2+oe.tz2+of.tz2) ) continue;
          int twoJmin = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jde-of.j2) );
          int twoJmax = std::min( 2*Jab+oc.j2, 2*Jde+of.j2 );
          for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
          {
            z0 += (twoJ+1) * ( X3.GetME_pn( Jab, Jde, twoJ, a,b,c,d,e,f) * Y3.GetME_pn( Jde, Jab, twoJ, d,e,f,a,b,c)
                             - Y3.GetME_pn( Jab, Jde, twoJ, a,b,c,d,e,f) * X3.GetME_pn( Jde, Jab, twoJ, d,e,f,a,b,c) );
//	    std::cout << "  " <<  X3.GetME_pn( Jab, Jde, twoJ, a,b,c,d,e,f) << "  " << Y3.GetME_pn( Jde, Jab, twoJ, d,e,f,a,b,c)
//	              << "  " <<  Y3.GetME_pn( Jab, Jde, twoJ, a,b,c,d,e,f) << "  " << X3.GetME_pn( Jde, Jab, twoJ, d,e,f,a,b,c) << std::endl;


          }
         }
        }
       }
      }
     }
    }
   }
  }
  std::cout << "Adding " << z0 << "  to Zero Body" << std::endl;
  Z.ZeroBody += z0;
}


//*****************************************************************************************
//                   |i
//    ~~~~~~~~~~~~~~~|        Uncoupled expression: 
//   /\      /\      |          Z_ij = 1/4 sum_abcde (nanb n`c n`dn`e) (X_abicde Y_cdeabj - Y_abicde X_cdeabj)
// a(  )c  b(  )d    |e 
//   \/      \/      |       Coupled expression:
//    ~~~~~~~~~~~~~~~|        Z_ij = 1/4 sum_abcde sum_J1,J2,J  (nanb n`c n`dn`e)  (2J+1)/(2ji+1)
//                   |j               (X^{J1J2J}_abicde Y^{J1J2J}_cdeabj - Y^{J1J2J}_abicde X^{J1J2J}_cdeabj)
//
void comm331ss( const Operator& X, const Operator& Y, Operator& Z )
{
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z1 = Z.OneBody;

  int norb = Z.modelspace->GetNumberOrbits();
  for (int i=0; i<norb; i++)
  {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    for ( auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      // TODO: At some point, use the hermiticity of Z and only do j<i
      double zij=0;
      for (auto a : Z.modelspace->holes ) // TODO: We can probably speed this up by looping over 2-body channels
      {
        Orbit& oa = Z.modelspace->GetOrbit(a);
        for (auto b : Z.modelspace->holes )
        {
          Orbit& ob = Z.modelspace->GetOrbit(b);
          int Jab_min = std::abs( oa.j2-ob.j2)/2;
          int Jab_max = ( oa.j2+ob.j2)/2;
          for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
          {
            for (auto c : Z.modelspace->particles)
            {
              Orbit& oc = Z.modelspace->GetOrbit(c);
              for (auto d : Z.modelspace->particles)
              {
                Orbit& od = Z.modelspace->GetOrbit(d);
                int Jcd_min = std::abs( oc.j2-od.j2)/2;
                int Jcd_max = ( oc.j2+od.j2)/2;
                for (int Jcd=Jcd_min; Jcd<=Jcd_max; Jcd++)
                {
                  for (auto e : Z.modelspace->particles)
                  {
                    Orbit& oe = Z.modelspace->GetOrbit(e);
                    int twoJ_min = std::max( std::abs(oi.j2-2*Jab), std::abs(oe.j2-2*Jcd) );
                    int twoJ_max = std::min( oi.j2+2*Jab, oe.j2+2*Jcd );
                    for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ++)
                    {
                      double prefactor = 0.25 * (twoJ+1)*(oi.j2+1) * oa.occ*ob.occ*(1-oc.occ)*(1-od.occ)*(1-oe.occ);
                      zij += prefactor * X3.GetME_pn( Jab, Jcd, twoJ, a,b,i,c,d,e) * Y3.GetME_pn( Jcd, Jab, twoJ, c,d,e,a,b,j);
                      zij -= prefactor * Y3.GetME_pn( Jab, Jcd, twoJ, a,b,i,c,d,e) * X3.GetME_pn( Jcd, Jab, twoJ, c,d,e,a,b,j);
                    }
                  }
                  for (auto e : Z.modelspace->holes) // same thing, just different occupation number factors.
                  {
                    Orbit& oe = Z.modelspace->GetOrbit(e);
                    int twoJ_min = std::max( std::abs(oi.j2-2*Jcd), std::abs(oe.j2-2*Jab) );
                    int twoJ_max = std::min( oi.j2+2*Jcd, oe.j2+2*Jab );
                    for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ++)
                    {
                      double prefactor = 0.25 * (twoJ+1)*(oi.j2+1) * oa.occ*ob.occ*oe.occ*(1-oc.occ)*(1-od.occ);
                      zij += prefactor * X3.GetME_pn( Jcd, Jab, twoJ, c,d,i,a,b,e) * Y3.GetME_pn( Jab, Jcd, twoJ, a,b,e,c,d,j);
                      zij -= prefactor * Y3.GetME_pn( Jcd, Jab, twoJ, c,d,i,a,b,e) * X3.GetME_pn( Jab, Jcd, twoJ, a,b,e,c,d,j);
                    }
                  }// for e
                }// for Jcd
              }// for d
            }// for c
          }// for Jab
        }// for b
      }// for a
      Z1(i,j) += zij;
    }// for j
  }// for i

}


//*****************************************************************************************
//                   |i
//    ~~~~~~~~       |        Uncoupled expression: 
//   /\      /\      |          Z_ij = 1/4 sum_abcd (nanb n`cn`d) (X_abcd Y_cdiabj - Y_abicdj X_cdab)
// a(  )c  b(  )d    | 
//   \/      \/      |       Coupled expression:
//    ~~~~~~~~~~~~~~~|        Z_ij = 1/4sum_abcd sum_J1J  (nanb n`c n`d)  (2J+1)/(2ji+1)
//                   |j               (X^{J}_abcd Y^{J1J1J}_cdiabj - Y^{J1J1J}_abicdj X^{J}_cdab)
//
void comm231ss( const Operator& X, const Operator& Y, Operator& Z )
{
  auto& X2 = X.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y2 = Y.TwoBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z1 = Z.OneBody;

  int norb = Z.modelspace->GetNumberOrbits();
  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  for (int i=0; i<norb; i++)
  {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    for ( auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      double zij=0;
      for (int ch=0; ch<nch; ch++)
      {
        auto tbc = Z.modelspace->GetTwoBodyChannel(ch);
        int J = tbc.J;
        for ( auto ibra : tbc.KetIndex_hh )
        {
          Ket& bra = tbc.GetKet(ibra);
          int a = bra.p;
          int b = bra.q;
          for ( auto iket : tbc.KetIndex_pp )
          {
            Ket& ket = tbc.GetKet(iket);
            int c = ket.p;
            int d = ket.q;
            double Xabcd = X2.GetTBME(ch,bra,ket);
            double Yabcd = Y2.GetTBME(ch,bra,ket);
            double Xcdab = X2.GetTBME(ch,ket,bra);
            double Ycdab = Y2.GetTBME(ch,ket,bra);
            double occfactor = (bra.op->occ)*(bra.oq->occ)*(1-ket.op->occ)*(1-ket.oq->occ);
            int twoJ_min = std::abs( 2*J - oi.j2);
            int twoJ_max = 2*J + oi.j2;
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; ++twoJ)
            {
              zij += occfactor* (twoJ+1)* ( Xabcd * Y3.GetME_pn(twoJ,J,J,c,d,i,a,b,j) - Y3.GetME_pn(twoJ,J,J,a,b,i,c,d,j) * Xcdab
                                           -Yabcd * X3.GetME_pn(twoJ,J,J,c,d,i,a,b,j) + X3.GetME_pn(twoJ,J,J,a,b,i,c,d,j) * Ycdab );
            }
          }
        }
      }
      Z1(i,j) += zij / (oi.j2+1.0);
    }// for j
  }// for i

}





//*****************************************************************************************
//
// i|  j|   .---X   Uncoupled expression:
//  |   |  / \          Z_ijkl = sum_ab (nan`b-n`anb) (X_ab * Y_ijbkla)
//  |   | (a  )b
//  |   |  \ /
//  |~~~Y~~~v       Coupled expression:
//  |   |              Z_{ijkl}^{J} = sum_ab (nan`b-n`anb) sum_J' (2J'+1)/(2J+1) ( X_ab * Y_{ijbkla}^{J,J,J'} )
// k|  l|                                          
//                           
//                              
//                  
//
void comm132ss( const Operator& X, const Operator& Y, Operator& Z )
{
  auto& X1 = X.OneBody;
  auto& X3 = X.ThreeBody;
  auto& Y1 = Y.OneBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;
  
  int norb = Z.modelspace->GetNumberOrbits();
  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  for (int ch=0; ch<nch; ch++)
  {
    auto& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;ibra++ ) // <ij| states
    {
      Ket& bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      Orbit& oi = Z.modelspace->GetOrbit(i);
//      Orbit& oj = Z.modelspace->GetOrbit(j);
      for (int iket=0;iket<nkets;iket++ ) // |kl> states
      {
        Ket& ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        double zijkl = 0;
        for (int a=0;a<norb;a++)
        {
          Orbit& oa = Z.modelspace->GetOrbit(a);
          for ( auto b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
          {
            Orbit& ob = Z.modelspace->GetOrbit(b);
            double occfactor = oa.occ - ob.occ;
            if (std::abs(occfactor)<1e-6) continue;
            int twoJ_min = std::abs( oi.j2 - 2*J );
            int twoJ_max = oi.j2 + 2*J;
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ++)
            {
              zijkl += occfactor * (twoJ+1)/(2*J+1) * ( X1(a,b) * Y3.GetME_pn(J,J,twoJ,i,j,b,k,l,a) - Y1(a,b) * X3.GetME_pn(J,J,twoJ,i,j,b,k,l,a) );
            }
          }
        }
        // normalize the tbme
        zijkl /= sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
        Z2.AddToTBME(ch, ch, bra,ket, zijkl );
      }
    }
  }

}






//*****************************************************************************************
//
//  |         |
// i|        j|     Uncoupled expression:
//  |~~~X~~*  |         Z_ijkl = -1/2 sum_abc (nanbn`c-n`an`bnc) ( (1-Pij) X_icab * Y_abjklc - (1-Pkl Yijcabl * Xabkc )
// a|  b|  c\ |
//  |   |    \|
//  |~~~Y~~~~~*      Coupled expression:
//  |   |              Z_{ijkl}^{J} = -1/2 sum_abc (nanbn`c-n`an`bnc) sum_J'J" (2J'+1)(2J"+1)/sqrt(2J+1) (-1)^{2J"+J'-J}
// k|  l|                           *  [   (1 - (-1)^{i+j-J}Pij) (-1)^{j-c} { j  J" J' } X_icab^{J'} * Y_{abjklc}^{J'JJ"}    
//                                                                          { c  i  J  }
//                           
//                                       -(1 - (-1)^{k+l-J}Pkl)  (-1)^{l-c} { l  J" J' } Y_{ijcabl}^{JJ'J"} * X_{abkc}^{J'}  ]
//                                                                          { c  k  J  }
//
void comm232ss( const Operator& X, const Operator& Y, Operator& Z )
{
  auto& X2 = X.TwoBody;
  auto& X3 = X.ThreeBody;
  auto& Y2 = Y.TwoBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z2 = Z.TwoBody;

  int nch = Z.modelspace->GetNumberTwoBodyChannels();
  for (int ch=0; ch<nch; ch++)
  {
    auto& tbc = Z.modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0; ibra<nkets; ibra++)
    {
      Ket& bra = tbc.GetKet(ibra);
      int i=bra.p;
      int j=bra.q;
      Orbit& oi = Z.modelspace->GetOrbit(i);
      Orbit& oj = Z.modelspace->GetOrbit(j);
      for (int iket=0; iket<nkets; iket++)
      {
        double zijkl = 0;
        Ket& ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        Orbit& ok = Z.modelspace->GetOrbit(k);
        Orbit& ol = Z.modelspace->GetOrbit(l);
        for (int ch2=0; ch2<nch; ch2++)
        {
          auto& tbc2 = X.modelspace->GetTwoBodyChannel(ch);
          int J2 = tbc2.J;
          for ( auto iket_ab : tbc2.GetKetIndex_pp() )
          {
            Ket& ket_ab = tbc2.GetKet(iket_ab);
            int a=ket_ab.p;
            int b=ket_ab.q;
            Orbit& oa = Z.modelspace->GetOrbit(a);
            Orbit& ob = Z.modelspace->GetOrbit(b);
            for (auto c : Z.modelspace->holes )
            {
              Orbit& oc = Z.modelspace->GetOrbit(c);
              double occfactor = oa.occ * ob.occ * (1-oc.occ);
              int twoJ_min = std::max( std::abs(oc.j2-2*J), std::abs( oj.j2-2*J2 ) );
              int twoJ_max = std::min( oc.j2+2*J,  oj.j2+2*J2 );
              for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ++)
              {
                double sixj = Z.modelspace->GetSixJ(0.5*oj.j2,0.5*twoJ,J2,0.5*oc.j2,0.5*oi.j2,J);
                zijkl +=  occfactor * Z.modelspace->phase((oj.j2-oc.j2)/2) * sixj * X2.GetTBME(ch,i,c,a,b) * Y3.GetME_pn(J2,J,twoJ,a,b,j,k,l,c);
                zijkl -=  occfactor * Z.modelspace->phase((oj.j2-oc.j2)/2) * sixj * Y2.GetTBME(ch,i,c,a,b) * X3.GetME_pn(J2,J,twoJ,a,b,j,k,l,c);
                sixj = Z.modelspace->GetSixJ(0.5*oi.j2,0.5*twoJ,J2,0.5*oc.j2,0.5*oj.j2,J);
                zijkl -=  occfactor * Z.modelspace->phase(oi.j2+(oj.j2-oc.j2)/2-J) * sixj * X2.GetTBME(ch,j,c,a,b) * Y3.GetME_pn(J2,J,twoJ,a,b,i,k,l,c);
                zijkl +=  occfactor * Z.modelspace->phase(oi.j2+(oj.j2-oc.j2)/2-J) * sixj * Y2.GetTBME(ch,j,c,a,b) * X3.GetME_pn(J2,J,twoJ,a,b,i,k,l,c);

                sixj = Z.modelspace->GetSixJ(0.5*ok.j2,0.5*twoJ,J2,0.5*oc.j2,0.5*ol.j2,J);
                zijkl -=  occfactor * Z.modelspace->phase((ol.j2-oc.j2)/2) * sixj * Y3.GetME_pn(J,J2,twoJ,i,j,c,a,b,l) * X2.GetTBME(ch,a,b,k,c) ;
                zijkl +=  occfactor * Z.modelspace->phase((ol.j2-oc.j2)/2) * sixj * X3.GetME_pn(J,J2,twoJ,i,j,c,a,b,l) * Y2.GetTBME(ch,a,b,k,c) ;
                sixj = Z.modelspace->GetSixJ(0.5*ol.j2,0.5*twoJ,J2,0.5*oc.j2,0.5*ok.j2,J);
                zijkl +=  occfactor * Z.modelspace->phase(ol.j2+(ok.j2-oc.j2)/2-J) * sixj * Y3.GetME_pn(J,J2,twoJ,i,j,c,a,b,k) * X2.GetTBME(ch,a,b,l,c) ;
                zijkl -=  occfactor * Z.modelspace->phase(ol.j2+(ok.j2-oc.j2)/2-J) * sixj * X3.GetME_pn(J,J2,twoJ,i,j,c,a,b,k) * Y2.GetTBME(ch,a,b,l,c) ;
              }
            }
          }
          for ( auto iket_ab : tbc2.GetKetIndex_hh() ) // same damn thing again, but with hhp rather than pph. There must be a better way...
          {
            Ket& ket_ab = tbc2.GetKet(iket_ab);
            int a=ket_ab.p;
            int b=ket_ab.q;
            Orbit& oa = Z.modelspace->GetOrbit(a);
            Orbit& ob = Z.modelspace->GetOrbit(b);
            for (auto c : Z.modelspace->particles )
            {
              Orbit& oc = Z.modelspace->GetOrbit(c);
              double occfactor = (1-oa.occ) * (1-ob.occ) * oc.occ;
              int twoJ_min = std::max( std::abs(oc.j2-2*J), std::abs( oj.j2-2*J2 ) );
              int twoJ_max = std::min( oc.j2+2*J,  oj.j2+2*J2 );
              for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ++)
              {
                double sixj = Z.modelspace->GetSixJ(0.5*oj.j2,0.5*twoJ,J2,0.5*oc.j2,0.5*oi.j2,J);
                zijkl +=  occfactor * Z.modelspace->phase((oj.j2-oc.j2)/2) * sixj * X2.GetTBME(ch,i,c,a,b) * Y3.GetME_pn(J2,J,twoJ,a,b,j,k,l,c);
                zijkl -=  occfactor * Z.modelspace->phase((oj.j2-oc.j2)/2) * sixj * Y2.GetTBME(ch,i,c,a,b) * X3.GetME_pn(J2,J,twoJ,a,b,j,k,l,c);
                sixj = Z.modelspace->GetSixJ(0.5*oi.j2,0.5*twoJ,J2,0.5*oc.j2,0.5*oj.j2,J);
                zijkl -=  occfactor * Z.modelspace->phase(oi.j2+(oj.j2-oc.j2)/2-J) * sixj * X2.GetTBME(ch,j,c,a,b) * Y3.GetME_pn(J2,J,twoJ,a,b,i,k,l,c);
                zijkl +=  occfactor * Z.modelspace->phase(oi.j2+(oj.j2-oc.j2)/2-J) * sixj * Y2.GetTBME(ch,j,c,a,b) * X3.GetME_pn(J2,J,twoJ,a,b,i,k,l,c);

                sixj = Z.modelspace->GetSixJ(0.5*ok.j2,0.5*twoJ,J2,0.5*oc.j2,0.5*ol.j2,J);
                zijkl -=  occfactor * Z.modelspace->phase((ol.j2-oc.j2)/2) * sixj * Y3.GetME_pn(J,J2,twoJ,i,j,c,a,b,l) * X2.GetTBME(ch,a,b,k,c) ;
                zijkl +=  occfactor * Z.modelspace->phase((ol.j2-oc.j2)/2) * sixj * X3.GetME_pn(J,J2,twoJ,i,j,c,a,b,l) * Y2.GetTBME(ch,a,b,k,c) ;
                sixj = Z.modelspace->GetSixJ(0.5*ol.j2,0.5*twoJ,J2,0.5*oc.j2,0.5*ok.j2,J);
                zijkl +=  occfactor * Z.modelspace->phase(ol.j2+(ok.j2-oc.j2)/2-J) * sixj * Y3.GetME_pn(J,J2,twoJ,i,j,c,a,b,k) * X2.GetTBME(ch,a,b,l,c) ;
                zijkl -=  occfactor * Z.modelspace->phase(ol.j2+(ok.j2-oc.j2)/2-J) * sixj * X3.GetME_pn(J,J2,twoJ,i,j,c,a,b,k) * Y2.GetTBME(ch,a,b,l,c) ;
              }
            }
          }
        }
        // normalize the tbme
        zijkl *= -0.5* sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
        Z2.AddToTBME(ch,bra,ket,zijkl);
      }// for iket
    }// for ibra
  }// for ch

}



//*****************************************************************************************
//
// i|  j|  k|       Uncoupled expression:
//  |   |   |---X      Z_ijklmn = sum_a P(ij/k) (X_ka * Y_ijalmn) - P(lm/n) (Y_ijklma * X_an)
//  |   |   |a
//  |~~~Y~~~|       Coupled expression:
//  |   |   |          Z_{ijklmn}^{J1,J2,J} = sum_a P(ij/k)^{J1,J} (X_ka * Y_{ijalmn}^{J1,J2,J})
// l|  m|  n|                                      -P(lm/n)^{J2,J} (Y_{ijklma}^{J1,J2,J} * X_an)
//
//                  This coupled expression is compact, but we can avoid recoupling by applying
//                  the permutations in the uncoupled expression, flipping some indices, and using
//                  the fact that the one-body scalar operator comes with a delta_jj.
//
void comm133ss( const Operator& X, const Operator& Y, Operator& Z )
{
  std::cout << "start comm133" << std::endl;
  int E3cut = 7;
  int e3maxcut = 6;
  auto& X3 = X.ThreeBody;
  auto& Y3 = Y.ThreeBody;
  auto& Z3 = Z.ThreeBody;
  auto& X1 = X.OneBody;
  auto& Y1 = Y.OneBody;
  int norbs = Z.modelspace->GetNumberOrbits();
  if (X3.Norm()<1e-6 and Y3.Norm()<1e-6 ) return;
  for (int i=0;i<norbs; i++)
  {
   Orbit& oi = Z.modelspace->GetOrbit(i);
   for (int j=0; j<=i; j++)
   {
    Orbit& oj = Z.modelspace->GetOrbit(j);
    int Jij_min = std::abs(oi.j2-oj.j2)/2;
    int Jij_max = (oi.j2+oj.j2)/2;
    for (int k=0; k<=j; k++)
    {
     Orbit& ok = Z.modelspace->GetOrbit(k);
     if (  2*( oi.n+oj.n+ok.n)+oi.l+oj.l+ok.l >e3maxcut ) continue;
     for (int l=0; l<=i; l++)
     {
      Orbit& ol = Z.modelspace->GetOrbit(l);
      for (int m=0; m<=l; m++)
      {
       if (l==i and m>j) continue;
       Orbit& om = Z.modelspace->GetOrbit(m);
       int Jlm_min=std::abs(ol.j2-om.j2)/2;
       int Jlm_max=(ol.j2+om.j2)/2;
       for (int n=0; n<=m; n++)
       {
        if (l==i and m==j and n>k) continue;
        Orbit& on = Z.modelspace->GetOrbit(n);
        // check isospin and parity
        if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2>0 ) continue;
        if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;
        if (  2*( ol.n+om.n+on.n)+ol.l+om.l+on.l >e3maxcut ) continue;
//        std::cout << "ijklmn = " << i << " " << j << " " << k << " "<< l << " " << m << " " << n << std::endl;
        for (int Jij=Jij_min; Jij<=Jij_max; Jij++)
        {
         if (i==j and Jij%2>0) continue;
         for (int Jlm=Jlm_min; Jlm<=Jlm_max; Jlm++)
         {
          if (l==m and Jlm%2>0) continue;
          int twoJmin = std::max( std::abs( ok.j2-Jij*2), std::abs( on.j2-Jlm*2) );
          int twoJmax = std::min( ok.j2+Jij*2 , on.j2+Jlm*2 );
          for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
          {
            double zsum =0;
            // First, connect on the bra side
            for (auto a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
              zsum += X1(i,a) * Y3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
              zsum -= Y1(i,a) * X3.GetME_pn(Jij, Jlm, twoJ, a, j, k, l, m, n);
            }
            for (auto a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
            {
              zsum += X1(j,a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
              zsum -= Y1(j,a) * X3.GetME_pn(Jij, Jlm, twoJ, i, a, k, l, m, n);
            }
            for (auto a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
              zsum += X1(k,a) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
              zsum -= Y1(k,a) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, a, l, m, n);
            }
            // Now connect on the ket side
            for (auto a : X.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
            {
              zsum -= X1(a,l) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
              zsum += Y1(a,l) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, a, m, n);
            }
            for (auto a : X.OneBodyChannels.at({om.l,om.j2,om.tz2}) )
            {
              zsum -= X1(a,m) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
              zsum += Y1(a,m) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, a, n);
            }
            for (auto a : X.OneBodyChannels.at({on.l,on.j2,on.tz2}) )
            {
              zsum -= X1(a,n) * Y3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
              zsum += Y1(a,n) * X3.GetME_pn(Jij, Jlm, twoJ, i, j, k, l, m, a);
            }
  
            Z3.AddToME_pn(Jij, Jlm, twoJ, i,j,k,l,m,n, zsum );
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



//*****************************************************************************************
//
//  i|  j|  k|   Uncoupled expression:  
//   |   |~~~|      Z_ijklmn  =  P(ij/k)P(lm/n) sum_a  (X_ijan * Y_aklm - Y_ijan * X_aklm)
//   |   |a  |
//   |~~~|   |   Coupled expression:
//  l|  m|  n|     Z_{ijklmn}^{J1,J2,J}  = -P(ij/k)^{J1,J} P(lm/n)^{J2,J} (-1)^(jn+jk+J1+J2)
//                       *(2J1+1)(2J2+1) * { jn ja J1 } ( X_{ijan}^{J1}-Y_{aklm}^{J2}  - YX)
//                                         { jk  J J2 }
//
void comm223ss( const Operator& X, const Operator& Y, Operator& Z )
{
  int emin = 0;
  int emax = 4;
  int e3maxcut = 6;
  int norbs = Z.modelspace->GetNumberOrbits();
  auto& Z3 = Z.ThreeBody;
  auto& X2 = X.TwoBody;
  auto& Y2 = Y.TwoBody;
  std::cout << "Begin the loop. Norm of X2 and Y2 " << X2.Norm() << " " << Y2.Norm() << std::endl;
  if ( std::abs( X2.Norm() * Y2.Norm() ) < 1e-6) return;
  for (int i=0;i<norbs; i++)
  {
   Orbit& oi = Z.modelspace->GetOrbit(i);
   for (int j=0; j<=i; j++)
   {
    Orbit& oj = Z.modelspace->GetOrbit(j);
    int Jij_min = std::abs(oi.j2-oj.j2)/2;
    int Jij_max = (oi.j2+oj.j2)/2;
    for (int k=0; k<=j; k++)
    {
     Orbit& ok = Z.modelspace->GetOrbit(k);
     if (  2*( oi.n+oj.n+ok.n)+oi.l+oj.l+ok.l >e3maxcut ) continue;
     for (int l=0; l<=i; l++)
     {
      Orbit& ol = Z.modelspace->GetOrbit(l);
      for (int m=0; m<=l; m++)
      {
       if (l==i and m>j) continue;
       Orbit& om = Z.modelspace->GetOrbit(m);
       int Jlm_min=std::abs(ol.j2-om.j2)/2;
       int Jlm_max=(ol.j2+om.j2)/2;
       for (int n=0; n<=m; n++)
       {
        if (l==i and m==j and n>k) continue;
        Orbit& on = Z.modelspace->GetOrbit(n);
        // check isospin and parity
        if ( (oi.l+oj.l+ok.l+ol.l+om.l+on.l)%2>0 ) continue;
        if ( (oi.tz2+oj.tz2+ok.tz2) != (ol.tz2+om.tz2+on.tz2) ) continue;
        if (  2*( ol.n+om.n+on.n)+ol.l+om.l+on.l >e3maxcut ) continue;
//        std::cout << "ijklmn = " << i << " " << j << " " << k << " "<< l << " " << m << " " << n << std::endl;
        for (int Jij=Jij_min; Jij<=Jij_max; Jij++)
        {
         if (i==j and Jij%2>0) continue;
         for (int Jlm=Jlm_min; Jlm<=Jlm_max; Jlm++)
         {
          if (l==m and Jlm%2>0) continue;
          int phasefactor = Z.modelspace->phase( (on.j2+ok.j2)/2 + Jij+Jlm);
          double hatfactor = sqrt( (2*Jij+1)*(2*Jlm+1) );

          int twoJmin = std::max( std::abs( ok.j2-Jij*2), std::abs( on.j2-Jlm*2) );
          int twoJmax = std::min( ok.j2+Jij*2 , on.j2+Jlm*2 );
          for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
          {
           double zdirect = 0;

//            std::cout << "    Jij Jlm twoJ = " << Jij << " " << Jlm << " " << twoJ << std::endl;
            for (int a=0; a<norbs; a++)
            {
             Orbit& oa = Z.modelspace->GetOrbit(a);
             if ( (oi.l+oj.l+oa.l+on.l)%2>0) continue;
             if ( (oa.l+ok.l+ol.l+om.l)%2>0) continue;
             if ( oi.tz2+oj.tz2 != oa.tz2+on.tz2) continue;
             if ( oa.tz2+ok.tz2 != ol.tz2+om.tz2) continue;
             if ( (std::abs(oa.j2-on.j2) > Jij) or (oa.j2+on.j2 < Jij) ) continue;
             if ( (std::abs(oa.j2-ok.j2) > Jlm) or (oa.j2+ok.j2 < Jlm) ) continue;
             double sixj = Z.modelspace->GetSixJ(on.j2*0.5, oa.j2*0.5, Jij, ok.j2*0.5, twoJ*0.5, Jlm );
             if (std::abs(sixj)<1e-6) continue;
             zdirect += -sixj * phasefactor * hatfactor * ( X2.GetTBME_J(Jij, i,j,a,n) * Y2.GetTBME_J(Jlm, a,k,l,m )
                                                         -  Y2.GetTBME_J(Jij, i,j,a,n) * X2.GetTBME_J(Jlm, a,k,l,m ) );
//             std:: cout << "a = " << a << "   sixj phase, hat = " << sixj << " " << phasefactor << " " << hatfactor
//                        << "   < " << oi.j2 << " " << oj.j2 << " " << Jij << " | " << oa.j2 << " " << on.j2 << "  ...  Jlm = " << Jlm 
//                        << "    X2 Y2 = " << std::setprecision(7) << std::setw(12) << X2.GetTBME_J(Jij, i,j,a,n)  << " " <<  Y2.GetTBME_J(Jlm, a,k,l,m )
//                        << "  Y2 X2 =  " <<  Y2.GetTBME_J(Jij, i,j,a,n)  << "  " <<  X2.GetTBME_J(Jlm, a,k,l,m )
//                        << std::endl;
//             if ( std::abs( X2.GetTBME_J(Jlm, a,k,l,m ) )>1e-6 or std::abs( X2.GetTBME_J(Jij, i,j,a,n) ) )
//             {
//               std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    " << zdirect  << std::endl;
//             }
            }

//            std::cout << " ijklmn = " << i << " " << j << " " << k << " " << l << " " << m << " " << n << "   zdirect = " << zdirect << std::endl;
            Z3.AddToME_pn( Jij, Jlm, twoJ, i,j,k,l,m,n, zdirect );
            // now all the permutations. Gulp.

            int Jik_min = std::max( std::abs( oi.j2-ok.j2 ), std::abs( oj.j2 -twoJ))/2;
            int Jik_max = std::min( oi.j2+ok.j2  , oj.j2+twoJ )/2;
            int Jjk_min = std::max( std::abs( oj.j2-ok.j2 ), std::abs( oi.j2 -twoJ))/2;
            int Jjk_max = std::min( oj.j2+ok.j2  , oi.j2+twoJ )/2;
            int Jln_min = std::max( std::abs( ol.j2-on.j2 ), std::abs( om.j2 -twoJ))/2;
            int Jln_max = std::min( ol.j2+on.j2  , om.j2+twoJ )/2;
            int Jmn_min = std::max( std::abs( om.j2-on.j2 ), std::abs( ol.j2 -twoJ))/2;
            int Jmn_max = std::min( om.j2+on.j2  , ol.j2+twoJ )/2;

            std::vector<int> Jik_list;
            std::vector<int> Jjk_list;
            std::vector<int> Jln_list;
            std::vector<int> Jmn_list;
            std::vector<double> recouple_ikj;
            std::vector<double> recouple_kji;
            std::vector<double> recouple_lnm;
            std::vector<double> recouple_nml;
            // loop through the possible recouplings and J values, storing the recoupling coefficients
            // because we need them for the double-recoupling case, like P_ik * P_ln.
//            std::cout << "Single flipsizes permutations..." << std::endl;
            for (int Jik=Jik_min; Jik<=Jik_max; Jik++)
            {
              if (i==k and Jik%2>0) continue;
              Jik_list.push_back(Jik);   // ijk -> ikj  couple ik to Jik
//              std::cout << "    recouple ikj" << std::endl;
              recouple_ikj.push_back( Z3.RecouplingCoefficient( ThreeBodyME::ACB, oi.j2*0.5, oj.j2*0.5, ok.j2*0.5, Jij, Jik, twoJ) );
              if (std::abs(recouple_ikj.back())<1e-6) continue;
              Z3.AddToME_pn( Jik, Jlm, twoJ, i,k,j,l,m,n, -recouple_ikj.back() * zdirect );
            }
            for (int Jjk=Jjk_min; Jjk<=Jjk_max; Jjk++)
            {
              if (j==k and Jjk%2>0) continue;
              Jjk_list.push_back(Jjk);  // ijk -> kji  couple kj to Jkj
//              std::cout << "    recouple kji" << std::endl;
              recouple_kji.push_back( Z3.RecouplingCoefficient( ThreeBodyME::CBA, oi.j2*0.5, oj.j2*0.5, ok.j2*0.5, Jij, Jjk, twoJ) );
              if (std::abs(recouple_kji.back())<1e-6) continue;
              Z3.AddToME_pn( Jjk, Jlm, twoJ, k,j,i,l,m,n, -recouple_kji.back() * zdirect );
            }
            for (int Jln=Jln_min; Jln<=Jln_max; Jln++)
            {
              if (l==n and Jln%2>0) continue;
              Jln_list.push_back(Jln);  // lmn -> lnm  couple ln to Jln
//              std::cout << "    recouple lnm" << std::endl;
              recouple_lnm.push_back( Z3.RecouplingCoefficient( ThreeBodyME::ACB, ol.j2*0.5, om.j2*0.5, on.j2*0.5, Jlm, Jln, twoJ) );
              Z3.AddToME_pn( Jij, Jln, twoJ, i,j,k,l,n,m, -recouple_lnm.back() * zdirect );
            }
            for (int Jmn=Jmn_min; Jmn<=Jmn_max; Jmn++)
            {
              if (m==n and Jmn%2>0) continue;
              Jmn_list.push_back(Jmn);  // lmn -> nml  couple mn to Jmn
//              std::cout << "    recouple nml" << std::endl;
              recouple_nml.push_back( Z3.RecouplingCoefficient( ThreeBodyME::CBA, ol.j2*0.5, om.j2*0.5, on.j2*0.5, Jlm, Jmn, twoJ) );
              Z3.AddToME_pn( Jij, Jmn, twoJ, i,j,k,n,m,l, -recouple_nml.back() * zdirect );
            }

//            std::cout << "Double flipsizes permutations..." << std::endl;
            // don't forget about the double-swapsies...
            for (size_t indx1=0;indx1<Jik_list.size();indx1++)
            {
              for (size_t indx2=0;indx2<Jln_list.size();indx2++)
              {
//                std::cout << "    recouple ikj,lnm : " << Jik_list[indx1] << " " << Jln_list[indx2] << std::endl;
                Z3.AddToME_pn( Jik_list[indx1], Jln_list[indx2], twoJ, i,k,j,l,n,m, recouple_ikj[indx1]*recouple_lnm[indx2] * zdirect );
              }
              for (size_t indx2=0;indx2<Jmn_list.size();indx2++)
              {
//                std::cout << "    recouple ikj,nml" << std::endl;
                Z3.AddToME_pn( Jik_list[indx1], Jmn_list[indx2], twoJ, i,k,j,n,m,l, recouple_ikj[indx1]*recouple_nml[indx2] * zdirect );
              }
            }
            for (size_t indx1=0;indx1<Jjk_list.size();indx1++)
            {
              for (size_t indx2=0;indx2<Jln_list.size();indx2++)
              {
//                std::cout << "    recouple kji,lnm" << std::endl;
                Z3.AddToME_pn( Jjk_list[indx1], Jln_list[indx2], twoJ, k,j,i,l,n,m, recouple_kji[indx1]*recouple_lnm[indx2] * zdirect );
              }
              for (size_t indx2=0;indx2<Jmn_list.size();indx2++)
              {
//                std::cout << "    recouple kji,nml" << std::endl;
                Z3.AddToME_pn( Jjk_list[indx1], Jmn_list[indx2], twoJ, k,j,i,n,m,l, recouple_kji[indx1]*recouple_nml[indx2] * zdirect );
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
  }
  std::cout << "Done with comm223ss. Norm of Z3 is " << Z.ThreeBodyNorm() << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
////////////   BEGIN SCALAR-TENSOR COMMUTATORS      //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |
//
// This is no different from the scalar-scalar version
void comm111st( const Operator & X, const Operator& Y, Operator& Z)
{
   double tstart = omp_get_wtime();
   comm111ss(X,Y,Z);
   X.profiler.timer["comm111st"] += omp_get_wtime() - tstart;
}


//*****************************************************************************************
//                                       |
//      i |              i |             |
//        |    ___.Y       |__X__        |
//        |___(_)    _     |   (_)__.    |  [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab 
//      j | X            j |        Y    |        * sum_J (2J+1) x_biaj^(J)  
//                                       |      
//---------------------------------------*        = 1/(2j+1) sum_a n_a sum_J (2J+1)
//                                                  * sum_b y_ab x_biaj - yba x_aibj
//
// X is scalar one-body, Y is tensor two-body
// There must be a better way to do this looping. 
//
//void Operator::comm121st( Operator& Y, Operator& Z) 
void comm121st( const Operator& X, const Operator& Y, Operator& Z) 
{

   double tstart = omp_get_wtime();
   int norbits = Z.modelspace->GetNumberOrbits();
   int Lambda = Z.GetJRank();
   Z.modelspace->PreCalculateSixJ();
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.rank_J))
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double ji = 0.5*oi.j2;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
      {
          Orbit &oj = Z.modelspace->GetOrbit(j);
          double jj = 0.5*oj.j2;
          if (j<i) continue; // only calculate upper triangle
          double& Zij = Z.OneBody(i,j);
          for (auto& a : Z.modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = Z.modelspace->GetOrbit(a);
             double ja = 0.5*oa.j2;
             for (auto& b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
             {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ);
                  int J1min = std::abs(ji-ja);
                  int J1max = ji + ja;
                  for (int J1=J1min; J1<=J1max; ++J1)
                  {
                    int phasefactor = Z.modelspace->phase(jj+ja+J1+Lambda);
                    int J2min = std::max(std::abs(Lambda - J1),std::abs(int(ja-jj)));
                    int J2max = std::min(Lambda + J1,int(ja+jj));
                    for (int J2=J2min; J2<=J2max; ++J2)
                    {
                      if ( ! ( J2>=std::abs(ja-jj) and J2<=ja+jj )) continue;
                      double prefactor = nanb*phasefactor * sqrt((2*J1+1)*(2*J2+1)) * Z.modelspace->GetSixJ(J1,J2,Lambda,jj,ji,ja);
                      Zij +=  prefactor * ( X.OneBody(a,b) * Y.TwoBody.GetTBME_J(J1,J2,b,i,a,j) - X.OneBody(b,a) * Y.TwoBody.GetTBME_J(J1,J2,a,i,b,j ));
                    }
                  }
             }
             // Now, X is scalar two-body and Y is tensor one-body
             for (auto& b : Y.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
             {

                Orbit &ob = Z.modelspace->GetOrbit(b);
                double jb = 0.5*ob.j2;
                if (std::abs(ob.occ-1) < OCC_CUT) continue;
                double nanb = oa.occ * (1-ob.occ);
                int J1min = std::max(std::abs(ji-jb),std::abs(jj-ja));
                int J1max = std::min(ji+jb,jj+ja);
                double zij = 0;
                for (int J1=J1min; J1<=J1max; ++J1)
                {
                  zij -= Z.modelspace->phase(ji+jb+J1) * (2*J1+1) * Z.modelspace->GetSixJ(ja,jb,Lambda,ji,jj,J1) * X.TwoBody.GetTBME_J(J1,J1,b,i,a,j);
                }

                J1min = std::max(std::abs(ji-ja),std::abs(jj-jb));
                J1max = std::min(ji+ja,jj+jb);
                for (int J1=J1min; J1<=J1max; ++J1)
                {
                  zij += Z.modelspace->phase(ji+jb+J1) * (2*J1+1) * Z.modelspace->GetSixJ(jb,ja,Lambda,ji,jj,J1) * X.TwoBody.GetTBME_J(J1,J1,a,i,b,j) ;
                }

                Zij += nanb * Y.OneBody(a,b) * zij;
             }

             
          }
      }
   }
   
   X.profiler.timer["comm121st"] += omp_get_wtime() - tstart;
}




//*****************************************************************************************
//
//    |     |               |      |           [X2,Y1](2) = sum_a ( Y_ia X_ajkl + Y_ja X_iakl - Y_ak X_ijal - Y_al X_ijka )
//    |     |___.Y          |__X___|         
//    |     |         _     |      |          
//    |_____|               |      |_____.Y        
//    |  X  |               |      |            
//
// -- AGREES WITH NATHAN'S RESULTS
// Right now, this is the slowest one...
// Agrees with previous code in the scalar-scalar limit
//void Operator::comm122st( Operator& Y, Operator& Z ) 
//void Operator::comm122st( const Operator& X, const Operator& Y ) 
void comm122st( const Operator& X, const Operator& Y , Operator& Z) 
{
   double tstart = omp_get_wtime();
   int Lambda = Z.rank_J;

    std::vector< int > bra_channels;
    std::vector< int > ket_channels;
    for ( auto& itmat : Z.TwoBody.MatEl )
    {
      bra_channels.push_back( itmat.first[0] );
      ket_channels.push_back( itmat.first[1] );
    }
    int nmat = bra_channels.size();
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.rank_J))
    for (int ii=0; ii<nmat; ++ii)
    {
     int ch_bra = bra_channels[ii];
     int ch_ket = ket_channels[ii];

      TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      double hatfactor = sqrt((2*J1+1)*(2*J2+1));
      arma::mat& Z2 = Z.TwoBody.GetMatrix(ch_bra,ch_ket);

      for (int ibra = 0;ibra<nbras; ++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = Z.modelspace->GetOrbit(i);
         Orbit& oj = Z.modelspace->GetOrbit(j);
         double ji = oi.j2/2.0;
         double jj = oj.j2/2.0;
         for (int iket=0;iket<nkets; ++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit& ok = Z.modelspace->GetOrbit(k);
            Orbit& ol = Z.modelspace->GetOrbit(l);
            double jk = ok.j2/2.0;
            double jl = ol.j2/2.0;

            double cijkl = 0;
            double c1 = 0;
            double c2 = 0;
            double c3 = 0;
            double c4 = 0;

            for ( int a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
              c1 += X.OneBody(i,a) * Y.TwoBody.GetTBME(ch_bra,ch_ket,a,j,k,l);
            }
            if (i==j)
            {
              c2 = c1; // there should be a phase here, but if the ket exists, it'd better be +1.
            }
            else
            {
              for ( int a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
              {
                 c2 += X.OneBody(j,a) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,a,k,l);
              }
            }
            for ( int a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
               c3 += X.OneBody(a,k) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,j,a,l);
            }
            if (k==l)
            {
              c4 = c3 ;
            }
            else
            {
              for ( int a : X.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
              {
                 c4 += X.OneBody(a,l) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,j,k,a);
              }
            }

            cijkl = c1 + c2 - c3 - c4;


            c1=0;
            c2=0;
            c3=0;
            c4=0;
            int phase1 = Z.modelspace->phase(ji+jj+J2+Lambda);
            int phase2 = Z.modelspace->phase(J1-J2+Lambda);
            int phase3 = Z.modelspace->phase(J1-J2+Lambda);
            int phase4 = Z.modelspace->phase(jk+jl-J1+Lambda);


            for ( int a : Y.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
               double ja = Z.modelspace->GetOrbit(a).j2*0.5;
               c1 -=   Z.modelspace->GetSixJ(J2,J1,Lambda,ji,ja,jj) * Y.OneBody(i,a) * X.TwoBody.GetTBME(ch_ket,ch_ket,a,j,k,l) ;
            }
            if (i==j)
            {
              c2 = -c1;
            }
            else
            {
              for ( int a : Y.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
              {
                 double ja = Z.modelspace->GetOrbit(a).j2*0.5;
                 c2 +=  Z.modelspace->GetSixJ(J2,J1,Lambda,jj,ja,ji) * Y.OneBody(j,a) * X.TwoBody.GetTBME(ch_ket,ch_ket,a,i,k,l);
              }
            }
            for ( int a : Y.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
               double ja = Z.modelspace->GetOrbit(a).j2*0.5;
               c3 -=  Z.modelspace->GetSixJ(J1,J2,Lambda,jk,ja,jl) * Y.OneBody(a,k) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,j,l,a) ;
            }
            if (k==l)
            {
              c4 = -c3;
            }
            else
            {
              for ( int a : Y.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
              {
                 double ja = Z.modelspace->GetOrbit(a).j2*0.5;
                 c4 +=  Z.modelspace->GetSixJ(J1,J2,Lambda,jl,ja,jk) * Y.OneBody(a,l) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,j,k,a) ;
              }
            }
            cijkl += hatfactor*(phase1*c1+phase2*c2+phase3*c3+phase4*c4);


            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            Z2(ibra,iket) += cijkl /norm;
         }
      }
   }
   X.profiler.timer["comm122st"] += omp_get_wtime() - tstart;
}





// Since comm222_pp_hh and comm211 both require the construction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
// X is a scalar, Y is a tensor
//void Operator::comm222_pp_hh_221st( Operator& Y, Operator& Z )  
//void Operator::comm222_pp_hh_221st( const Operator& X, const Operator& Y )  
void comm222_pp_hh_221st( const Operator& X, const Operator& Y, Operator& Z )  
{

   double tstart = omp_get_wtime();
   int Lambda = Z.GetJRank();

   TwoBodyME Mpp = Y.TwoBody;
   TwoBodyME Mhh = Y.TwoBody;
   TwoBodyME Mff = Y.TwoBody;

   std::vector<int> vch_bra;
   std::vector<int> vch_ket;
   std::vector<const arma::mat*> vmtx;
   for ( auto& itmat : Y.TwoBody.MatEl )
   {
     vch_bra.push_back(itmat.first[0]);
     vch_ket.push_back(itmat.first[1]);
     vmtx.push_back(&(itmat.second));
   }
   size_t nchan = vch_bra.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (size_t i=0;i<nchan; ++i)
   {
    int ch_bra = vch_bra[i];
    int ch_ket = vch_ket[i];

    TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);

    auto& LHS1 = X.TwoBody.GetMatrix(ch_bra,ch_bra);
    auto& LHS2 = X.TwoBody.GetMatrix(ch_ket,ch_ket);

    auto& RHS  =  *vmtx[i];

    arma::mat& Matrixpp =  Mpp.GetMatrix(ch_bra,ch_ket);
    arma::mat& Matrixhh =  Mhh.GetMatrix(ch_bra,ch_ket);
   
    const arma::uvec& bras_pp = tbc_bra.GetKetIndex_pp();
    const arma::uvec& bras_hh = tbc_bra.GetKetIndex_hh();
    const arma::uvec& bras_ph = tbc_bra.GetKetIndex_ph();
    const arma::uvec& kets_pp = tbc_ket.GetKetIndex_pp();
    const arma::uvec& kets_hh = tbc_ket.GetKetIndex_hh();
    const arma::uvec& kets_ph = tbc_ket.GetKetIndex_ph();

    // the complicated-looking construct after the % signs just multiply the matrix elements by the proper occupation numbers (nanb, etc.)

    arma::mat MLeft  = join_horiz( LHS1.cols(bras_hh) , -RHS.cols(kets_hh) );
    arma::mat MRight = join_vert( RHS.rows(bras_hh)  % tbc_bra.Ket_occ_hh.cols( arma::uvec(RHS.n_cols,arma::fill::zeros ) ),
                                 LHS2.rows(kets_hh)  % tbc_ket.Ket_occ_hh.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) ));

    Matrixhh = MLeft * MRight;


    MLeft  = join_horiz( LHS1.cols(join_vert(bras_pp,join_vert(bras_hh,bras_ph))), -RHS.cols(join_vert(kets_pp,join_vert(kets_hh,kets_ph))) );
    MRight = join_vert( join_vert(     RHS.rows(bras_pp), 
                          join_vert( RHS.rows(bras_hh)  % tbc_bra.Ket_unocc_hh.cols( arma::uvec(RHS.n_cols,arma::fill::zeros) )  ,
                                     RHS.rows(bras_ph)  % tbc_bra.Ket_unocc_ph.cols( arma::uvec(RHS.n_cols,arma::fill::zeros) ) )),
                      join_vert(     LHS2.rows(kets_pp),
                          join_vert( LHS2.rows(kets_hh) % tbc_ket.Ket_unocc_hh.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) ),
                                     LHS2.rows(kets_ph) % tbc_ket.Ket_unocc_ph.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) )))
                     );

    Matrixpp = MLeft * MRight;
                                

    if (Z.GetParticleRank()>1)
    {
      Z.TwoBody.GetMatrix(ch_bra,ch_ket) += Matrixpp - Matrixhh;
    }

   }// for itmat

      // The one body part takes some additional work

   int norbits = Z.modelspace->GetNumberOrbits();
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.rank_J))
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double ji = oi.j2/2.0;
      for (int j : Z.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         if (j<i) continue;
         Orbit &oj = Z.modelspace->GetOrbit(j);
         double jj = oj.j2/2.0;
         double cijJ = 0;
         // Sum c over holes and include the nbar_a * nbar_b terms
           for (auto& c : Z.modelspace->holes)
           {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2/2.0;
              int j1min = std::abs(jc-ji);
              int j1max = jc+ji;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
               int j2min = std::max( int(std::abs(jc-jj)), std::abs(Lambda-J1) );
               int j2max = std::min( int(jc+jj), J1+Lambda ); 
               for (int J2=j2min; J2<=j2max; ++J2)
               {
                double hatfactor = sqrt( (2*J1+1)*(2*J2+1) );
                double sixj = Z.modelspace->GetSixJ(J1, J2, Lambda, jj, ji, jc);
                cijJ += hatfactor * sixj * Z.modelspace->phase(jj + jc + J1 + Lambda)
                        * ( oc.occ * Mpp.GetTBME_J(J1,J2,c,i,c,j) + (1-oc.occ) * Mhh.GetTBME_J(J1,J2,c,i,c,j));
               }
              }
           // Sum c over particles and include the n_a * n_b terms
           }
           for (auto& c : Z.modelspace->particles)
           {
              Orbit &oc = Z.modelspace->GetOrbit(c);
              double jc = oc.j2/2.0;
              int j1min = std::abs(jc-ji);
              int j1max = jc+ji;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
               int j2min = std::max( int(std::abs(jc-jj)), std::abs(Lambda-J1) );
               int j2max = std::min( int(jc+jj), J1+Lambda ); 
               for (int J2=j2min; J2<=j2max; ++J2)
               {
                double hatfactor = sqrt( (2*J1+1)*(2*J2+1) );
                double sixj = Z.modelspace->GetSixJ(J1, J2, Lambda, jj, ji, jc);
                cijJ += hatfactor * sixj * Z.modelspace->phase(jj + jc + J1 + Lambda) * Mhh.GetTBME_J(J1,J2,c,i,c,j);
               }
              }
           }
//         #pragma omp critical
         Z.OneBody(i,j) += cijJ ;
      } // for j
    } // for i
    X.profiler.timer["comm222_pp_hh_221st"] += omp_get_wtime() - tstart;
}






//**************************************************************************
//
//  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
//                        { k l J'}
// TENSOR VARIETY
/// The scalar Pandya transformation is defined as
/// \f[
///  \bar{X}^{J}_{i\bar{j}k\bar{l}} = - \sum_{J'} (2J'+1)
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  X^{J}_{ilkj}
/// \f]
/// where the overbar indicates time-reversed orbits.
/// This function is designed for use with comm222_phss() and so it takes in
/// two arrays of matrices, one for hp terms and one for ph terms.
//void Operator::DoTensorPandyaTransformation(vector<arma::mat>& TwoBody_CC_hp, vector<arma::mat>& TwoBody_CC_ph)
//void Operator::DoTensorPandyaTransformation(map<array<int,2>,arma::mat>& TwoBody_CC_hp, map<array<int,2>,arma::mat>& TwoBody_CC_ph) const
void DoTensorPandyaTransformation( const Operator& Z, std::map<std::array<index_t,2>,arma::mat>& TwoBody_CC_ph)
{
   int Lambda = Z.rank_J;
   // loop over cross-coupled channels
   index_t nch = Z.modelspace->SortedTwoBodyChannels_CC.size();

   // Allocate map for matrices -- this needs to be serial.
   for ( index_t ch_bra_cc : Z.modelspace->SortedTwoBodyChannels_CC )
   {
      TwoBodyChannel& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
      index_t nph_bras = bras_ph.n_rows;
      for ( index_t ch_ket_cc : Z.modelspace->SortedTwoBodyChannels_CC )
      {
        TwoBodyChannel& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
        index_t nKets_cc = tbc_ket_cc.GetNumberKets();

         TwoBody_CC_ph[{ch_bra_cc,ch_ket_cc}] =  arma::mat(2*nph_bras,   nKets_cc, arma::fill::zeros);
      }
   }

   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Lambda))
   for (index_t ich=0;ich<nch;++ich)
   {
      index_t ch_bra_cc = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      TwoBodyChannel& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      int Jbra_cc = tbc_bra_cc.J;
      arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
//      arma::uvec& bras_ph = tbc_bra_cc.GetKetIndex_ph();
      index_t nph_bras = bras_ph.size();

      for ( index_t ch_ket_cc : Z.modelspace->SortedTwoBodyChannels_CC )
      {
        TwoBodyChannel& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
        int Jket_cc = tbc_ket_cc.J;
        if ( (Jbra_cc+Jket_cc < Z.GetJRank()) or std::abs(Jbra_cc-Jket_cc)>Z.GetJRank() ) continue;
        if ( (tbc_bra_cc.parity + tbc_ket_cc.parity + Z.GetParity())%2>0 ) continue;

        index_t nKets_cc = tbc_ket_cc.GetNumberKets();

//        arma::mat& MatCC_hp = TwoBody_CC_hp[{ch_bra_cc,ch_ket_cc}];
        arma::mat& MatCC_ph = TwoBody_CC_ph[{ch_bra_cc,ch_ket_cc}];
        // loop over ph bras <ad| in this channel
        for (index_t ibra=0; ibra<nph_bras; ++ibra)
        {
           Ket & bra_cc = tbc_bra_cc.GetKet( bras_ph[ibra] );
           index_t a = bra_cc.p;
           index_t b = bra_cc.q;
           Orbit & oa = Z.modelspace->GetOrbit(a);
           Orbit & ob = Z.modelspace->GetOrbit(b);
           double ja = oa.j2*0.5;
           double jb = ob.j2*0.5;

           // loop over kets |bc> in this channel
           index_t iket_max =  nKets_cc ;
           for (index_t iket_cc=0; iket_cc<iket_max; ++iket_cc)
           {
              Ket & ket_cc = tbc_ket_cc.GetKet(iket_cc%nKets_cc);
              index_t c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
              index_t d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
              Orbit & oc = Z.modelspace->GetOrbit(c);
              Orbit & od = Z.modelspace->GetOrbit(d);
              double jc = oc.j2*0.5;
              double jd = od.j2*0.5;


              int j1min = std::abs(ja-jd);
              int j1max = ja+jd;
              double sm = 0;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
                int j2min = std::max( int(std::abs(jc-jb)), std::abs(J1-Lambda) );
                int j2max = std::min(int(jc+jb),J1+Lambda);
                for (int J2=j2min; J2<=j2max; ++J2)
                {
                  double ninej = Z.modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
                  if (std::abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
                  double tbme = Z.TwoBody.GetTBME_J(J1,J2,a,d,c,b);
                  sm -= hatfactor * Z.modelspace->phase(jb+jd+Jket_cc+J2) * ninej * tbme ;
                }
              }
              MatCC_ph(ibra,iket_cc) = sm;

              // Exchange (a <-> b) to account for the (n_a - n_b) term
              // Get Tz,parity and range of J for <bd || ca > coupling
              j1min = std::abs(jb-jd);
              j1max = jb+jd;
              sm = 0;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
                int j2min = std::max( int(std::abs(jc-ja)), std::abs(J1-Lambda) );
                int j2max = std::min(int(jc+ja),J1+Lambda);
                for (int J2=j2min; J2<=j2max; ++J2)
                {
                  double ninej = Z.modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
                  if (std::abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
                  double tbme = Z.TwoBody.GetTBME_J(J1,J2,b,d,c,a);
                  sm -= hatfactor * Z.modelspace->phase(ja+jd+Jket_cc+J2) * ninej * tbme ;
                }
              }
              MatCC_ph(ibra+nph_bras,iket_cc) = sm;

           }
        }
    }
   }
}



// This happens inside an OMP loop, and so everything here needs to be thread safe
//
//void Operator::DoTensorPandyaTransformation_SingleChannel( arma::mat& TwoBody_CC_ph, int ch_bra_cc, int ch_ket_cc) const
//void Operator::DoTensorPandyaTransformation_SingleChannel( arma::mat& MatCC_ph, int ch_bra_cc, int ch_ket_cc) const
void DoTensorPandyaTransformation_SingleChannel( const Operator& Z, arma::mat& MatCC_ph, int ch_bra_cc, int ch_ket_cc)
{
   int Lambda = Z.rank_J;

   TwoBodyChannel& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
   arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
   int nph_bras = bras_ph.n_rows;

   TwoBodyChannel& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
   int nKets_cc = tbc_ket_cc.GetNumberKets();

   // The Pandya-transformed (formerly cross-coupled) particle-hole type matrix elements
   // (this is the output of this method)
   MatCC_ph =  arma::mat(2*nph_bras,   nKets_cc, arma::fill::zeros);

   int Jbra_cc = tbc_bra_cc.J;
   int Jket_cc = tbc_ket_cc.J;
   if ( (Jbra_cc+Jket_cc < Z.GetJRank() ) or std::abs(Jbra_cc-Jket_cc)>Z.GetJRank() ) return;
   if ( (tbc_bra_cc.parity + tbc_ket_cc.parity + Z.GetParity())%2>0 ) return;


   // loop over ph bras <ad| in this channel
   for (int ibra=0; ibra<nph_bras; ++ibra)
   {
      Ket & bra_cc = tbc_bra_cc.GetKet( bras_ph[ibra] );
      int a = bra_cc.p;
      int b = bra_cc.q;
      Orbit & oa = Z.modelspace->GetOrbit(a);
      Orbit & ob = Z.modelspace->GetOrbit(b);
      double ja = oa.j2*0.5;
      double jb = ob.j2*0.5;

      // loop over kets |bc> in this channel
      int iket_max =  nKets_cc ;
      for (int iket_cc=0; iket_cc<iket_max; ++iket_cc)
      {
         Ket & ket_cc = tbc_ket_cc.GetKet(iket_cc%nKets_cc);
         int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
         int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
         Orbit & oc = Z.modelspace->GetOrbit(c);
         Orbit & od = Z.modelspace->GetOrbit(d);
         double jc = oc.j2*0.5;
         double jd = od.j2*0.5;


         int j1min = std::abs(ja-jd);
         int j1max = ja+jd;
         double sm = 0;
         for (int J1=j1min; J1<=j1max; ++J1)
         {
           int j2min = std::max( int(std::abs(jc-jb)), std::abs(J1-Lambda) );
           int j2max = std::min( int(jc+jb), J1+Lambda) ;
           for (int J2=j2min; J2<=j2max; ++J2)
           {
             double ninej = Z.modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
             if (std::abs(ninej) < 1e-8) continue;
             double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
             double tbme = Z.TwoBody.GetTBME_J(J1,J2,a,d,c,b);
             sm -= hatfactor * Z.modelspace->phase(jb+jd+Jket_cc+J2) * ninej * tbme ;
           }
         }
         MatCC_ph(ibra,iket_cc) = sm;

         // Exchange (a <-> b) to account for the (n_a - n_b) term
         if (a==b)
         {
           MatCC_ph(ibra+nph_bras,iket_cc) = sm;
         }
         else
         {

           // Get Tz,parity and range of J for <bd || ca > coupling
           j1min = std::abs(jb-jd);
           j1max = jb+jd;
           sm = 0;
           for (int J1=j1min; J1<=j1max; ++J1)
           {
             int j2min = std::max( int(std::abs(jc-ja)), std::abs(J1-Lambda) );
             int j2max = std::min( int(jc+ja), J1+Lambda );
             for (int J2=j2min; J2<=j2max; ++J2)
             {
               double ninej = Z.modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
               if (std::abs(ninej) < 1e-8) continue;
               double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
               double tbme = Z.TwoBody.GetTBME_J(J1,J2,b,d,c,a);
               sm -= hatfactor * Z.modelspace->phase(ja+jd+Jket_cc+J2) * ninej * tbme ;
             }
           }
           MatCC_ph(ibra+nph_bras,iket_cc) = sm;
         }
      }
   }
}





// Take Pandya-transformed matrix Zbar for a single channel, invert the Pandya transformation and add the result to the current operator.
//void Operator::AddInverseTensorPandyaTransformation_SingleChannel(arma::mat& Zbar, int ch_bra_cc, int ch_ket_cc)
//void AddInverseTensorPandyaTransformation_SingleChannel(Operator& Z, arma::mat& Zbar, int ch_bra_cc, int ch_ket_cc)
void AddInverseTensorPandyaTransformation_SingleChannel(Operator& Z, arma::mat& Zbar, size_t ch_bra_cc, size_t ch_ket_cc)
{
    // Do the inverse Pandya transform
   if (ch_bra_cc > ch_ket_cc)  // hopefully this won't happen
   {
      std::cout << "WARNING: Called Operator::AddInverseTensorPandyaTransformation_SingleChannel with ch_bra_cc > ch_ket_cc : " << ch_bra_cc << " > " << ch_ket_cc << std::endl;
      std::cout << " Skipping this channel." << std::endl;
      return;
   }
//   int n_channels_kept = 0;
   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup(Z.GetJRank(), Z.GetTRank(),Z.GetParity())[{ch_bra_cc,ch_ket_cc}];
   int Lambda = Z.rank_J;
   int hZ = Z.IsHermitian() ? 1 : -1;
   TwoBodyChannel_CC& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
   TwoBodyChannel_CC& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
   int J3 = tbc_bra_cc.J;
   int J4 = tbc_ket_cc.J;
   int Tz_bra_cc = tbc_bra_cc.Tz;
   int Tz_ket_cc = tbc_ket_cc.Tz;
   int parity_bra_cc = tbc_bra_cc.parity;
   int parity_ket_cc = tbc_ket_cc.parity;
   int nkets_cc = tbc_ket_cc.GetNumberKets();

   int nchannels = pandya_lookup[0].size();
   for (int ich=0;ich<nchannels;++ich)
   {
      int ch_bra = pandya_lookup[0][ich];
      int ch_ket = pandya_lookup[1][ich];
      TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      int nBras = tbc_bra.GetNumberKets();
      int nKets = tbc_ket.GetNumberKets();
      arma::mat& Zijkl = Z.TwoBody.GetMatrix(ch_bra,ch_ket);
//      bool inner_loop = false;

      double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );

      for (int ibra=0; ibra<nBras; ++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = Z.modelspace->GetOrbit(i);
         Orbit & oj = Z.modelspace->GetOrbit(j);
         double ji = 0.5*oi.j2;
         double jj = 0.5*oj.j2;
         int ketmin = ch_bra==ch_ket ? ibra : 0;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = Z.modelspace->GetOrbit(k);
            Orbit & ol = Z.modelspace->GetOrbit(l);
            double jk = 0.5*ok.j2;
            double jl = 0.5*ol.j2;

            double commij = 0;
            double commji = 0;

            // Transform Z_ilkj
            int j3min = std::abs(int(ji-jl));
            int j3max = ji+jl;
            int j4min = std::abs(int(jk-jj));
            int j4max = int(jk+jj);

            if (   (oi.l+ol.l)%2==parity_bra_cc             and (ok.l+oj.l)%2==parity_ket_cc
                      and std::abs(oi.tz2+ol.tz2)==2*Tz_bra_cc   and std::abs(ok.tz2+oj.tz2)==2*Tz_ket_cc
                      and j3min<=J3 and J3<=j3max           and j4min<=J4 and J4<=j4max )
            {
               int indx_il = tbc_bra_cc.GetLocalIndex( std::min(i,l), std::max(i,l) );
               int indx_kj = tbc_ket_cc.GetLocalIndex( std::min(j,k), std::max(j,k) );

               double ninej = Z.modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
               double tbme = 0;

               if (i<=l) tbme = Zbar( indx_il ,indx_kj+(k>j?nkets_cc:0) );
               else      tbme = Zbar( indx_il ,indx_kj+(k>j?0:nkets_cc) ) * hZ * Z.modelspace->phase(J3-J4 + ji+jj+jk+jl);
               commij += hatfactor * Z.modelspace->phase(jj+jl+J2+J4) * ninej * tbme ;
//               inner_loop = true;
            }

            if (  (ch_bra_cc != ch_ket_cc)  and ((oi.l+ol.l)%2==parity_ket_cc)           and ((ok.l+oj.l)%2==parity_bra_cc)
                                            and (std::abs(oi.tz2+ol.tz2)==2*Tz_ket_cc)   and (std::abs(ok.tz2+oj.tz2)==2*Tz_bra_cc)
                                            and (j3min<=J4) and (J4<=j3max)  and  (j4min<=J3) and (J3<=j4max )    )
              {
                 int indx_kj = tbc_bra_cc.GetLocalIndex( std::min(j,k), std::max(j,k));
                 int indx_il = tbc_ket_cc.GetLocalIndex( std::min(i,l), std::max(i,l));

                 double ninej = Z.modelspace->GetNineJ(ji,jl,J4,jj,jk,J3,J1,J2,Lambda);
                 double tbme = 0;

                 if(k<=j) tbme = Zbar(indx_kj, indx_il+(i>l?nkets_cc:0) ) * hZ * Z.modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                 else     tbme = Zbar(indx_kj, indx_il+(i>l?0:nkets_cc) ) * Z.modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)

                 commij += hatfactor * Z.modelspace->phase(jj+jl+J2+J3) * ninej * tbme ;
//               inner_loop = true;
              }


            if (i==j)
            {
              commji = commij;
            }
            else
            {
              j3min = std::abs(int(jj-jl));
              j3max = jj+jl;
              j4min = std::abs(int(jk-ji));
              j4max = int(jk+ji);
              if (   ((oj.l+ol.l)%2==parity_bra_cc) and ((ok.l+oi.l)%2==parity_ket_cc)
                      and (std::abs(oj.tz2+ol.tz2)==2*Tz_bra_cc) and (std::abs(ok.tz2+oi.tz2)==2*Tz_ket_cc)
                      and (J3>=j3min) and (J3<=j3max) and (J4>=j4min) and (J4<=j4max)   )
              {
                // Transform Z_jlki
                int indx_jl = tbc_bra_cc.GetLocalIndex( std::min(j,l), std::max(j,l) );
                int indx_ki = tbc_ket_cc.GetLocalIndex( std::min(k,i), std::max(k,i) );
                double ninej = Z.modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);
                double tbme = 0;
                if (j<=l) tbme = Zbar( indx_jl ,indx_ki+(k>i?nkets_cc:0) );
                else      tbme = Zbar( indx_jl ,indx_ki+(k>i?0:nkets_cc) ) * hZ * Z.modelspace->phase(J3-J4 + ji+jj+jk+jl);
                commji += hatfactor * Z.modelspace->phase(ji+jl+J2+J4) * ninej * tbme ;
//               inner_loop = true;
              }
              if ( (ch_bra_cc!=ch_ket_cc) and ((oj.l+ol.l)%2==parity_ket_cc) and ((ok.l+oi.l)%2==parity_bra_cc)
                                       and (std::abs(oj.tz2+ol.tz2)==2*Tz_ket_cc) and (std::abs(ok.tz2+oi.tz2)==2*Tz_bra_cc)
                                       and (J4>=j3min) and (J4<=j3max) and (J3>=j4min) and (J3<=j4max)   )
              {
                int indx_jl = tbc_ket_cc.GetLocalIndex( std::min(j,l), std::max(j,l));
                int indx_ki = tbc_bra_cc.GetLocalIndex( std::min(k,i), std::max(k,i));
                double ninej = Z.modelspace->GetNineJ(jj,jl,J4,ji,jk,J3,J1,J2,Lambda);
                double tbme = 0;

                if(k<=i) tbme = Zbar(indx_ki, indx_jl+(j>l?nkets_cc:0) ) * hZ * Z.modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                else     tbme = Zbar(indx_ki, indx_jl+(j>l?0:nkets_cc) ) * Z.modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
                commji += hatfactor * Z.modelspace->phase(ji+jl+J2+J3) * ninej * tbme ;
//               inner_loop = true;
              }
            }


            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            #pragma omp atomic
            Zijkl(ibra,iket) +=  (commij - Z.modelspace->phase(ji+jj-J1)*commji) / norm;
            if (ch_bra==ch_ket) 
            {
               #pragma omp atomic write
               Zijkl(iket,ibra) = hZ * Zijkl(ibra,iket);
            }
         }
      }
//      if (inner_loop)  // if we made it to the inner loop, count it.
//      {
//        n_channels_kept++;
//      }
   }
}





//void Operator::AddInverseTensorPandyaTransformation( const std::map<std::array<index_t,2>,arma::mat>&  Zbar )
void AddInverseTensorPandyaTransformation( Operator& Z, const std::map<std::array<index_t,2>,arma::mat>&  Zbar )
{
    // Do the inverse Pandya transform
   int Lambda = Z.rank_J;
//   std::vector<std::map<std::array<int,2>,arma::mat>::iterator> iteratorlist;
   std::vector<std::map<std::array<size_t,2>,arma::mat>::iterator> iteratorlist;
//   for (std::map<std::array<int,2>,arma::mat>::iterator iter= Z.TwoBody.MatEl.begin(); iter!= Z.TwoBody.MatEl.end(); ++iter) iteratorlist.push_back(iter);
   for (auto iter= Z.TwoBody.MatEl.begin(); iter!= Z.TwoBody.MatEl.end(); ++iter) iteratorlist.push_back(iter);
   int niter = iteratorlist.size();
   int hZ = Z.IsHermitian() ? 1 : -1;
    // Only go parallel if we've previously calculated the SixJs/NineJs. Otherwise, it's not thread safe.
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Lambda))
   for (int i=0; i<niter; ++i)
   {
      const auto iter = iteratorlist[i];
      int ch_bra = iter->first[0];
      int ch_ket = iter->first[1];
      arma::mat& Zijkl = iter->second;
      const TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      const TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      index_t nBras = tbc_bra.GetNumberKets();
      index_t nKets = tbc_ket.GetNumberKets();

      for (index_t ibra=0; ibra<nBras; ++ibra)
      {
         const Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         const Orbit & oi = Z.modelspace->GetOrbit(i);
         const Orbit & oj = Z.modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;
         index_t ketmin = ch_bra==ch_ket ? ibra : 0;
         for (index_t iket=ketmin; iket<nKets; ++iket)
         {
            const Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            const Orbit & ok = Z.modelspace->GetOrbit(k);
            const Orbit & ol = Z.modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;

            double commij = 0;
            double commji = 0;

            // Transform Z_ilkj
            int parity_bra_cc = (oi.l+ol.l)%2;
            int parity_ket_cc = (ok.l+oj.l)%2;
            int Tz_bra_cc = std::abs(oi.tz2+ol.tz2)/2;
            int Tz_ket_cc = std::abs(ok.tz2+oj.tz2)/2;
            int j3min = std::abs(int(ji-jl));
            int j3max = ji+jl;
            for (int J3=j3min; J3<=j3max; ++J3)
            {
              index_t ch_bra_cc = Z.modelspace->GetTwoBodyChannelIndex(J3,parity_bra_cc,Tz_bra_cc);
              const TwoBodyChannel_CC& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
              index_t nbras = tbc_bra_cc.GetNumberKets();
              index_t indx_il = tbc_bra_cc.GetLocalIndex( std::min(i,l), std::max(i,l) );
              int j4min = std::max( std::abs(int(jk-jj)), std::abs(J3-Lambda) );
              int j4max = std::min( int(jk+jj), J3+Lambda );
              for (int J4=j4min; J4<=j4max; ++J4)
              {
                 index_t ch_ket_cc = Z.modelspace->GetTwoBodyChannelIndex(J4,parity_ket_cc,Tz_ket_cc);
                 const TwoBodyChannel_CC& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
                 index_t nkets = tbc_ket_cc.GetNumberKets();
                 index_t indx_kj = tbc_ket_cc.GetLocalIndex( std::min(j,k), std::max(j,k) );

                  double ninej = Z.modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
                  if (std::abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
                  double tbme = 0;
                  index_t ch_lo = std::min(ch_bra_cc,ch_ket_cc);
                  index_t ch_hi = std::max(ch_bra_cc,ch_ket_cc);
                  auto zbar_iter = Zbar.find({ch_lo,ch_hi});
                  if (zbar_iter == Zbar.end()) continue;
                  const auto& Zmat = zbar_iter->second;

                  if (ch_bra_cc <= ch_ket_cc)
                  {
                    if (i<=l) tbme = Zmat( indx_il ,indx_kj+(k>j?nkets:0) );
                    else      tbme = Zmat( indx_il ,indx_kj+(k>j?0:nkets) ) * hZ * Z.modelspace->phase(J3-J4 + ji+jj+jk+jl);
                  }
                  else
                  {
                      if(k<=j) tbme = Zmat(indx_kj, indx_il+(i>l?nbras:0) ) * hZ * Z.modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                      else     tbme = Zmat(indx_kj, indx_il+(i>l?0:nbras) ) * Z.modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
                  }
                  commij += hatfactor * Z.modelspace->phase(jj+jl+J2+J4) * ninej * tbme ;
              }
            }

            if (i==j)
            {
              commji = commij;
            }
            else
            {
              // Transform Z_jlki
              parity_bra_cc = (oj.l+ol.l)%2;
              parity_ket_cc = (ok.l+oi.l)%2;
              Tz_bra_cc = std::abs(oj.tz2+ol.tz2)/2;
              Tz_ket_cc = std::abs(ok.tz2+oi.tz2)/2;
              j3min = std::abs(int(jj-jl));
              j3max = jj+jl;
  
              for (int J3=j3min; J3<=j3max; ++J3)
              {
                int ch_bra_cc = Z.modelspace->GetTwoBodyChannelIndex(J3,parity_bra_cc,Tz_bra_cc);
                const TwoBodyChannel_CC& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
                int nbras = tbc_bra_cc.GetNumberKets();
                int indx_jl = tbc_bra_cc.GetLocalIndex( std::min(j,l), std::max(j,l) );
                int j4min = std::max( std::abs(int(jk-ji)), std::abs(J3-Lambda) );
                int j4max = std::min( int(jk+ji), J3+Lambda );
                for (int J4=j4min; J4<=j4max; ++J4)
                {
                   int ch_ket_cc = Z.modelspace->GetTwoBodyChannelIndex(J4,parity_ket_cc,Tz_ket_cc);
                   const TwoBodyChannel_CC& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
                   int nkets = tbc_ket_cc.GetNumberKets();
                   int indx_ki = tbc_ket_cc.GetLocalIndex( std::min(k,i), std::max(k,i) );
                    double ninej = Z.modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);
                    if (std::abs(ninej) < 1e-8) continue;
                    double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
  
                    index_t ch_lo = std::min(ch_bra_cc,ch_ket_cc);
                    index_t ch_hi = std::max(ch_bra_cc,ch_ket_cc);
                    auto zbar_iter = Zbar.find({ch_lo,ch_hi});
                    if (zbar_iter == Zbar.end()) continue;
                    const auto& Zmat = zbar_iter->second;
                    double tbme = 0;
                    if (ch_bra_cc <= ch_ket_cc)
                    {
                      if (j<=l) tbme = Zmat( indx_jl ,indx_ki+(k>i?nkets:0) );
                      else      tbme = Zmat( indx_jl ,indx_ki+(k>i?0:nkets) ) * hZ * Z.modelspace->phase(J3-J4 + ji+jj+jk+jl);
                    }
                    else
                    {
                        if(k<=i) tbme = Zmat(indx_ki, indx_jl+(j>l?nbras:0) ) * hZ * Z.modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                        else     tbme = Zmat(indx_ki, indx_jl+(j>l?0:nbras) ) * Z.modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
                    }
  
  
                      commji += hatfactor * Z.modelspace->phase(ji+jl+J2+J4) * ninej * tbme ;
                }
              }
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            Zijkl(ibra,iket) +=  (commij - Z.modelspace->phase(ji+jj-J1)*commji) / norm;
            if (ch_bra==ch_ket) Zijkl(iket,ibra) = hZ * Zijkl(ibra,iket);
         }
      }
   }
}




//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.     
//                                             
//   |          |      |          |           
//   |     __Y__|      |     __X__|            
//   |    /\    |      |    /\    |
//   |   (  )   |  _   |   (  )   |            
//   |____\/    |      |____\/    |            
//   |  X       |      |  Y       |            
//           
//            
// -- This appears to agree with Nathan's results
//
/// Calculates the part of \f$ [X_{(2)},\mathbb{Y}^{\Lambda}_{(2)}]_{ijkl} \f$ which involves ph intermediate states, here indicated by \f$ \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} \f$
/// \f[
/// \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} = \sum_{abJ_3J_4}(n_a-n_b) \hat{J_1}\hat{J_2}\hat{J_3}\hat{J_4}
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_l  &  J_3 \\
///  j_j  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
/// \left( \bar{X}^{J3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
///  -(-1)^{j_i+j_j-J_1}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_l  &  J_3 \\
///  j_i  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
/// \left( \bar{X}^{J_3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
/// \right]
/// \f]
/// This is implemented by defining an intermediate matrix
/// \f[
/// \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{i\bar{l}k\bar{j}} \equiv \sum_{ab}(n_a\bar{n}_b)
/// \left[ \left( \bar{X}^{J3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
/// -\left( \bar{X}^{J_3}_{i\bar{l}b\bar{a}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{b\bar{a}k\bar{j}} - 
///    \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}b\bar{a}}\bar{X}^{J_4}_{b\bar{a}k\bar{j}} \right)\right]
/// \f]
/// The Pandya-transformed matrix elements are obtained with DoTensorPandyaTransformation().
/// The matrices \f$ \bar{X}^{J_3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} \f$
/// and \f$ \bar{\mathbb{Y}}^{J_4J_3\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_3}_{a\bar{b}k\bar{j}} \f$
/// are related by a Hermitian conjugation, which saves two matrix multiplications, provided we
/// take into account the phase \f$ (-1)^{J_3-J_4} \f$ from conjugating the spherical tensor.
/// The commutator is then given by
/// \f[
/// \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} = \sum_{J_3J_4} \hat{J_1}\hat{J_2}\hat{J_3}\hat{J_4}
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_l  &  J_3 \\
///  j_j  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
///  \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{i\bar{l}k\bar{j}}
///  -(-1)^{j_i+j_j-J}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_l  &  J_3 \\
///  j_i  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
///  \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{j\bar{l}k\bar{i}}
///  \right]
///  \f]
///
//void Operator::comm222_phst( Operator& Y, Operator& Z ) 
//void Operator::comm222_phst( const Operator& X, const Operator& Y ) 
void comm222_phst( const Operator& X, const Operator& Y, Operator& Z ) 
{

   int hX = X.IsHermitian() ? 1 : -1;
   int hY = Y.IsHermitian() ? 1 : -1;

   double t_start = omp_get_wtime();
   // We reuse Xt_bar multiple times, so it makes sense to calculate them once and store them in a deque.
   std::deque<arma::mat> Xt_bar_ph = InitializePandya( Z, Z.nChannels, "transpose"); // We re-use the scalar part multiple times, so there's a significant speed gain for saving it
   std::map<std::array<index_t,2>,arma::mat> Y_bar_ph;
   DoPandyaTransformation(X, Xt_bar_ph, "transpose" );
   X.profiler.timer["DoTensorPandyaTransformation"] += omp_get_wtime() - t_start;


   t_start = omp_get_wtime();
   // Construct the intermediate matrix Z_bar
   // First, we initialize the map Z_bar with empty matrices
   // to avoid problems in the parallel loop -- (do we even want a parallel loop here?)
   std::map<std::array<index_t,2>,arma::mat> Z_bar;

   t_start = omp_get_wtime();
   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup(Z.GetJRank(), Z.GetTRank(), Z.GetParity() );
   X.profiler.timer["PandyaLookup"] += omp_get_wtime() - t_start;

   std::vector<index_t> ybras;
   std::vector<index_t> ykets;
   for (auto ich_bra : Z.modelspace->SortedTwoBodyChannels_CC)
   {
     int n_rows = Z.modelspace->GetTwoBodyChannel_CC(ich_bra).GetNumberKets();
     for (auto ich_ket : Z.modelspace->SortedTwoBodyChannels_CC)
     {
       if (ich_bra>ich_ket) continue;
//       if (pandya_lookup.at({(int)ich_bra,(int)ich_ket})[0].size()<1) continue;
       if (pandya_lookup.at({ich_bra,ich_ket})[0].size()<1) continue;
         int n_cols = 2*Z.modelspace->GetTwoBodyChannel_CC(ich_ket).GetNumberKets();
         ybras.push_back(ich_bra);
         ykets.push_back(ich_ket);
         Z_bar[{ich_bra,ich_ket}] = arma::mat(n_rows,n_cols);
     }
   }
   int counter = ybras.size();
//  std::cout << "Done allocating" << std::endl;


   X.profiler.timer["Allocate Z_bar_tensor"] += omp_get_wtime() - t_start;

   t_start = omp_get_wtime();

   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->tensor_transform_first_pass.at(Z.GetJRank()))
   #endif
   for(int i=0;i<counter;++i)
   {
      index_t ch_bra_cc = ybras[i];
      index_t ch_ket_cc = ykets[i];
//      const auto plookup = pandya_lookup.find({(int)ch_bra_cc,(int)ch_ket_cc});
      const auto plookup = pandya_lookup.find({ch_bra_cc,ch_ket_cc});
      if ( plookup == pandya_lookup.end() or plookup->second[0].size()<1 )
      {
       continue;
      }

      const auto& tbc_bra_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      const auto& tbc_ket_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
      int Jbra = tbc_bra_cc.J;
      int Jket = tbc_ket_cc.J;

      arma::mat YJ1J2;
      arma::mat YJ2J1;
      const auto& XJ1 = Xt_bar_ph[ch_bra_cc];
      const auto& XJ2 = Xt_bar_ph[ch_ket_cc];

      arma::uvec kets_ph = arma::join_cols( tbc_ket_cc.GetKetIndex_hh(), tbc_ket_cc.GetKetIndex_ph() );
      arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );

      DoTensorPandyaTransformation_SingleChannel(Y, YJ1J2, ch_bra_cc, ch_ket_cc);
      if (ch_bra_cc==ch_ket_cc)
      {
         YJ2J1 = YJ1J2;
      }
      else
      {
         DoTensorPandyaTransformation_SingleChannel(Y, YJ2J1, ch_ket_cc, ch_bra_cc);
      }

      int flipphaseY = hY * Z.modelspace->phase( Jbra - Jket ) ;
      // construct a matrix of phases (-1)^{k+j+p+h} used below to generate X_phkj for k>j
      arma::mat PhaseMatXJ2( tbc_ket_cc.GetNumberKets(), kets_ph.size(), arma::fill::ones) ;
      arma::mat PhaseMatYJ1J2( bras_ph.size(), tbc_ket_cc.GetNumberKets(), arma::fill::ones) ;
      for ( index_t iket=0;iket<(index_t)tbc_ket_cc.GetNumberKets();iket++)
      {
        const Ket& ket = tbc_ket_cc.GetKet(iket);
        if ( Z.modelspace->phase((ket.op->j2+ket.oq->j2)/2)<0)
        {
           PhaseMatXJ2.row(iket) *=-1;
           PhaseMatYJ1J2.col(iket) *=-1;
        }
      }
      for (index_t iph=0;iph<kets_ph.size();iph++)
      {
        const Ket& ket_ph = tbc_ket_cc.GetKet( kets_ph[iph] );
        if ( Z.modelspace->phase((ket_ph.op->j2+ket_ph.oq->j2)/2)<0)   PhaseMatXJ2.col(iph) *=-1;
      }
      for (index_t iph=0;iph<bras_ph.size();iph++)
      {
        const Ket& bra_ph = tbc_bra_cc.GetKet( bras_ph[iph] );
        if ( Z.modelspace->phase((bra_ph.op->j2+bra_ph.oq->j2)/2)<0)   PhaseMatYJ1J2.row(iph) *=-1;
      }
      PhaseMatYJ1J2 *= flipphaseY;



//                J2                       J1         J2                       J2          J2
//             k<=j     k>=j                hp  -ph    hp   ph                 k<=j       k<=j
//      J1   [       |       ]       J1   [           |          ]         [hp        |ph        ]
//     i<=j  [  Zbar | Zbar  ]  =   i<=j  [   Xbar    | -Ybar    ]   * J1  [   Ybar   |   Ybar'  ]
//           [       |       ]            [           |          ]         [ph        |hp        ]      where Ybar'_phkj = Ybar_hpkj * (-1)^{p+h+k+j}*(-1)^{J1-J2}*hY
//                                                                         [----------|----------]       and
//                                                                     J2  [hp        |ph        ]            Xbar'_phkj = Xbar_hpkj * (-1)^{p+h+k+j}*hX
//                                                                         [   Xbar   |   Xbar'  ]
//                                                                         [-ph       |-hp       ]
//    
//
      int halfncx2 = XJ2.n_cols/2;
      int halfnry12 = YJ1J2.n_rows/2;
      auto& Zmat = Z_bar.at({ch_bra_cc,ch_ket_cc});

      arma::mat Mleft = join_horiz( XJ1,  -flipphaseY * YJ2J1.t() );
      arma::mat Mright = join_vert( join_horiz( YJ1J2 ,  join_vert( YJ1J2.tail_rows(halfnry12)%PhaseMatYJ1J2 ,
                                                                    YJ1J2.head_rows(halfnry12)%PhaseMatYJ1J2  )   ), 
                                  hX*join_vert( XJ2,    join_horiz(   XJ2.tail_cols(halfncx2)%PhaseMatXJ2 ,
                                                                      XJ2.head_cols(halfncx2)%PhaseMatXJ2     )   ).t() );

      Zmat = Mleft * Mright;

   }
   X.profiler.timer["Build Z_bar_tensor"] += omp_get_wtime() - t_start;


   t_start = omp_get_wtime();
   AddInverseTensorPandyaTransformation(Z, Z_bar);

   X.profiler.timer["InverseTensorPandyaTransformation"] += omp_get_wtime() - t_start;

   Z.modelspace->tensor_transform_first_pass.at( Z.GetJRank() ) = false;

}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Begin Scalar-Dagger Commutator Terms////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////





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
// This could be modified to be more efficent since we only need one column from Y to go to one column in Z.
// But this is so far from being the bottleneck that it makes more sense to be clear.
void comm211sd( const Operator& X, const Operator& Y, Operator& Z )
{
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
// Sum_ab Sum_J (na n`b - n`a nb) X_ab Y_bia  hat(J)^2/hat(lambda)^2
// or
// Sum_ab Sum_J (2J+1)/(2lambda+1)  (na n`b) ( X_ab Y_bia - X_ba Y_aib )
//
// Adapted from
//    void Operator::comm121ss( const Operator& X, const Operator& Y) 
void comm231sd( const Operator& X, const Operator& Y, Operator& Z) 
{
   if (Y.legs<3) return;
   index_t norbits = Z.modelspace->GetNumberOrbits();
   index_t Q = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
//   for (index_t i=0;i<norbits;++i)
//   #pragma omp parallel for 
   for (auto i : Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
   {
          for (auto& a : Z.modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = Z.modelspace->GetOrbit(a);
             for (index_t b=0; b<norbits; ++b)
             {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ); // despite what the name suggests, nanb is n_a * (1-n_b).
                if (std::abs(nanb)<OCC_CUT) continue;
                  Z.OneBody(i,Q) += (ob.j2+1) * nanb *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,Q) ;  // Is this still the right way to do this? Do we need to worry about normalization? (It looks ok).
                  Z.OneBody(i,Q) -= (oa.j2+1) * nanb *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,Q) ;  // GetTBMEmonopole returns unnormalized TBME summed over J times (2J+1)/((2ji+1)*(2jj+1))
             }
          }
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
//  Zi = sum 1/2 1/(2j_i +1)  sum_abcJ (2J+1)  (na nb n_c - n_a n_b nc) AJ_ciab BJ_abc
//  note that in Hergert et al Phys. Rep. 621 165  this formula has an errant extra factor of 1/2
//
//  Adapted from
//    void Operator::comm221ss( const Operator& X, const Operator& Y) 
// THIS IS NOT USED
/*
void comm431sd( const Operator& X, const Operator& Y, Operator& Z)
{

   int norbits = Z.modelspace->GetNumberOrbits();

   // I think the static call was an attempt to reduce memory fragmentation, but this should be tested more thoroughly.
   static TwoBodyME Mpp = Z.TwoBody;  // Working directly with the matrix, not the setter/getters. The matrix should store the normalized matrix element.
   static TwoBodyME Mhh = Z.TwoBody;
//   TwoBodyME Mpp = Z.TwoBody;
//   TwoBodyME Mhh = Z.TwoBody;

   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo.
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
//   #ifndef OPENBLAS_NOUSEOMP
//   #pragma omp parallel for schedule(dynamic,1)
//   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

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
      // Fractional occupation of hole states means we need to do a bit more work.
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;

   } //for ch

   index_t Q = Z.GetQSpaceOrbit();
   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double ciQJ = 0;
      for (int ch=0;ch<Z.nChannels;++ch)
      {
         TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
         double Jfactor = (2*tbc.J+1.0);
         // Sum c over holes and include the nbar_a * nbar_b terms
         for (auto& c : Z.modelspace->holes)
         {
            Orbit& oc = Z.modelspace->GetOrbit(c);
            ciQJ += Jfactor * oc.occ     * Mpp.GetTBME(ch,c,i,c,Q);
            ciQJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,Q);
         }
         // Sum c over particles and include the n_a * n_b terms
         for (auto& c : Z.modelspace->particles)
         {
            ciQJ += Jfactor * Mhh.GetTBME(ch,c,i,c,Q);
         }
       }
       Z.OneBody(i,Q) += ciQJ /(oi.j2+1.0); 
   }
}
*/




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
void comm413_233sd( const Operator& X, const Operator& Y, Operator& Z) 
{
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;

   index_t Qorbit = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Qorbit);
   int norb = Z.modelspace->GetNumberOrbits();

   int n_nonzero = Z.modelspace->SortedTwoBodyChannels.size();
//   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = Z.modelspace->GetOrbit(i);
         Orbit& oj = Z.modelspace->GetOrbit(j);
         double norm_ij = (i==j) ? 1.0/SQRT2 : 1.0;
         for ( int k=0; k<norb; ++k )
         {
           Orbit& ok = Z.modelspace->GetOrbit(k);
           if (not tbc.CheckChannel_ket(&ok,&oQ)) continue;

           double norm_kQ = (k==Qorbit) ? 1.0/SQRT2 : 1.0;
           double cijk = 0;

           for ( int a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
           {
//             cijk += X1(i,a) * Y.TwoBody.GetTBME_norm(ch,ch,a,j,k,Qorbit);
             cijk += X1(i,a) * Y.TwoBody.GetTBME(ch,ch,a,j,k,Qorbit);
           }
           for ( int a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
           {
//             cijk += X1(j,a) * Y.TwoBody.GetTBME_norm(ch,ch,i,a,k,Qorbit);
             cijk += X1(j,a) * Y.TwoBody.GetTBME(ch,ch,i,a,k,Qorbit);
           }
           for ( int a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
           {
//             cijk += X1(a,k) * Y.TwoBody.GetTBME_norm(ch,ch,i,j,a,Qorbit); // should this have a minus sign?
//             cijk -= X1(a,k) * Y.TwoBody.GetTBME_norm(ch,ch,i,j,a,Qorbit); 
             cijk -= X1(a,k) * Y.TwoBody.GetTBME(ch,ch,i,j,a,Qorbit); 
           }

           for ( int a : Y.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
           {
//             cijk += Y1(a,Qorbit) * X.TwoBody.GetTBME_norm(ch,ch,i,j,k,a);   // This determines the normalization for the (adagger adagger a) term.
             cijk += Y1(a,Qorbit) * X.TwoBody.GetTBME(ch,ch,i,j,k,a);   // This determines the normalization for the (adagger adagger a) term.
           }

//           Z.TwoBody.SetTBME(ch,ch,i,j,k,Qorbit, cijk);   // SetTBME  directly sets the matrix element, i.e. it assumes a normalized matrix element.

           cijk *= (norm_ij * norm_kQ);  // We normalize like a scalar TBME so that the setter/getters make sense. We'll fix the wrong |kQ> normalization when writing to file.
           Z.TwoBody.SetTBME(ch,ch,i,j,k,Qorbit, cijk);   // SetTBME  directly sets the matrix element, i.e. it assumes a normalized matrix element.
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
// THIS IS NOT USED
/*
void comm433sd_pphh( const Operator& X, const Operator& Y, Operator& Z)  
{

//   static TwoBodyME Mpp = Z.TwoBody;
//   static TwoBodyME Mhh = Z.TwoBody;
   TwoBodyME Mpp = Z.TwoBody;
   TwoBodyME Mhh = Z.TwoBody;

   double t = omp_get_wtime();
   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo.
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

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


      // The two body part
      OUT += Matrixpp - Matrixhh;
   } //for ch
   Z.profiler.timer["sd_pphh TwoBody bit"] += omp_get_wtime() - t;
}
*/







void ConstructDaggerMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, TwoBodyME& Mpp, TwoBodyME& Mhh)
{
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);

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


   } //for ch

}


/// Since comm443_pp_hhsd() and comm431sd() both require the construction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
/// We do the same thing in the standard scalar-scalar commutators.
void comm433_pp_hh_431sd( const Operator& X, const Operator& Y, Operator& Z )  
{


//   static TwoBodyME Mpp = Z.TwoBody;
//   static TwoBodyME Mhh = Z.TwoBody;
   TwoBodyME Mpp = Z.TwoBody;
   TwoBodyME Mhh = Z.TwoBody;

   double t = omp_get_wtime();
   ConstructDaggerMpp_Mhh( X, Y, Z, Mpp, Mhh);

   Z.TwoBody += Mpp;
   Z.TwoBody -= Mhh;

   Z.profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;

   // Now, the one body part
   t = omp_get_wtime();
   int Q = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
   auto ilist = Z.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2});
   int ni = ilist.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int i_ind=0; i_ind<ni; ++i_ind)
   {
      int i = ilist[i_ind];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      double cijJ = 0;
      for (int ch=0;ch<Z.nChannels;++ch)
      {
         TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
         double Jfactor = (2*tbc.J+1.0);
         // Sum c over holes and include the nbar_a * nbar_b terms
         for (auto& c : Z.modelspace->holes)
         {
            Orbit& oc = Z.modelspace->GetOrbit(c);
            cijJ += Jfactor * oc.occ * Mpp.GetTBME(ch,c,i,c,Q);     // We use the GetTBME, which returns an unnormalized matrix element, as required.
            cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,Q); 
         }
         // Sum c over particles and include the n_a * n_b terms
         for (auto& c : Z.modelspace->particles)
         {
            cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,Q);
         }
      }
      Z.OneBody(i,Q) += cijJ /(oi.j2+1.0);   // The factor of 1/2 in the formula is absorbed by the fact that the mat-mult only sums a<=b.
   } // for i
   Z.profiler.timer["pphh One Body bit"] += omp_get_wtime() - t;
}




//*****************************************************************************************
// [X^(4),Y^(3)]^(3)]  ph piece, the slow way
//                                             
//   |           |       |           |           
//   |      _(Y)_|       |_(X)_      |            
//   |     /\            |    /\     |
//   |    (  )      __   |   (  )    |            
//   |_(X)_\/            |    \/_(Y)_|            
//   |                   |            
//
//  Straightfoward and very slow implementation, only used for unit testing the
//  faster mat-mult implementation.  The two implementations agree as of Nov 2, 2018 - SRS 
void comm433sd_ph_dumbway( const Operator& X, const Operator& Y, Operator& Z)
{

   int norb = Z.modelspace->GetNumberOrbits();
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
   index_t Q = Y.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
   double jQ = 0.5*oQ.j2;
   auto kets_ph = Z.modelspace->KetIndex_ph;
   for ( auto khh : Z.modelspace->KetIndex_hh ) // ph kets needs to include hh list in case there are fractionally occupied states
   {
     kets_ph.push_back(khh);
   }

   for (int ich=0; ich<nch; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      size_t nkets = tbc.GetNumberKets();
      for (size_t ibra=0; ibra<nkets; ibra++)
      {
        auto& bra = tbc.GetKet(ibra);
        std::vector<index_t> i_cases = { bra.p, bra.q };
        std::vector<index_t> j_cases = { bra.q, bra.p };
        std::vector<int> ijsign_cases = { +1, bra.Phase(J)};  // This accounts for the 1 - (-1)^(ji +jj-J) Pij  factor.
        double norm_ij = (bra.p==bra.q) ? 1.0 / SQRT2 : 1.0;
        for (int ijcase=0; ijcase<=1; ijcase++)
        {
          index_t i = i_cases[ijcase];
          index_t j = j_cases[ijcase];
          int ijsign = ijsign_cases[ijcase];
          Orbit& oi = Z.modelspace->GetOrbit(i);
          Orbit& oj = Z.modelspace->GetOrbit(j);
          double ji = 0.5* oi.j2;
          double jj = 0.5* oj.j2;
          for ( size_t k=0; k<norb; k++)
          {
            Orbit& ok = Z.modelspace->GetOrbit(k);
            if ( not tbc.CheckChannel_ket(&ok,&oQ) ) continue;
            double zijk = 0.;
            double norm_kQ = (k==Q) ? 1.0 / SQRT2 : 1.0;
            double jk = 0.5* ok.j2;


            int Jprime_min = std::max(  std::abs(ji-jQ), std::abs(jj-jk) );
            int Jprime_max = std::min(  ji+jQ, jj+jk  );
            for (int Jprime=Jprime_min; Jprime<=Jprime_max; ++Jprime)
            {
              double sixjprime =  Z.modelspace->GetSixJ(ji,jj,J,jk,jQ,Jprime);
              if (std::abs(sixjprime) < 1e-8) continue;

              double XYprod = 0;

              for (auto& iketab : kets_ph )
              {
                auto& ketab = Z.modelspace->GetKet(iketab);
                std::vector<index_t> ab_cases = { ketab.p, ketab.q };
                for (int abcase = 0; abcase<=1; abcase++)
                {
                  index_t a  = ab_cases[abcase];
                  index_t b  = ab_cases[1-abcase];

                  Orbit& oa = Z.modelspace->GetOrbit(a);
                  Orbit& ob = Z.modelspace->GetOrbit(b);
                  double ja = 0.5*oa.j2;
                  double jb = 0.5*ob.j2;
                  double nanb = oa.occ - ob.occ;


                  int JA_min = std::max(  std::abs(ja-jj), std::abs(jk-jb)  );
                  int JA_max = std::min( ja+jj,  jk+jb );

                  int JB_min = std::max(  std::abs(ja-jQ), std::abs(ji-jb)  );
                  int JB_max = std::min( ja+jQ,  ji+jb );
                  double matelX = 0;
                  double matelY = 0;
                  for (int JA = JA_min; JA<=JA_max; ++JA)
                  {
                    double sixjA = Z.modelspace->GetSixJ( ja, jb, Jprime, jk, jj, JA );
                    matelX -= (2*JA+1) * sixjA * X.TwoBody.GetTBME_J(JA,a,j,k,b);   // GetTBME_J returns an un-normalized matrix element
                  }
                  for (int JB = JB_min; JB<=JB_max; ++JB)
                  {
                    double sixjB = Z.modelspace->GetSixJ( ja, jb, Jprime, ji, jQ, JB );
                    matelY -= (2*JB+1) * sixjB * Y.TwoBody.GetTBME_J(JB,i,b,a,Q);   // GetTBME_J returns an un-normalized matrix element
                  }
                  XYprod += nanb * matelX * matelY;
                }
              } // loop over ab cases
               zijk -= ijsign * (2*Jprime+1) * sixjprime * XYprod;
            }  // loop over ab kets
 


            Z.TwoBody.AddToTBME_J(J,i,j,k,Q,  zijk * norm_ij*norm_kQ);   // AddToTBME_J assumes a normalized matrix element.
          }  // loop over k 
        } // loop over ij cases
      } // for ibra
   } // for ich


}




//*****************************************************************************************
// [X^(4),Y^(3)]^(3)]  ph piece
//                                             
//   |           |       |           |           
//   |      _(Y)_|       |_(X)_      |            
//   |     /\            |    /\     |
//   |    (  )      __   |   (  )    |            
//   |_(X)_\/            |    \/_(Y)_|            
//   |                   |            
//
// For formula and details of implementation, see Operator::comm222_phss. The only change
// made here to accommodate a dagger operator is to replace XY - YX  with -YX. This
// ensures that we don't include contributions from the Qspace orbit to X.
// Adapted from
//    void Operator::comm222_phss( const Operator& X, const Operator& Y ) 


void comm433sd_ph( const Operator& X, const Operator& Y, Operator& Z)  
{

   double t_start = omp_get_wtime();

   int hx = X.IsHermitian() ? 1 : -1; // assuming X is either hermitian or antihermitian. I'm sure this will come back to bite me some day.

   // Construct the intermediate matrix Z_bar and fill it with zeros.
   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup( Z.GetJRank(), Z.GetTRank(), Z.GetParity() );
   size_t nch = Z.modelspace->SortedTwoBodyChannels_CC.size();
   size_t norb = Z.modelspace->GetNumberOrbits();
   std::deque<arma::mat> Z_bar (Z.nChannels );
   std::vector<bool> lookup_empty(Z.nChannels,true);
   for (size_t ich=0;ich<nch;++ich)
   {
      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      index_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros( norb, 2*nKets_cc );  // Z_iQ`kj`   we want both orderings if k and j, but Q is fixed.
      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
   }


   // loop over channels for Z_bar
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   #endif
   for (size_t ich=0; ich<nch; ++ich )
   {
      if (lookup_empty.at(ich)) continue;
      size_t ch = Z.modelspace->SortedTwoBodyChannels_CC.at(ich);
      const TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      size_t nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

      arma::mat Y_bar_ph;
      arma::mat X_bar_ph;

      DoPandyaTransformation_SingleChannel_Dagger(Y,Y_bar_ph,ch);      // Generate  YbarJ'_iQ`ab` * (na-nb)      for a<=b and a>b
      DoPandyaTransformation_SingleChannel(X,X_bar_ph,ch,"normal");  // Generate XbarJ'_ab`kj`                  for a<=b and a>b
      auto& Zbar_ch = Z_bar.at(ch); // Marginally useful aliasing to avoid lookups...

     // Zbar will be
     // ZbarJ'_iQ`kj` = Sum_ab   (na-nb) <iQ`|YbarJ'|ab`> <ab`|XbarJ'|kj`>
     // So then we can perform an inverse Pandya transform (which is really just another Pandya transform) to get
     // ZJ_ijkQ = - Sum_J' (2J'+1) { ji  jj  J  }  ZbarJ'_iQ`kj`
     //                            { jk  jQ  J' }


     // The shape of Xt_bar_ph is  ( rows, columns) = ( 2*nph_kets,  nKets_cc)  =>    [     Xbar    ]  |  ab`, element of nph_kets, which has ph` and hp`
     //                                                                               [             ]  v
     //                                                                                    -> 
     //                                                                                    kj` is an element of kets_cc
     //
     // The shape of Y_bar_ph is   (rows, columns) = (norb, 2*nph_kets)         =>    [     Ybar    ]  | iQ`  runs over all orbits i, with Q fixed.
     //                                                                               [             ]  v
     //                                                                                    -> 
     //                                                                                    ab` is an element of nph_kets, here we have ph` and hp`
     //
     // So we should multiply Ybar * Xbar to sum over ab` and get out Zbar      =>      [     Zbar    ]  |  iQ`
     //                                                                                 [             ]  v
     //                                                                                      -> 
     //                                                                                      kj`
     //
     // Really, we need  <iQ`|Zbar|kj`> and <iQ`|Zbar|jk`> 



      if (Y_bar_ph.size()<1 or X_bar_ph.size()<1)   // for an armadillo matrix, .size()  gives the total number of elements, i.e. n_rows * n_cols
      {
        Zbar_ch = arma::zeros( norb, 2*nKets_cc );  // This seems unnecessary since we initialized things earlier... try getting rid of it and see if things break.
        continue;
      }

      // get the phases for taking the transpose
      arma::mat PhaseMat(nKets_cc, nKets_cc, arma::fill::ones );
      for (index_t iket=0;iket<nKets_cc;iket++)
      {
         const Ket& ket = tbc_cc.GetKet(iket);
         if ( Z.modelspace->phase( (ket.op->j2 + ket.oq->j2)/2 ) > 0) continue;
         PhaseMat.col( iket ) *= -1;
         PhaseMat.row( iket ) *= -1;
      }
      arma::uvec phkets = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
      auto PhaseMatX = PhaseMat.rows(phkets) * hx;

//                                           [      |     ]
//     create full Y matrix from the half:   [  Yhp | Y'ph]   where the prime indicates multiplication by (-1)^(i+j+k+l) h_x
//                                           [      |     ]   Flipping hp <-> ph and multiplying by the phase is equivalent to
//                                           [  Yph | Y'hp]   having kets |kj> with k>j.
      Zbar_ch =  Y_bar_ph * join_horiz(X_bar_ph, join_vert(   X_bar_ph.tail_rows(nph_kets)%PhaseMatX,
                                                              X_bar_ph.head_rows(nph_kets)%PhaseMatX) );

   }

   Z.profiler.timer["Build Z_bar_Dagger"] += omp_get_wtime() - t_start;
//   std::cout << "Done building Z_bar_dagger" << std::endl;

   // Perform inverse Pandya transform on Z_bar to get Z
   t_start = omp_get_wtime();
   AddInversePandyaTransformation_Dagger(Z_bar, Z);

//   Z.modelspace->scalar_transform_first_pass = false;
   Z.profiler.timer["InversePandyaTransformation_Dagger"] += omp_get_wtime() - t_start;

}






//**************************************************************************
// DAGGER VARIETY
//
//  X^J_iQ`ab` = - sum_J' { i Q J } (2J'+1) X^J'_ibaQ
//                        { a b J'}
/// A modification of the Pandya transform for the case of a dagger operator.
/// We embed the dagger operator as a scalar with a ficticious fourth index Q,
/// which is outside the Hilbert space. So our relevant transformation is
/// \f[
///  \bar{X}^{J}_{i\bar{Q}a\bar{b}} = - \sum_{J'} (2J'+1)
///  \left\{ \begin{array}{lll}
///  j_i  &  j_Q  &  J \\
///  j_a  &  j_b  &  J' \\
///  \end{array} \right\}
///  X^{J}_{ibaQ}
/// \f]
/// where the overbar indicates time-reversed orbits.
/// Since we use this for the commutator expression, the k and j orbits will be a and b,
/// i.e. the orbits we sum over with a factor (na-nb). We include the (na-nb) factor
/// in the output matrix.
void DoPandyaTransformation_SingleChannel_Dagger(const Operator& Z, arma::mat& TwoBody_CC_ph, int ch_cc)
{
   TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
   int nKets_cc = tbc_cc.GetNumberKets();
   size_t norb = Z.modelspace->GetNumberOrbits();
   arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
   int nph_kets = kets_ph.n_rows;
   int J_cc = tbc_cc.J;

//   TwoBody_CC_ph.zeros( nKets_cc, 2*nph_kets);
   TwoBody_CC_ph.zeros( norb, 2*nph_kets);   // factor 2 because we want both orderings ab and ba.

   index_t Q = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
   double jQ = oQ.j2*0.5;

   // loop over cross-coupled ph kets |ab> in this channel
   // (this is the side that gets summed over in the matrix multiplication)
   for (int iket=0; iket<nph_kets; ++iket)
   {
     Ket & ket_cc = tbc_cc.GetKet( kets_ph[iket] );
     int a = ket_cc.p;
     int b = ket_cc.q;
     Orbit & oa = Z.modelspace->GetOrbit(a);
     Orbit & ob = Z.modelspace->GetOrbit(b);
     double ja = oa.j2*0.5;
     double jb = ob.j2*0.5;

     // Here we make some lists so that we don't need to repeat code when switching (a,b) -> (b,a).
     // Instead we can just iterate over the two cases.
     std::vector<int> ab     = { a, b };
     std::vector<double>jab  = { ja, jb };
     std::vector<double>nanb = { (oa.occ-ob.occ),  (ob.occ-oa.occ) };

     // we loop over both orderings, a<b and a>b. Here, this is done by exchanging a<->b and taking a minus sign due to the (na-nb) factor.
     for ( int ab_case=0; ab_case<=1; ab_case++)
     {
       size_t ab1 = ab[ab_case];
       size_t ab2 = ab[1-ab_case];
       double jab1  = jab[ab_case];
       double jab2  = jab[1-ab_case];
       size_t indx_ab = iket + ab_case*nph_kets;

       // loop over orbits i, and check if <iQ| is in the desired channel.
       for (size_t i=0; i<norb; i++)
       {
          Orbit& oi = Z.modelspace->GetOrbit(i);
          if (not tbc_cc.CheckChannel_ket( &oi, &oQ ) ) continue;  //  <iQ|  isn't in this channel, so move along. (Note, for this check, the ordering of i,Q doesn't matter).
          size_t indx_iQ = i;   // since Q is fixed, we can label the bra <iQ| by the index i.
          double ji = oi.j2*0.5;

          int Jmin = std::max( std::abs(jab1-jQ), std::abs(ji-jab2) );
          int Jmax = std::min( jab1+jQ, ji+jab2 );
          double Xbar = 0;
          for (int J_std=Jmin; J_std<=Jmax; ++J_std)
          {
             double sixj = Z.modelspace->GetSixJ(ji,jQ,J_cc,jab1,jab2,J_std);
             if (std::abs(sixj) < 1e-8) continue;
             double tbme = Z.TwoBody.GetTBME_J(J_std,i,ab2,ab1,Q);
             Xbar -= (2*J_std+1) * sixj * tbme ;
          }
          TwoBody_CC_ph(indx_iQ, indx_ab) = Xbar * nanb[ab_case] ;
       }
     }
   }
}



// Same idea as the scalar commutator version, with enough modifications to make it no wrong (I hope).
// One big difference is that the Pandya-transformed operator Zbar has its bra indices numbered
// by orbit index rather than ket index, because it's   <iQ|Zbar|kj>  and Q is fixed.
// Otherwise, it should look pretty similar.
void AddInversePandyaTransformation_Dagger( const std::deque<arma::mat>& Zbar, Operator& Z )
{
    // Do the inverse Pandya transform
   int n_nonzeroChannels = Z.modelspace->SortedTwoBodyChannels.size();
   size_t norb = Z.modelspace->GetNumberOrbits();
   size_t Q = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Q);
   double jQ = 0.5 * oQ.j2;

    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
//   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   {
      size_t ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      size_t nKets = tbc.GetNumberKets();

      // the bra is the <ij| part of Z_ijkQ
      for (size_t ibra=0; ibra<nKets; ++ibra)
      {
        Ket & bra = tbc.GetKet(ibra);

        // Two-element lists for use below, because we need to antisymmetrize with 1-Pij, including a phase factor (-1)^{ji+jj-J)
        // The code would look like the exact same thing copy-pasted twice, with i and j exchanged in the second copy.
        std::vector<size_t> ij_switcheroo = {bra.p , bra.q};
        std::vector<int> phaseij = { +1, bra.Phase(J) }; 


        for (int ij_case=0; ij_case<=1; ij_case++)
        {
          size_t i = ij_switcheroo[ij_case];
          size_t j = ij_switcheroo[1-ij_case];
          int Pij = phaseij[ij_case];

          Orbit & oi = Z.modelspace->GetOrbit(i);
          Orbit & oj = Z.modelspace->GetOrbit(j);
          double ji = oi.j2/2.;
          double jj = oj.j2/2.;
          double norm_ij = (i==j) ? 1.0/SQRT2 : 1.0;


          for (size_t k=0; k<norb; k++)
          {
            Orbit & ok = Z.modelspace->GetOrbit(k);
            if ( not tbc.CheckChannel_ket(&ok, &oQ) ) continue; // if |kQ> doesn't live in this channel, move along.

            double jk = ok.j2/2.;
            double norm_kQ = (k==Q) ? 1.0/SQRT2  : 1.0;

            double commijk = 0;  // contribution of the 1 term 

            int parity_cc = (oi.l+oQ.l)%2;
            int Tz_cc = std::abs(oi.tz2+oQ.tz2)/2;
            int Jmin = std::max(std::abs(int(ji-jQ)),std::abs(int(jk-jj)));
            int Jmax =  std::min( ji+jQ, jk+jj );
            for (int Jprime=Jmin; Jprime<=Jmax; ++Jprime)
            {
              double sixj = Z.modelspace->GetSixJ(ji,jj,J,jk,jQ,Jprime);
              if (std::abs(sixj)<1e-8) continue;
              int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
              TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
              int nkets_cc = tbc_cc.GetNumberKets();
              size_t indx_iQ = i;
              size_t indx_kj = tbc_cc.GetLocalIndex(std::min(j,k),std::max(j,k)) +(k>j?nkets_cc:0);
              double me1 = Zbar.at(ch_cc)(indx_iQ,indx_kj);
              commijk -= (2*Jprime+1) * sixj * me1  * Pij;
            }

            commijk *= norm_ij * norm_kQ;

            Z.TwoBody.AddToTBME_J(J,i,j,k,Q,  commijk );   // AddToTBME_J assumes a normalized matrix element.
          }
        } // for ij_case
      } // for ibra, which is <ij|
   } // for ich, which is J etc
}


















} // namespace Commutator
