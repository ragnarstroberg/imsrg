
#include "Commutator.hh"
#include "Operator.hh"
#include "TwoBodyME.hh"
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



//Operator Operator::Commutator( Operator& opright)
/// Returns \f$ Z = [X,Y] \f$
Operator Commutator( const Operator& X, const Operator& Y)
{
  int jrank = std::max(X.rank_J,Y.rank_J);
  int trank = std::max(X.rank_T,Y.rank_T);
  int parity = (X.parity+Y.parity)%2;
  int particlerank = std::max(X.particle_rank,Y.particle_rank);
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
   comm433sd_ph( X, Y, Z ) ; 
   comm433_pp_hh_431sd( X, Y, Z ) ; 

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
            int ind2 = tbc.GetLocalIndex( std::min(a,j), std::max(a,j) );
            if (ind2<0 or ind2>=tbc.GetNumberKets()) continue;
            ind1_ia.push_back(a);
            ind2_aj.push_back(ind2);
            factor_ia.push_back( a>j ? flipphaseij : (a==j ? SQRT2 : 1));
         }
         if (i!=j)
         {
           for (int a : Z.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
           {
              int ind2 = tbc.GetLocalIndex( std::min(a,i), std::max(a,i) );
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
      int a = bra_cc.p;
      int b = bra_cc.q;
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
           TwoBody_CC_ph(ibra,iket_cc) = Xbar;
         }
         else // "transpose"
         {
           TwoBody_CC_ph(iket_cc,ibra) = herm * Xbar * na_nb_factor;
         }

         // Exchange (a <-> b) to account for the (n_a - n_b) term
         jmin = std::max(std::abs(jb-jd),std::abs(jc-ja));
         jmax = std::min(jb+jd,jc+ja);
         Xbar = 0;
         for (int J_std=jmin; J_std<=jmax; ++J_std)
         {
            double sixj = Z.modelspace->GetSixJ(jb,ja,J_cc,jc,jd,J_std);
            if (std::abs(sixj) < 1e-8) continue;
            double tbme = Z.TwoBody.GetTBME_J(J_std,b,d,c,a);
            Xbar -= (2*J_std+1) * sixj * tbme ;
         }
         if (orientation=="normal")
         {
           TwoBody_CC_ph(ibra+nph_kets,iket_cc) = Xbar;
         }
         else  // "transpose"
         {
           TwoBody_CC_ph(iket_cc,ibra+nph_kets) = herm * Xbar * -na_nb_factor;
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
               commij += Jhat2 * sixj * Zbar(indx_il,indx_kj) ;

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
                 commji += Jhat2 *  sixj *  Zbar(indx_ik, indx_lj) ;


              }
            }
            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            #pragma omp atomic
            Zmat(ibra,iket) += (commij - Z.modelspace->phase(jk+jl-J ) * commji) / norm;
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
            int jmin = std::max(std::abs(int(ji-jl)),std::abs(int(jk-jj)));
            int jmax = std::min(int(ji+jl),int(jk+jj));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = Z.modelspace->GetSixJ(ji,jj,J,jk,jl,Jprime);
               if (std::abs(sixj)<1e-8) continue;
               int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
               int nkets_cc = tbc_cc.GetNumberKets();
               int indx_il = tbc_cc.GetLocalIndex(std::min(i,l),std::max(i,l)) +(i>l?nkets_cc:0);
               int indx_kj = tbc_cc.GetLocalIndex(std::min(j,k),std::max(j,k)) +(k>j?nkets_cc:0);
               double me1 = Zbar.at(ch_cc)(indx_il,indx_kj);
               commij += (2*Jprime+1) * sixj * me1;
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
              jmin = std::max(std::abs(int(jj-jl)),std::abs(int(jk-ji)));
              jmax = std::min(int(jj+jl),int(jk+ji));
              for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
              {
                 double sixj = Z.modelspace->GetSixJ(jj,ji,J,jk,jl,Jprime);
                 if (std::abs(sixj)<1e-8) continue;
                 int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
                 TwoBodyChannel_CC& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
                 int nkets_cc = tbc_cc.GetNumberKets();
                 int indx_ik = tbc_cc.GetLocalIndex(std::min(i,k),std::max(i,k)) +(i>k?nkets_cc:0);
                 int indx_lj = tbc_cc.GetLocalIndex(std::min(l,j),std::max(l,j)) +(l>j?nkets_cc:0);
                 // we always have i<=k so we should always flip Z_jlki = (-1)^{i+j+k+l} Z_iklj
                 double me1 = Zbar.at(ch_cc)(indx_ik, indx_lj) ;//* modelspace->phase(ji+jj+jk+jl);
                 commji += (2*Jprime+1) *  sixj * me1;

              }
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            Z.TwoBody.GetMatrix(ch,ch)(ibra,iket) += (commij - Z.modelspace->phase(jk+jl-J ) * commji) / norm;
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
   int nch = Z.modelspace->SortedTwoBodyChannels_CC.size();
   t_start = omp_get_wtime();
   std::deque<arma::mat> Z_bar ( Z.nChannels );
   std::vector<bool> lookup_empty(Z.nChannels,true);
   for (int ich=0;ich<nch;++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      index_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros( nKets_cc, 2*nKets_cc );
      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
   }


   #ifndef OPENBLAS_NOUSEOMP
//   #pragma omp parallel for schedule(dynamic,1)
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   #endif
   for (int ich=0; ich<nch; ++ich )
   {
      if (lookup_empty.at(ich)) continue;
      int ch = Z.modelspace->SortedTwoBodyChannels_CC.at(ich);
      const TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

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
   AddInversePandyaTransformation(Z_bar, Z);

   Z.modelspace->scalar_transform_first_pass = false;
   X.profiler.timer["InversePandyaTransformation"] += omp_get_wtime() - t_start;

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
void AddInverseTensorPandyaTransformation_SingleChannel(Operator& Z, arma::mat& Zbar, int ch_bra_cc, int ch_ket_cc)
{
    // Do the inverse Pandya transform
   if (ch_bra_cc > ch_ket_cc)  // hopefully this won't happen
   {
      std::cout << "WARNING: Called Operator::AddInverseTensorPandyaTransformation_SingleChannel with ch_bra_cc > ch_ket_cc : " << ch_bra_cc << " > " << ch_ket_cc << std::endl;
      std::cout << " Skipping this channel." << std::endl;
      return;
   }
   int n_channels_kept = 0;
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
      bool inner_loop = false;

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
               inner_loop = true;
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
               inner_loop = true;
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
               inner_loop = true;
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
               inner_loop = true;
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
   std::vector<std::map<std::array<int,2>,arma::mat>::iterator> iteratorlist;
   for (std::map<std::array<int,2>,arma::mat>::iterator iter= Z.TwoBody.MatEl.begin(); iter!= Z.TwoBody.MatEl.end(); ++iter) iteratorlist.push_back(iter);
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
       if (pandya_lookup.at({(int)ich_bra,(int)ich_ket})[0].size()<1) continue;
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
      const auto plookup = pandya_lookup.find({(int)ch_bra_cc,(int)ch_ket_cc});
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
// This could be modified since we only need one column from Y to go to one column in Z.
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
// Sum_iab Sum_J X_ab Y_bia  hat(J)/hat(lambda)
//
//
// Adapted from
//    void Operator::comm121ss( const Operator& X, const Operator& Y) 
void comm231sd( const Operator& X, const Operator& Y, Operator& Z) 
{
   if (Y.legs<3) return;
   index_t norbits = Z.modelspace->GetNumberOrbits();
   index_t j = Z.GetQSpaceOrbit();
   Orbit& oj = Z.modelspace->GetOrbit(j);
//   for (index_t i=0;i<norbits;++i)
//   #pragma omp parallel for 
   for (auto i : Z.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
   {
//      Orbit &oi = modelspace->GetOrbit(i);
      // j is Q-space orbit
//      index_t jmin = Z.IsNonHermitian() ? 0 : i;
//      for (auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
//      {
//          if (j<jmin) continue; // only calculate upper triangle
          for (auto& a : Z.modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = Z.modelspace->GetOrbit(a);
             for (index_t b=0; b<norbits; ++b)
             {
                Orbit &ob = Z.modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ);
                if (std::abs(nanb)<OCC_CUT) continue;
//                if (Y.legs>1)
//                {
                  Z.OneBody(i,j) += (ob.j2+1) * nanb *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                  Z.OneBody(i,j) -= (oa.j2+1) * nanb *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,j) ;
//                }
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
void comm431sd( const Operator& X, const Operator& Y, Operator& Z)
{

   int norbits = Z.modelspace->GetNumberOrbits();

   // I think the static call was an attempt to reduce memory fragmentation, but this should be tested more thoroughly.
   static TwoBodyME Mpp = Z.TwoBody;
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
      Orbit &oi = Z.modelspace->GetOrbit(i);
//      int jmin = Z.IsNonHermitian() ? 0 : i;
//      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
//      {
//         if (j<jmin) continue;
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
void comm413_233sd( const Operator& X, const Operator& Y, Operator& Z) 
{

   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;
//   int hZ = Z.IsHermitian() ? 1 : -1;

   index_t Qorbit = Z.GetQSpaceOrbit();
   Orbit& oQ = Z.modelspace->GetOrbit(Qorbit);
   int norb = Z.modelspace->GetNumberOrbits();

   int n_nonzero = Z.modelspace->SortedTwoBodyChannels.size();
//   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = Z.modelspace->GetTwoBodyChannel(ch);
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
         Orbit& oi = Z.modelspace->GetOrbit(i);
         Orbit& oj = Z.modelspace->GetOrbit(j);
         for ( int k=0; k<norb; ++k )
         {
           Orbit& ok = Z.modelspace->GetOrbit(k);
           if (not tbc.CheckChannel_ket(&ok,&oQ)) continue;
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
   Z.profiler.timer["sd_pphh TwoBody bit"] += omp_get_wtime() - t;
}








void ConstructDaggerMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, TwoBodyME& Mpp, TwoBodyME& Mhh)
{
   int nch = Z.modelspace->SortedTwoBodyChannels.size();
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

   t = omp_get_wtime();
   // The one body part
   int j = Z.GetQSpaceOrbit();
   Orbit& oj = Z.modelspace->GetOrbit(j);
//   for (int i=0;i<norbits;++i)
   auto ilist = Z.OneBodyChannels.at({oj.l,oj.j2,oj.tz2});
   int ni = ilist.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int i_ind=0; i_ind<ni; ++i_ind)
   {
      int i = ilist[i_ind];
      Orbit &oi = Z.modelspace->GetOrbit(i);
//      int jmin = Z.IsNonHermitian() ? 0 : i;
//      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
//      {
//         if (j<jmin) continue;
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
//      } // for j
   } // for i
   Z.profiler.timer["pphh One Body bit"] += omp_get_wtime() - t;
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
void comm433sd_ph( const Operator& X, const Operator& Y, Operator& Z)  
{

   int hy = Y.IsHermitian() ? 1 : -1;
   double t_start = omp_get_wtime();


   // Construct the intermediate matrix Z_bar
   const auto& pandya_lookup = Z.modelspace->GetPandyaLookup( Z.GetJRank(), Z.GetTRank(), Z.GetParity() );
   int nch = Z.modelspace->SortedTwoBodyChannels_CC.size();
   t_start = omp_get_wtime();
   std::deque<arma::mat> Z_bar (Z.nChannels );
   std::vector<bool> lookup_empty(Z.nChannels,true);
   for (int ich=0;ich<nch;++ich)
   {
      int ch = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      index_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros( nKets_cc, 2*nKets_cc );
      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
   }


   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1) if (not Z.modelspace->scalar_transform_first_pass)
   #endif
   for (int ich=0; ich<nch; ++ich )
   {
      if (lookup_empty.at(ich)) continue;
      int ch = Z.modelspace->SortedTwoBodyChannels_CC.at(ich);
      const TwoBodyChannel& tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

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

   Z.profiler.timer["Build Z_bar"] += omp_get_wtime() - t_start;

   // Perform inverse Pandya transform on Z_bar to get Z
   t_start = omp_get_wtime();
   AddInversePandyaTransformation(Z_bar, Z);

   Z.modelspace->scalar_transform_first_pass = false;
   Z.profiler.timer["InversePandyaTransformation"] += omp_get_wtime() - t_start;

}







} // namespace Commutator
