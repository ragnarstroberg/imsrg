
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

   Z.comm211sd( X, Y );

   Z.comm413sd( X, Y );


   // NEED TO FIGURE OUT IF SYMMETRIZE IS A SANE THING TO DO WITH DAGGER OPERATORS
   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();

   profiler.timer["CommutatorScalarScalar"] += omp_get_wtime() - t_css;

}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////// BEGIN IMPLEMENTATION OF COMMUTATOR EXPRESSIONS //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////



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
//
//
void DaggerOperator::comm211sd( const Operator& X, const DaggerOperator& Y )
{
   DaggerOperator& Z = *this;
   arma::mat xy = X.OneBody * Y.OneBody;
   Z.OneBody = xy;
}


/////////////////////////////////////////////////////////////////////////////////////////


// [X^(2),Y^(3)]^(1)
//         |i             |i
//         |              |
//  (X)_ b |         a _(Y)
//    \_\  |   ---    /_/
//    a  (Y)         (X) b
//
// Sum_iab Sum_J X_ab Y_bia  hat(J)/hat(lambda)
//
void DaggerOperator::comm231sd( const Operator& X, const DaggerOperator& Y) 
{
   DaggerOperator& Z = *this;
   


}



/////////////////////////////////////////////////////////////////////////////////////////

void DaggerOperator::comm431sd( const Operator& X, const DaggerOperator& Y)
{



}



/////////////////////////////////////////////////////////////////////////////////////////

//
// [X^(4),Y^(1)]^(3)
//
// i|      j|         i|      j|
//  |--(A)--|    ---   |       | (B)
//  |       |a         |--(A)--|/a
// k|      (B)        k|       
//
// - Sum_ijk  Sum_a  X_ijak Y_a
//
//
void DaggerOperator::comm413sd( const Operator& X, const DaggerOperator& Y )
{
   DaggerOperator& Z = *this;

   int n_nonzero = modelspace->SortedTwoBodyChannels.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto& X2 = X.TwoBody;
//      auto& Y2 = Y.TwoBody.GetMatrix(ch,ch);
      auto& Z2 = Z.TwoBody;
   
      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
//         Orbit& oi = modelspace->GetOrbit(i);
//         Orbit& oj = modelspace->GetOrbit(j);

         for (int indx_kl = 0; indx_kl<npq; ++indx_kl)
         {
            Ket & ket = tbc.GetKet(indx_kl);
            int k = ket.p;
            int lam = ket.q;
//            Orbit& ok = modelspace->GetOrbit(k);
            Orbit& olam = modelspace->GetOrbit(lam);
            // HERE WOULD BE A REASONABLE PLACE TO CHECK THAT we have the appropriate lambda orbit
            // C_ijk^{Jlambda}
            double cijk = 0;
            for (int a : OneBodyChannels.at({olam.l,olam.j2,olam.tz2}) )
            {
              cijk -= X2.GetTBME(ch, i,j,k,a) * Y.OneBody(lam,a);
            }
            Z2.SetTBME(ch,i,j,k,lam,cijk);
         }

      }

   }

}


/////////////////////////////////////////////////////////////////////////////////////////

//  [X^(2),Y^(3)]^(3)
//  
//  |i     /j             |i    /j
// (X)    /               |    /
//  a\   /      ---       | (Y)
//    (Y)                 |/a |
//     |                 (X)  |
//     |k                     |k
//
//  Sum_iab Xab Yiba
//
  void comm233sd( const Operator& X, const DaggerOperator& Y) ; 




/////////////////////////////////////////////////////////////////////////////////////////


  void comm433sd_pphh( const Operator& X, const DaggerOperator& Y) ; 



/////////////////////////////////////////////////////////////////////////////////////////


  void comm433sd_ph( const Operator& X, const DaggerOperator& Y) ;


