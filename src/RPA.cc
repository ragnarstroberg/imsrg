

#include "RPA.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh"

RPA::RPA(ModelSpace& ms)
 : modelspace(&ms) 
{
};

RPA::RPA(Operator& h)
 : modelspace(h.modelspace), H(h) 
{
};

using PhysConst::SQRT2;


///  In case we want to construct the A matrix for a single channel
///  and it's more convenient to specify J,parity,Tz than the channel index.
void RPA::ConstructAMatrix(int J, int parity, int Tz, bool Isovector=false)
{
   size_t ich_CC = modelspace->GetTwoBodyChannelIndex( J, parity, Tz);
   ConstructAMatrix_byIndex(ich_CC, Isovector);
}


///  In case we want to construct the B matrix for a single channel
///  and it's more convenient to specify J,parity,Tz than the channel index.
void RPA::ConstructBMatrix(int J, int parity, int Tz, bool Isovector=false)
{
   size_t ich_CC = modelspace->GetTwoBodyChannelIndex( J, parity, Tz);
   ConstructBMatrix_byIndex(ich_CC, Isovector);
}


///
/// <ai|A|bj> = <ai^-1,J|H|bj^-1,J>
///           = deta_ab delta_ij eps_ai  + <ai^-1,J|V|bj^-1,J>
void RPA::ConstructAMatrix_byIndex(size_t ich_CC, bool Isovector=false)
{
     channel = ich_CC;
     TwoBodyChannel_CC& tbc_CC = modelspace->GetTwoBodyChannel_CC(ich_CC);
     size_t nkets_ph = tbc_CC.GetKetIndex_ph().size();
     A.zeros( nkets_ph, nkets_ph);
//     arma::mat Vbar_bjck( nkets_ph, nkets_ph, arma::fill::zeros );
     int Jph = tbc_CC.J;

     size_t I_ph = 0;
     for (auto iket_ai : tbc_CC.GetKetIndex_ph() )
     {
       Ket& ket_ai = tbc_CC.GetKet(iket_ai);
       index_t a = ket_ai.p;
       index_t i = ket_ai.q;
       std::cout << " a i = " << a << " " << i << "   spe: " << H.OneBody(a,a) << "   " << H.OneBody(i,i) << std::endl;
       double ja = 0.5*modelspace->GetOrbit(a).j2;
       double ji = 0.5*modelspace->GetOrbit(i).j2;

       int phase_ai = 1;
       int phase_ia = - AngMom::phase( ja+ji - Jph );
       if (Isovector)
       {
           phase_ai *= modelspace->GetOrbit(a).tz2;
           phase_ia *= modelspace->GetOrbit(a).tz2;
       }

       if ( ket_ai.op->occ > ket_ai.oq->occ ) // if orbit a is more occupied than orbit i, switch them
       {
         std::swap(a,i);
         std::swap(ja,ji);
         std::swap(phase_ai,phase_ia);
       }

       size_t II_ph = 0;
       for (auto iket_bj : tbc_CC.GetKetIndex_ph() )
       {
         Ket& ket_bj = tbc_CC.GetKet(iket_bj);
         index_t b = ket_bj.p;
         index_t j = ket_bj.q;

         double jb = 0.5*modelspace->GetOrbit(b).j2;
         double jj = 0.5*modelspace->GetOrbit(j).j2;

         int phase_bj = 1;
         int phase_jb = - AngMom::phase( jb+jj - Jph );
         if ( ket_bj.op->occ > ket_bj.oq->occ )
         {
           std::swap(b,j);
           std::swap(jb,jj);
           std::swap(phase_bj,phase_jb);
         }

//         double Delta_ijab = OneBody(i,i) + OneBody(j,j) - OneBody(a,a) - OneBody(b,b);
         double H1b = ( iket_ai==iket_bj ) ? H.OneBody(a,a) - H.OneBody(i,i)   :   0. ;

//         int J1min = std::max(std::abs(ja-jb),std::abs(ji-jj));
//         int J1max = std::min(ja+jb,ji+jj);
         int J1min = std::max(std::abs(ja-jj),std::abs(jb-ji));
         int J1max = std::min(ja+jj,jb+ji);
         double V_aibj = 0;
//         double tbme_bjck = 0;
//         double tbme_ckia = 0;
//
//       a\    j/        a\     /i
//        ^\_H_/^    =>   ^\_H_/v
//        ^/   \^         ^/   \v
//       b/    i\        b/     \j

         if ( AngMom::Triangle(jj,jb,Jph) and AngMom::Triangle(ji,ja,Jph))
         {
          for (int J1=J1min;J1<=J1max;++J1)  //Pandya 1: <ai`| V |bj`>_Jtot
          {
//            V_aibj -= modelspace->GetSixJ(ja,ji,Jph,jb,jj,J1)  * (2*J1 + 1) *  H.TwoBody.GetTBME_J_norm(J1,a,j,b,i);
            V_aibj -= modelspace->GetSixJ(ja,ji,Jph,jb,jj,J1)  * (2*J1 + 1) *  H.TwoBody.GetTBME_J(J1,a,j,b,i);
          }
         }
         A(I_ph,II_ph) = H1b  + V_aibj * phase_ai *phase_bj;
//         Vbar_iabj(I_ph,II_ph) = tbme_iabj * phase_ia * phase_bj / Delta_ijab;
         II_ph ++;
       }
       I_ph ++;
     }

}


///
/// <ai|B|bj> = <ai^-1bj^-1,J|H|0>
/// Suhonen (11.70)
//  B_aibj(J) = (-1)^{jb+jb+J} sqrt(1+delta_ab)(1+delta_ij)} sum_J' (-1)^J' (2J'+1)
//                 * {ja ji J } <ab; J' | V | ij; J'>
//                   {jb jj J'}
//  The sqrt factor out front can be incorporated by using un-normalized TBMEs.
void RPA::ConstructBMatrix_byIndex(size_t ich_CC, bool Isovector=false)
{
     
     channel = ich_CC;
     TwoBodyChannel_CC& tbc_CC = modelspace->GetTwoBodyChannel_CC(ich_CC);
     size_t nkets_ph = tbc_CC.GetKetIndex_ph().size();
     B.zeros( nkets_ph, nkets_ph);
//     arma::mat Vbar_bjck( nkets_ph, nkets_ph, arma::fill::zeros );
     int Jph = tbc_CC.J;

     size_t I_ph = 0;
     for (auto iket_ai : tbc_CC.GetKetIndex_ph() )
     {
       Ket& ket_ai = tbc_CC.GetKet(iket_ai);
       index_t a = ket_ai.p;
       index_t i = ket_ai.q;
       double ja = 0.5*modelspace->GetOrbit(a).j2;
       double ji = 0.5*modelspace->GetOrbit(i).j2;

       int phase_ai = 1;
       int phase_ia = - AngMom::phase( ja+ji - Jph );
       if (Isovector)
       {
           phase_ai *= modelspace->GetOrbit(a).tz2;
           phase_ia *= modelspace->GetOrbit(a).tz2;
       }

       if ( ket_ai.op->occ > ket_ai.oq->occ ) // if orbit a is more occupied than orbit i, switch them
       {
         std::swap(a,i);
         std::swap(ja,ji);
         std::swap(phase_ai,phase_ia);
       }

       size_t II_ph = 0;
       for (auto iket_bj : tbc_CC.GetKetIndex_ph() )
       {
         Ket& ket_bj = tbc_CC.GetKet(iket_bj);
         index_t b = ket_bj.p;
         index_t j = ket_bj.q;

         double jb = 0.5*modelspace->GetOrbit(b).j2;
         double jj = 0.5*modelspace->GetOrbit(j).j2;

         int phase_bj = 1;
         int phase_jb = - AngMom::phase( jb+jj - Jph );
         if ( ket_bj.op->occ > ket_bj.oq->occ )
         {
           std::swap(b,j);
           std::swap(jb,jj);
           std::swap(phase_bj,phase_jb);
         }


         int J1min = std::max(std::abs(ja-jb),std::abs(ji-jj));
         int J1max = std::min(ja+jb,ji+jj);
         double V_aibj = 0;
//
//       a\    b/        a\  /i  b\  /j
//        ^\_H_/^    =>   ^\/v_H__^\/v
//        ^/   \^         
//       i/    j\    .

	 int phase_ib = AngMom::phase( ji + jb + Jph);
         if ( AngMom::Triangle(jj,jb,Jph) and AngMom::Triangle(ji,ja,Jph))
         {
          for (int J1=J1min;J1<=J1max;++J1)  //Pandya 1: <ai`| V |bj`>_Jtot
          {
            V_aibj += AngMom::phase(J1) * (2*J1 + 1) *  modelspace->GetSixJ(ja,ji,Jph,jj,jb,J1)  * H.TwoBody.GetTBME_J(J1,a,b,i,j);
          }
         }
         B(I_ph,II_ph) =  V_aibj *phase_ib *  phase_ai *phase_bj;
         II_ph ++;
       }
       I_ph ++;
     }

}


void RPA::SolveCP()
{
   arma::vec eigvals;
   arma::mat eigvecs;
   arma::mat Adiag = arma::diagmat(A);
   arma::eig_sym(eigvals,eigvecs,Adiag);
   Energies = eigvals;
   X = eigvecs;
   Y = arma::zeros( arma::size(X) );
}

void RPA::SolveTDA()
{
   arma::vec eigvals;
   arma::mat eigvecs;
   arma::eig_sym(eigvals,eigvecs,A);
   Energies = eigvals;
   X = eigvecs;
   Y = arma::zeros( arma::size(X) );
}


void RPA::SolveRPA()
{
   arma::cx_vec eigvals;
   arma::cx_mat eigvecs;

   arma::mat AB = arma::join_vert(  arma::join_horiz( A, B) ,
		                    arma::join_horiz(-B,-A) );
   arma::eig_gen(eigvals,eigvecs,AB);

   // now deal with the fact that we get spurious copies of each eigenvector with a negative eigenvalue
   double normimag = arma::norm( arma::imag(eigvals),"fro");
   if (normimag>1e-3)
   {
      std::cout << "WARNING: non-zero imaginary part of eigenvalues! " << __FILE__ << " " << __func__ << " " << __LINE__ << std::endl;
   }
   arma::uvec positive_indices = arma::find( arma::real(eigvals) + arma::imag(eigvals) >0 );
   arma::vec tmp = arma::real( eigvals(positive_indices) );
   size_t len = positive_indices.n_rows;
   arma::mat good_vecs = arma::real ( eigvecs.cols( positive_indices) );
   arma::mat Xtmp = good_vecs.head_rows(len);
   arma::mat Ytmp = good_vecs.tail_rows(len);
   arma::vec Etmp = arma::real( eigvals(positive_indices) );
   arma::uvec ordered_indices = arma::sort_index(Etmp);
   Energies = Etmp(ordered_indices);
   X = Xtmp.cols(ordered_indices);
   Y = Ytmp.cols(ordered_indices);
   std::cout << "eigvals = " << std::endl << eigvals << std::endl;
   for (size_t mu=0; mu<len; mu++)
   {
//      std::cout << "X" << std::endl << X.col(mu) << std::endl;
      double xnorm = arma::norm( X.col(mu),"fro");
      double ynorm = arma::norm( Y.col(mu),"fro");
      double nxy = xnorm*xnorm - ynorm*ynorm ;
//      double wrongway = xnorm*xnorm + ynorm*ynorm;
      double xmax = X.col(mu).max();
      double xmin = X.col(mu).min();
      double ymax = Y.col(mu).max();
      double ymin = Y.col(mu).min();
//      double y2 = Y.col(mu).t()*Y.col(mu);
//      double nxy = x2 - y2;
      std::cout << " mu = " << mu << "  norm = " << nxy  << "   Xmin/max " << xmin << " " << xmax << "   Ymin/max " << ymin << " " << ymax << "  E = " << Energies(mu) << "  vs  "<< Etmp(mu)  << std::endl;
//      if ( std::abs(xmin)>std::abs(xmax)   and xmin<0)
//      {
//          X.col(mu) *= -1;
//          Y.col(mu) *= -1;
//      }
      X.col(mu) /= sqrt(nxy);
      Y.col(mu) /= sqrt(nxy);
   }

//   Energies = arma::real( eigvals(positive_indices) );
}



/// See Ring and Schuck 8.93
/// This has not been benchmarked in any sense.
double RPA::GetEgs()
{
  double Erpa = 0;

  size_t nch = modelspace->GetNumberTwoBodyChannels_CC();
  for ( size_t ch=0; ch<nch; ch++)
  {
     ConstructAMatrix_byIndex( ch );
     ConstructBMatrix_byIndex( ch );
     SolveRPA();
     int nbosons = Energies.size();
     for (int mu=0; mu<nbosons; mu++)
     {
//       Erpa -= Energies % arma::vecnorm( Y ) % arma::vecnorm( Y );
       Erpa -= Energies(mu) * arma::dot( Y.col(mu).t() ,  Y.col(mu) );
     }
  }

  return Erpa;

}


// Suhonen 11.223
double RPA::TransitionToGroundState( Operator& OpIn, size_t mu )
{

   double mat_el =0;
   
   int lambda = OpIn.GetJRank();
   size_t ich_CC = modelspace->GetTwoBodyChannelIndex( lambda, OpIn.GetParity(), OpIn.GetTRank());
   TwoBodyChannel_CC& tbc_CC = modelspace->GetTwoBodyChannel_CC(ich_CC);
   if ( ich_CC != channel )
   {
     std::cout << "Uh oh. OpIn has JPT " << lambda << " " << OpIn.GetParity() << " " << OpIn.GetTRank() << "   but A/B matrices computed for channel " << channel << std::endl;
     return 0;
   }
//   int Jph = tbc_CC.J;

   size_t I_mi=0;
   for (auto iket_mi : tbc_CC.GetKetIndex_ph() )
   {
      Ket& ket_mi = tbc_CC.GetKet(iket_mi);
      index_t m = ket_mi.p;
      index_t i = ket_mi.q;
//      double jm = 0.5*modelspace->GetOrbit(m).j2;
//      double ji = 0.5*modelspace->GetOrbit(i).j2;

      mat_el += OpIn.OneBody(m,i) * (  AngMom::phase(lambda) * X(I_mi,mu) + Y(I_mi,mu) );
//      std::cout << __func__ << "mu = " << mu << "   mi " << m << " " << i << "  X,Y,M " << X(I_mi,mu) << " " << Y(I_mi,mu) << " " << OpIn.OneBody(m,i) << std::endl;

      I_mi++;
   }


   return mat_el;
}


// See Ring and Schuck Section 9.3.5 "Effective Charges"
double RPA::PVCouplingEffectiveCharge( Operator& OpIn, size_t k, size_t l)
{
  size_t nbosons = X.n_rows;

  size_t ich_CC = modelspace->GetTwoBodyChannelIndex( OpIn.GetJRank(), OpIn.GetParity(), OpIn.GetTRank());
  TwoBodyChannel_CC& tbc_CC = modelspace->GetTwoBodyChannel_CC(ich_CC);
  size_t nkets_ph = tbc_CC.GetKetIndex_ph().size();

//      std::cout << "A is a " << A.n_rows << "x" << A.n_cols << " matrix.   "
//                << "X is a " << X.n_rows << "x" << X.n_cols << " matrix.   "
//                << "E is a " << Energies.n_rows << "x" << Energies.n_cols << " matrix.   " << std::endl;

//  size_t mu_collective = Energies.index_min();
   if ( ich_CC != channel )
   {
     std::cout << "YIKES! TROUBLE " << __FILE__ << " "<< __func__ << " " << __LINE__  << std::endl;
     std::cout << "Uh oh. OpIn has JPT " << OpIn.GetJRank() << " " << OpIn.GetParity() << " " << OpIn.GetTRank() << "   but A/B matrices computed for channel " << channel << std::endl;
     return 0;
   }

//  if ( nkets_ph != nbosons)
//  {
//      std::cout << "YIKES! TROUBLE " << __FILE__ << " "<< __func__ << " " << __LINE__ 
//                 << "  " << nkets_ph << " != " << nbosons << std::endl;
//  }
  int Jph = tbc_CC.J;

  double Teff = OpIn.OneBody(k,l);
  double ek = H.OneBody(k,k);
  double el = H.OneBody(l,l);

  double jk = 0.5*modelspace->GetOrbit(k).j2;
  double jl = 0.5*modelspace->GetOrbit(l).j2;

  for (size_t mu=0; mu<nbosons; mu++)
  {
//     if (mu == mu_collective) continue;
     double gamma_mu_kl = 0;
     double gamma_mu_lk = 0;
     double T_mu = 0;
     size_t I_mi = 0;
//     double e_mu = 0;
     for (auto iket_mi : tbc_CC.GetKetIndex_ph() )
     {
        Ket& ket_mi = tbc_CC.GetKet(iket_mi);
        index_t m = ket_mi.p;
        index_t i = ket_mi.q;
        double jm = 0.5*modelspace->GetOrbit(m).j2;
        double ji = 0.5*modelspace->GetOrbit(i).j2;
        double Omi = OpIn.OneBody(m,i);
        double Oim = OpIn.OneBody(i,m);

        int phase_mi = 1;
        int phase_im = - AngMom::phase( jm+ji - Jph );
//        int phase_im = 1;

        if ( ket_mi.op->occ > ket_mi.oq->occ ) // if orbit a is more occupied than orbit i, switch them
        {
          std::swap(m,i);
          std::swap(jm,ji);
          std::swap(phase_mi,phase_im);
//          std::swap(Omi,Oim);
        }


          // bj -> lk   and ai -> mi
         double VA_milk = 0;
         double VA_mikl = 0;

         int J1min = std::max(std::abs(jm-jk),std::abs(jl-ji));
         int J1max = std::min(jm+jk,jl+ji);

         if ( AngMom::Triangle(jk,jl,Jph) and AngMom::Triangle(ji,jm,Jph))
         {
          for (int J1=J1min;J1<=J1max;++J1)  //Pandya 1: <ai`| V |bj`>_Jtot
          {
//            VA_milk -= modelspace->GetSixJ(jm,ji,Jph,jl,jk,J1)  * (2*J1 + 1) *  H.TwoBody.GetTBME_J_norm(J1,m,k,l,i);
            VA_milk -= modelspace->GetSixJ(jm,ji,Jph,jl,jk,J1)  * (2*J1 + 1) *  H.TwoBody.GetTBME_J(J1,m,k,l,i);
//            if ( m==6 and i==0 and k==2 and l==2)
//            {
//              std::cout << " line " << __LINE__ << " mu = " << mu << "  J1=" << J1 << " :  " << AngMom::phase(J1) << " * " << (2*J1 + 1)
//                        << " * " <<  modelspace->GetSixJ(jm,ji,Jph,jl,jk,J1)<< " * "  << H.TwoBody.GetTBME_J_norm(J1,m,k,l,i)
//                        << "  => " << VA_milk << std::endl;
//            }
          }
         }
         VA_milk *= phase_mi ;


         J1min = std::max(std::abs(jm-jl),std::abs(jk-ji));
         J1max = std::min(jm+jl,jk+ji);

         if ( AngMom::Triangle(jl,jk,Jph) and AngMom::Triangle(ji,jm,Jph))
         {
          for (int J1=J1min;J1<=J1max;++J1)  //Pandya 1: <ai`| V |bj`>_Jtot
          {
//            VA_mikl -= modelspace->GetSixJ(jm,ji,Jph,jk,jl,J1)  * (2*J1 + 1) *  H.TwoBody.GetTBME_J_norm(J1,m,l,k,i);
            VA_mikl -= modelspace->GetSixJ(jm,ji,Jph,jk,jl,J1)  * (2*J1 + 1) *  H.TwoBody.GetTBME_J(J1,m,l,k,i);
          }
         }
         VA_mikl *= phase_mi ;


          // bj -> lk   and ai -> mi
         J1min = std::max(std::abs(jm-jl),std::abs(ji-jk));
         J1max = std::min(jm+jl,ji+jk);
         double VB_milk = 0;
         double VB_mikl = 0;
//
//       m\    l/        m\  /i  l\  /k
//        ^\___/^    =>   ^\/v____^\/v
//        ^/   \^         
//       i/    k\    .

	 int phase_il = AngMom::phase( ji + jl + Jph);
         if ( AngMom::Triangle(jk,jl,Jph) and AngMom::Triangle(ji,jm,Jph))
         {
          for (int J1=J1min;J1<=J1max;++J1)  //Pandya 1: <ai`| V |bj`>_Jtot
          {
            VB_milk += AngMom::phase(J1) * (2*J1 + 1) *  modelspace->GetSixJ(jm,ji,Jph,jk,jl,J1)  * H.TwoBody.GetTBME_J(J1,m,l,i,k);

          }
         }
         VB_milk *=  phase_il *  phase_mi ;

         J1min = std::max(std::abs(jm-jk),std::abs(ji-jl));
         J1max = std::min(jm+jk,ji+jl);
	 int phase_ik = AngMom::phase( ji + jk + Jph);
         if ( AngMom::Triangle(jl,jk,Jph) and AngMom::Triangle(ji,jm,Jph))
         {
          for (int J1=J1min;J1<=J1max;++J1)  //Pandya 1: <ai`| V |bj`>_Jtot
          {
            VB_mikl += AngMom::phase(J1) * (2*J1 + 1) *  modelspace->GetSixJ(jm,ji,Jph,jl,jk,J1)  * H.TwoBody.GetTBME_J(J1,m,k,i,l);
          }
         }
         VB_mikl *=  phase_ik *  phase_mi ;



////         B(ai,bj)
////	     int phase_ib = AngMom::phase( ji + jb + Jph);
////         if ( AngMom::Triangle(jj,jb,Jph) and AngMom::Triangle(ji,ja,Jph))
////         {
////          for (int J1=J1min;J1<=J1max;++J1)  //Pandya 1: <ai`| V |bj`>_Jtot
////          {
////            V_aibj += AngMom::phase(J1) * (2*J1 + 1) *  modelspace->GetSixJ(ja,ji,Jph,jj,jb,J1)  * H.TwoBody.GetTBME_J(J1,a,b,i,j);
////          }
////         }
////         B(I_ph,II_ph) =  V_aibj *phase_ib *  phase_ai *phase_bj;


//         // I Still don't fully understand these SQRT2 factors...
//         if (m==k)
//         {
//           VB_milk *= SQRT2;
//           VB_mikl *= SQRT2;
//         }
//         if (i==l)
//         {
//           VB_milk *= SQRT2;
//           VB_mikl *= SQRT2;
//         }


         // Ring and Schuck (9.151)
         gamma_mu_kl +=  X(I_mi,mu) * VA_milk  + Y(I_mi,mu) * VB_milk;
         gamma_mu_lk +=  X(I_mi,mu) * VA_mikl  + Y(I_mi,mu) * VB_mikl;

         // Ring and Schuck (9.155)
//         T_mu += OpIn.OneBody(m,i) * X(I_mi,mu) + OpIn.OneBody(i,m)* Y(I_mi,mu);
         T_mu += Omi * X(I_mi,mu) + Oim* Y(I_mi,mu);

//         if ( X(I_mi,mu)== )
//         if ( I_mi == X.col(mu).index_max() )
//         {
//            e_mu = H.OneBody(m,m) - H.OneBody(i,i);
////            std::cout << "mu=" << mu << " -> mi = " << m << " " << i << std::endl;
////            std::cout << " Tmu = " << Omi << " * " << X(I_mi,mu) << "  flipped is " << Oim << std::endl;
////            std::cout << " gamma_mu_kl = " << X(I_mi,mu) << " * " << VA_milk << std::endl;
////            std::cout << " gamma_mu_lk = " << X(I_mi,mu) << " * " << VA_mikl << std::endl;
//         }
         std::cout << "milk = " << m << " " << i << " " << l << " " << k << "   Omi = " << Omi << "  Vmilk = " << VA_milk << "  Vmikl = " << VA_mikl << std::endl;
         I_mi++;
     }// for iket_mi

     double Omega_mu = Energies(mu);
//     Omega_mu = e_mu;

     Teff += gamma_mu_kl * T_mu / ( el-ek-Omega_mu);
     Teff += gamma_mu_lk * T_mu / ( ek-el-Omega_mu);
     std::cout << "    Teff += " << gamma_mu_kl << " * " << T_mu << "  / ( " << el << " - " << ek << " - " << Omega_mu << " ) " << std::endl;
     std::cout << "    Teff += " << gamma_mu_kl << " * " << T_mu << "  / ( " << ek << " - " << el << " - " << Omega_mu << " ) " << std::endl;
     std::cout << "    => Teff = " << Teff << std::endl;
  }// for mu

  return Teff;
}



arma::vec RPA::GetX(size_t i)
{
   return X.col(i);
}

arma::vec RPA::GetY(size_t i)
{
   return Y.col(i);
}

arma::vec RPA::GetEnergies()
{
   return Energies;
}


