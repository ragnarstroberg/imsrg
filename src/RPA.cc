

#include "RPA.hh"
#include "AngMom.hh"


RPA::RPA(ModelSpace& ms)
 : modelspace(&ms) 
{
};

RPA::RPA(Operator& h)
 : modelspace(h.modelspace), H(h) 
{
};



///  In case we want to construct the A matrix for a single channel
///  and it's more convenient to specify J,parity,Tz than the channel index.
void RPA::ConstructAMatrix(int J, int parity, int Tz)
{
   size_t ich_CC = modelspace->GetTwoBodyChannelIndex( J, parity, Tz);
   ConstructAMatrix_byIndex(ich_CC);
}


///  In case we want to construct the B matrix for a single channel
///  and it's more convenient to specify J,parity,Tz than the channel index.
void RPA::ConstructBMatrix(int J, int parity, int Tz)
{
   size_t ich_CC = modelspace->GetTwoBodyChannelIndex( J, parity, Tz);
   ConstructBMatrix_byIndex(ich_CC);
}


///
/// <ai|A|bj> = <ai^-1,J|H|bj^-1,J>
///           = deta_ab delta_ij eps_ai  + <ai^-1,J|V|bj^-1,J>
void RPA::ConstructAMatrix_byIndex(size_t ich_CC)
{
     
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
       double ja = 0.5*modelspace->GetOrbit(a).j2;
       double ji = 0.5*modelspace->GetOrbit(i).j2;

       int phase_ai = 1;
       int phase_ia = - AngMom::phase( ja+ji - Jph );

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
            V_aibj -= modelspace->GetSixJ(ja,ji,Jph,jb,jj,J1)  * (2*J1 + 1) *  H.TwoBody.GetTBME_J_norm(J1,a,j,b,i);
//            V_aibj -= modelspace->GetSixJ(ja,ji,Jph,jb,jj,J1)  * (2*J1 + 1) *  H.TwoBody.GetTBME_J(J1,a,j,b,i);
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
void RPA::ConstructBMatrix_byIndex(size_t ich_CC)
{
     
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




arma::vec RPA::SolveTDA()
{
   arma::vec eigvals;
   arma::mat eigvecs;
   arma::eig_sym(eigvals,eigvecs,A);
//   std::cout << "eigenvalues:"
//	   << std::endl << eigvals << std::endl;
   return eigvals;
}


arma::vec RPA::SolveRPA()
{
   arma::cx_vec eigvals;
   arma::cx_mat eigvecs;

   arma::mat AB = arma::join_vert(  arma::join_horiz( A, B) ,
		                    arma::join_horiz(-B,-A) );
//   std::cout << "AB mat is " << std::endl << AB << std::endl;
   arma::eig_gen(eigvals,eigvecs,AB);
   std::cout << "eigenvalues:"
	   << std::endl << eigvals << std::endl;
   size_t len = eigvals.n_rows;
   double normimag = arma::norm( arma::imag(eigvals),"fro");
   if (normimag>1e-3)
   {
      std::cout << "WARNING: non-zero imaginary part of eigenvalues! " << __FILE__ << " " << __func__ << " " << __LINE__ << std::endl;
   }
   arma::vec good_eigenvals = arma::sort(arma::real( eigvals) );
   good_eigenvals = good_eigenvals.tail_rows(len/2);
   return good_eigenvals;
//   return arma::real(eigvals);
}




