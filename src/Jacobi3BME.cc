
#include "Jacobi3BME.hh"
#include "AngMom.hh"
#include "HartreeFock.hh" // needed for V3monHash and V3monUnHash
#include "IMSRGProfiler.hh"
#include <istream>
#include <set>
#include <sstream>
#include <iomanip>
#include <omp.h>




Jacobi3BME::Jacobi3BME()
{}

Jacobi3BME::Jacobi3BME( int nmax, int twojmin, int twojmax, int twotmin, int twotmax )
 : Nmax(nmax), twoJmin(twojmin), twoJmax(twojmax), twoTmin(twotmin), twoTmax(twotmax)
{
///  Allocate   // Don't call allocate upon construction because we need to read the relevant dimensions from file
  size_t hashmax_tjn = HashTJN(twotmax,twojmax,nmax);
  size_t hashmax_tjnn = HashTJNN(twotmax,twojmax,nmax,nmax);
  dimensionAS.resize(  hashmax_tjn+1, 0);
  dimensionNAS.resize( hashmax_tjn+1, 0);
  cfp_start_loc.resize(hashmax_tjn+1, 0); 
  NAS_jacobi_states.resize(hashmax_tjn+1); 
  start_locAS.resize(  hashmax_tjnn+1, 0);
  start_locNAS.resize( hashmax_tjnn+1, 0);
  std::cout << "Initialized Jacobi3BME, hashmax = " << hashmax_tjn << "  " << hashmax_tjnn << std::endl;
  CalculateJacobiLabels();
}
 

//////  Layout: {   (JTp0) [N0,N0]  [N0,N1]  [N0,N2] ... [N1,N1]  [N1,N2] ... |       (JTp1) ...    |     (JTp2) ... }
//////  Here [N0,N1] is a matrix with dimension dim(T,J,N0) x dim(T,J,N1) 
//////  We store the full symmetrix NN matrix. It's only a factor of 2, and it will avoid a bunch of comparisons later.
//////
void Jacobi3BME::Allocate()
{
//  size_t hashmax = HashTJNN(twoTmax, twoJmax, Nmax, Nmax);
//  start_locAS.resize(  hashmax+1, 0 ); // resize based on the largest hash we can have
//  start_locNAS.resize( hashmax+1 ,0 ); 
  size_t number_mat_el_AS = 0;
  size_t number_mat_el_NAS = 0;
  for (int t=twoTmin; t<=twoTmax; t+=2)
  {
    for (int j=twoJmin; j<=twoJmax; j+=2)
    {
      for (int parity=0; parity<=1; parity++)
      {
        for (int Nbra=parity; Nbra<=Nmax; Nbra+=2)
        {

          size_t dim_braAS = GetDimensionAS(t,j,parity,Nbra);
          size_t dim_braNAS = GetDimensionNAS(t,j,parity,Nbra);

          // set up the list relating non-antisymmetrized basis state index to jacobi1 and jacobi2 indices.
          size_t hashtjn = HashTJN(t,j,Nbra);
          for ( size_t i2=0; i2<jacobi_2.size(); i2++ )
          {
            auto& jac2 = jacobi_2.at(i2);
            for ( size_t i1=0; i1<jacobi_1.size(); i1++ )
            {
              auto& jac1 = jacobi_1.at(i1);
              if ( std::abs(jac1.t*2 - t) !=1 ) continue;  // isospin triangle
              if ( std::abs(jac1.j*2 - j) >jac2.j2 or (jac1.j*2+j)<jac2.j2 ) continue;  // angular momentum triangle
              if ( (2*jac1.n+2*jac2.n + jac1.l + jac2.l) != Nbra ) continue; // conservation of total oscillator quanta
              NAS_jacobi_states.at(hashtjn).push_back({i1,i2});
            }
          }
          if ( NAS_jacobi_states.at(hashtjn).size() != dim_braNAS )
          {
            std::cout << "TROUBLE. I don't understand how the NAS basis states are constructed. Dim = " << NAS_jacobi_states.at(hashtjn).size() << "  should be  " <<  dim_braNAS << std::endl;
          }

          // now allocate   start_locAS,  start_locNAS,  meAS,  meNAS

          for (int Nket=parity; Nket<=Nmax; Nket+=2)
          {
            size_t dim_ketAS = GetDimensionAS(t,j,parity,Nket);
            size_t dim_ketNAS = GetDimensionNAS(t,j,parity,Nket);
            size_t hashtjnn = HashTJNN(t,j,Nbra,Nket);
            start_locAS[ hashtjnn ] = number_mat_el_AS;
            start_locNAS[ hashtjnn ] = number_mat_el_NAS;

            number_mat_el_AS += dim_braAS * dim_ketAS;
            number_mat_el_NAS += dim_braNAS * dim_ketNAS;

          }
        }
      }
    }
  }
  meAS.resize(number_mat_el_AS, 0.);
  meNAS.resize(number_mat_el_NAS, 0.);

  std::cout << "Done with " << __func__ << std::endl;

}



// Perfect hash, assuming 2*T and 2*J can only be odd
// and N%2 = parity 
size_t Jacobi3BME::HashTJN( int twoT, int twoJ, int N)
{
  return   (twoT-twoTmin)/2 * ((twoJmax-twoJmin)/2*(Nmax+1)+Nmax+1 )
         + (twoJ-twoJmin)/2 * (Nmax+1)
         + (N);
}

// Perfect hash, assuming 2*T and 2*J can only be odd
// and N1%2 = N2%2 = parity 
size_t Jacobi3BME::HashTJNN( int twoT, int twoJ, int Nbra, int Nket)
{
  int Nround = (Nmax/2)*2;
////  std::cout << "1) " << (twoT-twoTmin)/2 * ( ((twoJmax-twoJmin)/2*(Nmax/2*(Nround+1))+ Nround+1 ) + Nmax/2*(Nround+1) + Nround + 1 ) << "  "
//  std::cout << "1) " <<  (twoT-twoTmin)/2 * ( (twoJmax-twoJmin)/2 * (Nmax/2*(Nround+1) + Nround+1)  + (Nmax/2)*(Nround+1) + Nround + 1 ) << "  "
//            << "2) " << (twoJ-twoJmin)/2 * (Nmax/2*(Nround+1) + Nround+1) << " "
//            << "3) " << (Nbra/2) * (Nround+1) << " "
//            << "4) " << (Nket/2)*2 + (Nket%2) << "   , Nround =  " << Nround
//            << std::endl;
//  return   (twoT-twoTmin)/2 * ( ((twoJmax-twoJmin)/2*(Nmax/2*(Nround+1))+ Nround+1 ) + Nmax/2*(Nround+1) + Nround + 1 )
  return   (twoT-twoTmin)/2 * ( (twoJmax-twoJmin)/2 * (Nmax/2*(Nround+1) + Nround+1)  + (Nmax/2)*(Nround+1) + Nround + 1  )
         + (twoJ-twoJmin)/2 * (Nmax/2*(Nround+1) + Nround+1)
         + (Nbra/2) * (Nround+1) 
         + (Nket/2)*2 + (Nket%2);   ///  N2%2 gets the parity, beyond that, whether N1 and N2 are even or odd is redundant information
}








// Setter / Getters for the jacobi matrix elements
// Access an antisymmetrized matrix element
// The matrices are stored in column-major order, meaning that sequential elements are (typically)
// from the same column. This means that sequential elements have the same |ket> state, with sequential <bra| states
double Jacobi3BME::GetMatElAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p)
{
//   size_t dimket = GetDimensionAS(twoT,twoJ,p,Nket); 
   size_t dimbra = GetDimensionAS(twoT,twoJ,p,Nbra); 
   size_t start_loc = GetStartLocAS(twoT,twoJ,Nbra,Nket);
//   return meAS.at( start_loc + ibra*dimket + iket );
   return meAS.at( start_loc + iket*dimbra + ibra );
}

// Access an antisymmetrized matrix element
void Jacobi3BME::SetMatElAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p, double matel)
{
//   size_t dimket = GetDimensionAS(twoT,twoJ,p,Nket); 
   size_t dimbra = GetDimensionAS(twoT,twoJ,p,Nbra); 
   size_t start_loc = GetStartLocAS(twoT,twoJ,Nbra,Nket);
//   meAS.at( start_loc + ibra*dimket + iket ) = matel;
   meAS.at( start_loc + iket*dimbra + ibra ) = matel;
}


double Jacobi3BME::GetMatElNAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p)
{
   auto hash = HashTJNN(twoT,twoJ,Nbra,Nket);
//   size_t dimket = GetDimensionNAS(twoT,twoJ,p,Nket); 
   size_t dimbra = GetDimensionNAS(twoT,twoJ,p,Nbra); 
   size_t start_loc = GetStartLocNAS(twoT,twoJ,Nbra,Nket);
//   return meNAS.at( start_loc + ibra*dimket + iket );
   return meNAS.at( start_loc + iket*dimbra + ibra );
}

void Jacobi3BME::SetMatElNAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p, double matel)
{
   auto hash = HashTJNN(twoT,twoJ,Nbra,Nket);
//   size_t dimket = GetDimensionNAS(twoT,twoJ,p,Nket); 
   size_t dimbra = GetDimensionNAS(twoT,twoJ,p,Nbra); 
   size_t start_loc = GetStartLocNAS(twoT,twoJ,Nbra,Nket);
   meNAS.at( start_loc + iket*dimbra + ibra ) = matel;
//   meNAS.at( start_loc + ibra*dimket + iket ) = matel;
}




void Jacobi3BME::GetJacobiStates( int twoT, int twoJ, int parity, int E12, int iNAS, jacobi1_state& jac1, jacobi2_state& jac2)
{
  size_t hash = HashTJN(twoT,twoJ,E12);
  auto& index_1_2 = NAS_jacobi_states.at(hash).at(iNAS);  // for each T,J,p,N block, this maps the NAS index to a pair of jacobi states (actually, to their index)
  jac1 = jacobi_1.at(index_1_2[0]);
  jac2 = jacobi_2.at(index_1_2[1]);

}




void Jacobi3BME::CalculateJacobiLabels()
{
  for (int N=0; N<=Nmax; N++)
  {
    for (int n=0; n<=Nmax/2; n++)
    {
      int l=N-n*2;
      for (int j2=std::abs(2*l-1); j2<=2*l+1; j2+=2)
      {
        jacobi_2.push_back( {n,l,j2} );
      }
    }
  }

  // the jacobi1 coordinate has a totally different ordering, for some reason...
  for (int t=0; t<=1; t++)
  {
    for (int j=0; j<=Nmax+1; j++)
    {
      for (int N=0; N<=Nmax; N++)
      {
        for (int n=0; n<=N/2; n++)
        {
          int l=N-2*n;
          int s = (l+t+1)%2;  // need to have antisymmetric 2-body wf
          if ( j<std::abs(l-s) or j>(l+s) ) continue;
          jacobi_1.push_back( {n,l,s,j,t} );
        }
      }
    }
  }

}







// just need to multiply the AS matrix elements by the CFPs
// CFPs are stored as <AS|cfp|NAS> in row-major order, or
// They're stored in an array as [  <NAS1|cfp|AS1> , <NAS2|cfp|AS1> , <NAS3|cfp|AS1> ... <NAS1|cfp|AS2> ... ]
// As an armadillo matrix, this is interpreted as
// CFP= [ <n1|a1>  <n1|a2> ... ]     so the number of rows is the NAS dimension
//      | <n2|a1>  <n2|a2> ... |
//      [ <n3|a1>  <n3|a2> ... ]
// So to get the NAS matrix, we should do
//  [ n11  n21  n31 ]    [ <n1|a1>  <n1|a2> ...]    [ a11  a21  a31 ]   [ <n1|a1>  <n2|a1>  <n3|a1> ]
//  | n12  n22  n32 | =  | <n2|a1>  <n2|a2> ...|  * | a12  a22  a32 | * | <n1|a2>  <n2|a2>  <n3|a2> |
//  [ n13  n23  n33 ]    [ <n3|a1>  <n3|a2> ...]    [ a13  a23  a33 ]   [ <n1|a3>  <n2|a3>  <n3|a3> ]
//
// NAS = CFP * AS * CFP.T
// <i|NAS|j> = sum_xy <i|CFP|x> <x|AS|y> <y|CFT.t|j>
// In general the flat matrix should be stored such that neighbboring entries typically have the same
// ket state, and sequential bra states.
//
// So I may have this backward...
// conversely, as <NAS|cfp|AS> in column-major order (I think?) so NAS=row,AS=col
// the AS matrix elements are stored as
// <bra|V|ket> in row-major order
// or <ket|V|bra> in column-major order , i.e. ket=row, bra=column
// So to get things in the NAS basis, we need
// <ketNAS|V|braNAS> = Cket <ket|V|bra> Cbra.T  where T is transpose
void Jacobi3BME::ComputeNAS_MatrixElements( )
{
  double t_start = omp_get_wtime();
//  std::cout << std::endl << " Computing NAS matrix elements " << std::endl << std::endl;
  // T2,J2,parity are conserved by V
  size_t num_AS=0;
  size_t num_NAS=0;
  for (int twoT=twoTmin; twoT<=twoTmax; twoT+=2)
  {
    for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
    {
      for (int parity=0; parity<=1; parity++)
      {
//        std::cout << "TJP = " << twoT << " " << twoJ << " " << parity << std::endl;

        for (int Nbra=parity; Nbra<=Nmax; Nbra+=2)  
        {
          size_t dim_braAS = GetDimensionAS(twoT,twoJ,parity,Nbra);
          size_t dim_braNAS = GetDimensionNAS(twoT,twoJ,parity,Nbra);
          size_t cfp_begin_bra = GetCFPStartLocation(twoT,twoJ,Nbra) ;
          if (dim_braAS==0) continue;
          // n_rows = dim_braNAS, n_cols = dim_braAS
          arma::mat cfp_bra( &(cfpvec[cfp_begin_bra]), dim_braNAS, dim_braAS, /*copy_aux_mem*/ true);

          for (int Nket=parity; Nket<=Nmax; Nket+=2)  
          {
//            std::cout << "Nbra,Nket = " << Nbra << " " << Nket << std::endl;
            // build matrix for AS jacobi MEs, CFPs, and  (as the output) NAS jacobi MEs
            size_t dim_ketAS = GetDimensionAS(twoT,twoJ,parity,Nket);
            if (dim_ketAS==0) continue;
         
            size_t dim_ketNAS = GetDimensionNAS(twoT,twoJ,parity,Nket);
            size_t cfp_begin_ket = GetCFPStartLocation(twoT,twoJ,Nket) ;
            arma::mat cfp_ket( &(cfpvec[cfp_begin_ket]), dim_ketNAS, dim_ketAS, /*copy_aux_mem*/ true);
 
            size_t startAS = GetStartLocAS(twoT, twoJ, Nbra, Nket);
            size_t startNAS = GetStartLocNAS(twoT, twoJ, Nbra, Nket);
           // signature is mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem=true, strict=false)
           // and n_rows=dim_bra, n_cols=dim_ket
            arma::mat ASmat( &meAS[startAS], dim_braAS, dim_ketAS, /*copy_aux_mem*/false);

//            arma::mat NASmat = cfp_ket * ASmat * cfp_bra.t();
            arma::mat NASmat = cfp_bra * ASmat * cfp_ket.t();
            num_AS += dim_braAS * dim_ketAS;
            num_NAS += dim_braNAS * dim_ketNAS;

//            std::cout << "ASmat: " << std::endl << ASmat << std::endl << std::endl;
//            std::cout << "NASmat = " << std::endl << NASmat << std::endl << std::endl;

            for (size_t iNAS=0; iNAS<dim_ketNAS*dim_braNAS; iNAS++)
            {
              // single index assumes a flat layout with column-major ordering
              // so we have  [ <bra1|NAS|ket1> , <bra2|NAS|ket1> , <bra3|NAS|ket1> ... <bra1|NAS|ket2> ...]
              meNAS[startNAS+iNAS] = NASmat(iNAS); 
            }
//            if (dim_braNAS>1 and dim_ketNAS>1)
//            {
//              std::cout << "Check: 1,1 element : " << ElementNAS(1,1,Nbra,Nket,T2,J2,parity) << std::endl;;
//            }

//            std::cout << "writing NAS starting at " << startNAS << "   size of meNAS = " << meNAS.size()
//                      << "  first element = " << meNAS[startNAS]
//                      << "  the relevant hash is " << HashTJNN(twoT,twoJ,Nbra,Nket ) << std::endl << std::endl;
          }
        }
      }
    }
  }

  std::cout << "Turned " << num_AS << "  AS elements into  " << num_NAS << "  NAS elements " << std::endl;
  IMSRGProfiler::timer[std::string(__func__)] += omp_get_wtime() - t_start;
}



/*

// Do it the slow stupid way.
double Jacobi3BME::GetLabMatEl( Ket3& bra, Ket3& ket, int Jab, int Jde, int twoJ, int Tab, int Tde, int twoT)
{
  double me_lab = 0;

  int Eabc = 2*(bra.op->n+bra.oq->n+bra.oR->n ) + (bra.op->l+bra.oq->l+bra.oR->l);
  int Edef = 2*(ket.op->n+ket.oq->n+ket.oR->n ) + (ket.op->l+ket.oq->l+ket.oR->l);

  int parity = Eabc%2;
  if ( Edef%2 != parity)
  {
   std::cout << __func__ << " TROUBLE WITH PARITY!!" << std::endl;
  }

  int Ecm_max = std::min( Eabc, Edef);
  for ( int Ecm=0; Ecm<=Ecm_max; Ecm++)
  {
    int E12_bra = Eabc - Ecm;
    int E12_ket = Edef - Ecm;
    int parity_12 = (parity+Ecm) %2;
    for ( int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
    {
      int Ncm = (Ecm-Lcm)/2;
      int twoJ12_min = std::abs(twoJ - 2*Lcm);
      int twoJ12_max = std::min( twoJ+2*Lcm, twoJmax) ;
      for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
      {
        size_t dim_NAS_bra = GetDimensionNAS( twoT, twoJ12, parity_12, E12_bra );
        size_t dim_NAS_ket = GetDimensionNAS( twoT, twoJ12, parity_12, E12_ket );
        for (size_t iNAS_bra=0; iNAS_bra<dim_NAS_bra; iNAS_bra++)
        {
          for (size_t iNAS_ket=0; iNAS_ket<dim_NAS_ket; iNAS_ket++)
          {
            jacobi1_state jac1_bra,jac1_ket;
            jacobi2_state jac2_bra,jac2_ket;
            GetJacobiStates( twoT, twoJ12, parity_12, E12_bra, iNAS_bra, jac1_bra, jac2_bra );
            GetJacobiStates( twoT, twoJ12, parity_12, E12_ket, iNAS_ket, jac1_ket, jac2_ket );
            if ( jac1_bra.t != Tab ) continue;
            if ( jac1_ket.t != Tde ) continue;
            
            double Tbra = Tcoeff_wrapper( bra, Jab, twoJ, jac1_bra, jac2_bra, twoJ12, Ncm, Lcm);
            double Tket = Tcoeff_wrapper( ket, Jde, twoJ, jac1_ket, jac2_ket, twoJ12, Ncm, Lcm);
            double matelNAS = GetMatElNAS( iNAS_bra, iNAS_ket, E12_bra, E12_ket, twoT, twoJ12, parity_12);
            size_t startNAS = GetStartLocNAS(twoT, twoJ, E12_bra, E12_ket);
            me_lab += 6 * Tbra * matelNAS * Tket;

//   size_t dimket = GetDimensionNAS(twoT,twoJ,p,Nket); 
//   size_t start_loc = GetStartLocNAS(twoT,twoJ,Nbra,Nket);
//   return meNAS.at( start_loc + ibra*dimket + iket );

//            std::cout << std::endl << " Ecm, Lcm, Ncm, twoJ12, E12_bra, E12_ket = " << Ecm << " " << Lcm << " " << Ncm << " " << twoJ12 << " " << E12_bra << " " << E12_ket
//                      << "     iNAS_bra,iNAS_ket = " << iNAS_bra << " " << iNAS_ket << std::endl;
//            std::cout << "  (( " << Tbra << "  *  " << matelNAS << "  *  " << Tket << " )) = " << Tbra * Tket * matelNAS << "   =>  " <<  Tbra * Tket * matelNAS*6 << "  sum = "  << me_lab << std::endl;
//            std::cout << "@@@@@ Looking at start location " << startNAS << "  with elements " << iNAS_bra << " , " << iNAS_ket << "  ->  offset = " << iNAS_bra*dim_NAS_ket + iNAS_ket
//                      << "  points to " << meNAS.at( startNAS + iNAS_bra*dim_NAS_ket + iNAS_ket )
//                      << "   first element = " << meNAS.at(startNAS) << std::endl;
          }
        }
      }
    }
  }
  return me_lab;
}

*/



/// Time to do some heavy lifting with T coefficients
double Jacobi3BME::GetLabMatEl( Ket3& bra, Ket3& ket, int Jab, int Jde, int twoJ, int Tab, int Tde, int twoT)
{
  double me_lab = 0;
//  int LabE3max = jacobi_basis.Nmax;
  int LabE3max = 16;

  if (bra.p==bra.q and bra.p==bra.r and bra.op->j2<3 and twoT>1) return 0;
  if (ket.p==ket.q and ket.p==ket.r and ket.op->j2<3 and twoT>1) return 0;

  int Ebra = 2*(bra.op->n+bra.oq->n+bra.oR->n ) + (bra.op->l+bra.oq->l+bra.oR->l);
  int Eket = 2*(ket.op->n+ket.oq->n+ket.oR->n ) + (ket.op->l+ket.oq->l+ket.oR->l);

  int parity_bra = Ebra%2;
  int parity_ket = Eket%2;

//  std::cout << "In GetLabMatEl. Ebra,Eket = " << Ebra << " " << Eket << std::endl;
  for (int Ecm=0; Ecm<=std::min(Ebra,Eket); Ecm++ )
  {
    int E12_bra = Ebra-Ecm; 
    int E12_ket = Eket-Ecm;
    if (E12_bra > Nmax or E12_ket > Nmax ) continue;
    int parity12_bra = E12_bra%2;
    int parity12_ket = E12_ket%2;
    if (parity12_bra != parity12_ket)
    {
      std::cout << "TROUBLE Things aren't looking good in terms of parity...." << std::endl;
      return 0.;
    }
//    for (int Lcm=(Ecm%2); Lcm<=Ecm; Lcm++)
    for (int Lcm=(Ecm%2); Lcm<=Ecm; Lcm+=2)
    {
      int Ncm = (Ecm-Lcm)/2; // Ncm is the number of radial nodes in the CM wave function. Ecm = 2Ncm + Lcm .
//      int twoJ12_min =  std::abs( twoJ - 2*Lcm );
//      int twoJ12_max =  twoJ + 2*Lcm ;
      int twoJ12_min = std::max(twoJmin-Lcm, std::abs( twoJ - 2*Lcm ));
//      int twoJ12_max = std::min(twoJmax+Lcm, twoJ + 2*Lcm) ;
      int twoJ12_max = std::min(twoJmax, twoJ + 2*Lcm) ;
//      int twoJ12_min = std::max(twoJmin, std::abs( twoJ - 2*Lcm ));
//      int twoJ12_max = std::min(twoJmax, twoJ + 2*Lcm) ;

//      std::cout << "Ecm,Lcm,Ncm = " << Ecm << " " << Lcm << " " << Ncm << std::endl;
      for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2) // the interaction conserves J12, parity, and T12=twoT
      {
//        std::cout << "Getting dimensions for J12 = " << twoJ12 << "/2" << "   twoJmax = " << twoJmax << std::endl;
        int NASdim_bra = GetDimensionNAS( twoT, twoJ12, parity12_bra, E12_bra ); 
        int NASdim_ket = GetDimensionNAS( twoT, twoJ12, parity12_ket, E12_ket ); 

        size_t ASdim_bra = GetDimensionAS( twoT, twoJ12, parity12_bra, E12_bra ); 
        size_t ASdim_ket = GetDimensionAS( twoT, twoJ12, parity12_ket, E12_ket ); 
        if (ASdim_bra==0 or ASdim_ket==0) continue;

        size_t cfp_begin_bra = GetCFPStartLocation(twoT,twoJ12,E12_bra);
        size_t cfp_begin_ket = GetCFPStartLocation(twoT,twoJ12,E12_ket);
        size_t startlocAS = GetStartLocAS( twoT, twoJ12, E12_bra, E12_ket);

        arma::mat matelAS( &meAS[startlocAS], ASdim_bra, ASdim_ket, false ); 
        
        arma::mat cfp_bra( &(cfpvec[cfp_begin_bra]), NASdim_bra, ASdim_bra, /*copy_aux_mem*/ true);
        arma::mat cfp_ket( &(cfpvec[cfp_begin_ket]), NASdim_ket, ASdim_ket, /*copy_aux_mem*/ true);


//        std::cout << " ^^^TJP E12 = " << twoT << " " << twoJ12 << " " << parity12_bra << " " << E12_bra << "," << E12_ket << "  Ncm,Lcm = " << Ncm << " " << Lcm << std::endl;
//        std::cout << "dimensions: " << NASdim_bra << " " << NASdim_ket << std::endl;
//           size_t dimbra = jacobi_basis.GetDimensionNAS(T,J,parity12_bra,E12_bra); 
//           std::cout << "dimbra = " <<  dimbra << std::endl;

//  compute the Tcoefficients that we'll need here
//        std::vector<double> Tcoeff_bra(NASdim_bra);
//        std::vector<double> Tcoeff_ket(NASdim_ket);
        arma::rowvec Tcoeff_bra(NASdim_bra, arma::fill::zeros);
        arma::vec Tcoeff_ket(NASdim_ket, arma::fill::zeros);
//        std::cout << "start loop over ibraNAS. NASdim_bra = " << NASdim_bra << std::endl;
        for (int ibraNAS=0; ibraNAS<NASdim_bra; ibraNAS++)
        {
          jacobi1_state jac1_bra;
          jacobi2_state jac2_bra;
//          std::cout << "GetJacobi bra state" << std::endl;
          GetJacobiStates( twoT, twoJ12, parity12_bra, E12_bra, ibraNAS, jac1_bra, jac2_bra);
          if (jac1_bra.t != Tab ) continue;
//          std::cout << " ibraNAS = " << ibraNAS << ": | " << jac1_bra.n << " " << jac1_bra.l << " " << jac1_bra.s << " " << jac1_bra.j << " " << jac1_bra.t << " > x | "
//                                     << jac2_bra.n << " " << jac2_bra.l << " " << jac2_bra.j2 << " > " << std::endl;
//          Tcoeff_bra.at(ibraNAS) = Tcoeff_wrapper( bra, Jab, twoJ, jac1_bra, jac2_bra, twoJ12, Ncm, Lcm) ;
          Tcoeff_bra[ibraNAS] = Tcoeff_wrapper( bra, Jab, twoJ, jac1_bra, jac2_bra, twoJ12, Ncm, Lcm) ;
        }
//        std::cout << "start loop over iketNAS. NASdim_ket = " << NASdim_ket << std::endl;
        for (int iketNAS=0; iketNAS<NASdim_ket; iketNAS++)
        {
          jacobi1_state jac1_ket;
          jacobi2_state jac2_ket;
//          std::cout << "GetJacobi ket state" << std::endl;
          GetJacobiStates( twoT, twoJ12, parity12_ket, E12_ket, iketNAS, jac1_ket, jac2_ket);
          if (jac1_ket.t != Tde ) continue;
//          std::cout << " iketNAS = " << iketNAS << ": | " << jac1_ket.n << " " << jac1_ket.l << " " << jac1_ket.s << " " << jac1_ket.j << " " << jac1_ket.t << " > x | "
//                                     << jac2_ket.n << " " << jac2_ket.l << " " << jac2_ket.j2 << " > " << std::endl;
//          Tcoeff_ket.at(iketNAS) = Tcoeff_wrapper( ket, Jde, twoJ, jac1_ket, jac2_ket, twoJ12, Ncm, Lcm) ;
          Tcoeff_ket[iketNAS] = Tcoeff_wrapper( ket, Jde, twoJ, jac1_ket, jac2_ket, twoJ12, Ncm, Lcm) ;
        }



        arma::mat result =  Tcoeff_bra * cfp_bra * matelAS * cfp_ket.t() * Tcoeff_ket;
        me_lab += result[0];

//        for (int ibraNAS=0; ibraNAS<NASdim_bra; ibraNAS++)
//        {
//          if (std::abs(Tcoeff_bra[ibraNAS])<1e-9) continue;
//          jacobi1_state jac1_bra,jac1_ket;
//          jacobi2_state jac2_bra,jac2_ket;
//
////          std::cout << "about to get jac1_bra and jac2_bra" << std::endl;
//          GetJacobiStates( twoT, twoJ12, parity12_bra, E12_bra, ibraNAS, jac1_bra, jac2_bra);
////          std::cout << "checking isospin bra: " << jac1_bra.t << " " << Tab << std::endl;
//          if ( jac1_bra.t != Tab ) continue;
////          double Tcoeff_bra = Tcoeff_wrapper( bra, Jab, twoJ, jac1_bra, jac2_bra, twoJ12, Ncm, Lcm);
////          if (std::abs(Tcoeff_bra)<1e-9) continue;
//
//          for (int iketNAS=0; iketNAS<NASdim_ket; iketNAS++)
//          {
//           if (std::abs(Tcoeff_ket[iketNAS])<1e-9) continue;
////           std::cout << "  ibraNAS,iketNAS = " << ibraNAS << " " << iketNAS << std::endl;
//           GetJacobiStates( twoT, twoJ12, parity12_ket, E12_ket, iketNAS, jac1_ket, jac2_ket);
////           std::cout << "     checking isospin ket: " << jac1_ket.t << " " << Tde << std::endl;
//           if ( jac1_ket.t != Tde ) continue;
////           double Tcoeff_ket = Tcoeff_wrapper( ket, Jde, twoJ, jac1_ket, jac2_ket, twoJ12, Ncm, Lcm);
//
////           std::cout << "    ** finding the NAS matrix element. T = " << T << "   J = " << J12 << "  p = " << parity12_bra << std::endl;
//
////           double meNAS = ElementNAS( ibraNAS, iketNAS, E12_bra, E12_ket, twoT, twoJ12, parity12_bra );
//           double meNAS = GetMatElNAS( ibraNAS, iketNAS, E12_bra, E12_ket, twoT, twoJ12, parity12_bra );
//
////           me_lab +=  Tcoeff_bra * meNAS * Tcoeff_ket  ;
//           me_lab +=  Tcoeff_bra[ibraNAS] * meNAS * Tcoeff_ket[iketNAS]  ;
////           std::cout << "Ecm,Lcm,twoJ12 = " << Ecm << " " << Lcm << " " << twoJ12 << "  Tab,Tde,twoT = " << Tab << " " << Tde << " " << twoT << std::endl;
////           std::cout << "jac1_bra " << jac1_bra.n << " " << jac1_bra.l << " " << jac1_bra.s << " " << jac1_bra.j << " " << jac1_bra.t << "  "
////                     << "jac2_bra " << jac2_bra.n << " " << jac2_bra.l << " " << jac2_bra.j2 <<  "  "
////                     << "jac1_ket " << jac1_ket.n << " " << jac1_ket.l << " " << jac1_ket.s << " " << jac1_ket.j << " " << jac1_ket.t << "  "
////                     << "jac2_bra " << jac2_bra.n << " " << jac2_bra.l << " " << jac2_bra.j2  << std::endl;
////           std::cout << " " << ibraNAS << " , " << iketNAS << ":  (( " << Tcoeff_bra[ibraNAS] << "  " << meNAS << "  " << Tcoeff_ket[iketNAS] << " ))  "
////                     << Tcoeff_bra[ibraNAS] * meNAS * Tcoeff_ket[iketNAS] <<  "  => sum_ME = " << me_lab << std::endl;
////           std::cout << "    < N1=" << jac1_bra.n << " L1=" << jac1_bra.l << " S1=" << jac1_bra.s << " J1=" << jac1_bra.j << " T1=" <<jac1_bra.t << ",  N2=" << jac2_bra.n << " L2=" << jac2_bra.l << " T=" << twoT << " J12=" << twoJ12 << " | ..." << std::endl;
//
//          }// for ibraNAS
//        } // for iketNAS
      } // for J12
    } // for Lcm
  } // for Ecm


  return 6 * me_lab;
}






/*
// It's not clear that this will ever be useful.
double Jacobi3BME::GetV3mon( size_t a, size_t b, size_t c, size_t d, size_t e, size_t f )
{


}

*/

// na,nb,and nc all should be <= emax/2, Jab will be <= 2*emax+1, twoJ will be <=6*emax+3, but only odd, Lcm will be <=3*emax, jac2 <= ~Nmax*Nmax,    jac1<= ~Nmax*Nmax*Nmax 
//  J lies in the range |Jab-jc| <= J <= Jab+jc,  so it has a span of no more than 2*emax+1.
// na + 10*nb + 100*nc fits in 10 bits (1023).  Jab fits in 6 bits (63).  twoJ/2 fits in 7 bits (127).  Lcm fits in 7 bits (127). jac2 fits in 12 bits (4095). jac1 fits in 18 bits (262143)
// This totals up to 10+6+7+7+7+12+18 = 67, which is too many.   For now, we will just use a string.
//uint64_t Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ, uint64_t jac1, uint64_t jac2, uint64_t twoJ12, uint64_t Lcm )
//std::string Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ, uint64_t jac1, uint64_t jac2, uint64_t twoJ12, uint64_t Lcm )
std::string Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ,  uint64_t twoJ12, uint64_t E12 )
{
  std::ostringstream oss;
//  oss << na << " " << nb << " " << nc << " " << Jab << " " << twoJ << " " << jac1 << " " << jac2 << " " << twoJ12 << " " << Lcm;
  oss << na << " " << nb << " " << nc << " " << Jab << " " << twoJ << " " <<  " " << twoJ12 << " " << E12;
  return oss.str();
}

//void Jacobi3BME::TcoeffUnHash(uint64_t key, uint64_t& na, uint64_t& nb, uint64_t& nc, uint64_t& Jab, uint64_t& twoJ, uint64_t& jac1, uint64_t& jac2, uint64_t& twoJ12, uint64_t& Lcm )
//void Jacobi3BME::TcoeffUnHash(std::string& key, int& na, int& nb, int& nc, int& Jab, int& twoJ, int& jac1, int& jac2, int& twoJ12, int& Lcm )
void Jacobi3BME::TcoeffUnHash(std::string& key, int& na, int& nb, int& nc, int& Jab, int& twoJ,  int& twoJ12, int& E12 )
{
  std::istringstream iss(key);
//  iss >> na >> nb >> nc >> Jab >> twoJ >> jac1 >> jac2 >> twoJ12 >> Lcm;
  iss >> na >> nb >> nc >> Jab >> twoJ >> twoJ12 >> E12;
}



// Find all the entries in HartreeFock::Vmon3 that we need to compute for a specific set of one-body channels (la,j2a, ...)
// We want to get all the free ones we can from symmetry (aside from swapping the third index, e.g. abc -> cab. I haven't mustered the gumption for that yet).
// This is a bit awkward because the Vmon3 and Vmon3_keys are just vectors of doubles and uint64_t, respectively.
// The strategy is to loop through Vmon_keys, and pick out all the ones that belong to this one-body channel, or to the one with ab->ba, de->ed.
// We store these in the std::set "need_to_compute". We then loop through those and figure out which ones are related by symmetry.
// We put a set of 8 indices into the structure "indices", where we explicitly compute entry [0], and get the others for free.
// We then add all 8 to the std::set "will_be_computed" and continue looping. One the next iteration, if the entry from "need_to_compute" is already in "will_be_computed"
// we skip it. Since a std::set is a container of unique entries, we avoid computing the same thing twice.
void Jacobi3BME::GetMonopoleIndices( int la, int j2a, int lb, int j2b, int lc, int j2c, HartreeFock& hf, std::vector<std::array<size_t,8>>& indices )
{
  std::set<size_t> need_to_compute;
  std::set<size_t> will_be_computed;

  bool verbose = false;

  if (verbose) std::cout << "---------------------------------------------------" << std::endl;
  if (verbose) std::cout << "lj channel : " << la << " " << j2a << " " << lb << " " << j2b << " " << lc << " " << j2c << std::endl;


  for (size_t imon=0; imon<hf.Vmon3_keys.size(); imon++)
  {
     auto key = hf.Vmon3_keys[imon];
     int a,b,c,d,e,f;
     hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we could call it without a class instance if we wanted...
     Orbit& oa = hf.modelspace->GetOrbit(a);
     Orbit& ob = hf.modelspace->GetOrbit(b);
     Orbit& oc = hf.modelspace->GetOrbit(c);
     if ( oc.l != lc or oc.j2 != j2c ) continue;
     if ( (oa.l==la and oa.j2==j2a and ob.l==lb and ob.j2==j2b) or (oa.l==lb and oa.j2==j2b and ob.l==la and ob.j2==j2a) )    need_to_compute.insert(imon);
  }

  for ( size_t imon : need_to_compute )
  {
     if ( will_be_computed.find( imon ) != will_be_computed.end() ) continue; // we've already got this one...
     auto key = hf.Vmon3_keys[imon];
     int a,b,c,d,e,f;
     hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we could call it without a class instance if we wanted...
     Orbit& oa = hf.modelspace->GetOrbit(a);
     Orbit& ob = hf.modelspace->GetOrbit(b);
     Orbit& oc = hf.modelspace->GetOrbit(c);
     Orbit& od = hf.modelspace->GetOrbit(d);
     Orbit& oe = hf.modelspace->GetOrbit(e);
     Orbit& of = hf.modelspace->GetOrbit(f);
     if ( not (oa.l==la and oa.j2==j2a and ob.l==lb and ob.j2==j2b) ) continue; // we'd like our [0] entry to be the direct term, not the flipped one.
     int aa = hf.modelspace->GetOrbitIndex( oa.n, oa.l, oa.j2, -oa.tz2 );
     int bb = hf.modelspace->GetOrbitIndex( ob.n, ob.l, ob.j2, -ob.tz2 );
     int cc = hf.modelspace->GetOrbitIndex( oc.n, oc.l, oc.j2, -oc.tz2 );
     int dd = hf.modelspace->GetOrbitIndex( od.n, od.l, od.j2, -od.tz2 );
     int ee = hf.modelspace->GetOrbitIndex( oe.n, oe.l, oe.j2, -oe.tz2 );
     int ff = hf.modelspace->GetOrbitIndex( of.n, of.l, of.j2, -of.tz2 );
     auto key_bacedf = hf.Vmon3Hash(b,a,c,e,d,f);
     auto key_defabc = hf.Vmon3Hash(d,e,f,a,b,c);
     auto key_edfbac = hf.Vmon3Hash(e,d,f,b,a,c);
     auto key_isoflip = hf.Vmon3Hash(aa,bb,cc,dd,ee,ff);
     auto key_bacedf_isoflip = hf.Vmon3Hash(bb,aa,cc,ee,dd,ff);
     auto key_defabc_isoflip = hf.Vmon3Hash(dd,ee,ff,aa,bb,cc);
     auto key_edfbac_isoflip = hf.Vmon3Hash(ee,dd,ff,bb,aa,cc);
     std::array<size_t,8> compute_these;
     for (int i=0;i<8;i++) compute_these[i] = imon;
     for ( auto imon1 : need_to_compute )
     {
       if ( hf.Vmon3_keys[imon1] == key_isoflip)         compute_these[1] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_defabc)          compute_these[2] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_defabc_isoflip)  compute_these[3] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_bacedf)          compute_these[4] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_edfbac)          compute_these[5] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_bacedf_isoflip)  compute_these[6] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_edfbac_isoflip)  compute_these[7] = imon1;
     }
     for (int i=0;i<8;i++) will_be_computed.insert( compute_these[i] );
     indices.push_back( compute_these );
  }

  if (verbose)
  {
    std::cout << "need_to_compute" << std::endl;
    for (auto i : need_to_compute ) std::cout << i << " ";
    std::cout << std::endl << std::endl << "will_compute" << std::endl;
    for (auto i : will_be_computed ) std::cout << i << " ";
    std::cout << std::endl << std::endl << "indices" << std::endl;
    for (auto& ilist : indices)
    {
     std::cout << "( ";
      for (int i=0;i<8;i++) std::cout << ilist[i] << " ";
     std::cout << "  ) " << std::endl;
    }
  }

}



/*

// Go through the Vmon3_keys from HartreeFock and find all the ones that live in the channel defined by la,j2a,lb,j2b,lc,j2c
// Then, because several of them are related by symmetry, we put 8 indices into each entry, each of which should have the same
// value for the 3-body monopole. That way, we compute one of them and set all 8.
void Jacobi3BME::GetMonopoleIndices( int la, int j2a, int lb, int j2b, int lc, int j2c, HartreeFock& hf, std::vector<std::array<size_t,8>>& indices )
{
   std::vector<size_t> imon_dir;
   std::vector<size_t> imon_exch;
   std::vector<size_t> imonlist;
   for (size_t imon=0; imon<hf.Vmon3_keys.size(); imon++)
   {
     // Check if the given key corresponds to the set of one-body channels we're currently working on.
     auto key = hf.Vmon3_keys[imon];
     int a,b,c,d,e,f;
     hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we could call it without a class instance if we wanted...
     Orbit& oa = hf.modelspace->GetOrbit(a);
     Orbit& ob = hf.modelspace->GetOrbit(b);
     Orbit& oc = hf.modelspace->GetOrbit(c);

     Orbit& od = hf.modelspace->GetOrbit(d);
     Orbit& oe = hf.modelspace->GetOrbit(e);
     Orbit& of = hf.modelspace->GetOrbit(f);
     if ( oc.l != lc or oc.j2 != j2c ) continue;

//     if ( oa.l==la and oa.j2==j2a and ob.l==lb and ob.j2==j2b)
//     {
//       imon_dir.push_back(imon);
//       if (a<d or (a==d and b<e) or (a==d and b==e and c<=f) ) // Choose a specific ordering which we will compute. Then get the others by symmetry.
//       {
////         if ( a<=b and d<=e and (oa.tz2+ob.tz2+oc.tz2)<0 )
//         if ( (oa.tz2+ob.tz2+oc.tz2)<0 )
//           imonlist.push_back(imon);
//       }
//     }
//     else if ( oa.l==lb and oa.j2==j2b and ob.l==la and ob.j2==j2a)
//     {
////       imon_exch.push_back(imon);
//       imon_dir.push_back(imon);
//     }
//   }

     if ( oa.l==la and oa.j2==j2a and ob.l==lb and ob.j2==j2b)
     {

       imon_dir.push_back(imon);

       if ( (oa.tz2+ob.tz2+oc.tz2)<0 )
       {
         if (la==lb and j2a==j2b and oa.n > ob.n ) continue;
//         if (oa.tz2 > ob.tz2 ) continue;
//         if ((oa.n > od.n)  ) continue;
//         if ((oa.n > od.n) or (oa.n==od.n and (ob.n+oc.n)>(oe.n+of.n)) ) continue;
//         if (a<d or (a==d and b<e) or (a==d and b==e and c<=f) )  imonlist.push_back(imon);
         imonlist.push_back(imon);
       }
//       if (a<d or (a==d and b<e) or (a==d and b==e and c<=f) ) // Choose a specific ordering which we will compute. Then get the others by symmetry.
//       {
////         if ( a<=b and d<=e and (oa.tz2+ob.tz2+oc.tz2)<0 )
//         if ( (oa.tz2+ob.tz2+oc.tz2)<0 )
//           imonlist.push_back(imon);
//       }
     }
     else if ( oa.l==lb and oa.j2==j2b and ob.l==la and ob.j2==j2a)
     {
//       imon_exch.push_back(imon);
       imon_dir.push_back(imon);
//       imonlist.push_back(imon);
     }
   }

   indices.resize( imonlist.size() );

   for ( size_t index=0; index<imonlist.size(); index++)
   {
     size_t imon = imonlist[index];
     for (int i=0;i<8;i++) indices[index][i] = imon;

     auto key = hf.Vmon3_keys[imon];
     int a,b,c,d,e,f;
     hf.Vmon3UnHash( key, a, b, c, d, e, f);  
     Orbit& oa = hf.modelspace->GetOrbit(a);
     Orbit& ob = hf.modelspace->GetOrbit(b);
     Orbit& oc = hf.modelspace->GetOrbit(c);
     Orbit& od = hf.modelspace->GetOrbit(d);
     Orbit& oe = hf.modelspace->GetOrbit(e);
     Orbit& of = hf.modelspace->GetOrbit(f);
     int aa = hf.modelspace->GetOrbitIndex( oa.n, oa.l, oa.j2, -oa.tz2 );
     int bb = hf.modelspace->GetOrbitIndex( ob.n, ob.l, ob.j2, -ob.tz2 );
     int cc = hf.modelspace->GetOrbitIndex( oc.n, oc.l, oc.j2, -oc.tz2 );
     int dd = hf.modelspace->GetOrbitIndex( od.n, od.l, od.j2, -od.tz2 );
     int ee = hf.modelspace->GetOrbitIndex( oe.n, oe.l, oe.j2, -oe.tz2 );
     int ff = hf.modelspace->GetOrbitIndex( of.n, of.l, of.j2, -of.tz2 );
     auto key_bacedf = hf.Vmon3Hash(b,a,c,e,d,f);
     auto key_defabc = hf.Vmon3Hash(d,e,f,a,b,c);
     auto key_edfbac = hf.Vmon3Hash(e,d,f,b,a,c);
     auto key_isoflip = hf.Vmon3Hash(aa,bb,cc,dd,ee,ff);
     auto key_bacedf_isoflip = hf.Vmon3Hash(bb,aa,cc,ee,dd,ff);
     auto key_defabc_isoflip = hf.Vmon3Hash(dd,ee,ff,aa,bb,cc);
     auto key_edfbac_isoflip = hf.Vmon3Hash(ee,dd,ff,bb,aa,cc);
     for ( auto imon1 : imon_dir )
     {
       if ( hf.Vmon3_keys[imon1] == key_isoflip) indices[index][1] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_defabc)  indices[index][2] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_defabc_isoflip)  indices[index][3] = imon1;
////       if ( hf.Vmon3_keys[imon1] == key_defabc_isoflip)  indices[index][3] = imon1;
//       if ( hf.Vmon3_keys[imon1] == key_defabc)  indices[index][3] = imon1;
//     }
////     for (auto imon2 : imon_exch )
//     for (auto imon2 : imon_dir )
//     {
//       if ( hf.Vmon3_keys[imon2] == key_bacedf) indices[index][4] = imon2;
//       if ( hf.Vmon3_keys[imon2] == key_edfbac) indices[index][5] = imon2;
//       if ( hf.Vmon3_keys[imon2] == key_bacedf_isoflip) indices[index][6] = imon2;
//       if ( hf.Vmon3_keys[imon2] == key_edfbac_isoflip) indices[index][7] = imon2;
       if ( hf.Vmon3_keys[imon1] == key_bacedf) indices[index][4] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_edfbac) indices[index][5] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_bacedf_isoflip) indices[index][6] = imon1;
       if ( hf.Vmon3_keys[imon1] == key_edfbac_isoflip) indices[index][7] = imon1;
     }
     
   }
}
*/


void Jacobi3BME::GetRelevantTcoeffs( int la, int j2a, int lb, int j2b, int lc, int j2c, HartreeFock& hf)   
{
    double t_internal = omp_get_wtime();
    bool verbose = false;
//    std::cout << " IN " << __func__ << std::endl;
//    if (la==0 and lc==0 and lb==1 and j2b==3) verbose = true;
//    if ((la==0 and lb==1 and lc==0 and j2b==3) or (la==1 and lb==0 and lc==0 and j2a==3)) verbose = true;
//    std::cout << "la, ..." << la << " " << j2a << " " << lb << " " << j2b << " " << lc << " " << j2c <<  "     " << verbose << std::endl;
    std::set<int> na_list, nb_list,nc_list; // a std::set is a sorted unique list of items
    for (auto& orb : hf.modelspace->Orbits )
    {
      if (orb.l==la and orb.j2==j2a) na_list.insert(orb.n); 
      if (orb.l==lb and orb.j2==j2b) nb_list.insert(orb.n); 
      if (orb.l==lc and orb.j2==j2c) nc_list.insert(orb.n); 
    }

    TcoeffLookup.clear();
    int tcoeff_count = 0;

    int Jab_min = std::abs(j2a-j2b)/2;
    int Jab_max = (j2a+j2b)/2;

    if (verbose) std::cout <<"Jab_min,Jab_max = " << Jab_min << " " << Jab_max << std::endl;

     for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
     {
       if (verbose) std::cout << " Jab = " << Jab << std::endl;
       int twoJ_min = std::abs(2*Jab-j2c);
       int twoJ_max = 2*Jab+j2c;
       for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
       {
         for ( auto na : na_list)
         {
           auto a = hf.modelspace->GetOrbitIndex( na, la, j2a, 1);
          for ( auto nb : nb_list )
          {
           auto b = hf.modelspace->GetOrbitIndex( nb, lb, j2b, 1);
//           if (a>b) continue;
           for (auto nc : nc_list )
           {
//    if (a==1 and b==2 and lc==0 and nc==0) verbose = true;
//    else verbose = false;
            int Eabc = 2*(na+nb+nc) + la+lb+lc;

            for (int Ecm=0; Ecm<=Eabc; Ecm++)
            {
             int E12 = Eabc-Ecm;
             if (E12 > Nmax) continue;
             int parity = E12%2;

             int twoJ12_min = std::max(1,twoJ-2*Ecm);
             int twoJ12_max = std::min( twoJ+2*Ecm, twoJmax);
             for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
             {
              size_t dimNAS_T1 = GetDimensionNAS(1,twoJ12, parity, E12 );
              size_t dimNAS_T3 = GetDimensionNAS(3,twoJ12, parity, E12 );
              // if we want, we could also check the AS dimension (make sure they're not both zero)
              auto hash_key = TcoeffHash(na,nb,nc,Jab,twoJ,twoJ12,E12);
              TcoeffLookup[ hash_key ] = tcoeff_count;
              for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
              {
                tcoeff_count += 2*dimNAS_T1 + dimNAS_T3 ; //  [T=1/2,Tab=0], [T=1/2,Tab=1],  [T=3/2,Tab=1]
                if (verbose) std::cout << "na,nb,nc  " <<  na << " " << nb << " " << nc << "   Jab,twoJ,twoJ12,E12,Lcm " << Jab << " " << twoJ  << " " << twoJ12 << " " << E12 <<  " " << Lcm << "   incremented to " << tcoeff_count << std::endl;
              }// for Lcm
             } // for twoJ12
            }// for Ecm
           } // for nc
          } // for nb
         } // for na
       } // for twoJ
     } // for Jab

     IMSRGProfiler::timer[std::string(__func__)+"_GenerateKeys"] += omp_get_wtime() - t_internal;
     t_internal = omp_get_wtime();

    if (verbose) std::cout << "Resizing TcoeffList to " << tcoeff_count << std::endl;
    TcoeffList.assign(tcoeff_count,0.0);
//    std::fill(TcoeffList.begin(),TcoeffList.end(),0.0);
    if (verbose)
    {
    std::cout << "size is " << TcoeffList.size() << std::endl;
    for (size_t i=0;i<TcoeffList.size(); i++) std::cout << TcoeffList[i] << " ";
    std::cout << std::endl;
    }
     int tcomputed = 0;
//     #pragma omp parallel  reduction(+ : tzero, time_zero, time_nonzero )
     #pragma omp parallel  reduction(+ : tcomputed)
     {
       size_t cnt = 0;
       int ithread = omp_get_thread_num();
       int nthreads = omp_get_num_threads();
       for(auto element = TcoeffLookup.begin(); element !=TcoeffLookup.end(); ++element, cnt++)
       {
         if (verbose) std::cout << "ithread, cnt, nthreads = " << ithread << " " << cnt << " " << nthreads << std::endl;
         if(cnt%nthreads != ithread) continue; // Check if this thread should compute this element
         auto hash_key = element->first;
         double localtime = omp_get_wtime();
         int na,nb,nc,Jab,twoJ,twoJ12,E12,Lcm;
         TcoeffUnHash(hash_key, na,nb,nc,Jab,twoJ,twoJ12,E12);
         bool ab_same = (la==lb and j2a==j2b and na==nb);

         size_t start_point = element->second;
         size_t offset = 0;
         if (verbose) std::cout << " ~ na,nb,nc " << na << " " << nb << " " << nc << "     Jab,twoJ,twoJ12: " << Jab << " " << twoJ << " " << twoJ12 << "   E12: " << E12 << std::endl;
         if (verbose) std::cout << "hash: " << element->first << "  ->  " << element->second << std::endl;
 

         int Eabc = 2*(na+nb+nc)+la+lb+lc;
         int Ecm = Eabc-E12;
         int parity = E12%2;

         size_t dimNAS_T1 = GetDimensionNAS(1,twoJ12, parity, E12 );
         size_t dimNAS_T3 = GetDimensionNAS(3,twoJ12, parity, E12 );

         if (verbose)
         {
           std::cout << "Workspace, starting at " << start_point << ",  size = " << (Ecm/2+1)*(2*dimNAS_T1+dimNAS_T3) << "  : " << std::endl;
           for (int i=0;i<(Ecm/2+1)*(2*dimNAS_T1 + dimNAS_T3 ); i++) std::cout << TcoeffList[start_point +i] << " ";
           std::cout << std::endl;
         }


         for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
         {
           if ( std::abs(twoJ-twoJ12)>2*Lcm or (twoJ+twoJ12)<2*Lcm) continue;
           int Ncm = (Ecm-Lcm)/2;
           offset = (Lcm-Ecm%2)/2 * (2*dimNAS_T1 + dimNAS_T3 );

           if (verbose)
           {
              std::cout << "Lcm = " << Lcm << "    workspace for this block is " << std::endl;
              for (int i=0; i<(2*dimNAS_T1+dimNAS_T3); i++) std::cout << TcoeffList[ start_point + offset + i] << " ";
              std::cout << std::endl;
           }
           // Layout is [ T=1/2,Tab=0] [T=1/2,Tab=1] [T=3/2,Tab=1] ... (repeat for each Lcm)
           //              dimNAS_T1     dimNAS_T1     dimNAS_T3
           for (size_t iNAS=0; iNAS<dimNAS_T1; iNAS++)
           {
             jacobi1_state jac1;
             jacobi2_state jac2;
             GetJacobiStates( 1, twoJ12, parity, E12, iNAS, jac1, jac2);
             if (ab_same and  (jac1.t+Jab)%2==0) continue;
//             double tcoef = ComputeTcoeff(hf, na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
             double tcoef = ComputeTcoeff( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
             size_t location = start_point + offset + iNAS + jac1.t*dimNAS_T1;
             if (verbose) std::cout << "T=1/2: Writing tcoefficent at location " << start_point << " + " << offset << " + " << iNAS << "  + " << jac1.t*dimNAS_T1 << " = " << location
                       << "    and setting it to " << tcoef << "    size of TcoeffList = " << TcoeffList.size() << std::endl;
//             TcoeffList[ start_point + offset + iNAS + jac1.t*dimNAS_T1 ] = tcoef; // add an offset for Tab=1.
             TcoeffList[ location ] = tcoef; // add an offset for Tab=1.
           }
           // Now the T=3/2 part
           for (size_t iNAS=0; iNAS<dimNAS_T3; iNAS++)
           {
             jacobi1_state jac1;
             jacobi2_state jac2;
             GetJacobiStates( 3, twoJ12, parity, E12, iNAS, jac1, jac2);
             if (ab_same and  (jac1.t+Jab)%2==0) continue;
//             double tcoef = ComputeTcoeff(hf, na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
             double tcoef = ComputeTcoeff( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
             if (verbose) std::cout << "T=3/2: Writing tcoefficent at location " << start_point << " + " << offset << " + " << 2*dimNAS_T1 << " + " << iNAS  << " = " << start_point + offset + 2*dimNAS_T1 + iNAS 
                       << "    and setting it to " << tcoef << std::endl;
             TcoeffList[ start_point + offset + 2*dimNAS_T1 + iNAS ] = tcoef; 
             tcomputed++;
//             offset++;
           }
         }
       }
     } 
     IMSRGProfiler::timer[std::string(__func__)+"_ComputeToeff"] += omp_get_wtime() - t_internal;
     IMSRGProfiler::counter["Tcoefficients"] += tcoeff_count;
     IMSRGProfiler::counter["Tcoeff computed"] += tcomputed;

}





//void Jacobi3BME::GetV3mon_all( std::vector<uint64_t>& keys, std::vector<double>& v3mon, ModelSpace& modelspace )
void Jacobi3BME::GetV3mon_all( HartreeFock& hf )
{
  double t_start = omp_get_wtime();
  // The keys are already computed elsewhere and are passed in as input
//  std::cout << "Begin " << __func__ << std::endl;
//  hf.modelspace->PreCalculateMoshinsky();
  PreComputeMoshinsky1();
  PreComputeMoshinsky2();
  PreComputeSixJ();
  PreComputeNineJ();
//  hf.modelspace->PreCalculateAdditionalSixJ();
  hf.Vmon3.resize( hf.Vmon3_keys.size(), 0.);
  struct ljChannel{
     int l; int j2;
     bool operator == (const ljChannel& rhs){return (rhs.l==l and rhs.j2==j2); };
     bool operator <  (const ljChannel& rhs) const {return (l<rhs.l or (l==rhs.l and j2>rhs.j2) ) ;};
  };

//  std::vector<ljChannel> ljchannels;
  std::set<ljChannel> ljchannels;

//  auto& onebodychan = hf.modelspace->OneBodyChannels;
//  typedef obc_type std::map<std::array<int,3>,std::vector<index_t> >;
//  for ( auto& obc : std::sort ( std::begin(onebodychan), std::end(onebodychan), [](obc_type& a, obc_type& b){return std::min_element(a.second}  )
  for ( auto& obc : hf.modelspace->OneBodyChannels )
  {
    ljChannel ljchan = {obc.first[0], obc.first[1] };
    
    ljchannels.insert( ljchan );
//    if ( std::find( ljchannels.begin(), ljchannels.end(), ljchan ) == ljchannels.end() ) ljchannels.push_back( ljchan) ;
  }

  int n_mon = 0;
  size_t num_lj = ljchannels.size();
//  for ( size_t ilj_a=0; ilj_a<num_lj; ilj_a++ )
  for ( auto& ljchan_a : ljchannels )
  {
    int la = ljchan_a.l;
    int j2a = ljchan_a.j2;
//    int la = ljchannels[ilj_a].l;
//    int j2a = ljchannels[ilj_a].j2;
//    for ( size_t ilj_b=ilj_a; ilj_b<num_lj; ilj_b++ ) 
    for ( auto& ljchan_b : ljchannels )
    {
      if (ljchan_b < ljchan_a ) continue;
      int lb = ljchan_b.l;
      int j2b = ljchan_b.j2;
//      int lb = ljchannels[ilj_b].l;
//      int j2b = ljchannels[ilj_b].j2;
      int Jab_min = std::abs(j2a-j2b)/2;
      int Jab_max = (j2a+j2b)/2;
//      for ( size_t ilj_c=0; ilj_c<num_lj; ilj_c++ ) // TODO: possibly also change lower limit from 0 to ilj_b, although this becomes more complicated, due to recoupling
      for ( auto& ljchan_c : ljchannels )
      {
        int lc = ljchan_c.l;
        int j2c = ljchan_c.j2;
//        int lc = ljchannels[ilj_c].l;
//        int j2c = ljchannels[ilj_c].j2;

        double t_internal = omp_get_wtime();

//        bool verbose = true;
        bool verbose = false;
        if (verbose) std::cout << "************************ lj channel: "  << la << " " << j2a << "  " << lb << " " << j2b << "  " << lc << " " << j2c << std::endl;

        std::vector<std::array<size_t,8>> imon_indices;
        GetMonopoleIndices(la, j2a, lb, j2b, lc, j2c, hf, imon_indices );

        IMSRGProfiler::timer[std::string(__func__)+"_FindMonKeys"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();

        if (verbose) std::cout << "done.  Now the loop over imon... Size of imonlist is " << imon_indices.size() << std::endl;
//        for (size_t i=0; i<imon_indices.size(); i++)
//        {
//          std::cout << i << " : ";
//          for (int k=0;k<8;k++) std::cout << imon_indices[i][k] << " ";
//          std::cout << std::endl;
//        }
        

        if (imon_indices.size()<1) continue;

//        std::unordered_map<std::string,double> T3bList;  // we will put things into a hash table.
//        GetRelevantTcoeffs(la, j2a, lb, j2b, lc, j2c, hf, T3bList); 
//        GetRelevantTcoeffs(la, j2a, lb, j2b, lc, j2c, hf); 
        GetRelevantTcoeffs(la, j2a, lb, j2b, lc, j2c, hf); 

//        std::unordered_map<std::string,bool> T3b_usedList;
//        for ( auto& iter : T3bList )   T3b_usedList[iter.first] = false;

        t_internal = omp_get_wtime();




        int tcoeff_counter = 0;
        int nonzero_vmon = 0;

//        #pragma omp parallel for schedule(dynamic,1) reduction(+ : tcoeff_counter,nonzero_vmon)
        for (size_t ilist=0; ilist<imon_indices.size(); ilist++)
        {
          size_t imon = imon_indices[ilist][0];

//          std::cout << ")))))))))))))))))))   imon = " << imon << std::endl;
//          if (imon==3 or imon==44) verbose = true;
//          else verbose = false;
//          std::cout << "verbose is " << verbose << std::endl;

          auto key = hf.Vmon3_keys[imon];
//          if (verbose) std::cout << " imon,imon_ab,imon_de,imon_abde = " << imon << " " << imon_ab << " " << imon_de << " " << imon_abde << std::endl;
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          // V3mon(a,b,c,d,e,f) = 1/(2jc+1) * sum_{Jab,J} (2J+1) <abc Jab,J| V |def Jab,J>
          //                    = 1/(2jc+1) * sum_{Jab,J} (2J+1) sum_{jac1,jac2,J12,jac1',jac2',J12',Ncm,Lcm}  Tcoef(abc,Jab,J;jac1,jac2,J12,Ncm,Lcm) * Tcoef(def,Jab,J;jac1',jac2',J12',Ncm,Lcm)
          //                                                                                 sum_{T12,T12',T}  * CG(ta,tb,T12)*CG(tc,T12,T) * CG(td,te,T12')*CG(df,T12',T)   <- Here T means isospin
          //                                                                                                   * < jac1,jac2,J12,T12 | V | jac1',jac2',J12',T12' >
          Orbit& oa = hf.modelspace->GetOrbit(a);
          Orbit& ob = hf.modelspace->GetOrbit(b);
          Orbit& oc = hf.modelspace->GetOrbit(c);
          Orbit& od = hf.modelspace->GetOrbit(d);
          Orbit& oe = hf.modelspace->GetOrbit(e);
          Orbit& of = hf.modelspace->GetOrbit(f);
          int na=oa.n;
          int nb=ob.n;
          int nc=oc.n;
          int nd=od.n;
          int ne=oe.n;
          int nf=of.n;
          int Eabc = 2*(oa.n+ob.n+oc.n) + oa.l+ob.l+oc.l;
          int Edef = 2*(od.n+oe.n+of.n) + od.l+oe.l+of.l;
          if (Eabc > E3max or Edef>E3max) continue;
//          if (Eabc>Nmax or Edef>Nmax) continue;
//          int parity = Eabc%2;

          int Tzab = (oa.tz2+ob.tz2)/2;
          int Tab_min = std::abs(Tzab);
          int twoTz = oa.tz2 + ob.tz2 + oc.tz2;
//          if (twoTz>0) continue;
          int twoT_min = std::abs( twoTz ); // this will either be 1 or 3
//          std::array<double,2> isospin2_Clebsch = {oa.tz2*sqrt(0.5), Tzab + std::abs(oa.tz2-ob.tz2)/2*sqrt(0.5) };
//          std::array<double,4> isospin3_Clebsch = {   1.0,   0.0,  AngMom::CG(1,Tzab,0.5,0.5*oc.tz2, 0.5, 0.5*twoTz),  AngMom::CG(1,Tzab,0.5,0.5*oc.tz2, 1.5,0.5*twoTz) };

          double v_monopole = 0;

          if (verbose) std::cout << std::endl << "===== imon = " << imon << " ==============" << std::endl <<  "abcdef : " << a << " " << b << " " << c << " " << d << " "<< e << " "<< f << std::endl;
          if (verbose) std::cout << "      (" << oa.n << " " <<  oa.l << " " << oa.j2 << " )  ( " << ob.n << " " << ob.l << " " << ob.j2 << " )  ( " << oc.n << " " << oc.l << " " << oc.j2 << std::endl; 

// start new attempt here

          int twoJ_min = std::max(1, j2c-j2a-j2b);
          int twoJ_max = j2a+j2b+j2c;

          for (int Ecm=0; Ecm<=std::min(Eabc,Edef); Ecm++)
          {
           int E12abc = Eabc-Ecm;
           int E12def = Edef-Ecm;
           int parity = E12abc%2;
           if (E12abc > Nmax  or E12def>Nmax) continue;
           for (int twoJ12=1; twoJ12<=2*(Nmax-Ecm)+3; twoJ12+=2)
           {
            // This is clunky, but we'll see how it goes for now...
            size_t dimNAS_abc_T1 = GetDimensionNAS( 1, twoJ12, parity, E12abc ); 
            size_t dimNAS_abc_T3 = GetDimensionNAS( 3, twoJ12, parity, E12abc ); 
            size_t dimNAS_def_T1 = GetDimensionNAS( 1, twoJ12, parity, E12def ); 
            size_t dimNAS_def_T3 = GetDimensionNAS( 3, twoJ12, parity, E12def ); 

            for (int twoT=twoT_min; twoT<=3; twoT+=2)
            {

              auto hashTJN_abc = HashTJN(twoT,twoJ12,E12abc);
              auto hashTJN_def = HashTJN(twoT,twoJ12,E12def);
  
              size_t dimNAS_abc = GetDimensionNAS( twoT, twoJ12, parity, E12abc ); 
              size_t dimNAS_def = GetDimensionNAS( twoT, twoJ12, parity, E12def ); 
              if (dimNAS_abc==0 or dimNAS_def==0) continue;
  
              size_t dimAS_abc = GetDimensionAS( twoT, twoJ12, parity, E12abc ); 
              size_t dimAS_def = GetDimensionAS( twoT, twoJ12, parity, E12def ); 
              if (dimAS_abc==0 or dimAS_def==0) continue;
  
              size_t startloc   = GetStartLocNAS(twoT, twoJ12, E12abc, E12def) ;
              size_t startlocAS = GetStartLocAS( twoT, twoJ12, E12abc, E12def);
  
              size_t cfp_begin_abc = GetCFPStartLocation(twoT,twoJ12,E12abc);
              size_t cfp_begin_def = GetCFPStartLocation(twoT,twoJ12,E12def);


              auto& jacobi_indices_abc = NAS_jacobi_states.at(hashTJN_abc);
              auto& jacobi_indices_def = NAS_jacobi_states.at(hashTJN_def);
  
              arma::mat matelAS( &meAS[startlocAS], dimAS_abc, dimAS_def, false ); 
              
              arma::mat cfp_abc( &(cfpvec[cfp_begin_abc]), dimNAS_abc, dimAS_abc, /*copy_aux_mem*/ false);
              arma::mat cfp_def( &(cfpvec[cfp_begin_def]), dimNAS_def, dimAS_def, /*copy_aux_mem*/ false);
  
              arma::mat matelNAS = cfp_abc * matelAS * cfp_def.t(); // Compute the non-antisymmetrized matrix elements 

//              arma::rowvec Tabc(dimNAS_abc, arma::fill::zeros);
//              arma::vec    Tdef(dimNAS_def, arma::fill::zeros);

              for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
              {
                int twoJ_min = std::abs(twoJ12-2*Lcm);
                int twoJ_max = twoJ12+2*Lcm;

                for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
                {
                  int Jab_min = std::max( std::abs(j2a-j2b), std::abs(twoJ-j2c) )/2;
                  int Jab_max = std::min( j2a+j2b, twoJ+j2c)/2;
                  int Jab_step = 1;
                  if (verbose) std::cout << "Lcm, twoJ12, twoJ, twoT " << Lcm << " " << twoJ12 << " "  << twoJ << " " << twoT << "   E12abc, E12def " << E12abc << " " << E12def << std::endl;


//                  if (verbose) std::cout << "Jabmin,Jabmax,Jabstep = " << Jab_min << " " << Jab_max << " " << Jab_step << std::endl;
                  for (int Jab=Jab_min; Jab<=Jab_max; Jab+=Jab_step)
                  {
                    if (verbose) std::cout << "Jab = " << Jab << "  Tab min = " << Tab_min << std::endl;


                    auto tcoeff_hash_abc = TcoeffHash(na,nb,nc,Jab,twoJ,twoJ12,E12abc);
                    auto tcoeff_hash_def = TcoeffHash(nd,ne,nf,Jab,twoJ,twoJ12,E12def);


                    size_t tcoeff_start_abc = TcoeffLookup[ tcoeff_hash_abc ];
                    size_t tcoeff_start_def = TcoeffLookup[ tcoeff_hash_def ];

                    if (verbose) std::cout << "hashabc : " << tcoeff_hash_abc << "   ->  " << tcoeff_start_abc << "         hashdef : " << tcoeff_hash_def << "    ->  " << tcoeff_start_def << std::endl;

                    size_t offset_abc = (Lcm-Ecm%2)/2 * (2*dimNAS_abc_T1 + dimNAS_abc_T3) + (twoT/2) * 2 * dimNAS_abc_T1;
                    size_t offset_def = (Lcm-Ecm%2)/2 * (2*dimNAS_def_T1 + dimNAS_def_T3) + (twoT/2) * 2 * dimNAS_def_T1;

//                    int rows = (2-Tab_min);
                    int rows = 2-(twoT/2); // T=1/2 -> Tab=0,1   T=3/2 -> Tab=1.

                    if (verbose) std::cout << "Taking Tabc from TcoeffList[ " << tcoeff_start_abc << " + " << offset_abc << " ]"
                                           << "   offset was computed as " << (Lcm-Ecm%2)/2 << " * ( 2* " << dimNAS_abc_T1 << " + " << dimNAS_abc_T3 << " )"
                                           << " + " << (twoT/2) * 2 * dimNAS_abc << std::endl;

//                    arma::mat Tabc( &TcoeffList[tcoeff_start_abc + offset_abc], rows, dimNAS_abc, false ); // false means don't copy the matrix from auxiliary memory to a new location, just use it in-place
//                    arma::mat Tdef( &TcoeffList[tcoeff_start_def + offset_def], rows, dimNAS_def, false );
                    arma::mat Tabc( &TcoeffList[tcoeff_start_abc + offset_abc],  dimNAS_abc, rows, false ); // false means don't copy the matrix from auxiliary memory to a new location, just use it in-place
                    arma::mat Tdef( &TcoeffList[tcoeff_start_def + offset_def],  dimNAS_def, rows, false );
                    arma::mat isospin_mat( rows, rows, arma::fill::eye );


                    if (verbose)
                    {
                      std::cout << "In flat layout, Tabc looks like " << std::endl;
                      for (int i=0; i<dimNAS_abc*rows; i++) std::cout << TcoeffList[tcoeff_start_abc + offset_abc + i] << " ";
                      std::cout << std::endl;
                    }

//                    arma::mat Tabc( rows, dimNAS_abc, arma::fill::zeros );
//                    arma::mat Tdef( dimNAS_def, rows, arma::fill::zeros );
//                    arma::mat isospin_mat( rows, rows, arma::fill::zeros );

//                    bool nonzero_abc = false;
//                    bool nonzero_def = false;
//                    for (size_t iNAS_abc = 0; iNAS_abc<dimNAS_abc; iNAS_abc++)
//                    {
//                      auto& index_1_2_abc = jacobi_indices_abc[iNAS_abc];
//                      auto& jac1= jacobi_1[index_1_2_abc[0]];
//                      auto& jac2= jacobi_2[index_1_2_abc[1]];
//                      int Ncm = (Ecm-Lcm)/2;
//                      int Tab = jacobi_1[index_1_2_abc[0]].t;
//                        double tc = ComputeTcoeff( hf, na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
//                        if (verbose) std::cout << "iNAS_abc = " << iNAS_abc << "  (Tab = " << jac1.t << ", Ncm,Lcm= " << Ncm << " " << Lcm << " )   tcoef = " << tc << std::endl;
//                        if (verbose) std::cout << "   called T( " << na << " " << la << " " << j2a << " " << nb << " " << lb << " " << j2b << " " << nc << " "<< lc << " " << j2c
//                                                                  << " ; " << Jab << " " << twoJ << " |  " << jac1.n << " "<< jac1.l << " " << jac1.s << " "<< jac1.j << ", "
//                                                                  << jac2.n << " " << jac2.l << " " << jac2.j2 << " ; " << twoJ12 << " ,  " << Ncm << " " << Lcm << " ) " << std::endl;
//                      if (oa.n==ob.n and la==lb and j2a==j2b and (Tab+Jab)%2<1) continue;
//                      if (Tab<Tab_min) continue;
//                      int rowT =  Tab-Tab_min;
//                      auto tcoeff_hash = TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,index_1_2_abc[0],index_1_2_abc[1],twoJ12,Lcm );
//                      Tabc(rowT,iNAS_abc) = T3bList.at(tcoeff_hash) ;
//                      if (verbose) std::cout << "   Tabc[" << rowT << "," << iNAS_abc << "]  = " << T3bList[ TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,index_1_2_abc[0],index_1_2_abc[1],twoJ12,Lcm )]
//                                               << std::endl;
//                      tcoeff_counter++;
//                      nonzero_abc = true;
//                    }
//                    for (size_t iNAS_def = 0; iNAS_def<dimNAS_def; iNAS_def++)
//                    {
//                      auto& index_1_2_def = jacobi_indices_def[iNAS_def];
//                      int Tde = jacobi_1[index_1_2_def[0]].t;
//                      if (od.n==oe.n and la==lb and j2a==j2b and (Tde+Jab)%2<1) continue;
//                      if (Tde<Tab_min) continue;
//                      int rowT =  Tde-Tab_min;
//                      auto tcoeff_hash = TcoeffHash(od.n,oe.n,of.n,Jab,twoJ,index_1_2_def[0],index_1_2_def[1],twoJ12,Lcm );
//                      if (verbose) std::cout << " Jab,Tde = " << Jab << " " << Tde << "  colJT = " << rowT << " iNAS_def = " << iNAS_def << "  dimNAS_def = " << dimNAS_def << std::endl;
//                      Tdef(iNAS_def,rowT) = T3bList.at(tcoeff_hash) ;//* (twoJ+1)  ;
//                      if (verbose) std::cout << "   Tdef[" << iNAS_def << "," << rowT << "]  = " << T3bList[ TcoeffHash(od.n,oe.n,of.n,Jab,twoJ,index_1_2_def[0],index_1_2_def[1],twoJ12,Lcm )]
//                                             <<  std::endl;
//                      tcoeff_counter++;
//                      nonzero_def = true;
//                    }
//                    if ( not (nonzero_abc and nonzero_def) ) continue;

//                   if (twoT==1 )
//                   {
//                    for (int Tab=Tab_min; Tab<=1; Tab++)
                    for (int Tab=twoT/2; Tab<=1; Tab++)
                    {
//                      Tabc.row(Tab) *= isospin2_Clebsch[Tab] * isospin3_Clebsch[2*Tab + twoT/2];
//                      Tdef.row(Tab) *= isospin2_Clebsch[Tab] * isospin3_Clebsch[2*Tab + twoT/2];
                      isospin_mat.row(Tab-twoT/2) *= AngMom::CG(0.5,0.5*oa.tz2,0.5,0.5*ob.tz2, Tab, Tzab) * AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2, 0.5*twoT, 0.5*twoTz) ;
//                      isospin_mat.row(Tab-twoT/2) *= isospin2_Clebsch[Tab] * isospin3_Clebsch[2*Tab + twoT/2];
//                      if (verbose) std::cout << "Tab = " << Tab << "  multiplying row " << Tab-twoT/2 << " by " << isospin2_Clebsch[Tab] << " * " << isospin3_Clebsch[2*Tab + twoT/2] << std::endl;
                      if (verbose) std::cout << "Tab = " << Tab << "  multiplying row " << Tab-twoT/2 << " by " <<  AngMom::CG(0.5,0.5*oa.tz2,0.5,0.5*ob.tz2, Tab, Tzab) << " * " << AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2, 0.5*twoT, 0.5*twoTz) << std::endl;
                      if (verbose) std::cout << " clebsch2 = " <<  AngMom::CG(0.5,0.5*oa.tz2,0.5,0.5*ob.tz2, Tab, Tzab) << std::endl;
                      if (verbose) std::cout << " clebsch3 = " <<  AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2, 0.5*twoT, 0.5*twoTz) << std::endl;
//                      Tabc.row(Tab-Tab_min) *= isospin2_Clebsch[Tab] * isospin3_Clebsch[2*Tab + twoT/2];
//                      Tdef.col(Tab-Tab_min) *= isospin2_Clebsch[Tab] * isospin3_Clebsch[2*Tab + twoT/2];
                    }
//                   }

//                  arma::mat result =  Tabc * matelNAS * Tdef ;
//                  arma::mat result =  Tabc * matelNAS * Tdef.t() ;
//                  arma::mat result =  isospin_mat * Tabc * matelNAS * Tdef.t() * isospin_mat  ;
                  arma::mat result =  isospin_mat * Tabc.t() * matelNAS * Tdef * isospin_mat  ;
                  v_monopole += 6* (twoJ+1) * arma::accu( result ) ;
                  if (verbose) std::cout << "  twoT, Ecm, twoJ12, Lcm, twoJ: " << twoT << " " << Ecm << " " << twoJ12 << "  " << Lcm << " " << twoJ << std::endl;
                  if (verbose) std::cout << "matrices:" << std::endl << std::endl << Tabc.t() << std::endl << std::endl << matelNAS << std::endl << std::endl << Tdef << std::endl << std::endl;
                  if (verbose) std::cout << "isospin matrix : " << std::endl << isospin_mat << std::endl;
//                  if (verbose) std::cout << std::endl << "intermediate: " << std::endl << (matelNAS * Tdef.t()) << std::endl << std::endl;
                  if (verbose) std::cout << std::endl << "intermediate: " << std::endl << (matelNAS * Tdef) << std::endl << std::endl;
                  if (verbose) std::cout << "accumulate : " << arma::accu(result ) << std::endl;
                  if (verbose) std::cout << " result, vmonopole = " << std::endl << result << std::endl  << v_monopole << std::endl;
                  } // for Jab
                } // for twoJ



              } // for Lcm

            } // for twoJ12
           } // for Ecm
          } // for twoT

          v_monopole *= 1./(j2c+1.);
          if (verbose) std::cout << "multiplying 6/2jc+1 : " << v_monopole << std::endl;
          if ( std::abs( v_monopole)>1e-8 ) nonzero_vmon++;
         // There are some symmetries we can exploit here to avoid redundant calculations
          for ( auto& imon_sym : imon_indices[ilist] )
          {
              hf.Vmon3[imon_sym] = v_monopole;
          }


        } // for imon

        n_mon += imon_indices.size();
//        std::cout << "done " << std::endl;
        if (verbose) std::cout << "Looked up a Tcoefficient " << tcoeff_counter << "  times " << std::endl;
//        std::cout << "Used " << iused << " different T coefficients " << std::endl;
        if (verbose) std::cout << "Found " << nonzero_vmon << "  nonzero monopoles" << "  out of " << imon_indices.size() << " terms" << std::endl;
        if (verbose) std::cout << " and it took " << omp_get_wtime() - t_internal << "  seconds" << std::endl;
//        if (tcoeff_counter < 1  and T3bList.size()>0)
//        {
//          std::cout << "@@@@@@@@@@@@  channel = (" << la << " " << j2a << " , " << lb << " " << j2b << " , " << lc << " " << j2c << std::endl;
//        }

        IMSRGProfiler::timer[std::string(__func__)+"_iMonLoop"] += omp_get_wtime() - t_internal;

      }// for ilj_c
    }// for ilj_b
  }// for ilj_a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  IMSRGProfiler::counter["V3Monopoles"] += n_mon;
}



////////////////////////////////////////////////////////////////////////////////////
// FOR REFERENCE, HERE's WHAT THE HF CODE IS DOING
// WE WANT TO PROVIDE THE MATRIX V3NO FOR A SINGLE 2-BODY CHANNEL
//
//    for (int ch=0;ch<nchan;++ch)
//    {
//       TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
//       int J = tbc.J;
//       int npq = tbc.GetNumberKets();
// 
//       arma::mat D(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
//       arma::mat V3NO(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
//
//       ....
// 
//       #pragma omp parallel for schedule(dynamic,1) // confirmed that this improves performance
//       for (int i=0; i<npq; ++i)    
//       {
//          Ket & bra = tbc.GetKet(i);
//          int e2bra = 2*bra.op->n + bra.op->l + 2*bra.oq->n + bra.oq->l;
//          for (int j=0; j<npq; ++j)
//          {
//            if (i>j) continue;
//            for ( auto a : modelspace->all_orbits )
//            {
//              Orbit & oa = modelspace->GetOrbit(a);
//              if ( 2*oa.n+oa.l+e2bra > Hbare.GetE3max() ) continue;
//              for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
//              {
//                Orbit & ob = modelspace->GetOrbit(b);
//                if ( 2*ob.n+ob.l+e2ket > Hbare.GetE3max() ) continue;
//                if ( std::abs(rho(a,b)) < 1e-8 ) continue; // Turns out this helps a bit (factor of 5 speed up in tests)
//                int J3min = std::abs(2*J-oa.j2);
//                int J3max = 2*J + oa.j2;
//                for (int J3=J3min; J3<=J3max; J3+=2)
//                {
//                  V3NO(i,j) += rho(a,b) * (J3+1) * Hbare.ThreeBody.GetME_pn(J,J,J3,bra.p,bra.q,a,ket.p,ket.q,b);
//                }
//              }
//            }
//            V3NO(i,j) /= (2*J+1);
//            if (bra.p==bra.q)  V3NO(i,j) /= SQRT2; 
//            if (ket.p==ket.q)  V3NO(i,j) /= SQRT2; 
//            V3NO(j,i) = V3NO(i,j);
//
///////////////////////////////////////////////////////////////////////////////
void Jacobi3BME::GetNO2b_single_channel( HartreeFock& hf, int ch, arma::mat& V3NO )
{
  std::cout << "Enter " <<__func__  << " with  ch = " << ch << "  number of threads: " << omp_get_num_threads() << std::endl;
  double t_start = omp_get_wtime();
  TwoBodyChannel& tbc = hf.modelspace->GetTwoBodyChannel(ch);
  int Jab = tbc.J;
  int nkets = tbc.GetNumberKets();

//  std::cout << "size of V3NO = " << V3NO.n_rows << " x " << V3NO.n_cols  << ".  nkets = " << nkets << ".  size of rho = " << hf.rho.n_rows << " x " << hf.rho.n_cols << std::endl;
  V3NO.zeros( nkets,nkets);

//  if (ch==1) std::cout << " Jab = " << Jab << std::endl;

  std::set<size_t> occupied_orbits;
  for (size_t irow=0; irow<hf.rho.n_rows; irow++)
  {
    if ( arma::norm(hf.rho.row(irow),"fro")>1e-8) occupied_orbits.insert(irow);
  }


  double t_internal = omp_get_wtime();

  std::unordered_map<std::string,arma::mat> TcoeffTable;
  std::unordered_map<std::string,bool> TcoeffSkip;

  // We twice pass through the list of T coefficients that we'll need for this calculation
  // On the first pass, which is single-threaded, we allocate the structure
  // On the second pass, which is multi-threaded, we compute the T coefficients to put them into the structure
  for (int loop_pass=0; loop_pass<2;loop_pass++)
  {
   #pragma omp parallel for  schedule(dynamic,1) if (loop_pass>0)   // don't go parallel on the first pass
   for (int ibra=0; ibra<nkets; ibra++)
   {
     Ket& bra = tbc.GetKet(ibra);
     int a = bra.p;
     int b = bra.q;
     Orbit& oa = hf.modelspace->GetOrbit(a);
     Orbit& ob = hf.modelspace->GetOrbit(b);
     if ( 2*(oa.n+ob.n)+oa.l+ob.l > E2max  ) continue;
     for (auto c : occupied_orbits )
     {
       Orbit& oc = hf.modelspace->GetOrbit(c);
      // loop over   twoJ,   Ecm,  twoJ12,  twoT,   Lcm
       int Eabc = 2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l;
       if (Eabc>E3max) continue;
       int twoTz = oa.tz2 + ob.tz2 + oc.tz2;
       int twoT_min = std::abs(twoTz);
       int twoJ_min = std::abs( 2*Jab-oc.j2);
       int twoJ_max = 2*Jab+oc.j2;
       for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
       {
         for (int Ecm=0; Ecm<=Eabc; Ecm++)
         {
           int E12abc = Eabc-Ecm;
           if (E12abc>Nmax) continue;
           int parity = E12abc%2;
           int twoJ12_min = std::max( 1, (twoJ - 2*Ecm) ); // maybe take a second look at these limits...
           int twoJ12_max = std::min( twoJmax, twoJ +2*Ecm);
           for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
           {
             for ( int twoT=twoT_min; twoT<=3; twoT+=2)
             {
               size_t dimNAS = GetDimensionNAS( twoT, twoJ12, parity, E12abc ); 
               if (dimNAS==0 ) continue;
               size_t dimAS = GetDimensionAS( twoT, twoJ12, parity, E12abc ); 
               if (dimAS==0 ) continue;
               auto hashTJN = HashTJN(twoT,twoJ12,E12abc);
               int rows = 2-(twoT/2); // T=1/2 -> Tab=0,1   T=3/2 -> Tab=1.
               auto& jacobi_indices = NAS_jacobi_states.at(hashTJN);
//               for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
               for (int Lcm=std::max(Ecm%2,std::abs(twoJ-twoJ12)/2); Lcm<=std::min(Ecm,(twoJ+twoJ12)/2); Lcm+=2)
               {
//                 if ( std::abs(twoJ-twoJ12)>2*Lcm or (twoJ+twoJ12)<2*Lcm) continue;
                 int Ncm=(Ecm-Lcm)/2;
                 std::ostringstream oss;
                 oss << ibra << " " << c << " " << twoJ << " " << Ecm << " " << twoJ12 << " " << twoT << " " << Lcm;
                 std::string tcoeff_hash = oss.str();

                 if (loop_pass==0) // first pass, not parallel, just make space
                 {
                   TcoeffTable.emplace( tcoeff_hash,   arma::mat( dimNAS, rows, arma::fill::zeros ) );
                   TcoeffSkip.emplace( tcoeff_hash,   false );
                 }
                 else  // second pass, in parallel, compute the T coefficients and place them in the data structure
                 {
                  auto& Tcoef_mat = TcoeffTable[tcoeff_hash];
                  for (int iNAS=0;iNAS<dimNAS; iNAS++)
                  {
                    auto& index_1_2 = jacobi_indices[iNAS];
                    auto& jac1= jacobi_1[index_1_2[0]];
                    auto& jac2= jacobi_2[index_1_2[1]];
                    double tcoef = ComputeTcoeff( oa.n, oa.l, oa.j2, ob.n, ob.l, ob.j2, oc.n, oc.l, oc.j2, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
                    Tcoef_mat( iNAS, jac1.t-twoT/2 ) = tcoef;
                  }
                  if (arma::norm( Tcoef_mat, "fro")<1e-9 ) TcoeffSkip[tcoeff_hash]=true;
                 }
               }
             }
           }
         }
       }
     }
   }
  }
  // all done precomputing
  IMSRGProfiler::timer[std::string(__func__)+"_computeT"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();

//  #pragma omp parallel for schedule(dynamic,1)
  for (int ibra=0; ibra<nkets; ibra++)
  {
    Ket& bra = tbc.GetKet(ibra);
//    if(ch==1) std::cout << "ibra = " << ibra << "  pq = " << bra.p << " " << bra.q << std::endl;
    int a = bra.p;
    int b = bra.q;
    Orbit& oa = hf.modelspace->GetOrbit(a);
    Orbit& ob = hf.modelspace->GetOrbit(b);
    if ( 2*(oa.n+ob.n)+oa.l+ob.l > E2max  ) continue;
    int Tzab = (oa.tz2+ob.tz2)/2;

    for (int iket=0; iket<=ibra; iket++ )
    {
      Ket& ket = tbc.GetKet(iket);
//      if(ch==1)std::cout << " iket = " << iket << std::endl;
      int d = ket.p;
      int e = ket.q;
      Orbit& od = hf.modelspace->GetOrbit(d);
      Orbit& oe = hf.modelspace->GetOrbit(e);
      if ( 2*(od.n+oe.n)+od.l+oe.l > E2max  ) continue;


//     for (auto c : hf.modelspace->holes )
     for (auto c : occupied_orbits )
     {
//      std::cout << "  c = " << c << std::endl;
      Orbit& oc = hf.modelspace->GetOrbit(c);
      int twoJ_min = std::abs( 2*Jab - oc.j2 );
      int twoJ_max = (2*Jab+oc.j2);
      int twoTz = oa.tz2 + ob.tz2 + oc.tz2;
      int twoT_min = std::abs(twoTz);
      if ( 2*(oa.n+oc.n)+oa.l+oc.l > E2max  ) continue;
      if ( 2*(oc.n+ob.n)+oc.l+ob.l > E2max  ) continue;
      int Eabc = 2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l;
      if (Eabc>E3max) continue;

//     if (ch==1)
//     {
//       std::cout << " Looping over f states :  ";
//       for (auto h : hf.modelspace->OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) ) std::cout << h << " ";
//       std::cout << std::endl;
//     }
       for (auto f : hf.modelspace->OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) )
       {
//         std::cout << "   f = " << f << std::endl;
//         if (ch==1) std::cout << " Jab, bra,ket  abcdef  " << Jab << " " << ibra << " " << iket << "  " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
         if ( std::abs( hf.rho(c,f)) <1e-8) continue;
         Orbit& of = hf.modelspace->GetOrbit(f);
         if ( 2*(od.n+of.n)+od.l+of.l > E2max  ) continue;
         if ( 2*(of.n+oe.n)+of.l+oe.l > E2max  ) continue;
         int Edef = 2*(od.n+oe.n+of.n)+od.l+oe.l+of.l;
         if (Edef>E3max) continue;

         double v_no2b_cf = 0;

         // Inner loops over jacobi stuff
         for (int Ecm=0; Ecm<=std::min(Eabc,Edef); Ecm++)
         {
           int E12abc = Eabc - Ecm;
           int E12def = Edef - Ecm;
           if (E12abc>Nmax or E12def>Nmax) continue;
           int parity = E12abc%2;
           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
           {
            double vsum_J = 0;
            int twoJ12_min = std::max( 1, (twoJ - 2*Ecm) ); // maybe take a second look at these limits...
            int twoJ12_max = std::min( twoJmax, twoJ +2*Ecm);
            for ( int twoT=twoT_min; twoT<=3; twoT+=2)
            {
                int rows = 2-(twoT/2); // T=1/2 -> Tab=0,1   T=3/2 -> Tab=1.

                arma::mat isospin_mat_abc( rows, rows, arma::fill::eye );
                arma::mat isospin_mat_def( rows, rows, arma::fill::eye );

                for (int Tab=twoT/2; Tab<=1; Tab++)
                {   // because c and f are in the same one-body channel,  Tzab = Tzde, although Tab need not be Tde, and tza need not be tzd, etc.
                  isospin_mat_abc(Tab-twoT/2,Tab-twoT/2) = AngMom::CG(0.5,0.5*oa.tz2,0.5,0.5*ob.tz2, Tab, Tzab) * AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2, 0.5*twoT, 0.5*twoTz) ;
                  isospin_mat_def(Tab-twoT/2,Tab-twoT/2) = AngMom::CG(0.5,0.5*od.tz2,0.5,0.5*oe.tz2, Tab, Tzab) * AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2, 0.5*twoT, 0.5*twoTz) ;
                }

              for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
              {
                
                auto hashTJN_abc = HashTJN(twoT,twoJ12,E12abc);
                auto hashTJN_def = HashTJN(twoT,twoJ12,E12def);
    
                size_t dimNAS_abc = GetDimensionNAS( twoT, twoJ12, parity, E12abc ); 
                size_t dimNAS_def = GetDimensionNAS( twoT, twoJ12, parity, E12def ); 
                if (dimNAS_abc==0 or dimNAS_def==0) continue;
    
                size_t dimAS_abc = GetDimensionAS( twoT, twoJ12, parity, E12abc ); 
                size_t dimAS_def = GetDimensionAS( twoT, twoJ12, parity, E12def ); 
                if (dimAS_abc==0 or dimAS_def==0) continue;
    
                size_t startloc   = GetStartLocNAS(twoT, twoJ12, E12abc, E12def) ;
                size_t startlocAS = GetStartLocAS( twoT, twoJ12, E12abc, E12def);
    
                size_t cfp_begin_abc = GetCFPStartLocation(twoT,twoJ12,E12abc);
                size_t cfp_begin_def = GetCFPStartLocation(twoT,twoJ12,E12def);
  
  
                auto& jacobi_indices_abc = NAS_jacobi_states.at(hashTJN_abc);
                auto& jacobi_indices_def = NAS_jacobi_states.at(hashTJN_def);
    
                arma::mat matelAS( &meAS[startlocAS], dimAS_abc, dimAS_def, false ); 
                
                arma::mat cfp_abc( &(cfpvec[cfp_begin_abc]), dimNAS_abc, dimAS_abc, /*copy_aux_mem*/ false);
                arma::mat cfp_def( &(cfpvec[cfp_begin_def]), dimNAS_def, dimAS_def, /*copy_aux_mem*/ false);
    
                arma::mat matelNAS = 6 * cfp_abc * matelAS * cfp_def.t(); // Compute the non-antisymmetrized matrix elements 
 

 
                for (int Lcm=std::max(Ecm%2,std::abs(twoJ-twoJ12)/2); Lcm<=std::min(Ecm,(twoJ+twoJ12)/2); Lcm+=2)
                {
                  int Ncm=(Ecm-Lcm)/2;
                  int rows = 2-(twoT/2); // T=1/2 -> Tab=0,1   T=3/2 -> Tab=1.
//                  arma::mat Tabc( dimNAS_abc, rows, arma::fill::zeros ); 
//                  arma::mat Tdef( dimNAS_def, rows, arma::fill::zeros );

                  std::ostringstream oss_abc, oss_def;
                  oss_abc << ibra << " " << c << " " << twoJ << " " << Ecm << " " << twoJ12 << " " << twoT << " " << Lcm;
                  oss_def << iket << " " << f << " " << twoJ << " " << Ecm << " " << twoJ12 << " " << twoT << " " << Lcm;
                  std::string t_hash_abc = oss_abc.str();
                  std::string t_hash_def = oss_def.str();
                  if (  TcoeffSkip[t_hash_abc] or TcoeffSkip[t_hash_def]) continue;  // it's not clear that this will help much

                  arma::mat& Tabc = TcoeffTable[ t_hash_abc ];
                  arma::mat& Tdef = TcoeffTable[ t_hash_def ];

                  
                  arma::mat result =  isospin_mat_abc * Tabc.t() * matelNAS * Tdef * isospin_mat_def  ;
//                  if(ch==1) std::cout << "J,Ecm,T,Lcm = " << twoJ << " " << Ecm << " " << twoT << " " << Lcm << "     matrices: "
//                                      << std::endl << Tabc.t() << std::endl << matelNAS << std::endl << Tdef << std::endl << result << std::endl << " isospin mats" << std::endl << isospin_mat_abc << std::endl << isospin_mat_def <<std::endl << std::endl;
                  
//                  std::cout << "     bra,ket: " << t_hash_abc << " | " << t_hash_def << "   :  " << arma::accu( result ) << std::endl;
//                  if (  arma::norm( result, "fro") < 1e-7 )
//                  {
//                     std::cout << Tabc << std::endl << Tdef << std::endl;
//                  }
                  vsum_J += arma::accu( result ) ;

                } // for Lcm

              } // for twoJ12
             } // for twoT
            
            v_no2b_cf +=  vsum_J * (twoJ+1);
           } // for twoJ
          } // for Ecm
//          if(ch==1) std::cout << "  v_no2b_cd +=" << vsum_J  << " * " << twoJ+1 << "  -> " << v_no2b_cf << std::endl;
//         V3NO(ibra,iket) += v_no2b_cf * hf.rho(c,f); 
         V3NO(ibra,iket) += v_no2b_cf * hf.rho(c,f); 
//         if(ch==1) std::cout << "   V3NO += " << v_no2b_cf << " * " << hf.rho(c,f) << "  ->  " << V3NO(ibra,iket) << std::endl;
       } // for f
     } // for c
     if (bra.p==bra.q) V3NO(ibra,iket) /= SQRT2;
     if (ket.p==ket.q) V3NO(ibra,iket) /= SQRT2;
     V3NO(iket,ibra) = V3NO(ibra,iket);
    } // for iket
  } // for ibra
  V3NO /= 2*Jab+1;
  std::cout << "Done with NO2B loop" << std::endl;


  IMSRGProfiler::timer[std::string(__func__)+"_computeNO2B"] += omp_get_wtime() - t_internal;
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
}







// This is just a convenience wrapper interface so that we don't have to pass so damn many parameters to Tcoeff
double Jacobi3BME::Tcoeff_wrapper( Ket3& ket, int Jab, int twoJ, jacobi1_state& jac1, jacobi2_state& jac2, int twoJ12, int Ncm, int Lcm)
{
  return AngMom::Tcoeff( ket.op->n, ket.op->l, ket.op->j2,
                         ket.oq->n, ket.oq->l, ket.oq->j2,
                         ket.oR->n, ket.oR->l, ket.oR->j2,  Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
}






void Jacobi3BME::ReadTcoeffNavratil( std::string fname, std::vector<double>& tcoeff, std::vector<labstate_nav>& labst, std::vector<relativestate_nav>& relst, std::vector<bookkeeping_nav>& bookkeeping)
{
  // First, we populate the lab and relative states
  labst.resize(0);
  relst.resize(0);
//  bookkeeping.resize( (Nmax+3) * Nmax+1 );
  bookkeeping.resize( Nmax*(Nmax+2) + Nmax+3 );
  std::cout << "Begin lab states. size of bookkeeping is " << bookkeeping.size() << std::endl;

/*
         DO jtot=1,2*nhom+3,2   ! la+lb+lc+3/2; jtot half integer
          DO ntot=max((jtot-3)/2,0),nhom      ! total N
           DO la=0,ntot         
            DO ja=abs(2*la-1),2*la+1,2 ! ja half integer
             DO lb=0,ntot-la
              DO jb=abs(2*lb-1),2*lb+1,2 ! jb half integer
               DO lc=mod(ntot+la+lb,2),ntot-la-lb,2 ! parity fixed by ntot
                DO jc=abs(2*lc-1),2*lc+1,2 ! jb half integer
                 DO bj12=abs((ja-jb))/2,(ja+jb)/2 ! bj12 is integer
                  IF(jtot.GE.abs(2*bj12-jc).AND.    jtot.LE.2*bj12+jc) THEN
                    DO na=0,(ntot-la-lb-lc)/2
                     DO nb=0,(ntot-la-lb-lc)/2
                      IF(ntot-la-lb-lc-2*na-2*nb  .GE.0) THEN
                        nc=(ntot-la-lb-lc-2*na-2*nb)/2
                        IF(2*na+la+2*nb+lb.LE.nhom2sp.AND. 2*nc+lc+2*nb+lb.LE.nhom2sp.AND. 2*nc+lc+2*na+la.LE.nhom2sp.AND. 2*na+la.LE.nhom1sp.AND. 2*nb+lb.LE.nhom1sp.AND. 2*nc+lc.LE.nhom1sp) THEN
                          IF(ntot.EQ.2*na+2*nb+2*nc+la+lb+lc) THEN
                            alphaspmax=alphaspmax+1
*/

  // lab states
  for (int jtot=1; jtot<=2*Nmax+3; jtot+=2 ) // jtot is twice the half-integer total J of the state
  {
   for (int ntot=std::max((jtot-3)/2,0); ntot<=Nmax; ntot++) // total harmonic oscillator quanta e1+e2+e3
   {
//    std::cout << "Trying to access bookkeeping[" << ntot*(Nmax+2) + jtot/2 << " ] " << std::endl;
    auto& bk = bookkeeping.at( ntot*(Nmax+2)+(jtot/2) );
    bk.nspstart = labst.size();
    for (int la=0;la<=ntot;la++)
    {
     for (int ja=std::abs(2*la-1); ja<=2*la+1; ja+=2 ) // ja twice the half-integer j
     {
      for (int lb=0; lb<=ntot-la; lb++)
      {
       for (int jb=std::abs(2*lb-1); jb<=2*lb+1; jb+=2)
       {
        for (int lc=(ntot+la+lb)%2; lc<=ntot-la-lb; lc+=2)
        {
         for (int jc=std::abs(2*lc-1); jc<=2*lc+1; jc+=2)
         {
          for (int bj12=std::abs(ja-jb)/2; bj12<=(ja+jb)/2; bj12++) // ja and jb couple to bj12, which is integer so we don't multiply by 2.
          {
           if ( std::abs(2*bj12-jc)>jtot or (2*bj12+jc)<jtot ) continue;
           for (int na=0; na<=(ntot-la-lb-lc)/2; na++)
           {
            for (int nb=0; nb<=(ntot-la-lb-lc)/2; nb++)
            {
             if ((ntot-la-lb-lc-2*na-2*nb)<0) continue;
             int nc=(ntot-2*na-2*nb-la-lb-lc)/2;
             if (nc<0) continue;
             if ( ntot != 2*(na+nb+nc)+la+lb+lc ) continue;
             if ( ((2*na+la) > emax) or ((2*nb+lb)>emax) or ((2*nc+lc)>emax) ) continue;
             if ( ((2*na+la+2*nb+lb)>E2max) or ((2*na+la+2*nc+lc)>E2max) or ((2*nb+lb+2*nc+lc)>E2max) ) continue;
             labst.push_back( {na,nb,nc,la,lb,lc,ja,jb,jc,bj12,jtot} );
             if (ntot==2 and jtot==1)
             {
               std::cout << std::fixed << std::setw(8) << labst.size() - bk.nspstart << " :  " << na << " " << la << " " << ja << "  " << nb << " " << lb << " " << jb << "  " << nc << " " << lc << " " << jc << "    " << bj12 << "  " << jtot << std::endl;
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
    bk.nspnum = labst.size() - bk.nspstart;
   }
  }

  std::cout << "Now the relative states" << std::endl;
  // now the relative states
  int nummat = 0;
  for (int jtot=1; jtot<=2*Nmax+3; jtot+=2 ) // jtot is twice the half-integer total J of the state
  {
   for (int ntot=std::max((jtot-3)/2,0); ntot<=Nmax; ntot++) // total harmonic oscillator quanta e1+e2+e3
   {
    auto& bk = bookkeeping.at( ntot*(Nmax+2)+(jtot/2) );
    bk.nrelstart = relst.size();
    bk.matstart = nummat;
    for (int ncm=0; ncm<=ntot/2; ncm++)
    {
     for (int lcm=0; lcm<=ntot-2*ncm; lcm++)
     {
      for (int j3=std::max(1,std::abs(jtot-2*lcm)); j3<=std::min(2*(ntot-2*ncm-lcm)+3,jtot+2*lcm); j3+=2) // j3 is twice the total relative angular momentum
      {
       for (int l3=0; l3<=ntot-2*ncm-lcm; l3++) // orbital angular momentum of the second jacobi coordinate
       {
        for (int I3=std::abs(2*l3-1); I3<=2*l3+1; I3+=2) // I3 is twice the half-integer value, angular momentum of l3 coupled with spin 1/2
        {
         for (int j12=std::abs(j3-I3)/2; j12<=(j3+I3)/2; j12++) // j12 is integer, so no factor of 2
         {
          for (int s12=0; s12<=1; s12++)
          {
           for (int l12=std::abs(s12-j12); l12<=s12+j12; l12++)
           {
            if ( (ntot+l12+l3+lcm)%2 >0 ) continue; // check parity
            for (int n12=0; n12<=(ntot-2*ncm-l12-l3-lcm)/2; n12++)
            {
             int n3=(ntot-l12-l3-lcm-2*n12-2*ncm)/2;
             if ( n3<0 or ( (2*n3+2*n12+2*ncm + l3+l12+lcm) != ntot ) ) continue; // check energy conservation
             relst.push_back( {n12,n3,ncm,l12,l3,lcm,j12,j3,I3,s12,jtot} );
            }
           }
          }
         }
        }
       }
      }
     }
    }
    bk.nrelnum = relst.size() - bk.nrelstart;
    nummat += bk.nrelnum * bk.nspnum;
   }
  }
  auto& bk = bookkeeping.at( Nmax*(Nmax+2) + (2*Nmax+3)/2) ;

  tcoeff.resize( bk.matstart + bk.nspnum * bk.nrelnum + 1 );


  // at last we can actually read the file
  std::ifstream t_file(fname, std::ios::binary );
  for (int jtot=1; jtot<=2*Nmax+3; jtot+=2 ) // jtot is twice the half-integer total J of the state
  {
   if (not t_file.good() ) break;
   for (int ntot=std::max((jtot-3)/2,0); ntot<=Nmax; ntot++) // total harmonic oscillator quanta e1+e2+e3
   {
     if (not t_file.good() ) break;
     uint32_t delimiter;
     uint32_t ntot_in, jtot_in; // these are declared as integer, which probably means 32bit?
     uint64_t dim, start; // these are declared as integer, which probably means 32bit?
     t_file.read((char*)&delimiter,    sizeof(delimiter));
     t_file.read((char*)&ntot_in,      sizeof(ntot_in));
     t_file.read((char*)&jtot_in,      sizeof(jtot_in));
     t_file.read((char*)&dim,          sizeof(dim));
     t_file.read((char*)&start,        sizeof(start));
     t_file.read((char*)&delimiter,    sizeof(delimiter));

     std::cout << " ntot,jtot = " << ntot << " " << jtot << std::endl;
     auto& bk = bookkeeping.at( ntot*(Nmax+2)+(jtot/2) );
     if (bk.matstart != start-1)
     {
       std::cout << __func__ << "  DANGER!! start != matstart : " << start << " != " << bk.matstart << std::endl;
     }
     if (dim != (bk.nspnum * bk.nrelnum) )
     {
       std::cout << __func__ << "  DANGER!! dim != nspnum*nrelnum : " << dim << " != " << bk.nspnum << " * " << bk.nrelnum << " = " << bk.nspnum*bk.nrelnum << std::endl;
     }

     start -=1;

     if ( tcoeff.size() < start + dim )
     { 
      std::cout << "resizing to " << start << " + " << dim << std::endl;
      tcoeff.resize( start + dim +1 );
     }

     // now we read dim coefficients
     float dummy_t;
     t_file.read((char*)&delimiter,    sizeof(delimiter));
     for (int i=0; i<dim; i++)
     {
       if (not t_file.good() ) break;
       t_file.read((char*)&(dummy_t),  sizeof(dummy_t));
       tcoeff[start+i] = dummy_t;
//       std::cout << "   tcoeff " << i << " =  " << tcoeff[start+i] << std::endl;
     }
     t_file.read((char*)&delimiter,    sizeof(delimiter));
     
   }
  }

  std::cout << "Cool. All done reading. The sizes of the lab and rel states are " << labst.size() << "  " << relst.size() << std::endl;
  std::cout << " last bookkeeping entry : " << bk.matstart << " " << bk.nspstart << " " << bk.nspnum << "   " << bk.nrelstart << " " << bk.nrelnum << std::endl;

}



void Jacobi3BME::TestReadTcoeffNavratil(std::string fname )
{

 std::vector<double> tcoeff;
 std::vector<labstate_nav> labst;
 std::vector<relativestate_nav> relst;
 std::vector<bookkeeping_nav> bookkeeping;

 ReadTcoeffNavratil(  fname, tcoeff,  labst, relst, bookkeeping );
 std::cout << "All done reading. Now what have we done?" << std::endl;

  for (int jtot=1; jtot<=2*Nmax+3; jtot+=2 ) // jtot is twice the half-integer total J of the state
  {
//   if (jtot>3) continue;
   for (int ntot=std::max((jtot-3)/2,0); ntot<=Nmax; ntot++) // total harmonic oscillator quanta e1+e2+e3
   {
//     if (ntot>1) continue;
     std::cout << " jtot,ntot = " << jtot << " " << ntot << std::endl;
     auto& bk = bookkeeping.at( ntot*(Nmax+2)+(jtot/2) );
     for (int isp=bk.nspstart; isp< bk.nspstart+bk.nspnum; isp++)
     {
       auto& lab = labst[isp];
       for (int irel=bk.nrelstart; irel< bk.nrelstart+bk.nrelnum; irel++)
       {
         auto& rel = relst[irel];
         int index = bk.matstart + bk.nrelnum * (isp-bk.nspstart) + (irel-bk.nrelstart) ;
         double tc = tcoeff.at(index);
         double mytc = AngMom::Tcoeff( lab.na,lab.la,lab.ja, lab.nb,lab.lb,lab.jb, lab.nc,lab.lc,lab.jc, lab.bj12, jtot,
                                       rel.n12,rel.l12,rel.s12, rel.j12, rel.n3, rel.l3, rel.I3, rel.j3, rel.ncm, rel.lcm);
         double mytc_bf = AngMom::Tcoeff_bruteforce( lab.na,lab.la,lab.ja, lab.nb,lab.lb,lab.jb, lab.nc,lab.lc,lab.jc, lab.bj12, jtot,
                                                  rel.n12,rel.l12,rel.s12, rel.j12, rel.n3, rel.l3, rel.I3, rel.j3, rel.ncm, rel.lcm);
         double mytc_fancy = AngMom::Tcoeff_fancy( lab.na,lab.la,lab.ja, lab.nb,lab.lb,lab.jb, lab.nc,lab.lc,lab.jc, lab.bj12, jtot,
                                       rel.n12,rel.l12,rel.s12, rel.j12, rel.n3, rel.l3, rel.I3, rel.j3, rel.ncm, rel.lcm);
//         double mytc_reorder = AngMom::Tcoeff( lab.na,lab.la,lab.ja, lab.nb,lab.lb,lab.jb, lab.nc,lab.lc,lab.jc, lab.bj12, jtot,
         double mytc_reorder = AngMom::Tcoeff_reorder( lab.na,lab.la,lab.ja, lab.nb,lab.lb,lab.jb, lab.nc,lab.lc,lab.jc, lab.bj12, jtot,
                                       rel.n12,rel.l12,rel.s12, rel.j12, rel.n3, rel.l3, rel.I3, rel.j3, rel.ncm, rel.lcm);
         std::string good1 = std::abs(mytc - tc)<1e-6 ? "good" : "BADval";
         if (good1=="BADval" and  std::abs(mytc + tc)<1e-6 )  good1 = "BADphase";
         std::string good2 = std::abs(mytc - mytc_bf)<1e-6 ? "good" : "BAD";
         std::string good3 = std::abs(mytc - mytc_fancy)<1e-6 ? "good" : "BADfancy";
         std::string good4 = std::abs(mytc - mytc_reorder)<1e-6 ? "good" : "BADreorder";
//   double Tcoeff( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm);
         std::cout << "( " << lab.na << " " << lab.la << " " << lab.ja << " " <<  lab.nb << " " << lab.lb << " " << lab.jb << " " <<  lab.nc << " " << lab.lc << " " << lab.jc << " " <<  lab.bj12 << " " <<  jtot << " "
                           <<  rel.n12 << " " << rel.l12 << " " << rel.s12 << " " <<  rel.j12 << " " <<  rel.n3 << " " <<  rel.l3 << " " <<  rel.I3 << " " <<  rel.j3 << " " <<  rel.ncm << " " <<  rel.lcm << " ) " << std::endl;
         std::cout << "   " << isp << " " << irel << "  " << index << " ( " << tcoeff.size() << " )  -> "
                   << std::fixed << std::setw(12) << tc << "   "
                   << std::fixed << std::setw(12) << mytc << "   "
                   << std::fixed << std::setw(12) << mytc_bf << "  "
                   << std::fixed << std::setw(12) << mytc_fancy << "  "
                   << std::fixed << std::setw(12) << mytc_reorder << "  "
                   << std::fixed << std::setw(8) << good1 
                   << std::fixed << std::setw(8) << good2
                   << std::fixed << std::setw(10) << good3
                   << std::fixed << std::setw(12) << good4
                   << std::endl;
       }
     }
   }
  }

  std::cout << "All done. Exiting nicely. " << std::endl;

}



//  double Jacobi3BME::ComputeTcoeff( HartreeFock& hf, int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm)
  double Jacobi3BME::ComputeTcoeff( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm)
  {
    double ja = 0.5*j2a;
    double sa = 0.5;
    double jb = 0.5*j2b;
    double sb = 0.5;
    double jc = 0.5*j2c;
    double sc = 0.5;
//
    double S2 = 0.5; 
//    std::cout << std::endl << " " << __func__ << " ( " << na << " " << la << " " << j2a << "  " << nb << " " << lb << " " << j2b << "  " << nc << " " << lc << " " << j2c
//              << "   " << Jab << " " << twoJ << "  " << N1 << " " << L1 << " " << S1 << "   " << N2 << " " << L2 << "  " << twoJ2 << " " << twoJ12 << "  " << Ncm << " " << Lcm << " ) " << std::endl;

    // Limits
    int Lab_min = std::max( std::abs(la-lb), std::abs(Jab-S1) );
    int Lab_max = std::min( la+lb , Jab+S1 );
    int twoS12_min = std::abs(2*S1-1);
    int twoS12_max = 2*S1+1;
    if ( Lab_max<Lab_min or  twoS12_max<twoS12_min ) return 0.0;
//    std::cout << "   Lab min/max = " << Lab_min << " " << Lab_max << std::endl;

    if ( std::abs(twoJ-twoJ12)>2*Lcm or (twoJ+twoJ12)<Lcm) return 0.0;

    double tcoeff=0.;
    for (int Lab=Lab_min; Lab<=Lab_max; Lab++)
    {
//      double ninejab_right = AngMom::phase(Lab) * (2*Lab+1) * AngMom::NineJ(la, lb, Lab, sa, sb, S1, ja, jb, Jab ); 
      double ninejab = AngMom::phase(Lab) * (2*Lab+1) * GetNineJ(2*la, 2*lb, 2*Lab, 1, 1, 2*S1, j2a, j2b, 2*Jab ); 
      if (std::abs(ninejab)<1e-8) continue; 

      // Moshinsky's are expensive, so let's do them in the outer loop (maybe even more outer than this??)
      int curlyL_min =  std::abs(L1-Lab) ;
      int curlyL_max =  L1+Lab ;
      if ( (curlyL_min + L1 + la + lb)%2>0 ) curlyL_min++; // parity for the first Moshinsky bracket
      if ( (curlyL_min + lc + Lcm + L2)%2>0 ) continue; // parity for the second Moshinsky bracket
      for (int curlyL=curlyL_min; curlyL<=curlyL_max; curlyL+=2)
      {
        int curlyN = na+nb-N1  +(la+lb-L1-curlyL)/2;
        if (curlyN<0) continue;
        if ( 2*curlyN+curlyL +2*nc+lc != 2*Ncm+Lcm + 2*N2+L2) continue; // energy conservation in second Moshinsky bracket
        int mosh1phase = AngMom::phase( curlyL+L1+Lab + 0*(curlyL+L1 - lb-la)/2 );
//        double mosh1 = mosh1phase  * AngMom::Moshinsky( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab,  1.0);
//        double mosh1 = mosh1phase  * hf.modelspace->GetMoshinsky( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab); // we've cached the d=1 brackets
        double mosh1 = mosh1phase  * GetMoshinsky1( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab); // we've cached the d=1 brackets
        if (std::abs(mosh1)<1e-9) continue;

        int Lambda_min = std::max( std::abs(Lcm-L2), std::abs(curlyL-lc) );
        int Lambda_max = std::min( Lcm+L2, curlyL+lc );
        for (int Lambda=Lambda_min; Lambda<=Lambda_max; Lambda++) 
        {
          // It will speed things up a bit to precompute these
          int moshphase2 = AngMom::phase( Lcm+L2+Lambda + 0*(Lcm+L2 - curlyL-lc)/2 ); // converting phase conventions (comment at the beginning of the function)
          double mosh2 = moshphase2 * (2*Lambda+1) * GetMoshinsky2(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda);
//          double mosh2 = moshphase2 * (2*Lambda+1) * AngMom::Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  2.0);
          if (std::abs(mosh2)<1e-9) continue;

          // This inner loop is all just summing over sixj's and a ninej, and can probably be re-expressed with some fancy recoupling
          double sum_L = 0;

//   Maybe try inserting the 12j stuff here instead of doin the loops
//          int phase_L = AngMom::phase( lc + Lab + (twoS12+twoJ+1-twoJ2)/2 + S1 + Jab);
//          double sixjL = 0

          int L_min = std::max( std::abs(Lambda-L1), std::abs(Lab-lc) );
          int L_max = std::min( Lambda+L1, Lab+lc );

          for (int L=std::max(std::abs(Lambda-L1),L_min); L<=std::min(L_max,Lambda+L1); L++)
          {
//            double sixj1 =  AngMom::phase(L+Lambda) * (2*L+1) * AngMom::SixJ( lc, curlyL, Lambda, L1,L,Lab);
            // { lc curlyL Lambda } Triangles: (lc,curlyL,Lambda), (L1,L,Lambda), (lc,L,Lab), (L1,curlyL,Lab)
            // { L1   L     Lab   }
            double sixj1 =  AngMom::phase(L+Lambda) * (2*L+1) * GetSixJ( 2*lc, 2*curlyL, 2*Lambda, 2*L1,2*L,2*Lab);
            if (std::abs(sixj1)<1e-9) continue;

            double sum_12 = 0;
            for (int twoS12=twoS12_min; twoS12<=twoS12_max; twoS12+=2)
            {
              if ( std::abs(2*L-twoS12)>twoJ or 2*L+twoS12<twoJ ) continue; // triangle condition for ninejL
              double ninejL =  GetNineJ( 2*lc, 2*L, 2*Lab,
                                         1, twoS12, 2*S1,
                                         j2c, twoJ, 2*Jab);
//              double ninejL =  AngMom::NineJ( Lab, lc, L,
//                                              S1,  sc, 0.5*twoS12,
//                                              Jab, jc, 0.5*twoJ);

              if ( std::abs(ninejL)<1e-9) continue;
              int L12_min = std::max( std::abs(L1-L2), std::max( std::abs(twoS12-twoJ12)/2, std::abs(Lcm-L)) ) ;
              int L12_max = std::min( L1+L2, std::min( (twoS12+twoJ12)/2, Lcm+L));
              for (int L12=L12_min; L12<=L12_max; L12++)
              {
                if ( std::abs(Lcm-L12)>L or (Lcm+L12)<L ) continue;
                if ( std::abs(twoS12-twoJ)>2*L or (twoS12+twoJ)<2*L ) continue;
                double ninej12 = AngMom::phase((twoS12+twoJ)/2)*(2*L12+1)*(twoS12+1)
                                              * GetNineJ(2*L2,  2*L12,  2*L1,
                                                         1,     twoS12, 2*S1,
                                                         twoJ2, twoJ12, 2*J1);
                if (std::abs(ninej12)<1e-9) continue;
//                double sixj2 = AngMom::SixJ( Lcm,L12,L,0.5*twoS12,0.5*twoJ,0.5*twoJ12);
                // { Lcm  L12  L   }  Triangles : (Lcm,L12,L),  (S12,J,L),  (Lcm,J,J12),  (S12,L12,J12)
                // { S12  J    J12 }
                double sixj2 = GetSixJ( 2*Lcm,2*L12,2*L, twoS12,twoJ,twoJ12);
//                double sixj3 = AngMom::SixJ( Lcm,L2,Lambda,L1,L,L12);
                double sixj3 = GetSixJ( 2*Lcm,2*L2,2*Lambda,2*L1,2*L,2*L12);

                sum_12 += ninejL * ninej12 * sixj2 * sixj3;
              } // for L12
            } // for twoS12
            sum_L += sixj1 * sum_12;
          } // for L
          tcoeff += ninejab * mosh1 * mosh2 * sum_L;
        } // for Lambda
      } // for curlyL
    } // for Lab

  // multiply by all the j-hat symbols.   
  // accounting for all the hat factors: 
  //  We also have (2S12+1)(2L12+1) included in ninej12,  (2L+1) included in sixj1,  (2Lab+1) included in ninejab, (2Lambda+1) included in mosh2, and (2curlyL+1) included in mosh1
    // hat factors:  sqrt( (2*ja+1)*(2*jb+1)*(2*jc+1)*(2*Jab+1)*(twoJ12+1)*(2*J1+1)*(twoJ2+1)*(2*S1+1) ) *   (2*Lab+1) * (twoS12+1) (2*L12+1) * (2*L+1) * (2*Lambda+1) 
    //                                                 global                                                (ninejab) * (    ninej12       ) * (sixj1)   (   mosh2  )
    // phases:    ( L12 + S12 + Lcm + J + L1 + L2 - L12 + Lcm + L2 + L1 + L  + curlyL + lc - Lambda + lc + Lab - L + lc + L1 + curlyL + L )
    // reduces to (       S12       + J                                               + lc - Lambda      + Lab - L      + L1              )
    //            (L1 + lc) * ( S12+J ) * (Lab) * (Lambda +L)
    //              global     ninej12   ninejab     sixj1
    //
    // phase( lc + Lab + L + L1 + (twoJ + twoS12)/2 ) * phase( Lambda)
    double jhats = sqrt( (2*ja+1)*(2*jb+1)*(2*jc+1)*(2*Jab+1)*(twoJ12+1)*(2*J1+1)*(twoJ2+1)*(2*S1+1) ) ;
    int globalphase = AngMom::phase(L1+lc);
    return tcoeff * jhats * globalphase;
  }















uint64_t Jacobi3BME::SixJHash(int j1, int j2, int j3, int J1, int J2, int J3)
{
   return   (((uint64_t)(j1)) << 50)
          + (((uint64_t)(j2)) << 40)
          + (((uint64_t)(j3)) << 30)
          + (((uint64_t)(J1)) << 20)
          + (((uint64_t)(J2)) << 10)
          +  ((uint64_t)(J3));
}

void Jacobi3BME::SixJUnHash(uint64_t key, uint64_t& j1, uint64_t& j2, uint64_t& j3, uint64_t& J1, uint64_t& J2, uint64_t& J3)
{
   j1 = (key >> 50) & 0x3FFL;
   j2 = (key >> 40) & 0x3FFL;
   j3 = (key >> 30) & 0x3FFL;
   J1 = (key >> 20) & 0x3FFL;
   J2 = (key >> 10) & 0x3FFL;
   J3 = (key      ) & 0x3FFL;
}

//uint64_t Jacobi3BME::MoshinskyHash(uint64_t N, uint64_t Lam, uint64_t n, uint64_t lam, uint64_t n1, uint64_t l1, uint64_t n2, uint64_t l2, uint64_t L, uint64_t d)
uint64_t Jacobi3BME::MoshinskyHash(uint64_t N, uint64_t Lam, uint64_t n, uint64_t lam, uint64_t n1, uint64_t l1, uint64_t n2, uint64_t l2, uint64_t L)
{
   return   (N   << 54)
          + (Lam << 47)
          + (n   << 41)
          + (lam << 34)
          + (n1  << 28)
          + (l1  << 21)
          + (n2  << 15)
          + (l2  << 8 )
          + (L        );
//          + (L   << 1 )
//          + (d%2);
}

//void Jacobi3BME::MoshinskyUnHash(uint64_t key,uint64_t& N,uint64_t& Lam,uint64_t& n,uint64_t& lam,uint64_t& n1,uint64_t& l1,uint64_t& n2,uint64_t& l2,uint64_t& L,uint64_t d)
void Jacobi3BME::MoshinskyUnHash(uint64_t key,uint64_t& N,uint64_t& Lam,uint64_t& n,uint64_t& lam,uint64_t& n1,uint64_t& l1,uint64_t& n2,uint64_t& l2,uint64_t& L)
{
   N   = (key >> 54) & 0x3FL;
   Lam = (key >> 47) & 0x7FL;
   n   = (key >> 41) & 0x3FL;
   lam = (key >> 34) & 0x7FL;
   n1  = (key >> 28) & 0x3FL;
   l1  = (key >> 21) & 0x7FL;
   n2  = (key >> 15) & 0x3FL;
   l2  = (key >> 8 ) & 0x7FL;
   L   = (key      ) & 0xFFL;
//   L   = (key >> 1 ) & 0x7FL;
//   d = 2-(key%2);
}

size_t Jacobi3BME::NineJHash( int twol1, int twol2, int twol3, int twos1, int twos2, int twos3, int twoj1, int twoj2, int twoj3)
{
//  return 6*(twoj3-twol3+3*twos3-2) + (twoj2-twol2+3*twos2-2) + (twoj1-twol1+3*twos1-2)/2;
  return 48*( (twol1/2*(Lmax_nj+1))*(twol1/2*(Lmax_nj+1)) + twol2/2*(Lmax_nj+1) + twol3)
       + 6*(twoj3-twol3+2*twos3) + (twoj2-twol2+3*twos2-2) + (twoj1-twol1+3*twos1-2)/2;
}


double Jacobi3BME::GetSixJ(int j1, int j2, int j3, int J1, int J2, int J3)
{
// { j1 j2 j3 }
// { J1 J2 J3 }
//   std::cout << " in " << __func__ << "  " << j1 << " " << j2 << " "<< j3 << "    " << J1 << " " << J2 << " " << J3 << std::endl;
   uint64_t key = SixJHash(j1,j2,j3,J1,J2,J3);

   const auto it = SixJList.find(key);
   double sixj=0.0;
   if (it != SixJList.end() )
   {
     sixj = it->second;
   }
   else
   {
    sixj = AngMom::SixJ(0.5*j1,0.5*j2,0.5*j3,0.5*J1,0.5*J2,0.5*J3);
//    if (omp_get_num_threads()<0)
    if (omp_get_num_threads()<2)
    {
      #pragma omp critical
      {
        SixJList[key] = sixj;
      }
    }
//    else
//    {
//      std::cout << "DANGER!!!!!!!  Updating SixJList inside a parellel loop breaks thread safety!" << std::endl;
//      std::cout << "  I shouldn't be here in GetSixJ("
//                << std::setprecision(1) << std::fixed << j1 << " " << std::setprecision(1) << std::fixed << j2 << " "
//                << std::setprecision(1) << std::fixed << j3 << " " << std::setprecision(1) << std::fixed << J1 << " "
//                << std::setprecision(1) << std::fixed << J2 << " " << std::setprecision(1) << std::fixed << J3 << "). key = "
//                << std::hex << key << "   sixj = " << std::dec << sixj << std::endl;
//      IMSRGProfiler::counter["N_CalcSixJ_in_Parallel_loop"] +=1;
//      exit(EXIT_FAILURE);
//    }
   }
   return sixj;
}

double Jacobi3BME::GetMoshinsky1( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L)
{
  int phase_mosh = 1;
  int switches = 10;
//  std::cout << "Enter GetMoshinsky " << N << " " << Lam << " " << n << " " << lam << " " << n1 << " " << l1 << " " << n2 << " " << l2 << " " << L << " " << d << std::endl;

  while (switches > 0)
  {
   switches = 0;
   if (n2>n1 or (n2==n1 and l2>l1))
   {
      std::swap(n1,n2);
      std::swap(l1,l2);
      phase_mosh *= AngMom::phase(Lam+L);
      ++switches;
   }
   if (n>N or (n==N and lam>Lam))
   {
      std::swap(n,N);
      std::swap(lam,Lam);
      phase_mosh *= AngMom::phase(l1 +L);
      ++switches;
   }

//   if (l1>Lam or (l1==Lam and n1>N) or (l1==Lam and n1==N and l2>lam) or (l1==Lam and n1==N and l2==lam and n2>n) )
   if (n1>N or (n1==N and l1>Lam) or (n1==N and l1==Lam and n2>n) or (n1==N and l1==Lam and n2==n and l2>lam) )
   {
      std::swap(n1,N);
      std::swap(l1,Lam);
      std::swap(n2,n);
      std::swap(l2,lam);
      ++switches;
//      phase_mosh *= phase(l2+lam); // This phase is given in Moshinsky and Brody, but with the current algorithm, it appears not to be required.
   }
  }

//   std::cout << "    after swapping, " <<  N << " " << Lam << " " << n << " " << lam << " " << n1 << " " << l1 << " " << n2 << " " << l2 << " " << L << " " << d << std::endl;
//   uint64_t key = MoshinskyHash(N,Lam,n,lam,n1,l1,n2,l2,L,d);
   uint64_t key = MoshinskyHash(N,Lam,n,lam,n1,l1,n2,l2,L);

   auto it = Moshinsky1List.find(key);

   if ( it != Moshinsky1List.end() )  return it->second * phase_mosh;
//   if (omp_get_num_threads()>1)
//   {
//     std::cout << "TROUBLE IN MOSHINSKY LAND!!!!!    <" << N << " " << Lam << " " << n << " " << lam << " | " << n1 << " " << l1 << " " << n2 << " " << l2 << ">_" << L  << " d = 1" << std::endl;
//   }

   // if we didn't find it, we need to calculate it.
   double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L,1);
//   #pragma omp atomic
   if (omp_get_num_threads()<2)
   {
     Moshinsky1List[key] = mosh;
   }
   return mosh * phase_mosh;

}


// The mass factor screws up the symmetries, so there's less to worry about here.
double Jacobi3BME::GetMoshinsky2( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L)
{
   uint64_t key = MoshinskyHash(N,Lam,n,lam,n1,l1,n2,l2,L);
   auto it = Moshinsky2List.find(key);
   if ( it != Moshinsky2List.end() )  return it->second ;
//   if (omp_get_num_threads()>1)
//   {
//     std::cout << "TROUBLE IN MOSHINSKY LAND!!!!!    <" << N << " " << Lam << " " << n << " " << lam << " | " << n1 << " " << l1 << " " << n2 << " " << l2 << ">_" << L  << " d = 2" << std::endl;
//   }
   // if we didn't find it, we need to calculate it.
   double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L,2);
//   #pragma omp atomic
   if (omp_get_num_threads()<2)
   {
     Moshinsky2List[key] = mosh;
   }
   return mosh;
}


double Jacobi3BME::GetNineJ( int twol1, int twol2, int twol3, int twos1, int twos2, int twos3, int twoj1, int twoj2, int twoj3)
{
//  std::cout << " GetNineJ: " << twol1 << " " << twol2 << " " << twol3 << " "
//                             << twos1 << " " << twos2 << " " << twos3 << " "
//                             << twoj1 << " " << twoj2 << " " << twoj3 << std::endl;
//  std::cout << " NineJHash = " << NineJHash( twol1, twol2, twol3, twos1, twos2, twos3, twoj1, twoj2, twoj3) << std::endl;
  return NineJList.at( NineJHash( twol1, twol2, twol3, twos1, twos2, twos3, twoj1, twoj2, twoj3) );
}






//  This is structured in essentially the same way as the analogous function in ModelSpace
//  We need { j1 j2 j3 }  for sixj1 : ( lc,  curlyL, Lambda,   L1,  L, Lab ) -- all integer,  0 < j < 3*max(lc)
//          { J4 J5 J6 }      sixj2 : ( Lcm, L12,    L,        S12, J, J12 )    -- top row integer, bottom row half-integer
//                            sixj3 : ( Lcm, L2,     Lambda,   L1,  L ,L12 )   -- all integer, with the same restrictions as sixj1
// 
// Limit on Lcm : must be <= E3max from energy cons.
// Limit on curlyL : must be <= E2max from energy cons.
// Limit on Lambda : Triangle.
// Limit on L2 : must be <= Nmax
// Limit on L : 3 * lmax.
  void Jacobi3BME::PreComputeSixJ()
  {
    double t_start = omp_get_wtime();
    std::cout << "Precalculating SixJ's" << std::endl;
    std::vector<uint64_t> KEYS;
    // first, we do the all-integer ones
    for (int j1=0; j1<=(6*emax+1); j1+=2)  // 2 * Lcm can go up to 6emax
    {
     for (int j2=0; j2<=(6*emax+1); j2+=2)
     {
      for (int j3=std::abs(j1-j2); j3<=(j1+j2); j3+=2)
      {
//       for (int J4=0; J4<=3*(2*emax+1); J4+=2)
       for (int J4=0; J4<=(2*Nmax+1); J4+=2)
       {
        for (int J5=std::abs(j3-J4); J5<=(j3+J4); J5+=2)
        {
         for (int J6=std::max(std::abs(j1-J5),std::abs(j2-J4)); J6<=std::min(j1+J5,j2+J4); J6+=2)
         {
           if ( (j1 + J6)>2*std::max(E3max,Nmax)) continue;
           uint64_t key = Jacobi3BME::SixJHash(j1,j2,j3,J4,J5,J6);
           if ( SixJList.count(key) == 0 ) 
           {
             KEYS.push_back(key);
             SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
           }
         } // for J6
        } // for J5
       } // for J4

       // Now we do the half-integer bottom row
       //       double sixj2 = GetSixJ( 2*Lcm,2*L12,2*L, twoS12,twoJ,twoJ12); twoJ <= 6emax+3, twoJ12<=Jmax, L could be up to 6emax+6
       // { j1  j2 j3  }
       // { S12 J  J12 }
       if ((j1+j2)>2*E3max) continue;
       for (int S12=1; S12<=3; S12+=2)
       {
        for (int J=std::abs(S12-j3); J<=(S12+j3); J+=2)
        {
         for (int J12=std::max(std::abs(j1-J),std::abs(j2-S12)); J12<=std::min((j1+J),(j2+S12)); J12+=2)
         {
           uint64_t key = Jacobi3BME::SixJHash(j1,j2,j3,S12,J,J12);
           if ( SixJList.count(key) == 0 ) 
           {
             KEYS.push_back(key);
             SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
           }
         } // for J12
        } // for J
       } // for S12
      } // for j3
     } // for j2
    } // for j1
    
    // now we actually compute them all.
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t i=0;i< KEYS.size(); ++i)
    {
      uint64_t j1,j2,j3,J1,J2,J3;
      uint64_t key = KEYS[i];
      SixJUnHash(key, j1,j2,j3,J1,J2,J3);
      SixJList[key] = AngMom::SixJ(0.5*j1,0.5*j2,0.5*j3,0.5*J1,0.5*J2,0.5*J3);
    }
//    sixj_has_been_precalculated = true;
    std::cout << "done calculating sixJs (" << KEYS.size() << " of them)" << std::endl;
    std::cout << "Hash table has " << SixJList.bucket_count() << " buckets and a load factor " << SixJList.load_factor() 
         << "  estimated storage ~ " << ((SixJList.bucket_count()+SixJList.size()) * (sizeof(size_t)+sizeof(void*))) / (1024.*1024.*1024.) << " GB" << std::endl;
    IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  }

  // All the relevant 9js have the structure
  // { l1 l2 l3 }  
  // { s1 s2 s3 }  <- s1=1/2,  s2=1/2,3/2,  s3=0,1
  // { j1 j2 j3 }
  //
  void Jacobi3BME::PreComputeNineJ()
  {
    double t_start = omp_get_wtime();
    Lmax_nj = std::max(Nmax,E3max);
    std::cout <<"Calculating NineJ with Lmax = " << Lmax_nj << std::endl;
//    NineJList.resize(48*Lmax_nj*Lmax_nj*Lmax_nj);
    NineJList.resize( NineJHash( 2*Lmax_nj, 2*Lmax_nj, 2*Lmax_nj, 1, 3, 2, 2*Lmax_nj+1, 2*Lmax_nj+3, 2*Lmax_nj+2) +1  );

    #pragma omp parallel for schedule(dynamic,1)
    for ( int l1=0; l1<=Lmax_nj; l1++)
    {
     int twoj1_min=std::max(2*l1-1,1);
     int twoj1_max=2*l1+1;
     for (int l2=0; l2<=Lmax_nj; l2++)
     {
      int l3min= std::abs(l1-l2);
      int l3max = std::min(l1+l2,Lmax_nj);
      for (int l3=l3min; l3<=l3max; l3++)
      {
       for (int twos2=1; twos2<=3; twos2+=2)
       {
        int twoj2_min=std::max(1,2*l2-twos2);
        int twoj2_max=2*l2+twos2;
        int s3_min = std::abs(twos2-1)/2;
        int s3_max = 1;
        for (int s3=s3_min; s3<=s3_max; s3++)
        {
         for (int twoj1=twoj1_min; twoj1<=twoj1_max; twoj1+=2)
         {
          for (int twoj2=twoj2_min; twoj2<=twoj2_max; twoj2+=2)
          {
           int j3min = std::max( std::abs(l3-s3), std::abs(twoj1-twoj2)/2);
           int j3max = std::min( l3+s3, (twoj1+twoj2)/2 );
           for (int j3=j3min; j3<=j3max; j3++)
           {
             size_t hash = NineJHash(2*l1,2*l2,2*l3, 1, twos2, 2*s3, twoj1, twoj2, 2*j3);
             if (hash > NineJList.size() )
             {
               std::cout << __func__ << " Computing:: hash = " << hash << "  which is > " << NineJList.size()
                         << "   :  " << l1 << " " << l2 << " " << l3
                         << "  "     << 1  << " " << twos2 << " " << s3
                         << "  "     << twoj1 << " "<< twoj2 << " " << j3 << std::endl;
//               std::cout << "      " << 6*(2*j3-2*l3+3*2*s3-2)
//               std::cout << "      " << 6*(2*j3-2*l3+2*2*s3)
//                         << "   " << (twoj2-2*l2+3*twos2-2)
//                         << "   " << (twoj1-2*l1+3*1-2)/2 << std::endl;

//  return 6*(twoj3-twol3+3*twos3-2) + (twoj2-twol2+3*twos2-2) + (twoj1-twol1+3*twos1-2)/2;
             }
             NineJList[hash] = AngMom::NineJ( l1,l2,l3, 0.5,0.5*twos2,s3, 0.5*twoj1, 0.5*twoj2, j3);
//             std::cout << "     -> Ninej = " << NineJList[hash] << std::endl;
           }
          }
         }
        }
       }
      }
     }
    }
    std::cout << "done calculating ninej (" << NineJList.size() << " elements)"
         << "  estimated storage ~ " << ((NineJList.size()) * (sizeof(size_t)+sizeof(void*))) / (1024.*1024.*1024.) << " GB" << std::endl;
    IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  }


//        double mosh1 = mosh1phase  * hf.modelspace->ambdaGetMoshinsky( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab); // we've cached the d=1 brackets
//          double mosh2 = moshphase2 * (2*Lambda+1) * AngMom::Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  2.0);

  // Two kinds of Moshinsky bracket to be computed
  // mosh1 = << curlyN,curlyL, N1,L1 |  na,la, nb,lb; Lab >>_1
  // The first one has E12 <= E2max;
  // The second one has a more complicated restriction. Ecurly <=  Eab - E1, and E1+E2<= Nmax, and Ecm+Nmax <= E3max
  // so Ecm + E2 <= E3max ?
  void Jacobi3BME::PreComputeMoshinsky1()
  {
    double t_start = omp_get_wtime();
    std::cout <<"Calculating moshinsky1 with Lmax = " << emax << std::endl;
    int Lmax = emax;
  
    // generating all the keys is fast, so we do this first without parallelization
  //  std::vector<unsigned long long int> KEYS;
    std::vector<uint64_t> KEYS;
    for (int N=0; N<=E2max/2; ++N)
    {
     for (int n=0; n<=std::min(N,E3max/2-N); ++n)
     {
//      std::cout << "N,n = " << N << " " << n << "   n is limited to " << std::min(N,E3max/2-N) << std::endl;
      int Lam_max = std::min( E2max-2*n-2*n, 2*Lmax ); // Here Lmax is the max L of the s.p. basis
      for (int Lam=0; Lam<=E2max-2*N-2*n; ++Lam)
      {
       int lam_max = N==n ? std::min(Lam,E2max-2*N-2*n-Lam) : E2max-2*N-2*n-Lam ; 
       for (int lam=0; lam<=lam_max; ++lam)
       {
        int e2 = 2*N+Lam + 2*n+lam;
        for (int L=std::abs(Lam-lam); L<=Lam+lam; ++L)
        {
         if (L>3*Lmax) continue;
         for (int n1=0; n1<=N; ++n1)
         {
          for (int n2=0; n2<=std::min(n1,e2/2-n1); ++n2)
          {
           int l1max = n1==N? std::min(Lam,e2-2*n1-2*n2) : e2-2*n1-2*n2;
           for (int l1=0; l1<=l1max; ++l1 )
           {
            int l2 = e2-2*n1-2*n2-l1;
            if ( (l1+l2+lam+Lam)%2 >0 ) continue; // parity conservation
            if ( l2<std::abs(L-l1) or l2>L+l1 ) continue; // triangle
  //          if (l1>Lmax or l2>Lmax) continue;
            if ( (l1>Lmax and l2>Lmax) and (lam>Lmax or Lam>Lmax)) continue; // maybe we don't want this?
            if ( (l1>Lmax or l2>Lmax) and (lam>Lmax and Lam>Lmax)) continue;
  
//            std::cout << "::: " << N << " " << Lam << " " << n << " " << lam << " " << n1 << " " << l1 << " " << n2 << " " << l2 << " " << L << std::endl;

            uint64_t key = MoshinskyHash(N,Lam,n,lam,n1,l1,n2,l2,L);
            KEYS.push_back(key);
            Moshinsky1List[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop

//            key = MoshinskyHash(N,Lam,n,lam,n1,l1,n2,l2,L,2);
//            KEYS.push_back(key);
//            MoshinskyList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
           }
          }
         }
        }
       }
      }
     }
    }
    // Now we calculate the Moshinsky brackets in parallel
  //  std::vector<double> mosh_vals( KEYS.size() );
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t i=0;i< KEYS.size(); ++i)
    {
      uint64_t key = KEYS[i];
//      uint64_t N,Lam,n,lam,n1,l1,n2,l2,L,d;
      uint64_t N,Lam,n,lam,n1,l1,n2,l2,L;
      MoshinskyUnHash(key,N,Lam,n,lam,n1,l1,n2,l2,L);
  
      Moshinsky1List[key] = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L,1);
    }
  
//    moshinsky_has_been_precalculated = true;
    std::cout << "done calculating moshinsky (" << KEYS.size() << " elements)" << std::endl;
    std::cout << "Hash table has " << Moshinsky1List.bucket_count() << " buckets and a load factor " << Moshinsky1List.load_factor() 
              << "  estimated storage ~ " << ((Moshinsky1List.bucket_count()+Moshinsky1List.size()) * (sizeof(size_t)+sizeof(void*))) / (1024.*1024.*1024.) << " GB" << std::endl;
    IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  }



 // mosh2 = << Ncm,Lcm,  N2,L2  |  curlyN,curlyL, nc,lc; Lambda >>_2
 //
 // ranges nc,lc given by emax.
 // Ncm,Lcm <= E3max or 3*emax
 // N2,L2 <= Nmax
 // curlyN,curlyL given by energy conservation.
 void Jacobi3BME::PreComputeMoshinsky2()
  {
    double t_start = omp_get_wtime();
    std::cout <<"Calculating moshinsky2 with E3max = " << E3max << std::endl;
    int Lmax = emax;
  
    // generating all the keys is fast, so we do this first without parallelization
  //  std::vector<unsigned long long int> KEYS;
    std::vector<uint64_t> KEYS;
    for (int Ncm=0; Ncm<=E3max/2; Ncm++)
    {
     for (int N2=0; N2<=Nmax/2; N2++)
     {
//      for (int Lcm=0; 2*Ncm+Lcm<=E3max/2; Lcm++) // this isn't all of them, but it (maybe?) keeps the number of terms manageable
      for (int Lcm=0; 2*Ncm+Lcm<=E3max; Lcm++) // this isn't all of them, but it (maybe?) keeps the number of terms manageable
      {
       for (int L2=0; 2*N2+L2<=Nmax; L2++)
       {
        for (int nc=0; nc<=emax/2; nc++)
        {
         for (int lc=0; 2*nc+lc<=emax; lc++)
         {
          for (int curlyN=0; 2*curlyN + 2*nc+lc <= 2*(Ncm+N2)+Lcm+L2; curlyN++)
          {
           int curlyL = 2*(Ncm+N2-nc-curlyN)+Lcm+L2-lc;
           if (curlyL<0) continue;
           if ( (Lcm+L2+lc+curlyL)%2>0 ) continue; // parity conservation
           int Lambda_min = std::max( std::abs(Lcm-L2), std::abs(lc-curlyL) );
           int Lambda_max = std::min( Lcm+L2, lc+curlyL);
           for (int Lambda=Lambda_min; Lambda<=Lambda_max; Lambda++)
           {
//       std::cout << __func__ << "  precomputing    <" << Ncm << " " << Lcm << " " << N2 << " " << L2 << " | " << curlyN << " " << curlyL << " " << nc << " " << lc << ">_" << Lambda  << " d = 2" << std::endl;
             uint64_t key = MoshinskyHash(Ncm,Lcm,N2,L2,curlyN,curlyL,nc,lc,Lambda);
             KEYS.push_back(key);
             Moshinsky2List[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
           }
          } // for curlyN
         } // for lc
        } // for nc
       } // for L2
      } // for Lcm
     } // for N2
    } // for Ncm
    // Now we calculate the Moshinsky brackets in parallel
  //  std::vector<double> mosh_vals( KEYS.size() );
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t i=0;i< KEYS.size(); ++i)
    {
      uint64_t key = KEYS[i];
//      uint64_t N,Lam,n,lam,n1,l1,n2,l2,L,d;
      uint64_t N,Lam,n,lam,n1,l1,n2,l2,L;
      MoshinskyUnHash(key,N,Lam,n,lam,n1,l1,n2,l2,L);
  
      Moshinsky2List[key] = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L,2);
    }
  
//    moshinsky_has_been_precalculated = true;
    std::cout << "done calculating moshinsky2 (" << KEYS.size() << " elements)" << std::endl;
    std::cout << "Hash table has " << Moshinsky2List.bucket_count() << " buckets and a load factor " << Moshinsky2List.load_factor() 
              << "  estimated storage ~ " << ((Moshinsky2List.bucket_count()+Moshinsky2List.size()) * (sizeof(size_t)+sizeof(void*))) / (1024.*1024.*1024.) << " GB" << std::endl;
    IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  }






