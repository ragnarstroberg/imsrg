
#include "Jacobi3BME.hh"
#include "AngMom.hh"
#include "HartreeFock.hh" // needed for V3monHash and V3monUnHash
#include "IMSRGProfiler.hh"
#include <istream>
#include <set>
#include <unordered_set>
#include <sstream>
#include <iomanip>
#include <omp.h>




Jacobi3BME::Jacobi3BME()
{}

Jacobi3BME::Jacobi3BME( int nmax, int twojmin, int twojmax, int twotmin, int twotmax )
 : Nmax(nmax), twoJmin(twojmin), twoJmax(twojmax), twoTmin(twotmin), twoTmax(twotmax), twoJmax_cut(twojmax), lmax_NO2b(1000)
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
            std::cout << "TROUBLE. I don't understand how the NAS basis states are constructed. Dim = " << NAS_jacobi_states.at(hashtjn).size() << "  should be  " <<  dim_braNAS <<  "    TJN =  " << t <<  " " << j << " " << Nbra << std::endl;
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

  int na=bra.op->n, nb=bra.oq->n, nc=bra.oR->n,  nd=ket.op->n, ne=ket.oq->n, nf=ket.oR->n;
  int j2a=bra.op->j2, j2b=bra.oq->j2, j2c=bra.oR->j2,  j2d=ket.op->j2, j2e=ket.oq->j2, j2f=ket.oR->j2;
  int la=bra.op->l, lb=bra.oq->l, lc=bra.oR->l,  ld=ket.op->l, le=ket.oq->l, lf=ket.oR->l;

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
//      int twoJ12_max = std::min(twoJmax, twoJ + 2*Lcm) ;
      int twoJ12_max = std::min(twoJmax_cut, twoJ + 2*Lcm) ;
//      int twoJ12_min = std::max(twoJmin, std::abs( twoJ - 2*Lcm ));
//      int twoJ12_max = std::min(twoJmax, twoJ + 2*Lcm) ;

//      std::cout << "Ecm,Lcm,Ncm = " << Ecm << " " << Lcm << " " << Ncm << "   twoJ12 runs from " << twoJ12_min << " to " << twoJ12_max << std::endl;
      for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2) // the interaction conserves J12, parity, and T12=twoT
      {
//        std::cout << "Getting dimensions for J12 = " << twoJ12 << "/2" << "   twoJmax = " << twoJmax << std::endl;
        int NASdim_bra = GetDimensionNAS( twoT, twoJ12, parity12_bra, E12_bra ); 
        int NASdim_ket = GetDimensionNAS( twoT, twoJ12, parity12_ket, E12_ket ); 

        size_t ASdim_bra = GetDimensionAS( twoT, twoJ12, parity12_bra, E12_bra ); 
        size_t ASdim_ket = GetDimensionAS( twoT, twoJ12, parity12_ket, E12_ket ); 
//        if (ASdim_bra==0 or ASdim_ket==0) continue;

        size_t cfp_begin_bra = GetCFPStartLocation(twoT,twoJ12,E12_bra);
        size_t cfp_begin_ket = GetCFPStartLocation(twoT,twoJ12,E12_ket);
        size_t startlocAS = GetStartLocAS( twoT, twoJ12, E12_bra, E12_ket);

        arma::mat matelAS( &meAS[startlocAS], ASdim_bra, ASdim_ket, false ); 
        
        arma::mat cfp_bra( &(cfpvec[cfp_begin_bra]), NASdim_bra, ASdim_bra, /*copy_aux_mem*/ true);
        arma::mat cfp_ket( &(cfpvec[cfp_begin_ket]), NASdim_ket, ASdim_ket, /*copy_aux_mem*/ true);


//        std::cout << " ^^^T J12 P E12 = " << twoT << " " << twoJ12 << " " << parity12_bra << " " << E12_bra << "," << E12_ket << "  Ncm,Lcm = " << Ncm << " " << Lcm << std::endl;
//        std::cout << "dimensions: " << NASdim_bra << " " << NASdim_ket << std::endl;

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
//          Tcoeff_bra[ibraNAS] = ComputeTcoeff( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jac1_bra.n, jac1_bra.l, jac1_bra.s, jac1_bra.j, jac2_bra.n, jac2_bra.l, jac2_bra.j2, twoJ12, Ncm, Lcm);
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
//          Tcoeff_ket[iketNAS] = ComputeTcoeff( nd, ld, j2d, ne, le, j2e, nf, lf, j2f, Jde, twoJ, jac1_ket.n, jac1_ket.l, jac1_ket.s, jac1_ket.j, jac2_ket.n, jac2_ket.l, jac2_ket.j2, twoJ12, Ncm, Lcm);
        }


        arma::mat result =  Tcoeff_bra * cfp_bra * matelAS * cfp_ket.t() * Tcoeff_ket;
        me_lab += result[0];


//        std::cout << "matrices:" << std::endl;
//        std::cout << "T bra " << std::endl << Tcoeff_bra << std::endl;
//        std::cout << "mAS : " << std::endl << matelAS << std::endl;
//        std::cout << "mNAS*6: " << std::endl << ( cfp_bra * matelAS * cfp_ket.t() )*6 << std::endl;
//        std::cout << "mAS: " << std::endl << ( matelAS ) << std::endl;
//        std::cout << "T ket " << std::endl << Tcoeff_ket << std::endl;
//        std::cout << "mNAS * Tket " << std::endl << ( cfp_bra * matelAS * cfp_ket.t() )*Tcoeff_ket  *6 << std::endl;
//        std::cout << "result " << result[0]  << " *6 = " << result[0]*6 << "  me_lab " << me_lab << "   *6= " << me_lab*6 << std::endl;


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




// These could probably be written in a much more concise, but totally incomprehensible templated/variadic function

std::array<unsigned short,3> Jacobi3BME::MakeUshort3( const std::array<int,3>& arr)
{
  return {  static_cast<unsigned short>(arr[0]), 
            static_cast<unsigned short>(arr[1]), 
            static_cast<unsigned short>(arr[2])  }; 
}

std::array<unsigned short,4> Jacobi3BME::MakeUshort4( const std::array<int,4>& arr)
{
  return {  static_cast<unsigned short>(arr[0]), 
            static_cast<unsigned short>(arr[1]), 
            static_cast<unsigned short>(arr[2]), 
            static_cast<unsigned short>(arr[3])  }; 
}

std::array<unsigned short,5> Jacobi3BME::MakeUshort5( const std::array<int,5>& arr)
{
  return {  static_cast<unsigned short>(arr[0]), 
            static_cast<unsigned short>(arr[1]), 
            static_cast<unsigned short>(arr[2]), 
            static_cast<unsigned short>(arr[3]), 
            static_cast<unsigned short>(arr[4])  }; 
}

std::array<unsigned short,6> Jacobi3BME::MakeUshort6( const std::array<int,6>& arr)
{
  return {  static_cast<unsigned short>(arr[0]), 
            static_cast<unsigned short>(arr[1]), 
            static_cast<unsigned short>(arr[2]), 
            static_cast<unsigned short>(arr[3]), 
            static_cast<unsigned short>(arr[4]), 
            static_cast<unsigned short>(arr[5])  }; 
}

std::array<unsigned short,7> Jacobi3BME::MakeUshort7( const std::array<int,7>& arr)
{
  return {  static_cast<unsigned short>(arr[0]), 
            static_cast<unsigned short>(arr[1]), 
            static_cast<unsigned short>(arr[2]), 
            static_cast<unsigned short>(arr[3]), 
            static_cast<unsigned short>(arr[4]), 
            static_cast<unsigned short>(arr[5]), 
            static_cast<unsigned short>(arr[6])  }; 
}

std::array<unsigned short,8> Jacobi3BME::MakeUshort8( const std::array<int,8>& arr)
{
  return {  static_cast<unsigned short>(arr[0]), 
            static_cast<unsigned short>(arr[1]), 
            static_cast<unsigned short>(arr[2]), 
            static_cast<unsigned short>(arr[3]), 
            static_cast<unsigned short>(arr[4]), 
            static_cast<unsigned short>(arr[5]), 
            static_cast<unsigned short>(arr[6]), 
            static_cast<unsigned short>(arr[7])  }; 
}

std::array<unsigned short,11> Jacobi3BME::MakeUshort11( const std::array<int,11>& arr)
{
  return {  static_cast<unsigned short>(arr[0]), 
            static_cast<unsigned short>(arr[1]), 
            static_cast<unsigned short>(arr[2]), 
            static_cast<unsigned short>(arr[3]), 
            static_cast<unsigned short>(arr[4]), 
            static_cast<unsigned short>(arr[5]), 
            static_cast<unsigned short>(arr[6]), 
            static_cast<unsigned short>(arr[7]), 
            static_cast<unsigned short>(arr[8]), 
            static_cast<unsigned short>(arr[9]), 
            static_cast<unsigned short>(arr[10])  }; 
}



size_t Jacobi3BME::array3_hash::operator() (const std::array<unsigned short,3>& key) const
 {
   return  ( static_cast<unsigned long>(key[0])       )
          +( static_cast<unsigned long>(key[1]) <<  8 )
          +( static_cast<unsigned long>(key[2]) << 16 );
 }

size_t Jacobi3BME::array4_hash::operator() (const std::array<unsigned short,4>& key) const
 {
   return  ( static_cast<unsigned long>(key[0])       )
          +( static_cast<unsigned long>(key[1]) <<  8 )
          +( static_cast<unsigned long>(key[2]) << 16 )
          +( static_cast<unsigned long>(key[3]) << 24 );
 }

size_t Jacobi3BME::array5_hash::operator() (const std::array<unsigned short,5>& key) const
 {
   return  ( static_cast<unsigned long>(key[0])       )
          +( static_cast<unsigned long>(key[1]) <<  8 )
          +( static_cast<unsigned long>(key[2]) << 16 )
          +( static_cast<unsigned long>(key[3]) << 24 )
          +( static_cast<unsigned long>(key[4]) << 32 );
 }

size_t Jacobi3BME::array6_hash::operator() (const std::array<unsigned short,6>& key) const
 {
   return  ( static_cast<unsigned long>(key[0])       )
          +( static_cast<unsigned long>(key[1]) <<  8 )
          +( static_cast<unsigned long>(key[2]) << 16 )
          +( static_cast<unsigned long>(key[3]) << 24 )
          +( static_cast<unsigned long>(key[4]) << 32 )
          +( static_cast<unsigned long>(key[5]) << 40 );
 }

size_t Jacobi3BME::array7_hash::operator() (const std::array<unsigned short,7>& key) const
 {
   return  ( static_cast<unsigned long>(key[0])       )
          +( static_cast<unsigned long>(key[1]) <<  8 )
          +( static_cast<unsigned long>(key[2]) << 16 )
          +( static_cast<unsigned long>(key[3]) << 24 )
          +( static_cast<unsigned long>(key[4]) << 32 )
          +( static_cast<unsigned long>(key[5]) << 40 )
          +( static_cast<unsigned long>(key[6]) << 48 );
 }

size_t Jacobi3BME::array8_hash::operator() (const std::array<unsigned short,8>& key) const
 {
   return  ( static_cast<unsigned long>(key[0])       )
          +( static_cast<unsigned long>(key[1]) <<  8 )
          +( static_cast<unsigned long>(key[2]) << 16 )
          +( static_cast<unsigned long>(key[3]) << 24 )
          +( static_cast<unsigned long>(key[4]) << 32 )
          +( static_cast<unsigned long>(key[5]) << 40 )
          +( static_cast<unsigned long>(key[6]) << 48 )
          +( static_cast<unsigned long>(key[7]) << 56 );
 }

//  Our 12j is of the form  
//  { J     jc     1/2    J2      }
//  {   Jab      lc     L     J12 }
//  { J1    curlyL Lambda  Lcm    }
size_t Jacobi3BME::array11_hash::operator() (const std::array<unsigned short,11>& key) const
 {
   return  ( static_cast<unsigned long>(key[0])       )
          +( static_cast<unsigned long>(key[1]+1-key[4]) <<  6 ) // jc and lc differ only by 1/2, so we can compress that information
          +( static_cast<unsigned long>(key[2]) << 8 )
          +( static_cast<unsigned long>(key[3]) << 14 )
          +( static_cast<unsigned long>(key[4]) << 20 )
          +( static_cast<unsigned long>(key[5]) << 26 )
          +( static_cast<unsigned long>(key[6]) << 32 )
          +( static_cast<unsigned long>(key[7]) << 38 )
          +( static_cast<unsigned long>(key[8]) << 44 )
          +( static_cast<unsigned long>(key[9]) << 50 )
          +( static_cast<unsigned long>(key[10]) << 56 );
 }


// na,nb,and nc all should be <= emax/2, Jab will be <= 2*emax+1, twoJ will be <=6*emax+3, but only odd, Lcm will be <=3*emax, jac2 <= ~Nmax*Nmax,    jac1<= ~Nmax*Nmax*Nmax 
//  J lies in the range |Jab-jc| <= J <= Jab+jc,  so it has a span of no more than 2*emax+1.
// na + 10*nb + 100*nc fits in 10 bits (1023).  Jab fits in 6 bits (63).  twoJ/2 fits in 7 bits (127).  Lcm fits in 7 bits (127). jac2 fits in 12 bits (4095). jac1 fits in 18 bits (262143)
// This totals up to 10+6+7+7+7+12+18 = 67, which is too many.   For now, we will just use a string.
//uint64_t Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ, uint64_t jac1, uint64_t jac2, uint64_t twoJ12, uint64_t Lcm )
//std::string Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ, uint64_t jac1, uint64_t jac2, uint64_t twoJ12, uint64_t Lcm )
//std::array<unsigned short,7> Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ,  uint64_t twoJ12, uint64_t E12 )
//std::array<unsigned short,7> Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t twoJ,  uint64_t twoJ12, uint64_t E12, uint64_t Lcm )
//std::array<unsigned short,8> Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t twoJ,  uint64_t twoJ12, uint64_t twoT, uint64_t E12, uint64_t Lcm )
std::array<unsigned short,8> Jacobi3BME::TcoeffHash(unsigned short na, unsigned short nb, unsigned short nc, unsigned short twoJ,  unsigned short twoJ12, unsigned short twoT, unsigned short E12, unsigned short Lcm )
{
//  return {na,nb,nc,Jab,twoJ,twoJ12,E12};
//  return {na,nb,nc,twoJ,twoJ12,E12,Lcm};
  return {na,nb,nc,twoJ,twoJ12,twoT,E12,Lcm};
}

//void Jacobi3BME::TcoeffUnHash(std::array<unsigned short,7>& key, int& na, int& nb, int& nc, int& Jab, int& twoJ,  int& twoJ12, int& E12 )
//void Jacobi3BME::TcoeffUnHash(std::array<unsigned short,7>& key, int& na, int& nb, int& nc, int& twoJ,  int& twoJ12, int& E12, int& Lcm )
void Jacobi3BME::TcoeffUnHash(std::array<unsigned short,8>& key, int& na, int& nb, int& nc, int& twoJ,  int& twoJ12, int& twoT, int& E12, int& Lcm )
{
  na=key[0]; nb=key[1]; nc=key[2]; twoJ=key[3]; twoJ12=key[4]; twoT=key[5]; E12=key[6]; Lcm=key[7];
}






//std::string Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ,  uint64_t twoJ12, uint64_t E12 )
//{
//  std::ostringstream oss;
////  oss << na << " " << nb << " " << nc << " " << Jab << " " << twoJ << " " << jac1 << " " << jac2 << " " << twoJ12 << " " << Lcm;
//  oss << na << " " << nb << " " << nc << " " << Jab << " " << twoJ << " " <<  " " << twoJ12 << " " << E12;
//  return oss.str();
//}

//void Jacobi3BME::TcoeffUnHash(uint64_t key, uint64_t& na, uint64_t& nb, uint64_t& nc, uint64_t& Jab, uint64_t& twoJ, uint64_t& jac1, uint64_t& jac2, uint64_t& twoJ12, uint64_t& Lcm )
//void Jacobi3BME::TcoeffUnHash(std::string& key, int& na, int& nb, int& nc, int& Jab, int& twoJ, int& jac1, int& jac2, int& twoJ12, int& Lcm )
//void Jacobi3BME::TcoeffUnHash(std::string& key, int& na, int& nb, int& nc, int& Jab, int& twoJ,  int& twoJ12, int& E12 )
//{
//  std::istringstream iss(key);
////  iss >> na >> nb >> nc >> Jab >> twoJ >> jac1 >> jac2 >> twoJ12 >> Lcm;
//  iss >> na >> nb >> nc >> Jab >> twoJ >> twoJ12 >> E12;
//}



// Find all the entries in HartreeFock::Vmon3 that we need to compute for a specific set of one-body channels (la,j2a, ...)
// We want to get all the free ones we can from symmetry (aside from swapping the third index, e.g. abc -> cab. I haven't mustered the gumption for that yet).
// This is a bit awkward because the Vmon3 and Vmon3_keys are just vectors of doubles and uint64_t, respectively.
// The strategy is to loop through Vmon_keys, and pick out all the ones that belong to this one-body channel, or to the one with ab->ba, de->ed.
// We store these in the std::set "need_to_compute". We then loop through those and figure out which ones are related by symmetry.
// We put a set of 8 indices into the structure "indices", where we explicitly compute entry [0], and get the others for free.
// We then add all 8 to the std::set "will_be_computed" and continue looping. One the next iteration, if the entry from "need_to_compute" is already in "will_be_computed"
// we skip it. Since a std::set is a container of unique entries, we avoid computing the same thing twice.
//void Jacobi3BME::GetMonopoleIndices( int la, int j2a, int lb, int j2b, int lc, int j2c, HartreeFock& hf, std::vector<std::array<size_t,8>>& indices )
//void Jacobi3BME::GetMonopoleIndices( int la, int j2a, int lb, int j2b, int lc, int j2c, HartreeFock& hf, std::vector<std::unordered_set<size_t>>& indices )
void Jacobi3BME::GetMonopoleIndices( int la, int j2a, int lb, int j2b, int lc, int j2c, HartreeFock& hf, std::vector<std::vector<size_t>>& indices )
{
  std::unordered_set<size_t> need_to_compute;
  std::unordered_set<size_t> will_be_computed;
//  std::set<size_t> need_to_compute;
//  std::set<size_t> will_be_computed;

  bool verbose = false;

  if (verbose) std::cout << "---------------------------------------------------" << std::endl;
  if (verbose) std::cout << "lj channel : " << la << " " << j2a << " " << lb << " " << j2b << " " << lc << " " << j2c << std::endl;


  double t_internal = omp_get_wtime();

//  std::unordered_map<std::array<unsigned short,6>,size_t,array6_hash> imon_lookup;
  std::unordered_map<size_t,size_t> imon_lookup;

  // this parenthesis is just so that the need_to_compute_this_thread structure goes out of scope
  {
  int nthreads = omp_get_max_threads();
//  std::vector<std::set<size_t>> need_to_compute_this_thread(nthreads);
  std::vector<std::unordered_set<size_t>> need_to_compute_this_thread(nthreads);
  #pragma omp parallel for
  for (size_t imon=0; imon<hf.Vmon3_keys.size(); imon++)
  {
     int ithread = omp_get_thread_num();
     auto key = hf.Vmon3_keys[imon];
     int a,b,c,d,e,f;
     hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we could call it without a class instance if we wanted...
     Orbit& oa = hf.modelspace->GetOrbit(a);
     Orbit& ob = hf.modelspace->GetOrbit(b);
     Orbit& oc = hf.modelspace->GetOrbit(c);
     if ( 2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l > E3max) continue;

     if (  (oa.l==la and oa.j2==j2a and ob.l==lb and ob.j2==j2b and oc.l==lc and oc.j2==j2c)
        or (oa.l==lb and oa.j2==j2b and ob.l==la and ob.j2==j2a and oc.l==lc and oc.j2==j2c)
        or (oa.l==la and oa.j2==j2a and ob.l==lc and ob.j2==j2c and oc.l==lb and oc.j2==j2b)
        or (oa.l==lc and oa.j2==j2c and ob.l==la and ob.j2==j2a and oc.l==lb and oc.j2==j2b)
        or (oa.l==lb and oa.j2==j2b and ob.l==lc and ob.j2==j2c and oc.l==la and oc.j2==j2a)
        or (oa.l==lc and oa.j2==j2c and ob.l==lb and ob.j2==j2b and oc.l==la and oc.j2==j2a)
        )
       {
//          need_to_compute.insert(imon);
//          imon_lookup[key] = imon;
          need_to_compute_this_thread[ithread].insert(imon);
       }
  }

  // Now join them all together
  for ( auto& need : need_to_compute_this_thread )
  {
    for (size_t imon : need )
    {
      need_to_compute.insert(imon);
    }
  }
  }

  for (size_t imon : need_to_compute )
  {
      imon_lookup[hf.Vmon3_keys[imon]] = imon;
  }

  IMSRGProfiler::timer[std::string(__func__)+"_firstloop"] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();
//  std::cout << "size of Vmon_keys = " << hf.Vmon3_keys.size() << "   size of need_to_compute = " << need_to_compute.size() << std::endl;

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
     if ( (oa.tz2 + ob.tz2 + oc.tz2)>0 ) continue;

     std::set<std::array<int,6>> permuted_indices = { {a,b,c,d,e,f}, {b,a,c,e,d,f}, {c,a,b,f,d,e}, {a,c,b,d,f,e}, {c,b,a,f,e,d}, {b,c,a,e,f,d}}; 
//     std::set<std::array<int,6>> permuted_indices = { {a,b,c,d,e,f}, {b,a,c,e,d,f}, {c,a,b,f,d,e}, {a,c,b,d,f,e}, {c,b,a,f,e,d}, {b,c,a,e,f,d}, 
//                                                      {d,e,f,a,b,c}, {e,d,f,b,a,c}, {f,d,e,c,a,b}, {d,f,e,a,c,b}, {f,e,d,c,b,a}, {e,f,d,b,c,a} };
     std::set<std::array<int,6>> more_permutations;
//     more_permutations.clear();
     for (auto& p : permuted_indices)
     {
       std::array<int,6> isoflip; // Because we have isospin symmetry, we can flip the tz of all the nucleons and the matrix element is the same
       for (int i=0;i<6;i++)
       {
         Orbit& oi = hf.modelspace->GetOrbit(p[i]);
         isoflip[i] = hf.modelspace->GetOrbitIndex( oi.n, oi.l, oi.j2, -oi.tz2);
       }
       more_permutations.insert(isoflip);
     }
     for (auto& p : more_permutations) permuted_indices.insert(p); // add them to the list

     more_permutations.clear();
     for (auto& p : permuted_indices) more_permutations.insert({p[3],p[4],p[5],p[0],p[1],p[2]}); //Hermitian conjucation  <abc|V|def> = <def|V|abc>
     for (auto& p : more_permutations) permuted_indices.insert(p); // add them to the list

     
     std::set<uint64_t> permuted_keys;
     for (auto& p : permuted_indices )  permuted_keys.insert( hf.Vmon3Hash(p[0],p[1],p[2],p[3],p[4],p[5] ) ); // get the corresponding key for each permutation

     std::vector<size_t> compute_these;


     compute_these.push_back(imon); 

     for ( auto& k : permuted_keys )
     {
       if (k==key) continue;
       auto imon_iter = imon_lookup.find(k);
       if (imon_iter != imon_lookup.end())  compute_these.push_back(  imon_iter->second );
     }

//     for ( auto imon1 : need_to_compute )
//     {
//       if (imon1==imon) continue;
//       for (auto& k : permuted_keys )
//       {
//        if ( hf.Vmon3_keys[imon1] == k )  compute_these.push_back(imon1);
//       }
//     }




     for ( auto& index : compute_these ) will_be_computed.insert(index);
     indices.push_back( compute_these );

  } // for imon in need_to_compute

  IMSRGProfiler::timer[std::string(__func__)+"_secondloop"] += omp_get_wtime() - t_internal;

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
      for (auto index : ilist) std::cout << index << " ";
     std::cout << "  ) " << std::endl;
    }
  }

}


//void Jacobi3BME::GetRelevantTcoeffs( int la, int j2a, int lb, int j2b, int lc, int j2c, std::vector<std::array<unsigned short,7>>& lab_kets)   
void Jacobi3BME::GetRelevantTcoeffs( int la, int j2a, int lb, int j2b, int lc, int j2c, std::vector<std::vector<std::array<unsigned short,7>>>& lab_kets)   
{
//  std::cout << "I'm here in " << __func__ << std::endl;
  double t_internal = omp_get_wtime();
  bool verbose = true;
  TcoeffLookup.clear();

//  size_t dim_lab = lab_kets.size();
  for (int nth_pass=0; nth_pass<=1; nth_pass++) // take two passes, one to allocate, one to calculate, with the calculate pass in parallel
  {
   #pragma omp parallel for schedule(dynamic,1) collapse(4) if (nth_pass==1)
   for (unsigned short E12=0; E12<=Nmax; E12++)
   {
    for (unsigned short twoT=1; twoT<=3; twoT+=2)
    {
//     for (unsigned short twoJ12=1; twoJ12<=twoJmax; twoJ12+=2)
     for (unsigned short twoJ12=1; twoJ12<=twoJmax_cut; twoJ12+=2)
     {
       for (unsigned short Ecm=0; Ecm<=E3max; Ecm++)
       {
        int Eabc = E12+Ecm;
        if ( (Ecm+E12+Eabc)%2>0 ) continue;
        if ( Eabc > E3max ) continue;
        size_t dim_lab = lab_kets[Eabc].size();
        for (unsigned short Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
        {
         int Ncm = (Ecm-Lcm)/2;
         int parity = E12%2;
         size_t dimNAS = GetDimensionNAS(twoT,twoJ12, parity, E12 );
         size_t dimAS = GetDimensionAS(twoT,twoJ12, parity, E12 );
         std::array<unsigned short,5> T_key = {E12,twoT,twoJ12,Ecm,Lcm};
         if (nth_pass==0) // first pass through, not in parallel so we can allocate the structure
         {
//            TcoeffLookup[ {E12,twoT,twoJ12,Ecm,Lcm} ] = arma::mat( dim_lab, dimNAS, arma::fill::zeros );
//            TcoeffLookup[ T_key ] = arma::mat( dim_lab, dimNAS, arma::fill::zeros );
            TcoeffLookup[ T_key ] = arma::mat( dim_lab, dimAS, arma::fill::zeros ); 
         }
         else // this is inside the parallel block. No more allocation, just calculation and assignment
         {
//           auto& TcoeffMat = TcoeffLookup[ {E12,twoT,twoJ12,Ecm,Lcm} ];
//           auto& TcoeffMat = TcoeffLookup[ T_key ];
           arma::mat Tcoeff_tmp( dim_lab, dimNAS,arma::fill::zeros);
           if (dimAS<1 or dimNAS<1 or dim_lab<1) continue;

           for (size_t iNAS=0; iNAS<dimNAS; iNAS++)
           {
             jacobi1_state jac1;
             jacobi2_state jac2;
             GetJacobiStates( twoT, twoJ12, parity, E12, iNAS, jac1, jac2);
             
             for (size_t ilab=0; ilab<dim_lab; ilab++)
             {
               auto lab_arr = lab_kets[Eabc][ilab];
               int na=lab_arr[0], nb=lab_arr[1], nc=lab_arr[2], Jab=lab_arr[3], Tab=lab_arr[4], twoJlab=lab_arr[5], twoTlab=lab_arr[6];
//               int Eabc = 2*(na+nb+nc)+la+lb+lc;
//               if ( (2*(na+nb+nc)+la+lb+lc) != Eabc ) continue; // it will probably pay to use the block-diagonal energy structure. for now, lots of zeros.
               if ( (twoTlab!=twoT) or  (Tab != jac1.t) ) continue;
               if ( std::abs(twoJlab-twoJ12) > 2*Lcm  or twoJlab+twoJ12<2*Lcm ) continue;
               double tcoef = ComputeTcoeff( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJlab, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
//               TcoeffMat( ilab, iNAS ) = tcoef;
               Tcoeff_tmp( ilab, iNAS ) = tcoef;
             } // for ilab
           } // for iNAS
           size_t cfp_begin = GetCFPStartLocation(twoT,twoJ12,E12);
           arma::mat cfp_mat( &(cfpvec[cfp_begin]), dimNAS, dimAS, /*copy_aux_mem*/ false);
           TcoeffLookup[T_key] = Tcoeff_tmp * cfp_mat * sqrt(6.0);
        } // else (nth_pass check)
       } // for Lcm
      } // for Ecm
     } //for twoJ12
    } // for twoT
   } // for E12
  } // for nth_pass
  
}




/*

void Jacobi3BME::GetRelevantTcoeffs( int la, int j2a, int lb, int j2b, int lc, int j2c, HartreeFock& hf)   
{
    double t_internal = omp_get_wtime();
    bool verbose = false;
    if (verbose) std::cout << "In " <<__func__ << std::endl;
    std::set<int> na_list, nb_list,nc_list; // a std::set is a sorted unique list of items
    for (auto& orb : hf.modelspace->Orbits ) // it doesn't matter if we run through +-tz, because it's a set, so only stores one copy
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

//     for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
//     {
     if (verbose) std::cout << "Start loop ober twoJ12" << std::endl;
         int twoJ12_min = 1;
         int twoJ12_max = twoJmax;
         for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
         {
         for (int twoT=1; twoT<=3; twoT+=2)
         {
         for (int E12=0; E12<=Nmax; E12++)
         {
         int parity = E12%2;
         for ( auto na : na_list)
         {
           auto a = hf.modelspace->GetOrbitIndex( na, la, j2a, 1);
          for ( auto nb : nb_list )
          {
           auto b = hf.modelspace->GetOrbitIndex( nb, lb, j2b, 1);
//           if (a>b) continue;
           for (auto nc : nc_list )
           {
            int Eabc = 2*(na+nb+nc) + la+lb+lc;
            if (Eabc>E3max) continue;
            int Ecm = Eabc - E12;
            if (Ecm<0) continue;

            int twoJ_min = 1;
            int twoJ_max = 2*(la+lb+lc)+3;
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
            {
//            for (int Ecm=0; Ecm<=Eabc; Ecm++)
//            {
//             int E12 = Eabc-Ecm;
//             if (E12 > Nmax) continue;

                size_t dimNAS = GetDimensionNAS(twoT,twoJ12, parity, E12 );
                if (dimNAS<1) continue;
                for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
                {
                 if ( std::abs(twoJ-twoJ12)>2*Lcm  or (twoJ+twoJ12)<2*Lcm ) continue;
//                 TcoeffLookup[ {na,nb,nc,twoJ,twoJ12,E12,Lcm} ] = tcoeff_count;
                 TcoeffLookup[ {na,nb,nc,twoJ,twoJ12,twoT,E12,Lcm} ] = tcoeff_count;
                 tcoeff_count += (Jab_max+1-Jab_min) * dimNAS;
//                 std::cout << "TcoeffLookup = " << TcoeffLookup[{na,nb,nc,twoJ,twoJ12,E12,Lcm}] << "   tcoeff_count = " << tcoeff_count
                 if (verbose) std::cout << "TcoeffLookup = " << TcoeffLookup[{na,nb,nc,twoJ,twoJ12,twoT,E12,Lcm}] << "   tcoeff_count = " << tcoeff_count
                           << "Jab_max-Jab_min = " << Jab_max << " - " << Jab_min
                           << "  GetDimensionNAS(" << twoT << " " << twoJ12 << " " << parity << " " << E12  << ") =  dimNAS = " << dimNAS << std::endl;
                }// for Lcm
//            }// for Ecm
           } // for nc
          } // for nb
         } // for na
         } // for E12
         } // for twoT
         } // for twoJ12
       } // for twoJ
//     } // for Jab
       if (verbose) std::cout << "Done looping over twoJ12" << std::endl;

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
    if (verbose) std::cout << "Start loop over TcoeffLookup" << std::endl;
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
         int na,nb,nc,Jab,twoJ,twoJ12,twoT,E12,Lcm;
//         TcoeffUnHash(hash_key, na,nb,nc,Jab,twoJ,twoJ12,E12);
         TcoeffUnHash(hash_key, na,nb,nc,twoJ,twoJ12,twoT,E12,Lcm);
//         TcoeffUnHash(hash_key, na,nb,nc,twoJ,twoJ12,E12,Lcm);
         bool ab_same = (la==lb and j2a==j2b and na==nb);

         size_t start_point = element->second;
         if (verbose) std::cout << "na,nb,nc,Jab,twoJ,twoJ12,twoT,E12,Lcm = " << na << " " << nb << " " << nc << " " << Jab << " " << twoJ << " " << twoJ12 << " " << twoT << " " << E12 << " " << Lcm << "   start_point = " << start_point << std::endl;

         int Eabc = 2*(na+nb+nc)+la+lb+lc;
         int Ecm = Eabc-E12;
         int parity = E12%2;

         if ( std::abs(twoJ-twoJ12)>2*Lcm or (twoJ+twoJ12)<2*Lcm) continue;
         int Ncm = (Ecm-Lcm)/2;


//         for (int twoT=1; twoT<=3; twoT+=2)
//         {
           size_t dimNAS = GetDimensionNAS(twoT,twoJ12, parity, E12 );
           if (verbose) std::cout << "GetDimensionNAS( " << twoT << " " << twoJ12 << " " << parity << " " << E12 << ")  DIMNAS = " << dimNAS << std::endl;
           for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
           {
            for (size_t iNAS=0; iNAS<dimNAS; iNAS++)
            {
              jacobi1_state jac1;
              jacobi2_state jac2;
              GetJacobiStates( twoT, twoJ12, parity, E12, iNAS, jac1, jac2);
              if (ab_same and  (jac1.t+Jab)%2==0) continue;
              double tcoef = ComputeTcoeff( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
//              size_t location = start_point +(twoT-1)/2 * (Jab-Jab_min) * dimNAS + iNAS;
//              size_t location = start_point + ( (twoT-1)/2 * (Jab_max-Jab_min) + (Jab-Jab_min) ) * dimNAS + iNAS;
              size_t location = start_point + (Jab-Jab_min) * dimNAS + iNAS;
              if (verbose) std::cout << "iNAS = " << iNAS << "  start_point = " << start_point << "   Jab-Jab_min = " << Jab << " - " << Jab_min << "   dimNAS = " << dimNAS << "  location = " << location << "    size = " << TcoeffList.size() << std::endl;
//              TcoeffList[ location ] = tcoef; 
              TcoeffList.at( location ) = tcoef; 
              tcomputed++;
            }
           }
//         }
       }
     } 
     if (verbose )std::cout << "Done with loop over TcoeffList" << std::endl;
     IMSRGProfiler::timer[std::string(__func__)+"_ComputeToeff"] += omp_get_wtime() - t_internal;
     IMSRGProfiler::counter["Tcoefficients"] += tcoeff_count;
     IMSRGProfiler::counter["Tcoeff computed"] += tcomputed;

}
*/



/*

void Jacobi3BME::GetRelevantTcoeffs( int la, int j2a, int lb, int j2b, int lc, int j2c, HartreeFock& hf)   
{
    double t_internal = omp_get_wtime();
    bool verbose = false;
    std::set<int> na_list, nb_list,nc_list; // a std::set is a sorted unique list of items
    for (auto& orb : hf.modelspace->Orbits ) // it doesn't matter if we run through +-tz, because it's a set, so only stores one copy
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
            if (Eabc>E3max) continue;

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
//         if (verbose) std::cout << "hash: " << element->first << "  ->  " << element->second << std::endl;
 

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




*/



void Jacobi3BME::GetV3mon_all( HartreeFock& hf )
{
  double t_start = omp_get_wtime();


//        bool verbose = true;
        bool verbose = false;

  // The keys are already computed elsewhere and are passed in as input
  PreComputeMoshinsky1();
  PreComputeMoshinsky2();
  PreComputeSixJ();
  PreComputeNineJ();
//  PreComputeTwelveJ();
  hf.Vmon3.resize( hf.Vmon3_keys.size(), 0.);
  struct ljChannel{
     int l; int j2;
     bool operator == (const ljChannel& rhs){return (rhs.l==l and rhs.j2==j2); };
     bool operator <  (const ljChannel& rhs) const {return (l<rhs.l or (l==rhs.l and j2>rhs.j2) ) ;};
  };

  std::set<ljChannel> ljchannels;


  for ( auto& obc : hf.modelspace->OneBodyChannels )
  {
    ljChannel ljchan = {obc.first[0], obc.first[1] };
    
    ljchannels.insert( ljchan );
  }

//  ljchannels.insert( {5,11} ); // TODO: DELETE THIS!!!

  int n_mon = 0;
  size_t num_lj = ljchannels.size();
  for ( auto& ljchan_a : ljchannels )
  {
    int la = ljchan_a.l;
    int j2a = ljchan_a.j2;
    for ( auto& ljchan_b : ljchannels )
    {
      if (ljchan_b < ljchan_a ) continue;
      int lb = ljchan_b.l;
      int j2b = ljchan_b.j2;
      int Jab_min = std::abs(j2a-j2b)/2;
      int Jab_max = (j2a+j2b)/2;
//      std::cout << "la,j2a  lb,j2b " << la << " " << j2a << " " << lb << " " << j2b << std::endl;
      for ( auto& ljchan_c : ljchannels )
      {
        if (ljchan_c < ljchan_b ) continue;
        int lc = ljchan_c.l;
        int j2c = ljchan_c.j2;
        if (verbose) std::cout << "lc,j2c = " << lc << " " << j2c << std::endl;

//        if (la==0 and lb==1 and lc==1 and j2b==3 and j2c==3) verbose = true;
//        else verbose = false;

        double t_internal = omp_get_wtime();


        std::vector<std::vector<size_t>> imon_indices;
        GetMonopoleIndices(la, j2a, lb, j2b, lc, j2c, hf, imon_indices );
        if (verbose) std::cout << "done calling GetMonopoleIndices" << std::endl;

        IMSRGProfiler::timer[std::string(__func__)+"_FindMonKeys"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();

        if (verbose) std::cout << "done.  Now the loop over imon... Size of imonlist is " << imon_indices.size() << std::endl;

        if (imon_indices.size()<1) continue;

//        GetRelevantTcoeffs(la, j2a, lb, j2b, lc, j2c, hf); 
//        if (verbose) std::cout << "done gettin tcoeffs" << std::endl; 


        t_internal = omp_get_wtime();


        std::set<int> na_list, nb_list,nc_list; // a std::set is a sorted unique list of items
        for (auto& orb : hf.modelspace->Orbits )
        {
          if (orb.l==la and orb.j2==j2a) na_list.insert(orb.n); 
          if (orb.l==lb and orb.j2==j2b) nb_list.insert(orb.n); 
          if (orb.l==lc and orb.j2==j2c) nc_list.insert(orb.n); 
        }
        int Jab_min = std::abs(j2a-j2b)/2;
        int Jab_max = (j2a+j2b)/2;
        // I want a list of the |abc,Jab,Tab,JT> kets, and a way to map abc,Jab,Tab,JT to an index
//        std::vector<std::array<unsigned short, 7>> lab_kets;
        std::vector<std::vector<std::array<unsigned short, 7>>> lab_kets;
        std::unordered_map<std::array<unsigned short,7>,size_t,array7_hash> lab_ket_lookup;
        lab_kets.resize(E3max+1);
        if (verbose) std::cout << "Start loop to fill lab_kets" << std::endl;

        for ( auto na : na_list )
        {
          for ( auto nb : nb_list )
          {
            for ( auto nc : nc_list )
            {
              int Eabc = 2*(na+nb+nc)+la+lb+lc;
              if ( Eabc > E3max ) continue;
//              if (la==lb and la==lc and j2a==j2b and j2a==j2c and na==nb and na==nc and j2a<3 ) continue; // cant fit 
              for ( int Jab=Jab_min; Jab<=Jab_max; Jab++ )
              {
              for ( int Tab=0; Tab<=1; Tab++ )
              {
               if (la==lb and j2a==j2b and na==nb and (Tab+Jab)%2==0 ) continue;
               for ( int twoT=1; twoT<=2*Tab+1; twoT+=2 )
               {
                 for ( int twoJ=std::abs(2*Jab-j2c); twoJ<=(2*Jab+j2c); twoJ+=2)
                 {
                    if (la==lb and la==lc and j2a==j2b and j2a==j2c and na==nb and na==nc and ( twoJ==3*j2a or twoT>j2a) ) continue;
                    auto key = MakeUshort7({na,nb,nc,Jab,Tab,twoJ,twoT});
                    lab_ket_lookup[key] = lab_kets[Eabc].size();
                    lab_kets[Eabc].push_back( key );
                  }
                }
               }
              }
            }
          }
        }
        if (verbose) std::cout << "done with that." << std::endl;
        arma::field<arma::mat> lab_mats(E3max+1,E3max+1);
        for (size_t Eabc=0; Eabc<=E3max; Eabc++)
        {
          for (size_t Edef=0; Edef<=E3max; Edef++)
          {
            lab_mats(Eabc,Edef).zeros( lab_kets[Eabc].size(), lab_kets[Edef].size() );
          }
        }

        IMSRGProfiler::timer[std::string(__func__)+"lab_kets_loop"] += omp_get_wtime() - t_internal;

        t_internal = omp_get_wtime();
        GetRelevantTcoeffs(la, j2a, lb, j2b, lc, j2c, lab_kets); 
        if (verbose) std::cout << "done gettin tcoeffs" << std::endl; 
        IMSRGProfiler::timer[std::string(__func__)+"getTcoeffs"] += omp_get_wtime() - t_internal;

        t_internal = omp_get_wtime();
        if (verbose) std::cout << "created and filled lab_ket_lookup and lab_kets  size =" << lab_kets.size() << std::endl;

        #pragma omp parallel for schedule(dynamic,1) collapse(3)
        for (int Ecm=0; Ecm<=E3max; Ecm++ )
        {
         for (int twoT=1; twoT<=3; twoT+=2)
         {
//          for (int twoJ12=1; twoJ12<=twoJmax; twoJ12+=2)
          for (int twoJ12=1; twoJ12<=twoJmax_cut; twoJ12+=2)
          {
           for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
           {
            if ( (j2a+j2b+j2c+2*Lcm) < twoJ12 ) continue;
            for (int E12abc=0; E12abc<=std::min(Nmax,E3max-Ecm); E12abc++)
            {
              int Eabc = E12abc + Ecm;
              if ( (E12abc + Ecm + Eabc)%2>0) continue;
              if (verbose) std::cout << "Ecm,Lcm,twoT,twoJ12,E12abc = " << Ecm << " " << Lcm << " " << twoT << " " << twoJ12 << " " << E12abc << std::endl;
              int parity=E12abc%2;
              size_t dimNAS_abc = GetDimensionNAS( twoT, twoJ12, parity, E12abc ); 
              size_t dimAS_abc = GetDimensionAS( twoT, twoJ12, parity, E12abc ); 
              if (dimNAS_abc==0 or dimAS_abc==0) continue;

              auto Tkey_abc = MakeUshort5({E12abc,twoT,twoJ12,Ecm,Lcm});
              arma::mat& Tabc = TcoeffLookup[ Tkey_abc ] ;
              if (arma::norm(Tabc,"fro")<1e-8) continue;


              for (int E12def=E12abc; E12def<=std::min(Nmax,E3max-Ecm); E12def+=2)
              {
                int Edef = E12def + Ecm;
                if ( (E12def + Ecm + Edef)%2>0) continue;
                size_t dimNAS_def = GetDimensionNAS( twoT, twoJ12, parity, E12def ); 
                size_t dimAS_def = GetDimensionAS( twoT, twoJ12, parity, E12def ); 
                if (dimNAS_def==0 or dimAS_def==0) continue;
                auto Tkey_def = MakeUshort5({E12def,twoT,twoJ12,Ecm,Lcm}); 
                arma::mat& Tdef = TcoeffLookup[ Tkey_def ] ;

                if (verbose) std::cout << " E12def,E12abc, twoJ12,twoT,Lcm,Ecm = " << E12def << " " << E12abc << " " << twoJ12 << " " << twoT << " " << Lcm << " " << Ecm << std::endl;

                size_t startlocAS = GetStartLocAS( twoT, twoJ12, E12abc, E12def);
                arma::mat matelAS( &meAS[startlocAS], dimAS_abc, dimAS_def, false ); 
//                arma::mat local_copy = Tabc * matelAS * cfp_def.t() * Tdef.t(); 
                arma::mat local_copy = Tabc * matelAS * Tdef.t(); 
                #pragma omp critical
                {
                  lab_mats(Eabc,Edef) += local_copy; // the add bit should be fast, so hopefully the blocking isn't too bad
                }

                if (verbose )
                {
                std::cout << "Tabc: " << std::endl << Tabc << std::endl << std::endl
                          << "matelAS: " << std::endl << matelAS << std::endl << std::endl
//                          << "cfp_def.t(): " << std::endl << cfp_def.t() << std::endl << std::endl
                          << "Tdef: " << std::endl << Tdef << std::endl << std::endl
                          << "lab_mat: " << std::endl << lab_mats << std::endl << std::endl << std::endl;
                }


              } // for E12def
            } // for E12abc
           } // for Lcm
          } // for twoJ12
         } // for twoT
        } // for Ecm
//        std::cout << "out of Ecm loop" << std::endl;

        // now fill out the lower-diagonal blocks of the matrix (or "field" of matrices)
        for (int Eabc=0;Eabc<=E3max;Eabc++)
        {
         for (int Edef=0; Edef<Eabc; Edef++)
         {
           lab_mats(Eabc,Edef) = lab_mats(Edef,Eabc).t() ;
         }
        }

        IMSRGProfiler::timer[std::string(__func__)+"matmult_loop"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();

        // Now we need to construct the monopole terms

        if (verbose) std::cout << "Start loop over imon" << std::endl;
        for (size_t ilist=0; ilist<imon_indices.size(); ilist++)
        {
          size_t imon = imon_indices[ilist][0];
          auto key = hf.Vmon3_keys[imon];
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          Orbit& oa = hf.modelspace->GetOrbit(a);
          Orbit& ob = hf.modelspace->GetOrbit(b);
          Orbit& oc = hf.modelspace->GetOrbit(c);
          Orbit& od = hf.modelspace->GetOrbit(d);
          Orbit& oe = hf.modelspace->GetOrbit(e);
          Orbit& of = hf.modelspace->GetOrbit(f);
          bool same_ab = ( (oa.n==ob.n) and (oa.l==ob.l) and (oa.j2==ob.j2) );
          bool same_ac = ( (oa.n==oc.n) and (oa.l==oc.l) and (oa.j2==oc.j2) );
          bool same_de = ( (od.n==oe.n) and (od.l==oe.l) and (od.j2==oe.j2) );
          bool same_df = ( (od.n==of.n) and (od.l==of.l) and (od.j2==of.j2) );
          int Eabc = 2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l;
          int Edef = 2*(od.n+oe.n+of.n)+od.l+oe.l+of.l;
          if (Eabc>E3max or Edef>E3max) continue;

          int Tzab = (oa.tz2 + ob.tz2)/2;
          int Tzde = (od.tz2 + oe.tz2)/2;
          int twoTz = 2*Tzab + oc.tz2;
          if (verbose) std::cout << std::endl << "abcdef " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;

          double v3mon = 0;
          int Jab_min = std::abs(oa.j2-ob.j2)/2;
          int Jab_max = (oa.j2+ob.j2)/2;
          for (int twoT=1; twoT<=3; twoT++)
          {
           for (int Tab=twoT/2; Tab<=1; Tab++)
           {
            double iso_clebsch_abc = AngMom::CG(0.5,0.5*oa.tz2,0.5,0.5*ob.tz2,Tab,Tzab) * AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2,0.5*twoT,0.5*twoTz);
            if (std::abs(iso_clebsch_abc)<1e-8) continue;
            for (int Tde=twoT/2; Tde<=1; Tde++)
            {
              double iso_clebsch_def = AngMom::CG(0.5,0.5*od.tz2,0.5,0.5*oe.tz2,Tde,Tzde) * AngMom::CG(Tde,Tzde,0.5,0.5*of.tz2,0.5*twoT,0.5*twoTz);
              if (std::abs(iso_clebsch_def)<1e-8) continue;
              for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
              {
                if (same_ab and (Tab+Jab)%2==0) continue;
                if (same_de and (Tde+Jab)%2==0) continue;
//                if (a==b and oa.j2==Jab) continue;
//                if (d==e and od.j2==Jab) continue;
                int twoJ_min = std::abs(2*Jab-oc.j2);
                int twoJ_max = 2*Jab+oc.j2;
                for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
                {
                  if ( ( (same_ab and same_ac) or (same_de and same_df) ) and  ( twoJ==3*oa.j2 or twoT>oa.j2)  ) continue;
                 if (verbose) std::cout << "twoT,Tab,Tde,Jab,twoJ " << twoT << " " << Tab << " " << Tde << " "<< Jab << " " << twoJ << std::endl;
                  size_t iket_abc = lab_ket_lookup[MakeUshort7({oa.n,ob.n,oc.n,Jab,Tab,twoJ,twoT})];
                  size_t iket_def = lab_ket_lookup[MakeUshort7({od.n,oe.n,of.n,Jab,Tde,twoJ,twoT})];
//                   if (verbose) std::cout << "iket_abc, iket_def = " << iket_abc << " " << iket_def << " dimensions of lab_mat " << lab_mat.n_rows << "x" << lab_mat.n_cols << std::endl;
//                  v3mon += (twoJ+1) * iso_clebsch_abc  * lab_mat(iket_abc,iket_def) * iso_clebsch_def ;
                  v3mon += (twoJ+1) * iso_clebsch_abc  * lab_mats(Eabc,Edef)(iket_abc,iket_def) * iso_clebsch_def ;
                  if (verbose) std::cout << "v3mon += (" << twoJ << "+1) * " << iso_clebsch_abc << " * " << lab_mats(Eabc,Edef)(iket_abc,iket_def) << " * " << iso_clebsch_def << " = " << v3mon << std::endl << std::endl;;
                }
              }
             }
           }
          }
          for ( auto& imon_sym : imon_indices[ilist] )
          {
              int aa,bb,cc,dd,ee,ff;
              auto key = hf.Vmon3_keys[imon_sym];
              hf.Vmon3UnHash(key, aa,bb,cc,dd,ee,ff);
              hf.Vmon3[imon_sym] = v3mon / (hf.modelspace->GetOrbit(cc).j2+1);  // don't forget that 2jc+1 factor...
              if (verbose) std::cout << "assigning imon_sym = " << imon_sym << "  = " << v3mon << " / " << (hf.modelspace->GetOrbit(cc).j2+1) << " = " << hf.Vmon3[imon_sym] << "  abcdef " << aa << " " << bb << " " << cc << " " << dd << " " << ee << " " << ff << std::endl;
          }
          if (verbose) std::cout << "done assigning" << std::endl;
        } // for ilist
        if (verbose) std::cout << "done with loop over ilist" << std::endl;

        IMSRGProfiler::timer[std::string(__func__)+"vmon_loop"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();


        int tcoeff_counter = 0;
        int nonzero_vmon = 0;

        n_mon += imon_indices.size();
        if (verbose) std::cout << "Looked up a Tcoefficient " << tcoeff_counter << "  times " << std::endl;
        if (verbose) std::cout << "Found " << nonzero_vmon << "  nonzero monopoles" << "  out of " << imon_indices.size() << " terms" << std::endl;
        if (verbose) std::cout << " and it took " << omp_get_wtime() - t_internal << "  seconds" << std::endl;

        IMSRGProfiler::timer[std::string(__func__)+"_iMonLoop"] += omp_get_wtime() - t_internal;

        if (verbose) std::cout << "finished the loop for this channel" << std::endl;
      }// for ilj_c
    }// for ilj_b
  }// for ilj_a
//  if (verbose) std::cout << "FINISHED THE LOOP FOR ALL CHANNELS" << std::endl;

//  std::cout << "all done" << std::endl;

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  IMSRGProfiler::counter["V3Monopoles"] += n_mon;
}





void Jacobi3BME::GetNO2b_all( ModelSpace& ms )
{

  struct ket3b { int a; int b; int c; int Tab;};
  bool verbose = false;
  PreComputeMoshinsky1();
  PreComputeMoshinsky2();
  PreComputeSixJ();
  PreComputeNineJ();

  // Make a list of orbits that could be in rho
  std::set<int> occupied_orbits; // really just potentially occupied orbits
  std::vector<std::pair<int,int>> rho_pairs;
  for ( auto& a : ms.all_orbits )
  {
    Orbit& oa = ms.GetOrbit(a);
    if (oa.l>lmax_NO2b) continue;
    occupied_orbits.insert(a);
    for ( auto b : ms.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
    {
      rho_pairs.push_back( std::make_pair(a,b) );
    }
  }

  // Allocate the vector of relevant 3bmes. This bit also demonstrates how to loop through and interpret the vector
  size_t n_elements = 0;
  int nch=ms.GetNumberTwoBodyChannels();
  std::vector<size_t> ch_begin(nch,0);
  for (int ch=0; ch<nch; ch++)
  {
    ch_begin.at(ch) = n_elements; // Location of first entry in matel for this channel
    TwoBodyChannel& tbc = ms.GetTwoBodyChannel(ch);
    int nkets = tbc.GetNumberKets();
    for (int ibra=0; ibra<nkets; ibra++)
    {
      Ket& bra = tbc.GetKet(ibra);
      int Ebra = 2*(bra.op->n+bra.oq->n) + bra.op->l + bra.oq->l;
      if (Ebra > E3max) continue;
      for (int iket=ibra; iket<nkets; iket++)
      {
        Ket& ket = tbc.GetKet(iket);
        int Eket = 2*(ket.op->n+ket.oq->n) + ket.op->l + ket.oq->l;
        if (Eket>E3max) continue;
        for (auto& ab : rho_pairs )
        {
          Orbit& oa = ms.GetOrbit(ab.first);
          Orbit& ob = ms.GetOrbit(ab.second);
          if ( 2*oa.n+oa.l + Ebra > E3max) continue;
          if ( 2*ob.n+ob.l + Eket > E3max) continue;
          n_elements++;
        }
      }
    }
  }
  std::cout << "Resizing matelNO2b to have " << n_elements << " elements" << std::endl;
  matelNO2b.resize(n_elements,0.0);

  // Loop over 2-body channels
  // Loop over possible density matrix entries <a|rho|b>
  for (int ch=0; ch<nch; ch++)
  {
    TwoBodyChannel& tbc = ms.GetTwoBodyChannel(ch);
    int Jab = tbc.J;
    int Tzab = tbc.Tz;
    int Tab_min = std::abs( Tzab );
    int Tab_max = 1;
    int nkets = tbc.GetNumberKets(); 
    // Make a list of relevant 3-body kets for this channel

    // hash table (twoJ twoT Eabc) -> {vector of lab-frame 3b kets in that channel}    --- local is probably a bad choice of terminology...
    std::unordered_map<std::array<unsigned short,3>,std::vector<ket3b>,array3_hash> local_3b_kets;
    // hash table (a,b,c,Tab,twoJ,twoT) -> index in the matrix where it lives (need to infer Eabc from a,b,c)
    std::unordered_map<std::array<unsigned short,6>,size_t,array6_hash > ket_3b_lookup;  // hash table for reverse lookups
    // structure for the lab-frame 3b matrix elements. (twoJ twoT parity) -> field (Eabc,Edef) -> matrix of 3-body matrix elements
    std::unordered_map<std::array<unsigned short,3>, arma::field<arma::mat>, array3_hash> V3Full;

    if (verbose) std::cout << "Begin loop enumerating 3-body states" << std::endl;
    double t_internal = omp_get_wtime();
    // Go through all the kets in this 2-body channel and figure out which 3-body states we'll need.
    for (int iket=0; iket<nkets; iket++)
    {
      Ket& ket = tbc.GetKet(iket);
      int a = ket.p;
      int b = ket.q;
      if (verbose) std::cout << "iket = " << iket << "  p,q = " << a << " " << b << std::endl;
      for ( int c : occupied_orbits ) // occupied orbits -> oscillator orbits which end up with a non-zero occupation in the HF reference
      {
        Orbit& oc = ms.GetOrbit(c);
        int Eabc = 2*(ket.op->n+ket.oq->n+oc.n)+ket.op->l+ket.oq->l+oc.l ;
        if ( Eabc > E3max ) continue;
        int parity = (ket.op->l + ket.oq->l + oc.l)%2;
        int twoTz = 2*Tzab + oc.tz2;
        int twoT_min = std::abs(twoTz);
        int twoT_max = 3;
        int twoJ_min = std::abs( oc.j2 - 2*Jab );
        int twoJ_max = oc.j2 + 2*Jab;
        for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
        {
         for (int twoT=twoT_min; twoT<=twoT_max; twoT+=2)
         {
          for (int Tab=Tab_min; Tab<=Tab_max; Tab++)
          {
           if ( (twoT-2*Tab)>1 ) continue;
           auto& channel_vec = local_3b_kets[ MakeUshort3({twoJ,twoT,Eabc})];  // find the J,T,E channel this 3b state belongs to
           ket_3b_lookup[ MakeUshort6({a,b,c,Tab,twoJ,twoT}) ] = channel_vec.size();  // this tells us where in the matrix this 3b state will live
           channel_vec.push_back( {a,b,c,Tab}); // add the 3b state to the vector of 3b states in this J,T,E channel
//           std::cout << "Added state |" << a << " " << b << " " << c << " " << Tab << ">  to channel  ( " << twoJ << " " << twoT << " " << Eabc << " ) " << std::endl;
           auto V3key = MakeUshort3({twoJ,twoT,parity});
           if ( V3Full.find( V3key ) == V3Full.end())
           {
             // Allocate the arma::field of matrices for this J,T,p channel.
             // We only use ether the even or the odd Eabc, but since we don't allocate the others I don't expect this to be much of an issue.
             V3Full[V3key] = arma::field<arma::mat>(E3max+1,E3max+1); 
           }
          }
         }
        }
      }
    }
    // Done enumerating 3-body states

    if (verbose) std::cout << "Allocate matrices in V3Full" << std::endl;
    // Allocate the <abc,Tab|def,Tde> matrices that live in the (Eabc,Edef) fields that live in the JTp hash table...Yikes.
    for ( auto& iter : V3Full )
    {
      int twoJ = iter.first[0];  int twoT = iter.first[1]; int parity = iter.first[2];
      for (int Eabc=parity; Eabc<=E3max; Eabc+=2)
      {
        size_t dim_abc = local_3b_kets[ MakeUshort3({twoJ,twoT,Eabc})].size();
        for (int Edef=parity; Edef<=E3max; Edef+=2)
        {
          size_t dim_def = local_3b_kets[ MakeUshort3({twoJ,twoT,Edef})].size();
          iter.second(Eabc,Edef).zeros(dim_abc, dim_def);
        }
      }
    }
    // done allocating the matrices in V3Full


    // Now we allocate (first pass) and compute (second pass) the relevant T coefficients
    std::unordered_map<std::array<unsigned short,6>, arma::mat, array6_hash> Tcoeffs;

    if (verbose) std::cout << "Begin loop allocating and computing Tcoefficients" << std::endl;
    for (int nth_pass=0; nth_pass<=1; nth_pass++)
    {
     if (verbose) std::cout << "  pass " << nth_pass << std::endl;
     #pragma omp parallel if(nth_pass==1)
     {
         size_t count = 0;
         int ithread = omp_get_thread_num();
         int nthreads = omp_get_num_threads();
//      for (auto& iterJTE : local_3b_kets)
      for (auto element = local_3b_kets.begin(); element !=local_3b_kets.end(); ++element, count++)
      {
        if(nth_pass==1 and count%nthreads != ithread) continue; // Check if this thread should compute this element
        int twoJ=element->first[0], twoT=element->first[1], Eabc=element->first[2];
        size_t dim_abc = element->second.size();
//        if (verbose) std::cout << "||labket channel JTE " << twoJ << " " << twoT << " " << Eabc << "  dimension " << dim_abc << std::endl;
        for (int E12=0; E12<=std::min(Nmax,Eabc); E12++)
        {
          int Ecm = Eabc-E12;
          for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
          {
//            if (verbose) std::cout << "..E12,Ecm,Lcm = " << E12 << " " << Ecm << " " << Lcm << std::endl;
            int Ncm = (Ecm-Lcm)/2;
            int twoJ12_min = std::abs(twoJ-2*Lcm);
            int twoJ12_max = std::min(twoJ+2*Lcm,twoJmax);
            for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
            {
//              if (verbose) std::cout << "checking dimensions for " << twoT << " " << twoJ12 << " " << E12%2 << " " << E12 << std::endl;
              size_t dimAS = GetDimensionAS( twoT, twoJ12, E12%2, E12 );
              size_t dimNAS = GetDimensionNAS( twoT, twoJ12, E12%2, E12 );
              if (dimNAS<1) continue;
              if (nth_pass==0) // on the first pass, we allocate
              {
//                Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12,Lcm}) ] = arma::mat( dim_abc, dimNAS, arma::fill::zeros);
                Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12,Lcm}) ] = arma::mat( dim_abc, dimAS, arma::fill::zeros);
//                if (verbose) std::cout << "allocated a tcoeff matrix for J,T,Eabc J12,E12,Lcm = " << twoJ << " " << twoT << " " << Eabc << "  " << twoJ12 << " " << E12 << " " << Lcm << "  with size " << dim_abc << " x " << dimNAS << std::endl;
                if (verbose) std::cout << "allocated a tcoeff matrix for J,T,Eabc J12,E12,Lcm = " << twoJ << " " << twoT << " " << Eabc << "  " << twoJ12 << " " << E12 << " " << Lcm << "  with size " << dim_abc << " x " << dimAS << std::endl;
              }
              if (nth_pass==1) // on the second pass, we compute all the matrix elements
              {
//                auto& TcoeffMat = Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12,Lcm}) ];
                arma::mat T_mat_tmp( dim_abc, dimNAS, arma::fill::zeros);
                if (verbose) std::cout << "start loop dimNAS = " << dimNAS << "  dim_abc = " << dim_abc << std::endl;
                for (size_t iNAS=0; iNAS<dimNAS; iNAS++)
                {
                  jacobi1_state jac1;
                  jacobi2_state jac2;
                  GetJacobiStates( twoT, twoJ12, E12%2, E12, iNAS, jac1, jac2);
                  if (verbose) std::cout << "iNAS = " << iNAS << std::endl;
                  for (size_t iabc=0; iabc<dim_abc; iabc++)
                  {
                    if (verbose) std::cout << "iabc = " << iabc << std::endl;
                    auto& ket_abc = element->second[iabc];
                    Orbit& oa = ms.GetOrbit(ket_abc.a);
                    Orbit& ob = ms.GetOrbit(ket_abc.b);
                    Orbit& oc = ms.GetOrbit(ket_abc.c);
                    if ( ket_abc.Tab != jac1.t ) continue;
                    if ( 2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l > E3max) continue;
                    if (verbose) std::cout << "calling ComputeTcoeff" << std::endl;
                    double tcoeff = ComputeTcoeff( oa.n, oa.l, oa.j2, ob.n, ob.l, ob.j2, oc.n, oc.l, oc.j2, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
//                    TcoeffMat(iabc,iNAS) = tcoeff;
                  if (verbose) std::cout << "assigning to " << iabc << " , " << iNAS << std::endl;
                    T_mat_tmp(iabc,iNAS) = tcoeff;
                  }
                }
                if (verbose) std::cout << "finished loop" << std::endl;
               //TODO CFPs here 
                size_t cfp_begin = GetCFPStartLocation(twoT,twoJ12,E12);
                if (verbose) std::cout << "assigning cfp_mat" << std::endl;
                arma::mat cfp_mat( &(cfpvec[cfp_begin]), dimNAS, dimAS, false); // false refers to copy_aux_mem
                if (verbose) std::cout << "assigning Tcoeffs[ ] " << std::endl;
                Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12,Lcm}) ] = sqrt(6) * T_mat_tmp * cfp_mat;
              }
            }
          }
        }
      }
     }
    IMSRGProfiler::timer[std::string(__func__)+"_Tcoeffs_"+std::to_string(nth_pass)] += omp_get_wtime() - t_internal;
    t_internal = omp_get_wtime();
    }
    if (verbose) std::cout << "Done computing Tcoefficients" << std::endl;

    int Eabc_min = std::max( Jab-1, 0); // for Jab > 1, we can't do it with  Eabc=0, so don't bother worrying about that
    
    if (verbose) std::cout << std::endl << "===============================================" << std::endl
                           << "start loops doing the mat mults. Eabc_min = " << Eabc_min << "based on Jab = " << Jab << std::endl;
   
   // Next we do the mat mult to perform the transformation
    for (int Ecm=0; Ecm<=E3max; Ecm++)
    {
     for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
     {
//       if (verbose) std::cout << " Ecm,Lcm = " << Ecm << " " << Lcm << std::endl;
       int Ncm = (Ecm-Lcm)/2;
//       if (verbose) std::cout << "J12 range " << std::max(1,2*Jab-2*emax-1-2*Lcm) << " to " << std::min(2*Jab+2*emax+1+2*Lcm,twoJmax) << "  from Jab=" << Jab << " emax = " << emax << "  twoJmax = " << twoJmax << "  and Lcm = " << Lcm << std::endl;
       for (int twoJ12=std::max(1,2*Jab-2*emax-1-2*Lcm); twoJ12<=std::min(2*Jab+2*emax+1+2*Lcm,twoJmax); twoJ12+=2)
       {
        for (int twoT=1; twoT<=3; twoT+=2)
        {
         for (int E12abc=0; E12abc<=std::min(Nmax,E3max-Ecm); E12abc++)
         {
//           if (verbose) std::cout << "~~  Ecm,Lcm,J12,T,E12abc  " << Ecm << " " << Lcm << " " << twoJ12 << " " << twoT << " " << E12abc << std::endl;
           size_t dimNAS_abc = GetDimensionNAS( twoT, twoJ12, E12abc%2, E12abc ); 
           size_t dimAS_abc = GetDimensionAS( twoT, twoJ12, E12abc%2, E12abc ); 
           if (dimNAS_abc==0 or dimAS_abc==0) continue;
           for (int E12def=E12abc; E12def<=std::min(Nmax,E3max-Ecm); E12def+=2 ) // we can probably make use of some hermiticity here
           {
             size_t dimNAS_def = GetDimensionNAS( twoT, twoJ12, E12def%2, E12def ); 
             size_t dimAS_def = GetDimensionAS( twoT, twoJ12, E12def%2, E12def ); 
             if (dimNAS_def==0 or dimAS_def==0) continue;
   
             size_t startlocAS = GetStartLocAS( twoT, twoJ12, E12abc, E12def);
             arma::mat matelAS( &meAS[startlocAS], dimAS_abc, dimAS_def, false ); 
   
//             if (verbose) std::cout << "...done" << std::endl;
   
             // Now we loop over the lab-frame matrix elements and update them with the contribution from this jacobi channel
             int twoJ_min = std::max( std::abs(twoJ12-2*Lcm), 2*Jab-2*emax-1);
             int twoJ_max = std::min( twoJ12+2*Lcm, 2*Jab+2*emax+1 );
//             if (verbose) std::cout << "running twoJ from " << twoJ_min << " to " << twoJ_max << "  based on Jab = " << Jab << " J12 = " << twoJ12 << " emax= " << emax << "  Lcm = " << Lcm << std::endl;
             int lab_parity = (E12abc+Ecm)%2;
             for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
             {
               auto& V3field = V3Full[MakeUshort3({twoJ,twoT,lab_parity})];
               int Eabc = E12abc + Ecm;
               int Edef = E12def + Ecm;
                 size_t dim_abc = local_3b_kets[ MakeUshort3({twoJ,twoT,Eabc})].size();
                 if (dim_abc<1) continue;
                 arma::mat matelAS_abc =  Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12abc,Lcm}) ] * matelAS;
   
                   size_t dim_def = local_3b_kets[ MakeUshort3({twoJ,twoT,Edef})].size();
                   if (dim_def<1) continue;
                   arma::mat local_lab_mat = matelAS_abc * Tcoeffs[ MakeUshort6({twoJ,twoT,Edef, twoJ12,E12def,Lcm}) ].t();
//                 if (verbose and twoJ==3 and Eabc==1 and Edef==1 and twoT==3)
//                  {
//                    std::cout << " local_lab_mat.  matelNAS_abc = " << std::endl << matelAS_abc << std::endl << "Tdef = " << std::endl
//                                         << Tcoeffs[ MakeUshort6({twoJ,twoT,Edef, twoJ12,E12def,Lcm}) ] << std::endl << "local_lab_mat = " << std::endl
//                                         << local_lab_mat << std::endl;
//                  }
//                   #pragma omp critical
//                   {
                     V3field(Eabc,Edef) += local_lab_mat;
                     if (Eabc != Edef)   V3field(Edef,Eabc) += local_lab_mat.t();
//                   }
               } // for twoJ
             } // for E12def
           } // for E12abc
         } // for twoT
       } // for twoJ12
     } // for Lcm
   } // for Ecm
   if (verbose) std::cout << "...done" << std::endl;
   
   IMSRGProfiler::timer[std::string(__func__)+"_matmult"] += omp_get_wtime() - t_internal;



  // Store it in the vector.
  // This appears to be very inefficient
  if (verbose) std::cout << "Now storing the vector" << std::endl;
  size_t istart = ch_begin.at(ch);
  size_t offset = 0;
  for (int ibra=0; ibra<nkets; ibra++)
  {
    if (verbose) std::cout << "ibra =  "<< ibra << std::endl;
    Ket& bra = tbc.GetKet(ibra);
    int Ebra = 2*(bra.op->n+bra.oq->n) + bra.op->l + bra.oq->l;
    int p = bra.p;
    int q = bra.q;
    if (Ebra > E3max) continue;
    for (int iket=ibra; iket<nkets; iket++)
    {
      if (verbose) std::cout << "  iket = " << iket << std::endl;
      Ket& ket = tbc.GetKet(iket);
      int r = ket.p;
      int s = ket.q;
      int Eket = 2*(ket.op->n+ket.oq->n) + ket.op->l + ket.oq->l;
      if (Eket>E3max) continue;
      for (auto& ab : rho_pairs )
      {
        int a = ab.first;
        int b = ab.second;
        if (verbose) std::cout << "a,b = " << a << " " << b << std::endl;
        Orbit& oa = ms.GetOrbit(a);
        Orbit& ob = ms.GetOrbit(b);
        int Epqa = Ebra + 2*oa.n+oa.l;
        int Ersb = Eket + 2*ob.n+ob.l;
        if ( Epqa > E3max) continue;
        if ( Ersb > E3max) continue;
        // Sum over twoJ, twoT, Tpq, Trs
        double me_no2b = 0;
        int twoJ_min = std::abs(2*Jab - oa.j2);
        int twoJ_max = 2*Jab + oa.j2;
        for (int Tpq=Tab_min; Tpq<=Tab_max; Tpq++)
        {
         double isoClebsch_pq = AngMom::CG(0.5,0.5*bra.op->tz2, 0.5,0.5*bra.oq->tz2, Tpq, Tzab);
         for (int Trs=Tab_min; Trs<=Tab_max; Trs++)
         {
          double isoClebsch_rs = AngMom::CG(0.5,0.5*ket.op->tz2, 0.5,0.5*ket.oq->tz2, Trs, Tzab);
          for (int twoT=std::abs(Tzab*2+oa.tz2); twoT<=std::min(Tpq,Trs)+1; twoT+=2)
          {
           double isoClebschTpq = AngMom::CG(Tpq,Tzab, 0.5,0.5*oa.tz2, 0.5*twoT, Tzab+0.5*oa.tz2);
           double isoClebschTrs = AngMom::CG(Trs,Tzab, 0.5,0.5*oa.tz2, 0.5*twoT, Tzab+0.5*oa.tz2);
           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
           {
             auto lab_pqa_key = MakeUshort6({p,q,a,Tpq,twoJ,twoT});
             size_t ipqa = ket_3b_lookup.at(lab_pqa_key);
             auto lab_rsb_key = MakeUshort6({r,s,b,Trs,twoJ,twoT});
             size_t irsb = ket_3b_lookup.at(lab_rsb_key);
             auto V3key = MakeUshort3({twoJ,twoT,Epqa%2});
             double me3b = V3Full.at(V3key)(Epqa,Ersb)(ipqa,irsb);
             me_no2b += (twoJ+1) * isoClebsch_pq * isoClebsch_rs * isoClebschTpq * isoClebschTrs * me3b;
           }
          }
         }
        }
        matelNO2b[istart + offset] = me_no2b;
        offset++;
      }
    }
  }
  if (verbose)std::cout << "Done." << std::endl;

  } // for ch
  // all done.

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
  bool verbose = false;
//  if (ch==2) verbose=true;
  if (verbose) std::cout << "Enter " <<__func__  << " with  ch = " << ch << "  number of threads: " << omp_get_num_threads() << std::endl;
  double t_start = omp_get_wtime();
  TwoBodyChannel& tbc = hf.modelspace->GetTwoBodyChannel(ch);
  int Jab = tbc.J;
  int Tzab = tbc.Tz;
  int Tab_min = std::abs( Tzab );
  int Tab_max = 1;
  int nkets = tbc.GetNumberKets();

  if (verbose) std::cout << "This channel: Jab,parity,Tzab = " << Jab << " " << tbc.parity << " " << Tzab << std::endl;

  std::set<size_t> occupied_orbits;
  for (size_t irow=0; irow<hf.rho.n_rows; irow++)
  {
    if ( arma::norm(hf.rho.row(irow),"fro")>1e-8) occupied_orbits.insert(irow);
  }

  // 3-body states that will contribute to the NO2B matrix elements in this channel. Tab is not conserved by V, so it doesn't go with the sub-block labels
  struct ket3b { int a; int b; int c; int Tab;};

  // hash table (twoJ twoT Eabc) -> {vector of lab-frame 3b kets in that channel}    --- local is probably a bad choice of terminology...
  std::unordered_map<std::array<unsigned short,3>,std::vector<ket3b>,array3_hash> local_3b_kets;
  // hash table (a,b,c,Tab,twoJ,twoT) -> index in the matrix where it lives (need to infer Eabc from a,b,c)
  std::unordered_map<std::array<unsigned short,6>,size_t,array6_hash > ket_3b_lookup;  // hash table for reverse lookups
  // structure for the lab-frame 3b matrix elements. (twoJ twoT parity) -> field (Eabc,Edef) -> matrix of 3-body matrix elements
  std::unordered_map<std::array<unsigned short,3>, arma::field<arma::mat>, array3_hash> V3Full;

  if (verbose) std::cout << "Begin loop enumerating 3-body states" << std::endl;
  double t_internal = omp_get_wtime();
  // Go through all the kets in this 2-body channel and figure out which 3-body states we'll need.
  for (int iket=0; iket<nkets; iket++)
  {
    Ket& ket = tbc.GetKet(iket);
//    int Tzab = (ket.op->tz2 + ket.oq->tz2)/2;
    int a = ket.p;
    int b = ket.q;
    if (verbose) std::cout << "iket = " << iket << "  p,q = " << a << " " << b << std::endl;
    for ( int c : occupied_orbits ) // occupied orbits -> oscillator orbits which end up with a non-zero occupation in the HF reference
    {
      Orbit& oc = hf.modelspace->GetOrbit(c);
      int Eabc = 2*(ket.op->n+ket.oq->n+oc.n)+ket.op->l+ket.oq->l+oc.l ;
      if ( Eabc > E3max ) continue;
      int parity = (ket.op->l + ket.oq->l + oc.l)%2;
      int twoTz = 2*Tzab + oc.tz2;
      int twoT_min = std::abs(twoTz);
      int twoT_max = 3;
      int twoJ_min = std::abs( oc.j2 - 2*Jab );
//      int twoJ_max = oc.j2 + 2*Jab;
      int twoJ_max = std::min( oc.j2 + 2*Jab, twoJmax_cut);
      for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
      {
       for (int twoT=twoT_min; twoT<=twoT_max; twoT+=2)
       {
        for (int Tab=Tab_min; Tab<=Tab_max; Tab++)
        {
         if ( (twoT-2*Tab)>1 ) continue;
         auto& channel_vec = local_3b_kets[ MakeUshort3({twoJ,twoT,Eabc})];  // find the J,T,E channel this 3b state belongs to
         ket_3b_lookup[ MakeUshort6({a,b,c,Tab,twoJ,twoT}) ] = channel_vec.size();  // this tells us where in the matrix this 3b state will live
         channel_vec.push_back( {a,b,c,Tab}); // add the 3b state to the vector of 3b states in this J,T,E channel
//         std::cout << "Added state |" << a << " " << b << " " << c << " " << Tab << ">  to channel  ( " << twoJ << " " << twoT << " " << Eabc << " ) " << std::endl;
         auto V3key = MakeUshort3({twoJ,twoT,parity});
         if ( V3Full.find( V3key ) == V3Full.end())
         {
           // Allocate the arma::field of matrices for this J,T,p channel.
           // We only use ether the even or the odd Eabc, but since we don't allocate the others I don't expect this to be much of an issue.
           V3Full[V3key] = arma::field<arma::mat>(E3max+1,E3max+1); 
         }
        }
       }
      }
    }
    
  }
  if (verbose) std::cout << "...done." << std::endl;
  IMSRGProfiler::timer[std::string(__func__) + "_enumkets"] += omp_get_wtime() - t_internal;

  t_internal = omp_get_wtime();
  if (verbose) std::cout << "Allocate matrices in V3Full" << std::endl;
  // Allocate the <abc,Tab|def,Tde> matrices that live in the (Eabc,Edef) fields that live in the JTp hash table...Yikes.
  for ( auto& iter : V3Full )
  {
    int twoJ = iter.first[0];  int twoT = iter.first[1]; int parity = iter.first[2];
    for (int Eabc=parity; Eabc<=E3max; Eabc+=2)
    {
      size_t dim_abc = local_3b_kets[ MakeUshort3({twoJ,twoT,Eabc})].size();
      for (int Edef=parity; Edef<=E3max; Edef+=2)
      {
        size_t dim_def = local_3b_kets[ MakeUshort3({twoJ,twoT,Edef})].size();
        iter.second(Eabc,Edef).zeros(dim_abc, dim_def);
      }
    }
  }
  if (verbose) std::cout << "...done." << std::endl;

  IMSRGProfiler::timer[std::string(__func__)+"V3FullAllocate"] += omp_get_wtime() - t_internal;

  // Now, we just need those T coefficients. We need a matrix for every lab frame (J T Eabc),  and every jacobi+CM (J12 T12 E12, Lcm), but T12=T so that's redundant

  t_internal = omp_get_wtime();

  std::unordered_map<std::array<unsigned short,6>, arma::mat, array6_hash> Tcoeffs;

  if (verbose) std::cout << "Begin loop allocating and computing Tcoefficients" << std::endl;
  for (int nth_pass=0; nth_pass<=1; nth_pass++)
  {
   if (verbose) std::cout << "  pass " << nth_pass << std::endl;
   #pragma omp parallel if(nth_pass==1)
   {
       size_t count = 0;
       int ithread = omp_get_thread_num();
       int nthreads = omp_get_num_threads();
//    for (auto& iterJTE : local_3b_kets)
    for (auto element = local_3b_kets.begin(); element !=local_3b_kets.end(); ++element, count++)
    {
      if(nth_pass==1 and count%nthreads != ithread) continue; // Check if this thread should compute this element
      int twoJ=element->first[0], twoT=element->first[1], Eabc=element->first[2];
      size_t dim_abc = element->second.size();
      if (verbose) std::cout << "||labket channel JTE " << twoJ << " " << twoT << " " << Eabc << "  dimension " << dim_abc << std::endl;
      for (int E12=0; E12<=std::min(Nmax,Eabc); E12++)
      {
        int Ecm = Eabc-E12;
        for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
        {
          if (verbose) std::cout << "..E12,Ecm,Lcm = " << E12 << " " << Ecm << " " << Lcm << std::endl;
          int Ncm = (Ecm-Lcm)/2;
          int twoJ12_min = std::abs(twoJ-2*Lcm);
//          int twoJ12_max = std::min(twoJ+2*Lcm,twoJmax);
          int twoJ12_max = std::min(twoJ+2*Lcm,twoJmax_cut);
          for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
          {
            if (verbose) std::cout << "checking dimensions for " << twoT << " " << twoJ12 << " " << E12%2 << " " << E12 << std::endl;
            size_t dimAS = GetDimensionAS( twoT, twoJ12, E12%2, E12 );
            size_t dimNAS = GetDimensionNAS( twoT, twoJ12, E12%2, E12 );
            if (dimNAS<1) continue;
            if (nth_pass==0) // on the first pass, we allocate
            {
//              Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12,Lcm}) ] = arma::mat( dim_abc, dimNAS, arma::fill::zeros);
              Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12,Lcm}) ] = arma::mat( dim_abc, dimAS, arma::fill::zeros);
//              if (verbose) std::cout << "allocated a tcoeff matrix for J,T,Eabc J12,E12,Lcm = " << twoJ << " " << twoT << " " << Eabc << "  " << twoJ12 << " " << E12 << " " << Lcm << "  with size " << dim_abc << " x " << dimNAS << std::endl;
              if (verbose) std::cout << "allocated a tcoeff matrix for J,T,Eabc J12,E12,Lcm = " << twoJ << " " << twoT << " " << Eabc << "  " << twoJ12 << " " << E12 << " " << Lcm << "  with size " << dim_abc << " x " << dimAS << std::endl;
            }
            if (nth_pass==1) // on the second pass, we compute all the matrix elements
            {
//              auto& TcoeffMat = Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12,Lcm}) ];
              arma::mat T_mat_tmp( dim_abc, dimNAS, arma::fill::zeros);
              for (size_t iNAS=0; iNAS<dimNAS; iNAS++)
              {
                jacobi1_state jac1;
                jacobi2_state jac2;
                GetJacobiStates( twoT, twoJ12, E12%2, E12, iNAS, jac1, jac2);
                for (size_t iabc=0; iabc<dim_abc; iabc++)
                {
                  auto& ket_abc = element->second[iabc];
                  Orbit& oa = hf.modelspace->GetOrbit(ket_abc.a);
                  Orbit& ob = hf.modelspace->GetOrbit(ket_abc.b);
                  Orbit& oc = hf.modelspace->GetOrbit(ket_abc.c);
                  if ( ket_abc.Tab != jac1.t ) continue;
                  if ( 2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l > E3max) continue;
                  double tcoeff = ComputeTcoeff( oa.n, oa.l, oa.j2, ob.n, ob.l, ob.j2, oc.n, oc.l, oc.j2, Jab, twoJ, jac1.n, jac1.l, jac1.s, jac1.j, jac2.n, jac2.l, jac2.j2, twoJ12, Ncm, Lcm);
//                  TcoeffMat(iabc,iNAS) = tcoeff;
                  T_mat_tmp(iabc,iNAS) = tcoeff;
                }
              }
             //TODO CFPs here 
              size_t cfp_begin = GetCFPStartLocation(twoT,twoJ12,E12);
              arma::mat cfp_mat( &(cfpvec[cfp_begin]), dimNAS, dimAS, false); // false refers to copy_aux_mem
              Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12,Lcm}) ] = sqrt(6) * T_mat_tmp * cfp_mat;
            }
          }
        }
      }
    }
   }
  IMSRGProfiler::timer[std::string(__func__)+"_Tcoeffs_"+std::to_string(nth_pass)] += omp_get_wtime() - t_internal;
  t_internal = omp_get_wtime();
  }
  if (verbose) std::cout << "...done" <<std::endl;

//  IMSRGProfiler::timer[std::string(__func__)+"computeTcoeffs"] += omp_get_wtime() - t_internal;
//  t_internal = omp_get_wtime();

  int Eabc_min = std::max( Jab-1, 0); // for Jab > 1, we can't do it with  Eabc=0, so don't bother worrying about that


   if (verbose) std::cout << std::endl << "===============================================" << std::endl
                          << "start loops doing the mat mults. Eabc_min = " << Eabc_min << "based on Jab = " << Jab << std::endl;

  // Next we do the mat mult to perform the transformation
//   for (int Ecm=Eabc_min; Ecm<=E3max; Ecm++)
//   # pragma omp parallel for collapse(2) schedule(dynamic,1) 
   for (int Ecm=0; Ecm<=E3max; Ecm++)
   {
//    for (int Lcm_step=0; Lcm_step<=E3max; Lcm_step+=2)
    for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
    {
//      int Lcm = Lcm_step + Ecm%2;
//      if (Lcm>Ecm) continue;
      if (verbose) std::cout << " Ecm,Lcm = " << Ecm << " " << Lcm << std::endl;
      int Ncm = (Ecm-Lcm)/2;
//      if (verbose) std::cout << "J12 range " << std::max(1,2*Jab-2*emax-1) << " to " << std::min(2*Jab+2*emax+1,twoJmax) << "  from Jab=" << Jab << " emax = " << emax << "  twoJmax = " << twoJmax << std::endl;
      if (verbose) std::cout << "J12 range " << std::max(1,2*Jab-2*emax-1-2*Lcm) << " to " << std::min(2*Jab+2*emax+1+2*Lcm,twoJmax) << "  from Jab=" << Jab << " emax = " << emax << "  twoJmax = " << twoJmax << "  and Lcm = " << Lcm << std::endl;
//      for (int twoJ12=std::max(1,2*Jab-2*emax-1); twoJ12<=std::min(2*Jab+2*emax+1,twoJmax); twoJ12+=2)
      for (int twoJ12=std::max(1,2*Jab-2*emax-1-2*Lcm); twoJ12<=std::min(2*Jab+2*emax+1+2*Lcm,twoJmax); twoJ12+=2)
      {
//       for (int twoT=1; twoT<=2*std::abs(Tzab)+1; twoT+=2)
       for (int twoT=1; twoT<=3; twoT+=2)
       {
//        if (ch==1 and twoT==3) verbose=true;
//        else verbose=false;
        for (int E12abc=0; E12abc<=std::min(Nmax,E3max-Ecm); E12abc++)
        {
          if (verbose) std::cout << "~~  Ecm,Lcm,J12,T,E12abc  " << Ecm << " " << Lcm << " " << twoJ12 << " " << twoT << " " << E12abc << std::endl;
          size_t dimNAS_abc = GetDimensionNAS( twoT, twoJ12, E12abc%2, E12abc ); 
          size_t dimAS_abc = GetDimensionAS( twoT, twoJ12, E12abc%2, E12abc ); 
          if (dimNAS_abc==0 or dimAS_abc==0) continue;
//          size_t cfp_begin_abc = GetCFPStartLocation(twoT,twoJ12,E12abc);
//          arma::mat cfp_abc( &(cfpvec[cfp_begin_abc]), dimNAS_abc, dimAS_abc, false); // false refers to copy_aux_mem

         
//          for (int E12def=E12abc%2; E12def<=std::min(Nmax,E3max-Ecm); E12def+=2 ) // we can probably make use of some hermiticity here
          for (int E12def=E12abc; E12def<=std::min(Nmax,E3max-Ecm); E12def+=2 ) // we can probably make use of some hermiticity here
          {
            size_t dimNAS_def = GetDimensionNAS( twoT, twoJ12, E12def%2, E12def ); 
            size_t dimAS_def = GetDimensionAS( twoT, twoJ12, E12def%2, E12def ); 
            if (dimNAS_def==0 or dimAS_def==0) continue;
//            size_t cfp_begin_def = GetCFPStartLocation(twoT,twoJ12,E12def);
//            arma::mat cfp_def( &(cfpvec[cfp_begin_def]), dimNAS_def, dimAS_def, false); // false refers to copy_aux_mem

            size_t startlocAS = GetStartLocAS( twoT, twoJ12, E12abc, E12def);
            arma::mat matelAS( &meAS[startlocAS], dimAS_abc, dimAS_def, false ); 
//            arma::mat local_copy = Tabc * matelAS * cfp_def.t() * Tdef.t(); 

//            if (verbose) std::cout << "matelNAS = cfp * matelAS * cfp" << std::endl;
//            if (verbose) std::cout << "matelAS = " << std::endl << matelAS << std::endl;
//            arma::mat matelNAS = 6*cfp_abc * matelAS * cfp_def.t();
            if (verbose) std::cout << "...done" << std::endl;

            // Now we loop over the lab-frame matrix elements and update them with the contribution from this jacobi channel
            int twoJ_min = std::max( std::abs(twoJ12-2*Lcm), 2*Jab-2*emax-1);
            int twoJ_max = std::min( twoJ12+2*Lcm, 2*Jab+2*emax+1 );
            if (verbose) std::cout << "running twoJ from " << twoJ_min << " to " << twoJ_max << "  based on Jab = " << Jab << " J12 = " << twoJ12 << " emax= " << emax << "  Lcm = " << Lcm << std::endl;
            int lab_parity = (E12abc+Ecm)%2;
//           # pragma omp parallel for schedule(dynamic,1)
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
            {
              auto& V3field = V3Full[MakeUshort3({twoJ,twoT,lab_parity})];
//              for (int Eabc=Eabc_min+(Eabc_min+lab_parity)%2; Eabc<=E3max; Eabc+=2)
              int Eabc = E12abc + Ecm;
              int Edef = E12def + Ecm;
//              {
//                if (std::abs(E12abc-Lcm)>Eabc or (E12abc+Lcm)<Eabc) continue;
                size_t dim_abc = local_3b_kets[ MakeUshort3({twoJ,twoT,Eabc})].size();
                if (dim_abc<1) continue;
//                auto Tabc = Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12abc,Lcm}) ];
//                if (verbose) std::cout << " matelNAS_abc = T * matelNAS Eabc,E12abc,Ecm,Lcm: " << Eabc << " " << E12abc << " " << Ecm << " " << Lcm << "  twoJ, twoJ12,twoT = " << twoJ << " " << twoJ12 << " " << twoT << " Jab = " << Jab << "  dim_abc = " << dim_abc << "   tbc.parity = " << tbc.parity << std::endl;
//                arma::mat matelNAS_abc =  Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12abc,Lcm}) ] * matelNAS;
                arma::mat matelAS_abc =  Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12abc,Lcm}) ] * matelAS;
//                if (verbose) std::cout << " ...done" << std::endl;

//                if (verbose and twoJ==1 and twoJ12==1 and Ecm==0 and E12abc==1 and E12def==1 and twoT==3 and Lcm==0)
//                if (verbose and twoJ==3 and Eabc==1 and Edef==1 and twoT==3 )
//                {
//                  std::cout << "  matmult abc: Tabc " << std::endl << Tcoeffs[ MakeUshort6({twoJ,twoT,Eabc, twoJ12,E12abc,Lcm}) ] << std::endl << "matelNAS " << std::endl << matelNAS << std::endl;
//                }

//                for (int Edef=Eabc; Edef<=E3max; Edef+=2)
//                {
//                  if ( std::abs(E12def-Lcm)>Edef or (E12def+Lcm)<Edef) continue;
                  size_t dim_def = local_3b_kets[ MakeUshort3({twoJ,twoT,Edef})].size();
                  if (dim_def<1) continue;
//                  if (verbose) std::cout << "  local_lab_mat = matelNAS_abc * T   Edef,E12def: " << Edef << " " << E12def << std::endl;
//                  arma::mat local_lab_mat = matelNAS_abc * Tcoeffs[ MakeUshort6({twoJ,twoT,Edef, twoJ12,E12def,Lcm}) ].t();
                  arma::mat local_lab_mat = matelAS_abc * Tcoeffs[ MakeUshort6({twoJ,twoT,Edef, twoJ12,E12def,Lcm}) ].t();
//                  if (verbose) std::cout << "  ...done  lab_mat dimensions: " << local_lab_mat.n_rows << " x " << local_lab_mat.n_cols
//                            << "   lab_parity = " << lab_parity << "  E12abc,Ecm" << E12abc << " " << Ecm  << std::endl;
//                 if (verbose)
//                if (verbose and twoJ==1 and twoJ12==1 and Ecm==0 and E12abc==1 and E12def==1 and twoT==3 and Lcm==0)
                if (verbose and twoJ==3 and Eabc==1 and Edef==1 and twoT==3)
                 {
//                   std::cout << " local_lab_mat.  matelNAS_abc = " << std::endl << matelNAS_abc << std::endl << "Tdef = " << std::endl
                   std::cout << " local_lab_mat.  matelNAS_abc = " << std::endl << matelAS_abc << std::endl << "Tdef = " << std::endl
                                        << Tcoeffs[ MakeUshort6({twoJ,twoT,Edef, twoJ12,E12def,Lcm}) ] << std::endl << "local_lab_mat = " << std::endl
                                        << local_lab_mat << std::endl;
                 }
                  #pragma omp critical
                  {
                    V3field(Eabc,Edef) += local_lab_mat;
                    if (Eabc != Edef)   V3field(Edef,Eabc) += local_lab_mat.t();
//                    if (Eabc > Edef)   V3field(Edef,Eabc) += local_lab_mat.t();
                  }
//                }// for Edef
//               } // for Eabc
              } // for twoJ
            } // for E12def
          } // for E12abc
        } // for twoT
      } // for twoJ12
    } // for Lcm
  } // for Ecm
  if (verbose) std::cout << "...done" << std::endl;

  IMSRGProfiler::timer[std::string(__func__)+"_matmult"] += omp_get_wtime() - t_internal;

  t_internal = omp_get_wtime();

  // and finally, we compute the NO2B part.

   if (verbose) std::cout << "Start loop over nkets, computing the NO2B part "<< std::endl;
   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0; i<nkets; ++i)    
   {
      Ket & bra = tbc.GetKet(i);
      int a = bra.p;
      int b = bra.q;
      int e2bra = 2*bra.op->n + bra.op->l + 2*bra.oq->n + bra.oq->l;
      for (int j=0; j<nkets; ++j)
      {
        if (i>j) continue;
        Ket & ket = tbc.GetKet(j); 
        int d = ket.p;
        int e = ket.q;
        int e2ket = 2*ket.op->n + ket.op->l + 2*ket.oq->n + ket.oq->l;
//        for ( auto a : modelspace->all_orbits )
        for ( int c : occupied_orbits )
        {
          Orbit & oc = hf.modelspace->GetOrbit(c);
          int Eabc = 2*oc.n+oc.l+e2bra;
          if ( Eabc > E3max ) continue;
          int twoTz = ket.op->tz2 + ket.oq->tz2 + oc.tz2;
          for (int f : hf.modelspace->OneBodyChannels.at({oc.l,oc.j2,oc.tz2}))
          {
            Orbit & of = hf.modelspace->GetOrbit(f);
            int Edef = 2*of.n+of.l+e2ket;
            if ( (Eabc+Edef)%2>0) continue;
            if ( Edef > E3max ) continue;
            if ( std::abs(hf.rho(c,f)) < 1e-8 ) continue; // Turns out this helps a bit (factor of 5 speed up in tests)
            int twoJ_min = std::abs(2*Jab-oc.j2);
            int twoJ_max = 2*Jab + oc.j2;
            if (verbose) std::cout << "Now running J from " << twoJ_min << " to " << twoJ_max << " based on Jab = " << Jab << " jc = " << oc.j2 << std::endl;
            for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
            {
//              double v3_pn = 0;
              for (int twoT=std::abs(twoTz); twoT<=3; twoT+=2)
              {
              double v3_pn = 0;
                for (int Tab=Tab_min; Tab<=Tab_max; Tab++)
                {
                  double iso_clebsch_abc =  AngMom::CG(0.5,0.5*bra.op->tz2,0.5,0.5*bra.oq->tz2,Tab,Tzab)
                                          * AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2,0.5*twoT,Tzab+0.5*oc.tz2);
                  if (std::abs(iso_clebsch_abc)<1e-9) continue;
                  auto lab_abc_key = MakeUshort6({a,b,c,Tab,twoJ,twoT});
//                  if ( ket_3b_lookup.find(lab_abc_key) == ket_3b_lookup.end() )  std::cout << "TROUBLE!!! abc ..." << a << " " << b << " " << c << " " << Tab << " " << twoJ << " " << twoT << std::endl;
                  size_t iabc = ket_3b_lookup.at(lab_abc_key);
//                  size_t iabc = ket_3b_lookup[MakeUshort6({a,b,c,Tab,twoJ,twoT})];
                  for ( int Tde=Tab_min; Tde<=Tab_max; Tde++)
                  {
                    double iso_clebsch_def =  AngMom::CG(0.5,0.5*ket.op->tz2,0.5,0.5*ket.oq->tz2,Tde,Tzab)
                                            * AngMom::CG(Tde,Tzab,0.5,0.5*of.tz2,0.5*twoT,Tzab+0.5*oc.tz2);
                    if (std::abs(iso_clebsch_def)<1e-9) continue;
                    auto lab_def_key = MakeUshort6({d,e,f,Tde,twoJ,twoT});
//                    if ( ket_3b_lookup.find(lab_def_key) == ket_3b_lookup.end() )  std::cout << "TROUBLE!!! def ..." << d << " " << e << " " << f << " " << Tde << " " << twoJ << " " << twoT << std::endl;
                    size_t idef = ket_3b_lookup.at(lab_def_key);
//                    size_t idef = ket_3b_lookup[MakeUshort6({d,e,f,Tde,twoJ,twoT})];
                    auto V3key = MakeUshort3({twoJ,twoT,Eabc%2});
//                    if ( V3Full.find(V3key) == V3Full.end() ) std::cout << "TROUBLE!!! LOOKING FOR " << twoJ << " " << twoT << " " << Eabc%2 << "   and came up empty" << std::endl;
                    v3_pn += iso_clebsch_abc * iso_clebsch_def *  V3Full.at(V3key)(Eabc,Edef)(iabc,idef);
//                    v3_pn += iso_clebsch_abc * iso_clebsch_def *  V3Full[MakeUshort3({twoJ,twoT,Eabc%2})](Eabc,Edef)(iabc,idef);
                    if (verbose ) std::cout << "J, T, Tab, Tde = " << twoJ << " " << twoT << " " << Tab << " " << Tde << " clebsch " << iso_clebsch_abc << " " << iso_clebsch_def
                                           << "  abcdef = " << a << " " << b << " " << c << " " << d << " " << e << " " << f << "   "
                              << V3Full[V3key](Eabc,Edef)(iabc,idef) << "   -> v3pn = " << v3_pn << std::endl;
//                              << V3Full[MakeUshort3({twoJ,twoT,Eabc%2})](Eabc,Edef)(iabc,idef) << "   -> v3pn = " << v3_pn << std::endl;
                  }
                }
//                if (verbose) std::cout << "T = " << twoT << "  J = " << twoJ <<  "   abcdef = " << a << " " << b << " " << c << " " << d << " " << e << " " << f
//                                       << "   v3_pn = " << v3_pn << std::endl;
              V3NO(i,j) += hf.rho(c,f) * (twoJ+1) * v3_pn;
              if (verbose) std::cout << "V3NO += " << hf.rho(c,f) << " * " << twoJ+1 << " * " << v3_pn << "   -> V3NO =  " << V3NO(i,j) << std::endl << std::endl;;
              }
              
//              V3NO(i,j) += hf.rho(c,f) * (twoJ+1) * v3_pn;

//              if (verbose) std::cout << "V3NO += " << hf.rho(c,f) << " * " << twoJ+1 << " * " << v3_pn << "   -> V3NO =  " << V3NO(i,j) << std::endl;
            }
          }
        }
//        V3NO(i,j) /= (2*Jab+1);
        if (bra.p==bra.q)  V3NO(i,j) /= SQRT2; 
        if (ket.p==ket.q)  V3NO(i,j) /= SQRT2; 
        V3NO(j,i) = V3NO(i,j);
      }
   }
   if (verbose) std::cout << "...done "<< std::endl;

  V3NO /= 2*Jab+1;
//  std::cout << "Done with NO2B loop" << std::endl;

  IMSRGProfiler::timer[std::string(__func__)+"_no2bLoop"] += omp_get_wtime() - t_internal;

//  IMSRGProfiler::timer[std::string(__func__)+"_computeNO2B"] += omp_get_wtime() - t_internal;
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
//  IMSRGProfiler::counter[std::string(__func__)+"matmultAS"] += mat_mult_AS;
//  IMSRGProfiler::timer[std::string(__func__)+"LcmLoop"] += LcmLoopTime;
//  IMSRGProfiler::timer[std::string(__func__)+"J12Loop"] += J12LoopTime;
}




/*
void Jacobi3BME::GetNO2b_single_channel( HartreeFock& hf, int ch, arma::mat& V3NO )
{
  std::cout << "Enter " <<__func__  << " with  ch = " << ch << "  number of threads: " << omp_get_num_threads() << std::endl;
  double t_start = omp_get_wtime();
  TwoBodyChannel& tbc = hf.modelspace->GetTwoBodyChannel(ch);
  int Jab = tbc.J;
  int nkets = tbc.GetNumberKets();

//  std::cout << "size of V3NO = " << V3NO.n_rows << " x " << V3NO.n_cols  << ".  nkets = " << nkets << ".  size of rho = " << hf.rho.n_rows << " x " << hf.rho.n_cols << std::endl;
  V3NO.zeros( nkets,nkets);

//  if (ch==6) std::cout << " Jab = " << Jab << std::endl;

  std::set<size_t> occupied_orbits;
  for (size_t irow=0; irow<hf.rho.n_rows; irow++)
  {
    if ( arma::norm(hf.rho.row(irow),"fro")>1e-8) occupied_orbits.insert(irow);
  }


  double t_internal = omp_get_wtime();

  std::unordered_map<std::array<unsigned short,7>,arma::mat,array7_hash> TcoeffTable;
  std::unordered_map<std::array<unsigned short,7>,bool,array7_hash> TcoeffSkip;
//  std::unordered_map<std::string,arma::mat> TcoeffTable;
//  std::unordered_map<std::string,bool> TcoeffSkip;

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
     for (int c : occupied_orbits )
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
//               for (int Lcm=std::max(Ecm%2,std::abs(twoJ-twoJ12)/2); Lcm<=std::min(Ecm,(twoJ+twoJ12)/2); Lcm+=2)
               for (int Lcm=Ecm%2; Lcm<=std::min(Ecm,(twoJ+twoJ12)/2); Lcm+=2)
               {
//                 if ( std::abs(twoJ-twoJ12)>2*Lcm or (twoJ+twoJ12)<2*Lcm) continue;
                 if ( std::abs(twoJ-twoJ12)>2*Lcm ) continue;
                 int Ncm=(Ecm-Lcm)/2;
//                 std::ostringstream oss;
//                 oss << ibra << " " << c << " " << twoJ << " " << Ecm << " " << twoJ12 << " " << twoT << " " << Lcm;
//                 std::string tcoeff_hash = oss.str();
//                 std::array<unsigned short,7> tcoeff_hash = {ibra,c,twoJ,Ecm,twoJ12,twoT,Lcm};
                 auto tcoeff_hash = MakeUshort7({ibra,c,twoJ,Ecm,twoJ12,twoT,Lcm});

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


  long int mat_mult_AS = 0;
  double LcmLoopTime = 0;
  double J12LoopTime = 0;
  // TODO: It may be more effective to make a matrix of 3-body matrix elements (organized into the appropriate J,pi,T,Eabc,Edef blocks), compute those, and then get the NO2B piece

  #pragma omp parallel for schedule(dynamic,1) reduction(+:mat_mult_AS,LcmLoopTime,J12LoopTime)
  for (int ibra=0; ibra<nkets; ibra++)
  {
    Ket& bra = tbc.GetKet(ibra);
    if(ch==1) std::cout << "ibra = " << ibra << "  pq = " << bra.p << " " << bra.q << std::endl;
    int a = bra.p;
    int b = bra.q;
    Orbit& oa = hf.modelspace->GetOrbit(a);
    Orbit& ob = hf.modelspace->GetOrbit(b);
    if ( 2*(oa.n+ob.n)+oa.l+ob.l > E2max  ) continue;
    int Tzab = (oa.tz2+ob.tz2)/2;

    for (int iket=0; iket<=ibra; iket++ )
    {
      Ket& ket = tbc.GetKet(iket);
//      if(ch==12)std::cout << " iket = " << iket << std::endl;
      int d = ket.p;
      int e = ket.q;
      Orbit& od = hf.modelspace->GetOrbit(d);
      Orbit& oe = hf.modelspace->GetOrbit(e);
      if ( 2*(od.n+oe.n)+od.l+oe.l > E2max  ) continue;


//     for (auto c : hf.modelspace->holes )
     for (int c : occupied_orbits )
     {
//      if (ch==12) std::cout << "  c = " << c << std::endl;
      Orbit& oc = hf.modelspace->GetOrbit(c);
      int twoJ_min = std::abs( 2*Jab - oc.j2 );
      int twoJ_max = (2*Jab+oc.j2);
      int twoTz = oa.tz2 + ob.tz2 + oc.tz2;
      int twoT_min = std::abs(twoTz);
      if ( 2*(oa.n+oc.n)+oa.l+oc.l > E2max  ) continue;
      if ( 2*(oc.n+ob.n)+oc.l+ob.l > E2max  ) continue;
      int Eabc = 2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l;
      if (Eabc>E3max) continue;

//     if (ch==12)
//     {
//       std::cout << " Looping over f states :  ";
//       for (auto h : hf.modelspace->OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) ) std::cout << h << " ";
//       std::cout << std::endl;
//     }
       for (int f : hf.modelspace->OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) )
       {
//         std::cout << "   f = " << f << std::endl;
//         if (ch==12 and ibra==0 and iket==0) std::cout << " Jab, bra,ket  abcdef  " << Jab << " " << ibra << " " << iket << "  " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
         if ( std::abs( hf.rho(c,f)) <1e-8) continue;
         Orbit& of = hf.modelspace->GetOrbit(f);
         if ( 2*(od.n+of.n)+od.l+of.l > E2max  ) continue;
         if ( 2*(of.n+oe.n)+of.l+oe.l > E2max  ) continue;
         int Edef = 2*(od.n+oe.n+of.n)+od.l+oe.l+of.l;
         if (Edef>E3max) continue;

         double v_no2b_cf = 0;

         // Inner loops over jacobi stuff

//           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
//           {
//            double vsum_J = 0;
            for ( int twoT=twoT_min; twoT<=3; twoT+=2)
            {
                int rows = 2-(twoT/2); // T=1/2 -> Tab=0,1   T=3/2 -> Tab=1.

                arma::mat isospin_mat_abc( rows, rows, arma::fill::eye );
                arma::mat isospin_mat_def( rows, rows, arma::fill::eye );

         for (int Ecm=0; Ecm<=std::min(Eabc,Edef); Ecm++)
         {
           int E12abc = Eabc - Ecm;
           int E12def = Edef - Ecm;
           if (E12abc>Nmax or E12def>Nmax) continue;
           int parity = E12abc%2;

//            int twoJ12_min = std::max( 1, (twoJ - 2*Ecm) ); // maybe take a second look at these limits...
            int twoJ12_min = std::max( 1, (twoJ_min - 2*Ecm) ); // maybe take a second look at these limits...
            int twoJ12_max = std::min( twoJmax, twoJ_max +2*Ecm);
//            if ( ch==12 and ibra==0 and iket==0 and twoJ==3 and twoT==1 and c==1 and f==1)
//            {
//             std::cout << "Ecm = " << Ecm << "   twoJ12_min/max = " << twoJ12_min << " " << twoJ12_max << std::endl;
//            }

                for (int Tab=twoT/2; Tab<=1; Tab++)
                {   // because c and f are in the same one-body channel,  Tzab = Tzde, although Tab need not be Tde, and tza need not be tzd, etc.
                  isospin_mat_abc(Tab-twoT/2,Tab-twoT/2) = AngMom::CG(0.5,0.5*oa.tz2,0.5,0.5*ob.tz2, Tab, Tzab) * AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2, 0.5*twoT, 0.5*twoTz) ;
                  isospin_mat_def(Tab-twoT/2,Tab-twoT/2) = AngMom::CG(0.5,0.5*od.tz2,0.5,0.5*oe.tz2, Tab, Tzab) * AngMom::CG(Tab,Tzab,0.5,0.5*oc.tz2, 0.5*twoT, 0.5*twoTz) ;
                }

              double t_J12_loop = omp_get_wtime();;
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
                
                arma::mat cfp_abc( &(cfpvec[cfp_begin_abc]), dimNAS_abc, dimAS_abc,  false);
                arma::mat cfp_def( &(cfpvec[cfp_begin_def]), dimNAS_def, dimAS_def,  false);
    
                mat_mult_AS++;
                arma::mat matelNAS = 6 * cfp_abc * matelAS * cfp_def.t(); // Compute the non-antisymmetrized matrix elements 


//            if ( ch==12 and ibra==0 and iket==0 and twoJ==3 and twoT==1 and c==1 and f==1)
//            {
//              std::cout << "Lcm will run from " << std::max(Ecm%2,std::abs(twoJ-twoJ12)/2) << "  to " << std::min(Ecm,(twoJ+twoJ12)/2) << "  in steps of 2" << std::endl;
//            }
//
           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
           {
            double vsum_J = 0;

                double t_LCM_loop = omp_get_wtime(); 
//                for (int Lcm=std::max(Ecm%2,std::abs(twoJ-twoJ12)/2); Lcm<=std::min(Ecm,(twoJ+twoJ12)/2); Lcm+=2)
                for (int Lcm=Ecm%2; Lcm<=std::min(Ecm,(twoJ+twoJ12)/2); Lcm+=2)
                {
                  if (2*Lcm < std::abs(twoJ-twoJ12)) continue;
//                  if (ch==6) std::cout << "Ecm = " << Ecm << "  Lcm  = " << Lcm << "  J,J12 " << twoJ << " " << twoJ12 << "   dimAS,dimNAS " << dimAS_abc << " " << dimAS_def << ", " << dimNAS_abc << " " << dimNAS_def << std::endl;
                  int Ncm=(Ecm-Lcm)/2;
                  int rows = 2-(twoT/2); // T=1/2 -> Tab=0,1   T=3/2 -> Tab=1.
//                  arma::mat Tabc( dimNAS_abc, rows, arma::fill::zeros ); 
//                  arma::mat Tdef( dimNAS_def, rows, arma::fill::zeros );

//                  std::ostringstream oss_abc, oss_def;
//                  oss_abc << ibra << " " << c << " " << twoJ << " " << Ecm << " " << twoJ12 << " " << twoT << " " << Lcm;
//                  oss_def << iket << " " << f << " " << twoJ << " " << Ecm << " " << twoJ12 << " " << twoT << " " << Lcm;
//                  std::string t_hash_abc = oss_abc.str();
//                  std::string t_hash_def = oss_def.str();
                  auto t_hash_abc = MakeUshort7({ibra,c,twoJ,Ecm,twoJ12,twoT,Lcm});
                  auto t_hash_def = MakeUshort7({iket,f,twoJ,Ecm,twoJ12,twoT,Lcm});
                  if (  TcoeffSkip[t_hash_abc] or TcoeffSkip[t_hash_def]) continue;  // looks like this maybe helps a little

                  arma::mat& Tabc = TcoeffTable[ t_hash_abc ];
                  arma::mat& Tdef = TcoeffTable[ t_hash_def ];

//                  if (ch==6) std::cout << "t_hash_abc = " << t_hash_abc << "    ,  t_hash_def = " << t_hash_def << std::endl;
//                  if (ch==6) std::cout << "About to multiply. dim Tabc " << Tabc.n_rows << "x" << Tabc.n_cols << "   Tdef " << Tdef.n_rows << "x" << Tdef.n_cols << std::endl;
                  
                  arma::mat result =  isospin_mat_abc * Tabc.t() * matelNAS * Tdef * isospin_mat_def  ;
//                  if (ch==6) std::cout << "Im ok. result = " << result[0] << std::endl;

                  vsum_J += arma::accu( result ) ;

//                  if(ch==12 and ibra==0 and iket==0 and twoJ==3 and twoT==1 and c==1 and f==1) std::cout << "J,J12,Ecm,E12abc,E12def,T,Lcm = " << twoJ << " " << twoJ12  << " " << Ecm << " " << E12abc << " " << E12def << " " << twoT << " " << Lcm 
//                  if(ch==1 and ibra==0 and iket==0 and twoJ==1 and twoT==1 and c==0 and f==0) std::cout << "J,J12,Ecm,E12abc,E12def,T,Lcm = " << twoJ << " " << twoJ12  << " " << Ecm << " " << E12abc << " " << E12def << " " << twoT << " " << Lcm 
//                  if(ch==0 and ibra==1 and iket==0 and twoJ==1 and twoT==1 and c==1 and f==1) std::cout << "J,J12,Ecm,E12abc,E12def,T,Lcm = " << twoJ << " " << twoJ12  << " " << Ecm << " " << E12abc << " " << E12def << " " << twoT << " " << Lcm 
//                                     << "     matrices: "
//                                      << std::endl << Tabc.t() << std::endl << matelNAS << std::endl << Tdef << std::endl << result << std::endl
//                                      << " isospin mats" << std::endl << isospin_mat_abc << std::endl << isospin_mat_def <<std::endl
//                                      << " mNAS * Tdef " << std::endl << matelNAS * Tdef << std::endl
//                                      << " result : " << std::endl << result << std::endl << "  vsum_J = " << vsum_J << std::endl << std::endl;
                  
//                  std::cout << "     bra,ket: " << t_hash_abc << " | " << t_hash_def << "   :  " << arma::accu( result ) << std::endl;
//                  if (  arma::norm( result, "fro") < 1e-7 )
//                  {
//                     std::cout << Tabc << std::endl << Tdef << std::endl;
//                  }

                } // for Lcm
                LcmLoopTime += omp_get_wtime() - t_LCM_loop;

            v_no2b_cf +=  vsum_J * (twoJ+1);
           } // for twoJ
              } // for twoJ12
              J12LoopTime += omp_get_wtime() - t_J12_loop;
          } // for Ecm
//            if(ch==12 and ibra==0 and iket==0 and twoJ==3 and twoT==1) std::cout << "  intermediate +=" << vsum_J  
//                                              << "  ibra,iket,twoJ,twoT,c,f = " << ibra << " " << iket << " "  << twoJ << " " << twoT << " " << c << " " << f  << std::endl;
             } // for twoT
            
//            v_no2b_cf +=  vsum_J * (twoJ+1);
//            if(ch==12 and ibra==0 and iket==0) std::cout << "  v_no2b_cf +=" << vsum_J  << " * " << twoJ+1 << "  -> " << v_no2b_cf
//                                              << "  ibra,iket,twoJ,c,f = " << ibra << " " << iket << " "  << twoJ << " " << c << " " << f  << std::endl;
//           } // for twoJ
//         V3NO(ibra,iket) += v_no2b_cf * hf.rho(c,f); 
         V3NO(ibra,iket) += v_no2b_cf * hf.rho(c,f); 
//         if(ch==0 and ibra==1 and iket==0) std::cout << "   V3NO += v_no2b_cf (" << v_no2b_cf << ") * rho (" << hf.rho(c,f) << ")  -> V3NO " << V3NO(ibra,iket) << std::endl;
       } // for f
     } // for c
     if (bra.p==bra.q) V3NO(ibra,iket) /= SQRT2;
     if (ket.p==ket.q) V3NO(ibra,iket) /= SQRT2;
     V3NO(iket,ibra) = V3NO(ibra,iket);
    } // for iket
  } // for ibra




  V3NO /= 2*Jab+1;
//  std::cout << "Done with NO2B loop" << std::endl;


  IMSRGProfiler::timer[std::string(__func__)+"_computeNO2B"] += omp_get_wtime() - t_internal;
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  IMSRGProfiler::counter[std::string(__func__)+"matmultAS"] += mat_mult_AS;
  IMSRGProfiler::timer[std::string(__func__)+"LcmLoop"] += LcmLoopTime;
  IMSRGProfiler::timer[std::string(__func__)+"J12Loop"] += J12LoopTime;
}

*/







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
//    double timestamp = omp_get_wtime();
    double ja = 0.5*j2a;
    double sa = 0.5;
    double jb = 0.5*j2b;
    double sb = 0.5;
    double jc = 0.5*j2c;
    double sc = 0.5;
    double S2 = 0.5; 
//    std::cout << std::endl << " " << __func__ << " ( " << na << " " << la << " " << j2a << "  " << nb << " " << lb << " " << j2b << "  " << nc << " " << lc << " " << j2c
//              << "   " << Jab << " " << twoJ << "  " << N1 << " " << L1 << " " << S1 << "   " << N2 << " " << L2 << "  " << twoJ2 << " " << twoJ12 << "  " << Ncm << " " << Lcm << " ) " << std::endl;

    if ( (la+lb+lc + L1+L2+Lcm)%2 >0 ) return 0; // parity conservation
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
      int curlyL_min =  std::max( std::abs(L1-Lab), std::abs(Jab-J1) ) ;
      int curlyL_max =  std::min( L1+Lab,  Jab+J1) ;
//      int curlyL_min =  std::abs(L1-Lab) ;
//      int curlyL_max =  L1+Lab ;
      if ( (curlyL_min + L1 + la + lb)%2>0 ) curlyL_min++; // parity for the first Moshinsky bracket
//      if ( (curlyL_min + lc + Lcm + L2)%2>0 ) continue; // parity for the second Moshinsky bracket
      for (int curlyL=curlyL_min; curlyL<=curlyL_max; curlyL+=2)
      {
//          if ( std::abs(curlyL-J1)>Jab   ) continue;
//          if ( (curlyL+J1)<Jab  ) continue;
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
//          if (Lambda>E3max) continue; // this is tentative...
          // It will speed things up a bit to precompute these
          int moshphase2 = AngMom::phase( Lcm+L2+Lambda + 0*(Lcm+L2 - curlyL-lc)/2 ); // converting phase conventions (comment at the beginning of the function)
          double mosh2 = moshphase2 * (2*Lambda+1) * GetMoshinsky2(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda);
//          double mosh2 = moshphase2 * (2*Lambda+1) * AngMom::Moshinsky(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda,  2.0);
          if (std::abs(mosh2)<1e-9) continue;

          // first, check some triangle conditions that will kill the 12j. TODO These should really be incorporated in the summation ranges...
//          if ( std::abs(curlyL-J1)>Jab  or std::abs(Lambda-curlyL)>lc or std::abs(Lcm-Lambda)>L2 ) continue;
//          if ( (curlyL+J1)<Jab or (Lambda+curlyL)<lc or (Lcm+Lambda)<L2 ) continue;
          double sum_L = 0;

//   Maybe try inserting the 12j stuff here instead of doin the loops
//          int phase_L = AngMom::phase( lc + Lab + (twoS12+twoJ+1-twoJ2)/2 + S1 + Jab);
//          double sixjL = 0

          // { J     jc     sc     J2     }
          // {   Jab    lc      L2   J12  } = sum_x (-1)^(S-x) { J      jc  Jab } { jc      sc   lc } { sc    J2   L2 } { J2  J1  J12 }
          // { J1  curlyL Lambda   Lcm    }                    {curlyL  J1  x   } { Lambda curlyL x } { Lcm Lambda x  } { J   Lcm  x  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // We rewrite the 6js in a way that better matches what we precompute
          // { J     jc     sc     J2     }
          // {   Jab    lc      L2   J12  } = sum_x (-1)^(S-x) {J1 curlyL   Jab } { curlyL Lambda  lc } { Lcm  Lambda  L2 } { Lcm J12 J  }
          // { J1  curlyL Lambda   Lcm    }                    {jc J        x   } { sc     jc       x } { sc   J2      x  } { J1  x   J2 }
          //
          // 
//          std::cout << "twelvej " << twoJ << " " << j2c << " " << 1 << " " << twoJ2
//                           << "  " << 2*Jab << " " << 2*lc << " " << 2*L2 << " " << twoJ12
//                           << "  " << 2*J1  << " " << 2*curlyL << " " << 2*Lambda << " " << 2*Lcm << std::endl;
//          double tstart = omp_get_wtime();
// The 12j obeys the tetragonal conditions for (J,J1,sc,Lambda) and (jc,curlyL,J2,Lcm).
// For a tetrad, each j should not be greater than the sum of the other 3.
          if (   ( twoJ> (2*J1+1+2*Lambda))      or ( (2*J1)>(twoJ+1+2*Lambda))     or ((2*Lambda)>(twoJ+1+2*J1))
              or ( j2c > (2*curlyL+twoJ2+2*Lcm)) or ( (2*curlyL)>(j2c+twoJ2+2*Lcm)) or (twoJ2>(j2c+2*curlyL+2*Lcm)) or ((2*Lcm)>(j2c+2*curlyL+twoJ2)) ) continue;


//          double twelvej = GetTwelveJ( twoJ, j2c, 1, twoJ2,
          double twelvej = ComputeTwelveJ( twoJ, j2c, 1, twoJ2,
                                            2*Jab, 2*lc, 2*L2, twoJ12,
                                           2*J1, 2*curlyL, 2*Lambda, 2*Lcm );
//          IMSRGProfiler::timer["ComputeTwelveJ"] += omp_get_wtime() - tstart;
//          double twelvej = 0;
//          int twox_min = std::max( { std::abs(2*curlyL-j2c),  std::abs(twoJ-2*J1), std::abs(2*Lambda-1), std::abs(2*Lcm-twoJ2)   });
//          int twox_max = std::min( { (2*curlyL+j2c),  (twoJ+2*J1), (2*Lambda+1),  (2*Lcm+twoJ2)   });
//          //                           sixj1,sixj2    sixj1,sixj4   sixj2,sixj3     sixj3,sixj4            
//          for (int twox = twox_min; twox<=twox_max; twox +=2)
//          {
//            twelvej += (twox+1) * AngMom::phase( (twoJ+j2c+1+twoJ2+2*Jab+2*lc+2*L2+twoJ12+2*J1+2*curlyL+2*Lambda+2*Lcm -twox )/2 ) 
//                                                   * GetSixJ( 2*J1, 2*curlyL,  2*Jab, j2c, twoJ,  twox)
//                                                   * GetSixJ( 2*curlyL , 2*Lambda,  2*lc, 1, j2c,  twox)
//                                                   * GetSixJ( 2*Lcm, 2*Lambda , 2*L2,  1,     twoJ2, twox)
//                                                   * GetSixJ( 2*Lcm, twoJ12, twoJ,   2*J1, twox,  twoJ2); 
////                                                   * AngMom::SixJ_int( 2*Lcm, twoJ12, twoJ,   2*J1, twox,  twoJ2); 
//          }

//          sum_L = AngMom::phase( Lcm + (1-twoJ2)/2 + S1 + curlyL + Jab + lc + Lab + Lambda  )  *  AngMom::SixJ_int( 2*J1,2*curlyL,2*Jab, 2*Lab, 2*S1, 2*L1) * twelvej;
          sum_L = AngMom::phase( Lcm + (1-twoJ2)/2 + S1 + curlyL + Jab + lc + Lab + Lambda  )  *  GetSixJ( 2*Jab,2*J1,2*curlyL, 2*L1, 2*Lab, 2*S1) * twelvej;



//          int L_min = std::max( std::abs(Lambda-L1), std::abs(Lab-lc) );
//          int L_max = std::min( Lambda+L1, Lab+lc );
//
//          for (int L=std::max(std::abs(Lambda-L1),L_min); L<=std::min(L_max,Lambda+L1); L++)
//          {
////            double sixj1 =  AngMom::phase(L+Lambda) * (2*L+1) * AngMom::SixJ( lc, curlyL, Lambda, L1,L,Lab);
//            // { lc curlyL Lambda } Triangles: (lc,curlyL,Lambda), (L1,L,Lambda), (lc,L,Lab), (L1,curlyL,Lab)
//            // { L1   L     Lab   }
//            double sixj1 =  AngMom::phase(L+Lambda) * (2*L+1) * GetSixJ( 2*lc, 2*curlyL, 2*Lambda, 2*L1,2*L,2*Lab);
//            if (std::abs(sixj1)<1e-9) continue;
//
//            double sum_12 = 0;
//            for (int twoS12=twoS12_min; twoS12<=twoS12_max; twoS12+=2)
//            {
//              if ( std::abs(2*L-twoS12)>twoJ or 2*L+twoS12<twoJ ) continue; // triangle condition for ninejL
//              double ninejL =  GetNineJ( 2*lc, 2*L, 2*Lab,
//                                         1, twoS12, 2*S1,
//                                         j2c, twoJ, 2*Jab);
////              double ninejL =  AngMom::NineJ( Lab, lc, L,
////                                              S1,  sc, 0.5*twoS12,
////                                              Jab, jc, 0.5*twoJ);
//
//              if ( std::abs(ninejL)<1e-9) continue;
//              int L12_min = std::max( std::abs(L1-L2), std::max( std::abs(twoS12-twoJ12)/2, std::abs(Lcm-L)) ) ;
//              int L12_max = std::min( L1+L2, std::min( (twoS12+twoJ12)/2, Lcm+L));
//              for (int L12=L12_min; L12<=L12_max; L12++)
//              {
//                if ( std::abs(Lcm-L12)>L or (Lcm+L12)<L ) continue;
//                if ( std::abs(twoS12-twoJ)>2*L or (twoS12+twoJ)<2*L ) continue;
//                double ninej12 = AngMom::phase((twoS12+twoJ)/2)*(2*L12+1)*(twoS12+1)
//                                              * GetNineJ(2*L2,  2*L12,  2*L1,
//                                                         1,     twoS12, 2*S1,
//                                                         twoJ2, twoJ12, 2*J1);
//                if (std::abs(ninej12)<1e-9) continue;
////                double sixj2 = AngMom::SixJ( Lcm,L12,L,0.5*twoS12,0.5*twoJ,0.5*twoJ12);
//                // { Lcm  L12  L   }  Triangles : (Lcm,L12,L),  (S12,J,L),  (Lcm,J,J12),  (S12,L12,J12)
//                // { S12  J    J12 }
//                double sixj2 = GetSixJ( 2*Lcm,2*L12,2*L, twoS12,twoJ,twoJ12);
////                double sixj3 = AngMom::SixJ( Lcm,L2,Lambda,L1,L,L12);
//                double sixj3 = GetSixJ( 2*Lcm,2*L2,2*Lambda,2*L1,2*L,2*L12);
//
//                sum_12 += ninejL * ninej12 * sixj2 * sixj3;
//              } // for L12
//            } // for twoS12
//            sum_L += sixj1 * sum_12;
//          } // for L


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
//    IMSRGProfiler::timer[__func__] += omp_get_wtime() - timestamp;
    return tcoeff * jhats * globalphase;
  }






          // { J     jc     sc     J2     }
          // {   Jab    lc      L2   J12  } = sum_x (-1)^(S-x) {J1 curlyL   Jab } { curlyL Lambda  lc } { Lcm  Lambda  L2 } { Lcm J12 J  }
          // { J1  curlyL Lambda   Lcm    }                    {jc J        x   } { sc     jc       x } { sc   J2      x  } { J1  x   J2 }

//  { a1   a2   a3    a4    }
//  {   b12  b23   b34  b41 }  = sum_x (-1)^(S-x)  { c1  c2  b12 } { c2  c3  b23 } { c4 c3 b34 } { c4 b41 a1 }
//  { c1   c2   c3   c4     }                      { a2  a1  x   } { a3  a2  x   } { a3 a4 x   } { c1 x   a4 }
//

//  { a1   a2   a3    a4    }
//  {   b12  b23   b34  b41 }  = sum_x (-1)^(S-x)  { c1  b12 c2  } { c2  b23 c3  } { c4 b34 c3 } { c4 b41 a1 }
//  { c1   c2   c3   c4     }                      { a2  x   a1  } { a3  x   a2  } { a3 x   a4 } { c1 x   a4 }
double Jacobi3BME::ComputeTwelveJ(int a1, int a2, int a3, int a4, int b12, int b23, int b34, int b41, int c1, int c2, int c3, int c4)
{
  double twelvej = 0;

// no need to check this if we do it in the precompute stage
//  if (   (a1>(a3+c1+c3)) or (a3>(a1+c1+c3)) or (c1>(a1+a3+c3)) or (c3>(a1+a3+c1))
//      or (a2>(a4+c2+c4)) or (a4>(a2+c2+c4)) or (c2>(a2+a4+c4)) or (c4>(a2+a4+c2)) ) return 0.0;

//  if (c4==0) // this check doesn't seem to do all that much
//  {
//    if ( (b41!=a1) or (c3!=b34) ) return 0;
//    if ( std::abs(c1-c2)>b12 or (c1+c2)<b12 or std::abs(a2-b41)>b12 or (a2+b41)<b12 or std::abs(c1-b41)>a4 or (c1+b41)<a4 or std::abs(c2-a2)>a4 or (c2+a2)<a4 ) return 0;
//    twelvej = AngMom::phase( (b12+b23+b34+b41-a2-c2)/2 ) * GetSixJ( c1,c2, b12, a2,b41,a4) * GetSixJ(c2,b34,b23,a3,a2,a4) /  sqrt( (b41+1.)*(b34+1) );
//    return twelvej;
//  }

  int x_min = std::max( { std::abs(c1-a1),  std::abs(c2-a2), std::abs(c3-a3), std::abs(c4-a4)   });
  int x_max = std::min( { a1+c1, a2+c2, a3+c3, a4+c4   });
  int S = a1+a2+a3+a4 + b12+b23+b34+b41 + c1+c2+c3+c4;
//  std::cout << "      x range " << x_min << " " << x_max << std::endl;
  for (int x = x_min; x<=x_max; x +=2)
  {
//    std::cout << "     sixJs  ( " << c1  << " " << b12 << " " << c2<< "  " << a2 << " " << x<< " " << a1  << " ) "
//                     << "     ( " << c2  << " " << b23 << " " << c3<< "  " << a3 << " " << x<< " " << a2  << " ) "
//                     << "     ( " << c4  << " " << b34 << " " << c3<< "  " << a3 << " " << x<< " " << a4  << " ) "
//                     << "     ( " << c4 << " " << b41 << " " << a1 << "  " << c1 << " " << x << " " << a4 << " ) "
//                     << std::endl;
//    double sixj1 = GetSixJ( c1,c2,b12, a2,a1,x);
//    double sixj2 = GetSixJ( c2,c3,b23, a3,a2,x);
//    double sixj3 = GetSixJ( c4,c3,b34, a3,a4,x);
    double sixj1 = GetSixJ( c1,b12,c2, a2,x,a1);
    double sixj2 = GetSixJ( c2,b23,c3, a3,x,a2);
    double sixj3 = GetSixJ( c4,b34,c3, a3,x,a4);
    double sixj4 = GetSixJ( c4,b41,a1, c1,x,a4);
    twelvej += (x+1) * AngMom::phase( (S+x )/2 ) * sixj1 * sixj2 * sixj3 * sixj4;
//    twelvej += (x+1) * AngMom::phase( (S+x )/2 ) 
//                                           * GetSixJ( c1,c2,b12, a2,a1,x)
//                                           * GetSixJ( c2,c3,b23, a3,a2,x)
//                                           * GetSixJ( c4,c3,b34, a3,a4,x)
//                                           * GetSixJ( c4,b41,a1, c1,x,a4); 

//    std::cout << "      evaluate to " << sixj1 << " " << sixj2 << " " << sixj3 << " " << sixj4 << "   =>   " << (x+1) << " " << AngMom::phase( (S+x)/2) << "  "
//              << std::setprecision(8) <<std::fixed << sixj1 * sixj2 * sixj3 * sixj4 << "   =  " << twelvej << std::endl;
  }
  return twelvej;

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
    if (omp_get_num_threads()<0)
//    if (omp_get_num_threads()<2)
    {
      #pragma omp critical
      {
        SixJList[key] = sixj;
      }
    }
    else
    {
      std::cout << "DANGER!!!!!!!  Updating SixJList inside a parellel loop breaks thread safety!" << std::endl;
      std::cout << "  I shouldn't be here in GetSixJ("
                << std::setprecision(1) << std::fixed << j1 << " " << std::setprecision(1) << std::fixed << j2 << " "
                << std::setprecision(1) << std::fixed << j3 << " " << std::setprecision(1) << std::fixed << J1 << " "
                << std::setprecision(1) << std::fixed << J2 << " " << std::setprecision(1) << std::fixed << J3 << "). key = "
                << std::hex << key << "   sixj = " << std::dec << sixj << std::endl;
      IMSRGProfiler::counter["N_CalcSixJ_in_Parallel_loop"] +=1;
//      exit(EXIT_FAILURE);
    }
   }
   return sixj;
}

double Jacobi3BME::GetTwelveJ(int a1, int a2, int a3, int a4,  int b12, int b23, int b34, int b41, int c1, int c2, int c3, int c4)
{
  if (a3==1)
  {
    auto key = MakeUshort11({a1,a2,a4,b12,b23,b34,b41,c1,c2,c3,c4});
    auto iter = TwelveJList.find(key);
    if ( iter != TwelveJList.end() )
    {
      return iter->second;
    }
  }
  double twelvej = ComputeTwelveJ(a1,a2,a3,a4,b12,b23,b34,b41,c1,c2,c3,c4);
  std::cout << "Asked for a TwelveJ and missed: " << a1  << " " << a2  << " " <<  a3 << " " << a4  << " "
                                                  << b12 << " " << b23 << " " << b34 << " " << b41 << " "
                                                  << c1  << " " << c2  << " " << c3  << " " << c4  << "   ->  " << twelvej << std::endl;
  return twelvej ;

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
//  we also have this one
//  { J1  curlyL  Jab }  =>  { Jab  curlyL  J1 }    // 2*Jab < 6emax+1 .. good.   ;   2*curlyL<=6emax+1 ... good.  ;  J1 Triangle.
//  { Lab   S1    L1  }      { L1  S1      Lab }   // 2*L1 <= 2*Nmax+1 ... good.  ; S1 triangle ... good ;   Lab  triangle ... good;      check that Jab + Lab < max(E3max,Nmax) ??
//
// Limit on Lcm : must be <= E3max from energy cons.
// Limit on curlyL : must be <= E2max from energy cons.
// Limit on Lambda : Triangle.
// Limit on L2 : must be <= Nmax
// Limit on L : 3 * lmax.
//
//   GetSixJ( 2*Jab,2*J1,2*curlyL, 2*L1, 2*Lab, 2*S1)
//
//           double twelvej = ComputeTwelveJ( twoJ, j2c, 1, twoJ2,
//                                            2*Jab, 2*lc, 2*L2, twoJ12,
//                                           2*J1, 2*curlyL, 2*Lambda, 2*Lcm );
//
//     twelvej += (x+1) * AngMom::phase( (S+x )/2 ) 
//                                           * GetSixJ( c1,c2,b12, a2,a1,x)
//                                           * GetSixJ( c2,c3,b23, a3,a2,x)
//                                           * GetSixJ( c4,c3,b34, a3,a4,x)
//                                           * GetSixJ( c4,b41,a1, c1,x,a4);
//  The 6J coefficients that we need for the Tcoefficients are:
//   1: { J1  curlyL  Jab }   2: { curlyL Lambda lc }   3: { Lcm Lambda L2 } 4: { Lcm J12 J }  5: { Jab  J1  curlyL }
//      { jc  J       x   }      { 1/2    jc     x  }      { 1/2   J2   x  }    { J1  x  J2 }     { L1   Lab  S1    }
//
// Different ordering
//   1: { J1 Jab curlyL }   2: { curlyL lc Lambda }   3: { Lcm  L2 Lambda } 4: { Lcm  J12  J }  5: { Jab  J1  curlyL }
//      { jc x   J      }      { 1/2    x  jc     }      { 1/2  x    J2   }    { J1   x   J2 }     { L1   Lab  S1    }
//
//  Limits J1 <= Nmax
//         Lcm <= E3max
//         Jab <= min(2emax+1,E3max+1)
//         Lambda is x+- 1/2, and 
//
// { Jab J1 curlyL }
// { L1  Lab  S1   }
//
//
  void Jacobi3BME::PreComputeSixJ()
  {
    double t_start = omp_get_wtime();
    std::cout << "Precalculating SixJ's" << std::endl;
    std::vector<uint64_t> KEYS;
//    std::cout << "j1 goes up to " << (2*E2max)+2 << std::endl;
//    std::cout << "j2 goes up to " <<  2*Nmax+2 << std::endl;
//    std::cout << "J4 goes up to " << 2*Nmax << std::endl;
//    std::cout << "limit on j1+J6 = " << 2*std::max(E3max,Nmax) << std::endl;
//
//  
// { Jab J1 curlyL }
// { L1  Lab  S1   }
    for (int twoJab=0; twoJab<=2*E2max+2; twoJab+=2)
    {
     for (int twoJ1=0; twoJ1<=2*Nmax+2; twoJ1+=2)
     {
      for (int twocurlyL=std::abs(twoJab-twoJ1); twocurlyL<=(twoJab+twoJ1); twocurlyL+=2)
      {
       for (int twoL1=0; twoL1<=2*Nmax; twoL1+=2)
       {
        for (int twoS1=0; twoS1<=2; twoS1+=2)
        {
         for (int twoLab=std::abs(twoJab-twoS1); twoLab<=twoJab+twoS1; twoLab+=2)
         {
          // maybe apply some other constraints here? Tetragons?
           uint64_t key = Jacobi3BME::SixJHash(twoJab,twoJ1,twocurlyL,twoL1,twoLab,twoS1);
           if ( SixJList.count(key) == 0 ) 
           {
             KEYS.push_back(key);
             SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
           }
         }
        }
       }
      }
     }
    }

// 
//        double mosh1 = mosh1phase  * GetMoshinsky1( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab); // we've cached the d=1 brackets
//          double mosh2 = moshphase2 * (2*Lambda+1) * GetMoshinsky2(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda);
//  mosh1 limits curlyL to E2max
//
//   1: { J1 Jab curlyL }   2: { curlyL lc Lambda }   3: { Lcm  L2 Lambda }  
//      { jc x   J      }      { 1/2    x  jc     }      { 1/2  x    J2   }  
    for (int twol1=0; twol1<=std::max({2*Nmax+2, 2*E2max, 2*E3max}); twol1+=2)
    {
     for (int twol2=0; twol2<=std::max({2*E2max+2,2*emax,2*Nmax}); twol2+=2)
     {
      for (int twol3=std::abs(twol1-twol2); twol3<=twol1+twol2; twol3+=2)
      {
       for (int twoj1=1; twoj1<=2*emax+1; twoj1+=2)
       {
        if (twoj1 + twol2 > 2*E3max+3) continue;
        for (int twoj2=std::abs(twoj1-twol3); twoj2<=twoj1+twol3; twoj2+=2)
        {
         for (int twoj3=std::max(std::abs(twol1-twoj2),std::abs(twol2-twoj1)); twoj3<=std::min(twol1+twoj2,twol2+twoj1); twoj3+=2)
         {
          // check some tetragons
           uint64_t key = Jacobi3BME::SixJHash(twol1,twol2,twol3,twoj1,twoj2,twoj3);
           if ( SixJList.count(key) == 0 ) 
           {
             KEYS.push_back(key);
             SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
           }
         }
        }
       }
      }
     }
    }

// third type of 6j { Lcm  J12  J }
//                  { J1   x   J2 }

   for (int twoLcm=0; twoLcm<=2*E3max; twoLcm+=2)
   {
    for (int twoJ12=1; twoJ12<=twoJmax; twoJ12+=2)
    {
     for (int twoJ=std::abs(twoLcm-twoJ12); twoJ<=std::min(twoLcm+twoJ12,2*E3max+3); twoJ+=2)
     {
      for (int twoJ1=0; twoJ1<=2*Nmax+2; twoJ1+=2)
      {
       for (int twox=std::abs(twoJ1-twoJ); twox<=twoJ1+twoJ; twox+=2)
       {
        for (int twoJ2=std::max(std::abs(twoJ12-twoJ1),std::abs(twoLcm-twox)); twoJ2<=std::min({twoJ12+twoJ1,twoLcm+twox,2*Nmax+3-twoJ1}); twoJ2+=2)
        {
         // more checking ?
           uint64_t key = Jacobi3BME::SixJHash(twoLcm,twoJ12,twoJ,twoJ1,twox,twoJ2);
           if ( SixJList.count(key) == 0 ) 
           {
             KEYS.push_back(key);
             SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
           }
        }
       }
      }
     }
    }
   }

/////////////////////////////////////////
//    // first, we do the all-integer ones
////    for (int j1=0; j1<=(6*emax+1); j1+=2)  // 2 * Lcm can go up to 6emax. Wait, no. It can only to up to 2*E3max
//    for (int j1=0; j1<=std::max(2*E3max,4*emax+2); j1+=2)  // 2 * Lcm can go up to 6emax. Wait, no. It can only to up to 2*E3max
//    {
////     for (int j2=0; j2<=(6*emax+1); j2+=2)
////     for (int j2=0; j2<=std::max({4*emax,2*Nmax+3,2*E3max+2}); j2+=2)
////                                   Jab       L2,J1      J12      J1
//     for (int j2=0; j2<=std::max({4*emax+2, 2*E3max+2-j1, twoJmax, 2*Nmax+2}); j2+=2)
//     {
//      for (int j3=std::abs(j1-j2); j3<=(j1+j2); j3+=2)
//      {
////       for (int J4=0; J4<=3*(2*emax+1); J4+=2)
//       for (int J4=0; J4<=(2*Nmax+0); J4+=2)
//       {
//        for (int J5=std::abs(j3-J4); J5<=(j3+J4); J5+=2)
//        {
//         for (int J6=std::max(std::abs(j1-J5),std::abs(j2-J4)); J6<=std::min(j1+J5,j2+J4); J6+=2)
//         {
//           if (J6>2) continue;
////           if ( (j1 + J6)>2*std::max(E3max,Nmax)) continue;   // why this???
//           uint64_t key = Jacobi3BME::SixJHash(j1,j2,j3,J4,J5,J6);
//           if ( SixJList.count(key) == 0 ) 
//           {
//             KEYS.push_back(key);
//             SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
//           }
//         } // for J6
//        } // for J5
//       } // for J4
//
//       // Now we do the half-integer bottom row
//       //       double sixj2 = GetSixJ( 2*Lcm,2*L12,2*L, twoS12,twoJ,twoJ12); twoJ <= 6emax+3, twoJ12<=Jmax, L could be up to 6emax+6
//       // { j1  j2 j3  }
//       // { S12 J  J12 }
//
//        // replaced by
//          // { J     jc     sc     J2     }
//          // {   Jab    lc      L2   J12  } = sum_x (-1)^(S-x) {curlyL J1   Jab } { Lambda curlyL lc } { Lcm  Lambda  L2 } { J2  J1  J12 }
//          // { J1  curlyL Lambda   Lcm    }                    {J      jc   x   } { jc     sc      x } { sc   J2      x  } { J   Lcm  x  }
////       if ((j1+j2)>2*E3max) continue;
////       for (int S12=1; S12<=3; S12+=2)
//// 
//       for (int S12=1; S12<=std::max(twoJmax,2*emax+1); S12+=2)
//       {
//        for (int J=std::abs(S12-j3); J<=(S12+j3); J+=2)
//        {
//         for (int J12=std::max(std::abs(j1-J),std::abs(j2-S12)); J12<=std::min((j1+J),(j2+S12)); J12+=2)
//         {
//           uint64_t key = Jacobi3BME::SixJHash(j1,j2,j3,S12,J,J12);
//           if ( SixJList.count(key) == 0 ) 
//           {
//             KEYS.push_back(key);
//             SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
//           }
//         } // for J12
//        } // for J
//       } // for S12
//      } // for j3
//     } // for j2
//    } // for j1
//
//    // one last form that doesn't fit with the others
//    // { Lcm  J12  J2 }
//    // { J1   x    J  }
//   // { J2  J1  J12 }  =>  { Lcm  J12  J  }
//   // { J   Lcm  x  }      { J1   x    J2 }
////          int twox_min = std::max( { std::abs(2*curlyL-j2c),  std::abs(twoJ-2*J1), std::abs(2*Lambda-1), std::abs(2*Lcm-twoJ2)   });
////          int twox_max = std::min( { (2*curlyL+j2c),  (twoJ+2*J1), (2*Lambda+1),  (2*Lcm+twoJ2)   });
////   x = Lambda += 1/2;  Triangle( Lambda,Lcm,L2), Triangle(Lambda,curlyL,lc),  Triangle(L2,J2,1),   so Lambda + Lcm <= L2 <= J2+1    so x <= Lambda+1 <= Lcm + L2 + 1 <= Lcm + J2 + 2
////  Triangle( x,sc,Lambda) 
//
//
//     for (int twoJ12=1; twoJ12<=twoJmax; twoJ12+=2)
//     {
////      for (int twoJ=std::abs(twoLcm-twoJ12); twoJ<=(twoLcm+twoJ12); twoJ+=2)
//      for (int twoJ=1; twoJ<=(2*E3max+3); twoJ+=2)
//      {
////      for (int twoLcm=0; twoLcm<=(6*emax+1); twoLcm+=2)  // 2 * Lcm can go up to 6emax
////      for (int twoLcm=std::abs(twoJ-twoJ12); twoLcm<=std::min(twoJ+twoJ12,6*emax); twoLcm+=2)  // 2 * Lcm can go up to 6emax
//      for (int twoLcm=std::abs(twoJ-twoJ12); twoLcm<=std::min({twoJ+twoJ12,2*E3max,6*emax,2*E3max+3-twoJ12}); twoLcm+=2)  // 2 * Lcm can go up to 6emax
//      {
//        for (int twoJ1=0; twoJ1<=2*(Nmax+1); twoJ1+=2)
//        {
//         for (int twox=std::abs(twoJ-twoJ1); twox<=(twoJ+twoJ1); twox+=2)
//         {
////          for (int twoJ2=std::max({ std::abs(twoLcm-twox), std::abs(twoJ1-twoJ12), twox-twoLcm-2  } ); twoJ2<=std::min( { (twoLcm+twox), (twoJ1+twoJ12), 2*Nmax+3-twoJ1 } ); twoJ2+=2)
//          for (int twoJ2=std::max({ std::abs(twoLcm-twox), std::abs(twoJ1-twoJ12)  } ); twoJ2<=std::min( { (twoLcm+twox), (twoJ1+twoJ12), 2*Nmax+3-twoJ1 } ); twoJ2+=2)
//          {
////           if ( twoJ2 < twox - twoLcm - 2 ) continue;
////           if ( (twoJ1+twoJ2) > 2*Nmax+3) continue;
//           if ( (twoJ1+twoJ2+twoLcm) > 2*(E3max+3) ) continue;
//           uint64_t key = Jacobi3BME::SixJHash(twoLcm,twoJ12,twoJ,twoJ1,twox,twoJ2);
//           if ( SixJList.count(key) == 0 ) 
//           {
//             KEYS.push_back(key);
//             SixJList[key] = 0.; // Make sure eveything's in there to avoid a rehash in the parallel loop
//           }
//           
//          }
//         }
//        }
//      }
//     }
//    }
    
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


//  Our 12j is of the form  
//  { J     jc     1/2    J2      }
//  {   Jab      lc     L     J12 }
//  { J1    curlyL Lambda  Lcm    }
//
// Triangles:  (Jab,jc,J)  (lc,1/2,jc)  ( L,1/2,J2)   (J1,J2,J12)    (J1,curlyL,Jab)   ( curlyL,Lambda,lc)  (Lambda,L,Lcm)  (Lcm,J12,J)
//
// Tetragonal conditions:  (J,J1,Lambda,1/2) , (jc,curlyL,J2,Lcm)
//  the tetragonal conditions (see Varshalovich ch 10) on (a,b,c,d) require  a+b+c+d = integer,   a<=b+c+d,  b<=a+c+d, c<=a+b+d, d<=a+b+c 
//
// Energy constraints: lc<=emax, Jab<=min(2*emax,E2max)+1, J2<=Nmax+1/2, J1<=Nmax+1, (J1+J2)<=Nmax+3/2,  J12<=Nmax+3/2,   Lcm<=E3max,   J<=E3max+3/2,  Lcm+J12 <= E3max+3/2
//
// 
//        double mosh1 = mosh1phase  * GetMoshinsky1( curlyN,curlyL,  N1,L1,      na,la,     nb,lb,     Lab); // we've cached the d=1 brackets
//          double mosh2 = moshphase2 * (2*Lambda+1) * GetMoshinsky2(    Ncm,Lcm,     N2,L2,  curlyN,curlyL, nc,lc,  Lambda);
//
  void Jacobi3BME::PreComputeTwelveJ()
  {
    double tstart = omp_get_wtime();
    std::vector< std::array<unsigned short,11>> keylist;
    for (int twoJ12=1; twoJ12<=twoJmax; twoJ12+=2)
    {
     for (int twoJ1=0; twoJ1<=2*Nmax+2; twoJ1+=2)
     {
      for (int twoJ2=std::abs(twoJ12-twoJ1); twoJ2<=std::min({2*Nmax+1, 2*Nmax+3-twoJ1, twoJ12+twoJ1}); twoJ2+=2)
      {
       for (int twolc=0; twolc<=2*emax; twolc+=2 )
       {
        for (int twojc=std::abs(twolc-1); twojc<=twolc+1; twojc+=2)
        {
         for (int twoLcm=0; twoLcm<=std::min(2*E3max, 2*E3max+3-twoJ12); twoLcm+=2)
         {
          if ( (twoJ1 + twoJ2 + twoLcm) > 2*E3max+3) continue;
          for (int twoJ=std::abs(twoJ12-twoLcm); twoJ<=std::min(twoJ12+twoLcm, 2*E3max+3); twoJ+=2)
          {
           for (int twoJab=std::abs(twoJ-twojc); twoJab<=std::min(twoJ+twojc, 2*std::min(E2max,2*emax)+2); twoJab+=2)
           {
            if ( (twoJab + twolc)>2*E3max+2 ) continue;
            for (int twoL=std::abs(twoJ2-1); twoL<=twoJ2+1; twoL+=2)
            {
             for (int twoLambda=std::abs(twoLcm-twoL); twoLambda<=twoLcm+twoL; twoLambda+=2)
             {
              if (twoJ>(twoLambda+twoJ1+1)) continue; // tetragonal conditions
              if (twoJ1>(twoLambda+twoJ+1)) continue;
              if (twoLambda>(twoJ+twoJ1+1)) continue;
              for (int twocurlyL=std::max( std::abs(twoJab-twoJ1),std::abs(twoLambda-twolc)); twocurlyL<=std::min( twoJab+twoJ1, twoLambda+twolc); twocurlyL+=2)
              {
               if (twocurlyL>2*E2max) continue; // energy conservation in moshinsky bracket
               if (twocurlyL+twoJ1>2*E2max+2) continue; // energy conservation in moshinsky bracket
               if (twojc>(twocurlyL+twoJ2+twoLcm)) continue;
               if (twocurlyL>(twojc+twoJ2+twoLcm)) continue;
               if (twoJ2>(twojc+twocurlyL+twoLcm)) continue;
               if (twoLcm>(twojc+twocurlyL+twoJ2)) continue;
               auto key = MakeUshort11({twoJ,twojc,twoJ2, twoJab,twolc,twoL,twoJ12, twoJ1,twocurlyL,twoLambda,twoLcm});
               keylist.push_back( key );
               TwelveJList[key] = 0.0;
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
    #pragma omp parallel for schedule(dynamic,1)
    for (int i=0;i<keylist.size(); i++)
    {
      auto& key = keylist[i];
      int a1=key[0];
      int a2=key[1];
      // we skip a3 because it's sc=1/2, so we don't need to store that information.
      int a4=key[2];
      int b12=key[3];
      int b23=key[4];
      int b34=key[5];
      int b41=key[6];
      int c1=key[7];
      int c2=key[8];
      int c3=key[9];
      int c4=key[10];
//      std::cout << std::endl << "Precomputing TwelveJ: " << a1 << " " << a2 << " " << 1 << " " << a4 << " "
//                                            << b12 << " " << b23 << " " << b34 << " " << b41 << " "
//                                            << c1  << " " << c2  << " " << c3  << " " << c4  << std::endl;
      TwelveJList[key] = ComputeTwelveJ(a1,a2,1,a4,b12,b23,b34,b41,c1,c2,c3,c4);
//      std::cout << "got " << std::fixed << std::setprecision(6) << TwelveJList[key] << std::endl;
    }
    
    std::cout << "Done with loops in " << __func__ << "  size of keylist is " << keylist.size() << std::endl;
    IMSRGProfiler::timer[__func__] += omp_get_wtime() - tstart;
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






