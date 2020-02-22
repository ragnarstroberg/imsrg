
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
double Jacobi3BME::GetMatElAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p)
{
   size_t dimket = GetDimensionAS(twoT,twoJ,p,Nket); 
   size_t start_loc = GetStartLocAS(twoT,twoJ,Nbra,Nket);
   return meAS.at( start_loc + ibra*dimket + iket );
}

// Access an antisymmetrized matrix element
void Jacobi3BME::SetMatElAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p, double matel)
{
   size_t dimket = GetDimensionAS(twoT,twoJ,p,Nket); 
   size_t start_loc = GetStartLocAS(twoT,twoJ,Nbra,Nket);
   meAS.at( start_loc + ibra*dimket + iket ) = matel;
}


double Jacobi3BME::GetMatElNAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p)
{
//   auto hash = HashTJNN(twoT,twoJ,Nbra,Nket);
   size_t dimket = GetDimensionNAS(twoT,twoJ,p,Nket); 
   size_t start_loc = GetStartLocNAS(twoT,twoJ,Nbra,Nket);
//   std::cout << "start_loc = " << start_loc << "  elements: " << std::endl;
//   size_t dimbra = GetDimensionNAS(twoT,twoJ,p,Nbra);
//   for (int i=0;i<dimbra*dimket; i++)
//   {
//     std::cout << meNAS.at( start_loc + i ) << " ";
//     if (i%10==0) std::cout << std::endl;
//   }
//   std::cout << std::endl;
   return meNAS.at( start_loc + ibra*dimket + iket );
}

void Jacobi3BME::SetMatElNAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p, double matel)
{
//   auto hash = HashTJNN(twoT,twoJ,Nbra,Nket);
   size_t dimket = GetDimensionNAS(twoT,twoJ,p,Nket); 
   size_t start_loc = GetStartLocNAS(twoT,twoJ,Nbra,Nket);
   meNAS.at( start_loc + ibra*dimket + iket ) = matel;
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
           // and rows=ket, cols=bra
            arma::mat ASmat( &meAS[startAS], dim_ketAS, dim_braAS, /*copy_aux_mem*/false);
            arma::mat NASmat = cfp_ket * ASmat * cfp_bra.t();
//            std::cout << "ASmat: " << std::endl << ASmat << std::endl << std::endl;
//            std::cout << "NASmat = " << std::endl << NASmat << std::endl << std::endl;

            for (size_t iNAS=0; iNAS<dim_ketNAS*dim_braNAS; iNAS++)
            {
              meNAS[startNAS+iNAS] = NASmat(iNAS); // single index assumes a flat layout with column-major ordering
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
//  int LabE3max = 16;

  if (bra.p==bra.q and bra.p==bra.r and bra.op->j2<3 and twoT>1) return 0;
  if (ket.p==ket.q and ket.p==ket.r and ket.op->j2<3 and twoT>1) return 0;

  int Ebra = 2*(bra.op->n+bra.oq->n+bra.oR->n ) + (bra.op->l+bra.oq->l+bra.oR->l);
  int Eket = 2*(ket.op->n+ket.oq->n+ket.oR->n ) + (ket.op->l+ket.oq->l+ket.oR->l);

//  int parity_bra = Ebra%2;
//  int parity_ket = Eket%2;

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
//        std::cout << " ^^^TJP E12 = " << twoT << " " << twoJ12 << " " << parity12_bra << " " << E12_bra << "," << E12_ket << "  Ncm,Lcm = " << Ncm << " " << Lcm << std::endl;
//        std::cout << "dimensions: " << NASdim_bra << " " << NASdim_ket << std::endl;
//           size_t dimbra = jacobi_basis.GetDimensionNAS(T,J,parity12_bra,E12_bra); 
//           std::cout << "dimbra = " <<  dimbra << std::endl;

//  compute the Tcoefficients that we'll need here
        std::vector<double> Tcoeff_bra(NASdim_bra);
        std::vector<double> Tcoeff_ket(NASdim_ket);
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
          Tcoeff_bra.at(ibraNAS) = Tcoeff_wrapper( bra, Jab, twoJ, jac1_bra, jac2_bra, twoJ12, Ncm, Lcm) ;
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
          Tcoeff_ket.at(iketNAS) = Tcoeff_wrapper( ket, Jde, twoJ, jac1_ket, jac2_ket, twoJ12, Ncm, Lcm) ;
        }


//        std::cout << "    Tcoef_bra: ";
//        for (auto tc : Tcoeff_bra ) std::cout << tc << " ";
//        std::cout << std::endl;
//        std::cout << "    Tcoef_ket: ";
//        for (auto tc : Tcoeff_ket ) std::cout << tc << " ";
//        std::cout << std::endl;
        for (int ibraNAS=0; ibraNAS<NASdim_bra; ibraNAS++)
        {
          if (std::abs(Tcoeff_bra[ibraNAS])<1e-9) continue;
          jacobi1_state jac1_bra,jac1_ket;
          jacobi2_state jac2_bra,jac2_ket;

//          std::cout << "about to get jac1_bra and jac2_bra" << std::endl;
          GetJacobiStates( twoT, twoJ12, parity12_bra, E12_bra, ibraNAS, jac1_bra, jac2_bra);
//          std::cout << "checking isospin bra: " << jac1_bra.t << " " << Tab << std::endl;
          if ( jac1_bra.t != Tab ) continue;
//          double Tcoeff_bra = Tcoeff_wrapper( bra, Jab, twoJ, jac1_bra, jac2_bra, twoJ12, Ncm, Lcm);
//          if (std::abs(Tcoeff_bra)<1e-9) continue;

          for (int iketNAS=0; iketNAS<NASdim_ket; iketNAS++)
          {
           if (std::abs(Tcoeff_ket[iketNAS])<1e-9) continue;
//           std::cout << "  ibraNAS,iketNAS = " << ibraNAS << " " << iketNAS << std::endl;
           GetJacobiStates( twoT, twoJ12, parity12_ket, E12_ket, iketNAS, jac1_ket, jac2_ket);
//           std::cout << "     checking isospin ket: " << jac1_ket.t << " " << Tde << std::endl;
           if ( jac1_ket.t != Tde ) continue;
//           double Tcoeff_ket = Tcoeff_wrapper( ket, Jde, twoJ, jac1_ket, jac2_ket, twoJ12, Ncm, Lcm);

//           std::cout << "    ** finding the NAS matrix element. T = " << T << "   J = " << J12 << "  p = " << parity12_bra << std::endl;

//           double meNAS = ElementNAS( ibraNAS, iketNAS, E12_bra, E12_ket, twoT, twoJ12, parity12_bra );
           double meNAS = GetMatElNAS( ibraNAS, iketNAS, E12_bra, E12_ket, twoT, twoJ12, parity12_bra );

//           me_lab +=  Tcoeff_bra * meNAS * Tcoeff_ket  ;
           me_lab +=  Tcoeff_bra[ibraNAS] * meNAS * Tcoeff_ket[iketNAS]  ;
//           std::cout << "Ecm,Lcm,twoJ12 = " << Ecm << " " << Lcm << " " << twoJ12 << "  Tab,Tde,twoT = " << Tab << " " << Tde << " " << twoT << std::endl;
//           std::cout << "jac1_bra " << jac1_bra.n << " " << jac1_bra.l << " " << jac1_bra.s << " " << jac1_bra.j << " " << jac1_bra.t << "  "
//                     << "jac2_bra " << jac2_bra.n << " " << jac2_bra.l << " " << jac2_bra.j2 <<  "  "
//                     << "jac1_ket " << jac1_ket.n << " " << jac1_ket.l << " " << jac1_ket.s << " " << jac1_ket.j << " " << jac1_ket.t << "  "
//                     << "jac2_bra " << jac2_bra.n << " " << jac2_bra.l << " " << jac2_bra.j2  << std::endl;
//           std::cout << " " << ibraNAS << " , " << iketNAS << ":  (( " << Tcoeff_bra[ibraNAS] << "  " << meNAS << "  " << Tcoeff_ket[iketNAS] << " ))  "
//                     << Tcoeff_bra[ibraNAS] * meNAS * Tcoeff_ket[iketNAS] <<  "  => sum_ME = " << me_lab << std::endl;
//           std::cout << "    < N1=" << jac1_bra.n << " L1=" << jac1_bra.l << " S1=" << jac1_bra.s << " J1=" << jac1_bra.j << " T1=" <<jac1_bra.t << ",  N2=" << jac2_bra.n << " L2=" << jac2_bra.l << " T=" << twoT << " J12=" << twoJ12 << " | ..." << std::endl;

          }// for ibraNAS
        } // for iketNAS
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
std::string Jacobi3BME::TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ, uint64_t jac1, uint64_t jac2, uint64_t twoJ12, uint64_t Lcm )
{
  std::ostringstream oss;
  oss << na << " " << nb << " " << nc << " " << Jab << " " << twoJ << " " << jac1 << " " << jac2 << " " << twoJ12 << " " << Lcm;
  return oss.str();
}

//void Jacobi3BME::TcoeffUnHash(uint64_t key, uint64_t& na, uint64_t& nb, uint64_t& nc, uint64_t& Jab, uint64_t& twoJ, uint64_t& jac1, uint64_t& jac2, uint64_t& twoJ12, uint64_t& Lcm )
void Jacobi3BME::TcoeffUnHash(std::string& key, int& na, int& nb, int& nc, int& Jab, int& twoJ, int& jac1, int& jac2, int& twoJ12, int& Lcm )
{
  std::istringstream iss(key);
  iss >> na >> nb >> nc >> Jab >> twoJ >> jac1 >> jac2 >> twoJ12 >> Lcm;
}


/*

//void Jacobi3BME::GetV3mon_all( std::vector<uint64_t>& keys, std::vector<double>& v3mon, ModelSpace& modelspace )
void Jacobi3BME::GetV3mon_all( HartreeFock& hf )
{
  double t_start = omp_get_wtime();
  // The keys are already computed elsewhere and are passed in as input
  std::cout << "Begin " << __func__ << std::endl;
  hf.Vmon3.resize( hf.Vmon3_keys.size(), 0.);
  struct ljChannel{ int l; int j2;};
  std::vector<ljChannel> ljchannels;
  for ( auto& obc : hf.modelspace->OneBodyChannels )
  {
    ljChannel ljchan(obc.first[0], obc.first[1] );
    if ( std::find( ljchan, ljchannels.begin(), ljchannels.end() ) == ljchannels.end() ) ljchannels.push_back( ljchan) ;
  }
//  v3mon.resize( keys.size(), 0.);
  int n_mon = 0;
  for (auto& obc_a : hf.modelspace->OneBodyChannels )
  {
    int la = obc_a.first[0];
    int j2a = obc_a.first[1];
    int tz2a = obc_a.first[2];
    for (auto& obc_b : hf.modelspace->OneBodyChannels )
    {
      int lb = obc_b.first[0];
      int j2b = obc_b.first[1];
      int tz2b = obc_b.first[2];
      int Jab_min = std::abs(j2a-j2b)/2;
      int Jab_max = (j2a+j2b)/2;
      int Tzab = (tz2a+tz2b)/2;
      int Tab_min = std::abs(Tzab);
       // possibilities  (tz2a,tz2b,Tab)=>CG :  (1,1,1)=>1 ,  (1,-1,1)=>sqrt(2), (-1,1,1)=>sqrt(2),  (1,-1,0)=>sqrt(2) , (-1,1,0)=>-sqrt(2)
      std::array<double,2> isospin2_Clebsch = {tz2a*sqrt(0.5), 0.5*(tz2a+tz2b) + std::abs(tz2a-tz2b)/2*sqrt(0.5) };
      for (auto& obc_c : hf.modelspace->OneBodyChannels )
      {
        int lc = obc_c.first[0];
        int j2c = obc_c.first[1];
        int tz2c = obc_c.first[2];
        std::cout << "OneBody Channels: (" << la << " " << j2a << " " << tz2a << ") , (" << lb << " " << j2b << " " << tz2b << ") , (" << lc << " " << j2c << " " << tz2c << ")"
                  << "   monopoles computed = " << n_mon << std::endl;
        int twoTz = tz2a + tz2b + tz2c;
        int twoT_min = std::abs( twoTz ); // this will either be 1 or 3
        // index by  2*Tab + T/2 ->             { (0,1), (0,3), (1,1), (1,3) }
        std::array<double,4> isospin3_Clebsch = {   1.0,   0.0,  AngMom::CG(1,Tzab,0.5,0.5*tz2c, 0.5, 0.5*twoTz),  AngMom::CG(1,Tzab,0.5,0.5*tz2c, 1.5,0.5*twoTz) };
        // possibilities (Tab,Tzab, 2tc, 2tzc| 2T)=>CG :  (0,0, 1,1|1)=>1 ,  (0,0, 1,-1|1)=>1,     (1,1, 1,1| 3,3)=>1,  (1,-1,1,-1|3,-3)=>1,   (1,0, 1,1|3,1)
        // compute all the relevant T coefficients.  The coefficients have the following indices:
        //  na, nb, nc, Jab, J, jacobi1, jacobi2, J12, Lcm   (Ncm can be inferred from energy conservation?)
//        std::unordered_map<uint64_t,double> T3bList;  // we will put things into a hash table.
        std::unordered_map<std::string,double> T3bList;  // we will put things into a hash table.

        double t_internal = omp_get_wtime();
        // First, generate all the keys that we'll need, but don't compute the T coefficients
//        std::cout << "  generate all the keys..." << std::endl;
        for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
        {
          int twoJ_min = std::abs(2*Jab-j2c);
          int twoJ_max = 2*Jab+j2c;
          for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
          {
            for (auto a : obc_a.second) // run over the list of orbit indices in channel obc_a
            {
             Orbit& oa = hf.modelspace->GetOrbit( a );
             for (auto b : obc_b.second)
             {
              if (a==b and Jab%2>0 ) continue;
              Orbit& ob = hf.modelspace->GetOrbit( b );
              for (auto c : obc_c.second)
              {
                if (a==b and a==c and j2a<3 ) continue;
                Orbit& oc = hf.modelspace->GetOrbit( c );
                int Eabc = 2*(oa.n+ob.n+oc.n) + oa.l+ob.l+oc.l;
                if (Eabc>Nmax) continue;

                for (int jac1=0; jac1<jacobi_1.size(); jac1++)
                {
                 int N1 = jacobi_1[jac1].n;
                 int L1 = jacobi_1[jac1].l;
                 int J1 = jacobi_1[jac1].j;
                 for (int jac2=0; jac2<jacobi_2.size(); jac2++)
                 {
                   int N2 = jacobi_2[jac2].n;
                   int L2 = jacobi_2[jac2].l;
                   int twoJ2 = jacobi_2[jac2].j2;
                   if ( (2*(N1+N2)+L1+L2) > Eabc ) continue;
                   int twoJ12_min = std::abs(2*J1-twoJ2);
                   int twoJ12_max = 2*J1+twoJ2;
                   for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
                   {
                    for (int Lcm=(L1+L2+Eabc)%2; Lcm<=(Eabc-2*(N1+N2)-L1-L2); Lcm++)
                    {
                      if ( (std::abs(2*Lcm - twoJ12)>twoJ) or (2*Lcm+twoJ12)<twoJ ) continue; // check triangle condition for J12,Lcm J
                      auto hash_key = TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,jac1,jac2,twoJ12,Lcm);
                      T3bList[hash_key] = 0;
                    }// for Lcm
                   } // for twoJ12
                 } // for jac2
               } // for jac1
              } // for c
             } // for b
            } // for a
          } // for twoJ
        } // for Jab

        IMSRGProfiler::timer[std::string(__func__)+"_GenerateKeys"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
        std::cout << "There are " << T3bList.size() << " Tcoefficients to compute in this channel" << std::endl;

//        std::cout << "done. Now compute the T coefficients" << std::endl;

        // Next, we compute all the T coefficients. We'll probably want to do this in parallel.
        // In the current construction, all threads loop through the entire set of elements, but
        // they only stop to compute when it's their turn.
        #pragma omp parallel
        {
          size_t cnt = 0;
          int ithread = omp_get_thread_num();
          int nthreads = omp_get_num_threads();
          for(auto element = T3bList.begin(); element !=T3bList.end(); ++element, cnt++)
          {
            if(cnt%nthreads != ithread) continue; // Check if this thread should compute this element
            auto hash_key = element->first;
            int na,nb,nc,Jab,twoJ,jac1,jac2,twoJ12,Lcm;
            TcoeffUnHash(hash_key, na,nb,nc,Jab,twoJ,jac1,jac2,twoJ12,Lcm);
            auto& jacobi1 = jacobi_1[jac1];
            auto& jacobi2 = jacobi_2[jac2];
            int Ncm = (2*(na+nb+nc)+la+lb+lc - 2*(jacobi1.n+jacobi2.n) - jacobi1.l - jacobi2.l - Lcm)/2;
            if (Ncm<0) continue;
            double tcoef = AngMom::Tcoeff( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jacobi1.n, jacobi1.l, jacobi1.s, jacobi1.j, jacobi2.n, jacobi2.l, jacobi2.j2, twoJ12, Ncm, Lcm);
            element->second = tcoef;
          }
        }

        IMSRGProfiler::timer[std::string(__func__)+"_ComputeToeff"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();

        std::vector<size_t> imonlist_all;
        std::vector<size_t> imonlist;
        std::vector<size_t> imonlist_swap_ab;
        std::vector<size_t> imonlist_swap_de;
        std::vector<size_t> imonlist_swap_abde;
//        for (size_t imon=0; imon<v3mon.size(); imon++)
        for (size_t imon=0; imon<hf.Vmon3_keys.size(); imon++)
        {
          // Check if the given key corresponds to the set of one-body channels we're currently working on.
          auto key = hf.Vmon3_keys[imon];
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          if ( std::find( obc_a.second.begin(), obc_a.second.end(), a ) == obc_a.second.end() ) continue;
          if ( std::find( obc_b.second.begin(), obc_b.second.end(), b ) == obc_b.second.end() ) continue;
          if ( std::find( obc_c.second.begin(), obc_c.second.end(), c ) == obc_c.second.end() ) continue;
//          imonlist.push_back( imon );
          imonlist_all.push_back( imon );
          if (a<=b and d<=e) imonlist.push_back(imon);
//          else if (a>=b and d<=e) imonlist_swap_ab.push_back(imon);
//          else if (a<=b and d>=e) imonlist_swap_de.push_back(imon);
//          else if (a>=b and d>=e) imonlist_swap_abde.push_back(imon);
//          std::cout << " imon = " << imon << "   =>  " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
        }
        // This is fast, so we can still afford to do this twice without parallelization
        // Also, we're now working with a greatly reduced subset of the monopole terms

        imonlist_swap_ab.resize( imonlist.size(),-1 );
        imonlist_swap_de.resize( imonlist.size(),-1 );
        imonlist_swap_abde.resize( imonlist.size(),-1 );
//        std::cout << "Begin imonlistloop. imonlist_all:" << std::endl;
//        for ( auto imon : imonlist_all ) std::cout << imon << " ";
//        std::cout << std::endl;
//        for (auto imon : imonlist )
        for (int ind=0; ind<imonlist.size(); ind++ )
        {
          auto imon = imonlist[ind];
          auto key = hf.Vmon3_keys[imon];
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          auto key_ab = hf.Vmon3Hash( b,a,c,d,e,f);
          auto key_de = hf.Vmon3Hash( a,b,c,e,d,f);
          auto key_abde = hf.Vmon3Hash( b,a,c,e,d,f);
          for ( auto iall : imonlist_all )
          {
            if (hf.Vmon3_keys[iall] == key_ab) imonlist_swap_ab[ind] = iall;
            if (hf.Vmon3_keys[iall] == key_de) imonlist_swap_de[ind] = iall;
            if (hf.Vmon3_keys[iall] == key_abde) imonlist_swap_abde[ind] = iall;
          }
//          std::cout << " ind = " << ind << " : " << imon << "  " << imonlist_swap_ab[ind] << "  " << imonlist_swap_de[ind] << "  " << imonlist_swap_abde[ind] << std::endl;
//          std::cout << "        " << a << " " << b << " " << c << "  " << d << " " << e << " " << f << std::endl;
//          std::cout << "        " << key << "  " << key_ab << " " << key_de << " " << key_abde << std::endl;
//          int aa,bb,cc,dd,ee,ff;
//          hf.Vmon3UnHash( key_ab ,aa,bb,cc,dd,ee,ff);
//          std::cout << "     " << aa << " " << bb << " " << cc << " " << dd << " " <<ee << " " << ff << "   ";
//          hf.Vmon3UnHash( key_de ,aa,bb,cc,dd,ee,ff);
//          std::cout << "     " << aa << " " << bb << " " << cc << " " << dd << " " <<ee << " " << ff << "   ";
//          hf.Vmon3UnHash( key_abde ,aa,bb,cc,dd,ee,ff);
//          std::cout << "     " << aa << " " << bb << " " << cc << " " << dd << " " <<ee << " " << ff << "   " << std::endl;
        }
        std::cout << "done. size of swapab = " << imonlist_swap_ab.size() << std::endl;

        IMSRGProfiler::timer[std::string(__func__)+"_FindMonKeys"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
        bool verbose = false;
//        if (lb==3 and lc==3) verbose = true;

        if (verbose) std::cout << "done.  Now the loop over imon... Size of imonlist is " << imonlist.size() << std::endl;
        // And at last, we compute all the monopole terms for this combination of one-body channels
        // This can also be done in parallel
//        for (size_t imon=0; imon<v3mon.size(); imon++)
//        #pragma omp parallel for schedule(dynamic,1)
        #pragma omp parallel for schedule(dynamic,1)
        for (size_t ilist=0; ilist<imonlist.size(); ilist++)
        {
          size_t imon = imonlist[ilist];
          size_t imon_ab = imonlist_swap_ab[ilist];
          size_t imon_de = imonlist_swap_de[ilist];
          size_t imon_abde = imonlist_swap_abde[ilist];
          auto key = hf.Vmon3_keys[imon];
//          if (verbose) std::cout << " imon,imon_ab,imon_de,imon_abde = " << imon << " " << imon_ab << " " << imon_de << " " << imon_abde << std::endl;
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          // V3mon(a,b,c,d,e,f) = 1/(2jc+1) * sum_{Jab,J} (2J+1) <abc Jab,J| V |def Jab,J>
          //                    = 1/(2jc+1) * sum_{Jab,J} (2J+1) sum_{jac1,jac2,J12,jac1',jac2',J12',Ncm,Lcm}  Tcoef(abc,Jab,J;jac1,jac2,J12,Ncm,Lcm) * Tcoef(def,Jab,J;jac1',jac2',J12',Ncm,Lcm)
          //                                                                                 sum_{T12,T12',T}  * CG(ta,tb,T12)*CG(tc,T12,T) * CG(td,te,T12')*CG(df,T12',T)   <- Here T means isospin
          //                                                                                                   * < jac1,jac2,J12,T12 | V | jac1',jac2',J12',T12' >
          if (a>b) continue;
          if (d>e) continue;
          Orbit& oa = hf.modelspace->GetOrbit(a);
          Orbit& ob = hf.modelspace->GetOrbit(b);
          Orbit& oc = hf.modelspace->GetOrbit(c);
          Orbit& od = hf.modelspace->GetOrbit(d);
          Orbit& oe = hf.modelspace->GetOrbit(e);
          Orbit& of = hf.modelspace->GetOrbit(f);
          int Eabc = 2*(oa.n+ob.n+oc.n) + oa.l+ob.l+oc.l;
          int Edef = 2*(od.n+oe.n+of.n) + od.l+oe.l+of.l;
          if (Eabc>Nmax or Edef>Nmax) continue;
          int parity = Eabc%2;
          double v_monopole = 0;
          double v_monopole_swap_ab = 0;
          double v_monopole_swap_de = 0;
          int Jab_step = (a==b or d==e)? 2 : 1;
          for (int Jab=Jab_min; Jab<=Jab_max; Jab+=Jab_step)
          {
           if ( (a==b and a==c and j2a<3) or (d==e and d==f and od.j2<3) ) continue; // can't fit 3 identical particles in a j=1/2 orbit
           int twoJ_min = std::abs(2*Jab-j2c);
           int twoJ_max = 2*Jab+j2c;
           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
           {
             for (int Tab=Tab_min; Tab<=1; Tab++)
             {
              int ab_swap_phase = - AngMom::phase( (j2a+j2b)/2 - Jab + 1 - Tab); // don't forget that minus sign out in front...
              double isoClebsch_ab = isospin2_Clebsch[Tab];
              for (int Tde=Tab_min; Tde<=1; Tde++)
              {
                double isoClebsch_de = isospin2_Clebsch[Tde];
                for (int twoT=twoT_min; twoT<=std::min(2*Tab+1,2*Tde+1); twoT+=2)
                {
                  int de_swap_phase = - AngMom::phase( (j2a+j2b)/2 - Jab + 1 - Tde); // don't forget that minus sign out in front...
                  double isoClebsch_c = isospin3_Clebsch[2*Tab + twoT/2];
                  double isoClebsch_f = isospin3_Clebsch[2*Tde + twoT/2];
                  if (std::abs(isoClebsch_c)<1e-6 or std::abs(isoClebsch_f)<1e-6) continue;
                  double v_sumJT = 0;
                  for (int Ecm=0; Ecm<=std::min(Eabc,Edef); Ecm++)
                  {
                   int E12abc = Eabc-Ecm;
                   int E12def = Edef-Ecm;
                   if (E12abc > Nmax  or E12def>Nmax) continue;
                   for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
                   {
                     int twoJ12_min = std::abs(twoJ-2*Lcm);
                     int twoJ12_max = std::min( twoJ+2*Lcm, twoJmax);
                     for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
                     {
                       if (verbose)   std::cout << "twoJ,Lcm,twoJ12,twoJmax = " << twoJ << " " << Lcm << " " << twoJ12 << " " << twoJmax << std::endl;
                       if (verbose)   std::cout << "Ecm,E12abc,E12def " << Ecm << " " << E12abc << " " << E12def  << std::endl;
                       auto hashTJN_abc = HashTJN(twoT,twoJ12,E12abc);
                       auto hashTJN_def = HashTJN(twoT,twoJ12,E12def);
                       if (verbose) std::cout << "hash_abc = " << hashTJN_abc << " hash_def = " << hashTJN_def << "  size of dimNAS = " << dimensionNAS.size() << std::endl;
                       size_t dimNAS_abc = GetDimensionNAS( twoT, twoJ12, parity, E12abc ); 
                       size_t dimNAS_def = GetDimensionNAS( twoT, twoJ12, parity, E12def ); 

                       size_t startloc = GetStartLocNAS(twoT, twoJ12, E12abc, E12def) ;
//                       size_t otherstartloc = GetStartLocNAS(twoT, twoJ12, E12def, E12abc) ;

//                       std::cout << "   dimensions: " << dimNAS_abc << " , " << dimNAS_def << std::endl;
                       // As a reminder, armadillo stores matrices   [ M11  M12  M13 ]
                       // in column-major order, as in               | M21  M22  M23 |
                       //                                            [ M31  M32  M33 ]
                       //
                       //  [ Tdef_1  Tdef_2 ... ]  *  [ <def|V|abc>  ... ]  *  [ Tabc_1 ]
                       //                             [   ...   ...  ... ]     [ Tabc_2 ]
                       //                             [   ...   ...  ... ]     [  ...   ]

                       arma::vec Tabc(dimNAS_abc, arma::fill::zeros);
                       arma::rowvec Tdef(dimNAS_def, arma::fill::zeros);
                       arma::mat matelNAS( &meNAS[startloc], dimNAS_def, dimNAS_abc, false );
//                       std::cout << "Getting matrix beginning at " << startloc << "  ( should it be " << otherstartloc << " ? ) " <<  std::endl;
                       if (verbose) std::cout << "Getting matrix beginning at " << startloc << "   size of vector = " << meNAS.size() <<  std::endl;
                       for (size_t iNAS_abc = 0; iNAS_abc<dimNAS_abc; iNAS_abc++)
                       {
                         auto& index_1_2_abc = NAS_jacobi_states.at(hashTJN_abc).at(iNAS_abc);
                         if ( jacobi_1.at(index_1_2_abc[0]).t != Tab ) continue;
                         Tabc[iNAS_abc] = T3bList[ TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,index_1_2_abc[0],index_1_2_abc[1],twoJ12,Lcm )];
                       }
                       for (size_t iNAS_def = 0; iNAS_def<dimNAS_def; iNAS_def++)
                       {
                         auto& index_1_2_def = NAS_jacobi_states.at(hashTJN_def).at(iNAS_def);
                         if ( jacobi_1.at(index_1_2_def[0]).t != Tde ) continue;
                         Tdef[iNAS_def] = T3bList[ TcoeffHash(od.n,oe.n,of.n,Jab,twoJ,index_1_2_def[0],index_1_2_def[1],twoJ12,Lcm )];
                       }

                       if (verbose) std::cout << "Doing the mat mult" << std::endl;
                       arma::mat result = Tdef * matelNAS * Tabc;
//                       arma::mat backresult = Tabc.t() * matelNAS * Tdef.t();
                       v_sumJT += result[0];
//                       std::cout << "matrix multiply version: " << std::endl << Tdef << std::endl << std::endl << matelNAS << std::endl << std::endl << Tabc << std::endl << std::endl << std::endl;
//                       std::cout << "matelNAS * Tabc :" << std::endl << ( matelNAS * Tabc ) << std::endl;
//                       std::cout << "result = " << result[0] << "    or " << backresult[0] << std::endl;
                  

                       
//                       for (size_t iNAS_abc = 0; iNAS_abc<dimNAS_abc; iNAS_abc++)
//                       {
//                         auto& index_1_2_abc = NAS_jacobi_states.at(hashTJN_abc).at(iNAS_abc);
//                         auto indx_jac1_abc = index_1_2_abc[0];
//                         auto indx_jac2_abc = index_1_2_abc[1];
//                         jacobi1_state& jac1_abc = jacobi_1.at(indx_jac1_abc);
//                         jacobi2_state& jac2_abc = jacobi_2.at(indx_jac2_abc);
////                         GetJacobiStates( twoT, twoJ, parity, E12abc, iNAS_abc, jac1_abc, jac2_abc);
//                         if (jac1_abc.t != Tab) continue;
//
//                         for (size_t iNAS_def = 0; iNAS_def<dimNAS_def; iNAS_def++)
//                         {
//                           auto& index_1_2_def = NAS_jacobi_states.at(hashTJN_def).at(iNAS_def);
//                           auto indx_jac1_def = index_1_2_def[0];
//                           auto indx_jac2_def = index_1_2_def[1];
//                           jacobi1_state& jac1_def = jacobi_1.at(indx_jac1_def);
//                           jacobi2_state& jac2_def = jacobi_2.at(indx_jac2_def);
//
////                           GetJacobiStates( twoT, twoJ, parity, E12def, iNAS_def, jac1_def, jac2_def);
//                           if (jac1_def.t != Tde) continue;
//                           if ( std::abs(twoJ12-2*Lcm)>twoJ or (twoJ12+2*Lcm)<twoJ ) continue;
//                           auto THash_abc = TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,indx_jac1_abc,indx_jac2_abc,twoJ12,Lcm );
//                           auto THash_def = TcoeffHash(od.n,oe.n,of.n,Jab,twoJ,indx_jac1_def,indx_jac2_def,twoJ12,Lcm );
//                           double Tcoef_abc = T3bList[THash_abc];
//                           double Tcoef_def = T3bList[THash_def];
//                           if ( std::abs(Tcoef_abc*Tcoef_def)>1e-6)
//                           {
////                             std::cout << "        Tcoefs = " << Tcoef_abc << " " << Tcoef_def << std::endl;
////                             std::cout << "       accesing " << iNAS_abc << " " << iNAS_def << " " << E12abc << " " << E12def << " " << twoT << " " << twoJ << " " << parity << std::endl;
////                             std::cout << "       Tab,Tde,twoT = " << Tab << " " << Tde << " " << twoT << "   isospin clebsch's " << isoClebsch_ab << " " << isoClebsch_c << "   " << isoClebsch_de << " " << isoClebsch_f << std::endl;
//                             double v_NAS = GetMatElNAS( iNAS_abc, iNAS_def,  E12abc,  E12def, twoT, twoJ12,  parity);
////                             v_sumJT +=  Tcoef_abc * Tcoef_def * v_NAS;
////                             std::cout << "  " << iNAS_abc << " " << iNAS_def << "  " << Tcoef_abc << " " << Tcoef_def << "  " << v_NAS << "   " <<  Tcoef_abc * Tcoef_def * v_NAS << "  sum = " << v_sumJT << std::endl;
////                             std::cout << "         v_NAS = " << v_NAS << "  v_sumJT = " << v_sumJT << std::endl;
//                           }
//                         } // for iNAS_def
//                       } // for iNAS_abc
                     } // for twoJ12
                   } // for Lcm
                  } // for Ecm
                  // v_sumJT is equal to <abc Jab Tab JT | V | def Jde Tde JT>
                  Ket3 bra = Ket3(oa,ob,oc);
                  Ket3 ket = Ket3(od,oe,of);
                  if (verbose) std::cout << "Assigning" << std::endl;
//                  double vcheck_JT =  GetLabMatEl( bra, ket, Jab, Jab, twoJ, Tab, Tde, twoT);
                  double vterm = 6 * (twoJ+1) * isoClebsch_ab * isoClebsch_c * isoClebsch_de * isoClebsch_f * v_sumJT;
                  v_monopole += vterm;
                  v_monopole_swap_ab += ab_swap_phase * vterm;
                  v_monopole_swap_de += de_swap_phase * vterm;
//                  std::cout << "    vmonopole: += " << twoJ+1 << " * " << isoClebsch_ab << " * " << isoClebsch_c << " * " << isoClebsch_de << " * " << isoClebsch_f << " * " << 6*v_sumJT
//                            << "  ( " << vcheck_JT << " )   " << "  => " << v_monopole << std::endl;
                } // for twoT
              } // for Tde
             } // for Tab
           } // for twoJ
          } // for Jab
          v_monopole /= j2c+1.;
          v_monopole_swap_ab /= j2c+1.;
          v_monopole_swap_de /= j2c+1.;
         // There are likely some symmetries we can exploit here to avoid redundant calculations
          hf.Vmon3[imon] = v_monopole; // put that bad boy into the output vector
          if (imon_abde !=-1 and imon_abde != imon)
            hf.Vmon3[imon_abde] = v_monopole; // put that bad boy into the output vector
          if (imon_ab !=-1)
            hf.Vmon3[imon_ab] = v_monopole_swap_ab; // put that bad boy into the output vector
          if (imon_de !=-1)
            hf.Vmon3[imon_de] = v_monopole_swap_de; // put that bad boy into the output vector
          if (verbose)
          {
            std::cout << "abcdef = " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
            std::cout << "SETTING " << imon << " -> " << v_monopole << std::endl << std::endl;
            if ( imon_abde !=-1 and imon_abde != imon )
            std::cout << "Also SETTING " << imon_abde << " -> " << v_monopole << std::endl;
            if ( imon_ab !=-1 and imon_ab != imon )
            std::cout << "Also SETTING " << imon_ab << " -> " << v_monopole_swap_ab << std::endl;
            if ( imon_de !=-1 and imon_de != imon )
            std::cout << "Also SETTING " << imon_de << " -> " << v_monopole_swap_de << std::endl;

          }
        } // for imon
        n_mon += imonlist.size();
//        std::cout << "done " << std::endl;

        IMSRGProfiler::timer[std::string(__func__)+"_iMonLoop"] += omp_get_wtime() - t_internal;
      } // for obc_c
    } // for obc_b
  } // for obc_a

  // all done!

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
}



*/

//void Jacobi3BME::GetV3mon_all( std::vector<uint64_t>& keys, std::vector<double>& v3mon, ModelSpace& modelspace )
void Jacobi3BME::GetV3mon_all( HartreeFock& hf )
{
  double t_start = omp_get_wtime();
  // The keys are already computed elsewhere and are passed in as input
  std::cout << "Begin " << __func__ << std::endl;
  hf.Vmon3.resize( hf.Vmon3_keys.size(), 0.);
  struct ljChannel{ int l; int j2; bool operator==(const ljChannel& rhs){return (rhs.l==l and rhs.j2==j2); };};
  std::vector<ljChannel> ljchannels;
  for ( auto& obc : hf.modelspace->OneBodyChannels )
  {
    ljChannel ljchan = {obc.first[0], obc.first[1] };
    if ( std::find( ljchannels.begin(), ljchannels.end(), ljchan ) == ljchannels.end() ) ljchannels.push_back( ljchan) ;
  }

  int n_mon = 0;
  size_t num_lj = ljchannels.size();
  for ( size_t ilj_a=0; ilj_a<num_lj; ilj_a++ )
  {
    int la = ljchannels[ilj_a].l;
    int j2a = ljchannels[ilj_a].j2;
    std::set<int> na_list; // a std::set is a sorted unique list of items
    for (auto& orb : hf.modelspace->Orbits ) if (orb.l==la and orb.j2==j2a) na_list.insert(orb.n); 
//    for ( size_t ilj_b=0; ilj_b<num_lj; ilj_b++ ) // TODO: change the lower limit from 0 to ilj_a
    for ( size_t ilj_b=ilj_a; ilj_b<num_lj; ilj_b++ ) // TODO: change the lower limit from 0 to ilj_a
    {
      int lb = ljchannels[ilj_b].l;
      int j2b = ljchannels[ilj_b].j2;
      std::set<int> nb_list; // a std::set is a sorted unique list of items
      for (auto& orb : hf.modelspace->Orbits ) if (orb.l==lb and orb.j2==j2b) nb_list.insert(orb.n); 
      int Jab_min = std::abs(j2a-j2b)/2;
      int Jab_max = (j2a+j2b)/2;
      for ( size_t ilj_c=0; ilj_c<num_lj; ilj_c++ ) // TODO: possibly also change lower limit from 0 to ilj_b, although this becomes more complicated, due to recoupling
      {
        int lc = ljchannels[ilj_c].l;
        int j2c = ljchannels[ilj_c].j2;
        std::set<int> nc_list; // a std::set is a sorted unique list of items
        for (auto& orb : hf.modelspace->Orbits ) if (orb.l==lc and orb.j2==j2c) nc_list.insert(orb.n); 

        // This is where we can compute the T coefficients
        double t_internal = omp_get_wtime();
        std::unordered_map<std::string,double> T3bList;  // we will put things into a hash table.

        for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
        {
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
              if (a>b) continue;
              for (auto nc : nc_list )
              {
               int Eabc = 2*(na+nb+nc) + la+lb+lc;

               for (int Ecm=0; Ecm<=Eabc; Ecm++)
               {
                int E12 = Eabc-Ecm;
                if (E12 > Nmax) continue;
                int parity = E12%2;
                for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
                {
                  int twoJ12_min = std::abs(twoJ-2*Lcm);
                  int twoJ12_max = std::min( twoJ+2*Lcm, twoJmax);
                  for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
                  {
                    for (int twoT=1; twoT<=3; twoT+=2)
                    {
                      size_t dim_NAS = GetDimensionNAS(twoT,twoJ12, parity, E12 ); 
                      size_t hash = HashTJN(twoT,twoJ12,E12);
                      auto& NAS_jacobi_statelist = NAS_jacobi_states.at(hash);
                      for ( size_t iNAS=0; iNAS<dim_NAS; iNAS++)
                      {
                        size_t jac1_index = NAS_jacobi_statelist[iNAS][0];
                        size_t jac2_index = NAS_jacobi_statelist[iNAS][1];

                        auto hash_key = TcoeffHash(na,nb,nc,Jab,twoJ,jac1_index,jac2_index,twoJ12,Lcm);
                        T3bList[hash_key] = 0;
                      }
                    } // for twoT
                  } // for twoJ12
                 }// for Lcm
               }// for Ecm
              } // for nc
             } // for nb
            } // for na
          } // for twoJ
        } // for Jab

        IMSRGProfiler::timer[std::string(__func__)+"_GenerateKeys"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
        std::cout << "There are " << T3bList.size() << " Tcoefficients to compute in this channel" << std::endl;

//        std::cout << "done. Now compute the T coefficients" << std::endl;

        // Next, we compute all the T coefficients. We'll probably want to do this in parallel.
        // In the current construction, all threads loop through the entire set of elements, but
        // they only stop to compute when it's their turn.
        // TODO: A lot of these end up being zero. Figure out what's going wrong.
        int tzero=0;
        #pragma omp parallel  reduction(+ : tzero )
        {
          int cnt = 0;
          int ithread = omp_get_thread_num();
          int nthreads = omp_get_num_threads();
          for(auto element = T3bList.begin(); element !=T3bList.end(); ++element, cnt++)
          {
            if(cnt%nthreads != ithread) continue; // Check if this thread should compute this element
            auto hash_key = element->first;
            int na,nb,nc,Jab,twoJ,jac1,jac2,twoJ12,Lcm;
            TcoeffUnHash(hash_key, na,nb,nc,Jab,twoJ,jac1,jac2,twoJ12,Lcm);
            auto& jacobi1 = jacobi_1[jac1];
            auto& jacobi2 = jacobi_2[jac2];
            int Ncm = (2*(na+nb+nc)+la+lb+lc - 2*(jacobi1.n+jacobi2.n) - jacobi1.l - jacobi2.l - Lcm)/2;
            if (Ncm<0) continue;
            double tcoef = AngMom::Tcoeff( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jacobi1.n, jacobi1.l, jacobi1.s, jacobi1.j, jacobi2.n, jacobi2.l, jacobi2.j2, twoJ12, Ncm, Lcm);
            if (std::abs(tcoef)<1e-9) tzero++;
//            double tcoef = AngMom::Tcoeff_bruteforce( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jacobi1.n, jacobi1.l, jacobi1.s, jacobi1.j, jacobi2.n, jacobi2.l, jacobi2.j2, twoJ12, Ncm, Lcm);
            element->second = tcoef;
          }
        }

        std::cout << "computed a zero T coefficient " << tzero << " times" << std::endl;
        IMSRGProfiler::timer[std::string(__func__)+"_ComputeToeff"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();



//        std::vector<size_t> imonlist_all;
        std::vector<size_t> imon_dir;
        std::vector<size_t> imon_exch;
        std::vector<size_t> imonlist;
//        std::vector<size_t> imonlist_swap_ab;
//        std::vector<size_t> imonlist_swap_de;
//        std::vector<size_t> imonlist_swap_abde;
        for (size_t imon=0; imon<hf.Vmon3_keys.size(); imon++)
        {
          // Check if the given key corresponds to the set of one-body channels we're currently working on.
          auto key = hf.Vmon3_keys[imon];
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          Orbit& oa = hf.modelspace->GetOrbit(a);
          Orbit& ob = hf.modelspace->GetOrbit(b);
          Orbit& oc = hf.modelspace->GetOrbit(c);
//          Orbit& od = hf.modelspace->GetOrbit(d);
//          Orbit& oe = hf.modelspace->GetOrbit(e);
//          Orbit& of = hf.modelspace->GetOrbit(f);
          if ( oc.l != lc or oc.j2 != j2c ) continue;

          if ( oa.l==la and oa.j2==j2a and ob.l==lb and ob.j2==j2b)
          {
            imon_dir.push_back(imon);
            if (a<d or (a==d and b<e) or (a==d and b==e and c<=f) ) // Choose a specific ordering which we will compute. Then get the others by symmetry.
            {
              if ( a<=b and d<=e and (oa.tz2+ob.tz2+oc.tz2)<0 )
                imonlist.push_back(imon);
            }
          }
          if ( oa.l==lb and oa.j2==j2b and ob.l==la and ob.j2==j2a)
          {
            imon_exch.push_back(imon);
          }
        }

//        std::cout << "sizes:  imon_dir " << imon_dir.size() << "   imonlist " << imonlist.size() << "   imon_exch " << imon_exch.size() << std::endl;

        std::vector<size_t> imonlist_defabc( imonlist.size(), -1);
        std::vector<size_t> imonlist_edfbac( imonlist.size(), -1);
        std::vector<size_t> imonlist_bacedf( imonlist.size(), -1);
        std::vector<size_t> imonlist_isoflip( imonlist.size(), -1);
        std::vector<size_t> imonlist_defabc_isoflip( imonlist.size(), -1);
        std::vector<size_t> imonlist_edfbac_isoflip( imonlist.size(), -1);
        std::vector<size_t> imonlist_bacedf_isoflip( imonlist.size(), -1);
        for ( size_t index=0; index<imonlist.size(); index++)
        {
          size_t imon = imonlist[index];
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
            if ( hf.Vmon3_keys[imon1] == key_isoflip) imonlist_isoflip[index] = imon1;
            if ( hf.Vmon3_keys[imon1] == key_defabc) imonlist_defabc[index] = imon1;
            if ( hf.Vmon3_keys[imon1] == key_defabc_isoflip) imonlist_defabc_isoflip[index] = imon1;
          }
          for (auto imon2 : imon_exch )
          {
            if ( hf.Vmon3_keys[imon2] == key_bacedf) imonlist_bacedf[index] = imon2;
            if ( hf.Vmon3_keys[imon2] == key_edfbac) imonlist_edfbac[index] = imon2;
            if ( hf.Vmon3_keys[imon2] == key_bacedf_isoflip) imonlist_bacedf_isoflip[index] = imon2;
            if ( hf.Vmon3_keys[imon2] == key_edfbac_isoflip) imonlist_edfbac_isoflip[index] = imon2;
          }
          
        }



        IMSRGProfiler::timer[std::string(__func__)+"_FindMonKeys"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
        bool verbose = false;
        if (verbose) std::cout << "done.  Now the loop over imon... Size of imonlist is " << imonlist.size() << std::endl;
        std::cout << "done.  Now the loop over imon... Size of imonlist is " << imonlist.size() << std::endl;

        int tcoeff_counter = 0;

        #pragma omp parallel for schedule(dynamic,1) reduction(+ : tcoeff_counter)
        for (size_t ilist=0; ilist<imonlist.size(); ilist++)
        {
          size_t imon = imonlist[ilist];
          size_t imon_bacedf = imonlist_bacedf[ilist];
          size_t imon_edfbac = imonlist_edfbac[ilist];
          size_t imon_defabc = imonlist_defabc[ilist];
          size_t imon_isoflip = imonlist_isoflip[ilist];
          size_t imon_bacedf_isoflip = imonlist_bacedf_isoflip[ilist];
          size_t imon_edfbac_isoflip = imonlist_edfbac_isoflip[ilist];
          size_t imon_defabc_isoflip = imonlist_defabc_isoflip[ilist];
//          if (imon==8) verbose = true;
//          else verbose = false;
          auto key = hf.Vmon3_keys[imon];
//          if (verbose) std::cout << " imon,imon_ab,imon_de,imon_abde = " << imon << " " << imon_ab << " " << imon_de << " " << imon_abde << std::endl;
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          // V3mon(a,b,c,d,e,f) = 1/(2jc+1) * sum_{Jab,J} (2J+1) <abc Jab,J| V |def Jab,J>
          //                    = 1/(2jc+1) * sum_{Jab,J} (2J+1) sum_{jac1,jac2,J12,jac1',jac2',J12',Ncm,Lcm}  Tcoef(abc,Jab,J;jac1,jac2,J12,Ncm,Lcm) * Tcoef(def,Jab,J;jac1',jac2',J12',Ncm,Lcm)
          //                                                                                 sum_{T12,T12',T}  * CG(ta,tb,T12)*CG(tc,T12,T) * CG(td,te,T12')*CG(df,T12',T)   <- Here T means isospin
          //                                                                                                   * < jac1,jac2,J12,T12 | V | jac1',jac2',J12',T12' >
//          if (a>d) continue; // we can use hermiticity to get the swapped order
//          if (a>b) continue;
//          if (d>e) continue;
          Orbit& oa = hf.modelspace->GetOrbit(a);
          Orbit& ob = hf.modelspace->GetOrbit(b);
          Orbit& oc = hf.modelspace->GetOrbit(c);
          Orbit& od = hf.modelspace->GetOrbit(d);
          Orbit& oe = hf.modelspace->GetOrbit(e);
          Orbit& of = hf.modelspace->GetOrbit(f);
          int Eabc = 2*(oa.n+ob.n+oc.n) + oa.l+ob.l+oc.l;
          int Edef = 2*(od.n+oe.n+of.n) + od.l+oe.l+of.l;
          if (Eabc>Nmax or Edef>Nmax) continue;
          int parity = Eabc%2;

          int Tzab = (oa.tz2+ob.tz2)/2;
          int Tab_min = std::abs(Tzab);
          int twoTz = oa.tz2 + ob.tz2 + oc.tz2;
          if (twoTz>0) continue;
          int twoT_min = std::abs( twoTz ); // this will either be 1 or 3
          std::array<double,2> isospin2_Clebsch = {oa.tz2*sqrt(0.5), Tzab + std::abs(oa.tz2-ob.tz2)/2*sqrt(0.5) };
          std::array<double,4> isospin3_Clebsch = {   1.0,   0.0,  AngMom::CG(1,Tzab,0.5,0.5*oc.tz2, 0.5, 0.5*twoTz),  AngMom::CG(1,Tzab,0.5,0.5*oc.tz2, 1.5,0.5*twoTz) };
          double v_monopole = 0;
//          double v_monopole_swap_abde = 0;
//          double v_monopole_swap_de = 0;
          int Jab_step = (a==b or d==e)? 2 : 1;
          for (int Jab=Jab_min; Jab<=Jab_max; Jab+=Jab_step)
          {
           if ( (a==b and a==c and j2a<3) or (d==e and d==f and od.j2<3) ) continue; // can't fit 3 identical particles in a j=1/2 orbit
           int twoJ_min = std::abs(2*Jab-j2c);
           int twoJ_max = 2*Jab+j2c;
           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
           {
             for (int Tab=Tab_min; Tab<=1; Tab++)
             {
//              int ab_swap_phase = - AngMom::phase( (j2a+j2b)/2 - Jab + 1 - Tab); // don't forget that minus sign out in front...
              double isoClebsch_ab = isospin2_Clebsch[Tab];
              for (int Tde=Tab_min; Tde<=1; Tde++)
              {
                double isoClebsch_de = isospin2_Clebsch[Tde];
                for (int twoT=twoT_min; twoT<=std::min(2*Tab+1,2*Tde+1); twoT+=2)
                {
//                  int de_swap_phase = - AngMom::phase( (j2a+j2b)/2 - Jab + 1 - Tde); // don't forget that minus sign out in front...
                  double isoClebsch_c = isospin3_Clebsch[2*Tab + twoT/2];
                  double isoClebsch_f = isospin3_Clebsch[2*Tde + twoT/2];
                  if (std::abs(isoClebsch_c)<1e-6 or std::abs(isoClebsch_f)<1e-6) continue;
                  double v_sumJT = 0;
                  for (int Ecm=0; Ecm<=std::min(Eabc,Edef); Ecm++)
                  {
                   int E12abc = Eabc-Ecm;
                   int E12def = Edef-Ecm;
                   if (E12abc > Nmax  or E12def>Nmax) continue;
                   for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
                   {
                     int twoJ12_min = std::abs(twoJ-2*Lcm);
                     int twoJ12_max = std::min( twoJ+2*Lcm, twoJmax);
                     for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
                     {
                       if (verbose)   std::cout << "twoJ,Lcm,twoJ12,twoJmax = " << twoJ << " " << Lcm << " " << twoJ12 << " " << twoJmax << std::endl;
                       if (verbose)   std::cout << "twoT,Tab,Tde = " << twoT << " " << Tab << " " << Tde << std::endl;
                       if (verbose)   std::cout << "Ecm,E12abc,E12def " << Ecm << " " << E12abc << " " << E12def  << std::endl;
                       auto hashTJN_abc = HashTJN(twoT,twoJ12,E12abc);
                       auto hashTJN_def = HashTJN(twoT,twoJ12,E12def);
//                       if (verbose) std::cout << "hash_abc = " << hashTJN_abc << " hash_def = " << hashTJN_def << "  size of dimNAS = " << dimensionNAS.size() << std::endl;
                       size_t dimNAS_abc = GetDimensionNAS( twoT, twoJ12, parity, E12abc ); 
                       size_t dimNAS_def = GetDimensionNAS( twoT, twoJ12, parity, E12def ); 

                       size_t startloc = GetStartLocNAS(twoT, twoJ12, E12abc, E12def) ;
//                       size_t otherstartloc = GetStartLocNAS(twoT, twoJ12, E12def, E12abc) ;

//                       std::cout << "   dimensions: " << dimNAS_abc << " , " << dimNAS_def << std::endl;
                       // As a reminder, armadillo stores matrices   [ M11  M12  M13 ]
                       // in column-major order, as in               | M21  M22  M23 |
                       //                                            [ M31  M32  M33 ]
                       //
                       //  [ Tdef_1  Tdef_2 ... ]  *  [ <def|V|abc>  ... ]  *  [ Tabc_1 ]
                       //                             [   ...   ...  ... ]     [ Tabc_2 ]
                       //                             [   ...   ...  ... ]     [  ...   ]

                       arma::vec Tabc(dimNAS_abc, arma::fill::zeros);
                       arma::rowvec Tdef(dimNAS_def, arma::fill::zeros);
                       arma::mat matelNAS( &meNAS[startloc], dimNAS_def, dimNAS_abc, false );
//                       std::cout << "Getting matrix beginning at " << startloc << "  ( should it be " << otherstartloc << " ? ) " <<  std::endl;
                       if (verbose) std::cout << "Getting matrix beginning at " << startloc << "   size of vector = " << meNAS.size() <<  std::endl;
                       for (size_t iNAS_abc = 0; iNAS_abc<dimNAS_abc; iNAS_abc++)
                       {
                         auto& index_1_2_abc = NAS_jacobi_states.at(hashTJN_abc).at(iNAS_abc);
                         if ( jacobi_1.at(index_1_2_abc[0]).t != Tab ) continue;
                         tcoeff_counter++;
                         Tabc[iNAS_abc] = T3bList[ TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,index_1_2_abc[0],index_1_2_abc[1],twoJ12,Lcm )];
                       }
                       for (size_t iNAS_def = 0; iNAS_def<dimNAS_def; iNAS_def++)
                       {
                         auto& index_1_2_def = NAS_jacobi_states.at(hashTJN_def).at(iNAS_def);
                         if ( jacobi_1.at(index_1_2_def[0]).t != Tde ) continue;
                         tcoeff_counter++;
                         Tdef[iNAS_def] = T3bList[ TcoeffHash(od.n,oe.n,of.n,Jab,twoJ,index_1_2_def[0],index_1_2_def[1],twoJ12,Lcm )];
                       }

                       if (verbose) std::cout << "Doing the mat mult" << std::endl;
                       arma::mat result = Tdef * matelNAS * Tabc;
                       v_sumJT += result[0];
                  

                       
                     } // for twoJ12
                   } // for Lcm
                  } // for Ecm
                  // v_sumJT is equal to <abc Jab Tab JT | V | def Jde Tde JT>  (and Jde=Jab)
//                  Ket3 bra = Ket3(oa,ob,oc);
//                  Ket3 ket = Ket3(od,oe,of);
//                  Ket3 bra_flip = Ket3(ob,oa,oc);
//                  Ket3 ket_flip = Ket3(oe,od,of);
                  if (verbose) std::cout << "Assigning" << std::endl;
//                  double vcheck_JT =  GetLabMatEl( bra, ket, Jab, Jab, twoJ, Tab, Tde, twoT);
//                  double vcheck_JT_flip =  GetLabMatEl( bra_flip, ket_flip, Jab, Jab, twoJ, Tab, Tde, twoT);
                  double vterm = 6 * (twoJ+1) * isoClebsch_ab * isoClebsch_c * isoClebsch_de * isoClebsch_f * v_sumJT;
                  v_monopole += vterm;
//                  v_monopole_swap_abde += ab_swap_phase * de_swap_phase * vterm;
//                  v_monopole_swap_abde += 1 * vterm;
//                  v_monopole_swap_de += de_swap_phase * vterm;
//                  if (verbose) std::cout << "    vmonopole: += " << twoJ+1 << " * " << isoClebsch_ab << " * " << isoClebsch_c << " * " << isoClebsch_de << " * " << isoClebsch_f << " * " << 6*v_sumJT
//                                         << "   swap phase = " << ab_swap_phase << " * " << de_swap_phase  << std::endl;
//                  if (verbose) std::cout << "vterm = " <<  vterm << "  v_sumJT = " << v_sumJT*6 << "   vswap = " << v_sumJT*6 * ab_swap_phase * de_swap_phase  << "  checks : " << vcheck_JT << "   " << vcheck_JT_flip << std::endl;
//                  std::cout << "    vmonopole: += " << twoJ+1 << " * " << isoClebsch_ab << " * " << isoClebsch_c << " * " << isoClebsch_de << " * " << isoClebsch_f << " * " << 6*v_sumJT
//                            << "  ( " << vcheck_JT << " )   " << "  => " << v_monopole << std::endl;
                } // for twoT
              } // for Tde
             } // for Tab
           } // for twoJ
          } // for Jab
          v_monopole /= j2c+1.;
//          v_monopole_swap_abde /= j2c+1.;
//          v_monopole_swap_de /= j2c+1.;
         // There are likely some symmetries we can exploit here to avoid redundant calculations
          hf.Vmon3[imon] = v_monopole; // put that bad boy into the output vector

//          int isymmm=0;
          for ( auto& imon_sym : { imon_defabc, imon_bacedf, imon_edfbac, imon_isoflip, imon_defabc_isoflip, imon_bacedf_isoflip, imon_edfbac_isoflip} )
          {
//            if ( hf.Vmon3_keys[imon_sym] == 1127000493261825L) std::cout << " in symmetry,  based on imon = " << imon << "   imon_sym = " << imon_sym <<  "    isymmm= " << isymmm << "   vmonopole = " << v_monopole << std::endl;
            if ( imon_sym != size_t(-1) )  hf.Vmon3[imon_sym] = v_monopole;
//            isymmm++;
          }

//          if (imon_defabc !=-1 )
//            hf.Vmon3[imon_defabc] = v_monopole; // put that bad boy into the output vector
//          if (imon_bacedf !=-1 )
//            hf.Vmon3[imon_bacedf] = v_monopole; // put that bad boy into the output vector
////          if (imon_edfbac !=-1 and imon_edfbac != imon)
//          if (imon_edfbac !=-1 )
//            hf.Vmon3[imon_edfbac] = v_monopole; // put that bad boy into the output vector
//
//          if (imon_defabc_isoflip !=-1 )
//            hf.Vmon3[imon_defabc] = v_monopole; // put that bad boy into the output vector
//          if (imon_bacedf_isoflip !=-1 )
//            hf.Vmon3[imon_bacedf] = v_monopole; // put that bad boy into the output vector
//          if (imon_edfbac_isoflip !=-1 )
//            hf.Vmon3[imon_edfbac_isoflip] = v_monopole; // put that bad boy into the output vector


          if (verbose)
          {
            std::cout << "abcdef = " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
            std::cout << "SETTING " << imon << " -> " << v_monopole << std::endl << std::endl;
            if ( imon_defabc !=size_t(-1) and imon_defabc != imon )
            std::cout << "1Also SETTING " << imon_defabc << " -> " << v_monopole << std::endl;
            if ( imon_bacedf !=size_t(-1) and imon_bacedf != imon )
            std::cout << "2Also SETTING " << imon_bacedf << " -> " << v_monopole << std::endl;
            if ( imon_edfbac !=size_t(-1) and imon_edfbac != imon )
            std::cout << "3Also SETTING " << imon_edfbac << " -> " << v_monopole << std::endl;

          }
        } // for imon
        n_mon += imonlist.size();
//        std::cout << "done " << std::endl;
        std::cout << "Looked up a Tcoefficient " << tcoeff_counter << "  times " << std::endl;
        if (tcoeff_counter < 1  and T3bList.size()>0)
        {
          std::cout << "@@@@@@@@@@@@  channel = (" << la << " " << j2a << " , " << lb << " " << j2b << " , " << lc << " " << j2c << std::endl;
        }

        IMSRGProfiler::timer[std::string(__func__)+"_iMonLoop"] += omp_get_wtime() - t_internal;




      }// for ilj_c
    }// for ilj_b
  }// for ilj_a















/*

  for (auto& obc_a : hf.modelspace->OneBodyChannels )
  {
    int la = obc_a.first[0];
    int j2a = obc_a.first[1];
    int tz2a = obc_a.first[2];
    for (auto& obc_b : hf.modelspace->OneBodyChannels )
    {
      int lb = obc_b.first[0];
      int j2b = obc_b.first[1];
      int tz2b = obc_b.first[2];
      int Jab_min = std::abs(j2a-j2b)/2;
      int Jab_max = (j2a+j2b)/2;
      int Tzab = (tz2a+tz2b)/2;
      int Tab_min = std::abs(Tzab);
       // possibilities  (tz2a,tz2b,Tab)=>CG :  (1,1,1)=>1 ,  (1,-1,1)=>sqrt(2), (-1,1,1)=>sqrt(2),  (1,-1,0)=>sqrt(2) , (-1,1,0)=>-sqrt(2)
      std::array<double,2> isospin2_Clebsch = {tz2a*sqrt(0.5), 0.5*(tz2a+tz2b) + std::abs(tz2a-tz2b)/2*sqrt(0.5) };
      for (auto& obc_c : hf.modelspace->OneBodyChannels )
      {
        int lc = obc_c.first[0];
        int j2c = obc_c.first[1];
        int tz2c = obc_c.first[2];
        std::cout << "OneBody Channels: (" << la << " " << j2a << " " << tz2a << ") , (" << lb << " " << j2b << " " << tz2b << ") , (" << lc << " " << j2c << " " << tz2c << ")"
                  << "   monopoles computed = " << n_mon << std::endl;
        int twoTz = tz2a + tz2b + tz2c;
        int twoT_min = std::abs( twoTz ); // this will either be 1 or 3
        // index by  2*Tab + T/2 ->             { (0,1), (0,3), (1,1), (1,3) }
        std::array<double,4> isospin3_Clebsch = {   1.0,   0.0,  AngMom::CG(1,Tzab,0.5,0.5*tz2c, 0.5, 0.5*twoTz),  AngMom::CG(1,Tzab,0.5,0.5*tz2c, 1.5,0.5*twoTz) };
        // possibilities (Tab,Tzab, 2tc, 2tzc| 2T)=>CG :  (0,0, 1,1|1)=>1 ,  (0,0, 1,-1|1)=>1,     (1,1, 1,1| 3,3)=>1,  (1,-1,1,-1|3,-3)=>1,   (1,0, 1,1|3,1)
        // compute all the relevant T coefficients.  The coefficients have the following indices:
        //  na, nb, nc, Jab, J, jacobi1, jacobi2, J12, Lcm   (Ncm can be inferred from energy conservation?)
//        std::unordered_map<uint64_t,double> T3bList;  // we will put things into a hash table.
        std::unordered_map<std::string,double> T3bList;  // we will put things into a hash table.

        double t_internal = omp_get_wtime();
        // First, generate all the keys that we'll need, but don't compute the T coefficients
        for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
        {
          int twoJ_min = std::abs(2*Jab-j2c);
          int twoJ_max = 2*Jab+j2c;
          for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
          {
            for (auto a : obc_a.second) // run over the list of orbit indices in channel obc_a
            {
             Orbit& oa = hf.modelspace->GetOrbit( a );
             for (auto b : obc_b.second)
             {
              if (a==b and Jab%2>0 ) continue;
              Orbit& ob = hf.modelspace->GetOrbit( b );
              for (auto c : obc_c.second)
              {
                if (a==b and a==c and j2a<3 ) continue;
                Orbit& oc = hf.modelspace->GetOrbit( c );
                int Eabc = 2*(oa.n+ob.n+oc.n) + oa.l+ob.l+oc.l;
                if (Eabc>Nmax) continue;

                for (int jac1=0; jac1<jacobi_1.size(); jac1++)
                {
                 int N1 = jacobi_1[jac1].n;
                 int L1 = jacobi_1[jac1].l;
                 int J1 = jacobi_1[jac1].j;
                 for (int jac2=0; jac2<jacobi_2.size(); jac2++)
                 {
                   int N2 = jacobi_2[jac2].n;
                   int L2 = jacobi_2[jac2].l;
                   int twoJ2 = jacobi_2[jac2].j2;
                   if ( (2*(N1+N2)+L1+L2) > Eabc ) continue;
                   int twoJ12_min = std::abs(2*J1-twoJ2);
                   int twoJ12_max = 2*J1+twoJ2;
                   for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
                   {
                    for (int Lcm=(L1+L2+Eabc)%2; Lcm<=(Eabc-2*(N1+N2)-L1-L2); Lcm++)
                    {
                      if ( (std::abs(2*Lcm - twoJ12)>twoJ) or (2*Lcm+twoJ12)<twoJ ) continue; // check triangle condition for J12,Lcm J
                      auto hash_key = TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,jac1,jac2,twoJ12,Lcm);
                      T3bList[hash_key] = 0;
                    }// for Lcm
                   } // for twoJ12
                 } // for jac2
               } // for jac1
              } // for c
             } // for b
            } // for a
          } // for twoJ
        } // for Jab

        IMSRGProfiler::timer[std::string(__func__)+"_GenerateKeys"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
        std::cout << "There are " << T3bList.size() << " Tcoefficients to compute in this channel" << std::endl;

//        std::cout << "done. Now compute the T coefficients" << std::endl;

        // Next, we compute all the T coefficients. We'll probably want to do this in parallel.
        // In the current construction, all threads loop through the entire set of elements, but
        // they only stop to compute when it's their turn.
        #pragma omp parallel
        {
          size_t cnt = 0;
          int ithread = omp_get_thread_num();
          int nthreads = omp_get_num_threads();
          for(auto element = T3bList.begin(); element !=T3bList.end(); ++element, cnt++)
          {
            if(cnt%nthreads != ithread) continue; // Check if this thread should compute this element
            auto hash_key = element->first;
            int na,nb,nc,Jab,twoJ,jac1,jac2,twoJ12,Lcm;
            TcoeffUnHash(hash_key, na,nb,nc,Jab,twoJ,jac1,jac2,twoJ12,Lcm);
            auto& jacobi1 = jacobi_1[jac1];
            auto& jacobi2 = jacobi_2[jac2];
            int Ncm = (2*(na+nb+nc)+la+lb+lc - 2*(jacobi1.n+jacobi2.n) - jacobi1.l - jacobi2.l - Lcm)/2;
            if (Ncm<0) continue;
            double tcoef = AngMom::Tcoeff( na, la, j2a, nb, lb, j2b, nc, lc, j2c, Jab, twoJ, jacobi1.n, jacobi1.l, jacobi1.s, jacobi1.j, jacobi2.n, jacobi2.l, jacobi2.j2, twoJ12, Ncm, Lcm);
            element->second = tcoef;
          }
        }

        IMSRGProfiler::timer[std::string(__func__)+"_ComputeToeff"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();

        std::vector<size_t> imonlist_all;
        std::vector<size_t> imonlist;
        std::vector<size_t> imonlist_swap_ab;
        std::vector<size_t> imonlist_swap_de;
        std::vector<size_t> imonlist_swap_abde;
//        for (size_t imon=0; imon<v3mon.size(); imon++)
        for (size_t imon=0; imon<hf.Vmon3_keys.size(); imon++)
        {
          // Check if the given key corresponds to the set of one-body channels we're currently working on.
          auto key = hf.Vmon3_keys[imon];
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          if ( std::find( obc_a.second.begin(), obc_a.second.end(), a ) == obc_a.second.end() ) continue;
          if ( std::find( obc_b.second.begin(), obc_b.second.end(), b ) == obc_b.second.end() ) continue;
          if ( std::find( obc_c.second.begin(), obc_c.second.end(), c ) == obc_c.second.end() ) continue;
//          imonlist.push_back( imon );
          imonlist_all.push_back( imon );
          if (a<=b and d<=e) imonlist.push_back(imon);
//          else if (a>=b and d<=e) imonlist_swap_ab.push_back(imon);
//          else if (a<=b and d>=e) imonlist_swap_de.push_back(imon);
//          else if (a>=b and d>=e) imonlist_swap_abde.push_back(imon);
//          std::cout << " imon = " << imon << "   =>  " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
        }
        // This is fast, so we can still afford to do this twice without parallelization
        // Also, we're now working with a greatly reduced subset of the monopole terms

        imonlist_swap_ab.resize( imonlist.size(),-1 );
        imonlist_swap_de.resize( imonlist.size(),-1 );
        imonlist_swap_abde.resize( imonlist.size(),-1 );
//        std::cout << "Begin imonlistloop. imonlist_all:" << std::endl;
//        for ( auto imon : imonlist_all ) std::cout << imon << " ";
//        std::cout << std::endl;
//        for (auto imon : imonlist )
        for (int ind=0; ind<imonlist.size(); ind++ )
        {
          auto imon = imonlist[ind];
          auto key = hf.Vmon3_keys[imon];
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          auto key_ab = hf.Vmon3Hash( b,a,c,d,e,f);
          auto key_de = hf.Vmon3Hash( a,b,c,e,d,f);
          auto key_abde = hf.Vmon3Hash( b,a,c,e,d,f);
          for ( auto iall : imonlist_all )
          {
            if (hf.Vmon3_keys[iall] == key_ab) imonlist_swap_ab[ind] = iall;
            if (hf.Vmon3_keys[iall] == key_de) imonlist_swap_de[ind] = iall;
            if (hf.Vmon3_keys[iall] == key_abde) imonlist_swap_abde[ind] = iall;
          }
//          std::cout << " ind = " << ind << " : " << imon << "  " << imonlist_swap_ab[ind] << "  " << imonlist_swap_de[ind] << "  " << imonlist_swap_abde[ind] << std::endl;
//          std::cout << "        " << a << " " << b << " " << c << "  " << d << " " << e << " " << f << std::endl;
//          std::cout << "        " << key << "  " << key_ab << " " << key_de << " " << key_abde << std::endl;
//          int aa,bb,cc,dd,ee,ff;
//          hf.Vmon3UnHash( key_ab ,aa,bb,cc,dd,ee,ff);
//          std::cout << "     " << aa << " " << bb << " " << cc << " " << dd << " " <<ee << " " << ff << "   ";
//          hf.Vmon3UnHash( key_de ,aa,bb,cc,dd,ee,ff);
//          std::cout << "     " << aa << " " << bb << " " << cc << " " << dd << " " <<ee << " " << ff << "   ";
//          hf.Vmon3UnHash( key_abde ,aa,bb,cc,dd,ee,ff);
//          std::cout << "     " << aa << " " << bb << " " << cc << " " << dd << " " <<ee << " " << ff << "   " << std::endl;
        }
        std::cout << "done. size of swapab = " << imonlist_swap_ab.size() << std::endl;

        IMSRGProfiler::timer[std::string(__func__)+"_FindMonKeys"] += omp_get_wtime() - t_internal;
        t_internal = omp_get_wtime();
        bool verbose = false;
//        if (lb==3 and lc==3) verbose = true;

        if (verbose) std::cout << "done.  Now the loop over imon... Size of imonlist is " << imonlist.size() << std::endl;
        // And at last, we compute all the monopole terms for this combination of one-body channels
        // This can also be done in parallel
//        for (size_t imon=0; imon<v3mon.size(); imon++)
//        #pragma omp parallel for schedule(dynamic,1)
        #pragma omp parallel for schedule(dynamic,1)
        for (size_t ilist=0; ilist<imonlist.size(); ilist++)
        {
          size_t imon = imonlist[ilist];
          size_t imon_ab = imonlist_swap_ab[ilist];
          size_t imon_de = imonlist_swap_de[ilist];
          size_t imon_abde = imonlist_swap_abde[ilist];
          auto key = hf.Vmon3_keys[imon];
//          if (verbose) std::cout << " imon,imon_ab,imon_de,imon_abde = " << imon << " " << imon_ab << " " << imon_de << " " << imon_abde << std::endl;
          int a,b,c,d,e,f;
          hf.Vmon3UnHash( key, a, b, c, d, e, f);  // static method so we coud call it without a class instance if we wanted...
          // V3mon(a,b,c,d,e,f) = 1/(2jc+1) * sum_{Jab,J} (2J+1) <abc Jab,J| V |def Jab,J>
          //                    = 1/(2jc+1) * sum_{Jab,J} (2J+1) sum_{jac1,jac2,J12,jac1',jac2',J12',Ncm,Lcm}  Tcoef(abc,Jab,J;jac1,jac2,J12,Ncm,Lcm) * Tcoef(def,Jab,J;jac1',jac2',J12',Ncm,Lcm)
          //                                                                                 sum_{T12,T12',T}  * CG(ta,tb,T12)*CG(tc,T12,T) * CG(td,te,T12')*CG(df,T12',T)   <- Here T means isospin
          //                                                                                                   * < jac1,jac2,J12,T12 | V | jac1',jac2',J12',T12' >
          if (a>b) continue;
          if (d>e) continue;
          Orbit& oa = hf.modelspace->GetOrbit(a);
          Orbit& ob = hf.modelspace->GetOrbit(b);
          Orbit& oc = hf.modelspace->GetOrbit(c);
          Orbit& od = hf.modelspace->GetOrbit(d);
          Orbit& oe = hf.modelspace->GetOrbit(e);
          Orbit& of = hf.modelspace->GetOrbit(f);
          int Eabc = 2*(oa.n+ob.n+oc.n) + oa.l+ob.l+oc.l;
          int Edef = 2*(od.n+oe.n+of.n) + od.l+oe.l+of.l;
          if (Eabc>Nmax or Edef>Nmax) continue;
          int parity = Eabc%2;
          double v_monopole = 0;
          double v_monopole_swap_ab = 0;
          double v_monopole_swap_de = 0;
          int Jab_step = (a==b or d==e)? 2 : 1;
          for (int Jab=Jab_min; Jab<=Jab_max; Jab+=Jab_step)
          {
           if ( (a==b and a==c and j2a<3) or (d==e and d==f and od.j2<3) ) continue; // can't fit 3 identical particles in a j=1/2 orbit
           int twoJ_min = std::abs(2*Jab-j2c);
           int twoJ_max = 2*Jab+j2c;
           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
           {
             for (int Tab=Tab_min; Tab<=1; Tab++)
             {
              int ab_swap_phase = - AngMom::phase( (j2a+j2b)/2 - Jab + 1 - Tab); // don't forget that minus sign out in front...
              double isoClebsch_ab = isospin2_Clebsch[Tab];
              for (int Tde=Tab_min; Tde<=1; Tde++)
              {
                double isoClebsch_de = isospin2_Clebsch[Tde];
                for (int twoT=twoT_min; twoT<=std::min(2*Tab+1,2*Tde+1); twoT+=2)
                {
                  int de_swap_phase = - AngMom::phase( (j2a+j2b)/2 - Jab + 1 - Tde); // don't forget that minus sign out in front...
                  double isoClebsch_c = isospin3_Clebsch[2*Tab + twoT/2];
                  double isoClebsch_f = isospin3_Clebsch[2*Tde + twoT/2];
                  if (std::abs(isoClebsch_c)<1e-6 or std::abs(isoClebsch_f)<1e-6) continue;
                  double v_sumJT = 0;
                  for (int Ecm=0; Ecm<=std::min(Eabc,Edef); Ecm++)
                  {
                   int E12abc = Eabc-Ecm;
                   int E12def = Edef-Ecm;
                   if (E12abc > Nmax  or E12def>Nmax) continue;
                   for (int Lcm=Ecm%2; Lcm<=Ecm; Lcm+=2)
                   {
                     int twoJ12_min = std::abs(twoJ-2*Lcm);
                     int twoJ12_max = std::min( twoJ+2*Lcm, twoJmax);
                     for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2)
                     {
                       if (verbose)   std::cout << "twoJ,Lcm,twoJ12,twoJmax = " << twoJ << " " << Lcm << " " << twoJ12 << " " << twoJmax << std::endl;
                       if (verbose)   std::cout << "Ecm,E12abc,E12def " << Ecm << " " << E12abc << " " << E12def  << std::endl;
                       auto hashTJN_abc = HashTJN(twoT,twoJ12,E12abc);
                       auto hashTJN_def = HashTJN(twoT,twoJ12,E12def);
                       if (verbose) std::cout << "hash_abc = " << hashTJN_abc << " hash_def = " << hashTJN_def << "  size of dimNAS = " << dimensionNAS.size() << std::endl;
                       size_t dimNAS_abc = GetDimensionNAS( twoT, twoJ12, parity, E12abc ); 
                       size_t dimNAS_def = GetDimensionNAS( twoT, twoJ12, parity, E12def ); 

                       size_t startloc = GetStartLocNAS(twoT, twoJ12, E12abc, E12def) ;
//                       size_t otherstartloc = GetStartLocNAS(twoT, twoJ12, E12def, E12abc) ;

//                       std::cout << "   dimensions: " << dimNAS_abc << " , " << dimNAS_def << std::endl;
                       // As a reminder, armadillo stores matrices   [ M11  M12  M13 ]
                       // in column-major order, as in               | M21  M22  M23 |
                       //                                            [ M31  M32  M33 ]
                       //
                       //  [ Tdef_1  Tdef_2 ... ]  *  [ <def|V|abc>  ... ]  *  [ Tabc_1 ]
                       //                             [   ...   ...  ... ]     [ Tabc_2 ]
                       //                             [   ...   ...  ... ]     [  ...   ]

                       arma::vec Tabc(dimNAS_abc, arma::fill::zeros);
                       arma::rowvec Tdef(dimNAS_def, arma::fill::zeros);
                       arma::mat matelNAS( &meNAS[startloc], dimNAS_def, dimNAS_abc, false );
//                       std::cout << "Getting matrix beginning at " << startloc << "  ( should it be " << otherstartloc << " ? ) " <<  std::endl;
                       if (verbose) std::cout << "Getting matrix beginning at " << startloc << "   size of vector = " << meNAS.size() <<  std::endl;
                       for (size_t iNAS_abc = 0; iNAS_abc<dimNAS_abc; iNAS_abc++)
                       {
                         auto& index_1_2_abc = NAS_jacobi_states.at(hashTJN_abc).at(iNAS_abc);
                         if ( jacobi_1.at(index_1_2_abc[0]).t != Tab ) continue;
                         Tabc[iNAS_abc] = T3bList[ TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,index_1_2_abc[0],index_1_2_abc[1],twoJ12,Lcm )];
                       }
                       for (size_t iNAS_def = 0; iNAS_def<dimNAS_def; iNAS_def++)
                       {
                         auto& index_1_2_def = NAS_jacobi_states.at(hashTJN_def).at(iNAS_def);
                         if ( jacobi_1.at(index_1_2_def[0]).t != Tde ) continue;
                         Tdef[iNAS_def] = T3bList[ TcoeffHash(od.n,oe.n,of.n,Jab,twoJ,index_1_2_def[0],index_1_2_def[1],twoJ12,Lcm )];
                       }

                       if (verbose) std::cout << "Doing the mat mult" << std::endl;
                       arma::mat result = Tdef * matelNAS * Tabc;
//                       arma::mat backresult = Tabc.t() * matelNAS * Tdef.t();
                       v_sumJT += result[0];
//                       std::cout << "matrix multiply version: " << std::endl << Tdef << std::endl << std::endl << matelNAS << std::endl << std::endl << Tabc << std::endl << std::endl << std::endl;
//                       std::cout << "matelNAS * Tabc :" << std::endl << ( matelNAS * Tabc ) << std::endl;
//                       std::cout << "result = " << result[0] << "    or " << backresult[0] << std::endl;
                  

                       
//                       for (size_t iNAS_abc = 0; iNAS_abc<dimNAS_abc; iNAS_abc++)
//                       {
//                         auto& index_1_2_abc = NAS_jacobi_states.at(hashTJN_abc).at(iNAS_abc);
//                         auto indx_jac1_abc = index_1_2_abc[0];
//                         auto indx_jac2_abc = index_1_2_abc[1];
//                         jacobi1_state& jac1_abc = jacobi_1.at(indx_jac1_abc);
//                         jacobi2_state& jac2_abc = jacobi_2.at(indx_jac2_abc);
////                         GetJacobiStates( twoT, twoJ, parity, E12abc, iNAS_abc, jac1_abc, jac2_abc);
//                         if (jac1_abc.t != Tab) continue;
//
//                         for (size_t iNAS_def = 0; iNAS_def<dimNAS_def; iNAS_def++)
//                         {
//                           auto& index_1_2_def = NAS_jacobi_states.at(hashTJN_def).at(iNAS_def);
//                           auto indx_jac1_def = index_1_2_def[0];
//                           auto indx_jac2_def = index_1_2_def[1];
//                           jacobi1_state& jac1_def = jacobi_1.at(indx_jac1_def);
//                           jacobi2_state& jac2_def = jacobi_2.at(indx_jac2_def);
//
////                           GetJacobiStates( twoT, twoJ, parity, E12def, iNAS_def, jac1_def, jac2_def);
//                           if (jac1_def.t != Tde) continue;
//                           if ( std::abs(twoJ12-2*Lcm)>twoJ or (twoJ12+2*Lcm)<twoJ ) continue;
//                           auto THash_abc = TcoeffHash(oa.n,ob.n,oc.n,Jab,twoJ,indx_jac1_abc,indx_jac2_abc,twoJ12,Lcm );
//                           auto THash_def = TcoeffHash(od.n,oe.n,of.n,Jab,twoJ,indx_jac1_def,indx_jac2_def,twoJ12,Lcm );
//                           double Tcoef_abc = T3bList[THash_abc];
//                           double Tcoef_def = T3bList[THash_def];
//                           if ( std::abs(Tcoef_abc*Tcoef_def)>1e-6)
//                           {
////                             std::cout << "        Tcoefs = " << Tcoef_abc << " " << Tcoef_def << std::endl;
////                             std::cout << "       accesing " << iNAS_abc << " " << iNAS_def << " " << E12abc << " " << E12def << " " << twoT << " " << twoJ << " " << parity << std::endl;
////                             std::cout << "       Tab,Tde,twoT = " << Tab << " " << Tde << " " << twoT << "   isospin clebsch's " << isoClebsch_ab << " " << isoClebsch_c << "   " << isoClebsch_de << " " << isoClebsch_f << std::endl;
//                             double v_NAS = GetMatElNAS( iNAS_abc, iNAS_def,  E12abc,  E12def, twoT, twoJ12,  parity);
////                             v_sumJT +=  Tcoef_abc * Tcoef_def * v_NAS;
////                             std::cout << "  " << iNAS_abc << " " << iNAS_def << "  " << Tcoef_abc << " " << Tcoef_def << "  " << v_NAS << "   " <<  Tcoef_abc * Tcoef_def * v_NAS << "  sum = " << v_sumJT << std::endl;
////                             std::cout << "         v_NAS = " << v_NAS << "  v_sumJT = " << v_sumJT << std::endl;
//                           }
//                         } // for iNAS_def
//                       } // for iNAS_abc
                     } // for twoJ12
                   } // for Lcm
                  } // for Ecm
                  // v_sumJT is equal to <abc Jab Tab JT | V | def Jde Tde JT>
                  Ket3 bra = Ket3(oa,ob,oc);
                  Ket3 ket = Ket3(od,oe,of);
                  if (verbose) std::cout << "Assigning" << std::endl;
//                  double vcheck_JT =  GetLabMatEl( bra, ket, Jab, Jab, twoJ, Tab, Tde, twoT);
                  double vterm = 6 * (twoJ+1) * isoClebsch_ab * isoClebsch_c * isoClebsch_de * isoClebsch_f * v_sumJT;
                  v_monopole += vterm;
                  v_monopole_swap_ab += ab_swap_phase * vterm;
                  v_monopole_swap_de += de_swap_phase * vterm;
//                  std::cout << "    vmonopole: += " << twoJ+1 << " * " << isoClebsch_ab << " * " << isoClebsch_c << " * " << isoClebsch_de << " * " << isoClebsch_f << " * " << 6*v_sumJT
//                            << "  ( " << vcheck_JT << " )   " << "  => " << v_monopole << std::endl;
                } // for twoT
              } // for Tde
             } // for Tab
           } // for twoJ
          } // for Jab
          v_monopole /= j2c+1.;
          v_monopole_swap_ab /= j2c+1.;
          v_monopole_swap_de /= j2c+1.;
         // There are likely some symmetries we can exploit here to avoid redundant calculations
          hf.Vmon3[imon] = v_monopole; // put that bad boy into the output vector
          if (imon_abde !=-1 and imon_abde != imon)
            hf.Vmon3[imon_abde] = v_monopole; // put that bad boy into the output vector
          if (imon_ab !=-1)
            hf.Vmon3[imon_ab] = v_monopole_swap_ab; // put that bad boy into the output vector
          if (imon_de !=-1)
            hf.Vmon3[imon_de] = v_monopole_swap_de; // put that bad boy into the output vector
          if (verbose)
          {
            std::cout << "abcdef = " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
            std::cout << "SETTING " << imon << " -> " << v_monopole << std::endl << std::endl;
            if ( imon_abde !=-1 and imon_abde != imon )
            std::cout << "Also SETTING " << imon_abde << " -> " << v_monopole << std::endl;
            if ( imon_ab !=-1 and imon_ab != imon )
            std::cout << "Also SETTING " << imon_ab << " -> " << v_monopole_swap_ab << std::endl;
            if ( imon_de !=-1 and imon_de != imon )
            std::cout << "Also SETTING " << imon_de << " -> " << v_monopole_swap_de << std::endl;

          }
        } // for imon
        n_mon += imonlist.size();
//        std::cout << "done " << std::endl;

        IMSRGProfiler::timer[std::string(__func__)+"_iMonLoop"] += omp_get_wtime() - t_internal;
      } // for obc_c
    } // for obc_b
  } // for obc_a

  // all done!
*/

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
     if ( size_t(bk.matstart) != start-1)
     {
       std::cout << __func__ << "  DANGER!! start != matstart : " << start << " != " << bk.matstart << std::endl;
     }
     if (dim != size_t(bk.nspnum * bk.nrelnum) )
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
     for (size_t i=0; i<dim; i++)
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
         std::string good1 = std::abs(mytc - tc)<1e-6 ? "good" : "BADval";
         if (good1=="BADval" and  std::abs(mytc + tc)<1e-6 )  good1 = "BADphase";
         std::string good2 = std::abs(mytc - mytc_bf)<1e-6 ? "good" : "BAD";
//   double Tcoeff( int na, int la, int j2a, int nb, int lb, int j2b, int nc, int lc, int j2c, int Jab, int twoJ, int N1, int L1, int S1, int J1, int N2, int L2, int twoJ2, int twoJ12, int Ncm, int Lcm);
         std::cout << "( " << lab.na << " " << lab.la << " " << lab.ja << " " <<  lab.nb << " " << lab.lb << " " << lab.jb << " " <<  lab.nc << " " << lab.lc << " " << lab.jc << " " <<  lab.bj12 << " " <<  jtot << " "
                           <<  rel.n12 << " " << rel.l12 << " " << rel.s12 << " " <<  rel.j12 << " " <<  rel.n3 << " " <<  rel.l3 << " " <<  rel.I3 << " " <<  rel.j3 << " " <<  rel.ncm << " " <<  rel.lcm << " ) " << std::endl;
         std::cout << "   " << isp << " " << irel << "  " << index << " ( " << tcoeff.size() << " )  -> "
                   << std::fixed << std::setw(12) << tc << "   "
                   << std::fixed << std::setw(12) << mytc << "   "
                   << std::fixed << std::setw(12) << mytc_bf << "  "
                   << std::fixed << std::setw(8) << good1
                   << std::fixed << std::setw(8) << good2
                   << std::endl;
       }
     }
   }
  }

  std::cout << "All done. Exiting nicely. " << std::endl;

}









