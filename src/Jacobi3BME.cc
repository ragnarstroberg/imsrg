
#include "Jacobi3BME.hh"
#include "AngMom.hh"
#include <istream>
#include <iomanip>




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
  return   (twoT-twoTmin)/2 * ( ((twoJmax-twoJmin)/2*(Nmax/2*(Nround+1))+ Nround+1 ) + Nmax/2*(Nround+1) + Nround + 1 )
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
   auto hash = HashTJNN(twoT,twoJ,Nbra,Nket);
   size_t dimket = GetDimensionNAS(twoT,twoJ,p,Nket); 
   size_t start_loc = GetStartLocNAS(twoT,twoJ,Nbra,Nket);
   return meNAS.at( start_loc + ibra*dimket + iket );
}

void Jacobi3BME::SetMatElNAS(size_t ibra, size_t iket, int Nbra, int Nket, int twoT, int twoJ, int p, double matel)
{
   auto hash = HashTJNN(twoT,twoJ,Nbra,Nket);
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
  std::cout << std::endl << " Computing NAS matrix elements " << std::endl << std::endl;
  // T2,J2,parity are conserved by V
  for (int twoT=twoTmin; twoT<=twoTmax; twoT+=2)
  {
    for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
    {
      for (int parity=0; parity<=1; parity++)
      {
        std::cout << "TJP = " << twoT << " " << twoJ << " " << parity << std::endl;

        for (int Nbra=parity; Nbra<=Nmax; Nbra+=2)  
        {
          size_t dim_braAS = GetDimensionAS(twoT,twoJ,parity,Nbra);
          size_t dim_braNAS = GetDimensionNAS(twoT,twoJ,parity,Nbra);
          size_t cfp_begin_bra = GetCFPStartLocation(twoT,twoJ,Nbra) ;
          if (dim_braAS==0) continue;
          arma::mat cfp_bra( &(cfpvec[cfp_begin_bra]), dim_braNAS, dim_braAS, /*copy_aux_mem*/ true);

          for (int Nket=parity; Nket<=Nmax; Nket+=2)  
          {
            std::cout << "   Nbra,Nket = " << Nbra << " " << Nket << std::endl;
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
            std::cout << "ASmat: " << std::endl << ASmat << std::endl << std::endl;
            std::cout << "NASmat = " << std::endl << NASmat << std::endl << std::endl;

            for (size_t iNAS=0; iNAS<dim_ketNAS*dim_braNAS; iNAS++)
            {
              meNAS[startNAS+iNAS] = NASmat(iNAS); // single index assumes a flat layout with column-major ordering
            }
//            if (dim_braNAS>1 and dim_ketNAS>1)
//            {
//              std::cout << "Check: 1,1 element : " << ElementNAS(1,1,Nbra,Nket,T2,J2,parity) << std::endl;;
//            }
          }
        }
      }
    }
  }

}







/// Time to do some heavy lifting with T coefficients
double Jacobi3BME::GetLabMatEl( Ket3& bra, Ket3& ket, int Jab, int Jde, int twoJ, int Tab, int Tde, int twoT)
{
  double me_lab = 0;
//  int LabE3max = jacobi_basis.Nmax;
  int LabE3max = 16;

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
 //          std::cout << "  ibraNAS,iketNAS = " << ibraNAS << " " << iketNAS << std::endl;
           GetJacobiStates( twoT, twoJ12, parity12_ket, E12_ket, iketNAS, jac1_ket, jac2_ket);
//           std::cout << "     checking isospin ket: " << jac1_ket.t << " " << Tde << std::endl;
           if ( jac1_ket.t != Tde ) continue;
//           double Tcoeff_ket = Tcoeff_wrapper( ket, Jde, twoJ, jac1_ket, jac2_ket, twoJ12, Ncm, Lcm);

//           std::cout << "    ** finding the NAS matrix element. T = " << T << "   J = " << J12 << "  p = " << parity12_bra << std::endl;

//           double meNAS = ElementNAS( ibraNAS, iketNAS, E12_bra, E12_ket, twoT, twoJ12, parity12_bra );
           double meNAS = GetMatElNAS( ibraNAS, iketNAS, E12_bra, E12_ket, twoT, twoJ12, parity12_bra );

//           me_lab +=  Tcoeff_bra * meNAS * Tcoeff_ket  ;
           me_lab +=  Tcoeff_bra[ibraNAS] * meNAS * Tcoeff_ket[iketNAS]  ;
//           std::cout << " " << ibraNAS << " , " << iketNAS << ":  (( " << Tcoeff_bra[ibraNAS] << "  " << meNAS << "  " << Tcoeff_ket[iketNAS] << " ))   => sum_ME = " << me_lab << std::endl;
//           std::cout << "    < N1=" << jac1_bra.n << " L1=" << jac1_bra.l << " S1=" << jac1_bra.s << " J1=" << jac1_bra.j << " T1=" <<jac1_bra.t << ",  N2=" << jac2_bra.n << " L2=" << jac2_bra.l << " T=" << twoT << " J12=" << twoJ12 << " | ..." << std::endl;

          }// for ibraNAS
        } // for iketNAS
      } // for J12
    } // for Lcm
  } // for Ecm


  return 6 * me_lab;
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











