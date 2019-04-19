
#include "Jacobi3BME.hh"
#include "AngMom.hh"




Jacobi3BME::Jacobi3BME()
{}

Jacobi3BME::Jacobi3BME( int nmax, int twojmin, int twojmax, int twotmin, int twotmax )
 : Nmax(nmax), twoJmin(twojmin), twoJmax(twojmax), twoTmin(twotmin), twoTmax(twotmax)
{
///  Allocate   // Don't call allocate upon construction because we need to read the relevant dimensions from file
}
 

//////  Layout: {   (JTp0) [N0,N0]  [N0,N1]  [N0,N2] ... [N1,N1]  [N1,N2] ... |       (JTp1) ...    |     (JTp2) ... }
//////  Here [N0,N1] is a matrix with dimension dim(T,J,N0) x dim(T,J,N1) 
//////  We store the full symmetrix NN matrix. It's only a factor of 2, and it will avoid a bunch of comparisons later.
//////
void Jacobi3BME::Allocate()
{
  size_t hashmax = HashTJNN(twoTmax, twoJmax, Nmax, Nmax);
  start_locAS.resize(  hashmax+1, 0 ); // resize based on the largest hash we can have
  start_locNAS.resize( hashmax+1 ,0 ); 
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
double Jacobi3BME::GetMatElAS(size_t ibra, size_t iket, int Nbra, int Nket, int T2, int J2, int p)
{
   size_t dimket = GetDimensionAS(T2,J2,p,Nket); 
   size_t start_loc = GetStartLocAS(T2,J2,Nbra,Nket);
   return meAS.at( start_loc + ibra*dimket + iket );
}

// Access an antisymmetrized matrix element
void Jacobi3BME::SetMatElAS(size_t ibra, size_t iket, int Nbra, int Nket, int T2, int J2, int p, double matel)
{
   size_t dimket = GetDimensionAS(T2,J2,p,Nket); 
   size_t start_loc = GetStartLocAS(T2,J2,Nbra,Nket);
   meAS.at( start_loc + ibra*dimket + iket ) = matel;
}


double Jacobi3BME::GetMatElNAS(size_t ibra, size_t iket, int Nbra, int Nket, int T2, int J2, int p)
{
   auto hash = HashTJNN(T2,J2,Nbra,Nket);
   size_t dimket = GetDimensionNAS(T2,J2,p,Nket); 
   size_t start_loc = GetStartLocNAS(T2,J2,Nbra,Nket);
   return meNAS.at( start_loc + ibra*dimket + iket );
}

void Jacobi3BME::SetMatElNAS(size_t ibra, size_t iket, int Nbra, int Nket, int T2, int J2, int p, double matel)
{
   auto hash = HashTJNN(T2,J2,Nbra,Nket);
   size_t dimket = GetDimensionNAS(T2,J2,p,Nket); 
   size_t start_loc = GetStartLocNAS(T2,J2,Nbra,Nket);
   meNAS.at( start_loc + ibra*dimket + iket ) = matel;
}




void Jacobi3BME::GetJacobiStates( int T2, int J2, int parity, int E12, int iNAS, jacobi1_state& jac1, jacobi2_state& jac2)
{
  size_t hash = HashTJN(T2,J2,E12);
  auto& index_1_2 = NAS_jacobi_states.at(hash).at(iNAS);  // for each T,J,p,N block, this maps the NAS index to a pair of jacobi states (actually, to their index)
  jac1 = jacobi_1.at(index_1_2[0]);
  jac2 = jacobi_2.at(index_1_2[1]);

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
  for (int T2=twoTmin; T2<=twoTmax; T2+=2)
  {
    for (int J2=twoJmin; J2<=twoJmax; J2+=2)
    {
      for (int parity=0; parity<=1; parity++)
      {
        std::cout << "TJP = " << T2 << " " << J2 << " " << parity << std::endl;

        for (int Nbra=parity; Nbra<=Nmax; Nbra+=2)  
        {
          size_t dim_braAS = GetDimensionAS(T2,J2,parity,Nbra);
          size_t dim_braNAS = GetDimensionNAS(T2,J2,parity,Nbra);
          size_t cfp_begin_bra = GetCFPStartLocation(T2,J2,Nbra) ;
          if (dim_braAS==0) continue;
          arma::mat cfp_bra( &(cfpvec[cfp_begin_bra]), dim_braNAS, dim_braAS, /*copy_aux_mem*/ true);

          for (int Nket=parity; Nket<=Nmax; Nket+=2)  
          {
            std::cout << "   Nbra,Nket = " << Nbra << " " << Nket << std::endl;
            // build matrix for AS jacobi MEs, CFPs, and  (as the output) NAS jacobi MEs
            size_t dim_ketAS = GetDimensionAS(T2,J2,parity,Nket);
            if (dim_ketAS==0) continue;
         
            size_t dim_ketNAS = GetDimensionNAS(T2,J2,parity,Nket);
            size_t cfp_begin_ket = GetCFPStartLocation(T2,J2,Nket) ;
            arma::mat cfp_ket( &(cfpvec[cfp_begin_ket]), dim_ketNAS, dim_ketAS, /*copy_aux_mem*/ true);
 
            size_t startAS = GetStartLocAS(T2, J2, Nbra, Nket);
            size_t startNAS = GetStartLocNAS(T2, J2, Nbra, Nket);
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

  std::cout << "In GetLabMatEl. Ebra,Eket = " << Ebra << " " << Eket << std::endl;
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
    for (int Lcm=(Ecm%2); Lcm<=Ecm; Lcm ++)
    {
      int Ncm = (Ecm-Lcm)/2; // Ncm is the number of radial nodes in the CM wave function. Ecm = 2Ncm + Lcm (but the variable Lcm is twice the actual value).
      int twoJ12_min = std::max(twoJmin, std::abs( twoJ - 2*Lcm ));
      int twoJ12_max = std::min(twoJmax, twoJ + 2*Lcm) ;

      for (int twoJ12=twoJ12_min; twoJ12<=twoJ12_max; twoJ12+=2) // the interaction conserves J12, parity, and T12=twoT
      {
        int NASdim_bra = GetDimensionNAS( twoT, twoJ12, parity12_bra, E12_bra ); 
        int NASdim_ket = GetDimensionNAS( twoT, twoJ12, parity12_ket, E12_ket ); 
        std::cout << "dimensions: " << NASdim_bra << " " << NASdim_ket << std::endl;
//           size_t dimbra = jacobi_basis.GetDimensionNAS(T,J,parity12_bra,E12_bra); 
//           std::cout << "dimbra = " <<  dimbra << std::endl;
        for (int ibraNAS=0; ibraNAS<NASdim_bra; ibraNAS++)
        {
          jacobi1_state jac1_bra,jac1_ket;
          jacobi2_state jac2_bra,jac2_ket;

//          std::cout << "about to get jac1_bra and jac2_bra" << std::endl;
          GetJacobiStates( twoT, twoJ12, parity12_bra, E12_bra, ibraNAS, jac1_bra, jac2_bra);
//          std::cout << "checking isospin bra: " << jac1_bra.t << " " << Tab << std::endl;
          if ( jac1_bra.t != Tab ) continue;
          double Tcoeff_bra = Tcoeff_wrapper( bra, Jab, twoJ, jac1_bra, jac2_bra, twoJ12, Ncm, Lcm);
          if (std::abs(Tcoeff_bra)<1e-9) continue;

          for (int iketNAS=0; iketNAS<NASdim_ket; iketNAS++)
          {
 //          std::cout << "  ibraNAS,iketNAS = " << ibraNAS << " " << iketNAS << std::endl;
           GetJacobiStates( twoT, twoJ12, parity12_ket, E12_ket, iketNAS, jac1_ket, jac2_ket);
//           std::cout << "     checking isospin ket: " << jac1_ket.t << " " << Tde << std::endl;
           if ( jac1_ket.t != Tde ) continue;
           double Tcoeff_ket = Tcoeff_wrapper( ket, Jde, twoJ, jac1_ket, jac2_ket, twoJ12, Ncm, Lcm);

//           std::cout << "    ** finding the NAS matrix element. T = " << T << "   J = " << J12 << "  p = " << parity12_bra << std::endl;

//           double meNAS = ElementNAS( ibraNAS, iketNAS, E12_bra, E12_ket, twoT, twoJ12, parity12_bra );
           double meNAS = GetMatElNAS( ibraNAS, iketNAS, E12_bra, E12_ket, twoT, twoJ12, parity12_bra );

           me_lab +=  Tcoeff_bra * meNAS * Tcoeff_ket  ;
           std::cout << " " << ibraNAS << " , " << iketNAS << ":  (( " << Tcoeff_bra << "  " << meNAS << "  " << Tcoeff_ket << " )) " << std::endl;
           std::cout << "    < N1=" << jac1_bra.n << " L1=" << jac1_bra.l << " S1=" << jac1_bra.s << " J1=" << jac1_bra.j << " T1=" <<jac1_bra.t << ",  N2=" << jac2_bra.n << " L2=" << jac2_bra.l << " T=" << twoT << " J12=" << twoJ12 << " | ..." << std::endl;

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



