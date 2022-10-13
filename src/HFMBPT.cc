#include <iomanip>
#include "HFMBPT.hh"
#include "HartreeFock.hh"
#include "PhysicalConstants.hh"

#include <omp.h>

HFMBPT::~HFMBPT()
{}

HFMBPT::HFMBPT(Operator& hbare)
  : HartreeFock(hbare),
    C_HO2NAT(modelspace->GetNumberOrbits(),modelspace->GetNumberOrbits(),arma::fill::eye),
    C_HF2NAT(modelspace->GetNumberOrbits(),modelspace->GetNumberOrbits(),arma::fill::eye),
    use_NAT_occupations(false), NAT_order("occupation")
{}

//*********************************************************************
// post Hartree-Fock method
// This method assumes that we've already performed a Hartree-Fock calculation, so that the matrix C maps
// the HO basis to the HF basis. The first step is to compute the 1b density matrix rho up to the MBPT2 level,
// using the NO2B Hamiltonian in the HF basis. The second step is to diagonalize the density matrix.
// An optional third step is to change the occupations to match the eigenvalues of the density matrix.
// This requires some thinking, since if all orbits are fractionally occupied, then count as both hole and particle
// states, and routines that benefit from only acting on hh or ph blocks now act on all orbits, and things slow down.
// As a compromise, we place a threshold on the occupations, and all orbits with occupations below that threshold
// (which will be most of them) are set to zero and labeled as "particle". To maintain the right particle number,
// we "back-fill" the hole states using a somewhat arbitrary but not totally outlandish procedure.
//*********************************************************************
void HFMBPT::GetNaturalOrbitals()
{

  int norbits = HartreeFock::modelspace->GetNumberOrbits();
  int A = HartreeFock::modelspace->GetTargetMass();
  Occ      = arma::vec(norbits,arma::fill::zeros);  // Occupations are the eigenvalues of the density matrix
  GetDensityMatrix();  // Include 2nd order MBPT corrections to rho
  DiagonalizeRho();  // Find the 1b transformation that diagonalizes rho, but don't apply it to anything yet.

  double AfromTr = 0.0;
  for(int i=0; i< norbits; ++i)
  {
    Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
    AfromTr += rho(i,i) * (oi.j2+1);
  }

  if(std::abs(AfromTr - A) > 1e-8)
  {
    std::cout << "Warning: Mass != Tr(rho)   " << A << " != " << AfromTr <<  std::endl;
    exit(0);
  }
  C_HO2NAT = C * C_HF2NAT;

  // set the occ_nat values
  for ( auto i : modelspace->all_orbits)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    oi.occ_nat = std::abs(Occ(i));  // it's possible that Occ(i) is negative, and for occ_nat, we don't want that.
  }

  if (use_NAT_occupations) // use fractional occupation
  {

    double keep_occ_threshold = 0.02; // anything with occupation less than keep_occ_threshold gets set to zero, otherwise everything is a hole.

    double NfromTr=0;
    double ZfromTr=0;
    std::cout << "Switching to occupation numbers obtained from 2nd order 1b density matrix." << std::endl;
    std::vector<index_t> holeorbs_tmp;
    std::vector<double> hole_occ_tmp;
    // Figure out how many particles are living in orbits with occupations above our threshold.
    // We do this separately for protons and neutrons.
    for(auto& i : modelspace->all_orbits)
    {
      Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
      if (Occ(i) > keep_occ_threshold)
      {
        holeorbs_tmp.push_back(i);
        hole_occ_tmp.push_back(Occ(i));
        NfromTr += (1+oi.tz2)/2 * Occ(i) * (oi.j2+1);
        ZfromTr += (1-oi.tz2)/2 * Occ(i) * (oi.j2+1);
      }
    }

    // Now do the back-filling.
    int Z = modelspace->GetZref();
    int N = A-Z;
    double aloquot = 0.005; // This is how much we increase an occupation in each back-filling pass.

   // Back-fill the protons
    while ( (Z-ZfromTr) > ModelSpace::OCC_CUT )
    {
      for (size_t i=0; i<holeorbs_tmp.size(); i++)
      {
        Orbit& oi = HartreeFock::modelspace->GetOrbit(holeorbs_tmp[i]);
        if ( oi.tz2 >0 ) continue;
        double occ_increase = std::min( { aloquot, (Z-ZfromTr)/(oi.j2+1.), 1.0-hole_occ_tmp[i] } );
        hole_occ_tmp[i] += occ_increase;
        ZfromTr += occ_increase * (oi.j2+1);
        if ( (Z-ZfromTr) < ModelSpace::OCC_CUT ) break;
      }
    }

   // Back-fill the neutrons
    while ( (N-NfromTr) > ModelSpace::OCC_CUT )
    {
      for (size_t i=0; i<holeorbs_tmp.size(); i++)
      {
        Orbit& oi = HartreeFock::modelspace->GetOrbit(holeorbs_tmp[i]);
        if ( oi.tz2 <0 ) continue;
        double occ_increase = std::min( { aloquot, (N-NfromTr)/(oi.j2+1.), 1.0-hole_occ_tmp[i] } );
        hole_occ_tmp[i] += occ_increase;
        NfromTr += occ_increase * (oi.j2+1);
        if ( (N-NfromTr) < ModelSpace::OCC_CUT ) break;
      }
    }

    holeorbs = arma::uvec( holeorbs_tmp );
    hole_occ = arma::rowvec( hole_occ_tmp );

    // Now we tell modelspace about the new occupations, and any needed reclassification of 'holes' and 'particles'
    UpdateReference();

   // Calling UpdateReference screws up the occ_nat values, so we need to set them again here.
    for ( auto i : HartreeFock::modelspace->all_orbits )
    {
      auto& oi = HartreeFock::modelspace->GetOrbit(i);
      oi.occ_nat = std::abs(Occ(i));  // it's possible that Occ(i) is negative, and for occ_nat, we don't want that.
    }

  } // if use_NAT_occupations

  // In principle, this could be iterated until self-consistency. But I don't do that for now. This would only have an impact
  // if the energy ordering of occupied and unoccupied levels gets flipped, which would be a fairly pathological case.
  arma::mat tmp = C_HO2NAT.cols(holeorbs);
  rho = (tmp.each_row() % hole_occ) * tmp.t(); // now rho is in the HO basis, with our prescribed occupations in the NAT basis
  UpdateF();  // Now F is in the HO basis, but with rho from filling in the NAT basis.

  ReorderHFMBPTCoefficients();
  C_HO2NAT = C * C_HF2NAT; // bug fix suggested by Emily Love
}

//*********************************************************************
// Diagonalize the 1b density matrix
//*********************************************************************
void HFMBPT::DiagonalizeRho()
{
  for (auto& it : Hbare.OneBodyChannels)
  {
//    arma::uvec orbvec(it.second);
    arma::uvec orbvec(std::vector<index_t>(it.second.begin(),it.second.end()));
//    arma::uvec orbvec_d = arma::sort(orbvec, "descend");
//    arma::uvec orbvec_d = sort(orbvec, "descend");
    arma::mat rho_ch = rho.submat(orbvec, orbvec);
    arma::mat vec;
    arma::vec eig;
    bool success = false;
    success = arma::eig_sym(eig, vec, rho_ch); // eigenvalues of rho (i.e. occupations) are in ascending order
    if(not success)
    {
      std::cout << "Error in diagonalization of density matrix" << std::endl;
      std::cout << "Density Matrix:" << std::endl;
      rho_ch.print();
      exit(0);
    }
//    Occ(orbvec_d) = eig; // assigning with descending indexes produces NAT orbits ordered by descending occupation.
//    C_HF2NAT.submat(orbvec, orbvec_d) = vec; // NAT orbits ordered by descending occupation

    Occ(orbvec) = arma::reverse(eig); // reverse so NAT orbits ordered by descending occupation.
    C_HF2NAT.submat(orbvec, orbvec) = arma::reverse(vec, 1); // "1" means reverse elements in each row
//    std::cout << "1-body channel ";
//    for ( auto x : it.first ) std::cout << x << " ";
//    std::cout << "  : " << std::endl << Occ(orbvec) << std::endl;
//    std::cout << " C submat: " << std::endl << C_HF2NAT.submat(orbvec, orbvec ) << std::endl;
  }
 // Choose ordering and phases so that C_HF2NAT looks as close to the identity as possible
 // NB: Rather than C_HF2NAT being close to the identity, we really want the orbits for a given (l,j,tz) ordered
 // according to decreasing occupation, so that a later emax_imsrg truncation is reasonable. -SRS Jan 2021.
//  ReorderHFMBPTCoefficients();
//  std::cout << " line " << __LINE__ << "   Occ = " << std::endl << Occ << std::endl;
}

//*********************************************************************
// Transform some operator from the HF basis to the NAT basis
//*********************************************************************
Operator HFMBPT::TransformHFToNATBasis( Operator& OpHF)
{
  Operator OpNAT(OpHF);
  OpNAT.OneBody = C_HF2NAT.t() * OpHF.OneBody * C_HF2NAT;

  for (auto& it : OpHF.TwoBody.MatEl )
  {
    int ch_bra = it.first[0];
    int ch_ket = it.first[1];
    TwoBodyChannel& tbc_bra = OpNAT.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = OpNAT.modelspace->GetTwoBodyChannel(ch_ket);
    int nbras = it.second.n_rows;
    int nkets = it.second.n_cols;
    arma::mat Dbra(nbras,nbras);
    arma::mat Dket(nkets,nkets);

    for (int i = 0; i<nkets; ++i)
    {
      Ket & ket_hf = tbc_ket.GetKet(i);
      for (int j=0; j<nkets; ++j)
      {
        Ket & ket_nat = tbc_ket.GetKet(j);
        Dket(i,j) = C_HF2NAT(ket_hf.p,ket_nat.p) * C_HF2NAT(ket_hf.q,ket_nat.q);
        if(ket_hf.p != ket_hf.q)
        {
          Dket(i,j) += C_HF2NAT(ket_hf.q, ket_nat.p) * C_HF2NAT(ket_hf.p, ket_nat.q) *
            ket_hf.Phase(tbc_ket.J);
        }
        if (ket_hf.p==ket_hf.q)    Dket(i,j) *= PhysConst::SQRT2;
        if (ket_nat.p==ket_nat.q)  Dket(i,j) /= PhysConst::SQRT2;
      }
    }
    if (ch_bra == ch_ket) {
      Dbra = Dket.t();
    }
    else
    {
      for (int i=0; i<nbras; ++i)
      {
        Ket & bra_nat = tbc_bra.GetKet(i);
        for (int j=0; j<nbras; ++j)
        {
          Ket & bra_hf = tbc_bra.GetKet(j);
          Dbra(i,j) = C_HF2NAT(bra_hf.p,bra_nat.p) * C_HF2NAT(bra_hf.q,bra_nat.q);
          if (bra_hf.p!=bra_hf.q)
          {
            Dbra(i,j) += C_HF2NAT(bra_hf.q, bra_nat.p) * C_HF2NAT(bra_hf.p, bra_nat.q)
              * bra_hf.Phase(tbc_bra.J);
          }
          if (bra_hf.p==bra_hf.q)    Dbra(i,j) *= PhysConst::SQRT2;
          if (bra_nat.p==bra_nat.q)  Dbra(i,j) /= PhysConst::SQRT2;
        }
      }
    }
    auto& IN  =  it.second;
    auto& OUT =  OpNAT.TwoBody.GetMatrix(ch_bra,ch_ket);
    OUT  =    Dbra * IN * Dket;
  }
  return OpNAT;
}

//*********************************************************************
// Transfrom some operator from the HO basis to the NAT basis
//*********************************************************************
Operator HFMBPT::TransformHOToNATBasis( Operator& OpHO)
{
  Operator OpNAT(OpHO);
  OpNAT.OneBody = C_HO2NAT.t() * OpHO.OneBody * C_HO2NAT;

  for (auto& it : OpHO.TwoBody.MatEl )
  {
    int ch_bra = it.first[0];
    int ch_ket = it.first[1];
    TwoBodyChannel& tbc_bra = OpNAT.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = OpNAT.modelspace->GetTwoBodyChannel(ch_ket);
    int nbras = it.second.n_rows;
    int nkets = it.second.n_cols;
    arma::mat Dbra(nbras,nbras);
    arma::mat Dket(nkets,nkets);

    for (int i = 0; i<nkets; ++i)
    {
      Ket & ket_ho = tbc_ket.GetKet(i);
      for (int j=0; j<nkets; ++j)
      {
        Ket & ket_nat = tbc_ket.GetKet(j);
        Dket(i,j) = C_HO2NAT(ket_ho.p,ket_nat.p) * C_HO2NAT(ket_ho.q,ket_nat.q);
        if(ket_ho.p != ket_ho.q)
        {
          Dket(i,j) += C_HO2NAT(ket_ho.q, ket_nat.p) * C_HO2NAT(ket_ho.p, ket_nat.q) *
            ket_ho.Phase(tbc_ket.J);

        }
        if (ket_ho.p==ket_ho.q)    Dket(i,j) *= PhysConst::SQRT2;
        if (ket_nat.p==ket_nat.q)  Dket(i,j) /= PhysConst::SQRT2;
      }
    }
    if (ch_bra == ch_ket) {
      Dbra = Dket.t();
    }
    else
    {
      for (int i=0; i<nbras; ++i)
      {
        Ket & bra_nat = tbc_bra.GetKet(i);
        for (int j=0; j<nbras; ++j)
        {
          Ket & bra_ho = tbc_bra.GetKet(j);
          Dbra(i,j) = C_HO2NAT(bra_ho.p,bra_nat.p) * C_HO2NAT(bra_ho.q,bra_nat.q);
          if (bra_ho.p!=bra_ho.q)
          {
            Dbra(i,j) += C_HO2NAT(bra_ho.q,bra_nat.p) * C_HO2NAT(bra_ho.p,bra_nat.q) *
              bra_ho.Phase(tbc_bra.J);
          }
          if (bra_ho.p==bra_ho.q)    Dbra(i,j) *= PhysConst::SQRT2;
          if (bra_nat.p==bra_nat.q)  Dbra(i,j) /= PhysConst::SQRT2;
        }
      }
    }

    auto& IN  =  it.second;
    auto& OUT =  OpNAT.TwoBody.GetMatrix(ch_bra,ch_ket);
    OUT  =    Dbra * IN * Dket;
   }
   return OpNAT;
}

//*********************************************************************
// Get the normal ordered Hamiltonian in the NAT basis (with residual 3N
// discarded). The methods UpdateF() and CalcEHF() use the density matrix
// rho which, if we're using the naive-filling occupations, is different
// from the 1b density matrix obtained from MBPT2. So we temporarily
// store rho in rho_swap and make rho correspond to the desired occupations
// in the NAT basis.
//*********************************************************************
//Operator HFMBPT::GetNormalOrderedHNAT()
Operator HFMBPT::GetNormalOrderedHNAT(int particle_rank)
{
  double start_time = omp_get_wtime();
  std::cout << "Getting normal-ordered H in NAT basis" << std::endl;
//  arma::mat rho_swap = rho;
//  arma::mat tmp = C_HO2NAT.cols(holeorbs);
//  rho = (tmp.each_row() % hole_occ) * tmp.t();

//  UpdateF();  // Now F is in the HO basis, but with rho from filling in the NAT basis.
  CalcEHF();
  std::cout << std::fixed <<  std::setprecision(7);
  std::cout << "e1Nat = " << e1hf << std::endl;
  std::cout << "e2Nat = " << e2hf << std::endl;
  std::cout << "e3Nat = " << e3hf << std::endl;
  std::cout << "E_Nat = "  << EHF  << std::endl;

//  Operator HNO = Operator(*HartreeFock::modelspace,0,0,0,2);
  Operator HNO = Operator(*HartreeFock::modelspace,0,0,0,particle_rank);
  HNO.ZeroBody = EHF;
  HNO.OneBody = C_HO2NAT.t() * F * C_HO2NAT;

  int nchan = HartreeFock::modelspace->GetNumberTwoBodyChannels();



    for (int ch=0;ch<nchan;++ch)
    {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int npq = tbc.GetNumberKets();

      arma::mat D(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
      arma::mat V3NO(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
#pragma omp parallel for schedule(dynamic,1)
      for (int i=0; i<npq; ++i)
      {
        Ket & bra = tbc.GetKet(i);
        int e2bra = 2*bra.op->n + bra.op->l + 2*bra.oq->n + bra.oq->l;
        for (int j=0; j<npq; ++j)
        {
          Ket & ket = tbc.GetKet(j);
          int e2ket = 2*ket.op->n + ket.op->l + 2*ket.oq->n + ket.oq->l;
          D(i,j) = C_HO2NAT(bra.p,ket.p) * C_HO2NAT(bra.q,ket.q);
          if (bra.p!=bra.q)
          {
            D(i,j) += C_HO2NAT(bra.q,ket.p) * C_HO2NAT(bra.p,ket.q) * bra.Phase(J);
          }
          if (bra.p==bra.q)    D(i,j) *= PhysConst::SQRT2;
          if (ket.p==ket.q)    D(i,j) /= PhysConst::SQRT2;

          // Now generate the NO2B part of the 3N interaction
          if (Hbare.GetParticleRank()<3) continue;
          if (i>j) continue;
          for ( auto a : modelspace->all_orbits )
          {
            Orbit & oa = modelspace->GetOrbit(a);
            if ( 2*oa.n+oa.l+e2bra > Hbare.GetE3max() ) continue;
            for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
            {
              if ( std::abs(rho(a,b)) < 1e-8 ) continue; // Turns out this helps a bit (factor of 5 speed up in tests)
              Orbit & ob = HartreeFock::modelspace->GetOrbit(b);
              if ( 2*ob.n+ob.l+e2ket > Hbare.GetE3max() ) continue;
//              V3NO(i,j) += rho(a,b) * GetVNO2B(bra.p, bra.q, a, ket.p, ket.q, b, J);
              V3NO(i,j) += rho(a,b) * Hbare.ThreeBody.GetME_pn_no2b(bra.p,bra.q,a,ket.p,ket.q,b,J);
            }
          }
          V3NO(i,j) /= (2*J+1);
          if (bra.p==bra.q)  V3NO(i,j) /= PhysConst::SQRT2;
          if (ket.p==ket.q)  V3NO(i,j) /= PhysConst::SQRT2;
          V3NO(j,i) = V3NO(i,j);
        }
      }

      auto& V2  =  Hbare.TwoBody.GetMatrix(ch);
      auto& OUT =  HNO.TwoBody.GetMatrix(ch);
      OUT  =    D.t() * (V2 + V3NO) * D;
    }


  if (particle_rank>2)
  {
//    HNO.ThreeBody = GetTransformed3B();
//    HNO.ThreeBody = GetTransformed3B( Hbare );
    HNO.ThreeBody = GetTransformed3B( Hbare, C_HO2NAT );
  }

//  rho = rho_swap;
  profiler.timer["HFMBT_GetNormalOrderedHNO"] += omp_get_wtime() - start_time;
  return HNO;
}

//*********************************************************************
// Compute the MBPT2 corrections to the 1b density matrix.
//*********************************************************************
void HFMBPT::GetDensityMatrix()
{
  Operator Hhf = HartreeFock::GetNormalOrderedH();
  Operator& H(Hhf);
  double t_start = omp_get_wtime();

  // After the HF step, rho is the density of harmonic oscillator states for a filled HF reference
  //  i.e. <a|rho|b> = sum_i <a|i> <b|i>  where i is an occupied HF state and a and b are HO basis states.
  // Now we switch to the HF basis, so rho should be diagonal before adding in the perturbative corrections.
  rho.zeros(); // This and the following line fixes bug found by Baishan Dec 2020.
  for (auto& i : HartreeFock::modelspace->holes)  rho(i,i) = HartreeFock::modelspace->GetOrbit(i).occ; // Set hole occupations to 1.

//  std::cout << std::endl << "before perturbative correction, rho is" << std::endl << rho << std::endl;
  // compute second order corrections to the density matrix
  DensityMatrixPP(H);
  DensityMatrixHH(H);
  DensityMatrixPH(H);
//  std::cout << std::endl << "after perturbative correction, rho is" << std::endl << rho << std::endl;
  profiler.timer["HFMBPT DensityMatrix"] += omp_get_wtime() - t_start;
}

//*********************************************************************
// Pretty self explanatory. Print the quantum numbers and occupations
// of all orbits, except that the occupation is the eigenvalue of the
// density matrix, not the value set in modelspace.
//*********************************************************************
void HFMBPT::PrintOccupation()
{

//  for (int i=0; i<HartreeFock::modelspace->GetNumberOrbits(); ++i)
  for (auto& i : modelspace->all_orbits)
  {
    Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
    std::cout << std::fixed << std::setw(4) << oi.n << std::setw(4) << oi.l <<
      std::setw(4) << oi.j2 << std::setw(4) << oi.tz2 << "   " <<
      std::setw(8) << Occ(i) << std::endl;
  }
}

//*********************************************************************
// Compute the MBPT2 contribution to rho due to < 1| rho |1 > where
//  |1> is the 1st order correction to the HF ground state.
// Here we treat the contribution to the particle-particle block.
//
//        *~~~~~~~*      <a|rho|b> = 1/2 sum_{cijJ} (2J+1)/(2j_a+1) <ac|V|ij><ij|V|bc> / Delta
// R0--- / \     / \a
//     c(  i)  j(  (RHO)      with Delta = (ea+ec -ei-ej)(eb+ec-ei-ej)
// R0--- \ /     \ /b
//        *~~~~~~~*        R0---  indicates the MBPT resolvent lines.
//                         ijk are holes, abc are particles
//
// In the case that a hole and particle level are very closely spaced in energy,
// the denominator can get small, and perturbation theory breaks down.
// In this case, the second order correction can produce a diagonal
// element larger than 1 (or less than zero for the hole orbits).
// To deal with this, we replace the V*V/Delta by the result of two
// level mixing
//   V*V / (e_acij*e_bcij) -> 1/2 * [sqrt( V*V + E*E) - E] / sqrt( V*V + E*E)
// where E*E = 1/4 (e_acij * e_bcij)
// In the limit V*V << E*E, this coincides with the perturbative expression,
// so we can make the replacement even if the levels aren't closely spaced.
//
// When using fraction occupations, including the factor (1-na)(1-nb) [ni nj(1-nc)]^2
// and an equivalent one in the HH term produces a density matrix with the correct
// particle number, encoded in the 2j+1 weighted trace of rho.
//*********************************************************************
void HFMBPT::DensityMatrixPP(Operator& H)
{
  arma::mat rho_pp2 = arma::zeros( arma::size(rho) );
//  for (auto& a : HartreeFock::modelspace->all_orbits)
  for (auto& a : HartreeFock::modelspace->particles)
  {
    double ea = H.OneBody(a,a);
    Orbit& oa = HartreeFock::modelspace->GetOrbit(a);
    if ( (1-oa.occ)<ModelSpace::OCC_CUT) continue;

    for (auto& b : modelspace->OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
    {
      if(b > a) continue;
      double eb = H.OneBody(b,b);
      Orbit& ob = HartreeFock::modelspace->GetOrbit(b);
      if((1-ob.occ) <ModelSpace::OCC_CUT) continue;

      double r = 0.0;
      for(auto& c : HartreeFock::modelspace->particles)
//      for(auto& c : HartreeFock::modelspace->all_orbits)
      {
        double ec = H.OneBody(c,c);
        Orbit& oc = HartreeFock::modelspace->GetOrbit(c);
        if ( (1-oc.occ)<ModelSpace::OCC_CUT) continue;

        for(auto& i : HartreeFock::modelspace->holes)
        {
          double ei = H.OneBody(i,i);
          Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
          if ( oi.occ < ModelSpace::OCC_CUT ) continue;

          for(auto& j : HartreeFock::modelspace->holes){
            double ej = H.OneBody(j,j);
            Orbit& oj = HartreeFock::modelspace->GetOrbit(j);
            if ( oj.occ < ModelSpace::OCC_CUT ) continue;

            double e_acij = ea + ec - ei - ej;
            double e_bcij = eb + ec - ei - ej;
            if(std::abs(e_acij*e_bcij) < 1.e-8) continue;
            int Jmin = std::max(std::abs(oa.j2-oc.j2), std::max(std::abs(oi.j2-oj.j2), std::abs(ob.j2-oc.j2)))/2;
            int Jmax = std::min(         oa.j2+oc.j2,  std::min(         oi.j2+oj.j2,           ob.j2+oc.j2))/2;

            double tbme = 0.0;
            for(int J = Jmin; J <= Jmax; ++J){
              tbme +=  (2*J+1) * H.TwoBody.GetTBME_J(J,a,c,i,j)
                               * H.TwoBody.GetTBME_J(J,i,j,b,c);
            }
            tbme *=  (1-oa.occ) * (1-ob.occ) *  (1-oc.occ)*(1-oc.occ) * oi.occ*oi.occ * oj.occ*oj.occ ;
            double epsilon = 0.5*sqrt(std::abs(e_acij * e_bcij));
            r += 0.5* ( sqrt(tbme + epsilon*epsilon) - epsilon ) / sqrt(tbme + epsilon*epsilon);
          }
        }
      }
      rho_pp2(a,b) = r * 0.5 / (oa.j2+1);
      rho_pp2(b,a) = r * 0.5 / (oa.j2+1);
    }
  }
  rho += rho_pp2;
}


//*********************************************************************
// Compute the MBPT2 contribution to rho due to < 1| rho |1 > where
//  |1> is the 1st order correction to the HF ground state.
// Here we treat the contribution to the hole-hole block.
//
//        *~~~~~~~*      <i|rho|j> = -1/2 sum_{abkJ} (2J+1)/(2j_i+1) <ab|V|ik><jk|V|ab> / Delta
// R0--- / \     / \j
//     a(  k)  b(  (RHO)      with Delta = (ea+eb -ei-ek)(ea+eb-ej-ek)
// R0--- \ /     \ /i
//        *~~~~~~~*        R0---  indicates the MBPT resolvent lines.
//                         ijk are holes, abc are particles
//
// See discussion above HFMBPT::DensityMatrixPP for more details.
//*********************************************************************
void HFMBPT::DensityMatrixHH(Operator& H)
{
  arma::mat rho_hh2 = arma::zeros( arma::size(rho) );
  for (auto& i : HartreeFock::modelspace->holes)
  {
    double ei = H.OneBody(i,i);
    Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
    if ( oi.occ < ModelSpace::OCC_CUT) continue;

    for (auto& j : HartreeFock::modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}))
    {
      if(j > i) continue;
      double ej = H.OneBody(j,j);
      Orbit& oj = HartreeFock::modelspace->GetOrbit(j);
      if ( oj.occ<ModelSpace::OCC_CUT) continue;

      double r = 0.0;
      for(auto& a : HartreeFock::modelspace->particles)
//      for(auto& a : HartreeFock::modelspace->all_orbits)
      {
        double ea = H.OneBody(a,a);
        Orbit& oa = HartreeFock::modelspace->GetOrbit(a);
        if ( (1-oa.occ)<ModelSpace::OCC_CUT) continue;

        for(auto& b : HartreeFock::modelspace->particles)
//        for(auto& b : HartreeFock::modelspace->all_orbits)
        {
          double eb = H.OneBody(b,b);
          Orbit& ob = HartreeFock::modelspace->GetOrbit(b);
          if ( (1-ob.occ)<ModelSpace::OCC_CUT) continue;

          for(auto& k : HartreeFock::modelspace->holes)
          {
            double ek = H.OneBody(k,k);
            Orbit& ok = HartreeFock::modelspace->GetOrbit(k);
            if ( ok.occ<ModelSpace::OCC_CUT) continue;

            double e_abik = ea + eb - ei - ek;
            double e_abjk = ea + eb - ek - ej;
            if( std::abs(e_abik*e_abjk) < 1.e-8) continue;
            int Jmin = std::max(std::abs(oa.j2-ob.j2), std::max(std::abs(oi.j2-ok.j2), std::abs(oj.j2-ok.j2)))/2;
            int Jmax = std::min(         oa.j2+ob.j2,  std::min(         oi.j2+ok.j2,           oj.j2+ok.j2))/2;

            double tbme = 0.0;
            for(int J = Jmin; J <= Jmax; ++J)
            {
              tbme += (2*J+1) * H.TwoBody.GetTBME_J(J,a,b,i,k)
                              * H.TwoBody.GetTBME_J(J,j,k,a,b);
            }

            tbme *=  (1-oa.occ)*(1-oa.occ) * (1-ob.occ)*(1-ob.occ) * ok.occ*ok.occ   *  oi.occ * oj.occ ;
            if (true)
            {
              double epsilon = 0.5*sqrt(std::abs(e_abik * e_abjk));
              r += 0.5* ( sqrt(tbme + epsilon*epsilon) - epsilon ) / sqrt(tbme + epsilon*epsilon) ;
            }
            else // the MBPT expression. We don't actually use this.
            {
              r += tbme / (e_abik * e_abjk);
            }
          }
        }
      }
      rho_hh2(i,j) = - r * 0.5 / (oi.j2+1);
      rho_hh2(j,i) = - r * 0.5 / (oi.j2+1);
    }
  }
  rho += rho_hh2;
}

//*********************************************************************
// Compute the MBPT2 contribution to rho due to <0|rho|2> + <2|rho|0>
//  where |0> is the HF ground state and |2> is the 2nd order correction.
//
//      (RHO)           <i|rho|a> = 1/2 sum_{bcjJ} (2J+1)/(2j_i+1) <aj|V|bc><bc|V|ij> / Delta
// R0--- / \a
//     i(   )~~~~*       with Delta = (ea-ei)(eb+ec-ei-ej)
// R0--- \ /b  j( )c
//        *~~~~~~*        R0---  indicates the MBPT resolvent lines.
//                         ijk are holes, abc are particles
// and
//

//      (RHO)           <i|rho|a> = -1/2 sum_{bkjJ} (2J+1)/(2j_i+1) <kj|V|ib><ab|V|kj> / Delta
// R0--- / \i
//     a(   )~~~~*       with Delta = (ea-ei)(ea+eb-ej-ek)
// R0--- \ /j  b( )k
//        *~~~~~~*
//
// Equivalent diagrams can be drawn with rho on the bottom, corresponding to <2|rho|0>,
// and the formulas are the same.
// In (limited) tests, a small gap between particle and hole levels did not appear
// to be a problem for these diagrams, so the MBPT2 expression is used directly.
// This may need to be revisited.
//*********************************************************************
void HFMBPT::DensityMatrixPH(Operator& H)
{

  arma::mat rho_ph2 = arma::zeros( arma::size(rho) );
  for (auto& i : HartreeFock::modelspace->holes)
  {
    double ei = H.OneBody(i,i);
    Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
    if ( oi.occ < ModelSpace::OCC_CUT) continue;

    for (auto& a : HartreeFock::modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}))
    {
      double ea = H.OneBody(a,a);
      Orbit& oa = HartreeFock::modelspace->GetOrbit(a);
      if ( (1-oa.occ)<ModelSpace::OCC_CUT) continue;

      double r = 0.0;
      for(auto& b : HartreeFock::modelspace->particles)
      {
        double eb = H.OneBody(b,b);
        Orbit& ob = HartreeFock::modelspace->GetOrbit(b);
        if ( (1-ob.occ)<ModelSpace::OCC_CUT) continue;

        for(auto& c : HartreeFock::modelspace->particles)
        {
          double ec = H.OneBody(c,c);
          Orbit& oc = HartreeFock::modelspace->GetOrbit(c);
          if ( (1-oc.occ)<ModelSpace::OCC_CUT) continue;

          for(auto& j : HartreeFock::modelspace->holes)
          {
            double ej = H.OneBody(j,j);
            Orbit& oj = HartreeFock::modelspace->GetOrbit(j);
            if ( oj.occ < ModelSpace::OCC_CUT) continue;

            double e_ai = ea - ei;
            double e_bcij = eb + ec - ei - ej;
            if(e_ai*e_bcij < 1.e-8) continue;
            int Jmin = std::max(std::abs(oa.j2-oj.j2), std::max(std::abs(ob.j2-oc.j2), std::abs(oi.j2-oj.j2)))/2;
            int Jmax = std::min(         oa.j2+oj.j2,  std::min(         ob.j2+oc.j2,           oi.j2+oj.j2))/2;

            double tbme = 0.0;
            for(int J = Jmin; J <= Jmax; ++J)
            {
              tbme += (2*J+1) * H.TwoBody.GetTBME_J(J,a,j,b,c)
                              * H.TwoBody.GetTBME_J(J,b,c,i,j);
            }

            tbme *= (1-oa.occ) * (1-ob.occ) * (1-oc.occ) * oi.occ * oj.occ ;
            r += tbme / (e_ai * e_bcij);
          }
        }
      }
      rho_ph2(a,i) += r * 0.5 / (oa.j2+1);
      rho_ph2(i,a) += r * 0.5 / (oa.j2+1);
    }
  }



  for (auto& i : HartreeFock::modelspace->holes)
  {
      double ei = H.OneBody(i,i);
      Orbit& oi = HartreeFock::modelspace->GetOrbit(i);

    for (auto& a : HartreeFock::modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}))
    {
      double ea = H.OneBody(a,a);
      Orbit& oa = HartreeFock::modelspace->GetOrbit(a);
      if ( (1-oa.occ)<ModelSpace::OCC_CUT) continue;


      double r = 0.0;
      for(auto& b : HartreeFock::modelspace->particles)
      {
        double eb = H.OneBody(b,b);
        Orbit& ob = HartreeFock::modelspace->GetOrbit(b);
        if ( (1-ob.occ)<ModelSpace::OCC_CUT) continue;

        for(auto& j : HartreeFock::modelspace->holes)
        {
          double ej = H.OneBody(j,j);
          Orbit& oj = HartreeFock::modelspace->GetOrbit(j);

          for(auto& k : HartreeFock::modelspace->holes)
          {
            double ek = H.OneBody(k,k);
            Orbit& ok = HartreeFock::modelspace->GetOrbit(k);

            double e_ai = ea - ei;
            double e_abkj = ea + eb - ek - ej;
            if(e_ai*e_abkj < 1.e-8) continue;
            int Jmin = std::max(std::abs(ok.j2-oj.j2), std::max(std::abs(oi.j2-ob.j2), std::abs(oa.j2-ob.j2)))/2;
            int Jmax = std::min(         ok.j2+oj.j2,  std::min(         oi.j2+ob.j2,           oa.j2+ob.j2))/2;

            double tbme = 0.0;
            for(int J = Jmin; J <= Jmax; ++J)
            {
              tbme += (2*J+1) * H.TwoBody.GetTBME_J(J,k,j,i,b)
                              * H.TwoBody.GetTBME_J(J,a,b,k,j);
            }
            tbme *= (1-oa.occ) * oi.occ * oj.occ * ok.occ * (1-ob.occ);
            r += tbme / (e_ai * e_abkj);
          }
        }
      }
      rho_ph2(a,i) -= r * 0.5 / (oa.j2+1);
      rho_ph2(i,a) -= r * 0.5 / (oa.j2+1);
    }
  }
  rho += rho_ph2;
}




//*********************************************************************
// Specialization of the HartreeFock version. This is because we want
// to print out the wave functions in terms of the HO components,
// and so we need to use C_HO2NAT rather than C.
// If no transformation to the NAT basis has been performed, then C_HF2NAT
// is just the identity and we get the Hartree-Fock wave functions.
//*********************************************************************
void HFMBPT::PrintSPEandWF()
{
  C_HO2NAT = C * C_HF2NAT;
  arma::mat F_natbasis = C_HO2NAT.t() * F * C_HO2NAT;


   
  std::cout << std::fixed << std::setw(3) << "i" << ": " << std::setw(3) << "n" << " " << std::setw(3) << "l" << " "
       << std::setw(3) << "2j" << " " << std::setw(3) << "2tz" << "   " << std::setw(12) << "SPE" << " " << std::setw(12) << "occ."
       << " " << std::setw(12) << "occNAT" << "   |   " << " overlaps" << std::endl;
  for ( auto i : modelspace->all_orbits )
  {
    Orbit& oi = modelspace->GetOrbit(i);
    std::cout << std::fixed << std::setw(3) << i << ": " << std::setw(3) << oi.n << " " << std::setw(3) << oi.l << " "
         << std::setw(3) << oi.j2 << " " << std::setw(3) << oi.tz2 << "   " << std::setw(12) << std::setprecision(6) << F_natbasis(i,i) << " " << std::setw(12) << oi.occ << " " << std::setw(12) << oi.occ_nat << "   | ";
    for (int j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) // j runs over HO states
    {
      std::cout << std::setw(9) << C_HO2NAT(j,i) << "  ";  // C is <HO|NAT>
    }
    std::cout << std::endl;
  }
}




//*********************************************************************
// The reordering business probably isn't necessary because it should
// be taken care of in the DiagonalizeRho() step. However, there may
// be some unsightly minus signs. While they have no impact on any observables,
// it will occasionally make our lives easier to get rid of those.
//*********************************************************************
void HFMBPT::ReorderHFMBPTCoefficients()
{
   for (index_t i=0;i<C_HF2NAT.n_rows;++i) // loop through original basis states
   {
      if (C_HF2NAT(i,i) < 0)  C_HF2NAT.col(i) *= -1;
   }

  // would enums and switch/case be prettier here?
  if (NAT_order == "energy")
  {
    // note that we do this in just one pass, and that in the unlikely case that one of the occupied orbits
    // gets swapped with an unoccupied orbit, we would in principle need to recompute the energies and orderings
    // but let's just hope that this doesn't happen.
    std::cout << "Ordering NAT orbits according to increasing energy..." << std::endl;
    for (auto& it : Hbare.OneBodyChannels)
    {
      arma::uvec orbvec(std::vector<index_t>(it.second.begin(),it.second.end()));
      arma::mat CNAT_chan = C_HF2NAT.submat(orbvec, orbvec);
      arma::mat CHF_chan = C.submat(orbvec, orbvec);
      arma::mat fHO_chan = F.submat(orbvec,orbvec);
      arma::mat fNAT_chan = CNAT_chan.t() * CHF_chan.t() * fHO_chan * CHF_chan * CNAT_chan;
      arma::vec spe_NAT = fNAT_chan.diag();
      arma::uvec sorted_indices = arma::sort_index( spe_NAT, "ascend");
      arma::uvec orbvec_sorted = orbvec(sorted_indices);
      C_HF2NAT.submat(orbvec, orbvec) = C_HF2NAT( orbvec, orbvec_sorted); // sort the column indices, <row|col> = <HF|NAT>
      Occ(orbvec) = Occ(orbvec_sorted); 
    }
    for ( auto i : HartreeFock::modelspace->all_orbits )
    {
      auto& oi = HartreeFock::modelspace->GetOrbit(i);
      oi.occ_nat = std::abs(Occ(i));  // it's possible that Occ(i) is negative, and for occ_nat, we don't want that.
    }
  }
  else if (NAT_order == "mp2")
  {
    std::cout << "Ordering NAT orbits according to second order energy impact..." << std::endl;
    Operator H_temp = GetNormalOrderedHNAT(2); // particle rank 2 - we don't care about threebody here
//    arma::vec impacts = H_temp.GetMP2_Impacts();
    arma::vec impacts = GetMP2_Impacts(H_temp);

    for (auto& it : Hbare.OneBodyChannels)
    {
      arma::uvec orbvec(std::vector<index_t>(it.second.begin(),it.second.end()));
      arma::vec impacts_chan = impacts(orbvec);
      arma::uvec sorted_indices = arma::sort_index( impacts_chan, "ascend");
      arma::uvec orbvec_sorted = orbvec(sorted_indices);
      C_HF2NAT.submat(orbvec, orbvec) = C_HF2NAT( orbvec, orbvec_sorted); // sort the column indices, <row|col> = <HF|NAT>
      Occ(orbvec) = Occ(orbvec_sorted); 
    }
    for ( auto i : HartreeFock::modelspace->all_orbits )
    {
      auto& oi = HartreeFock::modelspace->GetOrbit(i);
      oi.occ_nat = std::abs(Occ(i));  // it's possible that Occ(i) is negative, and for occ_nat, we don't want that.
    }
  }
}







// Get a single 3-body matrix element in the HartreeFock basis.
// This is the straightforward but inefficient way to do it.
//double HartreeFock::GetHF3bme( int Jab, int Jde, int J2,  size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
//double HFMBPT::GetTransformed3bme( int Jab, int Jde, int J2,  size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
double HFMBPT::GetTransformed3bme( Operator& OpIn, int Jab, int Jde, int J2,  size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  double V_nat = 0.;
  Orbit& oa = modelspace->GetOrbit(a);
  Orbit& ob = modelspace->GetOrbit(b);
  Orbit& oc = modelspace->GetOrbit(c);
  Orbit& od = modelspace->GetOrbit(d);
  Orbit& oe = modelspace->GetOrbit(e);
  Orbit& of = modelspace->GetOrbit(f);

  for (auto alpha : OpIn.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
  {
   if ( std::abs(C_HO2NAT(alpha,a)) < 1e-8 ) continue;
   for (auto beta : OpIn.OneBodyChannels.at({ob.l,ob.j2,ob.tz2}) )
   {
    if ( std::abs(C_HO2NAT(beta,b)) < 1e-8 ) continue;
    for (auto gamma : OpIn.OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) )
    {
     if ( std::abs(C_HO2NAT(gamma,c)) < 1e-8 ) continue;
     for (auto delta : OpIn.OneBodyChannels.at({od.l,od.j2,od.tz2}) )
     {
      if ( std::abs(C_HO2NAT(delta,d)) < 1e-8 ) continue;
      for (auto epsilon : OpIn.OneBodyChannels.at({oe.l,oe.j2,oe.tz2}) )
      {
       if ( std::abs(C_HO2NAT(epsilon,e)) < 1e-8 ) continue;
       for (auto phi : OpIn.OneBodyChannels.at({of.l,of.j2,of.tz2}) )
       {
         double V_ho = OpIn.ThreeBody.GetME_pn( Jab,  Jde,  J2,  alpha,  beta,  gamma,  delta,  epsilon,  phi);
//         double V_ho = OpIn.ThreeBody.GetME( Jab,  Jde,  J2,  tab,  tde,  T2,  alpha,  beta,  gamma,  delta,  epsilon,  phi);
         V_nat += V_ho * C_HO2NAT(alpha,a) * C_HO2NAT(beta,b) * C_HO2NAT(gamma,c) * C_HO2NAT(delta,d) * C_HO2NAT(epsilon,e) * C_HO2NAT(phi,f);
       } // for phi
      } // for epsilon
     } // for delta
    } // for gamma
   } // for beta
  } // for alpha
  return V_nat;
}




// Modified version of GetMP2_Energy. Determines each orbital's impact on the total MP2 energy 
// (i.e., by how much would EMP2 change if this single orbital were removed)
arma::vec HFMBPT::GetMP2_Impacts(Operator& OpIn) const
{
//   std::cout << "  a    b    i    j    J    na     nb    tbme    denom     dE" << std::endl;
   double t_start = omp_get_wtime();
   double de = 0;

   int nparticles = OpIn.modelspace->particles.size();
   arma::vec orbit_impacts(OpIn.modelspace->all_orbits.size(), arma::fill::zeros);

   std::vector<index_t> particles_vec(OpIn.modelspace->particles.begin(),OpIn.modelspace->particles.end()); // convert set to vector for OMP looping
//   for ( auto& i : modelspace->particles)
//   #pragma omp parallel for reduction(+:Emp2)
   for ( int ii=0;ii<nparticles;++ii)
   {
//     index_t i = modelspace->particles[ii];
     index_t i = particles_vec[ii];
     //  std::cout << " i = " << i << std::endl;

     double ei = OpIn.OneBody(i,i);
     Orbit& oi = OpIn.modelspace->GetOrbit(i);
     for (auto& a : OpIn.modelspace->holes)
     {
       Orbit& oa = OpIn.modelspace->GetOrbit(a);
       double ea = OpIn.OneBody(a,a);
       double fia = OpIn.OneBody(i,a);
//       if (abs(fia)>1e-6)
        de = (oa.j2+1) * oa.occ * fia*fia/(ea-ei);
        orbit_impacts(a) += de;
        orbit_impacts(i) += de;

       for (index_t j : OpIn.modelspace->particles)
       {
         if (j<i) continue;
         double ej = OpIn.OneBody(j,j);
         Orbit& oj = OpIn.modelspace->GetOrbit(j);
         for ( auto& b: OpIn.modelspace->holes)
         {
           if (b<a) continue;
           Orbit& ob = OpIn.modelspace->GetOrbit(b);
           if ( (oi.l+oj.l+oa.l+ob.l)%2 >0) continue;
           if ( (oi.tz2 + oj.tz2) != (oa.tz2 +ob.tz2) ) continue;
           double eb = OpIn.OneBody(b,b);
           double denom = ea+eb-ei-ej;
           int Jmin = std::max(std::abs(oi.j2-oj.j2),std::abs(oa.j2-ob.j2))/2;
           int Jmax = std::min(oi.j2+oj.j2,oa.j2+ob.j2)/2;
           int dJ = 1;
           if (a==b or i==j)
           {
             Jmin += Jmin%2;
             dJ=2;
           }
           for (int J=Jmin; J<=Jmax; J+=dJ)
           {
             double tbme = OpIn.TwoBody.GetTBME_J_norm(J,a,b,i,j);
             if (std::abs(tbme)>1e-6)
             {
              de = (2*J+1)* oa.occ * ob.occ * tbme*tbme/denom; // no factor 1/4 because of the restricted sum

              orbit_impacts(a) += de;
              if (a!=b) orbit_impacts(b) += de;

              orbit_impacts(i) += de;
              if (i!=j) orbit_impacts(j) += de;

              }
           }
         }
       }
     }
   }
   IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
   return orbit_impacts;
}


