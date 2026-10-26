#include <iomanip>
#include "HFMBPT.hh"
#include "HartreeFock.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"
#include "Commutator.hh"
#include "UnitTest.hh" // for debugging purposes to use Mscheme

#include <omp.h>

HFMBPT::~HFMBPT()
{}

HFMBPT::HFMBPT(Operator& hbare)
  : HartreeFock(hbare),
    C_HO2NAT(modelspace->GetNumberOrbits(),modelspace->GetNumberOrbits(),arma::fill::eye),
    C_HF2NAT(modelspace->GetNumberOrbits(),modelspace->GetNumberOrbits(),arma::fill::eye),
    use_NAT_occupations(false), NAT_order("occupation"),vv_approx("imsrg2")
{
   vv_gen.SetType("white");
   vv_gen.SetDenominatorPartitioning("Moller_Plesset");
}

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

  int nchan = modelspace->GetNumberTwoBodyChannels();



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
            if ( 2*oa.n+oa.l+e2bra > modelspace->GetE3max() ) continue;
            for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
            {
              if ( std::abs(rho(a,b)) < 1e-8 ) continue; // Turns out this helps a bit (factor of 5 speed up in tests)
              Orbit & ob = modelspace->GetOrbit(b);
              if ( 2*ob.n+ob.l+e2ket > modelspace->GetE3max() ) continue;
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



/// Calculate the second-order perturbation theory correction to the energy
/// \f[
/// E^{(2)} = \sum_{ia} (2 j_a +1) \frac{|f_{ia}|^2}{f_{aa}-f_{ii}}
/// +  \sum_{\substack{i\leq j // a\leq b}}\sum_{J} (2J+1)\frac{|\Gamma_{ijab}^{J}|^2}{f_{aa}+f_{bb}-f_{ii}-f_{jj}}
/// \f]
///
double HFMBPT::GetMP2_Energy(const Operator& H) const
{
  double t_start = omp_get_wtime();
  double Emp2 = 0;
  int nparticles = H.modelspace->particles.size();
  std::vector<index_t> particles_vec(H.modelspace->particles.begin(), H.modelspace->particles.end()); // convert set to vector for OMP looping
  #pragma omp parallel for reduction(+:Emp2)
  for (int ii = 0; ii < nparticles; ++ii)
  {
    index_t i = particles_vec[ii];
    double ei = H.OneBody(i, i);
    Orbit &oi = H.modelspace->GetOrbit(i);
    for (auto &a : H.modelspace->holes)
    {
      Orbit &oa = H.modelspace->GetOrbit(a);
      double ea = H.OneBody(a, a);
      if (abs(H.OneBody(i, a)) > 1e-9)
        Emp2 += (oa.j2 + 1) * oa.occ * H.OneBody(i, a) * H.OneBody(a, i) / (ea-ei) ;
      for (index_t j : H.modelspace->particles)
      {
        if (j < i)
          continue;
        double ej = H.OneBody(j, j);
        Orbit &oj = H.modelspace->GetOrbit(j);
        for (auto &b : H.modelspace->holes)
        {
          if (b < a)
            continue;
          Orbit &ob = H.modelspace->GetOrbit(b);
          if ((oi.l + oj.l + oa.l + ob.l) % 2 > 0)
            continue;
          if ((oi.tz2 + oj.tz2) != (oa.tz2 + ob.tz2))
            continue;
          double eb = H.OneBody(b, b);
          double denom = ea + eb - ei - ej;
//          int Jmin = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
//          int Jmax = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
          int Jmin = AngMom::Jmin({ {oi.j2,oj.j2}, {oa.j2,ob.j2} }) /2;
          int Jmax = AngMom::Jmax({ {oi.j2,oj.j2}, {oa.j2,ob.j2} }) /2;
          int dJ = 1;
          if (a == b or i == j)
          {
            Jmin += Jmin % 2;
            dJ = 2;
          }
          for (int J = Jmin; J <= Jmax; J += dJ)
          {
            double tbme = H.TwoBody.GetTBME_J_norm(J, a, b, i, j);
            if (std::abs(tbme) > 1e-9)
            {
              Emp2 += (2 * J + 1) * oa.occ * ob.occ * tbme * tbme / denom; // no factor 1/4 because of the restricted sum
            }
          }
        }
      }
    }
  }
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return Emp2;
}

double HFMBPT::GetMP3_pp(const Operator& H) const
{
  double t_start = omp_get_wtime();
  double Epp = 0;
  // This can certainly be optimized, but I'll wait until this is the bottleneck.
  int nch = modelspace->GetNumberTwoBodyChannels();

  //   #pragma omp parallel for  schedule(dynamic,1) reduction(+:Emp3)
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Epp)
  for (int ich = 0; ich < nch; ++ich)
  {
    TwoBodyChannel &tbc = H.modelspace->GetTwoBodyChannel(ich);
    auto &Mat = H.TwoBody.GetMatrix(ich, ich);

    size_t n_pp = tbc.GetKetIndex_pp().size();
    size_t n_hh = tbc.GetKetIndex_hh().size();
    arma::mat M_hhpp(n_hh, n_pp, arma::fill::zeros);
    arma::mat M_pppp(n_pp, n_pp, arma::fill::zeros);

    size_t I_hh = 0;
    for (auto iket_ij : tbc.GetKetIndex_hh())
    {
      Ket &ket_ij = tbc.GetKet(iket_ij);
      index_t i = ket_ij.p;
      index_t j = ket_ij.q;

      size_t II_pp = 0;
      for (auto iket_ab : tbc.GetKetIndex_pp())
      {
        Ket &ket_ab = tbc.GetKet(iket_ab);
        index_t a = ket_ab.p;
        index_t b = ket_ab.q;
        double Delta_ijab = H.OneBody(i, i) + H.OneBody(j, j) - H.OneBody(a, a) - H.OneBody(b, b);
        M_hhpp(I_hh, II_pp) = Mat(iket_ij, iket_ab) / Delta_ijab;
        II_pp++;
      }
      I_hh++;
    }

    size_t I_pp = 0;
    for (auto iket_ab : tbc.GetKetIndex_pp())
    {
      size_t II_pp = 0;
      for (auto iket_cd : tbc.GetKetIndex_pp())
      {
        M_pppp(I_pp, II_pp) = Mat(iket_ab, iket_cd);
        II_pp++;
      }
      I_pp++;
    }

    int J = tbc.J;
    Epp += (2 * J + 1) * arma::trace(M_hhpp * M_pppp * M_hhpp.t());

  } // for ich

  
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return Epp;

}


double HFMBPT::GetMP3_hh(const Operator& H) const
{
  double t_start = omp_get_wtime();
  double Ehh = 0;
  int nch = modelspace->GetNumberTwoBodyChannels();

#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Ehh)
  for (int ich = 0; ich < nch; ++ich)
  {
    TwoBodyChannel &tbc = H.modelspace->GetTwoBodyChannel(ich);
    auto &Mat = H.TwoBody.GetMatrix(ich, ich);

    size_t n_hh = tbc.GetKetIndex_hh().size();
    size_t n_pp = tbc.GetKetIndex_pp().size();
    arma::mat M_hhpp(n_hh, n_pp, arma::fill::zeros);
    arma::mat M_hhhh(n_hh, n_hh, arma::fill::zeros);

    size_t I_hh = 0;
    for (auto iket_ij : tbc.GetKetIndex_hh())
    {
      Ket &ket_ij = tbc.GetKet(iket_ij);
      index_t i = ket_ij.p;
      index_t j = ket_ij.q;

      size_t II_pp = 0;
      for (auto iket_ab : tbc.GetKetIndex_pp())
      {
        Ket &ket_ab = tbc.GetKet(iket_ab);
        index_t a = ket_ab.p;
        index_t b = ket_ab.q;
        double Delta_ijab = H.OneBody(i, i) + H.OneBody(j, j) - H.OneBody(a, a) - H.OneBody(b, b);
        M_hhpp(I_hh, II_pp) = Mat(iket_ij, iket_ab) / Delta_ijab;
        II_pp++;
      }
      size_t II_hh = 0;
      for (auto iket_kl : tbc.GetKetIndex_hh())
      {
        M_hhhh(I_hh, II_hh) = Mat(iket_ij, iket_kl);
        II_hh++;
      }
      I_hh++;
    }
    int J = tbc.J;
    Ehh += (2 * J + 1) * arma::trace(M_hhpp.t() * M_hhhh * M_hhpp);
  } // for ich


  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return Ehh;
}

double HFMBPT::GetMP3_ph(const Operator& H) const
{
  double t_start = omp_get_wtime();
  double Eph = 0;

  H.modelspace->PreCalculateSixJ();

  int nch_CC = modelspace->GetNumberTwoBodyChannels_CC();

#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Eph)
  for (int ich_CC = 0; ich_CC < nch_CC; ++ich_CC)
  {
    TwoBodyChannel_CC &tbc_CC = H.modelspace->GetTwoBodyChannel_CC(ich_CC);
    size_t nkets_ph = tbc_CC.GetKetIndex_ph().size();
    arma::mat Vbar_iabj(nkets_ph, nkets_ph, arma::fill::zeros);
    arma::mat Vbar_bjck(nkets_ph, nkets_ph, arma::fill::zeros);
    int Jph = tbc_CC.J;

    size_t I_ph = 0;
    for (auto iket_ai : tbc_CC.GetKetIndex_ph())
    {
      Ket &ket_ai = tbc_CC.GetKet(iket_ai);
      index_t a = ket_ai.p;
      index_t i = ket_ai.q;
      double ja = 0.5 * H.modelspace->GetOrbit(a).j2;
      double ji = 0.5 * H.modelspace->GetOrbit(i).j2;

      int phase_ai = 1;
      int phase_ia = -AngMom::phase(ja + ji - Jph);
      if (ket_ai.op->occ < ket_ai.oq->occ)
      {
        std::swap(a, i);
        std::swap(ja, ji);
        std::swap(phase_ai, phase_ia);
      }

      size_t II_ph = 0;
      for (auto iket_bj : tbc_CC.GetKetIndex_ph())
      {
        Ket &ket_bj = tbc_CC.GetKet(iket_bj);
        index_t b = ket_bj.p;
        index_t j = ket_bj.q;

        double jb = 0.5 * H.modelspace->GetOrbit(b).j2;
        double jj = 0.5 * H.modelspace->GetOrbit(j).j2;

        int phase_bj = 1;
        int phase_jb = -AngMom::phase(jb + jj - Jph);
        if (ket_bj.op->occ < ket_bj.oq->occ)
        {
          std::swap(b, j);
          std::swap(jb, jj);
          std::swap(phase_bj, phase_jb);
        }

        double Delta_ijab = H.OneBody(i, i) + H.OneBody(j, j) - H.OneBody(a, a) - H.OneBody(b, b);
        int J1min = std::max(std::abs(ja - jb), std::abs(ji - jj));
        int J1max = std::min(ja + jb, ji + jj);
        double tbme_iabj = 0;
        double tbme_bjck = 0;
        //         double tbme_ckia = 0;

        if (AngMom::Triangle(jj, jb, Jph) and AngMom::Triangle(ji, ja, Jph))
        {
          for (int J1 = J1min; J1 <= J1max; ++J1) // Pandya 1: <ai`| V |jb`>_Jtot
          {
            tbme_iabj -= H.modelspace->GetSixJ(ja, ji, Jph, jj, jb, J1) * (2 * J1 + 1) * H.TwoBody.GetTBME_J(J1, i, j, b, a);
          }
        }

        J1min = std::max(std::abs(ji - jb), std::abs(ja - jj));
        J1max = std::min(ji + jb, ja + jj);

        if (AngMom::Triangle(jj, jb, Jph) and AngMom::Triangle(ji, ja, Jph))
        {
          for (int J1 = J1min; J1 <= J1max; ++J1) // Pandya 1: <ai`| V |jb`>_Jtot
          {
            tbme_bjck -= H.modelspace->GetSixJ(jb, jj, Jph, ja, ji, J1) * (2 * J1 + 1) * H.TwoBody.GetTBME_J(J1, b, i, a, j);
          }
        }

        Vbar_iabj(I_ph, II_ph) = tbme_iabj * phase_ia * phase_bj / Delta_ijab;
        Vbar_bjck(II_ph, I_ph) = tbme_bjck * phase_bj * phase_ai;
        II_ph++;
      }
      I_ph++;
    }
    auto Vbar_ckia = Vbar_iabj.t();
    Eph += (2 * Jph + 1) * arma::trace(Vbar_iabj * Vbar_bjck * Vbar_ckia);

  } // for ich_CC


  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return Eph;
}


double HFMBPT::GetMP3_Energy(const Operator& H) const
{
  double Epp = GetMP3_pp(H);
  double Ehh = GetMP3_hh(H);
  double Eph = GetMP3_ph(H);
  return Epp + Ehh + Eph;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//// Implement Van Vleck Perturbation Theory  ////

void HFMBPT::SetVVApprox( std::string apr)
{
   vv_approx = apr;
}

void HFMBPT::SetVVGenerator_type( std::string g)
{
   vv_gen.SetType(g);
}
void HFMBPT::SetVVGenerator_partitioning( std::string p)
{
   vv_gen.SetDenominatorPartitioning(p);
}

std::vector<Operator>& HFMBPT::GetVVOmegas()
{
   return vv_omegas;
}


////////////////////////////////////////////////////////////////////////////
/// Solve order-by-order for the generator Omega
/// in Van Vleck perturbation theory.
/// cf eq (41) in  J. Chem. Phys. 73, 5711–5717 (1980)
///  https://doi-org/10.1063/1.440050
void HFMBPT::SolveVanVleck(Operator& HNO, int order)
{
   vv_omegas.resize(order);
   for (int i=0; i<order; i++)
   {
      vv_omegas[i] = 0*HNO;
      vv_omegas[i].SetAntiHermitian();
   }
   Operator H0 = Operator(HNO);
   Operator VD = Operator(HNO);
   H0.EraseTwoBody();
   H0.EraseThreeBody();
   VD.ZeroBody = 0;
   VD.EraseOneBody();
   Operator VX = vv_gen.GetHod(VD);
   VD -= VX;

   // 1st order
   // [O1,H0] = -VX
   vv_gen.UpdateGeneral( VX, H0, vv_omegas[0] ); // first order Omega

   if ( order <2 )   return ;

   if ( vv_approx == "imsrg3n7")
   {
       Commutator::SetUseIMSRG3(true);
       Commutator::SetUseIMSRG3N7(true);
       for (int i=1;i<order;i++)
       {
          vv_omegas[i].ThreeBody.SetMode("pn");
          vv_omegas[i].SetParticleRank(3);
       }
   }

   // 2nd order [for a 2-body V, this goes up to 3-body]
   // [O2,H0] = [O1, (VD+1/2Vx)]
   Operator O1V1 = Commutator::Commutator( vv_omegas[0], (VD+1./2*VX) );
   vv_gen.UpdateGeneral( O1V1, H0, vv_omegas[1] ); // second order Omega
   if ( order <3 )   return ;

   // 3rd order [for a 2-body V, this should go up to 4-body]
   // [O3,H0] = [O2,(VD+1/2VX)] + 1/3[O1,[O1,VX]]  -1/4[O1, [O1,VX]X]
   Operator O2V1 = Commutator::Commutator( vv_omegas[1], VD+1./2*VX);
   Operator O1VX = Commutator::Commutator( vv_omegas[0], VX);
   Operator O1O1VX = Commutator::Commutator( vv_omegas[0], O1VX );
   Operator O1O1VXx = Commutator::Commutator( vv_omegas[0], vv_gen.GetHod( O1VX));
   if ( vv_approx == "imsrg3f2")
   {
      Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[0],VX, O1O1VX);
      Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[0],VX, O1O1VX);
   }
   Operator Hod = O2V1 + +1./3*O1O1VX - 1./4*O1O1VXx;
   vv_gen.UpdateGeneral( Hod,  H0, vv_omegas[2] ); // third order Omega

   if ( order <4 )   return ;
/*
   // 4th order [for a 2-body V, this should go up to 5-body]
   Operator O3V1   = Commutator::Commutator(vv_omegas[2],V1);
   Operator O2O1V1 = Commutator::Commutator(vv_omegas[1],O1V1);
   O2O1V1 += Commutator::Commutator(vv_omegas[0],O2V1); // Add in [O1,[O2,V]]
   if ( vv_approx == "imsrg3f2")
   {
      // To capture the cross terms [O1,[O2,V]], we do [O1+O2,[O1+O2,V]]
      // = [O1,[O1,V]] + [O1,[O2,V]] + [O2,[O1,V]] + [O2,[O2,V]].
      // We want the middle 2 terms, and we can subtract off the 1st and 4th.
      // In principle, we've already computed the [O1,[O1,V]] so we could
      // reuse it. But for simplicity I didn't bother.
      Operator O1O2 = vv_omegas[0]+vv_omegas[1]; // Omega_1 + Omega_2
      Commutator::FactorizedDoubleCommutator::comm223_231(O1O2,V1, O2O1V1);
      Commutator::FactorizedDoubleCommutator::comm223_232(O1O2,V1, O2O1V1);

      Operator nV = -V1;
      Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[0],nV, O2O1V1); // -[O1,[O1,V]]
      Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[0],nV, O2O1V1);
      Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[1],nV, O2O1V1); // -[O2,[O2,V]]
      Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[1],nV, O2O1V1);
   }

   Hod = O3V1 + 1./3*O2O1V1;
   vv_gen.UpdateGeneral( Hod, H0, vv_omegas[3]); // 4th order Omega
   if ( order <5 )   return ;
*/
   if ( order >=4 )
   {
      std::cout << "Orders beyond 4 not yet implemented. So you get up to 4th order." << std::endl;
   }

}


////////////////////////////////////////////////////////////////////////////
// The Hamiltonian is special because parts of it are driven to be diagonal
// We return H2,H3,H4,...
// In principle, H2 has a 3-body piece, H3 has a 4-body piece, H4 has as 5-body piece.
// To get the 0b-part correct to 4th order, we need 3f2 corrections for O1O1V and O1O1O1V
// As well as the perturbative triples, which accounts for contributions from the 3b part of O2
// to the 2b part of O3, which then contribute to H4 here. We can't use 3f2 for that because of
// the presence of energy denominators which spoil factorization.
// cf eq (46) in  J. Chem. Phys. 73, 5711–5717 (1980)
///  https://doi-org/10.1063/1.440050
std::vector<Operator> HFMBPT::VV_TransformH( Operator& HNO, int order, bool singleref=true)
{
   std::vector<Operator> Hout;
   if ( vv_approx == "imsrg3n7")
   {
       Commutator::SetUseIMSRG3(true);
       Commutator::SetUseIMSRG3N7(true);
   }
   Operator VX = vv_gen.GetHod(HNO);
   Hout.push_back( 1./2* Commutator::Commutator(vv_omegas[0],VX) ); // Hout[0] = H2
   Hout.push_back( 1./2* Commutator::Commutator(vv_omegas[1],VX) ); // Hout[1] = H3
   Hout.push_back( 1./2* Commutator::Commutator(vv_omegas[2],VX) ); // Hout[2] = H4
   Operator O1O1VX = Commutator::Commutator( 2*vv_omegas[0],Hout[0] ); 
   if (vv_approx == "imsrg3f2")
   {
     Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[0],VX, O1O1VX);
     if (not singleref)
        Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[0],VX, O1O1VX);
     Hout[1] += 1./12 * O1O1VX; // SRS added this. It's not in the Shavitt Redmon paper.
   }
   Operator O1O1O1V1 = Commutator::Commutator(vv_omegas[0], O1O1VX );
   // This part doesn't matter for the 4th order energy
   if (vv_approx == "imsrg3f2" and (not singleref))
   {
     Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[0],Hout[0], O1O1O1V1);
     Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[0],Hout[0], O1O1O1V1);
   }
   Hout[2] -= 1./24 * O1O1O1V1; 
   
   if ( vv_approx == "imsrg3f2")
   {
      Operator Wbar(HNO);
      Wbar.ZeroBody =0;
      int emax = HNO.modelspace->GetEmax();
      int E3maxsave = HNO.modelspace->GetE3max();
      HNO.modelspace->SetE3max(3*emax);
      Commutator::perturbative_triples = true;
      Commutator::comm223ss(vv_omegas[0], HNO, Wbar);
      Commutator::perturbative_triples = false;
      HNO.modelspace->SetE3max(E3maxsave);
      std::cout << "4th order triples = " << Wbar.ZeroBody << std::endl;
      Hout[2].ZeroBody += Wbar.ZeroBody;
   }

   return Hout;
}

////////////////////////////////////////////////////////////////////////////
// We implement the BCH transformation A' = A + [O,A] + 1/2[O,[O,A]] + ... order by order
//
//
std::vector<Operator> HFMBPT::VV_Transform(Operator& Xin, int order)
{
  if (order > vv_omegas.size() )
  {
     std::cout << "Uh oh. I only have Omega calculated to order " << vv_omegas.size()-1 << "  and you're asking for " << order << std::endl;
     order = vv_omegas.size();
  }
  std::vector<Operator> Xout;


  // 1st order [O1,X]
  Operator O1X = Commutator::Commutator(vv_omegas[0], Xin );
  Xout.push_back( O1X );
  if (order<2) return Xout;

  // 2nd order [O2,X] + 1/2 [O1,[O1,X]]
  Operator O2X = Commutator::Commutator(vv_omegas[1], Xin );
  Operator O1O1X = Commutator::Commutator(vv_omegas[0],O1X );

  if (vv_approx == "imsrg3f2")
  {
      Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[0],Xin, O1O1X);
      Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[0],Xin, O1O1X);
  }
  Xout.push_back( O2X + 1./2*O1O1X );
  if (order<3) return Xout;

  // 3rd order [O3,X] + 1/2 [O2,[O1,X]] + 1/2[O1,[O2,X]] + 1/6[O1,[O1,[O1,X]]]
  Operator O3X = Commutator::Commutator(vv_omegas[2], Xin );
  Operator O2O1X = Commutator::Commutator(vv_omegas[1], O1X );
  Operator O1O2X = Commutator::Commutator(vv_omegas[0], O2X );
  Operator O1O1O1X = Commutator::Commutator(vv_omegas[0], O1O1X );

  if (vv_approx == "imsrg3f2")
  {
      Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[0],O1X, O1O1O1X);
      Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[0],O1X, O1O1O1X);
      Operator O1O2 = vv_omegas[0]+vv_omegas[1]; // Omega_1 + Omega_2
      Commutator::FactorizedDoubleCommutator::comm223_231(O1O2,Xin, O2O1X);
      Commutator::FactorizedDoubleCommutator::comm223_232(O1O2,Xin, O2O1X);

      Operator nX = -Xin;
      Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[0],nX, O2O1X); // -[O1,[O1,V]]
      Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[0],nX, O2O1X);
      Commutator::FactorizedDoubleCommutator::comm223_231(vv_omegas[1],nX, O2O1X); // -[O2,[O2,V]]
      Commutator::FactorizedDoubleCommutator::comm223_232(vv_omegas[1],nX, O2O1X);
  }

  Xout.push_back( O3X + 1./2*(O2O1X+O1O2X) + 1./6*O1O1O1X );

//  // 4th order. dont do this for now



  return Xout;

}

double HFMBPT::GetDenom(const Operator& H, const std::vector<index_t>& holes, const std::vector<index_t>& particles) const
{
   double denom = 0;
   for ( auto& h : holes )
   {
      denom += H.OneBody(h,h);
   }
   for ( auto& p : particles )
   {
      denom -= H.OneBody(p,p);
   }
   return denom;

}

double HFMBPT::GetMP4_term( const Operator& H, int diagram) const
{
   double E = 0;
   if      (diagram==1)   E = GetMP4_F1(H);
   else if (diagram==2)   E = GetMP4_F2(H);
   else if (diagram==3)   E = GetMP4_F3(H);
   else if (diagram==4)   E = GetMP4_F4(H);
   else if (diagram==5)   E = GetMP4_F5(H);
   else if (diagram==6)   E = GetMP4_F6(H);
   else if (diagram==7)   E = GetMP4_F7(H);
   else if (diagram==8)   E = GetMP4_F8(H);
   else if (diagram==9)   E = GetMP4_F9(H);
   else if (diagram==10)   E = GetMP4_F10(H);
   else if (diagram==11)   E = GetMP4_F11(H);
   else if (diagram==12)   E = GetMP4_F12(H);
   else if (diagram==13)   E = GetMP4_F13(H);
   else if (diagram==14)   E = GetMP4_F14(H);
   else if (diagram==15)   E = GetMP4_F15(H);
   else if (diagram==16)   E = GetMP4_F16(H);
   else if (diagram==17)   E = GetMP4_F17(H);
   else if (diagram==18)   E = GetMP4_F18(H);
   else if (diagram==19)   E = GetMP4_F19(H);
   else if (diagram==20)   E = GetMP4_F20(H);
   else if (diagram==21)   E = GetMP4_F21(H);
   else if (diagram==22)   E = GetMP4_F22(H);
   else if (diagram==23)   E = GetMP4_F23(H);
   else if (diagram==24)   E = GetMP4_F24(H);
   else if (diagram==25)   E = GetMP4_F25(H);
   else if (diagram==26)   E = GetMP4_F26(H);
   else if (diagram==27)   E = GetMP4_F27(H);
   else if (diagram==28)   E = GetMP4_F28(H);
   else if (diagram==29)   E = GetMP4_F29(H);
   else if (diagram==30)   E = GetMP4_F30(H);
   else if (diagram==31)   E = GetMP4_F31(H);
   else if (diagram==32)   E = GetMP4_F32(H);
   else if (diagram==33)   E = GetMP4_F33(H);
   else if (diagram==34)   E = GetMP4_F34(H);
   else if (diagram==35)   E = GetMP4_F35(H);
   else if (diagram==36)   E = GetMP4_F36(H);
   else if (diagram==37)   E = GetMP4_F37(H);
   else if (diagram==38)   E = GetMP4_F38(H);
   else if (diagram==39)   E = GetMP4_F39(H);
   else
   {
      std::cout << __func__ << " Term " << diagram << " not yet implemented." << std::endl;
   }
   return E;
}

// Diagram F1 (as numbered by ADG)   corresponds to diagram 4 from Shavitt & Bartlett
// mscheme expression: F1 = 1/4 sum_abcijklm (v_abij v_ijak v_kclm v_lmbc) / (eps_abij eps_bk eps_bclm)
// agrees.

// Diagram F1 (as numbered by ADG)
// Factorized version: build two one-body-like intermediate matrices,
//   chi_A(b,k) = sum_J0 Jhat0^2 sum_{a,i,j} V_abij^J0 V_ijak^J0 / e_ijab
//   chi_B(k,b) = sum_J1 Jhat1^2 sum_{c,l,m} V_kclm^J1 V_lmbc^J1 / e_lmbc
// then F1 = 1/4 * jhat_b^-2 * sum_{b,k} delta_{jb,jk} chi_A(b,k) chi_B(k,b) / e_kb
double HFMBPT::GetMP4_F1( const Operator& H) const
{
   double t_start = omp_get_wtime();

   size_t norb = H.modelspace->GetNumberOrbits();
   arma::mat chi_A(norb,norb,arma::fill::zeros);  // indexed chi_A(b,k)
   arma::mat chi_B(norb,norb,arma::fill::zeros);  // indexed chi_B(k,b)

   // Build chi_A(b,k)
   for ( auto b : H.modelspace->particles )
   {
     Orbit& ob = H.modelspace->GetOrbit(b);
     for ( auto k : H.modelspace->holes )
     {
       Orbit& ok = H.modelspace->GetOrbit(k);
       double chi = 0;
       for ( auto a : H.modelspace->particles )
       {
         Orbit& oa = H.modelspace->GetOrbit(a);
         for ( auto i : H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j : H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             double e_ijab = GetDenom(H,{i,j},{a,b});
             int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,ok.j2}} ) /2;
             int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,ok.j2}} ) /2;
             for (int J0=J0_min; J0<=J0_max; J0++)
             {
               double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
               double vijak = H.TwoBody.GetTBME_J(J0,J0,i,j,a,k);
               chi += (2*J0+1) * vabij * vijak / e_ijab;
             }// for J0
           }// for j
         }// for i
       }// for a
       chi_A(b,k) = chi;
     }// for k
   }// for b

   // Build chi_B(k,b)
   for ( auto k : H.modelspace->holes )
   {
     Orbit& ok = H.modelspace->GetOrbit(k);
     for ( auto b : H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       double chi = 0;
       for ( auto c : H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto l : H.modelspace->holes )
         {
           Orbit& ol = H.modelspace->GetOrbit(l);
           for ( auto m : H.modelspace->holes )
           {
             Orbit& om = H.modelspace->GetOrbit(m);
             double e_lmbc = GetDenom(H,{l,m},{b,c});
             int J1_min = AngMom::Jmin( {{ok.j2,oc.j2},{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
             int J1_max = AngMom::Jmax( {{ok.j2,oc.j2},{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
             for (int J1=J1_min; J1<=J1_max; J1++)
             {
               double vkclm = H.TwoBody.GetTBME_J(J1,J1,k,c,l,m);
               double vlmbc = H.TwoBody.GetTBME_J(J1,J1,l,m,b,c);
               chi += (2*J1+1) * vkclm * vlmbc / e_lmbc;
             }// for J1
           }// for m
         }// for l
       }// for c
       chi_B(k,b) = chi;
     }// for b
   }// for k

   // Contract the two intermediates over the shared (b,k) pair
   double F1 =0;
   for ( auto b : H.modelspace->particles )
   {
     Orbit& ob = H.modelspace->GetOrbit(b);
     for ( auto k : H.modelspace->holes )
     {
       Orbit& ok = H.modelspace->GetOrbit(k);
       if ( ok.j2 != ob.j2) continue;
       double e_kb = GetDenom(H,{k},{b});
       F1 += 1./4 * chi_A(b,k) * chi_B(k,b) / ( e_kb * (ob.j2+1) );
     }// for k
   }// for b

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F1;
}

// Factorized version:
//   chi_A(b,k) = sum_J0 Jhat0^2 sum_{a,i,j} V_abij^J0 V_ijak^J0 / e_ijab
//   chi_C(b,k) = sum_J1 Jhat1^2 sum_{c,d,l} V_cdbl^J1 V_klcd^J1 / e_klcd
// then F2 = -1/4 * jhat_b^-2 * sum_{b,k} delta_{jb,jk} chi_A(b,k) chi_C(b,k) / e_kb
double HFMBPT::GetMP4_F2( const Operator& H) const
{
   double t_start = omp_get_wtime();

   size_t norb = H.modelspace->GetNumberOrbits();
   arma::mat chi_A(norb,norb,arma::fill::zeros);  // indexed chi_A(b,k)
   arma::mat chi_C(norb,norb,arma::fill::zeros);  // indexed chi_C(b,k)

   // Build chi_A(b,k)
   for ( auto b : H.modelspace->particles )
   {
     Orbit& ob = H.modelspace->GetOrbit(b);
     for ( auto k : H.modelspace->holes )
     {
       Orbit& ok = H.modelspace->GetOrbit(k);
       double chi = 0;
       for ( auto a : H.modelspace->particles )
       {
         Orbit& oa = H.modelspace->GetOrbit(a);
         for ( auto i : H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j : H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             double e_ijab = GetDenom(H,{i,j},{a,b});
             int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,ok.j2}} ) /2;
             int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,ok.j2}} ) /2;
             for (int J0=J0_min; J0<=J0_max; J0++)
             {
               double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
               double vijak = H.TwoBody.GetTBME_J(J0,J0,i,j,a,k);
               chi += (2*J0+1) * vabij * vijak / e_ijab;
             }// for J0
           }// for j
         }// for i
       }// for a
       chi_A(b,k) = chi;
     }// for k
   }// for b

   // Build chi_C(b,k)
   for ( auto b : H.modelspace->particles )
   {
     Orbit& ob = H.modelspace->GetOrbit(b);
     for ( auto k : H.modelspace->holes )
     {
       Orbit& ok = H.modelspace->GetOrbit(k);
       double chi = 0;
       for ( auto c : H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d : H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto l : H.modelspace->holes )
           {
             Orbit& ol = H.modelspace->GetOrbit(l);
             double e_klcd = GetDenom(H,{k,l},{c,d});
             int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{ob.j2,ol.j2},{ok.j2,ol.j2}} ) /2;
             int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{ob.j2,ol.j2},{ok.j2,ol.j2}} ) /2;
             for (int J1=J1_min; J1<=J1_max; J1++)
             {
               double vcdbl = H.TwoBody.GetTBME_J(J1,J1,c,d,b,l);
               double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
               chi += (2*J1+1) * vcdbl * vklcd / e_klcd;
             }// for J1
           }// for l
         }// for d
       }// for c
       chi_C(b,k) = chi;
     }// for k
   }// for b

   // Contract the two intermediates over the shared (b,k) pair
   double F2 =0;
   for ( auto b : H.modelspace->particles )
   {
     Orbit& ob = H.modelspace->GetOrbit(b);
     for ( auto k : H.modelspace->holes )
     {
       Orbit& ok = H.modelspace->GetOrbit(k);
       if ( ok.j2 != ob.j2) continue;
       double e_kb = GetDenom(H,{k},{b});
       F2 += -1./4 * chi_A(b,k) * chi_C(b,k) / ( e_kb * (ob.j2+1) );
     }// for k
   }// for b

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F2;
}



// Diagram F3 (as numbered by ADG)   complex conjugate diagram: F2
// Factorized version:
//   chi_D(j,c) = sum_J0 Jhat0^2 sum_{a,b,i} V_abij^J0 V_icab^J0 / e_ijab
//   chi_E(j,c) = sum_J1 Jhat1^2 sum_{d,k,l} V_jdkl^J1 V_klcd^J1 / e_klcd
// then F3 = -1/4 * jhat_c^-2 * sum_{j,c} delta_{jj,jc} chi_D(j,c) chi_E(j,c) / e_jc
double HFMBPT::GetMP4_F3( const Operator& H) const
{
   double F3 = GetMP4_F2(H); // conjugate diagram
   return F3;
}

// Diagram F4 (as numbered by ADG)
// Factorized version:
//   chi_D(j,c) = sum_J0 Jhat0^2 sum_{a,b,i} V_abij^J0 V_icab^J0 / e_ijab
//   chi_F(j,c) = sum_J1 Jhat1^2 sum_{d,e,k} V_deck^J1 V_jkde^J1 / e_jkde
// then F4 = 1/4 * jhat_c^-2 * sum_{j,c} delta_{jj,jc} chi_D(j,c) chi_F(j,c) / e_jc
double HFMBPT::GetMP4_F4( const Operator& H) const
{
   double t_start = omp_get_wtime();

   size_t norb = H.modelspace->GetNumberOrbits();
   arma::mat chi_D(norb,norb,arma::fill::zeros);  // indexed chi_D(j,c)
   arma::mat chi_F(norb,norb,arma::fill::zeros);  // indexed chi_F(j,c)

   // Build chi_D(j,c)
   for ( auto j : H.modelspace->holes )
   {
     Orbit& oj = H.modelspace->GetOrbit(j);
     for ( auto c : H.modelspace->particles )
     {
       Orbit& oc = H.modelspace->GetOrbit(c);
       double chi = 0;
       for ( auto a : H.modelspace->particles )
       {
         Orbit& oa = H.modelspace->GetOrbit(a);
         for ( auto b : H.modelspace->particles )
         {
           Orbit& ob = H.modelspace->GetOrbit(b);
           for ( auto i : H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             double e_ijab = GetDenom(H,{i,j},{a,b});
             int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oi.j2,oc.j2}} ) /2;
             int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oi.j2,oc.j2}} ) /2;
             for (int J0=J0_min; J0<=J0_max; J0++)
             {
               double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
               double vicab = H.TwoBody.GetTBME_J(J0,J0,i,c,a,b);
               chi += (2*J0+1) * vabij * vicab / e_ijab;
             }// for J0
           }// for i
         }// for b
       }// for a
       chi_D(j,c) = chi;
     }// for c
   }// for j

   // Build chi_F(j,c)
   for ( auto j : H.modelspace->holes )
   {
     Orbit& oj = H.modelspace->GetOrbit(j);
     for ( auto c : H.modelspace->particles )
     {
       Orbit& oc = H.modelspace->GetOrbit(c);
       double chi = 0;
       for ( auto d : H.modelspace->particles )
       {
         Orbit& od = H.modelspace->GetOrbit(d);
         for ( auto e : H.modelspace->particles )
         {
           Orbit& oe = H.modelspace->GetOrbit(e);
           for ( auto k : H.modelspace->holes )
           {
             Orbit& ok = H.modelspace->GetOrbit(k);
             double e_jkde = GetDenom(H,{j,k},{d,e});
             int J1_min = AngMom::Jmin( {{od.j2,oe.j2},{oc.j2,ok.j2},{oj.j2,ok.j2}} ) /2;
             int J1_max = AngMom::Jmax( {{od.j2,oe.j2},{oc.j2,ok.j2},{oj.j2,ok.j2}} ) /2;
             for (int J1=J1_min; J1<=J1_max; J1++)
             {
               double vdeck = H.TwoBody.GetTBME_J(J1,J1,d,e,c,k);
               double vjkde = H.TwoBody.GetTBME_J(J1,J1,j,k,d,e);
               chi += (2*J1+1) * vdeck * vjkde / e_jkde;
             }// for J1
           }// for k
         }// for e
       }// for d
       chi_F(j,c) = chi;
     }// for c
   }// for j

   // Contract the two intermediates over the shared (j,c) pair
   double F4 =0;
   for ( auto j : H.modelspace->holes )
   {
     Orbit& oj = H.modelspace->GetOrbit(j);
     for ( auto c : H.modelspace->particles )
     {
       Orbit& oc = H.modelspace->GetOrbit(c);
       if ( oj.j2 != oc.j2) continue;
       double e_jc = GetDenom(H,{j},{c});
       F4 += 1./4 * chi_D(j,c) * chi_F(j,c) / ( e_jc * (oc.j2+1) );
     }// for c
   }// for j

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F4;
}






// Diagram F5 (as numbered by ADG)
// mscheme expression: F5 = 1/16 sum_abijklmn (v_abij v_ijkl v_klmn v_mnab) / (eps^ij_ab eps^kl_ab eps^mn_ab)
// agrees.
double HFMBPT::GetMP4_F5( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F5 =0;

   size_t nch = H.modelspace->GetNumberTwoBodyChannels();
   for (size_t ch=0; ch<nch; ch++)
   {
      TwoBodyChannel& tbc = H.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      arma::uvec hh = tbc.GetKetIndex_hh();
      arma::uvec pp = tbc.GetKetIndex_pp();
      const arma::mat& V = H.TwoBody.MatEl.at({ch,ch});
      arma::mat Vpphh = V.submat(pp,hh);
      arma::mat Vhhhh = V.submat(hh,hh);
      arma::mat S = Vpphh;
      size_t npp = pp.n_rows;
      size_t nhh = hh.n_rows;
      for ( size_t ipp=0;ipp<npp; ipp++ )
      {
        auto& ab = pp(ipp);
        Ket& ket_ab = tbc.GetKet(ab);
        size_t a = ket_ab.p;
        size_t b = ket_ab.q;
        for ( size_t ihh=0;ihh<nhh; ihh++ )
        {
           auto& ij = hh(ihh);
           Ket& ket_ij = tbc.GetKet(ij);
           size_t i = ket_ij.p;
           size_t j = ket_ij.q;
           double e_ijab = GetDenom(H,{i,j},{a,b});
           S(ipp,ihh) /= e_ijab;
        }
      }
      S *= Vhhhh;
      for ( size_t ipp=0;ipp<npp; ipp++ )
      {
        auto& ab = pp(ipp);
        Ket& ket_ab = tbc.GetKet(ab);
        size_t a = ket_ab.p;
        size_t b = ket_ab.q;
        for ( size_t ihh=0;ihh<nhh; ihh++ )
        {
           auto& ij = hh(ihh);
           Ket& ket_ij = tbc.GetKet(ij);
           size_t i = ket_ij.p;
           size_t j = ket_ij.q;
           double e_ijab = GetDenom(H,{i,j},{a,b});
           S(ipp,ihh) /= e_ijab;
        }
      }
      S *= Vhhhh;
      for ( size_t ipp=0;ipp<npp; ipp++ )
      {
        auto& ab = pp(ipp);
        Ket& ket_ab = tbc.GetKet(ab);
        size_t a = ket_ab.p;
        size_t b = ket_ab.q;
        for ( size_t ihh=0;ihh<nhh; ihh++ )
        {
           auto& ij = hh(ihh);
           Ket& ket_ij = tbc.GetKet(ij);
           size_t i = ket_ij.p;
           size_t j = ket_ij.q;
           double e_ijab = GetDenom(H,{i,j},{a,b});
           S(ipp,ihh) /= e_ijab;
        }
      }
      S *= Vpphh.t(); // multiply by transpose
      // no factor 1/16 becase we store normalized matrix elements
      // and we only store kets |pq> with p<=q
      F5 += (2*J+1) * arma::trace(S);   
   }
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F5;
}



// Diagram F6 (as numbered by ADG)   complex conjugate diagram: F8
// mscheme expression: F6 = 1/2 sum_abcijklm (v_abij v_ijkl v_kcam v_lmbc) / (eps^ij_ab eps^kl_ab eps^lm_bc)
// agrees.
double HFMBPT::GetMP4_F6( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F6 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double e_lmbc = GetDenom(H,{l,m},{b,c});
                   double denom = e_ijab * e_klab * e_lmbc;
                   int phase_exponent = (ob.j2+ol.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ok.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{ok.j2,oc.j2},{oa.j2,om.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{ok.j2,oc.j2},{oa.j2,om.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijkl = H.TwoBody.GetTBME_J(J0,J0,i,j,k,l);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vkcam = H.TwoBody.GetTBME_J(J1,J1,k,c,a,m);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vlmbc = H.TwoBody.GetTBME_J(J2,J2,l,m,b,c);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, ok.j2/2., J2, om.j2/2., ol.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F6 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vijkl * vkcam * vlmbc / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F6;
}

// Diagram F7 (as numbered by ADG)   complex conjugate diagram: F14
// Factorized version, following the same matrix-chain approach as F5.
// F7 = 1/16 sum_J0 Jhat0^2 sum_{ab,ij,kl,cd} V_abij V_ijkl V_cdab V_klcd
//                                            / (e_ijab e_klab e_klcd)
// Chain (spectator pp pair ab through the hole-ladder step, then a
// particle-ladder step ab->cd with kl as the new spectator):
//   S1(ab,ij) = V_abij / e_ijab
//   S2(ab,kl) = [sum_ij S1(ab,ij) V_ijkl] / e_klab
//   S3(cd,kl) = [sum_ab V_cdab S2(ab,kl)] / e_klcd
//   F7 = sum_J0 Jhat0^2 sum_{cd,kl} S3(cd,kl) V_klcd(kl,cd)
// Since diagram 14 is the complex conjugate of this one, no separate
// routine is needed for F14.
double HFMBPT::GetMP4_F7( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F7 =0;

   size_t nch = H.modelspace->GetNumberTwoBodyChannels();
   for (size_t ch=0; ch<nch; ch++)
   {
      TwoBodyChannel& tbc = H.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      arma::uvec hh = tbc.GetKetIndex_hh();
      arma::uvec pp = tbc.GetKetIndex_pp();
      const arma::mat& V = H.TwoBody.MatEl.at({ch,ch});
      arma::mat Vpphh = V.submat(pp,hh);
      arma::mat Vhhhh = V.submat(hh,hh);
      arma::mat Vpppp = V.submat(pp,pp);
      arma::mat S = Vpphh;
      size_t npp = pp.n_rows;
      size_t nhh = hh.n_rows;

      // Stage 1: S(ab,ij) = V_abij / e_ijab
      for ( size_t ipp=0;ipp<npp; ipp++ )
      {
        auto& ab = pp(ipp);
        Ket& ket_ab = tbc.GetKet(ab);
        size_t a = ket_ab.p;
        size_t b = ket_ab.q;
        for ( size_t ihh=0;ihh<nhh; ihh++ )
        {
           auto& ij = hh(ihh);
           Ket& ket_ij = tbc.GetKet(ij);
           size_t i = ket_ij.p;
           size_t j = ket_ij.q;
           double e_ijab = GetDenom(H,{i,j},{a,b});
           S(ipp,ihh) /= e_ijab;
        }
      }

      // Stage 2: S(ab,kl) = [sum_ij S(ab,ij) V_ijkl] / e_klab
      S *= Vhhhh;
      for ( size_t ipp=0;ipp<npp; ipp++ )
      {
        auto& ab = pp(ipp);
        Ket& ket_ab = tbc.GetKet(ab);
        size_t a = ket_ab.p;
        size_t b = ket_ab.q;
        for ( size_t ihh=0;ihh<nhh; ihh++ )
        {
           auto& kl = hh(ihh);
           Ket& ket_kl = tbc.GetKet(kl);
           size_t k = ket_kl.p;
           size_t l = ket_kl.q;
           double e_klab = GetDenom(H,{k,l},{a,b});
           S(ipp,ihh) /= e_klab;
        }
      }

      // Stage 3: S(cd,kl) = [sum_ab V_cdab S(ab,kl)] / e_klcd
      S = Vpppp * S;
      for ( size_t ipp=0;ipp<npp; ipp++ )
      {
        auto& cd = pp(ipp);
        Ket& ket_cd = tbc.GetKet(cd);
        size_t c = ket_cd.p;
        size_t d = ket_cd.q;
        for ( size_t ihh=0;ihh<nhh; ihh++ )
        {
           auto& kl = hh(ihh);
           Ket& ket_kl = tbc.GetKet(kl);
           size_t k = ket_kl.p;
           size_t l = ket_kl.q;
           double e_klcd = GetDenom(H,{k,l},{c,d});
           S(ipp,ihh) /= e_klcd;
        }
      }

      // Close the ring: sum_{cd,kl} S(cd,kl) V_klcd(kl,cd)
      S *= Vpphh.t();
      // no factor 1/16 because we store normalized matrix elements
      // and we only store kets |pq> with p<=q
      F7 += (2*J+1) * arma::trace(S);
   }
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F7;
}

// Diagram F8 (as numbered by ADG)   complex conjugate diagram: F6
// mscheme expression: F8 = 1/2 sum_abcijklm (v_abij v_icak v_jklm v_lmbc) / (eps^ij_ab eps^jk_bc eps^lm_bc)
// agrees.
double HFMBPT::GetMP4_F8( const Operator& H) const
{
  double F8 =GetMP4_F6(H); // conjugate diagram
  return F8;
}



// Diagram F9 (as numbered by ADG)
// mscheme expression: F9 = sum_abcdijkl (v_abik v_icaj v_jdcl v_klbd) / (eps^ik_ab eps^jk_cb eps^kl_bd)
// agrees.
double HFMBPT::GetMP4_F9( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F9 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F9)
   for ( auto a : part_vec )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: part_vec )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: part_vec )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_jkcb = GetDenom(H,{j,k},{c,b});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ikab * e_jkcb * e_klbd;
                   int phase_exponent = (ob.j2+oc.j2+oj.j2+ok.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,oj.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,oj.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{oc.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{oc.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   // J4 does not label any two-body matrix element; its range comes purely from
                   // the triangle conditions of the four 6j symbols it appears in.
                   int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{oc.j2,oj.j2},{od.j2,ol.j2}} ) /2;
                   int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{oc.j2,oj.j2},{od.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicaj = H.TwoBody.GetTBME_J(J1,J1,i,c,a,j);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdcl = H.TwoBody.GetTBME_J(J2,J2,j,d,c,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oc.j2/2., oj.j2/2., J4, oa.j2/2., oi.j2/2., J1 );
                             double sixj3 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, ok.j2/2., ob.j2/2., J3 );
                             double sixj4 = H.modelspace->GetSixJ( ol.j2/2., od.j2/2., J4, oj.j2/2., oc.j2/2., J2 );
                             F9 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                   * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vicaj * vjdcl * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F9;
}

/*
// Implementation by Claude. It's very wrong.
// Diagram F9 (as numbered by ADG)
//
// Factorized (Pandya-transformed) version. Comparing the AMC-generated
// "cc" (cross-coupled) form of this diagram,
//   F9cc = sum_{abcdijkl,J4} Jhat4^2 * Vtilde_baik^{J4} Vtilde_icja^{J4}
//                                    * Vtilde_jdlc^{J4} Vtilde_lkbd^{J4}
// against the standard-coupled form (obtained from the original F9 by
// swapping index order within each vertex and using antisymmetry,
// V_abik = -V_baik etc.) term by term fixes each Vtilde uniquely:
//
//   Vtilde_baik(J4) = sum_J0 (2J0+1) * sixj(jk,jb,J4;ja,ji,J0) * V_baik^{J0J0}
//   Vtilde_icja(J4) = sum_J1 (2J1+1) * sixj(jc,jj,J4;ja,ji,J1) * V_icja^{J1J1}
//   Vtilde_jdlc(J4) = sum_J2 (2J2+1) * sixj(jl,jd,J4;jj,jc,J2) * V_jdlc^{J2J2}
//   Vtilde_lkbd(J4) = sum_J3 (2J3+1) * sixj(jd,jl,J4;jk,jb,J3) * V_lkbd^{J3J3}
//
// Each Vtilde is a matrix in the ph-channel space of pairs coupled to
// J4: Vtilde_baik has rows=(b,k), cols=(a,i); Vtilde_icja has
// rows=(a,i), cols=(c,j); Vtilde_jdlc has rows=(c,j), cols=(d,l);
// Vtilde_lkbd has rows=(d,l), cols=(b,k), closing the ring. So
//   F9 = sum_J4 Jhat4^2 * trace( Vtilde_baik * Vtilde_icja * Vtilde_jdlc * Vtilde_lkbd )
double HFMBPT::GetMP4_F9( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F9 =0;

   H.modelspace->PreCalculateSixJ();

   int nch_CC = H.modelspace->GetNumberTwoBodyChannels_CC();

   for (int ich_CC=0; ich_CC<nch_CC; ich_CC++)
   {
      TwoBodyChannel_CC& tbc_CC = H.modelspace->GetTwoBodyChannel_CC(ich_CC);
      int J4 = tbc_CC.J;
      size_t nkets_ph = tbc_CC.GetKetIndex_ph().size();

      // Collect the (particle,hole) orbital for every ph-ket in this
      // channel once, so all four matrices below share one basis.
      std::vector<index_t> ph_particle(nkets_ph);
      std::vector<index_t> ph_hole(nkets_ph);
      size_t idx=0;
      for ( auto iket : tbc_CC.GetKetIndex_ph() )
      {
        Ket& ket = tbc_CC.GetKet(iket);
        index_t p = ket.p;
        index_t q = ket.q;
        // occ=0 for particles, occ=1 for holes -- pick out which is which.
        if ( H.modelspace->GetOrbit(p).occ < H.modelspace->GetOrbit(q).occ )
        {
          ph_particle[idx] = p;
          ph_hole[idx] = q;
        }
        else
        {
          ph_particle[idx] = q;
          ph_hole[idx] = p;
        }
        idx++;
      }

      arma::mat Vtilde_baik(nkets_ph,nkets_ph,arma::fill::zeros); // rows=(b,k), cols=(a,i)
      arma::mat Vtilde_icja(nkets_ph,nkets_ph,arma::fill::zeros); // rows=(a,i), cols=(c,j)
      arma::mat Vtilde_jdlc(nkets_ph,nkets_ph,arma::fill::zeros); // rows=(c,j), cols=(d,l)
      arma::mat Vtilde_lkbd(nkets_ph,nkets_ph,arma::fill::zeros); // rows=(d,l), cols=(b,k)

      for ( size_t irow=0; irow<nkets_ph; irow++ )
      {
        for ( size_t icol=0; icol<nkets_ph; icol++ )
        {
          // --- Vtilde_baik(b,k;a,i) ---
                   double denom = e_ikab * e_jkcb * e_klbd;
          {
            index_t b = ph_particle[irow];
            index_t k = ph_hole[irow];
            index_t a = ph_particle[icol];
            index_t i = ph_hole[icol];
            Orbit& oa = H.modelspace->GetOrbit(a);
            Orbit& ob = H.modelspace->GetOrbit(b);
            Orbit& oi = H.modelspace->GetOrbit(i);
            Orbit& ok = H.modelspace->GetOrbit(k);
            double ja=oa.j2/2., jb=ob.j2/2., ji=oi.j2/2., jk=ok.j2/2.;
            int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
            int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
            double val = 0;
            for (int J0=J0_min; J0<=J0_max; J0++)
            {
              double sixj = H.modelspace->GetSixJ(jk,jb,J4,ja,ji,J0);
              double v = H.TwoBody.GetTBME_J(J0,J0,b,a,i,k);
              val += (2*J0+1) * sixj * v;
            }
            Vtilde_baik(irow,icol) = val / e_ikab;
          }

          // --- Vtilde_icja(a,i;c,j) ---
          {
            index_t a = ph_particle[irow];
            index_t i = ph_hole[irow];
            index_t c = ph_particle[icol];
            index_t j = ph_hole[icol];
            Orbit& oa = H.modelspace->GetOrbit(a);
            Orbit& oi = H.modelspace->GetOrbit(i);
            Orbit& oc = H.modelspace->GetOrbit(c);
            Orbit& oj = H.modelspace->GetOrbit(j);
            double ja=oa.j2/2., ji=oi.j2/2., jc=oc.j2/2., jj=oj.j2/2.;
            int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,oa.j2}} ) /2;
            int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,oa.j2}} ) /2;
            double val = 0;
            for (int J1=J1_min; J1<=J1_max; J1++)
            {
              double sixj = H.modelspace->GetSixJ(jc,jj,J4,ja,ji,J1);
              double v = H.TwoBody.GetTBME_J(J1,J1,i,c,j,a);
              val += (2*J1+1) * sixj * v;
            }
            Vtilde_icja(irow,icol) = val;
          }

          // --- Vtilde_jdlc(c,j;d,l) ---
          {
            index_t c = ph_particle[irow];
            index_t j = ph_hole[irow];
            index_t d = ph_particle[icol];
            index_t l = ph_hole[icol];
            Orbit& oc = H.modelspace->GetOrbit(c);
            Orbit& oj = H.modelspace->GetOrbit(j);
            Orbit& od = H.modelspace->GetOrbit(d);
            Orbit& ol = H.modelspace->GetOrbit(l);
            double jc=oc.j2/2., jj=oj.j2/2., jd=od.j2/2., jl=ol.j2/2.;
            int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{ol.j2,oc.j2}} ) /2;
            int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{ol.j2,oc.j2}} ) /2;
            double val = 0;
            for (int J2=J2_min; J2<=J2_max; J2++)
            {
              double sixj = H.modelspace->GetSixJ(jl,jd,J4,jj,jc,J2);
              double v = H.TwoBody.GetTBME_J(J2,J2,j,d,l,c);
              val += (2*J2+1) * sixj * v;
            }
            Vtilde_jdlc(irow,icol) = val;
          }

                   double e_jkcb = GetDenom(H,{j,k},{c,b});
          // --- Vtilde_lkbd(d,l;b,k) ---
          {
            index_t d = ph_particle[irow];
            index_t l = ph_hole[irow];
            index_t b = ph_particle[icol];
            index_t k = ph_hole[icol];
            Orbit& od = H.modelspace->GetOrbit(d);
            Orbit& ol = H.modelspace->GetOrbit(l);
            Orbit& ob = H.modelspace->GetOrbit(b);
            Orbit& ok = H.modelspace->GetOrbit(k);
            double e_klbd = GetDenom(H,{k,l},{b,d});
            double jd=od.j2/2., jl=ol.j2/2., jb=ob.j2/2., jk=ok.j2/2.;
            int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
            int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
            double val = 0;
            for (int J3=J3_min; J3<=J3_max; J3++)
            {
              double sixj = H.modelspace->GetSixJ(jd,jl,J4,jk,jb,J3);
              double v = H.TwoBody.GetTBME_J(J3,J3,l,k,b,d);
              val += (2*J3+1) * sixj * v;
            }
            Vtilde_lkbd(irow,icol) = val / e_klbd;
          }

        }// for icol
      }// for irow
      std::cout << " ich_CC = " << ich_CC << "  four Pandya-transformed matrices: " << std::endl
                << Vtilde_baik << std::endl << std::endl
                << Vtilde_icja << std::endl << std::endl
                << Vtilde_jdlc << std::endl << std::endl
                << Vtilde_lkbd << std::endl;

      F9 += (2*J4+1) * arma::trace( Vtilde_baik * Vtilde_icja * Vtilde_jdlc * Vtilde_lkbd );

   }// for ich_CC

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F9;
}
*/


// Diagram F10 (as numbered by ADG)
// mscheme expression: F10 = -sum_abcdijkl (v_abij v_icak v_jdcl v_klbd) / (eps^ij_ab eps^jk_cb eps^kl_bd)
// POSSIBLE SIGN ERROR?
double HFMBPT::GetMP4_F10( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F10 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F10)
//   #pragma omp parallel for
   for ( auto a : part_vec )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: part_vec )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: part_vec )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jkcb = GetDenom(H,{j,k},{c,b});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ijab * e_jkcb * e_klbd;
                   int phase_exponent = (od.j2+ol.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{oc.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{oc.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicak = H.TwoBody.GetTBME_J(J1,J1,i,c,a,k);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdcl = H.TwoBody.GetTBME_J(J2,J2,j,d,c,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           // J4 does not label any two-body matrix element. Two of the 6j symbols
                           // couple it to single-particle j's, but the other two couple it to J2
                           // and J3 directly, so its range can only be fixed once J2,J3 are known.
                           int J4_min = AngMom::Jmin( {{oj.j2,ob.j2},{oc.j2,ok.j2},{2*J2,2*J3}} ) /2;
                           int J4_max = AngMom::Jmax( {{oj.j2,ob.j2},{oc.j2,ok.j2},{2*J2,2*J3}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( oj.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oc.j2/2., ok.j2/2., J4, oa.j2/2., oi.j2/2., J1 );
                             double sixj3 = H.modelspace->GetSixJ( J3, J4, J2, oj.j2/2., od.j2/2., ob.j2/2. );
                             double sixj4 = H.modelspace->GetSixJ( J4, J2, J3, ol.j2/2., ok.j2/2., oc.j2/2. );
                             F10 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabij * vicak * vjdcl * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F10;
}



// Diagram F11 (as numbered by ADG)
// mscheme expression: F11 = -sum_abcdijkl (v_abik v_icaj v_jdbl v_klcd) / (eps^ik_ab eps^jk_bc eps^kl_cd)
// POSSIBLE SIGN ERROR?
double HFMBPT::GetMP4_F11( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F11 =0;
//   #pragma omp parallel for
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F11)
   for ( auto a : part_vec )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: part_vec )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: part_vec )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_jkbc = GetDenom(H,{j,k},{b,c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ikab * e_jkbc * e_klcd;
                   int phase_exponent = (od.j2+ol.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,oj.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,oj.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{ob.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{ob.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicaj = H.TwoBody.GetTBME_J(J1,J1,i,c,a,j);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdbl = H.TwoBody.GetTBME_J(J2,J2,j,d,b,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklcd = H.TwoBody.GetTBME_J(J3,J3,k,l,c,d);
                           // As in F10, J4 is tied to J2,J3 directly through two of the 6j symbols.
                           int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{oc.j2,oj.j2},{2*J2,2*J3}} ) /2;
                           int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{oc.j2,oj.j2},{2*J2,2*J3}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oc.j2/2., oj.j2/2., J4, oa.j2/2., oi.j2/2., J1 );
                             double sixj3 = H.modelspace->GetSixJ( J3, J4, J2, ob.j2/2., ol.j2/2., ok.j2/2. );
                             double sixj4 = H.modelspace->GetSixJ( J3, J2, J4, oj.j2/2., oc.j2/2., od.j2/2. );
                             F11 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vicaj * vjdbl * vklcd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F11;
}



// Diagram F12 (as numbered by ADG)
// mscheme expression: F12 = sum_abcdijkl (v_abij v_icak v_jdbl v_klcd) / (eps^ij_ab eps^jk_bc eps^kl_cd)
// Fixed the triangle condidions. Now it works.
double HFMBPT::GetMP4_F12( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F12 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F12)
   for ( auto a : part_vec )
   {
     for ( auto b: part_vec )
     {
       for ( auto c: part_vec )
       {
         Orbit& ob = H.modelspace->GetOrbit(b);
         Orbit& oa = H.modelspace->GetOrbit(a);
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jkbc = GetDenom(H,{j,k},{b,c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_jkbc * e_klcd;
                   int phase_exponent = (ob.j2+oc.j2+oj.j2+ok.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{ob.j2,ol.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{ob.j2,ol.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   // Here all four 6j symbols tie J4 to single-particle j's only.
                   // SRS: J4 has triangle cond with (j,b) (a,i) (c,k) (l,d).  It appears Claude messed up the conditions
//                   int J4_min = AngMom::Jmin( {{oj.j2,ob.j2},{oc.j2,ok.j2},{ol.j2,od.j2},{ol.j2,ok.j2}} ) /2;
//                   int J4_max = AngMom::Jmax( {{oj.j2,ob.j2},{oc.j2,ok.j2},{ol.j2,od.j2},{ol.j2,ok.j2}} ) /2;
                   int J4_min = AngMom::Jmin( {{oj.j2,ob.j2},{oc.j2,ok.j2},{ol.j2,od.j2},{oa.j2,oi.j2}} ) /2;
                   int J4_max = AngMom::Jmax( {{oj.j2,ob.j2},{oc.j2,ok.j2},{ol.j2,od.j2},{oa.j2,oi.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicak = H.TwoBody.GetTBME_J(J1,J1,i,c,a,k);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdbl = H.TwoBody.GetTBME_J(J2,J2,j,d,b,l);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklcd = H.TwoBody.GetTBME_J(J3,J3,k,l,c,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( oj.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oc.j2/2., ok.j2/2., J4, oa.j2/2., oi.j2/2., J1 );
                             double sixj3 = H.modelspace->GetSixJ( ol.j2/2., od.j2/2., J4, oj.j2/2., ob.j2/2., J2 );
                             double sixj4 = H.modelspace->GetSixJ( oc.j2/2., od.j2/2., J3, ol.j2/2., ok.j2/2., J4 );
                             F12 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabij * vicak * vjdbl * vklcd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F12;
}



// Diagram F13 (as numbered by ADG)   complex conjugate diagram: F15
// mscheme expression: F13 = 1/2 sum_abcdeijk (v_abij v_icak v_debc v_jkde) / (eps^ij_ab eps^jk_bc eps^jk_de)
// agrees.
double HFMBPT::GetMP4_F13( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F13 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jkbc = GetDenom(H,{j,k},{b,c});
                   double e_jkde = GetDenom(H,{j,k},{d,e});
                   double denom = e_ijab * e_jkbc * e_jkde;
                   int phase_exponent = (ob.j2+oj.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{od.j2,oe.j2},{ob.j2,oc.j2},{oj.j2,ok.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{od.j2,oe.j2},{ob.j2,oc.j2},{oj.j2,ok.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicak = H.TwoBody.GetTBME_J(J1,J1,i,c,a,k);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vdebc = H.TwoBody.GetTBME_J(J2,J2,d,e,b,c);
                         double vjkde = H.TwoBody.GetTBME_J(J2,J2,j,k,d,e);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, oi.j2/2., J2, ok.j2/2., oj.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F13 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vicak * vdebc * vjkde / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F13;
}



// Diagram F14 (as numbered by ADG)   complex conjugate diagram: F7
// mscheme expression: F14 = 1/16 sum_abcdijkl (v_abij v_cdab v_ijkl v_klcd) / (eps^ij_ab eps^ij_cd eps^kl_cd)
// agrees.
double HFMBPT::GetMP4_F14( const Operator& H) const
{
   double F14 = GetMP4_F7(H); // conjugate diagram
   return F14;
}


// Diagram F15 (as numbered by ADG)   complex conjugate diagram: F13
// mscheme expression: F15 = 1/2 sum_abcdeijk (v_abij v_cdab v_ieck v_jkde) / (eps^ij_ab eps^ij_cd eps^jk_de)
// agrees.
double HFMBPT::GetMP4_F15( const Operator& H) const
{

   double F15 = GetMP4_F13(H); // conjugate diagram
   return F15;
}



// Diagram F16 (as numbered by ADG)
// Factorized version, following the same matrix-chain approach as F5,
// but with the hole pair ij as the spectator throughout, and a chain
// of three particle-particle ladder steps ab->cd->ef.
// F16 = 1/16 sum_J0 Jhat0^2 sum_{ab,cd,ef,ij} V_abij V_cdab V_efcd V_ijef
//                                             / (e_ijab e_ijcd e_ijef)
// Chain:
//   S1(ij,ab) = V_abij / e_ijab
//   S2(ij,cd) = [sum_ab S1(ij,ab) V_abcd] / e_ijcd
//   S3(ij,ef) = [sum_cd S2(ij,cd) V_cdef] / e_ijef
//   F16 = sum_J0 Jhat0^2 sum_{ij,ef} S3(ij,ef) V_ijef(ij,ef)
double HFMBPT::GetMP4_F16( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F16 =0;

   size_t nch = H.modelspace->GetNumberTwoBodyChannels();
   for (size_t ch=0; ch<nch; ch++)
   {
      TwoBodyChannel& tbc = H.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      arma::uvec hh = tbc.GetKetIndex_hh();
      arma::uvec pp = tbc.GetKetIndex_pp();
      const arma::mat& V = H.TwoBody.MatEl.at({ch,ch});
      arma::mat Vpphh = V.submat(pp,hh);
      arma::mat Vpppp = V.submat(pp,pp);
      arma::mat S = Vpphh.t();  // S(ij,ab) = V_abij
      size_t npp = pp.n_rows;
      size_t nhh = hh.n_rows;

      // Stage 1: S(ij,ab) = V_abij / e_ijab
      for ( size_t ihh=0;ihh<nhh; ihh++ )
      {
        auto& ij = hh(ihh);
        Ket& ket_ij = tbc.GetKet(ij);
        size_t i = ket_ij.p;
        size_t j = ket_ij.q;
        for ( size_t ipp=0;ipp<npp; ipp++ )
        {
           auto& ab = pp(ipp);
           Ket& ket_ab = tbc.GetKet(ab);
           size_t a = ket_ab.p;
           size_t b = ket_ab.q;
           double e_ijab = GetDenom(H,{i,j},{a,b});
           S(ihh,ipp) /= e_ijab;
        }
      }

      // Stage 2: S(ij,cd) = [sum_ab S(ij,ab) V_abcd] / e_ijcd
      S = S * Vpppp;
      for ( size_t ihh=0;ihh<nhh; ihh++ )
      {
        auto& ij = hh(ihh);
        Ket& ket_ij = tbc.GetKet(ij);
        size_t i = ket_ij.p;
        size_t j = ket_ij.q;
        for ( size_t ipp=0;ipp<npp; ipp++ )
        {
           auto& cd = pp(ipp);
           Ket& ket_cd = tbc.GetKet(cd);
           size_t c = ket_cd.p;
           size_t d = ket_cd.q;
           double e_ijcd = GetDenom(H,{i,j},{c,d});
           S(ihh,ipp) /= e_ijcd;
        }
      }

      // Stage 3: S(ij,ef) = [sum_cd S(ij,cd) V_cdef] / e_ijef
      S = S * Vpppp;
      for ( size_t ihh=0;ihh<nhh; ihh++ )
      {
        auto& ij = hh(ihh);
        Ket& ket_ij = tbc.GetKet(ij);
        size_t i = ket_ij.p;
        size_t j = ket_ij.q;
        for ( size_t ipp=0;ipp<npp; ipp++ )
        {
           auto& ef = pp(ipp);
           Ket& ket_ef = tbc.GetKet(ef);
           size_t e = ket_ef.p;
           size_t f = ket_ef.q;
           double e_ijef = GetDenom(H,{i,j},{e,f});
           S(ihh,ipp) /= e_ijef;
        }
      }

      // Close the ring: sum_{ij,ef} S(ij,ef) V_ijef(ij,ef)
      S = S * Vpphh;
      // no factor 1/16 because we store normalized matrix elements
      // and we only store kets |pq> with p<=q
      F16 += (2*J+1) * arma::trace(S);
   }
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F16;
}





// Diagram F17 (as numbered by ADG)
// mscheme expression: F17 = 1/4 sum_abcijklm (v_abil v_icjk v_jkcm v_lmab) / (eps^il_ab eps^jkl_cab eps^lm_ab)
// Claude was OFF BY A MINUS SIGN. I fixed it.
double HFMBPT::GetMP4_F17( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F17 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   if ( om.j2 != oi.j2) continue;
                   double e_ilab = GetDenom(H,{i,l},{a,b});
                   double e_jklcab = GetDenom(H,{j,k,l},{c,a,b});
                   double e_lmab = GetDenom(H,{l,m},{a,b});
                   double denom = e_ilab * e_jklcab * e_lmab;
                   int phase_exponent = (oc.j2+ol.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ol.j2},{ol.j2,om.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ol.j2},{ol.j2,om.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,ok.j2},{oc.j2,om.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,ok.j2},{oc.j2,om.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabil = H.TwoBody.GetTBME_J(J0,J0,a,b,i,l);
                     double vlmab = H.TwoBody.GetTBME_J(J0,J0,l,m,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicjk = H.TwoBody.GetTBME_J(J1,J1,i,c,j,k);
                       double vjkcm = H.TwoBody.GetTBME_J(J1,J1,j,k,c,m);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F17 -= 1./4 * phase * (2*J0+1) * (2*J1+1) / (oi.j2+1) * vabil * vicjk * vjkcm * vlmab / denom;
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F17;
}



// Diagram F18 (as numbered by ADG)
// mscheme expression: F18 = 1/2 sum_abcijklm (v_abij v_ickl v_jkcm v_lmab) / (eps^ij_ab eps^jkl_cab eps^lm_ab)
// agrees.
double HFMBPT::GetMP4_F18( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F18 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jklcab = GetDenom(H,{j,k,l},{c,a,b});
                   double e_lmab = GetDenom(H,{l,m},{a,b});
                   double denom = e_ijab * e_jklcab * e_lmab;
                   int phase_exponent = (oc.j2+oj.j2+ok.j2+ol.j2)/2;  // + J0, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ol.j2,om.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{ol.j2,om.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{ok.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,ok.j2},{oc.j2,om.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,ok.j2},{oc.j2,om.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vlmab = H.TwoBody.GetTBME_J(J0,J0,l,m,a,b);
                     double phase = ( (J0+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vickl = H.TwoBody.GetTBME_J(J1,J1,i,c,k,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjkcm = H.TwoBody.GetTBME_J(J2,J2,j,k,c,m);
                         double ninej = H.modelspace->GetNineJ( oj.j2/2., oi.j2/2., J0, ok.j2/2., J1, ol.j2/2., J2, oc.j2/2., om.j2/2. );
                         F18 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vickl * vjkcm * vlmab / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F18;
}



// Diagram F19 (as numbered by ADG)
// mscheme expression: F19 = 1/2 sum_abcijklm (v_abil v_icjk v_jkam v_lmbc) / (eps^il_ab eps^jkl_abc eps^lm_bc)
// agrees.
double HFMBPT::GetMP4_F19( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F19 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ilab = GetDenom(H,{i,l},{a,b});
                   double e_jklabc = GetDenom(H,{j,k,l},{a,b,c});
                   double e_lmbc = GetDenom(H,{l,m},{b,c});
                   double denom = e_ilab * e_jklabc * e_lmbc;
                   int phase_exponent = (ob.j2+ol.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ol.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,ok.j2},{oa.j2,om.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,ok.j2},{oa.j2,om.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabil = H.TwoBody.GetTBME_J(J0,J0,a,b,i,l);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicjk = H.TwoBody.GetTBME_J(J1,J1,i,c,j,k);
                       double vjkam = H.TwoBody.GetTBME_J(J1,J1,j,k,a,m);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vlmbc = H.TwoBody.GetTBME_J(J2,J2,l,m,b,c);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, oi.j2/2., J2, om.j2/2., ol.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F19 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabil * vicjk * vjkam * vlmbc / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F19;
}



// Diagram F20 (as numbered by ADG)
// mscheme expression: F20 = sum_abcijklm (v_abij v_ickl v_jkam v_lmbc) / (eps^ij_ab eps^jkl_abc eps^lm_bc)
// doesnt agree, and I don't know why
double HFMBPT::GetMP4_F20( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F20 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F20)
   for ( auto a : part_vec )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: part_vec )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: part_vec )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto i: H.modelspace->holes )
         {
           Orbit& oi = H.modelspace->GetOrbit(i);
           for ( auto j: H.modelspace->holes )
           {
             Orbit& oj = H.modelspace->GetOrbit(j);
             for ( auto k: H.modelspace->holes )
             {
               Orbit& ok = H.modelspace->GetOrbit(k);
               for ( auto l: H.modelspace->holes )
               {
                 Orbit& ol = H.modelspace->GetOrbit(l);
                 for ( auto m: H.modelspace->holes )
                 {
                   Orbit& om = H.modelspace->GetOrbit(m);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jklabc = GetDenom(H,{j,k,l},{a,b,c});
                   double e_lmbc = GetDenom(H,{l,m},{b,c});
                   double denom = e_ijab * e_jklabc * e_lmbc;
                   int phase_exponent = (oc.j2+oi.j2+oj.j2+om.j2)/2;  // + J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{ok.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,ok.j2},{oa.j2,om.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,ok.j2},{oa.j2,om.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ol.j2,om.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vickl = H.TwoBody.GetTBME_J(J1,J1,i,c,k,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjkam = H.TwoBody.GetTBME_J(J2,J2,j,k,a,m);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vlmbc = H.TwoBody.GetTBME_J(J3,J3,l,m,b,c);
                           double phase = ( (J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           // J4 does not label any two-body matrix element. Two of the 6j symbols
                           // couple it to single-particle j's; the other two couple it to J1 and J3.
                           int J4_min = AngMom::Jmin( {{oi.j2,ob.j2},{oj.j2,oa.j2},{om.j2,ok.j2},{2*J1,2*J3}} ) /2;
                           int J4_max = AngMom::Jmax( {{oi.j2,ob.j2},{oj.j2,oa.j2},{om.j2,ok.j2},{2*J1,2*J3}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( oj.j2/2., oa.j2/2., J4, ob.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( J1, J3, J4, ob.j2/2., oi.j2/2., oc.j2/2. );
                             double sixj3 = H.modelspace->GetSixJ( om.j2/2., ok.j2/2., J4, oj.j2/2., oa.j2/2., J2 );
                             double sixj4 = H.modelspace->GetSixJ( J1, J3, J4, om.j2/2., ok.j2/2., ol.j2/2. );
                             F20 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabij * vickl * vjkam * vlmbc / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for m
               }// for l
             }// for k
           }// for j
         }// for i
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F20;
}



// Diagram F21 (as numbered by ADG)   complex conjugate diagram: F25
// mscheme expression: F21 = -sum_abcdijkl (v_abik v_icjl v_jdac v_klbd) / (eps^ik_ab eps^jkl_acb eps^kl_bd)
// agrees.
double HFMBPT::GetMP4_F21( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F21 =0;
//   #pragma omp parallel for
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F21)
   for ( auto a : part_vec )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: part_vec )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: part_vec )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_jklacb = GetDenom(H,{j,k,l},{a,c,b});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ikab * e_jklacb * e_klbd;
                   int phase_exponent = (oa.j2+ob.j2+oc.j2+od.j2+oi.j2+oj.j2+ok.j2+ol.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{oa.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{oa.j2,oc.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicjl = H.TwoBody.GetTBME_J(J1,J1,i,c,j,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdac = H.TwoBody.GetTBME_J(J2,J2,j,d,a,c);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           // J4 is tied to two single-particle pairs and to J1,J2 directly.
                           int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{od.j2,ol.j2},{2*J1,2*J2}} ) /2;
                           int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{od.j2,ol.j2},{2*J1,2*J2}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( J1, J2, J4, oa.j2/2., oi.j2/2., oc.j2/2. );
                             double sixj3 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, ok.j2/2., ob.j2/2., J3 );
                             double sixj4 = H.modelspace->GetSixJ( J1, J4, J2, od.j2/2., oj.j2/2., ol.j2/2. );
                             F21 += -phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vicjl * vjdac * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F21;
}



// Diagram F22 (as numbered by ADG)   complex conjugate diagram: F26
// mscheme expression: F22 = 1/2 sum_abcdijkl (v_abij v_ickl v_jdac v_klbd) / (eps^ij_ab eps^jkl_acb eps^kl_bd)
// Missing minus sign? Yes, dropped the minus sign from AMC
double HFMBPT::GetMP4_F22( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F22 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jklacb = GetDenom(H,{j,k,l},{a,c,b});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ijab * e_jklacb * e_klbd;
                   int phase_exponent = (ob.j2+oj.j2)/2;  // + J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,od.j2},{oa.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,od.j2},{oa.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vickl = H.TwoBody.GetTBME_J(J1,J1,i,c,k,l);
                       double vklbd = H.TwoBody.GetTBME_J(J1,J1,k,l,b,d);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjdac = H.TwoBody.GetTBME_J(J2,J2,j,d,a,c);
                         double ninej = H.modelspace->GetNineJ( oa.j2/2., ob.j2/2., J0, oc.j2/2., J1, oi.j2/2., J2, od.j2/2., oj.j2/2. );
                         double phase = ( (J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F22 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vickl * vjdac * vklbd / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F22;
}



// Diagram F23 (as numbered by ADG)   complex conjugate diagram: F29
// mscheme expression: F23 = 1/2 sum_abcdijkl (v_abik v_icjl v_jdab v_klcd) / (eps^ik_ab eps^jkl_abc eps^kl_cd)
// missing minus sign? yep. fixed it.
double HFMBPT::GetMP4_F23( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F23 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_jklabc = GetDenom(H,{j,k,l},{a,b,c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ikab * e_jklabc * e_klcd;
                   int phase_exponent = (oc.j2+ok.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2},{oj.j2,od.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2},{oj.j2,od.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{oj.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{oj.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     double vjdab = H.TwoBody.GetTBME_J(J0,J0,j,d,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vicjl = H.TwoBody.GetTBME_J(J1,J1,i,c,j,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vklcd = H.TwoBody.GetTBME_J(J2,J2,k,l,c,d);
                         double ninej = H.modelspace->GetNineJ( ok.j2/2., oi.j2/2., J0, ol.j2/2., J1, oj.j2/2., J2, oc.j2/2., od.j2/2. );
                         double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F23 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabik * vicjl * vjdab * vklcd / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F23;
}



// Diagram F24 (as numbered by ADG)   complex conjugate diagram: F30
// mscheme expression: F24 = -1/4 sum_abcdijkl (v_abij v_ickl v_jdab v_klcd) / (eps^ij_ab eps^jkl_abc eps^kl_cd)
// missing minus sign. fixed it.
double HFMBPT::GetMP4_F24( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F24 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     double a_sum = 0;
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             if ( oi.j2 != od.j2) continue;
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_jklabc = GetDenom(H,{j,k,l},{a,b,c});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ijab * e_jklabc * e_klcd;
                   int phase_exponent = (oc.j2+oj.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oj.j2,od.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oj.j2,od.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oi.j2,oc.j2},{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oi.j2,oc.j2},{ok.j2,ol.j2},{oc.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vjdab = H.TwoBody.GetTBME_J(J0,J0,j,d,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vickl = H.TwoBody.GetTBME_J(J1,J1,i,c,k,l);
                       double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F24 += 1./4 * phase * (2*J0+1) * (2*J1+1) / (od.j2+1) * vabij * vickl * vjdab * vklcd / denom;
                       a_sum += 1./4 * phase * (2*J0+1) * (2*J1+1) / (od.j2+1) * vabij * vickl * vjdab * vklcd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
     std::cout << "    " << a << "  " << a_sum << "  " << oa.j2+1 << std::endl;
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F24;
}

double HFMBPT::GetMP4_F25( const Operator& H) const
{
   double F25 = GetMP4_F21(H); // conjugate diagram
   return F25;
}
double HFMBPT::GetMP4_F26( const Operator& H) const
{
   double F26 = GetMP4_F22(H); // conjugate diagram
   return F26;
}




// Diagram F27 (as numbered by ADG)
// mscheme expression: F27 = 1/4 sum_abcdeijk (v_abjk v_cdai v_iecd v_jkbe) / (eps^jk_ab eps^ijk_cdb eps^jk_be)
// missing minus sign?  yep. fixed it.
double HFMBPT::GetMP4_F27( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F27 =0;
//   #pragma omp parallel for
             for ( auto i: H.modelspace->holes )
             {
               double a_sum=0;
               Orbit& oi = H.modelspace->GetOrbit(i);
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             if ( oe.j2 != oa.j2) continue;
//             for ( auto i: H.modelspace->holes )
//             {
//               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_jkab = GetDenom(H,{j,k},{a,b});
                   double e_ijkcdb = GetDenom(H,{i,j,k},{c,d,b});
                   double e_jkbe = GetDenom(H,{j,k},{b,e});
                   double denom = e_jkab * e_ijkcdb * e_jkbe;
                   int phase_exponent = (ob.j2+oi.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oj.j2,ok.j2},{ob.j2,oe.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oj.j2,ok.j2},{ob.j2,oe.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,oi.j2},{oi.j2,oe.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,oi.j2},{oi.j2,oe.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabjk = H.TwoBody.GetTBME_J(J0,J0,a,b,j,k);
                     double vjkbe = H.TwoBody.GetTBME_J(J0,J0,j,k,b,e);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdai = H.TwoBody.GetTBME_J(J1,J1,c,d,a,i);
                       double viecd = H.TwoBody.GetTBME_J(J1,J1,i,e,c,d);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F27 -= 1./4 * phase * (2*J0+1) * (2*J1+1) / (oa.j2+1) * vabjk * vcdai * viecd * vjkbe / denom;
                       a_sum -= 1./4 * phase * (2*J0+1) * (2*J1+1) / (oa.j2+1) * vabjk * vcdai * viecd * vjkbe / denom;
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
     std::cout << "    " << i << "  " << a_sum << "  " << oi.j2+1 << std::endl;
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F27;
}



// Diagram F28 (as numbered by ADG)
// mscheme expression: F28 = 1/2 sum_abcdeijk (v_abij v_cdak v_iecd v_jkbe) / (eps^ij_ab eps^ijk_cdb eps^jk_be)
// agrees.
double HFMBPT::GetMP4_F28( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F28 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijkcdb = GetDenom(H,{i,j,k},{c,d,b});
                   double e_jkbe = GetDenom(H,{j,k},{b,e});
                   double denom = e_ijab * e_ijkcdb * e_jkbe;
                   int phase_exponent = (ob.j2+oj.j2)/2;  // + J0+J1+J2, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,ok.j2},{oi.j2,oe.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,ok.j2},{oi.j2,oe.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oj.j2,ok.j2},{ob.j2,oe.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oj.j2,ok.j2},{ob.j2,oe.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdak = H.TwoBody.GetTBME_J(J1,J1,c,d,a,k);
                       double viecd = H.TwoBody.GetTBME_J(J1,J1,i,e,c,d);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vjkbe = H.TwoBody.GetTBME_J(J2,J2,j,k,b,e);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oe.j2/2., J1, oi.j2/2., J2, ok.j2/2., oj.j2/2. );
                         double phase = ( (J0+J1+J2+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                         F28 += 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabij * vcdak * viecd * vjkbe / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F28;
}

double HFMBPT::GetMP4_F29( const Operator& H) const
{
   double F29 = GetMP4_F23(H); // conjugate diagram
   return F29;
}
double HFMBPT::GetMP4_F30( const Operator& H) const
{
   double F30 = GetMP4_F24(H); // conjugate diagram
   return F30;
}


// Diagram F31 (as numbered by ADG)
// mscheme expression: F31 = 1/2 sum_abcdeijk (v_abjk v_cdai v_iebc v_jkde) / (eps^jk_ab eps^ijk_bcd eps^jk_de)
// missing minus sign. fixed it.
double HFMBPT::GetMP4_F31( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F31 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_jkab = GetDenom(H,{j,k},{a,b});
                   double e_ijkbcd = GetDenom(H,{i,j,k},{b,c,d});
                   double e_jkde = GetDenom(H,{j,k},{d,e});
                   double denom = e_jkab * e_ijkbcd * e_jkde;
                   int phase_exponent = (ob.j2+oc.j2+od.j2+oi.j2)/2;  // + J0, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,oi.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,oi.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oi.j2,oe.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oi.j2,oe.j2},{ob.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabjk = H.TwoBody.GetTBME_J(J0,J0,a,b,j,k);
                     double vjkde = H.TwoBody.GetTBME_J(J0,J0,j,k,d,e);
                     double phase = ( (J0+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdai = H.TwoBody.GetTBME_J(J1,J1,c,d,a,i);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double viebc = H.TwoBody.GetTBME_J(J2,J2,i,e,b,c);
                         double ninej = H.modelspace->GetNineJ( ob.j2/2., oa.j2/2., J0, oc.j2/2., J1, od.j2/2., J2, oi.j2/2., oe.j2/2. );
                         F31 -= 1./2 * phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * ninej * vabjk * vcdai * viebc * vjkde / denom;
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F31;
}



// Diagram F32 (as numbered by ADG)
// mscheme expression: F32 = sum_abcdeijk (v_abij v_cdak v_iebc v_jkde) / (eps^ij_ab eps^ijk_bcd eps^jk_de)
// agrees.
double HFMBPT::GetMP4_F32( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F32 =0;
   std::vector<size_t> part_vec;
   for ( auto a : H.modelspace->particles) part_vec.push_back(a);
   #pragma omp parallel for  collapse(3) reduction(+:F32)
   for ( auto a : part_vec )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: part_vec )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: part_vec )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto e: H.modelspace->particles )
           {
             Orbit& oe = H.modelspace->GetOrbit(e);
             for ( auto i: H.modelspace->holes )
             {
               Orbit& oi = H.modelspace->GetOrbit(i);
               for ( auto j: H.modelspace->holes )
               {
                 Orbit& oj = H.modelspace->GetOrbit(j);
                 for ( auto k: H.modelspace->holes )
                 {
                   Orbit& ok = H.modelspace->GetOrbit(k);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijkbcd = GetDenom(H,{i,j,k},{b,c,d});
                   double e_jkde = GetDenom(H,{j,k},{d,e});
                   double denom = e_ijab * e_ijkbcd * e_jkde;
                   int phase_exponent = (oa.j2+ob.j2+oe.j2+ok.j2)/2;  // + J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oa.j2,ok.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oa.j2,ok.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oi.j2,oe.j2},{ob.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oi.j2,oe.j2},{ob.j2,oc.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{oj.j2,ok.j2},{od.j2,oe.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdak = H.TwoBody.GetTBME_J(J1,J1,c,d,a,k);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double viebc = H.TwoBody.GetTBME_J(J2,J2,i,e,b,c);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vjkde = H.TwoBody.GetTBME_J(J3,J3,j,k,d,e);
                           double phase = ( (J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           // J4 is tied to two single-particle pairs and to J1,J3 directly.
                           int J4_min = AngMom::Jmin( {{oi.j2,ob.j2},{oc.j2,oe.j2},{2*J3,2*J1}} ) /2;
                           int J4_max = AngMom::Jmax( {{oi.j2,ob.j2},{oc.j2,oe.j2},{2*J3,2*J1}} ) /2;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( oi.j2/2., ob.j2/2., J4, oa.j2/2., oj.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( J3, J1, J4, oa.j2/2., oj.j2/2., ok.j2/2. );
                             double sixj3 = H.modelspace->GetSixJ( oc.j2/2., oe.j2/2., J4, oi.j2/2., ob.j2/2., J2 );
                             double sixj4 = H.modelspace->GetSixJ( J3, J4, J1, oc.j2/2., od.j2/2., oe.j2/2. );
                             F32 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabij * vcdak * viebc * vjkde / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for k
               }// for j
             }// for i
           }//for e
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F32;
}

// Diagram F33 (as numbered by ADG)
// mscheme expression: F33 = -1/4 sum_abcdijkl (v_abik v_cdjl v_ijcd v_klab) / (eps^ik_ab eps^ijkl_cdab eps^kl_ab)
// minus sign error. Fixed.
double HFMBPT::GetMP4_F33( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F33 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   if ( ol.j2 != oi.j2) continue;
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_ijklcdab = GetDenom(H,{i,j,k,l},{c,d,a,b});
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double denom = e_ikab * e_ijklcdab * e_klab;
                   int phase_exponent = (oj.j2+ok.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2},{ok.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oj.j2,ol.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oj.j2,ol.j2},{oi.j2,oj.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     double vklab = H.TwoBody.GetTBME_J(J0,J0,k,l,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdjl = H.TwoBody.GetTBME_J(J1,J1,c,d,j,l);
                       double vijcd = H.TwoBody.GetTBME_J(J1,J1,i,j,c,d);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F33 += 1./4 * phase * (2*J0+1) * (2*J1+1) / (oi.j2+1) * vabik * vcdjl * vijcd * vklab / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F33;
}


// Diagram F34 (as numbered by ADG)
// mscheme expression: F34 = 1/16 sum_abcdijkl (v_abij v_cdkl v_ijcd v_klab) / (eps^ij_ab eps^ijkl_cdab eps^kl_ab)
// agrees.
double HFMBPT::GetMP4_F34( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F34 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijklcdab = GetDenom(H,{i,j,k,l},{c,d,a,b});
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double denom = e_ijab * e_ijklcdab * e_klab;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2},{ok.j2,ol.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oc.j2,od.j2},{ok.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vcdkl = H.TwoBody.GetTBME_J(J0,J0,c,d,k,l);
                     double vijcd = H.TwoBody.GetTBME_J(J0,J0,i,j,c,d);
                     double vklab = H.TwoBody.GetTBME_J(J0,J0,k,l,a,b);
                     F34 += 1./16 * (2*J0+1) * vabij * vcdkl * vijcd * vklab / denom;
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F34;
}


// Diagram F35 (as numbered by ADG)
// mscheme expression: F35 = -1/4 sum_abcdijkl (v_abkl v_cdij v_ijac v_klbd) / (eps^kl_ab eps^ijkl_acbd eps^kl_bd)
// minus sign error. Fixed.
double HFMBPT::GetMP4_F35( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F35 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     double a_sum = 0;
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           if ( od.j2 != oa.j2) continue;
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double e_ijklacbd = GetDenom(H,{i,j,k,l},{a,c,b,d});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_klab * e_ijklacbd * e_klbd;
                   int phase_exponent = (ob.j2+oc.j2)/2;  // + J0+J1, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabkl = H.TwoBody.GetTBME_J(J0,J0,a,b,k,l);
                     double vklbd = H.TwoBody.GetTBME_J(J0,J0,k,l,b,d);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdij = H.TwoBody.GetTBME_J(J1,J1,c,d,i,j);
                       double vijac = H.TwoBody.GetTBME_J(J1,J1,i,j,a,c);
                       double phase = ( (J0+J1+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                       F35 += 1./4 * phase * (2*J0+1) * (2*J1+1) / (oa.j2+1) * vabkl * vcdij * vijac * vklbd / denom;
                       a_sum += 1./4 * phase * (2*J0+1) * (2*J1+1) / (oa.j2+1) * vabkl * vcdij * vijac * vklbd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
     std::cout << "    " << a << "  " << a_sum << "  " << oa.j2+1 << std::endl;
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F35;
}



// Diagram F36 (as numbered by ADG)
// mscheme expression: F36 = sum_abcdijkl (v_abik v_cdjl v_ijac v_klbd) / (eps^ik_ab eps^ijkl_acbd eps^kl_bd)
// agrees.
double HFMBPT::GetMP4_F36( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F36 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_ijklacbd = GetDenom(H,{i,j,k,l},{a,c,b,d});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ikab * e_ijklacbd * e_klbd;
                   int phase_exponent = (ob.j2+oc.j2+oj.j2+ok.j2)/2;  // + J0+J1+J2+J3, added inside the loops below
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oj.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oj.j2,ol.j2}} ) /2;
                   int J2_min = AngMom::Jmin( {{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J2_max = AngMom::Jmax( {{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J3_min = AngMom::Jmin( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J3_max = AngMom::Jmax( {{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   // J4 does not label any two-body matrix element; its range comes purely from
                   // the (single-particle) triangle conditions of the four 6j symbols it appears in.
                   int J4_min = AngMom::Jmin( {{ok.j2,ob.j2},{oj.j2,oc.j2},{od.j2,ol.j2}} ) /2;
                   int J4_max = AngMom::Jmax( {{ok.j2,ob.j2},{oj.j2,oc.j2},{od.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdjl = H.TwoBody.GetTBME_J(J1,J1,c,d,j,l);
                       for (int J2=J2_min; J2<=J2_max; J2++)
                       {
                         double vijac = H.TwoBody.GetTBME_J(J2,J2,i,j,a,c);
                         for (int J3=J3_min; J3<=J3_max; J3++)
                         {
                           double vklbd = H.TwoBody.GetTBME_J(J3,J3,k,l,b,d);
                           double phase = ( (J0+J1+J2+J3+phase_exponent)%2==0 ) ? 1.0 : -1.0;
                           for (int J4=J4_min; J4<=J4_max; J4++)
                           {
                             double sixj1 = H.modelspace->GetSixJ( ok.j2/2., ob.j2/2., J4, oa.j2/2., oi.j2/2., J0 );
                             double sixj2 = H.modelspace->GetSixJ( oj.j2/2., oc.j2/2., J4, oa.j2/2., oi.j2/2., J2 );
                             double sixj3 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, ok.j2/2., ob.j2/2., J3 );
                             double sixj4 = H.modelspace->GetSixJ( od.j2/2., ol.j2/2., J4, oj.j2/2., oc.j2/2., J1 );
                             F36 += phase * (2*J0+1) * (2*J1+1) * (2*J2+1) * (2*J3+1) * (2*J4+1)
                                    * sixj1 * sixj2 * sixj3 * sixj4 * vabik * vcdjl * vijac * vklbd / denom;
                           }// for J4
                         }// for J3
                       }// for J2
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F36;
}



// Diagram F37 (as numbered by ADG)
// mscheme expression: F37 = -1/4 sum_abcdijkl (v_abij v_cdkl v_ijac v_klbd) / (eps^ij_ab eps^ijkl_acbd eps^kl_bd)
// agrees.
double HFMBPT::GetMP4_F37( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F37 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         if ( oc.j2 != ob.j2) continue;
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ijab = GetDenom(H,{i,j},{a,b});
                   double e_ijklacbd = GetDenom(H,{i,j,k,l},{a,c,b,d});
                   double e_klbd = GetDenom(H,{k,l},{b,d});
                   double denom = e_ijab * e_ijklacbd * e_klbd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,oj.j2},{oa.j2,oc.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{ok.j2,ol.j2},{ob.j2,od.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabij = H.TwoBody.GetTBME_J(J0,J0,a,b,i,j);
                     double vijac = H.TwoBody.GetTBME_J(J0,J0,i,j,a,c);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdkl = H.TwoBody.GetTBME_J(J1,J1,c,d,k,l);
                       double vklbd = H.TwoBody.GetTBME_J(J1,J1,k,l,b,d);
                       F37 += -1./4 * (2*J0+1) * (2*J1+1) / (ob.j2+1) * vabij * vcdkl * vijac * vklbd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F37;
}



// Diagram F38 (as numbered by ADG)
// mscheme expression: F38 = 1/16 sum_abcdijkl (v_abkl v_cdij v_ijab v_klcd) / (eps^kl_ab eps^ijkl_abcd eps^kl_cd)
// agrees.
double HFMBPT::GetMP4_F38( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F38 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_klab = GetDenom(H,{k,l},{a,b});
                   double e_ijklabcd = GetDenom(H,{i,j,k,l},{a,b,c,d});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_klab * e_ijklabcd * e_klcd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{ok.j2,ol.j2},{oc.j2,od.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{ok.j2,ol.j2},{oc.j2,od.j2},{oi.j2,oj.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabkl = H.TwoBody.GetTBME_J(J0,J0,a,b,k,l);
                     double vcdij = H.TwoBody.GetTBME_J(J0,J0,c,d,i,j);
                     double vijab = H.TwoBody.GetTBME_J(J0,J0,i,j,a,b);
                     double vklcd = H.TwoBody.GetTBME_J(J0,J0,k,l,c,d);
                     F38 += 1./16 * (2*J0+1) * vabkl * vcdij * vijab * vklcd / denom;
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F38;
}



// Diagram F39 (as numbered by ADG)
// mscheme expression: F39 = -1/4 sum_abcdijkl (v_abik v_cdjl v_ijab v_klcd) / (eps^ik_ab eps^ijkl_abcd eps^kl_cd)
// agrees.
double HFMBPT::GetMP4_F39( const Operator& H) const
{
   double t_start = omp_get_wtime();

   double F39 =0;
//   #pragma omp parallel for
   for ( auto a : H.modelspace->particles )
   {
     double a_sum =0;
     Orbit& oa = H.modelspace->GetOrbit(a);
     for ( auto b: H.modelspace->particles )
     {
       Orbit& ob = H.modelspace->GetOrbit(b);
       for ( auto c: H.modelspace->particles )
       {
         Orbit& oc = H.modelspace->GetOrbit(c);
         for ( auto d: H.modelspace->particles )
         {
           Orbit& od = H.modelspace->GetOrbit(d);
           for ( auto i: H.modelspace->holes )
           {
             Orbit& oi = H.modelspace->GetOrbit(i);
             for ( auto j: H.modelspace->holes )
             {
               Orbit& oj = H.modelspace->GetOrbit(j);
               for ( auto k: H.modelspace->holes )
               {
                 Orbit& ok = H.modelspace->GetOrbit(k);
                 if ( ok.j2 != oj.j2) continue;
                 for ( auto l: H.modelspace->holes )
                 {
                   Orbit& ol = H.modelspace->GetOrbit(l);
                   double e_ikab = GetDenom(H,{i,k},{a,b});
                   double e_ijklabcd = GetDenom(H,{i,j,k,l},{a,b,c,d});
                   double e_klcd = GetDenom(H,{k,l},{c,d});
                   double denom = e_ikab * e_ijklabcd * e_klcd;
                   int J0_min = AngMom::Jmin( {{oa.j2,ob.j2},{oi.j2,ok.j2},{oi.j2,oj.j2}} ) /2;
                   int J0_max = AngMom::Jmax( {{oa.j2,ob.j2},{oi.j2,ok.j2},{oi.j2,oj.j2}} ) /2;
                   int J1_min = AngMom::Jmin( {{oc.j2,od.j2},{oj.j2,ol.j2},{ok.j2,ol.j2}} ) /2;
                   int J1_max = AngMom::Jmax( {{oc.j2,od.j2},{oj.j2,ol.j2},{ok.j2,ol.j2}} ) /2;
                   for (int J0=J0_min; J0<=J0_max; J0++)
                   {
                     double vabik = H.TwoBody.GetTBME_J(J0,J0,a,b,i,k);
                     double vijab = H.TwoBody.GetTBME_J(J0,J0,i,j,a,b);
                     for (int J1=J1_min; J1<=J1_max; J1++)
                     {
                       double vcdjl = H.TwoBody.GetTBME_J(J1,J1,c,d,j,l);
                       double vklcd = H.TwoBody.GetTBME_J(J1,J1,k,l,c,d);
                       F39 += -1./4 * (2*J0+1) * (2*J1+1) / (oj.j2+1) * vabik * vcdjl * vijab * vklcd / denom;
                       a_sum += -1./4 * (2*J0+1) * (2*J1+1) / (oj.j2+1) * vabik * vcdjl * vijab * vklcd / denom;
                     }// for J1
                   }// for J0
                 }// for l
               }// for k
             }// for j
           }// for i
         }//for d
       }//for c
     }//for b
     std::cout << "    " << a << "  " << a_sum << "  " << oa.j2+1 << std::endl;
   }//for a

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return F39;
}



