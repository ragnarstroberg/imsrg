///////////////////////////////////////////////////////////////////////////////////
//    HartreeFock.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef HartreeFock_h
#define HartreeFock_h 1


#include "ModelSpace.hh"
#include "Operator.hh"
#include "IMSRGProfiler.hh"
#include "Jacobi3BME.hh"
#include "ThreeBodyME.hh"
#include <armadillo>
#include <vector>
#include <array>
#include <deque>

class Jacobi3BME; // Forward declaration

class HartreeFock
{
// using PhysConst::HBARC;
// using PhysConst::M_NUCLEON;
// using PhysConst::SQRT2;

 public:
   Operator& Hbare;         ///< Input bare Hamiltonian
   ModelSpace * modelspace; ///< Model Space of the Hamiltonian
   ModelSpace * ms_for_output_3N; ///< Model Space for the transformed 3N matrix elements (we may want a smaller space for this).
   arma::mat C;             ///< transformation coefficients, 1st index is ho basis, 2nd = HF basis
   arma::mat rho;           ///< density matrix rho_ij
   arma::mat KE;            ///< kinetic energy
   arma::mat Vij;           ///< 1 body piece of 2 body potential
   arma::mat V3ij;          ///< 1 body piece of 3 body potential
   arma::mat F;             ///< Fock matrix
   std::array< std::array< arma::mat,2>,3> Vmon;          ///< Monopole 2-body interaction
   std::array< std::array< arma::mat,2>,3> Vmon_exch;          ///< Monopole 2-body interaction
   arma::uvec holeorbs;     ///< list of hole orbits for generating density matrix
   arma::rowvec hole_occ;   /// occupations of hole orbits
   arma::vec energies;      ///< vector of single particle energies
   arma::vec prev_energies; ///< SPE's from last iteration
   double tolerance;        ///< tolerance for convergence
   double EHF;              ///< Hartree-Fock energy (Normal-ordered 0-body term)
   double e1hf;             ///< One-body contribution to EHF
   double e2hf;             ///< Two-body contribution to EHF
   double e3hf;             ///< Three-body contribution to EHF
   int iterations;          ///< iterations used in Solve()
   std::vector<uint64_t> Vmon3_keys;
   std::vector< double> Vmon3;
   IMSRGProfiler profiler;  ///< Profiler for timing, etc.
   std::deque<double> convergence_ediff; ///< Save last few convergence checks for diagnostics
   std::deque<double> convergence_EHF; ///< Save last few convergence checks for diagnostics
   std::deque<arma::mat> DIIS_density_mats; ///< Save density matrix from past iterations for DIIS
   std::deque<arma::mat> DIIS_error_mats; ///< Save error from past iterations for DIIS

   Jacobi3BME* jacobi3bme;  ///< pointer to 3-body matrix elements in jacobi basis, if we want to use that
   bool use_jacobi_3body;
   bool freeze_occupations;
   bool discard_NO2B_from_3N; // This should generally be false for any production runs.

// Methods
   HartreeFock(Operator&  hbare); ///< Constructor
   void BuildMonopoleV();         ///< Only the monopole part of V is needed, so construct it.
   void SetUpMonopoleV3Keys();    ///< Figure out the needed monopole terms.
   void BuildMonopoleV3();        ///< Only the monopole part of V3 is needed.
   void BuildMonopoleV3Jacobi();  ///< Construct directly from the Jacobi 3-body matrix elements
   void Diagonalize();            ///< Diagonalize the Fock matrix
   void UpdateF();                ///< Update the Fock matrix with the new transformation coefficients C
   void UpdateDensityMatrix();    ///< Update the density matrix with the new coefficients C
   void UpdateDensityMatrix_DIIS();    ///< Update the density matrix using Direct Inversion in Iterative Subspace
   void FillLowestOrbits();       ///< Get new occupations based on the current single-particle energies
   void UpdateReference();        ///< If we got new occupations in FillLowestOrbits, then we should update the hole states in the reference.
   bool CheckConvergence();       ///< Compare the current energies with those from the previous iteration
   void Solve();                  ///< Diagonalize and UpdateF until convergence
   void CalcEHF();                ///< Evaluate the Hartree Fock energy
   void PrintEHF();               ///< Print out the Hartree Fock energy
   void ReorderCoefficients();    ///< Reorder the coefficients in C to eliminate phases etc.
   void SetJacobi3BME( Jacobi3BME* jac ) {jacobi3bme = jac; use_jacobi_3body=true;}; ///< Setter.
   void DiscardNO2Bfrom3N() {discard_NO2B_from_3N = true;};
   Operator TransformToHFBasis( Operator& OpIn); ///< Transform an operator from oscillator basis to HF basis
//   Operator GetNormalOrderedH();  ///< Return the Hamiltonian in the HF basis at the normal-ordered 2body level.
   Operator GetNormalOrderedH(int particle_rank=2);  ///< Return the Hamiltonian in the HF basis at the normal-ordered 2body level.
   Operator GetNormalOrderedH(arma::mat& Cin, int particle_rank=2);  ///< Return the Hamiltonian in the HF basis at the normal-ordered 2body level.
   Operator GetOmega();           ///< Return a generator of the Hartree Fock transformation
   Operator GetHbare(){return Hbare;}; ///< Getter function for Hbare
   void PrintSPE(); ///< Print out the single-particle energies
   void PrintSPEandWF(); ///< Print out the single-particle energies and wave functions
   void FreeVmon();               ///< Free up the memory used to store Vmon3.
   void GetRadialWF(index_t index, std::vector<double>& R, std::vector<double>& PSI); ///< Return the radial wave function of an orbit in the HF basis
   double GetRadialWF_r(index_t index, double R); ///< Return the radial wave function of an orbit in the HF basis
   double GetHFPotential( size_t i, double r, double rprime);
   double GetAverageHFPotential( double r, double rprime);
   void FreezeOccupations(){freeze_occupations = true;};
   void UnFreezeOccupations(){freeze_occupations = false;};
   uint64_t Vmon3JacobiHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ, uint64_t E12, uint64_t alpha, uint64_t Ncm, uint64_t Lcm);
   void Vmon3JacobiUnHash(uint64_t key, uint64_t& na, uint64_t& nb, uint64_t& nc, uint64_t& Jab, uint64_t& twoJ, uint64_t& E12, uint64_t& alpha, uint64_t& Ncm, uint64_t& Lcm);
   static uint64_t Vmon3Hash(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e, uint64_t f);
   static void Vmon3UnHash(uint64_t key, int& a, int& b, int& c, int& d, int& e, int& f);
//   ThreeBodyMEpn GetTransformed3B(  );
//   ThreeBodyMEpn GetTransformed3B( Operator& OpIn );
//   ThreeBodyME GetTransformed3B( Operator& OpIn );
   ThreeBodyME GetTransformed3B( Operator& OpIn, arma::mat& C3b );
//   ThreeBodyMEpn GetValence3B( int emax, int E3max );
//   ThreeBodyMEpn GetValence3B( Operator& OpIn, int emax, int E3max );
   ThreeBodyME GetValence3B( Operator& OpIn, int emax, int E3max );
//   double GetTransformed3bme( int Jab, int Jde, int J2, size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);
   double GetTransformed3bme( Operator& OpIn, int Jab, int Jde, int J2, size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);
//   double GetHF3bme( int Jab, int Jde, int J2, size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);
//   double GetHF3bme( int Jab, int Jde, int J2, int tab, int tde, int T2, size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);
   void SetModelspaceForOutput3N( ModelSpace& ms );

/*
   void GetSecondOrderRho();
   void SwitchToNaturalOrbitals();
*/

};



#endif

