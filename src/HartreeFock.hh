#ifndef HartreeFock_h
#define HartreeFock_h

#include "Operator.hh"
#include <armadillo>

class HartreeFock
{
 public:
   Operator& Hbare;         // Input bare Hamiltonian
   arma::mat C;             // coefficients, 1st index is ho basis, 2nd = HF basis
   arma::mat rho;           // density matrix rho_ij
   arma::mat t;             // kinetic energy
   arma::mat Vab;           // 1 body potenial
   arma::mat H;             // 1 body Hamiltonian H_ab
   arma::mat Vmon;          // Monopole 2-body interaction
   arma::vec energies;      // vector of single particle energies
   arma::vec prev_energies; // SPE's from last iteration
   double tolerance;        // tolerance for convergence
   double EHF;              // Hartree-Fock energy (Normal-ordered 0-body term)
   

   HartreeFock(Operator&  hbare);
   void BuildMonopoleV();
   void Diagonalize();
   void UpdateH();
   void UpdateDensityMatrix();
   bool CheckConvergence();
   void Solve();
   void PrintOrbits();
   void BuildKineticEnergy(float hw);
   void PrintRho(int );
   void CalcEHF();
   void ReorderCoefficients();
   Operator TransformToHFBasis( Operator& OpIn);

   Operator GetHbare(){return Hbare;};

};

#endif

