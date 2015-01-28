#ifndef HartreeFock_h
#define HartreeFock_h

#include "Operator.hh"
#include "Operator3.hh"
#include <armadillo>
#include <map>

template<class OPERATOR>
class HartreeFock
{
 public:
//   Operator3& Hbare;         // Input bare Hamiltonian
//   Operator& Hbare;         // Input bare Hamiltonian
   OPERATOR& Hbare;         // Input bare Hamiltonian
   arma::mat C;             // coefficients, 1st index is ho basis, 2nd = HF basis
   arma::mat rho;           // density matrix rho_ij
   arma::mat t;             // kinetic energy
   arma::mat Vab;           // 1 body piece of 2 body potential
   arma::mat V3ab;          // 1 body piece of 3 body potential
   arma::mat H;             // 1 body Hamiltonian H_ab
   arma::mat Vmon;          // Monopole 2-body interaction
   arma::mat Vmon_exch;          // Monopole 2-body interaction
   arma::vec energies;      // vector of single particle energies
   arma::vec prev_energies; // SPE's from last iteration
   double tolerance;        // tolerance for convergence
   double EHF;              // Hartree-Fock energy (Normal-ordered 0-body term)
   int iterations;          // iterations used in Solve()

   map<unsigned long long int, double> Vmon3; // Monopole 3-body interaction
   

   HartreeFock(OPERATOR&  hbare);
//   HartreeFock(Operator&  hbare);
//   HartreeFock(Operator3&  hbare);
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
//   Operator TransformToHFBasis( Operator& OpIn);
   OPERATOR TransformToHFBasis( OPERATOR& OpIn);

   void BuildMonopoleV3();
//   void BuildMonopoleV3(Operator& Hbare);
//   void BuildMonopoleV3(Operator3& Hbare);

//   Operator GetHbare(){return Hbare;};
   OPERATOR GetHbare(){return Hbare;};

};



#endif

