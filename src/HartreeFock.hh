#ifndef HartreeFock_h
#define HartreeFock_h

#include "Operator.hh"
#include <armadillo>

class HartreeFock
{
 public:
   arma::mat C;             // coefficients, 1st index is ho basis, 2nd = HF basis
   arma::mat rho;           // density matrix rho_ij
   arma::mat t;             // kinetic energy
   arma::mat Vab;           // 1 body potenial
   arma::mat H;             // 1 body Hamiltonian H_ab
   arma::mat Vmon;          // Monopole 2-body interaction
   arma::vec energies;      // vector of single particle energies
   arma::vec prev_energies; // SPE's from last iteration
   Operator& Hbare;          // Pointer to (input) bare Hamiltonian
   double ediff;              // Energy difference from last iteration
   double tolerance;          // tolerance for convergence
   double EHF;
   

 //  HartreeFock(Operator * hbare=NULL, Operator *hhf=NULL);
   HartreeFock(Operator&  hbare);
   void BuildMonopoleV();
   void Diagonalize();
   void Diagonalize2();
   void UpdateH();
   void UpdateDensityMatrix();
   bool CheckConvergence();
   void Solve();
   void PrintOrbits();
   void UpdateHFOrbits();
   void BuildKineticEnergy(float hw);
   void PrintRho(int );
   void CalcEHF();
   Operator TransformToHFBasis( Operator& OpIn);
//   void TransformToHFBasis(Operator *OpIn, Operator *OpHF);
//   void TransformToHFBasis(Operator& OpIn, Operator& OpHF);
   void ReorderCoefficients();

};

#endif

