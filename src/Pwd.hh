#ifndef PWD_h
#define PWD_h 1

#include <iostream>
#include <unordered_map>
#include <gsl/gsl_integration.h>
#include <unordered_map>
#include <map>
#include <functional>
#include <cstdint>

/// Class containing the partial wave decomposition of any scalar two-body potential
/// as described in https://doi.org/10.1016/0375-9474(71)90279-X. The advantages
/// is that any 2B scalar potential can be splitted into 
/// \f[
///        V_(q) = \sum_i [\Psi^{(0)}_i(q)\tau_1\tau_2+\Psi^{(1)}_i(q)]\Omega_i
/// \f]
/// where q is the momentum transfer, psi is a scalar function and Omega is 
/// one of the six possibility for the spin operator. It is therefore very easy
/// to implement any scalar 2B operator using this class to generatre the relative
/// frame TBME. Example on how this is use can be found for the 0vbb operators in
/// M0nu.cc
class PWD
{
  private:
    /// Integration mesh. By making these attributes, we can easily share them
    /// when precomputing the relative space integrals over the input and output
    /// momenta
    int momentum_mesh_size;
    gsl_integration_glfixed_table *momentum_mesh;
    int angular_mesh_size;
    gsl_integration_glfixed_table *angular_mesh;
    double max_momentum;

    /// The six possible spin operators and the function associated with them
    /// they are initiated by the function setPotential
    const char *types[6] = {"central", "spin-spin", "spin-orbit", "sigma_L", "tensor", "sigma_k"};
    std::function<double(double, double, double)> potential_central;
    std::function<double(double, double, double)> potential_spin_spin;
    std::function<double(double, double, double)> potential_spin_orbit;
    std::function<double(double, double, double)> potential_sigma_L;
    std::function<double(double, double, double)> potential_tensor;
    std::function<double(double, double, double)> potential_sigma_k;
    std::map<std::string, std::function<double(double, double, double)>> potential_list = {{"central", potential_central},
                                                                                           {"spin-spin", potential_spin_spin},
                                                                                           {"spin-orbit", potential_spin_orbit},
                                                                                           {"sigma_L", potential_sigma_L},
                                                                                           {"tensor", potential_tensor},
                                                                                           {"sigma_k", potential_sigma_k}};
    bool central_on = false;
    bool spin_spin_on = false;
    bool spin_orbit_on = false;
    bool sigma_L_on = false;
    bool tensor_on = false;
    bool sigma_k_on = false;
    std::map<std::string, bool> potential_bool = {{"central", central_on},
                                                  {"spin-spin", spin_spin_on},
                                                  {"spin-orbit", spin_orbit_on},
                                                  {"sigma_L", sigma_L_on},
                                                  {"tensor", tensor_on},
                                                  {"sigma_k", sigma_k_on}};
    /// List of precalculated angular integrals. This allows for significant
    /// speed up as the angular integrals are reused multiple times in the 
    /// momentum integrals.
    std::unordered_map<uint64_t, double> AList;
    /// Regulator function. Can be set to either local or nonlocal by setRegulator.
    std::function<double(double, double, double)> regulator;

  public:
    void initializeMomentumMesh(int size);
    void setMomentumMesh(gsl_integration_glfixed_table *t, int size);
    int getMomentumMeshSize();
    gsl_integration_glfixed_table *getMomentumMesh();
    void freeMomentumMesh();
     
    void initializeAngularMesh(int size);
    void setAngularMesh(gsl_integration_glfixed_table *momentum_mesh, int size);
    int getAngularMeshSize();
    gsl_integration_glfixed_table *getAngularMesh();
    void freeAngularMesh();

    void setMaxMomentum(double max);
    double getMaxMomentum();

    void setPotential(std::function<double(double, double, double)> potential_func,  std::string potential_type);

    void setRegulator(double regulator_cutoff, int regulator_power, std::string type);

    void calcA(int e2max, int lmax);
    int getAsize();
    void clearA();

    double getW(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J);

    PWD();
};

// Functions to precalculate and hash the angular integrals.
double AIntegrand(double p, double pp, int J, int l, std::function<double(double,double,double)> potential);
double A(double p, double pp, int J, int l, std::function<double(double, double, double)> potential, gsl_integration_glfixed_table *t, int z_size);
uint64_t AHash(int index_p, int index_pp, int J, int l, int type);
void AUnHash(uint64_t key, uint64_t &index_p, uint64_t &index_pp, uint64_t &J, uint64_t &l, uint64_t &type);
void PrecalculateA(int e2max, std::function<double(double, double, double)> potential, int lmax, int type, double max_momentum, gsl_integration_glfixed_table *t_momentum, gsl_integration_glfixed_table *t_z, int n_momentum_points, int n_z_points, std::unordered_map<uint64_t, double>& AList);
double GetA(int index_p, int index_pp, int J, int l, int type, std::unordered_map<uint64_t, double> &AList);

/// The 6 different decomposition depending on the spin part of the operator
double central_force_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double> &AList);
double spin_spin_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double> &AList);
double spin_orbit_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double> &AList);
double sigma_L_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J,  std::unordered_map<uint64_t, double> &AList);
double tensor_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double> &AList);
double sigma_k_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double> &AList);


double regulator_local(double q, double regulator_cutoff, int regulator_power);
double regulator_nonlocal(double p, double pp, double regulator_cutoff, int regulator_power);

#endif
