#include "Pwd.hh"
#include <gsl/gsl_integration.h>
#include <unordered_map>
#include <vector>
#include <map>
#include <functional>
#include<iostream>
#include <gsl/gsl_sf_legendre.h>
#include "PhysicalConstants.hh"

using namespace PhysConst;

PWD::PWD()
{
  for (int i = 0; i < 6; i++)
  {
    setPotential([](double p, double pp, double z) { return 0; }, types[i]);
    potential_bool[types[i]] =  false;
    regulator = [](double p, double pp, double z) { return 1.0; };
  }
}

void PWD::initializeMomentumMesh(int size)
{
  momentum_mesh_size = size;
  momentum_mesh = gsl_integration_glfixed_table_alloc(momentum_mesh_size);
}

void PWD::setMomentumMesh(gsl_integration_glfixed_table *t, int size)
{
  momentum_mesh_size = size;
  momentum_mesh = t;
}

int PWD::getMomentumMeshSize()
{
  return momentum_mesh_size;
}

gsl_integration_glfixed_table* PWD::getMomentumMesh()
{
  return momentum_mesh;
}

void PWD::freeMomentumMesh()
{
  gsl_integration_glfixed_table_free(momentum_mesh);
}

void PWD::initializeAngularMesh(int size)
{
  angular_mesh_size = size;
  angular_mesh = gsl_integration_glfixed_table_alloc(angular_mesh_size);
}

void PWD::setAngularMesh(gsl_integration_glfixed_table *t, int size)
{
  angular_mesh_size = size;
  angular_mesh = t; 
}

int PWD::getAngularMeshSize()
{
  return angular_mesh_size;
}

gsl_integration_glfixed_table *PWD::getAngularMesh()
{
  return angular_mesh;
}

void PWD::freeAngularMesh()
{
  gsl_integration_glfixed_table_free(angular_mesh);
}

void PWD::setMaxMomentum(double max)
{
  max_momentum = max;
}


double PWD::getMaxMomentum()
{
  return max_momentum;
}

void PWD::setRegulator(double regulator_cutoff, int regulator_power, std::string type)
{
  if (type == "local")
  {
    regulator = [regulator_cutoff, regulator_power](double p, double pp, double z)
                {
                  double q = sqrt(p*p+pp*pp-2*p*pp*z);
                  return regulator_local(q, regulator_cutoff, regulator_power);
                 };
  }
  else if (type == "nonlocal")
  {
    regulator = [regulator_cutoff, regulator_power](double p, double pp, double z)
    {
      return regulator_nonlocal(p, pp, regulator_cutoff, regulator_power);
    };
  }
}

void PWD::setPotential(std::function<double(double, double, double)> potential_func, std::string potential_type)
{
  potential_list[potential_type] = potential_func;
  potential_bool[potential_type] = true;

}



void PWD::calcA(int e2max, int lmax)
{
  std::map<std::string, int> types_lmax = {{"central", 0},
                                           {"spin-spin", 0},
                                           {"spin-orbit", 0},
                                           {"sigma_L", 2},
                                           {"tensor", 1},
                                           {"sigma_k", 1}
                                          };
  for (int i=0; i<6; i++)
  {
    if (potential_bool[types[i]] == true)
    { 
      std::cout<<types[i]<<std::endl;
      int l = std::min(types_lmax[types[i]], lmax);
     
      auto pot = [&](double p, double pp, double z)
      {
        return potential_list[types[i]](p,pp,z)*regulator(p,pp,z);
      };
      PrecalculateA(e2max, pot, l, i, max_momentum, momentum_mesh, angular_mesh, momentum_mesh_size, angular_mesh_size, AList);
    } 
  }
}

void PWD::clearA()
{
  AList.clear();
}

int PWD::getAsize()
{
  return AList.size();
}


double PWD::getW(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J)
{
  double W = 0;
  std::map<std::string, std::function<double(double,double,int,int,int,int,int,int,std::unordered_map<uint64_t,double>&)>> Wtypes  = {{"central", central_force_decomposition},
                                                                                                                                      {"spin-spin", spin_spin_decomposition},
                                                                                                                                      {"spin-orbit", spin_orbit_decomposition},
                                                                                                                                      {"sigma_L", sigma_L_decomposition},
                                                                                                                                      {"tensor", tensor_decomposition},
                                                                                                                                      {"sigma_k", sigma_k_decomposition}
                                                                                                                                     };
  for (int i =  0; i<6; i++)
  {
    if (potential_bool[types[i]]  == true)
    {
      double Wval = Wtypes[types[i]](p, pp, index_p, index_pp, S, L, Lp, J, AList);
      W += Wval;
    }
  }
  return W;
}


double AIntegrand(double p, double pp, double z, int J, int l, std::function<double(double, double, double)> potential)
{
  return potential(p, pp, z) * pow(z, l) * gsl_sf_legendre_Pl(J, z);
}

double A(double p, double pp, int J, int l, std::function<double(double, double, double)> potential, gsl_integration_glfixed_table *t, int n_z_points)
{
  double Aval = 0;
  for (int i = 0; i < n_z_points; i++)
  {
    double zi, wi;
    gsl_integration_glfixed_point(-1, 1, i, &zi, &wi, t);
    double Aint = AIntegrand(p, pp, zi, J, l, potential);
    Aval += wi * Aint;
  }
  return PI * Aval;
}

uint64_t AHash(int index_p, int index_pp, int J, int l, int type)
{
  return (((uint64_t)(index_p)) << 40) + (((uint64_t)(index_pp)) << 30) + (((uint64_t)(J)) << 20) + ((uint64_t)(l) << 10) + ((uint64_t)(type));
}

void AUnHash(uint64_t key, uint64_t &index_p, uint64_t &index_pp, uint64_t &J, uint64_t &l, uint64_t &type)
{
  index_p = (key >> 40) & 0x3FFL;
  index_pp = (key >> 30) & 0x3FFL;
  J = (key >> 20) & 0x3FFL;
  l = (key >> 10) & 0x3FFL;
  type = (key) &0x3FFL;
}

void PrecalculateA(int e2max, std::function<double(double, double, double)> potential, int lmax, int type, double max_momentum, gsl_integration_glfixed_table *t_momentum, gsl_integration_glfixed_table *t_z, int n_momentum_points, int n_z_points, std::unordered_map<uint64_t, double>& AList)
{
  int Jmax = e2max + 1; // maximal coupling of 2 l=emax with S=1
  std::vector<uint64_t> KEYS;
  for (int J = 0; J <= Jmax; J++)
  {
    for (int index_p = 0; index_p < n_momentum_points; index_p++)
    {
      for (int index_pp = index_p; index_pp < n_momentum_points; index_pp++)
      {
        for (int l = 0; l <= lmax; l++)
        {
          uint64_t key = AHash(index_p, index_pp, J, l, type);
          KEYS.push_back(key);
          AList[key] = 0.0; // Everything needs to be in this loop to avoid re-hash in parralel
        }
      }
    }
  }
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t n = 0; n < KEYS.size(); n++)
  {
    uint64_t key = KEYS[n];
    uint64_t index_p, index_pp, J, l,  type;
    AUnHash(key, index_p, index_pp, J, l, type);
    double p, pp, wi, wj;
    gsl_integration_glfixed_point(0, max_momentum, index_p, &p, &wi, t_momentum);
    gsl_integration_glfixed_point(0, max_momentum, index_pp, &pp, &wj, t_momentum);
    double aval = A(p, pp, J, l, potential, t_z, n_z_points);
    AList[key] = aval;
  }
}

double GetA(int index_p, int index_pp, int J, int l, int type, std::unordered_map<uint64_t, double>& AList)
{
  double A;
  if (index_p > index_pp)
    std::swap(index_p, index_pp);
  uint64_t key = AHash(index_p, index_pp, J, l, type);
  auto it = AList.find(key);
  if (it != AList.end()) // return what we've found
  {
    A = it->second;
  }
  else // if we didn't find it, calculate it and add it to the list!
  {
    std::cout << "DANGER!!!!!!!  Updating AList inside a parellel loop breaks thread safety!" << std::endl;
    std::cout << "   I shouldn't be here in GetA(" << index_p << " , " << index_pp << " , " << J << " , " << l << "):   key = " << key << std::endl;
//    printf("DANGER!!!!!!!  Updating AList inside a parellel loop breaks thread safety!\n");
//    printf("   I shouldn't be here in GetA(%d, %d, %d, %d):   key =%lx", index_p, index_pp, J, l, key);
    exit(EXIT_FAILURE);
  }
  return A;
}

double central_force_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
{
  int type = 0;
  double W = 0;
  if ((S == 0 and L == J and Lp == J)
    or (S == 1 and (L == Lp) and ((L == J) or (L == J - 1) or (L == J + 1)))) // Singlet state and triplet state are the same
    {
      W = 2 * GetA(index_p, index_pp, L, 0, type, AList);
    }
  return W;
}

double spin_spin_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J,  std::unordered_map<uint64_t, double>& AList)
{
  int type = 1;
  double W = 0;
  if (S == 0 and L == J and Lp == J) // Singlet state
  {
    W = -6 * GetA(index_p, index_pp, J, 0, type, AList);
  }
  else if (S == 1) // Triplet state
  {
    // The fucntions are the same if L=Lp regardless of the value of L and is 0 otherwise
    if ((L == Lp) and (((L == J) or (L == J - 1) or (L == J + 1))))
    {
      W = 2 * GetA(index_p, index_pp, L, 0, type, AList);
    }
  }
  return W;
}

double spin_orbit_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
{
  double W = 0;
  int type = 2;
  if (S == 1)
  {
    if (L == Lp)
    {
      if (L == J)
      {
        return 2 * p * pp / (2 * J + 1) * (GetA(index_p, index_pp, J + 1, 0, type, AList) - GetA(index_p, index_pp, J - 1, 0, type, AList));
      }
      else if (L == J - 1)
      {
        return 2 * p * pp * (J - 1) / (2 * J + 1) * (GetA(index_p, index_pp, J - 2, 0, type, AList) - GetA(index_p, index_pp, J, 0, type, AList));
      }
      else if (L == J + 1)
      {
        return 2 * p * pp * (J + 2) / (2 * J + 1) * (GetA(index_p, index_pp, J + 2, 0, type,  AList) - GetA(index_p, index_pp, J, 0, type, AList));
      }
    }
  }
  return W;
}

double sigma_L_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
{
  int type = 3;
  double W = 0;
  if (S==0 and  L==J and Lp==J) //Singlet state
  {
    return 2*p*p*pp*pp*(GetA(index_p,index_pp,J,2,type,AList)-GetA(index_p,index_pp,J,0,type,AList));
  }
  else if(S==1) //Triplet state
  {
    if (L==Lp)
    {
      if (L==J)
      {
        W = 2*p*p*pp*pp*(-GetA(index_p,index_pp,J,0,3,AList)+(J-1)/(2*J+1)*GetA(index_p,index_pp,J+1,1,type,AList)+(J-1)/(2*J+1)*GetA(index_p,index_pp,J-1,1,type,AList));
      }
      else if (L==J-1)
      {
        W = 2*p*p*pp*pp*((2*J-1)/(2*J+1)*GetA(index_p,index_pp,J-1,0,type,AList)-2/(2*J+1)*GetA(index_p,index_pp,J,1,type,AList)-GetA(index_p,index_pp,J-1,1,type,AList));
      }
      else if (L==J+1)
      {
        W = 2*p*p*pp*pp*((2*J+3)/(2*J+1)*GetA(index_p,index_pp,J+1,0,type,AList)-2/(2*J+1)*GetA(index_p,index_pp,J,1,type,AList)-GetA(index_p,index_pp,J+1,1,type,AList));
      }
    }
    else if ((Lp==J-1 and L==J+1) or (Lp==J+1 and L==J-1))
    {
      W = 4*p*p*pp*pp*sqrt(J*(J+1))/(2*J+1)/(2*J+1)*(GetA(index_p,index_pp,J+1,0,type,AList)-GetA(index_p,index_pp,J-1,0,type,AList));      
    }
  }
  return W;
}

double tensor_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
{
  int type = 4;
  double W = 0;
  if (S == 0 and L == J and Lp == J) // Singlet state
  {
    W = 2 * (-(p * p + pp * pp) * GetA(index_p, index_pp, J, 0, type, AList) + 2 * p * pp * GetA(index_p, index_pp, J, 1, type, AList));
  }
  else if (S == 1) // Triple state
  {
    if (L == Lp)
    {
      if (L == J)
      {
        W = 2 * ( (pp * pp + p * p) * GetA(index_p, index_pp, J, 0, type, AList)
                   - 2 * p * pp / (2 * J + 1) * (
                                  J * GetA(index_p, index_pp, J + 1, 0, type, AList) 
                                  + (J + 1) * GetA(index_p, index_pp, J - 1, 0, type, AList)));
      }
      else if (L == J - 1)
      {
        W = 2 * ( (p * p + pp * pp) * GetA(index_p, index_pp, J - 1, 0, type, AList) - 2 * p * pp * GetA(index_p, index_pp, J, 0, type, AList)) / (2 * J + 1);

      }
      else if (L == J + 1)
      {
        W = 2 * (-(p * p + pp * pp) * GetA(index_p, index_pp, J + 1, 0, type, AList) + 2 * p * pp * GetA(index_p, index_pp, J, 0, type, AList)) / (2 * J + 1);
      }
    }
    else if (Lp == J + 1 and L == J - 1)
    {
      W = -4 * (sqrt(J * (J + 1)) / (2 * J + 1)) * (p * p * GetA(index_p, index_pp, J + 1, 0, type, AList) + pp * pp * GetA(index_p, index_pp, J - 1, 0, type, AList) - 2 * p * pp * GetA(index_p, index_pp, J, 0, type, AList));
    }
    else if (Lp == J - 1 and L == J + 1)
    {
      W = -4 * (sqrt(J * (J + 1)) / (2 * J + 1)) * (p * p * GetA(index_p, index_pp, J - 1, 0, type, AList) + pp * pp * GetA(index_p, index_pp, J + 1, 0, type, AList) - 2 * p * pp * GetA(index_p, index_pp, J, 0, type, AList));
    }
  }
  return W;
}

double sigma_k_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
{
  double W = 0;
  int type = 5;
  if (S == 0 and L == J and Lp == J) // Singlet state
  {
    W = 1 / 2 * (-(p * p + pp * pp) * GetA(index_p, index_pp, J, 0, type, AList) - 2 * p * pp * GetA(index_p, index_pp, J, 1, type, AList));
  }
  else if (S == 1) // Triplet state
  {
    if (L == Lp)
    {
      if (L == J)
      {
        W = 1 / 2 * ((pp * pp + p * p) * GetA(index_p, index_pp, J, 0, type, AList) + 2 * p * pp / (2 * J + 1) * (GetA(index_p, index_pp, J + 1, 0, type, AList) + (J + 1) * GetA(index_p, index_pp, J - 1, 0, type, AList)));
      }
      else if (L == J - 1)
      {
        W = 1 / 2 / (2 * J + 1) * ((p * p + pp * pp) * GetA(index_p, index_pp, J - 1, 0, type, AList) + 2 * p * pp * GetA(index_p, index_pp, J, 0, type, AList));
      }
      else if (L == J + 1)
      {
        W = 1 / 2 / (2 * J + 1) * (-(p * p + pp * pp) * GetA(index_p, index_pp, J + 1, 0, type, AList) - 2 * p * pp * GetA(index_p, index_pp, J, 0, type, AList));
      }
    }
    else if (Lp == J - 1 and L == J + 1)
    {
      W = sqrt(J * (J + 1)) / (2 * J + 1) * (p * p * GetA(index_p, index_pp, J + 1, 0, type, AList) + pp * pp * GetA(index_p, index_pp, J - 1, 0, type, AList) + 2 * p * pp * GetA(index_p, index_pp, J, 0, type, AList));
    }
    else if (Lp == J + 1 and L == J - 1)
    {
      W = sqrt(J * (J + 1)) / (2 * J + 1) * (p * p * GetA(index_p, index_pp, J - 1, 0, type, AList) + pp * pp * GetA(index_p, index_pp, J + 1, 0, type, AList) + 2 * p * pp * GetA(index_p, index_pp, J, 0, type, AList));
    }
  }
  return W;
}

double regulator_local(double q, double regulator_cutoff, int regulator_power)
{
  return exp(-pow(q * HBARC / regulator_cutoff, 2 * regulator_power));
}

double regulator_nonlocal(double p, double pp, double regulator_cutoff, int regulator_power)
{
  return exp(-pow(p * HBARC / regulator_cutoff, 2 * regulator_power)) * exp(-pow(pp * HBARC / regulator_cutoff, 2 * regulator_power));
}
