#include "Helicity.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "PhysicalConstants.hh"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include "AngMom.hh"
#include <map>
#include <vector>
#include <sstream>
#include <functional>

using namespace PhysConst;

namespace Helicity
{
  int phase(int x)
  {
    return x%2==0 ? 1 : -1;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///                                                                                     ///
  ///  Partial Wave Decomposition of 2-body scalar potentials in momentum space           ///
  ///  for  . Implementation of the derivation by Erkelenz et al.                         ///
  ///  found in K. Erkelenz, R. Alzetta, and K. Holinde, Nucl. Phys. A 176, 413 (1971).   ///
  ///  Implementd by AB                                                                   ///
  ///                                                                                     ///
  ///////////////////////////////////////////////////////////////////////////////////////////

  /// Integrand for the angular integration of the partial wave decomposition
  /// presented above eq 4.14 in Erkelenz et al. as the abbreviation A.
  /// Here p and pp are momenta of the two states and z = cos(theta).
  /// The potential needs to be written as a function of pp,p,z
  double AIntegrand(double p, double pp, double z, int J, int l, std::function<double(double,double,double)> potential)
  {
    return potential(p,pp,z)*pow(z,l)*gsl_sf_legendre_Pl(J,z);
  }

  // /// Integral of the angular part presented above eq 4.14 in Erkelenz et al. as the abbreviation A..
  // double A(double p, double pp, int J, int l, std::function<double(double,double,double)> potential, int n_z_points)
  // { 
  //   double A = 0;
  //   gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(n_z_points);
  //   for (int i=0;i<n_z_points;i++)
  //   {
  //     double zi, wi;
  //     gsl_integration_glfixed_point(-1,1,i,&zi,&wi,t);
  //     A += wi*AIntegrand(p,pp,zi,J,l,potential);
  //   }
  //   gsl_integration_glfixed_table_free(t);
  //   return PI*A;
  // }

  
  // /// Helicity decomposition of a scalar potential with a central force as 
  // /// described in eq 4.20 of Erkelenz et al.
  // double central_force_decomposition(double p, double pp, int S, int L, int Lp, int J, std::function<double(double, double, double)> potential, int n_z_points)
  // {
  //   double W  = 0;
  //   if (S==0 and L==J and Lp==J) //Singlet state
  //   {
  //     W = 2*A(p,pp,J,0,potential,n_z_points);
  //   } 
  //   else if(S==1) //Triplet state
  //   {
  //     //The fucntions are the same if L=Lp regardless of the value of L and is 0 otherwise
  //     if ( (L == Lp) and (((L == J) or (L == J - 1) or (L == J + 1))))
  //     {
  //       if ((L==J) or (L==J-1) or (L==J+1))
  //       {
  //         W  = 2*A(p,pp,L,0,potential,n_z_points);
  //       }
  //     }
  //   }
  //   return W;
  // }

  // /// Helicity decomposition of a scalar potential with a spin-spin force as 
  // /// described in eq 4.21 of Erkelenz et al. The value is  the same as the central_force_decomposition
  // /// with the exception of a factor of -3 for the singlet state. This is handle with the variable Seval which is
  // /// -3 for S = 0 and 1 for S = 1. Not that W =0 otherwise so all case are considered.
  // double spin_spin_decomposition(double p, double pp, int S, int L, int Lp, int J,  std::function<double(double,double,double)> potential, int n_z_points)
  // {
  //   int Seval = Seval = 2 * S * (S + 1) - 3;
  //   double  W = Seval*central_force_decomposition(p, pp, S, L, Lp, J, potential, n_z_points);
  //   return W;
  // }

  // /// Helicity decomposition of a scalar potential with a spin-orbit force as 
  // /// described in eq 4.22 of Erkelenz et al. The singlet state gives 0 and state with
  // /// L != Lp are also 0 so we only need to check for the triplet states with L=Lp
  // double spin_orbit_decomposition(double p, double pp, int S, int L, int Lp, int J,  std::function<double(double,double,double)> potential, int n_z_points)
  // {
  //   double W = 0;
  //   if(S==1)
  //   {
  //     if (L==Lp)
  //     {
  //       if (L==J)
  //       {
  //         return 2*p*pp/(2*J+1)*(A(p,pp,J+1,0,potential,n_z_points)-A(p,pp,J-1,0,potential,n_z_points));
  //       }
  //       else if (L==J-1)
  //       {
  //         return 2*p*pp*(J-1)/(2*J+1)*(A(p,pp,J-2,0,potential,n_z_points)-A(p,pp,J,0,potential,n_z_points));
  //       }
  //       else if (L==J+1)
  //       {
  //         return 2*p*pp*(J+2)/(2*J+1)*(A(p,pp,J+2,0,potential,n_z_points)-A(p,pp,J,0,potential,n_z_points));
  //       }
  //     }
  //   }
  //   return W;
  // }

  // /// Helicity decomposition of a scalar potential with a sigma-L force as 
  // /// described in eq 4.23 of Erkelenz et al.
  // double sigma_L_decomposition(double p, double pp, int S, int L, int Lp, int J,  std::function<double(double,double,double)> potential, int n_z_points)
  // {
  //   double W = 0;
  //   if (S==0 and  L==J and Lp==J) //Singlet state
  //   {
  //     return 2*p*p*pp*pp*(A(p,pp,J,2,potential,n_z_points)-A(p,pp,J,0,potential,n_z_points));
  //   }
  //   else if(S==1) //Triplet state
  //   {
  //     if (L==Lp)
  //     {
  //       if (L==J)
  //       {
  //         W = 2*p*p*pp*pp*(-A(p,pp,J,0,potential,n_z_points)+(J-1)/(2*J+1)*A(p,pp,J+1,1,potential,n_z_points)+(J-1)/(2*J+1)*A(p,pp,J-1,1,potential,n_z_points));
  //       }
  //       else if (L==J-1)
  //       {
  //         W = 2*p*p*pp*pp*((2*J-1)/(2*J+1)*A(p,pp,J-1,0,potential,n_z_points)-2/(2*J+1)*A(p,pp,J,1,potential,n_z_points)-A(p,pp,J-1,1,potential,n_z_points));
  //       }
  //       else if (L==J+1)
  //       {
  //         W = 2*p*p*pp*pp*((2*J+3)/(2*J+1)*A(p,pp,J+1,0,potential,n_z_points)-2/(2*J+1)*A(p,pp,J,1,potential,n_z_points)-A(p,pp,J+1,1,potential,n_z_points));
  //       }
  //     }
  //     else if ((Lp==J-1 and L==J+1) or (Lp==J+1 and L==J-1))
  //     {
  //       W = 4*p*p*pp*pp*sqrt(J*(J+1))/(2*J+1)/(2*J+1)*(A(p,pp,J+1,0,potential,n_z_points)-A(p,pp,J-1,0,potential,n_z_points));      
  //     }
  //   }
  //   return W;
  // }

  // /// Helicity decomposition of a scalar potential with a tensor force as 
  // /// described in eq 4.24 of Erkelenz et al.
  // double tensor_decomposition(double p, double pp, int S, int L, int Lp, int J,  std::function<double(double,double,double)> potential, int n_z_points)
  // {
  //   double W = 0;
  //   if (S==0  and L==J and Lp==J) //Singlet state
  //   {
  //     W = 2*(-(p*p+pp*pp)*A(p,pp,J,0,potential,n_z_points)+2*p*pp*A(p,pp,J,1,potential,n_z_points));
  //   }
  //   else if(S==1) //Triple state
  //   {
  //     if (L==Lp)
  //     {
  //       if (L==J)
  //       {
  //         W = 2*((pp*pp+p*p)*A(p,pp,J,0,potential,n_z_points)-2*p*pp/(2*J+1)*(A(p,pp,J+1,0,potential,n_z_points)+(J+1)*A(p,pp,J-1,0,potential,n_z_points)));
  //       }
  //       else if (L==J-1)
  //       {
  //         W =  2/(2*J+1)*((p*p+pp*pp)*A(p,pp,J-1,0,potential,n_z_points)-2*p*pp*A(p,pp,J,0,potential,n_z_points));
  //       }
  //       else if (L==J+1)
  //       {
  //         W = 2/(2*J+1)*(-(p*p+pp*pp)*A(p,pp,J+1,0,potential,n_z_points)+2*p*pp*A(p,pp,J,0,potential,n_z_points));
  //       }
  //     }
  //     else if (Lp==J-1 and L==J+1) 
  //     {
  //       W = 4*sqrt(J*(J+1))/(2*J+1)*(p*p*A(p,pp,J+1,0,potential,n_z_points)+pp*pp*A(p,pp,J-1,0,potential,n_z_points))-2*p*pp*A(p,pp,J,0,potential,n_z_points);      
  //     }
  //     else if (Lp==J+1 and L==J-1)
  //     {
  //       W = 4*sqrt(J*(J+1))/(2*J+1)*(p*p*A(p,pp,J-1,0,potential,n_z_points)+pp*pp*A(p,pp,J+1,0,potential,n_z_points))-2*p*pp*A(p,pp,J,0,potential,n_z_points);
  //     }
  //   }
  //   return W;
  // }

  // /// Helicity decomposition of a scalar potential with a spin-spin force as 
  // /// described in eq 4.25 of Erkelenz et al.
  // double sigma_k_decomposition(double p, double pp, int S, int L, int Lp, int J,  std::function<double(double,double,double)> potential, int n_z_points)
  // {
  //   double W = 0;
  //   if (S==0 and L==J and Lp==J) //Singlet state
  //   {
  //     W = 1/2*(-(p*p+pp*pp)*A(p,pp,J,0,potential,n_z_points)-2*p*pp*A(p,pp,J,1,potential,n_z_points));
  //   }
  //   else if(S==1) //Triplet state
  //   {
  //     if (L==Lp)
  //     {
  //       if (L==J)
  //       {
  //         W = 1/2*((pp*pp+p*p)*A(p,pp,J,0,potential,n_z_points)+2*p*pp/(2*J+1)*(A(p,pp,J+1,0,potential,n_z_points)+(J+1)*A(p,pp,J-1,0,potential,n_z_points)));
  //       }
  //       else if (L==J-1)
  //       {
  //         W = 1/2/(2*J+1)*((p*p+pp*pp)*A(p,pp,J-1,0,potential,n_z_points)+2*p*pp*A(p,pp,J,0,potential,n_z_points));
  //       }
  //       else if (L==J+1)
  //       {
  //         W = 1/2/(2*J+1)*(-(p*p+pp*pp)*A(p,pp,J+1,0,potential,n_z_points)-2*p*pp*A(p,pp,J,0,potential,n_z_points));
  //       }
  //     }
  //     else if (Lp==J-1 and L==J+1) 
  //     {
  //       W = sqrt(J*(J+1))/(2*J+1)*(p*p*A(p,pp,J+1,0,potential,n_z_points)+pp*pp*A(p,pp,J-1,0,potential,n_z_points))+2*p*pp*A(p,pp,J,0,potential,n_z_points);      
  //     }
  //     else if (Lp==J+1 and L==J-1)
  //     {
  //       W = sqrt(J*(J+1))/(2*J+1)*(p*p*A(p,pp,J-1,0,potential,n_z_points)+pp*pp*A(p,pp,J+1,0,potential,n_z_points))+2*p*pp*A(p,pp,J,0,potential,n_z_points);
  //     }
  //   }
  //   return W;
  // }

  // Functions to allow for precomputation for efficiency
  // overloaded function if we pass the integration table for efficiency when many A's need to be computed
  double A(double p, double pp, int J, int l, std::function<double(double, double, double)> potential, gsl_integration_glfixed_table *t, int n_z_points)
  {
    double A = 0;
    for (int i = 0; i < n_z_points; i++)
    {
      double zi, wi;
      gsl_integration_glfixed_point(-1, 1, i, &zi, &wi, t);
      A += wi * AIntegrand(p, pp, zi, J, l, potential);
    }
    return PI * A;
  }

  uint64_t AHash(int index_p, int index_pp, int J, int l)
  {
    return (((uint64_t)(index_p)) << 30) + (((uint64_t)(index_pp)) << 20) + (((uint64_t)(J)) << 10) + ((uint64_t)(l));
  }

  void AUnHash(uint64_t key, uint64_t &index_p, uint64_t &index_pp, uint64_t &J, uint64_t &l)
  {
    index_p = (key >> 30) & 0x3FFL;
    index_pp = (key >> 20) & 0x3FFL;
    J = (key >> 10) & 0x3FFL;
    l = (key)&0x3FFL;
  }

  std::unordered_map<uint64_t, double> PrecalculateA(int e2max, std::function<double(double, double, double)> potential, int lmax, double max_momentum, gsl_integration_glfixed_table *t_momentum, gsl_integration_glfixed_table *t_z, int n_momentum_points, int n_z_points)
  {
    std::unordered_map<uint64_t, double> AList;
    int Jmax = e2max + 1; // coupling of 2 l=emax with S=1
    std::vector<uint64_t> KEYS;
    for (int J = 0; J <= Jmax; J++)
    {
      for (int index_p = 0; index_p < n_momentum_points; index_p++)
      {
        for (int index_pp = index_p; index_pp < n_momentum_points; index_pp++)
        {
          for (int l = 0; l <= lmax; l++)
          {
            uint64_t key = AHash(index_p, index_pp, J, l);
            KEYS.push_back(key);
            AList[key] = 0.0; // Everything needs to be in this loop to avoid re-hash in parralel
          }
        }
      }
    }
    #pragma omp parallel for schedule(dynamic, 1) // this works as long as the gsl_function handle is within this for-loop
    for (size_t n = 0; n < KEYS.size(); n++)
    {
      uint64_t key = KEYS[n];
      uint64_t index_p, index_pp, J, l;
      AUnHash(key, index_p, index_pp, J, l);
      double p, pp, wi, wj;
      gsl_integration_glfixed_point(0, max_momentum, index_p, &p, &wi, t_momentum);
      gsl_integration_glfixed_point(0, max_momentum, index_pp, &pp, &wj, t_momentum);
      double aval = A(p, pp, J, l, potential, t_z, n_z_points);
      AList[key] =aval;
      // std::stringstream intvalue;
      // intvalue<<p<<" "<<pp<<" "<<J<<" "<<aval<<std::endl;
      // std::cout<<intvalue.str();
    }
    // gsl_integration_glfixed_table_free(t_momentum);
    // gsl_integration_glfixed_table_free(t_z);
    return AList;
  }

  double GetA(int index_p, int index_pp, int J, int l, std::unordered_map<uint64_t, double> &AList)
  {
    double A;
    if (index_p > index_pp)
      std::swap(index_p, index_pp);
    uint64_t key = AHash(index_p, index_pp, J, l);
    auto it = AList.find(key);
    if (it != AList.end()) // return what we've found
    {
      A = it->second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
      std::cout << "DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!" << std::endl;
      std::cout << "   I shouldn't be here in " << __func__ << "(" << index_p << " , " << index_pp << " , " << J << " , " << l << ")   key = " << key << std::endl;
//      printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
//      printf("   I shouldn't be here in GetA(%d, %d, %d, %d):   key =%llx", index_p, index_pp, J, l, key);
      exit(EXIT_FAILURE);
    }
    return A;
  }

  double central_force_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
  {
    double W = 0;
    if (S == 0 and L == J and Lp == J) // Singlet state
    {
      W = 2 * GetA(index_p, index_pp, J, 0, AList);
    }
    else if (S == 1) // Triplet state
    {
      // The fucntions are the same if L=Lp regardless of the value of L and is 0 otherwise
      if ((L == Lp) and (((L == J) or (L == J - 1) or (L == J + 1))))
      {
        if ((L == J) or (L == J - 1) or (L == J + 1))
        {
          W = 2 * GetA(index_p, index_pp, L, 0, AList);
        }
      }
    }
    return W;
  }

  double spin_spin_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
  {
    int Seval = Seval = 2 * S * (S + 1) - 3;
    double W = Seval * central_force_decomposition(p, pp, index_p, index_pp, S, L, Lp, J,AList); 
    return W;
  }

  double spin_orbit_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
  {
    double W = 0;
    if (S == 1)
    {
      if (L == Lp)
      {
        if (L == J)
        {
          return 2 * p * pp / (2 * J + 1) * (GetA(index_p, index_pp, J + 1, 0, AList) - GetA(index_p, index_pp, J - 1, 0, AList));
        }
        else if (L == J - 1)
        {
          return 2 * p * pp * (J - 1) / (2 * J + 1) * (GetA(index_p, index_pp, J - 2, 0, AList) - GetA(index_p, index_pp, J, 0, AList));
        }
        else if (L == J + 1)
        {
          return 2 * p * pp * (J + 2) / (2 * J + 1) * (GetA(index_p, index_pp, J + 2, 0, AList) - GetA(index_p, index_pp, J, 0, AList));
        }
      }
    }
    return W;
  }

  double sigma_L_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
  {
    double W = 0;
    if (S==0 and  L==J and Lp==J) //Singlet state
    {
      return 2*p*p*pp*pp*(GetA(index_p,index_pp,J,2,AList)-GetA(index_p,index_pp,J,0,AList));
    }
    else if(S==1) //Triplet state
    {
      if (L==Lp)
      {
        if (L==J)
        {
          W = 2*p*p*pp*pp*(-GetA(index_p,index_pp,J,0,AList)+(J-1)/(2*J+1)*GetA(index_p,index_pp,J+1,1,AList)+(J-1)/(2*J+1)*GetA(index_p,index_pp,J-1,1,AList));
        }
        else if (L==J-1)
        {
          W = 2*p*p*pp*pp*((2*J-1)/(2*J+1)*GetA(index_p,index_pp,J-1,0,AList)-2/(2*J+1)*GetA(index_p,index_pp,J,1,AList)-GetA(index_p,index_pp,J-1,1,AList));
        }
        else if (L==J+1)
        {
          W = 2*p*p*pp*pp*((2*J+3)/(2*J+1)*GetA(index_p,index_pp,J+1,0,AList)-2/(2*J+1)*GetA(index_p,index_pp,J,1,AList)-GetA(index_p,index_pp,J+1,1,AList));
        }
      }
      else if ((Lp==J-1 and L==J+1) or (Lp==J+1 and L==J-1))
      {
        W = 4*p*p*pp*pp*sqrt(J*(J+1))/(2*J+1)/(2*J+1)*(GetA(index_p,index_pp,J+1,0,AList)-GetA(index_p,index_pp,J-1,0,AList));      
      }
    }
    return W;
  }

  double tensor_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
  {
    double W = 0;
    if (S == 0 and L == J and Lp == J) // Singlet state
    {
      W = 2 * (-(p * p + pp * pp) * GetA(index_p, index_pp, J, 0, AList) + 2 * p * pp * GetA(index_p, index_pp, J, 1, AList));
    }
    else if (S == 1) // Triple state
    {
      if (L == Lp)
      {
        if (L == J)
        {
          W = 2 * ((pp * pp + p * p) * GetA(index_p, index_pp, J, 0, AList) - 2 * p * pp / (2 * J + 1) * (GetA(index_p, index_pp, J + 1, 0, AList) + (J + 1) * GetA(index_p, index_pp, J - 1, 0, AList)));
        }
        else if (L == J - 1)
        {
          W = 2 / (2 * J + 1) * ((p * p + pp * pp) * GetA(index_p, index_pp, J - 1, 0, AList) - 2 * p * pp * GetA(index_p, index_pp, J, 0, AList));
        }
        else if (L == J + 1)
        {
          W = 2 / (2 * J + 1) * (-(p * p + pp * pp) *GetA(index_p, index_pp, J + 1, 0, AList) + 2 * p * pp * GetA(index_p, index_pp, J, 0, AList));
        }
      }
      else if (Lp == J - 1 and L == J + 1)
      {
        W = 4 * sqrt(J * (J + 1)) / (2 * J + 1) * (p * p * GetA(index_p, index_pp, J + 1, 0, AList) + pp * pp * GetA(index_p, index_pp, J - 1, 0, AList)) - 2 * p * pp * GetA(index_p, index_pp, J, 0, AList);
      }
      else if (Lp == J + 1 and L == J - 1)
      {
        W = 4 * sqrt(J * (J + 1)) / (2 * J + 1) * (p * p * GetA(index_p, index_pp, J - 1, 0, AList) + pp * pp * GetA(index_p, index_pp, J + 1, 0, AList)) - 2 * p * pp * GetA(index_p, index_pp, J, 0, AList);
      }
    }
    return W;
  }

  double sigma_k_decomposition(double p, double pp, int index_p, int index_pp, int S, int L, int Lp, int J, std::unordered_map<uint64_t, double>& AList)
  {
    double W = 0;
    if (S == 0 and L == J and Lp == J) // Singlet state
    {
      W = 1 / 2 * (-(p * p + pp * pp) * GetA(index_p, index_pp, J, 0, AList) - 2 * p * pp * GetA(index_p, index_pp, J, 1, AList));
    }
    else if (S == 1) // Triplet state
    {
      if (L == Lp)
      {
        if (L == J)
        {
          W = 1 / 2 * ((pp * pp + p * p) * GetA(index_p, index_pp, J, 0, AList) + 2 * p * pp / (2 * J + 1) * (GetA(index_p, index_pp, J + 1, 0, AList) + (J + 1) * GetA(index_p, index_pp, J - 1, 0, AList)));
        }
        else if (L == J - 1)
        {
          W = 1 / 2 / (2 * J + 1) * ((p * p + pp * pp) * GetA(index_p, index_pp, J - 1, 0, AList) + 2 * p * pp * GetA(index_p, index_pp, J, 0, AList));
        }
        else if (L == J + 1)
        {
          W = 1 / 2 / (2 * J + 1) * (-(p * p + pp * pp) * GetA(index_p, index_pp, J + 1, 0, AList) - 2 * p * pp * GetA(index_p, index_pp, J, 0, AList));
        }
      }
      else if (Lp == J - 1 and L == J + 1)
      {
        W = sqrt(J * (J + 1)) / (2 * J + 1) * (p * p * GetA(index_p, index_pp, J + 1, 0, AList) + pp * pp * GetA(index_p, index_pp, J - 1, 0, AList)) + 2 * p * pp * GetA(index_p, index_pp, J, 0, AList);
      }
      else if (Lp == J + 1 and L == J - 1)
      {
        W = sqrt(J * (J + 1)) / (2 * J + 1) * (p * p * GetA(index_p, index_pp, J - 1, 0, AList) + pp * pp * GetA(index_p, index_pp, J + 1, 0, AList)) + 2 * p * pp * GetA(index_p, index_pp, J, 0, AList);
      }
    }
    return W;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///                                                                                     ///
  ///  Generalization of the above formula to tensor operators.                           ///
  ///  This follows the derivation I did which can be found at                            ///
  ///  https://drive.google.com/file/d/1QQrLwYNSbMaJrA2WRX0wlATcl0mL5MRp/view?usp=sharing ///
  ///                                                                                     ///
  ///////////////////////////////////////////////////////////////////////////////////////////

  double binomial_coefficient(int n, int k)
  {
    return gsl_sf_fact(n)/gsl_sf_fact(k)/gsl_sf_fact(n-k); 
  }

  // Jacobi polynomials as expressed in the Alternate expression for real arguments section
  // on Wikipedia
  double jacobi_polynomials(int n, int a, int b, double z)
  {
    double polynomial = 0;
    for(int s=0;s<=n;s++)
    {
      polynomial += binomial_coefficient(n+a, n-s)*binomial_coefficient(n+b,s)*pow((z-1)/2,s)*pow((z+1)/2, n-s); 
    } 
    return polynomial;
  }

  // Wigener small-d matrix implemented using the formula
  // using Jacobi polynomials on wikipedia
  double wigner_small_d_matrix(int j, int m1, int m2, double z)
  {
    int a, lam;
    int k = std::min(std::min(j+m2,j+m1),std::min(j-m2,j-m1));    
    if (k==j+m2)
    {
      a = m1-m2;
      lam = m1-m2;
    }
    else if (k==j-m2)
    {
      a = m2-m1;
      lam = 0;
    }
    else if (k==j+m1)
    {
      a = m2-m1;
      lam = 0;
    }
    else if (k==j-m1)
    {
      a = m1 - m2;
      lam = m1-m2;
    } 
    double b = 2*j - 2*k - a;
    return phase(lam)*sqrt(binomial_coefficient(2*j-k,k+a)/binomial_coefficient(k+b,b))*pow(sqrt(0.5*(1-z)),a)*pow(sqrt(0.5*(1+z)),b)*jacobi_polynomials(k,a,b,z);

  }

  // Performs the integral over cos(theta) of the wigner small-d matrices and the matrix element of the operator in helicity basis.
  // This corresponds to the integral in eq. 13 in my notes, found at https://drive.google.com/file/d/1QQrLwYNSbMaJrA2WRX0wlATcl0mL5MRp/view?usp=sharing
  double integrate_op_rotation (double p, double pp, std::function<double(double, double, double,int,int,double,double,double,double)> op, int Jp, int operator_index, int operator_rank, double lam1, double lam2, double lamp1, double lamp2, double Lam, double Lamp)
  {
    double integral = 0;
    int size = 100;
    gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(size);
    for (int i=0;i<size;i++)
    {
      double zi, wi;
      gsl_integration_glfixed_point(-1,1,i,&zi,&wi,t);
      integral += wi*op(p,pp, zi, operator_index, operator_rank, lam1, lam2, lamp1, lamp2)*wigner_small_d_matrix(Jp, operator_index+Lam, Lamp, zi);
    }
    gsl_integration_glfixed_table_free(t);
    return integral;
  }

  // Partial Wave Decomposition of a general tensor operator in the J basis using the helicity representation of the operator. 
  // Corresponds to equation 13 of my notes.
  double projection_to_j_basis(double p, double pp, int J, int Jp, double lam1, double lam2, double lamp1, double lamp2, double Lam, double Lamp,  std::function<double(double, double, double, int,int,double,double,double,double)> op, int operator_rank)
  {
    double tensor_representation = 0;
    if (abs(Lamp)>Jp) return 0;
    for (int u = -operator_rank; u <= operator_rank; u++)
    {
      if (abs(u+Lam)>Jp) continue;
      tensor_representation += AngMom::CG(J,Lam,operator_rank,u,Jp,u+Lam)*integrate_op_rotation(p, pp, op, Jp, u, operator_rank, lam1, lam2, lamp1, lamp2, Lam, Lamp);
    }
    tensor_representation *= 2.0*PI*sqrt(2.0*J+1.0);
    return tensor_representation;
  }

  // Compute the overlap of the states in heliciy basis and in lsj basis.
  // It is used at to convert the result in helicity basis to lsj basis in the 
  // full partial wave_decomposition. Corresponds to equation 14 of my notes.
  double overlap_lsj(int S, int L, int J, double lam1, double lam2, double Lam)
  {
    if (S < abs(Lam)) return 0;
    double overlap = sqrt((2.0*L+1.0)/(2.0*J+1.0))*AngMom::CG(L,0,S,Lam,J,Lam)*AngMom::CG(0.5,lam1,0.5,-lam2,S,Lam);
    return overlap;
  }

  // Full partial wave decomposition of a general tensor operator using helicity formalism.
  // Gives the matrix element of the opertor in lsj basis. Corresponds to equation 15 of my notes.
  double tensor_partial_wave_decomposition(double p, double pp, int L, int Lp, int S, int Sp, int J, int Jp,  std::function<double(double, double, double, int,int,double,double,double,double)> op, int operator_rank)
  {
    // There are 4 possible heliciies for each state, i.e. ++, +-. -+, --
    // The following arrrays allow to iterate over those helicites states.
    int helicity1[4] = {1,1,-1,-1};
    int helicity2[4] = {1,-1,1,-1};
    double result = 0;
    double lam1, lam2, Lam, lamp1, lamp2, Lamp;
    double overlap_bra, overlap_ket;
    for (int i=0; i<4; i++)
    {
      lam1 = 0.5*helicity1[i];
      lam2 = 0.5*helicity2[i];
      Lam = (lam1-lam2);
      if (S < abs(Lam)) continue;
      overlap_ket = overlap_lsj(S,L,J,lam1,lam2,Lam);
      if (abs(overlap_ket)<1e-16) continue;
      for (int j=0; j<4; j++)
      {
        lamp1 = 0.5*helicity1[j];
        lamp2 = 0.5*helicity2[j];
        Lamp = (lamp1-lamp2); 
        if (Sp < abs(Lamp)) continue;
        overlap_bra = overlap_lsj(Sp,Lp,Jp,lamp1,lamp2,Lamp);
        if (abs(overlap_bra)<1e-16) continue;
        double term = overlap_bra*projection_to_j_basis(p,pp,J,Jp,lam1,lam2,lamp1,lamp2,Lam,Lamp,op,operator_rank)*overlap_ket;
        result += term;
      }
    }
    if ((Lp-L)%2 == 0) result *= phase((Lp-L)/2); //comes from Fourier trans
    if ((Lp-L)%2 == 1) result *= phase((Lp-L-3)/2); //comes from Fourier trans, (-i) is factored out
    return result;
  }




  //////////////////////////////////////////////////////////////////////////////////////
  ///                                                                                ///
  ///  Helicity functions, useful to write the tensor operators in helicity basis.   ///
  ///  They are implemented using equation 4.11 and following in Erkelenz et al.     ///
  ///  All the following are written in terms of z=cos(theta).                       ///
  ///                                                                                ///
  //////////////////////////////////////////////////////////////////////////////////////

  double helicity_expectation_value_identity(double lam, double lamp,  double z)
  {
    return abs(lamp+lam) * sqrt(0.5*(1+z)) + (lamp-lam) * sqrt(0.5*(1-z));
  }

  double helicity_expectation_value_x(double lam, double lamp, double z)
  {
    return abs(lamp-lam) * sqrt(0.5*(1+z)) + (lamp+lam) * sqrt(0.5*(1-z));
  }

  //Factor of i is not taken into account
  double helicity_expectation_value_y(double lam, double lamp, double z)
  {
    return (lam-lamp) * sqrt(0.5*(1+z)) + abs(lamp+lam) * sqrt(0.5*(1-z));
  }

  double helicity_expectation_value_z(double lam, double lamp, double z)
  {
    return (lamp+lam) * sqrt(0.5*(1+z)) - abs(lamp-lam) * sqrt(0.5*(1-z));
  }

  double helicity_expectation_value_sigma(double lam, double lamp, double z, int m)
  {
    double result;
    if (m==1)
    {
      result = -(helicity_expectation_value_x(lam, lamp, z) - helicity_expectation_value_y(lam, lamp, z))/sqrt(2);
    }
    else if (m==0)
    {
      result = helicity_expectation_value_z(lam, lamp, z);
    }
    else if (m==-1)
    {
      result = (helicity_expectation_value_x(lam, lamp, z) + helicity_expectation_value_y(lam, lamp, z))/sqrt(2);
    }
    else
    {
      result = 0;
    }
    return result;
  }

  double identity_helicity_representaition(double lam1, double lam2, double lamp1, double lamp2, double z)
  {
    return helicity_expectation_value_identity(lam1, lamp1, z) * helicity_expectation_value_identity(-lam2, -lamp2, z);  
  }

  //s1 . s2
  double sigma_dot_sigma(double lam1, double lam2, double lamp1, double lamp2, double z)
  {
    return    helicity_expectation_value_x(lam1, lamp1, z)*helicity_expectation_value_x(-lam2, -lamp2, z)
            - helicity_expectation_value_y(lam1, lamp1, z)*helicity_expectation_value_y(-lam2, -lamp2, z)
            + helicity_expectation_value_z(lam1, lamp1, z)*helicity_expectation_value_z(-lam2, -lamp2, z);
  }

  // (s1 + s2) . q 
  double sigma_plus_sigma_dot_q(double lam1, double lam2, double lamp1, double lamp2, double p_out, double p_in, double z)
  {
    return   (p_out * (lamp1-lamp2) - p_in * (lam1-lam2)) 
           * helicity_expectation_value_identity(lam1, lamp1, z)
           * helicity_expectation_value_identity(-lam2, -lamp2, z);
  }

  // (s1 + s2) . (p_out + p_in)
  double sigma_plus_sigma_dot_pout_plus_pin(double lam1, double lam2, double lamp1, double lamp2, double p_out, double p_in, double z)
  {
    return   (p_out * (lamp1-lamp2) + p_in * (lam1-lam2)) 
           * helicity_expectation_value_identity(lam1, lamp1, z)
           * helicity_expectation_value_identity(-lam2, -lamp2, z);
  }

  // (s1 - s2) . q
  double sigma_minus_sigma_dot_q(double lam1, double lam2, double lamp1, double lamp2, double p_out, double p_in, double z)
  {
    return   (p_out * (lamp1+lamp2) - p_in * (lam1+lam2)) 
           * helicity_expectation_value_identity(lam1, lamp1, z)
           * helicity_expectation_value_identity(-lam2, -lamp2, z);
  }

  // (s1 - s2) . (p_out+p_in)
  double sigma_minus_sigma_dot_pout_plus_pin(double lam1, double lam2, double lamp1, double lamp2, double p_out, double p_in, double z)
  {
    return   (p_out * (lamp1+lamp2) + p_in * (lam1+lam2)) 
           * helicity_expectation_value_identity(lam1, lamp1, z)
           * helicity_expectation_value_identity(-lam2, -lamp2, z);
  }

  // (s1 x s2) . q
  double sigma_cross_sigma_dot_q(double lam1, double lam2, double lamp1, double lamp2, double p_out, double p_in, double z)
  {
    return    (p_out*z-p_in)
            * (
                 helicity_expectation_value_x(lam1, lamp1, z) * helicity_expectation_value_y(-lam2, -lamp2, z)
               - helicity_expectation_value_y(lam1, lamp1, z) * helicity_expectation_value_x(-lam2, -lamp2, z)
              )
            * p_out * sqrt(1-z*z)
            * (
                 helicity_expectation_value_y(lam1, lamp1, z) * helicity_expectation_value_z(-lam2, -lamp2, z)
               - helicity_expectation_value_z(lam1, lamp1, z) * helicity_expectation_value_y(-lam2, -lamp2, z)
              );
  }

  // (s1 x s2) . (pout+pin)
  double sigma_cross_sigma_dot_pout_plus_pin(double lam1, double lam2, double lamp1, double lamp2, double p_out, double p_in, double z)
  {
    return    (p_out*z+p_in)
            * (
                 helicity_expectation_value_x(lam1, lamp1, z) * helicity_expectation_value_y(-lam2, -lamp2, z)
               - helicity_expectation_value_y(lam1, lamp1, z) * helicity_expectation_value_x(-lam2, -lamp2, z)
              )
            * p_out * sqrt(1-z*z)
            * (
                 helicity_expectation_value_y(lam1, lamp1, z) * helicity_expectation_value_z(-lam2, -lamp2, z)
               - helicity_expectation_value_z(lam1, lamp1, z) * helicity_expectation_value_y(-lam2, -lamp2, z)
              );
  }

  // (s1 . q) * (s2 . q )
  double sigma_dot_q_sigma_dot_q(double lam1, double lam2, double lamp1, double lamp2, double p_out, double p_in, double z)
  {
    return - (p_out*lamp1 - p_in*lam1) * (p_out*lamp2 - p_in*lam2)
           * helicity_expectation_value_identity(lam1, lamp1, z)
           * helicity_expectation_value_identity(-lam2, -lamp2, z);
  }

 


  // Operator to benchmark the method [[s1*q]^2[s2*q]^2]^K_M.
  double test_operator_helicity_basis(double p, double pp, double z, int operator_index, int operator_rank, double lam1, double lam2, double lamp1, double lamp2)
  { 
    double orbitpart;
    double qsquared =  pp*pp+p*p-2*p*pp*z;
    double result = 0;
    for (int s1 = -1;s1 <= 1; s1++)
    {
      for (int m1 = -1; m1 <= 1; m1++)
      {
        int mk1 = m1+s1;
        int mk2 = operator_index-mk1;
        if (abs(mk2) > 2 ) continue;
        double coupling_ms1 = AngMom::CG(1,m1,1,s1,2,mk1);
        double coupling_mk = AngMom::CG(2,mk1,2,mk2,operator_rank,operator_index);
        for (int s2 = -1; s2 <= 1; s2++)
        {
          int m2 = operator_index - mk1 - s2;
          if (abs(m2) > 1) continue;
          double CGs = coupling_ms1 * AngMom::CG(1,m2,1,s2,2,mk2) * coupling_mk;
          if (abs(CGs) < 1e-16) continue;
          if ((abs(m1)==1) && (abs(m2)==1))
          {
            orbitpart = 0.5 * pp * pp * (-m1) * (-m2) * (1-z*z);
          }
          else if ((abs(m1)==1) && (m2==0))
          {
            orbitpart = (-m1)* pp * sqrt(1-z*z)*(pp * z - p )/sqrt(2);
          }
          else if ((m1==0) && (abs(m2)==1))
          {
             orbitpart = (-m2)* pp * sqrt(1-z*z)*(pp * z - p )/sqrt(2);
          }
          else
          {
             orbitpart = ( pp * pp * z*z - 2*pp * p * z + p * p );
          }
          result += CGs * orbitpart * helicity_expectation_value_sigma(lam1, lamp1, z, s1) * helicity_expectation_value_sigma(-lam2, -lamp2, z, s2)/(M_PION_NEUTRAL*M_PION_CHARGED+qsquared);
        } 
      }
    } 
    return result;
  } 
    

}
