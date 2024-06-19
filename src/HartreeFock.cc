
#include "HartreeFock.hh"
#include "ModelSpace.hh"
#include "PhysicalConstants.hh"
#include "AngMom.hh"
#include <iomanip>
#include <vector>
#include <array>
#include <map>
#include <utility> // for make_pair
//#include "gsl/gsl_sf_gamma.h" // for radial wave function
#include "gsl/gsl_sf_laguerre.h" // for radial wave function
//#include <gsl/gsl_math.h> // for M_SQRTPI
#include <omp.h>


//using namespace std;

HartreeFock::HartreeFock(Operator& hbare)
  : Hbare(hbare), modelspace(hbare.GetModelSpace()), ms_for_output_3N(hbare.GetModelSpace()), 
    KE(Hbare.OneBody), energies(Hbare.OneBody.diag()),
    tolerance(1e-8), convergence_ediff(7,0), convergence_EHF(7,0), freeze_occupations(true),discard_NO2B_from_3N(false)
{
   int norbits = modelspace->GetNumberOrbits();

   C             = arma::mat(norbits,norbits,arma::fill::eye);
   Vij           = arma::mat(norbits,norbits,arma::fill::zeros);
   V3ij          = arma::mat(norbits,norbits,arma::fill::zeros);
   F             = arma::mat(norbits,norbits);
   for (int Tz=-1;Tz<=1;++Tz)
   {
     for (int parity=0; parity<=1; ++parity)
     {
       int nKetsMon = modelspace->MonopoleKets[Tz+1][parity].size();
       Vmon[Tz+1][parity] = arma::mat(nKetsMon,nKetsMon,arma::fill::zeros);
       Vmon_exch[Tz+1][parity] = arma::mat(nKetsMon,nKetsMon,arma::fill::zeros);
     }
   }
   prev_energies = arma::vec(norbits,arma::fill::zeros);
   std::vector<double> occvec;
   for (auto& h : modelspace->holes) occvec.push_back(modelspace->GetOrbit(h).occ);
   holeorbs = arma::uvec( std::vector<index_t>(modelspace->holes.begin(),modelspace->holes.end()));
   hole_occ = arma::rowvec(occvec);


}


//*********************************************************************
/// Diagonalize and update the Fock matrix until convergence.
/// Then, call ReorderCoefficients() to make sure the index
/// ordering and phases are preserved in the transformation
/// from the original basis to the Hatree-Fock basis.
//*********************************************************************
void HartreeFock::Solve()
{
   double t_start = omp_get_wtime();
   iterations = 0; // counter so we don't go on forever
   int maxiter = 1000;
   double density_mixing_factor = 0.2;
   double field_mixing_factor = 0.0;
   int fill_modulus = 3;

   BuildMonopoleV();
   if (Hbare.GetParticleRank()>2)
   {
      BuildMonopoleV3();
   }
   UpdateDensityMatrix();
   UpdateF();


   for (iterations=0; iterations<maxiter; ++iterations)
   {
      Diagonalize();          // Diagonalize the Fock matrix
      ReorderCoefficients();  // Reorder columns of C so we can properly identify the hole orbits.

      // This bit is occationally the source of non-convergence of the iterations
      if (iterations > 10  and (iterations%fill_modulus)==0 and not freeze_occupations)
      {
         FillLowestOrbits();// if we don't freeze the occupations, then calculate the new ones.
         fill_modulus = fill_modulus %5;
         fill_modulus+=1; // to avoid cycles, we keep increasing the interval between checking if we should change the occupations.
      }

      if (iterations == 500)
      {
        density_mixing_factor = 0.7;
        field_mixing_factor = 0.5;
        std::cout << "Still not converged after 500 iterations. Setting density_mixing_factor => " << density_mixing_factor
                  << " field_mixing_factor => " << field_mixing_factor << std::endl;

      }
      if (iterations>600 and iterations%10 == 3 ) // if we've made it to 600, really put on the brakes with a big mixing factor, with random noise
      {
        field_mixing_factor = 1 - 0.005 * (std::rand() % 100);
        density_mixing_factor = 1- 0.005 * ( std::rand() % 100);
      }
     if (iterations==800) // Things are really not healthy at this point. Most likely we're stuck in a occupation changing circle.
     {
        std::cout << "Not converged after 800 iterations. That's it, I'm freezing the occupations." << std::endl;
        freeze_occupations = true;
     }

      if (iterations >100 and (DIIS_error_mats.size()<1 or arma::norm( DIIS_error_mats.back(),"fro")>0.01)  )
      {
         if (iterations==100)
         {
           std::cout << "Still not converged after 100 iterations. Switching to DIIS." << std::endl;
           C.eye(); // reset the density matrix to the starting point
           UpdateDensityMatrix();
           UpdateF();
         }
         UpdateDensityMatrix_DIIS();
         if ( arma::norm( DIIS_error_mats.back(),"fro")<0.01)
         {
           std::cout << "DIIS error matrix below 0.01, switching back to simpler SCF algorithm." << std::endl;
         }
      }
      else
      {
        arma::mat last_rho = rho;
        UpdateDensityMatrix();  // Update the 1 body density matrix, used in UpdateF()
        rho = (1.0 - density_mixing_factor)*rho + density_mixing_factor*last_rho;
      }
      arma::mat last_F = F;
      UpdateF();              // Update the Fock matrix
      F = (1.0 - field_mixing_factor)*F + field_mixing_factor * last_F;

      if ( CheckConvergence() ) break;
   }
   CalcEHF();

   std::cout << std::setw(15) << std::setprecision(10);
   if (iterations < maxiter)
   {
      std::cout << "HF converged after " << iterations << " iterations. " << std::endl;
   }
   else
   {
      std::cout << "!!!! Warning: Hartree-Fock calculation didn't converge after " << iterations << " iterations." << std::endl;
      std::cout << "!!!! Last " << convergence_ediff.size() << " points in convergence check:";
      for (auto& x : convergence_ediff ) std::cout << x << " ";
      std::cout << "  (tolerance = " << tolerance << ")" << std::endl;
      std::cout << "!!!! Last " << convergence_EHF.size() << "  EHF values: ";
      for (auto& x : convergence_EHF ) std::cout << x << " ";
      std::cout << std::endl;

   }
   UpdateReference();
   PrintEHF();
   profiler.timer["HF_Solve"] += omp_get_wtime() - t_start;
}


//*********************************************************************
/// Calculate the HF energy.
/// \f{eqnarray*} E_{HF} &=& \sum_{\alpha} t_{\alpha\alpha}
///                    + \frac{1}{2}\sum_{\alpha\beta} V_{\alpha\beta\alpha\beta}
///                    + \frac{1}{6}\sum_{\alpha\beta\gamma} V_{\alpha\beta\gamma\alpha\beta\gamma} \\
///    &=& \sum_{ij} (2j_i+1) \rho_{ij} ( t_{ij} +\frac{1}{2}\tilde{V}^{(2)}_{ij} + \frac{1}{6}\tilde{V}^{(3)}_{ij} )
/// \f}
/// Where the matrices \f{eqnarray*}
///  \tilde{V}^{(2)}_{ij} &=& \sum_{ab} \rho_{ab}\bar{V}^{(2)}_{iajb} \\
///  \tilde{V}^{(3)}_{ij} &=& \sum_{abcd} \rho_{ab}\rho_{cd} \bar{V}^{(3)}_{iacjbd} \\
/// \f}
/// have already been calculated by UpdateF().
//*********************************************************************
void HartreeFock::CalcEHF()
{
   EHF = 0;
   e1hf = 0;
   e2hf = 0;
   e3hf = 0;
   for (auto i : modelspace->all_orbits )
   {
      Orbit& oi = modelspace->GetOrbit(i);
      int jfactor = oi.j2 +1;
      for (int j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}))
      {
         e1hf += rho(i,j) * jfactor * KE(i,j);
         e2hf += rho(i,j) * jfactor * 0.5 * Vij(i,j);
         e3hf += rho(i,j) * jfactor * (1./6)*V3ij(i,j);
      }
   }
   EHF = e1hf + e2hf + e3hf;
}

//**************************************************************************************
/// Print out the Hartree Fock energy, and the 1-, 2-, and 3-body contributions to it.
//**************************************************************************************
void HartreeFock::PrintEHF()
{
   std::cout << std::fixed <<  std::setprecision(7);
   std::cout << "e1hf = " << e1hf << std::endl;
   std::cout << "e2hf = " << e2hf << std::endl;
   std::cout << "e3hf = " << e3hf << std::endl;
   std::cout << "EHF = "  << EHF  << std::endl;
}

//*********************************************************************
/// [See Suhonen eq. 4.85]
/// Diagonalize the fock matrix \f$ <a|F|b> \f$ and put the
/// eigenvectors in \f$C(i,\alpha) = <i|\alpha> \f$
/// and eigenvalues in the vector energies.
/// Save the last vector of energies to check
/// for convergence.
/// Submatrices corresponding to different channels are diagonalized independently.
/// This guarantees that J,Tz, and \f$ \pi \f$ remain good.
//*********************************************************************
void HartreeFock::Diagonalize()
{
//  std::cout << __func__ << " begin " << std::endl;
   prev_energies = energies;
   for ( auto& it : Hbare.OneBodyChannels )
   {
      arma::uvec orbvec(std::vector<index_t>(it.second.begin(),it.second.end()) ); // convert std::set to std::vector, and convert that to arma::uvec...
      arma::mat F_ch = F.submat(orbvec,orbvec);
      arma::mat C_ch;
      arma::vec E_ch;
      bool success = false;
      int diag_tries = 0;
      while ( not success)
      {
         success = arma::eig_sym(E_ch, C_ch, F_ch);
         ++diag_tries;
         if (diag_tries > 5)
         {
           std::cout << "Hartree-Fock: Failed to diagonalize the submatrix "
                << " on iteration # " << iterations << ". The submatrix looks like:" << std::endl;
           F_ch.print();
           break;
         }
      }
      // Update the full overlap matrix C and energy vector
      energies(orbvec) = E_ch;
      C.submat(orbvec,orbvec) = C_ch;

   }
}


//*********************************************************************
/// Construct an unnormalized two-body monopole interaction
/// \f[ \langle ab | \bar{V}^{(2)} | cd \rangle
///   = \sqrt{(1+\delta_{ab})(1+\delta_{cd})} \sum_{J} (2J+1) \langle ab | V^{(2)} | cd \rangle_{J} \f]
/// This method utilizes the operator method  TwoBodyME::GetTBMEmonopole()
///
//*********************************************************************
void HartreeFock::BuildMonopoleV()
{
   for (int Tz=-1; Tz<=1; ++Tz)
   {
     for (int parity=0; parity<=1; ++parity)
     {
        Vmon[Tz+1][parity].zeros();
        Vmon_exch[Tz+1][parity].zeros();
        for ( auto& itbra : modelspace->MonopoleKets[Tz+1][parity] )
        {
           Ket & bra = modelspace->GetKet(itbra.first);
           int a = bra.p;
           int b = bra.q;
           Orbit & oa = modelspace->GetOrbit(a);
           Orbit & ob = modelspace->GetOrbit(b);
           double norm = (oa.j2+1)*(ob.j2+1);
           for ( auto& itket : modelspace->MonopoleKets[Tz+1][parity] )
           {
              if (itket.second < itbra.second) continue;
              Ket & ket = modelspace->GetKet(itket.first);
              int c = ket.p;
              int d = ket.q;
              Vmon[Tz+1][parity](itbra.second,itket.second)       = Hbare.TwoBody.GetTBMEmonopole(a,b,c,d)*norm;
              Vmon_exch[Tz+1][parity](itbra.second,itket.second)  = Hbare.TwoBody.GetTBMEmonopole(a,b,d,c)*norm;
           }
        }
        Vmon[Tz+1][parity] = arma::symmatu(Vmon[Tz+1][parity]);
        Vmon_exch[Tz+1][parity] = arma::symmatu(Vmon_exch[Tz+1][parity]);
    }
  }
}



/// Set up the keys for the 3-body monopole terms
/// We do this separately because it's fast (so it doesn't need to be parallel)
/// and because it's the same whether we use lab frame 3N or jacobi 3N.
void HartreeFock::SetUpMonopoleV3Keys()
{
   double start_time = omp_get_wtime();
  // First, allocate. This is fast so don't parallelize.
  for ( auto i : modelspace->all_orbits)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    int ei = 2*oi.n + oi.l;
    for (uint64_t j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      if (j<i) continue;
      Orbit& oj = modelspace->GetOrbit(j);
      int ej = 2*oj.n + oj.l;

      for (auto a : modelspace->all_orbits)
      {
        Orbit& oa = modelspace->GetOrbit(a);
        int ea = 2*oa.n + oa.l;
        for (uint64_t b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
        {
          Orbit& ob = modelspace->GetOrbit(b);
          int eb = 2*ob.n + ob.l;

            for ( auto c : modelspace->all_orbits)
            {
              Orbit& oc = modelspace->GetOrbit(c);
              int ec = 2*oc.n + oc.l;
              if ( ea+ec+ei > Hbare.E3max ) continue;
              for (uint64_t d : Hbare.OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) )
              {
                Orbit& od = modelspace->GetOrbit(d);
                int ed = 2*od.n + od.l;

                if ( eb+ed+ej > Hbare.E3max ) continue;
                if ( (oi.l+oa.l+ob.l+oj.l+oc.l+od.l)%2 >0) continue;
                uint64_t key = Vmon3Hash(a,c,i,b,d,j);
                Vmon3_keys.push_back( key );
              }
            }
          }
        }
      }
    }

   profiler.timer[std::string(__func__)] += omp_get_wtime() - start_time;

}


//*********************************************************************
/// Construct an unnormalized three-body monopole interaction
/// \f[ \langle iab | \bar{V}^{(3)} | jcd \rangle =
///     \sum\limits_{J,J_{12}}\sum_{Tt_{12}}(2J+1)(2T+1)
///       \langle (ia)J_{12}t_{12};b JT| V^{(3)} | (jc)J_{12}t_{12}; d JT\rangle
/// \f]
//*********************************************************************
void HartreeFock::BuildMonopoleV3()
{

   SetUpMonopoleV3Keys();
   double start_time = omp_get_wtime();

   Vmon3.resize( Vmon3_keys.size(), 0. );
   profiler.timer["HF_BuildMonopoleV3_allocate"] += omp_get_wtime() - start_time;

   #pragma omp parallel for schedule(dynamic,1)
     for (size_t ind=0; ind<Vmon3.size(); ++ind)
     {
       double v=0;
       int a,b,c,d,i,j;
       Vmon3UnHash( Vmon3_keys[ind], a,c,i,b,d,j);

       int j2i = modelspace->GetOrbit(i).j2;

       v = Hbare.ThreeBody.GetME_pn_mono(a,c,i,b,d,j);

       v /= j2i+1.0;
       #pragma omp atomic write
       Vmon3[ind] = v ;
     }
   std::cout << "HartreeFock::BuildMonopoleV3  storing " << Vmon3.size() << " doubles for Vmon3 and "
             << Vmon3_keys.size() << " uint64's for Vmon3_keys." << std::endl;

   profiler.timer[std::string(__func__)] += omp_get_wtime() - start_time;
}



//*********************************************************************
/// Construct the V3 monopole matrix elements directly from the jacobi matrix elements.
/// So far this is completely untested and certainly wrong.
//*********************************************************************
void HartreeFock::BuildMonopoleV3Jacobi()
{
   SetUpMonopoleV3Keys();
   double start_time = omp_get_wtime();
   Vmon3.resize( Vmon3_keys.size(), 0. );

   for ( auto& obc_a : Hbare.OneBodyChannels )
   {
     int j2a = obc_a.first[1];
     for ( auto& obc_b : Hbare.OneBodyChannels )
     {
       int j2b = obc_b.first[1];
       int Jab_min = std::abs( j2a-j2b )/2;
       int Jab_max = ( j2a+j2b )/2;
       for ( auto& obc_c : Hbare.OneBodyChannels )
       {
         int j2c = obc_c.first[1];

         // Now we compute the relevalt T coefficients
         // These are only a small subset, and we discard them after use
         // so that the storage doesn't kill us.
         // The T coefficents in this sub-block can be labeled by: na,nb,nc, Jab, twoJ, E12, alpha, Ncm, Lcm
         std::unordered_map<size_t,double> Ttable;
         // We'll need one loop to allocate (not thread-safe), one loop to compute the Tcoefficients (can be parallel), and one to compute V3mon (also parallel).
         // first, allocate
         for ( int Jab=Jab_min; Jab<=Jab_max; Jab++)
         {
           int twoJ_min = std::abs(2*Jab-j2c);
           int twoJ_max = (2*Jab-j2c);
           for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
           {
            for (auto a : obc_a.second )
            {
             Orbit& oa = modelspace->GetOrbit(a);
             for (auto b : obc_b.second )
             {
              Orbit& ob = modelspace->GetOrbit(b);
              for (auto c : obc_c.second )
              {
               Orbit& oc = modelspace->GetOrbit(c);
               int Eabc = 2*(oa.n+ob.n+oc.n)+(oa.l+ob.l+oc.l);
               int Lcm_max = std::min(oa.l+ob.l+oc.l, Eabc);
               for (int Lcm=0; Lcm<=Lcm_max; Lcm++)
               {
                for ( int Ncm=0; 2*Ncm<=Eabc-Lcm; Ncm++)
                {
                  int E12 = Eabc - 2*Ncm - Lcm;
                  // TODO: Here's the tricky bit because of the isospin
                  int twoT = 1;
                  int alpha_dim = jacobi3bme->GetDimensionNAS( twoT, twoJ, E12%2, E12 );
                  for ( int alpha=0; alpha<alpha_dim; alpha++)
                  {
                    auto key = Vmon3JacobiHash(oa.n,ob.n,oc.n,Jab,twoJ,E12,alpha,Ncm,Lcm);
                    Ttable[ key ] = 0; // make sure there's somewhere to put the value when we calculate it.
                  } // for alpha
                } // for Ncm
               } // for Lcm
              } // for c
             } // for b
            } // for a
           } // for twoJ
         } // for Jab
         // Now the heavy lifting of computing the T coefficients
         // make a vector of keys so we can loop it in parallel
         std::vector<size_t> keyvector;
         keyvector.reserve( Ttable.size() );
         for ( auto& iter : Ttable ) keyvector.push_back(iter.first);
// parallelize here...
//         #pragma omp parallel for schedule(dynamic, 1)
         for (size_t ikey=0; ikey<keyvector.size(); ikey++)
         {
           auto key = keyvector[ikey];
           uint64_t na,nb,nc,Jab,twoJ,E12,alpha,Ncm,Lcm;
           Vmon3JacobiUnHash(key, na,nb,nc,Jab,twoJ,E12,alpha,Ncm,Lcm);
           jacobi1_state jac1;
           jacobi2_state jac2;
           int twoT = 1; // TODO: This is nonsense. Fix it.
           jacobi3bme->GetJacobiStates( twoT, twoJ, E12%2, E12, alpha, jac1, jac2);
           Orbit& oa = modelspace->GetOrbit( modelspace->GetOrbitIndex(na,obc_a.first[0],obc_a.first[1],obc_a.first[2]) );
           Orbit& ob = modelspace->GetOrbit( modelspace->GetOrbitIndex(nb,obc_b.first[0],obc_b.first[1],obc_b.first[2]) );
           Orbit& oc = modelspace->GetOrbit( modelspace->GetOrbitIndex(nc,obc_c.first[0],obc_c.first[1],obc_c.first[2]) );
           Ket3 ket3(oa,ob,oc);
           double Tcoef = jacobi3bme->Tcoeff_wrapper( ket3, Jab, twoJ, jac1, jac2, twoJ, Ncm, Lcm);
           Ttable[key] = Tcoef;
         }// for ikey
         // and now that we have all the relevant T coefficients computed, let's get all the 3-body monopole terms
         // we can again parallelize here
//         size_t asize = obc_a.second.size();
//         for ( auto a : obc_a.second )
//         #pragma omp parallel for schedule(dynamic,1)
//         for ( size_t ia=0; ia<asize; ia++)  // do it this way because omp doesn't play nice with range-based for loops
         for ( auto a : obc_a.second )
         {
//           auto a = obc_a.second[ia];
           Orbit& oa = modelspace->GetOrbit(a);
           for (auto b : obc_b.second )
           {
             Orbit& ob = modelspace->GetOrbit(b);
             for (auto c : obc_c.second )
             {
               Orbit& oc = modelspace->GetOrbit(c);
               int Eabc = 2*(oa.n+ob.n+oc.n) + (oa.l+ob.l+oc.l);
               for (auto d : obc_a.second )
               {
                 Orbit& od = modelspace->GetOrbit(d);
                 for (auto e : obc_b.second )
                 {
                   Orbit& oe = modelspace->GetOrbit(e);
                   for (auto f : obc_c.second )
                   {
                     if ( f < c ) continue;
                     Orbit& of = modelspace->GetOrbit(f);
                     int Edef = 2*(od.n+oe.n+of.n) + (od.l+oe.l+of.l);
                     double vmon3 = 0;
                     for ( int Jab=Jab_min; Jab<=Jab_max; Jab++ )
                     {
                       int twoJ_min = std::abs(2*Jab-j2c);
                       int twoJ_max = (2*Jab-j2c);
                       for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
                       {
                         for ( int Ncm=0; 2*Ncm<=std::min(Eabc,Edef); Ncm++)
                         {
                           for (int Lcm=0; Lcm<=std::min(Eabc,Edef)-2*Ncm; Lcm+=2)
                           {
                             int E12_abc = Eabc - 2*Ncm - Lcm;
                             int E12_def = Edef - 2*Ncm - Lcm;
                             if ( (E12_abc + Eabc +Lcm)%2 >0) continue; // check parity
                             int twoT = 1;
                             size_t alpha_dim_abc = jacobi3bme->GetDimensionNAS( twoT, twoJ, E12_abc%2, E12_abc );
                             size_t alpha_dim_def = jacobi3bme->GetDimensionNAS( twoT, twoJ, E12_def%2, E12_def );
                             for ( size_t alpha_abc=0; alpha_abc<alpha_dim_abc; alpha_abc++)
                             {
                               auto key_abc = Vmon3JacobiHash(oa.n, ob.n, oc.n, Jab, twoJ,  E12_abc, alpha_abc, Ncm, Lcm);
                               double T_abc = Ttable[key_abc];
                               for ( size_t alpha_def=0; alpha_def<alpha_dim_def; alpha_def++)
                               {
                                 auto key_def = Vmon3JacobiHash(od.n, oe.n, of.n, Jab, twoJ,  E12_def, alpha_def, Ncm, Lcm);
                                 double T_def = Ttable[key_def];
                                 double me_jacobi_NAS = jacobi3bme->GetMatElNAS( alpha_abc, alpha_def, E12_abc, E12_def, twoJ, twoT, E12_abc%2);
                                 vmon3 += (twoJ+1) * T_abc * me_jacobi_NAS * T_def;
                               }
                             }
                           }
                         }
                       }
                     }
                     // Probably should think of a more clever way to do the lookup here
                     auto key = Vmon3Hash(a,b,c,d,e,f);
                     auto ptr = std::find(Vmon3_keys.begin(),Vmon3_keys.end(),key) ;
                     Vmon3[ ptr- Vmon3_keys.begin() ] = vmon3;
                   }
                 }
               }
             }
           }
         }
       }// for obc_c
     }// for obc_b
   }// for obc_a

   profiler.timer[std::string(__func__)] += omp_get_wtime() - start_time;
}









//*********************************************************************
/// Hashing function for rolling six orbit indices into a single long unsigned int.
/// Each orbit gets 10 bits.
//*********************************************************************
uint64_t HartreeFock::Vmon3Hash(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e, uint64_t f)
{
  return a + (b<<10) + (c<<20) + (d<<30) + (e<<40) + (f<<50);
}


//*********************************************************************
/// Take a hashed key and extract the six orbit indices that went into it.
//*********************************************************************
//void HartreeFock::ParseVmon3HashKey(uint64_t key, int& a, int& b, int& c, int& d, int& e, int& f)
void HartreeFock::Vmon3UnHash(uint64_t key, int& a, int& b, int& c, int& d, int& e, int& f)
{
  a = (key    )&0x3FFL; // The L means long, and 0x3FF is 1's in the lowest 10 bits.
  b = (key>>10)&0x3FFL;
  c = (key>>20)&0x3FFL;
  d = (key>>30)&0x3FFL;
  e = (key>>40)&0x3FFL;
  f = (key>>50)&0x3FFL;
}




//*********************************************************************
/// Hashing function for rolling nine indices into a single long unsigned int.
/// Each n value get 5 bits, which means they should be less than 32 (emax=64, probably ok for the near future)
/// Jab gets 6 bits, twoJ gets 7 bits, Lcm gets 7 bits, Ncm gets 6 bits, E12 gets 7 bits
/// leaving 16 bits for alpha, which is hopefully (?) enough.
//*********************************************************************
uint64_t HartreeFock::Vmon3JacobiHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ, uint64_t E12, uint64_t alpha, uint64_t Ncm, uint64_t Lcm)
{
   return (  alpha + (na<<16) + (nb<<21) + (nc<<26) + (Jab<<31) + ((twoJ/2)<<37) + (E12<<44) + (Ncm<<51) + (Lcm<<57) );
}

//*********************************************************************
/// Take a hashed key and extract the six orbit indices that went into it.
//*********************************************************************
void HartreeFock::Vmon3JacobiUnHash(uint64_t key, uint64_t& na, uint64_t& nb, uint64_t& nc, uint64_t& Jab, uint64_t& twoJ, uint64_t& E12, uint64_t& alpha, uint64_t& Ncm, uint64_t& Lcm)
{
   alpha = (key    )&0xFFFF; // 16 bits
   na    = (key>>16)&0x1F;   // 5 bits
   nb    = (key>>21)&0x1F;   // 5 bits
   nc    = (key>>26)&0x1F;   // 5 bits
   Jab   = (key>>31)&0x3F;   // 6 bits
   twoJ  = (key>>37)&0x7F;   // 7 bits
   E12   = (key>>44)&0x7F;   // 7 bits
   Ncm   = (key>>51)&0x3F;   // 6 bits
   Lcm   = (key>>57)&0x7F;   // 7 bits
}




//*********************************************************************
/// one-body density matrix
/// \f$ <i|\rho|j> = \sum\limits_{\beta} n_{\beta} <i|\beta> <\beta|j> \f$
/// where \f$n_{\beta} \f$ ensures that beta runs over HF orbits in
/// the core (i.e. below the fermi surface)
//*********************************************************************
void HartreeFock::UpdateDensityMatrix()
{
  arma::mat tmp = C.cols(holeorbs);
  rho = (tmp.each_row() % hole_occ) * tmp.t();
}


// DIIS: Direct Inversion in the Iterative Subspace
// an approach for accelerating the convergence of the HF iterations
// See P Pulay J. Comp. Chem. 3(4) 556 (1982), and  Garza & Scuseria J. Chem. Phys. 137 054110 (2012)
void HartreeFock::UpdateDensityMatrix_DIIS()
{

  size_t N_MATS_STORED = 5; // How many past density matrices and error matrices to store

  // If we're at the solution, the Fock matrix F and rho will commute
  arma::mat error_mat = F * rho - rho * F;

  // Compute the new density matrix rho from the C matrix which diagonalizes F
  arma::mat tmp = C.cols(holeorbs);
  rho = (tmp.each_row() % hole_occ) * tmp.t();


  // save this error matrix in a list
  size_t nsave = DIIS_error_mats.size();
  if (nsave>N_MATS_STORED) DIIS_error_mats.pop_front();    // out with the old
  DIIS_error_mats.push_back(error_mat);                    // in with the new


  // save the new rho in the list of rhos
  if (DIIS_density_mats.size()>N_MATS_STORED)  DIIS_density_mats.pop_front();   // out with the old
  DIIS_density_mats.push_back( rho );                                           // in with the new

  // Now construct the B matrix which is inverted to find the combination of
  // previous density matrices which will minimize the error

  if (nsave<N_MATS_STORED) return; // check that have enough previous rhos and errors stored to do this

  arma::mat Bij( nsave+1, nsave+1 );
  Bij.row(nsave).ones();
  Bij.col(nsave).ones();
  Bij(nsave,nsave) = 0;

  for (size_t i=0; i<nsave; i++)
  {
   for (size_t j=0; j<nsave; j++)
   {
    Bij(i,j) = arma::norm( DIIS_error_mats[i] * DIIS_error_mats[j].t(), "fro" );
   }
  }
  arma::vec rhs(nsave+1, arma::fill::zeros);
  rhs(nsave) = 1;
  // Now invert that bad boy.
  arma::vec cvec = arma::solve( Bij, rhs );
  cvec = arma::abs(cvec);
  cvec /= arma::accu(cvec.head_rows(nsave));


  rho.zeros();
  // our next guess for rho is a linear combination of the previous ones
  for (size_t i=0; i<nsave; i++)   rho += cvec(i) * DIIS_density_mats[i];

}




//********************************************************************
/// Get new occupation numbers by filling the orbits in order of their
/// single-particle energies. The last proton/neutron orbit can have
/// a fractional filling, corresponding to ensemble normal ordering.
//********************************************************************
void HartreeFock::FillLowestOrbits()
{
  // vector of indices such that they point to elements of F(i,i)
  // in ascending order of energy
  arma::mat F_hfbasis = C.t() * F * C;
  arma::uvec sorted_indices = arma::stable_sort_index( F_hfbasis.diag() );
  double refereceZ = modelspace->GetZref();
  double refereceN = modelspace->GetAref() - refereceZ;
  double placedZ = 0;
  double placedN = 0;
  std::vector<index_t> holeorbs_tmp;
  std::vector<double> hole_occ_tmp;

  for (auto i : sorted_indices)
  {

    Orbit& oi = modelspace->GetOrbit(i);
    if (oi.tz2 < 0 and (placedZ<refereceZ))
    {
      holeorbs_tmp.push_back(i);
      hole_occ_tmp.push_back( std::min(1.0,double(refereceZ-placedZ)/(oi.j2+1.) ) );
      placedZ = std::min(placedZ+oi.j2+1.,refereceZ);
    }
    else if (oi.tz2 > 0 and (placedN<refereceN))
    {
      holeorbs_tmp.push_back(i);
      hole_occ_tmp.push_back( std::min(1.0,double(refereceN-placedN)/(oi.j2+1.) ) );
      placedN = std::min(placedN+oi.j2+1.,refereceN);
    }

    if((placedZ >= refereceZ) and (placedN >= refereceN) ) break;
  }

  std::set<index_t> newholes;
  std::set<index_t> oldholes;

  for (auto i : holeorbs_tmp ) newholes.insert(i);
  for (auto i : holeorbs ) oldholes.insert(i);
  if ( oldholes != newholes )
  {
     std::cout << "Changing hole orbits. New holes:" << std::endl;
     for (auto i : holeorbs_tmp ) std::cout << i << " ";
     std::cout << std::endl;
     std::cout << "                     Old holes:" << std::endl;
     for (auto i : holeorbs ) std::cout << i << " ";
     std::cout << std::endl;
  }

  holeorbs = arma::uvec( holeorbs_tmp );
  hole_occ = arma::rowvec( hole_occ_tmp );
}




//*********************************************************************
///  [See Suhonen eq 4.85]
/// \f[ F_{ij} = t_{ij}  +  \frac{1}{2j_i+1}\sum_{ab} \rho_{ab} \bar{V}^{(2)}_{iajb}
///               + \frac{1}{2(2j_i+1)}\sum_{abcd}\rho_{ab} \rho_{cd} \bar{V}^{(3)}_{iacjbd}  \f]
/// * \f$ F \f$ is the Fock matrix, to be diagonalized
/// * \f$ t \f$ is the kinetic energy
/// * \f$\rho\f$ is the density matrix defined in UpdateDensityMatrix()
/// * \f$ \bar{V}^{(2)} \f$ is the monopole component of the 2-body interaction defined in BuildMonopoleV().
/// * \f$ \bar{V}^{(3)} \f$ is the monopole component of the 3-body interaction devined in BuildMonopoleV3().
//*********************************************************************
void HartreeFock::UpdateF()
{
   double start_time = omp_get_wtime();
//   int norbits = modelspace->GetNumberOrbits();
   Vij.zeros();
   V3ij.zeros();


   // This loop isn't thread safe for some reason. Regardless, parallelizing it makes it slower.
   std::vector<std::array<int,3>> onebodychannels;
   for (auto it : Hbare.OneBodyChannels ) onebodychannels.push_back(it.first);

   #pragma omp parallel for schedule(dynamic,1)
   for (size_t obc=0; obc<onebodychannels.size(); obc++)
   {
    auto& chvec = Hbare.OneBodyChannels.at(onebodychannels[obc]);
   for (index_t i : chvec )
   {
      Orbit& oi = modelspace->GetOrbit(i);
      for (auto j : chvec )
      {
         if (j<i) continue;
         for (auto a : modelspace->all_orbits)
         {
            Orbit& oa = modelspace->GetOrbit(a);
            int Tz = (oi.tz2+oa.tz2)/2;
            int parity = (oi.l+oa.l)%2;
            int bra = modelspace->GetKetIndex(std::min(i,a),std::max(i,a));
            int local_bra = modelspace->MonopoleKets[Tz+1][parity][bra];
            for (auto b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
            {
               int ket = modelspace->GetKetIndex(std::min(j,b),std::max(j,b));
               int local_ket = modelspace->MonopoleKets[Tz+1][parity][ket];
               // 2body term <ai|V|bj>
               if ((a>i) xor (b>j))  // code needed some obfuscation, so threw an xor in there...
               {
                  Vij(i,j) += rho(a,b)*Vmon_exch[Tz+1][parity](local_bra,local_ket); // <a|rho|b> * <ai|Vmon|jb>
               }
               else
               {
                  Vij(i,j) += rho(a,b)*Vmon[Tz+1][parity](local_bra,local_ket); // <a|rho|b> * <ai|Vmon|bj>
               }
           }
         }
      }
      Vij.col(i) /= (oi.j2+1);
   }
   }

   if (Hbare.GetParticleRank()>=3)
   {
      // it's just a one-body matrix, so we can store different copy for each thread.
      std::vector<arma::mat> V3vec(omp_get_max_threads(),V3ij);
      #pragma omp parallel for schedule(dynamic,1)
      for (size_t ind=0;ind<Vmon3.size(); ++ind)
      {
        arma::mat& v3ij = V3vec[omp_get_thread_num()];
        int a,b,c,d,i,j;
        Vmon3UnHash(Vmon3_keys[ind], a,c,i,b,d,j);
        v3ij(i,j) += rho(a,b) * rho(c,d) * Vmon3[ind] ;
      }
      for (auto& v : V3vec) V3ij += v;
      V3vec.clear(); // free up the memory (is this at all necessary, since it goes out of scope?)
   }

   Vij  = arma::symmatu(Vij);
   V3ij = arma::symmatu(V3ij);

   F = KE + Vij + 0.5*V3ij;


   profiler.timer["HF_UpdateF"] += omp_get_wtime() - start_time;
}



//********************************************************
/// Check for convergence using difference in s.p. energies
/// between iterations.
/// Converged when
/// \f[ \delta_{e} \equiv \sqrt{ \sum_{i}(e_{i}^{(n)}-e_{i}^{(n-1)})^2} < \textrm{tolerance} \f]
/// where \f$ e_{i}^{(n)} \f$ is the \f$ i \f$th eigenvalue of the Fock matrix after \f$ n \f$ iterations.
///
//********************************************************
bool HartreeFock::CheckConvergence()
{
   CalcEHF();
   convergence_EHF.push_back(EHF);
   convergence_EHF.pop_front();
   double ediff = arma::norm(energies-prev_energies, "frob") / energies.size();
   convergence_ediff.push_back(ediff); // update list of convergence checks
   convergence_ediff.pop_front();
   return (ediff < tolerance);
}


//**********************************************************************
/// Eigenvectors/values come out of the diagonalization energy-ordered.
/// We want them ordered corresponding to the input ordering, i.e. we want
/// the l,j,tz sub-blockes of the matrix C to be energy-ordered and positive along the diagonal.
/// For a 3x3 matrix this would be something like (this needs to be updated)
/// \f[
/// \left(
/// \begin{array}{rrr}
///   -0.8 & 0.2 & -0.6 \\
///   -0.3 & 0.3 &  0.9 \\
///    0.2 & 0.9 & -0.4 \\
/// \end{array}\right)
/// \rightarrow
/// \left(\begin{array}{rrr}
///    0.8 & -0.6 & 0.2  \\
///    0.3 &  0.9 & 0.3  \\
///   -0.2 & -0.4 & 0.9  \\
/// \end{array}\right)
/// \f]
//**********************************************************************
void HartreeFock::ReorderCoefficients()
{
   for ( auto& it : Hbare.OneBodyChannels )
   {
     arma::uvec orbvec(std::vector<index_t>(it.second.begin(),it.second.end())); // convert from std::set to std::vector, and then to arma::uvec
     arma::vec E_ch = energies(orbvec);
     int nswaps = 10; // keep track of the number of swaps we had to do, iterate until nswaps==0
     // first, make the orbits with the same l,j,tz energy ordered
     while (nswaps>0) // loop until we don't have to make any more swaps
     {
       nswaps = 0;
       for (index_t i=0;i<E_ch.size()-1;i++)
       {
         if (E_ch[i] > E_ch[i+1])
         {
           E_ch.swap_rows(orbvec[i],orbvec[i+1]);
           C.swap_cols(orbvec[i],orbvec[i+1]);
           nswaps++;
         }
       }
      }

     // Make sure the diagonal terms are positive (to avoid confusion later).
     for (index_t i=0;i<C.n_cols;++i) // loop through original basis states
     {
        if (C(i,i) < 0)  C.col(i) *= -1; // C is <HO|HF>, armadillow ordering is C(row,col)  so column j is the jth HF state
     }
   }

}





//**************************************************************************
/// Takes in an operator expressed in the basis of the original Hamiltonian,
/// and returns that operator in the Hartree-Fock basis.
/// \f[ t_{HF} = C^{\dagger} t_{HO} C \f]
/// \f[ V_{HF}^{J} = D^{\dagger} V^{J}_{HO} D \f]
/// The matrix \f$ D \f$ is defined as
/// \f[ D_{ab\alpha\beta} \equiv \sqrt{ \frac{1+\delta_{ab}} {1+\delta_{\alpha\beta}} }  C_{a\alpha} C_{b\beta} \f]
/// The factor in the square root is due to the fact that we're using normalized TBME's.
/// Since only kets with \f$ a\leq b\f$ are stored, we can use the antisymmetry of the TBME's and define
/// \f[ D(J)_{ab\alpha\beta} \equiv \sqrt{ \frac{1+\delta_{ab}} {1+\delta_{\alpha\beta}} }
///      \left( C_{a\alpha} C_{b\beta} -(1-\delta_{ab})(-1)^{j_a+j_b-J} C_{b\alpha}C_{a\beta}\right) \f]
///
//**************************************************************************
Operator HartreeFock::TransformToHFBasis( Operator& OpHO)
{

   Operator OpHF(OpHO);
   // Easy part:
   //Update the one-body part by multiplying by the matrix C(i,a) = <i|a>
   // where |i> is the original basis and |a> is the HF basis
   if ( OpHF.legs%2== 0)
   {
     OpHF.OneBody = C.t() * OpHO.OneBody * C;
   }
   else
   {
     OpHF.OneBody = C.t() * OpHO.OneBody ;
   }

   // Moderately difficult part:
   // Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz. Most of the effort here is in constructing D.

   if ( OpHF.legs%2== 0 and OpHO.legs>3)
   {
     for ( auto& it : OpHO.TwoBody.MatEl )
     {
        int ch_bra = it.first[0];
        int ch_ket = it.first[1];
        TwoBodyChannel& tbc_bra = OpHF.GetModelSpace()->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel& tbc_ket = OpHF.GetModelSpace()->GetTwoBodyChannel(ch_ket);
        int nbras = it.second.n_rows;
        int nkets = it.second.n_cols;
        arma::mat Dbra(nbras,nbras);
        arma::mat Dket(nkets,nkets);
        // loop over all possible original basis configurations <pq| in this J,p,Tz channel.
        // and all possible HF configurations |p'q'> in this J,p,Tz channel
        // bra is in the original basis, ket is in the HF basis
        // i and j are the indices of the matrix D for this channel
        for (int i=0; i<nkets; ++i)
        {
           Ket & ket_ho = tbc_ket.GetKet(i);
           for (int j=0; j<nkets; ++j)
           {
              Ket & ket_hf = tbc_ket.GetKet(j);
              Dket(i,j) = C(ket_ho.p,ket_hf.p) * C(ket_ho.q,ket_hf.q);
              if (ket_ho.p!=ket_ho.q)
              {
                 Dket(i,j) += C(ket_ho.q, ket_hf.p) * C(ket_ho.p, ket_hf.q) * ket_ho.Phase(tbc_ket.J);
              }
              if (ket_ho.p==ket_ho.q)    Dket(i,j) *= PhysConst::SQRT2;
              if (ket_hf.p==ket_hf.q)    Dket(i,j) /= PhysConst::SQRT2;
           }
        }
        if (ch_bra == ch_ket)
        {
          Dbra = Dket.t();
        }
        else
        {
          for (int i=0; i<nbras; ++i)
          {
             Ket & bra_hf = tbc_bra.GetKet(i);
             for (int j=0; j<nbras; ++j)
             {
                Ket & bra_ho = tbc_bra.GetKet(j);
                Dbra(i,j) = C(bra_ho.p,bra_hf.p) * C(bra_ho.q,bra_hf.q);
                if (bra_ho.p!=bra_ho.q)
                {
                   Dbra(i,j) += C(bra_ho.q, bra_hf.p) * C(bra_ho.p, bra_hf.q) * bra_ho.Phase(tbc_bra.J);
                }
                if (bra_ho.p==bra_ho.q)    Dbra(i,j) *= PhysConst::SQRT2;
                if (bra_hf.p==bra_hf.q)    Dbra(i,j) /= PhysConst::SQRT2;
             }
          }
        }
        auto& IN  =  it.second;
        auto& OUT =  OpHF.TwoBody.GetMatrix(ch_bra,ch_ket);
        OUT  =    Dbra * IN * Dket;

     }
   }
   else
   {
     for ( auto& it : OpHO.ThreeLeg.MatEl )
     {
        int ch_bra = it.first;
        TwoBodyChannel& tbc_bra = OpHF.GetModelSpace()->GetTwoBodyChannel(ch_bra);
//        int nbras = it.second.n_rows;
        int nbras = it.second.n_rows;
        arma::mat Dbra(nbras,nbras);
        for (int i=0; i<nbras; ++i)
        {
           Ket & bra_hf = tbc_bra.GetKet(i);
           for (int j=0; j<nbras; ++j)
           {
              Ket & bra_ho = tbc_bra.GetKet(j);
              Dbra(i,j) = C(bra_ho.p,bra_hf.p) * C(bra_ho.q,bra_hf.q);
              if (bra_ho.p!=bra_ho.q)
              {
                 Dbra(i,j) += C(bra_ho.q, bra_hf.p) * C(bra_ho.p, bra_hf.q) * bra_ho.Phase(tbc_bra.J);
              }
              if (bra_ho.p==bra_ho.q)    Dbra(i,j) *= PhysConst::SQRT2;
              if (bra_hf.p==bra_hf.q)    Dbra(i,j) /= PhysConst::SQRT2;
           }
        }
        auto& IN  =  it.second;
        auto& OUT =  OpHF.ThreeLeg.GetMatrix(ch_bra);
        OUT  =    Dbra * IN * C;
     }
   }

   if ( OpHF.GetParticleRank() >= 3 )
   {
       OpHF.ThreeBody = GetTransformed3B( OpHO, C );
   }

   return OpHF;
}

//**************************************************************************
/// If the lowest orbits are different from our previous guess, we should update the reference.
//**************************************************************************
void HartreeFock::UpdateReference()
{

  bool changed_occupations = false;
  std::map<index_t,double> hole_map;
  for (index_t i=0;i<holeorbs.size();++i)
  {
     hole_map[holeorbs[i]] = hole_occ[i];
     if ( std::abs(modelspace->GetOrbit( holeorbs[i] ).occ - hole_occ[i]) > 1e-3)
     {
        changed_occupations = true;
        std::cout << "After HF, occupation of orbit " << holeorbs[i] << " has changed. Modelspace will be updated." << std::endl;
     }
  }
  if (changed_occupations)  modelspace->SetReference( hole_map );

}


//**************************************************************************
/// Returns the Hamiltonian normal ordered with respect to the reference
/// in the basis defined by Cin, where Cin is a transformation from the
/// original basis (typically the harmonic oscillator).
//**************************************************************************
Operator HartreeFock::GetNormalOrderedH(arma::mat& Cin, int particle_rank)
{
  C=Cin;
//  ReorderCoefficients();  // Reorder columns of C so we can properly identify the hole orbits.
  UpdateDensityMatrix();  // Update the 1 body density matrix, used in UpdateF()
  UpdateF();              // Update the Fock matrix
  CalcEHF();
  PrintEHF();

  return GetNormalOrderedH(particle_rank);
}

//**************************************************************************
/// Returns the normal-ordered Hamiltonian in the Hartree-Fock basis, neglecting the residual 3-body piece.
/// \f[ E_0 = E_{HF} \f]
/// \f[ f = C^{\dagger} F C \f]
/// \f[ \Gamma = D^{\dagger} \left(V^{(2)}+V^{(3\rightarrow 2)} \right) D \f]
/// \f[ V^{(2\rightarrow 3)J}_{ijkl} \equiv \frac{1}{\sqrt{(1+\delta_{ij})(1+\delta_{kl})}}\sum_{ab}\sum_{J_3}(2J_{3}+1)\rho_{ab}V^{JJJ_{3}}_{ijaklb} \f]
/// Where \f$ F\f$ is the Fock matrix obtained in UpdateF() and the matrix \f$ D\f$ is the same as the one defined in TransformToHFBasis().
///
//**************************************************************************
Operator HartreeFock::GetNormalOrderedH(int particle_rank)
{
   double start_time = omp_get_wtime();

   // First, check if we need to update the occupation numbers for the reference
   
   if (not freeze_occupations)
   {
     UpdateReference();
   }


   Operator HNO = Operator(*modelspace,0,0,0,particle_rank);
   HNO.ZeroBody = EHF;
   HNO.OneBody = C.t() * F * C;

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
     TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
     int J = tbc.J;
     int npq = tbc.GetNumberKets();

     arma::mat D(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
     arma::mat V3NO(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>

#pragma omp parallel for schedule(dynamic,1) // confirmed that this improves performance
     for (int i=0; i<npq; ++i)
     {
       Ket & bra = tbc.GetKet(i);
       int e2bra = 2*bra.op->n + bra.op->l + 2*bra.oq->n + bra.oq->l;
       for (int j=0; j<npq; ++j)
       {
         Ket & ket = tbc.GetKet(j);
         int e2ket = 2*ket.op->n + ket.op->l + 2*ket.oq->n + ket.oq->l;
         D(i,j) = C(bra.p,ket.p) * C(bra.q,ket.q);
         if (bra.p!=bra.q)
         {
           D(i,j) += C(bra.q,ket.p) * C(bra.p,ket.q) * bra.Phase(J);
         }
         if (bra.p==bra.q)    D(i,j) *= PhysConst::SQRT2;
         if (ket.p==ket.q)    D(i,j) /= PhysConst::SQRT2;

         // Now generate the NO2B part of the 3N interaction
         if (Hbare.GetParticleRank()<3) continue;
         if (i>j) continue;
         //            for (int a=0; a<norb; ++a)
         for ( auto a : modelspace->all_orbits )
         {
           Orbit & oa = modelspace->GetOrbit(a);
           if ( 2*oa.n+oa.l+e2bra > Hbare.GetE3max() ) continue;
           for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
           {
             Orbit & ob = modelspace->GetOrbit(b);
             if ( 2*ob.n+ob.l+e2ket > Hbare.GetE3max() ) continue;
             if ( std::abs(rho(a,b)) < 1e-8 ) continue; // Turns out this helps a bit (factor of 5 speed up in tests)
             V3NO(i,j) += rho(a,b) * Hbare.ThreeBody.GetME_pn_no2b(bra.p,bra.q,a, ket.p,ket.q,b, J);
           }
         }
         V3NO(i,j) /= (2*J+1);
         if (bra.p==bra.q)  V3NO(i,j) /= PhysConst::SQRT2;
         if (ket.p==ket.q)  V3NO(i,j) /= PhysConst::SQRT2;
         V3NO(j,i) = V3NO(i,j);
       }
     }

     // This option is just for running diagnostics. No reason to use it for production runs
     if (discard_NO2B_from_3N)
     {
        std::cout << " ====  Discarding NO2B contribution from 3N. Is that what you want? ===== " << std::endl;
        V3NO *= 0;
     }

     auto& V2  =  Hbare.TwoBody.GetMatrix(ch);
     auto& OUT =  HNO.TwoBody.GetMatrix(ch);
     OUT  =    D.t() * (V2 + V3NO) * D;
   }// for ch

   if (particle_rank>2)
   {
//     HNO.ThreeBody = GetTransformed3B( Hbare );
     HNO.ThreeBody = GetTransformed3B( Hbare, C );
   }


   profiler.timer["HF_GetNormalOrderedH"] += omp_get_wtime() - start_time;

   return HNO;

}


//**************************************************************************
/// Not sure if this is actually helpful. Memory management is a mystery to me.
//**************************************************************************
void HartreeFock::FreeVmon()
{
   // free up some memory
   std::array< std::array< arma::mat,2>,3>().swap(Vmon);
   std::array< std::array< arma::mat,2>,3>().swap(Vmon_exch);
   std::vector< double>().swap( Vmon3 );
   std::vector< uint64_t>().swap( Vmon3_keys );
}


//**************************************************************************
/// Get the one-body generator corresponding to the transformation to the HF basis.
/// Since the unitary transformation for HF is given by the \f$ U_{HF} = C^{\dagger} \f$ matrix, we have
/// \f$ e^{-\Omega} = C \Rightarrow \Omega = -\log(C) \f$.
/// The log is evaluated by diagonalizing the one-body submatrix and taking the log of the diagonal entries.
/// This is much slower than the other methods, but it might be useful for debugging/benchmarking.
//**************************************************************************
Operator HartreeFock::GetOmega()
{
   Operator Omega(*modelspace,0,0,0,1);
   Omega.SetAntiHermitian();
   for ( auto& it : Hbare.OneBodyChannels )
   {
      arma::uvec orbvec(std::vector<index_t>(it.second.begin(),it.second.end()));
      arma::mat C_ch = C.submat(orbvec,orbvec);
      arma::cx_mat eigvect;
      arma::cx_vec eigval;
      arma::eig_gen(eigval, eigvect, C_ch);
      Omega.OneBody.submat(orbvec,orbvec) = -arma::real( eigvect * arma::diagmat(arma::log(eigval)) * eigvect.t()) ;
   }
   return Omega;
}


//**************************************************************************
/// Print out the single particle orbits with their energies.
//**************************************************************************
void HartreeFock::PrintSPE()
{
  arma::mat F_hfbasis = C.t() * F * C;
  for ( auto i : modelspace->all_orbits )
  {
    Orbit& oi = modelspace->GetOrbit(i);
    std::cout << std::fixed << std::setw(3) << oi.n << " " << std::setw(3) << oi.l << " "
         << std::setw(3) << oi.j2 << " " << std::setw(3) << oi.tz2 << "   " << std::setw(10) << F_hfbasis(i,i) << std::endl;
  }

}


//**************************************************************************
/// Print out the single particle orbits with their energies, as well as 
/// the overlap with the input (oscillator) basis states
//**************************************************************************
void HartreeFock::PrintSPEandWF()
{
  arma::mat F_hfbasis = C.t() * F * C;
  std::cout << std::fixed << std::setw(3) << "i" << ": " << std::setw(3) << "n" << " " << std::setw(3) << "l" << " "
       << std::setw(3) << "2j" << " " << std::setw(3) << "2tz" << "   " << std::setw(12) << "SPE" << " " << std::setw(12) << "occ." << "   |   " << " overlaps" << std::endl;

  for ( auto i : modelspace->all_orbits )
  {
    Orbit& oi = modelspace->GetOrbit(i);
    std::cout << std::fixed << std::setw(3) << i << ": " << std::setw(3) << oi.n << " " << std::setw(3) << oi.l << " "
         << std::setw(3) << oi.j2 << " " << std::setw(3) << oi.tz2 << "   " << std::setw(12) << std::setprecision(6) << F_hfbasis(i,i) << " " << std::setw(12) << oi.occ << "   | ";
    for (int j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) // j loops over harmonic oscillator states
    {
      std::cout << std::setw(9) << C(j,i) << "  ";  // C is <HO|HF>
    }
    std::cout << std::endl;
  }
}


void HartreeFock::GetRadialWF(index_t index, std::vector<double>& R, std::vector<double>& PSI)
{
  PSI.resize(0);
  for (double r : R)
  {
    PSI.push_back( GetRadialWF_r(index, r) );
  }
}


//**************************************************************************
/// Get the radial wave function of an orbit in the HF basis at point R.
/// It is assumed that the input basis is the harmonic oscillator.
//**************************************************************************
double HartreeFock::GetRadialWF_r(index_t index, double R)
{
  using PhysConst::M_NUCLEON;
  using PhysConst::HBARC;
  using PhysConst::SQRTPI;
  double b = sqrt( (HBARC*HBARC) / (modelspace->GetHbarOmega() * M_NUCLEON) );
  Orbit& orb = modelspace->GetOrbit(index);
   double x = R/b;
   double psi = 0;
   for ( index_t j : Hbare.OneBodyChannels.at({orb.l,orb.j2,orb.tz2}) )
   {
     Orbit& oj = modelspace->GetOrbit(j);
     double Norm = 2*sqrt( AngMom::factorial(oj.n) * pow(2,oj.n+oj.l) / SQRTPI / AngMom::double_fact(2*oj.n+2*oj.l+1) * pow(b,-3.0) );
     psi += C(index,j) * Norm * pow(x,oj.l) * exp(-x*x*0.5) * gsl_sf_laguerre_n(oj.n,oj.l+0.5,x*x);
   }
   return psi;
}





//**************************************************************************
/// The effective one-body Hartree-Fock potential for an orbit with quantum
/// numbers l,j,tz. This is a non-local potential, i.e. it depends on r and r'.
//**************************************************************************
double HartreeFock::GetHFPotential( size_t i, double r, double rprime)
{
   double Uhf = 0;
   Orbit& oi = modelspace->GetOrbit(i);
   for (auto j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
   {
     for (auto k : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
     {
       Uhf +=  (Vij(j,k) + 0.5*V3ij(j,k) )  * GetRadialWF_r(j,r) * GetRadialWF_r(k,rprime);
     }
   }
   return Uhf;
}


double HartreeFock::GetAverageHFPotential( double r, double rprime)
{
  double Uhf = 0;
  size_t norb = modelspace->GetNumberOrbits();
  for (size_t i=0; i<norb; i++)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    if (oi.n == 0)
    {
      Uhf += GetHFPotential(i,r,rprime);
    }
  }
  Uhf /= 2*(2*modelspace->GetEmax()+1); // 2 for protons/neutrons, and 2*emax+1 is number of different l,j combinations
  return Uhf;
}


//**************************************************************************
/// In most standard calculations, we take the NO2B approximation after the HF
/// step, so we don't need this. However, if we want to use the residual
/// 3-body in any post-HF steps, we need to transform it to the HF basis.
//**************************************************************************
//ThreeBodyME HartreeFock::GetTransformed3B( Operator& OpIn )
ThreeBodyME HartreeFock::GetTransformed3B( Operator& OpIn, arma::mat& C_3b )
{
  double t_start = omp_get_wtime();

//  ThreeBodyME hf3bme(modelspace, OpIn.GetJRank(), OpIn.GetTRank(), OpIn.GetParity());
  ThreeBodyME hf3bme(ms_for_output_3N, OpIn.GetJRank(), OpIn.GetTRank(), OpIn.GetParity());
  hf3bme.SwitchToPN_and_discard();
//  std::cout << __func__ << "  Norm of input 3bme in oscillator basis = " << OpIn.ThreeBodyNorm() << std::endl;

  std::map<int,double> e_fermi = ms_for_output_3N->GetEFermi();

  std::vector<size_t> chbra_list,chket_list;
  for ( auto& itmat : hf3bme.Get_ch_start() )
  {
    chbra_list.push_back( itmat.first.ch_bra );
    chket_list.push_back( itmat.first.ch_ket );
  }
  size_t nch3 = chbra_list.size();
  #pragma omp parallel for schedule(dynamic,1)
  for ( size_t ich=0; ich<nch3; ich++ )
  {
    size_t ch_bra = chbra_list[ich];
    size_t ch_ket = chket_list[ich];
    ThreeBodyChannel& Tbc_bra_HO = modelspace->GetThreeBodyChannel(ch_bra);
    ThreeBodyChannel& Tbc_ket_HO = modelspace->GetThreeBodyChannel(ch_ket);
    ThreeBodyChannel& Tbc_bra_HF = ms_for_output_3N->GetThreeBodyChannel(ch_bra);
    ThreeBodyChannel& Tbc_ket_HF = ms_for_output_3N->GetThreeBodyChannel(ch_ket);
    int twoJ = Tbc_bra_HO.twoJ;
    size_t nbras_HO = Tbc_bra_HO.GetNumberKets();
    size_t nkets_HO = Tbc_ket_HO.GetNumberKets();
    size_t nbras_HF = Tbc_bra_HF.GetNumberKets();
    size_t nkets_HF = Tbc_ket_HF.GetNumberKets();
    std::vector<size_t> bras_kept;
    std::vector<size_t> kets_kept;
    std::map<size_t,size_t> kept_lookup_bra;
    std::map<size_t,size_t> kept_lookup_ket;


    for ( size_t ibra=0; ibra<nbras_HF; ibra++)
    {
      Ket3& bra = Tbc_bra_HF.GetKet(ibra);
      int ei = 2*bra.op->n + bra.op->l;
      int ej = 2*bra.oq->n + bra.op->l;
      int ek = 2*bra.oR->n + bra.op->l;
      int tz2i = bra.op->tz2;
      int tz2j = bra.oq->tz2;
      int tz2k = bra.oR->tz2;
      if (  ( std::abs(ei - e_fermi[tz2i]) + std::abs(ej-e_fermi[tz2j]) + std::abs(ek-e_fermi[tz2k])) > ms_for_output_3N->GetdE3max() ) continue;
      bras_kept.push_back( ibra );
      kept_lookup_bra[ibra] = bras_kept.size()-1;
    }
    size_t nbras_kept = bras_kept.size();


    if ( ch_bra==ch_ket )
    {
      kets_kept = bras_kept;
    }
    else
    {
      for ( size_t iket=0; iket<nkets_HF; iket++)
      {
        Ket3& ket = Tbc_ket_HF.GetKet(iket);
        int ei = 2*ket.op->n + ket.op->l;
        int ej = 2*ket.oq->n + ket.op->l;
        int ek = 2*ket.oR->n + ket.op->l;
        int tz2i = ket.op->tz2;
        int tz2j = ket.oq->tz2;
        int tz2k = ket.oR->tz2;
        if (  ( std::abs(ei - e_fermi[tz2i]) + std::abs(ej-e_fermi[tz2j]) + std::abs(ek-e_fermi[tz2k])) > ms_for_output_3N->GetdE3max() ) continue;
        kets_kept.push_back( iket );
        kept_lookup_ket[iket] = kets_kept.size()-1;
      }
    }
    size_t nkets_kept = kets_kept.size();



   arma::mat Dbra(nbras_kept, nbras_HO, arma::fill::zeros );
   arma::mat Dket(nkets_HO, nkets_kept, arma::fill::zeros );
//   arma::mat Dbra(nbras_kept, nbras_kept, arma::fill::zeros );
//   arma::mat Dket(nkets_kept, nkets_kept, arma::fill::zeros );
// Potentially worth exploring if any improvement is made by using sparse matrices
//   arma::sp_mat Dbra(nbras_kept, nbras_kept );
//   arma::sp_mat Dket(nkets_kept, nkets_kept );

//   arma::mat Vho( nbras_kept, nkets_kept, arma::fill::zeros );
   arma::mat Vho( nbras_HO, nkets_HO, arma::fill::zeros );


    for ( size_t indxHF=0; indxHF<nbras_kept; indxHF++ )
    {
       size_t ibra_HF = bras_kept[indxHF];
       Ket3& bra_HF = Tbc_bra_HF.GetKet(ibra_HF);
       // Get the index of the HF orbit in the original modelspace
       size_t iHF = modelspace->GetOrbitIndex( bra_HF.op->n, bra_HF.op->l, bra_HF.op->j2, bra_HF.op->tz2 );
       size_t jHF = modelspace->GetOrbitIndex( bra_HF.oq->n, bra_HF.oq->l, bra_HF.oq->j2, bra_HF.oq->tz2 );
       size_t kHF = modelspace->GetOrbitIndex( bra_HF.oR->n, bra_HF.oR->l, bra_HF.oR->j2, bra_HF.oR->tz2 );
//       size_t iHF = bra_HF.p;
//       size_t jHF = bra_HF.q;
//       size_t kHF = bra_HF.r;
       int JijHF = bra_HF.Jpq;
       for ( size_t indxHO=0; indxHO<nbras_HO; indxHO++ )
       {
//         size_t ibra_HO = bras_kept[indxHO];
//         Ket3& bra_HO = Tbc_bra.GetKet(ibra_HO);
         Ket3& bra_HO = Tbc_bra_HO.GetKet(indxHO);
         size_t iHO = bra_HO.p;
         size_t jHO = bra_HO.q;
         size_t kHO = bra_HO.r;
         Orbit& oiHO = modelspace->GetOrbit(iHO);
         Orbit& ojHO = modelspace->GetOrbit(jHO);
         Orbit& okHO = modelspace->GetOrbit(kHO);
         double ji = oiHO.j2 *0.5;
         double jj = ojHO.j2 *0.5;
         double jk = okHO.j2 *0.5;
         int JijHO = bra_HO.Jpq;
         double dbra = 0;

         for ( auto perm3b : OpIn.ThreeBody.UniquePermutations( iHO, jHO, kHO ) )
         {
           size_t iiHO,jjHO,kkHO;
           OpIn.ThreeBody.Permute(perm3b, iHO,jHO,kHO, iiHO,jjHO,kkHO);
           double overlap = C_3b(iiHO,iHF) * C_3b(jjHO,jHF) * C_3b(kkHO,kHF);
           if ( std::abs( overlap) > 1e-8 )
           {
             double phase = OpIn.ThreeBody.PermutationPhase(perm3b); // the fermionic sign from the permutation
             double recouple = OpIn.ThreeBody.RecouplingCoefficient( perm3b, ji,jj,jk, JijHF, JijHO, twoJ); // the angular momentum recoupling coefficient
             dbra += phase * recouple * overlap;
           }
         }

         Dbra( indxHF, indxHO ) = dbra ;
      }// for indxHO
    }// for indxHF

    if (ch_bra==ch_ket)
    {
      Dket = Dbra.t();
    }
    else
    {
//      for ( size_t indxHO=0; indxHO<nkets_kept; indxHO++ )
      for ( size_t indxHO=0; indxHO<nkets_HO; indxHO++ )
      {
//        size_t iket_HO = kets_kept[indxHO];
//        Ket3& ket_HO = Tbc_ket.GetKet(iket_HO);
        Ket3& ket_HO = Tbc_ket_HO.GetKet(indxHO);
        size_t iHO = ket_HO.p;
        size_t jHO = ket_HO.q;
        size_t kHO = ket_HO.r;
        int JijHO = ket_HO.Jpq;
        Orbit& oiHO = modelspace->GetOrbit(iHO);
        Orbit& ojHO = modelspace->GetOrbit(jHO);
        Orbit& okHO = modelspace->GetOrbit(kHO);
        double ji = oiHO.j2 *0.5;
        double jj = ojHO.j2 *0.5;
        double jk = okHO.j2 *0.5;

        
        for ( size_t indxHF=0; indxHF<nkets_kept; indxHF++ )
        {
           size_t iket_HF = kets_kept[indxHF];
           Ket3& ket_HF = Tbc_ket_HF.GetKet(iket_HF);

           size_t iHF = modelspace->GetOrbitIndex( ket_HF.op->n, ket_HF.op->l, ket_HF.op->j2, ket_HF.op->tz2 );
           size_t jHF = modelspace->GetOrbitIndex( ket_HF.oq->n, ket_HF.oq->l, ket_HF.oq->j2, ket_HF.oq->tz2 );
           size_t kHF = modelspace->GetOrbitIndex( ket_HF.oR->n, ket_HF.oR->l, ket_HF.oR->j2, ket_HF.oR->tz2 );
//           size_t iHF = ket_HF.p;
//           size_t jHF = ket_HF.q;
//           size_t kHF = ket_HF.r;
           int JijHF = ket_HF.Jpq;

           double dket = 0;

           for ( auto perm3b : OpIn.ThreeBody.UniquePermutations( iHO, jHO, kHO ) )
           {
             size_t iiHO,jjHO,kkHO;
             OpIn.ThreeBody.Permute(perm3b, iHO,jHO,kHO, iiHO,jjHO,kkHO);
             double overlap = C(iiHO,iHF) * C(jjHO,jHF) * C(kkHO,kHF);
             if ( std::abs( overlap) > 1e-8 )
             {
               double phase = OpIn.ThreeBody.PermutationPhase(perm3b); // the fermionic sign from the permutation
               double recouple = OpIn.ThreeBody.RecouplingCoefficient( perm3b, ji,jj,jk, JijHF, JijHO, twoJ); // the angular momentum recoupling coefficient
               dket += phase * recouple * overlap;
             }
           }

           Dket( indxHO, indxHF ) = dket ;

        }// for indxHF
      }// for indxHO
    }// else ch_bra != ch_ket



//    for ( size_t indx_bra=0; indx_bra<nbras_kept; indx_bra++ )
    for ( size_t indx_bra=0; indx_bra<nbras_HO; indx_bra++ )
    {
//      size_t ibra_HO = bras_kept[indx_bra];
//      Ket3& bra_HO = Tbc_bra.GetKet(ibra_HO);
      Ket3& bra_HO = Tbc_bra_HO.GetKet(indx_bra);
      size_t iHO = bra_HO.p;
      size_t jHO = bra_HO.q;
      size_t kHO = bra_HO.r;
      int JijHO = bra_HO.Jpq;

//      for ( size_t indx_ket=0; indx_ket<nkets_kept; indx_ket++ )
      for ( size_t indx_ket=0; indx_ket<nkets_HO; indx_ket++ )
      {
//        size_t iket_HO = kets_kept[indx_ket];
//        Ket3& ket_HO = Tbc_ket.GetKet(iket_HO);
        Ket3& ket_HO = Tbc_ket_HO.GetKet(indx_ket);
        size_t lHO = ket_HO.p;
        size_t mHO = ket_HO.q;
        size_t nHO = ket_HO.r;
        int JlmHO = ket_HO.Jpq;
        Vho(indx_bra,indx_ket) = OpIn.ThreeBody.GetME_pn(JijHO, JlmHO, twoJ, iHO, jHO, kHO, lHO, mHO, nHO );
      }
    }

    arma::mat Vhf = Dbra * Vho * Dket;

    for ( size_t indx_bra=0; indx_bra<nbras_kept; indx_bra++ )
    {
      size_t ibra_HF = bras_kept[indx_bra];
//      Ket3& bra_HF = Tbc_bra.GetKet(ibra_HF);
//      size_t iHF = bra_HF.p;
//      size_t jHF = bra_HF.q;
//      size_t kHF = bra_HF.r;
//      int JijHF = bra_HF.Jpq;

      for ( size_t indx_ket=0; indx_ket<nkets_kept; indx_ket++ )
      {
        size_t iket_HF = kets_kept[indx_ket];
//        Ket3& ket_HF = Tbc_ket.GetKet(iket_HF);
//        size_t lHF = ket_HF.p;
//        size_t mHF = ket_HF.q;
//        size_t nHF = ket_HF.r;
//        int JlmHF = ket_HF.Jpq;
        double VHF = Vhf(indx_bra, indx_ket);
        if (std::abs(VHF)<1e-9) continue;

        hf3bme.SetME_pn_ch( ch_bra,ch_ket, ibra_HF,iket_HF, VHF);

      }
    }




    // the slow and easy way...
 /*
//    for ( size_t ibra : kets_kept )
    for ( size_t ibra : bras_kept )
    {
      Ket3& bra = Tbc_bra.GetKet(ibra);
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      int Jij = bra.Jpq;
      for ( size_t iket : kets_kept )
      {
//        if (iket<ibra) continue;
        if ((ch_bra==ch_ket) and (iket<ibra)) continue;
        Ket3& ket = Tbc_ket.GetKet(iket);
        size_t l = ket.p;
        size_t m = ket.q;
        size_t n = ket.r;
        int Jlm = ket.Jpq;
//        double VHO = GetHF3bme(  Jij,  Jlm, twoJ,  i,j,k,l,m,n);
//        double VHO = GetTransformed3bme(  Jij,  Jlm, twoJ,  i,j,k,l,m,n);
        double VHF = GetTransformed3bme( OpIn, Jij,  Jlm, twoJ,  i,j,k,l,m,n);
//        if (ich < 3)
//        {        
//        double VHO = OpIn.ThreeBody.GetME_pn( Jij,  Jlm, twoJ,  i,j,k,l,m,n);
//          std::cout << " ich = " << ich << " ibra,iket " << ibra << " " << iket << "   VHO VHF = " << VHO << " " << VHF
//                    << "  ijklmn = " << i << " " << j << " " << k << " " << l << " " << m << " " << n << "   "
//                    << "Jij Jlm twoJ = " << Jij << " " << Jlm << " " << twoJ << std::endl;
//        }
////        hf3bme.SetME_pn_PN_ch( ch3,ch3, ibra,iket, VHO);
////        hf3bme.SetME_pn_PN_ch( ch_bra,ch_ket, ibra,iket, VHO);
        hf3bme.SetME_pn_ch( ch_bra,ch_ket, ibra,iket, VHF);
////        hf3bme.SetME_pn_PN_ch( ch_bra,ch_ket, ibra,iket, VHF);
      }
    }
*/


  }// for ich
//  std::cout << __func__ << "   After transforming norm is " << hf3bme.Norm() << std::endl;

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return hf3bme;

}

/// Set a different modelspace for the output 3N upon transforming to the HF basis.
/// This might be useful if we do the HF step with a large 3N space, and then want
/// to do IMSRG(3) in a smaller 3N space. This way we can avoid transforming 3N matrix elements
/// that we then immediately discard.
void HartreeFock::SetModelspaceForOutput3N( ModelSpace& ms )
{
   ms_for_output_3N = &ms;
}



// Get the 3-body interaction in the Hartree-Fock basis (valence space only).
ThreeBodyME HartreeFock::GetValence3B( Operator& OpIn, int emax, int E3max )
{
  double t_start = omp_get_wtime();
  ThreeBodyME hf3bme(modelspace, E3max);
  hf3bme.Setemax(emax);
  hf3bme.SwitchToPN_and_discard();

  // big loop over elements of hf3bme...
//  auto norbits = modelspace->GetNumberOrbits();
  size_t nch3 = modelspace->GetNumberThreeBodyChannels();
  for ( size_t ch3=0; ch3<nch3; ch3++)
  {
    ThreeBodyChannel& Tbc = modelspace->GetThreeBodyChannel(ch3);
    int twoJ = Tbc.twoJ;
    size_t nkets = Tbc.GetNumberKets();
    for (size_t ibra=0; ibra<nkets; ibra++)
    {
      Ket3& bra = Tbc.GetKet(ibra);
      if ( (bra.op->cvq != 1) or (bra.oq->cvq != 1) or (bra.oR->cvq!=1) ) continue;
      size_t a = bra.p;
      size_t b = bra.q;
      size_t c = bra.r;
      int Jab = bra.Jpq;
      for (size_t iket=0; iket<nkets; iket++)
      {
        Ket3& ket = Tbc.GetKet(iket);
        if ( (ket.op->cvq != 1) or (ket.oq->cvq != 1) or (ket.oR->cvq!=1) ) continue;
        size_t d = ket.p;
        size_t e = ket.q;
        size_t f = ket.r;
        int Jde = ket.Jpq;
        double V = GetTransformed3bme(OpIn, Jab, Jde, twoJ, a, b, c, d, e, f );
        hf3bme.SetME_pn_ch(ch3,ch3,ibra,iket, V);
      }// for iket
    }// for ibra
  }// for ch3
  IMSRGProfiler::timer["HartreeFock::GetValence3B"] += omp_get_wtime() - t_start;
  return hf3bme;

}




// Get a single 3-body matrix element in the HartreeFock basis.
// This is the straightforward but inefficient way to do it.
// Good for debugging/benchmarking.
double HartreeFock::GetTransformed3bme( Operator& OpIn, int Jab, int Jde, int J2,  size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  double V_hf = 0.;
//  double debugV = 0.;
  Orbit& oa = modelspace->GetOrbit(a);
  Orbit& ob = modelspace->GetOrbit(b);
  Orbit& oc = modelspace->GetOrbit(c);
  Orbit& od = modelspace->GetOrbit(d);
  Orbit& oe = modelspace->GetOrbit(e);
  Orbit& of = modelspace->GetOrbit(f);

  for (auto alpha : modelspace->OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
  {
   if ( std::abs(C(alpha,a)) < 1e-8 ) continue;
   for (auto beta : modelspace->OneBodyChannels.at({ob.l,ob.j2,ob.tz2}) )
   {
    if ( std::abs(C(beta,b)) < 1e-8 ) continue;
    for (auto gamma : modelspace->OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) )
    {
     if ( std::abs(C(gamma,c)) < 1e-8 ) continue;
     for (auto delta : modelspace->OneBodyChannels.at({od.l,od.j2,od.tz2}) )
     {
      if ( std::abs(C(delta,d)) < 1e-8 ) continue;
      for (auto epsilon : modelspace->OneBodyChannels.at({oe.l,oe.j2,oe.tz2}) )
      {
       if ( std::abs(C(epsilon,e)) < 1e-8 ) continue;
       for (auto phi : modelspace->OneBodyChannels.at({of.l,of.j2,of.tz2}) )
       {
         double V_ho = OpIn.ThreeBody.GetME_pn( Jab,  Jde,  J2,  alpha,  beta,  gamma,  delta,  epsilon,  phi);
         V_hf += V_ho * C(alpha,a) * C(beta,b) * C(gamma,c) * C(delta,d) * C(epsilon,e) * C(phi,f);

       } // for phi
      } // for epsilon
     } // for delta
    } // for gamma
   } // for beta
  } // for alpha
  return V_hf;
}


