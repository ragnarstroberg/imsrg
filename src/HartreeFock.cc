
#include "HartreeFock.hh"
#include "ModelSpace.hh"
#include <iomanip>
#include <vector>
#include <array>
#include <map>
#include <utility> // for make_pair
#include "gsl/gsl_sf_gamma.h" // for radial wave function
#include "gsl/gsl_sf_laguerre.h" // for radial wave function
#include <gsl/gsl_math.h> // for M_SQRTPI
#include <omp.h>

#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif
#define HBARC 197.3269718 // hc in MeV * fm
#define M_NUCLEON 938.9185 // average nucleon mass in MeV


//using namespace std;

HartreeFock::HartreeFock(Operator& hbare)
  : Hbare(hbare), modelspace(hbare.GetModelSpace()), 
    KE(Hbare.OneBody), energies(Hbare.OneBody.diag()),
    tolerance(1e-8), convergence_ediff(7,0), convergence_EHF(7,0), freeze_occupations(true)
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
       Vmon[Tz+1][parity] = arma::mat(nKetsMon,nKetsMon);
       Vmon_exch[Tz+1][parity] = arma::mat(nKetsMon,nKetsMon);
     }
   }
   prev_energies = arma::vec(norbits,arma::fill::zeros);
   std::vector<double> occvec;
   for (auto& h : modelspace->holes) occvec.push_back(modelspace->GetOrbit(h).occ);
   holeorbs = arma::uvec(modelspace->holes);
   hole_occ = arma::rowvec(occvec);
   BuildMonopoleV();
   if (hbare.GetParticleRank()>2)
   {
      BuildMonopoleV3();
   }
   UpdateDensityMatrix();
   UpdateF();

}


//*********************************************************************
/// Diagonalize and update the Fock matrix until convergence.
/// Then, call ReorderCoefficients() to make sure the index
/// ordering and phases are preserved in the transformation
/// from the original basis to the Hatree-Fock basis.
//*********************************************************************
void HartreeFock::Solve()
{
   iterations = 0; // counter so we don't go on forever
   int maxiter = 1000;

   for (iterations=0; iterations<maxiter; ++iterations)
   {
      Diagonalize();          // Diagonalize the Fock matrix
      ReorderCoefficients();  // Reorder columns of C so we can properly identify the hole orbits.
      if (not freeze_occupations) FillLowestOrbits(); // if we don't freeze the occupations, then calculate the new ones.
      UpdateDensityMatrix();  // Update the 1 body density matrix, used in UpdateF()
      UpdateF();              // Update the Fock matrix

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
   PrintEHF();
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
   int norbits = modelspace->GetNumberOrbits();
   for (int i=0;i<norbits;i++)
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
   prev_energies = energies;
   for ( auto& it : Hbare.OneBodyChannels )
   {
      arma::uvec orbvec(it.second);
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



//*********************************************************************
/// Construct an unnormalized three-body monopole interaction
/// \f[ \langle iab | \bar{V}^{(3)} | jcd \rangle =
///     \sum\limits_{J,J_{12}}\sum_{Tt_{12}}(2J+1)(2T+1) 
///       \langle (ia)J_{12}t_{12};b JT| V^{(3)} | (jc)J_{12}t_{12}; d JT\rangle
/// \f]
//*********************************************************************
void HartreeFock::BuildMonopoleV3()
{
   double start_time = omp_get_wtime();
  // First, allocate. This is fast so don't parallelize.
  size_t norbits = modelspace->GetNumberOrbits();
  for (uint64_t i=0; i<norbits; ++i)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    int ei = 2*oi.n + oi.l;
    for (uint64_t j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      if (j<i) continue;
      Orbit& oj = modelspace->GetOrbit(j);
      int ej = 2*oj.n + oj.l;


      for (uint64_t a=0; a<norbits; ++a)
      {
        Orbit& oa = modelspace->GetOrbit(a);
        int ea = 2*oa.n + oa.l;
        for (uint64_t b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
        {
          Orbit& ob = modelspace->GetOrbit(b);
          int eb = 2*ob.n + ob.l;
     
            for (uint64_t c=0; c<norbits; ++c)
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

   Vmon3.resize( Vmon3_keys.size(), 0. );

   #pragma omp parallel for schedule(dynamic,1) 
   for (size_t ind=0; ind<Vmon3.size(); ++ind)
   {
      double v=0;
      int a,b,c,d,i,j;
      Vmon3UnHash( Vmon3_keys[ind], a,c,i,b,d,j);

      int j2a = modelspace->GetOrbit(a).j2;
      int j2c = modelspace->GetOrbit(c).j2;
      int j2i = modelspace->GetOrbit(i).j2;
      int j2b = modelspace->GetOrbit(b).j2;
      int j2d = modelspace->GetOrbit(d).j2;
      int j2j = modelspace->GetOrbit(j).j2;
 
      int j2min = std::max( std::abs(j2a-j2c), std::abs(j2b-j2d) )/2;
      int j2max = std::min( j2a+j2c, j2b+j2d )/2;
      for (int j2=j2min; j2<=j2max; ++j2)
      {
        int Jmin = std::max( std::abs(2*j2-j2i), std::abs(2*j2-j2j) );
        int Jmax = 2*j2 + std::min(j2i, j2j);
        for (int J2=Jmin; J2<=Jmax; J2+=2)
        {
           v += Hbare.ThreeBody.GetME_pn(j2,j2,J2,a,c,i,b,d,j) * (J2+1);
        }
      }
      v /= j2i+1.0;
      #pragma omp atomic write
      Vmon3[ind] = v ;
   }
   std::cout << "HartreeFock::BuildMonopoleV3  storing " << Vmon3.size() << " doubles for Vmon3 and "
             << Vmon3_keys.size() << " uint64's for Vmon3_keys." << std::endl;

   profiler.timer["HF_BuildMonopoleV3"] += omp_get_wtime() - start_time;
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
  a = (key    )&0x3FFL;
  b = (key>>10)&0x3FFL;
  c = (key>>20)&0x3FFL;
  d = (key>>30)&0x3FFL;
  e = (key>>40)&0x3FFL;
  f = (key>>50)&0x3FFL;
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



//********************************************************************
/// Get new occupation numbers by filling the orbits in order of their
/// single-particle energies. The last proton/neutron orbit can have
/// a fractional filling, corresponding to ensemble normal ordering.
//********************************************************************
void HartreeFock::FillLowestOrbits()
{
  // vector of indices such that they point to elements of F(i,i)
  // in ascending order of energy
  arma::mat F_trans = C.t() * F * C;
//  arma::uvec sorted_indices = arma::stable_sort_index( F.diag() );
  arma::uvec sorted_indices = arma::stable_sort_index( F_trans.diag() );
  int targetZ = modelspace->GetZref();
  int targetN = modelspace->GetAref() - targetZ;
  int placedZ = 0;
  int placedN = 0;
  std::vector<index_t> holeorbs_tmp; 
  std::vector<double> hole_occ_tmp;

  for (auto i : sorted_indices)
  {

    Orbit& oi = modelspace->GetOrbit(i);
    if (oi.tz2 < 0 and (placedZ<targetZ))
    {
      holeorbs_tmp.push_back(i);
      hole_occ_tmp.push_back( std::min(1.0,double(targetZ-placedZ)/(oi.j2+1) ) );
      placedZ = std::min(placedZ+oi.j2+1,targetZ);
    }
    else if (oi.tz2 > 0 and (placedN<targetN))
    {
      holeorbs_tmp.push_back(i);
      hole_occ_tmp.push_back( std::min(1.0,double(targetN-placedN)/(oi.j2+1) ) );
      placedN = std::min(placedN+oi.j2+1,targetN);
    }

    if((placedZ >= targetZ) and (placedN >= targetN) ) break;
  }

  holeorbs = arma::uvec( holeorbs_tmp );
  hole_occ = arma::rowvec( hole_occ_tmp );
//  hole_occ = arma::rowvec( hole_occ );
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
   int norbits = modelspace->GetNumberOrbits();
   Vij.zeros();
   V3ij.zeros();


   // This loop isn't thread safe for some reason. Regardless, parallelizing it makes it slower. 
   for (int i=0;i<norbits;i++)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      for (int j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         if (j<i) continue;
         for (int a=0;a<norbits;++a)
         {
            Orbit& oa = modelspace->GetOrbit(a);
            int Tz = (oi.tz2+oa.tz2)/2;
            int parity = (oi.l+oa.l)%2;
            int bra = modelspace->GetKetIndex(std::min(i,a),std::max(i,a));
            int local_bra = modelspace->MonopoleKets[Tz+1][parity][bra];
            for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
            {
               int ket = modelspace->GetKetIndex(std::min(j,b),std::max(j,b));
               int local_ket = modelspace->MonopoleKets[Tz+1][parity][ket];
               // 2body term <ai|V|bj>
               if ((a>i) xor (b>j))  // code needed some obfuscation, so threw an xor in there...
                  Vij(i,j) += rho(a,b)*Vmon_exch[Tz+1][parity](local_bra,local_ket); // <a|rho|b> * <ai|Vmon|jb>
               else
                  Vij(i,j) += rho(a,b)*Vmon[Tz+1][parity](local_bra,local_ket); // <a|rho|b> * <ai|Vmon|bj>
           }
         }
      }
      Vij.col(i) /= (oi.j2+1);
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
      V3vec.clear(); // free up the memory
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
     arma::uvec orbvec(it.second);
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
     for (index_t i=0;i<C.n_rows;++i) // loop through original basis states
     {
        if (C(i,i) < 0)  C.col(i) *= -1;
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
   OpHF.OneBody = C.t() * OpHO.OneBody * C;


   // Moderately difficult part:
   // Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz. Most of the effort here is in constructing D.

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
            if (ket_ho.p==ket_ho.q)    Dket(i,j) *= SQRT2;
            if (ket_hf.p==ket_hf.q)    Dket(i,j) /= SQRT2;
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
              if (bra_ho.p==bra_ho.q)    Dbra(i,j) *= SQRT2;
              if (bra_hf.p==bra_hf.q)    Dbra(i,j) /= SQRT2;
           }
        }
      }
      auto& IN  =  it.second;
      auto& OUT =  OpHF.TwoBody.GetMatrix(ch_bra,ch_ket);
      OUT  =    Dbra * IN * Dket;

   }

   return OpHF;
}


// If the lowest orbits are different from our previous guess, we should update the reference.
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


Operator HartreeFock::GetNormalOrderedH(arma::mat& Cin) 
{
  C=Cin;
//  ReorderCoefficients();  // Reorder columns of C so we can properly identify the hole orbits.
  UpdateDensityMatrix();  // Update the 1 body density matrix, used in UpdateF()
  UpdateF();              // Update the Fock matrix
  CalcEHF();
  PrintEHF();

  return GetNormalOrderedH();
}
/// Returns the normal-ordered Hamiltonian in the Hartree-Fock basis, neglecting the residual 3-body piece.
/// \f[ E_0 = E_{HF} \f]
/// \f[ f = C^{\dagger} F C \f]
/// \f[ \Gamma = D^{\dagger} \left(V^{(2)}+V^{(3\rightarrow 2)} \right) D \f]
/// \f[ V^{(2\rightarrow 3)J}_{ijkl} \equiv \frac{1}{\sqrt{(1+\delta_{ij})(1+\delta_{kl})}}\sum_{ab}\sum_{J_3}(2J_{3}+1)\rho_{ab}V^{JJJ_{3}}_{ijaklb} \f]
/// Where \f$ F\f$ is the Fock matrix obtained in UpdateF() and the matrix \f$ D\f$ is the same as the one defined in TransformToHFBasis().
///
Operator HartreeFock::GetNormalOrderedH() 
{
   double start_time = omp_get_wtime();
   std::cout << "Getting normal-ordered H in HF basis" << std::endl;

   // First, check if we need to update the occupation numbers for the reference
   if (not freeze_occupations)
   {
     bool changed_occupations = false;
     std::map<index_t,double> hole_map;
     for (index_t i=0;i<holeorbs.size();++i)  
     {
        hole_map[holeorbs[i]] = hole_occ[i];
        if ( std::abs(modelspace->GetOrbit( holeorbs[i] ).occ - hole_occ[i]) > 1e-3)
        {
           changed_occupations = true;
           std::cout << "After HF, occupation of orbit " << i << " has changed. Modelspace will be updated." << std::endl;
        }
     }
     if (changed_occupations)  modelspace->SetReference( hole_map );
   }

   Operator HNO = Operator(*modelspace,0,0,0,2);
   HNO.ZeroBody = EHF;
   HNO.OneBody = C.t() * F * C;

   int nchan = modelspace->GetNumberTwoBodyChannels();
   int norb = modelspace->GetNumberOrbits();
//   #pragma omp parallel for schedule(dynamic,1) // have not yet confirmed that this improves performance ... no sign of significant improvement
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
            if (bra.p==bra.q)    D(i,j) *= SQRT2;
            if (ket.p==ket.q)    D(i,j) /= SQRT2;

            // Now generate the NO2B part of the 3N interaction
            if (Hbare.GetParticleRank()<3) continue;
            if (i>j) continue;
            for (int a=0; a<norb; ++a)
            {
              Orbit & oa = modelspace->GetOrbit(a);
              if ( 2*oa.n+oa.l+e2bra > Hbare.GetE3max() ) continue;
              for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
              {
                Orbit & ob = modelspace->GetOrbit(b);
                if ( 2*ob.n+ob.l+e2ket > Hbare.GetE3max() ) continue;
                if ( std::abs(rho(a,b)) < 1e-8 ) continue; // Turns out this helps a bit (factor of 5 speed up in tests)
                int J3min = std::abs(2*J-oa.j2);
                int J3max = 2*J + oa.j2;
                for (int J3=J3min; J3<=J3max; J3+=2)
                {
                  V3NO(i,j) += rho(a,b) * (J3+1) * Hbare.ThreeBody.GetME_pn(J,J,J3,bra.p,bra.q,a,ket.p,ket.q,b);
                }
              }
            }
            V3NO(i,j) /= (2*J+1);
            if (bra.p==bra.q)  V3NO(i,j) /= SQRT2; 
            if (ket.p==ket.q)  V3NO(i,j) /= SQRT2; 
            V3NO(j,i) = V3NO(i,j);
         }
      }

     auto& V2  =  Hbare.TwoBody.GetMatrix(ch);
     auto& OUT =  HNO.TwoBody.GetMatrix(ch);
     OUT  =    D.t() * (V2 + V3NO) * D;
   }
   
//   FreeVmon();

   profiler.timer["HF_GetNormalOrderedH"] += omp_get_wtime() - start_time;
   
   return HNO;

}


void HartreeFock::FreeVmon()
{
   // free up some memory
   std::array< std::array< arma::mat,2>,3>().swap(Vmon);
   std::array< std::array< arma::mat,2>,3>().swap(Vmon_exch);
//   vector< pair<const array<int,6>,double>>().swap( Vmon3 );
//   vector< pair<const uint64_t,double>>().swap( Vmon3 );
   std::vector< double>().swap( Vmon3 );
   std::vector< uint64_t>().swap( Vmon3_keys );
}


/// Get the one-body generator corresponding to the transformation to the HF basis.
/// Since the unitary transformation for HF is given by the \f$ U_{HF} = C^{\dagger} \f$ matrix, we have
/// \f$ e^{-\Omega} = C \Rightarrow \Omega = -\log(C) \f$.
/// The log is evaluated by diagonalizing the one-body submatrix and taking the log of the diagonal entries.
/// This is much slower than the other methods, but it might be useful.
Operator HartreeFock::GetOmega()
{
   Operator Omega(*modelspace,0,0,0,1);
   Omega.SetAntiHermitian();
   for ( auto& it : Hbare.OneBodyChannels )
   {
      arma::uvec orbvec(it.second);
      arma::mat C_ch = C.submat(orbvec,orbvec);
      arma::cx_mat eigvect;
      arma::cx_vec eigval;
      arma::eig_gen(eigval, eigvect, C_ch);
      Omega.OneBody.submat(orbvec,orbvec) = -arma::real( eigvect * arma::diagmat(arma::log(eigval)) * eigvect.t()) ;
   }
   return Omega;
}


void HartreeFock::PrintSPE()
{
  arma::mat F_trans = C.t() * F * C;
  for (int i=0;i<modelspace->GetNumberOrbits();++i)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    std::cout << std::fixed << std::setw(3) << oi.n << " " << std::setw(3) << oi.l << " "
         << std::setw(3) << oi.j2 << " " << std::setw(3) << oi.tz2 << "   " << std::setw(10) << F_trans(i,i) << std::endl;
//         << std::setw(3) << oi.j2 << " " << std::setw(3) << oi.tz2 << "   " << std::setw(10) << F(i,i) << std::endl;
  }

}

void HartreeFock::PrintSPEandWF()
{
  arma::mat F_trans = C.t() * F * C;
  std::cout << std::fixed << std::setw(3) << "i" << ": " << std::setw(3) << "n" << " " << std::setw(3) << "l" << " "
       << std::setw(3) << "2j" << " " << std::setw(3) << "2tz" << "   " << std::setw(12) << "SPE" << " " << std::setw(12) << "occ." << "   |   " << " overlaps" << std::endl;
  for (int i=0;i<modelspace->GetNumberOrbits();++i)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    std::cout << std::fixed << std::setw(3) << i << ": " << std::setw(3) << oi.n << " " << std::setw(3) << oi.l << " "
         << std::setw(3) << oi.j2 << " " << std::setw(3) << oi.tz2 << "   " << std::setw(12) << std::setprecision(6) << F_trans(i,i) << " " << std::setw(12) << oi.occ << "   | ";
//         << std::setw(3) << oi.j2 << " " << std::setw(3) << oi.tz2 << "   " << std::setw(12) << std::setprecision(6) << F(i,i) << " " << std::setw(12) << oi.occ << "   | ";
    for (int j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      std::cout << std::setw(9) << C(i,j) << "  ";
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

double HartreeFock::GetRadialWF_r(index_t index, double R)
{
  double b = sqrt( (HBARC*HBARC) / (modelspace->GetHbarOmega() * M_NUCLEON) );
  Orbit& orb = modelspace->GetOrbit(index);
   double x = R/b;
   double psi = 0;
   for ( index_t j : Hbare.OneBodyChannels.at({orb.l,orb.j2,orb.tz2}) )
   {
     Orbit& oj = modelspace->GetOrbit(j);
     double Norm = 2*sqrt( gsl_sf_fact(oj.n) * pow(2,oj.n+oj.l) / M_SQRTPI / gsl_sf_doublefact(2*oj.n+2*oj.l+1) * pow(b,-3.0) );
     psi += C(index,j) * Norm * pow(x,oj.l) * exp(-x*x*0.5) * gsl_sf_laguerre_n(oj.n,oj.l+0.5,x*x);
   }
   return psi;
}



