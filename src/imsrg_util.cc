
#include "imsrg_util.hh"
#include "AngMom.hh"
#include "Commutator.hh"
#include "GaussLaguerre.hh"
//#include "DarkMatterNREFT.hh"
#include "omp.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h> // to use bessel functions
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
//#include <boost/math/special_functions/gamma.hpp>
//#include <boost/math/special_functions/factorials.hpp>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <array>

#define LOG2 log(2.0)
//using namespace AngMom;

/// imsrg_util namespace. Used to define some helpful functions.
namespace imsrg_util
{

 std::vector<std::string> split_string(std::string s, std::string delimiter)
 {
  std::vector<std::string> out;
   size_t last = 0;
   size_t next = 0;
   while ((next = s.find(delimiter, last)) != std::string::npos)
   {
     out.push_back( s.substr(last, next-last) );
     last = next + 1;
   }
   out.push_back( s.substr(last) );
   return out;
 }

 Operator OperatorFromString(ModelSpace& modelspace, std::string opname )
 {
      std::vector<std::string> opnamesplit = split_string( opname, "_" );  // split std::string on _ into a vector of std::string so that, e.g. "R2_p1"  =>  {"R2", "p1"}

           if (opname == "R2_p1")         return R2_1body_Op(modelspace,"proton") ;
      else if (opname == "R2_p2")         return R2_2body_Op(modelspace,"proton") ;
      else if (opname == "R2_n1")         return R2_1body_Op(modelspace,"neutron") ;
      else if (opname == "R2_n2")         return R2_2body_Op(modelspace,"neutron") ;
      else if (opname == "Rp2")           return Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "Rn2")           return Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "Rm2")           return Rm2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "E1")            return ElectricMultipoleOp(modelspace,1) ;
      else if (opname == "E2")            return ElectricMultipoleOp(modelspace,2) ;
      else if (opname == "E3")            return ElectricMultipoleOp(modelspace,3) ;
      else if (opname == "E4")            return ElectricMultipoleOp(modelspace,4) ;
      else if (opname == "E5")            return ElectricMultipoleOp(modelspace,5) ;
      else if (opname == "E6")            return ElectricMultipoleOp(modelspace,6) ;
      else if (opname == "E2int")         return IntrinsicElectricMultipoleOp(modelspace,2) ; // Untested
      else if (opname == "nE2")           return NeutronElectricMultipoleOp(modelspace,2) ;
      else if (opname == "M1")            return MagneticMultipoleOp(modelspace,1) ;
      else if (opname == "M2")            return MagneticMultipoleOp(modelspace,2) ;
      else if (opname == "M3")            return MagneticMultipoleOp(modelspace,3) ;
      else if (opname == "M4")            return MagneticMultipoleOp(modelspace,4) ;
      else if (opname == "M5")            return MagneticMultipoleOp(modelspace,5) ;
      else if (opname == "M1p")           return MagneticMultipoleOp_pn(modelspace,1,"proton") ;
      else if (opname == "M1n")           return MagneticMultipoleOp_pn(modelspace,1,"neutron") ;
      else if (opname == "Fermi")         return AllowedFermi_Op(modelspace) ;
      else if (opname == "GamowTeller")   return AllowedGamowTeller_Op(modelspace) ;
      else if (opname == "Iso2")          return Isospin2_Op(modelspace) ;
      else if (opname == "R2CM")          return R2CM_Op(modelspace) ;
      else if (opname == "Trel")          return Trel_Op(modelspace) ;
      else if (opname == "TCM")           return TCM_Op(modelspace) ;
      else if (opname == "Rso")           return RpSpinOrbitCorrection(modelspace) ;
      else if (opname == "RadialOverlap") return RadialOverlap(modelspace); // Untested...
      else if (opname == "Sigma")         return Sigma_Op(modelspace);
      else if (opname == "Sigma_p")       return Sigma_Op_pn(modelspace,"proton");
      else if (opname == "Sigma_n")       return Sigma_Op_pn(modelspace,"neutron");
      else if (opname == "L2rel")         return L2rel_Op(modelspace); // Untested...
      else if (opname == "QdotQ")         return QdotQ_Op(modelspace); // Untested...
      else if (opname == "VCoul")         return VCoulomb_Op(modelspace); // Untested...
      else if (opname == "hfsNMS")         return atomic_hfs::NormalMassShift(modelspace, 1);
      else if (opname == "hfsSMS")         return atomic_hfs::SpecificMassShift(modelspace, 1);
      else if (opname == "VCentralCoul")         return VCentralCoulomb_Op(modelspace); // Untested...
      else if (opnamesplit[0] =="HCM")
      {
         if ( opnamesplit.size() == 1 ) return HCM_Op(modelspace);
         double hw_HCM; // frequency of trapping potential
         std::istringstream( opnamesplit[1] ) >> hw_HCM;
         int A = modelspace.GetTargetMass();
         return TCM_Op(modelspace) + 0.5*A*M_NUCLEON*hw_HCM*hw_HCM/HBARC/HBARC*R2CM_Op(modelspace);
      }
      else if (opnamesplit[0] == "VCM") // GetHCM with a different frequency, ie HCM_24 for hw=24
      {
         double hw_VCM; // frequency of trapping potential
         std::istringstream(opnamesplit[1]) >> hw_VCM;
         int A = modelspace.GetTargetMass();
         return 0.5*A*M_NUCLEON*hw_VCM*hw_VCM/HBARC/HBARC*R2CM_Op(modelspace); 
      }
      else if (opnamesplit[0] == "Rp2Z") // Get point proton radius for specified Z, e.g. Rp2Z_10 for neon
      {
        int Z_rp;
        std::istringstream(opnamesplit[1]) >> Z_rp;
        return Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) ;
      }
      else if (opnamesplit[0] == "Rp2AZ") // Get point proton radius for specified A and Z, e.g. Rp2AZ_20_10 for neon
      {
        int A_rp;
        int Z_rp;
        std::istringstream(opnamesplit[1]) >> A_rp;
        std::istringstream(opnamesplit[2]) >> Z_rp;
        return Rp2_corrected_Op(modelspace,A_rp,Z_rp) ;
      }
      else if (opnamesplit[0] == "Rn2Z") // Get point neutron radius for specified Z
      {
        int Z_rp;
        std::istringstream(opnamesplit[1]) >> Z_rp;
        return Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) ;
      }
      else if (opnamesplit[0] == "rhop") // point proton  density at position r, e.g. rhop_1.25
      {
        double rr;
        std::istringstream(opnamesplit[1]) >> rr;
        return ProtonDensityAtR(modelspace,rr);
      }
      else if (opnamesplit[0] == "rhon") // point neutron density at position r
      {
        double rr;
        std::istringstream(opnamesplit[1]) >> rr;
        NeutronDensityAtR(modelspace,rr);
      }
      else if (opnamesplit[0] == "OneOcc") // Get occupation of specified orbit, e.g. OneOccp_1p3
      {
         index_t ind = modelspace.String2Index( {  opnamesplit[1] } )[0];
         Orbit& oi = modelspace.GetOrbit(ind);
         return NumberOp(modelspace,oi.n,oi.l,oi.j2,oi.tz2) ;
      }
      else if (opnamesplit[0]== "AllOcc") // Get occupation of orbit, summed over all values of radial quantum number n, e.g. AllOccpp3
      {
         index_t ind = modelspace.String2Index( { "0"+ opnamesplit[1] } )[0];
         Orbit& oi = modelspace.GetOrbit(ind);
         return NumberOpAlln(modelspace,oi.l,oi.j2,oi.tz2) ;
      }
      else if (opnamesplit[0] == "protonFBC") // Fourier bessel coefficient of order nu
      {
         int nu;
         std::istringstream(opnamesplit[1]) >> nu;
         return FourierBesselCoeff( modelspace, nu, 8.0, modelspace.proton_orbits);
      }
      else if (opnamesplit[0] == "neutronFBC") // Fourier bessel coefficient of order nu
      {
         int nu;
         std::istringstream(opnamesplit[1]) >> nu;
         return FourierBesselCoeff( modelspace, nu, 8.0, modelspace.neutron_orbits) ;
      }
      else if (opnamesplit[0] == "M0nu" and opnamesplit[1] == "TBME") // 0\nu\beta\beta decay TBME, M0nu_TBME_${Nq}_${SRC} (CP)
      {
         int Nquad; // number of quadrature points
         std::string src; // chosen SRC parameters (none, Argonne, CD-Bonn, Miller/Spencer)
         std::istringstream(opnamesplit[2]) >> Nquad;
         std::istringstream(opnamesplit[3]) >> src;
         return M0nu_TBME_Op(modelspace,Nquad,src);
      }
//      else if (opnamesplit[0] == "DMNREFT") // point radius density at position r, e.g. rhop1.25
//      {
//        double q;
//        int J;
//        std::string dmopname = opnamesplit[1];
//        std::istringstream(opnamesplit[2]) >> q;
//        std::istringstream(opnamesplit[3]) >> J;
//
//        std::map<string, Operator (*)(ModelSpace&, int, double) > dmop = { {"M",       &DM_NREFT::M},
//                                                                           {"Sigma",   &DM_NREFT::Sigma},
//                                                                           {"Sigmap",  &DM_NREFT::Sigmap},
//                                                                           {"Sigmapp", &DM_NREFT::Sigmapp},
//                                                                           {"Delta",   &DM_NREFT::Delta},
//                                                                           {"Deltap",  &DM_NREFT::Deltap},
//                                                                           {"Phipp",   &DM_NREFT::Phipp},
//                                                                           {"Phitp",   &DM_NREFT::Phitp},
//                                                                           {"Omega",   &DM_NREFT::Omega},
//                                                                         };
//        if ( dmop.find(dmopname) != dmop.end() )
//        {
//        return dmop[dmopname](modelspace, J, q );
//        }
//      }
      else if (opnamesplit[0] == "Dagger" or opnamesplit[0] == "DaggerHF" )
      {
        index_t Q = modelspace.String2Index({opnamesplit[1]})[0];
        return Dagger_Op( modelspace, Q);
      }
      else //need to remove from the list
      {
         std::cout << "Unknown operator: " << opname << std::endl;
      }
      return Operator();
 
 }

 Operator NumberOp(ModelSpace& modelspace, int n, int l, int j2, int tz2)
 {
   Operator NumOp = Operator(modelspace);
   int indx = modelspace.Index1(n,l,j2,tz2);
   NumOp.ZeroBody = 0;
   NumOp.EraseOneBody();
   NumOp.EraseTwoBody();
   NumOp.OneBody(indx,indx) = 1;
   return NumOp;
 }

 Operator NumberOpAlln(ModelSpace& modelspace, int l, int j2, int tz2)
 {
   Operator NumOp = Operator(modelspace);
   for (int n=0;n<=(modelspace.GetEmax()-l)/2;++n)
   {
     int indx = modelspace.Index1(n,l,j2,tz2);
     NumOp.OneBody(indx,indx) = 1;
   }
   return NumOp;
 }

 double HO_density(int n, int l, double hw, double r)
 {
    double v = M_NUCLEON * hw / (HBARC*HBARC);
    double Norm = pow(v/2.,1.5+l) * M_SQRT2/M_SQRTPI * pow(2,n+2*l+3) * gsl_sf_fact(n) / gsl_sf_doublefact(2*n + 2*l + 1);
    double L = gsl_sf_laguerre_n(n, l+0.5, v*r*r);
    double rho = Norm * pow(r,2*l) * exp(-v * r*r) * L * L;
    return rho;
 }

double HO_Radial_psi(int n, int l, double hw, double r)
{
   double b = sqrt( (HBARC*HBARC) / (hw * M_NUCLEON) );
   double x = r/b;
   double Norm = 2*sqrt( gsl_sf_fact(n) * pow(2,n+l) / M_SQRTPI / gsl_sf_doublefact(2*n+2*l+1) / pow(b,3.0) );
   double L = gsl_sf_laguerre_n(n,l+0.5,x*x);
   double psi = Norm * pow(x,l) * exp(-x*x*0.5) * L;
   return psi;
}

 // Just do the HF transformation
 std::vector<double> GetOccupationsHF(HartreeFock& hf)
 {
    ModelSpace* modelspace = hf.Hbare.modelspace;
    int norb = modelspace->GetNumberOrbits();
    std::vector<double> occupation(norb);

    for (int i=0; i<norb; ++i)
    {
      Orbit & oi = modelspace->GetOrbit(i);
      // Get the number operator for orbit i
      Operator N_bare = NumberOp(*modelspace,oi.n,oi.l,oi.j2,oi.tz2);
      // Transform it to the normal-ordered HF basis
      Operator N_NO = hf.TransformToHFBasis(N_bare).DoNormalOrdering();
      occupation[i] = N_NO.ZeroBody;
      std::cout << oi.n << " " << oi.l << " " << oi.j2 << "/2 " << occupation[i] << std::endl;
    }
    return occupation;
 }

 // Do the full IMSRG transformation
 std::vector<double> GetOccupations(HartreeFock& hf, IMSRGSolver& imsrgsolver)
 {
    ModelSpace* modelspace = imsrgsolver.modelspace;
    int norb = modelspace->GetNumberOrbits();
    std::vector<double> occupation(norb,0);

    for (int i=0; i<norb; ++i)
    {
      Orbit & oi = modelspace->GetOrbit(i);
      // Get the number operator for orbit i
      Operator N_bare = NumberOp(*modelspace,oi.n,oi.l,oi.j2,oi.tz2);
      // Transform it to the normal-ordered HF basis
      Operator N_NO = hf.TransformToHFBasis(N_bare).DoNormalOrdering();
      // Transform to the imsrg evolved basis
      Operator N_final = imsrgsolver.Transform(N_NO);
      occupation[i] = N_final.ZeroBody;
    }
    return occupation;
 }

 std::vector<double> GetDensity( std::vector<double>& occupation, std::vector<double>& R, std::vector<int>& orbits, ModelSpace& modelspace )
 {
     int nr_steps = R.size();
     double hw = modelspace.GetHbarOmega();
     std::vector<double> dens(nr_steps,0);
     for (int& i : orbits)
     {
       Orbit & oi = modelspace.GetOrbit(i);
       for (int ir=0;ir<nr_steps;++ir)
       {
          dens[ir] += HO_density(oi.n,oi.l,hw,R[ir]) * occupation[i];
       }
     }
     return dens;
 }


 void WriteSPWaveFunctions( std::vector<std::string>& spwf, HartreeFock& hf, std::string intfile )
 {
   int n_radial_points = 200;
   double Rmax = 20.0;
   std::vector<index_t> spwf_indices = hf.modelspace->String2Index(spwf);
   std::vector<double> R(n_radial_points);
   std::vector<double> PSI(n_radial_points);
   for (int rstep=0;rstep<n_radial_points;++rstep)   R[rstep] = Rmax/n_radial_points * rstep;
   for ( index_t i=0; i< spwf.size(); ++i)
   {
     hf.GetRadialWF(spwf_indices[i], R, PSI);
     std::ofstream wf_file (intfile + "_spwf_" + spwf[i] + ".dat");
     for ( index_t rstep=0; rstep<R.size(); ++rstep)  wf_file << std::fixed << std::setw(10) << std::setprecision(7) << R[rstep] << "   " << std::setw(10) << std::setprecision(7) << PSI[rstep] << std::endl;
     wf_file.close();
   }
 }


 Operator Single_Ref_1B_Density_Matrix(ModelSpace& modelspace)
 {
    Operator DM(modelspace,0,0,0,2);
//    for (index_t a : modelspace.holes)
    for (auto& a : modelspace.holes)
    {
//       index_t a = it_a.first;
        Orbit& oa = modelspace.GetOrbit(a);
//       DM.OneBody(a,a) = 1.0;
       DM.OneBody(a,a) = oa.j2+1.0;
    }
    return DM;
 }


 double Get_Charge_Density( Operator& DM, double r)
 {
   ModelSpace* modelspace = DM.GetModelSpace();
   double hw = modelspace->GetHbarOmega();
   double rho=0;
   for (index_t i : modelspace->proton_orbits )
   {
      Orbit& oi = modelspace->GetOrbit(i);
      rho += (oi.j2+1) * DM.OneBody(i,i) * HO_density(oi.n, oi.l, hw, r);
   }
   return rho;
 }

/*
 double Get_Charge_Density( Operator& DM, double r)
 {
   ModelSpace* modelspace = DM.GetModelSpace();
   double hw = modelspace->GetHbarOmega();
   double rho=0;
   for (index_t i : modelspace->proton_orbits )
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (std::abs(DM.OneBody(i,i))<1e-7) continue;
//      std::cout << i << " " << (oi.j2+1)*DM.OneBody(i,i) << std::endl;
//      rho += (oi.j2+1) * DM.OneBody(i,i) * HO_density(oi.n, oi.l, hw, r);
      rho += (1) * DM.OneBody(i,i) * HO_density(oi.n, oi.l, hw, r);
   }
   return rho;
 }
*/


// double Get_Charge_Density( Operator& DM, double r)
// {
//   ModelSpace* modelspace = DM.GetModelSpace();
//   double hw = modelspace->GetHbarOmega();
//   double rho=0;
//   for (index_t i : modelspace->proton_orbits )
//   {
//      Orbit& oi = modelspace->GetOrbit(i);
//      for ( index_t j : DM.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
//      {
//        if (std::abs(DM.OneBody(i,j))<1e-7) continue;
//        Orbit& oj = modelspace->GetOrbit(j);
//        rho += DM.OneBody(i,j) * HO_Radial_psi(oi.n, oi.l, hw, r) * HO_Radial_psi(oj.n, oj.l, hw, r);
//      }
//   }
//   return rho;
// }

Operator KineticEnergy_Op(ModelSpace& modelspace)
{
   Operator T(modelspace,0,0,0,2);
//   int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
//   for (int a=0;a<norbits;++a)
   for ( auto a : modelspace.all_orbits )
   {
      Orbit & oa = modelspace.GetOrbit(a);
      T.OneBody(a,a) = 0.5 * hw * (2*oa.n + oa.l +3./2); 
      for ( int b : T.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
         if (b<=a) continue;
         Orbit & ob = modelspace.GetOrbit(b);
         if (oa.n == ob.n+1)
            T.OneBody(a,b) = 0.5 * hw * sqrt( (oa.n)*(oa.n + oa.l +1./2));
         else if (oa.n == ob.n-1)
            T.OneBody(a,b) = 0.5 * hw * sqrt( (ob.n)*(ob.n + ob.l +1./2));
         T.OneBody(b,a) = T.OneBody(a,b);
      }
   }
//   std::cout << "Done with kinetic energy op. One body part looks like" << std::endl << T.OneBody << std::endl << std::endl;
   return T;
}

/// Relative kinetic energy, includingthe hw/A factor
/// \f[
/// T_{rel} = T - T_{CM}
/// \f]
 Operator Trel_Op(ModelSpace& modelspace)
 {
   Operator Trel( KineticEnergy_Op(modelspace) );
   Trel -= TCM_Op(modelspace);
   return Trel;
 }


// Leading relativistic correction to the kinetic energy
// T = p^2/2m  - p^4/(8m^3c^2) + ...
// We use p^4/(4m^2) = (p^2/2m)^2
// and the fact that the nonrelativistic kinetic energy acts as a 'hopping operator'
// so acting p^2/2m on a state |nl> returns a linear combination of |nl> |(n+1)l> and |(n-1)l>.
// Specifically, we write this as (p^2/2m) |nl> = N0 |nl> + Np |(n+1)l> + Nm |(n-1)l>
// where N0 = (2n+l+3/2)hw/2,   Np = sqrt((n+1)(n+l+3/2))hw/2,   Nm = sqrt(n(n+l+1/2))hw/2.
// Then we act again with p^2/2m to get the action of p^4/4m.
// p^4/(4m^2) |n> = Npp|n+2> + (Np0+N0p)|n+1> + (N00 + Npm + Nmp)|n> + (Nm0+N0m)|n-1> + Nmm |n-2>
// where
// Np0 = sqrt[ (n+1)(n+l+3/2) ] (2n+l+3/2)             (hw/2)^2
// N0p = (2n+l+7/2) sqrt[ (n+1)(n+l+3/2) ]             (hw/2)^2
// N00 = (2n+2+3/2) (2n+l+3/2)                         (hw/2)^2
// Nm0 = sqrt[ n(n+l+1/2) ]  (2n+l+3/2)                (hw/2)^2
// N0m = (2n+l-1/2) sqrt[ n(n+l+1/2) ]                 (hw/2)^2
// Npp = sqrt[ (n+2)(n+l+5/2)] sqrt[ (n+1)(n+l+3/2) ]  (hw/2)^2
// Npm = sqrt[ n(n+l+1/2)]      sqrt[ n(n+l+1/2) ]     (hw/2)^2
// Nmp = sqrt[ (n+1)(n+l+3/2)]  sqrt[ (n+1)(n+l+3/2) ] (hw/2)^2
// Nmm = sqrt[ (n-1)(n+l-1/2) ] sqrt[ n(n+l+1/2) ]     (hw/2)^2
Operator KineticEnergy_RelativisticCorr(ModelSpace& modelspace)
{
   Operator Trc(modelspace);
   int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
   double coeff = 0.25*hw*hw;
   for (int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      int na = oa.n;
      int la = oa.l;
      // I have pulled off a global factor (1/2 hw)^2 above and called it coeff.
      double Npp = sqrt( (na+2)*(na+la+2.5) ) * sqrt( (na+1)*(na+la+1.5) ) ;
      double Np0 = sqrt( (na+1)*(na+la+1.5) ) * (2*na+la+1.2) ;
      double N0p = (2*na+la+3.5) * sqrt((na+1)*(na+la+1.5)) ;
      double N00 = (2*na+la+1.5) * (2*na*la+1.5) ;
      double Nmp = sqrt( (na+1)*(na+la+1.5) ) * sqrt( (na+1)*(na+la+1.5) ) ;
      double Npm = na<1 ? 0 : sqrt( na*(na+la+0.5) ) * sqrt( na*(na+la+0.5) ) ;
//      double Nm0 = sqrt( na*(na+la+0.5)) * (2*na+la+1.5) ;
//      double N0m = na<1 ? 0 : (2*na+la-0.5) * sqrt( na*(na+la+0.5)) ;
//      double Nmm = na<2 ? 0 : sqrt( (na-1)*(na+la-0.5) ) * sqrt( na*(na+la+0.5) );

      Trc.OneBody(a,a) =  N00 + Npm + Nmp ; 
      // off-diagonal terms
      for ( int b : Trc.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
         Orbit & ob = modelspace.GetOrbit(b);
         if ( ob.n<=oa.n or ob.n>oa.n+2 ) continue;
         if (ob.n == oa.n+1)
            Trc.OneBody(b,a) = Np0 + N0p;
         else if (ob.n == oa.n+2)
            Trc.OneBody(b,a) = Npp;
         // Make it hermitian:
         Trc.OneBody(a,b) = Trc.OneBody(b,a);
      }
   }  // 
   return -Trc*coeff /(2*M_NUCLEON);  // the minus sign is put in here
}



/// Center of mass kinetic energy, including the hw/A factor
/// \f[
/// T = \frac{\hbar\omega}{A}\sum_{ij} t_{ij} a^{\dagger}_{i} a_{j} + \frac{\hbar\omega}{A}\frac{1}{4}\sum_{ijkl} t_{ijkl} a^{\dagger}_{i}a^{\dagger}_{j}a_{l}a_{k}
/// \f]
/// with a one-body piece
/// \f[
/// t_{ij} = \frac{1}{\hbar\omega} \left\langle i | T_{12} | j \right\rangle = \frac{1}{2}(2n_i+\ell_i+3/2) \delta_{ij} + \frac{1}{2}\sqrt{n_j(n_j+\ell_j+\frac{1}{2})} \delta_{n_i,n_j-1}\delta_{k_i k_j}
/// \f]
/// where \f$k\f$ labels all quantum numbers other than \f$n\f$ and a two-body piece
/// \f[  
/// t_{ijkl} = \frac{1}{\hbar\omega} \left\langle ij | (T^{CM}_{12} - T^{rel}_{12}) | kl \right\rangle
/// \f]
 Operator TCM_Op(ModelSpace& modelspace)
 {
   double t_start = omp_get_wtime();
   int E2max = modelspace.GetE2max();
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   Operator TcmOp = Operator(modelspace);
   TcmOp.SetHermitian();
   // One body piece = p**2/(2mA)
   int norb = modelspace.GetNumberOrbits();
   for (int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      for (int j : TcmOp.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         Orbit & oj = modelspace.GetOrbit(j);
         if (oj.n<oi.n) continue;
         double tij = 0;
         if (oi.n == oj.n) tij = 0.5*(2*oi.n+oi.l + 1.5) * hw/A;
         else if (oi.n == oj.n-1) tij = 0.5*sqrt(oj.n*(oj.n+oj.l + 0.5)) * hw/A;
         TcmOp.OneBody(i,j) = tij;
         TcmOp.OneBody(j,i) = tij;
      }
   }

   // Two body piece = 2*p1*p2/(2mA) = (Tcm-Trel)/A
   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1)  // In order to make this parallel, need to precompute the Moshinsky brackets.
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         Orbit & oi = modelspace.GetOrbit(bra.p);
         Orbit & oj = modelspace.GetOrbit(bra.q);
         if ( 2*(oi.n+oj.n)+oi.l+oj.l > E2max) continue;
         for (int iket=ibra;iket<nkets;++iket)
         {
            
            Ket & ket = tbc.GetKet(iket);
            Orbit & ok = modelspace.GetOrbit(ket.p);
            Orbit & ol = modelspace.GetOrbit(ket.q);
            if ( 2*(ok.n+ol.n)+ok.l+ol.l > E2max) continue;
//            double p1p2 = Calculate_p1p2(modelspace,bra,ket,tbc.J) * hw/A;
            double p1p2 = Calculate_p1p2(modelspace,bra,ket,tbc.J) / A;
            if (std::abs(p1p2)>1e-7)
            {
              TcmOp.TwoBody.SetTBME(ch,ibra,iket,p1p2);
            }
         }
      }
   }
   TcmOp.profiler.timer["TCM_Op"] += omp_get_wtime() - t_start;
   return TcmOp;
 }


 // evaluate <bra| p1*p2/2 | ket> 
/// This returns the antisymmetrized J-coupled two body matrix element of \f$ \vec{p}_1 \cdot \vec{p}_2 / (m) \f$.
/// The formula is
/// \f{eqnarray*}{
/// \frac{1}{m}
/// \left \langle a b \right| \vec{p}_1 \cdot \vec{p}_2 \left| c d \right \rangle_J
/// = \frac{1}{\sqrt{(1+\delta_{ab})(1+\delta_{cd})}} &\sum\limits_{LS}
/// \left[ \begin{array}{ccc}
///  \ell_a & s_a & j_a \\
///  \ell_b & s_b & j_b \\
///  L      & S   & J
/// \end{array} \right]
/// \left[ \begin{array}{ccc}
///  \ell_c & s_c & j_c \\
///  \ell_d & s_d & j_d \\
///  L      & S   & J
/// \end{array} \right] & \\
/// & \times \sum\limits_{\substack{N_{ab}N_{cd} \Lambda \\ n_{ab} n_{cd} \lambda}} \mathcal{A}_{abcd}^{\lambda S} \times
/// \left\langle N_{ab}\Lambda n_{ab} \lambda | n_{a} \ell_{a} n_{b} \ell_{b} \right\rangle_{L} &
/// \left\langle N_{cd}\Lambda n_{cd} \lambda | n_{c} \ell_{c} n_{d} \ell_{d} \right\rangle_{L} \\ & \times
/// \left( \left\langle N_{ab}\Lambda | t_{cm} | N_{cd} \Lambda \right \rangle
/// -\left\langle n_{ab}\lambda | t_{rel} | n_{cd} \lambda \right \rangle \right)
/// \f}
/// The antisymmetrization factor \f$ \mathcal{A}_{abcd}^{\lambda S} \f$ ensures that
/// the relative wave function is antisymmetrized.
/// It is given by \f$ \mathcal{A}_{abcd}^{\lambda S} = \left|t_{za}+t_{zc}\right| + \left| t_{za} + t_{zd} \right| (-1)^{\lambda + S + \left|T_z\right|}\f$ .
///
/// The center-of-mass and relative kinetic energies can be found by the same equation as used in the one-body piece of TCM_Op()
 double Calculate_p1p2(ModelSpace& modelspace, Ket & bra, Ket & ket, int J)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
   int nb = ob.n;
   int nc = oc.n;
   int nd = od.n;

   int la = oa.l;
   int lb = ob.l;
   int lc = oc.l;
   int ld = od.l;

   double ja = oa.j2/2.0;
   double jb = ob.j2/2.0;
   double jc = oc.j2/2.0;
   double jd = od.j2/2.0;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   // p1*p2 only connects kets with delta N = 0,1 ==> delta E = 0,2
   if (std::abs(fab-fcd)>2 or std::abs(fab-fcd)%2 >0 ) return 0; 

   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double p1p2=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( std::abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<std::abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
              // factor to account for antisymmetrization

              int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (std::abs(mosh_ab)<1e-8) continue;

              for (int N_cd=std::max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd and N_ab != N_cd) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (std::abs(mosh_cd)<1e-8) continue;

                double tcm = 0;
                double trel = 0;
                if (n_ab == n_cd)
                {
                  if      (N_ab == N_cd)   tcm = (2*N_ab+Lam_ab+1.5);
                  else if (N_ab == N_cd+1) tcm = sqrt(N_ab*( N_ab+Lam_ab+0.5));
                  else if (N_ab == N_cd-1) tcm = sqrt(N_cd*( N_cd+Lam_ab+0.5));
                }
                if (N_ab == N_cd)
                {
                  if      (n_ab == n_cd)   trel = (2*n_ab+lam_ab+1.5);
                  else if (n_ab == n_cd+1) trel = sqrt(n_ab*( n_ab+lam_ab+0.5));
                  else if (n_ab == n_cd-1) trel = sqrt(n_cd*( n_cd+lam_cd+0.5));
                }
                double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                p1p2 += (tcm-trel) * prefactor ;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize. The 0.5 comes from t ~ 0.5 * (N+3/2) hw
   p1p2 *= 0.5*modelspace.GetHbarOmega() / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return p1p2 ;

 }


/// Correction to Trel due to the proton-neutron mass differences
/// \f[
/// T_{rel} = T - T_{CM}
/// \f]
///
/// \f[
/// \delta T = \sum_i \frac{p_i^2}{2m} \left( \frac{m-m_i}{m_i} \right)
///\f]
///
/// \f[
/// \delta T_{CM} = \left( \frac{Am}{Zm_p+Nm_n} \right) \frac{1}{2mA} P_{CM}^2
///\f]
 Operator Trel_Masscorrection_Op(ModelSpace& modelspace)
 {
   Operator dTrel( modelspace );
   int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
   double m_avg = 0.5*(M_PROTON+M_NEUTRON);
   int A = modelspace.GetTargetMass();
   int Z = modelspace.GetTargetZ();
   int N = A-Z;

   for (int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      double m_a = (oa.tz2 == -1) ? M_PROTON : M_NEUTRON ;
      double correction = (m_avg-m_a)/m_a;
      dTrel.OneBody(a,a) = correction * 0.5 * hw * (2*oa.n + oa.l +3./2); 
      for ( int b : dTrel.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
         if (b<=a) continue;
         Orbit & ob = modelspace.GetOrbit(b);
         if (oa.n == ob.n+1)
            dTrel.OneBody(a,b) = correction * 0.5 * hw * sqrt( (oa.n)*(oa.n + oa.l +1./2));
         else if (oa.n == ob.n-1)
            dTrel.OneBody(a,b) = correction * 0.5 * hw * sqrt( (ob.n)*(ob.n + ob.l +1./2));
         dTrel.OneBody(b,a) = dTrel.OneBody(a,b);
      }
   }

   double CM_correction = A*m_avg / (Z*M_PROTON + N*M_NEUTRON) - 1;

   dTrel -= CM_correction * TCM_Op(modelspace);
   return dTrel;
 }






////////////////////////////////////////////////////////////////////////////
/////////////  IN PROGRESS...  doesn't work yet, and for now it's slower.  /
////////////////////////////////////////////////////////////////////////////

 void Calculate_p1p2_all(Operator& OpIn)
 {
   ModelSpace* modelspace = OpIn.GetModelSpace();
//   modelspace->PreCalculateMoshinsky();
   for ( int ch : modelspace->SortedTwoBodyChannels )
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int parity = tbc.parity;
//      int Tz = tbc.Tz;
      arma::mat& MatJJ = OpIn.TwoBody.GetMatrix(ch);
      int nkets_JJ = tbc.GetNumberKets();

      // Find the maximum oscillator energy for this channel
      Ket& ketlast = tbc.GetKet( tbc.GetNumberKets()-1 );
      int emax_ket = 2*ketlast.op->n + 2*ketlast.oq->n + ketlast.op->l + ketlast.oq->l;

      std::vector<std::array<int,6>> JacobiBasis;  // L,S,N,Lambda,n,lambda
      for (int L=std::max(J-1,0); L<=J+1; ++L)
      {
       for ( int S=std::abs(J-L); S<=1; ++S)
       {
        for ( int N=0; N<=emax_ket/2; ++N )
        {
         for ( int Lambda=0; Lambda<=(emax_ket-2*N); ++Lambda)
         {
          for ( int lambda=std::abs(L-Lambda)+(L+parity)%2; lambda<=std::min(Lambda+L,emax_ket-2*N-Lambda); lambda+=2)
          {
           for ( int n =0; n<=(emax_ket-2*N-Lambda-lambda)/2; ++n)
           {
             JacobiBasis.push_back({L,S,N,Lambda,n,lambda});
           }
          }
         }
        }
       }
      }


      int nkets_Jacobi = JacobiBasis.size();
      arma::mat MatJacobi(nkets_Jacobi,nkets_Jacobi);
      arma::mat Trans(nkets_Jacobi,nkets_JJ);
      int n_nonzero = 0;
      for (int iJJ=0; iJJ<nkets_JJ; ++iJJ)
      {
        Ket & ket = tbc.GetKet(iJJ);
        int la = ket.op->l;
        int lb = ket.oq->l;
        int na = ket.op->n;
        int nb = ket.oq->n;
        double ja = ket.op->j2*0.5;
        double jb = ket.oq->j2*0.5;
//        int ta = ket.op->n;
//        int tb = ket.oq->n;
        for (int iJac=0; iJac<nkets_Jacobi; ++iJac)
        {
          int L      = JacobiBasis[iJac][0];
          int S      = JacobiBasis[iJac][1];
          int N      = JacobiBasis[iJac][2];
          int Lambda = JacobiBasis[iJac][3];
          int n      = JacobiBasis[iJac][4];
          int lambda = JacobiBasis[iJac][5];

//          int Asym = 1; // Fix this...
          double ninej = AngMom::NormNineJ(la,0.5,ja,lb,0.5,jb,L,S,J);
          if (std::abs(ninej)<1e-6) continue;
          double mosh = modelspace->GetMoshinsky(N,Lambda,n,lambda,na,la,nb,lb,L);
          if (std::abs(mosh)<1e-6) continue;
          Trans(iJac,iJJ) = ninej * mosh;
          n_nonzero += 1;
           
        }
      }
      for (int i=0; i<nkets_Jacobi; ++i)
      {
        int Li      = JacobiBasis[i][0];
        int Si      = JacobiBasis[i][1];
        int Ni      = JacobiBasis[i][2];
        int Lambdai = JacobiBasis[i][3];
        int ni      = JacobiBasis[i][4];
        int lambdai = JacobiBasis[i][5];
        for (int j=0; j<nkets_Jacobi; ++j)
        {
          int Lj      = JacobiBasis[j][0];
          int Sj      = JacobiBasis[j][1];
          int Nj      = JacobiBasis[j][2];
          int Lambdaj = JacobiBasis[j][3];
          int nj      = JacobiBasis[j][4];
          int lambdaj = JacobiBasis[j][5];
          if(Li!=Lj or Si!=Sj or Lambdai!=Lambdaj or lambdai!=lambdaj) continue;
          double tcm = 0;
          double trel = 0;
          if (ni == nj)
          {
            if      (Ni == Nj)   tcm = (2*Ni+Lambdai+1.5);
            else if (Ni == Nj+1) tcm = sqrt(Ni*( Ni+Lambdai+0.5));
            else if (Ni == Nj-1) tcm = sqrt(Nj*( Nj+Lambdai+0.5));
          }
          if (Ni == Nj)
          {
            if      (ni == nj)   tcm = (2*ni+lambdai+1.5);
            else if (ni == nj+1) tcm = sqrt(ni*( ni+lambdai+0.5));
            else if (ni == nj-1) tcm = sqrt(nj*( nj+lambdai+0.5));
          }
          MatJacobi(i,j) = tcm - trel;
        }
      }


      MatJJ = Trans.t() * MatJacobi * Trans;
      std::cout << "ch = " << ch << "   size of JJ basis = " << nkets_JJ << "  size of Jacobi Basis = " << nkets_Jacobi << "   nonzero matrix elements = " << n_nonzero << std::endl;

   }
 }

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////





// Center of mass R^2, in units of fm^2
/// Returns
/// \f[ 
/// R^{2}_{CM} = \left( \frac{1}{A}\sum_{i}\vec{r}_{i}\right)^2 =
/// \frac{1}{A^2} \left( \sum_{i}r_{i}^{2} + 2\sum_{i<j}\vec{r}_i\cdot\vec{r}_j  \right)
/// \f]
/// evaluated in the oscillator basis.


 Operator R2CM_Op(ModelSpace& modelspace)
 {
//   Operator R2cmOp = Operator(modelspace);
   Operator R2cmOp = RSquaredOp(modelspace);

//   unsigned int norb = modelspace.GetNumberOrbits();
//   double oscillator_b2 = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
//   for (unsigned int i=0; i<norb; ++i)
//   {
//      Orbit & oi = modelspace.GetOrbit(i);
//      for (auto j : R2cmOp.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
//      {
//         if (j<i) continue;
//         Orbit & oj = modelspace.GetOrbit(j);
//         double rij = 0;
//         if (oi.n == oj.n)        rij = (2*oi.n+oi.l + 1.5);
//         else if (oi.n == oj.n-1) rij = -sqrt(oj.n*(oj.n+oj.l + 0.5));
////         R2cmOp.OneBody(i,j) = rij;
////         R2cmOp.OneBody(j,i) = rij;
//         R2cmOp.OneBody(i,j) = rij * oscillator_b2;
//         R2cmOp.OneBody(j,i) = rij * oscillator_b2;
//      }
//   }

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            // factor of 2 comes from limiting sum to i<j. Otherwise it would be r1*r2 + r2*r1.
            double mat_el = 2*Calculate_r1r2(modelspace,bra,ket,tbc.J); 
             
            R2cmOp.TwoBody.SetTBME(ch,ibra,iket,mat_el);
            R2cmOp.TwoBody.SetTBME(ch,iket,ibra,mat_el);
         }
      }
   }
//   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
//   return R2cmOp * (HBARC*HBARC/M_NUCLEON/hw)/(A*A);
   return R2cmOp /(A*A);
 }



// Intrinsic point proton radius squared
/// Returns
/// \f[ 
/// R_p^{2} = \frac{1}{Z} \sum_{p}\left(\vec{r}_{p}-\vec{R}_{CM}\right)^2 =
/// R^2_{CM} + \frac{A-2}{AZ} \sum_{p}r_{p}^{2} - \frac{4}{AZ}\sum_{i<j}\vec{r}_i\cdot\vec{r}_j  
/// \f]
/// evaluated in the oscillator basis.
 Operator Rp2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   if (Z==0) return 0.0*KineticEnergy_Op(modelspace);
   return R2CM_Op(modelspace) + (A-2.0)/(A*Z)*R2_1body_Op(modelspace,"proton")
                                   - 4./(A*Z)*R2_2body_Op(modelspace,"proton")
                                   + 1./Z * RpSpinOrbitCorrection(modelspace);
 }

 Operator Rn2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   if (Z==A) return 0.0*KineticEnergy_Op(modelspace);
   return R2CM_Op(modelspace) + (A-2.0)/(A*(A-Z))*R2_1body_Op(modelspace,"neutron")
                                   - 4./(A*(A-Z))*R2_2body_Op(modelspace,"neutron");
 }

 Operator Rm2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   return (1./A)*RSquaredOp(modelspace) - R2CM_Op(modelspace)  ;
 }



 // Evaluate <bra | r1*r2 | ket>, omitting the factor (hbar * omega) /(m * omega^2)
/// Returns the normalized, anti-symmetrized, J-coupled, two-body matrix element of \f$  \vec{r}_1\cdot\vec{r}_2 \f$.
/// Calculational details are similar to Calculate_p1p2().
 double Calculate_r1r2(ModelSpace& modelspace, Ket & bra, Ket & ket, int J)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
   int nb = ob.n;
   int nc = oc.n;
   int nd = od.n;

   int la = oa.l;
   int lb = ob.l;
   int lc = oc.l;
   int ld = od.l;

   double ja = oa.j2*0.5;
   double jb = ob.j2*0.5;
   double jc = oc.j2*0.5;
   double jd = od.j2*0.5;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   if (std::abs(fab-fcd)%2 >0) return 0; // p1*p2 only connects kets with delta N = 0,1
   if (std::abs(fab-fcd)>2) return 0; // p1*p2 only connects kets with delta N = 0,1

   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double r1r2=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( std::abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<std::abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
              // factor to account for antisymmetrization

              int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (std::abs(mosh_ab)<1e-8) continue;

              for (int N_cd=std::max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd and N_ab != N_cd) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (std::abs(mosh_cd)<1e-8) continue;

                double r2cm = 0;
                double r2rel = 0;
                if (n_ab == n_cd)
                {
                  if      (N_ab == N_cd)   r2cm = (2*N_ab+Lam_ab+1.5);
                  else if (N_ab == N_cd+1) r2cm = -sqrt(N_ab*( N_ab+Lam_ab+0.5));
                  else if (N_ab == N_cd-1) r2cm = -sqrt(N_cd*( N_cd+Lam_ab+0.5));
                }
                if (N_ab == N_cd)
                {
                  if      (n_ab == n_cd)   r2rel = (2*n_ab+lam_ab+1.5);
                  else if (n_ab == n_cd+1) r2rel = -sqrt(n_ab*( n_ab+lam_ab+0.5));
                  else if (n_ab == n_cd-1) r2rel = -sqrt(n_cd*( n_cd+lam_cd+0.5));
                }
                double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                r1r2 += (r2cm-r2rel) * prefactor;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   double oscillator_b2 = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
   // normalize and give dimension of length
   // the 0.5 comes from the virial theorem, <V> = <T> = 1/2 (N+3/2)hw
   r1r2 *= oscillator_b2 *0.5 / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
//   r1r2 *= 1.0 / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return r1r2 ;

 }




/// Center of mass Hamiltonian
/// \f{align*}{
/// H_{CM} &= T_{CM} + \frac{1}{2} Am\omega^2 R^2 \\
///        &= T_{CM} + \frac{1}{2b^2} AR^2 \hbar\omega
/// \f}
 Operator HCM_Op(ModelSpace& modelspace)
 {
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
//   double oscillator_b2 = HBARC*HBARC/M_NUCLEON/hw;
//   Operator HcmOp = TCM_Op(modelspace) + R2CM_Op(modelspace) * (0.5*A * hw / oscillator_b2);
   Operator HcmOp = TCM_Op(modelspace) + 0.5*A*M_NUCLEON*hw*hw/HBARC/HBARC * R2CM_Op(modelspace) ;
   std::cout << "HcmOp: first 1b element = " << HcmOp.OneBody(0,0) << std::endl;
   return HcmOp;
 }



/// Returns
/// \f[ r^2 = \sum_{i} r_{i}^2 \f]
///
Operator RSquaredOp(ModelSpace& modelspace)
{
   Operator r2 = Operator(modelspace);
   r2.OneBody.zeros();
   unsigned int norbits = modelspace.GetNumberOrbits();
   double oscillator_b2 = HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega();
   for (unsigned int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      r2.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for ( unsigned int b : r2.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
      {
        if ( b < a ) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              r2.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              r2.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           r2.OneBody(b,a) = r2.OneBody(a,b);
        }
      }
   }
//   r2.OneBody *= (HBARC*HBARC/M_NUCLEON/hw);
   r2.OneBody *= oscillator_b2;
   return r2;
}


 Operator R2_p1_Op(ModelSpace& modelspace)
 {
   return R2_1body_Op(modelspace,"proton");
 }

/// One-body part of the proton charge radius operator.
/// Returns
/// \f[ 
/// \hat{R}^{2}_{p1} = \sum_{i} e_{i}{r}_i^2
/// \f]
 Operator R2_1body_Op(ModelSpace& modelspace,std::string option)
 {
   Operator r2(modelspace);
   double oscillator_b = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());

   auto orbitlist = modelspace.proton_orbits;
   if (option == "neutron") orbitlist = modelspace.neutron_orbits;
   else if (option == "matter")  orbitlist.insert(orbitlist.end(),modelspace.neutron_orbits.begin(),modelspace.neutron_orbits.end());
   else if (option != "proton") std::cout << "!!! WARNING. BAD OPTION "  << option << " FOR imsrg_util::R2_p1_Op !!!" << std::endl;
 
   for (unsigned int a : orbitlist )
   {
      Orbit & oa = modelspace.GetOrbit(a);
      r2.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for ( unsigned int b : r2.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
      {
        if ( b < a ) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              r2.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              r2.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           r2.OneBody(b,a) = r2.OneBody(a,b);
        }
      }
   }
   r2.OneBody *= oscillator_b;
   return r2;
 }

 Operator R2_p2_Op(ModelSpace& modelspace)
 {
   return R2_2body_Op(modelspace,"proton");
 }

/// Two-body part of the proton charge radius operator.
/// Returns
/// \f[ 
/// \hat{R}^{2}_{p2} = \sum_{i\neq j} e_{i}\vec{r}_i\cdot\vec{r}_j 
/// \f]
/// evaluated in the oscillator basis.
 Operator R2_2body_Op(ModelSpace& modelspace,std::string option)
 {
   Operator Rp2Op(modelspace,0,0,0,2);
//   double oscillator_b = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int Tz = tbc.Tz;
      if (option=="proton" and Tz > 0) continue; // don't bother with nn channel
      if (option=="neutron" and Tz < 0) continue; // don't bother with pp channel
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         double prefactor = 1; // factor to account for double counting in pn channel.
         if (Tz==0 and (option=="proton" or option=="neutron")) prefactor = 0.5;
//         if (option=="proton" and bra.op->tz2>0) continue;
//         else if (option=="neutron" and bra.op->tz2<0) continue;
         if (option!="matter" and option!="proton" and option!="neutron") std::cout << "!!! WARNING. BAD OPTION "  << option << " FOR imsrg_util::R2_p2_Op !!!" << std::endl;
         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket & ket = tbc.GetKet(iket);
//            double mat_el = Calculate_r1r2(modelspace,bra,ket,tbc.J) * oscillator_b ; 
            double mat_el = Calculate_r1r2(modelspace,bra,ket,tbc.J) * prefactor; 
            Rp2Op.TwoBody.SetTBME(ch,ibra,iket,mat_el);
            Rp2Op.TwoBody.SetTBME(ch,iket,ibra,mat_el);
         }
      }
   }
   return Rp2Op;
 }


Operator ProtonDensityAtR(ModelSpace& modelspace, double R)
{
  Operator Rho(modelspace,0,0,0,2);
  for ( auto i : modelspace.proton_orbits)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    Rho.OneBody(i,i) = HO_density(oi.n,oi.l,modelspace.GetHbarOmega(),R);
  }
  return Rho;
}

Operator NeutronDensityAtR(ModelSpace& modelspace, double R)
{
  Operator Rho(modelspace,0,0,0,2);
  for ( auto i : modelspace.neutron_orbits)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    Rho.OneBody(i,i) = HO_density(oi.n,oi.l,modelspace.GetHbarOmega(),R);
  }
  return Rho;
}


Operator RpSpinOrbitCorrection(ModelSpace& modelspace)
{
  Operator dr_so(modelspace,0,0,0,2);
  double M2 = M_NUCLEON*M_NUCLEON/(HBARC*HBARC);
  int norb = modelspace.GetNumberOrbits();
  for (int i=0;i<norb;i++)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    double mu_i = oi.tz2<0 ? 1.79 : -1.91;
    int kappa = oi.j2 < 2*oi.l ? oi.l : -(oi.l+1);
    dr_so.OneBody(i,i) = -mu_i/M2*(kappa+1);
  }
  return dr_so;
}

// Electric monopole operator
/// Returns
/// \f[ r_{e}^2 = \sum_{i} e_{i} r_{i}^2 \f]
///
Operator E0Op(ModelSpace& modelspace)
{
   Operator e0(modelspace);
   e0.EraseZeroBody();
   e0.OneBody.zeros();
//   unsigned int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
   for (unsigned int a : modelspace.proton_orbits)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      e0.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for (unsigned int b : e0.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
        if (b<=a) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              e0.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              e0.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           e0.OneBody(b,a) = e0.OneBody(a,b);
        }
      }
   }
   e0.OneBody *= (HBARC*HBARC/M_NUCLEON/hw);
   return e0;
}





//Operator FourierBesselCoeff(ModelSpace& modelspace, int nu, double R)
// providing an index list allows us to select e.g. just protons or just neutrons
Operator FourierBesselCoeff(ModelSpace& modelspace, int nu, double R, std::vector<index_t> index_list)
{
  Operator a_nu(modelspace,0,0,0,2);

  size_t npoints = 99;  // should be a multiple of 3 so we can use the Simpson 3/8 rule
  std::vector<double> RGRID(npoints);
  std::vector<double> BESSEL(npoints);
  std::vector<double> PSI_p(npoints);
  std::vector<double> INTEGRAND(npoints);
  double dr = R / npoints;
  double Q = nu*M_PI/R;

  for (size_t i=0; i<npoints; i++)
  {
    RGRID[i] = i*dr;
    BESSEL[i] =  gsl_sf_bessel_j0( Q * RGRID[i] );
  }

  double hw = modelspace.GetHbarOmega();

  for ( auto p : index_list )  
  {
    Orbit& op = modelspace.GetOrbit(p);
    for (size_t i=0; i<npoints; i++)
    {
      PSI_p[i] = HO_Radial_psi( op.n, op.l, hw, RGRID[i]);
    }

    // Use Simpson's rule to perform the integral
    // int  dr  r^2  psi_p(r)  j0(nu *pi * r /R ) psi_q(r)
//    size_t q = p;
    for (auto q : a_nu.OneBodyChannels.at({op.l,op.j2,op.tz2}) )
    {
      if (q>p) continue;
      Orbit& oq = modelspace.GetOrbit(q);
      for (size_t i=0; i<npoints; i++)
      {
        INTEGRAND[i] = HO_Radial_psi( oq.n, oq.l, hw, RGRID[i] ) * RGRID[i]*RGRID[i] * BESSEL[i] * PSI_p[i];
      }

      double integral = INTEGRAND[0] + INTEGRAND[npoints-1];
      for (size_t i=1; i<npoints; i+=3)
      {
        integral += 3*INTEGRAND[i] + 3*INTEGRAND[i+1] + 2*INTEGRAND[i+2];
      }
      integral *= 3./8 * dr;
      integral *= Q*Q;

      a_nu.OneBody(p,q) = integral;
      a_nu.OneBody(q,p) = integral;
    }
  }

  return a_nu; // There may need to be some additional normalization by pi or something...
}

//struct FBCIntegrandParameters{int n; int l; double hw;};
//
//double FBCIntegrand(double x, void *p)
//{
//  struct FBCIntegrandParameters * params = (struct FBCIntegrandParameters *)p;
//  return x*HO_density(params->n, params->l, params->hw, x);
//}
//
//Operator FourierBesselCoeff(ModelSpace& modelspace, int nu, double R, std::vector<index_t> index_list)
//{
//  Operator a_nu(modelspace,0,0,0,2);
//  double omega = nu * M_PI / R; // coefficient of sine function, i.e. sin(omega*x)
//  double L = R; // range of integration
//  size_t n = 20; // number of bisections
//  gsl_integration_qawo_table * table = gsl_integration_qawo_table_alloc (omega, L, GSL_INTEG_SINE, n);
//  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc ( n );
//
//  const double epsabs = 1e-5; // std::absolute error
//  const double epsrel = 1e-5; // relative error
//  const size_t limit = n; // maximum number of subintervals (maybe should be different?)
//  const double start = 0.0; // lower limit on integration range
//  double result;
//  double abserr;
//  gsl_function F;
//  F.function = &FBCIntegrand;
//
//  for (auto i : index_list )
//  {
//    Orbit& oi = modelspace.GetOrbit(i);
//    struct FBCIntegrandParameters params = {oi.n, oi.l, modelspace.GetHbarOmega()};
//    F.params = &params;
//    //int status = gsl_integration_qawo (&F, start, epsstd::abs, epsrel, limit, workspace, table, &result, &std::abserr);
//    gsl_integration_qawo (&F, start, epsabs, epsrel, limit, workspace, table, &result, &abserr);
//    a_nu.OneBody(i,i) = M_PI*M_PI/R/R/R * R/nu/M_PI*(result);
//    std::cout << "orbit,nu = " << i << "," << nu << "  => " << a_nu.OneBody(i,i) << "  from " << result << " (" << abserr << ")" << std::endl;
//  }
//  return a_nu;
//}



/// Returns the \f$ T^{2} \f$ operator
 Operator Isospin2_Op(ModelSpace& modelspace)
 {
   Operator T2 = Operator(modelspace,0,0,0,2);
   T2.OneBody.diag().fill(0.75);

   for (int ch=0; ch<T2.nChannels; ++ch)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
     arma::mat& TB = T2.TwoBody.GetMatrix(ch);
     // pp,nn:  2<t2.t1> = 1/(2(1+delta_ab)) along diagonal
     if (std::abs(tbc.Tz) == 1)
     {
        TB.diag().fill(0.5); // pp,nn TBME's
         // I managed to confuse myself.  t_{pp'pp'} should give 2 t1*t2. This gives reasonable results in nushellx.
//        for (int ibra=0;ibra<tbc.GetNumberKets(); ++ibra)
//        {
//           Ket& bra = tbc.GetKet(ibra);
//           if (bra.p == bra.q)
//           {
//             TB(ibra,ibra) /= 2.;
//           }
//        }
     }
     else if (tbc.Tz == 0)
     {
        for (size_t ibra=0;ibra<tbc.GetNumberKets(); ++ibra)
        {
           Ket& bra = tbc.GetKet(ibra);
           Orbit& oa = modelspace.GetOrbit(bra.p);
           Orbit& ob = modelspace.GetOrbit(bra.q);
           for (size_t iket=ibra;iket<tbc.GetNumberKets(); ++iket)
           {
             Ket& ket = tbc.GetKet(iket);
             Orbit& oc = modelspace.GetOrbit(ket.p);
             Orbit& od = modelspace.GetOrbit(ket.q);
             if (oa.j2==oc.j2 and oa.n==oc.n and oa.l==oc.l
               and ob.j2==od.j2 and ob.n==od.n and ob.l==od.l )
             {
               // tz1 tz2 case
               if( oa.tz2 == oc.tz2 and ob.tz2==od.tz2)
               {
                 TB(ibra,iket) -= 0.5;
               }
               // t+ t- case
               if( oa.tz2 == od.tz2 and ob.tz2==oc.tz2)
               {
                 TB(ibra,iket) += 1.0;
               }
               // if a==b==c==d, we need to consider the exchange term
               if (oa.j2==ob.j2 and oa.n==ob.n and oa.l==ob.l)
               {
                  int phase = bra.Phase(tbc.J);
                  // tz1 tz2 case
                  if( oa.tz2 == oc.tz2 and ob.tz2==od.tz2)
                  {
                    TB(ibra,iket) += phase * 1.0;
                  }
                  // t+ t- case
                  if( oa.tz2 == od.tz2 and ob.tz2==oc.tz2)
                  {
                    TB(ibra,iket) -= phase * 0.5;
                  }
               }
               TB(iket,ibra) = TB(ibra,iket); // hermitian
             }

           }
        }
     }
   }
   return T2;
 }



  /// Returns a reduced electric multipole operator with units \f$ e\f$ fm\f$^{\lambda} \f$
  /// See Suhonen eq. (6.23)
  Operator ElectricMultipoleOp(ModelSpace& modelspace, int L)
  {
    Operator EL(modelspace, L,0,L%2,2);
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*L); // b^L where b=sqrt(hbar/mw)
    for (int i : modelspace.proton_orbits)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,L) * bL ;
        EL.OneBody(i,j) = modelspace.phase(jj+L-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(ji,jj, L, 0.5, -0.5,0) * r2int;
        EL.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * EL.OneBody(i,j);
      }
    }
    return EL;
  }

  /// Returns a reduced electric multipole operator with units \f$ e\f$ fm\f$^{\lambda} \f$
  /// See Suhonen eq. (6.23)
  Operator NeutronElectricMultipoleOp(ModelSpace& modelspace, int L)
  {
    Operator EL(modelspace, L,0,L%2,2);
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*L); // b^L where b=sqrt(hbar/mw)
    for (int i : modelspace.neutron_orbits)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,L) * bL ;
        EL.OneBody(i,j) = modelspace.phase(jj+L-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(ji,jj, L, 0.5, -0.5,0) * r2int;
        EL.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * EL.OneBody(i,j);
      }
    }
    return EL;
  }

  /// Returns a reduced magnetic multipole operator with units \f$ \mu_{N}\f$ fm\f$ ^{\lambda-1} \f$
  Operator MagneticMultipoleOp(ModelSpace& modelspace, int L)
  {
    return MagneticMultipoleOp_pn(modelspace,L,"both");
  }

  /// Returns a reduced magnetic multipole operator with units \f$ \mu_{N}\f$ fm\f$ ^{\lambda-1} \f$
  /// This version allows for the selection of just proton or just neutron contributions, or both.
  /// See Suhonen eq. (6.24)
  Operator MagneticMultipoleOp_pn(ModelSpace& modelspace, int L, std::string pn)
  {
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*(L-1));
    Operator ML(modelspace, L,0,(L+1)%2,2);
    if (L<1)
    {
      std::cout << "A magnetic monopole operator??? Setting it to zero..." << std::endl;
      return ML;
    }
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      if (pn=="proton" and oi.tz2>0) continue;
      if (pn=="neutron" and oi.tz2<0) continue;
      double gl = oi.tz2<0 ? 1.0 : 0.0;
      double gs = oi.tz2<0 ? 5.586 : -3.826;
      double ji = 0.5*oi.j2;
      for ( int j : ML.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        // multiply radial integral by b^L-1 = (hbar/mw)^L-1/2
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,L-1) * bL;
        double kappa =  modelspace.phase(oi.l+ji+0.5) * (ji+0.5)  +  modelspace.phase(oj.l+jj+0.5) * (jj+0.5);
        ML.OneBody(i,j) = modelspace.phase(jj+L-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(ji, jj,L, 0.5,-0.5,0)
                        * (L - kappa) *(gl*(1+kappa/(L+1.))-0.5*gs )  * r2int ;
        ML.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * ML.OneBody(i,j);
      }
    }
    return ML;
  }



  /// Returns a reduced electric multipole operator with units \f$ e\f$ fm\f$^{\lambda} \f$
  /// See Suhonen eq. (6.23)
  Operator IntrinsicElectricMultipoleOp(ModelSpace& modelspace, int L)
  {
    Operator EL(modelspace, L,0,L%2,2);
    double bL = pow( HBARC*HBARC/(0.5*M_NUCLEON)/modelspace.GetHbarOmega(),0.5*L); // b^L where b=sqrt(hbar/mw)

    int emax = modelspace.GetEmax();
    int norb_rel = (2*emax+1)*(2*emax+2)/2;
    arma::mat ME_rel(norb_rel, norb_rel);
    for (int nrela=0; 2*nrela<2*emax; nrela+=1)
    {
     for (int lrela=0; lrela<2*emax-2*nrela; lrela+=1)
     {
       int ia = (2*nrela + lrela)*(2*nrela+lrela+1)/2 + lrela;
       for (int lrelb=std::abs(lrela-L); lrelb<=lrela+L; lrelb+=2)
       {
         for (int nrelb=0; 2*nrelb<=2*emax-lrelb; nrelb+=1)
         {
           int ib = (2*nrelb + lrelb)*(2*nrelb+lrelb+1)/2 + lrelb;
           double r2int = RadialIntegral(nrela,lrela,nrelb,lrelb,L) * bL ;
           double angint = modelspace.phase(lrela) * sqrt( (2*lrela+1)*(2*lrelb+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(lrela, L, lrelb, 0,0,0);
           ME_rel(ia,ib) = angint * r2int;
           printf("< %d %d || E2 || %d %d > = %f   angular: %f \n",nrela,lrela,nrelb,lrelb,ME_rel(ia,ib), angint );
         }
       }
     }
    }



    // Moshinsky transformation to lab frame bra and ket
    // lab frame: < na la ja  nb lb jb  Jab || EL || nc lc jc  nd ld jd  Jcd>
    for ( auto& itmat : EL.TwoBody.MatEl )
    {
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(itmat.first[0]);
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(itmat.first[1]);
      if (tbc_bra.Tz>=0) continue;
      int Jab = tbc_bra.J;
      int Jcd = tbc_ket.J;
      for (size_t ibra=0; ibra<tbc_bra.GetNumberKets(); ++ibra)
      {
       Ket& bra = tbc_bra.GetKet(ibra);
       int na = bra.op->n;
       int la = bra.op->l;
       double ja = bra.op->j2*0.5;
       int nb = bra.oq->n;
       int lb = bra.oq->l;
       double jb = bra.oq->j2*0.5;
       int rho_ab = 2*na+2*nb+la+lb;
       for (size_t iket=0; iket<tbc_ket.GetNumberKets(); ++iket)
       {
        Ket& ket = tbc_bra.GetKet(ibra);
        int nc = ket.op->n;
        int lc = ket.op->l;
        double jc = ket.op->j2*0.5;
        int nd = ket.oq->n;
        int ld = ket.oq->l;
        double jd = ket.oq->j2*0.5;
        int rho_cd = 2*nc+2*nd+lc+ld;

        double labME = 0;
        
        for (int Lab=std::max(std::abs(la-lb),std::abs(Jab-1)); Lab<=std::min(la+lb,Jab+1); Lab+=1)
        {
         for (int Lcd=std::max(std::abs(lc-ld),std::abs(Jcd-1))+(Lab+L)%2; Lcd<=std::min(lc+ld,Jcd+1); Lcd+=2)
         {
          for (int S=std::max(std::abs(Jab-Lab),std::abs(Jcd-Lcd)); S<=std::min(1,std::min(Jab+Lab,Jcd+Lcd)); S+=1 )
          {
            double njab = AngMom::NormNineJ( la, 0.5, ja, lb, 0.5, jb, Lab, S, Jab);
            double njcd = AngMom::NormNineJ( lc, 0.5, jc, ld, 0.5, jd, Lcd, S, Jcd);
            if (std::abs(njab)<1e-8 or std::abs(njcd)<1e-8) continue;
            for ( int nab=0; nab<=rho_ab; nab+=1)
            {
             for ( int lab=0; lab<=rho_ab-2*nab; lab+=1)
             {
              int iab = (2*nab + lab)*(2*nab+lab+1)/2 + lab;
              for (int Lam=std::abs(lab-Lab); Lam<=std::min(lab+Lab,rho_ab-2*nab-lab); Lam+=2)
              {
               int N = rho_ab - 2*nab - lab - Lam;
               double moshab = modelspace.GetMoshinsky( N, Lam, nab, lab, na,la,nb,lb,Lab);
               if (std::abs(moshab)<1e-8) continue;
               for ( int ncd=0; ncd<=rho_cd; ncd+=1)
               {
                for ( int lcd=std::abs(Lam-Lcd)+(Lcd+Lam)%2; lcd<=std::min(rho_cd-2*ncd,Lab+Lcd); lcd+=2)
                {
                 int icd = (2*ncd + lcd)*(2*ncd+lcd+1)/2 + lcd;
                 double moshcd = modelspace.GetMoshinsky( N, Lam, ncd, lcd, nc,lc,nd,ld,Lcd);
                 labME += njab * njcd * moshab * moshcd * ME_rel(iab,icd);
                }
               }
              }
             }
            }
          }
         }
        }
        itmat.second(ibra,iket) = labME;
       }
      }
    }
    
    std::cout << "done with intrinsic EL. one body = " <<  std::endl;
    return EL;
  }



/// Evaluate the radial integral \f[
/// \tilde{\mathcal{R}}^{\lambda}_{ab} = \int_{0}^{\infty} dx \tilde{g}_{n_a\ell_a}(x)x^{\lambda+2}\tilde{g}_{n_b\ell_b}(x)
/// \f]
/// where \f$ \tilde{g}(x) \f$ is the radial part of the harmonic oscillator wave function with unit oscillator length \f$ b=1 \f$
/// and \f$ x = r/b \f$.
/// To obtain the radial integral for some other oscillator length, multiply by \f$ b^{\lambda} \f$.
/// This implementation uses eq (6.41) from Suhonen.
/// Note this is only valid for \f$ \ell_a+\ell_b+\lambda\f$ = even.
/// If \f$ \ell_a+\ell_b+\lambda\f$ is odd, RadialIntegral_RpowK() is called.
  double RadialIntegral(int na, int la, int nb, int lb, int L)
  {
    if ((la+lb+L)%2!=0) return RadialIntegral_RpowK(na,la,nb,lb,L);
    int tau_a = std::max((lb-la+L)/2,0);
    int tau_b = std::max((la-lb+L)/2,0);
    int sigma_min = std::max(std::max(na-tau_a,nb-tau_b),0);
    int sigma_max = std::min(na,nb);
  
    double term1 = AngMom::phase(na+nb) * gsl_sf_fact(tau_a)*gsl_sf_fact(tau_b) * sqrt(gsl_sf_fact(na)*gsl_sf_fact(nb)
                   / (tgamma(na+la+1.5)*tgamma(nb+lb+1.5) ) );
//    double term1 = AngMom::phase(na+nb) * gsl_sf_fact(tau_a)*gsl_sf_fact(tau_b) * sqrt(gsl_sf_fact(na)*gsl_sf_fact(nb)
//                   / (gsl_sf_gamma(na+la+1.5)*gsl_sf_gamma(nb+lb+1.5) ) );
    double term2 = 0;
    for (int sigma=sigma_min; sigma<=sigma_max; ++sigma)
    {
      term2 += gsl_sf_gamma(0.5*(la+lb+L)+sigma+1.5) / (gsl_sf_fact(sigma)*gsl_sf_fact(na-sigma)*gsl_sf_fact(nb-sigma)*gsl_sf_fact(sigma+tau_a-na)*gsl_sf_fact(sigma+tau_b-nb) );
    }
    return term1*term2;
  
  }


 // the nomenclature on the variable k has sort of gotten out of hand here...
 double RadialIntegral_RpowK(int na, int la, int nb, int lb, int k)
 {
   long double I = 0;
   int pmin = (la+lb)/2;
   int pmax = pmin + na + nb;

   std::vector<double> Ip(pmax+1);
   Ip[pmin] = TalmiI(pmin,k);
   for (int p=pmin+1; p<=pmax;++p) Ip[p] = (2*p+1+k)/(2*p+1.) * Ip[p-1];

//   if (pmax < -10 )
   {
     for (int p=pmin;p<=pmax;++p)
     {
//        I += TalmiB(na,la,nb,lb,p) * Ip[p];
        I += TalmiB(na,la,nb,lb,p) * TalmiI(p,k);
//        I += AngMom::TalmiB(na,la,nb,lb,p) * TalmiI(p,k);
     }
   }
//   else // just do it with quadrature
   {
     double Iquad = 0;
     double Norm = 2*sqrt( gsl_sf_gamma(na+1)*gsl_sf_gamma(nb+1)/gsl_sf_gamma(na+la+1.5)/gsl_sf_gamma(nb+lb+1.5)) ;
//     double Norm1 = sqrt( 2*gsl_sf_gamma(na+1)/gsl_sf_gamma(na+la+1.5)) ;
//     double Norm2 = sqrt( 2*gsl_sf_gamma(nb+1)/gsl_sf_gamma(nb+lb+1.5)) ;
//     double Norm = 2.0 * exp( 0.5*( lgamma(na+1) + lgamma(nb+1) - lgamma(na+la+1.5) - lgamma(nb+lb+1.5)) ) ;
//     double Norm = 2.0 * exp( 0.5*( lgamma(na+1) + lgamma(nb+1) - lgamma(na+la+1.5) - lgamma(nb+lb+1.5)) ) ;
//     double Norm = 2*sqrt( gsl_sf_fact(na) * pow(2,na+la) / M_SQRTPI / gsl_sf_doublefact(2*na+2*la+1)  );
//     Norm *=       2*sqrt( gsl_sf_fact(nb) * pow(2,nb+lb) / M_SQRTPI / gsl_sf_doublefact(2*nb+2*lb+1)  );
//     int npoints = std::min(50,2*na+2*nb+la+lb+20);
     int npoints = 50;
//     double I1=0;
//     double I2 = 0;
     for (int i=0;i<npoints;i++)
     {
       double x_i = GaussLaguerre::gauss_laguerre_points[npoints][i][0];
       double w_i = GaussLaguerre::gauss_laguerre_points[npoints][i][1];
       double f_i = Norm * exp(-x_i*x_i) * gsl_sf_laguerre_n(na,la+0.5,x_i*x_i)  * gsl_sf_laguerre_n(nb,lb+0.5,x_i*x_i) * pow(x_i,la+lb+2+k)  ;
//       double f1 = Norm1*Norm1 * exp(-x_i*x_i) * gsl_sf_laguerre_n(na,la+0.5,x_i*x_i) * gsl_sf_laguerre_n(na,la+0.5,x_i*x_i) * pow(x_i,la+la+2);
//       double f2 = Norm2*Norm2 * exp(-x_i*x_i) * gsl_sf_laguerre_n(nb,lb+0.5,x_i*x_i) * gsl_sf_laguerre_n(nb,lb+0.5,x_i*x_i) * pow(x_i,lb+lb+2);
       Iquad += w_i * f_i;
//       I1 += w_i * f1;
//       I2 += w_i * f2;
     }
//     std::cout << " " << na << " " << la << " " << nb << " " << lb << " " << k << "   :  " << I << "  ,  " << Iquad  << "    " << Iquad / I << "  norms   " << I1 << " " << I2 << std::endl; 
   }
//   else // for large p, we need to be more careful. TODO: This doesn't work yet.
//   {
//     int q = (la+lb)/2;
//     int Kmin = std::max(0, pmin-q-nb);
//     int Kmax = std::min(na, pmax-q);
//     for (int K=Kmin;K<=Kmax;++K)
//     {
//       int inner_pmin = K+q;
//       int inner_pmax = K+q+nb;
//       for (int p=inner_pmin;p<=inner_pmax;p+=2)
//       {
//         if (p+1<=inner_pmax)
//         {
//           I += TalmiB_SingleTermPair(na,la,nb,lb,p,K,k) * Ip[p]; // this is still giving nonsense
//         }
//         else
//         {
//           I += TalmiB_SingleTerm(na,la,nb,lb,p,K) * Ip[p];
//         }
//       }
//     }
//   }

   return I;
 }

/// General Talmi integral for a potential r**k
/// 1/gamma(p+3/2) * 2*INT dr r**2 r**2p r**k exp(-r**2/b**2)
/// This is valid for (2p+3+k) > 0. The Gamma function diverges for non-positive integers.
 long double TalmiI(int p, double k)
 {
//   return gsl_sf_gamma(p+1.5+0.5*k) / gsl_sf_gamma(p+1.5);
//   return boost::math::tgamma_ratio(p+1.5+0.5*k, p+1.5);
     return exp( lgamma(p+1.5+0.5*k)-lgamma(p+1.5) ); // Tested that this agrees with the boost version to better than 1 in 10^12. Also stable against overflow.
 }

/// Calculate B coefficient for Talmi integral. Formula given in Brody and Moshinsky
/// "Tables of Transformation Brackets for Nuclear Shell-Model Calculations"
 long double TalmiB(int na, int la, int nb, int lb, int p)
 {
   if ( (la+lb)%2>0 ) return 0;
   
   int q = (la+lb)/2;

   if ( std::max(na+la+p, nb+lb+p) < 10 )
   {
     return AngMom::TalmiB( na, la, nb, lb, p);
   }
//   double B1 = AngMom::phase(p-q) * exp(lgamma(2*p+2)-lgamma(p+1)) / pow(2,(na+nb))
//               *  exp(0.5*(lgamma(na+1)+lgamma(nb+1)-lgamma(na+la+1)-lgamma(nb+lb+1) 
//   double B1 = AngMom::phase(p-q) * exp(lgamma(2*p+2)-lgamma(p+1) +0.5*(  lgamma(na+1)+lgamma(nb+1)-lgamma(na+la+1)-lgamma(nb+lb+1) + lgamma(2*na+2*la+2) + lgamma(2*nb+2*lb+2)) - (na+nb)*log(2)   ) ;
   long double logB1 = (lgamma(2*p+2)-lgamma(p+1) +0.5*(  lgamma(na+1)+lgamma(nb+1)-lgamma(na+la+1)-lgamma(nb+lb+1) + lgamma(2*na+2*la+2) + lgamma(2*nb+2*lb+2)) - (na+nb)*LOG2   ) ;
//               *  exp(0.5*(lgamma(na+1)+lgamma(nb+1)-lgamma(na+la+1)-lgamma(nb+lb+1) 
//                     + lgamma(2*na+2*la+2) + lgamma(2*nb+2*lb+2)) );
//               * sqrt(  exp(lgamma(na+1)) * exp(lgamma(nb+1)) / exp(lgamma(na+la+1))/ exp(lgamma(nb+lb+1))
//                     * exp(lgamma(2*na+2*la+2)) * exp(lgamma(2*nb+2*lb+2)) );
//              * sqrt( gsl_sf_fact(na)*gsl_sf_fact(nb)/gsl_sf_fact(na+la)/gsl_sf_fact(nb+lb) 
//                   * gsl_sf_fact(2*na+2*la+1) * gsl_sf_fact(2*nb+2*lb+1) );
   
   long double B2 = 0;
   int kmin = std::max(0, p-q-nb);
   int kmax = std::min(na, p-q);
   for (int k=kmin;k<=kmax;++k)
   {
//      B2  += exp(lgamma(la+k+1)-lgamma(k+1) +lgamma(p-(la-lb)/2-k+1) -lgamma(2*p-la+lb-2*k+2) )
//             / (  gsl_sf_fact(2*la+2*k+1) * gsl_sf_fact(na-k)  
//                * gsl_sf_fact(nb - p + q + k) * gsl_sf_fact(p-q-k) );
      B2  += exp(logB1+lgamma(la+k+1)-lgamma(k+1) +lgamma(p-(la-lb)/2-k+1) -lgamma(2*p-la+lb-2*k+2) 
                - lgamma(2*la+2*k+2) -lgamma(na-k+1)  
                - lgamma(nb - p + q + k+1) - lgamma(p-q-k+1) );
   }
   
   return AngMom::phase(p-q) *  B2;
//   return  B1 * B2;
 }


 long double TalmiB_SingleTerm(int na, int la, int nb, int lb, int p, int K)
 {

   if ( (la+lb)%2>0 ) return 0;
   
   int q = (la+lb)/2;
   long double logB1 = (lgamma(2*p+2)-lgamma(p+1) +0.5*(  lgamma(na+1)+lgamma(nb+1)-lgamma(na+la+1)-lgamma(nb+lb+1) + lgamma(2*na+2*la+2) + lgamma(2*nb+2*lb+2)) - (na+nb)*LOG2   ) ;
//   long double B2 = 0;
   int k = K;
//   int kmin = std::max(0, p-q-nb);
//   int kmax = std::min(na, p-q);
//   for (int k=kmin;k<=kmax;++k)
//   {
//      B2  += exp(lgamma(la+k+1)-lgamma(k+1) +lgamma(p-(la-lb)/2-k+1) -lgamma(2*p-la+lb-2*k+2) )
//             / (  gsl_sf_fact(2*la+2*k+1) * gsl_sf_fact(na-k)  
//                * gsl_sf_fact(nb - p + q + k) * gsl_sf_fact(p-q-k) );
   long double B2  = exp(logB1+lgamma(la+k+1)-lgamma(k+1) +lgamma(p-(la-lb)/2-k+1) -lgamma(2*p-la+lb-2*k+2) 
                - lgamma(2*la+2*k+2) -lgamma(na-k+1)  
                - lgamma(nb - p + q + k+1) - lgamma(p-q-k+1) );
//   }

   return AngMom::phase(p-q) *  B2;
 }


 double factorial_ratio( int a, int b )
 {
   double c = 1.0;
   for (int x=std::min(a,b)+1; x<=std::max(a,b); x++ ) c*=x;
   return (a>b) ? c : 1.0/c;
 }


 long double TalmiB_SingleTermPair(int na, int la, int nb, int lb, int p, int K, int nu)
 {

   if ( (la+lb)%2>0 ) return 0;
   
   int q = (la+lb)/2;
   int k = K;
//   long double logB1 = (lgamma(2*p+2)-lgamma(p+1) +0.5*(  lgamma(na+1)+lgamma(nb+1)-lgamma(na+la+1)-lgamma(nb+lb+1) + lgamma(2*na+2*la+2) + lgamma(2*nb+2*lb+2)) - (na+nb)*LOG2   ) ;
   long double logB1 = ( 0.5*(  lgamma(na+1)+lgamma(nb+1)-lgamma(na+la+1)-lgamma(nb+lb+1) + lgamma(2*na+2*la+2) + lgamma(2*nb+2*lb+2)) - (na+nb)*LOG2   ) ;
   logB1 +=  lgamma(la+K+1) - lgamma(k+1) - lgamma(2*la+2*k+2) - lgamma(na-k+1);
   logB1 +=  -lgamma(p+1) -lgamma(na+k+(la+lb)/2-p+1);

   long double B2 = factorial_ratio( 2*p+1, 2*p+1-la-lb-2*k)  *  factorial_ratio( p-k-(la-lb)/2,  p-k-(la+lb)/2);
   long double B3 = 1.0 - (2*p+3)*(nb+k+(la+lb)/2-p)/(2*p+3-la+lb-2*k)/(p+1-k-(la+lb)/2) * (2*p+3+nu)/(p+3);
//   long double B2 = 0;
//   int kmin = std::max(0, p-q-nb);
//   int kmax = std::min(na, p-q);
//   for (int k=kmin;k<=kmax;++k)
//   {
//      B2  += exp(lgamma(la+k+1)-lgamma(k+1) +lgamma(p-(la-lb)/2-k+1) -lgamma(2*p-la+lb-2*k+2) )
//             / (  gsl_sf_fact(2*la+2*k+1) * gsl_sf_fact(na-k)  
//                * gsl_sf_fact(nb - p + q + k) * gsl_sf_fact(p-q-k) );
//   long double B2  = exp(logB1+lgamma(la+k+1)-lgamma(k+1) +lgamma(p-(la-lb)/2-k+1) -lgamma(2*p-la+lb-2*k+2) 
//                - lgamma(2*la+2*k+2) -lgamma(na-k+1)  
//                - lgamma(nb - p + q + k+1) - lgamma(p-q-k+1) );
//   }

   return AngMom::phase(p-q) *  exp(logB1)*B2*B3;
 }



  Operator AllowedFermi_Op(ModelSpace& modelspace)
  {
    Operator Fermi(modelspace,0,1,0,2);
    Fermi.SetHermitian();
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      for (int j : Fermi.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        Orbit& oj = modelspace.GetOrbit(j);
        if (oi.n!=oj.n or oi.tz2 == oj.tz2) continue;
        Fermi.OneBody(i,j) = sqrt(oi.j2+1.0);  // Reduced matrix element
      }
    }
    return Fermi;
  }

/// Note that there is a literature convention to include the 1/sqrt(Lambda) factor
/// in the reduced matrix element rather than in the expression involving the sum
/// over one-body densities (see footnote on pg 165 of Suhonen).
/// I do not follow this convention, and instead produce the reduced matrix element
///  \f[ \langle f \| \sigma \tau_{\pm} \| i \rangle \f]
///
  Operator AllowedGamowTeller_Op(ModelSpace& modelspace)
  {
    Operator GT(modelspace,1,1,0,2);
    GT.SetHermitian();
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      for (int j : GT.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        Orbit& oj = modelspace.GetOrbit(j);
        if (oi.n!=oj.n or oi.l != oj.l or oi.tz2==oj.tz2) continue;
        double sixj = modelspace.GetSixJ(0.5,0.5,1.0,oj.j2/2.,oi.j2/2.,oi.l);
        double M_gt = 2 * modelspace.phase(oi.l+oi.j2/2.0+1.5) * sqrt((oi.j2+1)*(oj.j2+1)) * sqrt(1.5) * sixj;
        GT.OneBody(i,j) = M_gt;
      }
    }
    return GT;
  }



 /// Pauli spin operator \f[ \langle f \| \sigma \| i \rangle \f]
 Operator Sigma_Op(ModelSpace& modelspace)
 {
   return Sigma_Op_pn(modelspace,"both");
 }

 /// Pauli spin operator \f[ \langle f \| \sigma \| i \rangle \f]
 Operator Sigma_Op_pn(ModelSpace& modelspace, std::string pn)
 {
   Operator Sig(modelspace,1,0,0,2);
   Sig.SetHermitian();
   size_t norbits = modelspace.GetNumberOrbits();
   for (size_t i=0; i<norbits; ++i)
   {
     Orbit& oi = modelspace.GetOrbit(i);
      if (pn=="proton" and oi.tz2>0) continue;
      if (pn=="neutron" and oi.tz2<0) continue;
      for (auto j : Sig.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        Orbit& oj = modelspace.GetOrbit(j);
        if ((oi.n!=oj.n) or (oi.l != oj.l) or (oi.tz2!=oj.tz2)) continue;
        double sixj = modelspace.GetSixJ(0.5,0.5,1.0,oj.j2/2.,oi.j2/2.,oi.l);
        double M_sig = 2 * modelspace.phase(oi.l+oi.j2/2.0+1.5) * sqrt((oi.j2+1)*(oj.j2+1)) * sqrt(1.5) * sixj;
        Sig.OneBody(i,j) = M_sig;
      }
   } 
   return Sig;
 }

  void Reduce(Operator& X)
  {
    ModelSpace* modelspace = X.GetModelSpace();
    size_t norbits = modelspace->GetNumberOrbits();
    for (size_t i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace->GetOrbit(i);
      for ( auto j : X.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         X.OneBody(i,j) *= sqrt(oi.j2+1.);
      }
    }

    for ( auto& itmat : X.TwoBody.MatEl )
    {
      size_t ch_bra = itmat.first[0];
      int J = modelspace->GetTwoBodyChannel(ch_bra).J;
      itmat.second *= sqrt(2*J+1.);
    }
  }

  void UnReduce(Operator& X)
  {
    ModelSpace* modelspace = X.GetModelSpace();
    size_t norbits = modelspace->GetNumberOrbits();
    for (size_t i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace->GetOrbit(i);
      for ( auto j : X.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         X.OneBody(i,j) /= sqrt(oi.j2+1.);
      }
    }

    for ( auto& itmat : X.TwoBody.MatEl )
    {
      size_t ch_bra = itmat.first[0];
      int J = modelspace->GetTwoBodyChannel(ch_bra).J;
      itmat.second /= sqrt(2*J+1.);
    }

  }


  void SplitUp(Operator& OpIn, Operator& OpLow, Operator& OpHi, int ecut)
  {
    ModelSpace* modelspace = OpIn.GetModelSpace();
    OpLow = OpIn;
    size_t norbits = modelspace->GetNumberOrbits();
    size_t ncut = 0;
    for (size_t i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace->GetOrbit(i);
      if ( (2*oi.n + oi.l) > ecut)
      {
         ncut = i;
         break;
      }
    }

    for (size_t i=ncut; i<norbits; ++i)
    {
      for (size_t j=0; j<norbits; ++j)
      {
         OpLow.OneBody(i,j) = 0;
         OpLow.OneBody(j,i) = 0;
      }
    }

    for (auto& itmat : OpLow.TwoBody.MatEl)
    {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel( itmat.first[0] );
      size_t nkets = tbc.GetNumberKets();
      for (size_t ibra=0;ibra<nkets;++ibra)
      {
       Ket& bra = tbc.GetKet(ibra);
       for (size_t iket=ibra;iket<nkets;++iket)
       {
         Ket& ket = tbc.GetKet(iket);
         if ( bra.p>ncut or bra.q>ncut or ket.p>ncut or ket.q>ncut)
         {
           itmat.second(ibra,iket) = 0;
         }
       }
      }
    }
    OpHi = OpIn - OpLow;
  }


  Operator RadialOverlap(ModelSpace& modelspace)
  {
     Operator OVL(modelspace,0,1,0,1);
     index_t norb = modelspace.GetNumberOrbits();
     for (index_t i=0; i<norb; ++i)
     {
       Orbit& oi = modelspace.GetOrbit(i);
       for (index_t j=0; j<norb; ++j)
       {
         Orbit& oj = modelspace.GetOrbit(j);
         OVL.OneBody(i,j) = RadialIntegral(oi.n, oi.l, oj.n, oj.l, 0 ); // This is not quite right. Only works for li+lj=even.
       } 
     }
     return OVL;
  }


  Operator LdotS_Op(ModelSpace& modelspace)
  {
    Operator LdS(modelspace,0,0,0,2);
    index_t norbits = modelspace.GetNumberOrbits();
    for (index_t i=0;i<norbits;++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      int li = oi.l;
      LdS.OneBody(i,i) = 0.5* (ji*(ji+1) -li*(li+1)-0.75);
    }
    return LdS;
  }


 Operator LCM_Op(ModelSpace& modelspace)
 {

   Operator LCM(modelspace, 1,0,0,2);
/*

   unsigned int norb = modelspace.GetNumberOrbits();
   for (unsigned int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      for (auto j : LCM.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         if (j<i) continue;
         Orbit & oj = modelspace.GetOrbit(j);
         if (oi.n != oj.n or oi.l != oj.l) continue ;
         double Lij = modelspace.phase(oi.l + (oj.j2+3)/2) * sqrt(oi.j2+1)*sqrt(oj.j2+1)  * sqrt(oi.l*(oi.l+1)*(2*oi.l+1)) * modelspace.GetSixJ(oi.l,oj.l,1,oj.j2/2.,oi.j2/2.,0.5);
         LCM.OneBody(i,j) = Lij;
         LCM.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * Lij;
      }
   }

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();

*/
  return LCM;
 }


 // < ij J || Q * Q || kl J > where Q is the quadrupole operator (possibly up to overall factors like square roots of pi, etc...)
 // < ij J || Q*Q || kl J > = <i||Q||l> <j||Q||k> (2J+1)/sqrt(5) (-1)^(jk-jj) { i j J }
 //                                                                           { k l 2 }
 //
 Operator QdotQ_Op(ModelSpace& modelspace)
 {
    
//   // temporarily store <i||Q||j> in the one body part.
//   Operator QdotQ_op = ElectricMultipoleOp(modelspace,2);
   Operator QdotQ_op(modelspace,0,0,0,2);
   auto Qmat = ElectricMultipoleOp(modelspace,2).OneBody;
//   double b2 =  HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(); // b^2 = hbar/mw 
   int nchan = modelspace.GetNumberTwoBodyChannels();
//
//   // temporarily store <i||Q||j> in the one body part.
//   for (size_t i=0;i<modelspace.GetNumberOrbits();i++)
//   {
//     for (size_t j=0;j<=i;j++)
//     {
//       Orbit & oi = modelspace.GetOrbit(i);
//       Orbit & oj = modelspace.GetOrbit(j);
//       double ji = oi.j2*0.5;
//       double jj = oj.j2*0.5;
//       double r2_ij = RadialIntegral(oi.n,oi.l,oj.n,oj.l,2) * b2 ;
//       double Qij = modelspace.phase(jj+2-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*2+1)/4./3.1415926) * AngMom::ThreeJ(ji,jj, 2, 0.5, -0.5,0) * r2_ij;
//       QdotQ_op.OneBody(i,j) = Qij;
//       QdotQ_op.OneBody(j,i) = modelspace.phase( ji-jj ) * Qij;
//     }
//   }

   for (size_t i=0;i<modelspace.GetNumberOrbits();i++)
   {
     Orbit & oi = modelspace.GetOrbit(i);
     for (auto j : QdotQ_op.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
     {
       Orbit & oj = modelspace.GetOrbit(j);
//       for (auto k : QdotQ_op.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
       double Qij =0;
       for (size_t k=0;k<modelspace.GetNumberOrbits();k++)
       {
         Orbit & ok = modelspace.GetOrbit(k);
//         QdotQ_op.OneBody(i,j) += modelspace.phase( (oj.j2-oj.j2)/2 ) * Qmat(i,k) * Qmat(k,j) / sqrt(oi.j2+1);
         Qij += modelspace.phase( (oj.j2-oj.j2)/2 ) * Qmat(i,k) * Qmat(k,j) / sqrt(oi.j2+1);
         if (i==0)
         {
           std::cout << "i,j,k = " << i << " " << j << " " << k << "Qki = " << Qmat(i,k) << "  Qkj = " << Qmat(k,j) << "    Qij = " << Qij << std::endl;
         }
       }
       QdotQ_op.OneBody(i,j) = Qij;
     }
   }

   std::cout << QdotQ_op.OneBody << std::endl << std::endl; 

   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      int J = tbc.J;
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace.GetOrbit(i);
         Orbit & oj = modelspace.GetOrbit(j);
         double ji = oi.j2*0.5;
         double jj = oj.j2*0.5;


         for (int iket=ibra;iket<nkets;++iket)
         {
            
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace.GetOrbit(k);
            Orbit & ol = modelspace.GetOrbit(l);
            double jk = ok.j2*0.5;
            double jl = ol.j2*0.5 ;

//            double r2_il = RadialIntegral(oi.n,oi.l,ol.n,ol.l,2) * b2 ;
//            double Qil = modelspace.phase(jl+2-0.5) * sqrt( (2*ji+1)*(2*jl+1)*(2*2+1)/4./3.1415926) * AngMom::ThreeJ(ji,jl, 2, 0.5, -0.5,0) * r2_il;
//            double r2_jk = RadialIntegral(oj.n,oj.l,ok.n,ok.l,2) * b2 ;
//            double Qjk = modelspace.phase(jk+2-0.5) * sqrt( (2*jj+1)*(2*jk+1)*(2*2+1)/4./3.1415926) * AngMom::ThreeJ(jj,jk, 2, 0.5, -0.5,0) * r2_jk;

            double Qki = Qmat(k,i);
            double Qli = Qmat(l,i);
            double Qjk = Qmat(j,k);
            double Qjl = Qmat(j,l);
            double Qik = Qmat(i,k);
            double Qil = Qmat(i,l);

            // Formula just taken from Suhonen 8.55, 8.56
            double QdQ = modelspace.phase( ji+jj+J)     * modelspace.GetSixJ(ji,jj,J,jl,jk,2) * Qki * Qjl
                       - modelspace.phase( ji+jj+jk+jl) * modelspace.GetSixJ(ji,jj,J,jk,jl,2) * Qli * Qjk;

//            double QdQ = Qil * Qjk * (2*J+1)/sqrt(5.0) * modelspace.phase( jk-jj ) * modelspace.GetSixJ(ji,jj,J,jk,jl,2.0);
//            double QdQ = 0.5 * Qil * Qjk * (2*J+1)/sqrt(5.0) * modelspace.phase( jk+jj ) * modelspace.GetSixJ(ji,jj,J,jk,jl,2.0)
//                       - 0.5 * Qik * Qjl * (2*J+1)/sqrt(5.0) * modelspace.phase( jl+jj ) * modelspace.GetSixJ(ji,jj,J,jl,jk,2.0);
            if (i==j) QdQ /= sqrt(2.0);
            if (k==l) QdQ /= sqrt(2.0);
            QdotQ_op.TwoBody.SetTBME(ch,ibra,iket,QdQ);
         }
      }
   }
   return QdotQ_op;
 }




 Operator Dagger_Op( ModelSpace& modelspace, index_t Q )
 {
   Operator dag(modelspace);
   dag.SetNumberLegs(3);
   dag.SetQSpaceOrbit(Q);
   dag.OneBody(Q,Q)= 1.0;
   dag.SetNonHermitian();
   std::cout << "Making a dagger operator. I think Q = " << Q << std::endl;
   return dag;
 }



 Operator VCentralCoulomb_Op( ModelSpace& modelspace, int lmax ) // default lmax=99999
 {
//   std::cout << "Making VCentralCoulomb_Op. lmax is " << lmax  << std::endl;
   Operator VCoul(modelspace, 0,0,0,2);
   double oscillator_b = sqrt(HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
   double alpha_FS = 1.0 / 137.035999;
  
// First, the one-body piece <a|1/r|b>
//   int norb = modelspace.GetNumberOrbits();
//   for (int a=0; a<norb; a++)
   for (auto a : modelspace.all_orbits)
   {
     Orbit& oa = modelspace.GetOrbit(a);
     if (oa.tz2>0) continue; // protons only
     if (oa.l>lmax) continue;
     for (int b : VCoul.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
     {
       if (b<a) continue;
       Orbit& ob = modelspace.GetOrbit(b);
       double rad_int =  RadialIntegral_RpowK(oa.n, oa.l, ob.n, ob.l, -1) ;  
       VCoul.OneBody(a,b) = rad_int;
       VCoul.OneBody(b,a) = rad_int;
     }
   }
   VCoul.OneBody *= alpha_FS * HBARC / oscillator_b; // convert from oscillator units to fermi
//   std::cout << "Oscillator b = " << oscillator_b << std::endl;
//   std::cout << "One body part done. it looks like" << std::endl << VCoul.OneBody << std::endl;
   return VCoul;
 }




 Operator VCoulomb_Op( ModelSpace& modelspace, int lmax ) //default lmax=99999
 {
   std::cout << "Making VCoulomb_Op" << std::endl;
   double t_start = omp_get_wtime();
   Operator VCoul(modelspace, 0,0,0,2);
   double oscillator_b = sqrt(HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
   double alpha_FS = 1.0 / 137.035999;
  
//// First, the one-body piece <a|1/r|b>
//   int norb = modelspace.GetNumberOrbits();
//   for (int a=0; a<norb; a++)
//   {
//     Orbit& oa = modelspace.GetOrbit(a);
//     if (oa.tz2>0) continue; // protons only
//     for (int b : VCoul.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
//     {
//       Orbit& ob = modelspace.GetOrbit(b);
//       double rad_int =  RadialIntegral_RpowK(oa.n, oa.l, ob.n, ob.l, -1) ;  
//       VCoul.OneBody(a,b) = rad_int;
//       VCoul.OneBody(b,a) = rad_int;
//     }
//   }
//   VCoul.OneBody /= oscillator_b; // convert from oscillator units to fermi
   int nmax = modelspace.GetEmax();
//   int lrelmax = 2*modelspace.GetLmax();
   int lrelmax = modelspace.GetEmax();
//   std::vector<double> RadialIntegrals( (nmax+1)*(nmax+1)*(lrelmax+1) );
   std::unordered_map<size_t,double> RadialIntegrals;
// #pragma omp parallel for schedule(dynamic,1) // not sure this is even necessary...
   for (int na=0;na<=nmax;na++)
   {
    for (int nb=0;nb<=na; nb++)
    {
      for (int l=0;l<=lrelmax ; l++)
      {
        size_t hash      = na*(nmax+1)*(lrelmax+1) + nb*(lrelmax+1) + l;   
        size_t flip_hash = nb*(nmax+1)*(lrelmax+1) + na*(lrelmax+1) + l;   
        double rint = RadialIntegral_RpowK(na, l, nb, l, -1);
        RadialIntegrals[hash] = rint;
        RadialIntegrals[flip_hash ] = rint;
      }
    }
   }

   std::cout << "now the big loop... lrelmax = " << lrelmax << "  size of RadInt = " << RadialIntegrals.size()  << std::endl;

// Now the (antisymmetrized) two-body piece <ab| 1/r_rel |cd>
   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   std::cout << "Done Precalculating Moshinsky." << std::endl;
   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;
   #pragma omp parallel for schedule(dynamic,1)  // It would appear that something's not thread-safe in this routine...
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      if (tbc.Tz >= 0) continue; // 2-body coulomb only acts in pp channel
      int J = tbc.J;
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         Orbit & oa = modelspace.GetOrbit(bra.p);
         Orbit & ob = modelspace.GetOrbit(bra.q);
         int na = oa.n;
         int nb = ob.n;
         int la = oa.l;
         int lb = ob.l;
//         std::cout << "ibra = " << ibra << "  a,b   la,lb = " << bra.p << " " << bra.q << "    " << la << " " <<lb << std::endl;
         double ja = oa.j2*0.5;
         double jb = ob.j2*0.5;
         int fab = 2*na + 2*nb + la + lb;

         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket & ket = tbc.GetKet(iket);

            Orbit & oc = modelspace.GetOrbit(ket.p);
            Orbit & od = modelspace.GetOrbit(ket.q);
         
            int nc = oc.n;
            int nd = od.n;
         
            int lc = oc.l;
            int ld = od.l;
            if (la>lmax or lb>lmax or lc>lmax or ld>lmax) continue;
         
            double jc = oc.j2*0.5;
            double jd = od.j2*0.5;
            int fcd = 2*nc + 2*nd + lc + ld;
            if (std::abs(fab-fcd)%2 >0) continue; //  parity conservation


            double rinv=0;
            // Transform to LS coupling using 9j coefficients
            for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
            {
              for (int Sab=0; Sab<=1; ++Sab)
              {
                if ( std::abs(Lab-Sab)>J or Lab+Sab<J) continue;
         
                double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
                if (std::abs(njab) <1e-7) continue;
                int Scd = Sab;
                int Lcd = Lab;
                double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
                if (std::abs(njcd) <1e-7) continue;
                // Next, transform to rel / com coordinates with Moshinsky tranformation
                for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
                {
                  for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
                  {
                    int Lam_cd = Lam_ab; // 1/r conserves lam and Lam, ie relative and com orbital angular momentum
                    for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
                    {
                       if (Lab<std::abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;

                       // factor to account for antisymmetrization
                       int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
                       if ( asymm_factor ==0 ) continue;
         
                       int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
                       int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation
         
                       double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                       if (std::abs(mosh_ab)<1e-8) continue;
//                       std::cout << "I asked for moshab: <" << N_ab << " " << Lam_ab << " " << n_ab << " " << lam_ab << " | " << na << " " << la << " " << nb << " " << lb << " >_" << Lab << std::endl;
         
                       int N_cd = N_ab;
                       int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                       if (n_cd < 0) continue;
//                       if  (n_ab != n_cd and N_ab != N_cd) continue;
         
                       double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                       if (std::abs(mosh_cd)<1e-8) continue;

                       double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;

                       size_t hash      = n_ab*(nmax+1)*(lrelmax+1) + n_cd*(lrelmax+1) + lam_ab;   
                       if ( RadialIntegrals.find(hash) == RadialIntegrals.end() )
                       {
                         std::cout << "AAAHHH!!!  trying to access radial integral for " << n_ab << " " << n_cd << " " << lam_ab << "    and it's not there!!!!" << std::endl;
                       }
                       double rad_int =  RadialIntegrals.at(hash) ;  
//                       double rad_int_lookup =  RadialIntegrals[hash] ;  
//                       double rad_int =  RadialIntegral_RpowK(n_ab, lam_ab, n_cd, lam_cd, -1) ;  
//                       if (std::abs( rad_int_lookup - rad_int)>1e-6)
//                       {
//                         std::cout << "discrepancy!!!  " << n_ab << " " << n_cd << " " << lam_ab << "    |   " << rad_int << "   " << rad_int_lookup << "   (maxvals) " << nmax << " " << lrelmax << std::endl;
//                       }
//                       size_t a_eff = modelspace.GetOrbitIndex(n_ab, lam_ab, 2*lam_ab+1, -1 );
//                       size_t b_eff = modelspace.GetOrbitIndex(n_cd, lam_cd, 2*lam_cd+1, -1 );
//                       std::cout << "a_eff, b_eff = " << a_eff << " " << b_eff
//                                 << "   " << n_ab << " " << lam_ab << " " << 2*lam_ab+1 << "    " << fab << ","
//                                 << "   " << n_cd << " " << lam_cd << " " << 2*lam_cd+1 << "    " << fcd
//                                 << std::endl;
//                       rinv += prefactor * VCoul.OneBody(a_eff, b_eff); 
                       rinv += prefactor * rad_int; 
    
                    } // lam_ab
                  } // Lam_ab
                } // N_ab
         
              } // Sab
            } // Lab

            // In Moshinsky's convention, r_rel = (r1-r2)/sqrt(2).  We want 1/|r1-r2| = 1/sqrt(2) * 1/r_rel
            rinv *=  1 / sqrt(2*(1.0+bra.delta_pq())*(1.0+ket.delta_pq())); // normalize and account for sqrt(2) Moshinsky convention
//            std::cout << "setting " << ch << " " << ibra << " " << iket << "  " << rinv << std::endl;
            VCoul.TwoBody.SetTBME(ch,ibra,iket,rinv);
            VCoul.TwoBody.SetTBME(ch,iket,ibra,rinv);
                         
         }
      }
   }
//   VCoul.OneBody.zeros(); // We don't want the 1-body Coulomb (unless we're doing an atom...)
   VCoul *= alpha_FS * HBARC / oscillator_b;  // convert to MeV.  V = e^2/r = alpha*hc / r

   std::cout << "All done with VCoul." << std::endl;
   VCoul.profiler.timer["VCoulomb_Op"] += omp_get_wtime() - t_start;
   return VCoul ;

 }



 namespace atomic_fs
 { // operators related to fine structure
  
   Operator Darwin(ModelSpace& modelspace, int Z )
   {
     double alpha_FS = 1.0 / 137.035999;
     double constants = M_PI * Z * alpha_FS * HBARC*HBARC*HBARC / (2*M_ELECTRON*M_ELECTRON*1e6*1e6) ; // convert to eV. M_PI is 3.1415... not the pion mass
     Operator Hdarwin( modelspace,0,0,0,2);
     for (auto a : modelspace.all_orbits )
     {
       Orbit& oa = modelspace.GetOrbit(a);
       double wf0_a = imsrg_util::HO_Radial_psi(oa.n, oa.l, modelspace.GetHbarOmega(), 0.0);
       if ( oa.l!=0) continue;  // no spin-orbit in s-wave
       for ( auto b : Hdarwin.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
       {
         Orbit& ob = modelspace.GetOrbit(b);
         double wf0_b = imsrg_util::HO_Radial_psi(ob.n, ob.l, modelspace.GetHbarOmega(), 0.0);
         Hdarwin.OneBody(a,b) = constants * wf0_a * wf0_b;
         Hdarwin.OneBody(b,a) = Hdarwin.OneBody(a,b);
       }
     }
     return Hdarwin;
   }

   Operator RelativisticT(ModelSpace& modelspace )
   {
     Operator Hrel = imsrg_util::KineticEnergy_RelativisticCorr(modelspace) * 1e6*M_NUCLEON/(M_ELECTRON); // change to electron mass, as use eV rather than MeV.
     return Hrel;
   }

   Operator SpinOrbit( ModelSpace& modelspace, int Z )
   {
     Operator Hso( modelspace, 0,0,0,2);
     double oscillator_b = sqrt(HBARC*HBARC/(1e6*M_ELECTRON)/modelspace.GetHbarOmega()); // convert electron mass to eV
     double oscillator_b3 = pow(oscillator_b,3);
     double alpha_FS = 1.0 / 137.035999;
     double gspin = 2.002319; // electron spin g factor
     double constants = Z*alpha_FS * HBARC*HBARC * gspin / (M_ELECTRON*M_ELECTRON*1e6*1e6) / 32;  // it's 1/8, but we use 4 * LdotS, so 1/32.
     for (auto a : modelspace.all_orbits )
     {
       Orbit& oa = modelspace.GetOrbit(a);
       if ( oa.l==0) continue;  // no spin-orbit in s-wave
       int four_ldots = oa.j2*(oa.j2+2) - 4*oa.l*(oa.l+1) -3 ;
       for ( auto b : Hso.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
       {
         Orbit& ob = modelspace.GetOrbit(b);
         double r3inv = imsrg_util::RadialIntegral_RpowK(oa.n, oa.l, ob.n, ob.l, -3);
         Hso.OneBody(a,b) = constants * four_ldots * r3inv  / oscillator_b3;
         Hso.OneBody(b,a) = Hso.OneBody(a,b);
       }
     }
     return Hso;
   }


 }// namespace atomic_fs


 namespace atomic_hfs
 { // operators related to hyperfine structure

   Operator hQ(ModelSpace& modelspace )
   {
     Operator Hq( modelspace,0,0,0,2);
     return Hq;
   }

   // The magnetic dipole term consists of three contributions:
   // The orbit term, the tensor term, and the contact term
   // Hd = -0.5*alpha(hbarc)^3/(m_ec^2 m_pc^2) g_nuc I * [ r^-3 L  +1/2 g_s r^-3 ( 3(\vec{s}*\hat{r})\hat{r} - \vec{s} ) + 4pi/3 gs delta(r) \vec{s} )
   // we rewrite the tensor bit as
   //                                 3(s*r)r-s = -sqrt{2pi}[s^(2) x Y^(2)]^(1) 
   Operator hD(ModelSpace& modelspace )
   {
     Operator Hd( modelspace,1,0,0,2);  // J rank is 1, even parity.

     double oscillator_b = sqrt(HBARC*HBARC/(1e6*M_ELECTRON)/modelspace.GetHbarOmega()); // convert electron mass to eV
     double oscillator_b3 = pow(oscillator_b,3);
     double alpha_FS = 1.0 / 137.035999;
     double gspin = 2.002319; // electron spin g factor
     double constants = - 0.5*alpha_FS *HBARC*HBARC*HBARC/(M_ELECTRON*M_NUCLEON*1e12);  // convert both masses to eV
     for ( auto a : modelspace.all_orbits )
     {
       Orbit& oa = modelspace.GetOrbit(a);
       for (auto b : modelspace.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
       {
         Orbit& ob = modelspace.GetOrbit(b);
         if (oa.l==0 and ob.l==0)
         {
           double wf0_a = imsrg_util::HO_Radial_psi(oa.n, oa.l, modelspace.GetHbarOmega(), 0.0);
           double wf0_b = imsrg_util::HO_Radial_psi(ob.n, ob.l, modelspace.GetHbarOmega(), 0.0);
           // the reduced matrix element of s is <1/2|| s || 1/2> = sqrt(3/2)
           Hd.OneBody(a,b) = 4*M_PI/3 * gspin * sqrt(3./2) * wf0_a * wf0_b ;
           Hd.OneBody(b,a) = Hd.OneBody(a,b);
         }
         else
         {
           double r3inv = imsrg_util::RadialIntegral_RpowK(oa.n, oa.l, ob.n, ob.l, -3) / oscillator_b3;
           double L = oa.l!=ob.l ? 0 : sqrt((oa.j2+1.0)/(oa.j2*(oa.j2+2))) * (oa.j2*(oa.j2+2.)/4 +oa.l*(oa.l+1) -3./4);
           double T = modelspace.phase(oa.l) * 3*sqrt(5)*sqrt((oa.j2+1)*(ob.j2+1)*(2*oa.l+1)*(2*ob.l+1)) * AngMom::ThreeJ(oa.l,2,ob.l,0,0,0) * AngMom::NineJ(oa.l,0.5,0.5*oa.j2, ob.l,0.5,0.5*ob.j2, 2,1,1);
           Hd.OneBody(a,b) = constants * r3inv *( L - gspin/2 * T );
           Hd.OneBody(b,a) = Hd.OneBody(a,b);
         }
       }
     }

     return Hd;
   }



   // Kinetic energy T = T_el + T_nuc =  1/2m sum_i (p_i)^2 + 1/2M_nuc ( sum_i p_i )^2  =  (1/2m + 1/2M_n) sum_i (p_i)^2 + 1/2M_nuc sum_ij (p_i * p_j)
   // The first correction, the 1/2M_n one-body part is responsible for what is called the "Normal Mass Shift", while the second correction
   // which goes like p_i * p_j (* means a vector dot product here), is responsible for the "Specific Mass Shift".
   Operator NormalMassShift( ModelSpace& modelspace, int A )
   {
     Operator Hnms = (M_ELECTRON/A/M_NUCLEON) * imsrg_util::KineticEnergy_Op( modelspace ) ;  // kinetic energy is in units of hw, so no change needed
     if (A!=modelspace.GetTargetMass()) Hnms *= (modelspace.GetTargetMass()/double(A));
     return Hnms;
   }

   
   Operator SpecificMassShift( ModelSpace& modelspace, int A )
   {
     Operator Hsms = imsrg_util::TCM_Op( modelspace ) ;  // TCM_Op returns a 1-body piece, plus the 1-body part pi*pj/mA. We don't want the 1-body part.
     Hsms.OneBody.zeros(); // The specific shift is just the two-body part.
     if (A!=modelspace.GetTargetMass()) Hsms *= (modelspace.GetTargetMass()/double(A));
     return Hsms;
   }

   // Maybe we want it all in one operator
   Operator CombinedMassShift( ModelSpace& modelspace, int A )
   {
     Operator Hcms = imsrg_util::TCM_Op( modelspace ) ; 
     if (A!=modelspace.GetTargetMass()) Hcms *= (modelspace.GetTargetMass()/double(A));
     return Hcms;
   }


 }// namespace atomic_hfs


/// Get the first-order perturbative correction to a one-body operator
/// The function returns a one-body operator, i.e. the two-body part is not computed.
/// Formula is 
/// \f[
///  \langle p \| O^{\lambda} \| q \rangle = \sum_{ia} (1+P_{ia}) (n_i \bar{n}_a) \sum_{J}(2J+1)
/// \begin{Bmatrix}
/// j_i & j_a & \lambda \\
/// j_p & j_q & J
/// \end{Bmatrix}  \langle i \| O^{\lambda} \| a \rangle \frac{ \tilde{\Gamma}^{J}_{paij}}{\Delta_{paiq}}
/// /f]
//
//So what you get should look like this:
//      p|                p|______
//       |     ^~~~~X      |     ^
//       |   a( )i     +   |   a( )i
//       |_____v           |     v~~~~X
//      q|                q|     
//
//
 Operator FirstOrderCorr_1b( const Operator& OpIn, const Operator& H )
 {

//   Operator OpIn = OpInx;
   Operator OpOut = 0. * OpIn;
   int Lambda = OpOut.GetJRank();
   size_t norb = OpIn.modelspace->GetNumberOrbits();

//   arma::uvec holevec ( OpIn.modelspace->holes );
//   arma::uvec particlevec ( OpIn.modelspace->particles );
//   OpIn.OneBody.submat(holevec,holevec).zeros();
//   OpIn.OneBody.submat(particlevec,particlevec).zeros();

   for ( auto p : OpIn.modelspace->valence )
   {
     Orbit& op = OpIn.modelspace->GetOrbit(p);
     double jp = 0.5 * op.j2;
     for ( auto q : OpIn.OneBodyChannels.at({op.l, op.j2, op.tz2}) )
     {
       if (  std::find(  OpIn.modelspace->valence.begin(), OpIn.modelspace->valence.end(), q )  == OpIn.modelspace->valence.end()  ) continue;
       Orbit& oq = OpIn.modelspace->GetOrbit(q);
       double jq = 0.5 * oq.j2;
       double Opq = 0;
//       double Opq_alt = 0;
       for ( size_t i=0; i<norb; i++)
       {
         Orbit& oi = OpIn.modelspace->GetOrbit(i);
         double ji = 0.5 * oi.j2;
         double ni = oi.occ;
         for ( auto a : OpIn.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
         {
           Orbit& oa = OpIn.modelspace->GetOrbit(a);
           double ja = 0.5 * oa.j2;
           double na = oa.occ;
           double nanifactor = ni-na;
           double Oia = OpIn.OneBody(i,a);
           if (std::abs(nanifactor)<1e-6) continue;
           if (std::abs(Oia)<1e-6) continue;
//           std::cout << "Evaluating PT1.  proceeding with term a,i = " << a << " " << i << " : " << ja << " " << ji << "  Oia = " << Oia <<  std::endl;
           int Jmin = std::max(  std::max( std::abs(jp-ja), std::abs(ji-jq) ) ,   std::max( std::abs(jp-ji), std::abs(ja-jq) )  );
           int Jmax = std::max(  std::min( jp+ja , ji+jq ) ,   std::min( jp+ji , ja+jq )  );
           for (int J=Jmin; J<=Jmax; J++)
           {
             double sixj = OpIn.modelspace->GetSixJ( ji, ja, Lambda, jp, jq, J );
             double Gamma_paiq = H.TwoBody.GetTBME_J(J,p,a,i,q);
             double Delta_paiq = H.OneBody(p,p) + H.OneBody(a,a) - H.OneBody(i,i) - H.OneBody(q,q);
             Opq += nanifactor * (2*J+1) * sixj * Oia * Gamma_paiq / Delta_paiq;
           }           
         }
       }
       OpOut.OneBody(p,q) = Opq;
//       std::cout << "Set Opq  with p,q = " << p << " " << q << "   to " << Opq << std::endl;
     }
   }
   return OpOut;
 }





 Operator RPA_resummed_1b( const Operator& OpIn, const Operator& H, std::string mode )
 {

   // construct hp and ph kets,  as well as Oph and Ohp
   int Lambda = OpIn.GetJRank();
//   size_t ch_CC = OpIn.modelspace->GetTwoBodyChannelIndex( Lambda, OpIn.GetParity(), OpIn.GetTRank() );
//   std::cout << "going for a channel with J,p,Tz = " << Lambda << " " << OpIn.GetParity() << " " << OpIn.GetTRank() << std::endl;
//   TwoBodyChannel_CC& tbc_cc_ph = OpIn.modelspace->GetTwoBodyChannel_CC( ch_CC );
//   auto ketindex_ph = tbc_cc_ph.GetKetIndex_ph();
//   size_t nkets = ketindex_ph.size();

   std::vector<std::pair<size_t,size_t>> ph_kets;
   std::vector<std::pair<size_t,size_t>> hp_kets;

   // maybe we try the dumb way
   for ( auto h : OpIn.modelspace->holes )
   {
//     Orbit& oh = OpIn.modelspace->GetOrbit(h);
//     for (auto p : OpIn.OneBodyChannels.at({oh.l,oh.j2,oh.tz2}) )
     for (auto p : OpIn.modelspace->particles )
     {
       Orbit& op = OpIn.modelspace->GetOrbit(p);
       if (op.occ>0.01) continue;
//       Oph( ph_kets.size() ) = OpIn.OneBody(p,h);
//       Ohp( hp_kets.size() ) = OpIn.OneBody(h,p);
       ph_kets.push_back( std::make_pair(p,h) );
       hp_kets.push_back( std::make_pair(h,p) );
     }
   }
   size_t nkets = ph_kets.size();
   arma::vec Oph( nkets, arma::fill::zeros );
   arma::vec Ohp( nkets, arma::fill::zeros );
   for (size_t i=0; i<nkets; i++)
   {
     auto p = ph_kets[i].first;
     auto h = ph_kets[i].second;
     Oph( i ) = OpIn.OneBody(p,h);
     Ohp( i ) = OpIn.OneBody(h,p);
   }

   
   // Next, construct the Mphph etc matrices
   arma::mat Mphph = GetPH_transformed_Gamma( ph_kets, ph_kets, H, Lambda );
   arma::mat Mphhp = GetPH_transformed_Gamma( ph_kets, hp_kets, H, Lambda );
   arma::mat Mhpph = GetPH_transformed_Gamma( hp_kets, ph_kets, H, Lambda );
   arma::mat Mhphp = GetPH_transformed_Gamma( hp_kets, hp_kets, H, Lambda );

   // to get TDA, we just get rid of the off-diagonal blocks
   if (mode=="TDA")
   {
     Mphhp *=0;
     Mphhp *=0;
   }

   // for the full M matrix  M = [ Mphph   Mphhp ]
   //                            [ Mhpph   Mhphp ]
   arma::mat M = arma::join_vert( arma::join_horiz( Mphph, Mphhp) ,
                                  arma::join_horiz( Mhpph, Mhphp) );

   // make the base case denominator eq-ep for ket |pq>
   arma::mat Delta(arma::size(M), arma::fill::ones );
   for ( size_t i=0; i<nkets; i++)
   {
     size_t p = ph_kets[i].first;
     size_t h = ph_kets[i].second;
     double del_ph = H.OneBody(p,p) - H.OneBody(h,h);
     Delta.col(i) *= -del_ph;
     Delta.col(i+nkets) *= del_ph;
   }

   // combine Oph and Ohp into a single column vector
   arma::vec Ovec = arma::join_vert(  Oph, Ohp );

  // At this point, we can't get around the fact that the denominators depend on the initial and final state of the
  // entire operator, as in < p | Oeff | q >, and we need to add ep-eq to the denominator. So for every p,q pair
  // we need to do the matrix inversion and evaluation in 1st order PT.

   Operator OpOut = OpIn;
   OpOut.OneBody.zeros();

   for ( auto v1 : OpIn.modelspace->valence )
   {
     for ( auto v2 : OpIn.modelspace->valence )
     {

      arma::mat del12 = arma::ones(arma::size(Delta)) * (H.OneBody(v1,v1) - H.OneBody(v2,v2)) ;
      // do the fancy resummation of the series by matrix inversion
      // Minv  =  (I-M)^-1  where I is the identity matrix. The slash here means element-wise division.
       arma::mat Minv = arma::inv(  arma::eye(arma::size(M))  - M/( Delta+del12) );
       arma::vec ORPA = Minv * Ovec;


       // now we need to unpack ORPA into a matrix again.
       Operator OpRPA = OpIn;
       OpRPA.OneBody.zeros();
       for ( size_t i=0; i<nkets; i++ )
       {
         size_t p = ph_kets[i].first;
         size_t h = ph_kets[i].second;
         OpRPA.OneBody(p,h) = ORPA(i);
         OpRPA.OneBody(h,p) = ORPA(i+nkets);
       }
      // evaluate this in first order perturbation theory, and take the valence part that we're interested in
       double op12 = FirstOrderCorr_1b( OpRPA, H).OneBody(v1,v2);
       OpOut.OneBody(v1,v2) += op12;

      }
    }
   return OpOut;

 }




// Return Pandya-transformed matrix elements. Input a list of bras <pq`| a list of kets |rs`> and an angular momentum Lambda (as well as the Hamiltonian H
// and obtain  Gamma`_pq`rs`^Lambda = - sum_J  (2J+1)  { jp  jq  Lambda  }  Gamma_psrq^J
//                                                     { jr  js  J       }
//
// arma::mat GetPH_transformed_Gamma( std::vector<size_t> bras, std::vector<size_t> kets, Operator& H, int Lambda )
 arma::mat GetPH_transformed_Gamma( std::vector<std::pair<size_t,size_t>>& bras, std::vector<std::pair<size_t,size_t>>& kets, const Operator& H, int Lambda )
 {
   arma::mat Mat_pandya( bras.size(), kets.size(), arma::fill::zeros);
   for ( size_t ibra=0; ibra<bras.size(); ibra++)
   {
     size_t p = bras[ibra].first;
     size_t q = bras[ibra].second;
     Orbit& op = H.modelspace->GetOrbit(p);
     Orbit& oq = H.modelspace->GetOrbit(q);
     double jp = 0.5 * op.j2;
     double jq = 0.5 * oq.j2;
     for ( size_t iket=0; iket<kets.size(); iket++)
     {
       size_t r = kets[iket].first;
       size_t s = kets[iket].second;
       Orbit& oR = H.modelspace->GetOrbit(r);
       Orbit& os = H.modelspace->GetOrbit(s);
       double jr = 0.5 * oR.j2;
       double js = 0.5 * os.j2;

       double matel = 0;
       int Jmin = std::max( std::abs(jp-js) , std::abs(jq-jr) );
       int Jmax = std::min( jp+js ,  jq+jr );
       for (int J=Jmin; J<=Jmax; J++)
       {
         double sixj = H.modelspace->GetSixJ( jp, jq, Lambda, jr, js, J );
         matel -= (2*J+1) * sixj * H.TwoBody.GetTBME_J(J,p,s,r,q);
       }
       Mat_pandya( ibra,iket ) = matel;
     }
   }
   return Mat_pandya;
 }





 // Evaluate <bra | r1*r2 | ket>, omitting the factor (hbar * omega) /(m * omega^2)
/// Returns the normalized, anti-symmetrized, J-coupled, two-body matrix element of \f$ \frac{m\omega^2}{\hbar \omega} \vec{r}_1\cdot\vec{r}_2 \f$.
/// Calculational details are similar to Calculate_p1p2().
 double Calculate_r1xp2(ModelSpace& modelspace, Ket & bra, Ket & ket, int Jab, int Jcd)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
   int nb = ob.n;
   int nc = oc.n;
   int nd = od.n;

   int la = oa.l;
   int lb = ob.l;
   int lc = oc.l;
   int ld = od.l;

   double ja = oa.j2/2.0;
   double jb = ob.j2/2.0;
   double jc = oc.j2/2.0;
   double jd = od.j2/2.0;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   if (std::abs(fab-fcd)%2 >0) return 0; // p1*p2 only connects kets with delta N = 0,1
   if (std::abs(fab-fcd)>2) return 0; // p1*p2 only connects kets with delta N = 0,1


   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double r1xp2=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( std::abs(Lab-Sab)>Jab or Lab+Sab<Jab) continue;

       double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,Jab);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,Jcd);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // Lcm and Lrel conserve lam and Lam, 
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<std::abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
              // factor to account for antisymmetrization

              int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // Lcm and Lrel operators conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (std::abs(mosh_ab)<1e-8) continue;

//              for (int N_cd=std::max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
               int N_cd = N_ab;
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd ) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (std::abs(mosh_cd)<1e-8) continue;

//                double lcm = 0;
//                double lrel = 0;
                double sj_rel = modelspace.GetSixJ( lam_ab, lam_cd, 1, Lcd, Lab, Lam_ab );
                double sj_cm = modelspace.GetSixJ( Lam_ab, Lam_cd, 1, Lcd, Lab, lam_ab );
                double lrel = modelspace.phase( Lab + Lcd ) * ( (2*Lab+1)*(2*Lcd+1) ) * sj_rel * sj_rel  * sqrt(lam_ab*(lam_ab+1)*(2*lam_ab+1));
                double lcm = ( (2*Lab+1)*(2*Lcd+1) ) * sj_cm * sj_cm  * sqrt(Lam_ab*(Lam_ab+1)*(2*Lam_ab+1));
                r1xp2 += (lcm-lrel) * 0.5;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize.
   r1xp2 *= 1.0 / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return r1xp2 ;

 }


//////////// M0v functions written by Charlie Payne //////////////////////////

/// This is the M^{0\nu} TBME from Equation (1) of [PRC 87, 064315 (2013)]
/// it was coded up by me, ie) Charlie Payne (CP)
/// from my thesis, I employed Equations: 
  Operator M0nu_TBME_Op(ModelSpace& modelspace, int Nquad, std::string src)
  {
    // VVV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~VVV
    // VVV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~VVV
    // adjustable parameters, lines BELOW
    //double reltol = 2*pow(10,-4); // relative tolerance for GLQ integration convergence
    const double mpro = 938.27231; // the proton mass [MeV] for the g-factors
    const double mpion = 139.57; // the pion mass [MeV] for the g-factors
    const double magmom = 3.706; // the difference between the (anomolous?) magnetic moment of a proton and neutron (units of \mu_N)
    const double g0V = 1.0; // the vector g-factor at zero momentum
    const double g0A = 1.27; // the axial-vector g-factor at zero momentum
    const double cutoffV = 850.0; // the vector finite-size parameter [MeV]
    const double cutoffA = 1086.0; // the axial-vector finite-size parameter [MeV]
    const double Ebar = 5.0; // the ref-closure energy [MeV], should be roughly independent of this...
    const double Ei = 0.0; // the initial energy [MeV] in the closure energy (4.26798)
    const double Ef = 0.0; // " final " [MeV] " " " "
    //const double mnusq = pow(1.2,-14); // the (max) electron-neutrino mass squared [MeV^2] (this is primarily for testing)
    // adjustable parameters, lines ABOVE
    // ^^^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^
    // ^^^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^
    Operator M0nu_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    M0nu_TBME.SetHermitian(); // it should be Hermitian...
    const double Rnuc = (1.2/HBARC)*pow(modelspace.GetTargetMass(),1.0/3.0); // the nuclear radius [MeV^-1]
    const double prefact = 2*(2*Rnuc)/(PI*g0A*g0A); // factor in-front of M0nu TBME [MeV^-1], factor of 2 from sum_{ab} = 2sum_{a<b} to match with JE
    const double Ediff = Ebar - ((Ei + Ef)/2.0); // the closure energy [MeV]
    // set the SRC parameters a,b,c in the Jastrow-type correlation function
    double aa,bb; // [MeV^2]
    double cc; // [Unitless]
    std::string argstr = "Argonne";
    std::string cdbstr = "CD-Bonn";
    std::string masstr = "Miller-Spencer";
    if (src == argstr)
    {
      aa = 1.59*HBARC*HBARC;
      bb = 1.45*HBARC*HBARC;
      cc = 0.92;
    }
    else if (src == cdbstr)
    {
      aa = 1.52*HBARC*HBARC;
      bb = 1.88*HBARC*HBARC;
      cc = 0.46;
    }
    else if (src == masstr)
    {
      aa = 1.1*HBARC*HBARC;
      bb = 0.68*HBARC*HBARC;
      cc = 1.0;
    }
    else // src == "none" or otherwise
    {
      aa = 0;  bb = 0;  cc = 0; // set them to zero to avoid variable warnings
    }
    // set the GLQ roots and weights
    if (Nquad != 187)
    {
      std::cout<<"ERROR 187: please set Nquad to 187 by running the operator as M0nu_TBME_187_none"<<std::endl;
      std::cout<<"which I do to stay consistent with previous naming conventions."<<std::endl;
      exit(1);
    }
    double nodes[187][2] = { {0.00771093190434205, 0.01978880917191989},
                               {0.04062903529180086, 0.04606592889118855},
                               {0.09985365551375681, 0.07238518879149865},
                               {0.1854021036766674, 0.09871291032827842},
                               {0.2972819985528331, 0.1250482436418515},
                               {0.4355015848139718, 0.151392567186406},
                               {0.6000707059043154, 0.1777476120597593},
                               {0.7910009885390217, 0.2041151966450439},
                               {1.008305895218917, 0.2304971697691649},
                               {1.252000745415518, 0.2568953952378961},
                               {1.522102727728737, 0.2833117470989823},
                               {1.818630909285556, 0.3097481082650387},
                               {2.141606244467358, 0.3362063703580633},
                               {2.491051583758736, 0.3626884340417398},
                               {2.866991683053023, 0.389196209573469},
                               {3.269453213570937, 0.4157316174666621},
                               {3.698464772473222, 0.4422965892116879},
                               {4.154056894213615, 0.4688930680324417},
                               {4.636262062662657, 0.4955230096710074},
                               {5.145114724024852, 0.5221883831886691},
                               {5.680651300567963, 0.5488911717891486},
                               {6.242910205181583, 0.5756333736530594},
                               {6.831931856781384, 0.6024170027954671},
                               {7.447758696575838, 0.6292440899340502},
                               {8.09043520521219, 0.6561166833776932},
                               {8.760007920819465, 0.6830368499335451},
                               {9.456525457966828, 0.7100066758308241},
                               {10.18003852755672, 0.7370282676660339},
                               {10.93059995767296, 0.7641037533676859},
                               {11.70826471540523, 0.7912352831858124},
                               {12.5130899296724, 0.8184250306993499},
                               {13.34513491506818, 0.8456751938532073},
                               {14.20446119675403, 0.8729879960168208},
                               {15.09113253642537, 0.9003656870720292},
                               {16.00521495937847, 0.9278105445271643},
                               {16.94677678270679, 0.9553248746614917},
                               {17.91588864465706, 0.9829110137007926},
                               {18.91262353517667, 1.01057132902482},
                               {19.93705682768588, 1.038308220409159},
                               {20.98926631210941, 1.066124121303237},
                               {22.0693322292045, 1.094021500145207},
                               {23.17733730622346, 1.122002861716167},
                               {24.31336679395113, 1.15007074853609},
                               {25.47750850515957, 1.17822774230309},
                               {26.6698528545242, 1.20647646537708},
                               {27.89049290004778, 1.234819582312125},
                               {29.13952438604095, 1.263259801437972},
                               {30.4170457877105, 1.291799876493899},
                               {31.72315835740869, 1.320442608317681},
                               {33.05796617259991, 1.349190846591676},
                               {34.42157618560343, 1.378047491649795},
                               {35.81409827517415, 1.407015496346448},
                               {37.23564529998588, 1.436097867993353},
                               {38.68633315408532, 1.465297670364799},
                               {40.16628082438788, 1.494618025776432},
                               {41.67561045029016, 1.524062117241149},
                               {43.21444738547788, 1.553633190703724},
                               {44.78292026201128, 1.583334557362363},
                               {46.38116105677503, 1.613169596076766},
                               {48.00930516038326, 1.643141755871362},
                               {49.66749144863559, 1.67325455853451},
                               {51.35586235662428, 1.703511601322361},
                               {53.07456395559839, 1.733916559770406},
                               {54.82374603269579, 1.764473190618135},
                               {56.60356217365987, 1.795185334854757},
                               {58.41416984866382, 1.826056920890234},
                               {60.25573050137191, 1.85709196785874},
                               {62.12840964137386, 1.888294589062953},
                               {64.03237694013582, 1.919668995564201},
                               {65.96780633061908, 1.951219499929162},
                               {67.93487611072592, 1.98295052013875},
                               {69.93376905074039, 2.014866583670929},
                               {71.9646725049416, 2.046972331764233},
                               {74.02777852757615, 2.079272523874329},
                               {76.12328399338745, 2.111772042333272},
                               {78.25139072291041, 2.144475897222693},
                               {80.41230561275184, 2.177389231474337},
                               {82.6062407710895, 2.21051732620899},
                               {84.83341365863637, 2.243865606330362},
                               {87.09404723533066, 2.277439646385225},
                               {89.38837011302772, 2.311245176709602},
                               {91.71661671448609, 2.345288089874654},
                               {94.07902743895778, 2.379574447450897},
                               {96.47584883471066, 2.41411048711277},
                               {98.90733377883214, 2.448902630099688},
                               {101.3737416646834, 2.483957489059988},
                               {103.8753385973973, 2.519281876298236},
                               {106.4123975978368, 2.554882812454203},
                               {108.9851988154589, 2.590767535637621},
                               {111.5940297505536, 2.626943511051393},
                               {114.2391854863636, 2.663418441132556},
                               {116.9209689316168, 2.700200276245322},
                               {119.6396910740445, 2.73729722596323},
                               {122.3956712454915, 2.774717770979608},
                               {125.189237399269, 2.812470675688552},
                               {128.0207264004418, 2.850565001482234},
                               {130.8904843297925, 2.889010120814019},
                               {133.7988668022529, 2.927815732081122},
                               {136.7462393006521, 2.966991875383191},
                               {139.7329775256874, 3.006548949222304},
                               {142.759467763093, 3.046497728208849},
                               {145.8261072690484, 3.08684938184865},
                               {148.9333046749477, 3.127615494490461},
                               {152.0814804127326, 3.168808086522331},
                               {155.2710671620838, 3.210439636907599},
                               {158.5025103208623, 3.25252310716648},
                               {161.7762685003008, 3.29507196691414},
                               {165.092814046561, 3.338100221076514},
                               {168.4526335904016, 3.381622438918201},
                               {171.8562286268402, 3.425653785027162},
                               {175.3041161268443, 3.470210052418758},
                               {178.7968291832567, 3.515307697927283},
                               {182.3349176933386, 3.560963880085994},
                               {185.9189490805198, 3.607196499697949},
                               {189.5495090581641, 3.654024243333757},
                               {193.2272024383993, 3.701466630010565},
                               {196.9526539893333, 3.749544061332636},
                               {200.7265093442698, 3.79827787540431},
                               {204.5494359668643, 3.847690404859263},
                               {208.4221241765241, 3.897805039384193},
                               {212.3452882387537, 3.948646293160198},
                               {216.3196675255918, 4.000239877687697},
                               {220.3460277517803, 4.052612780516117},
                               {224.4251622928536, 4.105793350456135},
                               {228.5578935919482, 4.159811389922316},
                               {232.7450746628161, 4.214698255126153},
                               {236.9875906972909, 4.27048696492881},
                               {241.2863607863087, 4.327212319259251},
                               {245.6423397645486, 4.384911028115464},
                               {250.0565201898382, 4.443621852294881},
                               {254.5299344696856, 4.503385757145375},
                               {259.0636571486737, 4.564246080797834},
                               {263.6588073720075, 4.62624871853364},
                               {268.3165515422674, 4.689442325160785},
                               {273.0381061884222, 4.75387853753741},
                               {277.8247410684369, 4.819612219670748},
                               {282.6777825294121, 4.88670173317499},
                               {287.5986171521681, 4.95520923627262},
                               {292.5886957106057, 5.025201014998875},
                               {297.6495374801001, 5.096747850818959},
                               {302.7827349337207, 5.169925429530114},
                               {307.9899588703066, 5.24481479708569},
                               {313.2729640245108, 5.321502868882359},
                               {318.6335952160018, 5.400083000160843},
                               {324.0737941032713, 5.480655626428459},
                               {329.5956066171773, 5.563328984392797},
                               {335.2011911607244, 5.648219925738174},
                               {340.8928276750128, 5.735454838322556},
                               {346.67292768718, 5.825170692111637},
                               {352.5440454750781, 5.917516230461636},
                               {358.5088905060044, 6.012653331443834},
                               {364.5703413339042, 6.11075856888457},
                               {370.7314611721156, 6.212025008974579},
                               {376.9955153982821, 6.3166642859827},
                               {383.3659912962149, 6.424909010223284},
                               {389.8466203984501, 6.537015573515938},
                               {396.4414038658405, 6.653267432708673},
                               {403.1546414304768, 6.773978971401204},
                               {409.9909645404044, 6.899500065148414},
                               {416.9553744854671, 7.030221508018704},
                               {424.0532864618092, 7.166581500960829},
                               {431.2905807597857, 7.309073458588748},
                               {438.6736625522015, 7.458255465736875},
                               {446.2095321389017, 7.614761815679939},
                               {453.9058680004354, 7.779317198652938},
                               {461.7711256711422, 7.952754297479251},
                               {469.8146563226232, 8.136035809504971},
                               {478.0468501423313, 8.330282285006321},
                               {486.4793112320975, 8.536807704774064},
                               {495.1250730378306, 8.75716549697966},
                               {503.998866560597, 8.993208849137494},
                               {513.1174582695169, 9.247170924190588},
                               {522.5000815040654, 9.521773309539423},
                               {532.1689954739326, 9.820375353805916},
                               {542.1502218601211, 10.14718412234403},
                               {552.4745341908956, 10.50755665060706},
                               {563.1788163012253, 10.9084470923048},
                               {574.3079759021866, 11.35908954153151},
                               {585.9177225644977, 11.87208048166133},
                               {598.0787485964511, 12.4651733647271},
                               {610.8833037499799, 13.16442107333889},
                               {624.4561163494692, 14.01006813896018},
                               {638.9738566152255, 15.06861784527307},
                               {654.70324721585, 16.46066092751913},
                               {672.0863369696656, 18.43694351949129},
                               {691.9752748484226, 21.65051648474522},
                               {716.5834864211865, 28.8014926391697} }; // nodes[ row ][ 0 = roots, 1 = weights ] 
    // pre-calculate the needed form factors, for efficency
    double hF[Nquad]; // the Fermi form factor
    double hGT[Nquad]; // the Gamow-Teller form factor
    for (int i=0; i<Nquad; i++)
    {
      double q = nodes[i][0]; // set the transfer momentum [MeV] via the impending GLQ
      double qsq = q*q; // q squared
      double mprosq = mpro*mpro; // the proton mass squared
      double mpionsq = mpion*mpion; // the pion mass squared
      double gV = g0V/pow((1.0 + (qsq/(cutoffV*cutoffV))),2); // from Equation (4.3)
      double gA = g0A/pow((1.0 + (qsq/(cutoffA*cutoffA))),2); // " " "
      double gP = (2*mpro*gA)/(qsq + mpionsq); // " " "
      double gM = magmom*gV; // " " "
      hF[i] = -1*gV*gV; // from Equation (4.2)
      hGT[i] = (gA*gA) - ((gA*gP*qsq)/(3*mpro)) + (pow(gP*qsq,2)/(12*mprosq)) + ((gM*gM*qsq)/(6*mprosq)); // " " "
      //hGT[i] = (gA*gA)*( 1.0 - (2.0/3.0)*(qsq/(qsq + mpionsq)) + (1.0/3.0)*pow((qsq/(qsq + mpionsq)),2) ) + (pow(magmom,2)/(6.0*mprosq))*pow(gV,2)*qsq; // also works, le algebra
    }
    // create the TBME of M0nu
    for (auto& itmat : M0nu_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      double J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat = sqrt(2*J + 1); // the hat factor of J
      for (int ibra=0; ibra<nbras; ibra++)
      {
        Ket& bra = tbc_bra.GetKet(ibra); // get the final state = <ab| == <a_f,b_f|
        int ia = bra.p; // get the integer label a
        int ib = bra.q; // get the integer label b
        Orbit& oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit& ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        for (int iket=0; iket<nkets; iket++)
        {
          Ket& ket = tbc_ket.GetKet(iket); // get the initial state = |cd> == |a_i,b_i>
          int ic = ket.p; // get the integer label c
          int id = ket.q; // get the integer label d
          Orbit& oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit& od = modelspace.GetOrbit(id); // get the |d> state orbit
          // initialize the matrix element to zero
          double MF = 0;
          double MGT = 0;
          int na = oa.n; // this is just...
          int nb = ob.n;
          int nc = oc.n;
          int nd = od.n;
          int la = oa.l;
          int lb = ob.l;
          int lc = oc.l;
          int ld = od.l;
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience
std::cout<<J<<" | "<<ia<<", "<<ib<<", "<<ic<<", "<<id<<" || "
<<na<<", "<<nb<<", "<<nc<<", "<<nd<<" | "
<<la<<", "<<lb<<", "<<lc<<", "<<ld<<" | "
<<ja<<", "<<jb<<", "<<jc<<", "<<jd<<" |=|  ";
          int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
          int eps_cd = 2*nc + lc + 2*nd + ld; // ...likewise
          double sumglqF = 0; // for the GLQ integration below (F)
          double sumglqFas = 0; // ...anti-symmetric part (F)
          double sumglqGT = 0; // (GT)
          double sumglqGTas = 0; // (GT)
          //double sumGLQ = 0; // (F+GT)
          //double sumGLQas = 0; // (F+GT)
          for (int i=0; i<Nquad; i++)
          {
            double q = nodes[i][0]; // set the transfer momentum [MeV] via the impending GLQ
            double sumB = 0; // for the pure Bessel's Matrix Elemets (BMEs)
            double sumBas = 0; // ...anti-symmetric part
            double sumBSdS = 0; // for the scalar product of the Bessel's functions and Sigma dot Sigma (SdS)
            double sumBSdSas = 0; // ...anti-symmetric part
            for (int S=0; S<=1; S++) // sum over total spin...
            {
              for (int L = std::abs(la-lb); L <= la+lb; L++) // ...and sum over angular momentum coupled to l_{a_f} and l_{b_f}, NOTE: get same result if used l_{a_i} and l_{b_i} (good)
              {
                double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
                double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
                double nNJab = normab*modelspace.GetNineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
                double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
                double nNJcd = normcd*modelspace.GetNineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
                double nNJdc = normcd*modelspace.GetNineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // ...anti-symmetric part
                double sumMT = 0; // for the Moshinsky-transformed relative BMEs (RBMEs)
                double sumMTas = 0; // ...anti-symmetric part
                double tempmaxnr = floor((eps_ab - L)/2.0); // just for the limits below
                double tempmaxnpr = floor((eps_cd - L)/2.0); // " " " " "
                for (int nr = 0; nr <= tempmaxnr; nr++)
                {
                  double tempmaxNcom = tempmaxnr - nr; // just for the limits below
                  for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                  {
                    int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                    int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                    for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                    {
                      int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation () of my thesis
                      for (int npr = 0; npr <= tempmaxnpr; npr++) // npr = n'_r in my latex notation
                      {
                        double RBME; // calculate via functions below this Operator
                        if (src == argstr or src == cdbstr or src == masstr)
                        {
                          RBME = CPrbmeGen(modelspace,0,q,nr,lr,npr,lr,0,0)
                            - 2*cc*CPrbmeGen(modelspace,0,q,nr,lr,npr,lr,0,aa)
                            + 2*cc*bb*CPrbmeGen(modelspace,0,q,nr,lr,npr,lr,2,aa)
                            + cc*cc*CPrbmeGen(modelspace,0,q,nr,lr,npr,lr,0,2.0*aa)
                            - 2*bb*cc*cc*CPrbmeGen(modelspace,0,q,nr,lr,npr,lr,2,2.0*aa)
                            + cc*cc*bb*bb*CPrbmeGen(modelspace,0,q,nr,lr,npr,lr,4,2.0*aa); // testing...
                        }
                        else
                        {
                          RBME = CPrbmeGen(modelspace,0,q,nr,lr,npr,lr,0,0);
                          //RBME = CPrbmeGLQ(modelspace,Nquad,nodes,q,nr,lr,npr,lr,0,0);
                        }
                        double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                        double Di =  modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L);
                        double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // ...anti-symmetric part
                        double temphat = 1.0/sqrt(2*lr + 1); // just for the several lines below
                        double tempD = Df*Di; // " " " " " "
                        double tempDas = Df*asDi; // " " " " " "
                        sumMT += temphat*tempD*RBME; // perform the Talmi-Moshinsky transformation
                        sumMTas += temphat*tempDas*RBME; // ...anti-symmetric part
                      } // end of for-loop over: npr
                    } // end of for-loop over: lr
                  } // end of for-loop over: Ncom
                } // end of for-loop over: nr
                double temprod = nNJab*nNJcd*sumMT; // just for efficiency, only used in the lines below
                double temprodas = nNJab*nNJdc*sumMTas; // ...anti-symmetric part
                int tempSeval = 2*S*(S + 1) - 3; // eigenvalue of S, only used in the lines below
                sumB += temprod; // this completes Equation (4.64), modulo Jhat
                sumBas += temprodas; // ...anti-symmetric part
                sumBSdS += tempSeval*temprod; // this completes Equation (4.49), modulo Jhat
                sumBSdSas += tempSeval*temprodas; // ...anti-symmetric part
              } // end of for-loop over: L
            } // end of for-loop over: S
            //
            double integrandF = (q/(q + Ediff))*(hF[i]*sumB); // Fermi part only
            double integrandFas = (q/(q + Ediff))*(hF[i]*sumBas); // ...anti-symmetric part
            double integrandGT = (q/(q + Ediff))*(hGT[i]*sumBSdS); // Gamow-Teller part only
            double integrandGTas = (q/(q + Ediff))*(hGT[i]*sumBSdSas); // ...anti-symmetric part
            //
            /*
            double qsq = pow(q,2);
            double Enu = sqrt(qsq + mnusq);
            double integrandF = (qsq/(Enu + Ediff))*((hF[i]/Enu)*sumB);
            double integrandFas = (qsq/(Enu + Ediff))*((hF[i]/Enu)*sumBas);
            double integrandGT = (qsq/(Enu + Ediff))*((hGT[i]/Enu)*sumBSdS);
            double integrandGTas = (qsq/(Enu + Ediff))*((hGT[i]/Enu)*sumBSdSas);
            */
            //double integrand = (q/(q + Ediff))*(hF[i]*sumB + hGT[i]*sumBSdS); // F+GT full integrand
            //double integrandas = (q/(q + Ediff))*(hF[i]*sumBas + hGT[i]*sumBSdSas); // ...anti-symmetric part
            sumglqF += nodes[i][1]*integrandF; // perform the GLQ integration (F)
            sumglqFas += nodes[i][1]*integrandFas; // ...anti-symmetric part (F)
            sumglqGT += nodes[i][1]*integrandGT; // perform the GLQ integration (GT)
            sumglqGTas += nodes[i][1]*integrandGTas; // ...anti-symmetric part (GT)
            //sumGLQ += nodes[i][1]*integrand; // perform the GLQ integration (F+GT)
            //sumGLQas += nodes[i][1]*integrandas; // ...anti-symmetric part (F+GT)
          } // end of for-loop wrt: GLQ
          double scale = 1.0; // global factor to compare with JE
          double tempfact = scale*prefact*Jhat; // just for the lines below
          double tempnorm = cpNorm(ia,ib)*cpNorm(ic,id); // just for the lines below
          double temphase = modelspace.phase(jc + jd - J); // just for the lines below
          MF = tempnorm*tempfact*(sumglqF - temphase*sumglqFas); // see Equation () of my thesis, (F)
          MGT = tempnorm*tempfact*(sumglqGT - temphase*sumglqGTas); // see Equation () of my thesis, (GT)
          double Mtbme = MF + MGT; // see Equation () of my thesis
std::cout<<MF<<",  "<<MGT<<",  "<<Mtbme<<std::endl;
          M0nu_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to M_ab
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
    } // end of for-loop over: auto
    return M0nu_TBME;
  }
 

/// testing...
/// mm = 0,2,4
/// pp = 0,a,2a
  double CPrbmeGen(ModelSpace& modelspace, double rho, double x, int n, int l, int np, int lp, int mm, double pp)
  {
    x *= SQRT2; // this is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
    double b = 1.0/sqrt(M_NUCLEON*modelspace.GetHbarOmega()); // the oscillator length [MeV^-1], from Equation (1.95) of my thesis
    double bisq = 1.0/(b*b); // the inverse squared oscillator length
    double yolo = (x*x)/(4.0*(pp + bisq)); // commonly reoccuring argument in the RBMEs
    double normy = (1.0/pow(2,rho+1))*(1.0/pow(b,3))*sqrt((PI*gsl_sf_fact(n)*gsl_sf_fact(np))/(gsl_sf_gamma(n + l + 1.5)*gsl_sf_gamma(np + lp + 1.5))); // the normalization factor
    double sum = 0; // initialize the double sum over k, k'
    for (int k=0; k<=n; k++)
    {
      for (int kp=0; kp<=np; kp++)
      {
        int kappa = ((l + lp - rho + mm)/2.0) + k + kp; // rho should be divisible by 2! >:|
        double llkk = l + lp + 2*k + 2*kp;
        double hahaha = pow(b,-llkk)*pow((pp + bisq),-(llkk + 3 + mm)/2.0);
        double one = modelspace.phase(k + kp)/(gsl_sf_fact(k)*gsl_sf_fact(kp));
        double bci = gsl_sf_gamma(n + l + 1.5)/(gsl_sf_gamma(n - k + 1.0)*gsl_sf_gamma(l + k + 1.5)); // the first binomial coefficient, via taking ratios of Gamma functions
        double bcj = gsl_sf_gamma(np + lp + 1.5)/(gsl_sf_gamma(np - kp + 1.0)*gsl_sf_gamma(lp + kp + 1.5)); // the second binomial coefficient, " " " " " "
        double two = bci*bcj;
        double three = gsl_sf_fact(kappa);
        double four = gsl_sf_laguerre_n(kappa,rho+0.5,yolo); // from GSL -- L^a_n(x) = gsl_sf_laguerre_n(const int n, const double a, const double x)
        sum += hahaha*one*two*three*four; // perform the summation
      }
    }
    double rbme = normy*pow(b*x,rho)*exp(-yolo)*sum; // see Equation (4.61) in my thesis
    return rbme;
  }

/////////////// end of M0v functions from Charlie //////////////////////////





 // Orbital angular momentum squared L^2 in the relative coordinate.
 // This was written with the deuteron in mind. Not sure if it will be useful for other systems.
 Operator L2rel_Op(ModelSpace& modelspace)
 {
   Operator OpOut(modelspace, 0,0,0,2);

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
    TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;++ibra)
    {
     Ket & bra = tbc.GetKet(ibra);
     for (int iket=ibra;iket<nkets;++iket)
     {
       Ket & ket = tbc.GetKet(iket);

       Orbit & oa = modelspace.GetOrbit(bra.p);
       Orbit & ob = modelspace.GetOrbit(bra.q);
       Orbit & oc = modelspace.GetOrbit(ket.p);
       Orbit & od = modelspace.GetOrbit(ket.q);
    
       int na = oa.n;
       int nb = ob.n;
       int nc = oc.n;
       int nd = od.n;
    
       int la = oa.l;
       int lb = ob.l;
       int lc = oc.l;
       int ld = od.l;
    
       double ja = oa.j2/2.0;
       double jb = ob.j2/2.0;
       double jc = oc.j2/2.0;
       double jd = od.j2/2.0;
    
       int fab = 2*na + 2*nb + la + lb;
       int fcd = 2*nc + 2*nd + lc + ld;
       if (fab != fcd) continue; // L2 conserves all the quantum numbers
    
       double sa,sb,sc,sd;
       sa=sb=sc=sd=0.5;
    
       double L2rel=0;
    
       for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
       {
         for (int Sab=0; Sab<=1; ++Sab)
         {
           if ( std::abs(Lab-Sab)>J or Lab+Sab<J) continue;
    
           double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
           if (njab == 0) continue;
           int Scd = Sab;
           int Lcd = Lab;
           double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
           if (njcd == 0) continue;
    
           // Next, transform to rel / com coordinates with Moshinsky tranformation
           for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
           {
             for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
             {
               int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
               for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
               {
                  if (Lab<std::abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
                  // factor to account for antisymmetrization
    
                  int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
                  if ( asymm_factor ==0 ) continue;
    
                  int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
                  int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation
    
                  double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
    
                  if (std::abs(mosh_ab)<1e-8) continue;
    
                  for (int N_cd=std::max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
                  {
                    int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                    if (n_cd < 0) continue;
                    if  (n_ab != n_cd or N_ab != N_cd) continue;
    
                    double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                    if (std::abs(mosh_cd)<1e-8) continue;
    
                    double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                    L2rel += lam_ab*(lam_ab+1) * prefactor;
    
                  } // N_cd
               } // lam_ab
             } // Lam_ab
           } // N_ab
    
         } // Sab
       } // Lab


       OpOut.TwoBody.SetTBME(ch,ibra,iket,L2rel);
       OpOut.TwoBody.SetTBME(ch,iket,ibra,L2rel);

     }
   }
  }

  return OpOut;

 }



 Operator MinnesotaPotential( ModelSpace& modelspace )
 {
   Operator Vminnesota(modelspace, 0,0,0,2);

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
//   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
    TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;++ibra)
    {
     Ket& bra = tbc.GetKet(ibra);
     for (int iket=ibra;iket<nkets;iket++)
     {
       Ket& ket = tbc.GetKet(iket);
       double vminn = MinnesotaMatEl( modelspace, bra, ket, J);
       Vminnesota.TwoBody.SetTBME(ch,ibra,iket,vminn);
       Vminnesota.TwoBody.SetTBME(ch,iket,ibra,vminn);
     }
    }
   }
   return Vminnesota;
 }

// Minnesota potential of Thompson, Lemere and Tang,  Nuc. Phys.A 286 53 (1977)
// I have modified the triplet channel strength by multiplying by 0.2.
// Without this modification, it seems that the triplet channel does not have enough
// repulsion to support against collapse, which makes the potential not terribly
// useful for finite nuclei.
// This is mostly useful for benchmarking without needing large input files.
//
// On March 12 2019, I calculated O16 in a basis of hw=16 and emax=4, 6, and 8.
// I used IMSRG(2) with method=magnus  omega_norm_max=0.25, dsmax=0.5
// Results:
//  emax = 4                         emax = 6                  emax = 8
//  e1hf = 177.4709518              e1hf = 132.7170792          e1hf = 111.4959623
//  e2hf = -102.1762390             e2hf = -138.2451400         e2hf = -210.3298042
//  e3hf = 0.0000000                e3hf = 0.0000000            e3hf = 0.0000000    
//  EHF = 75.2947127                EHF  = -5.5280608           EHF = -98.8338419   
//  EIMSRG = -10.047908             EIMSRG = -69.170313         EIMSRG =  -144.9899081
//  Rp2HF = 9.5591609               Rp2HF = 15.1991899          Rp2HF =  20.8748671
//  Rp2IMSRG = 9.362856             Rp2IMSRG = 14.399832        Rp2IMSRG = 21.353270
// 
 double MinnesotaMatEl( ModelSpace& modelspace, Ket& bra, Ket& ket, int J )
 {
   double oscillator_b2 = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
   double VR = 200 ;
   double VT = -178. * 0.2; // scaled to make things not crazy for finite nuclei
   double VS = -91.85 ;
   double kR = 1.487 * oscillator_b2; // make things dimensionless
   double kT = 0.639 * oscillator_b2;
   double kS = 0.465 * oscillator_b2;

   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
   int nb = ob.n;
   int nc = oc.n;
   int nd = od.n;

   int la = oa.l;
   int lb = ob.l;
   int lc = oc.l;
   int ld = od.l;

   double ja = oa.j2*0.5;
   double jb = ob.j2*0.5;
   double jc = oc.j2*0.5;
   double jd = od.j2*0.5;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   if (std::abs(fab+fcd)%2 >0) return 0; // check parity conservation
//   if (std::abs(fab-fcd)>2) return 0; // p1*p2 only connects kets with delta N = 0,1

   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double Vminn=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
//       double isospin_factor = 1.0;
       if ( (oa.tz2 == ob.tz2) and Sab==1 ) continue;
//       if ( oa.tz2 != ob.tz2) isospin_factor = 0.5 * (oa.tz2 * oc.tz2);  // <pn|V|pn> = 1/2 (V(T=1) + V(T=0)).  <pn|V|np> = -<pn|V|pn>
       if ( std::abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<std::abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;

              // factor to account for antisymmetrization. I sure wish this were more transparent...
              //  antisymmetrized matrix elements: V = <ab|V|cd> - <ab|V|dc>
              //  here, we've gone to LS, coupling in CM/relative coordinates
              //  so exchange gives us a phase factor and we have V = <ab|V|cd> + (-1)**(Lcd+Scd) <ab|V|dc>
              //  case pp,nn -> 1 + (-1)**(L+S)
              //  pnpn -> 1 + 0
              //  pnnp -> 0 + (-1)**(L+S)
              int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // V conserves lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (std::abs(mosh_ab)<1e-8) continue;

              for (int N_cd=std::max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd and N_ab != N_cd) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (std::abs(mosh_cd)<1e-8) continue;

                double vrel = 0;
                int pmin = lam_ab;
                int pmax = lam_ab + n_ab + n_cd;
                for (int p=pmin; p<=pmax; p++)
                {
                   double Bp = TalmiB(n_ab,lam_ab,n_cd,lam_cd,p);  // It's probably inefficient doing this in the inner loop, but oh well...
                   double Ip = VR / pow(1.0+kR, p+1.5);  // Talmi integral of a Gaussian
                   if ( Sab==1 ) Ip += VT / pow(1.0+kT, p+1.5); // Triplet short-range potential
                   else Ip += VS / pow(1.0+kS, p+1.5); // Singlet short-range potential
                   vrel += Bp * Ip;
                }
                Vminn += vrel * njab * njcd * mosh_ab * mosh_cd * asymm_factor;


              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   return Vminn ;



 }






 Operator EKKShift( Operator& Hin, int Nlower, int Nupper)
 {
   ModelSpace* modelspace = Hin.GetModelSpace();
   Operator EKKShift(*modelspace,0,0,0,2);
   double elower = 0;
   double eupper = 0;
   std::vector<index_t> index_lower;
   std::vector<index_t> index_upper;
   for (size_t i=0; i<modelspace->GetNumberOrbits(); ++i)
   {
     Orbit& oi = modelspace->GetOrbit(i);
     int N = 2*oi.n + oi.l;
     if (N==Nlower)
     {
       elower += Hin.OneBody(i,i);
       index_lower.push_back(i);
     }
     if (N==Nupper)
     {
       eupper += Hin.OneBody(i,i);
       index_upper.push_back(i);
     }
   }
   elower /= index_lower.size();
   eupper /= index_upper.size();
   double emid = 0.5 * (elower + eupper);
   for ( auto i : index_lower ) EKKShift.OneBody(i,i) = emid - elower;
   for ( auto i : index_upper ) EKKShift.OneBody(i,i) = emid - eupper;
   std::cout << "In EKKShift, elower, eupper, emid = " << elower << " , " << eupper << " , " << emid << std::endl;
   return EKKShift;
 }


  std::map<index_t,double> GetSecondOrderOccupations(Operator& H, int emax)
  {
//    ModelSpace* modelspace = H.GetModelSpace();
    std::map<index_t,double> hole_list;
    std::cout << "GetSecondOrderOccupations : Not yet implemented" << std::endl;
    return hole_list;
  }


  /// Embeds the one-body operator of op1 in the two-body part, using mass number A in the embedding.
  /// Note that the embedded operator is added to the two-body part, rather than overwriting.
  /// The one-body part is left as-is.
  void Embed1BodyIn2Body(Operator& op1, int A)
  {
    if (A<2)
    {
      std::cout << "Embed1BodyIn2Body: A = " << A << ". You clearly didn't mean to do that..." << std::endl;
      return;
    }
    ModelSpace* modelspace = op1.GetModelSpace();
    int Lambda = op1.GetJRank();
    for (auto& itmat : op1.TwoBody.MatEl)
    {
      index_t ch_bra = itmat.first[0];
      index_t ch_ket = itmat.first[1];
      arma::mat& TBME = itmat.second;

      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      index_t nbras = tbc_bra.GetNumberKets();
      index_t nkets = tbc_ket.GetNumberKets();
      int Jbra = tbc_bra.J;
      int Jket = tbc_ket.J;
      for (index_t ibra=0;ibra<nbras;++ibra)
      {
        Ket& bra = tbc_bra.GetKet(ibra);
        index_t i = bra.p;
        index_t j = bra.q;
        for (index_t iket=0;iket<nkets;++iket)
        {
          Ket& ket = tbc_ket.GetKet(iket);
          index_t k = ket.p;
          index_t l = ket.q;
          double embedded_tbme = GetEmbeddedTBME(op1,i,j,k,l,Jbra,Jket,Lambda) / (A-1.0);
          TBME(ibra,iket) += embedded_tbme;
        }
      }
    }
  }

  /// Returns a normalized TBME formed by embedding the one body part of op1 in a two body operator. Assumes A=2 in the formula.
  /// For other values of A, divide by (A-1).
  double GetEmbeddedTBME(Operator& op1, index_t i, index_t j, index_t k, index_t l, int Jbra,int Jket, int Lambda)
  {
    double embedded_tbme = 0;
    ModelSpace* modelspace = op1.GetModelSpace();
    double ji = modelspace->GetOrbit(i).j2 * 0.5;
    double jj = modelspace->GetOrbit(j).j2 * 0.5;
    double jk = modelspace->GetOrbit(k).j2 * 0.5;
    double jl = modelspace->GetOrbit(l).j2 * 0.5;
    arma::mat& OB = op1.OneBody;
    if (Lambda==0) // scalar => no sixJs, tbmes are not reduced.
    {
       if (j==l)  embedded_tbme += OB(i,k);
       if (i==k)  embedded_tbme += OB(j,l);
       if (i==l)  embedded_tbme -= OB(j,k) * modelspace->phase(ji+jj-Jbra);
       if (j==k)  embedded_tbme -= OB(i,l) * modelspace->phase(ji+jj-Jbra);
    }
    else // Tensor => slightly more complicated, tbmes are reduced.
    {
       if (j==l)  embedded_tbme += OB(i,k) * modelspace->phase(ji+jj+Jket) * AngMom::SixJ(Jbra,Jket,Lambda,jk,ji,jj);
       if (i==k)  embedded_tbme += OB(j,l) * modelspace->phase(jk+jl-Jbra) * AngMom::SixJ(Jbra,Jket,Lambda,jl,jj,ji);
       if (j==k)  embedded_tbme -= OB(i,l) * modelspace->phase(ji+jj+jk+jl)* AngMom::SixJ(Jbra,Jket,Lambda,jl,ji,jj);
       if (i==l)  embedded_tbme -= OB(j,k) * modelspace->phase(Jket-Jbra)  * AngMom::SixJ(Jbra,Jket,Lambda,jk,jj,ji);
       embedded_tbme *= sqrt((2*Jbra+1)*(2*Jket+1)) * modelspace->phase(Lambda);
    }
    if (i==j) embedded_tbme /=SQRT2;
    if (k==l) embedded_tbme /=SQRT2;
    return embedded_tbme;
  }



//  double FCcoefFunction(int n1, int l1, double hw1, int n2, int l2, double hw2, double r)
  double FCcoefIntegrand(double r, void* params)
  {
    int n1     = ((double*) params)[0];
    int l1     = ((double*) params)[1];
    double hw1 = ((double*) params)[2];
    int n2     = ((double*) params)[3];
    int l2     = ((double*) params)[4];
    double hw2 = ((double*) params)[5];
    return r*r*HO_Radial_psi( n1,l1,hw1,r)*HO_Radial_psi( n2,l2,hw2,r);
  }

  double FrequencyConversionCoeff(int n1, int l1, double hw1, int n2, int l2, double hw2)
  {
    double parameters[6] = {(double)n1, (double)l1, hw1, (double)n2, (double)l2, hw2};
    gsl_function F;
    F.function = &FCcoefIntegrand;
    F.params = &(parameters[0]);
    double r0 = 0;
    double epsabs = 1e-7;
    double epsrel = 1e-7;
    size_t limit = 1000; // is this a reasonable number?
    gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(limit);
    double result = 0;
    double abserr;

    gsl_integration_qagiu(&F, r0, epsabs, epsrel, limit, workspace, &result, &abserr );

    gsl_integration_workspace_free(workspace);
    return result;

  }



//  vector<double> FindLaguerreRoots( int n )
//  {
//    // For a Laguerre polynominal of order n, we should find n roots.
//    // The kth derivative of L is given by d^k/dx^k L^{0}_{n}(x) = (-1)^k L^{k}_{n-k}(x)  if n>k and 0 otherwise
//    vector<double> roots;
//    double x0 = 0;
//    double x1 = 0;
//    double f1 = gsl_sf_laguerre_n( n, 0,  x0);
//    double fprime1 = n>1 ? : gsl_sf_laguerre_n( n-1,1, x0 ) 0 ;
//  }



  void CommutatorTest(Operator& X, Operator& Y)
  {
    Operator Zscalar(X);
    if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Zscalar.SetAntiHermitian();
    if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Zscalar.SetHermitian();
    Zscalar.Erase();
    Operator Ztensor(Zscalar);
    Operator Yred = Y;
    Reduce(Yred);

    std::cout << "operator norms: " << X.Norm() << "  " << Y.Norm() << std::endl;
//    X.comm111ss(Y,Zscalar);
//    X.comm111st(Yred,Ztensor);
    Commutator::comm111ss(X,Y,Zscalar);
    Commutator::comm111st(X,Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm111 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm111 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;

    Zscalar.Erase();
    Ztensor.Erase();
    std::cout << "121ss" << std::endl;
    Commutator::comm121ss(X,Y,Zscalar);
    std::cout << "121st" << std::endl;
    Commutator::comm121st(X,Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm121 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm121 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;

    Zscalar.Erase();
    Ztensor.Erase();
    Commutator::comm122ss(X,Y,Zscalar);
    Commutator::comm122st(X,Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm122 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm122 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;

    Zscalar.Erase();
    Ztensor.Erase();
    Commutator::comm222_pp_hh_221ss(X,Y,Zscalar);
    Commutator::comm222_pp_hh_221st(X,Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm222_pp_hh_221 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm222_pp_hh_221 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;

    Zscalar.Erase();
    Ztensor.Erase();
    Commutator::comm222_phss(Y,Zscalar,X);
//    Reduce(Y); // Not sure why I can't use Yred...
    Commutator::comm222_phss(X,Y,Zscalar);
    Commutator::comm222_phst(X,Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm222_ph norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm222_ph diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;


  }



/*
  void CommutatorTest(Operator& X, Operator& Y)
  {
    Operator Zscalar(X);
    if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Zscalar.SetAntiHermitian();
    if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Zscalar.SetHermitian();
    Zscalar.Erase();
    Operator Ztensor(Zscalar);
    Operator Yred = Y;
    Reduce(Yred);

    std::cout << "operator norms: " << X.Norm() << "  " << Y.Norm() << std::endl;
//    X.comm111ss(Y,Zscalar);
//    X.comm111st(Yred,Ztensor);
    Zscalar.comm111ss(X,Y);
    Ztensor.comm111st(X,Yred);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm111 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm111 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm121ss(Y,Zscalar);
    X.comm121st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm121 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm121 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm122ss(Y,Zscalar);
    X.comm122st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm122 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm122 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm222_pp_hh_221ss(Y,Zscalar);
    X.comm222_pp_hh_221st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm222_pp_hh_221 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm222_pp_hh_221 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm222_phss(Y,Zscalar);
    Reduce(Y); // Not sure why I can't use Yred...
    X.comm222_phst(Y,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    std::cout << "comm222_ph norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << std::endl;
    Zscalar -= Ztensor;
    std::cout << "comm222_ph diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << std::endl;


  }
*/




// template <typename T>
// T VectorUnion(T& v1)
// {
//   return v1;
// }
// 
// template <typename T, typename... Args>
// T VectorUnion(T& v1, T& v2, Args... args)
// {
//   T vec(v1.size()+v2.size());
//   copy(v1.begin(),v1.end(),vec.begin());
//   copy(v2.begin(),v2.end(),vec.begin()+v1.size());
//   return VectorUnion(vec, args...);
// }
 

}// namespace imsrg_util

