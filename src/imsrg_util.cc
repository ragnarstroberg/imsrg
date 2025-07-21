
#include "imsrg_util.hh"
#include "AngMom.hh"
#include "Commutator.hh"
#include "GaussLaguerre.hh"
#include "DarkMatterNREFT.hh"
#include "M0nu.hh"
#include "omp.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h> // to use bessel functions
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <array>


/// imsrg_util namespace. Used to define some helpful functions.
namespace imsrg_util
{
 using PhysConst::HBARC;
 using PhysConst::M_PROTON;
 using PhysConst::M_NEUTRON;
 using PhysConst::M_NUCLEON;
 using PhysConst::M_ELECTRON;
 using PhysConst::PROTON_SPIN_G;
 using PhysConst::NEUTRON_SPIN_G;
 using PhysConst::ELECTRON_SPIN_G;
 using PhysConst::F_PI;
 using PhysConst::ALPHA_FS;
 using PhysConst::PI;
 using PhysConst::SQRT2;
 using PhysConst::SQRTPI;
 using PhysConst::INVSQRT2;
 using PhysConst::LOG2;

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

 /// Convenient way to generate pre-coded operators
 Operator OperatorFromString(ModelSpace& modelspace, std::string opname )
 {
      double t_start = omp_get_wtime();
      Operator theop;
      std::vector<std::string> opnamesplit = split_string( opname, "_" );  // split std::string on _ into a vector of std::string so that, e.g. "R2_p1"  =>  {"R2", "p1"}

           if (opname == "R2_p1")         theop =  R2_1body_Op(modelspace,"proton") ;
      else if (opname == "R2_p2")         theop =  R2_2body_Op(modelspace,"proton") ;
      else if (opname == "R2_n1")         theop =  R2_1body_Op(modelspace,"neutron") ;
      else if (opname == "R2_n2")         theop =  R2_2body_Op(modelspace,"neutron") ;
      else if (opname == "Rp2")           theop =  Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "Rn2")           theop =  Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "Rm2")           theop =  Rm2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "Rm2lab")        theop =  RSquaredOp(modelspace) ;
      else if (opname == "ISM")           theop =  MultipoleResponseOp(modelspace,2,0,0) ; // Isoscalar monopole (see PRC97(2018)054306 ) --added by bhu
      else if (opname == "IVM")           theop =  MultipoleResponseOp(modelspace,2,0,1) ; // Isovector monopole
      else if (opname == "ISQ")           theop =  MultipoleResponseOp(modelspace,2,2,0) ; // Isoscalar quadrupole
      else if (opname == "IVQ")           theop =  MultipoleResponseOp(modelspace,2,2,1) ; // Isosvector quadrupole
      else if (opname == "ISO")           theop =  MultipoleResponseOp(modelspace,3,3,0) ; // Isoscalar octupole
      else if (opname == "IVO")           theop =  MultipoleResponseOp(modelspace,3,3,1) ; // Isosvector octupole
      else if (opname == "IVD")           theop =  IVDipoleOp(modelspace,1,1) ; // Isovector dipole (see PRC97(2018)054306 ) --added by bhu
      else if (opnamesplit[0] == "ISD")  // Isoscalar dipole, the CoM term are 5/3*Rms^2*r e.g., ISD_3.2111, here Rms=3.2111  --added by bhu
      {
          double Rms ; // center-of-mass correction terms
          std::istringstream( opnamesplit[1] ) >> Rms ;
          theop =  ISDipoleOp(modelspace,3,1,Rms) ;
      }
      else if (opname == "E0")            theop =  E0Op(modelspace) ;
      else if (opname == "E1")            theop =  ElectricMultipoleOp(modelspace,1) ;
      else if (opname == "E2")            theop =  ElectricMultipoleOp(modelspace,2) ;
      else if (opname == "E3")            theop =  ElectricMultipoleOp(modelspace,3) ;
      else if (opname == "E4")            theop =  ElectricMultipoleOp(modelspace,4) ;
      else if (opname == "E5")            theop =  ElectricMultipoleOp(modelspace,5) ;
      else if (opname == "E6")            theop =  ElectricMultipoleOp(modelspace,6) ;
      else if (opname == "E2int")         theop =  IntrinsicElectricMultipoleOp(modelspace,2) ; // Untested
      else if (opname == "nE2")           theop =  NeutronElectricMultipoleOp(modelspace,2) ;
      else if (opname == "M1")            theop =  MagneticMultipoleOp(modelspace,1) ;
      else if (opname == "M2")            theop =  MagneticMultipoleOp(modelspace,2) ;
      else if (opname == "M3")            theop =  MagneticMultipoleOp(modelspace,3) ;
      else if (opname == "M4")            theop =  MagneticMultipoleOp(modelspace,4) ;
      else if (opname == "M5")            theop =  MagneticMultipoleOp(modelspace,5) ;
      else if (opname == "M1p")           theop =  MagneticMultipoleOp_pn(modelspace,1,"proton") ;
      else if (opname == "M1n")           theop =  MagneticMultipoleOp_pn(modelspace,1,"neutron") ;
      else if (opname == "M1S")           theop =  MagneticMultipoleOp_pn(modelspace,1,"spin") ;
      else if (opname == "M1L")           theop =  MagneticMultipoleOp_pn(modelspace,1,"orbit") ;
      else if (opname == "Fermi")         theop =  AllowedFermi_Op(modelspace) ;
      else if (opname == "GamowTeller")   theop =  AllowedGamowTeller_Op(modelspace) ;
      else if (opname == "Iso2")          theop =  Isospin2_Op(modelspace) ;
      else if (opname == "Tz2")           theop =  TzSquared_Op(modelspace) ;
      else if (opname == "R2CM")          theop =  R2CM_Op(modelspace) ;
      else if (opname == "Trel")          theop =  Trel_Op(modelspace) ;
      else if (opname == "TCM")           theop =  TCM_Op(modelspace) ;
      else if (opname == "Tlab")          theop =  KineticEnergy_Op(modelspace) ;
      else if (opname == "TrelMassCorrection")          theop =  Trel_Masscorrection_Op(modelspace) ;
      else if (opname == "Rso")           theop =  RpSpinOrbitCorrection(modelspace) ;
      else if (opname == "RadialOverlap") theop =  RadialOverlap(modelspace); // Untested...
      else if (opname == "Sigma")         theop =  Sigma_Op(modelspace);
      else if (opname == "Sigma_p")       theop =  Sigma_Op_pn(modelspace,"proton");
      else if (opname == "Sigma_n")       theop =  Sigma_Op_pn(modelspace,"neutron");
      else if (opname == "L2rel")         theop =  L2rel_Op(modelspace); // Untested...
      else if (opname == "QdotQ")         theop =  QdotQ_Op(modelspace); // Untested...
      else if (opname == "VQQ")           theop =  VQQ_Op(modelspace); 
      else if (opname == "VCoul")         theop =  VCoulomb_Op(modelspace); // Untested...
      else if (opname == "hfsNMS")        theop =  atomic_hfs::NormalMassShift(modelspace, 1);
      else if (opname == "hfsSMS")        theop =  atomic_hfs::SpecificMassShift(modelspace, 1);
      else if (opname == "VCentralCoul")  theop =  VCentralCoulomb_Op(modelspace); 
      else if (opname == "AxialCharge")   theop =  AxialCharge_Op(modelspace); // Untested...
      else if (opname == "VMinnesota")    theop =  MinnesotaPotential( modelspace );
      else if (opname == "VBareDelta")    theop =  BareDelta( modelspace );
      else if (opname == "OccRef")        theop =  NumberOpRef( modelspace );
      else if (opname == "LdotS")         theop =  LdotS_Op( modelspace);
      else if (opname == "DGT")           theop = M0nu::DGT_Op(modelspace);
      else if (opname == "Anapole")       theop = AnapoleMoment(modelspace);
      else if (opname == "ChargeDensity") theop = Charge_Density_Op(modelspace);
      else if (opnamesplit[0] =="VGaus")
      {
         double sigma = 1.0;
         if ( opnamesplit.size() > 1 ) 
         {
           std::istringstream( opnamesplit[1] ) >> sigma;
         }
         theop =  GaussianPotential(modelspace,sigma);
      }
      else if (opnamesplit[0] =="VGausCSB")
      {
         double sigma = 1.0;
         if ( opnamesplit.size() > 1 ) 
         {
           std::istringstream( opnamesplit[1] ) >> sigma;
         }
         theop =  GaussianCSBPotential(modelspace,sigma);
      }
      else if (opnamesplit[0] =="VSDI")
      {
//        std::cout << "    " << opnamesplit[0] << " " << opnamesplit[1] << " " << opnamesplit[2] << std::endl;
         double V0 = 1.0;
         double R = 1.0;
         if ( opnamesplit.size() > 2 ) 
         {
           std::istringstream( opnamesplit[1] ) >> V0;
           std::istringstream( opnamesplit[2] ) >> R;
         }
         theop =  SurfaceDeltaInteraction(modelspace,V0,R);
      }
      else if (opnamesplit[0] =="HCM")
      {
//         if ( opnamesplit.size() == 1 ) theop =  HCM_Op(modelspace);
//         else
//         {
           double hw_HCM = modelspace.GetHbarOmega(); // frequency of trapping potential
           if (opnamesplit.size()>1)
           {
              std::istringstream( opnamesplit[1] ) >> hw_HCM;
           }
           int A = modelspace.GetTargetMass();
//           std::cout << "Calling HCM with hw = " << hw_HCM << " target mass = " << A << std::endl;
           theop =  TCM_Op(modelspace) + 0.5*A*M_NUCLEON*hw_HCM*hw_HCM/HBARC/HBARC*R2CM_Op(modelspace);
//         }
      }
      else if (opnamesplit[0] == "VCM") // GetHCM with a different frequency, ie HCM_24 for hw=24
      {
         double hw_VCM = modelspace.GetHbarOmega(); // frequency of trapping potential
         if ( opnamesplit.size() > 1 )
         {
            std::istringstream(opnamesplit[1]) >> hw_VCM;
         }
         int A = modelspace.GetTargetMass();
         theop =  0.5*A*M_NUCLEON*hw_VCM*hw_VCM/HBARC/HBARC*R2CM_Op(modelspace); 
      }
      else if (opnamesplit[0] == "Rp2Z") // Get point proton radius for specified Z, e.g. Rp2Z_10 for neon
      {
        int Z_rp;
        std::istringstream(opnamesplit[1]) >> Z_rp;
        theop =  Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) ;
      }
      else if (opnamesplit[0] == "Rp2AZ") // Get point proton radius for specified A and Z, e.g. Rp2AZ_20_10 for neon
      {
        int A_rp;
        int Z_rp;
        std::istringstream(opnamesplit[1]) >> A_rp;
        std::istringstream(opnamesplit[2]) >> Z_rp;
        theop =  Rp2_corrected_Op(modelspace,A_rp,Z_rp) ;
      }
      else if (opnamesplit[0] == "Rn2Z") // Get point neutron radius for specified Z
      {
        int Z_rp;
        std::istringstream(opnamesplit[1]) >> Z_rp;
        theop =  Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) ;
      }
      else if (opnamesplit[0] == "rhop") // point proton  density at position r, e.g. rhop_1.25
      {
        double rr;
        std::istringstream(opnamesplit[1]) >> rr;
        theop =  DensityAtR(modelspace,rr,"proton");
      }
      else if (opnamesplit[0] == "rhon") // point neutron density at position r
      {
        double rr;
        std::istringstream(opnamesplit[1]) >> rr;
        theop =  DensityAtR(modelspace,rr,"neutron"); // whoops... I'd forgotten the "theop = ". Thanks Johannes...
      }
      else if (opnamesplit[0] == "rhocentralp") // central point proton density, averaged over a sphere of radius rmax
      {
        double rmax;
        std::istringstream(opnamesplit[1]) >> rmax;
        int npoints = 10;
        double dr = rmax/npoints;
        theop = 0.5*0.5*DensityAtR(modelspace,dr/2,"proton");
        double norm = 0.5*0.5;
        for ( int ipoint=1; ipoint<npoints; ipoint++)  theop += (ipoint+0.5)*(ipoint+0.5) * DensityAtR(modelspace,dr*(ipoint+0.5),"proton");
        for ( int ipoint=1; ipoint<npoints; ipoint++)  norm  += (ipoint+0.5)*(ipoint+0.5) ;
        theop /= norm;
      }
      else if (opnamesplit[0] == "rhocentraln")// central point neutron density, averaged over a sphere of radius rmax
      {
        double rmax;
        std::istringstream(opnamesplit[1]) >> rmax;
        int npoints = 10;
        double dr = rmax/npoints;
        theop = 0.5*0.5*DensityAtR(modelspace,dr/2,"neutron");
        double norm = 0.5*0.5;
        for ( int ipoint=1; ipoint<npoints; ipoint++)  theop += (ipoint+0.5)*(ipoint+0.5) * DensityAtR(modelspace,dr*(ipoint+0.5),"neutron");
        for ( int ipoint=1; ipoint<npoints; ipoint++)  norm  += (ipoint+0.5)*(ipoint+0.5) ;
        theop /= norm;
      }
      else if (opnamesplit[0] == "FFp") // point proton  form factor F(q), e.g. FFp_1.25, where q is in fm^-1
      {
        double q;
        std::istringstream(opnamesplit[1]) >> q;
        theop =  FormfactorAtQ(modelspace,q,"proton");
      }
      else if (opnamesplit[0] == "FFn") // point neutron  form factor F(q), e.g. FFp_1.25, where q is in fm^-1
      {
        double q;
        std::istringstream(opnamesplit[1]) >> q;
        theop =  FormfactorAtQ(modelspace,q,"neutron");
      }
      else if (opnamesplit[0] == "OneOcc") // Get occupation of specified orbit, e.g. OneOccp_1p3
      {
         index_t ind = modelspace.String2Index( {  opnamesplit[1] } )[0];
         Orbit& oi = modelspace.GetOrbit(ind);
         theop =  NumberOp(modelspace,oi.n,oi.l,oi.j2,oi.tz2) ;
      }
      else if (opnamesplit[0]== "AllOcc") // Get occupation of orbit, summed over all values of radial quantum number n, e.g. AllOccpp3
      {
         index_t ind = modelspace.String2Index( { "0"+ opnamesplit[1] } )[0];
         Orbit& oi = modelspace.GetOrbit(ind);
         theop =  NumberOpAlln(modelspace,oi.l,oi.j2,oi.tz2) ;
      }
      else if (opnamesplit[0] == "OBD")
      {
         index_t i = modelspace.String2Index( { opnamesplit[1] } )[0];
         index_t j = modelspace.String2Index( { opnamesplit[2] } )[0];
         theop =  OneBodyDensity(modelspace,i,j);
      }
      else if (opnamesplit[0] == "protonFBC") // Fourier bessel coefficient of order nu
      {
         int nu;
         std::istringstream(opnamesplit[1]) >> nu;
         theop =  FourierBesselCoeff( modelspace, nu, 8.0, modelspace.proton_orbits);
      }
      else if (opnamesplit[0] == "neutronFBC") // Fourier bessel coefficient of order nu
      {
         int nu;
         std::istringstream(opnamesplit[1]) >> nu;
         theop =  FourierBesselCoeff( modelspace, nu, 8.0, modelspace.neutron_orbits) ;
      }
      else if (opnamesplit[0] == "M0nuCT" )
      {
        double R0;
        if (opnamesplit.size() < 2 ) std::cout << "ERROR!!!  " << __func__ << "    need to specify a cutoff for M0nuCT" << std::endl;
        std::istringstream( opnamesplit[1]) >> R0;
        theop =  M0nu_contact_Op(modelspace, R0);
      }

      else if (opnamesplit[0] == "DMNREFT") // Dark matter non-relativistic EFT operators: M+, M-, ...
      {
        std::string IsoSV; // "+" labels isoscalar and "-" is isovector
        double q;
        int J;
        std::string dmopname = opnamesplit[1].substr(0, opnamesplit[1].length() - 1);
        IsoSV = opnamesplit[1].back();
        std::istringstream(opnamesplit[2]) >> q;
        std::istringstream(opnamesplit[3]) >> J;

        std::map<std::string, Operator (*)(ModelSpace&, std::string, int, double) > dmop = {
              {"M",       &DM_NREFT::M},
              {"Sigma",   &DM_NREFT::Sigma},
              {"Sigmap",  &DM_NREFT::Sigmap},
              {"Sigmapp", &DM_NREFT::Sigmapp},
              {"Delta",   &DM_NREFT::Delta},
              {"Deltap",  &DM_NREFT::Deltap},
              {"Phip",    &DM_NREFT::Phip},
              {"Phipp",   &DM_NREFT::Phipp},
              {"Phitp",   &DM_NREFT::Phitp},
              {"Omega",   &DM_NREFT::Omega},
              {"Omegat",  &DM_NREFT::Omegat},
             };
        if ( dmop.find(dmopname) != dmop.end() )
        {
        theop =  dmop[dmopname](modelspace, IsoSV, J, q );
        }
      }
      else if (opnamesplit[0] == "Dagger" or opnamesplit[0] == "DaggerHF" )
      {
        index_t Q = modelspace.String2Index({opnamesplit[1]})[0];
        std::cout << "call Dagger_Op with Q = " << Q << std::endl;
        theop =  Dagger_Op( modelspace, Q);
      }
      else if (opnamesplit[0] == "DaggerAlln")
      {
        index_t Q = modelspace.String2Index({opnamesplit[1]})[0];
        theop =  DaggerAlln_Op( modelspace, Q);
      }
      else if (opnamesplit[0] == "M0nu") // Neutrinoless Double Beta Decay Operators   format e.g.  M0nu_GT_7.72_none or M0nu_F_12.6_AV18
      {
        
        std::string M0nuopname = opnamesplit[1];
        std::map<std::string, Operator (*)(ModelSpace&, double, std::string, std::function<double(double)>)> OPList{
            {"GT", &M0nu::GamowTeller},
            {"F", &M0nu::Fermi},
            {"T", &M0nu::Tensor}};

        if (M0nuopname != "C")
        {
          double Eclosure;
          std::istringstream(opnamesplit[2]) >> Eclosure;
          std::string src = opnamesplit[3];
          if (opnamesplit.size() == 4)
          {
            std::map<std::string, std::function<double(double)>> FormFactorList{
                {"GT",M0nu::GTFormFactor},
                {"F", M0nu::FermiFormFactor},
                {"T", M0nu::TensorFormFactor}
            };
            theop = OPList[M0nuopname](modelspace, Eclosure, src, FormFactorList[M0nuopname]);
          }
          else if (opnamesplit.size() == 5)
          { 
            std::string formfactor = opnamesplit[4];
            if (M0nuopname == "GT")
            {
              if (formfactor == "AA")
              {
                theop = M0nu::GamowTeller(modelspace, Eclosure, src, M0nu::hGT_AA);
              }
              else if (formfactor == "AP")
              {
                theop = M0nu::GamowTeller(modelspace, Eclosure, src, M0nu::hGT_AP);
              }
              else if (formfactor == "PP")
              {
                theop = M0nu::GamowTeller(modelspace, Eclosure, src, M0nu::hGT_PP);
              }
              else if (formfactor == "MM")
              {
                theop = M0nu::GamowTeller(modelspace, Eclosure, src, M0nu::hGT_MM);
              }
            }
            else if (M0nuopname == "T")
            {
              if (formfactor == "AA")
              {
                theop = M0nu::Tensor(modelspace, Eclosure, src, M0nu::hT_AA);
              }
              else if (formfactor == "AP")
              {
                theop = M0nu::Tensor(modelspace, Eclosure, src, M0nu::hT_AP);
              }
              else if (formfactor == "PP")
              {
                theop = M0nu::Tensor(modelspace, Eclosure, src, M0nu::hT_PP);
              }
              else if (formfactor == "MM")
              {
                theop = M0nu::Tensor(modelspace, Eclosure, src, M0nu::hT_MM);
              }
            }
          }
        }
        else
        {
          double regulator_cutoff;
          std::istringstream(opnamesplit[2]) >> regulator_cutoff;
          int regulator_power;
          std::istringstream(opnamesplit[3]) >> regulator_power;
          theop = M0nu::Contact(modelspace, regulator_cutoff, regulator_power);
        }
        
      }
      else if (opnamesplit[0] == "M0nuHeavy") // Neutrinoless Double Beta Decay Operators   format e.g.  M0nu_GT_7.72_none or M0nu_F_12.6_AV18
      {
        std::string M0nuopname = opnamesplit[1];
        std::string src = opnamesplit[2];
        std::string formfactor = opnamesplit[3];
        if (M0nuopname == "GT")
        {
          if (formfactor == "AA")
          {
            theop = M0nu::GamowTellerHeavy(modelspace, src, M0nu::hGT_AA);
          }
          else if (formfactor == "AP")
          {
            theop = M0nu::GamowTellerHeavy(modelspace, src, M0nu::hGT_AP);
          }
          else if (formfactor == "PP")
          {
            theop = M0nu::GamowTellerHeavy(modelspace, src, M0nu::hGT_PP);
          }
        }
        else if (M0nuopname == "F")
        {
          theop = M0nu::FermiHeavy(modelspace, src, M0nu::hF_VV);
        }
        else if (M0nuopname == "T")
        {

          if (formfactor == "AP")
          {
            theop = M0nu::TensorHeavy(modelspace, src, M0nu::hT_AP);
          }
          else if (formfactor == "PP")
          {
            theop = M0nu::TensorHeavy(modelspace, src, M0nu::hT_PP);
          }
        }
      }
      else if (opnamesplit[0] == "M0nuIntegrand") // Neutrinoless Double Beta Decay Operators   format e.g.  M0nu_GT_7.72_none or M0nu_F_12.6_AV18
      {
        std::string M0nuopname = opnamesplit[1];
        std::map<std::string, Operator (*)(ModelSpace &, double, std::string, std::function<double(double)>, int, int)> OPList{
            {"GT", &M0nu::GamowTeller_integrand},
            {"F", &M0nu::Fermi_integrand},
            {"T", &M0nu::Tensor_integrand}};
        if (M0nuopname != "C")
        {
          double Eclosure;
          int p;
          int pp;
          std::istringstream(opnamesplit[2]) >> Eclosure;
          std::string src = opnamesplit[3];
          std::istringstream(opnamesplit[4]) >> p;
          std::istringstream(opnamesplit[5]) >> pp;
          std::map<std::string, std::function<double(double)>> FormFactorList{
              {"GT", M0nu::GTFormFactor},
              {"F", M0nu::FermiFormFactor},
              {"T", M0nu::TensorFormFactor}};
          theop = OPList[M0nuopname](modelspace, Eclosure, src, FormFactorList[M0nuopname], p, pp); 
        }
        else
        {
          double regulator_cutoff;
          std::istringstream(opnamesplit[2]) >> regulator_cutoff;
          int regulator_power;
          std::istringstream(opnamesplit[3]) >> regulator_power;
          int p;
          int pp;
          std::istringstream(opnamesplit[4]) >> p;
          std::istringstream(opnamesplit[5]) >> pp;
          theop = M0nu::Contact_integrand(modelspace, regulator_cutoff, regulator_power, p, pp);
        }
      }

      else if (opnamesplit[0] == "VWS")
      {
         double V0,R0,a0;
         std::istringstream( opnamesplit[1] ) >> V0;
         std::istringstream( opnamesplit[2] ) >> R0;
         std::istringstream( opnamesplit[3] ) >> a0;
         theop = WoodsSaxon1b_Op( modelspace, V0, R0, a0);
      }
      else if (opnamesplit[0] == "HOtrap")
      {
         double hw_trap;
         std::istringstream( opnamesplit[1] ) >> hw_trap;
         theop = HOtrap_Op( modelspace, hw_trap);
      }
      else if (opnamesplit[0] == "Schiff")
      {
        double R;
        std::istringstream( opnamesplit[1]) >> R;
        theop =  SchiffOp( modelspace,3,1,R);
      }
      else if (opnamesplit[0] == "VPT" )
      {
        std::vector<double> LECs;
        for (size_t i=1; i<opnamesplit.size(); i++)
        {
          double tmp;
          std::istringstream( opnamesplit[i] ) >> tmp;
          LECs.push_back(tmp);
        }
        if (LECs.size() < 3)
        {
          std::cout << "Uh oh. You didn't specify enough LECs. I'm giving you an empy operator instead." << std::endl;
          return theop;
        }
        theop = TViolatingPotential_Op( modelspace, LECs );
      }
      else if (opnamesplit[0] == "M0nuR") // Radial dependance of GT part of M0nu
      {
        double Eclosure;
        double r12;
        std::istringstream(opnamesplit[1]) >> Eclosure;
        std::istringstream(opnamesplit[2]) >> r12;
        theop = M0nu::GamowTeller_R(modelspace, Eclosure, r12);
      }
      else if (opnamesplit[0] == "DGTR") //Radial depedance of DGT operator
      {
        double r12;
        std::istringstream(opnamesplit[1]) >> r12;
        theop = M0nu::DGT_R(modelspace, r12);
      }
      else if (opnamesplit[0] == "DGTRLocal") //Radial depedance of DGT operator with surface localization
      {
        double r12;
        std::istringstream(opnamesplit[1]) >> r12;
        theop = M0nu::DGT_R_SurfaceLocalization(modelspace, r12);
      }
      else //need to remove from the list
      {
         std::cout << "Unknown operator: " << opname << std::endl;
      }
//      theop =  Operator();
      theop.profiler.timer[opname] += omp_get_wtime() - t_start;
      return theop;
 
 }



 Operator NumberOp(ModelSpace& modelspace, int n, int l, int j2, int tz2)
 {
   Operator NumOp = Operator(modelspace);
   int indx = modelspace.Index1(n,l,j2,tz2);
   NumOp.ZeroBody = 0;
   NumOp.EraseOneBody();
   NumOp.EraseTwoBody();
   NumOp.OneBody(indx,indx) = 1;
   std::cout << "Number op indx = " << indx << " one body = " << NumOp.OneBody << std::endl;
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

 Operator NumberOpRef(ModelSpace& modelspace)
 {
   Operator NumOp = Operator(modelspace);
   for (int h : modelspace.holes)
   {
//     Orbit& oh = modelspace.GetOrbit(h);
     //NumOp.OneBody(h,h) = (oh.j2+1) * oh.occ ;
//     NumOp.OneBody(h,h) = 1 * oh.occ ;
     NumOp.OneBody(h,h) = 1  ;
   }
   return NumOp;
 }

 Operator OneBodyDensity( ModelSpace& modelspace, index_t i, index_t j)
 {
   Operator OBD = Operator(modelspace,0,0,0,2);
   OBD.OneBody(i,j) += 0.5;
   OBD.OneBody(j,i) += 0.5;
   std::cout << "OBD ij = " << i << " " << j << "  One body = " << OBD.OneBody << std::endl;
   return OBD;
 }

 /// Gives \f$ rho_{n,l}(r)\f$ , the norm squared of the harmonic oscillator radial wave function evaluated at r.
 /// This is normalized so that integrating \f$ rho(r) r^2 dr \f$ gives 1.
 /// This means that it is related to the acutal density by a factor of \f$ 4\pi \f$.
 /// This is equivalent to squaring the return value of HO_Radial_psi().
 double HO_density(int n, int l, double hw, double r)
 {
    double v = M_NUCLEON * hw / (HBARC*HBARC);
    double Norm = pow(v/2.,1.5+l) * SQRT2/SQRTPI * pow(2,n+2*l+3) * gsl_sf_fact(n) / gsl_sf_doublefact(2*n + 2*l + 1);
    double L = gsl_sf_laguerre_n(n, l+0.5, v*r*r);
    double rho = Norm * pow(r,2*l) * exp(-v * r*r) * L * L;
    return rho;
 }

 double HO_Radial_psi(int n, int l, double hw, double r)
 {
   double b = sqrt( (HBARC*HBARC) / (hw * M_NUCLEON) );
   double x = r/b;
   double Norm = 2*sqrt( gsl_sf_fact(n) * pow(2,n+l) / SQRTPI / gsl_sf_doublefact(2*n+2*l+1) / pow(b,3.0) );
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
   if ( spwf.size() <1 ) return;
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
    for (auto& a : modelspace.holes)
    {
        Orbit& oa = modelspace.GetOrbit(a);
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

 Operator Charge_Density_Op(ModelSpace& modelspace)
 {
   // in unit of e
   Operator rho(modelspace,0,0,0,2);
   for ( auto a : modelspace.proton_orbits )
   {
      // Orbit & oa = modelspace.GetOrbit(a);
      rho.OneBody(a,a) = 1.; 
   }
   return rho;
 }

/// Lab-frame kinetic energy in the oscillator basis
Operator KineticEnergy_Op(ModelSpace& modelspace)
{
   Operator T(modelspace,0,0,0,2);
   double hw = modelspace.GetHbarOmega();
   for ( auto a : modelspace.all_orbits )
   {
      Orbit & oa = modelspace.GetOrbit(a);
      T.OneBody(a,a) = 0.5 * hw * (2*oa.n + oa.l +3./2); 
      for ( auto b : T.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
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
// \f$ T = p^2/2m  - p^4/(8m^3c^2) + \ldots \f$
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
/// Leading relativistic correction to the kinetic energy
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
   }   
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
/// where \f$ k \f$ labels all quantum numbers other than \f$ n \f$ and a two-body piece
/// \f[  
/// t_{ijkl} = \frac{1}{\hbar\omega} \left\langle ij | (T^{CM}_{12} - T^{rel}_{12}) | kl \right\rangle
/// \f]
 Operator TCM_Op(ModelSpace& modelspace)
 {
   double t_start = omp_get_wtime();
   int E2max = modelspace.GetE2max();
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   if (A<1) std::cout << "DANGER!!! Calling " << __func__ << "  with A=" << A << std::endl;
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
//      std::cout << "== J = " << tbc.J << " == " << std::endl;
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
//            std::cout << "ijkl " << bra.p << " " << bra.q << " " << ket.p << " " << ket.q << "   p1p2 = " << p1p2 << "  E2max is " << E2max << std::endl;
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
              njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
              njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
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

              // factor to account for antisymmetrization. I sure wish this were more transparent...
              //  antisymmetrized matrix elements: V = <ab|V|cd> - <ab|V|dc>
              //  coupled to J, this becomes <abJ|V|cdJ> - (-1)^{jc+jd-J}<ab|V|dcJ>
              //  The phase we obtain for switching c<->d in the transformation from lab to cm/rel is
              //  due to the 9j symbol an the Moshinsky bracket. Together, they amount to -(-1)^{lam+S + jc+jd-J}
              //  so this becomes <abJ|V|cdJ> + (-1)^{lam+S}<abJ|V|cdJ> = (1+ (-1)^{lam+S})<abJ|V|cdJ>.
              //  so exchange gives us a phase factor and we have V = <ab|V|cd> + (-1)**(Lcd+Scd) <ab|V|dc>
              //  case pp,nn -> 1 + (-1)**(L+S)
              //  pnpn -> 1 + 0
              //  pnnp -> 0 + (-1)**(L+S)
              int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
//              double asymm_factor = 1 + modelspace.phase( lam_ab + Sab );
//              if ( bra.op->tz2 != bra.oq->tz2 ) asymm_factor /= 2;
              if ( asymm_factor ==0 ) continue;

//              int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
//              if ( asymm_factor ==0 ) continue;

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
/// \f]
///
/// \f[
/// \delta T_{CM} = \left( \frac{Am}{Zm_p+Nm_n} \right) \frac{1}{2mA} P_{CM}^2
/// \f]
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





/// Center of mass \f$ R^2 \f$, in units of fm \f$^2\f$
/// Returns
/// \f[ 
/// R^{2}_{CM} = \left( \frac{1}{A}\sum_{i}\vec{r}_{i}\right)^2 =
/// \frac{1}{A^2} \left( \sum_{i}r_{i}^{2} + 2\sum_{i<j}\vec{r}_i\cdot\vec{r}_j  \right)
/// \f]
/// evaluated in the oscillator basis.


 Operator R2CM_Op(ModelSpace& modelspace)
 {
   Operator R2cmOp = RSquaredOp(modelspace);

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
   int A = modelspace.GetTargetMass();
   return R2cmOp /(A*A);
 }



/// Intrinsic point proton radius squared
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
/// Calculational details are similar to imsrg_util::Calculate_p1p2().
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

//              std::cout << __func__ << "na la nb lb " << na << " " << la << " " << nb << " " << lb << std::endl;
              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (std::abs(mosh_ab)<1e-8) continue;

              for (int N_cd=std::max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd and N_ab != N_cd) continue;

//              std::cout << __func__ << "nc lc nd ld " << nc << " " << lc << " " << nd << " " << ld << std::endl;
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
   r2.OneBody *= oscillator_b2;
   return r2;
}


 /// One-body part of the \f$ R^2 \f$ operator for protons
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

   std::set<index_t> orbitlist;
   if (option == "neutron" or option == "matter")  orbitlist.insert(modelspace.neutron_orbits.begin(),modelspace.neutron_orbits.end());
   if (option == "proton" or option == "matter")  orbitlist.insert(modelspace.proton_orbits.begin(),modelspace.proton_orbits.end());
   else if (option != "proton" and option != "neutron" and option !="matter") std::cout << "!!! WARNING. " << __func__ << "  BAD OPTION "  << option << std::endl;
 
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

   std::cout << __func__ << " begin" << std::endl;
   int nchan = modelspace.GetNumberTwoBodyChannels();
   if (option!="matter" and option!="proton" and option!="neutron") std::cout << "!!! WARNING. " << __func__ << "  BAD OPTION "  << option << std::endl;
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
//          if ( (bra.op->l > modelspace.GetLmax()) or (bra.oq->l > modelspace.GetLmax()) )
//          {
//            std::cout << "ibra is " << ibra << " => " << bra.p << " " << bra.q << "  => " << bra.op->l << " " <<  bra.oq->l << std::endl;
//          }
         if (Tz==0 and (option=="proton" or option=="neutron")) prefactor = 0.5;
         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            double mat_el = Calculate_r1r2(modelspace,bra,ket,tbc.J) * prefactor; 
            Rp2Op.TwoBody.SetTBME(ch,ibra,iket,mat_el);
            Rp2Op.TwoBody.SetTBME(ch,iket,ibra,mat_el);
         }
      }
   }
   return Rp2Op;
 }


/// Operator whose expectation value gives the proton density at radius R
/// \f$ \hat{rho}(r) = \sum_{ij} \phi_i^{*}(r) \phi_j(r) a^{\dagger}_ia_j \f$
/// Normalized such that \f$ 4\pi \int dr r^2 \rho(r) = Z$.
//Operator ProtonDensityAtR(ModelSpace& modelspace, double R)
Operator DensityAtR(ModelSpace& modelspace, double R, std::string pn)
{
  Operator Rho(modelspace,0,0,0,1);
  double fourpi = 4*PI;
  double hw = modelspace.GetHbarOmega();

  auto pn_list = pn=="proton" ? modelspace.proton_orbits  : modelspace.neutron_orbits;
  for ( auto i : pn_list)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    for ( auto j : Rho.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
       Orbit& oj = modelspace.GetOrbit(j);
       Rho.OneBody(i,j) = HO_Radial_psi(oi.n,oi.l,hw,R) * HO_Radial_psi( oj.n, oj.l, hw, R) / fourpi;
    }
  }
  return Rho;
}


/// Get the point nucleon form factor at a given momentum transfer q.
/// \f$ F(q) = 4\pi \int_0^{\infty} dr r^2 j_0(qr) \rho(r) \f$
/// We absorb the \f$ 4pi \f$ by using radial wave functions which are normalized to 1.
/// The form factor is normalized such that for protons \f$ F(q=0) = Z \f$ and likewise for neutrons.
/// At low q, the form factor behaves as \f$ F(q\approx = 0) = Z(1 - r^2 q^2 / 6) \f$, where $r^2$ is the rms point proton radius.
/// NB, PREX-II measurement is at Q^2=0.0062 GeV^2 => Q = 78.74 MeV =~ 0.4 fm^-1.  PREX-I is at 0.475 fm^-1
Operator FormfactorAtQ(ModelSpace& modelspace, double q, std::string pn)
{
  Operator Fq(modelspace,0,0,0,1);
  double hw = modelspace.GetHbarOmega();
  double bosc = HBARC / sqrt( M_NUCLEON * hw );
  double k = q * bosc; // put q in oscillator units

  int npoints = 200;  // current options for npoints are 0-50, 100, and 200.
  std::vector<double> xvals(npoints,0.);
  std::vector<double> weights(npoints,0.);
  std::vector<double> bessel_j0(npoints,0.);
  for (int ix=0;ix<npoints;ix++)
  {
     xvals[ix]    = GaussLaguerre::gauss_laguerre_points_200[ix][0] ;
     weights[ix]  =   GaussLaguerre::gauss_laguerre_points_200[ix][1] * xvals[ix] * xvals[ix]  ; // includes x^2 from jacobian
     bessel_j0[ix] =   gsl_sf_bessel_j0( k * xvals[ix] ) ;
   }
 
 
  auto pn_list = pn=="proton" ? modelspace.proton_orbits  : modelspace.neutron_orbits;
  for ( auto i : pn_list)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    std::vector<double> phi_i(npoints,0.);

    // compare with HO_Radial_psi: no factor sqrt(1/b^3) because we have a b^3 out front from changing variables from r to r/b.
    double Norm_i = 2*sqrt( gsl_sf_fact(oi.n) * pow(2,oi.n+oi.l) / SQRTPI / gsl_sf_doublefact(2*oi.n+2*oi.l+1) );

    for (int ix=0;ix<npoints;ix++)
    {
       double x = xvals[ix];
       // extra +x in exponent comes from Gauss-Laguerre weighting.
       phi_i[ix] = Norm_i * pow(x,oi.l) * exp(-x*x*0.5 + x) * gsl_sf_laguerre_n(oi.n,oi.l+0.5,x*x) ;
    }
    for ( auto j : Fq.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
       if (i<j) continue;
       Orbit& oj = modelspace.GetOrbit(j);

       double Norm_j = 2*sqrt( gsl_sf_fact(oj.n) * pow(2,oj.n+oj.l) / SQRTPI / gsl_sf_doublefact(2*oj.n+2*oj.l+1) );
       std::vector<double> phi_j(npoints,0.);
       for (int ix=0;ix<npoints;ix++)
       {
          double x = xvals[ix];
          // Notice no +x in the exponent here
          phi_j[ix] = Norm_j * pow(x,oj.l) * exp(-x*x*0.5 ) * gsl_sf_laguerre_n(oj.n,oj.l+0.5,x*x) ;
       }
       double fq_ij = 0;
       for (int ix=0;ix<npoints;ix++)
       {
          fq_ij += weights[ix] *  bessel_j0[ix] * phi_i[ix] * phi_j[ix];
       }
       Fq.OneBody(i,j) = fq_ij;
       Fq.OneBody(j,i) = fq_ij;
    }
  }
  return Fq;
}


// Spin-orbit charge radius correction, originally described
// by Ong et al PRC 82 014320. That expression contained an error,
// as pointed out by Martin Hoferichter.
// The corrected expression is given in Heinz et al PRC 111, 0343411.
// The difference was mu_i - Q_i  ->  mu_i - Q_i/2.
Operator RpSpinOrbitCorrection(ModelSpace& modelspace)
{
  Operator dr_so(modelspace,0,0,0,2);
  double M2 = M_NUCLEON*M_NUCLEON/(HBARC*HBARC);
  int norb = modelspace.GetNumberOrbits();
  double mup = 2.793;
  double mun = -1.913;
  for (int i=0;i<norb;i++)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    double mu_i = oi.tz2<0 ? mup-1./2 : mun;
    int kappa = oi.j2 < 2*oi.l ? oi.l : -(oi.l+1);
    dr_so.OneBody(i,i) = -mu_i/M2*(kappa+1);
  }
  return dr_so;
}

/// Electric monopole operator
/// Returns
/// \f[ r_{e}^2 = \sum_{i} e_{i} r_{i}^2 \f]
///
Operator E0Op(ModelSpace& modelspace)
{
   Operator e0(modelspace);
   e0.EraseZeroBody();
   e0.OneBody.zeros();
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





// Operator FourierBesselCoeff(ModelSpace& modelspace, int nu, double R)
// providing an index list allows us to select e.g. just protons or just neutrons
// Operator FourierBesselCoeff(ModelSpace& modelspace, int nu, double R, std::vector<index_t> index_list)
Operator FourierBesselCoeff(ModelSpace& modelspace, int nu, double R, std::set<index_t> index_list)
{
  Operator a_nu(modelspace,0,0,0,2);

  size_t npoints = 99;  // should be a multiple of 3 so we can use the Simpson 3/8 rule
  std::vector<double> RGRID(npoints);
  std::vector<double> BESSEL(npoints);
  std::vector<double> PSI_p(npoints);
  std::vector<double> INTEGRAND(npoints);
  double dr = R / npoints;
  double Q = nu*PI/R;

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



/// Returns the \f$ T^{2} \f$ operator
 Operator Isospin2_Op(ModelSpace& modelspace)
 {
   Operator T2 = Operator(modelspace,0,0,0,2);
   T2.OneBody.diag().fill(0.75);

   for (int ch=0; ch<T2.nChannels; ++ch)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
     arma::mat& TB = T2.TwoBody.GetMatrix(ch);

     // 2 tz1 * tz2 = |Tz|-1/2
     TB.diag().fill( std::abs( tbc.Tz)-0.5 );

     if (tbc.Tz == 0)
     {
        for (size_t ibra=0;ibra<tbc.GetNumberKets(); ++ibra)
        {
           Ket& bra = tbc.GetKet(ibra);
           Orbit& oa = modelspace.GetOrbit(bra.p);
           Orbit& ob = modelspace.GetOrbit(bra.q);

           size_t abar = modelspace.GetOrbitIndex( oa.n, oa.l, oa.j2, -oa.tz2 );
           size_t bbar = modelspace.GetOrbitIndex( ob.n, ob.l, ob.j2, -ob.tz2 );
           size_t iket = tbc.GetLocalIndex( std::min(abar,bbar), std::max(abar,bbar) ); 
           Ket& ket = tbc.GetKet(iket);
           int flipphase = abar > bbar ? ket.Phase(tbc.J) : 1 ;
           TB(ibra,iket) += flipphase * 1.0;
//           TB(iket,ibra) = TB(ibra,iket);

        }// for ibra
     }// if Tz=-0
   }// for ch
   return T2;
 }
/*
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
             double matel2b = 0;
             if (  oa.j2==oc.j2 and oa.n==oc.n and oa.l==oc.l
               and ob.j2==od.j2 and ob.n==od.n and ob.l==od.l )
             {
               // tz1 tz2 case
               if( oa.tz2 == oc.tz2 and ob.tz2==od.tz2)
               {
                 matel2b -= 0.5;
               }
               // t+ t- case
               if( oa.tz2 == od.tz2 and ob.tz2==oc.tz2)
               {
                 matel2b += 1.0;
               }
             }
             if (  oa.j2==od.j2 and oa.n==od.n and oa.l==od.l
               and ob.j2==oc.j2 and ob.n==oc.n and ob.l==oc.l )
             {

                  int phase = bra.Phase(tbc.J); // this phase includes the fermionic minus sign
                  // tz1 tz2 case
                  if( oa.tz2 == od.tz2 and ob.tz2==oc.tz2)
                  {
                    matel2b -= phase * 0.5;
                  }
                  // t+ t- case
                  if( oa.tz2 == oc.tz2 and ob.tz2==od.tz2)
                  {
                    matel2b += phase * 1.0;
                  }
             }

             TB(ibra,iket) = matel2b;
             TB(iket,ibra) = TB(ibra,iket); // hermitian

           }
        }
     }
   }
   return T2;
 }
*/


/// Returns the \f$ T^{2} \f$ operator
 Operator TzSquared_Op(ModelSpace& modelspace)
 {
   Operator T2 = Operator(modelspace,0,0,0,2);
   T2.OneBody.diag().fill(0.25);

   for (int ch=0; ch<T2.nChannels; ++ch)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
     arma::mat& TB = T2.TwoBody.GetMatrix(ch);

     // 2 tz1 * tz2 = |Tz|-1/2
     TB.diag().fill( std::abs( tbc.Tz)-0.5 );

   }// for ch
   return T2;
 }

  /// Ananople moment operator -- added by A. B.
  Operator AnapoleMoment(ModelSpace& modelspace)
  {
    Operator As(modelspace, 1, 0, 1, 2);
    As.SetAntiHermitian();
    double bL = pow(HBARC * HBARC / M_NUCLEON / modelspace.GetHbarOmega(), 0.5 * 1); // b^L where b=sqrt(hbar/mw)
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit & oi = modelspace.GetOrbit(i);
      double ji = 0.5 * oi.j2;
      double magnetic_moment = oi.tz2 < 0 ? PROTON_SPIN_G/2 : NEUTRON_SPIN_G/2; // These are 5.586 and -3.826 for proton and neutron, respectively. Defined in PhysicalConstants.hh
      for (int j: As.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        if (oi.tz2 != oj.tz2) continue;
        double jj = 0.5 * oj.j2; 
        double r2int = RadialIntegral(oi.n, oi.l, oj.n, oj.l, 1) * bL ;
        double nineJ = AngMom::NineJ(oi.l, 0.5, ji, oj.l, 0.5, jj, 1, 1, 1);
        double hatfactors = sqrt((2 * ji + 1) * (2 * jj + 1) * (2 * oi.l + 1) * (2 * oj.l + 1));
        // Including the factor of sqrt(4pi/3) from the spherical harmonics
        As.OneBody(i, j) = 3 * sqrt(2) * magnetic_moment * r2int * hatfactors * nineJ * AngMom::ThreeJ(oj.l, 1, oi.l, 0, 0, 0);
        As.OneBody(j, i) = -  modelspace.phase((oi.j2-oj.j2)/2) * As.OneBody(i, j); // Operator is imagniary and therefore antisymmetric, reduced ME need a phase under hermitian conjugation
      }
    }
    return As;
  }

  /// Mutipole responses (except dipole) with units fm\f$ ^{rL}\f$ (see PRC97(2018)054306 ) --added by bhu
  Operator MultipoleResponseOp(ModelSpace& modelspace, int rL, int YL, int isospin)
  {
    Operator EL(modelspace, YL,0,YL%2,2);
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*rL); // b^L where b=sqrt(hbar/mw)
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,rL) * bL ;
        EL.OneBody(i,j) = modelspace.phase(jj+YL-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*YL+1)/4./3.1415926) * AngMom::ThreeJ(ji,jj,YL,0.5,-0.5,0) * r2int;
	      if( YL == 0 ) EL.OneBody(i,j) = EL.OneBody(i,j) * sqrt(4.*3.1415926);
        double fac = isospin ? oi.tz2 : 1.0;
        EL.OneBody(i,j) = EL.OneBody(i,j) * fac;
        EL.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * EL.OneBody(i,j);
      }
    }
    return EL;
  }

  /// Isovector dipole with units e * fm , see PRC 97 054306 (2018) --added by bhu
  Operator IVDipoleOp(ModelSpace& modelspace, int rL, int YL)
  {
    Operator EL(modelspace, YL,0,YL%2,2);
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*rL); // b^L where b=sqrt(hbar/mw)
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,rL) * bL ;
        EL.OneBody(i,j) = modelspace.phase(jj+YL-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*YL+1)/4./3.1415926) * AngMom::ThreeJ(ji,jj,YL,0.5,-0.5,0) * r2int;
        EL.OneBody(i,j) = EL.OneBody(i,j) * (0.5 - oi.tz2*0.5 - modelspace.GetTargetZ()*1.0/modelspace.GetTargetMass() );
	      EL.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * EL.OneBody(i,j);
      }
    }
    return EL;
  }

  /// Isoscalar dipole with units e * fm^3  --added by bhu
  Operator ISDipoleOp(ModelSpace& modelspace, int rL, int YL, double Rms)
  {
    Operator EL(modelspace, YL,0,YL%2,2);
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*rL); // b^L where b=sqrt(hbar/mw)
    double bLp = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*1);
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,rL) * bL ;
        double r2intp = RadialIntegral(oi.n,oi.l,oj.n,oj.l,1) * bLp ;
        r2int = r2int - 5.0/3.0*r2intp*Rms*Rms ;
        EL.OneBody(i,j) = modelspace.phase(jj+YL-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*YL+1)/4./3.1415926) * AngMom::ThreeJ(ji,jj,YL,0.5,-0.5,0) * r2int;
	      EL.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * EL.OneBody(i,j);
      }
    }
    return EL;
  }


  /// Schiff Moment = Isoscalar dipole / 10 (where the sum is over proton orbits only)  with units e * fm^3  --added by DK. Ref. PHYSICAL REVIEW C 89, 014335 (2014)
  Operator SchiffOp(ModelSpace& modelspace, int rL, int YL, double Rms)
  {   
    Operator EL(modelspace, YL,0,YL%2,2);
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*rL); // b^L where b=sqrt(hbar/mw)
    double bLp = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*1);
    // int norbits = modelspace.GetNumberOrbits();
    for (int i : modelspace.proton_orbits)
   // for (int i=0; i<norbits; ++i) //DK: modify for protons only
    {     
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,rL) * bL ;
        double r2intp = RadialIntegral(oi.n,oi.l,oj.n,oj.l,1) * bLp ;
        r2int = r2int - 5.0/3.0*r2intp*Rms*Rms ;
        EL.OneBody(i,j) = modelspace.phase(jj+YL-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*YL+1)/4./PI) * AngMom::ThreeJ(ji,jj,YL,0.5,-0.5,0) * r2int;
        EL.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * EL.OneBody(i,j);
      }
    }
    return EL/10;
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
    bool spin = true;
    bool orbit = true;
    if ( pn.find("spin") != std::string::npos ) orbit = false;
    if ( pn.find("orbit") != std::string::npos ) spin = false;
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      if (pn=="proton" and oi.tz2>0) continue;
      if (pn=="neutron" and oi.tz2<0) continue;
      double gl = oi.tz2<0 ? 1.0 : 0.0;
//      double gs = oi.tz2<0 ? 5.586 : -3.826;
      double gs = oi.tz2<0 ? PROTON_SPIN_G : NEUTRON_SPIN_G; // These are 5.586 and -3.826 for proton and neutron, respectively. Defined in PhysicalConstants.hh
      if (not spin) gs = 0;
      if (not orbit) gl = 0;
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
//           printf("< %d %d || E2 || %d %d > = %f   angular: %f \n",nrela,lrela,nrelb,lrelb,ME_rel(ia,ib), angint );
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
//      term2 += gsl_sf_gamma(0.5*(la+lb+L)+sigma+1.5) / (gsl_sf_fact(sigma)*gsl_sf_fact(na-sigma)*gsl_sf_fact(nb-sigma)*gsl_sf_fact(sigma+tau_a-na)*gsl_sf_fact(sigma+tau_b-nb) );
      term2 += tgamma(0.5*(la+lb+L)+sigma+1.5) / (gsl_sf_fact(sigma)*gsl_sf_fact(na-sigma)*gsl_sf_fact(nb-sigma)*gsl_sf_fact(sigma+tau_a-na)*gsl_sf_fact(sigma+tau_b-nb) );
    }
    return term1*term2;
  
  }


 // the nomenclature on the variable k has sort of gotten out of hand here...
 double RadialIntegral_RpowK(int na, int la, int nb, int lb, int k)
 {
   long double I = 0;
   int pmin = (la+lb)/2;
   int pmax = pmin + na + nb;

//   std::vector<double> Ip(pmax+1);
//   Ip[pmin] = TalmiI(pmin,k);
//   for (int p=pmin+1; p<=pmax;++p) Ip[p] = (2*p+1+k)/(2*p+1.) * Ip[p-1];

   if (pmax < 20 )
   {
     for (int p=pmin;p<=pmax;++p)
     {
//        I += TalmiB(na,la,nb,lb,p) * Ip[p];
        I += TalmiB(na,la,nb,lb,p) * TalmiI(p,k);
//        I += AngMom::TalmiB(na,la,nb,lb,p) * TalmiI(p,k);
     }
   }
   else // just do it with quadrature
   {
//     double Iquad = 0;
     double Norm = 2*sqrt( tgamma(na+1)*tgamma(nb+1)/tgamma(na+la+1.5)/tgamma(nb+lb+1.5)) ;
//     double Norm = 2*sqrt( gsl_sf_gamma(na+1)*gsl_sf_gamma(nb+1)/gsl_sf_gamma(na+la+1.5)/gsl_sf_gamma(nb+lb+1.5)) ;
//     int npoints = 50;
//     int npoints = std::min(50,na+nb+(la+lb+1+k)/2+10000);
//     int poly_order = na+nb+(la+lb+1+k)/2+1;
//     int npoints = std::min(50, 2*poly_order-1);
     int npoints = 200;  // current options for npoints are 0-50, 100, and 200.
//     double I1=0;
//     double I2 = 0;
     for (int i=0;i<npoints;i++)
     {
//       double x_i = GaussLaguerre::gauss_laguerre_points[npoints][i][0];
//       double w_i = GaussLaguerre::gauss_laguerre_points[npoints][i][1];
       double x_i = GaussLaguerre::gauss_laguerre_points_200[i][0];
       double w_i = GaussLaguerre::gauss_laguerre_points_200[i][1];
       double f_i = Norm * 0.5* gsl_sf_laguerre_n(na,la+0.5,x_i)  * gsl_sf_laguerre_n(nb,lb+0.5,x_i) * pow(x_i,0.5*(la+lb+1+k) )  ;
       I += w_i * f_i;
     }
   }

   return I;
//   return Iquad;
 }


 double RadialIntegral_Gauss( int na, int la, int nb, int lb, double sigma )
 {
   long double I = 0;
   int pmin = (la+lb)/2;
   int pmax = pmin + na + nb;
   double kappa = 1.0 / (2*sigma*sigma); // Gaussian exp( -x^2 / (2 sigma^2)) =>  exp(-kappa x^2)

//   if (pmax < 20 )
//   {
     for (int p=pmin;p<=pmax;++p)
     {
//        I += TalmiB(na,la,nb,lb,p) * TalmiI(p,k);
        I += TalmiB(na,la,nb,lb,p) * pow(1+kappa, -(p+1.5) )  ; // Talmi integral is analytic
     }
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
//    for (int i=0; i<norbits; ++i)
    for ( auto i : modelspace.proton_orbits )
    {
      Orbit& oi = modelspace.GetOrbit(i);
      int j = modelspace.GetOrbitIndex( oi.n, oi.l, oi.j2, -oi.tz2);
      Fermi.OneBody(i,j) = sqrt(oi.j2+1); // Reduced matrix element
      Fermi.OneBody(j,i) = Fermi.OneBody(i,j); // Hermitian
//      for (int j : Fermi.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
//      {
//        Orbit& oj = modelspace.GetOrbit(j);
//        if (oi.n!=oj.n or oi.tz2 == oj.tz2) continue;
//        Fermi.OneBody(i,j) = sqrt(oi.j2+1.0);  // Reduced matrix element
//      }
    }
    return Fermi;
  }

/// Note that there is a literature convention to include the 1/sqrt(Lambda) factor
/// in the reduced matrix element rather than in the expression involving the sum
/// over one-body densities (see footnote on pg 165 of Suhonen).
/// I do not follow this convention, and instead produce the reduced matrix element
///  \f[ \langle f \| \sigma \tau_{\pm} \| i \rangle \f]
/// Note that this matrix element does not include the axial coupling gA.
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
 // We use equation (A5) of Caurier et al Rev Mod Phys 2005
 // S
 Operator QdotQ_Op(ModelSpace& modelspace)
 {
    
//   // temporarily store <i||Q||j> in the one body part.
   Operator QdotQ_op(modelspace,0,0,0,2);
   Operator E2op = ElectricMultipoleOp(modelspace,2) + NeutronElectricMultipoleOp(modelspace,2);
   auto& Qmat = E2op.OneBody;

   std::cout << "Oops! This operator is still under construction! " << __FILE__ << "  line " << __LINE__ << std::endl;

   // The one-body piece should be
   // < i | Q*Q | j > = sqrt(2*2+1)/(2*ji+1) * sum_k <k||Q||i><k||Q||j>
//   QdotQ_op.OneBody =  Qmat*Qmat.t() / sqrt(2);
//   for (auto i : modelspace.all_orbits)
//   {
//      Orbit& oi = modelspace.GetOrbit(i);
//      QdotQ_op.OneBody.row(i) /= oi.j2+1;
//   }

//   // We subtract off the one-body piece acting in the two-body space
//   Embed1BodyIn2Body( QdotQ_op, 2);
//   QdotQ_op.TwoBody *= -1;

   int nchan = modelspace.GetNumberTwoBodyChannels();


   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      int J = tbc.J;
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int a = bra.p;
         int b = bra.q;
         Orbit & oa = modelspace.GetOrbit(a);
         Orbit & ob = modelspace.GetOrbit(b);
         double ja = oa.j2*0.5;
         double jb = ob.j2*0.5;

         for (int iket=ibra;iket<nkets;++iket)
         {
            
            Ket & ket = tbc.GetKet(iket);
            int c = ket.p;
            int d = ket.q;
            Orbit & oc = modelspace.GetOrbit(c);
            Orbit & od = modelspace.GetOrbit(d);
            double jc = oc.j2*0.5;
            double jd = od.j2*0.5 ;

            // Suhonen (8.56)
            double Aabcd = modelspace.phase(ja+jb+J) * modelspace.GetSixJ(ja,jb,J,jd,jc,2) * Qmat(c,a) * Qmat(b,d);
            double Aabdc = modelspace.phase(ja+jb+J) * modelspace.GetSixJ(ja,jb,J,jc,jd,2) * Qmat(d,a) * Qmat(b,c);
            // Suhonen (8.55),(8.57),(8.58)
            double QdQ;
            QdQ = Aabcd - modelspace.phase(jc+jd+J)*Aabdc;  // pppp or nnnn
            if (a==b) QdQ /= sqrt(2.0);
            if (c==d) QdQ /= sqrt(2.0);


//            QdotQ_op.TwoBody.SetTBME(ch,ibra,iket,QdQ);
            QdotQ_op.TwoBody.AddToTBME(ch,ibra,iket,QdQ);
         }

      }
   }

   return QdotQ_op;
 }



 // Quadrupole-quadrupole interaction
 // This is a pure two-body operator.
 Operator VQQ_Op(ModelSpace& modelspace)
 {
   Operator VQQ = Operator(modelspace,0,0,0,2);
   // temporarily store <i||Q||j> in a one-body operator.
   Operator Q1b = Operator(modelspace,2,0,0,1);
   Q1b.OneBody += ElectricMultipoleOp(modelspace,2).OneBody;
   Q1b.OneBody += NeutronElectricMultipoleOp(modelspace,2).OneBody;
   auto& Qmat = Q1b.OneBody;

   int nchan = modelspace.GetNumberTwoBodyChannels();

   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      int J = tbc.J;
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int a = bra.p;
         int b = bra.q;
         Orbit & oa = modelspace.GetOrbit(a);
         Orbit & ob = modelspace.GetOrbit(b);
         double ja = oa.j2*0.5;
         double jb = ob.j2*0.5;

         for (int iket=ibra;iket<nkets;++iket)
         {
            
            Ket & ket = tbc.GetKet(iket);
            int c = ket.p;
            int d = ket.q;
            Orbit & oc = modelspace.GetOrbit(c);
            Orbit & od = modelspace.GetOrbit(d);
            double jc = oc.j2*0.5;
            double jd = od.j2*0.5 ;

            // Suhonen (8.56)
            double Aabcd = modelspace.phase(ja+jb+J) * modelspace.GetSixJ(ja,jb,J,jd,jc,2) * Qmat(c,a) * Qmat(b,d);
            double Aabdc = modelspace.phase(ja+jb+J) * modelspace.GetSixJ(ja,jb,J,jc,jd,2) * Qmat(d,a) * Qmat(b,c);
            // Suhonen (8.55),(8.57),(8.58)
            double QdQ;
            // I think that the pppp case can actually be applied to the pn case too...
//            if ( tbc.Tz==0 ) // proton-neutron channel
//            {
//              if ( oa.tz2 == oc.tz2)   QdQ = Aabcd;  // pnpn, Suhonen (8.57)
//              else                     QdQ = -modelspace.phase(jc+jd+J) * Aabdc;  // pnnp
//            }
//            else  // like-nucleon channel,  Suhonen (8.55)
//            {
               QdQ = Aabcd - modelspace.phase(jc+jd+J)*Aabdc;  // pppp or nnnn
               if (a==b) QdQ /= sqrt(2.0);
               if (c==d) QdQ /= sqrt(2.0);
//            }

            VQQ.TwoBody.SetTBME(ch,ibra,iket,QdQ);
         }
      }
   }

   return VQQ;
 }



 Operator Dagger_Op( ModelSpace& modelspace, index_t Q )
 {
   Operator dag(modelspace);
   dag.SetNumberLegs(3);
   dag.SetQSpaceOrbit(Q);
   dag.OneBody(Q,0)= 1.0;
//   dag.OneBody(Q,Q)= 1.0;
   dag.SetNonHermitian();
   std::cout << "Making a dagger operator. I think Q = " << Q << std::endl;
   return dag;
 }

 Operator DaggerAlln_Op( ModelSpace& modelspace, index_t Q )
 {
   Orbit& oQ = modelspace.GetOrbit(Q);
   Operator dag(modelspace);
   dag.SetNumberLegs(3);
   dag.SetQSpaceOrbit(Q);
   dag.SetNonHermitian();
   for ( auto nQ : modelspace.OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
   {
     //dag.OneBody(nQ,Q) = 1.0;
     dag.OneBody(nQ,0) = 1.0;
   }
   return dag;
 }

 Operator VCentralCoulomb_Op( ModelSpace& modelspace, int lmax ) // default lmax=99999
 {
//   std::cout << "Making VCentralCoulomb_Op. lmax is " << lmax  << std::endl;
   Operator VCoul(modelspace, 0,0,0,2);
   double oscillator_b = sqrt(HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
//   double alpha_FS = 1.0 / 137.035999;
  
// First, the one-body piece <a|1/r|b>
//   int norb = modelspace.GetNumberOrbits();
//   for (int a=0; a<norb; a++)
   for (auto a : modelspace.all_orbits)
   {
     Orbit& oa = modelspace.GetOrbit(a);
     if (oa.tz2>0) continue; // protons only
     if (oa.l>lmax) continue;
     for (auto b : VCoul.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
     {
       if (b<a) continue;
       Orbit& ob = modelspace.GetOrbit(b);
       double rad_int =  RadialIntegral_RpowK(oa.n, oa.l, ob.n, ob.l, -1) ;  
       VCoul.OneBody(a,b) = rad_int;
       VCoul.OneBody(b,a) = rad_int;
     }
   }
   VCoul.OneBody *= ALPHA_FS * HBARC / oscillator_b; // convert from oscillator units to fermi. ALPHA_FS ~ 1/137 is the fine structure constant
//   VCoul.OneBody *= alpha_FS * HBARC / oscillator_b; // convert from oscillator units to fermi
//   std::cout << "Oscillator b = " << oscillator_b << std::endl;
//   std::cout << "One body part done. it looks like" << std::endl << VCoul.OneBody << std::endl;
   return VCoul;
 }




 Operator VCoulomb_Op( ModelSpace& modelspace, int lmax ) //default lmax=99999
 {
//   std::cout << "Making VCoulomb_Op" << std::endl;
   double t_start = omp_get_wtime();
   Operator VCoul(modelspace, 0,0,0,2);
   double oscillator_b = sqrt(HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());

   int nmax = modelspace.GetEmax();
   int lrelmax = modelspace.GetEmax();
   double tt_start = omp_get_wtime();
   std::unordered_map<size_t,double> RadialIntegrals;
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
   VCoul.profiler.timer["ComputeCoulombIntegrals"] += omp_get_wtime() - tt_start;

//   std::cout << "now the big loop... lrelmax = " << lrelmax << "  size of RadInt = " << RadialIntegrals.size()  << std::endl;

// Now the (antisymmetrized) two-body piece <ab| 1/r_rel |cd>
   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
//   std::cout << "Done Precalculating Moshinsky." << std::endl;
   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;
//   #pragma omp parallel for schedule(dynamic,1)  // It would appear that something's not thread-safe in this routine...
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
         
                       int N_cd = N_ab;
                       int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                       if (n_cd < 0) continue;
         
                       double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                       if (std::abs(mosh_cd)<1e-8) continue;

                       double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;

                       size_t hash      = n_ab*(nmax+1)*(lrelmax+1) + n_cd*(lrelmax+1) + lam_ab;   
                       if ( RadialIntegrals.find(hash) == RadialIntegrals.end() )
                       {
                         std::cout << "AAAHHH!!!  trying to access radial integral for " << n_ab << " " << n_cd << " " << lam_ab << "    and it's not there!!!!" << std::endl;
                       }
                       double rad_int =  RadialIntegrals.at(hash) ;  
 
                       rinv += prefactor * rad_int; 
    
                    } // lam_ab
                  } // Lam_ab
                } // N_ab
         
              } // Sab
            } // Lab

            // In Moshinsky's convention, r_rel = (r1-r2)/sqrt(2).  We want 1/|r1-r2| = 1/sqrt(2) * 1/r_rel
            rinv *=  1 / sqrt(2*(1.0+bra.delta_pq())*(1.0+ket.delta_pq())); // normalize and account for sqrt(2) Moshinsky convention
            VCoul.TwoBody.SetTBME(ch,ibra,iket,rinv);
            VCoul.TwoBody.SetTBME(ch,iket,ibra,rinv);
                         
         }
      }
   }
   VCoul *= ALPHA_FS * HBARC / oscillator_b;  // convert to MeV.  V = e^2/r = alpha*hc / r

//   std::cout << "All done with VCoul." << std::endl;
//   std::cout << "0s pp,pn,nn: " << VCoul.TwoBody.GetTBME_J(0,0, 0,0, 0,0) << " " << VCoul.TwoBody.GetTBME_J(0,0, 1,0, 1,0) << " " << VCoul.TwoBody.GetTBME_J(0,0, 1,1, 1,1) << std::endl;
   VCoul.profiler.timer[__func__] += omp_get_wtime() - t_start;
   return VCoul ;

 }




 Operator AxialCharge_Op( ModelSpace& modelspace )
 {
   Operator AxCh(modelspace, 0,1,1,2 );
   double oscillator_b = HBARC*HBARC/sqrt(M_NUCLEON * modelspace.GetHbarOmega());
   for (auto a : modelspace.all_orbits )
   {
     Orbit& oa = modelspace.GetOrbit(a);
     double prefactor = modelspace.phase((oa.j2+1)/2) * sqrt( (oa.j2+1.0)*6.0/4/M_PI )  ;
     std::cout << " a = " << a << " size of one body channels = " << AxCh.OneBodyChannels.size() << " pre: " << (oa.j2+1)/2 << "  " << sqrt((oa.j2+1.0)*6.0/4/M_PI) << "  " << prefactor << std::endl;
     std::cout << " OBC.at(a) : ";
     for (auto b : AxCh.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) std::cout << b << "  ";
     std::cout << std::endl;
     for (auto b : AxCh.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
     {
       Orbit& ob = modelspace.GetOrbit(b);
       double sixj = modelspace.GetSixJ(oa.l,0.5,0.5*oa.j2,0.5,ob.l,1);
       double threej = AngMom::ThreeJ(oa.l, 1, ob.l, 0,0,0);
       double radialint = imsrg_util::RadialIntegral_RpowK( oa.n, oa.l, ob.n, ob.l, 1 );
       std::cout << "prefactor " << prefactor << "   sixj " << sixj << "   threej " << threej << "   radialint " << radialint << "  b " << oscillator_b << std::endl;
       AxCh.OneBody(a,b) = prefactor * sixj * threej * radialint / oscillator_b;
     }
   }

   std::cout << "AxialCharge 1b looks like"<< std::endl << AxCh.OneBody << std::endl <<std::endl;

  return AxCh;

 }

 Operator TViolatingPotential_Op( ModelSpace& modelspace, std::vector<double> LECs )
 {
  //alpha0,alpha1,alpha2 are the inputs which are gg_0,2mpi/(8piMnucleon)
  double alpha0 = LECs[0];
  double alpha1 = LECs[1];
  double alpha2 = LECs[2];
  Operator VPT(modelspace, 0,0,1,2 ); // J=0, dTz=0, parity=odd, particle-rank=2

  // Isospin factors converted to proton-neutron basis
    // ATT'Tz = alpha0[4T-3]+2*alpha1*Tz + 2*alpha2*delta_T,1*[delta_|Tz|,1 -2delta_|Tz|,0
    // BTT' = 2 alpha1 delta_T+T',1
    // A111 = alpha0 + 2 alpha1 + 2 alpha2
    // A110 = alpha0 - 4 alpha2
    // A000 = -3 alpha0 -4 alpha2
    // B10 = B01 = 2 alpha1
    // |np> = (|1,0> + |0,0>)/sqrt(2),  |pn> = (|1,0> - |0,0>)/sqrt(2)
    // Apppp = A111  = alpha0 + 2 alpha1 + 2 alpha2
    // Annnn = A11-1 = alpha0 - 2 alpha1 + 2 alpha2
    // Apnpn = Anpnp = 1/2 (A110 + A000) = -alpha0 - 4 alpha2
    // Anppn = Apnnp = 1/2 (A110 - A000) = 2 alpha0
    // Bnpnp = -Bpnpn = 1/2 (B10 + B01) = 2 alpha1
    // Bnppn = -Bpnnp = 1/2(B01 - B10) = 0
  double Apppp = alpha0 + 2*alpha1 + 2*alpha2;
  double Annnn = alpha0 - 2*alpha1 + 2*alpha2;
  double Apnpn = -alpha0 - 4*alpha2;
  double Anpnp = Apnpn;
  double Anppn = 2 * alpha0;
  double Apnnp = Anppn;
  double Bnpnp = 2 * alpha1;
  double Bpnpn = -Bnpnp;

 //  auto isoA = [alpha0,alpha1,alpha2] (int T,int Tp,int Tz) { return alpha0*(4*T-3) + 2*alpha1*Tz + 2*alpha2*T*( 3*std::abs(Tz) - 2) * (T==Tp) ; } ;
  auto isoA = [alpha0,alpha1,alpha2] (int T,int Tp,int Tz) { return (alpha0*(4*T-3) + 2*alpha1*Tz + 2*alpha2*(T==1)*((std::abs(Tz)==1) - 2*(Tz==0))) * (T==Tp) ; } ;
  auto isoB = [alpha1] (int T,int Tp) { return 2*alpha1 * ((T+Tp)==1) ; };

  double hw = modelspace.GetHbarOmega();
  double bosc = HBARC/sqrt( hw * M_NUCLEON);
  double mpi = (2*PhysConst::M_PION_CHARGED + PhysConst::M_PION_NEUTRAL)/3;
  int emax = modelspace.Emax;

   std::map< std::array<int,4>,double> IntegralLookup;
   #pragma omp parallel for schedule(dynamic, 1)
   // Do the radial integral and cache it
   for (int n=0; n <= emax; n++)
   {
     for (int np=0; np <= emax; np++)
     {
       for ( int l=0; 2*n+l <= 2*emax; l++)
       {
         for ( int lp=0; 2*np+lp <= 2*emax; lp++)
         {
           double the_integral=0;
           double dr = 0.01;
           double rmax = 5*bosc;
           for ( double r=0.001; r<=rmax; r+=dr )
           {
             double psi_nl   = HO_Radial_psi(n,l,hw,r);
             double psi_nplp = HO_Radial_psi(np,lp,hw,r);
             double x = sqrt(2) * r;
             double fprime = -exp(-mpi * x/HBARC) / x * ( 1 + 1/(mpi*x/HBARC));
             the_integral += r*r * psi_nl * psi_nplp * fprime * dr;
           }
           IntegralLookup[{n,l,np,lp}] = the_integral;
         }
       }
     }
   }


   double sa,sb,sc,sd;
   sa=0.5; sb=0.5; sc=0.5; sd=0.5; // just for clarity...
   std::vector<int> ch_bra_list, ch_ket_list;
   for (auto& itmat : VPT.TwoBody.MatEl )
   {
     ch_bra_list.push_back(itmat.first[0]);
     ch_ket_list.push_back(itmat.first[1]);
   }
   int nbraket = ch_bra_list.size();
   #pragma omp parallel for schedule(dynamic, 1)
   for (int ibraket=0; ibraket<nbraket; ibraket++ )
   {
     int ch_bra = ch_bra_list[ibraket];
     int ch_ket = ch_ket_list[ibraket];
     TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(ch_bra);
     TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(ch_ket);
     int J = tbc_bra.J;
     int Tz = tbc_bra.Tz;
     int nbras = tbc_bra.GetNumberKets();
     int nkets = tbc_ket.GetNumberKets();
     for (int ibra=0; ibra<nbras; ibra++)
     {
       Ket& bra = tbc_bra.GetKet(ibra);
       int la = bra.op->l;
       int lb = bra.oq->l;
       int na = bra.op->n;
       int nb = bra.oq->n;
       double ja = 0.5*bra.op->j2;
       double jb = 0.5*bra.oq->j2;
       int ebra = 2*na+la + 2*nb+lb;

       for (int iket=0; iket<nkets; iket++)
       {
         Ket& ket = tbc_ket.GetKet(iket);
         int lc = ket.op->l;
         int ld = ket.oq->l;
         int nc = ket.op->n;
         int nd = ket.oq->n;
         double jc = 0.5*ket.op->j2;
         double jd = 0.5*ket.oq->j2;
         int eket = 2*nc+lc + 2*nd+ld;

         std::string pn_case = "";
         for ( auto tz : { bra.op->tz2, bra.oq->tz2 , ket.op->tz2, ket.oq->tz2} )
         {
            if (tz==-1) pn_case += "p";
            else pn_case += "n";
         }

         double vpt = 0;
         for (int S=0; S<=1; S++)
         {
           for (int Sp=0; Sp<=1; Sp++)
           {
             if (S+Sp==0)
                continue;
             int Lmin = std::max( J-S, 0);
             int Lmax = J+S;
             int Lpmin = std::max( J-Sp, 0);
             int Lpmax = J+Sp;
             for (int L=Lmin; L<=Lmax; L++)
             {
               double ninej1 = AngMom::NineJ( la,sa,ja,  lb,sb,jb, L,S,J);
               //std::cout << "la=" << la << "sa=" << sa << "ja=" << ja << "lb=" << lb << "sb="<< sb << "jb=" << jb << "L=" << L << "S=" << S << "J=" << J << "ninej1 = " << ninej1 << std::endl;
               for (int Lp=Lpmin; Lp<=Lpmax; Lp++)
               {
                 if ( std::abs(L-Lp)>1 ) continue;
                 double sixj1 =  AngMom::SixJ( L, S, J, Sp, Lp, 1);
                 //std::cout << "L=" << L << "S=" << S << "J=" << J << "Sp=" << Sp << "Lp=" << Lp << "sixj1 = " << sixj1 << std::endl;
                 //std::cout << "sixj1 = " << sixj1 << std::endl;
                 double ninej2 = AngMom::NineJ( lc,sc,jc, ld,sd,jd, Lp,Sp,J);
                 //std::cout << "lc=" << lc << "sc=" << sc << "jc=" << jc << "ld=" << ld << "sd="<< sd << "jd=" << jd << "Lp=" << Lp << "Sp=" << Sp << "J=" << J << "ninej2 = " << ninej2 << std::endl;
                 //std::cout << "ninej2 = " << ninej2 << std::endl;
                 // Now sum over relative/CM coordinates
                 for (int Ncm=0; 2*Ncm<=std::min(ebra,eket); Ncm++ )
                 {
                   for (int Lcm=0; 2*Ncm+Lcm <= std::min(ebra,eket); Lcm++)
                   {
                     for (int nrel=0; 2*nrel <= ebra-(2*Ncm+Lcm); nrel++)
                     {
                       int lrel = ebra - (2*Ncm+Lcm+2*nrel);
                       for (int nprel=0; 2*nprel <= eket-(2*Ncm+Lcm); nprel++)
                       {
                         int lprel = eket - (2*Ncm+Lcm+2*nprel);
                         double AplusB = 0;
                         for (int Tp=std::abs(Tz); Tp<=1; Tp++)
                         {
                            int asymmfactor = 1 - AngMom::phase(lprel+Sp+Tp);
                            for (int T=std::abs(Tz); T<=1; T++)
                            {
                               double iso_to_pn_factor = 1;  // This account for converting an isospin expression to pn formalism
                               if (pn_case == "npnp") iso_to_pn_factor = 0.5;
                               if (pn_case == "pnpn") iso_to_pn_factor = 0.5 * AngMom::phase( T+Tp ) ;
                               if (pn_case == "nppn") iso_to_pn_factor = 0.5 * AngMom::phase( Tp+1  ) ;
                               if (pn_case == "pnnp") iso_to_pn_factor = 0.5 * AngMom::phase( T+1 ) ;
                               
                               if ( (S+Sp)==1)           AplusB += iso_to_pn_factor * asymmfactor * isoA(T,Tp,Tz);
                               if ( (S==1) and (Sp==1) ) AplusB += iso_to_pn_factor * asymmfactor * isoB(T,Tp) * sqrt(2); //sqrt(2) missing here
                               //std::cout << "asymmfactor = " << asymmfactor << std::endl; 
                           }
                         }
                         //if ((lprel + Sp)%2==0) continue; // factor [1 - (-1)^(lprel+Sp) ]
                         //std::cout << "AplusB = " << AplusB << std::endl;
                         double threej = AngMom::ThreeJ( lrel,lprel,1, 0,0,0);
                         //std::cout << "threej = " << threej << std::endl;
                         //std::cout << "lrel=" << lrel << "lprel=" << lprel << "threej=" << threej << std::endl;
                         double sixj2 = AngMom::SixJ( lrel,L,Lcm,  Lp,lprel,1);
                         //std::cout << "sixj2 = " << sixj2 << std::endl;
                         //std::cout << "lrel=" << lrel << "L=" << L << "Lcm=" << Lcm << "Lp=" << Lp << "lprel=" << lprel << "sixj2 = " << sixj2 << std::endl;
                         double mosh1 = AngMom::Moshinsky(Ncm,Lcm,nrel,lrel,na,la,nb,lb,L);
                         //std::cout << "mosh1 = " << mosh1 << std::endl;
                         //std::cout << "Ncm=" << Ncm << "Lcm=" << Lcm << "nrel=" << nrel << "lrel=" << lrel << "na="<< na << "la=" << la << "nb=" << nb << "lb=" << lb << "L=" << L << "mosh1 = " << mosh1 << std::endl;
                         double mosh2 = AngMom::Moshinsky(Ncm,Lcm,nprel,lprel,nc,lc,nd,ld,Lp);
                         //std::cout << "mosh2 = " << mosh2 << std::endl;
                         //std::cout << "Ncm=" << Ncm << "Lcm=" << Lcm << "nprel=" << nprel << "lprel=" << lprel << "nc="<< nc << "lc=" << lc << "nd=" << nd << "ld=" << ld << "Lp=" << Lp << "mosh2 = " << mosh2 << std::endl;
                         double integral = IntegralLookup[{nrel,lrel,nprel,lprel}];
                         //std::cout << "integral = " << integral << std::endl;
                         vpt += sqrt((2*S+1)*(2*Sp+1)) * AplusB * (2*L+1) * (2*Lp+1) * sixj1 * ninej1 * ninej2
                               * AngMom::phase(Lcm) * sqrt((2*lrel+1)*(2*lprel+1)) * threej * sixj2 * mosh1*mosh2 * integral; // *2 factor is from the integral, was already included bt Ragnar  
                         //std::cout << "sixj1 = " << sixj1 << std::endl;
                         //std::cout << "ninej1 = " << ninej1 << std::endl;
                         //std::cout << "ninej2 = " << ninej2 << std::endl;
                         //std::cout << "AplusB = " << AplusB << std::endl;
                         //std::cout << "threej = " << threej << std::endl;
                         //std::cout << "sixj2 = " << sixj2 << std::endl;
                         //std::cout << "mosh1 = " << mosh1 << std::endl;
                         //std::cout << "mosh2 = " << mosh2 << std::endl;
                         //std::cout << "integral = " << integral << std::endl;
                      }//nprel
                     }//nrel
                   }//Lcm
                 }//Ncm
               }//Lp
             }//L
           }//Sp
         }//S
         vpt *= AngMom::phase(J) * 2*sqrt(3) * sqrt((2*ja+1)*(2*jb+1)*(2*jc+1)*(2*jd+1)) ;// * (2*J+1); //(2*J+1) should not be here?
         vpt *= sqrt(2*J+1.); //Because the operator violates parity, the code expects the operator to be reduced
         if (bra.p==bra.q) vpt /= SQRT2;
         if (ket.p==ket.q) vpt /= SQRT2;

         //std::cout << "Setting  channels " << ch_bra << " " << ch_ket << " " << ibra << " " << iket << "  = " << vpt << std::endl;

         //std::cout << "The Matrix: " << std::endl << VPT.TwoBody.GetMatrix(ch_bra,ch_ket) << std::endl;
         VPT.TwoBody.SetTBME(ch_bra,ch_ket,ibra,iket,vpt);
       }//iket
     }//ibra
     
   }//itmat

   return VPT;

 }//TViolatingPotential_Op





 Operator WoodsSaxon1b_Op( ModelSpace& modelspace, double V0, double R0, double a0)
 {
    Operator VWS_op( modelspace, 0,0,0,2 ); // This should just be a 1b operator...

    double hw = modelspace.GetHbarOmega();
    int ngrid = 500;
    double rmax = 3 * R0;
    double dr = rmax/ngrid;
    std::vector<double> rgrid(ngrid,0.0);
    std::vector<double> Vgrid(ngrid,0.0);
    std::vector<double> psi_a(ngrid,0.0);
    std::vector<double> psi_b(ngrid,0.0);
    for (int igrid=0;igrid<ngrid;igrid++)  rgrid[igrid] = dr * igrid;
    for (int igrid=0;igrid<ngrid;igrid++)  Vgrid[igrid] = V0 / (1.0 + exp( (rgrid[igrid]-R0)/a0 ) );

    for (auto a : modelspace.all_orbits)
    {
      Orbit& oa = modelspace.GetOrbit(a);
      for (int igrid=0;igrid<ngrid;igrid++) psi_a[igrid] =  HO_Radial_psi( oa.n, oa.l, hw, rgrid[igrid] );
      for (auto b: VWS_op.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
        if (b<a) continue;
        Orbit& ob = modelspace.GetOrbit(b);
        double Vab = 0;
        for (int igrid=0;igrid<ngrid;igrid++)
        {
           double psi_b_r =  HO_Radial_psi( ob.n, ob.l, hw, rgrid[igrid] );
           Vab += psi_a[igrid] * psi_b_r * Vgrid[igrid] * rgrid[igrid] * rgrid[igrid] * dr;
        }
        VWS_op.OneBody(a,b) = Vab;
        VWS_op.OneBody(b,a) = Vab;
      }
    }
    return VWS_op;
 }



 Operator HOtrap_Op( ModelSpace& modelspace, double hw_trap)
 {
    Operator VHO_op( modelspace, 0,0,0,2 ); // This should just be a 1b operator...

    double hw = modelspace.GetHbarOmega();
    double bosc2 = HBARC*HBARC/( M_NUCLEON * hw_trap );
    int ngrid = 500;
    double rmax = 4 * sqrt(bosc2);
    double dr = rmax/ngrid;
    std::vector<double> rgrid(ngrid,0.0);
    std::vector<double> Vgrid(ngrid,0.0);
    std::vector<double> psi_a(ngrid,0.0);
    std::vector<double> psi_b(ngrid,0.0);
    for (int igrid=0;igrid<ngrid;igrid++)  rgrid[igrid] = dr * igrid;
    for (int igrid=0;igrid<ngrid;igrid++)  Vgrid[igrid] = 0.5 * rgrid[igrid]*rgrid[igrid]/bosc2  * hw_trap;

    for (auto a : modelspace.all_orbits)
    {
      Orbit& oa = modelspace.GetOrbit(a);
      for (int igrid=0;igrid<ngrid;igrid++) psi_a[igrid] =  HO_Radial_psi( oa.n, oa.l, hw, rgrid[igrid] );
      for (auto b: VHO_op.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
        if (b<a) continue;
        Orbit& ob = modelspace.GetOrbit(b);
        double Vab = 0;
        for (int igrid=0;igrid<ngrid;igrid++)
        {
           double psi_b_r =  HO_Radial_psi( ob.n, ob.l, hw, rgrid[igrid] );
           Vab += psi_a[igrid] * psi_b_r * Vgrid[igrid] * rgrid[igrid] * rgrid[igrid] * dr;
        }
        VHO_op.OneBody(a,b) = Vab;
        VHO_op.OneBody(b,a) = Vab;
      }
    }
    return VHO_op;
 }


 // Second-order estimate of the spectroscopic factor
 // for removal of a nuleon in orbit p.
 //                          o               o
 //             o  /~~~/\    p\             p|
 //   o         p\/|  (  )     \/|~~~/\      |~~~/\
 //  p|    +     /\|a b\/i  +  /\|a b\/i  +  |i a\/j
 //   |        p/  ~~~~~     p/   ~~~~      p|~~~~
 // __|__    __/_____      __/______       __|_____
 //
 // The formula that I derived, and have not checked exhaustively, is
 // < A | adagger | A-1> = 1 + 1/2 sum_abi |V_piab|^2 * ( 1/(ep+ei-ea-eb)+1/(ep))/(2ep+ei-ea-eb)
 //                          - 1/2 sum_aij |V_apij|^2 * 1/( ep*(ei+ej-ea-eb) )
 //
 //          when doing the j coupled version, the Vpiab term gets a 2J+1/2ji+1  and the Vapij gets 2J+1/2ja+1
 //          The spectroscopic factor is then the reduced matrix element squared which gives
 //          SF = 2jp+1 |amplitude|^2
 //
 double MBPT2_SpectroscopicFactor( Operator H, index_t p)
 {
   Orbit& op = H.modelspace->GetOrbit(p);
   double amplitude = op.occ; // leading order amplitude
   double ep = H.OneBody(p,p);
   double ef = 0;
   std::cout << "SPEs: " << std::endl;
   for (auto p : H.modelspace->all_orbits) std::cout << p << " :  " << H.OneBody(p,p) << std::endl;

   for (auto a : H.modelspace->particles )
   {
     Orbit& oa = H.modelspace->GetOrbit(a);
     double ea = H.OneBody(a,a);
     for (auto i : H.modelspace->holes )
     {
       Orbit& oi = H.modelspace->GetOrbit(i);
       double ei = H.OneBody(i,i);
       for (auto b : H.modelspace->particles )
       {
         Orbit& ob = H.modelspace->GetOrbit(b);
         double eb = H.OneBody(b,b);
         int Jmin = std::max( std::abs(oa.j2-ob.j2), std::abs(oi.j2-op.j2) ) /2;
         int Jmax = std::min( oa.j2+ob.j2 ,  oi.j2+op.j2 )/2;
         double Jterm = 0;
         for (int J=Jmin; J<=Jmax; J++)
         {
           // GetTBME_J returns an unnormalized matrix element, which is what I want.
           double Vpiab = H.TwoBody.GetTBME_J(J,p,i,a,b);
           Jterm += 0.5 * (2*J+1.)/(oi.j2+1) * Vpiab*Vpiab * ( 1.0/(ep+ei-ea-eb-ef) + 1.0/ep)/(2*ep+ei-ea-eb);
           amplitude += 0.5 * (2*J+1.)/(oi.j2+1) * Vpiab*Vpiab * ( 1.0/(2*ep+ei-ea-eb-ef) + 1.0/ep)/(ep+ei-ea-eb);
         }
         std::cout << "i,a,b, " << i << " " << a << " " <<b << " spes  " << ei << " " << ea << " " << eb << "  denominator  " << ( 1.0/(2*ep+ei-ea-eb-ef) + 1.0/ep)/(ep+ei-ea-eb) << "  contribution  " << Jterm << std::endl;
       }
       for (auto j : H.modelspace->holes )
       {
         Orbit& oj = H.modelspace->GetOrbit(j);
         double ej = H.OneBody(j,j);
         int Jmin = std::max( std::abs(oa.j2-op.j2), std::abs(oi.j2-oj.j2) ) /2;
         int Jmax = std::min( oa.j2+op.j2 ,  oi.j2+oj.j2 )/2;
         double Jterm = 0;
         for (int J=Jmin; J<=Jmax; J++)
         {
           double Vpaij = H.TwoBody.GetTBME_J(J,p,a,i,j);
           Jterm -= 0.5 * (2*J+1.)/(oa.j2+1) * Vpaij*Vpaij * 1.0/((ep-ef)*(ei+ej-ea-ef));
           amplitude -= 0.5 * (2*J+1.)/(oa.j2+1) * Vpaij*Vpaij * 1.0/((ep-ef)*(ei+ej-ea-ef));
         }
         std::cout << "i,j,a " << i << " " << j << " " << a << " spes  " << ei << " " << ej << " " << ea << "  denominator  " << 1.0/((ep-ef)*(ei+ej-ea-ef)) << "   contribution  " << Jterm << std::endl;
       }
     }	 
   }
   std::cout << "Amplitude is " << amplitude << std::endl;
   double SF = (op.j2+1.0) * amplitude * amplitude;
   return SF;

 }


 namespace atomic_fs
 { // operators related to fine structure
  
   Operator Darwin(ModelSpace& modelspace, int Z )
   {
     double constants = PI * Z * ALPHA_FS * HBARC*HBARC*HBARC / (2*M_ELECTRON*M_ELECTRON*1e6*1e6) ; // convert to eV and nanometers.  
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
     double oscillator_b = sqrt(HBARC*HBARC/(1e6*M_ELECTRON)/modelspace.GetHbarOmega()); // convert energies to eV, and lengths to nanometers
     double oscillator_b3 = pow(oscillator_b,3);
//     double alpha_FS = 1.0 / 137.035999;
//     double gspin = 2.002319; // electron spin g factor
//     double gspin = ELECTRON_SPIN_G; // electron spin g factor
//     double constants = Z*ALPHA_FS * HBARC*HBARC * gspin / (M_ELECTRON*M_ELECTRON*1e6*1e6) / 32;  // it's 1/8, but we use 4 * LdotS, so 1/32.
     double constants = Z*ALPHA_FS * HBARC*HBARC * ELECTRON_SPIN_G / (M_ELECTRON*M_ELECTRON*1e6*1e6) / 32;  // it's 1/8, but we use 4 * LdotS, so 1/32.
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

     double oscillator_b = sqrt(HBARC*HBARC/(1e6*M_ELECTRON)/modelspace.GetHbarOmega()); // convert energies to eV, lengths to nanometers
     double oscillator_b3 = pow(oscillator_b,3);
//     double alpha_FS = 1.0 / 137.035999;
//     double gspin = 2.002319; // electron spin g factor
//     double constants = - 0.5*alpha_FS *HBARC*HBARC*HBARC/(M_ELECTRON*M_NUCLEON*1e12);  // convert both masses to eV
     double constants = - 0.5*ALPHA_FS *HBARC*HBARC*HBARC/(M_ELECTRON*M_NUCLEON*1e6*1e6);  // convert masses to eV, lengths to nanometers
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
//           Hd.OneBody(a,b) = 4*PI/3 * gspin * sqrt(3./2) * wf0_a * wf0_b ;
           Hd.OneBody(a,b) = 4*PI/3 * ELECTRON_SPIN_G * sqrt(3./2) * wf0_a * wf0_b ;
           Hd.OneBody(b,a) = Hd.OneBody(a,b);
         }
         else
         {
           double r3inv = imsrg_util::RadialIntegral_RpowK(oa.n, oa.l, ob.n, ob.l, -3) / oscillator_b3;
           double L = oa.l!=ob.l ? 0 : sqrt((oa.j2+1.0)/(oa.j2*(oa.j2+2))) * (oa.j2*(oa.j2+2.)/4 +oa.l*(oa.l+1) -3./4);
           double T = modelspace.phase(oa.l) * 3*sqrt(5)*sqrt((oa.j2+1)*(ob.j2+1)*(2*oa.l+1)*(2*ob.l+1)) * AngMom::ThreeJ(oa.l,2,ob.l,0,0,0) * AngMom::NineJ(oa.l,0.5,0.5*oa.j2, ob.l,0.5,0.5*ob.j2, 2,1,1);
//           Hd.OneBody(a,b) = constants * r3inv *( L - gspin/2 * T );
           Hd.OneBody(a,b) = constants * r3inv *( L - ELECTRON_SPIN_G/2 * T );
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
/// \f{equation}{
///  \langle p \| O^{\lambda} \| q \rangle = -\sum_{ia} (n_i -na) \sum_{J}(2J+1)
///  \begin{Bmatrix}
///    j_i & j_a & \lambda \\
///    j_p & j_q & J
///  \end{Bmatrix}
///   \langle i \| O^{\lambda} \| a \rangle \frac{ \tilde{\Gamma}^{J}_{paiq}}{\Delta_{paiq}}
/// \f}
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

//   std::cout << "BEGIN " << __func__ << std::endl;
//   Operator OpIn = OpInx;
   Operator OpOut = 0. * OpIn;
   int Lambda = OpOut.GetJRank();
   size_t norb = OpIn.modelspace->GetNumberOrbits();

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
//       std::cout << "p,q = " << p << " " << q << std::endl;
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
           int Jmin = std::max(  std::max( std::abs(jp-ja), std::abs(ji-jq) ) ,   std::max( std::abs(jp-ji), std::abs(ja-jq) )  );
           int Jmax = std::max(  std::min( jp+ja , ji+jq ) ,   std::min( jp+ji , ja+jq )  );
//           double Delta_paiq = H.OneBody(p,p) + H.OneBody(a,a) - H.OneBody(i,i) - H.OneBody(q,q);
           double Delta_qipa = H.OneBody(q,q) + H.OneBody(i,i) - H.OneBody(p,p) - H.OneBody(a,a);
           double gbar = 0; // gbar_pq`ia`
           for (int J=Jmin; J<=Jmax; J++)
           {
             double sixj = OpIn.modelspace->GetSixJ( ji, ja, Lambda, jp, jq, J );
             double Gamma_paiq = H.TwoBody.GetTBME_J(J,p,a,i,q);
//             double Gamma_paiq = H.TwoBody.GetTBME_J_norm(J,p,a,i,q);
             gbar -= (2*J+1) * sixj * Gamma_paiq;
//             if ( p==2 and q==2 and a==6 and i==0)
//             {
//                std::cout << " line " <<__LINE__ << "  J=" << J << "  :  " << (2*J+1) << " * " << sixj << " * " << Gamma_paiq << "  =>  " << gbar << std::endl;
//             }
         }           
//           Opq += nanifactor * (-gbar) / (Delta_paiq) * Oia;
           Opq += nanifactor * gbar / (Delta_qipa) * Oia;
//           if (p==2 and p==2)
//           {
//           std::cout << "  ai = " << a << " " << i << "  :  " << nanifactor << " * " <<gbar << " * " << Oia
//                     << " / ( " << H.OneBody(q,q) << " + " << H.OneBody(i,i) << " - " << H.OneBody(p,p) << " - " << H.OneBody(a,a) << " )   => Opq = " << Opq << std::endl;
//           }
         }
       }
       OpOut.OneBody(p,q) = Opq;
     }
   }
   return OpOut;
 }




 // Same idea as FirstOrderCorr_1b, except we do a TDA or RPA style resummation.
 // Diagrams like this:
 //                                                 
 //  |    _______                 |                 ^~~~~~X                               |
 //  |   ^      ^                 |   ____         / \                                    |
 //  |  / \    ( )                |  ^    ^       /   \                                   |         ~~~~X 
 //  | (   )    v~~~~X   + ... +  | / \  ( )     (     )     + ... etc.    TDA is just    |      __()
 //  |  \ /                       | \ /   V___    \   /                       this ->     |   __()
 //  |___v                        |__V       ( )   \ /                    type of diagram |__()
 //  |                            |           V_____v                                     |
 //
 //
 //   The n-th order TDA contribution is  TODO this is still not right in the comment. Work this out cleanly.
 //
 //   <p|| O(TDA,n) ||q> = sum_{a,b,...c,d}  sum_{i,j,...k,l}  [ G_pq`ia` / Del_qipa ] 
 //                                                         * [G_ia`jb` / Del_qjpb ] * ... *[G_kc`ld` / Del_qlpd] * <l||O||d>
 //                                                          + [ G_pq`ai` / Del_piqa ]
 //                                                         * [G_ia`jb` / Del_pjqb ] * ... *[G_kc`ld` / Del_plqd] * <d||O||l>
 //   G_ia`jb` is the Pandya-Transformed interaction,and Del is the energy denomonator Del_ijkl = ei + ej - ek -el
 //  
////// I Think this way is wrong. It assumes some symmetries that don't seem to hold.
/*
 Operator RPA_resummed_1b( const Operator& OpIn, const Operator& H, std::string mode )
 {

   // construct hp and ph kets,  as well as Oph and Ohp
   int Lambda = OpIn.GetJRank();
   std::vector<std::pair<size_t,size_t>> ph_kets;
   std::vector<std::pair<size_t,size_t>> hp_kets;
   std::vector<int> ph_phase;

   // maybe we try the dumb way
   for ( auto h : OpIn.modelspace->holes )
   {
     Orbit& oh = OpIn.modelspace->GetOrbit(h);
     for (auto p : OpIn.modelspace->particles )
     {
       Orbit& op = OpIn.modelspace->GetOrbit(p);
       if (op.occ>0.01) continue;
       if ( (oh.j2 + op.j2)<2*Lambda or std::abs(oh.j2-op.j2)>2*Lambda) continue;
       if ( (oh.l + op.l + OpIn.GetParity())%2 > 0 ) continue;
       if ( std::abs( oh.tz2 - op.tz2) != 2*std::abs(OpIn.GetTRank())) continue;
       ph_kets.push_back( std::make_pair(p,h) );
       hp_kets.push_back( std::make_pair(h,p) );
       ph_phase.push_back( AngMom::phase( (oh.j2-op.j2)/2) );
     }
   }
   size_t nkets = ph_kets.size();
   arma::vec Oph( nkets, arma::fill::zeros );
   for (size_t i=0; i<nkets; i++)
   {
     auto p = ph_kets[i].first;
     auto h = ph_kets[i].second;
     Oph( i ) = OpIn.OneBody(p,h);
   }

   
   // Next, construct the Mphph etc matrices
   arma::mat Mphph = GetPH_transformed_Gamma( ph_kets, ph_kets, H, Lambda );

   // The off-diagonal blocks are zero for TDA, but non-zero for RPA
   arma::mat Mphhp( arma::size(Mphph), arma::fill::zeros);

   if (mode=="RPA")
   {
     Mphhp = GetPH_transformed_Gamma( ph_kets, hp_kets, H, Lambda );

      // we need a phase because we multiply this by Oph, which differs from Ohp by (-1)^{jp-jh}.
      arma::Row<double> phases(nkets, arma::fill::ones);
      for ( size_t i=0; i<nkets; i++) phases(i) = ph_phase[i];
      Mphhp.each_row() %= phases;

   }
   if ( mode=="CP" ) // CP means first order core-polarization correction. M is zero because we don't want the iterations.
   {
     Mphph.zeros();
   }


   arma::mat Delta(arma::size(Mphph), arma::fill::ones );

   for ( size_t i=0; i<nkets; i++)
   {
     size_t p = ph_kets[i].first;
     size_t h = ph_kets[i].second;
     double del_ph = H.OneBody(h,h) - H.OneBody(p,p);
     Delta.col(i).fill(del_ph);  // denominator for the Oph bits
   }




  // At this point, we can't get around the fact that the denominators depend on the initial and final state of the
  // entire operator, as in < p | Oeff | q >, and we need to add ep-eq to the denominator. So for every p,q pair
  // we need to do the matrix inversion and evaluation in 1st order PT.

   Operator OpOut = OpIn;
   OpOut.OneBody.zeros();

   for ( auto v1 : OpIn.modelspace->valence )
   {
     for ( auto v2 : OpIn.modelspace->valence )
     {

      arma::mat del12(  arma::size(Delta) );
      del12.fill(  H.OneBody(v1,v1) - H.OneBody(v2,v2) );

      arma::mat M = Mphph / (Delta + del12 ) + Mphhp/(Delta - del12);


      // do the fancy resummation of the series by matrix inversion
      // Minv  =  (I-M)^-1  where I is the identity matrix.
       arma::mat ONE = arma::eye(arma::size(M));
       arma::mat Minv = arma::inv(  ONE  - M );

       arma::vec ORPA = Minv * Oph;

       // now we need to unpack ORPA into a matrix again.
       Operator OpRPA = OpIn;
       OpRPA.OneBody.zeros();
       for ( size_t i=0; i<nkets; i++ )
       {
         size_t p = ph_kets[i].first;
         size_t h = ph_kets[i].second;
         
         OpRPA.SetOneBody(p,h,  ORPA(i) ) ;
       }

      // evaluate this in first order perturbation theory, and take the valence part that we're interested in
       double op12 = FirstOrderCorr_1b( OpRPA, H).OneBody(v1,v2);
       OpOut.OneBody(v1,v2) += op12;

      }
    }
   return OpOut;

 }
*/

 Operator RPA_resummed_1b( const Operator& OpIn, const Operator& H, std::string mode )
 {

   // construct hp and ph kets,  as well as Oph and Ohp
   int Lambda = OpIn.GetJRank();
   std::vector<std::pair<size_t,size_t>> ph_kets;
   std::vector<std::pair<size_t,size_t>> hp_kets;
   std::vector<int> ph_phase;

   // maybe we try the dumb way
   for ( auto h : OpIn.modelspace->holes )
   {
     Orbit& oh = OpIn.modelspace->GetOrbit(h);
//     for (auto p : OpIn.OneBodyChannels.at({oh.l,oh.j2,oh.tz2}) )
     for (auto p : OpIn.modelspace->particles )
     {
       Orbit& op = OpIn.modelspace->GetOrbit(p);
       if (op.occ>0.01) continue;
//       Oph( ph_kets.size() ) = OpIn.OneBody(p,h);
//       Ohp( hp_kets.size() ) = OpIn.OneBody(h,p);
       if ( (oh.j2 + op.j2)<2*Lambda or std::abs(oh.j2-op.j2)>2*Lambda) continue;
       if ( (oh.l + op.l + OpIn.GetParity())%2 > 0 ) continue;
       if ( std::abs( oh.tz2 - op.tz2) != 2*std::abs(OpIn.GetTRank())) continue;
       ph_kets.push_back( std::make_pair(p,h) );
       hp_kets.push_back( std::make_pair(h,p) );
       ph_phase.push_back( AngMom::phase( (oh.j2-op.j2)/2) );
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
//     std::cout << " setting  i " << i << "  ph " << p << " " << h << "  Oph " << Oph(i) << "  Ohp " << Ohp(i) << std::endl;
   }

   
   // Next, construct the Mphph etc matrices
   arma::mat Mphph = GetPH_transformed_Gamma( ph_kets, ph_kets, H, Lambda );
   arma::mat Mhphp = GetPH_transformed_Gamma( hp_kets, hp_kets, H, Lambda );

   // The off-diagonal blocks are zero for TDA, but non-zero for RPA
   arma::mat Mphhp( arma::size(Mphph), arma::fill::zeros);
   arma::mat Mhpph( arma::size(Mphph), arma::fill::zeros);

//   arma::mat Mphhp = GetPH_transformed_Gamma( ph_kets, hp_kets, H, Lambda );
//   arma::mat Mhpph = GetPH_transformed_Gamma( hp_kets, ph_kets, H, Lambda );
   // to get TDA, we just get rid of the off-diagonal blocks
   if (mode=="RPA")
   {
     Mphhp = GetPH_transformed_Gamma( ph_kets, hp_kets, H, Lambda );
     Mhpph = GetPH_transformed_Gamma( hp_kets, ph_kets, H, Lambda );
   }
   if ( mode=="CP" ) // CP means first order core-polarization correction. M is zero because we don't want the iterations.
   {
     Mphph.zeros();
     Mhphp.zeros();
   }
//   if (mode=="TDA")
//   {
//     Mphhp *=0;
//     Mhpph *=0;
//   }

   // for the full M matrix  M = [ Mphph   Mphhp ]
   //                            [ Mhpph   Mhphp ]
   arma::mat M = arma::join_vert( arma::join_horiz( Mphph, Mphhp) ,
                                  arma::join_horiz( Mhpph, Mhphp) );

//   std::cout << "Size of Mphph is " << Mphph.n_rows << " x " << Mphph.n_cols << std::endl;
   // make the base case denominator eq-ep for ket |pq>
   arma::mat Delta(arma::size(M), arma::fill::ones );
   for ( size_t i=0; i<nkets; i++)
   {
     size_t p = ph_kets[i].first;
     size_t h = ph_kets[i].second;
//     double del_ph = H.OneBody(p,p) - H.OneBody(h,h);
     double del_ph = H.OneBody(h,h) - H.OneBody(p,p);
//     Delta.col(i) *= -del_ph;
//     Delta.col(i+nkets) *= del_ph;
//     Delta.col(i).fill(-del_ph);  // denominator for the Oph bits
//     Delta.col(i+nkets).fill(-del_ph);  // denominator for the Ohp bits
     Delta.col(i).fill(del_ph);  // denominator for the Oph bits
     Delta.col(i+nkets).fill(del_ph);  // denominator for the Ohp bits
   }

//   arma::mat tmp = Mphph/Delta.submat(0,0,nkets-1,nkets-1)*Oph;
//   std::cout << "Mphph is " << std::endl << Mphph << std::endl << std::endl << "Delta is " << std::endl << Delta << std::endl << "Oph is " << std::endl << Oph << std::endl << "Mphph*Oph is " << std::endl << tmp << std::endl;

   // combine Oph and Ohp into a single column vector
   arma::vec Ovec = arma::join_vert(  Oph, Ohp );



  // At this point, we can't get around the fact that the denominators depend on the initial and final state of the
  // entire operator, as in < p | Oeff | q >, and we need to add ep-eq to the denominator. So for every p,q pair
  // we need to do the matrix inversion and evaluation in 1st order PT.

   Operator OpOut = OpIn;
   OpOut.OneBody.zeros();
//   std::cout << "Mphph = " << std::endl << Mphph << std::endl;
//   std::cout << "Mphhp = " << std::endl << Mphhp << std::endl;
//   std::cout << "Mhpph = " << std::endl << Mhpph << std::endl;
//   std::cout << "Mhphp = " << std::endl << Mhphp << std::endl;

   for ( auto v1 : OpIn.modelspace->valence )
   {
     for ( auto v2 : OpIn.modelspace->valence )
     {

//      if (v1 != v2) continue; // GET RID OF THIS!!

//      arma::mat del12 = arma::ones(arma::size(Delta)) * (H.OneBody(v1,v1) - H.OneBody(v2,v2)) ;
      arma::mat del12(  arma::size(Delta) );
      del12.fill(  H.OneBody(v1,v1) - H.OneBody(v2,v2) );
      del12.cols(nkets,2*nkets-1) *=-1;
      // do the fancy resummation of the series by matrix inversion
      // Minv  =  (I-M)^-1  where I is the identity matrix. The slash here means element-wise division.
       arma::mat Minv = arma::inv(  arma::eye(arma::size(M))  - M/( Delta+del12) );

//       arma::mat Mdel = M/(Delta+del12); // FOR DEBUGGING. REMOVE THIS
//       Minv = arma::eye(arma::size(M)) + Mdel;// + Mdel*Mdel;   // FOR DEBUGGING. REMOVE THIS
       arma::vec ORPA = Minv * Ovec;

//       arma::mat mphph_tmp = Mphph/(Delta.submat(0,0,nkets-1,nkets-1) + del12.submat(0,0,nkets-1,nkets-1) );
//       arma::mat mhphp_tmp = Mhphp/(Delta.submat(nkets,nkets,2*nkets-1,2*nkets-1) + del12.submat(nkets,nkets,2*nkets-1,2*nkets-1) );
//       arma::mat mphhp_tmp = Mphhp/(Delta.submat(0,nkets,nkets-1,2*nkets-1) + del12.submat(0,nkets,nkets-1,2*nkets-1) );
//       arma::mat mphhp_alt = mphhp_tmp;
//       arma::Row<double> phases(nkets, arma::fill::ones);
//       for ( size_t i=0; i<nkets; i++) phases(i) = ph_phase[i];
//       mphhp_alt.each_row() %= phases;
//
////       std::cout << "mphhp_tmp = " << std::endl << mphhp_tmp << std::endl << std::endl << " mphhp_alt = " << std::endl << mphhp_alt << std::endl;
//       arma::mat mhpph_tmp = Mhpph/(Delta.submat(nkets,0,2*nkets-1,nkets-1) + del12.submat(nkets,0,2*nkets-1,nkets-1) );
//       arma::vec Orpa_ph = mphph_tmp * Oph + mphhp_tmp * Ohp;
//       arma::vec Orpa_ph_alt = (mphph_tmp + mphhp_alt) * Oph;
//       arma::vec Orpa_hp = mhphp_tmp * Ohp + mhpph_tmp * Oph;
    
//       std::cout << "====================================  " << v1 << " " << v2 << " ===============================" << std::endl;
//       for (size_t i=0; i<Oph.n_rows;i++)
//       {
//          std::cout << "i " << i << "    Oph, Ohp  " << Oph(i) << "   " << Ohp(i) << "    Orpa_ph,hp " << Orpa_ph(i) << " , " <<  Orpa_ph_alt(i) << "   " << Orpa_hp(i)
//                      << "  phase " << ph_phase[i] << "       " << Orpa_ph(i) - ph_phase[i]*Orpa_hp(i) << std::endl;
//       }
       // now we need to unpack ORPA into a matrix again.
       Operator OpRPA = OpIn;
       OpRPA.OneBody.zeros();
       for ( size_t i=0; i<nkets; i++ )
       {
         size_t p = ph_kets[i].first;
         size_t h = ph_kets[i].second;
        // This is risky business, if in the future we assume some symmetries about these matrix elements.
        // We have symmetry for  v1 <-> v2, but we don't have a symmetry for p <-> h for a fixed v1,v2.
         OpRPA.OneBody(p,h) = ORPA(i);
         OpRPA.OneBody(h,p) = ORPA(i+nkets);
//         OpRPA.SetOneBody(p,h,  ORPA(i) ) ;
//         OpRPA.OneBody(p,h) = ORPA(i) + ORPA(i+nkets);
//         OpRPA.OneBody(h,p) = OpRPA.OneBody(p,h);
//          std::cout << "i " << i << "   ph " << p << "  " << h << "   Oph " << OpRPA.OneBody(p,h) << "   Ohp " << OpRPA.OneBody(h,p) << std::endl;
//         OpRPA.SetOneBody(p,h,  ORPA(i) ) ;
//          std::cout << "  " << i << "   ph " << p << "  " << h << "   Oph " << OpRPA.OneBody(p,h) << "   Ohp " << OpRPA.OneBody(h,p) << std::endl;
       }

//       // BRUTE FORCE TDA2
//       Orbit ov1 = OpIn.modelspace->GetOrbit(v1);
//       Orbit ov2 = OpIn.modelspace->GetOrbit(v2);
//       double jv1 = ov1.j2*0.5;
//       double jv2 = ov2.j2*0.5;
//       double tda1 = 0;
//       double tda2 = 0;
//       for ( auto a : OpIn.modelspace->particles)
//       {
//         for ( auto i : OpIn.modelspace->holes)
//         {
//           Orbit& oa = OpIn.modelspace->GetOrbit(a);
//           Orbit& oi = OpIn.modelspace->GetOrbit(i);
//           double ja = oa.j2*0.5;
//           double ji = oi.j2*0.5;
//           int Jmin = std::max( std::abs(ov1.j2-oa.j2), std::abs(ov2.j2-oi.j2))/2;
//           int Jmax = std::min( (ov1.j2+oa.j2), (ov2.j2+oi.j2))/2;
//           double gbar1up = 0;
//           for (int J=Jmin; J<=Jmax; J++)
//           {
//             gbar1up -= (2*J+1) * OpIn.modelspace->GetSixJ( jv1, ja, J, ji, jv2, Lambda) * H.TwoBody.GetTBME_J(J,v1,a,i,v2);
//           }
//           double denom1up = H.OneBody(v2,v2) + H.OneBody(i,i) - H.OneBody(v1,v1) - H.OneBody(a,a);
//
//           Jmin = std::max( std::abs(ov1.j2-oi.j2), std::abs(ov2.j2-oa.j2))/2;
//           Jmax = std::min( (ov1.j2+oi.j2), (ov2.j2+oa.j2))/2;
//           double gbar1down = 0;
//           for (int J=Jmin; J<=Jmax; J++)
//           {
//             gbar1down -= (2*J+1) * OpIn.modelspace->GetSixJ( jv1, ji, J, ja, jv2, Lambda) * H.TwoBody.GetTBME_J(J,v1,i,a,v2);
//           }
//           double denom1down = H.OneBody(v1,v1) + H.OneBody(i,i) - H.OneBody(v2,v2) - H.OneBody(a,a);
//
//           tda1 += gbar1up/denom1up * OpIn.OneBody(i,a);
//           tda1 += gbar1down/denom1down * OpIn.OneBody(a,i);
//           if (v1==7 and std::abs(OpIn.OneBody(i,a))>1e-6)
//           {
//           std::cout << " i,a = " << i << " " << a << "   tda1up " << gbar1up << " / " << denom1up << " * " << OpIn.OneBody(i,a) << " = " << gbar1up/denom1up * OpIn.OneBody(i,a) << std::endl;
//           std::cout << "       " << i << " " << a << "   tda1down " << gbar1down << " / " << denom1down << " * " << OpIn.OneBody(a,i) << " = " << gbar1down/denom1down * OpIn.OneBody(a,i) << std::endl;
//            }
//
//           for ( auto b : OpIn.modelspace->particles)
//           {
//              for ( auto j : OpIn.modelspace->holes)
//              {
//                Orbit& ob = OpIn.modelspace->GetOrbit(b);
//                Orbit& oj = OpIn.modelspace->GetOrbit(j);
//                double jb = ob.j2*0.5;
//                double jj = oj.j2*0.5;
//                int Jpmin = std::max( std::abs(oi.j2-ob.j2), std::abs(oj.j2-oa.j2))/2;
//                int Jpmax = std::min( (oi.j2+ob.j2), (oj.j2+oa.j2))/2;
//                double gbar2up = 0;
//                for (int Jp=Jpmin; Jp<=Jpmax; Jp++)
//                {
//                  gbar2up -= (2*Jp+1) * OpIn.modelspace->GetSixJ( ji, jb, Jp, jj, ja, Lambda) * H.TwoBody.GetTBME_J(Jp,i,b,j,a);
//                }
//                double denom2up = H.OneBody(v2,v2) + H.OneBody(j,j) - H.OneBody(v1,v1) - H.OneBody(b,b);
//
//                Jpmin = std::max( std::abs(oi.j2-ob.j2), std::abs(oj.j2-oa.j2))/2;
//                Jpmax = std::min( (oi.j2+ob.j2), (oj.j2+oa.j2))/2;
//                double gbar2down = 0;
//                for (int Jp=Jpmin; Jp<=Jpmax; Jp++)
//                {
//                  gbar2down -= (2*Jp+1) * OpIn.modelspace->GetSixJ( ji, jb, Jp, jj, ja, Lambda) * H.TwoBody.GetTBME_J(Jp,i,b,j,a);
//                }
//                double denom2down = H.OneBody(v1,v1) + H.OneBody(j,j) - H.OneBody(v2,v2) - H.OneBody(b,b);
//
//                tda2 += (gbar1up)/denom1up * (gbar2up)/denom2up * OpIn.OneBody(j,b);
//                tda2 += (gbar1down)/denom1down * (gbar2down)/denom2down * OpIn.OneBody(b,j);
//              }
//            }
//         }
//       }
//       std::cout << " v1 v2 =  "<< v1 << " " << v2 << "   tda1, tda2 = " << tda1 << " " << tda2 << std::endl;

//       if (v1==2 and v2==2)
//       {
//          std::cout << "========== THIS ONE: ===========" << std::endl;
//       }
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





 // Evaluate <bra | r1xp2 | ket>, omitting the factor (hbar * omega) /(m * omega^2)
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



 // This is an attempt to implement the contact term in the neutrinoless double beta decay operator
 // presented in Cirigliano et al PRL 120 202001 (2018)
 // I implemented the cutoff-regulated version given in r-space by
 // V(r) = -2 g_nuNN(R) delta_R(r)
 // where delta_R is a Gaussian regulated delta function
 // delta_R(r) = 1/ (sqrt(pi)R)^3 * exp(-r^2/R^2)
 // for cutoff R.  This delta is normalized to 1 when integrated over x,y,x.
 // So in spherical coordinates, the angular integration picks up a factor 4pi and the radial integration comes with r^2dr.
 //  The low energy constant g_nuNN is estimated by Cirigliano et al to be order fpi^-2
 //  where fpi is the pion decay constant. To make things dimensionless, we should 
 // After discussing with Javier Menendez, the correct thing to do is to compute the dimensionless quantity
 // (r0 A^1/3) / mpi^2  * delta_R(r)
 // Since g_nuNN is order fpi^-2, the combination fpi^2 * g_nuNN is a dimensionless number of order 1
 // The usual neutrino exchange matrix element has units of MeV (or fm^-1), but is multiplied by the nuclear radius R=1.2 A^1/3
 // to make it dimensionless (and typically of order 1). This nuclear radius factor is compensated in the phase space factor.
 // So in the end, if the expectation value of R/fpi^2 * delta_R(r) is of order 1, it is non-negligible compared to the neutrino exchange operator.
 Operator M0nu_contact_Op(ModelSpace& modelspace, double R0 )
 {
   double t_start = omp_get_wtime();
   Operator M0nuCT(modelspace, 0,2,0,2);
   double oscillator_b = sqrt(HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
//   std::cout << "oscillator b = " << oscillator_b << std::endl;

   // In Moshinsky's convention, r_rel = (r1-r2)/sqrt(2).  We want ( |r1-r2|^2 / R0^2 )  =  ( 2 r_rel^2 / R0^2 ) => 2 sigma^2 = R0^2/2 => sigma = R0/2
   double sigma = R0/2.0 / oscillator_b ;  // we work in units of the oscillator length b
//   double normalization = 1.0 / pow( sqrt(M_PI * R0/HBARC), 3 );
   double normalization =4*PI *  pow( SQRTPI * R0, -3 ); // this has units fm^-3. The 4pi comes from integration over angles
   int A = modelspace.GetTargetMass(); // The mass of the thing we will be calculating
   double nuclear_radius = 1.2 * pow( A, 1./3);  // empirical estimate for nuclear radius, has units of fm
//   double mpi = 134.0; // pion mass, or whatever..
//   double mpi2 = mpi*mpi / (HBARC*HBARC); // pion mass squared, in units of fm^-2, (a.k.a inverse compton wavelength squared)
//   normalization *= nuclear_radius / mpi2; // now the normalization is dimensionless
   //
//   double fpi = 92.2; // pion decay constant in MeV
//   normalization *= nuclear_radius * HBARC*HBARC / (fpi*fpi); // now the normalization is dimensionless
   normalization *= nuclear_radius * HBARC*HBARC / (F_PI*F_PI); // now the normalization is dimensionless

   // Making the units work out. In the future, this all should be done consistently with the interaction
//   double Ctilde = -0.4 / (fpi*fpi);
//   double C1tilde = 2.0;  // from fitting LO scattering at R0=0.5  (see Cirigliano et al. PRL 120 202001 )
//   double C1 = pow(M_NUCLEON*Ctilde/4/M_PI, 2) * C1tilde;
////   double g_nuNN = C1;
   double g_nuNN = 1;
   normalization *= -2*g_nuNN;

   //std::cout << "Constructing normalization from " << pow( sqrt(3.14159) * R0, -3 ) << " * " << nuclear_radius << " * " << HBARC*HBARC/(fpi*fpi) << " * -2  = " << normalization << std::endl;

//   std::vector<double> rgrid;
//   for (double r=0; r<10.1; r+=0.1) rgrid.push_back(r);

//   std::vector<double> rgrid = { 0.10000 ,0.20000 ,0.30000 ,0.40000 ,0.50000 ,0.60000 ,
//                                 0.70000 ,0.80000 ,0.90000 ,1.00000 ,1.50000 ,2.00000 ,
//                                 2.50000 ,3.00000 ,3.50000 ,4.00000 ,4.50000 ,5.00000 ,
//                                 5.50000 ,6.00000 ,6.50000 , 7.00000 ,7.50000 ,8.00000 };


   modelspace.PreCalculateMoshinsky();
//   std::cout << "Done Precalculating Moshinsky." << std::endl;
   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   std::vector<std::array<size_t,2>> braketchannels;
   for ( auto channel : M0nuCT.TwoBody.MatEl ) braketchannels.push_back(channel.first);
   int nchan = braketchannels.size();
//   #pragma omp parallel for schedule(dynamic,1)  
   for (int ch=0; ch<nchan; ++ch)
   {
      int chbra = braketchannels[ch][0];
      int chket = braketchannels[ch][1];
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra);
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket);
//      if (tbc.Tz >= 0) continue; // 2-body coulomb only acts in pp channel
      int lmax = modelspace.GetLmax();
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra=0;ibra<nbras;++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         Orbit & oa = modelspace.GetOrbit(bra.p);
         Orbit & ob = modelspace.GetOrbit(bra.q);
         int na = oa.n;
         int nb = ob.n;
         int la = oa.l;
         int lb = ob.l;
         double ja = oa.j2*0.5;
         double jb = ob.j2*0.5;
         int fab = 2*na + 2*nb + la + lb;

         for (int iket=0;iket<nkets;++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);

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

            double mcont=0;
//            bool verbose = false;
//            if (J==0 and la==3 and lb==3 and lc==3 and ld==3 and na+nb+nc+nd==0 and oa.j2==7 and ob.j2==7 and oc.j2==7 and od.j2==7) verbose = true;
//            if (verbose) std::cout << "osc_b, R0, sigma = " << oscillator_b << " " << R0 << " " << sigma << "  nuc Radius, norm = " << nuclear_radius << " " << normalization  << std::endl;
//            std::vector<double> psirel( rgrid.size(), 0);
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
//                if (verbose) std::cout << "L S = " << Lab << " " << Sab << "   norm nine-j = " << njab << " , " << njcd << std::endl;
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
//                       int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
                       // matrix elements are < pp | M | nn >, so we need (-1)(lam + S) to be positive
                       int asymm_factor = (1 + 1*modelspace.phase( lam_ab + Sab ))/ 2;
                       if ( asymm_factor ==0 ) continue;
         
                       int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
                       int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation
         
                       double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                       if (std::abs(mosh_ab)<1e-8) continue;
         
                       int N_cd = N_ab;
                       int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                       if (n_cd < 0) continue;
//                       if  (n_ab != n_cd and N_ab != N_cd) continue;
         
                       double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                       if (std::abs(mosh_cd)<1e-8) continue;

                       double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;

                       double rad_int =  RadialIntegral_Gauss( n_ab,lam_ab, n_cd,lam_cd, sigma) ;  
 
                       mcont += prefactor * rad_int; 
//                       if (verbose) std::cout << "   Nab Lamab nab lamab = " << N_ab << " " << Lam_ab << " " << n_ab << " " << lam_ab
//                                              << "  Ncd Lamcd ncd lamcd = " << N_cd << " " << Lam_cd << " " << n_cd << " " << lam_cd
//                                              << " mosh_ab , mosh_cd = " << mosh_ab << " " << mosh_cd << "  asymm_factor = " << asymm_factor
//                                              << " rad_int = " << rad_int << "   mcont = " << mcont << std::endl;
//
//                      if (verbose)
//                      {
//                        for ( size_t i=0; i<rgrid.size(); i++ )
//                        {
//                          psirel[i] += prefactor * pow(2,-1.5)* HO_density(n_ab, lam_ab, modelspace.GetHbarOmega(), rgrid[i] /sqrt(2) );
//                          
//                        }
//                      }

    
                    } // lam_ab
                  } // Lam_ab
                } // N_ab
         
              } // Sab
            } // Lab
//            if (verbose)
//            {
//              for (size_t i=0; i<rgrid.size(); i++)
//              {
//                std::cout << rgrid[i] << "   " << psirel[i] << std::endl;
//              }
//            }

            mcont *=  normalization / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq())); // normalize 
//            if (verbose) std::cout << " with normalization " << normalization << " , setting matrix element to " << mcont << std::endl;
            M0nuCT.TwoBody.SetTBME(chbra,chket,ibra,iket,mcont);
                         
         }
      }
   }

   M0nuCT.profiler.timer["M0nu_contact_Op"] += omp_get_wtime() - t_start;
   return M0nuCT ;
 }


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

/// Minnesota potential of Thompson, Lemere and Tang,  Nuc. Phys.A 286 53 (1977)
/// I have modified the triplet channel strength by multiplying by 0.2.
/// Without this modification, it seems that the triplet channel does not have enough
/// repulsion to support against collapse, which makes the potential not terribly
/// useful for finite nuclei.
/// This is mostly useful for benchmarking without needing large input files.
///
/// On March 12 2019, I calculated O16 in a basis of hw=16 and emax=4, 6, and 8.
/// I used IMSRG(2) with method=magnus  omega_norm_max=0.25, dsmax=0.5
/// Results:
///  emax = 4                         emax = 6                  emax = 8
///  e1hf = 177.4709518              e1hf = 132.7170792          e1hf = 111.4959623
///  e2hf = -102.1762390             e2hf = -138.2451400         e2hf = -210.3298042
///  e3hf = 0.0000000                e3hf = 0.0000000            e3hf = 0.0000000    
///  EHF = 75.2947127                EHF  = -5.5280608           EHF = -98.8338419   
///  EIMSRG = -10.047908             EIMSRG = -69.170313         EIMSRG =  -144.9899081
///  Rp2HF = 9.5591609               Rp2HF = 15.1991899          Rp2HF =  20.8748671
///  Rp2IMSRG = 9.362856             Rp2IMSRG = 14.399832        Rp2IMSRG = 21.353270
// Note that the triplet channel strength VT is scaled by 0.2 to
// give non-crazy results for finite nuclei.
 Operator MinnesotaPotential( ModelSpace& modelspace )
 {
   double VR = 200 ; // MeV
////   double VT = -178. ; // MeV
//   double VT = -178. * 0.2; // scaled to make things not crazy for finite nuclei
   double VT = -178. * 0.4; // scaled to make things not crazy for finite nuclei
   double VS = -91.85 ; // MeV
   double kR = 1.487 ; // in fm^-2
   double kT = 0.639 ; // in fm^-2
   double kS = 0.465 ; // in fm^-2
   std::array<double,6> params = {VR,VT,VS, kR,kT,kS};
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
//       double vminn = MinnesotaMatEl( modelspace, bra, ket, J);
       double vminn = MinnesotaMatEl( modelspace, bra, ket, J, params);
       Vminnesota.TwoBody.SetTBME(ch,ibra,iket,vminn);
       Vminnesota.TwoBody.SetTBME(ch,iket,ibra,vminn);
     }
    }
   }
   return Vminnesota;
 }

 // A gaussian potential with strength 1.0
 // We just reuse the Minnesota potential code and set the singlet and triplet strengths to zero
 Operator GaussianPotential( ModelSpace& modelspace, double sigma )
 {
   double VR = 1.0 ;
////   double VT = -178. ;
//   double VT = -178. * 0.2; // scaled to make things not crazy for finite nuclei
//   double VS = -91.85 ;
   double kR = 1.0 / (2*sigma*sigma); // make things dimensionless
//   double kR = 1.487 * oscillator_b2; // make things dimensionless
//   double kT = 0.639 * oscillator_b2;
//   double kS = 0.465 * oscillator_b2;
   std::array<double,6> params = {VR,0.0,0.0, kR,1,1};
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
//       double vminn = MinnesotaMatEl( modelspace, bra, ket, J);
       double vminn = MinnesotaMatEl( modelspace, bra, ket, J, params);
       Vminnesota.TwoBody.SetTBME(ch,ibra,iket,vminn);
       Vminnesota.TwoBody.SetTBME(ch,iket,ibra,vminn);
     }
    }
   }
   return Vminnesota;
 }

 // A gaussian potential with strength 1.0
 // We just reuse the Minnesota potential code and set the singlet and triplet strengths to zero
 Operator GaussianCSBPotential( ModelSpace& modelspace, double sigma )
 {
   double VR = 1.0 ;
////   double VT = -178. ;
//   double VT = -178. * 0.2; // scaled to make things not crazy for finite nuclei
//   double VS = -91.85 ;
   double kR = 1.0 / (2*sigma*sigma); // make things dimensionless
//   double kR = 1.487 * oscillator_b2; // make things dimensionless
//   double kT = 0.639 * oscillator_b2;
//   double kS = 0.465 * oscillator_b2;
   std::array<double,6> params = {VR,0.0,0.0, kR,1,1};
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
//       double vminn = MinnesotaMatEl( modelspace, bra, ket, J);
       double vminn = -tbc.Tz * MinnesotaMatEl( modelspace, bra, ket, J, params);
       Vminnesota.TwoBody.SetTBME(ch,ibra,iket,vminn);
       Vminnesota.TwoBody.SetTBME(ch,iket,ibra,vminn);
     }
    }
   }
   return Vminnesota;
 }


 // A bare delta function of the relative coordinate (effectively regulated by the finite oscillator basis).
 // The delta is a 1D delta of the relative coordinate, as defined by Moshinsky: r = 1/sqrt(2) |r1-r2|.
 // To convert this to the usual definition, we have delta(|r1-r2|) = 1/sqrt(2) delta(r)
 // Evaluating the diagonal matrix element between two 0s states, we get
 // < 0s 0s | delta(r) | 0s 0s> = 4 / [sqrt(pi) b^3], where b is the oscillator length.
 Operator BareDelta( ModelSpace& modelspace )
 {
   Operator Vdel = Operator(modelspace,0,0,0,2);

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();

   double hw = modelspace.GetHbarOmega();
//   double b_osc = HBARC / sqrt( hw * M_NUCLEON);

   // it's a delta, so we only need the wave function evaluated at r=0. pre-store those
//   int nmax = modelspace.GetEmax() /2;
   int nmax = modelspace.GetEmax();
   std::vector<double> psi_at_zero(nmax+1, 0.) ;
   for (int n=0; n<=nmax; n++) psi_at_zero[n] = HO_Radial_psi( n, 0, hw, 0.0);
   
   // the spins are always 1/2
   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
    TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;++ibra)
    {
     Ket& bra = tbc.GetKet(ibra);
     Orbit & oa = modelspace.GetOrbit(bra.p);
     Orbit & ob = modelspace.GetOrbit(bra.q);
     int na = oa.n;
     int nb = ob.n;
     int la = oa.l;
     int lb = ob.l;
     double ja = oa.j2/2.0;
     double jb = ob.j2/2.0;
     int fab = 2*na + 2*nb + la + lb; // oscillator energy in |ab>


     for (int iket=ibra;iket<nkets;iket++)
     {
       Ket& ket = tbc.GetKet(iket);

       double vdelta =0;

       Orbit & oc = modelspace.GetOrbit(ket.p);
       Orbit & od = modelspace.GetOrbit(ket.q);
       int nc = oc.n;
       int nd = od.n;
       int lc = oc.l;
       int ld = od.l;
       double jc = oc.j2/2.0;
       double jd = od.j2/2.0;
       int fcd = 2*nc + 2*nd + lc + ld;
    

       // First, transform to LS coupling using 9j coefficients
       for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
       {
         if ( Lab%2 != fab%2 )  continue; // lam=0, so parity needs to come from COM
         for (int Sab=0; Sab<=1; ++Sab)
         {
           if ( std::abs(Lab-Sab)>J or Lab+Sab<J) continue;
    
           // the delta doesn't act on spin or angular components
           int Scd = Sab;
           int Lcd = Lab;
           double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
           if (std::abs(njab)<1e-8) continue;
           double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
           if (std::abs(njcd)<1e-8) continue;


          int lam_ab = 0; // radial delta is nonzero only for s-waves
          int lam_cd = lam_ab; // radial delta conserves orbital angular momentum
          int Lam_ab = Lab; // so the total Lab needs to be all COM motion.
          int Lam_cd = Lam_ab; // center of mass is untouched

          // factor to account for antisymmetrization   
          int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
          if ( asymm_factor ==0 ) continue;

          // Next, transform to rel / com coordinates with Moshinsky tranformation
          for (int N_ab=0; N_ab<=(fab-Lam_ab)/2; ++N_ab)  // N_ab = CoM n for a,b
          {

            int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

            int N_cd = N_ab; // com is untouched by the delta
            int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
            if (n_cd < 0) continue;

   
            double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
            if (std::abs(mosh_ab)<1e-8) continue;

            double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
            if (std::abs(mosh_cd)<1e-8) continue;

            double prefactor =  njab * njcd * mosh_ab * mosh_cd * asymm_factor;
            vdelta += psi_at_zero[n_ab] * psi_at_zero[n_cd] * prefactor;

          }// for N_ab
         }// for Sab
       }// for Lab

       // This should be normalized.
       if (bra.p==bra.q) vdelta /= SQRT2;
       if (ket.p==ket.q) vdelta /= SQRT2;
       Vdel.TwoBody.SetTBME(ch,ibra,iket,vdelta);
       Vdel.TwoBody.SetTBME(ch,iket,ibra,vdelta);
     }// for iket
    }// for ibra
   }// for ch

   return Vdel;
 }




/// 
/// A single matrix element of the Minnesota potential. See more details above MinnesotaPotential
 double MinnesotaMatEl( ModelSpace& modelspace, Ket& bra, Ket& ket, int J, const std::array<double,6>& params )
 {
   double oscillator_b2 = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
//   double VR = 200 ;
////   double VT = -178. ;
//   double VT = -178. * 0.2; // scaled to make things not crazy for finite nuclei
//   double VS = -91.85 ;
//   double kR = 1.487 * oscillator_b2; // make things dimensionless
//   double kT = 0.639 * oscillator_b2;
//   double kS = 0.465 * oscillator_b2;
   double VR = params[0];
   double VT = params[1];
   double VS = params[2];
   double kR = params[3]* oscillator_b2;
   double kT = params[4]* oscillator_b2;
   double kS = params[5]* oscillator_b2;

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
   if ((la+lb+lc+ld)%2 !=0) return 0; // check parity conservation
//   if (std::abs(fab+fcd)%2 >0) return 0; // check parity conservation
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
//       if ( (oa.tz2 == ob.tz2) and Sab==1 ) continue; // WRONG!!!
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


//              for (int N_cd=std::max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              for (int N_cd=N_ab; N_cd<=N_ab; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
//                if  (n_ab != n_cd and N_ab != N_cd) continue;

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

                // for debugging
//                if ( (n_ab==n_cd)  and (lam_ab==lam_cd)  and (N_ab==N_cd)  and (Lam_ab==Lam_cd) and (Sab==Scd) and (Lab==Lcd) )  vrel =1;
//                if ( (n_ab==n_cd) and  (N_ab==N_cd)  )  vrel =1;
//                if ( true )  vrel =1;
                Vminn += vrel * njab * njcd * mosh_ab * mosh_cd * asymm_factor;

//                if ( (J==0) and ((bra.op->l+bra.oq->l)%2==0)  and ((bra.op->tz2+bra.oq->tz2)==-2) )
//                {
//                   std::cout << "  abcd << " << bra.p << " " << bra.q << " " << ket.p << " " << ket.q 
//                             << "  nab lab  ncd lcd  " << n_ab << " " << lam_ab << " " << n_cd << " " << lam_cd
//                             << "  Nab Lab  Ncd Lcd  " << N_ab << " " << Lam_ab << " " << N_cd << " " << Lam_cd
//                             << " L S " << Lab << " " << Lcd << "  " << Sab << " " << Scd 
//                             << "    adding " << vrel << " * " << njab << " * " << njcd << " * " << mosh_ab << " * " << mosh_cd
//                             << " * " << asymm_factor << "  -> " << Vminn << std::endl;
//                }


              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // Dont forget to normalize!
   if (bra.p==bra.q) Vminn /= SQRT2;
   if (ket.p==ket.q) Vminn /= SQRT2;

   return Vminn ;


 }



 // Serber-type potential V = f(r) [ A + B sig*sig + C tau*tau + D (sig*sig)(tau*tau) ]
 // with f(r) = V_0 exp[-r^2/mu^2]
 // See Suhonen section 8.1.4
 Operator SerberTypePotential( ModelSpace& modelspace, double V0, double mu, double A, double B, double C, double D)
 {
//   double VR = 1.0 ;
   std::array<double,6> params = {V0,mu,A,B,C,D};
   Operator Vserber(modelspace, 0,0,0,2);

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
//       double vminn = MinnesotaMatEl( modelspace, bra, ket, J);
       double vserber = SerberMatEl( modelspace, bra, ket, J, params);
       Vserber.TwoBody.SetTBME(ch,ibra,iket,vserber);
       Vserber.TwoBody.SetTBME(ch,iket,ibra,vserber);
     }
    }
   }
   return Vserber;
 }



/// A single matrix element of the Serber potential
 double SerberMatEl( ModelSpace& modelspace, Ket& bra, Ket& ket, int J, const std::array<double,6>& params )
 {
   double oscillator_b2 = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
   double V0 = params[0];
   double mu = params[1];
   double A = params[2];
   double B = params[3];
   double C = params[4];
   double D = params[5];
   double kR = oscillator_b2/(mu*mu);
//   std::cout << "V0 " << V0 << " mu " << mu << "  A  " << A << "  B " << B << std::endl;

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

   double Vserber=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=std::abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
//       double isospin_factor = 1.0;
//       if ( (oa.tz2 == ob.tz2) and Sab==1 ) continue;
//       if ( oa.tz2 != ob.tz2) isospin_factor = 0.5 * (oa.tz2 * oc.tz2);  // <pn|V|pn> = 1/2 (V(T=1) + V(T=0)).  <pn|V|np> = -<pn|V|pn>
       if ( std::abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = AngMom::NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (std::abs(njab) <1e-7) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = AngMom::NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (std::abs(njcd) <1e-7) continue;

       double sig_dot_sig = 4*Sab - 3;

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
              //  This is achieved in a one-liner (  |tza+tzc| + |tza+tzd|*(-1)**(lab +Sab)  ), which gives either 1 or 0. The /2 is because tz2 is 2*tz.
              int asymm_factor = (std::abs(bra.op->tz2+ket.op->tz2) + std::abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // V conserves lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (std::abs(mosh_ab)<1e-8) continue;

              int tab = (Sab + lam_ab + 1)%2; // e.g. deuteron has S=1,L=0,T=0, S+L+T must be odd, so T = (L+S+1)%2
              if ( tab*2 < std::abs( oa.tz2+ob.tz2) ) continue; // If T < Tz, then this can't contribute.
              double tau_dot_tau = 4*tab - 3; // 1 if T=1, -3 if T=0.

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
                double RadialInt = 0;
                for (int p=pmin; p<=pmax; p++)
                {
                   double Bp = TalmiB(n_ab,lam_ab,n_cd,lam_cd,p);  // It's probably inefficient doing this in the inner loop, but oh well...
                   double Ip = 1.0 / pow(1.0+kR, p+1.5);  // Talmi integral of a Gaussian
                   RadialInt += Bp * Ip;
                }
                vrel = RadialInt * V0 * ( A + B*sig_dot_sig + C* tau_dot_tau + D*sig_dot_sig*tau_dot_tau );
                Vserber += vrel * njab * njcd * mosh_ab * mosh_cd * asymm_factor;


              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   return Vserber ;

 }


/// Direct implementation of equations (8.55), (8.62) and (8.65) in Suhonen.
 Operator SurfaceDeltaInteraction( ModelSpace& modelspace, double V0, double R)
 {

   Operator Vsdi(modelspace,0,0,0,2);

   int nchan = modelspace.GetNumberTwoBodyChannels();
   double hw = modelspace.GetHbarOmega();

//   modelspace.PreCalculateMoshinsky();
//   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
    TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;++ibra)
    {
     Ket& bra = tbc.GetKet(ibra);
     int a = bra.p;
     int b = bra.q;
     int na = bra.op->n;
     int la = bra.op->l;
     int j2a = bra.op->j2;
     int nb = bra.oq->n;
     int lb = bra.oq->l;
     int j2b = bra.oq->j2;
     double psi_a = HO_Radial_psi(na,la,hw,R);
     double psi_b = HO_Radial_psi(nb,lb,hw,R);
     for (int iket=ibra;iket<nkets;iket++)
     {
       Ket& ket = tbc.GetKet(iket);
       int c = ket.p;
       int d = ket.q;
       int nc = ket.op->n;
       int lc = ket.op->l;
       int j2c = ket.op->j2;
       int nd = ket.oq->n;
       int ld = ket.oq->l;
       int j2d = ket.oq->j2;
       double psi_c = HO_Radial_psi(nc,lc,hw,R);
       double psi_d = HO_Radial_psi(nd,ld,hw,R);

       double kappa_ac =  psi_a * psi_c * R;
       double kappa_bd =  psi_b * psi_d * R;

       double Kabcd = -V0 * kappa_ac * kappa_bd / (16*PI);
//       std::cout << " abcd " << a << " " << b << " " << c << " " << d << "   Kabcd " << Kabcd << std::endl;
       double AT = 1;
       Kabcd = -1./4 * AT * AngMom::phase(na+nb+nc+nd);

       double threej0ab = AngMom::ThreeJ(0.5*j2a, 0.5*j2b, J, 0.5,-0.5,0);
       double threej1ab = AngMom::ThreeJ(0.5*j2a, 0.5*j2b, J, 0.5, 0.5,-1);
       double threej0cd = AngMom::ThreeJ(0.5*j2c, 0.5*j2d, J, 0.5,-0.5,0);
       double threej1cd = AngMom::ThreeJ(0.5*j2c, 0.5*j2d, J, 0.5, 0.5,-1);

       double threej0dc = AngMom::ThreeJ(0.5*j2d, 0.5*j2c, J, 0.5,-0.5,0);
       double threej1dc = AngMom::ThreeJ(0.5*j2d, 0.5*j2c, J, 0.5, 0.5,-1);

       double hatfactor = sqrt((j2a+1)*(j2b+1)*(j2c+1)*(j2d+1));

       // Factor of 2 relative to Suhonen (8.65) because we just enforce parity conservation
       double Aabcd = 2*Kabcd *hatfactor* ( threej1ab*threej1cd  - AngMom::phase((2*la+2*lc+j2b+j2d)/2) * threej0ab*threej0cd);
       double Aabdc = 2*Kabcd *hatfactor* ( threej1ab*threej1dc  - AngMom::phase((2*la+2*ld+j2b+j2c)/2) * threej0ab*threej0dc);


       double vsdi = Aabcd;

       if (tbc.Tz !=0) // like nucleons
       {
	  vsdi -= AngMom::phase((j2c+j2d+2*J)/2) * Aabdc; // Suhonen (8.55)
       }
       else if ( bra.op->tz2 != ket.op->tz2 ) // if we have pnnp or nppn
       {
	  vsdi = AngMom::phase( (j2c + j2d + 2*J)/2 +1 ) * Aabdc; // Suhonen (8.56
       }

       if (a==b) vsdi /= SQRT2;
       if (c==d) vsdi /= SQRT2;

       Vsdi.TwoBody.SetTBME(ch,ibra,iket,vsdi);
       Vsdi.TwoBody.SetTBME(ch,iket,ibra,vsdi);
     }
    }
   }

    return Vsdi;
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
//    if (Lambda==0) // scalar => no sixJs, tbmes are not reduced.
    if (not op1.IsReduced() ) // scalar => no sixJs, tbmes are not reduced.
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



 

}// namespace imsrg_util
