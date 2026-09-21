
#include "imsrg_util.hh"
#include "AngMom.hh"
//#include "Commutator.hh"
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


 

}// namespace imsrg_util
