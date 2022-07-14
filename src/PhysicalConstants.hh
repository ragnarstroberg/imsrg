
#ifndef PhysicalConstants_h
#define PhysicalConstants_h

#include <cmath>


// This is where we define all the physical and mathematical constants that we will use

namespace PhysConst {
// I take all my physical constants from Wikipedia, because I'm a professional...

const double HBARC              = 197.3269718  ;                  // reduced Planck constant * speed of light in MeV * fm
const double M_PROTON           = 938.2720813  ;                  // proton mass in MeV/c^2
const double M_NEUTRON          = 939.5654133  ;                  // neutron mass in MeV/c^2
const double M_ELECTRON         = 0.5109989461 ;                  // electron mass in MeV/c^2
const double M_NUCLEON          = (M_PROTON+M_NEUTRON) / 2.0 ;    // avarage nucleon mass in MeV/c^2
const double M_PION_CHARGED     = 139.57018 ;                     // charged pion mass in MeV/c^2
const double M_PION_NEUTRAL     = 134.9766 ;                      // neutron pion mass in MeV/c^2
const double NUCLEON_VECTOR_G   = 1.00 ;                          // nucleon vector g factor, which is 1 due to CVC
const double NUCLEON_AXIAL_G    = 1.27 ;                          // nucleon axial g factor, aka gA. PDG lists 1.2723(23). From Lattice 1.271(13)
const double PROTON_SPIN_G      =  5.585690569;                   // proton spin g factor for magnetic moment from PDG 2020
const double NEUTRON_SPIN_G     = -3.826085 ;                     // neutron spin g factor for magnetic moment from PDG 2020
const double ELECTRON_SPIN_G    = 2.002319;                       // electron spin g factor
const double ALPHA_FS           = 1.0 / 137.035999;               // fine structure constant
const double F_PI               = 92.2 ;                          // pion decay constant
const double HARTREE            = 27.21138602;                    // 1 Hartree in eV
//const double PROTON_RCH2	= 0.770;			  // proton charge radius squared. This corresponds to Rp=0.8775 from CODATA 2010
const double PROTON_RCH2	= 0.707;			  // proton charge radius squared. This corresponds to Rp=0.8409 from PDG 2020
//const double NEUTRON_RCH2	= -0.1149;			  // neutron charge radius squared. From Angeli and Marinova (not sure how they got it)
const double NEUTRON_RCH2	= -0.1161;			  // neutron charge radius squared. From PDG 2020
//const double DARWIN_FOLDY	= 0.033;			  // Darwin-Foldy relativistic correction to charge radii  = 3hbar^2/4 mp^2 c^2
const double DARWIN_FOLDY	= 0.033172;			  // Darwin-Foldy relativistic correction to charge radii  = 3hbar^2/4 mp^2 c^2

// Math constants
const double SQRT2    = sqrt(2.0);
const double INVSQRT2 = 1.0 / SQRT2;
const double LOG2     = log(2.0)  ;                           // natural logarithm of 2
const double PI       = 4.0*atan(1.0) ;
const double SQRTPI   = sqrt(PI) ;
//const double SQRT2 =  1.41421356237309505                  // square root of 2
//const double INVSQRT2 = 0.7071067811865                  // 1 / sqrt(2)
//const double PI     = 3.14159265358979324                  // pi, the ratio of circumference to diameter
//const double SQRTPI = 1.7724538509055160                  // square root of pi

// I can't remember why I originally wanted these to be long doubles. It sort of seems like overkill...
//const long double SQRT2 = 1.4142135623730950488L         // square root of 2
//const long double INVSQRT2 = 0.70710678118654752440L     // 1 / sqrt(2)

} // namespace PhysConst




#endif
