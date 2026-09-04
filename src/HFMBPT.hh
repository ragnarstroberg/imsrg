///////////////////////////////////////////////////////////////////////////////////
//    HFMBPT.hh, part of  imsrg++
//    Copyright (C) 2019 Takayuki Miyagi
//      -- code modified by Ragnar Stroberg 2019 
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////
#ifndef NaturalOrbital_h
#define NaturalOrbital_h

#include "HartreeFock.hh"
#include <armadillo>
#include <map>
//#define OCC_CUT 1e-6

typedef unsigned long long int index_t;
class HFMBPT : public HartreeFock
{
  public:
    arma::mat C_HO2NAT; // transforamtion coefficients, 1st index ho, 2nd index NAT basis
    arma::mat C_HF2NAT; // transforamtion coefficients, 1st index hf, 2nd index NAT basis
    arma::vec Occ; // Occupation number

    bool use_NAT_occupations; // Option to use occupations from (if true) density matrix, or (if false) use naive filling.
    std::string NAT_order; // The default is to order by occupation

    ~HFMBPT();
    HFMBPT(Operator& hbare); // same as HartreeFock constructor
    void GetNaturalOrbitals();
    void GetDensityMatrix();
    void DensityMatrixPP(Operator& H);
    void DensityMatrixHH(Operator& H);
    void DensityMatrixPH(Operator& H);
    void DiagonalizeRho();
    void PrintOccupation();
    Operator TransformHFToNATBasis(Operator& OpIn);
    Operator TransformHOToNATBasis(Operator& OpIn);
    Operator GetNormalOrderedHNAT(int particle_rank=2);
//    Operator GetNormalOrderedHNAT();
//    double GetTransformed3bme( int Jab, int Jde, int J2, size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);
    double GetTransformed3bme( Operator& OpIn, int Jab, int Jde, int J2, size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);

    void PrintSPEandWF(); // Function override, since we want to express the SPWF in terms of HO states
    void ReorderHFMBPTCoefficients();
    void UseNATOccupations( bool tf=true ){ use_NAT_occupations=tf;}; // Choose whether to use occupations from rho.
    void OrderNATBy( std::string order ){ NAT_order = order;}; // Choose how to label orbits ("occupation", "energy", "mp2")

    double GetMP2_Energy(const Operator& H) const;
    double GetMP3_pp(const Operator& H) const;
    double GetMP3_hh(const Operator& H) const;
    double GetMP3_ph(const Operator& H)const ;
    double GetMP3_Energy(const Operator& H)const ;
    double GetMP4_Diagram(int id, const Operator& H) const ;

    arma::vec GetMP2_Impacts(Operator& OpIn) const;


    double GetDenom(const Operator& H, const std::vector<index_t>& holes, const std::vector<index_t>& particles) const;
    double GetMP4_F1( const Operator& H) const ;
    double GetMP4_F2( const Operator& H) const ;
    double GetMP4_F3( const Operator& H) const ;
    double GetMP4_F4( const Operator& H) const ;
    double GetMP4_F5( const Operator& H) const ;
    double GetMP4_F6( const Operator& H) const ;
    double GetMP4_F7( const Operator& H) const ;
    double GetMP4_F8( const Operator& H) const ;
    double GetMP4_F9( const Operator& H) const ;
    double GetMP4_F10( const Operator& H) const ;
    double GetMP4_F11( const Operator& H) const ;
    double GetMP4_F12( const Operator& H) const ;
    double GetMP4_F13( const Operator& H) const ;
    double GetMP4_F14( const Operator& H) const ;
    double GetMP4_F15( const Operator& H) const ;
    double GetMP4_F16( const Operator& H) const ;
    double GetMP4_F17( const Operator& H) const ;
    double GetMP4_F18( const Operator& H) const ;
    double GetMP4_F19( const Operator& H) const ;
    double GetMP4_F20( const Operator& H) const ;
    double GetMP4_F21( const Operator& H) const ;
    double GetMP4_F22( const Operator& H) const ;
    double GetMP4_F23( const Operator& H) const ;
    double GetMP4_F24( const Operator& H) const ;
    double GetMP4_F25( const Operator& H) const ;
    double GetMP4_F26( const Operator& H) const ;
    double GetMP4_F27( const Operator& H) const ;
    double GetMP4_F28( const Operator& H) const ;
    double GetMP4_F29( const Operator& H) const ;
    double GetMP4_F30( const Operator& H) const ;
    double GetMP4_F31( const Operator& H) const ;
    double GetMP4_F32( const Operator& H) const ;
    double GetMP4_F33( const Operator& H) const ;
    double GetMP4_F34( const Operator& H) const ;
    double GetMP4_F35( const Operator& H) const ;
    double GetMP4_F36( const Operator& H) const ;
    double GetMP4_F37( const Operator& H) const ;
    double GetMP4_F38( const Operator& H) const ;
    double GetMP4_F39( const Operator& H) const ;

};
#endif
