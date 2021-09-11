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

    arma::vec GetMP2_Impacts(Operator& OpIn) const;

};
#endif
