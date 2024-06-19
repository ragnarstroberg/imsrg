///////////////////////////////////////////////////////////////////////////////////
//    RPA.hh, part of  imsrg++
//    Copyright (C) 2021  Ragnar Stroberg
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

#ifndef RPA_hh
#define RPA_hh 1

#include <armadillo>

#include "ModelSpace.hh"
#include "Operator.hh"


class RPA
{
 public:
  // Fields
  ModelSpace * modelspace; ///< Pointer to the associated modelspace

  Operator H;
  arma::mat A;  
  arma::mat B;  
  arma::mat X;
  arma::mat Y;
  arma::vec Energies;
  size_t channel;


  // Methods
  RPA(); ///< Default constructor
  RPA(ModelSpace& ms);
  RPA(Operator& H);

  void ConstructAMatrix(int J, int parity, int Tz, bool Isovector);
  void ConstructBMatrix(int J, int parity, int Tz, bool Isovector);

  void ConstructAMatrix_byIndex(size_t ich_CC, bool Isovector);
  void ConstructBMatrix_byIndex(size_t ich_CC, bool Isovector);
//  void SolveTDA();
  void SolveCP(); // core polarization, i.e. 1st order approximation of TDA
  void SolveTDA();
  void SolveRPA();

  double GetEgs();

  double TransitionToGroundState( Operator& OpIn, size_t mu );
  double PVCouplingEffectiveCharge( Operator& OpIn, size_t k, size_t l);

  arma::vec GetX(size_t i);
  arma::vec GetY(size_t i);
  arma::vec GetEnergies();

};





#endif
