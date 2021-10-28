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


  // Methods
  RPA(); ///< Default constructor
  RPA(ModelSpace& ms);
  RPA(Operator& H);

  void ConstructAMatrix(int J, int parity, int Tz);
  void ConstructBMatrix(int J, int parity, int Tz);

  void ConstructAMatrix_byIndex(size_t ich_CC);
  void ConstructBMatrix_byIndex(size_t ich_CC);
//  void SolveTDA();
  arma::vec SolveTDA();
  arma::vec SolveRPA();

};





#endif
