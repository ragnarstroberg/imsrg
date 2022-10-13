/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                  ____                                         ///
///        _________________           _____________/   /\               _________________        ///
///       /____/_____/_____/|         /____/_____/ /___/  \             /____/_____/_____/|       ///
///      /____/_____/__G_ /||        /____/_____/|/   /\  /\           /____/_____/____ /||       ///
///     /____/_____/__+__/|||       /____/_____/|/ G /  \/  \         /____/_____/_____/|||       ///
///    |     |     |     ||||      |     |     |/___/   /\  /\       |     |     |     ||||       ///
///    |  I  |  M  |     ||/|      |  I  |  M  /   /\  /  \/  \      |  I  |  M  |     ||/|       ///
///    |_____|_____|_____|/||      |_____|____/ + /  \/   /\  /      |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |   /___/   /\  /  \/       |     |     |     ||||       ///
///    |  S  |  R  |     ||/|      |  S  |   \   \  /  \/   /        |  S  |  R  |  G  ||/|       ///
///    |_____|_____|_____|/||      |_____|____\ __\/   /\  /         |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |     \   \  /  \/          |     |     |     ||||       ///
///    |     |  +  |     ||/       |     |  +  |\ __\/   /           |     |  +  |  +  ||/        ///
///    |_____|_____|_____|/        |_____|_____|/\   \  /            |_____|_____|_____|/         ///
///                                               \___\/                                          ///
///                                                                                               ///
///           imsrg++ : Interface for performing standard IMSRG calculations.                     ///
///                     Usage is imsrg++  option1=value1 option2=value2 ...                       ///
///                     To get a list of options, type imsrg++ help                               ///
///                                                                                               ///
///                                                      - Ragnar Stroberg 2016                   ///
///                                                                                               ///
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
//    imsrg++.cc, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
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


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <omp.h>
#include "../IMSRG.hh"
#include "../Parameters.hh"
#include "../PhysicalConstants.hh"
#include "../version.hh"

struct OpFromFile {
   std::string file2name,file3name,opname;
   int j,p,t,r; // J rank, parity, dTz, particle rank
};

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
  std::cout << "######  imsrg++ build version: " << version::BuildVersion() << std::endl;

  Parameters parameters(argc,argv);
  if (parameters.help_mode) return 0;

  int emax = parameters.i("emax_imsrg");
  int emax_3body = parameters.i("emax_3body_imsrg");
  int e3max = parameters.i("e3max_imsrg");
  std::string ref = parameters.s("reference");
  // We coopt emax for number of commutators here
  int num_comms = parameters.i("emax");

  ModelSpace ms(emax, emax_3body, ref, ref);
  ms.SetE3max(e3max);

  Operator op1 = Operator(ms, 0, 0, 0, 3);
  op1.SetAntiHermitian();
  op1.ThreeBody.SwitchToPN_and_discard();
  Operator op2 = Operator(ms, 0, 0, 0, 3);
  op2.SetHermitian();
  op2.ThreeBody.SwitchToPN_and_discard();
  Operator op3 = Operator(ms, 0, 0, 0, 3);
  op3.SetHermitian();
  op3.ThreeBody.SwitchToPN_and_discard();

  for (int i = 0; i < num_comms; i+=1){
    std::cout << "Evaluating commutator " << i << std::endl;
    Commutator::comm232ss_mh_optimized(op1, op2, op3);
    std::cout << "Evaluated commutator " << i << std::endl;
  }


  op1.PrintTimes();

  return 0;
}