// Copyright (c) 2022 Matthias Heinz
// 
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include <iostream>

#include "Commutator.hh"
#include "IMSRG.hh"
#include "UnitTest.hh"

int main(void) {
    int emax = 4;
    int emax_3body = 4;
    int e3max = 6;
    std::string ref = "He4";

  ModelSpace ms(emax, emax_3body, ref, ref);
  ms.SetE3max(e3max);
  // ms.PreCalculateSixJ();

  UnitTest ut(ms);

  std::cout << "SanityCheck() " << ut.SanityCheck() << "\n";

  Commutator::SetUseIMSRG3(true);
  Commutator::SetUseIMSRG3N7(true);

  std::cout << "Commutators() " << ut.TestCommutators() << "\n";

  // Operator A = ut.RandomOp(ms, 0, 0, 0, 3, 1);
  // ut.Test3BodyAntisymmetry(A);
  // std::cout << A.IsHermitian() << "\n";
  // Operator B = ut.RandomOp(ms, 0, 0, 0, 3, 1);
  // std::cout << B.IsHermitian() << "\n";
  // ut.Test3BodyAntisymmetry(B);

  //   std::vector<std::string> skiplist = {"allN9", "allN8"};
  // ut.TestCommutators3(A, B, skiplist);
}