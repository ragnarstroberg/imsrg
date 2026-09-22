#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(2,'He6','He6')
ut = pyIMSRG.UnitTest(ms)

jA,tzA,pA,hermA = 0,0,0,-1
jB,tzB,pB,hermB = 0,0,0,+1
A = ut.RandomOp(ms,jA,tzA,pA,2,hermA)
B = ut.RandomOp(ms,jB,tzB,pB,2,hermB)
if A.IsReduced():
  A.MakeNotReduced()
if B.IsReduced():
  B.MakeNotReduced()

passed = ut.TestFactorizedDoubleCommutators(A,B)

jA,tzA,pA,hermA = 0,0,0,-1
jB,tzB,pB,hermB = 0,0,0,-1
A = ut.RandomOp(ms,jA,tzA,pA,2,hermA)
B = ut.RandomOp(ms,jB,tzB,pB,2,hermB)
if A.IsReduced():
  A.MakeNotReduced()
if B.IsReduced():
  B.MakeNotReduced()

passed &= ut.TestFactorizedDoubleCommutators(A,B)

print('passed? ',passed)

exit(not passed)
