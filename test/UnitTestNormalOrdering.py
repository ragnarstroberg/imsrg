#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(2,'He6','He6')
ms.PreCalculateSixJ()
ut = pyIMSRG.UnitTest(ms)

print('scalar 3b',flush=True)
Rando = ut.RandomOp(ms,0,0,0,3,+1)
passed = ut.TestNormalOrdering(Rando)

print('lambda=1 2b',flush=True)
Rando = ut.RandomOp(ms,1,0,0,2,+1)
passed &= ut.TestNormalOrdering(Rando)

print('parity-changing 2b',flush=True)
Rando = ut.RandomOp(ms,0,0,1,2,+1)
passed &= ut.TestNormalOrdering(Rando)

print('charge-changing 2b',flush=True)
Rando = ut.RandomOp(ms,0,1,0,2,+1)
passed &= ut.TestNormalOrdering(Rando)


ms2 = pyIMSRG.ModelSpace(2,'C12','p-shell')
ut = pyIMSRG.UnitTest(ms2)
Rando = ut.RandomOp(ms2,0,0,0,3,+1)

Rv1 = 1.0 * Rando
Rv2 = 1.0 * Rando
Rv1 = Rv1.UndoNormalOrdering()
Rv1.ThreeBody = Rando.ThreeBody
Rv1 = Rv1.DoNormalOrderingCore()
Rv2 = Rv2.ReNormalOrderCore()

diff = Rv1 - Rv2
print('Norm of diff = ',diff.Norm())



print('passed? ',passed)

exit(not passed)
