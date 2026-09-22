#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(2,'He6','He6')
ut = pyIMSRG.UnitTest(ms)
#pyIMSRG.Commutator.SetUseIMSRG3(True)
#pyIMSRG.Commutator.SetUseIMSRG3N7(True)
jx,tx,px,prx,hx = 0,0,0,2,-1
jy,ty,py,pry,hy = 0,0,0,2,+1
X = ut.RandomOp(ms, jx, tx, px, prx,hx)
Y = ut.RandomOp(ms, jy, ty, py, pry,hy)
X.MakeNotReduced()
Y.MakeNotReduced()
passed = ut.TestCommutators(X,Y)


print('passed? ',passed)

exit(not passed)
