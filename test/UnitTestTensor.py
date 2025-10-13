#!/usr/bin/env python3

import pyIMSRG

emax=  2
ms = pyIMSRG.ModelSpace(emax,'He6','He6')
ut = pyIMSRG.UnitTest(ms)
passed = True


#pyIMSRG.Commutator.TurnOnTerm('comm223st')
#pyIMSRG.Commutator.TurnOnTerm('comm331st')

Jx,px,Tx,rankx,hx = 0,0,0,3,+1
Jy,py,Ty,ranky,hy = 1,0,0,3,+1
X = ut.RandomOp(ms,Jx,Tx,px,rankx,hx)
Y = ut.RandomOp(ms,Jy,Ty,py,ranky,hy)
passed &= ut.TestCommutators_Tensor(X,Y)

Jx,px,Tx,rankx,hx = 0,1,0,3,-1
Jy,py,Ty,ranky,hy = 1,0,0,3,+1
X = ut.RandomOp(ms,Jx,Tx,px,rankx,hx)
Y = ut.RandomOp(ms,Jy,Ty,py,ranky,hy)
passed &= ut.TestCommutators_Tensor(X,Y)

Jx,px,Tx,rankx,hx = 0,0,0,3,-1
Jy,py,Ty,ranky,hy = 2,0,0,3,+1
X = ut.RandomOp(ms,Jx,Tx,px,rankx,hx)
Y = ut.RandomOp(ms,Jy,Ty,py,ranky,hy)
passed &= ut.TestCommutators_Tensor(X,Y)

print('passed? ',passed)

exit(not passed)
