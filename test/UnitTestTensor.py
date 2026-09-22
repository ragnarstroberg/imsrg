#!/usr/bin/env python3

import pyIMSRG

emax=  1
ms = pyIMSRG.ModelSpace(emax,'He6','He6')
ut = pyIMSRG.UnitTest(ms)
passed = True


#pyIMSRG.Commutator.TurnOnTerm("comm331st")
#pyIMSRG.Commutator.TurnOnTerm("comm223st")
#pyIMSRG.Commutator.TurnOnTerm("comm231st")
#pyIMSRG.Commutator.TurnOnTerm("comm232st")
#pyIMSRG.Commutator.TurnOnTerm("comm133st")
#pyIMSRG.Commutator.TurnOnTerm("comm132st")
#pyIMSRG.Commutator.TurnOnTerm('comm233_pp_hhst')
#pyIMSRG.Commutator.TurnOnTerm('comm233_phst')

#pyIMSRG.Commutator.TurnOffTerm('comm222_phst')

Jx,px,Tx,rankx,hx = 0,0,0,3,+1
Jy,py,Ty,ranky,hy = 2,0,0,3,+1
X = ut.RandomOp(ms,Jx,Tx,px,rankx,hx)
Y = ut.RandomOp(ms,Jy,Ty,py,ranky,hy)
passed &= ut.TestCommutators_Tensor(X,Y)

#Jx,px,Tx,rankx,hx = 0,1,0,3,-1
#Jy,py,Ty,ranky,hy = 1,0,0,3,+1
#X = ut.RandomOp(ms,Jx,Tx,px,rankx,hx)
#Y = ut.RandomOp(ms,Jy,Ty,py,ranky,hy)
#passed &= ut.TestCommutators_Tensor(X,Y)
#
#Jx,px,Tx,rankx,hx = 0,0,0,3,-1
#Jy,py,Ty,ranky,hy = 2,0,0,3,+1
#X = ut.RandomOp(ms,Jx,Tx,px,rankx,hx)
#Y = ut.RandomOp(ms,Jy,Ty,py,ranky,hy)
#passed &= ut.TestCommutators_Tensor(X,Y)

print('passed? ',passed)

prof = pyIMSRG.IMSRGProfiler()
prof.PrintTimes()

exit(not passed)
