#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(2,'He6','He6')
ut = pyIMSRG.UnitTest(ms)

jx,tx,px,prx,hx = 0,0,0,3,-1
jy,ty,py,pry,hy = 0,0,0,3,+1
X = ut.RandomOp(ms, jx, tx, px, prx,hx)
Y = ut.RandomOp(ms, jy, ty, py, pry,hy)

#Y.EraseThreeBody()

for t in ['comm110ss','comm220ss','comm111ss','comm121ss','comm221ss','comm122ss','comm222_pp_hhss','comm222_phss']:
   pyIMSRG.Commutator.TurnOffTerm(t)
for t in ['comm330ss','comm131ss','comm231ss','comm331ss','comm132ss','comm232ss','comm223ss','comm133ss']:
   pyIMSRG.Commutator.TurnOnTerm(t)

pyIMSRG.Commutator.SetUseIMSRG3(True)
pyIMSRG.Commutator.SetUseIMSRG3N7(True)
passed = ut.TestCommutators(X,Y)


ms = pyIMSRG.ModelSpace(1,'He6','He6')
ut = pyIMSRG.UnitTest(ms)

jx,tx,px,prx,hx = 0,0,0,3,-1
jy,ty,py,pry,hy = 0,0,0,3,+1
X = ut.RandomOp(ms, jx, tx, px, prx,hx)
Y = ut.RandomOp(ms, jy, ty, py, pry,hy)


for t in ['comm330ss','comm131ss','comm231ss','comm331ss','comm132ss','comm232ss','comm223ss','comm133ss']:
   pyIMSRG.Commutator.TurnOffTerm(t)
for t in ['comm332_ppph_hhhpss','comm332_pphhss','comm233_pp_hhss','comm233_phss','comm333_ppp_hhhss','comm333_pph_hhpss']:
   pyIMSRG.Commutator.TurnOnTerm(t)
passed &= ut.TestCommutators(X,Y)

print('passed? ',passed)

exit(not passed)
