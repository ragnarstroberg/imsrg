#!/usr/bin/env python3


from pyIMSRG import *


emax=2
hw = 12 # This is a guess.
ref = 'vacuum'
val = 'psd-shell'
fpsdmwk = 'PSDMWK.snt' # or whatever it is...


ms = ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)

rw = ReadWrite()

for i in range(20):
    kappa = -0.2*(i-10)
    fout = 'PSDMWK_QQ_kappa{}.snt'.format(round(kappa,2))
    H = Operator(ms,0,0,0,2)
    rw.ReadTokyo(fpsdmwk, H)
    H += kappa * OperatorFromString(ms, 'VQQ')
    rw.WriteTokyo(H, fout, '')



