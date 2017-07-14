#!/usr/bin/python

from pyIMSRG import *

hw = 24.0
valence_space = 'sd-shell'
emax=2

rw = ReadWrite()
modelspace = ModelSpace(emax,valence_space)
modelspace.SetHbarOmega(hw)

E2bare = ElectricMultipoleOp(modelspace,2)

rw.WriteTensorOneBody('E2bare_hw%.0f_%s_1b.op'%(hw,valence_space),E2bare,'E2')



