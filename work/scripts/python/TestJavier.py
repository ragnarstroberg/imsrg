#!/usr/bin/python

from pyIMSRG import *

hw = 40.0
valence_space = 'sp-shell'
emax=2
StateFile = 'jhernandez/2B_ME_Instrinsic_States.txt'
MatelFile = 'jhernandez/2B_ME_MEC_NLO_intrinsic.txt'

rw = ReadWrite()
modelspace = ModelSpace(emax,valence_space)
modelspace.SetHbarOmega(hw)

M1bare = MagneticMultipoleOp(modelspace,1)

rw.ReadRelCMOpFromJavier(StateFile, MatelFile, M1bare)


#rw.WriteTensorOneBody('M1bare_hw%.0f_%s_1b.op'%(hw,valence_space),M1bare,'M1')
rw.WriteTensorTwoBody('M1bare_hw%.0f_%s_2b.op'%(hw,valence_space),M1bare,'M1')
rw.WriteOperatorHuman(M1bare,'M1bare.op')


