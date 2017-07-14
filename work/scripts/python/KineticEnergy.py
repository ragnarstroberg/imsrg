#!/usr/bin/python

from pyIMSRG import *

emax=2
hw=1
valence_space='He4'

rw = ReadWrite()
modelspace = ModelSpace(emax,valence_space)
modelspace.SetHbarOmega(hw)

Trel = Trel_Op(modelspace)

rw.WriteOperatorHuman(Trel,'Trel_hw1.op')
