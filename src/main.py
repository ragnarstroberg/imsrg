#!/usr/bin/python

from pyIMSRG import *

inputsps = "../input/He4_lmax6.sps"
inputtbme = "../input/vsrg2.0_n3lo_noEmpot_Gaute_emax8_hw24.int"
hw = 24.0
targetMass = 4

rw = ReadWrite()

#rw.ReadSettingsFile("settings.inp")

modelspace = rw.ReadModelSpace(inputsps)
modelspace.SetHbarOmega(hw)
modelspace.SetTargetMass(targetMass)

Hbare = Operator(modelspace)
Hbare.SetHermitian()

Hbare.CalculateKineticEnergy()

rw.SetCoMCorr(False)
rw.ReadBareTBME(inputtbme, Hbare)
Tcm = TCM_Op(modelspace)

Hbare -= Tcm

hf = HartreeFock(Hbare)

hf.Solve()
print "EHF = ",hf.EHF
Hhf = hf.TransformToHFBasis(Hbare)

HhfNO = Hhf.DoNormalOrdering()

imsrgsolver = IMSRGSolver(HhfNO)

imsrgsolver.SetGenerator("atan")

imsrgsolver.SetSmax(1.0)

imsrgsolver.Solve()


