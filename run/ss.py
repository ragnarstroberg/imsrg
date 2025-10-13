#!/usr/bin/env python3
from pyIMSRG import *
from lanczos import *

emax =2         # maximum number of oscillator quanta in the model space
ref = 'He4'     # reference used for normal ordering
val = 'p-shell' # valence space
hw=16
ms = ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)

