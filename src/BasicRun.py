#!/usr/bin/env python3

from pyIMSRG import *
from math import pi, fmod
import numpy as np
#import matplotlib.pyplot as plt
#plt.style.use('seaborn-whitegrid')
from time import process_time,perf_counter
from multiprocessing import Pool

t1_start = perf_counter() #process_time()


#x = np.linspace(0, 10, 30)
#y = np.sin(x)

#plt.plot(x, y, 'o', color='black');
#plt.show()

#Smax = 501
#def parallel(smax):
emax = 4   # maximum number of oscillator quanta in the model space
ref = 'Be9'     # reference used for normal ordering
val = 'psd-shell' # valence space
A = 9

hw = 20    # harmonic oscillator basis frequency
core_generator = 'white'   # definition of generator eta for decoupling the core (could also use 'white')
valence_generator = 'shell-model'  # definition of generator for decoupling the valence space (could also use 'shell-model-white'
smax_core = 5       # limit of integration in flow parameter s for first stage of decoupling
smax_valence = 10   # limit of s for second stage of decoupling

#### Example format of how to read input interaction matrix elements from file (these are not included with the code)
#f2b = 'TwBME-HO_NN-only_N3LO_EM500_srg1000_hw16_emax4_e2max8.me2j.gz'
f2b = 'TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw20_emax14_e2max28.me2j.gz'
f2e1,f2e2,f2l = 14,28,14 #emax=spemaxtruncation, e2max=sumof2bspemaxtruncation, lmaxtruncationonl
f3b = 'NO2B_ThBME_srg2.0_ramp46_N3LO_EM500_JJmax13_c1_-0.81_c3_-3.2_c4_5.4_cD_0.7_cE_-0.06_LNL2_650_500_IS_hw20from30_ms14_28_18.stream.bin'
#f3b = 'ThBME_c1_-0.81_c3_-3.2_c4_5.4_cD_0.83_cE_-0.052_Local2_500_IS_hw16from16_ms4_8_12.me3j.gz'
#f3b = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax12_EMax12_hwHO020.me3j.gz'
f3e1,f3e2,f3e3 = 14,28,18
LECs = 'srg2.0'
format3n = 'no2b'


### Otherwise, we use the Minnesota NN potential
LECs = 'Minnesota'
f3b = 'none'


##########################################################################
###  END PARAMETER SETTING. BEGIN ACTUALLY DOING STUFF ##################
##########################################################################


### Create an instance of the ModelSpace class
ms = ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)

### the ReadWrite object handles reading and writing of files
rw = ReadWrite()

rank_j, parity, rank_Tz, particle_rank = 0,0,0,2
if f3b != 'none':
   particle_rank = 3

### Create an instance of the Operator class, representing the Hamiltonian
H0 = Operator(ms,rank_j, rank_Tz, parity, particle_rank)

mpi = (2*M_PION_CHARGED + M_PION_NEUTRAL)/3
a0 = mpi/(8*pi*M_NUCLEON)  #g0=1
a1 = -a0/2
a2 = a0

#VPT = UnitTest(ms).RandomOp(ms, 0, 0, 1, 2, 1)
#VPT = OperatorFromString( ms, "VPT_{}_0_0".format(a0))

#VPT = OperatorFromString( ms, "VPT_a0_0_0")

### Either generate the matrix elements of the Minnesota potential, or read in matrix elements from file
if LECs == 'Minnesota':
    H0 += OperatorFromString(ms,'VMinnesota')

else:
  rw.ReadBareTBME_Darmstadt(f2b,H0,f2e1,f2e2,f2l)
  if f3b != 'none':
     if format3n == 'no2b':
        H0.ThreeBody.SetMode(format3n)
        H0.ThreeBody.ReadFile( [f3b],  [f3e1,f3e2,f3e3] )
     else:
        rw.Read_Darmstadt_3body(f3b,H0,f3e1,f3e2,f3e3)

### Add the relative kinetic energy, so H = Trel + V
H0 += OperatorFromString(ms,'Trel')


### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H0)
hf.Solve()
hf.PrintSPEandWF()
#VPT = hf.TransformToHFBasis(VPT)
 
### Do normal ordering with respect to the HF basis
H0NO = hf.GetNormalOrderedH(2)
#VPT = VPT.DoNormalOrdering()
#ueta = UnitTest(ms).RandomOp(ms, 0, 0, 0, 1, -1 )
#uetapt=hf.TransformToHFBasis(UnitTest(ms).RandomOp(ms, 0,0,1,2,1))
#print('rank of ueta=',uetapt.GetParticleRank())


### Create an instance of the IMSRGSolver class, used for solving the IMSRG flow equations
imsrgsolver = IMSRGSolver(H0NO)
imsrgsolver.SetMethod('flow_RK4')  # Solve using the Magnus formulation. Could also be 'flow_RK4'

imsrgsolver.SetGenerator(core_generator)
imsrgsolver.SetSmax(smax_core)

### Do the first stage of integration to decouple the core
imsrgsolver.Solve()


### Now set the generator for the second stage to decouple the valence space
imsrgsolver.SetGenerator(valence_generator)
imsrgsolver.SetSmax(smax_valence)
imsrgsolver.Solve()

### Hs is the IMSRG-evolved Hamiltonian
Hs = imsrgsolver.GetH_s()

### Shell model codes assume the interaction is normal ordered with respect to the core
### and typically we choose a reference different from the core, so we need to re-normal-order
### with respect to the core
Hs = Hs.UndoNormalOrdering()

Hs = Hs.DoNormalOrderingCore()

### Write out the effective valence space interaction

valence_fnameH0 = 'output/{}_{}_{}_e{}_hw{}_H0.snt'.format(val,ref,LECs,emax,hw)
rw.WriteTokyo( H0, valence_fnameH0 ,'H0')
#rw.WriteOperatorHuman(VPTs,'VPT.op')
#rw.WriteTokyo(H0,'H0.snt','')




t1_stop = perf_counter() #process_time()
print("Elapsed time during the whole program in {} minutes".format((t1_stop-t1_start)/60))


