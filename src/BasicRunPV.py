#!/usr/bin/env python3

from errno import EEXIST
from os import makedirs,path,system
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
emax = 2   # maximum number of oscillator quanta in the model space
ref = 'Be9'     # reference used for normal ordering
val = 'psd-shell' # valence space
#val = 'C12,p0p1,p0d5,p1s1,n0p1,n0d5,n1s1'
A = 9

hw = 20    # harmonic oscillator basis frequency
core_generator = 'atan'   # definition of generator eta for decoupling the core (could also use 'white')
valence_generator = 'shell-model-atan'  # definition of generator for decoupling the valence space (could also use 'shell-model-white'

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
#Any file name enters with me2j or me3j is Darmstadt


### Otherwise, we use the Minnesota NN potential
#LECs = 'Minnesota'
#f3b = 'none'


##########################################################################
###  END PARAMETER SETTING. BEGIN ACTUALLY DOING STUFF ##################
##########################################################################


### Create an instance of the ModelSpace class
ms = ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)

"""
def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''
    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise 
#output_dir = "../Runs/{}".format(ref)
#mkdir_p(output_dir)
#Path=path('/Users/davidkekejian/imsrg/jobs/Runs/{}'.format(ref))
#system("cd output_dir")
#jobsfolder = open(output_dir, "w")
"""

### the ReadWrite object handles reading and writing of files
rw = ReadWrite()

rank_j, parity, rank_Tz, particle_rank = 0,0,0,2
if f3b != 'none':
   particle_rank = 3

### Create an instance of the Operator class, representing the Hamiltonian
H0 = Operator(ms,rank_j, rank_Tz, parity, particle_rank)

mpi = (2*M_PION_CHARGED + M_PION_NEUTRAL)/3
a0 = mpi*HBARC/(8*pi*M_NUCLEON)  #g0=1
a1 = -a0/2
a2 = a0

#VPT = UnitTest(ms).RandomOp(ms, 0, 0, 1, 2, 1)
VPT = OperatorFromString( ms, "VPT_{}_0_0".format(a0))

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

#H0 += -3 * OperatorFromString(ms,'LdotS')

### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H0)
hf.Solve()
hf.PrintSPEandWF()
VPT = hf.TransformToHFBasis(VPT)
 
### Do normal ordering with respect to the HF basis
H0 = hf.GetNormalOrderedH(2)
VPT = VPT.DoNormalOrdering()
#ueta = UnitTest(ms).RandomOp(ms, 0, 0, 0, 1, -1 )
#uetapt=hf.TransformToHFBasis(UnitTest(ms).RandomOp(ms, 0,0,1,2,1))
#print('rank of ueta=',uetapt.GetParticleRank())

RMS = 1.2*pow(A,1/3)
#H = Operator(ms,0,0,0,2)
#H += H0NO 
#H += VPT
#H.PrintOneBody()

#Same code as above without parallel counting
Schiff = SchiffOp(ms,3,1,RMS)
Schiff = hf.TransformToHFBasis(Schiff)
#print('parity=',Schiff.GetParity())
#print('particlerank=',Schiff.GetParticleRank())
#print('Jrank=',Schiff.GetJRank())

Schiffpp = Operator(ms,1,0,0,2)
Schiffpp = hf.TransformToHFBasis(Schiffpp)

#delta_s = 1
#Smax = 5000
#H0 = H0NO 
#VPTs = VPT #uetapt #VPT 
#H0.PrintOneBody()
#print('after')
#H0 += 5
#Schiffpp.PrintOneBody()


eta = Operator(ms, 0, 0, 0, 2)
etapv = Operator(ms, 0, 0, 1, 2)
eta.SetAntiHermitian()
etapv.SetAntiHermitian()
generator =  GeneratorPV()
#H0data = []


imsrgsolver = IMSRGSolverPV(H0,VPT,Schiff,Schiffpp)
imsrgsolver.SetDsmax(0.1)
imsrgsolver.SetSmax(5)
imsrgsolver.SetGeneratorPV(core_generator)
imsrgsolver.Solve_RK4()

imsrgsolver.SetSmax(10)
imsrgsolver.SetGeneratorPV(valence_generator)
imsrgsolver.Solve_RK4()

Hs = imsrgsolver.GetH_s()
VPTs = imsrgsolver.GetVPT_s()
Spps = imsrgsolver.GetSchiffpp_s()
Sp = imsrgsolver.GetSchiff_s()

#VPTs.PrintTwoBody()
#Hs.PrintOneBody()
#Hs.PrintTwoBody()

valence_fnameHs = '{}_{}_{}_e{}_hw{}_Hs.snt'.format(val,ref,LECs,emax,hw)
rw.WriteTokyo(Hs, valence_fnameHs, '')
valence_fnameVs = '{}_{}_{}_e{}_hw{}_Vs.snt'.format(val,ref,LECs,emax,hw)
valence_fnameSp = '{}_{}_{}_e{}_hw{}_Sp.snt'.format(val,ref,LECs,emax,hw)
#rw.WriteTensorTokyo(valence_fnameVs, VPTs)
#rw.WriteTokyo(VPTs, valence_fnameVs, '')
#rw.WriteTensorTokyo(valence_fnameVs, Spps)
rw.WriteTensorTokyo(valence_fnameVs, VPTs)
rw.WriteTensorTokyo(valence_fnameSp, Sp)



#jobsfolder.close()

### Finally, print out profiling information so we know why this took so dang long to run...
prof = IMSRGProfiler()
prof.PrintAll()



t1_stop = perf_counter() #process_time()
print("Elapsed time during the whole program in {} minutes".format((t1_stop-t1_start)/60))


