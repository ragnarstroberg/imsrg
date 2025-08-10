#!/usr/bin/env python3
from pyIMSRG import *
import numpy as np
from lanczos import *



emax =2         # maximum number of oscillator quanta in the model space
ref = 'He4'     # reference used for normal ordering
val = 'p-shell' # valence space

core_generator = 'atan'   # definition of generator eta for decoupling the core (could also use 'white')
valence_generator = 'shell-model-atan'  # definition of generator for decoupling the valence space (could also use 'shell-model-white'
smax_core = 50       # limit of integration in flow parameter s for first stage of decoupling
smax_valence = 50   # limit of s for second stage of decoupling

#### Example format of how to read input interaction matrix elements from file (these are not included with the code)
#f2b = 'input/chi2b_srg0800_eMax16_EMax16_hwHO020.me2j.gz'
#f2e1,f2e2,f2l = 16,16,16
#f3b = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax12_EMax12_hwHO020.me3j.gz'
#f3e1,f3e2,f3e3 = 12,24,12
#LECs = 'srg0800'

f2b='input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz'
f2e1,f2e2,f2l = 14,28,14
f3b='input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin'
f3e1,f3e2,f3e3 = 18,36,18
hw=16
mode3n = 'no2b'
LECs = 'EM1820'

#### Otherwise, we use the Minnesota NN potential
#LECs = 'Minnesota'
#f3b = 'none'
#hw = 20    # harmonic oscillator basis frequency

### name of file to write resulting shell model effective interaction.
### *.snt is the exension used with KSHELL
valence_fname = 'output/{}_{}_{}_e{}_hw{}.snt'.format(val,ref,LECs,emax,hw)


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
H = Operator(ms,rank_j, parity, rank_Tz, particle_rank)


### Either generate the matrix elements of the Minnesota potential, or read in matrix elements from file
if LECs == 'Minnesota':
    H += OperatorFromString(ms,'VMinnesota')


else:
  ### Read Two-body matrix elements
  rw.ReadBareTBME_Darmstadt(f2b,H,f2e1,f2e2,f2l)
  ### Read Three-body matrix elements
  if f3b != 'none':
     if mode3n == 'no2b':
        H.ThreeBody.SetMode('no2b')
        H.ThreeBody.ReadFile([f3b],[f3e1,f3e2,f3e3])
     else:
        rw.Read_Darmstadt_3body(f3b,H,f3e1,f3e2,f3e3)


### Add the relative kinetic energy, so H = Trel + V
H += OperatorFromString(ms,'Trel')

### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H)
hf.Solve()
hf.PrintSPEandWF()

### Do normal ordering with respect to the HF basis
HNO = hf.GetNormalOrderedH(2)

### Create an instance of the IMSRGSolver class, used for solving the IMSRG flow equations
imsrgsolver = IMSRGSolver(HNO)
imsrgsolver.SetMethod('magnus')  # Solve using the Magnus formulation. Could also be 'flow_RK4'

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


cm=Commutator
gm=Generator()


unt = UnitTest(ms)
rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1
 
h3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)

chi= gm.GetVSEOM_ladder(h3,0)
nm=Norm_vs(chi,chi)
chi=chi/np.sqrt(nm)
nm=Norm_vs(chi,chi)

gm.force_decouple(Hs)

#lvs=[]
#print('initial vector')
info=1
gm.PrintEOMVS_ladder(Hs,info);
Hs.PrintOneBody()
#
#
vec2=chi*0
vec2.SetAntiHermitian()
#cm.comm110ss(Hs, chi, vec2)
#cm.comm220ss(Hs, chi, vec2)
#cm.comm111ss(Hs, chi, vec2)
#cm.comm121ss(Hs, chi, vec2)
#cm.comm122ss(Hs, chi, vec2)
#cm.comm222_pp_hhss(Hs, chi, vec2)
#cm.comm221ss(Hs, chi, vec2)
#cm.comm222_phss(Hs, chi, vec2)
#vec2=cm.CommutatorScalarScalar(Hs, chi)
#nm=gm.GetVSEOM_Overlap(vec2)
#print(nm)

#vec2=htc_vs(Hs, chi)
#vec2=gm.GetVSEOM_ladder(vec2,0)
#gm.PrintEOMVS_ladder(vec2,info)
#
#ai=Norm_vs(chi,vec2)
#print(ai)


#v1=htc_vs(Hs,chi)

#a=Norm_vs(chi,v1)
#
#v1=v1-a*chi
#
#bb=Norm_vs(v1,v1)
#
#v1=v1/np.sqrt(bb)

## methd a

#w1=htc_vs(Hs,v1)
#w2=htc_vs(Hs,chi)
#print('Norm condition: ', Norm_vs(v1,chi), Norm_vs(chi,v1))
#print(Norm_vs(w1,chi))
#print(Norm_vs(w2,v1))


#e,v1,v2=lanczos_proc(htc_vs, Norm_vs, Hs, chi, 30, 2)
e,v1,v2=arnoldi_proc(htc_vs, Norm_vs_new, Hs, chi, 80, 5)

dim=len(v2)
fnm=np.zeros([dim,dim])

for i ,vi in enumerate(v2):
    for j ,vj in enumerate(v2):
        nm=Norm_vs_new(vi,vj)
        if(abs(nm) >0.000001):
            print(i,j, nm)

#va=v2[0]
#vb=v2[1]
#
#print(Norm_vs_new(va,vb))
#print(Norm_vs_new(vb,va))
#
#Hs.EraseOneBody()
#vs1=htc_vs(Hs, va)
#print(Norm_vs_new(vb,vs1))
#print(Norm_vs_new(vs1,vb))

#vs1=htc_vs(Hs, vb)
#print(Norm_vs_new(va,vs1))
