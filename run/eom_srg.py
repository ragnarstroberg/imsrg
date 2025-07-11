from pyIMSRG import *
import numpy as np

emax =3        # maximum number of oscillator quanta in the model space
ref = 'O16'     # reference used for normal ordering
val = ref # valence space

core_generator = 'atan'   # definition of generator eta for decoupling the core (could also use 'white')
smax_core = 10      # limit of integration in flow parameter s for first stage of decoupling
#smax_core = 0       # limit of integration in flow parameter s for first stage of decoupling

##### Example format of how to read input interaction matrix elements from file (these are not included with the code)
#f2b='input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz'
#f2e1,f2e2,f2l = 14,28,14
#f3b='input/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin'
#f3e1,f3e2,f3e3 = 18,36,18
#mode3n='no2b'
#LECs = 'EM1820'
#hw=16

#### Example format of how to read input interaction matrix elements from file (these are not included with the code)
f2b='input/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz'
f2e1,f2e2,f2l = 14,28,14
f3b='input/NO2B_ThBME_EM7.5_1.8_2.0_IS_hw16from16_ms14_28_18.me3j.gz'
f3e1,f3e2,f3e3 = 14,28,18
mode3n='no2b'
LECs = 'EM7.5_1820'
hw=16



#hw = 20    # harmonic oscillator basis frequency
#LECs='Minnesota'
#f3b='none'

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
print('after reading files, 3-body norm is',H.ThreeBodyNorm())

### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H)
hf.Solve()
hf.PrintSPEandWF()

### Do normal ordering with respect to the HF basis, and retain only up to 2-body operators
HNO = hf.GetNormalOrderedH(2)

### Create an instance of the IMSRGSolver class, used for solving the IMSRG flow equations
imsrgsolver = IMSRGSolver(HNO)
imsrgsolver.SetMethod('magnus')  # Solve using the Magnus formulation. Could also be 'flow_RK4'

imsrgsolver.SetGenerator(core_generator)
imsrgsolver.SetSmax(smax_core)

### Do the first stage of integration to decouple the core
imsrgsolver.Solve()


### Hs is the IMSRG-evolved Hamiltonian
Hs = imsrgsolver.GetH_s()


from  lanczos import *
cm=Commutator
gm=Generator()


ndim=10
unt = UnitTest(ms)
rank_j, parity, rank_Tz, particle_rank, herm= 2,0,0,2,0
h3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
#h3.MakeReduced()
chi= gm.GetEOM_ladder(h3,0)
cnorm=Norm(chi,chi)
chi=chi/np.sqrt(cnorm)
e,v,lvs = lanczos_proc( htc, Norm, Hs, chi, ndim)
e=np.sort(e)
print(e)
