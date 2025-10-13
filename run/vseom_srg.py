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
smax_valence = 100   # limit of s for second stage of decoupling

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
ips1 = ms.GetOrbitIndex(0,0,1,-1)
ins1 = ms.GetOrbitIndex(0,0,1,+1)
HNO.SetOneBody(ips1,ips1, HNO.GetOneBody(ips1,ips1)-10)
HNO.SetOneBody(ins1,ins1, HNO.GetOneBody(ins1,ins1)-10)
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
#rw.WriteTokyo( HNO, valence_fname,'')
Hs = imsrgsolver.GetH_s()

#rw.WriteTokyo( Hs, 'sd.snt','')
#Hs = Hs.DoNormalOrderingCore()

### Write out the effective valence space interaction

Hs.ZeroBody=0.



cm=Commutator
cm.SetIMSRG3Noqqq(True)
gm=Generator()

## read in the reduce transition matrix for reference 
tdm_op=read_tdm("he6.ref",ms)


unt = UnitTest(ms)
rank_j, parity, rank_Tz, particle_rank, herm= 0,0,0,2,1
 
h3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)

chi= gm.GetVSEOM_ladder(h3,0)

nm=Norm_vs_new2(chi,chi,tdm_op)

print('Norm of random Op: ', nm)
chi=chi/np.sqrt(nm)

nm=Norm_vs_new2(chi,chi,tdm_op)

gm.force_decouple(Hs)

print('compute reference energy to check hermitian')
nm=gm.GetVSEOM_Overlap_rd(Hs,tdm_op)
#nm=gm.GetVSEOM_Overlap(Hs)
print('reference energy: ', nm, 'vs [ -4.19234292  9.75174484 ] in nmax2 for he8_2',)


Hs=Hs*0.01


##e,v1,v2=lanczos_proc(htc_vs, Norm_vs_new, Hs, chi, 160, 2,ms)
#e,vs,v2=arnoldi_proc(htc_vs, Norm_vs_new2,Norm4, Hs, chi, 88,2,ms,tdm_op)
#unt.SetRandomSeed(1)
#h3= unt.RandomOp( ms, rank_j,  rank_Tz, parity, particle_rank,herm)
#chi= gm.GetVSEOM_ladder(h3,0)
#nm=Norm_vs_new2(chi,chi,tdm_op)
#chi=chi/np.sqrt(nm)
#
#runs=1
#for ard_step in range(runs):
#    e,vs,lvs=arnoldi_proc(htc_vs, Norm_vs_new2,Norm4, Hs, chi, 40,2,ms,tdm_op)
#    vec_a = lvs[0]*0.
#    for ii in range(2):
#        for jj in range(len(e)):
#            vec_a += lvs[jj]*vs[jj,ii]
#    chi=vec_a/np.sqrt(Norm_vs_new2(vec_a,vec_a,tdm_op))

#
#
#        
#dims= len(ard_vectors)
#
#unmat=np.zeros([dims,dims])
#uhmat=np.zeros([dims,dims])
#
#for ii in range(len(ard_vectors)):
#
#    for jj in range(len(ard_vectors)):
#        if(jj > ii):
#            continue
#        unmat[ii,jj]=Norm_vs_new(ard_vectors[ii],ard_vectors[jj],tdm_op)
#        unmat[jj,ii]=unmat[ii,jj]
#
#
#
##unmat[np.abs(unmat) < 0.0000000001] = 0
#
#enn, vnn = np.linalg.eig(unmat)
#print(enn)
#
#vn_trunc=[]
#
#for i, ei in enumerate(enn):
##   if(ei.real > 0.0001):
#    vn_trunc.append(vnn[:,i]/np.sqrt(ei.real))
#
#
#
#norm_vectors=[]
#
#for ii in range(len(vn_trunc)):
#    vec_a = ard_vectors[0]*0.
#    vni=vn_trunc[ii]
#    for jj in range(len(ard_vectors)):
#        vec_a += ard_vectors[jj] * vni[jj]
#    norm_vectors.append(vec_a)
#
#c= len(norm_vectors)
#norm_new=np.zeros([c,c])
#hmat_new=np.zeros([c,c])
#
#for i in range(c):
#    vk=norm_vectors[i]
#    for j in range(c):
#        if(j > i ):
#            continue
#        vl=norm_vectors[j]
#
#        vec_c=htc_vs(Hs, vl)
#
#        hmat_new[i,j]=Norm_vs_new(vk, vec_c,tdm_op)+Norm3(vk,vl,Hs,ms,tdm_op)
#        hmat_new[i,j]=hmat_new[j,i]
#
#        norm_new[i,j]=Norm_vs_new(vk, vl,tdm_op)
#        norm_new[i,j]=norm_new[j,i]
#
#ef, vf = np.linalg.eig(hmat_new)
#print(ef)
#
#print('new norm')
#print(norm_new)
#print('new haml')
#print(hmat_new)

#vs=print_op(chi,ms)

#vout=chi*0

#vout = htc_vs(Hs, chi)
#print(Norm_vs_new(chi,vout,tdm_op))

#vs=print_op(vout,ms)
#print(vs)

#print(Norm3(chi,chi,Hs,ms,tdm_op))

## Generate all configurations
## ppvv only J=0

cfs=[]
chiv=[]
jj=0
pp=0
tt=1
for jj in [0]:
    ch=ms.GetTwoBodyChannelIndex(jj,0,1)
    kcf=ms.GetTwoBodyChannel(ch)
    bras=kcf.GetKetIndex_qq()+kcf.GetKetIndex_qv()
    kets=kcf.GetKetIndex_vv()
    for ibra in bras:
        dbra=kcf.GetKet(ibra)
        for iket in kets:
            dket=kcf.GetKet(iket)
            cfs.append([dbra.p,dbra.q,dket.p,dket.q,jj,pp,tt])
            chiv.append(chi.GetTwoBody(ch,ch,ibra,iket))
### ppvh
nch=ms.GetNumberTwoBodyChannels()

for ich in range(nch):
    dch = ms.GetTwoBodyChannel(ich)
    jj=dch.J
    pp=dch.parity
    tt=dch.Tz
    ch=ich
    
#    ch=ms.GetTwoBodyChannelIndex(jj,pp,tt)
#    print(jj,pp,tt,ich,ch)
    kcf=ms.GetTwoBodyChannel(ich)
#    if(kcf.GetNumberKets() == 0):
#        continue
            #print(ch, kcf.GetNumberKets())
    bras=kcf.GetKetIndex_qq()+kcf.GetKetIndex_qv()+kcf.GetKetIndex_vv()
    kets=kcf.GetKetIndex_vc()
    for ibra in bras:
        dbra=kcf.GetKet(ibra)
        for iket in kets:
            dket=kcf.GetKet(iket)
        #    if(dket.q %2 == 0):
#                continue
            cfs.append([dbra.p,dbra.q,dket.p,dket.q,jj,pp,tt])
            chiv.append(chi.GetTwoBody(ch,ch,ibra,iket))
    print(jj,pp,tt,len(cfs))
dims=len(cfs)

print(dims)

##hmat=np.zeros([dims,dims])
#h3mat=np.zeros([dims,dims])
#nmat=np.zeros([dims,dims])
#print(cfs)
##print(chiv)
#for i,cfl in enumerate(cfs):
#    ql=Hs*0.
#    ql.SetTwoBody(cfl[4],cfl[5],cfl[6],cfl[4],cfl[5],cfl[6], cfl[0],cfl[1],cfl[2],cfl[3],1.)
#    print(i, "th row")
#    for j,cfr in enumerate(cfs):
#        qr=Hs*0.
#        qr.SetTwoBody(cfr[4],cfr[5],cfr[6],cfr[4],cfr[5],cfr[6], cfr[0],cfr[1],cfr[2],cfr[3],1.)
#        nmat[i,j] = Norm_vs_new2(ql,qr,tdm_op)
#        qv=Hs*0.
#        qv=htc_vs(Hs, qr)
#        h3mat[i,j]=Norm4(ql,qr,Hs,ms,tdm_op)
##
##
#np.save("nmat_he8_2_ppvv_reftest.npy", nmat)
#np.save("h3mat_he8_2_ppvv_reftest.npy",h3mat)
#
#hmatbare=[]
#for i,cfl in enumerate(cfs):
#    ql=Hs*0.
#    ql.SetTwoBody(cfl[4],cfl[5],cfl[6],cfl[4],cfl[5],cfl[6], cfl[0],cfl[1],cfl[2],cfl[3],1.)
#    qv=htc_vs(Hs, ql)
#    qout=print_op(qv,ms)
#    hmatbare.append(qout)
#hmatbare=np.array(hmatbare)
##
##print(hmatbare)
##
##
#np.save("hmat_new_he8_2_ppvv_reftest.npy", hmatbare)


