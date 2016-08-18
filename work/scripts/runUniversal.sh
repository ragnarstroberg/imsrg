#!/bin/bash

exe=imsrg++
#exe=../compiled/Universal
#exe="valgrind --track-origins=yes ../compiled/Universal"
#exe="valgrind --track-origins=yes imsrg++"

#vnn=/itch/exch/BlockGen/me2j/chi2b_srg0800_eMax12_lMax10_hwHO020.me2j.gz
#v3n=/itch/exch/me3j/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO020.me3j.gz

#vnn=../../input/chi2b_srg0800_eMax12_lMax10_hwHO020.me2j.gz
#v3n=../../input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax12_EMax12_hwHO020.me3j.gz
vnn=input/chi2b_srg0200_eMax04_hwHO024.me2j
v3n=none

#vnn=/itch/exch/BlockGen/me2j/chi2b_srg0800_eMax12_lMax10_hwHO024.me2j.gz
#v3n=/itch/exch/me3j/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO024.me3j.gz

#vnn=/itch/exch/N2LOSAT/me2j/chi2bSAT_srg0800_eMax14_lMax10_hwHO022.me2j.gz
#v3n=/itch/exch/N2LOSAT/me3j/chi2b3bSAT_J7666_4503c1-112c3-393c4377cD08cE-0040_hwconv036_JTP731_srg0800ho40J_eMax14_EMax14_hwHO022.me3j.gz
#v3n=none

#vnn=../../input/chi2b_srg0625_eMax12_hwHO020.me2j.gz
#v3n=../../input/jsTNF_Nmax_18_J12max_8_hbarOmega_20.00_Fit_cutoff_2.00_nexp_4_c1_-0.81_c3_-3.20_c4_5.40_cD_1.27_cE_-0.13_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_id_1_new_ant_E3_12_e_12.h5

BetaCM=0.0

#for ((A=14;A<=28;A++)); do
#for ((A=6;A<=6;A++)); do
A=22
#for reference in He6 C20 N29 Ti44 ; do
#for reference in vacuum ; do
for reference in O22 ; do

hw=20
#valence_space=O$A
#valence_space=sd-shell
#valence_space=O16,p0d5,n0d5,p0d3,n0d3,p1s1,n1s1
valence_space="customsd"

#valence_space=Ca40
#valence_space=Si34
#valence_space=psd-shell
#reference=default
#reference=Si34
#reference=Si28
#reference=O$A
#reference=He4
#A=18
#smax=200
smax=2

emax=4
e3max=12
#method=HF
#method=MP3
method=magnus
#method=flow
#method=flow-omega
omega_norm_max=0.25
#flowfile=output/BCH_SingleRef_${}_e${emax}.dat
file3='file3e1max=12 file3e2max=28 file3e3max=12'
Operators=
#Operators=Rp2,Rm2,HCM,HCM_28
#Operators=Rp2,rho0.0,rho0.2,rho0.5,rho1.0,rho1.2,rho1.5,rho2.0,rho2.5,rho3.0,rho3.5,rho4.5,rho6.0
Operators=E2
#Operators=protonFBC1,protonFBC2,protonFBC3,protonFBC4,protonFBC5,protonFBC6,protonFBC7,protonFBC8,protonFBC9
scratch=SCRATCH
other_arguments="SPWF=p0s1,p0p3"
other_arguments=${other_arguments}" custom_valence_space=O16,p0d5,n0d5,p0d3,n0d3,p1s1,n1s1"
#scratch=

OMP_NUM_THREADS=1

$exe 2bme=${vnn} 3bme=${v3n} emax=${emax} e3max=${e3max} method=${method} valence_space=${valence_space} hw=${hw} smax=${smax} ${file3} omega_norm_max=${omega_norm_max} reference=${reference} Operators=${Operators} scratch=${scratch} A=${A} use_brueckner_bch=false $other_arguments


done

