#!/bin/bash

#vnn=/itch/exch/BlockGen/me2j/chi2b_srg0800_eMax12_lMax10_hwHO020.me2j.gz
#v3n=/itch/exch/me3j/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO020.me3j.gz
vnn=/itch/exch/me2j/chi2b_srg0800_eMax12_lMax10_hwHO011.me2j.gz
v3n=/itch/exch/me3j/chi2b3b400c1-081c3-320c4540cD-02cE0098_hwconv036_srg0800ho40J_eMax14_EMax14_hwHO011.me3j.gz

#hw=20
hw=10.51
#valence_space=sd-shell
valence_space=fp-shell
smax=20
emax=5
e3max=8
delta=0.0
omega_norm_max=1.0
domega=0.25
flowfile=output/BCH_ValenceSpace_fp.dat
file3='file3e1max=14 file3e2max=28 file3e3max=14'
intfile=output/FP_e${emax}_d${delta}.int

../compiled/ValenceSpace 2bme=${vnn} 3bme=${v3n} emax=${emax} e3max=${e3max} valence_space=${valence_space} hw=${hw} flowfile=${flowfile} smax=${smax} ${file3} denominator_delta=${delta} intfile=${intfile} omega_norm_max=${omega_norm_max} domega=${domega}


