#!/bin/bash

#vnn=/itch/exch/BlockGen/me2j/chi2b_srg0800_eMax12_lMax10_hwHO020.me2j.gz
#v3n=/itch/exch/me3j/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO020.me3j.gz

vnn=/itch/exch/N2LOSAT/me2j/chi2bSAT_srg0800_eMax14_lMax10_hwHO022.me2j.gz
v3n=/itch/exch/N2LOSAT/me3j/chi2b3bSAT_J7666_4503c1-112c3-393c4377cD08cE-0040_hwconv036_JTP731_srg0800ho40J_eMax14_EMax14_hwHO022.me3j.gz

#vnn=../../input/chi2b_srg0625_eMax12_hwHO020.me2j.gz
#v3n=../../input/jsTNF_Nmax_18_J12max_8_hbarOmega_20.00_Fit_cutoff_2.00_nexp_4_c1_-0.81_c3_-3.20_c4_5.40_cD_1.27_cE_-0.13_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_id_1_new_ant_E3_12_e_12.h5

hw=22
nucl=O16
smax=25
emax=4
e3max=4
method=magnus
#method=flow-omega
omega_norm_max=0.25
flowfile=output/BCH_SingleRefjs_O16_e${emax}.dat
file3='file3e1max=14 file3e2max=28 file3e3max=14'

../compiled/SingleRef 2bme=${vnn} 3bme=${v3n} emax=${emax} e3max=${e3max} method=${method} nucleus=${nucl} hw=${hw} flowfile=${flowfile} smax=${smax} ${file3} omega_norm_max=${omega_norm_max}


