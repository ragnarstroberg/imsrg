#!/bin/bash

exe=$HOME/ragnar_imsrg/work/compiled/SingleRef
#exe='valgrind --tool=massif --detailed-freq=2 --time-unit=ms $HOME/ragnar_imsrg/work/compiled/SingleRef'
WORKDIR=$HOME/ragnar_imsrg/work/scripts
NTHREADS=12

#vnn='/itch/exch/N2LOSAT/me2j/chi2bSAT_srg0400_eMax14_lMax10_hwHO022.me2j.gz'
#v3n='/itch/exch/N2LOSAT/me3j/chi2b3bSAT_J7666_4503c1-112c3-393c4377cD08cE-0040_hwconv036_JTP731_srg0400ho40J_eMax14_EMax14_hwHO022.me3j.gz'
#vnn='/itch/exch/N2LOSAT/me2j/chi2bSAT_srg0800_eMax14_lMax10_hwHO022.me2j.gz'
#v3n='/itch/exch/N2LOSAT/me3j/chi2b3bSAT_J7666_4503c1-112c3-393c4377cD08cE-0040_hwconv036_JTP731_srg0800ho40J_eMax14_EMax14_hwHO022.me3j.gz'

#vnn=input/TBMENNLOsat_eMax12_lMax10_hwHO022.me2j
#v3n=input/v3trans_J3T3.int_NN3Nnnlosat_nu3_330_161615.22.me3j.gz
#vnn=/itch/exch/N2LOSAT/me2j/chi2bSAT_srg0000_eMax14_lMax10_hwHO022.me2j.gz
#v3n=/itch/exch/N2LOSAT/me3j/chi2b3bSAT_J7666_4503c1-112c3-393c4377cD08cE-0040_hwconv036_JT3Nfull73_srg0000ho0_eMax14_EMax14_hwHO022.me3j.bin


#vnn=/itch/exch/BlockGen/me2j/chi2b_srg0800_eMax12_lMax10_hwHO020.me2j.gz
#v3n=/itch/exch/me3j/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO020.me3j.gz

#vnn=/itch/exch/BlockGen/me2j/chi2b_srg0800_eMax12_lMax10_hwHO024.me2j.gz
#v3n=/itch/exch/me3j/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO024.me3j.gz

vnn=/itch/exch/BlockGen/me2j/chi2b_srg0625_eMax12_hwHO020.me2j.gz
v3n=/itch/exch/me3j/chi2b3b400cD-02cE0098_srg0625ho40C_eMax14_EMax14_hwHO020.me3j.gz
#v3n=none

#vnn=../../input/chi2b_srg0625_eMax12_hwHO020.me2j.gz
#v3n=../../input/jsTNF_Nmax_18_J12max_8_hbarOmega_20.00_Fit_cutoff_2.00_nexp_4_c1_-0.81_c3_-3.20_c4_5.40_cD_1.27_cE_-0.13_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_id_1_new_ant_E3_12_e_12.h5

hw=20
#hw=22
#hw=24
nucl=O24
#smax=20
smax=200
emax=4
#e3max=14
e3max=12
omega_norm_max=0.25
file3='file3e1max=14 file3e2max=28 file3e3max=14'
#file3='file3e1max=16 file3e2max=16 file3e3max=15'
method=magnus

for ((emax=4;emax<=12;emax+=4)); do
#for ((e3max=14;e3max<=14;e3max+=2)); do
#for nucl in O16 O22 O24 O28; do
for nucl in Si34 S32; do
for method in magnus; do
#for method in magnus restore_4th_order; do
#for method in  restore_4th_order; do

flowfile=$HOME/ragnar_imsrg/output/BCH_SR_hw${hw}_${nucl}_${method}_e${emax}_s${smax}.dat
#jobname=SR_${method}_sat400_e${emax}_E${e3max}
jobname=SR__${method}_${nucl}_e${emax}_E${e3max}
qsub -N ${jobname} -q batchmpi -d $PWD -l walltime=192:00:00 -l nodes=1:ppn=${NTHREADS} -l vmem=60gb -m ae -M sstroberg@triumf.ca -j oe -o pbslog/${jobname}.o.`date +"%g%m%d%H%M"` << END
cd $WORKDIR
export OMP_NUM_THREADS=${NTHREADS}
${exe} 2bme=${vnn} 3bme=${v3n} emax=${emax} e3max=${e3max} nucleus=${nucl} hw=${hw} flowfile=${flowfile} smax=${smax} ${file3} method=${method} omega_norm_max=${omega_norm_max}
END


done
done
done
#done
