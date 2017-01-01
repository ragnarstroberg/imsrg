#!/usr/bin/python

from os import path,environ,mkdir,remove
from subprocess import call
from time import sleep

NTHREADS=24
#exe = '%s/Research/ragnar_imsrg/work/compiled/Universal'%(environ['HOME'])
exe = '%s/bin/imsrg++'%(environ['HOME'])
#exe = '/work/hda21/hda212/imsrg++' # trying Jason's version of the executable
#batch_mode=False
batch_mode=True

ELEM = ['n','H','He','Li','Be','B','C','N',
       'O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',
       'Ca','Sc','Ti','V','Cr','Mn','Fe','Co',  'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']# ,'Bi','Po','At','Rn','Fr','Ra','Ac','Th','U','Np','Pu']

ARGS  =  {}

#ARGS['2bme'] = 'input/chi2bSAT_srg0000_eMax14_lMax10_hwHO022.me2j.gz'
#ARGS['3bme'] = 'input/chi2b3bSAT_J7666_4503c1-112c3-393c4377cD08cE-0040_hwconv036_JT3Nfull73_srg0400ho40J_eMax14_EMax14_hwHO022.me3j.gz'
#ARGS['3bme'] = 'input/chi2b3bSAT_J7666_hwconv022_JT3Nfull73_srg0000ho40C_eMax14_EMax14_hwHO022.me3j.gz'
#ARGS['hw'] = '22'
#ARGS['LECs'] = 'SAT'
#ARGS['LECs'] = 'srg0625'

ARGS['smax'] = '500'
ARGS['emax'] = '14'
#ARGS['e3max'] = '16'
ARGS['e3max'] = '14'
#ARGS['lmax3'] = '10' # for comparing with Heiko
ARGS['omega_norm_max'] = '0.25'
ARGS['ode_tolerance'] = '1e-5'
#ARGS['valence_space'] = 'sd-shell'
#ARGS['reference'] = 'O22'
#ARGS['A'] = '40'
#ARGS['nsteps'] = '2'
ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=10'
ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=14'
#ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
#ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=16'
#ARGS['scratch'] = 'SCRATCH' JDH

#ARGS['method'] = 'MP3'
ARGS['method'] = 'magnus'
#ARGS['method'] = 'brueckner'
#ARGS['method'] = 'flow'
#ARGS['method'] = 'HF'

FILECONTENT = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#SBATCH --output=slurm_log/%s.%%j
#SBATCH --time=%s
#SBATCH --mail-user=sstroberg@triumf.ca
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
echo NTHREADS = %d
export OMP_NUM_THREADS=%d
time srun %s
"""


if not path.exists('slurm_log'): mkdir('slurm_log')

# for A in range(44,64):
#for A in [46,47,48,49,53]:
#for A in [64,65,66,67,68,69,70]:
#for A in range(68,80):
#for A in [48]:
#for Z in range(8,20):
# for A in range( 2*Z-2, 2*Z+3 ):
#for Z in range(6,7):
#for betaCM in [0,1,3]:
for A in [8]:
#for A in [74,76,78,80]:
 Z=4
# ARGS['BetaCM'] = '%.1f'%betaCM
# for etacrit in ['1e-6','1e-8']:
#  ARGS['eta_criterion'] = etacrit
# Z=32
#  if (A-Z)<8 or (A-Z)>20: continue
#for A in range(51,57):
#for A in range(10,11):
#  e = 12
# for e in [4,6,8,10,12]:
# for e in [8,10,12]:
 for e in [10]:
# for e in [8]:
# for e in [8,10]:
# for e in [4]:
  for hw in [12]:
#  for hw in [16,24]:
#  for hw in [12,16,20,24,28]:
#  for hw in [24]:
#   if e==8 and A==24: continue
#   if  e==12 and hw==24: continue
#  for hw in [12,16,20,24,28]
# for e in [4,10,12]:
#  for method in ['magnus','flow']:
#   for reference in ['F%d'%(A)]:
#   for reference in ['Ca40', 'V%d'%(A)]:
#   for reference in ['Ca40', 'V%d'%(A)]:
#   for reference in ['Ge%d'%(A),'Se%d'%(A)]:
#   for reference in ['He%d'%(10)]:
#   for reference in ['Ge%d'%(A),'Se%d'%(A)]:
   for reference in ['%s%d'%(ELEM[Z],A)]:
#   for reference in {10:['B10','He4','He8'],22:['Na22','O16'],46:['V46','Ca40']}[A]:
#   for reference in ['He%d'%A,'Li%d'%A]:
#   for reference in ['O22']:
#  for method in ['flow']:
#   for hw in [14,22]:
#   for hw in [22]:
     ARGS['reference'] = reference
#     ARGS['valence_space'] = {10:'p-shell',22:'sd-shell',46:'fp-shell'}[A]
#     ARGS['valence_space'] = 'psdNR-shell'
#   ARGS['reference'] = 'O%d'%refmass
#     ARGS['2bme'] = 'input/chi2bSAT_srg0000_eMax14_lMax10_hwHO0%d.me2j.gz'%(hw)
#     ARGS['3bme'] = 'input/chi2b3bSAT_J7666_hwconv022_JT3Nfull73_srg0000ho40C_eMax14_EMax14_hwHO0%d.me3j.gz'%(hw)
#     ARGS['LECs'] = 'SAT'
     ARGS['2bme'] = 'input/chi2b_srg0625_eMax14_lMax10_hwHO0%d.me2j.gz'%(hw)
     ARGS['3bme'] = 'input/me3j/chi2b3b400cD-02cE0098_hwconv036_srg0625ho40J_eMax14_EMax14_hwHO0%d.me3j.gz'%(hw)
     ARGS['LECs'] = 'srg0625'
#     ARGS['3bme'] = 'none'
#     ARGS['LECs'] = 'NNonlysrg0625'
#     ARGS['2bme'] = 'input/chi2b_srg0800_eMax14_lMax10_hwHO024.me2j.gz'
#     ARGS['2bme'] = 'input/chi2b_srg0800_eMax14_lMax10_hwHO%03d.me2j.gz'%(hw)
#     ARGS['3bme'] = 'input/chi2b3b400cD-02cE0098_srg0800ho40J_hwconv36_eMax14_EMax14_hwHO024.me3j.gz'
#     ARGS['3bme'] = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO%03d.me3j.gz'%(hw)
#     ARGS['3bme'] = 'input/me3j/chi2b3b400cD-02cE0098_srg0800ho40J_hwconv36_eMax14_EMax14_hwHO%03d.me3j.gz'%(hw)
#     ARGS['LECs'] = 'srg0800'
#     ARGS['2bme'] = '/work/hda21/hda219/input/chi2b_srg0625_eMax14_lMax10_hwHO%03d.me2j.gz'%(hw)
#     ARGS['2bme'] = '/work/hda21/hda215/ME_share/vnn_hw%d.00_kvnn10_lambda1.80_mesh_kmax_7.0_100_pc_R15.00_N15.dat_to_me2j.gz'%(hw)
#     ARGS['3bme'] = '/work/hda21/hda215/ME_share/jsTNF_Nmax_16_J12max_8_hbarOmega_%d.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_new_E3_16_e_14_ant_EM1.8_2.0.h5_to_me3j.gz'%(hw)
#     ARGS['LECs'] = 'EM1.8_2.0'
#     ARGS['3bme'] = 'none'
#     ARGS['LECs'] = 'NNonlysrg0625'
#     ARGS['3bme'] = 'input/me3j/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO%03d.me3j.gz'%(hw)
#     ARGS['LECs'] = 'srg0800ho40C'
#     ARGS['3bme'] = 'input/chi2b3b_hwconv036_srg0800ho40J_eMax14_EMax14_hwHO024.me3j.gz'
#     ARGS['LECs'] = '500srg0800'
#     ARGS['3bme'] = 'input/me3j/chi2b3b_hwconv036_srg0800ho40J_eMax14_EMax14_hwHO%03d.me3j.gz'%(hw)
#     ARGS['LECs'] = '500srg0800ho40J'
#     ARGS['3bme'] = 'input/me3j/chi2b3b_hwconv036_srg0800ho40J_eMax14_EMax14_hwHO%03d.me3j.gz'%(hw)
#     ARGS['LECs'] = '500srg0800ho40C'
#     ARGS['3bme'] = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO024.me3j.gz'
#     ARGS['2bme'] = 'input/chi2b_srg0800_eMax14_lMax10_hwHO%03d.me2j.gz'%(hw)
#     ARGS['3bme'] = 'input/chi2b3b400cD-02cE0098_srg0800ho40J_hwconv36_eMax14_EMax14_hwHO%03d.me3j.gz'%(hw)
#     ARGS['3bme'] = 'input/chi2b3b_hwconv036_srg0800ho40J_eMax14_EMax14_hwHO024.me3j.gz'
#     ARGS['LECs'] = '500srg0800'
#     ARGS['3bme'] = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax14_EMax14_hwHO%03d.me3j.gz'%(hw)
#     ARGS['LECs'] = 'NFCsrg0800'
#     ARGS['2bme'] = 'input/chi2b_srg0800_eMax14_lMax10_hwHO024.me2j.gz'
#     ARGS['2bme'] = 'input/chi2b_srg0625_eMax14_lMax10_hwHO0%d.me2j.gz'%(hw)
#     ARGS['3bme'] = 'input/jsTNF_Nmax_18_J12max_8_hbarOmega_28.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_id_1_new_ant_E3_14_e_14_new.h5'
#     ARGS['LECs'] = 'PWA2.0_2.0'
#     ARGS['LECs'] = 'EM2.0_2.0'
     ARGS['hw'] = '%d'%hw
     ARGS['A'] = '%d'%A
#     ARGS['valence_space'] = reference
     ARGS['valence_space'] = '0hw-shell'
#     ARGS['valence_space'] = 'p-shell'
#     ARGS['valence_space'] = 'sd-shell'
#     ARGS['valence_space'] = 'sdfpNR-shell'
#     ARGS['valence_space'] = 'psdNR-shell'
#     ARGS['valence_space'] = 'fp-shell'
#     ARGS['valence_space'] = 'fpg9NR-shell'
#     ARGS['valence_space'] = 'fpgdsNR-shell'
#     ARGS['valence_space'] = 'gds-shell'
#     ARGS['valence_space'] = 'Cr%d'%A
#     if A<30:
#       ARGS['valence_space'] = 'sd-shell'
#     else:
#       ARGS['valence_space'] = 'sdfpNR-shell'
#     ARGS['core_generator'] = 'imaginary-time'
#     ARGS['valence_generator'] = 'shell-model-imaginary-time'
     ARGS['emax'] = '%d'%e
#     ARGS['method'] = method

     ARGS['Operators'] = ''
#     ARGS['Operators'] = 'Rp2'
#     ARGS['Operators'] = 'Rp2,Rn2'
#     ARGS['Operators'] = 'Rp2Z8,Rp2Z10'
#     ARGS['Operators'] = 'Rp2Z16'
#     ARGS['Operators'] = 'E2'
#     ARGS['Operators'] = 'E2,M1'
#     ARGS['Operators'] = 'E2,M1,GamowTeller'
     ARGS['Operators'] = 'M1p,M1n,Sigma_p,Sigma_n'
#     ARGS['Operators'] = 'Sigma_n'
#     ARGS['Operators'] = 'E2,M1,Rp2Z10'
#     ARGS['Operators'] = 'Rp2Z10'
#     ARGS['Operators'] = 'rhop0.0,R2CM'
#     ARGS['Operators'] = 'GamowTeller'

     time_request = '24:00:00'
     if e<5: time_request = '00:10:00'
     if e < 8 : time_request = '01:00:00'
     elif e < 10 : time_request = '04:00:00'
     elif e < 12 : time_request = '12:00:00'
     jobname  = '%s_%s_%s_%s_e%s_E%s_s%s_hw%s_A%s' %(ARGS['valence_space'], ARGS['LECs'],ARGS['method'],ARGS['reference'],ARGS['emax'],ARGS['e3max'],ARGS['smax'],ARGS['hw'],ARGS['A'])
     if 'lmax3' in ARGS:  jobname  += '_l%d'%(ARGS['lmax3'])
     if 'eta_criterion' in ARGS: jobname += '_eta%s'%(ARGS['eta_criterion'])
     if 'core_generator' in ARGS: jobname += '_' + ARGS['core_generator']
     if 'BetaCM' in ARGS: jobname += '_' + ARGS['BetaCM']
     ARGS['flowfile'] = 'output/BCH_' + jobname + '.dat'
     ARGS['intfile']  = 'output/' + jobname
     cmd = '%s %s'%(exe,' '.join(['%s=%s'%(x,ARGS[x]) for x in ARGS]) )
     if batch_mode==True:
       sfile = open(jobname,'w')
       sfile.write(FILECONTENT%(NTHREADS,jobname,time_request,NTHREADS,NTHREADS,cmd))
       sfile.close()
       call(['sbatch', jobname])
       remove(jobname) # delete the file
       sleep(0.1)
     else:
       call(cmd.split())  # Run in the terminal, rather than submitting


