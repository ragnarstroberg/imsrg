#!/usr/bin/env python

##########################################################################
##  goUniversal.py
##
##  A python script to run or submit jobs for the common use cases
##  of the IMSRG++ code. We check whether there is a pbs or slurm
##  scheduler, assign the relevant input parameters, set names
##  for the output files, and run or submit.
##  						-Ragnar Stroberg
##  						TRIUMF Nov 2016
######################################################################

from os import path,environ,mkdir,remove
from subprocess import call,PIPE
from time import time,sleep
from datetime import datetime

### Check to see what type of batch submission system we're dealing with
BATCHSYS = 'NONE'
if call('type '+'qsub', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'PBS'
elif call('type '+'srun', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'SLURM'

### The code uses OpenMP and benefits from up to at least 24 threads
NTHREADS=12
exe = '%s/bin/imsrg++'%(environ['HOME'])

### Flag to swith between submitting to the scheduler or running in the current shell
batch_mode=False
#batch_mode=True

### Don't forget to change this. I don't want emails about your calculations...
mail_address = 'sstroberg@triumf.ca'

### This comes in handy if you want to loop over Z
ELEM = ['n','H','He','Li','Be','B','C','N',
       'O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',
       'Ca','Sc','Ti','V','Cr','Mn','Fe','Co',  'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']# ,'Bi','Po','At','Rn','Fr','Ra','Ac','Th','U','Np','Pu']

### ARGS is a (string => string) dictionary of input variables that are passed to the main program
ARGS  =  {}

### Maximum value of s, and maximum step size ds
ARGS['smax'] = '500'
ARGS['dsmax'] = '0.5'

#ARGS['lmax3'] = '10' # for comparing with Heiko

### Norm of Omega at which we split off and start a new transformation
ARGS['omega_norm_max'] = '0.25'

### Model space parameters used for reading Darmstadt-style interaction files
ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=10'
ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=14'

### Name of a directory to write Omega operators so they don't need to be stored in memory. If not given, they'll just be stored in memory.
#ARGS['scratch'] = 'SCRATCH'    

### Generator for core decoupling, can be atan, white, imaginary-time.  (atan is default)
#ARGS['core_generator'] = 'imaginary-time' 
### Generator for valence deoupling, can be shell-model, shell-model-atan, shell-model-npnh, shell-model-imaginary-time (shell-model-atan is default)
#ARGS['valence_generator'] = 'shell-model-imaginary-time' 

### Solution method
ARGS['method'] = 'magnus'
#ARGS['method'] = 'brueckner'
#ARGS['method'] = 'flow'
#ARGS['method'] = 'HF'
#ARGS['method'] = 'MP3'

### Tolerance for ODE solver if using flow solution method
#ARGS['ode_tolerance'] = '1e-5'

if BATCHSYS == 'PBS':
  FILECONTENT = """#!/bin/bash
#PBS -N %s
#PBS -q batchmpi
#PBS -d %s
#PBS -l walltime=192:00:00
#PBS -l nodes=1:ppn=%d
#PBS -l vmem=60gb
#PBS -m ae
#PBS -M %s
#PBS -j oe
#PBS -o imsrg_log/%s.o
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=%d
%s
  """

elif BATCHSYS == 'SLURM':
  FILECONTENT = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#SBATCH --output=imsrg_log/%s.%%j
#SBATCH --time=%s
#SBATCH --mail-user=%s
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
echo NTHREADS = %d
export OMP_NUM_THREADS=%d
time srun %s
"""

### Make a directory for the log files, if it doesn't already exist
if not path.exists('imsrg_log'): mkdir('imsrg_log')

### Loop over multiple jobs to submit
for Z in range(16,17):
 A=32
 for reference in ['%s%d'%(ELEM[Z],A)]:
  ARGS['reference'] = reference
  print 'Z = ',Z
  for e in [2]:
   for hw in [13]:

     ARGS['emax'] = '%d'%e
     ARGS['e3max'] = '14'

     ARGS['2bme'] = 'input/chi2b_srg0625_eMax14_lMax10_hwHO0%d.me2j.gz'%(hw)
     ARGS['3bme'] = 'input/me3j/chi2b3b400cD-02cE0098_hwconv036_srg0625ho40J_eMax14_EMax14_hwHO0%d.me3j.gz'%(hw)
     ARGS['LECs'] = 'srg0625'

#     ARGS['2bme'] = 'input/usdbpn.int'
#     ARGS['3bme'] = 'none'
#     ARGS['LECs'] = 'usdb'
#     ARGS['fmt2'] = 'nushellx'
#     ARGS['basis'] = 'oscillator'

#     ARGS['2bme'] = '/work/hda21/hda215/ME_share/vnn_hw%d.00_kvnn10_lambda1.80_mesh_kmax_7.0_100_pc_R15.00_N15.dat_to_me2j.gz'%(hw)
#     ARGS['3bme'] = '/work/hda21/hda215/ME_share/jsTNF_Nmax_16_J12max_8_hbarOmega_%d.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_new_E3_16_e_14_ant_EM1.8_2.0.h5_to_me3j.gz'%(hw)
#     ARGS['2bme'] = '/itch/exch/me2j/chi2b_srg0625_eMax14_lMax10_hwHO0%d.me2j.gz'%(hw)
#     ARGS['3bme'] = '/itch/exch/me3j/new/chi2b3b400cD-02cE0098_hwconv036_srg0625ho40J_eMax14_EMax14_hwHO0%d.me3j.gz'%(hw)
#     ARGS['LECs'] = 'srg0625'

     ARGS['hw'] = '%d'%hw
     ARGS['A'] = '%d'%A
     ARGS['valence_space'] = reference
#     ARGS['valence_space'] = '0hw-shell'
#     ARGS['valence_space'] = 'Cr%d'%A
#     ARGS['core_generator'] = 'imaginary-time'
#     ARGS['valence_generator'] = 'shell-model-imaginary-time'
     ARGS['emax'] = '%d'%e
#     ARGS['method'] = method

     ARGS['Operators'] = ''    # Operators to consistenly transform, separated by commas.
#     ARGS['Operators'] = 'Rp2'
#     ARGS['Operators'] = 'Rp2,Rn2'
#     ARGS['Operators'] = 'E2'
#     ARGS['Operators'] = 'E2,M1'
#     ARGS['Operators'] = 'E2,M1,GamowTeller'
#     ARGS['Operators'] = 'M1p,M1n,Sigma_p,Sigma_n'
#     ARGS['Operators'] = 'GamowTeller'



    ### Make an estimate of how much time to request. Only used for slurm at the moment.
     time_request = '24:00:00'
     if   e <  5 : time_request = '00:10:00'
     elif e <  8 : time_request = '01:00:00'
     elif e < 10 : time_request = '04:00:00'
     elif e < 12 : time_request = '12:00:00'

     jobname  = '%s_%s_%s_%s_e%s_E%s_s%s_hw%s_A%s' %(ARGS['valence_space'], ARGS['LECs'],ARGS['method'],ARGS['reference'],ARGS['emax'],ARGS['e3max'],ARGS['smax'],ARGS['hw'],ARGS['A'])
     logname = jobname + datetime.fromtimestamp(time()).strftime('_%y%m%d%H%M.log')

  ### Some optional parameters that we probably want in the output name if we're using them
     if 'lmax3' in ARGS:  jobname  += '_l%d'%(ARGS['lmax3'])
     if 'eta_criterion' in ARGS: jobname += '_eta%s'%(ARGS['eta_criterion'])
     if 'core_generator' in ARGS: jobname += '_' + ARGS['core_generator']
     if 'BetaCM' in ARGS: jobname += '_' + ARGS['BetaCM']
     ARGS['flowfile'] = 'output/BCH_' + jobname + '.dat'
     ARGS['intfile']  = 'output/' + jobname

     cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])

  ### Submit the job if we're running in batch mode, otherwise just run in the current shell
     if batch_mode==True:
       sfile = open(jobname+'.batch','w')
       if BATCHSYS == 'PBS':
         sfile.write(FILECONTENT%(jobname,environ['PWD'],NTHREADS,mail_address,logname,NTHREADS,cmd))
         sfile.close()
         call(['qsub', jobname+'.batch'])
       elif BATCHSYS == 'SLURM':
         sfile.write(FILECONTENT%(NTHREADS,jobname,time_request,mail_address,NTHREADS,NTHREADS,cmd))
         sfile.close()
         call(['sbatch', jobname])
       remove(jobname+'.batch') # delete the file
       sleep(0.1)
     else:
       call(cmd.split())  # Run in the terminal, rather than submitting

