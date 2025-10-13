#!/usr/bin/env python

##########################################################################
# goUniversal.py
##
# A python script to run or submit jobs for the common use cases
# of the IMSRG++ code. We check whether there is a pbs or slurm
# scheduler, assign the relevant input parameters, set names
# for the output files, and run or submit.
# -Ragnar Stroberg
# TRIUMF Nov 2016
######################################################################

from os import path, environ, mkdir, remove, system
from sys import argv
from subprocess import call, PIPE
from time import time, sleep
from datetime import datetime

# Check to see what type of batch submission system we're dealing with
BATCHSYS = 'NONE'
if call('type '+'qsub', shell=True, stdout=PIPE, stderr=PIPE) == 0:
  BATCHSYS = 'PBS'
elif call('type '+'srun', shell=True, stdout=PIPE, stderr=PIPE) == 0:
  BATCHSYS = 'SLURM'
  system("source ~/.bashrc")

# The code uses OpenMP and benefits from up to at least 24 threads
NTHREADS = 48
exe = '/Users/antoinebelley/Documents/TRIUMF/imsrg/src/imsrg++'

# Flag to swith between submitting to the scheduler or running in the current shell
# batch_mode=False
batch_mode = False
if 'terminal' in argv[1:]:
  batch_mode = False

# Don't forget to change this. I don't want emails about your calculations...
mail_address = 'antoine.belley@mail.mcgill.ca'

# This comes in handy if you want to loop over Z
ELEM = ['n', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N',
        'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
        'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',  'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
        'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',  'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
        'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']  # ,'Bi','Po','At','Rn','Fr','Ra','Ac','Th','U','Np','Pu']

# ARGS is a (string => string) dictionary of input variables that are passed to the main program
ARGS = {}

# Maximum value of s, and maximum step size ds
ARGS['smax'] = '0'
ARGS['dsmax'] = '0.5'
ARGS['basis'] = "oscillator"
# ARGS['lmax3'] = '10' # for comparing with Heiko

# Norm of Omega at which we split off and start a new transformation
ARGS['omega_norm_max'] = '0.25'

# Model space parameters used for reading Darmstadt-style interaction files
ARGS['file2e1max'] = '18 file2e2max=36 file2lmax=18'
ARGS['file3e1max'] = '18 file3e2max=36 file3e3max=24'

# Name of a directory to write Omega operators so they don't need to be stored in memory. If not given, they'll just be stored in memory.
# ARGS['scratch'] = 'SCRATCH'

# Generator for core decoupling, can be atan, white, imaginary-time.  (atan is default)
# ARGS['core_generator'] = 'imaginary-time'
# Generator for valence deoupling, can be shell-model, shell-model-atan, shell-model-npnh, shell-model-imaginary-time (shell-model-atan is default)
# ARGS['valence_generator'] = 'shell-model-imaginary-time'

# Solution method
ARGS['method'] = 'magnus'
# ARGS['method'] = 'brueckner'
# ARGS['method'] = 'flow'
# ARGS['method'] = 'HF'
# ARGS['method'] = 'MP3'

# Tolerance for ODE solver if using flow solution method
# ARGS['ode_tolerance'] = '1e-5'

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
#SBATCH --account=rrg-holt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#SBATCH --output=/home/alextodd/projects/rrg-holt/alextodd/results/imsrg_log/%s.%%j
#SBATCH --time=%s
#SBATCH --mail-user=%s
#SBATCH --mail-type=END
#SBATCH --mem=187G
cd $SLURM_SUBMIT_DIR
echo NTHREADS = %d
export OMP_NUM_THREADS=%d
time srun %s
"""

# Make a directory for the log files, if it doesn't already exist
# if not path.exists('imsrg_log'): mkdir('imsrg_log')

# Loop over multiple jobs to submit
for Z in range(4, 5):
 A = 10
 for reference in ['%s%d' % (ELEM[Z], A)]:
  ARGS['reference'] = reference
  print('Z = ', Z)
  for e in [4]:
   for hw in [16]:

     ARGS['emax'] = '%d' % e
     ARGS['e3max'] = '12'

     # Delta
    # [2007, 2460, 3448, 3373, 1141, 3319, 3098, 1429, 3895, 1172, 90, 723, 2125, 2245, 750, 1469, 1177, 493, 500, 3621, 606, 3480, 3260, 3813, 3105, 2411, 3350, 774, 2618, 4117, 922, 1173, 1802, 3472]
     # ARGS['2bme'] = '/home/alextodd/projects/def-holt/shared/me2j/TwBME-HO_NN-only_DeltaGO394_sample90_bare_hw10_emax16_e2max32.me2j.gz'
     # ARGS['3bme'] = '/home/alextodd/projects/def-holt/shared/me3j/delta_samples/NO2B_half_ThBME_3NFJmax15_DN2LOGO394_sample90_IS_hw10_ms16_32_28.stream.bin'
     # ARGS['LECs'] = 'DeltaGO'

     # N3LO
     # ARGS['2bme'] = '/home/alextodd/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N3LO_EMN500_srg2.0_hw%d_emax16_e2max32.me2j.gz'%(hw)
     # ARGS['3bme'] = '/home/alextodd/projects/def-holt/shared/me3j/NO2B_ThBME_srg2.0_ramp46-9-44-15-42_N3LO_EM500_3NFJmax15_c1_-0.81_c3_-3.2_c4_5.4_cD_0.7_cE_-0.06_LNL2_650_500_IS_hw%dfrom30_ms16_32_22.stream.bin'%(hw)
     # ARGS['LECs'] = 'N3LO_EM500'

     # EM(1.8/2.0)
     ARGS['2bme'] = '/Users/antoinebelley/Documents/TRIUMF/Interactions/TwBME-HO_NN-only_N3LO_EM500_srg1.80_hw16_emax18_e2max36.me2j.gz'
     ARGS['3bme'] = '/Users/antoinebelley/Documents/TRIUMF/Interactions/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_24.stream.bin'
     ARGS['LECs'] = 'EM1.8_2.0'

     ARGS["3bme_type"] = "no2b"
     # ARGS["no2b_precision"] = "half"

     ARGS['hw'] = '%d' % hw
     ARGS['A'] = '%d' % A
     # ARGS['valence_space'] = 'fp-shell'
     # this is just a label when custom_valence_space is set
     ARGS['valence_space'] = 'p-shell'
     # ARGS['custom_valence_space'] = 'Ni56,p0f5,n0f5,p1p3,n1p3,p1p1,n1p1,p0g9,n0g9'  # AKA: Ni56 core with jj44pn
     # ARGS['valence_space'] = '3hw-shell'
#     ARGS['valence_space'] = 'Cr%d'%A
#     ARGS['core_generator'] = 'imaginary-time'
#     ARGS['valence_generator'] = 'shell-model-imaginary-time'
     # ARGS['emax'] = '%d'%e
#     ARGS['method'] = method

     # ARGS['Operators'] = 'GamowTeller'    # Operators to consistenly transform, separated by commas.
     ARGS['Operators'] = 'M0nu_F_3.54_none,M0nu_GT_3.54_none,M0nu_T_3.54_none'
#     ARGS['Operators'] = 'Rp2,Rn2'
#     ARGS['Operators'] = 'E2'
#     ARGS['Operators'] = 'E2,M1'
#     ARGS['Operators'] = 'E2,M1,GamowTeller'
#     ARGS['Operators'] = 'M1p,M1n,Sigma_p,Sigma_n'
#     ARGS['Operators'] = 'GamowTeller'
    #  ARGS['input_op_fmt'] = 'miyagi'

     # all 6 short range operators:
     # ARGS['OperatorsFromFile'] ='0vbbFermiVVSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbFermiVV_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz,0vbbGamowTellerAASR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbGamowTellerAA_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz,0vbbGamowTellerAPSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbGamowTellerAP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz,0vbbGamowTellerPPSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbGamowTellerPP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz,0vbbTensorAPSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbTensorAP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz,0bvvTensorPPSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbTensorPP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz'

     # two tensor operators:
     # ARGS['OperatorsFromFile'] ='0vbbTensorAPSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbTensorAP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz,0bvvTensorPPSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbTensorPP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz'

     # ARGS['OperatorsFromFile'] ='0vbbGamowTellerPPSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbGamowTellerPP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz'
     # ARGS['OperatorsFromFile'] ='0vbbFermiVVSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbFermiVV_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz'

     # ARGS['OperatorsFromFile']= '0vbbTensorPPSR^0_2_0_2^/home/alextodd/projects/rrg-holt/alextodd/NuHamil-public/exe/0vbbTensorPP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz'

    # Make an estimate of how much time to request. Only used for slurm at the moment.
     # if   e <  5 : time_request = '00:10:00'
     # elif e <  8 : time_request = '01:00:00'
     # elif e < 10 : time_request = '04:00:00'
     # elif e < 12 : time_request = '12:00:00'
     # time_request = '100:00:00'
     if e <= 4:
       time_request = '00:06:00'
     elif e <= 6:
       time_request = '01:00:00'
     elif e <= 8:
       time_request = '6:00:00'
     elif e <= 10:
       time_request = '12:00:00'
     elif e <= 12:
       time_request = '24:00:00'
     elif e <= 14:
       time_request = '35:00:00'
     # time_request = '40:00:00'
# /home/alextodd/projects/rrg-holt/alextodd/results
     jobname = '%s_%s_%s_%s_e%s_E%s_s%s_hw%s_A%s' % (
         ARGS['valence_space'], ARGS['LECs'], ARGS['method'], ARGS['reference'], ARGS['emax'], ARGS['e3max'], ARGS['smax'], ARGS['hw'], ARGS['A'])
     logname = jobname + \
         datetime.fromtimestamp(time()).strftime('_%y%m%d%H%M.log')
     print(jobname)
  # Some optional parameters that we probably want in the output name if we're using them
     if 'lmax3' in ARGS:
       jobname += '_l%d' % (ARGS['lmax3'])
     if 'eta_criterion' in ARGS:
       jobname += '_eta%s' % (ARGS['eta_criterion'])
     if 'core_generator' in ARGS:
       jobname += '_' + ARGS['core_generator']
     if 'BetaCM' in ARGS:
       jobname += '_' + ARGS['BetaCM']
     ARGS['flowfile'] = '/Users/antoinebelley/Documents/TRIUMF/results/BCH_' + jobname + '.dat'
     ARGS['intfile'] = '/Users/antoinebelley/Documents/TRIUMF/results/' + jobname

     cmd = ' '.join([exe] + ['%s=%s' % (x, ARGS[x]) for x in ARGS])
     print(cmd)

  # Submit the job if we're running in batch mode, otherwise just run in the current shell
     if batch_mode == True:
       sfile = open(jobname+'.batch', 'w')
       if BATCHSYS == 'PBS':
         sfile.write(FILECONTENT % (
             jobname, environ['PWD'], NTHREADS, mail_address, logname, NTHREADS, cmd))
         sfile.close()
         call(['qsub', jobname+'.batch'])
       elif BATCHSYS == 'SLURM':
         sfile.write(FILECONTENT % (NTHREADS, jobname, time_request,
                     mail_address, NTHREADS, NTHREADS, cmd))
         sfile.close()
         call(['sbatch', jobname+'.batch'])
       remove(jobname+'.batch')  # delete the file
       sleep(0.1)
     else:
       call(cmd.split())  # Run in the terminal, rather than submitting
