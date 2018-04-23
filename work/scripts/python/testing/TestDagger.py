#!/usr/bin/env python

###################################################################
#
#
###################################################################

from pyIMSRG import *
from os import path



################################# Inputs go here ##################################
emax = 4
E3max = 14
hw = 16
#file2N = '../../input/chi2b_srg0625_eMax14_lMax10_hwHO020.me2j.gz'
#file3N = '../input/chi2b3b_srg0800ho40C_eMax12_EMax12_hwHO020.me3j'
#file3N = 'none'
file2N = '/global/scratch/exch/me2j/fromJohannes/vnn_hw16.00_kvnn10_lambda1.80_mesh_kmax_7.0_100_pc_R15.00_N15.dat_to_me2j.gz'
file3N = '/global/scratch/exch/me3j/fromJohannes/jsTNF_Nmax_16_J12max_8_hbarOmega_16.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_new_E3_14_e_14_ant_EM1.8_2.0.h5_to_me3j.gz'


reference = 'O16'  # The reference for the normal ordering
core = 'O16'       # The core of the valence space. At the end, normal order wrt this
valence = 'sd-shell'
basis = 'IMSRG'       # basis can be either 'oscillator' or 'HF' (for Hartree-Fock) or 'IMSRG'

output_file_name = 'NormalOrdered_%s_%s_e%d_hw%d'%(reference,basis,emax,hw)

formatNN = 'darmstadt'
format3N = 'darmstadt'  # can be 'darmstadt' or 'navratil'

file2N_emax  = 14  # These parameter are the truncations of the interaction file
file2N_e2max = 28  # if the file is in Darmstad format. If these are not set
file2N_lmax  = 14  # properly, then I'll misread the file.

file3N_emax  = 14  #  Same deal for the 3N file
file3N_e2max = 28  # 
file3N_e3max = 14  #

daggerorbit = "p0d3"
###################################################################################
###################################################################################


def main():
  # Before doing too much, make sure the files are there...
  if not path.exists( file2N ):
    print "Uh, oh. I can't open file:",file2N
    return
  
  if (file3N not in ['','none']) and (not path.exists(file3N)):
    print "Uh, oh. I can't open file:",file3N
    return
  
  

  rw = ReadWrite()
  modelspace = ModelSpace(emax,reference,valence)
  modelspace.SetHbarOmega(hw)
  modelspace.SetE3max(E3max)

  Q = modelspace.GetOrbitIndex_fromString(daggerorbit)

  particle_rank = 3
  if file3N == '' or file3N == 'none': particle_rank = 2
  
  Hbare = Operator(modelspace,0,0,0,particle_rank)
  
  
  if formatNN == 'darmstadt':
    rw.ReadBareTBME_Darmstadt( file2N, Hbare, file2N_emax, file2N_e2max, file2N_lmax )
  
  
  # Due to unfortunate naming conventions, Read_Darmstadt_3body also handles
  #  Petr Navratil's 3N format. Sorry about that.
  if particle_rank > 2:
    if format3N == 'darmstadt' or format3N=='me3j':
      rw.Set3NFormat( 'me3j' )
    else:
      rw.Set3NFormat( format3N )
    rw.Read_Darmstadt_3body( file3N, Hbare, file3N_emax, file3N_e2max, file3N_e3max )
      
  
  dag = DaggerOperator(modelspace,Q)
  Hbare += Trel_Op(modelspace)

  ## If we want an oscillator reference, just do the normal ordering.
  if basis == 'oscillator':
    HNO = Hbare.DoNormalOrdering()
  
  ## If we want a Hartree Fock reference, first do the HF calculation, then normal order.
  elif basis == 'HF' or basis == 'IMSRG':
    hf = HartreeFock(Hbare)
    hf.Solve()
    HNO = hf.GetNormalOrderedH()

  ## There are lots of solver parameters you can tweak, but they really shouldn't
  ## have much of an effect, so let's just go with some defaults.
  if basis == 'IMSRG':
    imsrgsolver = IMSRGSolver(HNO)
#    imsrgsolver.SetSmax(500)
    imsrgsolver.SetSmax(5)
    imsrgsolver.SetReadWrite(rw)
    imsrgsolver.SetEtaCriterion(1e-6)
    imsrgsolver.SetMethod('magnus')
    imsrgsolver.SetOmegaNormMax(0.5)
    imsrgsolver.SetGenerator('atan')
    imsrgsolver.Solve()
    HNO = imsrgsolver.GetH_s()
    dagT = imsrgsolver.TransformDagger(dag)
    rw.WriteNuShellX_int(dagT,'dagger_%s'%(daggerorbit) + '_hw%d_e%d_%s_%s.int'%(hw,emax,reference,valence))
#  dagT.PrintOneBody()

 


  ## If the core of the eventual valence space is different from normal ordering reference
  ## then we should at this point switch to the core.
  if core != reference:
    ms_core = modelspace(emax,core,core)
    HNO = HNO.UndoNormalOrdering()
    HNO.SetModelSpace(ms_core)
    HNO = HNO.DoNormalOrdering()





  
  print dag.GetQSpaceOrbit()


  HNO.PrintTimes()

##########################
##########################

if __name__ == '__main__':
  main()
