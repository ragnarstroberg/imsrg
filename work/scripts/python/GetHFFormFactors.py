#!/usr/bin/env python

from pyIMSRG import *
from numpy import array,arange,linspace
from os import path

###################################################################################
################################# Inputs go here ##################################
emax = 4
E3max = 14
hw = 20
file2N = '../input/chi2b_srg0800_eMax12_lMax10_hwHO020.me2j.gz'
file3N = '../input/chi2b3b_srg0800ho40C_eMax12_EMax12_hwHO020.me3j'
#file3N = 'none'


reference = 'O16'  # The reference for the normal ordering
core = 'O16'       # The core of the valence space. At the end, normal order wrt this
basis = 'HF'       # basis can be either 'oscillator' or 'HF' (for Hartree-Fock) or 'IMSRG'

output_file_name = 'SP_wavefunction_%s_%s_e%d_hw%d'%(reference,basis,emax,hw)

formatNN = 'darmstadt'
format3N = 'darmstadt'  # can be 'darmstadt' or 'navratil'

file2N_emax  = 14  # These parameter are the truncations of the interaction file
file2N_e2max = 28  # if the file is in Darmstad format. If these are not set
file2N_lmax  = 10  # properly, then I'll misread the file.

file3N_emax  = 12  #  Same deal for the 3N file
file3N_e2max = 28  # 
file3N_e3max = 12  # 

orbit_labels = ['p0d5','p0d3','p1s1','n0d5','n0d3','n1s1']

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
  modelspace = ModelSpace(emax,reference,reference)
  modelspace.SetHbarOmega(hw)
  
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
      
  
  ## other inputs not yet supported in this script, because I'm lazy.
  
  Hbare += Trel_Op(modelspace)
  
  
  ## If we want an oscillator reference, just do the normal ordering.
  if basis == 'oscillator':
    HNO = Hbare.DoNormalOrdering()
  
  ## If we want a Hartree Fock reference, first do the HF calculation, then normal order.
  dr = 0.05
  R = arange(0,10,dr)
  if basis == 'HF' or basis == 'IMSRG':
    hf = HartreeFock(Hbare)
    hf.Solve()
    for orb in orbit_labels:
      index = modelspace.GetOrbitIndex_fromString(orb)
      print orb,index
      f = open(output_file_name + '_'+ orb + '.dat','w')
      Psi = []
      for r in R:
       Psi.append( hf.GetRadialWF_r(index,r) )
      for (r,psi) in zip(R,Psi):
        f.write('%12.6f  %12.6e\n'%( r,psi ))
      f.close()
      norm = 0
      for (r,psi) in zip(R,Psi): norm += r*r*dr*psi*psi
      print 'Norm = ',norm
#    HNO = hf.GetNormalOrderedH()

  ## There are lots of solver parameters you can tweak, but they really shouldn't
  ## have much of an effect, so let's just go with some defaults.
  if basis == 'IMSRG':
    imsrgsolver = IMSRGSolver(HNO)
    imsrgsolver.SetSmax(500)
    imsrgsolver.SetReadWrite(rw)
    imsrgsolver.SetEtaCriterion(1e-6)
    imsrgsolver.SetMethod('magnus')
    imsrgsolver.SetOmegaNormMax(0.5)
    imsrgsolver.SetGenerator('atan')
    imsrgsolver.Solve()
    HNO = imsrgsolver.GetH_s()
    
  
#  ## If the core of the eventual valence space is different from normal ordering reference
#  ## then we should at this point switch to the core.
#  if core != reference:
#    ms_core = modelspace(emax,core,core)
#    HNO = HNO.UndoNormalOrdering()
#    HNO.SetModelSpace(ms_core)
#    HNO = HNO.DoNormalOrdering()
  
  
  ## Now we write out the interaction in Oslo format for use with MBPT code
  
#  rw.WriteOneBody_Oslo( output_file_name + '_1N.dat', HNO)
#  rw.WriteTwoBody_Oslo( output_file_name + '_2N.dat', HNO)
  
  Hbare.PrintTimes()
  

##########################
##########################

if __name__ == '__main__':
  main()


