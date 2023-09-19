
from pyIMSRG import *

emax = 2
nuc = 'C12'
E3max = 3 * emax
particle_rank = 3
use_imsrg3 = True
use_imsrg3n7 = True

print('ref= {}\temax= {}\tE3max= {}\t IMSRG3= {}\t IMSRG3n7= {}'.format(nuc,emax,E3max, use_imsrg3,use_imsrg3n7))

ms = ModelSpace( emax, nuc, nuc )
ms.SetE3max(E3max)


print(len(ms.holes))



if particle_rank>2:
    Commutator.SetUseIMSRG3(use_imsrg3)
    #Commutator.SetUseIMSRG3N7(use_imsrg3n7)

utest = UnitTest(ms)
H   = utest.RandomOp(ms,0,0,0,particle_rank, +1)
Eta = utest.RandomOp(ms,0,0,0,particle_rank, -1)

H.ThreeBody.TransformToPN()
Eta.ThreeBody.TransformToPN()

H.ScaleOneBody(0.)
Eta.ScaleOneBody(0.)

H.ThreeBody.Erase()
Eta.ThreeBody.Erase()

Htemp = Operator(ms,0,0,0,particle_rank)
Hout1 = Operator(ms,0,0,0,particle_rank)
Hout2 = Operator(ms,0,0,0,particle_rank)
Hout3 = Operator(ms,0,0,0,particle_rank)

Htemp.ThreeBody.SetMode('pn')     
Hout1.ThreeBody.SetMode('pn')    
Hout2.ThreeBody.SetMode('pn')    
Hout3.ThreeBody.SetMode('pn')    

Htemp.ThreeBody.Erase()
Hout1.ThreeBody.Erase()
Hout2.ThreeBody.Erase()
Hout3.ThreeBody.Erase()

Htemp.TwoBody.Erase()

Htemp.ScaleOneBody(0.)
Hout1.ScaleOneBody(0.)
Hout2.ScaleOneBody(0.)
Hout3.ScaleOneBody(0.)

# double commutator
Commutator.comm223ss(Eta, H, Htemp)
# Htemp.TwoBody.Erase()
Commutator.comm231ss(Eta, Htemp, Hout1)

# direct term
# ReferenceImplementations.comm223_231_BruteForce(Eta, H, Hout2)

# Factorization
Commutator.comm223_231_Factorization(Eta, H, Hout3)

print("\n   One Body part: ")
print( "H1    {}    ".format(  Hout1.OneBodyNorm()) )
print( "H2    {}    ".format(  Hout2.OneBodyNorm() ))
print( "H3    {}    ".format(  Hout3.OneBodyNorm() ))
print("\n")





Hout2_1 = Operator(ms,0,0,0,particle_rank)
Hout2_1.ThreeBody.SetMode('pn')    
Hout2_1.ThreeBody.Erase()
Hout2_1.ScaleOneBody(0.)
Hout2_1.TwoBody.Erase()

Hout2_2 = Operator(ms,0,0,0,particle_rank)
Hout2_2.ThreeBody.SetMode('pn')    
Hout2_2.ThreeBody.Erase()
Hout2_2.ScaleOneBody(0.)
Hout2_2.TwoBody.Erase()



print("\n   Two Body part: ")

# double commutator
# Commutator.comm232ss(Eta, Htemp, Hout2_1)

Commutator.comm223_232_Test(Eta, H, Hout2_1)


'''
Hout2_4 = Operator(ms,0,0,0,particle_rank)
Hout2_4.ThreeBody.SetMode('pn')    
Hout2_4.ThreeBody.Erase()
Hout2_4.ScaleOneBody(0.)
Hout2_4.TwoBody.Erase()

ReferenceImplementations.diagram_DIVa(Eta, H, Hout2_4)
print( "DIVa    {}    ".format(  Hout2_4.TwoBodyNorm() ))
Hout2_4.TwoBody.Erase()
ReferenceImplementations.diagram_DIVb(Eta, H, Hout2_4)
print( "DIVb    {}    ".format(  Hout2_4.TwoBodyNorm() ))
Hout2_4.TwoBody.Erase()
ReferenceImplementations.diagram_DIa(Eta, H, Hout2_4)
print( "DIa     {}    ".format(  Hout2_4.TwoBodyNorm() ))
Hout2_4.TwoBody.Erase()
ReferenceImplementations.diagram_DIb(Eta, H, Hout2_4)
print( "DIb     {}    ".format(  Hout2_4.TwoBodyNorm() ))
Hout2_4.TwoBody.Erase()

print("\n\n")

ReferenceImplementations.diagram_DIVb_intermediate(Eta, H, Hout2_4)
print( "diagram DIV_ragnar     {}    ".format(  Hout2_4.TwoBodyNorm() ))
print("\n\n")
'''

'''
print("\n  One Body part: ")
ReferenceImplementations.diagram_CIa(Eta, H, Hout2_4)
ReferenceImplementations.diagram_CIb(Eta, H, Hout2_4)
ReferenceImplementations.diagram_CIIa(Eta, H, Hout2_4)
ReferenceImplementations.diagram_CIIb(Eta, H, Hout2_4)
ReferenceImplementations.diagram_CIIc(Eta, H, Hout2_4)
ReferenceImplementations.diagram_CIId(Eta, H, Hout2_4)
ReferenceImplementations.diagram_CIIIa(Eta, H, Hout2_4)
ReferenceImplementations.diagram_CIIIb(Eta, H, Hout2_4)

print( "H4    {}    ".format(  Hout2_4.OneBodyNorm() ))
'''

# direct term
#Commutator.comm223_322_BruteForce(Eta, H, Hout2_2)
#ReferenceImplementations.comm223_232_BruteForce(Eta, H, Hout2_2)

Commutator.comm223_232_Factorization(Eta, H, Hout2_2)

print( "H1    {}    ".format(  Hout2_1.TwoBodyNorm()) )
print( "H2    {}    ".format(  Hout2_2.TwoBodyNorm() ))

print("\n\n")

H.PrintTimes()

