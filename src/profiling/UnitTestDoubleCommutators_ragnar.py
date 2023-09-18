#!/usr/bin/env python3

from pyIMSRG import *


emax=2
ref='C12'
val=ref

ms = ModelSpace(emax,ref,val)
def Get_pphh(O):
   Opphh = O * 0
   nch = ms.GetNumberTwoBodyChannels()
   for ch in range(nch):
      tbc = ms.GetTwoBodyChannel(ch)
      for ibra in tbc.GetKetIndex_pp() : 
         for iket in tbc.GetKetIndex_hh() : 
            tbme = O.TwoBody.GetTBMEnorm_chij(ch,ch,ibra,iket)
            Opphh.TwoBody.SetTBME_chij(ch,ch,ibra,iket,tbme)
   return Opphh

def Get_pppp(O):
   Opphh = O * 0
   nch = ms.GetNumberTwoBodyChannels()
   for ch in range(nch):
      tbc = ms.GetTwoBodyChannel(ch)
      for ibra in tbc.GetKetIndex_pp() : 
         for iket in tbc.GetKetIndex_pp() : 
            tbme = O.TwoBody.GetTBMEnorm_chij(ch,ch,ibra,iket)
            Opphh.TwoBody.SetTBME_chij(ch,ch,ibra,iket,tbme)
   return Opphh

def Get_hhhh(O):
   Opphh = O * 0
   nch = ms.GetNumberTwoBodyChannels()
   for ch in range(nch):
      tbc = ms.GetTwoBodyChannel(ch)
      for ibra in tbc.GetKetIndex_hh() : 
         for iket in tbc.GetKetIndex_hh() : 
            tbme = O.TwoBody.GetTBMEnorm_chij(ch,ch,ibra,iket)
            Opphh.TwoBody.SetTBME_chij(ch,ch,ibra,iket,tbme)
   return Opphh

def Get_phph(O):
   Opphh = O * 0
   nch = ms.GetNumberTwoBodyChannels()
   for ch in range(nch):
      tbc = ms.GetTwoBodyChannel(ch)
      for ibra in tbc.GetKetIndex_ph() : 
         for iket in tbc.GetKetIndex_ph() : 
            tbme = O.TwoBody.GetTBMEnorm_chij(ch,ch,ibra,iket)
            Opphh.TwoBody.SetTBME_chij(ch,ch,ibra,iket,tbme)
   return Opphh

def Get_ppph(O):
   Opphh = O * 0
   nch = ms.GetNumberTwoBodyChannels()
   for ch in range(nch):
      tbc = ms.GetTwoBodyChannel(ch)
      for ibra in tbc.GetKetIndex_pp() : 
         for iket in tbc.GetKetIndex_ph() : 
            tbme = O.TwoBody.GetTBMEnorm_chij(ch,ch,ibra,iket)
            Opphh.TwoBody.SetTBME_chij(ch,ch,ibra,iket,tbme)
   return Opphh

def Get_phhh(O):
   Opphh = O * 0
   nch = ms.GetNumberTwoBodyChannels()
   for ch in range(nch):
      tbc = ms.GetTwoBodyChannel(ch)
      for ibra in tbc.GetKetIndex_ph() : 
         for iket in tbc.GetKetIndex_hh() : 
            tbme = O.TwoBody.GetTBMEnorm_chij(ch,ch,ibra,iket)
            Opphh.TwoBody.SetTBME_chij(ch,ch,ibra,iket,tbme)
   return Opphh

def SelectJ(O,J):
   Opphh = O * 0
   nch = ms.GetNumberTwoBodyChannels()
   for ch in range(nch):
      tbc = ms.GetTwoBodyChannel(ch)
      if tbc.J != J: continue
      for ibra in tbc.GetKetIndex_pp() : 
         for iket in tbc.GetKetIndex_hh() : 
            tbme = O.TwoBody.GetTBMEnorm_chij(ch,ch,ibra,iket)
            Opphh.TwoBody.SetTBME_chij(ch,ch,ibra,iket,tbme)
   return Opphh

ut = UnitTest(ms)

jrank,parity,trank, particle_rank = 0,0,0,3
hermitian = +1
antihermitian = -1
X = ut.RandomOp(ms,jrank,parity,trank,particle_rank, antihermitian)
Y = ut.RandomOp(ms,jrank,parity,trank,particle_rank, hermitian)

gen = Generator()
gen.SetType('white')
#X = gen.GetHod(X)

#Xpphh = Get_pphh(X)
#X -= Get_pphh(X)
#Xrest = X-Xpphh
X = Get_pphh(X) #+ Get_ppph(X)
#Y = Get_phph(Y)
#Y -= Get_pphh(Y)
Y = Get_pphh(Y)
#Y = SelectJ(Y,2)
#Y -= SelectJ(Y,1)

Z = Operator(ms, jrank, parity,trank,particle_rank )
Z.ThreeBody.SetMode("pn")
Z.SetHermitian()

Zfull = Z*1.0;
ms.PreCalculateSixJ()

#Commutator.comm223ss(X,Y,Zfull)
#Commutator.comm231ss(X,Zfull,Zfull)
#Commutator.comm232ss(X,Zfull,Zfull)
ReferenceImplementations.comm223ss(X,Y,Zfull)
ReferenceImplementations.comm232ss(X,Zfull,Zfull)
#Commutator.comm223ss(X,Y,Zfull)
#Commutator.comm231ss(X,Zfull,Zfull)


##  Testing notes:
##   * Restricting Y=Ypppp selects IIa and IIIa.   This agrees.
##   * Restricting Y=Yhhhh selects IId and IIIb.   This agrees.
##   * Restricting Y=Yppph and X=Xpphh selects Ia IIa IIb IIIa
##   * Restricting X=Xpphh and Zpp selects IIa IIc IIIa IIIb.
##      
##


#ReferenceImplementations.diagram_CIa(X,Y,Z)
#ReferenceImplementations.diagram_CIb(X,Y,Z)
##
#ReferenceImplementations.diagram_CIIa(X,Y,Z)
#ReferenceImplementations.diagram_CIIb(X,Y,Z)
#ReferenceImplementations.diagram_CIIc(X,Y,Z)
#ReferenceImplementations.diagram_CIId(X,Y,Z)
##
#ReferenceImplementations.diagram_CIIIa(X,Y,Z)
#ReferenceImplementations.diagram_CIIIb(X,Y,Z)



#ReferenceImplementations.diagram_DIa(X,Y,Z)
#ReferenceImplementations.diagram_DIb(X,Y,Z)
#
#ReferenceImplementations.diagram_DIVa(X,Y,Z)
#ReferenceImplementations.diagram_DIVb(X,Y,Z)
ReferenceImplementations.diagram_DIVb_intermediate(X,Y,Z)

Zfull = Get_pphh(Zfull)
Z = Get_pphh(Z)
### Checking the 1b part
#print('Z1b is\n',Z.OneBody)
#print('[231] =\n',Zfull.OneBody)
#diff = (Z.OneBody-Zfull.OneBody).Norm()


nch = ms.GetNumberTwoBodyChannels()
#for ch in range(nch):
#for ch in range(29):
for ch in [9]:
   tbc = ms.GetTwoBodyChannel(ch)
   print(40*'=')
   print('channel {} JPT= {} {} {}'.format(ch,tbc.J,tbc.parity,tbc.Tz))
   nkets = tbc.GetNumberKets()
   for iket in range(nkets):
      ket = tbc.GetKet(iket)
      op = ms.GetOrbit(ket.p)
      oq = ms.GetOrbit(ket.q)
      print('{} : ({},{}), ;   {}{}'.format(iket,ket.p,ket.q, op.cvq,oq.cvq ))
   mat_diagram = Z.TwoBody.GetChannelMatrix(tbc.J,tbc.parity,tbc.Tz)
   mat_full = Zfull.TwoBody.GetChannelMatrix(tbc.J,tbc.parity,tbc.Tz)
   print('diagrams:\n',mat_diagram)
   print('full:\n',mat_full)
   diff = (mat_diagram-mat_full).Norm()
   print('norms,',mat_diagram.Norm(),' ',mat_full.Norm(),' difference',diff )
   if (abs(diff)>1e-8): print('  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
#   print('BTW, X and Y look like this')
#   X.TwoBody.PrintMatrix(ch,ch)
#   print()
#   Y.TwoBody.PrintMatrix(ch,ch)

   print('\n')


print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
print('norms:',Z.TwoBody.Norm(),Zfull.TwoBody.Norm())
diff = (Z.TwoBody-Zfull.TwoBody).Norm()
print('difference = ',diff)
if abs(diff)<1e-8:
    print('HOORAY!!!')
else:
    print("NNNOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO!!!")



#mschemeresult = TestTerm(X,Y,ms)
#print('The mscheme result is ',mschemeresult,'->',sum(mschemeresult))

prof = IMSRGProfiler()
prof.PrintAll()
