#!/usr/bin/python

from sys import argv
from numpy import sqrt

def main():
  fname = argv[1]
  ParseFile(fname)

def GetIsospinIndex(orbit_pn):
  Norb = 2*orbit_pn['n'] + orbit_pn['l']
  j2orb = orbit_pn['j']
  index_iso = (Norb*(Norb+1)+j2orb+1)/2 
  return index_iso

def ParseFile(fname):
  f = open(fname)
  line = f.readline()
  title = ''
  while 'Rank_J' not in line:
   title += line
   line = f.readline()
#  title  = f.readline()
#  while 'Rank_J'
  rankj  = int(line.split()[3])
  rankt  = int(f.readline().split()[3])
  parity = int(f.readline().split()[3])
  Zerobody = float(f.readline().split()[4])
  Opdata = {'rankj':rankj, 'rankt':rankt, 'parity':parity}
  line = f.readline()
  while 'index' not in line: line = f.readline()
  line = f.readline()
  orbits_pn = []
  jmax = -1
  while len(line.split())==6:
    ldat = line.split()
    index,n,l,j,tz = [int(x) for x in ldat[1:]]
    orbits_pn.append({'index':index,'n':n,'l':l,'j':j,'tz':-tz})
    jmax = max(jmax,j)
    line = f.readline()
    
  Opdata['jmax'] = jmax
  orbits_iso = []
  for orbit in orbits_pn:
    if orbit['tz'] > 0:
      ind=GetIsospinIndex(orbit)
      orbits_iso.append({'index':ind,'n':orbit['n'],'l':orbit['l'],'j':orbit['j']})
    else:
      if {'index':(orbit['index']+len(orbits_pn)/2)%len(orbits_pn),'n':orbit['n'],'l':orbit['l'],'j':orbit['j'],'tz':-orbit['tz']} not in orbits_pn:
        print "Error: proton and neutron orbits do not match!!!"
        return
  print title.rstrip('\n')
  print '!  Converted to isospin formalism with python script written by Ragnar Stroberg, so if this is wrong, blame him.'
  print '!  Rank_J :  %d'%(rankj)
  print '!  Rank_T :  %d'%(rankt)
  print '!  Parity :  %+d'%(parity)
  print '!  Zero body term: %10.6f'%(Zerobody)  # this should be zero...
  print '!  index   n   l   2j'
  for orbit in orbits_iso:
    print '!   %3d  %3d %3d %3d'%(orbit['index'],orbit['n'],orbit['l'],orbit['j'])  
  print '!\n!  a   b   c   d  Jab  Jcd  Tab  Tcd  Tzab  Tzcd  <ab Jab ||O|| cd Jcd>'
  

  tbmes = ReadIn(f,orbits_pn)
#  for t in sorted(tbmes):
#    print t,tbmes[t]
#  print 30*'-'
  Convert(tbmes,orbits_pn,orbits_iso,Opdata)


def GetLabelAndPhase(orbits_pn,a,b,c,d,Jab,Jcd):
  phase = 0
  if a>b:
    a,b = b,a
    phase += (orbits_pn[a-1]['j']+orbits_pn[b-1]['j'])/2 + Jab +1
  if c>d:
    c,d = d,c
    phase += (orbits_pn[c-1]['j']+orbits_pn[d-1]['j'])/2 + Jcd +1
  if a>c or a==c and b>d:
    a,b,c,d,Jab,Jcd = c,d,a,b,Jcd,Jab
    phase += Jab-Jcd
  label = '%02d-%02d-%02d-%02d-%02d-%02d'%(a,b,c,d,Jab,Jcd)
  return label,phase

def ReadIn(f,orbits_pn):
  uncoupled_TBMEs = {}
  for line in f:
    if '!' in line: continue
    ldat = line.split()
    a,b,c,d,Jab,Jcd = [int(x) for x in ldat[:6]]
    Op = float(ldat[6])
    label,phase = GetLabelAndPhase(orbits_pn,a,b,c,d,Jab,Jcd)
    if a==b: Op*= sqrt(2)   # I don't think this is right... but maybe it is.
    if c==d: Op*= sqrt(2)
    uncoupled_TBMEs[label] = Op*(-1)**phase
  return uncoupled_TBMEs


def GetMatEl(orbits_pn,TBMEs,a,b,c,d,Jab,Jcd):
  label,phase = GetLabelAndPhase(orbits_pn,a,b,c,d,Jab,Jcd)
  tbme = 0.0
  if label in TBMEs:
   tbme = TBMEs[label]
  return tbme * (-1)**(phase)


# For Petr's code protons have tz +1/2
# |T=1,Tz=0> = ( |pn> + |np> ) / sqrt(2)
# |T=0,Tz=0> = ( |pn> - |np> ) / sqrt(2)
#
def GetShifts(T,Tz):
  if (T,Tz) == (1,-1): return [1],[1],[1]
  if (T,Tz) == (1, 1): return [0],[0],[1]
  if (T,Tz) == (1, 0): return [0,1],[1,0],[sqrt(0.5),sqrt(0.5)]
  if (T,Tz) == (0, 0): return [0,1],[1,0],[sqrt(0.5),-sqrt(0.5)]


def Convert(uncoupled_TBMEs, orbits_pn,orbits_iso, Opdata):
  for Jab in range(Opdata['jmax']+1):
   for Jcd in range(Jab,Opdata['jmax']+1):
    if Jab+Jcd < Opdata['rankj'] or abs(Jab-Jcd)>Opdata['rankj']: continue
    for a in range(len(orbits_iso)):
     aiso = orbits_iso[a]['index']
     ja = orbits_iso[a]['j'] * 0.5
     for b in range(a,len(orbits_iso)):
      if (orbits_iso[a]['j'] + orbits_iso[b]['j'])/2 < Jab: continue
      if abs(orbits_iso[a]['j'] - orbits_iso[b]['j'])/2 > Jab: continue
      biso = orbits_iso[b]['index']
      jb = orbits_iso[b]['j'] * 0.5
      for c in range(0,len(orbits_iso)):
       ciso = orbits_iso[c]['index']
       jc = orbits_iso[c]['j'] * 0.5
       dmin = c
#       if a==c: dmin = b
       for d in range(dmin,len(orbits_iso)):
        if (a,b)==(c,d) and Jcd<Jab: continue
        if (orbits_iso[c]['j'] + orbits_iso[d]['j'])/2 < Jcd: continue
        if abs(orbits_iso[c]['j'] - orbits_iso[d]['j'])/2 > Jcd: continue
        jd = orbits_iso[d]['j'] * 0.5
        diso = orbits_iso[d]['index']
        for Tab in range(2):
         if a==b and (Jab+Tab)%2 < 1: continue
         for Tcd in range(2):
          if (a,b,Jab)==(c,d,Jcd) and Tab>Tcd: continue
#          if (Tab + Tcd < Opdata['rankt']) or (abs(Tab-Tcd)>Opdata['rankt']): continue
          if c==d and (Jcd+Tcd)%2 < 1: continue
          for Tzab in range(-Tab,Tab+1):
           for Tzcd in range(-Tcd,Tcd+1):
            if abs(Tzab-Tzcd)!=Opdata['rankt']: continue
            if (a,b,Jab,Tab)==(c,d,Jcd,Tcd) and Tzcd<Tzab: continue
            ap,an = a+1,a+1+len(orbits_pn)/2
            bp,bn = b+1,b+1+len(orbits_pn)/2
            cp,cn = c+1,c+1+len(orbits_pn)/2
            dp,dn = d+1,d+1+len(orbits_pn)/2
#            opval = 0
            if Tzab== 1: (ap,bp)=(an,bn)
            if Tzab==-1: (an,bn)=(ap,bp)
            if Tzcd== 1: (cp,dp)=(cn,dn)
            if Tzcd==-1: (cn,dn)=(cp,dp)
            opval = 0.5*( GetMatEl(orbits_pn,uncoupled_TBMEs,ap,bn,cp,dn,Jab,Jcd)  \
                        + GetMatEl(orbits_pn,uncoupled_TBMEs,bp,an,cp,dn,Jab,Jcd) * (-1)**(Tab+ja+jb-Jab)  \
                        + GetMatEl(orbits_pn,uncoupled_TBMEs,ap,bn,dp,cn,Jab,Jcd) * (-1)**(Tcd+jc+jd-Jcd)  \
                        + GetMatEl(orbits_pn,uncoupled_TBMEs,bp,an,dp,cn,Jab,Jcd) * (-1)**(Tab+Tcd+ja+jb+jc+jd-Jab-Jcd) )
              
#            shifts_a,shifts_b,factors_ab = GetShifts(Tab,Tzab)
#            shifts_c,shifts_d,factors_cd = GetShifts(Tcd,Tzcd)
#            for s_a,s_b,f_ab in zip(shifts_a,shifts_b,factors_ab):
#             apn = ( a + 1 + s_a*len(orbits_iso) )%(len(orbits_iso)*2)
#             bpn = ( b + 1 + s_b*len(orbits_iso) )%(len(orbits_iso)*2)
#             for s_c,s_d,f_cd in zip(shifts_c,shifts_d,factors_cd):
#              cpn = ( c + 1 + s_c*(len(orbits_iso) ))%(len(orbits_iso)*2)
#              dpn = ( d + 1 + s_d*(len(orbits_iso) ))%(len(orbits_iso)*2)
#              opval += GetMatEl(orbits_pn,uncoupled_TBMEs,apn,bpn,cpn,dpn,Jab,Jcd)*f_ab*f_cd
            if abs(opval)<1e-9: continue
            if aiso==biso: opval /= sqrt(2)
            if ciso==diso: opval /= sqrt(2)
#            opval /=2
            label = ' %3d %3d %3d %3d  %3d %4d %4d %4d %5d %5d    %16.9f'%(aiso,biso,ciso,diso,Jab,Jcd,Tab,Tcd,Tzab,Tzcd,opval)
            print label
        


if __name__ == '__main__':
  main()
