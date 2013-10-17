#!/usr/bin/python
"""\
calculates pH dependent free energy

usage dgl.py xxx.pka 

This gives the same results at dg.py but without the integration
"""

import string,re,sys,os,math

def get_pk(filename):
    input = open(filename,'r')

    pkinfo = []
    for line in input.readlines():
        words = string.split(line)
        pkinfo.append(words)
    return pkinfo
        
def Qmod(pkinfo,pH):
    Qmod = 0.0
    
    for words in pkinfo:
       if len(words) != 5: resnam = words[0] 
       if len(words) == 5: resnam = words[4] 
       pKa = eval(words[3]) 
       Qmod = Qmod + rescharge(resnam,pH,pKa)
    return Qmod

def Qprot(pkinfo,pH):
    Qprot = 0.0
    
    for words in pkinfo:
       if len(words) != 5: resnam = words[0] 
       if len(words) == 5: resnam = words[4] 
       pKa = eval(words[2])
       Qprot = Qprot + rescharge(resnam,pH,pKa)
    return Qprot

def rescharge(resnam,pH,pKa):
   
    x =  pH-pKa
    y = 10**x
    rescharge = math.log10(1+y)
    return rescharge

def main():

   filename = sys.argv[1]

   pkinfo = get_pk(filename)

   print 'Free energy of unfolding (kcal/mol) as a function of pH'

   dg = 0.0 
   dgmin = 9999999.0
   pHmin = 0.0
   for ipH in range (0, 140):
       pH = ipH/10.0
       dq = -Qprot(pkinfo,pH) + Qmod(pkinfo,pH)
       dg = 1.36*dq
       if (pH%1 == 0): print '%5.2f %5.2f' % (pH, -dg)
       if dg < dgmin:
          dgmin = dg
          pHmin = pH

   print ' '
   print 'The pH of optimum stability is %4.1f for which the free energy is %5.1f kcal/mol at 298K' % (pHmin, -dgmin)

   pH = pHmin
   pHlow = pHmin
   pHlowrange = 0.0
   dg = 0.0
   dg80 = 0.8*abs(dgmin)
   for idpH in range (0, 100):
       pH = pH - 0.1
       if pH < 0.0: break
       dq = Qprot(pkinfo,pH)-Qmod(pkinfo,pH)
       dg = 1.36*dq
       if abs(dg) > dg80:
          pHlow = pH
       if dg > 0.0:
          pHlowrange = pH

   pH = pHmin
   pHhigh = pHmin
   pHhighrange = 14.0
   dg = 0.0
   for idpH in range (0, 100):
       pH = pH + 0.1
       if pH > 13.0: break
       dq = Qprot(pkinfo,pH)-Qmod(pkinfo,pH)
       dg = 1.36*dq
       if abs(dg) > dg80:
          pHhigh = pH
       if dg > 0.0:
          pHhighrange = pH


   print 'The free energy is within 80 percent of maximum range %4.1f - %4.1f' % (pHlow, pHhigh)

   print 'The free energy is positive in the range %4.1f - %4.1f' % (pHlowrange, pHhighrange)
   print ' '

if __name__ == '__main__': main() 
