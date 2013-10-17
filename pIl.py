#!/usr/bin/python
"""\
calculates the pI of the protein

usage pIl.py xxx.pka 

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
    gamma = {'ASP': -1.0, 'GLU':-1.0, 'C-':-1.0, 'HIS':1.0, 'CYS':-1.0, 'TYR':-1.0, 'N+':1.0, 'LYS':1.0, 'ARG':1.0, 'Npl':1.0, 'N3':1.0, 'Nar':1.0, 'Oco':-1.0} 
   
    x = gamma[resnam]*(pKa - pH)
    y = 10**x
    rescharge = gamma[resnam]*(y/(1.0+y))
    return rescharge

def getA(pkinfo):
    A = 0.0
    acid = {'ASP': 1.0, 'GLU':1.0, 'C-':1.0, 'HIS':0.0, 'CYS':0.0, 'TYR':0.0, 'N+':0.0, 'LYS':0.0, 'ARG':0.0, 'Npl':0.0, 'N3':0.0, 'Nar':0.0, 'Oco':1.0}

    for words in pkinfo:
       if len(words) != 5: resnam = words[0] 
       if len(words) == 5: resnam = words[4] 
       A = A + acid[resnam]
    return A

def getB(pkinfo):
    B = 0.0
    base = {'ASP': 0.0, 'GLU':0.0, 'C-':0.0, 'HIS':0.0, 'CYS':0.0, 'TYR':0.0, 'N+':1.0, 'LYS':1.0, 'ARG':1.0, 'Npl':1.0, 'N3':1.0, 'Nar':1.0, 'Oco':0.0}

    for words in pkinfo:
       if len(words) != 5: resnam = words[0] 
       if len(words) == 5: resnam = words[4] 
       B = B + base[resnam]
    return B
   

def main():

   filename = sys.argv[1]

   pkinfo = get_pk(filename)

   pHprot = 7.0
   pHmod = 7.0

   qmax = 9999.0
   print 'Protein charge of folded and unfolded state as a function of pH'
   for pH in range(14):
       qp = Qprot(pkinfo,pH)
       print '%5.2f %5.2f %5.2f' % (pH, qp, Qmod(pkinfo,pH) )
       if abs(qp) < qmax:
          qmax = qp
          pHprot = pH
   
   qmax = 9999.0
   for pH in range(14):
       qp = Qmod(pkinfo,pH)
       if abs(qp) < qmax:
          qmax = qp
          pHmod = pH

   pIprot = Qprot(pkinfo,pHprot)
   pImod = Qmod(pkinfo,pHmod)


   n = 4.0
   count = 0
   while abs(pIprot) > 0.01:
       count = count + 1
       if count == 200:
	  print 'Folded pI not converged in 200 steps'
	  print '%5.2f' % pIprot
	  break
       if abs(pIprot) > n:
          pIprot = pIprot*n/abs(pIprot)
          n = n/2.0
       pHprot = pHprot + pIprot
       pIprot = Qprot(pkinfo,pHprot) 

   n = 4.0
   count = 0
   while abs(pImod) > 0.01:
       count = count + 1
       if count == 200:
	  print 'Unfolded pI not converged in 200 steps'
	  break
       if abs(pImod) > n:
          pImod = pImod*n/abs(pImod)
          n = n/2.0
       pHmod = pHmod + pImod
       pImod = Qmod(pkinfo,pHmod) 

   print 'The pI is %5.2f (folded) and %5.2f (unfolded)' % (pHprot,pHmod)

   BFprot = (Qprot(pkinfo,pHprot+0.1) - Qprot(pkinfo,pHprot))/0.1
   BFmod = (Qmod(pkinfo,pHmod+0.1) - Qmod(pkinfo,pHmod))/0.1

   print 'The buffering capacity (dQ/dpH) is %5.2f (folded) and %5.2f (unfolded)' % (BFprot,BFmod)

   print 'The charge of the unfolded protein at folded pI is %5.2f' % (Qmod(pkinfo,pHprot))

   A = getA(pkinfo)
   B = getB(pkinfo)
   R = A/B

   pKmod = 4.2
   if R < 1.0:
      pKmod = 11.2

   pIR = pKmod - math.log10(R)

   print 'Patrickios/Yamasaki pI is %5.2f for %3i acidic and %3i basic residues' % (pIR, A, B) 
   print 'CS Patrickios, EN Yamasaki, Analytical Biochemistry 1995, 231, 82-81'

if __name__ == '__main__': main() 
