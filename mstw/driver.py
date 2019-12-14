#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
#import lhapdf
import hqlib as hq

iord  = 2        #--0=LO, 1=NLO, 2=NNLO
ipn   = 1        #--1 = proton, 2 = neutron, 3 = nucleus DSSZ, 4 = nucleus EPS09 
fr2   = 1.0      #--ratio of mu_f^2 to mu_r^
mur   = 1.       #--input mu_r in G
asmur = 0.49128  #--input value of alpha_s at mu_r
mc    = 1.4      #--charm quark mass
mb    = 4.75     #--bottom quark mass
mt    = 175.     #--top quark mas          


#print(dir(hq.mstwcommon))
#sys.exit()
hq.mstwcommon.fr2     = fr2     
hq.mstwcommon.mur     = mur     
hq.mstwcommon.asmur   = asmur 
hq.mstwcommon.mcharm  = mc      
hq.mstwcommon.mbottom = mb      

hq.iordcommon.iord=iord

hq.wate96()
hq.initalphas(iord, fr2, mur, asmur, mc, mb, mt)

fname = 'JAM19PDF_proton_nlo'
fname = 'EPPS16nlo_CT14nlo_He4'
#fname = 'CT10nlo'
iset  = 0 
hq.setup(fname.ljust(100),iset)


x  = 0.5
q2 = 10.0
q = np.sqrt(q2)
#,f2,f2c,f2b,fl,flc,flb
data=hq.mstwnc(x,q,ipn)
print('%20.10f'%hq.alphas(q))
for _ in data: print('%20.10f'%_)


























