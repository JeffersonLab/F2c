#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import pandas as pd

Input=sys.argv[1]
Output=sys.argv[2]

L=open(Input).readlines()
L=[l.strip() for l in L]
L=[l for l in L if l!='']

COLS=[l for l in L if l.startswith('#')]
TAB=[[float(x) for x in l.split()] for l in L if l.startswith('#')==False]
TAB=np.transpose(TAB)

D={}
D['X']=TAB[0]
D['Q2']=TAB[1]
D['relerr']=TAB[2]

D=pd.DataFrame(D)
writer = pd.ExcelWriter(Output)
D.to_excel(writer, index=False,header=True)

