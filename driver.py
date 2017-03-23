#!/usr/bin/env python
import sys,os
import numpy as np
from BAYESR import BAYESR
#################################

#--step 0: set the paths for the experimental data 
conf={}
conf['datasets']={}
conf['datasets'][10000]='expdata/10000.xlsx'
conf['gauss noise']=True


#--step 1: initalize the BAYESR class
br=BAYESR(conf)

#--step 2: compute F2C for the data sets, gen MC F2C and get weights
br.gen_F2C_input()
br.get_weights()

#--step 3: select kinematics to visualize impact
#--choose a name for this analysis
fname='test'

#--choose X
X1=10**np.linspace(-3,-1,100)
X2=np.linspace(0.1001,0.4,100)
X=np.append(X1,X2)
#--choose Q2
Q2=np.ones(X.size)*5.0
#--gen MC F2C and do reweighting
br.gen_priors_and_posteriors(X,Q2,fname)
#--make plot
br.plot_F2C(fname)
br.plot_gl(fname)



