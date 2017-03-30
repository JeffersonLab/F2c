#!/usr/bin/env python
import sys,os
import numpy as np
from BAYESR import BAYESR
#################################

#--step 0: set the paths for the experimental data 
conf={}
conf['datasets']={}
#conf['datasets'][20000]='expdata/20000.xlsx'
#conf['datasets'][20001]='expdata/20001.xlsx'
#conf['datasets'][20002]='expdata/20002.xlsx'
#conf['datasets'][20003]='expdata/20003.xlsx'

conf['datasets'][30000]='expdata/30000.xlsx'
conf['datasets'][30001]='expdata/30001.xlsx'
conf['datasets'][30002]='expdata/30002.xlsx'
conf['datasets'][30003]='expdata/30003.xlsx'
conf['gauss noise']=False

#-- case selector
#conf['case']=0  # data sets 2000X
conf['case']=1  # data sets 3000X


#--step 1: initalize the BAYESR class
br=BAYESR(conf)

#--step 2: compute F2C for the data sets, gen MC F2C and get weights
br.gen_F2C_input()
br.get_weights()

#--step 3: select kinematics to visualize impact
#--choose X
X1=10**np.linspace(-2,-1,100)
X2=np.linspace(0.1001,0.4,100)
X=np.append(X1,X2)
#--choose Q2
Q2=np.ones(X.size)*5.0
#--gen MC F2C and do reweighting
br.gen_priors_and_posteriors(X,'F2C')
br.gen_priors_and_posteriors(X,'gl')
#--make plot
br.plot_F2C()
br.plot_gl()



