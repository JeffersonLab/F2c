#!/usr/bin/env python
import sys,os
import numpy as np
import pandas as pd
import lhapdf

#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('text',usetex=True)
import pylab as py

#--local
from tools import load, save,checkdir
from mstw import hqlib as hq

fname = 'EPPS16nlo_CT14nlo_Fe56'

def gen_eic_table():
    t1=pd.read_excel('expdata/30000.xlsx')
    t2=pd.read_excel('expdata/30001.xlsx')
    t3=pd.read_excel('expdata/30002.xlsx')
    t4=pd.read_excel('expdata/30003.xlsx')
    t=pd.concat([t1,t2,t3,t4], ignore_index=True)
    t.to_excel('expdata/eic.xlsx') 

def set_mstw():

    iord  = 2        #--0=LO, 1=NLO, 2=NNLO
    ipn   = 1        #--1 = proton, 2 = neutron, 3 = nucleus DSSZ, 4 = nucleus EPS09 
    fr2   = 1.0      #--ratio of mu_f^2 to mu_r^
    mur   = 1.       #--input mu_r in G
    asmur = 0.49128  #--input value of alpha_s at mu_r
    mc    = 1.4      #--charm quark mass
    mb    = 4.75     #--bottom quark mass
    mt    = 175.     #--top quark mas          


    hq.mstwcommon.fr2     = fr2     
    hq.mstwcommon.mur     = mur     
    hq.mstwcommon.asmur   = asmur 
    hq.mstwcommon.mcharm  = mc      
    hq.mstwcommon.mbottom = mb      
    
    hq.iordcommon.iord=iord
    
    hq.wate96()
    hq.initalphas(iord, fr2, mur, asmur, mc, mb, mt)

def gen_predictions():

    #--setup nlo code
    set_mstw()

    #--setup pdfs
    hq.setup(fname.ljust(100),0)
    nsets=97    

    #--read kinematics
    t=pd.read_excel('expdata/eic.xlsx')
    t=t.to_dict(orient='list')
    npts=len(t['X'])  

    #--compute F2c accross kinematics
    for iset in range(nsets):
        hq.setup_pdfset(iset)

        t[iset]=[] #--create empty array to fill theory 

        for i in range(npts):
            x=t['X'][i]
            Q2=t['Q2'][i]
            Q = np.sqrt(Q2)
            #,f2,f2c,f2b,fl,flc,flb
            F2c=hq.mstwnc(x,Q,1)[1]
            t[iset].append(F2c)

    #--convert all arrays to numpy arrays
    for _ in t: t[_]=np.array(t[_])

    checkdir('data')
    save(t,'data/predictions.po')

def gen_hess_error():
    t=load('data/predictions.po')
    nsets=97    
    err2=np.zeros(t[0].size)
    for iset in range(1,nsets,2):
        err2+=(t[iset]-t[iset+1])**2/4.0
    t['hess']=err2**0.5
    save(t,'data/hess.po')

def gen_hess_glue(Q2=2.0):

    X=10**np.linspace(-4,-1,100)
    X=np.append(X,np.linspace(0.1,0.99,100))

    data={'X':X,'Q2':Q2}
    pdfs=lhapdf.mkPDFs(fname)
    for i in range(len(pdfs)):
        data[i]=np.array([pdfs[i].xfxQ2(21,x,Q2) for x in X])

    err2=np.zeros(data[0].size)
    for iset in range(1,len(pdfs),2):
        err2+=(data[iset]-data[iset+1])**2/4.0
    data['err']=err2**0.5
    save(data,'data/hess_glue_%.2f.po'%Q2)

def gen_rand(nrep=1000):
    """
    random pdfs will be constructed from eigen directions as
    f_k = f_0 + sum_i rand_i (f[i]-f[-i])/2 
    similarly the random  F2c will computed as
    F2c_k = F2c_0 + sum_i rand_i (F2c[i]-F2c[-i])/2 
    here we generate the rand_i
    """
    nsets=97
    rand=np.random.randn(nrep,(nsets-1)/2)
    save(rand,'data/rand.po')

def gen_mc_F2c():
    rnd=load('data/rand.po')
    t=load('data/predictions.po')

    #--gen MCF2c
    nsets=97    
    mcF2c=[]
    for k in range(len(rnd)):
        F2c=np.copy(t[0])
        cnt=0
        for iset in range(1,nsets,2):
            F2c+=rnd[k][cnt]*(t[iset]-t[iset+1])/2
            cnt+=1
        mcF2c.append(F2c)
    mcF2c=np.array(mcF2c)
    save(mcF2c,'data/mcF2c.po')

def gen_mc_glue(Q2=2.0):
    rnd=load('data/rand.po')
    t=load('data/hess_glue_%.2f.po'%Q2)

    #--gen MCF2c
    nsets=97    
    mcglue=[]
    for k in range(len(rnd)):
        g=np.copy(t[0])
        cnt=0
        for iset in range(1,nsets,2):
            g+=rnd[k][cnt]*(t[iset]-t[iset+1])/2
            cnt+=1
        mcglue.append(g)
    mcglue=np.array(mcglue)
    save(mcglue,'data/mcglue_%.2f.po'%Q2)

def gen_weights():
    rnd=load('data/rand.po')
    t=load('data/predictions.po')
    mcF2c=load('data/mcF2c.po')

    #--get absolute simulated errors
    t['alpha']=t[0]*t['relerr']

    #--gen chi2
    chi2=[]
    for F2c in  mcF2c:
        exp=np.copy(t[0])
        res=(exp-F2c)/t['alpha']
        chi2.append(np.sum(res**2))
    chi2=np.array(chi2)
    dchi2=chi2-t[0].size

    #--gen weights         
    weights=np.exp(-0.5*dchi2)
    norm=np.sum(weights)
    weights/=norm

    save(weights,'data/weights.po')
  
def plot1(Q2=2.0):  
    hess = load('data/hess_glue_%.2f.po'%Q2)
    mc   = load('data/mcglue_%.2f.po'%Q2)

    ncols,nrows=1,1
    py.figure(figsize=(ncols*5,nrows*4))
    ax=py.subplot(nrows,ncols,1)

    X=hess['X']

    #--plot current errorbands (hess version)
    ax.fill_between(X,(hess[0]-hess['err'])/hess[0]
                     ,(hess[0]+hess['err'])/hess[0]
                     ,color='Y',alpha=0.5
                     ,label=r'$\rm EPPS16$')

    #--plot current errorbands (mc version)
    #g  = np.mean(mc,axis=0)
    #dg = np.std(mc,axis=0)
    #ax.fill_between(X,(g-dg)/g
    #                 ,(g+dg)/g
    #                 ,facecolor='none',hatch='//',edgecolor='k')
        

    #--plot mc reweighted errorbands
    weights=load('data/weights.po')
    g=np.sum([weights[i]*mc[i] for i in range(weights.size)],axis=0)
    dg=np.sum([weights[i]*(mc[i]-g)**2 for i in range(weights.size)],axis=0)**0.5


    ax.fill_between(X,(g-dg)/g
                     ,(g+dg)/g
                     ,facecolor='r',alpha=0.8
                     ,label=r'$\rm EPPS16+EIC$')

    #--plot eic kinematics
    t=pd.read_excel('expdata/eic.xlsx')
    xmin=np.amin(t.X)
    xmax=np.amax(t.X)
    ax.plot([xmin,xmax],[0.6,0.6],'g-',lw=5)
    #ax.text(0.3,0.15,r'$\rm EIC~kinematics$',size=10,transform=ax.transAxes)

    #--grids prop
    ax.axhline(1,color='k',ls=':')
    ax.set_ylim(0.5,1.5)
    ax.set_xlim(1e-2,0.9)
    ax.semilogx()
    #ax.set_yticks([0.8,0.9,1,1.1,1.2])
    ax.set_xticks([0.01,0.1])
    ax.set_xticklabels([r'$0.01$',r'$0.1$'])


    ax.legend(loc=1,framealpha=1)

    #--axis labels
    ax.set_ylabel(r'$\delta g/ g$',size=20)

    ax.set_xlabel(r'$x$',size=20)
    ax.xaxis.set_label_coords(0.95,-0.03)

    ax.text(0.3,0.9,r'$\mu^2=2{~\rm GeV^2}$',size=15,transform=ax.transAxes)
    ax.text(0.3,0.8,r'$A=56$',size=15,transform=ax.transAxes)


    py.tight_layout()
    checkdir('gallery')
    #py.savefig('gallery/glue.pdf')
    py.savefig('gallery/glue.png')
    py.close()

def plot2(Q2=2.0):  

    hess = load('data/hess_glue_%.2f.po'%Q2)
    mc   = load('data/mcglue_%.2f.po'%Q2)


    ncols,nrows=1,1
    py.figure(figsize=(ncols*5,nrows*4))
    ax=py.subplot(nrows,ncols,1)

    X=hess['X']

    ct14  = lhapdf.mkPDFs('CT14nlo')
    denom = np.array([ct14[0].xfxQ2(21,x,Q2) for x in X])


    #--plot current errorbands (hess version)
    ax.fill_between(X,(hess[0]-hess['err'])/denom
                     ,(hess[0]+hess['err'])/denom
                     ,color='Y',alpha=0.5
                     ,label=r'$\rm EPPS16$')

    #--plot current errorbands (mc version)
    #g  = np.mean(mc,axis=0)
    #dg = np.std(mc,axis=0)
    #ax.fill_between(X,(g-dg)/denom
    #                 ,(g+dg)/denom
    #                 ,facecolor='none',hatch='//',edgecolor='k')
        

    #--plot mc reweighted errorbands
    weights=load('data/weights.po')
    g=np.sum([weights[i]*mc[i] for i in range(weights.size)],axis=0)
    dg=np.sum([weights[i]*(mc[i]-g)**2 for i in range(weights.size)],axis=0)**0.5


    ax.fill_between(X,(g-dg)/denom
                     ,(g+dg)/denom
                     ,facecolor='r',alpha=0.8
                     ,label=r'$\rm EPPS16+EIC$')

    ##--plot eic kinematics
    t=pd.read_excel('expdata/eic.xlsx')
    xmin=np.amin(t.X)
    xmax=np.amax(t.X)
    ax.plot([xmin,xmax],[0.25,0.25],'g-',lw=5)
    #ax.text(0.3,0.15,r'$\rm EIC~kinematics$',size=10,transform=ax.transAxes)

    #--grids prop
    ax.axhline(1,color='k',ls=':')
    ax.set_ylim(0,2)
    ax.set_xlim(1e-2,0.9)
    ax.semilogx()
    ax.set_yticks([0,0.5,1,1.5,2])
    ax.set_xticks([0.01,0.1])
    ax.set_xticklabels([r'$0.01$',r'$0.1$'])

    ax.text(0.05,0.9,r'$\mu^2=2{~\rm GeV^2}$',size=15,transform=ax.transAxes)
    ax.text(0.05,0.8,r'$A=56$',size=15,transform=ax.transAxes)

    ax.legend(loc=1,framealpha=1)

    #--axis labels
    ax.set_ylabel(r'$R_g$',size=20)

    ax.set_xlabel(r'$x$',size=20)
    ax.xaxis.set_label_coords(0.95,-0.03)



    py.tight_layout()
    checkdir('gallery')
    #py.savefig('gallery/Rg.pdf')
    py.savefig('gallery/Rg.png')
    py.close()

if __name__=="__main__":

    Q2=2.0

    #gen_predictions()
    #gen_hess_error()
    #gen_hess_glue(Q2)
    #gen_rand(nrep=10000)
    #gen_mc_F2c()
    #gen_mc_glue(Q2)
    #gen_weights()
    plot1(Q2)
    plot2(Q2)














