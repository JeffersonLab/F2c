#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import pandas as pd
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
#################################
from F2C.F2C import F2C
from tools import BAR,save,load,checkdir,load_config,tex
#################################

class BAYESR: 

  def __init__(self,conf):
    self.conf=conf
    self.f2c=F2C('F2C/')
    checkdir('data')
    checkdir('gallery')

  def gen_F2C(self,X,Q2,fname,isave=True):
    """
    Notes
    pdfset : 1     = central fit
             2,3   = +,- 
             ...   = ...
             30,31 = +,-
    isf: 1=F2H  2=FLH
    """
    isf=1
    D={}
    self.f2c.set_pdfset(1)
    D[0]=np.array([self.f2c.get_F2C(X[j],Q2[j],isf) for j in range(X.size)])
    bar=BAR('calc F2C',15)
    for i in range(1,16):
      k=2*i
      self.f2c.set_pdfset(k)
      D[i]=np.array([self.f2c.get_F2C(X[j],Q2[j],isf) for j in range(X.size)])
      k=2*i+1
      self.f2c.set_pdfset(k)
      D[-i]=np.array([self.f2c.get_F2C(X[j],Q2[j],isf) for j in range(X.size)])
      bar.next()
    bar.finish()
    data={'X':X,'Q2':Q2,'F2C':D}
    if isave: 
      save(data,'data/'+fname)
      print 'F2C has been saved at data/%s'%fname
    return data

  def gen_XPDF(self,X,Q2,flav):
    isf=1
    D={}
    self.f2c.set_pdfset(1)
    D[0]=np.array([self.f2c.get_xpdf(X[j],Q2[j],flav) for j in range(X.size)])
    bar=BAR('calc PDF=%s'%flav,15)
    for i in range(1,16):
      k=2*i
      self.f2c.set_pdfset(k)
      D[i]=np.array([self.f2c.get_xpdf(X[j],Q2[j],flav) for j in range(X.size)])
      k=2*i+1
      self.f2c.set_pdfset(k)
      D[-i]=np.array([self.f2c.get_xpdf(X[j],Q2[j],flav) for j in range(X.size)])
      bar.next()
    bar.finish()
    data={'X':X,'Q2':Q2,flav:D}
    return data

  def gen_F2C_input(self):
    conf=self.conf
    for k in conf['datasets']:
      tab=pd.read_excel(conf['datasets'][k])
      data=self.gen_F2C(tab.X.values,tab.Q2.values,str(k),isave=False)
      data['dF2C(exp)']=data['F2C'][0]*tab.relerr.values
      if conf['gauss noise']:
        data['F2C(exp)']=data['F2C'][0]+np.random.randn(tab.X.size)*data['dF2C(exp)']
      else:
        data['F2C(exp)']=np.copy(data['F2C'][0])
      save(data,'data/%s'%str(k))
      print 'F2C has been saved at data/%s'%str(k)

  def get_hess_errors(self,D,key):
    err2=np.zeros(D[key][0].size)
    for k in range(1,16):
      err2+=(D[key][k]-D[key][-k])**2/4.0
    return {'d%s(hess)'%key:err2**0.5}

  def print_F2C(self,fname):
    D=load(fname)
    self.get_hess_errors(D,'F2C')
    fmt='x=%5.3f  Q2=%10.2e  F2C=%10.3e  dF2C(hess)/F2C=%10.3e dF2C(exp)/F2C=%10.3e'
    for i in range(D['X'].size):
      print fmt%(D['X'][i],D['Q2'][i],D['F2C'][0][i],\
                 D['dF2C(hess)'][i]/D['F2C'][0][i],
                 D['dF2C(exp)'][i]/D['F2C'][0][i])

  def gen_rand(self,nrep):
    RAND={}
    for k in range(nrep):
      RAND[k]={}
      for s in range(1,16):
        RAND[k][s]=np.random.randn(1)
    return RAND

  def gen_mc(self,data,RAND,nrep=1000):
    MCF2C=[]
    for k in range(nrep):
      F2Ck=np.copy(data[0])
      for s in range(1,16):
        F2Ck+=(data[s]-data[-s])*RAND[k][s]/2
      MCF2C.append(F2Ck)
    return np.array(MCF2C)

  def _get_weights(self,CHI2):
    chi2min=np.amin(CHI2)
    dCHI2=CHI2-chi2min
    W=np.exp(-0.5*dCHI2)
    wnorm=np.sum(W)
    return W/wnorm

  def get_weights(self):
    nrep=10000
    D={}
    D['RAND']=self.gen_rand(nrep)
    conf=self.conf
    D['Widx']={}
    D['exp']={}
    CHI2TOT=np.zeros(nrep)
    bar=BAR('compute bayes weights',len(conf['datasets'].keys()))
    for idx in conf['datasets']:
      d=load('data/%d'%idx)
      MCF2C=self.gen_mc(d['F2C'],D['RAND'],nrep=nrep)
      CHI2=np.array([np.sum(((d['F2C(exp)']-F2CK)/d['dF2C(exp)'])**2) for F2CK in MCF2C])
      D['Widx'][idx]=self._get_weights(CHI2)
      D['exp'][idx]=d
      CHI2TOT+=CHI2
      bar.next()
    bar.finish()
    D['Wtot']=self._get_weights(CHI2TOT)
    print 'saving bayes weights...'
    save(D,'data/BayesWeights')

  def get_priors(self,D,key):
    mc=self.gen_mc(D[key],D['RAND'],nrep=D['Wtot'].size)
    mean=np.mean(mc,axis=0)
    std=np.std(mc,axis=0)
    return {'MC%s'%key:mc,'%s(pri)'%key:mean,'d%s(pri)'%key:std}

  def get_posteriors(self,D,key):
    mean=np.einsum('i,ij',D['Wtot'],D['MC%s'%key])
    std=np.einsum('i,ij',D['Wtot'],(D['MC%s'%key]-mean)**2)**0.5
    return {'%s(pos)'%key:mean,'d%s(pos)'%key:std}

  def gen_priors_and_posteriors(self,X,Q2,fname):
    D=load('data/BayesWeights')

    D.update(self.gen_F2C(X,Q2,fname,isave=False))
    D.update(self.get_priors(D,'F2C'))
    D.update(self.get_posteriors(D,'F2C'))
    D.update(self.get_hess_errors(D,'F2C'))

    D.update(self.gen_XPDF(X,Q2,'gl'))
    D.update(self.get_priors(D,'gl'))
    D.update(self.get_posteriors(D,'gl'))
    D.update(self.get_hess_errors(D,'gl'))

    for k in D:
      if k.startswith('MC'):
        D[k]=D[k][:50]

    print 'saving reweighted data at data/%s'%fname
    save(D,'data/%s'%fname)

  def plot_F2C(self,fname):  
    D=load('data/%s'%fname)
    for i in range(D['F2C'][0].size):
      if D['F2C'][0][i]==0: I=i-1;break
    X=D['X']
    norm=D['F2C'][0] 

    ax=py.subplot(111)
    Lhess,=ax.plot(X[:I],(D['F2C'][0]/norm)[:I],'r--')
    DO=(D['F2C'][0]-D['dF2C(hess)'])/norm
    UP=(D['F2C'][0]+D['dF2C(hess)'])/norm
    Bhess=ax.fill_between(X[:I],DO[:I],UP[:I],alpha=0.3,color='y',zorder=1)

    Lpri,=ax.plot(X[:I],(D['F2C(pri)']/norm)[:I],'k:')
    DO=(D['F2C(pri)']-D['dF2C(pri)'])/norm
    UP=(D['F2C(pri)']+D['dF2C(pri)'])/norm
    Bpri=ax.fill_between(X[:I],DO[:I],UP[:I],alpha=0.3,\
      facecolor='none',edgecolor='k',hatch='...',zorder=10)

    Lpos,=ax.plot(X[:I],(D['F2C(pos)']/norm)[:I],'b--')
    DO=(D['F2C(pos)']-D['dF2C(pos)'])/norm
    UP=(D['F2C(pos)']+D['dF2C(pos)'])/norm
    Bpos=ax.fill_between(X[:I],DO[:I],UP[:I],alpha=0.3,color='b',zorder=10)
  
    for k in range(30):
      Y=D['MCF2C'][k]/norm
      Lmc,=ax.plot(X[:I],Y[:I],'r-',alpha=0.3,zorder=0)

    for idx in D['exp']:
      d=D['exp'][idx]
      norm=d['F2C'][0]
      Lexp=ax.errorbar(d['X'],d['F2C(exp)']/norm,\
          yerr=d['dF2C(exp)']/norm,fmt='k.',capsize=0)
  
    L=[tex('Hess'),tex('pri'),tex('pri(MC)'),tex('pos'),tex('sim~dat')]
    H=[(Bhess,Lhess),(Bpri,Lpri),Lmc,(Bpos,Lpos),Lexp]
    ax.legend(H,L,loc=2,frameon=0,fontsize=15,ncol=2)

    ax.set_xlabel(r'$x$',size=20)
    ax.set_ylabel(r'$F_2^c/F_2^c({\rm mean~priors})$',size=20)
    #ax.semilogy()
    ax.set_ylim(0.2,1.8)
    #ax.set_xlim()
    ax.semilogx()
    py.tight_layout()
    py.savefig('gallery/%s-F2C.pdf'%fname.split('/')[-1])
    py.close()

  def plot_gl(self,fname):  
    D=load('data/%s'%fname)
    for i in range(D['gl'][0].size):
      if D['gl'][0][i]==0: I=i-1;break
    I=len(D['gl'][0])
    X=D['X']
    norm=D['gl'][0] 

    ax=py.subplot(111)
    Lhess,=ax.plot(X[:I],(D['gl'][0]/norm)[:I],'r--')
    DO=(D['gl'][0]-D['dgl(hess)'])/norm
    UP=(D['gl'][0]+D['dgl(hess)'])/norm
    Bhess=ax.fill_between(X[:I],DO[:I],UP[:I],alpha=0.3,color='y',zorder=1)

    Lpri,=ax.plot(X[:I],(D['gl(pri)']/norm)[:I],'k:')
    DO=(D['gl(pri)']-D['dgl(pri)'])/norm
    UP=(D['gl(pri)']+D['dgl(pri)'])/norm
    Bpri=ax.fill_between(X[:I],DO[:I],UP[:I],alpha=0.3,\
      facecolor='none',edgecolor='k',hatch='...',zorder=10)

    Lpos,=ax.plot(X[:I],(D['gl(pos)']/norm)[:I],'b--')
    DO=(D['gl(pos)']-D['dgl(pos)'])/norm
    UP=(D['gl(pos)']+D['dgl(pos)'])/norm
    Bpos=ax.fill_between(X[:I],DO[:I],UP[:I],alpha=0.3,color='b',zorder=10)
  
    for k in range(30):
      Y=D['MCgl'][k]/norm
      Lmc,=ax.plot(X[:I],Y[:I],'r-',alpha=0.3,zorder=0)

    L=[tex('Hess'),tex('pri'),tex('pri(MC)'),tex('pos')]
    H=[(Bhess,Lhess),(Bpri,Lpri),Lmc,(Bpos,Lpos)]
    ax.legend(H,L,loc=2,frameon=0,fontsize=15,ncol=2)

    ax.set_xlabel(r'$x$',size=20)
    ax.set_ylabel(r'$g/g({\rm mean~priors})$',size=20)
    #ax.semilogy()
    ax.set_ylim(0.2,1.8)
    #ax.set_xlim()
    ax.semilogx()
    py.tight_layout()
    py.savefig('gallery/%s-gl.pdf'%fname.split('/')[-1])
    py.close()


