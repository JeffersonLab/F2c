#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import fF2C

class F2C:

  def __init__(self,root='./'):

    fF2C.root.root=root.ljust(255)

  def set_pdfset(self,ipdf):
    fF2C.setipdf(ipdf)

  def get_F2C(self,x,Q2,isf):
    """
    isf == 0: F2H
    isf == 1: F2H
    """
    return fF2C.sfh(x,Q2,isf)


if __name__=='__main__':

  f2c=F2C('./')
  f2c.set_pdfset(1)
  x,Q2,isf=0.01,5.,1
  print f2c.get_F2C(x,Q2,isf)
  



