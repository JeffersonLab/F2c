#!/usr/bin/env python
import sys,os
import numpy as np
import time
import fnmatch
import cPickle 
import zlib
from operator import mul

def checkdir(path):
  if not os.path.exists(path): 
    os.makedirs(path)

def tex(x):
  return r'$\mathrm{'+x+'}$'

def save(data,name):  
  compressed=zlib.compress(cPickle.dumps(data))
  f=open(name,"wb")
  f.writelines(compressed)
  f.close()

def load(name): 
  compressed=open(name,"rb").read()
  data=cPickle.loads(zlib.decompress(compressed))
  return data

class BAR(object):

  def __init__(self,msg,size):
    self.msg=msg
    self.size=size
    self.cnt=0

  def next(self):
    sys.stdout.write('\r')
    percentage=int(self.cnt/float(self.size)*100)
    sys.stdout.write('%s [%d%%]' % (self.msg,percentage))
    sys.stdout.flush()
    self.cnt+=1

  def finish(self):
    self.next()
    print 

def load_config(fname):
  L=open(fname).readlines()
  for l in L: exec l
  return conf   

