#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
timing.py
---------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.constants import *
from planetplanet.photo import Star, Planet, System
import matplotlib.pyplot as pl
import numpy as np
import batman
import pysyzygy as ps
from tqdm import tqdm
import timeit, builtins

# System params
time = np.arange(-0.2, 0.2, 0.0001)
mstar = 1.
rstar = 1.
limbdark = [0.4, 0.26]
per = 5.
inc = 90.
r = 10.
t0 = 0.
w = 60.
ecc = 0.3

def run_pp():
  '''
  
  '''
  
  # planetplanet
  star = Star('A', m = mstar, r = rstar, nz = 11, limbdark = limbdark)
  b = Planet('b', m = 0., per = per, inc = inc, r = r, t0 = t0, 
             nz = 1, Omega = 0., w = w, ecc = ecc, phasecurve = False)
  system = System(star, b, R = 1, lambda1 = 5, lambda2 = 6, quiet = True)
  system.compute(time)
  flux_pp = system.A.flux[:,0]
  flux_pp /= flux_pp[0]
  return flux_pp
  
def run_bm():
  '''
  
  '''
  
  # batman
  params = batman.TransitParams()
  params.t0 = t0
  params.per = per
  params.inc = inc
  params.ecc = ecc
  params.w = w - 180.
  params.limb_dark = "quadratic"
  params.u = limbdark
  params.rp = r / (rstar * RSUNREARTH)
  params.a = ((per) ** 2 * GEARTH * (mstar * MSUNMEARTH) / (4 * np.pi ** 2)) ** (1. / 3.) / (rstar * RSUNREARTH)
  m = batman.TransitModel(params, time)
  flux_bm = m.light_curve(params)
  return flux_bm

def run_ps():
  '''
  
  '''
  
  # pysyzygy
  a = ((per) ** 2 * GEARTH * (mstar * MSUNMEARTH) / (4 * np.pi ** 2)) ** (1. / 3.) / (rstar * RSUNREARTH)
  trn = ps.Transit(t0 = t0, per = per, RpRs = r / (rstar * RSUNREARTH),
                   ecc = ecc, w = w * np.pi / 180 + np.pi, u1 = limbdark[0], u2 = limbdark[1],
                   bcirc = a * np.cos(inc * np.pi / 180), aRs = a)                 
  flux_ps = trn(time, 'unbinned')
  return flux_ps
  
# Time it!
builtins.__dict__.update(locals())
tpp = timeit.timeit('run_pp()', number = 10) / 10.
tbm = timeit.timeit('run_bm()', number = 10) / 10.
tps = timeit.timeit('run_ps()', number = 10) / 10.

print(tpp, tbm, tps)