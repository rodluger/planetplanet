#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
timing.py |github|
------------------

Timing tests for :py:obj:`planetplanet`. Here we compare the time it takes
to compute a transit light curve among :py:obj:`planetplanet`, :py:obj:`batman`,
and :py:obj:`pysyzygy`. :py:obj:`planetplanet` is currently about 10 times slower
than the other two for transits. I'm working on speeding it up, but the main source
of the speed difference is the fact that :py:obj:`planetplanet` computes 
wavelength-dependent light curves, which means the total stellar flux is not unity (as it
is in the other two codes); instead, we need to evaluate the Planck function repeatedly,
even when we set the spectral resolution `R = 1`. It's fairly straightforward to optimize
this, so stay tuned for updates!
 
  .. plot::
     :align: center
     
     from scripts import timing
     import matplotlib.pyplot as pl
     timing.plot()
     pl.show()
  
  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/timing.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet.constants import *
from planetplanet.photo import Star, Planet, System
import matplotlib.pyplot as pl
import numpy as np
import batman
import pysyzygy as ps
from tqdm import tqdm
import timeit, builtins

# System params

mstar = 1.
rstar = 1.
limbdark = [0.4, 0.26]
per = 5.
inc = 90.
r = 10.
t0 = 0.
w = 60.
ecc = 0.3

def run_pp(N = 1000):
  '''
  
  '''
  
  # planetplanet
  time = np.linspace(-0.1, 0.1, N)
  star = Star('A', m = mstar, r = rstar, nz = 11, limbdark = limbdark)
  b = Planet('b', m = 0., per = per, inc = inc, r = r, t0 = t0, 
             nz = 1, Omega = 0., w = w, ecc = ecc, phasecurve = False)
  system = System(star, b, R = 1, lambda1 = 5, lambda2 = 6, quiet = True, batmanopt = True, circleopt = True, nbody = False)
  system.compute(time)
  flux_pp = system.A.flux[:,0]
  flux_pp /= flux_pp[0]
  return flux_pp
  
def run_bm(N = 1000):
  '''
  
  '''
  
  # batman
  time = np.linspace(-0.2, 0.2, N)
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

def run_ps(N = 1000):
  '''
  
  '''
  
  # pysyzygy
  time = np.linspace(-0.2, 0.2, N)
  a = ((per) ** 2 * GEARTH * (mstar * MSUNMEARTH) / (4 * np.pi ** 2)) ** (1. / 3.) / (rstar * RSUNREARTH)
  trn = ps.Transit(t0 = t0, per = per, RpRs = r / (rstar * RSUNREARTH),
                   ecc = ecc, w = w * np.pi / 180 + np.pi, u1 = limbdark[0], u2 = limbdark[1],
                   bcirc = a * np.cos(inc * np.pi / 180), aRs = a)                 
  flux_ps = trn(time, 'unbinned')
  return flux_ps

 
def plot():
  '''
  
  '''
  
  # Register the functions
  builtins.__dict__.update(globals())
  
  # Loop over various dataset sizes
  Narr = np.logspace(0, 5, 5)
  tpp = np.zeros_like(Narr)
  tbm = np.zeros_like(Narr)
  tps = np.zeros_like(Narr)
  for i, N in enumerate(Narr):
    tpp[i] = timeit.timeit('run_pp(%d)' % N, number = 10) / 10.
    tbm[i] = timeit.timeit('run_bm(%d)' % N, number = 10) / 10.
    tps[i] = timeit.timeit('run_ps(%d)' % N, number = 10) / 10.
  
  pl.plot(Narr, tpp, '-o', label = 'planetplanet')
  pl.plot(Narr, tbm, '-o', label = 'batman')
  pl.plot(Narr, tps, '-o', label = 'pysyzygy')
  pl.legend()
  pl.yscale('log')
  pl.xscale('log')
  pl.ylabel('Time [seconds]', fontweight = 'bold')
  pl.xlabel('Number of datapoints', fontweight = 'bold')

if __name__ == '__main__':
  plot()
  pl.show()