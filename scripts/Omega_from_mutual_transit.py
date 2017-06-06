#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Omega_from_mutual_transit.py
----------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np
from tqdm import tqdm
np.random.seed(1234)
RSUN = 6.957e8
REARTH = 6.3781e6

# Draw from prior several times and plot
niter = 100
time = np.linspace(-0.05, 0.05, 1000)
dt = (time[1] - time[0]) * 1440.
Omega = np.linspace(0, 180, 100)
duration = np.zeros((niter, len(Omega)))
for n in tqdm(range(niter)):
  
  # First run is the mean model
  if n == 0:
    N = lambda mu, sigma: mu
  else:
    N = lambda mu, sigma: mu + sigma * np.random.randn()

  # Instantiate the star
  mstar = N(0.0802, 0.0073)
  rstar = N(0.117, 0.0036)
  star = Star('A', m = mstar, r = rstar, color = 'k')

  # Planet b
  per = N(1.51087081, 0.60e-6)
  m = 0
  while m <= 0:
    m = N(0.85, 0.72)
  inc = N(89.65, 0.245)
  if inc > 90:
    inc = 180 - inc
  w = 360. * np.random.rand()
  ecc = 1
  while (ecc < 0) or (ecc >= 1):
    ecc = N(0.0005, 0.0001)
  mu = np.sqrt(0.7266 / 100)
  sig = 0.5 * 0.0088 / 100 / mu
  RpRs = N(mu, sig)
  r = RpRs * rstar * RSUN / REARTH
  b = Planet('b', m = m, per = per, inc = inc, r = r, trn0 = 0., Omega = 0, w = w, ecc = ecc, color = 'r')

  # Planet c
  per = N(2.4218233, 0.17e-5)
  m = 0
  while m <= 0:
    m = N(1.38, 0.61)
  inc = N(89.67, 0.17)
  if inc > 90:
    inc = 180 - inc
  w = 360. * np.random.rand()
  ecc = 1
  while (ecc < 0) or (ecc >= 1):
    ecc = N(0.004, 0.001)
  mu = np.sqrt(0.687 / 100)
  sig = 0.5 * 0.010 / 100 / mu
  RpRs = N(mu, sig)
  r = RpRs * rstar * RSUN / REARTH
  c = Planet('c', m = m, per = per, inc = inc, r = r, trn0 = 0., Omega = 0, w = w, ecc = ecc, color = 'b')

  # The system
  system = System(star, b, c, ttvs = False, adaptive = True)

  # Compute
  for i in range(len(Omega)):
    system.bodies[2].Omega = Omega[i] * np.pi / 180
    system.compute_orbits(time)
    duration[n,i] = dt * len(np.where(system.bodies[1].occultor == 4)[0])
    
# Plot duration
fig = pl.figure()
for n in range(niter):
  pl.plot(duration[n], Omega, alpha = 0.1, color = 'k')
pl.plot(duration[0], Omega, 'r')
pl.xlabel('Duration of mutual transit [minutes]', fontsize = 12, fontweight = 'bold')
pl.ylabel(r'$\mathbf{\Delta\Omega}$ [degrees]', fontsize = 14, fontweight = 'bold')
pl.ylim(0, 180)
pl.xlim(0.1, 50)

pl.show()