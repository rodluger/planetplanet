#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Omega_from_mutual_transit.py
----------------------------

For random draws from the prior, computes the duration of a
mutual transit between TRAPPIST-1 b and c as a function of the
difference in their longitude of ascending nodes. With some
scatter, the difference in this angle is inversely proportional
to the duration of the mutual transit. Observing such an event
can place strong constraints on the longitudes of ascending node
of the TRAPPIST-1 planets.


'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np
from tqdm import tqdm
import corner
np.random.seed(1234)
RSUN = 6.957e8
REARTH = 6.3781e6

# Draw from prior several times and plot
niter = 10000
time = np.linspace(-0.05, 0.05, 1000)
dt = (time[1] - time[0]) * 1440.
Omega = np.linspace(0, 180, 100)
Omega_zoom = np.linspace(0, 3, 100)
duration = np.zeros((niter, len(Omega)))
duration_zoom = np.zeros((niter, len(Omega)))
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
  b = Planet('b', m = m, per = per, inc = inc, r = r, t0 = 0., Omega = 0, w = w, ecc = ecc, color = 'r')

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
  c = Planet('c', m = m, per = per, inc = inc, r = r, t0 = 0., Omega = 0, w = w, ecc = ecc, color = 'b')

  # The system
  system = System(star, b, c, nbody = True, adaptive = True, quiet = True)

  # Compute
  for i in range(len(Omega)):
    system.bodies[2].Omega = Omega[i]
    system.compute_orbits(time)
    duration[n,i] = dt * len(np.where(system.bodies[1].occultor == 4)[0])
  
  # Compute (zoom)
  for i in range(len(Omega_zoom)):
    system.bodies[2].Omega = Omega_zoom[i]
    system.compute_orbits(time)
    duration_zoom[n,i] = dt * len(np.where(system.bodies[1].occultor == 4)[0])
  
# Plot the heatmap
fig = pl.figure(figsize = (6, 8))
ax = [pl.subplot2grid((3,1), (0,0), rowspan = 2),
      pl.subplot2grid((3,1), (2,0))]

# Top plot
duration = duration.reshape(-1)
Omega = np.concatenate([Omega for i in range(niter)])
iszero = np.where(duration == 0)
duration = np.delete(duration, iszero)
Omega = np.delete(Omega, iszero)
corner.hist2d(duration, Omega, bins = 100, ax = ax[0])

# Zoom inset
duration_zoom = duration_zoom.reshape(-1)
Omega_zoom = np.concatenate([Omega_zoom for i in range(niter)])
iszero = np.where(duration_zoom == 0)
duration_zoom = np.delete(duration_zoom, iszero)
Omega_zoom = np.delete(Omega_zoom, iszero)
corner.hist2d(duration_zoom, Omega_zoom, bins = 100, ax = ax[1])

# Appearance
ax[0].set_ylabel(r'$\mathbf{\Delta\Omega}$ [degrees]', fontsize = 14, fontweight = 'bold')
ax[1].set_ylabel(r'$\mathbf{\Delta\Omega}$ [degrees]', fontsize = 14, fontweight = 'bold')
ax[1].set_xlabel('Duration of mutual transit [minutes]', fontsize = 14, fontweight = 'bold')
ax[0].set_ylim(0, 180)
ax[1].set_ylim(0, 3)

fig.savefig("../img/Omega.pdf")
pl.show()