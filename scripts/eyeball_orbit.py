#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball_orbit.py
----------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import eyeball, Star, Planet, System
import matplotlib.pyplot as pl
import numpy as np

# Orbital params
inc = 60.
Omega = 0.
ecc = 0.5
w = 0

# Compute two cases: one with no hotspot offset, one with an offset
fig = [None, None, None, None]
ax = [None, None, None, None]
phase = [None, None]
flux = [None, None]

for i, dpsi, dlambda in zip([0, 1], [0, 60], [0, 30]):

  # Plot the geometry
  fig[i], ax[i] = eyeball.DrawOrbit(inc = inc, Omega = Omega, ecc = ecc, w = w, size = 1.75, dpsi = dpsi, dlambda = dlambda,
                                    plot_phasecurve = False, label_phases = True, rasterize = True)

  # Compute the phase curve
  star = Star('A')
  b = Planet('b', per = 10., inc = inc, Omega = Omega, t0 = 0, ecc = ecc, w = w, 
             dlambda = dlambda, dpsi = dpsi, airless = True, phasecurve = True)
  system = System(star, b, mintheta = 0.001)
  time = np.linspace(-5, 5, 1000)
  system.compute(time)
  phase[i] = np.linspace(0, 1, len(b.time))
  flux[i] = np.array(b.flux[:,0])

# Plot the first phasecurve
fig[2], ax[2] = pl.subplots(1, figsize = (8, 2))
fig[2].subplots_adjust(bottom = 0.3)
ax[2].plot(phase[0], flux[0] / np.nanmax(flux[0]), 'k-')

# Plot both phasecurves
fig[3], ax[3] = pl.subplots(1, figsize = (8, 2))
fig[3].subplots_adjust(bottom = 0.3)
ax[3].plot(phase[0], flux[0] / np.nanmax(flux[1]), 'k-', alpha = 0.3)
ax[3].plot(phase[1], flux[1] / np.nanmax(flux[1]), 'k-')

# Adjust appearance
for axis in [ax[2], ax[3]]:
  axis.set_xlabel('Orbital phase', fontweight = 'bold', fontsize = 12)
  axis.set_ylabel('Relative flux', fontweight = 'bold', fontsize = 12)
  axis.spines['top'].set_visible(False)
  axis.spines['right'].set_visible(False)
  axis.margins(0.01, 0.2)

fig[0].savefig('../img/eyeball_orbit1.pdf', bbox_inches = 'tight', dpi = 600)
fig[1].savefig('../img/eyeball_orbit2.pdf', bbox_inches = 'tight', dpi = 600)
fig[2].savefig('../img/eyeball_phasecurve1.pdf', bbox_inches = 'tight')
fig[3].savefig('../img/eyeball_phasecurve2.pdf', bbox_inches = 'tight')

# Show
pl.show()