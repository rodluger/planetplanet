#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
spectrum.py
-----------

Computes and plots a light curve of the TRAPPIST-1 system over three days in the
wavelength range 5-15 microns, with orbital parameters drawn at random 
from their prior distributions. All transits, secondary eclipses, 
planet-planet occultations, and mutual transits are shown.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, ttvs = False, phasecurve = False, adaptive = True)

# Get the occultation light curves for the first 10 days
time = np.linspace(0., 10., 10000)
system.compute(time, lambda1 = 5, lambda2 = 15, R = 1000)

# Get the flux at the first, middle, and last wavelength
wavelength = []
flux = []
for index in [0, len(system.wavelength) // 2, -1]:
  wavelength.append(system.wavelength[index])
  flux.append(system.flux[:,index] / np.nanmedian(system.flux[:, index]))

# Plot
fig, ax = pl.subplots(3, figsize = (11, 7), sharex = True)
fig.subplots_adjust(bottom = 0.1, top = 0.95, left = 0.1, right = 0.95)
ax[0].plot(system.time, flux[0], color = 'b', lw = 1)
ax[1].plot(system.time, flux[1], color = 'g', lw = 1)
ax[2].plot(system.time, flux[2], color = 'r', lw = 1)

# Appearance
for n in range(3):
  ax[n].set_ylabel('%.1f $\mathbf{\mu}$m Flux' % wavelength[n], fontsize = 12, fontweight = 'bold')
ax[2].set_xlabel('Time [days]', fontsize = 14, fontweight = 'bold')
pl.setp(ax[0].get_xticklabels(), visible = False)
pl.setp(ax[1].get_xticklabels(), visible = False)
for axis in ax:
  axis.margins(0, 0.1)

pl.show()