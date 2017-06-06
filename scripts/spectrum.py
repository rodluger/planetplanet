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
system.compute(time, lambda1 = 5, lambda2 = 15, R = 3)

# Get the normalized fluxes at each wavelength
flux5 = system.flux[:,0]
flux5 /= np.nanmedian(flux5)
flux10 = system.flux[:,1]
flux10 /= np.nanmedian(flux10)
flux15 = system.flux[:,2]
flux15 /= np.nanmedian(flux15)

# Plot
fig, ax = pl.subplots(3, figsize = (11, 7), sharex = True)
fig.subplots_adjust(bottom = 0.1, top = 0.95, left = 0.1, right = 0.95)
ax[0].plot(system.time, flux5, color = 'b', lw = 1, label = '5 microns')
ax[1].plot(system.time, flux10, color = 'g', lw = 1,  label = '10 microns')
ax[2].plot(system.time, flux15, color = 'r', lw = 1,  label = '15 microns')

# Appearance
ax[0].set_ylabel('5 $\mathbf{\mu}$m Flux', fontsize = 12, fontweight = 'bold')
ax[1].set_ylabel('10 $\mathbf{\mu}$m Flux', fontsize = 12, fontweight = 'bold')
ax[2].set_ylabel('15 $\mathbf{\mu}$m Flux', fontsize = 12, fontweight = 'bold')
ax[2].set_xlabel('Time [days]', fontsize = 14, fontweight = 'bold')
pl.setp(ax[0].get_xticklabels(), visible = False)
pl.setp(ax[1].get_xticklabels(), visible = False)
for axis in ax:
  axis.margins(0, 0.1)
  
pl.show()