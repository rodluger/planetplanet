#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
contrast.py
-----------

Plots an occultation event in two different limits: the airless limit
and the thick atmosphere limit. The asymmetry of the light curve in the 
former case betrays a strong day/night temperature contrast on the occulted
planet.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo.ppo import Planet, Star, System
from planetplanet.constants import *
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np
np.random.seed(1234)
  
# The figure for the paper
fig, ax = pl.subplots(4, figsize = (8, 18)) #, sharex = True)
fig.subplots_adjust(top = 0.8)
  
# Plot both an airless and a limb-darkened planet
for color, airless, label, dt, df in zip(['b', 'g'], [False, True], ['Thick atmosphere', 'Airless'], [0, -0.0015], [1, 1.1728478216578166]):
  
  # Instantiate the star
  mstar = 0.0802
  rstar = 0.121
  teff = (0.000524 * LSUN / (4 * np.pi * (rstar * RSUN) ** 2 * SBOLTZ)) ** 0.25
  star = Star('A', m = mstar, r = rstar, teff = teff, color = 'k')
  
  # Instantiate `c`
  RpRs = np.sqrt(0.687 / 100)
  r = RpRs * rstar * RSUN / REARTH
  c = Planet('c', m = 1.38, per = 2.4218233, inc = 89.67 - 0.05, r = r, t0 = 0, 
             Omega = 0, w = 0, ecc = 0, color = 'coral', tnight = 40., albedo = 0.3, 
             airless = True, phasecurve = False, nz = 31)

  # Instantiate `d`
  RpRs = np.sqrt(0.367 / 100)
  r = RpRs * rstar * RSUN / REARTH    
  d = Planet('d', m = 0.41, per = 4.049610, inc = 89.75 + 0.16, r = r, t0 = 0, 
             Omega = 0, w = 0, ecc = 0, color = 'firebrick', tnight = 40., albedo = 0.3, 
             airless = True, phasecurve = False)

  # Instantiate the system
  system = System(star, c, d, distance = 12, oversample = 1)

  # There's an occultation of `c` at this time
  time = np.arange(-259.684 + 2 * 0.00025, -259.665, 0.01 * MINUTE)
  minutes = (time - np.nanmedian(time)) / MINUTE
  
  # System
  system = System(star, c, d)
  system.c.airless = airless
  system.compute(time, lambda2 = 15)
  flux = np.array(c.flux[:,-1])

  # Stellar baseline
  norm = np.nanmedian(star.flux[:,-1])
  tmid = len(time) // 2

  # Plot the light curves
  ax[0].plot(minutes, (norm + flux) / norm, '-', color = color, label = label)
  
  # Compute the shifted, normalized airless light curve
  if airless:
    fairless1 = (norm + flux) / norm
    fairless2 = np.interp(time, time + dt, 1 - df * (1 - (norm + flux) / norm))
  else:
    fthickatm = (norm + flux) / norm

# Plot the residuals
ax[1].plot(minutes, 1e6 * (fairless1 - fthickatm), '-', color = 'k', label = 'Residuals')

# Plot the shifted, normalized airless light curve and the thick atmosphere light curve for comparison
ax[2].plot(minutes, fthickatm, '-', color = 'b', label = 'Thick atmosphere')
ax[2].plot(minutes, fairless2, '-', color = 'g', label = 'Airless')
  
# Plot the residuals
ax[3].plot(minutes, 1e6 * (fairless2 - fthickatm), '-', color = 'k', label = 'Residuals')

# Plot the images
axim = [None for i in range(5)]
for i, t in enumerate([1333 - 600, 1333 - 300, 1333, 1333 + 300, 1333 + 600]):
  axim[i] = fig.add_axes([0.16 + 0.15 * i, 0.8, 0.1, 0.15])
  axim[i].set_aspect('equal')
  axim[i].set_xlim(-2.6, 2.6)
  axim[i].axis('off')
  system.plot_image(t, c, ax = axim[i], occultors = [2])

# Arrows
axim[0].annotate("", xy = (0, -1.2), xytext = (83, -79), textcoords = "offset points", clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
axim[1].annotate("", xy = (0, -1.2), xytext = (44, -79), textcoords = "offset points", clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
axim[2].annotate("", xy = (0, -1.2), xytext = (0, -79), textcoords = "offset points", clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
axim[3].annotate("", xy = (0, -1.2), xytext = (-44, -79), textcoords = "offset points", clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
axim[4].annotate("", xy = (0, -1.2), xytext = (-83, -79), textcoords = "offset points", clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

# Appeareance
ax[0].legend(loc = 'lower left', fontsize = 10, frameon = False)
for axis in ax:
  axis.get_yaxis().set_major_locator(MaxNLocator(4))
  axis.get_xaxis().set_major_locator(MaxNLocator(8))
  for tick in axis.get_xticklabels() + axis.get_yticklabels():
    tick.set_fontsize(12)
  axis.ticklabel_format(useOffset = False)
ax[0].set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 16)
ax[1].set_ylabel(r'Residuals [ppm]', fontweight = 'bold', fontsize = 16, labelpad = 29)
ax[2].set_ylabel(r'Shifted Flux', fontweight = 'bold', fontsize = 16)
ax[3].set_ylabel(r'Residuals [ppm]', fontweight = 'bold', fontsize = 16, labelpad = 38)
ax[0].margins(None, 0.1)
ax[1].set_ylim(-45,91)
ax[2].margins(None, 0.1)
ax[3].set_ylim(-3, 46)
ax[3].set_xlabel('Time [minutes]', fontweight = 'bold', fontsize = 16, labelpad = 15)
fig.savefig("../img/contrast.pdf", bbox_inches = 'tight')