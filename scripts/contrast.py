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
from planetplanet.photo.ppo import Planet, Star, System, DrawEyeball
from planetplanet.constants import *
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np
np.random.seed(1234)
  
# The figure for the paper
fig = pl.figure(figsize = (8, 14))
ax = [pl.subplot2grid((6, 1), (0, 0), colspan = 1, rowspan = 2),
      pl.subplot2grid((6, 1), (2, 0), colspan = 1, rowspan = 1),
      pl.subplot2grid((6, 1), (3, 0), colspan = 1, rowspan = 2),
      pl.subplot2grid((6, 1), (5, 0), colspan = 1, rowspan = 1)]
fig.subplots_adjust(top = 0.8)
  
# Plot both an airless and a limb-darkened planet
for color, airless, label, dt, df in zip(['b', 'g'], [False, True], ['Thick atmosphere', 'Airless'], [0, -0.00165], [1, 1.015 * 1.1728478216578166]):
  
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
  # I eyeballed the scaled depth and duration to get the best match
  # to egress. We can probably do a little better, but it's still going
  # to be a ~20 ppm signal in the residuals -- there's no way around the
  # asymmetry!
  if airless:
    fairless1 = (norm + flux) / norm
    fairless2 = np.interp(time, time + dt, 1 - df * (1 - (norm + flux) / norm))
    fairless3 = np.interp(minutes, 1.5 * minutes, fairless2)
  else:
    fthickatm = (norm + flux) / norm

# Plot the residuals
ax[1].plot(minutes, 1e6 * (fairless1 - fthickatm), '-', color = 'k', label = 'Residuals')

# Plot the shifted, normalized airless light curve and the thick atmosphere light curve for comparison
ax[2].plot(minutes, fthickatm, '-', color = 'b', label = 'Thick atmosphere')
ax[2].plot(minutes, fairless3, '-', color = 'g', label = 'Airless')
  
# Plot the residuals
ax[3].plot(minutes, 1e6 * (fairless3 - fthickatm), '-', color = 'k', label = 'Residuals')

# Plot the images
x0 = 0.5102
dx = 0.15
px = [x0 - 2 * dx, x0 - dx, x0, x0 + dx, x0 + 2 * dx]
for i, t in enumerate([1333 - 600, 1333 - 300, 1333, 1333 + 300, 1333 + 600]):
  rp = c._r
  x0 = c.x_hr[t]
  y0 = c.y_hr[t]
  z0 = c.z_hr[t]
  x = x0 * np.cos(c._Omega) + y0 * np.sin(c._Omega)
  y = y0 * np.cos(c._Omega) - x0 * np.sin(c._Omega)
  z = z0
  r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
  xprime = c._r * np.cos(c._Phi) * np.sin(c._Lambda)
  yprime = c._r * np.sin(c._Phi)
  zprime = r - c._r * np.cos(c._Phi) * np.cos(c._Lambda)
  rxz = np.sqrt(x ** 2 + z ** 2)
  xstar = ((z * r) * xprime - (x * y) * yprime + (x * rxz) * zprime) / (r * rxz)
  ystar = (rxz * yprime + y * zprime) / r
  zstar = (-(x * r) * xprime - (y * z) * yprime + (z * rxz) * zprime) / (r * rxz)
  xstar, ystar = xstar * np.cos(c._Omega) - ystar * np.sin(c._Omega), \
                 ystar * np.cos(c._Omega) + xstar * np.sin(c._Omega)
  x = x0
  y = y0
  dist = np.sqrt((xstar - x) ** 2 + (ystar - y) ** 2)
  gamma = np.arctan2(ystar - y, xstar - x) + np.pi
  if (zstar - z) <= 0:
    theta = np.arccos(dist / c._r)
  else:
    theta = -np.arccos(dist / c._r)
  occ_dict = [dict(x = d.x_hr[t] - x0, y = d.y_hr[t] - y0, r = d._r, zorder = i + 1, alpha = 1)]
  DrawEyeball(x0 = 0, y0 = 0, r = c._r, theta = theta, nz = 31, gamma = gamma, 
              draw_ellipses = False,
              occultors = occ_dict, cmap = 'inferno', fig = fig, 
              pos = [px[i], 0.85, 0.05, 0.05], rasterize = True)

# Arrows
ax[0].annotate("", xy = (minutes[1333 - 600], 1.000008), xycoords = "data", xytext = (-80, 40), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
ax[0].annotate("", xy = (minutes[1333 - 300], 1.000008), xycoords = "data", xytext = (-40, 40), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
ax[0].annotate("", xy = (minutes[1333], 1.000008), xycoords = "data", xytext = (0, 40), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
ax[0].annotate("", xy = (minutes[1333 + 300], 1.000008), xycoords = "data", xytext = (40, 40), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
ax[0].annotate("", xy = (minutes[1333 + 600], 1.000008), xycoords = "data", xytext = (80, 40), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))


# Appeareance
ax[0].legend(loc = 'lower left', fontsize = 10, frameon = False)
for i, axis in enumerate(ax):
  axis.get_yaxis().set_major_locator(MaxNLocator(4))
  axis.get_xaxis().set_major_locator(MaxNLocator(8))
  for tick in axis.get_xticklabels() + axis.get_yticklabels():
    tick.set_fontsize(12)
  axis.ticklabel_format(useOffset = False)
  if i < 3:
    axis.set_xticklabels([])
ax[0].set_ylabel(r'Flux', fontweight = 'bold', fontsize = 16)
ax[1].set_ylabel(r'$\Delta$ [ppm]', fontweight = 'bold', fontsize = 16, labelpad = 29)
ax[2].set_ylabel(r'Scaled Flux', fontweight = 'bold', fontsize = 16)
ax[3].set_ylabel(r'$\Delta$ [ppm]', fontweight = 'bold', fontsize = 16, labelpad = 29)
ax[0].margins(None, 0.1)
ax[1].set_ylim(-45,91)
ax[2].margins(None, 0.1)
ax[3].set_ylim(-15, 15)
ax[3].set_yticks([-10,0,10])
ax[-1].set_xlabel('Time [minutes]', fontweight = 'bold', fontsize = 16, labelpad = 15)
fig.savefig("../img/contrast.pdf", bbox_inches = 'tight', dpi = 600)