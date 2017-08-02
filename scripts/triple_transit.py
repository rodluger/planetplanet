#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
mutual_transit.py
-----------------

Computes and plots a hypothetical triple mutual transit event, where three large 
planets transit the star and occult each other simultaneously:

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System, DrawEyeball
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np
np.random.seed(1234)

def u1(lam):
  '''
  A really silly linear limb darkening law with a linear
  wavelength dependence.
  
  '''
  
  lam = np.atleast_1d(lam)
  result = 0.5 * (1 - (lam - 5) / 10) + 0.5
  result[lam < 5] = 1.
  result[lam > 15] = 0.5
  return result
  
# Instantiate the star
star = Star('Star', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [u1])

# Planet b
b = Planet('b', m = 1, per = 2, inc = 90.4, r = 2., t0 = 0, 
           nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# Planet c
c = Planet('c', m = 1, per = 8, inc = 90., r = 2., t0 = 0.0005, 
           nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# Planet d
d = Planet('d', m = 1, per = 32, inc = 89.94, r = 2., t0 = 0.002, 
           nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# System
system = System(star, b, c, d)

# Get the occultation light curves
time = np.linspace(-0.06, 0.06, 1000)
system.compute(time)

# Set up the figure
fig = pl.figure(figsize = (7, 7))
fig.subplots_adjust(left = 0.175)

# Plot three different wavelengths (first, mid, and last)
axlc = pl.subplot2grid((60, 5), (15, 0), colspan = 5, rowspan = 40)
axlc.plot(star.time * 1440, star.flux[:, 0] / star.flux[0, 0], 'b-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * star.wavelength[0])) + r"\ \mu\mathrm{m}$")
axlc.plot(star.time * 1440, star.flux[:, star.flux.shape[-1] // 2] / star.flux[0, star.flux.shape[-1] // 2], 'g-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * star.wavelength[star.flux.shape[-1] // 2])) + r"\ \mu\mathrm{m}$")
axlc.plot(star.time * 1440, star.flux[:, -1] / star.flux[0, -1], 'r-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * star.wavelength[-1])) + r"\ \mu\mathrm{m}$")
axlc.set_xlabel('Time [minutes]', fontweight = 'bold', fontsize = 14)
axlc.set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 14)
axlc.get_yaxis().set_major_locator(MaxNLocator(4))
axlc.get_xaxis().set_major_locator(MaxNLocator(8))
axlc.margins(0, None)
for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
  tick.set_fontsize(12)
axlc.legend(loc = 'lower right', fontsize = 12)

# Plot the images
t = [300, 400, 500, 600, 700]

x0 = 0.535
dx = 0.15
px = [x0 - 2 * dx, x0 - dx, x0, x0 + dx, x0 + 2 * dx]
for n in range(5):
  
  # Convert them into a list of dicts
  occ_dict = []
  for i, occultor in enumerate([b, c, d]):
    occ_dict.append(dict(x = occultor.x_hr[t[n]], y = occultor.y_hr[t[n]], r = occultor._r, zorder = i + 1, alpha = 1))
  
  # Draw the eyeball planet and the occultors
  DrawEyeball(x0 = 0, y0 = 0, r = star._r, theta = np.pi / 2, nz = 31, gamma = 0, 
              occultors = occ_dict, cmap = 'inferno', fig = fig, 
              draw_ellipses = False, pos = [px[n], 0.85, 0.085, 0.085])

# Arrows
axlc.annotate("", xy = (star.time[t[0]] * 1440, 1.006), xycoords = "data", xytext = (-80, 63), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

axlc.annotate("", xy = (star.time[t[1]] * 1440, 1.006), xycoords = "data", xytext = (-40, 63), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

axlc.annotate("", xy = (star.time[t[2]] * 1440, 1.006), xycoords = "data", xytext = (0, 63), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

axlc.annotate("", xy = (star.time[t[3]] * 1440, 1.006), xycoords = "data", xytext = (40, 63), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

axlc.annotate("", xy = (star.time[t[4]] * 1440, 1.006), xycoords = "data", xytext = (80, 63), textcoords = "offset points", 
              clip_on = False, arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
    
fig.savefig('../img/triple.pdf', bbox_inches = 'tight')
pl.show()
pl.close()

# Animate!
fig, axlc, axxz, axim = system.plot_occultation('A', 0.) #, gifname = 'triple')
pl.show()