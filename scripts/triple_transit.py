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
from planetplanet.photo import Planet, Star, System
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
b = Planet('b', m = 1, per = 2, inc = 90.4, r = 2., trn0 = 0, 
           nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# Planet c
c = Planet('c', m = 1, per = 8, inc = 90., r = 2., trn0 = 0.0005, 
           nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# Planet d
d = Planet('d', m = 1, per = 32, inc = 89.94, r = 2., trn0 = 0.002, 
           nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# System
system = System(star, b, c, d)

# Get the occultation light curves
time = np.linspace(-0.06, 0.06, 1000)
system.compute(time)

# Set up the figure
fig = pl.figure(figsize = (8, 5))
fig.subplots_adjust(left = 0.175)

# Plot three different wavelengths (first, mid, and last)
axlc = pl.subplot2grid((5, 5), (1, 0), colspan = 5, rowspan = 4)
axlc.plot(star.time * 1440, star.flux[:, 0] / star.flux[0, 0], 'b-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * star.wavelength[0])) + r"\ \mu\mathrm{m}$")
axlc.plot(star.time * 1440, star.flux[:, star.flux.shape[-1] // 2] / star.flux[0, star.flux.shape[-1] // 2], 'g-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * star.wavelength[star.flux.shape[-1] // 2])) + r"\ \mu\mathrm{m}$")
axlc.plot(star.time * 1440, star.flux[:, -1] / star.flux[0, -1], 'r-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * star.wavelength[-1])) + r"\ \mu\mathrm{m}$")
axlc.set_xlabel('Time [minutes]', fontweight = 'bold', fontsize = 10)
axlc.set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 10)
axlc.get_yaxis().set_major_locator(MaxNLocator(4))
axlc.get_xaxis().set_major_locator(MaxNLocator(8))
axlc.margins(0, None)
for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
  tick.set_fontsize(8)
axlc.legend(loc = 'lower right', fontsize = 8)

# Plot the images
axim = [pl.subplot2grid((5, 5), (0, n), colspan = 1, rowspan = 1) for n in range(5)]
t = [300, 400, 500, 600, 700]
for n in range(5):
  system.plot_image(t[n], star, occultors = [1,2,3], ax = axim[n])
  axim[n].axis('off')
  axim[n].set_aspect('equal')
  axim[n].annotate('%d' % (n + 1), xy = (0, 12), xycoords = 'data', ha = 'center', va = 'bottom', fontweight = 'bold', fontsize = 8, clip_on = False)
  axlc.annotate('%d' % (n + 1), xy = (star.time[t[n]] * 1440, 0.01 + star.flux[t[n], -1] / star.flux[0, -1]), xycoords = 'data', ha = 'center', fontweight = 'bold', fontsize = 8)
  
axim[0].set_xlim(-30,15)
axim[1].set_xlim(-26.25,18.75)
axim[2].set_xlim(-22.5,22.5)
axim[3].set_xlim(-18.75,26.25)
axim[4].set_xlim(-15,30)

fig.savefig('../img/triple.pdf', bbox_inches = 'tight')
pl.close()

# Animate!
fig, axlc, axxz, axim = system.plot_occultation('A', 0.) #, gifname = 'triple')
pl.show()