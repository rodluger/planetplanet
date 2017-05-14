#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
example.py
----------


'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet.eyeball import Planet, Occultor, Circle
import matplotlib.pyplot as pl
cmap = pl.get_cmap('RdBu_r')
import numpy as np
np.seterr(divide = 'ignore', invalid = 'ignore')
TOL = 1e-10

def style(lat):
  '''
  
  '''
  
  coslat = np.cos(lat)
  if np.abs(coslat) < TOL:
    return dict(color = 'k', ls = '--', lw = 1)
  else:
    return dict(color = cmap(0.5 * (coslat + 1)), ls = '-', lw = 1)

# Set up
theta = np.pi / 8
occultor = Occultor(0.5)
planet = Planet(0., -0.25, 1., theta, occultor, n = 31, noon = 0.3, midnight = 0.3)
x0 = -2
v = 4

# Plotting arrays
time = np.linspace(0, 1, 100)
flux = np.zeros_like(time)
xarr = np.linspace(-1, 1, 1000)
dummy = Circle(0, 0, 1)
ylower = dummy.val_lower(xarr)
yupper = dummy.val_upper(xarr)

# Set up the plot
fig = pl.figure(figsize = (12,6))
ax = pl.axes([0.125, 0.1, 0.825, 0.7])
ax.margins(0, None)
ax.set_xlabel('Time', fontsize = 16)
ax.set_ylabel('Flux', fontsize = 16)
aximg = pl.axes([0.125, 0.825, 0.825, 0.07])
aximg.set_xlim(0, 50)
aximg.axis('off')

# Plot light curves for different contrasts
for contrast, ls in zip([0, 0.5, 1], ['-', '--', '-.']):
  
  planet.midnight = 0.3 * contrast
  
  # Compute the light curve
  for i in range(100):
  
    # Get the flux
    planet.x0 = x0 + v * time[i]
    planet.compute()
    flux[i] = 1 + planet.delta_flux

  # Plot the light curve
  ax.plot(time, flux, color = 'k', ls = ls, label = contrast)

# Plot the planet images
for i in [5, 15, 25, 35, 45, 55, 65, 75, 85, 95]:
  planet.x0 = x0 + v * time[i]
  planet.compute()
  dx = 0.5 * i
  ellipses = planet.ellipses
  vertices = planet.vertices
  aximg.plot(dx + xarr, ylower, color = 'k', zorder = 98, lw = 1)
  aximg.plot(dx + xarr, yupper, color = 'k', zorder = 98, lw = 1)
  aximg.plot(dx + occultor.r * xarr - planet.x0, occultor.r * ylower - planet.y0, color = 'k', zorder = 99, lw = 1)
  aximg.plot(dx + occultor.r * xarr - planet.x0, occultor.r * yupper - planet.y0, color = 'k', zorder = 99, lw = 1)
  aximg.fill_between(dx + occultor.r * xarr - planet.x0, occultor.r * ylower - planet.y0, occultor.r * yupper - planet.y0, color = 'lightgray', zorder = 99)
  for ellipse in ellipses:
    x = np.linspace(ellipse.x0 - ellipse.b, ellipse.x0 + ellipse.b, 1000)
    if planet.theta > 0:
      x[x < ellipse.x0 - ellipse.xlimb] = np.nan
    else:
      x[x > ellipse.x0 - ellipse.xlimb] = np.nan
    aximg.plot(dx + x - planet.x0, ellipse.val_lower(x) - planet.y0, **style(ellipse.latitude))
    aximg.plot(dx + x - planet.x0, ellipse.val_upper(x) - planet.y0, **style(ellipse.latitude))

ax.legend(loc = 'lower left', title = 'Night/Day')
pl.show()

