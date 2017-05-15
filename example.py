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

def Style(lat):
  '''
  
  '''
  
  coslat = np.cos(lat)
  if np.abs(coslat) < TOL:
    return dict(color = 'k', ls = '--', lw = 1)
  else:
    return dict(color = cmap(0.5 * (coslat + 1)), ls = '-', lw = 1)

def Plot(y0 = -1.5, x0 = 2, vx = -4, vy = 3, ro = 0.85, theta = -np.pi / 6, n = 11):
  '''
  
  '''
  
  # The planet and the occultor
  occultor = Occultor(ro)
  planet = Planet(0., 0., 1., theta, occultor, n = n, noon = 1, midnight = 0.1)

  # Plotting arrays
  time = np.linspace(0, 1, 100)
  flux = np.zeros_like(time)
  xarr = np.linspace(-1, 1, 1000)
  dummy = Circle(0, 0, 1)
  ylower = dummy.val_lower(xarr)
  yupper = dummy.val_upper(xarr)

  # Set up the plot
  fig = pl.figure(figsize = (12,6))
  ax = pl.axes([0.125, 0.1, 0.825, 0.7], adjustable = 'box')
  ax.margins(0, None)
  ax.set_xlabel('Time', fontsize = 16)
  ax.set_ylabel('Normalized Flux', fontsize = 16)
  aximg = pl.axes([0.125, 0.825, 0.825, 0.15], sharex = ax)
  aximg.set_xlim(-0.01, 1.01)
  aximg.set_ylim(-0.05, 0.05)
  aximg.axis('off')

  # Plot light curves for different contrasts
  for contrast, ls in zip([0.0, 0.5, 1.0], ['-', '--', '-.']):
  
    planet.midnight = planet.noon * contrast
  
    # Compute the light curve
    for i in range(100):
  
      # Get the flux
      planet.x0 = x0 + vx * time[i]
      planet.y0 = y0 + vy * time[i]
      planet.compute_delta()
      flux[i] = planet.norm_flux

    # Plot the light curve
    ax.plot(time, flux, color = 'k', ls = ls, label = '%.1f' % contrast)

  # Plot the planet images
  stretch = 50
  for i in [5, 15, 25, 35, 45, 55, 65, 75, 85, 95]:
    dx = i / 100
    planet.x0 = x0 + vx * time[i]
    planet.y0 = y0 + vy * time[i]
    planet.compute_delta()
    ellipses = planet.ellipses
    vertices = planet.vertices
    aximg.plot(dx + xarr / stretch, ylower / stretch, color = 'k', zorder = 98, lw = 1)
    aximg.plot(dx + xarr / stretch, yupper / stretch, color = 'k', zorder = 98, lw = 1)
    aximg.plot(dx + (occultor.r * xarr - planet.x0) / stretch, (occultor.r * ylower - planet.y0) / stretch, color = 'k', zorder = 99, lw = 1)
    aximg.plot(dx + (occultor.r * xarr - planet.x0) / stretch, (occultor.r * yupper - planet.y0) / stretch, color = 'k', zorder = 99, lw = 1)
    aximg.fill_between(dx + (occultor.r * xarr - planet.x0) / stretch, (occultor.r * ylower - planet.y0) / stretch, (occultor.r * yupper - planet.y0) / stretch, color = 'lightgray', zorder = 99)
    for ellipse in ellipses:
      x = np.linspace(ellipse.x0 - ellipse.b, ellipse.x0 + ellipse.b, 1000)
      if planet.theta > 0:
        x[x < ellipse.x0 - ellipse.xlimb] = np.nan
      else:
        x[x > ellipse.x0 - ellipse.xlimb] = np.nan
      aximg.plot(dx + (x - planet.x0) / stretch, (ellipse.val_lower(x) - planet.y0) / stretch, **Style(ellipse.latitude))
      aximg.plot(dx + (x - planet.x0) / stretch, (ellipse.val_upper(x) - planet.y0) / stretch, **Style(ellipse.latitude))

  legend = ax.legend(loc = 'lower left', title = 'Night/Day\n Contrast')
  pl.setp(legend.get_title(), fontsize = 9, fontweight = 'bold')  
  pl.show()

Plot()