#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
orbits.py |github|
------------------

Plots the orbital path of each of the seven TRAPPIST-1 planets as seen
by an observer on Earth. The width of each path is the planet diameter.
Planet-planet occultations may occur anywhere where two orbits cross.

  .. plot::
     :align: center
     
     from scripts import orbits
     orbits._test()

  .. role:: raw-html(raw)
     :format: html
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/orbits.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import matplotlib.pyplot as pl

def _test():
  '''
  
  '''
  
  plot()
  pl.show()

def plot():
  '''
  
  '''
  
  # Setup
  AUREARTH = 23454.9271
  rstar = 12.758 / AUREARTH
  semis = np.array([11.11, 15.21, 21.44, 28.17, 37.1, 45.1, 59.3]) * 1e-3
  radii = np.array([1.086, 1.056, 0.772, 0.918, 1.045, 1.127, 0.755])
  incs = np.array([89.65, 89.67, 89.75, 89.86, 89.680, 89.710, 89.80]) * np.pi / 180
  colors = ['firebrick', 'coral', 'gold', 'mediumseagreen', 'turquoise', 'cornflowerblue', 'midnightblue']
  labels = ['b', 'c', 'd', 'e', 'f' , 'g', 'h']
  fig = pl.figure(figsize = (12,5))

  # Plot the star
  x = np.linspace(-rstar, rstar, 1000)
  pl.fill_between(x, -np.sqrt(rstar ** 2 - x ** 2) * AUREARTH, np.zeros_like(x), color = 'sandybrown')
  pl.plot(x, -np.sqrt(rstar ** 2 - x ** 2) * AUREARTH, color = 'k', lw = 0.5)
  pl.fill_between(x, np.sqrt(rstar ** 2 - x ** 2) * AUREARTH, np.zeros_like(x), color = 'sandybrown', zorder = 99)
  pl.plot(x, np.sqrt(rstar ** 2 - x ** 2) * AUREARTH, color = 'k', lw = 0.5, zorder = 99)

  # Plot the planet orbits
  for i in range(7):

    # Contours
    x = np.linspace(-semis[i], semis[i], 1000)
    a = semis[i]
    b = semis[i] * np.cos(incs[i])
    y = (b / a) * np.sqrt(a ** 2 - x ** 2) * AUREARTH
    pl.plot(x, y, lw = 1, color = colors[i])
    pl.plot(x, -y, lw = 1, color = colors[i])
  
    # Dummy line for legend
    pl.plot(x, y + 999, lw = 2, color = colors[i], label = labels[i])
  
    # Fill
    pl.fill_between(x, y - radii[i], y + radii[i], color = colors[i], alpha = 0.5)
    pl.fill_between(x, -y - radii[i], -y + radii[i], color = colors[i], alpha = 0.5)

  pl.xlabel('x [AU]', fontsize = 16, fontweight = 'bold')
  pl.ylabel(r'y [R$_\oplus$]', fontsize = 16, fontweight = 'bold')
  pl.ylim(-14, 14)
  pl.legend(ncol = 2, loc = 'lower left', frameon = False)
  
  return fig, pl.gca()
  
if __name__ == '__main__':
  plot()
  fig.savefig('../pdf/orbits.pdf', bbox_inches = 'tight')