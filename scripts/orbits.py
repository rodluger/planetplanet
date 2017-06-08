#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
orbits.py
---------

Plots the orbital path of each of the seven TRAPPIST-1 planets as seen
by an observer on Earth. The width of each path is the planet diameter.
Planet-planet occultations may occur anywhere where two orbits cross.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import matplotlib.pyplot as pl

AUREARTH = 23454.9271
rstar = 12.758 / AUREARTH
semis = np.array([11.11, 15.21, 21.44, 28.17, 37.1, 45.1, 60]) * 1e-3
radii = np.array([1.086, 1.056, 0.772, 0.918, 1.045, 1.127, 0.755])
incs = np.array([89.65, 89.67, 89.75, 89.86, 89.680, 89.710, 89.80]) * np.pi / 180
colors = ['#92C6FF', '#97F0AA', '#FF9F9A', '#D0BBFF', '#FFFEA3', '#B0E0E6', '#999999']

fig = pl.figure(figsize = (12,5))

# Plot the star
x = np.linspace(-rstar, rstar, 1000)
pl.fill_between(x, -np.sqrt(rstar ** 2 - x ** 2) * AUREARTH, np.zeros_like(x), color = 'orange')
pl.plot(x, -np.sqrt(rstar ** 2 - x ** 2) * AUREARTH, color = 'gray', lw = 0.5)
pl.fill_between(x, np.sqrt(rstar ** 2 - x ** 2) * AUREARTH, np.zeros_like(x), color = 'orange', zorder = 99)
pl.plot(x, np.sqrt(rstar ** 2 - x ** 2) * AUREARTH, color = 'gray', lw = 0.5, zorder = 99)

# Plot the planet orbits
for i in range(7):
  x = np.linspace(-semis[i], semis[i], 1000)
  a = semis[i]
  b = semis[i] * np.cos(incs[i])
  y = (b / a) * np.sqrt(a ** 2 - x ** 2) * AUREARTH
  
  pl.plot(x, y, lw = 1, color = colors[i])
  pl.plot(x, -y, lw = 1, color = colors[i])
  
  pl.fill_between(x, y - radii[i], y + radii[i], color = colors[i], alpha = 0.5)
  pl.fill_between(x, -y - radii[i], -y + radii[i], color = colors[i], alpha = 0.5)

pl.xlabel('x [AU]', fontsize = 16, fontweight = 'bold')
pl.ylabel(r'y [R$_\oplus$]', fontsize = 16, fontweight = 'bold')
pl.ylim(-0.0003 * AUREARTH, 0.0003 * AUREARTH)
pl.show()