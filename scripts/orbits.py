#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
orbits.py
---------

The observer's view of the planet orbits, highlighting
how they overlap on the sky.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import matplotlib.pyplot as pl

AUREARTH = 23454.9271
semis = np.array([11.11, 15.21, 21.44, 28.17, 37.1, 45.1, 60]) * 1e-3
radii = np.array([1.086, 1.056, 0.772, 0.918, 1.045, 1.127, 0.755]) / AUREARTH
incs = np.array([89.65, 89.67, 89.75, 89.86, 89.680, 89.710, 89.80]) * np.pi / 180
colors = ['#92C6FF', '#97F0AA', '#FF9F9A', '#D0BBFF', '#FFFEA3', '#B0E0E6', '#999999']

fig = pl.figure(figsize = (8,5))
fig.subplots_adjust(left = 0.2)

for i in range(7):
  x = np.linspace(-semis[i], semis[i], 1000)
  a = semis[i]
  b = semis[i] * np.cos(incs[i])
  y = (b / a) * np.sqrt(a ** 2 - x ** 2)
  
  pl.plot(x, y, lw = 1, color = colors[i])
  pl.plot(x, -y, lw = 1, color = colors[i])
  
  pl.fill_between(x, y - radii[i], y + radii[i], color = colors[i], alpha = 0.3)
  pl.fill_between(x, -y - radii[i], -y + radii[i], color = colors[i], alpha = 0.3)

pl.xlabel('x [AU]', fontsize = 16, fontweight = 'bold')
pl.ylabel('y [AU]', fontsize = 16, fontweight = 'bold')
pl.show()