#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
integration.py
--------------

Not yet ready!

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np
cmap = pl.get_cmap('Greys_r')

# Params
nl = 5
r = 1
theta = np.pi / 4
ro = 1.
xo = -0.75
yo = 0.85

vertices = [\
(-1.00000, 0.00000),
(-0.96593, 0.00000),
(-0.25882, 0.00000),
(-0.70711, 0.70711),
(-0.00000, 1.00000),
(0.23075, 0.94525),
(0.22449, 0.97447),
(0.16088, 0.453621)
]

lines = [-0.25882, -0.70711]

# Set up the plot
fig, ax = pl.subplots(1)

# Plot the occulted body
x = np.linspace(-r, r, 1000)
y = np.sqrt(r ** 2 - x ** 2)
ax.plot(x, y, color = 'k', zorder = 98, lw = 1)
ax.plot(x, -y, color = 'k', zorder = 98, lw = 1)

# Plot the latitude ellipses
for lat in np.linspace(0, np.pi, nl + 2)[1:-1]:

  # The ellipse
  a = r * np.abs(np.sin(lat))
  b = a * np.abs(np.sin(theta))
  xE = -r * np.cos(lat) * np.cos(theta)
  yE = 0
  xlimb = r * np.cos(lat) * np.sin(theta) * np.tan(theta)
  if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
    xmin = xE - b
  else:
    xmin = xE - xlimb
  if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
    xmax = xE + b
  else:
    xmax = xE - xlimb
        
  # Plot it
  x = np.linspace(xE - b, xE + b, 1000)
  if theta > 0:
    x[x < xE - xlimb] = np.nan
  else:
    x[x > xE - xlimb] = np.nan
  A = b ** 2 - (x - xE) ** 2
  A[A < 0] = 0
  y = (a / b) * np.sqrt(A)
  if np.abs(np.cos(lat)) < 1e-5:
    style = dict(color = 'k', ls = '--', lw = 1)
  else:
    style = dict(color = 'k', ls = '-', lw = 1)
  ax.plot(x, y, **style)
  ax.plot(x, -y, **style)
  color = cmap(0.3 + 0.3 * (np.cos(lat) + 1))
  ax.fill_between(x, -y, y, color = color, zorder = int(-100 * lat))
  x = np.linspace(-r, xE - xlimb, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.fill_between(x, -y, y, color = color, zorder = int(-100 * lat))

# Plot the vertices that are inside the occultor
for x, y in vertices:
  pl.plot(x, y, 'ro', zorder = 300, ms = 3)

for x in lines:
  pl.plot([x, x], [yo - np.sqrt(ro ** 2 - (x - xo) ** 2), np.sqrt(r ** 2 - x ** 2)], 'r-')
  
# Plot the occultor
x = np.linspace(xo - ro, xo + ro, 1000)
y = np.sqrt(ro ** 2 - (x - xo) ** 2)
ax.fill_between(x, yo - y, yo + y, 
                color = 'lightgray', zorder = 99, lw = 1,
                alpha = 0.1)
ax.plot(x, yo - y, 'k-')
ax.plot(x, yo + y, 'k-')

ax.set_aspect('equal')
pl.show()