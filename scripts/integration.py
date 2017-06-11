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
yo = 0.6

vertices = [\
(-1.00000, 0.00000),
(-0.96593, 0.00000),
(-0.334832, -0.309744),
(-0.258819, 0.00000),
(-0.70711, 0.70711),
(-0.00000, 1.00000),
(0.18023, 0.966971),
(0.172941, 0.984932),
(0.212737, 0.329584),
(-0.922941, -0.384932),
(-0.827296, -0.397008)
]

# Integration limits
x1, x2 = (-0.70711, -0.334832)

# Set up the plot
fig, ax = pl.subplots(1, figsize = (8, 8))

# Plot the occulted body
x = np.linspace(-r, r, 1000)
y = np.sqrt(r ** 2 - x ** 2)
ax.plot(x, y, color = 'k', zorder = 98, lw = 1)
ax.plot(x, -y, color = 'k', zorder = 98, lw = 1)

# Integration bounds
inds = np.where((x >= x1) & (x <= x2))
ax.plot(x[inds], y[inds], 'b-', zorder = 999, lw = 2)
ax.plot(x[inds], -y[inds], '-', zorder = 999, lw = 2, color = '#aaaaee')

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

  # Integration bounds
  inds = np.where((x >= x1) & (x <= x2))
  ax.plot(x[inds], y[inds], 'b-', zorder = 999, lw = 2)
  ax.plot(x[inds], -y[inds], '-', zorder = 999, lw = 2, color = '#aaaaee')
  
  # Fill
  color = cmap(0.3 + 0.3 * (np.cos(lat) + 1))
  ax.fill_between(x, -y, y, color = color, zorder = int(-100 * lat))
  x = np.linspace(-r, xE - xlimb, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.fill_between(x, -y, y, color = color, zorder = int(-100 * lat))

# Plot the vertices that are inside the occultor
for x, y in vertices:
  ax.plot(x, y, 'ro', zorder = 1000, ms = 5)

# Integration bounds
for x in (x1, x2):
  ax.plot([x, x], [yo + np.sqrt(ro ** 2 - (x - xo) ** 2), -np.sqrt(r ** 2 - x ** 2)], '-', lw = 2, color = '#aaaaff')
  ax.plot([x, x], [yo - np.sqrt(ro ** 2 - (x - xo) ** 2), np.sqrt(r ** 2 - x ** 2)], 'b-', lw = 2)
  
# Plot the occultor
x = np.linspace(xo - ro, xo + ro, 1000)
y = np.sqrt(ro ** 2 - (x - xo) ** 2)
ax.fill_between(x, yo - y, yo + y, 
                color = 'lightgray', zorder = 99, lw = 1,
                alpha = 0.1)
ax.plot(x, yo - y, 'k-')
ax.plot(x, yo + y, 'k-')

# Integration bounds
inds = np.where((x >= x1) & (x <= x2))
ax.plot(x[inds], yo - y[inds], 'b-', zorder = 999, lw = 2)
ax.plot(x[inds], yo + y[inds], '-', zorder = 999, lw = 2,  color = '#aaaaee')

# Label regions
ax.annotate(r'$A_1$', xy = (-0.54, 0), xycoords = 'data', ha = 'center', va = 'center', color = 'b', fontsize = 14)
ax.annotate(r'$A_2$', xy = (-0.54, 0.64), xycoords = 'data', ha = 'center', va = 'center', color = 'b', fontsize = 14)
ax.annotate(r'$A_3$', xy = (-0.54, 0.85), xycoords = 'data', ha = 'center', va = 'center', color = 'b', fontsize = 14,
            xytext = (0, 40), textcoords = 'offset points', arrowprops = dict(arrowstyle = '-|>', color = 'b'))
ax.annotate(r'$v_n$', xy = (-0.70711, 0.70711), xycoords = 'data', ha = 'center', va = 'center', color = 'r', 
            xytext = (-14, 8), textcoords = 'offset points', fontsize = 14)
ax.annotate(r'$v_{n+1}$', xy = (-0.334832, -0.309744), xycoords = 'data', ha = 'center', va = 'center', color = 'r', 
            xytext = (20, -8), textcoords = 'offset points', fontsize = 14)

ax.annotate(r'$Occultor$', xy = (-1.6, 1.5), xycoords = 'data', ha = 'center', va = 'center', color = 'k', fontsize = 18)
ax.annotate(r'$Occulted$', xy = (0, -1.2), xycoords = 'data', ha = 'center', va = 'center', color = 'k', fontsize = 18)

# Appearance
ax.set_aspect('equal')
ax.set_frame_on(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_ylim(-1.2, 1.7)

# Save!
fig.savefig('../img/integration.pdf', bbox_inches = 'tight')