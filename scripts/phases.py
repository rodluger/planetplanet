#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
phases.py
---------

*Still working on this.*

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import eyeball, Star, Planet, System
import numpy as np
import matplotlib.pyplot as pl

# Get orbital elements
star = Star('A')
b = Planet('b', per = 10., inc = 70., Omega = 30., t0 = 0, ecc = 0.5, w = 0.)
system = System(star, b)
time = np.linspace(-5, 5, 1000)
system.compute_orbits(time)

# Plot the orbit
fig, ax = pl.subplots(1, figsize = (8,8))
fig.subplots_adjust(left = 0, top = 1, bottom = 0, right = 1)
ax.plot(b.x, b.y, 'k-', alpha = 0.5)

# Adjust the plot range
xmin = min(b.x.min(), b.y.min())
xmax = max(b.x.max(), b.y.max())
dx = xmax - xmin
xmin -= 0.1 * dx
xmax += 0.1 * dx
ax.set_xlim(xmin, xmax)
ax.set_ylim(xmin, xmax)
x = (b.x - xmin) / (xmax - xmin)
y = (b.y - xmin) / (xmax - xmin)

# Get the indices of the images we'll plot, sorted by zorder
inds = np.array(list(range(0, 1000, 50)), dtype = int)
inds = inds[np.argsort([-b.z[i] for i in inds])]

# Plot images at different phases
for i in inds:
  
  # Get the xy coordinates of the substellar point
  r = np.sqrt(b.x[i] ** 2 + b.y[i] ** 2 + b.z[i] ** 2)
  xss = -b.x[i] * (b._r / r)
  yss = -b.y[i] * (b._r / r)
  
  # Find the rotation angle that would put it on the x axis
  rotation = -np.arctan(yss / xss)

  # Find the x coordinate in the rotated frame
  xss_r = xss * np.cos(rotation) - yss * np.sin(rotation)

  # Find the y coordinate in a frame where the orbital plane
  # runs parallel to the x axis
  yprime = b.y[i] * np.cos(-b._Omega) + b.x[i] * np.sin(-b._Omega)
  
  # This is the effective orbital phase
  if yprime > 0:
    theta = np.pi - np.arccos(xss_r / b._r)
  else:
    theta = np.pi + np.arccos(xss_r / b._r)
  
  
  ax.plot([0, b.x[i]], [0, b.y[i]], 'k-', alpha = 0.5, lw = 1)
  eyeball.Draw(theta = theta, rotation = rotation, fig = fig, pos = [x[i] - 0.015, y[i] - 0.015, 0.03, 0.03])

# Show
pl.show()
