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
b = Planet('b', per = 10., inc = 80., Omega = 0., t0 = 0)
system = System(star, b)
time = np.linspace(-5, 5, 1000)
system.compute_orbits(time)

# Get the eyeball angles
theta = np.arctan2(b.z, b.x)
dlambda = -np.arctan2(b.y, np.abs(b.x))

# Plot the orbit
fig, ax = pl.subplots(1, figsize = (8,8))
fig.subplots_adjust(left = 0, top = 1, bottom = 0, right = 1)
ax.plot(b.x, b.y, 'k-', alpha = 0.5)

# Adjust the plot range
xmin = min(b.x.min(), b.y.min())
xmax = max(b.x.max(), b.x.max())
dx = xmax - xmin
xmin -= 0.1 * dx
xmax += 0.1 * dx
ax.set_xlim(xmin, xmax)
ax.set_ylim(xmin, xmax)
x = (b.x - xmin) / (xmax - xmin)
y = (b.y - xmin) / (xmax - xmin)

# Plot images at different phases
for i in range(0, 1000, 50):
  ax.plot([0, b.x[i]], [0, b.y[i]], 'k-', alpha = 0.5, lw = 1)
  eyeball.Draw(theta = theta[i], dlambda = dlambda[i], fig = fig, pos = [x[i] - 0.015, y[i] - 0.015, 0.03, 0.03])

# Show
pl.show()
