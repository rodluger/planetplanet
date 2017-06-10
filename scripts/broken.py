#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
broken.py
---------

Numerical issues lead to NANs in the integral
computation. Here's an example. Need to fix this.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Star, Planet, System
import matplotlib.pyplot as pl
import numpy as np

vertices = [[-0.5321789172328408,
1.0000000000000000,
0.2339105413836380,
0.2339105413835211,
-0.2572677222029636,
-0.1620439617193123,
-0.1620439174799243,
],

[-0.5278434299325454,
1.0000000000000000,
0.2360782850337852,
0.2360782850336694,
-0.2572431900217150,
],

[-0.5235077282962379,
1.0000000000000000,
0.2382461358519384,
0.2382461358518236,
-0.2572186576746783,
-0.1650050903760042,
-0.1650050903760042,
],]

# Instantiate the system
star = Star('A')
b = Planet('b', per = 1., t0 = 0., nl = 1)
c = Planet('c', per = 2., t0 = 0.)
system = System(star, b, c)

# Get the occultation light curve
time = np.linspace(0.2055, 0.21, 1000)
#time = np.linspace(0.2085, 0.2089, 100)[22:25]
system.compute(time)

# Plot the image
fig, ax = pl.subplots(1,3, figsize = (12,6))
for t in range(3):
  system.plot_image(t, system.b, [2], occultor_alpha = 0.1, ax = ax[t])
  for v in vertices[t]:
    ax[t].axvline(v, alpha = 0.5, color = 'r')
  ax[t].set_aspect('equal')

# Plot all of the occultations
system.plot_lightcurve()
pl.show()