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
from planetplanet.photo import Star, Planet, System, Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

vertices = [[-4.1894798819512191,
-2.0437458961562811,
-3.9473799137185788,
-2.1231321097257405,
-3.9577115446967541,
-3.2993770706601921,
],

[-4.1055166255840643,
-1.9597826397891258,
-3.8968980723907398,
-2.0505096226210138,
-3.7420059879226288,
-3.3762553925035377,
],

[-4.0215533423894367,
-1.8758193565944983,
-3.8438571521412888,
-1.9782098284259473,
]]


# Instantiate the Trappist-1 system
system = Trappist1(sample = True, nbody = False, phasecurve = False, adaptive = True)
time = np.linspace(2.1, 2.24, 1000)#[783:786]

system.compute(time, lambda1 = 5, lambda2 = 15, R = 1000)

# Plot the image
fig, ax = pl.subplots(1,3, figsize = (12,6))
for t in range(3):
  system.plot_image(t, system.A, [3,5], occultor_alpha = 0.1, ax = ax[t])
  for v in vertices[t]:
    ax[t].axvline(v, alpha = 0.5, color = 'r')
  ax[t].set_xlim(-4.5,-1.5)
  ax[t].set_ylim(-7,-4)
  ax[t].set_aspect('equal')
  

# Plot all of the occultations
system.plot_lightcurve()
pl.show()