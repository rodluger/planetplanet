#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
transit.py
----------

A simple transit light curve.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np


# Instantiate the star
star = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [0.4, 0.26])

# Planet b
b = Planet('b', m = 1, per = 2, inc = 90.4, r = 10., t0 = 0, 
           nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# Compute the light curve, no optimization
system = System(star, b, batmanopt = False)
time = np.arange(-0.05, 0.05, MINUTE)
system.compute(time)
flux1 = system.A.flux[:,0] / system.A.flux[0,0]

# Compute the light curve w/ batman optimization
system = System(star, b, batmanopt = True)
system.compute(time)
flux2 = system.A.flux[:,0] / system.A.flux[0,0]

# Plot it
pl.plot(system.A.time, flux1, label = 'Standard')
pl.plot(system.A.time, flux2, '--', label = 'Batman')
pl.legend()
pl.show()