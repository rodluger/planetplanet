#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
exomoon.py
----------

*Not yet working.*

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, Moon, System
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

# Instantiate the star
star = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [0.4, 0.26])

# Planet b
b = Planet('b', m = 1., per = 10, inc = 90., r = 1., t0 = 0, 
           nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# Moon
m = Moon('m', 'b', m = 0.001, per = 1, inc = 90., r = 0.5, t0 = 0.5, 
         nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# Compute the light curve
system = System(star, b, m, nbody = True, integrator = 'ias15', timestep = MINUTE)
time = np.arange(-6, 6, 1 * MINUTE)
system.compute(time)

# Plot
system.plot_lightcurve()
pl.show()