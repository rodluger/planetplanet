#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
mutual_transit.py
-----------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the star
star = Star('A', m = 0.1, r = 0.1, nl = 31, color = 'k')

# Planet b
b = Planet('b', m = 1, per = 3, inc = 90., r = 3, trn0 = 0, 
           nl = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

# Planet c
c = Planet('c', m = 1, per = 100, inc = 89.99, r = 3., trn0 = 0., 
           nl = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'b')

# Planet c
system = System(star, b, c, ttvs = False, adaptive = False)

# Get the occultation light curves
time = np.linspace(-1, 1, 10000)
system.compute(time)

# Plot all of the occultations
system.plot_occultations('A')
pl.show()