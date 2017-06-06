#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
flower.py
---------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the star
star = Star('A', m = 0.1, r = 0.1, nl = 21, color = 'k')

# Planet b
b = Planet('b', m = 1, per = 3, inc = 89.6, r = 3., trn0 = 0, 
           nl = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

# Planet c
c = Planet('c', m = 1, per = 3, inc = 89.6, r = 3., trn0 = 0, 
           nl = 11, Omega = 90, w = 0., ecc = 0., phasecurve = False, color = 'b')

# Planet c
d = Planet('d', m = 1, per = 3, inc = 89.6, r = 3., trn0 = 0, 
           nl = 11, Omega = 180, w = 0., ecc = 0., phasecurve = False, color = 'b')
           
# Planet c
e = Planet('e', m = 1, per = 3, inc = 89.6, r = 3., trn0 = 0, 
           nl = 11, Omega = 270, w = 0., ecc = 0., phasecurve = False, color = 'b')

# System
system = System(star, b, c, d, e, ttvs = False, adaptive = True)

# Get the occultation light curves
time = np.linspace(-0.02, 0.02, 1000)
system.compute(time)
system.plot_occultations('A') #, gifname = 'flower')
pl.show()