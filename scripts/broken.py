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

# Instantiate the system
star = Star('A')
b = Planet('b', per = 1., t0 = 0.)
c = Planet('c', per = 2., t0 = 0.)
system = System(star, b, c)

# Get the occultation light curve
time = np.linspace(0.2055, 0.21, 10000)
system.compute(time)

# Plot all of the occultations
system.plot_lightcurve()
pl.show()