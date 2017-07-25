#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
occultation.py
--------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, phasecurve = True, airless = True, nbody = True, seed = 999)

# Give `c` a large latitudinal offset in its hotspot just for fun
system.c.dlambda = 30

# Compute an occultation by `b`
time = np.linspace(9552.9364, 9552.9564, 100)
system.compute(time)

# Plot the occultation
system.plot_occultation('c', 9552.95)
pl.show()
