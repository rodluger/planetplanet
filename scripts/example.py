#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
example.py
----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, ttvs = False, phasecurve = False, adaptive = True)

# Get the occultation light curves for the first 10 days
time = np.linspace(0, 1, 1000)
system.compute(time)

# Plot all of the occultations
system.plot_lightcurve()
a,b,c,d,e = system.plot_occultations('A')
pl.show()