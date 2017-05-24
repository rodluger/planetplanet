#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test.py
-------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.trappist1 import Trappist1
import matplotlib.pyplot as pl
import numpy as np

system = Trappist1(uncertainty = True)

# Get the light curves 
time = np.linspace(0, 10, 100000)
system.compute(time)

# Plot the occultations of `b`
system.plot_occultations('b')
pl.show()