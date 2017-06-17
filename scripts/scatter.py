#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
scatter.py
----------

Computes all occultations that occur in the TRAPPIST-1 system over a 
3 year time period for a random draw from the prior. Plots each
occultation as a circle in a top-view of the system; the circle size, 
transparency, and color indicate the duration, impact parameter, and 
occulting body, respectively.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, nbody = False)
system.scatter_plot(0, 365 * 3)
pl.show()