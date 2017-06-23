from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(99)

# Instantiate the Trappist-1 system
system = Trappist1(sample = True)

# Get the next occultation of b
t = system.next_occultation(100, system.b)

# Get the light curve around that point
time = np.arange(t - 0.1, t + 0.1, 0.01 / 1440.)
system.compute(time, lambda1 = 4, lambda2 = 30, R = 3000)
system.plot_lightcurve()
pl.show()