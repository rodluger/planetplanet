'''
trappist_eclipse_estimates.py
-----------------------------

Computes and plots the SNR on the planet when observed
in secondary eclipse with JWST/MIRI, considering photon
and thermal noise.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from planetplanet.detect import jwst
import matplotlib.pyplot as plt
import numpy as np
np.random.seed(1234)

jwst.estimate_eclipse_snr(tint = 36.4*60.,
                          nout = 4.0,
                          Tstar = 2560.,
                          Tplan = 400.,
                          Rs = 0.12,
                          Rp = 1.086,
                          d = 12.2,
                          verbose=True,
                          plot=True)
