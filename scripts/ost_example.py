#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
ost_example.py
---------------

Simulate an observation of a triple occultation of TRAPPIST-1 `c` by `b`
with OST "MIRI" at 15 microns.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.constants import *
from planetplanet.photo import Star, Planet, System
from planetplanet.detect import create_tophat_filter
import matplotlib.pyplot as pl
import numpy as np

def Triple_bc():
    '''
    '''

    # Instantiate the star
    mstar = 0.0802
    rstar = 0.121
    teff = (0.000524 * LSUN / (4 * np.pi * (rstar * RSUN) ** 2 * SBOLTZ)) ** 0.25
    star = Star('A', m = mstar, r = rstar, teff = teff, color = 'k')

    # Instantiate `b`
    RpRs = np.sqrt(0.7266 / 100)
    r = RpRs * rstar * RSUN / REARTH
    b = Planet('b', m = 0.85, per = 1.51087081, inc = 89.65, r = r, t0 = 0,
               Omega = 0, w = 0, ecc = 0, color = 'firebrick', tnight = 40., albedo = 0.3,
               airless = True, phasecurve = True)

    # Instantiate `c`
    RpRs = np.sqrt(0.687 / 100)
    r = RpRs * rstar * RSUN / REARTH
    c = Planet('c', m = 1.38, per = 2.4218233, inc = 89.67, r = r, t0 = 0,
               Omega = 0, w = 0, ecc = 0, color = 'coral', tnight = 40., albedo = 0.3,
               airless = True, phasecurve = True)

    # Instantiate the system
    system = System(star, b, c, distance = 12, oversample = 10)

    # There's a triple occultation of `c` at this time
    time = np.arange(252.75, 253.50, 10 * MINUTE)

    # Compute and plot the light curve
    system.compute(time, lambda1 = 10, lambda2 = 60)
    system.plot_lightcurve(50.)

    # Create custom filter
    f50 = create_tophat_filter(45., 55., dlam = 0.1, Tput = 0.3, name = r"50 $\pm$5 $\mu$m")

    # Observe it (one exposure)
    system.observe(stack = 1, filter = f50, instrument = 'ost')
    pl.show()

def SNR_v_wavelength():
    '''
    '''
    # Instantiate the star
    mstar = 0.0802
    rstar = 0.121
    teff = (0.000524 * LSUN / (4 * np.pi * (rstar * RSUN) ** 2 * SBOLTZ)) ** 0.25
    star = Star('A', m = mstar, r = rstar, teff = teff, color = 'k')

    # Instantiate `b`
    RpRs = np.sqrt(0.7266 / 100)
    r = RpRs * rstar * RSUN / REARTH
    b = Planet('b', m = 0.85, per = 1.51087081, inc = 89.65, r = r, t0 = 0,
               Omega = 0, w = 0, ecc = 0, color = 'firebrick', tnight = 40., albedo = 0.3,
               airless = True, phasecurve = True)

    # Instantiate `c`
    RpRs = np.sqrt(0.687 / 100)
    r = RpRs * rstar * RSUN / REARTH
    c = Planet('c', m = 1.38, per = 2.4218233, inc = 89.67, r = r, t0 = 0,
               Omega = 0, w = 0, ecc = 0, color = 'coral', tnight = 40., albedo = 0.3,
               airless = True, phasecurve = True)

    # Instantiate the system
    system = System(star, b, c, distance = 12, oversample = 10)

    # There's a triple occultation of `c` at this time
    time = np.arange(252.98, 253.05, 5 * MINUTE)

    # Compute and plot the light curve
    system.compute(time, lambda1 = 5, lambda2 = 100)
    #system.plot_lightcurve(50.)

    Tput = 0.3
    lams = np.linspace(5,80,30)
    dlam = 5.0
    SNRs = np.zeros_like(lams)

    for i in range(len(lams)):

        # Create custom filter
        f = create_tophat_filter(lams[i]-0.5*dlam, lams[i]+0.5*dlam, dlam = 0.1, Tput = Tput, name = "$%.1f \pm %.1f \mu$m" %(lams[i], 0.5*dlam))

        # Observe it (one exposure)
        system.observe(stack = 1, filter = f, instrument = 'ost')

        # get SNR on feature
        SNRs[i] = system.filter.lightcurve.event_SNRs[0]

        pl.clf()
        pl.close()

    fig, ax = pl.subplots(figsize=(12,6))
    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel("SNR")
    ax.plot(lams, SNRs, "o-", color="k")
    pl.show()

Triple_bc()
SNR_v_wavelength()
