#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
spitzer_example.py
------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.constants import *
from planetplanet.photo import Star, Planet, System
from planetplanet.detect import jwst
import matplotlib.pyplot as pl
import numpy as np

def Stacked_bc_all_filters(N = 100, albedo = 0.0, airless = True):
   '''
   `N` stacked exposures of `b` occulting `c` observed in both Warm Spitzer filters

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
              Omega = 0, w = 0, ecc = 0, color = 'firebrick', tnight = 40., albedo = albedo,
              airless = airless, phasecurve = False)

   # Instantiate `c`
   RpRs = np.sqrt(0.687 / 100)
   r = RpRs * rstar * RSUN / REARTH
   c = Planet('c', m = 1.38, per = 2.4218233, inc = 89.67, r = r, t0 = 0,
              Omega = 0, w = 0, ecc = 0, color = 'coral', tnight = 40., albedo = albedo,
              airless = airless, phasecurve = False)

   # Instantiate the system
   system = System(star, b, c, distance = 12, oversample = 10)

   # There's a triple occultation of `c` at this time
   time = np.arange(252.75, 253.50, 5 * MINUTE)

   # Compute and plot the light curve
   system.compute(time, lambda1 = 1, lambda2 = 15)

   # Hackishly remove the other occultations and secondary eclipses
   # so we can plot just the occultation of `c` by `b`
   a = np.argmax(time > 253.013 - 0.025)
   b = np.argmax(time > 253.013 + 0.025)
   system.c.flux[:a] = system.c.flux[0]
   system.b.flux[:a] = system.b.flux[0]
   system.c.flux[b:] = system.c.flux[0]
   system.b.flux[b:] = system.b.flux[0]
   system.c.occultor[:a] = 0
   system.c.occultor[b:] = 0
   system.b.occultor[:a] = 0
   system.b.occultor[b:] = 0
   a = np.argmax(system.time_hr > 253.013 - 0.025)
   b = np.argmax(system.time_hr > 253.013 + 0.025)
   system.c.flux_hr[:a] = system.c.flux_hr[0]
   system.b.flux_hr[:a] = system.b.flux_hr[0]
   system.c.flux_hr[b:] = system.c.flux_hr[0]
   system.b.flux_hr[b:] = system.b.flux_hr[0]
   system.c.occultor_hr[:a] = 0
   system.c.occultor_hr[b:] = 0
   system.b.occultor_hr[:a] = 0
   system.b.occultor_hr[b:] = 0

   # Let's re-center the time array for a prettier x axis
   system.A.time -= 253.013
   system.A.time_hr -= 253.013

   # Get all MIRI filters
   filters = jwst.get_spitzer_filter_wheel(warm = True)

   print("Planet c Teff : %.2f K" %system.c.teff)

   # Loop over all MIRI filters
   SNRs = np.zeros(len(filters))
   wls = np.zeros(len(filters))
   for i in range(len(filters)):

       # Observe it (ten exposures)
       np.random.seed(123)
       fig, ax = system.observe(stack = N, filter = filters[i], instrument = 'spitzer')

       SNRs[i] = system.filter.lightcurve.event_SNRs[0]
       wls[i] = system.filter.eff_wl
       print("SNR: %.3f" %SNRs[i])

       # Save
       #fig.savefig("../img/stacked_bc.pdf", bbox_inches = 'tight')
       #pl.show()
       pl.clf()
       pl.close()

   fig, ax = pl.subplots(figsize=(10,7))
   ax.plot(wls, SNRs, "-o", color="k")
   ax.set_xlabel(r"Wavelength [$\mu$m]", fontweight = 'bold', fontsize = 25)
   ax.set_ylabel("SNR", fontweight = 'bold', fontsize = 25)

   wl_filt, dwl_filt, tputs, names = jwst.readin_miri_filters()
   #jwst.plot_miri_filters(ax, wl_filt, tputs, names, ylim=[0.0,1.0], leg=True)
   #ax2 = fig.get_axes()[1]
   #ax2.set_ylabel("Throughput", fontweight = 'bold', fontsize = 25)
   ax.tick_params(labelsize=20)
   #ax2.tick_params(labelsize=20)
   #fig.savefig("../img/stacked_bc_all_filters.pdf", bbox_inches = 'tight')
   pl.show()

Stacked_bc_all_filters(N = 100000, airless = False, albedo = 0.0)
Stacked_bc_all_filters(N = 10000, airless = True, albedo = 0.0)
Stacked_bc_all_filters(N = 100, airless = False, albedo = -15)
