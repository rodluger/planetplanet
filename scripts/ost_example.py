#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
ost_example.py
--------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.constants import *
from planetplanet.photo import Star, Planet, System
from planetplanet.detect import create_tophat_filter
import matplotlib.pyplot as pl
import numpy as np

"""
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
    system.compute(time, lambda1 = 5, lambda2 = 35)
    #system.plot_lightcurve(50.)

    # Create custom filter
    f50 = create_tophat_filter(10., 30., dlam = 0.1, Tput = 0.3, name = r"20 $\pm$5 $\mu$m")

    # Observe it (one exposure)
    system.observe(stack = 1, filter = f50, instrument = 'ost')
    pl.show()
"""

def Triple_bc():
  '''
  Simulate an observation of a triple occultation of TRAPPIST-1 `c` by `b`
  with OST at 50 microns.

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
  time = np.arange(252.75, 253.50, 5 * MINUTE)

  # Compute the light curve
  system.compute(time, lambda1 = 5, lambda2 = 35)

  # Let's re-center the time array for a prettier x axis
  system.A.time_hr -= time[0]
  system.A.time -= time[0]

  # Create custom filter
  f50 = create_tophat_filter(10., 30., dlam = 0.1, Tput = 0.3, name = r"OST 10-30 $\mu$m")

  # Observe it (one exposure)
  np.random.seed(1234)
  fig, ax = system.observe(stack = 1, filter = f50, instrument = 'ost')
  fig.set_size_inches(12, 6)
  fig.subplots_adjust(top = 0.7)
  ax.set_title("")

  # Plot the orbits of all bodies
  colors = ['k', 'firebrick', 'coral']
  for rect, index in zip([[0.11, 0.75, 0.15, 0.15], [0.315, 0.75, 0.15, 0.15], [0.405, 0.75, 0.15, 0.15], [0.71, 0.75, 0.15, 0.15]], [np.argmax(system.A.time_hr > 0.06), np.argmax(system.A.time_hr > 0.26), np.argmax(system.A.time_hr > 0.32), np.argmax(system.A.time_hr > 0.632)]):
    axxz = fig.add_axes(rect)
    f = np.linspace(0, 2 * np.pi, 1000)
    for j, b in enumerate(system.bodies):
      r = b.a * (1 - b.ecc ** 2) / (1 + b.ecc * np.cos(f))
      x = r * np.cos(b._w + f) - r * np.sin(b._w + f) * np.cos(b._inc) * np.sin(b._Omega)
      z = r * np.sin(b._w + f) * np.sin(b._inc)
      axxz.plot(x, z, color = 'gray', lw = 1)
      axxz.plot(b.x_hr[index], b.z_hr[index], 'o', color = colors[j], alpha = 1, markeredgecolor = 'k', zorder = 99, clip_on=False)
      for i in range(index, index - 100, -1):
        axxz.plot(b.x_hr[i], b.z_hr[i], 'o', color = colors[j], alpha = 0.01, markeredgecolor = colors[j], zorder = 99)
    axxz.axis('off')
    axxz.set_aspect(1)
    # Label the first image
    if rect[0] == 0.11:
      axxz.annotate("b", xy = (system.b.x_hr[index], system.b.z_hr[index]), va = "center", ha = "center", xytext = (12, 4), textcoords = "offset points", color = "firebrick",
                    fontweight = 'bold', fontsize = 8)
      axxz.annotate("c", xy = (system.c.x_hr[index], system.c.z_hr[index]), va = "center", ha = "center", xytext = (12, 4), textcoords = "offset points", color = "coral",
                    fontweight = 'bold', fontsize = 8)

  fig.savefig("../img/triple_bc_ost.pdf", bbox_inches = 'tight')

def Stacked_bc():
  '''
  Ten stacked exposures of `b` occulting `c`

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
             airless = True, phasecurve = False)

  # Instantiate `c`
  RpRs = np.sqrt(0.687 / 100)
  r = RpRs * rstar * RSUN / REARTH
  c = Planet('c', m = 1.38, per = 2.4218233, inc = 89.67, r = r, t0 = 0,
             Omega = 0, w = 0, ecc = 0, color = 'coral', tnight = 40., albedo = 0.3,
             airless = True, phasecurve = False)

  # Instantiate the system
  system = System(star, b, c, distance = 12, oversample = 10)

  # There's a triple occultation of `c` at this time
  time = np.arange(252.75, 253.50, 1 * MINUTE)

  # Compute and plot the light curve
  system.compute(time, lambda1 = 5, lambda2 = 35)

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

  # Create custom filter
  f50 = create_tophat_filter(10., 30., dlam = 0.1, Tput = 0.3, name = r"OST 10-30 $\mu$m")

  # Observe it (ten exposures)
  np.random.seed(123)
  fig, ax = system.observe(stack = 60, filter = f50, instrument = 'ost')
  fig.set_size_inches(12, 6)
  fig.subplots_adjust(top = 0.7)

  # Center on the event
  ax.set_xlim(-0.1, 0.1)
  ax.set_title("")

  # Plot the orbits of all bodies
  colors = ['k', 'firebrick', 'coral']
  rect = [0.433, 0.75, 0.15, 0.15]
  index = np.argmax(system.A.time_hr > 0)
  axxz = fig.add_axes(rect)
  f = np.linspace(0, 2 * np.pi, 1000)
  for j, b in enumerate(system.bodies):
    r = b.a * (1 - b.ecc ** 2) / (1 + b.ecc * np.cos(f))
    x = r * np.cos(b._w + f) - r * np.sin(b._w + f) * np.cos(b._inc) * np.sin(b._Omega)
    z = r * np.sin(b._w + f) * np.sin(b._inc)
    axxz.plot(x, z, color = 'gray', lw = 1)
    axxz.plot(b.x_hr[index], b.z_hr[index], 'o', color = colors[j], alpha = 1, markeredgecolor = 'k', zorder = 99, clip_on=False)
    for i in range(index, index - 100, -1):
      axxz.plot(b.x_hr[i], b.z_hr[i], 'o', color = colors[j], alpha = 0.01, markeredgecolor = colors[j], zorder = 99)
  axxz.axis('off')
  axxz.set_aspect(1)
  axxz.annotate("b", xy = (system.b.x_hr[index], system.b.z_hr[index]), va = "center", ha = "center", xytext = (18, 3), textcoords = "offset points", color = "firebrick",
                fontweight = 'bold', fontsize = 8)
  axxz.annotate("c", xy = (system.c.x_hr[index], system.c.z_hr[index]), va = "center", ha = "center", xytext = (18, 3), textcoords = "offset points", color = "coral",
                fontweight = 'bold', fontsize = 8)

  # Save
  fig.savefig("../img/stacked_bc_ost.pdf", bbox_inches = 'tight')

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
    system.compute(time, lambda1 = 1, lambda2 = 100, R=1000)
    #system.plot_lightcurve(50.)

    Tput = 0.3
    lams = np.linspace(5,80,100)
    dlam = 5.0

    SNRs = np.zeros_like(lams)
    signal = np.zeros_like(lams)
    noise = np.zeros_like(lams)

    for i in range(len(lams)):

        # Create custom filter
        f = create_tophat_filter(lams[i]-0.5*dlam, lams[i]+0.5*dlam, dlam = 0.1, Tput = Tput, name = "$%.1f \pm %.1f \mu$m" %(lams[i], 0.5*dlam))

        # Observe it (one exposure)
        system.observe(stack = 1, filter = f, instrument = 'ost')

        # get SNR on feature
        SNRs[i] = system.filter.lightcurve.event_SNRs[0]

        # get ppm noise on lc
        noise[i] = system.filter.lightcurve.event_noise[0]

        # get ppm amplitude on event
        signal[i] = system.filter.lightcurve.event_signal[0]

        # let's not have all these plots show up
        pl.clf()
        pl.close()

    fig, ax = pl.subplots(figsize=(6,5))
    ax.set_xlabel(r"Wavelength [$\mu$m]", fontweight = 'bold', fontsize = 10)
    ax.set_ylabel("SNR", fontweight = 'bold', fontsize = 10)
    ax.plot(lams, SNRs, ls="-", color="k")

    #""" Plot the noise and signal amplitude as a function of wavelength
    ax2 = ax.twinx()
    #ax2.semilogy()
    ax2.set_ylabel("Signal [ppm]", fontweight = 'bold', fontsize = 10, rotation = 270, labelpad = 15)
    ax2.plot(lams, signal, ls="--", color="C0", label="signal")
    ax2.plot(lams, noise, ls="--", color="C3", label="noise")
    ax2.legend(loc=9)
    #"""

    # Save
    fig.savefig("../img/snr_bc_ost.pdf", bbox_inches = 'tight')

    pl.show()

#Triple_bc()
#Stacked_bc()
SNR_v_wavelength()