#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
validation.py |github|
----------------------

Validates the code against :py:obj:`batman` for transit light curves
and against a direct (brute force) integration of the radiance map for
planet-planet occultations.

  .. plot::
     :align: center
     
     from scripts import validation
     validation._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/validation.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet.constants import *
from planetplanet import Star, Planet, System
import matplotlib.pyplot as pl
import numpy as np
try:
  import batman
except:
  batman = None
from tqdm import tqdm

def _test():
  '''
  
  '''
  
  ValidateTransits()
  ValidateOccultations()
  pl.show()
  
def ZenithAngle(x, y, r, theta):
  '''
  Compute the zenith angle.
  
  '''
  
  # Normalize
  x = x / r
  y = y / r
  x2 = x * x
  y2 = y * y
  
  # This is a solution to a quadratic equation in z = sin(za) **  2 
  z = 0.5 * ((1 - 2 * x2 - y2) * np.cos(2 * theta) + 2 * x * np.sqrt(1 - x2 - y2) * np.sin(2 * theta) + y2 + 1)
  
  # Where are we relative to the terminator?
  xterm = np.sin(theta) * np.sqrt(np.abs(1 - y2))
  
  # Solve for the zenith angle
  if np.abs(theta) <= np.pi / 2:
    if (x <= xterm):
      return np.arcsin(np.sqrt(z)) 
    else:
      return np.pi - np.arcsin(np.sqrt(z))
  else:  
    if (x >= -xterm):
      return np.arcsin(np.sqrt(z)) 
    else:
      return np.pi - np.arcsin(np.sqrt(z))

def Radiance(z, irrad, lam = 15, albedo = 0.3, tnight = 40):
  '''
  Compute the radiance at a point on the planet's surface
  at a given wavelength.
  
  '''
  
  # Cosine law
  if (z < np.pi / 2):
    temp = ((irrad * np.cos(z) * (1 - albedo)) / SBOLTZ) ** 0.25
    if (temp < tnight): 
      temp = tnight
  else:
    temp = tnight

  # Planck's law
  lam /= 1e6
  a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
  b = HPLANCK * CLIGHT / (lam * KBOLTZ * temp)
  return a / (np.exp(b) - 1)

def ValidateTransits():
  '''
  
  '''
  
  if batman is None:
    print("Please install the `batman` package.")
    return
  
  # System params
  time = np.arange(-0.12, 0.12, 0.001)
  mstar = 1.
  rstar = 1.
  limbdark = [0.4, 0.26]
  per = 5.
  inc = 90.
  r = 10.
  t0 = 0.
  w = 60.
  ecc = 0.3

  # planetplanet
  star = Star('A', m = mstar, r = rstar, nz = 99, limbdark = limbdark)
  b = Planet('b', m = 0., per = per, inc = inc, r = r, t0 = t0, 
             nz = 1, Omega = 0., w = w, ecc = ecc, phasecurve = False)
  system = System(star, b, batmanopt = False, circleopt = False)
  system.compute(time)
  flux_pp = system.A.flux[:,0]
  flux_pp /= flux_pp[0]

  # batman
  params = batman.TransitParams()
  params.t0 = t0
  params.per = per
  params.inc = inc
  params.ecc = ecc
  params.w = w - 180.
  params.limb_dark = "quadratic"
  params.u = limbdark
  params.rp = r / (rstar * RSUNREARTH)
  params.a = ((b.per) ** 2 * GEARTH * (mstar * MSUNMEARTH) / (4 * np.pi ** 2)) ** (1. / 3.) / (rstar * RSUNREARTH)
  m = batman.TransitModel(params, time)
  flux_bm = m.light_curve(params)
  
  # Ensure the two are equal to within 0.1 ppm
  assert np.max(np.abs((flux_pp - flux_bm) * 1e6)) < 0.1, "Flux mismatch between `planetplanet` and `batman`."
  
  # Plot the comparison
  fig, ax = pl.subplots(2, sharex = True)
  ax[0].plot(time, flux_pp, color = 'b', label = 'planetplanet (pp)')
  ax[0].plot(time, flux_bm, color = 'g', ls = '--', label = 'batman (bm)')
  ax[1].plot(time, (flux_pp - flux_bm) * 1e6, color = 'k')
  ax[0].legend(loc = 9)
  ax[0].set_ylabel('Flux', fontweight = 'bold')
  ax[1].set_ylabel('pp - bm [ppm]', fontweight = 'bold')
  ax[1].set_xlabel('Time [days]', fontweight = 'bold')
  
  return fig, ax

def ValidateOccultations():
  '''
  
  '''
  
  # Instantiate the star
  mstar = 0.0802
  rstar = 0.121
  teff = (0.000524 * LSUN / (4 * np.pi * (rstar * RSUN) ** 2 * SBOLTZ)) ** 0.25
  star = Star('A', m = mstar, r = rstar, teff = teff, color = 'k')
  
  # Instantiate `c`
  RpRs = np.sqrt(0.687 / 100)
  r = RpRs * rstar * RSUN / REARTH
  c = Planet('c', m = 1.38, per = 2.4218233, inc = 89.67 - 0.05, r = r, t0 = 0, 
             Omega = 0, w = 0, ecc = 0, color = 'coral', tnight = 40., albedo = 0.3, 
             airless = True, phasecurve = False, nz = 31)

  # Instantiate `d`
  RpRs = np.sqrt(0.367 / 100)
  r = RpRs * rstar * RSUN / REARTH    
  d = Planet('d', m = 0.41, per = 4.049610, inc = 89.75 + 0.16, r = r, t0 = 0, 
             Omega = 0, w = 0, ecc = 0, color = 'firebrick', tnight = 40., albedo = 0.3, 
             airless = True, phasecurve = False)

  # Instantiate the system
  system = System(star, c, d, distance = 12, oversample = 1, nbody = False)

  # There's a triple occultation of `c` at this time
  time = np.arange(-259.684 + 2 * 0.00025, -259.665, 0.01 * MINUTE)
  
  # Compute the light curve using planetplanet
  system = System(star, c, d)
  system.compute(time, lambda2 = 15)
  flux_pp = 1 + np.array(c.flux[:,-1]) / c.total_flux[-1]
  
  # Rescale the time array
  time = (time - np.nanmedian(time)) / MINUTE
  
  # Now compute the light curve by brute force direct integration
  flux_bf = np.zeros_like(time)
  for t in tqdm(range(len(time))):
    
    # Grid up the planet
    xarr = c.x[t] + np.linspace(-c.r, c.r, 30)
    yarr = c.y[t] + np.linspace(-c.r, c.r, 30)
    
    rad = np.zeros((len(xarr), len(yarr)))
    for i, x in enumerate(xarr):
      for j, y in enumerate(yarr):
        
        # Are we outside the planet?
        if (x - c.x[t]) ** 2 + (y - c.y[t]) ** 2 >= c.r ** 2:
          continue
      
        # Are we outside the occultor?
        if (x - d.x[t]) ** 2 + (y - d.y[t]) ** 2 >= d.r ** 2:
          continue

        # Get the orbital phase
        theta = np.arctan2(c.z[t], c.x[t])
        
        # Zenith angle
        z = ZenithAngle(x - c.x[t], y - c.y[t], c.r, theta)
        
        # Get the irradiance on the planet
        d2 = c.x[t] ** 2 + c.y[t] ** 2 + c.z[t] ** 2
        irrad = star._r ** 2 * SBOLTZ * star.teff ** 4 / d2
        
        # Get the planet's radiance
        rad[i,j] = Radiance(z, irrad)

    flux_bf[t] = -np.sum(rad)

  # Normalize it
  dbf = -np.min(flux_bf)
  dpp = 1 - np.min(flux_pp)
  flux_bf = flux_bf * dpp / dbf + 1
  
  # Ensure the two are equal to within one percent
  assert np.max(np.abs((flux_pp - flux_bf))) < 0.01, "Flux mismatch between `planetplanet` and brute force integration."
  
  # Plot the light curve
  fig = pl.figure()
  pl.plot(time, flux_pp, '-', label = 'planetplanet')
  pl.plot(time, flux_bf, '-', label = 'brute force')
  pl.ylabel('Flux', fontweight = 'bold')
  pl.xlabel('Time [days]', fontweight = 'bold')
  pl.legend()
  
  return fig, pl.gca()

if __name__ == '__main__':
  ValidateTransits()
  ValidateOccultations()
  pl.show()