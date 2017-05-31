#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`ppo.py` - Python interface to C
----------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .ppo import Star, Planet, System
import numpy as np
import matplotlib.pyplot as pl
import os
import emcee, corner
from tqdm import tqdm
MSUN = 1.988416e30
RSUN = 6.957e8
G = 6.67428e-11
MEARTH = 5.9722e24
REARTH = 6.3781e6
DAYSEC = 86400.
AUM = 1.49598e11

__all__ = ['Trappist1']

def h(i, Om):
  '''
  The angular momentum unit vector given the inclination `i`
  and the longitude of ascending node `Om`
  
  '''
  
  hx = np.sin(i) * np.sin(Om)
  hy = -np.sin(i) * np.cos(Om)
  hz = np.cos(i)
  return np.array([hx, hy, hz])

def LnLikelihood(p, imu, isig, imin):
  '''
  The likelihood of observing all seven planets transit given
  a set of angular momentum vectors for each planet and a value for the
  standard deviation of the distribution of longitude of ascending
  nodes.
  
  '''
  
  # Get the parameters; set Om_b = 0
  sigOm = p[0]
  x = np.array(p[0::3]); x[0] = 0
  y = np.array(p[1::3])
  z = np.array(p[2::3])
  n = len(x)
  
  # Initialize
  lnlike = 0
  Om = np.zeros(n) * np.nan
  Om[0] = sigOm
  
  # Enforce positive standard deviation for omega
  if sigOm < 0:
    return -np.inf, Om
    
  # Loop over all planets
  for j in range(n):
    
    # Enforce normal priors with unit standard deviation on x, y, and z
    lnlike -= 0.5 * (x[j] ** 2 + y[j] ** 2 + z[j] ** 2)
    
    # Normalize to get the unit vector
    norm = np.sqrt(x[j] ** 2 + y[j] ** 2 + z[j] ** 2)
    x[j] /= norm
    z[j] /= norm
    
    # Compute the inclination in the range [0-90]
    i = np.arccos(np.abs(z[j]))
    
    # Compute the inclination if the orbit precessed 90 degrees
    # and reject this sample if the planet doesn't transit
    i90 = np.arccos(np.abs(x[j]))
    if i90 < imin[j]:
      return -np.inf, Om
    
    # Probability given the data for either
    # positive or negative impact parameter
    # (we pick the highest likelihood one)
    l1 = 0.5 * ((i - imu[j]) / isig[j]) ** 2
    l2 = 0.5 * ((i - (np.pi - imu[j])) / isig[j]) ** 2
    lnlike -= min(l1, l2)

    # Compute omega and its probability given sigOm
    if j > 0:
      Om[j] = np.arcsin(x[j] / np.sin(i))
      lnlike -= 0.5 * (Om[j] / sigOm) ** 2
      lnlike -= np.log(np.sqrt(2 * np.pi) * sigOm)

  return lnlike, Om

def OmegaPrior(nwalk = 42, nsteps = 50000, nburn = 10000, thin = 10, calc = False):
  '''
  Computes the prior standard deviation of the distribution of longitude of ascending
  nodes for the planets in the TRAPPSIT-1 system. This yields sigma_Om ~ 0.4 degrees.
  
  '''
  
  # Parameters from Gillon et al. (2017) and Luger et al. (2017)
  imu = np.array([89.65, 89.670, 89.750, 89.86, 89.680, 89.710, 89.80]) * np.pi / 180
  isig = np.array([0.245, 0.170, 0.160, 0.110, 0.034, 0.025, 0.075]) * np.pi / 180
  imin = np.arccos(1. / np.array([20.50, 28.08, 39.55, 51.97, 68.4, 83.2, 109.]))
  ndim = len(imu) * 3
  
  # Initial values
  def get_guess():
    x0 = []
    # Initial values for x, y, and z for each planet, computed from
    # the observed inclination distribution and assuming the Omegas
    # are close to zero with standard deviation 0.1 degrees
    for m,s in zip(imu, isig):
      if np.random.rand() < 0.5:
        # Positive impact parameter
        x, y, z = h(m + s * np.random.randn(), 0.1 * np.random.randn() * np.pi / 180)
      else:
        # Negative impact parameter
        x, y, z = h(np.pi - m + s * np.random.randn(), 0.1 * np.random.randn() * np.pi / 180)
      x0.extend([x,y,z])
    # Initial value for sigma omega -- about 0.1 degrees
    x0[0] = np.abs(0.1 * np.random.randn() * np.pi / 180)
    return x0
  x0 = np.array([get_guess() for w in range(nwalk)])
  _, blobs0 = LnLikelihood(x0[0], imu, isig, imin)
  
  # Run the chain
  if calc or not os.path.exists('omega.npz'):
    sampler = emcee.EnsembleSampler(nwalk, ndim, LnLikelihood, args = (imu, isig, imin))
    for i in tqdm(sampler.sample(x0, iterations = nsteps, thin = thin, blobs0 = blobs0), total = nsteps):
      pass
    chain = sampler.chain
    blobs = sampler.blobs
    np.savez('omega.npz', chain = chain, blobs = blobs)
  else:
    data = np.load('omega.npz')
    blobs = data['blobs']
    chain = data['chain']
    
  # Plot the corner plot, with values in degrees
  labels = [r'$\sigma_\Omega$'] + [r'$\Omega_\mathrm{%s}$' % p for p in ['c', 'd', 'e', 'f', 'g', 'h']]
  blobs = np.swapaxes(blobs, 0, 1)
  blobs = blobs[:,(nburn // thin):,:].reshape(nwalk * (nsteps - nburn) // thin, -1)
  blobs *= 180 / np.pi
  fig = corner.corner(blobs, labels = labels, range = [(0, 2)] + [1 for j in range(len(labels) - 1)])
  fig.set_size_inches(9,9)
  
  # Appearance tweaks
  for ax in fig.axes:
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
      tick.set_fontsize(8)
    ax.set_ylabel(ax.get_ylabel(), fontsize = 14, fontweight = 'bold', labelpad = 35)
    ax.set_xlabel(ax.get_xlabel(), fontsize = 14, fontweight = 'bold', labelpad = 35)

  # Show
  fig.savefig('omega.pdf', bbox_inches = 'tight')

def Trappist1(nl = 11, polyeps1 = 1e-8, polyeps2 = 1e-15, ttvs = True, uncertainty = True, phasecurve = False):
  '''
  
  '''
  
  # Account for the uncertainty?
  if not uncertainty:
    N = lambda mu, sigma: mu
  else: 
    N = lambda mu, sigma: mu + sigma * np.random.randn()
    
  # Instantiate the star
  mstar = N(0.0802, 0.0073)
  rstar = N(0.117, 0.0036)
  star = Star('A', m = mstar, r = rstar, color = 'k')
  
  # Parameters from Gillon et al. (2017) and Luger et al. (2017)
  # Mass for `h` is currently unconstrained, so basing it loosely on 
  # the mass distribution for `d`, which has a similar radius.
  planets = [None for i in range(7)]
  names = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
  periods = [(1.51087081, 0.60e-6), (2.4218233, 0.17e-5), (4.049610, 0.63e-4), (6.099615, 0.11e-4), 
             (9.206690, 0.15e-4), (12.35294, 0.12e-3), (18.767, 0.004)]
  transits = [(7322.51736, 0.00010), (7282.80728, 0.00019), (7670.14165, 0.00035), (7660.37859, 0.00038),
              (7671.39767, 0.00023), (7665.34937, 0.00021), (7662.55284, 0.00037)]
  masses = [(0.85, 0.72), (1.38, 0.61), (0.41, 0.27), (0.62, 0.58), (0.68, 0.18), (1.34, 0.88), (0.4, 0.3)]
  inclinations = [(89.65, 0.245), (89.67, 0.17), (89.75, 0.16), (89.86, 0.11), (89.680, 0.034),
                  (89.710, 0.025), (89.80, 0.075)]
  depths = [(0.7266, 0.0088), (0.687, 0.010), (0.367, 0.017), (0.519, 0.026), (0.673, 0.023), 
            (0.782, 0.027), (0.752, 0.032)]
  colors = ['firebrick', 'coral', 'gold', 'mediumseagreen', 'turquoise', 'cornflowerblue', 'midnightblue']
  
  # These are eyeballed from Supplementary Figure 6 in Luger et al. (2017).
  # These are likely quite biased and model-specific. Need to re-think them.
  eccentricities = [(0.0005, 0.0001), (0.004, 0.001), (0.0004, 0.0003), (0.007, 0.0005), (0.009, 0.001), (0.004, 0.001), (0.003, 0.001)]
  
  # Instantiate the planets
  for i in range(7):
  
    # Period and time of transit
    per = N(*periods[i])
    trn0 = N(*transits[i])
  
    # Positive mass
    m = 0
    while m <= 0:
      m = N(*masses[i])
  
    # Inclination in range [0, 90]
    inc = N(*inclinations[i])
    if inc > 90:
      inc = 180 - inc
    
    # Longitude of ascending node in degrees
    # This is computed from a distribution with
    # standard deviation 0.4 degrees, which we
    # obtain in `OmegaPrior()` above. For
    # definiteness and WLOG, we choose Omega_b = 0.
    if i == 0:
      Omega = 0
    else:
      Omega = 0.4 * np.random.randn()
    
    # Longitude of pericenter (uniform over [0-360 deg])
    w = 360. * np.random.rand()
    
    # Eccentricity
    ecc = 1
    while (ecc < 0) or (ecc >= 1):
      ecc = N(*eccentricities[i])
    
    # Semi-major axis in AU from Kepler's law
    a = ((per * DAYSEC) ** 2 * G * (mstar * MSUN + m * MEARTH) / (4 * np.pi ** 2)) ** (1. / 3.) / AUM
  
    # Radius from Rp / Rstar
    mu = np.sqrt(depths[i][0] / 100)
    sig = 0.5 * depths[i][1] / 100 / mu
    RpRs = N(mu, sig)
    r = RpRs * rstar * RSUN / REARTH
  
    # Instantiate!
    planets[i] = Planet(names[i], m = m, per = per, inc = inc, a = a, r = r, trn0 = trn0, 
                        nl = nl, Omega = Omega, w = w, ecc = ecc, phasecurve = phasecurve,
                        color = colors[i])

  # Return the system
  return System(star, *planets, polyeps1 = polyeps1, polyeps2 = polyeps2, ttvs = ttvs)