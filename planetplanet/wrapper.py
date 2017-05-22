#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`wrapper.py` - C wrapper
--------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import ctypes
import numpy as np
np.seterr(invalid = 'ignore')
import os
from numpy.ctypeslib import ndpointer, as_ctypes
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
rdbu = pl.get_cmap('RdBu_r')
greys = pl.get_cmap('Greys')
plasma = pl.get_cmap('plasma')

# Define constants
AUREARTH = 23454.9271
SEARTH = 1.361e3
MDFAST = 0
NEWTON = 1

DOUBLEPP = np.ctypeslib.ndpointer(dtype = np.uintp, ndim=1, flags='C')

# Load the library
try:
  ppo = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ppo.so'))
except:
  raise Exception("Can't find `ppo.so`; please run `make` to compile it.")

class Settings(ctypes.Structure):
  '''
  The class that contains the model settings.
  
  :param float keptol: Kepler solver tolerance. Default `1.e-15`
  :param int maxkepiter: Maximum number of Kepler solver iterations. Default `100`
  :param str kepsolver: Kepler solver (`newton` | `mdfast`). Default `newton`
  :param float polyeps1: Tolerance in the polynomial root-finding routine. Default `2.0e-6`
  :param float polyeps2: Tolerance in the polynomial root-finding routine. Default `6.0e-9`
  :param int maxpolyiter: Maximum number of root finding iterations. Default `100`
  
  '''
  
  _fields_ = [("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("kepsolver", ctypes.c_int),
              ("polyeps1", ctypes.c_double),
              ("polyeps2", ctypes.c_double),
              ("maxpolyiter", ctypes.c_int)]
  
  def __init__(self, **kwargs):
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = eval(kwargs.pop('kepsolver', 'newton').upper())
    self.polyeps1 = kwargs.pop('polyeps1', 2.0e-6)
    self.polyeps2 = kwargs.pop('polyeps2', 6.0e-9)
    self.maxpolyiter = kwargs.pop('maxpolyiter', 100)

class Planet(ctypes.Structure):
  '''
  The class containing all the input planet parameters.
  
  :param float per: Orbital period in days. Default `1.51087081`
  :param float inc: Orbital inclination in degrees. Default `89.65`
  :param float ecc: Orbital eccentricity. Default `0.`
  :param float w: Longitude of pericenter in degrees. `0.`
  :param float a: Semi-major axis in AU. Default `0.01111`
  :param float t0: Time of first transit in days. Default `7322.51736`
  :param float r: Planet radius in Earth radii. Default `1.086`
  :param float albedo: Planet albedo. Default `0.3`
  :param float irrad: Stellar irradiation at the planet's distance in units \
         of the solar constant (1370 W/m^2). Default `0.3`
  :param bool phasecurve: Compute the full phase curve? Default `False`
  :param int nl: Number of latitude slices. Default `11`
  
  '''
  
  _fields_ = [("per", ctypes.c_double),
              ("inc", ctypes.c_double),
              ("ecc", ctypes.c_double),
              ("w", ctypes.c_double),
              ("a", ctypes.c_double),
              ("t0", ctypes.c_double),
              ("r", ctypes.c_double),
              ("albedo", ctypes.c_double),
              ("irrad", ctypes.c_double),
              ("phasecurve", ctypes.c_int),
              ("nl", ctypes.c_int),
              ("t", ctypes.c_int),
              ("nt", ctypes.c_int),
              ("nw", ctypes.c_int),
              ("_time", ctypes.POINTER(ctypes.c_double)),
              ("_wavelength", ctypes.POINTER(ctypes.c_double)),
              ("_x", ctypes.POINTER(ctypes.c_double)),
              ("_y", ctypes.POINTER(ctypes.c_double)),
              ("_z", ctypes.POINTER(ctypes.c_double)),
              ("_occultor", ctypes.POINTER(ctypes.c_int)),
              ("_flux", ctypes.POINTER(ctypes.c_double))]
              
  def __init__(self, name, **kwargs):
  
    # User
    self.name = name
    self.per = kwargs.pop('per', 1.51087081)
    self.inc = kwargs.pop('inc', 89.65) * np.pi / 180.
    self.ecc = kwargs.pop('ecc', 0.)
    self.w = kwargs.pop('w', 0.) * np.pi / 180.
    self.a = kwargs.pop('a', 0.01111) * AUREARTH
    self.r = kwargs.pop('r', 1.086)
    self.albedo = kwargs.pop('albedo', 0.3)
    self.irrad = kwargs.pop('irrad', 4.25) * SEARTH
    self.phasecurve = int(kwargs.pop('phasecurve', False))
    self.nl = kwargs.pop('nl', 11)
    
    # System
    self.t = 0
    self.nt = 0
    self.nw = 0
    self._inds = []
    self._computed = False
    
    # Compute the time of pericenter passage (e.g. Shields et al. 2015) 
    fi = (3. * np.pi / 2.) - self.w + np.pi
    tperi0 = (self.per * np.sqrt(1. - self.ecc * self.ecc) / (2. * np.pi) * (self.ecc * np.sin(fi) / 
             (1. + self.ecc * np.cos(fi)) - 2. / np.sqrt(1. - self.ecc * self.ecc) * 
             np.arctan2(np.sqrt(1. - self.ecc * self.ecc) * np.tan(fi/2.), 1. + self.ecc)))
    
    # We define the mean anomaly to be zero at t = t0 = trn0 + tperi0
    self.t0 = kwargs.pop('trn0', 7322.51736) + tperi0
  
  def plot_lightcurve(self):
    '''
    
    '''
    
    fig, ax = pl.subplots(1, figsize = (12, 4))
    ax.plot(self.time, self.flux[:, 0], 'b-')
    ax.plot(self.time, self.flux[:, self.nw // 2], 'g-')
    ax.plot(self.time, self.flux[:, -1], 'r-')
    ax.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    ax.set_ylabel('Occulted Flux [W/sr]', fontweight = 'bold', fontsize = 10)
    ax.get_yaxis().set_major_locator(MaxNLocator(4))
    ax.get_xaxis().set_major_locator(MaxNLocator(4))
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
      tick.set_fontsize(8)
    
    return fig, ax

class System(object):

  def __init__(self, *planets, **kwargs):
    '''
    
    '''
    
    self.planets = planets
    self.settings = Settings(**kwargs)
    self._names = np.array([p.name for p in self.planets])
  
  def compute(self, time, lambda1 = 5, lambda2 = 15, R = 100):
    '''
    
    '''
    
    # Compute the wavelength grid
    wav = [lambda1]
    while(wav[-1] < lambda2):
      wav.append(wav[-1] + wav[-1] / R) 
    wavelength = np.array(wav) * 1e-6
  
    # Dimensions
    n = len(self.planets)
    nt = len(time)
    nw = len(wavelength)

    # Initialize the C interface
    Orbit = ppo.OrbitXYZ
    Orbit.restype = ctypes.c_int
    Orbit.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(Planet)), Settings]
    Flux = ppo.Flux
    Flux.restype = ctypes.c_int
    Flux.argtypes = [ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nt),
                     ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nw),
                     ctypes.c_int, ctypes.POINTER(ctypes.POINTER(Planet)),
                     Settings]
  
    # Allocate memory for all the arrays
    for planet in self.planets:
      planet.time = np.zeros(nt)
      planet._time = np.ctypeslib.as_ctypes(planet.time)
      planet.wavelength = np.zeros(nw)
      planet._wavelength = np.ctypeslib.as_ctypes(planet.wavelength)
      planet.x = np.zeros(nt)
      planet._x = np.ctypeslib.as_ctypes(planet.x)
      planet.y = np.zeros(nt)
      planet._y = np.ctypeslib.as_ctypes(planet.y)
      planet.z = np.zeros(nt)
      planet._z = np.ctypeslib.as_ctypes(planet.z)
      planet.occultor = np.zeros(nt, dtype = 'int32')
      planet._occultor = np.ctypeslib.as_ctypes(planet.occultor)
      # HACK: The fact that flux is 2d is a nightmare for ctypes. We will
      # treat it as a 1d array within C and keep track of the row/column
      # indices by hand...
      planet.flux = np.zeros((nt, nw))
      planet._flux1d = planet.flux.reshape(nt * nw)
      planet._flux = np.ctypeslib.as_ctypes(planet._flux1d)

    # A pointer to a pointer to `Planet`. This is an array of `n` `Planet` instances, 
    # passed by reference. The contents can all be accessed through `planets`
    ptr_planets = (ctypes.POINTER(Planet) * n)(*[ctypes.pointer(p) for p in self.planets])

    # Call the light curve routine
    Flux(nt, np.ctypeslib.as_ctypes(time), nw, np.ctypeslib.as_ctypes(wavelength), n, ptr_planets, self.settings)
  
    # Loop over all planets and store each occultation event as a separate attribute
    for planet in self.planets:
    
      # Set the flag
      planet._computed = True
    
      # Identify the different events
      inds = np.where(planet.occultor >= 0)[0]
      si = np.concatenate(([0], inds[np.where(np.diff(inds) > 1)] + 1, [nt]))

      # Loop over the events
      for i in range(len(si) - 1):
  
        # Split the light curve, trim it, and add a little padding
        t = time[si[i]:si[i+1]]
        f = planet.flux[si[i]:si[i+1]]
        inds = np.where(f)[0]
        if len(inds):        
          t = t[inds]
          tdur = t[-1] - t[0]
          a = np.argmin(np.abs(time - (t[0] - tdur)))
          b = np.argmin(np.abs(time - (t[-1] + tdur)))
          planet._inds.append(list(range(a,b)))
  
  def plot_occultations(self, planet):
    '''
    
    '''
    
    # Get the occulted planet
    p = np.argmax(self._names == planet)
    planet = self.planets[p]
    
    # Loop over all events
    for t in planet._inds:
      
      # Set up the figure
      fig = pl.figure(figsize = (5, 6))
      fig.subplots_adjust(left = 0.175)
  
      # Plot three different wavelengths (first, mid, and last)
      axlc = pl.subplot2grid((5, 3), (3, 0), colspan = 3, rowspan = 2)
      axlc.plot(planet.time[t], planet.flux[t, 0], 'b-')
      axlc.plot(planet.time[t], planet.flux[t, planet.flux.shape[-1] // 2], 'g-')
      axlc.plot(planet.time[t], planet.flux[t, -1], 'r-')
      axlc.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
      axlc.set_ylabel('Occulted Flux [W/sr]', fontweight = 'bold', fontsize = 10)
      axlc.get_yaxis().set_major_locator(MaxNLocator(4))
      axlc.get_xaxis().set_major_locator(MaxNLocator(4))
      for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
        tick.set_fontsize(8)
    
      # Get the times of ingress, midpoint, and egress
      tstart = t[0] + np.argmax(planet.flux[t,0] < 0)
      tend = t[0] + len(planet.time[t]) - np.argmax(planet.flux[t,0][::-1] < 0)
      tmid = (tstart + tend) // 2
      o = planet.occultor[tmid]
      occultor = self.planets[o]

      # Plot the orbits of the two planets and all interior ones
      axxz = pl.subplot2grid((5, 3), (0, 0), colspan = 3, rowspan = 2)
      for j, _ in enumerate(self.planets):
        if (j > p) and (j > o):
          break
        if j == p:
          style = dict(color = 'r', alpha = 1, ls = '-', lw = 1)
        elif j == o:
          style = dict(color = 'k', alpha = 1, ls = '-', lw = 1)
        else:
          style = dict(color = 'k', alpha = 0.1, ls = '--', lw = 1)
        tmax = np.argmin(np.abs(self.planets[j].time - (self.planets[j].time[0] + self.planets[j].per)))
        axxz.plot(self.planets[j].x[:tmax], self.planets[j].z[:tmax], **style)
      
      # Plot their current positions
      axxz.plot(planet.x[tmid], planet.z[tmid], 'o', color = 'r', alpha = 1, markeredgecolor = 'k', zorder = 99)
      axxz.plot(occultor.x[tmid], occultor.z[tmid], 'o', color = 'lightgrey', alpha = 1, markeredgecolor = 'k', zorder = 99)
      
      # Which object is moving faster in the sky plane?
      xp0 = planet.x[tmid]
      xo0 = occultor.x[tmid]
      xp1 = planet.x[tmid+1]
      xo1 = occultor.x[tmid+1]
      if np.abs(xp1 - xp0) > np.abs(xo1 - xo0):
        t0 = np.argmin(np.abs(planet.time - (planet.time[tmid] - planet.per / 8)))
        ti = np.array(sorted(list(set(np.array(np.linspace(t0, tmid, 30), dtype = int)))))
        axxz.plot(planet.x[ti], planet.z[ti], 'o', color = 'r', alpha = float(j) / 30.)
      else:
        t0 = np.argmin(np.abs(occultor.time - (occultor.time[tmid] - occultor.per / 8)))
        ti = np.array(sorted(list(set(np.array(np.linspace(t0, tmid, 30), dtype = int)))))
        axxz.plot(occultor.x[ti], occultor.z[ti], 'o', color = 'lightgrey', alpha = float(j) / 30.)
      
      # Appearance
      axxz.set_aspect('equal')
      axxz.axis('off')
    
      # Plot the images
      axim = [None, None, None]
      axim[0] = pl.subplot2grid((5, 3), (2, 0), colspan = 1, rowspan = 1)
      axim[1] = pl.subplot2grid((5, 3), (2, 1), colspan = 1, rowspan = 1)
      axim[2] = pl.subplot2grid((5, 3), (2, 2), colspan = 1, rowspan = 1)
      for j in range(3): 
        axim[j].axis('off')
        axim[j].set_aspect('equal')
  
      # Recall that in quadrants II and III we mirror the geometry,
      # so flip the plots around
      if planet.x[tmid] < 0:
        self.image(tend, planet, occultor, ax = axim[0])
        self.image(tmid, planet, occultor, ax = axim[1])
        self.image(tstart, planet, occultor, ax = axim[2])
      else:
        self.image(tstart, planet, occultor, ax = axim[0])
        self.image(tmid, planet, occultor, ax = axim[1])
        self.image(tend, planet, occultor, ax = axim[2])
  
      # The title
      axxz.annotate(planet.name, xy = (0.15, 1.1),
                    xycoords = "axes fraction", ha = 'right', va = 'center',
                    fontweight = 'bold', color = 'r', fontsize = 12)
      axxz.annotate(occultor.name, xy = (0.85, 1.1),
                    xycoords = "axes fraction", ha = 'left', va = 'center',
                    fontweight = 'bold', color = 'grey', fontsize = 12)
      axxz.annotate("occulted by", xy = (0.5, 1.1),
                    xycoords = "axes fraction", ha = 'center', va = 'center',
                    fontweight = 'bold', fontsize = 12)
    
    pl.show()
  
  def plot_orbits(self, t, ax = None):
    '''
    
    '''

    # Set up the figure
    if ax is None:
      fig, ax = pl.subplots(1, figsize = (5, 6))
      fig.subplots_adjust(left = 0.175)

    # Plot the orbits of the two planets and all interior ones
    for j, _ in enumerate(self.planets):
    
      # The full orbit
      tmax = np.argmin(np.abs(self.planets[j].time - (self.planets[j].time[0] + self.planets[j].per)))
      x = self.planets[j].x[:tmax]
      z = self.planets[j].z[:tmax]
      
      # Thin
      thin = max(1, len(x) // 100)
      x = np.append(x[::thin], x[-1])
      z = np.append(z[::thin], z[-1])
      
      # Plot
      for i in range(len(x) - 1):
        ax.plot(x[i:i+2], z[i:i+2], '-', lw = 1, color = greys(i / (len(x) - 1.)))
      
      # The current position
      ax.plot(self.planets[j].x[t], self.planets[j].z[t], 'o', color = plasma(1 - self.planets[j].per / self.planets[-1].per), alpha = 1, markeredgecolor = 'k', zorder = 99)
    
    # Appearance
    ax.set_aspect('equal')
    ax.axis('off')
  
    return ax
      
  def image(self, t, occulted, occultor, ax = None, pad = 2.5, **kwargs):
    '''
    Plots an image of the `occulted` planet and the `occultor` at a given `time`.
  
    '''
  
    # Set up the plot
    if ax is None:
      fig, ax = pl.subplots(1, figsize = (6,6))
  
    # Plot the occultor
    r = occultor.r
    x = np.linspace(-r, r, 1000)
    y = np.sqrt(r ** 2 - x ** 2)
    ax.plot(x, y, color = 'k', zorder = 99, lw = 1)
    ax.plot(x, -y, color = 'k', zorder = 99, lw = 1)
    ax.fill_between(x, -y, y, color = 'lightgray', zorder = 99, lw = 1)

    # Plot the occulted planet
    r = occulted.r
    x0 = occulted.x[t] - occultor.x[t]
    y0 = occulted.y[t] - occultor.y[t]
    theta = np.arctan(occulted.z[t] / np.abs(occulted.x[t]))
    x = np.linspace(x0 - r, x0 + r, 1000)
    y = np.sqrt(r ** 2 - (x - x0) ** 2)
    ax.plot(x, y0 + y, color = 'k', zorder = 98, lw = 1)
    ax.plot(x, y0 - y, color = 'k', zorder = 98, lw = 1)
  
    # Plot the latitude ellipses
    for lat in np.linspace(0, np.pi, occulted.nl + 2)[1:-1]:
      a = occulted.r * np.abs(np.sin(lat))
      b = a * np.abs(np.sin(theta))
      x0 = occulted.x[t] - occulted.r * np.cos(lat) * np.cos(theta)
      y0 = occulted.y[t]
      xlimb = occulted.r * np.cos(lat) * np.sin(theta) * np.tan(theta)
      if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
        xmin = x0 - b
      else:
        xmin = x0 - xlimb
      if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
        xmax = x0 + b
      else:
        xmax = x0 - xlimb
      x = np.linspace(x0 - b, x0 + b, 1000)
      if theta > 0:
        x[x < x0 - xlimb] = np.nan
      else:
        x[x > x0 - xlimb] = np.nan
      A = b ** 2 - (x - x0) ** 2
      A[A<0] = 0
      y = (a / b) * np.sqrt(A)
      if np.abs(np.cos(lat)) < 1e-5:
        style = dict(color = 'k', ls = '--', lw = 1)
      else:
        style = dict(color = rdbu(0.5 * (np.cos(lat) + 1)), ls = '-', lw = 1)
      ax.plot(x - occultor.x[t], occulted.y[t] - occultor.y[t] + y, **style)
      ax.plot(x - occultor.x[t], occulted.y[t] - occultor.y[t] - y, **style)

    # Limits
    ax.set_xlim(occulted.x[t] - occultor.x[t] - (pad + 1) * r, occulted.x[t] - occultor.x[t] + (pad + 1) * r)
    ax.set_ylim(occulted.y[t] - occultor.y[t] - (pad + 1) * r, occulted.y[t] - occultor.y[t] + (pad + 1) * r)
  
    # In quadrants II and III we mirror the problem to compute it!
    if occulted.x[t] < 0:
      ax.invert_xaxis()
  
    return ax

# Define the planets
b = Planet('b', per = 1.51087081, inc = 89.65, a = 0.01111, r = 1.086, trn0 = 7322.51736, nlat = 11)
c = Planet('c', per = 2.4218233, inc = 89.67, a = 0.01521, r = 1.056, trn0 = 7282.80728, nlat = 11)
d = Planet('d', per = 4.049610, inc = 89.75, a = 0.02144, r = 0.772, trn0 = 7670.14165, nlat = 11)
e = Planet('e', per = 6.099615, inc = 89.86, a = 0.02817, r = 0.918, trn0 = 7660.37859, nlat = 11)
f = Planet('f', per = 9.206690, inc = 89.68, a = 0.0371, r = 1.045, trn0 = 7671.39767, nlat = 11)
g = Planet('g', per = 12.35294, inc = 89.71, a = 0.0451, r = 1.127, trn0 = 7665.34937, nlat = 11)
h = Planet('h', per = 18.766, inc = 89.80, a = 0.06, r = 0.755, trn0 = 7662.55463, nlat = 11)
system = System(b, c, d, e, f, g, h)

# Get the light curves
time = np.linspace(0, 10, 100000)
system.compute(time)

system.plot_occultations('c')
pl.show()