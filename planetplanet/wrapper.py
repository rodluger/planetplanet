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
cmap = pl.get_cmap('RdBu_r')

# Define constants
AUREARTH = 23454.9271
SEARTH = 1.361e3
MDFAST = 0
NEWTON = 1

# Load the library
try:
  ppo = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ppo.so'))
except:
  raise Exception("Can't find `ppo.so`; please run `make` to compile it.")

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
  :param int nlat: Number of latitude slices. Default `11`
  
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
              ("nlat", ctypes.c_int),
              ("x", ctypes.c_double),
              ("y", ctypes.c_double),
              ("z", ctypes.c_double)]
        
  def __init__(self, name, **kwargs):
    self.name = name
    self.per = kwargs.pop('per', 1.51087081)
    self.inc = kwargs.pop('inc', 89.65) * np.pi / 180.
    self.ecc = kwargs.pop('ecc', 0.)
    self.w = kwargs.pop('w', 0.) * np.pi / 180.
    self.a = kwargs.pop('a', 0.01111) * AUREARTH
    self.t0 = kwargs.pop('t0', 7322.51736)
    self.r = kwargs.pop('r', 1.086)
    self.albedo = kwargs.pop('albedo', 0.3)
    self.irrad = kwargs.pop('irrad', 4.25) * SEARTH
    self.nlat = kwargs.pop('nlat', 11)
    self.x = np.nan
    self.y = np.nan
    self.z = np.nan
    
class Settings(ctypes.Structure):
  '''
  The class that contains the model settings.
  
  :param float keptol: Kepler solver tolerance. Default `1.e-15`
  :param int maxkepiter: Maximum number of Kepler solver iterations. Default `100`
  :param str kepsolver: Kepler solver (`newton` | `mdfast`). Default `newton`
  :param float polyeps1: Tolerance in the polynomial root-finding routine. Default `2.0e-6`
  :param float polyeps2: Tolerance in the polynomial root-finding routine. Default `6.0e-9`
  :param int maxpolyiter: Maximum number of root finding iterations. Default `100`
  :param bool phasecurve: Compute the full phase curve of the planets? Default `False`
  
  '''
  
  _fields_ = [("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("kepsolver", ctypes.c_int),
              ("polyeps1", ctypes.c_double),
              ("polyeps2", ctypes.c_double),
              ("maxpolyiter", ctypes.c_int),
              ("phasecurve", ctypes.c_int)]
  
  def __init__(self, **kwargs):
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = eval(kwargs.pop('kepsolver', 'newton').upper())
    self.polyeps1 = kwargs.pop('polyeps1', 2.0e-6)
    self.polyeps2 = kwargs.pop('polyeps2', 6.0e-9)
    self.maxpolyiter = kwargs.pop('maxpolyiter', 100)
    self.phasecurve = int(kwargs.pop('phasecurve', False))

class Occultation(object):
  '''
  
  '''
  
  def __init__(self, planets, x, y, z, occulted, occultor, time, flux, settings):
    '''
    
    '''
    
    self.planets = planets
    self.x = x
    self.y = y
    self.z = z
    self.occulted = occulted
    self.occultor = occultor
    self.time = time
    self.flux = flux
    self.settings = settings
  
  def plot(self):
    '''
    
    '''
    
    # Set up the figure
    fig = pl.figure(figsize = (5, 6))
    fig.subplots_adjust(left = 0.175)
    
    # Plot three different wavelengths (first, mid, and last)
    axlc = pl.subplot2grid((5, 3), (3, 0), colspan = 3, rowspan = 2)
    axlc.plot(self.time, self.flux[:, 0], 'b-')
    axlc.plot(self.time, self.flux[:, self.flux.shape[-1] // 2], 'g-')
    axlc.plot(self.time, self.flux[:, -1], 'r-')
    axlc.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    axlc.set_ylabel('Occulted Flux [W/sr]', fontweight = 'bold', fontsize = 10)
    axlc.get_yaxis().set_major_locator(MaxNLocator(4))
    axlc.get_xaxis().set_major_locator(MaxNLocator(4))
    for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
      tick.set_fontsize(8)
      
    # Get the times of ingress, midpoint, and egress
    tstart = self.time[np.argmax(self.flux[:,0] < 0)]
    tend = self.time[::-1][np.argmax(self.flux[:,0][::-1] < 0)]
    tmid = (tstart + tend) / 2.
    dt = tend - tstart
    
    # Plot the orbits of the two planets and all interior ones
    axxz = pl.subplot2grid((5, 3), (0, 0), colspan = 3, rowspan = 2)
    for j, _ in enumerate(self.planets):
      if (j > self.occulted) and (j > self.occultor):
        break
      if j == self.occulted:
        style = dict(color = 'r', alpha = 1, ls = '-', lw = 1)
      elif j == self.occultor:
        style = dict(color = 'k', alpha = 1, ls = '-', lw = 1)
      else:
        style = dict(color = 'k', alpha = 0.1, ls = '--', lw = 1)
      axxz.plot(self.x[j], self.z[j], **style)

    # Plot their current positions
    Orbit = ppo.OrbitXYZ
    Orbit.restype = ctypes.c_int
    Orbit.argtypes = [ctypes.c_double, ctypes.POINTER(Planet), Settings]
    Orbit(tmid, ctypes.byref(self.planets[self.occulted]), self.settings)
    Orbit(tmid, ctypes.byref(self.planets[self.occultor]), self.settings)
    axxz.plot(self.planets[self.occulted].x, self.planets[self.occulted].z, 'o', color = 'r', alpha = 1, markeredgecolor = 'k', zorder = 99)
    axxz.plot(self.planets[self.occultor].x, self.planets[self.occultor].z, 'o', color = 'lightgrey', alpha = 1, markeredgecolor = 'k', zorder = 99)
    
    # Which object is moving faster in the sky plane?
    xp0 = self.planets[self.occulted].x
    xo0 = self.planets[self.occultor].x
    Orbit(tend, ctypes.byref(self.planets[self.occulted]), self.settings)
    Orbit(tend, ctypes.byref(self.planets[self.occulted]), self.settings)
    xp1 = self.planets[self.occulted].x
    xo1 = self.planets[self.occultor].x
    if np.abs(xp1 - xp0) > np.abs(xo1 - xo0):
      
      # Occulted planet is moving past occultor
      for j, t in enumerate(np.linspace(tmid - 0.1 * self.planets[self.occulted].per, tmid, 25)):
        Orbit(t, ctypes.byref(self.planets[self.occulted]), self.settings)
        axxz.plot(self.planets[self.occulted].x, self.planets[self.occulted].z, 'o', color = 'r', alpha = float(j) / 250.)
      
    else:
    
      # Occultor is moving past occulted planet
      for j, t in enumerate(np.linspace(tmid - 0.1 * self.planets[self.occultor].per, tmid, 50)):
        Orbit(t, ctypes.byref(self.planets[self.occultor]), self.settings)
        axxz.plot(self.planets[self.occultor].x, self.planets[self.occultor].z, 'o', color = 'lightgrey', alpha = float(j) / 250.)
    
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
    if self.planets[self.occulted].x < 0:
      Image(tend, self.planets[self.occulted], self.planets[self.occultor], ax = axim[0])
      Image(tmid, self.planets[self.occulted], self.planets[self.occultor], ax = axim[1])
      Image(tstart, self.planets[self.occulted], self.planets[self.occultor], ax = axim[2])
    else:
      Image(tstart, self.planets[self.occulted], self.planets[self.occultor], ax = axim[0])
      Image(tmid, self.planets[self.occulted], self.planets[self.occultor], ax = axim[1])
      Image(tend, self.planets[self.occulted], self.planets[self.occultor], ax = axim[2])
    
    # The title
    axxz.annotate(self.planets[self.occulted].name, xy = (0.15, 1.1),
                  xycoords = "axes fraction", ha = 'right', va = 'center',
                  fontweight = 'bold', color = 'r', fontsize = 12)
    axxz.annotate(self.planets[self.occultor].name, xy = (0.85, 1.1),
                  xycoords = "axes fraction", ha = 'left', va = 'center',
                  fontweight = 'bold', color = 'grey', fontsize = 12)
    axxz.annotate("occulted by", xy = (0.5, 1.1),
                  xycoords = "axes fraction", ha = 'center', va = 'center',
                  fontweight = 'bold', fontsize = 12)
    
    return fig, [axxz, axim, axlc]
    
class Lightcurve(object):
  '''
  Computes the occultation lightcurve for a system of `planets` at the
  time array `time` and in the wavelength range (`lambda1`, `lambda2`)
  at resolution `R`.
  
  '''
  
  def __init__(self, time, planets, lambda1 = 5, lambda2 = 15, R = 100, **kwargs):
    '''
    
    '''
    
    # Make planets into a list if necessary
    if type(planets) is Planet:
      self.planets = [planet]
    else:
      self.planets = planets
    self.time = time
    
    # Compute the wavelength grid
    wav = [lambda1]
    while(wav[-1] < lambda2):
      wav.append(wav[-1] + wav[-1] / R) 
    self.wavelength = np.array(wav) * 1e-6
    
    # Dimensions
    n = len(self.planets)
    nt = len(self.time)
    nwav = len(self.wavelength)
  
    # Initialize the C interface
    Orbit = ppo.OrbitXYZ
    Orbit.restype = ctypes.c_int
    Orbit.argtypes = [ctypes.c_double, ctypes.POINTER(Planet), Settings]
    Flux = ppo.Flux
    Flux.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.POINTER(Planet), 
                     Settings, ctypes.ARRAY(ctypes.c_double, nwav), 
                     ctypes.ARRAY(ctypes.c_int, n),
                     ctypes.ARRAY(ctypes.ARRAY(ctypes.c_double, nwav), n)]
    p = (Planet * n)(*self.planets)
    l = ctypes.ARRAY(ctypes.c_double, nwav)(*self.wavelength)
    f = ctypes.ARRAY(ctypes.ARRAY(ctypes.c_double, nwav), n)()
    o = ctypes.ARRAY(ctypes.c_int, n)()
    s = Settings(**kwargs)
    
    # Compute the orbits
    self.x = np.zeros((n, 50))
    self.y = np.zeros((n, 50))
    self.z = np.zeros((n, 50))
    for i in range(n):
      for j, t in enumerate(np.linspace(0, self.planets[i].per, 50)):
        Orbit(t, ctypes.byref(self.planets[i]), s)
        self.x[i,j] = self.planets[i].x
        self.y[i,j] = self.planets[i].y
        self.z[i,j] = self.planets[i].z

    # Call the C light curve function in a loop
    self._flux = np.zeros((n, nt, nwav))
    self._occultor = np.zeros((n, nt), dtype = 'int32')
    for i, t in enumerate(self.time):
      Flux(t, n, nwav, p, s, l, o, f)
      self._flux[:,i] = np.array(f)
      self._occultor[:,i] = np.array(o)
  
    # Loop over all planets and store each occultation
    # event as a separate object
    self.occultations = []
    for n in range(len(self.planets)):
    
      # Identify the different events
      inds = np.where(self._flux[n,:,0])[0]
      si = np.concatenate(([0], inds[np.where(np.diff(inds) > 1)] + 1, [len(self.time)]))
      ne = len(si) - 1
      
      # Loop over the events
      for i in range(ne):
    
        # Split the light curve and trim it
        t = self.time[si[i]:si[i+1]]
        f = self._flux[n, si[i]:si[i+1]]
        o = self._occultor[n, si[i]:si[i+1]]
        inds = np.where(f)[0]
        if len(inds):        
          
          t = t[inds]
          f = f[inds,:]
          o = o[inds[0]]

          # Add padding to ingress and egress
          pad = np.zeros_like(f)
          f = np.concatenate((pad, f, pad))
          tstart = t[0]
          tmid = (t[-1] + t[0]) / 2
          tend = t[-1]
          dt = 0.5 * (tend - tstart)
          lpad = np.linspace(tstart - dt, tstart, len(t))
          rpad = np.linspace(tend, tend + dt, len(t))
          t = np.concatenate((lpad, t, rpad))
          
          # Save the event
          self.occultations.append(Occultation(self.planets, self.x, self.y, self.z, n, o, t, f, s))
          
  def lightcurve(self, planet = 'all'):
    '''
    Returns a light curve of dimensions (`ntime`, `nwav`)
    evaluated at times `self.time` and wavelengths `self.wavelength`
    for a given `planet` name.
    
    '''
    
    # Return the planet(s) flux:
    if planet == 'all':
      return np.sum(self._flux, axis = 0)
    else:
      n = np.argmax([p.name == planet for p in self.planets])
      return self._flux[n]
  
  def plot(self, planet = 'all'):
    '''
    
    '''
    
    f = self.lightcurve()
    
    fig, ax = pl.subplots(1, figsize = (12, 4))
    ax.plot(self.time, f[:, 0], 'b-')
    ax.plot(self.time, f[:, f.shape[-1] // 2], 'g-')
    ax.plot(self.time, f[:, -1], 'r-')
    ax.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    ax.set_ylabel('Occulted Flux [W/sr]', fontweight = 'bold', fontsize = 10)
    ax.get_yaxis().set_major_locator(MaxNLocator(4))
    ax.get_xaxis().set_major_locator(MaxNLocator(4))
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
      tick.set_fontsize(8)
    
    return fig, ax
    
  def xyz(self, planet):
    '''
    Returns arrays corresponding to the `x`, `y` and `z` positions
    for a given `planet` name.
    
    '''
    
    # Get the planet index
    n = np.argmax([p.name == planet for p in self.planets])
    return self.x[n], self.y[n], self.z[n]

def Image(time, occulted, occultor, ax = None, pad = 2.5, **kwargs):
  '''
  Plots an image of the `occulted` planet and the `occultor` at a given `time`.
  
  '''
  
  # Set up the plot
  if ax is None:
    fig, ax = pl.subplots(1, figsize = (6,6))
  
  # Set up the orbit function
  OrbitXYZ = ppo.OrbitXYZ
  OrbitXYZ.restype = ctypes.c_int
  OrbitXYZ.argtypes = [ctypes.c_double, ctypes.POINTER(Planet), Settings]
  settings = Settings(**kwargs)
  OrbitXYZ(time, ctypes.byref(occultor), settings)
  OrbitXYZ(time, ctypes.byref(occulted), settings)
  
  # Plot the occultor
  r = occultor.r
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x, y, color = 'k', zorder = 99, lw = 1)
  ax.plot(x, -y, color = 'k', zorder = 99, lw = 1)
  ax.fill_between(x, -y, y, color = 'lightgray', zorder = 99, lw = 1)

  # Plot the occulted planet
  r = occulted.r
  x0 = occulted.x - occultor.x
  y0 = occulted.y - occultor.y
  theta = np.arctan(occulted.z / np.abs(occulted.x))
  x = np.linspace(x0 - r, x0 + r, 1000)
  y = np.sqrt(r ** 2 - (x - x0) ** 2)
  ax.plot(x, y0 + y, color = 'k', zorder = 98, lw = 1)
  ax.plot(x, y0 - y, color = 'k', zorder = 98, lw = 1)
  
  # Plot the latitude ellipses
  for lat in np.linspace(0, np.pi, occulted.nlat + 2)[1:-1]:
    a = occulted.r * np.abs(np.sin(lat))
    b = a * np.abs(np.sin(theta))
    x0 = occulted.x - occulted.r * np.cos(lat) * np.cos(theta)
    y0 = occulted.y
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
      style = dict(color = cmap(0.5 * (np.cos(lat) + 1)), ls = '-', lw = 1)
    ax.plot(x - occultor.x, occulted.y - occultor.y + y, **style)
    ax.plot(x - occultor.x, occulted.y - occultor.y - y, **style)

  # Limits
  ax.set_xlim(occulted.x - occultor.x - (pad + 1) * r, occulted.x - occultor.x + (pad + 1) * r)
  ax.set_ylim(occulted.y - occultor.y - (pad + 1) * r, occulted.y - occultor.y + (pad + 1) * r)
  
  # In quadrants II and III we mirror the problem to compute it!
  if occulted.x < 0:
    ax.invert_xaxis()
  
  return ax


# Define the planets
b = Planet('b', per = 1.51087081, inc = 89.65, a = 0.01111, r = 1.086, t0 = 7322.51736, nlat = 11)
c = Planet('c', per = 2.4218233, inc = 89.67, a = 0.01521, r = 1.056, t0 = 7282.80728, nlat = 11)
d = Planet('d', per = 4.049610, inc = 89.75, a = 0.02144, r = 0.772, t0 = 7670.14165, nlat = 11)
e = Planet('e', per = 6.099615, inc = 89.86, a = 0.02817, r = 0.918, t0 = 7660.37859, nlat = 11)
f = Planet('f', per = 9.206690, inc = 89.68, a = 0.0371, r = 1.045, t0 = 7671.39767, nlat = 11)
g = Planet('g', per = 12.35294, inc = 89.71, a = 0.0451, r = 1.127, t0 = 7665.34937, nlat = 11)
h = Planet('h', per = 18.766, inc = 89.80, a = 0.06, r = 0.755, t0 = 7662.55463, nlat = 11)
planets = [b, c, d, e, f, g, h]

# Get the light curves
time = np.linspace(0, 5, 100000)
lc = Lightcurve(time, planets, phasecurve = False)

for o in lc.occultations:
  o.plot()
pl.show()