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
cmap = pl.get_cmap('RdBu_r')
AUREARTH = 23454.9271
MDFAST = 0
NEWTON = 1
try:
  ppo = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ppo.so'))
except:
  raise Exception("Can't find `ppo.so`; please run `make` to compile it.")

class Planet(ctypes.Structure):
  '''
  The class containing all the input planet parameters
  
  '''
  
  _fields_ = [("per", ctypes.c_double),
              ("inc", ctypes.c_double),
              ("ecc", ctypes.c_double),
              ("w", ctypes.c_double),
              ("a", ctypes.c_double),
              ("t0", ctypes.c_double),
              ("r", ctypes.c_double),
              ("noon", ctypes.c_double),
              ("midnight", ctypes.c_double),
              ("nlat", ctypes.c_int),
              ("x", ctypes.c_double),
              ("y", ctypes.c_double),
              ("z", ctypes.c_double)]
        
  def __init__(self, name, **kwargs):
    self.name = name
    self.per = kwargs.pop('per', 1.51087081)
    inc = kwargs.pop('inc', 89.65)
    self.inc = inc * np.pi / 180.
    self.ecc = kwargs.pop('ecc', 0.)
    self.w = kwargs.pop('w', 0.)
    self.a = kwargs.pop('a', 0.01111) * AUREARTH
    self.t0 = kwargs.pop('t0', 7322.51736)
    self.r = kwargs.pop('r', 1.086)
    self.noon = kwargs.pop('noon', 1.)
    self.midnight = kwargs.pop('midnight', 0.1)
    self.nlat = kwargs.pop('nlat', 11)
    self.x = np.nan
    self.y = np.nan
    self.z = np.nan
    
class Settings(ctypes.Structure):
  '''
  The class that contains the model settings
  
  '''
  
  _fields_ = [("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("kepsolver", ctypes.c_int),
              ("phasecurve", ctypes.c_int)]
  
  def __init__(self, **kwargs):
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = kwargs.pop('kepsolver', NEWTON)
    self.phasecurve = int(kwargs.pop('phasecurve', False))
    
def Image(time, occultor, occulted, ax = None, pad = 2.5, **kwargs):
  '''
  
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
  
  #if occulted.x < 0:
  #  ax.invert_xaxis()
  
  return ax

def Lightcurve(time, planets, **kwargs):
  '''
  
  '''
  
  # Initialize
  Flux = ppo.Flux
  Flux.restype = ctypes.c_double
  Flux.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.POINTER(Planet), 
                   Settings, ctypes.ARRAY(ctypes.c_double, len(planets)),
                   ctypes.ARRAY(ctypes.c_int, len(planets))]
  if type(planets) is Planet:
    planets = [planet]
  n = len(planets)
  p = (Planet * n)(*planets)
  f = ctypes.ARRAY(ctypes.c_double, n)()
  o = ctypes.ARRAY(ctypes.c_int, n)()
  s = Settings(**kwargs)
  
  # Call the C function in a loop
  flux = np.zeros((len(time), n))
  occultor = np.zeros((len(time), n), dtype = 'int32')
  for i, t in enumerate(time):
    Flux(t, n, p, s, f, o)
    flux[i] = np.array(f)
    occultor[i] = np.array(o)
  return flux, occultor

b = Planet('b', per = 1.51087081, inc = 89.65, a = 0.01111, r = 1.086, t0 = 7322.51736, midnight = 0, nlat = 31)
c = Planet('c', per = 2.4218233, inc = 89.67, a = 0.01521, r = 1.056, t0 = 7282.80728, midnight = 0, nlat = 31)
d = Planet('d', per = 4.049610, inc = 89.75, a = 0.02144, r = 0.772, t0 = 7670.14165, midnight = 0, nlat = 31)
e = Planet('e', per = 6.099615, inc = 89.86, a = 0.02817, r = 0.918, t0 = 7660.37859, midnight = 0, nlat = 31)
f = Planet('f', per = 9.206690, inc = 89.68, a = 0.0371, r = 1.045, t0 = 7671.39767, midnight = 0, nlat = 31)
g = Planet('g', per = 12.35294, inc = 89.71, a = 0.0451, r = 1.127, t0 = 7665.34937, midnight = 0, nlat = 31)
h = Planet('h', per = 18.766, inc = 89.80, a = 0.06, r = 0.755, t0 = 7662.55463, midnight = 0, nlat = 31)
planets = [b, c, d, e, f, g, h]

# Get the light curves
time = np.linspace(0, 10, 100000)
fluxes, occultors = Lightcurve(time, planets, phasecurve = False)

# Loop over each planet
for n, _ in enumerate(planets):

  # Get the flux
  flux = fluxes[:,n]
  occultor = occultors[:,n]

  # Identify the different events
  inds = np.where(flux)[0]
  s = np.concatenate(([0], inds[np.where(np.diff(inds) > 1)] + 1, [len(flux)]))
  ne = len(s) - 1
  
  # Set up the figure
  fig = pl.figure(figsize = (10, 4))
  axlc = [None for i in range(ne)]
  axim1 = [None for i in range(ne)]
  axim2 = [None for i in range(ne)]
  axim3 = [None for i in range(ne)]

  # Loop over the events
  for i in range(ne):
    
    # Split the light curve and trim it
    t = time[s[i]:s[i+1]]
    f = flux[s[i]:s[i+1]]
    o = occultor[s[i]:s[i+1]]
    inds = np.where(f)[0]
    t = t[inds]
    f = f[inds]
    o = o[inds[0]]

    # Add padding to ingress and egress
    pad = np.zeros(len(f))
    f = np.concatenate((pad, f, pad))
    tstart = t[0]
    tmid = (t[-1] + t[0]) / 2
    tend = t[-1]
    dt = 0.5 * (tend - tstart)
    lpad = np.linspace(tstart - dt, tstart, len(t))
    rpad = np.linspace(tend, tend + dt, len(t))
    t = np.concatenate((lpad, t, rpad))
    
    # Plot the occultation
    axlc[i] = pl.subplot2grid((4, 3 * ne), (1, 3 * i), colspan = 3, rowspan = 3, sharey = axlc[0] if i > 0 else None)
    if i > 0:
      axlc[i].set_yticklabels([])
    axlc[i].plot(t, f)
    
    # Plot the image
    axim1[i] = pl.subplot2grid((4, 3 * ne), (0, 3 * i), colspan = 1, rowspan = 1)
    axim2[i] = pl.subplot2grid((4, 3 * ne), (0, 3 * i + 1), colspan = 1, rowspan = 1)
    axim3[i] = pl.subplot2grid((4, 3 * ne), (0, 3 * i + 2), colspan = 1, rowspan = 1)
    Image(tstart, planets[n], planets[o], ax = axim1[i])
    Image(tmid, planets[n], planets[o], ax = axim2[i])
    Image(tend, planets[n], planets[o], ax = axim3[i])
    
    axim1[i].annotate(planets[o].name, xy = (0,0), ha = 'center', va = 'center', zorder = 100, fontsize = 8)
    axim2[i].annotate(planets[o].name, xy = (0,0), ha = 'center', va = 'center', zorder = 100, fontsize = 8)
    axim3[i].annotate(planets[o].name, xy = (0,0), ha = 'center', va = 'center', zorder = 100, fontsize = 8)
    axim1[i].axis('off')
    axim2[i].axis('off')
    axim3[i].axis('off')
    
    print(planets[n].nlat)
    
  pl.show()
  quit()

fig, ax = pl.subplots(1, figsize = (12, 4))
ax.plot(time, flux[:,0], lw = 1)
pl.show()