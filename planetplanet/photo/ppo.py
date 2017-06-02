#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`ppo.py` - Python interface to C
----------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import ctypes
import numpy as np
np.seterr(invalid = 'ignore')
import os, shutil
from numpy.ctypeslib import ndpointer, as_ctypes
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import corner
rdbu = pl.get_cmap('RdBu_r')
greys = pl.get_cmap('Greys')
plasma = pl.get_cmap('plasma')
AUREARTH = 23454.9271
MSUNMEARTH = 332968.308
RSUNREARTH = 109.045013
SEARTH = 1.361e3
MSUN = 1.988416e30
RSUN = 6.957e8
G = 6.67428e-11
MEARTH = 5.9722e24
REARTH = 6.3781e6
DAYSEC = 86400.
AUM = 1.49598e11
G = 6.67428e-11
MDFAST = 0
NEWTON = 1

__all__ = ['Settings', 'Star', 'Planet', 'Body', 'System']

# Load the library
try:
  cwd = os.getcwd()
  if cwd != os.path.dirname(os.path.abspath(__file__)):
    shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'librebound.so'), cwd)
  libppo = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libppo.so'))
except:
  import pdb; pdb.set_trace()
  raise Exception("Can't find `libppo.so`; please run `make` to compile it.")

class Settings(ctypes.Structure):
  '''
  The class that contains the model settings.
  
  :param bool ttvs: Allow for TTVs? Uses `REBOUND` N-body code to compute orbits. Default `False`
  :param float keptol: Kepler solver tolerance. Default `1.e-15`
  :param int maxkepiter: Maximum number of Kepler solver iterations. Default `100`
  :param str kepsolver: Kepler solver (`newton` | `mdfast`). Default `newton`
  :param float polyeps1: Tolerance in the polynomial root-finding routine. Default `2.0e-6`
  :param float polyeps2: Tolerance in the polynomial root-finding routine. Default `6.0e-9`
  :param int maxpolyiter: Maximum number of root finding iterations. Default `100`
  :param float dt: Maximum timestep in days for the N-body solver. Default `0.01`
  :param bool adaptive: Adaptive grid for limb-darkened bodies? Default `False`
  
  '''
  
  _fields_ = [("ttvs", ctypes.c_int),
              ("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("kepsolver", ctypes.c_int),
              ("polyeps1", ctypes.c_double),
              ("polyeps2", ctypes.c_double),
              ("maxpolyiter", ctypes.c_int),
              ("dt", ctypes.c_double),
              ("adaptive", ctypes.c_int)]
  
  def __init__(self, **kwargs):
    self.ttvs = int(kwargs.pop('ttvs', False))
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = eval(kwargs.pop('kepsolver', 'newton').upper())
    self.polyeps1 = kwargs.pop('polyeps1', 1.0e-8)  # was 2.0e-6
    self.polyeps2 = kwargs.pop('polyeps2', 1.0e-15) # was 6.0e-9
    self.maxpolyiter = kwargs.pop('maxpolyiter', 100)
    self.dt = kwargs.pop('dt', 0.01)
    self.adaptive = int(kwargs.pop('adaptive', False))

def Star(*args, **kwargs):
  '''
  
  '''
  
  # Effective temperature and limb darkening
  T = kwargs.get('T', 2559.)
  u = kwargs.get('u', np.array([1., -1.]))
  u *= T / u[0]
  
  # Number of layers
  nl = kwargs.get('nl', 31)
  
  kwargs.update(dict(m = kwargs.get('m', 0.0802) * MSUNMEARTH, 
                     r = kwargs.get('r', 0.117) * RSUNREARTH, 
                     per = 0., inc = 0., ecc = 0., w = 0., 
                     Omega = 0., a = 0., t0 = 0., irrad = 0.,
                     albedo = 0., phasecurve = False, u = u,
                     nl = nl))
                     
  return Body(*args, **kwargs)

def Planet(*args, **kwargs):
  '''
  
  '''
  
  return Body(*args, **kwargs)

class Body(ctypes.Structure):
  '''
  The class containing all the input planet/star parameters.
  
  :param float per: Body mass in Earth masses. Default `0.85`
  :param float per: Orbital period in days. Default `1.51087081`
  :param float inc: Orbital inclination in degrees. Default `89.65`
  :param float ecc: Orbital eccentricity. Default `0.`
  :param float w: Longitude of pericenter in degrees. `0.`
  :param float Omega: Longitude of ascending node in degrees. `0.`
  :param float t0: Time of first transit in days. Default `7322.51736`
  :param float r: Body radius in Earth radii. Default `1.086`
  :param float albedo: Body albedo. Default `0.3`
  :param float irrad: Stellar irradiation at the body's distance in units \
         of the solar constant (1370 W/m^2). Default `0.3`
  :param bool phasecurve: Compute the full phase curve? Default `False`
  :param int nl: Number of latitude slices. Default `11`
  
  '''
  
  _fields_ = [("m", ctypes.c_double),
              ("per", ctypes.c_double),
              ("inc", ctypes.c_double),
              ("ecc", ctypes.c_double),
              ("w", ctypes.c_double),
              ("Omega", ctypes.c_double),
              ("a", ctypes.c_double),
              ("t0", ctypes.c_double),
              ("r", ctypes.c_double),
              ("albedo", ctypes.c_double),
              ("irrad", ctypes.c_double),
              ("phasecurve", ctypes.c_int),
              ("nu", ctypes.c_int),
              ("nl", ctypes.c_int),
              ("nt", ctypes.c_int),
              ("nw", ctypes.c_int),
              ("_u", ctypes.POINTER(ctypes.c_double)),
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
    self.m = kwargs.pop('m', 0.85)
    self.per = kwargs.pop('per', 1.51087081)
    self.inc = kwargs.pop('inc', 89.65) * np.pi / 180.
    self.ecc = kwargs.pop('ecc', 0.)
    self.w = kwargs.pop('w', 0.) * np.pi / 180.
    self.Omega = kwargs.pop('Omega', 0.) * np.pi / 180.
    self.r = kwargs.pop('r', 1.086)
    self.albedo = kwargs.pop('albedo', 0.3)
    self.irrad = kwargs.pop('irrad', 4.25) * SEARTH
    self.phasecurve = int(kwargs.pop('phasecurve', False))
    self.nl = kwargs.pop('nl', 11)
    self.color = kwargs.pop('color', 'k')
    
    # C stuff
    self.nt = 0
    self.nw = 0
    self.u = kwargs.pop('u', np.array([], dtype = float))
    self.nu = len(self.u)
    
    # Semi-major axis computed in `System` class
    self.a = 0.
    
    # Python stuff
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
    ax.plot(self.time, 1e-12 * self.flux[:, 0], 'b-')
    ax.plot(self.time, 1e-12 * self.flux[:, self.nw // 2], 'g-')
    ax.plot(self.time, 1e-12 * self.flux[:, -1], 'r-')
    ax.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    ax.set_ylabel(r'Occulted Flux [TW/$\mathbf{\mu}$m/sr]', fontweight = 'bold', fontsize = 10)
    ax.get_yaxis().set_major_locator(MaxNLocator(4))
    ax.get_xaxis().set_major_locator(MaxNLocator(4))
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
      tick.set_fontsize(8)
    
    return fig, ax

class System(object):

  def __init__(self, *bodies, **kwargs):
    '''
    
    '''
    
    self.bodies = bodies
    self.star = self.bodies[0]
    self.colors = [b.color for b in self.bodies]
    self.settings = Settings(**kwargs)
    self._names = np.array([p.name for p in self.bodies])
  
    # Compute the semi-major axis for each planet (in Earth radii)
    for body in self.bodies:
      body.a = ((body.per * DAYSEC) ** 2 * G * (self.star.m + body.m) * MEARTH / (4 * np.pi ** 2)) ** (1. / 3.) / REARTH
      
  def scatter_plot(self, tstart, tend):
    '''
    
    '''
  
    # Dimensions
    n = len(self.bodies)
    time = np.arange(tstart, tend, self.settings.dt)
    nt = len(time)
    nw = 1

    # Initialize the C interface
    Orbits = libppo.Orbits
    Orbits.restype = ctypes.c_int
    Orbits.argtypes = [ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nt),
                       ctypes.c_int, ctypes.POINTER(ctypes.POINTER(Body)),
                       Settings]
  
    # Allocate memory for all the arrays
    for body in self.bodies:
      body._u = np.ctypeslib.as_ctypes(body.u)
      body.time = np.zeros(nt)
      body._time = np.ctypeslib.as_ctypes(body.time)
      body.wavelength = np.zeros(nw)
      body._wavelength = np.ctypeslib.as_ctypes(body.wavelength)
      body.x = np.zeros(nt)
      body._x = np.ctypeslib.as_ctypes(body.x)
      body.y = np.zeros(nt)
      body._y = np.ctypeslib.as_ctypes(body.y)
      body.z = np.zeros(nt)
      body._z = np.ctypeslib.as_ctypes(body.z)
      body.occultor = np.zeros(nt, dtype = 'int32')
      body._occultor = np.ctypeslib.as_ctypes(body.occultor)
      # HACK: The fact that flux is 2d is a nightmare for ctypes. We will
      # treat it as a 1d array within C and keep track of the row/column
      # indices by hand...
      body.flux = np.zeros((nt, nw))
      body._flux1d = body.flux.reshape(nt * nw)
      body._flux = np.ctypeslib.as_ctypes(body._flux1d)

    # A pointer to a pointer to `Body`. This is an array of `n` `Body` instances, 
    # passed by reference. The contents can all be accessed through `bodies`
    ptr_bodies = (ctypes.POINTER(Body) * n)(*[ctypes.pointer(p) for p in self.bodies])

    # Call the light curve routine
    Orbits(nt, np.ctypeslib.as_ctypes(time), n, ptr_bodies, self.settings)
    
    # Loop over all bodies and plot each occultation event as a circle
    figp, axp = pl.subplots(1, figsize = (8,8))
    figp.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.1, top = 0.9)
    axp.axis('off')
    for body in self.bodies:
  
      # Identify the different events
      inds = np.where(body.occultor >= 0)[0]
      difs = np.where(np.diff(inds) > 1)[0]
    
      # Plot the orbit outline
      f = np.linspace(0, 2 * np.pi, 1000)
      r = body.a * (1 - body.ecc ** 2) / (1 + body.ecc * np.cos(f))
      x = r * np.cos(body.w + f) - r * np.sin(body.w + f) * np.cos(body.inc) * np.sin(body.Omega)
      z = r * np.sin(body.w + f) * np.sin(body.inc)
      axp.plot(x, z, 'k-', lw = 1, alpha = 0.05)
      n = np.argmin(1e10 * (x < 0) + np.abs(z))
      axp.annotate(body.name, xy = (x[n], z[n]), color = 'k', alpha = 0.2, fontweight = 'bold',
                   fontsize = 8, zorder = -99, ha = 'center', va = 'center')
      
      # Loop over individual ones
      plot_secondary = True
      for i in inds[difs]:
        
        # The occultor index
        occ = body.occultor[i]
        
        # Note that `i` is the last index of the occultation
        duration = np.argmax(body.occultor[:i][::-1] != body.occultor[i])          
        
        # Compute the minimum impact parameter
        impact = np.min(np.sqrt((self.bodies[occ].x[i-duration:i+1] - body.x[i-duration:i+1]) ** 2 + (self.bodies[occ].y[i-duration:i+1] - body.y[i-duration:i+1]) ** 2)) / (self.bodies[occ].r + body.r)
        
        # Transparency proportional to the impact parameter
        alpha = 1 - impact
        
        # Size = duration in minutes / 3
        ms = duration * self.settings.dt * 1440 / 3
        
        # If the occultor is the star, plot it only once
        if (body.occultor[i] == 0):
          if plot_secondary:
            axp.plot(body.x[i], body.z[i], 'o', color = self.colors[body.occultor[i]], alpha = alpha, ms = ms, markeredgecolor = 'none')
            plot_secondary = False
        else:
          axp.plot(body.x[i], body.z[i], 'o', color = self.colors[body.occultor[i]], alpha = alpha, ms = ms, markeredgecolor = 'none')
    
    # Legend 1: Occultor names/colors
    axl1 = pl.axes([0.025, 0.775, 0.2, 0.2])
    axl1.axis('off')
    axl1.set_xlim(-0.5, 1.5)
    axl1.set_ylim(-len(self.bodies) // 2 - 1, 1.5)
    axl1.annotate('Occultations by', xy = (0.5, 1), ha = 'center', va = 'center', fontweight = 'bold')
    for j, body in enumerate(self.bodies):
      if j < len(self.bodies) // 2:
        x, y = (0, -j)
      else:
        x, y = (0.825, len(self.bodies) // 2 - j)
      axl1.plot(x, y, 'o', color = self.colors[j], ms = 6, alpha = 1, markeredgecolor = 'none')
      axl1.annotate(body.name, xy = (x + 0.1, y), xycoords = 'data', 
                    ha = 'left', va = 'center', color = self.colors[j])
    
    # Legend 2: Size/duration
    axl2 = pl.axes([0.775, 0.775, 0.2, 0.2])
    axl2.axis('off')
    axl2.set_xlim(-1, 1)
    axl2.set_ylim(-3, 1.5)
    axl2.annotate('Duration', xy = (0., 1), ha = 'center', va = 'center', fontweight = 'bold')
    for j, duration in enumerate([10, 30, 60]):
      ms = duration / 3.
      axl2.plot(-0.65, -0.75 * j + 0.2, 'o', color = 'k', ms = ms, alpha = 0.65, markeredgecolor = 'none')
      axl2.annotate('%d minutes' % duration, xy = (-0.3, -0.75 * j + 0.2), xycoords = 'data', 
                    ha = 'left', va = 'center', color = 'k')

    # Legend 3: Transparency/impact parameter
    axl3 = pl.axes([0.025, 0.025, 0.2, 0.2])
    axl3.axis('off')
    axl3.set_xlim(-0.5, 1.5)
    axl3.set_ylim(-2, 1.5)
    axl3.annotate('Impact parameter', xy = (0.5, 0.65), ha = 'center', va = 'center', fontweight = 'bold')
    for j, impact in enumerate([0, 0.2, 0.4]):
      alpha = 1 - impact
      axl3.plot(-0.15, -0.75 * j, 'o', color = 'k', ms = 8, alpha = alpha, markeredgecolor = 'none')
      axl3.annotate('%.1f' % impact, xy = (-0.05, -0.75 * j), xycoords = 'data', 
                    ha = 'left', va = 'center', color = 'k')
    for j, impact in enumerate([0.6, 0.8, 0.99]):
      alpha = 1 - impact
      axl3.plot(0.675, -0.75 * j, 'o', color = 'k', ms = 8, alpha = alpha, markeredgecolor = 'none')
      axl3.annotate('%.1f' % impact, xy = (0.775, -0.75 * j), xycoords = 'data', 
                    ha = 'left', va = 'center', color = 'k')                
    
    # Observer direction
    axp.annotate("To observer", xy = (0.5, -0.1), xycoords = "axes fraction", xytext = (0, 30),
                 ha = 'center', va = 'center', annotation_clip = False, color = 'cornflowerblue',
                 textcoords = "offset points", arrowprops=dict(arrowstyle = "-|>", color = 'cornflowerblue'))
  
    return figp
      
  def corner_plot(self, tstart, tend):
    '''
    
    '''
  
    # Dimensions
    n = len(self.bodies)
    time = np.arange(tstart, tend, self.settings.dt)
    nt = len(time)
    nw = 1

    # Initialize the C interface
    Orbits = libppo.Orbits
    Orbits.restype = ctypes.c_int
    Orbits.argtypes = [ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nt),
                       ctypes.c_int, ctypes.POINTER(ctypes.POINTER(Body)),
                       Settings]
  
    # Allocate memory for all the arrays
    for body in self.bodies:
      body._u = np.ctypeslib.as_ctypes(body.u)
      body.time = np.zeros(nt)
      body._time = np.ctypeslib.as_ctypes(body.time)
      body.wavelength = np.zeros(nw)
      body._wavelength = np.ctypeslib.as_ctypes(body.wavelength)
      body.x = np.zeros(nt)
      body._x = np.ctypeslib.as_ctypes(body.x)
      body.y = np.zeros(nt)
      body._y = np.ctypeslib.as_ctypes(body.y)
      body.z = np.zeros(nt)
      body._z = np.ctypeslib.as_ctypes(body.z)
      body.occultor = np.zeros(nt, dtype = 'int32')
      body._occultor = np.ctypeslib.as_ctypes(body.occultor)
      # HACK: The fact that flux is 2d is a nightmare for ctypes. We will
      # treat it as a 1d array within C and keep track of the row/column
      # indices by hand...
      body.flux = np.zeros((nt, nw))
      body._flux1d = body.flux.reshape(nt * nw)
      body._flux = np.ctypeslib.as_ctypes(body._flux1d)

    # A pointer to a pointer to `Body`. This is an array of `n` `Body` instances, 
    # passed by reference. The contents can all be accessed through `bodies`
    ptr_bodies = (ctypes.POINTER(Body) * n)(*[ctypes.pointer(p) for p in self.bodies])

    # Call the light curve routine
    Orbits(nt, np.ctypeslib.as_ctypes(time), n, ptr_bodies, self.settings)

    # A corner plot for each planet, showing
    # distribution of phases, impact parameters, and durations
    fig = [None for i in self.bodies[1:]]
    for k, body in enumerate(self.bodies[1:]):
      
      # Identify the different planet-planet events
      inds = np.where(body.occultor > 0)[0]
      difs = np.where(np.diff(inds) > 1)[0]
          
      # Loop over individual ones
      duration = np.zeros(len(difs), dtype = int)
      phase = np.zeros(len(difs))
      impact = np.zeros(len(difs))
      for j, i in enumerate(inds[difs]):
        occ = body.occultor[i]
        duration[j] = np.argmax(body.occultor[:i][::-1] != occ) 
        phase[j] = np.arctan2(body.z[i], body.x[i]) * 180 / np.pi
        impact[j] = np.min(np.sqrt((self.bodies[occ].x[i-duration[j]:i+1] - body.x[i-duration[j]:i+1]) ** 2 + (self.bodies[occ].y[i-duration[j]:i+1] - body.y[i-duration[j]:i+1]) ** 2)) / (self.bodies[occ].r + body.r)

      samples = np.array([phase, impact, np.log10(duration * self.settings.dt * 1440)]).T
      fig[k] = corner.corner(samples, range = [(-180,180), (0,1), (0, 3)])
    
    return fig
        
  def compute(self, time, lambda1 = 5, lambda2 = 15, R = 100):
    '''
    
    '''
    
    # Compute the wavelength grid
    wav = [lambda1]
    while(wav[-1] < lambda2):
      wav.append(wav[-1] + wav[-1] / R) 
    wavelength = np.array(wav) * 1e-6
  
    # Dimensions
    n = len(self.bodies)
    nt = len(time)
    nw = len(wavelength)

    # Initialize the C interface
    Flux = libppo.Flux
    Flux.restype = ctypes.c_int
    Flux.argtypes = [ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nt),
                     ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nw),
                     ctypes.c_int, ctypes.POINTER(ctypes.POINTER(Body)),
                     Settings]
  
    # Allocate memory for all the arrays
    for body in self.bodies:
      body._u = np.ctypeslib.as_ctypes(body.u)
      body.time = np.zeros(nt)
      body._time = np.ctypeslib.as_ctypes(body.time)
      body.wavelength = np.zeros(nw)
      body._wavelength = np.ctypeslib.as_ctypes(body.wavelength)
      body.x = np.zeros(nt)
      body._x = np.ctypeslib.as_ctypes(body.x)
      body.y = np.zeros(nt)
      body._y = np.ctypeslib.as_ctypes(body.y)
      body.z = np.zeros(nt)
      body._z = np.ctypeslib.as_ctypes(body.z)
      body.occultor = np.zeros(nt, dtype = 'int32')
      body._occultor = np.ctypeslib.as_ctypes(body.occultor)
      # HACK: The fact that flux is 2d is a nightmare for ctypes. We will
      # treat it as a 1d array within C and keep track of the row/column
      # indices by hand...
      body.flux = np.zeros((nt, nw))
      body._flux1d = body.flux.reshape(nt * nw)
      body._flux = np.ctypeslib.as_ctypes(body._flux1d)

    # A pointer to a pointer to `Body`. This is an array of `n` `Body` instances, 
    # passed by reference. The contents can all be accessed through `bodies`
    ptr_bodies = (ctypes.POINTER(Body) * n)(*[ctypes.pointer(p) for p in self.bodies])

    # Call the light curve routine
    Flux(nt, np.ctypeslib.as_ctypes(time), nw, np.ctypeslib.as_ctypes(wavelength), n, ptr_bodies, self.settings)
  
    # Loop over all bodies and store each occultation event as a separate attribute
    for body in self.bodies:
    
      # Set the flag
      body._computed = True
    
      # Identify the different events
      inds = np.where(body.occultor > 0)[0]
      si = np.concatenate(([0], inds[np.where(np.diff(inds) > 1)] + 1, [nt]))

      # Loop over the events
      for i in range(len(si) - 1):
  
        # Split the light curve, trim it, and add a little padding
        t = time[si[i]:si[i+1]]
        f = body.flux[si[i]:si[i+1]]
        o = body.occultor[si[i]:si[i+1]]
        inds = np.where(o > 0)[0]
        if len(inds):        
          t = t[inds]
          tdur = t[-1] - t[0]
          a = np.argmin(np.abs(time - (t[0] - tdur)))
          b = np.argmin(np.abs(time - (t[-1] + tdur)))
          if b > a:
            body._inds.append(list(range(a,b)))
  
  def plot_occultations(self, body):
    '''
    
    '''
    
    # Get the occulted body
    p = np.argmax(self._names == body)
    body = self.bodies[p]
    
    # Figure lists
    fig = [None for i in range(len(body._inds))]
    axlc = [None for i in range(len(body._inds))]
    axxz = [None for i in range(len(body._inds))]
    axim = [[None, None, None] for i in range(len(body._inds))]
    
    # Stellar flux (baseline)
    normb = np.nanmedian(self.bodies[0].flux[:,0])
    normg = np.nanmedian(self.bodies[0].flux[:,body.flux.shape[-1] // 2])
    normr = np.nanmedian(self.bodies[0].flux[:,-1])
    
    # Loop over all events
    for i, t in enumerate(body._inds):
      
      # Set up the figure
      fig[i] = pl.figure(figsize = (5, 6))
      fig[i].subplots_adjust(left = 0.175)
  
      # Plot three different wavelengths (first, mid, and last)
      axlc[i] = pl.subplot2grid((5, 3), (3, 0), colspan = 3, rowspan = 2)
      axlc[i].plot(body.time[t], (normb + body.flux[t, 0]) / normb, 'b-')
      axlc[i].plot(body.time[t], (normg + body.flux[t, body.flux.shape[-1] // 2]) / normg, 'g-')
      axlc[i].plot(body.time[t], (normr + body.flux[t, -1]) / normr, 'r-')
      axlc[i].set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
      axlc[i].set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 10)
      axlc[i].get_yaxis().set_major_locator(MaxNLocator(4))
      axlc[i].get_xaxis().set_major_locator(MaxNLocator(4))
      for tick in axlc[i].get_xticklabels() + axlc[i].get_yticklabels():
        tick.set_fontsize(8)
    
      # Get the times of ingress, midpoint, and egress
      tstart = t[0] + np.argmax(body.occultor[t] >= 0)
      tend = t[0] + len(body.time[t]) - np.argmax(body.occultor[t][::-1] >= 0)
      tmid = (tstart + tend) // 2
      occultors = []
      for b in range(len(self.bodies)):
        for ti in t:
          if (body.occultor[ti] & 2 ** b):
            occultors.append(b)
      occultors = list(set(occultors))

      # Plot the orbits of all interior bodies
      axxz[i] = pl.subplot2grid((5, 3), (0, 0), colspan = 3, rowspan = 2)
      for j, _ in enumerate(self.bodies):
        if (j > p) and (j > occultors[-1]):
          break
        if j == p:
          style = dict(color = 'r', alpha = 1, ls = '-', lw = 1)
        elif j in occultors:
          style = dict(color = 'k', alpha = 1, ls = '-', lw = 1)
        else:
          style = dict(color = 'k', alpha = 0.1, ls = '--', lw = 1)
        tmax = np.argmin(np.abs(self.bodies[j].time - (self.bodies[j].time[0] + self.bodies[j].per)))
        axxz[i].plot(self.bodies[j].x[:tmax], self.bodies[j].z[:tmax], **style)
      
      # Plot their current positions
      axxz[i].plot(body.x[tmid], body.z[tmid], 'o', color = 'r', alpha = 1, markeredgecolor = 'k', zorder = 99)
      for occultor in [self.bodies[o] for o in occultors]:
        axxz[i].plot(occultor.x[tmid], occultor.z[tmid], 'o', color = 'lightgrey', alpha = 1, markeredgecolor = 'k', zorder = 99)
      
      # Appearance
      axxz[i].set_ylim(-max(np.abs(axxz[i].get_ylim())), max(np.abs(axxz[i].get_ylim())))
      axxz[i].set_xlim(-max(np.abs(axxz[i].get_xlim())), max(np.abs(axxz[i].get_xlim())))
      axxz[i].set_aspect('equal')
      axxz[i].axis('off')
      
      # Plot the other planets
      for b in self.bodies:
        axxz[i].plot(b.x[tmid], b.z[tmid], 'o', color = '#cccccc', alpha = 1, markeredgecolor = 'none', zorder = 98)
      
      # Plot orbital motion
      t0 = np.argmin(np.abs(body.time - (body.time[tmid] - body.per / 8)))
      ti = np.array(sorted(list(set(np.array(np.linspace(t0, tmid, 30), dtype = int)))))
      axxz[i].plot(body.x[ti], body.z[ti], 'o', color = 'r', alpha = float(j) / 30.)
      for occultor in [self.bodies[o] for o in occultors]:
        t0 = np.argmin(np.abs(occultor.time - (occultor.time[tmid] - occultor.per / 8)))
        ti = np.array(sorted(list(set(np.array(np.linspace(t0, tmid, 30), dtype = int)))))
        axxz[i].plot(occultor.x[ti], occultor.z[ti], 'o', color = 'lightgrey', alpha = float(j) / 30.)

      # Plot the images
      axim[i][0] = pl.subplot2grid((5, 3), (2, 0), colspan = 1, rowspan = 1)
      axim[i][1] = pl.subplot2grid((5, 3), (2, 1), colspan = 1, rowspan = 1)
      axim[i][2] = pl.subplot2grid((5, 3), (2, 2), colspan = 1, rowspan = 1)
      for j in range(3): 
        axim[i][j].axis('off')
        axim[i][j].set_aspect('equal')
      self.image(tstart, body, occultors, ax = axim[i][0])
      self.image(tmid, body, occultors, ax = axim[i][1])
      self.image(tend, body, occultors, ax = axim[i][2])
  
      # The title
      if len(occultors) == 1:
        axxz[i].annotate("%s occulted by %s" % (body.name, self.bodies[occultors[0]].name), xy = (0.5, 1.25),
                         xycoords = "axes fraction", ha = 'center', va = 'center',
                         fontweight = 'bold', fontsize = 12)
      else:
        axxz[i].annotate("%s occulted by %s" % (body.name, 
                        ", ".join([occultor.name for occultor in [self.bodies[o] for o in occultors]])), 
                         xy = (0.5, 1.25),
                         xycoords = "axes fraction", ha = 'center', va = 'center',
                         fontweight = 'bold', fontsize = 12)
      axxz[i].annotate("Duration: %.2f minutes" % ((body.time[tend] - body.time[tstart]) * 1440.),
                       xy = (0.5, 1.1), ha = 'center', va = 'center', xycoords = 'axes fraction',
                       fontsize = 10, style = 'italic')
      
    return fig, axlc, axxz, axim
  
  def plot_orbits(self, t, ax = None):
    '''
    
    '''

    # Set up the figure
    if ax is None:
      fig, ax = pl.subplots(1, figsize = (5, 6))
      fig.subplots_adjust(left = 0.175)

    # Plot the orbits of the two bodies and all interior ones
    for j, _ in enumerate(self.bodies):
    
      # The full orbit
      tmax = np.argmin(np.abs(self.bodies[j].time - (self.bodies[j].time[0] + self.bodies[j].per)))
      x = self.bodies[j].x[:tmax]
      z = self.bodies[j].z[:tmax]
      
      # Thin
      thin = max(1, len(x) // 100)
      x = np.append(x[::thin], x[-1])
      z = np.append(z[::thin], z[-1])
      
      # Plot
      for i in range(len(x) - 1):
        ax.plot(x[i:i+2], z[i:i+2], '-', lw = 1, color = greys(i / (len(x) - 1.)))
      
      # The current position
      ax.plot(self.bodies[j].x[t], self.bodies[j].z[t], 'o', color = plasma(1 - self.bodies[j].per / self.bodies[-1].per), alpha = 1, markeredgecolor = 'k', zorder = 99)
    
    # Appearance
    ax.set_aspect('equal')
    ax.axis('off')
  
    return ax
      
  def image(self, t, occulted, occultors, ax = None, pad = 2.5, **kwargs):
    '''
    Plots an image of the `occulted` body and the `occultor` at a given `time`.
  
    '''
  
    # Set up the plot
    if ax is None:
      fig, ax = pl.subplots(1, figsize = (6,6))
  
    # Plot the occultor
    for occultor in [self.bodies[o] for o in occultors]: 
      r = occultor.r
      x = np.linspace(occultor.x[t] - r, occultor.x[t] + r, 1000)
      y = np.sqrt(r ** 2 - (x - occultor.x[t]) ** 2)
      ax.plot(x, occultor.y[t] + y, color = 'k', zorder = 99, lw = 1)
      ax.plot(x, occultor.y[t] - y, color = 'k', zorder = 99, lw = 1)
      ax.fill_between(x, occultor.y[t] - y, occultor.y[t] + y, color = 'lightgray', zorder = 99, lw = 1)

    # Plot the occulted body
    r = occulted.r
    x0 = occulted.x[t]
    y0 = occulted.y[t]
    if occulted.nu == 0:
      theta = np.arctan(occulted.z[t] / np.abs(occulted.x[t]))
    else:
      theta = np.pi / 2
    x = np.linspace(-r, r, 1000)
    y = np.sqrt(r ** 2 - x ** 2)
    ax.plot(x, y, color = 'k', zorder = 98, lw = 1)
    ax.plot(x, -y, color = 'k', zorder = 98, lw = 1)
  
    # Plot the latitude ellipses
    for lat in np.linspace(0, np.pi, occulted.nl + 2)[1:-1]:
      a = occulted.r * np.abs(np.sin(lat))
      b = a * np.abs(np.sin(theta))
      x0 = -occulted.r * np.cos(lat) * np.cos(theta)
      y0 = 0
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
      ax.plot(x, y, **style)
      ax.plot(x, -y, **style)

    return ax

  def plot_lightcurve(self, wavelength = 15.):
    '''
    
    '''
    
    # Plot
    fig, ax = pl.subplots(1, figsize = (12, 4))
    time = self.bodies[0].time
    assert (wavelength >= self.bodies[0].wavelength[0]) and (wavelength >= self.bodies[0].wavelength[-1]), "Wavelength value outside of computed grid."
    w = np.argmax(1e-6 * wavelength <= self.bodies[0].wavelength)
    flux = np.sum([b.flux[:,w] for b in self.bodies], axis = 0)
    flux /= np.nanmedian(flux)
    ax.plot(time, flux, 'k-', lw = 1)
    
    # Appearance
    ax.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    ax.set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 10)
    ax.get_yaxis().set_major_locator(MaxNLocator(4))
    ax.get_xaxis().set_major_locator(MaxNLocator(4))
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
      tick.set_fontsize(8)
    
    # Limits
    ymax = np.nanmax(flux)
    ymin = np.nanmin(flux)
    yrng = ymax - ymin
    ax.set_ylim(ymin - 0.2 * yrng, ymax + 0.2 * yrng)
    
    # Label all of the events
    for body in self.bodies:
      for i, t in enumerate(body._inds):
        tstart = t[0] + np.argmax(body.occultor[t] > 0)
        tend = t[0] + len(body.time[t]) - np.argmax(body.occultor[t][::-1] > 0)
        tmid = (tstart + tend) // 2
        occultors = []
        for b in range(len(self.bodies)):
          for ti in t:
            if (body.occultor[ti] & 2 ** b):
              occultors.append(b)
        occultors = list(set(occultors))
        time = body.time[tmid]
        for n, occultor in enumerate([self.bodies[o] for o in occultors]):
          ax.annotate("%s" % body.name, xy = (time, ymax + (0.1 + 0.05 * n) * yrng), ha = 'center',
                      va = 'center', color = occultor.color, fontweight = 'bold',
                      fontsize = 8)
    
    return fig, ax