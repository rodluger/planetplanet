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
import matplotlib.animation as animation
from matplotlib.ticker import MaxNLocator
from scipy.integrate import quad
from tqdm import tqdm
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
HPLANCK = 6.62607004e-34
CLIGHT = 2.998e8
KBOLTZ = 1.38064852e-23
MDFAST = 0
NEWTON = 1
MINUTE = 1. / 1440.

__all__ = ['Star', 'Planet', 'System']

# Load the library
libppo = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libppo.so'))

class Settings(ctypes.Structure):
  '''
  The class that contains the model settings. This class is used internally.
  
  :param bool nbody: Uses `REBOUND` N-body code to compute orbits. Default `False`
  :param float keptol: Kepler solver tolerance. Default `1.e-15`
  :param int maxkepiter: Maximum number of Kepler solver iterations. Default `100`
  :param str kepsolver: Kepler solver (`newton` | `mdfast`). Default `newton`
  :param float polyeps1: Tolerance in the polynomial root-finding routine. Default `1.0e-8`
  :param float polyeps2: Tolerance in the polynomial root-finding routine. Default `1.0e-15`
  :param int maxpolyiter: Maximum number of root finding iterations. Default `100`
  :param float timestep: Timestep in days for the N-body solver. Default `0.01`
  :param bool adaptive: Adaptive grid for limb-darkened bodies? Default `True`
  :param bool quiet: Suppress output? Default `False`
  :param float mintheta: Absolute value of the minimum phase angle in degrees. Below this \
         angle, elliptical boundaries of constant surface brightness on the planet surface are \
         treated as vertical lines. Default `1.`
  :param int maxvertices: Maximum number of vertices allowed in the area computation. Default `999`
  :param int maxfunctions: Maximum number of functions allowed in the area computation. Default `999`
  :param int oversample: Oversampling factor for each exposure. Default `1`
  :param float distance: Distance to the system in parsecs. Default `10.`
  
  '''
  
  _fields_ = [("_nbody", ctypes.c_int),
              ("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("_kepsolver", ctypes.c_int),
              ("polyeps1", ctypes.c_double),
              ("polyeps2", ctypes.c_double),
              ("maxpolyiter", ctypes.c_int),
              ("timestep", ctypes.c_double),
              ("_adaptive", ctypes.c_int),
              ("_quiet", ctypes.c_int),
              ("_mintheta", ctypes.c_double),
              ("maxvertices", ctypes.c_int),
              ("maxfunctions", ctypes.c_int),
              ("oversample", ctypes.c_int),
              ("distance", ctypes.c_double)]
  
  def __init__(self, **kwargs):
    self.nbody = kwargs.pop('nbody', False)
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = kwargs.pop('kepsolver', 'newton')
    self.polyeps1 = kwargs.pop('polyeps1', 1.0e-8)
    self.polyeps2 = kwargs.pop('polyeps2', 1.0e-15)
    self.maxpolyiter = kwargs.pop('maxpolyiter', 100)
    self.timestep = kwargs.pop('timestep', 0.01)
    self.adaptive = kwargs.pop('adaptive', True)
    self.quiet = kwargs.pop('quiet', False)
    self.mintheta = kwargs.pop('mintheta', 1.)
    self.maxvertices = kwargs.pop('maxvertices', 999)
    self.maxfunctions = kwargs.pop('maxfunctions', 999)
    self.oversample = max(1, kwargs.pop('oversample', 1))
    self.distance = kwargs.pop('distance', 10.)
  
  @property
  def params(self):
    return ['nbody', 'keptol', 'maxkepiter', 'kepsolver', 'polyeps1', 'polyeps2',
            'maxpolyiter', 'timestep', 'adaptive', 'quiet', 'mintheta', 'maxvertices',
            'maxfunctions', 'oversample', 'distance']
  
  @property
  def mintheta(self):
    return self._mintheta * 180 / np.pi
  
  @mintheta.setter
  def mintheta(self, val):
    self._mintheta = val * np.pi / 180.

  @property
  def nbody(self):
    return bool(self._nbody)
  
  @nbody.setter
  def nbody(self, val):
    self._nbody = int(val)

  @property
  def adaptive(self):
    return bool(self._adaptive)
  
  @adaptive.setter
  def adaptive(self, val):
    self._adaptive = int(val)

  @property
  def quiet(self):
    return bool(self._quiet)
  
  @quiet.setter
  def quiet(self, val):
    self._quiet = int(val)

  @property
  def kepsolver(self):
    return self._kepsolver
  
  @kepsolver.setter
  def kepsolver(self, val):
    self._kepsolver = eval(val.upper())

def Star(name, **kwargs):
  '''
  
  Returns a `Body` instance of type `star`.
  
  :param str name: A unique identifier for this star
  :param float m: Mass in solar masses. Default `1.`
  :param float r: Radius in solar radii. Default `1.`
  :param float per: Orbital period in days. Default `0.`
  :param float inc: Orbital inclination in degrees. Default `90.`
  :param float ecc: Orbital eccentricity. Default `0.`
  :param float w: Longitude of pericenter in degrees. `0.`
  :param float Omega: Longitude of ascending node in degrees. `0.`
  :param float t0: Time of primary eclipse in days. Default `0.`
  :param float teff: The effective temperature of the star in Kelvin. Default `5577.`
  :param array_like limbdark: The limb darkening coefficients. These are the coefficients \
         in the Taylor expansion of `(1 - mu)`, starting with the first order (linear) \
         coefficient, where `mu = cos(theta)` is the radial coordinate on the surface of \
         the star. Each coefficient may either be a scalar, in which case limb darkening is \
         assumed to be grey (the same at all wavelengths), or a callable whose single argument \
         is the wavelength array in microns. Default is `[1.0]`, a grey linear limb darkening law.
  :param int nz: Number of zenith angle slices. Default `31`
  :param str color: Object color (for plotting). Default `k`
  
  '''
  
  return Body(name, 'star', **kwargs)

def Planet(name, **kwargs):
  '''
  
  Returns a `Body` instance of type `planet`.
  
  :param str name: A unique identifier for this planet
  :param float m: Mass in Earth masses. Default `1.`
  :param float r: Radius in Earth radii. Default `1.`
  :param float per: Orbital period in days. Default `3.`
  :param float inc: Orbital inclination in degrees. Default `90.`
  :param float ecc: Orbital eccentricity. Default `0.`
  :param float w: Longitude of pericenter in degrees. `0.`
  :param float Omega: Longitude of ascending node in degrees. `0.`
  :param float t0: Time of primary eclipse in days. Default `0.`
  :param bool phasecurve: Compute the phasecurve for this planet? Default `False`
  :param bool airless: Treat this as an airless planet? If `True`, computes light curves \
         in the instant re-radiation limit, where the surface brightness is proportional \
         to the cosine of the zenith angle (the angle between the line connecting \
         the centers of the planet and the star and the line connecting the center of the \
         planet and a given point on its surface. A fixed nightside temperature may be specified \
         via the `tnight` kwarg. If `False`, treats the planet as a limb-darkened blackbody. \
         Default `True`
  :param float albedo: Planetary albedo (airless limit). Default `0.3`
  :param float tnight: Nightside temperature in Kelvin (airless limit). Default `40`
  :param array_like limbdark: The limb darkening coefficients (thick atmosphere limit). These are the coefficients \
         in the Taylor expansion of `(1 - mu)`, starting with the first order (linear) \
         coefficient, where `mu = cos(theta)` is the radial coordinate on the surface of \
         the star. Each coefficient may either be a scalar, in which case limb darkening is \
         assumed to be grey (the same at all wavelengths), or a callable whose single argument \
         is the wavelength in microns. Default is `[1.0]`, a grey linear limb darkening law.
  :param int nz: Number of zenith angle slices. Default `11`
  :param str color: Object color (for plotting). Default `r`
  
  '''
  
  return Body(name, 'planet', **kwargs)

class Body(ctypes.Structure):
  '''
  The class containing all the input planet/star parameters.

  '''
  
  _fields_ = [("_m", ctypes.c_double),
              ("per", ctypes.c_double),
              ("_inc", ctypes.c_double),
              ("ecc", ctypes.c_double),
              ("_w", ctypes.c_double),
              ("_Omega", ctypes.c_double),
              ("a", ctypes.c_double),
              ("t0", ctypes.c_double),
              ("_r", ctypes.c_double),
              ("albedo", ctypes.c_double),
              ("_teff", ctypes.c_double),
              ("tnight", ctypes.c_double),
              ("_phasecurve", ctypes.c_int),
              ("_blackbody", ctypes.c_int),
              ("nu", ctypes.c_int),
              ("nz", ctypes.c_int),
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
              
  def __init__(self, name, body_type, **kwargs):
    
    # Check
    self.name = name
    self.body_type = body_type
    assert body_type in ['planet', 'star'], "Argument `body_type` must be either `planet` or `star`."
    
    # User
    self.m = kwargs.pop('m', 1.)
    self.r = kwargs.pop('r', 1.)
    self.ecc = kwargs.pop('ecc', 0.)
    self.w = kwargs.pop('w', 0.)
    self.Omega = kwargs.pop('Omega', 0.)
    self.inc = kwargs.pop('inc', 90.)
    self.trn0 = kwargs.pop('trn0', 0.)
    
    # These defaults are different depending on body type
    if self.body_type == 'planet':
      self.airless = kwargs.pop('airless', True)
      self.nz = kwargs.pop('nz', 11)
      self.per = kwargs.pop('per', 3.)
      self.albedo = kwargs.pop('albedo', 0.3)
      self.teff = 0
      if self.airless:
        self.tnight = kwargs.pop('tnight', 40.)
        self.limbdark = []
      else:
        self.tnight = 0
        self.limbdark = kwargs.pop('limbdark', [])
      self.phasecurve = kwargs.pop('phasecurve', False)
      self.color = kwargs.pop('color', 'r')
    elif self.body_type == 'star':
      self.airless = False
      self.nz = kwargs.pop('nz', 31)
      self.per = kwargs.pop('per', 0.)
      self.albedo = 0.
      self.tnight = 0.
      self.teff = kwargs.pop('teff', 5577.)
      self.limbdark = kwargs.pop('limbdark', [1.])
      self.phasecurve = False
      self.color = kwargs.pop('color', 'k')

    # C stuff, computed in `System` class
    self.nt = 0
    self.nw = 0
    self.nu = 0
    self.t0 = 0.
    self.a = 0.
    
    # Python stuff
    self._inds = []
    self._computed = False
  
  @property
  def m(self):
    if self.body_type == 'planet':
      return self._m
    elif self.body_type == 'star':
      return self._m / MSUNMEARTH
      
  @m.setter
  def m(self, val):
    if self.body_type == 'planet':
      self._m = val
    elif self.body_type == 'star':
      self._m = val * MSUNMEARTH

  @property
  def r(self):
    if self.body_type == 'planet':
      return self._r
    elif self.body_type == 'star':
      return self._r / RSUNREARTH
      
  @r.setter
  def r(self, val):
    if self.body_type == 'planet':
      self._r = val
    elif self.body_type == 'star':
      self._r = val * RSUNREARTH  
      
  @property
  def limbdark(self):
    return self._limbdark
  
  @limbdark.setter
  def limbdark(self, val):
    assert hasattr(val, '__len__'), "Limb darkening coefficients must be provided as a list of scalars or as a list of functions."
    self._limbdark = val
    
  @property
  def phasecurve(self):
    return bool(self._phasecurve)
  
  @phasecurve.setter
  def phasecurve(self, val):
    self._phasecurve = int(val)

  @property
  def airless(self):
    return not bool(self._blackbody)
  
  @airless.setter
  def airless(self, val):
    self._blackbody = int(not bool(val))
  
  @property
  def inc(self):
    return self._inc * 180 / np.pi
  
  @inc.setter
  def inc(self, val):
    self._inc = val * np.pi / 180.

  @property
  def w(self):
    return self._w * 180 / np.pi
  
  @w.setter
  def w(self, val):
    self._w = val * np.pi / 180.

  @property
  def Omega(self):
    return self._Omega * 180 / np.pi
  
  @Omega.setter
  def Omega(self, val):
    self._Omega = val * np.pi / 180.
  
  @property
  def teff(self):
    return self._teff
  
  @teff.setter
  def teff(self, val):
    self._teff = val
    if self._teff == 0:
      self.u = []
  
class System(object):
  '''
  
  A planetary system class. Instantiate with all bodies in the system
  and the desired settings, passed as kwargs.
  
  ** Calculation settings **
  
  :param *bodies: Any number of `Planet` or `Star` instances \
         comprising all the bodies in the system. The first body is assumed to be the primary.
  :param bool nbody: Uses `REBOUND` N-body code to compute orbits. Default `False`
  :param float keptol: Kepler solver tolerance. Default `1.e-15`
  :param int maxkepiter: Maximum number of Kepler solver iterations. Default `100`
  :param str kepsolver: Kepler solver (`newton` | `mdfast`). Default `newton`
  :param float polyeps1: Tolerance in the polynomial root-finding routine. Default `1.0e-8`
  :param float polyeps2: Tolerance in the polynomial root-finding routine. Default `1.0e-15`
  :param int maxpolyiter: Maximum number of root finding iterations. Default `100`
  :param float timestep: Timestep in days for the N-body solver. Default `0.01`
  :param bool adaptive: Adaptive grid for limb-darkened bodies? Default `True`
  :param bool quiet: Suppress output? Default `False`
  :param float mintheta: Absolute value of the minimum phase angle in degrees. Below this \
         angle, elliptical boundaries of constant surface brightness on the planet surface are \
         treated as vertical lines. Default `1.`
  :param int maxvertices: Maximum number of vertices allowed in the area computation. Default `999`
  :param int maxfunctions: Maximum number of functions allowed in the area computation. Default `999`
  :param int oversample: Oversampling factor for each exposure. Default `1`
  :param float distance: Distance to the system in parsecs. Default `10.`

  '''
  
  def __init__(self, *bodies, **kwargs):
    '''

    '''
    
    self.bodies = bodies
    self.settings = Settings(**kwargs)
    self.reset()
    
  def reset(self):
    '''
    Resets the system and recomputes some orbital properties passed to the integrator.
    
    '''
    
    # Move params set by the user over to the settings class
    for param in self.settings.params:
      if hasattr(self, param):
        setattr(self.settings, param, getattr(self, param))
        delattr(self, param)
        
    # Make planets accessible as properties
    self.primary = self.bodies[0]
    for body in self.bodies:
      setattr(self, body.name, body)
    self._names = np.array([p.name for p in self.bodies])
    self.colors = [b.color for b in self.bodies]
    
    # Compute the semi-major axis for each planet (in Earth radii)
    # and the time of pericenter passage 
    for body in self.bodies:
    
      # From Kepler's law
      body.a = ((body.per * DAYSEC) ** 2 * G * (self.primary._m + body._m) * MEARTH / (4 * np.pi ** 2)) ** (1. / 3.) / REARTH
    
      # See, e.g., Shields et al. 2015
      fi = (3 * np.pi / 2.) - body._w
      tperi0 = (body.per * np.sqrt(1. - body.ecc * body.ecc) / (2. * np.pi) * (body.ecc * np.sin(fi) / 
               (1. + body.ecc * np.cos(fi)) - 2. / np.sqrt(1. - body.ecc * body.ecc) * 
               np.arctan2(np.sqrt(1. - body.ecc * body.ecc) * np.tan(fi/2.), 1. + body.ecc)))
    
      # We define the mean anomaly to be zero at t = t0 = trn0 + tperi0
      body.t0 = body.trn0 + tperi0
    
    # Reset animations
    self._animations = []
       
  def scatter_plot(self, tstart, tend, dt = 0.001):
    '''
    
    '''
  
    # Reset
    self.reset()
  
    # Dimensions
    n = len(self.bodies)
    time = np.arange(tstart, tend, dt)
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
      # HACK: Same for the limb darkening coefficients
      body._u1d = np.array([], dtype = float)
      body._u = np.ctypeslib.as_ctypes(body._u1d)
      # Dimensions
      body.nu = 0
      body.nt = nt
      body.nw = nw

    # A pointer to a pointer to `Body`. This is an array of `n` `Body` instances, 
    # passed by reference. The contents can all be accessed through `bodies`
    ptr_bodies = (ctypes.POINTER(Body) * n)(*[ctypes.pointer(p) for p in self.bodies])

    # Call the light curve routine
    Orbits(nt, np.ctypeslib.as_ctypes(time), n, ptr_bodies, self.settings)
    
    # Loop over all bodies and plot each occultation event as a circle
    figp, axp = pl.subplots(1, figsize = (8,8))
    figp.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.1, top = 0.9)
    axp.axis('off')
    for body in self.bodies[1:]:
  
      # Identify the different events
      inds = np.where(body.occultor > 0)[0]
      difs = np.where(np.diff(inds) > 1)[0]
    
      # Plot the orbit outline
      f = np.linspace(0, 2 * np.pi, 1000)
      r = body.a * (1 - body.ecc ** 2) / (1 + body.ecc * np.cos(f))
      x = r * np.cos(body._w + f) - r * np.sin(body._w + f) * np.cos(body._inc) * np.sin(body._Omega)
      z = r * np.sin(body._w + f) * np.sin(body._inc)
      axp.plot(x, z, 'k-', lw = 1, alpha = 0.05)
      n = np.argmin(1e10 * (x < 0) + np.abs(z))
      axp.annotate(body.name, xy = (x[n], z[n]), color = 'k', alpha = 0.2, fontweight = 'bold',
                   fontsize = 8, zorder = -99, ha = 'center', va = 'center')
      
      # Loop over individual ones
      plot_secondary = True
      for i in inds[difs]:
        
        # Loop over all possible occultors
        for occ in range(len(self.bodies)):
        
          # Is body `occ` occulting?
          if (body.occultor[i] & 2 ** occ):
            
            # Note that `i` is the last index of the occultation
            duration = np.argmax(body.occultor[:i][::-1] & 2 ** occ == 0)          

            # Compute the minimum impact parameter
            idx = range(i - duration, i + 1)
            impact = np.min(np.sqrt((self.bodies[occ].x[idx] - body.x[idx]) ** 2 + 
                                    (self.bodies[occ].y[idx] - body.y[idx]) ** 2)) / (self.bodies[occ]._r + body._r)

            # Transparency proportional to the impact parameter
            alpha = 0.8 * (1 - impact) + 0.01
        
            # Size = duration in minutes / 3
            ms = duration * dt * 1440 / 3
        
            # If the occultor is the star, plot it only once
            if (occ == 0):
              if plot_secondary:
                axp.plot(body.x[i], body.z[i], 'o', color = self.colors[occ], alpha = alpha, ms = ms, markeredgecolor = 'none')
                plot_secondary = False
            else:
              axp.plot(body.x[i], body.z[i], 'o', color = self.colors[occ], alpha = alpha, ms = ms, markeredgecolor = 'none')
          
        # Check for mutual transits
        if self.bodies[0].occultor[i]:
          
          # Get all bodies currently occulting the star
          occultors = []
          for occ in range(1, len(self.bodies)):
            if (self.bodies[0].occultor[i] & 2 ** occ):
              occultors.append(occ)
          
          # Check if any of these occult each other
          for occ1 in occultors:
            for occ2 in occultors:
              if self.bodies[occ1].occultor[i] & 2 ** occ2:
                axp.plot(self.bodies[occ1].x[i], self.bodies[occ1].z[i], 'x', 
                         color = self.colors[occ2], alpha = 1, zorder = 100, 
                         ms = 20)
                
          
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
      alpha = 0.8 * (1 - impact) + 0.01
      axl3.plot(-0.15, -0.75 * j, 'o', color = 'k', ms = 8, alpha = alpha, markeredgecolor = 'none')
      axl3.annotate('%.1f' % impact, xy = (-0.05, -0.75 * j), xycoords = 'data', 
                    ha = 'left', va = 'center', color = 'k')
    for j, impact in enumerate([0.6, 0.8, 1.]):
      alpha = 1 - impact
      axl3.plot(0.675, -0.75 * j, 'o', color = 'k', ms = 8, alpha = alpha, markeredgecolor = 'none')
      axl3.annotate('%.1f' % impact, xy = (0.775, -0.75 * j), xycoords = 'data', 
                    ha = 'left', va = 'center', color = 'k')                
    
    # Observer direction
    axp.annotate("To observer", xy = (0.5, -0.1), xycoords = "axes fraction", xytext = (0, 30),
                 ha = 'center', va = 'center', annotation_clip = False, color = 'cornflowerblue',
                 textcoords = "offset points", arrowprops=dict(arrowstyle = "-|>", color = 'cornflowerblue'))
  
    return figp
      
  def histogram(self, tstart, tend, dt = 0.0001):
    '''
    
    '''
  
    # Reset
    self.reset()
  
    # Dimensions
    n = len(self.bodies)
    time = np.arange(tstart, tend, dt)
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
      # HACK: Same for the limb darkening coefficients
      body._u1d = np.array([], dtype = float)
      body._u = np.ctypeslib.as_ctypes(body._u1d)
      # Dimensions
      body.nu = 0
      body.nt = nt
      body.nw = nw
      
    # A pointer to a pointer to `Body`. This is an array of `n` `Body` instances, 
    # passed by reference. The contents can all be accessed through `bodies`
    ptr_bodies = (ctypes.POINTER(Body) * n)(*[ctypes.pointer(p) for p in self.bodies])

    # Call the light curve routine
    Orbits(nt, np.ctypeslib.as_ctypes(time), n, ptr_bodies, self.settings)

    # A histogram of the distribution of phases, impact parameters, and durations
    hist = [[] for body in self.bodies[1:]]
    for k, body in enumerate(self.bodies[1:]):
      
      # Identify the different planet-planet events
      inds = np.where(body.occultor > 0)[0]
      difs = np.where(np.diff(inds) > 1)[0]
      
      # Loop over individual ones
      for i in inds[difs]:
        
        # Loop over possible occultors
        for occ in range(1, len(self.bodies)):
          
          # Is body `occ` occulting (but not behind the star)?
          if (body.occultor[i] & 2 ** occ) and (body.occultor[i] & 1 == 0):

            # Note that `i` is the last index of the occultation
            duration = np.argmax(body.occultor[:i][::-1] & 2 ** occ == 0)
            if duration > 0:
            
              # Orbital phase
              phase = np.arctan2(body.x[i], -body.z[i]) * 180 / np.pi
            
              # Compute the minimum impact parameter
              idx = range(i - duration, i + 1)
              impact = np.min(np.sqrt((self.bodies[occ].x[idx] - body.x[idx]) ** 2 + 
                                      (self.bodies[occ].y[idx] - body.y[idx]) ** 2)) / (self.bodies[occ]._r + body._r)
            
              # Convert duration to log
              duration = np.log10(duration * dt * 1440)
            
              # Running list
              hist[k].append((phase, impact, duration))
      
      # Make into array  
      hist[k] = np.array(hist[k])
    
    return hist
        
  def compute(self, time, lambda1 = 5, lambda2 = 15, R = 100):
    '''
    
    '''
    
    # Reset
    self.reset()

    # Compute the wavelength grid
    wav = [lambda1]
    while(wav[-1] < lambda2):
      wav.append(wav[-1] + wav[-1] / R) 
    wavelength = np.array(wav)
    
    # Compute all limb darkening coefficients
    for body in self.bodies:
      body.u = [None for ld in body.limbdark]
      for n, ld in enumerate(body.limbdark):
        if callable(ld):
          body.u[n] = ld(wavelength)
        elif not hasattr(ld, '__len__'):
          body.u[n] = ld * np.ones_like(wavelength)
        else:
          raise Exception("Limb darkening coefficients must be provided as a list of scalars or as a list of functions.")
      body.u = np.array(body.u)

    # Convert from microns to meters
    wavelength *= 1e-6
    
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
      # HACK: Same for the limb darkening coefficients
      body._u1d = body.u.reshape(-1)
      body._u = np.ctypeslib.as_ctypes(body._u1d)
      # Dimensions
      body.nu = len(body.u)
      body.nt = nt
      body.nw = nw

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
          a = np.argmin(np.abs(time - (t[0] - 0.25 * tdur)))
          b = np.argmin(np.abs(time - (t[-1] + 0.25 * tdur)))
          if b > a:
            body._inds.append(list(range(a,b)))
  
  def compute_orbits(self, time):
    '''
    
    '''
    
    # Reset
    self.reset()
    
    # Dimensions
    n = len(self.bodies)
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
      # HACK: Same for the limb darkening coefficients
      body._u1d = np.array([], dtype = float)
      body._u = np.ctypeslib.as_ctypes(body._u1d)
      # Dimensions
      body.nu = 0
      body.nt = nt
      body.nw = nw

    # A pointer to a pointer to `Body`. This is an array of `n` `Body` instances, 
    # passed by reference. The contents can all be accessed through `bodies`
    ptr_bodies = (ctypes.POINTER(Body) * n)(*[ctypes.pointer(p) for p in self.bodies])

    # Call the light curve routine
    Orbits(nt, np.ctypeslib.as_ctypes(time), n, ptr_bodies, self.settings)
  
  def next_occultation(self, tstart, occulted, min_duration = 10, max_impact = 0.5, occultor = None, maxruns = 100, dt = 0.001):
    '''
    
    '''
    
    # Quiet?
    quiet = self.settings.quiet
    self.settings.quiet = True
    if quiet:
      iterator = lambda x: range(x)
    else:
      iterator = lambda x: tqdm(range(x))
      
    # Convert occultors to indices if necessary
    if occultor is None:
      occultor = list(range(1, len(self.bodies)))
    elif not hasattr(occultor, '__len__'):
      occultor = [occultor]
    for i, occ in enumerate(occultor):
      if occ in self.bodies:
        occultor[i] = np.argmax([b == occ for b in self.bodies])

    # Loop until we find one
    for n in iterator(maxruns):
      
      # Compute the orbits, 100 days at a time
      time = np.arange(tstart + n * 100., tstart + (n + 1) * 100., dt)
      self.compute_orbits(time)
      
      # Identify the different planet-planet events
      inds = np.where(occulted.occultor > 0)[0]
      difs = np.where(np.diff(inds) > 1)[0]
      
      # Loop over individual ones
      for i in inds[difs]:
        
        # Loop over possible planet occultors
        for occ in occultor:
          
          # Is body `occ` occulting (but not behind the star)?
          if (occulted.occultor[i] & 2 ** occ) and ((occ == 0) or (occulted.occultor[i] & 1 == 0)):

            # Note that `i` is the last index of the occultation
            duration = np.argmax(occulted.occultor[:i][::-1] & 2 ** occ == 0)
            if duration * dt * 1440 >= min_duration:
            
              # Compute the minimum impact parameter
              idx = range(i - duration, i + 1)
              b = np.sqrt((self.bodies[occ].x[idx] - occulted.x[idx]) ** 2 + (self.bodies[occ].y[idx] - occulted.y[idx]) ** 2) / (self.bodies[occ]._r + occulted._r)
              ind = np.argmin(b)
              if b[ind] <= max_impact:
                if not quiet:
                  print("Next occultation of %s by %s is at t = %.2f days." % (occulted.name, self.bodies[occ].name, time[idx[ind]]))
                self.settings.quiet = quiet
                return time[idx[ind]]
    
    # Nothing found...
    if not quiet:
      print("No occultation %s by %s found." % (occulted.name, self.bodies[occ].name))
    self.settings.quiet = quiet
    return np.nan
    
  def observe(self):
    '''
    TODO: Jake
    
    '''
    
    from ..detect import jwst
    
    # Call Jake's code
    self.observation = 0
  
  def plot_occultation(self, body, time, interval = 50, gifname = None):
    '''
    
    '''
    
    if not self.settings.quiet:
      print("Plotting the occultation...")
    
    # Check file name
    if gifname is not None:
      if gifname.endswith(".gif"):
        gifname = gifname[:-4]    
    
    # Get the occulted body
    p = np.argmax(self._names == body)
    body = self.bodies[p]
    
    # Get the indices of the occultation
    tind = np.argmin(np.abs(body.time - time))
    iind = np.argmax([tind in inds for inds in body._inds])
    if (iind == 0) and not (tind in body._inds[0]):
      return None
    t = body._inds[iind]
    
    # Stellar flux (baseline)
    normb = np.nanmedian(self.bodies[0].flux[:,0])
    normg = np.nanmedian(self.bodies[0].flux[:,body.flux.shape[-1] // 2])
    normr = np.nanmedian(self.bodies[0].flux[:,-1])
    
    # Set up the figure
    fig = pl.figure(figsize = (7, 8))
    fig.subplots_adjust(left = 0.175)

    # Plot three different wavelengths (first, mid, and last)
    axlc = pl.subplot2grid((5, 3), (3, 0), colspan = 3, rowspan = 2)
    axlc.plot(body.time[t], (int(p > 0) * normb + body.flux[t, 0]) / normb, 'b-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * body.wavelength[0])) + r"\ \mu\mathrm{m}$")
    axlc.plot(body.time[t], (int(p > 0) * normg + body.flux[t, body.flux.shape[-1] // 2]) / normg, 'g-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * body.wavelength[body.flux.shape[-1] // 2])) + r"\ \mu\mathrm{m}$")
    axlc.plot(body.time[t], (int(p > 0) * normr + body.flux[t, -1]) / normr, 'r-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * body.wavelength[-1])) + r"\ \mu\mathrm{m}$")
    axlc.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    axlc.set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 10)
    axlc.get_yaxis().set_major_locator(MaxNLocator(4))
    axlc.get_xaxis().set_major_locator(MaxNLocator(4))
    tracker = axlc.axvline(body.time[t[0]], color = 'k', alpha = 0.5, lw = 1, ls = '--')
    for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
      tick.set_fontsize(8)
    axlc.legend(loc = 'lower right', fontsize = 8)
    axlc.ticklabel_format(useOffset = False)
    if body.time[t[0]] > 1e4:
      for label in axlc.get_xmajorticklabels():
        label.set_rotation(30)
    
    # Get the times of ingress, midpoint, and egress
    tstart = t[0] + np.argmax(body.occultor[t] > 0)
    tend = t[0] + len(body.time[t]) - np.argmax(body.occultor[t][::-1] > 0)
    tmid = (tstart + tend) // 2
    
    # Sort occultors by z-order (occultor closest to observer last)
    occultors = []
    for b in range(len(self.bodies)):
      for ti in t:
        if (body.occultor[ti] & 2 ** b):
          occultors.append(b)
    occultors = list(set(occultors))
    zorders = [-self.bodies[o].z[tmid] for o in occultors]
    occultors = [o for (z,o) in sorted(zip(zorders, occultors))]

    # Plot the orbits of all bodies
    axxz = pl.subplot2grid((5, 3), (0, 0), colspan = 3, rowspan = 2)
    f = np.linspace(0, 2 * np.pi, 1000)
    for j, b in enumerate(self.bodies):
      if j == p:
        style = dict(color = 'r', alpha = 1, ls = '-', lw = 1)
      elif j in occultors:
        style = dict(color = 'k', alpha = 1, ls = '-', lw = 1)
      else:
        style = dict(color = 'k', alpha = 0.1, ls = '--', lw = 1)
      r = b.a * (1 - b.ecc ** 2) / (1 + b.ecc * np.cos(f))
      x = r * np.cos(b._w + f) - r * np.sin(b._w + f) * np.cos(b._inc) * np.sin(b._Omega)
      z = r * np.sin(b._w + f) * np.sin(b._inc)
      axxz.plot(x, z, **style)

    # Plot the locations of the bodies
    ptb = [None for b in self.bodies]
    for bi, b in enumerate(self.bodies):
      if b == body:
        ptb[bi], = axxz.plot(b.x[tmid], b.z[tmid], 'o', color = 'r', alpha = 1, markeredgecolor = 'k', zorder = 99)
      elif bi in occultors:
        ptb[bi], = axxz.plot(b.x[tmid], b.z[tmid], 'o', color = 'lightgrey', alpha = 1, markeredgecolor = 'k', zorder = 99)
      else:
        ptb[bi], = axxz.plot(b.x[tmid], b.z[tmid], 'o', color = '#dddddd', alpha = 1, markeredgecolor = '#999999', zorder = 99)
      
    # Appearance
    axxz.set_ylim(-max(np.abs(axxz.get_ylim())), max(np.abs(axxz.get_ylim())))
    axxz.set_xlim(-max(np.abs(axxz.get_xlim())), max(np.abs(axxz.get_xlim())))
    axxz.set_aspect('equal')
    axxz.axis('off')

    # Plot the image
    axim = pl.subplot2grid((5, 3), (2, 0), colspan = 3, rowspan = 1) 
    _, pto = self.plot_image(tmid, body, occultors, ax = axim)
    xmin = min([self.bodies[o].x[tstart] - 3 * self.bodies[o]._r for o in occultors])
    xmax = max([self.bodies[o].x[tend] + 3 * self.bodies[o]._r for o in occultors])
    if xmin > xmax: xmin, xmax = xmax, xmin
    if (body.x[tmid] - xmin) > (xmax - body.x[tmid]):
      dx = body.x[tmid] - xmin
    else:
      dx = xmax - body.x[tmid]
    dx = max(dx, 1.5 * body._r)
    ymin = min([self.bodies[o].y[tstart] - 3 * self.bodies[o]._r for o in occultors])
    ymax = max([self.bodies[o].y[tend] + 3 * self.bodies[o]._r for o in occultors])
    if ymin > ymax: ymin, ymax = ymax, ymin
    if (body.y[tmid] - ymin) > (ymax - body.y[tmid]):
      dy = body.y[tmid] - ymin
    else:
      dy = ymax - body.y[tmid]
    dy = max(dy, 1.5 * body._r)
    axim.set_xlim(0 - dx, 0 + dx)
    axim.set_ylim(0 - dy, 0 + dy)
    axim.axis('off')
    axim.set_aspect('equal')
    
    # The title
    if len(occultors) == 1:
      axxz.annotate("%s occulted by %s" % (body.name, self.bodies[occultors[0]].name), xy = (0.5, 1.25),
                       xycoords = "axes fraction", ha = 'center', va = 'center',
                       fontweight = 'bold', fontsize = 12)
    else:
      axxz.annotate("%s occulted by %s" % (body.name, 
                      ", ".join([occultor.name for occultor in [self.bodies[o] for o in occultors]])), 
                       xy = (0.5, 1.25),
                       xycoords = "axes fraction", ha = 'center', va = 'center',
                       fontweight = 'bold', fontsize = 12)
    axxz.annotate("Duration: %.2f minutes" % ((body.time[tend] - body.time[tstart]) * 1440.),
                     xy = (0.5, 1.1), ha = 'center', va = 'center', xycoords = 'axes fraction',
                     fontsize = 10, style = 'italic')
    
    # Animate!
    if gifname is not None:
      tmp = '%s.%03d.gif' % (gifname, len(self._animations) + 1)
    else:
      tmp = None
    self._animations.append(Animation(t, fig, axim, tracker, pto, ptb, body, 
                            self.bodies, [self.bodies[o] for o in occultors],
                            interval = interval, gifname = tmp, quiet = self.settings.quiet))

    return fig, axlc, axxz, axim
      
  def plot_image(self, t, occulted, occultors = None, ax = None, pad = 2.5, occultor_alpha = 1, **kwargs):
    '''
    Plots an image of the `occulted` body and its occultors at a given index of the time array `t`.
  
    '''
  
    # Set up the plot
    if ax is None:
      fig, ax = pl.subplots(1, figsize = (6,6))
    
    # Get the occultors
    if occultors is None:
      occultors = []
      for b in range(len(self.bodies)):
        if (occulted.occultor[t] & 2 ** b):
          occultors.append(b)
      occultors = list(set(occultors))
    
    # Plot the occulted body
    r = occulted._r
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
  
    # Plot the zenith angle ellipses
    for za in np.linspace(0, np.pi, occulted.nz + 2)[1:-1]:
      a = occulted._r * np.abs(np.sin(za))
      b = a * np.abs(np.sin(theta))
      xE = -occulted._r * np.cos(za) * np.cos(theta)
      yE = 0
      xlimb = occulted._r * np.cos(za) * np.sin(theta) * np.tan(theta)
      if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
        xmin = xE - b
      else:
        xmin = xE - xlimb
      if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
        xmax = xE + b
      else:
        xmax = xE - xlimb
      x = np.linspace(xE - b, xE + b, 1000)
      if theta > 0:
        x[x < xE - xlimb] = np.nan
      else:
        x[x > xE - xlimb] = np.nan
      A = b ** 2 - (x - xE) ** 2
      A[A<0] = 0
      y = (a / b) * np.sqrt(A)
      if np.abs(np.cos(za)) < 1e-5:
        style = dict(color = 'k', ls = '--', lw = 1, alpha = 0.5)
      else:
        style = dict(color = rdbu(0.5 * (np.cos(za) + 1)), ls = '-', lw = 1, alpha = 0.5)
      if (occulted.x[t] < 0):
        ax.plot(-x, y, **style)
        ax.plot(-x, -y, **style)
      else:
        ax.plot(x, y, **style)
        ax.plot(x, -y, **style)
    
    # Plot the occultors
    pto = [None for o in occultors]
    for i, occultor in enumerate([self.bodies[o] for o in occultors]): 
      r = occultor._r
      x = np.linspace(occultor.x[t] - r, occultor.x[t] + r, 1000)
      y = np.sqrt(r ** 2 - (x - occultor.x[t]) ** 2)
      pto[i] = ax.fill_between(x - x0, occultor.y[t] - y - y0, occultor.y[t] + y - y0, 
                               color = 'lightgray', zorder = 99 + i, lw = 1,
                               alpha = occultor_alpha)
      pto[i].set_edgecolor('k')
    
    return ax, pto
  
  def plot_lightcurve(self, wavelength = 15.):
    '''
    
    '''
    
    if not self.settings.quiet:
      print("Plotting the light curve...")
    
    # Plot
    fig, ax = pl.subplots(1, figsize = (12, 4))
    time = self.bodies[0].time
    assert (wavelength >= self.bodies[0].wavelength[0]) and (wavelength >= self.bodies[0].wavelength[-1]), "Wavelength value outside of computed grid."
    w = np.argmax(1e-6 * wavelength <= self.bodies[0].wavelength)
    flux = np.sum([b.flux[:,w] for b in self.bodies], axis = 0)
    flux /= np.nanmedian(flux)
    curve, = ax.plot(time, flux, 'k-', lw = 1, picker = 10)
    fig.canvas.mpl_connect('pick_event', self._onpick)
    
    # Appearance
    ax.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    ax.set_ylabel(r'Normalized Flux @ %.1f$\mathbf{\mu}$m' % wavelength, fontweight = 'bold', fontsize = 10)
    ax.get_yaxis().set_major_locator(MaxNLocator(4))
    ax.get_xaxis().set_major_locator(MaxNLocator(4))
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
      tick.set_fontsize(8)
    ax.ticklabel_format(useOffset = False)
    
    # Limits
    ymax = np.nanmax(flux)
    ymin = np.nanmin(flux)
    yrng = ymax - ymin
    ax.set_ylim(ymin - 0.2 * yrng, ymax + 0.2 * yrng)
    ax.margins(0, None)
    
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
        tmin = body.time[t][np.argmin(flux[t])]
        fmin = np.min(flux[t])
        for n, occultor in enumerate([self.bodies[o] for o in occultors]):
          ax.annotate("%s" % body.name, xy = (tmin, fmin), ha = 'center',
                      va = 'center', color = occultor.color, fontweight = 'bold',
                      fontsize = 10, xytext = (0, -15), textcoords = 'offset points')
    
    return fig, ax
  
  def _onpick(self, event):
    '''
    
    '''
    
    index = event.ind[len(event.ind) // 2]
    for body in self.bodies:
      for occultation in body._inds:
        if index in occultation:
          self.plot_occultation(body.name, body.time[index])
    pl.show()
  
  @property
  def flux(self):
    '''
    The total flux of the system computed on a grid of time and wavelength.
    
    '''
    
    return np.sum([b.flux for b in self.bodies], axis = 0)
    
  @property
  def time(self):
    '''
    Time array in days.
    
    '''
    
    return self.bodies[0].time
  
  @property
  def wavelength(self):
    '''
    Wavelength array in microns.
    
    '''
    
    return self.bodies[0].wavelength * 1.e6

class Animation(object):
  '''
  An animation class for occultation movies.
  
  '''
  
  def __init__(self, t, fig, axim, tracker, pto, ptb, body, bodies, occultors, 
               interval = 50, gifname = None, quiet = False):
    '''
    
    '''
    
    self.t = t
    self.fig = fig
    self.axim = axim
    self.tracker = tracker
    self.pto = pto
    self.ptb = ptb
    self.body = body
    self.bodies = bodies
    self.occultors = occultors
    self.pause = True
    self.animation = animation.FuncAnimation(self.fig, self.animate, frames = 100, 
                                             interval = interval, repeat = True)
    self.fig.canvas.mpl_connect('button_press_event', self.toggle)
    
    # Save?
    if gifname is not None:
      self.pause = False
      if not gifname.endswith('.gif'):
        gifname += '.gif'
      if not quiet:
        print("Saving %s..." % gifname)
      self.animation.save(gifname, writer = 'imagemagick', fps = 20, dpi = 150)
      self.pause = True
      
  def toggle(self, event):
    '''
    
    '''
    
    self.pause ^= True
    
  def animate(self, j):
    '''
    
    '''
    
    if not self.pause:
      
      # Normalize the time index
      j = int(j * len(self.t) / 100.)
      
      # Time tracker
      self.tracker.set_xdata(self.bodies[0].time[self.t[j]])
      
      # Occultor images
      x0 = self.body.x[self.t[j]]
      y0 = self.body.y[self.t[j]]
      for k, occultor in enumerate(self.occultors): 
        r = occultor._r
        x = np.linspace(occultor.x[self.t[j]] - r, occultor.x[self.t[j]] + r, 1000)
        y = np.sqrt(r ** 2 - (x - occultor.x[self.t[j]]) ** 2)
        try:
          self.pto[k].remove()
        except:
          pass
        self.pto[k] = self.axim.fill_between(x - x0, occultor.y[self.t[j]] - y - y0, occultor.y[self.t[j]] + y - y0, color = 'lightgray', zorder = 99 + k, lw = 1)
        self.pto[k].set_edgecolor('k')
      
      # Body orbits
      for k, b in enumerate(self.bodies):
        self.ptb[k].set_xdata(b.x[self.t[j]])
        self.ptb[k].set_ydata(b.z[self.t[j]])