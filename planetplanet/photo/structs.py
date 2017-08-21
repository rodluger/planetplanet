#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
structs.py |github|
-------------------

  .. role:: raw-html(raw)
     :format: html

  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/photo/structs.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from ..constants import *
import ctypes
import numpy as np

__all__ = ['BODY', 'SETTINGS', 'Star', 'Planet', 'Moon']

class BODY(ctypes.Structure):
  '''
  The class containing all the input planet/star parameters. This is a :py:mod:`ctypes` interface
  to the C :py:obj:`BODY` struct. Users should instantiate these via the :py:func:`Star`,
  :py:func:`Planet`, and :py:func:`Moon` functions.

  '''

  #: All the fields
  _fields_ = [("_m", ctypes.c_double),
              ("_per", ctypes.c_double),
              ("_inc", ctypes.c_double),
              ("_ecc", ctypes.c_double),
              ("_w", ctypes.c_double),
              ("_Omega", ctypes.c_double),
              ("a", ctypes.c_double),
              ("tperi0", ctypes.c_double),
              ("_r", ctypes.c_double),
              ("albedo", ctypes.c_double),
              ("_teff", ctypes.c_double),
              ("tnight", ctypes.c_double),
              ("_phasecurve", ctypes.c_int),
              ("_blackbody", ctypes.c_int),
              ("_Lambda", ctypes.c_double),
              ("_Phi", ctypes.c_double),
              ("_host", ctypes.c_int),
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
              ("_flux", ctypes.POINTER(ctypes.c_double)),
              ("_total_flux", ctypes.POINTER(ctypes.c_double)),
              ("_custommap", ctypes.c_int),
              ("_radiancemap", ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)),
              ]

  def __init__(self, name, body_type, **kwargs):
    '''

    :param str name: The name of the body.
    :param str body_type: One of :py:obj:`planet`, :py:obj:`star`, or :py:obj:`moon`.
    :param kwargs: Any body-specific :py:obj:`kwargs`. See :py:func:`Star`, \
    :py:func:`Planet`, and :py:func:`Moon`.

    '''
    
    # Check
    self.name = name
    self.body_type = body_type
    assert body_type in ['planet', 'star', 'moon'], "Argument `body_type` must be `planet`, `moon`, or `star`."

    # User
    self.m = kwargs.pop('m', 1.)
    self.r = kwargs.pop('r', 1.)
    self.t0 = kwargs.pop('t0', 0.)
    self.ecc = kwargs.pop('ecc', 0.)
    self.w = kwargs.pop('w', 0.)
    self.Omega = kwargs.pop('Omega', 0.)
    self.inc = kwargs.pop('inc', 90.)

    # These defaults are different depending on body type
    if self.body_type in ['planet', 'moon']:
      self.airless = kwargs.pop('airless', True)
      self.nz = kwargs.pop('nz', 11)
      self.per = kwargs.pop('per', 3.)
      self.albedo = kwargs.pop('albedo', 0.3)
      self.teff = 0
      if self.airless:
        self.tnight = kwargs.pop('tnight', 40.)
        self.limbdark = []
        self.Lambda = kwargs.pop('Lambda', 0)
        self.Phi = kwargs.pop('Phi', 0)
      else:
        self.tnight = 0
        self.limbdark = kwargs.pop('limbdark', [])
        self.Lambda = 0
        self.Phi = 0
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
      self.Lambda = 0
      self.Phi = 0
      self.color = kwargs.pop('color', 'k')
      
    # User-defined radiance map
    self.radiancemap = kwargs.get('radiancemap', None)
    
    # Settings for moons
    if self.body_type == 'moon':
      self.host = kwargs.pop('host', None)
      if not type(self.host) is str:
        raise ValueError("Please specify the `host` kwarg for all moons.")
    else:
      self.host = 0

    # C stuff, computed in `System` class
    self.nt = 0
    self.nw = 0
    self.nu = 0
    self.tperi0 = 0.
    self.a = 0.

    # Python stuff
    self._inds = []
    self._computed = False

  @property
  def m(self):
    if self.body_type in ['planet', 'moon']:
      return self._m
    elif self.body_type == 'star':
      return self._m / MSUNMEARTH

  @m.setter
  def m(self, val):
    if self.body_type in ['planet', 'moon']:
      self._m = val
    elif self.body_type == 'star':
      self._m = val * MSUNMEARTH

  @property
  def r(self):
    if self.body_type in ['planet', 'moon']:
      return self._r
    elif self.body_type == 'star':
      return self._r / RSUNREARTH

  @r.setter
  def r(self, val):
    if self.body_type in ['planet', 'moon']:
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
    # Force an update to the time of pericenter passage
    self.ecc = self._ecc

  @property
  def Omega(self):
    return self._Omega * 180 / np.pi

  @Omega.setter
  def Omega(self, val):
    self._Omega = val * np.pi / 180.

  @property
  def Lambda(self):
    return self._Lambda * 180 / np.pi

  @Lambda.setter
  def Lambda(self, val):
    self._Lambda = val * np.pi / 180.

  @property
  def Phi(self):
    return self._Phi * 180 / np.pi

  @Phi.setter
  def Phi(self, val):
    self._Phi = val * np.pi / 180.

  @property
  def teff(self):
    return self._teff

  @teff.setter
  def teff(self, val):
    self._teff = val
    if self._teff == 0:
      self.u = []

  @property
  def t0(self):
    return self._t0

  @t0.setter
  def t0(self, val):
    self._t0 = val
    # Force an update to the time of pericenter passage
    self.ecc = self._ecc

  @property
  def per(self):
    return self._per

  @per.setter
  def per(self, val):
    self._per = val
    # Force an update to the time of pericenter passage
    self.ecc = self._ecc

  @property
  def ecc(self):
    return self._ecc

  @ecc.setter
  def ecc(self, val):
    '''
    We need to update the time of pericenter passage whenever the eccentricty,
    longitude of pericenter, period, or time of transit changes. See the appendix in
    Shields et al. (2016).

    '''

    self._ecc = val
    fi = (3 * np.pi / 2.) - self._w
    self.tperi0 = self.t0 + (self.per * np.sqrt(1. - self.ecc * self.ecc) / (2. * np.pi) * (self.ecc * np.sin(fi) /
             (1. + self.ecc * np.cos(fi)) - 2. / np.sqrt(1. - self.ecc * self.ecc) *
             np.arctan2(np.sqrt(1. - self.ecc * self.ecc) * np.tan(fi/2.), 1. + self.ecc)))

  def M(self, t):
    '''
    The mean anomaly at a time `t`.

    '''

    return 2. * np.pi / self.per * ((t - self.tperi0) % self.per)


class SETTINGS(ctypes.Structure):
  '''
  The class that contains the model settings. This class is used internally; settings
  should be specified as :py:obj:`kwargs` to :py:class:`System` or by assignment once
  a :py:class:`System` object has been instantiated.

  :param bool nbody: Uses the :py:obj:`REBOUND` N-body code to compute orbits. Default :py:obj:`False`
  :param float keptol: Kepler solver tolerance. Default `1.e-15`
  :param int maxkepiter: Maximum number of Kepler solver iterations. Default `100`
  :param str kepsolver: Kepler solver (`newton` | `mdfast`). Default `newton`
  :param float timestep: Timestep in days for the N-body solver. Default `0.01`
  :param bool adaptive: Adaptive grid for limb-darkened bodies? Default :py:obj:`True`
  :param bool quiet: Suppress output? Default :py:obj:`False`
  :param float mintheta: Absolute value of the minimum phase angle in degrees. Below this \
         angle, elliptical boundaries of constant surface brightness on the planet surface are \
         treated as vertical lines. Default `1.`
  :param int maxvertices: Maximum number of vertices allowed in the area computation. Default `999`
  :param int maxfunctions: Maximum number of functions allowed in the area computation. Default `999`
  :param int oversample: Oversampling factor for each exposure. Default `1`
  :param float distance: Distance to the system in parsecs. Default `10.`
  :param bool circleopt: Solve the simpler quadratic problem for circle-ellipse intersections when \
         the axes of the ellipse are equal to within :math:`10^{-10}`? Default :py:obj:`True`
  :param bool batmanopt: Use the :py:mod:`batman` algorithm to compute light curves of radially \
         symmetric bodies? This can significantly speed up the code. Default :py:obj:`True`
  :param str integrator: The N-body integrator (:py:obj:`whfast` | :py:obj:`ias15`) to use. Default :py:obj:`whfast`

  '''

  _fields_ = [("_nbody", ctypes.c_int),
              ("_integrator", ctypes.c_int),
              ("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("_kepsolver", ctypes.c_int),
              ("timestep", ctypes.c_double),
              ("_adaptive", ctypes.c_int),
              ("_circleopt", ctypes.c_int),
              ("_batmanopt", ctypes.c_int),
              ("_quarticsolver", ctypes.c_int),
              ("_quiet", ctypes.c_int),
              ("_mintheta", ctypes.c_double),
              ("maxvertices", ctypes.c_int),
              ("maxfunctions", ctypes.c_int),
              ("distance", ctypes.c_double)]

  def __init__(self, **kwargs):
    self.nbody = kwargs.pop('nbody', False)
    self.integrator = kwargs.pop('integrator', 'whfast')
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = kwargs.pop('kepsolver', 'newton')
    self.timestep = kwargs.pop('timestep', 0.01)
    self.adaptive = kwargs.pop('adaptive', True)
    self.circleopt = kwargs.pop('circleopt', True)
    self.batmanopt = kwargs.pop('batmanopt', True)
    self.quarticsolver = kwargs.pop('quarticsolver', 'gsl')
    self.quiet = kwargs.pop('quiet', False)
    self.mintheta = kwargs.pop('mintheta', 1.)
    self.maxvertices = kwargs.pop('maxvertices', 999)
    self.maxfunctions = kwargs.pop('maxfunctions', 999)
    self.oversample = max(1, kwargs.pop('oversample', 1))
    self.distance = kwargs.pop('distance', 10.)

  @property
  def params(self):
    return ['nbody', 'integrator', 'keptol', 'maxkepiter', 'kepsolver',
            'timestep', 'adaptive', 'circleopt', 'batmanopt', 'quarticsolver', 'quiet', 
            'mintheta', 'maxvertices', 'maxfunctions', 'oversample', 'distance']

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
  def circleopt(self):
    return bool(self._circleopt)

  @circleopt.setter
  def circleopt(self, val):
    self._circleopt = int(val)

  @property
  def batmanopt(self):
    return bool(self._batmanopt)

  @batmanopt.setter
  def batmanopt(self, val):
    self._batmanopt = int(val)

  @property
  def integrator(self):
    if self._integrator == REB_INTEGRATOR_WHFAST:
      return 'whfast'
    elif self._integrator == REB_INTEGRATOR_IAS15:
      return 'ias15'

  @integrator.setter
  def integrator(self, val):
    if val.lower() == 'whfast':
      self._integrator = REB_INTEGRATOR_WHFAST
    elif val.lower() == 'ias15':
      self._integrator = REB_INTEGRATOR_IAS15
    else:
      raise ValueError("Unsupported integrator.")
  
  @property
  def quarticsolver(self):
    if self._quarticsolver == QGSL:
      return 'gsl'
      
  @quarticsolver.setter
  def quarticsolver(self, val):
    if val.lower() == 'gsl':
      self._quarticsolver = QGSL
    else:		
      raise ValueError("Unsupported solver.")

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

  Returns a :py:class:`BODY` instance of type :py:obj:`star`.

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

  return BODY(name, 'star', **kwargs)

def Planet(name, **kwargs):
  '''

  Returns a :py:class:`BODY` instance of type :py:obj:`planet`.

  :param str name: A unique identifier for this planet
  :param float m: Mass in Earth masses. Default `1.`
  :param float r: Radius in Earth radii. Default `1.`
  :param float per: Orbital period in days. Default `3.`
  :param float inc: Orbital inclination in degrees. Default `90.`
  :param float ecc: Orbital eccentricity. Default `0.`
  :param float w: Longitude of pericenter in degrees. `0.`
  :param float Omega: Longitude of ascending node in degrees. `0.`
  :param float t0: Time of transit in days. Default `0.`
  :param bool phasecurve: Compute the phasecurve for this planet? Default :py:obj:`False`
  :param bool airless: Treat this as an airless planet? If :py:obj:`True`, computes light curves \
         in the instant re-radiation limit, where the surface brightness is proportional \
         to the cosine of the zenith angle (the angle between the line connecting \
         the centers of the planet and the star and the line connecting the center of the \
         planet and a given point on its surface. A fixed nightside temperature may be specified \
         via the `tnight` kwarg. If :py:obj:`False`, treats the planet as a limb-darkened blackbody. \
         Default :py:obj:`True`
  :param float albedo: Planetary albedo (airless limit). Default `0.3`
  :param float tnight: Nightside temperature in Kelvin (airless limit). Default `40`
  :param array_like limbdark: The limb darkening coefficients (thick atmosphere limit). These are the coefficients \
         in the Taylor expansion of `(1 - mu)`, starting with the first order (linear) \
         coefficient, where `mu = cos(theta)` is the radial coordinate on the surface of \
         the star. Each coefficient may either be a scalar, in which case limb darkening is \
         assumed to be grey (the same at all wavelengths), or a callable whose single argument \
         is the wavelength in microns. Default is `[1.0]`, a grey linear limb darkening law.
  :param float Lambda: Latitudinal hotspot offset in degrees, with positive values corresponding to a northward \
         shift. Airless bodies only. Default `0.`
  :param float Phi: Longitudinal hotspot offset in degrees, with positive values corresponding to an eastward \
         shift. Airless bodies only. Default `0.`
  :param int nz: Number of zenith angle slices. Default `11`
  :param str color: Object color (for plotting). Default `r`

  '''

  return BODY(name, 'planet', **kwargs)

def Moon(name, host, **kwargs):
  '''

  Returns a :py:class:`BODY` instance of type :py:obj:`moon`.

  :param str name: A unique identifier for this moon
  :param str host: The name of the moon's host planet
  :param float m: Mass in Earth masses. Default `1.`
  :param float r: Radius in Earth radii. Default `1.`
  :param float per: Orbital period in days. Default `3.`
  :param float inc: Orbital inclination in degrees. Default `90.`
  :param float ecc: Orbital eccentricity. Default `0.`
  :param float w: Longitude of pericenter in degrees. `0.`
  :param float Omega: Longitude of ascending node in degrees. `0.`
  :param float t0: Time of transit in days. Default `0.`
  :param bool phasecurve: Compute the phasecurve for this planet? Default :py:obj:`False`
  :param bool airless: Treat this as an airless planet? If :py:obj:`True`, computes light curves \
         in the instant re-radiation limit, where the surface brightness is proportional \
         to the cosine of the zenith angle (the angle between the line connecting \
         the centers of the planet and the star and the line connecting the center of the \
         planet and a given point on its surface. A fixed nightside temperature may be specified \
         via the `tnight` kwarg. If :py:obj:`False`, treats the planet as a limb-darkened blackbody. \
         Default :py:obj:`True`
  :param float albedo: Planetary albedo (airless limit). Default `0.3`
  :param float tnight: Nightside temperature in Kelvin (airless limit). Default `40`
  :param array_like limbdark: The limb darkening coefficients (thick atmosphere limit). These are the coefficients \
         in the Taylor expansion of `(1 - mu)`, starting with the first order (linear) \
         coefficient, where `mu = cos(theta)` is the radial coordinate on the surface of \
         the star. Each coefficient may either be a scalar, in which case limb darkening is \
         assumed to be grey (the same at all wavelengths), or a callable whose single argument \
         is the wavelength in microns. Default is `[1.0]`, a grey linear limb darkening law.
  :param float Lambda: Latitudinal hotspot offset in degrees, with positive values corresponding to a northward \
         shift. Airless bodies only. Default `0.`
  :param float Phi: Longitudinal hotspot offset in degrees, with positive values corresponding to an eastward \
         shift. Airless bodies only. Default `0.`
  :param int nz: Number of zenith angle slices. Default `11`
  :param str color: Object color (for plotting). Default `r`

  '''

  return BODY(name, 'moon', host = host, **kwargs)
