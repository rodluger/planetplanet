#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
example.py
----------


'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
import re
import numpy as np
np.seterr(invalid = 'ignore')
cmap = pl.get_cmap('RdBu_r')
TOL = 1e-10

def Roots(occultor, ellipse):
  '''

  '''
  
  # Get the params
  a = ellipse.a
  b = ellipse.b
  xE = ellipse.x0
  yE = ellipse.y0
  r = occultor.r
  
  # Get the coefficients
  A = a ** 2 / b ** 2 - 1
  B = -2 * a ** 2 * xE / b ** 2
  C = r ** 2 - yE ** 2 - a ** 2 + a ** 2 / b ** 2 * xE ** 2
  D = 4 * yE ** 2 * a ** 2 / b ** 2
  c4 = A ** 2
  c3 = 2 * A * B
  c2 = 2 * A * C + B ** 2 + D
  c1 = 2 * B * C - 2 * D * xE
  c0 = C ** 2 - (b ** 2 - xE ** 2) * D
  
  # Solve numerically
  roots = np.roots([c4, c3, c2, c1, c0])
  
  # Filter out imaginary roots
  good_roots = []
  for x in roots:
    if (x.imag == 0):
      good_roots.append(x.real)
      
  return good_roots

class Circle(object):
  '''
  
  '''
  
  def __init__(self, x0, y0, r):
    '''
    
    '''
    
    self.x0 = x0
    self.y0 = y0
    self.r = r
  
  @property
  def xmin(self):
    '''
    
    '''
    
    return self.x0 - self.r
    
  @property
  def xmax(self):
    '''
    
    '''
    
    return self.x0 + self.r
  
  @property
  def vertices(self):
    '''
    
    '''
    
    v1 = (self.x0 - self.r, self.y0)
    v2 = (self.x0 + self.r, self.y0)
    return [v1, v2]
    
  def val_upper(self, x):
    '''
    
    '''
    
    A = self.r ** 2 - (x - self.x0) ** 2
    if hasattr(x, '__len__'):
      A[np.abs(A) < TOL] = 0
      A[np.where((x > self.xmax) | (x < self.xmin))] = np.nan
    else:
      if np.abs(A) < TOL:
        A = 0
      if (x > self.xmax) or (x < self.xmin):
        return np.nan
    return self.y0 + np.sqrt(A)

  def val_lower(self, x):
    '''
    
    '''
    
    A = self.r ** 2 - (x - self.x0) ** 2
    if hasattr(x, '__len__'):
      A[np.abs(A) < TOL] = 0
      A[np.where((x > self.xmax) | (x < self.xmin))] = np.nan
    else:
      if np.abs(A) < TOL:
        A = 0
      if (x > self.xmax) or (x < self.xmin):
        return np.nan
    return self.y0 - np.sqrt(A)
  
  def int_upper(self, x0, x1):
    '''
    
    '''
    
    # Check if out of bounds
    if (x0 > self.xmax) or (x0 < self.xmin) or (x1 > self.xmax) or (x1 < self.xmin):
      return np.nan
    
    F = [0, 0]
    for i, x in enumerate((x0, x1)):
      y = x - self.x0
      z = np.sqrt((self.r - y) * (self.r + y))
      F[i] = (1 / 2.) * (y * z + self.r ** 2 * np.arctan(y / z)) + self.y0 * x
    return F[1] - F[0]

  def int_lower(self, x0, x1):
    '''
    
    '''
    
    # Check if out of bounds
    if (x0 > self.xmax) or (x0 < self.xmin) or (x1 > self.xmax) or (x1 < self.xmin):
      return np.nan
    
    F = [0, 0]
    for i, x in enumerate((x0, x1)):
      y = x - self.x0
      z = np.sqrt((self.r - y) * (self.r + y))
      F[i] = -(1 / 2.) * (y * z + self.r ** 2 * np.arctan(y / z)) + self.y0 * x
    return F[1] - F[0]
  
  @property
  def curves(self):
    '''
    
    '''
    
    return [self.val_lower, self.val_upper]
  
  @property
  def integrals(self):
    '''
    
    '''
    
    return [self.int_lower, self.int_upper]
      
class Ellipse(object):
  '''
  
  '''
  
  def __init__(self, planet, occultor, latitude, theta):
    '''
    
    '''

    # Generate the ellipse
    self.a = planet.r * np.abs(np.sin(latitude))
    self.b = self.a * np.abs(np.sin(theta))
    self.x0 = planet.x0 - planet.r * np.cos(latitude) * np.cos(theta)
    self.y0 = planet.y0
    self.latitude = latitude
    
    # The x position of the limb
    self.xlimb = planet.r * np.cos(latitude) * np.sin(theta) * np.tan(theta)
    
    # Compute the vertices / intersection points
    self.vertices = []
    self.intersections = []

    # Ellipse x minimum
    if ((theta > 0) and (self.b < self.xlimb)) or ((theta <= 0) and (self.b > self.xlimb)):
      self.vertices.append((self.x0 - self.b, self.y0))
      self.xmin = self.x0 - self.b
    else:
      self.xmin = self.x0 - self.xlimb
      
    # Ellipse x maximum
    if ((theta > 0) and (self.b > -self.xlimb)) or ((theta <= 0) and (self.b < -self.xlimb)):
      self.vertices.append((self.x0 + self.b, self.y0))
      self.xmax = self.x0 + self.b
    else:
      self.xmax = self.x0 - self.xlimb
    
    # Ellipse-planet vertices
    if self.b ** 2 >= self.xlimb ** 2:
      x = self.x0 - self.xlimb
      y = self.a / self.b * np.sqrt(self.b ** 2 - self.xlimb ** 2)
      self.vertices.append((x, self.y0 - y))
      self.vertices.append((x, self.y0 + y))
      
    # Ellipse-occultor intersections
    for root in Roots(occultor, self):
      if ((theta > 0) and (root >= self.x0 - self.xlimb - TOL)) or ((theta <= 0) and (root <= self.x0 - self.xlimb + TOL)):
        eup = self.val_upper(root)
        elo = self.val_lower(root)
        oup = occultor.val_upper(root)
        olo = occultor.val_lower(root)
        if (np.abs(eup - oup) < TOL):
          self.intersections.append((root, eup))
        elif (np.abs(eup - olo) < TOL):
          self.intersections.append((root, eup))
        elif (np.abs(elo - oup) < TOL):
          self.intersections.append((root, elo))
        elif (np.abs(elo - olo) < TOL):
          self.intersections.append((root, elo))
    
  def val_upper(self, x):
    '''
    
    '''
    
    A = self.b ** 2 - (x - self.x0) ** 2
    if hasattr(x, '__len__'):
      A[np.abs(A) < TOL] = 0
      A[np.where((x > self.xmax) | (x < self.xmin))] = np.nan
    else:
      if np.abs(A) < TOL:
        A = 0
      if (x > self.xmax) or (x < self.xmin):
        return np.nan
    return self.y0 + (self.a / self.b) * np.sqrt(A)

  def val_lower(self, x):
    '''
    
    '''
    
    A = self.b ** 2 - (x - self.x0) ** 2
    if hasattr(x, '__len__'):
      A[np.abs(A) < TOL] = 0
      A[np.where((x > self.xmax) | (x < self.xmin))] = np.nan
    else:
      if np.abs(A) < TOL:
        A = 0
      if (x > self.xmax) or (x < self.xmin):
        return np.nan
    return self.y0 - (self.a / self.b) * np.sqrt(A)
  
  def int_upper(self, x0, x1):
    '''
    
    '''
    
    # Check if out of bounds
    if (x0 > self.xmax) or (x0 < self.xmin) or (x1 > self.xmax) or (x1 < self.xmin):
      return np.nan
    
    F = [0, 0]
    for i, x in enumerate((x0, x1)):
      y = x - self.x0
      z = np.sqrt((self.b - y) * (self.b + y))
      F[i] = (self.a / (2 * self.b)) * (y * z + self.b ** 2 * np.arctan(y / z)) + self.y0 * x
    return F[1] - F[0]

  def int_lower(self, x0, x1):
    '''
    
    '''
    
    # Check if out of bounds
    if (x0 > self.xmax) or (x0 < self.xmin) or (x1 > self.xmax) or (x1 < self.xmin):
      return np.nan
    
    F = [0, 0]
    for i, x in enumerate((x0, x1)):
      y = x - self.x0
      z = np.sqrt((self.b - y) * (self.b + y))
      F[i] = -(self.a / (2 * self.b)) * (y * z + self.b ** 2 * np.arctan(y / z)) + self.y0 * x
    return F[1] - F[0]

  @property
  def curves(self):
    '''
    
    '''
    
    return [self.val_lower, self.val_upper]
  
  @property
  def integrals(self):
    '''
    
    '''
    
    return [self.int_lower, self.int_upper]

class Occultor(Circle):
  '''
  The occultor is simply a Circle centered at the origin.
  
  '''
  
  def __init__(self, r):
    '''
    
    '''
    
    super(Occultor, self).__init__(0, 0, r)

class Planet(Circle):
  '''
  
  '''
  
  def __init__(self, x0, y0, r, theta, occultor, n = 11, noon = 0.3, midnight = 0.01):
    '''
    
    '''
    
    self._r = r
    self._theta = theta
    self._n = n
    self._latitudes = np.linspace(0, np.pi, self._n + 2)[1:-1]
    self._noon = noon
    self._midnight = midnight    
    self.x0 = x0
    self.y0 = y0
    self.occultor = occultor
    self.compute_total()
    self.compute_delta()
  
  @property
  def r(self):
    return self._r
  
  @r.setter
  def r(self, val):
    self._r = val
    self.compute_total()

  @property
  def theta(self):
    return self._theta
  
  @theta.setter
  def theta(self, val):
    self._theta = val
    self.compute_total()

  @property
  def n(self):
    return self._n
  
  @n.setter
  def n(self, val):
    self._n = val
    self._latitudes = np.linspace(0, np.pi, self._n + 2)[1:-1]
    self.compute_total()

  @property
  def latitudes(self):
    return self._latitudes

  @property
  def noon(self):
    return self._noon
  
  @noon.setter
  def noon(self, val):
    self._noon = val
    self.compute_total()

  @property
  def midnight(self):
    return self._midnight
  
  @midnight.setter
  def midnight(self, val):
    self._midnight = val
    self.compute_total()
    
  @property
  def vertices(self):
    '''
    
    '''
    
    # The planet's own vertices
    vertices = [(self.x0 - self.r, self.y0), (self.x0 + self.r, self.y0)]
      
    # Add the occultor's vertices
    vertices.extend(self.occultor.vertices)
    
    # Add the ellipse vertices
    for ellipse in self.ellipses:
      vertices.extend(ellipse.vertices)
    
    # Now discard the vertices that are outside either the occultor or the planet
    vin = []
    for v in vertices:
      x, y = v
      if (x ** 2 + y ** 2 <= self.occultor.r ** 2 + TOL) and ((x - self.x0) ** 2 + (y - self.y0) ** 2 <= self.r ** 2 + TOL):
        vin.append((x, y))
    vertices = list(vin)
    
    # Compute vertices of intersection with the occultor
    # Adapted from http://mathworld.wolfram.com/Circle-CircleIntersection.html
    d = np.sqrt(self.x0 ** 2 + self.y0 ** 2)
    if d <= (self.occultor.r + self.r):
      A = (-d + self.r - self.occultor.r) * (-d - self.r + self.occultor.r) * (-d + self.r + self.occultor.r) * (d + self.r + self.occultor.r)
      if A >= 0:
        y = np.sqrt(A) / (2 * d)
        x = (d ** 2 - self.r ** 2 + self.occultor.r ** 2) / (2 * d)
        cost = -self.x0 / d
        sint = -self.y0 / d
        x1 = -x * cost + y * sint
        y1 = -x * sint - y * cost
        x2 = -x * cost - y * sint
        y2 = -x * sint + y * cost
        vertices.extend([(x1, y1), (x2, y2)])
    
    # Get the points of intersection between the ellipses and the planet
    for ellipse in self.ellipses:
      vertices.extend(ellipse.intersections)
    
    # Sort them
    vertices = sorted(list(set(vertices)))

    return vertices
  
  @property
  def curves(self):
    '''
    
    '''
    
    curves = [self.val_lower, self.val_upper]
    curves.extend(self.occultor.curves)
    for ellipse in self.ellipses:
      curves.extend(ellipse.curves)
    return curves
  
  @property
  def integrals(self):
    '''
    
    '''
    
    integrals = [self.int_lower, self.int_upper]
    integrals.extend(self.occultor.integrals)
    for ellipse in self.ellipses:
      integrals.extend(ellipse.integrals)
    return integrals
  
  @property
  def flux(self):
    '''
    
    '''
    
    return self.total_flux + self.delta_flux
  
  @property
  def norm_flux(self):
    '''
    
    '''
    
    return self.flux / self.total_flux
  
  def surf_brightness(self, lat):
    '''
    
    '''

    grid = np.concatenate(([0], self.latitudes, [np.pi + TOL]))
    if lat < np.pi / 2:
      i = np.argmax(lat < grid) - 1
    else:
      i = np.argmax(lat < grid)
    cos = np.cos(grid[i])
    return (0.5 * (self.noon - self.midnight) * (cos + 1) + self.midnight)

  def get_latitude(self, x, y):
    '''
  
    '''
    
    # Compute the latitude. We will solve
    # a quadratic equation for z = sin(lat) **  2  
    alpha = (x - self.x0) / (self.r * np.sin(self.theta))
    beta = 1 / np.tan(self.theta)
    gamma = (y - self.y0) / self.r
    c1 = 4 * alpha ** 2 * beta ** 2
    c2 = 1 + beta ** 2
    c3 = alpha ** 2 + gamma ** 2 + beta ** 2
    b = (c1 - 2 * c2 * c3) / c2 ** 2
    c = (c3 ** 2 - c1) / c2 ** 2
    z0 = (-b + np.sqrt(b ** 2 - 4 * c)) / 2
    z1 = (-b - np.sqrt(b ** 2 - 4 * c)) / 2
    
    # Where is the terminator for this value of `y`?
    if (self.theta <= 0):
      xterm = self.x0 - (np.abs(np.sin(self.theta))) * np.sqrt(self.r ** 2 - (y - self.y0) ** 2)
    else:
      xterm = self.x0 + (np.abs(np.sin(self.theta))) * np.sqrt(self.r ** 2 - (y - self.y0) ** 2)

    # Are we on the dayside?
    if x <= xterm:
    
      # We have two possible solutions. But only one is on the
      # observer's side of the planet. TODO: Speed this up?
      for z in (z0, z1):
        if (z >= 0) and (z <= 1):
          a = self.r * np.abs(np.sqrt(z))
          b = a * np.abs(np.sin(self.theta))
          x0 = self.x0 - self.r * np.sqrt(1 - z) * np.cos(self.theta)
          y0 = self.y0
          dx = (b / a) * np.sqrt(a ** 2 - (y - y0) ** 2)
          xlimb = self.r * np.sqrt(1 - z) * np.sin(self.theta) * np.tan(self.theta)
          if ((self.theta > 0) and (b < xlimb)) or ((self.theta <= 0) and (b > xlimb)):
            xmin = x0 - b
          else:
            xmin = x0 - xlimb
          if ((self.theta > 0) and (b > -xlimb)) or ((self.theta <= 0) and (b < -xlimb)):
            xmax = x0 + b
          else:
            xmax = x0 - xlimb
          if (x >= xmin - TOL) and (x <= xmax + TOL): 
            if ((np.abs(x - (x0 + dx)) < TOL) or (np.abs(x - (x0 - dx)) < TOL)):
              return np.arcsin(np.sqrt(z))
    
    # Or the nightside?
    else:
      
      # We have two possible solutions. But only one is on the
      # observer's side of the planet. TODO: Speed this up?
      for z in (z0, z1):
        if (z >= 0) and (z <= 1):
          a = self.r * np.abs(np.sqrt(z))
          b = a * np.abs(np.sin(self.theta))
          x0 = self.x0 + self.r * np.sqrt(1 - z) * np.cos(self.theta)
          y0 = self.y0
          dx = (b / a) * np.sqrt(a ** 2 - (y - y0) ** 2)
          xlimb = -self.r * np.sqrt(1 - z) * np.sin(self.theta) * np.tan(self.theta)
          if ((self.theta > 0) and (b < xlimb)) or ((self.theta <= 0) and (b > xlimb)):
            xmin = x0 - b
          else:
            xmin = x0 - xlimb
          if ((self.theta > 0) and (b > -xlimb)) or ((self.theta <= 0) and (b < -xlimb)):
            xmax = x0 + b
          else:
            xmax = x0 - xlimb
          if (x >= xmin - TOL) and (x <= xmax + TOL): 
            if ((np.abs(x - (x0 + dx)) < TOL) or (np.abs(x - (x0 - dx)) < TOL)):
              return np.pi - np.arcsin(np.sqrt(z))
    
    return np.nan

  def compute_total(self):
    '''
    
    '''
    
    # Save current values
    x0 = self.x0
    y0 = self.y0
    r = self.occultor.r
    
    # Compute missing flux when occultor completely
    # blocks the planet: this is the total flux.
    self.x0 = 0
    self.y0 = 0
    self.occultor.r = 1.5 * self.r
    self.compute_delta()
    self.total_flux = -self.delta_flux
    
    # Restore
    self.x0 = x0
    self.y0 = y0
    self.occultor.r = r
    
    # Finally, compute the new delta flux
    self.compute_delta()

  def compute_delta(self):
    '''
  
    '''
  
    # Avoid the singular point
    if np.abs(self.theta) < 1e-3:
      self.theta = 1e-3
    flux = []
    
    # Compute the ellipses
    self.ellipses = [Ellipse(self, self.occultor, lat, self.theta) for lat in self.latitudes]
    
    # Get all vertices, curves, and integrals
    vertices = self.vertices
    curves = self.curves
    integrals = self.integrals
    
    # Get tuples of integration limits
    limits = list(zip(vertices[:-1], vertices[1:]))

    # Loop over all the sub-regions and compute their fluxes
    for i, lim in enumerate(limits):
    
      # The integration limits
      (xleft, _), (xright, _) = lim
    
      # Check if they are identical
      if xright - TOL <= xleft + TOL:
        continue
      
      # Perturb them for numerical stability, as we
      # need to ensure the functions are finite-valued
      # at the limits.
      xleft += TOL
      xright -= TOL
    
      # Bisect the limits. Find the boundary functions that are
      # finite valued at this point and sort their integrals 
      # in order of increasing function value.
      x = (xleft + xright) / 2
      vals = []
      ints = []
      for curve, integral in zip(curves, integrals):
        y = curve(x)
        if np.isfinite(y) and (x ** 2 + y ** 2 <= self.occultor.r ** 2 + TOL) and ((x - self.x0) ** 2 + (y - self.y0) ** 2 <= self.r ** 2 + TOL):
          vals.append(y)
          ints.append(integral)
      order = np.argsort(np.array(vals))
      ints = [ints[j] for j in order]
      vals = [vals[j] for j in order]
    
      # The areas are just the difference of successive integrals
      for int1, int0, y1, y0 in zip(ints[1:], ints[:-1], vals[1:], vals[:-1]):
        area = int1(xleft, xright) - int0(xleft, xright)
      
        # Get the latitude of the midpoint and
        # the corresponding surface brightness
        y = (y1 + y0) / 2
        lat = self.get_latitude(x, y)
        f = self.surf_brightness(lat)
        flux.append(f * area)
      
        # DEBUG: Let's catch these NaNs
        if np.isnan(area):
          print("The integration returned NaN!")
          pl.close()
          fig = pl.figure(figsize = (6,6))
          x = np.linspace(self.x0 - self.r, self.x0 + self.r, 1000)
          pl.plot(x, self.y0 + np.sqrt(self.r ** 2 - (x - self.x0) ** 2), 'k-')
          pl.plot(x, self.y0 - np.sqrt(self.r ** 2 - (x - self.x0) ** 2), 'k-')
          x = np.linspace(-self.occultor.r, self.occultor.r, 1000)
          pl.plot(x, np.sqrt(self.occultor.r ** 2 - x ** 2), 'r-')
          pl.plot(x, -np.sqrt(self.occultor.r ** 2 - x ** 2), 'r-')
          for ellipse in self.ellipses:
            x = np.linspace(ellipse.x0 - ellipse.b, ellipse.x0 + ellipse.b, 1000)
            if self.theta > 0:
              x[x < ellipse.x0 - ellipse.xlimb] = np.nan
            else:
              x[x > ellipse.x0 - ellipse.xlimb] = np.nan
            pl.plot(x, ellipse.val_lower(x), 'k-')
            pl.plot(x, ellipse.val_upper(x), 'k-')          
          for v in self.vertices:
            pl.plot(v[0], v[1], 'o', color = 'k', ms = 3)
          pl.xlim(self.x0 - 1.5 * self.r, self.x0 + 1.5 * self.r)
          pl.ylim(self.y0 - 1.5 * self.r, self.y0 + 1.5 * self.r)
          pl.axvline(xleft, color = 'r', alpha = 0.3)
          pl.axvline(xright, color = 'r', alpha = 0.3)
          pl.show()
          import pdb; pdb.set_trace()
      
    # Total flux
    self.delta_flux = -np.sum(flux)

class Interact(object):

  def __init__(self, **kwargs):
    '''
  
    '''
    
    # Occultor (always at the origin)
    self.occultor = Occultor(0.5)
    
    # Planet (occulted)
    self.planet = Planet(1.35, -0.75, 1., np.pi / 8, self.occultor, **kwargs)
    
    # Set up the figure
    self.fig = pl.figure(figsize = (6,7))
    self.ax = pl.axes([0.125, 0.2, 0.8, 0.625])
    self.ax.set_xlim(-3,3)
    self.ax.set_ylim(-3,3)
    
    # Plotting arrays
    self.x = np.linspace(-1, 1, 1000)
    self.dummy = Circle(0, 0, 1)
    self.ylower = self.dummy.val_lower(self.x)
    self.yupper = self.dummy.val_upper(self.x)
    
    # The light curve
    self.axlc = pl.axes([0.125, 0.85, 0.8, 0.1])
    self.time = np.array(np.arange(100), dtype = float)
    self.flux = np.ones_like(self.time)
    self.lc, = self.axlc.plot(self.time, self.flux, 'k-', ms = 2, lw = 1, alpha = 0.75)
    self.lcr, = self.axlc.plot(self.time[-1], self.flux[-1], 'ro', ms = 3)
    self.lcl = self.axlc.annotate('%.3f' % self.flux[-1], xy = (101, self.flux[-1]), 
                                   fontsize = 6, color = 'r', ha = 'left', va = 'center')
    self.axlc.set_xlim(0, 110)
    self.axlc.set_xticks([])  
    
    # The theta slider
    self.axtheta = pl.axes([0.125, 0.035, 0.725, 0.02])
    self.sltheta = Slider(self.axtheta, r'$\theta$', -90, 90, valinit = self.planet.theta * 180 / np.pi)
    self.sltheta.on_changed(self.theta_changed)
    
    # The occultor radius slider
    self.axradius = pl.axes([0.125, 0.065, 0.725, 0.02])
    self.slradius = Slider(self.axradius, r'$R_o$', 0.01, 2, valinit = self.occultor.r)
    self.slradius.on_changed(self.radius_changed)
    
    # The occultor y0 slider
    self.axy0 = pl.axes([0.125, 0.095, 0.725, 0.02])
    self.sly0 = Slider(self.axy0, r'$y_o$', -3, 3, valinit = -self.planet.y0)
    self.sly0.on_changed(self.y0_changed)
    
    # The occultor x0 slider
    self.axx0 = pl.axes([0.125, 0.125, 0.725, 0.02])
    self.slx0 = Slider(self.axx0, r'$x_o$', -3, 3, valinit = -self.planet.x0)
    self.slx0.on_changed(self.x0_changed)
    
    # Show!
    self.update()
    pl.show()

  def style(self, lat):
    '''
  
    '''
  
    coslat = np.cos(lat)
    if np.abs(coslat) < TOL:
      return dict(color = 'k', ls = '--')
    else:
      return dict(color = cmap(0.5 * (coslat + 1)), ls = '-')
  
  def theta_changed(self, theta):
    '''
    
    '''
    
    self.planet.theta = theta * np.pi / 180
    self.update()

  def radius_changed(self, radius):
    '''
    
    '''
    
    self.occultor.r = radius
    self.update()

  def x0_changed(self, x0):
    '''
    
    '''
    
    self.planet.x0 = -x0
    self.update()

  def y0_changed(self, y0):
    '''
    
    '''
    
    self.planet.y0 = -y0
    self.update()
  
  def update(self):
    '''
  
    '''
    
    self.planet.compute_delta()
    ellipses = self.planet.ellipses
    vertices = self.planet.vertices
  
    # Set up plot
    self.ax.clear()
    self.ax.set_xlim(-2, 2)
    self.ax.set_ylim(-2, 2)
  
    # Plot the planet
    self.ax.plot(self.x, self.ylower, color = 'k', zorder = 98)
    self.ax.plot(self.x, self.yupper, color = 'k', zorder = 98)
  
    # Plot the occultor
    self.ax.plot(self.occultor.r * self.x - self.planet.x0, self.occultor.r * self.ylower - self.planet.y0, color = 'k', zorder = 99)
    self.ax.plot(self.occultor.r * self.x - self.planet.x0, self.occultor.r * self.yupper - self.planet.y0, color = 'k', zorder = 99)
    self.ax.fill_between(self.occultor.r * self.x - self.planet.x0, self.occultor.r * self.ylower - self.planet.y0, self.occultor.r * self.yupper - self.planet.y0, color = 'k', alpha = 0.1)
    
    # Plot the ellipses
    for ellipse in ellipses:

      # First identify and remove points behind the limb
      x = np.linspace(ellipse.x0 - ellipse.b, ellipse.x0 + ellipse.b, 1000)
      if self.planet.theta > 0:
        x[x < ellipse.x0 - ellipse.xlimb] = np.nan
      else:
        x[x > ellipse.x0 - ellipse.xlimb] = np.nan

      self.ax.plot(x - self.planet.x0, ellipse.val_lower(x) - self.planet.y0, **self.style(ellipse.latitude))
      self.ax.plot(x - self.planet.x0, ellipse.val_upper(x) - self.planet.y0, **self.style(ellipse.latitude))

    # Plot the vertices  
    for v in vertices:
      self.ax.plot(v[0] - self.planet.x0, v[1] - self.planet.y0, 'o', color = 'k', ms = 3)
    
    # Update the light curve
    self.flux = np.roll(self.flux, -1)
    self.flux[-1] = 1 + self.planet.delta_flux
    self.lc.set_ydata(self.flux)
    self.lcr.set_ydata(self.flux[-1])
    ymin = self.flux.min()
    ymax = self.flux.max()
    ypad = max(0.1, (ymax - ymin) * 0.1)
    self.axlc.set_ylim(ymin - ypad, ymax + ypad)
    self.lcl.set_text('%.3f' % self.flux[-1])
    self.lcl.set_y(self.flux[-1])
    
    # Re-draw!
    self.fig.canvas.draw_idle()