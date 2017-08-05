#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball.py
----------

'''

import numpy as np
np.seterr(invalid = 'ignore')
import matplotlib.pyplot as pl
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.widgets import Slider

def ZenithColor(z, cmap = 'inferno'):
  '''
  Color of a given zenith angle slice
  
  '''
  
  return pl.get_cmap(cmap)(0.2 + 0.4 * (np.cos(z) + 1))

def DrawEyeball(x0 = 0, y0 = 0, r = 1, theta = np.pi / 3, nz = 11, gamma = 0, occultors = [], cmap = 'inferno', fig = None, 
                pos = None, draw_terminator = True, draw_outline = True, draw_ellipses = True, rasterize = False):
  '''

  '''
  
  # The rotation transformation, Equation (E6) in the paper
  xy = lambda x, y: (x * np.cos(gamma) + y * np.sin(gamma), y * np.cos(gamma) - x * np.sin(gamma))
  
  # Set up the floating axis
  if fig is None:
    fig = pl.figure(figsize = (6,6))
  tr = Affine2D().rotate_deg(gamma * 180 / np.pi)
  x = 1. / (np.abs(np.cos(gamma)) + np.abs(np.sin(gamma)))
  scale = max([r] + [occultor['r'] for occultor in occultors])
  x *= scale
  grid_helper = floating_axes.GridHelperCurveLinear(tr, extremes=(-x, x, -x, x))
  ax_orig = floating_axes.FloatingSubplot(fig, 111, grid_helper = grid_helper)
  if pos is not None:
    pos[0] -= pos[2] / 2
    pos[1] -= pos[3] / 2
    ax_orig.set_position(pos)
  ax_orig.axis["bottom"].set_visible(False)
  ax_orig.axis["top"].set_visible(False)
  ax_orig.axis["left"].set_visible(False)
  ax_orig.axis["right"].set_visible(False)
  ax_orig.patch.set_alpha(0)
  fig.add_subplot(ax_orig)
  ax = ax_orig.get_aux_axes(tr)
  if rasterize:
    ax_orig.set_rasterization_zorder(9999)
    ax.set_rasterization_zorder(9999)
  
  # Plot the occultors. Note that we need to transform
  # their position vectors since we're in a rotated frame.
  occ = [None for i in occultors]
  for i, occultor in enumerate(occultors):
    xo = occultor['x']
    yo = occultor['y']
    ro = occultor['r']
    zo = occultor.get('zorder', 1)
    co = occultor.get('color', 'lightgray')
    ao = occultor.get('alpha', 1)
    xo, yo = xy(xo, yo)
    occ[i] = pl.Circle((xo, yo), ro, color = co, ec = 'k' if draw_outline else 'none', alpha = ao, zorder = zo, clip_on = False)
    ax.add_artist(occ[i])
    
  # Plot the occulted body
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  if draw_outline:
    ax.plot(x0 + x, y0 + y, color = 'k', zorder = 0, lw = 1, clip_on = False)
    ax.plot(x0 + x, y0 - y, color = 'k', zorder = 0, lw = 1, clip_on = False)

  # Plot the zenith angle ellipses
  zarr = np.linspace(0, np.pi, nz + 2)
  for i, z in enumerate(zarr[1:]):

    # The ellipse
    a = r * np.abs(np.sin(z))
    b = max(0.001, a * np.abs(np.sin(theta)))
    xE = -r * np.cos(z) * np.cos(theta)
    xlimb = r * np.cos(z) * np.sin(theta) * np.tan(theta)
    if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
      xmin = xE - b
    else:
      xmin = xE - xlimb
    if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
      xmax = xE + b
    else:
      xmax = xE - xlimb
        
    # Plot it
    x = np.linspace(xE - b, xE + b, 1000)
    if theta > 0:
      x[x < xE - xlimb] = np.nan
    elif theta > -np.pi / 2:
      x[x > xE - xlimb] = np.nan
    A = b ** 2 - (x - xE) ** 2
    A[A < 0] = 0
    y = (a / b) * np.sqrt(A)
    if np.abs(np.cos(z)) < 1e-5:
      # This is the terminator
      if draw_terminator and scale < 3 * r:
        style = dict(color = 'k', ls = '--', lw = 0.5, zorder = 0, clip_on = False)
        ax.plot(x0 + x, y0 + y, **style)
        ax.plot(x0 + x, y0 - y, **style)
    else:
      # These are the ellipse boundaries
      if draw_ellipses and scale < 3 * r:
        style = dict(color = 'k', ls = '-', lw = 0.5, zorder = 0, clip_on = False)
        ax.plot(x0 + x, y0 + y, **style)
        ax.plot(x0 + x, y0 - y, **style)
      
    # Fill the ellipses
    if theta < 0:
      ax.fill_between(x0 + x, y0 - y, y0 + y, color = ZenithColor(zarr[i+1], cmap = cmap), zorder = 0.5 * (z / np.pi - 1), clip_on = False)
    else:
      ax.fill_between(x0 + x, y0 - y, y0 + y, color = ZenithColor(zarr[i], cmap = cmap), zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
  
    # Fill the ellipses that are cut off by the limb
    if theta < 0:
      x_ = np.linspace(-r, xE - xlimb, 1000)
      y_ = np.sqrt(r ** 2 - x_ ** 2)
      ax.fill_between(x0 + x_, y0 - y_, y0 + y_, color = ZenithColor(zarr[i], cmap = cmap), zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
    else:
      x_ = np.linspace(-r, xE - xlimb, 1000)
      y_ = np.sqrt(r ** 2 - x_ ** 2)
      ax.fill_between(x0 + x_, y0 - y_, y0 + y_, color = ZenithColor(zarr[i], cmap = cmap), zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
  
  return fig, ax, occ, xy

def DrawOrbit(inc = 70., Omega = 0., ecc = 0., w = 0., Phi = 0., Lambda = 0., nphases = 20, size = 1, 
              draw_orbit = True, draw_orbital_vectors = True, plot_phasecurve = False, 
              label_phases = False, **kwargs):
  '''
  
  '''
      
  # Get the orbital elements over a full orbit of the planet
  # We are assuming a period of 10 days, but it doesn't matter for the plot!
  from . import Star, Planet, System
  star = Star('A')
  b = Planet('b', per = 10., inc = inc, Omega = Omega, t0 = 0, ecc = ecc, w = w, 
             Phi = Phi, Lambda = Lambda, airless = True, phasecurve = True)
  system = System(star, b, mintheta = 0.001)
  time = np.linspace(-5, 5, 1000)
  if plot_phasecurve:
    system.compute(time)
  else:
    system.compute_orbits(time)

  # Plot stuff
  fig, ax = pl.subplots(1, figsize = (8,8))
  ax.margins(0.1, 0.1)
  
  # Phase curve
  if plot_phasecurve:
    figphase, axphase = pl.subplots(1, figsize = (8, 2))
    figphase.subplots_adjust(bottom = 0.3)
    axphase.plot(np.linspace(0, 1, len(b.time)), b.flux[:,0] / (np.nanmax(b.flux[:,0])), 'k-')
    axphase.set_xlabel('Orbital phase', fontweight = 'bold', fontsize = 12)
    axphase.set_ylabel('Relative flux', fontweight = 'bold', fontsize = 12)
      
  # Orbit outline
  if draw_orbit:
    ax.plot(b.x, b.y, 'k-', alpha = 0.5)
  
  # Adjust the figure dimensions so the aspect ratio is unity
  left = 0.125
  right = 0.9
  xmin, xmax = ax.get_xlim()
  ymin, ymax = ax.get_ylim()
  width = right - left
  height = width * (ymin - ymax) / (xmin - xmax)
  bottom = 0.5 - height / 2
  top = 0.5 + height / 2
  fig.subplots_adjust(left = left, right = right, bottom = bottom, top = top)
  ax.axis('off')

  # Get the indices of the images we'll plot, sorted by zorder
  inds = np.array(list(range(0, 1000, 1000 // nphases)), dtype = int)
  inds = inds[np.argsort([-b.z[i] for i in inds])]

  # Plot images at different phases
  axes = [ax]
  for i in inds:
    
    # The position of the planet in the rotated sky plane
    x = b.x[i] * np.cos(b._Omega) + b.y[i] * np.sin(b._Omega)
    y = b.y[i] * np.cos(b._Omega) - b.x[i] * np.sin(b._Omega)
    z = b.z[i]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        
    # Coordinates of the hotspot in a frame where the planet is
    # at x, y, z = (0, 0, r), at full phase
    xprime = b._r * np.cos(b._Phi) * np.sin(b._Lambda)
    yprime = b._r * np.sin(b._Phi)
    zprime = r - b._r * np.cos(b._Phi) * np.cos(b._Lambda)

    # Transform to the rotated sky plane
    rxz = np.sqrt(x ** 2 + z ** 2)
    xstar = ((z * r) * xprime - (x * y) * yprime + (x * rxz) * zprime) / (r * rxz)
    ystar = (rxz * yprime + y * zprime) / r
    zstar = (-(x * r) * xprime - (y * z) * yprime + (z * rxz) * zprime) / (r * rxz)
        
    # Transform back to the true sky plane
    xstar, ystar = xstar * np.cos(b._Omega) - ystar * np.sin(b._Omega), \
                   ystar * np.cos(b._Omega) + xstar * np.sin(b._Omega)
    x = b.x[i]
    y = b.y[i]
    
    # Distance from planet center to hotspot
    d = np.sqrt((xstar - x) ** 2 + (ystar - y) ** 2)
    
    # Get the rotation and phase angles
    gamma = np.arctan2(ystar - y, xstar - x) + np.pi
    if (zstar - z) <= 0:
      theta = np.arccos(d / b._r)
    else:
      theta = -np.arccos(d / b._r)
    
    # Plot the radial vector
    if draw_orbital_vectors:
      ax.plot([0, x], [0, y], 'k-', alpha = 0.5, lw = 1)
  
    # Get the figure coordinates of the point
    disp_coords = ax.transData.transform((x, y))
    xf, yf = fig.transFigure.inverted().transform(disp_coords)
  
    # Draw the planet
    _, tmp, _, _ = DrawEyeball(theta = theta, gamma = gamma, fig = fig, pos = [xf, yf, 0.03 * size, 0.03 * size], **kwargs)
    axes.append(tmp)
    
    # Indicate the orbital phase
    if label_phases:
      dx = x / r
      dy = y / r
      dr = np.sqrt(dx ** 2 + dy ** 2)
      dx *= 16 * size / dr
      dy *= 16 * size / dr
      tmp.annotate("%.2f" % (i / 1000.), xy = (0, 0), xytext = (dx, dy), xycoords = 'data', 
                   textcoords = 'offset points', fontsize = 8, ha = 'center', va = 'center',
                   zorder = 10000)

  if plot_phasecurve:
    return fig, axes, figphase, axphase
  else:
    return fig, axes
 
class Interact(object):
  '''
  
  '''
  
  def __init__(self, **kwargs):
    '''
  
    '''
  
    self.kwargs = kwargs
    self.fig = pl.figure(figsize = (6,6))
    self.fig.subplots_adjust(bottom = 0.25)
    self.axtheta = pl.axes([0.3, 0.05, 0.44, 0.03])
    self.theta = Slider(self.axtheta, r'$\theta$', -180., 180., valinit = 90.)
    self.axpsi = pl.axes([0.3, 0.1, 0.44, 0.03])
    self.psi = Slider(self.axpsi, r'$\psi$', -90, 90., valinit = 0.)
    self.axlam = pl.axes([0.3, 0.15, 0.44, 0.03])
    self.lam = Slider(self.axlam, r'$\lambda$', -90, 90., valinit = 0.)
    self.theta.on_changed(self.update)
    self.psi.on_changed(self.update)
    self.lam.on_changed(self.update)
    self.update(90.)
    pl.show()
  
  def update(self, val):
    '''
  
    '''
  
    # Remove all axes except sliders
    for ax in self.fig.get_axes():
      if ax not in [self.axtheta, self.axpsi, self.axlam]:
        ax.remove()
  
    # Get the angles
    theta = self.theta.val * np.pi / 180
    Lambda = self.psi.val * np.pi / 180
    Phi = self.lam.val * np.pi / 180
  
    # The coordinates of the substellar point
    x_ss = -np.cos(theta + Lambda) * np.cos(Phi)
    y_ss = np.sin(Phi)

    # The rotation angle that makes the planet symmetric about the x-axis
    gamma = -np.arctan2(y_ss, -x_ss)

    # Compute the new effective theta
    if theta + Lambda < -np.pi:
      theta = 2 * np.pi - np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
    elif theta + Lambda < 0:
      theta = -np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
    elif theta + Lambda > np.pi:
      theta = -np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
    else:
      theta = np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
    
    # Plot the planet
    DrawEyeball(fig = self.fig, theta = theta, gamma = gamma, **self.kwargs)