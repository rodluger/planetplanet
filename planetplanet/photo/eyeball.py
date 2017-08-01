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
                pos = None, draw_terminator = True, draw_outline = True, draw_ellipses = True):
  '''

  '''

  # If theta is in quadrants II and III, let's
  # shift it by 180 degrees in theta and alpha
  # so that we don't have to change the plotting
  # stuff below. The equations would still work, but
  # figuring out the z-order is a pain, so this is
  # simpler.
  if theta > np.pi / 2:
    theta = np.pi - theta
    gamma += np.pi
  elif theta < -np.pi / 2:
    theta = -np.pi - theta
    gamma += np.pi

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

def DrawOrbit(inc = 70., Omega = 0., ecc = 0., w = 0., dlambda = 0., dpsi = 0., nphases = 20, size = 1, 
              draw_orbit = True, draw_orbital_vectors = True, **kwargs):
  '''
  
  '''
  
  # Get the orbital elements over a full orbit of the planet
  # We are assuming a period of 10 days, but it doesn't matter for the plot!
  from . import Star, Planet, System
  star = Star('A')
  b = Planet('b', per = 10., inc = inc, Omega = Omega, t0 = 0, ecc = ecc, w = w, dlambda = dlambda, dpsi = dpsi)
  system = System(star, b)
  time = np.linspace(-5, 5, 1000)
  system.compute_orbits(time)

  # Plot the orbit
  fig, ax = pl.subplots(1, figsize = (8,8))
  if draw_orbit:
    ax.plot(b.x, b.y, 'k-', alpha = 0.5)

  # Adjust the plot range
  xmin = min(b.x.min(), b.y.min())
  xmax = max(b.x.max(), b.y.max())
  dx = xmax - xmin
  xmin -= 0.1 * dx
  xmax += 0.1 * dx
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(xmin, xmax)
  ax.axis('off')

  # Get the indices of the images we'll plot, sorted by zorder
  inds = np.array(list(range(0, 1000, 1000 // nphases)), dtype = int)
  inds = inds[np.argsort([-b.z[i] for i in inds])]

  # Plot images at different phases
  ax_eye = []
  for i in inds:
  
    # The position of the planet in the original frame (frame #0)
    x0 = b.x[i]
    y0 = b.y[i]
    z0 = b.z[i]
    r = np.sqrt(x0 ** 2 + y0 ** 2 + z0 ** 2)
    
    # Rotate so that the orbital plane is parallel to the x axis (frame #1)
    x1 = x0 * np.cos(b._Omega) + y0 * np.sin(b._Omega)
    y1 = y0 * np.cos(b._Omega) - x0 * np.sin(b._Omega)
    z1 = z0
  
    # Now rotate so that it's also parallel to the z axis (frame #2)
    x2 = x1
    y2 = 0
    z2 = z1 * np.sin(b._inc) + y1 * np.cos(b._inc)
  
    # Now rotate so that the substellar point faces the observer (frame #3)
    # This is a counter-clockwise rotation through an angle `xi`: 
    xi = np.pi / 2 - np.arctan2(z2, x2)
    x3 = 0
    y3 = 0
    z3 = z2 * np.cos(xi) + x2 * np.sin(xi)
    
    # Finally, translate the axis so that the planet is centered at the origin (frame #4)
    x4 = 0
    y4 = 0
    z4 = 0
  
    # Get the coordinates of the hot spot (frame #4)
    xh4 = 0
    yh4 = 0
    zh4 = -z3 * (b._r / r)
  
    # Add the offset in latitude. This is a *clockwise* rotation in the
    # zy plane by an angle `dlambda` about the planet center
    zh4, yh4 = zh4 * np.cos(b._dlambda) + yh4 * np.sin(b._dlambda), \
               yh4 * np.cos(b._dlambda) - zh4 * np.sin(b._dlambda)
  
    # Add the offset in longitude. This is a counterclockwise rotation
    # in the xz plane by an angle `dpsi` about the planet center
    xh4, zh4 = -zh4 * np.sin(b._dpsi), \
                zh4 * np.cos(b._dpsi)
  
    # Get the coordinates relative to the star
    xh3 = x3 + xh4
    yh3 = y3 + yh4
    zh3 = z3 + zh4
  
    # Rotate back to frame #2
    xh2 = xh3 * np.cos(xi) + zh3 * np.sin(xi)
    yh2 = yh3
    zh2 = zh3 * np.cos(xi) - xh3 * np.sin(xi)
  
    # Rotate back to frame #1
    xh1 = xh2
    yh1 = yh2 * np.sin(b._inc) + zh2 * np.cos(b._inc)
    zh1 = zh2 * np.sin(b._inc) - yh2 * np.cos(b._inc)
  
    # Finally, rotate back to the original frame (frame #0)
    xh0 = xh1 * np.cos(b._Omega) - yh1 * np.sin(b._Omega)
    yh0 = yh1 * np.cos(b._Omega) + xh1 * np.sin(b._Omega)
    zh0 = zh1

    # Compute the effective rotation angle
    gamma = np.arctan((yh0 - y0) / (xh0 - x0))
  
    # Find the x coordinate of the hotspot in the rotated frame (r)
    # relative to the planet center
    dxhr = (xh0 - x0) * np.cos(gamma) + (yh0 - y0) * np.sin(gamma)
  
    # Prevent tiny numerical errors from yielding NaNs in the arccos below
    if dxhr < -1: 
      dxhr = -1
    elif dxhr > 1: 
      dxhr = 1
  
    # Find the effective phase angle
    if (zh0 - z0) <= 0:
      # Quadrants I & II
      theta = np.pi - np.arccos(dxhr / b._r)
    else:
      # Quadrants III & IV
      theta = np.arccos(dxhr / b._r) - np.pi
    
    # Plot the radial vector
    if draw_orbital_vectors:
      ax.plot([0, x0], [0, y0], 'k-', alpha = 0.5, lw = 1)
  
    # Get the figure coordinates of the point
    disp_coords = ax.transData.transform((x0, y0))
    xf, yf = fig.transFigure.inverted().transform(disp_coords)
  
    # Draw the planet
    _, tmp, _, _ = DrawEyeball(theta = theta, gamma = gamma, fig = fig, pos = [xf, yf, 0.03 * size, 0.03 * size], **kwargs)
    ax_eye.append(tmp)
    
  return fig, ax, ax_eye
 
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
    dpsi = self.psi.val * np.pi / 180
    dlambda = self.lam.val * np.pi / 180
  
    # The coordinates of the substellar point
    x_ss = -np.cos(theta + dpsi) * np.cos(dlambda)
    y_ss = np.sin(dlambda)

    # The rotation angle that makes the planet symmetric about the x-axis
    gamma = -np.arctan2(y_ss, -x_ss)

    # Compute the new effective theta
    if theta + dpsi < -np.pi:
      theta = 2 * np.pi - np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
    elif theta + dpsi < 0:
      theta = -np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
    elif theta + dpsi > np.pi:
      theta = -np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
    else:
      theta = np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
    
    # Plot the planet
    DrawEyeball(fig = self.fig, theta = theta, gamma = gamma, **self.kwargs)