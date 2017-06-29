#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
daynight.py
-----------

Computes and plots a hypothetical mutual transit event, where two large 
planets transit the star and occult each other simultaneously.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo.ppo import Planet, Star, System
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import matplotlib.animation as animation
import numpy as np
np.random.seed(1234)
  
class Animation(object):
  '''
  An animation class for occultation movies.
  
  '''
  
  def __init__(self, t, fig, axim, axlc, pto, ptb, body, bodies, occultors, 
               time, flux, interval = 50, color = 'r', gifname = None, quiet = False):
    '''
    
    '''
    
    self.t = t
    self.fig = fig
    self.axim = axim
    self.pto = pto
    self.ptb = ptb
    self.body = body
    self.bodies = bodies
    self.occultors = occultors
    self.time = time
    self.flux = flux
    self.curve, = axlc.plot([], [], '-', color = color)
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
      
      # Curve
      self.curve.set_data(self.time[:self.t[j]], self.flux[:self.t[j]])
      
      # Body orbits
      for k, b in enumerate(self.bodies):
        self.ptb[k].set_xdata(b.x[self.t[j]])
        self.ptb[k].set_ydata(b.z[self.t[j]])

# Instantiate the star
star = Star('A', m = 0.1, r = 0.1, teff = 3000, nz = 31, color = 'k')

# Planet b
b = Planet('b', m = 1, per = 3, inc = 89.8, r = 3., t0 = 1.5, 
           nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

# Planet c
c = Planet('c', m = 1, per = 6, inc = 90., r = 3., t0 = 0., 
           nz = 51, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'b')

# System
system = System(star, b, c)
time = np.linspace(2.35, 2.4, 1000)

# Airless body
color = 'b'
gifname = '../img/occultation_limbdark.gif'
system.c.airless = False
system.compute(time, lambda2 = 15)
flux_airless = np.array(c.flux[:,-1])

# Stellar baseline
norm = np.nanmedian(star.flux[:,-1])
tmid = len(time) // 2

# Set up the figure
fig = pl.figure(figsize = (7, 8))
fig.subplots_adjust(left = 0.175)

# Plot the light curves
axlc = pl.subplot2grid((5, 3), (3, 0), colspan = 3, rowspan = 2)
axlc.plot(time, (norm + flux_airless) / norm, '-', color = color, label = 'Airless', alpha = 0.1)
axlc.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
axlc.set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 10)
axlc.get_yaxis().set_major_locator(MaxNLocator(4))
axlc.get_xaxis().set_major_locator(MaxNLocator(4))
for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
  tick.set_fontsize(8)
axlc.ticklabel_format(useOffset = False)
if system.time[0] > 1e4:
  for label in axlc.get_xmajorticklabels():
    label.set_rotation(30)

# Plot the orbits of all bodies
axxz = pl.subplot2grid((5, 3), (0, 0), colspan = 3, rowspan = 2)
f = np.linspace(0, 2 * np.pi, 1000)
for j, body in enumerate(system.bodies):
  if body == c:
    style = dict(color = 'r', alpha = 1, ls = '-', lw = 1)
  elif body == b:
    style = dict(color = 'k', alpha = 1, ls = '-', lw = 1)
  else:
    style = dict(color = 'k', alpha = 0.1, ls = '--', lw = 1)
  r = body.a * (1 - body.ecc ** 2) / (1 + body.ecc * np.cos(f))
  x = r * np.cos(body._w + f) - r * np.sin(body._w + f) * np.cos(body._inc) * np.sin(body._Omega)
  z = r * np.sin(body._w + f) * np.sin(body._inc)
  axxz.plot(x, z, **style)

# Plot the locations of the bodies
ptb = [None for body in system.bodies]
for bi, body in enumerate(system.bodies):
  if body == c:
    ptb[bi], = axxz.plot(body.x[0], body.z[0], 'o', color = 'r', alpha = 1, markeredgecolor = 'k', zorder = 99)
  elif body == b:
    ptb[bi], = axxz.plot(body.x[0], body.z[0], 'o', color = 'lightgrey', alpha = 1, markeredgecolor = 'k', zorder = 99)
  else:
    ptb[bi], = axxz.plot(body.x[0], body.z[0], 'o', color = '#dddddd', alpha = 1, markeredgecolor = '#999999', zorder = 99)
  
# Appearance
axxz.set_ylim(-max(np.abs(axxz.get_ylim())), max(np.abs(axxz.get_ylim())))
axxz.set_xlim(-max(np.abs(axxz.get_xlim())), max(np.abs(axxz.get_xlim())))
axxz.set_aspect('equal')
axxz.axis('off')

# Plot the image
axim = pl.subplot2grid((5, 3), (2, 0), colspan = 3, rowspan = 1)
_, pto = system.plot_image(0, c, ax = axim, occultors = [1])
axim.set_xlim(-15, 15)
axim.axis('off')
axim.set_aspect('equal')

# Animate!
animation = Animation(range(len(time)), fig, axim, axlc, pto, ptb, c, system.bodies, [b], time, (norm + flux_airless) / norm, gifname = gifname, color = color)
system._animations.append(animation)
pl.show()