#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
geometry.py
-----------

Not yet ready!

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np

def Observer(ax, r = 1, theta = np.pi / 8, lat = np.pi / 4):
  '''
  Observer view
  
  '''
  
  # Plot the planet
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x, y, color = 'k', lw = 1, zorder = 100)
  ax.plot(x, -y, color = 'k', lw = 1, zorder = 100)

  # Compute the ellipse
  a = r * np.abs(np.sin(lat))
  b = a * np.abs(np.sin(theta))
  x0 = -r * np.cos(lat) * np.cos(theta)
  y0 = 0
  xlimb = r * np.cos(lat) * np.sin(theta) * np.tan(theta)
  if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
    xmin = x0 - b
  else:
    xmin = x0 - xlimb
  if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
    xmax = x0 + b
  else:
    xmax = x0 - xlimb
    
  # Plot it: invisble
  x = np.linspace(x0 - b, x0 + b, 1000)
  if theta > 0:
    x[x > x0 - xlimb] = np.nan
  else:
    x[x < x0 - xlimb] = np.nan
  A = b ** 2 - (x - x0) ** 2
  A[A < 0] = 0
  y = (a / b) * np.sqrt(A)
  ax.plot(x, y, color = 'k', ls = '--', lw = 1)
  ax.plot(x, -y, color = 'k', ls = '--', lw = 1)

  # Plot it: visible
  x = np.linspace(x0 - b, x0 + b, 1000)
  if theta > 0:
    x[x < x0 - xlimb] = np.nan
  else:
    x[x > x0 - xlimb] = np.nan
  A = b ** 2 - (x - x0) ** 2
  A[A < 0] = 0
  y = (a / b) * np.sqrt(A)
  ax.plot(x, y, color = 'k', ls = '-', lw = 1)
  ax.plot(x, -y, color = 'k', ls = '-', lw = 1)

  # Ellipse axes
  ax.plot([x0, x0], [0, a], color = 'k', lw = 1, zorder = 100)
  ax.plot([x0, x0 + b], [0, 0], color = 'k', lw = 1, zorder = 100)
  ax.annotate(r'$a$', xy = (x0, a / 2), ha = 'right', va = 'center', xytext = (-10, 0), textcoords = 'offset points', fontsize = 16)
  ax.annotate(r'$b$', xy = (x0 + b / 2, 0), ha = 'center', va = 'top', xytext = (0, -10), textcoords = 'offset points', fontsize = 16)
  ax.plot(x0, y0, 'ko')
  ax.annotate(r'$(x_0, y_0)$', xy = (x0, y0), ha = 'right', va = 'top', xytext = (-5, -5), textcoords = 'offset points', fontsize = 12)
  ylimb = (a / b) * np.sqrt(b ** 2 - (x0 - xlimb - x0) ** 2)
  ax.plot(x0 - xlimb, -ylimb, 'ko')
  ax.plot(x0 - xlimb, ylimb, 'ko')
  ax.annotate(r'$(x_\mathrm{limb}, y_0 - y_\mathrm{limb})$', xy = (x0 - xlimb, -ylimb), ha = 'right', 
              va = 'top', xytext = (-5, -5), textcoords = 'offset points', fontsize = 12)
  ax.annotate(r'$(x_\mathrm{limb}, y_0 + y_\mathrm{limb})$', xy = (x0 - xlimb, ylimb), ha = 'right', 
              va = 'bottom', xytext = (-5, 5), textcoords = 'offset points', fontsize = 12)
  ax.plot(0, 0, 'ko')
  ax.annotate(r'$(0, 0)$', xy = (0, 0), ha = 'left', va = 'top', xytext = (5, -5), textcoords = 'offset points', fontsize = 12)

  # Fill
  ax.fill_between(x, -y, y, color = 'lightgray', zorder = 0)
  x = np.linspace(-r, x0 - xlimb, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.fill_between(x, -y, y, color = 'lightgray', zorder = 0)
    
  # Appearance
  ax.set_aspect('equal')
  ax.set_frame_on(False)
  ax.set_xticks([])
  ax.set_yticks([])

def Side(ax, r = 1, theta = np.pi / 8, lat = np.pi / 4):
  '''
  Side view
  
  '''
  
  # Plot the planet
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x, y, color = 'k', lw = 1, zorder = 100)
  ax.plot(x, -y, color = 'k', lw = 1, zorder = 100)

  # Plot the ellipse
  d = r * np.cos(lat)
  a = r * np.sin(lat)
  
  x = np.linspace(-r, -d, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.fill_between(x, -y, y, color = 'lightgray', zorder = 0)
  
  ax.plot([-d, 0], [0, 0], 'k-', lw = 1)
  ax.plot([-d, 0], [a, 0], 'k-', lw = 1)
  ax.plot([-d, -d], [a, -a], 'k-', lw = 1)
  ax.plot(0, 0, 'ko')
  ax.annotate(r'$(0, 0)$', xy = (0, 0), ha = 'left', va = 'top', xytext = (5, -5), textcoords = 'offset points', fontsize = 12)
  ax.annotate(r'$a$', xy = (-d, a / 2), ha = 'right', va = 'top', xytext = (-10, 0), textcoords = 'offset points', fontsize = 16)
  ax.annotate(r'$d$', xy = (-d / 2, 0), ha = 'center', va = 'top', xytext = (0, -10), textcoords = 'offset points', fontsize = 16)
  ax.annotate(r'$r$', xy = (-d / 2, a / 2), ha = 'center', va = 'top', xytext = (15, 15), textcoords = 'offset points', fontsize = 16)
  x = np.linspace(-0.1, -0.0707107, 10000)
  y = np.sqrt(0.01 - x ** 2)
  ax.plot(x, y, 'k-', lw = 1)
  l = ax.annotate(r'l', xy = (-0.15, 0.1), ha = 'center', va = 'top', xytext = (0, 0), 
                  textcoords = 'offset points', fontsize = 18, family='PilGi')
  
  ax.annotate("To star", xy = (-r, 0), xytext = (-70, 0),
              ha = 'center', va = 'center', annotation_clip = False, color = 'k',
              textcoords = "offset points", arrowprops=dict(arrowstyle = "<|-", color = 'k'))
  
  # Appearance
  ax.set_aspect('equal')
  ax.set_frame_on(False)
  ax.set_xticks([])
  ax.set_yticks([])

def Top(ax, r = 1, theta = np.pi / 8, lat = np.pi / 4):
  '''
  Top view
  
  '''
  
  # Plot the planet
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x, y, color = 'k', lw = 1, zorder = 100)
  ax.plot(x, -y, color = 'k', lw = 1, zorder = 100)

  # Compute the ellipse
  a = r * np.abs(np.sin(lat))
  b = a * np.abs(np.sin(theta))
  d = r * np.cos(lat)
  x0 = -r * np.cos(lat) * np.cos(theta)
  y0 = 0
  xlimb = r * np.cos(lat) * np.sin(theta) * np.tan(theta)
  
  
  # The ellipse (a line from this vantage point)
  ax.plot([x0 - b, x0 - xlimb], [np.sqrt(r ** 2 - (x0 - b) ** 2), 0], 'k--', lw = 1)
  ax.plot([x0 - xlimb, x0 + b], [0, -np.sqrt(r ** 2 - (x0 + b) ** 2)], 'k-', lw = 1)
  
  # To observer
  ax.plot([0, 0], [0, -r], 'k-', lw = 1)
  
  # To star
  ax.plot([0, -r * np.cos(theta)], [0, -r * np.sin(theta)], 'k-', lw = 1)
  
  # Terminator
  ax.plot([0, r * np.sin(theta)], [0, -r * np.cos(theta)], 'k-', lw = 1)
  
  # Fill it
  x = np.linspace(-r, x0 + b, 1000)
  x1, y1 = (x0 - b, np.sqrt(r ** 2 - (x0 - b) ** 2))
  x2, y2 = (x0 + b, -np.sqrt(r ** 2 - (x0 + b) ** 2))  
  m = (y2 - y1) / (x2 - x1)
  c = y1 - m * x1
  ytop = m * x + c
  ybot = -np.sqrt(r ** 2 - x ** 2)
  ax.fill_between(x, ybot, ytop, color = 'lightgray', zorder = 0, where = (ytop < np.sqrt(r ** 2 - x ** 2)))
  ax.fill_between(x, -ybot, ybot, color = 'lightgray', zorder = 0, where = (x < x0 - b))

  # Origin
  ax.plot(0, 0, 'ko')
  ax.annotate(r'$(0, 0)$', xy = (0, 0), ha = 'left', va = 'bottom', xytext = (5, 5), textcoords = 'offset points', fontsize = 12)
  
  # x0
  ax.plot([x0, 0], [0, 0], 'k-', lw = 1)
  ax.plot([x0, x0], [0, -r * np.cos(lat) * np.sin(theta)], 'k-', lw = 1)
  
  # limb
  ax.plot([x0 - xlimb, x0], [0, 0], 'r-', lw = 1)
  ax.plot(x0 - xlimb, 0, 'ro')
  
  # a, b
  ax.plot([x0 + b, x0 + b], [-r * np.cos(lat) * np.sin(theta), -np.sqrt(r ** 2 - (x0 + b) ** 2)], 'k-', lw = 1)
  ax.plot([x0, x0 + b], [-r * np.cos(lat) * np.sin(theta), -r * np.cos(lat) * np.sin(theta)], 'k-', lw = 1)
  
  # Labels  
  ax.annotate("To observer", xy = (0, -0.99 * r), xytext = (0, -50),
              ha = 'center', va = 'center', annotation_clip = False, color = 'k',
              textcoords = "offset points", arrowprops=dict(arrowstyle = "<|-", color = 'k'))  
  
  ax.annotate("To star", xy = (-r * np.cos(theta), -r * np.sin(theta)), xytext = (-66 * np.cos(theta), -66 * np.sin(theta)),
              ha = 'center', va = 'center', annotation_clip = False, color = 'k',
              textcoords = "offset points", arrowprops=dict(arrowstyle = "<|-", color = 'k'))
  
  ax.annotate(r'$d$', xy = (-0.3375, -0.18), ha = 'center', va = 'center', xytext = (2, -2), textcoords = 'offset points', fontsize = 14)
  ax.annotate(r'$x_0$', xy = (x0 / 2, 0), ha = 'center', va = 'center', xytext = (0, 10), textcoords = 'offset points', fontsize = 14)
  ax.annotate(r'$b$', xy = (x0 + b / 2, -0.32), ha = 'center', va = 'center', xytext = (0, -2), textcoords = 'offset points', fontsize = 14)
  ax.annotate(r'$x_\mathrm{limb}$', xy = (-0.707, 0.043), ha = 'center', va = 'center', xytext = (0, 5), textcoords = 'offset points', fontsize = 12, color = 'r')
  ax.annotate(r'$r$', xy = (0.28, -0.5), ha = 'center', va = 'center', xytext = (0, 0), textcoords = 'offset points', fontsize = 14)

  
  ax.annotate(r'$\theta$', xy = (0.0464, -0.23), ha = 'center', va = 'center', xytext = (0, -5), textcoords = 'offset points', fontsize = 14)
  ax.annotate(r'$\theta$', xy = (-0.225, -0.048), ha = 'center', va = 'center', xytext = (0, 0), textcoords = 'offset points', fontsize = 12)
  ax.annotate(r'$\theta$', xy = (-0.6786, -0.1286), ha = 'center', va = 'center', xytext = (0, 0), textcoords = 'offset points', fontsize = 10)

  
  # Angles
  x = np.linspace(0, 0.066, 10000)
  y = -np.sqrt(0.03 - x ** 2)
  ax.plot(x, y, 'k-', lw = 1)
  
  x = np.linspace(-0.1732, -0.16, 10000)
  y = -np.sqrt(0.03 - x ** 2)
  ax.plot(x, y, 'k-', lw = 1)
  
  x = np.linspace(-0.6882, x0, 10000)
  y = -0.28 + np.sqrt(0.01 - (x - x0) ** 2)
  ax.plot(x, y, 'k-', lw = 1)
  
  ax.plot([x0 - 0.03, x0], [-0.03, -0.03], 'r-', lw = 1) 
  ax.plot([x0 - 0.03, x0 - 0.03], [-0.03, 0], 'r-', lw = 1) 
  
  ax.plot([x0 + b - 0.03, x0 + b], [-r * np.cos(lat) * np.sin(theta) - 0.03, -r * np.cos(lat) * np.sin(theta) - 0.03], 'k-', lw = 1) 
  ax.plot([x0 + b - 0.03, x0 + b - 0.03], [-r * np.cos(lat) * np.sin(theta) - 0.03, -r * np.cos(lat) * np.sin(theta)], 'k-', lw = 1) 
    
  # Appearance
  ax.set_aspect('equal')
  ax.set_frame_on(False)
  ax.set_xticks([])
  ax.set_yticks([])

def Front(ax, r = 1, theta = np.pi / 8, lat = np.pi / 4):
  '''
  Front view
  
  '''
  
  # Plot the planet
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x, y, color = 'k', lw = 1, zorder = 100)
  ax.plot(x, -y, color = 'k', lw = 1, zorder = 100)

  # Compute the ellipse
  a = r * np.abs(np.sin(lat))
  b = a * np.abs(np.sin(theta))
  
  # Plot it
  x = np.linspace(-a, a, 1000)
  y = np.sqrt(a ** 2 - x ** 2)
  ax.plot(x, -y, 'k-')
  ax.plot(x, y, 'k-')
  ax.fill_between(x, -y, y, color = 'lightgray', zorder = 0)
  
  ax.plot([0, 0], [0, a], 'k-', lw = 1)
  ax.plot([-r, 0], [0, 0], 'k-', lw = 1)
  ax.annotate(r'$a$', xy = (0, a / 2), ha = 'left', va = 'center', xytext = (5, 0), textcoords = 'offset points', fontsize = 16)
  ax.annotate(r'$r$', xy = (-r / 2, 0), ha = 'center', va = 'top', xytext = (0, -5), textcoords = 'offset points', fontsize = 16)

  ax.plot(0,0,'ko')
  ax.annotate(r'$(0, 0)$', xy = (0, 0), ha = 'left', va = 'top', xytext = (5, -5), textcoords = 'offset points', fontsize = 12)
  
  # Appearance
  ax.set_aspect('equal')
  ax.set_frame_on(False)
  ax.set_xticks([])
  ax.set_yticks([])

fig, ax = pl.subplots(2,2, figsize = (16, 16))
Observer(ax[0,0]); ax[0,0].set_title('Observer', fontweight = 'bold', fontsize = 18)
Side(ax[0,1]); ax[0,1].set_title('Side', fontweight = 'bold', fontsize = 18)
Top(ax[1,0]); ax[1,0].set_title('Top', fontweight = 'bold', fontsize = 18)
Front(ax[1,1]); ax[1,1].set_title('Front', fontweight = 'bold', fontsize = 18)

fig.savefig('geometry.pdf', bbox_inches = 'tight')