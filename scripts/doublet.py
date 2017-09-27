#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
doublet.py |github|
-------------------

An example of a PPO doublet in TRAPPIST-1. To see an animation, check out
:doc:`next_occultation.py </scripts/next_occultation>`.

  .. plot::
     :align: center
     
     from scripts import next_occultation
     next_occultation._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/doublet.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from planetplanet import Trappist1
from planetplanet.photo import GetAngles, DrawEyeball
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

def _test():
    '''

    '''

    pass

def plot():
    '''

    '''

    # Time of the doublet (I found this with `next_occultation.py`)
    t = 7918.73

    # Get the light curve up to that point plus a little bit
    system = Trappist1(sample = True, nbody = True, seed = 1234)
    b = system.bodies[1]
    c = system.bodies[2]
    time = np.arange(7900., t + 0.1, MINUTE)
    system.compute(time)
    
    # Get the indices of the occultation
    i0 = -300
    time = system.time[i0:]
    time -= time[0]
    time /= MINUTE
    #time -= (i0 // 2)
    flux = c.flux[i0:,-1]
    flux /= flux[0]
    c_x = c.x[i0:]
    c_y = c.y[i0:]
    c_z = c.z[i0:]
    c_vx = c.vx[i0:]
    c_vy = c.vy[i0:]
    c_vz = c.vz[i0:]
    b_x = b.x[i0:]
    b_y = b.y[i0:]
    b_vx = b.vx[i0:]
    b_vy = b.vy[i0:]
    
    # The figure for the paper
    fig, ax = pl.subplots(1, figsize = (7, 5))
    fig.subplots_adjust(top = 0.75, bottom = 0.25)
    
    # Plot the light curve
    ax.plot(time, flux, 'k-')
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax)
    ax.set_xticklabels([])
    ax.set_ylabel('Normalized Flux', fontweight = 'bold')
    ax.tick_params(direction = 'in')
    
    # Indices of ingress/egress
    ingress = 80
    egress = 230
    center = (ingress + egress) // 2
    inds = np.array(np.linspace(ingress, egress, 10, endpoint = True), dtype = int)
    topinds = inds[::2]
    botinds = inds[1::2]
    
    # Relative velocity vector normalization
    vnorm = 0.03 * np.mean(np.abs(c_vx - b_vx))
    
    # Plot the images
    x0 = 0.5102
    dx = 0.15
    pxinds = [x0 - 2 * dx, x0 - dx, x0, x0 + dx, x0 + 2 * dx]
    for px, i, j in zip(pxinds, topinds, botinds):
        
        # Top
        theta, gamma = GetAngles(c_x[i], c_y[i], c_z[i], c_vx[i], c_vy[i], c_vz[i])
        occ_dict = [dict(x = (b_x[i] - c_x[i]) / c._r, 
                         y = (b_y[i] - c_y[i]) / c._r, 
                         r = b._r / c._r, zorder = 99, alpha = 1)]
        DrawEyeball(px, 0.9, 0.025, theta = theta, nz = 31, gamma = gamma, 
                    draw_ellipses = False, radiancemap = c.radiancemap,
                    occultors = occ_dict, cmap = 'inferno', 
                    fig = fig, rasterize = True)
        ax.annotate("", xy = (time[i], ymax), 
                       xycoords = "data", xytext = (px, 0.9), 
                       textcoords = "figure fraction", 
                       clip_on = False, 
                       arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
        ax.annotate(r"$%d$" % time[i], xy = (px, 0.96), xycoords = "figure fraction",
                    ha = 'center', va = 'center', fontsize = 8)
        vx = (c_vx[i] - b_vx[i]) / vnorm
        if np.abs(vx) > 5:
            ax.annotate("", xy = (time[i], flux[i]), textcoords = 'offset points',
                        xytext = (-vx, 0), arrowprops = dict(arrowstyle = '<|-', color = 'k', alpha = 1))
        else:
            ax.plot(time[i], flux[i], 'ko', ms = 3)
            
        # Bottom
        theta, gamma = GetAngles(c_x[j], c_y[j], c_z[j], c_vx[j], c_vy[j], c_vz[j])
        occ_dict = [dict(x = (b_x[j] - c_x[j]) / c._r, 
                         y = (b_y[j] - c_y[j]) / c._r, 
                         r = b._r / c._r, zorder = 99, alpha = 1)]
        DrawEyeball(px, 0.1, 0.025, theta = theta, nz = 31, gamma = gamma, 
                    draw_ellipses = False, radiancemap = c.radiancemap,
                    occultors = occ_dict, cmap = 'inferno', 
                    fig = fig, rasterize = True)
        ax.annotate("", xy = (time[j], ymin), 
                       xycoords = "data", xytext = (px, 0.1), 
                       textcoords = "figure fraction", 
                       clip_on = False, 
                       arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
        ax.annotate(r"$%d$" % time[j], xy = (px, 0.04), xycoords = "figure fraction",
                    ha = 'center', va = 'center', fontsize = 8)
        vx = (c_vx[j] - b_vx[j]) / vnorm
        if np.abs(vx) > 5:
            ax.annotate("", xy = (time[j], flux[j]), textcoords = 'offset points',
                        xytext = (-vx, 0), arrowprops = dict(arrowstyle = '<|-', color = 'k', alpha = 1))
        else:
            ax.plot(time[j], flux[j], 'ko', ms = 3)
            
        ax.plot((time[i],time[i]), (ymax, flux[i]), color = 'k', alpha = 0.25, lw = 1)
        ax.plot((time[j],time[j]), (ymin, flux[j]), color = 'k', alpha = 0.25, lw = 1)

    fig.savefig('retro.pdf', bbox_inches = 'tight')
    
if __name__ == '__main__':
  plot()