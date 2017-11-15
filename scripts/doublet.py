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
    time -= 163
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
    fig, ax = pl.subplots(1, figsize = (6, 4))
    fig.subplots_adjust(top = 0.75, bottom = 0.15, left = 0.1, right = 0.95)

    # Plot the light curve
    ax.plot(time, flux, 'k-')
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('Time [minutes]', fontweight = 'bold')
    ax.set_ylabel('Normalized Flux', fontweight = 'bold')
    ax.tick_params(direction = 'in')

    # Indices of ingress/egress
    ingress = 80
    egress = 230
    center = (ingress + egress) // 2
    inds = np.array(np.linspace(ingress, egress, 7, endpoint = True), dtype = int)

    # Relative velocity vector normalization
    vnorm = 0.03 * np.mean(np.abs(c_vx - b_vx))

    # Plot the images
    x0 = 0.5102
    dx = 0.14
    pxinds = [x0 - 3 * dx, x0 - 2 * dx, x0 - dx, x0, x0 + dx, x0 + 2 * dx, x0 + 3 * dx]
    for px, i in zip(pxinds, inds):

        # Top
        theta, gamma = GetAngles(c_x[i], c_y[i], c_z[i], c_vx[i], c_vy[i], c_vz[i])
        occ_dict = [dict(x = (b_x[i] - c_x[i]) / c._r,
                         y = (b_y[i] - c_y[i]) / c._r,
                         r = b._r / c._r, zorder = 99, alpha = 1,
                         lw = 1)]
        DrawEyeball(px, 0.9, 0.03, theta = theta, nz = 31, gamma = gamma,
                    draw_ellipses = False, radiancemap = c.radiancemap,
                    occultors = occ_dict, cmap = 'inferno',
                    fig = fig, rasterize = True, lw = 1, color = 'k')
        ax.annotate("", xy = (time[i], ymax),
                       xycoords = "data", xytext = (px, 0.9),
                       textcoords = "figure fraction",
                       clip_on = False,
                       arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))
        vx = (c_vx[i] - b_vx[i]) / vnorm
        if np.abs(vx) > 5:
            ax.annotate("", xy = (time[i], flux[i]), textcoords = 'offset points',
                        xytext = (-vx, 0), arrowprops = dict(arrowstyle = '<|-', color = 'k', alpha = 1))
        else:
            ax.plot(time[i], flux[i], 'ko', ms = 3)

        ax.plot((time[i], time[i]), (ymax, flux[i]), color = 'k', alpha = 0.25, lw = 1)

    fig.savefig('retro.pdf', dpi = 800)
    pl.show()

if __name__ == '__main__':
  plot()
