#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball_orbit.py
----------------

Visualize an "eyeball" planet's orbit with and without hot spot offsets. 
See :py:mod:`planetplanet.photo.eyeball`. |github|

  .. plot::
     :align: center
     
     from planetplanet.photo.eyeball import DrawOrbit
     import matplotlib.pyplot as pl
     DrawOrbit(inc = 60., Omega = 0., ecc = 0.5, size = 2, figsize = (6, 6))
     pl.show()

  .. role:: raw-html(raw)
     :format: html
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/eyeball_orbit.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet.photo import eyeball, Star, Planet, System
import matplotlib.pyplot as pl
import numpy as np

if __name__ == '__main__':

  # Orbital params
  inc = 60.
  Omega = 0.
  ecc = 0.5
  w = 0

  # Compute two cases: one with no hotspot offset, one with an offset
  fig = [None, None, None, None]
  ax = [None, None, None, None]
  phase = [None, None]
  flux = [None, None]

  for i, Lambda, Phi in zip([0, 1], [0, 60], [0, 30]):

    # Plot the geometry
    fig[i], ax[i] = eyeball.DrawOrbit(inc = inc, Omega = Omega, ecc = ecc, w = w, size = 1.75, Lambda = Lambda, Phi = Phi,
                                      plot_phasecurve = False, label_phases = True, rasterize = True)

    # Compute the phase curve
    star = Star('A')
    b = Planet('b', per = 10., inc = inc, Omega = Omega, t0 = 0, ecc = ecc, w = w, 
               Phi = Phi, Lambda = Lambda, airless = True, phasecurve = True)
    system = System(star, b, mintheta = 0.001)
    time = np.linspace(-5, 5, 1000)
    system.compute(time)
    phase[i] = np.linspace(0, 1, len(b.time))
    flux[i] = np.array(b.flux[:,0])

  # Plot the first phasecurve
  fig[2], ax[2] = pl.subplots(1, figsize = (8, 2))
  fig[2].subplots_adjust(bottom = 0.3)
  ax[2].plot(phase[0], flux[0] / np.nanmax(flux[0]), 'k-')

  # Plot both phasecurves
  fig[3], ax[3] = pl.subplots(1, figsize = (8, 2))
  fig[3].subplots_adjust(bottom = 0.3)
  ax[3].plot(phase[0], flux[0] / np.nanmax(flux[1]), 'k-', alpha = 0.3)
  ax[3].plot(phase[1], flux[1] / np.nanmax(flux[1]), 'k-')

  # Adjust appearance
  for axis in [ax[2], ax[3]]:
    axis.set_xlabel('Orbital phase', fontweight = 'bold', fontsize = 12)
    axis.set_ylabel('Relative flux', fontweight = 'bold', fontsize = 12)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.margins(0.01, 0.2)

  fig[0].savefig('eyeball_orbit1.pdf', bbox_inches = 'tight', dpi = 600)
  fig[1].savefig('eyeball_orbit2.pdf', bbox_inches = 'tight', dpi = 600)
  fig[2].savefig('eyeball_phasecurve1.pdf', bbox_inches = 'tight')
  fig[3].savefig('eyeball_phasecurve2.pdf', bbox_inches = 'tight')

  # Show
  pl.show()