#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
timestep.py |github|
--------------------

Tests the minimum timestep we need in the N-body code.
We integrate for one year and look at the errors as a fraction
of the planet radius. Looks like a timestep of 1 hour leads to
negligible (< 1 percent) error over 1 year.

  .. plot::
     :align: center
     
     from scripts import timestep
     import matplotlib.pyplot as pl
     timestep.plot()
     pl.show()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/timestep.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Trappist1
import matplotlib.pyplot as pl
import numpy as np

def plot():
  '''
  
  '''
  
  fig, ax = pl.subplots(3, figsize = (6,8))
  fig.subplots_adjust(left = 0.2)
  time = np.linspace(0, 365, 10000)

  # Tiny timestep (1 minute)
  np.random.seed(1234)
  system1 = Trappist1(nbody = True, timestep = 1. / 1440)
  system1.compute(time)

  # Large timestep (1 hour)
  np.random.seed(1234)
  system2 = Trappist1(nbody = True, timestep = 1. / 24)
  system2.compute(time)

  for body1, body2 in zip(system1.bodies[1:], system2.bodies[1:]):
    ax[0].plot(system1.time, 100 * (body1.x - body2.x) / body1._r, color = body1.color)
    ax[1].plot(system1.time, 100 * (body1.y - body2.y) / body1._r, color = body1.color)
    ax[2].plot(system1.time, 100 * (body1.z - body2.z) / body1._r, color = body1.color)

  ax[0].set_xticklabels([])
  ax[1].set_xticklabels([])
  ax[2].set_xlabel('Time [days]', fontsize = 16, fontweight = 'bold')
  ax[0].set_ylabel('x error (% of radius)', fontsize = 8, fontweight = 'bold')
  ax[1].set_ylabel('y error (% of radius)', fontsize = 8, fontweight = 'bold')
  ax[2].set_ylabel('z error (% of radius)', fontsize = 8, fontweight = 'bold')

  return fig, ax

if __name__ == '__main__':
  fig, ax = plot()
  pl.show()