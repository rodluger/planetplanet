#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
mcmc_example.py |github|
------------------------

An example of how to perform joint modeling of multiple
occultations to recover certain system parameters from the
PPOs.

.. todo:: This script is not yet working. Stay tuned!

  .. plot::
     :align: center
     
     #from scripts import mcmc_example
     #import matplotlib.pyplot as pl
     #mcmc_example.plot()
     #pl.show()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/mcmc_example.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Trappist1
from planetplanet import LimbDarkenedMap
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

def plot():
  '''
  
  '''
  
  pass

def run():
  '''
  
  '''
  
  pass

def simulate():
  '''
  
  '''
  
  # Instantiate the system. Let's assume the planets are all
  # blackbodies with albedo = 0.3 (default for `Trappist1()`)
  # Let's also assume no limb darkening.
  system = Trappist1(sample = True, nbody = True, seed = 1234,
                     radiancemap = LimbDarkenedMap(),
                     limbdark = [0.])
  
  # NOTE: Let's make the planets massless so we can use the Kepler
  # solver in the modeling. This isn't realistic, but will speed up 
  # the code in this example.
  for planet in system.bodies[1:]:
    planet.m = 0.
  
  # Let's pick October 8, 2016 (t − 2,450,000 = 7670.)
  tstart = OCTOBER_08_2016
  
  # Find up to 100 PPOs of `c` by `b` in the next year
  times, _, durations = system.next_occultation('c', occultors = 'b', tstart = tstart,
                                                tend = tstart + 100, noccultations = 100)
  
  # Get all the occultations longer than 10 minutes
  times = np.array(times)[np.array(durations) > 10]
  
  # Get the light curve, this time using the Kepler solver, on a cadence
  # of five minutes and with an oversampling factor of 5.
  system = Trappist1(sample = True, nbody = False, seed = 1234,
                     radiancemap = LimbDarkenedMap(),
                     limbdark = [0.], oversample = 5)
  time = np.arange(tstart, tstart + 100, 5 * MINUTE)
  system.compute(time)
  
  # Plot each of the 6 occultations
  fig, ax = pl.subplots(3, 2)
  ax = ax.flatten()
  for i, axis in enumerate(ax):
    a = np.argmax(system.time > times[i] - 0.1)
    b = np.argmax(system.time > times[i] + 0.1)
    axis.plot(system.time[a:b], system.c.flux[a:b,-1] + system.c.total_flux[-1])
  pl.show()
  
if __name__ == '__main__':
  raise NotImplementedError("This script is not yet working!")