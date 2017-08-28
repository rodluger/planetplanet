#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
next_occultation.py |github|
----------------------------

Compute the time of the next occultation of a given planet and plot its light curve. 

  .. plot::
     :align: center
     
     from scripts import next_occultation
     next_occultation._test()

This is a **double** occultation of `c`, as `b` goes into retrograde halfway through
the event! The duration is 300 minutes, or about 5 hours (!)

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/next_occultation.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Trappist1
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

def _test():
  '''
  
  '''
  
  plot()

def plot():
  '''
  
  '''
  
  # Instantiate the Trappist-1 system
  system = Trappist1(sample = True, nbody = True, seed = 1234)

  # Get the next 10 occultations of c
  # I'm starting the integration on May 26, 2017. The ephemerides aren't really accurate
  # given the transit times from Gillon et al. (2017), which are from October 2016.
  # But that's ok for this example.
  times, _, durations = system.next_occultation(system.c, occultors = system.b, noccultations = 10, tstart = 7900., tend = 8000.)
  
  # Grab the longest one
  t = times[np.argmax(durations)]
  
  # Now let's plot the light curve to check it out. Note that we need to re-instantiate
  # the system (with the same seed) since the integration already carried us past the occultation. 
  # We should also integrate it from the same `tstart` we used above to get the exact same orbital solution.
  
  # Get the light curve up to that point plus a little bit
  system = Trappist1(sample = True, nbody = True, seed = 1234)
  time = np.arange(7900., t + 0.1, MINUTE)
  system.compute(time)
  
  # Now plot just the occultation
  system.plot_occultation('c', t)
  pl.show()
  
if __name__ == '__main__':
  plot()