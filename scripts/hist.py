#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
hist.py
-------

Histograms of the occultation events as a function of phase, duration, and
impact parameter for each of the seven TRAPPIST-1 planets.

```
screen -dm python -c "import hist; hist.Compute(nsamp = 100)"
```

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib
import matplotlib.pyplot as pl
import numpy as np
import corner
from tqdm import tqdm
from zipfile import BadZipFile
from scipy.stats import norm

def Compute(nsamp = 3000, mind = 10., maxb = 0.5, nbody = True):
  '''
  
  '''
  
  # Draw samples from the prior
  hist = [[] for k in range(7)]
  count = [np.zeros(nsamp) for k in range(7)]
  for n in tqdm(range(nsamp)):

    # Instantiate the Trappist-1 system
    system = Trappist1(sample = True, nbody = nbody, quiet = True)
    system.settings.timestep = 1. / 24.
    try:
      h = system.histogram(0, 365)
    except:
      continue
      
    # Loop over the planets
    for k in range(7):
  
      # Count the number of events longer than `mind`
      # and with impact parameter below `minb`
      if len(h[k]):
        count[k][n] = len(np.where((h[k][:,1] < maxb) & (h[k][:,2] > np.log10(mind)))[0])
    
      # Append to cumulative histogram
      hist[k].extend(list(h[k]))
  
  # Convert to numpy arrays
  for k in range(7):
    hist[k] = np.array(hist[k])
  
  # Save
  n = 0
  while os.path.exists('hist%03d.npz' % n): 
    n += 1
  np.savez('hist%03d.npz' % n, hist = hist, count = count)

def Plot():
  '''
  
  '''
  
  # Load
  print("Loading...")
  for n in tqdm(range(1000)):
    try:
      data = np.load('hist%03d.npz' % n)
      data['hist'][0]
    except FileNotFoundError:
      if n == 0:
        raise Exception("Please run `Compute()` first.")
      break
    except BadZipFile:
      print("Bad zip file: %d." % n)
      continue
    if n == 0:
      hist = data['hist']
      count = data['count']
    else:
      for k in range(7):
        hist[k] = np.vstack((hist[k], data['hist'][k]))
      count = np.hstack((count, data['count']))
  
  # Corner plot
  for k, planet in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):  
    samples = np.array(hist[k])
    fig = corner.corner(samples, data_kwargs = {'alpha': 0.005}, range = [(-180,180), (0,1), (0, 3)], labels = ["Phase [deg]", "Impact parameter", "Duration [min]"])
    for i, ax in enumerate(fig.axes):
      ax.set_xlabel(ax.get_xlabel(), fontsize = 14, fontweight = 'bold')
      ax.set_ylabel(ax.get_ylabel(), fontsize = 14, fontweight = 'bold')
      for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(12)
    fig.axes[0].set_xticks([-180, -90, 0, 90, 180])
    fig.axes[3].set_xticks([-180, -90, 0, 90, 180])
    fig.axes[6].set_xticks([-180, -90, 0, 90, 180])
    fig.axes[3].set_yticks([0.2, 0.4, 0.6, 0.8])
    fig.axes[7].set_xticks([0.2, 0.4, 0.6, 0.8])
    fig.axes[8].set_xticks([0, 1, 2, 3])
    fig.axes[8].set_xticklabels([1, 10, 100, 1000])
    fig.axes[6].set_yticks([0, 1, 2, 3])
    fig.axes[6].set_yticklabels([1, 10, 100, 1000])
    fig.savefig('../img/%s.corner.pdf' % planet, bbox_inches = 'tight')
    pl.close()
  
  # Frequency histogram
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  fig, ax = pl.subplots(1, figsize = (8, 6))
  system = Trappist1(sample = True, nbody = False, quiet = True)
  for k, planet in enumerate(system.bodies[1:]): 
      
    # Fit a gaussians
    mu, sig = norm.fit(count[k])
    mu = '%.1f' % mu
    if len(mu) == 3: 
      label = r"$\mathbf{%s}: \ \ %s \pm %3.1f\ \mathrm{yr}^{-1}$" % (planet.name, mu, sig)
    else:
      label = r"$\mathbf{%s}: %s \pm %3.1f\ \mathrm{yr}^{-1}$" % (planet.name, mu, sig)
    
    # Plot
    ax.hist(count[k], color = planet.color, edgecolor = 'none', alpha = 0.25, histtype = 'stepfilled', normed = True, range = (0,50), zorder = 0, label = label, bins = 50)
    ax.hist(count[k], color = planet.color, histtype = 'step', normed = True, range = (0,50), zorder = 1, lw = 2, bins = 50)
  
  leg = ax.legend(loc = 'upper right', fontsize = 15)
  ax.set_xlabel('Occultations per year', fontsize = 16, fontweight = 'bold')
  ax.set_ylabel('Probability', fontsize = 16, fontweight = 'bold')
  for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(12)
  fig.savefig('../img/hist.pdf', bbox_inches = 'tight')
  pl.close()
  
if __name__ == '__main__':
  Compute()
  Plot()