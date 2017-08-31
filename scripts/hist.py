#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
hist.py |github|
----------------

Histograms of the occultation events as a function of mean longitude, duration, and
impact parameter for each of the seven TRAPPIST-1 planets, as well as a marginalized
histogram of the total number of potentially detectable planet-planet occultations in
one Earth year.

.. note:: When I sample too many times from the prior in a single run, the code often \
          hangs. There is likely a memory leak somewhere, but I haven't been able to find \
          it yet. If you want to run large ensembles, I recommend parallelizing smaller \
          batches. A brain-dead way of doing it is to instantiate a bunch of **screen** \
          sessions: `screen -dm python -c "import hist; hist.Compute(nsamp = 100)"`
          
.. plot::
   :align: center
   
   from scripts import hist
   hist._test()

.. role:: raw-html(raw)
   :format: html
.. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/hist.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os
import planetplanet
from planetplanet import Trappist1
from planetplanet.constants import *
import matplotlib
import matplotlib.pyplot as pl
import numpy as np
import corner
from tqdm import tqdm
from scipy.stats import norm
datapath = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(planetplanet.__file__))), 'scripts', 'data')
if not os.path.exists(datapath):
  os.makedirs(datapath)

def _test():
  '''
  
  '''
  
  if not os.path.exists(os.path.join(datapath, 'hist000.npz')):
    Compute(nsamp = 1)
  figs = Plot()
  for fig in figs[:-1]:
    pl.close(fig)
  pl.show()

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
      h = system.histogram(OCTOBER_08_2016, OCTOBER_08_2016 + 365)
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
  
  while os.path.exists(os.path.join(datapath, 'hist%03d.npz' % n)): 
    n += 1
  np.savez(os.path.join(datapath, 'hist%03d.npz' % n), hist = hist, count = count)

def Plot():
  '''
  
  '''
  
  # Load
  print("Loading...")
  for n in tqdm(range(1000)):
    if os.path.exists(os.path.join(datapath, 'hist%03d.npz' % n)):
      data = np.load(os.path.join(datapath, 'hist%03d.npz' % n))
      data['hist'][0]
    else:
      if n == 0:
        raise Exception("Please run `Compute()` first.")
      break
    if n == 0:
      hist = data['hist']
      count = data['count']
    else:
      for k in range(7):
        hist[k] = np.vstack((hist[k], data['hist'][k]))
      count = np.hstack((count, data['count']))
  
  # For reference, total number of systems instantiated (samples)
  print("Total number of samples: %d" % len(count[0]))

  # For reference, the average number of occultations *per day* is
  occ_day = np.sum([hist[n].shape[0] for n in range(7)]) / count.shape[1] / 365
  print("Average number of occultations per day: %.2f" % occ_day)
  # I get 1.1 (!) These are occultations at all impact parameters and durations,
  # so most are grazing / not really detectable.
  
  # We will plot 8 different figures
  figs = [None for i in range(8)]
  
  # Corner plot
  for k, planet in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):  
    samples = np.array(hist[k])
    
    # But first check if we have enough samples
    if samples.shape[0] <= samples.shape[1]:
      figs[k] = pl.figure()
      continue
    
    figs[k] = corner.corner(samples, data_kwargs = {'alpha': 0.005}, range = [(-180,180), (0,1), (0, 3)], labels = ["Longitude [deg]", "Impact parameter", "Duration [min]"], bins = 30)
    for i, ax in enumerate(figs[k].axes):
      ax.set_xlabel(ax.get_xlabel(), fontsize = 14, fontweight = 'bold')
      ax.set_ylabel(ax.get_ylabel(), fontsize = 14, fontweight = 'bold')
      for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(12)
    for i in [0,3,6]:
      # IMPORTANT: The `histogram()` method returns the **orbital phase angle**, which is
      # measured from *transit* (phase = 0 deg). The mean longitude is measured from
      # *quadrature*, so there's a 90 deg offset we must apply.
      # Order is secondary eclipse, quadrature left, transit, quadrature right, secondary eclipse
      figs[k].axes[i].set_xticks([-180, -90, 0, 90, 180])
    figs[k].axes[6].set_xticklabels([r"$+$90", r"$\pm$180", r"$-$90", "0", r"$+$90"])
    figs[k].axes[3].set_yticks([0.2, 0.4, 0.6, 0.8])
    figs[k].axes[7].set_xticks([0.2, 0.4, 0.6, 0.8])
    figs[k].axes[8].set_xticks([0, 1, 2, 3])
    figs[k].axes[8].set_xticklabels([1, 10, 100, 1000])
    figs[k].axes[6].set_yticks([0, 1, 2, 3])
    figs[k].axes[6].set_yticklabels([1, 10, 100, 1000])
  
  # Frequency histogram
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  figs[-1] = pl.figure(figsize = (7, 8))
  figs[-1].subplots_adjust(hspace = 0.075)
  ax = pl.subplot2grid((5, 1), (1, 0), rowspan = 4)
  axt = pl.subplot2grid((5, 1), (0, 0), rowspan = 1, sharex = ax, zorder = -99)
  system = Trappist1(sample = True, nbody = False, quiet = True)
  for k, planet in enumerate(system.bodies[1:]): 
      
    # Fit a gaussian
    mu, sig = norm.fit(count[k])
    mu = '%.1f' % mu
    if len(mu) == 3: 
      label = r"$\mathbf{%s}: \ \ %s \pm %3.1f\ \mathrm{yr}^{-1}$" % (planet.name, mu, sig)
    else:
      label = r"$\mathbf{%s}: %s \pm %3.1f\ \mathrm{yr}^{-1}$" % (planet.name, mu, sig)
    
    # Plot
    ax.hist(count[k], color = planet.color, edgecolor = 'none', alpha = 0.25, histtype = 'stepfilled', normed = True, range = (0,50), zorder = 0, label = label, bins = 50)
    ax.hist(count[k], color = planet.color, histtype = 'step', normed = True, range = (0,50), zorder = 1, lw = 2, bins = 50)
    
    # HACK: Force a broken axis for planet `f`
    if planet.name == 'f':
      axt.hist(count[k], color = planet.color, edgecolor = 'none', alpha = 0.25, histtype = 'stepfilled', normed = True, range = (0,50), zorder = 0, label = label, bins = 50)
      axt.hist(count[k], color = planet.color, histtype = 'step', normed = True, range = (0,50), zorder = 1, lw = 2, bins = 50)
    
  leg = ax.legend(loc = 'upper right', fontsize = 15, bbox_to_anchor=(0.89, 0.865), bbox_transform = figs[-1].transFigure)
  ax.set_xlabel('Occultations per year', fontsize = 16, fontweight = 'bold')
  ax.set_ylabel('Probability', fontsize = 16, fontweight = 'bold')
  ax.yaxis.set_label_coords(-0.1, 0.6)
  
  for tick in ax.get_xticklabels() + ax.get_yticklabels() + axt.get_yticklabels():
    tick.set_fontsize(12)
  
  # HACK: Force a broken axis for planet `f`
  ax.set_ylim(0, 0.32)
  axt.set_ylim(0.4725, 0.55)
  axt.spines['bottom'].set_visible(False)
  ax.spines['top'].set_visible(False)
  axt.tick_params(bottom='off',labelbottom='off')
  axt.set_yticks([0.5, 0.55])
  d = .015
  kwargs = dict(transform=axt.transAxes, color='k', clip_on=False, lw = 1)
  axt.plot((-d, +d), (-d, +d), **kwargs)
  axt.plot((1 - d, 1 + d), (-d, +d), **kwargs)
  kwargs.update(transform=ax.transAxes)
  ax.plot((-d, +d), (1 - 0.25 * d, 1 + 0.25 * d), **kwargs)
  ax.plot((1 - d, 1 + d), (1 - 0.25 * d, 1 + 0.25 * d), **kwargs)
  
  return figs

if __name__ == '__main__':
  if not os.path.exists(os.path.join(datapath, 'hist000.npz')):
    Compute()
  figs = Plot()
  for k, planet in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):
    figs[k].savefig('%s.corner.pdf' % planet, bbox_inches = 'tight')
  figs[-1].savefig('hist.pdf', bbox_inches = 'tight')
  pl.close()