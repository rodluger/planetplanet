#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
hist_mutual.py |github|
-----------------------

Histograms of the mutual transit events in TRAPPIST-1. Shows histograms
of the fractional depth and duration of these events for all pairs of
planets.

TRAPPIST-1b
~~~~~~~~~~~
          
.. image:: /b_mutual.jpg
   :width: 400px
   :align: center

TRAPPIST-1c
~~~~~~~~~~~
          
.. image:: /c_mutual.jpg
   :width: 400px
   :align: center

TRAPPIST-1d
~~~~~~~~~~~
          
.. image:: /d_mutual.jpg
   :width: 400px
   :align: center

TRAPPIST-1e
~~~~~~~~~~~
          
.. image:: /e_mutual.jpg
   :width: 400px
   :align: center

TRAPPIST-1f
~~~~~~~~~~~
          
.. image:: /f_mutual.jpg
   :width: 400px
   :align: center

TRAPPIST-1g
~~~~~~~~~~~
          
.. image:: /g_mutual.jpg
   :width: 400px
   :align: center

TRAPPIST-1h
~~~~~~~~~~~
          
.. image:: /h_mutual.jpg
   :width: 400px
   :align: center

.. role:: raw-html(raw)
   :format: html
.. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/hist_mutual.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
import os
import subprocess
import planetplanet
from planetplanet import jwst
from planetplanet import Trappist1
from planetplanet.constants import *
from planetplanet.pool import Pool
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.ticker import FuncFormatter
import numpy as np
import corner
from tqdm import tqdm
from scipy.stats import norm
datapath = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.abspath(planetplanet.__file__))), 
                        'scripts', 'data')
histpath = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.abspath(planetplanet.__file__))), 
                        'scripts')
if not os.path.exists(datapath):
    os.makedirs(datapath)
                        
def _test():
    '''
    This routine is too expensive to test on Travis, so I'm
    bypassing it for now.
    
    '''
    
    pass

def Submit(queue = None, email = None, walltime = 8, nodes = 5, ppn = 12,
           mpn = None, nsamp = 50000, batch_size = 30, nproc = None):
    '''
    Submits a PBS cluster job to run :py:func:`Compute` in parallel.

    :param str queue: The name of the queue to submit to. \
           Default :py:obj:`None`
    :param str email: The email to send job status notifications to. \
           Default :py:obj:`None`
    :param int walltime: The number of hours to request. Default `8`
    :param int nodes: The number of nodes to request. Default `5`
    :param int ppn: The number of processors per node to request. Default `12`
    :param int nsamp: The number of prior samples to draw. Default `50,000`
    :param int batch_size: Size of each batch used in the parallelization. \
           Default `100`
    :param int mpn: Memory per node in gb to request. Default no setting.
    :param int nproc: Number of processes to spawn. Default is the number of \
           core.
    '''
        
    if nproc is None:
        nproc = ppn * nodes
    str_w = 'walltime=%d:00:00' % walltime
    if mpn is not None:
        str_n = 'nodes=%d:ppn=%d,feature=%dcore,mem=%dgb' % \
                (nodes, ppn, ppn, mpn * nodes)
    else:
        str_n = 'nodes=%d:ppn=%d,feature=%dcore' % (nodes, ppn, ppn)    
    str_v = 'NPROC=%d,HISTPATH=%s,NSAMP=%d,BATCHSZ=%d' % \
            (nproc, histpath, nsamp, batch_size)
    str_name = 'planetplanet'
    str_out = 'hist_mutual.log'
    qsub_args = ['qsub', 'hist_mutual.pbs', 
                 '-v', str_v, 
                 '-o', str_out,
                 '-j', 'oe', 
                 '-N', str_name,
                 '-l', str_n,
                 '-l', str_w]
    if email is not None: 
        qsub_args.append(['-M', email, '-m', 'ae'])
    if queue is not None:
        qsub_args += ['-q', queue] 
    print("Submitting the job...")
    subprocess.call(qsub_args)

class _FunctionWrapper(object):
    '''
    A simple function wrapper class. Stores :py:obj:`args` and :py:obj:`kwargs` 
    and allows an arbitrary function to be called via :py:func:`map`.
    Used internally.
    
    '''
    
    def __init__(self, f, *args, **kwargs):
        '''
        
        '''
        
        self.f = f
        self.args = args
        self.kwargs = kwargs
    
    def __call__(self, x):
        '''
        
        '''
        
        return self.f(*self.args, **self.kwargs)

def _Parallelize(nsamp, batch_size):
    '''
    Runs the actual parallelized computations. Used internally.
    
    '''

    # Get our function wrapper
    m = _FunctionWrapper(Compute, nsamp = batch_size,
                         progress_bar = False)
    
    # Parallelize. We will run `N` iterations
    N = int(np.ceil(nsamp / batch_size))
    with Pool() as pool:
        pool.map(m, range(N))

def histogram(system, tstart, tend, dt = 0.0001):
    '''
    Computes statistical properties of mutual events (PPOs occuring on the face
    of the star).

    :param system: A system instance.
    :type system: :py:obj:`planetplanet.structs.System`
    :param float tstart: The integration start time (BJD − 2,450,000)
    :param float tend: The integration end time (BJD − 2,450,000)
    :param float dt: The time resolution in days. Occultations shorter \
           than this will not be registered.

    '''

    # Compute the orbits
    time = np.arange(tstart, tend, dt)
    system.compute_orbits(time)
    pairs = []
    durs = []
    depths = []
    for bi, body in enumerate(system.bodies[1:]):
        
        # Loop over all times w/ occultations of the star
        inds = np.where(system.bodies[0].occultor)[0]
        for i in inds:

            # Get all bodies currently occulting the star
            occultors = []
            for occ in range(1, len(system.bodies)):
                if (system.bodies[0].occultor[i] & 2 ** occ):
                    occultors.append(occ)

            # Check if any of these occult each other
            if len(occultors) > 1:
                for occ1 in occultors:
                    for occ2 in occultors:
                        if system.bodies[occ1].occultor[i] & 2 ** occ2:
                            
                            # Sort the planet indices
                            occ1, occ2 = sorted([occ1, occ2])
                            
                            # Is this a new occultation?
                            if (len(pairs)==0) or (pairs[-1] != [occ1, occ2]):
                                pairs.append([occ1, occ2])
                                durs.append(0.)
                                depths.append(0.)
                            
                            # Update the maximum depth. Based on
                            # http://mathworld.wolfram.com/
                            # Circle-CircleIntersection.html
                            dx2 = (system.bodies[occ1].x[i] - 
                                   system.bodies[occ2].x[i]) ** 2
                            dy2 = (system.bodies[occ1].y[i] - 
                                   system.bodies[occ2].y[i]) ** 2
                            d = np.sqrt(dx2 + dy2)
                            r1 = system.bodies[occ1]._r
                            r2 = system.bodies[occ2]._r
                            A = r1 ** 2 * np.arccos((d ** 2 + r1 ** 2 - r2 ** 2) 
                                                   / (2 * d * r1)) \
                              + r2 ** 2 * np.arccos((d ** 2 + r2 ** 2 - r1 ** 2) 
                                                   / (2 * d * r2)) \
                              - 0.5 * np.sqrt((-d + r1 + r2) *
                                               (d + r1 - r2) *
                                               (d - r1 + r2) *
                                               (d + r1 + r2))
                            Astar = np.pi * (system.bodies[0]._r) ** 2
                            depth = A / Astar
                            depths[-1] = max(depths[-1], depth)
                            
                            # Update the duration
                            durs[-1] += dt
    
    return np.array(pairs, dtype = int), \
           np.array(durs, dtype = float), \
           np.array(depths, dtype = float)
    
def Compute(nsamp = 100, nbody = True, progress_bar = True, **kwargs):
    '''
    Runs the simulations.
    
    :param int nsamp: The number of prior samples to draw. Default `300`
    :param bool nbody: Use the N-Body solver? Default :py:obj:`True`
    :param bool progress_bar: Display a progress bar? Default :py:obj:`True`

    '''
       
    # Draw samples from the prior
    pairs = np.empty([0, 2], dtype = int)
    durs = np.array([], dtype = float)
    depths = np.array([], dtype = float)
    if progress_bar:
        wrap = tqdm
    else:
        wrap = lambda x: x
    for n in wrap(range(nsamp)):
        
        # Instantiate the Trappist-1 system
        system = Trappist1(sample = True, nbody = nbody, 
                           quiet = True, **kwargs)
        system.settings.timestep = 1. / 24.
        
        # Run!
        try:
            p, t, d = histogram(system, OCTOBER_08_2016, OCTOBER_08_2016 + 365)
        except:
            print("ERROR in routine `hist.Compute()`")
            continue
        
        if len(p):
            pairs = np.vstack([pairs, p])
            durs = np.append(durs, t)
            depths = np.append(depths, d)

    # Save
    n = 0
    while os.path.exists(os.path.join(datapath, 'hist_mutual%03d.npz' % n)): 
        n += 1
    np.savez(os.path.join(datapath, 'hist_mutual%03d.npz' % n), 
             pairs = pairs, durs = durs, depths = depths)

def MergeFiles():
    '''
    Merge all the `npz` savesets into a single one for faster plotting.
    
    '''
    
    # Load
    pairs = np.empty([0, 2], dtype = int)
    durs = np.array([], dtype = float)
    depths = np.array([], dtype = float)
    print("Loading...")
    for n in tqdm(range(1000)):
        if os.path.exists(os.path.join(datapath, 'hist_mutual%03d.npz' % n)):
            
            # Skip corrupt files
            try:
                data = np.load(os.path.join(datapath, 'hist_mutual%03d.npz' % n))
                os.remove(os.path.join(datapath, 'hist_mutual%03d.npz' % n))
                data['pairs'][0]
                data['durs'][0]
                data['depths'][0]
            except:
                continue
        
        else:
            break
        
        pairs = np.vstack([pairs, data['pairs']])
        durs = np.append(durs, data['durs'])
        depths = np.append(depths, data['depths'])
    
    # Save as one big file
    if n > 0:
        print("Saving...")
        np.savez(os.path.join(datapath,'hist_mutual000.npz'), 
                 pairs = pairs, durs = durs, depths = depths)

def Plot():
    '''
    
    '''
    
    # Load
    pairs = np.empty([0, 2], dtype = int)
    durs = np.array([], dtype = float)
    depths = np.array([], dtype = float)
    print("Loading...")
    for n in tqdm(range(1000)):
        if os.path.exists(os.path.join(datapath, 'hist_mutual%03d.npz' % n)):
            
            # Skip corrupt files
            try:
                data = np.load(os.path.join(datapath, 'hist_mutual%03d.npz' % n))
                data['pairs'][0]
                data['durs'][0]
                data['depths'][0]
            except:
                continue
        
        else:
            if n == 0:
                raise Exception("Please run `Compute()` first.")
            break
        
        pairs = np.vstack([pairs, data['pairs']])
        durs = np.append(durs, data['durs'])
        depths = np.append(depths, data['depths'])
    
    # Dummy system to get colors
    system = Trappist1()
    colors = [system.bodies[n].color for n in range(1, 8)]
    
    # For the paper, we ran 30,000 simulations, so the average
    # number of mutual transits per year is...
    print("Mutual transits per year: %.3f" % (len(pairs) / 30000.))
    
    # Loop over all planets
    for k in range(1, 8):
        
        # Indices of events involving this planet
        inds = np.where((pairs[:,0] == k) | (pairs[:,1] == k))[0]
        
        # Again, for the 30,000 simulations we ran...
        print("%s: %.3f" % (system.bodies[k].name, len(pairs[inds]) / 30000.))
        
        # Duration
        dt = durs[inds] / MINUTE
        
        # Depth
        d = depths[inds] * 1e2
        
        # Corner plot
        samples = np.vstack((dt, d)).T
        fig = corner.corner(samples, plot_datapoints = False,
                            range = [(0, 60), (0, 1)], 
                            labels = ["Duration [min]", 
                                      "Depth [%]"], 
                            bins = 30,
                            hist_kwargs = {'color': 'w'})
        
        # Indices of events involving each of the planets
        pinds = [[] for j in range(1, 8)]
        for j in range(1, 8):
            if j != k:
                pinds[j - 1] = np.where((pairs[inds,0] == j) | (pairs[inds,1] == j))[0]

        # Duration stacked histogram
        n, _, _ = fig.axes[0].hist([dt[p] for p in pinds], bins = 30,
                                   range = (0, 60), 
                                   stacked = True,
                                   normed = 1,
                                   color = colors,
                                   alpha = 0.5)
        maxn = np.max(np.array(n)[-1])
        fig.axes[0].hist(dt, bins = 30, range = (0, 60), normed = 1,
                         color = 'k', alpha = 1, histtype = 'step')
        fig.axes[0].set_ylim(0, 1.1 * maxn)
        
        # Depth stacked histogram
        n, _, _ = fig.axes[3].hist([d[p] for p in pinds], bins = 30,
                                   range = (0, 1), 
                                   stacked = True,
                                   normed = 1,
                                   color = colors,
                                   alpha = 0.5)
        maxn = np.max(np.array(n)[-1])
        fig.axes[3].hist(d, bins = 30, range = (0, 1), normed = 1,
                         color = 'k', alpha = 1, histtype = 'step')
        fig.axes[3].set_ylim(0, 1.1 * maxn)
        
        # Tweak appearance
        for i, ax in enumerate(fig.axes):
            ax.set_xlabel(ax.get_xlabel(), fontsize = 14, fontweight = 'bold')
            ax.set_ylabel(ax.get_ylabel(), fontsize = 14, fontweight = 'bold')
            for tick in ax.get_xticklabels() + ax.get_yticklabels():
                tick.set_fontsize(12)
        
        # Save!
        fig.savefig('%s_mutual.pdf' % system.bodies[k].name, 
                    bbox_inches = 'tight')