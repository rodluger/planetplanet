#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
hist.py |github|
----------------

Histograms of the occultation events as a function of mean longitude, duration, 
and impact parameter for each of the seven TRAPPIST-1 planets, as well as a 
marginalized histogram of the total number of potentially detectable 
planet-planet occultations in one Earth year.

.. note:: When I sample too many times from the prior in a single run, the \
          code often hangs. There may be a memory leak somewhere, but I \
          haven't been able to find it yet. If you want to run large \
          ensembles, I recommend using the parallelization scheme I \
          implemented below. Alternatively, a brain-dead \
          way of doing it is to instantiate a bunch of **screen** \
          sessions: \
          :py:obj:`screen -dm python -c "import hist; hist.Compute(nsamp = 100)"`

TRAPPIST-1b
~~~~~~~~~~~
          
.. image:: /b.corner.jpg
   :width: 400px
   :align: center

TRAPPIST-1c
~~~~~~~~~~~
          
.. image:: /c.corner.jpg
   :width: 400px
   :align: center

TRAPPIST-1d
~~~~~~~~~~~
          
.. image:: /d.corner.jpg
   :width: 400px
   :align: center

TRAPPIST-1e
~~~~~~~~~~~
          
.. image:: /e.corner.jpg
   :width: 400px
   :align: center

TRAPPIST-1f
~~~~~~~~~~~
          
.. image:: /f.corner.jpg
   :width: 400px
   :align: center

TRAPPIST-1g
~~~~~~~~~~~
          
.. image:: /g.corner.jpg
   :width: 400px
   :align: center

TRAPPIST-1h
~~~~~~~~~~~
          
.. image:: /h.corner.jpg
   :width: 400px
   :align: center
   
Marginal distributions
~~~~~~~~~~~~~~~~~~~~~~
          
.. image:: /hist.jpg
   :width: 400px
   :align: center

.. role:: raw-html(raw)
   :format: html
.. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/hist.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

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
    
    '''
    
    if not os.path.exists(os.path.join(datapath, 'hist_d000.npz')):
        Compute(nsamp = 1, photo = False)
    fig_corner, _, fig_hist = Plot(photo = False)
    for fig in fig_corner:
        pl.close(fig)
    pl.show()

def Submit(queue = None, email = None, walltime = 8, nodes = 5, ppn = 12,
           mpn = None, nsamp = 50000, photo = True, eyeball = True, 
           batch_size = 100, nproc = None):
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
    :param bool photo: Run the full photodynamical model? \
           Default :py:obj:`True`
    :param bool eyeball: Use the radiative equilibrium surface map? \
           Default :py:obj:`True`
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
    str_v = 'NPROC=%d,HISTPATH=%s,NSAMP=%d,PHOTO=%d,EYEBALL=%d,BATCHSZ=%d' % \
            (nproc, histpath, nsamp, int(photo), int(eyeball), batch_size)
    str_name = 'planetplanet'
    str_out = 'hist.log'
    qsub_args = ['qsub', 'hist.pbs', 
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

def _Parallelize(nsamp, photo, eyeball, batch_size):
    '''
    Runs the actual parallelized computations. Used internally.
    
    '''

    # Get our function wrapper
    m = _FunctionWrapper(Compute, nsamp = batch_size, photo = bool(photo), 
                         eyeball = eyeball, progress_bar = False)
    
    # Parallelize. We will run `N` iterations
    N = int(np.ceil(nsamp / batch_size))
    with Pool() as pool:
        pool.map(m, range(N))

def histogram(system, tstart, tend, dt = 0.0001, photo = False):
    '''
    Computes statistical properties of planet-planet occultations that 
    occur over a given interval. Computes the frequency of occultations as 
    a function of orbital phase, duration, and impact parameter, as well 
    as the fully marginalized occultation frequency for each planet in the 
    system. Occultations by the star are not included, nor are occultations
    occuring behind the star, which are not visible to the observer.

    :param system: A system instance.
    :type system: :py:obj:`planetplanet.structs.System`
    :param float tstart: The integration start time (BJD − 2,450,000)
    :param float tend: The integration end time (BJD − 2,450,000)
    :param float dt: The time resolution in days. Occultations shorter \
           than this will not be registered.
    :param bool photo: Run the full photodynamical routine? Default \
           :py:obj:`False`, in which case only dynamical properties \
           (durations, impact parameters, and phases) are computed. \
           If :py:obj:`True`, computes photometric statistics, including \
           occultation depths and SNRs. Note that calling the full \
           photodynamical code is much more computationally expensive.
           
    :returns: A list of \
              :py:obj:`(phase angle, impact parameter, duration)` tuples \
              for each planet in the system. The phase angle is measured \
              in degrees and the duration is measured in days. If `photo` \
              is :py:obj:`True`, the tuples contain three additional floats: \
              :py:obj:`signal`, :py:obj:`noise`, :py:obj:`snr`; the first two \
              are measured in ppm.

    .. warning:: This routine computes the **orbital phase angle**, which \
                 is measured from **transit**. This is different from the \
                 mean longitude by :math:`\pi/2`

    '''

    # Reset
    system._reset()
    time = np.arange(tstart, tend, dt)
    
    # Dynamics only
    if not photo:
    
        for body in system.bodies:
            body.u = np.array([], dtype = float)
            
        system._malloc(len(time), 1)

        # Call the orbit routine
        err = system._Orbits(len(time), np.ctypeslib.as_ctypes(time), 
                           len(system.bodies), system._ptr_bodies, 
                           system.settings)
        assert err <= 0, "Error in C routine `Orbits` (%d)." % err
    
    # Full photodynamical model
    else:
    
        # Compute the wavelength grid. We are hard-coding the
        # 15 micron JWST MIRI band here.
        lambda1 = 12.5
        lambda2 = 17.5
        R = 50
        wav = [lambda1]
        while(wav[-1] < lambda2):
            wav.append(wav[-1] + wav[-1] / R)
        wavelength = np.array(wav)
        
        # Compute all limb darkening coefficients
        for body in system.bodies:
            body.u = [None for ld in body.limbdark]
            for n, ld in enumerate(body.limbdark):
                if callable(ld):
                    body.u[n] = ld(wavelength)
                elif not hasattr(ld, '__len__'):
                    body.u[n] = ld * np.ones_like(wavelength)
                else:
                    raise Exception("Limb darkening coefficients must be "
                                    + "provided as a list of scalars or "
                                    + "as a list of functions.")
            body.u = np.array(body.u)
            
            # HACK: Disable phase curves. The code would take *way*
            # too long to run, and they don't affect these statistics.
            body.phasecurve = False

        # Convert from microns to meters
        wavelength *= 1e-6
                    
        # No oversampling
        time_hr = np.array(time)
        
        # Continuum flux
        system._continuum = np.zeros(len(time_hr) * len(wavelength))
        
        # Allocate memory
        system._malloc(len(time_hr), len(wavelength))

        # Call the light curve routine
        err = system._Flux(len(time_hr), np.ctypeslib.as_ctypes(time_hr), 
                         len(wavelength), 
                         np.ctypeslib.as_ctypes(wavelength), 
                         np.ctypeslib.as_ctypes(system._continuum), 
                         len(system.bodies), system._ptr_bodies, 
                         system.settings)
        assert err <= 0, "Error in C routine `Flux` (%d)." % err
        
    # A histogram of the distribution of phases, 
    # impact parameters, and durations
    hist = [[] for body in system.bodies[1:]]
    for k, body in enumerate(system.bodies[1:]):
        
        # Full photo model?
        if photo:
        
            # Simulate an observation w/ JWST at 15 microns
            # Same syntax as in `observe()`
            w = jwst.get_miri_filter_wheel()
            filter = w[np.argmax([f.name.lower() == 'f1500w' for f in w])]
            filter.compute_lightcurve(time, body.flux, 
                                      system.continuum, 
                                      system.wavelength, 
                                      stack = 1,
                                      atel = 25., 
                                      thermal = True,
                                      quiet = True)
        
        # Identify the different planet-planet events
        inds = np.where(body.occultor > 0)[0]
        difs = np.where(np.diff(inds) > 1)[0]

        # Loop over individual ones
        for i in inds[difs]:
            
            # Loop over possible occultors
            for occ in range(1, len(system.bodies)):

                # Is body `occ` occulting (but not behind the star)?
                if (body.occultor[i] & 2 ** occ) and \
                   (body.occultor[i] & 1 == 0):

                    # Note that `i` is the last index of the occultation
                    duration = np.argmax(body.occultor[:i][::-1] 
                                         & 2 ** occ == 0)
                    if duration > 0:

                        # Orbital phase, **measured from transit**
                        # At transit, phase = 0; at secondary, phase = 180.
                        phase = np.arctan2(body.x[i], 
                                          -body.z[i]) * 180 / np.pi

                        # Indices of the occultation
                        idx = range(i - duration, i + 1)
                        
                        # Compute the minimum impact parameter
                        impact = np.min(np.sqrt((system.bodies[occ].x[idx] 
                                                 - body.x[idx]) ** 2 +
                                                (system.bodies[occ].y[idx] 
                                                - body.y[idx]) ** 2)) \
                                                / (system.bodies[occ]._r 
                                                   + body._r)

                        # Convert duration to log
                        duration = np.log10(duration * dt * 1440)
                        
                        # Are we done?
                        if not photo:
                        
                            # Running list
                            hist[k].append((phase, impact, duration))
                        
                        # No. Let's get some photometric info
                        else:

                            # Planet, background, and star photons
                            Nplan = filter.lightcurve.Nsys[idx]
                            Nback = filter.lightcurve.Nback[idx]
                            Nstar = filter.lightcurve.Ncont[idx]
                            
                            # Compute the number of photons *missing*
                            Nplan = np.nanmedian(Nplan) - Nplan
                            
                            # Compute signal of and noise on the event
                            # in parts per million
                            norm = 1.e6 / np.sum(Nstar + Nback)
                            signal = norm * np.sum(np.fabs(Nplan)) 
                            noise = norm * np.sqrt(np.sum(Nstar + Nback))
                            
                            # Compute the actual SNR on event. Note that this 
                            # is NOT the sum of the signals divided by the sum 
                            # of the noises! We need to add the SNR of each 
                            # *datapoint* individually in quadrature.
                            snr = np.sqrt(np.sum((Nplan) ** 2 
                                          / (Nstar + Nback)))
                            
                            # Running list
                            hist[k].append((phase, impact, duration,
                                            signal, noise, snr))
        # Make into array
        hist[k] = np.array(hist[k])
        
    return hist
    
def Compute(nsamp = 300, mind = 10., maxb = 0.5, nbody = True, photo = True,
            progress_bar = True, eyeball = True, **kwargs):
    '''
    Compute occultation histograms by drawing `nsamp` draws from the 
    system prior. Saves the results to `data/histXXX.npz`.
    
    :param int nsamp: The number of prior samples to draw. Default `300`
    :param float mind: The minimum occultation duration in minutes to include\
           in the tally. Default `10.`
    :param float maxb: The maximum occultation impact parameter to include\
           in the tally. Default `0.5`
    :param bool nbody: Use the N-Body solver? Default :py:obj:`True`
    :param bool nbody: Run the full photodynamical code? Default :py:obj:`True`
    :param bool progress_bar: Display a progress bar? Default :py:obj:`True`
    :param bool eyeball: Assume eyeball planets? Default :py:obj:`True`. If
           :py:obj:`False`, uses the limb darkened surface map.
    
    '''
       
    # The dataset name
    name = 'hist_' + ('p' + ('e' if eyeball else 'l') if photo else 'd')   
        
    # Draw samples from the prior
    hist = [[] for k in range(7)]
    count = [np.zeros(nsamp) for k in range(7)]
    if progress_bar:
        wrap = tqdm
    else:
        wrap = lambda x: x
    for n in wrap(range(nsamp)):

        # Instantiate the Trappist-1 system
        system = Trappist1(sample = True, nbody = nbody, 
                           quiet = True, **kwargs)
        system.settings.timestep = 1. / 24.
        
        # Choose the correct radiance map
        if eyeball:
            radiancemap = planetplanet.RadiativeEquilibriumMap()
        else:
            radiancemap = planetplanet.LimbDarkenedMap()
        
        # Run!
        try:
            h = histogram(system, OCTOBER_08_2016, OCTOBER_08_2016 + 365, 
                          photo = photo)
        except:
            print("ERROR in routine `hist.Compute()`")
            continue
            
        # Loop over the planets
        for k in range(7):
    
            # Count the number of events longer than `mind`
            # and with impact parameter below `minb`
            if len(h[k]):
                count[k][n] = len(np.where((h[k][:,1] < maxb) & 
                                 (h[k][:,2] > np.log10(mind)))[0])
        
            # Append to cumulative histogram
            hist[k].extend(list(h[k]))
    
    # Convert to numpy arrays
    for k in range(7):
        hist[k] = np.array(hist[k])
    
    # Save
    n = 0
    
    while os.path.exists(os.path.join(datapath, '%s%03d.npz' % (name, n))): 
        n += 1
    np.savez(os.path.join(datapath, '%s%03d.npz' % (name, n)), 
             hist = hist, count = count)

def Plot(photo = True, eyeball = True):
    '''
    Plots the results of a `Compute()` run and returns several figures.
    
    '''
    
    # The dataset name
    name = 'hist_' + ('p' + ('e' if eyeball else 'l') if photo else 'd')  
    
    # Load
    print("Loading...")
    for n in tqdm(range(1000)):
        if os.path.exists(os.path.join(datapath, '%s%03d.npz' % (name, n))):
            data = np.load(os.path.join(datapath, '%s%03d.npz' % (name, n)))
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
    occ_day = np.sum([hist[n].shape[0] for n in range(7)]) \
              / count.shape[1] / 365
    print("Average number of occultations per day: %.2f" % occ_day)
    # I get 1.1 (!) These are occultations at all impact parameters 
    # and durations, so most are grazing / not really detectable.
        
    # All the figures we'll plot
    fig_corner = [None for i in range(7)]
    fig_snr = [None for i in range(7)]
    fig_hist = None

    # Corner plot
    for k, planet in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):    
        samples = np.array(hist[k])
        
        # But first check if we have enough samples
        if not len(samples):
            fig_snr[k] = pl.figure()
            continue
        
        # Did we run the full photodynamical model?
        if samples.shape[1] == 6:
            snr = samples[:,5]
            fig_snr[k], ax = pl.subplots(1)
            color = 'cornflowerblue'
            ax.hist(snr, 
                    weights = np.ones_like(snr) / len(count[0]),
                    color = color, edgecolor = 'none', 
                    alpha = 0.5, histtype = 'stepfilled', 
                    bins = 50)
            ax.hist(snr, 
                    weights = np.ones_like(snr) / len(count[0]),
                    color = color, histtype = 'step', 
                    lw = 2, bins = 50)
            ax.set_xlabel('SNR', fontsize = 14, fontweight = 'bold')
            ax.set_ylabel(r'Occultations [yr$^{-1}$]', fontsize = 14, 
                          fontweight = 'bold')
            samples = samples[:,:3]
                 
        # But first check if we have enough samples
        if samples.shape[0] <= samples.shape[1]:
            fig_corner[k] = pl.figure()
            continue
        
        fig_corner[k] = corner.corner(samples, data_kwargs = {'alpha': 0.005}, 
                                range = [(-180,180), (0,1), (0, 3)], 
                                labels = ["Longitude [deg]", 
                                          "Impact parameter", 
                                          "Duration [min]"], bins = 30)
        for i, ax in enumerate(fig_corner[k].axes):
            ax.set_xlabel(ax.get_xlabel(), fontsize = 14, fontweight = 'bold')
            ax.set_ylabel(ax.get_ylabel(), fontsize = 14, fontweight = 'bold')
            for tick in ax.get_xticklabels() + ax.get_yticklabels():
                tick.set_fontsize(12)
        for i in [0,3,6]:
            # IMPORTANT: The `histogram()` method returns the 
            # **orbital phase angle**, which is
            # measured from *transit* (phase = 0 deg). The mean longitude 
            # is measured from *quadrature*, so there's a 90 deg offset we 
            # must apply. Order is secondary eclipse, quadrature left, 
            # transit, quadrature right, secondary eclipse
            fig_corner[k].axes[i].set_xticks([-180, -90, 0, 90, 180])
        fig_corner[k].axes[6].set_xticklabels([r"$+$90", r"$\pm$180", 
                                               r"$-$90", "0", r"$+$90"])
        fig_corner[k].axes[3].set_yticks([0.2, 0.4, 0.6, 0.8])
        fig_corner[k].axes[7].set_xticks([0.2, 0.4, 0.6, 0.8])
        fig_corner[k].axes[8].set_xticks([0, 1, 2, 3])
        fig_corner[k].axes[8].set_xticklabels([1, 10, 100, 1000])
        fig_corner[k].axes[6].set_yticks([0, 1, 2, 3])
        fig_corner[k].axes[6].set_yticklabels([1, 10, 100, 1000])
    
    # Frequency histogram
    matplotlib.rcParams['mathtext.fontset'] = 'cm'
    fig_hist = pl.figure(figsize = (7, 8))
    fig_hist.subplots_adjust(hspace = 0.075)
    ax = pl.subplot2grid((5, 1), (1, 0), rowspan = 4)
    axt = pl.subplot2grid((5, 1), (0, 0), rowspan = 1, 
                          sharex = ax, zorder = -99)
    system = Trappist1(sample = True, nbody = False, quiet = True)
    for k, planet in enumerate(system.bodies[1:]): 
            
        # Fit a gaussian
        mu, sig = norm.fit(count[k])
        mu = '%.1f' % mu
        if len(mu) == 3: 
            label = r"$\mathbf{%s}: \ \ %s \pm %3.1f\ \mathrm{yr}^{-1}$" \
                    % (planet.name, mu, sig)
        else:
            label = r"$\mathbf{%s}: %s \pm %3.1f\ \mathrm{yr}^{-1}$" \
                    % (planet.name, mu, sig)
        
        # Plot
        ax.hist(count[k], color = planet.color, edgecolor = 'none', 
                 alpha = 0.25, histtype = 'stepfilled', normed = True, 
                 range = (0,50), zorder = 0, label = label, bins = 50)
        ax.hist(count[k], color = planet.color, histtype = 'step', 
                normed = True, range = (0,50), zorder = 1, lw = 2, bins = 50)
        
        # HACK: Force a broken axis for planet `f`
        if planet.name == 'f':
            axt.hist(count[k], color = planet.color, edgecolor = 'none', 
                     alpha = 0.25, histtype = 'stepfilled', normed = True, 
                     range = (0,50), zorder = 0, label = label, bins = 50)
            axt.hist(count[k], color = planet.color, histtype = 'step', 
                     normed = True, range = (0,50), zorder = 1, 
                     lw = 2, bins = 50)
        
    leg = ax.legend(loc = 'upper right', fontsize = 15, 
                    bbox_to_anchor = (0.89, 0.865), 
                    bbox_transform = fig_hist.transFigure)
    ax.set_xlabel('Occultations per year', fontsize = 16, fontweight = 'bold')
    ax.set_ylabel('Probability', fontsize = 16, fontweight = 'bold')
    ax.yaxis.set_label_coords(-0.1, 0.6)
    
    for tick in ax.get_xticklabels() + ax.get_yticklabels() \
              + axt.get_yticklabels():
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
    
    return fig_corner, fig_snr, fig_hist

def MakeFigures(photo = True, eyeball = True):
    '''
    Plots all histogram figures for the paper.
    
    '''
    
    # Plot all figures   
    fig_corner, fig_snr, fig_hist = Plot(photo = photo, eyeball = eyeball)
    for k, planet in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):
        fig_corner[k].savefig('%s.corner.pdf' % planet, bbox_inches = 'tight')
        if fig_snr is not None:
            if eyeball:
                fig_snr[k].savefig('%s.snr.eyeball.pdf' % planet, 
                                   bbox_inches = 'tight')
            else:
                fig_snr[k].savefig('%s.snr.limbdark.pdf' % planet, 
                                   bbox_inches = 'tight')
    fig_hist.savefig('hist.pdf', bbox_inches = 'tight')
    pl.close()