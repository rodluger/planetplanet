dynamics |github|
=================

The analytic dynamics module of :py:obj:`planetplanet`. These are a collection of
:py:obj:`Julia` scripts and are separate from the rest of the code.

offset_diagram2.jl
~~~~~~~~~~~~~~~~~~

|github1| Creates a diagram showing the analytic means to compute a timing offset.  
Requires installation of the `ExoJulia package <https://github.com/jlustigy/ExoJulia>`_
for carrying out Kepler solver 
before calling. Also requires :py:obj:`PyPlot`, which calls matplotlib, and :py:obj:`CGS.jl`
and :py:obj:`regress.jl`.  For example, from Julia prompt:

.. code-block:: julia

   julia> include("/PathToExoJulia/exojulia.jl")

where `PathToExoJulia` indicates the path where it is installed.

This is followed by:

.. code-block:: julia
   
   julia> include("offset_diagram2.jl")

timing_offset2.jl
~~~~~~~~~~~~~~~~~

|github2| Creates simulated timing offsets for planet-planet
events, and then fits with analytic formulae to recover the eccentricities
and make a plot for the paper.  Same dependencies as :py:obj:`offset_diagram2.jl`.

monte_carlo_circular.jl
~~~~~~~~~~~~~~~~~~~~~~~

|github3| Computes the probability of density and
coplanarity of TRAPPIST-1 based upon the measured transit periods and
the observed durations plus uncertainties (as well as the measured
planet-star radius ratios).  Uses :py:obj:`CGS.jl` and  
`JLD <https://github.com/JuliaIO/JLD.jl>`_, which saves variables in HDF files.

plot_circular.jl
~~~~~~~~~~~~~~~~

|github4| Creates plots of the posterior distribution
of density and coplanarity parameter for TRAPPIST-1.  Uses the
`Optim package <https://github.com/JuliaNLSolvers/Optim.jl>`_ for optimizing a fit to
the distribution, as well as :py:obj:`PyPlot` and :py:obj:`JLD`.  Requires first running 
:py:obj:`monte_carlo_circular.jl` and then reading in the :py:obj:`.jld` file (the name
of the :py:obj:`.jld` file contains the number of simulations carried out;  this
needs to be modified in :py:obj:`plot_circular.jld`).

.. role:: raw-html(raw)
   :format: html

.. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/dynamics/"><i class="fa fa-github" aria-hidden="true"></i></a>`
.. |github1| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/dynamics/offset_diagram2.jl"><i class="fa fa-github" aria-hidden="true"></i></a>`
.. |github2| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/dynamics/timing_offset2.jl"><i class="fa fa-github" aria-hidden="true"></i></a>`
.. |github3| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/dynamics/monte_carlo_circular.jl"><i class="fa fa-github" aria-hidden="true"></i></a>`
.. |github4| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/dynamics/plot_circular.jl"><i class="fa fa-github" aria-hidden="true"></i></a>`
