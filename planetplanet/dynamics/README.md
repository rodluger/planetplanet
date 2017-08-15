dynamics
========

1).  offset_diagram2.jl:  Creates a diagram showing the analytic means
to compute a timing offset.  Requires installation of ExoJulia package
(https://github.com/jlustigy/ExoJulia) for carrying out Kepler solver 
before calling.  Also requires PyPlot, which calls matplotlib, and CGS.jl 
& regress.jl.  For example, from Julia prompt:

julia> include("/PathToExoJulia/exojulia.jl")

where "PathToExoJulia" indicates the path where it is installed.

This is followed by:

julia> include("offset_diagram2.jl")

2).  timing_offset2.jl: Creates simulated timing offsets for planet-planet
events, and then fits with analytic formulae to recover the eccentricities
and make a plot for the paper.  Same dependencies as offset_diagram2.jl.

3).  monte_carlo_circular.jl:  Computes the probability of density and
coplanarity of TRAPPIST-1 based upon the measured transit periods and
the observed durations plus uncertainties (as well as the measured
planet-star radius ratios).  Uses CGS.jl and JLD (https://github.com/JuliaIO/JLD.jl
which saves variables in HDF files).

4).  plot_circular.jl:  Creates plots of the posterior distribution
of density and coplanarity parameter for TRAPPIST-1.  Uses Optim package
(https://github.com/JuliaNLSolvers/Optim.jl) for optimizing a fit to
the distribution, as well as PyPlot and JLD.  Requires first running 
monte_carlo_circular.jl and then reading in the .jld file (the name
of the .jld file contains the number of simulations carried out;  this
needs to be modified in plot_circular.jl).
