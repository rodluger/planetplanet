def retrograde_bc():
  '''
  A retrograde `c` occults `b` for 180 minutes.
  
  '''
  
  # Instantiate the Trappist-1 system
  system = Trappist1(sample = False, oversample = 10, airless = True, seed = 1234)

  # Get the next occultation
  t = system.next_occultation(2000, system.b, occultor = system.c)
  system.b.limbdark = [0]
  time = np.arange(t - 0.15, t + 0.15, 5. / 1440.)

  # Compute and plot the light curve
  system.compute(time)
  system.plot_lightcurve(15.)

  # Observe it (one exposure)
  system.observe(stack = 1)
  pl.show()