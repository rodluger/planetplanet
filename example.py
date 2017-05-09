import planet_planet

def Test():
  '''
  
  '''

  # Defaults
  ro = 1
  xo = 0
  yo = 0
  rp = 1.25
  xp = np.linspace(-3,3,1000)
  xp0 = 0.
  yp = -1.5
  theta = -np.pi / 2
  fday = 1
  fnight = 0.25
  
  # Plot
  fig, ax = pl.subplots(1, figsize = (6, 6))
  fig.subplots_adjust(left = 0.15)
  angles = [-np.pi / 2, -np.pi / 4, 0, np.pi / 4, np.pi / 2]
  colors = ['r', 'g', 'b', 'y', 'k']
  
  for n, theta in enumerate(angles):
    
    # Plot the light curve
    df = np.zeros(1000)
    for i in range(1000):
      df[i] = dFlux(xp[i], yp, rp, xo, yo, ro, theta = theta, fday = fday, fnight = fnight)
    ax.plot(xp, df, color = colors[n], alpha = 0.75)
  
    # Plot the image
    axl = fig.add_axes([0.725, 0.65 - n * 0.125, 0.15, 0.1])
    axl.xaxis.set_visible(False)
    axl.yaxis.set_visible(False)
    axl.axis('off')

    # Draw the planet
    x = np.linspace(xp0 - rp, xp0 + rp, 10000)
    ytop = yp + np.sqrt(rp ** 2 - (x - xp0) ** 2)
    ybot = yp - np.sqrt(rp ** 2 - (x - xp0) ** 2)
    axl.plot(x, ytop, color = 'k', lw = 1)
    axl.plot(x, ybot, color = 'k', lw = 1)
    axl.fill_between(x, ybot, ytop, color = 'lightgray')

    # Draw the night side
    x = np.linspace(xp0, xp0 + rp, 10000)
    ytop = yp + np.sqrt(rp ** 2 - (x - xp0) ** 2)
    ybot = yp - np.sqrt(rp ** 2 - (x - xp0) ** 2)
    axl.plot(x, ytop, color = 'k', lw = 1)
    axl.plot(x, ybot, color = 'k', lw = 1)
    axl.fill_between(x, ybot, ytop, color = 'gray')

    # Draw the terminator
    b = rp * np.sin(theta)
    if b >= 0:
      color = 'gray'
    else:
      color = 'lightgray'
    x = np.linspace(xp0 - rp * b, xp0, 10000)
    ytop = yp + rp * np.sqrt(1 - (x - xp0) ** 2 / b ** 2)
    ybot = yp - rp * np.sqrt(1 - (x - xp0) ** 2 / b ** 2)
    axl.fill_between(x, ybot, ytop, color = color)

    # Draw the occultor
    x = np.linspace(xo - ro, xo + ro, 10000)
    ytop = yo + np.sqrt(ro ** 2 - (x - xo) ** 2)
    ybot = yo - np.sqrt(ro ** 2 - (x - xo) ** 2)
    axl.plot(x, ytop, color = 'k', lw = 1)
    axl.plot(x, ybot, color = 'k', lw = 1)
    axl.fill_between(x, ybot, ytop, color = 'k', alpha = 0.75)
    
    # Limits
    axl.plot([-3.75, -2.5],[-1.25, -1.25], lw = 2, color = colors[n])
    axl.set_xlim(-4.0, 2.0)
    axl.set_ylim(-2.9, 1.1)
  
  ax.set_ylabel('Relative Flux', fontsize = 16)
  ax.set_xlabel('Time', fontsize = 16)
  pl.show()