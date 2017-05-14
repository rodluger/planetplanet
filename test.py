import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
import re
import numpy as np
np.seterr(invalid = 'ignore')
color = pl.get_cmap('plasma')
TOL = 1e-10

def Roots(planet, ellipse):
  '''

  '''
  
  # Get the params
  a = ellipse.a
  b = ellipse.b
  xE = ellipse.x0
  yE = ellipse.y0
  r = planet.r
  
  # Get the coefficients
  A = a ** 2 / b ** 2 - 1
  B = -2 * a ** 2 * xE / b ** 2
  C = r ** 2 - yE ** 2 - a ** 2 + a ** 2 / b ** 2 * xE ** 2
  D = 4 * yE ** 2 * a ** 2 / b ** 2
  c4 = A ** 2
  c3 = 2 * A * B
  c2 = 2 * A * C + B ** 2 + D
  c1 = 2 * B * C - 2 * D * xE
  c0 = C ** 2 - (b ** 2 - xE ** 2) * D
  
  # Solve numerically
  roots = np.roots([c4, c3, c2, c1, c0])
  
  # Filter out imaginary roots
  good_roots = []
  for x in roots:
    if (x.imag == 0):
      good_roots.append(x.real)
      
  return good_roots

class Planet(object):
  '''
  
  '''
  
  def __init__(self, x0, y0, r):
    '''
    
    '''
    
    self.x0 = x0
    self.y0 = y0
    self.r = r
    v1 = (self.x0 - self.r, self.y0)
    v2 = (self.x0 + self.r, self.y0)
    self.vertices = [v1, v2]
    self.xmin = self.x0 - self.r
    self.xmax = self.x0 + self.r
  
  def get_latitude(self, x, y, theta):
    '''
  
    '''
    
    # Where is the terminator for this value of `y`?
    if (theta <= 0):
      xterm = self.x0 - (np.abs(np.sin(theta))) * np.sqrt(self.r ** 2 - (y - self.y0) ** 2)
    else:
      xterm = self.x0 + (np.abs(np.sin(theta))) * np.sqrt(self.r ** 2 - (y - self.y0) ** 2)
    
    # If we are beyond the terminator, return the latitude of the pole
    # TODO: Allow nightside temperature gradients
    if (x > xterm):
      return np.pi / 2
    
    # Otherwise, compute the latitude. We will solve
    # a quadratic equation for z = sin(lat) **  2  
    alpha = (x - self.x0) / (self.r * np.sin(theta))
    beta = 1 / np.tan(theta)
    gamma = (y - self.y0) / self.r
    c1 = 4 * alpha ** 2 * beta ** 2
    c2 = 1 + beta ** 2
    c3 = alpha ** 2 + gamma ** 2 + beta ** 2
    b = (c1 - 2 * c2 * c3) / c2 ** 2
    c = (c3 ** 2 - c1) / c2 ** 2
    z0 = (-b + np.sqrt(b ** 2 - 4 * c)) / 2
    z1 = (-b - np.sqrt(b ** 2 - 4 * c)) / 2
    
    # We have two possible solutions. But only one is on the
    # observer's side of the planet. TODO: Speed this up?
    for z in (z0, z1):
      if (z >= 0) and (z <= 1):
        a = self.r * np.abs(np.sqrt(z))
        b = a * np.abs(np.sin(theta))
        x0 = self.x0 - self.r * np.sqrt(1 - z) * np.cos(theta)
        y0 = self.y0
        dx = (b / a) * np.sqrt(a ** 2 - (y - y0) ** 2)
        xlimb = self.r * np.sqrt(1 - z) * np.sin(theta) * np.tan(theta)
        if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
          xmin = x0 - b
        else:
          xmin = x0 - xlimb
        if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
          xmax = x0 + b
        else:
          xmax = x0 - xlimb
        if (x >= xmin - TOL) and (x <= xmax + TOL): 
          if ((np.abs(x - (x0 + dx)) < TOL) or (np.abs(x - (x0 - dx)) < TOL)):
            return np.arcsin(np.sqrt(z))
    
    return np.nan
  
  def val_upper(self, x):
    '''
    
    '''
    
    A = self.r ** 2 - (x - self.x0) ** 2
    if hasattr(x, '__len__'):
      A[np.abs(A) < TOL] = 0
      A[np.where((x > self.xmax) | (x < self.xmin))] = np.nan
    else:
      if np.abs(A) < TOL:
        A = 0
      if (x > self.xmax) or (x < self.xmin):
        return np.nan
    return self.y0 + np.sqrt(A)

  def val_lower(self, x):
    '''
    
    '''
    
    A = self.r ** 2 - (x - self.x0) ** 2
    if hasattr(x, '__len__'):
      A[np.abs(A) < TOL] = 0
      A[np.where((x > self.xmax) | (x < self.xmin))] = np.nan
    else:
      if np.abs(A) < TOL:
        A = 0
      if (x > self.xmax) or (x < self.xmin):
        return np.nan
    return self.y0 - np.sqrt(A)
  
  def int_upper(self, x0, x1):
    '''
    
    '''
    
    # Check if out of bounds
    if (x0 > self.xmax) or (x0 < self.xmin) or (x1 > self.xmax) or (x1 < self.xmin):
      return np.nan
    
    F = [0, 0]
    for i, x in enumerate((x0, x1)):
      y = x - self.x0
      z = np.sqrt((self.r - y) * (self.r + y))
      F[i] = (1 / 2.) * (y * z + self.r ** 2 * np.arctan(y / z)) + self.y0 * x
    return F[1] - F[0]

  def int_lower(self, x0, x1):
    '''
    
    '''
    
    # Check if out of bounds
    if (x0 > self.xmax) or (x0 < self.xmin) or (x1 > self.xmax) or (x1 < self.xmin):
      return np.nan
    
    F = [0, 0]
    for i, x in enumerate((x0, x1)):
      y = x - self.x0
      z = np.sqrt((self.r - y) * (self.r + y))
      F[i] = -(1 / 2.) * (y * z + self.r ** 2 * np.arctan(y / z)) + self.y0 * x
    return F[1] - F[0]
  
  @property
  def curves(self):
    '''
    
    '''
    
    return [self.val_lower, self.val_upper]
  
  @property
  def integrals(self):
    '''
    
    '''
    
    return [self.int_lower, self.int_upper]
    
class Ellipse(object):
  '''
  
  '''
  
  def __init__(self, planet, occultor, latitude, theta):
    '''
    
    '''

    # Generate the ellipse
    self.a = planet.r * np.abs(np.sin(latitude))
    self.b = self.a * np.abs(np.sin(theta))
    self.x0 = planet.x0 - planet.r * np.cos(latitude) * np.cos(theta)
    self.y0 = planet.y0
    self.latitude = latitude
    
    # The x position of the limb
    self.xlimb = planet.r * np.cos(latitude) * np.sin(theta) * np.tan(theta)
    
    # Compute the vertices
    self.vertices = []
    
    # Ellipse-planet intersections
    if self.b ** 2 >= self.xlimb ** 2:
      x = self.x0 - self.xlimb
      y = self.a / self.b * np.sqrt(self.b ** 2 - self.xlimb ** 2)
      self.vertices.append((x, self.y0 - y))
      self.vertices.append((x, self.y0 + y))
  
    # Ellipse x minimum
    if ((theta > 0) and (self.b < self.xlimb)) or ((theta <= 0) and (self.b > self.xlimb)):
      self.vertices.append((self.x0 - self.b, self.y0))
      self.xmin = self.x0 - self.b
    else:
      self.xmin = self.x0 - self.xlimb
      
    # Ellipse x maximum
    if ((theta > 0) and (self.b > -self.xlimb)) or ((theta <= 0) and (self.b < -self.xlimb)):
      self.vertices.append((self.x0 + self.b, self.y0))
      self.xmax = self.x0 + self.b
    else:
      self.xmax = self.x0 - self.xlimb
      
    # Ellipse-occultor intersections
    for root in Roots(occultor, self):
      if ((theta > 0) and (root > self.x0 - self.xlimb)) or ((theta <= 0) and (root < self.x0 - self.xlimb)):
        eup = self.val_upper(root)
        elo = self.val_lower(root)
        oup = occultor.val_upper(root)
        olo = occultor.val_lower(root)
        if (np.abs(eup - oup) < 1e-10):
          self.vertices.append((root, eup))
        elif (np.abs(eup - olo) < 1e-10):
          self.vertices.append((root, eup))
        elif (np.abs(elo - oup) < 1e-10):
          self.vertices.append((root, elo))
        elif (np.abs(elo - olo) < 1e-10):
          self.vertices.append((root, elo))
    
  def val_upper(self, x):
    '''
    
    '''
    
    A = self.b ** 2 - (x - self.x0) ** 2
    if hasattr(x, '__len__'):
      A[np.abs(A) < TOL] = 0
      A[np.where((x > self.xmax) | (x < self.xmin))] = np.nan
    else:
      if np.abs(A) < TOL:
        A = 0
      if (x > self.xmax) or (x < self.xmin):
        return np.nan
    return self.y0 + (self.a / self.b) * np.sqrt(A)

  def val_lower(self, x):
    '''
    
    '''
    
    A = self.b ** 2 - (x - self.x0) ** 2
    if hasattr(x, '__len__'):
      A[np.abs(A) < TOL] = 0
      A[np.where((x > self.xmax) | (x < self.xmin))] = np.nan
    else:
      if np.abs(A) < TOL:
        A = 0
      if (x > self.xmax) or (x < self.xmin):
        return np.nan
    return self.y0 - (self.a / self.b) * np.sqrt(A)
  
  def int_upper(self, x0, x1):
    '''
    
    '''
    
    # Check if out of bounds
    if (x0 > self.xmax) or (x0 < self.xmin) or (x1 > self.xmax) or (x1 < self.xmin):
      return np.nan
    
    F = [0, 0]
    for i, x in enumerate((x0, x1)):
      y = x - self.x0
      z = np.sqrt((self.b - y) * (self.b + y))
      F[i] = (self.a / (2 * self.b)) * (y * z + self.b ** 2 * np.arctan(y / z)) + self.y0 * x
    return F[1] - F[0]

  def int_lower(self, x0, x1):
    '''
    
    '''
    
    # Check if out of bounds
    if (x0 > self.xmax) or (x0 < self.xmin) or (x1 > self.xmax) or (x1 < self.xmin):
      return np.nan
    
    F = [0, 0]
    for i, x in enumerate((x0, x1)):
      y = x - self.x0
      z = np.sqrt((self.b - y) * (self.b + y))
      F[i] = -(self.a / (2 * self.b)) * (y * z + self.b ** 2 * np.arctan(y / z)) + self.y0 * x
    return F[1] - F[0]

  @property
  def curves(self):
    '''
    
    '''
    
    return [self.val_lower, self.val_upper]
  
  @property
  def integrals(self):
    '''
    
    '''
    
    return [self.int_lower, self.int_upper]

class Occultor(Planet):
  '''
  
  '''
  
  def __init__(self, r, planet):
    '''
    
    '''
    
    # Initialize
    super(Occultor, self).__init__(0, 0, r)
    
    # Compute vertices of intersection with the planet
    # Adapted from http://mathworld.wolfram.com/Circle-CircleIntersection.html
    d = np.sqrt(planet.x0 ** 2 + planet.y0 ** 2)
    if d <= (self.r + planet.r):
      y = np.sqrt((-d + planet.r - self.r) * (-d - planet.r + self.r) * (-d + planet.r + self.r) * (d + planet.r + self.r)) / (2 * d)
      x = (d ** 2 - planet.r ** 2 + self.r ** 2) / (2 * d)
      cost = -planet.x0 / d
      sint = -planet.y0 / d
      x1 = -x * cost + y * sint
      y1 = -x * sint - y * cost
      x2 = -x * cost - y * sint
      y2 = -x * sint + y * cost

      # We're done!
      self.vertices.append((x1, y1))
      self.vertices.append((x2, y2))


# Set up the figure
fig, ax = pl.subplots(1, figsize = (7,7))
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)
axslider = pl.axes([0.125, 0.035, 0.725, 0.03])
slider = Slider(axslider, r'$\theta$', -np.pi / 2, np.pi / 2, valinit = -np.pi / 8)

# The latitude grid
latitude = np.linspace(0, np.pi / 2, 5)[1:]
def surf_brightness(lat):
  '''
  
  '''
  
  grid = np.concatenate(([0], latitude, [np.pi]))
  i = np.argmax(lat < grid) - 1
  return np.cos(grid[i])

# DEBUG
#pl.close()
#x = np.linspace(0, np.pi / 2 + 0.01, 100)
#pl.plot(x, [surf_brightness(xi) for xi in x])
#pl.show()
#quit()


# Planet (occulted)
planet = Planet(0.5, -1.25, 1.25)

# Occultor (always at the origin)
occultor = Occultor(3., planet)

# Arrays for plotting
xp = np.linspace(planet.x0 - planet.r, planet.x0 + planet.r, 1000)
yp0 = planet.val_lower(xp)
yp1 = planet.val_upper(xp)
xo = np.linspace(-occultor.r, occultor.r, 1000)
yo0 = occultor.val_lower(xo)
yo1 = occultor.val_upper(xo)

def compute(theta):
  '''
  
  '''
  
  flux = []
  
  # Avoid the singular point
  if theta == 0:
    theta = TOL
  
  # Compute the ellipses
  ellipses = [Ellipse(planet, occultor, lat, theta) for lat in latitude]
  shapes = [planet, occultor] + ellipses
  
  # Get the curves and their integrals
  curves = []
  integrals = []
  for shape in shapes:
    curves.extend(shape.curves)
    integrals.extend(shape.integrals)
  
  # Compute the vertices inside both the occultor and the planet and sort them
  vertices = []
  for shape in shapes:
    vertices.extend(shape.vertices)
  vin = []
  for v in vertices:
    x, y = v
    if (x ** 2 + y ** 2 <= occultor.r ** 2 + TOL) and ((x - planet.x0) ** 2 + (y - planet.y0) ** 2 <= planet.r ** 2 + TOL) :
      vin.append((x, y))
  vertices = sorted(list(set(vin)))
  
  # Get tuples of integration limits
  limits = list(zip(vertices[:-1], vertices[1:]))

  # Loop over all the sub-regions and compute their fluxes
  for i, lim in enumerate(limits):
    
    # The integration limits
    (xleft, _), (xright, _) = lim
    
    # Check if they are identical
    if xright - TOL <= xleft + TOL:
      continue
      
    # Perturb them for numerical stability, as we
    # need to ensure the functions are finite valued
    # at the limits.
    xleft += TOL
    xright -= TOL
    
    # Bisect the limits. Find the boundary functions that are
    # finite valued at this point and sort their integrals 
    # in order of increasing function value.
    x = (xleft + xright) / 2
    vals = []
    ints = []
    for curve, integral in zip(curves, integrals):
      y = curve(x)
      if np.isfinite(y) and (x ** 2 + y ** 2 <= occultor.r ** 2 + TOL) and ((x - planet.x0) ** 2 + (y - planet.y0) ** 2 <= planet.r ** 2 + TOL):
        vals.append(y)
        ints.append(integral)
    order = np.argsort(np.array(vals))
    ints = [ints[j] for j in order]
    vals = [vals[j] for j in order]
    
    # The areas are just the difference of successive integrals
    for int1, int0, y1, y0 in zip(ints[1:], ints[:-1], vals[1:], vals[:-1]):
      area = int1(xleft, xright) - int0(xleft, xright)
      
      # Get the latitude of the midpoint and
      # the corresponding surface brightness
      y = (y1 + y0) / 2
      lat = planet.get_latitude(x, y, theta)
      f = surf_brightness(lat)
      flux.append(f * area)

  # Total flux
  flux = np.sum(flux)
    
  return flux, ellipses, vertices

def update(theta):
  '''
  
  '''
  
  flux, ellipses, vertices = compute(theta)
  
  # Set up plot
  ax.clear()
  ax.set_xlim(planet.x0 - 1.25 * planet.r, planet.x0 + 1.25 * planet.r)
  ax.set_ylim(planet.y0 - 1.25 * planet.r, planet.y0 + 1.25 * planet.r)
  
  # Plot the planet
  ax.plot(xp, yp0, color = 'k')
  ax.plot(xp, yp1, color = 'k')
  
  # Plot the occultor
  ax.plot(xo, yo0, color = 'k', zorder = 99)
  ax.plot(xo, yo1, color = 'k', zorder = 99)
  
  # Plot the ellipses
  for ellipse in ellipses:

    # First identify and remove points behind the limb
    x = np.linspace(ellipse.x0 - ellipse.b, ellipse.x0 + ellipse.b, 1000)
    if theta > 0:
      x[x < ellipse.x0 - ellipse.xlimb] = np.nan
    else:
      x[x > ellipse.x0 - ellipse.xlimb] = np.nan

    ax.plot(x, ellipse.val_lower(x), lw = 2, color = color(np.cos(ellipse.latitude)))
    ax.plot(x, ellipse.val_upper(x), lw = 2, color = color(np.cos(ellipse.latitude)))

  # Plot the vertices  
  for v in vertices:
    ax.plot(v[0], v[1], 'o', color = 'k', ms = 4)
  
  # Report the flux
  ax.set_title('%.4f' % flux)
    
  # Re-draw!
  fig.canvas.draw_idle()
  
slider.on_changed(update)
update(-np.pi / 8)

pl.show()
