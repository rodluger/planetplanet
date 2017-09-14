/**
@file eyeball.c
@brief Occultation light curve routines.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>

/**
Computes the blackbody intensity (W / m^2 / m / sr).

@param lambda The wavelength in m
@param T The effective temperature in K
@return The intensity given by Planck's Law

*/
double Blackbody(double lambda, double T) {
  
  // Planck's Law
  double a = 2 * HPLANCK * CLIGHT * CLIGHT / (lambda * lambda * lambda * lambda * lambda);
  double b = HPLANCK * CLIGHT / (lambda * KBOLTZ * T);
  return a / (exp(b) - 1);
}

/**
An implementation of the BATMAN algorithm (Kreidberg 2016) to compute
occultations of radially-symmetric bodies. This applies to transits
across stars and planet-planet occultations in the limb-darkened limit
or of planets at full or new phase.

@param r The occulted body's radius in Earth radii
@param x0 The x coordinate of the occultor relative to the occulted body in Earth radii
@param y0 The y coordinate of the occultor relative to the occulted body in Earth radii
@param ro The radius of the occultor in Earth radii
@param teff The effective temperature of the occulted body
@param distance The distance to the system in parsecs
@param nu The number of limb darkening coefficients
@param nz The number of zenith angle slices
@param nw The size of the wavelength grid
@param u The limb darkening coefficient grid: \a nu coefficients at each wavelength
@param lambda The wavelength grid in m
@param flux The wavelength-dependent light curve computed by this function

*/
void BatmanFlux(double r, double x0, double y0, double ro, double teff, double distance, int nu, int nz, int nw, double u[nu * nw], double lambda[nw], double flux[nw]) {
  
  int i, j, k;  
  double d, d2, ro2, x2;
  double U, V, W;
  double x[nz], A[nz];
  double xavg;
  double area, norm, cosza, y;
  double z, zmin, zmax;
  double B[nw], B0[nw];
  double distance2 = distance * distance * PARSEC * PARSEC;
  
  // Normalize distances by the occulted planet's radius
  x0 /= r;
  y0 /= r;
  ro /= r;
  
  // Pre-compute some stuff
  ro2 = ro * ro;
  d2 = x0 * x0 + y0 * y0;
  d = sqrt(d2);
  
  // Pre-normalize the polynomial limb darkening
  // Equation (E5) in Luger et al. (2017)
  for (j = 0; j < nw; j++) {
    norm = 0;
    for (k = 0; k < nu; k++)
      norm += u[nw * k + j] / ((k + 2) * (k + 3));
    norm = 1 - 2 * norm;
    B0[j] = Blackbody(lambda[j], teff) / norm;
  }
  
  // Integration grid, constant in zenith angle
  if (d - ro <= 0) zmin = 0;
  else zmin = asin(d - ro);
  if (d + ro >= 1) zmax = PI / 2;
  else zmax = asin(d + ro);
  for (i = 0; i < nz; i++) {
    z = (i + 0.) / (nz - 1.) * (zmax - SMALL - zmin) + zmin + SMALL;
    x[i] = sin(z);
  }
  
  // Iterate over each zenith angle slice
  for (i = 0; i < nz; i++) {
      
    // Compute the area below
    if (x[i] <= ro - d) {
      x2 = x[i] * x[i];
      A[i] = PI * x2;
    } else if (x[i] >= ro + d) {
      A[i] = PI * ro2;
    } else {
      x2 = x[i] * x[i];
      U = (d2 + x2 - ro2) / (2 * d * x[i]);
      V = (d2 + ro2 - x2) / (2 * d * ro);
      W = (-d + x[i] + ro) * (d + x[i] - ro) * (d - x[i] + ro) * (d + x[i] + ro); 
      // HACK: Prevent numerical issues
      if (W < 0) W = 0;
      if (U < -1) U = -1;
      else if (U > 1) U = 1;
      if (V < -1 ) V = -1;
      else if (V > 1) V = 1;
      A[i] = x2 * acos(U) + ro2 * acos(V) - 0.5 * sqrt(W);
    }

  }

  // Integrate to compute the flux
  for (i = 1; i < nz; i++) {
  
    // The area of the current segment
    area = (A[i] - A[i - 1]) * r * r;
    
    // Pre-compute mu, the cosine of the zenith angle
    // Note that cos(arcsin(x)) = sqrt(1 - x^2)
    xavg = 0.5 * (x[i] + x[i - 1]);
    cosza = sqrt(1 - xavg * xavg);
    
    // Initialize the intensity grid
    for (j = 0; j < nw; j++)
      B[j] = B0[j];
    
    // Loop over the coefficient order
    for (k = 0; k < nu; k++) {
    
      // The Taylor expansion is in y = 1 - mu
      y = pow(1 - cosza, k + 1);
      
      // Compute the wavelength-dependent intensity
      for (j = 0; j < nw; j++) {
        B[j] -= u[nw * k + j] * B0[j] * y;
      } 
      
    }
    
    // Finally, flux is radiance times area
    for (j = 0; j < nw; j++)
      flux[j] += B[j] * area;
      
  }
  
  // Scale flux by REARTH**2, convert from m^-1 to micron^-1, and
  // divide by the square of the distance to get a radiance (W / m^2 / um)
  for (j = 0; j < nw; j++) 
    flux[j] *= REARTH * REARTH * MICRON / distance2;
  
}

/**
Compares two doubles. Used with qsort.

@param a A double
@param b A double
@return 0 if a == b, 1 if a < b, -1 otherwise

*/
int dblcomp( const void* a, const void* b) {
  double dbl_a = * ( (double*) a );
  double dbl_b = * ( (double*) b );
  if ( dbl_a == dbl_b ) return 0;
  else if ( dbl_a < dbl_b ) return -1;
  else return 1;
}

/**
Compares the y values of two FUNCTIONs. Used with qsort.

@param a A FUNCTION instance
@param b A FUNCTION instance
@return 0 if a.y == b.y, 1 if a.y < b.y, -1 otherwise

*/
int funcomp( const void* a, const void* b) {
  FUNCTION *fun_a = (FUNCTION *)a;
  FUNCTION *fun_b = (FUNCTION *)b;
  if ( fun_a->y == fun_b->y ) return 0;
  else if ( fun_a->y < fun_b->y ) return -1;
  else return 1;
}

/**
Computes the flux at a given zenith angle on an eyeball planet.

@param lambda The wavelength in m
@param za The zenith angle in radians, measured from the hotspot
@param teff The body's equilibrium temperature
@param tnight The night side temperature in K
@return The emitted spectral flux density from that point on the surface

*/
double RadiativeEquilibriumMap(double lambda, double za, double teff, double tnight) {
  
  double temp;
  
  if (za < PI / 2) {
    temp = teff * pow(4 * cos(za), 0.25);
    if (temp < tnight)
      temp = tnight;
  } else {
    temp = tnight;
  }
  
  // Return the flux
  return Blackbody(lambda, temp);
  
}

/**
Computes the flux at a given zenith angle on a limb-darkened planet.

@param j The index of the wavelength grid
@param za The zenith angle in radians, measured from the hotspot
@param teff The body's effective temperature
@param nu The number of limb-darkening coefficients per wavelength bin
@param nw The number of wavelength bins
@param u The limb-darkening coefficient grid
@param B0 The normalization constant for the radiance at this wavelength
@return The emitted spectral flux density from that point on the surface

*/
double LimbDarkenedMap(int j, double za, double teff, int nu, int nw, double u[nu * nw], double B0) {
  
  int k;
  double x;
  double B = B0;
  double cosza = cos(za);
  
  // Loop over the coefficient order
  for (k = 0; k < nu; k++) {
  
    // The Taylor expansion is in (1 - mu)
    x = pow(1 - cosza, k + 1);
    B -= u[nw * k + j] * B0 * x;
    
  } 
  
  // Return the flux
  return B;
  
}

/**
Computes the real roots of the circle-ellipse intersection quartic equation. These roots
are the x coordinates of the intersection point(s) relative to the center of the circle.

@param a The semi-major axis of the ellipse, which is aligned with the y axis
@param b The semi-minor axis of the ellipse, which is aligned with the x axis
@param xE The x coordinate of the center of the ellipse
@param yE The y coordinate of the center of the ellipse
@param xC The x coordinate of the center of the circle
@param yC The y coordinate of the center of the ellipse
@param r The radius of the circle
@param roots The x coordinates of the point(s) of intersection

*/
void GetRootsGSL(double a, double b, double xE, double yE, double xC, double yC, double r, double roots[4]) {
  /*

  */
  
  int i, j;
  double A, B, C, D;
  double z[8];
  double r2 = r * r;
  double a2 = a * a;
  double b2 = b * b;
  double a2b2 = a2 / b2;
  double x0 = xE - xC;
  double y0 = yE - yC;
  double y2 = y0 * y0;
  double x2 = x0 * x0;
  int iErr;
  
  // Initialize the roots
  roots[0] = NAN;
  roots[1] = NAN;
  roots[2] = NAN;
  roots[3] = NAN;
   
  // Get the quartic coefficients
  A = a2b2 - 1.;
  B = -2. * x0 * a2b2;
  C = r2 - y2 - a2 + a2b2 * x2;
  D = 4. * y2 * a2b2;
  double c[5] = { 
                  C * C - (b2 - x2) * D, 
                  2. * B * C - 2. * D * x0, 
                  2. * A * C + B * B + D,
                  2. * A * B, 
                  A * A 
                };
  
  // Turn off the default error handler
  gsl_set_error_handler_off();
  
  // Allocate memory
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(5);
  
  // Solve the polynomial
  iErr = gsl_poly_complex_solve (c, 5, w, z);
  
  // TODO: These errors pop up here and there. I'm issuing a warning
  // and ignoring them so we don't cause large jobs to fail. This needs
  // to be investigated in the future. I suspect it's due to occultations
  // of planets that are very close to quadrature.
  if (iErr == GSL_EFAILED) {
    printf("WARNING: Circle/ellipse intersection point finder failed to converge.\n");
    return;
  }
  
  // Free the memory
  gsl_poly_complex_workspace_free (w);
  
  // Get the real roots (up to 4)
  j = 0;
  for (i = 0; i < 8; i+=2) {
    if (!(isnan(z[i])) && (fabs(z[i + 1]) < MAXIM)) {
      roots[j] = z[i] + xC;
      j += 1;
    }
  }

}

/**
Computes the zenith angle of a point (x, y) on the projected disk of a planet.

@param x The x coordinate of the point
@param y The y coordinate of the point
@param r The radius of the planet
@param theta The phase angle of the eyeball
@return The zenith angle in radians

*/
double ZenithAngle(double x, double y, double r, double theta) {
    
  // Normalize
  x = x / r;
  y = y / r;
  double x2 = x * x;
  double y2 = y * y;
  double xterm, z;
  
  // Are we dealing with circles?
  if ((fabs(fabs(theta) - PI / 2)) < SMALL) {
    if (theta > 0)
      return asin(sqrt(x2 + y2));
    else
      return PI - asin(sqrt(x2 + y2));
  }
  
  // This is a solution to a quadratic equation in z = sin(za) **  2 
  z = 0.5 * ((1 - 2 * x2 - y2) * cos(2 * theta) + 2 * x * sqrt(1 - x2 - y2) * sin(2 * theta) + y2 + 1);
  xterm = sin(theta) * sqrt(fabs(1 - y2));
  if (fabs(theta) <= PI / 2) {
    if (x <= xterm)
      return asin(sqrt(z));
    else
      return PI - asin(sqrt(z));
  } else {
    if (x >= -xterm)
      return asin(sqrt(z));
    else
      return PI - asin(sqrt(z));
  }

}

/**
Computes the blackbody intensity at the center of each zenith_angle slice 
evaluated at a given array of wavelengths.

@param tnight The night side temperature in K (eyeball limit)
@param teff The effective temperature of the planet (blackbody limit)
@param maptype The code for the radiance map
@param radiancemap A pointer to the RADIANCEMAP function
@param nz The number of zenith angle slices
@param zenithgrid The zenith angle grid
@param nw The number of wavelength points
@param lambda The wavelength array
@param nu The number of limb darkening coefficients
@param u The limb darkening coefficient grid: \a nu coefficients at each wavelength
@param B0 The wavelength-dependent normalization for the limb-darkened case
@param B The blackbody intensity grid

*/
void SurfaceIntensity(double tnight, double teff, int maptype, RADIANCEMAP radiancemap, int nz, double zenithgrid[nz], int nw, double lambda[nw], int nu, double u[nu * nw], double B0[nw], double B[nz + 1][nw]) {
  
  int i, j;
  double zenith_angle;
  
  // Loop over each slice
  for (i = 0; i < nz + 1; i++) {
    
    // Get the zenith_angle halfway between the `zenithgrid` points.
    // In the regular grid, B[0] is the intensity of the region between the 
    // sub-stellar point and zenithgrid[0] (the location of the first zenith_angle ellipse); 
    // B[1] is the intensity between zenithgrid[1] and zenithgrid[2], and so forth.
    // B[nz] is the intensity between zenithgrid[nz - 1] and the anti-stellar point.
    if (i == 0) {
      zenith_angle = zenithgrid[0] - 0.5 * (zenithgrid[1] - zenithgrid[0]);
      if (zenith_angle < 0) zenith_angle = 0;
    } else if (i == nz) {
      zenith_angle = zenithgrid[nz - 1] + 0.5 * (zenithgrid[nz - 1] - zenithgrid[nz - 2]);
      if (zenith_angle > PI) zenith_angle = PI;
    } else
      zenith_angle = 0.5 * (zenithgrid[i] + zenithgrid[i - 1]);

    // Get the flux      
    for (j = 0; j < nw; j++) {
      if (maptype == MAP_ELLIPTICAL_DEFAULT)
        B[i][j] = RadiativeEquilibriumMap(lambda[j], zenith_angle, teff, tnight);
      else if (maptype == MAP_RADIAL_DEFAULT)
        B[i][j] = LimbDarkenedMap(j, zenith_angle, teff, nu, nw, u, B0[j]);
      else
        B[i][j] = radiancemap(lambda[j], zenith_angle);
    }
  } 
}

/**
Evaluates a FUNCTION instance.

@param x The point at which to evaluate the function
@param function A FUNCTION instance
@return The value of the function at \a x

*/
double curve(double x, FUNCTION function) {
  return function.curve(x, function.ellipse);
}

/**
Evaluates the definite integral of a FUNCTION instance.

@param x0 The lower integration limit
@param x1 The upper integration limit
@param function A FUNCTION instance
@param oob Flag set to 1 if the integration limits are out of bounds
@return The definite integral of the function from \a x0 to \a x1

*/
double integral(double x0, double x1, FUNCTION function, int *oob) {
  return function.integral(x0, x1, function.ellipse, oob);
}

/**
The upper curve of an ellipse.

@param x The point at which to evaluate the curve
@param ellipse An ELLIPSE instance
@return The value of the upper curve of the ellipse at \a x.

*/
double fupper(double x, ELLIPSE *ellipse) {
  double A;
  if (ellipse->circle) {
    A = ellipse->r * ellipse->r - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < TINY) A = 0;
    else if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 + sqrt(A);
  } else {
    A = ellipse->b * ellipse->b - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < TINY) A = 0;
    else if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 + (ellipse->a / ellipse->b) * sqrt(A);
  }
}

/**
The lower curve of an ellipse.

@param x The point at which to evaluate the curve
@param ellipse An ELLIPSE instance
@return The value of the lower curve of the ellipse at \a x.

*/
double flower(double x, ELLIPSE *ellipse) {
  double A;
  if (ellipse->circle) {
    A = ellipse->r * ellipse->r - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < TINY) A = 0;
    else if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 - sqrt(A);
  } else {
    A = ellipse->b * ellipse->b - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < TINY) A = 0;
    else if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 - (ellipse->a / ellipse->b) * sqrt(A);
  }
}

/**
The area under the upper curve of an ellipse.

@param xL The left (lower) integration limit
@param xR The right (upper) integration limit
@param ellipse An ELLIPSE instance
@param oob Flag set to 1 if the integration limits are out of bounds
@return The value of the definite integral of the upper curve of the ellipse

*/
double iupper(double xL, double xR, ELLIPSE *ellipse, int *oob) {
  /*
  
  */
  
  double yL,yR,zL,zR,FL,FR;
  
  // Catch NaNs silently
  if ((xL > ellipse->xmax) || (xL < ellipse->xmin)) *oob = 1;
  if ((xR > ellipse->xmax) || (xR < ellipse->xmin)) *oob = 1;
  
  yL = xL - ellipse->x0;
  yR = xR - ellipse->x0;
  if (ellipse->circle) {
    zL = sqrt((ellipse->r - yL) * (ellipse->r + yL));
    if (isnan(zL)) zL = 0;
    zR = sqrt((ellipse->r - yR) * (ellipse->r + yR));
    if (isnan(zR)) zR = 0;
    FL = 0.5 * (yL * zL + ellipse->r * ellipse->r * atan(yL / zL)) + ellipse->y0 * xL;
    FR = 0.5 * (yR * zR + ellipse->r * ellipse->r * atan(yR / zR)) + ellipse->y0 * xR;
  } else {
    zL = sqrt((ellipse->b - yL) * (ellipse->b + yL));
    if (isnan(zL)) zL = 0;
    zR = sqrt((ellipse->b - yR) * (ellipse->b + yR));
    if (isnan(zR)) zR = 0;
    FL = (ellipse->a / (2 * ellipse->b)) * (yL * zL + ellipse->b * ellipse->b * atan(yL / zL)) + ellipse->y0 * xL;
    FR = (ellipse->a / (2 * ellipse->b)) * (yR * zR + ellipse->b * ellipse->b * atan(yR / zR)) + ellipse->y0 * xR;
  }
  return FR - FL;
}

/**
The area under the lower curve of an ellipse.

@param xL The left (lower) integration limit
@param xR The right (upper) integration limit
@param ellipse An ELLIPSE instance
@param oob Flag set to 1 if the integration limits are out of bounds
@return The value of the definite integral of the lower curve of the ellipse

*/
double ilower(double xL, double xR, ELLIPSE *ellipse, int *oob) {
  /*
  
  */
  
  double yL,yR,zL,zR,FL,FR;
  
  // Catch NaNs silently
  if ((xL > ellipse->xmax) || (xL < ellipse->xmin)) *oob = 1;
  if ((xR > ellipse->xmax) || (xR < ellipse->xmin)) *oob = 1;
    
  yL = xL - ellipse->x0;
  yR = xR - ellipse->x0;
  if (ellipse->circle) {
    zL = sqrt((ellipse->r - yL) * (ellipse->r + yL));
    if (isnan(zL)) zL = 0;
    zR = sqrt((ellipse->r - yR) * (ellipse->r + yR));
    if (isnan(zR)) zR = 0;
    FL = -0.5 * (yL * zL + ellipse->r * ellipse->r * atan(yL / zL)) + ellipse->y0 * xL;
    FR = -0.5 * (yR * zR + ellipse->r * ellipse->r * atan(yR / zR)) + ellipse->y0 * xR;
  } else {
    zL = sqrt((ellipse->b - yL) * (ellipse->b + yL));
    if (isnan(zL)) zL = 0;
    zR = sqrt((ellipse->b - yR) * (ellipse->b + yR));
    if (isnan(zR)) zR = 0;
    FL = -(ellipse->a / (2 * ellipse->b)) * (yL * zL + ellipse->b * ellipse->b * atan(yL / zL)) + ellipse->y0 * xL;
    FR = -(ellipse->a / (2 * ellipse->b)) * (yR * zR + ellipse->b * ellipse->b * atan(yR / zR)) + ellipse->y0 * xR;
  }
  return FR - FL;
}

/**
Computes the full geometry for a single zenith angle ellipse. Computes the 
functional form of the ellipse, its integral, and all points of intersection
with the limb of the occulted planet and the limb of the occultor. Adds the
vertices and functions to the running lists.

@param zenith_angle The zenith angle in radians
@param r The radius of the occulted planet
@param no The number of occulting bodies
@param x0 The x coordinate of the center of each occulting body
@param y0 The y coordinate of the center of each occulting body
@param ro The radius of each occulting body
@param quarticsolver Which quartic solver to use?
@param theta The phase angle of the eyeball
@param maxvertices Maximum number of vertices in the problem
@param maxfunctions Maximum number of functions in the problem
@param vertices The list of vertices
@param v The number of vertices so far
@param functions The list of functions
@param f The number of functions so far

*/
void AddZenithAngleEllipse(double zenith_angle, double r, int no, double x0[no], double y0[no], double ro[no], double theta, int quarticsolver, int maxvertices, int maxfunctions, double *vertices, int *v, FUNCTION *functions, int *f){
    
  int i, j;
  double xlimb, x, y;
  double roots[4];
  double ro2[no];
  for (i = 0; i < no; i++) {
    ro2[i] = ro[i] * ro[i];
    ro2[i] += TINY * ro2[i];
  }
  
  // Check if this ellipse is at all visible
  double zmin, zmax;
  if (theta <= 0) {
    zmax = PI;
    if (theta <= -PI / 2)
      zmin = PI + theta;
    else
      zmin = -theta;
  } else {
    zmin = 0;
    if (theta >= PI / 2)
      zmax = theta;
    else
      zmax = PI - theta;
  }
  if ((zenith_angle < zmin) || (zenith_angle > zmax)) return;
  
  // Construct the ellipse
  ELLIPSE *ellipse;
  ellipse = malloc(sizeof(ELLIPSE)); 
  ellipse->circle = 0;
  ellipse->a = r * fabs(sin(zenith_angle));
  ellipse->b = ellipse->a * fabs(sin(theta));
  ellipse->x0 = -r * cos(zenith_angle) * cos(theta);
  ellipse->y0 = 0;
  
  // The x position of the intersection with the occulted planet
  // limb, relative to the ellipse center
  xlimb = -r * cos(zenith_angle) * sin(theta) * tan(theta);

  // Ellipse x minimum
  if (((theta >= 0) && (zenith_angle <= theta)) || ((theta <= 0) && (zenith_angle >= -theta))) {
    ellipse->xmin = ellipse->x0 - ellipse->b;
    // Is this point inside __at least__ one occultor?
    for (i = 0; i < no; i++) {
      if ((ellipse->xmin - x0[i]) * (ellipse->xmin - x0[i]) + y0[i] * y0[i] < ro2[i]) {
        vertices[(*v)++] = ellipse->xmin;
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
        break;
      }
    }
  } else {
    ellipse->xmin = ellipse->x0 + xlimb;
  }
    
  // Ellipse x maximum
  if (((theta >= 0) && (zenith_angle <= PI - theta)) || ((theta <= 0) && (zenith_angle >= PI + theta))) {
    ellipse->xmax = ellipse->x0 + ellipse->b;
    // Is this point inside __at least__ one occultor?
    for (i = 0; i < no; i++) {
      if ((ellipse->xmax - x0[i]) * (ellipse->xmax - x0[i]) + y0[i] * y0[i] < ro2[i]) {
        vertices[(*v)++] = ellipse->xmax;
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
        break;    
      }
    }
  } else {
    ellipse->xmax = ellipse->x0 + xlimb;
  }
  
  // Are the limb vertices inside __at least__ one occultor? 
  if (ellipse->b >= fabs(xlimb)) {
    x = ellipse->x0 + xlimb;
    y = fupper(x, ellipse);
    for (i = 0; i < no; i++) {
      if ((x - x0[i]) * (x - x0[i]) + (y - y0[i]) * (y - y0[i]) < ro2[i]) {
        vertices[(*v)++] = x;
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
        break;    
      }
    }
    y = flower(x, ellipse);
    for (i = 0; i < no; i++) {
      if ((x - x0[i]) * (x - x0[i]) + (y - y0[i]) * (y - y0[i]) < ro2[i]) {
        vertices[(*v)++] = x;
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
        break;    
      }
    }
  }
  
  // Ellipse-occultor vertices
  for (i = 0; i < no; i++) {
    
    // Is the occulted planet entirely within the occultor?
    if (sqrt(x0[i] * x0[i] + y0[i] * y0[i]) + r + SMALL < ro[i]) continue;
    
    // Solve the quartic
    if (quarticsolver == QGSL)
      GetRootsGSL(ellipse->a, ellipse->b, ellipse->x0, ellipse->y0, x0[i], y0[i], ro[i], roots);
    else {
      // TODO: Implement more solvers
      printf("ERROR: Invalid quartic solver.\n");
      abort();
    }
    
    // Keep only real roots within the visible bounds of the ellipse
    for (j = 0; j < 4; j++) {
      if (!(isnan(roots[j])) && (roots[j] > ellipse->xmin) && (roots[j] < ellipse->xmax)) {
        vertices[(*v)++] = roots[j];
        if (*v > maxvertices) {
            printf("ERROR: Maximum number of vertices exceeded.\n");
            abort();
          }
      }
    }
    
  }
  
  // Finally, the boundary curves and their integrals
  functions[*f].ellipse = ellipse;
  functions[*f].curve = fupper;
  functions[(*f)++].integral = iupper;
  functions[*f].ellipse = ellipse;
  functions[*f].curve = flower;
  functions[(*f)++].integral = ilower;
  if (*f > maxfunctions) {
    printf("ERROR: Maximum number of functions exceeded.\n");
    abort();
  }
  
}

/**
Computes the full geometry for a single zenith angle circle. This is the limiting
case of \a AddZenithAngleEllipse for limb-darkened planets or those at full/new 
phase. Computes the functional form of the circle, its integral, and all points of intersection
with the limb of the occultor. Adds the vertices and functions to the running lists.

@param zenith_angle The zenith angle in radians
@param r The radius of the occulted planet
@param no The number of occulting bodies
@param x0 The x coordinate of the center of each occulting body
@param y0 The y coordinate of the center of each occulting body
@param ro The radius of each occulting body
@param maxvertices Maximum number of vertices in the problem
@param maxfunctions Maximum number of functions in the problem
@param vertices The list of vertices
@param v The number of vertices so far
@param functions The list of functions
@param f The number of functions so far

*/
void AddZenithAngleCircle(double zenith_angle, double r, int no, double x0[no], double y0[no], double ro[no], int maxvertices, int maxfunctions, double *vertices, int *v, FUNCTION *functions, int *f){

  int i;
  double ro2[no];
  double d, A, x, y, frac, cost, sint;
  for (i = 0; i < no; i++) {
    ro2[i] = ro[i] * ro[i];
    ro2[i] += TINY * ro2[i];
  }
  ELLIPSE *ellipse;
  ellipse = malloc(sizeof(ELLIPSE)); 
  
  // Initialize the circle
  ellipse->circle = 1;
  ellipse->r = r * fabs(sin(zenith_angle));
  ellipse->x0 = 0;
  ellipse->y0 = 0;
  
  // Circle x minimum
  ellipse->xmin = ellipse->x0 - ellipse->r;
  // Is this point inside __at least__ one occultor?
  for (i = 0; i < no; i++) {
    if ((ellipse->xmin - x0[i]) * (ellipse->xmin - x0[i]) + y0[i] * y0[i] < ro2[i]) {
      vertices[(*v)++] = ellipse->xmin;
      if (*v > maxvertices) {
        printf("ERROR: Maximum number of vertices exceeded.\n");
        abort();
      }
      break;    
    }
  }

  // Circle x maximum
  ellipse->xmax = ellipse->x0 + ellipse->r;
  // Is this point inside __at least__ one occultor?
  for (i = 0; i < no; i++) {
    if ((ellipse->xmax - x0[i]) * (ellipse->xmax - x0[i]) + y0[i] * y0[i] < ro2[i]) {
      vertices[(*v)++] = ellipse->xmax;
      if (*v > maxvertices) {
        printf("ERROR: Maximum number of vertices exceeded.\n");
        abort();
      }
      break;    
    }
  }
    
  // Circle-occultor vertices
  for (i = 0; i < no; i++) {
    
    // Is the occulted planet entirely within the occultor?
    if (sqrt(x0[i] * x0[i] + y0[i] * y0[i]) + r + SMALL < ro[i]) continue;
    
    // Adapted from http://mathworld.wolfram.com/Circle-CircleIntersection.html
    d = sqrt(x0[i] * x0[i] + y0[i] * y0[i]);
    if (d < (ro[i] + ellipse->r)) {
      A = (-d + ellipse->r - ro[i]) * (-d - ellipse->r + ro[i]) * (-d + ellipse->r + ro[i]) * (d + ellipse->r + ro[i]);
      if (A >= 0) {
        y = sqrt(A) / (2 * d);
        x = sqrt(ellipse->r * ellipse->r - y * y);
        frac = y0[i] / x0[i];
        cost = 1 / sqrt(frac * frac + 1);
        sint = frac * cost;
        if (x0[i] < 0) {
          cost *= -1;
          sint *= -1;
        }
        vertices[(*v)++] = x * cost + y * sint;
        vertices[(*v)++] = x * cost - y * sint;
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
      }
    }
      
  }
  
  // Finally, the boundary curves and their integrals
  functions[*f].ellipse = ellipse;
  functions[*f].curve = fupper;
  functions[(*f)++].integral = iupper;
  functions[*f].ellipse = ellipse;
  functions[*f].curve = flower;
  functions[(*f)++].integral = ilower;
  if (*f > maxfunctions) {
    printf("ERROR: Maximum number of functions exceeded.\n");
    abort();
  }
  
}

/**
Computes the geometry for all occulting bodies. Adds relevant vertices
and functions to the running lists.

@param r The radius of the occulted planet
@param no The number of occulting bodies
@param x0 The x coordinate of the center of each occulting body
@param y0 The y coordinate of the center of each occulting body
@param ro The radius of each occulting body
@param maxvertices Maximum number of vertices in the problem
@param maxfunctions Maximum number of functions in the problem
@param vertices The list of vertices
@param v The number of vertices so far
@param functions The list of functions
@param f The number of functions so far

*/
void AddOccultors(double r, int no, double x0[no], double y0[no], double ro[no], int maxvertices, int maxfunctions, double *vertices, int *v, FUNCTION *functions, int *f) {

  int i;
  double r2 = r * r;
  r2 += TINY * r2;
  ELLIPSE *occultor[no];
   
  for (i = 0; i < no; i++) {
    
    // Allocate memory
    occultor[i] = malloc(sizeof(ELLIPSE));
    
    // Initialize the occultor
    occultor[i]->r = ro[i];
    occultor[i]->x0 = x0[i];
    occultor[i]->y0 = y0[i];
    occultor[i]->circle = 1;
    
    // Occultor x minimum
    occultor[i]->xmin = x0[i] - occultor[i]->r;
    if (occultor[i]->xmin * occultor[i]->xmin + occultor[i]->y0 * occultor[i]->y0 < r2) {
      vertices[(*v)++] = occultor[i]->xmin;  
      if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }  
    }
    
    // Occultor x maximum
    occultor[i]->xmax = x0[i] + occultor[i]->r;
    if (occultor[i]->xmax * occultor[i]->xmax + occultor[i]->y0 * occultor[i]->y0 < r2) {
      vertices[(*v)++] = occultor[i]->xmax;    
      if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
    }

    // Finally, the boundary curves and their integrals
    functions[*f].ellipse = occultor[i];
    functions[*f].curve = fupper;
    functions[(*f)++].integral = iupper;
    functions[*f].ellipse = occultor[i];
    functions[*f].curve = flower;
    functions[(*f)++].integral = ilower;
    if (*f > maxfunctions) {
      printf("ERROR: Maximum number of functions exceeded.\n");
      abort();
    }
  
  }
  
}

/**
Computes the geometry for the occulted body. Adds the vertices of intersection
with each of the occultors and the functional form of the planet disk to the running
lists.

@param r The radius of the occulted planet
@param no The number of occulting bodies
@param x0 The x coordinate of the center of each occulting body
@param y0 The y coordinate of the center of each occulting body
@param ro The radius of each occulting body
@param maxvertices Maximum number of vertices in the problem
@param maxfunctions Maximum number of functions in the problem
@param vertices The list of vertices
@param v The number of vertices so far
@param functions The list of functions
@param f The number of functions so far

*/
void AddOcculted(double r, int no, double x0[no], double y0[no], double ro[no], int maxvertices, int maxfunctions, double *vertices, int *v, FUNCTION *functions, int *f) {
  /*
  
  */ 

  int i;
  double ro2[no];
  for (i = 0; i < no; i++) {
    ro2[i] = ro[i] * ro[i];
    ro2[i] += TINY * ro2[i];
  }
  double d, A, x, y, frac, cost, sint;
  ELLIPSE *occulted;
  occulted = malloc(sizeof(ELLIPSE));
  
  // Initialize the occulted planet at the origin
  occulted->r = r;
  occulted->x0 = 0.;
  occulted->y0 = 0.;
  occulted->circle = 1;
  
  // Occulted planet x minimum
  occulted->xmin = -r;
  // Is this point inside __at least__ one occultor?
  for (i = 0; i < no; i++) {
    if ((occulted->xmin - x0[i]) * (occulted->xmin - x0[i]) + y0[i] * y0[i] < ro2[i]) {
      vertices[(*v)++] = occulted->xmin;
      if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
      break;    
    }
  }  
  
  // Occulted planet x maximum
  occulted->xmax = r;
  // Is this point inside __at least__ one occultor?
  for (i = 0; i < no; i++) {
    if ((occulted->xmax - x0[i]) * (occulted->xmax - x0[i]) + y0[i] * y0[i] < ro2[i]) {
      vertices[(*v)++] = occulted->xmax;
      if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
      break;    
    }
  }  

  // Vertices of intersection with the occultors
  // Adapted from http://mathworld.wolfram.com/Circle-CircleIntersection.html
  for (i = 0; i < no; i++) {
    d = sqrt(x0[i] * x0[i] + y0[i] * y0[i]);
    if (d < (ro[i] + r)) {
      A = (-d + r - ro[i]) * (-d - r + ro[i]) * (-d + r + ro[i]) * (d + r + ro[i]);
      if (A >= 0) {
        y = sqrt(A) / (2 * d);
        x = sqrt(r * r - y * y);
        frac = y0[i] / x0[i];
        cost = 1 / sqrt(frac * frac + 1);
        sint = frac * cost;
        if (x0[i] < 0) {
          cost *= -1;
          sint *= -1;
        }
        vertices[(*v)++] = x * cost + y * sint;
        vertices[(*v)++] = x * cost - y * sint;        
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
      }
    }
  }
  
  // Finally, the boundary curves and their integrals
  functions[*f].ellipse = occulted;
  functions[*f].curve = fupper;
  functions[(*f)++].integral = iupper;
  functions[*f].ellipse = occulted;
  functions[*f].curve = flower;
  functions[(*f)++].integral = ilower;
  if (*f > maxfunctions) {
    printf("ERROR: Maximum number of functions exceeded.\n");
    abort();
  }
  
}

/**
Computes the flux of an eyeball planet occulted by one or more bodies
over a grid of wavelengths.

@param r The radius of the occulted planet
@param no The number of occulting bodies
@param x0 The x coordinate of the center of each occulting body
@param y0 The y coordinate of the center of each occulting body
@param ro The radius of each occulting body
@param theta The eyeball phase angle
@param tnight The planet's night side temperature in K (airless limit only)
@param teff The planet's effective temperature in K (blackbody limit only)
@param distance The distance to the system in parsecs
@param mintheta The minimum absolute value of the phase angle, 
       below which it is assumed constant to prevent numerical errors
@param maxvertices Maximum number of vertices in the problem
@param maxfunctions Maximum number of functions in the problem
@param adaptive Adaptive zenith angle grid? Limb-darkened limit only
@param circleopt Treat ellipses as circles for limb-darkened bodies? There's
       no reason *not* to do this! Option left in mainly for testing purposes.
@param batmanopt Use the BATMAN algorithm to speed up limb-darkened occultations?
@param quarticsolver Which quartic solver to use?
@param nu The number of limb darkening coefficients
@param nz The number of zenith angle slices
@param nw The size of the wavelength grid
@param u The limb darkening coefficient grid: \a nu coefficients at each wavelength
@param lambda The wavelength grid in m
@param flux The wavelength-dependent light curve computed by this function
@param maptype The code for the radiance map
@param radiancemap A pointer to the RADIANCEMAP function
@param quiet Suppress output?
@param iErr Flag set if an error occurs

*/
void OccultedFlux(double r, int no, double x0[no], double y0[no], double ro[no], double theta, double tnight, double teff, double distance, double mintheta, int maxvertices, int maxfunctions, int adaptive, int circleopt, int batmanopt, int quarticsolver, int nu, int nz, int nw,  double u[nu * nw], double lambda[nw], double flux[nw], int maptype, RADIANCEMAP radiancemap, int quiet, int *iErr) {

  int i, j, k, m;
  int oob;
  int b = 0;
  int v = 0;
  int f = 0;
  int good;
  double lmin, lmax, za;
  double xL, xR, x, y, area;
  double r2 = r * r;
  r2 += SMALL * r2;
  double d, dmin, dmax;
  double ro2[no];
  for (i = 0; i < no; i++) {
    ro2[i] = ro[i] * ro[i];
    ro2[i] += SMALL * ro2[i];
  }
  double norm;
  double B0[nw];
  double d2 = distance * distance * PARSEC * PARSEC;
  *iErr = ERR_NONE;
  
  // Zero out the flux
  for (m = 0; m < nw; m++) 
    flux[m] = 0.;
  
  // Can we optimize this with the *batman* algorithm?
  // TODO: Adapt this to work with custom radial maps?
  if ((no == 1) && (batmanopt) && (maptype == MAP_RADIAL_DEFAULT)) {
    BatmanFlux(r, x0[0], y0[0], ro[0], teff, distance, nu, nz, nw, u, lambda, flux);
    return;
  }
    
  // Handle numerical issues
  if (theta < MINCRESCENT) {
  
    // Numerical issues pop up when the crescent is too thin and
    // the nightside is dark, so let's increase the number of grid slices
    if (nz < CRESCENTNZ)
      nz = CRESCENTNZ;
      
  } else if (fabs(theta) < mintheta) {
  
    // Avoid the singular point
    theta = mintheta;
    
  }
  
  // Allocate the expensive stuff
  double vertices[maxvertices];
  FUNCTION functions[maxfunctions];
  FUNCTION boundaries[maxfunctions];
  double B[nz * no + 1][nw];
  double zenithgrid[nz * no];
  
  // Generate all the shapes and get their vertices and curves
  AddOccultors(r, no, x0, y0, ro, maxvertices, maxfunctions, vertices, &v, functions, &f);     
  AddOcculted(r, no, x0, y0, ro, maxvertices, maxfunctions, vertices, &v, functions, &f); 
  
  // Compute the zenith_angle grid
  if (adaptive && ((maptype == MAP_RADIAL_DEFAULT) || (maptype == MAP_RADIAL_CUSTOM))) {
  
    // Adaptive zenith_angle grid
    // Loop over each occultor
    for (i = 0; i < no; i++) {
    
      // Distance to occultor center
      d = sqrt(x0[i] * x0[i] + y0[i] * y0[i]);
      
      // Distance to far side
      dmax = (d + ro[i]) / r + SMALL;
      if (dmax >= 1 - SMALL) dmax = 1;
      lmax = asin(dmax);
          
      // Distance to near side
      dmin = (d - ro[i]) / r - SMALL;
      if (dmin <= SMALL) dmin = 0;
      lmin = asin(dmin);

      // Add to zenith_angle grid
      for (j = 0; j < nz; j++) {
        zenithgrid[nz * i + j] = (j + 1.) / (nz + 1.) * (lmax - lmin) + lmin;
      }
      
      // Adjust the boundaries
      zenithgrid[nz * i] = lmin;
      zenithgrid[nz * (i + 1) - 1] = lmax;
      
    }

    // Sort the grid
    qsort(zenithgrid, nz * no, sizeof(double), dblcomp);
    
  } else {
  
    // Linearly spaced zenith_angle grid
    lmin = 0 - SMALL;
    lmax = PI + SMALL;
    for (i = 0; i < nz * no; i++) {
      zenithgrid[i] = (i + 1.) / (nz * no + 1.) * (lmax - lmin) + lmin;
    }
    
  }

  // Add the ellipses  
  for (i = 0; i < nz * no; i++) {

    if ((circleopt || (quarticsolver == QGSL)) && fabs(cos(theta)) < SMALL)
      AddZenithAngleCircle(zenithgrid[i], r, no, x0, y0, ro, maxvertices, maxfunctions, vertices, &v, functions, &f);
    else
      AddZenithAngleEllipse(zenithgrid[i], r, no, x0, y0, ro, theta, quarticsolver, maxvertices, maxfunctions, vertices, &v, functions, &f);
  
  }
   
  // Pre-compute the blackbody normalization if necessary
  if (maptype == MAP_RADIAL_DEFAULT) {
    for (j = 0; j < nw; j++) {
      norm = 0;
      for (k = 0; k < nu; k++) {
        norm += u[nw * k + j] / ((k + 2) * (k + 3));
      }
      norm = 1 - 2 * norm;
      B0[j] = Blackbody(lambda[j], teff) / norm;
    }
  }
  
  // Pre-compute the surface intensity in each zenith_angle slice
  SurfaceIntensity(tnight, teff, maptype, radiancemap, nz * no, zenithgrid, nw, lambda, nu, u, B0, B);
  
  // Sort the vertices
  qsort(vertices, v, sizeof(double), dblcomp);
      
  // Loop over all vertices
  for (i = 0; i < v - 1; i++) {

    // Get the integration limits, with a
    // small perturbation for stability
    xL = vertices[i] + TINY;
    xR = vertices[i+1] - TINY;
    
    // Check if they are identical
    if (xR <= xL + TINY)
      continue;
        
    // Bisect the limits. Find the boundary functions that are
    // finite valued at this point and sort their integrals 
    // in order of increasing function value.
    x = 0.5 * (xL + xR);
    
    // Loop over all functions
    b = 0;
    for (j = 0; j < f; j++) {
    
      y = curve(x, functions[j]);
      
      // Check that it's inside the planet
      if (!(isnan(y)) && (x * x + y * y < r2)) {
        
        // Check that it's inside __at least__ one occultor
        for (k = 0; k < no; k++) {
          if ((x - x0[k]) * (x - x0[k]) + (y - y0[k]) * (y - y0[k]) < ro2[k]) {
            
            // Add the function to the list of boundaries
            boundaries[b] = functions[j];
            boundaries[b++].y = y;
            break;
          
          }
        
        }
      
      }
      
    }

    // Sort boundaries based on y
    qsort(boundaries, b, sizeof(FUNCTION), funcomp);
    
    // Loop over all regions bounded by xL and xR
    for (j = 0; j < b - 1; j++) {
    
      // Get the midpoint
      y = 0.5 * (boundaries[j + 1].y + boundaries[j].y);
      
      // Check again: is it in the planet?
      if (x * x + y * y > r2) continue;
      
      // Check again: is it in at least one occultor?
      good = 0;
      for (k = 0; k < no; k++) {
        if ((x - x0[k]) * (x - x0[k]) + (y - y0[k]) * (y - y0[k]) < ro2[k]) {
          good = 1;
          break;
        }
      }
      if (!good) continue;
      
      // The area of each region is just the difference of successive integrals
      oob = 0;
      area = integral(xL, xR, boundaries[j + 1], &oob) - integral(xL, xR, boundaries[j], &oob);      
      if (oob) *iErr = ERR_OOB;
      
      // Get the zenith_angle of the midpoint
      za = ZenithAngle(x, y, r, theta);
      
      // Get the index `k` of the zenith_angle grid *above* this zenith_angle.
      // B[k] is the intensity of this region.
      for (k = 0; k < nz * no; k++) {
        if (zenithgrid[k] > za)
          break;
      }
      
      // Multiply by the area of the zenith_angle slice to get the flux at each wavelength
      for (m = 0; m < nw; m++)
        flux[m] += B[k][m] * area;
    } 
    
  }
  
  // Scale flux by REARTH**2, convert from m^-1 to micron^-1, and
  // divide by the square of the distance to get a radiance (W / m^2 / um)
  for (m = 0; m < nw; m++) 
    flux[m] *= REARTH * REARTH * MICRON / d2;
  
  // Free all of the ellipse instances
  for (j = 0; j < f; j+=2) {
    free(functions[j].ellipse);
  }  
   
}

/**
Computes the total flux of an eyeball over a grid of wavelengths. This routine
(hackishly) calls \a OccultedFlux with an imaginary occultor covering the entire disk
of the body.

@param r The radius of the occulted planet
@param theta The eyeball phase angle
@param tnight The planet's night side temperature in K (airless limit only)
@param teff The planet's effective temperature in K (blackbody limit only)
@param distance The distance to the system in parsecs
@param mintheta The minimum absolute value of the phase angle, 
       below which it is assumed constant to prevent numerical errors
@param maxvertices Maximum number of vertices in the problem
@param maxfunctions Maximum number of functions in the problem
@param adaptive Adaptive zenith angle grid? Limb-darkened limit only
@param circleopt Treat ellipses as circles for limb-darkened bodies? There's
       no reason *not* to do this! Option left in mainly for testing purposes.
@param batmanopt Use the BATMAN algorithm to speed up limb-darkened occultations?
@param quarticsolver Which quartic solver to use?
@param nu The number of limb darkening coefficients
@param nz The number of zenith angle slices
@param nw The size of the wavelength grid
@param u The limb darkening coefficient grid: \a nu coefficients at each wavelength
@param lambda The wavelength grid in m
@param flux The wavelength-dependent light curve computed by this function
@param maptype The code for the radiance map
@param radiancemap A pointer to the RADIANCEMAP function
@param quiet Suppress output?
@param iErr Flag set if an error occurs

*/
void UnoccultedFlux(double r, double theta, double tnight, double teff, double distance, double mintheta, int maxvertices, int maxfunctions, int adaptive, int circleopt,  int batmanopt, int quarticsolver, int nu, int nz, int nw, double u[nu * nw], double lambda[nw], double flux[nw], int maptype, RADIANCEMAP radiancemap, int quiet, int *iErr) {

  double x0[1] = {0};
  double y0[1] = {0};
  double ro[1] = {2 * r};
  
  // Hack: compute the occulted flux with a single huge occultor
  OccultedFlux(r, 1, x0, y0, ro, theta, tnight, teff, distance, mintheta, maxvertices, maxfunctions, adaptive, circleopt, batmanopt, quarticsolver, nu, nz, nw, u, lambda, flux, maptype, radiancemap, quiet, iErr);
    
}