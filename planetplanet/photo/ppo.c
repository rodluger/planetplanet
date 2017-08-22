/**
@file ppo.c
@brief Main interface to the photodynamical routines.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"

/**
Computes the phase angle \a theta and the rotation angle \a gamma at a given
point in a planet's orbit.

@param x0 The x coordinate of the planet's position on the sky in Earth radii
@param y0 The y coordinate of the planet's position on the sky in Earth radii
@param z0 The z coordinate of the planet's position on the sky in Earth radii
@param rp The radius of the planet in Earth radii
@param Omega The longitude of ascending node in radians
@param Lambda The latitudinal hotspot offset in radians (north positive)
@param Phi The longitudinal hotspot offset in radians (east positive)
@param theta The eyeball phase angle
@param gamma The eyeball rotation angle

*/
void GetAngles(double x0, double y0, double z0, double rp, double Omega, double Lambda, double Phi, double *theta, double *gamma) {
 
  double r = sqrt(x0 * x0 + y0 * y0 + z0 * z0);
  double d, xstar, ystar, zstar;
  
  // The "easy" case, with no hotspot offset
  if ((Lambda == 0) && (Phi == 0)) {
  
    // The coordinates of the sub-stellar point on the sky
    xstar = x0 * (1 - rp / r);
    ystar = y0 * (1 - rp / r);
    zstar = z0 * (1 - rp / r);
  
  // The general case   
  } else {
    
    double x, y, z;
    double xprime, yprime, zprime;
    double rxz;
    double tmpx, tmpy;
    double cosO = cos(Omega);
    double sinO = sin(Omega);
    double cosl = cos(Lambda);
    double sinl = sin(Lambda);
    double cosp = cos(Phi);
    double sinp = sin(Phi);
    
    // The position of the planet in the rotated sky plane
    x = x0 * cosO + y0 * sinO;
    y = y0 * cosO - x0 * sinO;
    z = z0;
        
    // Coordinates of the hotspot in a frame where the planet is
    // at x, y, z = (0, 0, r), at full phase
    xprime = rp * cosl * sinp;
    yprime = rp * sinl;
    zprime = r - rp * cosl * cosp;

    // Transform to the rotated sky plane
    rxz = sqrt(x * x + z * z);
    xstar = ((z * r) * xprime - (x * y) * yprime + (x * rxz) * zprime) / (r * rxz);
    ystar = (rxz * yprime + y * zprime) / r;
    zstar = (-(x * r) * xprime - (y * z) * yprime + (z * rxz) * zprime) / (r * rxz);
        
    // Transform back to the true sky plane
    tmpx = xstar * cosO - ystar * sinO;
    tmpy = ystar * cosO + xstar * sinO;
    xstar = tmpx;
    ystar = tmpy;
                   
  }
  
  // The rotation angle
  *gamma = PI + atan2(ystar - y0, xstar - x0);
  
  // The distance from the center of the planet disk to the substellar point
  d = sqrt((xstar - x0) * (xstar - x0) + (ystar - y0) * (ystar - y0));
  
  // The phase angle
  if (zstar - z0 <= 0)
    *theta = acos(d / rp);
  else
    *theta = -acos(d / rp);
  
}

/**
Raises an integer to an integer power.

@param base The base
@param exp The exponent
@return The result of base^exponent

*/
int ipow(int base, int exp){
  
  int result = 1;
  while (exp) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
  
}

/**
Evolves the orbital positions of all bodies forward in time with either a Keplerian
or an N-body solver.

@param nt The size of the time array
@param time The array of times at which to evaluate the orbits
@param np The number of particles (bodies) in the system
@param body An array of BODY pointers corresponding to all the bodies in the system
@param settings An instance of the SETTINGS class containing all settings

*/
int Orbits(int nt, double time[nt], int np, BODY **body, SETTINGS settings){

  int t, p, o;
  double dx, dy, d;
  int iErr = ERR_NONE;
  
  // Initialize the arrays for each body
  for (p = 0; p < np; p++) {
    body[p]->nt = nt;
    for (t = 0; t < nt; t++)
      body[p]->time[t] = time[t];
  }
  
  // Solve for the orbits
  if (settings.nbody)
    iErr = NBody(np, body, settings);
  else
    iErr = Kepler(np, body, settings);
  if (iErr != ERR_NONE) return iErr;
  
  // Log
  if (!settings.quiet)
    printf("Computing occultation events...\n");
  
  // Loop over the time array
  for (t = 0; t < nt; t++) {
    
    // Loop over each body
    for (p = 0; p < np; p++) {
      
      // Default is no occultation
      body[p]->occultor[t] = 0;
      
      // Loop over all possible occultors, including the star
      // and check whether an occultation is occurring
      for (o = 0; o < np; o++) {
      
        // Skip self
        if (o == p) continue;
      
        // Compute the body-body separation between body `o` and body `p`
        dx = (body[o]->x[t] - body[p]->x[t]);
        dy = (body[o]->y[t] - body[p]->y[t]);
        d = sqrt(dx * dx + dy * dy);
        
        // Is body `o` occulting body `p`?
        if ((d <= body[p]->r + body[o]->r) && (body[p]->z[t] > body[o]->z[t])){

          // Yes! Add to the occultor flag
          body[p]->occultor[t] += ipow(2, o);

        }
        
      }
  
    }
  
  }
  
  // Log
  if (!settings.quiet)
    printf("Done!\n");

  return iErr;

}

/**
Computes the full light curve for all bodies in the system.

@param nt The size of the time array
@param time The array of times at which to evaluate the orbits
@param nw The size of the wavelength grid
@param wavelength The wavelength grid in microns
@param np The number of particles (bodies) in the system
@param body An array of BODY pointers corresponding to all the bodies in the system
@param settings An instance of the SETTINGS class containing all settings

*/
int Flux(int nt, double time[nt], int nw, double wavelength[nw], int np, BODY **body, SETTINGS settings){

  double d, dx, dy, dz, d2;
  double lum, irrad;
  int no;
  double xo[np-1], yo[np-1], ro[np-1];
  double tmp[nw];
  double tmpflux[nw];
  double theta, gamma;
  double tmpx, tmpy;
  int t, p, o, w;
  int iErr = ERR_NONE;
  double norm, sflx;
    
  // Initialize the arrays for each body
  for (p = 0; p < np; p++) {
    body[p]->nw = nw;
    for (w = 0; w < nw; w++)
      body[p]->wavelength[w] = wavelength[w];
    body[p]->nt = nt;
    for (t = 0; t < nt; t++) {
      body[p]->time[t] = time[t];
    }    
  }
    
  // Solve for the orbits
  if (settings.nbody)
    iErr = NBody(np, body, settings);
  else
    iErr = Kepler(np, body, settings);
  if (iErr != ERR_NONE) {
    printf("ERROR: Kepler solver failure (%d).\n", iErr);
    abort();
  }
  
  // Compute the stellar radiance from Planck's law
  norm = PI * body[0]->r * REARTH * body[0]->r * REARTH * MICRON / (settings.distance * settings.distance * PARSEC * PARSEC);
  for (w = 0; w < nw; w++) {
    sflx = Blackbody(wavelength[w], body[0]->teff) * norm;
    body[0]->total_flux[w] = sflx;
    for (t = 0; t < nt; t++) {
      body[0]->flux[nw * t + w] = sflx;
    }
  }
    
  // Pre-compute the stellar luminosity per solid angle
  lum = (body[0]->r * body[0]->r) * SBOLTZ * (body[0]->teff * body[0]->teff * body[0]->teff * body[0]->teff);
  
  // Compute the total flux from each of the planets at full phase
  // and store it in the `total_flux` attribute
  for (p = 1; p < np; p++) {
    
    // The planet effective temperature from radiation balance
    dx = (body[0]->x[0] - body[p]->x[0]);
    dy = (body[0]->y[0] - body[p]->y[0]);
    dz = (body[0]->z[0] - body[p]->z[0]);
    d2 = dx * dx + dy * dy + dz * dz;
    irrad = lum / d2;
    body[p]->teff = pow(irrad * (1 - body[p]->albedo) / (4 * SBOLTZ), 0.25);
    
    // Call the eyeball routine at *full phase*
    UnoccultedFlux(body[p]->r, PI / 2, body[p]->tnight, body[p]->teff,
                   settings.distance, settings.mintheta, settings.maxvertices, settings.maxfunctions, 
                   settings.adaptive, settings.circleopt, settings.batmanopt, settings.quarticsolver, 
                   body[p]->nu, body[p]->nz, nw, body[p]->u, wavelength, 
                   tmp, body[p]->maptype, body[p]->radiancemap, settings.quiet, &iErr);
    for (w = 0; w < nw; w++) {
      body[p]->total_flux[w] = tmp[w];  
              
    }    
  }
  
  // Log
  if (!settings.quiet)
    printf("Computing occultation light curves...\n");

  // Loop over the time array
  for (t = 0; t < nt; t++) {
                
    // Compute the light curve for each body
    for (p = 0; p < np; p++) {
      
      // Compute the phase curve for this body?
      if ((p > 0) && (body[p]->phasecurve) && ((body[p]->maptype == MAP_ELLIPTICAL) || (body[p]->maptype == MAP_ELLIPTICAL_CUSTOM))) {
        
        // TODO: Interpolate here to save time!
        
        // Get the eyeball angles `theta` and `gamma`
        // Note that `gamma` does not matter for phase curves.
        GetAngles(body[p]->x[t], body[p]->y[t], body[p]->z[t], body[p]->r, body[p]->Omega, body[p]->Lambda, body[p]->Phi, &theta, &gamma);
        
        // The planet effective temperature from radiation balance
        dx = (body[0]->x[t] - body[p]->x[t]);
        dy = (body[0]->y[t] - body[p]->y[t]);
        dz = (body[0]->z[t] - body[p]->z[t]);
        d2 = dx * dx + dy * dy + dz * dz;
        irrad = (body[0]->r * body[0]->r) * SBOLTZ * (body[0]->teff * body[0]->teff * body[0]->teff * body[0]->teff) / d2;
        body[p]->teff = pow(irrad * (1 - body[p]->albedo) / (4 * SBOLTZ), 0.25);
                
        // Call the eyeball routine
        UnoccultedFlux(body[p]->r, theta, body[p]->tnight, body[p]->teff,
                       settings.distance, settings.mintheta, settings.maxvertices, settings.maxfunctions, 
                       settings.adaptive, settings.circleopt, settings.batmanopt, settings.quarticsolver, 
                       body[p]->nu, body[p]->nz, nw, body[p]->u, wavelength, 
                       tmp, body[p]->maptype, body[p]->radiancemap, settings.quiet, &iErr);
        for (w = 0; w < nw; w++) {
          body[p]->flux[nw * t + w] = tmp[w];          
        }
      
      } else if ((body[p]->maptype == MAP_RADIAL) || (body[p]->maptype == MAP_RADIAL_CUSTOM)) {
        
        // The flux is constant because the planet emission is radially
        // symmetrical about the center of the planet disk
        for (w = 0; w < nw; w++) {
          body[p]->flux[nw * t + w] = body[p]->total_flux[w];
        }
        
      } else if (p > 0) {
        
        // Initialize to zero at all wavelengths
        for (w = 0; w < nw; w++) {
          body[p]->flux[nw * t + w] = 0;
        }
        
      }
            
      // Default is no occultation
      no = 0;
      body[p]->occultor[t] = 0;
      for (w = 0; w < nw; w++)
        tmpflux[w] = 0;
      
      // Loop over all possible occultors, including the star
      // and check whether an occultation is occurring
      for (o = 0; o < np; o++) {
  
        // Skip self
        if (o == p) continue;
  
        // Compute the body-body separation between body `o` and body `p`
        dx = body[o]->x[t] - body[p]->x[t];
        dy = body[o]->y[t] - body[p]->y[t];
        dz = body[o]->z[t] - body[p]->z[t];
        d = sqrt(dx * dx + dy * dy);
        
        // Is body `o` occulting body `p`?
        if ((d <= body[p]->r + body[o]->r) && (dz < 0)){

          // Yes! Add to the occultor flag if it's not already there
          if (!(body[p]->occultor[t] & ipow(2, o)))
            body[p]->occultor[t] += ipow(2, o);
      
          // Relative position of the occultor
          xo[no] = dx;          
          yo[no] = dy;
          ro[no++] = body[o]->r;

        }
    
      }
              
      // Now compute the light curve for this planet
      if (no > 0) {
        
        if ((p > 0) && ((body[p]->maptype == MAP_ELLIPTICAL) || (body[p]->maptype == MAP_ELLIPTICAL_CUSTOM))) {
        
          // Get the eyeball angles `theta` and `gamma`
          GetAngles(body[p]->x[t], body[p]->y[t], body[p]->z[t], body[p]->r, body[p]->Omega, body[p]->Lambda, body[p]->Phi, &theta, &gamma);
        
          // Rotate the occultors to a frame in which the ellipses are symmetric about the x axis
          for (o = 0; o < no; o++) {
            tmpx = xo[o] * cos(gamma) + yo[o] * sin(gamma);
            tmpy = yo[o] * cos(gamma) - xo[o] * sin(gamma);
            xo[o] = tmpx;
            yo[o] = tmpy;          
          }
        
        } else {
          
          // The body is radially symmetric
          theta = PI / 2.;
          
        }
        
        // The planet effective temperature from radiation balance
        if (p > 0) {
          dx = (body[0]->x[t] - body[p]->x[t]);
          dy = (body[0]->y[t] - body[p]->y[t]);
          dz = (body[0]->z[t] - body[p]->z[t]);
          d2 = dx * dx + dy * dy + dz * dz;
          irrad = lum / d2;
          body[p]->teff = pow(irrad * (1 - body[p]->albedo) / (4 * SBOLTZ), 0.25);
        }
        
        // Call the eyeball routine
        OccultedFlux(body[p]->r, no, xo, yo, ro, theta, body[p]->tnight,  body[p]->teff, settings.distance, 
                     settings.mintheta, settings.maxvertices, settings.maxfunctions, settings.adaptive, 
                     settings.circleopt, settings.batmanopt, settings.quarticsolver, 
                     body[p]->nu, body[p]->nz, nw, body[p]->u, wavelength, tmp, body[p]->maptype, body[p]->radiancemap, settings.quiet, &iErr);
        
        // Update the body light curve
        for (w = 0; w < nw; w++)
          body[p]->flux[nw * t + w] -= tmp[w];
                    
      }
      
    }    
  
  }

  // Log
  if (!settings.quiet) {
    if (iErr == ERR_OOB) {
      printf("WARNING: Precision loss detected in integration.\n");
      printf("Consider increasing `nz` or increasing `mintheta`.\n");
    }
    printf("Done!\n");
  }
  
  return iErr;

}