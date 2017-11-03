/**
@file ppo.c
@brief Main interface to the photodynamical routines.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"

/**
The Rodrigues rotation formula in 3D, given a vector \a v, a unit vector \a k
normal to the plane of rotation, and the angle of rotation \a theta in radians.
The rotation is applied in-place to the vector \a v.

@param vx The x coordinate of the vector \a v
@param vy The y coordinate of the vector \a v
@param vz The z coordinate of the vector \a v
@param kx The x coordinate of the unit vector \a k
@param ky The y coordinate of the unit vector \a k
@param kz The z coordinate of the unit vector \a k
@param theta The angle of rotation in radians

*/
void Rodrigues(double *vx, double *vy, double *vz, double kx, double ky, double kz, double theta) {
  
  double ct = cos(theta);
  double st = sin(theta);
  double rx, ry, rz;
  double kdotv = (kx * *vx) + (ky * *vy) + (kz * *vz);
  
  rx = *vx * ct + (ky * *vz - kz * *vy) * st + kx * kdotv * (1 - ct);
  ry = *vy * ct + (kz * *vx - kx * *vz) * st + ky * kdotv * (1 - ct);
  rz = *vz * ct + (kx * *vy - ky * *vx) * st + kz * kdotv * (1 - ct);
  
  *vx = rx;
  *vy = ry;
  *vz = rz;
  
}

/**
Computes the phase angle \a theta and the rotation angle \a gamma at a given
point in a planet's orbit.

@param x The x coordinate of the planet's position on the sky 
@param y The y coordinate of the planet's position on the sky 
@param z The z coordinate of the planet's position on the sky 
@param vx The x coordinate of the planet's velocity on the sky 
@param vy The y coordinate of the planet's velocity on the sky 
@param vz The z coordinate of the planet's velocity on the sky
@param Lambda The latitudinal hotspot offset in radians (north positive)
@param Phi The longitudinal hotspot offset in radians (east positive)
@param theta The eyeball phase angle
@param gamma The eyeball rotation angle

*/
void GetAngles(double x, double y, double z, double vx, double vy, double vz, double Lambda, double Phi, double *theta, double *gamma) {
 
  double r, xstar, ystar, zstar, v, d;
  
  // The orbital radius
  r = sqrt(x * x + y * y + z * z);
  
  // The orbital speed
  v = sqrt(vx * vx + vy * vy + vz * vz);
  
  // Normalize the vectors
  x /= r;
  y /= r;
  z /= r;
  vx /= v;
  vy /= v;
  vz /= v;
  
  // The position vector of the substellar point,
  // relative to the planet center, normalized to 1
  xstar = -x;
  ystar = -y;
  zstar = -z;
  
  // Do we need to apply a hotspot offset?
  if ((Lambda != 0) || (Phi != 0)) {
  
    double vlon, vlonx, vlony, vlonz;
    double vlat, vlatx, vlaty, vlatz;

    // Unit vector normal to the longitudinal plane
    vlonx = (y * vz - z * vy);
    vlony = (z * vx - x * vz);
    vlonz = (x * vy - y * vx);
    vlon = sqrt(vlonx * vlonx + vlony * vlony + vlonz * vlonz);
    vlonx /= vlon;
    vlony /= vlon;
    vlonz /= vlon;    
    
    // Apply the longitudinal offset
    Rodrigues(&xstar, &ystar, &zstar, vlonx, vlony, vlonz, Lambda);

    // Unit vector normal to the latitudinal plane
    vlatx = (vlony * zstar - vlonz * ystar);
    vlaty = (vlonz * xstar - vlonx * zstar);
    vlatz = (vlonx * ystar - vlony * xstar);
    vlat = sqrt(vlatx * vlatx + vlaty * vlaty + vlatz * vlatz);
    vlatx /= vlat;
    vlaty /= vlat;
    vlatz /= vlat;
  
    // Apply the latitudinal offset
    Rodrigues(&xstar, &ystar, &zstar, vlatx, vlaty, vlatz, Phi);

  }
  
  // Projected distance from planet center to hotspot
  d = sqrt(xstar * xstar + ystar * ystar);
  
  // Get the rotation and phase angles
  *gamma = atan2(ystar, xstar) + PI;
  if (zstar <= 0)
    *theta = acos(d);
  else
    *theta = -acos(d);
  
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
Computes the time(s) of the next occultation(s) of the body of index
`occulted` by the body(ies) of index(ices) `occultors`.

@param nt The size of the time array
@param time The array of times at which to evaluate the orbits
@param np The number of particles (bodies) in the system
@param body An array of BODY pointers corresponding to all the bodies in the system
@param settings An instance of the SETTINGS class containing all settings
@param occulted The index of the body being occulted
@param noccultors Size of occultor array
@param occultors The indices of the occultors of `occulted`
@param noccultations The number of occultations to look for
@param occultation_times The times of each of the occultations, computed by this function
@param occultation_inds The indices of the occultors corresponding to each of the occultations identified by this function
@param occultation_durs The duration of each of the occultations identified by this function
@return The error code

*/
int NextOccultation(int nt, double time[nt], int np, BODY **body, SETTINGS settings, int occulted, int noccultors, int occultors[noccultors], int noccultations, double occultation_times[noccultations], int occultation_inds[noccultations], double occultation_durs[noccultations]){
    
  int t, p;
  
  // Initialize the arrays for each body
  for (p = 0; p < np; p++) {
    body[p]->nt = nt;
    for (t = 0; t < nt; t++)
      body[p]->time[t] = time[t];
  }
  
  // Find the occultations
  return NBody(np, body, settings, 1, occulted, noccultors, occultors, noccultations, occultation_times, occultation_inds, occultation_durs);

}

/**
Evolves the orbital positions of all bodies forward in time with either a Keplerian
or an N-body solver.

@param nt The size of the time array
@param time The array of times at which to evaluate the orbits
@param np The number of particles (bodies) in the system
@param body An array of BODY pointers corresponding to all the bodies in the system
@param settings An instance of the SETTINGS class containing all settings
@return The error code

*/
int Orbits(int nt, double time[nt], int np, BODY **body, SETTINGS settings){

  int t, p, o;
  double dx, dy, d;
  int iErr = ERR_NONE;
  int dummyInt[0];
  double dummyDouble[0];
  
  // Initialize the arrays for each body
  for (p = 0; p < np; p++) {
    body[p]->nt = nt;
    for (t = 0; t < nt; t++)
      body[p]->time[t] = time[t];
  }
  
  // Solve for the orbits
  if (settings.nbody)
    iErr = NBody(np, body, settings, 0, 0, 0, dummyInt, 0, dummyDouble, dummyInt, dummyDouble);
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
@param continuum The continuum flux (i.e., the total flux without occultations) of the system on a time/wavelength grid
@param np The number of particles (bodies) in the system
@param body An array of BODY pointers corresponding to all the bodies in the system
@param settings An instance of the SETTINGS class containing all settings
@return The error code

*/
int Flux(int nt, double time[nt], int nw, double wavelength[nw], double continuum[nt * nw], int np, BODY **body, SETTINGS settings){

  double d, dx, dy, dz, d2;
  double lum[settings.nstars];
  double irrad;
  int no;
  double xo[np-1], yo[np-1], ro[np-1];
  double tmp[nw];
  double tmpflux[nw];
  double theta, gamma;
  double tmpx, tmpy;
  int t, p, q, o, w;
  int iErr = ERR_NONE;
  int dummyInt[0];
  double dummyDouble[0];
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
  
  // Initialize the continuum flux
  for (t = 0; t < nt; t++) {
    for (w = 0; w < nw; w++) {
      continuum[nw * t + w] = 0;
    }
  }
  
  // Solve for the orbits
  if (settings.nbody)
    iErr = NBody(np, body, settings, 0, 0, 0, dummyInt, 0, dummyDouble, dummyInt, dummyDouble);
  else
    iErr = Kepler(np, body, settings);
  if (iErr != ERR_NONE) {
    printf("ERROR: Kepler solver failure (%d).\n", iErr);
    abort();
  }
  
  // Compute the stellar radiance(s) from Planck's law
  for (p = 0; p < settings.nstars; p++) {
    norm = PI * body[p]->r * REARTH * body[p]->r * REARTH * MICRON / (settings.distance * settings.distance * PARSEC * PARSEC);
    for (w = 0; w < nw; w++) {
      sflx = Blackbody(wavelength[w], body[p]->teff) * norm;
      body[p]->total_flux[w] = sflx;
      for (t = 0; t < nt; t++) {
        body[p]->flux[nw * t + w] = sflx;
      }
    }
  
    // Pre-compute the stellar luminosity per solid angle
    lum[p] = (body[p]->r * body[p]->r) * SBOLTZ * (body[p]->teff * body[p]->teff * body[p]->teff * body[p]->teff);
  
  }

  // Compute the total flux from each of the planets at full phase
  // and store it in the `total_flux` attribute
  for (p = settings.nstars; p < np; p++) {
    
    // The planet effective temperature from radiation balance
    irrad = 0;
    for (q = 0; q < settings.nstars; q++) {
      dx = (body[q]->x[0] - body[p]->x[0]);
      dy = (body[q]->y[0] - body[p]->y[0]);
      dz = (body[q]->z[0] - body[p]->z[0]);
      d2 = dx * dx + dy * dy + dz * dz;
      irrad += lum[q] / d2;
    }
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
      if ((p > settings.nstars - 1) && (body[p]->phasecurve) && ((body[p]->maptype == MAP_ELLIPTICAL_DEFAULT) || (body[p]->maptype == MAP_ELLIPTICAL_CUSTOM))) {
        
        // TODO: Interpolate here to save time!
        
        // Get the eyeball angles `theta` and `gamma`
        // Note that `gamma` does not matter for phase curves.
        GetAngles(body[p]->x[t], body[p]->y[t], body[p]->z[t], body[p]->vx[t], body[p]->vy[t], body[p]->vz[t], body[p]->Lambda, body[p]->Phi, &theta, &gamma);
        
        // The planet effective temperature from radiation balance
        irrad = 0;
        for (q = 0; q < settings.nstars; q++) {
          dx = (body[q]->x[t] - body[p]->x[t]);
          dy = (body[q]->y[t] - body[p]->y[t]);
          dz = (body[q]->z[t] - body[p]->z[t]);
          d2 = dx * dx + dy * dy + dz * dz;
          irrad += lum[q] / d2;
        }
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
      
      } else if ((body[p]->phasecurve) && ((body[p]->maptype == MAP_RADIAL_DEFAULT) || (body[p]->maptype == MAP_RADIAL_CUSTOM))) {
        
        // The flux is constant because the planet emission is radially
        // symmetrical about the center of the planet disk
        for (w = 0; w < nw; w++) {
          body[p]->flux[nw * t + w] = body[p]->total_flux[w];
        }
        
      } else if (p > settings.nstars - 1) {
        
        // NOTE: If phase curves are turned off, I set the body's
        // flux to be the flux at full phase. This is obviously
        // incorrect for an eyeball planet but has a negligible effect on 
        // the total light curve. Recall that occultations are differential 
        // measurements and the continuum is set by the star. If this is
        // ever an issue, just turn phase curves on!
        for (w = 0; w < nw; w++) {
          body[p]->flux[nw * t + w] = body[p]->total_flux[w];
        }
        
      }
        
      // Add to the continuum flux
      for (w = 0; w < nw; w++) {
        continuum[nw * t + w] += body[p]->flux[nw * t + w];
      }
             
      // Default is no occultation
      no = 0;
      body[p]->occultor[t] = 0;
      for (w = 0; w < nw; w++)
        tmpflux[w] = 0;
      
      // Loop over all possible occultors, including the star(s)
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
        
        if ((p > settings.nstars - 1) && ((body[p]->maptype == MAP_ELLIPTICAL_DEFAULT) || (body[p]->maptype == MAP_ELLIPTICAL_CUSTOM))) {
        
          // Get the eyeball angles `theta` and `gamma`
          GetAngles(body[p]->x[t], body[p]->y[t], body[p]->z[t], body[p]->vx[t], body[p]->vy[t], body[p]->vz[t], body[p]->Lambda, body[p]->Phi, &theta, &gamma);

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
        if (p > settings.nstars - 1) {
          // The planet effective temperature from radiation balance
          irrad = 0;
          for (q = 0; q < settings.nstars; q++) {
            dx = (body[q]->x[t] - body[p]->x[t]);
            dy = (body[q]->y[t] - body[p]->y[t]);
            dz = (body[q]->z[t] - body[p]->z[t]);
            d2 = dx * dx + dy * dy + dz * dz;
            irrad += lum[q] / d2;
          }
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