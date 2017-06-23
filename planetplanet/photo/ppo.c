#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"

double QuadInterp(int n, double x[n], int t, double i) {
  /*
  
  */
  
  double xL, xC, xR;
  
  if ((t >= n) || (t < 0)) {
    printf("ERROR: Interpolation out of bounds.\n");
    abort();
  }
  
  // Left
  if (t == 0)  xL = x[t] - 0.5 * (x[t+1] - x[t]);
  else xL = x[t-1];
  
  // Right
  if (t == n - 1) xR = x[t] + 0.5 * (x[t] - x[t-1]);
  else xR = x[t+1];

  // Center
  xC = x[t];

  return (2 * (xL + xR) - 4 * xC) * i * i + (4 * xC - 3 * xL - xR) * i + xL;

}

int ipow(int base, int exp){
  /*
  
  */
  
  int result = 1;
  while (exp) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
  
}

int Orbits(int nt, double time[nt], int np, BODY **body, SETTINGS settings){
  /*
  
  Compute just the orbits and whether occultations occur.
  
  */

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

int Flux(int nt, double time[nt], int nw, double wavelength[nw], int np, BODY **body, SETTINGS settings){
  /*
  
  */

  double d, dx, dy, dz, d2;
  double xp, yp, zp;
  double lum, irrad;
  int no;
  double xo[np-1], yo[np-1], ro[np-1];
  double tmp[nw];
  double tmpflux[nw];
  double theta;
  int i, t, p, o, w;
  int iErr = ERR_NONE;
  double norm, sflx;
  
  // Initialize the arrays for each body
  for (p = 0; p < np; p++) {
    body[p]->nw = nw;
    for (w = 0; w < nw; w++)
      body[p]->wavelength[w] = wavelength[w];
    body[p]->nt = nt;
    for (t = 0; t < nt; t++)
      body[p]->time[t] = time[t];
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
    for (t = 0; t < nt; t++)
      body[0]->flux[nw * t + w] = sflx;
  }
    
  // Pre-compute the stellar luminosity per solid angle
  lum = (body[0]->r * body[0]->r) * SBOLTZ * (body[0]->teff * body[0]->teff * body[0]->teff * body[0]->teff);
  
  // Log
  if (!settings.quiet)
    printf("Computing occultation light curves...\n");
  
  // Loop over the time array
  for (t = 0; t < nt; t++) {
                
    // Compute the light curve for each body
    for (p = 0; p < np; p++) {
      
      // Compute the phase curve for this body?
      if ((p > 0) && (body[p]->phasecurve)) {
        
        // The orbital phase (edge-on limit!)
        theta = atan(body[p]->z[t] / fabs(body[p]->x[t]));
        
        // The irradiation
        dx = (body[0]->x[t] - body[p]->x[t]);
        dy = (body[0]->y[t] - body[p]->y[t]);
        dz = (body[0]->z[t] - body[p]->z[t]);
        d2 = dx * dx + dy * dy + dz * dz;
        irrad = (body[0]->r * body[0]->r) * SBOLTZ * (body[0]->teff * body[0]->teff * body[0]->teff * body[0]->teff) / d2;
        
        // The planet effective temperature from radiation balance
        if (body[p]->blackbody) {
          body[p]->teff = pow(irrad * (1 - body[p]->albedo) / (4 * SBOLTZ), 0.25);
        }
        
        // Call the eyeball routine
        UnoccultedFlux(body[p]->r, theta, body[p]->albedo, irrad, body[p]->tnight, body[p]->teff,
                       settings.distance, settings.polyeps1, settings.polyeps2, settings.maxpolyiter,
                       settings.mintheta, settings.maxvertices, settings.maxfunctions, 
                       settings.adaptive, body[p]->nu, body[p]->nz, nw, body[p]->u, wavelength, 
                       tmp, settings.quiet, &iErr);
        for (w = 0; w < nw; w++)
          body[p]->flux[nw * t + w] = tmp[w];
      
      } else if (p > 0) {
        
        // Initialize to zero at all wavelengths
        for (w = 0; w < nw; w++) {
          body[p]->flux[nw * t + w] = 0;
        }
        
      }
            
      // Default is no occultation
      body[p]->occultor[t] = 0;
      for (w = 0; w < nw; w++)
        tmpflux[w] = 0;

      // Oversample the light curve
      for (i = 0; i < settings.oversample; i++) { 
          
        // Reset the number of occultors
        no = 0;
        
        // Quadratically interpolate to get the xyz coordinates of the occulted body
        // at this oversampled time grid point
        if (settings.oversample > 1) {
          xp = QuadInterp(nt, body[p]->x, t, ((float) i) / (settings.oversample - 1.));
          yp = QuadInterp(nt, body[p]->y, t, ((float) i) / (settings.oversample - 1.));
          zp = QuadInterp(nt, body[p]->z, t, ((float) i) / (settings.oversample - 1.));
        } else {
          xp = body[p]->x[t];
          yp = body[p]->y[t];
          zp = body[p]->z[t];
        }
                
        // Loop over all possible occultors, including the star
        // and check whether an occultation is occurring
        for (o = 0; o < np; o++) {
    
          // Skip self
          if (o == p) continue;
    
          // Compute the body-body separation between body `o` and body `p` by quadratic interpolation
          if (settings.oversample > 1) {
            dx = QuadInterp(nt, body[o]->x, t, ((float) i) / (settings.oversample - 1.)) - xp;
            dy = QuadInterp(nt, body[o]->y, t, ((float) i) / (settings.oversample - 1.)) - yp;
            dz = QuadInterp(nt, body[o]->z, t, ((float) i) / (settings.oversample - 1.)) - zp;
          } else {
            dx = body[o]->x[t] - xp;
            dy = body[o]->y[t] - yp;
            dz = body[o]->z[t] - zp;
          }
          d = sqrt(dx * dx + dy * dy);
          
          // Is body `o` occulting body `p`?
          if ((d <= body[p]->r + body[o]->r) && (dz < 0)){

            // Yes! Add to the occultor flag if it's not already there
            if (!(body[p]->occultor[t] & ipow(2, o)))
              body[p]->occultor[t] += ipow(2, o);
        
            // If the body is in quadrants II or III, we need to mirror
            // the problem, since `OccultedFlux` assumes the star is always
            // to the *left* of the body.
            if (xp < 0) 
              xo[no] = -dx;
            else 
              xo[no] = dx;
            yo[no] = dy;
            ro[no++] = body[o]->r;

          }
      
        }
                
        // Now compute the light curve for this planet
        if (no > 0) {

          // The orbital phase (edge-on limit!)
          theta = atan(zp / fabs(xp));
    
          // The irradiation on the planets; assume the star
          // hasn't moved across this oversampled time step
          if (p > 0) {
            dx = (body[0]->x[t] - xp);
            dy = (body[0]->y[t] - yp);
            dz = (body[0]->z[t] - zp);
            d2 = dx * dx + dy * dy + dz * dz;
            irrad = lum / d2;
            
            // The planet effective temperature from radiation balance
            if (body[p]->blackbody) {
              body[p]->teff = pow(irrad * (1 - body[p]->albedo) / (4 * SBOLTZ), 0.25);
            }
    
          } else {
            irrad = 0.;
          }
      
          // Call the eyeball routine
          OccultedFlux(body[p]->r, no, xo, yo, ro, theta, body[p]->albedo, 
                       irrad, body[p]->tnight,  body[p]->teff, settings.distance, settings.polyeps1, settings.polyeps2, 
                       settings.maxpolyiter, settings.mintheta, settings.maxvertices,
                       settings.maxfunctions, settings.adaptive, body[p]->nu, body[p]->nz, nw, 
                       body[p]->u, wavelength, tmp, settings.quiet, &iErr);
              
          // Update the body light curve
          for (w = 0; w < nw; w++)
            tmpflux[w] -= tmp[w];
    
        }
      }
      
      // Finally, take the average of the oversampled flux
      for (w = 0; w < nw; w++)
        body[p]->flux[nw * t + w] += tmpflux[w] / settings.oversample;
      
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