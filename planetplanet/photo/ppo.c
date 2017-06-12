#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"

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
  double irrad;
  int no;
  double xo[np-1], yo[np-1], ro[np-1];
  double tmp[nw];
  double theta;
  int t, p, o, w;
  int iErr = ERR_NONE;
  
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
    
  // Compute the stellar flux
  UnoccultedFlux(body[0]->r, PI / 2., 0., 0., 0., body[0]->teff, settings.polyeps1, settings.polyeps2, 
                 settings.maxpolyiter, settings.mintheta, settings.maxvertices, settings.maxfunctions,
                 settings.adaptive, body[0]->nu, body[0]->nl, nw, body[0]->u, 
                 wavelength, tmp, settings.quiet, &iErr);
  for (t = 0; t < nt; t++) {
    for (w = 0; w < nw; w++)
      body[0]->flux[nw * t + w] = tmp[w];
  }
  
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
                       settings.polyeps1, settings.polyeps2, settings.maxpolyiter,
                       settings.mintheta, settings.maxvertices, settings.maxfunctions, 
                       settings.adaptive, body[p]->nu, body[p]->nl, nw, body[p]->u, wavelength, 
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
      no = 0;
      
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
          
          // If the body is in quadrants II or III, we need to mirror
          // the problem, since `OccultedFlux` assumes the star is always
          // to the *left* of the body.
          if (body[p]->x[t] < 0) 
            xo[no] = -dx;
          else 
            xo[no] = dx;
          yo[no] = dy;
          ro[no++] = body[o]->r;

        }
        
      }
      
      // Finally, compute the light curve for this planet
      if (no > 0) {
      
        // The orbital phase (edge-on limit!)
        theta = atan(body[p]->z[t] / fabs(body[p]->x[t]));
        
        // The irradiation on the planets
        if (p > 0) {
          dx = (body[0]->x[t] - body[p]->x[t]);
          dy = (body[0]->y[t] - body[p]->y[t]);
          dz = (body[0]->z[t] - body[p]->z[t]);
          d2 = dx * dx + dy * dy + dz * dz;
          irrad = (body[0]->r * body[0]->r) * SBOLTZ * (body[0]->teff * body[0]->teff * body[0]->teff * body[0]->teff) / d2;
        
          // The planet effective temperature from radiation balance
          if (body[p]->blackbody) {
            body[p]->teff = pow(irrad * (1 - body[p]->albedo) / (4 * SBOLTZ), 0.25);
          }
        
        } else 
          irrad = 0.;
        
        // Call the eyeball routine
        OccultedFlux(body[p]->r, no, xo, yo, ro, theta, body[p]->albedo, 
                     irrad, body[p]->tnight,  body[p]->teff, settings.polyeps1, settings.polyeps2, 
                     settings.maxpolyiter, settings.mintheta, settings.maxvertices,
                     settings.maxfunctions, settings.adaptive, body[p]->nu, body[p]->nl, nw, 
                     body[p]->u, wavelength, tmp, settings.quiet, &iErr);
      
        // Update the body light curve
        for (w = 0; w < nw; w++)
          body[p]->flux[nw * t + w] -= tmp[w];
      
      }
    }
  
  }

  // Log
  if (!settings.quiet) {
    if (iErr == ERR_OOB)
      printf("WARNING: Precision loss detected in integration.\n");
    printf("Done!\n");
  }
  
  return iErr;

}