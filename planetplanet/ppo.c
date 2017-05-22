#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"

int Flux(int nt, double time[nt], int nw, double wavelength[nw], int np, PLANET **planet, SETTINGS settings){
  /*
  
  NOTE: This does not yet handle the pathological case of *multiple* simultaneous 
  planet-planet occultations involving the same planet!
  
  */

  double d, dx, dy, x0;
  double tmp[nw];
  double theta;
  int t, p, o, w;
  int iErr = ERR_NONE;
  
  // Initialize the arrays for each planet
  for (p = 0; p < np; p++) {
    for (w = 0; w < nw; w++)
      planet[p]->wavelength[w] = wavelength[w];
    for (t = 0; t < nt; t++)
      planet[p]->time[t] = time[t];
  }

  // Loop over the time array
  for (t = 0; t < nt; t++) {
    
    // Update the time array index in each planet struct
    for (p = 0; p < np; p++)
      planet[p]->t = t;

    // Compute the instantaneous orbital positions of all the planets
    iErr = OrbitXYZ(np, planet, settings);
    if (iErr != ERR_NONE) return iErr;
    
    // Compute the light curve for each planet
    for (p = 0; p < np; p++) {
    
      // Compute the phase curve for this planet?
      if (planet[p]->phasecurve) {
        
        // The orbital phase (edge-on limit!)
        theta = atan(planet[p]->z[t] / fabs(planet[p]->x[t]));
        
        // Call the eyeball routine
        UnoccultedFlux(planet[p]->r, theta, planet[p]->albedo, planet[p]->irrad, 
                       settings.polyeps1, settings.polyeps2, settings.maxpolyiter, 
                       planet[p]->nl, nw, wavelength, tmp);
        for (w = 0; w < nw; w++)
          planet[p]->flux[nw * t + w] = tmp[w];
        
      } else {
        
        // Initialize to zero at all wavelengths
        for (w = 0; w < nw; w++) {
          planet[p]->flux[nw * t + w] = 0;
        }
        
      }

      // Default is no occultation
      planet[p]->occultor[t] = -1;
      
      // Loop over all possible occultors
      for (o = 0; o < np; o++) {
      
        // Skip self
        if (o == p) continue;
      
        // Compute the planet-planet separation between planet `p` and planet `o`
        dx = (planet[p]->x[t] - planet[o]->x[t]);
        dy = (planet[p]->y[t] - planet[o]->y[t]);
        d = sqrt(dx * dx + dy * dy);
        
        // Is planet `o` occulting planet `p`?
        if ((d <= planet[p]->r + planet[o]->r) && (planet[p]->z[t] > planet[o]->z[t])){
          
          // Yes!
          planet[p]->occultor[t] = o;
          
          // If the planet is in quadrants II or III, we need to mirror
          // the problem, since `OccultedFlux` assumes the star is always
          // to the *left* of the planet.
          if (planet[p]->x[t] < 0) 
            x0 = -dx;
          else 
            x0 = dx;
          
          // The orbital phase (edge-on limit!)
          theta = atan(planet[p]->z[t] / fabs(planet[p]->x[t]));
          
          // Call the eyeball routine
          OccultedFlux(planet[p]->r, x0, dy, planet[o]->r, theta, planet[p]->albedo, 
                       planet[p]->irrad, settings.polyeps1, settings.polyeps2, 
                       settings.maxpolyiter, planet[p]->nl, nw, wavelength, tmp);
          
          // Update the planet light curve
          for (w = 0; w < nw; w++)
            planet[p]->flux[nw * t + w] -= tmp[w];
        
        }
        
      }
      
    }
  
  }

  return iErr;

}