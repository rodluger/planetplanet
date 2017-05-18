#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"

void Flux(double time, int n, PLANET planet[n], SETTINGS settings, double flux[n], int occultor[n]){
  /*
  
  NOTE: This does not properly handle the pathological case of *multiple* simultaneous 
  planet-planet occultations involving the same planet!
  
  */

  double d, dx, dy;
  double x0, theta;
  double total;
  int i, j, nlat;
  
  // Compute the instantaneous orbital positions
  // of all the planets
  for (i = 0; i < n; i++)
    OrbitXYZ(time, &planet[i], settings);
  
  // Loop over all planets
  for (i = 0; i < n; i++) {
    
    // Compute the phase curve?
    if ((settings.phasecurve) && (planet[i].noon != planet[i].midnight)) {
      theta = atan(planet[i].z / fabs(planet[i].x));
      flux[i] = UnoccultedFlux(planet[i].r, theta, planet[i].noon, planet[i].midnight, planet[i].nlat);
    } else {
      flux[i] = 0;
    }
    
    // Loop over all possible occultors
    occultor[i] = 0;
    for (j = 0; j < n; j++) {
      
      // Skip self
      if (i == j) continue;
      
      // Compute the planet-planet separation between
      // the ith planet and all the others and check
      // if the first planet is occulted.
      dx = (planet[i].x - planet[j].x);
      dy = (planet[i].y - planet[j].y);
      d = sqrt(dx * dx + dy * dy);
      if ((d <= planet[i].r + planet[j].r) && (planet[i].z < planet[j].z)){
        occultor[i] = j;
        if (planet[i].x < i) x0 = -dx;
        else x0 = dx;
        theta = atan(planet[i].z / fabs(planet[i].x));
        flux[i] -= OccultedFlux(planet[i].r, x0, dy, planet[j].r, theta, planet[i].noon, planet[i].midnight, planet[i].nlat);
      
      }
    }
  }
}