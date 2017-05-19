#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"

void Flux(double time, int n, int nlam, PLANET planet[n], SETTINGS settings, double lambda[nlam], int occultor[n], double flux[n][nlam]){
  /*
  
  NOTE: This does not properly handle the pathological case of *multiple* simultaneous 
  planet-planet occultations involving the same planet!
  
  */

  double d, dx, dy;
  double x0, theta;
  double tmp[nlam];
  int i, j, m;
  
  // Compute the instantaneous orbital positions
  // of all the planets
  for (i = 0; i < n; i++)
    OrbitXYZ(time, &planet[i], settings);
  
  // Loop over all planets
  for (i = 0; i < n; i++) {
    
    // Compute the phase curve?
    if (settings.phasecurve) {
      theta = atan(planet[i].z / fabs(planet[i].x));
      UnoccultedFlux(planet[i].r, theta, planet[i].albedo, planet[i].irrad, planet[i].nlat, nlam, lambda, flux[i]);
    } else {
      for (m = 0; m < nlam; m++)
        flux[i][m] = 0;
    }

    // Loop over all possible occultors
    occultor[i] = -1;
    for (j = 0; j < n; j++) {
      
      // Skip self
      if (i == j) continue;
      
      // Compute the planet-planet separation between
      // the ith planet and all the others and check
      // if the first planet is occulted.
      dx = (planet[i].x - planet[j].x);
      dy = (planet[i].y - planet[j].y);
      d = sqrt(dx * dx + dy * dy);
      if ((d <= planet[i].r + planet[j].r) && (planet[i].z > planet[j].z)){
        occultor[i] = j;
        if (planet[i].x < i) x0 = -dx;
        else x0 = dx;
        theta = atan(planet[i].z / fabs(planet[i].x));
        OccultedFlux(planet[i].r, x0, dy, planet[j].r, theta, planet[i].albedo, planet[i].irrad, planet[i].nlat, nlam, lambda, tmp);
        for (m = 0; m < nlam; m++)
          flux[i][m] -= tmp[m];
        
      }
    }
  }
}