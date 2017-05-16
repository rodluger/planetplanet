#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"

double Flux(double time, PLANET *planet1, PLANET *planet2, SETTINGS *settings){
  /*
  
  */
  
  double d, dx, dy;
  double r, x0, y0, ro, theta, noon, midnight;
  double total;
  int nlat;
  
  // Compute the instantaneous orbital positions
  OrbitXYZ(time, planet1, settings);
  OrbitXYZ(time, planet2, settings);
  
  // Compute the planet-planet separation
  dx = (planet1->x - planet2->x);
  dy = (planet1->y - planet2->y);
  d = sqrt(dx * dx + dy * dy);
  
  // Do they occult?
  if (d <= planet1->r + planet2->r) {
    if (planet1->z < planet2->z) {
      // Planet 1 is occulted
      r = planet1->r;
      if (planet1->x < 0)
        x0 = -dx;
      else
        x0 = dx;
      y0 = dy;
      ro = planet2->r;
      theta = atan(planet1->z / fabs(planet1->x));
      noon = planet1->noon;
      midnight = planet1->midnight;
      nlat = planet1->nlat;
    } else {
      // Planet 2 is occulted
      r = planet2->r;
      if (planet1->x < 0)
        x0 = dx;
      else
        x0 = -dx;
      y0 = -dy;
      ro = planet1->r;
      theta = atan(planet2->z / fabs(planet2->x));
      noon = planet2->noon;
      midnight = planet2->midnight;
      nlat = planet2->nlat;
    }

    return -OccultedFlux(r, x0, y0, ro, theta, noon, midnight, nlat);
    
  }
  
  return 0;
  
}