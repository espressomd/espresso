#include "thermostat.h"
#include "global.h"

/** Friction coefficient gamma. */
double friction_gamma = 0;

void friction_thermo()
{
  int i, j;

  for(i=0;i<n_particles+n_ghosts;i++)
    for(j=0;j<3;j++)
	particles[i].f[j] = -friction_gamma*particles[i].v[j];
}
