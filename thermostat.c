/** \file thermostat.c
    Implementation of \ref thermostat.h "thermostat.h"
 */
#include <math.h>
#include "global.h"
#include "thermostat.h"
#include "particle_data.h"
#include "communication.h"
#include "random.h"

/** Friction coefficient gamma. */
double friction_gamma = 0;
double time_step;

int gamma_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;

  if (data < 0) {
    Tcl_AppendResult(interp, "Gamma must be non negativ.", (char *) NULL);
    return (TCL_ERROR);
  }
  friction_gamma = data;

  mpi_bcast_parameter(FIELD_GAMMA);

  return (TCL_OK);
}

void friction_thermo()
{
  int i, j;
  double temperature=1.0;
  for(i=0;i<n_particles;i++)
    for(j=0;j<3;j++)
            particles[i].f[j] = - friction_gamma/time_step*particles[i].v[j] 
                        + sqrt(12.0 * 2.0*temperature*friction_gamma/time_step)*(d_random()-0.5);
}
