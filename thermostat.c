/** \file thermostat.c
    Implementation of \ref thermostat.h "thermostat.h"
 */
#include <math.h>
#include "global.h"
#include "thermostat.h"
#include "particle_data.h"
#include "communication.h"
#include "random.h"
#include "integrate.h"
#include "cells.h"
#include "debug.h"

/** Friction coefficient gamma. */
double friction_gamma = 0.0;
/** Temperature */
double temperature = -1.0;

static double pref1;
static double pref2;

int gamma_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;

  if (data < 0) {
    Tcl_AppendResult(interp, "Gamma must be non negativ.", (char *) NULL);
    return (TCL_ERROR);
  }
  friction_gamma = data;

  mpi_bcast_parameter(FIELD_GAMMA);
  mpi_bcast_event(PARAMETER_CHANGED);

  return (TCL_OK);
}

int temp_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;

  if (data < 0) {
    Tcl_AppendResult(interp, "Temperature must be non-negative.", (char *) NULL);
    return (TCL_ERROR);
  }
  temperature = data;

  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_event(PARAMETER_CHANGED);

  return (TCL_OK);
}

void thermo_init()
{
  pref1 = -friction_gamma/time_step;
  pref2 = sqrt(24.0*temperature*friction_gamma/time_step);
  /* fprintf(stderr,"%d: pref1=%f, pref2=%f\n",this_node,pref1,pref2); */
}

void friction_thermo(Particle *p)
{
#ifdef EXTERNAL_FORCES
  if(p->ext_flag != PARTICLE_FIXED) 
#endif
    {
      p->f[0] = pref1*p->v[0] + pref2*(d_random()-0.5);
      p->f[1] = pref1*p->v[1] + pref2*(d_random()-0.5);
      p->f[2] = pref1*p->v[2] + pref2*(d_random()-0.5);

      ONEPART_TRACE(if(p->r.identity==check_id) fprintf(stderr,"%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",this_node,p->f[0],p->f[1],p->f[2]));

      THERMO_TRACE(fprintf(stderr,"%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",this_node,p->r.identity,p->f[0],p->f[1],p->f[2]));
    }
}

#ifdef ROTATION
void friction_thermo_rotation(Particle *p)
{
#ifdef EXTERNAL_FORCES
  if(p->ext_flag != PARTICLE_FIXED) 
#endif
    {
      p->torque[0] = -friction_gamma*p->omega[0] + pref2*(d_random()-0.5);
      p->torque[1] = -friction_gamma*p->omega[1] + pref2*(d_random()-0.5);
      p->torque[2] = -friction_gamma*p->omega[2] + pref2*(d_random()-0.5);

      ONEPART_TRACE(if(p->r.identity==check_id) fprintf(stderr,"%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",this_node,p->f[0],p->f[1],p->f[2]));

      THERMO_TRACE(fprintf(stderr,"%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",this_node,p->r.identity,p->f[0],p->f[1],p->f[2]));
    }
}
#endif

