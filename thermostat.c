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

/** Friction coefficient gamma. */
double friction_gamma = 0.;
/** Temperature */
double temperature = 1.8;

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

  return (TCL_OK);
}

int temp_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;

  if (data < 0) {
    Tcl_AppendResult(interp, "Temperature must be non negativ.", (char *) NULL);
    return (TCL_ERROR);
  }
  temperature = data;

  mpi_bcast_parameter(FIELD_TEMPERATURE);

  return (TCL_OK);
}

void thermo_init()
{
  pref1 = -friction_gamma/time_step;
  pref2 = sqrt(12.0 * 2.0*temperature*friction_gamma/time_step);
}

void friction_thermo(Particle *p)
{
  p->f[0] = pref1*p->v[0] + pref2*(d_random()-0.5);
  p->f[1] = pref1*p->v[1] + pref2*(d_random()-0.5);
  p->f[2] = pref1*p->v[2] + pref2*(d_random()-0.5);
}
