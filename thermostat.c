// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
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
#ifdef NPT
#include "pressure.h"
#endif

/** Friction coefficient gamma. */
double friction_gamma = 0.0;
double friction_g0 = 0.0;
double friction_gv = 0.0;

/** Temperature */
double temperature = -1.0;

static double pref1;
static double pref2;
#ifdef NPT
static double pref3;
static double pref4;
#endif

#ifdef ROTATION
static double friction_gamma_rotation;
static double pref2_rotation;
#endif

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

int g0_callback(Tcl_Interp *interp, void *_data) {
  double data = *(double *)_data;
  if (data < 0) { Tcl_AppendResult(interp, "Gamma_0 must be non negativ.", (char *) NULL); return (TCL_ERROR); }
  friction_g0 = data;
  mpi_bcast_parameter(FIELD_FRICTION_G0);
  return (TCL_OK);
}
int gv_callback(Tcl_Interp *interp, void *_data) {
  double data = *(double *)_data;
  if (data < 0) { Tcl_AppendResult(interp, "Gamma_V must be non negativ.", (char *) NULL); return (TCL_ERROR); }
  friction_gv = data;
  mpi_bcast_parameter(FIELD_FRICTION_GV);
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

  return (TCL_OK);
}

void thermo_init()
{
  pref1 = -friction_gamma/time_step;
  pref2 = sqrt(24.0*temperature*friction_gamma/time_step);
#ifdef NPT
  if (piston != 0.0) {
    pref1 = -0.5*friction_g0;
    pref2 = sqrt(12.0*temperature*friction_g0*time_step);
    pref3 = -0.5*friction_gv*inv_piston*0.5*time_step;
    pref4 = sqrt(12.0*temperature*friction_gv*time_step);
  }
#endif
#ifdef ROTATION 
  friction_gamma_rotation = friction_gamma/3;
  pref2_rotation = sqrt(24.0*temperature*friction_gamma_rotation/time_step);
#endif
  /* fprintf(stderr,"%d: pref1=%f, pref2=%f\n",this_node,pref1,pref2); */
}

void friction_thermo(Particle *p)
{
  int j;
  for ( j = 0 ; j < 3 ; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p->l.ext_flag & COORD_FIXED(j)))
#endif
      p->f.f[j] = pref1*p->m.v[j] + pref2*(d_random()-0.5);
  }

  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
  THERMO_TRACE(fprintf(stderr,"%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",this_node,p->p.identity,p->f.f[0],p->f.f[1],p->f.f[2]));
}

#ifdef NPT
double friction_thermo_NpT(void) {
  return (pref3*p_diff + pref4*(d_random()-0.5) );
}
#endif


#ifdef ROTATION
void friction_thermo_rotation(Particle *p)
{
  int j;
#ifdef EXTERNAL_FORCES
  if(!(p->l.ext_flag & COORDS_FIX_MASK))
#endif
    {
      for ( j = 0 ; j < 3 ; j++)
	p->f.torque[j] = -friction_gamma*p->m.omega[j] + pref2*(d_random()-0.5);

      ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
      THERMO_TRACE(fprintf(stderr,"%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",this_node,p->p.identity,p->f.f[0],p->f.f[1],p->f.f[2]));
    }
}
#endif

