// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#include <sys/time.h>
#include <sys/resource.h>
#include <tcl.h>
#include "communication.h"
#include "debug.h"

int timing_samples = 0;

/* timing helper variables */
static struct rusage time1, time2;

void markTime()
{
  time2 = time1;
  getrusage(RUSAGE_SELF, &time1);
}

double diffTime()
{
  return 1e-3*(time1.ru_utime.tv_usec - time2.ru_utime.tv_usec) +
    1e3*(time1.ru_utime.tv_sec - time2.ru_utime.tv_sec);
}

int timings_callback(Tcl_Interp *interp, void *data)
{
  if (*(int *)data <= 0)
    timing_samples = 0;
  else 
    timing_samples = *(int *)data;
  return TCL_OK;
}

double time_force_calc(int default_samples)
{
  int rds = timing_samples > 0 ? timing_samples : default_samples;
  int i;

  mpi_integrate(0);

  /* perform force calculation test */
  markTime();
  for (i = 0; i < rds; i++) {
    mpi_bcast_event(INVALIDATE_SYSTEM);
    mpi_integrate(0);
  }
  markTime();
  return diffTime()/rds;
}
