#include <sys/resource.h>
#include <tcl.h>
#include "communication.h"

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
    mpi_bcast_event(PARTICLE_CHANGED);		
    mpi_integrate(0);
  }
  markTime();
  return diffTime()/rds;
}
