/*
  Copyright (C) 2010,2012 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file tuning.c
    Implementation of tuning.h .
*/
#include <sys/time.h>
#include <sys/resource.h>
#include "utils.h"
#include "communication.h"
#include "errorhandling.h"

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

double time_force_calc(int default_samples)
{
  int rds = timing_samples > 0 ? timing_samples : default_samples;
  int i;

  if (mpi_integrate(0))
    return -1;

  /* perform force calculation test */
  markTime();
  for (i = 0; i < rds; i++) {
    mpi_bcast_event(INVALIDATE_SYSTEM);
    if (mpi_integrate(0))
      return -1;
  }
  markTime();
  return diffTime()/rds;
}
