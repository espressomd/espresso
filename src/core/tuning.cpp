/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file tuning.cpp
    Implementation of tuning.hpp .
*/
#include <sys/time.h>
#include <sys/resource.h>
#include "utils.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrate.hpp"
#include "global.hpp"
#include <limits>

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

  if (mpi_integrate(0, 0))
    return -1;

  /* perform force calculation test */
  markTime();
  for (i = 0; i < rds; i++) {
    if (mpi_integrate(0, -1))
      return -1;
  }
  markTime();
  return diffTime()/rds;
}

static double time_calc(int rds)
{
  if (mpi_integrate(0, 0))
    return -1;

  /* perform force calculation test */
  markTime();
    if (mpi_integrate(rds, -1))
      return -1;
  markTime();
  return diffTime()/rds;
}


void tune_skin(double min, double max, double tol, int steps) {
  skin_set = true;

  double a = min;
  double b = max;
  double time_a, time_b;

  while(fabs(a - b) > tol) {
    skin = a;
    mpi_bcast_parameter(FIELD_SKIN);    
    time_a = time_calc(steps);

    skin = b;
    mpi_bcast_parameter(FIELD_SKIN);    
    time_b = time_calc(steps);

    if(time_a > time_b) {
      a = 0.5*(a + b);
    } else {
      b = 0.5*(a + b);
    }
  }
  skin = 0.5*(a+b);
  mpi_bcast_parameter(FIELD_SKIN);    
}
