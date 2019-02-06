/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
    Implementation of tuning.hpp .
*/
#include "communication.hpp"
#include "errorhandling.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "utils.hpp"
#include "utils/statistics/RunningAverage.hpp"
#include <limits>
#include <sys/resource.h>
#include <sys/time.h>

int timing_samples = 10;

/**
 * \brief Time the force calculation.
 * This times the force calculation without
 * propagating the system. It therefore does
 * not include e.g. Verlet list updates.
 *
 * @return Time per integration in ms.
 */
double time_force_calc(int default_samples) {
  int rds = timing_samples > 0 ? timing_samples : default_samples;
  int i;
  Utils::Statistics::RunningAverage<double> running_average;

  if (mpi_integrate(0, 0))
    return -1;

  /* perform force calculation test */
  for (i = 0; i < rds; i++) {
    const double tick = MPI_Wtime();

    if (mpi_integrate(0, -1))
      return -1;

    const double tock = MPI_Wtime();
    running_average.add_sample((tock - tick));
  }

  if (running_average.avg() <= 5 * MPI_Wtick()) {
    runtimeWarning("Clock resolution is too low to reliably time integration.");
  }

  if (running_average.sig() >= 0.1 * running_average.avg()) {
    runtimeWarning("Statistics of tuning samples is very bad.");
  }

  /* MPI returns s, return value should be in ms. */
  return 1000. * running_average.avg();
}

/**
 * \brief Time the integration.
 * This times the integration and
 * propagates the system.
 *
 * @param rds Number of steps to integrate.
 * @return Time per integration in ms.
 */
static double time_calc(int rds) {
  if (mpi_integrate(0, 0))
    return -1;

  /* perform force calculation test */
  const double tick = MPI_Wtime();
  if (mpi_integrate(rds, -1))
    return -1;
  const double tock = MPI_Wtime();

  /* MPI returns s, return value should be in ms. */
  return 1000. * (tock - tick) / rds;
}

void tune_skin(double min, double max, double tol, int steps) {
  skin_set = true;

  double a = min;
  double b = max;
  double time_a, time_b;

  while (fabs(a - b) > tol) {
    skin = a;
    mpi_bcast_parameter(FIELD_SKIN);
    time_a = time_calc(steps);

    skin = b;
    mpi_bcast_parameter(FIELD_SKIN);
    time_b = time_calc(steps);

    if (time_a > time_b) {
      a = 0.5 * (a + b);
    } else {
      b = 0.5 * (a + b);
    }
  }
  skin = 0.5 * (a + b);
  mpi_bcast_parameter(FIELD_SKIN);
}
