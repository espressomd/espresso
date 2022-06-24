/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "tuning.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <utils/statistics/RunningAverage.hpp>

#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/min_element.hpp>

#include <mpi.h>

#include <algorithm>
#include <functional>
#include <string>

std::string TuningFailed::get_first_error() const {
  using namespace ErrorHandling;
  auto const queued_warnings = mpi_gather_runtime_errors_all(this_node == 0);
  auto message = std::string("tuning failed: an exception was thrown while "
                             "benchmarking the integration loop");
  for (auto const &warning : queued_warnings) {
    if (warning.level() == RuntimeError::ErrorLevel::ERROR) {
      message += " (" + warning.what() + ")";
      break;
    }
  }
  return message;
}

static void check_statistics(Utils::Statistics::RunningAverage<double> &acc) {
  if (acc.avg() <= 5 * MPI_Wtick()) {
    runtimeWarningMsg()
        << "Clock resolution is too low to reliably time integration.";
  }
  if (acc.sig() >= 0.1 * acc.avg()) {
    runtimeWarningMsg() << "Statistics of tuning samples is very bad.";
  }
}

static void throw_on_error() {
  auto const n_errors = check_runtime_errors_local();
  if (boost::mpi::all_reduce(comm_cart, n_errors, std::plus<>()) != 0) {
    throw TuningFailed{};
  }
}

double benchmark_integration_step(int int_steps) {
  Utils::Statistics::RunningAverage<double> running_average;

  integrate(0, 0);
  throw_on_error();

  /* perform force calculation test */
  for (int i = 0; i < int_steps; i++) {
    auto const tick = MPI_Wtime();
    integrate(0, -1);
    auto const tock = MPI_Wtime();
    running_average.add_sample((tock - tick));
    throw_on_error();
  }

  if (this_node == 0) {
    check_statistics(running_average);
  }

  /* MPI returns in seconds, returned value should be in ms. */
  auto retval = 1000. * running_average.avg();
  boost::mpi::broadcast(comm_cart, retval, 0);
  return retval;
}

/**
 * \brief Time the integration.
 * This times the integration and
 * propagates the system.
 *
 * @param int_steps Number of steps to integrate.
 * @return Time per integration in ms.
 */
static double time_calc(int int_steps) {
  if (mpi_integrate(0, 0))
    return -1;

  /* perform force calculation test */
  const double tick = MPI_Wtime();
  if (mpi_integrate(int_steps, -1))
    return -1;
  const double tock = MPI_Wtime();

  /* MPI returns s, return value should be in ms. */
  return 1000. * (tock - tick) / int_steps;
}

void tune_skin(double min_skin, double max_skin, double tol, int int_steps,
               bool adjust_max_skin) {

  double a = min_skin;
  double b = max_skin;

  /* The maximal skin is the remainder from the required cutoff to
   * the maximal range that can be supported by the cell system, but
   * never larger than half the box size. */
  double const max_permissible_skin =
      std::min(*boost::min_element(cell_structure.max_cutoff()) -
                   maximal_cutoff(n_nodes),
               0.5 * *boost::max_element(box_geo.length()));

  if (adjust_max_skin and max_skin > max_permissible_skin)
    b = max_permissible_skin;

  while (fabs(a - b) > tol) {
    mpi_set_skin(a);
    auto const time_a = time_calc(int_steps);

    mpi_set_skin(b);
    auto const time_b = time_calc(int_steps);

    if (time_a > time_b) {
      a = 0.5 * (a + b);
    } else {
      b = 0.5 * (a + b);
    }
  }
  auto const new_skin = 0.5 * (a + b);
  mpi_set_skin(new_skin);
}
