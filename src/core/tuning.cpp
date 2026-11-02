/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"

#include <utils/statistics/RunningAverage.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

std::string TuningFailed::get_first_error() {
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

static void run_integration(System::System &system, int n_steps,
                            int reuse_forces) {
  auto const error_code = system.integrate(n_steps, reuse_forces);
  if (error_code == INTEG_ERROR_RUNTIME) {
    throw TuningFailed{};
  }
}

static void run_full_force_calc(System::System &system, int reuse_forces) {
  run_integration(system, 0, reuse_forces);
}

double benchmark_integration_step(System::System &system, int int_steps,
                                  bool with_integration) {
  Utils::Statistics::RunningAverage<double> running_average;

  // check if the system can be integrated with the current parameters;
  // this also computes the forces for the current positions
  run_full_force_calc(system, INTEG_REUSE_FORCES_CONDITIONALLY);

  if (with_integration) {
    // measure the time of full integration steps: the particles move, so
    // Verlet list updates and particle propagation are included; the forces
    // from the previous step are reused at the start of each step
    for (int i = 0; i < int_steps; i++) {
      auto const tick = MPI_Wtime();
      run_integration(system, 1, INTEG_REUSE_FORCES_CONDITIONALLY);
      auto const tock = MPI_Wtime();
      running_average.add_sample((tock - tick));
    }
  } else {
    // measure force calculation time on fixed particle positions
    for (int i = 0; i < int_steps; i++) {
      auto const tick = MPI_Wtime();
      run_full_force_calc(system, INTEG_REUSE_FORCES_NEVER);
      auto const tock = MPI_Wtime();
      running_average.add_sample((tock - tick));
    }
  }

  if (this_node == 0) {
    check_statistics(running_average);
  }

  /* MPI returns in seconds, returned value should be in ms. */
  auto retval = 1000. * running_average.avg();
  boost::mpi::broadcast(::comm_cart, retval, 0);
  return retval;
}

/**
 * \brief Time the integration.
 * This times the integration and
 * propagates the system.
 *
 * @param system The system to tune.
 * @param int_steps Number of steps to integrate.
 * @return Time per integration in ms.
 */
static double time_calc(System::System &system, int int_steps) {
  auto const error_code_init =
      system.integrate(0, INTEG_REUSE_FORCES_CONDITIONALLY);
  if (error_code_init == INTEG_ERROR_RUNTIME) {
    return -1;
  }

  /* perform force calculation test */
  auto const tick = MPI_Wtime();
  auto const error_code = system.integrate(int_steps, INTEG_REUSE_FORCES_NEVER);
  auto const tock = MPI_Wtime();
  if (error_code == INTEG_ERROR_RUNTIME) {
    return -1;
  }

  /* MPI returns in seconds, returned value should be in ms. */
  return 1000. * (tock - tick) / int_steps;
}

namespace System {
double System::get_max_permissible_skin() const {
  /* The maximal skin is the remainder from the required cutoff to
   * the maximal range that can be supported by the cell system, but
   * never larger than half the box size. */
  return std::min(
      std::ranges::min(cell_structure->max_cutoff()) - maximal_cutoff(),
      0.5 * std::ranges::max(box_geo->length()));
}

void System::tune_verlet_skin(double min_skin, double max_skin, double tol,
                              int int_steps, bool adjust_max_skin) {

  double a = min_skin;
  double b = max_skin;

  auto const max_permissible_skin = get_max_permissible_skin();
  if (adjust_max_skin and max_skin > max_permissible_skin)
    b = max_permissible_skin;

  if (b < a) {
    throw std::runtime_error(
        "Cannot tune the Verlet skin: the maximal permissible skin " +
        std::to_string(b) + " is smaller than the minimal skin " +
        std::to_string(a));
  }

  /* Minimize the time per step: coarse scan of interior skins, then
   * repeated halving of the bracket around the best sample. The best
   * measured skin is kept. The timings of the head node are used on all
   * ranks, so that they all take the same decisions. */
  auto best_skin = 0.5 * (a + b);
  auto best_time = std::numeric_limits<double>::max();
  auto const evaluate = [&](double skin) {
    cell_structure->set_verlet_skin(skin);
    auto time = time_calc(*this, int_steps);
    boost::mpi::broadcast(::comm_cart, time, 0);
    if (time >= 0. and time < best_time) {
      best_time = time;
      best_skin = skin;
    }
  };

  if (b - a <= tol) {
    evaluate(best_skin);
  } else {
    // up to six interior points, but not finer than the tolerance
    auto const n_coarse =
        std::clamp(static_cast<int>(std::ceil((b - a) / tol)) - 1, 1, 6);
    auto const h = (b - a) / (n_coarse + 1);
    for (int k = 1; k <= n_coarse; ++k) {
      evaluate(a + k * h);
    }
    auto lo = std::max(a, best_skin - h);
    auto hi = std::min(b, best_skin + h);
    while (hi - lo > tol) {
      auto const centre = best_skin;
      auto const m1 = 0.5 * (lo + centre);
      auto const m2 = 0.5 * (centre + hi);
      evaluate(m1);
      evaluate(m2);
      if (best_skin == m1) {
        hi = centre;
      } else if (best_skin == m2) {
        lo = centre;
      } else {
        lo = m1;
        hi = m2;
      }
    }
  }
  if (this_node == 0 and best_time == std::numeric_limits<double>::max()) {
    runtimeWarningMsg() << "Verlet skin tuning failed: no integration could "
                           "be timed, the skin was set to the middle of the "
                           "interval without a measurement";
  }
  cell_structure->set_verlet_skin(best_skin);
}
} // namespace System
