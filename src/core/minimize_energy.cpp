/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

#include "minimize_energy.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "rotation.hpp"
#include "utils.hpp"
#include <algorithm>
#include <limits>

#ifdef MINIMIZE_ENERGY_DEBUG
#define MINIMIZE_ENERGY_TRACE(A) A
#else
#define MINIMIZE_ENERGY_TRACE(A)
#endif

struct MinimizeEnergyParameters {
  double f_max;
  double gamma;
  int max_steps;
  double max_displacement;
};

static MinimizeEnergyParameters *params = 0;

/* Signum of val */
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

bool steepest_descent_step(void) {
  // Maximal force encountered on node
  double f_max = -std::numeric_limits<double>::max();
  // and globally
  double f_max_global;

  // Square of force,torque on particle
  double f, t;

  // Positional increments
  double dp, dp2, dp2_max = -std::numeric_limits<double>::max();

  // Iteration over all local particles

  for (auto &p : local_cells.particles()) {
    f = 0.0;
    t = 0.0;
    dp2 = 0.0;
    // For all Cartesian coordinates
    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      // Skip, if coordinate is fixed
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
#ifdef VIRTUAL_SITES
        // Skip positional increments of virtual particles
        if (!p.p.isVirtual)
#endif
        {
          // Square of force on particle
          f += SQR(p.f.f[j]);

          // Positional increment
          dp = params->gamma * p.f.f[j];
          if (fabs(dp) > params->max_displacement)
            // Crop to maximum allowed by user
            dp = sgn<double>(dp) * params->max_displacement;
          dp2 += SQR(dp);

          // Move particle
          p.r.p[j] += dp;
          MINIMIZE_ENERGY_TRACE(printf("part %d dim %d dp %e gamma*f %e\n", i,
                                       j, dp, params->gamma * p.f.f[j]));
        }
    }
#ifdef ROTATION
    // Rotational increment
    double dq[3]; // Vector parallel to torque

    for (int j = 0; j < 3; j++) {
      dq[j] = 0;
      // Square of torque
      t += SQR(p.f.torque[j]);

      // Rotational increment
      dq[j] = params->gamma * p.f.torque[j];
    }
    // Normalize rotation axis and compute amount of rotation
    double l = normr(dq);
    if (l > 0.0) {
      for (int j = 0; j < 3; j++)
        dq[j] /= l;

      if (fabs(l) > params->max_displacement)
        // Crop to maximum allowed by user
        l = sgn(l) * params->max_displacement;

      // Rotate the particle around axis dq by amount l
      rotate_particle(&(p), dq, l);
    }
#endif

    // Note maximum force/torque encountered
    f_max = std::max(f_max, f);
    f_max = std::max(f_max, t);
    dp2_max = std::max(dp2_max, dp2);
  }

  set_resort_particles(Cells::RESORT_LOCAL);

  MINIMIZE_ENERGY_TRACE(
      printf("f_max %e resort_particles %d\n", f_max, resort_particles));
  announce_resort_particles();

  // Synchronize maximum force/torque encountered
  MPI_Allreduce(&f_max, &f_max_global, 1, MPI_DOUBLE, MPI_MAX, comm_cart);

  // Return true, if the maximum force/torque encountered is below the user
  // limit.
  return (sqrt(f_max_global) < params->f_max);
}

void minimize_energy_init(const double f_max, const double gamma,
                          const int max_steps, const double max_displacement) {
  if (!params)
    params = new MinimizeEnergyParameters;

  params->f_max = f_max;
  params->gamma = gamma;
  params->max_steps = max_steps;
  params->max_displacement = max_displacement;
}

bool minimize_energy(void) {
  if (!params)
    params = new MinimizeEnergyParameters;

  MPI_Bcast(params, sizeof(MinimizeEnergyParameters), MPI_BYTE, 0, comm_cart);
  int integ_switch_old = integ_switch;
  integ_switch = INTEG_METHOD_STEEPEST_DESCENT;
  integrate_vv(params->max_steps, -1);
  integ_switch = integ_switch_old;

  return true;
}
