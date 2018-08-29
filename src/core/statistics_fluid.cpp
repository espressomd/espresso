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
/** \file statistics_fluid.cpp
 *
 * Fluid related analysis functions.
 * Implementation of \ref statistics_fluid.hpp.
 *
 */

#include "statistics_fluid.hpp"

#ifdef LB

#include "communication.hpp"
#include "lb.hpp"
#include "lbboundaries.hpp"
#include "utils.hpp"
#include <mpi.h>

/** Calculate mass of the LB fluid.
 * \param result Fluid mass
 */
void lb_calc_fluid_mass(double *result) {
  int x, y, z, index;
  double sum_rho = 0.0, rho = 0.0;

  for (x = 1; x <= lblattice.grid[0]; x++) {
    for (y = 1; y <= lblattice.grid[1]; y++) {
      for (z = 1; z <= lblattice.grid[2]; z++) {
        index = get_linear_index(x, y, z, lblattice.halo_grid);

        lb_calc_local_rho(index, &rho);
        // fprintf(stderr,"(%d,%d,%d) %e\n",x,y,z,rho);
        sum_rho += rho;
      }
    }
  }

  MPI_Reduce(&sum_rho, result, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}

/** Calculate momentum of the LB fluid.
 * \param result Fluid momentum
 */
void lb_calc_fluid_momentum(double *result) {

  int x, y, z, index;
  double j[3], momentum[3] = {0.0, 0.0, 0.0};

  for (x = 1; x <= lblattice.grid[0]; x++) {
    for (y = 1; y <= lblattice.grid[1]; y++) {
      for (z = 1; z <= lblattice.grid[2]; z++) {
        index = get_linear_index(x, y, z, lblattice.halo_grid);

        lb_calc_local_j(index, j);
        momentum[0] += j[0] + lbfields[index].force_density[0];
        momentum[1] += j[1] + lbfields[index].force_density[1];
        momentum[2] += j[2] + lbfields[index].force_density[2];
      }
    }
  }

  momentum[0] *= lbpar.agrid / lbpar.tau;
  momentum[1] *= lbpar.agrid / lbpar.tau;
  momentum[2] *= lbpar.agrid / lbpar.tau;

  MPI_Reduce(momentum, result, 3, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}

/** Calculate temperature of the LB fluid.
 * \param result Fluid temperature
 */
void lb_calc_fluid_temp(double *result) {
  int x, y, z, index;
  double rho, j[3];
  double temp = 0.0;
  int number_of_non_boundary_nodes = 0;

  for (x = 1; x <= lblattice.grid[0]; x++) {
    for (y = 1; y <= lblattice.grid[1]; y++) {
      for (z = 1; z <= lblattice.grid[2]; z++) {

        index = get_linear_index(x, y, z, lblattice.halo_grid);

#ifdef LB_BOUNDARIES
        if (!lbfields[index].boundary)
#endif
        {
          lb_calc_local_fields(index, &rho, j, nullptr);
          temp += scalar(j, j);
          number_of_non_boundary_nodes++;
        }
      }
    }
  }

  // @Todo: lblattice.agrid is 3d. What to use here?
  temp *= 1. /
          (3. * lbpar.rho * number_of_non_boundary_nodes * lbpar.tau *
           lbpar.tau * lbpar.agrid) /
          n_nodes;

  MPI_Reduce(&temp, result, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}

void lb_collect_boundary_forces(double *result) {
#ifdef LB_BOUNDARIES
  int n_lb_boundaries = LBBoundaries::lbboundaries.size();
  double *boundary_forces =
      (double *)Utils::malloc(3 * n_lb_boundaries * sizeof(double));

  int i = 0;
  for (auto it = LBBoundaries::lbboundaries.begin();
       it != LBBoundaries::lbboundaries.end(); ++it, i++)
    for (int j = 0; j < 3; j++)
      boundary_forces[3 * i + j] = (**it).force()[j];

  MPI_Reduce(boundary_forces, result, 3 * n_lb_boundaries, MPI_DOUBLE, MPI_SUM,
             0, comm_cart);
  free(boundary_forces);
#endif
}

#endif /* LB */
