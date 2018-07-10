/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group,

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
/** \file lb-boundaries.cpp
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb-boundaries.hpp.
 *
 */

#include "config.hpp"

#include <limits>

#include "communication.hpp"
#include "constraints.hpp"
#include "electrokinetics.hpp"
#include "electrokinetics_pdb_parse.hpp"
#include "interaction_data.hpp"
#include "lb.hpp"
#include "lbboundaries.hpp"
#include "lbboundaries/LBBoundary.hpp"
#include "lbgpu.hpp"
#include "utils.hpp"
#include <memory>
#include <vector>

namespace LBBoundaries {

std::vector<std::shared_ptr<LBBoundary>> lbboundaries;
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)

void lbboundary_mindist_position(const Vector3d& pos, double *mindist,
                                 double distvec[3], int *no) {

  double vec[3] = {std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};
  double dist = std::numeric_limits<double>::infinity();
  *mindist = std::numeric_limits<double>::infinity();

  int n = 0;
  for (auto lbb = lbboundaries.begin(); lbb != lbboundaries.end(); ++lbb, n++) {
    (**lbb).calc_dist(pos.data(), &dist, vec);

    if (dist < *mindist || lbb == lbboundaries.begin()) {
      *no = n;
      *mindist = dist;
      distvec[0] = vec[0];
      distvec[1] = vec[1];
      distvec[2] = vec[2];
    }
  }
}

/** Initialize boundary conditions for all constraints in the system. */
void lb_init_boundaries() {

  double pos[3];
  if (lattice_switch & LATTICE_LB_GPU) {
#if defined(LB_GPU) && defined(LB_BOUNDARIES_GPU)
    int number_of_boundnodes = 0;
    int *host_boundary_node_list = (int *)Utils::malloc(sizeof(int));
    int *host_boundary_index_list = (int *)Utils::malloc(sizeof(int));
    size_t size_of_index;
    int boundary_number =
        -1; // the number the boundary will actually belong to.

#ifdef EK_BOUNDARIES
    ekfloat *host_wallcharge_species_density = nullptr;
    float node_wallcharge = 0.0f;
    int wallcharge_species = -1, charged_boundaries = 0;
    int node_charged = 0;

    for (auto lbb = lbboundaries.begin(); lbb != lbboundaries.end(); ++lbb) {
      (**lbb).set_net_charge(0.0);
    }

    if (ek_initialized) {
      host_wallcharge_species_density = (ekfloat *)Utils::malloc(
          ek_parameters.number_of_nodes * sizeof(ekfloat));
      for (auto lbb = lbboundaries.begin(); lbb != lbboundaries.end(); ++lbb) {
        if ((**lbb).charge_density() != 0.0) {
          charged_boundaries = 1;
          break;
        }
      }
      if (pdb_charge_lattice) {
        charged_boundaries = 1;
      }

      for (int n = 0; n < int(ek_parameters.number_of_species); n++)
        if (ek_parameters.valency[n] != 0.0) {
          wallcharge_species = n;
          break;
        }

      if (wallcharge_species == -1 && charged_boundaries) {
        runtimeErrorMsg()
            << "no charged species available to create wall charge\n";
      }
    }
#endif
    for (int z = 0; z < int(lbpar_gpu.dim_z); z++) {
      for (int y = 0; y < int(lbpar_gpu.dim_y); y++) {
        for (int x = 0; x < int(lbpar_gpu.dim_x); x++) {
          pos[0] = (x + 0.5) * lbpar_gpu.agrid;
          pos[1] = (y + 0.5) * lbpar_gpu.agrid;
          pos[2] = (z + 0.5) * lbpar_gpu.agrid;

          double dist = 1e99;
          double dist_tmp = 0.0;
          double dist_vec[3];

#ifdef EK_BOUNDARIES
          if (ek_initialized) {
            host_wallcharge_species_density[ek_parameters.dim_y *
                                                ek_parameters.dim_x * z +
                                            ek_parameters.dim_x * y + x] = 0.0f;
            node_charged = 0;
            node_wallcharge = 0.0f;
          }
#endif

          int n = 0;
          for (auto lbb = lbboundaries.begin(); lbb != lbboundaries.end();
               ++lbb, n++) {
            (**lbb).calc_dist(pos, &dist_tmp, dist_vec);

            if (dist > dist_tmp || n == 0) {
              dist = dist_tmp;
              boundary_number = n;
            }
#ifdef EK_BOUNDARIES
            if (ek_initialized) {
              if (dist_tmp <= 0 && (**lbb).charge_density() != 0.0f) {
                node_charged = 1;
                node_wallcharge += (**lbb).charge_density() *
                                   ek_parameters.agrid * ek_parameters.agrid *
                                   ek_parameters.agrid;
                (**lbb).set_net_charge(
                    (**lbb).net_charge() +
                    (**lbb).charge_density() * ek_parameters.agrid *
                        ek_parameters.agrid * ek_parameters.agrid);
              }
            }
#endif
          }

#ifdef EK_BOUNDARIES
          if (pdb_boundary_lattice &&
              pdb_boundary_lattice[ek_parameters.dim_y * ek_parameters.dim_x *
                                       z +
                                   ek_parameters.dim_x * y + x]) {
            dist = -1;
            boundary_number = lbboundaries.size(); // Makes sure that
            // boundary_number is not used by
            // a constraint
          }
#endif
          if (dist <= 0 && boundary_number >= 0 &&
              (lbboundaries.size() > 0 || pdb_boundary_lattice)) {
            size_of_index = (number_of_boundnodes + 1) * sizeof(int);
            host_boundary_node_list =
                Utils::realloc(host_boundary_node_list, size_of_index);
            host_boundary_index_list =
                Utils::realloc(host_boundary_index_list, size_of_index);
            host_boundary_node_list[number_of_boundnodes] =
                x + lbpar_gpu.dim_x * y + lbpar_gpu.dim_x * lbpar_gpu.dim_y * z;
            host_boundary_index_list[number_of_boundnodes] =
                boundary_number + 1;
            number_of_boundnodes++;
          }

          lbpar_gpu.number_of_boundnodes = number_of_boundnodes;

#ifdef EK_BOUNDARIES
          if (ek_initialized) {
            ek_parameters.number_of_boundary_nodes = number_of_boundnodes;

            if (wallcharge_species != -1) {
              if (pdb_charge_lattice &&
                  pdb_charge_lattice[ek_parameters.dim_y * ek_parameters.dim_x *
                                         z +
                                     ek_parameters.dim_x * y + x] != 0.0f) {
                node_charged = 1;
                node_wallcharge +=
                    pdb_charge_lattice[ek_parameters.dim_y *
                                           ek_parameters.dim_x * z +
                                       ek_parameters.dim_x * y + x];
              }
              if (node_charged)
                host_wallcharge_species_density[ek_parameters.dim_y *
                                                    ek_parameters.dim_x * z +
                                                ek_parameters.dim_x * y + x] =
                    node_wallcharge / ek_parameters.valency[wallcharge_species];
              else if (dist <= 0)
                host_wallcharge_species_density[ek_parameters.dim_y *
                                                    ek_parameters.dim_x * z +
                                                ek_parameters.dim_x * y + x] =
                    0.0f;
              else
                host_wallcharge_species_density[ek_parameters.dim_y *
                                                    ek_parameters.dim_x * z +
                                                ek_parameters.dim_x * y + x] =
                    ek_parameters.density[wallcharge_species] *
                    ek_parameters.agrid * ek_parameters.agrid *
                    ek_parameters.agrid;
            }
          }
#endif
        }
      }
    }

    /**call of cuda fkt*/
    float *boundary_velocity =
        (float *)Utils::malloc(3 * (lbboundaries.size() + 1) * sizeof(float));
    int n = 0;
    for (auto lbb = lbboundaries.begin(); lbb != lbboundaries.end();
         ++lbb, n++) {
      boundary_velocity[3 * n + 0] = (**lbb).velocity()[0];
      boundary_velocity[3 * n + 1] = (**lbb).velocity()[1];
      boundary_velocity[3 * n + 2] = (**lbb).velocity()[2];
    }

    boundary_velocity[3 * lbboundaries.size() + 0] = 0.0f;
    boundary_velocity[3 * lbboundaries.size() + 1] = 0.0f;
    boundary_velocity[3 * lbboundaries.size() + 2] = 0.0f;

    if (lbboundaries.size() || pdb_boundary_lattice) {
      lb_init_boundaries_GPU(lbboundaries.size(), number_of_boundnodes,
                             host_boundary_node_list, host_boundary_index_list,
                             boundary_velocity);
    }
    free(boundary_velocity);
    free(host_boundary_node_list);
    free(host_boundary_index_list);

#ifdef EK_BOUNDARIES
    if (ek_initialized) {
      ek_init_species_density_wallcharge(host_wallcharge_species_density,
                                         wallcharge_species);
      free(host_wallcharge_species_density);
    }
#endif

#endif /* defined (LB_GPU) && defined (LB_BOUNDARIES_GPU) */
  } else {
#if defined(LB) && defined(LB_BOUNDARIES)
    int node_domain_position[3], offset[3];
    int the_boundary = -1;
    map_node_array(this_node, node_domain_position);

    offset[0] = node_domain_position[0] * lblattice.grid[0];
    offset[1] = node_domain_position[1] * lblattice.grid[1];
    offset[2] = node_domain_position[2] * lblattice.grid[2];

    for (int n = 0; n < lblattice.halo_grid_volume; n++) {
      lbfields[n].boundary = 0;
    }

    if (lblattice.halo_grid_volume == 0)
      return;

    for (int z = 0; z < lblattice.grid[2] + 2; z++) {
      for (int y = 0; y < lblattice.grid[1] + 2; y++) {
        for (int x = 0; x < lblattice.grid[0] + 2; x++) {
          pos[0] = (offset[0] + (x - 0.5)) * lblattice.agrid[0];
          pos[1] = (offset[1] + (y - 0.5)) * lblattice.agrid[1];
          pos[2] = (offset[2] + (z - 0.5)) * lblattice.agrid[2];

          double dist = 1e99;
          double dist_tmp = 0.0;
          double dist_vec[3];

          int n = 0;
          for (auto it = lbboundaries.begin(); it != lbboundaries.end();
               ++it, ++n) {
            (**it).calc_dist(pos, &dist_tmp, dist_vec);

            if (dist_tmp < dist || n == 0) {
              dist = dist_tmp;
              the_boundary = n;
            }
          }

          if (dist <= 0 && the_boundary >= 0 &&
              LBBoundaries::lbboundaries.size() > 0) {
            lbfields[get_linear_index(x, y, z, lblattice.halo_grid)].boundary =
                the_boundary + 1;
          } else {
            lbfields[get_linear_index(x, y, z, lblattice.halo_grid)].boundary =
                0;
          }
        }
      }
    }
#endif
  }
}

// TODO dirty hack. please someone get rid of void*
int lbboundary_get_force(void *lbb, double *f) {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)

  int no = 0;
  for (auto it = lbboundaries.begin(); it != lbboundaries.end(); ++it, ++no) {
    if (&(**it) == lbb)
      break;
  }
  if (no == lbboundaries.size())
    throw std::runtime_error("You probably tried to get the force of an "
                             "lbboundary that was not added to "
                             "system.lbboundaries.");

  std::vector<double> forces(3*lbboundaries.size());

  if (lattice_switch & LATTICE_LB_GPU) {
#if defined(LB_BOUNDARIES_GPU) && defined(LB_GPU)
    lb_gpu_get_boundary_forces(forces.data());

    f[0] = -forces[3 * no + 0];
    f[1] = -forces[3 * no + 1];
    f[2] = -forces[3 * no + 2];
#else
    return ES_ERROR;
#endif
  } else {
#if defined(LB_BOUNDARIES) && defined(LB)
    mpi_gather_stats(8, forces.data(), nullptr, nullptr, nullptr);

    f[0] = forces[3 * no + 0] * lbpar.agrid /
           lbpar.tau; // lbpar.tau; TODO this makes the units wrong and
    f[1] = forces[3 * no + 1] * lbpar.agrid /
           lbpar.tau; // lbpar.tau; the result correct. But it's 3.13AM
    f[2] = forces[3 * no + 2] * lbpar.agrid /
           lbpar.tau; // lbpar.tau; on a Saturday at the ICP. Someone fix.
#else
    return ES_ERROR;
#endif
  }

#endif
  return 0;
}

#endif /* LB_BOUNDARIES or LB_BOUNDARIES_GPU */

#ifdef LB_BOUNDARIES

void lb_bounce_back() {

#ifdef D3Q19
#ifndef PULL
  int k, i, l;
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0] * lblattice.halo_grid[1];
  int next[19];
  int x, y, z;
  double population_shift;
  double modes[19];
  next[0] = 0;                     // ( 0, 0, 0) =
  next[1] = 1;                     // ( 1, 0, 0) +
  next[2] = -1;                    // (-1, 0, 0)
  next[3] = yperiod;               // ( 0, 1, 0) +
  next[4] = -yperiod;              // ( 0,-1, 0)
  next[5] = zperiod;               // ( 0, 0, 1) +
  next[6] = -zperiod;              // ( 0, 0,-1)
  next[7] = (1 + yperiod);         // ( 1, 1, 0) +
  next[8] = -(1 + yperiod);        // (-1,-1, 0)
  next[9] = (1 - yperiod);         // ( 1,-1, 0)
  next[10] = -(1 - yperiod);       // (-1, 1, 0) +
  next[11] = (1 + zperiod);        // ( 1, 0, 1) +
  next[12] = -(1 + zperiod);       // (-1, 0,-1)
  next[13] = (1 - zperiod);        // ( 1, 0,-1)
  next[14] = -(1 - zperiod);       // (-1, 0, 1) +
  next[15] = (yperiod + zperiod);  // ( 0, 1, 1) +
  next[16] = -(yperiod + zperiod); // ( 0,-1,-1)
  next[17] = (yperiod - zperiod);  // ( 0, 1,-1)
  next[18] = -(yperiod - zperiod); // ( 0,-1, 1) +
  int reverse[] = {0, 2,  1,  4,  3,  6,  5,  8,  7, 10,
                   9, 12, 11, 14, 13, 16, 15, 18, 17};

  /* bottom-up sweep */
  //  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume;k++)
  for (z = 0; z < lblattice.grid[2] + 2; z++) {
    for (y = 0; y < lblattice.grid[1] + 2; y++) {
      for (x = 0; x < lblattice.grid[0] + 2; x++) {
        k = get_linear_index(x, y, z, lblattice.halo_grid);

        if (lbfields[k].boundary) {
          lb_calc_modes(k, modes);

          for (i = 0; i < 19; i++) {
            population_shift = 0;
            for (l = 0; l < 3; l++) {
              population_shift -=
                  lbpar.agrid * lbpar.agrid * lbpar.agrid * lbpar.rho * 2 *
                  lbmodel.c[i][l] * lbmodel.w[i] *
                  (*LBBoundaries::lbboundaries[lbfields[k].boundary - 1])
                      .velocity()[l] *
                  (lbpar.tau / lbpar.agrid) / // TODO
                  lbmodel.c_sound_sq;
            }

            if (x - lbmodel.c[i][0] > 0 &&
                x - lbmodel.c[i][0] < lblattice.grid[0] + 1 &&
                y - lbmodel.c[i][1] > 0 &&
                y - lbmodel.c[i][1] < lblattice.grid[1] + 1 &&
                z - lbmodel.c[i][2] > 0 &&
                z - lbmodel.c[i][2] < lblattice.grid[2] + 1) {
              if (!lbfields[k - next[i]].boundary) {
                for (l = 0; l < 3; l++) {
                  (*LBBoundaries::lbboundaries[lbfields[k].boundary - 1])
                      .force()[l] += // TODO
                      (2 * lbfluid[1][i][k] + population_shift) *
                      lbmodel.c[i][l];
                }
                lbfluid[1][reverse[i]][k - next[i]] =
                    lbfluid[1][i][k] + population_shift;
              } else {
                lbfluid[1][reverse[i]][k - next[i]] = lbfluid[1][i][k] = 0.0;
              }
            }
          }
        }
      }
    }
  }
#else
#error Bounce back boundary conditions are only implemented for PUSH scheme!
#endif
#else
#error Bounce back boundary conditions are only implemented for D3Q19!
#endif
}

#endif
} // namespace
