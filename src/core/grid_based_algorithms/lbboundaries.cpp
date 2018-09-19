/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#include "communication.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/electrokinetics_pdb_parse.hpp"
#include "grid_based_algorithms/lb.hpp"
#include "grid_based_algorithms/lbboundaries.hpp"
#include "grid_based_algorithms/lbgpu.hpp"
#include "initialize.hpp"
#include "lbboundaries/LBBoundary.hpp"

#include <boost/range/algorithm.hpp>

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

namespace LBBoundaries {
constexpr const int no_boundary = -1;

std::vector<int> raster(std::vector<const LBBoundary *> const &boundaries,
                        Vector3d const &offset) {
  std::vector<int> ids(lblattice.halo_grid_volume, no_boundary);

  for (int z = 0; z < lblattice.grid[2] + 2; z++) {
    for (int y = 0; y < lblattice.grid[1] + 2; y++) {
      for (int x = 0; x < lblattice.grid[0] + 2; x++) {
        double pos[3];
        pos[0] = x * lblattice.agrid[0] - offset[0];
        pos[1] = y * lblattice.agrid[1] - offset[1];
        pos[2] = z * lblattice.agrid[2] - offset[2];

        auto it = boost::find_if(boundaries, [&pos](const LBBoundary *b) {
          double dist;
          double dist_vec[3];

          b->calc_dist(pos, &dist, dist_vec);

          return dist <= 0.;
        });

        if (it != boundaries.end()) {
          ids[get_linear_index(x, y, z, lblattice.halo_grid)] =
              std::distance(boundaries.begin(), it);
        }
      }
    }
  }

  return ids;
}
}

namespace LBBoundaries {

std::vector<std::shared_ptr<LBBoundary>> lbboundaries;
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)

void add(const std::shared_ptr<LBBoundary> &b) {
  lbboundaries.emplace_back(b);

  on_lbboundary_change();
}

void remove(const std::shared_ptr<LBBoundary> &b) {
  auto &lbb = lbboundaries;

  lbboundaries.erase(std::remove(lbb.begin(), lbb.end(), b), lbb.end());

  on_lbboundary_change();
}

void lbboundary_mindist_position(const Vector3d &pos, double *mindist,
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

    lb_init_boundaries_GPU(lbboundaries.size(), number_of_boundnodes,
                           host_boundary_node_list, host_boundary_index_list,
                           boundary_velocity);

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
    std::vector<const LBBoundary *> boundaries(lbboundaries.size());
    boost::transform(
        lbboundaries, boundaries.begin(),
        [](std::shared_ptr<LBBoundary> const &b) { return b.get(); });

    Vector3d offset;
    for (int i = 0; i < 3; i++)
      offset[i] = (-node_pos[i] * lblattice.grid[i] + 0.5) * lblattice.agrid[i];

    auto const ids = raster(boundaries, offset);
    boost::transform(ids, lbfields, lbfields.begin(),
                     [&](int id, LB_FluidNode node) {
                       if (id == no_boundary) {
                         node.boundary = 0;
                         node.boundary_velocity = Vector3d{};
                       } else {
                         node.boundary = id + 1;
                         node.boundary_velocity = (lbpar.tau / lbpar.agrid) *
                                                  boundaries[id]->velocity();
                       }
                       return node;
                     });
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

  std::vector<double> forces(3 * lbboundaries.size());

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

namespace {
static constexpr const std::array<std::array<int, 3>, 19> c = {{{{0, 0, 0}},
                                                                {{1, 0, 0}},
                                                                {{-1, 0, 0}},
                                                                {{0, 1, 0}},
                                                                {{0, -1, 0}},
                                                                {{0, 0, 1}},
                                                                {{0, 0, -1}},
                                                                {{1, 1, 0}},
                                                                {{-1, -1, 0}},
                                                                {{1, -1, 0}},
                                                                {{-1, 1, 0}},
                                                                {{1, 0, 1}},
                                                                {{-1, 0, -1}},
                                                                {{1, 0, -1}},
                                                                {{-1, 0, 1}},
                                                                {{0, 1, 1}},
                                                                {{0, -1, -1}},
                                                                {{0, 1, -1}},
                                                                {{0, -1, 1}}}};

constexpr std::array<int, 19> push_stencil(int yperiod, int zperiod) {
  return {0,                     // ( 0, 0, 0) =
          1,                     // ( 1, 0, 0) +
          -1,                    // (-1, 0, 0)
          yperiod,               // ( 0, 1, 0) +
          -yperiod,              // ( 0,-1, 0)
          zperiod,               // ( 0, 0, 1) +
          -zperiod,              // ( 0, 0,-1)
          (1 + yperiod),         // ( 1, 1, 0) +
          -(1 + yperiod),        // (-1,-1, 0)
          (1 - yperiod),         // ( 1,-1, 0)
          -(1 - yperiod),        // (-1, 1, 0) +
          (1 + zperiod),         // ( 1, 0, 1) +
          -(1 + zperiod),        // (-1, 0,-1)
          (1 - zperiod),         // ( 1, 0,-1)
          -(1 - zperiod),        // (-1, 0, 1) +
          (yperiod + zperiod),   // ( 0, 1, 1) +
          -(yperiod + zperiod),  // ( 0,-1,-1)
          (yperiod - zperiod),   // ( 0, 1,-1)
          -(yperiod - zperiod)}; // ( 0,-1, 1) +
}
}

void lb_bounce_back(LB_Fluid &lbfluid) {
  static constexpr int reverse[] = {0, 2,  1,  4,  3,  6,  5,  8,  7, 10,
                                    9, 12, 11, 14, 13, 16, 15, 18, 17};
  const int yperiod = lblattice.halo_grid[0];
  const int zperiod = lblattice.halo_grid[0] * lblattice.halo_grid[1];
  auto const next = push_stencil(yperiod, zperiod);

  std::array<Vector3d, 19> pop_shifts;
  boost::transform(lbmodel.c, lbmodel.w, pop_shifts.begin(),
                   [](std::array<double, 3> const &c_i, double w_i) {
                     return (lbpar.agrid * lbpar.agrid * lbpar.agrid *
                             lbpar.rho * 2. / lbmodel.c_sound_sq) *
                            w_i * Vector3d{c_i[0], c_i[1], c_i[2]};
                   });

  /* bottom-up sweep */
  for (int z = 0; z < lblattice.halo_grid[2]; z++) {
    for (int y = 0; y < lblattice.halo_grid[1]; y++) {
      for (int x = 0; x < lblattice.halo_grid[0]; x++) {
        auto const k = get_linear_index(x, y, z, lblattice.halo_grid);
        if (lbfields[k].boundary) {
          for (int i = 0; i < 19; i++) {
            if (x - c[i][0] > 0 && x - c[i][0] < lblattice.grid[0] + 1 &&
                y - c[i][1] > 0 && y - c[i][1] < lblattice.grid[1] + 1 &&
                z - c[i][2] > 0 && z - c[i][2] < lblattice.grid[2] + 1) {
              auto const stream_to = k - next[i];
              if (!lbfields[stream_to].boundary) {
                auto const population_shift =
                    -(pop_shifts[i] * lbfields[k].boundary_velocity);

                lbfluid[reverse[i]][stream_to] =
                    lbfluid[i][k] + population_shift;
              } else {
                lbfluid[reverse[i]][stream_to] = lbfluid[i][k] = 0.0;
              }
            }
          }
        }
      }
    }
  }
}

#endif
} // namespace LBBoundaries
