/*
 * Copyright (C) 2019-2023 The ESPResSo project
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
#define BOOST_TEST_MODULE Walberla interpolation test
#define BOOST_TEST_DYN_LINK
#include "config/config.hpp"

#ifdef WALBERLA

#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "tests_common_lb.hpp"

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <numeric>

using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

static LBTestParameters params; // populated in main()

BOOST_DATA_TEST_CASE(force_interpolation_bspline, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(params);

  /* Check the bspline weights sum up to 1 in each direction.
   * The position for the interpolation is sampled uniformly
   * in the range (-0.5, +0.5) around the LB node mid point.
   */

  constexpr auto dx = 0.02;
  auto const f = Vector3d{{-1.0, 0.5, 1.5}};
  auto offset = Vector3d::broadcast(-0.5 + dx);
  int index = 0;
  for (auto const &n : local_nodes_incl_ghosts(lb->get_lattice(), false)) {
    if (lb->get_lattice().node_in_local_halo(n)) {
      index = (index + 1) % 3;
      offset[index] = std::fmod(offset[index] + 0.5, 1. - dx) - 0.5 + dx;
      auto const pos = n + offset;
      lb->add_force_at_pos(pos, f);
      // Check neighboring nodes for bspline weights
      Vector3d sum{};
      for (int x : {0, 1}) {
        for (int y : {0, 1}) {
          for (int z : {0, 1}) {
            Vector3i const check_node{{n[0] - x, n[1] - y, n[2] - z}};
            if (lb->get_lattice().node_in_local_halo(check_node)) {
              auto const res = lb->get_node_force_to_be_applied(check_node);
              sum += *res;
            }
          }
        }
      }
      BOOST_CHECK_SMALL((sum - f).norm(), 1E-10);
      // Apply counter force to clear force field
      lb->add_force_at_pos(pos, -f);
    }
  }
}

BOOST_DATA_TEST_CASE(velocity_interpolation_bspline, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(params);

  /* Check linear interpolation of the velocity. LB cells can couple
   * to particles that are at most 1 agrid away from the cell mid point.
   * The test assigns a velocity to every third node on a simple cubic
   * lattice with lattice constant l0 = 3 * agrid. A particle moving on a
   * line along the x-, y- or z-axis should experience a coupling whose
   * profile is a series of peaks with formula:
   *   f(x) = sum_i v_i * max(0, 1 - abs(x - x_i))
   * where x_i are the peak centers and v_i their velocity.
   * In ASCII art: _/\_/\_/\_/\_
   */

  // make sure the lattice constant is commensurate with the box dimensions
  assert(params.grid_dimensions[0] % 3 == 0 and
         params.grid_dimensions[1] % 3 == 0 and
         params.grid_dimensions[2] % 3 == 0);

  // set node velocities on a simple cubic lattice
  auto const vel = Vector3d{{-1., 0.5, 1.5}};
  for (auto const &n : local_nodes_incl_ghosts(lb->get_lattice(), false)) {
    if (lb->get_lattice().node_in_local_domain(n)) {
      if ((n[0] + 2) % 3 == 0 and (n[1] + 2) % 3 == 0 and (n[2] + 2) % 3 == 0) {
        BOOST_CHECK(lb->set_node_velocity(n, vel));
      }
    }
  }

  lb->ghost_communication();

  for (double x = 0.0; x < params.box_dimensions[0]; x += 0.3) {
    for (double y = 0.1; y < params.box_dimensions[1]; y += 0.3) {
      for (double z = 0.2; z < params.box_dimensions[2]; z += 0.3) {
        Vector3d const pos{x, y, z};
        if (lb->get_lattice().pos_in_local_domain(pos)) {
          auto const factor = std::accumulate(
              pos.begin(), pos.end(), 1., [](double a, double x) {
                return a * std::max(0., 1. - std::fabs(std::fmod(x, 3.) - 1.5));
              });
          auto const ref = factor * vel;
          auto const res = lb->get_velocity_at_pos(pos, true);
          BOOST_CHECK(res); // locally available
          BOOST_CHECK_SMALL((*res - ref).norm(), 1E-10);
        }
      }
    }
  }
}

// TODO: check last applied force on a ghost node, i.e. when two forces
// are applied at (agrid/2, 0, 0) and (box_l - agrid/2, 0, 0)

int main(int argc, char **argv) {
  int n_nodes;
  Vector3i mpi_shape{};

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());
  walberla::mpi_init();

  params.seed = 0u;
  params.kT = 1.3E-4;
  params.viscosity = 0.003;
  params.density = 1.4;
  params.grid_dimensions = Vector3i{12, 6, 9};
  params.box_dimensions = Vector3d{12, 6, 9};
  params.lattice =
      std::make_shared<LatticeWalberla>(params.grid_dimensions, mpi_shape, 1u);

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // WALBERLA
int main(int argc, char **argv) {}
#endif
