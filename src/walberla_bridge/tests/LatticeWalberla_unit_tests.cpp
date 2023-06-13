/*
 * Copyright (C) 2020-2023 The ESPResSo project
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
#define BOOST_TEST_MODULE LatticeWalberla tests
#define BOOST_TEST_DYN_LINK
#include "config/config.hpp"

#ifdef WALBERLA

#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "tests_common.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <functional>
#include <stdexcept>
#include <type_traits>

using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

static LatticeTestParameters params; // populated in main()
static Vector3i mpi_shape;           // populated in main

BOOST_DATA_TEST_CASE(domain_and_halo, bdata::xrange(3u), n_ghost_layers) {
  auto const lattice =
      LatticeWalberla(params.grid_dimensions, mpi_shape, n_ghost_layers);
  auto const [my_left, my_right] = lattice.get_local_domain();

  for (auto const &n : all_nodes_incl_ghosts(lattice)) {
    auto const pos = n + Vector3d::broadcast(.5);
    int is_local = 0;
    // Nodes in local domain
    if (Vector3d(n) >= my_left and Vector3d(n) < my_right) {
      BOOST_CHECK(lattice.node_in_local_domain(n));
      BOOST_CHECK(lattice.node_in_local_halo(n));

      BOOST_CHECK(lattice.pos_in_local_domain(pos));
      BOOST_CHECK(lattice.pos_in_local_halo(pos));
      is_local = 1;
    } else {
      // in local halo?
      if ((n + Vector3d::broadcast(n_ghost_layers)) >= my_left and
          (n - Vector3d::broadcast(n_ghost_layers)) < my_right) {
        BOOST_CHECK(!lattice.node_in_local_domain(n));
        BOOST_CHECK(lattice.node_in_local_halo(n));

        BOOST_CHECK(!lattice.pos_in_local_domain(pos));
        BOOST_CHECK(lattice.pos_in_local_halo(pos));
      } else {
        // neither in domain nor in halo
        BOOST_CHECK(!lattice.node_in_local_domain(n));
        BOOST_CHECK(!lattice.node_in_local_halo(n));

        BOOST_CHECK(!lattice.pos_in_local_domain(pos));
        BOOST_CHECK(!lattice.pos_in_local_halo(pos));
      }
    }

    // If the cell is in the global physical domain
    // check that only one mpi rank said the node was local
    constexpr auto origin = Vector3i{0, 0, 0};
    if (n >= origin and n < params.grid_dimensions) {
      boost::mpi::communicator world;
      auto const is_local_sum =
          boost::mpi::all_reduce(world, is_local, std::plus<int>());
      BOOST_CHECK(is_local_sum == 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(exceptions) {
  for (int i : {0, 1, 2}) {
    auto node_grid = Vector3i::broadcast(1);
    auto grid_dims = Vector3i::broadcast(1);
    grid_dims[i] = 3;
    node_grid[i] = 2;
    BOOST_CHECK_THROW(LatticeWalberla(grid_dims, node_grid, 1u),
                      std::runtime_error);
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int n_nodes;

  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());

  params.grid_dimensions = Vector3i{12, 12, 18};
  params.box_dimensions = Vector3d{12, 12, 18};

  walberla::mpi_init();
  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // WALBERLA
int main(int argc, char **argv) {}
#endif
