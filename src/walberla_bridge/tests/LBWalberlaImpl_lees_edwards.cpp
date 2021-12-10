/*
 * Copyright (C) 2019-2020 The ESPResSo project
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
#define BOOST_TEST_MODULE Walberla point force test
#define BOOST_TEST_DYN_LINK
#include "config.hpp"

#ifdef LB_WALBERLA

#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "tests_common.hpp"

#include "LBWalberlaImpl.hpp"
#include <LBWalberlaBase.hpp>
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

Vector3i mpi_shape{};
BOOST_AUTO_TEST_CASE(test_lees_edwards) {
  auto lattice =
      std::make_shared<LatticeWalberla>(Vector3i{64, 64, 64}, mpi_shape, 1);
  auto lb = walberla::LBWalberlaImpl<double>(lattice, 0.001, 1);
  double v0 = 0.064;
  lb.set_collision_model(LeesEdwardsPack(
      0, 1, [&]() { return 0.0; }, [&]() { return v0; }));
  for (int i = 0; i < 2500; i++) {
    lb.integrate();
    for (int j : {-1, 0, 1, 2, 3, 4, 5, 63, 64}) {
      std::cout << (*(lb.get_node_velocity(Vector3i{32, j, 32}, true)))[0]
                << " ";
    }
    std::cout << std::endl;
  }
  for (int i = -1; i <= 64; i++)
    std::cout << (*(lb.get_node_velocity(Vector3i{32, i, 32}, true)))[0]
              << std::endl;
}

int main(int argc, char **argv) {
  int n_nodes;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());
  walberla_mpi_init();

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
