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
#include <math.h>
#include <memory>
#include <vector>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;
constexpr const double v0 = 0.064;

double u_expected(double x, double t, double nu, double v_0, double h,
                  int k_max = 100) {
  double u = x / h - 0.5;
  for (int k = 1; k <= k_max; k++) {
    u += 1.0 / (M_PI * k) * exp(-4 * M_PI * M_PI * nu * k * k / (h * h) * t) *
         sin(2 * M_PI / h * k * x);
  }
  return v_0 * u;
}
Vector3i mpi_shape{};
BOOST_AUTO_TEST_CASE(test_lees_edwards) {
  double density = 1;
  double viscosity = 1. / 7.;
  auto lattice =
      std::make_shared<LatticeWalberla>(Vector3i{8, 64, 8}, mpi_shape, 1);
  auto lb = walberla::LBWalberlaImpl<double>(lattice, viscosity, density);
  auto le_pack = std::make_unique<LeesEdwardsPack>(
      0, 1, []() { return 0.0; }, [=]() { return v0; });
  lb.set_collision_model(std::move(le_pack));
  auto const grid_size_y = lattice->get_grid_dimensions()[1];
  for (int i = 0; i < 200; i++) {
    lb.integrate();
    if (i < grid_size_y / 2.)
      continue;
    for (double y : {-0.5, 0., 0.5, 5.5, 6.5, 7.5, 14.0, 14.5,
                     0.23 * grid_size_y, 0.7 * grid_size_y, 1. * grid_size_y}) {
      auto u = lb.get_velocity_at_pos(Vector3d{4, y, 4}, true);
      auto expected = u_expected(y, i, viscosity, v0, grid_size_y);
      std::cout << y << " " << *u << " " << expected << std::endl;
    }
    std::cout << std::endl;
  }
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
