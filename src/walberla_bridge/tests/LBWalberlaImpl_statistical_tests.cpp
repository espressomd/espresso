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
#define BOOST_TEST_MODULE Walberla node setters and getters test
#define BOOST_TEST_DYN_LINK
#include "config.hpp"

#ifdef LB_WALBERLA

#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <LBWalberlaBase.hpp>
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include "tests_common.hpp"

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

LBTestParameters params; // populated in main()
Vector3i mpi_shape;      // populated in main

BOOST_DATA_TEST_CASE(velocity_fluctuation, bdata::make(thermalized_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);

  // Warmup
  for (int i = 0; i < 200; i++)
    lb->integrate();

  // Sample
  int steps = 800;

  auto my_left = lb->get_local_domain().first;
  auto my_right = lb->get_local_domain().second;

  Vector3d sum_v{}, sum_v_square{};
  auto const count =
      static_cast<int>((my_right[0] - my_left[0]) * (my_right[1] - my_left[1]) *
                       (my_right[2] - my_left[2]));

  for (int i = 0; i < steps; i++) {
    Vector3d step_v{}, step_v_square{};
    for (int x = static_cast<int>(my_left[0]);
         x < static_cast<int>(my_right[0]); x++) {
      for (int y = static_cast<int>(my_left[1]);
           y < static_cast<int>(my_right[1]); y++) {
        for (int z = static_cast<int>(my_left[2]);
             z < static_cast<int>(my_right[2]); z++) {
          const Vector3i node{{x, y, z}};
          auto v = *(lb->get_node_velocity(node));
          auto rho = *(lb->get_node_density(node));
          step_v += v * rho;
          step_v_square += rho * hadamard_product(v, v);
        }
      }
    }
    step_v /= double(count);
    step_v_square /= double(count);

    sum_v += step_v;
    sum_v_square += step_v_square;
    auto tmp = sum_v_square / double(i + 1);
    printf("%g %g %g\n", tmp[0], tmp[1], tmp[2]);

    lb->integrate();
    lb->integrate();
    lb->integrate();
  }

  // aggregate
  MPI_Allreduce(MPI_IN_PLACE, sum_v.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, sum_v_square.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  const int num_ranks = mpi_shape[0] * mpi_shape[1] * mpi_shape[2];
  sum_v /= double(num_ranks);
  sum_v_square /= double(num_ranks);
  // check
  const double tol_v = 3E-6;

  BOOST_CHECK_SMALL(std::abs(sum_v[0] / steps), tol_v * 100); // boost oddity
  BOOST_CHECK_SMALL(std::abs(sum_v[1] / steps), tol_v * 100);
  BOOST_CHECK_SMALL(std::abs(sum_v[2] / steps), tol_v * 100);

  const double tol_kT = 5; // this is in percent ...
  BOOST_CHECK_CLOSE(sum_v_square[0] / steps, params.kT, tol_kT);
  BOOST_CHECK_CLOSE(sum_v_square[1] / steps, params.kT, tol_kT);
  BOOST_CHECK_CLOSE(sum_v_square[2] / steps, params.kT, tol_kT);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int n_nodes;

  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());

  params.viscosity = 0.02;
  params.kT = 1.1E-4;
  params.density = 1.4;
  params.grid_dimensions = Vector3i{12, 12, 18};
  params.box_dimensions = Vector3d{6, 6, 9};

  walberla_mpi_init();
  auto res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
