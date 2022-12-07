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
#define BOOST_TEST_MODULE Walberla statistical tests
#define BOOST_TEST_DYN_LINK
#include "config/config.hpp"

#ifdef WALBERLA

#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "tests_common.hpp"

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <cmath>
#include <functional>
#include <iostream>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

static LBTestParameters params; // populated in main()

BOOST_DATA_TEST_CASE(velocity_fluctuation, bdata::make(thermalized_lbs()),
                     lb_generator) {
  auto lb = lb_generator(params);

  // Warmup
  for (int i = 0; i < 200; i++)
    lb->integrate();

  // Sample
  int steps = 800;

  auto const [my_left, my_right] = lb->get_lattice().get_local_domain();
  auto const denominator = Utils::product(my_right - my_left);

  Vector3d sum_v_local{}, sum_v_square_local{};

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
    step_v /= denominator;
    step_v_square /= denominator;

    sum_v_local += step_v;
    sum_v_square_local += step_v_square;
    std::cout << sum_v_square_local / static_cast<double>(i + 1) << std::endl;

    lb->integrate();
    lb->integrate();
    lb->integrate();
  }

  // aggregate
  boost::mpi::communicator world;
  auto sum_v =
      boost::mpi::all_reduce(world, sum_v_local, std::plus<Vector3d>());
  auto sum_v_square =
      boost::mpi::all_reduce(world, sum_v_square_local, std::plus<Vector3d>());
  sum_v /= static_cast<double>(world.size());
  sum_v_square /= static_cast<double>(world.size());

  // check
  auto const tol_v = 3E-6;
  BOOST_CHECK_SMALL(std::abs(sum_v[0] / steps), tol_v * 100); // boost oddity
  BOOST_CHECK_SMALL(std::abs(sum_v[1] / steps), tol_v * 100);
  BOOST_CHECK_SMALL(std::abs(sum_v[2] / steps), tol_v * 100);

  const double tol_kT = 5; // this is in percent ...
  BOOST_CHECK_CLOSE(sum_v_square[0] / steps, params.kT, tol_kT);
  BOOST_CHECK_CLOSE(sum_v_square[1] / steps, params.kT, tol_kT);
  BOOST_CHECK_CLOSE(sum_v_square[2] / steps, params.kT, tol_kT);
}

int main(int argc, char **argv) {
  int n_nodes;
  Vector3i mpi_shape{};

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());
  walberla_mpi_init();

  params.seed = 1u;
  params.kT = 1.1E-4;
  params.viscosity = 0.02;
  params.density = 1.4;
  params.grid_dimensions = Vector3i{12, 12, 18};
  params.box_dimensions = Vector3d{6, 6, 9};
  params.lattice =
      std::make_shared<LatticeWalberla>(params.grid_dimensions, mpi_shape, 1u);

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // WALBERLA
int main(int argc, char **argv) {}
#endif
