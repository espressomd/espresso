
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

BOOST_DATA_TEST_CASE(integrate_with_point_force_thermalized,
                     bdata::make(thermalized_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);

  // Check that momentum stays zero after initial integration
  lb->integrate();
  lb->integrate();
  Vector3d mom = lb->get_momentum();
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  BOOST_CHECK_SMALL(mom.norm(), 1E-10);

  // Check that momentum changes as expected when applying forces
  // auto f = Vector3d{0.15, 0.25, -0.22};
  // auto f = Vector3d{0.0006, -0.0013, 0.000528};
  auto f = Vector3d{-0., 0., 0.};
  auto f2 = Vector3d{0.1, 0.2, -0.3};
  lb->set_external_force(f);
  Vector3d force_location{1.5, 1.5, 1.5};
  Vector3i force_node{int(force_location[0]), int(force_location[1]),
                      int(force_location[2])};

  lb->add_force_at_pos(force_location, f2);
  lb->integrate();
  for (auto n : all_nodes_incl_ghosts(params.grid_dimensions, 1)) {
    if (lb->node_in_local_halo(n)) {
      if (n == force_node) {
        BOOST_CHECK_SMALL(
            (*(lb->get_node_last_applied_force(n, true)) - f - f2).norm(),
            1E-10);
      } else {
        BOOST_CHECK_SMALL(
            (*(lb->get_node_last_applied_force(n, true)) - f).norm(), 1E-10);
      }
    }
  }
  mom = lb->get_momentum();

  // Expected momentum = momentum added in prev. time step
  // + f/2 from velocity shift due to last applied forces
  auto mom_exp = 1.5 * f * params.grid_dimensions[0] *
                     params.grid_dimensions[1] * params.grid_dimensions[2] +
                 1.5 * f2;
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  auto d = mom - mom_exp;
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
  printf("thermalized: %g %g %g | %g %g %g | %g %g %g\n", mom[0], mom[1],
         mom[2], mom_exp[0], mom_exp[1], mom_exp[2], d[0], d[1], d[2]);

  // check that momentum doesn't drift when no force is applied again
  lb->set_external_force(Vector3d{});
  // The expected moment is just that applied during a single time step
  // No f/2 correction, since no force was applied in last time step
  mom_exp = 1 * f * params.grid_dimensions[0] * params.grid_dimensions[1] *
                params.grid_dimensions[2] +
            1 * f2;
  lb->integrate();
  mom = lb->get_momentum();
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
}

// this can be merged with the thermalized test, once that passes
BOOST_DATA_TEST_CASE(integrate_with_point_force_unthermalized,
                     bdata::make(unthermalized_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);

  // Check that momentum stays zero after initial integration
  lb->integrate();
  BOOST_CHECK_SMALL(lb->get_momentum().norm(), 1E-10);

  // Check that momentum changes as expected when applying forces
  // auto f = Vector3d{0.0006, -0.0013, 0.000528};
  Vector3d f{};
  auto f2 = Vector3d{0.095, 0.23, -0.52};
  lb->set_external_force(f);
  lb->add_force_at_pos(Utils::Vector3d{2, 2, 2}, f2);
  lb->integrate();
  auto mom = lb->get_momentum();

  // Expected momentum = momentum added in prev. time step
  // + f/2 from velocity shift due to last applied forces
  auto mom_exp = 1.5 * f * params.grid_dimensions[0] *
                     params.grid_dimensions[1] * params.grid_dimensions[2] +
                 1.5 * f2;
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  auto d = mom - mom_exp;
  printf("%g %g %g | %g %g %g | %g %g %g\n", mom[0], mom[1], mom[2], mom_exp[0],
         mom_exp[1], mom_exp[2], d[0], d[1], d[2]);
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);

  // check that momentum doesn't drift when no force is applied again
  lb->set_external_force(Vector3d{});
  lb->integrate();
  // The expected moment is just that applied during a single time step
  // No f/2 correction, since no force was applied in last time step
  mom_exp = 1 * f * params.grid_dimensions[0] * params.grid_dimensions[1] *
                params.grid_dimensions[2] +
            1 * f2;
  mom = lb->get_momentum();
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
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
