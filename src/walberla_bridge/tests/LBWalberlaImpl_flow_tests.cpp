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
#define BOOST_TEST_MODULE Walberla point force test
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

static LBTestParameters params; // populated in main()

BOOST_DATA_TEST_CASE(integrate_with_point_force_thermalized,
                     bdata::make(thermalized_lbs()), lb_generator) {
  auto lb = lb_generator(params);
  boost::mpi::communicator world;

  // Check that momentum stays zero after initial integration
  lb->integrate();
  lb->integrate();
  auto mom_local = lb->get_momentum();
  auto mom = boost::mpi::all_reduce(world, mom_local, std::plus<Vector3d>());
  BOOST_CHECK_SMALL(mom.norm(), 1E-10);

  // Check that momentum changes as expected when applying forces
  // auto f = Vector3d{0.15, 0.25, -0.22};
  // auto f = Vector3d{0.0006, -0.0013, 0.000528};
  auto const f1 = Vector3d{0., 0., 0.};
  auto const f2 = Vector3d{0.1, 0.2, -0.3};
  lb->set_external_force(f1);
  auto const force_node = Vector3i{{1, 1, 1}};
  lb->add_force_at_pos(force_node + Vector3d::broadcast(.5), f2);
  lb->integrate();
  for (auto const &n : all_nodes_incl_ghosts(lb->get_lattice())) {
    if (lb->get_lattice().node_in_local_halo(n)) {
      auto const laf = *(lb->get_node_last_applied_force(n, true));
      if (n == force_node) {
        BOOST_CHECK_SMALL((laf - f1 - f2).norm(), 1E-10);
      } else {
        BOOST_CHECK_SMALL((laf - f1).norm(), 1E-10);
      }
    }
  }
  mom_local = lb->get_momentum();
  mom = boost::mpi::all_reduce(world, mom_local, std::plus<Vector3d>());

  // Expected momentum = momentum added in prev. time step
  // + f/2 from velocity shift due to last applied forces
  auto mom_exp = 1.5 * f1 * Utils::product(params.grid_dimensions) + 1.5 * f2;
  auto d = mom - mom_exp;
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
  std::cout << "thermalized: " << mom << " | " << mom_exp << " | " << d << "\n";

  // check that momentum doesn't drift when no force is applied again
  lb->set_external_force(Vector3d{});
  // The expected moment is just that applied during a single time step
  // No f/2 correction, since no force was applied in last time step
  mom_exp = 1.0 * f1 * Utils::product(params.grid_dimensions) + 1.0 * f2;
  lb->integrate();
  mom_local = lb->get_momentum();
  mom = boost::mpi::all_reduce(world, mom_local, std::plus<Vector3d>());
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
}

// this can be merged with the thermalized test, once that passes
BOOST_DATA_TEST_CASE(integrate_with_point_force_unthermalized,
                     bdata::make(unthermalized_lbs()), lb_generator) {
  auto lb = lb_generator(params);
  boost::mpi::communicator world;

  // Check that momentum stays zero after initial integration
  lb->integrate();
  BOOST_CHECK_SMALL(lb->get_momentum().norm(), 1E-10);

  // Check that momentum changes as expected when applying forces
  // auto f = Vector3d{0.0006, -0.0013, 0.000528};
  auto const f1 = Vector3d{0., 0., 0.};
  auto const f2 = Vector3d{0.095, 0.23, -0.52};
  lb->set_external_force(f1);
  lb->add_force_at_pos(Utils::Vector3d{2, 2, 2}, f2);
  lb->integrate();

  auto mom_local = lb->get_momentum();
  auto mom = boost::mpi::all_reduce(world, mom_local, std::plus<Vector3d>());

  // Expected momentum = momentum added in prev. time step
  // + f/2 from velocity shift due to last applied forces
  auto mom_exp = 1.5 * f1 * Utils::product(params.grid_dimensions) + 1.5 * f2;
  auto d = mom - mom_exp;
  std::cout << mom << " | " << mom_exp << " | " << d << "\n";
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);

  // check that momentum doesn't drift when no force is applied again
  lb->set_external_force(Vector3d{});
  lb->integrate();
  // The expected moment is just that applied during a single time step
  // No f/2 correction, since no force was applied in last time step
  mom_exp = 1.0 * f1 * Utils::product(params.grid_dimensions) + 1.0 * f2;
  mom_local = lb->get_momentum();
  mom = boost::mpi::all_reduce(world, mom_local, std::plus<Vector3d>());
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
}

int main(int argc, char **argv) {
  int n_nodes;
  Vector3i mpi_shape{};

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());
  walberla::mpi_init();

  params.seed = 0u;
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
