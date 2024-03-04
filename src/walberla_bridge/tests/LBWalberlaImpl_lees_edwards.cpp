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

#include "../src/lattice_boltzmann/LBWalberlaImpl.hpp"

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

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
auto constexpr v0 = 0.064;
static Vector3i mpi_shape{};

static double u_expected(double x, double t, double nu, double v_0, double h,
                         int k_max = 100) {
  auto u = x / h - 0.5;
  for (int k = 1; k <= k_max; k++) {
    u += 1.0 / (M_PI * k) * exp(-4 * M_PI * M_PI * nu * k * k / (h * h) * t) *
         sin(2 * M_PI / h * k * x);
  }
  return v_0 * u;
}

BOOST_AUTO_TEST_CASE(test_transient_shear) {
  using LBImplementation = walberla::LBWalberlaImpl<double, lbmpy::Arch::CPU>;
  double density = 1;
  double viscosity = 1. / 7.;
  auto lattice =
      std::make_shared<LatticeWalberla>(Vector3i{8, 64, 8}, mpi_shape, 1);
  auto lb = LBImplementation(lattice, viscosity, density);
  auto le_pack = std::make_unique<LeesEdwardsPack>(
      0u, 1u, []() { return 0.0; }, [=]() { return v0; });
  lb.set_collision_model(std::move(le_pack));
  lb.ghost_communication();
  auto const grid_size_y = lattice->get_grid_dimensions()[1];
  for (int i = 0; i < 200; i++) {
    lb.integrate();
    if (i < grid_size_y / 2.)
      continue;
    for (double y :
         {0., 0.13 * grid_size_y, 0.7 * grid_size_y, 1. * grid_size_y}) {
      auto u = lb.get_velocity_at_pos(Vector3d{4, y, 4}, true);
      auto expected = u_expected(y, i, viscosity, v0, grid_size_y);
      BOOST_CHECK_SMALL((*u)[0] - expected, 3E-5);
    }
  }
}

static auto setup_lb_with_offset(double offset, double vel_shift = 0.) {
  using LBImplementation = walberla::LBWalberlaImpl<double, lbmpy::Arch::CPU>;
  auto density = 1.;
  auto viscosity = 1. / 7.;
  auto lattice =
      std::make_shared<LatticeWalberla>(Vector3i{10, 10, 10}, mpi_shape, 1);
  auto lb = std::make_shared<LBImplementation>(lattice, viscosity, density);
  auto le_pack = std::make_unique<LeesEdwardsPack>(
      0u, 1u, [=]() { return offset; }, [=]() { return vel_shift; });
  lb->set_collision_model(std::move(le_pack));
  lb->ghost_communication();
  return lb;
}

BOOST_AUTO_TEST_CASE(test_interpolation_force) {
  auto const offset = 2;
  auto lb = setup_lb_with_offset(offset);
  auto const shape = lb->get_lattice().get_grid_dimensions();
  auto const xz = shape[0] / 2;
  auto const y_max = shape[1] - 1;

  auto const force_pos = Vector3d{xz + 0.5, y_max + 0.5, xz + 0.5};
  auto const force_node = Vector3i{xz, y_max, xz};
  auto const f1 = Vector3d{0.3, -0.2, 0.3};
  lb->add_force_at_pos(force_pos, f1);

  lb->integrate();

  auto const ghost_node = Vector3i{force_node[0] - offset, -1, force_node[2]};
  auto const laf = *(lb->get_node_last_applied_force(ghost_node, true));
  BOOST_CHECK_SMALL((laf - f1).norm(), 1E-10);
}

int fold(int i, int size) {
  while (i < 0)
    i += size;
  while (i >= size)
    i -= size;
  return i;
}

std::vector<Utils::Vector3i>
get_mirroring_xz_ghosts(const Utils::Vector3i &n, const Utils::Vector3i shape) {
  std::vector<Utils::Vector3i> res;

  auto in_ghost_layer = [](int i, int shape) {
    return (i == -1 or i == shape);
  };
  auto in_box_or_ghost_layer = [](int i, int shape) {
    return (i >= -1 and i <= shape);
  };

  for (int dx : {-shape[0], 0, shape[0]}) {
    for (int dz : {-shape[2], 0, shape[2]}) {
      if (dx == 0 and dz == 0) // no shift
        continue;

      auto shifted = Utils::Vector3i{n[0] + dx, n[1], n[2] + dz};
      if ((in_ghost_layer(shifted[0], shape[0]) and
           in_box_or_ghost_layer(shifted[2], shape[2])) or
          (in_ghost_layer(shifted[2], shape[2]) and
           in_box_or_ghost_layer(shifted[0], shape[0]))) {
        res.push_back(shifted);
      }
    }
  }
  return res;
}
template <typename LB>
void check_mirroring_ghost_vel(const LB &lb, Vector3i orig_ghost_node) {
  auto const ref_vel = *(lb->get_node_velocity(orig_ghost_node, true));
  auto const shape = lb->get_lattice().get_grid_dimensions();
  //    std::cerr << "mirrors of "<<orig_ghost_node<<" v="<<ref_vel<<" lb
  //    shape="<<shape<<std::endl;
  for (auto const mirror_node :
       get_mirroring_xz_ghosts(orig_ghost_node, shape)) {
    auto const mirror_vel = *(lb->get_node_velocity(mirror_node, true));
    //    std::cerr << "node "<<mirror_node<<" v="<<mirror_vel<<std::endl;
    BOOST_CHECK_SMALL((mirror_vel - ref_vel).norm(), 1E-10);
  }
}

BOOST_AUTO_TEST_CASE(test_interpolation_velocity_int_offset) {
  auto const offset = 12;
  auto vel_shift = 1.5;
  Vector3d vel_shift_vec{vel_shift, 0., 0.};
  auto lb = setup_lb_with_offset(offset, vel_shift);
  auto const shape = lb->get_lattice().get_grid_dimensions();
  for (int x = 0; x < shape[0]; x++) {
    for (int z : {0, 3, shape[2] - 1}) {
      auto const y_max = shape[1] - 1;
      auto const v_upper = Vector3d{0.3 + x, -0.2 + y_max, 0.3 + z};
      auto const v_lower = Vector3d{0.5 + x, -0.7, 0.9 + z};

      auto const upper_source_node = Vector3i{x, y_max, z};
      auto const lower_source_node = Vector3i{x, 0, z};
      auto const ghost_of_upper_node =
          Vector3i{fold(upper_source_node[0] - offset, shape[0]), -1,
                   upper_source_node[2]};
      auto const ghost_of_lower_node =
          Vector3i{fold(lower_source_node[0] + offset, shape[0]), y_max + 1,
                   lower_source_node[2]};

      lb->set_node_velocity(upper_source_node, v_upper);
      lb->set_node_velocity(lower_source_node, v_lower);

      lb->ghost_communication();

      auto const ghost_vel_for_upper =
          *(lb->get_node_velocity(ghost_of_upper_node, true));
      auto const ghost_vel_for_lower =
          *(lb->get_node_velocity(ghost_of_lower_node, true));

      BOOST_CHECK_SMALL((ghost_vel_for_upper + vel_shift_vec - v_upper).norm(),
                        1E-10);
      BOOST_CHECK_SMALL((ghost_vel_for_lower - vel_shift_vec - v_lower).norm(),
                        1E-10);
      check_mirroring_ghost_vel(lb, ghost_of_upper_node);
      check_mirroring_ghost_vel(lb, ghost_of_lower_node);
    }
  }
}


Vector3d pos_from_node(Vector3i n) {
  return {0.5+n[0],0.5+n[1],0.5+n[2]};
}

BOOST_AUTO_TEST_CASE(test_interpolation_velocity_non_int_offset) {
  auto const offset = .5;
  auto vel_shift = 0;
  Vector3d vel_shift_vec{vel_shift, 0., 0.};
  auto lb = setup_lb_with_offset(offset, vel_shift);
  auto const shape = lb->get_lattice().get_grid_dimensions();
  for (int x = 0; x < shape[0]; x++) {
    for (int z : {0, 3, shape[2] - 1}) {
      auto const y_max = shape[1] - 1;
      auto const v_upper = Vector3d{x, -0.2 + y_max, 0.3 + z};
      auto const v_lower = Vector3d{x, -0.7, 0.9 + z};

      auto const upper_source_node = Vector3i{x, y_max, z};
      auto const lower_source_node = Vector3i{x, 0, z};
      lb->set_node_velocity(upper_source_node, v_upper);
      lb->set_node_velocity(lower_source_node, v_lower);
      auto upper_source_pos = pos_from_node(upper_source_node);
      auto lower_source_pos = pos_from_node(lower_source_node);
      BOOST_CHECK_SMALL((*(lb->get_velocity_at_pos(lower_source_pos)) - v_lower).norm(), 1E-8);
      BOOST_CHECK_SMALL((*(lb->get_velocity_at_pos(upper_source_pos)) - v_upper).norm(), 1E-8);

    }
  }
  
  lb->ghost_communication();
  
  for (int y: {0,10,9,-1}) {
    std::cout << "y: "<<y<<"| ";
    for (int x = -1; x <= shape[0]; x++) {
       Vector3i node{x,y,0};
       std::cout << (*(lb->get_node_velocity(node,true)))[0] << " ";
    }
    std::cout << std::endl;
  }
  return;
  for (int x = 0; x < shape[0]; x++) {
    int z=0;
    Vector3i node = {x,-1,z};
    auto vel = *(lb->get_node_velocity(node,true));
    BOOST_CHECK_SMALL((*(lb->get_velocity_at_pos(pos_from_node(node),true)) - vel).norm(), 1E-8);
    check_mirroring_ghost_vel(lb,node);

    node = {x,10,z};
    vel = *(lb->get_node_velocity(node,true));
    BOOST_CHECK_SMALL((*(lb->get_velocity_at_pos(pos_from_node(node), true)) - vel).norm(), 1E-8);
    check_mirroring_ghost_vel(lb,node);
  }
}

BOOST_AUTO_TEST_CASE(test_interpolation_pdf) {
  auto const offset = 2;
  auto lb = setup_lb_with_offset(offset);
  auto const shape = lb->get_lattice().get_grid_dimensions();
  auto const xz = shape[0] / 2;
  auto const y_max = shape[1] - 1;

  auto const source_node = Vector3i{xz, y_max, xz};

  std::vector<double> source_pop(19);
  auto x = -1.;
  std::for_each(source_pop.begin(), source_pop.end(), [&x](auto &v) {
    v = x;
    x += .1;
  });
  lb->set_node_population(source_node, source_pop);
  lb->ghost_communication();

  auto const ghost_node = Vector3i{source_node[0] - offset, -1, source_node[2]};
  auto const ghost_pop = *(lb->get_node_population(ghost_node, true));
  for (unsigned int i = 0u; i < source_pop.size(); ++i) {
    BOOST_CHECK_EQUAL(source_pop[i], ghost_pop[i]);
  }
}

int main(int argc, char **argv) {
  int n_nodes;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());
  walberla::mpi_init();

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // WALBERLA
int main(int argc, char **argv) {}
#endif
