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

#include <mpi.h>

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

// Disable printing of type which does not support it

BOOST_DATA_TEST_CASE(dimensions, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  BOOST_CHECK(lb->get_grid_dimensions() == params.grid_dimensions);

  auto local_domain = lb->get_local_domain();
  auto my_size = local_domain.second - local_domain.first;
  for (int i : {0, 1, 2}) {
    BOOST_CHECK(my_size[i] > 0);
    BOOST_CHECK(local_domain.first[i] >= 0);
    BOOST_CHECK(local_domain.second[i] <= params.grid_dimensions[i]);
  }
}

BOOST_DATA_TEST_CASE(set_viscosity, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  double new_viscosity = 2.0;
  lb->set_viscosity(new_viscosity);
  BOOST_CHECK_CLOSE(lb->get_viscosity(), new_viscosity, 1E-11);
}

BOOST_DATA_TEST_CASE(initial_state, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  for (auto &node :
       local_nodes_incl_ghosts(lb->get_local_domain(), lb->n_ghost_layers())) {
    bool consider_ghosts = !lb->node_in_local_domain(node);
    BOOST_CHECK(!(*lb->get_node_is_boundary(node, consider_ghosts)));
    if (lb->node_in_local_domain(node)) {
      BOOST_CHECK((*lb->get_node_force_to_be_applied(node)) == Vector3d{});
      BOOST_CHECK((*lb->get_node_last_applied_force(node)) == Vector3d{});
      BOOST_CHECK_CLOSE((*lb->get_node_density(node)), params.density, 1E-10);
    }
    // Todo: add initial pressure tensor check
  }

  BOOST_CHECK(lb->get_momentum() == Vector3d{});
  BOOST_CHECK_CLOSE(lb->get_viscosity(), params.viscosity, 1E-11);
}

BOOST_DATA_TEST_CASE(kT_unthermalized, bdata::make(unthermalized_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  BOOST_CHECK_EQUAL(lb->get_kT(), 0.0);
}

BOOST_DATA_TEST_CASE(kT_thermalized, bdata::make(thermalized_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  BOOST_CHECK_EQUAL(lb->get_kT(), params.kT);
}

BOOST_DATA_TEST_CASE(boundary, bdata::make(all_lbs()), lb_generator) {
  Vector3d vel = {0.2, 3.8, 4.2};
  auto lb = lb_generator(mpi_shape, params);
  for (auto const &node : std::vector<Vector3i>{
           {-lb->n_ghost_layers(), 0, 0}, {0, 0, 0}, {0, 1, 2}, {9, 9, 9}}) {
    if (lb->node_in_local_halo(node)) {
      BOOST_CHECK(lb->set_node_velocity_at_boundary(node, vel));
      auto vel_check = lb->get_node_velocity_at_boundary(node);
      // Do we have a value
      BOOST_CHECK(vel_check);
      // Check the value
      BOOST_CHECK_SMALL((*vel_check - vel).norm(), 1E-12);
      BOOST_CHECK(lb->remove_node_from_boundary(node));
      auto res = lb->get_node_is_boundary(node, true);
      // Did we get a value?
      BOOST_CHECK(res);
      // Should not be a boundary node
      BOOST_CHECK(*res == false);
    } else {
      // Not in the local halo.
      BOOST_CHECK(!lb->set_node_velocity_at_boundary(node, vel));
      BOOST_CHECK(!lb->get_node_velocity_at_boundary(node));
      BOOST_CHECK(!lb->remove_node_from_boundary(node));
      BOOST_CHECK(!lb->get_node_is_boundary(node));
    }
  }

  lb->clear_boundaries();
  for (auto &node :
       local_nodes_incl_ghosts(lb->get_local_domain(), lb->n_ghost_layers())) {
    BOOST_CHECK(!(*lb->get_node_is_boundary(node, true)));
  }
}

BOOST_DATA_TEST_CASE(domain_and_halo, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  auto const n_ghost_layers = lb->n_ghost_layers();

  auto my_left = lb->get_local_domain().first;
  auto my_right = lb->get_local_domain().second;

  for (auto const &n :
       all_nodes_incl_ghosts(params.grid_dimensions, n_ghost_layers)) {
    const Vector3d pos = n + Vector3d::broadcast(.5);
    int is_local = 0;
    // Nodes in local domain
    if ((n[0] >= my_left[0] and n[1] >= my_left[1] and n[2] >= my_left[2]) and
        (n[0] < my_right[0] and n[1] < my_right[1] and n[2] < my_right[2])) {
      BOOST_CHECK(lb->node_in_local_domain(n));
      BOOST_CHECK(lb->node_in_local_halo(n));

      BOOST_CHECK(lb->pos_in_local_domain(pos));
      BOOST_CHECK(lb->pos_in_local_halo(pos));
      is_local = 1;
    } else {
      // in local halo?
      if ((n[0] + n_ghost_layers >= my_left[0] and
           n[1] + n_ghost_layers >= my_left[1] and
           n[2] + n_ghost_layers >= my_left[2]) and
          (n[0] - n_ghost_layers < my_right[0] and
           n[1] - n_ghost_layers < my_right[1] and
           n[2] - n_ghost_layers < my_right[2])) {
        BOOST_CHECK(!lb->node_in_local_domain(n));
        BOOST_CHECK(lb->node_in_local_halo(n));

        BOOST_CHECK(!lb->pos_in_local_domain(pos));
        BOOST_CHECK(lb->pos_in_local_halo(pos));
      } else {
        // neither in domain nor in halo
        BOOST_CHECK(!lb->node_in_local_domain(n));
        BOOST_CHECK(!lb->node_in_local_halo(n));

        BOOST_CHECK(!lb->pos_in_local_domain(pos));
        BOOST_CHECK(!lb->pos_in_local_halo(pos));
      }
    }

    // If the cell is in the global physical domain
    // check that only one mpi rank said the node was local
    if (n[0] >= 0 and n[1] >= 0 and n[2] >= 0 and
        n[0] < params.grid_dimensions[0] and
        n[1] < params.grid_dimensions[1] and n[2] < params.grid_dimensions[2]) {
      MPI_Allreduce(MPI_IN_PLACE, &is_local, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      BOOST_CHECK(is_local == 1);
    }
  }
}

BOOST_DATA_TEST_CASE(velocity_at_node_and_pos, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  auto const n_ghost_layers = lb->n_ghost_layers();

  // Values
  auto n_pos = [](Vector3i node) { return node + Vector3d::broadcast(.5); };

  auto fold = [&](Vector3i n) {
    for (int i = 0; i < 3; i++) {
      if (n[i] < 0)
        n[i] += params.grid_dimensions[i];
      else if (n[i] >= params.grid_dimensions[i])
        n[i] -= params.grid_dimensions[i];
    }
    return n;
  };

  auto n_vel = [&fold](Vector3i node) {
    return fold(node) + Vector3d{{1, 2, -.5}};
  };

  // Assign velocities
  for (auto const &node :
       all_nodes_incl_ghosts(params.grid_dimensions, n_ghost_layers)) {
    if (lb->node_in_local_domain(node)) {
      BOOST_CHECK(lb->set_node_velocity(node, n_vel(node)));
    } else {
      // Check that access to node velocity is not possible
      BOOST_CHECK(!lb->set_node_velocity(node, Vector3d{}));
    }
  }

  lb->ghost_communication();

  // check velocities
  for (auto const &node :
       all_nodes_incl_ghosts(params.grid_dimensions, n_ghost_layers)) {
    double eps = 1E-8;

    if (lb->node_in_local_halo(node)) {
      bool consider_ghosts = !lb->node_in_local_domain(node);
      auto res = lb->get_node_velocity(node, consider_ghosts);
      BOOST_CHECK(res);                               // value available
      BOOST_CHECK((*res - n_vel(node)).norm() < eps); // correct value?
      // Check that the interpolated velocity at the node pos equals the node
      // vel
      res = lb->get_velocity_at_pos(n_pos(node), consider_ghosts);
      BOOST_CHECK(res);                                    // locally available
      BOOST_CHECK_SMALL((*res - n_vel(node)).norm(), eps); // value correct?
    } else {
      // Check that access to node velocity is not possible
      BOOST_CHECK(!lb->get_node_velocity(node));
      BOOST_CHECK(!lb->get_velocity_at_pos(n_pos(node), true));
    }
  }
}

BOOST_DATA_TEST_CASE(total_momentum, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  Vector3i n1{1, 2, 3};
  Vector3i n2{9, 2, 10};
  auto v1 = Vector3d{1.5, 2.5, -2.2};
  auto v2 = Vector3d{-.5, 3.5, -.2};
  if (lb->node_in_local_domain(n1)) {
    lb->set_node_velocity(n1, v1);
  };
  if (lb->node_in_local_domain(n2)) {
    lb->set_node_velocity(n2, v2);
  };

  auto mom = lb->get_momentum();
  auto mom_exp = params.density * (v1 + v2);
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
}

BOOST_DATA_TEST_CASE(forces_interpolation, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);

  // todo: check a less symmetrical situation, where the force is applied not
  // in the middle between the nodes

  for (Vector3i n : all_nodes_incl_ghosts(params.grid_dimensions, 1)) {
    if (lb->node_in_local_halo(n)) {
      Vector3d pos = 1. * n; // Mid point between nodes
      Vector3d f = {1, 2, -3.5};
      lb->add_force_at_pos(pos, f);
      // Check neighboring nodes for force to be applied
      for (int x : {0, 1})
        for (int y : {0, 1})
          for (int z : {0, 1}) {
            Vector3i check_node{{n[0] - x, n[1] - y, n[2] - z}};
            if (lb->node_in_local_halo(check_node)) {
              auto res = lb->get_node_force_to_be_applied(check_node);
              BOOST_CHECK_SMALL(((*res) - f / 8.0).norm(), 1E-10);
            }
          }
      // Apply counter force to clear force field
      lb->add_force_at_pos(pos, -f);
    }
  }
}

BOOST_DATA_TEST_CASE(forces_book_keeping, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);

  // Forces added go to force_to_be_applied. After integration, they should be
  // in last_applied_force, where they are used for velocity calculation

  Vector3i origin{};
  Vector3i middle = params.grid_dimensions / 2;
  Vector3i right = params.grid_dimensions - Vector3i{{1, 1, 1}};

  Vector3d f{{1, -2, 3.1}};

  for (auto n : {origin, middle, right}) {
    // Add force to node position
    if (lb->node_in_local_domain(n)) {
      lb->add_force_at_pos(Vector3d{{n[0] + .5, n[1] + .5, n[2] + .5}}, f);
      BOOST_CHECK_SMALL((*(lb->get_node_force_to_be_applied(n)) - f).norm(),
                        1E-10);
    }
    lb->integrate();
    // Check nodes incl some of the ghosts
    for (auto cn : {n, n + params.grid_dimensions, n - params.grid_dimensions,
                    n + Vector3i{{params.grid_dimensions[0], 0, 0}}}) {
      if (lb->node_in_local_halo(cn)) {
        BOOST_CHECK_SMALL(
            (*(lb->get_node_last_applied_force(cn, true)) - f).norm(), 1E-10);
        BOOST_CHECK_SMALL((*(lb->get_node_force_to_be_applied(cn))).norm(),
                          1E-10);
      }
    }
    lb->integrate();
    for (auto cn : {n, n + params.grid_dimensions, n - params.grid_dimensions,
                    n + Vector3i{{params.grid_dimensions[0], 0, 0}}}) {
      if (lb->node_in_local_halo(cn)) {
        BOOST_CHECK_SMALL((*(lb->get_node_last_applied_force(cn, true))).norm(),
                          1E-10);
        BOOST_CHECK_SMALL((*(lb->get_node_force_to_be_applied(cn))).norm(),
                          1E-10);
      }
    }
  }
}

BOOST_DATA_TEST_CASE(force_in_corner, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);

  // Add forces in all box corners. If domain boundaries are treated correclty
  // each corner node should get 1/8 of the force.

  auto const l = params.box_dimensions;
  const Vector3d f{0.1, .02, -0.3};
  for (double x : {0., l[0]})
    for (double y : {0., l[1]})
      for (double z : {0., l[2]}) {
        const Vector3d pos{x, y, z};
        bool res = lb->add_force_at_pos(pos, f);
      }

  // check forces to be applied
  // Each corner node should have 1/8 of the force
  const double tol = 1E-10;
  int count = 0;
  for (auto c : corner_nodes(params.grid_dimensions)) {
    auto res = lb->get_node_force_to_be_applied(c);
    if (res) {
      BOOST_CHECK_SMALL(((*res) - f / 8.0).norm(), tol);
      count += 1;
    }
  };
  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  BOOST_CHECK_EQUAL(count, 8);

  lb->integrate();

  // check applied forces from last integration step
  count = 0;
  for (auto c : corner_nodes(params.grid_dimensions)) {
    auto res = lb->get_node_last_applied_force(c);
    if (res) {
      BOOST_CHECK_SMALL(((*res) - f / 8.0).norm(), tol);
      count += 1;
    }
  };
  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  BOOST_CHECK_EQUAL(count, 8);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int n_nodes;

  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());

  params.viscosity = 0.003;
  params.kT = 1.3E-4;
  params.density = 1.4;
  params.grid_dimensions = Vector3i{12, 12, 18};
  params.box_dimensions = Vector3d{12, 12, 18};

  walberla_mpi_init();
  auto res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
