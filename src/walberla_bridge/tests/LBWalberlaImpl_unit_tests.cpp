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
#include <LBWalberlaD3Q19FluctuatingMRT.hpp>
#include <LBWalberlaD3Q19MRT.hpp>
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

using walberla::LBWalberlaD3Q19FluctuatingMRT;
using walberla::LBWalberlaD3Q19MRT;

constexpr double viscosity = 0.4;
constexpr Vector3d box_dimensions = {6, 6, 9};
constexpr double agrid = 0.5;
constexpr Vector3i grid_dimensions{int(box_dimensions[0] / agrid),
                                   int(box_dimensions[1] / agrid),
                                   int(box_dimensions[2] / agrid)};
constexpr double tau = 0.34;
constexpr double density = 2.5;
Vector3i mpi_shape;
constexpr unsigned int seed = 3;
constexpr double kT = 0.0014;

namespace bdata = boost::unit_test::data;

using LbGeneratorVector =
    std::vector<std::function<std::shared_ptr<LBWalberlaBase>(Vector3i)>>;

// Add all LBs with kT=0 to be tested, here
LbGeneratorVector unthermalized_lbs() {
  LbGeneratorVector lbs;
  // Unthermalized D3Q19 MRT
  lbs.push_back([](Vector3i mpi_shape) {
    return std::make_shared<LBWalberlaD3Q19MRT>(viscosity, density, agrid, tau,
                                                box_dimensions, mpi_shape, 1);
  });

  // Thermalized D3Q19 MRT with kT set to 0
  lbs.push_back([](Vector3i mpi_shape) {
    return std::make_shared<LBWalberlaD3Q19FluctuatingMRT>(
        viscosity, density, agrid, tau, box_dimensions, mpi_shape, 1, 0.0,
        seed);
  });
  return lbs;
}

// Add all LBs with thermalization to be tested, here
LbGeneratorVector thermalized_lbs() {
  LbGeneratorVector lbs;

  // Thermalized D3Q19 MRT with kT set to 0
  lbs.push_back([](Vector3i mpi_shape) {
    return std::make_shared<LBWalberlaD3Q19FluctuatingMRT>(
        viscosity, density, agrid, tau, box_dimensions, mpi_shape, 1, kT, seed);
  });
  return lbs;
}

LbGeneratorVector all_lbs() {
  auto lbs = unthermalized_lbs();
  auto thermalized = thermalized_lbs();
  lbs.insert(lbs.end(), thermalized.begin(), thermalized.end());
  return lbs;
}

// Disable printing of type which does not support it
BOOST_TEST_DONT_PRINT_LOG_VALUE(LbGeneratorVector::value_type)

BOOST_DATA_TEST_CASE(dimensions, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape);
  BOOST_CHECK(lb->get_grid_dimensions() == grid_dimensions);

  auto local_domain = lb->get_local_domain();
  auto my_size = local_domain.second - local_domain.first;
  for (int i : {0, 1, 2}) {
    BOOST_CHECK(my_size[i] > 0);
    BOOST_CHECK(local_domain.first[i] >= 0);
    BOOST_CHECK(local_domain.second[i] < grid_dimensions[i]);
  }

  Vector3d total_size;
  MPI_Allreduce(my_size.data(), total_size.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  for (int i : {0, 1, 2}) {
    BOOST_CHECK(total_size[i] == lb->get_grid_dimensions()[i]);
  }
}

BOOST_DATA_TEST_CASE(set_viscosity, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape);
  double new_viscosity = 2.0;
  lb->set_viscosity(new_viscosity);
  BOOST_CHECK_CLOSE(lb->get_viscosity(), new_viscosity, 1E-11);
}

std::vector<Vector3i> all_nodes_incl_ghosts(int n_ghost_layers) {
  std::vector<Vector3i> res;
  for (int x = -n_ghost_layers; x < grid_dimensions[0] + n_ghost_layers; x++) {
    for (int y = -n_ghost_layers; y < grid_dimensions[1] + n_ghost_layers;
         y++) {
      for (int z = -n_ghost_layers; z < grid_dimensions[2] + n_ghost_layers;
           z++) {
        res.push_back(Vector3i{x, y, z});
      }
    }
  }
  return res;
}

std::vector<Vector3i>
local_nodes_incl_ghosts(std::pair<Vector3d, Vector3d> local_domain,
                        int n_ghost_layers) {
  std::vector<Vector3i> res;
  auto const left = local_domain.first;
  auto const right = local_domain.second;
  for (int x = left[0] - n_ghost_layers; x < right[0] + n_ghost_layers; x++) {
    for (int y = left[1] - n_ghost_layers; y < right[1] + n_ghost_layers; y++) {
      for (int z = left[2] - n_ghost_layers; z < right[2] + n_ghost_layers;
           z++) {
        res.push_back(Vector3i{x, y, z});
      }
    }
  }
  return res;
}

BOOST_DATA_TEST_CASE(initial_state, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape);
  for (auto &node :
       local_nodes_incl_ghosts(lb->get_local_domain(), lb->n_ghost_layers())) {
    bool consider_ghosts = !lb->node_in_local_domain(node);
    BOOST_CHECK(!(*lb->get_node_is_boundary(node, consider_ghosts)));
    if (lb->node_in_local_domain(node)) {
      BOOST_CHECK((*lb->get_node_force_to_be_applied(node)) == Vector3d{});
      BOOST_CHECK((*lb->get_node_last_applied_force(node)) == Vector3d{});
      BOOST_CHECK_CLOSE((*lb->get_node_density(node)), density, 1E-10);
    }
    // Todo: add initial pressure tensor check
  }

  BOOST_CHECK(lb->get_momentum() == Vector3d{});
  BOOST_CHECK_CLOSE(lb->get_viscosity(), viscosity, 1E-11);
}

BOOST_DATA_TEST_CASE(kT_unthermalized, bdata::make(unthermalized_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape);
  BOOST_CHECK_EQUAL(lb->get_kT(), 0.0);
}

BOOST_DATA_TEST_CASE(kT_thermalized, bdata::make(thermalized_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape);
  BOOST_CHECK_EQUAL(lb->get_kT(), kT);
}

BOOST_DATA_TEST_CASE(boundary, bdata::make(all_lbs()), lb_generator) {
  Vector3d vel = {0.2, 3.8, 4.2};
  auto lb = lb_generator(mpi_shape);
  for (Vector3i node :

       std::vector<Vector3i>{
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
  auto lb = lb_generator(mpi_shape);
  int n_ghost_layers = 1;

  auto my_left = lb->get_local_domain().first;
  auto my_right = lb->get_local_domain().second;

  for (auto const &n : all_nodes_incl_ghosts(n_ghost_layers)) {
    const Vector3d pos =
        Vector3d{{double(n[0] + .5), double(n[1] + .5), double(n[2] + .5)}};
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
    if (n[0] >= 0 and n[1] >= 0 and n[2] >= 0 and n[0] < grid_dimensions[0] and
        n[1] < grid_dimensions[1] and n[2] < grid_dimensions[2]) {
      MPI_Allreduce(MPI_IN_PLACE, &is_local, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      BOOST_CHECK(is_local == 1);
    }
  }
}

BOOST_DATA_TEST_CASE(velocity_at_node_and_pos, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape);

  // Values
  auto n_pos = [](Vector3i node) {
    return Vector3d{
        {double(node[0] + .5), double(node[1] + .5), double(node[2] + .5)}};
  };

  auto fold = [&](Vector3i n) {
    for (int i = 0; i < 3; i++) {
      if (n[i] < 0)
        n[i] += grid_dimensions[i];
      else if (n[i] >= grid_dimensions[i])
        n[i] -= grid_dimensions[i];
    }
    return n;
  };

  auto n_vel = [&fold](Vector3i node) {
    return fold(node) + Vector3d{{1, 2, -.5}};
  };

  // Assign velocities
  int n_ghost_layers = 1;
  for (auto const &node : all_nodes_incl_ghosts(n_ghost_layers)) {
    if (lb->node_in_local_domain(node)) {
      BOOST_CHECK(lb->set_node_velocity(node, n_vel(node)));
    } else {
      // Check that access to node velocity is not possible
      BOOST_CHECK(!lb->set_node_velocity(node, Vector3d{}));
    }
  }

  lb->ghost_communication();

  // check velocities
  for (auto const &node : all_nodes_incl_ghosts(n_ghost_layers)) {
    double eps = 1E-8;

    if (lb->node_in_local_halo(node)) {
      auto res = lb->get_node_velocity(node);
      BOOST_CHECK(res);                               // value available
      BOOST_CHECK((*res - n_vel(node)).norm() < eps); // correct value?
      // Check that the interpolated velocity at the node pos equals the node
      // vel
      res = lb->get_velocity_at_pos(n_pos(node), true);
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
  auto lb = lb_generator(mpi_shape);
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
  auto mom_exp = density * (v1 + v2);
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
}

BOOST_DATA_TEST_CASE(forces_interpolation, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape);

  // todo: check a less symmetrical situation, wher the force is applied not in
  // the middle between the ndoes

  for (Vector3i n : all_nodes_incl_ghosts(1)) {
    if (lb->node_in_local_halo(n)) {
      Vector3d pos{double(n[0]), double(n[1]),
                   double(n[2])}; // Mid point between nodes
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
  auto lb = lb_generator(mpi_shape);

  // Forces added go to force_to_be_applied. After integration, they should be
  // in last_applied_force, where they are used for velocity calculation

  Vector3i origin{};
  Vector3i middle{
      {grid_dimensions[0] / 2, grid_dimensions[1] / 2, grid_dimensions[2] / 2}};
  Vector3i right = grid_dimensions - Vector3i{{1, 1, 1}};

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
    for (auto cn : {n, n + grid_dimensions, n - grid_dimensions,
                    n + Vector3i{{grid_dimensions[0], 0, 0}}}) {
      if (lb->node_in_local_halo(cn)) {
        BOOST_CHECK_SMALL(
            (*(lb->get_node_last_applied_force(cn, true)) - f).norm(), 1E-10);
        BOOST_CHECK_SMALL((*(lb->get_node_force_to_be_applied(cn))).norm(),
                          1E-10);
      }
    }
    lb->integrate();
    for (auto cn : {n, n + grid_dimensions, n - grid_dimensions,
                    n + Vector3i{{grid_dimensions[0], 0, 0}}}) {
      if (lb->node_in_local_halo(cn)) {
        BOOST_CHECK_SMALL((*(lb->get_node_last_applied_force(cn, true))).norm(),
                          1E-10);
        BOOST_CHECK_SMALL((*(lb->get_node_force_to_be_applied(cn))).norm(),
                          1E-10);
      }
    }
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int n_nodes;

  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());

  walberla_mpi_init();
  auto res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

BOOST_DATA_TEST_CASE(velocity_fluctuation, bdata::make(thermalized_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape);

  // Warmup
  for (int i = 0; i < 100; i++)
    lb->integrate();

  // Sample
  int steps = 3000;

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
          if (lb->node_in_local_domain(node)) {
            auto v = *(lb->get_node_velocity(node));
            auto rho = *(lb->get_node_density(node));
            step_v += v * rho;
            step_v_square += rho * hadamard_product(v, v);
          }
        }
      }
    }
    step_v /= double(count);
    step_v_square /= double(count);

    sum_v += step_v;
    sum_v_square += step_v_square;

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
  BOOST_CHECK_CLOSE(sum_v_square[0] / steps, kT, tol_kT);
  BOOST_CHECK_CLOSE(sum_v_square[1] / steps, kT, tol_kT);
  BOOST_CHECK_CLOSE(sum_v_square[2] / steps, kT, tol_kT);
}
#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
