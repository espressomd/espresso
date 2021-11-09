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

#include "tests_common.hpp"

#include <LBWalberlaBase.hpp>
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/inplace.hpp>
#include <boost/multi_array.hpp>

#include <mpi.h>

#include <functional>
#include <type_traits>
#include <vector>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

LBTestParameters params; // populated in main()
Vector3i mpi_shape;      // populated in main

BOOST_DATA_TEST_CASE(dimensions, bdata::make(all_lbs()), lb_generator) {
  using boost::test_tools::per_element;
  auto lb = lb_generator(mpi_shape, params);
  auto constexpr zero = Vector3i{0, 0, 0};

  auto const grid_dim = lb->get_grid_dimensions();
  BOOST_TEST(grid_dim == params.grid_dimensions, per_element());

  auto const [my_left, my_right] = lb->get_local_domain();
  auto const my_size = my_right - my_left;
  BOOST_TEST(my_size > zero, per_element());
  BOOST_TEST(my_left >= zero, per_element());
  BOOST_TEST(my_right <= params.grid_dimensions, per_element());
}

BOOST_DATA_TEST_CASE(set_viscosity, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  double new_viscosity = 2.0;
  lb->set_viscosity(new_viscosity);
  BOOST_CHECK_CLOSE(lb->get_viscosity(), new_viscosity, 1E-11);
}

BOOST_DATA_TEST_CASE(initial_state, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  for (auto const &node :
       local_nodes_incl_ghosts(lb->get_local_domain(), lb->n_ghost_layers())) {
    bool const consider_ghosts = !lb->node_in_local_domain(node);
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

BOOST_DATA_TEST_CASE(per_node_boundary, bdata::make(all_lbs()), lb_generator) {
  auto const vel = Vector3d{{0.2, 3.8, 4.2}};
  auto lb = lb_generator(mpi_shape, params);
  for (auto const &node : std::vector<Vector3i>{
           {-lb->n_ghost_layers(), 0, 0}, {0, 0, 0}, {0, 1, 2}, {9, 9, 9}}) {
    if (lb->node_in_local_halo(node)) {
      {
        auto const res = lb->get_node_is_boundary(node, true);
        // Did we get a value?
        BOOST_REQUIRE(res);
        // Should not be a boundary node
        BOOST_CHECK(*res == false);
      }
      {
        BOOST_CHECK(lb->set_node_velocity_at_boundary(node, vel, true));
        auto const res = lb->get_node_is_boundary(node, true);
        // Did we get a value?
        BOOST_REQUIRE(res);
        // Should be a boundary node
        BOOST_CHECK(*res == true);
      }
      {
        auto const vel_check = lb->get_node_velocity_at_boundary(node);
        // Do we have a value
        BOOST_REQUIRE(vel_check);
        // Check the value
        BOOST_CHECK_SMALL((*vel_check - vel).norm(), 1E-12);
      }
      {
        BOOST_CHECK(lb->remove_node_from_boundary(node, true));
        auto const res = lb->get_node_is_boundary(node, true);
        // Did we get a value?
        BOOST_REQUIRE(res);
        // Should not be a boundary node
        BOOST_CHECK(*res == false);
      }
    } else {
      // Not in the local halo.
      BOOST_CHECK(!lb->set_node_velocity_at_boundary(node, vel, true));
      BOOST_CHECK(!lb->get_node_velocity_at_boundary(node));
      BOOST_CHECK(!lb->remove_node_from_boundary(node, true));
      BOOST_CHECK(!lb->get_node_is_boundary(node));
    }
  }

  lb->clear_boundaries();
  for (auto const &node :
       local_nodes_incl_ghosts(lb->get_local_domain(), lb->n_ghost_layers())) {
    BOOST_CHECK(!(*lb->get_node_is_boundary(node, true)));
  }
}

BOOST_DATA_TEST_CASE(update_boundary_from_shape, bdata::make(all_lbs()),
                     lb_generator) {
  auto const vel = Vector3d{{0.2, 3.8, 4.2}};
  auto lb = lb_generator(mpi_shape, params);

  auto const vec3to4 = [](Utils::Vector<int, 3> const &d, int v) {
    return Utils::Vector<int, 4>{{d[0], d[1], d[2], v}};
  };

  auto const nodes = std::vector<Vector3i>{
      {-lb->n_ghost_layers(), 0, 0}, {0, 0, 0}, {0, 1, 2}, {9, 9, 9}};
  // set up boundary
  {
    auto const n_grid_points = Utils::product(params.grid_dimensions);
    boost::multi_array<int, 3> raster_3d(params.grid_dimensions);
    boost::multi_array<double, 4> vel_3d(vec3to4(params.grid_dimensions, 3));
    BOOST_CHECK_EQUAL(raster_3d.num_elements(), n_grid_points);
    for (auto const &node : nodes) {
      auto const idx = (node + params.grid_dimensions) % params.grid_dimensions;
      raster_3d(idx) = 1;
      for (auto const i : {0, 1, 2}) {
        vel_3d(vec3to4(idx, i)) = vel[i];
      }
    }
    std::vector<int> raster_flat(raster_3d.data(),
                                 raster_3d.data() + raster_3d.num_elements());
    std::vector<double> vel_flat(vel_3d.data(),
                                 vel_3d.data() + vel_3d.num_elements());
    lb->update_boundary_from_shape(raster_flat, vel_flat);
  }

  for (auto const &node : nodes) {
    if (lb->node_in_local_halo(node)) {
      {
        auto const res = lb->get_node_is_boundary(node, true);
        // Did we get a value?
        BOOST_REQUIRE(res);
        // Should be a boundary node
        BOOST_CHECK(*res == true);
      }
      {
        auto const vel_check = lb->get_node_velocity_at_boundary(node);
        // Do we have a value
        BOOST_REQUIRE(vel_check);
        // Check the value
        BOOST_CHECK_SMALL((*vel_check - vel).norm(), 1E-12);
      }
    } else {
      // Not in the local halo.
      BOOST_CHECK(!lb->get_node_velocity_at_boundary(node));
    }
  }

  lb->clear_boundaries();
  lb->ghost_communication();
  for (auto const &node :
       local_nodes_incl_ghosts(lb->get_local_domain(), lb->n_ghost_layers())) {
    BOOST_CHECK(!(*lb->get_node_is_boundary(node, true)));
  }
}

BOOST_DATA_TEST_CASE(update_boundary_from_list, bdata::make(all_lbs()),
                     lb_generator) {
  auto const vel = Vector3d{{0.2, 3.8, 4.2}};
  auto lb = lb_generator(mpi_shape, params);

  auto const nodes = std::vector<Vector3i>{
      {-lb->n_ghost_layers(), 0, 0}, {0, 0, 0}, {0, 1, 2}, {9, 9, 9}};
  // set up boundary
  {
    std::vector<double> velocities_flat;
    std::vector<int> indices_flat;
    for (auto const &node : nodes) {
      for (auto const i : {0, 1, 2}) {
        indices_flat.emplace_back(node[i]);
        velocities_flat.emplace_back(vel[i]);
      }
    }
    lb->update_boundary_from_list(indices_flat, velocities_flat);
  }

  for (auto const &node : nodes) {
    if (lb->node_in_local_halo(node)) {
      {
        auto const res = lb->get_node_is_boundary(node, true);
        // Did we get a value?
        BOOST_REQUIRE(res);
        // Should be a boundary node
        BOOST_CHECK(*res == true);
      }
      {
        auto const vel_check = lb->get_node_velocity_at_boundary(node);
        // Do we have a value
        BOOST_REQUIRE(vel_check);
        // Check the value
        BOOST_CHECK_SMALL((*vel_check - vel).norm(), 1E-12);
      }
    } else {
      // Not in the local halo.
      BOOST_CHECK(!lb->get_node_velocity_at_boundary(node));
    }
  }

  lb->clear_boundaries();
  lb->ghost_communication();
  for (auto const &node :
       local_nodes_incl_ghosts(lb->get_local_domain(), lb->n_ghost_layers())) {
    BOOST_CHECK(!(*lb->get_node_is_boundary(node, true)));
  }
}

BOOST_DATA_TEST_CASE(domain_and_halo, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  auto const n_ghost_layers = lb->n_ghost_layers();

  auto const my_left = lb->get_local_domain().first;
  auto const my_right = lb->get_local_domain().second;

  for (auto const &n :
       all_nodes_incl_ghosts(params.grid_dimensions, n_ghost_layers)) {
    auto const pos = n + Vector3d::broadcast(.5);
    int is_local = 0;
    // Nodes in local domain
    if (Vector3d(n) >= my_left and Vector3d(n) < my_right) {
      BOOST_CHECK(lb->node_in_local_domain(n));
      BOOST_CHECK(lb->node_in_local_halo(n));

      BOOST_CHECK(lb->pos_in_local_domain(pos));
      BOOST_CHECK(lb->pos_in_local_halo(pos));
      is_local = 1;
    } else {
      // in local halo?
      if ((n + Vector3d::broadcast(n_ghost_layers)) >= my_left and
          (n - Vector3d::broadcast(n_ghost_layers)) < my_right) {
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
    constexpr auto origin = Vector3i{0, 0, 0};
    if (n >= origin and n < params.grid_dimensions) {
      boost::mpi::communicator world;
      boost::mpi::all_reduce(world, boost::mpi::inplace(is_local),
                             std::plus<int>());
      BOOST_CHECK(is_local == 1);
    }
  }
}

auto fold_node(Vector3i n) {
  for (int i = 0; i < 3; i++) {
    if (n[i] < 0)
      n[i] += params.grid_dimensions[i];
    else if (n[i] >= params.grid_dimensions[i])
      n[i] -= params.grid_dimensions[i];
  }
  return n;
}

BOOST_DATA_TEST_CASE(velocity_at_node_and_pos, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  auto const n_ghost_layers = lb->n_ghost_layers();

  // Values
  auto n_pos = [](Vector3i const &n) { return n + Vector3d::broadcast(.5); };

  auto n_vel = [](Vector3i const &node) {
    return fold_node(node) + Vector3d{{1., 2., -.5}};
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
    auto constexpr eps = 1E-8;
    if (lb->node_in_local_halo(node)) {
      bool const consider_ghosts = !lb->node_in_local_domain(node);
      auto res = lb->get_node_velocity(node, consider_ghosts);
      BOOST_CHECK(res);                                    // value available
      BOOST_CHECK_SMALL((*res - n_vel(node)).norm(), eps); // value correct?
      // Check that the interpolated velocity at the node pos equals the node
      // vel
      res = lb->get_velocity_at_pos(n_pos(node), consider_ghosts);
      BOOST_CHECK(res);                                    // value available
      BOOST_CHECK_SMALL((*res - n_vel(node)).norm(), eps); // value correct?
    } else {
      // Check that access to node velocity is not possible
      BOOST_CHECK(!lb->get_node_velocity(node));
      BOOST_CHECK(!lb->get_velocity_at_pos(n_pos(node), true));
    }
  }
}

BOOST_DATA_TEST_CASE(interpolated_density_at_pos, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  auto const n_ghost_layers = lb->n_ghost_layers();

  // Values
  auto n_pos = [](Vector3i const &n) { return n + Vector3d::broadcast(.5); };

  auto n_dens = [](Vector3i const &node) {
    return 1.0 + static_cast<double>(Utils::product(fold_node(node))) * 1e-6;
  };

  // Assign densities
  for (auto const &node :
       all_nodes_incl_ghosts(params.grid_dimensions, n_ghost_layers)) {
    if (lb->node_in_local_domain(node)) {
      BOOST_CHECK(lb->set_node_density(node, n_dens(node)));
    } else {
      // Check that access to node density is not possible
      BOOST_CHECK(!lb->set_node_density(node, 0.));
    }
  }

  lb->ghost_communication();

  // check densities
  for (auto const &node :
       all_nodes_incl_ghosts(params.grid_dimensions, n_ghost_layers)) {
    auto constexpr eps = 1E-8;
    if (lb->node_in_local_halo(node)) {
      if (lb->node_in_local_domain(node)) {
        auto res = lb->get_node_density(node);
        BOOST_CHECK(res);                            // value available
        BOOST_CHECK_SMALL(*res - n_dens(node), eps); // value correct?
        // Check that the interpolated density at the node pos equals the node
        // density
        res = lb->get_interpolated_density_at_pos(n_pos(node));
        BOOST_CHECK(res);                            // value available
        BOOST_CHECK_SMALL(*res - n_dens(node), eps); // value correct?
      } else {
        BOOST_CHECK(!lb->get_node_density(node));
        BOOST_CHECK(!lb->get_interpolated_density_at_pos(n_pos(node), false));
        if (node == Vector3i{-1, -1, -1}) {
          BOOST_CHECK_THROW(
              lb->get_interpolated_density_at_pos(n_pos(node), true),
              std::runtime_error);
        }
      }
    } else {
      // Check that access to node density is not possible
      BOOST_CHECK(!lb->get_node_density(node));
      BOOST_CHECK(!lb->get_interpolated_density_at_pos(n_pos(node), true));
    }
  }
}

BOOST_DATA_TEST_CASE(total_momentum, bdata::make(all_lbs()), lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  auto const n1 = Vector3i{{1, 2, 3}};
  auto const n2 = Vector3i{{9, 2, 10}};
  auto const v1 = Vector3d{{1.5, 2.5, -2.2}};
  auto const v2 = Vector3d{{-.5, 3.5, -.2}};
  if (lb->node_in_local_domain(n1)) {
    lb->set_node_velocity(n1, v1);
  }
  if (lb->node_in_local_domain(n2)) {
    lb->set_node_velocity(n2, v2);
  }

  boost::mpi::communicator world;
  auto const mom_local = lb->get_momentum();
  auto const mom_exp = params.density * (v1 + v2);
  auto const mom =
      boost::mpi::all_reduce(world, mom_local, std::plus<Vector3d>());
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
}

BOOST_DATA_TEST_CASE(forces_interpolation, bdata::make(all_lbs()),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);

  // todo: check a less symmetrical situation, where the force is applied not
  // in the middle between the nodes

  for (Vector3i n : all_nodes_incl_ghosts(params.grid_dimensions, 1)) {
    if (lb->node_in_local_halo(n)) {
      auto const pos = 1. * n; // Mid point between nodes
      auto const f = Vector3d{{1., 2., -3.5}};
      lb->add_force_at_pos(pos, f);
      // Check neighboring nodes for force to be applied
      for (int x : {0, 1})
        for (int y : {0, 1})
          for (int z : {0, 1}) {
            auto const check_node = Vector3i{{n[0] - x, n[1] - y, n[2] - z}};
            if (lb->node_in_local_halo(check_node)) {
              auto const res = lb->get_node_force_to_be_applied(check_node);
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

  Vector3i const origin{};
  Vector3i const middle = params.grid_dimensions / 2;
  Vector3i const right = params.grid_dimensions - Vector3i{{1, 1, 1}};

  Vector3d const f{{1., -2., 3.1}};

  for (auto n : {origin, middle, right}) {
    // Add force to node position
    if (lb->node_in_local_domain(n)) {
      lb->add_force_at_pos(n + Vector3d::broadcast(.5), f);
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
  boost::mpi::communicator world;

  // Add forces in all box corners. If domain boundaries are treated correctly
  // each corner node should get 1/8 of the force.

  auto const l = params.box_dimensions;
  auto const f = Vector3d{{0.1, .02, -0.3}};
  for (double x : {0., l[0]}) {
    for (double y : {0., l[1]}) {
      for (double z : {0., l[2]}) {
        auto const pos = Vector3d{x, y, z};
        static_cast<void>(lb->add_force_at_pos(pos, f));
      }
    }
  }

  // check forces to be applied
  // Each corner node should have 1/8 of the force
  auto const tol = 1E-10;
  int count = 0;
  for (auto const &c : corner_nodes(params.grid_dimensions)) {
    auto const res = lb->get_node_force_to_be_applied(c);
    if (res) {
      BOOST_CHECK_SMALL(((*res) - f / 8.0).norm(), tol);
      count += 1;
    }
  };
  boost::mpi::all_reduce(world, boost::mpi::inplace(count), std::plus<int>());
  BOOST_CHECK_EQUAL(count, 8);

  lb->integrate();

  // check applied forces from last integration step
  count = 0;
  for (auto const &c : corner_nodes(params.grid_dimensions)) {
    auto const res = lb->get_node_last_applied_force(c);
    if (res) {
      BOOST_CHECK_SMALL(((*res) - f / 8.0).norm(), tol);
      count += 1;
    }
  };
  boost::mpi::all_reduce(world, boost::mpi::inplace(count), std::plus<int>());
  BOOST_CHECK_EQUAL(count, 8);
}

BOOST_DATA_TEST_CASE(vtk_exceptions,
                     bdata::make(LbGeneratorVector{unthermalized_lbs()[0]}),
                     lb_generator) {
  auto lb = lb_generator(mpi_shape, params);
  auto const flag =
      static_cast<std::underlying_type_t<OutputVTK>>(OutputVTK::density);
  // cannot create the same observable twice
  lb->create_vtk(1u, 0u, flag, "density", "vtk_out", "step");
  BOOST_CHECK_THROW(lb->create_vtk(1u, 0u, flag, "density", "vtk_out", "step"),
                    std::runtime_error);
  // cannot manually call an automatic observable
  lb->create_vtk(1u, 0u, flag, "auto", "vtk_out", "step");
  BOOST_CHECK_THROW(lb->write_vtk("vtk_out/auto"), std::runtime_error);
  // cannot activate a manual observable
  lb->create_vtk(0u, 0u, flag, "manual", "vtk_out", "step");
  BOOST_CHECK_THROW(lb->switch_vtk("vtk_out/manual", 0), std::runtime_error);
  // cannot call or activate observables that haven't been registered yet
  BOOST_CHECK_THROW(lb->write_vtk("unknown"), std::runtime_error);
  BOOST_CHECK_THROW(lb->switch_vtk("unknown", 0), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(macroscopic_accessor_equilibrium_distribution) {
  using namespace walberla::stencil;
  namespace EqDistr = walberla::lbm::accessor::EquilibriumDistribution;
  auto const x = std::sqrt(1. / 3.);
  auto const u = Utils::Vector3d::broadcast(x);
  auto const u_f = walberla::to_vector3<float>(u);
  auto const u_d = walberla::to_vector3<double>(u);
  auto const rho_f = 0.2f;
  auto const rho_d = 0.2;
  auto const tol_f = 100.f * 5e-7f;
  auto const tol_d = 100. * 5e-9;

  {
    auto const direction = Direction::C;
    auto const ref_d = rho_d * (1. / 3. - 1. / 2. * u.norm2());
    auto const ref_f = static_cast<float>(ref_d);
    auto const pop_f = EqDistr::get(direction, u_f, rho_f);
    auto const pop_d = EqDistr::get(direction, u_d, rho_d);
    BOOST_CHECK_CLOSE(pop_f, ref_f, tol_f);
    BOOST_CHECK_CLOSE(pop_d, ref_d, tol_d);
  }
  {
    auto const ref_d = rho_d * (1. / 18. - 1. / 12. * u.norm2() +
                                3. / 12. * x * x + 1. / 6. * x);
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::N, Direction::E, Direction::T}) {
      auto const pop_f = EqDistr::get(direction, u_f, rho_f);
      auto const pop_d = EqDistr::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(pop_f, ref_f, tol_f);
      BOOST_CHECK_CLOSE(pop_d, ref_d, tol_d);
    }
  }
  {
    auto const ref_d = -rho_d * (1. / 18. - 1. / 12. * u.norm2() -
                                 1. / 12. * x * x + 1. / 6. * x);
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::S, Direction::W, Direction::B}) {
      auto const pop_f = EqDistr::get(direction, u_f, rho_f);
      auto const pop_d = EqDistr::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(pop_f, ref_f, tol_f);
      BOOST_CHECK_CLOSE(pop_d, ref_d, tol_d);
    }
  }
  {
    auto const ref_d = -rho_d * (1. / 36. + 2. / 24. * u.norm2() -
                                 1. / 24. * x * x - 1. / 4. * x * x);
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::NW, Direction::SE, Direction::TS,
                                 Direction::TW, Direction::BN, Direction::BE}) {
      auto const pop_f = EqDistr::get(direction, u_f, rho_f);
      auto const pop_d = EqDistr::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(pop_f, ref_f, tol_f);
      BOOST_CHECK_CLOSE(pop_d, ref_d, tol_d);
    }
  }
  {
    auto const ref_d =
        rho_d * (1. / 36. + 2. / 24. * u.norm2() - 3. / 24. * x * x +
                 1. / 4. * x * x + 2. / 12. * x);
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::NE, Direction::TN, Direction::TE}) {
      auto const pop_f = EqDistr::get(direction, u_f, rho_f);
      auto const pop_d = EqDistr::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(pop_f, ref_f, tol_f);
      BOOST_CHECK_CLOSE(pop_d, ref_d, tol_d);
    }
  }
  {
    auto const ref_d =
        rho_d * (1. / 36. + 2. / 24. * u.norm2() - 3. / 24. * x * x +
                 1. / 4. * x * x - 2. / 12. * x);
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::SW, Direction::BS, Direction::BW}) {
      auto const pop_f = EqDistr::get(direction, u_f, rho_f);
      auto const pop_d = EqDistr::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(pop_f, ref_f, tol_f);
      BOOST_CHECK_CLOSE(pop_d, ref_d, tol_d);
    }
  }
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
  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
