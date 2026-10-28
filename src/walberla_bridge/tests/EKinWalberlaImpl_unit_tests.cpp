/*
 * Copyright (C) 2019-2026 The ESPResSo project
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
#define BOOST_TEST_MODULE EK walberla node setters and getters test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "tests_common_ek.hpp"

#include "../src/electrokinetics/PoissonSolverNone.hpp"
#include <walberla_bridge/VTKHandle.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/multi_array.hpp>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

static EKTestParameters params; // populated in main()
static Vector3i mpi_shape{};    // populated in main()

BOOST_DATA_TEST_CASE(dimensions, bdata::make(all_eks()), ek_generator) {
  using boost::test_tools::per_element;
  auto ek = ek_generator(params);
  auto constexpr zero = Vector3i{0, 0, 0};

  auto const grid_dim = ek->get_lattice().get_grid_dimensions();
  BOOST_TEST(grid_dim == params.grid_dimensions, per_element());

  auto const [my_left, my_right] = ek->get_lattice().get_local_domain();
  auto const my_size = my_right - my_left;
  BOOST_TEST(my_size > zero, per_element());
  BOOST_TEST(my_left >= zero, per_element());
  BOOST_TEST(my_right <= params.grid_dimensions, per_element());
}

BOOST_AUTO_TEST_CASE(stencil_size) {
  auto constexpr stencil_size = std::size_t{13u};
  auto ek = std::make_shared<walberla::EKinWalberlaImpl<stencil_size, float>>(
      params.lattice, params.diffusion, 0., params.valency, params.ext_efield,
      params.density, params.advection, params.friction_coupling, false, 0u);
  BOOST_CHECK_EQUAL(ek->stencil_size(), stencil_size);
}

BOOST_DATA_TEST_CASE(set_diffusion, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  auto new_diffusion = 0.005;
  ek->set_diffusion(new_diffusion);
  BOOST_CHECK_CLOSE(ek->get_diffusion(), new_diffusion, 1E-11);
}

BOOST_DATA_TEST_CASE(set_valency, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  auto new_valency = 2.;
  ek->set_valency(new_valency);
  BOOST_CHECK_CLOSE(ek->get_valency(), new_valency, 1E-11);
}

BOOST_DATA_TEST_CASE(set_kT, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  auto new_kT = 2.;
  ek->set_kT(new_kT);
  BOOST_CHECK_CLOSE(ek->get_kT(), new_kT, 1E-11);
}

BOOST_DATA_TEST_CASE(set_advection, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  auto new_advection = false;
  ek->set_advection(new_advection);
  BOOST_CHECK_EQUAL(ek->get_advection(), new_advection);
}

BOOST_DATA_TEST_CASE(set_coupling, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  auto new_friction_coupling = false;
  ek->set_friction_coupling(new_friction_coupling);
  BOOST_CHECK_EQUAL(ek->get_friction_coupling(), new_friction_coupling);
}

BOOST_DATA_TEST_CASE(initial_state, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  for (auto const &node : local_nodes_incl_ghosts(ek->get_lattice())) {
    auto const consider_ghosts = !ek->get_lattice().node_in_local_domain(node);
    BOOST_CHECK(!(*ek->get_node_is_boundary(node, consider_ghosts)));
    if (ek->get_lattice().node_in_local_domain(node)) {
      BOOST_CHECK_CLOSE((*ek->get_node_density(node)), params.density, 1E-10);
    }
  }

  BOOST_CHECK_CLOSE(ek->get_diffusion(), params.diffusion, 1E-11);
  BOOST_CHECK_CLOSE(ek->get_valency(), params.valency, 1E-11);
  BOOST_CHECK_EQUAL(ek->get_advection(), params.advection);
  BOOST_CHECK_EQUAL(ek->get_friction_coupling(), params.friction_coupling);
}

BOOST_DATA_TEST_CASE(kT_unthermalized, bdata::make(unthermalized_eks()),
                     ek_generator) {
  auto ek = ek_generator(params);
  BOOST_CHECK_EQUAL(ek->get_kT(), 0.);
}

BOOST_DATA_TEST_CASE(kT_thermalized, bdata::make(thermalized_eks()),
                     ek_generator) {
  auto ek = ek_generator(params);
  BOOST_CHECK_EQUAL(ek->get_kT(), params.kT);
}

BOOST_DATA_TEST_CASE(node_flux_boundary, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  auto const flux = Vector3d{{0.2, 3.8, 4.2}};
  auto const n_ghost_layers =
      static_cast<int>(ek->get_lattice().get_ghost_layers());
  for (auto const &node : std::vector<Vector3i>{
           {-n_ghost_layers, 0, 0}, {0, 0, 0}, {0, 1, 2}, {9, 9, 9}}) {
    if (ek->get_lattice().node_in_local_halo(node)) {
      {
        auto const res = ek->get_node_is_boundary(node, true);
        // Did we get a value?
        BOOST_REQUIRE(res);
        // Should not be a boundary node
        BOOST_CHECK(*res == false);
      }
      {
        BOOST_CHECK(ek->set_node_flux_boundary(node, flux));
        ek->ghost_communication();
        {
          auto const res = ek->get_node_is_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == true);
        }
        {
          auto const res = ek->get_node_is_flux_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == true);
        }
        {
          auto const res = ek->get_node_is_density_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == false);
        }
      }
      {
        auto const flux_check = ek->get_node_flux_at_boundary(node, true);
        // Do we have a value
        BOOST_REQUIRE(flux_check);
        // Check the value
        BOOST_CHECK_SMALL((*flux_check - flux).norm(), 1E-12);
      }
      {
        BOOST_CHECK(ek->remove_node_from_flux_boundary(node));
        ek->ghost_communication();
        {
          auto const res = ek->get_node_is_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == false);
        }
        {
          auto const res = ek->get_node_is_flux_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == false);
        }
        {
          auto const res = ek->get_node_is_density_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == false);
        }
      }
    } else {
      // Not in the local halo.
      BOOST_CHECK(!ek->set_node_flux_boundary(node, flux));
      ek->ghost_communication();
      BOOST_CHECK(!ek->get_node_flux_at_boundary(node));
      BOOST_CHECK(!ek->get_node_flux_vector(node));
      BOOST_CHECK(!ek->remove_node_from_flux_boundary(node));
      ek->ghost_communication();
      BOOST_CHECK(!ek->get_node_is_flux_boundary(node));
    }
  }

  ek->clear_flux_boundaries();
  ek->ghost_communication();
  for (auto const &node : local_nodes_incl_ghosts(ek->get_lattice())) {
    BOOST_CHECK(!(*ek->get_node_is_flux_boundary(node, true)));
  }
}

BOOST_DATA_TEST_CASE(node_dens_boundary, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  auto const density = 0.2;
  auto const n_ghost_layers =
      static_cast<int>(ek->get_lattice().get_ghost_layers());
  for (auto const &node : std::vector<Vector3i>{
           {-n_ghost_layers, 0, 0}, {0, 0, 0}, {0, 1, 2}, {9, 9, 9}}) {
    if (ek->get_lattice().node_in_local_halo(node)) {
      {
        auto const res = ek->get_node_is_boundary(node, true);
        // Did we get a value?
        BOOST_REQUIRE(res);
        // Should not be a boundary node
        BOOST_CHECK(*res == false);
      }
      {
        BOOST_CHECK(ek->set_node_density_boundary(node, density));
        {
          auto const res = ek->get_node_is_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == true);
        }
        {
          auto const res = ek->get_node_is_density_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == true);
        }
        {
          auto const res = ek->get_node_is_flux_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == false);
        }
      }
      {
        auto const density_check = ek->get_node_density_at_boundary(node, true);
        // Do we have a value
        BOOST_REQUIRE(density_check);
        // Check the value
        BOOST_CHECK_SMALL(std::abs(*density_check - density), 1E-12);
      }
      {
        BOOST_CHECK(ek->remove_node_from_density_boundary(node));
        {
          auto const res = ek->get_node_is_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == false);
        }
        {
          auto const res = ek->get_node_is_density_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == false);
        }
        {
          auto const res = ek->get_node_is_flux_boundary(node, true);
          BOOST_REQUIRE(res);
          BOOST_CHECK(*res == false);
        }
      }
    } else {
      // Not in the local halo.
      BOOST_CHECK(!ek->set_node_density_boundary(node, density));
      BOOST_CHECK(!ek->get_node_density_at_boundary(node));
      BOOST_CHECK(!ek->remove_node_from_density_boundary(node));
      BOOST_CHECK(!ek->get_node_is_density_boundary(node));
    }
  }

  ek->clear_density_boundaries();
  for (auto const &node : local_nodes_incl_ghosts(ek->get_lattice())) {
    BOOST_CHECK(!(*ek->get_node_is_density_boundary(node, true)));
  }
}

BOOST_DATA_TEST_CASE(update_flux_boundary_from_shape, bdata::make(all_eks()),
                     ek_generator) {
  auto ek = ek_generator(params);
  auto const n_ghost_layers =
      static_cast<int>(ek->get_lattice().get_ghost_layers());
  auto const flux = Vector3d{{0.2, 3.8, 4.2}};

  auto const vec3to4 = [](Utils::Vector<int, 3> const &d, int v) {
    return Utils::Vector<int, 4>{{d[0], d[1], d[2], v}};
  };

  auto const nodes = std::vector<Vector3i>{
      {-n_ghost_layers, 0, 0}, {0, 0, 0}, {0, 1, 2}, {9, 9, 9}};
  // set up boundary
  {
    auto const n_grid_points = Utils::product(params.grid_dimensions);
    boost::multi_array<int, 3> raster_3d(params.grid_dimensions);
    boost::multi_array<double, 4> flux_3d(vec3to4(params.grid_dimensions, 3));
    BOOST_CHECK_EQUAL(raster_3d.num_elements(), n_grid_points);
    for (auto const &node : nodes) {
      auto const idx = (node + params.grid_dimensions) % params.grid_dimensions;
      raster_3d(idx) = 1;
      for (auto const i : {0, 1, 2}) {
        flux_3d(vec3to4(idx, i)) = flux[i];
      }
    }
    std::vector<int> raster_flat(raster_3d.data(),
                                 raster_3d.data() + raster_3d.num_elements());
    std::vector<double> flux_flat(flux_3d.data(),
                                  flux_3d.data() + flux_3d.num_elements());
    ek->update_flux_boundary_from_shape(raster_flat, flux_flat);
    ek->ghost_communication();
  }

  for (auto const &node : nodes) {
    if (ek->get_lattice().node_in_local_halo(node)) {
      {
        auto const res = ek->get_node_is_boundary(node, true);
        BOOST_REQUIRE(res);
        BOOST_CHECK(*res == true);
      }
      {
        auto const res = ek->get_node_is_flux_boundary(node, true);
        BOOST_REQUIRE(res);
        BOOST_CHECK(*res == true);
      }
      {
        auto const res = ek->get_node_is_density_boundary(node, true);
        BOOST_REQUIRE(res);
        BOOST_CHECK(*res == false);
      }
      {
        auto const flux_check = ek->get_node_flux_at_boundary(node, true);
        // Do we have a value
        BOOST_REQUIRE(flux_check);
        // Check the value
        BOOST_CHECK_SMALL((*flux_check - flux).norm(), 1E-12);
      }
    } else {
      // Not in the local halo.
      BOOST_CHECK(!ek->get_node_flux_at_boundary(node));
    }
  }

  ek->clear_flux_boundaries();
  ek->ghost_communication();
  for (auto const &node : local_nodes_incl_ghosts(ek->get_lattice())) {
    BOOST_CHECK(!(*ek->get_node_is_flux_boundary(node, true)));
  }
}

BOOST_DATA_TEST_CASE(update_density_boundary_from_shape, bdata::make(all_eks()),
                     ek_generator) {
  auto ek = ek_generator(params);
  auto const n_ghost_layers =
      static_cast<int>(ek->get_lattice().get_ghost_layers());
  auto const density = 0.2;

  auto const nodes = std::vector<Vector3i>{
      {-n_ghost_layers, 0, 0}, {0, 0, 0}, {0, 1, 2}, {9, 9, 9}};
  // set up boundary
  {
    auto const n_grid_points = Utils::product(params.grid_dimensions);
    boost::multi_array<int, 3> raster_3d(params.grid_dimensions);
    boost::multi_array<double, 3> dens_3d(params.grid_dimensions);
    BOOST_CHECK_EQUAL(raster_3d.num_elements(), n_grid_points);
    for (auto const &node : nodes) {
      auto const idx = (node + params.grid_dimensions) % params.grid_dimensions;
      raster_3d(idx) = 1;
      dens_3d(idx) = density;
    }
    std::vector<int> raster_flat(raster_3d.data(),
                                 raster_3d.data() + raster_3d.num_elements());
    std::vector<double> dens_flat(dens_3d.data(),
                                  dens_3d.data() + dens_3d.num_elements());
    ek->update_density_boundary_from_shape(raster_flat, dens_flat);
  }

  for (auto const &node : nodes) {
    if (ek->get_lattice().node_in_local_halo(node)) {
      {
        auto const res = ek->get_node_is_boundary(node, true);
        BOOST_REQUIRE(res);
        BOOST_CHECK(*res == true);
      }
      {
        auto const res = ek->get_node_is_density_boundary(node, true);
        BOOST_REQUIRE(res);
        BOOST_CHECK(*res == true);
      }
      {
        auto const res = ek->get_node_is_flux_boundary(node, true);
        BOOST_REQUIRE(res);
        BOOST_CHECK(*res == false);
      }
      {
        auto const density_check = ek->get_node_density_at_boundary(node, true);
        // Do we have a value
        BOOST_REQUIRE(density_check);
        // Check the value
        BOOST_CHECK_SMALL(std::abs(*density_check - density), 1E-12);
      }
    } else {
      // Not in the local halo.
      BOOST_CHECK(!ek->get_node_density_at_boundary(node));
    }
  }

  ek->clear_density_boundaries();
  ek->ghost_communication();
  for (auto const &node : local_nodes_incl_ghosts(ek->get_lattice())) {
    BOOST_CHECK(!(*ek->get_node_is_density_boundary(node, true)));
  }
}

BOOST_DATA_TEST_CASE(domain_and_halo, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);
  auto const n_ghost_layers = ek->get_lattice().get_ghost_layers();
  auto const [my_left, my_right] = ek->get_lattice().get_local_domain();

  for (auto const &n : all_nodes_incl_ghosts(ek->get_lattice())) {
    auto const pos = n + Vector3d::broadcast(.5);
    int is_local = 0;
    // Nodes in local domain
    if (Vector3d(n) >= my_left and Vector3d(n) < my_right) {
      BOOST_CHECK(ek->get_lattice().node_in_local_domain(n));
      BOOST_CHECK(ek->get_lattice().node_in_local_halo(n));

      BOOST_CHECK(ek->get_lattice().pos_in_local_domain(pos));
      BOOST_CHECK(ek->get_lattice().pos_in_local_halo(pos));
      is_local = 1;
    } else {
      // in local halo?
      if ((n + Vector3d::broadcast(n_ghost_layers)) >= my_left and
          (n - Vector3d::broadcast(n_ghost_layers)) < my_right) {
        BOOST_CHECK(!ek->get_lattice().node_in_local_domain(n));
        BOOST_CHECK(ek->get_lattice().node_in_local_halo(n));

        BOOST_CHECK(!ek->get_lattice().pos_in_local_domain(pos));
        BOOST_CHECK(ek->get_lattice().pos_in_local_halo(pos));
      } else {
        // neither in domain nor in halo
        BOOST_CHECK(!ek->get_lattice().node_in_local_domain(n));
        BOOST_CHECK(!ek->get_lattice().node_in_local_halo(n));

        BOOST_CHECK(!ek->get_lattice().pos_in_local_domain(pos));
        BOOST_CHECK(!ek->get_lattice().pos_in_local_halo(pos));
      }
    }

    // If the cell is in the global physical domain
    // check that only one mpi rank said the node was local
    auto constexpr origin = Vector3i{0, 0, 0};
    if (n >= origin and n < params.grid_dimensions) {
      boost::mpi::communicator world;
      auto const is_local_sum =
          boost::mpi::all_reduce(world, is_local, std::plus<int>());
      BOOST_CHECK(is_local_sum == 1);
    }
  }
}

static auto fold_node(Vector3i n) {
  for (unsigned int i = 0; i < 3; i++) {
    if (n[i] < 0) {
      n[i] += params.grid_dimensions[i];
    } else if (n[i] >= params.grid_dimensions[i]) {
      n[i] -= params.grid_dimensions[i];
    }
  }
  return n;
}

BOOST_DATA_TEST_CASE(set_node_density, bdata::make(all_eks()), ek_generator) {
  auto ek = ek_generator(params);

  auto n_dens = [](Vector3i const &node) {
    return 1. + static_cast<double>(Utils::product(fold_node(node))) * 1e-6;
  };

  // Assign densities
  for (auto const &node : all_nodes_incl_ghosts(ek->get_lattice())) {
    if (ek->get_lattice().node_in_local_domain(node)) {
      BOOST_CHECK(ek->set_node_density(node, n_dens(node)));
    } else {
      // Check that access to node density is not possible
      BOOST_CHECK(!ek->set_node_density(node, 0.));
    }
  }

  ek->ghost_communication();

  // check densities
  for (auto const &node : all_nodes_incl_ghosts(ek->get_lattice())) {
    auto constexpr eps = 1E-8;
    if (ek->get_lattice().node_in_local_halo(node)) {
      auto const consider_ghosts =
          !ek->get_lattice().node_in_local_domain(node);
      auto res = ek->get_node_density(node, consider_ghosts);
      BOOST_REQUIRE(res);                          // value available?
      BOOST_CHECK_SMALL(*res - n_dens(node), eps); // value correct?
    } else {
      // Check that access to node density is not possible
      BOOST_CHECK(!ek->get_node_density(node));
    }
  }
}

BOOST_DATA_TEST_CASE(vtk_exceptions,
                     bdata::make(EkGeneratorVector{unthermalized_eks()[0]}),
                     ek_generator) {
  std::unordered_map<std::string, double> const units = {{"density", 1.}};
  auto ek = ek_generator(params);
  auto const flag =
      static_cast<std::underlying_type_t<OutputVTK>>(OutputVTK::density);
  // cannot create the same observable twice
  ek->create_vtk(1u, 0u, flag, units, "density", "vtk_out", "step", false);
  BOOST_CHECK_THROW(
      ek->create_vtk(1u, 0u, flag, units, "density", "vtk_out", "step", false),
      std::runtime_error);
  // cannot manually call an automatic observable
  ek->create_vtk(1u, 0u, flag, units, "auto", "vtk_out", "step", false);
  BOOST_CHECK_THROW(ek->write_vtk("vtk_out/auto"), std::runtime_error);
  // cannot activate a manual observable
  ek->create_vtk(0u, 0u, flag, units, "manual", "vtk_out", "step", false);
  BOOST_CHECK_THROW(ek->switch_vtk("vtk_out/manual", 0), std::runtime_error);
  // cannot call or activate observables that haven't been registered yet
  BOOST_CHECK_THROW(ek->write_vtk("unknown"), std::runtime_error);
  BOOST_CHECK_THROW(ek->switch_vtk("unknown", 0), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ek_exceptions) {
  auto ek = std::make_shared<walberla::EKinWalberlaImpl<>>(
      params.lattice, params.diffusion, 0., params.valency, params.ext_efield,
      params.density, params.advection, params.friction_coupling, false, 0u);
  BOOST_CHECK_THROW(
      ek->integrate(std::size_t{}, std::size_t{}, std::size_t{}, 0.),
      std::runtime_error);
  // no diffusion leads to early exit
  ek->set_diffusion(0.);
  ek->integrate(std::size_t{}, std::size_t{}, std::size_t{}, 0.);
}

BOOST_AUTO_TEST_CASE(ek_poisson_solver_none) {
  auto ek_solver =
      walberla::PoissonSolverNone<double, lbmpy::Arch::CPU>(params.lattice);
  // no-op
  ek_solver.add_charge_to_field(std::size_t{}, 0.);
  ek_solver.reset_charge_field();
  ek_solver.solve();
  BOOST_CHECK(not ek_solver.get_node_potential({9999, 9999, 9999}, true));
  BOOST_CHECK(
      ek_solver.get_slice_potential({9999, 9999, 9999}, {10000, 10000, 10000})
          .empty());
}

namespace {

/** EK of the species of the tests without advection and friction coupling. */
auto make_ek(double diffusion, double kT, double valency,
             Vector3d const &ext_efield, bool thermalized,
             std::shared_ptr<LatticeWalberla> lattice = params.lattice) {
  return std::make_shared<walberla::EKinWalberlaImpl<>>(
      std::move(lattice), diffusion, kT, valency, ext_efield, params.density,
      false, false, thermalized, 42u);
}

/** Density of the nodes of the lattice of all MPI ranks, in all processes. */
std::vector<double> gather_density(EKinWalberlaBase const &ek) {
  auto const &grid = params.grid_dimensions;
  std::vector<double> local(Utils::product(grid), 0.);
  std::size_t index = 0;
  for (int x = 0; x < grid[0]; ++x) {
    for (int y = 0; y < grid[1]; ++y) {
      for (int z = 0; z < grid[2]; ++z, ++index) {
        if (auto const density = ek.get_node_density({x, y, z})) {
          local[index] = *density;
        }
      }
    }
  }
  std::vector<double> total(local.size());
  boost::mpi::communicator world;
  boost::mpi::all_reduce(world, local.data(), static_cast<int>(local.size()),
                         total.data(), std::plus<double>());
  BOOST_REQUIRE(std::ranges::all_of(
      total, [](double const value) { return std::isfinite(value); }));
  return total;
}

/** Set the density of the local nodes. */
template <typename Function>
void set_density(EKinWalberlaBase &ek, Function &&density) {
  for (auto const &node : local_nodes_incl_ghosts(ek.get_lattice(), false)) {
    ek.set_node_density(node, density(node));
  }
  ek.ghost_communication();
}

/** Integrate with the communication of the density of each time step. */
void integrate(EKinWalberlaBase &ek, int n_steps,
               std::size_t potential_id = 0u) {
  for (int step = 0; step < n_steps; ++step) {
    ek.ghost_communication();
    ek.integrate(potential_id, 0u, 0u, 1.);
  }
}

/**
 * Density that is large enough for the fluctuations of a few time steps to
 * leave it positive (the amplitude of the noise is proportional to the square
 * root of the density).
 */
double smooth_density(Vector3i const &node) {
  auto const &grid = params.grid_dimensions;
  return 20. +
         3. * std::sin(2. * M_PI * node[0] / grid[0]) *
             std::cos(2. * M_PI * node[2] / grid[2]) +
         std::sin(2. * M_PI * node[1] / grid[1]);
}

double smooth_potential(Vector3i const &node) {
  auto const &grid = params.grid_dimensions;
  return 0.02 * std::cos(2. * M_PI * node[0] / grid[0]) +
         0.01 * std::sin(2. * M_PI * node[1] / grid[1]) *
             std::cos(2. * M_PI * node[2] / grid[2]);
}

/** Potential field of a solver without charges, set in the local domain. */
std::shared_ptr<walberla::PoissonSolverNone<double, lbmpy::Arch::CPU>>
make_potential(double scale) {
  auto solver =
      std::make_shared<walberla::PoissonSolverNone<double, lbmpy::Arch::CPU>>(
          params.lattice);
  for (auto const &node : local_nodes_incl_ghosts(*params.lattice, false)) {
    solver->set_node_potential(node, scale * smooth_potential(node));
  }
  return solver;
}

} // namespace

BOOST_AUTO_TEST_CASE(diffusion_conserves_mass) {
  auto ek = make_ek(0.1, 0., 0., Vector3d{}, false);
  auto const spike = Vector3i{5, 6, 9};
  set_density(*ek,
              [&](Vector3i const &node) { return node == spike ? 2.4 : 1.4; });
  auto const total = [&]() {
    auto const density = gather_density(*ek);
    return std::accumulate(density.begin(), density.end(), 0.);
  };
  auto const mass = total();
  BOOST_CHECK_CLOSE(
      mass, 1.4 * double(Utils::product(params.grid_dimensions)) + 1., 1e-10);
  integrate(*ek, 4);
  BOOST_CHECK_CLOSE(total(), mass, 1e-10);
  // the spike spreads to the 26 neighbors, which are all the same distance
  // apart in the faces they share: the neighbors in the opposite directions
  // have the same density
  auto const density = gather_density(*ek);
  auto const index = [&](Vector3i const &node) {
    return static_cast<std::size_t>(
        (node[0] * params.grid_dimensions[1] + node[1]) *
            params.grid_dimensions[2] +
        node[2]);
  };
  BOOST_CHECK_LT(density[index(spike)], 2.4);
  for (auto const &d : {Vector3i{1, 0, 0}, Vector3i{0, 1, 0}, Vector3i{0, 0, 1},
                        Vector3i{1, 1, 0}, Vector3i{1, 0, -1},
                        Vector3i{1, 1, 1}, Vector3i{-1, 1, 1}}) {
    BOOST_CHECK_GT(density[index(spike + d)], 1.4);
    BOOST_CHECK_CLOSE(density[index(spike + d)], density[index(spike - d)],
                      1e-10);
  }
}

BOOST_AUTO_TEST_CASE(set_diffusion_updates_the_kernel) {
  // the persistent kernels of the EK are updated by the setters of their
  // parameters: a species whose diffusion coefficient was changed behaves
  // like a species that was created with it
  auto const reference = make_ek(0.02, 0., 0., Vector3d{}, false);
  auto const changed = make_ek(0.3, 0., 0., Vector3d{}, false);
  changed->set_diffusion(0.02);
  for (auto const &ek : {reference, changed}) {
    set_density(*ek, smooth_density);
    integrate(*ek, 3);
  }
  auto const expected = gather_density(*reference);
  BOOST_CHECK(gather_density(*changed) == expected);
  // the first and last cells of the profile differ from the initial profile
  auto const initial = make_ek(0.02, 0., 0., Vector3d{}, false);
  set_density(*initial, smooth_density);
  BOOST_CHECK(gather_density(*initial) != expected);
  // no more diffusion: no change
  changed->set_diffusion(0.);
  integrate(*changed, 2);
  BOOST_CHECK(gather_density(*changed) == expected);
}

BOOST_AUTO_TEST_CASE(electrostatic_kernel_parameters) {
  auto const solver = make_potential(1.);
  auto const phi = solver->get_potential_field_id();
  auto const ext_efield = Vector3d{0.01, -0.02, 0.03};
  auto const reference = make_ek(0.02, 0.5, 1., ext_efield, false);
  auto const changed = make_ek(0.4, 2., 3., Vector3d{0.2, 0.1, 0.3}, false);
  // the electrostatic kernel is constructed with the first integration, with
  // the parameters that the species has at that time: change them afterwards
  set_density(*changed, smooth_density);
  integrate(*changed, 1, phi);
  changed->set_diffusion(0.02);
  changed->set_kT(0.5);
  changed->set_valency(1.);
  changed->set_ext_efield(ext_efield);
  for (auto const &ek : {reference, changed}) {
    set_density(*ek, smooth_density);
    integrate(*ek, 3, phi);
  }
  auto const expected = gather_density(*reference);
  BOOST_CHECK(gather_density(*changed) == expected);
  BOOST_CHECK_CLOSE(changed->get_kT(), 0.5, 1e-10);
  BOOST_CHECK_CLOSE(changed->get_valency(), 1., 1e-10);
  BOOST_CHECK_CLOSE(changed->get_ext_efield()[2], 0.03, 1e-10);
  // the species drifts in the electric field
  auto const only_diffusion = make_ek(0.02, 0.5, 0., Vector3d{}, false);
  set_density(*only_diffusion, smooth_density);
  integrate(*only_diffusion, 3);
  BOOST_CHECK(gather_density(*only_diffusion) != expected);
}

BOOST_AUTO_TEST_CASE(potential_field_can_change) {
  // the ID of the potential field is read-only in a generated kernel: the
  // kernel is created again when the potential field changes
  auto const solver_a = make_potential(1.);
  auto const solver_b = make_potential(-3.);
  auto const phi_a = solver_a->get_potential_field_id();
  auto const phi_b = solver_b->get_potential_field_id();
  BOOST_REQUIRE_NE(phi_a, phi_b);
  auto const run = [&](std::vector<std::size_t> const &sequence) {
    auto ek = make_ek(0.02, 0.5, 1., Vector3d{0.01, 0., 0.}, false);
    set_density(*ek, smooth_density);
    for (auto const phi : sequence) {
      integrate(*ek, 1, phi);
    }
    return gather_density(*ek);
  };
  auto const a = run({phi_a, phi_a});
  auto const b = run({phi_b, phi_b});
  BOOST_CHECK(a != b);
  // the first step is the same as the one of the species that only
  // knows the first potential
  auto const ab = run({phi_a, phi_b});
  auto const ba = run({phi_b, phi_a});
  BOOST_CHECK(ab != a and ab != b and ab != ba);
  BOOST_CHECK(run({phi_a, phi_b}) == ab);
  BOOST_CHECK(run({phi_a, phi_b, phi_a, phi_a}) ==
              run({phi_a, phi_b, phi_a, phi_a}));
}

BOOST_AUTO_TEST_CASE(thermalized_flux_is_independent_of_the_blocks) {
  // the RNG of a cell depends on its global coordinates, and the flux field
  // is computed in the ghost layers: the densities don't depend on the
  // number of blocks, here one and eight blocks for each MPI rank
  auto const run = [](std::shared_ptr<LatticeWalberla> lattice) {
    auto const ek = make_ek(0.02, 0.5, 0., Vector3d{}, true, lattice);
    set_density(*ek, smooth_density);
    integrate(*ek, 3);
    return gather_density(*ek);
  };
  auto const lattice_of_eight_blocks = std::make_shared<LatticeWalberla>(
      params.grid_dimensions, mpi_shape, 2 * mpi_shape, 1u);
  auto const reference = run(params.lattice);
  BOOST_CHECK(run(lattice_of_eight_blocks) == reference);
  // the noise is there: the density is not the one without fluctuations
  auto const ek = make_ek(0.02, 0.5, 0., Vector3d{}, false);
  set_density(*ek, smooth_density);
  integrate(*ek, 3);
  BOOST_CHECK(gather_density(*ek) != reference);
}

BOOST_AUTO_TEST_CASE(rng_state) {
  auto const make = [](double valency) {
    auto ek = make_ek(0.02, 0.5, valency, Vector3d{}, true);
    set_density(*ek, smooth_density);
    return ek;
  };
  auto const ek = make(0.);
  BOOST_REQUIRE(ek->get_rng_state());
  BOOST_CHECK_EQUAL(*ek->get_rng_state(), 0u);
  integrate(*ek, 3);
  BOOST_CHECK_EQUAL(*ek->get_rng_state(), 3u);
  ek->set_rng_state(10u);
  BOOST_CHECK_EQUAL(*ek->get_rng_state(), 10u);
  integrate(*ek, 1);
  BOOST_CHECK_EQUAL(*ek->get_rng_state(), 11u);
  // the noise depends on the counter, and the species at the same counter
  // is reproducible
  auto const fluctuations = [&](std::uint64_t counter) {
    auto ek = make(0.);
    ek->set_rng_state(counter);
    integrate(*ek, 1);
    return gather_density(*ek);
  };
  BOOST_CHECK(fluctuations(7u) == fluctuations(7u));
  BOOST_CHECK(fluctuations(7u) != fluctuations(8u));
  // the counter advances in the same way with electrostatics
  auto const solver = make_potential(1.);
  auto const ek_electrostatic = make(1.);
  integrate(*ek_electrostatic, 2, solver->get_potential_field_id());
  BOOST_CHECK_EQUAL(*ek_electrostatic->get_rng_state(), 2u);
  // and the noise is reproducible with electrostatics too
  auto const electrostatic = [&](std::uint64_t counter) {
    auto ek = make(1.);
    ek->set_rng_state(counter);
    integrate(*ek, 1, solver->get_potential_field_id());
    return gather_density(*ek);
  };
  BOOST_CHECK(electrostatic(4u) == electrostatic(4u));
  BOOST_CHECK(electrostatic(4u) != electrostatic(5u));
  // an EK without fluctuations has no RNG
  auto const unthermalized = make_ek(0.02, 0.5, 0., Vector3d{}, false);
  BOOST_CHECK(not unthermalized->get_rng_state());
  BOOST_CHECK_THROW(unthermalized->set_rng_state(1u), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ek_without_ghost_layers) {
  auto const lattice = std::make_shared<LatticeWalberla>(
      params.grid_dimensions, mpi_shape, mpi_shape, 0u);
  BOOST_CHECK_THROW(walberla::EKinWalberlaImpl<>(
                        lattice, params.diffusion, 0., 0., params.ext_efield,
                        params.density, false, false, false, 0u),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ek_block_too_large_for_32bit_indices) {
  // blocks of 600^3 cells: the flux field would have 13 * 602^3 > 2^31
  // elements (creating the lattice doesn't allocate any field, and the EK
  // constructor must throw before it does)
  auto const grid_dimensions =
      hadamard_product(Vector3i::broadcast(600), mpi_shape);
  auto const lattice = std::make_shared<LatticeWalberla>(
      grid_dimensions, mpi_shape, mpi_shape, 1u);
  auto const check_message = [](std::runtime_error const &error) {
    std::string const message(error.what());
    return message.find("flux field") != std::string::npos and
           message.find("2^31") != std::string::npos;
  };
  BOOST_CHECK_EXCEPTION(walberla::EKinWalberlaImpl<>(lattice, params.diffusion,
                                                     0., 0., params.ext_efield,
                                                     params.density, false,
                                                     false, false, 0u),
                        std::runtime_error, check_message);
}

int main(int argc, char **argv) { // NOLINT(bugprone-exception-escape)
  int n_nodes;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());
  walberla::mpi_init();

  params.seed = 0u;
  params.kT = 1.3E-4;
  params.density = 1.4;
  params.diffusion = 0.003;
  params.valency = 1.;
  params.advection = true;
  params.friction_coupling = true;
  params.ext_efield = Vector3d{0.01, 0.02, 0.03};
  params.grid_dimensions = Vector3i{12, 12, 18};
  params.box_dimensions = Vector3d{12, 12, 18};
  params.lattice = std::make_shared<LatticeWalberla>(params.grid_dimensions,
                                                     mpi_shape, mpi_shape, 1u);

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}
