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
#define BOOST_TEST_MODULE EK walberla node setters and getters test
#define BOOST_TEST_DYN_LINK
#include "config/config.hpp"

#ifdef WALBERLA

#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "tests_common_ek.hpp"

#include <walberla_bridge/VTKHandle.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/multi_array.hpp>

#include <mpi.h>

#include <cmath>
#include <functional>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

static EKTestParameters params; // populated in main()

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
  auto constexpr stencil_size = std::size_t{9u};
  auto ek = std::make_shared<walberla::EKinWalberlaImpl<stencil_size, float>>(
      params.lattice, params.diffusion, 0., params.valency, params.ext_efield,
      params.density, params.advection, params.friction_coupling);
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
      BOOST_CHECK(!ek->get_node_flux_at_boundary(node));
      BOOST_CHECK(!ek->remove_node_from_flux_boundary(node));
      BOOST_CHECK(!ek->get_node_is_flux_boundary(node));
    }
  }

  ek->clear_flux_boundaries();
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
  ek->create_vtk(1u, 0u, flag, units, "density", "vtk_out", "step");
  BOOST_CHECK_THROW(
      ek->create_vtk(1u, 0u, flag, units, "density", "vtk_out", "step"),
      std::runtime_error);
  // cannot manually call an automatic observable
  ek->create_vtk(1u, 0u, flag, units, "auto", "vtk_out", "step");
  BOOST_CHECK_THROW(ek->write_vtk("vtk_out/auto"), std::runtime_error);
  // cannot activate a manual observable
  ek->create_vtk(0u, 0u, flag, units, "manual", "vtk_out", "step");
  BOOST_CHECK_THROW(ek->switch_vtk("vtk_out/manual", 0), std::runtime_error);
  // cannot call or activate observables that haven't been registered yet
  BOOST_CHECK_THROW(ek->write_vtk("unknown"), std::runtime_error);
  BOOST_CHECK_THROW(ek->switch_vtk("unknown", 0), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ek_exceptions) {
  auto ek = std::make_shared<walberla::EKinWalberlaImpl<>>(
      params.lattice, params.diffusion, 0., params.valency, params.ext_efield,
      params.density, params.advection, params.friction_coupling);
  BOOST_CHECK_THROW(ek->integrate(std::size_t{}, std::size_t{}, std::size_t{}),
                    std::runtime_error);
  // no diffusion leads to early exit
  ek->set_diffusion(0.);
  ek->integrate(std::size_t{}, std::size_t{}, std::size_t{});
}

int main(int argc, char **argv) {
  int n_nodes;
  Vector3i mpi_shape{};

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
  params.lattice =
      std::make_shared<LatticeWalberla>(params.grid_dimensions, mpi_shape, 1u);

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // WALBERLA
int main(int argc, char **argv) {}
#endif
