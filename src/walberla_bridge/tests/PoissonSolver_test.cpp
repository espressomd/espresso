/*
 * Copyright (C) 2025 The ESPResSo project
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
#define BOOST_TEST_MODULE "Poisson Solver FFT test"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "tests_common_ek.hpp"

#include "../src/electrokinetics/PoissonSolverNone.hpp"
#include <walberla_bridge/VTKHandle.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <waLBerlaDefinitions.h>

#include <boost/mpl/list.hpp>

#if __has_include(<heffte.h>)
#define HAS_HEFFTE
#endif

#if defined(WALBERLA_BUILD_WITH_CUDA)
#include <cuda_runtime_api.h>
#endif

#include <mpi.h>

#include <cmath>
#include <functional>
#include <initializer_list>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;

static EKTestParameters params; // populated in main()

boost::test_tools::assertion_result has_gpu(boost::unit_test::test_unit_id) {
  bool has_compatible_device = false;
#if defined(WALBERLA_BUILD_WITH_CUDA)
  int n_devices = 0;
  cudaGetDeviceCount(&n_devices);
  if (n_devices > 0) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    if (prop.major >= 3) {
      has_compatible_device = true;
    }
  }
#endif
  return has_compatible_device;
}

template <lbmpy::Arch Architecture = lbmpy::Arch::CPU> struct Factory {
  static auto constexpr new_ek_species = walberla::new_ek_walberla_cpu;
  static auto constexpr new_ek_solver = walberla::new_ek_poisson_fft;
};

#if defined(WALBERLA_BUILD_WITH_CUDA)
template <> struct Factory<lbmpy::Arch::GPU> {
  static auto constexpr new_ek_species = walberla::new_ek_walberla_gpu;
  static auto constexpr new_ek_solver = walberla::new_ek_poisson_fft_cuda;
};
#endif

template <typename FloatType, lbmpy::Arch Architecture> struct Fixture {

  void runTest() {
    auto constexpr single_precision = std::is_same_v<FloatType, float>;
    auto constexpr eps = single_precision ? 6e-4 : 1e-8;
    auto ek_species_handle = Factory<Architecture>::new_ek_species(
        params.lattice, params.diffusion, 0., params.valency, params.ext_efield,
        params.density, params.advection, params.friction_coupling,
        single_precision, false, 0u);
    auto ek_solver_handle = Factory<Architecture>::new_ek_solver(
        params.lattice, 0.01, single_precision);
    ek_solver_handle->setup_fft(false);
    auto &ek_species = *ek_species_handle;
    auto &ek_solver = *ek_solver_handle;
    BOOST_CHECK_EQUAL(ek_solver.is_double_precision(), not single_precision);
    BOOST_CHECK_EQUAL(ek_solver.is_gpu(), Architecture == lbmpy::Arch::GPU);
    std::vector<double> density{};
    std::vector<double> potential_ref{};
    auto const &[lc, uc] = params.lattice->get_local_grid_range();
    auto const kw = 2. * std::numbers::pi / params.box_dimensions; // wavevector
    auto const phi0 = single_precision ? 15.2319936752 : 15.2319875119;
    for (int i = lc[0]; i < uc[0]; ++i) {
      for (int j = lc[1]; j < uc[1]; ++j) {
        for (int k = lc[2]; k < uc[2]; ++k) {
          auto const signal =
              std::cos(kw[0] * i) * std::cos(kw[1] * j) * std::cos(kw[2] * k);
          density.emplace_back(2. + signal);
          potential_ref.emplace_back(phi0 * signal);
        }
      }
    }
    ek_species.set_slice_density(lc, uc, density);
    ek_solver.add_charge_to_field(ek_species.get_density_id(), 0.1, false);
    ek_solver.solve();
    auto const potential_calc_local = ek_solver.get_slice_potential(lc, uc);
    for (std::size_t i = 0; i < potential_calc_local.size(); ++i) {
      BOOST_CHECK_SMALL(potential_ref[i] - potential_calc_local[i], eps);
    }
    // check node access
    auto const dim = uc - lc;
    for (auto const &node : all_nodes_incl_ghosts(*params.lattice)) {
      if (params.lattice->node_in_local_halo(node)) {
        if (params.lattice->node_in_local_domain(node)) {
          auto const found = ek_solver.get_node_potential(node, false);
          BOOST_REQUIRE(found);
          auto const node_pot = *found;
          auto const coord = node - lc;
          auto const i = (coord[0] * dim[1] + coord[1]) * dim[2] + coord[2];
          BOOST_CHECK_SMALL(potential_ref[i] - node_pot, eps);
        } else {
          auto const found = ek_solver.get_node_potential(node, false);
          BOOST_REQUIRE(not found);
          auto const found_in_halo = ek_solver.get_node_potential(node, true);
          BOOST_REQUIRE(found_in_halo);
        }
      } else {
        auto const found = ek_solver.get_node_potential(node, true);
        BOOST_REQUIRE(not found);
      }
    }
  }
};

BOOST_AUTO_TEST_SUITE(suite)

BOOST_AUTO_TEST_CASE(ek_poisson_solver_none) {
  auto ek_solver = walberla::PoissonSolverNone<double>(params.lattice);
  ek_solver.setup_fft(false);
  BOOST_CHECK(ek_solver.is_double_precision());
  BOOST_CHECK(not ek_solver.is_gpu());
  // no-op
  ek_solver.add_charge_to_field(std::size_t{}, 0., false);
  ek_solver.reset_charge_field();
  ek_solver.solve();
  // exceptions
  BOOST_CHECK_THROW(ek_solver.get_node_potential({0, 0, 0}, true),
                    std::runtime_error);
  BOOST_CHECK_THROW(ek_solver.get_slice_potential({0, 0, 0}, {1, 1, 1}),
                    std::runtime_error);
}

using test_types = boost::mpl::list<float, double>;

BOOST_AUTO_TEST_CASE_TEMPLATE(ek_poisson_solver_fft, FT, test_types) {
#if defined(HAS_HEFFTE)
  Fixture<FT, lbmpy::Arch::CPU>().runTest();
#endif
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(suite_cuda, *boost::unit_test::precondition(has_gpu))

using test_types = boost::mpl::list<float, double>;

BOOST_AUTO_TEST_CASE_TEMPLATE(ek_poisson_solver_fft_cuda, FT, test_types) {
#if defined(HAS_HEFFTE) and defined(WALBERLA_BUILD_WITH_CUDA)
  Fixture<FT, lbmpy::Arch::GPU>().runTest();
#endif
}

BOOST_AUTO_TEST_SUITE_END()

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
  params.lattice = std::make_shared<LatticeWalberla>(params.grid_dimensions,
                                                     mpi_shape, mpi_shape, 1u);

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  params.lattice.reset(); // free GPU fields before MPI_Finalize
  MPI_Finalize();
  return res;
}
