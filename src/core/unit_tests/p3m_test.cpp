/*
 * Copyright (C) 2020-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE "P3M solvers"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ParticleFactory.hpp"
#include "particle_management.hpp"

#include "actor/registration.hpp"
#include "electrostatics/p3m.hpp"
#include "electrostatics/p3m.impl.definitions.hpp"
#include "magnetostatics/dp3m.hpp"
#include "magnetostatics/dp3m_heffte.impl.hpp"
#include "p3m/FFTBackendLegacy.hpp"
#include "p3m/FFTBuffersLegacy.hpp"
#include "p3m/common.hpp"

#include "EspressoCoreGlobalConfig.hpp"
#include "Particle.hpp"
#include "communication.hpp"
#include "energy_inline.hpp"
#include "integrate.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/demangle.hpp>
#include <utils/index.hpp>

#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

BOOST_AUTO_TEST_CASE(calc_meshift_false) {
  std::array<std::vector<int>, 3> const ref = {
      {std::vector<int>{0}, std::vector<int>{0, 1, -2, -1},
       std::vector<int>{0, 1, 2, 3, -3, -2, -1}}};

  auto const mesh = Utils::Vector3i{{1, 4, 7}};
  auto const val = calc_p3m_mesh_shift(mesh, false);

  for (std::size_t i = 0u; i < 3u; ++i) {
    for (std::size_t j = 0u; j < ref[i].size(); ++j) {
      BOOST_CHECK_EQUAL(val[i][j], ref[i][j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(calc_meshift_true) {
  std::array<std::vector<int>, 3> const ref = {
      {std::vector<int>{0}, std::vector<int>{0, 1, 0, -1},
       std::vector<int>{0, 1, 2, 0, -3, -2, -1}}};

  auto const mesh = Utils::Vector3i{{1, 4, 7}};
  auto const val = calc_p3m_mesh_shift(mesh, true);

  for (std::size_t i = 0u; i < 3u; ++i) {
    for (std::size_t j = 0u; j < ref[i].size(); ++j) {
      BOOST_CHECK_EQUAL(val[i][j], ref[i][j]);
    }
  }
}

#if defined(ESPRESSO_P3M) or defined(ESPRESSO_DP3M)

namespace espresso {
// ESPResSo system instance
static std::shared_ptr<System::System> system;
} // namespace espresso

struct GlobalConfig : public EspressoCoreGlobalConfig {
  GlobalConfig() {
    espresso::system = System::System::create();
    espresso::system->set_cell_structure_topology(CellStructureType::REGULAR);
    ::System::set_system(espresso::system);
  }
  ~GlobalConfig() {
    espresso::system.reset();
    ::System::reset_system();
  }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);
BOOST_AUTO_TEST_SUITE(suite)

#ifdef ESPRESSO_P3M
template <auto... Pack> void test_all_p3m_fft_configs() {
  auto constexpr nest_level = sizeof...(Pack);
  if constexpr (nest_level == 0) {
    test_all_p3m_fft_configs<Pack..., Utils::MemoryOrder::ROW_MAJOR>();
    test_all_p3m_fft_configs<Pack..., Utils::MemoryOrder::COLUMN_MAJOR>();
  }
  if constexpr (nest_level == 1) {
    test_all_p3m_fft_configs<Pack..., Utils::MemoryOrder::ROW_MAJOR>();
    test_all_p3m_fft_configs<Pack..., Utils::MemoryOrder::COLUMN_MAJOR>();
  }
  if constexpr (nest_level == 2) {
    test_all_p3m_fft_configs<Pack..., true>();
    test_all_p3m_fft_configs<Pack..., false>();
  }
  if constexpr (nest_level == 3) {
    // in the complex-to-complex backend, short dimension is irrelevant,
    // assuming the flat index is properly incremented in the energy and
    // pressure kernels (i.e. outside the short dimension conditional!!)
    test_all_p3m_fft_configs<Pack..., 0>();
    test_all_p3m_fft_configs<Pack..., 1>();
    test_all_p3m_fft_configs<Pack..., 2>();
  }
  if constexpr (nest_level == 4) {
    using FFTConfig = P3MFFTConfig<Pack...>;
    auto constexpr Hardware = Arch::CPU;
    auto const comm = boost::mpi::communicator();
    auto const rank = comm.rank();
    auto &system = *espresso::system;
    // set up P3M
    auto const prefactor = 2.;
    auto tuning = TuningParameters{1, {std::nullopt, std::nullopt}, false};
    auto solver = std::make_shared<CoulombP3MImpl<double, Hardware, FFTConfig>>(
        std::make_unique<CoulombP3MState<double, Hardware, FFTConfig>>(
            P3MParameters{false, 0.0, 3.5, Utils::Vector3i::broadcast(12),
                          Utils::Vector3d::broadcast(0.5), 5, 0.615, 1e-3}),
        tuning, prefactor);
    add_actor(comm, espresso::system, system.coulomb.impl->solver, solver,
              [&system]() { system.on_coulomb_change(); });
    auto const &obs_energy = system.calculate_energy();
    system.integrate(0, INTEG_REUSE_FORCES_NEVER);
    if (rank == 0) {
      BOOST_TEST_CONTEXT(Utils::demangle<FFTConfig>()) {
        auto constexpr energy_ref = -0.114744451;
        auto const energy_k_space = obs_energy.coulomb[1];
        BOOST_CHECK_CLOSE(energy_k_space, energy_ref, 1e-6);
      }
    }
    // check FFT properties
    auto const &fft = *solver->p3m.fft;
    auto const ref_rs_local_size = 12 / ::communicator.node_grid;
    BOOST_CHECK_EQUAL(fft.rs_local_size(), ref_rs_local_size);
    // the k-space boxes follow the last pencil stage of the heFFTe plan on
    // several ranks: check that they partition the k-space mesh and that
    // every non-empty box is a pencil (spans the full mesh in one direction)
    auto ref_ks_size = Utils::Vector3i::broadcast(12);
    if constexpr (FFTConfig::use_r2c) {
      ref_ks_size[FFTConfig::r2c_dir] = 12 / 2 + 1;
    }
    auto const ks_ld = fft.ks_local_ld_index();
    auto const ks_ur = fft.ks_local_ur_index();
    std::array<int, 6> const local_box = {ks_ld[0], ks_ld[1], ks_ld[2],
                                          ks_ur[0], ks_ur[1], ks_ur[2]};
    std::vector<int> flat_boxes;
    boost::mpi::all_gather(comm, local_box.data(), 6, flat_boxes);
    BOOST_REQUIRE_EQUAL(flat_boxes.size(), 6u * comm.size());
    auto const get_box = [&flat_boxes](int i) {
      auto const it = flat_boxes.begin() + 6 * i;
      return std::make_pair(Utils::Vector3i(it, it + 3),
                            Utils::Vector3i(it + 3, it + 6));
    };
    auto const volume = [](Utils::Vector3i const &ld,
                           Utils::Vector3i const &ur) {
      auto v = 1;
      for (auto d : {0, 1, 2}) {
        v *= std::max(0, ur[d] - ld[d]);
      }
      return v;
    };
    auto total_volume = 0;
    for (int i = 0; i < comm.size(); ++i) {
      BOOST_TEST_CONTEXT(Utils::demangle<FFTConfig>() << " box " << i) {
        auto const [ld, ur] = get_box(i);
        for (auto d : {0, 1, 2}) {
          BOOST_CHECK_GE(ld[d], 0);
          BOOST_CHECK_LE(ld[d], ur[d]);
          BOOST_CHECK_LE(ur[d], ref_ks_size[d]);
        }
        auto const v = volume(ld, ur);
        total_volume += v;
        if (v != 0) {
          auto is_pencil = false;
          for (auto d : {0, 1, 2}) {
            is_pencil |= (ld[d] == 0 and ur[d] == ref_ks_size[d]);
          }
          BOOST_CHECK(is_pencil);
        }
        for (int j = 0; j < i; ++j) {
          auto const [ld2, ur2] = get_box(j);
          auto overlap = 1;
          for (auto d : {0, 1, 2}) {
            overlap *=
                std::max(0, std::min(ur[d], ur2[d]) - std::max(ld[d], ld2[d]));
          }
          BOOST_CHECK_EQUAL(overlap, 0);
        }
      }
    }
    BOOST_CHECK_EQUAL(total_volume, Utils::product(ref_ks_size));
    if (::communicator.size == 1) {
      BOOST_CHECK_EQUAL(fft.ks_local_size(), ref_ks_size);
    }
    // deactivate actor
    solver->detach_system(espresso::system);
    system.coulomb.impl->solver = std::nullopt;
    system.on_coulomb_change();
  }
}
#endif // ESPRESSO_P3M

#ifdef ESPRESSO_DP3M
template <auto... Pack> void test_all_dp3m_fft_configs() {
  auto constexpr nest_level = sizeof...(Pack);
  if constexpr (nest_level == 0) {
    test_all_dp3m_fft_configs<Pack..., Utils::MemoryOrder::ROW_MAJOR>();
  }
  if constexpr (nest_level == 1) {
    test_all_dp3m_fft_configs<Pack..., Utils::MemoryOrder::ROW_MAJOR>();
  }
  if constexpr (nest_level == 2) {
    test_all_dp3m_fft_configs<Pack..., true>();
    test_all_dp3m_fft_configs<Pack..., false>();
  }
  if constexpr (nest_level == 3) {
    // in the complex-to-complex backend, short dimension is irrelevant,
    // assuming the flat index is properly incremented in the energy
    // kernel (i.e. outside the short dimension conditional!!)
    test_all_dp3m_fft_configs<Pack..., 0>();
    test_all_dp3m_fft_configs<Pack..., 1>();
    test_all_dp3m_fft_configs<Pack..., 2>();
  }
  if constexpr (nest_level == 4) {
    using FFTConfig = P3MFFTConfig<Pack...>;
    auto constexpr Hardware = Arch::CPU;
    auto const comm = boost::mpi::communicator();
    auto const rank = comm.rank();
    auto &system = *espresso::system;
    // set up P3M
    auto const prefactor = 2.;
    auto tuning = TuningParameters{1, {std::nullopt, std::nullopt}, false};
    auto solver =
        std::make_shared<DipolarP3MHeffte<double, Hardware, FFTConfig>>(
            std::make_unique<DipolarP3MState<double, FFTConfig>>(
                P3MParameters{false, 0.0, 3.5, Utils::Vector3i::broadcast(12),
                              Utils::Vector3d::broadcast(0.5), 5, 0.615, 1e-3}),
            tuning, prefactor);
    solver->dp3m.template make_mesh_instance<FFTBuffersLegacy<double>>();
    solver->dp3m.template make_fft_instance<FFTBackendLegacy<double>>();
    add_actor(comm, espresso::system, system.dipoles.impl->solver, solver,
              [&system]() { system.on_dipoles_change(); });
    auto const &obs_energy = system.calculate_energy();
    system.integrate(0, INTEG_REUSE_FORCES_NEVER);
    if (rank == 0) {
      BOOST_TEST_CONTEXT(Utils::demangle<FFTConfig>()) {
        auto constexpr energy_ref = -0.0052257342364;
        auto const energy_k_space = obs_energy.dipolar[1];
        BOOST_CHECK_CLOSE(energy_k_space, energy_ref, 1e-6);
      }
    }
    // deactivate actor
    solver->detach_system(espresso::system);
    system.dipoles.impl->solver = std::nullopt;
    system.on_dipoles_change();
  }
}
#endif // ESPRESSO_DP3M

BOOST_FIXTURE_TEST_CASE(p3m_solvers, ParticleFactory) {
  auto const box_l = 12.;
  auto const time_step = 0.001;
  auto const skin = 0.4;
  auto &system = *espresso::system;
  system.set_box_l(Utils::Vector3d::broadcast(box_l));
  system.set_time_step(time_step);
  system.cell_structure->set_verlet_skin(skin);

  // particle properties
  auto const pid1 = 9;
  auto const pid2 = 2;
  system.nonbonded_ias->make_particle_type_exist(0);

  // place particles farther apart than P3M real-space cutoff
  create_particle({box_l / 2. + 0.1, 0.1, 0.1}, pid1, 0);
  create_particle({0.1, 0.1, 0.1}, pid2, 0);

#ifdef ESPRESSO_P3M
  set_particle_property(pid1, &Particle::q, +0.5);
  set_particle_property(pid2, &Particle::q, -0.5);
  test_all_p3m_fft_configs<>();
#endif // ESPRESSO_P3M

#ifdef ESPRESSO_DP3M
  set_particle_property(pid1, &Particle::dipm, +0.5);
  set_particle_property(pid2, &Particle::dipm, -0.5);
  test_all_dp3m_fft_configs<>();
#endif // ESPRESSO_DP3M
}

BOOST_AUTO_TEST_SUITE_END()

#endif // defined(ESPRESSO_P3M) or defined(ESPRESSO_DP3M)
