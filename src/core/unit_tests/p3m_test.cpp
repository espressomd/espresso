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

#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
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
/** @brief Mesh and charge assignment order of a P3M test. */
struct P3MTestMesh {
  Utils::Vector3i mesh;
  int cao;
};

/** @brief Run the force calculation with P3M in the given precision and
 *  FFT configuration.
 *  @param setup Mesh and charge assignment order.
 *  @param check Check of the k-space energy and the particle forces after
 *  the force calculation, called on all ranks.
 */
template <typename FloatType, class FFTConfig>
void run_p3m_config(P3MTestMesh const &setup, auto const &check) {
  auto constexpr Hardware = Arch::CPU;
  auto const comm = boost::mpi::communicator();
  auto &system = *espresso::system;
  // set up P3M
  auto const prefactor = 2.;
  auto tuning = TuningParameters{1, {std::nullopt, std::nullopt}, false};
  auto solver =
      std::make_shared<CoulombP3MImpl<FloatType, Hardware, FFTConfig>>(
          std::make_unique<CoulombP3MState<FloatType, Hardware, FFTConfig>>(
              P3MParameters{false, 0.0, 3.5, setup.mesh,
                            Utils::Vector3d::broadcast(0.5), setup.cao, 0.615,
                            1e-3}),
          tuning, prefactor);
  add_actor(comm, espresso::system, system.coulomb.impl->solver, solver,
            [&system]() { system.on_coulomb_change(); });
  auto const &obs_energy = system.calculate_energy();
  system.integrate(0, INTEG_REUSE_FORCES_NEVER);
  BOOST_TEST_CONTEXT(Utils::demangle<FloatType>()
                     << ", " << Utils::demangle<FFTConfig>()) {
    check(obs_energy.coulomb[1]);
  }
  // check FFT properties
  if (setup.mesh == Utils::Vector3i::broadcast(12)) {
    auto const &fft = *solver->p3m.fft;
    auto const ref_rs_local_size = 12 / ::communicator.node_grid;
    BOOST_CHECK_EQUAL(fft.rs_local_size(), ref_rs_local_size);
    if (::communicator.size <= 3) {
      auto const node_index = ::communicator.calc_node_index();
      auto ref_ks_local_size = 12 / ::communicator.node_grid;
      if constexpr (FFTConfig::use_r2c) {
        ref_ks_local_size[FFTConfig::r2c_dir] /= 2;
        if (node_index[FFTConfig::r2c_dir] == 0) {
          ref_ks_local_size[FFTConfig::r2c_dir] += 1;
        }
      }
      BOOST_CHECK_EQUAL(fft.ks_local_size(), ref_ks_local_size);
      BOOST_CHECK_EQUAL(fft.ks_local_ur_index() - fft.ks_local_ld_index(),
                        ref_ks_local_size);
    }
  }
  // deactivate actor
  solver->detach_system(espresso::system);
  system.coulomb.impl->solver = std::nullopt;
  system.on_coulomb_change();
}

/** @brief Run @ref run_p3m_config in double precision for all FFT
 *  configurations.
 */
template <auto... Pack>
void test_all_p3m_fft_configs(P3MTestMesh const &setup, auto const &check) {
  auto constexpr nest_level = sizeof...(Pack);
  if constexpr (nest_level == 0) {
    test_all_p3m_fft_configs<Pack..., Utils::MemoryOrder::ROW_MAJOR>(setup,
                                                                     check);
    test_all_p3m_fft_configs<Pack..., Utils::MemoryOrder::COLUMN_MAJOR>(setup,
                                                                        check);
  }
  if constexpr (nest_level == 1) {
    test_all_p3m_fft_configs<Pack..., Utils::MemoryOrder::ROW_MAJOR>(setup,
                                                                     check);
    test_all_p3m_fft_configs<Pack..., Utils::MemoryOrder::COLUMN_MAJOR>(setup,
                                                                        check);
  }
  if constexpr (nest_level == 2) {
    test_all_p3m_fft_configs<Pack..., true>(setup, check);
    test_all_p3m_fft_configs<Pack..., false>(setup, check);
  }
  if constexpr (nest_level == 3) {
    // in the complex-to-complex backend, short dimension is irrelevant,
    // assuming the flat index is properly incremented in the energy and
    // pressure kernels (i.e. outside the short dimension conditional!!)
    test_all_p3m_fft_configs<Pack..., 0>(setup, check);
    test_all_p3m_fft_configs<Pack..., 1>(setup, check);
    test_all_p3m_fft_configs<Pack..., 2>(setup, check);
  }
  if constexpr (nest_level == 4) {
    run_p3m_config<double, P3MFFTConfig<Pack...>>(setup, check);
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
  test_all_p3m_fft_configs<>(
      {Utils::Vector3i::broadcast(12), 5}, [](double energy_k_space) {
        if (boost::mpi::communicator().rank() == 0) {
          auto constexpr energy_ref = -0.114744451;
          BOOST_CHECK_CLOSE(energy_k_space, energy_ref, 1e-6);
        }
      });
#endif // ESPRESSO_P3M

#ifdef ESPRESSO_DP3M
  set_particle_property(pid1, &Particle::dipm, +0.5);
  set_particle_property(pid2, &Particle::dipm, -0.5);
  test_all_dp3m_fft_configs<>();
#endif // ESPRESSO_DP3M
}

#ifdef ESPRESSO_P3M
/** @brief Create charged particles and set up the system for the P3M
 *  force tests.
 */
static void setup_force_test(ParticleFactory &factory,
                             std::span<Utils::Vector3d const> positions,
                             std::span<double const> charges) {
  auto const box_l = 12.;
  auto &system = *espresso::system;
  system.set_box_l(Utils::Vector3d::broadcast(box_l));
  system.set_time_step(0.001);
  system.cell_structure->set_verlet_skin(0.4);
  system.nonbonded_ias->make_particle_type_exist(0);
  for (int pid = 0; pid < static_cast<int>(positions.size()); ++pid) {
    factory.create_particle(positions[pid], pid, 0);
    factory.set_particle_property(pid, &Particle::q, charges[pid]);
  }
}

/** @brief Check the particle forces against references.
 *  @param forces_ref Reference forces, in particle id order.
 *  @param tol        Tolerance relative to the largest reference force
 *                    component, used as an absolute tolerance for every
 *                    component.
 */
static void check_forces(std::span<Utils::Vector3d const> forces_ref,
                         double tol) {
  auto const comm = boost::mpi::communicator();
  auto &system = *espresso::system;
  auto f_max = 0.;
  for (auto const &f_ref : forces_ref) {
    for (auto d : {0u, 1u, 2u}) {
      f_max = std::max(f_max, std::abs(f_ref[d]));
    }
  }
  for (int pid = 0; pid < static_cast<int>(forces_ref.size()); ++pid) {
    auto const p_opt = copy_particle_to_head_node(comm, system, pid);
    if (comm.rank() == 0) {
      auto const &force = p_opt->force();
      for (auto d : {0u, 1u, 2u}) {
        BOOST_TEST_CONTEXT("particle " << pid << ", direction " << d) {
          BOOST_CHECK_SMALL(force[d] - forces_ref[pid][d], tol * f_max);
        }
      }
    }
  }
}

BOOST_FIXTURE_TEST_CASE(p3m_forces, ParticleFactory) {
  // particles farther apart than the P3M real-space cutoff, the forces are
  // the k-space forces; the stencils of the particles close to the box
  // boundaries reach into the lower and upper halo layers in all directions
  auto const positions = std::array<Utils::Vector3d, 4>{
      {{1.3, 2.7, 0.4}, {7.9, 3.1, 5.6}, {2.2, 8.6, 7.3}, {10.9, 11.6, 11.3}}};
  auto const charges = std::array<double, 4>{{+1., -1., +0.5, -0.5}};
  // reference forces of the halo exchange path
  auto const forces_ref = std::array<Utils::Vector3d, 4>{{
      {-3.948462637768e-02, -3.664156665389e-02, -2.366230322212e-03},
      {-1.796541823041e-03, 7.079949821642e-03, -5.895725754431e-03},
      {-1.684119697017e-03, 3.474991277469e-03, -2.461433231618e-03},
      {4.296528789774e-02, 2.608662555477e-02, 1.072338930826e-02},
  }};
  setup_force_test(*this, positions, charges);
  test_all_p3m_fft_configs<>({Utils::Vector3i::broadcast(12), 5},
                             [&](double) { check_forces(forces_ref, 1e-10); });
}

BOOST_FIXTURE_TEST_CASE(p3m_forces_non_cubic_mesh, ParticleFactory) {
  // even charge assignment order, i.e. different lower and upper halo
  // widths, and a different mesh size in every direction; every box face
  // has a particle whose stencil reaches into the halo layer behind it
  auto const positions = std::array<Utils::Vector3d, 4>{
      {{0.3, 0.2, 6.1}, {5.9, 11.8, 0.3}, {6.0, 5.8, 6.4}, {11.7, 6.2, 11.8}}};
  auto const charges = std::array<double, 4>{{+1., -1., +0.5, -0.5}};
  // (a mesh size divisible by 16 would select a k-space slab decomposition
  // along x on 8 ranks, which leaves ranks without k-space data when the
  // r2c dimension is x)
  auto const setup = P3MTestMesh{{12, 14, 10}, 4};
  // reference forces of the padded halo exchange path (upstream, 1 rank)
  auto const forces_ref = std::array<Utils::Vector3d, 4>{{
      {8.821271830162e-04, -4.067718746322e-03, -1.333226550091e-03},
      {-4.657770983589e-03, 4.042236216161e-03, 2.805867909286e-03},
      {2.009429791044e-03, 2.511179681604e-03, 2.752882713058e-03},
      {1.766214009529e-03, -2.485697151442e-03, -4.225524072253e-03},
  }};
  setup_force_test(*this, positions, charges);
  test_all_p3m_fft_configs<>(setup,
                             [&](double) { check_forces(forces_ref, 1e-10); });
  // single precision, default FFT configuration (differs from the double
  // precision reference by about 1e-6 of the largest force)
  using FFTConfig = P3MFFTConfig<Utils::MemoryOrder::ROW_MAJOR,
                                 Utils::MemoryOrder::ROW_MAJOR, true, 2>;
  run_p3m_config<float, FFTConfig>(
      setup, [&](double) { check_forces(forces_ref, 1e-5); });
}
#endif // ESPRESSO_P3M

BOOST_AUTO_TEST_SUITE_END()

#endif // defined(ESPRESSO_P3M) or defined(ESPRESSO_DP3M)
