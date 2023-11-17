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

#define BOOST_TEST_MODULE LB particle coupling test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
namespace bdata = boost::unit_test::data;
namespace utf = boost::unit_test;

#include "config/config.hpp"

#ifdef WALBERLA

#include "ParticleFactory.hpp"
#include "particle_management.hpp"

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "global_ghost_flags.hpp"
#include "lb/LBNone.hpp"
#include "lb/LBWalberla.hpp"
#include "lb/particle_coupling.hpp"
#include "particle_node.hpp"
#include "random.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi.hpp>
#include <boost/serialization/optional.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// multiply by 6 to account for error accumulation
auto constexpr eps = 6. * std::numeric_limits<double>::epsilon();

static struct {
  unsigned int seed = 23u;
  double kT = 0.;
  double viscosity = 1e-3;
  double density = 0.5;
  double tau = 0.01;
  double time_step = 0.01;
  double agrid = 1.;
  double skin = 0.5;
  Utils::Vector3d box_dimensions = Utils::Vector3d::broadcast(8.);
  Utils::Vector3i grid_dimensions = Utils::Vector3i::broadcast(8);
  auto force_md_to_lb(Utils::Vector3d const &md_force) const {
    return (this->time_step * this->tau / this->agrid) * md_force;
  }
} params;

/** Boost unit test dataset */
std::vector<double> const kTs{0., 1E-4};

namespace espresso {
// ESPResSo system instance
static std::shared_ptr<System::System> system;
// ESPResSo actors
static std::shared_ptr<LB::LBWalberlaParams> lb_params;
static std::shared_ptr<LatticeWalberla> lb_lattice;
static std::shared_ptr<LBWalberlaBase> lb_fluid;

static auto make_lb_actor() {
  auto constexpr n_ghost_layers = 1u;
  auto constexpr single_precision = false;
  lb_params = std::make_shared<LB::LBWalberlaParams>(params.agrid, params.tau);
  lb_lattice = std::make_shared<LatticeWalberla>(
      params.grid_dimensions, ::communicator.node_grid, n_ghost_layers);
  lb_fluid = new_lb_walberla(lb_lattice, params.viscosity, params.density,
                             single_precision);
  lb_fluid->set_collision_model(params.kT, params.seed);
  lb_fluid->ghost_communication();
}

static void add_lb_actor() {
  espresso::system->lb.set<::LB::LBWalberla>(lb_fluid, lb_params);
}

static void remove_lb_actor() { espresso::system->lb.reset(); }

static void set_lb_kT(double kT) {
  lb_fluid->set_collision_model(kT, params.seed);
}
} // namespace espresso

namespace LB {
static auto get_force_to_be_applied(Utils::Vector3d const &pos) {
  auto const agrid = espresso::lb_params->get_agrid();
  auto const ind = Utils::Vector3i{static_cast<int>(pos[0] / agrid),
                                   static_cast<int>(pos[1] / agrid),
                                   static_cast<int>(pos[2] / agrid)};
  auto const res = espresso::lb_fluid->get_node_force_to_be_applied(ind);
  if (!res) {
    auto const comm = boost::mpi::communicator();
    std::stringstream what;
    what << "Force to be applied could not be obtained from Walberla "
         << "on MPI rank " << comm.rank() << ": position = [" << pos << "]";
    throw std::runtime_error(what.str());
  }
  return *res;
}
} // namespace LB

/** Fixture to manage the lifetime of the LB actor. */
struct CleanupActorLB : public ParticleFactory {
  CleanupActorLB() : ParticleFactory() {
    params.kT = 0.;
    espresso::make_lb_actor();
    espresso::add_lb_actor();
  }

  // NOLINTNEXTLINE(bugprone-exception-escape)
  ~CleanupActorLB() { espresso::remove_lb_actor(); }
};

BOOST_FIXTURE_TEST_SUITE(suite, CleanupActorLB)

static void lb_lbcoupling_broadcast() {
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, lb_particle_coupling, 0);
}

BOOST_AUTO_TEST_CASE(activate) {
  lb_lbcoupling_deactivate();
  lb_lbcoupling_broadcast();
  lb_lbcoupling_activate();
  lb_lbcoupling_broadcast();
  BOOST_CHECK(lb_particle_coupling.couple_to_md);
}

BOOST_AUTO_TEST_CASE(de_activate) {
  lb_lbcoupling_activate();
  lb_lbcoupling_broadcast();
  lb_lbcoupling_deactivate();
  lb_lbcoupling_broadcast();
  BOOST_CHECK(not lb_particle_coupling.couple_to_md);
}

BOOST_AUTO_TEST_CASE(rng_initial_state) {
  BOOST_CHECK(lb_lbcoupling_is_seed_required());
  BOOST_CHECK(!lb_particle_coupling.rng_counter_coupling);
}

BOOST_AUTO_TEST_CASE(rng) {
  lb_lbcoupling_set_rng_state(17);
  lb_particle_coupling.gamma = 0.2;
  auto &lb = espresso::system->lb;

  LB::ParticleCoupling coupling{lb, true, params.time_step, 1.};
  BOOST_REQUIRE(lb_particle_coupling.rng_counter_coupling);
  BOOST_CHECK_EQUAL(lb_lbcoupling_get_rng_state(), 17);
  BOOST_CHECK(not lb_lbcoupling_is_seed_required());
  Particle test_partcl_1{};
  test_partcl_1.id() = 1;
  Particle test_partcl_2{};
  test_partcl_2.id() = 4;
  auto const step1_random1 = coupling.get_noise_term(test_partcl_1);
  auto const step1_random2 = coupling.get_noise_term(test_partcl_2);
  auto const step1_random2_try2 = coupling.get_noise_term(test_partcl_2);
  BOOST_REQUIRE(step1_random1.norm() > 1e-10);
  BOOST_CHECK(step1_random1 != step1_random2);
  BOOST_CHECK(step1_random2 == step1_random2_try2);

  // Propagation queries kT from LB, so LB needs to be initialized
  espresso::set_lb_kT(1E-4);
  lb_lbcoupling_propagate();

  BOOST_REQUIRE(lb_particle_coupling.rng_counter_coupling);
  BOOST_CHECK_EQUAL(lb_lbcoupling_get_rng_state(), 18);
  auto const step2_random1 = coupling.get_noise_term(test_partcl_1);
  auto const step2_random2 = coupling.get_noise_term(test_partcl_2);
  BOOST_CHECK(step1_random1 != step2_random1);
  BOOST_CHECK(step1_random1 != step2_random2);

  LB::ParticleCoupling coupling_unthermalized{lb, true, params.time_step, 0.};
  auto const step3_norandom =
      coupling_unthermalized.get_noise_term(test_partcl_2);
  BOOST_CHECK((step3_norandom == Utils::Vector3d{0., 0., 0.}));
}

BOOST_AUTO_TEST_CASE(drift_vel_offset) {
  Particle p{};
  auto &lb = espresso::system->lb;
  LB::ParticleCoupling coupling{lb, false, params.time_step};
  BOOST_CHECK_EQUAL(coupling.lb_drift_velocity_offset(p).norm(), 0);
  Utils::Vector3d expected{};
#ifdef LB_ELECTROHYDRODYNAMICS
  p.mu_E() = Utils::Vector3d{-2., 1.5, 1.};
  expected += p.mu_E();
#endif
  BOOST_CHECK_SMALL((coupling.lb_drift_velocity_offset(p) - expected).norm(),
                    eps);
}

BOOST_DATA_TEST_CASE(drag_force, bdata::make(kTs), kT) {
  espresso::set_lb_kT(kT);
  auto &lb = espresso::system->lb;
  Particle p{};
  p.v() = {-2.5, 1.5, 2.};
  p.pos() = espresso::lb_fluid->get_lattice().get_local_domain().first;
  lb_lbcoupling_set_gamma(0.2);
  Utils::Vector3d drift_offset{-1., 1., 1.};

  // Drag force in quiescent fluid
  {
    auto const observed = lb_drag_force(lb, p, p.pos(), drift_offset);
    const Utils::Vector3d expected{0.3, -0.1, -.2};
    BOOST_CHECK_SMALL((observed - expected).norm(), eps);
  }
}

#ifdef ENGINE
BOOST_DATA_TEST_CASE(swimmer_force, bdata::make(kTs), kT) {
  espresso::set_lb_kT(kT);
  auto &lb = espresso::system->lb;
  auto const first_lb_node =
      espresso::lb_fluid->get_lattice().get_local_domain().first;
  Particle p{};
  p.swimming().swimming = true;
  p.swimming().f_swim = 2.;
  p.swimming().is_engine_force_on_fluid = true;
  p.pos() = first_lb_node + Utils::Vector3d::broadcast(0.5);

  // swimmer coupling
  {
    if (in_local_halo(p.pos())) {
      LB::ParticleCoupling coupling{lb, true, params.time_step};
      coupling.kernel(p);
      auto const interpolated = LB::get_force_to_be_applied(p.pos());
      auto const expected =
          params.force_md_to_lb(Utils::Vector3d{0., 0., p.swimming().f_swim});

      // interpolation happened on the expected LB cell
      BOOST_CHECK_SMALL((interpolated - expected).norm(), eps);
    }

    // all other LB cells have no force
    for (int i = 0; i < params.grid_dimensions[0]; ++i) {
      for (int j = 0; j < params.grid_dimensions[1]; ++j) {
        for (int k = 0; k < params.grid_dimensions[2]; ++k) {
          auto const pos = Utils::Vector3d{
              0.5 + static_cast<double>(i) * params.agrid,
              0.5 + static_cast<double>(j) * params.agrid,
              0.5 + static_cast<double>(k) * params.agrid,
          };
          if ((pos - p.pos()).norm() < 1e-6)
            continue;
          if (in_local_halo(pos)) {
            auto const interpolated = LB::get_force_to_be_applied(pos);
            BOOST_CHECK_SMALL(interpolated.norm(), eps);
          }
        }
      }
    }
  }

  // remove force of the particle from the fluid
  {
    if (in_local_halo(p.pos())) {
      add_md_force(lb, p.pos(), -Utils::Vector3d{0., 0., p.swimming().f_swim},
                   params.time_step);
      auto const reset = LB::get_force_to_be_applied(p.pos());
      BOOST_REQUIRE_SMALL(reset.norm(), eps);
    }
  }
}
#endif // ENGINE

BOOST_DATA_TEST_CASE(particle_coupling, bdata::make(kTs), kT) {
  espresso::set_lb_kT(kT);
  lb_lbcoupling_set_rng_state(17);
  auto &lb = espresso::system->lb;
  auto const first_lb_node =
      espresso::lb_fluid->get_lattice().get_local_domain().first;
  auto const gamma = 0.2;
  lb_lbcoupling_set_gamma(gamma);
  Particle p{};
  LB::ParticleCoupling coupling{lb, false, params.time_step};
  auto expected = coupling.get_noise_term(p);
#ifdef LB_ELECTROHYDRODYNAMICS
  p.mu_E() = Utils::Vector3d{-2., 1.5, 1.};
  expected += gamma * p.mu_E();
#endif
  p.pos() = first_lb_node + Utils::Vector3d::broadcast(0.5);

  // coupling
  {
    if (in_local_halo(p.pos())) {
      coupling.kernel(p);
      BOOST_CHECK_SMALL((p.force() - expected).norm(), eps);

      auto const interpolated = -LB::get_force_to_be_applied(p.pos());
      BOOST_CHECK_SMALL((interpolated - params.force_md_to_lb(expected)).norm(),
                        eps);
    }
  }

  // remove force of the particle from the fluid
  {
    if (in_local_halo(p.pos())) {
      add_md_force(lb, p.pos(), -expected, params.time_step);
    }
  }
}

BOOST_DATA_TEST_CASE_F(CleanupActorLB, coupling_particle_lattice_ia,
                       bdata::make(kTs), kT) {
  auto const comm = boost::mpi::communicator();
  auto const rank = comm.rank();
  espresso::set_lb_kT(kT);
  lb_lbcoupling_set_rng_state(17);
  auto &system = *espresso::system;
  auto &cell_structure = *system.cell_structure;
  auto &lb = system.lb;
  auto const &box_geo = *system.box_geo;
  auto const first_lb_node =
      espresso::lb_fluid->get_lattice().get_local_domain().first;
  auto const gamma = 0.2;
  auto const pid = 0;
  auto const skin = params.skin;
  auto const &box_l = params.box_dimensions;
  create_particle({box_l[0] / 2. - skin * 2., skin * 2., skin * 2.}, 0, 0);

  // sanity checks
  BOOST_REQUIRE_EQUAL(get_particle_node_parallel(pid), rank ? -1 : 0);
  BOOST_REQUIRE_EQUAL(
      ErrorHandling::mpi_gather_runtime_errors_all(rank == 0).size(), 0);

#ifdef ENGINE
  set_particle_property(pid, &Particle::swimming,
                        ParticleParametersSwimming{2., true, false});
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
  set_particle_property(pid, &Particle::mu_E, Utils::Vector3d{-2., 1.5, 1.});
#endif

  LB::ParticleCoupling coupling{lb, thermo_virtual, params.time_step};
  auto const p_opt = copy_particle_to_head_node(comm, system, pid);
  auto expected = Utils::Vector3d{};
  if (rank == 0) {
    auto const &p = *p_opt;
    expected += coupling.get_noise_term(p);
#ifdef LB_ELECTROHYDRODYNAMICS
    expected += gamma * p.mu_E();
#endif
  }
  boost::mpi::broadcast(comm, expected, 0);
  auto const p_pos = first_lb_node + Utils::Vector3d::broadcast(0.5);
  set_particle_pos(pid, p_pos);
  lb_lbcoupling_set_gamma(gamma);

  for (bool with_ghosts : {false, true}) {
    {
      if (with_ghosts) {
        cell_structure.update_ghosts_and_resort_particle(global_ghost_flags());
      }
      if (rank == 0) {
        auto const particles = cell_structure.local_particles();
        auto const ghost_particles = cell_structure.ghost_particles();
        BOOST_REQUIRE_GE(particles.size(), 1);
        BOOST_REQUIRE_GE(ghost_particles.size(), static_cast<int>(with_ghosts));
      }
    }

    // check box shifts
    if (rank == 0) {
      auto constexpr reference_shifts =
          std::array<Utils::Vector3i, 8>{{{{0, 0, 0}},
                                          {{0, 0, 8}},
                                          {{0, 8, 0}},
                                          {{0, 8, 8}},
                                          {{8, 0, 0}},
                                          {{8, 0, 8}},
                                          {{8, 8, 0}},
                                          {{8, 8, 8}}}};
      boost::mpi::communicator world;
      assert(world.size() <= 4);
      auto const cutoff = 8 / world.size();
      auto const agrid = params.agrid;
      {
        auto const shifts = positions_in_halo({0., 0., 0.}, box_geo, agrid);
        BOOST_REQUIRE_EQUAL(shifts.size(), cutoff);
        for (std::size_t i = 0; i < shifts.size(); ++i) {
          BOOST_REQUIRE_EQUAL(shifts[i], reference_shifts[i]);
        }
      }
      {
        auto const reference_shift = Utils::Vector3d{1., 1., 1.};
        auto const shifts = positions_in_halo({1., 1., 1.}, box_geo, agrid);
        BOOST_REQUIRE_EQUAL(shifts.size(), 1);
        BOOST_REQUIRE_EQUAL(shifts[0], reference_shift);
      }
      {
        auto const reference_origin = Utils::Vector3d{1., 2., 0.};
        auto const reference_shift = Utils::Vector3d{1., 2., 8.};
        auto const shifts = positions_in_halo({1., 2., 0.}, box_geo, agrid);
        BOOST_REQUIRE_EQUAL(shifts.size(), 2);
        BOOST_REQUIRE_EQUAL(shifts[0], reference_origin);
        BOOST_REQUIRE_EQUAL(shifts[1], reference_shift);
      }
    }

    // check without LB coupling
    {
      lb_lbcoupling_deactivate();
      lb_lbcoupling_broadcast();
      auto const particles = cell_structure.local_particles();
      auto const ghost_particles = cell_structure.ghost_particles();
      LB::couple_particles(thermo_virtual, particles, ghost_particles,
                           params.time_step);
      auto const p_opt = copy_particle_to_head_node(comm, system, pid);
      if (rank == 0) {
        auto const &p = *p_opt;
        BOOST_CHECK_EQUAL(p.force().norm(), 0.);
      }
    }

    // check with LB coupling
    {
      lb_lbcoupling_activate();
      lb_lbcoupling_broadcast();
      auto const particles = cell_structure.local_particles();
      auto const ghost_particles = cell_structure.ghost_particles();
      Utils::Vector3d lb_before{};
      {
        auto const p_opt = copy_particle_to_head_node(comm, system, pid);
        if (rank == 0) {
          auto const &p = *p_opt;
          // get original LB force
          lb_before = -LB::get_force_to_be_applied(p.pos());
        }
      }
      // couple particle to LB
      LB::couple_particles(thermo_virtual, particles, ghost_particles,
                           params.time_step);
      {
        auto const p_opt = copy_particle_to_head_node(comm, system, pid);
        if (rank == 0) {
          auto const &p = *p_opt;
          // check particle force
          BOOST_CHECK_SMALL((p.force() - expected).norm(), eps);
          // check LB force
          auto const lb_after = -LB::get_force_to_be_applied(p.pos());
          auto const lb_expected = params.force_md_to_lb(expected) + lb_before;
          BOOST_CHECK_SMALL((lb_after - lb_expected).norm(), eps);
        }
      }
      // remove force of the particle from the fluid
      set_particle_property(pid, &Particle::force, Utils::Vector3d{});
      add_md_force(lb, p_pos, -expected, params.time_step);
    }
  }

  // clean-up and sanity checks
  {
    boost::mpi::communicator world;
    auto const error_message_ref = std::string(
        "Recalculating forces, so the LB coupling forces are not included in "
        "the particle force the first time step. This only matters if it "
        "happens frequently during sampling.");
    auto const error_messages =
        ErrorHandling::mpi_gather_runtime_errors_all(world.rank() == 0);
    for (auto const &error_message : error_messages) {
      BOOST_CHECK_EQUAL(error_message.what(), error_message_ref);
    }
    lb_particle_coupling.rng_counter_coupling = {};
    if (kT == 0.) {
      BOOST_CHECK_EQUAL(coupling.get_noise_term(Particle{}).norm(), 0.);
    } else {
      BOOST_CHECK_THROW(coupling.get_noise_term(Particle{}),
                        std::runtime_error);
    }
  }
}

BOOST_AUTO_TEST_CASE(runtime_exceptions) {
  boost::mpi::communicator world;
  auto &lb = espresso::system->lb;
  // LB prevents changing most of the system state
  {
    BOOST_CHECK_THROW(lb.on_boxl_change(), std::runtime_error);
    BOOST_CHECK_THROW(lb.on_timestep_change(), std::runtime_error);
    BOOST_CHECK_THROW(lb.on_temperature_change(), std::runtime_error);
    BOOST_CHECK_THROW(lb.on_node_grid_change(), std::runtime_error);
  }
  // check access out of the LB local domain
  {
    if (world.rank() != 0) {
      BOOST_CHECK(not lb.get_interpolated_density({0., 0., 0.}).has_value());
      BOOST_CHECK(not lb.get_interpolated_velocity({0., 0., 0.}).has_value());
      BOOST_CHECK_THROW(lb.add_force_density({0., 0., 0.}, {0., 0., 0.}),
                        std::runtime_error);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

bool test_lb_domain_mismatch_local() {
  boost::mpi::communicator world;
  auto const node_grid_original = ::communicator.node_grid;
  auto const node_grid_reversed =
      Utils::Vector3i{{::communicator.node_grid[2], ::communicator.node_grid[1],
                       ::communicator.node_grid[0]}};
  auto const n_ghost_layers = 1u;
  auto const params = std::make_shared<LB::LBWalberlaParams>(0.5, 0.01);
  ::communicator.node_grid = node_grid_reversed;
  auto const lattice = std::make_shared<LatticeWalberla>(
      Utils::Vector3i{12, 12, 12}, node_grid_original, n_ghost_layers);
  auto const ptr = new_lb_walberla(lattice, 1.0, 1.0, false);
  ptr->set_collision_model(0.0, 0);
  ::communicator.node_grid = node_grid_original;
  auto lb_instance = std::make_shared<LB::LBWalberla>(ptr, params);
  if (world.rank() == 0) {
    try {
      lb_instance->sanity_checks(*espresso::system);
    } catch (std::runtime_error const &err) {
      auto const what_ref = std::string("waLBerla and ESPResSo disagree "
                                        "about domain decomposition.");
      return err.what() == what_ref;
    }
  }
  return false;
}

BOOST_AUTO_TEST_CASE(lb_exceptions) {
  boost::mpi::communicator world;
  auto &lb = espresso::system->lb;
  // LB exceptions mechanism
  {
    using std::exception;
    // getters and setters
    BOOST_CHECK_THROW(lb.get_agrid(), exception);
    BOOST_CHECK_THROW(lb.get_tau(), exception);
    BOOST_CHECK_THROW(lb.get_kT(), exception);
    BOOST_CHECK_THROW(lb.get_pressure_tensor(), exception);
    BOOST_CHECK_THROW(LB::get_force_to_be_applied({-10., -10., -10.}),
                      std::runtime_error);
    // coupling, interpolation, boundaries
    BOOST_CHECK_THROW(lb_lbcoupling_get_rng_state(), std::runtime_error);
    BOOST_CHECK_THROW(lb_lbcoupling_set_rng_state(0ul), std::runtime_error);
    BOOST_CHECK_THROW(lb.get_momentum(), exception);
  }

  // waLBerla and ESPResSo must agree on domain decomposition
  {
    auto const has_thrown_correct_exception = test_lb_domain_mismatch_local();
    auto const n_errors = check_runtime_errors_local();
    auto const error_queue =
        ErrorHandling::mpi_gather_runtime_errors_all(world.rank() == 0);
    if (world.rank() == 0) {
      BOOST_TEST_REQUIRE(has_thrown_correct_exception);
      BOOST_REQUIRE_EQUAL(n_errors, 1);
      BOOST_REQUIRE_EQUAL(error_queue.size(), 1);
      auto const what_ref = std::string("MPI rank 0: left ESPResSo: [0, 0, 0], "
                                        "left waLBerla: [0, 0, 0]");
      for (auto const &error : error_queue) {
        auto const error_what = error.what().substr(1, what_ref.size());
        BOOST_CHECK_EQUAL(error_what, what_ref);
      }
    } else {
      BOOST_TEST_REQUIRE(not has_thrown_correct_exception);
    }
  }

  // LBNone exceptions mechanism
  {
    using LB::NoLBActive;
    auto const vec = Utils::Vector3d{};
    auto lb_impl = std::make_shared<LB::LBNone>();
    lb.set<LB::LBNone>(lb_impl);
    BOOST_CHECK_THROW(lb.get_agrid(), NoLBActive);
    BOOST_CHECK_THROW(lb.get_tau(), NoLBActive);
    BOOST_CHECK_THROW(lb.get_kT(), NoLBActive);
    BOOST_CHECK_THROW(lb.get_pressure_tensor(), NoLBActive);
    BOOST_CHECK_THROW(lb.get_momentum(), NoLBActive);
    BOOST_CHECK_THROW(lb.sanity_checks(), NoLBActive);
    BOOST_CHECK_THROW(lb.veto_time_step(0.), NoLBActive);
    BOOST_CHECK_THROW(lb.lebc_sanity_checks(0u, 1u), NoLBActive);
    BOOST_CHECK_THROW(lb.propagate(), NoLBActive);
    BOOST_CHECK_THROW(lb.on_cell_structure_change(), NoLBActive);
    BOOST_CHECK_THROW(lb.on_boxl_change(), NoLBActive);
    BOOST_CHECK_THROW(lb.on_node_grid_change(), NoLBActive);
    BOOST_CHECK_THROW(lb.on_timestep_change(), NoLBActive);
    BOOST_CHECK_THROW(lb.on_temperature_change(), NoLBActive);
    BOOST_CHECK_THROW(lb_impl->get_density_at_pos(vec, true), NoLBActive);
    BOOST_CHECK_THROW(lb_impl->get_velocity_at_pos(vec, true), NoLBActive);
    BOOST_CHECK_THROW(lb_impl->add_force_at_pos(vec, vec), NoLBActive);
    lb.reset();
  }
}

int main(int argc, char **argv) {
  mpi_init_stand_alone(argc, argv);
  espresso::system = System::System::create();
  espresso::system->set_box_l(params.box_dimensions);
  espresso::system->set_time_step(params.time_step);
  espresso::system->set_cell_structure_topology(CellStructureType::REGULAR);
  espresso::system->cell_structure->set_verlet_skin(params.skin);
  ::System::set_system(espresso::system);

  boost::mpi::communicator world;
  assert(world.size() <= 2);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

#else // WALBERLA
int main(int argc, char **argv) {}
#endif
