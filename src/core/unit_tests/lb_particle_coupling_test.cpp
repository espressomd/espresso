/*
 * Copyright (C) 2019 The ESPResSo project
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

/* Unit tests for LB particle coupling. */

#define BOOST_TEST_MODULE LB particle coupling test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
namespace bdata = boost::unit_test::data;
namespace utf = boost::unit_test;

#include "config.hpp"

#ifdef LB_WALBERLA

#include "ParticleFactory.hpp"

#include "EspressoSystemStandAlone.hpp"
#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include "random.hpp"
#include "thermostat.hpp"

#include <LBWalberlaBase.hpp>
#include <LatticeWalberla.hpp>

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// multiply by 100 because BOOST_CHECK_CLOSE takes a percentage tolerance,
// and by 6 to account for error accumulation
constexpr auto tol = 6 * 100 * std::numeric_limits<double>::epsilon();

struct LBTestParameters {
  unsigned int seed;
  double kT;
  double viscosity;
  double density;
  double tau;
  double time_step;
  double agrid;
  double skin;
  Utils::Vector3d box_dimensions;
  Utils::Vector3i grid_dimensions;
  Utils::Vector3d force_md_to_lb(Utils::Vector3d const &md_force) {
    return (-this->time_step * this->tau / this->agrid) * md_force;
  }
};

LBTestParameters params{23u,
                        0.,
                        1e-3,
                        0.5,
                        0.01,
                        0.01,
                        1.,
                        0.5,
                        Utils::Vector3d::broadcast(8.),
                        Utils::Vector3i::broadcast(8)};

/** Boost unit test dataset */
std::vector<double> const kTs{0, 1E-4};

namespace espresso {
// ESPResSo system instance
static std::unique_ptr<EspressoSystemStandAlone> system;
// ESPResSo actors
static std::shared_ptr<LBWalberlaParams> lb_params;
static std::shared_ptr<LatticeWalberla> lb_lattice;
static std::shared_ptr<LBWalberlaBase> lb_fluid;

void add_lb_actor_local(double kT) {
  constexpr unsigned int n_ghost_layers = 1u;
  params.kT = kT;
  lb_params = std::make_shared<LBWalberlaParams>(params.agrid, params.tau);
  lb_lattice = std::make_shared<LatticeWalberla>(params.grid_dimensions,
                                                 node_grid, n_ghost_layers);
  lb_fluid =
      mpi_init_lb_walberla_local(*lb_lattice, *lb_params, params.viscosity,
                                 params.density, params.kT, params.seed, false);
  mpi_activate_lb_walberla_local(lb_fluid, lb_params);
}

void set_lb_kT_local(double kT) {
  lb_fluid->set_collision_model(kT, params.seed);
}

void remove_lb_actor_local() { mpi_deactivate_lb_walberla_local(); }

void cells_update_ghosts_local() { cells_update_ghosts(global_ghost_flags()); }

REGISTER_CALLBACK(add_lb_actor_local)
REGISTER_CALLBACK(set_lb_kT_local)
REGISTER_CALLBACK(remove_lb_actor_local)
REGISTER_CALLBACK(cells_update_ghosts_local)
} // namespace espresso

/** Decorator to run a unit test only on the head node. */
struct if_head_node {
  boost::test_tools::assertion_result operator()(utf::test_unit_id) {
    return world.rank() == 0;
  }

private:
  boost::mpi::communicator world;
};

/** Fixture to manage the lifetime of the LB actor. */
struct CleanupActorLB : public ParticleFactory {
  CleanupActorLB() : ParticleFactory() {
    mpi_call_all(espresso::add_lb_actor_local, 0.);
  }
  ~CleanupActorLB() { mpi_call_all(espresso::remove_lb_actor_local); }
};

BOOST_FIXTURE_TEST_SUITE(suite, CleanupActorLB,
                         *utf::precondition(if_head_node()))

BOOST_AUTO_TEST_CASE(activate) {
  lb_lbcoupling_deactivate();
  mpi_bcast_lb_particle_coupling();
  lb_lbcoupling_activate();
  mpi_bcast_lb_particle_coupling();
  BOOST_CHECK(lb_particle_coupling.couple_to_md);
}

BOOST_AUTO_TEST_CASE(de_activate) {
  lb_lbcoupling_activate();
  mpi_bcast_lb_particle_coupling();
  lb_lbcoupling_deactivate();
  mpi_bcast_lb_particle_coupling();
  BOOST_CHECK(not lb_particle_coupling.couple_to_md);
}

BOOST_AUTO_TEST_CASE(rng_initial_state) {
  BOOST_CHECK(lb_lbcoupling_is_seed_required());
  BOOST_CHECK(!lb_particle_coupling.rng_counter_coupling);
}

BOOST_AUTO_TEST_CASE(rng) {
  lb_lbcoupling_set_rng_state(17);
  BOOST_REQUIRE(lb_particle_coupling.rng_counter_coupling);
  BOOST_CHECK_EQUAL(lb_lbcoupling_get_rng_state(), 17);
  BOOST_CHECK(not lb_lbcoupling_is_seed_required());
  auto step1_random1 = lb_particle_coupling_noise(
      true, 1, lb_particle_coupling.rng_counter_coupling);
  auto step1_random2 = lb_particle_coupling_noise(
      true, 4, lb_particle_coupling.rng_counter_coupling);
  auto step1_random2_try2 = lb_particle_coupling_noise(
      true, 4, lb_particle_coupling.rng_counter_coupling);
  BOOST_CHECK(step1_random1 != step1_random2);
  BOOST_CHECK(step1_random2 == step1_random2_try2);

  // Propagation queries kT from LB, so LB needs to be initialized
  mpi_call_all(espresso::set_lb_kT_local, 1E-4);
  lb_lbcoupling_propagate();

  BOOST_REQUIRE(lb_particle_coupling.rng_counter_coupling);
  BOOST_CHECK_EQUAL(lb_lbcoupling_get_rng_state(), 18);
  auto step2_random1 = lb_particle_coupling_noise(
      true, 1, lb_particle_coupling.rng_counter_coupling);
  auto step2_random2 = lb_particle_coupling_noise(
      true, 4, lb_particle_coupling.rng_counter_coupling);
  BOOST_CHECK(step1_random1 != step2_random1);
  BOOST_CHECK(step1_random1 != step2_random2);

  auto step3_norandom = lb_particle_coupling_noise(
      false, 4, lb_particle_coupling.rng_counter_coupling);
  BOOST_CHECK(step3_norandom == Utils::Vector3d{});
}

BOOST_AUTO_TEST_CASE(drift_vel_offset) {
  Particle p{};
  BOOST_CHECK_EQUAL(lb_particle_coupling_drift_vel_offset(p).norm(), 0);
  Utils::Vector3d expected{};
#ifdef ENGINE
  p.p.swim.swimming = true;
  p.p.swim.v_swim = 2.;
  expected += p.p.swim.v_swim * p.r.calc_director();
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
  p.p.mu_E = Utils::Vector3d{-2, 1.5, 1};
  expected += p.p.mu_E;
#endif
  BOOST_CHECK_SMALL(
      (lb_particle_coupling_drift_vel_offset(p) - expected).norm(), tol);
}

BOOST_DATA_TEST_CASE(drag_force, bdata::make(kTs), kT) {
  mpi_call_all(espresso::set_lb_kT_local, kT);
  Particle p{};
  p.m.v = {-2.5, 1.5, 2};
  p.r.p = lb_walberla()->get_local_domain().first;
  lb_lbcoupling_set_gamma(0.2);
  Utils::Vector3d drift_offset{-1, 1, 1};

  // Drag force in quiescent fluid
  {
    auto const observed = lb_drag_force(p, drift_offset);
    const Utils::Vector3d expected{0.3, -0.1, -.2};
    BOOST_CHECK_SMALL((observed - expected).norm(), tol);
  }
}

#ifdef ENGINE
BOOST_DATA_TEST_CASE(swimmer_force, bdata::make(kTs), kT) {
  mpi_call_all(espresso::set_lb_kT_local, kT);
  auto const first_lb_node = lb_walberla()->get_local_domain().first;
  Particle p{};
  p.p.swim.swimming = true;
  p.p.swim.f_swim = 2.;
  p.p.swim.dipole_length = 3.;
  p.p.swim.push_pull = 1;
  p.r.p = first_lb_node + Utils::Vector3d::broadcast(0.5);

  auto const coupling_pos =
      p.r.p + Utils::Vector3d{0., 0., p.p.swim.dipole_length / params.agrid};

  // swimmer coupling
  {
    if (in_local_halo(p.r.p)) {
      add_swimmer_force(p, params.time_step, false);
    }
    if (in_local_halo(coupling_pos)) {
      auto const interpolated =
          lb_lbfluid_get_force_to_be_applied(coupling_pos);
      auto const expected =
          params.force_md_to_lb(Utils::Vector3d{0., 0., p.p.swim.f_swim});

      // interpolation happened on the expected LB cell
      BOOST_CHECK_SMALL((interpolated - expected).norm(), tol);
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
          if ((pos - coupling_pos).norm() < 1e-6)
            continue;
          if (in_local_halo(pos)) {
            auto const interpolated = lb_lbfluid_get_force_to_be_applied(pos);
            BOOST_CHECK_SMALL(interpolated.norm(), tol);
          }
        }
      }
    }
  }

  // remove force of the particle from the fluid
  {
    if (in_local_halo(coupling_pos)) {
      add_md_force(coupling_pos, -Utils::Vector3d{0., 0., p.p.swim.f_swim},
                   params.time_step);
      auto const reset = lb_lbfluid_get_force_to_be_applied(coupling_pos);
      BOOST_REQUIRE_SMALL(reset.norm(), tol);
    }
  }
}
#endif // ENGINE

BOOST_DATA_TEST_CASE(particle_coupling, bdata::make(kTs), kT) {
  mpi_call_all(espresso::set_lb_kT_local, kT);
  lb_lbcoupling_set_rng_state(17);
  auto const first_lb_node = lb_walberla()->get_local_domain().first;
  auto const gamma = 0.2;
  auto const noise =
      (kT > 0.) ? std::sqrt(24. * gamma * kT / params.time_step) : 0.0;
  auto &rng = lb_particle_coupling.rng_counter_coupling;
  Particle p{};
  Utils::Vector3d expected = noise * Random::noise_uniform<RNGSalt::PARTICLES>(
                                         rng->value(), 0, p.identity());
#ifdef ENGINE
  p.p.swim.swimming = true;
  p.p.swim.v_swim = 2.;
  p.p.swim.push_pull = 1;
  expected += gamma * p.p.swim.v_swim * p.r.calc_director();
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
  p.p.mu_E = Utils::Vector3d{-2, 1.5, 1};
  expected += gamma * p.p.mu_E;
#endif
  p.r.p = first_lb_node + Utils::Vector3d::broadcast(0.5);
  lb_lbcoupling_set_gamma(gamma);

  // coupling
  {
    if (in_local_halo(p.r.p)) {
      couple_particle(p, false, noise, rng, params.time_step, false);
      BOOST_CHECK_SMALL((p.f.f - expected).norm(), tol);

      auto const interpolated = lb_lbfluid_get_force_to_be_applied(p.r.p);
      BOOST_CHECK_SMALL((interpolated - params.force_md_to_lb(expected)).norm(),
                        tol);
    }
  }

  // remove force of the particle from the fluid
  {
    if (in_local_halo(p.r.p)) {
      add_md_force(p.r.p, -expected, params.time_step);
    }
  }
}

BOOST_DATA_TEST_CASE_F(CleanupActorLB, coupling_particle_lattice_ia,
                       bdata::make(kTs), kT) {
  mpi_call_all(espresso::set_lb_kT_local, kT);
  lb_lbcoupling_set_rng_state(17);
  auto const first_lb_node = lb_walberla()->get_local_domain().first;
  auto const gamma = 0.2;
  auto const noise = std::sqrt(24. * gamma * kT / params.time_step *
                               Utils::sqr(params.agrid / params.tau));
  auto &rng = lb_particle_coupling.rng_counter_coupling;

  auto const pid = 0;
  auto const skin = params.skin;
  auto const &box_l = params.box_dimensions;
  create_particle({box_l[0] / 2. - skin * 2., skin * 2., skin * 2.});

  // sanity checks
  BOOST_REQUIRE_EQUAL(get_particle_node(pid), 0);
  BOOST_REQUIRE_EQUAL(ErrorHandling::mpi_gather_runtime_errors().size(), 0);

  auto const &p = get_particle_data(pid);
  auto expected =
      noise * Random::noise_uniform<RNGSalt::PARTICLES>(rng->value(), 0, pid);
#ifdef ENGINE
  set_particle_swimming(pid, ParticleParametersSwimming{true, 0., 2., 1, 3.});
  expected += gamma * p.p.swim.v_swim * p.r.calc_director();
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
  set_particle_mu_E(pid, {-2., 1.5, 1.});
  expected += gamma * p.p.mu_E;
#endif
  place_particle(pid, first_lb_node + Utils::Vector3d::broadcast(0.5));
  lb_lbcoupling_set_gamma(gamma);

  for (bool with_ghosts : {false, true}) {
    {
      if (with_ghosts) {
        mpi_call_all(espresso::cells_update_ghosts_local);
      }
      auto const particles = cell_structure.local_particles();
      auto const ghost_particles = cell_structure.ghost_particles();
      BOOST_REQUIRE_GE(particles.size(), 1);
      BOOST_REQUIRE_GE(ghost_particles.size(), with_ghosts);
    }

    // check box shifts
    {
      auto constexpr reference_shifts =
          std::array<Utils::Vector3i, 7>{{{{0, 0, 8}},
                                          {{0, 8, 0}},
                                          {{0, 8, 8}},
                                          {{8, 0, 0}},
                                          {{8, 0, 8}},
                                          {{8, 8, 0}},
                                          {{8, 8, 8}}}};
      boost::mpi::communicator world;
      assert(world.size() <= 4);
      auto const cutoff = 8 / world.size() - 1;
      {
        auto const shifts = shifted_positions({0., 0., 0.}, box_geo);
        BOOST_REQUIRE_EQUAL(shifts.size(), cutoff);
        for (std::size_t i = 0; i < shifts.size(); ++i) {
          BOOST_REQUIRE_EQUAL(shifts[i], reference_shifts[i]);
        }
      }
      {
        auto const shifts = shifted_positions({1., 1., 1.}, box_geo);
        BOOST_REQUIRE_EQUAL(shifts.size(), 0);
      }
      {
        auto const reference_shift = Utils::Vector3d{1., 2., 8.};
        auto const shifts = shifted_positions({1., 2., 0.}, box_geo);
        BOOST_REQUIRE_EQUAL(shifts.size(), 1);
        BOOST_REQUIRE_EQUAL(shifts[0], reference_shift);
      }
    }

    // check without LB coupling
    {
      lb_lbcoupling_deactivate();
      mpi_bcast_lb_particle_coupling();
      auto const particles = cell_structure.local_particles();
      auto const ghost_particles = cell_structure.ghost_particles();
      lb_lbcoupling_calc_particle_lattice_ia(thermo_virtual, particles,
                                             ghost_particles, params.time_step);
      auto const &p = get_particle_data(pid);
      BOOST_CHECK_EQUAL(p.f.f.norm(), 0.);
    }

    // check with LB coupling
    {
      lb_lbcoupling_activate();
      mpi_bcast_lb_particle_coupling();
      auto const particles = cell_structure.local_particles();
      auto const ghost_particles = cell_structure.ghost_particles();
      // get original LB force
      auto const lb_before = lb_lbfluid_get_force_to_be_applied(p.r.p);
      // couple particle to LB
      lb_lbcoupling_calc_particle_lattice_ia(thermo_virtual, particles,
                                             ghost_particles, params.time_step);
      // check particle force
      auto const &p = get_particle_data(pid);
      BOOST_CHECK_SMALL((p.f.f - expected).norm(), tol);
      // check LB force
      auto const lb_after = lb_lbfluid_get_force_to_be_applied(p.r.p);
      auto const lb_expected = params.force_md_to_lb(expected) + lb_before;
      BOOST_CHECK_SMALL((lb_after - lb_expected).norm(), tol);
      // remove force of the particle from the fluid
      set_particle_f(pid, {});
      add_md_force(p.r.p, -expected, params.time_step);
    }
  }

  // clean-up and sanity checks
  {
    auto const error_message_ref = std::string(
        "Recalculating forces, so the LB coupling forces are not included in "
        "the particle force the first time step. This only matters if it "
        "happens frequently during sampling.");
    auto const error_messages = ErrorHandling::mpi_gather_runtime_errors();
    for (auto const &error_message : error_messages) {
      BOOST_CHECK_EQUAL(error_message.what(), error_message_ref);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

int main(int argc, char **argv) {
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);
  espresso::system->set_box_l(params.box_dimensions);
  espresso::system->set_time_step(params.time_step);
  espresso::system->set_skin(params.skin);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
