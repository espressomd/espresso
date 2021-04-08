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

#include "config.hpp"

#ifdef LB_WALBERLA

#include <lb_walberla_init.hpp>

#include "Particle.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include "random.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <limits>
#include <vector>

// multiply by 100 because BOOST_CHECK_CLOSE takes a percentage tolerance,
// and by 6 to account for error accumulation
constexpr auto tol = 6 * 100 * std::numeric_limits<double>::epsilon();

void setup_lb(double kT) {
  init_lb_walberla(1E-3,                           // viscosity
                   0.5,                            // density
                   1,                              // agrid
                   0.01,                           // tau
                   Utils::Vector3i::broadcast(12), // grid
                   node_grid,
                   kT,  // kT
                   23); // seed
}

BOOST_AUTO_TEST_CASE(activate) {
  lb_lbcoupling_deactivate();
  lb_lbcoupling_activate();
  BOOST_CHECK(lb_particle_coupling.couple_to_md);
}
BOOST_AUTO_TEST_CASE(de_activate) {
  lb_lbcoupling_activate();
  lb_lbcoupling_deactivate();
  BOOST_CHECK(not lb_particle_coupling.couple_to_md);
}

BOOST_AUTO_TEST_CASE(rng_initial_state) {
  lattice_switch = ActiveLB::WALBERLA;
  BOOST_CHECK(lb_lbcoupling_is_seed_required());
  BOOST_CHECK_THROW(lb_lbcoupling_get_rng_state(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(rng) {
  lattice_switch = ActiveLB::WALBERLA;
  lb_lbcoupling_set_rng_state(17);
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

  // Propagation queries kT from lb, so lb needs to be initialized
  setup_lb(1E-4); // kT
  lb_lbcoupling_propagate();

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
#ifdef ENGINE
  p.p.swim.swimming = true;
  p.p.swim.v_swim = 2;
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
  p.p.mu_E = Utils::Vector3d{-2, 1.5, 1};
#endif
  Utils::Vector3d expected =
      Utils::Vector3d{-2, 1.5, 1} + 2 * p.r.calc_director();
  BOOST_CHECK_SMALL(
      (lb_particle_coupling_drift_vel_offset(p) - expected).norm(), tol);
}
const std::vector<double> kTs{0, 1E-4};

BOOST_DATA_TEST_CASE(drag_force, bdata::make(kTs), kT) {
  setup_lb(kT); // in LB units.
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

// todo: test remaining functionality from lb_particle_coupling.hpp

int main(int argc, char **argv) {
  auto mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  Communication::init(mpi_env);
  walberla_mpi_init();

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
