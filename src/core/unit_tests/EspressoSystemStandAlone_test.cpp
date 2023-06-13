/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE EspressoSystemStandAlone test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
namespace utf = boost::unit_test;

#include "ParticleFactory.hpp"
#include "particle_management.hpp"

#include "EspressoSystemStandAlone.hpp"
#include "Particle.hpp"
#include "accumulators/TimeSeries.hpp"
#include "bonded_interactions/bonded_interaction_utils.hpp"
#include "bonded_interactions/fene.hpp"
#include "bonded_interactions/harmonic.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "electrostatics/p3m.hpp"
#include "electrostatics/registration.hpp"
#include "energy.hpp"
#include "event.hpp"
#include "galilei/Galilei.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/lj.hpp"
#include "observables/ParticleVelocities.hpp"
#include "particle_node.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/range/numeric.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace espresso {
// ESPResSo system instance
static std::unique_ptr<EspressoSystemStandAlone> system;
} // namespace espresso

static void remove_translational_motion() {
  Galilei{}.kill_particle_motion(false);
}

BOOST_FIXTURE_TEST_CASE(espresso_system_stand_alone, ParticleFactory) {
  auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();
  auto const comm = boost::mpi::communicator();
  auto const rank = comm.rank();
  auto const n_nodes = comm.size();

  auto const box_l = 12.;
  auto const box_center = box_l / 2.;
  auto const time_step = 0.001;
  auto const skin = 0.4;
  espresso::system->set_box_l(Utils::Vector3d::broadcast(box_l));
  espresso::system->set_time_step(time_step);
  espresso::system->set_skin(skin);

  // particle properties
  auto const pid1 = 9;
  auto const pid2 = 2;
  auto const pid3 = 5;
  auto const type_a = 1;
  auto const type_b = 2;

  // we need at least 2 MPI ranks to test the communication logic, therefore
  // place particles close to the interface between different MPI domains
  auto const start_positions = std::unordered_map<int, Utils::Vector3d>{
      {pid1, {box_center - 0.1, box_center - 0.1, 1.0}},
      {pid2, {box_center + 0.1, box_center - 0.1, 1.0}},
      {pid3, {box_center + 0.1, box_center + 0.1, 1.0}}};
  create_particle(start_positions.at(pid1), pid1, type_a);
  create_particle(start_positions.at(pid2), pid2, type_b);
  create_particle(start_positions.at(pid3), pid3, type_b);
  if (n_nodes % 2 == 0) {
    BOOST_REQUIRE_EQUAL(get_particle_node_parallel(pid1), rank ? -1 : 0);
    BOOST_REQUIRE_GE(get_particle_node_parallel(pid2), rank ? -1 : 1);
    BOOST_REQUIRE_GE(get_particle_node_parallel(pid3), rank ? -1 : 1);
  }

  auto const reset_particle_positions = [&start_positions]() {
    for (auto const &kv : start_positions) {
      set_particle_pos(kv.first, kv.second);
    }
  };

  // check accumulators
  if (n_nodes == 1) {
    auto const pids = std::vector<int>{pid2};
    auto obs = std::make_shared<Observables::ParticleVelocities>(pids);
    auto acc = Accumulators::TimeSeries(obs, 1);

    auto const obs_shape = obs->shape();
    auto const ref_shape = std::vector<std::size_t>{pids.size(), 3u};
    BOOST_REQUIRE_EQUAL_COLLECTIONS(obs_shape.begin(), obs_shape.end(),
                                    ref_shape.begin(), ref_shape.end());

    remove_translational_motion();
    for (int i = 0; i < 5; ++i) {
      set_particle_v(pid2, {static_cast<double>(i), 0., 0.});

      acc.update();
      auto const time_series = acc.time_series();
      BOOST_REQUIRE_EQUAL(time_series.size(), i + 1);

      auto const acc_value = time_series.back();
      auto const obs_value = (*obs)();
      auto const &p = get_particle_data(pid2);
      BOOST_TEST(obs_value == p.v(), boost::test_tools::per_element());
      BOOST_TEST(acc_value == p.v(), boost::test_tools::per_element());
    }
  }

  // check kinetic energy
  {
    remove_translational_motion();
    for (int i = 0; i < 5; ++i) {
      set_particle_v(pid2, {static_cast<double>(i), 0., 0.});
      auto const obs_energy = calculate_energy();
      auto const p_opt = copy_particle_to_head_node(comm, pid2);
      if (rank == 0) {
        auto const &p = *p_opt;
        auto const kinetic_energy = 0.5 * p.mass() * p.v().norm2();
        BOOST_CHECK_CLOSE(obs_energy->kinetic[0], kinetic_energy, tol);
      }
    }
  }

  // check non-bonded energies
#ifdef LENNARD_JONES
  {
    // distance between particles
    auto const dist = 0.2;
    // set up LJ potential
    auto const eps = 1.0;
    auto const sig = 0.11;
    auto const shift = 0.0;
    auto const offset = 0.1;
    auto const min = 0.0;
    auto const r_off = dist - offset;
    auto const cut = r_off + 1e-3; // LJ for only 2 pairs: AB BB
    make_particle_type_exist(std::max(type_a, type_b));
    auto const key_ab = get_ia_param_key(type_a, type_b);
    auto const key_bb = get_ia_param_key(type_b, type_b);
    LJ_Parameters lj{eps, sig, cut, offset, min, shift};
    ::nonbonded_ia_params[key_ab]->lj = lj;
    ::nonbonded_ia_params[key_bb]->lj = lj;
    on_non_bonded_ia_change();

    // matrix indices and reference energy value
    auto const max_type = type_b + 1;
    auto const n_pairs = Utils::upper_triangular(type_b, type_b, max_type) + 1;
    auto const lj_pair_ab = Utils::upper_triangular(type_a, type_b, max_type);
    auto const lj_pair_bb = Utils::upper_triangular(type_b, type_b, max_type);
    auto const frac6 = Utils::int_pow<6>(sig / r_off);
    auto const lj_energy = 4.0 * eps * (Utils::sqr(frac6) - frac6 + shift);

    // measure energies
    auto const obs_energy = calculate_energy();
    if (rank == 0) {
      for (int i = 0; i < n_pairs; ++i) {
        auto const ref_inter = (i == lj_pair_ab) ? lj_energy : 0.;
        auto const ref_intra = (i == lj_pair_bb) ? lj_energy : 0.;
        BOOST_CHECK_CLOSE(obs_energy->non_bonded_inter[i], ref_inter, 1e-10);
        BOOST_CHECK_CLOSE(obs_energy->non_bonded_intra[i], ref_intra, 1e-10);
      }
    }
  }
#endif // LENNARD_JONES

  // check bonded energies
  {
    // distance between particles
    auto const dist = 0.2;
    // set up a harmonic bond and a FENE bond, with a gap
    auto const harm_bond_id = 0;
    auto const none_bond_id = 1;
    auto const fene_bond_id = 2;
    {
      // set up a harmonic bond
      auto const bond = HarmonicBond(200.0, 0.3, 1.0);
      auto const bond_ia = std::make_shared<Bonded_IA_Parameters>(bond);
      bonded_ia_params.insert(harm_bond_id, bond_ia);
    }
    {
      // set up a FENE bond
      auto const bond = FeneBond(300.0, 1.0, 0.3);
      auto const bond_ia = std::make_shared<Bonded_IA_Parameters>(bond);
      bonded_ia_params.insert(fene_bond_id, bond_ia);
    }
    auto const &harm_bond =
        *boost::get<HarmonicBond>(bonded_ia_params.at(harm_bond_id).get());
    auto const &fene_bond =
        *boost::get<FeneBond>(bonded_ia_params.at(fene_bond_id).get());
    insert_particle_bond(pid2, harm_bond_id, {pid1});
    insert_particle_bond(pid2, fene_bond_id, {pid3});

    // measure energies
    auto const obs_energy = calculate_energy();
    if (rank == 0) {
      auto const none_energy = 0.0;
      auto const harm_energy =
          0.5 * harm_bond.k * Utils::sqr(harm_bond.r - dist);
      auto const fene_energy =
          -0.5 * fene_bond.k * Utils::sqr(fene_bond.drmax) *
          std::log(1.0 - Utils::sqr((dist - fene_bond.r0) / fene_bond.drmax));
      BOOST_CHECK_CLOSE(obs_energy->bonded[none_bond_id], none_energy, 0.0);
      BOOST_CHECK_CLOSE(obs_energy->bonded[harm_bond_id], harm_energy, 1e-10);
      BOOST_CHECK_CLOSE(obs_energy->bonded[fene_bond_id], fene_energy, 1e-10);
    }
  }

  // check electrostatics
#ifdef P3M
  {
    // add charges
    set_particle_property(pid1, &Particle::q, +1.);
    set_particle_property(pid2, &Particle::q, -1.);

    // set up P3M
    auto const prefactor = 2.;
    auto p3m = P3MParameters{false,
                             0.0,
                             3.5,
                             Utils::Vector3i::broadcast(12),
                             Utils::Vector3d::broadcast(0.5),
                             5,
                             0.615,
                             1e-3};
    auto solver =
        std::make_shared<CoulombP3M>(std::move(p3m), prefactor, 1, false, true);
    ::Coulomb::add_actor(solver);

    // measure energies
    auto const step = 0.02;
    auto const pos1 = start_positions.at(pid1);
    Utils::Vector3d pos2{box_center, box_center - 0.1, 1.0};
    for (int i = 0; i < 10; ++i) {
      // move particle
      pos2[0] += step;
      set_particle_pos(pid2, pos2);
      auto const r = (pos2 - pos1).norm();
      // check P3M energy
      auto const obs_energy = calculate_energy();
      if (rank == 0) {
        // at very short distances, the real-space contribution to
        // the energy is much larger than the k-space contribution
        auto const energy_ref = -prefactor / r;
        auto const energy_p3m = obs_energy->coulomb[0] + obs_energy->coulomb[1];
        BOOST_CHECK_CLOSE(energy_p3m, energy_ref, 0.002);
      }
    }
  }
#endif // P3M

  // check integration
  {
    // set up velocity-Verlet integrator
    set_integ_switch(INTEG_METHOD_NVT);

    // reset system
    remove_translational_motion();
    reset_particle_positions();

    // recalculate forces without propagating the system
    integrate(0, INTEG_REUSE_FORCES_CONDITIONALLY);

    // particles are arranged in a right triangle
    auto const p1_opt = copy_particle_to_head_node(comm, pid1);
    auto const p2_opt = copy_particle_to_head_node(comm, pid2);
    auto const p3_opt = copy_particle_to_head_node(comm, pid3);
    if (rank == 0) {
      auto const &p1 = *p1_opt;
      auto const &p2 = *p2_opt;
      auto const &p3 = *p3_opt;
      // forces are symmetric
      BOOST_CHECK_CLOSE(p1.force()[0], -p2.force()[0], tol);
      BOOST_CHECK_CLOSE(p3.force()[1], -p2.force()[1], tol);
      // periodic image contributions to the electrostatic force are negligible
      BOOST_CHECK_LE(std::abs(p1.force()[1]), tol);
      BOOST_CHECK_LE(std::abs(p1.force()[2]), tol);
      BOOST_CHECK_LE(std::abs(p2.force()[2]), tol);
      // zero long-range contribution for uncharged particles
      BOOST_CHECK_EQUAL(p3.force()[0], 0.);
      BOOST_CHECK_EQUAL(p3.force()[2], 0.);
      // velocities are not propagated
      BOOST_CHECK_EQUAL(p1.v().norm(), 0.);
      BOOST_CHECK_EQUAL(p2.v().norm(), 0.);
      BOOST_CHECK_EQUAL(p3.v().norm(), 0.);
    }

    // check integrated trajectory; the time step is chosen
    // small enough so that particles don't travel too far
#ifndef NDEBUG
    auto const pos_com = Utils::Vector3d{box_center, box_center, 1.0};
#endif
    auto const pids = std::vector<int>{pid1, pid2, pid3};
    for (int i = 0; i < 10; ++i) {
      std::unordered_map<int, Utils::Vector3d> expected;
      for (auto pid : pids) {
        auto p_opt = copy_particle_to_head_node(comm, pid);
        if (rank == 0) {
          auto &p = *p_opt;
          p.v() += 0.5 * time_step * p.force() / p.mass();
          p.pos() += time_step * p.v();
          expected[pid] = p.pos();
        }
      }
      integrate(1, INTEG_REUSE_FORCES_CONDITIONALLY);
      for (auto pid : pids) {
        auto const p_opt = copy_particle_to_head_node(comm, pid);
        if (rank == 0) {
          auto &p = *p_opt;
          BOOST_CHECK_LE((p.pos() - expected[pid]).norm(), tol);
          assert((p.pos() - pos_com).norm() < 0.5);
        }
      }
    }
  }

  // check exceptions
  {
    BOOST_CHECK_THROW(get_particle_node_parallel(-1), std::domain_error);
    if (rank == 0) {
      BOOST_CHECK_THROW(get_particle_node_parallel(12345), std::runtime_error);
    } else {
      get_particle_node_parallel(12345);
    }
  }
}

int main(int argc, char **argv) {
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
