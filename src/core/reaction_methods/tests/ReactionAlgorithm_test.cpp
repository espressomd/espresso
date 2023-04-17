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
#define BOOST_TEST_MODULE ReactionAlgorithm test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config/config.hpp"

#include "reaction_methods/ReactionAlgorithm.hpp"
#include "reaction_methods/SingleReaction.hpp"

#include "EspressoSystemStandAlone.hpp"
#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "particle_node.hpp"
#include "unit_tests/ParticleFactory.hpp"

#include <boost/mpi.hpp>

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace espresso {
// ESPResSo system instance
static std::unique_ptr<EspressoSystemStandAlone> system;
} // namespace espresso

namespace Testing {
class ReactionAlgorithm : public ReactionMethods::ReactionAlgorithm {
public:
  using Base = ReactionMethods::ReactionAlgorithm;
  using Base::clear_old_system_state;
  using Base::displacement_mc_move;
  using Base::get_old_system_state;
  using Base::get_random_position_in_box;
  using Base::make_displacement_mc_move_attempt;
  using Base::ReactionAlgorithm;
};
} // namespace Testing

// Check the base class for all Monte Carlo algorithms.
BOOST_FIXTURE_TEST_CASE(ReactionAlgorithm_test, ParticleFactory) {
  using ReactionMethods::SingleReaction;
  auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();
  auto const comm = boost::mpi::communicator();

  // check acceptance rate
  auto r_algo = Testing::ReactionAlgorithm(comm, 42, 1., 0., {});
  for (int tried_moves = 1; tried_moves < 5; ++tried_moves) {
    for (int accepted_moves = 0; accepted_moves < 5; ++accepted_moves) {
      r_algo.m_tried_configurational_MC_moves = tried_moves;
      r_algo.m_accepted_configurational_MC_moves = accepted_moves;
      auto const ref_rate = static_cast<double>(accepted_moves) /
                            static_cast<double>(tried_moves);
      BOOST_CHECK_CLOSE(r_algo.get_acceptance_rate_configurational_moves(),
                        ref_rate, tol);
    }
  }

  // create a reaction A -> 3 B + 4 C
  int const type_A = 0;
  int const type_B = 1;
  int const type_C = 2;
  SingleReaction const reaction(2., {type_A}, {1}, {type_B, type_C}, {3, 4});
  r_algo.non_interacting_type = 5;

  // check reaction addition
  {
    r_algo.add_reaction(std::make_shared<SingleReaction>(
        reaction.gamma, reaction.reactant_types, reaction.reactant_coefficients,
        reaction.product_types, reaction.product_coefficients));
    BOOST_REQUIRE_EQUAL(r_algo.reactions.size(), 1ul);
    auto const &value = *r_algo.reactions[0];
    BOOST_TEST(value.reactant_types == reaction.reactant_types,
               boost::test_tools::per_element());
    BOOST_TEST(value.reactant_coefficients == reaction.reactant_coefficients,
               boost::test_tools::per_element());
    BOOST_TEST(value.product_types == reaction.product_types,
               boost::test_tools::per_element());
    BOOST_TEST(value.product_coefficients == reaction.product_coefficients,
               boost::test_tools::per_element());
    BOOST_CHECK_EQUAL(value.gamma, reaction.gamma);
  }

  // check reaction removal
  {
    auto const new_gamma = 5.;
    auto const new_reaction = std::make_shared<SingleReaction>(
        new_gamma, std::vector<int>{type_B}, std::vector<int>{1},
        std::vector<int>{type_C}, std::vector<int>{1});
    r_algo.add_reaction(new_reaction);
    BOOST_REQUIRE_EQUAL(r_algo.reactions.size(), 2ul);
    BOOST_CHECK_EQUAL(r_algo.reactions[1]->gamma, new_gamma);
    r_algo.delete_reaction(1);
    BOOST_REQUIRE_EQUAL(r_algo.reactions.size(), 1ul);
    BOOST_CHECK_EQUAL(r_algo.reactions[0]->gamma, reaction.gamma);
    r_algo.delete_reaction(0);
    BOOST_REQUIRE_EQUAL(r_algo.reactions.size(), 0ul);
  }

  // exception if deleting a non-existent particle
  BOOST_CHECK_THROW(r_algo.delete_particle(5), std::runtime_error);

  // check particle moves
  {
    // set up particles
    auto const box_l = 1.;
    std::vector<std::pair<Utils::Vector3d, Utils::Vector3d>> ref_positions{
        {{0.1, 0.2, 0.3}, {+10., +20., +30.}},
        {{0.4, 0.5, 0.6}, {-10., -20., -30.}}};
    ::make_new_particle(0, ref_positions[0].first);
    ::make_new_particle(1, ref_positions[1].first);
    set_particle_v(0, ref_positions[0].second);
    set_particle_v(1, ref_positions[1].second);
    set_particle_type(0, type_A);
    set_particle_type(1, type_A);
    // update particle positions and velocities
    BOOST_CHECK(!r_algo.particle_inside_exclusion_range_touched);
    r_algo.particle_inside_exclusion_range_touched = false;
    r_algo.exclusion_range = box_l;
    r_algo.displacement_mc_move(0, 2);
    auto const &bookkeeping = r_algo.get_old_system_state();
    BOOST_CHECK(r_algo.particle_inside_exclusion_range_touched);
    // check moves and bookkeeping
    for (auto const &[pid, old_pos, old_vel] : bookkeeping.moved) {
      BOOST_REQUIRE(pid == 0 or pid == 1);
      auto const ref_old_pos = ref_positions[pid].first;
      auto const ref_old_vel = ref_positions[pid].second;
      if (auto const p = ::cell_structure.get_local_particle(pid)) {
        auto const &new_pos = p->pos();
        auto const &new_vel = p->v();
        BOOST_CHECK_EQUAL(old_pos, ref_old_pos);
        BOOST_CHECK_EQUAL(old_vel, ref_old_vel);
        BOOST_CHECK_GE(new_pos, Utils::Vector3d::broadcast(0.));
        BOOST_CHECK_LE(new_pos, Utils::Vector3d::broadcast(box_l));
        BOOST_CHECK_GE((new_pos - ref_old_pos).norm(), 0.1);
        BOOST_CHECK_GE((new_vel - ref_old_vel).norm(), 10.);
      }
    }
    // cleanup
    r_algo.clear_old_system_state();
    remove_particle(0);
    remove_particle(1);
  }

  // check Monte Carlo moves
  {
    // set up particles
    auto const box_l = 1.;
    std::vector<Utils::Vector3d> ref_positions{{0.1, 0.2, 0.3},
                                               {0.6, 0.7, 0.8}};
    ::make_new_particle(0, ref_positions[0]);
    ::make_new_particle(1, ref_positions[1]);
    set_particle_type(0, type_A);
    set_particle_type(1, type_A);
    // check runtime errors when a MC move cannot be physically performed
    auto displacement_move = std::bind(
        &Testing::ReactionAlgorithm::make_displacement_mc_move_attempt, &r_algo,
        std::placeholders::_1, std::placeholders::_2);
    BOOST_REQUIRE(!displacement_move(type_C, 1));
    BOOST_REQUIRE(!displacement_move(type_B, 2));
    BOOST_REQUIRE(!displacement_move(type_A, 0));
    BOOST_CHECK_THROW(displacement_move(type_A, -2), std::domain_error);
    BOOST_CHECK_THROW(displacement_move(-2, 1), std::domain_error);
    // force all MC moves to be rejected by picking particles inside
    // their exclusion radius
    r_algo.exclusion_range = box_l;
    BOOST_REQUIRE(!displacement_move(type_A, 2));
    // check none of the particles moved
    for (auto const pid : {0, 1}) {
      auto const ref_old_pos = ref_positions[pid];
      if (auto const p = ::cell_structure.get_local_particle(pid)) {
        auto const &new_pos = p->pos();
        BOOST_CHECK_LE((new_pos - ref_old_pos).norm(), tol);
      }
    }
    // force a MC move to be accepted by using a constant Hamiltonian
    r_algo.exclusion_range = 0.;
    BOOST_REQUIRE(displacement_move(type_A, 1));
    std::vector<double> distances(2);
    // check that only one particle moved
    for (auto const pid : {0, 1}) {
      if (auto const p = ::cell_structure.get_local_particle(pid)) {
        distances[pid] = (ref_positions[pid] - p->pos()).norm();
      }
    }
    BOOST_CHECK_LE(std::min(distances[0], distances[1]), tol);
    BOOST_CHECK_GE(std::max(distances[0], distances[1]), 0.1);
    // cleanup
    remove_particle(0);
    remove_particle(1);
  }

  // check random positions generator
  {
    // setup box
    auto const box_l = Utils::Vector3d{0.5, 0.4, 0.7};
    auto const origin = Utils::Vector3d::broadcast(0.);
    espresso::system->set_box_l(box_l);
    // cubic case
    r_algo.remove_constraint();
    for (int i = 0; i < 100; ++i) {
      auto const pos = r_algo.get_random_position_in_box();
      BOOST_TEST(pos <= box_l, boost::test_tools::per_element());
      BOOST_TEST(pos >= origin, boost::test_tools::per_element());
      // nodes must agree on random numbers
      auto pos_head_node = Utils::Vector3d::broadcast(0.);
      if (comm.rank() == 0) {
        pos_head_node = pos;
      }
      boost::mpi::broadcast(comm, pos_head_node, 0);
      BOOST_TEST(pos == pos_head_node, boost::test_tools::per_element());
    }
    // slab case
    auto const start_z{0.2}, end_z{0.6};
    auto const slab_lower = Utils::Vector3d{0., 0., start_z};
    auto const slab_upper = Utils::Vector3d{box_l[0], box_l[1], end_z};
    r_algo.set_slab_constraint(start_z, end_z);
    for (int i = 0; i < 100; ++i) {
      auto const pos = r_algo.get_random_position_in_box();
      BOOST_TEST(pos <= slab_upper, boost::test_tools::per_element());
      BOOST_TEST(pos >= slab_lower, boost::test_tools::per_element());
    }
    auto const slab_params = r_algo.get_slab_constraint_parameters();
    BOOST_CHECK_CLOSE(slab_params[0], start_z, tol);
    BOOST_CHECK_CLOSE(slab_params[1], end_z, tol);
    // cylindrical case
    auto const cyl_x{0.2}, cyl_y{0.1}, radius{0.2};
    r_algo.set_cyl_constraint(cyl_x, cyl_y, radius);
    for (int i = 0; i < 400; ++i) {
      auto const pos = r_algo.get_random_position_in_box();
      auto const z = pos[2];
      auto const r = Utils::Vector2d{pos[0] - cyl_x, pos[1] - cyl_y}.norm();
      BOOST_REQUIRE_LE(r, radius);
      BOOST_REQUIRE_LE(z, box_l[2]);
      BOOST_REQUIRE_GE(z, 0.);
    }
    // restore box geometry
    espresso::system->set_box_l(Utils::Vector3d::broadcast(1.));
    // check exception mechanism
    using exception = std::domain_error;
    BOOST_CHECK_THROW(r_algo.set_slab_constraint(-1., 0.5), exception);
    BOOST_CHECK_THROW(r_algo.set_slab_constraint(0.5, 1.5), exception);
    BOOST_CHECK_THROW(r_algo.set_slab_constraint(0.5, 0.2), exception);
    BOOST_CHECK_THROW(r_algo.set_cyl_constraint(-1., 0.5, 0.5), exception);
    BOOST_CHECK_THROW(r_algo.set_cyl_constraint(1.5, 0.5, 0.5), exception);
    BOOST_CHECK_THROW(r_algo.set_cyl_constraint(0.5, -1., 0.5), exception);
    BOOST_CHECK_THROW(r_algo.set_cyl_constraint(0.5, 1.5, 0.5), exception);
    BOOST_CHECK_THROW(r_algo.set_cyl_constraint(0.5, 0.5, -1.0), exception);
    r_algo.remove_constraint();
  }

  {
    // domain error if negative exclusion_range is provided
    BOOST_CHECK_THROW(Testing::ReactionAlgorithm(comm, 40, 1., -1, {}),
                      std::domain_error);

    // domain error if a negative value is provided in exclusion_radius_per_type
    std::unordered_map<int, double> exclusion_radius_per_type;
    exclusion_radius_per_type[type_A] = 1;
    exclusion_radius_per_type[type_B] = -1;

    BOOST_CHECK_THROW(
        Testing::ReactionAlgorithm(comm, 40, 1., 1, exclusion_radius_per_type),
        std::domain_error);

    espresso::system->set_box_l({1., 1., 1.});

    // set up particle
    ::make_new_particle(0, {0.5, 0.5, 0.5});
    ::make_new_particle(1, {0.7, 0.7, 0.7});
    set_particle_type(0, type_A);
    set_particle_type(1, type_B);
    exclusion_radius_per_type[type_A] = 0.1;
    exclusion_radius_per_type[type_B] = 1;
    auto r_algo =
        Testing::ReactionAlgorithm(comm, 40, 1., 0, exclusion_radius_per_type);

    // the new position will always be in the excluded range since the sum of
    // the radii of both particle types is larger than box length. The exclusion
    // range value should be ignored
    r_algo.displacement_mc_move(type_B, 1);
    BOOST_REQUIRE(r_algo.particle_inside_exclusion_range_touched);
    r_algo.clear_old_system_state();

    // the new position will never be in the excluded range because the
    // exclusion_radius of the particle is 0
    r_algo.exclusion_radius_per_type[type_B] = 0;
    r_algo.particle_inside_exclusion_range_touched = false;
    r_algo.displacement_mc_move(type_B, 1);
    BOOST_REQUIRE(!r_algo.particle_inside_exclusion_range_touched);
    r_algo.clear_old_system_state();

    // the new position will never accepted since the value in exclusion_range
    // will be used if the particle does not have a defined excluded radius
    r_algo.exclusion_range = 1;
    r_algo.exclusion_radius_per_type = {{type_A, 0}};
    r_algo.displacement_mc_move(type_B, 1);
    BOOST_REQUIRE(r_algo.particle_inside_exclusion_range_touched);
    r_algo.clear_old_system_state();
  }
}

int main(int argc, char **argv) {
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
