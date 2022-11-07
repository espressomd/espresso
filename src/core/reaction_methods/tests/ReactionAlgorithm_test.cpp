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

/* Unit tests for the Reaction Algorithm. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ReactionAlgorithm test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "reaction_methods/ReactionAlgorithm.hpp"

#include "EspressoSystemStandAlone.hpp"
#include "Particle.hpp"
#include "communication.hpp"
#include "particle_data.hpp"
#include "particle_node.hpp"

#include <boost/mpi.hpp>

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace espresso {
// ESPResSo system instance
static std::unique_ptr<EspressoSystemStandAlone> system;
} // namespace espresso

// Check the base class for all Monte Carlo algorithms.
BOOST_AUTO_TEST_CASE(ReactionAlgorithm_test) {
  using namespace ReactionMethods;
  class ReactionAlgorithmTest : public ReactionAlgorithm {
  public:
    using ReactionAlgorithm::calculate_acceptance_probability;
    using ReactionAlgorithm::get_random_position_in_box;
    using ReactionAlgorithm::ReactionAlgorithm;
  };
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // check acceptance rate
  ReactionAlgorithmTest r_algo(42, 1., 0., {});
  for (int tried_moves = 1; tried_moves < 5; ++tried_moves) {
    for (int accepted_moves = 0; accepted_moves < 5; ++accepted_moves) {
      r_algo.N_trial_particle_displacement_MC_moves = tried_moves;
      r_algo.N_accepted_particle_displacement_MC_moves = accepted_moves;
      auto const ref_rate = static_cast<double>(accepted_moves) /
                            static_cast<double>(tried_moves);
      BOOST_CHECK_CLOSE(
          r_algo.get_acceptance_rate_particle_displacement_MC_moves(), ref_rate,
          tol);
    }
  }

  // exception if no reaction was added
  BOOST_CHECK_THROW(r_algo.check_reaction_method(), std::runtime_error);

  // create a reaction A -> 3 B + 4 C
  int const type_A = 0;
  int const type_B = 1;
  int const type_C = 2;
  SingleReaction const reaction(2., {type_A}, {1}, {type_B, type_C}, {3, 4});

  // track particles
  init_type_map(type_A);
  init_type_map(type_B);
  init_type_map(type_C);

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

  // check acceptance probability
  {
    double probability =
        r_algo.calculate_acceptance_probability(reaction, -1., -1., {{1, 2}});
    BOOST_CHECK_EQUAL(probability, -10.);
  }

#ifdef ELECTROSTATICS
  // exception if reactant types have no charge information
  BOOST_CHECK_THROW(r_algo.check_reaction_method(), std::runtime_error);
  r_algo.charges_of_types[0] = 1;
  // exception if product types have no charge information
  BOOST_CHECK_THROW(r_algo.check_reaction_method(), std::runtime_error);
  r_algo.charges_of_types[1] = 1;
  r_algo.charges_of_types[2] = 0;
#endif

  // sanity checks should now pass
  r_algo.check_reaction_method();

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

  // check particle displacement Monte Carlo moves
  {

    // set up one particle

    auto const box_l = 1.;
    std::vector<Utils::Vector3d> ref_positions{{0.1, 0.2, 0.3},
                                               {0.4, 0.5, 0.6}};
    place_particle(0, ref_positions[0]);
    set_particle_type(0, type_A);

    // check that the particle moves if there is no restriction
    r_algo.do_particle_displacement_MC_move(1, {});
    BOOST_CHECK_GE((get_particle_data(0).pos() - ref_positions[0]).norm(), tol);

    // check that only particles with types in particle_types_to_move are
    // allowed to move

    place_particle(0, ref_positions[0]);
    place_particle(1, ref_positions[1]);

    set_particle_type(1, type_B);
    r_algo.do_particle_displacement_MC_move(10, {type_B});

    BOOST_CHECK_EQUAL((get_particle_data(0).pos() - ref_positions[0]).norm(),
                      0);
    BOOST_CHECK_GE((get_particle_data(1).pos() - ref_positions[1]).norm(), tol);

    // check that no particle moves if they do not match  types in
    // particle_types_to_move

    place_particle(0, ref_positions[0]);
    place_particle(1, ref_positions[1]);
    r_algo.do_particle_displacement_MC_move(10, {3});
    BOOST_CHECK_EQUAL((get_particle_data(0).pos() - ref_positions[0]).norm(),
                      0);
    BOOST_CHECK_EQUAL((get_particle_data(1).pos() - ref_positions[1]).norm(),
                      0);
    // force all MC moves to be rejected by picking particles inside
    // their exclusion radius

    r_algo.exclusion_range = box_l;
    r_algo.particle_inside_exclusion_range_touched = false;

    place_particle(0, ref_positions[0]);
    place_particle(1, ref_positions[1]);

    r_algo.do_particle_displacement_MC_move(10, {type_A, type_B});

    BOOST_CHECK_EQUAL((get_particle_data(0).pos() - ref_positions[0]).norm(),
                      0);
    BOOST_CHECK_EQUAL((get_particle_data(1).pos() - ref_positions[1]).norm(),
                      0);
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

    BOOST_CHECK_THROW(ReactionAlgorithmTest r_algo(40, 1., -1, {}),
                      std::domain_error);

    // domain error if a negative value is provided in exclusion_radius_per_type
    std::unordered_map<int, double> exclusion_radius_per_type;
    exclusion_radius_per_type[type_A] = 1;
    exclusion_radius_per_type[type_B] = -1;

    BOOST_CHECK_THROW(
        ReactionAlgorithmTest r_algo(40, 1., 1, exclusion_radius_per_type),
        std::domain_error);

    auto const box_l = Utils::Vector3d{1, 1, 1};
    espresso::system->set_box_l(box_l);

    // set up particle
    place_particle(0, {0.5, 0.5, 0.5});
    set_particle_type(0, type_A);
    place_particle(1, {0.7, 0.7, 0.7});
    set_particle_type(1, type_B);
    exclusion_radius_per_type[type_A] = 0.1;
    exclusion_radius_per_type[type_B] = 1;
    ReactionAlgorithmTest r_algo(40, 1., 0, exclusion_radius_per_type);
  }
}

int main(int argc, char **argv) {
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
