/*
 * Copyright (C) 2021 The ESPResSo project
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

#include "communication.hpp"

#include <boost/mpi.hpp>

#include <limits>
#include <memory>
#include <stdexcept>

// Check the base class for all Monte Carlo algorithms.
BOOST_AUTO_TEST_CASE(ReactionAlgorithm_test) {
  using namespace ReactionMethods;
  class ReactionAlgorithmTest : public ReactionAlgorithm {
  public:
    using ReactionAlgorithm::calculate_acceptance_probability;
    using ReactionAlgorithm::ReactionAlgorithm;
  };
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // check acceptance rate
  ReactionAlgorithmTest r_algo(42);
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

  // exception if no reaction was added
  BOOST_CHECK_THROW(r_algo.check_reaction_method(), std::runtime_error);

  // create a reaction A -> 3 B + 4 C
  int const type_A = 0;
  int const type_B = 1;
  int const type_C = 2;
  SingleReaction const reaction(2., {type_A}, {1}, {type_B, type_C}, {3, 4});

  // check reaction addition
  {
    r_algo.add_reaction(reaction.gamma, reaction.reactant_types,
                        reaction.reactant_coefficients, reaction.product_types,
                        reaction.product_coefficients);
    BOOST_REQUIRE_EQUAL(r_algo.reactions.size(), 1ul);
    auto const &value = r_algo.reactions[0];
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
    double probability = r_algo.calculate_acceptance_probability(
        reaction, -1., -1., {{1, 2}}, -1, -1, false);
    BOOST_CHECK_EQUAL(probability, -10.);
  }

  // exception if temperature is negative
  BOOST_CHECK_THROW(r_algo.check_reaction_method(), std::runtime_error);

  // set temperature
  r_algo.temperature = 10;

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
    SingleReaction const new_reaction(5., {type_B}, {1}, {type_C}, {1});
    r_algo.add_reaction(new_reaction.gamma, new_reaction.reactant_types,
                        new_reaction.reactant_coefficients,
                        new_reaction.product_types,
                        new_reaction.product_coefficients);
    BOOST_REQUIRE_EQUAL(r_algo.reactions.size(), 2ul);
    BOOST_CHECK_EQUAL(r_algo.reactions[1].gamma, new_reaction.gamma);
    r_algo.delete_reaction(1);
    BOOST_REQUIRE_EQUAL(r_algo.reactions.size(), 1ul);
    BOOST_CHECK_EQUAL(r_algo.reactions[0].gamma, reaction.gamma);
    r_algo.delete_reaction(0);
    BOOST_REQUIRE_EQUAL(r_algo.reactions.size(), 0ul);
  }

  // exception if deleting a non-existent particle
  BOOST_CHECK_THROW(r_algo.delete_particle(5), std::runtime_error);
}

int main(int argc, char **argv) {
  auto mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  Communication::init(mpi_env);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
