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

/* Unit tests for the Reaction Ensemble. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE Reaction Ensemble test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "reaction_methods/ReactionEnsemble.hpp"
#include "reaction_methods/utils.hpp"

#include "EspressoSystemStandAlone.hpp"
#include "communication.hpp"
#include "particle_data.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>

namespace espresso {
// ESPResSo system instance
std::unique_ptr<EspressoSystemStandAlone> system;
} // namespace espresso

// Check the Monte Carlo algorithm where moves depend on the system
// configuration and energy.
BOOST_AUTO_TEST_CASE(ReactionEnsemble_test) {
  using namespace ReactionMethods;
  class ReactionEnsembleTest : public ReactionEnsemble {
  public:
    using ReactionEnsemble::calculate_acceptance_probability;
    using ReactionEnsemble::generic_oneway_reaction;
    using ReactionEnsemble::ReactionEnsemble;
  };
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // check basic interface
  {
    ReactionEnsembleTest r_algo(42);
    r_algo.volume = 10.;
    r_algo.kT = 20.;

    // exception if no reaction was added
    BOOST_CHECK_THROW(r_algo.check_reaction_method(), std::runtime_error);

    // create a reaction A -> 3 B + 4 C
    int const type_A = 0;
    int const type_B = 1;
    int const type_C = 2;
    SingleReaction const reaction(2., {type_A}, {1}, {type_B, type_C}, {3, 4});

    // check acceptance probability
    constexpr auto g = factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
          // system contains i x A, j x B, and k x C
          auto const p_numbers =
              std::map<int, int>{{type_A, i}, {type_B, j}, {type_C, k}};
          auto const energy = static_cast<double>(i + 1);
          auto const f_expr =
              calculate_factorial_expression(reaction, p_numbers);
          // acceptance = V^{nu_bar} * gamma * f_expr * exp(- E / T)
          auto const acceptance_ref = std::pow(r_algo.volume, reaction.nu_bar) *
                                      reaction.gamma * f_expr *
                                      std::exp(energy / r_algo.kT);
          auto const acceptance = r_algo.calculate_acceptance_probability(
              reaction, energy, 0., p_numbers);
          BOOST_CHECK_CLOSE(acceptance, acceptance_ref, 5 * tol);
        }
      }
    }
  }

  // check that the system energy is updated after a succesful reaction
  {
    ReactionEnsembleTest test_reaction(42);
    test_reaction.volume = 1.;
    test_reaction.kT = 1.;
    test_reaction.exclusion_radius = 0;

    // create a generic identity exchange reaction D <-> E
    int const type_D = 0;
    int const type_E = 1;

    test_reaction.charges_of_types[type_D] = 0;
    test_reaction.charges_of_types[type_E] = 0;

    // track particles
    init_type_map(type_D);
    init_type_map(type_E);

    auto const gamma = 1e100; // reaction completly shifted to product creation
    ReactionMethods::SingleReaction reaction(gamma, {type_D}, {1}, {type_E},
                                             {1});

    // resize system box
    auto const box_l = Utils::Vector3d{0.5, 0.4, 0.7};
    espresso::system->set_box_l(box_l);

    // create a D particle in the system
    auto const pid = 0;
    auto const ref_position = Utils::Vector3d{0.1, 0.2, 0.3};
    place_particle(pid, ref_position);
    set_particle_type(pid, type_D);

    // sentinel value to check energies get updated
    double energy = 1.0;

    // for an ideal system with gamma ~ inf, the reaction is always accepted
    test_reaction.generic_oneway_reaction(reaction, energy);

    // the potential energy of the new state has to be 0 since it is an ideal
    // system
    double const energy_ref = 0.0;
    BOOST_CHECK_CLOSE(energy, energy_ref, tol);

    // the reaction was updated
    BOOST_CHECK_EQUAL(reaction.tried_moves, 1);
    BOOST_CHECK_EQUAL(reaction.accepted_moves, 1);
  }
}

int main(int argc, char **argv) {
  auto mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);
  Communication::init(mpi_env);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
