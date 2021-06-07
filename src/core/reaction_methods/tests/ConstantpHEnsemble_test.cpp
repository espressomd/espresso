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

/* Unit tests for the Constant pH Ensemble. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE Constant pH Ensemble test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "reaction_methods/ConstantpHEnsemble.hpp"
#include "reaction_methods/utils.hpp"

#include "communication.hpp"
#include "particle_data.hpp"

#include <boost/mpi.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>

// Check the Monte Carlo algorithm where moves depend on the system
// configuration, energy and pH.
BOOST_AUTO_TEST_CASE(ConstantpHEnsemble_test) {
  using namespace ReactionMethods;
  class ConstantpHEnsembleTest : public ConstantpHEnsemble {
  public:
    using ConstantpHEnsemble::calculate_acceptance_probability;
    using ConstantpHEnsemble::ConstantpHEnsemble;
  };
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  ConstantpHEnsembleTest r_algo(42);
  r_algo.temperature = 20.;
  r_algo.m_constant_pH = 1.;
  r_algo.m_conversion_factor_from_mol_per_l_to_1_div_sigma_cubed = 1.0;

  // exception if no reaction was added
  BOOST_CHECK_THROW(r_algo.check_reaction_method(), std::runtime_error);

  // create a reaction A -> B + C
  int const type_A = 0;
  int const type_B = 1;
  int const type_C = 2;
  SingleReaction const reaction(2., {type_A}, {1}, {type_B, type_C}, {1, 1}, 1);

  // check acceptance probability
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        // system contains i x A, j x B, and k x C
        auto const p_numbers =
            std::map<int, int>{{type_A, i}, {type_B, j}, {type_C, k}};
        auto const energy = static_cast<double>(i + 1);
        auto const f_expr =
            calculate_factorial_expression_cpH(reaction, p_numbers);
        // acceptance = f_expr * exp(- E / T + nu_bar * log(10) * (pH - nu_bar *
        // pKa))
        auto const acceptance_ref =
            f_expr * std::exp(energy / r_algo.temperature +
                              std::log(10.) * (r_algo.m_constant_pH +
                                               std::log10(reaction.gamma)));
        auto const acceptance = r_algo.calculate_acceptance_probability(
            reaction, energy, 0., p_numbers, -1, -1, false);
        BOOST_CHECK_CLOSE(acceptance, acceptance_ref, 5 * tol);
      }
    }
  }
}

int main(int argc, char **argv) {
  auto mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  Communication::init(mpi_env);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
