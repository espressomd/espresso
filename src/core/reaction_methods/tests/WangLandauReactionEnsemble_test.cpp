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

/* Unit tests for the Wang-Landau Reaction Ensemble method. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE Wang Landau Reaction Ensemble test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "unit_tests/ParticleFactory.hpp"

#include "reaction_methods/WangLandauReactionEnsemble.hpp"
#include "reaction_methods/utils.hpp"

#include "communication.hpp"
#include "particle_data.hpp"

#include <boost/mpi.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/adaptor/strided.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>

#include <algorithm>
#include <cmath>
#include <istream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>

// Check the collective variable that tracks the degree of dissociation
// of a polybasic acid. The state measures the concentration of the base
// (i.e. fully deprotonated species) divided by the concentration of the
// base and all of its the conjugated acids.
BOOST_FIXTURE_TEST_CASE(DegreeOfAssociationCollectiveVariable_test,
                        ParticleFactory) {
  using namespace ReactionMethods;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // particle types
  int const type_A = 0;
  int const type_AH = 1;
  int const type_AH2 = 2;
  init_type_map(type_A);
  init_type_map(type_AH);
  init_type_map(type_AH2);

  // collective variable
  DegreeOfAssociationCollectiveVariable doa_cv{};
  doa_cv.corresponding_acid_types = {type_A, type_AH, type_AH2};
  doa_cv.associated_type = type_A;

  // exception if no base is present
  BOOST_CHECK_THROW(doa_cv.determine_current_state(), std::runtime_error);

  // add base
  create_particle({}, 0, type_A);
  BOOST_CHECK_CLOSE(doa_cv.determine_current_state(), 1.0, tol);

  // add acid
  create_particle({}, 1, type_AH);
  BOOST_CHECK_CLOSE(doa_cv.determine_current_state(), 0.5, tol);

  // add acid
  create_particle({}, 2, type_AH);
  create_particle({}, 3, type_AH2);
  BOOST_CHECK_CLOSE(doa_cv.determine_current_state(), 0.25, tol);
}

// Check the Monte Carlo algorithm where moves depend on the system
// configuration and/or energy.
BOOST_AUTO_TEST_CASE(WangLandauReactionEnsemble_test) {
  using namespace ReactionMethods;
  class WangLandauReactionEnsembleTest : public WangLandauReactionEnsemble {
  public:
    using WangLandauReactionEnsemble::calculate_acceptance_probability;
    using WangLandauReactionEnsemble::format_wang_landau_results;
    using WangLandauReactionEnsemble::WangLandauReactionEnsemble;
  };
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  WangLandauReactionEnsembleTest r_algo(42);
  r_algo.volume = 10.;
  r_algo.temperature = 20.;

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
        auto const f_expr = calculate_factorial_expression(reaction, p_numbers);
        auto const prefactor = 2e6 * f_expr;
        auto const boltzmann = std::exp(energy / 20.);
        double acceptance;
        r_algo.do_energy_reweighting = false;
        acceptance = r_algo.calculate_acceptance_probability(
            reaction, energy, 0., p_numbers, -1, -1, false);
        BOOST_CHECK_CLOSE(acceptance, prefactor * boltzmann, 5 * tol);
        acceptance = r_algo.calculate_acceptance_probability(
            reaction, energy, 0., p_numbers, -1, -1, true);
        BOOST_CHECK_CLOSE(acceptance, boltzmann, 5 * tol);
        r_algo.do_energy_reweighting = true;
        acceptance = r_algo.calculate_acceptance_probability(
            reaction, energy, 0., p_numbers, -1, -1, false);
        BOOST_CHECK_CLOSE(acceptance, prefactor, 5 * tol);
        acceptance = r_algo.calculate_acceptance_probability(
            reaction, energy, 0., p_numbers, -1, -1, true);
        BOOST_CHECK_CLOSE(acceptance, 1., 5 * tol);
        acceptance = r_algo.calculate_acceptance_probability(
            reaction, energy, 0., p_numbers, -1, 0, true);
        BOOST_CHECK_CLOSE(acceptance, 10., 5 * tol);
        acceptance = r_algo.calculate_acceptance_probability(
            reaction, energy, 0., p_numbers, 0, -1, true);
        BOOST_CHECK_CLOSE(acceptance, -10., 5 * tol);
      }
    }
  }

  // check collective variables
  {
    using namespace boost::adaptors;
    // setup input
    auto const delta_CV = 0.5;
    std::basic_stringstream<char> wl_file;
    wl_file << "# header 1.0  0.5\n";
    wl_file << "0  1.0  4.0\n";
    wl_file << "1\t2.0\t5.0\n";
    wl_file.flush();
    r_algo.do_energy_reweighting = false;
    // add collective variable
    r_algo.add_new_CV_potential_energy(wl_file, delta_CV);
    BOOST_REQUIRE_EQUAL(r_algo.collective_variables.size(), 1ul);
    BOOST_REQUIRE_EQUAL(r_algo.min_boundaries_energies.size(), 2ul);
    BOOST_REQUIRE_EQUAL(r_algo.max_boundaries_energies.size(), 2ul);
    BOOST_CHECK_EQUAL(r_algo.do_energy_reweighting, true);
    BOOST_CHECK_EQUAL(r_algo.collective_variables[0]->delta_CV, delta_CV);
    BOOST_CHECK_EQUAL(r_algo.collective_variables[0]->CV_minimum, 1.);
    BOOST_CHECK_EQUAL(r_algo.collective_variables[0]->CV_maximum, 5.);
    BOOST_CHECK_CLOSE(r_algo.min_boundaries_energies[0], 1., tol);
    BOOST_CHECK_CLOSE(r_algo.min_boundaries_energies[1], 2., tol);
    BOOST_CHECK_CLOSE(r_algo.max_boundaries_energies[0], 4., tol);
    BOOST_CHECK_CLOSE(r_algo.max_boundaries_energies[1], 5., tol);
    // check Wang-Landau histogram
    {
      std::basic_stringstream<char> outfile;
      r_algo.format_wang_landau_results(outfile);
      std::istream_iterator<double> start(outfile), end;
      std::vector<double> flat_array(start, end);
      BOOST_REQUIRE_EQUAL(flat_array.size(), 14ul);
      std::vector<double> bin_edges;
      std::vector<double> histogram;
      boost::range::push_back(bin_edges, flat_array | strided(2));
      boost::range::push_back(
          histogram, flat_array | sliced(1, flat_array.size()) | strided(2));
      std::vector<double> bin_edges_ref(7u, 0.);
      std::vector<double> histogram_ref(7u, 0.);
      std::generate(bin_edges_ref.begin(), bin_edges_ref.end(),
                    [n = .5]() mutable { return n += .5; });
      BOOST_TEST(bin_edges == bin_edges_ref, boost::test_tools::per_element());
      BOOST_TEST(histogram == histogram_ref, boost::test_tools::per_element());
    }
    // exception if file doesn't exist
    BOOST_CHECK_THROW(r_algo.add_new_CV_potential_energy("unknown", 0.),
                      std::runtime_error);
  }
}

int main(int argc, char **argv) {
  auto mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  Communication::init(mpi_env);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
