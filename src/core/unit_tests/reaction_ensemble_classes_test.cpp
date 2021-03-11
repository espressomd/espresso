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

/* Unit tests for the ReactionEnsemble classes. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ReactionEnsemble classes test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "reaction_ensemble.hpp"

#include "communication.hpp"
#include "particle_data.hpp"

#include <boost/mpi.hpp>

#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

/** Fixture to create particles during a test and remove them at the end. */
struct ParticleFactory {
  ParticleFactory() = default;
  ~ParticleFactory() {
    for (auto pid : particle_cache) {
      remove_particle(pid);
    }
  }
  void create_particle(int pid, int type) {
    double pos[3] = {0., 0., 0.};
    place_particle(pid, pos);
    set_particle_type(pid, type);
    particle_cache.emplace_back(pid);
  }

private:
  std::vector<int> particle_cache;
};

BOOST_FIXTURE_TEST_CASE(particle_type_map_test, ParticleFactory) {
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // particle properties
  int const type = 10;
  int const pid = 1;

  // exception for untracked particle ids
  BOOST_CHECK_THROW(number_of_particles_with_type(type), std::runtime_error);

  // exception for negative particle ids
  BOOST_CHECK_THROW(init_type_map(-10), std::runtime_error);

  // check particle counting
  init_type_map(type);
  BOOST_CHECK_EQUAL(number_of_particles_with_type(type), 0);
  create_particle(pid, type);
  BOOST_CHECK_EQUAL(number_of_particles_with_type(type), 1);

  // exception for random index that exceeds the number of particles
  BOOST_CHECK_THROW(get_random_p_id(type, 10), std::runtime_error);

  // check particle selection
  BOOST_CHECK_EQUAL(get_random_p_id(type, 0), pid);
}

BOOST_FIXTURE_TEST_CASE(DegreeOfAssociationCollectiveVariable_test,
                        ParticleFactory) {
  using namespace ReactionEnsemble;
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
  create_particle(0, type_A);
  BOOST_CHECK_CLOSE(doa_cv.determine_current_state(), 1.0, tol);

  // add acid
  create_particle(1, type_AH);
  BOOST_CHECK_CLOSE(doa_cv.determine_current_state(), 0.5, tol);

  // add acid
  create_particle(2, type_AH);
  create_particle(3, type_AH2);
  BOOST_CHECK_CLOSE(doa_cv.determine_current_state(), 0.25, tol);
}

BOOST_AUTO_TEST_CASE(SingleReaction_test) {
  using namespace ReactionEnsemble;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // check derived parameter
  SingleReaction reaction(2., {0}, {1}, {1, 2}, {3, 4});
  BOOST_CHECK_EQUAL(reaction.nu_bar, 6);

  // check acceptance rate
  for (int tried_moves = 1; tried_moves < 5; ++tried_moves) {
    for (int accepted_moves = 0; accepted_moves < 5; ++accepted_moves) {
      reaction.tried_moves = tried_moves;
      reaction.accepted_moves = accepted_moves;
      auto const ref_rate = static_cast<double>(accepted_moves) /
                            static_cast<double>(tried_moves);
      BOOST_CHECK_CLOSE(reaction.get_acceptance_rate(), ref_rate, tol);
    }
  }
}

BOOST_AUTO_TEST_CASE(ReactionAlgorithm_test) {
  using namespace ReactionEnsemble;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // check acceptance rate
  ReactionAlgorithm r_algo(42);
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
  BOOST_CHECK_THROW(r_algo.check_reaction_ensemble(), std::runtime_error);

  // check reaction addition
  {
    SingleReaction const ref(2., {0}, {1}, {1, 2}, {3, 4});
    r_algo.add_reaction(ref.gamma, ref.reactant_types,
                        ref.reactant_coefficients, ref.product_types,
                        ref.product_coefficients);
    BOOST_CHECK_EQUAL(r_algo.reactions.size(), 1ul);
    auto const &reaction = r_algo.reactions[0];
    BOOST_TEST(reaction.reactant_types == ref.reactant_types,
               boost::test_tools::per_element());
    BOOST_TEST(reaction.reactant_coefficients == ref.reactant_coefficients,
               boost::test_tools::per_element());
    BOOST_TEST(reaction.product_types == ref.product_types,
               boost::test_tools::per_element());
    BOOST_TEST(reaction.product_coefficients == ref.product_coefficients,
               boost::test_tools::per_element());
    BOOST_CHECK_EQUAL(reaction.gamma, ref.gamma);
  }

  // exception if temperature is negative
  BOOST_CHECK_THROW(r_algo.check_reaction_ensemble(), std::runtime_error);

  // set temperature
  r_algo.temperature = 10;

#ifdef ELECTROSTATICS
  // exception if reactant types have no charge information
  BOOST_CHECK_THROW(r_algo.check_reaction_ensemble(), std::runtime_error);
  r_algo.charges_of_types[0] = 1;
  // exception if product types have no charge information
  BOOST_CHECK_THROW(r_algo.check_reaction_ensemble(), std::runtime_error);
  r_algo.charges_of_types[1] = 1;
  r_algo.charges_of_types[2] = 0;
#endif

  // sanity checks should now pass
  r_algo.check_reaction_ensemble();

  // check reaction removal
  {
    r_algo.add_reaction(5, {1}, {1}, {2}, {1});
    BOOST_CHECK_EQUAL(r_algo.reactions.size(), 2ul);
    BOOST_CHECK_EQUAL(r_algo.reactions[1].gamma, 5);
    r_algo.delete_reaction(1);
    BOOST_CHECK_EQUAL(r_algo.reactions.size(), 1ul);
    BOOST_CHECK_EQUAL(r_algo.reactions[0].gamma, 2);
    r_algo.delete_reaction(0);
    BOOST_CHECK_EQUAL(r_algo.reactions.size(), 0ul);
  }

  // exception if deleting a non-existent particle
  BOOST_CHECK_THROW(r_algo.delete_particle(5), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(WidomInsertion_test) {
  using namespace ReactionEnsemble;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // check acceptance rate
  WidomInsertion r_algo(42);

  SingleReaction const ref(2., {0}, {1}, {1, 2}, {3, 4});
  r_algo.add_reaction(ref.gamma, ref.reactant_types, ref.reactant_coefficients,
                      ref.product_types, ref.product_coefficients);

  // exception if not enough particles
  BOOST_CHECK_THROW(r_algo.measure_excess_chemical_potential(0),
                    std::runtime_error);
}

int main(int argc, char **argv) {
  auto mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  Communication::init(mpi_env);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
