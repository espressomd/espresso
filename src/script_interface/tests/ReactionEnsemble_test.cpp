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
#define BOOST_TEST_MODULE ReactionEnsemble test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/reaction_methods/ReactionEnsemble.hpp"
#include "core/reaction_methods/SingleReaction.hpp"

#include "script_interface/LocalContext.hpp"
#include "script_interface/None.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/reaction_methods/ReactionEnsemble.hpp"

#include "core/Particle.hpp"
#include "core/cell_system/CellStructureType.hpp"
#include "core/communication.hpp"
#include "core/particle_node.hpp"
#include "core/unit_tests/ParticleFactory.hpp"

#include <boost/mpi.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using ::ReactionMethods::SingleReaction;

namespace espresso {
// ESPResSo system instance
static std::shared_ptr<System::System> system;
} // namespace espresso

namespace ScriptInterface::Testing {
class ReactionEnsemble
    : public ScriptInterface::ReactionMethods::ReactionEnsemble {
public:
  auto calculate_acceptance_probability(
      SingleReaction const &reaction, double E_pot_diff,
      std::unordered_map<int, int> const &old_particle_numbers) const {
    auto const factorial_expr =
        ::ReactionMethods::calculate_factorial_expression(reaction,
                                                          old_particle_numbers);
    return std::pow(RE()->get_volume(), reaction.nu_bar) * reaction.gamma *
           factorial_expr * std::exp(-E_pot_diff / RE()->kT);
  }
};
} // namespace ScriptInterface::Testing

// Check the Monte Carlo algorithm where moves depend on the system
// configuration and energy.
BOOST_FIXTURE_TEST_CASE(ReactionEnsemble_test, ParticleFactory) {
  using namespace ScriptInterface;
  auto constexpr tol = 100. * std::numeric_limits<double>::epsilon();

  Utils::Factory<ScriptInterface::ObjectHandle> factory;
  factory.register_new<Testing::ReactionEnsemble>("Testing::ReactionEnsemble");
  factory.register_new<ScriptInterface::ReactionMethods::SingleReaction>(
      "SingleReaction");

  auto const comm = boost::mpi::communicator();
  auto const make_algo = [&factory,
                          &comm](int seed, double kT, double exclusion_range,
                                 std::unordered_map<int, double> const &radii) {
    using namespace ScriptInterface;
    auto ctx = std::make_shared<ScriptInterface::LocalContext>(factory, comm);
    VariantMap params{};
    params["seed"] = seed;
    params["kT"] = kT;
    params["exclusion_range"] = exclusion_range;
    params["exclusion_radius_per_type"] = make_unordered_map_of_variants(radii);
    auto &&sp = ctx->make_shared_local("Testing::ReactionEnsemble", params);
    return std::dynamic_pointer_cast<Testing::ReactionEnsemble>(sp);
  };
  auto const make_reaction =
      [&factory, &comm](double gamma, std::vector<int> const &reactant_types,
                        std::vector<int> const &reactant_coefficients,
                        std::vector<int> const &product_types,
                        std::vector<int> const &product_coefficients) {
        using namespace ScriptInterface;
        auto ctx =
            std::make_shared<ScriptInterface::LocalContext>(factory, comm);
        VariantMap params{};
        params["gamma"] = gamma;
        params["reactant_types"] = make_vector_of_variants(reactant_types);
        params["reactant_coefficients"] =
            make_vector_of_variants(reactant_coefficients);
        params["product_types"] = make_vector_of_variants(product_types);
        params["product_coefficients"] =
            make_vector_of_variants(product_coefficients);
        auto &&si_obj = ctx->make_shared_local("SingleReaction", params);
        return std::dynamic_pointer_cast<
            ScriptInterface::ReactionMethods::SingleReaction>(si_obj);
      };

  // check basic interface
  {
    auto &&r_algo_si = make_algo(42, 20., 0., {});
    auto &r_algo = *r_algo_si->RE();
    r_algo.set_volume(10.);
    r_algo.non_interacting_type = 5;

    // create a reaction A -> 3 B + 4 C
    int const type_A = 0;
    int const type_B = 1;
    int const type_C = 2;
    SingleReaction const reaction(2., {type_A}, {1}, {type_B, type_C}, {3, 4});

    // check acceptance probability
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
          // system contains i x A, j x B, and k x C
          auto const p_numbers = std::unordered_map<int, int>{
              {type_A, i}, {type_B, j}, {type_C, k}};
          auto const energy = static_cast<double>(i + 1);
          auto const f_expr =
              calculate_factorial_expression(reaction, p_numbers);
          // acceptance = V^{nu_bar} * gamma * f_expr * exp(- E / T)
          auto const acceptance_ref = std::pow(r_algo.volume, reaction.nu_bar) *
                                      reaction.gamma * f_expr *
                                      std::exp(-energy / r_algo.kT);
          auto const acceptance = r_algo_si->calculate_acceptance_probability(
              reaction, energy, p_numbers);
          BOOST_CHECK_CLOSE(acceptance, acceptance_ref, 5 * tol);
        }
      }
    }
  }

  // check that the system energy is updated after a successful reaction
  {
    auto &&r_algo_si = make_algo(42, 1., 0., {});
    auto &r_algo = *r_algo_si->RE();
    r_algo.set_volume(1.);
    r_algo.non_interacting_type = 5;

    // create a generic identity exchange reaction D <-> E
    int const type_D = 0;
    int const type_E = 1;

    r_algo.charges_of_types[type_D] = 0;
    r_algo.charges_of_types[type_E] = 0;

    // track particles
    init_type_map(type_D);
    init_type_map(type_E);

    // set up a reaction completely shifted to product creation
    auto const gamma = 1e100;
    auto &&reaction_si = make_reaction(gamma, {type_D}, {1}, {type_E}, {1});
    auto &reaction = *reaction_si->get_reaction();
    r_algo_si->do_call_method("add_reaction", {{"reaction", reaction_si}});
    auto const reaction_id = 0;

    // resize system box
    auto const box_l = Utils::Vector3d{0.5, 0.4, 0.7};
    espresso::system->set_box_l(box_l);

    // without reactants, no reaction will take place
    auto const result = r_algo.create_new_trial_state(reaction_id);
    BOOST_REQUIRE(not result.has_value());

    // the reaction was updated
    BOOST_CHECK_EQUAL(reaction.tried_moves, 1);
    BOOST_CHECK_EQUAL(reaction.accepted_moves, 0);

    // add reactants in the system
    auto const ref_position = Utils::Vector3d{0.1, 0.2, 0.3};
    create_particle(ref_position, 1, type_D);
    create_particle(ref_position, 2, type_D);

    {
      // for an ideal system with gamma ~ inf, the reaction is always accepted;
      // the potential energy of the new state has to be 0 since it is an ideal
      // system
      auto const energy_ref = 0.;

      auto const result = r_algo.create_new_trial_state(reaction_id);
      BOOST_REQUIRE(result.has_value());
      auto const energy_move = *result;

      // verify bookkeeping
      auto const &bookkeeping = r_algo.get_old_system_state();
      BOOST_REQUIRE_EQUAL(bookkeeping.changed.size(), 1ul);
      BOOST_REQUIRE_EQUAL(bookkeeping.created.size(), 0ul);
      BOOST_REQUIRE_EQUAL(bookkeeping.hidden.size(), 0ul);
      BOOST_REQUIRE_EQUAL(bookkeeping.moved.size(), 0ul);
      BOOST_REQUIRE_EQUAL(bookkeeping.reaction_id, reaction_id);
      BOOST_REQUIRE_EQUAL(bookkeeping.old_particle_numbers.at(type_D), 2);
      BOOST_REQUIRE_EQUAL(bookkeeping.old_particle_numbers.at(type_E), 0);
      BOOST_REQUIRE_EQUAL(std::get<1>(bookkeeping.changed[0]), type_D);

      auto const bf = r_algo_si->calculate_acceptance_probability(
          reaction, energy_move, {{type_D, 1}, {type_E, 0}});

      auto const energy_end = r_algo.make_reaction_mc_move_attempt(
          reaction_id, bf, 0., energy_move);
      BOOST_CHECK_CLOSE(energy_end, energy_ref, tol);

      // verify bookkeeping was cleared
      BOOST_CHECK_THROW(r_algo.get_old_system_state(), std::runtime_error);

      // the reaction was updated
      BOOST_CHECK_EQUAL(reaction.tried_moves, 2);
      BOOST_CHECK_EQUAL(reaction.accepted_moves, 1);
    }
    {
      // attempt a second reaction
      auto const result = r_algo.create_new_trial_state(reaction_id);
      BOOST_REQUIRE(result.has_value());

      // verify bookkeeping
      auto const &bookkeeping = r_algo.get_old_system_state();
      BOOST_REQUIRE_EQUAL(bookkeeping.changed.size(), 1ul);
      BOOST_REQUIRE_EQUAL(bookkeeping.created.size(), 0ul);
      BOOST_REQUIRE_EQUAL(bookkeeping.hidden.size(), 0ul);
      BOOST_REQUIRE_EQUAL(bookkeeping.moved.size(), 0ul);
      BOOST_REQUIRE_EQUAL(bookkeeping.reaction_id, reaction_id);
      BOOST_REQUIRE_EQUAL(bookkeeping.old_particle_numbers.at(type_D), 1);
      BOOST_REQUIRE_EQUAL(bookkeeping.old_particle_numbers.at(type_E), 1);
      BOOST_REQUIRE_EQUAL(std::get<1>(bookkeeping.changed[0]), type_D);

      // force move to be rejected
      auto const energy_reject =
          r_algo.make_reaction_mc_move_attempt(reaction_id, 0., 0.2, 0.1);
      BOOST_CHECK_CLOSE(energy_reject, 0.2, tol);

      // the reaction was updated
      BOOST_CHECK_EQUAL(reaction.tried_moves, 3);
      BOOST_CHECK_EQUAL(reaction.accepted_moves, 1);

      // verify bookkeeping was cleared
      BOOST_CHECK_THROW(r_algo.get_old_system_state(), std::runtime_error);
    }
  }
}

int main(int argc, char **argv) {
  auto const mpi_handle = MpiContainerUnitTest(argc, argv);
  espresso::system = System::System::create();
  espresso::system->set_cell_structure_topology(CellStructureType::REGULAR);
  ::System::set_system(espresso::system);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
