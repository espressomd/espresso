/*
 * Copyright (C) 2021-2023 The ESPResSo project
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
#define BOOST_TEST_MODULE ConstantpHEnsemble test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/reaction_methods/ConstantpHEnsemble.hpp"
#include "core/reaction_methods/SingleReaction.hpp"
#include "script_interface/LocalContext.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/reaction_methods/ConstantpHEnsemble.hpp"
#include "script_interface/reaction_methods/ExclusionRadius.hpp"

#include "core/EspressoSystemStandAlone.hpp"
#include "core/Particle.hpp"
#include "core/particle_node.hpp"
#include "core/unit_tests/ParticleFactory.hpp"

#include <boost/mpi.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <unordered_map>

using ::ReactionMethods::SingleReaction;

namespace espresso {
// ESPResSo system instance
static std::unique_ptr<EspressoSystemStandAlone> system;
} // namespace espresso

namespace ScriptInterface::Testing {
class ConstantpHEnsemble
    : public ScriptInterface::ReactionMethods::ConstantpHEnsemble {
public:
  auto calculate_acceptance_probability(
      SingleReaction const &reaction, double E_pot_diff,
      std::unordered_map<int, int> const &old_particle_numbers) const {
    auto const factorial_expr =
        ::ReactionMethods::calculate_factorial_expression_cpH(
            reaction, old_particle_numbers);
    auto const pH =
        std::dynamic_pointer_cast<::ReactionMethods::ConstantpHEnsemble>(RE())
            ->m_constant_pH;
    auto const ln_bf =
        E_pot_diff - reaction.nu_bar * RE()->kT * std::log(10.) *
                         (pH + reaction.nu_bar * std::log10(reaction.gamma));
    return factorial_expr * std::exp(-ln_bf / RE()->kT);
  }
};
} // namespace ScriptInterface::Testing

// Check the Monte Carlo algorithm where moves depend on the system
// configuration, energy and pH.
BOOST_FIXTURE_TEST_CASE(ConstantpHEnsemble_test, ParticleFactory) {
  using namespace ScriptInterface;
  auto constexpr tol = 100. * std::numeric_limits<double>::epsilon();

  Utils::Factory<ScriptInterface::ObjectHandle> factory;
  factory.register_new<Testing::ConstantpHEnsemble>(
      "Testing::ConstantpHEnsemble");
  factory.register_new<ScriptInterface::ReactionMethods::ExclusionRadius>(
      "ExclusionRadius");

  auto const comm = boost::mpi::communicator();
  auto const make_algo = [&factory,
                          &comm](int seed, double kT, double pH,
                                 double exclusion_range,
                                 std::unordered_map<int, double> const &radii) {
    auto ctx = std::make_shared<ScriptInterface::LocalContext>(factory, comm);
    VariantMap params{};
    params["seed"] = seed;
    params["kT"] = kT;
    params["constant_pH"] = pH;
    params["exclusion_range"] = exclusion_range;
    params["exclusion_radius_per_type"] = make_unordered_map_of_variants(radii);
    params["exclusion"] = ctx->make_shared_local("ExclusionRadius", params);
    auto &&sp = ctx->make_shared_local("Testing::ConstantpHEnsemble", params);
    return std::dynamic_pointer_cast<Testing::ConstantpHEnsemble>(sp);
  };

  auto const constant_pH = 1.;
  auto &&r_algo_si = make_algo(42, 20., constant_pH, 0., {});
  auto &r_algo = *r_algo_si->RE();
  r_algo.non_interacting_type = 5;

  // create a reaction A -> B + C
  int const type_A = 0;
  int const type_B = 1;
  int const type_C = 2;
  SingleReaction const reaction(2., {type_A}, {1}, {type_B, type_C}, {1, 1});

  // check acceptance probability
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        // system contains i x A, j x B, and k x C
        auto const p_numbers =
            std::unordered_map<int, int>{{type_A, i}, {type_B, j}, {type_C, k}};
        auto const energy = -static_cast<double>(i + 1);
        auto const f_expr =
            calculate_factorial_expression_cpH(reaction, p_numbers);
        // bf = f_expr * exp(- E / kT + nu_bar * log(10) * (pH - nu_bar * pKa))
        auto const acceptance_ref =
            f_expr * std::exp(-energy / r_algo.kT +
                              std::log(10.) *
                                  (constant_pH + std::log10(reaction.gamma)));
        auto const acceptance = r_algo_si->calculate_acceptance_probability(
            reaction, energy, p_numbers);
        BOOST_CHECK_CLOSE(acceptance, acceptance_ref, 5. * tol);
      }
    }
  }
}

int main(int argc, char **argv) {
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
