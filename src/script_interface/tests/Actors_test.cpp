/*
 * Copyright (C) 2022 The ESPResSo project
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

#define BOOST_TEST_MODULE "Long-range actors test"
#define BOOST_TEST_DYN_LINK

#include "config/config.hpp"

#if defined(ELECTROSTATICS) or defined(DIPOLES) or defined(SCAFACOS) or        \
    defined(SCAFACOS_DIPOLES)

#include <boost/test/unit_test.hpp>

#include "script_interface/electrostatics/Actor.hpp"
#include "script_interface/electrostatics/Actor_impl.hpp"
#include "script_interface/magnetostatics/Actor.hpp"
#include "script_interface/magnetostatics/Actor_impl.hpp"
#include "script_interface/scafacos/scafacos.hpp"

#include "core/actor/visitors.hpp"
#include "core/electrostatics/coulomb.hpp"
#include "core/electrostatics/debye_hueckel.hpp"
#include "core/electrostatics/registration.hpp"
#include "core/magnetostatics/dipolar_direct_sum.hpp"
#include "core/magnetostatics/dipoles.hpp"
#include "core/magnetostatics/registration.hpp"

#include "core/communication.hpp"

#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

#include <utils/Vector.hpp>

#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace ScriptInterface {
namespace Coulomb {

#ifdef ELECTROSTATICS
struct MockDebyeHueckel : public Actor<MockDebyeHueckel, ::DebyeHueckel> {
  MockDebyeHueckel() = default;

  void do_construct(VariantMap const &params) override {
    m_actor = std::make_shared<CoreActorClass>(
        get_value<double>(params, "prefactor"),
        get_value<double>(params, "kappa"), get_value<double>(params, "r_cut"));
  }
};
#endif // ELECTROSTATICS

} // namespace Coulomb

namespace Dipoles {

#ifdef DIPOLES
struct MockDipolarDirectSum
    : public Actor<MockDipolarDirectSum, ::DipolarDirectSum> {
  MockDipolarDirectSum() = default;

  void do_construct(VariantMap const &params) override {
    m_actor = std::make_shared<CoreActorClass>(
        get_value<double>(params, "prefactor"),
        get_value<double>(params, "n_replicas"));
  }
};
#endif // DIPOLES

} // namespace Dipoles
} // namespace ScriptInterface

#ifdef ELECTROSTATICS
BOOST_AUTO_TEST_CASE(coulomb_actor) {
  auto constexpr tol = 100. * std::numeric_limits<double>::epsilon();
  ScriptInterface::Coulomb::MockDebyeHueckel actor;
  actor.do_construct({{"prefactor", 2.}, {"kappa", 3.}, {"r_cut", 4.}});
  // check const and non-const access
  BOOST_CHECK_CLOSE(actor.actor()->prefactor, 2., tol);
  BOOST_CHECK_CLOSE(std::as_const(actor).actor()->prefactor, 2., tol);
  // check visitors
  BOOST_CHECK(has_actor_of_type<::DebyeHueckel>(
      boost::optional<ElectrostaticsActor>(actor.actor())));
  BOOST_CHECK(not has_actor_of_type<decltype(actor)>(
      boost::optional<ElectrostaticsActor>(actor.actor())));
  BOOST_CHECK(is_already_stored(
      actor.actor(), boost::optional<ElectrostaticsActor>(actor.actor())));
  BOOST_CHECK(not is_already_stored(
      std::shared_ptr<::DebyeHueckel>{},
      boost::optional<ElectrostaticsActor>(actor.actor())));
  BOOST_CHECK(not is_already_stored(
      std::shared_ptr<decltype(actor)>{},
      boost::optional<ElectrostaticsActor>(actor.actor())));
}
#endif // ELECTROSTATICS

#ifdef DIPOLES
BOOST_AUTO_TEST_CASE(dipoles_actor) {
  auto constexpr tol = 100. * std::numeric_limits<double>::epsilon();
  n_nodes = 1;
  ScriptInterface::Dipoles::MockDipolarDirectSum actor;
  actor.do_construct({{"prefactor", 2.}, {"n_replicas", 3}});
  // check const and non-const access
  BOOST_CHECK_EQUAL(actor.actor()->n_replicas, 3);
  BOOST_CHECK_CLOSE(actor.actor()->prefactor, 2., tol);
  BOOST_CHECK_CLOSE(std::as_const(actor).actor()->prefactor, 2., tol);
}
#endif // DIPOLES

#if defined(SCAFACOS) or defined(SCAFACOS_DIPOLES)

BOOST_AUTO_TEST_CASE(scafacos_parameters_serialization) {
  using ScriptInterface::Variant;
  {
    auto const vec_int = std::vector<int>{1, 2, 3};
    auto const vec_var = std::vector<Variant>{4, 5, std::string("ab")};
    auto const parameters = ScriptInterface::Scafacos::serialize_parameters(
        Variant{std::unordered_map<std::string, Variant>{
            {"key", Variant{2}},
            {"vec_int", Variant{vec_int}},
            {"vec_var", Variant{vec_var}}}});
    auto const components =
        std::vector<std::string>{"key,2", "vec_int,1,2,3", "vec_var,4,5,ab"};
    for (auto const &item : components) {
      BOOST_CHECK(parameters.find(item) != item.npos);
    }
  }
}

BOOST_AUTO_TEST_CASE(scafacos_parameters_serialization_exceptions) {
  using ScriptInterface::Variant;
  {
    auto caught = false;
    try {
      ScriptInterface::Scafacos::serialize_parameters(Variant{5.});
    } catch (ScriptInterface::Exception const &) {
      caught = true;
    }
    BOOST_CHECK(caught);
  }
  {
    auto caught = false;
    try {
      auto const empty_map = std::unordered_map<int, Variant>{};
      ScriptInterface::Scafacos::serialize_parameters(Variant{empty_map});
    } catch (std::invalid_argument const &) {
      caught = true;
    }
    BOOST_CHECK(caught);
  }
  {
    auto caught = false;
    try {
      auto const invalid_map = std::unordered_map<std::string, Variant>{
          {"k", Variant{std::unordered_map<std::string, Variant>{}}}};
      ScriptInterface::Scafacos::serialize_parameters(Variant{invalid_map});
    } catch (std::runtime_error const &) {
      caught = true;
    }
    BOOST_CHECK(caught);
  }
}
#endif // SCAFACOS or SCAFACOS_DIPOLES

#else
int main(int argc, char **argv) {}
#endif // ELECTROSTATICS or DIPOLES or SCAFACOS or SCAFACOS_DIPOLES
