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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_ELECTROSTATICS_COULOMB_SCAFACOS_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_ELECTROSTATICS_COULOMB_SCAFACOS_HPP

#include "config.hpp"

#ifdef SCAFACOS

#include "Actor.hpp"

#include "core/MpiCallbacks.hpp"
#include "core/communication.hpp"
#include "core/electrostatics/scafacos.hpp"
#include "core/scafacos/ScafacosContextBase.hpp"

#include "script_interface/get_value.hpp"
#include "script_interface/scafacos/scafacos.hpp"

#include <iomanip>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Coulomb {

class CoulombScafacos : public Actor<CoulombScafacos, ::CoulombScafacos> {
  std::shared_ptr<boost::mpi::environment> m_mpi_env_lock;

public:
  CoulombScafacos() {
    add_parameters({
        {"method_name", AutoParameter::read_only,
         [this]() { return actor()->get_method(); }},
        {"method_params", AutoParameter::read_only,
         [this]() {
           auto const m_tuned_methods =
               std::set<std::string>{"ewald", "p2nfft", "p3m"};
           auto parameters_string = actor()->get_parameters();
           auto const method_name = actor()->get_method();
           auto const delegate = actor()->get_near_field_delegation();
           if (delegate and m_tuned_methods.count(method_name)) {
             auto const tuned_r_cut = actor()->get_r_cut();
             auto const field_name = method_name + "_r_cut";
             std::ostringstream serializer;
             serializer << std::scientific << std::setprecision(17);
             serializer << tuned_r_cut;
             auto const tuned_r_cut_string = serializer.str();
             auto const replacement =
                 "," + field_name + "," + tuned_r_cut_string;
             if (parameters_string.find(field_name) == field_name.npos) {
               parameters_string += replacement;
             } else {
               auto const field_pattern =
                   std::regex("," + field_name + ",[0-9eE\\-\\+\\.]+");
               parameters_string = std::regex_replace(
                   parameters_string, std::regex(field_pattern), replacement);
             }
           }
           return Scafacos::deserialize_parameters(parameters_string);
         }},
    });
  }

  ~CoulombScafacos() override {
    m_actor.reset();
    m_mpi_env_lock.reset();
  }

  void do_construct(VariantMap const &params) override {
    auto const method_name = get_value<std::string>(params, "method_name");
    auto const param_list = params.at("method_params");
    auto const prefactor = get_value<double>(params, "prefactor");

    context()->parallel_try_catch([&]() {
      ScafacosContextBase::sanity_check_method(method_name);
      auto const method_params = Scafacos::serialize_parameters(param_list);
      m_actor = make_coulomb_scafacos(method_name, method_params);
      actor()->set_prefactor(prefactor);
    });
    set_charge_neutrality_tolerance(params);
    // MPI communicator is needed to destroy the FFT plans
    m_mpi_env_lock = ::Communication::mpiCallbacksHandle()->share_mpi_env();
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "set_near_field_delegation") {
      auto const delegate = get_value<bool>(params, "delegate");
      context()->parallel_try_catch(
          [&]() { actor()->set_near_field_delegation(delegate); });
      return {};
    }
    if (name == "get_near_field_delegation") {
      return actor()->get_near_field_delegation();
    }
    if (name == "get_available_methods") {
      return make_vector_of_variants(Scafacos::available_methods());
    }
    return Actor<SIActorClass, CoreActorClass>::do_call_method(name, params);
  }
};

} // namespace Coulomb
} // namespace ScriptInterface

#endif // SCAFACOS
#endif
