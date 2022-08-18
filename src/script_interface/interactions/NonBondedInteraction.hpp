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

/** @file
 *  The ScriptInterface counterparts of the non-bonded interactions parameters
 *  structs from the core are defined here.
 */

#ifndef SCRIPT_INTERFACE_INTERACTIONS_NONBONDED_INTERACTION_HPP
#define SCRIPT_INTERFACE_INTERACTIONS_NONBONDED_INTERACTION_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"

#include "core/event.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Interactions {

template <class CoreIA>
class InteractionPotentialInterface
    : public AutoParameters<InteractionPotentialInterface<CoreIA>> {
  std::array<int, 2> m_types = {-1, -1};

public:
  using CoreInteraction = CoreIA;

protected:
  using AutoParameters<InteractionPotentialInterface<CoreIA>>::context;
  std::shared_ptr<CoreInteraction> m_ia_si;
  virtual CoreInteraction IA_parameters::*get_ptr_offset() const = 0;
  virtual void make_new_instance(VariantMap const &params) = 0;

  template <typename T>
  auto make_autoparameter(T CoreInteraction::*ptr, char const *name) {
    return AutoParameter{name, AutoParameter::read_only,
                         [this, ptr]() { return m_ia_si.get()->*ptr; }};
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "set_params") {
      context()->parallel_try_catch(
          [this, &params]() { make_new_instance(params); });
      if (m_types[0] != -1) {
        copy_si_to_core();
        on_non_bonded_ia_change();
      }
      return {};
    }
    if (name == "bind_types") {
      auto types = get_value<std::vector<int>>(params, "_types");
      if (types[0] > types[1]) {
        std::swap(types[0], types[1]);
      }
      if (m_types[0] == -1 or
          (m_types[0] == types[0] and m_types[1] == types[1])) {
        m_types[0] = types[0];
        m_types[1] = types[1];
      } else {
        context()->parallel_try_catch([this]() {
          throw std::runtime_error(
              "Non-bonded interaction is already bound to interaction pair [" +
              std::to_string(m_types[0]) + ", " + std::to_string(m_types[1]) +
              "]");
        });
      }
      return {};
    }
    return {};
  }

  void do_construct(VariantMap const &params) final {
    if (params.empty()) {
      m_ia_si = std::make_shared<CoreInteraction>();
    } else if (params.count("_types") != 0) {
      do_call_method("bind_types", params);
      m_ia_si = std::make_shared<CoreInteraction>();
      copy_core_to_si();
    } else {
      context()->parallel_try_catch(
          [this, &params]() { make_new_instance(params); });
    }
  }

  auto const &get_ia() const { return *m_ia_si; }

  void copy_si_to_core() {
    assert(m_ia_si != nullptr);
    auto const key = get_ia_param_key(m_types[0], m_types[1]);
    assert(key < ::nonbonded_ia_params.size());
    ::nonbonded_ia_params[key].get()->*get_ptr_offset() = *m_ia_si;
    ::old_nonbonded_ia_params[key].*get_ptr_offset() = *m_ia_si;
  }

  void copy_core_to_si() {
    assert(m_ia_si != nullptr);
    auto const key = get_ia_param_key(m_types[0], m_types[1]);
    assert(key < ::nonbonded_ia_params.size());
    *m_ia_si = ::nonbonded_ia_params[key].get()->*get_ptr_offset();
  }
};

class NonBondedInteractionHandle
    : public AutoParameters<NonBondedInteractionHandle> {
  std::array<int, 2> m_types = {-1, -1};
  std::shared_ptr<::IA_parameters> m_interaction;

  template <class T>
  auto make_autoparameter(std::shared_ptr<T> &member, const char *key) const {
    auto const setter = [this, &member](Variant const &v) {
      member = get_value<std::shared_ptr<T>>(v);
      if (m_types[0] != -1) {
        auto const types = Variant{std::vector<int>{{m_types[0], m_types[1]}}};
        member->do_call_method("bind_types", VariantMap{{"_types", types}});
        member->copy_si_to_core();
        on_non_bonded_ia_change();
      }
    };
    return AutoParameter{key, setter, [&member]() { return member; }};
  }

public:
  NonBondedInteractionHandle() {
    add_parameters({
    });
  }

private:
  template <class T>
  void set_member(std::shared_ptr<T> &member, std::string key,
                  std::string so_name, VariantMap const &params) {
    auto const ia_types = VariantMap{{"_types", params.at("_types")}};
    if (params.count(key) != 0) {
      member = get_value<std::shared_ptr<T>>(params.at(key));
      member->do_call_method("bind_types", ia_types);
      member->copy_si_to_core();
    } else {
      auto so_object = context()->make_shared_local(so_name, ia_types);
      member = std::dynamic_pointer_cast<T>(so_object);
      member->copy_core_to_si();
    }
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_types") {
      return std::vector<int>{{m_types[0], m_types[1]}};
    }
    return {};
  }

  void do_construct(VariantMap const &params) override {
    assert(params.count("_types") != 0);
    auto const types = get_value<std::vector<int>>(params.at("_types"));
    m_types[0] = std::min(types[0], types[1]);
    m_types[1] = std::max(types[0], types[1]);
    // make type exist
    mpi_realloc_ia_params_local(m_types[1] + 1);
    // create interface objects
    auto const key = get_ia_param_key(m_types[0], m_types[1]);
    m_interaction = ::nonbonded_ia_params[key];
  }

  auto get_ia() const { return m_interaction; }
};

} // namespace Interactions
} // namespace ScriptInterface

#endif
