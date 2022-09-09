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

#ifndef SCRIPT_INTERFACE_INTERACTIONS_NONBONDED_INTERACTIONS_HPP
#define SCRIPT_INTERFACE_INTERACTIONS_NONBONDED_INTERACTIONS_HPP

#include "NonBondedInteraction.hpp"

#include "core/event.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "script_interface/ObjectMap.hpp"
#include "script_interface/ScriptInterface.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Interactions {

class NonBondedInteractions : public ObjectHandle {
  using container_type =
      std::unordered_map<int, std::shared_ptr<NonBondedInteractionHandle>>;

  auto make_interaction(int i, int j) {
    assert(i <= j);
    auto const types = std::vector<int>{{i, j}};
    return std::dynamic_pointer_cast<NonBondedInteractionHandle>(
        context()->make_shared_local("Interactions::NonBondedInteractionHandle",
                                     {{"_types", Variant{types}}}));
  }

public:
  using key_type = typename container_type::key_type;
  using mapped_type = typename container_type::mapped_type;

  void reset() {
    auto const size = ::max_seen_particle_type;
    for (int i = 0; i < size; i++) {
      for (int j = i; j < size; j++) {
        auto const key = Utils::upper_triangular(i, j, size);
        ::nonbonded_ia_params[i] = std::make_shared<::IA_parameters>();
        m_nonbonded_ia_params[key] = make_interaction(i, j);
      }
    }
    on_non_bonded_ia_change();
  }

  void do_construct(VariantMap const &params) override {
    auto const size = ::max_seen_particle_type;
    {
      // when reloading from a checkpoint file, need to resize IA lists
      auto const new_size = ::max_seen_particle_type;
      auto const n_pairs = new_size * (new_size + 1) / 2;
      ::nonbonded_ia_params.resize(n_pairs);
      for (auto &ia_params : ::nonbonded_ia_params) {
        if (ia_params == nullptr) {
          ia_params = std::make_shared<::IA_parameters>();
        }
      }
    }
    for (int i = 0; i < size; i++) {
      for (int j = i; j < size; j++) {
        auto const key = Utils::upper_triangular(i, j, size);
        m_nonbonded_ia_params[key] = make_interaction(i, j);
      }
    }
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_n_types") {
      return Variant{::max_seen_particle_type};
    }
    if (name == "reset") {
      reset();
      return {};
    }
    if (name == "insert") {
      auto const types = get_value<std::vector<int>>(params.at("key"));
      auto const key = get_ia_param_key(std::min(types[0], types[1]),
                                        std::max(types[0], types[1]));
      auto obj_ptr = get_value<std::shared_ptr<NonBondedInteractionHandle>>(
          params.at("object"));
      ::nonbonded_ia_params[key] = obj_ptr->get_ia();
      m_nonbonded_ia_params[key] = obj_ptr;
      on_non_bonded_ia_change();
      return {};
    }

    return {};
  }

private:
  // disable serialization: pickling done by the python interface
  std::string get_internal_state() const override { return {}; }
  void set_internal_state(std::string const &state) override {}
  container_type m_nonbonded_ia_params;
};

} // namespace Interactions
} // namespace ScriptInterface

#endif
