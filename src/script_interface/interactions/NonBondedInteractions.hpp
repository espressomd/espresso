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
      std::unordered_map<unsigned int,
                         std::shared_ptr<NonBondedInteractionHandle>>;

public:
  using key_type = typename container_type::key_type;
  using mapped_type = typename container_type::mapped_type;

private:
  container_type m_nonbonded_ia_params;

  auto make_interaction(int i, int j) {
    assert(i <= j);
    auto const types = std::vector<int>{{i, j}};
    return std::dynamic_pointer_cast<NonBondedInteractionHandle>(
        context()->make_shared_local("Interactions::NonBondedInteractionHandle",
                                     {{"_types", Variant{types}}}));
  }

  void reset() {
    auto const size = ::max_seen_particle_type;
    for (int i = 0; i < size; i++) {
      for (int j = i; j < size; j++) {
        auto const key = get_ia_param_key(i, j);
        ::nonbonded_ia_params[i] = std::make_shared<::IA_parameters>();
        m_nonbonded_ia_params[key] = make_interaction(i, j);
      }
    }
    on_non_bonded_ia_change();
  }

  void do_construct(VariantMap const &params) override {
    auto const size = ::max_seen_particle_type;
    make_particle_type_exist(size);
    for (int i = 0; i < size; i++) {
      for (int j = i; j < size; j++) {
        auto const key = get_ia_param_key(i, j);
        m_nonbonded_ia_params[key] = make_interaction(i, j);
      }
    }
  }

  std::pair<int, int> get_key(Variant const &key) const {
    try {
      auto const types = get_value<std::vector<int>>(key);
      if (types.size() != 2ul or types[0] < 0 or types[1] < 0) {
        throw Exception("need two particle types");
      }
      return {std::min(types[0], types[1]), std::max(types[0], types[1])};
    } catch (...) {
      if (context()->is_head_node()) {
        throw std::invalid_argument(
            "NonBondedInteractions[] expects two particle types as indices");
      }
      throw;
    }
  }

protected:
  Variant do_call_method(std::string const &method,
                         VariantMap const &params) override {
    if (method == "get_n_types") {
      return Variant{::max_seen_particle_type};
    }
    if (method == "reset") {
      reset();
      return {};
    }
    if (method == "insert") {
      auto const [type_min, type_max] = get_key(params.at("key"));
      make_particle_type_exist(type_max);
      auto const key = get_ia_param_key(type_min, type_max);
      auto obj_ptr = get_value<std::shared_ptr<NonBondedInteractionHandle>>(
          params.at("object"));
      ::nonbonded_ia_params[key] = obj_ptr->get_ia();
      m_nonbonded_ia_params[key] = obj_ptr;
      on_non_bonded_ia_change();
      return {};
    }
    if (method == "check_key") {
      static_cast<void>(get_key(params.at("key")));
      return {};
    }

    return {};
  }
};

} // namespace Interactions
} // namespace ScriptInterface

#endif
