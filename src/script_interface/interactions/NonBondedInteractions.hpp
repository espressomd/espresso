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

#pragma once

#include "NonBondedInteraction.hpp"

#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/system/System.hpp"

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
  decltype(System::System::nonbonded_ias) m_nonbonded_ias;

  auto make_interaction(int i, int j) {
    assert(j <= i);
    auto const types = std::vector<int>{{i, j}};
    return std::dynamic_pointer_cast<NonBondedInteractionHandle>(
        context()->make_shared_local("Interactions::NonBondedInteractionHandle",
                                     {{"_types", Variant{types}}}));
  }

  void reset() {
    auto const max_type = m_nonbonded_ias->get_max_seen_particle_type();
    for (int i = 0; i <= max_type; i++) {
      for (int j = 0; j <= i; j++) {
        auto const key = m_nonbonded_ias->get_ia_param_key(i, j);
        m_nonbonded_ias->set_ia_param(i, j,
                                      std::make_shared<::IA_parameters>());
        m_nonbonded_ia_params[key] = make_interaction(i, j);
      }
    }
    System::get_system().on_non_bonded_ia_change();
  }

  void do_construct(VariantMap const &) override {
    m_nonbonded_ias = System::get_system().nonbonded_ias;
    auto const max_type = m_nonbonded_ias->get_max_seen_particle_type();
    for (int i = 0; i <= max_type; i++) {
      for (int j = 0; j <= i; j++) {
        auto const key = m_nonbonded_ias->get_ia_param_key(i, j);
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
      return {std::max(types[0], types[1]), std::min(types[0], types[1])};
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
      return Variant{m_nonbonded_ias->get_max_seen_particle_type() + 1};
    }
    if (method == "reset") {
      reset();
      return {};
    }
    if (method == "insert") {
      auto const [type_min, type_max] = get_key(params.at("key"));
      m_nonbonded_ias->make_particle_type_exist(type_max);
      auto const key = m_nonbonded_ias->get_ia_param_key(type_min, type_max);
      auto obj_ptr = get_value<std::shared_ptr<NonBondedInteractionHandle>>(
          params.at("object"));
      m_nonbonded_ias->set_ia_param(type_min, type_max, obj_ptr->get_ia());
      m_nonbonded_ia_params[key] = obj_ptr;
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
