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
#include "script_interface/system/Leaf.hpp"

#include <utils/serialization/pack.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Interactions {

class NonBondedInteractions : public System::Leaf {
  using container_type =
      std::unordered_map<unsigned int,
                         std::shared_ptr<NonBondedInteractionHandle>>;

public:
  using key_type = typename container_type::key_type;
  using mapped_type = typename container_type::mapped_type;

private:
  container_type m_nonbonded_ia_params;
  std::shared_ptr<::InteractionsNonBonded> m_handle;
  std::shared_ptr<std::function<void()>> m_notify_cutoff_change;

  void do_construct(VariantMap const &params) override {
    m_handle = std::make_shared<::InteractionsNonBonded>();
    m_notify_cutoff_change = std::make_shared<std::function<void()>>([]() {});
  }

  void on_bind_system(::System::System &system) override {
    auto const max_type = m_handle->get_max_seen_particle_type();
    system.nonbonded_ias = m_handle;
    m_handle->make_particle_type_exist(max_type);
    m_handle->bind_system(m_system.lock());
    m_handle->on_non_bonded_ia_change();
    *m_notify_cutoff_change = [this]() {
      if (m_handle) {
        m_handle->on_non_bonded_ia_change();
      }
    };
  }

  void on_detach_system(::System::System &) override {
    *m_notify_cutoff_change = []() {};
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
    if (method == "reset") {
      if (not context()->is_head_node()) {
        return {};
      }
      auto const max_type = m_handle->get_max_seen_particle_type();
      auto const obj_params = VariantMap{{"notify", false}};
      for (int i = 0; i <= max_type; i++) {
        for (int j = 0; j <= i; j++) {
          auto const key = m_handle->get_ia_param_key(i, j);
          if (m_nonbonded_ia_params.contains(key)) {
            m_nonbonded_ia_params.at(key)->call_method("reset", obj_params);
          }
        }
      }
      get_system().on_non_bonded_ia_change();
      return {};
    }
    if (method == "get_handle") {
      auto ctx = context();
      auto const [type_min, type_max] = get_key(params.at("key"));
      if (type_max > m_handle->get_max_seen_particle_type()) {
        m_handle->make_particle_type_exist(type_max);
      }
      if (not ctx->is_head_node()) {
        return {};
      }
      auto const key = m_handle->get_ia_param_key(type_min, type_max);
      if (m_nonbonded_ia_params.contains(key)) {
        return m_nonbonded_ia_params.at(key);
      }
      auto so = std::dynamic_pointer_cast<NonBondedInteractionHandle>(
          ctx->make_shared("Interactions::NonBondedInteractionHandle", {}));
      m_nonbonded_ia_params[key] = so;
      call_method("internal_attach", {{"key", params.at("key")}, {"obj", so}});
      return so;
    }
    if (method == "internal_set_max_type") {
      m_handle->make_particle_type_exist(get_value<int>(params, "max_type"));
      return {};
    }
    if (method == "internal_attach") {
      auto so = std::dynamic_pointer_cast<NonBondedInteractionHandle>(
          get_value<ObjectRef>(params, "obj"));
      auto const [i, j] = get_key(params.at("key"));
      auto const cb_register =
          [this, i, j](std::shared_ptr<::IA_parameters> const &core_ia) {
            m_handle->set_ia_param(i, j, core_ia);
          };
      so->attach(cb_register, m_notify_cutoff_change);
      return {};
    }

    return {};
  }

  std::string get_internal_state() const override {
    auto const max_type = m_handle->get_max_seen_particle_type();
    std::vector<std::string> object_states;
    object_states.emplace_back(Utils::pack(max_type));
    for (int i = 0; i <= max_type; i++) {
      for (int j = 0; j <= i; j++) {
        auto const key = m_handle->get_ia_param_key(i, j);
        if (m_nonbonded_ia_params.contains(key)) {
          object_states.emplace_back(
              m_nonbonded_ia_params.at(key)->serialize());
        } else {
          object_states.emplace_back("");
        }
      }
    }

    return Utils::pack(object_states);
  }

  void set_internal_state(std::string const &state) override {
    auto const object_states = Utils::unpack<std::vector<std::string>>(state);
    auto const max_type = Utils::unpack<int>(object_states.front());
    call_method("internal_set_max_type", {{"max_type", max_type}});
    auto const end = object_states.end();
    auto it = object_states.begin() + 1;
    for (int i = 0; i <= max_type; i++) {
      for (int j = 0; j <= i; j++) {
        auto const key = m_handle->get_ia_param_key(i, j);
        auto const &buffer = *it;
        if (not buffer.empty()) {
          auto so = std::dynamic_pointer_cast<NonBondedInteractionHandle>(
              ObjectHandle::deserialize(buffer, *context()));
          m_nonbonded_ia_params[key] = so;
          call_method("internal_attach",
                      {{"key", std::vector<int>{{j, i}}}, {"obj", so}});
        }
        ++it;
        if (it == end) {
          break;
        }
      }
    }
  }
};

} // namespace Interactions
} // namespace ScriptInterface
