/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "config/config.hpp"

#ifdef COLLISION_DETECTION

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/interactions/BondedInteraction.hpp"
#include "script_interface/interactions/BondedInteractions.hpp"
#include "script_interface/system/Leaf.hpp"
#include "script_interface/system/System.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/collision.hpp"

#include <cassert>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace ScriptInterface {
namespace CollisionDetection {

class CollisionDetection
    : public AutoParameters<CollisionDetection, System::Leaf> {
  std::shared_ptr<::CollisionDetection> m_handle;
  std::unique_ptr<VariantMap> m_params;
  std::weak_ptr<Interactions::BondedInteractions> m_bonded_ias;

  std::unordered_map<CollisionModeType, std::string> const cd_mode_to_name = {
      {CollisionModeType::OFF, "off"},
      {CollisionModeType::BIND_CENTERS, "bind_centers"},
      {CollisionModeType::BIND_VS, "bind_at_point_of_collision"},
      {CollisionModeType::GLUE_TO_SURF, "glue_to_surface"},
  };
  std::unordered_map<std::string, CollisionModeType> cd_name_to_mode;
  std::unordered_map<CollisionModeType,
                     std::vector<std::string>> const cd_mode_to_parameters = {
      {CollisionModeType::OFF, {"mode"}},
      {CollisionModeType::BIND_CENTERS, {"mode", "bond_centers", "distance"}},
      {CollisionModeType::BIND_VS,
       {"mode", "bond_centers", "bond_vs", "part_type_vs", "distance",
        "vs_placement"}},
      {CollisionModeType::GLUE_TO_SURF,
       {"mode", "bond_centers", "bond_vs", "part_type_vs",
        "part_type_to_be_glued", "part_type_to_attach_vs_to",
        "part_type_after_glueing", "distance",
        "distance_glued_particle_to_vs"}},
  };

  auto find_bond_id(Variant const &v) const {
    auto &system = get_system();
    if (is_type<int>(v)) {
      auto const bond_id = get_value<int>(v);
      std::optional<int> retval = std::nullopt;
      if (system.bonded_ias->contains(bond_id)) {
        retval = bond_id;
      }
      return retval;
    }
    auto obj = get_value<std::shared_ptr<Interactions::BondedInteraction>>(v);
    return system.bonded_ias->find_bond_id(obj->bonded_ia());
  }

public:
  CollisionDetection() {
    for (auto const &kv : cd_mode_to_name) {
      cd_name_to_mode[kv.second] = kv.first;
    }
    add_parameters({
        {"mode",
         [this](Variant const &v) {
           auto const name = get_value<std::string>(v);
           check_mode_name(name);
           m_handle->mode = cd_name_to_mode.at(name);
         },
         [this]() { return cd_mode_to_name.at(m_handle->mode); }},
        {"bond_centers",
         [this](Variant const &v) {
           auto const bond_id = find_bond_id(v);
           if (not bond_id) {
             throw std::invalid_argument("Bond in parameter 'bond_centers' was "
                                         "not added to the system");
           }
           m_handle->bond_centers = bond_id.value();
         },
         [this]() { return m_handle->bond_centers; }},
        {"bond_vs",
         [this](Variant const &v) {
           auto const bond_id = find_bond_id(v);
           if (not bond_id) {
             throw std::invalid_argument(
                 "Bond in parameter 'bond_vs' was not added to the system");
           }
           m_handle->bond_vs = bond_id.value();
         },
         [this]() { return m_handle->bond_vs; }},
        {"distance",
         [this](Variant const &v) {
           m_handle->distance = get_value<double>(v);
         },
         [this]() { return m_handle->distance; }},
        {"distance_glued_particle_to_vs",
         [this](Variant const &v) {
           m_handle->dist_glued_part_to_vs = get_value<double>(v);
         },
         [this]() { return m_handle->dist_glued_part_to_vs; }},
        {"vs_placement",
         [this](Variant const &v) {
           m_handle->vs_placement = get_value<double>(v);
         },
         [this]() { return m_handle->vs_placement; }},
        {"part_type_vs",
         [this](Variant const &v) {
           m_handle->vs_particle_type = get_value<int>(v);
         },
         [this]() { return m_handle->vs_particle_type; }},
        {"part_type_to_be_glued",
         [this](Variant const &v) {
           m_handle->part_type_to_be_glued = get_value<int>(v);
         },
         [this]() { return m_handle->part_type_to_be_glued; }},
        {"part_type_to_attach_vs_to",
         [this](Variant const &v) {
           m_handle->part_type_to_attach_vs_to = get_value<int>(v);
         },
         [this]() { return m_handle->part_type_to_attach_vs_to; }},
        {"part_type_after_glueing",
         [this](Variant const &v) {
           m_handle->part_type_after_glueing = get_value<int>(v);
         },
         [this]() { return m_handle->part_type_after_glueing; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
    if (params.empty()) {
      (*m_params)["mode"] = std::string("off");
    } else {
      // Assume we are reloading from a checkpoint file.
      // This else branch can be removed once issue #4483 is fixed.
      m_params = std::make_unique<VariantMap>();
      auto const name = get_value<std::string>(params, "mode");
      check_mode_name(name);
      for (auto const &name :
           cd_mode_to_parameters.at(cd_name_to_mode.at(name))) {
        (*m_params)[name] = params.at(name);
      }
      (*m_params)["mode"] = params.at("mode");
    }
  }

  Variant do_call_method(const std::string &name,
                         const VariantMap &params) override {
    if (name == "set_params") {
      context()->parallel_try_catch([this, &params]() {
        auto const backup = std::make_shared<::CollisionDetection>(*m_handle);
        auto &system = get_system();
        try {
          // check provided parameters
          check_input_parameters(params);
          // set parameters
          for (auto const &kv : params) {
            do_set_parameter(get_value<std::string>(kv.first), kv.second);
          }
          // sanitize parameters and calculate derived parameters
          m_handle->initialize();
          return none;
        } catch (...) {
          // restore original parameters and re-throw exception
          m_handle = system.collision_detection = backup;
          m_handle->initialize();
          throw;
        }
      });
    }
    if (name == "params_for_mode") {
      auto const name = get_value<std::string>(params, "mode");
      check_mode_name(name);
      auto const mode = cd_name_to_mode.at(name);
      return make_vector_of_variants(cd_mode_to_parameters.at(mode));
    }
    if (name == "get_bond_by_id") {
      if (not context()->is_head_node()) {
        return {};
      }
      return m_bonded_ias.lock()->call_method("get_bond", params);
    }
    return none;
  }

  void attach(std::weak_ptr<Interactions::BondedInteractions> bonded_ias) {
    m_bonded_ias = bonded_ias;
  }

private:
  void check_mode_name(std::string const &name) const {
    if (not cd_name_to_mode.contains(name)) {
      throw std::invalid_argument("Unknown collision mode '" + name + "'");
    }
  }

  void on_bind_system(::System::System &system) override {
    m_handle = system.collision_detection;
    m_handle->bind_system(m_system.lock());
    if (not m_params->empty()) {
      do_call_method("set_params", *m_params);
    }
    m_params.reset();
  }

  void check_input_parameters(VariantMap const &params) const {
    auto const name = get_value<std::string>(params, "mode");
    check_mode_name(name);
    auto const mode = cd_name_to_mode.at(name);
    auto const expected_params = cd_mode_to_parameters.at(mode);
    auto const expected_param_names =
        std::set<std::string>{expected_params.begin(), expected_params.end()};
    std::set<std::string> input_parameter_names = {};
    for (auto const &kv : params) {
      auto const &param_name = kv.first;
      if (not expected_param_names.contains(param_name)) {
        // throw if parameter is unknown
        std::ignore = get_parameter(param_name);
        // throw if parameter is known but doesn't match the selected mode
        throw std::runtime_error("Parameter '" + param_name + "' is not " +
                                 "required for mode '" + name + "'");
      }
      input_parameter_names.insert(param_name);
    }
    for (auto const &param_name : expected_param_names) {
      if (not input_parameter_names.contains(param_name)) {
        throw std::runtime_error("Parameter '" + param_name + "' is " +
                                 "required for mode '" + name + "'");
      }
    }
  }
};

} /* namespace CollisionDetection */
} /* namespace ScriptInterface */

#endif // COLLISION_DETECTION
