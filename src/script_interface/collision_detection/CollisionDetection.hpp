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

#ifndef SCRIPT_INTERFACE_COLLISION_DETECTION_COLLISION_DETECTION_HPP
#define SCRIPT_INTERFACE_COLLISION_DETECTION_COLLISION_DETECTION_HPP

#include "config/config.hpp"

#ifdef COLLISION_DETECTION

#include "script_interface/ScriptInterface.hpp"

#include "core/collision.hpp"

#include <cassert>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace ScriptInterface {
namespace CollisionDetection {

class CollisionDetection : public AutoParameters<CollisionDetection> {
  std::unordered_map<CollisionModeType, std::string> const cd_mode_to_name = {
      {CollisionModeType::OFF, "off"},
      {CollisionModeType::BIND_CENTERS, "bind_centers"},
      {CollisionModeType::BIND_VS, "bind_at_point_of_collision"},
      {CollisionModeType::GLUE_TO_SURF, "glue_to_surface"},
      {CollisionModeType::BIND_THREE_PARTICLES, "bind_three_particles"},
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
      {CollisionModeType::BIND_THREE_PARTICLES,
       {"mode", "bond_centers", "distance", "bond_three_particles",
        "three_particle_binding_angle_resolution"}},
  };

public:
  CollisionDetection() {
    for (auto const &kv : cd_mode_to_name) {
      cd_name_to_mode[kv.second] = kv.first;
    }
    add_parameters(
        {{"mode",
          [this](Variant const &v) {
            auto const name = get_value<std::string>(v);
            check_mode_name(name);
            collision_params.mode = cd_name_to_mode.at(name);
          },
          [this]() { return cd_mode_to_name.at(collision_params.mode); }},

         {"bond_centers", collision_params.bond_centers},
         {"bond_vs", collision_params.bond_vs},
         {"bond_three_particles", collision_params.bond_three_particles},
         {"three_particle_binding_angle_resolution",
          collision_params.three_particle_angle_resolution},

         {"distance", collision_params.distance},
         {"distance_glued_particle_to_vs",
          collision_params.dist_glued_part_to_vs},
         {"vs_placement", collision_params.vs_placement},

         {"part_type_vs", collision_params.vs_particle_type},
         {"part_type_to_be_glued", collision_params.part_type_to_be_glued},
         {"part_type_to_attach_vs_to",
          collision_params.part_type_to_attach_vs_to},
         {"part_type_after_glueing",
          collision_params.part_type_after_glueing}});
  }

  Variant do_call_method(const std::string &name,
                         const VariantMap &params) override {
    if (name == "instantiate") {
      context()->parallel_try_catch([this, &params]() {
        auto collision_params_backup = ::collision_params;
        try {
          // check provided parameters
          check_input_parameters(params);
          // set parameters
          ::collision_params = Collision_parameters();
          for (auto const &kv : params) {
            do_set_parameter(get_value<std::string>(kv.first), kv.second);
          }
          // sanitize parameters and calculate derived parameters
          ::collision_params.initialize();
          return none;
        } catch (...) {
          // restore original parameters and re-throw exception
          ::collision_params = collision_params_backup;
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
    return none;
  }

private:
  void check_mode_name(std::string const &name) const {
    if (cd_name_to_mode.count(name) == 0) {
      throw std::invalid_argument("Unknown collision mode '" + name + "'");
    }
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
      auto const param_name = get_value<std::string>(kv.first);
      if (expected_param_names.count(param_name) == 0) {
        // throw if parameter is unknown
        std::ignore = get_parameter(param_name);
        // throw if parameter is known but doesn't match the selected mode
        throw std::runtime_error("Parameter '" + param_name + "' is not " +
                                 "required for mode '" + name + "'");
      }
      input_parameter_names.insert(param_name);
    }
    for (auto const &param_name : expected_param_names) {
      if (input_parameter_names.count(param_name) == 0) {
        throw std::runtime_error("Parameter '" + param_name + "' is " +
                                 "required for mode '" + name + "'");
      }
    }
  }
};

} /* namespace CollisionDetection */
} /* namespace ScriptInterface */

#endif // COLLISION_DETECTION
#endif
