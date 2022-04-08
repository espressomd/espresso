/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "config.hpp"

#ifdef COLLISION_DETECTION

#include "script_interface/ScriptInterface.hpp"

#include "core/collision.hpp"

#include <cassert>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace CollisionDetection {

class CollisionDetection : public AutoParameters<CollisionDetection> {
public:
  CollisionDetection() {
    add_parameters(
        {{"mode",
          [](Variant const &v) {
            auto const mode = get_value<std::string>(v);
            if (mode == "off") {
              collision_params.mode = COLLISION_MODE_OFF;
            } else if (mode == "bind_centers") {
              collision_params.mode = COLLISION_MODE_BOND;
            } else if (mode == "bind_at_point_of_collision") {
              collision_params.mode = COLLISION_MODE_VS;
            } else if (mode == "glue_to_surface") {
              collision_params.mode = COLLISION_MODE_GLUE_TO_SURF;
            } else if (mode == "bind_three_particles") {
              collision_params.mode = COLLISION_MODE_BIND_THREE_PARTICLES;
            } else {
              throw std::invalid_argument("Unknown collision mode '" + mode +
                                          "'");
            }
          },
          []() {
            std::string mode = "";
            if (collision_params.mode == COLLISION_MODE_OFF) {
              mode = "off";
            } else if (collision_params.mode == COLLISION_MODE_BOND) {
              mode = "bind_centers";
            } else if (collision_params.mode == COLLISION_MODE_VS) {
              mode = "bind_at_point_of_collision";
            } else if (collision_params.mode == COLLISION_MODE_GLUE_TO_SURF) {
              mode = "glue_to_surface";
            } else if (collision_params.mode ==
                       COLLISION_MODE_BIND_THREE_PARTICLES) {
              mode = "bind_three_particles";
            }
            assert(not mode.empty());
            return Variant{mode};
          }},

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
      auto collision_params_backup = ::collision_params;
      try {
        // check provided parameters
        check_input_parameters(params);
        // set parameters
        ::collision_params = Collision_parameters();
        for (auto const &kv : params) {
          set_parameter(get_value<std::string>(kv.first), kv.second);
        }
        // sanitize parameters and calculate derived parameters
        ::collision_params.initialize();
        return none;
      } catch (...) {
        // restore original parameters and re-throw exception
        ::collision_params = collision_params_backup;
        if (context()->is_head_node()) {
          throw;
        }
        throw Exception("");
      }
    }
    if (name == "params_for_mode") {
      auto const mode = get_value<std::string>(params, "mode");
      return make_vector_of_variants(get_param_names(mode));
    }
    return none;
  }

private:
  std::vector<std::string> get_param_names(std::string const &mode) const {
    if (mode == "off") {
      return {"mode"};
    }
    if (mode == "bind_centers") {
      return {"mode", "bond_centers", "distance"};
    }
    if (mode == "bind_at_point_of_collision") {
      return {"mode",         "bond_centers", "bond_vs",
              "part_type_vs", "distance",     "vs_placement"};
    }
    if (mode == "glue_to_surface") {
      return {"mode",
              "bond_centers",
              "bond_vs",
              "part_type_vs",
              "part_type_to_be_glued",
              "part_type_to_attach_vs_to",
              "part_type_after_glueing",
              "distance",
              "distance_glued_particle_to_vs"};
    }
    if (mode == "bind_three_particles") {
      return {"mode", "bond_centers", "distance", "bond_three_particles",
              "three_particle_binding_angle_resolution"};
    }
    throw std::invalid_argument("Unknown collision mode '" + mode + "'");
  }

  void check_input_parameters(VariantMap const &params) const {
    auto const mode = get_value<std::string>(params, "mode");
    auto const expected_params = get_param_names(mode);
    auto const expected_param_names =
        std::set<std::string>{expected_params.begin(), expected_params.end()};
    std::set<std::string> input_parameter_names = {};
    for (auto const &kv : params) {
      auto const param_name = get_value<std::string>(kv.first);
      if (expected_param_names.count(param_name) == 0) {
        // throw if parameter is unknown
        static_cast<void>(get_parameter(param_name));
        // throw if parameter is known but doesn't match the selected mode
        throw std::runtime_error("Parameter '" + param_name + "' is not " +
                                 "required for mode '" + mode + "'");
      }
      input_parameter_names.insert(param_name);
    }
    for (auto const &param_name : expected_param_names) {
      if (input_parameter_names.count(param_name) == 0) {
        throw std::runtime_error("Parameter '" + param_name + "' is " +
                                 "required for mode '" + mode + "'");
      }
    }
  }
};

} /* namespace CollisionDetection */
} /* namespace ScriptInterface */

#endif // COLLISION_DETECTION
#endif
