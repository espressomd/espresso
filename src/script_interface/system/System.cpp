/*
 * Copyright (C) 2013-2022 The ESPResSo project
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

#include "System.hpp"

#include "config/config.hpp"

#include "core/grid.hpp"
#include "core/object-in-fluid/oif_global_forces.hpp"
#include "core/particle_data.hpp"
#include "core/particle_node.hpp"
#include "core/rotate_system.hpp"

#include <utils/Vector.hpp>

#include <string>
#include <vector>

namespace ScriptInterface {
namespace System {

static bool system_created = false;

System::System() {
  add_parameters({
      {"max_oif_objects", ::max_oif_objects},
  });
}

Variant System::do_call_method(std::string const &name,
                               VariantMap const &parameters) {
  if (name == "is_system_created") {
    return system_created;
  }
  if (name == "lock_system_creation") {
    system_created = true;
    return {};
  }
  if (name == "rescale_boxl") {
    auto const coord = get_value<int>(parameters, "coord");
    auto const length = get_value<double>(parameters, "length");
    auto const scale = (coord == 3) ? length * ::box_geo.length_inv()[0]
                                    : length * ::box_geo.length_inv()[coord];
    context()->parallel_try_catch([&]() {
      if (length <= 0.) {
        throw std::domain_error("Parameter 'd_new' be > 0");
      }
    });
    auto new_value = Utils::Vector3d{};
    if (coord == 3) {
      new_value = Utils::Vector3d::broadcast(length);
    } else {
      new_value = ::box_geo.length();
      new_value[coord] = length;
    }
    // when shrinking, rescale the particles first
    if (scale <= 1.) {
      rescale_particles(coord, scale);
    }
    set_box_length(new_value);
    if (scale > 1.) {
      rescale_particles(coord, scale);
    }
    return {};
  }
#ifdef EXCLUSIONS
  if (name == "auto_exclusions") {
    if (context()->is_head_node()) {
      auto const distance = get_value<int>(parameters, "distance");
      auto_exclusions(distance);
    }
    return {};
  }
#endif
  if (name == "setup_type_map") {
    if (context()->is_head_node()) {
      auto const types = get_value<std::vector<int>>(parameters, "type_list");
      for (auto const type : types) {
        init_type_map(type);
      }
    }
    return {};
  }
  if (name == "number_of_particles") {
    if (context()->is_head_node()) {
      auto const type = get_value<int>(parameters, "type");
      return number_of_particles_with_type(type);
    }
    return {};
  }
  if (name == "velocity_difference") {
    auto const pos1 = get_value<Utils::Vector3d>(parameters, "pos1");
    auto const pos2 = get_value<Utils::Vector3d>(parameters, "pos2");
    auto const v1 = get_value<Utils::Vector3d>(parameters, "v1");
    auto const v2 = get_value<Utils::Vector3d>(parameters, "v2");
    return ::box_geo.velocity_difference(pos2, pos1, v2, v1);
  }
  if (name == "distance_vec") {
    auto const pos1 = get_value<Utils::Vector3d>(parameters, "pos1");
    auto const pos2 = get_value<Utils::Vector3d>(parameters, "pos2");
    return ::box_geo.get_mi_vector(pos2, pos1);
  }
  if (name == "rotate_system") {
    rotate_system(get_value<double>(parameters, "phi"),
                  get_value<double>(parameters, "theta"),
                  get_value<double>(parameters, "alpha"));
    return {};
  }
  return {};
}

} // namespace System
} // namespace ScriptInterface
