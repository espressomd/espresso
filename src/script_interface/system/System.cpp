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

#include "core/cells.hpp"
#include "core/event.hpp"
#include "core/grid.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/object-in-fluid/oif_global_forces.hpp"
#include "core/particle_node.hpp"
#include "core/rotate_system.hpp"
#include "core/system/System.hpp"

#include <utils/Vector.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace System {

static bool system_created = false;

static Utils::Vector3b get_periodicity() {
  return {::box_geo.periodic(0), ::box_geo.periodic(1), ::box_geo.periodic(2)};
}

static void set_periodicity(Utils::Vector3b const &value) {
  for (int i = 0; i < 3; ++i) {
    ::box_geo.set_periodic(i, value[i]);
  }
  on_periodicity_change();
}

System::System() {
  add_parameters({
      {"box_l",
       [this](Variant const &v) {
         context()->parallel_try_catch([&]() {
           auto const new_value = get_value<Utils::Vector3d>(v);
           if (not(new_value > Utils::Vector3d::broadcast(0.))) {
             throw std::domain_error("Attribute 'box_l' must be > 0");
           }
           box_geo.set_length(new_value);
           on_boxl_change();
         });
       },
       []() { return ::box_geo.length(); }},
      {"min_global_cut",
       [this](Variant const &v) {
         context()->parallel_try_catch([&]() {
           auto const new_value = get_value<double>(v);
           if (new_value < 0. and new_value != INACTIVE_CUTOFF) {
             throw std::domain_error("Attribute 'min_global_cut' must be >= 0");
           }
           set_min_global_cut(new_value);
         });
       },
       []() { return ::get_min_global_cut(); }},
      {"periodicity",
       [this](Variant const &v) {
         context()->parallel_try_catch(
             [&]() { set_periodicity(get_value<Utils::Vector3b>(v)); });
       },
       []() { return get_periodicity(); }},
      {"max_oif_objects", ::max_oif_objects},
  });
}

void System::do_construct(VariantMap const &params) {
  /* When reloading the system state from a checkpoint file,
   * the order of global variables instantiation matters.
   * The @c box_l must be set before any other global variable.
   * All these globals re-initialize the cell system, and we
   * cannot use the default-constructed @c box_geo when e.g.
   * long-range interactions exist in the system, otherwise
   * runtime errors about the local geometry being smaller
   * than the interaction range would be raised.
   */
  auto const valid_params = get_valid_parameters();
  std::set<std::string> const skip{{"box_l", "periodicity", "min_global_cut"}};
  std::set<std::string> const allowed{valid_params.begin(), valid_params.end()};

  m_instance = std::make_shared<::System::System>();
  ::System::set_system(m_instance);
  m_instance->init();

  do_set_parameter("box_l", params.at("box_l"));
  if (params.count("periodicity")) {
    do_set_parameter("periodicity", params.at("periodicity"));
  }
  if (params.count("min_global_cut")) {
    do_set_parameter("min_global_cut", params.at("min_global_cut"));
  }
  for (auto const &[name, value] : params) {
    if (skip.count(name) == 0 and allowed.count(name) != 0) {
      do_set_parameter(name, value);
    }
  }
}

/** Rescale all particle positions in direction @p dir by a factor @p scale.
 *  @param dir   direction to scale (0/1/2 = x/y/z, 3 = x+y+z isotropically)
 *  @param scale factor by which to rescale (>1: stretch, <1: contract)
 */
static void rescale_particles(int dir, double scale) {
  for (auto &p : ::cell_structure.local_particles()) {
    if (dir < 3)
      p.pos()[dir] *= scale;
    else {
      p.pos() *= scale;
    }
  }
  on_particle_change();
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
  if (name == "setup_type_map") {
    auto const types = get_value<std::vector<int>>(parameters, "type_list");
    for (auto const type : types) {
      ::init_type_map(type);
    }
    return {};
  }
  if (name == "number_of_particles") {
    auto const type = get_value<int>(parameters, "type");
    return ::number_of_particles_with_type(type);
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
  if (name == "session_shutdown") {
    if (m_instance) {
      if (&::System::get_system() == m_instance.get()) {
        ::System::reset_system();
      }
      assert(m_instance.use_count() == 1l);
      m_instance.reset();
    }
    return {};
  }
  return {};
}

} // namespace System
} // namespace ScriptInterface
