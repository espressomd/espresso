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

#include "Globals.hpp"

#include "core/event.hpp"
#include "core/grid.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <utils/Vector.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace System {

static Utils::Vector3b get_periodicity() {
  return {::box_geo.periodic(0), ::box_geo.periodic(1), ::box_geo.periodic(2)};
}

static void set_periodicity(Utils::Vector3b const &value) {
  for (int i = 0; i < 3; ++i) {
    ::box_geo.set_periodic(i, value[i]);
  }
  on_periodicity_change();
}

Globals::Globals() {
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
  });
}

} // namespace System
} // namespace ScriptInterface
