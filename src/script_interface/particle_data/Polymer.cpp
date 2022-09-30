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

#include "Polymer.hpp"

#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

#include "core/partCfg_global.hpp"
#include "core/polymer.hpp"

#include <utils/Vector.hpp>

#include <string>
#include <vector>

namespace ScriptInterface {
namespace Particles {

Variant Polymer::do_call_method(std::string const &name,
                                VariantMap const &parameters) {
  if (name == "draw_polymer_positions") {
    auto const positions = draw_polymer_positions(
        partCfg(), get_value<int>(parameters, "n_polymers"),
        get_value<int>(parameters, "beads_per_chain"),
        get_value<double>(parameters, "bond_length"),
        get_value<std::vector<Utils::Vector3d>>(parameters, "start_positions"),
        get_value<double>(parameters, "min_distance"),
        get_value<int>(parameters, "max_tries"),
        get_value<bool>(parameters, "use_bond_angle"),
        get_value_or<double>(parameters, "bond_angle", 0.),
        get_value<bool>(parameters, "respect_constraints"),
        get_value<int>(parameters, "seed"));
    std::vector<Variant> pack;
    for (auto const &chain : positions) {
      pack.emplace_back(make_vector_of_variants(chain));
    }
    return pack;
  }
  return {};
}

} // namespace Particles
} // namespace ScriptInterface
