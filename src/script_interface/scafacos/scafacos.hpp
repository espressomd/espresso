/*
 * Copyright (C) 2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_SCAFACOS_SCAFACOS_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_SCAFACOS_SCAFACOS_HPP

#include "config/config.hpp"

#if defined(SCAFACOS) or defined(SCAFACOS_DIPOLAR)

#include "script_interface/Variant.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace ScriptInterface {
namespace Scafacos {

/** @brief Fetch list of methods compiled in ScaFaCoS. */
std::vector<std::string> available_methods();

/** @brief Flatten a parameter map. */
std::string serialize_parameters(Variant const &pack);

/** @brief Convert flattened parameters to a map. */
std::unordered_map<std::string, Variant>
deserialize_parameters(std::string const &parameters);

} // namespace Scafacos
} // namespace ScriptInterface

#endif // SCAFACOS or SCAFACOS_DIPOLAR
#endif
