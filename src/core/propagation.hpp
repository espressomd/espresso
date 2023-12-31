/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "PropagationMode.hpp"

#include <string>
#include <unordered_map>

/**
 * @brief Check allowlist of valid propagation modes combinations.
 */
bool is_valid_propagation_combination(int propagation);

/**
 * @brief Convert @ref PropagationMode::PropagationMode to name/value pairs.
 */
std::unordered_map<std::string, int> propagation_flags_map();

/**
 * @brief Convert a propagation modes bitmask to a string.
 */
std::string propagation_bitmask_to_string(int propagation);
