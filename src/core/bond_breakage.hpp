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

#pragma once

#include <memory>
#include <unordered_map>

namespace BondBreakage {

enum class ActionType {
  NONE = 0,
  DELETE_BOND = 1,
  REVERT_BIND_AT_POINT_OF_COLLISION = 2
};

struct BreakageSpec {
  double breakage_length;
  ActionType action_type;
};

extern std::unordered_map<int, std::shared_ptr<BreakageSpec>> breakage_specs;

/** @brief Check if the bond between the particles should break, if yes, queue
 *  it.
 */
bool check_and_handle_breakage(int particle_id, int bond_partner_id,
                               int bond_type, double distance);

void clear_queue();

void process_queue();

} // namespace BondBreakage
