/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef _SCRIPT_INTERFACE_INTERACTIONS_BONDED_HPP
#define _SCRIPT_INTERFACE_INTERACTIONS_BONDED_HPP

/** @file
 *  Functions to interface with the core boost::variant.
 */

#include "core/bonded_interactions/bonded_interaction_data.hpp"

#include <boost/variant.hpp>

/** Return the 0-based type number of the specified bond. */
inline int bonded_ia_params_zero_based_type(int bond_id) {
  if (bonded_ia_params.contains(bond_id)) {
    return (*bonded_ia_params.at(bond_id)).which();
  }
  return 0;
}

/** Return the total number of bonds. */
inline int bonded_ia_params_size() { return bonded_ia_params.size(); }

/** Return the next key. */
inline int bonded_ia_params_next_key() {
  return bonded_ia_params.get_next_key();
}

#endif
