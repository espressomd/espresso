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

#include <stdexcept>

/** Return the specified bonded parameters of the given type
 *  if the type matches.
 */
template <typename BondType> BondType bonded_ia_params_at(int bond_id) {
  if (bond_id < 0 or bond_id >= bonded_ia_params.size()) {
    throw std::out_of_range("Access out of bounds");
  }
  return boost::get<BondType>(bonded_ia_params[bond_id]);
}

/** Check if the specified bond is of a specific type. */
template <typename BondType> bool bonded_ia_params_is_type(int bond_id) {
  if (bond_id < 0 or bond_id >= bonded_ia_params.size()) {
    throw std::out_of_range("Access out of bounds");
  }
  return boost::get<BondType>(&bonded_ia_params[bond_id]) != nullptr;
}

/** Return the number of bonded partners for the specified bond. */
inline int bonded_ia_params_num_partners(int bond_id) {
  if (bond_id < 0 or bond_id >= bonded_ia_params.size()) {
    throw std::out_of_range("Access out of bounds");
  }
  return number_of_partners(bonded_ia_params[bond_id]);
}

/** Return the 0-based type number of the specified bond. */
inline int bonded_ia_params_zero_based_type(int bond_id) {
  return bonded_ia_params[bond_id].which();
}

/** Return the total number of bonds. */
inline int bonded_ia_params_size() { return bonded_ia_params.size(); }

#endif
