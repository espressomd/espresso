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

#ifndef ESPRESSO_SRC_WALBERLA_BRIDGE_ELECTROKINETICS_REACTIONS_UTILS_HPP
#define ESPRESSO_SRC_WALBERLA_BRIDGE_ELECTROKINETICS_REACTIONS_UTILS_HPP

#include "EKReactant.hpp"
#include <domain_decomposition/BlockDataID.h>

namespace walberla::detail {
template <typename FloatType>
auto get_reaction_details(
    const std::shared_ptr<EKReactant<FloatType>> &reactant) {
  const auto order = reactant->get_order();
  const auto stoech_coeff = reactant->get_stoech_coeff();
  const auto density_id =
      walberla::BlockDataID(reactant->get_species()->get_density_id());

  return std::make_tuple(density_id, order, stoech_coeff);
}
} // namespace walberla::detail

#endif // ESPRESSO_SRC_WALBERLA_BRIDGE_ELECTROKINETICS_REACTIONS_UTILS_HPP
