/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef _BONDED_INTERACTION_UTILS_HPP
#define _BONDED_INTERACTION_UTILS_HPP

#include "bonded_interaction_data.hpp"

#include "BondList.hpp"
#include "Particle.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/variant.hpp>

/** @brief Checks both particles for a specific bond, even on ghost particles.
 *
 *  @param p           particle to check for the bond
 *  @param p_partner   possible bond partner
 *  @tparam BondType   Bond type to check for. Must be of one of the types in
 *                     @ref Bonded_IA_Parameters.
 */
template <typename BondType>
inline bool pair_bond_enum_exists_on(Particle const &p,
                                     Particle const &p_partner) {
  return boost::algorithm::any_of(
      p.bonds(), [partner_id = p_partner.identity()](BondView const &bond) {
        auto const &bond_ptr = bonded_ia_params.at(bond.bond_id());
        return (boost::get<BondType>(bond_ptr.get()) != nullptr) and
               (bond.partner_ids()[0] == partner_id);
      });
}

/** @brief Checks both particles for a specific bond, even on ghost particles.
 *
 *  @param p1     particle on which the bond may be stored
 *  @param p2     particle on which the bond may be stored
 *  @tparam BondType Bond type to check for. Must be of one of the types in
 *                @ref Bonded_IA_Parameters.
 */
template <typename BondType>
inline bool pair_bond_enum_exists_between(Particle const &p1,
                                          Particle const &p2) {
  if (&p1 == &p2)
    return false;

  // Check if particles have bonds and search for the bond of interest.
  // Could be saved on both sides (and both could have other bonds), so
  // we need to check both.
  return pair_bond_enum_exists_on<BondType>(p1, p2) or
         pair_bond_enum_exists_on<BondType>(p2, p1);
}

#endif
