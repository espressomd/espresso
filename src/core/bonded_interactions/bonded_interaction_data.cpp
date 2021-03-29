/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "bonded_interaction_data.hpp"
#include "interactions.hpp"

#include <boost/range/numeric.hpp>
#include <boost/variant.hpp>

#include <utils/constants.hpp>

#include <algorithm>
#include <cstddef>
#include <vector>

std::vector<Bonded_IA_Parameters> bonded_ia_params;

/** Visitor to get the bond cutoff from the bond parameter variant */
class BondCutoff : public boost::static_visitor<double> {
public:
  template <typename T> double operator()(T const &bond) const {
    return bond.cutoff();
  }
};

double maximal_cutoff_bonded() {
  auto const max_cut_bonded = boost::accumulate(
      bonded_ia_params, BONDED_INACTIVE_CUTOFF,
      [](auto max_cut, Bonded_IA_Parameters const &bond) {
        return std::max(max_cut, boost::apply_visitor(BondCutoff(), bond));
      });

  /* Check if there are dihedrals */
  auto const any_dihedrals = std::any_of(
      bonded_ia_params.begin(), bonded_ia_params.end(), [](auto const &bond) {
        return (boost::get<DihedralBond>(&bond) ||
                boost::get<TabulatedDihedralBond>(&bond));
      });

  /* dihedrals: the central particle is indirectly connected to the fourth
   * particle via the third particle, so we have to double the cutoff */
  return (any_dihedrals) ? 2 * max_cut_bonded : max_cut_bonded;
}

void make_bond_type_exist(int type) {
  std::size_t ns = static_cast<std::size_t>(type) + 1;
  auto const old_size = bonded_ia_params.size();
  if (ns <= old_size) {
    return;
  }
  /* else allocate new memory */
  bonded_ia_params.resize(ns, NoneBond());
}
