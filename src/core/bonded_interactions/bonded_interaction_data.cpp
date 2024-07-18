/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#include "rigid_bond.hpp"
#include "system/System.hpp"
#include "thermalized_bond.hpp"
#include "thermostat.hpp"

#include <boost/variant.hpp>

#include <algorithm>
#include <numeric>

/** Visitor to get the bond cutoff from the bond parameter variant */
class BondCutoff : public boost::static_visitor<double> {
public:
  template <typename T> double operator()(T const &bond) const {
    return bond.cutoff();
  }
};

double BondedInteractionsMap::maximal_cutoff() const {
  auto const max_cut_bonded = std::accumulate(
      begin(), end(), BONDED_INACTIVE_CUTOFF, [](auto max_cut, auto const &kv) {
        return std::max(max_cut,
                        boost::apply_visitor(BondCutoff(), *kv.second));
      });

  /* Check if there are dihedrals */
  auto const any_dihedrals = std::any_of(begin(), end(), [](auto const &kv) {
    return (boost::get<DihedralBond>(&(*kv.second)) ||
            boost::get<TabulatedDihedralBond>(&(*kv.second)));
  });

  /* dihedrals: the central particle is indirectly connected to the fourth
   * particle via the third particle, so we have to double the cutoff */
  return (any_dihedrals) ? 2 * max_cut_bonded : max_cut_bonded;
}

void BondedInteractionsMap::on_ia_change() {
  n_thermalized_bonds = 0;
#ifdef BOND_CONSTRAINT
  n_rigid_bonds = 0;
#endif
  for (auto &kv : *this) {
    if (boost::get<ThermalizedBond>(&(*kv.second)) != nullptr) {
      ++n_thermalized_bonds;
    }
#ifdef BOND_CONSTRAINT
    if (boost::get<RigidBond>(&(*kv.second)) != nullptr) {
      ++n_rigid_bonds;
    }
#endif
  }
  if (auto system = m_system.lock()) {
    system->on_short_range_ia_change();
    system->on_thermostat_param_change(); // thermalized bonds
  }
}

void BondedInteractionsMap::activate_bond(mapped_type const &ptr) {
  auto &system = get_system();
  if (auto bond = boost::get<ThermalizedBond>(ptr.get())) {
    bond->set_thermostat_view(system.thermostat);
  }
}

void BondedInteractionsMap::deactivate_bond(mapped_type const &ptr) {
  if (auto bond = boost::get<ThermalizedBond>(ptr.get())) {
    bond->unset_thermostat_view();
  }
}
