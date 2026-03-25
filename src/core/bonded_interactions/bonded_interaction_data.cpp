/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include "immersed_boundary/ImmersedBoundaries.hpp"
#include "rigid_bond.hpp"
#include "system/System.hpp"
#include "thermalized_bond.hpp"

#include <algorithm>
#include <numeric>
#include <ranges>
#include <variant>

double BondedInteractionsMap::maximal_cutoff() const {
  auto const max_cut_bonded = std::accumulate(
      begin(), end(), bonded_inactive_cutoff, [](auto max_cut, auto const &kv) {
        auto constexpr visitor = [](auto const &bond) { return bond.cutoff(); };
        return std::max(max_cut, std::visit(visitor, *kv.second));
      });

  /* Check if there are dihedrals */
  auto const any_dihedrals = std::ranges::any_of(*this, [](auto const &kv) {
    return (std::holds_alternative<DihedralBond>(*kv.second) or
            std::holds_alternative<TabulatedDihedralBond>(*kv.second));
  });

  /* dihedrals: the central particle is indirectly connected to the fourth
   * particle via the third particle, so we have to double the cutoff */
  return (any_dihedrals) ? 2 * max_cut_bonded : max_cut_bonded;
}

void BondedInteractionsMap::on_ia_change() {
  n_thermalized_bonds = 0;
#ifdef ESPRESSO_BOND_CONSTRAINT
  n_rigid_bonds = 0;
#endif
  for (auto const &bond : std::views::elements<1>(*this)) {
    if (std::holds_alternative<ThermalizedBond>(*bond)) {
      ++n_thermalized_bonds;
    }
#ifdef ESPRESSO_BOND_CONSTRAINT
    if (std::holds_alternative<RigidBond>(*bond)) {
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
  if (auto bond = std::get_if<ThermalizedBond>(ptr.get())) {
    bond->set_thermostat_view(system.thermostat);
  }
  if (auto bond = std::get_if<IBMVolCons>(ptr.get())) {
    system.immersed_boundaries->register_softID(*bond);
  }
  if (auto bond = std::get_if<IBMTriel>(ptr.get())) {
    bond->initialize(*system.box_geo, *system.cell_structure);
  }
  if (auto bond = std::get_if<IBMTribend>(ptr.get())) {
    bond->initialize(*system.box_geo, *system.cell_structure);
  }
}

void BondedInteractionsMap::deactivate_bond(mapped_type const &ptr) {
  if (auto bond = std::get_if<ThermalizedBond>(ptr.get())) {
    bond->unset_thermostat_view();
    n_thermalized_bonds = -1;
  }
  if (auto bond = std::get_if<IBMVolCons>(ptr.get())) {
    bond->unset_volumes_view();
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (std::get_if<RigidBond>(ptr.get())) {
    n_rigid_bonds = -1;
  }
#endif
}
