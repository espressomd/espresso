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
#ifndef ESPRESSO_SRC_CORE_PAIR_CRITERIA_BOND_CRITERION_HPP
#define ESPRESSO_SRC_CORE_PAIR_CRITERIA_BOND_CRITERION_HPP

#include "pair_criteria/PairCriterion.hpp"

#include "bonds.hpp"

namespace PairCriteria {
/**
 * @brief True if a bond of given type exists between two particles.
 * Head node only: the bonds are fetched from the ranks of the particles.
 */
class BondCriterion : public PairCriterion {
public:
  bool decide(Particle const &p1, Particle const &p2) const override {
    auto const connects = [this](int owner_id, int partner_id) {
      for (auto const &[bond_id, ids] : get_owned_bonds(owner_id)) {
        if (bond_id == m_bond_type and ids.size() == 2u and
            ids[1] == partner_id) {
          return true;
        }
      }
      return false;
    };
    return p1.id() != p2.id() and
           (connects(p1.id(), p2.id()) or connects(p2.id(), p1.id()));
  }
  int get_bond_type() { return m_bond_type; }
  void set_bond_type(int t) { m_bond_type = t; }

private:
  int m_bond_type;
};
} // namespace PairCriteria

#endif
