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
#ifndef ESPRESSO_SRC_CORE_PAIR_CRITERIA_BOND_CRITERION_HPP
#define ESPRESSO_SRC_CORE_PAIR_CRITERIA_BOND_CRITERION_HPP

#include "pair_criteria/PairCriterion.hpp"

#include "BondList.hpp"
#include "BoxGeometry.hpp"
#include "system/System.hpp"

namespace PairCriteria {
/** @brief True if a bond of given type exists between two particles. */
class BondCriterion : public PairCriterion {
public:
  bool decide(Particle const &p1, Particle const &p2) const override {

    auto const &box_geo = *System::get_system().box_geo;
    auto const d = box_geo.get_mi_vector(p1.pos(), p2.pos()).norm();

    return (pair_bond_exists_on(p1.bonds(), p2.id(), m_bond_type) ||
            pair_bond_exists_on(p2.bonds(), p1.id(), m_bond_type)) &&
           d <= m_cut_off;
  }
  int get_bond_type() { return m_bond_type; }
  void set_bond_type(int t) { m_bond_type = t; }

  double get_cut_off() { return m_cut_off; }
  void set_cut_off(double c) { m_cut_off = c; }

private:
  int m_bond_type;
  double m_cut_off;
};
} // namespace PairCriteria

#endif
