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
#ifndef ESPRESSO_SRC_CORE_PAIR_CRITERIA_DISTANCE_CRITERION_HPP
#define ESPRESSO_SRC_CORE_PAIR_CRITERIA_DISTANCE_CRITERION_HPP

#include "pair_criteria/PairCriterion.hpp"

#include "grid.hpp"

namespace PairCriteria {
/**
 * @brief True if two particles are closer than a cut off distance,
 * respecting minimum image convention.
 */
class DistanceCriterion : public PairCriterion {
public:
  bool decide(const Particle &p1, const Particle &p2) const override {
    return box_geo.get_mi_vector(p1.pos(), p2.pos()).norm() <= m_cut_off;
  }
  double get_cut_off() { return m_cut_off; }
  void set_cut_off(double c) { m_cut_off = c; }

private:
  double m_cut_off;
};

} // namespace PairCriteria

#endif
