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
#ifndef ESPRESSO_SRC_CORE_PAIR_CRITERIA_ENERGY_CRITERION_HPP
#define ESPRESSO_SRC_CORE_PAIR_CRITERIA_ENERGY_CRITERION_HPP

#include "pair_criteria/PairCriterion.hpp"

#include "energy_inline.hpp"

namespace PairCriteria {
/**
 * @brief True if the short-range energy is larger than a cutoff value.
 */
class EnergyCriterion : public PairCriterion {
public:
  bool decide(Particle const &p1, Particle const &p2) const override {
    // Distance between particles
    auto const d = box_geo.get_mi_vector(p1.pos(), p2.pos());

    // Interaction parameters for particle types
    auto const &ia_params = *get_ia_param(p1.type(), p2.type());
    auto const coulomb_kernel = Coulomb::pair_energy_kernel();

    auto const energy = calc_non_bonded_pair_energy(
        p1, p2, ia_params, d, d.norm(), coulomb_kernel.get_ptr());

    return energy >= m_cut_off;
  }
  double get_cut_off() { return m_cut_off; }
  void set_cut_off(double c) { m_cut_off = c; }

private:
  double m_cut_off;
};
} // namespace PairCriteria

#endif
