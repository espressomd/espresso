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

#pragma once

#include "pair_criteria/PairCriterion.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "energy_inline.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"

namespace PairCriteria {
/**
 * @brief True if the short-range energy is larger than a cutoff value.
 */
class EnergyCriterion : public PairCriterion {
public:
  EnergyCriterion(System::System const &system) : m_system{system} {}
  bool decide(Particle const &p1, Particle const &p2) const override {
    // Distance between particles
    auto const d = m_system.box_geo->get_mi_vector(p1.pos(), p2.pos());

    // Interaction parameters for particle types
    auto const &ia_params =
        m_system.nonbonded_ias->get_ia_param(p1.type(), p2.type());
    auto const coulomb_kernel = m_system.coulomb.pair_energy_kernel();

    auto const energy = calc_non_bonded_pair_energy(
        p1, p2, ia_params, d, d.norm(), *m_system.bonded_ias,
        get_ptr(coulomb_kernel));

    return energy >= e_cut_off && d.norm() <= m_cut_off;
  }
  double get_cut_off() { return m_cut_off; }
  void set_cut_off(double c) { m_cut_off = c; }

  double get_e_cut_off() { return e_cut_off; }
  void set_e_cut_off(double e) { e_cut_off = e; }

private:
  double m_cut_off, e_cut_off;
  System::System const &m_system;
};
} // namespace PairCriteria
