/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include <config/config.hpp>

#include "Particle.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"

#include <utils/index.hpp>
#include <utils/math/sqr.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

struct GetNonbondedCutoff {
  GetNonbondedCutoff(System::System const &system) : m_system{system} {}
  auto operator()(int type_i, int type_j) const {
    return m_system.nonbonded_ias->get_ia_param(type_i, type_j).max_cut;
  }

private:
  System::System const &m_system;
};

/** Returns true if the particles are to be considered for short range
 *  interactions.
 *
 *  @tparam ShortRangeOnly  When true, the electrostatics / dipolar / collision
 *    early-accept branches are compiled out of the call operator. Selected by
 *    `update_verlet_state` only when none of those cutoffs are active, so the
 *    removed branches would never have accepted a pair -- the result is
 *    unchanged, but the per-candidate build loop drops the dead comparisons.
 */
template <typename CutoffGetter = GetNonbondedCutoff,
          bool ShortRangeOnly = false>
class VerletCriterion {
  const double m_skin;
  const double m_eff_max_cut2;
  const double m_eff_coulomb_cut2 = 0.;
  const double m_eff_dipolar_cut2 = 0.;
  const double m_collision_cut2 = 0.;
  double eff_cutoff_sqr(double x) const {
    if (x == inactive_cutoff)
      return inactive_cutoff;
    return Utils::sqr(x + m_skin);
  }
  /** Dense row-major table of squared effective (cutoff + skin) values per
   *  type pair, @ref inactive_cutoff for inactive pairs. The per-type-pair
   *  cutoff query runs once per candidate pair in the Verlet-list build, so
   *  it must be a plain load instead of a walk through the
   *  @ref InteractionsNonBonded pointer table.
   */
  std::vector<double> m_eff_cut2_table;
  int m_n_types;

public:
  VerletCriterion(System::System const &system, double skin, double max_cut,
                  double coulomb_cut = 0., double dipolar_cut = 0.,
                  double collision_detection_cutoff = 0.)
      : m_skin(skin), m_eff_max_cut2(eff_cutoff_sqr(max_cut)),
        m_eff_coulomb_cut2(eff_cutoff_sqr(coulomb_cut)),
        m_eff_dipolar_cut2(eff_cutoff_sqr(dipolar_cut)),
        m_collision_cut2(eff_cutoff_sqr(collision_detection_cutoff)) {
    CutoffGetter const get_nonbonded_cutoff(system);
    auto const max_type = system.nonbonded_ias->get_max_seen_particle_type();
    m_n_types = std::max(max_type + 1, 1);
    m_eff_cut2_table.assign(static_cast<std::size_t>(m_n_types) *
                                static_cast<std::size_t>(m_n_types),
                            inactive_cutoff);
    for (int type_i = 0; type_i <= max_type; ++type_i) {
      for (int type_j = type_i; type_j <= max_type; ++type_j) {
        auto const eff_cut2 =
            eff_cutoff_sqr(get_nonbonded_cutoff(type_i, type_j));
        m_eff_cut2_table[static_cast<std::size_t>(type_i) *
                             static_cast<std::size_t>(m_n_types) +
                         static_cast<std::size_t>(type_j)] = eff_cut2;
        m_eff_cut2_table[static_cast<std::size_t>(type_j) *
                             static_cast<std::size_t>(m_n_types) +
                         static_cast<std::size_t>(type_i)] = eff_cut2;
      }
    }
  }

  bool operator()(const Particle &p1, const Particle &p2, double dist2) const {
    if (dist2 > m_eff_max_cut2)
      return false;

    if constexpr (not ShortRangeOnly) {
#ifdef ESPRESSO_ELECTROSTATICS
      // Within real space cutoff of electrostatics and both are charged
      if (dist2 <= m_eff_coulomb_cut2 and p1.q() != 0. and p2.q() != 0.)
        return true;
#endif

#ifdef ESPRESSO_DIPOLES
      // Within dipolar cutoff and both carry magnetic moments
      if (dist2 <= m_eff_dipolar_cut2 and p1.dipm() != 0. and p2.dipm() != 0.)
        return true;
#endif

#ifdef ESPRESSO_COLLISION_DETECTION
      // Collision detection
      if (dist2 <= m_collision_cut2)
        return true;
#endif
    }

    // Within short-range distance (including dpd and the like). Inactive
    // pairs hold inactive_cutoff (negative) in the table, so the comparison
    // rejects them without a separate activity check.
    auto const type_i = p1.type();
    auto const type_j = p2.type();
    assert(type_i >= 0 and type_i < m_n_types);
    assert(type_j >= 0 and type_j < m_n_types);
    return dist2 <= m_eff_cut2_table[static_cast<std::size_t>(type_i) *
                                         static_cast<std::size_t>(m_n_types) +
                                     static_cast<std::size_t>(type_j)];
  }
};
