/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef CORE_NB_IA_VERLETCRITERION_HPP
#define CORE_NB_IA_VERLETCRITERION_HPP

#include "Particle.hpp"
#include "config.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <utils/index.hpp>
#include <utils/math/sqr.hpp>

struct GetNonbondedCutoff {
  auto operator()(int type_i, int type_j) const {
    return get_ia_param(type_i, type_j)->max_cut;
  }
};

/** Returns true if the particles are to be considered for short range
 *  interactions.
 */
template <typename CutoffGetter = GetNonbondedCutoff> class VerletCriterion {
  const double m_skin;
  const double m_eff_max_cut2;
  const double m_eff_coulomb_cut2 = 0.;
  const double m_eff_dipolar_cut2 = 0.;
  const double m_collision_cut2 = 0.;
  double eff_cutoff_sqr(double x) const {
    if (x == INACTIVE_CUTOFF)
      return INACTIVE_CUTOFF;
    return Utils::sqr(x + m_skin);
  }
  CutoffGetter get_nonbonded_cutoff;

public:
  VerletCriterion(double skin, double max_cut, double coulomb_cut = 0.,
                  double dipolar_cut = 0.,
                  double collision_detection_cutoff = 0.)
      : m_skin(skin), m_eff_max_cut2(eff_cutoff_sqr(max_cut)),
        m_eff_coulomb_cut2(eff_cutoff_sqr(coulomb_cut)),
        m_eff_dipolar_cut2(eff_cutoff_sqr(dipolar_cut)),
        m_collision_cut2(eff_cutoff_sqr(collision_detection_cutoff)) {}

  template <typename Distance>
  bool operator()(const Particle &p1, const Particle &p2,
                  Distance const &dist) const {
    auto const &dist2 = dist.dist2;
    if (dist2 > m_eff_max_cut2)
      return false;

#ifdef ELECTROSTATICS
    // Within real space cutoff of electrostatics and both charged
    if ((dist2 <= m_eff_coulomb_cut2) && (p1.p.q != 0) && (p2.p.q != 0))
      return true;
#endif

#ifdef DIPOLES
    // Within dipolar cutoff and both carry magnetic moments
    if ((dist2 <= m_eff_dipolar_cut2) && (p1.p.dipm != 0) && (p2.p.dipm != 0))
      return true;
#endif

#ifdef COLLISION_DETECTION
    // Collision detection
    if (dist2 <= m_collision_cut2)
      return true;
#endif

    // Within short-range distance (including dpd and the like)
    auto const ia_cut = get_nonbonded_cutoff(p1.p.type, p2.p.type);
    return (ia_cut != INACTIVE_CUTOFF) &&
           (dist2 <= Utils::sqr(ia_cut + m_skin));
  }
};
#endif
