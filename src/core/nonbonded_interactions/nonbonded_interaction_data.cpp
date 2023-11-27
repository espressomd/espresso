/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
/** \file
 *  Implementation of nonbonded_interaction_data.hpp
 */
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "electrostatics/coulomb.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

static double recalc_maximal_cutoff(IA_parameters const &data) {
  auto max_cut_current = INACTIVE_CUTOFF;

#ifdef LENNARD_JONES
  max_cut_current = std::max(max_cut_current, data.lj.max_cutoff());
#endif

#ifdef WCA
  max_cut_current = std::max(max_cut_current, data.wca.max_cutoff());
#endif

#ifdef DPD
  max_cut_current = std::max(max_cut_current, data.dpd.max_cutoff());
#endif

#ifdef LENNARD_JONES_GENERIC
  max_cut_current = std::max(max_cut_current, data.ljgen.max_cutoff());
#endif

#ifdef SMOOTH_STEP
  max_cut_current = std::max(max_cut_current, data.smooth_step.max_cutoff());
#endif

#ifdef HERTZIAN
  max_cut_current = std::max(max_cut_current, data.hertzian.max_cutoff());
#endif

#ifdef GAUSSIAN
  max_cut_current = std::max(max_cut_current, data.gaussian.max_cutoff());
#endif

#ifdef BMHTF_NACL
  max_cut_current = std::max(max_cut_current, data.bmhtf.max_cutoff());
#endif

#ifdef MORSE
  max_cut_current = std::max(max_cut_current, data.morse.max_cutoff());
#endif

#ifdef BUCKINGHAM
  max_cut_current = std::max(max_cut_current, data.buckingham.max_cutoff());
#endif

#ifdef SOFT_SPHERE
  max_cut_current = std::max(max_cut_current, data.soft_sphere.max_cutoff());
#endif

#ifdef HAT
  max_cut_current = std::max(max_cut_current, data.hat.max_cutoff());
#endif

#ifdef LJCOS
  max_cut_current = std::max(max_cut_current, data.ljcos.max_cutoff());
#endif

#ifdef LJCOS2
  max_cut_current = std::max(max_cut_current, data.ljcos2.max_cutoff());
#endif

#ifdef GAY_BERNE
  max_cut_current = std::max(max_cut_current, data.gay_berne.max_cutoff());
#endif

#ifdef TABULATED
  max_cut_current = std::max(max_cut_current, data.tab.cutoff());
#endif

#ifdef THOLE
  // If THOLE is active, use p3m cutoff
  if (data.thole.scaling_coeff != 0.)
    max_cut_current =
        std::max(max_cut_current, Coulomb::get_coulomb().cutoff());
#endif

  return max_cut_current;
}

void InteractionsNonBonded::recalc_maximal_cutoffs() {
  for (auto &data : m_nonbonded_ia_params) {
    data->max_cut = recalc_maximal_cutoff(*data);
  }
}

double InteractionsNonBonded::maximal_cutoff() const {
  auto max_cut_nonbonded = INACTIVE_CUTOFF;
  for (auto &data : m_nonbonded_ia_params) {
    max_cut_nonbonded = std::max(max_cut_nonbonded, data->max_cut);
  }
  return max_cut_nonbonded;
}
