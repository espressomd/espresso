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
/** \file
 *  Implementation of nonbonded_interaction_data.hpp
 */

#include <config/config.hpp>

#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "electrostatics/coulomb.hpp"
#include "system/System.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

// Set mask bit and fold cutoff into max_cut if the potential's max_cutoff()
// exceeds the current maximum. The "active" criterion matches each kernel's
// own guard: a kernel can only contribute when its max_cutoff() > 0 (since
// dist >= 0 and inactive_cutoff = -1).
[[maybe_unused]] static void consider(double &max_cut_current, unsigned &mask,
                                      PairPotential p, double sub_cutoff) {
  if (sub_cutoff > 0.) {
    mask |= pair_potential_bit(p);
    if (sub_cutoff > max_cut_current) {
      max_cut_current = sub_cutoff;
    }
  }
}

static std::pair<double, unsigned>
recalc_maximal_cutoff(IA_parameters const &data,
                      [[maybe_unused]] System::System const &system) {
  auto max_cut_current = inactive_cutoff;
  auto mask = 0u;

#ifdef ESPRESSO_LENNARD_JONES
  consider(max_cut_current, mask, PairPotential::LennardJones,
           data.lj.max_cutoff());
#endif

#ifdef ESPRESSO_WCA
  consider(max_cut_current, mask, PairPotential::WCA, data.wca.max_cutoff());
#endif

#ifdef ESPRESSO_DPD
  consider(max_cut_current, mask, PairPotential::DPD, data.dpd.max_cutoff());
#endif

#ifdef ESPRESSO_LENNARD_JONES_GENERIC
  consider(max_cut_current, mask, PairPotential::LennardJonesGeneric,
           data.ljgen.max_cutoff());
#endif

#ifdef ESPRESSO_SMOOTH_STEP
  consider(max_cut_current, mask, PairPotential::SmoothStep,
           data.smooth_step.max_cutoff());
#endif

#ifdef ESPRESSO_HERTZIAN
  consider(max_cut_current, mask, PairPotential::Hertzian,
           data.hertzian.max_cutoff());
#endif

#ifdef ESPRESSO_GAUSSIAN
  consider(max_cut_current, mask, PairPotential::Gaussian,
           data.gaussian.max_cutoff());
#endif

#ifdef ESPRESSO_BMHTF_NACL
  consider(max_cut_current, mask, PairPotential::BMHTF,
           data.bmhtf.max_cutoff());
#endif

#ifdef ESPRESSO_MORSE
  consider(max_cut_current, mask, PairPotential::Morse,
           data.morse.max_cutoff());
#endif

#ifdef ESPRESSO_BUCKINGHAM
  consider(max_cut_current, mask, PairPotential::Buckingham,
           data.buckingham.max_cutoff());
#endif

#ifdef ESPRESSO_SOFT_SPHERE
  consider(max_cut_current, mask, PairPotential::SoftSphere,
           data.soft_sphere.max_cutoff());
#endif

#ifdef ESPRESSO_HAT
  consider(max_cut_current, mask, PairPotential::Hat, data.hat.max_cutoff());
#endif

#ifdef ESPRESSO_LJCOS
  consider(max_cut_current, mask, PairPotential::LJCos,
           data.ljcos.max_cutoff());
#endif

#ifdef ESPRESSO_LJCOS2
  consider(max_cut_current, mask, PairPotential::LJCos2,
           data.ljcos2.max_cutoff());
#endif

#ifdef ESPRESSO_GAY_BERNE
  consider(max_cut_current, mask, PairPotential::GayBerne,
           data.gay_berne.max_cutoff());
#endif

#ifdef ESPRESSO_TABULATED
  consider(max_cut_current, mask, PairPotential::Tabulated, data.tab.cutoff());
#endif

#ifdef ESPRESSO_THOLE
  // If THOLE is active, use p3m cutoff
  if (data.thole.scaling_coeff != 0.)
    max_cut_current = std::max(max_cut_current, system.coulomb.cutoff());
#endif

  return {max_cut_current, mask};
}

void InteractionsNonBonded::recalc_maximal_cutoffs() {
  auto const &system = get_system();
  auto combined_mask = 0u;
#ifdef ESPRESSO_THOLE
  auto any_thole_configured = false;
#endif
  for (auto &data : m_nonbonded_ia_params) {
    auto const [mc, mask] = recalc_maximal_cutoff(*data, system);
    data->max_cut = mc;
    data->active_pair_mask = mask;
    combined_mask |= mask;
#ifdef ESPRESSO_THOLE
    any_thole_configured |=
        (data->thole.scaling_coeff != 0. and data->thole.q1q2 != 0.);
#endif
  }
  m_combined_active_pair_mask = combined_mask;
#ifdef ESPRESSO_THOLE
  m_any_thole_configured = any_thole_configured;
#endif
}

double InteractionsNonBonded::maximal_cutoff() const {
  auto max_cut_nonbonded = inactive_cutoff;
  for (auto &data : m_nonbonded_ia_params) {
    max_cut_nonbonded = std::max(max_cut_nonbonded, data->max_cut);
  }
  return max_cut_nonbonded;
}

void InteractionsNonBonded::on_non_bonded_ia_change() const {
  get_system().on_non_bonded_ia_change();
}
