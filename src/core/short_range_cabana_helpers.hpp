/*
 * Copyright (C) 2026 The ESPResSo project
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

#include "PropagationMode.hpp"
#include "aosoa_pack.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <Kokkos_Core.hpp>

// Helper functions to check if specific algorithms are active
#ifdef ESPRESSO_GAY_BERNE
// Configured-for-type-pair: mask bit set iff gay_berne has a real cutoff.
// In-range check is deferred to where it matters.
ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool
gay_berne_configured(IA_parameters const &ia_params) {
  return (ia_params.active_pair_mask &
          pair_potential_bit(PairPotential::GayBerne)) != 0u;
}
ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool
gay_berne_active(double dist, IA_parameters const &ia_params) {
  return gay_berne_configured(ia_params) and dist < ia_params.gay_berne.cut;
}
#endif

#ifdef ESPRESSO_DPD
// True iff the DPD thermostat is on AND this type pair has DPD configured.
ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool
dpd_active(IA_parameters const &ia_params, int thermo_switch) {
  return (thermo_switch & THERMO_DPD) and
         ((ia_params.active_pair_mask &
           pair_potential_bit(PairPotential::DPD)) != 0u);
}
#endif

#ifdef ESPRESSO_THOLE
KOKKOS_INLINE_FUNCTION bool thole_active(IA_parameters const &ia_params,
                                         bool has_coulomb_kernel) {
  return (ia_params.thole.scaling_coeff != 0. and ia_params.thole.q1q2 != 0. and
          has_coulomb_kernel);
}
#endif
