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

#include "aosoa_pack.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <Kokkos_Core.hpp>

// Helper functions to check if specific algorithms are active
#ifdef ESPRESSO_GAY_BERNE
KOKKOS_INLINE_FUNCTION bool gay_berne_active(double dist,
                                             IA_parameters const &ia_params) {
  return dist < ia_params.gay_berne.cut;
}
#endif

#ifdef ESPRESSO_THOLE
KOKKOS_INLINE_FUNCTION bool thole_active(IA_parameters const &ia_params,
                                         bool has_coulomb_kernel) {
  return (ia_params.thole.scaling_coeff != 0. and ia_params.thole.q1q2 != 0. and
          has_coulomb_kernel);
}
#endif

struct PairDataFlags {
  bool need_directors = false;
  bool need_particle_pointers = false;
};

KOKKOS_INLINE_FUNCTION PairDataFlags compute_pair_data_flags(
    [[maybe_unused]] double dist,
    [[maybe_unused]] IA_parameters const &ia_params,
    [[maybe_unused]] bool has_coulomb, [[maybe_unused]] bool has_dipoles,
    [[maybe_unused]] auto const &aosoa, [[maybe_unused]] std::size_t i,
    [[maybe_unused]] std::size_t j) {
  PairDataFlags flags;
#ifdef ESPRESSO_GAY_BERNE
  flags.need_directors |= gay_berne_active(dist, ia_params);
#endif
#ifdef ESPRESSO_DIPOLES
  flags.need_directors |= has_dipoles;
#endif
#ifdef ESPRESSO_EXCLUSIONS
  flags.need_particle_pointers |=
      aosoa.has_exclusion(i) or aosoa.has_exclusion(j);
#endif
#ifdef ESPRESSO_THOLE
  flags.need_particle_pointers |= thole_active(ia_params, has_coulomb);
#endif
  return flags;
}
