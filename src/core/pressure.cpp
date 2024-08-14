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

#include "config/config.hpp"

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "electrostatics/coulomb.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "pressure_inline.hpp"
#include "short_range_loop.hpp"
#include "system/System.hpp"
#include "virtual_sites/relative.hpp"

#include <utils/Vector.hpp>
#include <utils/flatten.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <span>

namespace System {
std::shared_ptr<Observable_stat> System::calculate_pressure() {

  auto obs_pressure_ptr = std::make_shared<Observable_stat>(
      9ul, static_cast<std::size_t>(bonded_ias->get_next_key()),
      nonbonded_ias->get_max_seen_particle_type());

  if (long_range_interactions_sanity_checks()) {
    return obs_pressure_ptr;
  }

  auto &obs_pressure = *obs_pressure_ptr;

  on_observable_calc();

  auto const volume = box_geo->volume();
  auto const local_parts = cell_structure->local_particles();

  for (auto const &p : local_parts) {
    add_kinetic_virials(p, obs_pressure);
  }

  auto const coulomb_force_kernel = coulomb.pair_force_kernel();
  auto const coulomb_pressure_kernel = coulomb.pair_pressure_kernel();

  short_range_loop(
      [this, coulomb_force_kernel_ptr = get_ptr(coulomb_force_kernel),
       &obs_pressure](Particle const &p1, int bond_id,
                      std::span<Particle *> partners) {
        auto const &iaparams = *bonded_ias->at(bond_id);
        auto const result = calc_bonded_pressure_tensor(
            iaparams, p1, partners, *box_geo, coulomb_force_kernel_ptr);
        if (result) {
          auto const &tensor = result.value();
          /* pressure tensor part */
          for (std::size_t k = 0u; k < 3u; k++)
            for (std::size_t l = 0u; l < 3u; l++)
              obs_pressure.bonded_contribution(bond_id)[k * 3u + l] +=
                  tensor(k, l);

          return false;
        }
        return true;
      },
      [coulomb_force_kernel_ptr = get_ptr(coulomb_force_kernel),
       coulomb_pressure_kernel_ptr = get_ptr(coulomb_pressure_kernel), this,
       &obs_pressure](Particle const &p1, Particle const &p2,
                      Distance const &d) {
        auto const &ia_params =
            nonbonded_ias->get_ia_param(p1.type(), p2.type());
        add_non_bonded_pair_virials(p1, p2, d.vec21, sqrt(d.dist2), ia_params,
                                    *bonded_ias, coulomb_force_kernel_ptr,
                                    coulomb_pressure_kernel_ptr, obs_pressure);
      },
      *cell_structure, maximal_cutoff(), bonded_ias->maximal_cutoff());

#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  auto const coulomb_pressure = coulomb.calc_pressure_long_range(local_parts);
  std::ranges::copy(coulomb_pressure, obs_pressure.coulomb.begin() + 9u);
#endif
#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  Dipoles::get_dipoles().calc_pressure_long_range();
#endif

#ifdef VIRTUAL_SITES_RELATIVE
  if (!obs_pressure.virtual_sites.empty()) {
    auto const vs_pressure = vs_relative_pressure_tensor(*cell_structure);
    std::ranges::copy(Utils::flatten(vs_pressure),
                      obs_pressure.virtual_sites.begin());
  }
#endif

  obs_pressure.rescale(volume);

  obs_pressure.mpi_reduce();
  return obs_pressure_ptr;
}
} // namespace System
