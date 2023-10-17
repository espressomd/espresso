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
 *  Implementation of pressure.hpp.
 */

#include "pressure.hpp"
#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "electrostatics/coulomb.hpp"
#include "event.hpp"
#include "interactions.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "pressure_inline.hpp"
#include "short_range_loop.hpp"
#include "system/System.hpp"
#include "virtual_sites.hpp"

#include "config/config.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/range/algorithm/copy.hpp>

#include <memory>

std::shared_ptr<Observable_stat> calculate_pressure() {
  auto &system = System::get_system();
  auto &cell_structure = *system.cell_structure;
  auto const &box_geo = *system.box_geo;
  auto const &nonbonded_ias = *system.nonbonded_ias;

  auto obs_pressure_ptr = std::make_shared<Observable_stat>(
      9ul, static_cast<std::size_t>(::bonded_ia_params.get_next_key()),
      max_seen_particle_type);

  if (long_range_interactions_sanity_checks()) {
    return obs_pressure_ptr;
  }

  auto &obs_pressure = *obs_pressure_ptr;

  on_observable_calc();

  auto const volume = box_geo.volume();
  auto const local_parts = cell_structure.local_particles();

  for (auto const &p : local_parts) {
    add_kinetic_virials(p, obs_pressure);
  }

  auto const &coulomb = system.coulomb;
  auto const coulomb_force_kernel = coulomb.pair_force_kernel();
  auto const coulomb_pressure_kernel = coulomb.pair_pressure_kernel();

  short_range_loop(
      [coulomb_force_kernel_ptr = get_ptr(coulomb_force_kernel), &obs_pressure,
       &box_geo](Particle const &p1, int bond_id,
                 Utils::Span<Particle *> partners) {
        auto const &iaparams = *bonded_ia_params.at(bond_id);
        auto const result = calc_bonded_pressure_tensor(
            iaparams, p1, partners, box_geo, coulomb_force_kernel_ptr);
        if (result) {
          auto const &tensor = result.get();
          /* pressure tensor part */
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
              obs_pressure.bonded_contribution(bond_id)[k * 3 + l] +=
                  tensor(k, l);

          return false;
        }
        return true;
      },
      [coulomb_force_kernel_ptr = get_ptr(coulomb_force_kernel),
       coulomb_pressure_kernel_ptr = get_ptr(coulomb_pressure_kernel),
       &obs_pressure, &nonbonded_ias](Particle const &p1, Particle const &p2,
                                      Distance const &d) {
        auto const &ia_params =
            nonbonded_ias.get_ia_param(p1.type(), p2.type());
        add_non_bonded_pair_virials(p1, p2, d.vec21, sqrt(d.dist2), ia_params,
                                    coulomb_force_kernel_ptr,
                                    coulomb_pressure_kernel_ptr, obs_pressure);
      },
      cell_structure, maximal_cutoff(), maximal_cutoff_bonded());

#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  auto const coulomb_pressure = coulomb.calc_pressure_long_range(local_parts);
  boost::copy(coulomb_pressure, obs_pressure.coulomb.begin() + 9);
#endif
#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  Dipoles::get_dipoles().calc_pressure_long_range();
#endif

#ifdef VIRTUAL_SITES
  if (!obs_pressure.virtual_sites.empty()) {
    auto const vs_pressure = virtual_sites()->pressure_tensor();
    boost::copy(flatten(vs_pressure), obs_pressure.virtual_sites.begin());
  }
#endif

  obs_pressure.rescale(volume);

  obs_pressure.mpi_reduce();
  return obs_pressure_ptr;
}
