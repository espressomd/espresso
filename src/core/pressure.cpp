/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "interactions.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "pressure_inline.hpp"
#include "virtual_sites.hpp"

#include "short_range_loop.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/range/algorithm/copy.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>

static std::shared_ptr<Observable_stat> calculate_pressure_local() {

  auto obs_pressure_ptr = std::make_shared<Observable_stat>(9);

  if (long_range_interactions_sanity_checks())
    return obs_pressure_ptr;

  auto &obs_pressure = *obs_pressure_ptr;

  on_observable_calc();

  auto const volume = box_geo.volume();
  auto const local_parts = cell_structure.local_particles();

  for (auto const &p : local_parts) {
    add_kinetic_virials(p, obs_pressure);
  }

  short_range_loop(
      [&obs_pressure](Particle const &p1, int bond_id,
                      Utils::Span<Particle *> partners) {
        auto const &iaparams = *bonded_ia_params.at(bond_id);
        auto const result = calc_bonded_pressure_tensor(iaparams, p1, partners);
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
      [&obs_pressure](Particle const &p1, Particle const &p2,
                      Distance const &d) {
        add_non_bonded_pair_virials(p1, p2, d.vec21, sqrt(d.dist2),
                                    obs_pressure);
      },
      maximal_cutoff(), maximal_cutoff_bonded());

#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  auto const coulomb_pressure = Coulomb::calc_pressure_long_range(local_parts);
  boost::copy(coulomb_pressure, obs_pressure.coulomb.begin() + 9);
#endif
#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  Dipole::calc_pressure_long_range();
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

REGISTER_CALLBACK_MAIN_RANK(calculate_pressure_local)

std::shared_ptr<Observable_stat> calculate_pressure() {
  return mpi_call(Communication::Result::main_rank, calculate_pressure_local);
}

Utils::Vector9d observable_compute_pressure_tensor() {
  auto const obs_pressure = calculate_pressure();
  Utils::Vector9d pressure_tensor{};
  for (std::size_t j = 0; j < 9; j++) {
    pressure_tensor[j] = obs_pressure->accumulate(0, j);
  }
  return pressure_tensor;
}
