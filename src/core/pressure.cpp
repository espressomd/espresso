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

#include "Observable_stat.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "pressure_inline.hpp"
#include "virtual_sites.hpp"
#include <boost/range/algorithm/copy.hpp>

#include "short_range_loop.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

/** Scalar pressure of the system */
Observable_stat obs_scalar_pressure{1};
/** Pressure tensor of the system */
Observable_stat obs_pressure_tensor{9};

nptiso_struct nptiso = {0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        0,
                        {NPTGEOM_XDIR, NPTGEOM_YDIR, NPTGEOM_ZDIR},
                        0,
                        false,
                        0};

/** Calculate long-range virials (P3M, ...). */
void calc_long_range_virials(const ParticleRange &particles) {
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  Coulomb::calc_pressure_long_range(obs_scalar_pressure, obs_pressure_tensor,
                                    particles);
#endif
#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  Dipole::calc_pressure_long_range();
#endif
}

static void add_single_particle_virials(Particle &p) {
  cell_structure.execute_bond_handler(p, add_bonded_pressure_tensor);
}

void pressure_calc() {
  auto const volume = box_geo.volume();

  if (!interactions_sanity_checks())
    return;

  obs_scalar_pressure.resize();
  obs_pressure_tensor.resize();

  on_observable_calc();

  for (auto const &p : cell_structure.local_particles()) {
    add_kinetic_virials(p);
  }

  short_range_loop([](Particle &p) { add_single_particle_virials(p); },
                   [](Particle &p1, Particle &p2, Distance const &d) {
                     add_non_bonded_pair_virials(p1, p2, d.vec21,
                                                 sqrt(d.dist2));
                   });

  calc_long_range_virials(cell_structure.local_particles());

#ifdef VIRTUAL_SITES
  if (!obs_scalar_pressure.virtual_sites.empty()) {
    auto const vs_pressure_tensor = virtual_sites()->pressure_tensor();

    obs_scalar_pressure.virtual_sites[0] += trace(vs_pressure_tensor);
    boost::copy(flatten(vs_pressure_tensor),
                obs_pressure_tensor.virtual_sites.begin());
  }
#endif

  /* rescale kinetic energy (=ideal contribution) */
  obs_scalar_pressure.rescale(3.0 * volume);

  obs_pressure_tensor.rescale(volume);

  /* gather data */
  obs_scalar_pressure.reduce(comm_cart);
  obs_pressure_tensor.reduce(comm_cart);
}

void update_pressure() { mpi_gather_stats(GatherStats::pressure); }

Utils::Vector9d observable_compute_pressure_tensor() {
  update_pressure();
  Utils::Vector9d pressure_tensor{};
  for (size_t j = 0; j < 9; j++) {
    pressure_tensor[j] = obs_pressure_tensor.accumulate(0, j);
  }
  return pressure_tensor;
}
