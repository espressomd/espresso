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
 *  Energy calculation.
 */

#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "cuda_interface.hpp"
#include "electrostatics_magnetostatics/magnetic_non_p3m_methods.hpp"
#include "energy_inline.hpp"
#include "event.hpp"
#include "forces.hpp"
#include <boost/range/numeric.hpp>

#include "short_range_loop.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

ActorList energyActors;

Observable_stat energy{};
Observable_stat total_energy{};

/************************************************************/

void init_energies(Observable_stat *stat) {
  size_t n_coulomb(0), n_dipolar(0);

#ifdef ELECTROSTATICS
  n_coulomb = Coulomb::energy_n();
#endif
#ifdef DIPOLES
  n_dipolar = Dipole::energy_n();
#endif

  stat->realloc_and_clear(n_coulomb, n_dipolar, 0, 1);
  stat->init_status = 0;
}

/************************************************************/

void master_energy_calc() {
  mpi_gather_stats(GatherStats::energy, reinterpret_cast<void *>(&total_energy),
                   nullptr, nullptr, nullptr);

  total_energy.init_status = 1;
}

/************************************************************/

void energy_calc(Observable_stat *result, const double time) {
  if (!interactions_sanity_checks())
    return;

  init_energies(&energy);

#ifdef CUDA
  clear_energy_on_GPU();
#endif

  EspressoSystemInterface::Instance().update();

  // Compute the energies from the energyActors
  for (auto &energyActor : energyActors)
    energyActor->computeEnergy(espressoSystemInterface);

  on_observable_calc();

  for (auto const &p : cell_structure.local_particles()) {
    add_kinetic_energy(p);
  }

  short_range_loop(
      [](Particle &p) { add_single_particle_energy(p); },
      [](Particle const &p1, Particle const &p2, Distance const &d) {
        add_non_bonded_pair_energy(p1, p2, d.vec21, sqrt(d.dist2), d.dist2);
      });

  calc_long_range_energies(cell_structure.local_particles());

  auto local_parts = cell_structure.local_particles();
  Constraints::constraints.add_energy(local_parts, time, energy);

#ifdef CUDA
  copy_energy_from_GPU();
#endif

  /* gather data */
  energy.reduce(result);
}

/************************************************************/

void calc_long_range_energies(const ParticleRange &particles) {
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  Coulomb::calc_energy_long_range(energy, particles);
#endif /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  Dipole::calc_energy_long_range(energy, particles);
#endif /* ifdef DIPOLES */
}

double calculate_current_potential_energy_of_system() {
  // calculate potential energy
  if (total_energy.init_status == 0) {
    init_energies(&total_energy);
    master_energy_calc();
  }

  return boost::accumulate(total_energy.data, -total_energy.data[0]);
}
