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
  int n_pre, n_non_bonded, n_dipolar(0);

  n_pre = 1;
  n_non_bonded = (max_seen_particle_type * (max_seen_particle_type + 1)) / 2;

#ifdef ELECTROSTATICS
  auto const n_coulomb = Coulomb::energy_n();
#else
  int const n_coulomb = 0;
#endif
#ifdef DIPOLES
  Dipole::energy_n(n_dipolar);
#endif

  obsstat_realloc_and_clear(stat, n_pre, bonded_ia_params.size(), n_non_bonded,
                            n_coulomb, n_dipolar, 0, 1);
  stat->init_status = 0;
}

/************************************************************/

void master_energy_calc() {
  mpi_gather_stats(1, total_energy.data.data(), nullptr, nullptr, nullptr);

  total_energy.init_status = 1;
}

/************************************************************/

void energy_calc(double *result, const double time) {
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
  MPI_Reduce(energy.data.data(), result, energy.data.size(), MPI_DOUBLE,
             MPI_SUM, 0, comm_cart);
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
