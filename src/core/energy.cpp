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
#include "Observable_stat.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "cuda_interface.hpp"
#include "electrostatics_magnetostatics/magnetic_non_p3m_methods.hpp"
#include "energy_inline.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "short_range_loop.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#include <memory>

ActorList energyActors;

static std::shared_ptr<Observable_stat> calculate_energy_local() {

  auto obs_energy_ptr = std::make_shared<Observable_stat>(1);

  if (long_range_interactions_sanity_checks())
    return obs_energy_ptr;

  auto &obs_energy = *obs_energy_ptr;

#ifdef CUDA
  clear_energy_on_GPU();
#endif

  EspressoSystemInterface::Instance().update();

  // Compute the energies from the energyActors
  for (auto &energyActor : energyActors)
    energyActor->computeEnergy(espressoSystemInterface);

  on_observable_calc();

  auto const local_parts = cell_structure.local_particles();

  for (auto const &p : local_parts) {
    obs_energy.kinetic[0] += calc_kinetic_energy(p);
  }

  short_range_loop(
      [&obs_energy](Particle const &p1, int bond_id,
                    Utils::Span<Particle *> partners) {
        auto const &iaparams = *bonded_ia_params.at(bond_id);
        auto const result = calc_bonded_energy(iaparams, p1, partners);
        if (result) {
          obs_energy.bonded_contribution(bond_id)[0] += result.get();
          return false;
        }
        return true;
      },
      [&obs_energy](Particle const &p1, Particle const &p2, Distance const &d) {
        add_non_bonded_pair_energy(p1, p2, d.vec21, sqrt(d.dist2), d.dist2,
                                   obs_energy);
      },
      maximal_cutoff(), maximal_cutoff_bonded());

#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  obs_energy.coulomb[1] = Coulomb::calc_energy_long_range(local_parts);
#endif

#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  obs_energy.dipolar[1] = Dipole::calc_energy_long_range(local_parts);
#endif

  Constraints::constraints.add_energy(local_parts, get_sim_time(), obs_energy);

#ifdef CUDA
  auto const energy_host = copy_energy_from_GPU();
  if (!obs_energy.coulomb.empty())
    obs_energy.coulomb[1] += energy_host.coulomb;
  if (!obs_energy.dipolar.empty())
    obs_energy.dipolar[1] += energy_host.dipolar;
#endif

  obs_energy.mpi_reduce();
  return obs_energy_ptr;
}

REGISTER_CALLBACK_MAIN_RANK(calculate_energy_local)

std::shared_ptr<Observable_stat> calculate_energy() {
  return mpi_call(Communication::Result::main_rank, calculate_energy_local);
}

double calculate_current_potential_energy_of_system() {
  auto const obs_energy = calculate_energy();
  return obs_energy->accumulate(-obs_energy->kinetic[0]);
}

double observable_compute_energy() {
  auto const obs_energy = calculate_energy();
  return obs_energy->accumulate(0);
}
