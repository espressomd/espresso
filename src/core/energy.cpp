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
 *  Energy calculation.
 */

#include "Observable_stat.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "energy_inline.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "short_range_loop.hpp"
#include "system/System.hpp"

#include "electrostatics/coulomb.hpp"
#include "magnetostatics/dipoles.hpp"

#include <utils/Span.hpp>
#include <utils/mpi/iall_gatherv.hpp>

#include <memory>

std::shared_ptr<Observable_stat> calculate_energy() {

  auto obs_energy_ptr = std::make_shared<Observable_stat>(1);

  if (long_range_interactions_sanity_checks()) {
    return obs_energy_ptr;
  }

  auto &obs_energy = *obs_energy_ptr;
#if defined(CUDA) and (defined(ELECTROSTATICS) or defined(DIPOLES))
  auto &gpu_particle_data = System::get_system().gpu;
  gpu_particle_data.clear_energy_on_device();
  gpu_particle_data.update();
#endif
  on_observable_calc();

  auto const local_parts = cell_structure.local_particles();

  for (auto const &p : local_parts) {
    obs_energy.kinetic[0] += calc_kinetic_energy(p);
  }

  auto const coulomb_kernel = Coulomb::pair_energy_kernel();
  auto const dipoles_kernel = Dipoles::pair_energy_kernel();

  short_range_loop(
      [&obs_energy, coulomb_kernel_ptr = coulomb_kernel.get_ptr()](
          Particle const &p1, int bond_id, Utils::Span<Particle *> partners) {
        auto const &iaparams = *bonded_ia_params.at(bond_id);
        auto const result =
            calc_bonded_energy(iaparams, p1, partners, coulomb_kernel_ptr);
        if (result) {
          obs_energy.bonded_contribution(bond_id)[0] += result.get();
          return false;
        }
        return true;
      },
      [&obs_energy, coulomb_kernel_ptr = coulomb_kernel.get_ptr(),
       dipoles_kernel_ptr = dipoles_kernel.get_ptr()](
          Particle const &p1, Particle const &p2, Distance const &d) {
        add_non_bonded_pair_energy(p1, p2, d.vec21, sqrt(d.dist2), d.dist2,
                                   coulomb_kernel_ptr, dipoles_kernel_ptr,
                                   obs_energy);
      },
      maximal_cutoff(n_nodes), maximal_cutoff_bonded());

#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  obs_energy.coulomb[1] = Coulomb::calc_energy_long_range(local_parts);
#endif

#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  obs_energy.dipolar[1] = Dipoles::calc_energy_long_range(local_parts);
#endif

  Constraints::constraints.add_energy(local_parts, get_sim_time(), obs_energy);

#if defined(CUDA) and (defined(ELECTROSTATICS) or defined(DIPOLES))
  auto const energy_host = gpu_particle_data.copy_energy_to_host();
  if (!obs_energy.coulomb.empty())
    obs_energy.coulomb[1] += static_cast<double>(energy_host.coulomb);
  if (!obs_energy.dipolar.empty())
    obs_energy.dipolar[1] += static_cast<double>(energy_host.dipolar);
#endif

  obs_energy.mpi_reduce();
  return obs_energy_ptr;
  // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
}

REGISTER_CALLBACK_MAIN_RANK(calculate_energy)

double mpi_calculate_potential_energy() {
  auto const obs = mpi_call(Communication::Result::main_rank, calculate_energy);
  return obs->accumulate(-obs->kinetic[0]);
}

double mpi_observable_compute_energy() {
  auto const obs = mpi_call(Communication::Result::main_rank, calculate_energy);
  return obs->accumulate(0);
}

double particle_short_range_energy_contribution(int pid) {
  double ret = 0.0;

  if (cell_structure.get_resort_particles()) {
    cells_update_ghosts(global_ghost_flags());
  }

  if (auto const p = cell_structure.get_local_particle(pid)) {
    auto const coulomb_kernel = Coulomb::pair_energy_kernel();
    auto kernel = [&ret, coulomb_kernel_ptr = coulomb_kernel.get_ptr()](
                      Particle const &p, Particle const &p1,
                      Utils::Vector3d const &vec) {
#ifdef EXCLUSIONS
      if (not do_nonbonded(p, p1))
        return;
#endif
      auto const &ia_params = get_ia_param(p.type(), p1.type());
      // Add energy for current particle pair to result
      ret += calc_non_bonded_pair_energy(p, p1, ia_params, vec, vec.norm(),
                                         coulomb_kernel_ptr);
    };
    cell_structure.run_on_particle_short_range_neighbors(*p, kernel);
  }
  return ret;
}

#ifdef DIPOLE_FIELD_TRACKING
void calc_long_range_fields() {
  auto particles = cell_structure.local_particles();
  Dipoles::calc_long_range_field(particles);
}

REGISTER_CALLBACK(calc_long_range_fields)

void mpi_calc_long_range_fields() { mpi_call_all(calc_long_range_fields); }
#endif
