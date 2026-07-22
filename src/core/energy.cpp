/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "bond_energy_kokkos.hpp"
#include "cell_system/CellStructure.hpp"
#include "constraints/Constraints.hpp"
#include "energy_cabana.hpp"
#include "energy_inline.hpp"
#include "integrators/Propagation.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "short_range_cabana.hpp"
#include "short_range_loop.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"

#include "electrostatics/coulomb.hpp"
#include "magnetostatics/dipoles.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace System {

std::shared_ptr<Observable_stat> System::calculate_energy() {

  auto obs_energy_ptr = std::make_shared<Observable_stat>(
      1ul, static_cast<std::size_t>(bonded_ias->get_next_key()),
      nonbonded_ias->get_max_seen_particle_type());

  if (long_range_interactions_sanity_checks()) {
    return obs_energy_ptr;
  }

  auto &obs_energy = *obs_energy_ptr;
#if defined(ESPRESSO_CUDA) and                                                 \
    (defined(ESPRESSO_ELECTROSTATICS) or defined(ESPRESSO_DIPOLES))
  gpu->clear_energy_on_device();
  gpu->update();
#endif
  on_observable_calc();

  auto const local_parts = cell_structure->local_particles();

  for (auto const &p : local_parts) {
    obs_energy.kinetic_lin[0] += translational_kinetic_energy(p);
    obs_energy.kinetic_rot[0] += rotational_kinetic_energy(p);
  }

  auto const coulomb_kernel = coulomb.pair_energy_kernel();
  auto const dipoles_kernel = dipoles.pair_energy_kernel();

#ifdef ESPRESSO_CALIPER
  CALI_MARK_BEGIN("cabana_short_range");
#endif
  VerletCriterion<> const verlet_criterion{*this,
                                           cell_structure->get_verlet_skin(),
                                           get_interaction_range(),
                                           coulomb.cutoff(),
                                           dipoles.cutoff(),
                                           inactive_cutoff};
  update_cabana_state(*cell_structure, verlet_criterion,
                      get_interaction_range(), propagation->integ_switch);

  EnergyBinLayout layout{
      static_cast<std::size_t>(bonded_ias->get_next_key()),
      std::size_t(nonbonded_ias->get_max_seen_particle_type() + 1)};

  using execution_space = Kokkos::DefaultHostExecutionSpace;
  Kokkos::View<double **, Kokkos::LayoutRight, execution_space> local_energy(
      "local_energy", execution_space().concurrency(), layout.total);
  auto const &unique_particles = cell_structure->get_unique_particles();
  auto const n_particles = unique_particles.size();
  Kokkos::View<int *, Kokkos::LayoutRight, execution_space> mol_id("mol_id",
                                                                   n_particles);
  for (std::size_t i = 0; i < n_particles; ++i) {
    mol_id(i) = unique_particles[i]->mol_id();
  }

  // Non Bonded energies
  EnergyKernel pair_e_kernel{*bonded_ias,
                             *nonbonded_ias,
                             coulomb,
                             get_ptr(coulomb_kernel),
                             get_ptr(dipoles_kernel),
                             *box_geo,
                             cell_structure->get_unique_particles(),
                             local_energy,
                             layout,
                             cell_structure->get_aosoa(),
                             mol_id,
                             maximal_cutoff()};

  // Bonded energies: write a BondsEnergyKernelData + *BondsEnergyKernel
  auto &bs = cell_structure->bond_state();
  BondsEnergyKernelData bonds_e_data{*bonded_ias, *box_geo, local_energy,
                                     layout, cell_structure->get_aosoa()};
  PairBondsEnergyKernel pair_be_kernel{bonds_e_data, bs.pair_list, bs.pair_ids,
                                       get_ptr(coulomb_kernel)};
  AngleBondsEnergyKernel angle_be_kernel{bonds_e_data, bs.angle_list,
                                         bs.angle_ids};
  DihedralBondsEnergyKernel dih_be_kernel{bonds_e_data, bs.dihedral_list,
                                          bs.dihedral_ids};

  cabana_short_range(pair_be_kernel, angle_be_kernel, dih_be_kernel,
                     pair_e_kernel, *cell_structure, get_interaction_range(),
                     bonded_ias->maximal_cutoff(), verlet_criterion,
                     propagation->integ_switch);

  reduce_cabana_energy(local_energy, layout, obs_energy, *bonded_ias,
                       nonbonded_ias->get_max_seen_particle_type() + 1);
#ifdef ESPRESSO_CALIPER
  CALI_MARK_END("cabana_short_range");
#endif

#ifdef ESPRESSO_ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  obs_energy.coulomb[1] = coulomb.calc_energy_long_range();
#endif

#ifdef ESPRESSO_DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  obs_energy.dipolar[1] = dipoles.calc_energy_long_range();
#endif

  constraints->add_energy(local_parts, get_sim_time(), obs_energy);

#if defined(ESPRESSO_CUDA) and                                                 \
    (defined(ESPRESSO_ELECTROSTATICS) or defined(ESPRESSO_DIPOLES))
  auto const energy_host = gpu->copy_energy_to_host();
  if (!obs_energy.coulomb.empty())
    obs_energy.coulomb[1] += static_cast<double>(energy_host.coulomb);
  if (!obs_energy.dipolar.empty())
    obs_energy.dipolar[1] += static_cast<double>(energy_host.dipolar);
#endif

  obs_energy.mpi_reduce();
  return obs_energy_ptr;
  // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
}

double System::particle_short_range_energy_contribution(int pid) {
  if (cell_structure->get_resort_particles()) {
    cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());
  }

  auto ret = 0.0;
  if (auto const p = cell_structure->get_local_particle(pid)) {
    auto const coulomb_kernel = coulomb.pair_energy_kernel();
    auto kernel = [&ret, this](Particle const &p, Particle const &p1,
                               Utils::Vector3d const &vec) {
#ifdef ESPRESSO_EXCLUSIONS
      if (not do_nonbonded(p, p1))
        return;
#endif
      auto const &ia_params = nonbonded_ias->get_ia_param(p.type(), p1.type());
      // Add energy for current particle pair to result
      ret += calc_non_bonded_pair_energy(p, p1, ia_params, vec, vec.norm(),
                                         *bonded_ias, coulomb, nullptr);
    };
    cell_structure->run_on_particle_short_range_neighbors(*p, kernel);
  }
  return ret;
}

std::optional<double> System::particle_bond_energy(int pid, int bond_id,
                                                   std::vector<int> partners) {
  if (cell_structure->get_resort_particles()) {
    cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());
  }
  Particle const *p = cell_structure->get_local_particle(pid);
  if (not p or p->is_ghost())
    return {}; // not available on this MPI rank or ghost
  auto const &iaparams = *bonded_ias->at(bond_id);
  try {
    auto resolved_partners = cell_structure->resolve_bond_partners(partners);
    auto const coulomb_kernel = coulomb.pair_energy_kernel();
    return calc_bonded_energy(
        iaparams, *p,
        std::span(resolved_partners.data(), resolved_partners.size()), *box_geo,
        get_ptr(coulomb_kernel));
  } catch (const BondResolutionError &) {
    bond_broken_error(p->id(), partners);
    return {};
  }
}

} // namespace System
