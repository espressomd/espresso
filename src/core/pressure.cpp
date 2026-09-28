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
#include "bond_pressure_kokkos.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/CellStructure.hpp"
#include "electrostatics/coulomb.hpp"
#include "integrators/Propagation.hpp"
#include "kokkos_helpers.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_reduction.hpp"
#include "pressure_cabana.hpp"
#include "pressure_inline.hpp"
#include "short_range_cabana.hpp"
#include "short_range_verlet.hpp"
#include "system/System.hpp"
#include "virtual_sites/relative.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>
#include <utils/matrix.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>

struct PressureObservable {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  std::unique_ptr<Observable_stat> observable;
  Kokkos::View<double **, Kokkos::LayoutRight, execution_space> local;
};

namespace System {

Observable_stat const &System::calculate_pressure() {

  if (not m_obs_pressure) {
    auto const n_threads = PressureObservable::execution_space{}.concurrency();
    m_obs_pressure = std::make_shared<PressureObservable>();
    m_obs_pressure->observable = std::make_unique<Observable_stat>(9ul, 0ul, 0);
    m_obs_pressure->local =
        decltype(m_obs_pressure->local)("local_pressure", n_threads, 9ul);
  }

  auto &local_pressure = m_obs_pressure->local;
  auto &obs_pressure = *m_obs_pressure->observable;
  obs_pressure.reset(static_cast<std::size_t>(bonded_ias->get_next_key()),
                     nonbonded_ias->get_max_seen_particle_type());

  if (long_range_interactions_sanity_checks()) {
    return obs_pressure;
  }

  on_observable_calc();

  auto const volume = box_geo->volume();

  // Kinetic virial — reduction over local particles
  auto const kinetic = reduce_over_local_particles<Utils::Matrix<double, 3, 3>>(
      *cell_structure,
      [](Utils::Matrix<double, 3, 3> &acc, Particle const &p) {
        if (!p.is_virtual()) {
          auto const vel = Utils::Vector3d(p.v());
          acc += Utils::tensor_product(vel, p.mass() * vel);
        }
      },
      [](auto &a, auto const &b) { a += b; });
  std::ranges::copy(Utils::flatten(kinetic), obs_pressure.kinetic_lin.begin());

  auto const coulomb_force_kernel = coulomb.pair_force_kernel();
  auto const coulomb_pressure_kernel = coulomb.pair_pressure_kernel();
  auto const dipoles_pressure_kernel = dipoles.pair_pressure_kernel();

  // Factory instead of an eager criterion: construction fills an O(n_types^2)
  // cutoff table, so it only runs on the link-cell fallback path.
  auto const make_verlet_criterion = [&] {
    return VerletCriterion<>{*this,
                             cell_structure->get_verlet_skin(),
                             get_interaction_range(),
                             coulomb.cutoff(),
                             dipoles.cutoff(),
                             inactive_cutoff};
  };
  update_verlet_state(*this, inactive_cutoff);
#ifdef ESPRESSO_ELECTROSTATICS
  // Refresh the pack-owned charge column once, guarded by an active coulomb
  // actor (the pressure pair kernel reads it contiguously).
  if (coulomb.impl->solver) {
    refresh_pack_charges(*cell_structure);
  }
#endif // ESPRESSO_ELECTROSTATICS
#ifdef ESPRESSO_DIPOLES
  if (dipoles.impl->solver) {
    refresh_pack_dipm(*cell_structure);
  }
#endif // ESPRESSO_DIPOLES

  PressureBinLayout layout{
      static_cast<std::size_t>(bonded_ias->get_next_key()),
      std::size_t(nonbonded_ias->get_max_seen_particle_type() + 1)};

  using exec = Kokkos::DefaultHostExecutionSpace;
  if (local_pressure.extent(1) != layout.total * 9ul) {
    Kokkos::realloc(Kokkos::WithoutInitializing, local_pressure,
                    exec{}.concurrency(), layout.total * 9ul);
  }
  kokkos_deep_copy(exec{}, local_pressure, 0.);

  auto const &unique_particles = cell_structure->get_unique_particles();
  auto const n_particles = unique_particles.size();
  Kokkos::View<int *, Kokkos::LayoutRight, exec> mol_id("mol_id", n_particles);
  for (std::size_t i = 0; i < n_particles; ++i) {
    mol_id(i) = unique_particles[i]->mol_id();
  }

  PressureKernel pair_p_kernel{*bonded_ias,
                               *nonbonded_ias,
                               coulomb,
                               get_ptr(coulomb_force_kernel),
                               get_ptr(coulomb_pressure_kernel),
                               get_ptr(dipoles_pressure_kernel),
                               *box_geo,
#ifdef ESPRESSO_DPD
                               thermostat->dpd.get(),
#endif
                               cell_structure->get_unique_particles(),
                               local_pressure,
                               layout,
                               cell_structure->get_aosoa(),
                               mol_id,
                               maximal_cutoff(),
                               thermostat->thermo_switch};

  auto &bs = cell_structure->bond_state();
  BondsPressureKernelData bonds_p_data{*bonded_ias, *box_geo, local_pressure,
                                       layout, cell_structure->get_aosoa()};
  PairBondsPressureKernel pair_bp_kernel{
      bonds_p_data, bs.pair_list, bs.pair_ids, get_ptr(coulomb_force_kernel)};
  AngleBondsPressureKernel angle_bp_kernel{bonds_p_data, bs.angle_list,
                                           bs.angle_ids};
  DihedralBondsPressureKernel dih_bp_kernel{bonds_p_data, bs.dihedral_list,
                                            bs.dihedral_ids};

  cabana_short_range(pair_bp_kernel, angle_bp_kernel, dih_bp_kernel,
                     pair_p_kernel, *cell_structure, get_interaction_range(),
                     bonded_ias->maximal_cutoff(), make_verlet_criterion,
                     propagation->integ_switch);

  reduce_cabana_pressure(local_pressure, layout, obs_pressure, *bonded_ias,
                         nonbonded_ias->get_max_seen_particle_type() + 1);

#ifdef ESPRESSO_ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  auto const coulomb_pressure = coulomb.calc_pressure_long_range();
  std::ranges::copy(coulomb_pressure, obs_pressure.coulomb.begin() + 9u);
#endif
#ifdef ESPRESSO_DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  auto const dipoles_pressure = dipoles.calc_pressure_long_range();
  std::ranges::copy(dipoles_pressure, obs_pressure.dipolar.begin() + 9u);
#endif

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  if (!obs_pressure.virtual_sites.empty()) {
    // vs_relative_pressure_tensor reads virtual-site particle forces; ensure
    // every particle has a valid ParticleStore row (on_observable_calc above
    // may have resorted). O(1) when clean; rank-local.
    cell_structure->ensure_particle_store_synchronized();
    auto const vs_pressure = vs_relative_pressure_tensor(*cell_structure);
    std::ranges::copy(Utils::flatten(vs_pressure),
                      obs_pressure.virtual_sites.begin());
  }
#endif

#ifdef ESPRESSO_BOND_CONSTRAINT
  if (propagation->is_inertial() and bonded_ias->get_n_rigid_bonds() >= 1) {
    // rigid_bond_virial was accumulated bond-by-bond inside
    // correct_position_shake() during the last integration step, so it is
    // already correct regardless of how many rigid bonds a particle
    // participates in; only the deferred 1/dt^2 factor is applied here.
    auto const sq_dt = Utils::sqr(get_time_step());
    auto const &rigid_bond_virial = bonded_ias->rigid_bond_virial;
    for (std::size_t bond_id = 0; bond_id < rigid_bond_virial.size();
         ++bond_id) {
      auto const stress = rigid_bond_virial[bond_id] / sq_dt;
      auto dest = obs_pressure.bonded_contribution(static_cast<int>(bond_id));
      for (std::size_t k = 0; k < 9u; ++k)
        dest[k] += stress[k];
    }
  }
#endif // ESPRESSO_BOND_CONSTRAINT

  obs_pressure.rescale(volume);

  obs_pressure.mpi_reduce();
  return obs_pressure;
}

} // namespace System
