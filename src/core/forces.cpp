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
#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/for_each_particle.hpp"
#include "cells.hpp"
#include "collision_detection/CollisionDetection.hpp"
#include "communication.hpp"
#include "constraints/Constraints.hpp"
#include "electrostatics/icc.hpp"
#include "forces_init.hpp"
#include "forces_inline.hpp"
#include "galilei/ComFixed.hpp"
#include "ghosts.hpp"
#include "immersed_boundary/ImmersedBoundaries.hpp"
#include "integrators/Propagation.hpp"
#include "lb/particle_coupling.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "rotation.hpp"
#include "short_range_cabana.hpp"
#include "short_range_loop.hpp"
#include "short_range_verlet.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "virtual_sites/com.hpp"
#include "virtual_sites/relative.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#ifdef ESPRESSO_CALIPER
#include "caliper_utils.hpp"
#endif

#include <Cabana_Core.hpp>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <span>
#include <variant>

/**
 * @brief Eligibility check for the split-phase ghost force reduction.
 *
 * The overlap is safe only when ALL of the following hold:
 *
 *  1. More than one MPI rank: at 1 rank the reduce is a local-only copy,
 *     there is nothing to hide, so the split adds overhead for no gain.
 *  2. ComFixed is inactive: ComFixed applies a per-type force correction that
 *     must see the final (post-reduce) forces, which arrive on the blocking
 *     path.
 *  3. force_cap == 0: force capping requires the final forces (post-reduce).
 *  4. Integrator is VV or symplectic Euler: only these inertial methods use
 *     step_2; steepest-descent, BD, SD, and NPT are not eligible (NPT virial
 *     is not accumulated during step_2 overlap; BD/SD have no step_2
 *     half-kick).
 *  5. NPT propagation not in use: the NPT step_2 updates the box pressure and
 *     box length using the virial gathered from forces; it must run after the
 *     full reduce.  NPT also implies integ_switch != NVT/SE, so condition 4
 *     already excludes it, but we check explicitly for clarity.
 *  6. Dipole field tracking not in use
 *  7. LB-tracer arm not active: TRANS_LB_TRACER particles have forces
 *     transferred to the LB fluid in the integrate loop *after*
 * calculate_forces and before step_2; with the split-phase path the fluid
 * coupling would run while the reduction is in flight, leading to incorrect
 * force accumulation.
 */
static bool ghost_reduce_overlap_eligible(System::System const &system) {
  // 1. Must have more than one MPI rank.
  if (::comm_cart.size() <= 1)
    return false;
  // 2. ComFixed must be inactive.
  if (not system.comfixed->get_fixed_types().empty())
    return false;
  // 3. Force capping must be off.
  if (system.get_force_cap() != 0.)
    return false;
  // 4. NPT propagation not in use.
#ifdef ESPRESSO_NPT
  if (system.propagation->used_propagations &
      PropagationMode::TRANS_LANGEVIN_NPT)
    return false;
#endif
  // 5. LB-tracer arm not active.
#ifdef ESPRESSO_VIRTUAL_SITES_INERTIALESS_TRACERS
  if (system.propagation->used_propagations & PropagationMode::TRANS_LB_TRACER)
    return false;
#endif
  // 6. Dipole field tracking not active.
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (system.dipoles.impl->solver.has_value()) {
    return false;
  }
#endif
  // 7. Integrator must be VV or symplectic Euler (inertial with step_2
  // half-kick).
  auto const integ = system.propagation->integ_switch;
  return integ == INTEG_METHOD_NVT or integ == INTEG_METHOD_SYMPLECTIC_EULER;
}

static void force_capping(CellStructure &cell_structure, double force_cap) {
  if (force_cap > 0.) {
    auto const force_cap_sq = Utils::sqr(force_cap);
    cell_structure.for_each_local_particle(
        [&force_cap, &force_cap_sq](Particle &p) {
          auto const force_sq = p.force().norm2();
          if (force_sq > force_cap_sq) {
            p.force() *= force_cap / std::sqrt(force_sq);
          }
        });
  }
}

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
static void reinit_dip_fld(CellStructure const &cell_structure) {
  cell_structure.for_each_local_particle(
      [](Particle &p) { p.dip_fld() = {0., 0., 0.}; });
}
#endif

static BondsKernelData
create_kokkos_bonds_kernel_data(System::System const &system) {
  auto scatter_force = system.cell_structure->get_scatter_force();
#ifdef ESPRESSO_NPT
  auto scatter_virial = system.cell_structure->get_scatter_virial();
#endif
  auto const &aosoa = system.cell_structure->get_aosoa();
  return /* BondsKernelData */ {*system.bonded_ias,
                                *system.bond_breakage,
                                *system.box_geo,
                                scatter_force,
#ifdef ESPRESSO_NPT
                                scatter_virial,
#endif
                                aosoa,
                                !system.bond_breakage->breakage_specs.empty()};
}

static ForcesKernel create_cabana_neighbor_kernel(
    System::System const &system, Utils::Vector3d *virial,
    auto const &elc_kernel, auto const &coulomb_kernel,
    auto const &dipoles_kernel, auto const &coulomb_u_kernel) {

  auto const &unique_particles = system.cell_structure->get_unique_particles();

  auto scatter_force = system.cell_structure->get_scatter_force();
#ifdef ESPRESSO_ROTATION
  auto scatter_torque = system.cell_structure->get_scatter_torque();
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  auto scatter_dip_fld = system.cell_structure->get_scatter_dip_fld();
#endif
#ifdef ESPRESSO_NPT
  auto scatter_virial = system.cell_structure->get_scatter_virial();
#endif
  auto const &aosoa = system.cell_structure->get_aosoa();

  return /* ForcesKernel */ {*system.bonded_ias,
                             *system.nonbonded_ias,
                             get_ptr(coulomb_kernel),
                             get_ptr(dipoles_kernel),
                             get_ptr(elc_kernel),
                             get_ptr(coulomb_u_kernel),
                             system.coulomb,
                             *system.thermostat,
                             *system.box_geo,
                             unique_particles,
                             scatter_force,
#ifdef ESPRESSO_ROTATION
                             scatter_torque,
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
                             scatter_dip_fld,
#endif
#ifdef ESPRESSO_NPT
                             virial,
                             scatter_virial,
#endif
                             aosoa,
                             system.maximal_cutoff()};
}

// Single construction-and-launch site for SpecializedForcesKernel, shared by
// the with- and without-coulomb branches of the dispatch below.
template <bool HasCoulomb>
static ShortRangeVerletPairLoop make_specialized_verlet_pair_loop(
    InteractionsNonBonded const &nonbonded_ias,
    CellStructure::AoSoA_pack const &aosoa,
    CellStructure::ScatterForce scatter_force,
    CuboidMinimumImage const &minimum_image, double const max_cutoff_sq
#ifdef ESPRESSO_ELECTROSTATICS
    ,
    Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_ptr = nullptr
#ifdef ESPRESSO_P3M
    ,
    CoulombP3M const *p3m = nullptr
#endif
#endif
) {
  return [=, &nonbonded_ias, &aosoa](CellStructure::ListType const &verlet_list,
                                     std::size_t const n) {
    SpecializedForcesKernel<HasCoulomb> const kernel{nonbonded_ias,
                                                     aosoa,
                                                     scatter_force,
                                                     verlet_list.counts,
                                                     verlet_list.neighbors,
                                                     minimum_image,
                                                     max_cutoff_sq
#ifdef ESPRESSO_ELECTROSTATICS
                                                     ,
                                                     coulomb_ptr
#ifdef ESPRESSO_P3M
                                                     ,
                                                     p3m
#endif
#endif
    };
    Kokkos::parallel_for("specialized_nonbonded_pairs",
                         Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(
                             std::size_t{0}, n),
                         kernel);
  };
}

// Build the compile-time-specialized Verlet pair loop when the active feature
// set is covered by SpecializedForcesKernel: cuboid box, no NPT virial, no
// dipolar or ELC kernel, only allowlisted (central-radial) pair potentials, no
// Thole pair, and no particle with an exclusion. Returns an empty
// ShortRangeVerletPairLoop otherwise, leaving cabana_short_range on the
// generic ForcesKernel path. The specialized kernel is bitwise-identical to
// the generic one on these systems.
static ShortRangeVerletPairLoop create_specialized_verlet_pair_loop(
    System::System const &system,
    [[maybe_unused]] Utils::Vector3d const *virial,
    [[maybe_unused]] auto const &elc_kernel, auto const &coulomb_kernel,
    [[maybe_unused]] auto const &dipoles_kernel) {
  if (system.box_geo->type() != BoxType::CUBOID)
    return {};
#ifdef ESPRESSO_NPT
  if (virial != nullptr)
    return {};
#endif
#ifdef ESPRESSO_DIPOLES
  if (get_ptr(dipoles_kernel) != nullptr)
    return {};
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  if (get_ptr(elc_kernel) != nullptr)
    return {};
#endif
  auto const &nonbonded_ias = *system.nonbonded_ias;
  // Allowlist over the aggregated pair-potential mask (O(1), maintained by
  // recalc_maximal_cutoffs). Any type pair with a potential the specialized
  // kernel does not compute -- Gay-Berne, DPD (in which case the DPD
  // thermostat could act on the pair), or any future addition -- falls back
  // to the generic kernel by default.
  if ((nonbonded_ias.combined_active_pair_mask() &
       ~specialized_kernel_pair_mask) != 0u)
    return {};
#ifdef ESPRESSO_THOLE
  // Thole damping is not in the pair-potential mask; check its own aggregate.
  if (nonbonded_ias.any_thole_configured())
    return {};
#endif
  auto &cell_structure = *system.cell_structure;
  auto const &aosoa = cell_structure.get_aosoa();
#ifdef ESPRESSO_EXCLUSIONS
  // The specialized kernel has no exclusion handling. The commit sweep (run by
  // update_verlet_state earlier in this same force call) accumulates whether
  // any packed particle carries an exclusion, so this is an O(1) read of the
  // same population the old per-particle sweep covered (local + ghosts).
  if (aosoa.has_any_exclusion())
    return {};
#endif

  auto scatter_force = cell_structure.get_scatter_force();
  auto const minimum_image = system.box_geo->cuboid_minimum_image();
  auto const max_cutoff_sq = Utils::sqr(system.maximal_cutoff());

#ifdef ESPRESSO_ELECTROSTATICS
  if (auto const *coulomb_ptr = get_ptr(coulomb_kernel);
      coulomb_ptr != nullptr) {
    return make_specialized_verlet_pair_loop<true>(
        nonbonded_ias, aosoa, scatter_force, minimum_image, max_cutoff_sq,
        coulomb_ptr
#ifdef ESPRESSO_P3M
        ,
        get_toplevel_p3m_solver(system.coulomb)
#endif
    );
  }
#else
  static_cast<void>(coulomb_kernel);
#endif

  return make_specialized_verlet_pair_loop<false>(
      nonbonded_ias, aosoa, scatter_force, minimum_image, max_cutoff_sq);
}

static void reduce_cabana_forces_and_torques(System::System const &system,
                                             Utils::Vector3d *virial) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif

  auto const &unique_particles = system.cell_structure->get_unique_particles();
  auto &local_force = system.cell_structure->get_local_force();
  auto scatter_force = system.cell_structure->get_scatter_force();
  Kokkos::Experimental::contribute(local_force, scatter_force);
#ifdef ESPRESSO_ROTATION
  // when no kernel scattered into the torque buffers this pass, they are
  // all-zero and both the reduction over replicas and the per-particle
  // accumulation would only add zeros
  auto const reduce_torque = system.cell_structure->torque_replicas_dirty();
  auto &local_torque = system.cell_structure->get_local_torque();
  if (reduce_torque) {
    auto scatter_torque = system.cell_structure->get_scatter_torque();
    Kokkos::Experimental::contribute(local_torque, scatter_torque);
  }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  auto const reduce_dip_fld = system.cell_structure->dip_fld_replicas_dirty();
  auto &local_dip_fld = system.cell_structure->get_local_dip_fld();
  if (reduce_dip_fld) {
    auto scatter_dip_fld = system.cell_structure->get_scatter_dip_fld();
    Kokkos::Experimental::contribute(local_dip_fld, scatter_dip_fld);
  }
#endif
#ifdef ESPRESSO_NPT
  auto &local_virial = system.cell_structure->get_local_virial();
  if (system.cell_structure->virial_replicas_dirty()) {
    auto scatter_virial = system.cell_structure->get_scatter_virial();
    Kokkos::Experimental::contribute(local_virial, scatter_virial);
  }
#endif

  using execution_space = Kokkos::DefaultHostExecutionSpace;
  Kokkos::RangePolicy<execution_space> policy(std::size_t{0},
                                              unique_particles.size());
  Kokkos::parallel_for("reduction", policy,
                       [&local_force,
#ifdef ESPRESSO_ROTATION
                        &local_torque, reduce_torque,
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
                        &local_dip_fld, reduce_dip_fld,
#endif
                        &unique_particles](std::size_t const i) {
                         Utils::Vector3d force{};
                         force[0] += local_force(i, 0);
                         force[1] += local_force(i, 1);
                         force[2] += local_force(i, 2);
                         unique_particles.at(i)->force() += force;
#ifdef ESPRESSO_ROTATION
                         if (reduce_torque) {
                           Utils::Vector3d torque{};
                           torque[0] += local_torque(i, 0);
                           torque[1] += local_torque(i, 1);
                           torque[2] += local_torque(i, 2);
                           unique_particles.at(i)->torque() += torque;
                         }
#endif // ESPRESSO_ROTATION
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
                         if (reduce_dip_fld) {
                           Utils::Vector3d dip_fld{};
                           dip_fld[0] += local_dip_fld(i, 0);
                           dip_fld[1] += local_dip_fld(i, 1);
                           dip_fld[2] += local_dip_fld(i, 2);
                           unique_particles.at(i)->dip_fld() += dip_fld;
                         }
#endif // ESPRESSO_DIPOLE_FIELD_TRACKING
                       });
  Kokkos::fence();

#ifdef ESPRESSO_NPT
  if (virial) {
    (*virial)[0] += local_virial(0);
    (*virial)[1] += local_virial(1);
    (*virial)[2] += local_virial(2);
  }
#endif
}

void System::System::calculate_forces() {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
#ifdef ESPRESSO_CUDA
  {
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("copy_particles_to_GPU");
#endif
    gpu->update();
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("copy_particles_to_GPU");
#endif
  }
#endif // ESPRESSO_CUDA

#ifdef ESPRESSO_COLLISION_DETECTION
  collision_detection->clear_queue();
  auto const collision_detection_cutoff = collision_detection->cutoff();
#else
  auto const collision_detection_cutoff = inactive_cutoff;
#endif
  bond_breakage->clear_queue();
  auto particles = cell_structure->local_particles();
#ifdef ESPRESSO_NPT
  if (propagation->used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) {
    // reset virial part of instantaneous pressure
    npt_inst_pressure->p_vir = Utils::Vector3d{};
  }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (dipoles.impl->solver.has_value()) {
    reinit_dip_fld(*cell_structure);
  }
#endif

  // Use combined function instead of two separate calls

  auto const elc_kernel = coulomb.pair_force_elc_kernel();
  auto const coulomb_kernel = coulomb.pair_force_kernel();
  auto const dipoles_kernel = dipoles.pair_force_kernel();
  auto const coulomb_u_kernel = coulomb.pair_energy_kernel();
  auto *const virial = get_npt_virial();

  // Factory instead of an eager criterion: construction fills an O(n_types^2)
  // cutoff table, so it only runs where a criterion is actually consumed (the
  // link-cell fallback and the collision-detection loop below).
  auto const make_verlet_criterion = [&] {
    return VerletCriterion<>{*this,
                             cell_structure->get_verlet_skin(),
                             get_interaction_range(),
                             coulomb.cutoff(),
                             dipoles.cutoff(),
                             collision_detection_cutoff};
  };

  update_verlet_state(*this, collision_detection_cutoff);
#ifdef ESPRESSO_ELECTROSTATICS
  if (coulomb.impl->extension) {
    update_icc_particles();
  }
#endif // ESPRESSO_ELECTROSTATICS
  init_forces_and_thermostat(*this);
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_BEGIN("calc_long_range_forces");
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.calc_long_range_force();
#endif
#ifdef ESPRESSO_DIPOLES
  dipoles.calc_long_range_force();
#endif
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_END("calc_long_range_forces");
#endif

#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_BEGIN("cabana_short_range");
#endif
  auto &bs = cell_structure->bond_state();
  auto bonds_kernel_data = create_kokkos_bonds_kernel_data(*this);
  auto pair_bonds_kernel = PairBondsKernel{
      bonds_kernel_data, bs.pair_list, bs.pair_ids, get_ptr(coulomb_kernel)};
  auto angle_bonds_kernel =
      AngleBondsKernel{bonds_kernel_data, bs.angle_list, bs.angle_ids};
  auto dihedral_bonds_kernel =
      DihedralBondsKernel{bonds_kernel_data, bs.dihedral_list, bs.dihedral_ids};

  auto first_neighbor_kernel =
      create_cabana_neighbor_kernel(*this, virial, elc_kernel, coulomb_kernel,
                                    dipoles_kernel, coulomb_u_kernel);

  auto const specialized_pair_loop = create_specialized_verlet_pair_loop(
      *this, virial, elc_kernel, coulomb_kernel, dipoles_kernel);

  // The specialized pair kernel scatters only into the force view. Record
  // before the launch whether this pass can write the torque/virial scatter
  // buffers at all, so their O(n_threads * N) zeroing and reduction can be
  // skipped otherwise (see CellStructure::reset_local_properties). All other
  // torque sources (Langevin rotation, virtual sites back-transfer, the
  // dipolar solvers except dp3m, constraints) write Particle::torque()
  // directly and never touch the scatter buffers; dp3m marks the flag at its
  // own scatter site.
  [[maybe_unused]] auto const generic_pair_path =
      get_interaction_range() > 0. and
      (not specialized_pair_loop or
       propagation->integ_switch == INTEG_METHOD_STEEPEST_DESCENT or
       not cell_structure->use_verlet_list);
#ifdef ESPRESSO_ROTATION
  // Within the generic kernel, only orientation-dependent pair potentials
  // (Gay-Berne) and the dipolar pair kernel scatter into the torque view.
  auto const gay_berne_active =
      nonbonded_ias->pair_potential_active(PairPotential::GayBerne);
  auto const dipolar_pair_kernel_active =
#ifdef ESPRESSO_DIPOLES
      get_ptr(dipoles_kernel) != nullptr;
#else
      false;
#endif
  if (generic_pair_path and (gay_berne_active or dipolar_pair_kernel_active)) {
    cell_structure->mark_torque_replicas_dirty();
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    cell_structure->mark_dip_fld_replicas_dirty();
#endif
    // Pair-kernel torques also land on ghost particles and only get home via
    // the TORQUE ghost reduce. Each torque-scattering kernel therefore needs
    // a matching arm in orientation_ghosts_needed() (System.cpp).
    assert(get_force_reduce_ghost_flags() & GHOSTTRANS_TORQUE);
  }
#endif // ESPRESSO_ROTATION
#ifdef ESPRESSO_NPT
  auto const have_bonds = cell_structure->get_local_pair_bond_numbers() > 0 or
                          cell_structure->get_local_angle_bond_numbers() > 0 or
                          cell_structure->get_local_dihedral_bond_numbers() > 0;
  if (generic_pair_path or
      (bonded_ias->maximal_cutoff() >= 0. and have_bonds)) {
    cell_structure->mark_virial_replicas_dirty();
  }
#endif // ESPRESSO_NPT

  cabana_short_range(pair_bonds_kernel, angle_bonds_kernel,
                     dihedral_bonds_kernel, first_neighbor_kernel,
                     *cell_structure, get_interaction_range(),
                     bonded_ias->maximal_cutoff(), make_verlet_criterion,
                     propagation->integ_switch, specialized_pair_loop);

  // Force and Torque reduction
  reduce_cabana_forces_and_torques(*this, virial);

#ifdef ESPRESSO_COLLISION_DETECTION
  auto collision_kernel = [&collision_detection = *collision_detection](
                              Particle const &p1, Particle const &p2,
                              Distance const &d) {
    collision_detection.detect_collision(p1, p2, d.dist2);
  };
  if (not collision_detection->is_off()) {
    auto const verlet_criterion = make_verlet_criterion();
    cell_structure->non_bonded_loop(collision_kernel, verlet_criterion);
  }
#endif // ESPRESSO_COLLISION_DETECTION

#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_END("cabana_short_range");
#endif

  constraints->add_forces(particles, get_sim_time());
  oif_global->calculate_forces();

  // Must be done here. Forces need to be ghost-communicated
  immersed_boundaries->volume_conservation(*cell_structure);

  if (thermostat->lb and (propagation->used_propagations &
                          PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE)) {
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("lb_particle_coupling");
#endif
    lb_couple_particles();
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("lb_particle_coupling");
#endif
  }

#ifdef ESPRESSO_CUDA
  {
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("copy_forces_from_GPU");
#endif
    gpu->copy_forces_to_host(particles, this_node);

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    gpu->copy_dip_fld_to_host(particles, this_node);
#endif

#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("copy_forces_from_GPU");
#endif
  }
#endif // ESPRESSO_CUDA

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  if (propagation->used_propagations &
      (PropagationMode::TRANS_VS_RELATIVE | PropagationMode::ROT_VS_RELATIVE |
       PropagationMode::ROT_VS_INDEPENDENT)) {
    vs_relative_back_transfer_forces_and_torques(*cell_structure);
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
  if (propagation->used_propagations &
      (PropagationMode::TRANS_VS_CENTER_OF_MASS)) {
    vs_com_back_transfer_forces_and_torques(*cell_structure);
  }
#endif

  if (ghost_reduce_overlap_eligible(*this)) {
    // Split-phase path: start the ghost force reduction now (non-blocking);
    // the integrator will run interior-cell step_2, finish the reduce, then
    // run boundary-cell step_2.  comfixed and force_capping are inactive on
    // this path.
    assert(comfixed->get_fixed_types().empty() &&
           "ghost_reduce_overlap: comfixed must be inactive on the eligible "
           "path");
    assert(force_cap == 0. &&
           "ghost_reduce_overlap: force_cap must be 0 on the eligible path");
    cell_structure->ghosts_reduce_forces_start();
  } else {
    // Blocking path: finish the reduce here, then apply comfixed/capping.
    cell_structure->ghosts_reduce_forces();
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    if (dipoles.impl->solver.has_value()) {
      cell_structure->ghosts_reduce_dipole_field();
    }
#endif

    // should be pretty late, since it needs to zero out the total force
    comfixed->apply(particles);

    // Needs to be the last one to be effective
    force_capping(*cell_structure, force_cap);
  }

  // mark that forces are now up-to-date
  propagation->recalc_forces = false;
}
