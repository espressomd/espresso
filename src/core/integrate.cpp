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

/** \file
 *  Molecular dynamics integrator.
 *
 *  For more information about the integrator
 *  see \ref integrate.hpp "integrate.hpp".
 */

#include "integrate.hpp"
#include "integrators/Propagation.hpp"
#include "integrators/brownian_inline.hpp"
#include "integrators/steepest_descent.hpp"
#include "integrators/stokesian_dynamics_inline.hpp"
#include "integrators/symplectic_euler_inline.hpp"
#include "integrators/velocity_verlet_inline.hpp"
#include "integrators/velocity_verlet_npt.hpp"

#include "BoxGeometry.hpp"
#include "PropagationMode.hpp"
#include "accumulators/AutoUpdateAccumulators.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/for_each_particle.hpp"
#include "cells.hpp"
#include "collision_detection/CollisionDetection.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "lb/particle_coupling.hpp"
#include "lb/utils.hpp"
#include "lees_edwards/lees_edwards.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "rattle.hpp"
#include "rotation.hpp"
#include "signalhandling.hpp"
#include "stokesian_dynamics/sd_interface.hpp"
#include "system/System.hpp"
#include "system/System.impl.hpp"
#include "thermostat.hpp"
#include "thermostats/langevin_inline.hpp"
#include "virtual_sites/com.hpp"
#include "virtual_sites/lb_tracers.hpp"
#include "virtual_sites/relative.hpp"

#include <instrumentation/fe_trap.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#ifdef ESPRESSO_CALIPER
#include "caliper_utils.hpp"
#endif

#ifdef ESPRESSO_VALGRIND
#include <callgrind.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef ESPRESSO_WALBERLA
#ifdef ESPRESSO_WALBERLA_STATIC_ASSERT
#error "waLberla headers should not be visible to the ESPResSo core"
#endif
#endif

namespace {
volatile std::sig_atomic_t ctrl_C = 0;
} // namespace

namespace LeesEdwards {

/**
 * @brief Update the Lees-Edwards parameters of the box geometry
 * for the current simulation time.
 */
void LeesEdwards::update_box_params(BoxGeometry &box_geo, double sim_time) {
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    assert(m_protocol != nullptr);
    box_geo.lees_edwards_update(get_pos_offset(sim_time, *m_protocol),
                                get_shear_velocity(sim_time, *m_protocol));
  }
}

void LeesEdwards::set_protocol(std::shared_ptr<ActiveProtocol> protocol) {
  auto &system = get_system();
  auto &cell_structure = *system.cell_structure;
  auto &box_geo = *system.box_geo;
  box_geo.set_type(BoxType::LEES_EDWARDS);
  m_protocol = std::move(protocol);
  update_box_params(box_geo, system.get_sim_time());
  system.propagation->recalc_forces = true;
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
}

void LeesEdwards::unset_protocol() {
  auto &system = get_system();
  auto &cell_structure = *system.cell_structure;
  auto &box_geo = *system.box_geo;
  m_protocol = nullptr;
  box_geo.set_type(BoxType::CUBOID);
  system.propagation->recalc_forces = true;
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
}

} // namespace LeesEdwards

void Propagation::update_default_propagation(int thermo_switch) {
  switch (integ_switch) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    default_propagation = PropagationMode::NONE;
    break;
  case INTEG_METHOD_NVT:
  case INTEG_METHOD_SYMPLECTIC_EULER: {
    // NOLINTNEXTLINE(bugprone-branch-clone)
    if ((thermo_switch & THERMO_LB) and (thermo_switch & THERMO_LANGEVIN)) {
      default_propagation = PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE;
#ifdef ESPRESSO_ROTATION
      default_propagation |= PropagationMode::ROT_LANGEVIN;
#endif
    } else if (thermo_switch & THERMO_LB) {
      default_propagation = PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE;
#ifdef ESPRESSO_ROTATION
      default_propagation |= PropagationMode::ROT_EULER;
#endif
    } else if (thermo_switch & THERMO_LANGEVIN) {
      default_propagation = PropagationMode::TRANS_LANGEVIN;
#ifdef ESPRESSO_ROTATION
      default_propagation |= PropagationMode::ROT_LANGEVIN;
#endif
    } else {
      default_propagation = PropagationMode::TRANS_NEWTON;
#ifdef ESPRESSO_ROTATION
      default_propagation |= PropagationMode::ROT_EULER;
#endif
    }
    break;
  }
#ifdef ESPRESSO_NPT
  case INTEG_METHOD_NPT_ISO_AND:
  case INTEG_METHOD_NPT_ISO_MTK:
    default_propagation = PropagationMode::TRANS_LANGEVIN_NPT;
    break;
#endif
  case INTEG_METHOD_BD:
    default_propagation = PropagationMode::TRANS_BROWNIAN;
#ifdef ESPRESSO_ROTATION
    default_propagation |= PropagationMode::ROT_BROWNIAN;
#endif
    break;
#ifdef ESPRESSO_STOKESIAN_DYNAMICS
  case INTEG_METHOD_SD:
    default_propagation = PropagationMode::TRANS_STOKESIAN;
    break;
#endif // ESPRESSO_STOKESIAN_DYNAMICS
  default:
    throw std::runtime_error("Unknown value for integ_switch");
  }
}

void System::System::update_used_propagations() {
  int used_propagations = PropagationMode::NONE;
  for (auto &p : cell_structure->local_particles()) {
    used_propagations |= p.propagation();
  }
  if (used_propagations & PropagationMode::SYSTEM_DEFAULT) {
    used_propagations |= propagation->default_propagation;
  }
  used_propagations = boost::mpi::all_reduce(::comm_cart, used_propagations,
                                             std::bit_or<int>());
  propagation->used_propagations = used_propagations;
  propagation->recalc_used_propagations = false;
}

void System::System::integrator_sanity_checks() const {
  auto const thermo_switch = thermostat->thermo_switch;
  if (time_step <= 0.) {
    runtimeErrorMsg() << "time_step not set";
  }
  if (propagation->integ_switch == INTEG_METHOD_STEEPEST_DESCENT) {
    if (thermo_switch != THERMO_OFF) {
      runtimeErrorMsg()
          << "The steepest descent integrator is incompatible with thermostats";
    }
  }
  if (propagation->integ_switch == INTEG_METHOD_NVT) {
    if (thermo_switch & (THERMO_NPT_ISO | THERMO_BROWNIAN | THERMO_SD)) {
      runtimeErrorMsg() << "The VV integrator is incompatible with the "
                           "currently active combination of thermostats";
    }
  }
#ifdef ESPRESSO_NPT
  if (propagation->used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) {
    if (thermo_switch != THERMO_NPT_ISO) {
      runtimeErrorMsg() << "The NpT integrator requires the NpT thermostat";
    }
    if (box_geo->type() == BoxType::LEES_EDWARDS) {
      runtimeErrorMsg() << "The NpT integrator cannot use Lees-Edwards";
    }
    try {
      nptiso->coulomb_dipole_sanity_checks(*this);
    } catch (std::runtime_error const &err) {
      runtimeErrorMsg() << err.what();
    }
  }
#endif
  if (propagation->used_propagations & PropagationMode::TRANS_BROWNIAN) {
    if (thermo_switch != THERMO_BROWNIAN) {
      runtimeErrorMsg() << "The BD integrator requires the BD thermostat";
    }
  }
  if (propagation->used_propagations & PropagationMode::TRANS_STOKESIAN) {
#ifdef ESPRESSO_STOKESIAN_DYNAMICS
    if (thermo_switch != THERMO_SD) {
      runtimeErrorMsg() << "The SD integrator requires the SD thermostat";
    }
#endif
  }
  if (lb.is_solver_set() and (propagation->used_propagations &
                              (PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE |
                               PropagationMode::TRANS_LB_TRACER))) {
    if (thermostat->lb == nullptr) {
      runtimeErrorMsg() << "The LB integrator requires the LB thermostat";
    }
  }
  if (bonded_ias->get_n_thermalized_bonds() >= 1 and
      (thermostat->thermalized_bond == nullptr or
       (thermo_switch & THERMO_BOND) == 0)) {
    runtimeErrorMsg()
        << "Thermalized bonds require the thermalized_bond thermostat";
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (bonded_ias->get_n_rigid_bonds() >= 1) {
    if (not propagation->is_inertial()) {
      runtimeErrorMsg()
          << "Rigid bonds (RATTLE) require an inertial integrator "
             "(VV or symplectic Euler); BD and SD are not supported";
    }
  }
#endif

#ifdef ESPRESSO_STOKESIAN_DYNAMICS
  if ((propagation->used_propagations & PropagationMode::TRANS_STOKESIAN) and
      (propagation->default_propagation & PropagationMode::TRANS_STOKESIAN)) {
    auto pred = PropagationPredicateStokesian(propagation->default_propagation);
    stokesian_dynamics->sanity_checks(
        cell_structure->local_particles().filter(pred));
  }
#endif // ESPRESSO_STOKESIAN_DYNAMICS

#ifdef ESPRESSO_ROTATION
  for (auto const &p : cell_structure->local_particles()) {
    using namespace PropagationMode;
    if (p.can_rotate() and not p.is_virtual() and
        (p.propagation() & (SYSTEM_DEFAULT | ROT_EULER | ROT_LANGEVIN |
                            ROT_BROWNIAN | ROT_STOKESIAN)) == 0) {
      runtimeErrorMsg()
          << "Rotating particles must have a rotation propagation mode enabled";
      break;
    }
  }
#endif // ESPRESSO_ROTATION

#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
#ifdef ESPRESSO_EXTERNAL_FORCES
  if (propagation->used_propagations &
      PropagationMode::TRANS_VS_CENTER_OF_MASS) {
    for (auto const &p : cell_structure->local_particles()) {
      using namespace PropagationMode;
      if ((p.propagation() & TRANS_VS_CENTER_OF_MASS) and
          p.has_fixed_coordinates()) {
        runtimeErrorMsg() << "VS COM particles cannot be fixed in space";
        break;
      }
    }
  }
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (bonded_ias->get_n_rigid_bonds()) {
    using namespace PropagationMode;
    for (auto const &p : cell_structure->local_particles()) {
      if (p.propagation() & TRANS_VS_CENTER_OF_MASS) {
        for (auto const bond : p.bonds()) {
          if (std::holds_alternative<RigidBond>(
                  *bonded_ias->at(bond.bond_id()))) {
            runtimeErrorMsg() << "VS COM particles cannot use rigid bonds";
            break;
          }
        }
      }
    }
  }
#endif // ESPRESSO_BOND_CONSTRAINT
#endif // ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  if ((thermo_switch & THERMO_LANGEVIN) == 0) {
    for (auto const &p : cell_structure->local_particles()) {
      if (p.stoner_wohlfarth_is_enabled()) {
        runtimeErrorMsg() << "The thermal Stoner-Wohlfarth model requires the "
                             "Langevin thermostat";
        break;
      }
    }
  }
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH
}

#ifdef ESPRESSO_WALBERLA
void walberla_tau_sanity_checks(std::string const &method, double tau,
                                double time_step) {
  if (time_step <= 0.) {
    return;
  }
  // use float epsilon since tau may be a float
  auto const eps = static_cast<double>(std::numeric_limits<float>::epsilon());
  if ((tau - time_step) / (tau + time_step) < -eps)
    throw std::invalid_argument(method + " tau (" + std::to_string(tau) +
                                ") must be >= MD time_step (" +
                                std::to_string(time_step) + ")");
  auto const factor = tau / time_step;
  if (std::fabs(std::round(factor) - factor) / factor > eps)
    throw std::invalid_argument(method + " tau (" + std::to_string(tau) +
                                ") must be an integer multiple of the "
                                "MD time_step (" +
                                std::to_string(time_step) + "). Factor is " +
                                std::to_string(factor));
}

void walberla_agrid_sanity_checks(std::string const &method,
                                  Utils::Vector3d const &geo_left,
                                  Utils::Vector3d const &geo_right,
                                  Utils::Vector3d const &lattice_left,
                                  Utils::Vector3d const &lattice_right,
                                  double agrid) {
  // waLBerla and ESPResSo must agree on domain decomposition
  auto const tol = agrid / 1E6;
  if ((lattice_left - geo_left).norm2() > tol or
      (lattice_right - geo_right).norm2() > tol) {
    std::stringstream error_msg;
    error_msg << "waLBerla and ESPResSo disagree about domain decomposition"
              << "\nMPI rank " << ::this_node << ": "
              << "left ESPResSo: [" << geo_left << "], "
              << "left waLBerla: [" << lattice_left << "]"
              << "\nMPI rank " << ::this_node << ": "
              << "right ESPResSo: [" << geo_right << "], "
              << "right waLBerla: [" << lattice_right << "]"
              << "\nfor method: " << method;
    throw std::runtime_error(error_msg.str());
  }
}
#endif // ESPRESSO_WALBERLA

static void resort_particles_if_needed(System::System &system) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  auto &cell_structure = *system.cell_structure;
  auto const offset = LeesEdwards::verlet_list_offset(
      *system.box_geo, cell_structure.get_le_pos_offset_at_last_resort());
  if (cell_structure.check_resort_required(offset)) {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

/** @brief Calls the hook for propagation kernels before the force calculation
 *  @return whether or not to stop the integration loop early.
 */
static bool integrator_step_1(CellStructure &cell_structure,
                              Propagation const &propagation,
                              System::System &system, double time_step) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  // steepest decent
  if (propagation.integ_switch == INTEG_METHOD_STEEPEST_DESCENT)
    return system.steepest_descent->propagate(cell_structure);

  auto const &thermostat = *system.thermostat;
  auto const kT = thermostat.kT;
  // Hoist the velocity-Verlet translation column-view handles ONCE outside the
  // parallel_for (the handle copy here is fine; per-element view rebinding is
  // not). The still-view-path operations (symplectic Euler, rotation,
  // Brownian) rebind a Particle lazily inside their branches only.
  auto &store = cell_structure.particle_store();
  auto vel_view = store.velocity_view();
  auto pos_view = store.position_view();
  auto force_view = store.force_view();
  auto id_view = store.id_view();
  // Disabled-feature handles are typed zero-extent dummy views (correct column
  // type, extent 0): they exist only to fix the shared kernel lambda signature
  // and are never indexed on a compiled-out path. See ParticleStore::dummy_*.
#ifdef ESPRESSO_MASS
  auto mass_view = store.mass_view();
#else
  auto mass_view = store.dummy_scalar_view();
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  auto ext_flag_view = store.ext_flag_view();
#else
  auto ext_flag_view = store.dummy_uint8_view();
#endif
#ifdef ESPRESSO_ROTATION
  auto quat_view = store.quaternion_view();
  auto omega_view = store.angular_velocity_view();
  auto torque_view = store.torque_view();
  auto rotation_view = store.rotation_view();
#else
  auto quat_view = store.dummy_quaternion_view();
  auto omega_view = store.dummy_vector_view();
  auto torque_view = store.dummy_vector_view();
  auto rotation_view = store.dummy_uint8_view();
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  auto rinertia_view = store.rinertia_view();
#else
  auto rinertia_view = store.dummy_vector_view();
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto gamma_view = store.gamma_view();
#else
  auto gamma_view = store.dummy_scalar_view();
#endif
#if defined(ESPRESSO_THERMOSTAT_PER_PARTICLE) && defined(ESPRESSO_ROTATION)
  auto gamma_rot_view = store.gamma_rot_view();
#else
  auto gamma_rot_view = store.dummy_scalar_view();
#endif
  cell_structure.for_each_local_particle_row([&](int const row) {
    Particle p;
    p.attach_to_store(store, row);
#ifdef ESPRESSO_VIRTUAL_SITES
    // virtual sites are updated later in the integration loop
    if (p.is_virtual())
      return;
#endif
    // Read the propagation bitfield ONCE per particle; every mode query below
    // reuses it instead of re-reading the ParticleStore propagation column.
    int const prop = p.propagation();
    if (propagation.integ_switch == INTEG_METHOD_SYMPLECTIC_EULER) {
      if (propagation.should_propagate_with(
              prop, PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE))
        symplectic_euler_propagator_1(p, time_step);
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_NEWTON))
        symplectic_euler_propagator_1(p, time_step);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop, PropagationMode::ROT_EULER))
        symplectic_euler_rotator_1(p, time_step);
#endif
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_LANGEVIN))
        symplectic_euler_propagator_1(p, time_step);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::ROT_LANGEVIN))
        symplectic_euler_rotator_1(p, time_step);
#endif
    } else {
      // Fixed-coordinate bitfield / mass read ONCE per particle from the
      // hoisted views (compile-time fallbacks when the feature is off, matching
      // Particle::fixed_flags_byte() / Particle::mass()).
#ifdef ESPRESSO_MASS
      double const mass = mass_view(row);
#else
      double const mass = 1.0;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
      std::uint8_t const fixed = ext_flag_view(row);
#else
      auto const fixed = static_cast<std::uint8_t>(0u);
#endif
      // Resolve the momentum row bases ONCE (one pointer + stride each); the
      // VV propagator then does stride-1 pointer arithmetic per axis. Cheaper
      // than re-subscripting a 2D column view per component (perf-iterate).
      auto vel = store.velocity_reference(row);
      auto pos = store.position_reference(row);
      auto const force = store.force_reference(row);
      if (propagation.should_propagate_with(
              prop, PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE))
        velocity_verlet_propagator_1(vel, pos, force, mass, fixed, time_step);
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_NEWTON))
        velocity_verlet_propagator_1(vel, pos, force, mass, fixed, time_step);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop, PropagationMode::ROT_EULER))
        velocity_verlet_rotator_1(quat_view, omega_view, rinertia_view,
                                  torque_view, rotation_view, row, time_step);
#endif
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_LANGEVIN))
        velocity_verlet_propagator_1(vel, pos, force, mass, fixed, time_step);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::ROT_LANGEVIN))
        velocity_verlet_rotator_1(quat_view, omega_view, rinertia_view,
                                  torque_view, rotation_view, row, time_step);
#endif
    }
    if (propagation.should_propagate_with(prop,
                                          PropagationMode::TRANS_BROWNIAN))
      brownian_dynamics_propagator(
          *thermostat.brownian,
          make_brownian_row_view(pos_view, vel_view, force_view, torque_view,
                                 quat_view, omega_view, rinertia_view, id_view,
                                 mass_view, rotation_view, ext_flag_view,
                                 gamma_view, gamma_rot_view, row),
          time_step, kT);
#ifdef ESPRESSO_ROTATION
    if (propagation.should_propagate_with(prop, PropagationMode::ROT_BROWNIAN))
      brownian_dynamics_rotator(
          *thermostat.brownian,
          make_brownian_row_view(pos_view, vel_view, force_view, torque_view,
                                 quat_view, omega_view, rinertia_view, id_view,
                                 mass_view, rotation_view, ext_flag_view,
                                 gamma_view, gamma_rot_view, row),
          time_step, kT);
#endif
  });

#ifdef ESPRESSO_NPT
  if ((propagation.used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) and
      (propagation.default_propagation & PropagationMode::TRANS_LANGEVIN_NPT)) {
    auto pred = PropagationPredicateNPT(propagation.default_propagation);
    if (propagation.integ_switch == INTEG_METHOD_NPT_ISO_AND) {
      velocity_verlet_npt_Andersen_step_1(
          cell_structure.local_particles().filter(pred), *thermostat.npt_iso,
          time_step, system);
    } else if (propagation.integ_switch == INTEG_METHOD_NPT_ISO_MTK) {
      velocity_verlet_npt_MTK_step_1(
          cell_structure.local_particles().filter(pred), *thermostat.npt_iso,
          time_step, system);
    }
  }
#endif

#ifdef ESPRESSO_STOKESIAN_DYNAMICS
  if ((propagation.used_propagations & PropagationMode::TRANS_STOKESIAN) and
      (propagation.default_propagation & PropagationMode::TRANS_STOKESIAN)) {
    auto pred = PropagationPredicateStokesian(propagation.default_propagation);
    stokesian_dynamics_step_1(cell_structure.local_particles().filter(pred),
                              *system.stokesian_dynamics, *thermostat.stokesian,
                              time_step, kT);
  }
#endif // ESPRESSO_STOKESIAN_DYNAMICS

  return false;
}

/**
 * @brief Build the per-particle half-kick callable for step_2.
 *
 * Returns a lambda that applies the velocity (and torque, if ROTATION is
 * enabled) update for a single particle. Virtual sites are skipped.
 *
 * Shared verbatim by @ref integrator_step_2 (full pass) and
 * @ref integrator_step_2_filtered (interior / boundary passes): the lambda is
 * constructed once, then handed to @c for_each_local_particle,
 * @c for_each_interior_particle, or @c for_each_boundary_particle.
 *
 * The per-particle scalar/rotation column-view handles are resolved ONCE here
 * (see integrator_step_1) and captured by the returned lambda, which then
 * indexes them by the particle's store row. The VV translation velocity/force
 * rows are resolved per row via velocity_reference/force_reference
 * (VectorReference), so no 2D velocity/force view handle is hoisted.
 *
 * NPT particles are intentionally absent: the NPT arm must run only on the
 * ineligible (full-reduce-then-step_2) path and is handled separately inside
 * @ref integrator_step_2.
 */
static auto make_step2_particle_kernel(CellStructure &cell_structure,
                                       Propagation const &propagation,
                                       double time_step) {
  auto &store = cell_structure.particle_store();
#ifdef ESPRESSO_MASS
  auto mass_view = store.mass_view();
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  auto ext_flag_view = store.ext_flag_view();
#endif
#ifdef ESPRESSO_ROTATION
  auto quat_view = store.quaternion_view();
  auto omega_view = store.angular_velocity_view();
  auto torque_view = store.torque_view();
  auto rotation_view = store.rotation_view();
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  auto rinertia_view = store.rinertia_view();
#else
  // Correct-typed zero-extent dummy: rinertia_view is never indexed when
  // rotational inertia is off (the rotator kernels use {1,1,1}); see
  // ParticleStore::dummy_vector_view.
  auto rinertia_view = store.dummy_vector_view();
#endif
#endif
  return [&propagation, time_step, &store
#ifdef ESPRESSO_MASS
          ,
          mass_view
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
          ,
          ext_flag_view
#endif
#ifdef ESPRESSO_ROTATION
          ,
          quat_view, omega_view, torque_view, rotation_view, rinertia_view
#endif
  ](Particle &p) {
    auto const row = p.store_row();
#ifdef ESPRESSO_VIRTUAL_SITES
    // virtual sites are updated later in the integration loop
    if (p.is_virtual())
      return;
#endif
    // Read the propagation bitfield ONCE per particle (see integrator_step_1).
    int const prop = p.propagation();
    if (propagation.integ_switch == INTEG_METHOD_SYMPLECTIC_EULER) {
      if (propagation.should_propagate_with(
              prop, PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE))
        symplectic_euler_propagator_2(p, time_step);
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_NEWTON))
        symplectic_euler_propagator_2(p, time_step);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop, PropagationMode::ROT_EULER))
        symplectic_euler_rotator_2(p, time_step);
#endif
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_LANGEVIN))
        symplectic_euler_propagator_2(p, time_step);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::ROT_LANGEVIN))
        symplectic_euler_rotator_2(p, time_step);
#endif
    } else {
#ifdef ESPRESSO_MASS
      double const mass = mass_view(row);
#else
      double const mass = 1.0;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
      std::uint8_t const fixed = ext_flag_view(row);
#else
      auto const fixed = static_cast<std::uint8_t>(0u);
#endif
      // Resolve the momentum row bases ONCE (see integrator_step_1).
      auto vel = store.velocity_reference(row);
      auto const force = store.force_reference(row);
      if (propagation.should_propagate_with(
              prop, PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE))
        velocity_verlet_propagator_2(vel, force, mass, fixed, time_step);
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_NEWTON))
        velocity_verlet_propagator_2(vel, force, mass, fixed, time_step);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop, PropagationMode::ROT_EULER))
        velocity_verlet_rotator_2(quat_view, omega_view, rinertia_view,
                                  torque_view, rotation_view, row, time_step);
#endif
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::TRANS_LANGEVIN))
        velocity_verlet_propagator_2(vel, force, mass, fixed, time_step);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(prop,
                                            PropagationMode::ROT_LANGEVIN))
        velocity_verlet_rotator_2(quat_view, omega_view, rinertia_view,
                                  torque_view, rotation_view, row, time_step);
#endif
    }
  };
}

static void integrator_step_2(CellStructure &cell_structure,
                              Propagation const &propagation,
                              [[maybe_unused]] System::System &system,
                              double time_step) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  if (propagation.integ_switch == INTEG_METHOD_STEEPEST_DESCENT)
    return;

  cell_structure.for_each_local_particle(
      make_step2_particle_kernel(cell_structure, propagation, time_step));

#ifdef ESPRESSO_NPT
  if ((propagation.used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) and
      (propagation.default_propagation & PropagationMode::TRANS_LANGEVIN_NPT)) {
    auto pred = PropagationPredicateNPT(propagation.default_propagation);
    if (propagation.integ_switch == INTEG_METHOD_NPT_ISO_AND) {
      velocity_verlet_npt_Andersen_step_2(
          cell_structure.local_particles().filter(pred), time_step, system);
    } else if (propagation.integ_switch == INTEG_METHOD_NPT_ISO_MTK) {
      velocity_verlet_npt_MTK_step_2(
          cell_structure.local_particles().filter(pred), time_step, system);
    }
  }
#endif
}

/**
 * @brief Restricted step_2 for the split-phase ghost-force-reduce overlap.
 *
 * Called twice per step: once with @p interior_pass = true (while the ghost
 * reduce is in flight) and once with @p interior_pass = false (after
 * @ref CellStructure::ghosts_reduce_forces_finish completes).  The per-particle
 * kernel is identical to the one used by @ref integrator_step_2 — shared via
 * @ref make_step2_particle_kernel.
 *
 * NPT is absent here: the eligibility check in calculate_forces() guarantees
 * TRANS_LANGEVIN_NPT is never active when this path runs.
 */
static void integrator_step_2_filtered(CellStructure &cell_structure,
                                       Propagation const &propagation,
                                       double time_step, bool interior_pass) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  // steepest_descent and NPT are excluded by the eligibility check; both
  // asserts are belt-and-suspenders cross-checks.
  assert(propagation.integ_switch != INTEG_METHOD_STEEPEST_DESCENT &&
         "integrator_step_2_filtered: steepest-descent is ineligible");
#ifdef ESPRESSO_NPT
  assert((propagation.used_propagations &
          PropagationMode::TRANS_LANGEVIN_NPT) == 0 &&
         "integrator_step_2_filtered: NPT propagation is ineligible");
#endif

  auto const kernel =
      make_step2_particle_kernel(cell_structure, propagation, time_step);
  if (interior_pass) {
    cell_structure.for_each_interior_particle(kernel);
  } else {
    cell_structure.for_each_boundary_particle(kernel);
  }
}

int System::System::integrate(int n_steps, int reuse_forces) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  auto &propagation = *this->propagation;
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  auto const has_vs_rel = [&propagation]() {
    return propagation.used_propagations &
           (PropagationMode::ROT_VS_RELATIVE |
            PropagationMode::ROT_VS_INDEPENDENT |
            PropagationMode::TRANS_VS_RELATIVE);
  };
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
  auto const has_vs_com = [&propagation]() {
    return propagation.used_propagations &
           (PropagationMode::TRANS_VS_CENTER_OF_MASS);
  };
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto const n_rigid_bonds = bonded_ias->get_n_rigid_bonds();
#endif

  // Prepare particle structure and run sanity checks of all active algorithms
  propagation.update_default_propagation(thermostat->thermo_switch);
  update_used_propagations();
  on_integration_start();

  // If any method vetoes (e.g. P3M not initialized), immediately bail out
  if (check_runtime_errors(comm_cart))
    return INTEG_ERROR_RUNTIME;

  // Additional preparations for the first integration step
  if (reuse_forces == INTEG_REUSE_FORCES_NEVER or
      ((reuse_forces != INTEG_REUSE_FORCES_ALWAYS) and
       propagation.recalc_forces)) {
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("Initial Force Calculation");
#endif
    thermostat->lb_coupling_deactivate();

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    if (has_vs_rel()) {
      vs_relative_update_particles(*cell_structure, *box_geo);
    }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
    if (has_vs_com()) {
      vs_com_update_particles(*cell_structure, *box_geo);
    }
#endif

    // Communication step: distribute ghost positions
    cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());

    calculate_forces();

    // If calculate_forces started a split-phase ghost reduce, finish it now:
    // the initial-force path has no step_2 to overlap with (n_steps may be 0,
    // or the forces are for the previous step's record), so we just complete
    // the reduce immediately.
    if (cell_structure->has_pending_ghost_reduce()) {
      cell_structure->ghosts_reduce_forces_finish();
    }

    if (propagation.integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
#ifdef ESPRESSO_ROTATION
      convert_initial_torques(cell_structure->local_particles());
#endif
    }

#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("Initial Force Calculation");
#endif
  }

  thermostat->lb_coupling_activate();

  if (check_runtime_errors(comm_cart))
    return INTEG_ERROR_RUNTIME;

  // Keep track of the number of Verlet updates (i.e. particle resorts)
  int n_verlet_updates = 0;

  // Keep track of whether an interrupt signal was caught (only in singleton
  // mode, since signal handlers are unreliable with more than 1 MPI rank)
  auto const singleton_mode = comm_cart.size() == 1;
  auto caught_sigint = false;
  auto caught_error = false;

  auto lb_active = false;
  auto ek_active = false;
  if (propagation.integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
    lb_active = lb.is_solver_set();
    ek_active = ek.is_ready_for_propagation();
  }
  auto const calc_md_steps_per_tau = [this](double tau) {
    return static_cast<int>(std::round(tau / time_step));
  };

#ifdef ESPRESSO_VALGRIND
  CALLGRIND_START_INSTRUMENTATION;
#endif
  // Integration loop
#ifdef ESPRESSO_CALIPER
  EspressoCaliLoop espresso_cali_integration_loop("Integration loop");
#endif
  int integrated_steps = 0;
  for (int step = 0; step < n_steps; step++) {
#ifdef ESPRESSO_CALIPER
    auto espresso_cali_integration_iter =
        espresso_cali_integration_loop.iteration(step);
#endif

    // Ensure every local/ghost particle has a valid ParticleStore row before
    // integrator_step_1 reads previous-step forces. Mid-step particle creation
    // (collision handling, bond breakage) at the end of the previous iteration
    // would otherwise leave new particles rowless. O(1) when the store is
    // clean; rank-local.
    cell_structure->ensure_particle_store_synchronized();

#ifdef ESPRESSO_BOND_CONSTRAINT
    if (n_rigid_bonds)
      save_old_position(cell_structure->local_particles(),
                        cell_structure->ghost_particles());
#endif

    lees_edwards->update_box_params(*box_geo, sim_time);
    bool early_exit =
        integrator_step_1(*cell_structure, propagation, *this, time_step);
    if (early_exit)
      break;

    sim_time += time_step;
    if (box_geo->type() == BoxType::LEES_EDWARDS) {
      auto const kernel = LeesEdwards::Push{*box_geo};
      cell_structure->for_each_local_particle(
          [&kernel](Particle &p) { kernel(p); });
    }

#ifdef ESPRESSO_NPT
    if (not has_npt_enabled())
#endif
    {
      resort_particles_if_needed(*this);
    }
    // Propagate philox RNG counters
    thermostat->philox_counter_increment();

#ifdef ESPRESSO_BOND_CONSTRAINT
    // Correct particle positions that participate in a rigid/constrained bond
    if (n_rigid_bonds) {
      correct_position_shake(*cell_structure, *box_geo, *bonded_ias);
    }
#endif

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    if (has_vs_rel()) {
#ifdef ESPRESSO_NPT
      if (has_npt_enabled()) {
        cell_structure->update_ghosts_and_resort_particle(
            Cells::DATA_PART_PROPERTIES);
      }
#endif // ESPRESSO_NPT
      vs_relative_update_particles(*cell_structure, *box_geo);
    }
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
    if (has_vs_com()) {
#ifdef ESPRESSO_NPT
      if (has_npt_enabled()) {
        cell_structure->update_ghosts_and_resort_particle(
            Cells::DATA_PART_PROPERTIES);
      }
#endif // ESPRESSO_NPT
      vs_com_update_particles(*cell_structure, *box_geo);
    }
#endif // ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS

    if (cell_structure->get_resort_particles() >= Cells::RESORT_LOCAL)
      n_verlet_updates++;

    // Communication step: distribute ghost positions
    cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    integrate_magnetodynamics();
#endif

    calculate_forces();

#ifdef ESPRESSO_VIRTUAL_SITES_INERTIALESS_TRACERS
    if (thermostat->lb and
        (propagation.used_propagations & PropagationMode::TRANS_LB_TRACER)) {
      // LB-tracer arm is ineligible for the split path; this block only runs
      // on the blocking path where has_pending_ghost_reduce() is false.
      assert(not cell_structure->has_pending_ghost_reduce() &&
             "LB-tracer arm must be inactive on the split-phase path");
      lb_tracers_add_particle_force_to_fluid(*cell_structure, *box_geo,
                                             *local_geo, lb);
    }
#endif
    if (cell_structure->has_pending_ghost_reduce()) {
      // Split-phase path: interior half-kick runs while the ghost reduce is
      // in flight; boundary half-kick runs after the reduce finishes.
      // NPT, BD, steepest-descent, and LB-tracer are ineligible and never
      // reach this branch (asserted inside integrator_step_2_filtered).
      //
      // RAII guard: if the interior pass throws (bad_alloc, Kokkos error, user
      // callback), finish the pending reduce on unwind so MPI requests are not
      // left dangling and the next start-assert does not fire.
      // ESPResSo normally propagates errors via runtimeErrorMsg rather than
      // exceptions, so this guard fires only in exceptional circumstances; the
      // normal path calls finish() explicitly and the guard becomes a no-op.
      struct ReduceGuard {
        CellStructure *cs;
        bool active;
        ~ReduceGuard() {
          if (active and cs->has_pending_ghost_reduce()) {
            try {
              cs->ghosts_reduce_forces_finish();
            } catch (...) { // NOLINT(bugprone-empty-catch)
              // The guard only runs during unwind from another exception;
              // letting a second one escape this (implicitly noexcept)
              // destructor would call std::terminate. Keep the original.
            }
          }
        }
      } guard{cell_structure.get(), true};

      integrator_step_2_filtered(*cell_structure, propagation, time_step,
                                 /*interior_pass=*/true);
      cell_structure->ghosts_reduce_forces_finish();
      guard.active = false; // normal path: finish already called above
      integrator_step_2_filtered(*cell_structure, propagation, time_step,
                                 /*interior_pass=*/false);
    } else {
      integrator_step_2(*cell_structure, propagation, *this, time_step);
    }
    if (propagation.integ_switch == INTEG_METHOD_BD) {
      resort_particles_if_needed(*this);
    }
    if (box_geo->type() == BoxType::LEES_EDWARDS) {
      auto const kernel = LeesEdwards::UpdateOffset{*box_geo};
      cell_structure->for_each_local_particle(
          [&kernel](Particle &p) { kernel(p); });
    }
#ifdef ESPRESSO_BOND_CONSTRAINT
    if (n_rigid_bonds) {
      correct_velocity_shake(*cell_structure, *box_geo, *bonded_ias);
    }
#endif

    // propagate one-step functionalities
    if (propagation.integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      if (lb_active and ek_active) {
        // assume that they are coupled, which is not necessarily true
        auto const md_steps_per_lb_step = calc_md_steps_per_tau(lb.get_tau());
        auto const md_steps_per_ek_step = calc_md_steps_per_tau(ek.get_tau());

        if (md_steps_per_lb_step != md_steps_per_ek_step) {
          runtimeErrorMsg()
              << "LB and EK are active but with different time steps.";
        }

        assert(lb.is_gpu() == ek.is_gpu());
        assert(propagation.lb_skipped_md_steps ==
               propagation.ek_skipped_md_steps);

        propagation.lb_skipped_md_steps += 1;
        propagation.ek_skipped_md_steps += 1;
        if (propagation.lb_skipped_md_steps >= md_steps_per_lb_step) {
          propagation.lb_skipped_md_steps = 0;
          propagation.ek_skipped_md_steps = 0;
#ifdef ESPRESSO_CALIPER
          ESPRESSO_CALI_MARK_BEGIN("lb_propagation");
#endif
          lb.propagate();
          lb.ghost_communication_vel();
#ifdef ESPRESSO_CALIPER
          ESPRESSO_CALI_MARK_END("lb_propagation");
#endif
#ifdef ESPRESSO_CALIPER
          ESPRESSO_CALI_MARK_BEGIN("ek_propagation");
#endif
          ek.propagate();
#ifdef ESPRESSO_CALIPER
          ESPRESSO_CALI_MARK_END("ek_propagation");
#endif
        }
      } else if (lb_active) {
        auto const md_steps_per_lb_step = calc_md_steps_per_tau(lb.get_tau());
        propagation.lb_skipped_md_steps += 1;
        if (propagation.lb_skipped_md_steps >= md_steps_per_lb_step) {
          propagation.lb_skipped_md_steps = 0;
#ifdef ESPRESSO_CALIPER
          ESPRESSO_CALI_MARK_BEGIN("lb_propagation");
#endif
          lb.propagate();
#ifdef ESPRESSO_CALIPER
          ESPRESSO_CALI_MARK_END("lb_propagation");
#endif
        }
      } else if (ek_active) {
        auto const md_steps_per_ek_step = calc_md_steps_per_tau(ek.get_tau());
        propagation.ek_skipped_md_steps += 1;
        if (propagation.ek_skipped_md_steps >= md_steps_per_ek_step) {
          propagation.ek_skipped_md_steps = 0;
#ifdef ESPRESSO_CALIPER
          ESPRESSO_CALI_MARK_BEGIN("ek_propagation");
#endif
          ek.propagate();
#ifdef ESPRESSO_CALIPER
          ESPRESSO_CALI_MARK_END("ek_propagation");
#endif
        }
      }
      if (lb_active and (propagation.used_propagations &
                         PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE)) {
        thermostat->lb->rng_increment();
      }

#ifdef ESPRESSO_VIRTUAL_SITES_INERTIALESS_TRACERS
      if (thermostat->lb and
          (propagation.used_propagations & PropagationMode::TRANS_LB_TRACER)) {
#ifdef ESPRESSO_CALIPER
        ESPRESSO_CALI_MARK_BEGIN("lb_tracers_propagation");
#endif
        if (lb_active) {
          lb.ghost_communication_vel();
        }
        lb_tracers_propagate(*cell_structure, lb, time_step);
#ifdef ESPRESSO_CALIPER
        ESPRESSO_CALI_MARK_END("lb_tracers_propagation");
#endif
      }
#endif

#ifdef ESPRESSO_COLLISION_DETECTION
      cell_structure->clear_new_bonds();
      collision_detection->handle_collisions();
      cell_structure->rebuild_bond_list();
#endif
      bond_breakage->process_queue(*this);
    }

    integrated_steps++;

    if (check_runtime_errors(comm_cart)) {
      caught_error = true;
      break;
    }

    // Check if SIGINT has been caught.
    if (singleton_mode and ctrl_C == 1) {
      caught_sigint = true;
      break;
    }

  } // for-loop over integration steps
  if (lb_active) {
    lb.ghost_communication();
  }
  lees_edwards->update_box_params(*box_geo, sim_time);
  // espresso_cali_integration_loop destructor ends the Caliper loop region.

#ifdef ESPRESSO_VALGRIND
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  if (has_vs_rel()) {
    vs_relative_update_particles(*cell_structure, *box_geo);
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
  if (has_vs_com()) {
    vs_com_update_particles(*cell_structure, *box_geo);
  }
#endif

  // Verlet list statistics
  cell_structure->update_verlet_stats(n_steps, n_verlet_updates);

#ifdef ESPRESSO_NPT
  if (has_npt_enabled()) {
    synchronize_npt_state();
  }
#endif
  if (caught_sigint) {
    ctrl_C = 0;
    return INTEG_ERROR_SIGINT;
  }
  if (caught_error) {
    return INTEG_ERROR_RUNTIME;
  }
  if (boost::mpi::all_reduce(::comm_cart, not cell_structure->use_verlet_list,
                             std::logical_or<>())) {
    cell_structure->use_verlet_list = false;
  }
  return integrated_steps;
}

int System::System::integrate_with_signal_handler(int n_steps, int reuse_forces,
                                                  bool update_accumulators) {
  assert(n_steps >= 0);

  // Override the signal handler so that the integrator obeys Ctrl+C
  SignalHandler sa(SIGINT, [](int) { ctrl_C = 1; });

  /* if skin wasn't set, do an educated guess now */
  if (not cell_structure->is_verlet_skin_set()) {
    try {
      cell_structure->set_verlet_skin_heuristic();
    } catch (...) {
      if (comm_cart.rank() == 0) {
        throw;
      }
      return INTEG_ERROR_RUNTIME;
    }
  }

  if (not update_accumulators or n_steps == 0) {
    return integrate(n_steps, reuse_forces);
  }

  for (int i = 0; i < n_steps;) {
    /* Integrate to either the next accumulator update, or the
     * end, depending on what comes first. */
    auto const steps =
        std::min((n_steps - i), auto_update_accumulators->next_update());

    auto const local_retval = integrate(steps, reuse_forces);

    // make sure all ranks exit when one rank fails
    std::remove_const_t<decltype(local_retval)> global_retval;
    boost::mpi::all_reduce(comm_cart, local_retval, global_retval,
                           std::plus<int>());
    if (global_retval < 0) {
      return global_retval; // propagate error code
    }

    reuse_forces = INTEG_REUSE_FORCES_ALWAYS;

    (*auto_update_accumulators)(comm_cart, steps);

    i += steps;
  }

  return 0;
}

void System::System::set_sim_time(double value) {
  sim_time = value;
  propagation->recalc_forces = true;
  lees_edwards->update_box_params(*box_geo, sim_time);
}
