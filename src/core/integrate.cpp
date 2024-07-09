/*
 * Copyright (C) 2010-2023 The ESPResSo project
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
#include "integrators/velocity_verlet_inline.hpp"
#include "integrators/velocity_verlet_npt.hpp"

#include "BoxGeometry.hpp"
#include "ParticleRange.hpp"
#include "PropagationMode.hpp"
#include "accumulators/AutoUpdateAccumulators.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "lb/particle_coupling.hpp"
#include "lb/utils.hpp"
#include "lees_edwards/lees_edwards.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "rattle.hpp"
#include "rotation.hpp"
#include "signalhandling.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "thermostats/langevin_inline.hpp"
#include "virtual_sites/lb_tracers.hpp"
#include "virtual_sites/relative.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>

#ifdef CALIPER
#include <caliper/cali.h>
#endif

#ifdef VALGRIND
#include <callgrind.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <csignal>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef WALBERLA
#ifdef WALBERLA_STATIC_ASSERT
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
  case INTEG_METHOD_NVT: {
    // NOLINTNEXTLINE(bugprone-branch-clone)
    if ((thermo_switch & THERMO_LB) and (thermo_switch & THERMO_LANGEVIN)) {
      default_propagation = PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE;
#ifdef ROTATION
      default_propagation |= PropagationMode::ROT_LANGEVIN;
#endif
    } else if (thermo_switch & THERMO_LB) {
      default_propagation = PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE;
#ifdef ROTATION
      default_propagation |= PropagationMode::ROT_EULER;
#endif
    } else if (thermo_switch & THERMO_LANGEVIN) {
      default_propagation = PropagationMode::TRANS_LANGEVIN;
#ifdef ROTATION
      default_propagation |= PropagationMode::ROT_LANGEVIN;
#endif
    } else {
      default_propagation = PropagationMode::TRANS_NEWTON;
#ifdef ROTATION
      default_propagation |= PropagationMode::ROT_EULER;
#endif
    }
    break;
  }
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    default_propagation = PropagationMode::TRANS_LANGEVIN_NPT;
    break;
#endif
  case INTEG_METHOD_BD:
    default_propagation = PropagationMode::TRANS_BROWNIAN;
#ifdef ROTATION
    default_propagation |= PropagationMode::ROT_BROWNIAN;
#endif
    break;
#ifdef STOKESIAN_DYNAMICS
  case INTEG_METHOD_SD:
    default_propagation = PropagationMode::TRANS_STOKESIAN;
    break;
#endif // STOKESIAN_DYNAMICS
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
#ifdef NPT
  if (propagation->used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) {
    if (thermo_switch != THERMO_NPT_ISO) {
      runtimeErrorMsg() << "The NpT integrator requires the NpT thermostat";
    }
    if (box_geo->type() == BoxType::LEES_EDWARDS) {
      runtimeErrorMsg() << "The NpT integrator cannot use Lees-Edwards";
    }
  }
#endif
  if (propagation->used_propagations & PropagationMode::TRANS_BROWNIAN) {
    if (thermo_switch != THERMO_BROWNIAN) {
      runtimeErrorMsg() << "The BD integrator requires the BD thermostat";
    }
  }
  if (propagation->used_propagations & PropagationMode::TRANS_STOKESIAN) {
#ifdef STOKESIAN_DYNAMICS
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
  if (::bonded_ia_params.get_n_thermalized_bonds() >= 1 and
      (thermostat->thermalized_bond == nullptr or
       (thermo_switch & THERMO_BOND) == 0)) {
    runtimeErrorMsg()
        << "Thermalized bonds require the thermalized_bond thermostat";
  }

#ifdef ROTATION
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
#endif
}

#ifdef WALBERLA
void walberla_tau_sanity_checks(std::string method, double tau,
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

void walberla_agrid_sanity_checks(std::string method,
                                  Utils::Vector3d const &geo_left,
                                  Utils::Vector3d const &geo_right,
                                  Utils::Vector3d const &lattice_left,
                                  Utils::Vector3d const &lattice_right,
                                  double agrid) {
  // waLBerla and ESPResSo must agree on domain decomposition
  auto const tol = agrid / 1E6;
  if ((lattice_left - geo_left).norm2() > tol or
      (lattice_right - geo_right).norm2() > tol) {
    runtimeErrorMsg() << "\nMPI rank " << ::this_node << ": "
                      << "left ESPResSo: [" << geo_left << "], "
                      << "left waLBerla: [" << lattice_left << "]"
                      << "\nMPI rank " << ::this_node << ": "
                      << "right ESPResSo: [" << geo_right << "], "
                      << "right waLBerla: [" << lattice_right << "]"
                      << "\nfor method: " << method;
    throw std::runtime_error(
        "waLBerla and ESPResSo disagree about domain decomposition.");
  }
}
#endif // WALBERLA

static void resort_particles_if_needed(System::System &system) {
  auto &cell_structure = *system.cell_structure;
  auto const offset = LeesEdwards::verlet_list_offset(
      *system.box_geo, cell_structure.get_le_pos_offset_at_last_resort());
  if (cell_structure.check_resort_required(offset)) {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

void System::System::thermostat_force_init() {
  auto const &propagation = *this->propagation;
  if ((not thermostat->langevin) or ((propagation.used_propagations &
                                      (PropagationMode::TRANS_LANGEVIN |
                                       PropagationMode::ROT_LANGEVIN)) == 0)) {
    return;
  }
  auto const &langevin = *thermostat->langevin;
  auto const kT = thermostat->kT;
  for (auto &p : cell_structure->local_particles()) {
    if (propagation.should_propagate_with(p, PropagationMode::TRANS_LANGEVIN))
      p.force() += friction_thermo_langevin(langevin, p, time_step, kT);
#ifdef ROTATION
    if (propagation.should_propagate_with(p, PropagationMode::ROT_LANGEVIN))
      p.torque() += convert_vector_body_to_space(
          p, friction_thermo_langevin_rotation(langevin, p, time_step, kT));
#endif
  }
}

/** @brief Calls the hook for propagation kernels before the force calculation
 *  @return whether or not to stop the integration loop early.
 */
static bool integrator_step_1(ParticleRange const &particles,
                              Propagation const &propagation,
                              System::System &system, double time_step) {
  // steepest decent
  if (propagation.integ_switch == INTEG_METHOD_STEEPEST_DESCENT)
    return steepest_descent_step(particles);

  auto const &thermostat = *system.thermostat;
  auto const kT = thermostat.kT;
  for (auto &p : particles) {
#ifdef VIRTUAL_SITES
    // virtual sites are updated later in the integration loop
    if (p.is_virtual())
      continue;
#endif
    if (propagation.should_propagate_with(
            p, PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE))
      velocity_verlet_propagator_1(p, time_step);
    if (propagation.should_propagate_with(p, PropagationMode::TRANS_NEWTON))
      velocity_verlet_propagator_1(p, time_step);
#ifdef ROTATION
    if (propagation.should_propagate_with(p, PropagationMode::ROT_EULER))
      velocity_verlet_rotator_1(p, time_step);
#endif
    if (propagation.should_propagate_with(p, PropagationMode::TRANS_LANGEVIN))
      velocity_verlet_propagator_1(p, time_step);
#ifdef ROTATION
    if (propagation.should_propagate_with(p, PropagationMode::ROT_LANGEVIN))
      velocity_verlet_rotator_1(p, time_step);
#endif
    if (propagation.should_propagate_with(p, PropagationMode::TRANS_BROWNIAN))
      brownian_dynamics_propagator(*thermostat.brownian, p, time_step, kT);
#ifdef ROTATION
    if (propagation.should_propagate_with(p, PropagationMode::ROT_BROWNIAN))
      brownian_dynamics_rotator(*thermostat.brownian, p, time_step, kT);
#endif
  }

#ifdef NPT
  if ((propagation.used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) and
      (propagation.default_propagation & PropagationMode::TRANS_LANGEVIN_NPT)) {
    auto pred = PropagationPredicateNPT(propagation.default_propagation);
    velocity_verlet_npt_step_1(particles.filter(pred), *thermostat.npt_iso,
                               time_step, system);
  }
#endif

#ifdef STOKESIAN_DYNAMICS
  if ((propagation.used_propagations & PropagationMode::TRANS_STOKESIAN) and
      (propagation.default_propagation & PropagationMode::TRANS_STOKESIAN)) {
    auto pred = PropagationPredicateStokesian(propagation.default_propagation);
    stokesian_dynamics_step_1(particles.filter(pred), *thermostat.stokesian,
                              time_step, kT);
  }
#endif // STOKESIAN_DYNAMICS

  return false;
}

static void integrator_step_2(ParticleRange const &particles,
                              Propagation const &propagation,
                              Thermostat::Thermostat const &thermostat,
                              double time_step) {
  if (propagation.integ_switch == INTEG_METHOD_STEEPEST_DESCENT)
    return;

  for (auto &p : particles) {
#ifdef VIRTUAL_SITES
    // virtual sites are updated later in the integration loop
    if (p.is_virtual())
      continue;
#endif
    if (propagation.should_propagate_with(
            p, PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE))
      velocity_verlet_propagator_2(p, time_step);
    if (propagation.should_propagate_with(p, PropagationMode::TRANS_NEWTON))
      velocity_verlet_propagator_2(p, time_step);
#ifdef ROTATION
    if (propagation.should_propagate_with(p, PropagationMode::ROT_EULER))
      velocity_verlet_rotator_2(p, time_step);
#endif
    if (propagation.should_propagate_with(p, PropagationMode::TRANS_LANGEVIN))
      velocity_verlet_propagator_2(p, time_step);
#ifdef ROTATION
    if (propagation.should_propagate_with(p, PropagationMode::ROT_LANGEVIN))
      velocity_verlet_rotator_2(p, time_step);
#endif
  }

#ifdef NPT
  if ((propagation.used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) and
      (propagation.default_propagation & PropagationMode::TRANS_LANGEVIN_NPT)) {
    auto pred = PropagationPredicateNPT(propagation.default_propagation);
    velocity_verlet_npt_step_2(particles.filter(pred), *thermostat.npt_iso,
                               time_step);
  }
#endif
}

int System::System::integrate(int n_steps, int reuse_forces) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  auto &propagation = *this->propagation;
#ifdef VIRTUAL_SITES_RELATIVE
  auto const has_vs_rel = [&propagation]() {
    return propagation.used_propagations & (PropagationMode::ROT_VS_RELATIVE |
                                            PropagationMode::TRANS_VS_RELATIVE);
  };
#endif
#ifdef BOND_CONSTRAINT
  auto const n_rigid_bonds = ::bonded_ia_params.get_n_rigid_bonds();
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
#ifdef CALIPER
    CALI_MARK_BEGIN("Initial Force Calculation");
#endif
    thermostat->lb_coupling_deactivate();

#ifdef VIRTUAL_SITES_RELATIVE
    if (has_vs_rel()) {
      vs_relative_update_particles(*cell_structure, *box_geo);
    }
#endif

    // Communication step: distribute ghost positions
    cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());

    calculate_forces();

    if (propagation.integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
#ifdef ROTATION
      convert_initial_torques(cell_structure->local_particles());
#endif
    }

#ifdef CALIPER
    CALI_MARK_END("Initial Force Calculation");
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

#ifdef VALGRIND
  CALLGRIND_START_INSTRUMENTATION;
#endif
  // Integration loop
#ifdef CALIPER
  CALI_CXX_MARK_LOOP_BEGIN(integration_loop, "Integration loop");
#endif
  int integrated_steps = 0;
  for (int step = 0; step < n_steps; step++) {
#ifdef CALIPER
    CALI_CXX_MARK_LOOP_ITERATION(integration_loop, step);
#endif

    auto particles = cell_structure->local_particles();

#ifdef BOND_CONSTRAINT
    if (n_rigid_bonds)
      save_old_position(particles, cell_structure->ghost_particles());
#endif

    lees_edwards->update_box_params(*box_geo, sim_time);
    bool early_exit =
        integrator_step_1(particles, propagation, *this, time_step);
    if (early_exit)
      break;

    sim_time += time_step;
    if (box_geo->type() == BoxType::LEES_EDWARDS) {
      auto const kernel = LeesEdwards::Push{*box_geo};
      for (auto &p : particles) {
        kernel(p);
      }
    }

#ifdef NPT
    if (propagation.integ_switch != INTEG_METHOD_NPT_ISO)
#endif
    {
      resort_particles_if_needed(*this);
    }

    // Propagate philox RNG counters
    thermostat->philox_counter_increment();

#ifdef BOND_CONSTRAINT
    // Correct particle positions that participate in a rigid/constrained bond
    if (n_rigid_bonds) {
      correct_position_shake(*cell_structure, *box_geo);
    }
#endif

#ifdef VIRTUAL_SITES_RELATIVE
    if (has_vs_rel()) {
#ifdef NPT
      if (propagation.integ_switch == INTEG_METHOD_NPT_ISO) {
        cell_structure->update_ghosts_and_resort_particle(
            Cells::DATA_PART_PROPERTIES);
      }
#endif // NPT
      vs_relative_update_particles(*cell_structure, *box_geo);
    }
#endif // VIRTUAL_SITES_RELATIVE

    if (cell_structure->get_resort_particles() >= Cells::RESORT_LOCAL)
      n_verlet_updates++;

    // Communication step: distribute ghost positions
    cell_structure->update_ghosts_and_resort_particle(get_global_ghost_flags());

    particles = cell_structure->local_particles();

    calculate_forces();

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
    if (thermostat->lb and
        (propagation.used_propagations & PropagationMode::TRANS_LB_TRACER)) {
      lb_tracers_add_particle_force_to_fluid(*cell_structure, *box_geo,
                                             *local_geo, lb);
    }
#endif
    integrator_step_2(particles, propagation, *thermostat, time_step);
    if (propagation.integ_switch == INTEG_METHOD_BD) {
      resort_particles_if_needed(*this);
    }
    if (box_geo->type() == BoxType::LEES_EDWARDS) {
      auto const kernel = LeesEdwards::UpdateOffset{*box_geo};
      for (auto &p : particles) {
        kernel(p);
      }
    }
#ifdef BOND_CONSTRAINT
    // SHAKE velocity updates
    if (n_rigid_bonds) {
      correct_velocity_shake(*cell_structure, *box_geo);
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

        assert(propagation.lb_skipped_md_steps ==
               propagation.ek_skipped_md_steps);

        propagation.lb_skipped_md_steps += 1;
        propagation.ek_skipped_md_steps += 1;
        if (propagation.lb_skipped_md_steps >= md_steps_per_lb_step) {
          propagation.lb_skipped_md_steps = 0;
          propagation.ek_skipped_md_steps = 0;
          lb.propagate();
          ek.propagate();
        }
      } else if (lb_active) {
        auto const md_steps_per_lb_step = calc_md_steps_per_tau(lb.get_tau());
        propagation.lb_skipped_md_steps += 1;
        if (propagation.lb_skipped_md_steps >= md_steps_per_lb_step) {
          propagation.lb_skipped_md_steps = 0;
          lb.propagate();
        }
      } else if (ek_active) {
        auto const md_steps_per_ek_step = calc_md_steps_per_tau(ek.get_tau());
        propagation.ek_skipped_md_steps += 1;
        if (propagation.ek_skipped_md_steps >= md_steps_per_ek_step) {
          propagation.ek_skipped_md_steps = 0;
          ek.propagate();
        }
      }
      if (lb_active and (propagation.used_propagations &
                         PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE)) {
        thermostat->lb->rng_increment();
      }

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
      if (thermostat->lb and
          (propagation.used_propagations & PropagationMode::TRANS_LB_TRACER)) {
        lb_tracers_propagate(*cell_structure, lb, time_step);
      }
#endif

#ifdef COLLISION_DETECTION
      handle_collisions(*cell_structure);
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
  lees_edwards->update_box_params(*box_geo, sim_time);
#ifdef CALIPER
  CALI_CXX_MARK_LOOP_END(integration_loop);
#endif

#ifdef VALGRIND
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

#ifdef VIRTUAL_SITES_RELATIVE
  if (has_vs_rel()) {
    vs_relative_update_particles(*cell_structure, *box_geo);
  }
#endif

  // Verlet list statistics
  cell_structure->update_verlet_stats(n_steps, n_verlet_updates);

#ifdef NPT
  if (propagation.integ_switch == INTEG_METHOD_NPT_ISO) {
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
