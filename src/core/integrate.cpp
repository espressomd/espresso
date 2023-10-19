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
 *  Molecular dynamics integrator.
 *
 *  For more information about the integrator
 *  see \ref integrate.hpp "integrate.hpp".
 */

#include "integrate.hpp"
#include "integrators/brownian_inline.hpp"
#include "integrators/steepest_descent.hpp"
#include "integrators/stokesian_dynamics_inline.hpp"
#include "integrators/velocity_verlet_inline.hpp"
#include "integrators/velocity_verlet_npt.hpp"

#include "BoxGeometry.hpp"
#include "ParticleRange.hpp"
#include "accumulators.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "bonded_interactions/rigid_bond.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "global_ghost_flags.hpp"
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
#include "virtual_sites.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/range/algorithm/min_element.hpp>

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

int integ_switch = INTEG_METHOD_NVT;

/** Actual simulation time. */
static double sim_time = 0.0;

bool recalc_forces = true;

static int lb_skipped_md_steps = 0;
static int ek_skipped_md_steps = 0;

namespace {
volatile std::sig_atomic_t ctrl_C = 0;
} // namespace

namespace LeesEdwards {
/** @brief Currently active Lees-Edwards protocol. */
static std::shared_ptr<ActiveProtocol> protocol = nullptr;

std::weak_ptr<ActiveProtocol> get_protocol() { return protocol; }

/**
 * @brief Update the Lees-Edwards parameters of the box geometry
 * for the current simulation time.
 */
static void update_box_params(BoxGeometry &box_geo) {
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    assert(protocol != nullptr);
    box_geo.lees_edwards_update(get_pos_offset(sim_time, *protocol),
                                get_shear_velocity(sim_time, *protocol));
  }
}

void set_protocol(std::shared_ptr<ActiveProtocol> new_protocol) {
  auto &system = System::get_system();
  auto &cell_structure = *system.cell_structure;
  auto &box_geo = *system.box_geo;
  box_geo.set_type(BoxType::LEES_EDWARDS);
  protocol = std::move(new_protocol);
  LeesEdwards::update_box_params(box_geo);
  ::recalc_forces = true;
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
}

void unset_protocol() {
  auto &system = System::get_system();
  auto &cell_structure = *system.cell_structure;
  auto &box_geo = *system.box_geo;
  protocol = nullptr;
  box_geo.set_type(BoxType::CUBOID);
  ::recalc_forces = true;
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
}

template <class Kernel> void run_kernel(BoxGeometry const &box_geo) {
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    auto &system = System::get_system();
    auto &cell_structure = *system.cell_structure;
    auto const kernel = Kernel{box_geo};
    auto const particles = cell_structure.local_particles();
    std::for_each(particles.begin(), particles.end(),
                  [&kernel](auto &p) { kernel(p); });
  }
}
} // namespace LeesEdwards

void integrator_sanity_checks() {
  if (get_time_step() < 0.0) {
    runtimeErrorMsg() << "time_step not set";
  }
  switch (integ_switch) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    if (thermo_switch != THERMO_OFF)
      runtimeErrorMsg()
          << "The steepest descent integrator is incompatible with thermostats";
    break;
  case INTEG_METHOD_NVT:
    if (thermo_switch & (THERMO_NPT_ISO | THERMO_BROWNIAN | THERMO_SD))
      runtimeErrorMsg() << "The VV integrator is incompatible with the "
                           "currently active combination of thermostats";
    break;
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    if (thermo_switch != THERMO_OFF and thermo_switch != THERMO_NPT_ISO)
      runtimeErrorMsg() << "The NpT integrator requires the NpT thermostat";
    if (System::get_system().box_geo->type() == BoxType::LEES_EDWARDS)
      runtimeErrorMsg() << "The NpT integrator cannot use Lees-Edwards";
    break;
#endif
  case INTEG_METHOD_BD:
    if (thermo_switch != THERMO_BROWNIAN)
      runtimeErrorMsg() << "The BD integrator requires the BD thermostat";
    break;
#ifdef STOKESIAN_DYNAMICS
  case INTEG_METHOD_SD:
    if (thermo_switch != THERMO_OFF and thermo_switch != THERMO_SD)
      runtimeErrorMsg() << "The SD integrator requires the SD thermostat";
    break;
#endif
  default:
    runtimeErrorMsg() << "Unknown value for integ_switch";
  }
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

void walberla_tau_sanity_checks(std::string method, double tau) {
  walberla_tau_sanity_checks(std::move(method), tau, get_time_step());
}

void walberla_agrid_sanity_checks(std::string method,
                                  Utils::Vector3d const &lattice_left,
                                  Utils::Vector3d const &lattice_right,
                                  double agrid) {
  // waLBerla and ESPResSo must agree on domain decomposition
  auto const &system = System::get_system();
  auto const &my_left = system.local_geo->my_left();
  auto const &my_right = system.local_geo->my_right();
  auto const tol = agrid / 1E6;
  if ((lattice_left - my_left).norm2() > tol or
      (lattice_right - my_right).norm2() > tol) {
    runtimeErrorMsg() << "\nMPI rank " << ::this_node << ": "
                      << "left ESPResSo: [" << my_left << "], "
                      << "left waLBerla: [" << lattice_left << "]"
                      << "\nMPI rank " << ::this_node << ": "
                      << "right ESPResSo: [" << my_right << "], "
                      << "right waLBerla: [" << lattice_right << "]";
    throw std::runtime_error(
        "waLBerla and ESPResSo disagree about domain decomposition.");
  }
}
#endif

static auto calc_md_steps_per_tau(double tau) {
  return static_cast<int>(std::round(tau / get_time_step()));
}

static void resort_particles_if_needed(ParticleRange const &particles) {
  auto &system = System::get_system();
  auto &cell_structure = *system.cell_structure;
  auto const verlet_skin = system.get_verlet_skin();
  auto const offset = LeesEdwards::verlet_list_offset(
      *system.box_geo, cell_structure.get_le_pos_offset_at_last_resort());
  if (cell_structure.check_resort_required(particles, verlet_skin, offset)) {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

/** @brief Calls the hook for propagation kernels before the force calculation
 *  @return whether or not to stop the integration loop early.
 */
static bool integrator_step_1(ParticleRange const &particles,
                              double time_step) {
  bool early_exit = false;
  switch (integ_switch) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    early_exit = steepest_descent_step(particles);
    break;
  case INTEG_METHOD_NVT:
    velocity_verlet_step_1(particles, time_step);
    break;
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    velocity_verlet_npt_step_1(particles, time_step);
    break;
#endif
  case INTEG_METHOD_BD:
    // the Ermak-McCammon's Brownian Dynamics requires a single step
    // so, just skip here
    break;
#ifdef STOKESIAN_DYNAMICS
  case INTEG_METHOD_SD:
    stokesian_dynamics_step_1(particles, time_step);
    break;
#endif // STOKESIAN_DYNAMICS
  default:
    throw std::runtime_error("Unknown value for integ_switch");
  }
  return early_exit;
}

/** Calls the hook of the propagation kernels after force calculation */
static void integrator_step_2(ParticleRange const &particles, double kT,
                              double time_step) {
  switch (integ_switch) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    // Nothing
    break;
  case INTEG_METHOD_NVT:
    velocity_verlet_step_2(particles, time_step);
    break;
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    velocity_verlet_npt_step_2(particles, time_step);
    break;
#endif
  case INTEG_METHOD_BD:
    // the Ermak-McCammon's Brownian Dynamics requires a single step
    brownian_dynamics_propagator(brownian, particles, time_step, kT);
    resort_particles_if_needed(particles);
    break;
#ifdef STOKESIAN_DYNAMICS
  case INTEG_METHOD_SD:
    // Nothing
    break;
#endif // STOKESIAN_DYNAMICS
  default:
    throw std::runtime_error("Unknown value for INTEG_SWITCH");
  }
}

int integrate(System::System &system, int n_steps, int reuse_forces) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  auto &cell_structure = *system.cell_structure;
  auto &box_geo = *system.box_geo;
  auto const time_step = system.get_time_step();

  // Prepare particle structure and run sanity checks of all active algorithms
  system.on_integration_start();

  // If any method vetoes (e.g. P3M not initialized), immediately bail out
  if (check_runtime_errors(comm_cart))
    return INTEG_ERROR_RUNTIME;

  // Additional preparations for the first integration step
  if (reuse_forces == INTEG_REUSE_FORCES_NEVER or
      (recalc_forces and reuse_forces != INTEG_REUSE_FORCES_ALWAYS)) {
#ifdef CALIPER
    CALI_MARK_BEGIN("Initial Force Calculation");
#endif
    lb_lbcoupling_deactivate();

#ifdef VIRTUAL_SITES
    virtual_sites()->update();
#endif

    // Communication step: distribute ghost positions
    cells_update_ghosts(global_ghost_flags());

    force_calc(system, time_step, temperature);

    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
#ifdef ROTATION
      convert_initial_torques(cell_structure.local_particles());
#endif
    }

#ifdef CALIPER
    CALI_MARK_END("Initial Force Calculation");
#endif
  }

  lb_lbcoupling_activate();

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
  if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
    lb_active = system.lb.is_solver_set();
    ek_active = system.ek.is_ready_for_propagation();
  }

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

    auto particles = cell_structure.local_particles();

#ifdef BOND_CONSTRAINT
    if (n_rigidbonds)
      save_old_position(particles, cell_structure.ghost_particles());
#endif

    LeesEdwards::update_box_params(box_geo);
    bool early_exit = integrator_step_1(particles, time_step);
    if (early_exit)
      break;

    LeesEdwards::run_kernel<LeesEdwards::Push>(box_geo);

#ifdef NPT
    if (integ_switch != INTEG_METHOD_NPT_ISO)
#endif
    {
      resort_particles_if_needed(particles);
    }

    // Propagate philox RNG counters
    philox_counter_increment();

#ifdef BOND_CONSTRAINT
    // Correct particle positions that participate in a rigid/constrained bond
    if (n_rigidbonds) {
      correct_position_shake(cell_structure, box_geo);
    }
#endif

#ifdef VIRTUAL_SITES
    virtual_sites()->update();
#endif

    if (cell_structure.get_resort_particles() >= Cells::RESORT_LOCAL)
      n_verlet_updates++;

    // Communication step: distribute ghost positions
    cells_update_ghosts(global_ghost_flags());

    particles = cell_structure.local_particles();

    force_calc(system, time_step, temperature);

#ifdef VIRTUAL_SITES
    virtual_sites()->after_force_calc(time_step);
#endif
    integrator_step_2(particles, temperature, time_step);
    LeesEdwards::run_kernel<LeesEdwards::UpdateOffset>(box_geo);
#ifdef BOND_CONSTRAINT
    // SHAKE velocity updates
    if (n_rigidbonds) {
      correct_velocity_shake(cell_structure, box_geo);
    }
#endif

    // propagate one-step functionalities
    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      if (lb_active and ek_active) {
        auto &lb = system.lb;
        auto &ek = system.ek;
        // assume that they are coupled, which is not necessarily true
        auto const md_steps_per_lb_step = calc_md_steps_per_tau(lb.get_tau());
        auto const md_steps_per_ek_step = calc_md_steps_per_tau(ek.get_tau());

        if (md_steps_per_lb_step != md_steps_per_ek_step) {
          runtimeErrorMsg()
              << "LB and EK are active but with different time steps.";
        }

        assert(lb_skipped_md_steps == ek_skipped_md_steps);

        lb_skipped_md_steps += 1;
        ek_skipped_md_steps += 1;
        if (lb_skipped_md_steps >= md_steps_per_lb_step) {
          lb_skipped_md_steps = 0;
          ek_skipped_md_steps = 0;
          lb.propagate();
          ek.propagate();
        }
        lb_lbcoupling_propagate();
      } else if (lb_active) {
        auto &lb = system.lb;
        auto const md_steps_per_lb_step = calc_md_steps_per_tau(lb.get_tau());
        lb_skipped_md_steps += 1;
        if (lb_skipped_md_steps >= md_steps_per_lb_step) {
          lb_skipped_md_steps = 0;
          lb.propagate();
        }
        lb_lbcoupling_propagate();
      } else if (ek_active) {
        auto &ek = system.ek;
        auto const md_steps_per_ek_step = calc_md_steps_per_tau(ek.get_tau());
        ek_skipped_md_steps += 1;
        if (ek_skipped_md_steps >= md_steps_per_ek_step) {
          ek_skipped_md_steps = 0;
          ek.propagate();
        }
      }

#ifdef VIRTUAL_SITES
      virtual_sites()->after_lb_propagation(time_step);
#endif

#ifdef COLLISION_DETECTION
      handle_collisions(cell_structure);
#endif
      BondBreakage::process_queue();
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
  LeesEdwards::update_box_params(box_geo);
#ifdef CALIPER
  CALI_CXX_MARK_LOOP_END(integration_loop);
#endif

#ifdef VALGRIND
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

#ifdef VIRTUAL_SITES
  virtual_sites()->update();
#endif

  // Verlet list statistics
  system.update_verlet_stats(n_steps, n_verlet_updates);

#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
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

int integrate_with_signal_handler(System::System &system, int n_steps,
                                  int reuse_forces, bool update_accumulators) {
  assert(n_steps >= 0);

  // Override the signal handler so that the integrator obeys Ctrl+C
  SignalHandler sa(SIGINT, [](int) { ctrl_C = 1; });

  /* if skin wasn't set, do an educated guess now */
  if (not system.is_verlet_skin_set()) {
    auto const max_cut = system.maximal_cutoff();
    if (max_cut <= 0.0) {
      if (comm_cart.rank() == 0) {
        throw std::runtime_error(
            "cannot automatically determine skin, please set it manually");
      }
      return INTEG_ERROR_RUNTIME;
    }
    auto &cell_structure = *system.cell_structure;
    /* maximal skin that can be used without resorting is the maximal
     * range of the cell system minus what is needed for interactions. */
    auto const max_range = *boost::min_element(cell_structure.max_cutoff());
    auto const new_skin = std::min(0.4 * max_cut, max_range - max_cut);
    system.set_verlet_skin(new_skin);
  }

  if (not update_accumulators or n_steps == 0) {
    return integrate(system, n_steps, reuse_forces);
  }

  using Accumulators::auto_update;
  using Accumulators::auto_update_next_update;

  for (int i = 0; i < n_steps;) {
    /* Integrate to either the next accumulator update, or the
     * end, depending on what comes first. */
    auto const steps = std::min((n_steps - i), auto_update_next_update());

    auto const local_retval = integrate(system, steps, reuse_forces);

    // make sure all ranks exit when one rank fails
    std::remove_const_t<decltype(local_retval)> global_retval;
    boost::mpi::all_reduce(comm_cart, local_retval, global_retval,
                           std::plus<int>());
    if (global_retval < 0) {
      return global_retval; // propagate error code
    }

    reuse_forces = INTEG_REUSE_FORCES_ALWAYS;

    auto_update(comm_cart, steps);

    i += steps;
  }

  return 0;
}

double get_time_step() { return System::get_system().get_time_step(); }

double get_sim_time() { return sim_time; }

void increment_sim_time(double amount) { sim_time += amount; }

void set_time(double value) {
  ::sim_time = value;
  ::recalc_forces = true;
  auto &system = System::get_system();
  auto &box_geo = *system.box_geo;
  LeesEdwards::update_box_params(box_geo);
}

void set_integ_switch(int value) {
  ::integ_switch = value;
  ::recalc_forces = true;
}
