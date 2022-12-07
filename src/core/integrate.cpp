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

#include "ParticleRange.hpp"
#include "accumulators.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "bonded_interactions/rigid_bond.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/ek_container.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "interactions.hpp"
#include "lees_edwards/lees_edwards.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "rattle.hpp"
#include "rotation.hpp"
#include "signalhandling.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#include <profiler/profiler.hpp>

#include <boost/mpi/collectives/reduce.hpp>
#include <boost/range/algorithm/min_element.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <csignal>
#include <functional>
#include <stdexcept>
#include <utility>

#ifdef VALGRIND_MARKERS
#include <callgrind.h>
#endif

#ifdef WALBERLA
#ifdef WALBERLA_STATIC_ASSERT
#error "waLberla headers should not be visible to the ESPResSo core"
#endif
#endif

int integ_switch = INTEG_METHOD_NVT;

/** Time step for the integration. */
static double time_step = -1.0;

/** Actual simulation time. */
static double sim_time = 0.0;

double skin = 0.0;

/** True iff the user has changed the skin setting. */
static bool skin_set = false;

bool recalc_forces = true;

/** Average number of integration steps the Verlet list has been re-using. */
static double verlet_reuse = 0.0;

static int fluid_step = 0;
static int ek_step = 0;

bool set_py_interrupt = false;
namespace {
volatile std::sig_atomic_t ctrl_C = 0;

void notify_sig_int() {
  ctrl_C = 0;              // reset
  set_py_interrupt = true; // global to notify Python
}
} // namespace

namespace LeesEdwards {
/** @brief Currently active Lees-Edwards protocol. */
static std::shared_ptr<ActiveProtocol> protocol = nullptr;

std::weak_ptr<ActiveProtocol> get_protocol() { return protocol; }

/**
 * @brief Update the Lees-Edwards parameters of the box geometry
 * for the current simulation time.
 */
static void update_box_params() {
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    assert(protocol != nullptr);
    box_geo.lees_edwards_update(get_pos_offset(sim_time, *protocol),
                                get_shear_velocity(sim_time, *protocol));
  }
}

void set_protocol(std::shared_ptr<ActiveProtocol> new_protocol) {
  box_geo.set_type(BoxType::LEES_EDWARDS);
  protocol = std::move(new_protocol);
  LeesEdwards::update_box_params();
  ::recalc_forces = true;
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
}

void unset_protocol() {
  protocol = nullptr;
  box_geo.set_type(BoxType::CUBOID);
  ::recalc_forces = true;
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
}

template <class Kernel> void run_kernel() {
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    auto const kernel = Kernel{box_geo};
    auto const particles = cell_structure.local_particles();
    std::for_each(particles.begin(), particles.end(),
                  [&kernel](auto &p) { kernel(p); });
  }
}
} // namespace LeesEdwards

void integrator_sanity_checks() {
  if (time_step < 0.0) {
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
    if (box_geo.type() == BoxType::LEES_EDWARDS)
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

static void resort_particles_if_needed(ParticleRange const &particles) {
  auto const offset = LeesEdwards::verlet_list_offset(
      box_geo, cell_structure.get_le_pos_offset_at_last_resort());
  if (cell_structure.check_resort_required(particles, skin, offset)) {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

/** @brief Calls the hook for propagation kernels before the force calculation
 *  @return whether or not to stop the integration loop early.
 */
static bool integrator_step_1(ParticleRange const &particles) {
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
static void integrator_step_2(ParticleRange const &particles, double kT) {
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

int integrate(int n_steps, int reuse_forces) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  // Prepare particle structure and run sanity checks of all active algorithms
  on_integration_start(time_step);

  // If any method vetoes (e.g. P3M not initialized), immediately bail out
  if (check_runtime_errors(comm_cart))
    return 0;

  // Additional preparations for the first integration step
  if (reuse_forces == -1 || (recalc_forces && reuse_forces != 1)) {
    ESPRESSO_PROFILER_MARK_BEGIN("Initial Force Calculation");
    lb_lbcoupling_deactivate();

#ifdef VIRTUAL_SITES
    virtual_sites()->update();
#endif

    // Communication step: distribute ghost positions
    cells_update_ghosts(global_ghost_flags());

    force_calc(cell_structure, time_step, temperature);

    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
#ifdef ROTATION
      convert_initial_torques(cell_structure.local_particles());
#endif
    }

    ESPRESSO_PROFILER_MARK_END("Initial Force Calculation");
  }

  lb_lbcoupling_activate();

  if (check_runtime_errors(comm_cart))
    return 0;

  // Keep track of the number of Verlet updates (i.e. particle resorts)
  int n_verlet_updates = 0;

#ifdef VALGRIND_MARKERS
  CALLGRIND_START_INSTRUMENTATION;
#endif
  // Integration loop
  ESPRESSO_PROFILER_CXX_MARK_LOOP_BEGIN(integration_loop, "Integration loop");
  int integrated_steps = 0;
  for (int step = 0; step < n_steps; step++) {
    ESPRESSO_PROFILER_CXX_MARK_LOOP_ITERATION(integration_loop, step);

    auto particles = cell_structure.local_particles();

#ifdef BOND_CONSTRAINT
    if (n_rigidbonds)
      save_old_position(particles, cell_structure.ghost_particles());
#endif

    LeesEdwards::update_box_params();
    bool early_exit = integrator_step_1(particles);
    if (early_exit)
      break;

    LeesEdwards::run_kernel<LeesEdwards::Push>();

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
      correct_position_shake(cell_structure);
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

    force_calc(cell_structure, time_step, temperature);

#ifdef VIRTUAL_SITES
    virtual_sites()->after_force_calc(time_step);
#endif
    integrator_step_2(particles, temperature);
    LeesEdwards::run_kernel<LeesEdwards::UpdateOffset>();
#ifdef BOND_CONSTRAINT
    // SHAKE velocity updates
    if (n_rigidbonds) {
      correct_velocity_shake(cell_structure);
    }
#endif

    // propagate one-step functionalities
    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      auto const lb_active = LB::get_lattice_switch() != ActiveLB::NONE;
#ifdef WALBERLA
      auto const ek_active = not EK::ek_container.empty();
#else
      auto constexpr ek_active = false;
#endif

      if (lb_active and ek_active) {
        // assume that they are coupled, which is not necessarily true
        auto const lb_steps_per_md_step = LB::get_steps_per_md_step(time_step);
        auto const ek_steps_per_md_step = EK::get_steps_per_md_step(time_step);

        if (lb_steps_per_md_step != ek_steps_per_md_step) {
          runtimeErrorMsg()
              << "LB and EK are active but with different timesteps.";
        }

        // only use fluid_step in this case
        assert(fluid_step == ek_step);

        fluid_step += 1;
        if (fluid_step >= lb_steps_per_md_step) {
          fluid_step = 0;
          LB::propagate();
          EK::propagate();
        }
        lb_lbcoupling_propagate();
      } else if (lb_active) {
        auto const lb_steps_per_md_step = LB::get_steps_per_md_step(time_step);
        fluid_step += 1;
        if (fluid_step >= lb_steps_per_md_step) {
          fluid_step = 0;
          LB::propagate();
        }
        lb_lbcoupling_propagate();
      } else if (ek_active) {
        auto const ek_steps_per_md_step = EK::get_steps_per_md_step(time_step);
        ek_step += 1;
        if (ek_step >= ek_steps_per_md_step) {
          ek_step = 0;
          EK::propagate();
        }
      }

#ifdef VIRTUAL_SITES
      virtual_sites()->after_lb_propagation(time_step);
#endif

#ifdef COLLISION_DETECTION
      handle_collisions();
#endif
      BondBreakage::process_queue();
    }

    integrated_steps++;

    if (check_runtime_errors(comm_cart))
      break;

    // Check if SIGINT has been caught.
    if (ctrl_C == 1) {
      notify_sig_int();
      break;
    }

  } // for-loop over integration steps
  LeesEdwards::update_box_params();
  ESPRESSO_PROFILER_CXX_MARK_LOOP_END(integration_loop);

#ifdef VALGRIND_MARKERS
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

#ifdef VIRTUAL_SITES
  virtual_sites()->update();
#endif

  // Verlet list statistics
  if (n_verlet_updates > 0)
    verlet_reuse = n_steps / static_cast<double>(n_verlet_updates);
  else
    verlet_reuse = 0;

#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    synchronize_npt_state();
  }
#endif
  return integrated_steps;
}

int python_integrate(int n_steps, bool recalc_forces_par,
                     bool reuse_forces_par) {

  assert(n_steps >= 0);

  // Override the signal handler so that the integrator obeys Ctrl+C
  SignalHandler sa(SIGINT, [](int) { ctrl_C = 1; });

  int reuse_forces = reuse_forces_par;

  if (recalc_forces_par) {
    if (reuse_forces) {
      runtimeErrorMsg() << "cannot reuse old forces and recalculate forces";
    }
    reuse_forces = -1;
  }

  /* if skin wasn't set, do an educated guess now */
  if (!skin_set) {
    auto const max_cut = maximal_cutoff(n_nodes);
    if (max_cut <= 0.0) {
      runtimeErrorMsg()
          << "cannot automatically determine skin, please set it manually";
      return ES_ERROR;
    }
    /* maximal skin that can be used without resorting is the maximal
     * range of the cell system minus what is needed for interactions. */
    auto const new_skin =
        std::min(0.4 * max_cut,
                 *boost::min_element(cell_structure.max_cutoff()) - max_cut);
    mpi_call_all(mpi_set_skin_local, new_skin);
  }

  using Accumulators::auto_update;
  using Accumulators::auto_update_next_update;

  for (int i = 0; i < n_steps;) {
    /* Integrate to either the next accumulator update, or the
     * end, depending on what comes first. */
    auto const steps = std::min((n_steps - i), auto_update_next_update());
    if (mpi_integrate(steps, reuse_forces))
      return ES_ERROR;

    reuse_forces = 1;

    auto_update(steps);

    i += steps;
  }

  if (n_steps == 0) {
    if (mpi_integrate(0, reuse_forces))
      return ES_ERROR;
  }

  return ES_OK;
}

static int mpi_steepest_descent_local(int steps) {
  return integrate(steps, -1);
}

REGISTER_CALLBACK_MAIN_RANK(mpi_steepest_descent_local)

int mpi_steepest_descent(int steps) {
  return mpi_call(Communication::Result::main_rank, mpi_steepest_descent_local,
                  steps);
}

static int mpi_integrate_local(int n_steps, int reuse_forces) {
  integrate(n_steps, reuse_forces);

  return check_runtime_errors_local();
}

REGISTER_CALLBACK_REDUCTION(mpi_integrate_local, std::plus<int>())

int mpi_integrate(int n_steps, int reuse_forces) {
  return mpi_call(Communication::Result::reduction, std::plus<int>(),
                  mpi_integrate_local, n_steps, reuse_forces);
}

double interaction_range() {
  /* Consider skin only if there are actually interactions */
  auto const max_cut = maximal_cutoff(n_nodes == 1);
  return (max_cut > 0.) ? max_cut + skin : INACTIVE_CUTOFF;
}

double get_verlet_reuse() { return verlet_reuse; }

double get_time_step() { return time_step; }

double get_sim_time() { return sim_time; }

void increment_sim_time(double amount) { sim_time += amount; }

void set_time_step(double value) {
  if (value <= 0.)
    throw std::domain_error("time_step must be > 0.");
  if (LB::get_lattice_switch() != ActiveLB::NONE) {
    LB::check_tau_time_step_consistency(LB::get_tau(), value);
  }
  ::time_step = value;
  on_timestep_change();
}

void mpi_set_skin_local(double value) {
  ::skin = value;
  skin_set = true;
  on_skin_change();
}

REGISTER_CALLBACK(mpi_set_skin_local)

void set_time(double value) {
  ::sim_time = value;
  ::recalc_forces = true;
  LeesEdwards::update_box_params();
}

void set_integ_switch(int value) {
  ::integ_switch = value;
  ::recalc_forces = true;
}
