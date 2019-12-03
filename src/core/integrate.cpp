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
 *  Molecular dynamics integrator.
 *
 *  For more information about the integrator
 *  see \ref integrate.hpp "integrate.hpp".
 */

#include "integrate.hpp"
#include "accumulators.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "domain_decomposition.hpp"
#include "dpd.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "ghosts.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "immersed_boundaries.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include "rattle.hpp"
#include "rotation.hpp"
#include "signalhandling.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#include "integrators/steepest_descent.hpp"
#include "integrators/velocity_verlet_inline.hpp"
#include "integrators/velocity_verlet_npt.hpp"

#include <profiler/profiler.hpp>
#include <utils/constants.hpp>

#include <boost/range/algorithm/min_element.hpp>
#include <cmath>
#include <cstdio>
#include <mpi.h>

#ifdef VALGRIND_INSTRUMENTATION
#include <callgrind.h>
#endif

/*******************  variables  *******************/

int integ_switch = INTEG_METHOD_NVT;

int n_verlet_updates = 0;

double time_step = -1.0;
double time_step_half = -1.0;
double time_step_squared = -1.0;
double time_step_squared_half = -1.0;

double sim_time = 0.0;
double skin = 0.0;
bool skin_set = false;

bool recalc_forces = true;

double verlet_reuse = 0.0;

bool set_py_interrupt = false;
namespace {
volatile std::sig_atomic_t ctrl_C = 0;

void notify_sig_int() {
  ctrl_C = 0;              // reset
  set_py_interrupt = true; // global to notify Python
}
} // namespace

/** \name Private Functions */
/************************************************************/
/*@{*/

/** Thermostats increment the RNG counter here. */
void philox_counter_increment();

/*@}*/

void integrator_sanity_checks() {
  if (time_step < 0.0) {
    runtimeErrorMsg() << "time_step not set";
  }
}

/************************************************************/

/** @brief Calls the hook for propagation kernels before the force calculation
 *  @return whether or not to stop the integration loop early.
 */
bool integrator_step_1(ParticleRange &particles) {
  switch (integ_switch) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    if (steepest_descent_step(particles))
      return true; // early exit
    break;
  case INTEG_METHOD_NVT:
    velocity_verlet_step_1(particles);
    break;
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    velocity_verlet_npt_step_1(particles);
    break;
#endif
  default:
    throw std::runtime_error("Unknown value for integ_switch");
  }
  return false;
}

/** Calls the hook of the propagation kernels after force calculation */
void integrator_step_2(ParticleRange &particles) {
  switch (integ_switch) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    // Nothing
    break;
  case INTEG_METHOD_NVT:
    velocity_verlet_step_2(particles);
    break;
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    velocity_verlet_npt_step_2(particles);
    break;
#endif
  default:
    throw std::runtime_error("Unknown value for INTEG_SWITCH");
  }
}

void integrate_vv(int n_steps, int reuse_forces) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  /* Prepare the integrator */
  on_integration_start();

  /* if any method vetoes (P3M not initialized), immediately bail out */
  if (check_runtime_errors(comm_cart))
    return;

  /* Verlet list criterion */

  /* Integration Step: Preparation for first integration step:
   * Calculate forces F(t) as function of positions x(t) (and velocities v(t))
   */
  if (reuse_forces == -1 || (recalc_forces && reuse_forces != 1)) {
    ESPRESSO_PROFILER_MARK_BEGIN("Initial Force Calculation");
    lb_lbcoupling_deactivate();

#ifdef VIRTUAL_SITES
    virtual_sites()->update();
#endif

    // Communication step: distribute ghost positions
    cells_update_ghosts();

    force_calc(cell_structure);

    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
#ifdef ROTATION
      convert_initial_torques(cell_structure.local_cells().particles());
#endif
    }

    ESPRESSO_PROFILER_MARK_END("Initial Force Calculation");
  }

  if (n_part > 0) {
    lb_lbcoupling_activate();
  }

  if (check_runtime_errors(comm_cart))
    return;

  n_verlet_updates = 0;

#ifdef VALGRIND_INSTRUMENTATION
  CALLGRIND_START_INSTRUMENTATION;
#endif

  /* Integration loop */
  ESPRESSO_PROFILER_CXX_MARK_LOOP_BEGIN(integration_loop, "Integration loop");
  for (int step = 0; step < n_steps; step++) {
    ESPRESSO_PROFILER_CXX_MARK_LOOP_ITERATION(integration_loop, step);

    auto particles = cell_structure.local_cells().particles();

#ifdef BOND_CONSTRAINT
    if (n_rigidbonds)
      save_old_pos(particles, ghost_cells.particles());
#endif

    bool early_exit = integrator_step_1(particles);
    if (early_exit)
      break;

    /* Propagate philox rng counters */
    philox_counter_increment();

#ifdef BOND_CONSTRAINT
    /* Correct those particle positions that participate in a rigid/constrained
     * bond */
    if (n_rigidbonds) {
      correct_pos_shake(cell_structure.local_cells().particles());
    }
#endif

#ifdef VIRTUAL_SITES
    // VIRTUAL_SITES pos (and vel for DPD) update for security reason!!!
    virtual_sites()->update();
#endif

    // Communication step: distribute ghost positions
    cells_update_ghosts();

    particles = cell_structure.local_cells().particles();

    force_calc(cell_structure);

#ifdef VIRTUAL_SITES
    virtual_sites()->after_force_calc();
#endif
    integrator_step_2(particles);
#ifdef BOND_CONSTRAINT
    // SHAKE velocity updates
    if (n_rigidbonds) {
      correct_vel_shake();
    }
#endif

    // propagate one-step functionalities
    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      lb_lbfluid_propagate();
      lb_lbcoupling_propagate();

#ifdef VIRTUAL_SITES
      virtual_sites()->after_lb_propagation();
#endif

#ifdef COLLISION_DETECTION
      handle_collisions();
#endif
    }

    if (check_runtime_errors(comm_cart))
      break;

    // Check if SIGINT has been caught.
    if (ctrl_C == 1) {
      notify_sig_int();
      break;
    }

  } // for-loop over integration steps
  ESPRESSO_PROFILER_CXX_MARK_LOOP_END(integration_loop);
#ifdef VALGRIND_INSTRUMENTATION
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

#ifdef VIRTUAL_SITES
  // VIRTUAL_SITES update vel
  virtual_sites()->update(false); // Recalc positions = false
#endif

  /* verlet list statistics */
  if (n_verlet_updates > 0)
    verlet_reuse = n_steps / (double)n_verlet_updates;
  else
    verlet_reuse = 0;

#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    synchronize_npt_state(n_steps);
  }
#endif
}

/************************************************************/

/* Private functions */
/************************************************************/

void philox_counter_increment() {
  if (thermo_switch & THERMO_LANGEVIN) {
    langevin_rng_counter_increment();
  }
  if (thermo_switch & THERMO_DPD) {
#ifdef DPD
    dpd_rng_counter_increment();
#endif
  }
  if (n_thermalized_bonds)
    thermalized_bond_rng_counter_increment();
}

int python_integrate(int n_steps, bool recalc_forces, bool reuse_forces_par) {
  // Override the signal handler so that the integrator obeys Ctrl+C
  SignalHandler sa(SIGINT, [](int) { ctrl_C = 1; });

  int reuse_forces = reuse_forces_par;

  if (recalc_forces) {
    if (reuse_forces) {
      runtimeErrorMsg() << "cannot reuse old forces and recalculate forces";
    }
    reuse_forces = -1;
  }

  /* go on with integrate <n_steps> */
  if (n_steps < 0) {
    runtimeErrorMsg() << "illegal number of steps (must be >0)";
    return ES_ERROR;
  }

  /* if skin wasn't set, do an educated guess now */
  if (!skin_set) {
    if (max_cut <= 0.0) {
      runtimeErrorMsg()
          << "cannot automatically determine skin, please set it manually";
      return ES_ERROR;
    }
    /* maximal skin that can be used without resorting is the maximal
     * range of the cell system minus what is needed for interactions. */
    skin = std::min(0.4 * max_cut,
                    *boost::min_element(cell_structure.max_range) - max_cut);
    mpi_bcast_parameter(FIELD_SKIN);
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

void integrate_set_nvt() {
  integ_switch = INTEG_METHOD_NVT;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
}

#ifdef NPT
int integrate_set_npt_isotropic(double ext_pressure, double piston,
                                bool xdir_rescale, bool ydir_rescale,
                                bool zdir_rescale, bool cubic_box) {
  nptiso.cubic_box = cubic_box;
  nptiso.p_ext = ext_pressure;
  nptiso.piston = piston;

  if (nptiso.piston <= 0.0) {
    runtimeErrorMsg() << "You must set <piston> as well before you can use "
                         "this integrator!\n";
    return ES_ERROR;
  }
  /* set the NpT geometry */
  nptiso.geometry = 0;
  nptiso.dimension = 0;
  nptiso.non_const_dim = -1;
  if (xdir_rescale) {
    nptiso.geometry |= NPTGEOM_XDIR;
    nptiso.dimension += 1;
    nptiso.non_const_dim = 0;
  }
  if (ydir_rescale) {
    nptiso.geometry |= NPTGEOM_YDIR;
    nptiso.dimension += 1;
    nptiso.non_const_dim = 1;
  }
  if (zdir_rescale) {
    nptiso.geometry |= NPTGEOM_ZDIR;
    nptiso.dimension += 1;
    nptiso.non_const_dim = 2;
  }

  /* Sanity Checks */
#ifdef ELECTROSTATICS
  if (nptiso.dimension < 3 && !nptiso.cubic_box && coulomb.prefactor > 0) {
    runtimeErrorMsg() << "WARNING: If electrostatics is being used you must "
                         "use the cubic box npt.";
    return ES_ERROR;
  }
#endif

#ifdef DIPOLES
  if (nptiso.dimension < 3 && !nptiso.cubic_box && dipole.prefactor > 0) {
    runtimeErrorMsg() << "WARNING: If magnetostatics is being used you must "
                         "use the cubic box npt.";
    return ES_ERROR;
  }
#endif

  if (nptiso.dimension == 0 || nptiso.non_const_dim == -1) {
    runtimeErrorMsg() << "You must enable at least one of the x y z components "
                         "as fluctuating dimension(s) for box length motion!";
    return ES_ERROR;
  }

  /* set integrator switch */
  integ_switch = INTEG_METHOD_NPT_ISO;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
  mpi_bcast_parameter(FIELD_NPTISO_PISTON);
  mpi_bcast_parameter(FIELD_NPTISO_PEXT);

  /* broadcast NpT geometry information to all nodes */
  mpi_bcast_nptiso_geom();
  return ES_OK;
}
#endif
