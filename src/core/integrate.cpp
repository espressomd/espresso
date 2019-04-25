/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "domain_decomposition.hpp"
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
#include "minimize_energy.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include "rattle.hpp"
#include "rotation.hpp"
#include "signalhandling.hpp"
#include "swimmer_reaction.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#include <profiler/profiler.hpp>
#include <utils/constants.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
double skin2 = 0.0;
bool skin_set = false;

int recalc_forces = 1;

double verlet_reuse = 0.0;

#ifdef ADDITIONAL_CHECKS
double db_max_force = 0.0, db_max_vel = 0.0;
int db_maxf_id = 0, db_maxv_id = 0;
#endif

bool set_py_interrupt = false;
namespace {
volatile std::sig_atomic_t ctrl_C = 0;
}

/** \name Private Functions */
/************************************************************/
/*@{*/

/** Propagate the velocities. Integration step 1 of the Velocity Verlet
   integrator:<br>
    \f[ v(t+0.5 \Delta t) = v(t) + 0.5 \Delta t f(t)/m \f] */
void propagate_vel();
/** Propagate the positions. Integration step 2 of the Velocity
   Verletintegrator:<br>
    \f[ p(t+\Delta t) = p(t) + \Delta t  v(t+0.5 \Delta t) \f] */
void propagate_pos();
/** Propagate the velocities and positions. Integration step 1 and 2
    of the Velocity Verlet integrator: <br>
    \f[ v(t+0.5 \Delta t) = v(t) + 0.5 \Delta t f(t)/m \f] <br>
    \f[ p(t+\Delta t) = p(t) + \Delta t  v(t+0.5 \Delta t) \f] */
void propagate_vel_pos();
/** Integration step 4 of the Velocity Verletintegrator and finalize
    instantaneous pressure calculation:<br>
    \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t)/m \f] */
void propagate_vel_finalize_p_inst();

/** Integrator stability check (see compile flag ADDITIONAL_CHECKS). */
void force_and_velocity_display();

void finalize_p_inst_npt();

/*@}*/

void integrator_sanity_checks() {
  // char *errtext;

  if (time_step < 0.0) {
    runtimeErrorMsg() << "time_step not set";
  }
}

#ifdef NPT

void integrator_npt_sanity_checks() {
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    if (nptiso.piston <= 0.0) {
      runtimeErrorMsg() << "npt on, but piston mass not set";
    }

#ifdef ELECTROSTATICS
    Coulomb::integrate_sanity_check();
#endif /*ELECTROSTATICS*/

#ifdef DIPOLES
    Dipole::integrate_sanity_check();
#endif /* ifdef DIPOLES */
  }
}
#endif /*NPT*/

/************************************************************/
void integrate_ensemble_init() {
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    /* prepare NpT-integration */
    nptiso.inv_piston = 1 / (1.0 * nptiso.piston);
    nptiso.p_inst_av = 0.0;
    if (nptiso.dimension == 0) {
      fprintf(stderr,
              "%d: INTERNAL ERROR: npt integrator was called but "
              "dimension not yet set. this should not happen. ",
              this_node);
      errexit();
    }

    nptiso.volume = pow(box_l[nptiso.non_const_dim], nptiso.dimension);

    if (recalc_forces) {
      nptiso.p_inst = 0.0;
      nptiso.p_vir[0] = nptiso.p_vir[1] = nptiso.p_vir[2] = 0.0;
      nptiso.p_vel[0] = nptiso.p_vel[1] = nptiso.p_vel[2] = 0.0;
    }
  }
#endif
}

/************************************************************/

void integrate_vv(int n_steps, int reuse_forces) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  /* Prepare the Integrator */
  on_integration_start();

#ifdef IMMERSED_BOUNDARY
  // Here we initialize volume conservation
  // This function checks if the reference volumes have been set and if
  // necessary calculates them
  immersed_boundaries.init_volume_conservation();
#endif

  /* if any method vetoes (P3M not initialized), immediately bail out */
  if (check_runtime_errors())
    return;

  /* Verlet list criterion */
  skin2 = Utils::sqr(0.5 * skin);

  INTEG_TRACE(fprintf(
      stderr, "%d: integrate_vv: integrating %d steps (recalc_forces=%d)\n",
      this_node, n_steps, recalc_forces));

  /* Integration Step: Preparation for first integration step:
     Calculate forces f(t) as function of positions p(t) ( and velocities v(t) )
     */
  /* reuse_forces logic:
     -1: recalculate forces unconditionally, mostly used for timing
      0: recalculate forces if recalc_forces is set, meaning it is probably
     necessary
      1: do not recalculate forces. Mostly when reading checkpoints with forces
   */
  if (reuse_forces == -1 || (recalc_forces && reuse_forces != 1)) {
    ESPRESSO_PROFILER_MARK_BEGIN("Initial Force Calculation");
    thermo_heat_up();

#if defined(LB) || defined(LB_GPU)
    lb_lbcoupling_deactivate();
#endif

    // Communication step: distribute ghost positions
    cells_update_ghosts();

// VIRTUAL_SITES pos (and vel for DPD) update for security reason !!!
#ifdef VIRTUAL_SITES
    virtual_sites()->update();
    if (virtual_sites()->need_ghost_comm_after_pos_update()) {
      ghost_communicator(&cell_structure.update_ghost_pos_comm);
    }
#endif

    // Langevin philox rng counter
    if (n_steps > 0) {
      langevin_rng_counter_increment();
    }

    force_calc();

    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
#ifdef ROTATION
      convert_initial_torques();
#endif
    }

    thermo_cool_down();

#ifdef COLLISION_DETECTION
    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      handle_collisions();
    }
#endif
    ESPRESSO_PROFILER_MARK_END("Initial Force Calculation");
  }

  if (check_runtime_errors())
    return;

  n_verlet_updates = 0;

#ifdef VALGRIND_INSTRUMENTATION
  CALLGRIND_START_INSTRUMENTATION;
#endif

  /* Integration loop */
  ESPRESSO_PROFILER_CXX_MARK_LOOP_BEGIN(integration_loop, "Integration loop");
  for (int step = 0; step < n_steps; step++) {
    ESPRESSO_PROFILER_CXX_MARK_LOOP_ITERATION(integration_loop, step);
    INTEG_TRACE(fprintf(stderr, "%d: STEP %d\n", this_node, step));

#ifdef BOND_CONSTRAINT
    if (n_rigidbonds)
      save_old_pos();

#endif

    /* Integration Steps: Step 1 and 2 of Velocity Verlet scheme:
       v(t+0.5*dt) = v(t) + 0.5*dt * a(t)
       p(t + dt)   = p(t) + dt * v(t+0.5*dt)
       NOTE: Depending on the integration method Step 1 and Step 2
       cannot be combined for the translation.
    */
    if (integ_switch == INTEG_METHOD_NPT_ISO) {
      propagate_vel();
      propagate_pos();

      /* Propagate time: t = t+dt */
      sim_time += time_step;
    } else if (integ_switch == INTEG_METHOD_STEEPEST_DESCENT) {
      if (steepest_descent_step())
        break;
    } else {
      propagate_vel_pos();

      /* Propagate time: t = t+dt */
      sim_time += time_step;
    }

#ifdef BOND_CONSTRAINT
    /**Correct those particle positions that participate in a rigid/constrained
     * bond */
    if (n_rigidbonds) {
      cells_update_ghosts();

      correct_pos_shake();
    }
#endif

    /* Integration Step: Step 3 of Velocity Verlet scheme:
       Calculate f(t+dt) as function of positions p(t+dt) ( and velocities
       v(t+0.5*dt) ) */

#if defined(LB) || defined(LB_GPU)
    if (n_part > 0)
      lb_lbcoupling_activate();
#endif

    // Communication step: distribute ghost positions
    cells_update_ghosts();

// VIRTUAL_SITES pos (and vel for DPD) update for security reason !!!
#ifdef VIRTUAL_SITES
    virtual_sites()->update();
    if (virtual_sites()->need_ghost_comm_after_pos_update()) {
      ghost_communicator(&cell_structure.update_ghost_pos_comm);
    }
#endif

    // Propagate langevin philox rng counter
    langevin_rng_counter_increment();

    force_calc();

#ifdef VIRTUAL_SITES
    virtual_sites()->after_force_calc();
#endif

#ifdef SWIMMER_REACTIONS
    integrate_reaction();
#endif

    /* Integration Step: Step 4 of Velocity Verlet scheme:
       v(t+dt) = v(t+0.5*dt) + 0.5*dt * f(t+dt) */
    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      propagate_vel_finalize_p_inst();
#ifdef ROTATION
      convert_torques_propagate_omega();
#endif
    }
// SHAKE velocity updates
#ifdef BOND_CONSTRAINT
    if (n_rigidbonds) {
      ghost_communicator(&cell_structure.update_ghost_pos_comm);
      correct_vel_shake();
    }
#endif

    // propagate one-step functionalities

    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
#if defined(LB) || defined(LB_GPU)
      lb_lbfluid_propagate();
      lb_lbcoupling_propagate();
#endif // LB || LB_GPU

#ifdef VIRTUAL_SITES
      virtual_sites()->after_lb_propagation();
#endif
    }

#ifdef NPT
    if ((this_node == 0) && (integ_switch == INTEG_METHOD_NPT_ISO))
      nptiso.p_inst_av += nptiso.p_inst;
#endif

    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
#ifdef COLLISION_DETECTION
      handle_collisions();
#endif
    }

    if (check_runtime_errors())
      break;

    // Check if SIGINT has been caught.
    if (ctrl_C == 1) {
      // Reset the flag.
      ctrl_C = 0;

      // Set global to notify Python of signal.
      set_py_interrupt = true;

      // Break the integration loop
      break;
    }
  }
  ESPRESSO_PROFILER_CXX_MARK_LOOP_END(integration_loop);
// VIRTUAL_SITES update vel
#ifdef VIRTUAL_SITES
  if (virtual_sites()->need_ghost_comm_before_vel_update()) {
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
  }
  virtual_sites()->update(false); // Recalc positions = false
#endif

#ifdef VALGRIND_INSTRUMENTATION
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

  /* verlet list statistics */
  if (n_verlet_updates > 0)
    verlet_reuse = n_steps / (double)n_verlet_updates;
  else
    verlet_reuse = 0;

#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    nptiso.invalidate_p_vel = 0;
    MPI_Bcast(&nptiso.p_inst, 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&nptiso.p_diff, 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&nptiso.volume, 1, MPI_DOUBLE, 0, comm_cart);
    if (this_node == 0)
      nptiso.p_inst_av /= 1.0 * n_steps;
    MPI_Bcast(&nptiso.p_inst_av, 1, MPI_DOUBLE, 0, comm_cart);
  }
#endif
}

/************************************************************/

/* Private functions */
/************************************************************/

void propagate_vel_finalize_p_inst() {
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    nptiso.p_vel[0] = nptiso.p_vel[1] = nptiso.p_vel[2] = 0.0;
  }
#endif

  INTEG_TRACE(
      fprintf(stderr, "%d: propagate_vel_finalize_p_inst:\n", this_node));

  for (auto &p : local_cells.particles()) {
    ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
        stderr, "%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",
        this_node, p.f.f[0], p.f.f[1], p.f.f[2], p.m.v[0], p.m.v[1], p.m.v[2]));
#ifdef VIRTUAL_SITES
    // Virtual sites are not propagated during integration
    if (p.p.is_virtual)
      continue;
#endif
    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
#endif
#ifdef NPT
        if (integ_switch == INTEG_METHOD_NPT_ISO &&
            (nptiso.geometry & nptiso.nptgeom_dir[j])) {
          nptiso.p_vel[j] += Utils::sqr(p.m.v[j] * time_step) * p.p.mass;
          p.m.v[j] += 0.5 * time_step / p.p.mass * p.f.f[j] +
                      friction_therm0_nptiso(p.m.v[j]) / p.p.mass;
        } else
#endif
          /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt) */
          p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;
#ifdef EXTERNAL_FORCES
      }
#endif
    }

    ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
        stderr, "%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n", this_node, p.m.v[0],
        p.m.v[1], p.m.v[2]));
  }

#ifdef NPT
  finalize_p_inst_npt();
#endif
}

void finalize_p_inst_npt() {
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    double p_tmp = 0.0;
    int i;
    /* finalize derivation of p_inst */
    nptiso.p_inst = 0.0;
    for (i = 0; i < 3; i++) {
      if (nptiso.geometry & nptiso.nptgeom_dir[i]) {
        nptiso.p_vel[i] /= Utils::sqr(time_step);
        nptiso.p_inst += nptiso.p_vir[i] + nptiso.p_vel[i];
      }
    }

    MPI_Reduce(&nptiso.p_inst, &p_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
    if (this_node == 0) {
      nptiso.p_inst = p_tmp / (nptiso.dimension * nptiso.volume);
      nptiso.p_diff = nptiso.p_diff +
                      (nptiso.p_inst - nptiso.p_ext) * 0.5 * time_step +
                      friction_thermV_nptiso(nptiso.p_diff);
    }
  }
#endif
}

void propagate_press_box_pos_and_rescale_npt() {
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    double scal[3] = {0., 0., 0.}, L_new = 0.0;

    /* finalize derivation of p_inst */
    finalize_p_inst_npt();

    /* adjust \ref nptiso_struct::nptiso.volume; prepare pos- and
     * vel-rescaling
     */
    if (this_node == 0) {
      nptiso.volume += nptiso.inv_piston * nptiso.p_diff * 0.5 * time_step;
      scal[2] = Utils::sqr(box_l[nptiso.non_const_dim]) /
                pow(nptiso.volume, 2.0 / nptiso.dimension);
      nptiso.volume += nptiso.inv_piston * nptiso.p_diff * 0.5 * time_step;
      if (nptiso.volume < 0.0) {

        runtimeErrorMsg()
            << "your choice of piston= " << nptiso.piston
            << ", dt= " << time_step << ", p_diff= " << nptiso.p_diff
            << " just caused the volume to become negative, decrease dt";
        nptiso.volume = box_l[0] * box_l[1] * box_l[2];
        scal[2] = 1;
      }

      L_new = pow(nptiso.volume, 1.0 / nptiso.dimension);
      // printf("Lnew, %f: volume, %f: dim, %f: press, %f \n", L_new,
      // nptiso.volume, nptiso.dimension,nptiso.p_inst );
      // fflush(stdout);

      scal[1] = L_new / box_l[nptiso.non_const_dim];
      scal[0] = 1 / scal[1];
    }
    MPI_Bcast(scal, 3, MPI_DOUBLE, 0, comm_cart);

    /* propagate positions while rescaling positions and velocities */
    for (auto &p : local_cells.particles()) {
#ifdef VIRTUAL_SITES
      if (p.p.is_virtual)
        continue;
#endif
      for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
        if (!(p.p.ext_flag & COORD_FIXED(j))) {
#endif
          if (nptiso.geometry & nptiso.nptgeom_dir[j]) {
            {
              p.r.p[j] = scal[1] * (p.r.p[j] + scal[2] * p.m.v[j] * time_step);
              p.l.p_old[j] *= scal[1];
              p.m.v[j] *= scal[0];
            }
          } else {
            p.r.p[j] += p.m.v[j] * time_step;
          }
#ifdef EXTERNAL_FORCES
        }
#endif
      }
      ONEPART_TRACE(if (p.p.identity == check_id)
                        fprintf(stderr, "%d: OPT:PV_1 v_new=(%.3e,%.3e,%.3e)\n",
                                this_node, p.m.v[0], p.m.v[1], p.m.v[2]));
      ONEPART_TRACE(if (p.p.identity == check_id)
                        fprintf(stderr, "%d: OPT:PPOS p=(%.3f,%.3f,%.3f)\n",
                                this_node, p.r.p[0], p.r.p[1], p.r.p[2]));
    }

    set_resort_particles(Cells::RESORT_LOCAL);

    /* Apply new volume to the box-length, communicate it, and account for
     * necessary adjustments to the cell geometry */
    if (this_node == 0) {
      for (int i = 0; i < 3; i++) {
        if (nptiso.geometry & nptiso.nptgeom_dir[i]) {
          box_l[i] = L_new;
        } else if (nptiso.cubic_box) {
          box_l[i] = L_new;
        }
      }
    }
    MPI_Bcast(box_l.data(), 3, MPI_DOUBLE, 0, comm_cart);

    /* fast box length update */
    grid_changed_box_l();
    recalc_maximal_cutoff();
    cells_on_geometry_change(CELL_FLAG_FAST);
  }
#endif
}

void propagate_vel() {
#ifdef NPT
  nptiso.p_vel[0] = nptiso.p_vel[1] = nptiso.p_vel[2] = 0.0;
#endif

  INTEG_TRACE(fprintf(stderr, "%d: propagate_vel:\n", this_node));

  for (auto &p : local_cells.particles()) {
#ifdef ROTATION
    propagate_omega_quat_particle(&p);
#endif

// Don't propagate translational degrees of freedom of vs
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      continue;
#endif
    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
#ifdef NPT
        if (integ_switch == INTEG_METHOD_NPT_ISO &&
            (nptiso.geometry & nptiso.nptgeom_dir[j])) {
          p.m.v[j] += p.f.f[j] * 0.5 * time_step / p.p.mass +
                      friction_therm0_nptiso(p.m.v[j]) / p.p.mass;
          nptiso.p_vel[j] += Utils::sqr(p.m.v[j] * time_step) * p.p.mass;
        } else
#endif
          /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * a(t) */
          p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;
      }

      ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
          stderr, "%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n", this_node,
          p.m.v[0], p.m.v[1], p.m.v[2]));
    }
  }

#ifdef ADDITIONAL_CHECKS
  force_and_velocity_display();
#endif
}

void propagate_pos() {
  INTEG_TRACE(fprintf(stderr, "%d: propagate_pos:\n", this_node));
  if (integ_switch == INTEG_METHOD_NPT_ISO)
    /* Special propagator for NPT ISOTROPIC */
    /* Propagate pressure, box_length (2 times) and positions, rescale
       positions and velocities and check Verlet list criterion (only NPT) */
    propagate_press_box_pos_and_rescale_npt();
  else {
    for (auto &p : local_cells.particles()) {
#ifdef VIRTUAL_SITES
      if (p.p.is_virtual)
        continue;
#endif
      for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
        if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
        {
          /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt *
           * v(t+0.5*dt) */
          p.r.p[j] += time_step * p.m.v[j];
        }
      }
      /* Verlet criterion check */
      if ((p.r.p - p.l.p_old).norm2() > skin2)
        set_resort_particles(Cells::RESORT_LOCAL);
    }
  }
  announce_resort_particles();
}

void propagate_vel_pos() {
  INTEG_TRACE(fprintf(stderr, "%d: propagate_vel_pos:\n", this_node));

#ifdef ADDITIONAL_CHECKS
  db_max_force = db_max_vel = 0;
  db_maxf_id = db_maxv_id = -1;
#endif

  for (auto &p : local_cells.particles()) {
#ifdef ROTATION
    propagate_omega_quat_particle(&p);
#endif

// Don't propagate translational degrees of freedom of vs
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      continue;
#endif
    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5 * dt * a(t) */
        p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;

        /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt *
         * v(t+0.5*dt) */
        p.r.p[j] += time_step * p.m.v[j];
      }
    }

    ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
        stderr, "%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n", this_node, p.m.v[0],
        p.m.v[1], p.m.v[2]));
    ONEPART_TRACE(if (p.p.identity == check_id)
                      fprintf(stderr, "%d: OPT: PPOS p = (%.3e,%.3e,%.3e)\n",
                              this_node, p.r.p[0], p.r.p[1], p.r.p[2]));

    /* Verlet criterion check*/
    if (Utils::sqr(p.r.p[0] - p.l.p_old[0]) +
            Utils::sqr(p.r.p[1] - p.l.p_old[1]) +
            Utils::sqr(p.r.p[2] - p.l.p_old[2]) >
        skin2)
      set_resort_particles(Cells::RESORT_LOCAL);
  }

  announce_resort_particles();

#ifdef ADDITIONAL_CHECKS
  force_and_velocity_display();
#endif
}

void force_and_velocity_display() {
#ifdef ADDITIONAL_CHECKS
  if (db_max_force > skin2)
    fprintf(stderr, "%d: max_force=%e, part=%d f=(%e,%e,%e)\n", this_node,
            sqrt(db_max_force), db_maxf_id, local_particles[db_maxf_id]->f.f[0],
            local_particles[db_maxf_id]->f.f[1],
            local_particles[db_maxf_id]->f.f[2]);
  if (db_max_vel > skin2)
    fprintf(stderr, "%d: max_vel=%e, part=%d v=(%e,%e,%e)\n", this_node,
            sqrt(db_max_vel), db_maxv_id, local_particles[db_maxv_id]->m.v[0],
            local_particles[db_maxv_id]->m.v[1],
            local_particles[db_maxv_id]->m.v[2]);
#endif
}

/** @todo This needs to go!! */

int python_integrate(int n_steps, bool recalc_forces, bool reuse_forces_par) {
  // Override the signal handler so that the integrator obeys Ctrl+C
  SignalHandler sa(SIGINT, [](int) { ctrl_C = 1; });

  int reuse_forces = 0;
  reuse_forces = reuse_forces_par;
  INTEG_TRACE(fprintf(stderr, "%d: integrate:\n", this_node));

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
    if (max_cut == 0.0) {
      runtimeErrorMsg()
          << "cannot automatically determine skin, please set it manually";
      return ES_ERROR;
    }
    skin = std::min(0.4 * max_cut, max_skin);
    mpi_bcast_parameter(FIELD_SKIN);
  }

  /* perform integration */
  if (!Accumulators::auto_update_enabled()) {
    if (mpi_integrate(n_steps, reuse_forces))
      return ES_ERROR;
  } else {
    for (int i = 0; i < n_steps; i++) {
      if (mpi_integrate(1, reuse_forces))
        return ES_ERROR;
      reuse_forces = 1;
      Accumulators::auto_update();
    }
    if (n_steps == 0) {
      if (mpi_integrate(0, reuse_forces))
        return ES_ERROR;
    }
  }
  return ES_OK;
}

void integrate_set_nvt() {
  integ_switch = INTEG_METHOD_NVT;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
}

/** Parse integrate npt_isotropic command */
int integrate_set_npt_isotropic(double ext_pressure, double piston, int xdir,
                                int ydir, int zdir, bool cubic_box) {
  nptiso.cubic_box = 0;
  nptiso.p_ext = ext_pressure;
  nptiso.piston = piston;

  if (nptiso.piston <= 0.0) {
    runtimeErrorMsg() << "You must set <piston> as well before you can use "
                         "this integrator!\n";
    return ES_ERROR;
  }
  if (xdir || ydir || zdir) {
    /* set the geometry to include rescaling specified directions only*/
    nptiso.geometry = 0;
    nptiso.dimension = 0;
    nptiso.non_const_dim = -1;
    if (xdir) {
      nptiso.geometry = (nptiso.geometry | NPTGEOM_XDIR);
      nptiso.dimension += 1;
      nptiso.non_const_dim = 0;
    }
    if (ydir) {
      nptiso.geometry = (nptiso.geometry | NPTGEOM_YDIR);
      nptiso.dimension += 1;
      nptiso.non_const_dim = 1;
    }
    if (zdir) {
      nptiso.geometry = (nptiso.geometry | NPTGEOM_ZDIR);
      nptiso.dimension += 1;
      nptiso.non_const_dim = 2;
    }
  } else {
    /* set the geometry to include rescaling in all directions; the default*/
    nptiso.geometry = 0;
    nptiso.geometry = (nptiso.geometry | NPTGEOM_XDIR);
    nptiso.geometry = (nptiso.geometry | NPTGEOM_YDIR);
    nptiso.geometry = (nptiso.geometry | NPTGEOM_ZDIR);
    nptiso.dimension = 3;
    nptiso.non_const_dim = 2;
  }

  if (cubic_box) {
    /* enable if the volume fluctuations should also apply to dimensions which
   are switched off by the above flags
   and which do not contribute to the pressure (3D) / tension (2D, 1D) */
    nptiso.cubic_box = 1;
  }

/* Sanity Checks */
#ifdef ELECTROSTATICS
  if (nptiso.dimension < 3 && !nptiso.cubic_box && coulomb.prefactor > 0) {
    runtimeErrorMsg() << "WARNING: If electrostatics is being used you must "
                         "use the the cubic box npt.";
    integ_switch = INTEG_METHOD_NVT;
    mpi_bcast_parameter(FIELD_INTEG_SWITCH);
    return ES_ERROR;
  }
#endif

#ifdef DIPOLES
  if (nptiso.dimension < 3 && !nptiso.cubic_box && dipole.prefactor > 0) {
    runtimeErrorMsg() << "WARNING: If magnetostatics is being used you must "
                         "use the the cubic box npt.";
    integ_switch = INTEG_METHOD_NVT;
    mpi_bcast_parameter(FIELD_INTEG_SWITCH);
    return ES_ERROR;
  }
#endif

  if (nptiso.dimension == 0 || nptiso.non_const_dim == -1) {
    runtimeErrorMsg() << "You must enable at least one of the x y z components "
                         "as fluctuating dimension(s) for box length motion!";
    runtimeErrorMsg() << "Cannot proceed with npt_isotropic, reverting to nvt "
                         "integration... \n";
    integ_switch = INTEG_METHOD_NVT;
    mpi_bcast_parameter(FIELD_INTEG_SWITCH);
    return (ES_ERROR);
  }

  /* set integrator switch */
  integ_switch = INTEG_METHOD_NPT_ISO;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
  mpi_bcast_parameter(FIELD_NPTISO_PISTON);
  mpi_bcast_parameter(FIELD_NPTISO_PEXT);

  /* broadcast npt geometry information to all nodes */
  mpi_bcast_nptiso_geom();
  return (ES_OK);
}
