/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

/** \file integrate.cpp   Molecular dynamics integrator.
 *
 *  For more information about the integrator 
 *  see \ref integrate.hpp "integrate.hpp".
*/

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "lees_edwards.hpp"
#include "utils.hpp"
#include "integrate.hpp"
#include "reaction.hpp"
#include "electrokinetics.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "cells.hpp"
#include "verlet.hpp"
#include "rotation.hpp"
#include "ghosts.hpp"
#include "pressure.hpp"
#include "p3m.hpp"
#include "maggs.hpp"
#include "thermostat.hpp"
#include "initialize.hpp"
#include "forces.hpp"
#include "nsquare.hpp"
#include "domain_decomposition.hpp"
#include "layered.hpp"
#include "nemd.hpp"
#include "rattle.hpp"
#include "errorhandling.hpp"
#include "lattice.hpp"
#include "lb.hpp"
#include "virtual_sites.hpp"
#include "statistics_correlation.hpp"
#include "ghmc.hpp"
#include "immersed_boundary/ibm_main.hpp"
#include "immersed_boundary/ibm_volume_conservation.hpp"
#include "minimize_energy.hpp"

/************************************************
 * DEFINES
 ************************************************/

/** Tag for communication in verlet fix: propagate_positions()  */
#define REQ_INT_VERLET   400

/*******************  variables  *******************/

int    integ_switch     = INTEG_METHOD_NVT;

int n_verlet_updates    = 0;

double time_step        = -1.0;
double time_step_half   = -1.0;
double time_step_squared= -1.0;
double time_step_squared_half = -1.0;

double sim_time         = 0.0;
double skin             = 0.0;
double skin2            = 0.0;
bool   skin_set         = false;

int    resort_particles = 1;
int    recalc_forces    = 1;

double verlet_reuse     = 0.0;

double smaller_time_step          = -1.0;
#ifdef MULTI_TIMESTEP
int    current_time_step_is_small = 0;
int    mts_index                  = 0;
int    mts_max                    = 0;
#ifdef NPT
double scal_store[3]              = {0.,0.,0.};
double virial_store[3]            = {0.,0.,0.};
#endif
#endif

/** For configurational temperature only */
double configtemp[2]              = {0.,0.};

#ifdef ADDITIONAL_CHECKS
double db_max_force = 0.0, db_max_vel = 0.0;
int    db_maxf_id   = 0,   db_maxv_id = 0;
#endif

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Rescale all particle forces with \f[ 0.5 \Delta t^2 \f]. */
void rescale_forces();
/** Propagate the velocities. Integration step 1 of the Velocity Verlet integrator:<br>
    \f[ v(t+0.5 \Delta t) = v(t) + 0.5 \Delta t f(t) \f] */
void propagate_vel();
/** Propagate the positions. Integration step 2 of the Velocity Verletintegrator:<br>
    \f[ p(t+\Delta t) = p(t) + \Delta t  v(t+0.5 \Delta t) \f] */
void propagate_pos();
/** Propagate the velocities and positions. Integration step 1 and 2
    of the Velocity Verlet integrator: <br>
    \f[ v(t+0.5 \Delta t) = v(t) + 0.5 \Delta t f(t) \f] <br>
    \f[ p(t+\Delta t) = p(t) + \Delta t  v(t+0.5 \Delta t) \f] */
void propagate_vel_pos();
/** Rescale all particle forces with \f[ 0.5 \Delta t^2 \f] and propagate the velocities.
    Integration step 4 of the Velocity Verletintegrator:<br> 
    \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t) \f] */
void rescale_forces_propagate_vel();

/** Integrator stability check (see compile flag ADDITIONAL_CHECKS). */
void force_and_velocity_check(Particle *p); 
/** Integrator stability check (see compile flag ADDITIONAL_CHECKS). */
void force_and_velocity_display();
 
void finalize_p_inst_npt();

/*@}*/

void integrator_sanity_checks()
{
  //char *errtext;

  if ( time_step < 0.0 ) {
      ostringstream msg;
      msg <<"time_step not set";
      runtimeError(msg);
  }
}

#ifdef NPT

void integrator_npt_sanity_checks()
{  
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    if (nptiso.piston <= 0.0) {
        ostringstream msg;
        msg <<"npt on, but piston mass not set";
        runtimeError(msg);
    }

#ifdef ELECTROSTATICS

    switch(coulomb.method) {
      case COULOMB_NONE:  break;
      case COULOMB_DH:    break;
      case COULOMB_RF:    break;
#ifdef P3M
      case COULOMB_P3M:   break;
#endif /*P3M*/
      default: {
        ostringstream msg;
        msg <<"npt only works with P3M, Debye-Huckel or reaction field";
        runtimeError(msg);
      }
    }
#endif /*ELECTROSTATICS*/

#ifdef DIPOLES

    switch (coulomb.Dmethod) {
      case DIPOLAR_NONE: break;
#ifdef DP3M
      case DIPOLAR_P3M: break;
#endif /* DP3M */
      default: {
        ostringstream msg;
        msg <<"NpT does not work with your dipolar method, please use P3M.";
        runtimeError(msg);
      }
    }
#endif  /* ifdef DIPOLES */
  }
}
#endif /*NPT*/

/************************************************************/
void integrate_ensemble_init()
{
#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    /* prepare NpT-integration */
    nptiso.inv_piston = 1/(1.0*nptiso.piston);
    nptiso.p_inst_av = 0.0;
    if ( nptiso.dimension == 0 ) {
      fprintf(stderr,"%d: INTERNAL ERROR: npt integrator was called but dimension not yet set. this should not happen. ", this_node);
      errexit();
    }

    nptiso.volume = pow(box_l[nptiso.non_const_dim],nptiso.dimension);

    if (recalc_forces) { 
      nptiso.p_inst = 0.0;  
      nptiso.p_vir[0] = nptiso.p_vir[1] = nptiso.p_vir[2] = 0.0;
      nptiso.p_vel[0] = nptiso.p_vel[1] = nptiso.p_vel[2] = 0.0;
#ifdef CONFIGTEMP
      configtemp[0] = configtemp[1] = 0.0;
#endif
    }
  }
#endif
}

/************************************************************/

void integrate_vv(int n_steps, int reuse_forces)
{
  /* Prepare the Integrator */
  on_integration_start();
 
  #ifdef IMMERSED_BOUNDARY
    // Here we initialize volume conservation
    // This function checks if the reference volumes have been set and if necessary calculates them
    IBM_InitVolumeConservation();
  #endif

  /* if any method vetoes (P3M not initialized), immediately bail out */
  if (check_runtime_errors())
    return;

#ifdef MULTI_TIMESTEP
  if (smaller_time_step > 0.){
    mts_max = time_step/smaller_time_step;
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO) {
      current_time_step_is_small = 1;
      // Compute forces for small timestep -> get virial contribution.
      if (recalc_forces)
        thermo_heat_up();
      force_calc();
      thermo_cool_down();
      ghost_communicator(&cell_structure.collect_ghost_force_comm);
      current_time_step_is_small = 0;
      // Store virial
      for(int j=0;j<3;++j)
        virial_store[j] = nptiso.p_vir[j];
      thermo_heat_up();
      force_calc();
      thermo_cool_down();
      ghost_communicator(&cell_structure.collect_ghost_force_comm);
      rescale_forces();
    }
#endif
  }
#endif

  /* Verlet list criterion */
  skin2 = SQR(0.5 * skin);

  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv: integrating %d steps (recalc_forces=%d)\n",
                      this_node, n_steps, recalc_forces));
   
  /* Integration Step: Preparation for first integration step:
     Calculate forces f(t) as function of positions p(t) ( and velocities v(t) ) */
  /* reuse_forces logic:
     -1: recalculate forces unconditionally, mostly used for timing
      0: recalculate forces if recalc_forces is set, meaning it is probably necessary
      1: do not recalculate forces. Mostly when reading checkpoints with forces
   */
  if (reuse_forces == -1 || (recalc_forces && reuse_forces != 1)) {
    thermo_heat_up();

#ifdef LB
    transfer_momentum = 0;
    if (lattice_switch & LATTICE_LB && this_node == 0)
      if (warnings) fprintf (stderr, "Warning: Recalculating forces, so the LB coupling forces are not included in the particle force the first time step. This only matters if it happens frequently during sampling.\n");
#endif
#ifdef LB_GPU
    transfer_momentum_gpu = 0;
    if (lattice_switch & LATTICE_LB_GPU && this_node == 0)
      if (warnings) fprintf (stderr, "Warning: Recalculating forces, so the LB coupling forces are not included in the particle force the first time step. This only matters if it happens frequently during sampling.\n");
#endif

    force_calc();

    if(integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      rescale_forces();
#ifdef ROTATION
      convert_initial_torques();
#endif
    }

    thermo_cool_down();

#ifdef MULTI_TIMESTEP
#ifdef NPT
    if (smaller_time_step > 0. && integ_switch == INTEG_METHOD_NPT_ISO) 
      for(int j=0;j<3;++j)
        nptiso.p_vir[j] += virial_store[j];
#endif
#endif

#ifdef COLLISION_DETECTION
    handle_collisions();
#endif
  }

#ifdef GHMC
  if(thermo_switch & THERMO_GHMC)
    ghmc_init();
#endif

  if (thermo_switch & THERMO_CPU)
    mpi_thermalize_cpu(temperature);

  if (check_runtime_errors())
    return;

  n_verlet_updates = 0;

  /* Integration loop */
  for (int step=0; step<n_steps; step++) {
    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n", this_node, step));

#ifdef BOND_CONSTRAINT
    save_old_pos();
#endif

#ifdef GHMC
    if(thermo_switch & THERMO_GHMC) {
      if (step % ghmc_nmd == 0)
        ghmc_momentum_update();
    }
#endif

#ifdef SD
    if (thermo_switch & THERMO_SD) {
      fprintf(stderr,"Warning: Use integrate_sd to use Stokesian Dynamics Thermalizer.");
    }
    if (thermo_switch & THERMO_BD) {
      fprintf(stderr,"Warning: Use integrate_sd to use Brownian Dynamics Thermalizer.");
    }
#endif

    /* Integration Steps: Step 1 and 2 of Velocity Verlet scheme:
       v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
       p(t + dt)   = p(t) + dt * v(t+0.5*dt)
       NOTE 1: Prefactors do not occur in formulas since we use 
       rescaled forces and velocities. 
       NOTE 2: Depending on the integration method Step 1 and Step 2 
       cannot be combined for the translation. 
    */
    if (integ_switch == INTEG_METHOD_NPT_ISO || nemd_method != NEMD_METHOD_OFF) {
      propagate_vel();
      propagate_pos(); 
    } else if(integ_switch == INTEG_METHOD_STEEPEST_DESCENT) {
      if(steepest_descent_step())
	break;
    } else { 
      propagate_vel_pos();
    }

#ifdef BOND_CONSTRAINT
    /**Correct those particle positions that participate in a rigid/constrained bond */
    cells_update_ghosts();

    correct_pos_shake();
#endif

#ifdef ELECTROSTATICS
    if(coulomb.method == COULOMB_MAGGS) {
      maggs_propagate_B_field(0.5*time_step); 
    }
#endif

#ifdef NPT
    if (check_runtime_errors())
      break;
#endif

#ifdef MULTI_TIMESTEP
    if (smaller_time_step > 0){
      current_time_step_is_small = 1;
      /* Calculate the forces */
      thermo_heat_up();
      force_calc();
      thermo_cool_down();
      ghost_communicator(&cell_structure.collect_ghost_force_comm);
      rescale_forces();
      for (mts_index=0;mts_index<mts_max;++mts_index) {
        /* Small integration steps */
        /* Propagate velocities and positions */
        /* Assumes: not NEMD_METHOD_OFF; NPT not updated during small steps */
        if(integ_switch == INTEG_METHOD_NPT_ISO || nemd_method != NEMD_METHOD_OFF) {
          propagate_vel();  propagate_pos(); }
        else
          propagate_vel_pos();        
        cells_update_ghosts();
        force_calc();
        ghost_communicator(&cell_structure.collect_ghost_force_comm);
#ifdef NPT
        // Store virial
        for(int j=0;j<3;++j)
          virial_store[j] = nptiso.p_vir[j];
#endif
        rescale_forces_propagate_vel();
      }
      current_time_step_is_small = 0;             
      thermo_heat_up();
      force_calc();
      thermo_cool_down();
      ghost_communicator(&cell_structure.collect_ghost_force_comm);
      rescale_forces();
      recalc_forces=0;
    }
#endif


    /* Integration Step: Step 3 of Velocity Verlet scheme:
       Calculate f(t+dt) as function of positions p(t+dt) ( and velocities v(t+0.5*dt) ) */

#ifdef LB
    transfer_momentum = 1;
#endif
#ifdef LB_GPU
    transfer_momentum_gpu = 1;
#endif

    force_calc();
    
// IMMERSED_BOUNDARY
#ifdef IMMERSED_BOUNDARY
    // Now the forces are computed and need to go into the LB fluid
    if (lattice_switch & LATTICE_LB) IBM_ForcesIntoFluid_CPU();
#ifdef LB_GPU
    if (lattice_switch & LATTICE_LB_GPU) IBM_ForcesIntoFluid_GPU();
#endif
#endif

#ifdef CATALYTIC_REACTIONS
    integrate_reaction();
#endif

    if (check_runtime_errors())
      break;

#ifdef MULTI_TIMESTEP
#ifdef NPT
    if (smaller_time_step > 0. && integ_switch == INTEG_METHOD_NPT_ISO) 
      for(int j=0;j<3;++j)
        nptiso.p_vir[j] += virial_store[j];
#endif
#endif
    /* Integration Step: Step 4 of Velocity Verlet scheme:
       v(t+dt) = v(t+0.5*dt) + 0.5*dt * f(t+dt) */
    if(integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      rescale_forces_propagate_vel();
#ifdef ROTATION
    convert_torques_propagate_omega();
#endif
    }
    // SHAKE velocity updates
#ifdef BOND_CONSTRAINT
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    correct_vel_shake();
#endif
    //VIRTUAL_SITES update vel
#ifdef VIRTUAL_SITES
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    update_mol_vel();
    if (check_runtime_errors()) break;
#endif

    // progagate one-step functionalities
#ifdef LB
    if (lattice_switch & LATTICE_LB)
      lattice_boltzmann_update();
      
    if (check_runtime_errors())
      break;
#endif

#ifdef LB_GPU
    if(this_node == 0){
#ifdef ELECTROKINETICS
      if (ek_initialized) {
        ek_integrate();
      }
      else {
#endif
        if (lattice_switch & LATTICE_LB_GPU)
          lattice_boltzmann_update_gpu();
#ifdef ELECTROKINETICS
      }
#endif
    }
#endif //LB_GPU
    
// IMMERSED_BOUNDARY
#ifdef IMMERSED_BOUNDARY
    
    IBM_UpdateParticlePositions();
    // We reset all since otherwise the halo nodes may not be reset
    // NB: the normal Espresso reset is also done after applying the forces
//    if (lattice_switch & LATTICE_LB) IBM_ResetLBForces_CPU();
#ifdef LB_GPU
    //if (lattice_switch & LATTICE_LB_GPU) IBM_ResetLBForces_GPU();
#endif
    
    if (check_runtime_errors()) break;
    
    // Ghost positions are now out-of-date
    // We should update.
    // Actually we seem to get the same results whether we do this here or not, but it is safer to do it
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    
#endif // IMMERSED_BOUNDARY

#ifdef ELECTROSTATICS
    if(coulomb.method == COULOMB_MAGGS) {
      maggs_propagate_B_field(0.5*time_step); 
    }
#endif

#ifdef NPT
    if((this_node==0) && (integ_switch == INTEG_METHOD_NPT_ISO))
      nptiso.p_inst_av += nptiso.p_inst;
#endif

#ifdef GHMC
    if(thermo_switch & THERMO_GHMC) {
      if (step % ghmc_nmd == ghmc_nmd-1)
        ghmc_mc();
    }
#endif

    if(integ_switch != INTEG_METHOD_STEEPEST_DESCENT) {
      /* Propagate time: t = t+dt */
      sim_time += time_step;
    }
  }

  /* verlet list statistics */
  if(n_verlet_updates>0) verlet_reuse = n_steps/(double) n_verlet_updates;
  else verlet_reuse = 0;

#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    nptiso.invalidate_p_vel = 0;
    MPI_Bcast(&nptiso.p_inst, 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&nptiso.p_diff, 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&nptiso.volume, 1, MPI_DOUBLE, 0, comm_cart);
    if(this_node==0) nptiso.p_inst_av /= 1.0*n_steps;
    MPI_Bcast(&nptiso.p_inst_av, 1, MPI_DOUBLE, 0, comm_cart);
  }
#endif

#ifdef GHMC
  if(thermo_switch & THERMO_GHMC)
    ghmc_close();
#endif

}

/************************************************************/

void rescale_velocities(double scale) 
{
  Particle *p;
  int i, np, c;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      p[i].m.v[0] *= scale;
      p[i].m.v[1] *= scale;
      p[i].m.v[2] *= scale;
    }
  }
}


/* Privat functions */
/************************************************************/

void rescale_forces()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  double scale ;

  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces:\n",this_node));

  scale = 0.5 * time_step * time_step;
#ifdef MULTI_TIMESTEP
  if (smaller_time_step > 0.) {
    if (current_time_step_is_small)
      scale = 0.5 * smaller_time_step * smaller_time_step;
    else
      scale = 0.5 * smaller_time_step *         time_step;
  }
#endif
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      check_particle_force(&p[i]);
      p[i].f.f[0] *= scale/PMASS(p[i]);
      p[i].f.f[1] *= scale/PMASS(p[i]);
      p[i].f.f[2] *= scale/PMASS(p[i]);

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));

    }
  }
}

void rescale_forces_propagate_vel()
{
  Cell *cell;
  Particle *p;
  int i, j, np, c;
  double scale;

#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO){
    nptiso.p_vel[0] = nptiso.p_vel[1] = nptiso.p_vel[2] = 0.0;}
#endif

  scale = 0.5 * time_step * time_step;
#ifdef MULTI_TIMESTEP
  if (smaller_time_step > 0.) {
    if (current_time_step_is_small)
      scale = 0.5 * smaller_time_step * smaller_time_step;
    else
      scale = 0.5 * smaller_time_step *         time_step;
  }
#endif
  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces_propagate_vel:\n",this_node));

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      check_particle_force(&p[i]);
      /* Rescale forces: f_rescaled = 0.5*dt*dt * f_calculated * (1/mass) */
      p[i].f.f[0] *= scale/PMASS(p[i]);
      p[i].f.f[1] *= scale/PMASS(p[i]);
      p[i].f.f[2] *= scale/PMASS(p[i]);

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
#ifdef VIRTUAL_SITES
       // Virtual sites are not propagated during integration
       if (ifParticleIsVirtual(&p[i])) continue; 
#endif
      for(j = 0; j < 3 ; j++) {
#ifdef EXTERNAL_FORCES
	if (!(p[i].p.ext_flag & COORD_FIXED(j))) {
#endif
#ifdef NPT
          if(integ_switch == INTEG_METHOD_NPT_ISO && ( nptiso.geometry & nptiso.nptgeom_dir[j] )) {
            nptiso.p_vel[j] += SQR(p[i].m.v[j])*PMASS(p[i]);
#ifdef MULTI_TIMESTEP
            if (smaller_time_step > 0. && current_time_step_is_small==1)
              p[i].m.v[j] += p[i].f.f[j];
            else
#endif
              p[i].m.v[j] += p[i].f.f[j] + friction_therm0_nptiso(p[i].m.v[j])/PMASS(p[i]);
          }
          else
#endif
            /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * f(t+dt) */
            p[i].m.v[j] += p[i].f.f[j]; 
#ifdef EXTERNAL_FORCES
        }
#endif
      }

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
    }
  }
#ifdef NPT
#ifdef MULTI_TIMESTEP
  if (smaller_time_step < 0. || current_time_step_is_small == 0)
#endif
    finalize_p_inst_npt();
#endif
}

void finalize_p_inst_npt()
{
#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    double p_tmp=0.0;
    int i;
    /* finalize derivation of p_inst */
    nptiso.p_inst = 0.0;
    for ( i = 0 ; i < 3 ; i++ ) {
      if( nptiso.geometry & nptiso.nptgeom_dir[i] ) {
#ifdef MULTI_TIMESTEP
        if (smaller_time_step > 0.)
          nptiso.p_vel[i] /= SQR(smaller_time_step);
        else
#endif
          nptiso.p_vel[i] /= SQR(time_step);
        nptiso.p_inst += nptiso.p_vir[i] + nptiso.p_vel[i];
      }
    }

    MPI_Reduce(&nptiso.p_inst, &p_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
    if (this_node == 0) {
      nptiso.p_inst = p_tmp/(nptiso.dimension*nptiso.volume);
      nptiso.p_diff = nptiso.p_diff  +  (nptiso.p_inst-nptiso.p_ext)*0.5*time_step + friction_thermV_nptiso(nptiso.p_diff);
    }
  }
#endif
}

void propagate_press_box_pos_and_rescale_npt()
{
#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    Cell *cell;
    Particle *p;
    int i, j, np, c;
    double scal[3]={0.,0.,0.}, L_new=0.0;

    /* finalize derivation of p_inst */
#ifdef MULTI_TIMESTEP
    if (smaller_time_step < 0. || current_time_step_is_small == 0)
#endif
      finalize_p_inst_npt();

    /* adjust \ref nptiso_struct::nptiso.volume; prepare pos- and vel-rescaling */
    if (this_node == 0) {
#ifdef MULTI_TIMESTEP
      if (smaller_time_step < 0. || current_time_step_is_small == 0)
#endif      
        nptiso.volume += nptiso.inv_piston*nptiso.p_diff*0.5*time_step;
      scal[2] = SQR(box_l[nptiso.non_const_dim])/pow(nptiso.volume,2.0/nptiso.dimension);
#ifdef MULTI_TIMESTEP
      if (smaller_time_step < 0. || current_time_step_is_small == 0)
#endif
        nptiso.volume += nptiso.inv_piston*nptiso.p_diff*0.5*time_step;
      if (nptiso.volume < 0.0) {

          ostringstream msg;
          msg << "your choice of piston= "<< nptiso.piston << ", dt= " << time_step << ", p_diff= " << nptiso.p_diff
                 << " just caused the volume to become negative, decrease dt";
          runtimeError(msg);
	nptiso.volume = box_l[0]*box_l[1]*box_l[2];
	scal[2] = 1;
      }

      L_new = pow(nptiso.volume,1.0/nptiso.dimension);
      // printf("Lnew, %f: volume, %f: dim, %f: press, %f \n", L_new, nptiso.volume, nptiso.dimension,nptiso.p_inst );
      // fflush(stdout);

      scal[1] = L_new/box_l[nptiso.non_const_dim];
      scal[0] = 1/scal[1];
#ifdef MULTI_TIMESTEP
      if (smaller_time_step > 0.) {
        if (current_time_step_is_small == 1) {
          // load scal variable
          scal[0] = scal_store[0];
          scal[1] = scal_store[1];
          scal[2] = scal_store[2];
        } else {
          // save scal variable
          scal_store[0] = scal[0];
          scal_store[1] = scal[1];
          scal_store[2] = scal[2];
        }
      }
#endif
    }
    MPI_Bcast(scal,  3, MPI_DOUBLE, 0, comm_cart);
    

    /* propagate positions while rescaling positions and velocities */
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c]; p  = cell->part; np = cell->n;
      for(i = 0; i < np; i++) { 
#ifdef VIRTUAL_SITES
        if (ifParticleIsVirtual(&p[i])) continue;
#endif
        for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
          if (!(p[i].p.ext_flag & COORD_FIXED(j))) {
#endif
            if(nptiso.geometry & nptiso.nptgeom_dir[j]) {
#ifdef MULTI_TIMESTEP
              if (smaller_time_step > 0.) { 
                if (current_time_step_is_small == 1) {
                  if (mts_index==mts_max-1) {
                    p[i].r.p[j]      = scal[1]*(p[i].r.p[j] + scal[2]*p[i].m.v[j]);
                    p[i].l.p_old[j] *= scal[1];
                    p[i].m.v[j]     *= scal[0];
                  }
                  else
                    p[i].r.p[j]     += p[i].m.v[j];
                }
              } else
#endif
              {
                p[i].r.p[j]      = scal[1]*(p[i].r.p[j] + scal[2]*p[i].m.v[j]);
                p[i].l.p_old[j] *= scal[1];
                p[i].m.v[j]     *= scal[0];          
              }
            } else {
#ifdef MULTI_TIMESTEP
              if (smaller_time_step < 0. || current_time_step_is_small == 1) 
#endif              
                p[i].r.p[j] += p[i].m.v[j];
            }


#ifdef EXTERNAL_FORCES
          }
#endif
        }
        ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT:PV_1 v_new=(%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
        ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT:PPOS p=(%.3f,%.3f,%.3f)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2])); 
#ifdef ADDITIONAL_CHECKS
        force_and_velocity_check(&p[i]); 
#endif
      }
    }
    resort_particles = 1; 

    /* Apply new volume to the box-length, communicate it, and account for necessary adjustments to the cell geometry */
    if (this_node == 0) {
      for ( i = 0 ; i < 3 ; i++ ){ 
        if ( nptiso.geometry & nptiso.nptgeom_dir[i] ) {
          box_l[i] = L_new;
        } else if ( nptiso.cubic_box ) {
          box_l[i] = L_new;
        }
      }
    }
    MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, comm_cart);

    /* fast box length update */
    grid_changed_box_l();
    recalc_maximal_cutoff();
    cells_on_geometry_change(CELL_FLAG_FAST);
  }
#endif
}

void propagate_vel()
{
  Cell *cell;
  Particle *p;
  int c, i, j, np;
#ifdef NPT
  nptiso.p_vel[0] = nptiso.p_vel[1] = nptiso.p_vel[2] = 0.0;
#endif

  INTEG_TRACE(fprintf(stderr,"%d: propagate_vel:\n",this_node));

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
#ifdef ROTATION
     propagate_omega_quat_particle(&p[i]);
#endif

        // Don't propagate translational degrees of freedom of vs
#ifdef VIRTUAL_SITES
       if (ifParticleIsVirtual(&p[i])) continue;
#endif
      for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	if (!(p[i].p.ext_flag & COORD_FIXED(j)))	
#endif
          {
#ifdef NPT
            if(integ_switch == INTEG_METHOD_NPT_ISO && (nptiso.geometry & nptiso.nptgeom_dir[j] )) {
#ifdef MULTI_TIMESTEP
              if (smaller_time_step > 0. && current_time_step_is_small==1)
                p[i].m.v[j] += p[i].f.f[j];
              else
#endif
                p[i].m.v[j] += p[i].f.f[j] + friction_therm0_nptiso(p[i].m.v[j])/PMASS(p[i]);
              nptiso.p_vel[j] += SQR(p[i].m.v[j])*PMASS(p[i]);
            }
            else
#endif
              /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
              p[i].m.v[j] += p[i].f.f[j];

            /* SPECIAL TASKS in particle loop */
#ifdef NEMD
            if(j==0) nemd_get_velocity(p[i]);
#endif
          }

        ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
#ifdef ADDITIONAL_CHECKS
      force_and_velocity_check(&p[i]);
#endif
      }
    }
  }
#ifdef ADDITIONAL_CHECKS
  force_and_velocity_display();
#endif

  /* SPECIAL TASKS after velocity propagation */
#ifdef NEMD
  nemd_change_momentum();
  nemd_store_velocity_profile();
#endif
}

void propagate_pos()
{
  INTEG_TRACE(fprintf(stderr,"%d: propagate_pos:\n",this_node));
  if(integ_switch == INTEG_METHOD_NPT_ISO)
    /* Special propagator for NPT ISOTROPIC */
    /* Propagate pressure, box_length (2 times) and positions, rescale
       positions and velocities and check verlet list criterion (only NPT) */
    propagate_press_box_pos_and_rescale_npt();
  else {
    Cell *cell;
    Particle *p;
    int c, i, j, np;

    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
#ifdef VIRTUAL_SITES
       if (ifParticleIsVirtual(&p[i])) continue;
#endif
        for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	  if (!(p[i].p.ext_flag & COORD_FIXED(j)))
#endif
            {
#ifdef NEMD
              /* change momentum of each particle in top and bottom slab */
              if(j==0) nemd_add_velocity(&p[i]);
#endif
              /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
              p[i].r.p[j] += p[i].m.v[j];
            }
        }
        /* Verlet criterion check */
        if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) resort_particles = 1;
      }
    }
  }
  announce_resort_particles();
}

void propagate_vel_pos()
{
  Cell *cell;
  Particle *p;
  int c, i, j, np;

  INTEG_TRACE(fprintf(stderr,"%d: propagate_vel_pos:\n",this_node));

#ifdef ADDITIONAL_CHECKS
  db_max_force = db_max_vel = 0;
  db_maxf_id = db_maxv_id = -1;
#endif

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {

#ifdef ROTATION
      propagate_omega_quat_particle(&p[i]);
#endif

       // Don't propagate translational degrees of freedom of vs
#ifdef VIRTUAL_SITES
       if (ifParticleIsVirtual(&p[i])) continue;
#endif
     for(j=0; j < 3; j++){   
#ifdef EXTERNAL_FORCES
        if (!(p[i].p.ext_flag & COORD_FIXED(j)))
#endif
        {
          /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
          p[i].m.v[j] += p[i].f.f[j];

#ifdef MULTI_TIMESTEP
        if (smaller_time_step < 0. || current_time_step_is_small==1)
#endif  
            /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
            p[i].r.p[j] += p[i].m.v[j];
        }

      }

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PPOS p = (%.3e,%.3e,%.3e)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2]));

#ifdef ADDITIONAL_CHECKS
      force_and_velocity_check(&p[i]);
#endif

#ifdef LEES_EDWARDS
      /* test for crossing of a y-pbc: requires adjustment of velocity.*/
      {
                    int   b1, delta_box;
                    b1           = (int)floor( p[i].r.p[1]*box_l_i[1]);
                    if( b1 != 0 ){
                         delta_box    = b1 - (int)floor(( p[i].r.p[1] - p[i].m.v[1])*box_l_i[1] );
                         if( abs(delta_box) > 1 ){
                             fprintf(stderr, "Error! Particle moved more than one box length in 1 step\n");
                             errexit();
                         } 
                         p[i].m.v[0]     -= delta_box * lees_edwards_rate;   
                         p[i].r.p[0]     -= delta_box * lees_edwards_offset; 
                         p[i].r.p[1]     -= delta_box * box_l[1];
                         p[i].l.i[1]     += delta_box; 
                         while( p[i].r.p[1] >  box_l[1] ) {p[i].r.p[1] -= box_l[1]; p[i].l.i[1]++;}
                         while( p[i].r.p[1] <  0.0 )      {p[i].r.p[1] += box_l[1]; p[i].l.i[1]--;}
                         resort_particles = 1;
                    }
                    /* Branch prediction on most systems should mean there is minimal cost here */ 
                    while( p[i].r.p[0] >  box_l[0] ) {p[i].r.p[0] -= box_l[0]; p[i].l.i[0]++;}
                    while( p[i].r.p[0] <  0.0 )      {p[i].r.p[0] += box_l[0]; p[i].l.i[0]--;}
                    while( p[i].r.p[2] >  box_l[2] ) {p[i].r.p[2] -= box_l[2]; p[i].l.i[2]++;}
                    while( p[i].r.p[2] <  0.0 )      {p[i].r.p[2] += box_l[2]; p[i].l.i[2]--;}
      }
#endif

      /* Verlet criterion check*/
      if(SQR(p[i].r.p[0]-p[i].l.p_old[0]) 
        +SQR(p[i].r.p[1]-p[i].l.p_old[1])
        +SQR(p[i].r.p[2]-p[i].l.p_old[2]) > skin2) 
            resort_particles=1;


    }
  }

#ifdef LEES_EDWARDS /* would be nice to be more refined about this */
  resort_particles = 1;
#endif

  announce_resort_particles();

#ifdef ADDITIONAL_CHECKS
  force_and_velocity_display();
#endif
}

void force_and_velocity_check(Particle *p)
{
#ifdef ADDITIONAL_CHECKS
  int i;
  double db_force,db_vel;
  /* distance_check */
  for (i = 0; i < 3; i++)
    if(fabs(p->r.p[i] - p->l.p_old[i]) > local_box_l[i]) {
      fprintf(stderr, "%d: particle %d moved further than local box length by %lf %lf %lf\n",
              this_node, p->p.identity, p->r.p[0] - p->l.p_old[0], p->r.p[1] - p->l.p_old[1],
              p->r.p[2] - p->l.p_old[2]);
    }

  /* force check */
  db_force = SQR(p->f.f[0])+SQR(p->f.f[1])+SQR(p->f.f[2]);
  if(db_force > skin2)
    fprintf(stderr,"%d: Part %d has force %f (%f,%f,%f)\n",
            this_node,p->p.identity,sqrt(db_force),
            p->f.f[0],p->f.f[1],p->f.f[2]);
  if(db_force > db_max_force) { db_max_force=db_force; db_maxf_id=p->p.identity; }
  /* velocity check */
  db_vel   = SQR(p->m.v[0])+SQR(p->m.v[1])+SQR(p->m.v[2]);
  if(db_vel > skin2)
    fprintf(stderr,"%d: Part %d has velocity %f (%f,%f,%f)\n",
            this_node,p->p.identity,sqrt(db_vel),
            p->m.v[0],p->m.v[1],p->m.v[2]);
  if(db_vel > db_max_vel) { db_max_vel=db_vel; db_maxv_id=p->p.identity; }
#endif
}

void force_and_velocity_display()
{
#ifdef ADDITIONAL_CHECKS
  if(db_max_force > skin2)
    fprintf(stderr,"%d: max_force=%e, part=%d f=(%e,%e,%e)\n",this_node,
            sqrt(db_max_force),db_maxf_id,local_particles[db_maxf_id]->f.f[0],
            local_particles[db_maxf_id]->f.f[1],local_particles[db_maxf_id]->f.f[2]);
  if(db_max_vel > skin2)
    fprintf(stderr,"%d: max_vel=%e, part=%d v=(%e,%e,%e)\n",this_node,
            sqrt(db_max_vel),db_maxv_id,local_particles[db_maxv_id]->m.v[0],
            local_particles[db_maxv_id]->m.v[1],local_particles[db_maxv_id]->m.v[2]);
#endif
}

int python_integrate(int n_steps, bool recalc_forces, bool reuse_forces)
{


  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if ( recalc_forces ) {
  	if ( reuse_forces ) {
      std::ostringstream msg;
      msg <<"cannot reuse old forces and recalculate forces";
      runtimeError(msg);
  	}
  	reuse_forces = -1;
  }

  /* go on with integrate <n_steps> */
  if(n_steps < 0) {
    std::ostringstream msg;
    msg <<"illegal number of steps (must be >0)";
    runtimeError(msg);
    return ES_ERROR;
  }

  /* if skin wasn't set, do an educated guess now */
  if (!skin_set) {
    if (max_cut == 0.0) {
      std::ostringstream msg;
      msg <<"cannot automatically determine skin, please set it manually";
      runtimeError(msg);
      return ES_ERROR;
    }
    skin = 0.4*max_cut;
    mpi_bcast_parameter(FIELD_SKIN);
  }

  /* perform integration */
  if (!correlations_autoupdate && !observables_autoupdate) {
    if (mpi_integrate(n_steps, reuse_forces))
      return ES_ERROR;
  } else  {
    for (int i=0; i<n_steps; i++) {
      if (mpi_integrate(1, reuse_forces))
        return ES_ERROR;
      reuse_forces=1;
      autoupdate_observables();
      autoupdate_correlations();
    }
    if (n_steps == 0){
      if (mpi_integrate(0, reuse_forces))
        return ES_ERROR;
    }
  }
  return ES_OK;
}

void integrate_set_nvt()
{
  integ_switch = INTEG_METHOD_NVT;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
}

/** Parse integrate npt_isotropic command */
int integrate_set_npt_isotropic(double ext_pressure, double piston, int xdir, int ydir, int zdir, bool cubic_box)
{
  nptiso.cubic_box = 0;
  nptiso.p_ext = ext_pressure;
  nptiso.piston = piston;

  if ( nptiso.piston <= 0.0 ) {
    std::ostringstream msg;
    msg <<"You must set <piston> as well before you can use this integrator!\n";
    runtimeError(msg);
    return ES_ERROR;
  }
  if ( xdir || ydir || zdir ) {
  		/* set the geometry to include rescaling specified directions only*/
  		nptiso.geometry = 0; nptiso.dimension = 0; nptiso.non_const_dim = -1;
  		if ( xdir ) {
  			nptiso.geometry = ( nptiso.geometry | NPTGEOM_XDIR );
  			nptiso.dimension += 1;
  			nptiso.non_const_dim = 0;
  		}
  		if ( ydir ) {
  			nptiso.geometry = ( nptiso.geometry | NPTGEOM_YDIR );
  			nptiso.dimension += 1;
  			nptiso.non_const_dim = 1;
  		}
  		if ( zdir ) {
  			nptiso.geometry = ( nptiso.geometry | NPTGEOM_ZDIR );
  			nptiso.dimension += 1;
  			nptiso.non_const_dim = 2;
  		}
  }
  else {
  	/* set the geometry to include rescaling in all directions; the default*/
  	nptiso.geometry = 0;
  	nptiso.geometry = ( nptiso.geometry | NPTGEOM_XDIR );
  	nptiso.geometry = ( nptiso.geometry | NPTGEOM_YDIR );
  	nptiso.geometry = ( nptiso.geometry | NPTGEOM_ZDIR );
  	nptiso.dimension = 3; nptiso.non_const_dim = 2;
  }

  if ( cubic_box ) {
  	/* enable if the volume fluctuations should also apply to dimensions which are switched off by the above flags
       and which do not contribute to the pressure (3D) / tension (2D, 1D) */
  	nptiso.cubic_box = 1;
  }

  /* Sanity Checks */
#ifdef ELECTROSTATICS
  if ( nptiso.dimension < 3 && !nptiso.cubic_box && coulomb.bjerrum > 0 ){
  	fprintf(stderr,"WARNING: If electrostatics is being used you must use the -cubic_box option!\n");
  	fprintf(stderr,"Automatically reverting to a cubic box for npt integration.\n");
  	fprintf(stderr,"Be aware though that all of the coulombic pressure is added to the x-direction only!\n");
  	nptiso.cubic_box = 1;
  }
#endif

#ifdef DIPOLES
  if ( nptiso.dimension < 3 && !nptiso.cubic_box && coulomb.Dbjerrum > 0 ){
  	fprintf(stderr,"WARNING: If magnetostatics is being used you must use the -cubic_box option!\n");
  	fprintf(stderr,"Automatically reverting to a cubic box for npt integration.\n");
  	fprintf(stderr,"Be aware though that all of the magnetostatic pressure is added to the x-direction only!\n");
  	nptiso.cubic_box = 1;
  }
#endif


  if( nptiso.dimension == 0 || nptiso.non_const_dim == -1) {
    std::ostringstream msg;
    msg <<"You must enable at least one of the x y z components as fluctuating dimension(s) for box length motion!";
    msg <<"Cannot proceed with npt_isotropic, reverting to nvt integration... \n";
    runtimeError(msg);
    integ_switch = INTEG_METHOD_NVT;
  	mpi_bcast_parameter(FIELD_INTEG_SWITCH);
  	return (ES_ERROR);
  }

  /* set integrator switch */
  integ_switch = INTEG_METHOD_NPT_ISO;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);

  /* broadcast npt geometry information to all nodes */
  mpi_bcast_nptiso_geom();
  return (ES_OK);
}
