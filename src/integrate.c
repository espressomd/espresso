/*
  Copyright (C) 2010,2011 The ESPResSo project
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

/** \file integrate.c   Molecular dynamics integrator.
 *
 *  For more information about the integrator 
 *  see \ref integrate.h "integrate.h".
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "integrate.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "communication.h"
#include "grid.h"
#include "cells.h"
#include "verlet.h"
#include "rotation.h"
#include "ghosts.h"
#include "pressure.h"
#include "p3m.h"
#include "maggs.h"
#include "thermostat.h"
#include "initialize.h"
#include "forces.h"
#include "nsquare.h"
#include "domain_decomposition.h"
#include "layered.h"
#include "nemd.h"
#include "rattle.h"
#include "errorhandling.h"
#include "lattice.h"
#include "lb.h"
#include "virtual_sites.h"
#include "adresso.h"
#include "lbgpu.h"

/************************************************
 * DEFINES
 ************************************************/

/** Tag for communication in verlet fix: propagate_positions()  */
#define REQ_INT_VERLET   400

/*******************  variables  *******************/

int    integ_switch     = INTEG_METHOD_NVT;

int n_verlet_updates    = 0;

double time_step        = -1.0;
double sim_time         = 0.0;
double skin             = -1.0;
double skin2;
double max_range        = -1.0;
double max_range2       = -1.0;
double max_range_non_bonded  = 0.0;
double max_range_non_bonded2 = 0.0;

int    resort_particles = 1;
int    recalc_forces    = 1;

double verlet_reuse     = 0.0;

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

/************************************************************/

int tclcommand_invalidate_system(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  mpi_bcast_event(INVALIDATE_SYSTEM);
  return TCL_OK;
}

void local_invalidate_system()
{
  resort_particles = 1;
  invalidate_obs();
}

/**  Hand over integrate usage information to tcl interpreter. */
int tclcommand_integrate_print_usage(Tcl_Interp *interp) 
{
  Tcl_AppendResult(interp, "Usage of tcl-command integrate:\n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate <INT n steps>' for integrating n steps \n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate set' for printing integrator status \n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate set nvt' for enabling NVT integration or \n" , (char *)NULL);
#ifdef NPT
  Tcl_AppendResult(interp, "'integrate set npt_isotropic <DOUBLE p_ext> [<DOUBLE piston>] [<INT, INT, INT system_geometry>] [-cubic_box]' for enabling isotropic NPT integration \n" , (char *)NULL);
#endif
  return (TCL_ERROR);
}

/** Hand over integrate status information to tcl interpreter. */
int tclcommand_integrate_print_status(Tcl_Interp *interp) 
{
  int i;
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];
  switch (integ_switch) {
  case INTEG_METHOD_NVT:
    Tcl_AppendResult(interp, "{ set nvt }", (char *)NULL);
    return (TCL_OK);
  case INTEG_METHOD_NPT_ISO:
    Tcl_PrintDouble(interp, nptiso.p_ext, buffer);
    Tcl_AppendResult(interp, "{ set npt_isotropic ", buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nptiso.piston, buffer);
    Tcl_AppendResult(interp, " ",buffer, (char *)NULL);
    for ( i = 0 ; i < 3 ; i++){
      if ( nptiso.geometry & nptiso.nptgeom_dir[i] ) {
	sprintf(buffer, " %d", 1 );
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      } else {
	sprintf(buffer, " %d", 0 );
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
    }
    if ( nptiso.cubic_box ) {
      Tcl_AppendResult(interp, " -cubic_box", (char *)NULL);
    }
    Tcl_AppendResult(interp, " } ", (char *)NULL);

    return (TCL_OK);
  }
  return (TCL_ERROR);
}

/** Parse integrate nvt command */
int tclcommand_integrate_set_nvt(Tcl_Interp *interp, int argc, char **argv)
{
  integ_switch = INTEG_METHOD_NVT;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
  return (TCL_OK);
}

/** Parse integrate npt_isotropic command */
int tclcommand_integrate_set_npt_isotropic(Tcl_Interp *interp, int argc, char **argv)
{
  int xdir, ydir, zdir;
  xdir = ydir = zdir = nptiso.cubic_box = 0;

  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args: \n", (char *)NULL);
    return tclcommand_integrate_print_usage(interp);
  }  
  /* set parameters p_ext and piston */
  if ( !ARG_IS_D(3, nptiso.p_ext) )  return tclcommand_integrate_print_usage(interp);
  tclcallback_p_ext(interp, &nptiso.p_ext);
  if ( argc > 4 ) { 
    if(!ARG_IS_D(4, nptiso.piston) ) return tclcommand_integrate_print_usage(interp);
    tclcallback_npt_piston(interp, &nptiso.piston); }
  else if ( nptiso.piston <= 0.0 ) {
    Tcl_AppendResult(interp, "You must set <piston> as well before you can use this integrator! \n", (char *)NULL);
    return tclcommand_integrate_print_usage(interp);
  }

  if ( argc > 5 ) {
    if (!ARG_IS_I(5,xdir) || !ARG_IS_I(6,ydir) || !ARG_IS_I(7,zdir) ) {
      return tclcommand_integrate_print_usage(interp);}
    else {
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
  } else {
    /* set the geometry to include rescaling in all directions; the default*/
    nptiso.geometry = 0;
    nptiso.geometry = ( nptiso.geometry | NPTGEOM_XDIR );
    nptiso.geometry = ( nptiso.geometry | NPTGEOM_YDIR );
    nptiso.geometry = ( nptiso.geometry | NPTGEOM_ZDIR );
    nptiso.dimension = 3; nptiso.non_const_dim = 2;
  }

  if ( argc > 8 ) {
    /* enable if the volume fluctuations should also apply to dimensions which are switched off by the above flags
       and which do not contribute to the pressure (3D) / tension (2D, 1D) */
    if (!ARG_IS_S(8,"-cubic_box")) {
      return tclcommand_integrate_print_usage(interp);
    } else {
      nptiso.cubic_box = 1;
    }
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
    Tcl_AppendResult(interp, "You must enable at least one of the x y z components as fluctuating dimension(s) for box length motion!", (char *)NULL);
    Tcl_AppendResult(interp, "Cannot proceed with npt_isotropic, reverting to nvt integration... \n", (char *)NULL);
    integ_switch = INTEG_METHOD_NVT;
    mpi_bcast_parameter(FIELD_INTEG_SWITCH);
    return (TCL_ERROR);
  }

  /* set integrator switch */
  integ_switch = INTEG_METHOD_NPT_ISO;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);

  /* broadcast npt geometry information to all nodes */
  mpi_bcast_nptiso_geom();
  return (TCL_OK);
}

int tclcommand_integrate(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int  n_steps;
  
  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if (argc < 1) {
    Tcl_AppendResult(interp, "wrong # args: \n\"", (char *) NULL);
    return tclcommand_integrate_print_usage(interp);  }
  else if (argc < 2) {                    return tclcommand_integrate_print_status(interp); }

  if (ARG1_IS_S("set")) {
    if      (argc < 3)                    return tclcommand_integrate_print_status(interp);
    if      (ARG_IS_S(2,"nvt"))           return tclcommand_integrate_set_nvt(interp, argc, argv);
#ifdef NPT
    else if (ARG_IS_S(2,"npt_isotropic")) return tclcommand_integrate_set_npt_isotropic(interp, argc, argv);
#endif
    else {
      Tcl_AppendResult(interp, "unknown integrator method:\n", (char *)NULL);
      return tclcommand_integrate_print_usage(interp);
    }
  } else if ( !ARG_IS_I(1,n_steps) ) return tclcommand_integrate_print_usage(interp);

  /* go on with integrate <n_steps> */
  if(n_steps < 0) {
    Tcl_AppendResult(interp, "illegal number of steps (must be >0) \n", (char *) NULL);
    return tclcommand_integrate_print_usage(interp);;
  }
  /* perform integration */
  if (mpi_integrate(n_steps))
    return mpi_gather_runtime_errors(interp, TCL_OK);
  return TCL_OK;
}

/************************************************************/

void integrate_vv_recalc_maxrange()
{
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_recalc_maxrange:\n",this_node));

  /* maximal interaction cutoff */
  calc_maximal_cutoff();
  if (max_cut <= 0.0) {
    max_range  = -1.0;
    max_range2 = -1.0;
    return;
  }
  max_range            = max_cut;
  max_range_non_bonded = max_cut_non_bonded;
  /* at beginning be nice */
  if (skin > 0.0) {
    max_range            += skin;
    max_range_non_bonded += skin;
  }
  max_range2            = SQR(max_range);
  max_range_non_bonded2 = SQR(max_range_non_bonded);
}

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
    }
  }
#endif
}

/************************************************************/

void integrate_vv(int n_steps)
{
  int i;

  /* Prepare the Integrator */
  on_integration_start();

  /* if any method vetoes (P3M not initialized), immediately bail out */
  if (check_runtime_errors())
    return;

  /* Verlet list criterion */
  skin2 = SQR(0.5 * skin);

  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv: integrating %d steps (recalc_forces=%d)\n",
		      this_node, n_steps, recalc_forces));
   
  /* Integration Step: Preparation for first integration step:
     Calculate forces f(t) as function of positions p(t) ( and velocities v(t) ) */
  if (recalc_forces) {
    thermo_heat_up();
#ifdef LB
    transfer_momentum = 0;
#endif
#ifdef LB_GPU
    transfer_momentum_gpu = 0;
#endif
//VIRTUAL_SITES pos (and vel for DPD) update for security reason !!!
#ifdef VIRTUAL_SITES
    update_mol_vel_pos();
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    if (check_runtime_errors()) return;
#ifdef ADRESS
    //    adress_update_weights();
   if (check_runtime_errors()) return;
#endif
#endif

   
   force_calc();

   
   //VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
   ghost_communicator(&cell_structure.collect_ghost_force_comm);
   init_forces_ghosts();
   distribute_mol_force();
   if (check_runtime_errors()) return;
#endif

ghost_communicator(&cell_structure.collect_ghost_force_comm);

#ifdef ROTATION
    convert_initial_torques();
#endif

    thermo_cool_down();

    /* Communication Step: ghost forces */


    /*apply trap forces to trapped molecules*/
#ifdef MOLFORCES         
    calc_and_apply_mol_constraints();
#endif

    rescale_forces();
    recalc_forces = 0;

  }

  if (check_runtime_errors())
    return;

  n_verlet_updates = 0;

  /* Integration loop */
  for(i=0;i<n_steps;i++) {
    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,i));



#ifdef BOND_CONSTRAINT
    save_old_pos();
#endif

    /* Integration Steps: Step 1 and 2 of Velocity Verlet scheme:
       v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
       p(t + dt)   = p(t) + dt * v(t+0.5*dt)
       NOTE 1: Prefactors do not occur in formulas since we use 
               rescaled forces and velocities. 
       NOTE 2: Depending on the integration method Step 1 and Step 2 
               cannot be combined for the translation. 
    */
    if(integ_switch == INTEG_METHOD_NPT_ISO || nemd_method != NEMD_METHOD_OFF) {
      propagate_vel();  propagate_pos(); }
    else
      propagate_vel_pos();
#ifdef ROTATION
    propagate_omega_quat();
#endif

#ifdef BOND_CONSTRAINT
    /**Correct those particle positions that participate in a rigid/constrained bond */
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
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

    cells_update_ghosts();

//VIRTUAL_SITES update pos and vel (for DPD)
#ifdef VIRTUAL_SITES
   update_mol_vel_pos();
   ghost_communicator(&cell_structure.update_ghost_pos_comm);

   if (check_runtime_errors()) break;
#ifdef VIRTUAL_SITES_RELATIVE && LB 
   // This is on a workaround stage: 
   // When using virtual sites relative and LB at the same time, it is necessary 
   // to reassemble the cell lists after all position updates, also of virtual
   // particles. 
    cells_update_ghosts();
#endif

#ifdef ADRESS
   //adress_update_weights();
   if (check_runtime_errors()) break;
#endif
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

//VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
   ghost_communicator(&cell_structure.collect_ghost_force_comm);
   init_forces_ghosts();
   distribute_mol_force();
   if (check_runtime_errors()) break;
#endif

    /* Communication step: ghost forces */
    ghost_communicator(&cell_structure.collect_ghost_force_comm);

    /*apply trap forces to trapped molecules*/
#ifdef MOLFORCES         
    calc_and_apply_mol_constraints();
#endif

    if (check_runtime_errors())
      break;

    /* Integration Step: Step 4 of Velocity Verlet scheme:
       v(t+dt) = v(t+0.5*dt) + 0.5*dt * f(t+dt) */
    rescale_forces_propagate_vel();

#ifdef LB
  if (lattice_switch & LATTICE_LB) lattice_boltzmann_update();
  if (check_runtime_errors()) break;
#endif

#ifdef LB_GPU
  if(this_node == 0){
	if (lattice_switch & LATTICE_LB_GPU) lattice_boltzmann_update_gpu();
  }
#endif

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

#ifdef ELECTROSTATICS
    if(coulomb.method == COULOMB_MAGGS) {
      maggs_propagate_B_field(0.5*time_step); 
    }
#endif

#ifdef ROTATION
    convert_torques_propagate_omega();
#endif
#ifdef NPT
    if((this_node==0) && (integ_switch == INTEG_METHOD_NPT_ISO))
      nptiso.p_inst_av += nptiso.p_inst;
#endif

    /* Propagate time: t = t+dt */
    sim_time += time_step;
  }

  /* after simulating the forces are necessarily set. Necessary since
     resort_particles sets recalc_forces to 1 */
  recalc_forces = 0;

  /* verlet list statistics */
  if(n_verlet_updates>0) verlet_reuse = n_steps/(double) n_verlet_updates;
  else verlet_reuse = 0;

#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    nptiso.invalidate_p_vel = 0;
    MPI_Bcast(&nptiso.p_inst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nptiso.p_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nptiso.volume, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(this_node==0) nptiso.p_inst_av /= 1.0*n_steps;
    MPI_Bcast(&nptiso.p_inst_av, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
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

/* Callback functions */
/************************************************************/


int tclcallback_skin(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "skin must be positive.", (char *) NULL);
    return (TCL_ERROR);
  }
  skin = data;
  mpi_bcast_parameter(FIELD_SKIN);
  return (TCL_OK);
}

int tclcallback_time_step(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
#ifdef LB_GPU
  float ts = (float)data;
#endif
  if (data < 0.0) {
    Tcl_AppendResult(interp, "time step must be positive.", (char *) NULL);
    return (TCL_ERROR);
  }
#ifdef LB
  else if ((lbpar.tau >= 0.0) && (data > lbpar.tau)) {
    Tcl_AppendResult(interp, "MD time step must be smaller than LB time step.", (char *)NULL);
    return (TCL_ERROR);
  }
#endif
#ifdef LB_GPU	
  else if ((lbpar_gpu.tau >= 0.0) && (ts > lbpar_gpu.tau)) { 
    Tcl_AppendResult(interp, "MD time step must be smaller than LB time step.", (char *)NULL);
    return (TCL_ERROR);
  }
#endif
  mpi_set_time_step(data);

  return (TCL_OK);
}

int tclcallback_time(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  sim_time = data;
  mpi_bcast_parameter(FIELD_SIMTIME);
  return (TCL_OK);
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
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
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
  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces_propagate_vel:\n",this_node));

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
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
	if (!(p[i].l.ext_flag & COORD_FIXED(j))) {
#endif
#ifdef NPT
	  if(integ_switch == INTEG_METHOD_NPT_ISO && ( nptiso.geometry & nptiso.nptgeom_dir[j] )) {
	    nptiso.p_vel[j] += SQR(p[i].m.v[j])*PMASS(p[i]);
	    p[i].m.v[j] += p[i].f.f[j] + friction_therm0_nptiso(p[i].m.v[j])/PMASS(p[i]);
	  }
	  else
#endif
	    /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * f(t+dt) */
	    { p[i].m.v[j] += p[i].f.f[j]; }
#ifdef EXTERNAL_FORCES
	}
#endif
      }

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
    }
  }
#ifdef NPT
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
	nptiso.p_vel[i] /= SQR(time_step);
	nptiso.p_inst += nptiso.p_vir[i] + nptiso.p_vel[i];
      }
    }

    MPI_Reduce(&nptiso.p_inst, &p_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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

    rebuild_verletlist = 0;

    /* finalize derivation of p_inst */
    finalize_p_inst_npt();

    /* adjust \ref nptiso_struct::nptiso.volume; prepare pos- and vel-rescaling */
    if (this_node == 0) {
      nptiso.volume += nptiso.inv_piston*nptiso.p_diff*0.5*time_step;
      scal[2] = SQR(box_l[nptiso.non_const_dim])/pow(nptiso.volume,2.0/nptiso.dimension);
      nptiso.volume += nptiso.inv_piston*nptiso.p_diff*0.5*time_step;
      if (nptiso.volume < 0.0) {
	char *errtxt = runtime_error(128 + 3*TCL_DOUBLE_SPACE);
        ERROR_SPRINTF(errtxt, "{015 your choice of piston=%g, dt=%g, p_diff=%g just caused the volume to become negative, decrease dt} ",
                nptiso.piston,time_step,nptiso.p_diff);
	nptiso.volume = box_l[0]*box_l[1]*box_l[2];
	scal[2] = 1;
      }

      L_new = pow(nptiso.volume,1.0/nptiso.dimension);
      //      printf(stdout,"Lnew, %f: volume, %f: dim, %f: press, %f \n", L_new, nptiso.volume, nptiso.dimension,nptiso.p_inst );
      //    fflush(stdout);

      scal[1] = L_new/box_l[nptiso.non_const_dim];
      scal[0] = 1/scal[1];
    }
    MPI_Bcast(scal,  3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* propagate positions while rescaling positions and velocities */
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c]; p  = cell->part; np = cell->n;
      for(i = 0; i < np; i++) {	
#ifdef VIRTUAL_SITES
       if (ifParticleIsVirtual(&p[i])) continue;
#endif
       for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	if (!(p[i].l.ext_flag & COORD_FIXED(j))) {
#endif	    
	  if(nptiso.geometry & nptiso.nptgeom_dir[j]) {
	    p[i].r.p[j]      = scal[1]*(p[i].r.p[j] + scal[2]*p[i].m.v[j]);
	    p[i].l.p_old[j] *= scal[1];
	    p[i].m.v[j]     *= scal[0];
	  } else {
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
      /* Verlet criterion check */
      //if(distance2(p[i].r.p,p[i].l.p_old) > skin2 )
      rebuild_verletlist = 1; 
      }
    }

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
    MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    on_NpT_boxl_change(scal[1]);
    //on_NpT_boxl_change(0.0);

    /* communicate verlet criterion */
    //announce_rebuild_vlist();
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
#ifdef VIRTUAL_SITES
       if (ifParticleIsVirtual(&p[i])) continue;
#endif
      for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	if (!(p[i].l.ext_flag & COORD_FIXED(j)))	
#endif
	  {
#ifdef NPT
	    if(integ_switch == INTEG_METHOD_NPT_ISO && (nptiso.geometry & nptiso.nptgeom_dir[j] )) {
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

    rebuild_verletlist = 0;

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
	  if (!(p[i].l.ext_flag & COORD_FIXED(j)))
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
	if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) rebuild_verletlist = 1;
      }
    }
  }
  announce_rebuild_vlist();
}

void propagate_vel_pos()
{
  Cell *cell;
  Particle *p;
  int c, i, j, np;

  INTEG_TRACE(fprintf(stderr,"%d: propagate_vel_pos:\n",this_node));

  rebuild_verletlist = 0;

#ifdef ADDITIONAL_CHECKS
  db_max_force = db_max_vel = 0;
  db_maxf_id = db_maxv_id = -1;
#endif

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
	if (!(p[i].l.ext_flag & COORD_FIXED(j)))
#endif
	  {
	    /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
	    p[i].m.v[j] += p[i].f.f[j];

	    /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
	    p[i].r.p[j] += p[i].m.v[j];
	  }
      }

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PPOS p = (%.3f,%.3f,%.3f)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2]));

#ifdef ADDITIONAL_CHECKS
      force_and_velocity_check(&p[i]);
#endif

      /* Verlet criterion check */
      if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) rebuild_verletlist = 1;
    }
  }

  if(dd.use_vList) announce_rebuild_vlist();

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
      /* put any cellsystems here which do not rely on rebuild_verletlist */
      if ((cell_structure.type != CELL_STRUCTURE_NSQUARE) && (dd.use_vList)) {
	fprintf(stderr, "%d: particle %d moved further than local box length by %lf %lf %lf\n",
		this_node, p->p.identity, p->r.p[0] - p->l.p_old[0], p->r.p[1] - p->l.p_old[1],
		p->r.p[2] - p->l.p_old[2]);
      }
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
