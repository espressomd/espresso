/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
/** \file initialize.c
    Implementation of \ref initialize.h "initialize.h"
*/
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "utils.h"
#include "initialize.h"
#include "global.h"
#include "tcl/global_tcl.h"
#include "particle_data.h"
#include "tcl/particle_data_tcl.h"
#include "interaction_data.h"
#include "tcl/interaction_data_tcl.h"
#include "tcl/binary_file_tcl.h"
#include "integrate.h"
#include "statistics.h"
#include "energy.h"
#include "pressure.h"
#include "polymer.h"
#include "imd.h"
#include "tcl/imd_tcl.h"
#include "random.h"
#include "tcl/random_tcl.h"
#include "communication.h"
#include "cells.h"
#include "grid.h"
#include "tcl/grid_tcl.h"
#include "thermostat.h"
#include "tcl/thermostat_tcl.h"
#include "rotation.h"
#include "p3m.h"
#include "p3m-dipolar.h"
#include "ewald.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "maggs.h"
#include "elc.h"
#include "lb.h"
#include "ghosts.h"
#include "debye_hueckel.h"
#include "reaction_field.h"
#include "forces.h"
#include "tcl/uwerr_tcl.h"
#include "nsquare.h"
#include "nemd.h"
#include "tcl/nemd_tcl.h"
#include "domain_decomposition.h"
#include "errorhandling.h"
#include "rattle.h"
#include "bin.h"
#include "tcl/bin_tcl.h"
#include "lattice.h"
#include "iccp3m.h" /* -iccp3m- */
#include "tcl/iccp3m_tcl.h" 
#include "adresso.h"
#include "tcl/adresso_tcl.h"
#include "metadynamics.h"
#include "tcl/integrate_tcl.h"
#include "tcl/cells_tcl.h"
#include "tcl/statistics_tcl.h"
#include "tcl/blockfile_tcl.h"
#include "tcl/iccp3m_tcl.h"
#include "tcl/polymer_tcl.h"
#include "tcl/lb-boundaries_tcl.h"
#include "lb-boundaries.h"
#include "tcl/initialize_interpreter.h"

#ifdef CUDA
#include "tcl/cuda_init_tcl.h"
#endif

// import function from scriptsdir.c
char *get_default_scriptsdir();

/** whether the thermostat has to be reinitialized before integration */
static int reinit_thermo = 1;
static int reinit_electrostatics = 0;
static int reinit_magnetostatics = 0;
#ifdef LB_GPU
static int lb_reinit_particles_gpu = 1;
#endif

static void init_tcl(Tcl_Interp *interp);

int on_program_start(Tcl_Interp *interp)
{
  EVENT_TRACE(fprintf(stderr, "%d: on_program_start\n", this_node));

  /*
    call the initialization of the modules here
  */
  init_random();
  init_bit_random();

  setup_node_grid();
  /* calculate initial minimimal number of cells (see tclcallback_min_num_cells) */
  min_num_cells = calc_processor_min_num_cells();

  cells_pre_init();
  ghost_init();
  /* Initialise force and energy tables */
  force_and_energy_tables_init();
#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
  adress_force_and_energy_tables_init();
#endif
  /* #ifdef THERMODYNAMIC_FORCE */
  tf_tables_init();
  /* #endif */
#endif
#ifdef P3M
  p3m_pre_init();
#endif
#ifdef DP3M
  dp3m_pre_init();
#endif

#ifdef LB_GPU
  if(this_node == 0){
    //lb_pre_init_gpu();
  }
#endif
#ifdef LB
  lb_pre_init();
#endif

  /*
    call all initializations to do only on the master node here.
  */
  if (this_node == 0) {
    /* interaction_data.c: make sure 0<->0 ia always exists */
    make_particle_type_exist(0);
    
    init_tcl(interp);
  }
  return TCL_OK;

}


void on_integration_start()
{
  char *errtext;

  EVENT_TRACE(fprintf(stderr, "%d: on_integration_start\n", this_node));
  INTEG_TRACE(fprintf(stderr,"%d: on_integration_start: reinit_thermo = %d, resort_particles=%d\n",
		      this_node,reinit_thermo,resort_particles));

  /********************************************/
  /* sanity checks                            */
  /********************************************/

  if ( time_step < 0.0 ) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext, "{010 time_step not set} ");
  }
  if ( skin < 0.0 ) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{011 skin not set} ");
  }
  if ( temperature < 0.0 ) {
    errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{012 thermostat not initialized} ");
  }
  
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    if (nptiso.piston <= 0.0) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{014 npt on, but piston mass not set} ");
    }
    if(cell_structure.type!=CELL_STRUCTURE_DOMDEC) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{014 npt requires domain decomposition cellsystem} ");
    }
  
#ifdef ELECTROSTATICS

    switch(coulomb.method) {
      case COULOMB_NONE:  break;
#ifdef P3M
      case COULOMB_P3M:   break;
#endif /*P3M*/
      case COULOMB_EWALD: break;
      default: {
        char *errtext = runtime_error(128);
        ERROR_SPRINTF(errtext,"{014 npt only works with Ewald sum or P3M} ");
      }
    }
#endif /*ELECTROSTATICS*/
  }
  
#endif /*NPT*/

  if (!check_obs_calc_initialized()) return;

#ifdef LB
  if(lattice_switch & LATTICE_LB) {
    if (lbpar.agrid <= 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{098 Lattice Boltzmann agrid not set} ");
    }
    if (lbpar.tau <= 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{099 Lattice Boltzmann time step not set} ");
    }
    if (lbpar.rho <= 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{100 Lattice Boltzmann fluid density not set} ");
    }
    if (lbpar.viscosity <= 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{101 Lattice Boltzmann fluid viscosity not set} ");
    }
    if (dd.use_vList && skin>=lbpar.agrid/2.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext, "{104 LB requires either no Verlet lists or that the skin of the verlet list to be less than half of lattice-Boltzmann grid spacing.} ");
    }
  }
#endif
#ifdef LB_GPU
if(this_node == 0){
  if(lattice_switch & LATTICE_LB_GPU) {
    if (lbpar_gpu.agrid < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{098 Lattice Boltzmann agrid not set} ");
    }
    if (lbpar_gpu.tau < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{099 Lattice Boltzmann time step not set} ");
    }
    if (lbpar_gpu.rho < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{100 Lattice Boltzmann fluid density not set} ");
    }
    if (lbpar_gpu.viscosity < 0.0) {
      errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{101 Lattice Boltzmann fluid viscosity not set} ");
    }
    if (lb_reinit_particles_gpu) {
	lb_realloc_particles_gpu();
	lb_reinit_particles_gpu = 0;
    }
  }
}

#endif

#ifdef METADYNAMICS
  meta_init();
#endif

  /********************************************/
  /* end sanity checks                        */
  /********************************************/


  /* Prepare the thermostat */
  if (reinit_thermo) {
    thermo_init();
    reinit_thermo = 0;
    recalc_forces = 1;
  }

  /* Ensemble preparation: NVT or NPT */
  integrate_ensemble_init();

  /* Update particle and observable information for routines in statistics.c */
  invalidate_obs();
  freePartCfg();

  on_observable_calc();
}

void on_observable_calc()
{
  /* Prepare particle structure: Communication step: number of ghosts and ghost information */

  if(resort_particles)
    cells_resort_particles(CELL_GLOBAL_EXCHANGE);

#ifdef ELECTROSTATICS  
  if(reinit_electrostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_electrostatics\n", this_node));
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_ELC_P3M:
    case COULOMB_P3M:
      p3m_count_charged_particles();
      break;
#endif
    case COULOMB_EWALD:
      EWALD_count_charged_particles();
      break;
    case COULOMB_MAGGS: 
      maggs_init(); 
      break;
    default: break;
    }
    reinit_electrostatics = 0;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  if(reinit_magnetostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_magnetostatics\n", this_node));
    switch (coulomb.Dmethod) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:
    case DIPOLAR_P3M:
      dp3m_count_magnetic_particles();
      break;
#endif
    default: break;
    }
    reinit_magnetostatics = 0;
  }
#endif /*ifdef ELECTROSTATICS */
}

void on_particle_change()
{
  // EVENT_TRACE(fprintf(stderr, "%d: on_particle_change\n", this_node));
  resort_particles = 1;
  reinit_electrostatics = 1;
  reinit_magnetostatics = 1;

#ifdef LB_GPU
  lb_reinit_particles_gpu = 1;
#endif

  invalidate_obs();

  /* the particle information is no longer valid */
  freePartCfg();
}

void on_max_cut_change()
{
  cells_on_max_cut_change(1);
}

void on_coulomb_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_coulomb_change\n", this_node));
  invalidate_obs();

#ifdef ELECTROSTATICS
  if(temperature > 0.0)
    coulomb.prefactor = coulomb.bjerrum * temperature; 
  else
    coulomb.prefactor = coulomb.bjerrum;

  switch (coulomb.method) {
  case COULOMB_DH:
    on_max_cut_change();
    break;    
#ifdef P3M
  case COULOMB_ELC_P3M:
    ELC_init();
    // fall through
  case COULOMB_P3M:
    p3m_init();
    on_max_cut_change();
    break;
#endif
  case COULOMB_EWALD:
    EWALD_init();
    on_max_cut_change();
    break;
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
  case COULOMB_MMM2D:
    MMM2D_init();
    break;
  case COULOMB_MAGGS: 
    maggs_init();
    on_max_cut_change();
    /* Maggs electrostatics needs ghost velocities */
    on_ghost_flags_change();
    break;

    break;
  default: break;
  }
#endif  /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  if(temperature > 0.0)
    coulomb.Dprefactor = coulomb.Dbjerrum * temperature; 
  else
    coulomb.Dprefactor = coulomb.Dbjerrum;
  
  switch (coulomb.Dmethod) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:
       // fall through
  case DIPOLAR_P3M:
    dp3m_init();
    on_max_cut_change();
    break;
#endif
  default: break;
  }
#endif  /* ifdef DIPOLES */

  recalc_forces = 1;
}

void on_short_range_ia_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_short_range_ia_changes\n", this_node));
  invalidate_obs();
  on_max_cut_change();
  recalc_forces = 1;
}

void on_constraint_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_constraint_change\n", this_node));
  invalidate_obs();
  recalc_forces = 1;
}

void on_lbboundary_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_lbboundary_change\n", this_node));
  invalidate_obs();

#ifdef LB_BOUNDARIES
  if(lattice_switch & LATTICE_LB) {
    lb_init_boundaries();
  }
#endif
#ifdef LB_BOUNDARIES_GPU
  if(this_node == 0){
    if(lattice_switch & LATTICE_LB_GPU) {
      lb_init_boundaries();
    }
  }
#endif

  recalc_forces = 1;
}

void on_cell_structure_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_cell_structure_change\n", this_node));
  on_coulomb_change();
  /* 
#ifdef LB
  if (!lb_sanity_checks()) return;
#endif
  */
}

void on_resort_particles()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_resort_particles\n", this_node));
#ifdef ELECTROSTATICS
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    ELC_on_resort_particles();
    break;
#endif
  case COULOMB_MMM2D:
    MMM2D_on_resort_particles();
    break;
  case COULOMB_EWALD:
    EWALD_on_resort_particles();
    break;
  default: break;
  }
#endif /* ifdef ELECTROSTATICS */
  
  /* DIPOLAR interactions so far don't need this */

  recalc_forces = 1;
}

#ifdef NPT
void on_NpT_boxl_change(double scal1) {
  grid_changed_box_l();
  
#ifdef ELECTROSTATICS
  switch(coulomb.method) {
  case COULOMB_NONE:
    break;
#ifdef P3M
  case COULOMB_P3M:
    p3m_scaleby_box_l();
    cells_on_max_cut_change(0);
    break;
#endif
  case COULOMB_EWALD:
    EWALD_scaleby_box_l();
    cells_on_max_cut_change(0);
    break;
  default: {
    char *errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"NpT does not work with your electrostatics method, please use P3M.");
  }
  }
#endif

#ifdef DIPOLES
  switch(coulomb.Dmethod) {
  case DIPOLAR_NONE:
    break;
#ifdef DP3M
  case DIPOLAR_P3M:
    dp3m_scaleby_box_l();
    cells_on_max_cut_change(0);
    break;
#endif
  default: {
    char *errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"NpT does not work with your dipolar method, please use P3M.");
  }
  }
#endif

  if (cell_structure.type == CELL_STRUCTURE_DOMDEC)
    dd_change_boxl();
}
#endif

void on_parameter_change(int field)
{
  /* to prevent two on_coulomb_change calls. If this value is set to 1,
     we really need to reinitialize the whole electrostatics. */
#if defined(ELECTROSTATICS) || defined(DIPOLES)
  int cc = 0;
#endif

  EVENT_TRACE(fprintf(stderr, "%d: on_parameter_change %s\n", this_node, fields[field].name));

  if (field == FIELD_SKIN) {
    cells_on_max_cut_change(1);
  }

  if (field == FIELD_NODEGRID)
    grid_changed_n_nodes();
  if (field == FIELD_BOXL || field == FIELD_NODEGRID)
    grid_changed_box_l();
  if (field == FIELD_TIMESTEP || field == FIELD_TEMPERATURE || field == FIELD_LANGEVIN_GAMMA || field == FIELD_DPD_TGAMMA
      || field == FIELD_DPD_GAMMA || field == FIELD_NPTISO_G0 || field == FIELD_NPTISO_GV || field == FIELD_NPTISO_PISTON )
    reinit_thermo = 1;

#ifdef NPT
  if ((field == FIELD_INTEG_SWITCH) && (integ_switch != INTEG_METHOD_NPT_ISO))
    nptiso.invalidate_p_vel = 1;  
#endif

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    if (field == FIELD_TEMPERATURE || field == FIELD_BOXL)
      cc = 1;
    // fall through
  case COULOMB_P3M:
    if (field == FIELD_TEMPERATURE || field == FIELD_NODEGRID || field == FIELD_SKIN)
      cc = 1;
    else if (field == FIELD_BOXL) {
      p3m_scaleby_box_l();
      on_max_cut_change(); 
    }
    break;
#endif
  case COULOMB_EWALD:
    if (field == FIELD_TEMPERATURE || field == FIELD_SKIN)
      cc = 1;
    else if (field == FIELD_BOXL) {
      EWALD_scaleby_box_l();
      on_max_cut_change(); 
    }
    break;
  case COULOMB_DH:
    if (field == FIELD_TEMPERATURE)
      cc = 1;
    break;
  case COULOMB_RF:
  case COULOMB_INTER_RF:
    if (field == FIELD_TEMPERATURE)
      cc = 1;
    break;
  case COULOMB_MMM1D:
    if (field == FIELD_TEMPERATURE || field == FIELD_BOXL)
      cc = 1;
    break;
  case COULOMB_MMM2D:
    if (field == FIELD_TEMPERATURE || field == FIELD_BOXL || field == FIELD_NLAYERS)
      cc = 1;
    break;
  default: break;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  switch (coulomb.Dmethod) {
   #ifdef DP3M
    case DIPOLAR_MDLC_P3M:
     if (field == FIELD_TEMPERATURE || field == FIELD_BOXL)
       cc = 1;
      // fall through
    case DIPOLAR_P3M:
      if (field == FIELD_TEMPERATURE || field == FIELD_NODEGRID || field == FIELD_SKIN)
        cc = 1;
      else if (field == FIELD_BOXL) {
        dp3m_scaleby_box_l();
        on_max_cut_change(); 
      }
      break;
#endif
  default: break;
  }
#endif /*ifdef DIPOLES */
  
#if defined(ELECTROSTATICS) || defined(DIPOLES)
  if (cc)
    on_coulomb_change();
#endif

  /* DPD needs ghost velocities, other thermostats not */
  if (field == FIELD_THERMO_SWITCH) {
    on_ghost_flags_change();
  }

  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    if (field == FIELD_NODEGRID) {
      if (node_grid[0] != 1 || node_grid[1] != 1) {
	char *errtext = runtime_error(128);
	ERROR_SPRINTF(errtext, "{091 layered cellsystem requires 1 1 n node grid} ");
      }
    }
    if (field == FIELD_BOXL || field == FIELD_THERMO_SWITCH)
      cells_re_init(CELL_STRUCTURE_LAYERED);
    break;
  case CELL_STRUCTURE_DOMDEC:
    if (field == FIELD_BOXL || field == FIELD_NODEGRID ||
	field == FIELD_MINNUMCELLS || field == FIELD_MAXNUMCELLS ||
	field == FIELD_THERMO_SWITCH)
      cells_re_init(CELL_STRUCTURE_DOMDEC);
    break;
  }

#ifdef LB
  /* LB needs ghost velocities */
  if (field == FIELD_LATTICE_SWITCH) {
    on_ghost_flags_change();
  }

  if (lattice_switch & LATTICE_LB) {
    if (field == FIELD_TEMPERATURE) {
      lb_reinit_parameters();
    }

    if (field == FIELD_BOXL || field == FIELD_CELLGRID || field == FIELD_NNODES || field == FIELD_NODEGRID) {
      lb_init();
    }
  }
#endif

#ifdef LB_GPU
if(this_node == 0){
  if (lattice_switch & LATTICE_LB_GPU) {
    if (field == FIELD_TEMPERATURE) lb_reinit_parameters_gpu();
  }
}
#endif
}

#ifdef LB
void on_lb_params_change(int field) {

  if (field == LBPAR_AGRID) {
    lb_init();
  }
  if (field == LBPAR_DENSITY) {
    lb_reinit_fluid();
  }

  lb_reinit_parameters();

}
#endif

#if defined (LB) || defined (LB_GPU)
void on_lb_params_change_gpu(int field) {
#ifdef LB_GPU
  if (field == LBPAR_AGRID) {
    lb_init_gpu();
#ifdef LB_BOUNDARIES_GPU
    lb_init_boundaries();
#endif
  }
  if (field == LBPAR_DENSITY) {
    lb_reinit_fluid_gpu();
  }

  lb_reinit_parameters_gpu();
#endif
}
#endif

void on_ghost_flags_change()
{
  /* that's all we change here */
  extern int ghosts_have_v;

  int old_have_v = ghosts_have_v;

  ghosts_have_v = 0;
  
  /* DPD and LB need also ghost velocities */
  if (thermo_switch & THERMO_DPD)
    ghosts_have_v = 1;
#ifdef LB
  if (lattice_switch & LATTICE_LB)
    ghosts_have_v = 1;
#endif
#ifdef BOND_CONSTRAINT
  else if (n_rigidbonds)
    ghosts_have_v = 1;
#endif
#ifdef ELECTROSTATICS
  /* Maggs electrostatics needs ghost velocities too */
  else if(coulomb.method == COULOMB_MAGGS)
      ghosts_have_v = 1;
#endif
#ifdef INTER_DPD
  //maybe we have to add a new global to differ between compile in and acctual use.
  ghosts_have_v = 1;
#endif
#ifdef VIRTUAL_SITES 
  //VIRUTAL_SITES need v to update v of virtual sites
  ghosts_have_v = 1;
#endif

  if (old_have_v != ghosts_have_v)
    cells_re_init(CELL_STRUCTURE_CURRENT);    
}

static void init_tcl(Tcl_Interp *interp)
{
  char cwd[1024];
  char *scriptdir;

  /*
    installation of tcl commands
  */
  register_tcl_commands(interp);


  /* evaluate the Tcl initialization script */
  scriptdir = getenv("ESPRESSO_SCRIPTS");
  if (!scriptdir)
    scriptdir = get_default_scriptsdir();
  
  /*  fprintf(stderr,"Script directory: %s\n", scriptdir);*/

  if ((getcwd(cwd, 1024) == NULL) || (chdir(scriptdir) != 0)) {
    fprintf(stderr,
	    "\n\ncould not change to script dir %s, please check ESPRESSO_SCRIPTS.\n\n\n",
	    scriptdir);
    exit(1);
  }
  if (Tcl_EvalFile(interp, "init.tcl") == TCL_ERROR) {
    fprintf(stderr, "\n\nerror in initialization script: %s\n\n\n",
	    Tcl_GetStringResult(interp));
    exit(1);
  }
  if (chdir(cwd) != 0) {
    fprintf(stderr,
	    "\n\ncould not change back to execution dir %s ????\n\n\n",
	    cwd);
    exit(1);
  }
}
