// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file initialize.c
    Implementation of \ref initialize.h "initialize.h"
*/
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "config.h"
#include "debug.h"
#include "initialize.h"
#include "global.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "binary_file.h"
#include "integrate.h"
#include "statistics.h"
#include "energy.h"
#include "pressure.h"
#include "polymer.h"
#include "imd.h"
#include "random.h"
#include "communication.h"
#include "blockfile_tcl.h"
#include "cells.h"
#include "grid.h"
#include "thermostat.h"
#include "rotation.h"
#include "p3m.h"
#include "fft.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "maggs.h"
#include "elc.h"
#include "ghosts.h"
#include "debye_hueckel.h"
#include "forces.h"
#include "uwerr.h"
#include "utils.h"
#include "nsquare.h"
#include "nemd.h"
#include "domain_decomposition.h"
#include "errorhandling.h"
#include "rattle.h"

/** whether before integration the thermostat has to be reinitialized */
static int reinit_thermo = 1;
static int reinit_electrostatics = 0;

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
  cells_pre_init();
  ghost_init();
  /* Initialise force and energy tables */
  force_and_energy_tables_init();

#ifdef ELECTROSTATICS
  fft_pre_init();
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
  int i;

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

  for (i = 0; i < 3; i++)
    if (local_box_l[i] < max_range) {
      errtext = runtime_error(128 + TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtext,"{013 box_l in direction %d is still too small} ", i);
    }
  
#ifdef NPT
  if((integ_switch == INTEG_METHOD_NPT_ISO) && (nptiso.piston <= 0.0)) {
    char *errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{014 npt on, but piston mass not set} ");
  }
#endif

  if (!check_obs_calc_initialized()) return;

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

  if(resort_particles) {
    cells_resort_particles(CELL_GLOBAL_EXCHANGE);
    resort_particles = 0;
  }
  
  if(reinit_electrostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_electrostatics\n", this_node));
    switch (coulomb.method) {
    case COULOMB_P3M:
      P3M_count_charged_particles();
      break;
    case COULOMB_MAGGS: 
      Maggs_init(); 
      break;
    default: break;
    }
    reinit_electrostatics = 0;
  }
}

void on_particle_change()
{
  // EVENT_TRACE(fprintf(stderr, "%d: on_particle_change\n", this_node));
  resort_particles = 1;
  reinit_electrostatics = 1;
  rebuild_verletlist = 1;

  invalidate_obs();

  /* the particle information is no longer valid */
  freePartCfg();
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
  case COULOMB_P3M:
    P3M_init();
    integrate_vv_recalc_maxrange();
    on_parameter_change(FIELD_MAXRANGE);
    break;
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
  case COULOMB_MMM2D:
    MMM2D_init();
    break;
  case COULOMB_MAGGS: 
    Maggs_init(); 
    break;
  default: break;
  }
  if (coulomb.use_elc)
    ELC_init();

  recalc_forces = 1;
#endif
}

void on_short_range_ia_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_short_range_ia_changes\n", this_node));
  invalidate_obs();

  integrate_vv_recalc_maxrange();
  on_parameter_change(FIELD_MAXRANGE);

  recalc_forces = 1;
}

void on_constraint_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_constraint_change\n", this_node));
  invalidate_obs();

  recalc_forces = 1;  
}

void on_cell_structure_change()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_cell_structure_change\n", this_node));
  on_coulomb_change();
}

void on_resort_particles()
{
  EVENT_TRACE(fprintf(stderr, "%d: on_resort_particles\n", this_node));
#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_MMM2D:
    MMM2D_on_resort_particles();
    break;
  default: break;
  }
  if (coulomb.use_elc)
    ELC_on_resort_particles();
#endif
}

#ifdef NPT
void on_NpT_boxl_change(double scal1) {
  grid_changed_box_l();

#ifdef ELECTROSTATICS
  if(coulomb.method == COULOMB_P3M) {
    P3M_scaleby_box_l();
    integrate_vv_recalc_maxrange();
  }
#endif

  if(cell_structure.type==CELL_STRUCTURE_DOMDEC)
    dd_NpT_update_cell_grid(scal1);
}
#endif

void on_parameter_change(int field)
{
  /* to prevent two on_coulomb_change */
#ifdef ELECTROSTATICS
  int cc;
#endif

  EVENT_TRACE(fprintf(stderr, "%d: on_parameter_change %s\n", this_node, fields[field].name));

  if (field == FIELD_SKIN) {
    integrate_vv_recalc_maxrange();
    on_parameter_change(FIELD_MAXRANGE);
  }

  if (field == FIELD_NODEGRID)
    grid_changed_n_nodes();
  if (field == FIELD_BOXL || field == FIELD_NODEGRID)
    grid_changed_box_l();
    
  if (field == FIELD_TIMESTEP || field == FIELD_TEMPERATURE || field == FIELD_LANGEVIN_GAMMA
      || FIELD_DPD_GAMMA || field == FIELD_NPTISO_G0 || field == FIELD_NPTISO_GV || field == FIELD_NPTISO_PISTON )
    reinit_thermo = 1;

#ifdef NPT
  if ((field == FIELD_INTEG_SWITCH) && (integ_switch != INTEG_METHOD_NPT_ISO))
    nptiso.invalidate_p_vel = 1;  
#endif

#ifdef ELECTROSTATICS
  cc = 0;
  switch (coulomb.method) {
  case COULOMB_P3M:
    if (field == FIELD_TEMPERATURE || field == FIELD_NODEGRID)
      cc = 1;
    else if (field == FIELD_BOXL) {
      P3M_scaleby_box_l();
      integrate_vv_recalc_maxrange(); 
    }
    break;
  case COULOMB_DH:
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
  case COULOMB_MAGGS:
    /* Maggs electrostatics needs ghost velocities */
    on_ghost_flags_change();
    cells_re_init(CELL_STRUCTURE_CURRENT);    
    break;
  default: break;
  }
  if (coulomb.use_elc && (field == FIELD_TEMPERATURE || field == FIELD_BOXL))
    cc = 1;
  if (cc)
    on_coulomb_change();
#endif

  /* DPD needs ghost velocities, other thermostats not */
  if (field == FIELD_THERMO_SWITCH) {
    on_ghost_flags_change();
    cells_re_init(CELL_STRUCTURE_CURRENT);
  }

  if (field == FIELD_MAXRANGE)
    rebuild_verletlist = 1;

  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    if (field == FIELD_BOXL || field == FIELD_MAXRANGE || field == FIELD_THERMO_SWITCH)
      cells_re_init(CELL_STRUCTURE_LAYERED);
    break;
  case CELL_STRUCTURE_DOMDEC:
    if (field == FIELD_BOXL || field == FIELD_NODEGRID || field == FIELD_MAXRANGE ||
	field == FIELD_MINNUMCELLS || field == FIELD_MAXNUMCELLS || field == FIELD_THERMO_SWITCH)
      cells_re_init(CELL_STRUCTURE_DOMDEC);
    break;
  }
}

void on_ghost_flags_change()
{
  /* that's all we change here */
  extern int ghosts_have_v;
  ghosts_have_v = 0;
  
  if (thermo_switch & THERMO_DPD)
    /* DPD needs also ghost velocities */
    ghosts_have_v = 1;
#ifdef BOND_CONSTRAINT
  else if (n_rigidbonds)
    ghosts_have_v = 1;
#endif
#ifdef ELECTROSTATICS
  /* Maggs electrostatics needs ghost velocities too */
  else if(coulomb.method == COULOMB_MAGGS)
      ghosts_have_v = 1;
#endif
}

static void init_tcl(Tcl_Interp *interp)
{
  char cwd[1024];
  char *scriptdir;

  /*
    installation of tcl commands
  */

  /* in cells.c */
  Tcl_CreateCommand(interp, "cellsystem", (Tcl_CmdProc *)cellsystem, 0, NULL);
  /* in integrate.c */
  Tcl_CreateCommand(interp, "invalidate_system", (Tcl_CmdProc *)invalidate_system, 0, NULL);
  Tcl_CreateCommand(interp, "integrate", (Tcl_CmdProc *)integrate, 0, NULL);
  /* in global.c */
  Tcl_CreateCommand(interp, "setmd", (Tcl_CmdProc *)setmd, 0, NULL);
  /* in grid.c */
  Tcl_CreateCommand(interp, "change_volume", (Tcl_CmdProc *)change_volume, 0, NULL);
  /* in interaction_data.c */
  Tcl_CreateCommand(interp, "code_info", (Tcl_CmdProc *)code_info, 0, NULL);
  /* in interaction_data.c */
  Tcl_CreateCommand(interp, "inter", (Tcl_CmdProc *)inter, 0, NULL);
  /* in particle_data.c */
  Tcl_CreateCommand(interp, "part", (Tcl_CmdProc *)part, 0, NULL);
  /* in file binaryfile.c */
  Tcl_CreateCommand(interp, "writemd", (Tcl_CmdProc *)writemd, 0, NULL);
  Tcl_CreateCommand(interp, "readmd", (Tcl_CmdProc *)readmd, 0, NULL);
  /* in file statistics.c */
  Tcl_CreateCommand(interp, "analyze", (Tcl_CmdProc *)analyze, 0, NULL);
  /* in file polymer.c */
  Tcl_CreateCommand(interp, "polymer", (Tcl_CmdProc *)polymer, 0, NULL);
  Tcl_CreateCommand(interp, "counterions", (Tcl_CmdProc *)counterions, 0, NULL);
  Tcl_CreateCommand(interp, "salt", (Tcl_CmdProc *)salt, 0, NULL);
  Tcl_CreateCommand(interp, "velocities", (Tcl_CmdProc *)velocities, 0, NULL);
  Tcl_CreateCommand(interp, "maxwell_velocities", (Tcl_CmdProc *)maxwell_velocities, 0, NULL);
  Tcl_CreateCommand(interp, "crosslink", (Tcl_CmdProc *)crosslink, 0, NULL);
  Tcl_CreateCommand(interp, "diamond", (Tcl_CmdProc *)diamond, 0, NULL);
  Tcl_CreateCommand(interp, "icosaeder", (Tcl_CmdProc *)icosaeder, 0, NULL);
   /* in file imd.c */
  Tcl_CreateCommand(interp, "imd", (Tcl_CmdProc *)imd, 0, NULL);
  /* in file random.c */
  Tcl_CreateCommand(interp, "t_random", (Tcl_CmdProc *)t_random, 0, NULL);
  Tcl_CreateCommand(interp, "bit_random", (Tcl_CmdProc *)bit_random, 0, NULL);
  /* in file blockfile_tcl.c */
  Tcl_CreateCommand(interp, "blockfile", (Tcl_CmdProc *)blockfile, 0, NULL);
  /* in constraint.c */
  Tcl_CreateCommand(interp, "constraint", (Tcl_CmdProc *)constraint, 0, NULL);
  /* in uwerr.c */
  Tcl_CreateCommand(interp, "uwerr", (Tcl_CmdProc *)uwerr, 0, NULL);
  /* in nemd.c */
  Tcl_CreateCommand(interp, "nemd", (Tcl_CmdProc *)nemd, 0, NULL);
  /* in thermostat.c */
  Tcl_CreateCommand(interp, "thermostat", (Tcl_CmdProc *)thermostat, 0, NULL);

  /* evaluate the Tcl initialization script */
  scriptdir = getenv("ESPRESSO_SCRIPTS");
  fprintf(stderr,"%d: Script directory: %s\n",this_node,scriptdir);
  if (!scriptdir)
    scriptdir= "scripts";


  if ((getcwd(cwd, 1024) == NULL) || (chdir(scriptdir) != 0)) {
    fprintf(stderr,
	    "\n\ncould not change to script dir %s, please check ESPRESSO_SCRIPTS.\n\n\n",
	    scriptdir);
    exit(-1);
  }
  if (Tcl_EvalFile(interp, "init.tcl") == TCL_ERROR) {
    fprintf(stderr, "\n\nerror in initialization script: %s\n\n\n",
	    Tcl_GetStringResult(interp));
    exit(-1);
  }
  if (chdir(cwd) != 0) {
    fprintf(stderr,
	    "\n\ncould not change back to execution dir %s ????\n\n\n",
	    cwd);
    exit(-1);
  }
}
