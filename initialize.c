// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file initialize.c
    Implementation of \ref initialize.h "initialize.h"
*/
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
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
#include "ghosts.h"
#include "debye_hueckel.h"
#include "mmm1d.h"
#include "forces.h"
#include "uwerr.h"
#include "nsquare.h"

/** wether before integration the thermostat has to be reinitialized */
static int reinit_thermo = 1;
/** wether to recalculate the maximal interaction range before integration */
static int recalc_maxrange = 1;
/** wether the coulomb method needs to be reset before integration */
static int reinit_coulomb = 1;

static void init_tcl(Tcl_Interp *interp);

int on_program_start(Tcl_Interp *interp)
{
  /*
    call the initialization of the modules here
  */
  init_random();
  init_bit_random();

  setup_node_grid();
  cells_pre_init();
  ghost_init();

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
  
  /* happens if 0 particles... */
  if (!node_grid_is_set())
    setup_node_grid();

  particle_invalidate_part_node();

  if (recalc_maxrange) {
    integrate_vv_recalc_maxrange();
    on_parameter_change(FIELD_MAXRANGE);
    recalc_maxrange = 0;
    recalc_forces = 1;
  }

  if (reinit_thermo) {
    thermo_init();
    reinit_thermo = 0;
    recalc_forces = 1;
  }

#ifdef ELECTROSTATICS
  if (reinit_coulomb) {
    if(temperature > 0.0)
      coulomb.prefactor = coulomb.bjerrum * temperature; 
    else
      coulomb.prefactor = coulomb.bjerrum;
    switch (coulomb.method) {
    case COULOMB_P3M:
      P3M_init();
      break;
    case COULOMB_MMM1D:
      MMM1D_init();
      break;
    default: break;
    }
    reinit_coulomb = 0;
    recalc_forces = 1;
  }
#endif

  invalidate_obs();

  /* the particle information is no longer valid */
  free(partCfg); partCfg=NULL;
}

void on_particle_change()
{
  resort_particles = 1;
  
  invalidate_obs();

  /* the particle information is no longer valid */
  free(partCfg); partCfg=NULL;
}

void on_coulomb_change()
{
  invalidate_obs();

  reinit_coulomb = 1;
}

void on_short_range_ia_change()
{
  invalidate_obs();

  recalc_maxrange = 1;  
}

void on_constraint_change()
{
  invalidate_obs();

  recalc_forces = 1;  
}

void on_cell_structure_change()
{
  reinit_coulomb = 1;
}

void on_parameter_change(int field)
{
  if (field == FIELD_BOXL || field == FIELD_NODEGRID)
    grid_changed_topology();

  if (field == FIELD_TIMESTEP || field == FIELD_GAMMA || field == FIELD_TEMPERATURE)
    reinit_thermo = 1;

  if (field == FIELD_TEMPERATURE)
    reinit_coulomb = 1;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_P3M:
    if (field == FIELD_BOXL || field == FIELD_NODEGRID)
      reinit_coulomb = 1;
    break;
  case COULOMB_MMM1D:
    if (field == FIELD_BOXL)
      reinit_coulomb = 1;
    break;
  default: break;
  }
#endif

  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    if (field == FIELD_BOXL || field == FIELD_NODEGRID || field == FIELD_MAXRANGE ||
	field == FIELD_MAXNUMCELLS || field == FIELD_SKIN)
      cells_re_init(CELL_STRUCTURE_DOMDEC);
  }
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
  /* in interaction_data.c */
  Tcl_CreateCommand(interp, "constraint", (Tcl_CmdProc *)constraint, 0, NULL);
  /* in uwerr.c */
  Tcl_CreateCommand(interp, "uwerr", (Tcl_CmdProc *)uwerr, 0, NULL);

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
