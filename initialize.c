/** \file initialize.c
    Implementation of \ref initialize.h "initialize.h"
*/
#include <unistd.h>
#include <stdlib.h>
#include "initialize.h"
#include "global.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "binary_file.h"
#include "integrate.h"
#include "statistics.h"
#include "imd.h"
#include "random.h"
#include "communication.h"
#include "blockfile_tcl.h"
#include "cells.h"
#include "grid.h"
#include "forces.h"
#include "thermostat.h"
#include "p3m.h"
#include "fft.h"
#include "ghosts.h"

static void init_tcl(Tcl_Interp *interp);

int on_program_start(Tcl_Interp *interp)
{
  /*
    call the initialization of the modules here
  */

  init_random();
  cells_pre_init();
  ghost_pre_init();
  fft_pre_init();


  /*
    call all initializations to don only on the master node here.
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

  fprintf(stderr,"%d: on_int_start: para_ch = %d inter_ch = %d top_ch %d part_ch =%d\n",this_node,
	  parameter_changed,interactions_changed,topology_changed,particle_changed);

  particle_invalidate_part_node();

  if (parameter_changed || interactions_changed || topology_changed) {
    integrate_vv_recalc_maxrange();
    invalidate_ghosts();
    exchange_and_sort_part();
    cells_re_init();
    ghost_init();
    force_init();
  }

  if (parameter_changed || topology_changed) {
    thermo_init();
  }

  if (interactions_changed || topology_changed) {
    P3M_init();
  }

  if (parameter_changed || particle_changed || interactions_changed || topology_changed) {
    rebuild_verletlist = 1;
  }

  /* the particle information is no longer valid */
  free(partCfg); partCfg=NULL;
}

void on_particle_change()
{
  particle_changed = 1;

  /* the particle information is no longer valid */
  free(partCfg); partCfg=NULL;
}

void on_topology_change()
{
  grid_changed_topology();
  cells_changed_topology();
  topology_changed = 1;
}

void on_ia_change()
{
  interactions_changed = 1;
}

void on_parameter_change()
{
  parameter_changed = 1;
}

static void init_tcl(Tcl_Interp *interp)
{
  char cwd[1024];
  char *scriptdir;

  /*
    installation of tcl commands
  */

  /* in integrate.c */
  Tcl_CreateCommand(interp, "integrate", (Tcl_CmdProc *)integrate, 0, NULL);
  /* in global.c */
  Tcl_CreateCommand(interp, "setmd", (Tcl_CmdProc *)setmd, 0, NULL);
  /* in interaction_data.c */
  Tcl_CreateCommand(interp, "inter", (Tcl_CmdProc *)inter, 0, NULL);
  /* in particle_data.c */
  Tcl_CreateCommand(interp, "part", (Tcl_CmdProc *)part, 0, NULL);
  /* in file binaryfile.c */
  Tcl_CreateCommand(interp, "writemd", (Tcl_CmdProc *)writemd, 0, NULL);
  Tcl_CreateCommand(interp, "readmd", (Tcl_CmdProc *)readmd, 0, NULL);
  /* in file statistics.c */
  Tcl_CreateCommand(interp, "mindist", (Tcl_CmdProc *)mindist, 0, NULL);
  Tcl_CreateCommand(interp, "analyze", (Tcl_CmdProc *)analyze, 0, NULL);
  /* in file imd.c */
  Tcl_CreateCommand(interp, "imd", (Tcl_CmdProc *)imd, 0, NULL);
  /* in file random.c */
  Tcl_CreateCommand(interp, "tcl_rand", (Tcl_CmdProc *)tcl_rand, 0, NULL);
  /* in file blockfile_tcl.c */
  Tcl_CreateCommand(interp, "blockfile", (Tcl_CmdProc *)blockfile, 0, NULL);

  /* evaluate the Tcl initialization script */
  scriptdir = getenv("TCLMD_SCRIPTS");
  if (!scriptdir)
    scriptdir= "scripts";

  if ((getcwd(cwd, 1024) == NULL) || (chdir(scriptdir) != 0)) {
    fprintf(stderr,
	    "\n\ncould not change to script dir %s, please check TCLMD_SCRIPTS.\n\n\n",
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
