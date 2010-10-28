/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
/** \file global.c
    Implementation of \ref global.h "global.h".
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "global.h"
/* from these modules we modify variables: */
#include "communication.h"
#include "cells.h"
#include "grid.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "integrate.h"
#include "thermostat.h"
#include "forces.h"
#include "verlet.h"
#include "p3m.h"
#include "imd.h"
#include "tuning.h"
#include "domain_decomposition.h"
#include "layered.h"
#include "pressure.h"
#include "rattle.h"
#include "lattice.h"
#include "adresso.h"

/**********************************************
 * description of variables
 * callbacks please define where the variables
 * comes from.
 **********************************************/

/** Read-only callback for \ref #fields.
    If you choose this, the variable cannot be
    changed by Tcl script code. */
int ro_callback(Tcl_Interp *interp, void *data);

/// List of all Tcl accessible global variables
const Datafield fields[] = {
  {box_l,            TYPE_DOUBLE, 3, "box_l",         boxl_callback,  1 },         /* 0  from grid.c */
  {dd.cell_grid,        TYPE_INT, 3, "cell_grid",     ro_callback,    6 },         /* 1  from cells.c */
  {dd.cell_size,     TYPE_DOUBLE, 3, "cell_size",     ro_callback,    6 },         /* 2  from cells.c */
  {&dpd_gamma,       TYPE_DOUBLE, 1, "dpd_gamma",     ro_callback,    5 },         /* 3  from thermostat.c */
  {&dpd_r_cut,       TYPE_DOUBLE, 1, "dpd_r_cut",     ro_callback,    5 },         /* 4  from thermostat.c */
  {&langevin_gamma,  TYPE_DOUBLE, 1, "gamma",         thermo_ro_callback, 1 },     /* 5  from thermostat.c */
  {&integ_switch,       TYPE_INT, 1, "integ_switch",  ro_callback,    1 },         /* 6  from integrate.c */
  {local_box_l,      TYPE_DOUBLE, 3, "local_box_l",   ro_callback,    2 },         /* 7  from global.c */
  {&max_cut,         TYPE_DOUBLE, 1, "max_cut",       ro_callback,    5 },         /* 8  from interaction_data.c */
  {&max_num_cells,      TYPE_INT, 1, "max_num_cells", max_num_cells_callback, 5 }, /* 9 from cells.c */
  {&max_seen_particle,  TYPE_INT, 1, "max_part",      ro_callback,    5 },         /* 10 from particle_data.c */
  {&max_range,       TYPE_DOUBLE, 1, "max_range",     ro_callback,    5 },         /* 11 from integrate.c */
  {&max_skin,        TYPE_DOUBLE, 1, "max_skin",      ro_callback,    5 },         /* 12 from integrate.c */
  {&min_num_cells,      TYPE_INT, 1, "min_num_cells", min_num_cells_callback, 5 }, /* 13  from cells.c */
  {&n_layers,           TYPE_INT, 1, "n_layers",      ro_callback,    3 },         /* 14 from layered.c */
  {&n_nodes,            TYPE_INT, 1, "n_nodes",       ro_callback,    3 },         /* 15 from communication.c */
  {&n_total_particles,  TYPE_INT, 1, "n_part",        ro_callback,    6 },         /* 16 from particle.c */
  {&n_particle_types,   TYPE_INT, 1, "n_part_types",  ro_callback,    8 },         /* 17 from interaction_data.c */
  {&n_rigidbonds,       TYPE_INT, 1, "n_rigidbonds",  ro_callback,    5 },         /* 18 from rattle.c */
  {node_grid,           TYPE_INT, 3, "node_grid",     node_grid_callback, 2 },     /* 19 from grid.c */
  {&nptiso_gamma0,   TYPE_DOUBLE, 1, "nptiso_gamma0", ro_callback,    13 },        /* 20 from thermostat.c */
  {&nptiso_gammav,   TYPE_DOUBLE, 1, "nptiso_gammav", ro_callback,    13 },        /* 21 from thermostat.c */
  {&nptiso.p_ext,    TYPE_DOUBLE, 1, "npt_p_ext",     ro_callback,     7 },        /* 22 from pressure.c */
  {&nptiso.p_inst,   TYPE_DOUBLE, 1, "npt_p_inst",    ro_callback,    10 },        /* 23 from pressure.c */
  {&nptiso.p_inst_av,TYPE_DOUBLE, 1, "npt_p_inst_av", ro_callback,    10 },        /* 24 from pressure.c */
  {&nptiso.p_diff,   TYPE_DOUBLE, 1, "npt_p_diff",    p_diff_callback, 7 },        /* 25 from pressure.c */
  {&nptiso.piston,   TYPE_DOUBLE, 1, "npt_piston",    piston_callback, 6 },        /* 26 from pressure.c */
  {&periodic,          TYPE_BOOL, 3, "periodicity",   per_callback,    1 },        /* 27 from grid.c */
  {&skin,            TYPE_DOUBLE, 1, "skin",          skin_callback,   2 },        /* 28 from integrate.c */
  {&temperature,     TYPE_DOUBLE, 1, "temperature",   thermo_ro_callback,   2 },   /* 29 from thermostat.c */
  {&thermo_switch,      TYPE_INT, 1, "thermo_switch", ro_callback,     2 },        /* 30 from thermostat.c */
  {&sim_time,        TYPE_DOUBLE, 1, "time",          time_callback,   4 },        /* 31 from integrate.c */
  {&time_step,       TYPE_DOUBLE, 1, "time_step",     time_step_callback, 5 },     /* 32 from integrate.c */
  {&timing_samples,     TYPE_INT, 1, "timings",       timings_callback, 4 },       /* 33 from tuning.c */
  {&transfer_rate,      TYPE_INT, 1, "transfer_rate", ro_callback,     2 },        /* 34 from imd.c */
  {&rebuild_verletlist,TYPE_BOOL, 1, "verlet_flag",   ro_callback,     8 },        /* 35 from verlet.c */
  {&verlet_reuse,    TYPE_DOUBLE, 1, "verlet_reuse",  ro_callback,     8 },        /* 36 from integrate.c */
  {&lattice_switch,     TYPE_INT, 1, "lattice_switch", ro_callback,    2 },          /* 37 from lattice.c */
  {&dpd_tgamma,      TYPE_DOUBLE, 1, "dpd_tgamma",    ro_callback,     6 },         /* 38 from thermostat.c */
  {&dpd_tr_cut,      TYPE_DOUBLE, 1, "dpd_tr_cut",    ro_callback,     6 },         /* 39 from thermostat.c */
  {&dpd_twf,            TYPE_INT, 1, "dpd_twf",    ro_callback,     6 },         /* 40 from thermostat.c */
  {&dpd_wf,             TYPE_INT, 1, "dpd_wf",    ro_callback,     5 },         /* 41 from thermostat.c */
  {adress_vars,      TYPE_DOUBLE, 7, "adress_vars",ro_callback,  1 },         /* 42  from adresso.c */
  { NULL, 0, 0, NULL, NULL, 0 }
};


/**********************************************
 * functions
 **********************************************/


int ro_callback(Tcl_Interp *interp, void *data)
{
  Tcl_AppendResult(interp, "variable is readonly", (char *)NULL);
  return (TCL_ERROR);
}

int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  char databuf[MAX_DIMENSION*(sizeof(int) + sizeof(double))];
  char buffer[TCL_DOUBLE_SPACE + 5];
  int i, j;
  int all = (argc == 1), writing = (argc >= 3);

  /* loop over all global variables. Has two purposes:
     either we write al variables or search for the one
     to write */
  for (i = 0; fields[i].data != NULL; i++) {
    if (all || !strncmp(argv[1], fields[i].name, strlen(argv[1]))) {
      if (!all) {
	if ((int)strlen(argv[1]) < fields[i].min_length) {
	  Tcl_AppendResult(interp, "Argument \"",argv[1],"\" not long ", (char *) NULL);
	  Tcl_AppendResult(interp, "enough to identify a setmd variable!", (char *) NULL);
	  return (TCL_ERROR);
	}
	if (writing) {
	  /* set */
	  /* parse in data */
	  if (argc != 2 + fields[i].dimension) {
	    sprintf(buffer, "%d", fields[i].dimension);	  
	    Tcl_AppendResult(interp, "\"", argv[1],
			     "\" has dimension ",
			     buffer, (char *) NULL);
	    sprintf(buffer, " not %d", argc - 2);	  
	    Tcl_AppendResult(interp, buffer, (char *) NULL);
	    return (TCL_ERROR);
	  }

	  /* get new value */
	  for (j = 0; j < fields[i].dimension; j++) {
	    switch (fields[i].type) {
	    case TYPE_INT:
	      if (Tcl_GetInt(interp, argv[2 + j], (int *)databuf + j) == TCL_ERROR)
		return (TCL_ERROR);
	      break;
	    case TYPE_BOOL: {
	      int dta;
	      if (Tcl_GetInt(interp, argv[2 + j], &dta))
		return (TCL_ERROR);
	      if (dta)
		*(int *)databuf |= (1L << j);
	      else
		*(int *)databuf &= ~(1L << j);
	      break;
	    }
	    case TYPE_DOUBLE:
	      if (Tcl_GetDouble(interp, argv[2 + j], (double *)databuf + j))
		return (TCL_ERROR);
	      break;
	    default: ;
	    }
	  }

	  if (fields[i].changeproc(interp, databuf) != TCL_OK)
	    return mpi_gather_runtime_errors(interp, TCL_ERROR);
	  /* fall through to write out the set value immediately again */
	}
      }

      /* get */
      if (all) {
	if (i != 0)
	  Tcl_AppendResult(interp, " ", (char *)NULL);
	Tcl_AppendResult(interp, "{", fields[i].name, " ", (char *)NULL);
      }
      for (j = 0; j < fields[i].dimension; j++) {
	switch (fields[i].type) {
	case TYPE_INT:
	  sprintf(buffer, "%d", ((int *)fields[i].data)[j]);
	  break;
	case TYPE_BOOL: {
	  if ((*(int *)fields[i].data) & (1L << j))
	    strcpy(buffer, "1");
	  else
	    strcpy(buffer, "0");
	  break;
	}
	case TYPE_DOUBLE:
	  Tcl_PrintDouble(interp, ((double *)fields[i].data)[j], buffer);
	  break;
	default: ;
	}
	Tcl_AppendResult(interp, buffer, (char *) NULL);
	if (j < fields[i].dimension - 1)
	  Tcl_AppendResult(interp, " ", (char *) NULL);
      }
      
      if (all)
	Tcl_AppendResult(interp, "}", (char *)NULL);
      /* wrote out one value, so skip rest */
      if (!all) {
	if (writing)
	  return mpi_gather_runtime_errors(interp, TCL_OK);
	else
	  return (TCL_OK);
      }
    }
  }
  if (all)
    return TCL_OK;

  Tcl_AppendResult(interp, "unknown md variable \"",
		   argv[1], "\"", (char *) NULL);
  return (TCL_ERROR);
}

int code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  if (argc < 2) {
    version_callback(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    compilation_callback(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    debug_callback(interp);
  }
  else {
    if(!strncmp(argv[1], "version" , strlen(argv[1]))) {
      version_callback(interp);
    }
    else if(!strncmp(argv[1], "compilation" , strlen(argv[1]))) {
      compilation_callback(interp);
    }
    else if(!strncmp(argv[1], "debug" , strlen(argv[1]))) {
      debug_callback(interp);
    }
    else {
      Tcl_AppendResult(interp, "info ",argv[1]," not known!", (char *) NULL);
    }
  }
  return (TCL_OK);
}
