/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file ghmc_tcl.c

    For more information see \ref ghmc_tcl.hpp
 */

#include <tcl.h>
#include <cmath>
#include "utils.hpp"
#include "parser.hpp"

#include "ghmc_tcl.hpp"
#include "ghmc.hpp"
#include "communication.hpp"
#include "thermostat.hpp"
#include "cells.hpp"
#include "statistics.hpp"

/************************************************************/
/* local prototypes                                         */
/************************************************************/

int tclcommand_ghmc_print_usage(Tcl_Interp *interp);
int tclcommand_ghmc_print_status(Tcl_Interp *interp);
int tclcommand_ghmc_print_statistics(Tcl_Interp *interp);

/************************************************************/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

int tclcommand_ghmc(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
#ifdef GHMC
  int status = TCL_OK;

  THERMO_TRACE(fprintf(stderr,"%d: ghmc:\n",this_node));
  Tcl_ResetResult(interp);

  /* print ghmc status */
  if(argc == 1) {
    status = tclcommand_ghmc_print_status(interp) ;
  }
  else if (ARG1_IS_S("statistics")) {
    status = tclcommand_ghmc_print_statistics(interp);
  }  
  else {
    Tcl_AppendResult(interp, "Unknown keyword: \n", (char *)NULL);
    status = tclcommand_ghmc_print_usage(interp);
  }

  return status;

#else

  INTEG_TRACE(fprintf(stderr,"%d: call to ghmc but not compiled in!\n",this_node));
  return tclcommand_ghmc_print_usage(interp);

#endif
}

int tclcommand_save_state(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
#ifdef GHMC

  Tcl_ResetResult(interp);

  /* save system state */
  save_last_state();

  return TCL_OK;

#else

  INTEG_TRACE(fprintf(stderr,"%d: call to ghmc but not compiled in!\n",this_node));
  return tclcommand_ghmc_print_usage(interp);

#endif
}

int tclcommand_load_state(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
#ifdef GHMC

  Tcl_ResetResult(interp);

  /* load last saved system state */
  load_last_state();
  cells_resort_particles(CELL_GLOBAL_EXCHANGE);
  invalidate_obs();

  return TCL_OK;

#else

  INTEG_TRACE(fprintf(stderr,"%d: call to ghmc but not compiled in!\n",this_node));
  return tclcommand_ghmc_print_usage(interp);

#endif
}

/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Hand over ghmc usage information to tcl interpreter. */
int tclcommand_ghmc_print_usage(Tcl_Interp *interp) 
{
#ifdef GHMC
  Tcl_AppendResult(interp, "Usage of tcl command ghmc:\n", (char *)NULL);
  Tcl_AppendResult(interp, "\"ghmc\" for returning the status or \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"ghmc statistics\" for returning the monte carlo statistics \n", (char *)NULL);
  return TCL_ERROR;
#else
  Tcl_AppendResult(interp, "ghmc not compiled in!", (char *)NULL);
  return TCL_ERROR;
#endif
}

#ifdef GHMC

/** Hand over ghmc status information to tcl interpreter. */
int tclcommand_ghmc_print_status(Tcl_Interp *interp) 
{
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];
  if (thermo_switch & THERMO_GHMC) {
		sprintf(buffer, "%d", ghmcdata.att);
    Tcl_AppendResult(interp, "mc moves attempted ",buffer, (char *)NULL);
		sprintf(buffer, "%d", ghmcdata.acc);
    Tcl_AppendResult(interp, ", mc moves accepted ",buffer, (char *)NULL);
    return TCL_OK;
  }
  else {
    Tcl_AppendResult(interp, "off", (char *)NULL);
    return TCL_OK;
  }
}

/** Hand over ghmc statistics information to tcl interpreter. */
int tclcommand_ghmc_print_statistics(Tcl_Interp *interp) 
{
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];
  if (thermo_switch & THERMO_GHMC) {
		sprintf(buffer, "%d", ghmc_att);
    Tcl_AppendResult(interp,buffer, (char *)NULL);
		sprintf(buffer, "%d", ghmc_acc);
    Tcl_AppendResult(interp, " ",buffer, (char *)NULL);
		sprintf(buffer, "%f", ghmc_acc == 0 ? 0.0 : (double) ghmc_acc/ghmc_att);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    return TCL_OK;
  }
  else {
    Tcl_AppendResult(interp, "{0 0 0}", (char *)NULL);
    return TCL_OK;
  }
}

/*@}*/

#endif
