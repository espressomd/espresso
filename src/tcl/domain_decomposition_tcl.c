/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file domain_decomposition.c
 *
 *  This file contains everything related to the cell system: domain decomposition.
 *  See also \ref domain_decomposition.h
 */
#include "utils.h"
#include "parser.h"

#include "domain_decomposition.h"

/** half the number of cell neighbors in 3 Dimensions. */
#define CELLS_MAX_NEIGHBORS 14

/*@}*/

/************************************************/
/** \name Variables */
/************************************************/
/*@{*/

extern int max_num_cells ;
extern int min_num_cells ;
extern double max_skin  ; 

/*@}*/


int tclcallback_max_num_cells(Tcl_Interp *interp, void *_data)
{
  int data = *(int *)_data;
  if (data < min_num_cells) {
    Tcl_AppendResult(interp, "max_num_cells cannot be smaller than min_num_cells", (char *) NULL);
    return (TCL_ERROR);
  }
  max_num_cells = data;
  mpi_bcast_parameter(FIELD_MAXNUMCELLS);
  return (TCL_OK);
}

int tclcallback_min_num_cells(Tcl_Interp *interp, void *_data)
{
  char buf[TCL_INTEGER_SPACE];
  int data = *(int *)_data;
  int min = calc_processor_min_num_cells();

  if (data < min) {
    sprintf(buf, "%d", min);
    Tcl_AppendResult(interp, "min_num_cells must be at least ", buf, (char *) NULL);
    return (TCL_ERROR);
  }
  if (data > max_num_cells) {
    Tcl_AppendResult(interp, "min_num_cells cannot be larger than max_num_cells", (char *) NULL);
    return (TCL_ERROR);
  }
  min_num_cells = data;
  mpi_bcast_parameter(FIELD_MINNUMCELLS);
  return (TCL_OK);
}


/************************************************************/
