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
/** \file grid.c   Domain decomposition for parallel computing.
 *
 *  For more information on the domain decomposition, 
 *  see \ref grid.h "grid.h". 
*/
#include "utils.h"
#include "parser.h"
#include "communication.h"
#include "grid.h"

int tclcallback_node_grid(Tcl_Interp *interp, void *_data)
{
  int *data = (int *)_data;
  if ((data[0] < 0) || (data[1] < 0) || (data[2] < 0)) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }

  if (data[0]*data[1]*data[2] != n_nodes) {
    Tcl_AppendResult(interp, "node grid does not fit n_nodes",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  /* outsourced to 
     sort_int_array(data,3); */

  node_grid[0] = data[0];
  node_grid[1] = data[1];
  node_grid[2] = data[2];

  mpi_bcast_parameter(FIELD_NODEGRID);

  return (TCL_OK);
}

#ifdef PARTIAL_PERIODIC
int tclcallback_periodicity(Tcl_Interp *interp, void *_data)
{
  periodic = *(int *)_data;

  mpi_bcast_parameter(FIELD_PERIODIC);

  return (TCL_OK);
}
#else

int tclcallback_periodicity(Tcl_Interp *interp, void *_data)
{
  int tmp_periodic;
  tmp_periodic = *(int *)_data;
  if ((tmp_periodic & 7) == 7)
    return (TCL_OK);

  Tcl_AppendResult(interp, "periodic cannot be set since PARTIAL_PERIODIC not configured.", (char *)NULL);  
  return (TCL_ERROR);
}

#endif

int tclcallback_box_l(Tcl_Interp *interp, void *_data)
{
  double *data = _data;

  if ((data[0] <= 0) || (data[1] <= 0) || (data[2] <= 0)) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }

  box_l[0] = data[0];
  box_l[1] = data[1];
  box_l[2] = data[2];

  mpi_bcast_parameter(FIELD_BOXL);

  return (TCL_OK);
}

int tclcommand_change_volume(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  char buffer[50 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  char *mode;
  double d_new = box_l[0]; 
  int dir = -1;

  if (argc < 2) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: change_volume { <V_new> | <L_new> { x | y | z | xyz } }", (char *)NULL); return (TCL_ERROR);
  }
  if (Tcl_GetDouble(interp, argv[1], &d_new) == TCL_ERROR) return (TCL_ERROR);
  if (argc == 3) { 
    mode = argv[2];
    if (!strncmp(mode, "x", strlen(mode))) dir = 0;
    else if (!strncmp(mode, "y", strlen(mode))) dir = 1;
    else if (!strncmp(mode, "z", strlen(mode))) dir = 2;
    else if (!strncmp(mode, "xyz", strlen(mode))) dir = 3;
  }
  else if (argc > 3) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: change_volume { <V_new> | <L_new> { x | y | z | xyz } }", (char *)NULL); return (TCL_ERROR);
  }

  if (dir < 0) {
    d_new = pow(d_new,1./3.);
    rescale_boxl(3,d_new);
  }
  else { 
    rescale_boxl(dir,d_new); 
  }
  sprintf(buffer, "%f", box_l[0]*box_l[1]*box_l[2]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return gather_runtime_errors(interp, TCL_OK);
}

