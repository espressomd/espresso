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

#include "immersed-boundary/vvolume_tcl.hpp"
#include "immersed-boundary/vvolume.hpp"
#include "communication.hpp"
#include "parser.hpp"

int tclcallback_vescnum(Tcl_Interp *interp, void *_data) {
  int data = *(int *)_data;
  
  if (data < 0 || data > 200) {
    Tcl_AppendResult(interp, "vescnum must be positive and smaller than 200.", (char *) NULL);
    return (TCL_ERROR);
  }
  vescnum = data;
  mpi_bcast_parameter(FIELD_VESCNUM);
  return (TCL_OK);
}

int tclcallback_vvolo(Tcl_Interp *interp, void *_data) {
  double *data = (double*)_data;
  int i;
	
  for(i=0; i<200; i++) {
		
    if(data[i]<0.0) {
      Tcl_AppendResult(interp, "illegal value, Volume must be positive", (char *) NULL);
      return (TCL_ERROR);
    }
		
    VVolo[i]=data[i];
  }
	
  mpi_bcast_parameter(FIELD_VVOLO);
	
  return (TCL_OK);
}
