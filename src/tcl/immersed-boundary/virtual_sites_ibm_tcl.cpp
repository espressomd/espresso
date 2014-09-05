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

#include "immersed-boundary/virtual_sites_ibm_tcl.hpp"
#include "global.hpp"
#ifdef VIRTUAL_SITES_IMMERSED_BOUNDARY
#include "immersed-boundary/virtual_sites_ibm.hpp"
#include "communication.hpp"

int tclcallback_integration_rule_ibm(Tcl_Interp *interp, void *_data) {
  int data = *(int *)_data;
  
  if (data < 0 || data > 1) {
    Tcl_AppendResult(interp, "integration rule for immersed boundary method is either 0 (euler) or 1 (2 step adams-bashforth).", (char *) NULL);
    return (TCL_ERROR);
  }
  integration_rule_ibm = data;
  mpi_bcast_parameter(FIELD_INTEGRATION_RULE_IBM);
  return (TCL_OK);
}

#endif
