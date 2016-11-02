/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

#include "minimize_energy_tcl.hpp"
#include "minimize_energy.hpp"
#include "communication.hpp"

static int usage(Tcl_Interp *interp) {
  Tcl_AppendResult(interp, "minimize_energy <f_max> <n_steps> <gamma> <max_displacement>\"\n", (char *)NULL);
  return TCL_ERROR;
}

int tclcommand_minimize_energy(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int  max_steps;
  double f_max, gamma, max_displacement;
  
  if (argc != 5) {
    Tcl_AppendResult(interp, "wrong # args: \n\"", (char *) NULL);
    return usage(interp);
  }
  else {
    if(!ARG_IS_D(1,f_max)) {
      return usage(interp);      
    }
    if(!ARG_IS_I(2,max_steps)) {
      return usage(interp);      
    }
    if(!ARG_IS_D(3,gamma)) {
      return usage(interp);      
    }
    if(!ARG_IS_D(4,max_displacement)) {
      return usage(interp);      
    }
  }

  minimize_energy_init(f_max, gamma, max_steps, max_displacement);
  mpi_minimize_energy();

  return TCL_OK;
}

