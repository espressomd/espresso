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
#ifndef GHMC_TCL_H
#define GHMC_TCL_H
/** \file ghmc_tcl.hpp

    This file contains the tcl interface for the GHMC (Generalized 
    Hybrid MOnte Carlo) thermostat.
    
 */

#include <tcl.h>

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** tcl procedure to get ghmc statistics.
    USAGE: ghmc \<statistics\>
    see also \ref tclcommand_ghmc
*/
int tclcommand_ghmc(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

/** tcl procedure to save the system state.
    USAGE: save_state
    see also \ref tclcommand_save_state
*/

int tclcommand_save_state(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

/** tcl procedure to load last saved system state.
    USAGE: load_state
    see also \ref tclcommand_load_state
*/

int tclcommand_load_state(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

/*@}*/

/* endif GHMC_TCL_H */
#endif
