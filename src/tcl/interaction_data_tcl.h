/*
  Copyright (C) 2010,2011 The ESPResSo project
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
#ifndef IA_DATA_TCL_H
#define IA_DATA_TCL_H
/** \file interaction_data.h
    Various procedures concerning interactions between particles.

    interaction_data.h contains the parser \ref tclcommand_inter for the
    Tcl command "inter". Therefore the parsing of bonded and nonbonded
    interaction definition commands both is done here. It also contains
    procedures for low-level interactions setup and some helper functions.
    Moreover it contains code for the treatments of constraints.
*/

#include <tcl.h>


/************************************************
 * exported functions
 ************************************************/

/** Implementation of the tcl command "inter". This function
    allows the interaction parameters to be modified.
 */
int tclcommand_inter(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

/** Implementation of the Tcl function constraint. This function
    allows to set and delete constraints.
 */
int tclcommand_constraint(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv);

#endif
