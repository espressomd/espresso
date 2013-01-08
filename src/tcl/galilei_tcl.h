/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
#ifndef _GALILEI_TCL_H
#define _GALILEI_TCL_H
#include "parser.h"

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** tcl procedure to set up reactions
*/

int tclcommand_kill_particle_motion(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);
int tclcommand_kill_particle_forces(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);
int tclcommand_system_CMS(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);
int tclcommand_system_CMS_velocity(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);
int tclcommand_galilei_transform(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);

/*@}*/

#endif
