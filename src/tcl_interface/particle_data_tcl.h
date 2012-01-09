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
#ifndef _PARTICLE_DATA_TCL_H
#define _PARTICLE_DATA_TCL_H

#include <tcl.h>
#include "utils.h"
#include "grid.h"
#include "global.h"

/************************************************
 * defines
 ************************************************/

/**  bonds_flag "bonds_flag" value for updating particle config without bonding information */
#define WITHOUT_BONDS 0
/**  bonds_flag "bonds_flag" value for updating particle config with bonding information */
#define WITH_BONDS 1


#ifdef EXTERNAL_FORCES
/** \ref ParticleLocal::ext_flag "ext_flag" value for particle subject to an external force. */
#define PARTICLE_EXT_FORCE 1
/** \ref ParticleLocal::ext_flag "ext_flag" value for fixed coordinate coord. */
#define COORD_FIXED(coord) (2L << coord)
/** \ref ParticleLocal::ext_flag "ext_flag" mask to check wether any of the coordinates is fixed. */
#define COORDS_FIX_MASK     (COORD_FIXED(0) | COORD_FIXED(1) | COORD_FIXED(2))
#endif


/************************************************
 * Functions
 ************************************************/

/** Implementation of the tcl command \ref tclcommand_part. This command allows to
    modify particle data. */
int tclcommand_part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

#endif
