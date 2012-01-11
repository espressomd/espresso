/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef NEMD_TCL_H
#define NEMD_TCL_H
/** \file nemd.h

    This file contains the implementation of the NEMD (Non Equilibrium
    Molecular Dynamics) algorithm. It allows one to shear a system
    with help of an unphysical momentum change in two slabs in the
    system.
 */

#include "nemd.h"
#include "particle_data.h"

#define NEMD_METHOD_OFF       0
#define NEMD_METHOD_EXCHANGE  1
#define NEMD_METHOD_SHEARRATE 2


/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** tcl procedure for nemd steering.
    USAGE: nemd \<n_slabs\> \<n_exchange\>   
    see also \ref tclcommand_nemd
*/
int tclcommand_nemd(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

#endif
