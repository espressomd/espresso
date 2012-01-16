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
/** \file lb-boundaries.h
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb-boundaries.c.
 *
 * In the current version only simple bounce back walls are implemented. Thus
 * after the streaming step, in all wall nodes all populations are bounced
 * back from where they came from. Ulf Schiller spent a lot of time
 * working on more powerful alternatives, they are to be found in the
 * lb_testing branch of espresso until the end of 2010. Now we stripped
 * down the code to a minimum, as most of it was not sufficiently understandable.
 *
 * Anyone who wants to revive these, please look into the git.
 *
 */

#ifndef LB_BOUNDARIES_TCL_H
#define LB_BOUNDARIES_TCL_H
#include <tcl.h>
#include "utils.h"
#include "halo.h"
#include "constraint.h"
#include "config.h"


// TCL Parser functions
extern int tclcommand_lbboundary(ClientData _data, Tcl_Interp *interp, int argc, char **argv);

#endif /* LB_BOUNDARIES_H */
