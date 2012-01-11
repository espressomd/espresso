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
#ifndef _LJ_TCL_H
#define _LJ_TCL_H

/** \file lj.h
 *  Routines to calculate the lennard jones energy and/or  force 
 *  for a particle pair.
 *  \ref forces.c
*/

#ifdef LENNARD_JONES

/* These headers are needed to define types used in this header, hence
 * they are included here.  */
#include "particle_data.h"
#include "interaction_data.h"
#include "mol_cut.h"

int tclprint_to_result_ljIA(Tcl_Interp *interp, int i, int j);

/** set the force cap for the LJ interaction.
    @param ljforcecap the maximal force, 0 to disable, -1 for individual cutoff
    for each of the interactions.
*/


/// parser for the forcecap
int tclcommand_inter_parse_ljforcecap(Tcl_Interp * interp, int argc, char ** argv);

int tclcommand_inter_parse_lj(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv);

#endif /* ifdef LENNARD_JONES */
#endif

