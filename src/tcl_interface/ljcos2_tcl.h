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
#ifndef _LJCOS2_TCL_H
#define _LJCOS2_TCL_H

/** \file ljcos2.h
 *  Routines to calculate the lennard-jones with cosine tail energy and/or  force 
 *  for a particle pair.  Cosine tail is different from that in ljcos.h
 *  Used for attractive tail/tail interactions in lipid bilayer calculations
 *  \ref forces.c
*/

#ifdef LJCOS2
#include <math.h>

/* These headers are needed to define types used in this header, hence
 * they are included here.  */
#include "particle_data.h"
#include "interaction_data.h"

int tclprint_to_result_ljcos2IA(Tcl_Interp *interp, int i, int j);


int tclcommand_inter_parse_ljcos2(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv);



#endif /* ifdef LJCOS2 */
#endif
