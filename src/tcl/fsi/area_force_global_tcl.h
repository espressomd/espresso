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
#ifndef AREA_FORCE_GLOBAL_TCL_H
#define AREA_FORCE_GLOBAL_TCL_H
/** \file area_force_global.h
 *  Routines to calculate the AREA_FORCE_GLOBAL energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.c
*/


#include "tcl/parser.h"
#include "interaction_data.h"

int tclprint_to_result_areaforceglobalIA(Tcl_Interp *interp, Bonded_ia_parameters *params);

/// parse parameters for the area_force_global potential
int tclcommand_inter_parse_area_force_global(Tcl_Interp *interp, int bond_type, int argc, char **argv);

#endif 
