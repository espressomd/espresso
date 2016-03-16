/*
  Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
  
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
#ifndef _OBJECT_IN_FLUID_OIF_GLOBAL_FORCES_TCL_H
#define _OBJECT_IN_FLUID_OIF_GLOBAL_FORCES_TCL_H
/** \file oif_global_forces.hpp
 *  Routines to calculate the OIF_GLOBAL_FORCES energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/


#include "tcl/parser.hpp"
#include "interaction_data.hpp"

int tclprint_to_result_oifglobalforcesIA(Tcl_Interp *interp, Bonded_ia_parameters *params);

/// parse parameters for the oif_global_forces potential
int tclcommand_inter_parse_oif_global_forces(Tcl_Interp *interp, int bond_type, int argc, char **argv);

#endif 
