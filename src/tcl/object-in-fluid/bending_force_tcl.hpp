/*
  Copyright (C) 2010,2013 The ESPResSo project
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
#ifndef BENDING_FORCE_TCL_H
#define BENDING_FORCE_TCL_H
/** \file bending_force.hpp Routines to calculate the bending_force energy or/and
 *  and force for a particle quadruple (two triangles that have 2 particles in common)
*/

#include "../../tcl/parser.hpp"
#include "interaction_data.hpp"


/// parse parameters for the bending_force potential
int tclcommand_inter_parse_bending_force(Tcl_Interp *interp, int bond_type, int argc, char **argv);

int tclprint_to_result_bendingforceIA(Tcl_Interp *interp, Bonded_ia_parameters *params);

#endif
