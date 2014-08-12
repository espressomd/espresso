/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#ifndef CG_DNA_TCL_H
#define CG_DNA_TCL_H
/** \file fene_tcl.hpp
 * Tcl interface for \ref fene.hpp
 */

#include "parser.hpp"
#include "interaction_data.hpp"
#include "cg_dna.hpp"

#ifdef CG_DNA

int tclcommand_inter_parse_cg_dna_basepair(Tcl_Interp *interp, int bond_type, int argc, char **argv);
int tclcommand_inter_parse_cg_dna_stacking(Tcl_Interp *interp, int bond_type, int argc, char **argv);

#endif

#endif
