/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file mmm2d.h MMM2D algorithm for long range coulomb interaction
    in 2d+h geometries.  Implementation of the MMM2D method for the
    calculation of the electrostatic interaction for two dimensionally
    periodic systems. For details on the method see MMM general. The
    MMM2D method works only with the layered or nsquared \ref
    tclcommand_cellsystem "cell system". The tuning is not automated,
    since the only tunable parameter is the cell size, which can be
    changed easily in Tcl. Moreover, only a few values make sense to
    be tested, since in general the number of cells will be between 5
    and 20.
 */
#ifndef MMM2D_TCL_H
#define MMM2D_TCL_H
#include "parser.hpp"

#ifdef ELECTROSTATICS

/// print the mmm2d parameters to the interpreters result
int tclprint_to_result_MMM2D(Tcl_Interp *interp);

/// parse the mmm2d parameters
int tclcommand_inter_coulomb_parse_mmm2d(Tcl_Interp * interp, int argc, char ** argv);

#endif

#endif
