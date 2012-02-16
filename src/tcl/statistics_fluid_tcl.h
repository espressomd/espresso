/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
/** \file statistics_fluid.h
 *
 * Fluid related analysis functions.
 * Header file for \ref statistics_fluid.c.
 *
 */

#ifndef STATISTICS_FLUID_TCL_H
#define STATISTICS_FLUID_TCL_H
#include "parser.h"

#ifdef LB
/** Parser for fluid related analysis functions. */
int tclcommand_analyze_parse_fluid_cpu(Tcl_Interp *interp, int argc, char **argv);
#endif /* LB */

#ifdef LB_GPU
/** Parser for fluid related analysis functions. */
int tclcommand_analyze_parse_fluid_gpu(Tcl_Interp *interp, int argc, char **argv);
#endif

/** Parser for fluid related analysis functions. */
int tclcommand_analyze_parse_fluid(Tcl_Interp *interp, int argc, char **argv);
#endif /* STATISTICS_FLUID_TCL_H */
