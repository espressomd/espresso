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
#ifndef STATISTICS_TCL_H
#define STATISTICS_TCL_H
#include "parser.h"

/** Implements the Tcl command \ref tclcommand_analyze. This allows for basic system analysis,
    both online and offline.
*/
int tclcommand_analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv);


/** return the approx diffusion constant of a special type of particle FIXME: this is not a smart way to compute D (very error-prone)!
 *  \param interp  TCL interpreter handle
 *  \param type_m  type of the particle, -1 for all
 *  \param n_time_steps number of timestep between saved configurations
 *  \param n_conf  number of saved contributions taken into account
 */
double tclcommand_analyze_print_MSD(Tcl_Interp *interp,int type_m, int n_time_steps,int n_conf);

#endif
