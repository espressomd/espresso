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
/** \file mmm1d.h MMM1D algorithm for long range coulomb interactions.
    Implementation of the MMM1D method for the calculation of the
    electrostatic interaction in one dimensionally periodic
    systems. For details on the method see MMM in general. The MMM1D
    method works only with the nsquared \ref tclcommand_cellsystem
    "cell system", since neither the near nor far formula can be
    decomposed. However, this implementation is reasonably fast, so
    that one can use up to 200 charges easily in a simulation.  */
#ifndef MMM1D_TCL_H
#define MMM1D_TCL_H

#include "parser.h"

#ifdef ELECTROSTATICS

/// print the mmm1d parameters to the interpreters result
int tclprint_to_result_MMM1D(Tcl_Interp *interp);

/// parse the mmm1d parameters
int tclcommand_inter_coulomb_parse_mmm1d(Tcl_Interp *interp, int argc, char **argv);

#endif
#endif
