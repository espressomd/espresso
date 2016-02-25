/*
  Copyright (C) 2014,2015,2016 The ESPResSo project
  
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
#ifndef MMM1D_GPU_TCL_H
#define MMM1D_GPU_TCL_H

#include "parser.hpp"
#include "config.hpp"

#ifdef MMM1D_GPU

#ifndef ELECTROSTATICS
#error MMM1D_GPU requires ELECTROSTATICS
#endif

/// print the mmm1d parameters to the interpreters result
int tclprint_to_result_MMM1DGPU (Tcl_Interp * interp);

/// parse the mmm1d parameters
int tclcommand_inter_coulomb_parse_mmm1dgpu (Tcl_Interp * interp, int argc, char **argv);

#endif
#endif
