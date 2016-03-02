/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _VIRTUAL_SITES_COM_TCL_H
#define _VIRTUAL_SITES_COM_TCL_H
#include "parser.hpp"

#ifdef VIRTUAL_SITES_COM

// Analyze the pressure on the molecule level
int tclcommand_analyze_parse_and_print_pressure_mol(Tcl_Interp *interp,int argc, char **argv);
// Analyze kinetic energy of the molecules
int tclcommand_analyze_parse_and_print_energy_kinetic_mol(Tcl_Interp *interp,int argc, char **argv);
// Sanity checks the positions of virtual sites
int tclcommand_analyze_parse_and_print_check_mol(Tcl_Interp *interp,int argc, char **argv);
// Analyze dipole moment on molecular basis
int tclcommand_analyze_parse_and_print_dipmom_mol(Tcl_Interp *interp,int argc, char **argv);
#endif

#endif
