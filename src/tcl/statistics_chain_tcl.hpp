/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef STATISTICS_CHAIN_TCL_H
#define STATISTICS_CHAIN_TCL_H
#include "parser.hpp"

int tclcommand_analyze_parse_set_chains(Tcl_Interp *interp, int argc, char **argv);
int tclcommand_analyze_set_parse_chain_topology(Tcl_Interp *interp, int argc, char **argv);

///
int tclcommand_analyze_parse_re(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_rg(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_rh(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_internal_dist(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_bond_l(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_bond_dist(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_g123(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_g_av(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_formfactor(Tcl_Interp *interp, int average, int argc, char **argv);
///
int tclcommand_analyze_parse_rdfchain(Tcl_Interp *interp, int argc, char **argv);
#ifdef ELECTROSTATICS
///
int tclcommand_analyze_parse_cwvac(Tcl_Interp *interp, int argc, char **argv);
#endif

#endif
