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

#include "utils.hpp"

#ifdef TABULATED

#include "tab_tcl.hpp"
#include "tab.hpp"
#include "parser.hpp"
#include "forcecap_tcl.hpp"

/// parse parameters for the tabulated bonded potential
int tclcommand_inter_parse_tabulated_bonded(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  int tab_type = TAB_UNKNOWN;

  if (argc < 3 ) {
    Tcl_AppendResult(interp, "tabulated needs two string parameter: "
		     "<type> <filename>", (char *) NULL);
    return (TCL_ERROR);
  }  

  if (ARG_IS_S(1,"bond"))     tab_type = TAB_BOND_LENGTH;
  if (ARG_IS_S(1,"angle"))    tab_type = TAB_BOND_ANGLE;
  if (ARG_IS_S(1,"dihedral")) tab_type = TAB_BOND_DIHEDRAL;
  if (tab_type == TAB_UNKNOWN) {
    Tcl_AppendResult(interp, "Unknown type of bonded tabulated interaction. Should be: "
		     "\"bond\" or \"angle\" or \"dihedral\"", (char *) NULL);
    return (TCL_ERROR);
  }

  switch (tabulated_bonded_set_params(bond_type, tab_type, argv[2])) {
  case 1:
    Tcl_AppendResult(interp, "illegal bond type", (char *)NULL);
    return TCL_ERROR;
  case 3:
    Tcl_AppendResult(interp, "cannot open \"", argv[2], "\"", (char *)NULL);
    return TCL_ERROR;
  case 4:
    Tcl_AppendResult(interp, "attempt to read file \"", argv[2], "\" failed. "
		     "Could not find start the start token <#>", (char *)NULL);
    return TCL_ERROR;
  case 5:
    Tcl_AppendResult(interp, "attempt to read file \"", argv[2], "\" failed. "
		     "Could not understand some numbers", (char *)NULL);
    return TCL_ERROR;
  case 6:
    if (tab_type == TAB_BOND_ANGLE) {
      Tcl_AppendResult(interp, "bond angle potential has to be defined in the interval 0 to pi", (char *)NULL);
    } else {
      Tcl_AppendResult(interp, "dihedral potential has to be defined in the interval 0 to 2pi", (char *)NULL);
    }
    return TCL_ERROR;
  default:
    return TCL_OK;
  }
}

/// parser for the force cap
int tclcommand_inter_parse_tabforcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  if (argc==1) {
    fprintf(stderr, "WARNING: \"inter tabforcecap\" is deprecated "
	    "and will be removed in some further version. "
	    "Use \"inter forcecap\" instead.\n");
  }
  return tclcommand_inter_parse_forcecap(interp, argc, argv);
}

int tclcommand_inter_parse_tab(Tcl_Interp * interp,
		int part_type_a, int part_type_b,
		int argc, char ** argv)
{
char *filename = NULL;

/* tabulated interactions should supply a file name for a file containing
both force and energy profiles as well as number of points, max
values etc.
*/
if (argc < 2) {
Tcl_AppendResult(interp, "tabulated potentials require a filename: "
	     "<filename>",
	     (char *) NULL);
return TCL_ERROR;
}

/* copy tabulated parameters */
filename = argv[1];

switch (tabulated_set_params(part_type_a, part_type_b, filename)) {
case 1:
Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
return 0;
case 2:
Tcl_AppendResult(interp, "the length of the filename must be less than 256 characters,"
	     "but is \"", filename, "\"", (char *)NULL);
return 0;
case 3:
Tcl_AppendResult(interp, "cannot open \"", filename, "\"", (char *)NULL);
return 0;
case 4:
Tcl_AppendResult(interp, "attempt to read file \"", filename, "\" failed. "
	     "Could not find start the start token <#>", (char *)NULL);
return 0;
case 5:
Tcl_AppendResult(interp, "attempt to read file \"", filename, "\" failed. "
	     "Could not understand some numbers", (char *)NULL);
return TCL_ERROR;
case 6:
Tcl_AppendResult(interp, "number of data points does not match the existing table", (char *)NULL);
return 0;

}
return 2;
}

int tclprint_to_result_tabulated_bondedIA(Tcl_Interp *interp,
					  Bonded_ia_parameters *params)
{
  switch (params->p.tab.type) {
  case TAB_BOND_LENGTH:
    Tcl_AppendResult(interp, "tabulated bond \"",params->p.tab.filename,"\"",(char *) NULL);
    return TCL_OK;
  case TAB_BOND_ANGLE:
    Tcl_AppendResult(interp, "tabulated angle \"",params->p.tab.filename,"\"",(char *) NULL);
    return TCL_OK;
  case TAB_BOND_DIHEDRAL:
    Tcl_AppendResult(interp, "tabulated dihedral \"",params->p.tab.filename,"\"",(char *) NULL);
    return TCL_OK;
  }
  Tcl_AppendResult(interp, "unknown type of tabulated bonded interaction encountered",(char *) NULL);
  return TCL_ERROR;
}

#endif

