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
/** \file mol_cut_tcl.c
 *
 *  Implementation of \ref mol_cut_tcl.h
 */
#include "parser.h"

#ifdef MOL_CUT
#include "mol_cut.h"
#include "interaction_data.h"

int tclprint_to_result_molcutIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer,"%i",data->mol_cut_type);
  Tcl_AppendResult(interp, "molcut ", buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, data->mol_cut_cutoff, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  return TCL_OK;
}

int tclcommand_inter_parse_molcut(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for molcut */
  int mol_cut_type;
  double mol_cut_cutoff;
  int change;

  /* get interaction type */
  if (argc < 3) {
    Tcl_AppendResult(interp, "molcut needs 2 parameter: "
		     "<type> <cutoff>",
		     (char *) NULL);
    return 0;
  }

  /* copy parameters */
  if (! ARG_IS_I(1, mol_cut_type)) {
    Tcl_AppendResult(interp, "type must be int",
		     (char *) NULL);
    return 0;
  }
  if (! ARG_IS_D(2, mol_cut_cutoff)) {
    Tcl_AppendResult(interp, "cutoff must be double",
		     (char *) NULL);
    return 0;
  }
  change = 3;
	
  if (! ((mol_cut_type==0) || (mol_cut_type==1)) ) {
    Tcl_AppendResult(interp, "type must be 0 or 1", (char *) NULL);
    return 0;
  }
  if (molcut_set_params(part_type_a, part_type_b,mol_cut_type,mol_cut_cutoff) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif
