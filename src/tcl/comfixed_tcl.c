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
/** \file comfixed_tcl.c
 *
 *  Implementation of \ref comfixed_tcl.h
 */
#include "comfixed_tcl.h"
#include "parser.h"

#ifdef COMFIXED
#include "comfixed.h"
#include "interaction_data.h"

int tclprint_to_result_comfixedIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[ES_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer,"%d",data->COMFIXED_flag);
  Tcl_AppendResult(interp, "comfixed ", buffer, " ", (char *) NULL);
  
  return TCL_OK;
}

int tclcommand_inter_parse_comfixed(Tcl_Interp * interp,
				    int part_type_a, int part_type_b,
				    int argc, char ** argv)
{
  int flagc;
  	
  if (argc != 2) {
    Tcl_AppendResult(interp, "comfixed needs 1 parameters: "
		     "<comfixed_flag> ", (char *) NULL);
    return 1;
  }
	 
  if (part_type_a != part_type_b) {
    Tcl_AppendResult(interp, "comfixed must be among same type interactions", (char *) NULL);
    return 1;
  }

  /* copy comfixed parameters */
  if ((! ARG_IS_I(1, flagc)) ) {
    Tcl_AppendResult(interp, "comfixed needs 1 INTEGER parameter: "
		     "<comfixed_flag>", (char *) NULL);
    return 2;
  }

  switch (comfixed_set_params(part_type_a, part_type_b, flagc)) {
  case 1:
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 2;
  case 2:
    Tcl_AppendResult(interp, "works only with a single CPU", (char *) NULL);
    return 2;
  case 3:
    Tcl_AppendResult(interp, "works only with non periodic BC", (char *) NULL);
    return 2;
  }

   
   return 2;
}


#endif

