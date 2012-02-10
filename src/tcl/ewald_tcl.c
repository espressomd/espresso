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
/** \file ewald_tcl.c
 *
 *  Implementation of \ref ewald_tcl.h
 */
#include "ewald_tcl.h"
#include "ewald.h"
#include "grid.h"

#ifdef ELECTROSTATICS

int tclprint_to_result_EWALD(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  double b=ewald.kmax;

  Tcl_PrintDouble(interp, ewald.r_cut, buffer);
  Tcl_AppendResult(interp, "ewald ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, ewald.alpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, b, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

int tclcommand_inter_coulomb_parse_ewald(Tcl_Interp * interp, int argc, char ** argv)
{
  double r_cut, alpha;
  int i, kmax;

  coulomb.method = COULOMB_EWALD;
    
#ifdef PARTIAL_PERIODIC
  if(PERIODIC(0) == 0 ||
     PERIODIC(1) == 0 ||
     PERIODIC(2) == 0)
    {
      Tcl_AppendResult(interp, "Need periodicity (1,1,1) with Coulomb EWALD",
		       (char *) NULL);
      return TCL_ERROR;  
    }
#endif

  if (argc < 2) {
    Tcl_AppendResult(interp, "expected: inter coulomb <bjerrum> ewald <r_cut> <alpha> <kmax>",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(! ARG0_IS_D(r_cut))
    return TCL_ERROR;  

  if(argc != 3) {
    Tcl_AppendResult(interp, "wrong # arguments: inter coulomb <bjerrum> ewald <r_cut> <alpha> <kmax>",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(! ARG_IS_D(1, alpha))
    return TCL_ERROR;

  if(! ARG_IS_I(2, kmax))
    return TCL_ERROR;

  if ((i = ewald_set_params(r_cut, alpha, kmax)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "r_cut must be positive", (char *) NULL);
      break;
    case -4:
      Tcl_AppendResult(interp, "alpha must be positive", (char *) NULL);
      break;
    case -5:
      Tcl_AppendResult(interp, "kmax must be greater than zero", (char *) NULL);
    default:;
      Tcl_AppendResult(interp, "unspecified error", (char *) NULL);
    }

    return TCL_ERROR;

  }

  return TCL_OK;
}

#endif
