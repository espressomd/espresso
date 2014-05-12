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
/** \file reaction_field_tcl.cpp
 *
 *  Implementation of \ref reaction_field_tcl.hpp
 */
#include "reaction_field_tcl.hpp"

#ifdef ELECTROSTATICS
#include "reaction_field.hpp"

int tclprint_to_result_rf(Tcl_Interp *interp, const char *name)
{
  char buffer[TCL_DOUBLE_SPACE];
  sprintf(buffer,"%s",name);
  Tcl_AppendResult(interp, buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, rf_params.kappa, buffer);
  Tcl_AppendResult(interp, buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, rf_params.epsilon1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, rf_params.epsilon2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, rf_params.r_cut, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

int tclcommand_inter_coulomb_parse_rf(Tcl_Interp * interp,
				      int argc, char ** argv,int method)
{
  double kappa,epsilon1,epsilon2, r_cut;
  int i;

  if(argc < 4) {
    Tcl_AppendResult(interp, "rf needs 4 parameters: "
                               "<kappa> <epsilon1> <epsilon2> <r_cut>",(char *) NULL);
    return TCL_ERROR;
  }

  coulomb.method = method;

  if ((! ARG_IS_D(0, kappa))      ||
      (! ARG_IS_D(1, epsilon1))   ||
      (! ARG_IS_D(2, epsilon2))   ||
      (! ARG_IS_D(3, r_cut)        )) {
      Tcl_AppendResult(interp, "rf needs 4 parameters: "
                               "<kappa> <epsilon1> <epsilon2> <r_cut>",(char *) NULL);
       return TCL_ERROR;
  }

  if ( (i = rf_set_params(kappa,epsilon1,epsilon2,r_cut)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "rf eps must be positive.",(char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "rf r_cut must be positive.",(char *) NULL);
      break;
    default:
      Tcl_AppendResult(interp, "unspecified error",(char *) NULL);
    }
    
    return TCL_ERROR;
  }

  return TCL_OK;
}

#ifdef INTER_RF

int tclcommand_inter_parse_interrf(Tcl_Interp * interp,
				   int part_type_a, int part_type_b,
				   int argc, char ** argv)
{
  /* parameters needed for RF */
  int rf_on;
  int change;

  /* get reaction_field interaction type */
  if (argc < 2) {
    Tcl_AppendResult(interp, "inter_rf needs 1 parameter: "
		     "<rf_on>",
		     (char *) NULL);
    return 0;
  }

  /* copy reaction_field parameters */
  if (! ARG_IS_I(1, rf_on)) {
    Tcl_AppendResult(interp, "<rf_on> must be int",
		     (char *) NULL);
    return 0;
  }
  change = 2;
	
  if (! ((rf_on==0) || (rf_on==1)) ) {
    Tcl_AppendResult(interp, "rf_on must be 0 or 1", (char *) NULL);
    return 0;
  }
  if (interrf_set_params(part_type_a, part_type_b,rf_on) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

int tclprint_to_result_interrfIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer,"%i",data->rf_on);
  Tcl_AppendResult(interp, "inter_rf ", buffer, " ",(char *) NULL);
  return TCL_OK;
}

#endif

#endif
