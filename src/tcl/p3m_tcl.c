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
#include "parser.h"

#ifdef P3M
#include "p3m_tcl.h"
#include "p3m.h"

int tclcommand_inter_coulomb_parse_p3m_tune(Tcl_Interp * interp, int argc, char ** argv, int adaptive)
{
  int mesh = -1, cao = -1, n_interpol = -1;
  double r_cut = -1, accuracy = -1;

  while(argc > 0) {
    if(ARG0_IS_S("r_cut")) {
      if (! (argc > 1 && ARG1_IS_D(r_cut) && r_cut >= -1)) {
	Tcl_AppendResult(interp, "r_cut expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("mesh")) {
      if(! (argc > 1 && ARG1_IS_I(mesh) && mesh >= -1)) {
	Tcl_AppendResult(interp, "mesh expects an integer >= -1",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("cao")) {
      if(! (argc > 1 && ARG1_IS_I(cao) && cao >= -1 && cao <= 7)) {
	Tcl_AppendResult(interp, "cao expects an integer between -1 and 7",
			 (char *) NULL);
	return TCL_ERROR;
      } 

    } else if(ARG0_IS_S("accuracy")) {
      if(! (argc > 1 && ARG1_IS_D(accuracy) && accuracy > 0)) {
	Tcl_AppendResult(interp, "accuracy expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }

    } else if (ARG0_IS_S("n_interpol")) {
      if (! (argc > 1 && ARG1_IS_I(n_interpol) && n_interpol >= 0)) {
	Tcl_AppendResult(interp, "n_interpol expects an nonnegative integer", (char *) NULL);
	return TCL_ERROR;
      }
    }
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }
  p3m_set_tune_params(r_cut, mesh, cao, -1.0, accuracy, n_interpol);

  /* check for optional parameters */
  if (argc > 0) {
    if (tclcommand_inter_coulomb_parse_p3m_opt_params(interp, argc, argv) == TCL_ERROR)
      return TCL_ERROR;
  }

  /* do the tuning */
  char *log = NULL;
  if (p3m_adaptive_tune(&log) == ES_ERROR) {  
    Tcl_AppendResult(interp, log, "\nfailed to tune P3M parameters to required accuracy", (char *) NULL);
    if (log)
      free(log);
    return TCL_ERROR;
  }
  
  /* Tell the user about the tuning outcome */
  Tcl_AppendResult(interp, log, (char *) NULL);

  if (log)
    free(log);

  return TCL_OK;
}

int tclcommand_inter_coulomb_parse_p3m(Tcl_Interp * interp, int argc, char ** argv)
{
  double r_cut, alpha, accuracy = -1.0;
  int mesh[3], cao, i;
  IntList il;
  init_intlist(&il);

  if (argc < 1) {
    Tcl_AppendResult(interp, "expected: inter coulomb <bjerrum> p3m tune | <r_cut> { <mesh> | \\{ <mesh_x> <mesh_y> <mesh_z> \\} } <cao> [<alpha> [<accuracy>]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if (ARG0_IS_S("tune"))
    return tclcommand_inter_coulomb_parse_p3m_tune(interp, argc-1, argv+1, 0);

  if (ARG0_IS_S("tunev2"))
    return tclcommand_inter_coulomb_parse_p3m_tune(interp, argc-1, argv+1, 1);
      
  if(! ARG0_IS_D(r_cut))
    return TCL_ERROR;  

  if(argc < 3 || argc > 5) {
    Tcl_AppendResult(interp, "wrong # arguments: inter coulomb <bjerrum> p3m <r_cut> { <mesh> | \\{ <mesh_x> <mesh_y> <mesh_z> \\} } <cao> [<alpha> [<accuracy>]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(! ARG_IS_I(1, mesh[0])) {
    if( ! ARG_IS_INTLIST(1, il) || !(il.n == 3) ) {
      Tcl_AppendResult(interp, "integer or interger list of length 3 expected", (char *) NULL);
      return TCL_ERROR;
    } else {
      mesh[0] = il.e[0];
      mesh[1] = il.e[1];
      mesh[2] = il.e[2];
    }
  } else {
    mesh[1] = mesh[2] = mesh[0];
  }

  if(! ARG_IS_I(2, cao)) {
    Tcl_AppendResult(interp, "integer expected", (char *) NULL);
    return TCL_ERROR;
  }
	
  if(argc > 3) {
    if(! ARG_IS_D(3, alpha))
      return TCL_ERROR;
  }
  else {
    Tcl_AppendResult(interp, "Automatic p3m tuning not implemented.",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(argc > 4) {
    if(! ARG_IS_D(4, accuracy)) {
      Tcl_AppendResult(interp, "double expected", (char *) NULL);
      return TCL_ERROR;
    }
  }

  if ((i = p3m_set_params(r_cut, mesh, cao, alpha, accuracy)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "r_cut must be positive", (char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "mesh must be positive", (char *) NULL);
      break;
    case -3:
      Tcl_AppendResult(interp, "cao must be between 1 and 7 and less than mesh",
		       (char *) NULL);
      break;
    case -4:
      Tcl_AppendResult(interp, "alpha must be positive", (char *) NULL);
      break;
    case -5:
      Tcl_AppendResult(interp, "accuracy must be positive", (char *) NULL);
      break;
    default:;
      Tcl_AppendResult(interp, "unspecified error", (char *) NULL);
    }

    return TCL_ERROR;
  }

  return TCL_OK;
}


int tclcommand_inter_coulomb_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv)
{
  int i; double d1, d2, d3;

  Tcl_ResetResult(interp);

  while (argc > 0) {
    /* p3m parameter: inter */
    if (ARG0_IS_S("n_interpol")) {
      
      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG1_IS_I(i)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 INTEGER parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (p3m_set_ninterpol(i) == TCL_ERROR) {
	Tcl_AppendResult(interp, argv[0], " argument must be positive",
			 (char *) NULL);
	return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;
    }
    
    /* p3m parameter: mesh_off */
    else if (ARG0_IS_S("mesh_off")) {
      
      if(argc < 4) {
	Tcl_AppendResult(interp, argv[0], " needs 3 parameters",
			 (char *) NULL);
	return TCL_ERROR;
      }
	
      if ((! ARG_IS_D(1, d1)) ||
	  (! ARG_IS_D(2, d2)) ||
	  (! ARG_IS_D(3, d3)))
	{
	  Tcl_AppendResult(interp, argv[0], " needs 3 DOUBLE parameters",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      if (p3m_set_mesh_offset(d1, d2 ,d3) == TCL_ERROR)
	{
	  Tcl_AppendResult(interp, argv[0], " parameters have to be between 0.0 an 1.0",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      argc -= 4;
      argv += 4;
    }
    
    /* p3m parameter: epsilon */
    else if(ARG0_IS_S( "epsilon")) {

      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }

      if (ARG1_IS_S("metallic")) {
	d1 = P3M_EPSILON_METALLIC;
      }
      else if (! ARG1_IS_D(d1)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 DOUBLE parameter or \"metallic\"",
	                 (char *) NULL);
	return TCL_ERROR;
      }
	
      if (p3m_set_eps(d1) == TCL_ERROR) {
        Tcl_AppendResult(interp, argv[0], " There is no error msg yet!",
                         (char *) NULL);
        return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;	    
    }
    else {
      Tcl_AppendResult(interp, "Unknown coulomb p3m parameter: \"",argv[0],"\"",(char *) NULL);
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}

/*********************** miscelanea of functions *************************************/

int tclprint_to_result_p3m(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, p3m.params.r_cut, buffer);
  Tcl_AppendResult(interp, "p3m ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.params.mesh[0]);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.params.cao);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.params.alpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.params.accuracy, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  Tcl_AppendResult(interp, "} {coulomb epsilon ", (char *) NULL);
  if (p3m.params.epsilon == P3M_EPSILON_METALLIC)
    Tcl_AppendResult(interp, " metallic ", (char *) NULL);
  else {
    Tcl_PrintDouble(interp, p3m.params.epsilon, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  }
  sprintf(buffer,"%d",p3m.params.inter);
  Tcl_AppendResult(interp, "n_interpol ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.params.mesh_off[0], buffer);
  Tcl_AppendResult(interp, "mesh_off ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.params.mesh_off[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.params.mesh_off[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

/************************************************
 * Debug functions printing p3m structures 
 ************************************************/

#endif /* of P3M */

