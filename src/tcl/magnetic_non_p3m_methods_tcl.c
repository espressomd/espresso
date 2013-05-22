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
/** \file magnetic_non_p3m_methods.c  All 3d non P3M methods to deal with the magnetic dipoles
 *   
 *  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *   Handling of a system of dipoles where no replicas 
 *   Assumes minimum image convention for those axis in which the 
 *   system is periodic as defined by setmd.
 *
 *   MDDS => Calculates dipole-dipole interaction of a perioidic system
 *   by explicitly summing the dipole-dipole interaction over several copies of the system
 *   Uses spherical summation order
 *
 */

#include "parser.h"

#ifdef DIPOLES
#include "magnetic_non_p3m_methods_tcl.h"
#include "magnetic_non_p3m_methods.h"

/* =============================================================================
                  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA                
   =============================================================================
*/

int  tclprint_to_result_DAWAANR(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, " dawaanr ", (char *) NULL);
  return TCL_OK;
}

/************************************************************/

int tclcommand_inter_magnetic_parse_dawaanr(Tcl_Interp * interp, int argc, char ** argv)
{
  if (dawaanr_set_params() != ES_OK) {
    Tcl_AppendResult(interp, "sorry: DAWAANR only works with 1 cpu", (char *) NULL);
    return TCL_ERROR;
  }
  return TCL_OK;
}

/************************************************************/



/************************************************************/

/*=================== */
/*=================== */
/*=================== */
/*=================== */
/*=================== */
/*=================== */

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS               
   =============================================================================
*/

extern int Ncut_off_magnetic_dipolar_direct_sum;

int tclprint_to_result_Magnetic_dipolar_direct_sum_(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_AppendResult(interp, " mdds ", buffer, (char *) NULL);
  Tcl_PrintDouble(interp,Ncut_off_magnetic_dipolar_direct_sum , buffer);
  Tcl_AppendResult(interp, " ", buffer, (char *) NULL);

  return TCL_OK;
}

/************************************************************/

int tclcommand_inter_magnetic_parse_mdds(Tcl_Interp * interp, int argc, char ** argv)
{
  int  n_cut=-1;
   
  while(argc > 0) {
    if (ARG0_IS_S("n_cut")) {
      if (! (argc > 1 && ARG1_IS_I(n_cut) && n_cut >= 0)) {
	Tcl_AppendResult(interp, "n_cut expects an nonnegative integer",
			 (char *) NULL);
	return TCL_ERROR;
      }
    } else { /* unknown parameter*/
      Tcl_AppendResult(interp, "unknown parameter/s for the magnetic dipolar direct sum, the only one accepted is:  n_cut  positive_integer", (char *) NULL);
      return TCL_ERROR;
    }
    
    argc -= 2;
    argv += 2;
  }

  if (mdds_set_params(n_cut) != ES_OK) {
    Tcl_AppendResult(interp, "sorry: magnetic dipolar direct sum only works with 1 cpu", (char *) NULL);
    return TCL_ERROR;
  }
  return TCL_OK;
}


#endif
