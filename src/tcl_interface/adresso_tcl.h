/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef ADRESSO_TCL_H
#define ADRESSO_TCL_H
/** \file adresso.h
    This is the place for adaptive resolution scheme (adress)
    Implementation of adresso.h

    For more details about adress see:
    - M. Praprotnik, L. Delle Site and K. Kremer, JCP 123, 224106, 2005. 
    - M. Praprotnik, L. Delle Site and K. Kremer, Annu. Rev. Phys. Chem. 59, 545-571, 2008. 
    - S. Poblete, M. Praprotnik, K. Kremer and L. Delle Site, J. Chem. Phys. 132, 114101, 2010. 

    For more detail about the implementation here see:
    - C. Junghans and S. Poblete, Comp. Phys. Comm. 181, 1449, 2010.
*/

#include <tcl.h>
#include "particle_data.h"
#include "virtual_sites.h"
#include "interaction_data.h"
#include "communication.h"



/** \name Exported Variables */
/************************************************************/
/*@{*/
extern double adress_vars[7];
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/
/** Implements the Tcl command "adress". This allows for seetings for adress
*/
int tclcommand_adress(ClientData data, Tcl_Interp *interp, int argc, char **argv);

int tclcommand_update_adress_weights(ClientData _data, Tcl_Interp * interp, int argc, char ** argv);

#ifdef ADRESS
// This code requires the "center of mass" implementation of virtual sites
#ifndef VIRTUAL_SITES_COM
 #error Adress requires the "center of mass"-implementation  of virtual sites. Please activate it in myconfig.h
#endif
/* #ifdef THERMODYNAMIC_FORCE */
int tclcommand_thermodynamic_force_parse_opt(Tcl_Interp * interp, int type, double prefactor, int argc, char ** argv);
int tclcommand_thermodynamic_force(ClientData _data, Tcl_Interp * interp, int argc, char ** argv);
/* #endif */



MDINLINE int adress_tab_parser(Tcl_Interp * interp,
			int part_type_a, int part_type_b,
			int argc, char ** argv)
{
  char *filename = NULL;

  /* adress_tab interactions should supply a file name for a file containing
     both force and energy profiles as well as number of points, max
     values etc.
  */
  if (argc < 2) {
    Tcl_AppendResult(interp, "tabulated potentials require a filename: "
		     "<filename>",
		     (char *) NULL);
    return 0;
  }

  /* copy tabulated parameters */
  filename = argv[1];

  switch (adress_tab_set_params(part_type_a, part_type_b, filename)) {
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
    Tcl_AppendResult(interp, "attempt to read file \"", filename,
		     "\" failed, could not find start the start token <#>", (char *)NULL);
    return 0;
  }
  return 2;
}



/* #endif */

#endif
/*@}*/
#endif
