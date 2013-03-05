/*
  Copyright (C) 2010,2011 Florian Fahrenberger
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

/** \file maggs.c
 *  Maxwell Equations Molecular Dynamics (MEMD) method for electrostatic
 *  interactions.
 *
 *  We use a local update scheme to propagate artificial B-fields on a
 *  lattice in the system. In principal, the algorithm simulates full
 *  electrodynamics, but with a tunable speed of light.
 *
 *  The method is very usable for large particle numbers or highly
 *  parallel architectures, since it is local and scales linearly.
 *  It is not suited for high-precision calculation of forces, since
 *  the simple interpolation scheme produces errors in the order of
 *  10^-5 in the force.
 *
 *  The chosen mesh should roughly be of the size of the particles.
 *
 *  Further reading on the algorithm:
 *  <ul>
 *  <li> I. Pasichnyk and B. Dunweg, Coulomb interaction via local dynamics: a molecular-dynamics algorithm. J. Phys: Condens. Matter, 16 ,p. 1399-4020, (2004).
 *  </ul>
 *  
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "interaction_data.h"
#include "maggs.h"
#include "parser.h"

#ifdef ELECTROSTATICS

/** parse TCL command.
    number of parameters is checked and maggs_set_parameters function is called.
    @return zero if successful
    @param interp  TCL interpreter handle
    @param argc    number of arguments given
    @param argv    array of arguments given
*/
int tclcommand_inter_coulomb_parse_maggs(Tcl_Interp * interp, int argc, char ** argv)
{
    int mesh;
    double f_mass;
    double epsilon = 1.0;
    int finite_epsilon_flag = 1;
	
    if(argc < 2) {
        Tcl_AppendResult(interp, "Not enough parameters: inter coulomb memd <f_mass> <mesh>", (char *) NULL);
        return TCL_ERROR;
    }
	
    if(! ARG_IS_D(0, f_mass))
        return TCL_ERROR;
	
    if(! ARG_IS_I(1, mesh)) {
        Tcl_AppendResult(interp, "integer expected", (char *) NULL);
        return TCL_ERROR;
    }
	
    if(argc > 4) {
        Tcl_AppendResult(interp, "Too many parameters: inter coulomb memd <f_mass> <mesh> [epsilon <eps>]", (char *) NULL);
        return TCL_ERROR;
    }
    if(argc == 3) {
        Tcl_AppendResult(interp, "Usage: inter coulomb memd <f_mass> <mesh> [epsilon <eps>]", (char *) NULL);
        return TCL_ERROR;
    }
    if(argc == 4) {
        if (ARG_IS_S(2, "epsilon")) {
            if(! (ARG_IS_D(3, epsilon) && epsilon > 0.0)) {
                Tcl_AppendResult(interp, "epsilon expects a positive double",
                             (char *) NULL);
                return TCL_ERROR;
            }
        }
    } else finite_epsilon_flag=0;

  coulomb.method = COULOMB_MAGGS;
	
  int res = maggs_set_parameters(coulomb.bjerrum, f_mass, mesh,
                                 finite_epsilon_flag, epsilon);
  switch (res) {
  case -1:
    Tcl_AppendResult(interp, "mass of the field is negative", (char *)NULL);
    return TCL_ERROR;
  case -2:
    Tcl_AppendResult(interp, "mesh must be positive", (char *) NULL);
    return TCL_ERROR;
  case ES_OK:
    return TCL_OK;
  }
  Tcl_AppendResult(interp, "unknown error", (char *) NULL);
  return TCL_ERROR;
}

int tclprint_to_result_Maggs(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
	
  Tcl_PrintDouble(interp, maggs.f_mass, buffer);
  Tcl_AppendResult(interp, "maggs ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",maggs.mesh);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL); 
	
  return TCL_OK;
}

#endif // ELECTROSTATICS
