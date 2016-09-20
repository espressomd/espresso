/*
  Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
  
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

/** \file oif_local_forces.hpp routines to calculate the local elastic forces
 for a particle quadruple (two triangles that have 1 edge in common)
 */

#include "tcl/parser.hpp"
#include "object-in-fluid/oif_local_forces.hpp"
#include "oif_local_forces_tcl.hpp"
#include "interaction_data.hpp"

// parse parameters for the oif_local_forces potential
int tclcommand_inter_parse_oif_local_forces(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
    double r0, ks, kslin, phi0, kb, A01, A02, kal;
    
    if (argc != 9) {
        Tcl_AppendResult(interp, "oif_global_forces needs 8 parameters: "
                         "<r0> <ks> <kslin> <phi0> <kb> <A01> <A02> <kal>", (char *) NULL);
        return (TCL_ERROR);
    }
    
    if ((! ARG_IS_D(1, r0)) || (! ARG_IS_D(2, ks)) || (! ARG_IS_D(3, kslin)) || (! ARG_IS_D(4, phi0)) || (! ARG_IS_D(5, kb)) || (! ARG_IS_D(6, A01)) || (! ARG_IS_D(7, A02)) || (! ARG_IS_D(8, kal)))
    {
        Tcl_AppendResult(interp, "oif_global_forces needs 7 DOUBLE parameters: "
                         "<r0> <ks> <kslin> <phi0> <kb> <A01> <A02> <kal>", (char *) NULL);
        return TCL_ERROR;
    }
    
    CHECK_VALUE(oif_local_forces_set_params(bond_type, r0, ks, kslin, phi0, kb, A01, A02, kal), "bond type must be nonnegative");
}

int tclprint_to_result_oiflocalforcesIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
    char buffer[TCL_DOUBLE_SPACE];
    Tcl_PrintDouble(interp, params->p.oif_local_forces.r0, buffer);
    Tcl_AppendResult(interp, "OIF_LOCAL_FORCES ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_local_forces.ks, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_local_forces.kslin, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_local_forces.phi0, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_local_forces.kb, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_local_forces.A01, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_local_forces.A02, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_local_forces.kal, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    return (TCL_OK);
}

