/*
  Copyright (C) 2010,2011,2012,2013,2016 The ESPResSo project
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
/** \file affinity_tcl.cpp
 *
 *  Implementation of \ref affinity_tcl.hpp
 */
#include "affinity_tcl.hpp"
#include "../../core/object-in-fluid/affinity.hpp"

#ifdef AFFINITY

int tclprint_to_result_affinityIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer, "%d", data->affinity_type);
  Tcl_AppendResult(interp, "affinity ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->affinity_kappa, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->affinity_r0, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->affinity_Kon, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->affinity_Koff, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->affinity_maxBond, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->affinity_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  return TCL_OK;
}

int tclcommand_inter_parse_affinity(Tcl_Interp * interp,
				int part_type_a, int part_type_b,
				int argc, char ** argv)
{
  /* parameters needed for affinity */
  double kappa, r0, Kon, Koff, maxBond, cut;
  int change;
  int type;

  /* get affinity interaction type */
  if (argc < 8) {
    Tcl_AppendResult(interp, "affinity potential needs 7 parameters: "
		     "<affinity_type> <affinity_kappa> <affinity_r0> <affinity_Kon> <affinity_Koff> <affinity_maxBond> <affinity_cut>",
		     (char *) NULL);
    return 0;
  }

  /* copy affinity parameters */
  if ((! ARG_IS_I(1, type))	||
      (! ARG_IS_D(2, kappa))		||
      (! ARG_IS_D(3, r0))		||
      (! ARG_IS_D(4, Kon))		||
      (! ARG_IS_D(5, Koff))		||
      (! ARG_IS_D(6, maxBond))	||
      (! ARG_IS_D(7, cut)   )) {
    Tcl_AppendResult(interp, "affinity potential needs 7 parameters: "
		     "<affinity_type> <affinity_kappa> <affinity_r0> <affinity_Kon> <affinity_Koff> <affinity_maxBond> <affinity_cut>",
		     (char *) NULL);
    return 0;
  }
  change = 8;

  Tcl_ResetResult(interp);
  if ( maxBond <= r0 && (type == 1 || type == 3)) {
    Tcl_AppendResult(interp, "tcl affinity parser: Caution, affinity_maxBond should be greater than affinity_r0", (char *) NULL);
    return 0;
  }

  if ( cut <= maxBond && (type == 1 || type == 3)) {
    Tcl_AppendResult(interp, "tcl affinity parser: affinity_cut must be greater than affinity_maxBond. It is reccommended to use cut >= maxBond + epsilon, epsilon depends on the spead of the system, 0.5 should be ok.", (char *) NULL);
    return 0;
  }

  if (affinity_set_params(part_type_a, part_type_b,
                             type, kappa, r0, Kon, Koff, maxBond, cut) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif

