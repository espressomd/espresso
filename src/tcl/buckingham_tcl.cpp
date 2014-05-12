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
/** \file buckingham_tcl.cpp
 *
 *  Implementation of \ref buckingham_tcl.hpp
 */
#include "buckingham_tcl.hpp"
#include "buckingham.hpp"
#include "forcecap_tcl.hpp"

#ifdef BUCKINGHAM

int tclprint_to_result_buckIA(Tcl_Interp *interp, int i, int j)
{
    char buffer[TCL_DOUBLE_SPACE];
    IA_parameters *data = get_ia_param(i, j);

    Tcl_PrintDouble(interp, data->BUCK_A, buffer);
    Tcl_AppendResult(interp, "buckingham ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_B, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_C, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_D, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_discont, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_shift, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_capradius, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_F1, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_F2, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
 
    return TCL_OK;
}

int tclcommand_inter_parse_buckingham(Tcl_Interp * interp,
				      int part_type_a, int part_type_b,
				      int argc, char ** argv)
{
  /* parameters needed for buckingham */
  double buck_A,buck_B,buck_C,buck_D,buck_cut,buck_discont,buck_shift,buck_cap_radius,F1,F2;
  int change;


  /* get buckingham interaction type */
  if (argc < 8) {
     Tcl_AppendResult(interp, "buckingham needs 7 parameters: "
		      "<buck_A> <buck_B> <buck_C> <buck_D> <buck_cut> <buck_discontinuity> <buck_shift> ",
		      (char *) NULL);
     return 0;
  }
  /* copy buckingham parameters */
  if ((! ARG_IS_D(1, buck_A))   ||
      (! ARG_IS_D(2, buck_B))   ||
      (! ARG_IS_D(3, buck_C))   ||
      (! ARG_IS_D(4, buck_D))   ||
      (! ARG_IS_D(5, buck_cut)) ||
      (! ARG_IS_D(6, buck_discont)) ||
      (! ARG_IS_D(7, buck_shift))) {
    Tcl_AppendResult(interp, "buckingham needs 7 DOUBLE parameters: "
           "<buck_A> <buck_B> <buck_C> <buck_D> <buck_cut> <buck_discontinuity> <buck_shift> ",
           (char *) NULL);
      return 0;
    }
    change = 8;
    buck_cap_radius = -1.0;
    F1 = 0.0;
    F2 = 0.0;
    /* check whether there are additional doubles, cap radius, F1 and F2*/
    if (argc >= 9 && ARG_IS_D(8, buck_cap_radius))
    {
      change++;
      if(argc >= 10 && ARG_IS_D(9, F1))
      {
         change++;
         if(argc >= 11 && ARG_IS_D(10, F1))
         change++;
      }
    }
    else
      Tcl_ResetResult(interp);
    if(buck_discont>buck_cut)
    {
       Tcl_AppendResult(interp, "ERROR: <buck_cut> must be greater than <buck_discontinuity>",
            (char *) NULL);
       return 0;
    }
    if (buckingham_set_params(part_type_a, part_type_b,
                  buck_A, buck_B, buck_C, buck_D,
            buck_cut, buck_discont, buck_shift,
            buck_cap_radius, F1, F2) == ES_ERROR) {
       Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
       return 0;
    }
    return change;
  }

  int tclcommand_inter_parse_buckforcecap(Tcl_Interp * interp, int argc, char ** argv)
  {
    if (argc==1) {
      fprintf(stderr, "WARNING: \"inter buckforcecap\" is deprecated "
	      "and will be removed in some further version. "
	      "Use \"inter forcecap\" instead.\n");
    }
  return tclcommand_inter_parse_forcecap(interp, argc, argv);
}

#endif
