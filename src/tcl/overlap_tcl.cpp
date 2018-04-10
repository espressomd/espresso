/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file overlap_tcl.cpp
 *
 *  Implementation of \ref overlap_tcl.hpp
 */
#include "utils.hpp"
#include "interaction_data.hpp"
#include "overlap_tcl.hpp"
#include "overlap.hpp"

#ifdef OVERLAPPED

int tclcommand_inter_parse_overlapped_bonded(Tcl_Interp *interp, int bond_type,
					     int argc, char **argv)
{
  OverlappedBondedInteraction overlap_type = OVERLAP_UNKNOWN;

  if (argc < 3 ) {
    Tcl_AppendResult(interp, "overlappedd needs two string parameter: "
		     "<type> <filename>", (char *) NULL);
    return (TCL_ERROR);
  }  

  if (ARG_IS_S(1,"bond"))     overlap_type = OVERLAP_BOND_LENGTH;
  if (ARG_IS_S(1,"angle"))    overlap_type = OVERLAP_BOND_ANGLE;
  if (ARG_IS_S(1,"dihedral")) overlap_type = OVERLAP_BOND_DIHEDRAL;
  if (overlap_type == OVERLAP_UNKNOWN) {
    Tcl_AppendResult(interp, "Unknown type of bonded overlapped interaction. Should be: "
		     "\"bond\" or \"angle\" or \"dihedral\"", (char *) NULL);
    return (TCL_ERROR);
  }

  switch (overlapped_bonded_set_params(bond_type, overlap_type, argv[2])) {
  case 1:
    Tcl_AppendResult(interp, "illegal bond type", (char *)NULL);
    return TCL_ERROR;
  case 2:
    Tcl_AppendResult(interp, "cannot open \"", argv[2], "\"", (char *)NULL);
    return TCL_ERROR;
  case 3:
    Tcl_AppendResult(interp, "the number of parameters is wrong in the file \"", argv[2], "\"", (char *)NULL);
    return TCL_ERROR;
  default:
    return TCL_OK;
  }
}

int tclprint_to_result_overlapIA(Tcl_Interp *interp,
				 Bonded_ia_parameters *params)
{
  switch (params->p.overlap.type) {
  case OVERLAP_BOND_LENGTH:
    Tcl_AppendResult(interp, "overlapped bond \"",params->p.overlap.filename,"\"",(char *) NULL);
    return TCL_OK;
  case OVERLAP_BOND_ANGLE:
    Tcl_AppendResult(interp, "overlapped angle \"",params->p.overlap.filename,"\"",(char *) NULL);
    return TCL_OK;
  case OVERLAP_BOND_DIHEDRAL:
    Tcl_AppendResult(interp, "overlapped dihedral \"",params->p.overlap.filename,"\"",(char *) NULL);
    return TCL_OK;
  default:
    fprintf(stderr, "INTERNAL ERROR: unexpected overlap bond type!\n");
    errexit();
  }
  return TCL_OK;
}

#endif

