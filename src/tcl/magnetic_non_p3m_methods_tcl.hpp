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
#ifndef MAG_NON_P3M_TCL_H
#define MAG_NON_P3M_TCL_H
#include "parser.hpp"

#ifdef DIPOLES

/* =============================================================================
                  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA                
   =============================================================================
*/

/*  Information about the status of the method */
int tclprint_to_result_DAWAANR(Tcl_Interp *interp);
         
/* Parsing function for the dawaanr method*/
int tclcommand_inter_magnetic_parse_dawaanr(Tcl_Interp * interp, int argc, char ** argv);

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS               
   =============================================================================
*/

/*  Information about the status of the method */
int tclprint_to_result_Magnetic_dipolar_direct_sum_(Tcl_Interp *interp);


/* Parsing function for the magnetic dipolar direct sum method*/
int tclcommand_inter_magnetic_parse_mdds(Tcl_Interp * interp, int argc, char ** argv);

/* Sanity checks for the magnetic dipolar direct sum*/
int magnetic_dipolar_direct_sum_sanity_checks();

#endif /*of ifdef DIPOLES  */
#endif /* of ifndef  MAG_NON_P3M_H */
