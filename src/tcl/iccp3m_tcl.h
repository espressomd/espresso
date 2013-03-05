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
#ifndef _ICCP3M_TCL_H
#define _ICCP3M_TCL_H

#if defined(ELECTROSTATICS)
#include "parser.h"

/** Implementation of the tcl-command <br>
    iccp3m  { \<last_ind_id\> \<e1\> \<num_iteration\> \<convergence\> \<relaxation\> \<area\> \<normal_components\> \<e_in/e_out\>  [\<ext_field\>] | iterate } 
    ICC sets up and calculates induced charges on dielectric surfaces. At the beginning of every simulation run particles on the surface boundary 
    have to be set up (before any real particle) together with the list of areas, normal vectors and dielectric constant associated with them. 
    After that the iterate flag can be used during the simulation to update the value of the induced charges.
    
    Parameters: <br>
                 \<last_ind_id\> ID of the last surface charge. Note that the IDs of the surface charges must range from 0 to \<last_ind_id\>
                 \<e1\>          = Dielectric Constant of the Bulk accessible to free particles 
                 \<num_iteration\> = Maximum number of ICCP3M iterations calculating the induced charges. 
                 \<relaxation\> = Relaxaxion parameter \f$omega\f$ for the successive over-relaxation scheme. 
                 \<area\>       = List of the areas of each surface element.
                 \<normal_components\> = List of normal vectors of each surface element. 3n list entries. Do not have to be normalized.
                 \<e_in/e_out\> = Ratio of dielectric co


                 iterate         = Indicates that a previous surface discretization shall be used. T
*/
int tclcommand_iccp3m(ClientData data, Tcl_Interp *interp, int argc, char **argv);

#endif /* ELECTROSTATICS */

#endif /* ICCP3M_TCL_H */
