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
#ifndef _INTEGRATE_SD_TCL_H
#define _INTEGRATE_SD_TCL_H
#include "parser.hpp"

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** tcl procedure for integrator steering. For documentation,
    see \ref tclcommand_integrate
*/
int tclcommand_integrate_sd(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);

/** tcl procedure to set particles to an non overlapping state
 */
int tclcommand_sd_set_particles_apart(ClientData data,
				     Tcl_Interp *interp,
				     int argc, char **argv);

/** tcl procedure to test some internal stuff in sd integrater
 */
int tclcommand_sd_test(ClientData data, Tcl_Interp *interp,
		       int argc, char **argv);

/** Callback for the sd particle radius (0.0 < sd_radius)
 */
int tclcallback_sd_radius(Tcl_Interp *interp, void *_data);

/** Callback for the sd viscosity (0.0 < sd_viscosity)
 */
int tclcallback_sd_viscosity(Tcl_Interp *interp, void *_data);

/** Callback for the sd random precision (0.0 < sd_random_precision)
 */
int tclcallback_sd_random_precision(Tcl_Interp *interp, void *_data);

/** Callback for the sd random state
 */
int tclcallback_sd_random_state(Tcl_Interp *interp, void *_data);

/** Callback for the sd random seed
 */
int tclcallback_sd_seed(Tcl_Interp *interp, void *_data);
/*@}*/

#endif
