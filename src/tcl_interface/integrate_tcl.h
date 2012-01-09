/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef INTEGRATE_TCL_H
#define INTEGRATE_TCL_H
/** \file integrate.h    Molecular dynamics integrator.
 *
 *  For more information see \ref integrate.c "integrate.c".
*/   
#include <tcl.h>

#define INTEG_METHOD_NPT_ISO   0
#define INTEG_METHOD_NVT       1

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** tcl procedure for integrator steering. For documentation,
    see \ref tclcommand_integrate
*/
int tclcommand_integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);

/** Calculate maximal interaction range. 
    Uses \ref calc_maximal_cutoff.
    \ref max_range  = \ref max_cut + \ref #skin;
 */
int tclcallback_skin(Tcl_Interp *interp, void *_data);

/** Callback for integration time_step (0.0 <= time_step).
    \return TCL status.
*/
int tclcallback_time_step(Tcl_Interp *interp, void *_data);

/** Callback for current time in the integration.
    If no value is set the integration starts at time = 0.0.
    \return TCL status.
*/
int tclcallback_time(Tcl_Interp *interp, void *_data);

/** Implements the tcl-command 'invalidate_system' which forces a system re-init. 
    For more information, see \ref tclcommand_invalidate_system. */
int tclcommand_invalidate_system(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
