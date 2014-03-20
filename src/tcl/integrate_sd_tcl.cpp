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

/** \file integrate.cpp   Molecular dynamics integrator.
 *
 *  For more information about the integrator 
 *  see \ref integrate.hpp "integrate.h".
*/

#include "integrate_sd.hpp"
#include "npt.hpp"
#include "interaction_data.hpp"
#include "lb.hpp"
#include "pressure_tcl.hpp"
#include "communication.hpp"
#include "parser.hpp"
#include "statistics_correlation.hpp"

/**  Hand over integrate usage information to tcl interpreter. */
int tclcommand_integrate_sd_print_usage(Tcl_Interp *interp) 
{
  Tcl_AppendResult(interp, "Usage of tcl-command integrate:\n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate <INT n steps>' for integrating n steps \n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate set' for printing integrator status \n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate set nvt' for enabling NVT integration or \n" , (char *)NULL);
#ifdef NPT
  Tcl_AppendResult(interp, "'integrate set npt_isotropic <DOUBLE p_ext> [<DOUBLE piston>] [<INT, INT, INT system_geometry>] [-cubic_box]' for enabling isotropic NPT integration \n" , (char *)NULL);
#endif
  return (TCL_ERROR);
}

/** Hand over integrate status information to tcl interpreter. */
int tclcommand_integrate_sd_print_status(Tcl_Interp *interp) 
{
  return (TCL_OK);
  return (TCL_ERROR);
}


int tclcommand_integrate_sd(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int  n_steps;
  
  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if (argc < 1) {
    Tcl_AppendResult(interp, "wrong # args: \n\"", (char *) NULL);
    return tclcommand_integrate_sd_print_usage(interp);  }
  else if (argc < 2) {                    return tclcommand_integrate_sd_print_status(interp); }

  if (ARG1_IS_S("set")) {
    if      (argc < 3)                    return tclcommand_integrate_sd_print_status(interp);
    
    Tcl_AppendResult(interp, "unknown integrator method:\n", (char *)NULL);
    return tclcommand_integrate_sd_print_usage(interp);
      
  } else if ( !ARG_IS_I(1,n_steps) ) return tclcommand_integrate_sd_print_usage(interp);

  /* go on with integrate <n_steps> */
  if(n_steps < 0) {
    Tcl_AppendResult(interp, "illegal number of steps (must be >0) \n", (char *) NULL);
    return tclcommand_integrate_sd_print_usage(interp);;
  }
  /* perform integration */
  integrate_sd(n_steps);
  return gather_runtime_errors(interp, TCL_OK);
}

/************************************************************/
/* Callback functions */
 
