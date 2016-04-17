/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file nemd.cpp

    For more information see \ref nemd.hpp
 */
#include "nemd.hpp"
#include <cstdio>
#include "integrate.hpp"
#include "communication.hpp"
#include "parser.hpp"
#include "grid.hpp"

/************************************************************/


/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Hand over nemd usage information to tcl interpreter. */
int tclcommand_nemd_print_usage(Tcl_Interp *interp) 
{
#ifdef NEMD
  Tcl_AppendResult(interp, "Usage of tcl command nemd:\n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd\" for returning the status or \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd off\" \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd exchange <INT n_slabs> <INT n_exchange>\" \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd shearrate <INT n_slabs> <DOUBLE shearrate>\" \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd profile\" for returning the velocity profile \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd viscosity\" for returning the viscosity \n", (char *)NULL);
  return (TCL_ERROR);
#else
  Tcl_AppendResult(interp, "nemd not compiled in!", (char *)NULL);
  return (TCL_ERROR);
#endif
}

#ifdef NEMD

/** Hand over nemd status information to tcl interpreter. */
int tclcommand_nemd_print_status(Tcl_Interp *interp) 
{
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];
  switch (nemd_method) {
  case NEMD_METHOD_OFF:
    Tcl_AppendResult(interp, "off", (char *)NULL);
    return (TCL_OK);
    break;
  case NEMD_METHOD_EXCHANGE:
    sprintf(buffer, "%d", nemddata.n_slabs);
    Tcl_AppendResult(interp, "exchange ",buffer, (char *)NULL);
    sprintf(buffer, "%d", nemddata.n_exchange);
    Tcl_AppendResult(interp, " ",buffer, (char *)NULL);
    return (TCL_OK);
    break;
  case NEMD_METHOD_SHEARRATE:
    sprintf(buffer, "%d", nemddata.n_slabs);
    Tcl_AppendResult(interp, "shearrate ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nemddata.shear_rate, buffer);
    Tcl_AppendResult(interp, " ",buffer, (char *)NULL);
    return (TCL_OK);
    break;
  default:
    return (TCL_ERROR);
  }
  return (TCL_ERROR);
}

/** Set nemd method to exchange and set nemd parameters */
int tclcommand_nemd_parse_exchange(Tcl_Interp *interp, int argc, char **argv) 
{
  int n_slabs, n_exchange;
  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args:  ", (char *)NULL);
    return tclcommand_nemd_print_usage(interp);
  }
  if ( !ARG_IS_I(2, n_slabs) || !ARG_IS_I(3, n_exchange) ) {
    Tcl_AppendResult(interp, "wrong argument type:  ", (char *)NULL);
    return tclcommand_nemd_print_usage(interp);
  }
  /* parameter sanity */
  if ( n_slabs<0 || n_slabs%2!=0 ) {
    Tcl_AppendResult(interp, "nemd <n_slabs> must be non negative and even!",(char *)NULL);
    return (TCL_ERROR);
  }  
  if ( n_slabs > 0 && n_exchange < 0 ) {
    Tcl_AppendResult(interp, "nemd <n_exchange> must be positive!",(char *)NULL);
    return (TCL_ERROR);
  } 

  nemd_method = NEMD_METHOD_EXCHANGE;
  nemd_init(n_slabs, n_exchange, 0.0);
  return (TCL_OK);
}

/** Set nemd method to shearrate and set nemd parameters */
int tclcommand_nemd_parse_shearrate(Tcl_Interp *interp, int argc, char **argv) 
{
  int n_slabs;
  double shearrate;
  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args:  ", (char *)NULL);
    return tclcommand_nemd_print_usage(interp);
  }
  if ( !ARG_IS_I(2, n_slabs) || !ARG_IS_D(3, shearrate) ) {
    Tcl_AppendResult(interp, "wrong argument type:  ", (char *)NULL);
    return tclcommand_nemd_print_usage(interp);
  }
 /* parameter sanity */
  if ( n_slabs<0 || n_slabs%2!=0 ) {
    Tcl_AppendResult(interp, "nemd <n_slabs> must be non negative and even!",(char *)NULL);
    return (TCL_ERROR);
  }  
  nemd_method = NEMD_METHOD_SHEARRATE;
  nemd_init(n_slabs, 0, shearrate);
  return (TCL_OK);
}

/** Hand over velocity profile to tcl interpreter */
int tclcommand_nemd_parse_and_print_profile(Tcl_Interp *interp)
{
  int i;
  double val;
  char buffer[TCL_DOUBLE_SPACE];
  
  INTEG_TRACE(fprintf(stderr,"%d: tclcommand_nemd_parse_and_print_profile:\n",this_node));
  if(nemd_method == NEMD_METHOD_OFF) {
    Tcl_AppendResult(interp, "nemd is off", (char *)NULL);
    return (TCL_OK);
  }
  
  for(i=0;i<nemddata.n_slabs;i++) {
    /* note: output velocities as usual have to be resacled by 1/time_step! */
    val = nemddata.velocity_profile[i]/(nemddata.profile_norm*time_step);
    Tcl_PrintDouble(interp, val, buffer);
    Tcl_AppendResult(interp," ", buffer, (char *)NULL);
    
    nemddata.velocity_profile[i] = 0.0;
  }
  
  nemddata.profile_norm = 0;
  return (TCL_OK);
}

/** Gives back the calculated viscosity */
int tclcommand_nemd_parse_and_print_viscosity(Tcl_Interp *interp)
{
  double shear_rate=0.0, mean_force, viscosity;
  char buffer[TCL_DOUBLE_SPACE];
  INTEG_TRACE(fprintf(stderr,"%d: tclcommand_nemd_parse_and_print_viscosity:\n",this_node));

  /* calculate shear_rate */
  switch (nemd_method) {
  case NEMD_METHOD_OFF:
    Tcl_AppendResult(interp, "nemd is off", (char *)NULL);
    return (TCL_OK);
  case NEMD_METHOD_EXCHANGE:
    shear_rate = 1.0;
  case NEMD_METHOD_SHEARRATE:
    shear_rate = nemddata.shear_rate;
  }
  /* rescale momentum exchange (vel_internal = time_step * vel) */
  nemddata.momentum /= time_step;
  /* Calculate average Force := Momentum transfer per time unit */
  mean_force = nemddata.momentum / (nemddata.momentum_norm*time_step);
  /* Calculate viscosity := mean_force / (shearrate * area) */
  viscosity = mean_force / (shear_rate*4.0*box_l[0]*box_l[1]);
  Tcl_PrintDouble(interp, viscosity, buffer);
  Tcl_AppendResult(interp,buffer, (char *)NULL);

  nemddata.momentum = 0.0;
  nemddata.momentum_norm = 0;
  return (TCL_OK);
}
#endif

/*@}*/

/************************************************************
 *            Exported Functions                            *
 ************************************************************/

int tclcommand_nemd(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
#ifdef NEMD
  int status = TCL_OK;

  INTEG_TRACE(fprintf(stderr,"%d: nemd:\n",this_node));
  Tcl_ResetResult(interp);

  /* print nemd status */
  if(argc == 1) {
    status = tclcommand_nemd_print_status(interp) ;
  }
  else if (ARG1_IS_S("off")) {
    nemd_method = NEMD_METHOD_OFF;
    status = nemd_free();
  }  
  else if (ARG1_IS_S("exchange")) {
    status = tclcommand_nemd_parse_exchange(interp,argc,argv);
  } 
  else if (ARG1_IS_S("shearrate")) {
    status = tclcommand_nemd_parse_shearrate(interp,argc,argv);
  } 
  else if (ARG1_IS_S("profile")) {
    status = tclcommand_nemd_parse_and_print_profile(interp);
  } 
  else if (ARG1_IS_S("viscosity")) {
    status = tclcommand_nemd_parse_and_print_viscosity(interp);
  } 
  else {
    Tcl_AppendResult(interp, "Unkwnown keyword: \n", (char *)NULL);
    return tclcommand_nemd_print_usage(interp);
  }

  return gather_runtime_errors(interp, status);

#endif
  INTEG_TRACE(fprintf(stderr,"%d: call to nemd but not compiled in!\n",this_node));
  return tclcommand_nemd_print_usage(interp);
}


