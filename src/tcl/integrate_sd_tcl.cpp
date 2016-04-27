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
#include "global.hpp"

/**  Hand over integrate usage information to tcl interpreter. */
int tclcommand_integrate_sd_print_usage(Tcl_Interp *interp) 
{
  Tcl_AppendResult(interp, "Usage of tcl-command integrate_sd:\n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate_sd <INT n steps>' for integrating n steps \n", (char *)NULL);
  //Tcl_AppendResult(interp, "'integrate set' for printing integrator status \n", (char *)NULL);
  //Tcl_AppendResult(interp, "'integrate set nvt' for enabling NVT integration or \n" , (char *)NULL);
#ifdef NPT
  //Tcl_AppendResult(interp, "'integrate set npt_isotropic <DOUBLE p_ext> [<DOUBLE piston>] [<INT, INT, INT system_geometry>] [-cubic_box]' for enabling isotropic NPT integration \n" , (char *)NULL);
#endif
  return (TCL_ERROR);
}

int tclcommand_sd_set_particle_apart_print_usage(Tcl_Interp * interp)
{
  Tcl_AppendResult(interp, "Usage of tcl-command sd_set_particles_apart:\n", (char *)NULL);
  Tcl_AppendResult(interp, "  The command does not need any argument. It just \n\
  sets the particles at least 2*sd_radius away from each\n\
  other, so they dont overlap. Make sure to have enogh space\n\
  to reach this state fast ...\n",(char *) NULL);
  return (TCL_ERROR);
}

/** Hand over integrate status information to tcl interpreter. */
int tclcommand_integrate_sd_print_status(Tcl_Interp *interp) 
{
  return (TCL_OK);
  return (TCL_ERROR);
}

int tclcommand_sd_set_particles_apart(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  #ifndef SD
  return TCL_ERROR;
  #else
  if (argc != 1){
    Tcl_AppendResult(interp, "wrong # args: \n\"", (char *) NULL);
    return tclcommand_sd_set_particle_apart_print_usage(interp);
  }
  if (sd_radius < 0){
    Tcl_AppendResult(interp, "The particle radius for SD was not set\n (set with setmd sd_radius <value>).\n", (char *)NULL);
    return tclcommand_sd_set_particle_apart_print_usage(interp);
  }
  int status = sd_set_particles_apart();
  if (status == -5){
    Tcl_AppendResult(interp, "The Volumefraction is above the highest possible with hcp.\n\
Therefore this function cannot suceed ...\n", (char *) NULL);
    return (TCL_ERROR);
  }
  if (status == -4){
    Tcl_AppendResult(interp, "The Volumefraction is above 0.7 (and 0.7408 is the limit with hcp)\n\
Therefore this function will not suceed ...\n", (char *) NULL);
    return (TCL_ERROR);
  }
  if (status == 0){
    return (TCL_OK);
  }
  else {
    //Tcl_AppendResult(interp, "Unknown return code.\n", (char *) NULL);
    return (TCL_ERROR);
  }
  #endif
}

int tclcommand_sd_test(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  if (argc < 3){
    Tcl_AppendResult(interp, "wrong # args: \n\"", (char *) NULL);
    return tclcommand_sd_set_particle_apart_print_usage(interp);
  }
  int type, size;
  int status=0;
  ARG_IS_I(1,size);
  ARG_IS_I(2,type);
  printf("size: %d  type: %d\n", size, type);
#ifdef SD
  status = sd_test(size, type);
  printf("test returned: %d\n", status);
#else
  fprintf(stderr,"I didn't find SD - not compiled in\n");
#endif
  return status;
}

int tclcommand_integrate_sd(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  #if defined(SD) || defined(BD)
  int  n_steps;
  
  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if (argc < 1) {
    Tcl_AppendResult(interp, "wrong # args: \n\"", (char *) NULL);
    return tclcommand_integrate_sd_print_usage(interp);
  }
  else if (argc < 2) {
    return tclcommand_integrate_sd_print_status(interp);
    }

  if (ARG1_IS_S("set")) {
    if      (argc < 3){
      return tclcommand_integrate_sd_print_status(interp);
    }
    
    Tcl_AppendResult(interp, "unknown integrator method:\n", (char *)NULL);
    return tclcommand_integrate_sd_print_usage(interp);
      
  } else if ( !ARG_IS_I(1,n_steps) ) return tclcommand_integrate_sd_print_usage(interp);
  
  /* go on with integrate <n_steps> */
  if(n_steps <= 0) {
    Tcl_AppendResult(interp, "illegal number of steps (must be >0) \n", (char *) NULL);
    return tclcommand_integrate_sd_print_usage(interp);
  }
  /* perform integration */
  if (!correlations_autoupdate && !observables_autoupdate) {
    integrate_sd(n_steps);
  } else  {
    for (int i=0; i<n_steps; i++) {
      integrate_sd(1);
      autoupdate_observables();
      autoupdate_correlations();
    }
    if (n_steps == 0){
      integrate_sd(0);
    }
  }
  ;
  return gather_runtime_errors(interp, TCL_OK);
  #else
  return TCL_ERROR;
  #endif
}

/************************************************************/
/* Callback functions */
 

int tclcallback_sd_radius(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data <= 0) {
    if (-1.001 < data && data < -.999){
#if defined(SD) || defined(BD)
      Tcl_AppendResult(interp, "SD radius was set to -1. It must be set to a reasonable value bevor the integrate_sd can be called.", (char *) NULL);
#endif
      sd_radius=-1;
      mpi_bcast_parameter(FIELD_SD_RADIUS);
      return (TCL_OK);
    }
    Tcl_AppendResult(interp, "radius must be > 0.", (char *) NULL);
    return (TCL_ERROR);
  }
  sd_radius = data;
  mpi_bcast_parameter(FIELD_SD_RADIUS);
  return (TCL_OK);
}

int tclcallback_sd_viscosity(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data <= 0) {
    Tcl_AppendResult(interp, "viscosity must be > 0.", (char *) NULL);
    return (TCL_ERROR);
  }
  sd_viscosity = data;
  mpi_bcast_parameter(FIELD_SD_VISCOSITY);
  return (TCL_OK);
}

int tclcallback_sd_seed(Tcl_Interp *interp, void *_data)
{
  int * data = (int *)_data;
  sd_seed[0]=data[0];
  sd_seed[1]=data[1];
  mpi_bcast_parameter(FIELD_SD_SEED);
  return (TCL_OK);
}

int tclcallback_sd_random_state(Tcl_Interp *interp, void *_data)
{
  int * data = (int *)_data;
  sd_random_state[0]=data[0];
  sd_random_state[1]=data[1];
  mpi_bcast_parameter(FIELD_SD_RANDOM_STATE);
  return (TCL_OK);
}

int tclcallback_sd_random_precision(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data <= 0) {
    Tcl_AppendResult(interp, "precision_random must be > 0.", (char *) NULL);
    return (TCL_ERROR);
  }
  sd_random_precision = data;
  mpi_bcast_parameter(FIELD_SD_RANDOM_PRECISION);
  return (TCL_OK);
}
