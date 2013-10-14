/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file thermostat.c
    Implementation of \ref thermostat.h "thermostat.h"
 */
#include <cmath>
#include "utils.hpp"

#include "communication.hpp"
#include "lattice.hpp"
#include "npt.hpp"
#include "ghmc.hpp"

#include "particle_data.hpp"
#include "parser.hpp"
#include "random.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "cells.hpp"
#include "lb.hpp"
#include "dpd.hpp"
#include "dpd_tcl.hpp"
#include "virtual_sites.hpp"
#include "thermostat_tcl.hpp"


int tclcommand_thermostat_parse_off(Tcl_Interp *interp, int argc, char **argv) 
{
  /* set temperature to zero */
  temperature = 0;
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  /* langevin thermostat */
  langevin_gamma = 0;
  mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA);
  /* dpd thermostat */
#ifdef DPD
  dpd_switch_off();
#endif
#ifdef INTER_DPD
  inter_dpd_switch_off();
#endif
#ifdef NPT
  /* npt isotropic thermostat */
  nptiso_gamma0 = 0;
  mpi_bcast_parameter(FIELD_NPTISO_G0);
  nptiso_gammav = 0;
  mpi_bcast_parameter(FIELD_NPTISO_GV);
#endif
#ifdef GHMC
  /* ghmc thermostat */
  ghmc_nmd = 1;
  mpi_bcast_parameter(FIELD_GHMC_NMD);
  ghmc_phi = 1;
  mpi_bcast_parameter(FIELD_GHMC_PHI);
  ghmc_mflip = GHMC_MFLIP_OFF;
  mpi_bcast_parameter(FIELD_GHMC_FLIP);
  ghmc_tscale = GHMC_TSCALE_OFF;
  mpi_bcast_parameter(FIELD_GHMC_SCALE);
#endif
  /* switch thermostat off */
  thermo_switch = THERMO_OFF;
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  return (TCL_OK);
}

int tclcommand_thermostat_parse_langevin(Tcl_Interp *interp, int argc, char **argv) 
{
  double temp, gamma;

  /* check number of arguments */
  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
		     argv[0]," ",argv[1]," <temp> <gamma>\"", (char *)NULL);
    return (TCL_ERROR);
  }

  /* check argument types */
  if ( !ARG_IS_D(2, temp) || !ARG_IS_D(3, gamma)) {
    Tcl_AppendResult(interp, argv[0]," ",argv[1]," needs two DOUBLES", (char *)NULL);
    return (TCL_ERROR);
  }

  if (temp < 0 || gamma < 0) {
    Tcl_AppendResult(interp, "temperature and gamma must be positive", (char *)NULL);
    return (TCL_ERROR);
  }

  /* broadcast parameters */
  temperature = temp;
  langevin_gamma = gamma;
  thermo_switch = ( thermo_switch | THERMO_LANGEVIN );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA);
  return (TCL_OK);
}

#ifdef NPT
int tclcommand_thermostat_parse_npt_isotropic(Tcl_Interp *interp, int argc, char **argv) 
{
  double temp, gamma0, gammav;
  /* check number of arguments */
  if (argc < 5) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
		     argv[0]," set ",argv[1]," <temp> <gamma0> <gammav>\"", (char *)NULL);
    return (TCL_ERROR);
  }
  /* check argument types */
  if ( !ARG_IS_D(2, temp) || !ARG_IS_D(3, gamma0) || !ARG_IS_D(4, gammav) ) {
    Tcl_AppendResult(interp, argv[0]," ",argv[1]," needs four DOUBLES", (char *)NULL);
    return (TCL_ERROR);
  }
  /* broadcast parameters */
  temperature = temp;
  nptiso_gamma0 = gamma0;
  nptiso_gammav = gammav;

  thermo_switch = ( thermo_switch | THERMO_NPT_ISO );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_parameter(FIELD_NPTISO_G0);
  mpi_bcast_parameter(FIELD_NPTISO_GV);
  return (TCL_OK);
}
#endif

#ifdef GHMC
int tclcommand_thermostat_parse_ghmc(Tcl_Interp *interp, int argc, char **argv) 
{
  double temp, phi;
	int n_md_steps;

  /* check number of arguments */
  if (argc < 5) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
		     argv[0]," ",argv[1]," <temp> <md_steps> <phi> {flip|no_flip|random_flip}\"", (char *)NULL);
    return (TCL_ERROR);
  }
  /* check argument types */
  if ( !ARG_IS_D(2, temp) || !ARG_IS_I(3, n_md_steps) || !ARG_IS_D(4, phi) ) {
    Tcl_AppendResult(interp, argv[2]," ",argv[3]," needs two DOUBLES and one INTEGER", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if (temp < 0 || n_md_steps <= 0) {
    Tcl_AppendResult(interp, "temperature and number of MD steps must be positive", (char *)NULL);
    return (TCL_ERROR);
  }

  if (phi < 0 || phi > 1) {
    Tcl_AppendResult(interp, "the parameter phi must be between zero and one", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if(argc > 5) {
    if (ARG_IS_S(5,"-flip")) {
      ghmc_mflip = GHMC_MFLIP_ON;
      
    } 
    else if (ARG_IS_S(5,"-no_flip")) {
      ghmc_mflip = GHMC_MFLIP_OFF;
    } 
    else if (ARG_IS_S(5,"-random_flip")) {
      ghmc_mflip = GHMC_MFLIP_RAND;
    } 
  } else {
    ghmc_mflip = GHMC_MFLIP_OFF;
  }
  
	if(argc > 6) {
    if (ARG_IS_S(6,"-scale")) {
      ghmc_tscale = GHMC_TSCALE_ON;
    } 
    else if (ARG_IS_S(6,"-no_scale")) {
      ghmc_tscale = GHMC_TSCALE_OFF;
    } 
  } else {
    ghmc_tscale = GHMC_TSCALE_OFF;
  }
  
  /* broadcast parameters */
  temperature = temp;
  ghmc_nmd = n_md_steps;
  ghmc_phi = phi;

  thermo_switch = ( thermo_switch | THERMO_GHMC );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_parameter(FIELD_GHMC_NMD);
  mpi_bcast_parameter(FIELD_GHMC_PHI);
  mpi_bcast_parameter(FIELD_GHMC_FLIP);
  mpi_bcast_parameter(FIELD_GHMC_SCALE);
  return (TCL_OK);
}
#endif

int tclcommand_thermostat_print_all(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  /* thermostat not initialized */
  if(temperature == -1.0) {
    Tcl_AppendResult(interp,"{ not initialized } ", (char *)NULL);
    return (TCL_OK);
  }

  /* no thermostat on */
  if(thermo_switch == THERMO_OFF) {
    Tcl_AppendResult(interp,"{ off } ", (char *)NULL);
    return (TCL_OK);
  }

  /* langevin */
  if(thermo_switch & THERMO_LANGEVIN ) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ langevin ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, langevin_gamma, buffer);
    Tcl_AppendResult(interp," ",buffer," } ", (char *)NULL);
  }
    
#ifdef DPD
 /* dpd */
  if(thermo_switch & THERMO_DPD) { tclcommand_thermostat_parse_and_print_dpd(interp);}
#endif

#ifdef NPT
  /* npt_isotropic */
  if(thermo_switch & THERMO_NPT_ISO) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ npt_isotropic ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nptiso_gamma0, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nptiso_gammav, buffer);
    Tcl_AppendResult(interp," ",buffer, " } ", (char *)NULL);
  }
#endif

#ifdef GHMC
  /* ghmc */
  if(thermo_switch & THERMO_GHMC) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ ghmc ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, ghmc_nmd, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, ghmc_phi, buffer);
    Tcl_AppendResult(interp," ",buffer, " ", (char *)NULL);
    switch (ghmc_mflip) {
      case GHMC_MFLIP_OFF:
        Tcl_AppendResult(interp, "-no_flip ", (char *)NULL);
				break;
      case GHMC_MFLIP_ON:
        Tcl_AppendResult(interp, "-flip ", (char *)NULL);
				break;
      case GHMC_MFLIP_RAND:
        Tcl_AppendResult(interp, "-random_flip ", (char *)NULL);
				break;
    }
    if (ghmc_tscale == GHMC_TSCALE_ON)
			Tcl_AppendResult(interp, "-scale }", (char *)NULL);
		else
			Tcl_AppendResult(interp, "-no_scale }", (char *)NULL);
	}
#endif

#if defined(LB) || defined(LB_GPU)
 /* lb */
  if(thermo_switch & THERMO_LB) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ lb ",buffer, " } ", (char *)NULL);
  }
#endif

#ifdef INTER_DPD
 /* inter_dpd */
  if(thermo_switch & THERMO_INTER_DPD) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ inter_dpd ",buffer, " } ", (char *)NULL);
  }
#endif
  return (TCL_OK);
}

int tclcommand_thermostat_print_usage(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "Usage of tcl-command thermostat:\n", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], "' for status return or \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " set off' to deactivate it (=> NVE-ensemble) \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " set langevin <temp> <gamma>' or \n ", (char *)NULL);
#ifdef DPD
  tclcommand_thermostat_print_usage_dpd(interp,argc,argv);
#endif
#ifdef NPT
  Tcl_AppendResult(interp, "'", argv[0], " set npt_isotropic <temp> <gamma0> <gammav>' ", (char *)NULL);
#endif
#ifdef LB
  Tcl_AppendResult(interp, "'", argv[0], " set lb <temperature>" , (char *)NULL);
#endif
#ifdef LB_GPU
  Tcl_AppendResult(interp, "'", argv[0], " set lb_gpu <temperature>" , (char *)NULL);
#endif
  return (TCL_ERROR);
}

int tclcommand_thermostat(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int err = TCL_OK;
  THERMO_TRACE(fprintf(stderr,"%d: thermostat:\n",this_node));

  /* print thermostat status */
  if(argc == 1) return tclcommand_thermostat_print_all(interp);
  
  if ( ARG1_IS_S("set") )          {
    argc--;
    argv++;

    if (argc == 1) {
      Tcl_AppendResult(interp, "wrong # args: \n", (char *)NULL);
      return tclcommand_thermostat_print_usage(interp, argc, argv);
    }
  }
  if ( ARG1_IS_S("off") )
    err = tclcommand_thermostat_parse_off(interp, argc, argv);
  else if ( ARG1_IS_S("langevin"))
    err = tclcommand_thermostat_parse_langevin(interp, argc, argv);
#ifdef DPD
  else if ( ARG1_IS_S("dpd") )
    err = tclcommand_thermostat_parse_dpd(interp, argc, argv);
#endif
#ifdef INTER_DPD
  else if ( ARG1_IS_S("inter_dpd") )
    err = tclcommand_thermostat_parse_inter_dpd(interp, argc, argv);
#endif
#ifdef NPT
  else if ( ARG1_IS_S("npt_isotropic") )
    err = tclcommand_thermostat_parse_npt_isotropic(interp, argc, argv);
#endif
#if defined(LB) || defined(LB_GPU)
  else if ( ARG1_IS_S("lb"))
    err = tclcommand_thermostat_parse_lb(interp, argc-1, argv+1);
#endif
#ifdef GHMC
  else if ( ARG1_IS_S("ghmc") )
    err = tclcommand_thermostat_parse_ghmc(interp, argc, argv);
#endif
  else {
    Tcl_AppendResult(interp, "Unknown thermostat ", argv[1], "\n", (char *)NULL);
    return tclcommand_thermostat_print_usage(interp, argc, argv);
  }
  return gather_runtime_errors(interp, err);
}

int tclcommand_thermostat_parse_lb(Tcl_Interp *interp, int argc, char ** argv)
{

#if defined(LB) || defined(LB_GPU)
  double temp;

  /* get lb interaction type */
  if (argc < 2) {
    Tcl_AppendResult(interp, "lattice-Boltzmann needs 1 parameter: "
		     "<temperature>",
		     (char *) NULL);
    return TCL_ERROR;
  }

  /* copy lattice-boltzmann parameters */
  if (! ARG_IS_D(1, temp)) { return TCL_ERROR; }

  if ( temp < 0.0 ) {
    Tcl_AppendResult(interp, "temperature must be non-negative", (char *) NULL);
    return TCL_ERROR;
  }
  temperature = temp;
  thermo_switch = ( thermo_switch | THERMO_LB );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  
#endif
  return TCL_OK;
}
