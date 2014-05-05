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
/** \file statistics_fluid.cpp
 *
 * Fluid related analysis functions.
 * Implementation of \ref statistics_fluid.hpp.
 *
 */

#include <mpi.h>
#include <tcl.h>
//#include "utils.hpp"
#include "parser.hpp"
#include "communication.hpp"
#include "lb.hpp"
#include "lb-boundaries.hpp"
#include "statistics_fluid.hpp"

#ifdef LB

static void lb_master_calc_densprof(double *profile, int pdir, int x1, int x2) {
    
    /* this is a quick and dirty hack to issue slave calls with parameters */
    
    int params[3] = { pdir, x1, x2 };
    
    mpi_gather_stats(9, profile, params, NULL, NULL);
    
}

static void lb_master_calc_velprof(double *result, int vcomp, int pdir, int x1, int x2) {
    
    /* this is a quick and dirty hack to issue slave calls with parameters */
    
    int params[4];
    
    params[0] = vcomp;
    params[1] = pdir;
    params[2] = x1;
    params[3] = x2;
    
    mpi_gather_stats(8, result, params, NULL, NULL);
    
}


static int tclcommand_analyze_fluid_parse_mass(Tcl_Interp *interp, int argc, char** argv) {
  char buffer[TCL_DOUBLE_SPACE];
  double mass;

  mpi_gather_stats(5, &mass, NULL, NULL, NULL);

  Tcl_PrintDouble(interp, mass, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  return TCL_OK;
}

static int tclcommand_analyze_fluid_parse_momentum(Tcl_Interp* interp, int argc, char *argv[]) {
  char buffer[TCL_DOUBLE_SPACE];
  double mom[3];

  mpi_gather_stats(6, mom, NULL, NULL, NULL);
  
  Tcl_PrintDouble(interp, mom[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, mom[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, mom[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  return TCL_OK;
}

static int tclcommand_analyze_fluid_parse_temp(Tcl_Interp *interp, int argc, char *argv[]) {
  char buffer[TCL_DOUBLE_SPACE];
  double temp;

  mpi_gather_stats(7, &temp, NULL, NULL, NULL);

  Tcl_PrintDouble(interp, temp, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  return TCL_OK;
}

static int tclcommand_analyze_fluid_parse_densprof(Tcl_Interp *interp, int argc, char **argv) {
  
  int i, pdir, x1, x2;
  char buffer[TCL_DOUBLE_SPACE];
  double *profile;

  if (argc <3) {
    Tcl_AppendResult(interp, "usage: analyze fluid density <p_dir> <x1> <x2>", (char *)NULL);
    return TCL_ERROR;
  }

  if (!ARG_IS_I(0,pdir)) return TCL_ERROR;
  if (!ARG_IS_I(1,x1)) return TCL_ERROR;
  if (!ARG_IS_I(2,x2)) return TCL_ERROR;

  if (pdir != 2) {
    Tcl_AppendResult(interp, "analyze fluid density is only implemented for pdir=2 yet!", (char *)NULL);
    return TCL_ERROR;
  }

  profile = (double*) malloc(lblattice.grid[pdir]*node_grid[pdir]*sizeof(double));

  lb_master_calc_densprof(profile, pdir, x1, x2);

  for (i=0; i<lblattice.grid[pdir]*node_grid[pdir]; i++) {
    Tcl_PrintDouble(interp, i*lbpar.agrid, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, profile[i], buffer);
    Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);
  }
  free(profile);
  
  return TCL_OK;

}

static int tclcommand_analyze_fluid_parse_velprof(Tcl_Interp *interp, int argc, char **argv) {
    int i, pdir, vcomp, x1, x2;
    char buffer[TCL_DOUBLE_SPACE];
    double *velprof;

    //fprintf(stderr, "NOTE: analyze fluid velprof is not completely implemented by now.\n      The calling interface might still change without backwards compatibility!\n");

    /*
    if (n_nodes > 1) {
	Tcl_AppendResult(interp, "velocity profile not yet implemented for parallel execution!", (char *)NULL);
	return TCL_ERROR;
    }
    */

    if (argc < 4) {
	Tcl_AppendResult(interp, "usage: analyze fluid velprof <v_comp> <p_dir> <x1> <x2>", (char *)NULL);
	return TCL_ERROR;
    }

    if (!ARG_IS_I(0,vcomp)) return TCL_ERROR;
    if (!ARG_IS_I(1,pdir)) return TCL_ERROR;
    if (!ARG_IS_I(2,x1)) return TCL_ERROR;
    if (!ARG_IS_I(3,x2)) return TCL_ERROR;

    if (pdir != 2) {
	Tcl_AppendResult(interp, "analyze fluid velprof is only implemented for pdir=2 yet!", (char *)NULL);
	return TCL_ERROR;
    }
    if (vcomp != 0) {
	Tcl_AppendResult(interp, "analyze fluid velprof is only implemented for vdir=0 yet", (char *)NULL);
	return TCL_ERROR;
    }

    velprof = (double*) malloc(box_l[pdir]/lblattice.agrid[pdir]*sizeof(double));

    lb_master_calc_velprof(velprof, vcomp, pdir, x1, x2);

    for (i=0; i<box_l[pdir]/lbpar.agrid; i++) {
	Tcl_PrintDouble(interp, i*lbpar.agrid, buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	Tcl_PrintDouble(interp, velprof[i], buffer);
	Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);
    }	

    free(velprof);

    return TCL_OK;

}
#endif /* LB */

/** Parser for fluid related analysis functions. */
int tclcommand_analyze_parse_fluid_cpu(Tcl_Interp *interp, int argc, char **argv) {
#ifdef LB
    int err = TCL_ERROR;

    if (argc==0) {
	Tcl_AppendResult(interp, "usage: analyze fluid <what>", (char *)NULL);
	return TCL_ERROR;
    } 

    if (ARG0_IS_S("mass"))
      err = tclcommand_analyze_fluid_parse_mass(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("momentum"))
      err = tclcommand_analyze_fluid_parse_momentum(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("temperature"))
      err = tclcommand_analyze_fluid_parse_temp(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("density"))
      err = tclcommand_analyze_fluid_parse_densprof(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("velprof"))
      err = tclcommand_analyze_fluid_parse_velprof(interp, argc - 1, argv + 1);
    else {
	Tcl_AppendResult(interp, "unkown feature \"", argv[0], "\" of analyze fluid", (char *)NULL);
	return TCL_ERROR;
    }

    return err;
#else /* !defined LB */
  Tcl_AppendResult(interp, "LB is not compiled in!?", NULL);
  return TCL_ERROR;
#endif
}

#ifdef LB_GPU
static int tclcommand_analyze_fluid_mass_gpu(Tcl_Interp* interp, int argc, char *argv[]){
  char buffer[TCL_DOUBLE_SPACE];
  double host_mass[1];
  
  lb_calc_fluid_mass_GPU(host_mass);

  Tcl_PrintDouble(interp, host_mass[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);

  return TCL_OK;

}

static int tclcommand_analyze_fluid_parse_momentum_gpu(Tcl_Interp* interp, int argc, char *argv[]) {
  char buffer[TCL_DOUBLE_SPACE];
  double host_mom[3];

  lb_calc_fluid_momentum_GPU(host_mom);
  
  Tcl_PrintDouble(interp, host_mom[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, host_mom[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, host_mom[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  return TCL_OK;
}

static int tclcommand_analyze_fluid_parse_temperature_gpu(Tcl_Interp* interp, int argc, char *argv[]) {
  char buffer[TCL_DOUBLE_SPACE];
  double host_temp[1];

  lb_calc_fluid_temperature_GPU(host_temp);
  
  Tcl_PrintDouble(interp, host_temp[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);

  return TCL_OK;
}
#endif

int tclcommand_analyze_parse_fluid_gpu(Tcl_Interp *interp, int argc, char **argv) {
#ifdef LB_GPU
    int err = TCL_ERROR;

    if (argc==0) {
	Tcl_AppendResult(interp, "usage: analyze fluid <what>", (char *)NULL);
	return TCL_ERROR;
    } 

    if (ARG0_IS_S("mass"))
		//fprintf(stderr, "sry not implemented yet");
      err = tclcommand_analyze_fluid_mass_gpu(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("momentum"))
      err = tclcommand_analyze_fluid_parse_momentum_gpu(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("temperature"))
      err = tclcommand_analyze_fluid_parse_temperature_gpu(interp, argc - 1, argv + 1);

    else if (ARG0_IS_S("velprof"))
		fprintf(stderr, "sry not implemented yet");
      //err = parse_analyze_fluid_velprof(interp, argc - 1, argv + 1);
    else {
	Tcl_AppendResult(interp, "unkown feature \"", argv[0], "\" of analyze fluid", (char *)NULL);
	return TCL_ERROR;
    }

    return err;
#else /* !defined LB_GPU */
  Tcl_AppendResult(interp, "LB_GPU is not compiled in!", NULL);
  return TCL_ERROR;
#endif
}

int tclcommand_analyze_parse_fluid(Tcl_Interp *interp, int argc, char **argv) {

  if (lattice_switch & LATTICE_LB_GPU){
      return tclcommand_analyze_parse_fluid_gpu(interp, argc, argv);
  } else 
      return tclcommand_analyze_parse_fluid_cpu(interp, argc, argv);

}
