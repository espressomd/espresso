/*
  Copyright (C) 2010,2011 The ESPResSo project
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
/** \file statistics_fluid.c
 *
 * Fluid related analysis functions.
 * Implementation of \ref statistics_fluid.h.
 *
 */

#include <mpi.h>
#include "utils.h"
#include "parser.h"
#include "communication.h"
#include "lb.h"
#include "statistics_fluid.h"
#include "lbgpu.h"

#ifdef LB

#include <fftw3.h>

/** Caclulate mass of the LB fluid.
 * \param result Fluid mass
 */
void lb_calc_fluid_mass(double *result) {
  int x, y, z, index;
  double mass = 0.0;

  for (x=1; x<=lblattice.grid[0]; x++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
      for (z=1; z<=lblattice.grid[2]; z++) {
	index = get_linear_index(x,y,z,lblattice.halo_grid);

	lb_calc_local_rho(&lbfluid[index]);
	mass += *lbfluid[index].rho;

      }
    }
  }

  MPI_Reduce(&mass, result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

/** Calculate momentum of the LB fluid.
 * \param result Fluid momentum
 */
void lb_calc_fluid_momentum(double *result) {

    int x, y, z, index;
    double momentum[3] = { 0.0, 0.0, 0.0 };

    for (x=1; x<=lblattice.grid[0]; x++) {
	for (y=1; y<=lblattice.grid[1]; y++) {
	    for (z=1; z<=lblattice.grid[2]; z++) {
		index = get_linear_index(x,y,z,lblattice.halo_grid);

		lb_calc_local_j(&lbfluid[index]);
		momentum[0] += lbfluid[index].j[0];
		momentum[1] += lbfluid[index].j[1];
		momentum[2] += lbfluid[index].j[2];

	    }
	}
    }

    momentum[0] *= lblattice.agrid/lbpar.tau;
    momentum[1] *= lblattice.agrid/lbpar.tau;
    momentum[2] *= lblattice.agrid/lbpar.tau;

    MPI_Reduce(momentum, result, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
}

/** Calculate temperature of the LB fluid.
 * \param result Fluid temperature
 */
void lb_calc_fluid_temp(double *result) {
  int x, y, z, index;
  double local_rho, local_j2;
  double temp = 0.0;

  for (x=1; x<=lblattice.grid[0]; x++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
      for (z=1; z<=lblattice.grid[2]; z++) {
	index = get_linear_index(x,y,z,lblattice.halo_grid);
	
	lb_calc_local_j(&lbfluid[index]);
	lb_calc_local_rho(&lbfluid[index]);

	local_rho = *lbfluid[index].rho;
	local_j2  = scalar(lbfluid[index].j,lbfluid[index].j);

	temp += local_j2;
      }
    }
  }

  temp *= 1./(lbpar.rho*lblattice.grid_volume*lbpar.tau*lbpar.tau*pow(lblattice.agrid,4));

  MPI_Reduce(&temp, result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

/** Calculate a velocity profile for the LB fluid. */	
void lb_calc_velocity_profile(double *velprof, int vcomp, int pdir, int x1, int x2) {

  int index, dir[3];
  double local_rho, local_j;
  /* \todo generalize and parallelize */

  dir[(pdir+1)%3] = x1;
  dir[(pdir+2)%3] = x2;
  for (dir[pdir]=1;dir[pdir]<=lblattice.grid[pdir];dir[pdir]++) {

      index = get_linear_index(dir[0],dir[1],dir[2],lblattice.halo_grid);
      lb_calc_local_j(&lbfluid[index]);
      lb_calc_local_rho(&lbfluid[index]);
      local_rho = *lbfluid[index].rho;
      local_j = lbfluid[index].j[vcomp]; 
      if (local_j == 0) {
	velprof[dir[pdir]-1] = 0.0;
      } else {
	velprof[dir[pdir]-1] = local_j/local_rho * lblattice.agrid/lbpar.tau;
      }

  }

}
/* TODO: This function is not used anywhere. To be removed?  */
/** Fourier transform the stress tensor into k-space using FFTW */
static void lb_calc_fourier_pi() {

  static fftw_plan plan;
  static int initialized = 0;
  static double *data;
  //static fftw_complex *result;

  /* prepare plan for FFTW */
  if (!initialized) {
    fftw_destroy_plan(plan);
    fftw_free(data);
    //data = fftw_malloc(volume*nfields*sizeof(double));
    //result = fftw_malloc((volume/2+1)*nfields*sizeof(fftw_complex));
    //plan = fftw_plan_many_dft_r2c(3,grid,6,data,NULL,nfields,1,result,NULL,nfields,1,FFTW_PATIENT);
  }

  /* prepare data: collect lattice from all nodes */

  /* compute 3d real-to-complex Fourier transform */
  fftw_execute(plan);

  /* */

}

/***********************************************************************/

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

static int tclcommand_analyze_fluid_parse_velprof(Tcl_Interp *interp, int argc, char **argv) {
    int i, pdir, vcomp, x1, x2;
    char buffer[TCL_DOUBLE_SPACE];
    double *velprof;

    fprintf(stderr, "NOTE: analyze fluid velprof is not completely implemented by now.\n      The calling interface might still change without backwards compatibility!\n");

    if (n_nodes > 1) {
	Tcl_AppendResult(interp, "velocity profil not yet implemented for parallel execution!", (char *)NULL);
	return TCL_ERROR;
    }

    if (argc < 4) {
	Tcl_AppendResult(interp, "usage: analyze fluid velprof <v_comp> <p_dir> <x1> <x2> <v_comp>", (char *)NULL);
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

    velprof = malloc(lblattice.grid[pdir]*sizeof(double));
    lb_calc_velocity_profile(velprof, vcomp, pdir, x1, x2);    
    for (i=0; i<lblattice.grid[pdir]; i++) {
	Tcl_PrintDouble(interp, my_left[pdir]+i*lblattice.agrid, buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	Tcl_PrintDouble(interp, velprof[i], buffer);
	Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);
    }	
    free(velprof);

    return TCL_OK;

}

/** Parser for fluid related analysis functions. */
int tclcommand_analyze_parse_fluid(Tcl_Interp *interp, int argc, char **argv) {
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
    else if (ARG0_IS_S("velprof"))
      err = tclcommand_analyze_fluid_parse_velprof(interp, argc - 1, argv + 1);
    else {
	Tcl_AppendResult(interp, "unkown feature \"", argv[0], "\" of analyze fluid", (char *)NULL);
	return TCL_ERROR;
    }

    return err;
}

#endif /* LB */

#ifdef LB_GPU

int tclcommand_analyze_parse_fluid(Tcl_Interp *interp, int argc, char **argv) {
    int err = TCL_ERROR;

    if (argc==0) {
	Tcl_AppendResult(interp, "usage: analyze fluid gpu <what>", (char *)NULL);
	return TCL_ERROR;
    } 

    if (ARG0_IS_S("mass"))
		fprintf(stderr, "sry not implemented yet");
      //err = parse_analyze_fluid_mass(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("momentum"))
		fprintf(stderr, "sry not implemented yet");
      //err = parse_analyze_fluid_momentum(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("temperature"))
		fprintf(stderr, "sry not implemented yet");
      //err = parse_analyze_fluid_temp(interp, argc - 1, argv + 1);
    else if (ARG0_IS_S("velprof"))
		fprintf(stderr, "sry not implemented yet");
      //err = parse_analyze_fluid_velprof(interp, argc - 1, argv + 1);
    else {
	Tcl_AppendResult(interp, "unkown feature \"", argv[0], "\" of analyze fluid", (char *)NULL);
	return TCL_ERROR;
    }

    return err;
}

#endif /* LB_GPU */ 
