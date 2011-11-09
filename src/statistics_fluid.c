/*
  Copyright (C) 2010,2011 The ESPResSo project
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
#include "lb-boundaries.h"
#include "statistics_fluid.h"
#include "lbgpu.h"

#ifdef LB

//#include <fftw3.h>

/** Caclulate mass of the LB fluid.
 * \param result Fluid mass
 */
void lb_calc_fluid_mass(double *result) {
  int x, y, z, index;
  double sum_rho=0.0, rho=0.0;

  for (x=1; x<=lblattice.grid[0]; x++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
      for (z=1; z<=lblattice.grid[2]; z++) {
	index = get_linear_index(x,y,z,lblattice.halo_grid);

	lb_calc_local_rho(index,&rho);
	//fprintf(stderr,"(%d,%d,%d) %e\n",x,y,z,rho);
	sum_rho += rho;

      }
    }
  }

  MPI_Reduce(&sum_rho, result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

/** Calculate momentum of the LB fluid.
 * \param result Fluid momentum
 */
void lb_calc_fluid_momentum(double *result) {

    int x, y, z, index;
    double j[3], momentum[3] = { 0.0, 0.0, 0.0 };

    for (x=1; x<=lblattice.grid[0]; x++) {
	for (y=1; y<=lblattice.grid[1]; y++) {
	    for (z=1; z<=lblattice.grid[2]; z++) {
		index = get_linear_index(x,y,z,lblattice.halo_grid);

		lb_calc_local_j(index,j);
		momentum[0] += j[0] + lbfields[index].force[0];
		momentum[1] += j[1] + lbfields[index].force[1];
		momentum[2] += j[2] + lbfields[index].force[2];

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
  double rho, j[3];
  double temp = 0.0;

  for (x=1; x<=lblattice.grid[0]; x++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
      for (z=1; z<=lblattice.grid[2]; z++) {
	index = get_linear_index(x,y,z,lblattice.halo_grid);
	
	lb_calc_local_fields(index, &rho, j, NULL);

	temp += scalar(j,j);
      }
    }
  }

  temp *= 1./(3.*lbpar.rho*lblattice.grid_volume*lbpar.tau*lbpar.tau*lblattice.agrid)/n_nodes;

  MPI_Reduce(&temp, result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void lb_collect_boundary_forces(double *result) {
#ifdef LB_BOUNDARIES
  double* boundary_forces=malloc(3*n_lb_boundaries*sizeof(double));

  for (int i = 0; i < n_lb_boundaries; i++) 
    for (int j = 0; j < 3; j++)
      boundary_forces[3*i+j]=lb_boundaries[i].force[j];

  MPI_Reduce(boundary_forces, result, 3*n_lb_boundaries, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  free(boundary_forces);
#endif
}
#define REQ_DENSPROF 702

/** Calculate a density profile of the fluid. */
void lb_calc_densprof(double *result, int *params) {

  int index, dir[3], grid[3];
  int newroot=0, subrank, involved=0;
  double *profile;
  MPI_Comm slice_comm;
  MPI_Status status;

  if (this_node !=0) params = malloc(3*sizeof(int));

  MPI_Bcast(params, 3, MPI_INT, 0, MPI_COMM_WORLD);

  int pdir = params[0];
  int x1   = params[1];
  int x2   = params[2];

  dir[pdir] = 0;
  dir[(pdir+1)%3] = x1;
  dir[(pdir+2)%3] = x2;

  newroot = map_lattice_to_node(&lblattice, dir, grid);
  map_node_array(this_node, node_pos);

  if (     (grid[(pdir+1)%3] == node_pos[(pdir+1)%3])
	&& (grid[(pdir+2)%3] == node_pos[(pdir+2)%3]) ) {
    involved = 1;
  }

  MPI_Comm_split(MPI_COMM_WORLD, involved, this_node, &slice_comm);
  MPI_Comm_rank(slice_comm, &subrank);

  if (this_node == newroot)
    result = realloc(result,box_l[pdir]/lblattice.agrid*sizeof(double));

  if (involved) {

    profile = malloc(lblattice.grid[pdir]*sizeof(double));
      
    //dir[(pdir+1)%3] += 1;
    //dir[(pdir+2)%3] += 1;
    for (dir[pdir]=1;dir[pdir]<=lblattice.grid[pdir];dir[pdir]++) {

      index = get_linear_index(dir[0],dir[1],dir[2],lblattice.halo_grid);
      lb_calc_local_rho(index,&profile[dir[pdir]-1]);
      //profile[dir[pdir]-1] = *lbfluid[index].rho;

      //if (dir[pdir]==lblattice.grid[pdir]) {
      //	int i;
      //	fprintf(stderr,"(%d,%d,%d)\n",dir[0],dir[1],dir[2]);
      //	fprintf(stderr,"%d\n",lbfluid[index].boundary);
      //	for (i=0;i<lbmodel.n_veloc;i++) fprintf(stderr,"local_n[%p][%d]=%.12e\n",lbfluid[index].n,i,lbfluid[index].n[i]+lbmodel.coeff[i][0]*lbpar.rho);
      //	fprintf(stderr,"local_rho=%e\n",*lbfluid[index].rho);
      //}

  }

    MPI_Gather(profile, lblattice.grid[pdir], MPI_DOUBLE, result, lblattice.grid[pdir], MPI_DOUBLE, 0, slice_comm);
    
    free(profile);

  }

  MPI_Comm_free(&slice_comm);

  if (newroot != 0) {
    if (this_node == newroot) {
      MPI_Send(result, lblattice.grid[pdir]*node_grid[pdir], MPI_DOUBLE, 0, REQ_DENSPROF, MPI_COMM_WORLD);
      free(result);
    }
    if (this_node == 0) {
      MPI_Recv(result, lblattice.grid[pdir]*node_grid[pdir], MPI_DOUBLE, newroot, REQ_DENSPROF, MPI_COMM_WORLD, &status);
    }
  }

  if (this_node != 0) free(params);

}

static void lb_master_calc_densprof(double *profile, int pdir, int x1, int x2) {

  /* this is a quick and dirty hack to issue slave calls with parameters */

  int params[3] = { pdir, x1, x2 };
  
  mpi_gather_stats(9, profile, params, NULL, NULL);
  
}

#define REQ_VELPROF 701

/** Calculate a velocity profile for the LB fluid. */	
void lb_calc_velprof(double *result, int *params) {

  int index, dir[3], grid[3];
  int newroot=0, subrank, involved=0;
  double rho, j[3], *velprof;
  MPI_Comm slice_comm;
  MPI_Status status;

  if (this_node != 0) params = malloc(4*sizeof(int));

  MPI_Bcast(params, 4, MPI_INT, 0, MPI_COMM_WORLD);

  int vcomp = params[0];
  int pdir  = params[1];
  int x1    = params[2];
  int x2    = params[3];

  dir[pdir] = 0;
  dir[(pdir+1)%3] = x1;
  dir[(pdir+2)%3] = x2;

  //fprintf(stderr,"%d: (%d,%d,%d)\n",this_node,dir[0],dir[1],dir[2]);

  newroot = map_lattice_to_node(&lblattice, dir, grid);
  map_node_array(this_node, node_pos);

  //fprintf(stderr,"%d: newroot=%d (%d,%d,%d)\n",this_node,newroot,grid[0],grid[1],grid[2]);

  if (    (grid[(pdir+1)%3] == node_pos[(pdir+1)%3])
       && (grid[(pdir+2)%3] == node_pos[(pdir+2)%3]) ) {
    involved = 1;
  }

  MPI_Comm_split(MPI_COMM_WORLD, involved, this_node, &slice_comm);
  MPI_Comm_rank(slice_comm, &subrank);

  if (this_node == newroot) 
    result = realloc(result,box_l[pdir]/lblattice.agrid*sizeof(double));

  //fprintf(stderr,"%d (%d,%d): result=%p vcomp=%d pdir=%d x1=%d x2=%d involved=%d\n",this_node,subrank,newroot,result,vcomp,pdir,x1,x2,involved);

  if (involved) {

    velprof = malloc(lblattice.grid[pdir]*sizeof(double));

    //dir[(pdir+1)%3] += 1;
    //dir[(pdir+2)%3] += 1;
    for (dir[pdir]=1;dir[pdir]<=lblattice.grid[pdir];dir[pdir]++) {
      
      index = get_linear_index(dir[0],dir[1],dir[2],lblattice.halo_grid);
      lb_calc_local_fields(index, &rho, j, NULL);
      
      //fprintf(stderr,"%p %d %.12e %.12e %d\n",lbfluid[0],index,rho,j[0],vcomp);

      if (rho < ROUND_ERROR_PREC) {
	velprof[dir[pdir]-1] = 0.0;
      } else {
	//velprof[dir[pdir]-1] = local_j / (SQR(lbpar.agrid)*lbpar.tau);
	velprof[dir[pdir]-1] = j[vcomp]/rho * lblattice.agrid/lbpar.tau;
	//fprintf(stderr,"%f %f %f\n",velprof[dir[pdir]-1],local_j,local_rho);
      }

      //if (dir[pdir]==lblattice.grid[pdir]) {
      //	int i;
      //	fprintf(stderr,"(%d,%d,%d)\n",dir[0],dir[1],dir[2]);
      //	fprintf(stderr,"%d\n",lbfluid[index].boundary);
      //	for (i=0;i<lbmodel.n_veloc;i++) fprintf(stderr,"local_n[%p][%d]=%.12e\n",lbfluid[index].n,i,lbfluid[index].n[i]+lbmodel.coeff[i][0]*lbpar.rho);
      //	fprintf(stderr,"local_rho=%e\n",local_rho);
      //	fprintf(stderr,"local_j=%e\n",local_j);
      //}

    }
    
    MPI_Gather(velprof, lblattice.grid[pdir], MPI_DOUBLE, result, lblattice.grid[pdir], MPI_DOUBLE, newroot, slice_comm);

    free(velprof);

  } 

  MPI_Comm_free(&slice_comm);

  if (newroot != 0) {
    if (this_node == newroot) {
      MPI_Send(result, lblattice.grid[pdir]*node_grid[pdir], MPI_DOUBLE, 0, REQ_VELPROF, MPI_COMM_WORLD);
      free(result);
    }
    if (this_node == 0) {
      //fprintf(stderr,"%d: I'm just here!\n",this_node);
      MPI_Recv(result, lblattice.grid[pdir]*node_grid[pdir], MPI_DOUBLE, newroot, REQ_VELPROF, MPI_COMM_WORLD, &status);
      //fprintf(stderr,"%d: And now I'm here!\n",this_node);
    }
  }

  if (this_node !=0) free(params);

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

  profile = malloc(lblattice.grid[pdir]*node_grid[pdir]*sizeof(double));

  lb_master_calc_densprof(profile, pdir, x1, x2);

  for (i=0; i<lblattice.grid[pdir]*node_grid[pdir]; i++) {
    Tcl_PrintDouble(interp, i*lblattice.agrid, buffer);
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

    velprof = malloc(box_l[pdir]/lblattice.agrid*sizeof(double));

    lb_master_calc_velprof(velprof, vcomp, pdir, x1, x2);

    for (i=0; i<box_l[pdir]/lblattice.agrid; i++) {
	Tcl_PrintDouble(interp, i*lblattice.agrid, buffer);
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
		fprintf(stderr, "sry not implemented yet");
      //err = parse_analyze_fluid_mass(interp, argc - 1, argv + 1);
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
