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
#include "utils.hpp"
#include "communication.hpp"
#include "lb.hpp"
#include "lb-boundaries.hpp"
#include "statistics_fluid.hpp"

#ifdef LB

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

  MPI_Reduce(&sum_rho, result, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
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

    momentum[0] *= lbpar.agrid/lbpar.tau;
    momentum[1] *= lbpar.agrid/lbpar.tau;
    momentum[2] *= lbpar.agrid/lbpar.tau;

    MPI_Reduce(momentum, result, 3, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
    
}

/** Calculate temperature of the LB fluid.
 * \param result Fluid temperature
 */
void lb_calc_fluid_temp(double *result) {
  int x, y, z, index;
  double rho, j[3];
  double temp = 0.0;
  int number_of_non_boundary_nodes = 0;

  for (x=1; x<=lblattice.grid[0]; x++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
      for (z=1; z<=lblattice.grid[2]; z++) {

	index = get_linear_index(x,y,z,lblattice.halo_grid);

#ifdef LB_BOUNDARIES
	if ( !lbfields[index].boundary )
#endif
	  {
	    lb_calc_local_fields(index, &rho, j, NULL);
	    temp += scalar(j,j);
	    number_of_non_boundary_nodes++;
	  }
      }
    }
  }

  // @Todo: lblattice.agrid is 3d. What to use here?
  temp *= 1./(3.*lbpar.rho[0]*number_of_non_boundary_nodes*
              lbpar.tau*lbpar.tau*lbpar.agrid)/n_nodes;

  MPI_Reduce(&temp, result, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}

void lb_collect_boundary_forces(double *result) {
#ifdef LB_BOUNDARIES
  double* boundary_forces = (double*) malloc(3*n_lb_boundaries*sizeof(double));

  for (int i = 0; i < n_lb_boundaries; i++) 
    for (int j = 0; j < 3; j++)
      boundary_forces[3*i+j]=lb_boundaries[i].force[j];

  MPI_Reduce(boundary_forces, result, 3*n_lb_boundaries, 
             MPI_DOUBLE, MPI_SUM, 0, comm_cart);
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

  if (this_node !=0) params = (int*) malloc(3*sizeof(int));

  MPI_Bcast(params, 3, MPI_INT, 0, comm_cart);

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

  MPI_Comm_split(comm_cart, involved, this_node, &slice_comm);
  MPI_Comm_rank(slice_comm, &subrank);

  if (this_node == newroot)
    result = (double*) realloc(result,box_l[pdir]/lblattice.agrid[pdir]*sizeof(double));

  if (involved) {

    profile = (double*) malloc(lblattice.grid[pdir]*sizeof(double));
      
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
      MPI_Send(result, lblattice.grid[pdir]*node_grid[pdir], MPI_DOUBLE, 0, REQ_DENSPROF, comm_cart);
      free(result);
    }
    if (this_node == 0) {
      MPI_Recv(result, lblattice.grid[pdir]*node_grid[pdir], MPI_DOUBLE, newroot, REQ_DENSPROF, comm_cart, &status);
    }
  }

  if (this_node != 0) free(params);

}

#define REQ_VELPROF 701

/** Calculate a velocity profile for the LB fluid. */	
void lb_calc_velprof(double *result, int *params) {

  int index, dir[3], grid[3];
  int newroot=0, subrank, involved=0;
  double rho, j[3], *velprof;
  MPI_Comm slice_comm;
  MPI_Status status;

  if (this_node != 0) params = (int*) malloc(4*sizeof(int));

  MPI_Bcast(params, 4, MPI_INT, 0, comm_cart);

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

  MPI_Comm_split(comm_cart, involved, this_node, &slice_comm);
  MPI_Comm_rank(slice_comm, &subrank);

  if (this_node == newroot) 
    result = (double*) realloc(result,box_l[pdir]/lblattice.agrid[pdir]*sizeof(double));

  //fprintf(stderr,"%d (%d,%d): result=%p vcomp=%d pdir=%d x1=%d x2=%d involved=%d\n",this_node,subrank,newroot,result,vcomp,pdir,x1,x2,involved);

  if (involved) {

    velprof = (double*) malloc(lblattice.grid[pdir]*sizeof(double));

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
	velprof[dir[pdir]-1] = j[vcomp]/rho * lbpar.agrid/lbpar.tau;
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
      MPI_Send(result, lblattice.grid[pdir]*node_grid[pdir], MPI_DOUBLE, 0, REQ_VELPROF, comm_cart);
      free(result);
    }
    if (this_node == 0) {
      //fprintf(stderr,"%d: I'm just here!\n",this_node);
      MPI_Recv(result, lblattice.grid[pdir]*node_grid[pdir], MPI_DOUBLE, newroot, REQ_VELPROF, comm_cart, &status);
      //fprintf(stderr,"%d: And now I'm here!\n",this_node);
    }
  }

  if (this_node !=0) free(params);

}

#endif /* LB */

