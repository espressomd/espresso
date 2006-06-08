/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you
 * accepted upon receiving the distribution and by which you are
 * legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * You should have received a copy of that license along with this
 * program; if not, refer to http://www.espresso.mpg.de/license.html
 * where its current version can be found, or write to
 * Max-Planck-Institute for Polymer Research, Theory Group, 
 * PO Box 3148, 55021 Mainz, Germany. 
 * Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
 */

/** \file lb.c
 *
 * Lattice Boltzmann algorithm for hydrodynamic degrees of freedom.
 *
 * Includes fluctuating LB and coupling to MD particles via frictional 
 * momentum transfer.
 *
 */

#include <mpi.h>
#include <tcl.h>
#include <stdio.h>
#include <fftw3.h>
#include "utils.h"
#include "parser.h"
#include "communication.h"
#include "grid.h"
#include "domain_decomposition.h"
#include "interaction_data.h"
#include "thermostat.h"
#include "lattice.h"
#include "halo.h"
#include "lb-boundaries.h"
#include "lb.h"

#ifdef LB

/** Flag indicating momentum exchange between particles and fluid */
int transfer_momentum = 0;

/** Struct holding the Lattice Boltzmann parameters */
LB_Parameters lbpar = { -1.0, -1.0, -1.0, -1.0, 0.0, { 0.0, 0.0, 0.0} };

/** The underlying lattice structure */
Lattice lblattice = { {0,0,0}, {0,0,0}, 0, 0, 0, 0, -1.0, -1.0, NULL, NULL };

/** Pointer to the fluid nodes
 * This variable is used for convenience instead of having to type lattice.fields everywhere */
LB_FluidNode *lbfluid=NULL;

/** Communicator for halo exchange between processors */
static HaloCommunicator update_halo_comm = { 0, NULL };

/** Velocity sub-lattice of the D3Q18 model */
static double d3q18_lattice[18][3] = { {  1.,  0.,  0. },
          			       { -1.,  0.,  0. },
          		               {  0.,  1.,  0. }, 
          			       {  0., -1.,  0. },
          			       {  0.,  0.,  1. }, 
          			       {  0.,  0., -1. },
          			       {  1.,  1.,  0. }, 
          			       { -1., -1.,  0. },
          			       {  1., -1.,  0. },
          			       { -1.,  1.,  0. },
          			       {  1.,  0.,  1. },
          			       { -1.,  0., -1. },
          			       {  1.,  0., -1. },
          			       { -1.,  0.,  1. },
          			       {  0.,  1.,  1. },
          			       {  0., -1., -1. },
          			       {  0.,  1., -1. },
          			       {  0., -1.,  1. } } ;

/** Coefficients for pseudo-equilibrium distribution of the D3Q18 model */
static double d3q18_coefficients[18][4] = { { 1./12.,  1./6., 1./4., -1./6. },
                                            { 1./12.,  1./6., 1./4., -1./6. },
                                            { 1./12.,  1./6., 1./4., -1./6. },
                                            { 1./12.,  1./6., 1./4., -1./6. },
                                            { 1./12.,  1./6., 1./4., -1./6. },
                                            { 1./12.,  1./6., 1./4., -1./6. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. },
				            { 1./24., 1./12., 1./8., 1./12. } } ;

/** The model used is D3Q18 by default. */
static LB_Model lbmodel = { 18, d3q18_lattice, d3q18_coefficients } ;

/** The number of velocities of the LB model.
 * Shortcut for lbmodel.n_veloc */
static int n_veloc;

/** The number of field variables on a local lattice site (counted in doubles). */
static int n_fields;

/** Flag indicating whether fluctuations are present. */
static int fluct;

/** \name Derived parameters */
/*@{*/
/** eigenvalue of the collision operator for relaxation of shear modes */
static double lblambda;
/** amplitude of the fluctuations in the fluid stress tensor */
static double lb_fluct_pref;
/** amplitude of the fluctuations in the viscous coupling */
static double lb_coupl_pref;
/*@}*/

/** Lattice spacing.
 * This variable is used for convenience instead of having to type lattice.agrid everywhere. */
static double agrid;
/** Lattice Boltzmann time step
 * This variable is used for convenience instead of having to type lattice.tau everywhere. */
static double tau;

/** measures the MD time since the last fluid update */
static double fluidstep=0.0;

/***********************************************************************/

/** Performs basic sanity checks. */
static int lb_sanity_checks() {

  char *errtxt;
  int ret = 0;

    if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
      errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{103 LB requires domain-decomposition cellsystem} ");
      ret = -1;
    } 
    else if (dd.use_vList) {
      errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{104 LB requires no Verlet Lists} ");
      ret = -1;
    }    

    return ret;

}

static void compare_buffers(double *buf1, double *buf2, int size) {
  if (memcmp(buf1,buf2,size)) {
    char *errtxt;
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{102 Halo buffers are not identical} ");
  }
}

/** Checks consistency of the halo regions (ADDITIONAL_CHECKS)
 * This function can be used as an additional check. It test whether the 
 * halo regions have been exchanged correctly. */
static void lb_check_halo_regions() {

  int x,y,z, index, s_node, r_node;
  double *s_buffer, *r_buffer;
  MPI_Status status[2];

  r_buffer = malloc(n_fields*sizeof(double));

  for (y=0;y<lblattice.halo_grid[1];++y) {
    for (z=0;z<lblattice.halo_grid[2];++z) {

      index  = get_linear_index(0,y,z,lblattice.halo_grid);
      s_buffer = lbfluid[index].n;
      s_node = node_neighbors[1];
      r_node = node_neighbors[0];
      MPI_Sendrecv(s_buffer, n_fields, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		   r_buffer, n_fields, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		   MPI_COMM_WORLD, status);
      index = get_linear_index(lblattice.grid[0],y,z,lblattice.halo_grid);
      compare_buffers(lbfluid[index].n,r_buffer,n_fields*sizeof(double));

      index = get_linear_index(lblattice.grid[0]+1,y,z,lblattice.halo_grid);
      s_buffer = lbfluid[index].n;
      s_node = node_neighbors[0];
      r_node = node_neighbors[1];
      MPI_Sendrecv(s_buffer, n_fields, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		   r_buffer, n_fields, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		   MPI_COMM_WORLD, status);
      index = get_linear_index(1,y,z,lblattice.halo_grid);
      compare_buffers(lbfluid[index].n,r_buffer,n_fields*sizeof(double));

    }
  }

  for (x=0;x<lblattice.grid[0];++x) {
    for (z=0;z<lblattice.grid[2];++z) {

      index = get_linear_index(x,0,z,lblattice.halo_grid);
      s_buffer = lbfluid[index].n;
      s_node = node_neighbors[3];
      r_node = node_neighbors[2];
      MPI_Sendrecv(s_buffer, n_fields, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		   r_buffer, n_fields, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		   MPI_COMM_WORLD, status);
      index = get_linear_index(x,lblattice.grid[1],z,lblattice.halo_grid);
      compare_buffers(lbfluid[index].n,r_buffer,n_fields*sizeof(double));

      index = get_linear_index(x,lblattice.grid[1]+1,z,lblattice.halo_grid);
      s_buffer = lbfluid[index].n;
      s_node = node_neighbors[2];
      r_node = node_neighbors[3];
      MPI_Sendrecv(s_buffer, n_fields, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		   r_buffer, n_fields, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		   MPI_COMM_WORLD, status);
      index = get_linear_index(x,1,z,lblattice.halo_grid);
      compare_buffers(lbfluid[index].n,r_buffer,n_fields*sizeof(double));

    }
  }

  for (x=0;x<lblattice.grid[0];++x) {
    for (y=0;y<lblattice.grid[1];++y) {

      index = get_linear_index(x,y,0,lblattice.halo_grid);
      s_buffer = lbfluid[index].n;
      s_node = node_neighbors[5];
      r_node = node_neighbors[4];
      MPI_Sendrecv(s_buffer, n_fields, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		   r_buffer, n_fields, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		   MPI_COMM_WORLD, status);
      index = get_linear_index(x,y,lblattice.grid[2],lblattice.halo_grid);
      compare_buffers(lbfluid[index].n,r_buffer,n_fields*sizeof(double));

      index = get_linear_index(x,y,lblattice.grid[2]+1,lblattice.halo_grid);
      s_buffer = lbfluid[index].n;
      s_node = node_neighbors[4];
      r_node = node_neighbors[5];
      MPI_Sendrecv(s_buffer, n_fields, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
		   r_buffer, n_fields, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
		   MPI_COMM_WORLD, status);
      index = get_linear_index(x,y,1,lblattice.halo_grid);
      compare_buffers(lbfluid[index].n,r_buffer,n_fields*sizeof(double));

    }
  }

  free(r_buffer);

}    

/***********************************************************************/

/** (Re-)allocate memory for the fluid and initialize pointers. */
static void lb_create_fluid() {

  int index;

  lblattice.fields = realloc(lblattice.fields,lblattice.halo_grid_volume*sizeof(LB_FluidNode));
  lblattice.data = realloc(lblattice.data,lblattice.halo_grid_volume*n_fields*sizeof(double));

  lbfluid = lblattice.fields;
  lbfluid[0].n = lblattice.data;
  
  for (index=0; index<lblattice.halo_grid_volume; index++) {
    lbfluid[index].n   = lbfluid[0].n + index*n_fields;
    lbfluid[index].rho = lbfluid[index].n + n_veloc;
    lbfluid[index].j   = lbfluid[index].rho + 1;
    lbfluid[index].pi  = lbfluid[index].j + SPACE_DIM;
#ifndef D3Q18
    lbfluid[index].n_new = lbfluid[index].pi + SPACE_DIM*(SPACE_DIM+1)/2;
#endif
  }

}

/** Sets up the structures for exchange of the halo regions.
 *  See also \ref halo.c */
static void lb_prepare_communication() {

    /* create types for lattice data layout */
    int lens[2] = { n_veloc, 1 };
    int disps[2] = { 0, n_fields*sizeof(double) };
    MPI_Aint adisps[2] = { 0, n_fields*sizeof(double) };
    MPI_Datatype types[2] = { MPI_DOUBLE, MPI_UB };
    MPI_Datatype datatype;
    MPI_Type_struct(2, lens, adisps, types, &datatype);
    MPI_Type_commit(&datatype);
    lens[0] *= sizeof(double);
    Fieldtype fieldtype;
    halo_create_fieldtype(1, lens, disps, disps[1], &fieldtype);

    /* setup the halo communication */
    prepare_halo_communication(&update_halo_comm,&lblattice,fieldtype,datatype);
 
    MPI_Type_free(&datatype);
    halo_free_fieldtype(&fieldtype);

}

/** Release the fluid. */
static void lb_release_fluid() {
  free(lbfluid[0].n);
  free(lbfluid);
}

/** (Re-)initializes the fluid. */
void lb_reinit_parameters() {

  agrid = lbpar.agrid;
  tau   = lbpar.tau;

  n_veloc = lbmodel.n_veloc;

  /* number of double entries in the data fields
   * velocity populations, density, momentum, stress tensor */
  n_fields = n_veloc + 1 + SPACE_DIM + SPACE_DIM*(SPACE_DIM+1)/2;
#ifndef D3Q18
  n_fields += n_veloc; /* temporary velocity populations */
#endif

  /* Eq. (3) Ahlrichs and Duenweg, JCP 111(17):8225 (1999). */
  lblambda = -2./(6.*lbpar.viscosity*tau/(agrid*agrid)+1.) ;
    
  if (temperature > 0.0) {  /* fluctuating hydrodynamics ? */

    fluct = 1 ;

    /* lb_fluct_pref is stored in lattice units (pressure)
     * Eq. (7) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
     * The factor 12 comes from the fact that we use random numbers
     * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
     */
    lb_fluct_pref = sqrt(12.*2.*lbpar.viscosity*lbpar.rho*temperature*SQR(lblambda)*tau*tau*tau/agrid);

    LB_TRACE(fprintf(stderr,"%d: lb_fluct_pref=%f (temp=%f, lambda=%f, tau=%f, agrid=%f)\n",this_node,lb_fluct_pref,temperature,lblambda,tau,agrid));

    /* lb_coupl_pref is stored in MD units (force)
     * Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
     * The factor 12 comes from the fact that we use random numbers
     * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
     * time_step comes from the discretization.
     */
    lb_coupl_pref = sqrt(12.*2.*lbpar.friction*temperature/time_step); 

    LB_TRACE(fprintf(stderr,"%d: lb_coupl_pref=%f (temp=%f, friction=%f, time_step=%f agrid=%f)\n",this_node,lb_coupl_pref,temperature,lbpar.friction,time_step,agrid));

  } else {
    /* no fluctuations at zero temperature */
    fluct = 0 ;
    lb_fluct_pref = 0.0;
    lb_coupl_pref = 0.0;
  }
  
  
}

/** (Re-)initializes the fluid according to the given value of rho. */
void lb_reinit_fluid() {

    int k ;

    /* default values for fields in lattice units */
    double rho = lbpar.rho*agrid*agrid*agrid ;
    double v[3] = { 0., 0., 0. };

    for (k=0;k<lblattice.halo_grid_volume;k++)

      if (lbfluid[k].boundary==0) {
	lb_set_local_fields(k,rho,v) ;
      } else {
	lb_set_local_fields(k,0.0,v);
      }

}

/** Performs a full initialization of
 *  the Lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values. */
void lb_init() {

  if (lb_sanity_checks()) return;

  /* initialize derived parameters */
  lb_reinit_parameters();

  /* initialize the local lattice domain */
  init_lattice(&lblattice,agrid,tau);  

  if (check_runtime_errors()) return;

  /* allocate memory for data structures */
  lb_create_fluid();

#ifdef CONSTRAINTS
  /* setup boundaries of constraints */
  lb_init_constraints();
#endif

  /* setup the initial particle velocity distribution */
  lb_reinit_fluid();

  /* prepare the halo communication */
  lb_prepare_communication();

  /* communicate the initial halo regions */
  halo_communication(&update_halo_comm) ;

}

/** Release fluid and communication. */
void lb_release() {
  release_halo_communication(&update_halo_comm);
  lb_release_fluid();
}

/***********************************************************************/

/** Calculate the average density of the fluid in the system.
 * This function has to be called after changing the density of
 * a local lattice site in order to set lbpar.rho consistently. */
void lb_calc_average_rho() {

  int x, y, z, index;
  double rho, sum_rho;

  rho = 0.0;
  for (x=1; x<=lblattice.grid[0]; x++) {
      for (y=1; y<=lblattice.grid[1]; y++) {
	  for (z=1; z<=lblattice.grid[2]; z++) {
	      index = get_linear_index(x,y,z,lblattice.halo_grid);

	      lb_calc_local_rho(&lbfluid[index]);
	      rho += *lbfluid[index].rho;

	  }
      }
  }

  MPI_Allreduce(&rho, &sum_rho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* calculate average density in MD units */
  lbpar.rho = sum_rho / (box_l[0]*box_l[1]*box_l[2]);

}

/** Returns the mass, momentum and stress of a local lattice site.
 * @param index The index of the lattice site within the local domain (Input)
 * @param rho   Local density of the fluid (Output)
 * @param j     Local momentum of the fluid (Output)
 * @param pi    Local stress tensor of the fluid (Output)
 */
void lb_get_local_fields(int index, double *rho, double *j, double *pi) {

  int i,k,m;

  double *local_rho = lbfluid[index].rho;
  double *local_j   = lbfluid[index].j;
  double *local_pi  = lbfluid[index].pi;

  lb_calc_local_fields(&lbfluid[index],1);

  *rho = *local_rho;
  m = 0;
  for (i=0;i<3;i++) {
    j[i] = local_j[i];
    for (k=0;k<i;k++) {
      pi[m] = local_pi[m];
      m++;
    }
  }

}

/** Sets the density and momentum on a local lattice site.
 * @param index The index of the lattice site within the local domain (Input)
 * @param rho   Local density of the fluid (Input)
 * @param v     Local velocity of the fluid (Input)
 */
void lb_set_local_fields(int index, const double rho, double *v) {

  int i ;

  double *local_n = lbfluid[index].n;
  double (*c)[3] = lbmodel.c ;
  double (*coeff)[4] = lbmodel.coeff ;

  *(lbfluid[index].rho) = rho;
  lbfluid[index].j[0] = rho * v[0];
  lbfluid[index].j[1] = rho * v[1];
  lbfluid[index].j[2] = rho * v[2];

  /* Eq. (2) Ahlrichs and Duenweg, JCP 111(17):8225 (1999). */
  for (i=0;i<n_veloc;i++) {

    local_n[i] = coeff[i][0] * rho ;

    local_n[i] += coeff[i][1] * rho * scalar(v,c[i]) ;

    local_n[i] += coeff[i][2] * rho * (scalar(v,c[i])*scalar(v,c[i])-scalar(v,v)*scalar(c[i],c[i])/SPACE_DIM) ;

    local_n[i] += coeff[i][3] * rho * scalar(v,v) ;

  }

}

/***********************************************************************/
/** \name External forces */
/***********************************************************************/
/*@{*/

/** Apply external forces to the fluid.
 * Eq. (28) Ladd and Verberg, J. Stat. Phys. 104(5/6):1191 (2001).
 * Note that the second moment of the force is neglected.
 */
MDINLINE void lb_external_forces() {

#ifdef D3Q18
  int x, y, z, index;
  double delta_j[3] = { 0.0, 0.0, 0.0 };
  double *local_n;
  
  /* calculate momentum due to ext_force in lattice units */
  delta_j[0] = lbpar.ext_force[0]*tau*tau/agrid;
  delta_j[1] = lbpar.ext_force[1]*tau*tau/agrid;
  delta_j[2] = lbpar.ext_force[2]*tau*tau/agrid;
  
  for (x=1; x<=lblattice.grid[0]; x++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
      for (z=1; z<=lblattice.grid[2]; z++) {
	
	index = get_linear_index(x,y,z,lblattice.halo_grid);
	
	if (lbfluid[index].boundary==0) {
	  local_n = lbfluid[index].n;
	  
	  local_n[0]  = local_n[0] + 1./6.*delta_j[0] ;
	  local_n[1]  = local_n[1] - 1./6.*delta_j[0] ;
	  local_n[2]  = local_n[2] + 1./6.*delta_j[1] ;
	  local_n[3]  = local_n[3] - 1./6.*delta_j[1] ;
	  local_n[4]  = local_n[4] + 1./6.*delta_j[2] ;
	  local_n[5]  = local_n[5] - 1./6.*delta_j[2] ;
	  local_n[6]  = local_n[6] + 1./12.*(delta_j[0]+delta_j[1]) ;
	  local_n[7]  = local_n[7] - 1./12.*(delta_j[0]+delta_j[1]) ;
	  local_n[8]  = local_n[8] + 1./12.*(delta_j[0]-delta_j[1]) ;
	  local_n[9]  = local_n[9] - 1./12.*(delta_j[0]-delta_j[1]) ;
	  local_n[10] = local_n[10] + 1./12.*(delta_j[0]+delta_j[2]) ;
	  local_n[11] = local_n[11] - 1./12.*(delta_j[0]+delta_j[2]) ;
	  local_n[12] = local_n[12] + 1./12.*(delta_j[0]-delta_j[2]) ;
	  local_n[13] = local_n[13] - 1./12.*(delta_j[0]-delta_j[2]) ;
	  local_n[14] = local_n[14] + 1./12.*(delta_j[1]+delta_j[2]) ;
	  local_n[15] = local_n[15] - 1./12.*(delta_j[1]+delta_j[2]) ;
	  local_n[16] = local_n[16] + 1./12.*(delta_j[1]-delta_j[2]) ;
	  local_n[17] = local_n[17] - 1./12.*(delta_j[1]-delta_j[2]) ;
	    
	}

      }
    }
  }
#else
#error External forces are only implemented for D3Q18!
#endif
}

/*@}*/

/***********************************************************************/
/** \name Collision step */
/***********************************************************************/
/*@{*/

/** Collision update of the stress tensor.
 * PI is updated according to PI_eq + (1-lambda)*(PI_neq-trace).
 * Eq. (2.14) of Ladd, J. Fluid Mech. 271, 285-309 (1994)
 * <br><em>Remarks:</em>
 * <ul>
 * <li>PI is not made traceless here. This will be taken care of in \ref lb_update_local_n.</li>
 *<li>All terms rhoc_sq (from PI_eq) can be dropped because they cancel when PI is made traceless.</li>
 * </ul>
 *
 * @param local_node Pointer to the local lattice site (Input).
 * @param trace      Trace of local stress tensor (Output).
 * @param trace_eq   Trace of equilibrium part of local stress tensor (Output).
 */
MDINLINE void lb_update_local_pi(LB_FluidNode *local_node, double *trace, double *trace_eq) {

  int k,l,m ;
  double onepluslambda = 1.0 + lblambda ;

  double *local_pi = local_node->pi;

#ifdef CREEPINGFLOW
  *trace = 0.0 ;
  m = 0 ;
  for (k=0;k<SPACE_DIM;k++) {
    for (l=0;l<k;l++) {
      /* non-diagonal elements */
      local_pi[m] = onepluslambda * local_pi[m] ;
      m++ ;
    }
    /* diagonal elements */
    /* The full formula would be:
     * local_pi[m] = rhoc_sq + onepluslambda * (local_pi[m] - rhoc_sq - trace)
     * We can drop rhoc_sq because the traceless part is zero.
     * The trace will be taken care of in lb_update_n.
     */
    local_pi[m] = onepluslambda * local_pi[m] ;
    /* calculate the trace on the fly */
    *trace += local_pi[m] ;
    m++ ;
  }
#else
  double local_pi_eq[6] ;

  trace = 0.0 ;
  trace_eq = 0.0 ;
  m = 0 ;
  for (k=0;k<SPACE_DIM;k++) {
    tmp = local_j[k]/local_rho ;
    for (l=0;l<k;l++) {
      /* non-diagonal elements */
      local_pi_eq[m] = tmp * local_j[l] ;
      local_pi[m] = local_pi_eq[m] + onepluslambda * (local_pi[m] - local_pi_eq[m]) ;
      m++ ;
    }
    /* diagonal elements */
    /* The full formulas would be:
     * local_pi_eq[m] = rhoc_sq + tmp * local_j[k]
     * local_pi[m] = local_pi_eq[m] + onepluslambda * (local_pi[m] - local_pi_eq[m] - trace)
     * We can drop rhoc_sq because the traceless part is zero.
     * The trace will be taken care of in lb_update_n.
     */
    local_pi_eq[m] = tmp * local_j[k] ;
    local_pi[m] = local_pi_eq[m] + onepluslambda * (local_pi[m] - local_pi_eq[m]) ;
    /* calculate the traces on the fly */
    *trace_eq += local_pi_eq[m] ;
    *trace += local_pi[m] ;
    m++ ;
  }
#endif

}

/** Add fluctuating part to the stress tensor and update the populations.
 *
 * Ladd, J. Fluid Mech. 271, 285-309 (1994).<br>
 * Ahlrichs, PhD-Thesis (2000).
 *   
 * @param local_node Pointer to the local lattice site.
 * @param trace      Trace of local stress tensor (Input).
 * @param trace_eq   Trace of equilibrium part of the stress tensor (Input).
 * @param badrandoms Flag/Counter for the occurence of negative populations (Output).
 */
MDINLINE void lb_add_fluct_update_local_n(LB_FluidNode *local_node, const double trace, const double trace_eq, int *badrandoms) {

  int i,k,l,m ;
  double tmp[3],sum ;

  double *local_n   = local_node->n;
  double *local_rho = local_node->rho;
  double *local_j   = local_node->j;
  double *local_pi  = local_node->pi;

  sum = 0.0 ;
  m = 0 ;
  for (k=0;k<SPACE_DIM;k++) {
    tmp[k] = sqrt(2) * lb_fluct_pref * (d_random()-0.5) ;
    sum -= tmp[k] ;
    for (l=0;l<k;l++) {
      /* non-diagonal elements */
      local_pi[m] -= lb_fluct_pref * (d_random()-0.5) ;
      m++ ;     
    }
    /* diagonal elements */
    local_pi[m] -= tmp[k] ;
    m++ ;
  }

  /* subtract the sum from the diagonal elements */
  sum /= SPACE_DIM ;
  m = 0 ;
  for (k=0;k<SPACE_DIM;k++) {
    local_pi[k+m] -= sum ;
    m += k+1 ;
  }

#ifdef D3Q18
  double tmp1,tmp2;

  /* For the D3Q18 model the c_i and the coefficients are known.
   * We use the explicit expressions in order to save
   * multiplications with 0 or 1.
   * The expressions can be derived from defaultc_g.
   */

  double trace_by_spacedim = trace/SPACE_DIM ;

  /* First we update the |c_i|=1 sublattice for which
   * coeff[i][0] = 1./12.
   * coeff[i][1] = 1./6.
   * coeff[i][2] = 1./4.
   * coeff[i][3] = -1./6.
   */

  double local_rho_times_coeff = 1./12. * *local_rho ;

  /* Take care to make local_pi traceless. */
  local_n[0] = local_rho_times_coeff + 1./6.*local_j[0] + 1./4.*(local_pi[0]-trace_by_spacedim) ;
  local_n[1] = local_rho_times_coeff - 1./6.*local_j[0] + 1./4.*(local_pi[0]-trace_by_spacedim) ;
  local_n[2] = local_rho_times_coeff + 1./6.*local_j[1] + 1./4.*(local_pi[2]-trace_by_spacedim) ;
  local_n[3] = local_rho_times_coeff - 1./6.*local_j[1] + 1./4.*(local_pi[2]-trace_by_spacedim) ;
  local_n[4] = local_rho_times_coeff + 1./6.*local_j[2] + 1./4.*(local_pi[5]-trace_by_spacedim) ;
  local_n[5] = local_rho_times_coeff - 1./6.*local_j[2] + 1./4.*(local_pi[5]-trace_by_spacedim) ;


#ifndef CREEPINGFLOW
  for (i=0;i<6;i++) {
    /* The full formula would be:
     * n[1] += -1./6.*(trace-3.0*rhoc_sq)
     * The nonequilibrium part of PI is traceless.
     * In the equilibrium part the rhoc_sq terms have been dropped,
     * which would cancel here with -3.0*rhoc_sq
     * so the latter can be dropped here.
     */
    local_n[i] += -1./6.*trace_eq ;
  }
#endif

  /* Check for negative populations */
  for (i=0;i<6;i++) {
    if (local_n[i]<0.0) {
      (*badrandoms)++;
      LB_TRACE(fprintf(stderr,"%d: population %d negative %f (local_rho=%.3f, local_j=(%.3f,%.3f,%.3f) local_pi=(%.3f,%.3f,%.3f,%.3f,%.3f,%.3f)\n",this_node,i,local_n[i],*local_rho,local_j[0],local_j[1],local_j[2],local_pi[0],local_pi[1],local_pi[2],local_pi[3],local_pi[4],local_pi[5]));
      return ;
    }
  }

  /* Second we update the |c_i|=sqrt(2) sublattice for which
   * coeff[i][0] = 1./24.
   * coeff[i][1] = 1./12.
   * coeff[i][2] = 1./8.
   * coeff[i][3] = 1./12.
   */

  local_rho_times_coeff = 1./24. * *local_rho ;

  /* Take care to make local_pi traceless. */
  tmp1 = local_pi[0]-trace_by_spacedim + local_pi[2]-trace_by_spacedim ; 
  tmp2 = 2.0 * local_pi[1] ;

  local_n[6] = local_rho_times_coeff + 1./12.*(local_j[0]+local_j[1]) + 1./8.*(tmp1+tmp2) ;
  local_n[7] = local_rho_times_coeff - 1./12.*(local_j[0]+local_j[1]) + 1./8.*(tmp1+tmp2) ;
  local_n[8] = local_rho_times_coeff + 1./12.*(local_j[0]-local_j[1]) + 1./8.*(tmp1-tmp2) ;
  local_n[9] = local_rho_times_coeff - 1./12.*(local_j[0]-local_j[1]) + 1./8.*(tmp1-tmp2) ;

  /* Take care to make local_pi traceless. */
  tmp1 = local_pi[0]-trace_by_spacedim + local_pi[5]-trace_by_spacedim ;
  tmp2 = 2.0 * local_pi[3] ;

  local_n[10] = local_rho_times_coeff + 1./12.*(local_j[0]+local_j[2]) + 1./8.*(tmp1+tmp2) ;
  local_n[11] = local_rho_times_coeff - 1./12.*(local_j[0]+local_j[2]) + 1./8.*(tmp1+tmp2) ;
  local_n[12] = local_rho_times_coeff + 1./12.*(local_j[0]-local_j[2]) + 1./8.*(tmp1-tmp2) ;
  local_n[13] = local_rho_times_coeff - 1./12.*(local_j[0]-local_j[2]) + 1./8.*(tmp1-tmp2) ;

  /* Take care to make local_pi traceless. */
  tmp1 = local_pi[2]-trace_by_spacedim + local_pi[5]-trace_by_spacedim ;
  tmp2 = 2.0 * local_pi[4] ;

  local_n[14] = local_rho_times_coeff + 1./12.*(local_j[1]+local_j[2]) + 1./8.*(tmp1+tmp2) ;
  local_n[15] = local_rho_times_coeff - 1./12.*(local_j[1]+local_j[2]) + 1./8.*(tmp1+tmp2) ;
  local_n[16] = local_rho_times_coeff + 1./12.*(local_j[1]-local_j[2]) + 1./8.*(tmp1-tmp2) ;
  local_n[17] = local_rho_times_coeff - 1./12.*(local_j[1]-local_j[2]) + 1./8.*(tmp1-tmp2) ;

#ifndef CREEPINGFLOW
  for (i=6;i<18;i++) {
    /* The full formula would be:
     * n[1] += 1./12.*(trace-3.0*rhoc_sq)
     * The nonequilibrium part of PI is traceless.
     * In the equilibrium part the rhoc_sq terms have been dropped,
     * which would cancel here with -3.0*rhoc_sq
     * so the latter can be dropped here.
     */
    local_n[i] += 1./12.*trace_eq ;
  }
#endif

  /* Check for negative populations. */
  for (i=6;i<18;i++) {
    if (local_n[i]<0.0) {
      (*badrandoms)++ ;
      LB_TRACE(fprintf(stderr,"%d: population %d negative %f (local_rho=%.3f, local_j=(%.3f,%.3f,%.3f) local_pi=(%.3f,%.3f,%.3f,%.3f,%.3f,%.3f)\n",this_node,i,local_n[i],*local_rho,local_j[0],local_j[1],local_j[2],local_pi[0],local_pi[1],local_pi[2],local_pi[3],local_pi[4],local_pi[5]));
      return ;
    }
  }

#else
  l = 0 ;
  for (i=0;i<n_abs_veloc;i++) {
    for (j=0;j<n_abs_which[i];j++) {
      for (k=0;k<4;k++) {
	coeff[l][k] = defaultcoef_eq[4*i+k] ;
      }
      l++ ;
    }
  }

  for (i=0;i<n_veloc;i++) {

    local_n[i] = coeff[i][0] * local_rho ;

    tmp = 0.0 ;
    for (k=0;k<SPACE_DIM;k++) {
      tmp += local_j[k] * c[i][k] ;
    }
    local_n[i] += coeff[i][1] * tmp ;

    tmp = 0.0 ;
    m = 0 ;
    for (k=0;k<SPACE_DIM;k++) {
      for (l=0;l<k;k++) {
	/* non-diagonal elements */
	tmp += 2.0 * local_pi[m] * c[i][k] * c[i][l] ; 
	m++ ;
      }
      /* diagonal elements */
      /* Here finally we take care to make local_pi traceless. */
      tmp += (local_pi[m]-trace_by_spacedim) * c[i][k] * c[i][k] ;
      m++ ;
    }
    local_n[i] += coeff[i][2] * tmp ;

#ifndef CREEPINGFLOW
    /* The full formula would be:
     * n[i] += coeff[i][3] * (trace - 3.0*rhoc_sq)
     * The nonequilibrium part of PI is traceless.
     * In the equilibrium part the rhoc_sq terms have been dropped,
     * which would cancel here with - 3.0*rhoc_sq
     * so the latter can be dropped here.
     */
    local_n[i] += coeff[i][3] * trace_eq ;
#endif

  }

  /* Check for negative populations */
  for (i=0;i<n_veloc;i++) {
    if (local_n[i]<0.0) {
      (*badrandoms)++ ;
      LB_TRACE(fprintf(stderr,"%d: population %d negative %f (local_rho=%.3f, local_j=(%.3f,%.3f,%.3f) local_pi=(%.3f,%.3f,%.3f,%.3f,%.3f,%.3f)\n",this_node,i,local_n[i],*local_rho,local_j[0],local_j[1],local_j[2],local_pi[0],local_pi[1],local_pi[2],local_pi[3],local_pi[4],local_pi[5]));
      return ;
    }
  }

#endif

}


/** Update the local populations without fluctuations.
 *
 * Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
 *
 * In the double tensor contraction, it does not matter 
 * which of the two tensors is made traceless.
 * Hence instead of using the traceless part of c_ic_i,
 * the traceless part of PI is contracted with the full c_ic_i.
 * <em>Remember: PI is not traceless yet and rhoc_sq terms have been dropped (see \ref lb_update_local_pi).</em>
 *
 * @param local_node Pointer to the local lattice site (Input).
 * @param trace      Trace of the local stress tensor (Input).
 * @param trace_eq   Trace of equilibriumd part of local stress tensor (Input).
 */
MDINLINE void lb_update_local_n(LB_FluidNode *local_node, const double trace, const double trace_eq) {

  double *local_n   = local_node->n;
  double *local_rho = local_node->rho;
  double *local_j   = local_node->j;
  double *local_pi  = local_node->pi;

  double trace_by_spacedim = trace/SPACE_DIM ;

#ifdef D3Q18
  double tmp1,tmp2;

  /* First we update the |c_i|=1 sublattice for which
   * coeff[0] = 1./12.
   * coeff[1] = 1./6.
   * coeff[2] = 1./4.
   * coeff[3] = -1./6.
   */
  double local_rho_times_coeff = 1./12. * *local_rho ;

  local_n[0] = local_rho_times_coeff ;
  local_n[1] = local_rho_times_coeff ;
  local_n[2] = local_rho_times_coeff ;
  local_n[3] = local_rho_times_coeff ;
  local_n[4] = local_rho_times_coeff ;
  local_n[5] = local_rho_times_coeff ;

  /* Take care to make local_pi traceless. */
  local_n[0] += 1./4.*(local_pi[0]-trace_by_spacedim) + 1./6.*local_j[0] ;
  local_n[1] += 1./4.*(local_pi[0]-trace_by_spacedim) - 1./6.*local_j[0] ;
  local_n[2] += 1./4.*(local_pi[2]-trace_by_spacedim) + 1./6.*local_j[1] ;
  local_n[3] += 1./4.*(local_pi[2]-trace_by_spacedim) - 1./6.*local_j[1] ;
  local_n[4] += 1./4.*(local_pi[5]-trace_by_spacedim) + 1./6.*local_j[2] ;
  local_n[5] += 1./4.*(local_pi[5]-trace_by_spacedim) - 1./6.*local_j[2] ;

#ifndef CREEPINGFLOW
  int i;
  for (i=0;i<6;i++) {
    /* The full formula would be:
     * n[1] += -1./6.*(trace-3.0*rhoc_sq)
     * The nonequilibrium part of PI is traceless.
     * In the equilibrium part the rhoc_sq terms have been dropped,
     * which would cancel here with -3.0*rhoc_sq
     * so the latter can be dropped here.
     */
    local_n[i] += -1./6.*trace_eq ;
  }
#endif

  /* Second we update the |c_i|=sqrt(2) sublattice for which
   * coeff[0] = 1./24.
   * coeff[1] = 1./12.
   * coeff[2] = 1./8.
   * coeff[3] = 1./12.
   */
  local_rho_times_coeff = 1./24. * *local_rho ;

  /* Take care to make local_pi traceless. */
  tmp1 = local_pi[0]-trace_by_spacedim + local_pi[2]-trace_by_spacedim ;
  tmp2 = 2.0*local_pi[1] ;

  local_n[6] = local_rho_times_coeff ;
  local_n[7] = local_rho_times_coeff ; 
  local_n[8] = local_rho_times_coeff ; 
  local_n[9] = local_rho_times_coeff ;
 
  local_n[6] += 1./8.*(tmp1+tmp2) + 1./12.*(local_j[0]+local_j[1]) ;
  local_n[7] += 1./8.*(tmp1+tmp2) - 1./12.*(local_j[0]+local_j[1]) ;
  local_n[8] += 1./8.*(tmp1-tmp2) + 1./12.*(local_j[0]-local_j[1]) ;
  local_n[9] += 1./8.*(tmp1-tmp2) - 1./12.*(local_j[0]-local_j[1]) ;

  /* Take care to make local_pi traceless. */
  tmp1 = local_pi[0]-trace_by_spacedim + local_pi[5]-trace_by_spacedim ;
  tmp2 = 2.0*local_pi[3] ;

  local_n[10] = local_rho_times_coeff ;
  local_n[11] = local_rho_times_coeff ; 
  local_n[12] = local_rho_times_coeff ; 
  local_n[13] = local_rho_times_coeff ;
 
  local_n[10] += 1./8.*(tmp1+tmp2) + 1./12.*(local_j[0]+local_j[2]) ;
  local_n[11] += 1./8.*(tmp1+tmp2) - 1./12.*(local_j[0]+local_j[2]) ;
  local_n[12] += 1./8.*(tmp1-tmp2) + 1./12.*(local_j[0]-local_j[2]) ;
  local_n[13] += 1./8.*(tmp1-tmp2) - 1./12.*(local_j[0]-local_j[2]) ;

  /* Take care to make local_pi traceless. */
  tmp1 = local_pi[2]-trace_by_spacedim + local_pi[5]-trace_by_spacedim ;
  tmp2 = 2.0*local_pi[4] ;

  local_n[14] = local_rho_times_coeff ;
  local_n[15] = local_rho_times_coeff ; 
  local_n[16] = local_rho_times_coeff ; 
  local_n[17] = local_rho_times_coeff ;
 
  local_n[14] += 1./8.*(tmp1+tmp2) + 1./12.*(local_j[1]+local_j[2]) ;
  local_n[15] += 1./8.*(tmp1+tmp2) - 1./12.*(local_j[1]+local_j[2]) ;
  local_n[16] += 1./8.*(tmp1-tmp2) + 1./12.*(local_j[1]-local_j[2]) ;
  local_n[17] += 1./8.*(tmp1-tmp2) - 1./12.*(local_j[1]-local_j[2]) ;

#ifndef CREEPINGFLOW
  for (i=6;i<18;i++) {
    /* The full formula would be:
     * n[1] += 1./12.*(trace-3.0*rhoc_sq)
     * The nonequilibrium part of PI is traceless.
     * In the equilibrium part the rhoc_sq terms have been dropped,
     * which would cancel here with -3.0*rhoc_sq
     * so the latter can be dropped here.
     */
    local_n[i] += 1./12.*trace_eq ;
  }
#endif

#else
  int i;
  for (i=0;i<n_veloc;i++) {

    local_n[i] = coeff[0] * local_rho ;

    tmp = 0.0 ;
    for (k=0;k<SPACE_DIM;k++) {
      tmp += local_j[k] * c[i][k] ;
    }
    local_n[i] += coeff[1] * tmp ;

    tmp = 0.0 ;
    m = 0 ;
    for (k=0;k<SPACE_DIM;k++) {
      for (l=0;l<k;k++) {
	/* non-diagonal elements */
	tmp += 2.0 * local_pi[m] * c[i][k] * c[i][l] ; 
	m++ ;
      }
      /* diagonal elements */
      /* Here finally we take care to make local_pi traceless. */
      tmp += (local_pi[m]-trace/SPACE_DIM) * c[i][k] * c[i][k] ;
      m++ ;
    }
    local_n[i] += coeff[2] * tmp ;

#ifndef CREEPINGFLOW
    /* The full formula would be:
     * n[i] += coeff[3] * (trace - 3.0*rhoc_sq)
     * The nonequilibrium part of PI is traceless.
     * In the equilibrium part the rhoc_sq terms have been dropped,
     * which would cancel here with - 3.0*rhoc_sq
     * so the latter can be dropped here.
     */
    local_n[i] += coeff[3] * trace_eq ;
#endif

  }

#endif

#ifdef ADDITIONAL_CHECKS
  int j;
  for (j=0;j<18;j++) {
    if (local_n[j]<0.0) {
      char *errtxt;
      errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt,"{105 Unexpected negative population} ");
    }
  }
#endif

}

/** The Lattice Boltzmann collision step.
 * Loop over all lattice sites and perform the collision update.
 * If fluctuations are present, the fluctuating part of the stress tensor
 * is added. The update is only accepted then, if no negative populations
 * occur.
 */
MDINLINE void lb_calc_collisions() {

  int x, y, z ;
  int i;
  int index ;
  int badrandoms ;
  double *local_n , *local_rho, *local_j, *local_pi ;
  double trace, trace_eq ;
  double save_local_n[n_veloc], save_local_pi[6] ; 

  /* loop over all nodes (halo excluded) */
  for (z=1;z<=lblattice.grid[2];z++) {
    for (y=1;y<=lblattice.grid[1];y++) {
      for (x=1;x<=lblattice.grid[0];x++) {

        index = get_linear_index(x,y,z,lblattice.halo_grid) ;
	local_n   = lbfluid[index].n;
	local_rho = lbfluid[index].rho;
	local_j   = lbfluid[index].j;
	local_pi  = lbfluid[index].pi;

	lb_calc_local_fields(&lbfluid[index],1) ;

#ifdef ADDITIONAL_CHECKS
	double old_rho = *local_rho;
#endif

	lb_update_local_pi(&lbfluid[index],&trace,&trace_eq) ;
	
	if (fluct) {

	  /* save the local population */
	  for (i=0;i<n_veloc;i++) {
	    save_local_n[i] = local_n[i] ;
	  }
	  /* save the local pressure tensor */
	  for (i=0;i<6;i++) {
	    save_local_pi[i] = local_pi[i] ;
	  }

	  do { /* try random numbers until no negative populations occur */
	    
	    badrandoms = 0 ;

	    lb_add_fluct_update_local_n(&lbfluid[index],trace,trace_eq,&badrandoms) ;

	    if (badrandoms>0) {
	      fprintf(stderr,"%d: Negative population (badrandoms=%d). Check your parameters if this happens too often!\n",this_node,badrandoms);
	      LB_TRACE(fprintf(stderr,"negative population at site (%d,%d,%d) %d (badrandoms=%d)\n",x,y,z,index,badrandoms));
	      /* restore the local population */
	      for (i=0;i<n_veloc;i++) {
		local_n[i] = save_local_n[i] ;
	      }
	      /* restore the local pressure tensor */
	      for (i=0;i<6;i++) {
		local_pi[i] = save_local_pi[i] ;
	      }
	    }
	    
	  } while (badrandoms>0) ;

	} else {

	  lb_update_local_n(&lbfluid[index],trace,trace_eq) ;

#ifdef ADDITIONAL_CHECKS
	  for (i=0;i<18;i++) {
	    if (local_n[i] < 0.0) {
	      char *errtxt;
	      errtxt = runtime_error(128);
	      ERROR_SPRINTF(errtxt,"{106 Unexpected negative population} ");
	    }
	  }
#endif

	}

#ifdef ADDITIONAL_CHECKS
	lb_calc_local_rho(&lbfluid[index]);
	if (fabs(*local_rho-old_rho) > ROUND_ERROR_PREC) {
	  char *errtxt = runtime_error(128 + TCL_DOUBLE_SPACE + 3*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt,"{107 Mass loss/gain %le in lb_calc_collisions on site (%d,%d,%d)} ",*local_rho-old_rho,x,y,z);
	}
#endif

      }
    }
  }

}

/*@}*/

/***********************************************************************/
/** \name Streaming step */
/***********************************************************************/
/*@{*/

/** The Lattice Boltzmann streaming step.
 * The populations are moved to the neighbouring lattice sites
 * according to the velocity sublattice. This can be done in two ways:
 * First, one can use a temporary field to store the updated configuration.
 * Second, one can order the updates such that only populations are
 * overwritten which have already been propagated. The halo region
 * serves as a buffer. This requires two sweeps through the lattice,
 * one bottom up and one bottom down. One has to be careful if the
 * velocities are upgoing or downgoing. This can be a real bugfest!
 */
MDINLINE void lb_propagate_n() {

#ifdef D3Q18

  /* For the D3Q18 model, we can precalculate the index shifts
   * and then propagate top down (bottom up) for shifts to 
   * higher (lower) indices. The halo region is used as buffer.
   */

  int k;
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next0  =   1;                       // ( 1, 0, 0) +
  int next1  = - 1;                       // (-1, 0, 0)  
  int next2  =   yperiod;                 // ( 0, 1, 0) +
  int next3  = - yperiod;                 // ( 0,-1, 0)	 
  int next4  =   zperiod;                 // ( 0, 0, 1) +
  int next5  = - zperiod;                 // ( 0, 0,-1)	 
  int next6  =   (yperiod+1);             // ( 1, 1, 0) +
  int next7  = - (yperiod+1);             // (-1,-1, 0)	 
  int next8  =   (-yperiod+1);            // ( 1,-1, 0)  
  int next9  = - (-yperiod+1);            // (-1, 1, 0) +
  int next10 =   (zperiod+1);             // ( 1, 0, 1) +
  int next11 = - (zperiod+1);             // (-1, 0,-1)	 
  int next12 =   (-zperiod+1);            // ( 1, 0,-1)	 
  int next13 = - (-zperiod+1);            // (-1, 0, 1) +
  int next14 =   (yperiod+zperiod);       // ( 0, 1, 1) +
  int next15 = - (yperiod+zperiod);       // ( 0,-1,-1)	 
  int next16 =   (yperiod-zperiod);       // ( 0, 1,-1)	 
  int next17 = - (yperiod-zperiod);       // ( 0,-1, 1) +

  for (k=(lblattice.halo_grid_volume-lblattice.halo_offset);k>=0;k-=1) {

    /* top down propagation of populations to higher indices */
    lbfluid[k+next0].n[0]   = lbfluid[k].n[0] ;
    lbfluid[k+next2].n[2]   = lbfluid[k].n[2] ;
    lbfluid[k+next4].n[4]   = lbfluid[k].n[4] ;
    lbfluid[k+next6].n[6]   = lbfluid[k].n[6] ;
    lbfluid[k+next9].n[9]   = lbfluid[k].n[9] ;
    lbfluid[k+next10].n[10] = lbfluid[k].n[10] ;
    lbfluid[k+next13].n[13] = lbfluid[k].n[13] ;
    lbfluid[k+next14].n[14] = lbfluid[k].n[14] ;
    lbfluid[k+next17].n[17] = lbfluid[k].n[17] ;

  }

  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume;k+=1) {

    /* bottom up propagation of populations to lower indices */
    lbfluid[k+next1].n[1]   = lbfluid[k].n[1] ;
    lbfluid[k+next3].n[3]   = lbfluid[k].n[3] ;
    lbfluid[k+next5].n[5]   = lbfluid[k].n[5] ;
    lbfluid[k+next7].n[7]   = lbfluid[k].n[7] ;
    lbfluid[k+next8].n[8]   = lbfluid[k].n[8] ;
    lbfluid[k+next11].n[11] = lbfluid[k].n[11] ;
    lbfluid[k+next12].n[12] = lbfluid[k].n[12] ;
    lbfluid[k+next15].n[15] = lbfluid[k].n[15] ;
    lbfluid[k+next16].n[16] = lbfluid[k].n[16] ;

  }

#else

  double next, *tmp ;

  /* In the general case, we don't know a priori which 
   * velocities propagate to higher or lower indices.
   * So we use a complete new array as buffer and
   * swap the pointers afterwards.
   */ 

  /* calculate the index shift for all velocities */
  for (i=0;i<n_veloc;i++) {
    next[i] = (c[i][0]+yperiod*c[i][1]+zperiod*c[i][2]) ;
  }

  /* propagate the populations */
  /* on the surface we have to check that shifts 
   * don't lead out of the cell's node range */
  for (k=0;k<gridbegin;k+=1) {
    for (i=0;i<n_veloc;i++) {
      next = k+next[i];
      if (next>=0) {
	lbfluid[next].n_new[i] = lbfluid[k].n[i];
      }
    }
  }
  for (k=gridbegin;k<(xyzcube+gridsurface-gridbegin);k+=1) {
    for (i=0;i<n_veloc;i++) {
      lbfluid[k+next[i]]n_new[i] = lbfluid[k].n[i] ;
    }
  }
  for (k=(xyzcube+gridsurface-gridbegin);k<(xyzcube+gridsurface);k+=1) {
    for (i=0;i<n_veloc;i++) {
      next = k+next[i] ;
      if (next<(xyzcube+gridsurface)) {
	lbfluid[next].n_new[i] = lbfluid[k].n[i] ;
      }
    }
  }
  
  /* swap the pointers to n and n_new */
  tmp = n ;
  n = n_new ;
  n_new = tmp ;

#endif

}

/*@}*/

/** Propagate the Lattice Boltzmann dynamics.
 * This function is called from the integrator. Since the time step
 * for the lattice dynamics can be coarser than the MD time step,
 * we monitor the time since the last lattice update.
 */
void lb_propagate() {

  fluidstep+=time_step ;

  if (fluidstep>=tau) {

    fluidstep=0.0 ;

    /* apply external forces */
    lb_external_forces();

    /* collision step */
    lb_calc_collisions() ;

    /* exchange halo regions */
    halo_communication(&update_halo_comm) ;

    /* streaming step */
    lb_propagate_n() ;

#ifdef CONSTRAINTS
    /* boundary conditions */
    lb_boundary_conditions();
#endif

    /* exchange halo regions */
    halo_communication(&update_halo_comm) ;

  }

}

/***********************************************************************/
/** \name Coupling part */
/***********************************************************************/
/*@{*/

/** Transfer a certain amount of momentum to a elementray cell of fluid.
 * 
 * Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
 *
 * @param momentum   Momentum to be transfered to the fluid (lattice
 *                   units) (Input).
 * @param node_index Indices of the sites of the elementary lattice
 *                   cell (Input).
 * @param delta      Weights for the assignment to the single lattice
 *                   sites (Input).
 * @param badrandoms Flag/Counter for the occurrence negative
 *                   populations (Output).
 */
MDINLINE void lb_transfer_momentum(const double momentum[3], const int node_index[8], const double delta[6], int *badrandoms) {

  int i,x,y,z,dir ;
  double delta_j[3] ;
  double *local_n, *n_new ;

  /* We don't need to save the local populations because 
   * we use a trick for their restoration:
   * We substract the old force from the new one,
   * hence the previous change in the local populations
   * is automatically revoked during the recalculation.
   * Note that this makes it necessary to actually apply 
   * all changes and forbids to return immediately when negative
   * populations occur.
   */

  for (x=0;x<2;x++) {
    for (y=0;y<2;y++) {
      for (z=0;z<2;z++) {

	n_new = local_n = lbfluid[node_index[z*4+y*2+x]].n;
  
	for (dir=0;dir<3;dir++) {
	  delta_j[dir] = delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*momentum[dir] ;
	}

#ifdef D3Q18

	n_new[0]  = local_n[0] + 1./6.*delta_j[0] ;
	n_new[1]  = local_n[1] - 1./6.*delta_j[0] ;
	n_new[2]  = local_n[2] + 1./6.*delta_j[1] ;
	n_new[3]  = local_n[3] - 1./6.*delta_j[1] ;
	n_new[4]  = local_n[4] + 1./6.*delta_j[2] ;
	n_new[5]  = local_n[5] - 1./6.*delta_j[2] ;
	n_new[6]  = local_n[6] + 1./12.*(delta_j[0]+delta_j[1]) ;
	n_new[7]  = local_n[7] - 1./12.*(delta_j[0]+delta_j[1]) ;
	n_new[8]  = local_n[8] + 1./12.*(delta_j[0]-delta_j[1]) ;
	n_new[9]  = local_n[9] - 1./12.*(delta_j[0]-delta_j[1]) ;
	n_new[10] = local_n[10] + 1./12.*(delta_j[0]+delta_j[2]) ;
	n_new[11] = local_n[11] - 1./12.*(delta_j[0]+delta_j[2]) ;
	n_new[12] = local_n[12] + 1./12.*(delta_j[0]-delta_j[2]) ;
	n_new[13] = local_n[13] - 1./12.*(delta_j[0]-delta_j[2]) ;
	n_new[14] = local_n[14] + 1./12.*(delta_j[1]+delta_j[2]) ;
	n_new[15] = local_n[15] - 1./12.*(delta_j[1]+delta_j[2]) ;
	n_new[16] = local_n[16] + 1./12.*(delta_j[1]-delta_j[2]) ;
	n_new[17] = local_n[17] - 1./12.*(delta_j[1]-delta_j[2]) ;

	for (i=0;i<18;i++) {
	  if (n_new[i]<0.0) {
	    (*badrandoms)++;
	    LB_TRACE(fprintf(stderr,"%d: (%d,%d,%d) negative population %d (badrandoms=%d)\n",this_node,x,y,z,i,*badrandoms));
	    /* DO NOT break and return immediately here! */
	  }
	}

#else

	for (i=0;i<n_veloc;i++) {
	  tmp = 0.0 ;
	  for (k=0;k<SPACE_DIM;k++) {
	    tmp += delta_j[k] * c_g_d[i][k] ;
	  }
	  n_new[i] = local_n[i] + coeff[i][1] * tmp ;
	}

	for (i=0;i<n_veloc;i++) {
	  if (n_new[i]<0.0) {
	    (*badrandoms)++ ;
	    LB_TRACE(fprintf(stderr,"%d: (%d,%d,%d) negative population %d (badrandoms=%d)\n",this_node,x,y,z,i,*badrandoms));
	    /* DO NOT break and return immediately here! */
	  }
	}

#endif

      }
    }
  }

}

/** Coupling of a particle to viscous fluid with Stokesian friction.
 * 
 * Section II.C. Ahlrichs and Duenweg, JCP 111(17):8225 (1999)
 *
 * @param p          The coupled particle (Input).
 * @param badrandoms Flag/Counter for occurence of negative
 *                   populations (Output).
 * @param p_is_ghost Flag indicating whether the particle is a ghost
 *                   particle. Ghost particles must not have the force
 *                   added since it is already included in the real
 *                   image. However, ghosts must be treated to
 *                   transfer momentum to sites on different processors.
 */
MDINLINE void lb_viscous_momentum_exchange(Particle *p, int *badrandoms, int p_is_ghost) {

  int x,y,z,dir ;
  int node_index[8] ;
  double delta[6] ;
  double *local_rho, *local_j, interpolated_u[3], delta_j[3];
#ifdef ADDITIONAL_CHECKS
  double old_rho[8];
#endif

  /* determine elementary lattice cell surrounding the particle 
     and the relative position of the particle in this cell */ 
  map_position_to_lattice(&lblattice,p->r.p,node_index,delta) ;

  ONEPART_TRACE(if(p->p.identity==check_id && !p_is_ghost) fprintf(stderr,"%d: OPT: LB delta=(%.3f,%.3f,%.3f,%.3f,%.3f,%.3f) pos=(%.3f,%.3f,%.3f)\n",this_node,delta[0],delta[1],delta[2],delta[3],delta[4],delta[5],p->r.p[0],p->r.p[1],p->r.p[2]));

  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  interpolated_u[0] = interpolated_u[1] = interpolated_u[2] = 0.0 ;
  for (x=0;x<2;x++) {
    for (y=0;y<2;y++) {
      for (z=0;z<2;z++) {

	local_rho = lbfluid[node_index[z*4+y*2+x]].rho;
#ifdef ADDITIONAL_CHECKS
	old_rho[z*4+y*2+x] = *local_rho;
#endif
	local_j = lbfluid[node_index[z*4+y*2+x]].j ;
	ONEPART_TRACE(if(p->p.identity==check_id && !p_is_ghost) fprintf(stderr,"%d: OPT: LB fluid (%d,%d,%d local_rho=%.3f local_j=(%.3e,%.3f,%.3f) weight=%.3f\n",this_node,x,y,z,*local_rho,local_j[0],local_j[1],local_j[2],delta[3*x+0]*delta[3*y+1]*delta[3*z+2]));

	for (dir=0;dir<3;dir++) {
	  interpolated_u[dir] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[dir]/(*local_rho) ;
	}

      }
    }
  }
  
  ONEPART_TRACE(if(p->p.identity==check_id && !p_is_ghost) fprintf(stderr,"%d: OPT: LB u = (%.3e,%.3e,%.3e) v = (%.16e,%.3e,%.3e)\n",this_node,interpolated_u[0],interpolated_u[1],interpolated_u[2],p->m.v[0],p->m.v[1],p->m.v[2]));

  /* calculate viscous force and add to random force
     (take care to rescale the velocities with the time_step
     and transform fluid velocity to MD units) 
     (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  for (dir=0;dir<3;dir++) {
      p->t.f_random[dir] += - lbpar.friction * (p->m.v[dir]/time_step - interpolated_u[dir]*agrid/tau);
  }

  /* exchange momentum */
  if (transfer_momentum) {

    /* add force to particle if not ghost */
    if (!p_is_ghost) {
      p->f.f[0] += p->t.f_random[0];
      p->f.f[1] += p->t.f_random[1];
      p->f.f[2] += p->t.f_random[2];

      ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f = (%.3e,%.3e,%.3e)\n",this_node,p->t.f_random[0],p->t.f_random[1],p->t.f_random[2]));
      ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f = (%.9e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
    }

    /* transform momentum transfer to lattice units
       (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
    for (dir=0;dir<3;dir++) {
      delta_j[dir] = - p->t.f_random[dir]*time_step*tau/agrid;
    }

    lb_transfer_momentum(delta_j,node_index,delta,badrandoms);

    ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB delta_j = (%.3e,%.3e,%.3e)\n",this_node,delta_j[0],delta_j[1],delta_j[2]));
  }

#ifdef ADDITIONAL_CHECKS
  int i;
  for (i=0;i<8;i++) {
    lb_calc_local_rho(&lbfluid[node_index[i]]);
    local_rho = lbfluid[node_index[i]].rho;
    if (fabs(*local_rho-old_rho[i]) > ROUND_ERROR_PREC) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt,"{108 Mass loss/gain %le in lb_viscous_momentum_exchange} ",*local_rho-old_rho[i]);
    }
  }
#endif

}

/** Calculate particle lattice interactions.
 * So far, only viscous coupling with Stokesian friction is
 * implemented.
 * Include all particle-lattice forces in this function.
 * The function is called from \ref force_calc.
 *
 * Parallelizing the fluid particle coupling is not straightforward
 * because drawing of random numbers makes the whole thing nonlocal.
 * One way to do it is to treat every particle only on one node, i.e.
 * the random numbers need not be communicated. The particles that are 
 * not fully inside the local lattice are taken into account via their
 * ghost images on the neighbouring nodes. But this requires that the 
 * correct values of the surrounding lattice nodes are available on 
 * the respective node, which means that we have to communicate the 
 * halo regions before treating the ghost particles. Moreover, after 
 * determining the ghost couplings, we have to communicate back the 
 * halo region such that all local lattice nodes have the correct values.
 * Thus two communication phases are involved which will most likely be 
 * the bottleneck of the computation.
 *
 * Another way of dealing with the particle lattice coupling is to 
 * treat a particle and all of it's images explicitly. This requires the
 * communication of the random numbers used in the calculation of the 
 * coupling force. The problem is now that, if random numbers have to 
 * be redrawn, we cannot efficiently determine which particles and which 
 * images have to be re-calculated. We therefore go back to the outset
 * and go through the whole system again until no failure occurs during
 * such a sweep. In the worst case, this is very inefficient because
 * many things are recalculated although they actually don't need.
 * But we can assume that this happens extremely rarely and then we have
 * on average only one communication phase for the random numbers, which
 * probably makes this method preferable compared to the above one.
 */
void calc_particle_lattice_ia() {

  int i, k, dir, c, np, badrandoms, allbadrandoms ;
  Cell *cell ;
  Particle *p ;

  for (k=0;k<lblattice.halo_grid_volume;k++) {
      lb_calc_local_fields(&lbfluid[k],0);
  }

  /* draw random numbers for local particles 
   * the old random numbers are subtracted in order to
   * remove the previously made update */
  for (c=0;c<local_cells.n;c++) {
    cell = local_cells.cell[c] ;
    p = cell->part ;
    np = cell->n ;
    for (i=0;i<np;i++) {
      double x[3];
      for (dir=0;dir<3;dir++) {
	x[dir] = d_random()-0.5;
	p[i].t.f_random[dir] = -lb_coupl_pref*x[dir] ;
      }
      ONEPART_TRACE(if (p[i].p.identity==check_id) fprintf(stderr, "%d: OPT: LB f_random = (%.3e,%.3e,%.3e) (%.3e,%.3e,%.3e)\n",this_node,p[i].t.f_random[0],p[i].t.f_random[1],p[i].t.f_random[2],2.*x[0],2.*x[1],2.*x[2]));
    }
  }

  /* try random numbers until no failure occurs during a whole sweep */
  do {

    allbadrandoms = 0;
    badrandoms = 0;

    /* communicate the random numbers */
    ghost_communicator(&cell_structure.ghost_temp_comm) ;
    
    /* local cells */
    for (c=0;c<local_cells.n;c++) {
      cell = local_cells.cell[c] ;
      p = cell->part ;
      np = cell->n ;

      for (i=0;i<np;i++) {

	lb_viscous_momentum_exchange(&p[i],&badrandoms,0) ;

      }

    }

    /* ghost cells */
    for (c=0;c<ghost_cells.n;c++) {
      cell = ghost_cells.cell[c] ;
      p = cell->part ;
      np = cell->n ;

      for (i=0;i<np;i++) {

	  /* for ghost particles we have to check if they lie
	   * in the range of the local lattice nodes */
	if (p[i].r.p[0] >= my_left[0]-lblattice.agrid && p[i].r.p[0] < my_right[0]
	    && p[i].r.p[1] >= my_left[1]-lblattice.agrid && p[i].r.p[1] < my_right[1]
	    && p[i].r.p[2] >= my_left[2]-lblattice.agrid && p[i].r.p[2] < my_right[2]) {

	  ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: coupling of ghost\n",this_node));
	  lb_viscous_momentum_exchange(&p[i],&badrandoms,1) ;

	}
      }
    }

    MPI_Allreduce(&badrandoms, &allbadrandoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;

    if (allbadrandoms>0) {

        fprintf(stderr, "%d: Momentum error (badrandoms=%d, allbadrandoms=%d). Check your parameters if this happens too often!\n",this_node,badrandoms,allbadrandoms);

	for (c=0;c<local_cells.n;c++) {
	  cell = local_cells.cell[c] ;
	  p = cell->part ;
	  np = cell->n ;
	  for (i=0;i<np;i++) {
	    for (dir=0;dir<3;dir++) {
	      p[i].t.f_random[dir] = lb_coupl_pref*(d_random()-0.5) - p[i].t.f_random[dir] ;
	    }
	  }
	}

    }

  } while (allbadrandoms>0) ;

  halo_communication(&update_halo_comm) ;

}

/*@}*/

/***********************************************************************/
/** \name TCL stuff */
/***********************************************************************/

static int lb_parse_set_fields(Tcl_Interp *interp, int argc, char **argv, int *change, int *ind) {

  int k, node, index ;
  double rho, j[3] ;

  *change = 4 ;
  if (argc < 4) return TCL_ERROR ;
  if (!ARG0_IS_D(rho)) return TCL_ERROR ;
  for (k=0;k<3;k++) {
    if (!ARG_IS_D(k+1,j[k])) return TCL_ERROR ;
  }
    
  node = map_lattice_to_node(&lblattice,ind);
  index = get_linear_index(ind[0],ind[1],ind[2],lblattice.halo_grid);

  /* transform to lattice units */
  rho  *= agrid*agrid*agrid;
  j[0] *= tau/agrid;
  j[1] *= tau/agrid;
  j[2] *= tau/agrid;

  mpi_send_fluid(node,index,rho,j) ;

  lb_calc_average_rho();
  lb_reinit_parameters();

  return TCL_OK ;

}

static int lb_print_local_fields(Tcl_Interp *interp, int argc, char **argv, int *change, int *ind) {

  char buffer[256+4*TCL_DOUBLE_SPACE+3*TCL_INTEGER_SPACE];
  int node, index;
  double rho, j[3];

  *change = 0;

  sprintf(buffer, "%d", ind[0]) ;
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  sprintf(buffer, "%d", ind[1]) ;
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  sprintf(buffer, "%d", ind[2]) ;
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  node = map_lattice_to_node(&lblattice,ind);
  index = get_linear_index(ind[0],ind[1],ind[2],lblattice.halo_grid);
  
  mpi_recv_fluid(node,index,&rho,j) ;

  /* transform to MD units */
  rho  *= 1./(agrid*agrid*agrid);
  j[0] *= agrid/tau;
  j[1] *= agrid/tau;
  j[2] *= agrid/tau;

  Tcl_PrintDouble(interp, rho, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  Tcl_PrintDouble(interp, j[0], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  Tcl_PrintDouble(interp, j[1], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  Tcl_PrintDouble(interp, j[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
    
  return TCL_OK ;

}

static int lbfluid_parse_tau(Tcl_Interp *interp, int argc, char *argv[], int *change) {
    double tau;

    if (argc < 1) {
	Tcl_AppendResult(interp, "tau requires 1 argument", NULL);
	return TCL_ERROR;
    }
    if (!ARG0_IS_D(tau)) {
	Tcl_AppendResult(interp, "wrong  argument for tau", (char *)NULL);
	return TCL_ERROR;
    }
    if (tau < 0.0) {
	Tcl_AppendResult(interp, "tau must be positive", (char *)NULL);
	return TCL_ERROR;
    }
    else if ((time_step >= 0.0) && (tau < time_step)) {
      Tcl_AppendResult(interp, "tau must be larger than MD time_step", (char *)NULL);
      return TCL_ERROR;
    }

    *change = 1;
    lbpar.tau = tau;

    mpi_bcast_lb_params(LBPAR_TAU);
 
    return TCL_OK;
}

static int lbfluid_parse_agrid(Tcl_Interp *interp, int argc, char *argv[], int *change) {
    double agrid;

    if (argc < 1) {
	Tcl_AppendResult(interp, "agrid requires 1 argument", (char *)NULL);
	return TCL_ERROR;
    }
    if (!ARG0_IS_D(agrid)) {
	Tcl_AppendResult(interp, "wrong argument for agrid", (char *)NULL);
	return TCL_ERROR;
    }
    if (agrid <= 0.0) {
	Tcl_AppendResult(interp, "agrid must be positive", (char *)NULL);
	return TCL_ERROR;
    }

    *change = 1;
    lbpar.agrid = agrid;

    mpi_bcast_lb_params(LBPAR_AGRID);
 
    return TCL_OK;
}

static int lbfluid_parse_density(Tcl_Interp *interp, int argc, char *argv[], int *change) {
    double density;

    if (argc < 1) {
	Tcl_AppendResult(interp, "density requires 1 argument", (char *)NULL);
	return TCL_ERROR;
    }
    if (!ARG0_IS_D(density)) {
	Tcl_AppendResult(interp, "wrong argument for density", (char *)NULL);
	return TCL_ERROR;
    }
    if (density <= 0.0) {
	Tcl_AppendResult(interp, "density must be positive", (char *)NULL);
	return TCL_ERROR;
    }

    *change = 1;
    lbpar.rho = density;

    mpi_bcast_lb_params(LBPAR_DENSITY);
 
    return TCL_OK;
}

static int lbfluid_parse_viscosity(Tcl_Interp *interp, int argc, char *argv[], int *change) {
    double viscosity;

    if (argc < 1) {
	Tcl_AppendResult(interp, "viscosity requires 1 argument", (char *)NULL);
	return TCL_ERROR;
    }
    if (!ARG0_IS_D(viscosity)) {
	Tcl_AppendResult(interp, "wrong argument for viscosity", (char *)NULL);
	return TCL_ERROR;
    }
    if (viscosity <= 0.0) {
	Tcl_AppendResult(interp, "viscosity must be positive", (char *)NULL);
	return TCL_ERROR;
    }

    *change = 1;
    lbpar.viscosity = viscosity;

    mpi_bcast_lb_params(LBPAR_VISCOSITY);
 
    return TCL_OK;
}

static int lbfluid_parse_friction(Tcl_Interp *interp, int argc, char *argv[], int *change) {
    double friction;

    if (argc < 1) {
	Tcl_AppendResult(interp, "friction requires 1 argument", (char *)NULL);
	return TCL_ERROR;
    }
    if (!ARG0_IS_D(friction)) {
	Tcl_AppendResult(interp, "wrong argument for friction", (char *)NULL);
	return TCL_ERROR;
    }
    if (friction <= 0.0) {
	Tcl_AppendResult(interp, "friction must be positive", (char *)NULL);
	return TCL_ERROR;
    }

    *change = 1;
    lbpar.friction = friction;

    mpi_bcast_lb_params(LBPAR_FRICTION);
 
    return TCL_OK;
}

static int lbfluid_parse_ext_force(Tcl_Interp *interp, int argc, char *argv[], int *change) {
    double ext_f[3];
    if (argc < 3) {
	Tcl_AppendResult(interp, "ext_force requires 3 arguments", (char *)NULL);
	return TCL_ERROR;
    }
    else {
 	if (!ARG_IS_D(0, ext_f[0])) return TCL_ERROR;
	if (!ARG_IS_D(1, ext_f[1])) return TCL_ERROR;
	if (!ARG_IS_D(2, ext_f[2])) return TCL_ERROR;
    }
    
    *change = 3;
    lbpar.ext_force[0] = ext_f[0];
    lbpar.ext_force[1] = ext_f[1];
    lbpar.ext_force[2] = ext_f[2];
    
    mpi_bcast_lb_params(LBPAR_EXTFORCE);
 
    return TCL_OK;
}

/** Parser for the \ref lbfluid command. */
int lbfluid_cmd(ClientData data, Tcl_Interp *interp, int argc, char **argv) {

  int err = TCL_OK;
  int change = 0;
  
  argc--; argv++;

  if (argc < 1) {
      Tcl_AppendResult(interp, "too few arguments to \"lbfluid\"", (char *)NULL);
      err = TCL_ERROR;
  }
  else if (ARG0_IS_S("off")) {
    err = TCL_ERROR;
  }
  else if (ARG0_IS_S("init")) {
    err = TCL_ERROR;
  }
  else while (argc > 0) {
      if (ARG0_IS_S("agrid"))
	  err = lbfluid_parse_agrid(interp, argc-1, argv+1, &change);
      else if (ARG0_IS_S("tau"))
	  err = lbfluid_parse_tau(interp, argc-1, argv+1, &change);
      else if (ARG0_IS_S("density"))
	  err = lbfluid_parse_density(interp, argc-1, argv+1, &change);
      else if (ARG0_IS_S("viscosity"))
	  err = lbfluid_parse_viscosity(interp, argc-1, argv+1, &change);
      else if (ARG0_IS_S("friction"))
	  err = lbfluid_parse_friction(interp, argc-1, argv+1, &change);
      else if (ARG0_IS_S("ext_force"))
	  err = lbfluid_parse_ext_force(interp, argc-1, argv+1, &change);
      else {
	  Tcl_AppendResult(interp, "unknown feature \"", argv[1],"\" of lbfluid", (char *)NULL);
	  err = TCL_ERROR ;
      }

      if ((err = mpi_gather_runtime_errors(interp, err))) break;
      
      argc -= (change + 1);
      argv += (change + 1);
  }

  lattice_switch = (lattice_switch | LATTICE_LB) ;
  mpi_bcast_parameter(FIELD_LATTICE_SWITCH) ;

  /* thermo_switch is retained for backwards compatibility */
  thermo_switch = (thermo_switch | THERMO_LB);
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);

  return err;    
}

/*@}*/

#endif /* LB */
