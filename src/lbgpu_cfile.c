/*
  Copyright (C) 2010,2011 The ESPResSo project

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
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "lbgpu.h"
#include "utils.h"
#include "parser.h"
#include "communication.h"
#include "thermostat.h"
#include "grid.h"
#include "domain_decomposition.h"
#include "integrate.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "global.h"
#include "lb_boundaries_gpu.h"

#ifdef LB_GPU

/** Action number for \ref mpi_get_particles. */
#define REQ_GETPARTS  16
#ifndef D3Q19
#error The implementation only works for D3Q19 so far!
#endif

/** Struct holding the Lattice Boltzmann parameters */
LB_parameters_gpu lbpar_gpu = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0 ,0.0, -1.0, 0, 0, 0, 0, 0, 0, 1, 0, {0.0, 0.0, 0.0}, 12345, 0.0, 0};

LB_values_gpu *host_values = NULL;
LB_nodes_gpu *host_nodes = NULL;
LB_particle_force_gpu *host_forces = NULL;
LB_particle_gpu *host_data = NULL;
/** Flag indicating momentum exchange between particles and fluid */
int transfer_momentum_gpu = 0;

static int max_ran = 1000000;
/*@}*/
//static double tau;

/** measures the MD time since the last fluid update */
static int fluidstep = 0;

/** c_sound_square in LB units*/
static float c_sound_sq = 1.f/3.f;

//clock_t start, end;
int i;

static FILE *datei;
//static char file[300];
static void mpi_get_particles_lb(LB_particle_gpu *host_result);
static void mpi_get_particles_slave_lb();
static void mpi_send_forces_lb(LB_particle_force_gpu *host_forces);
static void mpi_send_forces_slave_lb();

int n_extern_nodeforces = 0;
LB_extern_nodeforce_gpu *host_extern_nodeforces = NULL;

/*-----------------------------------------------------------*/
/** main of lb_gpu_programm */
/*-----------------------------------------------------------*/
/** lattice boltzmann update gpu called from integrate.c
*/
void lattice_boltzmann_update_gpu() {

  int factor = (int)round(lbpar_gpu.tau/time_step);

  fluidstep += 1;

  if (fluidstep>=factor) {
    fluidstep=0;

    lb_integrate_GPU();

    LB_TRACE (fprintf(stderr,"lb_integrate_GPU \n"));

  }
}

/** Calculate particle lattice interactions called from forces.c
*/
void lb_calc_particle_lattice_ia_gpu() {

  if (transfer_momentum_gpu) {
    mpi_get_particles_lb(host_data);

    if(this_node == 0){
#if 0
      for (i=0;i<n_total_particles;i++) {
        fprintf(stderr, "%i particle posi: , %f %f %f\n", i, host_data[i].p[0], host_data[i].p[1], host_data[i].p[2]);
      }
#endif

    if(lbpar_gpu.number_of_particles) lb_particle_GPU(host_data);

    LB_TRACE (fprintf(stderr,"lb_calc_particle_lattice_ia_gpu \n"));

    }
  }
}

/**copy forces from gpu to cpu and call mpi routines to add forces to particles
*/
void lb_send_forces_gpu(){

  if (transfer_momentum_gpu) {
    if(this_node == 0){
      if (lbpar_gpu.number_of_particles) lb_copy_forces_GPU(host_forces);

      LB_TRACE (fprintf(stderr,"lb_send_forces_gpu \n"));
#if 0
        for (i=0;i<n_total_particles;i++) {
          fprintf(stderr, "%i particle forces , %f %f %f \n", i, host_forces[i].f[0], host_forces[i].f[1], host_forces[i].f[2]);
        }
#endif
    }
    mpi_send_forces_lb(host_forces);
  }
}
/** allocation of the needed memory for phys. values and particle data residing in the cpu memory
*/
void lb_pre_init_gpu() {
	 
  lbpar_gpu.number_of_particles = 0;

  LB_TRACE (fprintf(stderr,"#nodes \t %u \n", lbpar_gpu.number_of_nodes));

  /*-----------------------------------------------------*/
  /* allocating of the needed memory for several structs */
  /*-----------------------------------------------------*/
  
  /* Struct holding calc phys values rho, j, phi of every node */
  size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_values_gpu);
  host_values = (LB_values_gpu*)malloc(size_of_values);

  LB_TRACE (fprintf(stderr,"lb_pre_init_gpu \n"));
}

/** (re-)allocation of the memory needed for the phys. values and if needed memory for the nodes located in the cpu memory
*/ 
static void lb_realloc_fluid_gpu() {
	 
  LB_TRACE (printf("#nodes \t %u \n", lbpar_gpu.number_of_nodes));

  /**-----------------------------------------------------*/
  /** reallocating of the needed memory for several structs */
  /**-----------------------------------------------------*/

  /**Struct holding calc phys values rho, j, phi of every node*/
  size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_values_gpu);
  host_values = realloc(host_values, size_of_values);

  LB_TRACE (fprintf(stderr,"lb_realloc_fluid_gpu \n"));
}
/** (re-) allocation of the memory need for the particles (cpu part)*/
void lb_realloc_particles_gpu(){

  lbpar_gpu.number_of_particles = n_total_particles;
  LB_TRACE (printf("#particles realloc\t %u \n", lbpar_gpu.number_of_particles));
  /**-----------------------------------------------------*/
  /** allocating of the needed memory for several structs */
  /**-----------------------------------------------------*/
  /**Allocate struct for particle forces */
  size_t size_of_forces = lbpar_gpu.number_of_particles * sizeof(LB_particle_force_gpu);
  host_forces = realloc(host_forces, size_of_forces);

  lbpar_gpu.your_seed = (unsigned int)i_random(max_ran);

  LB_TRACE (fprintf(stderr,"test your_seed %u \n", lbpar_gpu.your_seed));
  lb_realloc_particle_GPU(&lbpar_gpu, &host_data);
}
/** (Re-)initializes the fluid according to the given value of rho. */
void lb_reinit_fluid_gpu() {

  //lbpar_gpu.your_seed = (unsigned int)i_random(max_ran);
  lb_reinit_parameters_gpu();
  if(lbpar_gpu.number_of_nodes != 0){
    lb_reinit_GPU(&lbpar_gpu);
    lbpar_gpu.reinit = 1;
  }

  LB_TRACE (fprintf(stderr,"lb_reinit_fluid_gpu \n"));
}

/** Release the fluid. */
/*not needed in Espresso but still not deleted*/
void lb_release_gpu(){

  free(host_nodes);
  free(host_values);
  free(host_forces);
  free(host_data);
}
/** (Re-)initializes the fluid. */
void lb_reinit_parameters_gpu() {

  lbpar_gpu.mu = 0.0;
  lbpar_gpu.time_step = (float)time_step;

  if (lbpar_gpu.viscosity > 0.0) {
    /* Eq. (80) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    lbpar_gpu.gamma_shear = 1. - 2./(6.*lbpar_gpu.viscosity*lbpar_gpu.tau/(lbpar_gpu.agrid*lbpar_gpu.agrid) + 1.);   
  }

  if (lbpar_gpu.bulk_viscosity > 0.0) {
    /* Eq. (81) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    lbpar_gpu.gamma_bulk = 1. - 2./(9.*lbpar_gpu.bulk_viscosity*lbpar_gpu.tau/(lbpar_gpu.agrid*lbpar_gpu.agrid) + 1.);
  }

  if (temperature > 0.0) {  /* fluctuating hydrodynamics ? */

    lbpar_gpu.fluct = 1;
	LB_TRACE (fprintf(stderr, "fluct on \n"));
    /* Eq. (51) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007).*/
    /* Note that the modes are not normalized as in the paper here! */

    lbpar_gpu.mu = (float)temperature/c_sound_sq*lbpar_gpu.tau*lbpar_gpu.tau/(lbpar_gpu.agrid*lbpar_gpu.agrid);
    //lbpar_gpu->mu *= agrid*agrid*agrid;  // Marcello's conjecture

    /* lb_coupl_pref is stored in MD units (force)
     * Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
     * The factor 12 comes from the fact that we use random numbers
     * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
     * time_step comes from the discretization.
     */

    lbpar_gpu.lb_coupl_pref = sqrt(12.f*2.f*lbpar_gpu.friction*(float)temperature/lbpar_gpu.time_step);
    lbpar_gpu.lb_coupl_pref2 = sqrt(2.f*lbpar_gpu.friction*(float)temperature/lbpar_gpu.time_step);

  } else {
    /* no fluctuations at zero temperature */
    lbpar_gpu.fluct = 0;
    lbpar_gpu.lb_coupl_pref = 0.0;
    lbpar_gpu.lb_coupl_pref2 = 0.0;
  }
	LB_TRACE (fprintf(stderr,"lb_reinit_prarameters_gpu \n"));

  reinit_parameters_GPU(&lbpar_gpu);
}

/** Performs a full initialization of
 *  the Lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values. */
void lb_init_gpu() {

  //lb_free_GPU();
  /** set parameters for transfer to gpu */
  lb_reinit_parameters_gpu();

  lb_realloc_particles_gpu();
	
  lb_realloc_fluid_gpu();

  lb_init_GPU(&lbpar_gpu);

  LB_TRACE (fprintf(stderr,"lb_init_gpu \n"));
}

/*@}*/

/***********************************************************************/
/** \name MPI stuff */
/***********************************************************************/

/*************** REQ_GETPARTS ************/
static void mpi_get_particles_lb(LB_particle_gpu *host_data)
{
  int n_part;
  int g, pnode;
  Cell *cell;
  int c;
  MPI_Status status;

  int i;	
  int *sizes;
  sizes = malloc(sizeof(int)*n_nodes);

  n_part = cells_get_n_particles();

  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* just check if the number of particles is correct */
  if(this_node > 0){
    /* call slave functions to provide the slave datas */
    mpi_get_particles_slave_lb();
  }
  else {
    /* master: fetch particle informations into 'result' */
    g = 0;
    for (pnode = 0; pnode < n_nodes; pnode++) {
      if (sizes[pnode] > 0) {
        if (pnode == 0) {
          for (c = 0; c < local_cells.n; c++) {
            Particle *part;
            int npart;	
            int dummy[3] = {0,0,0};
            double pos[3];
            cell = local_cells.cell[c];
            part = cell->part;
            npart = cell->n;
            for (i=0;i<npart;i++) {
              memcpy(pos, part[i].r.p, 3*sizeof(double));
              fold_position(pos, dummy);
              host_data[i+g].p[0] = (float)pos[0];
              host_data[i+g].p[1] = (float)pos[1];
              host_data[i+g].p[2] = (float)pos[2];
								
              host_data[i+g].v[0] = (float)part[i].m.v[0];
              host_data[i+g].v[1] = (float)part[i].m.v[1];
              host_data[i+g].v[2] = (float)part[i].m.v[2];
#ifdef LB_ELECTROHYDRODYNAMICS
              host_data[i+g].mu_E[0] = (float)part[i].p.mu_E[0];
              host_data[i+g].mu_E[1] = (float)part[i].p.mu_E[1];
              host_data[i+g].mu_E[2] = (float)part[i].p.mu_E[2];
#endif
            }  
            g += npart;
          }  
        }
        else {
          MPI_Recv(&host_data[g], sizes[pnode]*sizeof(LB_particle_gpu), MPI_BYTE, pnode, REQ_GETPARTS,
          MPI_COMM_WORLD, &status);
          g += sizes[pnode];
        }
      }
    }
  }
  COMM_TRACE(fprintf(stderr, "%d: finished get\n", this_node));
  free(sizes);
}

static void mpi_get_particles_slave_lb(){
 
  int n_part;
  int g;
  LB_particle_gpu *host_data_sl;
  Cell *cell;
  int c, i;

  n_part = cells_get_n_particles();

  COMM_TRACE(fprintf(stderr, "%d: get_particles_slave, %d particles\n", this_node, n_part));

  if (n_part > 0) {
    /* get (unsorted) particle informations as an array of type 'particle' */
    /* then get the particle information */
    host_data_sl = malloc(n_part*sizeof(LB_particle_gpu));
    
    g = 0;
    for (c = 0; c < local_cells.n; c++) {
      Particle *part;
      int npart;
      int dummy[3] = {0,0,0};
      double pos[3];
      cell = local_cells.cell[c];
      part = cell->part;
      npart = cell->n;

      for (i=0;i<npart;i++) {
        memcpy(pos, part[i].r.p, 3*sizeof(double));
        fold_position(pos, dummy);	
			
        host_data_sl[i+g].p[0] = (float)pos[0];
        host_data_sl[i+g].p[1] = (float)pos[1];
        host_data_sl[i+g].p[2] = (float)pos[2];

        host_data_sl[i+g].v[0] = (float)part[i].m.v[0];
        host_data_sl[i+g].v[1] = (float)part[i].m.v[1];
        host_data_sl[i+g].v[2] = (float)part[i].m.v[2];
#ifdef LB_ELECTROHYDRODYNAMICS
        host_data_sl[i+g].mu_E[0] = (float)part[i].p.mu_E[0];
        host_data_sl[i+g].mu_E[1] = (float)part[i].p.mu_E[1];
        host_data_sl[i+g].mu_E[2] = (float)part[i].p.mu_E[2];
#endif
      }
      g+=npart;
    }
    /* and send it back to the master node */
    MPI_Send(host_data_sl, n_part*sizeof(LB_particle_gpu), MPI_BYTE, 0, REQ_GETPARTS, MPI_COMM_WORLD);
    free(host_data_sl);
  }  
}

static void mpi_send_forces_lb(LB_particle_force_gpu *host_forces){
	
  int n_part;
  int g, pnode;
  Cell *cell;
  int c;
  int i;	
  int *sizes;
  sizes = malloc(sizeof(int)*n_nodes);
  n_part = cells_get_n_particles();
  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* call slave functions to provide the slave datas */
  if(this_node > 0) {
    mpi_send_forces_slave_lb();
  }
  else{
  /* fetch particle informations into 'result' */
  g = 0;
    for (pnode = 0; pnode < n_nodes; pnode++) {
      if (sizes[pnode] > 0) {
        if (pnode == 0) {
          for (c = 0; c < local_cells.n; c++) {
            int npart;	
            cell = local_cells.cell[c];
            npart = cell->n;
            for (i=0;i<npart;i++) {
              cell->part[i].f.f[0] += (double)host_forces[i+g].f[0];
              cell->part[i].f.f[1] += (double)host_forces[i+g].f[1];
              cell->part[i].f.f[2] += (double)host_forces[i+g].f[2];
            }
 	    g += npart;
          }
        }
        else {
        /* and send it back to the slave node */
        MPI_Send(&host_forces[g], sizes[pnode]*sizeof(LB_particle_force_gpu), MPI_BYTE, pnode, REQ_GETPARTS, MPI_COMM_WORLD);			
        g += sizes[pnode];
        }
      }
    }
  }
  COMM_TRACE(fprintf(stderr, "%d: finished send\n", this_node));

  free(sizes);
}

static void mpi_send_forces_slave_lb(){

  int n_part;
  int g;
  LB_particle_force_gpu *host_forces_sl;
  Cell *cell;
  int c, i;
  MPI_Status status;

  n_part = cells_get_n_particles();

  COMM_TRACE(fprintf(stderr, "%d: send_particles_slave, %d particles\n", this_node, n_part));


  if (n_part > 0) {
    /* get (unsorted) particle informations as an array of type 'particle' */
    /* then get the particle information */
    host_forces_sl = malloc(n_part*sizeof(LB_particle_force_gpu));
    MPI_Recv(host_forces_sl, n_part*sizeof(LB_particle_force_gpu), MPI_BYTE, 0, REQ_GETPARTS,
    MPI_COMM_WORLD, &status);
    for (c = 0; c < local_cells.n; c++) {
      int npart;	
      cell = local_cells.cell[c];
      npart = cell->n;
      for (i=0;i<npart;i++) {
        cell->part[i].f.f[0] += (double)host_forces_sl[i+g].f[0];
        cell->part[i].f.f[1] += (double)host_forces_sl[i+g].f[1];
        cell->part[i].f.f[2] += (double)host_forces_sl[i+g].f[2];
      }
      g += npart;
    }
    free(host_forces_sl);
  } 
}
/*@}*/

/***********************************************************************/
/** \name TCL stuff */
/***********************************************************************/


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
    fprintf(stderr,"tau %f \n", lbpar_gpu.tau);
    fprintf(stderr,"time_step %f \n", time_step);
    Tcl_AppendResult(interp, "tau must be larger than MD time_step", (char *)NULL);
    return TCL_ERROR;
  }

  *change = 1;
  lbpar_gpu.tau = (float)tau;

  on_lb_params_change_gpu(0);

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
  lbpar_gpu.agrid = (float)agrid;

  lbpar_gpu.dim_x = (unsigned int)floor(box_l[0]/agrid);
  lbpar_gpu.dim_y = (unsigned int)floor(box_l[1]/agrid);
  lbpar_gpu.dim_z = (unsigned int)floor(box_l[2]/agrid);

  unsigned int tmp[3];
  tmp[0] = lbpar_gpu.dim_x;
  tmp[1] = lbpar_gpu.dim_y;
  tmp[2] = lbpar_gpu.dim_z;
  /* sanity checks */
  int dir;
  for (dir=0;dir<3;dir++) {
  /* check if box_l is compatible with lattice spacing */
    if (fabs(box_l[dir]-tmp[dir]*agrid) > ROUND_ERROR_PREC) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{097 Lattice spacing agrid=%f is incompatible with box_l[%i]=%f} ",agrid,dir,box_l[dir]);
    }
  }

  lbpar_gpu.number_of_nodes = lbpar_gpu.dim_x * lbpar_gpu.dim_y * lbpar_gpu.dim_z;

  on_lb_params_change_gpu(LBPAR_AGRID);

  LB_TRACE (printf("#nodes \t %u \n", lbpar_gpu.number_of_nodes));
 
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
  lbpar_gpu.rho = (float)density;

  on_lb_params_change_gpu(LBPAR_DENSITY);

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
  lbpar_gpu.viscosity = (float)viscosity;

  on_lb_params_change_gpu(LBPAR_VISCOSITY);
 
  return TCL_OK;
}

static int lbfluid_parse_bulk_visc(Tcl_Interp *interp, int argc, char *argv[], int *change) {
  double bulk_visc;

  if (argc < 1) {
    Tcl_AppendResult(interp, "bulk_viscosity requires 1 argument", (char *)NULL);
    return TCL_ERROR;
  }
  if (!ARG0_IS_D(bulk_visc)) {
    Tcl_AppendResult(interp, "wrong argument for bulk_viscosity", (char *)NULL);
    return TCL_ERROR;
  }
  if (bulk_visc < 0.0) {
    Tcl_AppendResult(interp, "bulk_viscosity must be positive", (char *)NULL);
    return TCL_ERROR;
  }

  *change =1;
  lbpar_gpu.bulk_viscosity = (float)bulk_visc;

  on_lb_params_change_gpu(LBPAR_BULKVISC);

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
  lbpar_gpu.friction = (float)friction;

  on_lb_params_change_gpu(LBPAR_FRICTION);

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

  /* external force density is stored in MD units */
  lbpar_gpu.ext_force[0] = (float)ext_f[0];
  lbpar_gpu.ext_force[1] = (float)ext_f[1];
  lbpar_gpu.ext_force[2] = (float)ext_f[2];

  lbpar_gpu.external_force = 1;

  lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);
    
  return TCL_OK;
}

static int lbfluid_parse_gamma_odd(Tcl_Interp *interp, int argc, char *argv[], int *change) {

  double g;

  if (argc < 1) {
    Tcl_AppendResult(interp, "gamma_odd requires 1 argument", (char *)NULL);
    return TCL_ERROR;
  }
  if (!ARG0_IS_D(g)) {
    Tcl_AppendResult(interp, "wrong argument for gamma_odd", (char *)NULL);
    return TCL_ERROR;
  }
  if (fabs( g > 1.0)) {
    Tcl_AppendResult(interp, "fabs(gamma_odd) must be > 1.", (char *)NULL);
    return TCL_ERROR;
  }

  *change = 1;
  lbpar_gpu.gamma_odd = (float)g;

  on_lb_params_change_gpu(0);

  return TCL_OK;
}

static int lbfluid_parse_gamma_even(Tcl_Interp *interp, int argc, char *argv[], int *change) {

  double g;

  if (argc < 1) {
    Tcl_AppendResult(interp, "gamma_even requires 1 argument", (char *)NULL);
    return TCL_ERROR;
  }
  if (!ARG0_IS_D(g)) {
    Tcl_AppendResult(interp, "wrong argument for gamma_even", (char *)NULL);
    return TCL_ERROR;
  }
  if (fabs( g > 1.0)) {
    Tcl_AppendResult(interp, "fabs(gamma_even) must be > 1.", (char *)NULL);
    return TCL_ERROR;
  }

  *change = 1;
  lbpar_gpu.gamma_even = (float)g;

  on_lb_params_change_gpu(0);

  return TCL_OK;
}

static int lbprint_parse_velocity(Tcl_Interp *interp, int argc, char *argv[], int *change, int vtk){

  if (argc < 1) {
    Tcl_AppendResult(interp, "file requires at least 1 argument", (char *)NULL);
    return TCL_ERROR;
  }

  *change = 1;

  datei=fopen(argv[0],"w");
  if(datei == NULL){
    fprintf(stderr, "couldn't open datafile! \n");
    exit(1);
  }
  lb_get_values_GPU(host_values);
  if(vtk == 1){
    fprintf(datei, "# vtk DataFile Version 2.0\ntest\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %u %u %u\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA %u\nSCALARS OutArray  floats 3\nLOOKUP_TABLE default\n", lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z, lbpar_gpu.number_of_nodes);
  }	
  int j;	
  for(j=0; j<lbpar_gpu.number_of_nodes; ++j){
    /** print of the calculated phys values */
    fprintf(datei, " %f \t %f \t %f \n", host_values[j].v[0], host_values[j].v[1], host_values[j].v[2]);
  }
  fclose(datei);
  return TCL_OK;
}
static int lbprint_parse_density(Tcl_Interp *interp, int argc, char *argv[], int *change) {

  if (argc < 1) {
    Tcl_AppendResult(interp, "file requires at least 1 argument", (char *)NULL);
    return TCL_ERROR;
  }

  *change = 1;

  datei=fopen(argv[0],"w");
  if(datei == NULL){
    fprintf(stderr, "couldn't open datafile! \n");
    exit(1);
  }
  lb_get_values_GPU(host_values);
  int j;	
  for(j=0; j<lbpar_gpu.number_of_nodes; ++j){
    /** print of the calculated phys values */
    fprintf(datei, " %f \n", host_values[j].rho);
  }

  return TCL_OK;
}
#if 0
static int lbprint_parse_stresstensor(Tcl_Interp *interp, int argc, char *argv[], int *change) {

    if (argc < 1) {
	Tcl_AppendResult(interp, "file requires at least 1 argument", (char *)NULL);
	return TCL_ERROR;
    }

    *change = 1;

	datei=fopen(argv[0],"w");
		if(datei == NULL){
			fprintf(stderr, "couldn't open datafile! \n");
			exit(1);
		}
	lb_get_values_GPU(host_values);
	int j;	
	for(j=0; j<lbpar_gpu.number_of_nodes; ++j){
	/** print of the calculated phys values */
		fprintf(datei, " %f \t %f \t %f \t %f \t %f \t %f \n", host_values[j].pi[0], host_values[j].pi[1], host_values[j].pi[2],
 															   host_values[j].pi[3], host_values[j].pi[4], host_values[j].pi[5]);
	}

    return TCL_OK;

}
#endif
#endif /* LB_GPU */

#ifdef LB_GPU
static int lbnode_parse_set(Tcl_Interp *interp, int argc, char **argv, int *ind) {
  unsigned int index;
  double f[3];
  size_t size_of_extforces;
  int change = 0;

  if ( ind[0] >=  lbpar_gpu.dim_x ||  ind[1] >= lbpar_gpu.dim_y ||  ind[2] >= lbpar_gpu.dim_z ) {
    Tcl_AppendResult(interp, "position is not in the LB lattice", (char *)NULL);
    return TCL_ERROR;
  }

  index = ind[0] + ind[1]*lbpar_gpu.dim_x + ind[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;
  while (argc > 0) {
    if (change==1) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "Error in lbnode_extforce force. You can only change one field at the same time.", (char *)NULL);
      return TCL_ERROR;
    }
  if(ARG0_IS_S("force")){
    if (ARG1_IS_D(f[0])) {
      argc--;
      argv++;
    } else return TCL_ERROR;
    if (ARG1_IS_D(f[1])) { 
      argc--;
      argv++;
    } else return TCL_ERROR;
    if (ARG1_IS_D(f[2])) {
      argc--;
      argv++;
    } else return TCL_ERROR;
    change=1;
  }
  size_of_extforces = (n_extern_nodeforces+1)*sizeof(LB_extern_nodeforce_gpu);
  host_extern_nodeforces = realloc(host_extern_nodeforces, size_of_extforces);
 
  host_extern_nodeforces[n_extern_nodeforces].force[0] = (float)f[0];
  host_extern_nodeforces[n_extern_nodeforces].force[1] = (float)f[1];
  host_extern_nodeforces[n_extern_nodeforces].force[2] = (float)f[2];

  host_extern_nodeforces[n_extern_nodeforces].index = index;
  n_extern_nodeforces++;
  
  if(lbpar_gpu.external_force == 0)lbpar_gpu.external_force = 1;

  --argc; ++argv;

  lb_init_extern_nodeforces_GPU(n_extern_nodeforces, host_extern_nodeforces, &lbpar_gpu);
  }

  return TCL_OK;
}
#endif /* LB_GPU */
/** Parser for the lbnode TCL command. 
*/
int tclcommand_lbnode_gpu(Tcl_Interp *interp, int argc, char **argv) {
#ifdef LB_GPU

  int coord[3];
  int counter;
  char double_buffer[TCL_DOUBLE_SPACE];
  LB_values_gpu *host_print_values;
  host_print_values = malloc(sizeof(LB_values_gpu));	
  int single_nodeindex;
  --argc; ++argv;
  if (argc < 3) {
    Tcl_AppendResult(interp, "too few arguments for lbnode", (char *)NULL);
    return TCL_ERROR;
  }

  if (!ARG_IS_I(0,coord[0]) || !ARG_IS_I(1,coord[1]) || !ARG_IS_I(2,coord[2])) {
    Tcl_AppendResult(interp, "wrong arguments for lbnode", (char *)NULL);
    return TCL_ERROR;
  } 
  argc-=3; argv+=3;
   
  if (argc == 0 ) { 
    Tcl_AppendResult(interp, "lbnode syntax: lbnode X Y Z [ print ] [ rho | u ]", (char *)NULL);
    return TCL_ERROR;
  }
  single_nodeindex = coord[0] + coord[1]*lbpar_gpu.dim_x + coord[2]*lbpar_gpu.dim_x*lbpar_gpu.dim_y;

  if (ARG0_IS_S("print")) {
    argc--; argv++;
    if (argc == 0 ) { 
      Tcl_AppendResult(interp, "lbnode syntax: lbnode X Y Z [ print ] [ rho | u ]", (char *)NULL);
      return TCL_ERROR;
    }
    while (argc > 0) {
      if (ARG0_IS_S("rho") || ARG0_IS_S("density")) {

      lb_print_node_GPU(single_nodeindex, host_print_values);
      Tcl_PrintDouble(interp, (double)host_print_values[0].rho, double_buffer);
      Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
      argc--; argv++;
      }
      else if (ARG0_IS_S("u") || ARG0_IS_S("v") || ARG0_IS_S("velocity")) { 
        lb_print_node_GPU(single_nodeindex, host_print_values);
        for (counter = 0; counter < 3; counter++) {
          Tcl_PrintDouble(interp, (double)host_print_values[0].v[counter], double_buffer);
          Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        }
        argc--; argv++;
      }
      else if (ARG0_IS_S("ux") || ARG0_IS_S("vx")) {
        lb_print_node_GPU(single_nodeindex, host_print_values);
        Tcl_PrintDouble(interp, (double)host_print_values[0].v[0], double_buffer);
        Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        argc--; argv++;
      }
      else if (ARG0_IS_S("uy") || ARG0_IS_S("vy")) {
        lb_print_node_GPU(single_nodeindex, host_print_values);
        Tcl_PrintDouble(interp, (double)host_print_values[0].v[1], double_buffer);
        Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        argc--; argv++;
      }
      else if (ARG0_IS_S("uz") || ARG0_IS_S("vz")) {
        lb_print_node_GPU(single_nodeindex, host_print_values);
        Tcl_PrintDouble(interp, (double)host_print_values[0].v[2], double_buffer);
        Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
        argc--; argv++;
      }
      else {
        Tcl_ResetResult(interp);
        Tcl_AppendResult(interp, "unknown fluid data \"", argv[0], "\" requested", (char *)NULL);
        return TCL_ERROR;
      }
    }
  }
  else {
    Tcl_AppendResult(interp, "unknown feature \"", argv[0], "\" of lbnode", (char *)NULL);
    return  TCL_ERROR;
  }     
  return TCL_OK;

#else /* !defined LB_GPU */
  Tcl_AppendResult(interp, "LB_GPU is not compiled in!", NULL);
  return TCL_ERROR;
#endif /* LB_GPU */
}
#ifdef LB_GPU
/** Parser for the \ref tclcommand_lbnode_extforce_gpu command. Can be used in future to set more values like rho,u e.g.
*/
int tclcommand_lbnode_extforce_gpu(ClientData data, Tcl_Interp *interp, int argc, char **argv) {

  int err=TCL_ERROR;
  int coord[3];

  --argc; ++argv;
  
  if (argc < 3) {
    Tcl_AppendResult(interp, "too few arguments for lbnode_extforce", (char *)NULL);
    return TCL_ERROR;
  }

  if (!ARG_IS_I(0,coord[0]) || !ARG_IS_I(1,coord[1]) || !ARG_IS_I(2,coord[2])) {
    Tcl_AppendResult(interp, "wrong arguments for lbnode", (char *)NULL);
    return TCL_ERROR;
  } 
  argc-=3; argv+=3;

  if (argc == 0 ) { 
    Tcl_AppendResult(interp, "lbnode_extforce syntax: lbnode_extforce X Y Z [ print | set ] [ F(X) | F(Y) | F(Z) ]", (char *)NULL);
    return TCL_ERROR;
  }

  if (ARG0_IS_S("set")) 
    err = lbnode_parse_set(interp, argc-1, argv+1, coord);
    else {
    Tcl_AppendResult(interp, "unknown feature \"", argv[0], "\" of lbnode_extforce", (char *)NULL);
    return  TCL_ERROR;
    }     
  return err;
}
#endif /* LB_GPU */

/** Parser for the \ref lbfluid command gpu.
*/
int tclcommand_lbfluid_gpu(Tcl_Interp *interp, int argc, char **argv) {
#ifdef LB_GPU
  int err = TCL_OK;
  int change = 0;

  while (argc > 0) {
    if (ARG0_IS_S("grid") || ARG0_IS_S("agrid"))
      err = lbfluid_parse_agrid(interp, argc-1, argv+1, &change);
    else if (ARG0_IS_S("tau"))
      err = lbfluid_parse_tau(interp, argc-1, argv+1, &change);
    else if (ARG0_IS_S("density") || ARG0_IS_S("dens"))
      err = lbfluid_parse_density(interp, argc-1, argv+1, &change);
    else if (ARG0_IS_S("viscosity") || ARG0_IS_S("visc"))
      err = lbfluid_parse_viscosity(interp, argc-1, argv+1, &change);
    else if (ARG0_IS_S("bulk_viscosity") || ARG0_IS_S("b_visc"))
      err = lbfluid_parse_bulk_visc(interp, argc-1, argv+1, &change);
    else if (ARG0_IS_S("friction") || ARG0_IS_S("coupling"))
      err = lbfluid_parse_friction(interp, argc-1, argv+1, &change);
    else if (ARG0_IS_S("ext_force"))
      err = lbfluid_parse_ext_force(interp, argc-1, argv+1, &change);
    else if (ARG0_IS_S("gamma_odd"))
      err = lbfluid_parse_gamma_odd(interp, argc-1, argv+1, &change);
    else if (ARG0_IS_S("gamma_even"))
      err = lbfluid_parse_gamma_even(interp, argc-1, argv+1, &change);
    else {
      Tcl_AppendResult(interp, "unknown feature \"", argv[0],"\" of lbfluid", (char *)NULL);
      err = TCL_ERROR ;
    }
    if (err == TCL_ERROR) return TCL_ERROR;
      argc -= (change + 1);
      argv += (change + 1);
  }

  mpi_bcast_parameter(FIELD_LATTICE_SWITCH) ;

  /* thermo_switch is retained for backwards compatibility */
  thermo_switch = (thermo_switch | THERMO_LB);
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
	
  LB_TRACE (fprintf(stderr,"tclcommand_lbfluid_gpu parser ok \n"));

  return err;
#else /* !defined LB_GPU */
  Tcl_AppendResult(interp, "LB_GPU is not compiled in!", NULL);
  return TCL_ERROR;
#endif
}

#ifdef LB_GPU
/** printing the hole fluid field to file with order x+y*dim_x+z*dim_x*dim_y  */
int tclcommand_lbprint_gpu(ClientData data, Tcl_Interp *interp, int argc, char **argv) {

  int err = TCL_OK;
  int change = 0;
  int vtk = 0;

  argc--; argv++;

  if (argc < 1) {
    Tcl_AppendResult(interp, "too few arguments to \"lbprint\"", (char *)NULL);
    err = TCL_ERROR;
  }
  else while (argc > 0) {
    if (ARG0_IS_S("u") || ARG0_IS_S("velocity") || ARG0_IS_S("v")){
      argc--; argv++;
      if (ARG0_IS_S("vtk")){
        vtk = 1;
	err = lbprint_parse_velocity(interp, argc-1, argv+1, &change, vtk); 
      }
      else
	err = lbprint_parse_velocity(interp, argc, argv, &change, vtk);
    }
    else if (ARG0_IS_S("rho") || ARG0_IS_S("density"))
      err = lbprint_parse_density(interp, argc-1, argv+1, &change);   
    else if (ARG0_IS_S("stresstensor")){
      //err = lbprint_parse_stresstensor(interp, argc-1, argv+1, &change); 
      Tcl_AppendResult(interp, "\"lbprint stresstensor\" is not available by default due to memory saving, pls ensure availablity of pi[6] (see lbgpu.h) and lbprint_parse_stresstensor()", (char *)NULL);
    err = TCL_ERROR;
    }  
    else {
      Tcl_AppendResult(interp, "unknown feature \"", argv[0],"\" of lbprint", (char *)NULL);
      err = TCL_ERROR ;
    }
    argc -= (change + 1);
    argv += (change + 1);

    LB_TRACE (fprintf(stderr,"tclcommand_lbprint_gpu parser ok \n"));
  }
  return err;    
}
#endif/* LB_GPU */
