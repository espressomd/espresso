/* $Id: lbgpu_cfile.c $
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
 * Copyright (c) 2002-2007; all rights reserved unless otherwise stated.
 */

/** \file lbgpu_cfile.c
 *
 * C file for the Lattice Boltzmann implementation on GPUs.
 * Header file for \ref lbgpu.h.
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
#include "lb-boundaries.h"
#ifdef LB_GPU

/** Action number for \ref mpi_get_particles. */
#define REQ_GETPARTS  16
#ifndef D3Q19
#error The implementation only works for D3Q19 so far!
#endif

/** Struct holding the Lattice Boltzmann parameters */
LB_parameters_gpu lbpar_gpu = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0 ,0.0, -1.0, 0, 0, 0, 0, 0, 0, 1, 0, {0.0, 0.0, 0.0}, 12345, 0};
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

/** (re-) allocation of the memory need for the particles (cpu part)*/
void lb_realloc_particles_gpu(){

  lbpar_gpu.number_of_particles = n_total_particles;
  LB_TRACE (printf("#particles realloc\t %u \n", lbpar_gpu.number_of_particles));
  //fprintf(stderr, "%u \t \n", lbpar_gpu.number_of_particles);
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

  LB_TRACE(printf("Begin initialzing fluid on GPU\n"));
  /** set parameters for transfer to gpu */
  lb_reinit_parameters_gpu();

  lb_realloc_particles_gpu();
	
  lb_init_GPU(&lbpar_gpu);

  LB_TRACE(printf("Initialzing fluid on GPU successful\n"));
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
  LB_particle_force_gpu *host_forces_sl;
  Cell *cell;
  int c, i;
  MPI_Status status;

  n_part = cells_get_n_particles();

  COMM_TRACE(fprintf(stderr, "%d: send_particles_slave, %d particles\n", this_node, n_part));


  if (n_part > 0) {
    int g = 0;
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
#endif/* LB_GPU */
