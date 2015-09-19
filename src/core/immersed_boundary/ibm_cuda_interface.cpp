

// *******
// This is an internal file of the IMMERSED BOUNDARY implementation
// It should not be included by any main Espresso routines
// Functions to be exported for Espresso are in ibm_main.hpp

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include <mpi.h>
#include "cells.hpp"
#include "grid.hpp"
#include "communication.hpp"
#include "particle_data.hpp"
#include "integrate.hpp"
#include "immersed_boundary/ibm_cuda_interface.hpp"

/// MPI tags for sending velocities and receiving particles
#define REQ_CUDAIBMSENDVEL 0xcc03
#define REQ_CUDAIBMGETPART 0xcc04

// Variables for communication
IBM_CUDA_ParticleDataInput *IBM_ParticleDataInput_host = NULL;
IBM_CUDA_ParticleDataOutput *IBM_ParticleDataOutput_host = NULL;

/*****************
   IBM_cuda_mpi_get_particles
Gather particle positions on the master node in order to communicate them to GPU
We transfer all particles (real and virtual), but acutally we would only need the virtual ones
Room for improvement...
 *****************/

void IBM_cuda_mpi_get_particles_slave();

// Analogous to the usual cuda_mpi_get_particles function
void IBM_cuda_mpi_get_particles()
{
  int g, pnode;
  Cell *cell;
  int c;
  MPI_Status status;
  
  int *sizes = (int*) Utils::malloc(sizeof(int)*n_nodes);
  int n_part = cells_get_n_particles();
  
  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);
  
  if(this_node > 0){
    /* call slave functions to provide the slave datas */
    IBM_cuda_mpi_get_particles_slave();
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
            for (int i=0;i<npart;i++) {
              memmove(pos, part[i].r.p, 3*sizeof(double));
              fold_position(pos, dummy);
              IBM_ParticleDataInput_host[i+g].pos[0] = (float)pos[0];
              IBM_ParticleDataInput_host[i+g].pos[1] = (float)pos[1];
              IBM_ParticleDataInput_host[i+g].pos[2] = (float)pos[2];
              
              IBM_ParticleDataInput_host[i+g].f[0] = (float)part[i].f.f[0];
              IBM_ParticleDataInput_host[i+g].f[1] = (float)part[i].f.f[1];
              IBM_ParticleDataInput_host[i+g].f[2] = (float)part[i].f.f[2];
              
              IBM_ParticleDataInput_host[i+g].isVirtual = part[i].p.isVirtual;
              
            }
            g += npart;
          }
        }
        else {
          MPI_Recv(&IBM_ParticleDataInput_host[g], sizes[pnode]*sizeof(IBM_CUDA_ParticleDataInput), MPI_BYTE, pnode, REQ_CUDAIBMGETPART,
                   comm_cart, &status);
          g += sizes[pnode];
        }
      }
    }
  }
  COMM_TRACE(fprintf(stderr, "%d: finished get\n", this_node));
  free(sizes);
}

void IBM_cuda_mpi_get_particles_slave()
{
  
  int n_part;
  int g;
  IBM_CUDA_ParticleDataInput *particle_input_sl;
  Cell *cell;
  int c, i;
  
  n_part = cells_get_n_particles();
  
  COMM_TRACE(fprintf(stderr, "%d: get_particles_slave, %d particles\n", this_node, n_part));
  
  if (n_part > 0) {
    /* get (unsorted) particle informations as an array of type 'particle' */
    /* then get the particle information */
    //        particle_data_host_sl = (IBM_CUDA_ParticleDataInput*) Utils::malloc(n_part*sizeof(IBM_CUDA_ParticleData));
    particle_input_sl = new IBM_CUDA_ParticleDataInput[n_part];
    
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
        memmove(pos, part[i].r.p, 3*sizeof(double));
        fold_position(pos, dummy);
        
        particle_input_sl[i+g].pos[0] = (float)pos[0];
        particle_input_sl[i+g].pos[1] = (float)pos[1];
        particle_input_sl[i+g].pos[2] = (float)pos[2];
        
        particle_input_sl[i+g].f[0] = (float)part[i].f.f[0];
        particle_input_sl[i+g].f[1] = (float)part[i].f.f[1];
        particle_input_sl[i+g].f[2] = (float)part[i].f.f[2];
        
        particle_input_sl[i+g].isVirtual = part[i].p.isVirtual;
      }
      g+=npart;
    }
    /* and send it back to the master node */
    MPI_Send(particle_input_sl, n_part*sizeof(IBM_CUDA_ParticleDataInput), MPI_BYTE, 0, REQ_CUDAIBMGETPART, comm_cart);
    delete []particle_input_sl;
  }
}

/*****************
   IBM_cuda_mpi_send_velocities
Particle velocities have been communicated from GPU, now transmit to all nodes
 ******************/
// Analogous to cuda_mpi_send_forces

void IBM_cuda_mpi_send_velocities_slave();

void IBM_cuda_mpi_send_velocities()
{
  int n_part;
  int g, pnode;
  Cell *cell;
  int c;
  int *sizes;
  sizes = (int *) Utils::malloc(sizeof(int)*n_nodes);
  n_part = cells_get_n_particles();
  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);
  
  /* call slave functions to provide the slave datas */
  if(this_node > 0) {
    IBM_cuda_mpi_send_velocities_slave();
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
            for (int i=0;i<npart;i++)
            {
              if ( cell->part[i].p.isVirtual )
              {
                cell->part[i].m.v[0] = (double)IBM_ParticleDataOutput_host[(i+g)].v[0];
                cell->part[i].m.v[1] = (double)IBM_ParticleDataOutput_host[(i+g)].v[1];
                cell->part[i].m.v[2] = (double)IBM_ParticleDataOutput_host[(i+g)].v[2];
              }
            }
            g += npart;
          }
        }
        else {
          /* and send it back to the slave node */
          MPI_Send(&IBM_ParticleDataOutput_host[g], sizes[pnode]*sizeof(IBM_CUDA_ParticleDataOutput), MPI_BYTE, pnode, REQ_CUDAIBMSENDVEL, comm_cart);
          g += sizes[pnode];
        }
      }
    }
  }
  COMM_TRACE(fprintf(stderr, "%d: finished send\n", this_node));
  
  free(sizes);
}

void IBM_cuda_mpi_send_velocities_slave()
{
  
  Cell *cell;
  int c, i;
  MPI_Status status;
  
  const int n_part = cells_get_n_particles();
  
  COMM_TRACE(fprintf(stderr, "%d: send_particles_slave, %d particles\n", this_node, n_part));
  
  
  if (n_part > 0) {
    int g = 0;
    IBM_CUDA_ParticleDataOutput *output_sl = new IBM_CUDA_ParticleDataOutput[n_part];
    MPI_Recv(output_sl, n_part*sizeof(IBM_CUDA_ParticleDataOutput), MPI_BYTE, 0, REQ_CUDAIBMSENDVEL,
             comm_cart, &status);
    for (c = 0; c < local_cells.n; c++) {
      int npart;
      cell = local_cells.cell[c];
      npart = cell->n;
      for (i=0;i<npart;i++)
      {
        if ( cell->part[i].p.isVirtual )
        {
          cell->part[i].m.v[0] = (double)output_sl[(i+g)].v[0];
          cell->part[i].m.v[1] = (double)output_sl[(i+g)].v[1];
          cell->part[i].m.v[2] = (double)output_sl[(i+g)].v[2];
        }
      }
      g += npart;
    }
    delete []output_sl;
  }
}


#endif