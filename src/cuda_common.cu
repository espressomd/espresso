/* 
   Copyright (C) 2010,2011,2012,2013 The ESPResSo project

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


#include "cuda_common.h"
extern "C" {

#include "config.h"
#include "random.h"
#include "particle_data.h"

  static void mpi_get_particles_lb(LB_particle_gpu *host_result);
  static void mpi_get_particles_slave_lb();
  static void mpi_send_forces_lb(LB_particle_force_gpu *host_forces);
  static void mpi_send_forces_slave_lb();
  
  static int max_ran = 1000000;
  static GPU_global_part_vars global_part_vars = {0,0,0};
  static __device__ __constant__ GPU_global_part_vars global_part_vars_gpu;
  
  /** struct for particle force */
  static LB_particle_force_gpu *particle_force = NULL;
  /** struct for particle position and veloctiy */
  static LB_particle_gpu *particle_data = NULL;
  /** struct for storing particle rn seed */
  static LB_particle_seed_gpu *part = NULL;

  LB_particle_gpu *host_data = NULL;
  
  /**cuda streams for parallel computing on cpu and gpu */
  extern cudaStream_t stream[1];

  extern cudaError_t err;
  extern cudaError_t _err;
  
}


__device__ unsigned int getThreadIndex() {

  return blockIdx.y * gridDim.x * blockDim.x +
         blockDim.x * blockIdx.x +
         threadIdx.x;
}

/** kernel for the initalisation of the particle force array
 * @param *particle_force	Pointer to local particle force (Output)
 * @param *part			Pointer to the particle rn seed storearray (Output)
*/
__global__ void init_particle_force(LB_particle_force_gpu *particle_force, LB_particle_seed_gpu *part){

  unsigned int part_index = getThreadIndex();

  if(part_index<global_part_vars_gpu.number_of_particles){
    particle_force[part_index].f[0] = 0.0f;
    particle_force[part_index].f[1] = 0.0f;
    particle_force[part_index].f[2] = 0.0f;

    part[part_index].seed = global_part_vars_gpu.seed + part_index;
  }

}


/** kernel for the initalisation of the partikel force array
 * @param *particle_force	pointer to local particle force (Input)
*/
__global__ void reset_particle_force(LB_particle_force_gpu *particle_force){
	
  unsigned int part_index = getThreadIndex();
	
  if(part_index<global_part_vars_gpu.number_of_particles){
    particle_force[part_index].f[0] = 0.0f;
    particle_force[part_index].f[1] = 0.0f;
    particle_force[part_index].f[2] = 0.0f;
  }			
}


extern "C" {

  void cuda_enable_particle_communication() {
    global_part_vars.number_of_particles = n_total_particles;
    //TODO
  }

  /**setup and call particle reallocation from the host
   * @param *lbpar_gpu	Pointer to parameters to setup the lb field
   * @param **host_data	Pointer to host information data
  */
  void gpu_init_particle_comm( LB_particle_gpu **host_data ) {
    
    //we only run the function if there are new particles which have been created since the last call of this function
    if ( global_part_vars.number_of_particles == n_total_particles ) {
      
      global_part_vars.seed = (unsigned int)i_random(max_ran);
      global_part_vars.number_of_particles = n_total_particles;
      global_part_vars.communication_enabled = 1;

      cuda_safe_mem(cudaMemcpyToSymbol(global_part_vars_gpu, &global_part_vars, sizeof(GPU_global_part_vars)));

      if ( host_data )      cudaFreeHost(*host_data);
      if ( particle_force ) cudaFree(particle_force);
      if ( particle_data )  cudaFree(particle_data);
      if ( part )           cudaFree(part);

    #if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200
      /**pinned memory mode - use special function to get OS-pinned memory*/
      cudaHostAlloc((void**)host_data, global_part_vars.number_of_particles * sizeof(LB_particle_gpu), cudaHostAllocWriteCombined);
    #else
      cudaMallocHost((void**)host_data, global_part_vars.number_of_particles * sizeof(LB_particle_gpu));
    #endif

      cuda_safe_mem(cudaMalloc((void**)&particle_force, global_part_vars.number_of_particles * sizeof(LB_particle_force_gpu)));
      cuda_safe_mem(cudaMalloc((void**)&particle_data, global_part_vars.number_of_particles * sizeof(LB_particle_gpu)));
      cuda_safe_mem(cudaMalloc((void**)&part, global_part_vars.number_of_particles * sizeof(LB_particle_seed_gpu)));
      
      /** values for the particle kernel */
      int threads_per_block_particles = 64;
      int blocks_per_grid_particles_y = 4;
      int blocks_per_grid_particles_x = (global_part_vars.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1)/(threads_per_block_particles * blocks_per_grid_particles_y);
      dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

      if ( global_part_vars.number_of_particles )
        KERNELCALL(init_particle_force, dim_grid_particles, threads_per_block_particles, (particle_force, part));

    }
  }

  void lb_get_particle_pointer(LB_particle_gpu** pointeradress) {
    *pointeradress = particle_data;
  }

  void lb_get_particle_force_pointer(LB_particle_force_gpu** pointeradress) {
    *pointeradress = particle_force;
  }

  void copy_part_data_to_gpu() {

    mpi_get_particles_lb(host_data);

    /** get espresso md particle values*/
    cudaMemcpyAsync(particle_data, host_data, global_part_vars.number_of_particles * sizeof(LB_particle_gpu), cudaMemcpyHostToDevice, stream[0]);
  }


  /** setup and call kernel to copy particle forces to host
   * @param *host_forces contains the particle force computed on the GPU
  */
  void lb_copy_forces_GPU(LB_particle_force_gpu *host_forces){

    /** Copy result from device memory to host memory*/
    cudaMemcpy(host_forces, particle_force, global_part_vars.number_of_particles * sizeof(LB_particle_force_gpu), cudaMemcpyDeviceToHost);

      /** values for the particle kernel */
    int threads_per_block_particles = 64;
    int blocks_per_grid_particles_y = 4;
    int blocks_per_grid_particles_x = (global_part_vars.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1)/(threads_per_block_particles * blocks_per_grid_particles_y);
    dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

    /** reset part forces with zero*/
    KERNELCALL(reset_particle_force, dim_grid_particles, threads_per_block_particles, (particle_force));
	
    cudaThreadSynchronize();
  }


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
    MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);

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

  #ifdef ELECTROSTATICS
                if (coulomb.method == COULOMB_P3M_GPU) {
                  host_data[i+g].q = (float)part[i].p.q;
                }
  #endif
              }  
              g += npart;
            }  
          }
          else {
            MPI_Recv(&host_data[g], sizes[pnode]*sizeof(LB_particle_gpu), MPI_BYTE, pnode, REQ_GETPARTS,
            comm_cart, &status);
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

  #ifdef ELECTROSTATICS
          if (coulomb.method == COULOMB_P3M_GPU) {
            host_data_sl[i+g].q = (float)part[i].p.q;
          }
  #endif
        }
        g+=npart;
      }
      /* and send it back to the master node */
      MPI_Send(host_data_sl, n_part*sizeof(LB_particle_gpu), MPI_BYTE, 0, REQ_GETPARTS, comm_cart);
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
    MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);

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
          MPI_Send(&host_forces[g], sizes[pnode]*sizeof(LB_particle_force_gpu), MPI_BYTE, pnode, REQ_GETPARTS, comm_cart);			
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
      comm_cart, &status);
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
  
}
