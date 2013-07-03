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


#include "cuda_common.h" //I can't go in extern C


#include "grid.h"
extern "C" {

#include "config.h"
#include "random.h"
#include "particle_data.h"
#include "interaction_data.h"


  
  static int max_ran = 1000000;
  static CUDA_global_part_vars global_part_vars_host = {0,0,0};
  static __device__ __constant__ CUDA_global_part_vars global_part_vars_device;
  
  /** struct for particle force */
  static CUDA_particle_force *particle_forces_device = NULL;
  /** struct for particle position and veloctiy */
  static CUDA_particle_data *particle_data_device = NULL;
  /** struct for storing particle rn seed */
  static CUDA_particle_seed *particle_seeds_device = NULL;

  CUDA_particle_data *particle_data_host = NULL;
  CUDA_particle_force *particle_forces_host = NULL;

  /**cuda streams for parallel computing on cpu and gpu */
  cudaStream_t stream[1];

  cudaError_t err;
  cudaError_t _err;
  
}

void _cuda_safe_mem(cudaError_t err, char *file, unsigned int line){
  if( cudaSuccess != err) {                                             
    fprintf(stderr, "Cuda Memory error at %s:%u.\n", file, line);
    printf("CUDA error: %s\n", cudaGetErrorString(err));
    if ( err == cudaErrorInvalidValue )
      fprintf(stderr, "You may have tried to allocate zero memory at %s:%u.\n", file, line);
    exit(EXIT_FAILURE);
  } else {
    _err=cudaGetLastError();
    if (_err != cudaSuccess) {
      fprintf(stderr, "Error found during memory operation. Possibly however from an failed operation before. %s:%u.\n", file, line);
      printf("CUDA error: %s\n", cudaGetErrorString(err));
      if ( _err == cudaErrorInvalidValue )
	fprintf(stderr, "You may have tried to allocate zero memory before %s:%u.\n", file, line);
      exit(EXIT_FAILURE);
    }
  }
}


__device__ unsigned int getThreadIndex() {

  return blockIdx.y * gridDim.x * blockDim.x +
         blockDim.x * blockIdx.x +
         threadIdx.x;
}



/** kernel for the initalisation of the particle force array
 * @param *particle_forces_device	    Pointer to local particle force (Output)
 * @param *particle_seeds_device			Pointer to the particle rn seed storearray (Output)
*/
__global__ void init_particle_force(CUDA_particle_force *particle_forces_device, CUDA_particle_seed *particle_seeds_device){

  unsigned int part_index = getThreadIndex();

  if(part_index<global_part_vars_device.number_of_particles){
    particle_forces_device[part_index].f[0] = 0.0f;
    particle_forces_device[part_index].f[1] = 0.0f;
    particle_forces_device[part_index].f[2] = 0.0f;

    particle_seeds_device[part_index].seed = global_part_vars_device.seed + part_index;
  }

}


/** kernel for the initalisation of the partikel force array
 * @param *particle_forces_device	pointer to local particle force (Input)
*/
__global__ void reset_particle_force(CUDA_particle_force *particle_forces_device){
	
  unsigned int part_index = getThreadIndex();
	
  if(part_index<global_part_vars_device.number_of_particles){
    particle_forces_device[part_index].f[0] = 0.0f;
    particle_forces_device[part_index].f[1] = 0.0f;
    particle_forces_device[part_index].f[2] = 0.0f;
  }			
}


extern "C" {


  /** change number of particles to be communicated to the GPU
   *  Note that in addition to calling this function the parameters must be broadcast with either:
   * 1) cuda_bcast_global_part_params(); (when just being executed on the master node) or
   * 2) MPI_Bcast(gpu_get_global_particle_vars_pointer_host(), sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart); (when executed on all nodes)
   */
  void gpu_change_number_of_part_to_comm() {
    //we only run the function if there are new particles which have been created since the last call of this function

    if ( global_part_vars_host.number_of_particles != n_total_particles && global_part_vars_host.communication_enabled == 1 && this_node == 0) {
      
      global_part_vars_host.seed = (unsigned int)i_random(max_ran);
      global_part_vars_host.number_of_particles = n_total_particles;

      cuda_safe_mem(cudaMemcpyToSymbol(global_part_vars_device, &global_part_vars_host, sizeof(CUDA_global_part_vars)));

      if ( particle_forces_host )    cudaFreeHost(particle_forces_host); //if the arrays exists free them to prevent memory leaks
      if ( particle_data_host )      cudaFreeHost(particle_data_host);
      if ( particle_forces_device )  cudaFree(particle_forces_device);
      if ( particle_data_device )    cudaFree(particle_data_device);
      if ( particle_seeds_device )   cudaFree(particle_seeds_device);


      if ( global_part_vars_host.number_of_particles ) {

    #if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200
        /**pinned memory mode - use special function to get OS-pinned memory*/
        cudaHostAlloc((void**)&particle_data_host, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_data), cudaHostAllocWriteCombined);
        cudaHostAlloc((void**)&particle_forces_host, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_force), cudaHostAllocWriteCombined);
    #else
        cudaMallocHost((void**)&particle_data_host, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_data));
        cudaMallocHost((void**)&particle_forces_host, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_force));
    #endif

        cuda_safe_mem(cudaMalloc((void**)&particle_forces_device, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_force)));
        cuda_safe_mem(cudaMalloc((void**)&particle_data_device, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_data)));
        cuda_safe_mem(cudaMalloc((void**)&particle_seeds_device, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_seed)));
        
        /** values for the particle kernel */
        int threads_per_block_particles = 64;
        int blocks_per_grid_particles_y = 4;
        int blocks_per_grid_particles_x = (global_part_vars_host.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1)/(threads_per_block_particles * blocks_per_grid_particles_y);
        dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

        KERNELCALL(init_particle_force, dim_grid_particles, threads_per_block_particles, (particle_forces_device, particle_seeds_device));
      }

    }

  }

  /** setup and call particle reallocation from the host
   *  Note that in addition to calling this function the parameters must be broadcast with either:
   * 1) cuda_bcast_global_part_params(); (when just being executed on the master node) or
   * 2) MPI_Bcast(gpu_get_global_particle_vars_pointer_host(), sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart); (when executed on all nodes)
   */
  void gpu_init_particle_comm() {
    if ( this_node == 0  && global_part_vars_host.communication_enabled == 0 ) {
      if( cuda_get_n_gpus() == -1 ) {
        fprintf(stderr, "Unable to initialize CUDA as no sufficient GPU is available.\n");
        exit(0);
      }
      if (cuda_get_n_gpus()>1) {
        fprintf (stderr, "More than one GPU detected, please note Espresso uses device 0 by default regardless of usage or capability\n");
        fprintf (stderr, "Note that the GPU to be used can be modified using cuda setdevice <int>\n");
        if (cuda_check_gpu(0)!=ES_OK) {
          fprintf (stderr, "WARNING!  CUDA device 0 is not capable of running Espresso but is used by default.  Espresso has detected a CUDA capable card but it is not the one used by Espresso by default\n");
          fprintf (stderr, "Please set the GPU to use with the cuda setdevice <int> command.\n");
          fprintf (stderr, "A list of available GPUs can be accessed using cuda list.\n");
        }
      }
    }
    global_part_vars_host.communication_enabled = 1;

    gpu_change_number_of_part_to_comm();

  }

  CUDA_particle_data* gpu_get_particle_pointer() {
    return particle_data_device;
  }
  CUDA_global_part_vars* gpu_get_global_particle_vars_pointer_host() {
    return &global_part_vars_host;
  }  
  CUDA_global_part_vars* gpu_get_global_particle_vars_pointer() {
    return &global_part_vars_device;
  }
  CUDA_particle_force* gpu_get_particle_force_pointer() {
    return particle_forces_device;
  }

  CUDA_particle_seed* gpu_get_particle_seed_pointer() {
    return particle_seeds_device;
  }

  void copy_part_data_to_gpu() {

    if ( global_part_vars_host.communication_enabled == 1 && global_part_vars_host.number_of_particles ) {
     
      cuda_mpi_get_particles(particle_data_host);

      /** get espresso md particle values*/
      if ( this_node == 0 ) cudaMemcpyAsync(particle_data_device, particle_data_host, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_data), cudaMemcpyHostToDevice, stream[0]);

    }
  }



  /** setup and call kernel to copy particle forces to host
  */
  void copy_forces_from_GPU() {

    if ( global_part_vars_host.communication_enabled == 1 && global_part_vars_host.number_of_particles ) {

      /** Copy result from device memory to host memory*/
      if ( this_node == 0 ) {
        cuda_safe_mem (cudaMemcpy(particle_forces_host, particle_forces_device, global_part_vars_host.number_of_particles * sizeof(CUDA_particle_force), cudaMemcpyDeviceToHost));


        /** values for the particle kernel */
        int threads_per_block_particles = 64;
        int blocks_per_grid_particles_y = 4;
        int blocks_per_grid_particles_x = (global_part_vars_host.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1)/(threads_per_block_particles * blocks_per_grid_particles_y);
        dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

        /** reset part forces with zero*/

        KERNELCALL(reset_particle_force, dim_grid_particles, threads_per_block_particles, (particle_forces_device));
        cudaThreadSynchronize();
      }
      cuda_mpi_send_forces(particle_forces_host);
    }
  }

  
}
