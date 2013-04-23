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
  static GPU_global_part_vars global_part_vars = {0,0,0};
  static __device__ __constant__ GPU_global_part_vars global_part_vars_gpu;
  
  /** struct for particle force */
  static LB_particle_force_gpu *particle_force = NULL;
  /** struct for particle position and veloctiy */
  static LB_particle_gpu *particle_data = NULL;
  /** struct for storing particle rn seed */
  static LB_particle_seed_gpu *part = NULL;

  LB_particle_gpu *host_data = NULL;
  LB_particle_force_gpu *host_forces = NULL;

  /**cuda streams for parallel computing on cpu and gpu */
  cudaStream_t stream[1];

  cudaError_t err;
  cudaError_t _err;
  
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

  /**setup and call particle reallocation from the host
   * @param *lbpar_gpu	Pointer to parameters to setup the lb field
   * @param **host_data	Pointer to host information data
  */
  void gpu_init_particle_comm() {
    
    //we only run the function if there are new particles which have been created since the last call of this function
    if ( global_part_vars.number_of_particles != n_total_particles ) {
      
      global_part_vars.seed = (unsigned int)i_random(max_ran);
      global_part_vars.number_of_particles = n_total_particles;
      global_part_vars.communication_enabled = 0; //assume no particles to begin with

      cuda_safe_mem(cudaMemcpyToSymbol(global_part_vars_gpu, &global_part_vars, sizeof(GPU_global_part_vars)));

      if ( host_forces )    cudaFreeHost(host_forces); //if the arrays exists free them to prevent memory leaks
      if ( host_data )      cudaFreeHost(host_data);
      if ( particle_force ) cudaFree(particle_force);
      if ( particle_data )  cudaFree(particle_data);
      if ( part )           cudaFree(part);

      if ( global_part_vars.number_of_particles ) {

        global_part_vars.communication_enabled = 1;
        
    #if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200
        /**pinned memory mode - use special function to get OS-pinned memory*/
        cudaHostAlloc((void**)&host_data, global_part_vars.number_of_particles * sizeof(LB_particle_gpu), cudaHostAllocWriteCombined);
        cudaHostAlloc((void**)&host_forces, global_part_vars.number_of_particles * sizeof(LB_particle_force_gpu), cudaHostAllocWriteCombined);
    #else
        cudaMallocHost((void**)&host_data, global_part_vars.number_of_particles * sizeof(LB_particle_gpu));
        cudaMallocHost((void**)&host_forces, global_part_vars.number_of_particles * sizeof(LB_particle_force_gpu));
    #endif

        cuda_safe_mem(cudaMalloc((void**)&particle_force, global_part_vars.number_of_particles * sizeof(LB_particle_force_gpu)));
        cuda_safe_mem(cudaMalloc((void**)&particle_data, global_part_vars.number_of_particles * sizeof(LB_particle_gpu)));
        cuda_safe_mem(cudaMalloc((void**)&part, global_part_vars.number_of_particles * sizeof(LB_particle_seed_gpu)));
        
        /** values for the particle kernel */
        int threads_per_block_particles = 64;
        int blocks_per_grid_particles_y = 4;
        int blocks_per_grid_particles_x = (global_part_vars.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1)/(threads_per_block_particles * blocks_per_grid_particles_y);
        dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

        KERNELCALL(init_particle_force, dim_grid_particles, threads_per_block_particles, (particle_force, part));
      }
    }
  }

  LB_particle_gpu* gpu_get_particle_pointer() {
    return particle_data;
  }
  
  GPU_global_part_vars* gpu_get_global_particle_vars_pointer() {
    return &global_part_vars;
  }
  LB_particle_force_gpu* gpu_get_particle_force_pointer() {
    return particle_force;
  }

  LB_particle_seed_gpu* gpu_get_particle_seed_pointer() {
    return part;
  }

  void copy_part_data_to_gpu() {

    if ( global_part_vars.communication_enabled == 1 ) {
      
      mpi_get_particles_lb(host_data);

      /** get espresso md particle values*/
      cudaMemcpyAsync(particle_data, host_data, global_part_vars.number_of_particles * sizeof(LB_particle_gpu), cudaMemcpyHostToDevice, stream[0]);

    }
  }



  /** setup and call kernel to copy particle forces to host
   * @param *host_forces contains the particle force computed on the GPU
  */
  void copy_forces_from_GPU() {

    if ( global_part_vars.communication_enabled == 1 ) {

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
      mpi_send_forces_lb(host_forces);
    }
  }

  
}
