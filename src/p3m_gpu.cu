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

/** \file lbgpu.cu
 *
 * Cuda (.cu) file for the Lattice Boltzmann implementation on GPUs.
 * Header file for \ref lbgpu.h.
 */

#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include "cuda_common.h"

extern "C" {
#include "p3m_gpu.h"
#include "lbgpu.h"
#include "config.h"

//extern cudaStream_t stream[1];
//extern cudaError_t _err;
}

__device__ unsigned int getThreadIndexP3M() { //rename is dumb but can't import same fnc from cuda_common

  return blockIdx.y * gridDim.x * blockDim.x +
         blockDim.x * blockIdx.x +
         threadIdx.x;
}



__global__ void add_p3m_farfield_force_gpu( LB_parameters_gpu* lb_parameters_gpu,
                                            LB_particle_gpu* lb_particle_gpu,
                                            LB_particle_force_gpu* lb_particle_force_gpu
                                          ) {

  unsigned int index = getThreadIndexP3M();

  if( index < lb_parameters_gpu->number_of_particles ) {
    
    lb_particle_force_gpu[ index ].f[0] = 1.0f;
    lb_particle_force_gpu[ index ].f[1] = 2.0f;
    lb_particle_force_gpu[ index ].f[2] = 3.0f;
  }
}


extern "C" {
  
void p3m_gpu_add_farfield_force() {
  
  LB_parameters_gpu* lb_parameters;
  LB_parameters_gpu* lb_parameters_gpu;
  LB_particle_gpu* lb_particle_gpu;
  LB_particle_force_gpu* lb_particle_force_gpu;
  
  lb_get_lbpar_pointer( &lb_parameters );
  lb_get_para_pointer( &lb_parameters_gpu );
  lb_particle_gpu = gpu_get_particle_pointer();
  lb_particle_force_gpu = gpu_get_particle_force_pointer();

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
    ( lb_parameters->number_of_particles + threads_per_block * blocks_per_grid_y - 1 ) /
    ( threads_per_block * blocks_per_grid_y );
  dim3 dim_grid = make_uint3( blocks_per_grid_x, blocks_per_grid_y, 1 );
  
  KERNELCALL( add_p3m_farfield_force_gpu, dim_grid, threads_per_block, ( lb_parameters_gpu, lb_particle_gpu, lb_particle_force_gpu ) );
}

}
