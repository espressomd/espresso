/*
  Copyright (C) 2013 The ESPResSo project
  
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
#ifndef _CUDA_COMMON_H
#define _CUDA_COMMON_H

#include <cstdio>
#include "config.hpp" //this is required so that the ifdefs are actually defined

#ifdef CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "cuda_init.hpp"

extern "C" {

/** Action number for \ref mpi_get_particles. */
#define REQ_GETPARTS  16

/**cuda streams for parallel computing on cpu and gpu */
extern cudaStream_t stream[1];

extern cudaError_t err;
extern cudaError_t _err;

/** data which must be copied from the GPU at each step run on the GPU */
typedef struct {
  /** force on the particle given to md part */
  float f[3];

} CUDA_particle_force;

/** data structure which must be copied to the GPU at each step run on the GPU */
typedef struct {
  /** particle position given from md part*/
  float p[3];
  /** particle momentum struct velocity p.m->v*/
  float v[3];

#ifdef SHANCHEN
  float solvation[2*LB_COMPONENTS];
#endif 
#ifdef LB_ELECTROHYDRODYNAMICS
  float mu_E[3];
#endif
#ifdef ELECTROSTATICS
  float q;
#endif
  unsigned int fixed;

} CUDA_particle_data;


/** Note the particle's seed gets its own struct since it doesn't get copied back and forth from the GPU */
typedef struct {

  unsigned int seed;

} CUDA_particle_seed;

extern CUDA_particle_data *particle_data_host;


/** This structure contains global variables associated with all of the particles and not with one individual particle */
typedef struct {
  
  /**  This is for seeding the particles' individual seeds and is initialized using irandom, beware if using for other purposes */
  unsigned int seed;
  
  unsigned int number_of_particles; 
  
  /** a boolean variable to indicate if particle info should be communicated between the cpu and gpu */
  unsigned int communication_enabled;
} CUDA_global_part_vars;


void copy_forces_from_GPU();
CUDA_global_part_vars* gpu_get_global_particle_vars_pointer_host();
CUDA_global_part_vars* gpu_get_global_particle_vars_pointer();
CUDA_particle_data* gpu_get_particle_pointer();
CUDA_particle_force* gpu_get_particle_force_pointer();
CUDA_particle_seed* gpu_get_particle_seed_pointer();
void gpu_change_number_of_part_to_comm();
void gpu_init_particle_comm();
void cuda_mpi_get_particles(CUDA_particle_data *host_result);
void copy_part_data_to_gpu();
void cuda_mpi_send_forces(CUDA_particle_force *host_forces);
void cuda_bcast_global_part_params();

/**erroroutput for memory allocation and memory copy 
 * @param err cuda error code
 * @param *file .cu file were the error took place
 * @param line line of the file were the error took place
*/

void _cuda_safe_mem(cudaError_t err, char *file, unsigned int line);

#define cuda_safe_mem(a) _cuda_safe_mem((a), __FILE__, __LINE__)

#define KERNELCALL(_f, _a, _b, _params) \
_f<<<_a, _b, 0, stream[0]>>>_params; \
_err=cudaGetLastError(); \
if (_err!=cudaSuccess){ \
  printf("CUDA error: %s\n", cudaGetErrorString(_err)); \
  fprintf(stderr, "error calling %s with dim %d %d %d in %s:%u\n", #_f, _a.x, _a.y, _a.z, __FILE__, __LINE__); \
  exit(EXIT_FAILURE); \
 }
}
#endif /* ifdef CUDA */

#endif /* ifdef CUDA_COMMON_H */
