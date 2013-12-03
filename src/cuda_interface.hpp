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
#ifndef _CUDA_INTERFACE_HPP
#define _CUDA_INTERFACE_HPP

#include "config.hpp" //this is required so that the ifdefs are actually defined

#include "SystemInterface.hpp"

#ifdef CUDA

/** data which must be copied from the GPU at each step run on the GPU */
typedef struct {

  /** force on the particle given to md part */
  float f[3];

} CUDA_particle_force;


typedef struct {
  /** fluid composition at the particle given to md part */
  float weight[LB_COMPONENTS];

} CUDA_fluid_composition;


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
void copy_composition_from_GPU();
CUDA_global_part_vars* gpu_get_global_particle_vars_pointer_host();
CUDA_global_part_vars* gpu_get_global_particle_vars_pointer();
CUDA_particle_data* gpu_get_particle_pointer();
CUDA_particle_force* gpu_get_particle_force_pointer();
CUDA_fluid_composition* gpu_get_fluid_composition_pointer();
CUDA_particle_seed* gpu_get_particle_seed_pointer();
void gpu_change_number_of_part_to_comm();
void gpu_init_particle_comm();
void cuda_mpi_get_particles(CUDA_particle_data *host_result);
void copy_part_data_to_gpu();
void cuda_mpi_send_forces(CUDA_particle_force *host_forces,CUDA_fluid_composition * host_fluid_composition);
void cuda_bcast_global_part_params();
void cuda_copy_to_device(void *host_data, void *device_data, size_t n);
void cuda_copy_to_host(void *host_device, void *device_host, size_t n);
#endif /* ifdef CUDA */

#endif /* ifdef CUDA_INTERFACE_HPP */
