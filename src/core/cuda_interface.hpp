/*
  Copyright (C) 2013,2014,2015,2016 The ESPResSo project
  
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

#ifdef ENGINE
// velocities which need to be copied from the GPU to the CPU to calculate a torque
typedef struct {

  // center and source velocity of the md part
  float v_cs[6];

} CUDA_v_cs;
#endif





typedef struct {
  /** fluid composition at the particle given to md part */
  float weight[LB_COMPONENTS];

} CUDA_fluid_composition;

// Parameters for swimmers
#ifdef ENGINE
typedef struct {
  // v_cs has to stay in the front for memmove reasons
  float v_cs[6];
  float v_swim;
  float f_swim;
  float quatu[3];
  int push_pull;
  float dipole_length;
  bool swimming;
} CUDA_ParticleParametersSwimming;
#endif

/** data structure which must be copied to the GPU at each step run on the GPU */
typedef struct {

//   // This has to stay in front of the struct for memmove reasons
#ifdef ENGINE
  CUDA_ParticleParametersSwimming swim;
#endif
  
  /** particle position given from md part*/
  float p[3];
  /** particle momentum struct velocity p.m->v*/
  float v[3];

#ifdef ROTATION
  float quatu[3];
#endif

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
  
#ifdef IMMERSED_BOUNDARY
  bool isVirtual;
#endif

#ifdef DIPOLES
  float dip[3];
#endif

} CUDA_particle_data;

/** data structure for the different kinds of energies */
typedef struct {
  float bonded, non_bonded, coulomb, dipolar;
} CUDA_energy;

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
void copy_energy_from_GPU();
void copy_CUDA_energy_to_energy(CUDA_energy energy_host);
void clear_energy_on_GPU();
void copy_composition_from_GPU();

CUDA_global_part_vars* gpu_get_global_particle_vars_pointer_host();
CUDA_global_part_vars* gpu_get_global_particle_vars_pointer();
CUDA_particle_data* gpu_get_particle_pointer();
float* gpu_get_particle_force_pointer();
#ifdef ROTATION
float* gpu_get_particle_torque_pointer();
#endif

CUDA_energy* gpu_get_energy_pointer();
float* gpu_get_particle_torque_pointer();
CUDA_fluid_composition* gpu_get_fluid_composition_pointer();
CUDA_particle_seed* gpu_get_particle_seed_pointer();
void gpu_change_number_of_part_to_comm();
void gpu_init_particle_comm();
void cuda_mpi_get_particles(CUDA_particle_data *host_result);
void copy_part_data_to_gpu();
void cuda_mpi_send_forces(float* host_forces,float* host_torques,CUDA_fluid_composition * host_fluid_composition);
void cuda_bcast_global_part_params();
void cuda_copy_to_device(void *host_data, void *device_data, size_t n);
void cuda_copy_to_host(void *host_device, void *device_host, size_t n);

#ifdef ENGINE
void copy_v_cs_from_GPU();
void cuda_mpi_send_v_cs(CUDA_v_cs *host_v_cs);
#endif

#endif /* ifdef CUDA */

#endif /* ifdef CUDA_INTERFACE_HPP */
