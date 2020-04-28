/*
 * Copyright (C) 2013-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CORE_CUDA_INTERFACE_HPP
#define CORE_CUDA_INTERFACE_HPP

#include "config.hpp"

#ifdef CUDA

#include "CudaHostAllocator.hpp"
#include "ParticleRange.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

// Parameters for swimmers
#ifdef ENGINE
struct CUDA_ParticleParametersSwimming {
  using Vector3f = Utils::Vector3f;

  // v_cs has to stay in the front for memmove reasons
  float v_cs[6];
  float v_swim;
  float f_swim;
  Vector3f director;
  int push_pull;
  float dipole_length;
  bool swimming;
};
#endif

/** data structure which must be copied to the GPU at each step run on the GPU
 */
struct CUDA_particle_data {
  using Vector3f = Utils::Vector3f;
//   // This has to stay in front of the struct for memmove reasons
#ifdef ENGINE
  CUDA_ParticleParametersSwimming swim;
#endif

  /** particle position given from md part*/
  Vector3f p;

#if defined(CUDA)
  /** particle id */
  int identity;
#ifdef VIRTUAL_SITES
  bool is_virtual;
#else
  static constexpr const bool is_virtual = false;
#endif

  /** particle momentum struct velocity p.m->v*/
  Vector3f v;
#endif

#ifdef ROTATION
  Vector3f director;
#endif

#if defined(LB_ELECTROHYDRODYNAMICS) && defined(CUDA)
  Vector3f mu_E;
#endif

#ifdef ELECTROSTATICS
  float q;
#endif

#ifdef MASS
  float mass;
#endif

#ifdef DIPOLES
  Vector3f dip;
#endif
};

/** data structure for the different kinds of energies */
typedef struct {
  float bonded, non_bonded, coulomb, dipolar;
} CUDA_energy;

/** Global variables associated with all of the particles and not with
 *  one individual particle.
 */
typedef struct {
  /** Boolean flag to indicate if particle info should be communicated
   *  between the cpu and gpu
   */
  unsigned int communication_enabled;
} CUDA_global_part_vars;

void copy_forces_from_GPU(ParticleRange &particles);
void copy_energy_from_GPU();
void copy_CUDA_energy_to_energy(CUDA_energy energy_host);
void clear_energy_on_GPU();

CUDA_global_part_vars *gpu_get_global_particle_vars_pointer_host();
Utils::Span<CUDA_particle_data> gpu_get_particle_pointer();
float *gpu_get_particle_force_pointer();
#ifdef ROTATION
float *gpu_get_particle_torque_pointer();
#endif

CUDA_energy *gpu_get_energy_pointer();
float *gpu_get_particle_torque_pointer();
void gpu_init_particle_comm();

void cuda_mpi_get_particles(
    const ParticleRange &particles,
    pinned_vector<CUDA_particle_data> &particle_data_host);
void copy_part_data_to_gpu(ParticleRange particles);

/**
 * @brief Distribute forces to the slaves, and add them to the particles.
 *
 * @param particles    The particles for which the forces (and torques) should
 *                     be added to.
 * @param host_forces  The forces as flat array of size 3 * particles.size(),
 *                     only relevant on the master.
 * @param host_torques The torques as flat array of size 3 * particles.size(),
 *                     this is only touched if ROTATION is active. Only
 *                     relevant on the master.
 *
 * This is a collective call.
 */
void cuda_mpi_send_forces(const ParticleRange &particles,
                          Utils::Span<float> host_forces,
                          Utils::Span<float> host_torques);
void cuda_bcast_global_part_params();
#endif /* ifdef CUDA */

#endif /* ifdef CUDA_INTERFACE_HPP */
