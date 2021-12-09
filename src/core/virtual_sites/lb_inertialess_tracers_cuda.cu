/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

// This is an internal file of the IMMERSED BOUNDARY implementation
// It should not be included by any main ESPResSo routines
// Functions to be exported for ESPResSo are in ibm_main.hpp

#include "config.hpp"

#if defined(VIRTUAL_SITES_INERTIALESS_TRACERS) && defined(CUDA)

#include "virtual_sites/lb_inertialess_tracers.hpp"
#include "virtual_sites/lb_inertialess_tracers_cuda_interface.hpp"

#include "Particle.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.cuh"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lbgpu.cuh"
#include "grid_based_algorithms/lbgpu.hpp"

#include <cuda.h>

#include <cstddef>

// Other functions for internal use
void InitCUDA_IBM(std::size_t numParticles);

// Our own global variables
IBM_CUDA_ParticleDataInput *IBM_ParticleDataInput_device = nullptr;
IBM_CUDA_ParticleDataOutput *IBM_ParticleDataOutput_device = nullptr;
bool IBM_initialized = false;
std::size_t IBM_numParticlesCache = 0; // To detect a change in particle number
                                       // which requires reallocation of memory

// These variables are defined in lbgpu_cuda.cu, but we also want them here
extern LB_node_force_density_gpu node_f;
extern LB_nodes_gpu *current_nodes;

// These variables are static in lbgpu_cuda.cu, so we need to duplicate them
// here. They are initialized in ForcesIntoFluid. The pointers are on the host,
// but point into device memory.
LB_parameters_gpu *para_gpu = nullptr;
float *lb_boundary_velocity_IBM = nullptr;

static constexpr unsigned int threads_per_block = 64;

__global__ void
ForcesIntoFluid_Kernel(const IBM_CUDA_ParticleDataInput *const particle_input,
                       std::size_t number_of_particles,
                       LB_node_force_density_gpu node_f,
                       const LB_parameters_gpu *const paraP) {
  const unsigned int particleIndex = blockIdx.y * gridDim.x * blockDim.x +
                                     blockDim.x * blockIdx.x + threadIdx.x;
  const LB_parameters_gpu &para = *paraP;

  if (particleIndex < number_of_particles &&
      particle_input[particleIndex].is_virtual) {
    // MD to LB units: mass is not affected, length are scaled by agrid, times
    // by para.tau
    const float factor = 1 / para.agrid * para.tau * para.tau;
    const float particleForce[3] = {particle_input[particleIndex].f[0] * factor,
                                    particle_input[particleIndex].f[1] * factor,
                                    particle_input[particleIndex].f[2] *
                                        factor};
    const float pos[3] = {particle_input[particleIndex].pos[0],
                          particle_input[particleIndex].pos[1],
                          particle_input[particleIndex].pos[2]};

    // First part is the same as for interpolation --> merge into a single
    // function
    float temp_delta[6];
    float delta[8];
    int my_left[3];
    unsigned int node_index[8];
    for (int i = 0; i < 3; ++i) {
      const float scaledpos = pos[i] / para.agrid - 0.5f;
      my_left[i] = static_cast<int>(floorf(scaledpos));
      temp_delta[3 + i] = scaledpos - static_cast<float>(my_left[i]);
      temp_delta[i] = 1.f - temp_delta[3 + i];
    }

    delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
    delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
    delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
    delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
    delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
    delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
    delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
    delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

    // modulo for negative numbers is strange at best, shift to make sure we are
    // positive
    auto const x = static_cast<unsigned int>(my_left[0] + para.dim[0]);
    auto const y = static_cast<unsigned int>(my_left[1] + para.dim[1]);
    auto const z = static_cast<unsigned int>(my_left[2] + para.dim[2]);

    node_index[0] = x % para.dim[0] + para.dim[0] * (y % para.dim[1]) +
                    para.dim[0] * para.dim[1] * (z % para.dim[2]);
    node_index[1] = (x + 1) % para.dim[0] + para.dim[0] * (y % para.dim[1]) +
                    para.dim[0] * para.dim[1] * (z % para.dim[2]);
    node_index[2] = x % para.dim[0] + para.dim[0] * ((y + 1) % para.dim[1]) +
                    para.dim[0] * para.dim[1] * (z % para.dim[2]);
    node_index[3] = (x + 1) % para.dim[0] +
                    para.dim[0] * ((y + 1) % para.dim[1]) +
                    para.dim[0] * para.dim[1] * (z % para.dim[2]);
    node_index[4] = x % para.dim[0] + para.dim[0] * (y % para.dim[1]) +
                    para.dim[0] * para.dim[1] * ((z + 1) % para.dim[2]);
    node_index[5] = (x + 1) % para.dim[0] + para.dim[0] * (y % para.dim[1]) +
                    para.dim[0] * para.dim[1] * ((z + 1) % para.dim[2]);
    node_index[6] = x % para.dim[0] + para.dim[0] * ((y + 1) % para.dim[1]) +
                    para.dim[0] * para.dim[1] * ((z + 1) % para.dim[2]);
    node_index[7] = (x + 1) % para.dim[0] +
                    para.dim[0] * ((y + 1) % para.dim[1]) +
                    para.dim[0] * para.dim[1] * ((z + 1) % para.dim[2]);

    for (int i = 0; i < 8; ++i) {
      // Atomic add is essential because this runs in parallel!
      atomicAdd(&(node_f.force_density[node_index[i]][0]),
                (particleForce[0] * delta[i]));
      atomicAdd(&(node_f.force_density[node_index[i]][1]),
                (particleForce[1] * delta[i]));
      atomicAdd(&(node_f.force_density[node_index[i]][2]),
                (particleForce[2] * delta[i]));
    }
  }
}

__global__ void ParticleVelocitiesFromLB_Kernel(
    LB_nodes_gpu n_curr,
    const IBM_CUDA_ParticleDataInput *const particles_input,
    std::size_t number_of_particles,
    IBM_CUDA_ParticleDataOutput *const particles_output,
    LB_node_force_density_gpu node_f, const float *const lb_boundary_velocity,
    const LB_parameters_gpu *const paraP) {

  const unsigned int particleIndex = blockIdx.y * gridDim.x * blockDim.x +
                                     blockDim.x * blockIdx.x + threadIdx.x;

  const LB_parameters_gpu &para = *paraP;

  if (particleIndex < number_of_particles &&
      particles_input[particleIndex].is_virtual) {

    // Get position
    float pos[3] = {particles_input[particleIndex].pos[0],
                    particles_input[particleIndex].pos[1],
                    particles_input[particleIndex].pos[2]};
    float v[3] = {0};

    // This part is copied from get_interpolated_velocity
    // + we add the force + we consider boundaries

    float temp_delta[6];
    float delta[8];
    int my_left[3];
    unsigned int node_index[8];
    Utils::Array<float, 4> mode;
#pragma unroll
    for (int i = 0; i < 3; ++i) {
      const float scaledpos = pos[i] / para.agrid - 0.5f;
      my_left[i] = static_cast<int>(floorf(scaledpos));
      temp_delta[3 + i] = scaledpos - static_cast<float>(my_left[i]);
      temp_delta[i] = 1.f - temp_delta[3 + i];
    }

    delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
    delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
    delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
    delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
    delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
    delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
    delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
    delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

    // modulo for negative numbers is strange at best, shift to make sure we are
    // positive
    auto const x = static_cast<unsigned int>(my_left[0] + para.dim[0]);
    auto const y = static_cast<unsigned int>(my_left[1] + para.dim[1]);
    auto const z = static_cast<unsigned int>(my_left[2] + para.dim[2]);

    node_index[0] = x % para.dim[0] + para.dim[0] * (y % para.dim[1]) +
                    para.dim[0] * para.dim[1] * (z % para.dim[2]);
    node_index[1] = (x + 1) % para.dim[0] + para.dim[0] * (y % para.dim[1]) +
                    para.dim[0] * para.dim[1] * (z % para.dim[2]);
    node_index[2] = x % para.dim[0] + para.dim[0] * ((y + 1) % para.dim[1]) +
                    para.dim[0] * para.dim[1] * (z % para.dim[2]);
    node_index[3] = (x + 1) % para.dim[0] +
                    para.dim[0] * ((y + 1) % para.dim[1]) +
                    para.dim[0] * para.dim[1] * (z % para.dim[2]);
    node_index[4] = x % para.dim[0] + para.dim[0] * (y % para.dim[1]) +
                    para.dim[0] * para.dim[1] * ((z + 1) % para.dim[2]);
    node_index[5] = (x + 1) % para.dim[0] + para.dim[0] * (y % para.dim[1]) +
                    para.dim[0] * para.dim[1] * ((z + 1) % para.dim[2]);
    node_index[6] = x % para.dim[0] + para.dim[0] * ((y + 1) % para.dim[1]) +
                    para.dim[0] * para.dim[1] * ((z + 1) % para.dim[2]);
    node_index[7] = (x + 1) % para.dim[0] +
                    para.dim[0] * ((y + 1) % para.dim[1]) +
                    para.dim[0] * para.dim[1] * ((z + 1) % para.dim[2]);

    for (int i = 0; i < 8; ++i) {
      double local_rho;
      double local_j[3];
#ifdef LB_BOUNDARIES_GPU
      if (n_curr.boundary[node_index[i]]) {
        // Boundary node
        auto const boundary_index =
            static_cast<int>(n_curr.boundary[node_index[i]]);

        // lb_boundary_velocity is given in MD units --> convert to LB and
        // reconvert back at the end of this function
        local_rho = para.rho;
        local_j[0] =
            para.rho * lb_boundary_velocity[3 * (boundary_index - 1) + 0];
        local_j[1] =
            para.rho * lb_boundary_velocity[3 * (boundary_index - 1) + 1];
        local_j[2] =
            para.rho * lb_boundary_velocity[3 * (boundary_index - 1) + 2];

      } else
#endif
      {
        calc_mass_and_momentum_mode(mode, n_curr, node_index[i]);
        local_rho = para.rho + mode[0];

        // Add the +f/2 contribution!!
        local_j[0] = mode[1] + node_f.force_density_buf[node_index[i]][0] / 2.f;
        local_j[1] = mode[2] + node_f.force_density_buf[node_index[i]][1] / 2.f;
        local_j[2] = mode[3] + node_f.force_density_buf[node_index[i]][2] / 2.f;
      }

      // Interpolate velocity
      v[0] += static_cast<float>(delta[i] * local_j[0] / local_rho);
      v[1] += static_cast<float>(delta[i] * local_j[1] / local_rho);
      v[2] += static_cast<float>(delta[i] * local_j[2] / local_rho);
    }

    // Rescale and store output
    particles_output[particleIndex].v[0] = v[0] * para.agrid / para.tau;
    particles_output[particleIndex].v[1] = v[1] * para.agrid / para.tau;
    particles_output[particleIndex].v[2] = v[2] * para.agrid / para.tau;
  }
}

__global__ void ResetLBForces_Kernel(LB_node_force_density_gpu node_f,
                                     const LB_parameters_gpu *const paraP) {

  const std::size_t index = blockIdx.y * gridDim.x * blockDim.x +
                            blockDim.x * blockIdx.x + threadIdx.x;
  const LB_parameters_gpu &para = *paraP;

  if (index < para.number_of_nodes) {
    const float force_factor = powf(para.agrid, 2) * para.tau * para.tau;
    if (para.external_force_density) {
      node_f.force_density[index][0] = para.ext_force_density[0] * force_factor;
      node_f.force_density[index][1] = para.ext_force_density[1] * force_factor;
      node_f.force_density[index][2] = para.ext_force_density[2] * force_factor;
    } else {
      node_f.force_density[index] = {};
    }
  }
}

/** Transfer particle forces into the LB fluid.
 *  Called from @ref integrate.
 *  This must be the first CUDA-IBM function to be called because it also does
 *  some initialization.
 */
void IBM_ForcesIntoFluid_GPU(ParticleRange const &particles, int this_node) {
  // This function does
  // (1) Gather forces from all particles via MPI
  // (2) Copy forces to the GPU
  // (3) interpolate on the LBM grid and spread forces

  auto const numParticles = gpu_get_particle_pointer().size();

  // Storage only needed on head node
  if (this_node == 0 &&
      (IBM_ParticleDataInput_host.empty() || !IBM_initialized ||
       numParticles != IBM_numParticlesCache))
    InitCUDA_IBM(numParticles);

  // We gather particle positions and forces from all nodes
  IBM_cuda_mpi_get_particles(particles);

  // GPU only on head node
  if (this_node == 0 && numParticles > 0) {

    // Copy data to device
    cuda_safe_mem(cudaMemcpy(IBM_ParticleDataInput_device,
                             IBM_ParticleDataInput_host.data(),
                             numParticles * sizeof(IBM_CUDA_ParticleDataInput),
                             cudaMemcpyHostToDevice));

    // Kernel call for spreading the forces on the LB grid
    dim3 dim_grid = calculate_dim_grid(static_cast<unsigned>(numParticles), 4,
                                       threads_per_block);
    KERNELCALL(ForcesIntoFluid_Kernel, dim_grid, threads_per_block,
               IBM_ParticleDataInput_device, numParticles, node_f, para_gpu);
  }
}

void InitCUDA_IBM(std::size_t const numParticles) {

  // Check if we have to delete
  if (!IBM_ParticleDataInput_host.empty()) {
    IBM_ParticleDataInput_host.clear();
    IBM_ParticleDataOutput_host.clear();
    cuda_safe_mem(cudaFree(IBM_ParticleDataInput_device));
    cuda_safe_mem(cudaFree(IBM_ParticleDataOutput_device));
    cuda_safe_mem(cudaFree(lb_boundary_velocity_IBM));
  }

  // Back and forth communication of positions and velocities
  IBM_ParticleDataInput_host.resize(numParticles);
  IBM_ParticleDataOutput_host.resize(numParticles);
  cuda_safe_mem(cudaMalloc((void **)&IBM_ParticleDataInput_device,
                           numParticles * sizeof(IBM_CUDA_ParticleDataInput)));
  cuda_safe_mem(cudaMalloc((void **)&IBM_ParticleDataOutput_device,
                           numParticles * sizeof(IBM_CUDA_ParticleDataOutput)));

  // Use LB parameters
  lb_get_para_pointer(&para_gpu);

  // Copy boundary velocities to the GPU
  // First put them into correct format
#ifdef LB_BOUNDARIES_GPU
  auto *host_lb_boundary_velocity =
      new float[3 * (LBBoundaries::lbboundaries.size() + 1)];

  for (int n = 0; n < LBBoundaries::lbboundaries.size(); n++) {
    host_lb_boundary_velocity[3 * n + 0] =
        static_cast<float>(LBBoundaries::lbboundaries[n]->velocity()[0]);
    host_lb_boundary_velocity[3 * n + 1] =
        static_cast<float>(LBBoundaries::lbboundaries[n]->velocity()[1]);
    host_lb_boundary_velocity[3 * n + 2] =
        static_cast<float>(LBBoundaries::lbboundaries[n]->velocity()[2]);
  }

  host_lb_boundary_velocity[3 * LBBoundaries::lbboundaries.size() + 0] = 0.0f;
  host_lb_boundary_velocity[3 * LBBoundaries::lbboundaries.size() + 1] = 0.0f;
  host_lb_boundary_velocity[3 * LBBoundaries::lbboundaries.size() + 2] = 0.0f;

  cuda_safe_mem(
      cudaMalloc((void **)&lb_boundary_velocity_IBM,
                 3 * LBBoundaries::lbboundaries.size() * sizeof(float)));
  cuda_safe_mem(
      cudaMemcpy(lb_boundary_velocity_IBM, host_lb_boundary_velocity,
                 3 * LBBoundaries::lbboundaries.size() * sizeof(float),
                 cudaMemcpyHostToDevice));

  delete[] host_lb_boundary_velocity;
#endif

  IBM_numParticlesCache = numParticles;
  IBM_initialized = true;
}

/** Call a kernel function to interpolate the velocity at each IBM particle's
 *  position. Store velocity in the particle data structure.
 */
void ParticleVelocitiesFromLB_GPU(ParticleRange const &particles,
                                  int this_node) {
  // This function performs three steps:
  // (1) interpolate velocities on GPU
  // (2) transfer velocities back to CPU
  // (3) spread velocities to local cells via MPI

  auto const numParticles = gpu_get_particle_pointer().size();

  // GPU only on head node
  if (this_node == 0 && numParticles > 0) {
    // Kernel call
    dim3 dim_grid = calculate_dim_grid(static_cast<unsigned>(numParticles), 4,
                                       threads_per_block);
    KERNELCALL(ParticleVelocitiesFromLB_Kernel, dim_grid, threads_per_block,
               *current_nodes, IBM_ParticleDataInput_device, numParticles,
               IBM_ParticleDataOutput_device, node_f, lb_boundary_velocity_IBM,
               para_gpu);

    // Copy velocities from device to host
    cuda_safe_mem(cudaMemcpy(IBM_ParticleDataOutput_host.data(),
                             IBM_ParticleDataOutput_device,
                             numParticles * sizeof(IBM_CUDA_ParticleDataOutput),
                             cudaMemcpyDeviceToHost));
  }

  // Scatter to all nodes
  IBM_cuda_mpi_send_velocities(particles);
}

#endif
