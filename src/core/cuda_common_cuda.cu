/*
   Copyright (C) 2010-2018 The ESPResSo project

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

#include "cuda_wrapper.hpp"

#include "config.hpp"
#include "debug.hpp"

#include "cuda_init.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "errorhandling.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <utils/constants.hpp>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

static CUDA_global_part_vars global_part_vars_host = {};
__device__ __constant__ CUDA_global_part_vars global_part_vars_device[1];

/** struct for particle force */
static float *particle_forces_device = nullptr;
static float *particle_torques_device = nullptr;

/** struct for particle position and velocity */
static CUDA_particle_data *particle_data_device = nullptr;
/** struct for energies */
static CUDA_energy *energy_device = nullptr;

CUDA_particle_data *particle_data_host = nullptr;
std::vector<float> particle_forces_host;
CUDA_energy energy_host;

std::vector<float> particle_torques_host;
#ifdef ENGINE
std::vector<CUDA_v_cs> host_v_cs;
#endif

/**cuda streams for parallel computing on cpu and gpu */
cudaStream_t stream[1];

cudaError_t _err;
cudaError_t CU_err;

void _cuda_safe_mem(cudaError_t CU_err, const char *file, unsigned int line) {
  if (cudaSuccess != CU_err) {
    fprintf(stderr, "Cuda Memory error at %s:%u.\n", file, line);
    printf("CUDA error: %s\n", cudaGetErrorString(CU_err));
    if (CU_err == cudaErrorInvalidValue)
      fprintf(stderr, "You may have tried to allocate zero memory at %s:%u.\n",
              file, line);
    errexit();
  } else {
    CU_err = cudaGetLastError();
    if (CU_err != cudaSuccess) {
      fprintf(stderr,
              "Error found during memory operation. Possibly however "
              "from an failed operation before. %s:%u.\n",
              file, line);
      printf("CUDA error: %s\n", cudaGetErrorString(CU_err));
      if (CU_err == cudaErrorInvalidValue)
        fprintf(stderr,
                "You may have tried to allocate zero memory before %s:%u.\n",
                file, line);
      errexit();
    }
  }
}

void _cuda_check_errors(const dim3 &block, const dim3 &grid,
                        const char *function, const char *file,
                        unsigned int line) {
/** If debugging is enabled, wait for Kernels to terminate before checking for
 * errors. This removes parallelism between host and device and should only be
 * enabled while debugging. */
#ifdef CUDA_DEBUG
  cudaDeviceSynchronize();
#endif
  CU_err = cudaGetLastError();
  if (CU_err != cudaSuccess) {
    fprintf(stderr,
            "%d: error \"%s\" calling %s with dim %d %d %d, grid %d %d "
            "%d in %s:%u\n",
            this_node, cudaGetErrorString(CU_err), function, block.x, block.y,
            block.z, grid.x, grid.y, grid.z, file, line);
    errexit();
  }
}

__device__ unsigned int getThreadIndex() {

  return blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x +
         threadIdx.x;
}

/** Kernel for the initialisation of the particle force array
 * @param[out] particle_forces_device    Local particle force
 * @param[out] particle_torques_device   Local particle torque
 */
__global__ void init_particle_force(float *particle_forces_device,
                                    float *particle_torques_device) {

  unsigned int part_index = getThreadIndex();

  if (part_index < global_part_vars_device->number_of_particles) {
    particle_forces_device[3 * part_index + 0] = 0.0f;
    particle_forces_device[3 * part_index + 1] = 0.0f;
    particle_forces_device[3 * part_index + 2] = 0.0f;

#ifdef ROTATION
    particle_torques_device[3 * part_index] = 0.0f;
    particle_torques_device[3 * part_index + 1] = 0.0f;
    particle_torques_device[3 * part_index + 2] = 0.0f;
#endif
  }
}

/** Kernel for the initialisation of the particle force array
 * @param[out] particle_forces_device    Local particle force
 * @param[out] particle_torques_device   Local particle torque
 */
__global__ void reset_particle_force(float *particle_forces_device,
                                     float *particle_torques_device) {

  unsigned int part_index = getThreadIndex();

  if (part_index < global_part_vars_device->number_of_particles) {
    particle_forces_device[3 * part_index + 0] = 0.0f;
    particle_forces_device[3 * part_index + 1] = 0.0f;
    particle_forces_device[3 * part_index + 2] = 0.0f;
#ifdef ROTATION
    particle_torques_device[3 * part_index + 0] = 0.0f;
    particle_torques_device[3 * part_index + 1] = 0.0f;
    particle_torques_device[3 * part_index + 2] = 0.0f;
#endif
  }
}

/** change number of particles to be communicated to the GPU
 *  Note that in addition to calling this function the parameters must be
 * broadcast with either:
 * 1) cuda_bcast_global_part_params(); (when just being executed on the master
 * node) or
 * 2) MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
 * sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart); (when executed on all
 * nodes)
 */
void gpu_change_number_of_part_to_comm() {
  // we only run the function if there are new particles which have been created
  // since the last call of this function

  if (global_part_vars_host.number_of_particles != n_part &&
      global_part_vars_host.communication_enabled == 1 && this_node == 0) {

    global_part_vars_host.number_of_particles = n_part;

    cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(global_part_vars_device),
                                     &global_part_vars_host,
                                     sizeof(CUDA_global_part_vars)));

    // if the arrays exists free them to prevent memory leaks
    particle_forces_host.clear();
    if (particle_data_host) {
      cuda_safe_mem(cudaFreeHost(particle_data_host));
      particle_data_host = nullptr;
    }
    if (particle_forces_device) {
      cudaFree(particle_forces_device);
      particle_forces_device = nullptr;
    }
    if (particle_data_device) {
      cudaFree(particle_data_device);
      particle_data_device = nullptr;
    }

#ifdef ENGINE
    host_v_cs.clear();
#endif
#ifdef ROTATION
    particle_torques_host.clear();
#endif

#ifdef ROTATION
    if (particle_torques_device) {
      cuda_safe_mem(cudaFree(particle_torques_device));
      particle_torques_device = nullptr;
    }
#endif

    if (global_part_vars_host.number_of_particles) {

      /**pinned memory mode - use special function to get OS-pinned memory*/
      cuda_safe_mem(cudaHostAlloc((void **)&particle_data_host,
                                  global_part_vars_host.number_of_particles *
                                      sizeof(CUDA_particle_data),
                                  cudaHostAllocWriteCombined));
      particle_forces_host.resize(3 *
                                  global_part_vars_host.number_of_particles);
#ifdef ENGINE
      host_v_cs.resize(global_part_vars_host.number_of_particles);
#endif
#if (defined DIPOLES || defined ROTATION)
      particle_torques_host.resize(3 *
                                   global_part_vars_host.number_of_particles);
#endif

      cuda_safe_mem(cudaMalloc((void **)&particle_forces_device,
                               3 * global_part_vars_host.number_of_particles *
                                   sizeof(float)));
#ifdef ROTATION
      cuda_safe_mem(cudaMalloc((void **)&particle_torques_device,
                               3 * global_part_vars_host.number_of_particles *
                                   sizeof(float)));
#endif

      cuda_safe_mem(cudaMalloc((void **)&particle_data_device,
                               global_part_vars_host.number_of_particles *
                                   sizeof(CUDA_particle_data)));

      /** values for the particle kernel */
      int threads_per_block_particles = 64;
      int blocks_per_grid_particles_y = 4;
      int blocks_per_grid_particles_x =
          (global_part_vars_host.number_of_particles +
           threads_per_block_particles * blocks_per_grid_particles_y - 1) /
          (threads_per_block_particles * blocks_per_grid_particles_y);
      dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x,
                                           blocks_per_grid_particles_y, 1);

      KERNELCALL(init_particle_force, dim_grid_particles,
                 threads_per_block_particles, particle_forces_device,
                 particle_torques_device);
    }
  }
}

/** setup and call particle reallocation from the host
 *  Note that in addition to calling this function the parameters must be
 * broadcast with either:
 * 1) cuda_bcast_global_part_params(); (when just being executed on the master
 * node) or
 * 2) MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
 * sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart); (when executed on all
 * nodes)
 */
void gpu_init_particle_comm() {
  if (this_node == 0 && global_part_vars_host.communication_enabled == 0) {
    if (cuda_get_n_gpus() == -1) {
      runtimeErrorMsg()
          << "Unable to initialize CUDA as no sufficient GPU is available.";
      errexit();
    }
    if (cuda_get_n_gpus() > 1) {
      runtimeWarningMsg() << "More than one GPU detected, please note ESPResSo "
                             "uses device 0 by default regardless of usage or "
                             "capability. The GPU to be used can be modified "
                             "by setting System.cuda_init_handle.device.";
      if (cuda_check_gpu(0) != ES_OK) {
        runtimeWarningMsg()
            << "CUDA device 0 is not capable of running ESPResSo but is used "
               "by default. Espresso has detected a CUDA capable card but it "
               "is not the one used by ESPResSo by default. Please set the "
               "GPU to use by setting System.cuda_init_handle.device. A list "
               "of avalable GPUs is available through "
               "System.cuda_init_handle.device_list.";
      }
    }
  }
  global_part_vars_host.communication_enabled = 1;
  gpu_change_number_of_part_to_comm();
}

CUDA_particle_data *gpu_get_particle_pointer() { return particle_data_device; }
CUDA_global_part_vars *gpu_get_global_particle_vars_pointer_host() {
  return &global_part_vars_host;
}
CUDA_global_part_vars *gpu_get_global_particle_vars_pointer() {
  return global_part_vars_device;
}
float *gpu_get_particle_force_pointer() { return particle_forces_device; }
CUDA_energy *gpu_get_energy_pointer() { return energy_device; }
float *gpu_get_particle_torque_pointer() { return particle_torques_device; }

void copy_part_data_to_gpu(ParticleRange particles) {
  COMM_TRACE(printf("global_part_vars_host.communication_enabled = %d && "
                    "global_part_vars_host.number_of_particles = %d\n",
                    global_part_vars_host.communication_enabled,
                    global_part_vars_host.number_of_particles));
  if (global_part_vars_host.communication_enabled == 1 &&
      global_part_vars_host.number_of_particles) {
    cuda_mpi_get_particles(particles, particle_data_host);

    /** get espresso md particle values*/
    if (this_node == 0)
      cudaMemcpyAsync(particle_data_device, particle_data_host,
                      global_part_vars_host.number_of_particles *
                          sizeof(CUDA_particle_data),
                      cudaMemcpyHostToDevice, stream[0]);
  }
}

/** setup and call kernel to copy particle forces to host
 */
void copy_forces_from_GPU(ParticleRange particles) {

  if (global_part_vars_host.communication_enabled == 1 &&
      global_part_vars_host.number_of_particles) {

    /** Copy result from device memory to host memory*/
    if (this_node == 0) {
      cuda_safe_mem(cudaMemcpy(
          &(particle_forces_host[0]), particle_forces_device,
          3 * global_part_vars_host.number_of_particles * sizeof(float),
          cudaMemcpyDeviceToHost));
#ifdef ROTATION
      cuda_safe_mem(cudaMemcpy(
          &(particle_torques_host[0]), particle_torques_device,
          global_part_vars_host.number_of_particles * 3 * sizeof(float),
          cudaMemcpyDeviceToHost));
#endif

      /** values for the particle kernel */
      int threads_per_block_particles = 64;
      int blocks_per_grid_particles_y = 4;
      int blocks_per_grid_particles_x =
          (global_part_vars_host.number_of_particles +
           threads_per_block_particles * blocks_per_grid_particles_y - 1) /
          (threads_per_block_particles * blocks_per_grid_particles_y);
      dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x,
                                           blocks_per_grid_particles_y, 1);

      /** reset part forces with zero*/

      KERNELCALL(reset_particle_force, dim_grid_particles,
                 threads_per_block_particles, particle_forces_device,
                 particle_torques_device);
      cudaDeviceSynchronize();
    }

    cuda_mpi_send_forces(particles, particle_forces_host,
                         particle_torques_host);
  }
}

#if defined(ENGINE) && defined(LB_GPU)
// setup and call kernel to copy v_cs to host
void copy_v_cs_from_GPU(ParticleRange particles) {
  if (global_part_vars_host.communication_enabled == 1 &&
      global_part_vars_host.number_of_particles) {
    // Copy result from device memory to host memory
    if (this_node == 0) {
      cuda_safe_mem(cudaMemcpy2D(
          host_v_cs.data(), sizeof(CUDA_v_cs), particle_data_device,
          sizeof(CUDA_particle_data), sizeof(CUDA_v_cs),
          global_part_vars_host.number_of_particles, cudaMemcpyDeviceToHost));
    }
    cuda_mpi_send_v_cs(particles, host_v_cs);
  }
}
#endif

void clear_energy_on_GPU() {
  if (!global_part_vars_host.communication_enabled)
    // || !global_part_vars_host.number_of_particles )
    return;
  if (energy_device == nullptr)
    cuda_safe_mem(cudaMalloc((void **)&energy_device, sizeof(CUDA_energy)));
  cuda_safe_mem(cudaMemset(energy_device, 0, sizeof(CUDA_energy)));
}

void copy_energy_from_GPU() {
  if (!global_part_vars_host.communication_enabled ||
      !global_part_vars_host.number_of_particles)
    return;
  cuda_safe_mem(cudaMemcpy(&energy_host, energy_device, sizeof(CUDA_energy),
                           cudaMemcpyDeviceToHost));
  copy_CUDA_energy_to_energy(energy_host);
}

/** Generic copy functions from an to device **/

void cuda_copy_to_device(void *host_data, void *device_data, size_t n) {
  cuda_safe_mem(cudaMemcpy(host_data, device_data, n, cudaMemcpyHostToDevice));
}

void cuda_copy_to_host(void *host_device, void *device_host, size_t n) {
  cuda_safe_mem(
      cudaMemcpy(host_device, device_host, n, cudaMemcpyDeviceToHost));
}
