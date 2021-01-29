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
#include "config.hpp"

#include "ParticleRange.hpp"
#include "cuda_init.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.cuh"
#include "errorhandling.hpp"

#include "CudaDeviceAllocator.hpp"
#include "CudaHostAllocator.hpp"

#include <utils/constants.hpp>

#include <thrust/device_vector.h>

#include <cuda.h>

#include <cstddef>
#include <cstdio>

extern int this_node;

template <class T>
using device_vector = thrust::device_vector<T, CudaDeviceAllocator<T>>;

static CUDA_global_part_vars global_part_vars_host = {};

template <class T, class A>
T *raw_data_pointer(thrust::device_vector<T, A> &vec) {
  return thrust::raw_pointer_cast(vec.data());
}

template <class SpanLike> size_t byte_size(SpanLike const &v) {
  return v.size() * sizeof(typename SpanLike::value_type);
}

/** struct for particle force */
static device_vector<float> particle_forces_device;
static device_vector<float> particle_torques_device;

/** struct for particle position and velocity */
static device_vector<CUDA_particle_data> particle_data_device;
/** struct for energies */
static CUDA_energy *energy_device = nullptr;

pinned_vector<CUDA_particle_data> particle_data_host;
pinned_vector<float> particle_forces_host;
CUDA_energy energy_host;

pinned_vector<float> particle_torques_host;

cudaStream_t stream[1];

void cuda_check_errors_exit(const dim3 &block, const dim3 &grid,
                            const char *function, const char *file,
                            unsigned int line) {
  cudaError_t CU_err = cudaGetLastError();
  if (CU_err != cudaSuccess) {
    fprintf(stderr,
            "%d: error \"%s\" calling %s with dim %d %d %d, grid %d %d "
            "%d in %s:%u\n",
            this_node, cudaGetErrorString(CU_err), function, block.x, block.y,
            block.z, grid.x, grid.y, grid.z, file, line);
    errexit();
  }
}

/**
 * @brief Resize a @ref device_vector.
 *
 * Due to a bug in thrust (https://github.com/thrust/thrust/issues/939),
 * resizing or appending to default constructed containers causes undefined
 * behavior by dereferencing a null-pointer for certain types. This
 * function is used instead of the resize member function to side-step
 * the problem. This is done by replacing the existing vector by a new
 * one constructed with the desired size if resizing from capacity zero.
 * Behaves as-if vec.resize(n) was called.
 *
 * @tparam T Type contained in the vector.
 * @param vec Vector To resize.
 * @param n Desired new size of the element.
 */
template <class T> void resize_or_replace(device_vector<T> &vec, size_t n) {
  if (vec.capacity() == 0) {
    vec = device_vector<T>(n);
  } else {
    vec.resize(n);
  }
}

void resize_buffers(size_t number_of_particles) {
  particle_data_host.resize(number_of_particles);
  resize_or_replace(particle_data_device, number_of_particles);

  particle_forces_host.resize(3 * number_of_particles);
  resize_or_replace(particle_forces_device, 3 * number_of_particles);

#ifdef ROTATION
  particle_torques_host.resize(3 * number_of_particles);
  resize_or_replace(particle_torques_device, 3 * number_of_particles);
#endif
}

/**
 * @brief Setup and call particle reallocation from the host.
 * Note that in addition to calling this function the parameters must be
 * broadcast with either:
 * 1. @ref cuda_bcast_global_part_params() (when just being executed on the
 *    head node) or
 * 2. `MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
 *    sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart)` (when executed
 *    on all nodes)
 */
void gpu_init_particle_comm() {
  if (this_node == 0 && global_part_vars_host.communication_enabled == 0) {
    try {
      if (cuda_get_n_gpus() == 0) {
        fprintf(stderr, "ERROR: No GPU was found.\n");
        errexit();
      }
      auto const devID = cuda_get_device();
      auto const compute_capability = cuda_check_gpu_compute_capability(devID);
      auto const communication_test = cuda_test_device_access();
      if (compute_capability != ES_OK or communication_test != ES_OK) {
        fprintf(stderr,
                "ERROR: CUDA device %i is not capable of running ESPResSo.\n",
                devID);
        errexit();
      }
    } catch (cuda_runtime_error const &err) {
      fprintf(stderr, "ERROR: %s\n", err.what());
      errexit();
    }
  }
  global_part_vars_host.communication_enabled = 1;
}

Utils::Span<CUDA_particle_data> gpu_get_particle_pointer() {
  return {raw_data_pointer(particle_data_device), particle_data_device.size()};
}
CUDA_global_part_vars *gpu_get_global_particle_vars_pointer_host() {
  return &global_part_vars_host;
}
float *gpu_get_particle_force_pointer() {
  return raw_data_pointer(particle_forces_device);
}
CUDA_energy *gpu_get_energy_pointer() { return energy_device; }
float *gpu_get_particle_torque_pointer() {
  return raw_data_pointer(particle_torques_device);
}

void copy_part_data_to_gpu(ParticleRange particles) {
  if (global_part_vars_host.communication_enabled == 1) {
    cuda_mpi_get_particles(particles, particle_data_host);

    resize_buffers(particle_data_host.size());

    /* get espressomd particle values */
    if (this_node == 0) {
      cudaMemsetAsync(raw_data_pointer(particle_forces_device), 0x0,
                      byte_size(particle_forces_device), stream[0]);
#ifdef ROTATION
      cudaMemsetAsync(raw_data_pointer(particle_torques_device), 0x0,
                      byte_size(particle_torques_device), stream[0]);
#endif
      cudaMemcpyAsync(raw_data_pointer(particle_data_device),
                      particle_data_host.data(), byte_size(particle_data_host),
                      cudaMemcpyHostToDevice, stream[0]);
    }
  }
}

/** setup and call kernel to copy particle forces to host
 */
void copy_forces_from_GPU(ParticleRange &particles) {
  if (global_part_vars_host.communication_enabled == 1) {
    /* Copy result from device memory to host memory*/
    if (this_node == 0 && (not particle_forces_device.empty())) {
      thrust::copy(particle_forces_device.begin(), particle_forces_device.end(),
                   particle_forces_host.begin());
#ifdef ROTATION
      thrust::copy(particle_torques_device.begin(),
                   particle_torques_device.end(),
                   particle_torques_host.begin());
#endif
    }

    cuda_mpi_send_forces(
        particles, {particle_forces_host.data(), particle_forces_host.size()},
        {particle_torques_host.data(), particle_torques_host.size()});
  }
}

void clear_energy_on_GPU() {
  if (!global_part_vars_host.communication_enabled)
    return;
  if (energy_device == nullptr)
    cuda_safe_mem(cudaMalloc((void **)&energy_device, sizeof(CUDA_energy)));
  cuda_safe_mem(cudaMemset(energy_device, 0, sizeof(CUDA_energy)));
}

CUDA_energy copy_energy_from_GPU() {
  if (!global_part_vars_host.communication_enabled)
    return {};
  cuda_safe_mem(cudaMemcpy(&energy_host, energy_device, sizeof(CUDA_energy),
                           cudaMemcpyDeviceToHost));
  return energy_host;
}

void cuda_safe_mem_exit(cudaError_t CU_err, const char *file,
                        unsigned int line) {
  if (CU_err != cudaSuccess) {
    fprintf(stderr, "CUDA Memory error at %s:%u.\n", file, line);
    fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(CU_err));
    if (CU_err == cudaErrorInvalidValue)
      fprintf(stderr, "You may have tried to allocate zero memory at %s:%u.\n",
              file, line);
    errexit();
  } else {
    CU_err = cudaGetLastError();
    if (CU_err != cudaSuccess) {
      fprintf(stderr,
              "Error found during memory operation. Possibly however "
              "from a failed operation before. %s:%u.\n",
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
