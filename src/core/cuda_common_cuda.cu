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

#include "cuda_wrapper.hpp"

#include "config.hpp"
#include "debug.hpp"

#include "cuda_init.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "errorhandling.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

#include <thrust/device_allocator.h>
#include <thrust/device_vector.h>
#include <utils/constants.hpp>

/**
 * @brief Allocator that uses CUDA to allocate CPU memory.
 *
 * Using the CUDA allocator can have performance benefits,
 * because it returns pinned memory that is suitable for
 * DMA.
 *
 * @tparam T Type to allocate memory for.
 */
template <class T> struct cuda_host_allocator {
  using value_type = T;
  using pointer = T *;
  using reference = T &;
  using const_reference = std::add_const_t<reference>;

  cuda_host_allocator() noexcept = default;
  template <class U>
  explicit cuda_host_allocator(const cuda_host_allocator<U> &) {}
  template <class U> bool operator==(const cuda_host_allocator<U> &) const {
    return true;
  }
  template <class U> bool operator!=(const cuda_host_allocator<U> &) const {
    return false;
  }

  T *allocate(const size_t n) const {
    T *result(0);
    cudaError_t error = cudaMallocHost(reinterpret_cast<void **>(&result),
                                       n * sizeof(value_type));

    if (error) {
      throw std::bad_alloc();
    }

    return result;
  }
  void deallocate(T *const p, size_t) const noexcept { cudaFreeHost(p); }
};

/**
 * @brief Wrapper around thrust::device_allocator.
 *
 * This is a thin wrapper around the thrust default device
 * allocator, which does ignore errors on deallocation. This
 * allows device containers with static lifetime, which might
 * be destroyed after the cuda api is finalized, which causes
 * the thrust allocator to throw. The implementation catches
 * all thrust::system::system_error during deallocation, otherwise
 * it works exactly as thrust::device_allocator.
 *
 * @tparam T Type to allocate memory for.
 */
template <class T>
struct cuda_device_allocator : public thrust::device_allocator<T> {
  using base_type = thrust::device_allocator<T>;
  using pointer = typename base_type::pointer;
  using size_type = typename base_type::size_type;
  using const_pointer = typename base_type::const_pointer;

  using base_type::address;
  using base_type::allocate;
  using base_type::base_type;
  using base_type::max_size;
  using base_type::rebind;

  void deallocate(pointer p, size_type cnt) {
    try {
      base_type::deallocate(p, cnt);
    } catch (thrust::system::system_error const &) {
      ;
    }
  }
};

template <class T>
using host_vector = thrust::host_vector<T, cuda_host_allocator<T>>;

template <class T>
using device_vector = thrust::device_vector<T, cuda_device_allocator<T>>;

static CUDA_global_part_vars global_part_vars_host = {};

template <class T, class A>
T *raw_data_pointer(thrust::device_vector<T, A> &vec) {
  return thrust::raw_pointer_cast(vec.data());
}

template <class T, class A>
size_t byte_size(thrust::device_vector<T, A> const &vec) {
  return vec.size() * sizeof(T);
}

/** struct for particle force */
static device_vector<float> particle_forces_device;
static device_vector<float> particle_torques_device;

/** struct for particle position and velocity */
static device_vector<CUDA_particle_data> particle_data_device;
/** struct for energies */
static CUDA_energy *energy_device = nullptr;

host_vector<CUDA_particle_data> particle_data_host;
std::vector<float> particle_forces_host;
CUDA_energy energy_host;

std::vector<float> particle_torques_host;

/**cuda streams for parallel computing on cpu and gpu */
cudaStream_t stream[1];

cudaError_t _err;
cudaError_t CU_err;

void _cuda_check_errors(const dim3 &block, const dim3 &grid,
                        const char *function, const char *file,
                        unsigned int line) {
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

/** change number of particles to be communicated to the GPU
 *  Note that in addition to calling this function the parameters must be
 * broadcast with either:
 * 1) cuda_bcast_global_part_params(); (when just being executed on the
 * master node) or 2) MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
 * sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart); (when executed on
 * all nodes)
 */
void gpu_change_number_of_part_to_comm() {
  // we only run the function if there are new particles which have been created
  // since the last call of this function

  if (global_part_vars_host.number_of_particles != n_part &&
      global_part_vars_host.communication_enabled == 1 && this_node == 0) {

    global_part_vars_host.number_of_particles = n_part;

    // if the arrays exists free them to prevent memory leaks
    particle_forces_host.clear();
    particle_forces_device.clear();
    particle_torques_device.clear();
    particle_data_device.clear();
#ifdef ROTATION
    particle_torques_host.clear();
#endif

    if (global_part_vars_host.number_of_particles) {
      /* pinned memory mode - use special function to get OS-pinned memory*/
      particle_data_host.resize(global_part_vars_host.number_of_particles);
      particle_forces_host.resize(3 *
                                  global_part_vars_host.number_of_particles);
#if (defined DIPOLES || defined ROTATION)
      particle_torques_host.resize(3 *
                                   global_part_vars_host.number_of_particles);
#endif

      particle_forces_device.resize(3 *
                                    global_part_vars_host.number_of_particles);
#ifdef ROTATION
      particle_torques_device.resize(3 *
                                     global_part_vars_host.number_of_particles);
#endif
      particle_data_device.resize(global_part_vars_host.number_of_particles);
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
               "by default. ESPResSo has detected a CUDA capable card but it "
               "is not the one used by ESPResSo by default. Please set the "
               "GPU to use by setting System.cuda_init_handle.device. A list "
               "of available GPUs is available through "
               "System.cuda_init_handle.device_list.";
      }
    }
  }
  global_part_vars_host.communication_enabled = 1;
  gpu_change_number_of_part_to_comm();
}

CUDA_particle_data *gpu_get_particle_pointer() {
  return raw_data_pointer(particle_data_device);
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
  if (global_part_vars_host.communication_enabled == 1 &&
      global_part_vars_host.number_of_particles) {
    cuda_mpi_get_particles(particles, particle_data_host.data());

    /* get espressomd particle values */
    if (this_node == 0) {
      cudaMemsetAsync(raw_data_pointer(particle_forces_device), 0x0,
                      byte_size(particle_forces_device), stream[0]);
      cudaMemsetAsync(raw_data_pointer(particle_torques_device), 0x0,
                      byte_size(particle_torques_device), stream[0]);
      cudaMemcpyAsync(raw_data_pointer(particle_data_device),
                      particle_data_host.data(),
                      particle_data_host.size() * sizeof(CUDA_particle_data),
                      cudaMemcpyHostToDevice, stream[0]);
    }
  }
}

/** setup and call kernel to copy particle forces to host
 */
void copy_forces_from_GPU(ParticleRange particles) {
  if (global_part_vars_host.communication_enabled == 1 &&
      global_part_vars_host.number_of_particles) {

    /* Copy result from device memory to host memory*/
    if (this_node == 0) {
      cuda_safe_mem(cudaMemcpy(
          &(particle_forces_host[0]), raw_data_pointer(particle_forces_device),
          3 * global_part_vars_host.number_of_particles * sizeof(float),
          cudaMemcpyDeviceToHost));
#ifdef ROTATION
      cuda_safe_mem(cudaMemcpy(&(particle_torques_host[0]),
                               raw_data_pointer(particle_torques_device),
                               global_part_vars_host.number_of_particles * 3 *
                                   sizeof(float),
                               cudaMemcpyDeviceToHost));
#endif
      cudaDeviceSynchronize();
    }
    using Utils::make_span;
    cuda_mpi_send_forces(particles, make_span(particle_forces_host),
                         make_span(particle_torques_host));
  }
}

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

/** @name Generic copy functions from and to device */

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

/*@}*/
