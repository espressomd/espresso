/*
 * Copyright (C) 2014-2022 The ESPResSo project
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
/**
 * @file
 * CUDA kernels to convert the particles AoS to a SoA on the device.
 */

#include "config/config.hpp"

#include "GpuParticleData.hpp"
#include "ResourceCleanup.hpp"
#include "System.hpp"

#include "ParticleRange.hpp"
#include "errorhandling.hpp"

#include "cuda/init.hpp"
#include "cuda/utils.cuh"

#include <utils/Span.hpp>

#include <thrust/copy.h>
#include <thrust/device_vector.h>

#include <cuda.h>

#include <cstddef>
#include <cstdio>
#include <memory>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

template <class T> T *raw_data_pointer(thrust::device_vector<T> &vec) {
  return thrust::raw_pointer_cast(vec.data());
}

template <class SpanLike> std::size_t byte_size(SpanLike const &v) {
  return v.size() * sizeof(typename SpanLike::value_type);
}

/**
 * @brief Resize a @c thrust::device_vector.
 *
 * Due to a bug in thrust (https://github.com/thrust/thrust/issues/939),
 * resizing or appending to default constructed containers causes undefined
 * behavior by dereferencing a null-pointer for certain types. This
 * function is used instead of the resize member function to side-step
 * the problem. This is done by replacing the existing vector by a new
 * one constructed with the desired size if resizing from capacity zero.
 * Behaves as-if @c vec.resize(n) was called.
 * This is fixed in Thrust 1.11, shipped in CUDA 11.3
 * (https://github.com/NVIDIA/thrust/commit/1c4f25d9).
 *
 * @tparam T Type contained in the vector.
 * @param vec Vector to resize.
 * @param n Desired new size of the vector.
 */
template <class T>
void resize_or_replace(thrust::device_vector<T> &vec, std::size_t n) {
  if (vec.capacity() == 0) {
    vec = thrust::device_vector<T>(n);
  } else {
    vec.resize(n);
  }
}

template <typename T> void free_device_vector(thrust::device_vector<T> &vec) {
  vec.clear();
  thrust::device_vector<T>().swap(vec);
}

/** @brief Host and device containers for particle data. */
class GpuParticleData::Storage {
  void free_device_memory();
  using DeviceMemory = ResourceCleanup::Attorney<&Storage::free_device_memory>;
  friend DeviceMemory;

public:
  /** @brief Which particle properties are needed by GPU methods. */
  GpuParticleData::prop::bitset m_need;
  GpuParticleData::GpuEnergy *energy_device = nullptr;
  std::size_t current_size = 0ul;
  pinned_vector<GpuParticle> particle_data_host;
  thrust::device_vector<GpuParticle> particle_data_device;
  pinned_vector<float> particle_forces_host;
  thrust::device_vector<float> particle_forces_device;
#ifdef ROTATION
  pinned_vector<float> particle_torques_host;
  thrust::device_vector<float> particle_torques_device;
#endif
  float *particle_pos_device = nullptr;
#ifdef DIPOLES
  float *particle_dip_device = nullptr;
#endif
#ifdef ELECTROSTATICS
  float *particle_q_device = nullptr;
#endif

  static auto make_shared() {
    auto obj = std::make_shared<GpuParticleData::Storage>();
    System::get_system().cleanup_queue.push<DeviceMemory>(obj);
    return obj;
  }

  ~Storage() { free_device_memory(); }
  void realloc_device_memory();
  void split_particle_struct();
  void copy_particles_to_device();
  void copy_particle_forces_to_host() {
    if (not particle_forces_device.empty()) {
      thrust::copy(particle_forces_device.begin(), particle_forces_device.end(),
                   particle_forces_host.begin());
    }
  }
#ifdef ROTATION
  void copy_particle_torques_to_host() {
    if (not particle_torques_device.empty()) {
      thrust::copy(particle_torques_device.begin(),
                   particle_torques_device.end(),
                   particle_torques_host.begin());
    }
  }
#endif
  Utils::Span<float> get_particle_forces_host_span() {
    return {particle_forces_host.data(), particle_forces_host.size()};
  }
#ifdef ROTATION
  Utils::Span<float> get_particle_torques_host_span() {
    return {particle_torques_host.data(), particle_torques_host.size()};
  }
#endif
};

void GpuParticleData::init() {
  m_data = GpuParticleData::Storage::make_shared();
}

GpuParticleData::~GpuParticleData() {}

std::size_t GpuParticleData::n_particles() const {
  return m_data->particle_data_device.size();
}

float *GpuParticleData::get_particle_positions_device() const {
  return m_data->particle_pos_device;
}

float *GpuParticleData::get_particle_forces_device() const {
  return raw_data_pointer(m_data->particle_forces_device);
}

#ifdef ROTATION
float *GpuParticleData::get_particle_torques_device() const {
  return raw_data_pointer(m_data->particle_torques_device);
}
#endif

#ifdef DIPOLES
float *GpuParticleData::get_particle_dipoles_device() const {
  return m_data->particle_dip_device;
}
#endif

#ifdef ELECTROSTATICS
float *GpuParticleData::get_particle_charges_device() const {
  return m_data->particle_q_device;
}
#endif

GpuParticleData::GpuEnergy *GpuParticleData::get_energy_device() const {
  return m_data->energy_device;
}

void GpuParticleData::enable_property(std::size_t property) {
  m_need_particles_update = true;
  m_data->m_need[property] = true;
  if (property != prop::force and property != prop::torque) {
    m_split_particle_struct = true;
  }
  enable_particle_transfer();
}

bool GpuParticleData::has_compatible_device_impl() const {
  auto result = true;
  try {
    cuda_check_device();
  } catch (cuda_runtime_error const &err) {
    result = false;
  }
  return result;
}

/**
 * @brief Setup and call particle reallocation from the host.
 */
void GpuParticleData::gpu_init_particle_comm() {
  try {
    cuda_check_device();
  } catch (cuda_runtime_error const &err) {
    fprintf(stderr, "ERROR: %s\n", err.what());
    errexit();
  }
  m_data->realloc_device_memory();
}

void GpuParticleData::Storage::copy_particles_to_device() {
  // resize buffers
  auto const n_part = particle_data_host.size();
  resize_or_replace(particle_data_device, n_part);
  particle_forces_host.resize(3ul * n_part);
  resize_or_replace(particle_forces_device, 3ul * n_part);
#ifdef ROTATION
  particle_torques_host.resize(3ul * n_part);
  resize_or_replace(particle_torques_device, 3ul * n_part);
#endif

  // zero out device memory for forces and torques
  cudaMemsetAsync(raw_data_pointer(particle_forces_device), 0x0,
                  byte_size(particle_forces_device), stream[0]);
#ifdef ROTATION
  cudaMemsetAsync(raw_data_pointer(particle_torques_device), 0x0,
                  byte_size(particle_torques_device), stream[0]);
#endif

  // copy particles to device
  cudaMemcpyAsync(raw_data_pointer(particle_data_device),
                  particle_data_host.data(), byte_size(particle_data_host),
                  cudaMemcpyHostToDevice, stream[0]);
}

void GpuParticleData::copy_particles_to_device(ParticleRange const &particles,
                                               int this_node) {
  if (m_communication_enabled) {
    gather_particle_data(particles, m_data->particle_data_host, this_node);
    if (this_node == 0) {
      m_data->copy_particles_to_device();
      if (m_split_particle_struct) {
        m_data->realloc_device_memory();
        m_data->split_particle_struct();
      }
    }
  }
}

void GpuParticleData::copy_forces_to_host(ParticleRange const &particles,
                                          int this_node) {
  if (m_communication_enabled) {
    // copy results from device memory to host memory
    if (this_node == 0) {
      m_data->copy_particle_forces_to_host();
#ifdef ROTATION
      m_data->copy_particle_torques_to_host();
#endif
    }

    auto forces_buffer = m_data->get_particle_forces_host_span();
#ifdef ROTATION
    auto torques_buffer = m_data->get_particle_torques_host_span();
#else
    auto torques_buffer = Utils::Span<float>{nullptr, std::size_t{0ul}};
#endif

    // add forces and torques to the particles
    particles_scatter_forces(particles, forces_buffer, torques_buffer);
  }
}

void GpuParticleData::clear_energy_on_device() {
  if (m_communication_enabled) {
    if (m_data->energy_device == nullptr) {
      cuda_safe_mem(cudaMalloc(&m_data->energy_device, sizeof(GpuEnergy)));
    }
    cuda_safe_mem(cudaMemset(m_data->energy_device, 0, sizeof(GpuEnergy)));
  }
}

GpuParticleData::GpuEnergy GpuParticleData::copy_energy_to_host() const {
  GpuEnergy energy_host{};
  if (m_communication_enabled) {
    cuda_safe_mem(cudaMemcpy(&energy_host, m_data->energy_device,
                             sizeof(GpuEnergy), cudaMemcpyDeviceToHost));
  }
  return energy_host;
}

// Position only
__global__ void split_kernel_r(GpuParticleData::GpuParticle *particles,
                               float *r, std::size_t n) {
  auto idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= n)
    return;

  auto const &p = particles[idx];
  idx *= 3u;
  r[idx + 0u] = p.p[0u];
  r[idx + 1u] = p.p[1u];
  r[idx + 2u] = p.p[2u];
}

#ifdef ELECTROSTATICS
// Position and charge
__global__ void split_kernel_rq(GpuParticleData::GpuParticle *particles,
                                float *r, float *q, std::size_t n) {
  auto const idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= n)
    return;

  auto const &p = particles[idx];
  r[3u * idx + 0u] = p.p[0u];
  r[3u * idx + 1u] = p.p[1u];
  r[3u * idx + 2u] = p.p[2u];
  q[idx] = p.q;
}

// Charge only
__global__ void split_kernel_q(GpuParticleData::GpuParticle *particles,
                               float *q, std::size_t n) {
  auto const idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= n)
    return;

  auto const &p = particles[idx];
  q[idx] = p.q;
}
#endif

#ifdef DIPOLES
// Dipole moment
__global__ void split_kernel_dip(GpuParticleData::GpuParticle *particles,
                                 float *dip, std::size_t n) {
  auto idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= n)
    return;

  auto const &p = particles[idx];

  idx *= 3u;

  dip[idx + 0u] = p.dip[0u];
  dip[idx + 1u] = p.dip[1u];
  dip[idx + 2u] = p.dip[2u];
}
#endif

void GpuParticleData::Storage::split_particle_struct() {
  auto const n_part = particle_data_device.size();
  if (n_part == 0ul)
    return;

  using prop = GpuParticleData::prop;
  dim3 const threadsPerBlock{512u, 1u, 1u};
  dim3 const numBlocks{static_cast<unsigned>(n_part / threadsPerBlock.x + 1ul)};

#ifdef ELECTROSTATICS
  if (m_need[prop::q] and m_need[prop::pos]) {
    split_kernel_rq<<<numBlocks, threadsPerBlock, 0, nullptr>>>(
        raw_data_pointer(particle_data_device), particle_pos_device,
        particle_q_device, n_part);
  } else if (m_need[prop::q]) {
    split_kernel_q<<<numBlocks, threadsPerBlock, 0, nullptr>>>(
        raw_data_pointer(particle_data_device), particle_q_device, n_part);
  } else
#endif
      if (m_need[prop::pos]) {
    split_kernel_r<<<numBlocks, threadsPerBlock, 0, nullptr>>>(
        raw_data_pointer(particle_data_device), particle_pos_device, n_part);
  }
#ifdef DIPOLES
  if (m_need[prop::dip]) {
    split_kernel_dip<<<numBlocks, threadsPerBlock, 0, nullptr>>>(
        raw_data_pointer(particle_data_device), particle_dip_device, n_part);
  }
#endif
}

void GpuParticleData::Storage::realloc_device_memory() {
  using prop = GpuParticleData::prop;
  auto const new_size = particle_data_device.size();
  auto const resize_needed = new_size != current_size;
  if (m_need[prop::pos] and (resize_needed or particle_pos_device == nullptr)) {
    if (particle_pos_device != nullptr) {
      cuda_safe_mem(cudaFree(particle_pos_device));
    }
    cuda_safe_mem(
        cudaMalloc(&particle_pos_device, 3ul * new_size * sizeof(float)));
  }
#ifdef DIPOLES
  if (m_need[prop::dip] and (resize_needed or particle_dip_device == nullptr)) {
    if (particle_dip_device != nullptr) {
      cuda_safe_mem(cudaFree(particle_dip_device));
    }
    cuda_safe_mem(
        cudaMalloc(&particle_dip_device, 3ul * new_size * sizeof(float)));
  }
#endif
#ifdef ELECTROSTATICS
  if (m_need[prop::q] and (resize_needed or particle_q_device == nullptr)) {
    if (particle_q_device != nullptr) {
      cuda_safe_mem(cudaFree(particle_q_device));
    }
    cuda_safe_mem(cudaMalloc(&particle_q_device, new_size * sizeof(float)));
  }
#endif
  current_size = new_size;
}

void GpuParticleData::Storage::free_device_memory() {
  auto const free_device_pointer = [](float *&ptr) {
    if (ptr != nullptr) {
      cuda_safe_mem(cudaFree(reinterpret_cast<void *>(ptr)));
      ptr = nullptr;
    }
  };
  free_device_vector(particle_data_device);
  free_device_vector(particle_forces_device);
#ifdef ROTATION
  free_device_vector(particle_torques_device);
#endif
  free_device_pointer(particle_pos_device);
#ifdef DIPOLES
  free_device_pointer(particle_dip_device);
#endif
#ifdef ELECTROSTATICS
  free_device_pointer(particle_q_device);
#endif
}
