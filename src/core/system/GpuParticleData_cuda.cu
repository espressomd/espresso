/*
 * Copyright (C) 2014-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "GpuParticleData.hpp"
#include "System.hpp"

#include "ParticleRange.hpp"
#include "errorhandling.hpp"

#include "cuda/init.hpp"
#include "cuda/utils.cuh"

#include <thrust/copy.h>
#include <thrust/device_vector.h>

#include <cuda.h>

#include <cstddef>
#include <cstdio>
#include <memory>
#include <span>

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
 * But Clang 20's UBSAN still triggers a runtime error as of CUDA 12.0
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
public:
  /** @brief Which particle properties are needed by GPU methods. */
  GpuParticleData::prop::bitset m_need;
  GpuParticleData::GpuEnergy *energy_device = nullptr;
  std::size_t current_size = 0ul;
  // Authoritative local particle count for n_particles() (consumed by the
  // P3M/DDS solvers). Set on both transfer paths; see copy_particles_to_device.
  std::size_t n_particles = 0ul;
  pinned_vector<GpuParticle> particle_data_host;
  thrust::device_vector<GpuParticle> particle_data_device;
  pinned_vector<float> particle_forces_host;
  thrust::device_vector<float> particle_forces_device;
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  pinned_vector<float> particle_dip_fld_host;
  thrust::device_vector<float> particle_dip_fld_device;
#endif
#ifdef ESPRESSO_ROTATION
  pinned_vector<float> particle_torques_host;
  thrust::device_vector<float> particle_torques_device;
#endif
  float *particle_pos_device = nullptr;
#ifdef ESPRESSO_DIPOLES
  float *particle_dip_device = nullptr;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  float *particle_q_device = nullptr;
#endif

  // -- single-rank per-field SoA staging buffers ---------------------------
  // Pinned host buffers in SoA layout (particle-major float3 for pos/dip, one
  // float per particle for q), filled directly from the ParticleStore columns
  // on the comm.size()==1 fast path. Their contents are byte-identical to what
  // the device split kernels produce, so they are cudaMemcpy'd straight into
  // the SoA device buffers, bypassing the AoS host pack + device AoS buffer +
  // split_kernel_* de-interleave (see copy_particles_soa_to_device below).
  pinned_vector<float> particle_pos_host_soa;
#ifdef ESPRESSO_DIPOLES
  pinned_vector<float> particle_dip_host_soa;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  pinned_vector<float> particle_q_host_soa;
#endif

  ~Storage();
  void realloc_device_memory(std::size_t n_part);
  void split_particle_struct();
  void resize_and_zero_return_buffers(std::size_t n_part);
  void resize_soa_staging_buffers(std::size_t n_part);
  void copy_particles_soa_to_device();
  std::span<float> get_particle_pos_host_soa_span() {
    return {particle_pos_host_soa.data(), particle_pos_host_soa.size()};
  }
#ifdef ESPRESSO_DIPOLES
  std::span<float> get_particle_dip_host_soa_span() {
    return {particle_dip_host_soa.data(), particle_dip_host_soa.size()};
  }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  std::span<float> get_particle_q_host_soa_span() {
    return {particle_q_host_soa.data(), particle_q_host_soa.size()};
  }
#endif
  void copy_particles_to_device();
  void copy_particle_forces_to_host() {
    if (not particle_forces_device.empty()) {
      thrust::copy(particle_forces_device.begin(), particle_forces_device.end(),
                   particle_forces_host.begin());
    }
  }

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  void copy_particle_dip_fld_to_host() {
    if (not particle_dip_fld_device.empty()) {
      thrust::copy(particle_dip_fld_device.begin(),
                   particle_dip_fld_device.end(),
                   particle_dip_fld_host.begin());
    }
  }
#endif

#ifdef ESPRESSO_ROTATION
  void copy_particle_torques_to_host() {
    if (not particle_torques_device.empty()) {
      thrust::copy(particle_torques_device.begin(),
                   particle_torques_device.end(),
                   particle_torques_host.begin());
    }
  }
#endif
  std::span<float> get_particle_forces_host_span() {
    return {particle_forces_host.data(), particle_forces_host.size()};
  }

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  std::span<float> get_particle_dip_fld_host_span() {
    return {particle_dip_fld_host.data(), particle_dip_fld_host.size()};
  }
#endif

#ifdef ESPRESSO_ROTATION
  std::span<float> get_particle_torques_host_span() {
    return {particle_torques_host.data(), particle_torques_host.size()};
  }
#endif
};

// default ctor/dtor definitions must appear out-of-line due to the
// forward-declaration of the Storage class
GpuParticleData::GpuParticleData() = default;
GpuParticleData::~GpuParticleData() = default;

void GpuParticleData::initialize() {
  m_data = std::make_unique<GpuParticleData::Storage>();
  get_system().cleanup_queue.push<DeviceMemory>(shared_from_this());
}

void GpuParticleData::deinitialize() noexcept { m_data.reset(); }

std::size_t GpuParticleData::n_particles() const { return m_data->n_particles; }

float *GpuParticleData::get_particle_positions_device() const {
  return m_data->particle_pos_device;
}

float *GpuParticleData::get_particle_forces_device() const {
  return raw_data_pointer(m_data->particle_forces_device);
}
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
float *GpuParticleData::get_particle_dip_fld_device() const {
  return raw_data_pointer(m_data->particle_dip_fld_device);
}
#endif

#ifdef ESPRESSO_ROTATION
float *GpuParticleData::get_particle_torques_device() const {
  return raw_data_pointer(m_data->particle_torques_device);
}
#endif

#ifdef ESPRESSO_DIPOLES
float *GpuParticleData::get_particle_dipoles_device() const {
  return m_data->particle_dip_device;
}
#endif

#ifdef ESPRESSO_ELECTROSTATICS
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
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (property != prop::force and property != prop::torque and
      property != prop::dip_fld) {
    m_split_particle_struct = true;
  }
#else
  if (property != prop::force and property != prop::torque) {
    m_split_particle_struct = true;
  }
#endif
  enable_particle_transfer();
}

bool GpuParticleData::has_compatible_device_impl() const {
  auto result = false;
  invoke_skip_cuda_exceptions([&result]() {
    cuda_check_device();
    result = true;
  });
  return result;
}

/**
 * @brief Setup and call particle reallocation from the host.
 */
void GpuParticleData::gpu_init_particle_comm() {
  try {
    cuda_check_device();
  } catch (cuda_runtime_error const &err) {
    throw cuda_fatal_error(err.what());
  }
  m_data->realloc_device_memory(m_data->n_particles);
}

void GpuParticleData::Storage::resize_and_zero_return_buffers(
    std::size_t const n_part) {
  // Size the force/torque/dip_fld return-path buffers and zero the device side
  // (the GPU solvers accumulate into these). Shared by the AoS split path and
  // the single-rank SoA staging path.
  particle_forces_host.resize(3ul * n_part);
  resize_or_replace(particle_forces_device, 3ul * n_part);

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  particle_dip_fld_host.resize(3ul * n_part);
  resize_or_replace(particle_dip_fld_device, 3ul * n_part);
#endif

#ifdef ESPRESSO_ROTATION
  particle_torques_host.resize(3ul * n_part);
  resize_or_replace(particle_torques_device, 3ul * n_part);
#endif

  // zero out device memory for forces and torques
  cudaMemsetAsync(raw_data_pointer(particle_forces_device), 0x0,
                  byte_size(particle_forces_device), stream[0]);

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  cudaMemsetAsync(raw_data_pointer(particle_dip_fld_device), 0x0,
                  byte_size(particle_dip_fld_device), stream[0]);
#endif

#ifdef ESPRESSO_ROTATION
  cudaMemsetAsync(raw_data_pointer(particle_torques_device), 0x0,
                  byte_size(particle_torques_device), stream[0]);
#endif
}

void GpuParticleData::Storage::copy_particles_to_device() {
  // resize buffers
  auto const n_part = particle_data_host.size();
  resize_or_replace(particle_data_device, n_part);
  resize_and_zero_return_buffers(n_part);

  // copy particles to device
  cudaMemcpyAsync(raw_data_pointer(particle_data_device),
                  particle_data_host.data(), byte_size(particle_data_host),
                  cudaMemcpyHostToDevice, stream[0]);
}

/**
 * @brief Size the pinned per-field SoA staging buffers for @p n_part particles
 * (single-rank fast path). Only the buffers for enabled properties are grown;
 * the others stay empty (and are skipped by the copy).
 */
void GpuParticleData::Storage::resize_soa_staging_buffers(
    std::size_t const n_part) {
  using prop = GpuParticleData::prop;
  if (m_need[prop::pos]) {
    particle_pos_host_soa.resize(3ul * n_part);
  }
#ifdef ESPRESSO_ELECTROSTATICS
  if (m_need[prop::q]) {
    particle_q_host_soa.resize(n_part);
  }
#endif
#ifdef ESPRESSO_DIPOLES
  if (m_need[prop::dip]) {
    particle_dip_host_soa.resize(3ul * n_part);
  }
#endif
}

/**
 * @brief Copy the pinned per-field SoA staging buffers straight into the SoA
 * device buffers (single-rank fast path).
 *
 * The staging buffers are filled directly from the ParticleStore columns in
 * SoA layout with the same f64->f32 casts and per-field element order that
 * @ref split_particle_struct produces from the AoS pack. Copying them here
 * therefore leaves the device SoA buffers (@ref particle_pos_device,
 * @ref particle_q_device, @ref particle_dip_device) bit-identical to the AoS
 * pack + split path, but without the device AoS buffer or the split kernels.
 */
void GpuParticleData::Storage::copy_particles_soa_to_device() {
  using prop = GpuParticleData::prop;
  if (m_need[prop::pos] and particle_pos_device != nullptr) {
    cudaMemcpyAsync(particle_pos_device, particle_pos_host_soa.data(),
                    byte_size(particle_pos_host_soa), cudaMemcpyHostToDevice,
                    stream[0]);
  }
#ifdef ESPRESSO_ELECTROSTATICS
  if (m_need[prop::q] and particle_q_device != nullptr) {
    cudaMemcpyAsync(particle_q_device, particle_q_host_soa.data(),
                    byte_size(particle_q_host_soa), cudaMemcpyHostToDevice,
                    stream[0]);
  }
#endif
#ifdef ESPRESSO_DIPOLES
  if (m_need[prop::dip] and particle_dip_device != nullptr) {
    cudaMemcpyAsync(particle_dip_device, particle_dip_host_soa.data(),
                    byte_size(particle_dip_host_soa), cudaMemcpyHostToDevice,
                    stream[0]);
  }
#endif
}

void GpuParticleData::copy_particles_to_device(ParticleRange const &particles,
                                               int this_node,
                                               bool single_rank) {
  if (not m_communication_enabled) {
    return;
  }

  // Single-rank fast path: fill per-field SoA host staging buffers directly
  // from the ParticleStore columns and cudaMemcpy each enabled field straight
  // into the device SoA buffers, bypassing the AoS host pack, the device AoS
  // buffer and the split kernels. This is bit-identical to the AoS+split path
  // (same f64->f32 casts, same per-field SoA layout / order; see
  // pack_particles_soa). The only split-relevant properties are pos (always
  // column-backed), q (ELECTROSTATICS column) and dip (DIPOLES columns via
  // calc_dip). On a single rank this_node == 0.
  // force/torque/dip_fld are return-path arrays, unaffected.
  if (single_rank and m_split_particle_struct) {
    auto const n_part = particles.size();
    // n_particles() reads m_data->n_particles (a plain counter); the device
    // AoS buffer is entirely unused on this path (no AoS host copy, no split
    // kernel).
    m_data->n_particles = n_part;
    m_data->resize_and_zero_return_buffers(n_part);
    m_data->realloc_device_memory(n_part);
    m_data->resize_soa_staging_buffers(n_part);
#ifdef ESPRESSO_ELECTROSTATICS
    auto q_span = m_data->get_particle_q_host_soa_span();
#else
    auto q_span = std::span<float>();
#endif
#ifdef ESPRESSO_DIPOLES
    auto dip_span = m_data->get_particle_dip_host_soa_span();
#else
    auto dip_span = std::span<float>();
#endif
    pack_particles_soa(particles, m_data->get_particle_pos_host_soa_span(),
                       q_span, dip_span);
    m_data->copy_particles_soa_to_device();
    return;
  }

  // Multi-rank (or no split needed): head-node MPI gather + AoS pack + split
  // path.
  gather_particle_data(particles, m_data->particle_data_host, this_node);
  if (this_node == 0) {
    m_data->copy_particles_to_device();
    m_data->n_particles = m_data->particle_data_device.size();
    if (m_split_particle_struct) {
      m_data->realloc_device_memory(m_data->n_particles);
      m_data->split_particle_struct();
    }
  }
}

void GpuParticleData::copy_forces_to_host(ParticleRange const &particles,
                                          int this_node) {
  if (m_communication_enabled) {
    // copy results from device memory to host memory
    if (this_node == 0) {
      m_data->copy_particle_forces_to_host();
#ifdef ESPRESSO_ROTATION
      m_data->copy_particle_torques_to_host();
#endif
    }

    auto forces_buffer = m_data->get_particle_forces_host_span();
#ifdef ESPRESSO_ROTATION
    auto torques_buffer = m_data->get_particle_torques_host_span();
#else
    auto torques_buffer = std::span<float>();
#endif

    // add forces and torques to the particles
    particles_scatter_forces(particles, forces_buffer, torques_buffer);
  }
}
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
void GpuParticleData::copy_dip_fld_to_host(ParticleRange const &particles,
                                           int this_node) {
  if (m_communication_enabled) {
    // copy results from device memory to host memory
    if (this_node == 0) {
      m_data->copy_particle_dip_fld_to_host();
    }

    auto dipole_field_buffer = m_data->get_particle_dip_fld_host_span();

    // add dip_fld to the particles
    particles_scatter_dip_fld(particles, dipole_field_buffer);
  }
}
#endif

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

#ifdef ESPRESSO_ELECTROSTATICS
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

#ifdef ESPRESSO_DIPOLES
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

#ifdef ESPRESSO_ELECTROSTATICS
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
#ifdef ESPRESSO_DIPOLES
  if (m_need[prop::dip]) {
    split_kernel_dip<<<numBlocks, threadsPerBlock, 0, nullptr>>>(
        raw_data_pointer(particle_data_device), particle_dip_device, n_part);
  }
#endif
}

void GpuParticleData::Storage::realloc_device_memory(std::size_t const n_part) {
  using prop = GpuParticleData::prop;
  auto const new_size = n_part;
  auto const resize_needed = new_size != current_size;
  if (m_need[prop::pos] and (resize_needed or particle_pos_device == nullptr)) {
    if (particle_pos_device != nullptr) {
      cuda_safe_mem(cudaFree(particle_pos_device));
    }
    cuda_safe_mem(
        cudaMalloc(&particle_pos_device, 3ul * new_size * sizeof(float)));
  }
#ifdef ESPRESSO_DIPOLES
  if (m_need[prop::dip] and (resize_needed or particle_dip_device == nullptr)) {
    if (particle_dip_device != nullptr) {
      cuda_safe_mem(cudaFree(particle_dip_device));
    }
    cuda_safe_mem(
        cudaMalloc(&particle_dip_device, 3ul * new_size * sizeof(float)));
  }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  if (m_need[prop::q] and (resize_needed or particle_q_device == nullptr)) {
    if (particle_q_device != nullptr) {
      cuda_safe_mem(cudaFree(particle_q_device));
    }
    cuda_safe_mem(cudaMalloc(&particle_q_device, new_size * sizeof(float)));
  }
#endif
  current_size = new_size;
}

GpuParticleData::Storage::~Storage() {
  auto const free_device_pointer = [](auto *&ptr) noexcept {
    if (ptr != nullptr) {
      cudaFree(reinterpret_cast<void *>(ptr));
      ptr = nullptr;
    }
  };
  free_device_pointer(particle_pos_device);
#ifdef ESPRESSO_DIPOLES
  free_device_pointer(particle_dip_device);
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  free_device_pointer(particle_q_device);
#endif
  free_device_pointer(energy_device);
}
