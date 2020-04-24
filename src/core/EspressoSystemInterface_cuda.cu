/*
 * Copyright (C) 2014-2019 The ESPResSo project
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

#include "EspressoSystemInterface.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "errorhandling.hpp"

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

// These functions will split the particle data structure into individual arrays
// for each property

// Position and charge
__global__ void split_kernel_rq(CUDA_particle_data *particles, float *r,
                                float *q, int n) {
  auto idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
  if (idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  r[3 * idx + 0] = p.p[0];
  r[3 * idx + 1] = p.p[1];
  r[3 * idx + 2] = p.p[2];
#ifdef ELECTROSTATICS
  q[idx] = p.q;
#endif
}

// Charge only
__global__ void split_kernel_q(CUDA_particle_data *particles, float *q, int n) {
  auto idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
  if (idx >= n)
    return;

#ifdef ELECTROSTATICS
  CUDA_particle_data p = particles[idx];

  q[idx] = p.q;
#endif
}

// Position only
__global__ void split_kernel_r(CUDA_particle_data *particles, float *r, int n) {
  auto idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
  if (idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  r[idx + 0] = p.p[0];
  r[idx + 1] = p.p[1];
  r[idx + 2] = p.p[2];
}

#ifdef CUDA
// Velocity
__global__ void split_kernel_v(CUDA_particle_data *particles, float *v, int n) {
  auto idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
  if (idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  v[idx + 0] = p.v[0];
  v[idx + 1] = p.v[1];
  v[idx + 2] = p.v[2];
}
#endif

#ifdef DIPOLES
// Dipole moment
__global__ void split_kernel_dip(CUDA_particle_data *particles, float *dip,
                                 int n) {
  auto idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
  if (idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  dip[idx + 0] = p.dip[0];
  dip[idx + 1] = p.dip[1];
  dip[idx + 2] = p.dip[2];
}
#endif

__global__ void split_kernel_director(CUDA_particle_data *particles,
                                      float *director, int n) {
#ifdef ROTATION
  auto idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
  if (idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  director[idx + 0] = p.director[0];
  director[idx + 1] = p.director[1];
  director[idx + 2] = p.director[2];
#endif
}

void EspressoSystemInterface::reallocDeviceMemory(int n) {
  if (m_needsRGpu && ((n != m_gpu_npart) || (m_r_gpu_begin == nullptr))) {
    if (m_r_gpu_begin != nullptr)
      cuda_safe_mem(cudaFree(m_r_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_r_gpu_begin, 3 * n * sizeof(float)));
    m_r_gpu_end = m_r_gpu_begin + 3 * n;
  }
#ifdef DIPOLES
  if (m_needsDipGpu && ((n != m_gpu_npart) || (m_dip_gpu_begin == nullptr))) {
    if (m_dip_gpu_begin != nullptr)
      cuda_safe_mem(cudaFree(m_dip_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_dip_gpu_begin, 3 * n * sizeof(float)));
    m_dip_gpu_end = m_dip_gpu_begin + 3 * n;
  }
#endif
  if (m_needsVGpu && ((n != m_gpu_npart) || (m_v_gpu_begin == nullptr))) {
    if (m_v_gpu_begin != nullptr)
      cuda_safe_mem(cudaFree(m_v_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_v_gpu_begin, 3 * n * sizeof(float)));
    m_v_gpu_end = m_v_gpu_begin + 3 * n;
  }

  if (m_needsQGpu && ((n != m_gpu_npart) || (m_q_gpu_begin == nullptr))) {
    if (m_q_gpu_begin != nullptr)
      cuda_safe_mem(cudaFree(m_q_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_q_gpu_begin, 3 * n * sizeof(float)));
    m_q_gpu_end = m_q_gpu_begin + 3 * n;
  }

  if (m_needsDirectorGpu &&
      ((n != m_gpu_npart) || (m_director_gpu_begin == nullptr))) {
    if (m_director_gpu_begin != nullptr)
      cuda_safe_mem(cudaFree(m_director_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_director_gpu_begin, 3 * n * sizeof(float)));
    m_director_gpu_end = m_director_gpu_begin + 3 * n;
  }

  m_gpu_npart = n;
}

void EspressoSystemInterface::split_particle_struct() {
  auto device_particles = gpu_get_particle_pointer();
  int n = device_particles.size();
  if (n == 0)
    return;

  dim3 grid(n / 512 + 1, 1, 1);
  dim3 block(512, 1, 1);

  if (m_needsQGpu && m_needsRGpu)
    hipLaunchKernelGGL(split_kernel_rq, dim3(grid), dim3(block), 0, nullptr,
                       device_particles.data(), m_r_gpu_begin, m_q_gpu_begin,
                       n);
  if (m_needsQGpu && !m_needsRGpu)
    hipLaunchKernelGGL(split_kernel_q, dim3(grid), dim3(block), 0, nullptr,
                       device_particles.data(), m_q_gpu_begin, n);
  if (!m_needsQGpu && m_needsRGpu)
    hipLaunchKernelGGL(split_kernel_r, dim3(grid), dim3(block), 0, nullptr,
                       device_particles.data(), m_r_gpu_begin, n);
#ifdef CUDA
  if (m_needsVGpu)
    hipLaunchKernelGGL(split_kernel_v, dim3(grid), dim3(block), 0, nullptr,
                       device_particles.data(), m_v_gpu_begin, n);
#endif
#ifdef DIPOLES
  if (m_needsDipGpu)
    hipLaunchKernelGGL(split_kernel_dip, dim3(grid), dim3(block), 0, nullptr,
                       device_particles.data(), m_dip_gpu_begin, n);

#endif

  if (m_needsDirectorGpu)
    hipLaunchKernelGGL(split_kernel_director, dim3(grid), dim3(block), 0,
                       nullptr, device_particles.data(), m_director_gpu_begin,
                       n);
}
