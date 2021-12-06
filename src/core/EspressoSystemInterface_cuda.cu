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
/**
 * @file
 * CUDA kernels to convert the particles AoS to a SoA on the device.
 */

#include "EspressoSystemInterface.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.cuh"
#include "errorhandling.hpp"

#include <cstddef>

#include <cuda.h>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

// Position and charge
__global__ void split_kernel_rq(CUDA_particle_data *particles, float *r,
                                float *q, unsigned int n) {
  auto const idx = blockDim.x * blockIdx.x + threadIdx.x;
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
__global__ void split_kernel_q(CUDA_particle_data *particles, float *q,
                               unsigned int n) {
  auto const idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= n)
    return;

#ifdef ELECTROSTATICS
  CUDA_particle_data p = particles[idx];

  q[idx] = p.q;
#endif
}

// Position only
__global__ void split_kernel_r(CUDA_particle_data *particles, float *r,
                               unsigned int n) {
  auto idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  r[idx + 0] = p.p[0];
  r[idx + 1] = p.p[1];
  r[idx + 2] = p.p[2];
}

#ifdef DIPOLES
// Dipole moment
__global__ void split_kernel_dip(CUDA_particle_data *particles, float *dip,
                                 unsigned int n) {
  auto idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  dip[idx + 0] = p.dip[0];
  dip[idx + 1] = p.dip[1];
  dip[idx + 2] = p.dip[2];
}
#endif

void EspressoSystemInterface::reallocDeviceMemory(std::size_t n) {
  if (m_needsRGpu && ((n != m_gpu_npart) || (m_r_gpu_begin == nullptr))) {
    if (m_r_gpu_begin != nullptr)
      cuda_safe_mem(cudaFree(m_r_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_r_gpu_begin, 3 * n * sizeof(float)));
  }
#ifdef DIPOLES
  if (m_needsDipGpu && ((n != m_gpu_npart) || (m_dip_gpu_begin == nullptr))) {
    if (m_dip_gpu_begin != nullptr)
      cuda_safe_mem(cudaFree(m_dip_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_dip_gpu_begin, 3 * n * sizeof(float)));
  }
#endif
#ifdef ELECTROSTATICS
  if (m_needsQGpu && ((n != m_gpu_npart) || (m_q_gpu_begin == nullptr))) {
    if (m_q_gpu_begin != nullptr)
      cuda_safe_mem(cudaFree(m_q_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_q_gpu_begin, 3 * n * sizeof(float)));
  }
#endif

  m_gpu_npart = n;
}

void EspressoSystemInterface::split_particle_struct() {
  auto const device_particles = gpu_get_particle_pointer();
  auto const n = static_cast<unsigned int>(device_particles.size());
  if (n == 0)
    return;

  dim3 grid(n / 512 + 1, 1, 1);
  dim3 block(512, 1, 1);

  if (m_needsQGpu && m_needsRGpu)
    split_kernel_rq<<<dim3(grid), dim3(block), 0, nullptr>>>(
        device_particles.data(), m_r_gpu_begin, m_q_gpu_begin, n);
  else if (m_needsQGpu)
    split_kernel_q<<<dim3(grid), dim3(block), 0, nullptr>>>(
        device_particles.data(), m_q_gpu_begin, n);
  else if (m_needsRGpu)
    split_kernel_r<<<dim3(grid), dim3(block), 0, nullptr>>>(
        device_particles.data(), m_r_gpu_begin, n);
#ifdef DIPOLES
  if (m_needsDipGpu)
    split_kernel_dip<<<dim3(grid), dim3(block), 0, nullptr>>>(
        device_particles.data(), m_dip_gpu_begin, n);
#endif
}
