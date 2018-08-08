/*
  Copyright (C) 2014,2015,2016 The ESPResSo project
  
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

#include "EspressoSystemInterface.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "cuda_init.hpp"
#include "errorhandling.hpp"

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

// These functions will split the paritlce data structure into individual arrays for each property

// Position and charge
__global__ void split_kernel_rq(CUDA_particle_data *particles, float *r, float *q, int n) {
  int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  r[3*idx + 0] = p.p[0];
  r[3*idx + 1] = p.p[1];
  r[3*idx + 2] = p.p[2];
  #ifdef ELECTROSTATICS
  q[idx] = p.q;
  #endif
}

// Charge only
__global__ void split_kernel_q(CUDA_particle_data *particles,float *q, int n) {
  int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx >= n)
    return;

#ifdef ELECTROSTATICS
  CUDA_particle_data p = particles[idx];

  q[idx] = p.q;
#endif
}

// Position only
__global__ void split_kernel_r(CUDA_particle_data *particles, float *r, int n) {
  int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  r[idx + 0] = p.p[0];
  r[idx + 1] = p.p[1];
  r[idx + 2] = p.p[2];
}

#ifdef LB_GPU
// Velocity
__global__ void split_kernel_v(CUDA_particle_data *particles, float *v, int n) {
  int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx >= n)
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
__global__ void split_kernel_dip(CUDA_particle_data *particles, float *dip, int n) {
  int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  dip[idx + 0] = p.dip[0];
  dip[idx + 1] = p.dip[1];
  dip[idx + 2] = p.dip[2];
}
#endif

__global__ void split_kernel_quatu(CUDA_particle_data *particles, float *quatu, int n) {
#ifdef ROTATION
  int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  idx *= 3;

  quatu[idx + 0] = p.quatu[0];
  quatu[idx + 1] = p.quatu[1];
  quatu[idx + 2] = p.quatu[2];
#endif
}

void EspressoSystemInterface::reallocDeviceMemory(int n) {
  if(m_needsRGpu && ((n != m_gpu_npart) || (m_r_gpu_begin == 0))) {
    if(m_r_gpu_begin != 0)
      cuda_safe_mem(cudaFree(m_r_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_r_gpu_begin, 3*n*sizeof(float)));
    m_r_gpu_end = m_r_gpu_begin + 3*n;
  }
#ifdef DIPOLES
  if(m_needsDipGpu && ((n != m_gpu_npart) || (m_dip_gpu_begin == 0))) {
    if(m_dip_gpu_begin != 0)
      cuda_safe_mem(cudaFree(m_dip_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_dip_gpu_begin, 3*n*sizeof(float)));
    m_dip_gpu_end = m_dip_gpu_begin + 3*n;
  }
#endif
  if(m_needsVGpu && ((n != m_gpu_npart) || (m_v_gpu_begin == 0))) {
    if(m_v_gpu_begin != 0)
      cuda_safe_mem(cudaFree(m_v_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_v_gpu_begin, 3*n*sizeof(float)));
    m_v_gpu_end = m_v_gpu_begin + 3*n;
  }

  if(m_needsQGpu && ((n != m_gpu_npart) || (m_q_gpu_begin == 0))) {
    if(m_q_gpu_begin != 0)
      cuda_safe_mem(cudaFree(m_q_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_q_gpu_begin, 3*n*sizeof(float)));
    m_q_gpu_end = m_q_gpu_begin + 3*n;
  }

  if(m_needsQuatuGpu && ((n != m_gpu_npart) || (m_quatu_gpu_begin == 0))) {
    if(m_quatu_gpu_begin != 0)
      cuda_safe_mem(cudaFree(m_quatu_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_quatu_gpu_begin, 3*n*sizeof(float)));
    m_quatu_gpu_end = m_quatu_gpu_begin + 3*n;
  }

  m_gpu_npart = n;
}

void EspressoSystemInterface::split_particle_struct() {
  int n = gpu_get_global_particle_vars_pointer_host()->number_of_particles;
  if(n == 0) 
    return;

  ESIF_TRACE(printf("n = %d, m_gpu_npart = %d\n", n, m_gpu_npart));
    
  dim3 grid(n/512+1,1,1);
  dim3 block(512,1,1);

  if(m_needsQGpu && m_needsRGpu)
    split_kernel_rq<<<grid,block>>>(gpu_get_particle_pointer(), m_r_gpu_begin,m_q_gpu_begin,n);
  if(m_needsQGpu && !m_needsRGpu)
    split_kernel_q<<<grid,block>>>(gpu_get_particle_pointer(), m_q_gpu_begin,n);
  if(!m_needsQGpu && m_needsRGpu)
    split_kernel_r<<<grid,block>>>(gpu_get_particle_pointer(), m_r_gpu_begin,n);
#ifdef LB_GPU
  if(m_needsVGpu)
    split_kernel_v<<<grid,block>>>(gpu_get_particle_pointer(), m_v_gpu_begin,n);
#endif
#ifdef DIPOLES
  if(m_needsDipGpu)
    split_kernel_dip<<<grid,block>>>(gpu_get_particle_pointer(), m_dip_gpu_begin,n);

#endif

  if(m_needsQuatuGpu)
    split_kernel_quatu<<<grid,block>>>(gpu_get_particle_pointer(), m_quatu_gpu_begin,n);
}
