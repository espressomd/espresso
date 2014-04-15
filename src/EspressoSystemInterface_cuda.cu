
#include "EspressoSystemInterface.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"

__global__ void split_kernel_rq(CUDA_particle_data *particles, float *r, float *q, int n) {
  int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx >= n)
    return;

  CUDA_particle_data p = particles[idx];

  r[idx + 0] = p.p[0];
  r[idx + 1] = p.p[1];
  r[idx + 2] = p.p[2];
  #ifdef ELECTROSTATICS
  q[idx] = p.q;
  #endif
}

__global__ void split_kernel_q(CUDA_particle_data *particles,float *q, int n) {
  int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx >= n)
    return;

#ifdef ELECTROSTRATICS
  CUDA_particle_data p = particles[idx];

  q[idx] = p.q;
#endif
}

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

void EspressoSystemInterface::reallocDeviceMemory(int n) {

  if(m_needsRGpu && ((n != m_gpu_npart) || (m_r_gpu_begin == 0))) {
    if(m_r_gpu_begin != 0)
      cuda_safe_mem(cudaFree(m_r_gpu_begin));
    cuda_safe_mem(cudaMalloc(&m_r_gpu_begin, 3*n*sizeof(float)));
    m_r_gpu_end = m_r_gpu_begin + 3*n;
  }

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
  if(m_needsVGpu)
    split_kernel_v<<<grid,block>>>(gpu_get_particle_pointer(), m_v_gpu_begin,n);
}
