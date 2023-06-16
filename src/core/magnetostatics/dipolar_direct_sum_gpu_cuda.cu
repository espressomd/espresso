/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef DIPOLAR_DIRECT_SUM

#include "magnetostatics/dipolar_direct_sum_gpu_cuda.cuh"

#include "cuda/utils.cuh"

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#include <cuda.h>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

__device__ inline float scalar_product(float const *a, float const *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

__device__ inline void vector_product(float const *a, float const *b,
                                      float *out) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

__device__ inline void get_mi_vector_dds(float res[3], float const a[3],
                                         float const b[3], float const box_l[3],
                                         int const periodic[3]) {
  for (int i = 0; i < 3; i++) {
    res[i] = a[i] - b[i];
    if (periodic[i])
      res[i] -= floor(res[i] / box_l[i] + 0.5f) * box_l[i];
  }
}

__device__ void dipole_ia_force(float pf, float const *r1, float const *r2,
                                float const *dip1, float const *dip2, float *f1,
                                float *torque1, float *torque2, float box_l[3],
                                int periodic[3]) {
  // Distance between particles
  float dr[3];
  get_mi_vector_dds(dr, r1, r2, box_l, periodic);

  // Powers of distance
  auto const r_sq = scalar_product(dr, dr);
  auto const r_sq_inv = 1.0f / r_sq;
  auto const r_inv = rsqrtf(r_sq);
  auto const r3_inv = 1.0f / r_sq * r_inv;
  auto const r5_inv = r3_inv * r_sq_inv;
  auto const r7_inv = r5_inv * r_sq_inv;

  // Dot products
  auto const pe1 = scalar_product(dip1, dip2);
  auto const pe2 = scalar_product(dip1, dr);
  auto const pe3 = scalar_product(dip2, dr);
  auto const pe4 = 3.0f * r5_inv;

  // Force
  auto const aa = pe4 * pe1;
  auto const bb = -15.0f * pe2 * pe3 * r7_inv;
  auto const ab = aa + bb;
  auto const cc = pe4 * pe3;
  auto const dd = pe4 * pe2;

  f1[0] = (pf * (ab * dr[0] + cc * dip1[0] + dd * dip2[0]));
  f1[1] = (pf * (ab * dr[1] + cc * dip1[1] + dd * dip2[1]));
  f1[2] = (pf * (ab * dr[2] + cc * dip1[2] + dd * dip2[2]));

#ifdef ROTATION
  // Torques
  float a[3];
  vector_product(dip1, dip2, a);

  float b[3];
  vector_product(dip1, dr, b);

  torque1[0] = pf * (-a[0] * r3_inv + b[0] * cc);
  torque1[1] = pf * (-a[1] * r3_inv + b[1] * cc);
  torque1[2] = pf * (-a[2] * r3_inv + b[2] * cc);

  vector_product(dip2, dr, b);

  torque2[0] = pf * (a[0] * r3_inv + b[0] * dd);
  torque2[1] = pf * (a[1] * r3_inv + b[1] * dd);
  torque2[2] = pf * (a[2] * r3_inv + b[2] * dd);
#endif
}

__device__ float dipole_ia_energy(float pf, float const *r1, float const *r2,
                                  float const *dip1, float const *dip2,
                                  float box_l[3], int periodic[3]) {
  // Distance between particles
  float dr[3];
  get_mi_vector_dds(dr, r1, r2, box_l, periodic);

  // Powers of distance
  auto const r_sq = scalar_product(dr, dr);
  auto const r_sq_inv = 1.0f / r_sq;
  auto const r_inv = rsqrtf(r_sq);
  auto const r3_inv = 1.0f / r_sq * r_inv;
  auto const r5_inv = r3_inv * r_sq_inv;

  // Dot products
  auto const pe1 = scalar_product(dip1, dip2);
  auto const pe2 = scalar_product(dip1, dr);
  auto const pe3 = scalar_product(dip2, dr);
  auto const pe4 = 3.0f * r5_inv;

  // Energy
  return pf * (pe1 * r3_inv - pe4 * pe2 * pe3);
}

__global__ void DipolarDirectSum_kernel_force(float pf, unsigned int n,
                                              float *pos, float *dip, float *f,
                                              float *torque, float box_l[3],
                                              int periodic[3]) {

  auto const i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= n)
    return;

  // Kahan summation based on the wikipedia article
  // Force
  float fi[3], fsum[3], tj[3];

  // Torque
  float ti[3], tsum[3];

  // There is one thread per particle. Each thread computes interactions
  // with particles whose id is smaller than the thread id.
  // The force and torque of all the interaction partners of the current thread
  // is atomically added to global results at once.
  // The result for the particle id equal to the thread id is atomically added
  // to global memory at the end.

  // Clear summation vars
  for (unsigned int j = 0; j < 3; j++) {
    // Force
    fsum[j] = 0;
    // Torque
    tsum[j] = 0;
  }

  for (unsigned int j = i + 1; j < n; j++) {
    dipole_ia_force(pf, pos + 3 * i, pos + 3 * j, dip + 3 * i, dip + 3 * j, fi,
                    ti, tj, box_l, periodic);
    for (unsigned int k = 0; k < 3; k++) {
      // Add rhs to global memory
      atomicAdd(f + 3 * j + k, -fi[k]);
      atomicAdd((torque + 3 * j + k), tj[k]);
      tsum[k] += ti[k];
      fsum[k] += fi[k];
    }
  }

  // Add the left hand side result to global memory
  for (int j = 0; j < 3; j++) {
    atomicAdd(f + 3 * i + j, fsum[j]);
    atomicAdd(torque + 3 * i + j, tsum[j]);
  }
}

__device__ void dds_sumReduction(float *input, float *sum) {
  auto const tid = static_cast<int>(threadIdx.x);
  for (auto i = static_cast<int>(blockDim.x); i > 1; i /= 2) {
    __syncthreads();
    if (tid < i / 2)
      input[tid] += input[i / 2 + tid];
    if ((i % 2 == 1) && (tid == 0))
      input[tid] += input[i - 1];
  }
  __syncthreads();
  if (tid == 0) {
    sum[0] = input[0];
  }
}

__global__ void DipolarDirectSum_kernel_energy(float pf, unsigned int n,
                                               float *pos, float *dip,
                                               float box_l[3], int periodic[3],
                                               float *energySum) {

  auto const i = blockIdx.x * blockDim.x + threadIdx.x;
  float sum = 0.0;
  extern __shared__ float res[];

  // There is one thread per particle. Each thread computes interactions
  // with particles whose id is larger than the thread id.
  // The result for the particle id equal to the thread id is added
  // to global memory at the end.

  if (i < n) {
    // Summation for particle i
    for (unsigned int j = i + 1; j < n; j++) {
      sum += dipole_ia_energy(pf, pos + 3 * i, pos + 3 * j, dip + 3 * i,
                              dip + 3 * j, box_l, periodic);
    }

    // Save per thread result into block shared mem
    res[threadIdx.x] = sum;
  } else
    res[threadIdx.x] = 0;

  // Sum results within a block
  __syncthreads(); // Wait till all threads in block are done
  dds_sumReduction(res, &(energySum[blockIdx.x]));
}

inline void copy_box_data(float **box_l_gpu, int **periodic_gpu,
                          float const *box_l, int const *periodic) {
  auto const s_box = 3u * sizeof(float);
  auto const s_per = 3u * sizeof(int);
  cuda_safe_mem(cudaMalloc(reinterpret_cast<void **>(box_l_gpu), s_box));
  cuda_safe_mem(cudaMalloc(reinterpret_cast<void **>(periodic_gpu), s_per));
  cuda_safe_mem(cudaMemcpy(*box_l_gpu, box_l, s_box, cudaMemcpyHostToDevice));
  cuda_safe_mem(
      cudaMemcpy(*periodic_gpu, periodic, s_per, cudaMemcpyHostToDevice));
}

void DipolarDirectSum_kernel_wrapper_force(float k, unsigned int n, float *pos,
                                           float *dip, float *f, float *torque,
                                           float box_l[3], int periodic[3]) {

  unsigned int const bs = 64;
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  if (n == 0)
    return;

  if (n <= bs) {
    grid.x = 1;
    block.x = n;
  } else {
    grid.x = n / bs + 1;
    block.x = bs;
  }

  float *box_l_gpu;
  int *periodic_gpu;
  copy_box_data(&box_l_gpu, &periodic_gpu, box_l, periodic);

  KERNELCALL(DipolarDirectSum_kernel_force, grid, block, k, n, pos, dip, f,
             torque, box_l_gpu, periodic_gpu);
  cudaFree(box_l_gpu);
  cudaFree(periodic_gpu);
}

void DipolarDirectSum_kernel_wrapper_energy(float k, unsigned int n, float *pos,
                                            float *dip, float box_l[3],
                                            int periodic[3], float *E) {

  unsigned int const bs = 512;
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  if (n == 0)
    return;

  if (n <= bs) {
    grid.x = 1;
    block.x = n;
  } else {
    grid.x = n / bs + 1;
    block.x = bs;
  }

  float *box_l_gpu;
  int *periodic_gpu;
  copy_box_data(&box_l_gpu, &periodic_gpu, box_l, periodic);

  float *energySum;
  cuda_safe_mem(cudaMalloc(&energySum, sizeof(float) * grid.x));

  // This will sum the energies up to the block level
  KERNELCALL_shared(DipolarDirectSum_kernel_energy, grid, block,
                    bs * sizeof(float), k, n, pos, dip, box_l_gpu, periodic_gpu,
                    energySum);

  // Sum the results of all blocks
  // One thread per block in the prev kernel
  // KERNELCALL(sumKernel,1,1,energySum,block.x,E);
  thrust::device_ptr<float> t(energySum);
  float x = thrust::reduce(t, t + grid.x);
  cuda_safe_mem(cudaMemcpy(E, &x, sizeof(float), cudaMemcpyHostToDevice));

  cuda_safe_mem(cudaFree(energySum));
  cuda_safe_mem(cudaFree(box_l_gpu));
  cuda_safe_mem(cudaFree(periodic_gpu));
}

#endif
