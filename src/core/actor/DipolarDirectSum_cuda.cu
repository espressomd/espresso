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
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#ifdef DIPOLAR_DIRECT_SUM

#include "cuda_utils.hpp"
#include <cstdio>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

// This needs to be done in the .cu file too!!!!!
typedef float dds_float;

__device__ inline void get_mi_vector_dds(dds_float res[3], dds_float const a[3],
                                         dds_float const b[3],
                                         dds_float const box_l[3],
                                         int const periodic[3]) {
  for (int i = 0; i < 3; i++) {
    res[i] = a[i] - b[i];
    if (periodic[i])
      res[i] -= floor(res[i] / box_l[i] + 0.5) * box_l[i];
  }
}

#define scalar(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

__device__ void dipole_ia_force(int id, dds_float pf, float const *r1,
                                float const *r2, float const *dip1,
                                float const *dip2, dds_float *f1,
                                dds_float *torque1, dds_float *torque2,
                                dds_float box_l[3], int periodic[3]) {
  dds_float r_inv, pe1, pe2, pe3, pe4, r_sq, r3_inv, r5_inv, r_sq_inv, r7_inv,
      a, b, cc, d, ab;
#ifdef ROTATION
  dds_float bx, by, bz, ax, ay, az;
#endif

  dds_float dr[3];
  dds_float _r1[3], _r2[3], _dip1[3], _dip2[3];

  for (int i = 0; i < 3; i++) {
    _r1[i] = r1[i];
    _r2[i] = r2[i];
    _dip1[i] = dip1[i];
    _dip2[i] = dip2[i];
  }

  // Distance between particles
  get_mi_vector_dds(dr, _r1, _r2, box_l, periodic);

  // Powers of distance
  r_sq = scalar(dr, dr);
  r_sq_inv = 1 / r_sq;
  r_inv = rsqrtf(r_sq);
  r3_inv = 1 / r_sq * r_inv;
  r5_inv = r3_inv * r_sq_inv;
  r7_inv = r5_inv * r_sq_inv;

  // Dot products
  pe1 = scalar(_dip1, _dip2);
  pe2 = scalar(_dip1, dr);
  pe3 = scalar(_dip2, dr);
  pe4 = 3.0f * r5_inv;

  // Force
  a = pe4 * pe1;
  b = -15.0f * pe2 * pe3 * r7_inv;
  ab = a + b;
  cc = pe4 * pe3;
  d = pe4 * pe2;

  f1[0] = (pf * (ab * dr[0] + cc * _dip1[0] + d * _dip2[0]));
  f1[1] = (pf * (ab * dr[1] + cc * _dip1[1] + d * _dip2[1]));
  f1[2] = (pf * (ab * dr[2] + cc * _dip1[2] + d * _dip2[2]));

#ifdef ROTATION
  // Torques
  ax = _dip1[1] * _dip2[2] - _dip2[1] * _dip1[2];
  ay = _dip2[0] * _dip1[2] - _dip1[0] * _dip2[2];
  az = _dip1[0] * _dip2[1] - _dip2[0] * _dip1[1];

  bx = _dip1[1] * dr[2] - dr[1] * _dip1[2];
  by = dr[0] * _dip1[2] - _dip1[0] * dr[2];
  bz = _dip1[0] * dr[1] - dr[0] * _dip1[1];

  torque1[0] = (pf * (-ax * r3_inv + bx * cc));
  torque1[1] = (pf * (-ay * r3_inv + by * cc));
  torque1[2] = (pf * (-az * r3_inv + bz * cc));

  bx = _dip2[1] * dr[2] - dr[1] * _dip2[2];
  by = dr[0] * _dip2[2] - _dip2[0] * dr[2];
  bz = _dip2[0] * dr[1] - dr[0] * _dip2[1];

  torque2[0] = pf * (ax * r3_inv + bx * d);
  torque2[1] = pf * (ay * r3_inv + by * d);
  torque2[2] = pf * (az * r3_inv + bz * d);
#endif
}

__device__ dds_float dipole_ia_energy(int id, dds_float pf, float const *r1,
                                      float const *r2, float const *dip1,
                                      float const *dip2, dds_float box_l[3],
                                      int periodic[3]) {
  dds_float r_inv, pe1, pe2, pe3, pe4, r_sq, r3_inv, r5_inv, r_sq_inv;
  dds_float dr[3];
  dds_float _r1[3], _r2[3], _dip1[3], _dip2[3];

  for (int i = 0; i < 3; i++) {
    _r1[i] = r1[i];
    _r2[i] = r2[i];
    _dip1[i] = dip1[i];
    _dip2[i] = dip2[i];
  }

  // Distance between particles
  get_mi_vector_dds(dr, _r1, _r2, box_l, periodic);

  // Powers of distance
  r_sq = scalar(dr, dr);
  r_sq_inv = 1 / r_sq;
  r_inv = rsqrtf(r_sq);
  r3_inv = 1 / r_sq * r_inv;
  r5_inv = r3_inv * r_sq_inv;

  // Dot products
  pe1 = scalar(_dip1, _dip2);
  pe2 = scalar(_dip1, dr);
  pe3 = scalar(_dip2, dr);
  pe4 = 3.0f * r5_inv;

  // Energy
  return pf * (pe1 * r3_inv - pe4 * pe2 * pe3);
}

__global__ void DipolarDirectSum_kernel_force(dds_float pf, int n, float *pos,
                                              float *dip, float *f,
                                              float *torque, dds_float box_l[3],
                                              int periodic[3]) {

  auto i = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);

  if (i >= n)
    return;

  // Kahan summation based on the wikipedia article
  // Force
  dds_float fi[3], fsum[3], tj[3];

  // Torque
  dds_float ti[3], tsum[3];

  // There is one thread per particle. Each thread computes interactions
  // with particles whose id is smaller than the thread id.
  // The force and torque of all the interaction partners of the current thread
  // is atomically added to global results ad once.
  // The result for the particle id equal to the thread id is atomically added
  // to global memory at the end.

  // Clear summation vars
  for (int j = 0; j < 3; j++) {
    // Force
    fsum[j] = 0;
    // Torque
    tsum[j] = 0;
  }

  for (int j = i + 1; j < n; j++) {
    dipole_ia_force(i, pf, pos + 3 * i, pos + 3 * j, dip + 3 * i, dip + 3 * j,
                    fi, ti, tj, box_l, periodic);
    for (int k = 0; k < 3; k++) {
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

__device__ void dds_sumReduction(dds_float *input, dds_float *sum) {
  auto tid = static_cast<int>(threadIdx.x);
  for (auto i = static_cast<int>(blockDim.x); i > 1; i /= 2) {
    __syncthreads();
    if (tid < i / 2)
      input[tid] += input[i / 2 + tid];
    if ((i % 2 == 1) && (tid == 0))
      input[tid] += input[i - 1];
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    sum[0] = input[0];
  }
}

__global__ void DipolarDirectSum_kernel_energy(dds_float pf, int n, float *pos,
                                               float *dip, dds_float box_l[3],
                                               int periodic[3],
                                               dds_float *energySum) {

  auto i = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  dds_float sum = 0.0;
  HIP_DYNAMIC_SHARED(dds_float, res)

  // There is one thread per particle. Each thread computes interactions
  // with particles whose id is larger than the thread id.
  // The result for the particle id equal to the thread id is added
  // to global memory at the end.

  if (i < n) {
    // Summation for particle i
    for (int j = i + 1; j < n; j++) {
      sum += dipole_ia_energy(i, pf, pos + 3 * i, pos + 3 * j, dip + 3 * i,
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

void DipolarDirectSum_kernel_wrapper_force(dds_float k, int n, float *pos,
                                           float *dip, float *f, float *torque,
                                           dds_float box_l[3],
                                           int periodic[3]) {

  const int bs = 64;
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

  dds_float *box_l_gpu;
  int *periodic_gpu;
  cuda_safe_mem(cudaMalloc((void **)&box_l_gpu, 3 * sizeof(dds_float)));
  cuda_safe_mem(cudaMalloc((void **)&periodic_gpu, 3 * sizeof(int)));
  cuda_safe_mem(cudaMemcpy(box_l_gpu, box_l, 3 * sizeof(dds_float),
                           cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(periodic_gpu, periodic, 3 * sizeof(int),
                           cudaMemcpyHostToDevice));

  KERNELCALL(DipolarDirectSum_kernel_force, grid, block, k, n, pos, dip, f,
             torque, box_l_gpu, periodic_gpu);
  cudaFree(box_l_gpu);
  cudaFree(periodic_gpu);
}

void DipolarDirectSum_kernel_wrapper_energy(dds_float k, int n, float *pos,
                                            float *dip, dds_float box_l[3],
                                            int periodic[3], float *E) {

  const int bs = 512;
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

  dds_float *box_l_gpu;
  int *periodic_gpu;
  cuda_safe_mem(cudaMalloc((void **)&box_l_gpu, 3 * sizeof(dds_float)));
  cuda_safe_mem(cudaMalloc((void **)&periodic_gpu, 3 * sizeof(int)));
  cuda_safe_mem(
      cudaMemcpy(box_l_gpu, box_l, 3 * sizeof(float), cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(periodic_gpu, periodic, 3 * sizeof(int),
                           cudaMemcpyHostToDevice));

  dds_float *energySum;
  cuda_safe_mem(cudaMalloc(&energySum, (int)(sizeof(dds_float) * grid.x)));

  // This will sum the energies up to the block level
  KERNELCALL_shared(DipolarDirectSum_kernel_energy, grid, block,
                    bs * sizeof(dds_float), k, n, pos, dip, box_l_gpu,
                    periodic_gpu, energySum);

  // Sum the results of all blocks
  // One thread per block in the prev kernel
  // KERNELCALL(sumKernel,1,1,energySum,block.x,E);
  thrust::device_ptr<dds_float> t(energySum);
  float x = thrust::reduce(t, t + grid.x);
  cuda_safe_mem(cudaMemcpy(E, &x, sizeof(float), cudaMemcpyHostToDevice));

  cuda_safe_mem(cudaFree(energySum));
  cuda_safe_mem(cudaFree(box_l_gpu));
  cuda_safe_mem(cudaFree(periodic_gpu));
}

#endif
