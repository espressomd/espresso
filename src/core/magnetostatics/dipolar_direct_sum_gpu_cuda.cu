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

// LCOV_EXCL_START
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

__device__ void dipole_ia_force(float pf, float const dr[3], float const *dip1,
                                float const *dip2, float *f1, float *torque1,
                                float *torque2
#ifdef DIPOLE_FIELD_TRACKING
                                ,
                                float *dip_fld1, float *dip_fld2
#endif
) {
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

#ifdef DIPOLE_FIELD_TRACKING
  auto const rep1 = pe3 * pe4;
  auto const rep2 = pe2 * pe4;
  dip_fld1[0] = pf * (rep1 * dr[0] - dip2[0] * r3_inv);
  dip_fld1[1] = pf * (rep1 * dr[1] - dip2[1] * r3_inv);
  dip_fld1[2] = pf * (rep1 * dr[2] - dip2[2] * r3_inv);

  dip_fld2[0] = pf * (rep2 * dr[0] - dip1[0] * r3_inv);
  dip_fld2[1] = pf * (rep2 * dr[1] - dip1[1] * r3_inv);
  dip_fld2[2] = pf * (rep2 * dr[2] - dip1[2] * r3_inv);
#endif

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
// LCOV_EXCL_STOP

__global__ void DipolarDirectSum_kernel_force(float pf, unsigned int n,
                                              const float *pos, const float *dip
#ifdef DIPOLE_FIELD_TRACKING
                                              ,
                                              float *dip_fld
#endif
                                              ,
                                              float *f, float *torque,
                                              const float box_l[3],
                                              const int periodic[3],
                                              int n_replicas) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;

  // Load particle i
  float xi = pos[3 * i + 0], yi = pos[3 * i + 1], zi = pos[3 * i + 2];

  const float *mi = dip + 3 * i;

  // Per‐thread accumulators
  float fsum[3] = {0.0f, 0.0f, 0.0f};
  float tsum[3] = {0.0f, 0.0f, 0.0f};
#ifdef DIPOLE_FIELD_TRACKING
  float dfsum[3] = {0.0f, 0.0f, 0.0f};
#endif

  float fi[3], ti1[3], ti2[3];
#ifdef DIPOLE_FIELD_TRACKING
  float dfi[3], dfj[3];
#endif

  // --- Minimal‐image branch (no replicas) ---
  if (n_replicas == 0) {
    for (unsigned int j = i + 1; j < n; ++j) {
      // Compute MIC vector
      float d0[3];
      get_mi_vector_dds(d0, (float[]){xi, yi, zi}, pos + 3 * j, box_l,
                        periodic);
      // Compute force/torque
      dipole_ia_force(pf, d0, mi, dip + 3 * j, fi, ti1, ti2
#ifdef DIPOLE_FIELD_TRACKING
                      ,
                      dfi, dfj
#endif
      );

      // Deposit on j (atomic)
      for (int k = 0; k < 3; ++k) {
        atomicAdd(f + 3 * j + k, -fi[k]);
        atomicAdd(torque + 3 * j + k, ti2[k]);
#ifdef DIPOLE_FIELD_TRACKING
        atomicAdd(dip_fld + 3 * j + k, dfj[k]);
#endif
      }
      // Accumulate on i
      for (int k = 0; k < 3; ++k) {
        fsum[k] += fi[k];
        tsum[k] += ti1[k];
#ifdef DIPOLE_FIELD_TRACKING
        dfsum[k] += dfi[k];
#endif
      }
    }
    // Flush i’s result
    for (int k = 0; k < 3; ++k) {
      atomicAdd(f + 3 * i + k, fsum[k]);
      atomicAdd(torque + 3 * i + k, tsum[k]);
#ifdef DIPOLE_FIELD_TRACKING
      atomicAdd(dip_fld + 3 * i + k, dfsum[k]);
#endif
    }
    return;
  }
  // --- Replica branch (n_replicas > 0) ---
  int nx = periodic[0] * n_replicas;
  int ny = periodic[1] * n_replicas;
  int nz = periodic[2] * n_replicas;

  // Self‐images of i (exclude primary)
  for (int dx = -nx; dx <= nx; ++dx) {
    for (int dy = -ny; dy <= ny; ++dy) {
      for (int dz = -nz; dz <= nz; ++dz) {
        if (dx == 0 && dy == 0 && dz == 0)
          continue;
        float dr[3] = {dx * box_l[0], dy * box_l[1], dz * box_l[2]};
        dipole_ia_force(pf, dr, mi, mi, fi, ti1, ti2
#ifdef DIPOLE_FIELD_TRACKING
                        ,
                        dfi, dfj
#endif
        );
        for (int k = 0; k < 3; ++k) {
          fsum[k] += fi[k];
          tsum[k] += ti1[k];
#ifdef DIPOLE_FIELD_TRACKING
          dfsum[k] += dfi[k];
#endif
        }
      }
    }
  }

  // Pairwise i <-> j over all images
  for (unsigned int j = i + 1; j < n; ++j) {
    float xj = pos[3 * j + 0], yj = pos[3 * j + 1], zj = pos[3 * j + 2];
    const float *mj = dip + 3 * j;

    for (int dx = -nx; dx <= nx; ++dx) {
      for (int dy = -ny; dy <= ny; ++dy) {
        for (int dz = -nz; dz <= nz; ++dz) {
          float dr[3] = {(xi - xj) + dx * box_l[0], (yi - yj) + dy * box_l[1],
                         (zi - zj) + dz * box_l[2]};
          dipole_ia_force(pf, dr, mi, mj, fi, ti1, ti2
#ifdef DIPOLE_FIELD_TRACKING
                          ,
                          dfi, dfj
#endif
          );

          // Deposit on j (atomic)
          for (int k = 0; k < 3; ++k) {
            atomicAdd(f + 3 * j + k, -fi[k]);
            atomicAdd(torque + 3 * j + k, ti2[k]);
#ifdef DIPOLE_FIELD_TRACKING
            atomicAdd(dip_fld + 3 * j + k, dfj[k]);
#endif
            fsum[k] += fi[k];
            tsum[k] += ti1[k];
#ifdef DIPOLE_FIELD_TRACKING
            dfsum[k] += dfi[k];
#endif
          }
        }
      }
    }
  }

  // Add i’s total to global memory **atomically**,
  // in case any other asynchronous update might race
  for (int k = 0; k < 3; ++k) {
    atomicAdd(f + 3 * i + k, fsum[k]);
    atomicAdd(torque + 3 * i + k, tsum[k]);
#ifdef DIPOLE_FIELD_TRACKING
    atomicAdd(dip_fld + 3 * i + k, dfsum[k]);
#endif
  }
}
// LCOV_EXCL_START
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
// LCOV_EXCL_STOP

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
                                           float *dip

#ifdef DIPOLE_FIELD_TRACKING
                                           ,
                                           float *dip_fld
#endif
                                           ,
                                           float *f, float *torque,
                                           float box_l[3], int periodic[3],
                                           int n_replicas) {

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

  KERNELCALL(DipolarDirectSum_kernel_force, grid, block, k, n, pos, dip
#ifdef DIPOLE_FIELD_TRACKING
             ,
             dip_fld
#endif
             ,
             f, torque, box_l_gpu, periodic_gpu, n_replicas);
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
