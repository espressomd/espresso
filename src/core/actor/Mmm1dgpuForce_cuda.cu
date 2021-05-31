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
/** \file
 *  This file contains the code for the polygamma expansions used for the
 *  near formulas of MMM1D on GPU, as well as the force kernels.
 */

#include "config.hpp"

#ifdef MMM1D_GPU

#include "EspressoSystemInterface.hpp"
#include "actor/Mmm1dgpuForce.hpp"
#include "actor/specfunc_cuda.hpp"
#include "cuda_utils.cuh"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/mmm-modpsi.hpp"
#include "electrostatics_magnetostatics/mmm1d.hpp"

#include <utils/constants.hpp>

#include <cuda.h>

#include <cstddef>
#include <iostream>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

// the code is mostly multi-GPU capable, but ESPResSo is not yet
const int deviceCount = 1;

#undef cudaSetDevice
#define cudaSetDevice(d)

__constant__ float far_switch_radius_2[1] = {0.05f * 0.05f};
__constant__ float boxz[1];
__constant__ float uz[1];
__constant__ float coulomb_prefactor[1] = {1.0f};
__constant__ int bessel_cutoff[1] = {5};
__constant__ float maxPWerror[1] = {1e-5f};

// order hardcoded. mmm1d_recalcTables() typically does order less than 30.
// As the coefficients are stored in __constant__ memory, the array needs to be
// sized in advance. We don't know exactly how many coefficients per order, so
// we size plentiful.
const int modpsi_order = 30;
const int modpsi_constant_size = modpsi_order * modpsi_order * 2;

// linearized array on device
__constant__ int device_n_modPsi[1] = {0};
__constant__ int device_linModPsi_offsets[2 * modpsi_order],
    device_linModPsi_lengths[2 * modpsi_order];
__constant__ float device_linModPsi[modpsi_constant_size];

__device__ float dev_mod_psi_even(int n, float x) {
  return evaluateAsTaylorSeriesAt(
      &device_linModPsi[device_linModPsi_offsets[2 * n]],
      device_linModPsi_lengths[2 * n], x * x);
}

__device__ float dev_mod_psi_odd(int n, float x) {
  return x * evaluateAsTaylorSeriesAt(
                 &device_linModPsi[device_linModPsi_offsets[2 * n + 1]],
                 device_linModPsi_lengths[2 * n + 1], x * x);
}

int modpsi_init() {
  create_mod_psi_up_to(modpsi_order);

  // linearized array on host
  std::vector<int> linModPsi_offsets(modPsi.size());
  std::vector<int> linModPsi_lengths(modPsi.size());
  for (size_t i = 0; i < modPsi.size(); i++) {
    if (i)
      linModPsi_offsets[i] =
          linModPsi_offsets[i - 1] + linModPsi_lengths[i - 1];
    linModPsi_lengths[i] = modPsi[i].size();
  }

  // linearize the coefficients array
  std::vector<float> linModPsi(linModPsi_offsets[modPsi.size() - 1] +
                               linModPsi_lengths[modPsi.size() - 1]);
  for (size_t i = 0; i < modPsi.size(); i++) {
    for (size_t j = 0; j < modPsi[i].size(); j++) {
      linModPsi[linModPsi_offsets[i] + j] = static_cast<float>(modPsi[i][j]);
    }
  }

  for (int d = 0; d < deviceCount; d++) {
    cudaSetDevice(d);

    // copy to GPU
    int linModPsiSize = linModPsi_offsets[modPsi.size() - 1] +
                        linModPsi_lengths[modPsi.size() - 1];
    if (linModPsiSize > modpsi_constant_size) {
      printf("ERROR: __constant__ device_linModPsi[] is not large enough\n");
      std::abort();
    }
    cuda_safe_mem(cudaMemcpyToSymbol(device_linModPsi_offsets,
                                     linModPsi_offsets.data(),
                                     modPsi.size() * sizeof(int)));
    cuda_safe_mem(cudaMemcpyToSymbol(device_linModPsi_lengths,
                                     linModPsi_lengths.data(),
                                     modPsi.size() * sizeof(int)));
    cuda_safe_mem(cudaMemcpyToSymbol(device_linModPsi, linModPsi.data(),
                                     linModPsiSize * sizeof(float)));
    auto const n_modPsi = static_cast<int>(modPsi.size() >> 1);
    cuda_safe_mem(cudaMemcpyToSymbol(device_n_modPsi, &n_modPsi, sizeof(int)));
  }

  return 0;
}

Mmm1dgpuForce::Mmm1dgpuForce(SystemInterface &s, float _coulomb_prefactor,
                             float _maxPWerror, float _far_switch_radius,
                             int _bessel_cutoff)
    : numThreads(64), host_boxz(0), host_npart(0), need_tune(true), pairs(-1),
      dev_forcePairs(nullptr), dev_energyBlocks(nullptr),
      coulomb_prefactor(_coulomb_prefactor), maxPWerror(_maxPWerror),
      far_switch_radius(_far_switch_radius), bessel_cutoff(_bessel_cutoff) {
  // interface sanity checks
  if (!s.requestFGpu())
    throw std::runtime_error("Mmm1dgpuForce needs access to forces on GPU!");

  if (!s.requestRGpu())
    throw std::runtime_error("Mmm1dgpuForce needs access to positions on GPU!");

  if (!s.requestQGpu())
    throw std::runtime_error("Mmm1dgpuForce needs access to charges on GPU!");

  // system sanity checks
  sanity_checks();

  modpsi_init();
}

void Mmm1dgpuForce::setup(SystemInterface &s) {
  if (s.box()[2] <= 0) {
    throw std::runtime_error(
        "Error: Please set box length before initializing MMM1D!");
  }
  if (need_tune && s.npart_gpu() > 0) {
    set_params(static_cast<float>(s.box()[2]),
               static_cast<float>(coulomb.prefactor), maxPWerror,
               far_switch_radius, bessel_cutoff);
    tune(s, maxPWerror, far_switch_radius, bessel_cutoff);
  }
  if (s.box()[2] != host_boxz) {
    set_params(static_cast<float>(s.box()[2]), 0, -1, -1, -1);
  }
  if (s.npart_gpu() == host_npart) // unchanged
  {
    return;
  }

  // For all but the largest systems, it is faster to store force pairs and then
  // sum them up. Atomics are just so slow: so unless we're limited by memory,
  // do the latter.
  pairs = 2;
  for (int d = 0; d < deviceCount; d++) {
    cudaSetDevice(d);

    size_t freeMem, totalMem;
    cudaMemGetInfo(&freeMem, &totalMem);
    if (freeMem / 2 <
        3 * s.npart_gpu() * s.npart_gpu() *
            sizeof(float)) // don't use more than half the device's memory
    {
      std::cerr << "Switching to atomicAdd due to memory constraints."
                << std::endl;
      pairs = 0;
      break;
    }
  }
  if (dev_forcePairs)
    cudaFree(dev_forcePairs);
  if (pairs) // we need memory to store force pairs
  {
    cuda_safe_mem(
        cudaMalloc((void **)&dev_forcePairs,
                   3 * s.npart_gpu() * s.npart_gpu() * sizeof(float)));
  }
  if (dev_energyBlocks)
    cudaFree(dev_energyBlocks);
  cuda_safe_mem(
      cudaMalloc((void **)&dev_energyBlocks, numBlocks(s) * sizeof(float)));
  host_npart = static_cast<int>(s.npart_gpu());
}

unsigned int Mmm1dgpuForce::numBlocks(SystemInterface &s) {
  auto b = static_cast<int>(s.npart_gpu() * s.npart_gpu() / numThreads) + 1;
  if (b > 65535)
    b = 65535;
  return b;
}

Mmm1dgpuForce::~Mmm1dgpuForce() { cudaFree(dev_forcePairs); }

__forceinline__ __device__ float sqpow(float x) { return x * x; }
__forceinline__ __device__ float cbpow(float x) { return x * x * x; }

__device__ void sumReduction(float *input, float *sum) {
  auto tid = static_cast<int>(threadIdx.x);
  for (auto i = static_cast<int>(blockDim.x) / 2; i > 0; i /= 2) {
    __syncthreads();
    if (tid < i)
      input[tid] += input[i + tid];
  }
  __syncthreads();
  if (tid == 0)
    sum[0] = input[0];
}

__global__ void sumKernel(float *data, int N) {
  extern __shared__ float partialsums[];
  if (blockIdx.x != 0)
    return;
  auto tid = static_cast<int>(threadIdx.x);
  float result = 0;

  for (int i = 0; i < N; i += static_cast<int>(blockDim.x)) {
    if (i + tid >= N)
      partialsums[tid] = 0;
    else
      partialsums[tid] = data[i + tid];

    sumReduction(partialsums, &result);
    if (tid == 0) {
      if (i == 0)
        data[0] = 0;
      data[0] += result;
    }
  }
}

__global__ void besselTuneKernel(int *result, float far_switch_radius,
                                 int maxCut) {
  const float c_2pif = 2 * Utils::pi<float>();
  float arg = c_2pif * *uz * far_switch_radius;
  float pref = 4 * *uz * max(1.0f, c_2pif * *uz);
  float err;
  int P = 1;
  do {
    err = pref * dev_K1(arg * static_cast<float>(P)) * exp(arg) / arg *
          (static_cast<float>(P) - 1 + 1 / arg);
    P++;
  } while (err > *maxPWerror && P <= maxCut);
  P--;

  result[0] = P;
}

void Mmm1dgpuForce::tune(SystemInterface &s, float _maxPWerror,
                         float _far_switch_radius, int _bessel_cutoff) {
  float far_switch_radius = _far_switch_radius;
  int bessel_cutoff = _bessel_cutoff;
  float maxrad = host_boxz;

  if (_far_switch_radius < 0 && _bessel_cutoff < 0)
  // autodetermine switching radius and Bessel cutoff
  {
    float bestrad = 0, besttime = INFINITY;

    // NOLINTNEXTLINE(clang-analyzer-security.FloatLoopCounter)
    for (far_switch_radius = 0.05f * maxrad; far_switch_radius < maxrad;
         far_switch_radius += 0.05f * maxrad) {
      set_params(0, 0, _maxPWerror, far_switch_radius, bessel_cutoff);
      tune(s, _maxPWerror, far_switch_radius, -2); // tune Bessel cutoff
      auto runtime = force_benchmark(s);
      if (runtime < besttime) {
        besttime = runtime;
        bestrad = far_switch_radius;
      }
    }
    far_switch_radius = bestrad;

    set_params(0, 0, _maxPWerror, far_switch_radius, bessel_cutoff);
    tune(s, _maxPWerror, far_switch_radius, -2); // tune Bessel cutoff
  }

  else if (_bessel_cutoff < 0)
  // autodetermine Bessel cutoff
  {
    int *dev_cutoff;
    int maxCut = 30;
    cuda_safe_mem(cudaMalloc((void **)&dev_cutoff, sizeof(int)));
    besselTuneKernel<<<dim3(1), dim3(1), 0, nullptr>>>(
        dev_cutoff, far_switch_radius, maxCut);
    cuda_safe_mem(cudaMemcpy(&bessel_cutoff, dev_cutoff, sizeof(int),
                             cudaMemcpyDeviceToHost));
    cudaFree(dev_cutoff);
    if (_bessel_cutoff != -2 &&
        bessel_cutoff >=
            maxCut) // we already have our switching radius and only need to
                    // determine the cutoff, i.e. this is the final tuning round
    {
      throw std::runtime_error(
          "No reasonable Bessel cutoff could be determined.");
    }

    set_params(0, 0, _maxPWerror, far_switch_radius, bessel_cutoff);
  }
}

void Mmm1dgpuForce::set_params(float _boxz, float _coulomb_prefactor,
                               float _maxPWerror, float _far_switch_radius,
                               int _bessel_cutoff, bool manual) {
  if (_boxz > 0 && _far_switch_radius > _boxz) {
    throw std::runtime_error(
        "switching radius must not be larger than box length");
  }
  float _far_switch_radius_2 = _far_switch_radius * _far_switch_radius;
  float _uz = 1.0f / _boxz;
  for (int d = 0; d < deviceCount; d++) {
    // double colons are needed to access the constant memory variables because
    // they are file globals and we have identically named class variables
    cudaSetDevice(d);
    if (manual) // tuning needs to be performed again
    {
      far_switch_radius = _far_switch_radius;
      bessel_cutoff = _bessel_cutoff;
    }
    if (_far_switch_radius >= 0) {
      mmm1d_params.far_switch_radius_2 =
          _far_switch_radius * _far_switch_radius;
      cuda_safe_mem(cudaMemcpyToSymbol(::far_switch_radius_2,
                                       &_far_switch_radius_2, sizeof(float)));
      far_switch_radius = _far_switch_radius;
    }
    if (_boxz > 0) {
      host_boxz = _boxz;
      cuda_safe_mem(cudaMemcpyToSymbol(::boxz, &_boxz, sizeof(float)));
      cuda_safe_mem(cudaMemcpyToSymbol(::uz, &_uz, sizeof(float)));
    }
    if (_coulomb_prefactor != 0) {
      cuda_safe_mem(cudaMemcpyToSymbol(::coulomb_prefactor, &_coulomb_prefactor,
                                       sizeof(float)));
      coulomb_prefactor = _coulomb_prefactor;
    }
    if (_bessel_cutoff > 0) {
      mmm1d_params.bessel_cutoff = _bessel_cutoff;
      cuda_safe_mem(
          cudaMemcpyToSymbol(::bessel_cutoff, &_bessel_cutoff, sizeof(int)));
      bessel_cutoff = _bessel_cutoff;
    }
    if (_maxPWerror > 0) {
      mmm1d_params.maxPWerror = _maxPWerror;
      cuda_safe_mem(
          cudaMemcpyToSymbol(::maxPWerror, &_maxPWerror, sizeof(float)));
      maxPWerror = _maxPWerror;
    }
  }
  need_tune = true;

  // The changed parameters in mmm1d_params do not need to be broadcast: they
  // are only accessed by the TCL print function (on node 0) when you call inter
  // coulomb. The CUDA code only runs on node 0, so other nodes do not need the
  // parameters. We couldn't broadcast from here anyway because set_params()
  // might be called from inside computeForces() which is not a time at which
  // the MPI loop on the slave nodes is waiting for broadcasts.
}

__global__ void forcesKernel(const float *__restrict__ r,
                             const float *__restrict__ q,
                             float *__restrict__ force, int N, int pairs,
                             int tStart, int tStop) {
  if (tStop < 0)
    tStop = N * N;

  const float c_2pif = 2 * Utils::pi<float>();

  for (int tid =
           static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x) + tStart;
       tid < tStop; tid += static_cast<int>(blockDim.x * gridDim.x)) {
    int p1 = tid % N, p2 = tid / N;
    float x = r[3 * p2] - r[3 * p1], y = r[3 * p2 + 1] - r[3 * p1 + 1],
          z = r[3 * p2 + 2] - r[3 * p1 + 2];
    float rxy2 = sqpow(x) + sqpow(y);
    float rxy = sqrt(rxy2);
    float sum_r = 0, sum_z = 0;

    // if (*boxz <= 0.0) return; // in case we are not initialized yet

    while (fabs(z) > *boxz / 2) // make sure we take the shortest distance
      z -= (z > 0 ? 1.f : -1.f) * *boxz;

    if (p1 == p2) // particle exerts no force on itself
    {
      rxy = 1; // so the division at the end doesn't fail with NaN (sum_r is 0
               // anyway)
    } else if (rxy2 <= *far_switch_radius_2) // near formula
    {
      float uzz = *uz * z;
      float uzr = *uz * rxy;
      sum_z = dev_mod_psi_odd(0, uzz);
      float uzrpow = uzr;
      for (int n = 1; n < *device_n_modPsi; n++) {
        float sum_r_old = sum_r;
        float mpe = dev_mod_psi_even(n, uzz);
        float mpo = dev_mod_psi_odd(n, uzz);

        sum_r += 2 * static_cast<float>(n) * mpe * uzrpow;
        uzrpow *= uzr;
        sum_z += mpo * uzrpow;
        uzrpow *= uzr;

        if (fabs(sum_r_old - sum_r) < *maxPWerror)
          break;
      }

      sum_r *= sqpow(*uz);
      sum_z *= sqpow(*uz);

      sum_r += rxy * cbpow(rsqrt(rxy2 + sqpow(z)));
      sum_r += rxy * cbpow(rsqrt(rxy2 + sqpow(z + *boxz)));
      sum_r += rxy * cbpow(rsqrt(rxy2 + sqpow(z - *boxz)));

      sum_z += z * cbpow(rsqrt(rxy2 + sqpow(z)));
      sum_z += (z + *boxz) * cbpow(rsqrt(rxy2 + sqpow(z + *boxz)));
      sum_z += (z - *boxz) * cbpow(rsqrt(rxy2 + sqpow(z - *boxz)));

      if (rxy == 0) // particles at the same radial position only exert a force
                    // in z direction
      {
        rxy = 1; // so the division at the end doesn't fail with NaN (sum_r is 0
                 // anyway)
      }
    } else // far formula
    {
      for (int p = 1; p < *bessel_cutoff; p++) {
        float arg = c_2pif * *uz * static_cast<float>(p);
        sum_r += static_cast<float>(p) * dev_K1(arg * rxy) * cos(arg * z);
        sum_z += static_cast<float>(p) * dev_K0(arg * rxy) * sin(arg * z);
      }
      sum_r *= sqpow(*uz) * 4 * c_2pif;
      sum_z *= sqpow(*uz) * 4 * c_2pif;
      sum_r += 2 * *uz / rxy;
    }

    float pref = *coulomb_prefactor * q[p1] * q[p2];
    if (pairs) {
      force[3 * (p1 + p2 * N - tStart)] = pref * sum_r / rxy * x;
      force[3 * (p1 + p2 * N - tStart) + 1] = pref * sum_r / rxy * y;
      force[3 * (p1 + p2 * N - tStart) + 2] = pref * sum_z;
    } else {
      atomicAdd(&force[3 * p2], pref * sum_r / rxy * x);
      atomicAdd(&force[3 * p2 + 1], pref * sum_r / rxy * y);
      atomicAdd(&force[3 * p2 + 2], pref * sum_z);
    }
  }
}

__global__ void energiesKernel(const float *__restrict__ r,
                               const float *__restrict__ q,
                               float *__restrict__ energy, int N, int pairs,
                               int tStart, int tStop) {
  if (tStop < 0)
    tStop = N * N;

  auto const c_2pif = 2 * Utils::pi<float>();
  auto const c_gammaf = Utils::gamma<float>();

  extern __shared__ float partialsums[];
  if (!pairs) {
    partialsums[threadIdx.x] = 0;
    __syncthreads();
  }
  for (int tid =
           static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x) + tStart;
       tid < tStop; tid += static_cast<int>(blockDim.x * gridDim.x)) {
    int p1 = tid % N, p2 = tid / N;
    float z = r[3 * p2 + 2] - r[3 * p1 + 2];
    float rxy2 =
        sqpow(r[3 * p2] - r[3 * p1]) + sqpow(r[3 * p2 + 1] - r[3 * p1 + 1]);
    float rxy = sqrt(rxy2);
    float sum_e = 0;

    // if (*boxz <= 0.0) return; // in case we are not initialized yet

    while (fabs(z) > *boxz / 2) // make sure we take the shortest distance
      z -= (z > 0 ? 1.f : -1.f) * *boxz;

    if (p1 == p2) // particle exerts no force on itself
    {
    } else if (rxy2 <= *far_switch_radius_2) // near formula
    {
      float uzz = *uz * z;
      float uzr2 = sqpow(*uz * rxy);
      float uzrpow = uzr2;
      sum_e = dev_mod_psi_even(0, uzz);
      for (int n = 1; n < *device_n_modPsi; n++) {
        float sum_e_old = sum_e;
        float mpe = dev_mod_psi_even(n, uzz);
        sum_e += mpe * uzrpow;
        uzrpow *= uzr2;

        if (fabs(sum_e_old - sum_e) < *maxPWerror)
          break;
      }

      sum_e *= -1 * *uz;
      sum_e -= 2 * *uz * c_gammaf;
      sum_e += rsqrt(rxy2 + sqpow(z));
      sum_e += rsqrt(rxy2 + sqpow(z + *boxz));
      sum_e += rsqrt(rxy2 + sqpow(z - *boxz));
    } else // far formula
    {
      sum_e = -(log(rxy * *uz / 2) + c_gammaf) / 2;
      for (int p = 1; p < *bessel_cutoff; p++) {
        float arg = c_2pif * *uz * static_cast<float>(p);
        sum_e += dev_K0(arg * rxy) * cos(arg * z);
      }
      sum_e *= *uz * 4;
    }

    if (pairs) {
      energy[p1 + p2 * N - tStart] = *coulomb_prefactor * q[p1] * q[p2] * sum_e;
    } else {
      partialsums[threadIdx.x] += *coulomb_prefactor * q[p1] * q[p2] * sum_e;
    }
  }
  if (!pairs) {
    sumReduction(partialsums, &energy[blockIdx.x]);
  }
}

__global__ void vectorReductionKernel(float const *src, float *dst, int N,
                                      int tStart, int tStop) {
  if (tStop < 0)
    tStop = N * N;

  for (auto tid = static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x);
       tid < N; tid += static_cast<int>(blockDim.x * gridDim.x)) {
    int offset = ((tid + (tStart % N)) % N);

    for (int i = 0; tid + i * N < (tStop - tStart); i++) {
#pragma unroll 3
      for (int d = 0; d < 3; d++) {
        dst[3 * offset + d] -= src[3 * (tid + i * N) + d];
      }
    }
  }
}

void Mmm1dgpuForce::computeForces(SystemInterface &s) {
  if (coulomb.method != COULOMB_MMM1D_GPU) {
    std::cerr << "MMM1D: coulomb.method has been changed, skipping calculation"
              << std::endl;
    return;
  }
  setup(s);

  if (pairs < 0) {
    throw std::runtime_error("MMM1D was not initialized correctly");
  }

  if (pairs) // if we calculate force pairs, we need to reduce them to forces
  {
    auto blocksRed = static_cast<int>(s.npart_gpu() / numThreads) + 1;
    KERNELCALL(forcesKernel, numBlocks(s), numThreads, s.rGpuBegin(),
               s.qGpuBegin(), dev_forcePairs, s.npart_gpu(), pairs, 0, -1)
    KERNELCALL(vectorReductionKernel, blocksRed, numThreads, dev_forcePairs,
               s.fGpuBegin(), s.npart_gpu(), 0, -1)
  } else {
    KERNELCALL(forcesKernel, numBlocks(s), numThreads, s.rGpuBegin(),
               s.qGpuBegin(), s.fGpuBegin(), s.npart_gpu(), pairs, 0, -1)
  }
}

__global__ void scaleAndAddKernel(float *dst, float const *src, int N,
                                  float factor) {
  for (auto tid = static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x);
       tid < N; tid += static_cast<int>(blockDim.x * gridDim.x)) {
    dst[tid] += src[tid] * factor;
  }
}

void Mmm1dgpuForce::computeEnergy(SystemInterface &s) {
  if (coulomb.method != COULOMB_MMM1D_GPU) {
    std::cerr << "MMM1D: coulomb.method has been changed, skipping calculation"
              << std::endl;
    return;
  }
  setup(s);

  if (pairs < 0) {
    throw std::runtime_error("MMM1D was not initialized correctly");
  }
  auto shared = static_cast<int>(numThreads * sizeof(float));

  KERNELCALL_shared(energiesKernel, numBlocks(s), numThreads, shared,
                    s.rGpuBegin(), s.qGpuBegin(), dev_energyBlocks,
                    s.npart_gpu(), 0, 0, -1);
  KERNELCALL_shared(sumKernel, 1, numThreads, shared, dev_energyBlocks,
                    numBlocks(s));
  KERNELCALL(scaleAndAddKernel, 1, 1, &(((CUDA_energy *)s.eGpu())->coulomb),
             &dev_energyBlocks[0], 1,
             0.5); // we have counted every interaction twice, so halve the
                   // total energy
}

float Mmm1dgpuForce::force_benchmark(SystemInterface &s) {
  cudaEvent_t eventStart, eventStop;
  float elapsedTime;
  float *dev_f_benchmark;

  cuda_safe_mem(
      cudaMalloc((void **)&dev_f_benchmark, 3 * s.npart_gpu() * sizeof(float)));
  cuda_safe_mem(cudaEventCreate(&eventStart));
  cuda_safe_mem(cudaEventCreate(&eventStop));
  cuda_safe_mem(cudaEventRecord(eventStart, stream[0]));
  KERNELCALL(forcesKernel, numBlocks(s), numThreads, s.rGpuBegin(),
             s.qGpuBegin(), dev_f_benchmark, s.npart_gpu(), 0, 0, -1)
  cuda_safe_mem(cudaEventRecord(eventStop, stream[0]));
  cuda_safe_mem(cudaEventSynchronize(eventStop));
  cuda_safe_mem(cudaEventElapsedTime(&elapsedTime, eventStart, eventStop));
  cuda_safe_mem(cudaEventDestroy(eventStart));
  cuda_safe_mem(cudaEventDestroy(eventStop));
  cuda_safe_mem(cudaFree(dev_f_benchmark));

  return elapsedTime;
}

#endif /* MMM1D_GPU */
