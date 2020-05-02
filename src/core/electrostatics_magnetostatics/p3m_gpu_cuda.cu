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

/** \file
 *
 *  P3M electrostatics on GPU.
 *
 *  The corresponding header file is p3m_gpu.hpp.
 */

#include "config.hpp"

#ifdef ELECTROSTATICS

#define _P3M_GPU_FLOAT
//#define _P3M_GPU_REAL_DOUBLE

#ifdef _P3M_GPU_FLOAT
#define REAL_TYPE float
#define FFT_TYPE_COMPLEX cufftComplex
#define FFT_FORW_FFT cufftExecR2C
#define FFT_BACK_FFT cufftExecC2R
#define FFT_PLAN_FORW_FLAG CUFFT_R2C
#define FFT_PLAN_BACK_FLAG CUFFT_C2R
#endif

#ifdef _P3M_GPU_REAL_DOUBLE
#define REAL_TYPE double
#define FFT_TYPE_COMPLEX cufftDoubleComplex
#define FFT_FORW_FFT cufftExecD2Z
#define FFT_BACK_FFT cufftExecZ2D
#define FFT_PLAN_FORW_FLAG CUFFT_D2Z
#define FFT_PLAN_BACK_FLAG CUFFT_Z2D
#endif

#include "BoxGeometry.hpp"

#include <cstdio>
#include <cstdlib>

#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "cufft_wrapper.hpp"

#include "electrostatics_magnetostatics/p3m_gpu.hpp"

#include "EspressoSystemInterface.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "global.hpp"

#include <utils/math/int_pow.hpp>
using Utils::int_pow;
#include <utils/math/sinc.hpp>
#include <utils/math/sqr.hpp>
using Utils::sqr;
#include <utils/math/bspline.hpp>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

extern BoxGeometry box_geo;

struct P3MGpuData {
  /** Charge mesh */
  FFT_TYPE_COMPLEX *charge_mesh;
  /** Force meshes */
  FFT_TYPE_COMPLEX *force_mesh_x;
  FFT_TYPE_COMPLEX *force_mesh_y;
  FFT_TYPE_COMPLEX *force_mesh_z;
  /** Influence Function */
  REAL_TYPE *G_hat;
  /** Charge assignment order */
  int cao;
  /** Total number of mesh points (including padding) */
  int mesh_size;
  /** Ewald parameter */
  REAL_TYPE alpha;
  /** Number of particles */
  int n_part;
  /** Box size */
  REAL_TYPE box[3];
  /** Mesh dimensions */
  int mesh[3];
  /** Padded size */
  int mesh_z_padded;
  /** Inverse mesh spacing */
  REAL_TYPE hi[3];
  /** Position shift */
  REAL_TYPE pos_shift;
};

P3MGpuData p3m_gpu_data;

struct p3m_gpu_fft_plans_t {
  /** FFT plans */
  cufftHandle forw_plan;
  cufftHandle back_plan;
} p3m_gpu_fft_plans;

static char p3m_gpu_data_initialized = 0;

template <int cao>
__device__ void static Aliasing_sums_ik(const P3MGpuData p, int NX, int NY,
                                        int NZ, REAL_TYPE *Zaehler,
                                        REAL_TYPE *Nenner) {
  REAL_TYPE S1, S2, S3;
  REAL_TYPE zwi;
  int MX, MY, MZ;
  REAL_TYPE NMX, NMY, NMZ;
  REAL_TYPE NM2;
  REAL_TYPE TE;
  REAL_TYPE Leni[3];
  REAL_TYPE Meshi[3];
  for (int i = 0; i < 3; ++i) {
    Leni[i] = 1.0f / p.box[i];
    Meshi[i] = 1.0f / static_cast<REAL_TYPE>(p.mesh[i]);
  }

  Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner = 0.0;

  for (MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++) {
    NMX = static_cast<REAL_TYPE>(((NX > p.mesh[0] / 2) ? NX - p.mesh[0] : NX) +
                                 p.mesh[0] * MX);
    S1 = int_pow<2 * cao>(Utils::sinc(Meshi[0] * NMX));
    for (MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++) {
      NMY = static_cast<REAL_TYPE>(
          ((NY > p.mesh[1] / 2) ? NY - p.mesh[1] : NY) + p.mesh[1] * MY);
      S2 = S1 * int_pow<2 * cao>(Utils::sinc(Meshi[1] * NMY));
      for (MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++) {
        NMZ = static_cast<REAL_TYPE>(
            ((NZ > p.mesh[2] / 2) ? NZ - p.mesh[2] : NZ) + p.mesh[2] * MZ);
        S3 = S2 * int_pow<2 * cao>(Utils::sinc(Meshi[2] * NMZ));

        NM2 = sqr(NMX * Leni[0]) + sqr(NMY * Leni[1]) + sqr(NMZ * Leni[2]);
        *Nenner += S3;

        TE = exp(-sqr(Utils::pi<float>() / (p.alpha)) * NM2);
        zwi = S3 * TE / NM2;
        Zaehler[0] += NMX * zwi * Leni[0];
        Zaehler[1] += NMY * zwi * Leni[1];
        Zaehler[2] += NMZ * zwi * Leni[2];
      }
    }
  }
}

/* Calculate influence function */
template <int cao>
__global__ void calculate_influence_function_device(const P3MGpuData p) {

  const auto NX = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
  const auto NY = static_cast<int>(blockDim.y * blockIdx.y + threadIdx.y);
  const auto NZ = static_cast<int>(blockDim.z * blockIdx.z + threadIdx.z);
  REAL_TYPE Dnx, Dny, Dnz;
  REAL_TYPE Zaehler[3] = {0.0, 0.0, 0.0}, Nenner = 0.0;
  REAL_TYPE zwi;
  int ind = 0;
  REAL_TYPE Leni[3];
  for (int i = 0; i < 3; ++i)
    Leni[i] = 1.0f / p.box[i];

  if ((NX >= p.mesh[0]) || (NY >= p.mesh[1]) || (NZ >= (p.mesh[2] / 2 + 1)))
    return;

  ind = NX * p.mesh[1] * (p.mesh[2] / 2 + 1) + NY * (p.mesh[2] / 2 + 1) + NZ;

  if (((NX == 0) && (NY == 0) && (NZ == 0)) ||
      ((NX % (p.mesh[0] / 2) == 0) && (NY % (p.mesh[1] / 2) == 0) &&
       (NZ % (p.mesh[2] / 2) == 0)))
    p.G_hat[ind] = 0.0;
  else {
    Aliasing_sums_ik<cao>(p, NX, NY, NZ, Zaehler, &Nenner);

    Dnx = static_cast<REAL_TYPE>((NX > p.mesh[0] / 2) ? NX - p.mesh[0] : NX);
    Dny = static_cast<REAL_TYPE>((NY > p.mesh[1] / 2) ? NY - p.mesh[1] : NY);
    Dnz = static_cast<REAL_TYPE>((NZ > p.mesh[2] / 2) ? NZ - p.mesh[2] : NZ);

    zwi = Dnx * Zaehler[0] * Leni[0] + Dny * Zaehler[1] * Leni[1] +
          Dnz * Zaehler[2] * Leni[2];
    zwi /= ((sqr(Dnx * Leni[0]) + sqr(Dny * Leni[1]) + sqr(Dnz * Leni[2])) *
            sqr(Nenner));
    p.G_hat[ind] = 2.0f * zwi / Utils::pi<float>();
  }
}

#ifdef _P3M_GPU_REAL_DOUBLE
__device__ double atomicAdd(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

namespace {
__device__ inline int linear_index_r(P3MGpuData const &p, int i, int j, int k) {
  return p.mesh[1] * p.mesh_z_padded * i + p.mesh_z_padded * j + k;
}

__device__ inline int linear_index_k(P3MGpuData const &p, int i, int j, int k) {
  return p.mesh[1] * (p.mesh[2] / 2 + 1) * i + (p.mesh[2] / 2 + 1) * j + k;
}
} // namespace

__global__ void apply_diff_op(const P3MGpuData p) {
  const int linear_index = linear_index_k(p, static_cast<int>(blockIdx.x),
                                          static_cast<int>(blockIdx.y),
                                          static_cast<int>(threadIdx.x));

  auto const nx = static_cast<int>(
      (blockIdx.x > p.mesh[0] / 2) ? blockIdx.x - p.mesh[0] : blockIdx.x);
  auto const ny = static_cast<int>(
      (blockIdx.y > p.mesh[1] / 2) ? blockIdx.y - p.mesh[1] : blockIdx.y);
  auto const nz = static_cast<int>(threadIdx.x);

  const FFT_TYPE_COMPLEX meshw = p.charge_mesh[linear_index];
  FFT_TYPE_COMPLEX buf;
  buf.x = -2.0f * Utils::pi<float>() * meshw.y;
  buf.y = 2.0f * Utils::pi<float>() * meshw.x;

  p.force_mesh_x[linear_index].x =
      static_cast<decltype(FFT_TYPE_COMPLEX::x)>(nx) * buf.x / p.box[0];
  p.force_mesh_x[linear_index].y =
      static_cast<decltype(FFT_TYPE_COMPLEX::x)>(nx) * buf.y / p.box[0];

  p.force_mesh_y[linear_index].x =
      static_cast<decltype(FFT_TYPE_COMPLEX::x)>(ny) * buf.x / p.box[1];
  p.force_mesh_y[linear_index].y =
      static_cast<decltype(FFT_TYPE_COMPLEX::x)>(ny) * buf.y / p.box[1];

  p.force_mesh_z[linear_index].x =
      static_cast<decltype(FFT_TYPE_COMPLEX::x)>(nz) * buf.x / p.box[2];
  p.force_mesh_z[linear_index].y =
      static_cast<decltype(FFT_TYPE_COMPLEX::x)>(nz) * buf.y / p.box[2];
}

__device__ inline int wrap_index(const int ind, const int mesh) {
  if (ind < 0)
    return ind + mesh;
  if (ind >= mesh)
    return ind - mesh;
  return ind;
}

__global__ void apply_influence_function(const P3MGpuData p) {
  const int linear_index = linear_index_k(p, static_cast<int>(blockIdx.x),
                                          static_cast<int>(blockIdx.y),
                                          static_cast<int>(threadIdx.x));

  p.charge_mesh[linear_index].x *= p.G_hat[linear_index];
  p.charge_mesh[linear_index].y *= p.G_hat[linear_index];
}

template <int cao, bool shared>
__global__ void assign_charge_kernel(const CUDA_particle_data *const pdata,
                                     const P3MGpuData par,
                                     const int parts_per_block) {
  const int part_in_block = threadIdx.x / cao;
  const int cao_id_x = threadIdx.x - part_in_block * cao;
  /* id of the particle */
  int id =
      parts_per_block * (blockIdx.x * gridDim.y + blockIdx.y) + part_in_block;
  if (id >= par.n_part)
    return;
  /* position relative to the closest gird point */
  REAL_TYPE m_pos[3];
  /* index of the nearest mesh point */
  int nmp_x, nmp_y, nmp_z;

  const CUDA_particle_data p = pdata[id];
  auto *charge_mesh = (REAL_TYPE *)par.charge_mesh;

  m_pos[0] = p.p[0] * par.hi[0] - par.pos_shift;
  m_pos[1] = p.p[1] * par.hi[1] - par.pos_shift;
  m_pos[2] = p.p[2] * par.hi[2] - par.pos_shift;

  nmp_x = static_cast<int>(floorf(m_pos[0] + 0.5f));
  nmp_y = static_cast<int>(floorf(m_pos[1] + 0.5f));
  nmp_z = static_cast<int>(floorf(m_pos[2] + 0.5f));

  m_pos[0] -= static_cast<REAL_TYPE>(nmp_x);
  m_pos[1] -= static_cast<REAL_TYPE>(nmp_y);
  m_pos[2] -= static_cast<REAL_TYPE>(nmp_z);

  nmp_x = wrap_index(nmp_x + cao_id_x, par.mesh[0]);
  nmp_y = wrap_index(nmp_y + static_cast<int>(threadIdx.y), par.mesh[1]);
  nmp_z = wrap_index(nmp_z + static_cast<int>(threadIdx.z), par.mesh[2]);

  const int ind = linear_index_r(par, nmp_x, nmp_y, nmp_z);

  HIP_DYNAMIC_SHARED(float, weights);

  if (shared) {
    if ((threadIdx.y < 3) && (threadIdx.z == 0)) {
      weights[3 * cao * part_in_block + 3 * cao_id_x + threadIdx.y] =
          Utils::bspline<cao>(cao_id_x, m_pos[threadIdx.y]);
    }

    __syncthreads();

    atomicAdd(&(charge_mesh[ind]),
              weights[3 * cao * part_in_block + 3 * cao_id_x + 0] *
                  weights[3 * cao * part_in_block + 3 * threadIdx.y + 1] *
                  weights[3 * cao * part_in_block + 3 * threadIdx.z + 2] * p.q);

  } else {
    atomicAdd(&(charge_mesh[ind]),
              Utils::bspline<cao>(cao_id_x, m_pos[0]) *
                  Utils::bspline<cao>(threadIdx.y, m_pos[1]) *
                  Utils::bspline<cao>(threadIdx.z, m_pos[2]) * p.q);
  }
}

void assign_charges(const CUDA_particle_data *const pdata, const P3MGpuData p) {
  dim3 grid, block;
  grid.z = 1;
  const int cao3 = p.cao * p.cao * p.cao;
  const int cao = p.cao;
  int parts_per_block = 1, n_blocks = 1;

  while ((parts_per_block + 1) * cao3 <= 1024) {
    parts_per_block++;
  }
  if ((p.n_part % parts_per_block) == 0)
    n_blocks = std::max<int>(1, p.n_part / parts_per_block);
  else
    n_blocks = p.n_part / parts_per_block + 1;

  grid.x = n_blocks;
  grid.y = 1;
  while (grid.x > 65536) {
    grid.y++;
    if ((n_blocks % grid.y) == 0)
      grid.x = std::max<int>(1, n_blocks / grid.y);
    else
      grid.x = n_blocks / grid.y + 1;
  }

  block.x = parts_per_block * cao;
  block.y = cao;
  block.z = cao;

  switch (cao) {
  case 1:
    hipLaunchKernelGGL((assign_charge_kernel<1, false>), dim3(grid),
                       dim3(block), 0, nullptr, pdata, p, parts_per_block);
    break;
  case 2:
    hipLaunchKernelGGL((assign_charge_kernel<2, false>), dim3(grid),
                       dim3(block), 0, nullptr, pdata, p, parts_per_block);
    break;
  case 3:
    hipLaunchKernelGGL((assign_charge_kernel<3, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), nullptr,
                       pdata, p, parts_per_block);
    break;
  case 4:
    hipLaunchKernelGGL((assign_charge_kernel<4, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), nullptr,
                       pdata, p, parts_per_block);
    break;
  case 5:
    hipLaunchKernelGGL((assign_charge_kernel<5, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), nullptr,
                       pdata, p, parts_per_block);
    break;
  case 6:
    hipLaunchKernelGGL((assign_charge_kernel<6, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), nullptr,
                       pdata, p, parts_per_block);
    break;
  case 7:
    hipLaunchKernelGGL((assign_charge_kernel<7, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), nullptr,
                       pdata, p, parts_per_block);
    break;
  default:
    break;
  }
  _cuda_check_errors(block, grid, "assign_charge", __FILE__, __LINE__);
}

template <int cao, bool shared>
__global__ void assign_forces_kernel(const CUDA_particle_data *const pdata,
                                     const P3MGpuData par,
                                     float *lb_particle_force_gpu,
                                     REAL_TYPE prefactor, int parts_per_block) {
  const int part_in_block = threadIdx.x / cao;
  const int cao_id_x = threadIdx.x - part_in_block * cao;
  /* id of the particle */
  int id =
      parts_per_block * (blockIdx.x * gridDim.y + blockIdx.y) + part_in_block;
  if (id >= par.n_part)
    return;
  /* position relative to the closest grid point */
  REAL_TYPE m_pos[3];
  /* index of the nearest mesh point */
  int nmp_x, nmp_y, nmp_z;

  const CUDA_particle_data p = pdata[id];

  m_pos[0] = p.p[0] * par.hi[0] - par.pos_shift;
  m_pos[1] = p.p[1] * par.hi[1] - par.pos_shift;
  m_pos[2] = p.p[2] * par.hi[2] - par.pos_shift;

  nmp_x = static_cast<int>(floorf(m_pos[0] + 0.5f));
  nmp_y = static_cast<int>(floorf(m_pos[1] + 0.5f));
  nmp_z = static_cast<int>(floorf(m_pos[2] + 0.5f));

  m_pos[0] -= static_cast<REAL_TYPE>(nmp_x);
  m_pos[1] -= static_cast<REAL_TYPE>(nmp_y);
  m_pos[2] -= static_cast<REAL_TYPE>(nmp_z);

  nmp_x = wrap_index(nmp_x + cao_id_x, par.mesh[0]);
  nmp_y = wrap_index(nmp_y + static_cast<int>(threadIdx.y), par.mesh[1]);
  nmp_z = wrap_index(nmp_z + static_cast<int>(threadIdx.z), par.mesh[2]);

  REAL_TYPE c;
  const int index = linear_index_r(par, nmp_x, nmp_y, nmp_z);

  HIP_DYNAMIC_SHARED(float, weights);

  if (shared) {
    if ((threadIdx.y < 3) && (threadIdx.z == 0)) {
      weights[3 * cao * part_in_block + 3 * cao_id_x + threadIdx.y] =
          Utils::bspline<cao>(cao_id_x, m_pos[threadIdx.y]);
    }

    __syncthreads();

    c = -prefactor * weights[3 * cao * part_in_block + 3 * cao_id_x + 0] *
        weights[3 * cao * part_in_block + 3 * threadIdx.y + 1] *
        weights[3 * cao * part_in_block + 3 * threadIdx.z + 2] * p.q;
  } else {
    c = -prefactor * Utils::bspline<cao>(cao_id_x, m_pos[0]) *
        Utils::bspline<cao>(threadIdx.y, m_pos[1]) *
        Utils::bspline<cao>(threadIdx.z, m_pos[2]) * p.q;
  }

  const REAL_TYPE *force_mesh_x = (REAL_TYPE *)par.force_mesh_x;
  const REAL_TYPE *force_mesh_y = (REAL_TYPE *)par.force_mesh_y;
  const REAL_TYPE *force_mesh_z = (REAL_TYPE *)par.force_mesh_z;

  atomicAdd(&(lb_particle_force_gpu[3 * id + 0]), c * force_mesh_x[index]);
  atomicAdd(&(lb_particle_force_gpu[3 * id + 1]), c * force_mesh_y[index]);
  atomicAdd(&(lb_particle_force_gpu[3 * id + 2]), c * force_mesh_z[index]);
}

void assign_forces(const CUDA_particle_data *const pdata, const P3MGpuData p,
                   float *lb_particle_force_gpu, REAL_TYPE prefactor) {
  dim3 grid, block;
  grid.z = 1;

  const int cao = p.cao;
  const int cao3 = cao * cao * cao;
  int parts_per_block = 1, n_blocks = 1;

  while ((parts_per_block + 1) * cao3 <= 1024) {
    parts_per_block++;
  }

  if ((p3m_gpu_data.n_part % parts_per_block) == 0)
    n_blocks = std::max<int>(1, p3m_gpu_data.n_part / parts_per_block);
  else
    n_blocks = p3m_gpu_data.n_part / parts_per_block + 1;

  grid.x = n_blocks;
  grid.y = 1;
  while (grid.x > 65536) {
    grid.y++;
    if ((n_blocks % grid.y) == 0)
      grid.x = std::max<int>(1, n_blocks / grid.y);
    else
      grid.x = n_blocks / grid.y + 1;
  }

  block.x = parts_per_block * cao;
  block.y = cao;
  block.z = cao;

  /* Switch for assignment templates, the shared version only is faster for cao
   * > 2 */
  switch (cao) {
  case 1:
    hipLaunchKernelGGL((assign_forces_kernel<1, false>), dim3(grid),
                       dim3(block), 0, nullptr, pdata, p, lb_particle_force_gpu,
                       prefactor, parts_per_block);
    break;
  case 2:
    hipLaunchKernelGGL((assign_forces_kernel<2, false>), dim3(grid),
                       dim3(block), 0, nullptr, pdata, p, lb_particle_force_gpu,
                       prefactor, parts_per_block);
    break;
  case 3:
    hipLaunchKernelGGL((assign_forces_kernel<3, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), nullptr,
                       pdata, p, lb_particle_force_gpu, prefactor,
                       parts_per_block);
    break;
  case 4:
    hipLaunchKernelGGL((assign_forces_kernel<4, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), nullptr,
                       pdata, p, lb_particle_force_gpu, prefactor,
                       parts_per_block);
    break;
  case 5:
    hipLaunchKernelGGL((assign_forces_kernel<5, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), nullptr,
                       pdata, p, lb_particle_force_gpu, prefactor,
                       parts_per_block);
    break;
  case 6:
    hipLaunchKernelGGL((assign_forces_kernel<6, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), nullptr,
                       pdata, p, lb_particle_force_gpu, prefactor,
                       parts_per_block);
    break;
  case 7:
    hipLaunchKernelGGL((assign_forces_kernel<7, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), nullptr,
                       pdata, p, lb_particle_force_gpu, prefactor,
                       parts_per_block);
    break;
  default:
    break;
  }
  _cuda_check_errors(block, grid, "assign_forces", __FILE__, __LINE__);
}

/* Init the internal data structures of the P3M GPU.
 * Mainly allocation on the device and influence function calculation.
 * Be advised: this needs mesh^3*5*sizeof(REAL_TYPE) of device memory.
 * We use real to complex FFTs, so the size of the reciprocal mesh
 * is (cuFFT convention) Nx x Ny x [ Nz /2 + 1 ].
 */
void p3m_gpu_init(int cao, const int mesh[3], double alpha) {
  espressoSystemInterface.requestParticleStructGpu();

  bool reinit_if = false, mesh_changed = false;
  p3m_gpu_data.n_part = gpu_get_particle_pointer().size();

  if ((p3m_gpu_data_initialized == 0) || (p3m_gpu_data.alpha != alpha)) {
    p3m_gpu_data.alpha = static_cast<REAL_TYPE>(alpha);
    reinit_if = true;
  }

  if ((p3m_gpu_data_initialized == 0) || (p3m_gpu_data.cao != cao)) {
    p3m_gpu_data.cao = cao;
    // NOLINTNEXTLINE(bugprone-integer-division)
    p3m_gpu_data.pos_shift = static_cast<REAL_TYPE>((p3m_gpu_data.cao - 1) / 2);
    reinit_if = true;
  }

  if ((p3m_gpu_data_initialized == 0) || (p3m_gpu_data.mesh[0] != mesh[0]) ||
      (p3m_gpu_data.mesh[1] != mesh[1]) || (p3m_gpu_data.mesh[2] != mesh[2])) {
    std::copy(mesh, mesh + 3, p3m_gpu_data.mesh);
    mesh_changed = true;
    reinit_if = true;
  }

  auto const box_l = box_geo.length();

  if ((p3m_gpu_data_initialized == 0) || (p3m_gpu_data.box[0] != box_l[0]) ||
      (p3m_gpu_data.box[1] != box_l[1]) || (p3m_gpu_data.box[2] != box_l[2])) {
    std::copy(box_l.begin(), box_l.end(), p3m_gpu_data.box);
    reinit_if = true;
  }

  p3m_gpu_data.mesh_z_padded = (mesh[2] / 2 + 1) * 2;
  p3m_gpu_data.mesh_size = mesh[0] * mesh[1] * p3m_gpu_data.mesh_z_padded;

  for (int i = 0; i < 3; i++) {
    p3m_gpu_data.hi[i] =
        static_cast<REAL_TYPE>(p3m_gpu_data.mesh[i]) / p3m_gpu_data.box[i];
  }

  if ((p3m_gpu_data_initialized == 1) && mesh_changed) {
    cuda_safe_mem(cudaFree(p3m_gpu_data.charge_mesh));
    p3m_gpu_data.charge_mesh = nullptr;
    cuda_safe_mem(cudaFree(p3m_gpu_data.force_mesh_x));
    p3m_gpu_data.force_mesh_x = nullptr;
    cuda_safe_mem(cudaFree(p3m_gpu_data.force_mesh_y));
    p3m_gpu_data.force_mesh_y = nullptr;
    cuda_safe_mem(cudaFree(p3m_gpu_data.force_mesh_z));
    p3m_gpu_data.force_mesh_z = nullptr;
    cuda_safe_mem(cudaFree(p3m_gpu_data.G_hat));
    p3m_gpu_data.G_hat = nullptr;

    cufftDestroy(p3m_gpu_fft_plans.forw_plan);
    cufftDestroy(p3m_gpu_fft_plans.back_plan);

    p3m_gpu_data_initialized = 0;
  }

  if ((p3m_gpu_data_initialized == 0) && (p3m_gpu_data.mesh_size > 0)) {
#if defined(__HIPCC__) and not defined(__CUDACC__)
    // rocFFT gives wrong P3M results for mesh sizes whose prime factors are
    // not 2, 3 or 5. So we check for the other supported prime factors and
    // throw an error if they are used with HIP.
    if (mesh[0] % 7 == 0 || mesh[0] % 11 == 0 || mesh[0] % 13 == 0 ||
        mesh[1] % 7 == 0 || mesh[1] % 11 == 0 || mesh[1] % 13 == 0 ||
        mesh[2] % 7 == 0 || mesh[2] % 11 == 0 || mesh[2] % 13 == 0) {
      throw std::runtime_error("Mesh size not supported");
    }
#endif

    /* Size of the complex mesh Nx * Ny * ( Nz / 2 + 1 ) */
    const int cmesh_size = p3m_gpu_data.mesh[0] * p3m_gpu_data.mesh[1] *
                           (p3m_gpu_data.mesh[2] / 2 + 1);

    cuda_safe_mem(cudaMalloc((void **)&(p3m_gpu_data.charge_mesh),
                             cmesh_size * sizeof(FFT_TYPE_COMPLEX)));
    cuda_safe_mem(cudaMalloc((void **)&(p3m_gpu_data.force_mesh_x),
                             cmesh_size * sizeof(FFT_TYPE_COMPLEX)));
    cuda_safe_mem(cudaMalloc((void **)&(p3m_gpu_data.force_mesh_y),
                             cmesh_size * sizeof(FFT_TYPE_COMPLEX)));
    cuda_safe_mem(cudaMalloc((void **)&(p3m_gpu_data.force_mesh_z),
                             cmesh_size * sizeof(FFT_TYPE_COMPLEX)));
    cuda_safe_mem(cudaMalloc((void **)&(p3m_gpu_data.G_hat),
                             cmesh_size * sizeof(REAL_TYPE)));

    if (cufftPlan3d(&(p3m_gpu_fft_plans.forw_plan), mesh[0], mesh[1], mesh[2],
                    FFT_PLAN_FORW_FLAG) != CUFFT_SUCCESS ||
        cufftPlan3d(&(p3m_gpu_fft_plans.back_plan), mesh[0], mesh[1], mesh[2],
                    FFT_PLAN_BACK_FLAG) != CUFFT_SUCCESS) {
      throw std::runtime_error("Unable to create fft plan");
    }
  }

  if ((reinit_if || (p3m_gpu_data_initialized == 0)) &&
      (p3m_gpu_data.mesh_size > 0)) {
    dim3 grid(1, 1, 1);
    dim3 block(1, 1, 1);
    block.x = 512 / mesh[0] + 1;
    block.y = mesh[1];
    block.z = 1;
    grid.x = mesh[0] / block.x + 1;
    grid.z = mesh[2] / 2 + 1;

    switch (p3m_gpu_data.cao) {
    case 1:
      KERNELCALL(calculate_influence_function_device<1>, grid, block,
                 p3m_gpu_data);
      break;
    case 2:
      KERNELCALL(calculate_influence_function_device<2>, grid, block,
                 p3m_gpu_data);
      break;
    case 3:
      KERNELCALL(calculate_influence_function_device<3>, grid, block,
                 p3m_gpu_data);
      break;
    case 4:
      KERNELCALL(calculate_influence_function_device<4>, grid, block,
                 p3m_gpu_data);
      break;
    case 5:
      KERNELCALL(calculate_influence_function_device<5>, grid, block,
                 p3m_gpu_data);
      break;
    case 6:
      KERNELCALL(calculate_influence_function_device<6>, grid, block,
                 p3m_gpu_data);
      break;
    case 7:
      KERNELCALL(calculate_influence_function_device<7>, grid, block,
                 p3m_gpu_data);
      break;
    }
  }
  if (p3m_gpu_data.mesh_size > 0)
    p3m_gpu_data_initialized = 1;
}

/**
 *  \brief The long range part of the P3M algorithm.
 */
void p3m_gpu_add_farfield_force() {
  auto device_particles = gpu_get_particle_pointer();
  CUDA_particle_data *lb_particle_gpu = device_particles.data();
  p3m_gpu_data.n_part = device_particles.size();

  float *lb_particle_force_gpu = gpu_get_particle_force_pointer();

  if (p3m_gpu_data.n_part == 0)
    return;

  dim3 gridConv(p3m_gpu_data.mesh[0], p3m_gpu_data.mesh[1], 1);
  dim3 threadsConv(p3m_gpu_data.mesh[2] / 2 + 1, 1, 1);

  REAL_TYPE prefactor =
      coulomb.prefactor / // NOLINT(bugprone-narrowing-conversions)
      (p3m_gpu_data.box[0] * p3m_gpu_data.box[1] * p3m_gpu_data.box[2] * 2.0);

  cuda_safe_mem(cudaMemset(p3m_gpu_data.charge_mesh, 0,
                           p3m_gpu_data.mesh_size * sizeof(REAL_TYPE)));

  /* Interpolate the charges to the mesh */
  assign_charges(lb_particle_gpu, p3m_gpu_data);

  /* Do forward FFT of the charge mesh */
  if (FFT_FORW_FFT(p3m_gpu_fft_plans.forw_plan,
                   (REAL_TYPE *)p3m_gpu_data.charge_mesh,
                   p3m_gpu_data.charge_mesh) != CUFFT_SUCCESS) {
    fprintf(stderr, "CUFFT error: Forward FFT failed\n");
    return;
  }

  /* Do convolution */
  KERNELCALL(apply_influence_function, gridConv, threadsConv, p3m_gpu_data);

  /* Take derivative */
  KERNELCALL(apply_diff_op, gridConv, threadsConv, p3m_gpu_data);

  /* Transform the components of the electric field back */
  FFT_BACK_FFT(p3m_gpu_fft_plans.back_plan, p3m_gpu_data.force_mesh_x,
               (REAL_TYPE *)p3m_gpu_data.force_mesh_x);
  FFT_BACK_FFT(p3m_gpu_fft_plans.back_plan, p3m_gpu_data.force_mesh_y,
               (REAL_TYPE *)p3m_gpu_data.force_mesh_y);
  FFT_BACK_FFT(p3m_gpu_fft_plans.back_plan, p3m_gpu_data.force_mesh_z,
               (REAL_TYPE *)p3m_gpu_data.force_mesh_z);

  /* Assign the forces from the mesh back to the particles */
  assign_forces(lb_particle_gpu, p3m_gpu_data, lb_particle_force_gpu,
                prefactor);
}

#endif /* ELECTROSTATICS */
