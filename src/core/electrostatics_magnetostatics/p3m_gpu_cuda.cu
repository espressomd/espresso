/*
   Copyright (C) 2010-2018 The ESPResSo project

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

#include "cuda_wrapper.hpp"

/** \file
 *
 *  P3M electrostatics on GPU.
 *
 *  The corresponding header file is p3m_gpu.hpp.
 */

#include "config.hpp"

#ifdef ELECTROSTATICS

#ifdef P3M_GPU_DEBUG
#define P3M_GPU_TRACE(A) A
#else
#define P3M_GPU_TRACE(A)
#endif

#include <stdio.h>
#include <stdlib.h>

#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "cufft_wrapper.hpp"

#include "electrostatics_magnetostatics/p3m_gpu.hpp"

#include "EspressoSystemInterface.hpp"
#include "global.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "p3m_gpu_common.hpp"

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

extern double box_l[3];

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

template <int cao_value, typename T> __device__ T caf(int i, T x) {
  switch (cao_value) {
  case 1:
    return 1.0;
  case 2: {
    switch (i) {
    case 0:
      return 0.5 - x;
    case 1:
      return 0.5 + x;
    default:
      return 0.0;
    }
  }
  case 3: {
    switch (i) {
    case 0:
      return 0.5 * sqr(0.5 - x);
    case 1:
      return 0.75 - sqr(x);
    case 2:
      return 0.5 * sqr(0.5 + x);
    default:
      return 0.0;
    }
  case 4: {
    switch (i) {
    case 0:
      return (1.0 + x * (-6.0 + x * (12.0 - x * 8.0))) / 48.0;
    case 1:
      return (23.0 + x * (-30.0 + x * (-12.0 + x * 24.0))) / 48.0;
    case 2:
      return (23.0 + x * (30.0 + x * (-12.0 - x * 24.0))) / 48.0;
    case 3:
      return (1.0 + x * (6.0 + x * (12.0 + x * 8.0))) / 48.0;
    default:
      return 0.0;
    }
  }
  case 5: {
    switch (i) {
    case 0:
      return (1.0 + x * (-8.0 + x * (24.0 + x * (-32.0 + x * 16.0)))) / 384.0;
    case 1:
      return (19.0 + x * (-44.0 + x * (24.0 + x * (16.0 - x * 16.0)))) / 96.0;
    case 2:
      return (115.0 + x * x * (-120.0 + x * x * 48.0)) / 192.0;
    case 3:
      return (19.0 + x * (44.0 + x * (24.0 + x * (-16.0 - x * 16.0)))) / 96.0;
    case 4:
      return (1.0 + x * (8.0 + x * (24.0 + x * (32.0 + x * 16.0)))) / 384.0;
    default:
      return 0.0;
    }
  }
  case 6: {
    switch (i) {
    case 0:
      return (1.0 +
              x * (-10.0 + x * (40.0 + x * (-80.0 + x * (80.0 - x * 32.0))))) /
             3840.0;
    case 1:
      return (237.0 +
              x * (-750.0 +
                   x * (840.0 + x * (-240.0 + x * (-240.0 + x * 160.0))))) /
             3840.0;
    case 2:
      return (841.0 +
              x * (-770.0 +
                   x * (-440.0 + x * (560.0 + x * (80.0 - x * 160.0))))) /
             1920.0;
    case 3:
      return (841.0 +
              x * (+770.0 +
                   x * (-440.0 + x * (-560.0 + x * (80.0 + x * 160.0))))) /
             1920.0;
    case 4:
      return (237.0 +
              x * (750.0 +
                   x * (840.0 + x * (240.0 + x * (-240.0 - x * 160.0))))) /
             3840.0;
    case 5:
      return (1.0 +
              x * (10.0 + x * (40.0 + x * (80.0 + x * (80.0 + x * 32.0))))) /
             3840.0;
    default:
      return 0.0;
    }
  }
  case 7: {
    switch (i) {
    case 0:
      return (1.0 +
              x * (-12.0 +
                   x * (60.0 + x * (-160.0 +
                                    x * (240.0 + x * (-192.0 + x * 64.0)))))) /
             46080.0;
    case 1:
      return (361.0 + x * (-1416.0 +
                           x * (2220.0 +
                                x * (-1600.0 +
                                     x * (240.0 + x * (384.0 - x * 192.0)))))) /
             23040.0;
    case 2:
      return (10543.0 +
              x * (-17340.0 +
                   x * (4740.0 +
                        x * (6880.0 +
                             x * (-4080.0 + x * (-960.0 + x * 960.0)))))) /
             46080.0;
    case 3:
      return (5887.0 + x * x * (-4620.0 + x * x * (1680.0 - x * x * 320.0))) /
             11520.0;
    case 4:
      return (10543.0 +
              x * (17340.0 +
                   x * (4740.0 +
                        x * (-6880.0 +
                             x * (-4080.0 + x * (960.0 + x * 960.0)))))) /
             46080.0;
    case 5:
      return (361.0 +
              x * (1416.0 +
                   x * (2220.0 +
                        x * (1600.0 +
                             x * (240.0 + x * (-384.0 - x * 192.0)))))) /
             23040.0;
    case 6:
      return (1.0 +
              x * (12.0 +
                   x * (60.0 +
                        x * (160.0 + x * (240.0 + x * (192.0 + x * 64.0)))))) /
             46080.0;
    default:
      return 0.0;
    }
  }
  }
  }
  return 0.0;
}

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
    Leni[i] = 1.0 / p.box[i];
    Meshi[i] = 1.0 / p.mesh[i];
  }

  Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner = 0.0;

  for (MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++) {
    NMX = ((NX > p.mesh[0] / 2) ? NX - p.mesh[0] : NX) + p.mesh[0] * MX;
    S1 = int_pow<2 * cao>(csinc(Meshi[0] * NMX));
    for (MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++) {
      NMY = ((NY > p.mesh[1] / 2) ? NY - p.mesh[1] : NY) + p.mesh[1] * MY;
      S2 = S1 * int_pow<2 * cao>(csinc(Meshi[1] * NMY));
      for (MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++) {
        NMZ = ((NZ > p.mesh[2] / 2) ? NZ - p.mesh[2] : NZ) + p.mesh[2] * MZ;
        S3 = S2 * int_pow<2 * cao>(csinc(Meshi[2] * NMZ));

        NM2 = sqr(NMX * Leni[0]) + sqr(NMY * Leni[1]) + sqr(NMZ * Leni[2]);
        *Nenner += S3;

        TE = exp(-sqr(PI / (p.alpha)) * NM2);
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

  const int NX = blockDim.x * blockIdx.x + threadIdx.x;
  const int NY = blockDim.y * blockIdx.y + threadIdx.y;
  const int NZ = blockDim.z * blockIdx.z + threadIdx.z;
  REAL_TYPE Dnx, Dny, Dnz;
  REAL_TYPE Zaehler[3] = {0.0, 0.0, 0.0}, Nenner = 0.0;
  REAL_TYPE zwi;
  int ind = 0;
  REAL_TYPE Leni[3];
  for (int i = 0; i < 3; ++i)
    Leni[i] = 1.0 / p.box[i];

  if ((NX >= p.mesh[0]) || (NY >= p.mesh[1]) || (NZ >= (p.mesh[2] / 2 + 1)))
    return;

  ind = NX * p.mesh[1] * (p.mesh[2] / 2 + 1) + NY * (p.mesh[2] / 2 + 1) + NZ;

  if ((NX == 0) && (NY == 0) && (NZ == 0))
    p.G_hat[ind] = 0.0;
  else if ((NX % (p.mesh[0] / 2) == 0) && (NY % (p.mesh[1] / 2) == 0) &&
           (NZ % (p.mesh[2] / 2) == 0))
    p.G_hat[ind] = 0.0;
  else {
    Aliasing_sums_ik<cao>(p, NX, NY, NZ, Zaehler, &Nenner);

    Dnx = (NX > p.mesh[0] / 2) ? NX - p.mesh[0] : NX;
    Dny = (NY > p.mesh[1] / 2) ? NY - p.mesh[1] : NY;
    Dnz = (NZ > p.mesh[2] / 2) ? NZ - p.mesh[2] : NZ;

    zwi = Dnx * Zaehler[0] * Leni[0] + Dny * Zaehler[1] * Leni[1] +
          Dnz * Zaehler[2] * Leni[2];
    zwi /= ((sqr(Dnx * Leni[0]) + sqr(Dny * Leni[1]) + sqr(Dnz * Leni[2])) *
            sqr(Nenner));
    p.G_hat[ind] = 2.0 * zwi / PI;
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
  const int linear_index =
      linear_index_k(p, blockIdx.x, blockIdx.y, threadIdx.x);

  const int nx =
      (blockIdx.x > p.mesh[0] / 2) ? blockIdx.x - p.mesh[0] : blockIdx.x;
  const int ny =
      (blockIdx.y > p.mesh[1] / 2) ? blockIdx.y - p.mesh[1] : blockIdx.y;
  const int nz = threadIdx.x;

  const FFT_TYPE_COMPLEX meshw = p.charge_mesh[linear_index];
  FFT_TYPE_COMPLEX buf;
  buf.x = -2.0 * PI * meshw.y;
  buf.y = 2.0 * PI * meshw.x;

  p.force_mesh_x[linear_index].x = nx * buf.x / p.box[0];
  p.force_mesh_x[linear_index].y = nx * buf.y / p.box[0];

  p.force_mesh_y[linear_index].x = ny * buf.x / p.box[1];
  p.force_mesh_y[linear_index].y = ny * buf.y / p.box[1];

  p.force_mesh_z[linear_index].x = nz * buf.x / p.box[2];
  p.force_mesh_z[linear_index].y = nz * buf.y / p.box[2];
}

__device__ inline int wrap_index(const int ind, const int mesh) {
  if (ind < 0)
    return ind + mesh;
  else if (ind >= mesh)
    return ind - mesh;
  else
    return ind;
}

__global__ void apply_influence_function(const P3MGpuData p) {
  const int linear_index =
      linear_index_k(p, blockIdx.x, blockIdx.y, threadIdx.x);

  p.charge_mesh[linear_index].x *= p.G_hat[linear_index];
  p.charge_mesh[linear_index].y *= p.G_hat[linear_index];
}

template <int cao, bool shared>
__global__ void assign_charge_kernel(const CUDA_particle_data *const pdata,
                                     const P3MGpuData par,
                                     const int parts_per_block) {
  const int part_in_block = threadIdx.x / cao;
  const int cao_id_x = threadIdx.x - part_in_block * cao;
  /** id of the particle **/
  int id =
      parts_per_block * (blockIdx.x * gridDim.y + blockIdx.y) + part_in_block;
  if (id >= par.n_part)
    return;
  /** position relative to the closest gird point **/
  REAL_TYPE m_pos[3];
  /** index of the nearest mesh point **/
  int nmp_x, nmp_y, nmp_z;

  const CUDA_particle_data p = pdata[id];
  REAL_TYPE *charge_mesh = (REAL_TYPE *)par.charge_mesh;

  m_pos[0] = p.p[0] * par.hi[0] - par.pos_shift;
  m_pos[1] = p.p[1] * par.hi[1] - par.pos_shift;
  m_pos[2] = p.p[2] * par.hi[2] - par.pos_shift;

  nmp_x = (int)floorf(m_pos[0] + 0.5f);
  nmp_y = (int)floorf(m_pos[1] + 0.5f);
  nmp_z = (int)floorf(m_pos[2] + 0.5f);

  m_pos[0] -= nmp_x;
  m_pos[1] -= nmp_y;
  m_pos[2] -= nmp_z;

  nmp_x = wrap_index(nmp_x + cao_id_x, par.mesh[0]);
  nmp_y = wrap_index(nmp_y + threadIdx.y, par.mesh[1]);
  nmp_z = wrap_index(nmp_z + threadIdx.z, par.mesh[2]);

  const int ind = linear_index_r(par, nmp_x, nmp_y, nmp_z);

  HIP_DYNAMIC_SHARED(float, weights);

  if (shared) {
    if ((threadIdx.y < 3) && (threadIdx.z == 0)) {
      weights[3 * cao * part_in_block + 3 * cao_id_x + threadIdx.y] =
          caf<cao>(cao_id_x, m_pos[threadIdx.y]);
    }

    __syncthreads();

    atomicAdd(&(charge_mesh[ind]),
              weights[3 * cao * part_in_block + 3 * cao_id_x + 0] *
                  weights[3 * cao * part_in_block + 3 * threadIdx.y + 1] *
                  weights[3 * cao * part_in_block + 3 * threadIdx.z + 2] * p.q);

  } else {
    atomicAdd(&(charge_mesh[ind]), caf<cao>(cao_id_x, m_pos[0]) *
                                       caf<cao>(threadIdx.y, m_pos[1]) *
                                       caf<cao>(threadIdx.z, m_pos[2]) * p.q);
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
                       dim3(block), 0, 0, pdata, p, parts_per_block);
    break;
  case 2:
    hipLaunchKernelGGL((assign_charge_kernel<2, false>), dim3(grid),
                       dim3(block), 0, 0, pdata, p, parts_per_block);
    break;
  case 3:
    hipLaunchKernelGGL((assign_charge_kernel<3, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), 0, pdata,
                       p, parts_per_block);
    break;
  case 4:
    hipLaunchKernelGGL((assign_charge_kernel<4, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), 0, pdata,
                       p, parts_per_block);
    break;
  case 5:
    hipLaunchKernelGGL((assign_charge_kernel<5, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), 0, pdata,
                       p, parts_per_block);
    break;
  case 6:
    hipLaunchKernelGGL((assign_charge_kernel<6, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), 0, pdata,
                       p, parts_per_block);
    break;
  case 7:
    hipLaunchKernelGGL((assign_charge_kernel<7, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(REAL_TYPE), 0, pdata,
                       p, parts_per_block);
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
  /** id of the particle **/
  int id =
      parts_per_block * (blockIdx.x * gridDim.y + blockIdx.y) + part_in_block;
  if (id >= par.n_part)
    return;
  /** position relative to the closest grid point **/
  REAL_TYPE m_pos[3];
  /** index of the nearest mesh point **/
  int nmp_x, nmp_y, nmp_z;

  const CUDA_particle_data p = pdata[id];

  m_pos[0] = p.p[0] * par.hi[0] - par.pos_shift;
  m_pos[1] = p.p[1] * par.hi[1] - par.pos_shift;
  m_pos[2] = p.p[2] * par.hi[2] - par.pos_shift;

  nmp_x = (int)floorf(m_pos[0] + 0.5f);
  nmp_y = (int)floorf(m_pos[1] + 0.5f);
  nmp_z = (int)floorf(m_pos[2] + 0.5f);

  m_pos[0] -= nmp_x;
  m_pos[1] -= nmp_y;
  m_pos[2] -= nmp_z;

  nmp_x = wrap_index(nmp_x + cao_id_x, par.mesh[0]);
  nmp_y = wrap_index(nmp_y + threadIdx.y, par.mesh[1]);
  nmp_z = wrap_index(nmp_z + threadIdx.z, par.mesh[2]);

  REAL_TYPE c;
  const int index = linear_index_r(par, nmp_x, nmp_y, nmp_z);

  HIP_DYNAMIC_SHARED(float, weights);

  if (shared) {
    if ((threadIdx.y < 3) && (threadIdx.z == 0)) {
      weights[3 * cao * part_in_block + 3 * cao_id_x + threadIdx.y] =
          caf<cao>(cao_id_x, m_pos[threadIdx.y]);
    }

    __syncthreads();

    c = -prefactor * weights[3 * cao * part_in_block + 3 * cao_id_x + 0] *
        weights[3 * cao * part_in_block + 3 * threadIdx.y + 1] *
        weights[3 * cao * part_in_block + 3 * threadIdx.z + 2] * p.q;
  } else {
    c = -prefactor * caf<cao>(cao_id_x, m_pos[0]) *
        caf<cao>(threadIdx.y, m_pos[1]) * caf<cao>(threadIdx.z, m_pos[2]) * p.q;
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

  if ((n_part % parts_per_block) == 0)
    n_blocks = std::max<int>(1, n_part / parts_per_block);
  else
    n_blocks = n_part / parts_per_block + 1;

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

  /** Switch for assignment templates, the shared version only is faster for cao
   * > 2 */
  switch (cao) {
  case 1:
    hipLaunchKernelGGL((assign_forces_kernel<1, false>), dim3(grid),
                       dim3(block), 0, 0, pdata, p, lb_particle_force_gpu,
                       prefactor, parts_per_block);
    break;
  case 2:
    hipLaunchKernelGGL((assign_forces_kernel<2, false>), dim3(grid),
                       dim3(block), 0, 0, pdata, p, lb_particle_force_gpu,
                       prefactor, parts_per_block);
    break;
  case 3:
    hipLaunchKernelGGL((assign_forces_kernel<3, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), 0, pdata, p,
                       lb_particle_force_gpu, prefactor, parts_per_block);
    break;
  case 4:
    hipLaunchKernelGGL((assign_forces_kernel<4, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), 0, pdata, p,
                       lb_particle_force_gpu, prefactor, parts_per_block);
    break;
  case 5:
    hipLaunchKernelGGL((assign_forces_kernel<5, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), 0, pdata, p,
                       lb_particle_force_gpu, prefactor, parts_per_block);
    break;
  case 6:
    hipLaunchKernelGGL((assign_forces_kernel<6, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), 0, pdata, p,
                       lb_particle_force_gpu, prefactor, parts_per_block);
    break;
  case 7:
    hipLaunchKernelGGL((assign_forces_kernel<7, true>), dim3(grid), dim3(block),
                       3 * parts_per_block * cao * sizeof(float), 0, pdata, p,
                       lb_particle_force_gpu, prefactor, parts_per_block);
    break;
  default:
    break;
  }
  _cuda_check_errors(block, grid, "assign_forces", __FILE__, __LINE__);
}

/* Init the internal datastructures of the P3M GPU.
 * Mainly allocation on the device and influence function calculation.
 * Be advised: this needs mesh^3*5*sizeof(REAL_TYPE) of device memory.
 * We use real to complex FFTs, so the size of the reciprocal mesh
 * is (cuFFT convention) Nx x Ny x [ Nz /2 + 1 ].
 */

void p3m_gpu_init(int cao, const int mesh[3], double alpha) {
  P3M_GPU_TRACE(printf("cao %d mesh %d %d %d, alpha %e, box (%e %e %e)\n", cao,
                       mesh[0], mesh[1], mesh[2], alpha, box_l[0], box_l[1],
                       box_l[2]));
  if (this_node == 0) {
    espressoSystemInterface.requestParticleStructGpu();

    int reinit_if = 0, mesh_changed = 0;
    p3m_gpu_data.n_part =
        gpu_get_global_particle_vars_pointer_host()->number_of_particles;

    if ((p3m_gpu_data_initialized == 0) || (p3m_gpu_data.alpha != alpha)) {
      p3m_gpu_data.alpha = alpha;
      reinit_if = 1;
    }

    if ((p3m_gpu_data_initialized == 0) || (p3m_gpu_data.cao != cao)) {
      p3m_gpu_data.cao = cao;
      p3m_gpu_data.pos_shift = (REAL_TYPE)((p3m_gpu_data.cao - 1) / 2);
      reinit_if = 1;
    }

    if ((p3m_gpu_data_initialized == 0) || (p3m_gpu_data.mesh[0] != mesh[0]) ||
        (p3m_gpu_data.mesh[1] != mesh[1]) ||
        (p3m_gpu_data.mesh[2] != mesh[2])) {
      std::copy(mesh, mesh + 3, p3m_gpu_data.mesh);
      mesh_changed = 1;
      reinit_if = 1;
    }

    if ((p3m_gpu_data_initialized == 0) || (p3m_gpu_data.box[0] != box_l[0]) ||
        (p3m_gpu_data.box[1] != box_l[1]) ||
        (p3m_gpu_data.box[2] != box_l[2])) {
      std::copy(box_l, box_l + 3, p3m_gpu_data.box);
      reinit_if = 1;
    }

    p3m_gpu_data.mesh_z_padded = (mesh[2] / 2 + 1) * 2;
    p3m_gpu_data.mesh_size = mesh[0] * mesh[1] * p3m_gpu_data.mesh_z_padded;

    for (int i = 0; i < 3; i++) {
      p3m_gpu_data.hi[i] = p3m_gpu_data.mesh[i] / p3m_gpu_data.box[i];
    }

    if ((p3m_gpu_data_initialized == 1) && (mesh_changed == 1)) {
      cuda_safe_mem(cudaFree(p3m_gpu_data.charge_mesh));
      p3m_gpu_data.charge_mesh = 0;
      cuda_safe_mem(cudaFree(p3m_gpu_data.force_mesh_x));
      p3m_gpu_data.force_mesh_x = 0;
      cuda_safe_mem(cudaFree(p3m_gpu_data.force_mesh_y));
      p3m_gpu_data.force_mesh_y = 0;
      cuda_safe_mem(cudaFree(p3m_gpu_data.force_mesh_z));
      p3m_gpu_data.force_mesh_z = 0;
      cuda_safe_mem(cudaFree(p3m_gpu_data.G_hat));
      p3m_gpu_data.G_hat = 0;

      cufftDestroy(p3m_gpu_fft_plans.forw_plan);
      cufftDestroy(p3m_gpu_fft_plans.back_plan);

      p3m_gpu_data_initialized = 0;
    }

    if ((p3m_gpu_data_initialized == 0) && (p3m_gpu_data.mesh_size > 0)) {
      /** Size of the complex mesh Nx * Ny * ( Nz / 2 + 1 ) */
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
        throw std::string("Unable to create fft plan");
      }
    }

    if (((reinit_if == 1) || (p3m_gpu_data_initialized == 0)) &&
        (p3m_gpu_data.mesh_size > 0)) {
      dim3 grid(1, 1, 1);
      dim3 block(1, 1, 1);
      block.x = 512 / mesh[0] + 1;
      block.y = mesh[1];
      block.z = 1;
      grid.x = mesh[0] / block.x + 1;
      grid.z = mesh[2] / 2 + 1;

      P3M_GPU_TRACE(printf("mesh %d %d %d, grid (%d %d %d), block (%d %d %d)\n",
                           mesh[0], mesh[1], mesh[2], grid.x, grid.y, grid.z,
                           block.x, block.y, block.z));

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
}

/**
 *  \brief The long range part of the P3M algorithm.
 */
void p3m_gpu_add_farfield_force() {
  CUDA_particle_data *lb_particle_gpu;
  float *lb_particle_force_gpu;

  lb_particle_gpu = gpu_get_particle_pointer();
  lb_particle_force_gpu = gpu_get_particle_force_pointer();

  p3m_gpu_data.n_part =
      gpu_get_global_particle_vars_pointer_host()->number_of_particles;

  if (p3m_gpu_data.n_part == 0)
    return;

  dim3 gridConv(p3m_gpu_data.mesh[0], p3m_gpu_data.mesh[1], 1);
  dim3 threadsConv(p3m_gpu_data.mesh[2] / 2 + 1, 1, 1);

  REAL_TYPE prefactor =
      coulomb.prefactor /
      (p3m_gpu_data.box[0] * p3m_gpu_data.box[1] * p3m_gpu_data.box[2] * 2.0);

  cuda_safe_mem(cudaMemset(p3m_gpu_data.charge_mesh, 0,
                           p3m_gpu_data.mesh_size * sizeof(REAL_TYPE)));

  /** Interpolate the charges to the mesh */
  assign_charges(lb_particle_gpu, p3m_gpu_data);

  /** Do forward FFT of the charge mesh */
  if (FFT_FORW_FFT(p3m_gpu_fft_plans.forw_plan,
                   (REAL_TYPE *)p3m_gpu_data.charge_mesh,
                   p3m_gpu_data.charge_mesh) != CUFFT_SUCCESS) {
    fprintf(stderr, "CUFFT error: Forward FFT failed\n");
    return;
  }

  /** Do convolution */
  KERNELCALL(apply_influence_function, gridConv, threadsConv, p3m_gpu_data);

  /** Take derivative */
  KERNELCALL(apply_diff_op, gridConv, threadsConv, p3m_gpu_data);

  /** Transform the components of the electric field back */
  FFT_BACK_FFT(p3m_gpu_fft_plans.back_plan, p3m_gpu_data.force_mesh_x,
               (REAL_TYPE *)p3m_gpu_data.force_mesh_x);
  FFT_BACK_FFT(p3m_gpu_fft_plans.back_plan, p3m_gpu_data.force_mesh_y,
               (REAL_TYPE *)p3m_gpu_data.force_mesh_y);
  FFT_BACK_FFT(p3m_gpu_fft_plans.back_plan, p3m_gpu_data.force_mesh_z,
               (REAL_TYPE *)p3m_gpu_data.force_mesh_z);

  /** Assign the forces from the mesh back to the particles */
  assign_forces(lb_particle_gpu, p3m_gpu_data, lb_particle_force_gpu,
                prefactor);
}

#endif /* ELECTROSTATICS */
