/*
 * Copyright (C) 2015-2022 The ESPResSo project
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
/** @file
 *
 *  P3M electrostatics on GPU.
 *
 *  The corresponding header file is p3m_gpu_error.hpp.
 */

#include "p3m_gpu_error.hpp"

#include "cuda/utils.cuh"

#include <utils/math/int_pow.hpp>
#include <utils/math/sinc.hpp>
#include <utils/math/sqr.hpp>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>

#include <cuda.h>

#include <algorithm>
#include <cstddef>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

using Utils::int_pow;
using Utils::sqr;

/** @todo Extend to higher order. This comes from some 1/sin expansion in
 *  @cite hockney88a
 */

template <int cao>
__device__ static double p3m_analytic_cotangent_sum(int n, double mesh_i) {
  const double c = sqr(cos(Utils::pi() * mesh_i * n));

  switch (cao) {
  case 1:
    return 1;
  case 2:
    return (1.0 + c * 2.0) / 3.0;
  case 3:
    return (2.0 + c * (11.0 + c * 2.0)) / 15.0;
  case 4:
    return (17.0 + c * (180.0 + c * (114.0 + c * 4.0))) / 315.0;
  case 5:
    return (62.0 + c * (1072.0 + c * (1452.0 + c * (247.0 + c * 2.0)))) /
           2835.0;
  case 6:
    return (1382.0 +
            c * (35396.0 +
                 c * (83021.0 + c * (34096.0 + c * (2026.0 + c * 4.0))))) /
           155925.0;
  case 7:
    return (21844.0 +
            c * (776661.0 +
                 c * (2801040.0 +
                      c * (2123860.0 +
                           c * (349500.0 + c * (8166.0 + c * 4.0)))))) /
           6081075.0;
  }

  return 0.0;
}

template <int cao>
__global__ void p3m_k_space_error_gpu_kernel_ik(int3 mesh, double3 meshi,
                                                double alpha_L, double *he_q) {
  const int nx =
      -mesh.x / 2 + static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
  const int ny =
      -mesh.y / 2 + static_cast<int>(blockDim.y * blockIdx.y + threadIdx.y);
  const int nz =
      -mesh.z / 2 + static_cast<int>(blockDim.z * blockIdx.z + threadIdx.z);

  if ((nx >= mesh.x / 2) || (ny >= mesh.y / 2) || (nz >= mesh.z / 2))
    return;

  const int lind = (nx + mesh.x / 2) * mesh.y * mesh.z +
                   (ny + mesh.y / 2) * mesh.z + (nz + mesh.z / 2);

  if ((nx != 0) || (ny != 0) || (nz != 0)) {
    const double alpha_L_i = 1. / alpha_L;
    const double n2 = sqr(nx) + sqr(ny) + sqr(nz);
    const double cs = p3m_analytic_cotangent_sum<cao>(nz, meshi.z) *
                      p3m_analytic_cotangent_sum<cao>(nx, meshi.x) *
                      p3m_analytic_cotangent_sum<cao>(ny, meshi.y);
    const double ex = exp(-sqr(Utils::pi() * alpha_L_i) * n2);
    const double ex2 = sqr(ex);
    const double U2 =
        int_pow<2 * cao>(Utils::sinc(meshi.x * nx) * Utils::sinc(meshi.y * ny) *
                         Utils::sinc(meshi.z * nz));
    auto const alias1 = ex2 / n2;
    auto const d = alias1 - sqr(U2 * ex / cs) / n2;

    if (d > 0 && (d / alias1 > ROUND_ERROR_PREC))
      he_q[lind] = d;
  } else {
    he_q[lind] = 0;
  }
}

double p3m_k_space_error_gpu(double prefactor, const int *mesh, int cao,
                             int npart, double sum_q2, double alpha_L,
                             const double *box) {
  static thrust::device_vector<double> he_q;
  he_q.resize(static_cast<unsigned>(mesh[0] * mesh[1] * mesh[2]));

  dim3 grid(std::max<unsigned>(1u, static_cast<unsigned>(mesh[0] / 8 + 1)),
            std::max<unsigned>(1u, static_cast<unsigned>(mesh[1] / 8 + 1)),
            std::max<unsigned>(1u, static_cast<unsigned>(mesh[2] / 8 + 1)));

  dim3 block(8, 8, 8);

  int3 mesh3;
  mesh3.x = mesh[0];
  mesh3.y = mesh[1];
  mesh3.z = mesh[2];

  double3 meshi;
  meshi.x = 1. / mesh[0];
  meshi.y = 1. / mesh[1];
  meshi.z = 1. / mesh[2];

  switch (cao) {
  case 1:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<1>, grid, block, mesh3, meshi,
               alpha_L, thrust::raw_pointer_cast(he_q.data()));
    break;
  case 2:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<2>, grid, block, mesh3, meshi,
               alpha_L, thrust::raw_pointer_cast(he_q.data()));
    break;
  case 3:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<3>, grid, block, mesh3, meshi,
               alpha_L, thrust::raw_pointer_cast(he_q.data()));
    break;
  case 4:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<4>, grid, block, mesh3, meshi,
               alpha_L, thrust::raw_pointer_cast(he_q.data()));
    break;
  case 5:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<5>, grid, block, mesh3, meshi,
               alpha_L, thrust::raw_pointer_cast(he_q.data()));
    break;
  case 6:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<6>, grid, block, mesh3, meshi,
               alpha_L, thrust::raw_pointer_cast(he_q.data()));
    break;
  case 7:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<7>, grid, block, mesh3, meshi,
               alpha_L, thrust::raw_pointer_cast(he_q.data()));
    break;
  }

  auto const he_q_final = thrust::reduce(he_q.begin(), he_q.end());

  return 2 * prefactor * sum_q2 * sqrt(he_q_final / npart) / (box[1] * box[2]);
}
