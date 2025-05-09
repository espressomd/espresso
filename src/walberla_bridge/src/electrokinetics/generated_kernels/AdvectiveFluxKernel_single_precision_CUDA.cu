//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file AdvectiveFluxKernel_single_precision_CUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 0aab9c0af2335b1f6fec75deae06e514ccb233ab

#include <cmath>

#include "AdvectiveFluxKernel_single_precision_CUDA.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX __global__

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning push
#pragma warning(disable : 1599)
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_advectivefluxkernel_single_precision_cuda_advectivefluxkernel_single_precision_cuda {
static FUNC_PREFIX __launch_bounds__(256) void advectivefluxkernel_single_precision_cuda_advectivefluxkernel_single_precision_cuda(float *RESTRICT const _data_j, float *RESTRICT const _data_rho, float *RESTRICT const _data_u, int64_t const _size_j_0, int64_t const _size_j_1, int64_t const _size_j_2, int64_t const _stride_j_0, int64_t const _stride_j_1, int64_t const _stride_j_2, int64_t const _stride_j_3, int64_t const _stride_rho_0, int64_t const _stride_rho_1, int64_t const _stride_rho_2, int64_t const _stride_u_0, int64_t const _stride_u_1, int64_t const _stride_u_2, int64_t const _stride_u_3) {
  if (blockDim.y * blockIdx.y + threadIdx.y < _size_j_1 && blockDim.z * blockIdx.z + threadIdx.z < _size_j_2 && blockDim.x * blockIdx.x + threadIdx.x + 1 < _size_j_0) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x + 1;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z;
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_1 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] = -xi_1 * xi_64 * (1.0f - fabs(xi_25)) * (1.0f - fabs(xi_5)) * ((xi_1 < 0.0f) ? (1.0f) : (0.0f)) + xi_39 - xi_40 * xi_50 * (1.0f - fabs(xi_29)) * (1.0f - fabs(xi_7)) * ((xi_40 > 0.0f) ? (1.0f) : (0.0f));
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3] = -xi_21 * xi_62 * (1.0f - fabs(xi_17)) * (1.0f - fabs(xi_45)) * ((xi_62 > 0.0f) ? (1.0f) : (0.0f)) - xi_5 * xi_64 * (1.0f - fabs(xi_1)) * (1.0f - fabs(xi_25)) * ((xi_5 < 0.0f) ? (1.0f) : (0.0f)) + xi_53;
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_0 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3] = -xi_22 * xi_6 * (1.0f - fabs(xi_23)) * (1.0f - fabs(xi_24)) * ((xi_6 > 0.0f) ? (1.0f) : (0.0f)) - xi_25 * xi_64 * (1.0f - fabs(xi_1)) * (1.0f - fabs(xi_5)) * ((xi_25 < 0.0f) ? (1.0f) : (0.0f)) + xi_59;
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_2 < _size_j_2 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = xi_1 * xi_5 * xi_64 * (1.0f - fabs(xi_25)) * ((xi_1 < 0.0f && xi_5 < 0.0f) ? (1.0f) : (0.0f)) - xi_14 * xi_51 * xi_65 * (1.0f - fabs(xi_2)) * ((xi_14 > 0.0f && xi_65 > 0.0f) ? (1.0f) : (0.0f)) + xi_8;
    }
    if (ctr_2 > 0 && ctr_1 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = -xi_1 * xi_5 * xi_64 * (1.0f - fabs(xi_25)) * ((xi_5 > 0.0f && xi_1 < 0.0f) ? (1.0f) : (0.0f)) + xi_27 + xi_41 * xi_61 * xi_63 * (1.0f - fabs(xi_44)) * ((xi_61 > 0.0f && xi_63 < 0.0f) ? (1.0f) : (0.0f));
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_1 < _size_j_1 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] = xi_0 + xi_1 * xi_25 * xi_64 * (1.0f - fabs(xi_5)) * ((xi_1 < 0.0f && xi_25 < 0.0f) ? (1.0f) : (0.0f)) - xi_15 * xi_38 * xi_67 * (1.0f - fabs(xi_28)) * ((xi_15 > 0.0f && xi_38 > 0.0f) ? (1.0f) : (0.0f));
    }
    if (ctr_1 > 0 && ctr_1 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] = -xi_1 * xi_25 * xi_64 * (1.0f - fabs(xi_5)) * ((xi_25 > 0.0f && xi_1 < 0.0f) ? (1.0f) : (0.0f)) + xi_16 * xi_3 * xi_66 * (1.0f - fabs(xi_60)) * ((xi_16 > 0.0f && xi_66 < 0.0f) ? (1.0f) : (0.0f)) + xi_48;
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_0 < _size_j_0 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = xi_12 + xi_25 * xi_5 * xi_64 * (1.0f - fabs(xi_1)) * ((xi_25 < 0.0f && xi_5 < 0.0f) ? (1.0f) : (0.0f)) - xi_33 * xi_34 * xi_42 * (1.0f - fabs(xi_30)) * ((xi_33 > 0.0f && xi_34 > 0.0f) ? (1.0f) : (0.0f));
    }
    if (ctr_1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = -xi_25 * xi_5 * xi_64 * (1.0f - fabs(xi_1)) * ((xi_25 > 0.0f && xi_5 < 0.0f) ? (1.0f) : (0.0f)) + xi_47 * xi_52 * xi_58 * (1.0f - fabs(xi_46)) * ((xi_52 > 0.0f && xi_58 < 0.0f) ? (1.0f) : (0.0f)) + xi_54;
    }
    if (ctr_1 > 0 && ctr_2 > 0) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = -xi_1 * xi_25 * xi_5 * xi_64 * ((xi_1 < 0.0f && xi_25 < 0.0f && xi_5 < 0.0f) ? (1.0f) : (0.0f)) + xi_18 - xi_19 * xi_31 * xi_32 * xi_49 * ((xi_19 > 0.0f && xi_32 > 0.0f && xi_49 > 0.0f) ? (1.0f) : (0.0f));
    }
    if (ctr_1 > 0 && ctr_2 < _size_j_2 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = xi_1 * xi_25 * xi_5 * xi_64 * ((xi_25 > 0.0f && xi_1 < 0.0f && xi_5 < 0.0f) ? (1.0f) : (0.0f)) + xi_10 * xi_20 * xi_36 * xi_68 * ((xi_36 > 0.0f && xi_68 > 0.0f && xi_20 < 0.0f) ? (1.0f) : (0.0f)) + xi_9;
    }
    if (ctr_2 > 0 && ctr_1 < _size_j_1 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = xi_1 * xi_25 * xi_5 * xi_64 * ((xi_5 > 0.0f && xi_1 < 0.0f && xi_25 < 0.0f) ? (1.0f) : (0.0f)) + xi_26 * xi_35 * xi_37 * xi_56 * ((xi_26 > 0.0f && xi_56 > 0.0f && xi_35 < 0.0f) ? (1.0f) : (0.0f)) + xi_43;
    }
    if (ctr_1 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
      const float xi_0 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const float xi_1 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_3 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_4 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_5 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_6 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const float xi_9 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const float xi_10 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_11 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const float xi_13 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_15 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_16 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_17 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const float xi_19 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_20 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_21 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_22 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_23 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_24 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_25 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_26 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_27 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const float xi_28 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_30 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_31 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_33 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_34 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_35 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const float xi_36 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_37 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_38 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_39 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const float xi_40 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const float xi_41 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_42 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_43 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const float xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const float xi_46 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const float xi_47 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const float xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const float xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const float xi_50 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_51 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const float xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_53 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const float xi_54 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const float xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_56 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const float xi_57 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const float xi_58 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_59 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const float xi_60 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const float xi_61 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const float xi_62 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_63 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_64 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const float xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const float xi_66 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const float xi_67 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const float xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = -xi_1 * xi_25 * xi_5 * xi_64 * ((xi_25 > 0.0f && xi_5 > 0.0f && xi_1 < 0.0f) ? (1.0f) : (0.0f)) - xi_11 * xi_13 * xi_4 * xi_55 * ((xi_4 > 0.0f && xi_11 < 0.0f && xi_55 < 0.0f) ? (1.0f) : (0.0f)) + xi_57;
    }
  }
}
} // namespace internal_advectivefluxkernel_single_precision_cuda_advectivefluxkernel_single_precision_cuda

void AdvectiveFluxKernel_single_precision_CUDA::run(IBlock *block, gpuStream_t stream) {

  auto rho = block->getData<gpu::GPUField<float>>(rhoID);
  auto u = block->getData<gpu::GPUField<float>>(uID);
  auto j = block->getData<gpu::GPUField<float>>(jID);

  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()))
  float *RESTRICT const _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho->nrOfGhostLayers()))
  float *RESTRICT const _data_rho = rho->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(u->nrOfGhostLayers()))
  float *RESTRICT const _data_u = u->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(u->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(int64_c(j->xSize()) + 2))
  const int64_t _size_j_0 = int64_t(int64_c(j->xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(int64_c(j->ySize()) + 2))
  const int64_t _size_j_1 = int64_t(int64_c(j->ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(int64_c(j->zSize()) + 2))
  const int64_t _size_j_2 = int64_t(int64_c(j->zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  const int64_t _stride_u_0 = int64_t(u->xStride());
  const int64_t _stride_u_1 = int64_t(u->yStride());
  const int64_t _stride_u_2 = int64_t(u->zStride());
  const int64_t _stride_u_3 = int64_t(1 * int64_t(u->fStride()));
  dim3 _block(uint32_c(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)), uint32_c(((1024 < ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))) ? 1024 : ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))), uint32_c(((64 < ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))) ? 64 : ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))))));
  dim3 _grid(uint32_c(((_size_j_0 - 1) % (((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)) == 0 ? (int64_t)(_size_j_0 - 1) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)) : ((int64_t)(_size_j_0 - 1) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))) + 1)), uint32_c(((_size_j_1) % (((1024 < ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))) ? 1024 : ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))) == 0 ? (int64_t)(_size_j_1) / (int64_t)(((1024 < ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))) ? 1024 : ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))) : ((int64_t)(_size_j_1) / (int64_t)(((1024 < ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))) ? 1024 : ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) + 1)), uint32_c(((_size_j_2) % (((64 < ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))) ? 64 : ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))))) == 0 ? (int64_t)(_size_j_2) / (int64_t)(((64 < ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))) ? 64 : ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))))) : ((int64_t)(_size_j_2) / (int64_t)(((64 < ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))) ? 64 : ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))))) + 1)));
  internal_advectivefluxkernel_single_precision_cuda_advectivefluxkernel_single_precision_cuda::advectivefluxkernel_single_precision_cuda_advectivefluxkernel_single_precision_cuda<<<_grid, _block, 0, stream>>>(_data_j, _data_rho, _data_u, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2, _stride_u_0, _stride_u_1, _stride_u_2, _stride_u_3);
}

void AdvectiveFluxKernel_single_precision_CUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block, gpuStream_t stream) {

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto rho = block->getData<gpu::GPUField<float>>(rhoID);
  auto u = block->getData<gpu::GPUField<float>>(uID);
  auto j = block->getData<gpu::GPUField<float>>(jID);

  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()))
  float *RESTRICT const _data_j = j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho->nrOfGhostLayers()))
  float *RESTRICT const _data_rho = rho->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(u->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(u->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(u->nrOfGhostLayers()))
  float *RESTRICT const _data_u = u->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(u->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
  const int64_t _size_j_0 = int64_t(int64_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
  const int64_t _size_j_1 = int64_t(int64_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
  const int64_t _size_j_2 = int64_t(int64_c(ci.zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  const int64_t _stride_u_0 = int64_t(u->xStride());
  const int64_t _stride_u_1 = int64_t(u->yStride());
  const int64_t _stride_u_2 = int64_t(u->zStride());
  const int64_t _stride_u_3 = int64_t(1 * int64_t(u->fStride()));
  dim3 _block(uint32_c(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)), uint32_c(((1024 < ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))) ? 1024 : ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))), uint32_c(((64 < ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))) ? 64 : ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))))));
  dim3 _grid(uint32_c(((_size_j_0 - 1) % (((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)) == 0 ? (int64_t)(_size_j_0 - 1) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)) : ((int64_t)(_size_j_0 - 1) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))) + 1)), uint32_c(((_size_j_1) % (((1024 < ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))) ? 1024 : ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))) == 0 ? (int64_t)(_size_j_1) / (int64_t)(((1024 < ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))) ? 1024 : ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))) : ((int64_t)(_size_j_1) / (int64_t)(((1024 < ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))) ? 1024 : ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) + 1)), uint32_c(((_size_j_2) % (((64 < ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))) ? 64 : ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))))) == 0 ? (int64_t)(_size_j_2) / (int64_t)(((64 < ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))) ? 64 : ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))))) : ((int64_t)(_size_j_2) / (int64_t)(((64 < ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))) ? 64 : ((_size_j_2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1))))))) ? _size_j_2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1) * ((_size_j_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))) ? _size_j_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 1) ? 128 : _size_j_0 - 1)))))))))) + 1)));
  internal_advectivefluxkernel_single_precision_cuda_advectivefluxkernel_single_precision_cuda::advectivefluxkernel_single_precision_cuda_advectivefluxkernel_single_precision_cuda<<<_grid, _block, 0, stream>>>(_data_j, _data_rho, _data_u, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2, _stride_u_0, _stride_u_1, _stride_u_2, _stride_u_3);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
