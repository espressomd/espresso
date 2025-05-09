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
//! \\file AdvectiveFluxKernel_double_precision_CUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 0aab9c0af2335b1f6fec75deae06e514ccb233ab

#include <cmath>

#include "AdvectiveFluxKernel_double_precision_CUDA.h"
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

namespace internal_advectivefluxkernel_double_precision_cuda_advectivefluxkernel_double_precision_cuda {
static FUNC_PREFIX __launch_bounds__(256) void advectivefluxkernel_double_precision_cuda_advectivefluxkernel_double_precision_cuda(double *RESTRICT const _data_j, double *RESTRICT const _data_rho, double *RESTRICT const _data_u, int64_t const _size_j_0, int64_t const _size_j_1, int64_t const _size_j_2, int64_t const _stride_j_0, int64_t const _stride_j_1, int64_t const _stride_j_2, int64_t const _stride_j_3, int64_t const _stride_rho_0, int64_t const _stride_rho_1, int64_t const _stride_rho_2, int64_t const _stride_u_0, int64_t const _stride_u_1, int64_t const _stride_u_2, int64_t const _stride_u_3) {
  if (blockDim.y * blockIdx.y + threadIdx.y < _size_j_1 && blockDim.z * blockIdx.z + threadIdx.z < _size_j_2 && blockDim.x * blockIdx.x + threadIdx.x + 1 < _size_j_0) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x + 1;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z;
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_1 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] = -xi_17 * xi_7 * (1.0 - fabs(xi_62)) * (1.0 - fabs(xi_65)) * ((xi_7 > 0.0) ? (1.0) : (0.0)) - xi_19 * xi_28 * (1.0 - fabs(xi_3)) * (1.0 - fabs(xi_40)) * ((xi_28 < 0.0) ? (1.0) : (0.0)) + xi_67;
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3] = -xi_16 * xi_50 * (1.0 - fabs(xi_56)) * (1.0 - fabs(xi_63)) * ((xi_50 > 0.0) ? (1.0) : (0.0)) - xi_19 * xi_40 * (1.0 - fabs(xi_28)) * (1.0 - fabs(xi_3)) * ((xi_40 < 0.0) ? (1.0) : (0.0)) + xi_38;
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_0 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3] = xi_15 - xi_19 * xi_3 * (1.0 - fabs(xi_28)) * (1.0 - fabs(xi_40)) * ((xi_3 < 0.0) ? (1.0) : (0.0)) - xi_25 * xi_52 * (1.0 - fabs(xi_60)) * (1.0 - fabs(xi_66)) * ((xi_52 > 0.0) ? (1.0) : (0.0));
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_2 < _size_j_2 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = xi_19 * xi_28 * xi_40 * (1.0 - fabs(xi_3)) * ((xi_28 < 0.0 && xi_40 < 0.0) ? (1.0) : (0.0)) - xi_29 * xi_39 * xi_51 * (1.0 - fabs(xi_53)) * ((xi_29 > 0.0 && xi_51 > 0.0) ? (1.0) : (0.0)) + xi_4;
    }
    if (ctr_2 > 0 && ctr_1 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = -xi_19 * xi_28 * xi_40 * (1.0 - fabs(xi_3)) * ((xi_40 > 0.0 && xi_28 < 0.0) ? (1.0) : (0.0)) + xi_20 * xi_30 * xi_33 * (1.0 - fabs(xi_34)) * ((xi_30 > 0.0 && xi_33 < 0.0) ? (1.0) : (0.0)) + xi_8;
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_1 < _size_j_1 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] = -xi_13 * xi_32 * xi_5 * (1.0 - fabs(xi_46)) * ((xi_13 > 0.0 && xi_32 > 0.0) ? (1.0) : (0.0)) + xi_19 * xi_28 * xi_3 * (1.0 - fabs(xi_40)) * ((xi_28 < 0.0 && xi_3 < 0.0) ? (1.0) : (0.0)) + xi_35;
    }
    if (ctr_1 > 0 && ctr_1 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] = xi_12 - xi_19 * xi_28 * xi_3 * (1.0 - fabs(xi_40)) * ((xi_3 > 0.0 && xi_28 < 0.0) ? (1.0) : (0.0)) + xi_22 * xi_31 * xi_37 * (1.0 - fabs(xi_43)) * ((xi_31 > 0.0 && xi_37 < 0.0) ? (1.0) : (0.0));
    }
    if (ctr_1 > 0 && ctr_2 > 0 && ctr_0 < _size_j_0 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = -xi_10 * xi_41 * xi_54 * (1.0 - fabs(xi_45)) * ((xi_10 > 0.0 && xi_54 > 0.0) ? (1.0) : (0.0)) + xi_18 + xi_19 * xi_3 * xi_40 * (1.0 - fabs(xi_28)) * ((xi_3 < 0.0 && xi_40 < 0.0) ? (1.0) : (0.0));
    }
    if (ctr_1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = -xi_19 * xi_3 * xi_40 * (1.0 - fabs(xi_28)) * ((xi_3 > 0.0 && xi_40 < 0.0) ? (1.0) : (0.0)) + xi_21 * xi_24 * xi_9 * (1.0 - fabs(xi_57)) * ((xi_21 > 0.0 && xi_9 < 0.0) ? (1.0) : (0.0)) + xi_61;
    }
    if (ctr_1 > 0 && ctr_2 > 0) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = xi_1 - xi_19 * xi_28 * xi_3 * xi_40 * ((xi_28 < 0.0 && xi_3 < 0.0 && xi_40 < 0.0) ? (1.0) : (0.0)) - xi_36 * xi_44 * xi_47 * xi_68 * ((xi_44 > 0.0 && xi_47 > 0.0 && xi_68 > 0.0) ? (1.0) : (0.0));
    }
    if (ctr_1 > 0 && ctr_2 < _size_j_2 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = xi_11 * xi_14 * xi_23 * xi_42 * ((xi_14 > 0.0 && xi_23 > 0.0 && xi_42 < 0.0) ? (1.0) : (0.0)) + xi_19 * xi_28 * xi_3 * xi_40 * ((xi_3 > 0.0 && xi_28 < 0.0 && xi_40 < 0.0) ? (1.0) : (0.0)) + xi_26;
    }
    if (ctr_2 > 0 && ctr_1 < _size_j_1 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = xi_0 * xi_55 * xi_59 * xi_64 * ((xi_0 > 0.0 && xi_55 > 0.0 && xi_59 < 0.0) ? (1.0) : (0.0)) + xi_19 * xi_28 * xi_3 * xi_40 * ((xi_40 > 0.0 && xi_28 < 0.0 && xi_3 < 0.0) ? (1.0) : (0.0)) + xi_48;
    }
    if (ctr_1 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
      const double xi_0 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_1 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3];
      const double xi_2 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_3 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_4 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3];
      const double xi_5 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_6 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3];
      const double xi_7 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_8 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3];
      const double xi_9 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_10 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_11 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_12 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3];
      const double xi_13 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_14 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_15 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3];
      const double xi_16 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_17 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_18 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3];
      const double xi_19 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
      const double xi_20 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_21 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_22 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_23 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_24 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_25 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_26 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3];
      const double xi_27 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_28 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2];
      const double xi_29 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_30 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_31 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_32 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_33 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_34 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_35 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3];
      const double xi_36 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_37 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_38 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3];
      const double xi_39 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2];
      const double xi_40 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_41 = _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_42 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_43 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_2 + _stride_u_3];
      const double xi_44 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_45 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_46 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_47 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_48 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3];
      const double xi_49 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2 + 2 * _stride_u_3];
      const double xi_50 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_51 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_52 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_53 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_54 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_55 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + 2 * _stride_u_3];
      const double xi_56 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_57 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 + _stride_u_2];
      const double xi_58 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2];
      const double xi_59 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_60 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2 + _stride_u_3];
      const double xi_61 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3];
      const double xi_62 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + 2 * _stride_u_3];
      const double xi_63 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2];
      const double xi_64 = _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2];
      const double xi_65 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 + _stride_u_3];
      const double xi_66 = _data_u[_stride_u_0 * ctr_0 + _stride_u_1 * ctr_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      const double xi_67 = _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2];
      const double xi_68 = _data_u[_stride_u_0 * ctr_0 - _stride_u_0 + _stride_u_1 * ctr_1 - _stride_u_1 + _stride_u_2 * ctr_2 - _stride_u_2];
      _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = -xi_19 * xi_28 * xi_3 * xi_40 * ((xi_3 > 0.0 && xi_40 > 0.0 && xi_28 < 0.0) ? (1.0) : (0.0)) - xi_2 * xi_27 * xi_49 * xi_58 * ((xi_2 > 0.0 && xi_27 < 0.0 && xi_49 < 0.0) ? (1.0) : (0.0)) + xi_6;
    }
  }
}
} // namespace internal_advectivefluxkernel_double_precision_cuda_advectivefluxkernel_double_precision_cuda

void AdvectiveFluxKernel_double_precision_CUDA::run(IBlock *block, gpuStream_t stream) {

  auto j = block->getData<gpu::GPUField<double>>(jID);
  auto u = block->getData<gpu::GPUField<double>>(uID);
  auto rho = block->getData<gpu::GPUField<double>>(rhoID);

  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()))
  double *RESTRICT const _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho->nrOfGhostLayers()))
  double *RESTRICT const _data_rho = rho->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(u->nrOfGhostLayers()))
  double *RESTRICT const _data_u = u->dataAt(-1, -1, -1, 0);
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
  internal_advectivefluxkernel_double_precision_cuda_advectivefluxkernel_double_precision_cuda::advectivefluxkernel_double_precision_cuda_advectivefluxkernel_double_precision_cuda<<<_grid, _block, 0, stream>>>(_data_j, _data_rho, _data_u, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2, _stride_u_0, _stride_u_1, _stride_u_2, _stride_u_3);
}

void AdvectiveFluxKernel_double_precision_CUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block, gpuStream_t stream) {

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto j = block->getData<gpu::GPUField<double>>(jID);
  auto u = block->getData<gpu::GPUField<double>>(uID);
  auto rho = block->getData<gpu::GPUField<double>>(rhoID);

  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()))
  double *RESTRICT const _data_j = j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho->nrOfGhostLayers()))
  double *RESTRICT const _data_rho = rho->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(u->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(u->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(u->nrOfGhostLayers()))
  double *RESTRICT const _data_u = u->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
  internal_advectivefluxkernel_double_precision_cuda_advectivefluxkernel_double_precision_cuda::advectivefluxkernel_double_precision_cuda_advectivefluxkernel_double_precision_cuda<<<_grid, _block, 0, stream>>>(_data_j, _data_rho, _data_u, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2, _stride_u_0, _stride_u_1, _stride_u_2, _stride_u_3);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
