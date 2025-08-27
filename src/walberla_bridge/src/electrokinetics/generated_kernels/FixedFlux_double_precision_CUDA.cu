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
//! \\file FixedFlux_double_precision_CUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7+13.gdfd203a, lbmpy v1.3.7+10.gd3f6236, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit c69cb11d6a95d32b2280544d3d9abde1fe5fdbb5

#include "FixedFlux_double_precision_CUDA.h"
#include "core/DataTypes.h"
#include "core/Macros.h"
#include "gpu/ErrorChecking.h"

#define FUNC_PREFIX __global__

using namespace std;

namespace walberla {
namespace pystencils {

#if defined(__NVCC__)
#define RESTRICT __restrict__
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // unused variable
#else
#pragma push
#pragma diag_suppress 177 // unused variable
#endif                    // defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstrict-aliasing"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wsign-compare"
#else
// clang compiling CUDA code in host mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstrict-aliasing"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wsign-compare"
#endif // defined(__CUDA_ARCH__)
#endif // defined(__CUDA__)
#elif defined(__GNUC__) or defined(__GNUG__)
#define RESTRICT __restrict__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wconversion"
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

// NOLINTBEGIN(readability-non-const-parameter*)
namespace internal_fixedflux_double_precision_cuda_boundary_FixedFlux_double_precision_CUDA {
static FUNC_PREFIX __launch_bounds__(256) void fixedflux_double_precision_cuda_boundary_FixedFlux_double_precision_CUDA(double *RESTRICT const _data_flux, uint8_t *RESTRICT const _data_indexVector, int64_t const _stride_flux_0, int64_t const _stride_flux_1, int64_t const _stride_flux_2, int64_t const _stride_flux_3, int32_t indexVectorSize) {
  if (blockDim.x * blockIdx.x + threadIdx.x < indexVectorSize) {
    uint8_t *RESTRICT _data_indexVector_10 = _data_indexVector;
    const int32_t x = *((int32_t *)(&_data_indexVector_10[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
    uint8_t *RESTRICT _data_indexVector_14 = _data_indexVector + 4;
    const int32_t y = *((int32_t *)(&_data_indexVector_14[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
    uint8_t *RESTRICT _data_indexVector_18 = _data_indexVector + 8;
    const int32_t z = *((int32_t *)(&_data_indexVector_18[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));

    const int32_t cx[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1};
    const int32_t cy[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1};
    const int32_t cz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1};
    const int32_t invdir[] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 16, 15, 18, 17, 12, 11, 14, 13, 26, 25, 24, 23, 22, 21, 20, 19};

    uint8_t *RESTRICT _data_indexVector_112 = _data_indexVector + 12;
    const int32_t dir = *((int32_t *)(&_data_indexVector_112[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
    if (((dir) == (26))) {
      double *RESTRICT _data_flux_10_20_39 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 9 * _stride_flux_3;
      uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
      uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
      uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
      _data_flux_10_20_39[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
    } else {
      if (((dir) == (25))) {
        double *RESTRICT _data_flux_1m1_2m1_312 = _data_flux + _stride_flux_1 * y - _stride_flux_1 + _stride_flux_2 * z - _stride_flux_2 + 12 * _stride_flux_3;
        uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
        uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
        uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
        _data_flux_1m1_2m1_312[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
      } else {
        if (((dir) == (24))) {
          double *RESTRICT _data_flux_10_20_311 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 11 * _stride_flux_3;
          uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
          uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
          uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
          _data_flux_10_20_311[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
        } else {
          if (((dir) == (23))) {
            double *RESTRICT _data_flux_11_2m1_310 = _data_flux + _stride_flux_1 * y + _stride_flux_1 + _stride_flux_2 * z - _stride_flux_2 + 10 * _stride_flux_3;
            uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
            uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
            uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
            _data_flux_11_2m1_310[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
          } else {
            if (((dir) == (22))) {
              double *RESTRICT _data_flux_10_20_310 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 10 * _stride_flux_3;
              uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
              uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
              uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
              _data_flux_10_20_310[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
            } else {
              if (((dir) == (21))) {
                double *RESTRICT _data_flux_1m1_21_311 = _data_flux + _stride_flux_1 * y - _stride_flux_1 + _stride_flux_2 * z + _stride_flux_2 + 11 * _stride_flux_3;
                uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                _data_flux_1m1_21_311[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
              } else {
                if (((dir) == (20))) {
                  double *RESTRICT _data_flux_10_20_312 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 12 * _stride_flux_3;
                  uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                  uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                  uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                  _data_flux_10_20_312[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                } else {
                  if (((dir) == (19))) {
                    double *RESTRICT _data_flux_11_21_39 = _data_flux + _stride_flux_1 * y + _stride_flux_1 + _stride_flux_2 * z + _stride_flux_2 + 9 * _stride_flux_3;
                    uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                    uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                    uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                    _data_flux_11_21_39[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                  } else {
                    if (((dir) == (18))) {
                      double *RESTRICT _data_flux_10_2m1_36 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z - _stride_flux_2 + 6 * _stride_flux_3;
                      uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                      uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                      _data_flux_10_2m1_36[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                    } else {
                      if (((dir) == (17))) {
                        double *RESTRICT _data_flux_10_20_35 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 5 * _stride_flux_3;
                        uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                        uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                        _data_flux_10_20_35[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                      } else {
                        if (((dir) == (16))) {
                          double *RESTRICT _data_flux_10_20_37 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 7 * _stride_flux_3;
                          uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                          uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                          _data_flux_10_20_37[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                        } else {
                          if (((dir) == (15))) {
                            double *RESTRICT _data_flux_11_2m1_38 = _data_flux + _stride_flux_1 * y + _stride_flux_1 + _stride_flux_2 * z - _stride_flux_2 + 8 * _stride_flux_3;
                            uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                            uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                            _data_flux_11_2m1_38[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                          } else {
                            if (((dir) == (14))) {
                              double *RESTRICT _data_flux_10_21_35 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + _stride_flux_2 + 5 * _stride_flux_3;
                              uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                              uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                              _data_flux_10_21_35[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                            } else {
                              if (((dir) == (13))) {
                                double *RESTRICT _data_flux_10_20_36 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 6 * _stride_flux_3;
                                uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                                uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                                _data_flux_10_20_36[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                              } else {
                                if (((dir) == (12))) {
                                  double *RESTRICT _data_flux_10_20_38 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 8 * _stride_flux_3;
                                  uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                                  uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                                  _data_flux_10_20_38[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                } else {
                                  if (((dir) == (11))) {
                                    double *RESTRICT _data_flux_11_21_37 = _data_flux + _stride_flux_1 * y + _stride_flux_1 + _stride_flux_2 * z + _stride_flux_2 + 7 * _stride_flux_3;
                                    uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                                    uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                                    _data_flux_11_21_37[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                  } else {
                                    if (((dir) == (10))) {
                                      double *RESTRICT _data_flux_1m1_20_34 = _data_flux + _stride_flux_1 * y - _stride_flux_1 + _stride_flux_2 * z + 4 * _stride_flux_3;
                                      uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                                      uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                                      _data_flux_1m1_20_34[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                    } else {
                                      if (((dir) == (9))) {
                                        double *RESTRICT _data_flux_10_20_33 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 3 * _stride_flux_3;
                                        uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                                        uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                                        _data_flux_10_20_33[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                      } else {
                                        if (((dir) == (8))) {
                                          double *RESTRICT _data_flux_11_20_33 = _data_flux + _stride_flux_1 * y + _stride_flux_1 + _stride_flux_2 * z + 3 * _stride_flux_3;
                                          uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                                          uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                                          _data_flux_11_20_33[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) - 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                        } else {
                                          if (((dir) == (7))) {
                                            double *RESTRICT _data_flux_10_20_34 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 4 * _stride_flux_3;
                                            uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                                            uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                                            _data_flux_10_20_34[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                          } else {
                                            if (((dir) == (6))) {
                                              double *RESTRICT _data_flux_10_20_32 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + 2 * _stride_flux_3;
                                              uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                                              _data_flux_10_20_32[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                            } else {
                                              if (((dir) == (5))) {
                                                double *RESTRICT _data_flux_10_21_32 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + _stride_flux_2 + 2 * _stride_flux_3;
                                                uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
                                                _data_flux_10_21_32[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                              } else {
                                                if (((dir) == (4))) {
                                                  double *RESTRICT _data_flux_10_20_30 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z;
                                                  uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                                                  _data_flux_10_20_30[_stride_flux_0 * x + _stride_flux_0] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                                } else {
                                                  if (((dir) == (3))) {
                                                    double *RESTRICT _data_flux_10_20_30 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z;
                                                    uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
                                                    _data_flux_10_20_30[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                                  } else {
                                                    if (((dir) == (2))) {
                                                      double *RESTRICT _data_flux_10_20_31 = _data_flux + _stride_flux_1 * y + _stride_flux_2 * z + _stride_flux_3;
                                                      uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                                                      _data_flux_10_20_31[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                                    } else {
                                                      if (((dir) == (1))) {
                                                        double *RESTRICT _data_flux_11_20_31 = _data_flux + _stride_flux_1 * y + _stride_flux_1 + _stride_flux_2 * z + _stride_flux_3;
                                                        uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
                                                        _data_flux_11_20_31[_stride_flux_0 * x] = -0.1111111111111111 * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
} // namespace internal_fixedflux_double_precision_cuda_boundary_FixedFlux_double_precision_CUDA

// NOLINTEND(readability-non-const-parameter*)

#if defined(__NVCC__)
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic pop
#else
#pragma pop
#endif // defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
#pragma clang diagnostic pop
#else
// clang compiling CUDA code in host mode
#pragma clang diagnostic pop
#endif // defined(__CUDA_ARCH__)
#endif // defined(__CUDA__)
#elif defined(__GNUC__) or defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

void FixedFlux_double_precision_CUDA::run_impl(IBlock *block, IndexVectors::Type type, gpuStream_t stream) {
  auto *indexVectors = block->getData<IndexVectors>(indexVectorID);
  int32_t indexVectorSize = int32_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerGpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto flux = block->getData<gpu::GPUField<double>>(fluxID);

  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(flux->nrOfGhostLayers()))
  double *RESTRICT const _data_flux = flux->dataAt(0, 0, 0, 0);
  const int64_t _stride_flux_0 = int64_t(flux->xStride());
  const int64_t _stride_flux_1 = int64_t(flux->yStride());
  const int64_t _stride_flux_2 = int64_t(flux->zStride());
  const int64_t _stride_flux_3 = int64_t(1 * int64_t(flux->fStride()));
  dim3 _block(uint32_c(((256 < indexVectorSize) ? 256 : indexVectorSize)), uint32_c(1), uint32_c(1));
  dim3 _grid(uint32_c(((indexVectorSize) % (((256 < indexVectorSize) ? 256 : indexVectorSize)) == 0 ? (int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize)) : ((int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize))) + 1)), uint32_c(1), uint32_c(1));
  internal_fixedflux_double_precision_cuda_boundary_FixedFlux_double_precision_CUDA::fixedflux_double_precision_cuda_boundary_FixedFlux_double_precision_CUDA<<<_grid, _block, 0, stream>>>(_data_flux, _data_indexVector, _stride_flux_0, _stride_flux_1, _stride_flux_2, _stride_flux_3, indexVectorSize);
}

void FixedFlux_double_precision_CUDA::run(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::ALL, stream);
}

void FixedFlux_double_precision_CUDA::inner(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::INNER, stream);
}

void FixedFlux_double_precision_CUDA::outer(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::OUTER, stream);
}

} // namespace pystencils
} // namespace walberla
