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
//! \\file DynamicUBBDoublePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7+13.gdfd203a, lbmpy v1.3.7+5.gd7100a3, sympy v1.10, lbmpy_walberla/pystencils_walberla from waLBerla commit c69cb11d6a95d32b2280544d3d9abde1fe5fdbb5

#include "DynamicUBBDoublePrecisionCUDA.h"
#include "core/DataTypes.h"
#include "core/Macros.h"
#include "gpu/ErrorChecking.h"

#define FUNC_PREFIX __global__

using namespace std;

namespace walberla {
namespace lbm {

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
namespace internal_dynamicubbdoubleprecisioncuda_boundary_DynamicUBBDoublePrecisionCUDA {
static FUNC_PREFIX __launch_bounds__(256) void dynamicubbdoubleprecisioncuda_boundary_DynamicUBBDoublePrecisionCUDA(uint8_t *RESTRICT _data_forceVector, uint8_t *RESTRICT const _data_indexVector, double *RESTRICT _data_pdfs, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int32_t forceVectorSize) {

  const int32_t f_in_inv_dir_idx[] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 16, 15, 18, 17, 12, 11, 14, 13};
  const int32_t f_in_inv_offsets_x[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1};
  const int32_t f_in_inv_offsets_y[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0};
  const int32_t f_in_inv_offsets_z[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};

  const double weights[] = {((double)(0.33333333333333333)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778))};

  const int32_t neighbour_offset_x[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1};
  const int32_t neighbour_offset_y[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0};
  const int32_t neighbour_offset_z[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};

  if (blockDim.x * blockIdx.x + threadIdx.x < forceVectorSize) {
    uint8_t *RESTRICT _data_indexVector_10 = _data_indexVector;
    const int32_t x = *((int32_t *)(&_data_indexVector_10[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
    uint8_t *RESTRICT _data_indexVector_14 = _data_indexVector + 4;
    const int32_t y = *((int32_t *)(&_data_indexVector_14[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
    uint8_t *RESTRICT _data_indexVector_18 = _data_indexVector + 8;
    const int32_t z = *((int32_t *)(&_data_indexVector_18[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
    uint8_t *RESTRICT _data_indexVector_112 = _data_indexVector + 12;
    const int32_t dir = *((int32_t *)(&_data_indexVector_112[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]));
    double *RESTRICT _data_pdfs_10_20_310 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_314 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_318 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_34 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_38 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 8 * _stride_pdfs_3;
    const double vel0Term = _data_pdfs_10_20_310[_stride_pdfs_0 * x] + _data_pdfs_10_20_314[_stride_pdfs_0 * x] + _data_pdfs_10_20_318[_stride_pdfs_0 * x] + _data_pdfs_10_20_34[_stride_pdfs_0 * x] + _data_pdfs_10_20_38[_stride_pdfs_0 * x];
    double *RESTRICT _data_pdfs_10_20_31 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_311 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_315 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_37 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 7 * _stride_pdfs_3;
    const double vel1Term = _data_pdfs_10_20_311[_stride_pdfs_0 * x] + _data_pdfs_10_20_315[_stride_pdfs_0 * x] + _data_pdfs_10_20_31[_stride_pdfs_0 * x] + _data_pdfs_10_20_37[_stride_pdfs_0 * x];
    double *RESTRICT _data_pdfs_10_20_312 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_313 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_35 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 5 * _stride_pdfs_3;
    const double vel2Term = _data_pdfs_10_20_312[_stride_pdfs_0 * x] + _data_pdfs_10_20_313[_stride_pdfs_0 * x] + _data_pdfs_10_20_35[_stride_pdfs_0 * x];
    double *RESTRICT _data_pdfs_10_20_30 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z;
    double *RESTRICT _data_pdfs_10_20_316 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_317 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_32 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_33 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_36 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_10_20_39 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 9 * _stride_pdfs_3;
    const double delta_rho = vel0Term + vel1Term + vel2Term + _data_pdfs_10_20_30[_stride_pdfs_0 * x] + _data_pdfs_10_20_316[_stride_pdfs_0 * x] + _data_pdfs_10_20_317[_stride_pdfs_0 * x] + _data_pdfs_10_20_32[_stride_pdfs_0 * x] + _data_pdfs_10_20_33[_stride_pdfs_0 * x] + _data_pdfs_10_20_36[_stride_pdfs_0 * x] + _data_pdfs_10_20_39[_stride_pdfs_0 * x];
    const double rho = delta_rho + 1.0;
    double *RESTRICT _data_pdfs8a76ae0ea77ad8b8 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_1 * f_in_inv_offsets_y[dir] + _stride_pdfs_2 * z + _stride_pdfs_2 * f_in_inv_offsets_z[dir] + _stride_pdfs_3 * f_in_inv_dir_idx[dir];
    uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
    uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
    uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
    double *RESTRICT _data_pdfs_10_204b8170353ebeb30e = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + _stride_pdfs_3 * dir;
    _data_pdfs8a76ae0ea77ad8b8[_stride_pdfs_0 * x + _stride_pdfs_0 * f_in_inv_offsets_x[dir]] = -rho * (6.0 * ((double)(neighbour_offset_x[dir])) * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 6.0 * ((double)(neighbour_offset_y[dir])) * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 6.0 * ((double)(neighbour_offset_z[dir])) * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]))) * weights[dir] + _data_pdfs_10_204b8170353ebeb30e[_stride_pdfs_0 * x];
    const double f = -rho * (6.0 * ((double)(neighbour_offset_x[dir])) * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 6.0 * ((double)(neighbour_offset_y[dir])) * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 6.0 * ((double)(neighbour_offset_z[dir])) * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]))) * weights[dir] + 2.0 * _data_pdfs_10_204b8170353ebeb30e[_stride_pdfs_0 * x];
    uint8_t *RESTRICT _data_forceVector_10 = _data_forceVector;
    *((double *)(&_data_forceVector_10[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x])) = f * ((double)(neighbour_offset_x[dir]));
    uint8_t *RESTRICT _data_forceVector_18 = _data_forceVector + 8;
    *((double *)(&_data_forceVector_18[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x])) = f * ((double)(neighbour_offset_y[dir]));
    uint8_t *RESTRICT _data_forceVector_116 = _data_forceVector + 16;
    *((double *)(&_data_forceVector_116[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x])) = f * ((double)(neighbour_offset_z[dir]));
  }
}
} // namespace internal_dynamicubbdoubleprecisioncuda_boundary_DynamicUBBDoublePrecisionCUDA

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

void DynamicUBBDoublePrecisionCUDA::run_impl(IBlock *block, IndexVectors::Type type, gpuStream_t stream) {
  auto *indexVectors = block->getData<IndexVectors>(indexVectorID);
  int32_t indexVectorSize = int32_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerGpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto *forceVector = block->getData<ForceVector>(forceVectorID);
  WALBERLA_ASSERT_EQUAL(indexVectorSize, int32_c(forceVector->forceVector().size()))

  auto forcePointer = forceVector->pointerGpu();
  int32_t forceVectorSize = int32_c(forceVector->forceVector().size());

  uint8_t *_data_forceVector = reinterpret_cast<uint8_t *>(forcePointer);

  auto pdfs = block->getData<gpu::GPUField<double>>(pdfsID);

  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  dim3 _block(uint32_c(((256 < indexVectorSize) ? 256 : indexVectorSize)), uint32_c(1), uint32_c(1));
  dim3 _grid(uint32_c(((indexVectorSize) % (((256 < indexVectorSize) ? 256 : indexVectorSize)) == 0 ? (int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize)) : ((int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize))) + 1)), uint32_c(1), uint32_c(1));
  internal_dynamicubbdoubleprecisioncuda_boundary_DynamicUBBDoublePrecisionCUDA::dynamicubbdoubleprecisioncuda_boundary_DynamicUBBDoublePrecisionCUDA<<<_grid, _block, 0, stream>>>(_data_forceVector, _data_indexVector, _data_pdfs, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, forceVectorSize);
}

void DynamicUBBDoublePrecisionCUDA::run(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::ALL, stream);
}

void DynamicUBBDoublePrecisionCUDA::inner(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::INNER, stream);
}

void DynamicUBBDoublePrecisionCUDA::outer(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::OUTER, stream);
}

} // namespace lbm
} // namespace walberla
