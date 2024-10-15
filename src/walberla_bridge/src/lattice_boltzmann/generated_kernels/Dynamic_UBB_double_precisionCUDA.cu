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
//! \\file Dynamic_UBB_double_precisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3, lbmpy_walberla/pystencils_walberla from waLBerla commit 04f4adbdfc0af983e2d9b72e244d775f37d77034

#include "Dynamic_UBB_double_precisionCUDA.h"
#include "core/DataTypes.h"
#include "core/Macros.h"
#include "gpu/ErrorChecking.h"

#define FUNC_PREFIX __global__

using namespace std;

namespace walberla {
namespace lbm {

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wconversion"
#endif

#ifdef __CUDACC__

#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diag_suppress 177
#else

#endif
#endif
// NOLINTBEGIN(readability-non-const-parameter*)
namespace internal_dynamic_ubb_double_precisioncuda_boundary_Dynamic_UBB_double_precisionCUDA {
static FUNC_PREFIX __launch_bounds__(256) void dynamic_ubb_double_precisioncuda_boundary_Dynamic_UBB_double_precisionCUDA(uint8_t *RESTRICT const _data_indexVector, double *RESTRICT _data_pdfs, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int32_t indexVectorSize) {

  const int32_t f_in_inv_dir_idx[] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 16, 15, 18, 17, 12, 11, 14, 13};
  const int32_t f_in_inv_offsets_x[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1};
  const int32_t f_in_inv_offsets_y[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0};
  const int32_t f_in_inv_offsets_z[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};

  const double weights[] = {0.33333333333333333, 0.055555555555555556, 0.055555555555555556, 0.055555555555555556, 0.055555555555555556, 0.055555555555555556, 0.055555555555555556, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778, 0.027777777777777778};

  const int32_t neighbour_offset_x[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1};
  const int32_t neighbour_offset_y[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0};
  const int32_t neighbour_offset_z[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};

  if (blockDim.x * blockIdx.x + threadIdx.x < indexVectorSize) {
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
    const double rho = vel0Term + vel1Term + vel2Term + _data_pdfs_10_20_30[_stride_pdfs_0 * x] + _data_pdfs_10_20_316[_stride_pdfs_0 * x] + _data_pdfs_10_20_317[_stride_pdfs_0 * x] + _data_pdfs_10_20_32[_stride_pdfs_0 * x] + _data_pdfs_10_20_33[_stride_pdfs_0 * x] + _data_pdfs_10_20_36[_stride_pdfs_0 * x] + _data_pdfs_10_20_39[_stride_pdfs_0 * x];
    double *RESTRICT _data_pdfs00178f3386915a72 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_1 * f_in_inv_offsets_y[dir] + _stride_pdfs_2 * z + _stride_pdfs_2 * f_in_inv_offsets_z[dir] + _stride_pdfs_3 * f_in_inv_dir_idx[dir];
    uint8_t *RESTRICT _data_indexVector_116 = _data_indexVector + 16;
    uint8_t *RESTRICT _data_indexVector_124 = _data_indexVector + 24;
    uint8_t *RESTRICT _data_indexVector_132 = _data_indexVector + 32;
    double *RESTRICT _data_pdfs_10_20c107984fcf09b503 = _data_pdfs + _stride_pdfs_1 * y + _stride_pdfs_2 * z + _stride_pdfs_3 * dir;
    _data_pdfs00178f3386915a72[_stride_pdfs_0 * x + _stride_pdfs_0 * f_in_inv_offsets_x[dir]] = -rho * (6.0 * ((double)(neighbour_offset_x[dir])) * *((double *)(&_data_indexVector_116[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 6.0 * ((double)(neighbour_offset_y[dir])) * *((double *)(&_data_indexVector_124[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x])) + 6.0 * ((double)(neighbour_offset_z[dir])) * *((double *)(&_data_indexVector_132[40 * blockDim.x * blockIdx.x + 40 * threadIdx.x]))) * weights[dir] + _data_pdfs_10_20c107984fcf09b503[_stride_pdfs_0 * x];
  }
}
} // namespace internal_dynamic_ubb_double_precisioncuda_boundary_Dynamic_UBB_double_precisionCUDA

// NOLINTEND(readability-non-const-parameter*)
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__

#endif

void Dynamic_UBB_double_precisionCUDA::run_impl(IBlock *block, IndexVectors::Type type, gpuStream_t stream) {
  auto *indexVectors = block->getData<IndexVectors>(indexVectorID);
  int32_t indexVectorSize = int32_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerGpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto pdfs = block->getData<gpu::GPUField<double>>(pdfsID);

  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  dim3 _block(uint32_c(((256 < indexVectorSize) ? 256 : indexVectorSize)), uint32_c(1), uint32_c(1));
  dim3 _grid(uint32_c(((indexVectorSize) % (((256 < indexVectorSize) ? 256 : indexVectorSize)) == 0 ? (int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize)) : ((int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize))) + 1)), uint32_c(1), uint32_c(1));
  internal_dynamic_ubb_double_precisioncuda_boundary_Dynamic_UBB_double_precisionCUDA::dynamic_ubb_double_precisioncuda_boundary_Dynamic_UBB_double_precisionCUDA<<<_grid, _block, 0, stream>>>(_data_indexVector, _data_pdfs, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, indexVectorSize);
}

void Dynamic_UBB_double_precisionCUDA::run(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::ALL, stream);
}

void Dynamic_UBB_double_precisionCUDA::inner(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::INNER, stream);
}

void Dynamic_UBB_double_precisionCUDA::outer(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::OUTER, stream);
}

} // namespace lbm
} // namespace walberla
