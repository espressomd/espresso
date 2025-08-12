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
//! \\file ReactionKernelIndexed_1_single_precision_CUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7+13.gdfd203a, lbmpy v1.3.7+10.gd3f6236, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit c69cb11d6a95d32b2280544d3d9abde1fe5fdbb5

#include "ReactionKernelIndexed_1_single_precision_CUDA.h"
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
#pragma clang diagnostic ignored "-Wunused-variable"
#else
// clang compiling CUDA code in host mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#endif // defined(__CUDA_ARCH__)
#endif // defined(__CUDA__)
#elif defined(__GNUC__) or defined(__GNUG__)
#define RESTRICT __restrict__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

// NOLINTBEGIN(readability-non-const-parameter*)
namespace internal_reactionkernelindexed_1_single_precision_cuda_boundary_ReactionKernelIndexed_1_single_precision_CUDA {
static FUNC_PREFIX __launch_bounds__(256) void reactionkernelindexed_1_single_precision_cuda_boundary_ReactionKernelIndexed_1_single_precision_CUDA(uint8_t *RESTRICT const _data_indexVector, float *RESTRICT _data_rho_0, int64_t const _stride_rho_0_0, int64_t const _stride_rho_0_1, int64_t const _stride_rho_0_2, int32_t indexVectorSize, float order_0, float rate_coefficient, float stoech_0) {
  if (blockDim.x * blockIdx.x + threadIdx.x < indexVectorSize) {
    uint8_t *RESTRICT _data_indexVector_10 = _data_indexVector;
    const int32_t x = *((int32_t *)(&_data_indexVector_10[12 * blockDim.x * blockIdx.x + 12 * threadIdx.x]));
    uint8_t *RESTRICT _data_indexVector_14 = _data_indexVector + 4;
    const int32_t y = *((int32_t *)(&_data_indexVector_14[12 * blockDim.x * blockIdx.x + 12 * threadIdx.x]));
    uint8_t *RESTRICT _data_indexVector_18 = _data_indexVector + 8;
    const int32_t z = *((int32_t *)(&_data_indexVector_18[12 * blockDim.x * blockIdx.x + 12 * threadIdx.x]));

    const int32_t cx[] = {0};
    const int32_t cy[] = {0};
    const int32_t cz[] = {0};
    const int32_t invdir[] = {0};

    float *RESTRICT _data_rho_0_10_20 = _data_rho_0 + _stride_rho_0_1 * y + _stride_rho_0_2 * z;
    const float local_rho_0 = _data_rho_0_10_20[_stride_rho_0_0 * x];
    const float rate_factor = rate_coefficient * powf(local_rho_0, order_0);
    _data_rho_0_10_20[_stride_rho_0_0 * x] = local_rho_0 + rate_factor * stoech_0;
  }
}
} // namespace internal_reactionkernelindexed_1_single_precision_cuda_boundary_ReactionKernelIndexed_1_single_precision_CUDA

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

void ReactionKernelIndexed_1_single_precision_CUDA::run_impl(IBlock *block, IndexVectors::Type type, gpuStream_t stream) {
  auto *indexVectors = block->uncheckedFastGetData<IndexVectors>(indexVectorID);
  int32_t indexVectorSize = int32_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerGpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto rho_0 = block->getData<gpu::GPUField<float>>(rho_0ID);

  auto &order_0 = order_0_;
  auto &stoech_0 = stoech_0_;
  auto &rate_coefficient = rate_coefficient_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_0->nrOfGhostLayers()))
  float *RESTRICT _data_rho_0 = rho_0->dataAt(0, 0, 0, 0);
  const int64_t _stride_rho_0_0 = int64_t(rho_0->xStride());
  const int64_t _stride_rho_0_1 = int64_t(rho_0->yStride());
  const int64_t _stride_rho_0_2 = int64_t(rho_0->zStride());
  dim3 _block(uint32_c(((256 < indexVectorSize) ? 256 : indexVectorSize)), uint32_c(1), uint32_c(1));
  dim3 _grid(uint32_c(((indexVectorSize) % (((256 < indexVectorSize) ? 256 : indexVectorSize)) == 0 ? (int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize)) : ((int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize))) + 1)), uint32_c(1), uint32_c(1));
  internal_reactionkernelindexed_1_single_precision_cuda_boundary_ReactionKernelIndexed_1_single_precision_CUDA::reactionkernelindexed_1_single_precision_cuda_boundary_ReactionKernelIndexed_1_single_precision_CUDA<<<_grid, _block, 0, stream>>>(_data_indexVector, _data_rho_0, _stride_rho_0_0, _stride_rho_0_1, _stride_rho_0_2, indexVectorSize, order_0, rate_coefficient, stoech_0);
}

void ReactionKernelIndexed_1_single_precision_CUDA::run(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::ALL, stream);
}

void ReactionKernelIndexed_1_single_precision_CUDA::inner(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::INNER, stream);
}

void ReactionKernelIndexed_1_single_precision_CUDA::outer(IBlock *block, gpuStream_t stream) {
  run_impl(block, IndexVectors::OUTER, stream);
}

} // namespace pystencils
} // namespace walberla
