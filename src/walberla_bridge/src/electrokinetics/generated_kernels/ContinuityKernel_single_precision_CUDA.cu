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
//! \\file ContinuityKernel_single_precision_CUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 007e77e077ad9d22b5eed6f3d3118240993e553c

#include <cmath>

#include "ContinuityKernel_single_precision_CUDA.h"
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

namespace internal_continuitykernel_single_precision_cuda_continuitykernel_single_precision_cuda {
static FUNC_PREFIX __launch_bounds__(256) void continuitykernel_single_precision_cuda_continuitykernel_single_precision_cuda(float *RESTRICT const _data_j, float *RESTRICT _data_rho, int64_t const _size_j_0, int64_t const _size_j_1, int64_t const _size_j_2, int64_t const _stride_j_0, int64_t const _stride_j_1, int64_t const _stride_j_2, int64_t const _stride_j_3, int64_t const _stride_rho_0, int64_t const _stride_rho_1, int64_t const _stride_rho_2) {
  if (blockDim.x * blockIdx.x + threadIdx.x + 1 < _size_j_0 - 1 && blockDim.y * blockIdx.y + threadIdx.y + 1 < _size_j_1 - 1 && blockDim.z * blockIdx.z + threadIdx.z + 1 < _size_j_2 - 1) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x + 1;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y + 1;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z + 1;
    _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] = _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_1 + _stride_j_2 * ctr_2 + _stride_j_2 + 9 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_1 + _stride_j_2 * ctr_2 - _stride_j_2 + 10 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_2 + 5 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 - _stride_j_2 + 6 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] + _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 - _stride_j_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 - _stride_j_1 + _stride_j_2 * ctr_2 + _stride_j_2 + 11 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_0 + _stride_j_1 * ctr_1 - _stride_j_1 + _stride_j_2 * ctr_2 - _stride_j_2 + 12 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_1 + _stride_j_2 * ctr_2 + _stride_j_2 + 7 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_1 + _stride_j_2 * ctr_2 + _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_1 + _stride_j_2 * ctr_2 - _stride_j_2 + 8 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] + _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_2 + 2 * _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3] - _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2];
  }
}
} // namespace internal_continuitykernel_single_precision_cuda_continuitykernel_single_precision_cuda

void ContinuityKernel_single_precision_CUDA::run(IBlock *block, gpuStream_t stream) {

  auto rho = block->getData<gpu::GPUField<float>>(rhoID);
  auto j = block->getData<gpu::GPUField<float>>(jID);

  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()))
  float *RESTRICT const _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho->nrOfGhostLayers()))
  float *RESTRICT _data_rho = rho->dataAt(-1, -1, -1, 0);
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
  dim3 _block(uint32_c(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)), uint32_c(((1024 < ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))) ? 1024 : ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))), uint32_c(((64 < ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))) ? 64 : ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))))));
  dim3 _grid(uint32_c(((_size_j_0 - 2) % (((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)) == 0 ? (int64_t)(_size_j_0 - 2) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)) : ((int64_t)(_size_j_0 - 2) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))) + 1)), uint32_c(((_size_j_1 - 2) % (((1024 < ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))) ? 1024 : ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))) == 0 ? (int64_t)(_size_j_1 - 2) / (int64_t)(((1024 < ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))) ? 1024 : ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))) : ((int64_t)(_size_j_1 - 2) / (int64_t)(((1024 < ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))) ? 1024 : ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) + 1)), uint32_c(((_size_j_2 - 2) % (((64 < ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))) ? 64 : ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))))) == 0 ? (int64_t)(_size_j_2 - 2) / (int64_t)(((64 < ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))) ? 64 : ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))))) : ((int64_t)(_size_j_2 - 2) / (int64_t)(((64 < ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))) ? 64 : ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))))) + 1)));
  internal_continuitykernel_single_precision_cuda_continuitykernel_single_precision_cuda::continuitykernel_single_precision_cuda_continuitykernel_single_precision_cuda<<<_grid, _block, 0, stream>>>(_data_j, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2);
}

void ContinuityKernel_single_precision_CUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block, gpuStream_t stream) {

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto rho = block->getData<gpu::GPUField<float>>(rhoID);
  auto j = block->getData<gpu::GPUField<float>>(jID);

  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()))
  float *RESTRICT const _data_j = j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(j->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho->nrOfGhostLayers()))
  float *RESTRICT _data_rho = rho->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
  dim3 _block(uint32_c(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)), uint32_c(((1024 < ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))) ? 1024 : ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))), uint32_c(((64 < ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))) ? 64 : ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))))));
  dim3 _grid(uint32_c(((_size_j_0 - 2) % (((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)) == 0 ? (int64_t)(_size_j_0 - 2) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)) : ((int64_t)(_size_j_0 - 2) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))) + 1)), uint32_c(((_size_j_1 - 2) % (((1024 < ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))) ? 1024 : ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))) == 0 ? (int64_t)(_size_j_1 - 2) / (int64_t)(((1024 < ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))) ? 1024 : ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))) : ((int64_t)(_size_j_1 - 2) / (int64_t)(((1024 < ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))) ? 1024 : ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) + 1)), uint32_c(((_size_j_2 - 2) % (((64 < ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))) ? 64 : ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))))) == 0 ? (int64_t)(_size_j_2 - 2) / (int64_t)(((64 < ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))) ? 64 : ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))))) : ((int64_t)(_size_j_2 - 2) / (int64_t)(((64 < ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))) ? 64 : ((_size_j_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2))))))) ? _size_j_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2) * ((_size_j_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))) ? _size_j_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_j_0 - 2) ? 128 : _size_j_0 - 2)))))))))) + 1)));
  internal_continuitykernel_single_precision_cuda_continuitykernel_single_precision_cuda::continuitykernel_single_precision_cuda_continuitykernel_single_precision_cuda<<<_grid, _block, 0, stream>>>(_data_j, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
