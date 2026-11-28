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
//! \\file PackInfoVecSinglePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 007e77e077ad9d22b5eed6f3d3118240993e553c

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "domain_decomposition/IBlock.h"
#include "stencil/Directions.h"

#include "PackInfoVecSinglePrecisionCUDA.h"

#define FUNC_PREFIX __global__

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

namespace walberla {
namespace pystencils {

using walberla::cell::CellInterval;
using walberla::stencil::Direction;

namespace internal_pack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE {
static FUNC_PREFIX __launch_bounds__(256) void pack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE(float *RESTRICT _data_buffer, float *RESTRICT const _data_field, int64_t const _size_field_0, int64_t const _size_field_1, int64_t const _size_field_2, int64_t const _stride_field_0, int64_t const _stride_field_1, int64_t const _stride_field_2, int64_t const _stride_field_3) {
  if (blockDim.x * blockIdx.x + threadIdx.x < _size_field_0 && blockDim.y * blockIdx.y + threadIdx.y < _size_field_1 && blockDim.z * blockIdx.z + threadIdx.z < _size_field_2) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z;
    _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0] = _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2];
    _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0 + 1] = _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2 + _stride_field_3];
    _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0 + 2] = _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2 + 2 * _stride_field_3];
  }
}
} // namespace internal_pack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE

namespace internal_unpack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE {
static FUNC_PREFIX __launch_bounds__(256) void unpack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE(float *RESTRICT const _data_buffer, float *RESTRICT _data_field, int64_t const _size_field_0, int64_t const _size_field_1, int64_t const _size_field_2, int64_t const _stride_field_0, int64_t const _stride_field_1, int64_t const _stride_field_2, int64_t const _stride_field_3) {
  if (blockDim.x * blockIdx.x + threadIdx.x < _size_field_0 && blockDim.y * blockIdx.y + threadIdx.y < _size_field_1 && blockDim.z * blockIdx.z + threadIdx.z < _size_field_2) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z;
    _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2] = _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0];
    _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2 + _stride_field_3] = _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0 + 1];
    _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2 + 2 * _stride_field_3] = _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0 + 2];
  }
}
} // namespace internal_unpack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE

void PackInfoVecSinglePrecisionCUDA::pack(Direction dir, unsigned char *byte_buffer, IBlock *block, gpuStream_t stream) {
  float *buffer = reinterpret_cast<float *>(byte_buffer);

  auto field = block->getData<gpu::GPUField<float>>(fieldID);

  CellInterval ci;
  field->getSliceBeforeGhostLayer(dir, ci, 1, false);

  switch (dir) {
  case stencil::SW:
  case stencil::BW:
  case stencil::W:
  case stencil::TW:
  case stencil::NW:
  case stencil::BS:
  case stencil::S:
  case stencil::TS:
  case stencil::B:
  case stencil::C:
  case stencil::T:
  case stencil::BN:
  case stencil::N:
  case stencil::TN:
  case stencil::SE:
  case stencil::BE:
  case stencil::E:
  case stencil::TE:
  case stencil::NE: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(field->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(field->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(field->nrOfGhostLayers()))
    float *RESTRICT const _data_field = field->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_field_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_field_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_field_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_field_0 = int64_t(field->xStride());
    const int64_t _stride_field_1 = int64_t(field->yStride());
    const int64_t _stride_field_2 = int64_t(field->zStride());
    const int64_t _stride_field_3 = int64_t(1 * int64_t(field->fStride()));
    dim3 _block(uint32_c(((128 < _size_field_0) ? 128 : _size_field_0)), uint32_c(((1024 < ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))) ? 1024 : ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))), uint32_c(((64 < ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))) ? 64 : ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))))));
    dim3 _grid(uint32_c(((_size_field_0) % (((128 < _size_field_0) ? 128 : _size_field_0)) == 0 ? (int64_t)(_size_field_0) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)) : ((int64_t)(_size_field_0) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))) + 1)), uint32_c(((_size_field_1) % (((1024 < ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))) ? 1024 : ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))) == 0 ? (int64_t)(_size_field_1) / (int64_t)(((1024 < ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))) ? 1024 : ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))) : ((int64_t)(_size_field_1) / (int64_t)(((1024 < ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))) ? 1024 : ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) + 1)), uint32_c(((_size_field_2) % (((64 < ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))) ? 64 : ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))))) == 0 ? (int64_t)(_size_field_2) / (int64_t)(((64 < ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))) ? 64 : ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))))) : ((int64_t)(_size_field_2) / (int64_t)(((64 < ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))) ? 64 : ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))))) + 1)));
    internal_pack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE::pack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_field, _size_field_0, _size_field_1, _size_field_2, _stride_field_0, _stride_field_1, _stride_field_2, _stride_field_3);
    break;
  }

  default:
    return;
  }
}

void PackInfoVecSinglePrecisionCUDA::unpack(Direction dir, unsigned char *byte_buffer, IBlock *block, gpuStream_t stream) {
  float *buffer = reinterpret_cast<float *>(byte_buffer);

  auto field = block->getData<gpu::GPUField<float>>(fieldID);

  CellInterval ci;
  field->getGhostRegion(dir, ci, 1, false);
  auto communciationDirection = stencil::inverseDir[dir];

  switch (communciationDirection) {
  case stencil::SW:
  case stencil::BW:
  case stencil::W:
  case stencil::TW:
  case stencil::NW:
  case stencil::BS:
  case stencil::S:
  case stencil::TS:
  case stencil::B:
  case stencil::C:
  case stencil::T:
  case stencil::BN:
  case stencil::N:
  case stencil::TN:
  case stencil::SE:
  case stencil::BE:
  case stencil::E:
  case stencil::TE:
  case stencil::NE: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(field->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(field->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(field->nrOfGhostLayers()))
    float *RESTRICT _data_field = field->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_field_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_field_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_field_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_field_0 = int64_t(field->xStride());
    const int64_t _stride_field_1 = int64_t(field->yStride());
    const int64_t _stride_field_2 = int64_t(field->zStride());
    const int64_t _stride_field_3 = int64_t(1 * int64_t(field->fStride()));
    dim3 _block(uint32_c(((128 < _size_field_0) ? 128 : _size_field_0)), uint32_c(((1024 < ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))) ? 1024 : ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))), uint32_c(((64 < ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))) ? 64 : ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))))));
    dim3 _grid(uint32_c(((_size_field_0) % (((128 < _size_field_0) ? 128 : _size_field_0)) == 0 ? (int64_t)(_size_field_0) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)) : ((int64_t)(_size_field_0) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))) + 1)), uint32_c(((_size_field_1) % (((1024 < ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))) ? 1024 : ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))) == 0 ? (int64_t)(_size_field_1) / (int64_t)(((1024 < ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))) ? 1024 : ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))) : ((int64_t)(_size_field_1) / (int64_t)(((1024 < ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))) ? 1024 : ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) + 1)), uint32_c(((_size_field_2) % (((64 < ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))) ? 64 : ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))))) == 0 ? (int64_t)(_size_field_2) / (int64_t)(((64 < ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))) ? 64 : ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))))) : ((int64_t)(_size_field_2) / (int64_t)(((64 < ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))) ? 64 : ((_size_field_2 < ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0))))))) ? _size_field_2 : ((int64_t)(256) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0) * ((_size_field_1 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))) ? _size_field_1 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_field_0) ? 128 : _size_field_0)))))))))) + 1)));
    internal_unpack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE::unpack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_field, _size_field_0, _size_field_1, _size_field_2, _stride_field_0, _stride_field_1, _stride_field_2, _stride_field_3);
    break;
  }

  default:
    return;
  }
}

uint_t PackInfoVecSinglePrecisionCUDA::size(stencil::Direction dir, IBlock *block) {
  auto field = block->getData<gpu::GPUField<float>>(fieldID);

  CellInterval ci;
  field->getGhostRegion(dir, ci, 1, false);

  uint_t elementsPerCell = 0;

  switch (dir) {
  case stencil::SW:
  case stencil::BW:
  case stencil::W:
  case stencil::TW:
  case stencil::NW:
  case stencil::BS:
  case stencil::S:
  case stencil::TS:
  case stencil::B:
  case stencil::C:
  case stencil::T:
  case stencil::BN:
  case stencil::N:
  case stencil::TN:
  case stencil::SE:
  case stencil::BE:
  case stencil::E:
  case stencil::TE:
  case stencil::NE:
    elementsPerCell = 3;
    break;

  default:
    elementsPerCell = 0;
  }
  return ci.numCells() * elementsPerCell * sizeof(float);
}

} // namespace pystencils
} // namespace walberla
