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
//! \\file UpdateVelFromPDFSinglePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7+4.gc7d65a7, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 0aab9c0af2335b1f6fec75deae06e514ccb233ab

#include <cmath>

#include "UpdateVelFromPDFSinglePrecisionCUDA.h"
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

namespace internal_updatevelfrompdfsingleprecisioncuda_updatevelfrompdfsingleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void updatevelfrompdfsingleprecisioncuda_updatevelfrompdfsingleprecisioncuda(float *RESTRICT const _data_force, float *RESTRICT const _data_pdfs, float *RESTRICT _data_velocity, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3) {
  if (blockDim.x * blockIdx.x + threadIdx.x + 1 < _size_force_0 - 1 && blockDim.y * blockIdx.y + threadIdx.y + 1 < _size_force_1 - 1 && blockDim.z * blockIdx.z + threadIdx.z + 1 < _size_force_2 - 1) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x + 1;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y + 1;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z + 1;
    const float vel0Term = _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    const float momdensity_0 = vel0Term - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3];
    const float vel1Term = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    const float momdensity_1 = vel1Term - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    const float vel2Term = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3];
    const float delta_rho = vel0Term + vel1Term + vel2Term + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2];
    const float momdensity_2 = vel2Term - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const float rho = delta_rho + 1.0f;
    const float u_0 = momdensity_0 * ((1.0f) / (rho)) + 0.5f * ((1.0f) / (rho)) * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float u_1 = momdensity_1 * ((1.0f) / (rho)) + 0.5f * ((1.0f) / (rho)) * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float u_2 = momdensity_2 * ((1.0f) / (rho)) + 0.5f * ((1.0f) / (rho)) * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    _data_velocity[_stride_velocity_0 * ctr_0 + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2] = u_0;
    _data_velocity[_stride_velocity_0 * ctr_0 + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2 + _stride_velocity_3] = u_1;
    _data_velocity[_stride_velocity_0 * ctr_0 + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3] = u_2;
  }
}
} // namespace internal_updatevelfrompdfsingleprecisioncuda_updatevelfrompdfsingleprecisioncuda

void UpdateVelFromPDFSinglePrecisionCUDA::run(IBlock *block, gpuStream_t stream) {

  auto force = block->getData<gpu::GPUField<float>>(forceID);
  auto pdfs = block->getData<gpu::GPUField<float>>(pdfsID);
  auto velocity = block->getData<gpu::GPUField<float>>(velocityID);

  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force->nrOfGhostLayers()))
  float *RESTRICT const _data_force = force->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()))
  float *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
  float *RESTRICT _data_velocity = velocity->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(force->xSize()) + 2))
  const int64_t _size_force_0 = int64_t(int64_c(force->xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(force->ySize()) + 2))
  const int64_t _size_force_1 = int64_t(int64_c(force->ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(force->zSize()) + 2))
  const int64_t _size_force_2 = int64_t(int64_c(force->zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  dim3 _block(uint32_c(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)), uint32_c(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))), uint32_c(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))));
  dim3 _grid(uint32_c(((_size_force_0 - 2) % (((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)) == 0 ? (int64_t)(_size_force_0 - 2) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)) : ((int64_t)(_size_force_0 - 2) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))) + 1)), uint32_c(((_size_force_1 - 2) % (((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))) == 0 ? (int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))) : ((int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) + 1)), uint32_c(((_size_force_2 - 2) % (((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))) == 0 ? (int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))) : ((int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))))) + 1)));
  internal_updatevelfrompdfsingleprecisioncuda_updatevelfrompdfsingleprecisioncuda::updatevelfrompdfsingleprecisioncuda_updatevelfrompdfsingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
}

void UpdateVelFromPDFSinglePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block, gpuStream_t stream) {

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto force = block->getData<gpu::GPUField<float>>(forceID);
  auto pdfs = block->getData<gpu::GPUField<float>>(pdfsID);
  auto velocity = block->getData<gpu::GPUField<float>>(velocityID);

  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force->nrOfGhostLayers()))
  float *RESTRICT const _data_force = force->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
  float *RESTRICT _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
  const int64_t _size_force_0 = int64_t(int64_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
  const int64_t _size_force_1 = int64_t(int64_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
  const int64_t _size_force_2 = int64_t(int64_c(ci.zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  dim3 _block(uint32_c(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)), uint32_c(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))), uint32_c(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))));
  dim3 _grid(uint32_c(((_size_force_0 - 2) % (((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)) == 0 ? (int64_t)(_size_force_0 - 2) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)) : ((int64_t)(_size_force_0 - 2) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))) + 1)), uint32_c(((_size_force_1 - 2) % (((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))) == 0 ? (int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))) : ((int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) + 1)), uint32_c(((_size_force_2 - 2) % (((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))) == 0 ? (int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))) : ((int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))))) + 1)));
  internal_updatevelfrompdfsingleprecisioncuda_updatevelfrompdfsingleprecisioncuda::updatevelfrompdfsingleprecisioncuda_updatevelfrompdfsingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
