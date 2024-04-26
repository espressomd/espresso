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
//! \\file StreamSweepDoublePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.2, lbmpy_walberla/pystencils_walberla from waLBerla commit 0c8b4b926c6979288fd8a6846d02ec0870e1fe41

#include <cmath>

#include "StreamSweepDoublePrecisionCUDA.h"
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

namespace internal_streamsweepdoubleprecisioncuda_streamsweepdoubleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void streamsweepdoubleprecisioncuda_streamsweepdoubleprecisioncuda(double *RESTRICT const _data_force, double *RESTRICT const _data_pdfs, double *RESTRICT _data_pdfs_tmp, double *RESTRICT _data_velocity, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3) {
  if (blockDim.x * blockIdx.x + threadIdx.x + 1 < _size_force_0 - 1 && blockDim.y * blockIdx.y + threadIdx.y + 1 < _size_force_1 - 1 && blockDim.z * blockIdx.z + threadIdx.z + 1 < _size_force_2 - 1) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x + 1;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y + 1;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z + 1;
    double *RESTRICT _data_pdfs_10_20_30 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2;
    const double streamed_0 = _data_pdfs_10_20_30[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_1m1_20_31 = _data_pdfs + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    const double streamed_1 = _data_pdfs_1m1_20_31[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_11_20_32 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    const double streamed_2 = _data_pdfs_11_20_32[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_10_20_33 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    const double streamed_3 = _data_pdfs_10_20_33[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
    double *RESTRICT _data_pdfs_10_20_34 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    const double streamed_4 = _data_pdfs_10_20_34[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
    double *RESTRICT _data_pdfs_10_2m1_35 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3;
    const double streamed_5 = _data_pdfs_10_2m1_35[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_10_21_36 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3;
    const double streamed_6 = _data_pdfs_10_21_36[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_1m1_20_37 = _data_pdfs + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    const double streamed_7 = _data_pdfs_1m1_20_37[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
    double *RESTRICT _data_pdfs_1m1_20_38 = _data_pdfs + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    const double streamed_8 = _data_pdfs_1m1_20_38[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
    double *RESTRICT _data_pdfs_11_20_39 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    const double streamed_9 = _data_pdfs_11_20_39[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
    double *RESTRICT _data_pdfs_11_20_310 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    const double streamed_10 = _data_pdfs_11_20_310[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
    double *RESTRICT _data_pdfs_1m1_2m1_311 = _data_pdfs + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3;
    const double streamed_11 = _data_pdfs_1m1_2m1_311[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_11_2m1_312 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3;
    const double streamed_12 = _data_pdfs_11_2m1_312[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_10_2m1_313 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3;
    const double streamed_13 = _data_pdfs_10_2m1_313[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
    double *RESTRICT _data_pdfs_10_2m1_314 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3;
    const double streamed_14 = _data_pdfs_10_2m1_314[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
    double *RESTRICT _data_pdfs_1m1_21_315 = _data_pdfs + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3;
    const double streamed_15 = _data_pdfs_1m1_21_315[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_11_21_316 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3;
    const double streamed_16 = _data_pdfs_11_21_316[_stride_pdfs_0 * ctr_0];
    double *RESTRICT _data_pdfs_10_21_317 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3;
    const double streamed_17 = _data_pdfs_10_21_317[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
    double *RESTRICT _data_pdfs_10_21_318 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3;
    const double streamed_18 = _data_pdfs_10_21_318[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
    const double vel0Term = streamed_10 + streamed_14 + streamed_18 + streamed_4 + streamed_8;
    const double momdensity_0 = streamed_13 * -1.0 + streamed_17 * -1.0 + streamed_3 * -1.0 + streamed_7 * -1.0 + streamed_9 * -1.0 + vel0Term;
    const double vel1Term = streamed_1 + streamed_11 + streamed_15 + streamed_7;
    const double momdensity_1 = streamed_10 * -1.0 + streamed_12 * -1.0 + streamed_16 * -1.0 + streamed_2 * -1.0 + streamed_8 + streamed_9 * -1.0 + vel1Term;
    const double vel2Term = streamed_12 + streamed_13 + streamed_5;
    const double rho = streamed_0 + streamed_16 + streamed_17 + streamed_2 + streamed_3 + streamed_6 + streamed_9 + vel0Term + vel1Term + vel2Term;
    const double momdensity_2 = streamed_11 + streamed_14 + streamed_15 * -1.0 + streamed_16 * -1.0 + streamed_17 * -1.0 + streamed_18 * -1.0 + streamed_6 * -1.0 + vel2Term;
    double *RESTRICT _data_force_10_20_30 = _data_force + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2;
    const double u_0 = momdensity_0 * ((1.0) / (rho)) + 0.5 * ((1.0) / (rho)) * _data_force_10_20_30[_stride_force_0 * ctr_0];
    double *RESTRICT _data_force_10_20_31 = _data_force + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3;
    const double u_1 = momdensity_1 * ((1.0) / (rho)) + 0.5 * ((1.0) / (rho)) * _data_force_10_20_31[_stride_force_0 * ctr_0];
    double *RESTRICT _data_force_10_20_32 = _data_force + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    const double u_2 = momdensity_2 * ((1.0) / (rho)) + 0.5 * ((1.0) / (rho)) * _data_force_10_20_32[_stride_force_0 * ctr_0];
    double *RESTRICT _data_velocity_10_20_30 = _data_velocity + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2;
    _data_velocity_10_20_30[_stride_velocity_0 * ctr_0] = u_0;
    double *RESTRICT _data_velocity_10_20_31 = _data_velocity + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    _data_velocity_10_20_31[_stride_velocity_0 * ctr_0] = u_1;
    double *RESTRICT _data_velocity_10_20_32 = _data_velocity + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    _data_velocity_10_20_32[_stride_velocity_0 * ctr_0] = u_2;
    double *RESTRICT _data_pdfs_tmp_10_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2;
    _data_pdfs_tmp_10_20_30[_stride_pdfs_tmp_0 * ctr_0] = streamed_0;
    double *RESTRICT _data_pdfs_tmp_10_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_31[_stride_pdfs_tmp_0 * ctr_0] = streamed_1;
    double *RESTRICT _data_pdfs_tmp_10_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_32[_stride_pdfs_tmp_0 * ctr_0] = streamed_2;
    double *RESTRICT _data_pdfs_tmp_10_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_33[_stride_pdfs_tmp_0 * ctr_0] = streamed_3;
    double *RESTRICT _data_pdfs_tmp_10_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_34[_stride_pdfs_tmp_0 * ctr_0] = streamed_4;
    double *RESTRICT _data_pdfs_tmp_10_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_35[_stride_pdfs_tmp_0 * ctr_0] = streamed_5;
    double *RESTRICT _data_pdfs_tmp_10_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_36[_stride_pdfs_tmp_0 * ctr_0] = streamed_6;
    double *RESTRICT _data_pdfs_tmp_10_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_37[_stride_pdfs_tmp_0 * ctr_0] = streamed_7;
    double *RESTRICT _data_pdfs_tmp_10_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_38[_stride_pdfs_tmp_0 * ctr_0] = streamed_8;
    double *RESTRICT _data_pdfs_tmp_10_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_39[_stride_pdfs_tmp_0 * ctr_0] = streamed_9;
    double *RESTRICT _data_pdfs_tmp_10_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_310[_stride_pdfs_tmp_0 * ctr_0] = streamed_10;
    double *RESTRICT _data_pdfs_tmp_10_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_311[_stride_pdfs_tmp_0 * ctr_0] = streamed_11;
    double *RESTRICT _data_pdfs_tmp_10_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_312[_stride_pdfs_tmp_0 * ctr_0] = streamed_12;
    double *RESTRICT _data_pdfs_tmp_10_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_313[_stride_pdfs_tmp_0 * ctr_0] = streamed_13;
    double *RESTRICT _data_pdfs_tmp_10_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_314[_stride_pdfs_tmp_0 * ctr_0] = streamed_14;
    double *RESTRICT _data_pdfs_tmp_10_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_315[_stride_pdfs_tmp_0 * ctr_0] = streamed_15;
    double *RESTRICT _data_pdfs_tmp_10_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_316[_stride_pdfs_tmp_0 * ctr_0] = streamed_16;
    double *RESTRICT _data_pdfs_tmp_10_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_317[_stride_pdfs_tmp_0 * ctr_0] = streamed_17;
    double *RESTRICT _data_pdfs_tmp_10_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3;
    _data_pdfs_tmp_10_20_318[_stride_pdfs_tmp_0 * ctr_0] = streamed_18;
  }
}
} // namespace internal_streamsweepdoubleprecisioncuda_streamsweepdoubleprecisioncuda

void StreamSweepDoublePrecisionCUDA::run(IBlock *block, gpuStream_t stream) {
  auto force = block->getData<gpu::GPUField<double>>(forceID);
  auto velocity = block->getData<gpu::GPUField<double>>(velocityID);
  auto pdfs = block->getData<gpu::GPUField<double>>(pdfsID);
  gpu::GPUField<double> *pdfs_tmp;
  {
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find(pdfs);
    if (it != cache_pdfs_.end()) {
      pdfs_tmp = *it;
    } else {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_.insert(pdfs_tmp);
    }
  }

  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force->nrOfGhostLayers()))
  double *RESTRICT const _data_force = force->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
  double *RESTRICT _data_velocity = velocity->dataAt(-1, -1, -1, 0);
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
  const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  dim3 _block(uint32_t(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)), uint32_t(((1024 < ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))), uint32_t(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))))));
  dim3 _grid(uint32_t(((_size_force_0 - 2) % (((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)) == 0 ? (int64_t)(_size_force_0 - 2) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)) : ((int64_t)(_size_force_0 - 2) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))) + 1)), uint32_t(((_size_force_1 - 2) % (((1024 < ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))) == 0 ? (int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))) : ((int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) + 1)), uint32_t(((_size_force_2 - 2) % (((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))))) == 0 ? (int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))))) : ((int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))))) + 1)));
  internal_streamsweepdoubleprecisioncuda_streamsweepdoubleprecisioncuda::streamsweepdoubleprecisioncuda_streamsweepdoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_pdfs_tmp, _data_velocity, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
  pdfs->swapDataPointers(pdfs_tmp);
}

void StreamSweepDoublePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block, gpuStream_t stream) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto force = block->getData<gpu::GPUField<double>>(forceID);
  auto velocity = block->getData<gpu::GPUField<double>>(velocityID);
  auto pdfs = block->getData<gpu::GPUField<double>>(pdfsID);
  gpu::GPUField<double> *pdfs_tmp;
  {
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find(pdfs);
    if (it != cache_pdfs_.end()) {
      pdfs_tmp = *it;
    } else {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_.insert(pdfs_tmp);
    }
  }

  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force->nrOfGhostLayers()))
  double *RESTRICT const _data_force = force->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
  double *RESTRICT _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
  const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  dim3 _block(uint32_t(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)), uint32_t(((1024 < ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))), uint32_t(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))))));
  dim3 _grid(uint32_t(((_size_force_0 - 2) % (((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)) == 0 ? (int64_t)(_size_force_0 - 2) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)) : ((int64_t)(_size_force_0 - 2) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))) + 1)), uint32_t(((_size_force_1 - 2) % (((1024 < ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))) == 0 ? (int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))) : ((int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) + 1)), uint32_t(((_size_force_2 - 2) % (((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))))) == 0 ? (int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))))) : ((int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2) * ((_size_force_1 - 2 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0 - 2) ? 16 : _size_force_0 - 2)))))))))) + 1)));
  internal_streamsweepdoubleprecisioncuda_streamsweepdoubleprecisioncuda::streamsweepdoubleprecisioncuda_streamsweepdoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_pdfs_tmp, _data_velocity, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
  pdfs->swapDataPointers(pdfs_tmp);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
