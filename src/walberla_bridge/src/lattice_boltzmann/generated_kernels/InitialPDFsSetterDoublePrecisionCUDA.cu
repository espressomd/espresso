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
//! \\file InitialPDFsSetterDoublePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.2, lbmpy_walberla/pystencils_walberla from waLBerla commit 0c8b4b926c6979288fd8a6846d02ec0870e1fe41

#include <cmath>

#include "InitialPDFsSetterDoublePrecisionCUDA.h"
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

namespace internal_initialpdfssetterdoubleprecisioncuda_initialpdfssetterdoubleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void initialpdfssetterdoubleprecisioncuda_initialpdfssetterdoubleprecisioncuda(double *RESTRICT const _data_force, double *RESTRICT _data_pdfs, double *RESTRICT const _data_velocity, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, double rho_0) {
  if (blockDim.x * blockIdx.x + threadIdx.x < _size_force_0 && blockDim.y * blockIdx.y + threadIdx.y < _size_force_1 && blockDim.z * blockIdx.z + threadIdx.z < _size_force_2) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z;
    const double rho = rho_0;
    double *RESTRICT _data_force_10_20_30 = _data_force + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2;
    double *RESTRICT _data_velocity_10_20_30 = _data_velocity + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2;
    const double u_0 = -0.5 * ((1.0) / (rho)) * _data_force_10_20_30[_stride_force_0 * ctr_0] + _data_velocity_10_20_30[_stride_velocity_0 * ctr_0];
    double *RESTRICT _data_force_10_20_31 = _data_force + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_velocity_10_20_31 = _data_velocity + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    const double u_1 = -0.5 * ((1.0) / (rho)) * _data_force_10_20_31[_stride_force_0 * ctr_0] + _data_velocity_10_20_31[_stride_velocity_0 * ctr_0];
    double *RESTRICT _data_force_10_20_32 = _data_force + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_velocity_10_20_32 = _data_velocity + _stride_velocity_1 * ctr_1 + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    const double u_2 = -0.5 * ((1.0) / (rho)) * _data_force_10_20_32[_stride_force_0 * ctr_0] + _data_velocity_10_20_32[_stride_velocity_0 * ctr_0];
    double *RESTRICT _data_pdfs_10_20_30 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2;
    _data_pdfs_10_20_30[_stride_pdfs_0 * ctr_0] = rho * -0.33333333333333331 * (u_0 * u_0) + rho * -0.33333333333333331 * (u_1 * u_1) + rho * -0.33333333333333331 * (u_2 * u_2) + rho * 0.33333333333333331;
    double *RESTRICT _data_pdfs_10_20_31 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    _data_pdfs_10_20_31[_stride_pdfs_0 * ctr_0] = rho * u_1 * 0.16666666666666666 + rho * -0.16666666666666666 * (u_0 * u_0) + rho * -0.16666666666666666 * (u_2 * u_2) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u_1 * u_1);
    double *RESTRICT _data_pdfs_10_20_32 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    _data_pdfs_10_20_32[_stride_pdfs_0 * ctr_0] = rho * u_1 * -0.16666666666666666 + rho * -0.16666666666666666 * (u_0 * u_0) + rho * -0.16666666666666666 * (u_2 * u_2) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u_1 * u_1);
    double *RESTRICT _data_pdfs_10_20_33 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    _data_pdfs_10_20_33[_stride_pdfs_0 * ctr_0] = rho * u_0 * -0.16666666666666666 + rho * -0.16666666666666666 * (u_1 * u_1) + rho * -0.16666666666666666 * (u_2 * u_2) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u_0 * u_0);
    double *RESTRICT _data_pdfs_10_20_34 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    _data_pdfs_10_20_34[_stride_pdfs_0 * ctr_0] = rho * u_0 * 0.16666666666666666 + rho * -0.16666666666666666 * (u_1 * u_1) + rho * -0.16666666666666666 * (u_2 * u_2) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u_0 * u_0);
    double *RESTRICT _data_pdfs_10_20_35 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    _data_pdfs_10_20_35[_stride_pdfs_0 * ctr_0] = rho * u_2 * 0.16666666666666666 + rho * -0.16666666666666666 * (u_0 * u_0) + rho * -0.16666666666666666 * (u_1 * u_1) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_36 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    _data_pdfs_10_20_36[_stride_pdfs_0 * ctr_0] = rho * u_2 * -0.16666666666666666 + rho * -0.16666666666666666 * (u_0 * u_0) + rho * -0.16666666666666666 * (u_1 * u_1) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_37 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    _data_pdfs_10_20_37[_stride_pdfs_0 * ctr_0] = rho * u_0 * u_1 * -0.25 + rho * u_0 * -0.083333333333333329 + rho * u_1 * 0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_0 * u_0) + rho * 0.083333333333333329 * (u_1 * u_1);
    double *RESTRICT _data_pdfs_10_20_38 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    _data_pdfs_10_20_38[_stride_pdfs_0 * ctr_0] = rho * u_0 * u_1 * 0.25 + rho * u_0 * 0.083333333333333329 + rho * u_1 * 0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_0 * u_0) + rho * 0.083333333333333329 * (u_1 * u_1);
    double *RESTRICT _data_pdfs_10_20_39 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    _data_pdfs_10_20_39[_stride_pdfs_0 * ctr_0] = rho * u_0 * u_1 * 0.25 + rho * u_0 * -0.083333333333333329 + rho * u_1 * -0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_0 * u_0) + rho * 0.083333333333333329 * (u_1 * u_1);
    double *RESTRICT _data_pdfs_10_20_310 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    _data_pdfs_10_20_310[_stride_pdfs_0 * ctr_0] = rho * u_0 * u_1 * -0.25 + rho * u_0 * 0.083333333333333329 + rho * u_1 * -0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_0 * u_0) + rho * 0.083333333333333329 * (u_1 * u_1);
    double *RESTRICT _data_pdfs_10_20_311 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    _data_pdfs_10_20_311[_stride_pdfs_0 * ctr_0] = rho * u_1 * u_2 * 0.25 + rho * u_1 * 0.083333333333333329 + rho * u_2 * 0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_1 * u_1) + rho * 0.083333333333333329 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_312 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    _data_pdfs_10_20_312[_stride_pdfs_0 * ctr_0] = rho * u_1 * u_2 * -0.25 + rho * u_1 * -0.083333333333333329 + rho * u_2 * 0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_1 * u_1) + rho * 0.083333333333333329 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_313 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    _data_pdfs_10_20_313[_stride_pdfs_0 * ctr_0] = rho * u_0 * u_2 * -0.25 + rho * u_0 * -0.083333333333333329 + rho * u_2 * 0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_0 * u_0) + rho * 0.083333333333333329 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_314 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    _data_pdfs_10_20_314[_stride_pdfs_0 * ctr_0] = rho * u_0 * u_2 * 0.25 + rho * u_0 * 0.083333333333333329 + rho * u_2 * 0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_0 * u_0) + rho * 0.083333333333333329 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_315 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    _data_pdfs_10_20_315[_stride_pdfs_0 * ctr_0] = rho * u_1 * u_2 * -0.25 + rho * u_1 * 0.083333333333333329 + rho * u_2 * -0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_1 * u_1) + rho * 0.083333333333333329 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_316 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    _data_pdfs_10_20_316[_stride_pdfs_0 * ctr_0] = rho * u_1 * u_2 * 0.25 + rho * u_1 * -0.083333333333333329 + rho * u_2 * -0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_1 * u_1) + rho * 0.083333333333333329 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_317 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    _data_pdfs_10_20_317[_stride_pdfs_0 * ctr_0] = rho * u_0 * u_2 * 0.25 + rho * u_0 * -0.083333333333333329 + rho * u_2 * -0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_0 * u_0) + rho * 0.083333333333333329 * (u_2 * u_2);
    double *RESTRICT _data_pdfs_10_20_318 = _data_pdfs + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    _data_pdfs_10_20_318[_stride_pdfs_0 * ctr_0] = rho * u_0 * u_2 * -0.25 + rho * u_0 * 0.083333333333333329 + rho * u_2 * -0.083333333333333329 + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u_0 * u_0) + rho * 0.083333333333333329 * (u_2 * u_2);
  }
}
} // namespace internal_initialpdfssetterdoubleprecisioncuda_initialpdfssetterdoubleprecisioncuda

void InitialPDFsSetterDoublePrecisionCUDA::run(IBlock *block, gpuStream_t stream) {
  auto force = block->getData<gpu::GPUField<double>>(forceID);
  auto velocity = block->getData<gpu::GPUField<double>>(velocityID);
  auto pdfs = block->getData<gpu::GPUField<double>>(pdfsID);

  auto &rho_0 = this->rho_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()))
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()))
  double *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(force->xSize()) + 0))
  const int64_t _size_force_0 = int64_t(int64_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(force->ySize()) + 0))
  const int64_t _size_force_1 = int64_t(int64_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(force->zSize()) + 0))
  const int64_t _size_force_2 = int64_t(int64_c(force->zSize()) + 0);
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
  dim3 _block(uint32_t(((16 < _size_force_0) ? 16 : _size_force_0)), uint32_t(((1024 < ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))) ? 1024 : ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))), uint32_t(((64 < ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))) ? 64 : ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))))));
  dim3 _grid(uint32_t(((_size_force_0) % (((16 < _size_force_0) ? 16 : _size_force_0)) == 0 ? (int64_t)(_size_force_0) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)) : ((int64_t)(_size_force_0) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))) + 1)), uint32_t(((_size_force_1) % (((1024 < ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))) ? 1024 : ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))) == 0 ? (int64_t)(_size_force_1) / (int64_t)(((1024 < ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))) ? 1024 : ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))) : ((int64_t)(_size_force_1) / (int64_t)(((1024 < ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))) ? 1024 : ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) + 1)), uint32_t(((_size_force_2) % (((64 < ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))) ? 64 : ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))))) == 0 ? (int64_t)(_size_force_2) / (int64_t)(((64 < ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))) ? 64 : ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))))) : ((int64_t)(_size_force_2) / (int64_t)(((64 < ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))) ? 64 : ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))))) + 1)));
  internal_initialpdfssetterdoubleprecisioncuda_initialpdfssetterdoubleprecisioncuda::initialpdfssetterdoubleprecisioncuda_initialpdfssetterdoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, rho_0);
}

void InitialPDFsSetterDoublePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block, gpuStream_t stream) {
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

  auto &rho_0 = this->rho_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()))
  double *RESTRICT const _data_force = force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()))
  double *RESTRICT const _data_velocity = velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
  const int64_t _size_force_0 = int64_t(int64_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
  const int64_t _size_force_1 = int64_t(int64_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
  const int64_t _size_force_2 = int64_t(int64_c(ci.zSize()) + 0);
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
  dim3 _block(uint32_t(((16 < _size_force_0) ? 16 : _size_force_0)), uint32_t(((1024 < ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))) ? 1024 : ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))), uint32_t(((64 < ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))) ? 64 : ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))))));
  dim3 _grid(uint32_t(((_size_force_0) % (((16 < _size_force_0) ? 16 : _size_force_0)) == 0 ? (int64_t)(_size_force_0) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)) : ((int64_t)(_size_force_0) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))) + 1)), uint32_t(((_size_force_1) % (((1024 < ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))) ? 1024 : ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))) == 0 ? (int64_t)(_size_force_1) / (int64_t)(((1024 < ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))) ? 1024 : ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))) : ((int64_t)(_size_force_1) / (int64_t)(((1024 < ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))) ? 1024 : ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) + 1)), uint32_t(((_size_force_2) % (((64 < ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))) ? 64 : ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))))) == 0 ? (int64_t)(_size_force_2) / (int64_t)(((64 < ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))) ? 64 : ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))))) : ((int64_t)(_size_force_2) / (int64_t)(((64 < ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))) ? 64 : ((_size_force_2 < ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0))))))) ? _size_force_2 : ((int64_t)(256) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0) * ((_size_force_1 < 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))) ? _size_force_1 : 16 * ((int64_t)(16) / (int64_t)(((16 < _size_force_0) ? 16 : _size_force_0)))))))))) + 1)));
  internal_initialpdfssetterdoubleprecisioncuda_initialpdfssetterdoubleprecisioncuda::initialpdfssetterdoubleprecisioncuda_initialpdfssetterdoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, rho_0);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
