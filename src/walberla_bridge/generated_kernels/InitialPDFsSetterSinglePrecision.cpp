// kernel generated with pystencils v0.4.3+12.g29e0e84, lbmpy v0.4.3+2.g0e17e61,
// lbmpy_walberla/pystencils_walberla from commit
// 55e6cf598e7e55f496ffaecd40bde632de76930e

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
//! \\file InitialPDFsSetterSinglePrecision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "InitialPDFsSetterSinglePrecision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
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

namespace internal_b8085d63d6b7e842485134abbac511e8 {
static FUNC_PREFIX void
initialpdfssettersingleprecision_initialpdfssettersingleprecision(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    float *RESTRICT const _data_velocity, int64_t const _size_force_0,
    int64_t const _size_force_1, int64_t const _size_force_2,
    int64_t const _stride_force_0, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t const _stride_velocity_0, int64_t const _stride_velocity_1,
    int64_t const _stride_velocity_2, int64_t const _stride_velocity_3,
    float rho_0) {
  const float rho = rho_0;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const float u_0 =
            _data_velocity_20_30_10[_stride_velocity_0 * ctr_0] -
            0.5f * _data_force_20_30_10[_stride_force_0 * ctr_0] / rho_0;
        const float u_1 =
            _data_velocity_20_31_10[_stride_velocity_0 * ctr_0] -
            0.5f * _data_force_20_31_10[_stride_force_0 * ctr_0] / rho_0;
        const float u_2 =
            _data_velocity_20_32_10[_stride_velocity_0 * ctr_0] -
            0.5f * _data_force_20_32_10[_stride_force_0 * ctr_0] / rho_0;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.5f + rho * (u_1 * u_1) * -0.5f +
            rho * (u_2 * u_2) * -0.5f + rho * 0.333333333333333f;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.0833333333333333f +
            rho * (u_1 * u_1) * 0.166666666666667f +
            rho * u_1 * 0.166666666666667f +
            rho * (u_2 * u_2) * -0.0833333333333333f +
            rho * 0.0555555555555556f;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.0833333333333333f +
            rho * (u_1 * u_1) * 0.166666666666667f +
            rho * u_1 * -0.166666666666667f +
            rho * (u_2 * u_2) * -0.0833333333333333f +
            rho * 0.0555555555555556f;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.166666666666667f +
            rho * u_0 * -0.166666666666667f +
            rho * (u_1 * u_1) * -0.0833333333333333f +
            rho * (u_2 * u_2) * -0.0833333333333333f +
            rho * 0.0555555555555556f;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.166666666666667f +
            rho * u_0 * 0.166666666666667f +
            rho * (u_1 * u_1) * -0.0833333333333333f +
            rho * (u_2 * u_2) * -0.0833333333333333f +
            rho * 0.0555555555555556f;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.0833333333333333f +
            rho * (u_1 * u_1) * -0.0833333333333333f +
            rho * (u_2 * u_2) * 0.166666666666667f +
            rho * u_2 * 0.166666666666667f + rho * 0.0555555555555556f;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.0833333333333333f +
            rho * (u_1 * u_1) * -0.0833333333333333f +
            rho * (u_2 * u_2) * 0.166666666666667f +
            rho * u_2 * -0.166666666666667f + rho * 0.0555555555555556f;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.0833333333333333f + rho * u_0 * u_1 * -0.25f +
            rho * u_0 * -0.0833333333333333f +
            rho * (u_1 * u_1) * 0.0833333333333333f +
            rho * u_1 * 0.0833333333333333f +
            rho * (u_2 * u_2) * -0.0416666666666667f +
            rho * 0.0277777777777778f;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.0833333333333333f + rho * u_0 * u_1 * 0.25f +
            rho * u_0 * 0.0833333333333333f +
            rho * (u_1 * u_1) * 0.0833333333333333f +
            rho * u_1 * 0.0833333333333333f +
            rho * (u_2 * u_2) * -0.0416666666666667f +
            rho * 0.0277777777777778f;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.0833333333333333f + rho * u_0 * u_1 * 0.25f +
            rho * u_0 * -0.0833333333333333f +
            rho * (u_1 * u_1) * 0.0833333333333333f +
            rho * u_1 * -0.0833333333333333f +
            rho * (u_2 * u_2) * -0.0416666666666667f +
            rho * 0.0277777777777778f;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.0833333333333333f + rho * u_0 * u_1 * -0.25f +
            rho * u_0 * 0.0833333333333333f +
            rho * (u_1 * u_1) * 0.0833333333333333f +
            rho * u_1 * -0.0833333333333333f +
            rho * (u_2 * u_2) * -0.0416666666666667f +
            rho * 0.0277777777777778f;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.0416666666666667f +
            rho * (u_1 * u_1) * 0.0833333333333333f + rho * u_1 * u_2 * 0.25f +
            rho * u_1 * 0.0833333333333333f +
            rho * (u_2 * u_2) * 0.0833333333333333f +
            rho * u_2 * 0.0833333333333333f + rho * 0.0277777777777778f;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.0416666666666667f +
            rho * (u_1 * u_1) * 0.0833333333333333f + rho * u_1 * u_2 * -0.25f +
            rho * u_1 * -0.0833333333333333f +
            rho * (u_2 * u_2) * 0.0833333333333333f +
            rho * u_2 * 0.0833333333333333f + rho * 0.0277777777777778f;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.0833333333333333f + rho * u_0 * u_2 * -0.25f +
            rho * u_0 * -0.0833333333333333f +
            rho * (u_1 * u_1) * -0.0416666666666667f +
            rho * (u_2 * u_2) * 0.0833333333333333f +
            rho * u_2 * 0.0833333333333333f + rho * 0.0277777777777778f;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.0833333333333333f + rho * u_0 * u_2 * 0.25f +
            rho * u_0 * 0.0833333333333333f +
            rho * (u_1 * u_1) * -0.0416666666666667f +
            rho * (u_2 * u_2) * 0.0833333333333333f +
            rho * u_2 * 0.0833333333333333f + rho * 0.0277777777777778f;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.0416666666666667f +
            rho * (u_1 * u_1) * 0.0833333333333333f + rho * u_1 * u_2 * -0.25f +
            rho * u_1 * 0.0833333333333333f +
            rho * (u_2 * u_2) * 0.0833333333333333f +
            rho * u_2 * -0.0833333333333333f + rho * 0.0277777777777778f;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * -0.0416666666666667f +
            rho * (u_1 * u_1) * 0.0833333333333333f + rho * u_1 * u_2 * 0.25f +
            rho * u_1 * -0.0833333333333333f +
            rho * (u_2 * u_2) * 0.0833333333333333f +
            rho * u_2 * -0.0833333333333333f + rho * 0.0277777777777778f;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.0833333333333333f + rho * u_0 * u_2 * 0.25f +
            rho * u_0 * -0.0833333333333333f +
            rho * (u_1 * u_1) * -0.0416666666666667f +
            rho * (u_2 * u_2) * 0.0833333333333333f +
            rho * u_2 * -0.0833333333333333f + rho * 0.0277777777777778f;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            rho * (u_0 * u_0) * 0.0833333333333333f + rho * u_0 * u_2 * -0.25f +
            rho * u_0 * 0.0833333333333333f +
            rho * (u_1 * u_1) * -0.0416666666666667f +
            rho * (u_2 * u_2) * 0.0833333333333333f +
            rho * u_2 * -0.0833333333333333f + rho * 0.0277777777777778f;
      }
    }
  }
}
} // namespace internal_b8085d63d6b7e842485134abbac511e8

void InitialPDFsSetterSinglePrecision::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);

  auto &rho_0 = this->rho_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
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
  internal_b8085d63d6b7e842485134abbac511e8::
      initialpdfssettersingleprecision_initialpdfssettersingleprecision(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
          _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
          _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, rho_0);
}

void InitialPDFsSetterSinglePrecision::runOnCellInterval(
    const shared_ptr<StructuredBlockStorage> &blocks,
    const CellInterval &globalCellInterval, cell_idx_t ghostLayers,
    IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);

  auto &rho_0 = this->rho_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity =
      velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
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
  internal_b8085d63d6b7e842485134abbac511e8::
      initialpdfssettersingleprecision_initialpdfssettersingleprecision(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
          _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
          _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, rho_0);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif