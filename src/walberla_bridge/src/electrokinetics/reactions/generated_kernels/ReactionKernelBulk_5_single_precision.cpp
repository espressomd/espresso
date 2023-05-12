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
//! \\file ReactionKernelBulk_5_single_precision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.2, lbmpy_walberla/pystencils_walberla from waLBerla commit ref: a839fac6ef7d0c58e7710e4d50490e9dd7146b4a

#include <cmath>

#include "ReactionKernelBulk_5_single_precision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

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

namespace internal_737d6904f7c65dcfc320d619f189641b {
static FUNC_PREFIX void reactionkernelbulk_5_single_precision_reactionkernelbulk_5_single_precision(float *RESTRICT _data_rho_0, float *RESTRICT _data_rho_1, float *RESTRICT _data_rho_2, float *RESTRICT _data_rho_3, float *RESTRICT _data_rho_4, int64_t const _size_rho_0_0, int64_t const _size_rho_0_1, int64_t const _size_rho_0_2, int64_t const _stride_rho_0_0, int64_t const _stride_rho_0_1, int64_t const _stride_rho_0_2, int64_t const _stride_rho_1_0, int64_t const _stride_rho_1_1, int64_t const _stride_rho_1_2, int64_t const _stride_rho_2_0, int64_t const _stride_rho_2_1, int64_t const _stride_rho_2_2, int64_t const _stride_rho_3_0, int64_t const _stride_rho_3_1, int64_t const _stride_rho_3_2, int64_t const _stride_rho_4_0, int64_t const _stride_rho_4_1, int64_t const _stride_rho_4_2, float order_0, float order_1, float order_2, float order_3, float order_4, float rate_coefficient, float stoech_0, float stoech_1, float stoech_2, float stoech_3, float stoech_4) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_rho_0_2; ctr_2 += 1) {
    float *RESTRICT _data_rho_0_20 = _data_rho_0 + _stride_rho_0_2 * ctr_2;
    float *RESTRICT _data_rho_1_20 = _data_rho_1 + _stride_rho_1_2 * ctr_2;
    float *RESTRICT _data_rho_2_20 = _data_rho_2 + _stride_rho_2_2 * ctr_2;
    float *RESTRICT _data_rho_3_20 = _data_rho_3 + _stride_rho_3_2 * ctr_2;
    float *RESTRICT _data_rho_4_20 = _data_rho_4 + _stride_rho_4_2 * ctr_2;
    for (int64_t ctr_1 = 0; ctr_1 < _size_rho_0_1; ctr_1 += 1) {
      float *RESTRICT _data_rho_0_20_10 = _stride_rho_0_1 * ctr_1 + _data_rho_0_20;
      float *RESTRICT _data_rho_1_20_10 = _stride_rho_1_1 * ctr_1 + _data_rho_1_20;
      float *RESTRICT _data_rho_2_20_10 = _stride_rho_2_1 * ctr_1 + _data_rho_2_20;
      float *RESTRICT _data_rho_3_20_10 = _stride_rho_3_1 * ctr_1 + _data_rho_3_20;
      float *RESTRICT _data_rho_4_20_10 = _stride_rho_4_1 * ctr_1 + _data_rho_4_20;
      for (int64_t ctr_0 = 0; ctr_0 < _size_rho_0_0; ctr_0 += 1) {
        const float local_rho_0 = _data_rho_0_20_10[_stride_rho_0_0 * ctr_0];
        const float local_rho_1 = _data_rho_1_20_10[_stride_rho_1_0 * ctr_0];
        const float local_rho_2 = _data_rho_2_20_10[_stride_rho_2_0 * ctr_0];
        const float local_rho_3 = _data_rho_3_20_10[_stride_rho_3_0 * ctr_0];
        const float local_rho_4 = _data_rho_4_20_10[_stride_rho_4_0 * ctr_0];
        const float rate_factor = rate_coefficient * powf(local_rho_0, order_0) * powf(local_rho_1, order_1) * powf(local_rho_2, order_2) * powf(local_rho_3, order_3) * powf(local_rho_4, order_4);
        _data_rho_0_20_10[_stride_rho_0_0 * ctr_0] = local_rho_0 + rate_factor * stoech_0;
        _data_rho_1_20_10[_stride_rho_1_0 * ctr_0] = local_rho_1 + rate_factor * stoech_1;
        _data_rho_2_20_10[_stride_rho_2_0 * ctr_0] = local_rho_2 + rate_factor * stoech_2;
        _data_rho_3_20_10[_stride_rho_3_0 * ctr_0] = local_rho_3 + rate_factor * stoech_3;
        _data_rho_4_20_10[_stride_rho_4_0 * ctr_0] = local_rho_4 + rate_factor * stoech_4;
      }
    }
  }
}
} // namespace internal_737d6904f7c65dcfc320d619f189641b

void ReactionKernelBulk_5_single_precision::run(IBlock *block) {
  auto rho_0 = block->getData<field::GhostLayerField<float, 1>>(rho_0ID);
  auto rho_2 = block->getData<field::GhostLayerField<float, 1>>(rho_2ID);
  auto rho_3 = block->getData<field::GhostLayerField<float, 1>>(rho_3ID);
  auto rho_4 = block->getData<field::GhostLayerField<float, 1>>(rho_4ID);
  auto rho_1 = block->getData<field::GhostLayerField<float, 1>>(rho_1ID);

  auto &order_1 = this->order_1_;
  auto &order_4 = this->order_4_;
  auto &order_0 = this->order_0_;
  auto &stoech_3 = this->stoech_3_;
  auto &stoech_4 = this->stoech_4_;
  auto &stoech_2 = this->stoech_2_;
  auto &stoech_0 = this->stoech_0_;
  auto &order_2 = this->order_2_;
  auto &stoech_1 = this->stoech_1_;
  auto &rate_coefficient = this->rate_coefficient_;
  auto &order_3 = this->order_3_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_0->nrOfGhostLayers()));
  float *RESTRICT _data_rho_0 = rho_0->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_1->nrOfGhostLayers()));
  float *RESTRICT _data_rho_1 = rho_1->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_2->nrOfGhostLayers()));
  float *RESTRICT _data_rho_2 = rho_2->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_3->nrOfGhostLayers()));
  float *RESTRICT _data_rho_3 = rho_3->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_4->nrOfGhostLayers()));
  float *RESTRICT _data_rho_4 = rho_4->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->xSizeWithGhostLayer(), int64_t(cell_idx_c(rho_0->xSize()) + 0));
  const int64_t _size_rho_0_0 = int64_t(cell_idx_c(rho_0->xSize()) + 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->ySizeWithGhostLayer(), int64_t(cell_idx_c(rho_0->ySize()) + 0));
  const int64_t _size_rho_0_1 = int64_t(cell_idx_c(rho_0->ySize()) + 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->zSizeWithGhostLayer(), int64_t(cell_idx_c(rho_0->zSize()) + 0));
  const int64_t _size_rho_0_2 = int64_t(cell_idx_c(rho_0->zSize()) + 0);
  const int64_t _stride_rho_0_0 = int64_t(rho_0->xStride());
  const int64_t _stride_rho_0_1 = int64_t(rho_0->yStride());
  const int64_t _stride_rho_0_2 = int64_t(rho_0->zStride());
  const int64_t _stride_rho_1_0 = int64_t(rho_1->xStride());
  const int64_t _stride_rho_1_1 = int64_t(rho_1->yStride());
  const int64_t _stride_rho_1_2 = int64_t(rho_1->zStride());
  const int64_t _stride_rho_2_0 = int64_t(rho_2->xStride());
  const int64_t _stride_rho_2_1 = int64_t(rho_2->yStride());
  const int64_t _stride_rho_2_2 = int64_t(rho_2->zStride());
  const int64_t _stride_rho_3_0 = int64_t(rho_3->xStride());
  const int64_t _stride_rho_3_1 = int64_t(rho_3->yStride());
  const int64_t _stride_rho_3_2 = int64_t(rho_3->zStride());
  const int64_t _stride_rho_4_0 = int64_t(rho_4->xStride());
  const int64_t _stride_rho_4_1 = int64_t(rho_4->yStride());
  const int64_t _stride_rho_4_2 = int64_t(rho_4->zStride());
  internal_737d6904f7c65dcfc320d619f189641b::reactionkernelbulk_5_single_precision_reactionkernelbulk_5_single_precision(_data_rho_0, _data_rho_1, _data_rho_2, _data_rho_3, _data_rho_4, _size_rho_0_0, _size_rho_0_1, _size_rho_0_2, _stride_rho_0_0, _stride_rho_0_1, _stride_rho_0_2, _stride_rho_1_0, _stride_rho_1_1, _stride_rho_1_2, _stride_rho_2_0, _stride_rho_2_1, _stride_rho_2_2, _stride_rho_3_0, _stride_rho_3_1, _stride_rho_3_2, _stride_rho_4_0, _stride_rho_4_1, _stride_rho_4_2, order_0, order_1, order_2, order_3, order_4, rate_coefficient, stoech_0, stoech_1, stoech_2, stoech_3, stoech_4);
}

void ReactionKernelBulk_5_single_precision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto rho_0 = block->getData<field::GhostLayerField<float, 1>>(rho_0ID);
  auto rho_2 = block->getData<field::GhostLayerField<float, 1>>(rho_2ID);
  auto rho_3 = block->getData<field::GhostLayerField<float, 1>>(rho_3ID);
  auto rho_4 = block->getData<field::GhostLayerField<float, 1>>(rho_4ID);
  auto rho_1 = block->getData<field::GhostLayerField<float, 1>>(rho_1ID);

  auto &order_1 = this->order_1_;
  auto &order_4 = this->order_4_;
  auto &order_0 = this->order_0_;
  auto &stoech_3 = this->stoech_3_;
  auto &stoech_4 = this->stoech_4_;
  auto &stoech_2 = this->stoech_2_;
  auto &stoech_0 = this->stoech_0_;
  auto &order_2 = this->order_2_;
  auto &stoech_1 = this->stoech_1_;
  auto &rate_coefficient = this->rate_coefficient_;
  auto &order_3 = this->order_3_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_0->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_0->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_0->nrOfGhostLayers()));
  float *RESTRICT _data_rho_0 = rho_0->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_1->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_1->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_1->nrOfGhostLayers()));
  float *RESTRICT _data_rho_1 = rho_1->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_2->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_2->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_2->nrOfGhostLayers()));
  float *RESTRICT _data_rho_2 = rho_2->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_3->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_3->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_3->nrOfGhostLayers()));
  float *RESTRICT _data_rho_3 = rho_3->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_4->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_4->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_4->nrOfGhostLayers()));
  float *RESTRICT _data_rho_4 = rho_4->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_rho_0_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_rho_0_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
  const int64_t _size_rho_0_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
  const int64_t _stride_rho_0_0 = int64_t(rho_0->xStride());
  const int64_t _stride_rho_0_1 = int64_t(rho_0->yStride());
  const int64_t _stride_rho_0_2 = int64_t(rho_0->zStride());
  const int64_t _stride_rho_1_0 = int64_t(rho_1->xStride());
  const int64_t _stride_rho_1_1 = int64_t(rho_1->yStride());
  const int64_t _stride_rho_1_2 = int64_t(rho_1->zStride());
  const int64_t _stride_rho_2_0 = int64_t(rho_2->xStride());
  const int64_t _stride_rho_2_1 = int64_t(rho_2->yStride());
  const int64_t _stride_rho_2_2 = int64_t(rho_2->zStride());
  const int64_t _stride_rho_3_0 = int64_t(rho_3->xStride());
  const int64_t _stride_rho_3_1 = int64_t(rho_3->yStride());
  const int64_t _stride_rho_3_2 = int64_t(rho_3->zStride());
  const int64_t _stride_rho_4_0 = int64_t(rho_4->xStride());
  const int64_t _stride_rho_4_1 = int64_t(rho_4->yStride());
  const int64_t _stride_rho_4_2 = int64_t(rho_4->zStride());
  internal_737d6904f7c65dcfc320d619f189641b::reactionkernelbulk_5_single_precision_reactionkernelbulk_5_single_precision(_data_rho_0, _data_rho_1, _data_rho_2, _data_rho_3, _data_rho_4, _size_rho_0_0, _size_rho_0_1, _size_rho_0_2, _stride_rho_0_0, _stride_rho_0_1, _stride_rho_0_2, _stride_rho_1_0, _stride_rho_1_1, _stride_rho_1_2, _stride_rho_2_0, _stride_rho_2_1, _stride_rho_2_2, _stride_rho_3_0, _stride_rho_3_1, _stride_rho_3_2, _stride_rho_4_0, _stride_rho_4_1, _stride_rho_4_2, order_0, order_1, order_2, order_3, order_4, rate_coefficient, stoech_0, stoech_1, stoech_2, stoech_3, stoech_4);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif